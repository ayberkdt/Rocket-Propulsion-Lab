"""Gas-side nozzle heat transfer using a source-labelled Bartz correlation."""

from __future__ import annotations

from dataclasses import dataclass
from math import isfinite, pi, sqrt

from rocket_propulsion.core.errors import DomainError
from rocket_propulsion.core.metadata import WarningMessage
from rocket_propulsion.geometry import NozzleContour


@dataclass(frozen=True, slots=True)
class HeatFluxStation:
    """Thermal state and accumulated heat at one contour station."""

    x_m: float
    radius_m: float
    mach: float
    adiabatic_wall_temperature_k: float
    gas_side_coefficient_w_m2_k: float
    heat_flux_w_m2: float
    cumulative_heat_w: float


@dataclass(frozen=True, slots=True)
class BartzHeatTransferResult:
    """Wall heat-flux distribution and correlation provenance."""

    stations: tuple[HeatFluxStation, ...]
    peak_heat_flux_w_m2: float
    peak_location_m: float
    total_heat_load_w: float
    energy_balance_residual_w: float
    correlation: str
    property_model: str
    validity: tuple[str, ...]
    warnings: tuple[WarningMessage, ...]
    source: str
    source_url: str


def calculate_bartz_heat_transfer(
    *,
    contour: NozzleContour,
    chamber_pressure_pa: float,
    chamber_temperature_k: float,
    wall_temperature_k: float,
    characteristic_velocity_m_s: float,
    dynamic_viscosity_pa_s: float,
    specific_heat_j_kg_k: float,
    prandtl_number: float,
    gamma: float,
    throat_radius_of_curvature_m: float,
) -> BartzHeatTransferResult:
    """Evaluate the Bartz gas-side coefficient over a sampled nozzle contour."""

    positive = {
        "chamber_pressure_pa": chamber_pressure_pa,
        "chamber_temperature_k": chamber_temperature_k,
        "wall_temperature_k": wall_temperature_k,
        "characteristic_velocity_m_s": characteristic_velocity_m_s,
        "dynamic_viscosity_pa_s": dynamic_viscosity_pa_s,
        "specific_heat_j_kg_k": specific_heat_j_kg_k,
        "prandtl_number": prandtl_number,
        "throat_radius_of_curvature_m": throat_radius_of_curvature_m,
    }
    for field, value in positive.items():
        if not isfinite(value) or value <= 0.0:
            raise DomainError(f"{field} must be finite and positive.", field=field)
    if not isfinite(gamma) or gamma <= 1.0:
        raise DomainError("Gamma must exceed one.", field="gamma")
    if wall_temperature_k >= chamber_temperature_k:
        raise DomainError(
            "Wall temperature must remain below chamber temperature.",
            field="wall_temperature_k",
        )

    throat_diameter = 2.0 * contour.throat_radius_m
    curvature_ratio = throat_diameter / throat_radius_of_curvature_m
    base_coefficient = (
        0.026
        * dynamic_viscosity_pa_s**0.2
        * specific_heat_j_kg_k
        / (prandtl_number**0.6 * throat_diameter**0.2)
        * (chamber_pressure_pa / characteristic_velocity_m_s) ** 0.8
        * curvature_ratio**0.1
    )
    recovery_factor = prandtl_number ** (1.0 / 3.0)
    records: list[HeatFluxStation] = []
    cumulative = 0.0
    previous_flux = 0.0
    for index, station in enumerate(contour.stations):
        compressibility = 1.0 + 0.5 * (gamma - 1.0) * station.mach**2
        static_temperature = chamber_temperature_k / compressibility
        adiabatic_wall_temperature = static_temperature * (
            1.0 + recovery_factor * 0.5 * (gamma - 1.0) * station.mach**2
        )
        sigma = (
            0.5 * wall_temperature_k / chamber_temperature_k * compressibility + 0.5
        ) ** -0.68 * compressibility**-0.12
        coefficient = base_coefficient * (1.0 / station.area_ratio) ** 0.9 * sigma
        heat_flux = max(0.0, coefficient * (adiabatic_wall_temperature - wall_temperature_k))
        if index:
            previous = contour.stations[index - 1]
            segment_length = sqrt(
                (station.x_m - previous.x_m) ** 2
                + (station.radius_m - previous.radius_m) ** 2
            )
            surface_area = pi * (station.radius_m + previous.radius_m) * segment_length
            cumulative += 0.5 * (previous_flux + heat_flux) * surface_area
        records.append(
            HeatFluxStation(
                x_m=station.x_m,
                radius_m=station.radius_m,
                mach=station.mach,
                adiabatic_wall_temperature_k=adiabatic_wall_temperature,
                gas_side_coefficient_w_m2_k=coefficient,
                heat_flux_w_m2=heat_flux,
                cumulative_heat_w=cumulative,
            )
        )
        previous_flux = heat_flux
    peak = max(records, key=lambda record: record.heat_flux_w_m2)
    warnings: list[WarningMessage] = []
    if not 0.01 <= throat_diameter <= 1.0:
        warnings.append(
            WarningMessage(
                code="bartz_diameter_extrapolation",
                message="Throat diameter is outside the declared preliminary range.",
                severity="warning",
            )
        )
    if chamber_pressure_pa > 30.0e6:
        warnings.append(
            WarningMessage(
                code="bartz_pressure_extrapolation",
                message="Chamber pressure exceeds the declared preliminary range.",
                severity="warning",
            )
        )
    return BartzHeatTransferResult(
        stations=tuple(records),
        peak_heat_flux_w_m2=peak.heat_flux_w_m2,
        peak_location_m=peak.x_m,
        total_heat_load_w=cumulative,
        energy_balance_residual_w=abs(cumulative - records[-1].cumulative_heat_w),
        correlation="Bartz gas-side heat-transfer correlation",
        property_model="user-supplied frozen bulk gas properties",
        validity=(
            "steady turbulent boundary layer",
            "0.01 m <= throat diameter <= 1 m (declared application band)",
            "chamber pressure <= 30 MPa (declared application band)",
        ),
        warnings=tuple(warnings),
        source="Bartz correlation as documented in NASA SP-125",
        source_url="https://ntrs.nasa.gov/citations/19710019929",
    )
