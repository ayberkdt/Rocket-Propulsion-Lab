"""Closed-system polytropic ideal-gas process paths and energy accounting."""

from dataclasses import dataclass
from math import isfinite, log

from rocket_propulsion.core.errors import DomainError


@dataclass(frozen=True, slots=True)
class PolytropicPoint:
    """One point along a quasistatic polytropic process path."""

    fraction: float
    pressure_pa: float
    temperature_k: float
    specific_volume_m3_kg: float
    entropy_change_j_kg_k: float


@dataclass(frozen=True, slots=True)
class PolytropicResult:
    """End state, first-law terms, and sampled path for ``p v^n = const``."""

    process: str
    exponent: float
    inlet_pressure_pa: float
    outlet_pressure_pa: float
    inlet_temperature_k: float
    outlet_temperature_k: float
    inlet_specific_volume_m3_kg: float
    outlet_specific_volume_m3_kg: float
    boundary_work_by_gas_j_kg: float
    internal_energy_change_j_kg: float
    enthalpy_change_j_kg: float
    heat_into_gas_j_kg: float
    entropy_change_j_kg_k: float
    path: tuple[PolytropicPoint, ...]


def calculate_polytropic_process(
    *,
    inlet_pressure_pa: float,
    outlet_pressure_pa: float,
    inlet_temperature_k: float,
    exponent: float,
    gamma: float = 1.4,
    gas_constant_j_kg_k: float = 287.05,
    point_count: int = 61,
) -> PolytropicResult:
    """Calculate a quasistatic closed-system ideal-gas polytropic process.

    Uses ``p v^n = constant`` and the sign convention ``Q - W_by = delta U``.
    ``n = 1`` is handled as an isothermal path. The model is intended to connect
    thermodynamic diagrams and energy terms, not to replace component efficiency
    models for compressors or turbines.
    """

    for value, label in (
        (inlet_pressure_pa, "Inlet pressure"),
        (outlet_pressure_pa, "Outlet pressure"),
        (inlet_temperature_k, "Inlet temperature"),
        (exponent, "Polytropic exponent"),
        (gas_constant_j_kg_k, "Specific gas constant"),
    ):
        if not isfinite(value) or value <= 0.0:
            raise DomainError(f"{label} must be finite and greater than zero.")
    if inlet_pressure_pa == outlet_pressure_pa:
        raise DomainError("Inlet and outlet pressures must differ.")
    if gamma <= 1.0:
        raise DomainError("Gamma must be greater than one.")
    if not 21 <= point_count <= 401:
        raise DomainError("Point count must be between 21 and 401.")

    pressure_ratio = outlet_pressure_pa / inlet_pressure_pa
    outlet_temperature = inlet_temperature_k * pressure_ratio ** ((exponent - 1.0) / exponent)
    inlet_volume = gas_constant_j_kg_k * inlet_temperature_k / inlet_pressure_pa
    outlet_volume = gas_constant_j_kg_k * outlet_temperature / outlet_pressure_pa
    if abs(exponent - 1.0) < 1.0e-10:
        boundary_work = gas_constant_j_kg_k * inlet_temperature_k * log(
            outlet_volume / inlet_volume
        )
    else:
        boundary_work = (
            outlet_pressure_pa * outlet_volume - inlet_pressure_pa * inlet_volume
        ) / (1.0 - exponent)

    cv = gas_constant_j_kg_k / (gamma - 1.0)
    cp = gamma * cv
    delta_u = cv * (outlet_temperature - inlet_temperature_k)
    delta_h = cp * (outlet_temperature - inlet_temperature_k)
    delta_s = cp * log(outlet_temperature / inlet_temperature_k) - gas_constant_j_kg_k * log(
        pressure_ratio
    )
    heat = delta_u + boundary_work

    path: list[PolytropicPoint] = []
    for index in range(point_count):
        fraction = index / (point_count - 1)
        pressure = inlet_pressure_pa * pressure_ratio**fraction
        temperature = inlet_temperature_k * (pressure / inlet_pressure_pa) ** (
            (exponent - 1.0) / exponent
        )
        volume = gas_constant_j_kg_k * temperature / pressure
        entropy = cp * log(temperature / inlet_temperature_k) - gas_constant_j_kg_k * log(
            pressure / inlet_pressure_pa
        )
        path.append(
            PolytropicPoint(
                fraction=fraction,
                pressure_pa=pressure,
                temperature_k=temperature,
                specific_volume_m3_kg=volume,
                entropy_change_j_kg_k=entropy,
            )
        )

    return PolytropicResult(
        process="compression" if pressure_ratio > 1.0 else "expansion",
        exponent=exponent,
        inlet_pressure_pa=inlet_pressure_pa,
        outlet_pressure_pa=outlet_pressure_pa,
        inlet_temperature_k=inlet_temperature_k,
        outlet_temperature_k=outlet_temperature,
        inlet_specific_volume_m3_kg=inlet_volume,
        outlet_specific_volume_m3_kg=outlet_volume,
        boundary_work_by_gas_j_kg=boundary_work,
        internal_energy_change_j_kg=delta_u,
        enthalpy_change_j_kg=delta_h,
        heat_into_gas_j_kg=heat,
        entropy_change_j_kg_k=delta_s,
        path=tuple(path),
    )

