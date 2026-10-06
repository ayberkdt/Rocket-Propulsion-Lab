"""Experimental one-dimensional regenerative cooling channel balance."""

from __future__ import annotations

from dataclasses import dataclass
from math import isfinite, log10

from rocket_propulsion.core.errors import DomainError
from rocket_propulsion.core.metadata import WarningMessage


@dataclass(frozen=True, slots=True)
class RegenerativeCoolingResult:
    """Bulk coolant energy rise and Darcy-Weisbach pressure loss."""

    inlet_temperature_k: float
    outlet_temperature_k: float
    temperature_rise_k: float
    absorbed_heat_w: float
    coolant_capacity_w_k: float
    channel_velocity_m_s: float
    reynolds_number: float
    darcy_friction_factor: float
    pressure_drop_pa: float
    energy_balance_residual_w: float
    experimental: bool
    warnings: tuple[WarningMessage, ...]
    assumptions: tuple[str, ...]


def calculate_regenerative_cooling(
    *,
    absorbed_heat_w: float,
    coolant_mass_flow_kg_s: float,
    coolant_specific_heat_j_kg_k: float,
    coolant_inlet_temperature_k: float,
    coolant_density_kg_m3: float,
    coolant_dynamic_viscosity_pa_s: float,
    channel_count: int,
    channel_flow_area_m2: float,
    hydraulic_diameter_m: float,
    channel_length_m: float,
    roughness_m: float,
    maximum_coolant_temperature_k: float,
) -> RegenerativeCoolingResult:
    """Close a bulk energy balance and Darcy-Weisbach channel pressure loss."""

    values = {
        "absorbed_heat_w": absorbed_heat_w,
        "coolant_mass_flow_kg_s": coolant_mass_flow_kg_s,
        "coolant_specific_heat_j_kg_k": coolant_specific_heat_j_kg_k,
        "coolant_inlet_temperature_k": coolant_inlet_temperature_k,
        "coolant_density_kg_m3": coolant_density_kg_m3,
        "coolant_dynamic_viscosity_pa_s": coolant_dynamic_viscosity_pa_s,
        "channel_flow_area_m2": channel_flow_area_m2,
        "hydraulic_diameter_m": hydraulic_diameter_m,
        "channel_length_m": channel_length_m,
        "maximum_coolant_temperature_k": maximum_coolant_temperature_k,
    }
    for field, value in values.items():
        if not isfinite(value) or value <= 0.0:
            raise DomainError(f"{field} must be finite and positive.", field=field)
    if not isinstance(channel_count, int) or isinstance(channel_count, bool) or channel_count < 1:
        raise DomainError("Channel count must be a positive integer.", field="channel_count")
    if not isfinite(roughness_m) or roughness_m < 0.0:
        raise DomainError("Roughness must be finite and non-negative.", field="roughness_m")
    if maximum_coolant_temperature_k <= coolant_inlet_temperature_k:
        raise DomainError("Maximum coolant temperature must exceed the inlet temperature.")

    capacity = coolant_mass_flow_kg_s * coolant_specific_heat_j_kg_k
    temperature_rise = absorbed_heat_w / capacity
    outlet_temperature = coolant_inlet_temperature_k + temperature_rise
    total_area = channel_count * channel_flow_area_m2
    velocity = coolant_mass_flow_kg_s / (coolant_density_kg_m3 * total_area)
    reynolds = (
        coolant_density_kg_m3
        * velocity
        * hydraulic_diameter_m
        / coolant_dynamic_viscosity_pa_s
    )
    if reynolds < 2300.0:
        friction_factor = 64.0 / reynolds
    else:
        friction_factor = 0.25 / (
            log10(
                roughness_m / (3.7 * hydraulic_diameter_m)
                + 5.74 / reynolds**0.9
            )
            ** 2
        )
    pressure_drop = (
        friction_factor
        * channel_length_m
        / hydraulic_diameter_m
        * 0.5
        * coolant_density_kg_m3
        * velocity**2
    )
    warnings: list[WarningMessage] = []
    if outlet_temperature > maximum_coolant_temperature_k:
        warnings.append(
            WarningMessage(
                code="coolant_temperature_limit_exceeded",
                message="Calculated bulk outlet temperature exceeds the supplied coolant limit.",
                severity="critical",
            )
        )
    if 2300.0 <= reynolds < 4000.0:
        warnings.append(
            WarningMessage(
                code="coolant_transition_regime",
                message="Channel Reynolds number is in the laminar-transition band.",
                severity="warning",
            )
        )
    recovered_heat = capacity * (outlet_temperature - coolant_inlet_temperature_k)
    return RegenerativeCoolingResult(
        inlet_temperature_k=coolant_inlet_temperature_k,
        outlet_temperature_k=outlet_temperature,
        temperature_rise_k=temperature_rise,
        absorbed_heat_w=absorbed_heat_w,
        coolant_capacity_w_k=capacity,
        channel_velocity_m_s=velocity,
        reynolds_number=reynolds,
        darcy_friction_factor=friction_factor,
        pressure_drop_pa=pressure_drop,
        energy_balance_residual_w=abs(absorbed_heat_w - recovered_heat),
        experimental=True,
        warnings=tuple(warnings),
        assumptions=(
            "single-phase constant coolant properties",
            "uniform parallel channels",
            "bulk-energy balance without wall conduction",
            "Darcy-Weisbach pressure loss with Swamee-Jain turbulent friction",
        ),
    )
