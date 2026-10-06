"""Preliminary rocket engine-cycle component energy balances."""

from __future__ import annotations

from dataclasses import dataclass
from math import isfinite

from rocket_propulsion.core.errors import DomainError
from rocket_propulsion.core.metadata import WarningMessage


@dataclass(frozen=True, slots=True)
class EngineCycleResult:
    cycle: str
    chamber_pressure_pa: float
    required_feed_pressure_pa: float
    pump_pressure_rise_pa: float
    pump_power_w: float
    turbine_specific_work_j_kg: float
    turbine_mass_flow_kg_s: float
    turbine_flow_fraction: float
    turbine_power_w: float
    power_balance_residual_w: float
    chamber_mass_flow_kg_s: float
    total_propellant_flow_kg_s: float
    delivered_mass_fraction: float
    warnings: tuple[WarningMessage, ...]
    assumptions: tuple[str, ...]


def calculate_engine_cycle(
    *,
    cycle: str,
    chamber_pressure_pa: float,
    chamber_mass_flow_kg_s: float,
    propellant_density_kg_m3: float,
    inlet_pressure_pa: float,
    injector_pressure_drop_fraction: float,
    pump_efficiency: float,
    turbine_efficiency: float,
    turbine_inlet_temperature_k: float,
    turbine_specific_heat_j_kg_k: float,
    turbine_gamma: float,
    turbine_pressure_ratio: float,
) -> EngineCycleResult:
    """Close feed/pump/turbine power for pressure-fed, GG, or staged combustion."""

    normalized = cycle.strip().lower().replace("_", "-")
    if normalized not in {"pressure-fed", "gas-generator", "staged-combustion"}:
        raise DomainError("Cycle must be pressure-fed, gas-generator, or staged-combustion.")
    positive = {
        "chamber_pressure_pa": chamber_pressure_pa,
        "chamber_mass_flow_kg_s": chamber_mass_flow_kg_s,
        "propellant_density_kg_m3": propellant_density_kg_m3,
        "turbine_inlet_temperature_k": turbine_inlet_temperature_k,
        "turbine_specific_heat_j_kg_k": turbine_specific_heat_j_kg_k,
        "turbine_pressure_ratio": turbine_pressure_ratio,
    }
    for field, value in positive.items():
        if not isfinite(value) or value <= 0.0:
            raise DomainError(f"{field} must be finite and positive.", field=field)
    if not isfinite(inlet_pressure_pa) or inlet_pressure_pa < 0.0:
        raise DomainError("Inlet pressure must be finite and non-negative.")
    if not 0.0 < injector_pressure_drop_fraction <= 1.0:
        raise DomainError("Injector pressure-drop fraction must be in (0, 1].")
    for field, value in {"pump_efficiency": pump_efficiency, "turbine_efficiency": turbine_efficiency}.items():
        if not isfinite(value) or not 0.0 < value <= 1.0:
            raise DomainError(f"{field} must be in (0, 1].", field=field)
    if not isfinite(turbine_gamma) or turbine_gamma <= 1.0:
        raise DomainError("Turbine gamma must exceed one.")
    if turbine_pressure_ratio <= 1.0:
        raise DomainError("Turbine pressure ratio must exceed one.")

    required_feed_pressure = chamber_pressure_pa * (1.0 + injector_pressure_drop_fraction)
    if normalized == "pressure-fed":
        return EngineCycleResult(
            cycle=normalized,
            chamber_pressure_pa=chamber_pressure_pa,
            required_feed_pressure_pa=required_feed_pressure,
            pump_pressure_rise_pa=0.0,
            pump_power_w=0.0,
            turbine_specific_work_j_kg=0.0,
            turbine_mass_flow_kg_s=0.0,
            turbine_flow_fraction=0.0,
            turbine_power_w=0.0,
            power_balance_residual_w=0.0,
            chamber_mass_flow_kg_s=chamber_mass_flow_kg_s,
            total_propellant_flow_kg_s=chamber_mass_flow_kg_s,
            delivered_mass_fraction=1.0,
            warnings=(),
            assumptions=(
                "tank pressure supplies injector pressure directly",
                "feed-line pressure loss excluded except declared injector fraction",
            ),
        )
    pressure_rise = max(0.0, required_feed_pressure - inlet_pressure_pa)
    pump_power = (
        chamber_mass_flow_kg_s
        * pressure_rise
        / (propellant_density_kg_m3 * pump_efficiency)
    )
    exponent = (turbine_gamma - 1.0) / turbine_gamma
    turbine_specific_work = (
        turbine_efficiency
        * turbine_specific_heat_j_kg_k
        * turbine_inlet_temperature_k
        * (1.0 - (1.0 / turbine_pressure_ratio) ** exponent)
    )
    if turbine_specific_work <= 0.0:
        raise DomainError("Turbine conditions do not yield positive specific work.")
    turbine_flow = pump_power / turbine_specific_work
    if normalized == "gas-generator":
        total_flow = chamber_mass_flow_kg_s + turbine_flow
        delivered_fraction = chamber_mass_flow_kg_s / total_flow
    else:
        total_flow = chamber_mass_flow_kg_s
        delivered_fraction = 1.0
        if turbine_flow > chamber_mass_flow_kg_s:
            raise DomainError(
                "Required staged-combustion turbine flow exceeds total chamber flow.",
                code="infeasible_cycle_power_balance",
            )
    turbine_power = turbine_flow * turbine_specific_work
    warnings: list[WarningMessage] = []
    if turbine_flow / total_flow > 0.08 and normalized == "gas-generator":
        warnings.append(
            WarningMessage(
                code="high_gas_generator_flow_fraction",
                message="Required gas-generator flow exceeds 8% of total propellant flow.",
            )
        )
    return EngineCycleResult(
        cycle=normalized,
        chamber_pressure_pa=chamber_pressure_pa,
        required_feed_pressure_pa=required_feed_pressure,
        pump_pressure_rise_pa=pressure_rise,
        pump_power_w=pump_power,
        turbine_specific_work_j_kg=turbine_specific_work,
        turbine_mass_flow_kg_s=turbine_flow,
        turbine_flow_fraction=turbine_flow / total_flow,
        turbine_power_w=turbine_power,
        power_balance_residual_w=abs(turbine_power - pump_power),
        chamber_mass_flow_kg_s=chamber_mass_flow_kg_s,
        total_propellant_flow_kg_s=total_flow,
        delivered_mass_fraction=delivered_fraction,
        warnings=tuple(warnings),
        assumptions=(
            "single equivalent propellant density",
            "constant pump and turbine efficiencies",
            "no pump map, cavitation, gearbox, or bearing loss model",
            "gas-generator turbine flow is dumped; staged-combustion flow returns to chamber",
        ),
    )
