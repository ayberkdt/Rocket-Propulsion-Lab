"""Ideal choked converging-diverging rocket-nozzle performance."""

from dataclasses import dataclass
from math import isfinite, sqrt

from rocket_propulsion.core.errors import DomainError

from .isentropic import calculate_isentropic, mach_from_isentropic

STANDARD_GRAVITY_M_S2 = 9.80665


@dataclass(frozen=True, slots=True)
class NozzleResult:
    """One-dimensional design-point state for an ideal choked nozzle."""

    throat_area_m2: float
    exit_area_m2: float
    area_ratio: float
    exit_mach: float
    exit_temperature_k: float
    exit_pressure_pa: float
    exit_velocity_m_s: float
    mass_flow_kg_s: float
    momentum_thrust_n: float
    pressure_thrust_n: float
    thrust_n: float
    characteristic_velocity_m_s: float
    thrust_coefficient: float
    specific_impulse_s: float
    pressure_regime: str


def _positive(value: float, label: str) -> None:
    if not isfinite(value) or value <= 0.0:
        raise DomainError(f"{label} must be finite and greater than zero.")


def calculate_ideal_nozzle(
    *,
    chamber_pressure_pa: float,
    chamber_temperature_k: float,
    throat_area_m2: float,
    area_ratio: float,
    ambient_pressure_pa: float,
    gamma: float = 1.22,
    gas_constant_j_kg_k: float = 355.0,
) -> NozzleResult:
    """Calculate ideal, frozen-property, supersonic C-D nozzle performance.

    Chamber properties are treated as stagnation conditions. The throat is
    choked and the requested exit area ratio uses the supersonic isentropic
    branch. Boundary-layer, divergence, chemistry, separation, and discharge
    losses are not included.
    """

    for value, label in (
        (chamber_pressure_pa, "Chamber pressure"),
        (chamber_temperature_k, "Chamber temperature"),
        (throat_area_m2, "Throat area"),
        (area_ratio, "Area ratio"),
        (gas_constant_j_kg_k, "Specific gas constant"),
    ):
        _positive(value, label)
    if area_ratio < 1.0:
        raise DomainError("Exit-to-throat area ratio must be at least one.")
    if not isfinite(ambient_pressure_pa) or ambient_pressure_pa < 0.0:
        raise DomainError("Ambient pressure must be finite and non-negative.")
    if not isfinite(gamma) or gamma <= 1.0:
        raise DomainError("Gamma must be finite and greater than one.")

    exit_mach = mach_from_isentropic("area_ratio_supersonic", area_ratio, gamma)
    exit_ratios = calculate_isentropic(exit_mach, gamma)
    exit_temperature = chamber_temperature_k * exit_ratios.temperature_total_ratio
    exit_pressure = chamber_pressure_pa * exit_ratios.pressure_total_ratio
    exit_velocity = exit_mach * sqrt(gamma * gas_constant_j_kg_k * exit_temperature)
    mass_flow = (
        chamber_pressure_pa
        * throat_area_m2
        / sqrt(chamber_temperature_k)
        * sqrt(gamma / gas_constant_j_kg_k)
        * (2.0 / (gamma + 1.0)) ** ((gamma + 1.0) / (2.0 * (gamma - 1.0)))
    )
    exit_area = area_ratio * throat_area_m2
    momentum_thrust = mass_flow * exit_velocity
    pressure_thrust = (exit_pressure - ambient_pressure_pa) * exit_area
    thrust = momentum_thrust + pressure_thrust
    if thrust <= 0.0:
        raise DomainError("The ideal nozzle state produces non-positive thrust.")
    characteristic_velocity = chamber_pressure_pa * throat_area_m2 / mass_flow
    if abs(exit_pressure - ambient_pressure_pa) <= 0.01 * max(exit_pressure, 1.0):
        regime = "adapted"
    elif exit_pressure > ambient_pressure_pa:
        regime = "underexpanded"
    else:
        regime = "overexpanded"

    return NozzleResult(
        throat_area_m2=throat_area_m2,
        exit_area_m2=exit_area,
        area_ratio=area_ratio,
        exit_mach=exit_mach,
        exit_temperature_k=exit_temperature,
        exit_pressure_pa=exit_pressure,
        exit_velocity_m_s=exit_velocity,
        mass_flow_kg_s=mass_flow,
        momentum_thrust_n=momentum_thrust,
        pressure_thrust_n=pressure_thrust,
        thrust_n=thrust,
        characteristic_velocity_m_s=characteristic_velocity,
        thrust_coefficient=thrust / (chamber_pressure_pa * throat_area_m2),
        specific_impulse_s=thrust / (mass_flow * STANDARD_GRAVITY_M_S2),
        pressure_regime=regime,
    )

