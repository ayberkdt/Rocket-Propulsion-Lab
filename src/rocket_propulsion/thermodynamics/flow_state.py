"""Connections between static, total, and process thermodynamic states."""

from dataclasses import dataclass
from math import isfinite, log, sqrt

from rocket_propulsion.core.errors import DomainError


@dataclass(frozen=True, slots=True)
class StagnationState:
    """Static and stagnation properties for adiabatic perfect-gas flow."""

    mach: float
    static_pressure_pa: float
    static_temperature_k: float
    static_density_kg_m3: float
    stagnation_pressure_pa: float
    stagnation_temperature_k: float
    stagnation_density_kg_m3: float
    speed_of_sound_m_s: float
    velocity_m_s: float
    dynamic_pressure_pa: float
    static_enthalpy_j_kg: float
    stagnation_enthalpy_j_kg: float


@dataclass(frozen=True, slots=True)
class IsentropicProcessResult:
    """Ideal and efficiency-corrected compression or expansion state."""

    process: str
    inlet_pressure_pa: float
    outlet_pressure_pa: float
    inlet_temperature_k: float
    ideal_outlet_temperature_k: float
    actual_outlet_temperature_k: float
    efficiency: float
    specific_work_j_kg: float
    entropy_change_j_kg_k: float


def _positive(value: float, label: str) -> None:
    if not isfinite(value) or value <= 0.0:
        raise DomainError(f"{label} must be finite and greater than zero.")


def calculate_stagnation_state(
    *,
    mach: float,
    static_pressure_pa: float,
    static_temperature_k: float,
    gamma: float = 1.4,
    gas_constant_j_kg_k: float = 287.05,
) -> StagnationState:
    """Connect a static perfect-gas flow state to its stagnation state.

    Uses adiabatic, reversible relations and ``h0 = h + V²/2``. The dynamic
    pressure is the conventional ``rho V²/2`` and is not equal to stagnation
    pressure minus static pressure at finite Mach number.
    """

    for value, label in (
        (mach, "Mach number"),
        (static_pressure_pa, "Static pressure"),
        (static_temperature_k, "Static temperature"),
        (gas_constant_j_kg_k, "Specific gas constant"),
    ):
        _positive(value, label)
    if not isfinite(gamma) or gamma <= 1.0:
        raise DomainError("Gamma must be finite and greater than one.")

    factor = 1.0 + 0.5 * (gamma - 1.0) * mach**2
    static_density = static_pressure_pa / (gas_constant_j_kg_k * static_temperature_k)
    stagnation_temperature = static_temperature_k * factor
    stagnation_pressure = static_pressure_pa * factor ** (gamma / (gamma - 1.0))
    stagnation_density = static_density * factor ** (1.0 / (gamma - 1.0))
    speed_of_sound = sqrt(gamma * gas_constant_j_kg_k * static_temperature_k)
    velocity = mach * speed_of_sound
    cp = gamma * gas_constant_j_kg_k / (gamma - 1.0)
    return StagnationState(
        mach=mach,
        static_pressure_pa=static_pressure_pa,
        static_temperature_k=static_temperature_k,
        static_density_kg_m3=static_density,
        stagnation_pressure_pa=stagnation_pressure,
        stagnation_temperature_k=stagnation_temperature,
        stagnation_density_kg_m3=stagnation_density,
        speed_of_sound_m_s=speed_of_sound,
        velocity_m_s=velocity,
        dynamic_pressure_pa=0.5 * static_density * velocity**2,
        static_enthalpy_j_kg=cp * static_temperature_k,
        stagnation_enthalpy_j_kg=cp * stagnation_temperature,
    )


def calculate_isentropic_process(
    *,
    inlet_pressure_pa: float,
    outlet_pressure_pa: float,
    inlet_temperature_k: float,
    efficiency: float = 1.0,
    gamma: float = 1.4,
    gas_constant_j_kg_k: float = 287.05,
) -> IsentropicProcessResult:
    """Calculate ideal and efficiency-corrected compression or expansion.

    Compressor efficiency is ``(T2s-T1)/(T2-T1)``. Turbine/nozzle efficiency is
    ``(T1-T2)/(T1-T2s)``. Specific work is positive for compression and negative
    for expansion under the work-input sign convention.
    """

    for value, label in (
        (inlet_pressure_pa, "Inlet pressure"),
        (outlet_pressure_pa, "Outlet pressure"),
        (inlet_temperature_k, "Inlet temperature"),
        (gas_constant_j_kg_k, "Specific gas constant"),
    ):
        _positive(value, label)
    if not isfinite(gamma) or gamma <= 1.0:
        raise DomainError("Gamma must be finite and greater than one.")
    if not isfinite(efficiency) or not 0.0 < efficiency <= 1.0:
        raise DomainError("Efficiency must be in the interval (0, 1].")
    if outlet_pressure_pa == inlet_pressure_pa:
        raise DomainError("Inlet and outlet pressures must differ.")

    pressure_ratio = outlet_pressure_pa / inlet_pressure_pa
    ideal_temperature = inlet_temperature_k * pressure_ratio ** ((gamma - 1.0) / gamma)
    if pressure_ratio > 1.0:
        process = "compression"
        actual_temperature = inlet_temperature_k + (
            ideal_temperature - inlet_temperature_k
        ) / efficiency
    else:
        process = "expansion"
        actual_temperature = inlet_temperature_k - efficiency * (
            inlet_temperature_k - ideal_temperature
        )

    cp = gamma * gas_constant_j_kg_k / (gamma - 1.0)
    return IsentropicProcessResult(
        process=process,
        inlet_pressure_pa=inlet_pressure_pa,
        outlet_pressure_pa=outlet_pressure_pa,
        inlet_temperature_k=inlet_temperature_k,
        ideal_outlet_temperature_k=ideal_temperature,
        actual_outlet_temperature_k=actual_temperature,
        efficiency=efficiency,
        specific_work_j_kg=cp * (actual_temperature - inlet_temperature_k),
        entropy_change_j_kg_k=cp * log(actual_temperature / inlet_temperature_k)
        - gas_constant_j_kg_k * log(pressure_ratio),
    )

