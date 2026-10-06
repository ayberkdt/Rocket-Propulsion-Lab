"""System-level rocket performance relations."""

from dataclasses import dataclass
from math import exp, isfinite, log

from rocket_propulsion.core.errors import DomainError

STANDARD_GRAVITY_M_S2 = 9.80665


@dataclass(frozen=True, slots=True)
class RocketEquationResult:
    """Mass and ideal delta-v quantities from the Tsiolkovsky equation."""

    initial_mass_kg: float
    final_mass_kg: float
    propellant_mass_kg: float
    mass_ratio: float
    propellant_mass_fraction: float
    specific_impulse_s: float
    effective_exhaust_velocity_m_s: float
    delta_v_m_s: float


@dataclass(frozen=True, slots=True)
class ThrustResult:
    """Momentum and pressure contributions to steady rocket thrust."""

    mass_flow_kg_s: float
    exit_velocity_m_s: float
    exit_area_m2: float
    exit_pressure_pa: float
    ambient_pressure_pa: float
    momentum_thrust_n: float
    pressure_thrust_n: float
    total_thrust_n: float
    effective_exhaust_velocity_m_s: float
    specific_impulse_s: float
    total_impulse_n_s: float | None


def _positive(value: float, label: str) -> None:
    if not isfinite(value) or value <= 0.0:
        raise DomainError(f"{label} must be finite and greater than zero.")


def calculate_rocket_equation(
    *,
    initial_mass_kg: float | None = None,
    final_mass_kg: float | None = None,
    delta_v_m_s: float | None = None,
    specific_impulse_s: float | None = None,
    effective_exhaust_velocity_m_s: float | None = None,
) -> RocketEquationResult:
    """Solve the ideal Tsiolkovsky rocket equation in forward or inverse form.

    Provide exactly one performance measure (`specific_impulse_s` or
    `effective_exhaust_velocity_m_s`). Then provide either both masses to solve
    delta-v, or `initial_mass_kg` and `delta_v_m_s` to solve final mass.

    The model uses ``delta_v = v_e ln(m_0/m_f)`` and neglects gravity, drag,
    steering, residual propellant, and finite-burn losses.
    """

    provided_performance = [specific_impulse_s is not None, effective_exhaust_velocity_m_s is not None]
    if sum(provided_performance) != 1:
        raise DomainError("Provide exactly one of specific impulse or effective exhaust velocity.")

    if specific_impulse_s is not None:
        _positive(specific_impulse_s, "Specific impulse")
        effective_velocity = specific_impulse_s * STANDARD_GRAVITY_M_S2
    else:
        assert effective_exhaust_velocity_m_s is not None
        _positive(effective_exhaust_velocity_m_s, "Effective exhaust velocity")
        effective_velocity = effective_exhaust_velocity_m_s
        specific_impulse_s = effective_velocity / STANDARD_GRAVITY_M_S2

    if initial_mass_kg is None:
        raise DomainError("Initial mass is required.")
    _positive(initial_mass_kg, "Initial mass")

    if final_mass_kg is not None and delta_v_m_s is not None:
        raise DomainError("Provide final mass or delta-v, not both.")
    if final_mass_kg is None and delta_v_m_s is None:
        raise DomainError("Provide final mass or delta-v.")

    if final_mass_kg is not None:
        _positive(final_mass_kg, "Final mass")
        if final_mass_kg >= initial_mass_kg:
            raise DomainError("Final mass must be smaller than initial mass.")
        mass_ratio = initial_mass_kg / final_mass_kg
        delta_v_m_s = effective_velocity * log(mass_ratio)
    else:
        assert delta_v_m_s is not None
        _positive(delta_v_m_s, "Delta-v")
        mass_ratio = exp(delta_v_m_s / effective_velocity)
        final_mass_kg = initial_mass_kg / mass_ratio

    propellant_mass = initial_mass_kg - final_mass_kg
    return RocketEquationResult(
        initial_mass_kg=initial_mass_kg,
        final_mass_kg=final_mass_kg,
        propellant_mass_kg=propellant_mass,
        mass_ratio=mass_ratio,
        propellant_mass_fraction=propellant_mass / initial_mass_kg,
        specific_impulse_s=specific_impulse_s,
        effective_exhaust_velocity_m_s=effective_velocity,
        delta_v_m_s=delta_v_m_s,
    )


def calculate_thrust(
    *,
    mass_flow_kg_s: float,
    exit_velocity_m_s: float,
    exit_area_m2: float,
    exit_pressure_pa: float,
    ambient_pressure_pa: float,
    burn_time_s: float | None = None,
) -> ThrustResult:
    """Calculate steady rocket thrust from momentum and pressure terms.

    Implements ``F = mdot V_e + A_e (p_e - p_a)``. A negative pressure term is
    allowed for an overexpanded nozzle, but total thrust must remain positive.
    """

    for value, label in (
        (mass_flow_kg_s, "Mass flow"),
        (exit_velocity_m_s, "Exit velocity"),
        (exit_area_m2, "Exit area"),
        (exit_pressure_pa, "Exit pressure"),
    ):
        _positive(value, label)
    if not isfinite(ambient_pressure_pa) or ambient_pressure_pa < 0.0:
        raise DomainError("Ambient pressure must be finite and non-negative.")
    if burn_time_s is not None:
        _positive(burn_time_s, "Burn time")

    momentum = mass_flow_kg_s * exit_velocity_m_s
    pressure = exit_area_m2 * (exit_pressure_pa - ambient_pressure_pa)
    total = momentum + pressure
    if total <= 0.0:
        raise DomainError("The specified state produces non-positive total thrust.")
    effective_velocity = total / mass_flow_kg_s
    return ThrustResult(
        mass_flow_kg_s=mass_flow_kg_s,
        exit_velocity_m_s=exit_velocity_m_s,
        exit_area_m2=exit_area_m2,
        exit_pressure_pa=exit_pressure_pa,
        ambient_pressure_pa=ambient_pressure_pa,
        momentum_thrust_n=momentum,
        pressure_thrust_n=pressure,
        total_thrust_n=total,
        effective_exhaust_velocity_m_s=effective_velocity,
        specific_impulse_s=effective_velocity / STANDARD_GRAVITY_M_S2,
        total_impulse_n_s=total * burn_time_s if burn_time_s is not None else None,
    )

