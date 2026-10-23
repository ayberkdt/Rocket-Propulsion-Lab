"""Ideal rocket and steady thrust models with checked retained arithmetic.

Numerical consistency does not validate either model against hardware.
"""
from dataclasses import dataclass
from fractions import Fraction

from rocket_propulsion.core.constants import STANDARD_GRAVITY_M_S2
from rocket_propulsion.core.errors import DomainError
from rocket_propulsion.core.validation import coerce_fields, real
from .burns.numeric_outputs import _finite_rational
from .performance_numerics import _PERFORMANCE_POLICY, _inverse_mass_outputs, _scaled_log


@dataclass(frozen=True, slots=True)
class RocketEquationResult:
    """Independent exports of an ideal solve, retaining its source mode and input."""
    initial_mass_kg: float
    final_mass_kg: float
    propellant_mass_kg: float
    mass_ratio: float
    propellant_mass_fraction: float
    specific_impulse_s: float
    effective_exhaust_velocity_m_s: float
    delta_v_m_s: float
    solve_mode: str = 'forward'
    performance_input: str = 'specific_impulse_s'
    numerical_policy: str = _PERFORMANCE_POLICY

    def __post_init__(self):
        coerce_fields(self, 'Rocket equation result')
        if self.numerical_policy != _PERFORMANCE_POLICY:
            raise DomainError('Unknown performance numerical policy.', field='Numerical policy')
        if self.solve_mode not in ('forward', 'inverse'):
            raise DomainError('Unknown rocket equation solve mode.', field='Solve mode')
        if self.performance_input not in ('specific_impulse_s', 'effective_exhaust_velocity_m_s'):
            raise DomainError('Unknown rocket performance input.', field='Performance input')
        expected = _rocket_values(self.initial_mass_kg, self.final_mass_kg, self.delta_v_m_s,
                                  getattr(self, self.performance_input), self.performance_input,
                                  self.solve_mode)
        if any(getattr(self, name) != value for name, value in expected.items()):
            raise DomainError('Rocket outputs do not reproduce their retained source inputs.', field='Rocket outputs')


@dataclass(frozen=True, slots=True)
class ThrustResult:
    """Steady pressure/momentum ledger and duration for replaying impulse."""
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
    burn_time_s: float | None = None
    numerical_policy: str = _PERFORMANCE_POLICY

    def __post_init__(self):
        coerce_fields(self, 'Thrust result')
        if self.numerical_policy != _PERFORMANCE_POLICY:
            raise DomainError('Unknown performance numerical policy.', field='Numerical policy')
        expected = _thrust_values(self.mass_flow_kg_s, self.exit_velocity_m_s, self.exit_area_m2,
                                  self.exit_pressure_pa, self.ambient_pressure_pa, self.burn_time_s)
        if any(getattr(self, name) != value for name, value in expected.items()):
            raise DomainError('Thrust outputs do not reproduce their retained source inputs.', field='Thrust outputs')


def _positive(value, label):
    value = real(value, label)
    if value <= 0:
        raise DomainError(f'{label} must be greater than zero.', field=label)
    return value


def _performance_values(value, kind):
    gravity = Fraction(STANDARD_GRAVITY_M_S2)
    if kind == 'specific_impulse_s':
        isp = _positive(value, 'Specific impulse')
        velocity = Fraction(isp)*gravity
        return isp, _finite_rational(velocity, 'Effective exhaust velocity'), velocity
    velocity = _positive(value, 'Effective exhaust velocity')
    return _finite_rational(Fraction(velocity)/gravity, 'Specific impulse'), velocity, Fraction(velocity)


def _rocket_values(initial, final, delta, performance, performance_input, mode):
    initial = _positive(initial, 'Initial mass')
    isp, velocity, exact_velocity = _performance_values(performance, performance_input)
    initial_exact = Fraction(initial)
    if mode == 'forward':
        final = _positive(final, 'Final mass')
        if final >= initial:
            raise DomainError('Final mass must be smaller than initial mass.', field='Final mass')
        final_exact = Fraction(final)
        ratio_exact = initial_exact/final_exact
        ratio = _finite_rational(ratio_exact, 'Mass ratio')
        consumed_exact = initial_exact-final_exact
        consumed = _finite_rational(consumed_exact, 'Propellant mass')
        fraction = _finite_rational(consumed_exact/initial_exact, 'Propellant mass fraction')
        delta = _scaled_log(exact_velocity, ratio_exact, 'Delta-v')
    else:
        delta = _positive(delta, 'Delta-v')
        exponent = Fraction(delta)/exact_velocity
        ratio, final, consumed, fraction = _inverse_mass_outputs(initial_exact, exponent)
    return dict(initial_mass_kg=initial, final_mass_kg=final, propellant_mass_kg=consumed,
                mass_ratio=ratio, propellant_mass_fraction=fraction, specific_impulse_s=isp,
                effective_exhaust_velocity_m_s=velocity, delta_v_m_s=delta)


def calculate_rocket_equation(
    *, initial_mass_kg: float | None = None, final_mass_kg: float | None = None,
    delta_v_m_s: float | None = None, specific_impulse_s: float | None = None,
    effective_exhaust_velocity_m_s: float | None = None,
) -> RocketEquationResult:
    """Solve ideal delta-v or final mass, with constant exhaust velocity.

    Give exactly one of Isp/exhaust velocity, initial mass, and exactly one of
    final mass/delta-v. Gravity, drag, steering and finite-burn losses are omitted.
    Mandatory unrepresentable outputs fail typed; unresolved enclosures are
    distinct convergence failures. Inverse rounded masses may coincide while
    their independently exported positive propellant consumption remains nonzero.
    """
    if (specific_impulse_s is None) == (effective_exhaust_velocity_m_s is None):
        raise DomainError('Provide exactly one of specific impulse or effective exhaust velocity.',
                          field='Performance input')
    if (final_mass_kg is None) == (delta_v_m_s is None):
        raise DomainError('Provide exactly one of final mass or delta-v.', field='Solve mode')
    performance_input = ('specific_impulse_s' if specific_impulse_s is not None
                         else 'effective_exhaust_velocity_m_s')
    performance = specific_impulse_s if specific_impulse_s is not None else effective_exhaust_velocity_m_s
    mode = 'forward' if final_mass_kg is not None else 'inverse'
    values = _rocket_values(initial_mass_kg, final_mass_kg, delta_v_m_s, performance,
                            performance_input, mode)
    return RocketEquationResult(**values, solve_mode=mode, performance_input=performance_input)


def _thrust_values(flow, velocity, area, exit_pressure, ambient_pressure, duration):
    flow = _positive(flow, 'Mass flow')
    velocity = _positive(velocity, 'Exit velocity')
    area = _positive(area, 'Exit area')
    exit_pressure = _positive(exit_pressure, 'Exit pressure')
    ambient_pressure = real(ambient_pressure, 'Ambient pressure')
    if ambient_pressure < 0:
        raise DomainError('Ambient pressure must be non-negative.', field='Ambient pressure')
    duration = None if duration is None else _positive(duration, 'Burn time')
    flow_exact = Fraction(flow)
    momentum = flow_exact*Fraction(velocity)
    pressure = Fraction(area)*(Fraction(exit_pressure)-Fraction(ambient_pressure))
    total = momentum+pressure
    if total <= 0:
        raise DomainError('The specified state produces non-positive total thrust.', field='Total thrust')
    return dict(mass_flow_kg_s=flow, exit_velocity_m_s=velocity, exit_area_m2=area,
                exit_pressure_pa=exit_pressure, ambient_pressure_pa=ambient_pressure,
                momentum_thrust_n=_finite_rational(momentum, 'Momentum thrust'),
                pressure_thrust_n=_finite_rational(pressure, 'Pressure thrust'),
                total_thrust_n=_finite_rational(total, 'Total thrust'),
                effective_exhaust_velocity_m_s=_finite_rational(total/flow_exact, 'Effective exhaust velocity'),
                specific_impulse_s=_finite_rational(total/(flow_exact*Fraction(STANDARD_GRAVITY_M_S2)), 'Specific impulse'),
                total_impulse_n_s=None if duration is None else _finite_rational(total*Fraction(duration), 'Total impulse'),
                burn_time_s=duration)


def calculate_thrust(
    *, mass_flow_kg_s: float, exit_velocity_m_s: float, exit_area_m2: float,
    exit_pressure_pa: float, ambient_pressure_pa: float, burn_time_s: float | None = None,
) -> ThrustResult:
    """F=mdot*V+A*(Pe-Pa), with signed pressure and strictly positive net thrust.

    Every derived output is rounded independently from the exact binary-input
    ledger. A provided duration is retained to verify constant-thrust impulse.
    """
    return ThrustResult(**_thrust_values(mass_flow_kg_s, exit_velocity_m_s, exit_area_m2,
                                        exit_pressure_pa, ambient_pressure_pa, burn_time_s))
