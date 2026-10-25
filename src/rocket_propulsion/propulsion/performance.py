"""Ideal rocket-equation and steady thrust-equation models with replayable results.

Both models are evaluated on an exact rational ledger of the binary inputs and every
exported quantity is rounded to binary64 once and independently. Result records
re-derive their fields from the retained inputs on construction. Numerical
consistency does not validate either model against hardware.

Notes
-----
Standard gravity is ``rocket_propulsion.core.constants.STANDARD_GRAVITY_M_S2``.

References
----------
.. [1] NASA Glenn Research Center, "Ideal Rocket Equation", Beginner's Guide to
       Aeronautics (Tsiolkovsky equation, mass ratio, effective exhaust velocity).
.. [2] NASA Glenn Research Center, "Rocket Thrust Equation", Beginner's Guide to
       Aeronautics (momentum thrust, pressure thrust, exit-area effects).
.. [3] NASA Glenn Research Center, "Specific Impulse", Beginner's Guide to
       Aeronautics (equivalent velocity, total impulse and Isp).
.. [4] NASA-STEADY-ROCKET-THRUST. NASA Glenn Research Center, momentum and pressure
       thrust equation for steady rocket flow.
.. [5] SUTTON-BIBLARZ-2017. G. P. Sutton and O. Biblarz, "Rocket Propulsion
       Elements", 9th ed., Wiley, 2017.
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
    """Result of an ideal rocket-equation solve, replayable from its retained inputs.

    Every numeric field is an independent export of the exact solve. On construction
    the record solves the equation again from ``initial_mass_kg``, the retained
    performance input and either ``final_mass_kg`` (forward mode) or ``delta_v_m_s``
    (inverse mode), and refuses any field that the solve does not reproduce.

    Attributes
    ----------
    initial_mass_kg : float
        Initial vehicle mass m0, kg.
    final_mass_kg : float
        Final vehicle mass mf, kg.
    propellant_mass_kg : float
        Consumed propellant mass m0 - mf, kg.
    mass_ratio : float
        Mass ratio m0 / mf.
    propellant_mass_fraction : float
        Consumed fraction (m0 - mf) / m0.
    specific_impulse_s : float
        Specific impulse, s.
    effective_exhaust_velocity_m_s : float
        Effective exhaust velocity c = g0 * Isp, m/s.
    delta_v_m_s : float
        Ideal velocity change, m/s.
    solve_mode : {'forward', 'inverse'}
        ``'forward'`` solves delta-v from the final mass; ``'inverse'`` solves the
        masses from delta-v.
    performance_input : {'specific_impulse_s', 'effective_exhaust_velocity_m_s'}
        Field supplied by the caller; the other performance field is derived.
    numerical_policy : str
        Identifier of the arithmetic and rounding policy.

    Raises
    ------
    DomainError
        If the policy, solve mode or performance input is unknown, or if a field does
        not reproduce the retained inputs.

    Notes
    -----
    The governing relation is the ideal rocket equation delta_v = c * ln(m0 / mf)
    [1]_ [2]_ with c = g0 * Isp [3]_.

    References
    ----------
    .. [1] NASA Glenn Research Center, "Ideal Rocket Equation", Beginner's Guide to
           Aeronautics (Tsiolkovsky equation, mass ratio, effective exhaust velocity).
    .. [2] SUTTON-BIBLARZ-2017. G. P. Sutton and O. Biblarz, "Rocket Propulsion
           Elements", 9th ed., Wiley, 2017.
    .. [3] NASA Glenn Research Center, "Specific Impulse", Beginner's Guide to
           Aeronautics (equivalent velocity, total impulse and Isp).
    """
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
    """Steady momentum/pressure thrust ledger, replayable from its retained inputs.

    On construction the record evaluates the thrust equation again from the retained
    flow, exit state, ambient pressure and optional burn time, and refuses any field
    that the evaluation does not reproduce.

    Attributes
    ----------
    mass_flow_kg_s : float
        Nozzle mass flow, kg/s.
    exit_velocity_m_s : float
        Exit-plane velocity, m/s.
    exit_area_m2 : float
        Nozzle exit area, m^2.
    exit_pressure_pa : float
        Exit-plane static pressure, Pa.
    ambient_pressure_pa : float
        Ambient pressure, Pa.
    momentum_thrust_n : float
        Momentum thrust mdot * Ve, N.
    pressure_thrust_n : float
        Signed pressure thrust Ae * (pe - pa), N.
    total_thrust_n : float
        Net thrust, N; always positive.
    effective_exhaust_velocity_m_s : float
        Effective exhaust velocity F / mdot, m/s.
    specific_impulse_s : float
        Specific impulse F / (mdot * g0), s.
    total_impulse_n_s : float or None
        Constant-thrust impulse F * t when a burn time is retained, N*s.
    burn_time_s : float or None
        Retained burn duration, s.
    numerical_policy : str
        Identifier of the arithmetic and rounding policy.

    Raises
    ------
    DomainError
        If the policy is unknown or a field does not reproduce the retained inputs.

    Notes
    -----
    F = mdot * Ve + Ae * (pe - pa) [1]_ [2]_; Isp = F / (mdot * g0) [3]_.

    References
    ----------
    .. [1] NASA-STEADY-ROCKET-THRUST. NASA Glenn Research Center, momentum and pressure
           thrust equation for steady rocket flow.
    .. [2] NASA Glenn Research Center, "Rocket Thrust Equation", Beginner's Guide to
           Aeronautics (momentum thrust, pressure thrust, exit-area effects).
    .. [3] NASA Glenn Research Center, "Specific Impulse", Beginner's Guide to
           Aeronautics (equivalent velocity, total impulse and Isp).
    """
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
    """Return ``value`` as a finite real, refusing zero and negative values.

    Parameters
    ----------
    value : object
        Candidate real number.
    label : str
        Field label used in the error.

    Returns
    -------
    float
        The validated value.

    Raises
    ------
    DomainError
        If ``value`` is not a finite real or is not greater than zero.
    """
    value = real(value, label)
    if value <= 0:
        raise DomainError(f'{label} must be greater than zero.', field=label)
    return value


def _performance_values(value, kind):
    """Convert the supplied Isp or exhaust velocity into the full performance triple.

    Parameters
    ----------
    value : float
        Specific impulse in s or effective exhaust velocity in m/s.
    kind : {'specific_impulse_s', 'effective_exhaust_velocity_m_s'}
        Which quantity ``value`` is.

    Returns
    -------
    tuple
        ``(isp_s, velocity_m_s, exact_velocity)``: the exported Isp and velocity
        and the exact rational velocity used by the solve.

    Notes
    -----
    c = g0 * Isp [1]_, with g0 the registered conventional standard gravity.

    References
    ----------
    .. [1] NASA Glenn Research Center, "Specific Impulse", Beginner's Guide to
           Aeronautics (equivalent velocity, total impulse and Isp).
    """
    gravity = Fraction(STANDARD_GRAVITY_M_S2)
    if kind == 'specific_impulse_s':
        isp = _positive(value, 'Specific impulse')
        velocity = Fraction(isp)*gravity
        return isp, _finite_rational(velocity, 'Effective exhaust velocity'), velocity
    velocity = _positive(value, 'Effective exhaust velocity')
    return _finite_rational(Fraction(velocity)/gravity, 'Specific impulse'), velocity, Fraction(velocity)


def _rocket_values(initial, final, delta, performance, performance_input, mode):
    """Solve the ideal rocket equation forward or inverse and export every field once.

    Parameters
    ----------
    initial, final, delta : float or None
        Initial mass (kg), final mass (kg, forward mode) and delta-v (m/s, inverse
        mode).
    performance : float
        Specific impulse (s) or effective exhaust velocity (m/s).
    performance_input : str
        Which performance quantity was supplied.
    mode : {'forward', 'inverse'}
        Solve direction.

    Returns
    -------
    dict
        Field values for :class:`RocketEquationResult`.

    Raises
    ------
    DomainError
        If an input is out of domain or an output is not representable.
    ConvergenceError
        If an enclosed log/exp output cannot be rounded within the precision budget.

    Notes
    -----
    Forward: delta_v = c * ln(m0 / mf). Inverse: mf = m0 * exp(-delta_v / c) [1]_.

    References
    ----------
    .. [1] NASA Glenn Research Center, "Ideal Rocket Equation", Beginner's Guide to
           Aeronautics (Tsiolkovsky equation, mass ratio, effective exhaust velocity).
    """
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
    return {'initial_mass_kg': initial, 'final_mass_kg': final, 'propellant_mass_kg': consumed,
                'mass_ratio': ratio, 'propellant_mass_fraction': fraction, 'specific_impulse_s': isp,
                'effective_exhaust_velocity_m_s': velocity, 'delta_v_m_s': delta}


def calculate_rocket_equation(
    *, initial_mass_kg: float | None = None, final_mass_kg: float | None = None,
    delta_v_m_s: float | None = None, specific_impulse_s: float | None = None,
    effective_exhaust_velocity_m_s: float | None = None,
) -> RocketEquationResult:
    """Solve the ideal rocket equation for delta-v or for the final mass.

    Supply ``initial_mass_kg``, exactly one of ``specific_impulse_s`` and
    ``effective_exhaust_velocity_m_s``, and exactly one of ``final_mass_kg`` (forward
    solve) and ``delta_v_m_s`` (inverse solve).

    Parameters
    ----------
    initial_mass_kg : float
        Initial vehicle mass, kg; positive.
    final_mass_kg : float, optional
        Final vehicle mass, kg; positive and smaller than the initial mass.
    delta_v_m_s : float, optional
        Required velocity change, m/s; positive.
    specific_impulse_s : float, optional
        Specific impulse, s; positive.
    effective_exhaust_velocity_m_s : float, optional
        Effective exhaust velocity, m/s; positive.

    Returns
    -------
    RocketEquationResult
        Masses, mass ratio, propellant fraction, performance and delta-v.

    Raises
    ------
    DomainError
        If the input combination is ambiguous, an input is out of domain, or a
        mandatory output is not representable as a finite binary64 value.
    ConvergenceError
        If an enclosed logarithm or exponential cannot be rounded within the
        precision budget.

    Notes
    -----
    delta_v = c * ln(m0 / mf) with constant effective exhaust velocity c = g0 * Isp
    [1]_ [2]_ [3]_. Gravity, drag, steering and finite-burn losses are omitted. In
    inverse mode the rounded initial and final masses may coincide while the
    independently exported propellant consumption remains positive.

    References
    ----------
    .. [1] NASA Glenn Research Center, "Ideal Rocket Equation", Beginner's Guide to
           Aeronautics (Tsiolkovsky equation, mass ratio, effective exhaust velocity).
    .. [2] NASA Glenn Research Center, "Specific Impulse", Beginner's Guide to
           Aeronautics (equivalent velocity, total impulse and Isp).
    .. [3] SUTTON-BIBLARZ-2017. G. P. Sutton and O. Biblarz, "Rocket Propulsion
           Elements", 9th ed., Wiley, 2017.
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
    """Evaluate the steady thrust equation on exact rationals and export each field once.

    Parameters
    ----------
    flow, velocity, area, exit_pressure : float
        Mass flow (kg/s), exit velocity (m/s), exit area (m^2) and exit pressure (Pa);
        all positive.
    ambient_pressure : float
        Ambient pressure, Pa; non-negative.
    duration : float or None
        Optional positive burn time, s.

    Returns
    -------
    dict
        Field values for :class:`ThrustResult`.

    Raises
    ------
    DomainError
        If an input is out of domain, the net thrust is not positive, or an output is
        not representable.

    Notes
    -----
    F = mdot * Ve + Ae * (pe - pa) [1]_; Isp = F / (mdot * g0) [2]_.

    References
    ----------
    .. [1] NASA-STEADY-ROCKET-THRUST. NASA Glenn Research Center, momentum and pressure
           thrust equation for steady rocket flow.
    .. [2] NASA Glenn Research Center, "Specific Impulse", Beginner's Guide to
           Aeronautics (equivalent velocity, total impulse and Isp).
    """
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
    return {'mass_flow_kg_s': flow, 'exit_velocity_m_s': velocity, 'exit_area_m2': area,
                'exit_pressure_pa': exit_pressure, 'ambient_pressure_pa': ambient_pressure,
                'momentum_thrust_n': _finite_rational(momentum, 'Momentum thrust'),
                'pressure_thrust_n': _finite_rational(pressure, 'Pressure thrust'),
                'total_thrust_n': _finite_rational(total, 'Total thrust'),
                'effective_exhaust_velocity_m_s': _finite_rational(total/flow_exact, 'Effective exhaust velocity'),
                'specific_impulse_s': _finite_rational(total/(flow_exact*Fraction(STANDARD_GRAVITY_M_S2)), 'Specific impulse'),
                'total_impulse_n_s': None if duration is None else _finite_rational(total*Fraction(duration), 'Total impulse'),
                'burn_time_s': duration}


def calculate_thrust(
    *, mass_flow_kg_s: float, exit_velocity_m_s: float, exit_area_m2: float,
    exit_pressure_pa: float, ambient_pressure_pa: float, burn_time_s: float | None = None,
) -> ThrustResult:
    """Evaluate steady rocket thrust from the exit-plane state.

    Parameters
    ----------
    mass_flow_kg_s : float
        Nozzle mass flow, kg/s; positive.
    exit_velocity_m_s : float
        Exit-plane velocity, m/s; positive.
    exit_area_m2 : float
        Nozzle exit area, m^2; positive.
    exit_pressure_pa : float
        Exit-plane static pressure, Pa; positive.
    ambient_pressure_pa : float
        Ambient pressure, Pa; non-negative.
    burn_time_s : float, optional
        Burn duration, s. When given, the constant-thrust total impulse is exported
        and the duration is retained so the impulse can be replayed.

    Returns
    -------
    ThrustResult
        Momentum, pressure and total thrust with effective exhaust velocity, Isp and
        optional total impulse.

    Raises
    ------
    DomainError
        If an input is out of domain, the net thrust is not strictly positive, or an
        output is not representable.

    Notes
    -----
    F = mdot * Ve + Ae * (pe - pa) [1]_ [2]_. The pressure term is signed, so an
    over-expanded nozzle reduces thrust; only a strictly positive net thrust is
    accepted. Isp = F / (mdot * g0) [3]_. Every derived output is rounded independently
    from the exact binary-input ledger.

    References
    ----------
    .. [1] NASA-STEADY-ROCKET-THRUST. NASA Glenn Research Center, momentum and pressure
           thrust equation for steady rocket flow.
    .. [2] NASA Glenn Research Center, "Rocket Thrust Equation", Beginner's Guide to
           Aeronautics (momentum thrust, pressure thrust, exit-area effects).
    .. [3] NASA Glenn Research Center, "Specific Impulse", Beginner's Guide to
           Aeronautics (equivalent velocity, total impulse and Isp).
    """
    return ThrustResult(**_thrust_values(mass_flow_kg_s, exit_velocity_m_s, exit_area_m2,
                                        exit_pressure_pa, ambient_pressure_pa, burn_time_s))
