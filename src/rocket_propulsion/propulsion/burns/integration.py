"""Closed-form solver for constant-performance propulsion burns.

With constant delivered thrust ``F`` and total tank flow ``mdot`` the vehicle
mass is ``m(t) = m0 - mdot t``, the delivered impulse ``I = F t`` and the
ideal velocity increment ``dv = (F_ax / mdot) ln(m0 / m(t))`` [1]_, [2]_.
Inventory ledgers use exact rational arithmetic; the exponential and
logarithm are evaluated in a local 1100-digit decimal context [3]_ before the
single rounding to binary64.

References
----------
.. [1] NASA Glenn Research Center, "Ideal Rocket Equation".
       https://www1.grc.nasa.gov/beginners-guide-to-aeronautics/ideal-rocket-equation/
.. [2] NASA Glenn Research Center, "Specific Impulse".
       https://www1.grc.nasa.gov/beginners-guide-to-aeronautics/specific-impulse/
.. [3] PYTHON-3.12-DECIMAL. Python 3.12, "decimal" module.
       https://docs.python.org/3.12/library/decimal.html
"""

from __future__ import annotations

from dataclasses import replace
from decimal import Context, Decimal, localcontext
from fractions import Fraction

from rocket_propulsion.core.errors import DomainError
from rocket_propulsion.propulsion.performance import STANDARD_GRAVITY_M_S2

from .inventory import validate_tank_coverage
from .models import (
    BurnDefinition,
    BurnProfilePoint,
    BurnResult,
    BurnSummary,
    BurnTargetKind,
    CutoffReason,
    EngineOperatingPoint,
    IntegrationEvidence,
    PropulsionState,
)
from .numeric_outputs import _finite_rational, _fraction
from .serialization import burn_input_hash, burn_result_hash


def _negative_exponential(argument: Fraction) -> Fraction:
    """Return ``exp(-argument)`` as an exact fraction of a 1100-digit value.

    Parameters
    ----------
    argument : Fraction
        Non-negative exponent ``x``.

    Returns
    -------
    Fraction
        ``exp(-x)`` evaluated in a local 1100-digit decimal context [1]_, or
        exactly zero for ``x > 4096``.

    Notes
    -----
    1100 digits cover the exponent range of binary64 inputs, including the
    three-factor product formed when the rocket equation is inverted, so a
    tiny positive exponent is not lost.  Above 4096 the surviving mass
    fraction cannot produce a non-zero binary64 final mass even from the
    largest initial mass; such a target is unreachable, and the burn ends at
    whichever finite cutoff limit binds first.

    References
    ----------
    .. [1] PYTHON-3.12-DECIMAL. Python 3.12, "decimal" module.
           https://docs.python.org/3.12/library/decimal.html
    """
    if argument > 4096:
        return Fraction(0)
    with localcontext(Context(prec=1100)):
        value = Decimal(argument.numerator) / Decimal(argument.denominator)
        return Fraction((-value).exp())


def _mass_log(initial: Fraction, final: Fraction) -> Fraction:
    """Return ``ln(initial / final)`` from the exact mass ratio.

    The logarithm of the rocket-equation mass ratio [1]_ is evaluated in a
    local 1100-digit decimal context [2]_ before any rounding, so the ratio is
    never rounded to binary64 first.

    References
    ----------
    .. [1] NASA Glenn Research Center, "Ideal Rocket Equation".
           https://www1.grc.nasa.gov/beginners-guide-to-aeronautics/ideal-rocket-equation/
    .. [2] PYTHON-3.12-DECIMAL. Python 3.12, "decimal" module.
           https://docs.python.org/3.12/library/decimal.html
    """
    with localcontext(Context(prec=1100)):
        ratio = initial / final
        return Fraction((Decimal(ratio.numerator) / Decimal(ratio.denominator)).ln())


def _target_duration(
    definition: BurnDefinition, operating_point: EngineOperatingPoint
) -> Fraction:
    """Return the exact burn duration that meets the requested target.

    ``DURATION`` gives the value itself; ``PROPELLANT_MASS`` gives
    ``m_p / mdot``; ``TOTAL_IMPULSE`` gives ``I / F`` with the delivered
    thrust [2]_; ``IDEAL_DELTA_V`` inverts the rocket equation [1]_ with the
    axial thrust ``F_ax = F * cant_efficiency``,
    ``t = m0 (1 - exp(-dv mdot / F_ax)) / mdot``.

    Raises
    ------
    DomainError
        If the target kind is not supported.

    References
    ----------
    .. [1] NASA Glenn Research Center, "Ideal Rocket Equation".
           https://www1.grc.nasa.gov/beginners-guide-to-aeronautics/ideal-rocket-equation/
    .. [2] NASA Glenn Research Center, "Specific Impulse".
           https://www1.grc.nasa.gov/beginners-guide-to-aeronautics/specific-impulse/
    """
    target = definition.target
    mass_flow = _fraction(operating_point.total_tank_flow_kg_s)
    axial_thrust = _fraction(operating_point.delivered_thrust_n) * _fraction(definition.cant_efficiency)
    value = _fraction(target.value)
    if target.kind is BurnTargetKind.DURATION:
        return value
    if target.kind is BurnTargetKind.PROPELLANT_MASS:
        return value / mass_flow
    if target.kind is BurnTargetKind.TOTAL_IMPULSE:
        return value / _fraction(operating_point.delivered_thrust_n)
    if target.kind is BurnTargetKind.IDEAL_DELTA_V:
        exponent = value * mass_flow / axial_thrust
        return _fraction(definition.initial_mass_kg) * (1 - _negative_exponential(exponent)) / mass_flow
    raise DomainError(f"Unsupported burn target: {target.kind}.")


def simulate_constant_burn(
    definition: BurnDefinition,
    operating_point: EngineOperatingPoint,
) -> BurnResult:
    """Evaluate an L0 constant-performance burn in closed form.

    Parameters
    ----------
    definition : BurnDefinition
        Initial mass, tanks, target, cant efficiency and cutoff limits.
    operating_point : EngineOperatingPoint
        Constant delivered thrust, tank streams and specific impulse.

    Returns
    -------
    BurnResult
        Initial and final states, a two-point profile, summary, integration
        evidence, warnings and input/result hashes.

    Raises
    ------
    DomainError
        If the tank coverage is invalid, the protected dry mass exceeds the
        initial mass, the target kind is unsupported, the burn would leave a
        non-positive vehicle mass, or an exported value is not representable
        as a finite binary64 value.

    Notes
    -----
    The burn ends at the earliest of the target duration, each tank's reserve
    limit, the protected dry-mass floor and the maximum duration; ties resolve
    in that order.  The method integrates only propulsion inventory, total
    impulse ``I = F t`` [2]_ and the ideal one-dimensional velocity increment
    ``integral(F_ax / m dt) = (F_ax / mdot) ln(m0 / m1)`` [1]_.  It does not
    model gravity, drag, steering, attitude or orbital motion.  Equivalent
    specific impulses are ``F / (g0 mdot)`` and ``F_ax / (g0 mdot)`` with
    ``g0`` the registered ``STANDARD_GRAVITY_M_S2``.  Durations are inverted
    with a bounded-precision (1100-digit) exponential [3]_; the reported
    delta-v closure error compares the velocity increment with
    ``Isp g0 cant_efficiency ln(m0 / m1)``.

    References
    ----------
    .. [1] NASA Glenn Research Center, "Ideal Rocket Equation".
           https://www1.grc.nasa.gov/beginners-guide-to-aeronautics/ideal-rocket-equation/
    .. [2] NASA Glenn Research Center, "Specific Impulse".
           https://www1.grc.nasa.gov/beginners-guide-to-aeronautics/specific-impulse/
    .. [3] PYTHON-3.12-DECIMAL. Python 3.12, "decimal" module.
           https://docs.python.org/3.12/library/decimal.html
    """

    validate_tank_coverage(definition.tanks, operating_point)
    requested_duration = _target_duration(definition, operating_point)
    mass_flow = _fraction(operating_point.total_tank_flow_kg_s)
    force = _fraction(operating_point.delivered_thrust_n)
    axial = force * _fraction(definition.cant_efficiency)
    axial_thrust = _finite_rational(axial, "Axial thrust")
    initial_mass = _fraction(definition.initial_mass_kg)
    rates = {}
    for stream in operating_point.streams:
        if stream.tank_depleting:
            rates[stream.tank_id] = rates.get(stream.tank_id, Fraction(0)) + _fraction(stream.mass_flow_kg_s)
    tank_limits = tuple(
        ((_fraction(tank.loaded_mass_kg) - _fraction(tank.reserve_mass_kg)) / rates[tank.tank_id],
         CutoffReason.TANK_RESERVE, tank.tank_id)
        for tank in definition.tanks if rates.get(tank.tank_id, 0) > 0)
    dry_mass_limit = (initial_mass - _fraction(definition.protected_dry_mass_kg)) / mass_flow
    if dry_mass_limit < 0:
        raise DomainError("Protected dry mass exceeds initial vehicle mass.", field="Protected dry mass")
    limits = (
        (requested_duration, CutoffReason.TARGET_REACHED, None),
        *tank_limits,
        (dry_mass_limit, CutoffReason.DRY_MASS_FLOOR, "protected_dry_mass_kg"),
        (_fraction(definition.maximum_duration_s), CutoffReason.MAXIMUM_DURATION, "maximum_duration_s"),
    )
    exact_duration, cutoff_reason, binding_constraint = min(limits, key=lambda item: item[0])
    target_achieved = cutoff_reason is CutoffReason.TARGET_REACHED
    duration = _finite_rational(exact_duration, "Burn duration")
    consumption = mass_flow * exact_duration
    consumed = _finite_rational(consumption, "Consumed propellant")
    remaining = initial_mass - consumption
    if remaining <= 0:
        raise DomainError("Burn would produce a non-positive vehicle mass.", field="Final vehicle mass")
    final_mass = _finite_rational(remaining, "Final vehicle mass")
    initial_tanks = tuple((tank.tank_id, tank.loaded_mass_kg) for tank in definition.tanks)
    tank_consumption = tuple((tank.tank_id, rates.get(tank.tank_id, Fraction(0)) * exact_duration)
                             for tank in definition.tanks)
    final_tanks = tuple((tank_id, _finite_rational(_fraction(loaded) - used, f"Tank {tank_id!r} mass"))
                        for (tank_id, loaded), (_, used) in zip(initial_tanks, tank_consumption))
    consumed_by_tank = tuple((name, _finite_rational(used, f"Tank {name!r} consumption"))
                             for name, used in tank_consumption)
    delivered_impulse = _finite_rational(force * exact_duration, "Delivered impulse")
    axial_impulse = _finite_rational(axial * exact_duration, "Axial impulse")
    mass_ratio_log = _mass_log(initial_mass, remaining)
    velocity = axial / mass_flow * mass_ratio_log
    ideal_delta_v = _finite_rational(velocity, "Ideal delta-v")
    gravity = _fraction(STANDARD_GRAVITY_M_S2)
    system_isp = _finite_rational(force / (gravity * mass_flow), "System equivalent specific impulse")
    rocket_equivalent_isp = _finite_rational(axial / (gravity * mass_flow), "Rocket equivalent specific impulse")

    initial_state = PropulsionState(0.0, definition.initial_mass_kg, initial_tanks)
    final_state = PropulsionState(duration, final_mass, final_tanks)
    profile = (
        BurnProfilePoint(
            time_s=0.0,
            delivered_thrust_n=operating_point.delivered_thrust_n,
            axial_thrust_n=axial_thrust,
            system_specific_impulse_s=operating_point.system_specific_impulse_s,
            total_mass_flow_kg_s=operating_point.total_tank_flow_kg_s,
            vehicle_mass_kg=definition.initial_mass_kg,
            tank_masses_kg=initial_tanks,
            phase="steady",
        ),
        BurnProfilePoint(
            time_s=duration,
            delivered_thrust_n=operating_point.delivered_thrust_n,
            axial_thrust_n=axial_thrust,
            system_specific_impulse_s=operating_point.system_specific_impulse_s,
            total_mass_flow_kg_s=operating_point.total_tank_flow_kg_s,
            vehicle_mass_kg=final_mass,
            tank_masses_kg=final_tanks,
            phase="cutoff",
        ),
    )
    summary = BurnSummary(
        requested_target_kind=definition.target.kind,
        requested_target_value=definition.target.value,
        target_achieved=target_achieved,
        cutoff_reason=cutoff_reason,
        binding_constraint=binding_constraint,
        duration_s=duration,
        initial_mass_kg=definition.initial_mass_kg,
        final_mass_kg=final_mass,
        consumed_propellant_kg=consumed,
        consumed_by_tank_kg=consumed_by_tank,
        delivered_total_impulse_n_s=delivered_impulse,
        axial_total_impulse_n_s=axial_impulse,
        average_delivered_thrust_n=operating_point.delivered_thrust_n if exact_duration else 0.,
        system_equivalent_specific_impulse_s=system_isp,
        ideal_delta_v_m_s=ideal_delta_v,
        rocket_equivalent_specific_impulse_s=rocket_equivalent_isp,
        thrust_centroid_time_s=_finite_rational(exact_duration / 2, "Thrust centroid time"),
        mass_closure_error_kg=_finite_rational(abs(consumption - sum(
            (used for _, used in tank_consumption), Fraction(0))), "Mass closure error"),
        impulse_closure_error_n_s=0.,
        delta_v_closure_error_m_s=_finite_rational(abs(velocity - _fraction(
            operating_point.system_specific_impulse_s) * gravity * _fraction(
                definition.cant_efficiency) * mass_ratio_log), "Delta-v closure error"),
        exact_constant_sidera_export_available=True,
    )
    evidence = IntegrationEvidence(
        method="analytic_constant_binary_ledgers_decimal1100_v2",
        absolute_tolerance=1e-12,
        relative_tolerance=1e-12,
        segment_count=1,
        accepted_steps=1,
        rejected_steps=0,
    )
    provisional = BurnResult(
        schema="rocket_propulsion_burn_v1",
        definition=definition,
        operating_point=operating_point,
        initial_state=initial_state,
        final_state=final_state,
        profile=profile,
        summary=summary,
        evidence=evidence,
        warnings=operating_point.warnings,
        input_hash=burn_input_hash(definition, operating_point),
        result_hash="",
    )
    return replace(provisional, result_hash=burn_result_hash(provisional))

