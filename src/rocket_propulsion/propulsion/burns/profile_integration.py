"""Propulsion-only L1 integration of piecewise-linear throttle schedules.

Thrust and every tank stream scale with the throttle command, so consumed
propellant is ``mdot_rated E(t)`` and delivered impulse ``F_rated E(t)``,
where ``E`` is the exact throttle exposure [2]_.  The ideal velocity increment
follows the rocket equation [1]_.  Exposure and cutoff-event ledgers are kept
in exact rational arithmetic and rounded once on export.

References
----------
.. [1] NASA Glenn Research Center, "Ideal Rocket Equation".
       https://www1.grc.nasa.gov/beginners-guide-to-aeronautics/ideal-rocket-equation/
.. [2] NASA Glenn Research Center, "Specific Impulse".
       https://www1.grc.nasa.gov/beginners-guide-to-aeronautics/specific-impulse/
"""

from __future__ import annotations

from dataclasses import replace
from fractions import Fraction
from math import isfinite, nextafter

from rocket_propulsion.core.errors import DomainError
from rocket_propulsion.propulsion.performance import STANDARD_GRAVITY_M_S2

from .integration import _mass_log, _target_duration
from .inventory import validate_tank_coverage
from .models import (
    BurnDefinition,
    BurnSummary,
    BurnTargetKind,
    CutoffReason,
    EngineOperatingPoint,
    IntegrationEvidence,
    ProfiledBurnPoint,
    ProfiledBurnResult,
    PropulsionState,
    ThrottleSchedule,
)
from .numeric_outputs import _finite_rational, _fraction
from .profiles import _inverse_exposure, _schedule_integrals, _schedule_slices, _schedule_totals
from .serialization import burn_input_hash, burn_result_hash


def _required_exposure(definition, operating_point):
    """Return the exposure that meets a non-duration target, or ``None``.

    The constant-burn duration at full command equals the required exposure,
    because consumption, impulse and the rocket-equation mass all depend on
    the command only through ``E``.  Duration targets return ``None``; they
    are handled as time limits.
    """
    return (
        None
        if definition.target.kind is BurnTargetKind.DURATION
        else _target_duration(definition, operating_point)
    )


def _tank_rates(operating_point):
    """Return the exact rated tank-depleting flow per tank identifier."""
    rates: dict[str, Fraction] = {}
    for stream in operating_point.streams:
        if stream.tank_depleting:
            rates[stream.tank_id] = rates.get(stream.tank_id, Fraction(0)) + _fraction(
                stream.mass_flow_kg_s
            )
    return rates


def _tank_masses_at_exposure(definition, operating_point, exposure_s):
    """Return ``(tank_id, loaded - rate * exposure)`` for every tank, rounded once."""
    rates = _tank_rates(operating_point)
    exposure = _fraction(exposure_s)
    return tuple(
        (
            tank.tank_id,
            _finite_rational(
                _fraction(tank.loaded_mass_kg) - rates.get(tank.tank_id, Fraction(0)) * exposure,
                f"Tank {tank.tank_id!r} mass",
            ),
        )
        for tank in definition.tanks
    )


def _residual_upper(value, label):
    """Round the magnitude of a residual upward to a binary64 bound.

    The exact magnitude is converted to binary64; when the conversion rounds
    down, the next representable value toward infinity [1]_ is returned, so
    the exported value bounds the exact residual, including in the subnormal
    range.

    Raises
    ------
    DomainError
        If the bound is not representable as a finite binary64 value.

    References
    ----------
    .. [1] PYTHON-3.12-NEXTAFTER. Python 3.12, "math.nextafter".
           https://docs.python.org/3.12/library/math.html#math.nextafter
    """
    value = abs(value)
    try:
        result = float(value)
    except OverflowError as error:
        raise DomainError("Residual bound exceeds finite output range.", field=label) from error
    if not isfinite(result):
        raise DomainError("Residual bound exceeds finite output range.", field=label)
    if _fraction(result) < value:
        result = nextafter(result, float("inf"))
    if not isfinite(result):
        raise DomainError("Residual bound exceeds finite output range.", field=label)
    return result


def _point(*, time_s, throttle, phase, exposure_s, definition, operating_point):
    """Build one profile point at a command level and accumulated exposure.

    Thrust, axial thrust, total and per-propellant flows and chamber pressure
    scale linearly with the command; system specific impulse is constant
    [1]_.  Vehicle and tank masses follow from the exact exposure.

    References
    ----------
    .. [1] NASA Glenn Research Center, "Specific Impulse".
           https://www1.grc.nasa.gov/beginners-guide-to-aeronautics/specific-impulse/
    """
    throttle, exposure = _fraction(throttle), _fraction(exposure_s)
    force = _fraction(operating_point.delivered_thrust_n)
    return ProfiledBurnPoint(
        time_s=_finite_rational(_fraction(time_s), "Profile point time"),
        commanded_throttle=_finite_rational(throttle, "Profile command"),
        realized_throttle=_finite_rational(throttle, "Profile throttle"),
        delivered_thrust_n=_finite_rational(force * throttle, "Profile thrust"),
        axial_thrust_n=_finite_rational(
            force * _fraction(definition.cant_efficiency) * throttle, "Profile axial thrust"
        ),
        system_specific_impulse_s=operating_point.system_specific_impulse_s,
        total_mass_flow_kg_s=_finite_rational(
            _fraction(operating_point.total_tank_flow_kg_s) * throttle, "Profile mass flow"
        ),
        oxidizer_mass_flow_kg_s=_finite_rational(
            _fraction(operating_point.oxidizer_mass_flow_kg_s) * throttle, "Profile oxidizer flow"
        ),
        fuel_mass_flow_kg_s=_finite_rational(
            _fraction(operating_point.fuel_mass_flow_kg_s) * throttle, "Profile fuel flow"
        ),
        vehicle_mass_kg=_finite_rational(
            _fraction(definition.initial_mass_kg)
            - _fraction(operating_point.total_tank_flow_kg_s) * exposure,
            "Profile vehicle mass",
        ),
        tank_masses_kg=_tank_masses_at_exposure(definition, operating_point, exposure),
        chamber_pressure_pa=_finite_rational(
            _fraction(operating_point.chamber_pressure_pa) * throttle, "Profile chamber pressure"
        ),
        phase=phase,
    )


def _build_profile(
    definition, operating_point, schedule, duration_s, cutoff_reason, final_exposure
):
    """Build the profile points from the schedule slices up to cutoff.

    A point is emitted at each slice end, plus at a slice start where the
    command is discontinuous.  The last point carries the cutoff phase label
    and the exact cutoff exposure.  A zero-duration burn yields one point.
    """
    slices = _schedule_slices(schedule, duration_s)
    if not slices:
        return (
            _point(
                time_s=Fraction(0),
                throttle=_fraction(schedule.segments[0].start_throttle),
                phase=f"{schedule.segments[0].phase}:cutoff",
                exposure_s=Fraction(0),
                definition=definition,
                operating_point=operating_point,
            ),
        )
    points = []
    exposure = Fraction(0)
    previous_throttle = None
    for index, s in enumerate(slices):
        if previous_throttle != s.a:
            points.append(
                _point(
                    time_s=s.start,
                    throttle=s.a,
                    phase=s.phase,
                    exposure_s=exposure,
                    definition=definition,
                    operating_point=operating_point,
                )
            )
        exposure += s.exposure
        last = index == len(slices) - 1
        points.append(
            _point(
                time_s=s.start + s.duration,
                throttle=s.b,
                phase=f"{s.phase}:cutoff:{cutoff_reason.value}" if last else s.phase,
                exposure_s=final_exposure if last else exposure,
                definition=definition,
                operating_point=operating_point,
            )
        )
        previous_throttle = s.b
    return tuple(points)


def _profile_cutoff(definition, operating_point, schedule, rates):
    """Return ``(time, exposure, reason, binding_constraint)`` of the cutoff.

    Limits carry a priority rank: target 0, tank reserve 1, protected
    dry-mass floor 2, maximum duration 3, schedule end 4.  Time limits
    (duration target, maximum duration, schedule end) and exposure limits
    (non-duration target, tank reserves, dry-mass floor) are each reduced to
    their earliest member, ties going to the lower rank.  Exposures are
    compared exactly before any root is inverted; the earliest exposure limit
    is inverted to a time only if it is reached no later than the earliest
    time limit.  On an exact tie in both exposure and time the lower rank
    wins.

    Raises
    ------
    DomainError
        If the protected dry mass exceeds the initial mass.
    """
    total_time, _ = _schedule_totals(schedule)
    requested = _required_exposure(definition, operating_point)
    times: list[tuple[Fraction, int, CutoffReason, str | None]] = [
        (
            _fraction(definition.maximum_duration_s),
            3,
            CutoffReason.MAXIMUM_DURATION,
            "maximum_duration_s",
        ),
        (total_time, 4, CutoffReason.SCHEDULE_END, "throttle_schedule"),
    ]
    if requested is None:
        times.append((_fraction(definition.target.value), 0, CutoffReason.TARGET_REACHED, None))
    cap, time_rank, time_reason, time_binding = min(times, key=lambda row: (row[0], row[1]))
    cap_exposure = _schedule_integrals(schedule, cap)[0]
    exposures: list[tuple[Fraction, int, CutoffReason, str | None]] = []
    if requested is not None:
        exposures.append((requested, 0, CutoffReason.TARGET_REACHED, None))
    for tank in definition.tanks:
        rate = rates.get(tank.tank_id, 0)
        if rate > 0:
            exposures.append(
                (
                    (_fraction(tank.loaded_mass_kg) - _fraction(tank.reserve_mass_kg)) / rate,
                    1,
                    CutoffReason.TANK_RESERVE,
                    tank.tank_id,
                )
            )
    dry = _fraction(definition.initial_mass_kg) - _fraction(definition.protected_dry_mass_kg)
    if dry < 0:
        raise DomainError("Protected dry mass exceeds initial mass.", field="Protected dry mass")
    exposures.append(
        (
            dry / _fraction(operating_point.total_tank_flow_kg_s),
            2,
            CutoffReason.DRY_MASS_FLOOR,
            "protected_dry_mass_kg",
        )
    )
    exposure, rank, reason, binding = min(exposures, key=lambda row: (row[0], row[1]))
    if exposure > cap_exposure:
        return cap, cap_exposure, time_reason, time_binding
    time = _inverse_exposure(schedule, exposure)
    if exposure == cap_exposure and time == cap and time_rank < rank:
        return cap, cap_exposure, time_reason, time_binding
    return min(time, cap), exposure, reason, binding


def simulate_profiled_burn(
    definition: BurnDefinition, operating_point: EngineOperatingPoint, schedule: ThrottleSchedule
) -> ProfiledBurnResult:
    """Integrate an L1 burn under a piecewise-linear throttle schedule.

    Parameters
    ----------
    definition : BurnDefinition
        Initial mass, tanks, target, cant efficiency and cutoff limits.
    operating_point : EngineOperatingPoint
        Rated operating point at full command.
    schedule : ThrottleSchedule
        Piecewise-linear throttle command.

    Returns
    -------
    ProfiledBurnResult
        Initial and final states, profile points, summary, integration
        evidence, warnings and input/result hashes.

    Raises
    ------
    DomainError
        If the tank coverage is invalid, the protected dry mass exceeds the
        initial mass, the burn would leave a non-positive vehicle mass, or an
        exported value is not representable as a finite binary64 value.

    Notes
    -----
    Thrust and all streams scale in proportion to the command and system
    specific impulse is constant [2]_; this is a declared modelling
    assumption, and the exact numerical accounting does not validate it.  The
    command polynomials are integrated exactly; the inverse of the exposure
    and the rocket-equation logarithm [1]_ are evaluated in a local
    1100-digit decimal context [3]_.  Reported closure errors are upward
    rounded bounds on the exact residuals.

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
    rates = _tank_rates(operating_point)
    time, exposure, cutoff_reason, binding_constraint = _profile_cutoff(
        definition, operating_point, schedule, rates
    )
    duration = _finite_rational(time, "Burn duration")
    flow = _fraction(operating_point.total_tank_flow_kg_s)
    force = _fraction(operating_point.delivered_thrust_n)
    cant = _fraction(definition.cant_efficiency)
    consumed_exact = flow * exposure
    consumed = _finite_rational(consumed_exact, "Consumed propellant")
    remaining = _fraction(definition.initial_mass_kg) - consumed_exact
    if remaining <= 0:
        raise DomainError(
            "Burn profile would produce a non-positive vehicle mass.", field="Final vehicle mass"
        )
    final_mass = _finite_rational(remaining, "Final vehicle mass")
    initial_tanks = tuple((tank.tank_id, tank.loaded_mass_kg) for tank in definition.tanks)
    final_tanks = _tank_masses_at_exposure(definition, operating_point, exposure)
    tank_consumption = tuple(
        (tank.tank_id, rates.get(tank.tank_id, Fraction(0)) * exposure) for tank in definition.tanks
    )
    consumed_by_tank = tuple(
        (name, _finite_rational(amount, f"Tank {name!r} consumption"))
        for name, amount in tank_consumption
    )
    impulse = force * exposure
    delivered_impulse = _finite_rational(impulse, "Delivered impulse")
    axial_impulse = _finite_rational(impulse * cant, "Axial impulse")
    mass_log = _mass_log(_fraction(definition.initial_mass_kg), remaining)
    velocity = force * cant / flow * mass_log
    ideal_delta_v = _finite_rational(velocity, "Ideal delta-v")
    actual_exposure, first_moment = _schedule_integrals(schedule, time)
    displayed_exposure = _schedule_integrals(schedule, _fraction(duration))[0]
    slices = _schedule_slices(schedule, time)
    profile = _build_profile(definition, operating_point, schedule, time, cutoff_reason, exposure)
    throttles = tuple(value for s in slices for value in (s.a, s.b))
    exact_constant_export = bool(throttles) and min(throttles) == max(throttles)
    gravity = _fraction(STANDARD_GRAVITY_M_S2)
    summary = BurnSummary(
        requested_target_kind=definition.target.kind,
        requested_target_value=definition.target.value,
        target_achieved=cutoff_reason is CutoffReason.TARGET_REACHED,
        cutoff_reason=cutoff_reason,
        binding_constraint=binding_constraint,
        duration_s=duration,
        initial_mass_kg=definition.initial_mass_kg,
        final_mass_kg=final_mass,
        consumed_propellant_kg=consumed,
        consumed_by_tank_kg=consumed_by_tank,
        delivered_total_impulse_n_s=delivered_impulse,
        axial_total_impulse_n_s=axial_impulse,
        average_delivered_thrust_n=_finite_rational(
            impulse / time if time else Fraction(0), "Average thrust"
        ),
        system_equivalent_specific_impulse_s=_finite_rational(
            force / (gravity * flow), "System equivalent specific impulse"
        ),
        ideal_delta_v_m_s=ideal_delta_v,
        rocket_equivalent_specific_impulse_s=_finite_rational(
            force * cant / (gravity * flow), "Rocket equivalent specific impulse"
        ),
        thrust_centroid_time_s=_finite_rational(
            first_moment / exposure if exposure else Fraction(0), "Thrust centroid time"
        ),
        mass_closure_error_kg=_residual_upper(
            consumed_exact - sum((used for _, used in tank_consumption), Fraction(0)),
            "Mass closure bound",
        ),
        impulse_closure_error_n_s=_residual_upper(
            max(
                abs(force * (exposure - actual_exposure)),
                abs(force * (exposure - displayed_exposure)),
            ),
            "Impulse closure bound",
        ),
        delta_v_closure_error_m_s=_residual_upper(
            velocity
            - _fraction(operating_point.system_specific_impulse_s) * gravity * cant * mass_log,
            "Delta-v closure bound",
        ),
        exact_constant_sidera_export_available=exact_constant_export,
    )
    evidence = IntegrationEvidence(
        method="analytic_profile_binary_ledgers_decimal1100_v2",
        absolute_tolerance=1e-12,
        relative_tolerance=1e-12,
        segment_count=len(slices),
        accepted_steps=len(slices),
        rejected_steps=0,
    )
    warnings = operating_point.warnings + (
        "Piecewise-linear throttle scales thrust and all tank streams proportionally; system specific impulse is constant.",
    )
    if not exact_constant_export:
        warnings += (
            "The current Sidera contract cannot receive this transient exactly; exact export is disabled.",
        )
    provisional = ProfiledBurnResult(
        schema="rocket_propulsion_burn_v2",
        definition=definition,
        operating_point=operating_point,
        schedule=schedule,
        initial_state=PropulsionState(0.0, definition.initial_mass_kg, initial_tanks),
        final_state=PropulsionState(duration, final_mass, final_tanks),
        profile=profile,
        summary=summary,
        evidence=evidence,
        warnings=warnings,
        input_hash=burn_input_hash(definition, operating_point, schedule),
        result_hash="",
    )
    return replace(provisional, result_hash=burn_result_hash(provisional))
