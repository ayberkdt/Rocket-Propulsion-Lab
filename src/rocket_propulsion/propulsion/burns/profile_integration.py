"""Analytic propulsion-only integration of piecewise-linear L1 burn profiles."""

from __future__ import annotations

from dataclasses import replace
from math import exp, inf, log

from rocket_propulsion.core.errors import DomainError
from rocket_propulsion.propulsion.performance import STANDARD_GRAVITY_M_S2

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
from .profiles import (
    integrated_throttle,
    resolved_schedule_intervals,
    throttle_first_moment,
    time_for_integrated_throttle,
)
from .serialization import burn_input_hash, burn_result_hash


def _required_exposure(
    definition: BurnDefinition, operating_point: EngineOperatingPoint
) -> float | None:
    target = definition.target
    mass_flow = operating_point.total_tank_flow_kg_s
    axial_thrust = operating_point.delivered_thrust_n * definition.cant_efficiency
    if target.kind is BurnTargetKind.DURATION:
        return None
    if target.kind is BurnTargetKind.PROPELLANT_MASS:
        return target.value / mass_flow
    if target.kind is BurnTargetKind.TOTAL_IMPULSE:
        return target.value / operating_point.delivered_thrust_n
    if target.kind is BurnTargetKind.IDEAL_DELTA_V:
        final_mass = definition.initial_mass_kg * exp(-target.value * mass_flow / axial_thrust)
        return (definition.initial_mass_kg - final_mass) / mass_flow
    raise DomainError(f"Unsupported burn target: {target.kind}.")


def _tank_rates(operating_point: EngineOperatingPoint) -> dict[str, float]:
    rates: dict[str, float] = {}
    for stream in operating_point.streams:
        if stream.tank_depleting:
            rates[stream.tank_id] = rates.get(stream.tank_id, 0.0) + stream.mass_flow_kg_s
    return rates


def _tank_masses_at_exposure(
    definition: BurnDefinition,
    operating_point: EngineOperatingPoint,
    exposure_s: float,
) -> tuple[tuple[str, float], ...]:
    rates = _tank_rates(operating_point)
    return tuple(
        (tank.tank_id, tank.loaded_mass_kg - rates.get(tank.tank_id, 0.0) * exposure_s)
        for tank in definition.tanks
    )


def _point(
    *,
    time_s: float,
    throttle: float,
    phase: str,
    exposure_s: float,
    definition: BurnDefinition,
    operating_point: EngineOperatingPoint,
) -> ProfiledBurnPoint:
    return ProfiledBurnPoint(
        time_s=time_s,
        commanded_throttle=throttle,
        realized_throttle=throttle,
        delivered_thrust_n=operating_point.delivered_thrust_n * throttle,
        axial_thrust_n=(
            operating_point.delivered_thrust_n * definition.cant_efficiency * throttle
        ),
        system_specific_impulse_s=operating_point.system_specific_impulse_s,
        total_mass_flow_kg_s=operating_point.total_tank_flow_kg_s * throttle,
        oxidizer_mass_flow_kg_s=operating_point.oxidizer_mass_flow_kg_s * throttle,
        fuel_mass_flow_kg_s=operating_point.fuel_mass_flow_kg_s * throttle,
        vehicle_mass_kg=(
            definition.initial_mass_kg
            - operating_point.total_tank_flow_kg_s * exposure_s
        ),
        tank_masses_kg=_tank_masses_at_exposure(
            definition, operating_point, exposure_s
        ),
        chamber_pressure_pa=operating_point.chamber_pressure_pa * throttle,
        phase=phase,
    )


def _build_profile(
    definition: BurnDefinition,
    operating_point: EngineOperatingPoint,
    schedule: ThrottleSchedule,
    duration_s: float,
    cutoff_reason: CutoffReason,
) -> tuple[ProfiledBurnPoint, ...]:
    intervals = resolved_schedule_intervals(schedule, stop_time_s=duration_s)
    if not intervals:
        first = schedule.segments[0]
        return (
            _point(
                time_s=0.0,
                throttle=first.start_throttle,
                phase=f"{first.phase}:cutoff",
                exposure_s=0.0,
                definition=definition,
                operating_point=operating_point,
            ),
        )
    points: list[ProfiledBurnPoint] = []
    exposure = 0.0
    for interval in intervals:
        if (
            not points
            or points[-1].time_s != interval.start_time_s
            or points[-1].commanded_throttle != interval.start_throttle
        ):
            points.append(
                _point(
                    time_s=interval.start_time_s,
                    throttle=interval.start_throttle,
                    phase=interval.phase,
                    exposure_s=exposure,
                    definition=definition,
                    operating_point=operating_point,
                )
            )
        exposure += interval.exposure_s
        is_final = interval is intervals[-1]
        phase = (
            f"{interval.phase}:cutoff:{cutoff_reason.value}"
            if is_final
            else interval.phase
        )
        points.append(
            _point(
                time_s=interval.end_time_s,
                throttle=interval.end_throttle,
                phase=phase,
                exposure_s=exposure,
                definition=definition,
                operating_point=operating_point,
            )
        )
    return tuple(points)


def simulate_profiled_burn(
    definition: BurnDefinition,
    operating_point: EngineOperatingPoint,
    schedule: ThrottleSchedule,
) -> ProfiledBurnResult:
    """Integrate an L1 schedule exactly without trajectory-state propagation.

    Delivered thrust and every propellant stream scale linearly with throttle;
    system specific impulse therefore remains constant. Every segment integral
    and cutoff inversion is analytic.
    """

    validate_tank_coverage(definition.tanks, operating_point)
    requested_exposure = _required_exposure(definition, operating_point)
    if definition.target.kind is BurnTargetKind.DURATION:
        target_time = definition.target.value
    else:
        resolved_target = time_for_integrated_throttle(schedule, requested_exposure or 0.0)
        target_time = inf if resolved_target is None else resolved_target

    rates = _tank_rates(operating_point)
    tank_limits = tuple(
        (tank.expendable_mass_kg / rates[tank.tank_id], tank.tank_id)
        for tank in definition.tanks
        if rates.get(tank.tank_id, 0.0) > 0.0
    )
    tank_exposure, binding_tank = min(tank_limits, default=(inf, None))
    tank_time_value = time_for_integrated_throttle(schedule, tank_exposure)
    tank_time = inf if tank_time_value is None else tank_time_value
    dry_exposure = (
        definition.initial_mass_kg - definition.protected_dry_mass_kg
    ) / operating_point.total_tank_flow_kg_s
    dry_time_value = time_for_integrated_throttle(schedule, dry_exposure)
    dry_time = inf if dry_time_value is None else dry_time_value

    limits = (
        (target_time, CutoffReason.TARGET_REACHED, None),
        (tank_time, CutoffReason.TANK_RESERVE, binding_tank),
        (dry_time, CutoffReason.DRY_MASS_FLOOR, "protected_dry_mass_kg"),
        (
            definition.maximum_duration_s,
            CutoffReason.MAXIMUM_DURATION,
            "maximum_duration_s",
        ),
        (schedule.total_duration_s, CutoffReason.SCHEDULE_END, "throttle_schedule"),
    )
    duration, cutoff_reason, binding_constraint = min(limits, key=lambda item: item[0])
    target_tolerance = max(1e-12, abs(target_time) * 1e-12) if target_time != inf else 1e-12
    target_achieved = target_time != inf and duration >= target_time - target_tolerance
    if target_achieved:
        duration = target_time
        cutoff_reason = CutoffReason.TARGET_REACHED
        binding_constraint = None

    exposure = integrated_throttle(schedule, duration)
    consumed = operating_point.total_tank_flow_kg_s * exposure
    final_mass = definition.initial_mass_kg - consumed
    if final_mass <= 0.0:
        raise DomainError("Burn profile would produce a non-positive vehicle mass.")
    initial_tanks = tuple((tank.tank_id, tank.loaded_mass_kg) for tank in definition.tanks)
    final_tanks = _tank_masses_at_exposure(definition, operating_point, exposure)
    final_tank_map = dict(final_tanks)
    consumed_by_tank = tuple(
        (tank_id, initial_mass - final_tank_map[tank_id])
        for tank_id, initial_mass in initial_tanks
    )

    delivered_impulse = operating_point.delivered_thrust_n * exposure
    axial_impulse = delivered_impulse * definition.cant_efficiency
    mass_ratio_log = log(definition.initial_mass_kg / final_mass)
    ideal_delta_v = (
        operating_point.delivered_thrust_n
        * definition.cant_efficiency
        / operating_point.total_tank_flow_kg_s
        * mass_ratio_log
    )
    first_moment = throttle_first_moment(schedule, duration)
    profile = _build_profile(
        definition, operating_point, schedule, duration, cutoff_reason
    )
    resolved_intervals = resolved_schedule_intervals(schedule, stop_time_s=duration)
    interval_impulse = sum(
        operating_point.delivered_thrust_n * interval.exposure_s
        for interval in resolved_intervals
    )
    throttle_values = tuple(
        value
        for interval in resolved_intervals
        for value in (interval.start_throttle, interval.end_throttle)
    )
    exact_constant_export = bool(throttle_values) and (
        max(throttle_values) - min(throttle_values) <= 1e-12
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
        average_delivered_thrust_n=(
            delivered_impulse / duration if duration > 0.0 else 0.0
        ),
        system_equivalent_specific_impulse_s=operating_point.system_specific_impulse_s,
        ideal_delta_v_m_s=ideal_delta_v,
        rocket_equivalent_specific_impulse_s=(
            operating_point.system_specific_impulse_s * definition.cant_efficiency
        ),
        thrust_centroid_time_s=(first_moment / exposure if exposure > 0.0 else 0.0),
        mass_closure_error_kg=abs(
            consumed - sum(amount for _, amount in consumed_by_tank)
        ),
        impulse_closure_error_n_s=abs(delivered_impulse - interval_impulse),
        delta_v_closure_error_m_s=abs(
            ideal_delta_v
            - operating_point.system_specific_impulse_s
            * STANDARD_GRAVITY_M_S2
            * definition.cant_efficiency
            * mass_ratio_log
        ),
        exact_constant_sidera_export_available=exact_constant_export,
    )
    evidence = IntegrationEvidence(
        method="analytic_piecewise_linear_throttle_v1",
        absolute_tolerance=1e-12,
        relative_tolerance=1e-12,
        segment_count=len(resolved_intervals),
        accepted_steps=len(resolved_intervals),
        rejected_steps=0,
    )
    warnings = operating_point.warnings + (
        (
            "Piecewise-linear throttle scales thrust and all tank streams proportionally; "
            "system specific impulse is constant."
        ),
    )
    if not exact_constant_export:
        warnings += (
            (
                "The current Sidera contract cannot receive this transient exactly; exact "
                "export is disabled."
            ),
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
