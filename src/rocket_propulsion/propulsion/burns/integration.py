"""Closed-form standalone solver for constant propulsion burns."""

from __future__ import annotations

from dataclasses import replace
from math import exp, log

from rocket_propulsion.core.errors import DomainError
from rocket_propulsion.propulsion.performance import STANDARD_GRAVITY_M_S2

from .inventory import tank_duration_limit, tank_masses_at, validate_tank_coverage
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
from .serialization import burn_input_hash, burn_result_hash


def _target_duration(
    definition: BurnDefinition, operating_point: EngineOperatingPoint
) -> float:
    target = definition.target
    mass_flow = operating_point.total_tank_flow_kg_s
    axial_thrust = operating_point.delivered_thrust_n * definition.cant_efficiency
    if target.kind is BurnTargetKind.DURATION:
        return target.value
    if target.kind is BurnTargetKind.PROPELLANT_MASS:
        return target.value / mass_flow
    if target.kind is BurnTargetKind.TOTAL_IMPULSE:
        return target.value / operating_point.delivered_thrust_n
    if target.kind is BurnTargetKind.IDEAL_DELTA_V:
        exponent = -target.value * mass_flow / axial_thrust
        return definition.initial_mass_kg * (1.0 - exp(exponent)) / mass_flow
    raise DomainError(f"Unsupported burn target: {target.kind}.")


def simulate_constant_burn(
    definition: BurnDefinition,
    operating_point: EngineOperatingPoint,
) -> BurnResult:
    """Evaluate an L0 constant burn with exact target and cutoff solutions.

    The method integrates only propulsion inventory, total impulse, and the
    ideal one-dimensional velocity increment ``integral(F_axial / m dt)``.
    It does not model gravity, drag, steering, attitude, or orbital motion.
    """

    validate_tank_coverage(definition.tanks, operating_point)
    requested_duration = _target_duration(definition, operating_point)
    mass_flow = operating_point.total_tank_flow_kg_s
    axial_thrust = operating_point.delivered_thrust_n * definition.cant_efficiency

    tank_limit, binding_tank = tank_duration_limit(definition.tanks, operating_point)
    dry_mass_limit = (definition.initial_mass_kg - definition.protected_dry_mass_kg) / mass_flow
    limits = (
        (requested_duration, CutoffReason.TARGET_REACHED, None),
        (tank_limit, CutoffReason.TANK_RESERVE, binding_tank),
        (dry_mass_limit, CutoffReason.DRY_MASS_FLOOR, "protected_dry_mass_kg"),
        (definition.maximum_duration_s, CutoffReason.MAXIMUM_DURATION, "maximum_duration_s"),
    )
    duration, cutoff_reason, binding_constraint = min(limits, key=lambda item: item[0])
    target_tolerance = max(1e-12, requested_duration * 1e-12)
    target_achieved = duration >= requested_duration - target_tolerance
    if target_achieved:
        duration = requested_duration
        cutoff_reason = CutoffReason.TARGET_REACHED
        binding_constraint = None

    final_mass = definition.initial_mass_kg - mass_flow * duration
    if final_mass <= 0.0:
        raise DomainError("Burn would produce a non-positive vehicle mass.")
    initial_tanks = tuple((tank.tank_id, tank.loaded_mass_kg) for tank in definition.tanks)
    final_tanks = tank_masses_at(definition.tanks, operating_point, duration)
    consumed_by_tank = tuple(
        (tank_id, initial_mass - dict(final_tanks)[tank_id])
        for tank_id, initial_mass in initial_tanks
    )
    consumed = mass_flow * duration
    delivered_impulse = operating_point.delivered_thrust_n * duration
    axial_impulse = axial_thrust * duration
    ideal_delta_v = axial_thrust / mass_flow * log(definition.initial_mass_kg / final_mass)
    mass_ratio_log = log(definition.initial_mass_kg / final_mass)
    rocket_equivalent_isp = ideal_delta_v / (STANDARD_GRAVITY_M_S2 * mass_ratio_log)

    initial_state = PropulsionState(0.0, definition.initial_mass_kg, initial_tanks)
    final_state = PropulsionState(duration, final_mass, final_tanks)
    profile = (
        BurnProfilePoint(
            time_s=0.0,
            delivered_thrust_n=operating_point.delivered_thrust_n,
            axial_thrust_n=axial_thrust,
            system_specific_impulse_s=operating_point.system_specific_impulse_s,
            total_mass_flow_kg_s=mass_flow,
            vehicle_mass_kg=definition.initial_mass_kg,
            tank_masses_kg=initial_tanks,
            phase="steady",
        ),
        BurnProfilePoint(
            time_s=duration,
            delivered_thrust_n=operating_point.delivered_thrust_n,
            axial_thrust_n=axial_thrust,
            system_specific_impulse_s=operating_point.system_specific_impulse_s,
            total_mass_flow_kg_s=mass_flow,
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
        average_delivered_thrust_n=operating_point.delivered_thrust_n,
        system_equivalent_specific_impulse_s=(
            delivered_impulse / (STANDARD_GRAVITY_M_S2 * consumed)
        ),
        ideal_delta_v_m_s=ideal_delta_v,
        rocket_equivalent_specific_impulse_s=rocket_equivalent_isp,
        thrust_centroid_time_s=duration / 2.0,
        mass_closure_error_kg=abs(
            consumed - sum(amount for _, amount in consumed_by_tank)
        ),
        impulse_closure_error_n_s=abs(
            delivered_impulse - operating_point.delivered_thrust_n * duration
        ),
        delta_v_closure_error_m_s=abs(
            ideal_delta_v
            - operating_point.system_specific_impulse_s
            * STANDARD_GRAVITY_M_S2
            * definition.cant_efficiency
            * mass_ratio_log
        ),
        exact_constant_sidera_export_available=True,
    )
    evidence = IntegrationEvidence(
        method="analytic_constant_v1",
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

