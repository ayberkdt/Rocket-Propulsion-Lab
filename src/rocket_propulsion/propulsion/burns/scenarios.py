"""Discrete propulsion execution and failure scenarios for engine clusters.

Scenarios are mutually exclusive, explicitly sourced, and probability weighted.
The module never invents component failure rates or assumes independence. Each
scenario produces a complete target-neutral thrust/mass-flow artifact for later
trajectory propagation by Sidera.

References
----------
NASA-SP-2011-3421: https://ntrs.nasa.gov/citations/20120001369
NASA-20050207429: https://ntrs.nasa.gov/citations/20050207429
NASA-20160007073: https://ntrs.nasa.gov/citations/20160007073
"""

from __future__ import annotations

import hashlib
from dataclasses import dataclass, replace
from enum import Enum
from math import isfinite

from rocket_propulsion.core.errors import DomainError

from .blowdown import _scaled_operating_point
from .cluster import (
    ClusterSamplingEvidence,
    EngineClusterMember,
    EngineTransientCommand,
    build_engine_cluster_profile,
)
from .models import ThrottleSchedule, ThrottleSegment
from .profiles import resolved_schedule_intervals
from .serialization import canonical_json_bytes
from .tabulated import TabulatedBurnArtifact, TabulatedBurnSegment

SCENARIO_REFERENCE_IDS = (
    "NASA-SP-2011-3421",
    "NASA-20050207429",
    "NASA-20160007073",
    "NASA-TM-107318",
)


class ScenarioEventType(str, Enum):
    """Supported discrete propulsion execution modifications.

    References
    ----------
    NASA-20050207429: https://ntrs.nasa.gov/citations/20050207429
    NASA-20160007073: https://ntrs.nasa.gov/citations/20160007073
    """

    ENGINE_DISABLED = "engine-disabled"
    IGNITION_DELAY = "ignition-delay"
    EARLY_CUTOFF = "early-cutoff"
    THRUST_SCALE = "thrust-scale"
    SPECIFIC_IMPULSE_SCALE = "specific-impulse-scale"
    RESPONSE_TIME_SCALE = "response-time-scale"


@dataclass(frozen=True, slots=True)
class ScenarioEvent:
    """One sourced modification applied to one or more cluster engines.

    Multiple engine identifiers in one event express an explicitly declared
    dependent/common-cause outcome. The code does not infer that dependence
    from identical hardware.

    References
    ----------
    NASA-SP-2011-3421: https://ntrs.nasa.gov/citations/20120001369
    NASA-20160007073: https://ntrs.nasa.gov/citations/20160007073
    """

    event_type: ScenarioEventType
    engine_ids: tuple[str, ...]
    value: float | None
    source: str
    rationale: str

    def __post_init__(self) -> None:
        if not self.engine_ids or any(not engine_id.strip() for engine_id in self.engine_ids):
            raise DomainError("Scenario event requires non-empty engine identifiers.")
        if len(self.engine_ids) != len(set(self.engine_ids)):
            raise DomainError("Scenario event engine identifiers must be unique.")
        if not self.source.strip() or not self.rationale.strip():
            raise DomainError("Scenario event source and rationale are required.")
        if self.event_type is ScenarioEventType.ENGINE_DISABLED:
            if self.value is not None:
                raise DomainError("Engine-disabled event does not accept a numeric value.")
            return
        if self.value is None or not isfinite(self.value):
            raise DomainError("Scenario event requires a finite numeric value.")
        if self.event_type in {
            ScenarioEventType.THRUST_SCALE,
            ScenarioEventType.SPECIFIC_IMPULSE_SCALE,
            ScenarioEventType.RESPONSE_TIME_SCALE,
        } and self.value <= 0.0:
            raise DomainError("Scenario scale factors must be greater than zero.")
        if self.event_type in {
            ScenarioEventType.IGNITION_DELAY,
            ScenarioEventType.EARLY_CUTOFF,
        } and self.value < 0.0:
            raise DomainError("Scenario timing values must be non-negative.")


@dataclass(frozen=True, slots=True)
class PropulsionScenario:
    """One mutually exclusive propulsion outcome with declared probability.

    A nominal scenario has an empty event tuple. Probability is supplied by the
    analyst from reliability/test evidence; it is never derived from event
    count or an implicit independence assumption.

    References
    ----------
    NASA-SP-2011-3421: https://ntrs.nasa.gov/citations/20120001369
    NASA-20050207429: https://ntrs.nasa.gov/citations/20050207429
    """

    scenario_id: str
    probability: float
    events: tuple[ScenarioEvent, ...]
    source: str
    rationale: str

    def __post_init__(self) -> None:
        if not self.scenario_id.strip():
            raise DomainError("Propulsion scenario identifier must not be empty.")
        if not isfinite(self.probability) or not 0.0 < self.probability <= 1.0:
            raise DomainError("Propulsion scenario probability must be in (0, 1].")
        if not self.source.strip() or not self.rationale.strip():
            raise DomainError("Propulsion scenario source and rationale are required.")


@dataclass(frozen=True, slots=True)
class PropulsionScenarioSet:
    """Exhaustive set of mutually exclusive propulsion scenarios.

    References
    ----------
    NASA-SP-2011-3421: https://ntrs.nasa.gov/citations/20120001369
    NASA-20160007073: https://ntrs.nasa.gov/citations/20160007073
    """

    scenarios: tuple[PropulsionScenario, ...]
    objective: str

    def __post_init__(self) -> None:
        if not self.scenarios:
            raise DomainError("A propulsion scenario set must not be empty.")
        if not self.objective.strip():
            raise DomainError("Propulsion scenario objective must not be empty.")
        identifiers = tuple(scenario.scenario_id for scenario in self.scenarios)
        if len(identifiers) != len(set(identifiers)):
            raise DomainError("Propulsion scenario identifiers must be unique.")
        probability = sum(scenario.probability for scenario in self.scenarios)
        if abs(probability - 1.0) > 1e-12:
            raise DomainError("Mutually exclusive scenario probabilities must sum to one.")


@dataclass(frozen=True, slots=True)
class ScenarioMetricStatistics:
    """Probability-weighted discrete statistics for one propulsion metric.

    References
    ----------
    NASA-SP-2011-3421: https://ntrs.nasa.gov/citations/20120001369
    """

    metric: str
    expected_value: float
    minimum: float
    p05: float
    p50: float
    p95: float
    maximum: float


@dataclass(frozen=True, slots=True)
class PropulsionScenarioRealization:
    """One probability-labelled cluster artifact and convergence evidence.

    ``evidence`` is absent only for the exact all-engines-disabled zero-thrust
    case, which needs no transient interpolation.

    References
    ----------
    NASA-20050207429: https://ntrs.nasa.gov/citations/20050207429
    NASA-TM-107318: https://ntrs.nasa.gov/citations/19970010379
    """

    scenario_id: str
    probability: float
    artifact: TabulatedBurnArtifact
    evidence: ClusterSamplingEvidence | None


@dataclass(frozen=True, slots=True)
class PropulsionScenarioEnsemble:
    """Discrete propulsion scenario ensemble for later Sidera propagation.

    References
    ----------
    NASA-SP-2011-3421: https://ntrs.nasa.gov/citations/20120001369
    NASA-20050207429: https://ntrs.nasa.gov/citations/20050207429
    NASA-20160007073: https://ntrs.nasa.gov/citations/20160007073
    """

    objective: str
    scenarios: tuple[PropulsionScenario, ...]
    realizations: tuple[PropulsionScenarioRealization, ...]
    statistics: tuple[ScenarioMetricStatistics, ...]
    zero_thrust_probability: float
    reference_ids: tuple[str, ...]
    warnings: tuple[str, ...]


def _schedule_after_cutoff(
    schedule: ThrottleSchedule, cutoff_time_s: float
) -> ThrottleSchedule | None:
    segments: list[ThrottleSegment] = []
    for interval in resolved_schedule_intervals(schedule):
        if interval.start_time_s >= cutoff_time_s:
            segments.append(
                ThrottleSegment(interval.duration_s, 0.0, 0.0, "failure-cutoff-off")
            )
            continue
        if interval.end_time_s <= cutoff_time_s:
            segments.append(
                ThrottleSegment(
                    interval.duration_s,
                    interval.start_throttle,
                    interval.end_throttle,
                    interval.phase,
                )
            )
            continue
        active_duration = cutoff_time_s - interval.start_time_s
        if active_duration > 0.0:
            fraction = active_duration / interval.duration_s
            cutoff_throttle = interval.start_throttle + fraction * (
                interval.end_throttle - interval.start_throttle
            )
            segments.append(
                ThrottleSegment(
                    active_duration,
                    interval.start_throttle,
                    cutoff_throttle,
                    f"{interval.phase}-until-failure-cutoff",
                )
            )
        remaining = interval.end_time_s - cutoff_time_s
        if remaining > 0.0:
            segments.append(ThrottleSegment(remaining, 0.0, 0.0, "failure-cutoff-off"))
    if not segments or sum(segment.exposure_s for segment in segments) <= 0.0:
        return None
    return ThrottleSchedule(tuple(segments), name=f"{schedule.name} with early cutoff")


def _weighted_percentile(
    values: tuple[tuple[float, float], ...], probability: float
) -> float:
    cumulative = 0.0
    for value, weight in sorted(values):
        cumulative += weight
        if cumulative + 1e-15 >= probability:
            return value
    return max(value for value, _ in values)


def _metric_statistics(
    metric: str, values: tuple[tuple[float, float], ...]
) -> ScenarioMetricStatistics:
    return ScenarioMetricStatistics(
        metric=metric,
        expected_value=sum(value * probability for value, probability in values),
        minimum=min(value for value, _ in values),
        p05=_weighted_percentile(values, 0.05),
        p50=_weighted_percentile(values, 0.50),
        p95=_weighted_percentile(values, 0.95),
        maximum=max(value for value, _ in values),
    )


def _zero_artifact(
    duration_s: float,
    *,
    members: tuple[EngineClusterMember, ...],
    commands: tuple[EngineTransientCommand, ...],
    scenario: PropulsionScenario,
) -> TabulatedBurnArtifact:
    source_hash = hashlib.sha256(
        canonical_json_bytes(
            {"members": members, "commands": commands, "scenario": scenario}
        )
    ).hexdigest()
    return TabulatedBurnArtifact(
        segments=(
            TabulatedBurnSegment(
                start_time_s=0.0,
                end_time_s=duration_s,
                start_delivered_thrust_n=0.0,
                end_delivered_thrust_n=0.0,
                start_axial_thrust_n=0.0,
                end_axial_thrust_n=0.0,
                start_total_mass_flow_kg_s=0.0,
                end_total_mass_flow_kg_s=0.0,
                phase="all-engines-disabled",
            ),
        ),
        source_id=f"propulsion-scenario:{scenario.scenario_id}",
        source_sha256=source_hash,
        reference_ids=SCENARIO_REFERENCE_IDS,
        warnings=(
            "All commanded engines are disabled in this explicitly declared scenario.",
        ),
    )


def build_propulsion_scenario_ensemble(
    members: tuple[EngineClusterMember, ...],
    commands: tuple[EngineTransientCommand, ...],
    scenario_set: PropulsionScenarioSet,
    *,
    maximum_step_s: float = 0.25,
    relative_tolerance: float = 1e-7,
    maximum_refinements: int = 12,
) -> PropulsionScenarioEnsemble:
    """Build exact probability-labelled cluster histories for discrete outcomes.

    Events are applied to the same nominal cluster definition. Multiple engine
    IDs in a single event remain one dependent outcome, allowing common-cause
    scenarios without multiplying independent component probabilities.

    References
    ----------
    NASA-SP-2011-3421: https://ntrs.nasa.gov/citations/20120001369
    NASA-20050207429: https://ntrs.nasa.gov/citations/20050207429
    NASA-20160007073: https://ntrs.nasa.gov/citations/20160007073
    NASA-TM-107318: https://ntrs.nasa.gov/citations/19970010379
    """

    if not members or not commands:
        raise DomainError("Scenario analysis requires cluster members and commands.")
    member_map = {member.engine_id: member for member in members}
    command_map = {command.engine_id: command for command in commands}
    if len(member_map) != len(members) or len(command_map) != len(commands):
        raise DomainError("Scenario cluster identifiers must be unique.")
    enabled_ids = {member.engine_id for member in members if member.enabled}
    if enabled_ids != command_map.keys():
        raise DomainError("Every nominally enabled scenario engine needs one command.")
    nominal_duration = max(
        command.ignition_delay_s
        + command.schedule.total_duration_s
        + command.settling_duration_s
        for command in commands
    )
    realizations: list[PropulsionScenarioRealization] = []
    for scenario in scenario_set.scenarios:
        disabled = {member.engine_id for member in members if not member.enabled}
        delays = {engine_id: 0.0 for engine_id in enabled_ids}
        cutoff: dict[str, float] = {}
        thrust_scale = {engine_id: 1.0 for engine_id in enabled_ids}
        isp_scale = {engine_id: 1.0 for engine_id in enabled_ids}
        response_scale = {engine_id: 1.0 for engine_id in enabled_ids}
        applied: set[tuple[ScenarioEventType, str]] = set()
        for event in scenario.events:
            unknown = set(event.engine_ids) - member_map.keys()
            if unknown:
                raise DomainError(
                    f"Scenario {scenario.scenario_id!r} references unknown engines: "
                    f"{', '.join(sorted(unknown))}."
                )
            for engine_id in event.engine_ids:
                if (
                    event.event_type is not ScenarioEventType.ENGINE_DISABLED
                    and engine_id not in enabled_ids
                ):
                    raise DomainError("Scenario modification targets a nominally disabled engine.")
                key = (event.event_type, engine_id)
                if key in applied:
                    raise DomainError("A scenario cannot repeat an event type for one engine.")
                applied.add(key)
                if event.event_type is ScenarioEventType.ENGINE_DISABLED:
                    disabled.add(engine_id)
                elif event.event_type is ScenarioEventType.IGNITION_DELAY:
                    assert event.value is not None
                    delays[engine_id] = event.value
                elif event.event_type is ScenarioEventType.EARLY_CUTOFF:
                    assert event.value is not None
                    if event.value >= command_map[engine_id].schedule.total_duration_s:
                        raise DomainError(
                            "Early-cutoff time must be before the command schedule ends."
                        )
                    cutoff[engine_id] = event.value
                elif event.event_type is ScenarioEventType.THRUST_SCALE:
                    assert event.value is not None
                    thrust_scale[engine_id] = event.value
                elif event.event_type is ScenarioEventType.SPECIFIC_IMPULSE_SCALE:
                    assert event.value is not None
                    isp_scale[engine_id] = event.value
                else:
                    assert event.value is not None
                    response_scale[engine_id] = event.value

        for engine_id in disabled:
            if any(
                applied_type is not ScenarioEventType.ENGINE_DISABLED
                and applied_engine == engine_id
                for applied_type, applied_engine in applied
            ):
                raise DomainError(
                    "A disabled engine cannot carry additional modifications in one scenario."
                )

        scenario_members: list[EngineClusterMember] = []
        scenario_commands: list[EngineTransientCommand] = []
        for engine_id, member in member_map.items():
            if engine_id in disabled:
                scenario_members.append(replace(member, enabled=False))
                continue
            command = command_map[engine_id]
            schedule = command.schedule
            if engine_id in cutoff:
                modified = _schedule_after_cutoff(schedule, cutoff[engine_id])
                if modified is None:
                    scenario_members.append(replace(member, enabled=False))
                    disabled.add(engine_id)
                    continue
                schedule = modified
            point = member.operating_point
            scaled_point = _scaled_operating_point(
                point,
                supply_pressure_pa=point.chamber_pressure_pa,
                delivered_thrust_n=point.delivered_thrust_n * thrust_scale[engine_id],
                system_specific_impulse_s=(
                    point.system_specific_impulse_s * isp_scale[engine_id]
                ),
                chamber_pressure_pa=point.chamber_pressure_pa * thrust_scale[engine_id],
                model_id="discrete_propulsion_scenario_v1",
                source=f"{point.source}; scenario {scenario.scenario_id}",
                warning="Operating point contains an explicitly declared scenario modification.",
            )
            scenario_members.append(replace(member, operating_point=scaled_point))
            scenario_commands.append(
                replace(
                    command,
                    schedule=schedule,
                    ignition_delay_s=command.ignition_delay_s + delays[engine_id],
                    response_time_constant_s=(
                        command.response_time_constant_s * response_scale[engine_id]
                    ),
                )
            )

        active_members = tuple(member for member in scenario_members if member.enabled)
        if not active_members:
            artifact = _zero_artifact(
                nominal_duration,
                members=members,
                commands=commands,
                scenario=scenario,
            )
            evidence = None
        else:
            profile = build_engine_cluster_profile(
                tuple(scenario_members),
                tuple(scenario_commands),
                maximum_step_s=maximum_step_s,
                relative_tolerance=relative_tolerance,
                maximum_refinements=maximum_refinements,
            )
            scenario_source_hash = hashlib.sha256(
                canonical_json_bytes(
                    {
                        "cluster_source_sha256": profile.artifact.source_sha256,
                        "scenario": scenario,
                    }
                )
            ).hexdigest()
            artifact = replace(
                profile.artifact,
                source_id=f"propulsion-scenario:{scenario.scenario_id}",
                source_sha256=scenario_source_hash,
                reference_ids=tuple(
                    dict.fromkeys(profile.artifact.reference_ids + SCENARIO_REFERENCE_IDS)
                ),
                warnings=profile.artifact.warnings
                + ("Artifact belongs to an explicitly probability-labelled scenario.",),
            )
            evidence = profile.evidence
        realizations.append(
            PropulsionScenarioRealization(
                scenario_id=scenario.scenario_id,
                probability=scenario.probability,
                artifact=artifact,
                evidence=evidence,
            )
        )

    metric_values = {
        "duration_s": tuple(
            (item.artifact.duration_s, item.probability) for item in realizations
        ),
        "delivered_impulse_n_s": tuple(
            (item.artifact.delivered_total_impulse_n_s, item.probability)
            for item in realizations
        ),
        "axial_impulse_n_s": tuple(
            (item.artifact.axial_total_impulse_n_s, item.probability)
            for item in realizations
        ),
        "consumed_propellant_kg": tuple(
            (item.artifact.consumed_propellant_kg, item.probability)
            for item in realizations
        ),
        "thrust_centroid_time_s": tuple(
            (item.artifact.thrust_centroid_time_s, item.probability)
            for item in realizations
        ),
    }
    return PropulsionScenarioEnsemble(
        objective=scenario_set.objective,
        scenarios=scenario_set.scenarios,
        realizations=tuple(realizations),
        statistics=tuple(
            _metric_statistics(metric, values) for metric, values in metric_values.items()
        ),
        zero_thrust_probability=sum(
            item.probability
            for item in realizations
            if item.artifact.delivered_total_impulse_n_s == 0.0
        ),
        reference_ids=SCENARIO_REFERENCE_IDS,
        warnings=(
            "Scenario probabilities are analyst-supplied and are not inferred by this library.",
            "Scenario metrics are propulsion consequences; trajectory consequences require Sidera.",
        ),
    )
