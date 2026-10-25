"""Discrete propulsion execution and failure scenarios for engine clusters.

Scenarios are mutually exclusive, explicitly sourced and probability weighted,
following the discrete-outcome treatment of probabilistic risk assessment
[1]_.  The module never invents component failure rates and never assumes
independence: engine-out, start and cutoff outcomes [2]_ and dependent
(common-cause) failures [3]_ are declared by the analyst.  Each scenario
produces a complete, target-neutral thrust/mass-flow artifact that a
trajectory propagator such as Sidera consumes; trajectory consequences are
outside this module.

References
----------
.. [1] NASA/SP-2011-3421. NASA, "Probabilistic Risk Assessment Procedures
       Guide for NASA Managers and Practitioners",
       https://ntrs.nasa.gov/citations/20120001369.
.. [2] NASA-20050207429. NASA, "Key Reliability Drivers of Liquid
       Propulsion Engines and a Reliability Model for Sensitivity Analysis",
       https://ntrs.nasa.gov/citations/20050207429.
.. [3] NASA-20160007073. NASA, "Common Cause Failure Modeling in Space
       Launch Vehicles", https://ntrs.nasa.gov/citations/20160007073.
"""

from __future__ import annotations

import hashlib
from dataclasses import dataclass, replace
from enum import Enum
from fractions import Fraction
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
from .numeric_outputs import _finite_rational, _fraction
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

    Members cover engine-out, ignition delay and early cutoff outcomes, which
    are reliability drivers of liquid engines [1]_, and multiplicative
    thrust, specific-impulse and response-time deviations.  Several engines
    in one event express a dependent outcome [2]_.

    References
    ----------
    .. [1] NASA-20050207429. NASA, "Key Reliability Drivers of Liquid
           Propulsion Engines and a Reliability Model for Sensitivity Analysis",
           https://ntrs.nasa.gov/citations/20050207429.
    .. [2] NASA-20160007073. NASA, "Common Cause Failure Modeling in Space
           Launch Vehicles", https://ntrs.nasa.gov/citations/20160007073.
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

    Parameters
    ----------
    event_type : ScenarioEventType
        Kind of modification.
    engine_ids : tuple of str
        Non-empty, unique engine identifiers affected by the event.
    value : float or None
        ``None`` for ``ENGINE_DISABLED``.  Otherwise finite: a positive factor
        for the scale events, or a non-negative time in s for
        ``IGNITION_DELAY`` (added delay) and ``EARLY_CUTOFF`` (cutoff time on
        the command schedule).
    source : str
        Non-empty provenance of the event.
    rationale : str
        Non-empty engineering rationale.

    Raises
    ------
    DomainError
        If identifiers are blank or duplicated, provenance is missing, or the
        value violates the rule for its event type.

    Notes
    -----
    Multiple engine identifiers in one event express an explicitly declared
    dependent (common-cause) outcome [2]_ within a discrete scenario
    model [1]_.  Dependence is never inferred from identical hardware.

    References
    ----------
    .. [1] NASA/SP-2011-3421. NASA, "Probabilistic Risk Assessment Procedures
           Guide for NASA Managers and Practitioners",
           https://ntrs.nasa.gov/citations/20120001369.
    .. [2] NASA-20160007073. NASA, "Common Cause Failure Modeling in Space
           Launch Vehicles", https://ntrs.nasa.gov/citations/20160007073.
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

    Parameters
    ----------
    scenario_id : str
        Non-empty identifier, unique within its scenario set.
    probability : float
        Outcome probability in ``(0, 1]``.
    events : tuple of ScenarioEvent
        Modifications; an empty tuple is the nominal scenario.
    source : str
        Non-empty provenance of the probability.
    rationale : str
        Non-empty engineering rationale.

    Raises
    ------
    DomainError
        If the identifier, source or rationale is blank, or the probability is
        outside ``(0, 1]``.

    Notes
    -----
    The probability is supplied by the analyst from reliability or test
    evidence [1]_, [2]_; it is never derived from the event count or from an
    implicit independence assumption.

    References
    ----------
    .. [1] NASA/SP-2011-3421. NASA, "Probabilistic Risk Assessment Procedures
           Guide for NASA Managers and Practitioners",
           https://ntrs.nasa.gov/citations/20120001369.
    .. [2] NASA-20050207429. NASA, "Key Reliability Drivers of Liquid
           Propulsion Engines and a Reliability Model for Sensitivity Analysis",
           https://ntrs.nasa.gov/citations/20050207429.
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

    Parameters
    ----------
    scenarios : tuple of PropulsionScenario
        Non-empty scenarios with unique identifiers whose probabilities sum to
        one within an absolute tolerance of ``1e-12``.
    objective : str
        Non-empty statement of the analysis objective.

    Raises
    ------
    DomainError
        If the set is empty, the objective is blank, identifiers repeat, or
        the probabilities do not sum to one.

    Notes
    -----
    Exhaustive, mutually exclusive outcomes form a discrete probability
    distribution [1]_; dependent failures are represented as their own
    outcomes rather than as products of independent probabilities [2]_.

    References
    ----------
    .. [1] NASA/SP-2011-3421. NASA, "Probabilistic Risk Assessment Procedures
           Guide for NASA Managers and Practitioners",
           https://ntrs.nasa.gov/citations/20120001369.
    .. [2] NASA-20160007073. NASA, "Common Cause Failure Modeling in Space
           Launch Vehicles", https://ntrs.nasa.gov/citations/20160007073.
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

    Parameters
    ----------
    metric : str
        Metric name, for example ``"delivered_impulse_n_s"``.
    expected_value : float
        Normalized probability-weighted mean [2]_.
    minimum, maximum : float
        Smallest and largest metric value with positive probability.
    p05, p50, p95 : float
        Discrete quantiles: the smallest value whose cumulative probability
        reaches 0.05, 0.5 and 0.95.

    Notes
    -----
    The statistics describe the analyst-declared discrete distribution [1]_;
    they are not calibrated probability evidence.

    References
    ----------
    .. [1] NASA/SP-2011-3421. NASA, "Probabilistic Risk Assessment Procedures
           Guide for NASA Managers and Practitioners",
           https://ntrs.nasa.gov/citations/20120001369.
    .. [2] NIST-WEIGHTED-MEAN. NIST Dataplot, "Weighted Mean", eq. 2-20,
           https://www.itl.nist.gov/div898/software/dataplot/refman2/ch2/weigmean.pdf.
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
    """One probability-labelled cluster artifact with its convergence evidence.

    Parameters
    ----------
    scenario_id : str
        Identifier of the realized scenario.
    probability : float
        Declared scenario probability.
    artifact : TabulatedBurnArtifact
        Cluster thrust and mass-flow history under the scenario's outcomes
        [1]_, assembled from per-engine transient models [2]_.
    evidence : ClusterSamplingEvidence or None
        Sampling-convergence evidence; ``None`` only for the exact
        all-engines-disabled zero-thrust case, which needs no transient
        interpolation.

    References
    ----------
    .. [1] NASA-20050207429. NASA, "Key Reliability Drivers of Liquid
           Propulsion Engines and a Reliability Model for Sensitivity Analysis",
           https://ntrs.nasa.gov/citations/20050207429.
    .. [2] NASA-TM-107318. NASA, "RL10A-3-3A Rocket Engine Modeling Project",
           https://ntrs.nasa.gov/citations/19970010379.
    """

    scenario_id: str
    probability: float
    artifact: TabulatedBurnArtifact
    evidence: ClusterSamplingEvidence | None


@dataclass(frozen=True, slots=True)
class PropulsionScenarioEnsemble:
    """Discrete propulsion scenario ensemble ready for trajectory propagation.

    Parameters
    ----------
    objective : str
        Analysis objective copied from the scenario set.
    scenarios : tuple of PropulsionScenario
        The declared scenarios.
    realizations : tuple of PropulsionScenarioRealization
        One artifact per scenario, in scenario order.
    statistics : tuple of ScenarioMetricStatistics
        Probability-weighted statistics of duration, delivered and axial
        impulse, consumed propellant and thrust-centroid time.
    zero_thrust_probability : float
        Normalized probability of the scenarios with zero delivered impulse.
    reference_ids : tuple of str
        Registered reference identifiers of the scenario method.
    warnings : tuple of str
        Interpretation warnings.

    Notes
    -----
    The ensemble holds propulsion consequences only; a trajectory propagator
    such as Sidera evaluates the trajectory consequences of each artifact.
    Probabilities follow the discrete-outcome treatment of [1]_, with
    engine-out and transient outcomes [2]_ and dependent failures [3]_
    declared by the analyst.

    References
    ----------
    .. [1] NASA/SP-2011-3421. NASA, "Probabilistic Risk Assessment Procedures
           Guide for NASA Managers and Practitioners",
           https://ntrs.nasa.gov/citations/20120001369.
    .. [2] NASA-20050207429. NASA, "Key Reliability Drivers of Liquid
           Propulsion Engines and a Reliability Model for Sensitivity Analysis",
           https://ntrs.nasa.gov/citations/20050207429.
    .. [3] NASA-20160007073. NASA, "Common Cause Failure Modeling in Space
           Launch Vehicles", https://ntrs.nasa.gov/citations/20160007073.
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
    """Return ``schedule`` with zero throttle from ``cutoff_time_s`` onward.

    The interval containing the cutoff is split, with its throttle linearly
    interpolated at the cutoff time.  ``None`` means no positive exposure
    remains, so the engine is effectively disabled.
    """
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
    """Return the smallest value whose cumulative weight reaches ``probability``."""
    from .probability_numerics import _weight_ledger
    return _weight_ledger(values).quantile(probability)


def _metric_statistics(
    metric: str, values: tuple[tuple[float, float], ...]
) -> ScenarioMetricStatistics:
    """Return the weighted mean [1]_, extremes and discrete quantiles of ``values``.

    ``values`` holds ``(metric_value, probability)`` pairs; weights are
    normalized with exact rational arithmetic before export.

    References
    ----------
    .. [1] NIST-WEIGHTED-MEAN. NIST Dataplot, "Weighted Mean", eq. 2-20,
           https://www.itl.nist.gov/div898/software/dataplot/refman2/ch2/weigmean.pdf.
    """
    from .numeric_outputs import _finite_rational
    from .probability_numerics import _quantiles, _weight_ledger
    ledger = _weight_ledger(values)
    p05, p50, p95 = _quantiles(ledger)
    return ScenarioMetricStatistics(
        metric=metric,
        expected_value=_finite_rational(ledger.mean, "Discrete weighted mean"),
        minimum=_finite_rational(ledger.rows[0][0], "Discrete weighted minimum"),
        p05=p05, p50=p50, p95=p95,
        maximum=_finite_rational(ledger.rows[-1][0], "Discrete weighted maximum"),
    )


def _zero_artifact(
    duration_s: float,
    *,
    members: tuple[EngineClusterMember, ...],
    commands: tuple[EngineTransientCommand, ...],
    scenario: PropulsionScenario,
) -> TabulatedBurnArtifact:
    """Return one zero-thrust, zero-flow segment for an all-disabled scenario.

    The source hash binds the nominal members, commands and scenario.
    """
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
    """Build probability-labelled cluster histories for discrete outcomes.

    Parameters
    ----------
    members : tuple of EngineClusterMember
        Nominal cluster members with unique engine identifiers.
    commands : tuple of EngineTransientCommand
        Exactly one command per nominally enabled member.
    scenario_set : PropulsionScenarioSet
        Exhaustive, mutually exclusive scenarios.
    maximum_step_s : float, optional
        Maximum sampling step of each cluster profile, in s.
    relative_tolerance : float, optional
        Relative sampling-convergence tolerance of each cluster profile.
    maximum_refinements : int, optional
        Maximum number of sampling refinements per cluster profile.

    Returns
    -------
    PropulsionScenarioEnsemble
        One artifact per scenario plus probability-weighted metric statistics.

    Raises
    ------
    DomainError
        If members or commands are missing or duplicated, a command is
        missing for an enabled engine, an event references an unknown or
        nominally disabled engine, an event type repeats for one engine, a
        disabled engine carries other modifications, or an early cutoff is not
        before the end of its schedule.

    Notes
    -----
    Every scenario modifies the same nominal cluster definition.  The thrust
    scale multiplies delivered thrust and chamber pressure together; the
    specific-impulse scale multiplies the system specific impulse; the
    response scale multiplies the response time constant; an early cutoff
    replaces the remaining schedule with zero throttle.  An engine whose
    cutoff leaves no positive exposure is treated as disabled.  Multiple
    engine identifiers in one event remain one dependent outcome [3]_, so
    common-cause scenarios never multiply independent component
    probabilities [1]_.  Engine-out and start/cutoff outcomes follow the
    reliability drivers of [2]_; each realized cluster is assembled from
    per-engine transient models in the manner of [4]_.  The model is not
    validated against hardware test data.

    References
    ----------
    .. [1] NASA/SP-2011-3421. NASA, "Probabilistic Risk Assessment Procedures
           Guide for NASA Managers and Practitioners",
           https://ntrs.nasa.gov/citations/20120001369.
    .. [2] NASA-20050207429. NASA, "Key Reliability Drivers of Liquid
           Propulsion Engines and a Reliability Model for Sensitivity Analysis",
           https://ntrs.nasa.gov/citations/20050207429.
    .. [3] NASA-20160007073. NASA, "Common Cause Failure Modeling in Space
           Launch Vehicles", https://ntrs.nasa.gov/citations/20160007073.
    .. [4] NASA-TM-107318. NASA, "RL10A-3-3A Rocket Engine Modeling Project",
           https://ntrs.nasa.gov/citations/19970010379.
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
        zero_thrust_probability=_finite_rational(sum((
            _fraction(item.probability)
            for item in realizations
            if item.artifact.delivered_total_impulse_n_s == 0.0
        ), Fraction(0))/sum((_fraction(item.probability) for item in realizations), Fraction(0)), "Discrete zero-thrust probability"),
        reference_ids=SCENARIO_REFERENCE_IDS,
        warnings=(
            "Scenario probabilities are analyst-supplied and are not inferred by this library.",
            "Scenario metrics are propulsion consequences; trajectory consequences require Sidera.",
        ),
    )
