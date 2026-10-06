"""Engine-cluster composition and converged transient exchange sampling."""

from __future__ import annotations

import hashlib
from dataclasses import dataclass
from itertools import pairwise
from math import isfinite

from rocket_propulsion.core.errors import DomainError

from .models import EngineOperatingPoint, ThrottleSchedule
from .response import RealizedThrottleHistory, realize_first_order_schedule
from .serialization import canonical_json_bytes
from .tabulated import TabulatedBurnArtifact, TabulatedBurnSegment


@dataclass(frozen=True, slots=True)
class EngineClusterMember:
    """One independently commanded engine or resolved engine group.

    References
    ----------
    NASA-SP-125: https://ntrs.nasa.gov/citations/19710019929
    """

    engine_id: str
    operating_point: EngineOperatingPoint
    cant_efficiency: float = 1.0
    enabled: bool = True

    def __post_init__(self) -> None:
        if not self.engine_id.strip():
            raise DomainError("Cluster engine identifier must not be empty.")
        if not isfinite(self.cant_efficiency) or not 0.0 < self.cant_efficiency <= 1.0:
            raise DomainError("Cluster member cant efficiency must be in (0, 1].")


@dataclass(frozen=True, slots=True)
class EngineTransientCommand:
    """Command and response policy for one cluster member.

    References
    ----------
    NASA-TM-107318: https://ntrs.nasa.gov/citations/19970010379
    """

    engine_id: str
    schedule: ThrottleSchedule
    response_time_constant_s: float = 0.0
    ignition_delay_s: float = 0.0
    settling_duration_s: float = 0.0

    def __post_init__(self) -> None:
        if not self.engine_id.strip():
            raise DomainError("Transient command engine identifier must not be empty.")
        for value, label in (
            (self.response_time_constant_s, "Response time constant"),
            (self.ignition_delay_s, "Ignition delay"),
            (self.settling_duration_s, "Settling duration"),
        ):
            if not isfinite(value) or value < 0.0:
                raise DomainError(f"{label} must be finite and non-negative.")


@dataclass(frozen=True, slots=True)
class ClusterSamplingEvidence:
    """Exact-response versus tabulated-interpolation convergence evidence.

    References
    ----------
    NASA transient-modeling overview:
    https://ntrs.nasa.gov/citations/20040000363
    """

    maximum_step_s: float
    refinement_count: int
    relative_tolerance: float
    exact_delivered_impulse_n_s: float
    exact_axial_impulse_n_s: float
    exact_propellant_mass_kg: float
    exact_thrust_centroid_time_s: float
    impulse_error_n_s: float
    axial_impulse_error_n_s: float
    propellant_error_kg: float
    centroid_error_s: float


@dataclass(frozen=True, slots=True)
class EngineClusterProfile:
    """Target-neutral transient cluster result and convergence evidence.

    References
    ----------
    NASA-TM-107318: https://ntrs.nasa.gov/citations/19970010379
    """

    artifact: TabulatedBurnArtifact
    evidence: ClusterSamplingEvidence
    histories: tuple[tuple[str, RealizedThrottleHistory], ...]


def _timeline(
    histories: dict[str, RealizedThrottleHistory], maximum_step_s: float
) -> tuple[float, ...]:
    times = {0.0, max(history.duration_s for history in histories.values())}
    for history in histories.values():
        times.update(history.sample_times(maximum_step_s))
    return tuple(sorted(times))


def _state(
    time_s: float,
    *,
    side: str,
    members: dict[str, EngineClusterMember],
    histories: dict[str, RealizedThrottleHistory],
) -> tuple[float, float, float, tuple[tuple[str, float], ...], tuple[str, ...]]:
    delivered = 0.0
    axial = 0.0
    total_flow = 0.0
    tank_flows: dict[str, float] = {}
    active: list[str] = []
    for engine_id, member in members.items():
        throttle = histories[engine_id].realized_throttle(time_s, side=side)
        if throttle > 1e-14:
            active.append(engine_id)
        point = member.operating_point
        delivered += point.delivered_thrust_n * throttle
        axial += point.delivered_thrust_n * member.cant_efficiency * throttle
        total_flow += point.total_tank_flow_kg_s * throttle
        for stream in point.streams:
            if stream.tank_depleting:
                tank_flows[stream.tank_id] = (
                    tank_flows.get(stream.tank_id, 0.0) + stream.mass_flow_kg_s * throttle
                )
    return delivered, axial, total_flow, tuple(sorted(tank_flows.items())), tuple(active)


def _artifact_for_step(
    *,
    members: dict[str, EngineClusterMember],
    histories: dict[str, RealizedThrottleHistory],
    maximum_step_s: float,
    source_hash: str,
) -> TabulatedBurnArtifact:
    times = _timeline(histories, maximum_step_s)
    segments: list[TabulatedBurnSegment] = []
    for start_time, end_time in pairwise(times):
        start = _state(
            start_time, side="right", members=members, histories=histories
        )
        end = _state(end_time, side="left", members=members, histories=histories)
        active = sorted(set(start[4]) | set(end[4]))
        phase = "cluster-off" if not active else f"cluster-active:{'+'.join(active)}"
        segments.append(
            TabulatedBurnSegment(
                start_time_s=start_time,
                end_time_s=end_time,
                start_delivered_thrust_n=start[0],
                end_delivered_thrust_n=end[0],
                start_axial_thrust_n=start[1],
                end_axial_thrust_n=end[1],
                start_total_mass_flow_kg_s=start[2],
                end_total_mass_flow_kg_s=end[2],
                start_stream_mass_flows_kg_s=start[3],
                end_stream_mass_flows_kg_s=end[3],
                phase=phase,
            )
        )
    return TabulatedBurnArtifact(
        segments=tuple(segments),
        source_id="engine-cluster-response",
        source_sha256=source_hash,
        reference_ids=(
            "NASA-TM-107318",
            "NASA-20040000363",
            "NASA-GRC-THRUST-EQUATION",
        ),
        warnings=(
            (
                "First-order engine response is exported as a convergence-qualified "
                "piecewise-linear propulsion history."
            ),
        ),
    )


def build_engine_cluster_profile(
    members: tuple[EngineClusterMember, ...],
    commands: tuple[EngineTransientCommand, ...],
    *,
    maximum_step_s: float = 0.25,
    relative_tolerance: float = 1e-7,
    maximum_refinements: int = 12,
) -> EngineClusterProfile:
    """Build a converged tabulation while retaining exact response integrals.

    References
    ----------
    NASA-TM-107318: https://ntrs.nasa.gov/citations/19970010379
    NASA transient-modeling overview:
    https://ntrs.nasa.gov/citations/20040000363
    """

    if not members or not commands:
        raise DomainError("Cluster members and transient commands must not be empty.")
    member_map = {member.engine_id: member for member in members}
    enabled_members = {
        engine_id: member for engine_id, member in member_map.items() if member.enabled
    }
    command_map = {command.engine_id: command for command in commands}
    if len(member_map) != len(members) or len(command_map) != len(commands):
        raise DomainError("Cluster member and command engine identifiers must be unique.")
    if not enabled_members:
        raise DomainError("An engine cluster requires at least one enabled member.")
    if enabled_members.keys() != command_map.keys():
        raise DomainError("Every enabled cluster member requires exactly one transient command.")
    if not isfinite(maximum_step_s) or maximum_step_s <= 0.0:
        raise DomainError("Maximum cluster sample step must be positive and finite.")
    if not isfinite(relative_tolerance) or relative_tolerance <= 0.0:
        raise DomainError("Cluster relative tolerance must be positive and finite.")
    if maximum_refinements < 0:
        raise DomainError("Maximum cluster refinements must be non-negative.")

    histories = {
        engine_id: realize_first_order_schedule(
            command.schedule,
            time_constant_s=command.response_time_constant_s,
            ignition_delay_s=command.ignition_delay_s,
            settling_duration_s=command.settling_duration_s,
        )
        for engine_id, command in command_map.items()
    }
    exact_delivered = sum(
        enabled_members[engine_id].operating_point.delivered_thrust_n * history.exposure_s
        for engine_id, history in histories.items()
    )
    exact_axial = sum(
        enabled_members[engine_id].operating_point.delivered_thrust_n
        * enabled_members[engine_id].cant_efficiency
        * history.exposure_s
        for engine_id, history in histories.items()
    )
    exact_mass = sum(
        enabled_members[engine_id].operating_point.total_tank_flow_kg_s
        * history.exposure_s
        for engine_id, history in histories.items()
    )
    exact_moment = sum(
        enabled_members[engine_id].operating_point.delivered_thrust_n
        * history.throttle_first_moment_s2
        for engine_id, history in histories.items()
    )
    exact_centroid = exact_moment / exact_delivered if exact_delivered > 0.0 else 0.0
    source_hash = hashlib.sha256(
        canonical_json_bytes({"members": members, "commands": commands})
    ).hexdigest()

    step = maximum_step_s
    artifact: TabulatedBurnArtifact | None = None
    errors: tuple[float, float, float, float] | None = None
    refinement = 0
    for refinement in range(maximum_refinements + 1):
        artifact = _artifact_for_step(
            members=enabled_members,
            histories=histories,
            maximum_step_s=step,
            source_hash=source_hash,
        )
        errors = (
            abs(artifact.delivered_total_impulse_n_s - exact_delivered),
            abs(artifact.axial_total_impulse_n_s - exact_axial),
            abs(artifact.consumed_propellant_kg - exact_mass),
            abs(artifact.thrust_centroid_time_s - exact_centroid),
        )
        scales = (exact_delivered, exact_axial, exact_mass, exact_centroid)
        if all(
            error <= max(1e-11, abs(scale) * relative_tolerance)
            for error, scale in zip(errors, scales, strict=True)
        ):
            break
        step *= 0.5
    else:
        raise DomainError("Cluster transient sampling did not meet the requested tolerance.")
    assert artifact is not None and errors is not None
    evidence = ClusterSamplingEvidence(
        maximum_step_s=step,
        refinement_count=refinement,
        relative_tolerance=relative_tolerance,
        exact_delivered_impulse_n_s=exact_delivered,
        exact_axial_impulse_n_s=exact_axial,
        exact_propellant_mass_kg=exact_mass,
        exact_thrust_centroid_time_s=exact_centroid,
        impulse_error_n_s=errors[0],
        axial_impulse_error_n_s=errors[1],
        propellant_error_kg=errors[2],
        centroid_error_s=errors[3],
    )
    return EngineClusterProfile(
        artifact=artifact,
        evidence=evidence,
        histories=tuple(sorted(histories.items())),
    )
