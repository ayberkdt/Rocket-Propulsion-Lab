"""Engine-cluster composition and convergence-qualified transient sampling.

Each enabled cluster member contributes its rated thrust and stream flows
scaled by a realized first-order throttle history.  The summed history is
exported as a piecewise-linear exchange artifact whose time step is halved
until its impulse, axial impulse, propellant mass and thrust centroid agree
with the exact response integrals.  [1]_ documents component transient-model
assembly and comparison with engine data; [2]_ gives the transient-modelling
context.

References
----------
.. [1] NASA TM-107318. NASA, "RL10A-3-3A Rocket Engine Modeling Project".
       https://ntrs.nasa.gov/citations/19970010379
.. [2] NASA, "Transient Mathematical Modeling for Liquid Rocket Engine
       Systems". https://ntrs.nasa.gov/citations/20040000363
"""

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

    Parameters
    ----------
    engine_id : str
        Non-blank identifier, unique within the cluster.
    operating_point : EngineOperatingPoint
        Rated (full-command) operating point of the member.
    cant_efficiency : float, default 1.0
        Axial fraction of delivered thrust, in ``(0, 1]``.
    enabled : bool, default True
        Whether the member participates in the cluster.

    Raises
    ------
    DomainError
        If the identifier is blank or the cant efficiency is outside
        ``(0, 1]``.

    Notes
    -----
    Liquid-engine system design context is given in [1]_.

    References
    ----------
    .. [1] NASA SP-125. NASA, "Design of Liquid Propellant Rocket Engines".
           https://ntrs.nasa.gov/citations/19710019929
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
    """Throttle schedule and first-order response policy for one member.

    Parameters
    ----------
    engine_id : str
        Non-blank identifier of the commanded cluster member.
    schedule : ThrottleSchedule
        Commanded throttle schedule.
    response_time_constant_s : float, default 0.0
        First-order response time constant, s.
    ignition_delay_s : float, default 0.0
        Delay between command and response onset, s.
    settling_duration_s : float, default 0.0
        Settling duration passed to the response realization, s.

    Raises
    ------
    DomainError
        If the identifier is blank or a duration is non-finite or negative.

    Notes
    -----
    [1]_ documents start, steady-state and shutdown transients of a modelled
    engine against test data.

    References
    ----------
    .. [1] NASA TM-107318. NASA, "RL10A-3-3A Rocket Engine Modeling Project".
           https://ntrs.nasa.gov/citations/19970010379
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
    """Convergence evidence of the tabulation against exact response integrals.

    The errors document the numerical closure of the tabulation; [1]_ gives
    the transient-modelling context.

    Attributes
    ----------
    maximum_step_s : float
        Sampling step of the accepted tabulation, s.
    refinement_count : int
        Number of step halvings performed.
    relative_tolerance : float
        Requested relative tolerance.
    exact_delivered_impulse_n_s : float
        Exact delivered impulse of the realized histories, N s.
    exact_axial_impulse_n_s : float
        Exact axial impulse, N s.
    exact_propellant_mass_kg : float
        Exact tank-propellant consumption, kg.
    exact_thrust_centroid_time_s : float
        Exact delivered-thrust centroid time, s.
    impulse_error_n_s : float
        Absolute delivered-impulse error of the tabulation, N s.
    axial_impulse_error_n_s : float
        Absolute axial-impulse error, N s.
    propellant_error_kg : float
        Absolute propellant-mass error, kg.
    centroid_error_s : float
        Absolute thrust-centroid error, s.

    References
    ----------
    .. [1] NASA, "Transient Mathematical Modeling for Liquid Rocket Engine
           Systems". https://ntrs.nasa.gov/citations/20040000363
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
    """Target-neutral transient cluster result with convergence evidence.

    Engine-level transient histories are combined into one cluster history,
    in the spirit of the component transient-model assembly of [1]_.

    Attributes
    ----------
    artifact : TabulatedBurnArtifact
        Piecewise-linear summed force and mass-flow history.
    evidence : ClusterSamplingEvidence
        Convergence evidence against the exact response integrals.
    histories : tuple of (str, RealizedThrottleHistory)
        Realized throttle history per engine, sorted by identifier.

    References
    ----------
    .. [1] NASA TM-107318. NASA, "RL10A-3-3A Rocket Engine Modeling Project".
           https://ntrs.nasa.gov/citations/19970010379
    """

    artifact: TabulatedBurnArtifact
    evidence: ClusterSamplingEvidence
    histories: tuple[tuple[str, RealizedThrottleHistory], ...]


def _timeline(
    histories: dict[str, RealizedThrottleHistory], maximum_step_s: float
) -> tuple[float, ...]:
    """Merge every history's sample times, the origin and the latest end time."""
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
    """Sum the cluster state at one instant on the requested one-sided limit.

    Each member's rated delivered thrust, axial thrust and tank-depleting
    stream flows are multiplied by its realized throttle; members whose
    throttle exceeds ``1e-14`` are reported active.  Thrust and mass flow scale
    with the same factor, so each member's specific impulse is constant along
    its response.

    Returns
    -------
    tuple
        ``(delivered_n, axial_n, total_flow_kg_s, tank_flows, active_ids)``.
    """
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
    """Tabulate the summed cluster history on one sampling step.

    Each segment interpolates linearly between the right-hand limit at its
    start and the left-hand limit at its end, so throttle discontinuities at
    sample boundaries are represented without smearing.
    """
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
    """Build a converged cluster tabulation while retaining exact integrals.

    Parameters
    ----------
    members : tuple of EngineClusterMember
        Cluster members with unique identifiers; at least one enabled.
    commands : tuple of EngineTransientCommand
        Exactly one command per enabled member, with unique identifiers.
    maximum_step_s : float, default 0.25
        Initial sampling step, s.
    relative_tolerance : float, default 1e-7
        Relative tolerance on impulse, axial impulse, propellant mass and
        thrust centroid (absolute floor ``1e-11``).
    maximum_refinements : int, default 12
        Maximum number of step halvings.

    Returns
    -------
    EngineClusterProfile
        Tabulated artifact, convergence evidence and realized histories.

    Raises
    ------
    DomainError
        If members or commands are empty or duplicated, no member is enabled,
        the enabled members and commands do not match one to one, a numerical
        argument is outside its domain, or sampling does not converge within
        ``maximum_refinements`` halvings.

    Notes
    -----
    Each command is realized as a first-order throttle response with ignition
    delay.  With ``E_i`` the time integral and ``M_i`` the first time moment of
    member ``i``'s realized throttle, ``F_i`` its rated delivered thrust,
    ``eta_i`` its cant efficiency and ``mdot_i`` its rated tank flow, the exact
    delivered impulse is ``sum(F_i * E_i)``, the axial impulse
    ``sum(F_i * eta_i * E_i)``, the propellant mass ``sum(mdot_i * E_i)`` and
    the thrust centroid ``sum(F_i * M_i) / sum(F_i * E_i)``.  The sampling
    step is halved until the tabulation reproduces all four quantities within
    tolerance.  [1]_ and [2]_ give the engine transient-modelling context.

    References
    ----------
    .. [1] NASA TM-107318. NASA, "RL10A-3-3A Rocket Engine Modeling Project".
           https://ntrs.nasa.gov/citations/19970010379
    .. [2] NASA, "Transient Mathematical Modeling for Liquid Rocket Engine
           Systems". https://ntrs.nasa.gov/citations/20040000363
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
