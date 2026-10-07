"""Optional, fail-closed adapter to the Lunaris lunar-orbit propagation framework.

Lunaris (``lunaris.core.propagation``) propagates a 6- or 7-element state
(mass at ``y[6]``) and supports *impulsive* maneuvers
(``ImpulsiveManeuver(t_burn_s, dv_mps, frame, isp_s)`` in a strictly
time-ordered ``ManeuverPlan``) whose mass drop is
``m exp(-|dv| / (g0 Isp))``.  It accepts any dynamics object with
``build_rhs()``.  This adapter offers two exact-where-possible paths without
modifying Lunaris:

1. **Impulsive equivalents.**  Each burn of a :class:`PropulsionTimeline`
   becomes one impulse with the burn's free-space delta-v, placed at the
   acceleration-weighted centroid ``t*``, and with
   ``isp_s = c_eff / g0`` where ``c_eff = dv / ln(m0 / mf)``.  Lunaris'
   ideal-rocket mass law then reproduces the finite burn's final mass
   exactly.  Only the time distribution of the impulse is approximated; the
   burn duration is reported with each impulse.
2. **Finite burns.**  :class:`LunarisFiniteBurnDynamics` wraps any Lunaris
   dynamics and adds the timeline's thrust and mass rate to its right-hand
   side.  :func:`lunaris_segmentation_plan` places zero-delta-v impulses at
   every burn discontinuity, so Lunaris restarts integration exactly there
   (Lunaris segments propagation at maneuver times).

Lunaris carries vehicle mass but no tank states, so finite-burn bridges must
be mass-only (no stream bindings).  Ambient back-pressure, performance
surfaces and jettison events have no Lunaris counterpart and are refused.
"""

from __future__ import annotations

from collections.abc import Callable, Iterator, Mapping
from dataclasses import dataclass
from importlib import import_module
from itertools import pairwise
from math import sqrt
from typing import Any

from rocket_propulsion.core.errors import DomainError, FeatureUnavailableError
from rocket_propulsion.core.validation import real, text, vector3
from rocket_propulsion.propulsion.burns import (
    PropulsionTimeline,
    free_space_burn_metrics,
)
from rocket_propulsion.propulsion.performance import STANDARD_GRAVITY_M_S2

LUNARIS_FRAMES = ("inertial", "ric")
_BOUNDARY_TOLERANCE_S = 1e-9


def _lunaris_plans() -> Any:
    try:
        return import_module("lunaris.core.propagation.plans")
    except ImportError as error:
        raise FeatureUnavailableError(
            "Lunaris is not installed; impulsive equivalents are still available as data."
        ) from error


@dataclass(frozen=True, slots=True)
class LunarisImpulsiveEquivalent:
    """One finite burn reduced to a Lunaris impulsive maneuver, with evidence.

    ``initial_mass_kg``/``final_mass_kg`` are identical for the finite burn and
    for the impulse (``isp_s`` is chosen for that).  ``finite_burn_duration_s``
    states how much time the impulse compresses; the positional effect of that
    compression (finite-burn / gravity loss) is not represented.

    References
    ----------
    NIST DLMF 3.5(v): https://dlmf.nist.gov/3.5#v
    NASA-20040000363: https://ntrs.nasa.gov/citations/20040000363
    """

    burn_id: str
    t_burn_s: float
    dv_mps: tuple[float, float, float]
    frame: str
    isp_s: float
    delta_v_m_s: float
    initial_mass_kg: float
    final_mass_kg: float
    ignition_time_s: float
    finite_burn_duration_s: float
    quadrature_error_estimate_m_s: float
    bridge_model_hash: str

    def to_lunaris(self) -> Any:
        """Return ``lunaris.core.propagation.plans.ImpulsiveManeuver``."""

        plans = _lunaris_plans()
        return plans.ImpulsiveManeuver(
            t_burn_s=self.t_burn_s, dv_mps=self.dv_mps, frame=self.frame, isp_s=self.isp_s
        )


def lunaris_impulsive_equivalents(
    timeline: PropulsionTimeline,
    *,
    initial_mass_kg: float,
    directions: Mapping[str, tuple[tuple[float, float, float], str]],
    parameter_overrides: Mapping[str, Mapping[str, float]] | None = None,
) -> tuple[LunarisImpulsiveEquivalent, ...]:
    """Reduce every timeline burn to an exact-mass impulsive equivalent.

    ``directions`` maps each burn id to ``(vector, frame)`` with frame
    ``"inertial"`` or ``"ric"`` (Lunaris' radial/in-track/cross-track order).
    Burns are processed in time order; each burn starts from the previous
    burn's final mass.

    References
    ----------
    NIST DLMF 3.5(v): https://dlmf.nist.gov/3.5#v
    """

    if not isinstance(timeline, PropulsionTimeline):
        raise DomainError("A PropulsionTimeline is required.")
    if timeline.jettisons:
        raise DomainError(
            "Lunaris impulses cannot drop mass without delta-v; remove jettisons or "
            "apply them to the Lunaris state yourself between segments."
        )
    if timeline.allow_overlapping_burns:
        raise DomainError("Overlapping burns cannot be reduced to ordered impulses.")
    if not isinstance(directions, Mapping):
        raise DomainError("Directions must map burn ids to (vector, frame).")
    overrides = dict(parameter_overrides or {})
    mass = real(initial_mass_kg, "Initial mass")
    if mass <= 0.0:
        raise DomainError("Initial mass must be greater than zero.")
    results: list[LunarisImpulsiveEquivalent] = []
    for burn in timeline.burns:
        bridge = burn.bridge
        if bridge.performance_surface is not None or bridge.nozzle_exit_area_m2 is not None:
            raise DomainError(
                f"Burn {burn.burn_id!r} depends on operating conditions or ambient pressure, "
                "which an impulsive equivalent cannot represent."
            )
        if burn.burn_id not in directions:
            raise DomainError(f"No direction for burn {burn.burn_id!r}.")
        entry = directions[burn.burn_id]
        if not isinstance(entry, tuple) or len(entry) != 2:
            raise DomainError("Each direction entry must be (vector, frame).")
        vector = vector3(entry[0], f"Direction of {burn.burn_id!r}")
        frame = text(entry[1], "Direction frame").strip().lower()
        if frame not in LUNARIS_FRAMES:
            raise DomainError(f"Lunaris supports frames {LUNARIS_FRAMES}, not {frame!r}.")
        norm = sqrt(sum(value * value for value in vector))
        if norm <= 0.0:
            raise DomainError("Direction must be non-zero.")
        drivers = bridge.resolved_parameters(overrides.get(burn.burn_id))
        metrics = free_space_burn_metrics(
            bridge.artifact,
            mass,
            thrust_scale=drivers["thrust_scale"],
            mass_flow_scale=drivers["mass_flow_scale"],
        )
        if metrics.effective_exhaust_velocity_m_s is None or metrics.delta_v_m_s <= 0.0:
            raise DomainError(f"Burn {burn.burn_id!r} produces no delta-v.")
        dry = bridge.protected_dry_mass_kg
        if dry is not None and metrics.final_mass_kg < dry * (1.0 - 1e-12):
            raise DomainError(
                f"Burn {burn.burn_id!r} would cut off at its protected dry mass "
                f"({dry!r} kg) before the end of its profile; a truncated burn has no "
                "impulsive equivalent here. Use LunarisFiniteBurnDynamics instead."
            )
        ignition = burn.ignition_time_s + drivers["ignition_time_bias_s"]
        results.append(
            LunarisImpulsiveEquivalent(
                burn_id=burn.burn_id,
                t_burn_s=ignition + metrics.impulsive_equivalent_time_s,
                dv_mps=tuple(metrics.delta_v_m_s * value / norm for value in vector),  # type: ignore[arg-type]
                frame=frame,
                isp_s=metrics.effective_exhaust_velocity_m_s / STANDARD_GRAVITY_M_S2,
                delta_v_m_s=metrics.delta_v_m_s,
                initial_mass_kg=mass,
                final_mass_kg=metrics.final_mass_kg,
                ignition_time_s=ignition,
                finite_burn_duration_s=bridge.artifact.duration_s,
                quadrature_error_estimate_m_s=metrics.delta_v_error_estimate_m_s,
                bridge_model_hash=bridge.model_hash,
            )
        )
        mass = metrics.final_mass_kg
    times = [item.t_burn_s for item in results]
    if any(later <= earlier for earlier, later in pairwise(times)):
        raise DomainError("Impulse times must be strictly increasing for a Lunaris plan.")
    return tuple(results)


def to_lunaris_maneuver_plan(equivalents: tuple[LunarisImpulsiveEquivalent, ...]) -> Any:
    """Return a Lunaris ``ManeuverPlan`` from impulsive equivalents.

    References
    ----------
    NIST DLMF 3.5(v): https://dlmf.nist.gov/3.5#v
    """

    plans = _lunaris_plans()
    return plans.ManeuverPlan(tuple(item.to_lunaris() for item in equivalents))


@dataclass(frozen=True, slots=True)
class LunarisBurnWindow:
    """Effective firing window of one burn, as the finite-burn wrapper sees it.

    ``ignition_time_s`` includes the resolved ``ignition_time_bias_s`` driver.
    ``depletion_time_s`` is the exact time at which the vehicle reaches the
    burn's ``protected_dry_mass_kg`` (``None`` when it never does or when the
    burn has no protected dry mass); the burn is off from that time on.
    ``mass_at_ignition_kg`` is ``None`` when no initial mass was given.

    References
    ----------
    NASA-SP-125: https://ntrs.nasa.gov/citations/19710019929
    """

    burn_id: str
    ignition_time_s: float
    cutoff_time_s: float
    mass_at_ignition_kg: float | None
    depletion_time_s: float | None
    time_events_s: tuple[float, ...]


def _consumed_by(bridge: Any, mass_flow_scale: float, relative_time_s: float) -> float:
    """Propellant used from ignition to ``relative_time_s`` (piecewise-linear flow)."""

    used = 0.0
    for segment in bridge.artifact.segments:
        if relative_time_s <= segment.start_time_s:
            break
        duration = segment.duration_s
        tau = min(relative_time_s, segment.end_time_s) - segment.start_time_s
        q0 = segment.start_total_mass_flow_kg_s
        q1 = segment.end_total_mass_flow_kg_s
        used += q0 * tau + 0.5 * (q1 - q0) * tau * tau / duration
    return used * mass_flow_scale


def _time_to_consume(bridge: Any, mass_flow_scale: float, amount_kg: float) -> float | None:
    """Return the relative time at which ``amount_kg`` has been consumed, or ``None``."""

    remaining = amount_kg
    if remaining <= 0.0:
        return 0.0
    for segment in bridge.artifact.segments:
        duration = segment.duration_s
        q0 = segment.start_total_mass_flow_kg_s * mass_flow_scale
        q1 = segment.end_total_mass_flow_kg_s * mass_flow_scale
        available = 0.5 * (q0 + q1) * duration
        if available < remaining:
            remaining -= available
            continue
        # Solve q0 tau + (q1 - q0) tau^2 / (2 d) = remaining with the stable root.
        a = 0.5 * (q1 - q0) / duration
        root = sqrt(max(q0 * q0 + 4.0 * a * remaining, 0.0))
        denominator = q0 + root
        tau = duration if denominator <= 0.0 else 2.0 * remaining / denominator
        return segment.start_time_s + min(max(tau, 0.0), duration)
    return None


def lunaris_burn_windows(
    timeline: PropulsionTimeline,
    *,
    parameter_overrides: Mapping[str, Mapping[str, float]] | None = None,
    initial_mass_kg: float | None = None,
) -> tuple[LunarisBurnWindow, ...]:
    """Return every burn's effective window: bias-shifted events and dry-mass cut.

    For a mass-only burn without operating-condition surfaces the mass flow is
    a function of time alone, so with a known ``initial_mass_kg`` the time at
    which the vehicle reaches a burn's ``protected_dry_mass_kg`` is computed
    exactly (piecewise-quadratic mass, stable quadratic root).  Lunaris event
    functions only see the first six state elements, so this is how a
    dry-mass cut becomes an exact segment boundary.  Bridges with a protected
    dry mass therefore need ``initial_mass_kg``; overlapping timelines are
    refused because a cut of one burn changes the other's mass history.

    References
    ----------
    NASA-SP-125: https://ntrs.nasa.gov/citations/19710019929
    Orekit EventDetector 13.1.5 API:
    https://www.orekit.org/static/apidocs/org/orekit/propagation/events/EventDetector.html
    """

    if not isinstance(timeline, PropulsionTimeline):
        raise DomainError("A PropulsionTimeline is required.")
    if parameter_overrides is not None and not isinstance(parameter_overrides, Mapping):
        raise DomainError("Parameter overrides must map burn ids to driver values.")
    overrides = dict(parameter_overrides or {})
    unknown = set(overrides) - {burn.burn_id for burn in timeline.burns}
    if unknown:
        raise DomainError("Parameter overrides for unknown burns: " + ", ".join(sorted(unknown)))
    mass: float | None = None
    if initial_mass_kg is not None:
        mass = real(initial_mass_kg, "Initial mass")
        if mass <= 0.0:
            raise DomainError("Initial mass must be greater than zero.")
    protected = any(burn.bridge.protected_dry_mass_kg is not None for burn in timeline.burns)
    if protected and mass is None:
        raise DomainError(
            "A burn has protected_dry_mass_kg: pass initial_mass_kg so the dry-mass cut "
            "becomes an exact Lunaris segment boundary."
        )
    if protected and timeline.allow_overlapping_burns:
        raise DomainError("Dry-mass cuts cannot be precomputed for overlapping burns.")
    resolved = []
    for burn in timeline.burns:
        drivers = burn.bridge.resolved_parameters(overrides.get(burn.burn_id))
        resolved.append((burn.ignition_time_s + drivers["ignition_time_bias_s"], burn, drivers))
    resolved.sort(key=lambda item: item[0])
    windows: list[LunarisBurnWindow] = []
    previous_cutoff = None
    for ignition, burn, drivers in resolved:
        bridge = burn.bridge
        if (
            not timeline.allow_overlapping_burns
            and previous_cutoff is not None
            and ignition < previous_cutoff - _BOUNDARY_TOLERANCE_S
        ):
            raise DomainError(
                f"Burn {burn.burn_id!r} starts before the previous burn ends once the "
                "ignition-time biases are applied."
            )
        flow_scale = drivers["mass_flow_scale"]
        depletion = None
        cutoff_relative = bridge.artifact.duration_s
        if mass is not None and bridge.protected_dry_mass_kg is not None:
            reached = _time_to_consume(bridge, flow_scale, mass - bridge.protected_dry_mass_kg)
            if reached is not None and reached < bridge.artifact.duration_s:
                depletion = ignition + reached
                cutoff_relative = reached
        events = tuple(
            ignition + surface.relative_time_s
            for surface in bridge.event_surfaces()
            if surface.relative_time_s is not None
            and surface.relative_time_s <= cutoff_relative + _BOUNDARY_TOLERANCE_S
        ) + ((depletion,) if depletion is not None else ())
        windows.append(
            LunarisBurnWindow(
                burn_id=burn.burn_id,
                ignition_time_s=ignition,
                cutoff_time_s=ignition + cutoff_relative,
                mass_at_ignition_kg=mass,
                depletion_time_s=depletion,
                time_events_s=tuple(sorted(events)),
            )
        )
        if mass is not None:
            mass -= _consumed_by(bridge, flow_scale, cutoff_relative)
        previous_cutoff = ignition + cutoff_relative
    return tuple(windows)


def lunaris_segmentation_times(
    timeline: PropulsionTimeline,
    *,
    parameter_overrides: Mapping[str, Mapping[str, float]] | None = None,
    initial_mass_kg: float | None = None,
) -> tuple[float, ...]:
    """Return every burn discontinuity time, strictly increasing.

    Times include each burn's resolved ``ignition_time_bias_s`` and, with
    ``initial_mass_kg``, the exact dry-mass cut times (see
    :func:`lunaris_burn_windows`).  Pass the same ``parameter_overrides`` as
    to :class:`LunarisFiniteBurnDynamics`, or use its ``segmentation_plan()``.

    References
    ----------
    Orekit EventDetector 13.1.5 API:
    https://www.orekit.org/static/apidocs/org/orekit/propagation/events/EventDetector.html
    """

    windows = lunaris_burn_windows(
        timeline, parameter_overrides=parameter_overrides, initial_mass_kg=initial_mass_kg
    )
    return _merged_times(value for window in windows for value in window.time_events_s)


def _merged_times(values: Iterator[float] | Any) -> tuple[float, ...]:
    merged: list[float] = []
    for value in sorted(values):
        if not merged or value - merged[-1] > _BOUNDARY_TOLERANCE_S:
            merged.append(value)
    return tuple(merged)


def _segmentation_plan(times: tuple[float, ...]) -> Any:
    plans = _lunaris_plans()
    return plans.ManeuverPlan(
        tuple(plans.ImpulsiveManeuver(t_burn_s=value, dv_mps=(0.0, 0.0, 0.0)) for value in times)
    )


def lunaris_segmentation_plan(
    timeline: PropulsionTimeline,
    *,
    parameter_overrides: Mapping[str, Mapping[str, float]] | None = None,
    initial_mass_kg: float | None = None,
) -> Any:
    """Return a Lunaris ``ManeuverPlan`` of zero-delta-v markers at discontinuities.

    Lunaris restarts integration at every maneuver time; zero impulses with no
    ``isp_s`` change neither velocity nor mass, so the finite-burn right-hand
    side is integrated piecewise-smoothly.

    References
    ----------
    Orekit EventDetector 13.1.5 API:
    https://www.orekit.org/static/apidocs/org/orekit/propagation/events/EventDetector.html
    """

    return _segmentation_plan(
        lunaris_segmentation_times(
            timeline, parameter_overrides=parameter_overrides, initial_mass_kg=initial_mass_kg
        )
    )


class _LazyDirections(Mapping):
    """Directions computed only for the burns the timeline asks for."""

    def __init__(self, law: Callable[[str, float, Any], tuple[float, float, float]],
                 burn_ids: tuple[str, ...], time_s: float, state: Any) -> None:
        self._law = law
        self._ids = burn_ids
        self._time = time_s
        self._state = state

    def __getitem__(self, burn_id: str) -> tuple[float, float, float]:
        if burn_id not in self._ids:
            raise KeyError(burn_id)
        return self._law(burn_id, self._time, self._state)

    def __iter__(self) -> Iterator[str]:
        return iter(self._ids)

    def __len__(self) -> int:
        return len(self._ids)


def velocity_direction(burn_id: str, time_s: float, state: Any) -> tuple[float, float, float]:
    """Direction law: thrust along the inertial velocity (prograde).

    References
    ----------
    NASA-20040000363: https://ntrs.nasa.gov/citations/20040000363
    """

    return (float(state[3]), float(state[4]), float(state[5]))


class LunarisFiniteBurnDynamics:
    """Lunaris dynamics wrapper that adds finite-burn thrust and mass flow.

    ``direction_law(burn_id, t, y)`` returns the inertial thrust direction.
    Attributes other than the right-hand-side builders are delegated to the
    wrapped dynamics, so Lunaris still finds ``grav``, ``flags`` and the rest.
    Both ``build_rhs`` and Lunaris' trusted ``_build_solver_rhs`` are
    overridden; without the second, Lunaris would bypass the propulsion.

    Each Lunaris segment builds a fresh right-hand side.  The first time it is
    evaluated is the segment start, where a discontinuity takes its right-hand
    limit; any later evaluation exactly on a burn boundary is the segment end
    and takes its left-hand limit.  Propagate with :meth:`segmentation_plan`,
    which uses the same overrides and initial mass, so every discontinuity
    (bias-shifted profile events and exact dry-mass cuts) is a segment
    boundary.

    ``apply_inventory_limits=None`` (default) turns the bridge's internal
    dry-mass clamp off when ``initial_mass_kg`` is given, because the cut is
    then applied exactly at its precomputed time; a clamp evaluated inside
    Runge-Kutta stages next to the cut would bias the final mass.  Without
    ``initial_mass_kg`` the clamp stays on.

    With ``initial_mass_kg`` the wrapper also checks, at the start of every
    segment, that the state mass equals the mass predicted from the timeline
    (relative tolerance ``mass_tolerance``).  A mismatch means the propagation
    was started with a different mass or another maneuver changed the mass,
    and the precomputed dry-mass cut would be wrong, so it fails closed.

    References
    ----------
    Orekit PropulsionModel 13.1.5 API:
    https://www.orekit.org/static/apidocs/org/orekit/forces/maneuvers/propulsion/PropulsionModel.html
    NASA-SP-125: https://ntrs.nasa.gov/citations/19710019929
    """

    def __init__(
        self,
        base_dynamics: Any,
        timeline: PropulsionTimeline,
        direction_law: Callable[[str, float, Any], tuple[float, float, float]] = velocity_direction,
        *,
        parameter_overrides: Mapping[str, Mapping[str, float]] | None = None,
        apply_inventory_limits: bool | None = None,
        initial_mass_kg: float | None = None,
        mass_tolerance: float = 1e-9,
    ) -> None:
        if not isinstance(timeline, PropulsionTimeline):
            raise DomainError("A PropulsionTimeline is required.")
        if not callable(direction_law):
            raise DomainError("direction_law must be callable(burn_id, t, y).")
        if not hasattr(base_dynamics, "build_rhs"):
            raise DomainError("Base dynamics must provide build_rhs().")
        if apply_inventory_limits is not None and not isinstance(apply_inventory_limits, bool):
            raise DomainError("apply_inventory_limits must be a boolean or None.")
        for burn in timeline.burns:
            bridge = burn.bridge
            if bridge.stream_bindings or bridge.state_constraints:
                raise DomainError(
                    f"Burn {burn.burn_id!r} binds tank states; Lunaris carries only vehicle "
                    "mass. Build the bridge with stream_bindings=() and "
                    "require_complete_stream_binding=False (protected_dry_mass_kg still works)."
                )
            if bridge.performance_surface is not None or bridge.nozzle_exit_area_m2 is not None:
                raise DomainError(
                    f"Burn {burn.burn_id!r} needs operating conditions or ambient pressure, "
                    "which Lunaris does not supply."
                )
        if timeline.jettisons:
            raise DomainError("Apply jettisons to the Lunaris state between segments.")
        tolerance = real(mass_tolerance, "Mass tolerance")
        if not 0.0 < tolerance < 1.0:
            raise DomainError("Mass tolerance must be in (0, 1).")
        windows = lunaris_burn_windows(
            timeline, parameter_overrides=parameter_overrides, initial_mass_kg=initial_mass_kg
        )
        self._base = base_dynamics
        self._timeline = timeline
        self._law = direction_law
        self._overrides = parameter_overrides
        # With a known initial mass the dry-mass cut is an exact boundary, and
        # the bridge's internal clamp would only bias the integrator's stages.
        self._limits = (
            initial_mass_kg is None if apply_inventory_limits is None else apply_inventory_limits
        )
        self._ids = tuple(burn.burn_id for burn in timeline.burns)
        self._windows = windows
        self._boundaries = _merged_times(
            value for window in windows for value in window.time_events_s
        )
        self._initial_mass = None if initial_mass_kg is None else real(initial_mass_kg, "Mass")
        self._mass_tolerance = tolerance
        self._burns = {burn.burn_id: burn for burn in timeline.burns}

    def __getattr__(self, name: str) -> Any:
        if name.startswith("__") or name == "_base":
            raise AttributeError(name)
        return getattr(self._base, name)

    @property
    def burn_windows(self) -> tuple[LunarisBurnWindow, ...]:
        """Effective burn windows used by this wrapper."""

        return self._windows

    @property
    def segmentation_times(self) -> tuple[float, ...]:
        """Every discontinuity time of this wrapper, strictly increasing."""

        return self._boundaries

    def segmentation_plan(self) -> Any:
        """Lunaris ``ManeuverPlan`` of zero-delta-v markers at :attr:`segmentation_times`."""

        return _segmentation_plan(self._boundaries)

    def predicted_mass_kg(self, time_s: float) -> float:
        """Vehicle mass implied by the timeline at ``time_s`` (needs ``initial_mass_kg``)."""

        if self._initial_mass is None:
            raise DomainError("Predicted mass needs initial_mass_kg.")
        time_s = real(time_s, "Time")
        mass = self._initial_mass
        for window in self._windows:
            if time_s <= window.ignition_time_s:
                break
            burn = self._burns[window.burn_id]
            drivers = burn.bridge.resolved_parameters(
                (self._overrides or {}).get(window.burn_id)
            )
            relative = min(time_s, window.cutoff_time_s) - window.ignition_time_s
            mass -= _consumed_by(burn.bridge, drivers["mass_flow_scale"], relative)
        return mass

    def _stopped(self, time_s: float, side: str) -> frozenset[str]:
        stopped = set()
        for window in self._windows:
            cut = window.depletion_time_s
            if cut is None:
                continue
            if time_s > cut + _BOUNDARY_TOLERANCE_S or (
                abs(time_s - cut) <= _BOUNDARY_TOLERANCE_S and side == "right"
            ):
                stopped.add(window.burn_id)
        return frozenset(stopped)

    def _wrap(self, base_rhs: Callable[[float, Any], Any]) -> Callable[[float, Any], Any]:
        first_time: list[float] = []
        boundaries = self._boundaries
        timeline = self._timeline

        def rhs(t: float, y: Any) -> Any:
            if len(y) != 7:
                raise DomainError("Finite burns need the 7-element Lunaris state (mass at y[6]).")
            out = base_rhs(t, y)
            time_s = float(t)
            if not first_time:
                first_time.append(time_s)
                if self._initial_mass is not None:
                    expected = self.predicted_mass_kg(time_s)
                    if abs(float(y[6]) - expected) > self._mass_tolerance * self._initial_mass:
                        raise DomainError(
                            f"State mass {float(y[6])!r} kg at t={time_s!r} s differs from the "
                            f"timeline prediction {expected!r} kg; the dry-mass cut would be "
                            "wrong. Start from initial_mass_kg and add no other mass changes."
                        )
            on_boundary = any(abs(time_s - value) <= _BOUNDARY_TOLERANCE_S for value in boundaries)
            side = (
                "left"
                if on_boundary and abs(time_s - first_time[0]) > _BOUNDARY_TOLERANCE_S
                else "right"
            )
            value = timeline.evaluate(
                time_s=time_s,
                vehicle_mass_kg=float(y[6]),
                directions=_LazyDirections(self._law, self._ids, time_s, y),
                parameter_overrides=self._overrides,
                side=side,
                stopped_burns=self._stopped(time_s, side),
                apply_inventory_limits=self._limits,
            )
            out[3] += value.acceleration_m_s2[0]
            out[4] += value.acceleration_m_s2[1]
            out[5] += value.acceleration_m_s2[2]
            out[6] += value.mass_derivative_kg_s
            return out

        return rhs

    def build_rhs(self) -> Callable[[float, Any], Any]:
        """Return the wrapped right-hand side (validated base builder)."""

        return self._wrap(self._base.build_rhs())

    def _build_solver_rhs(self) -> Callable[[float, Any], Any]:
        """Return the wrapped trusted right-hand side used by Lunaris' propagator."""

        trusted = getattr(self._base, "_build_solver_rhs", None)
        return self._wrap(trusted() if callable(trusted) else self._base.build_rhs())
