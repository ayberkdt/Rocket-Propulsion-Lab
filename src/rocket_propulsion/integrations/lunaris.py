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


def lunaris_segmentation_times(timeline: PropulsionTimeline) -> tuple[float, ...]:
    """Return every nominal burn discontinuity time, strictly increasing.

    References
    ----------
    Orekit EventDetector 13.1.5 API:
    https://www.orekit.org/static/apidocs/org/orekit/propagation/events/EventDetector.html
    """

    if not isinstance(timeline, PropulsionTimeline):
        raise DomainError("A PropulsionTimeline is required.")
    times = sorted({event.time_s for event in timeline.events() if event.time_s is not None})
    merged: list[float] = []
    for value in times:
        if not merged or value - merged[-1] > _BOUNDARY_TOLERANCE_S:
            merged.append(value)
    return tuple(merged)


def lunaris_segmentation_plan(timeline: PropulsionTimeline) -> Any:
    """Return a Lunaris ``ManeuverPlan`` of zero-delta-v markers at discontinuities.

    Lunaris restarts integration at every maneuver time; zero impulses with no
    ``isp_s`` change neither velocity nor mass, so the finite-burn right-hand
    side is integrated piecewise-smoothly.

    References
    ----------
    Orekit EventDetector 13.1.5 API:
    https://www.orekit.org/static/apidocs/org/orekit/propagation/events/EventDetector.html
    """

    plans = _lunaris_plans()
    return plans.ManeuverPlan(
        tuple(
            plans.ImpulsiveManeuver(t_burn_s=value, dv_mps=(0.0, 0.0, 0.0))
            for value in lunaris_segmentation_times(timeline)
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
    and takes its left-hand limit.  Use :func:`lunaris_segmentation_plan` so
    every discontinuity is a segment boundary.

    References
    ----------
    Orekit PropulsionModel 13.1.5 API:
    https://www.orekit.org/static/apidocs/org/orekit/forces/maneuvers/propulsion/PropulsionModel.html
    """

    def __init__(
        self,
        base_dynamics: Any,
        timeline: PropulsionTimeline,
        direction_law: Callable[[str, float, Any], tuple[float, float, float]] = velocity_direction,
        *,
        parameter_overrides: Mapping[str, Mapping[str, float]] | None = None,
        apply_inventory_limits: bool = True,
    ) -> None:
        if not isinstance(timeline, PropulsionTimeline):
            raise DomainError("A PropulsionTimeline is required.")
        if not callable(direction_law):
            raise DomainError("direction_law must be callable(burn_id, t, y).")
        if not hasattr(base_dynamics, "build_rhs"):
            raise DomainError("Base dynamics must provide build_rhs().")
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
        self._base = base_dynamics
        self._timeline = timeline
        self._law = direction_law
        self._overrides = parameter_overrides
        self._limits = apply_inventory_limits
        self._ids = tuple(burn.burn_id for burn in timeline.burns)
        self._boundaries = lunaris_segmentation_times(timeline)

    def __getattr__(self, name: str) -> Any:
        if name.startswith("__"):
            raise AttributeError(name)
        return getattr(self._base, name)

    def _wrap(self, base_rhs: Callable[[float, Any], Any]) -> Callable[[float, Any], Any]:
        first_time: list[float] = []
        boundaries = self._boundaries
        timeline = self._timeline

        def rhs(t: float, y: Any) -> Any:
            out = base_rhs(t, y)
            if len(y) != 7:
                raise DomainError("Finite burns need the 7-element Lunaris state (mass at y[6]).")
            time_s = float(t)
            if not first_time:
                first_time.append(time_s)
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
