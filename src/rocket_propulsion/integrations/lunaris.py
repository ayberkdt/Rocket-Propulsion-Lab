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
from fractions import Fraction
from importlib import import_module
from itertools import pairwise
from typing import Any

from rocket_propulsion.core.errors import DomainError, FeatureUnavailableError
from rocket_propulsion.core.validation import real, text, vector3
from rocket_propulsion.propulsion.burns import (
    PropulsionTimeline,
    free_space_burn_metrics,
)
from rocket_propulsion.propulsion.burns.numeric_outputs import _finite_rational, _fraction
from rocket_propulsion.propulsion.burns.propagator import INVENTORY_LIMIT_TOLERANCE
from rocket_propulsion.propulsion.burns.time_coordinates import (
    _clock_table,
    _epoch_profile_time,
    _profile_boundaries,
)
from rocket_propulsion.propulsion.performance import STANDARD_GRAVITY_M_S2

from .direction import _directed_vector
from .native_mass import _NATIVE_MASS_POLICY, _consumption, _cutoff, _scale, _upper_bound

LUNARIS_FRAMES = ("inertial", "ric")


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
    mass = real(initial_mass_kg, "Initial mass")
    if mass <= 0.0:
        raise DomainError("Initial mass must be greater than zero.")
    # Same override validation, bias-shifted ordering and overlap refusal as
    # the finite-burn path, so both paths accept exactly the same inputs.
    windows = lunaris_burn_windows(
        timeline, parameter_overrides=parameter_overrides, initial_mass_kg=mass
    )
    overrides = dict(parameter_overrides or {})
    by_id = {burn.burn_id: burn for burn in timeline.burns}
    results: list[LunarisImpulsiveEquivalent] = []
    for window in windows:
        burn = by_id[window.burn_id]
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
                dv_mps=_directed_vector(vector, metrics.delta_v_m_s, label="Lunaris impulse delta-v"),
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
    ``depletion_time_s`` is a guarded, represented cutoff for the burn's
    protected dry mass and declared inventory margin (``None`` when no cut
    is needed). The burn is off from that date on. The reported numerical
    shortfall bound measures unused mass relative to that margin-adjusted
    target; it is not a physical uncertainty bound.
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
    cutoff_profile_time_ratio: tuple[int, int] | None = None
    mass_at_cutoff_kg: float | None = None
    depletion_mass_shortfall_bound_kg: float | None = None
    cutoff_numerics: str = 'unspecified'

    def __post_init__(self):
        text(self.burn_id, 'Native burn identifier')
        for name in ('ignition_time_s','cutoff_time_s'):
            object.__setattr__(self,name,real(getattr(self,name),name))
        if self.cutoff_time_s < self.ignition_time_s:
            raise DomainError('Native cutoff cannot precede ignition.')
        for name in ('mass_at_ignition_kg','mass_at_cutoff_kg','depletion_time_s',
                     'depletion_mass_shortfall_bound_kg'):
            value=getattr(self,name)
            if value is not None:
                object.__setattr__(self,name,real(value,name))
        for name in ('mass_at_ignition_kg','mass_at_cutoff_kg'):
            if getattr(self,name) is not None and getattr(self,name) <= 0:
                raise DomainError('Native window masses must be positive.')
        if self.depletion_time_s is not None and self.depletion_time_s != self.cutoff_time_s:
            raise DomainError('A native depletion date must be the cutoff date.')
        if self.depletion_mass_shortfall_bound_kg is not None and (
                self.depletion_mass_shortfall_bound_kg < 0 or self.depletion_time_s is None):
            raise DomainError('A native depletion deficit requires a nonnegative bound and cut.')
        if isinstance(self.time_events_s,(str,bytes,Mapping)):
            raise DomainError('Native event dates must be a sequence of real scalars.')
        try:
            dates=tuple(real(t,'Native event date') for t in self.time_events_s)
        except TypeError as error:
            raise DomainError('Native event dates must be a sequence.') from error
        if not dates or dates[0]!=self.ignition_time_s or dates[-1]!=self.cutoff_time_s or any(
                b<=a for a,b in pairwise(dates)):
            raise DomainError('Native event dates must include both endpoints and increase strictly.')
        object.__setattr__(self,'time_events_s',dates)
        ratio=self.cutoff_profile_time_ratio
        if ratio is not None:
            if not isinstance(ratio,(tuple,list)) or len(ratio)!=2 or any(type(n) is not int for n in ratio):
                raise DomainError('Native cutoff coordinate requires an integer ratio.')
            if ratio[0]<0 or ratio[1]<=0:
                raise DomainError('Native cutoff coordinate ratio must be nonnegative.')
            exact=Fraction(*ratio)
            object.__setattr__(self,'cutoff_profile_time_ratio',(exact.numerator,exact.denominator))
        if self.cutoff_numerics not in ('unspecified',_NATIVE_MASS_POLICY):
            raise DomainError('Unsupported native cutoff numerical policy.')


def _consumed_by(bridge: Any, mass_flow_scale: float, relative_time_s: float) -> float:
    """Propellant used from ignition to ``relative_time_s`` (piecewise-linear flow)."""

    return _finite_rational(_consumption(bridge.artifact.segments, _scale(mass_flow_scale),
                                         _fraction(relative_time_s)), 'Native consumed mass')


def _time_to_consume(bridge: Any, mass_flow_scale: float, amount_kg: float) -> float | None:
    """Return the relative time at which ``amount_kg`` has been consumed, or ``None``.

    The returned scalar is the greatest safe represented relative time, with
    decisions made against the exact mass integral. The existing inventory
    margin remains a caller/model policy, not measurement uncertainty.
    """

    # The bridge's clamp treats a mass within INVENTORY_LIMIT_TOLERANCE of the
    # limit as reached; the cut uses the same margin, otherwise rounding can
    # carry it past a zero-flow plateau that already sits at the limit.
    target = _fraction(amount_kg)-_fraction(INVENTORY_LIMIT_TOLERANCE)
    cut = _cutoff(bridge.artifact.segments, _scale(mass_flow_scale), target, _fraction(0.))
    return None if cut is None else cut.date


def lunaris_burn_windows(
    timeline: PropulsionTimeline,
    *,
    parameter_overrides: Mapping[str, Mapping[str, float]] | None = None,
    initial_mass_kg: float | None = None,
) -> tuple[LunarisBurnWindow, ...]:
    """Return every burn's effective window: bias-shifted events and dry-mass cut.

    For a mass-only burn without operating-condition surfaces the mass flow is
    a function of time alone. With a known ``initial_mass_kg``, exact supplied-
    binary mass integrals select a safe date on the represented epoch clock.
    The date can precede the analytical crossing; the window reports an
    outward-rounded unused-mass bound. Lunaris event functions only see the
    first six state elements, so the cut is exposed as a segment boundary.
    Bridges with a protected
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
    initial = None
    if initial_mass_kg is not None:
        initial = _fraction(real(initial_mass_kg, "Initial mass"))
        if initial <= 0:
            raise DomainError("Initial mass must be greater than zero.")
    protected = any(burn.bridge.protected_dry_mass_kg is not None for burn in timeline.burns)
    if protected and initial is None:
        raise DomainError("Protected native burns require initial_mass_kg.")
    if protected and timeline.allow_overlapping_burns:
        raise DomainError("Dry-mass cuts cannot be precomputed for overlapping burns.")
    if initial is not None and timeline.jettisons:
        raise DomainError("Native mass prediction requires jettisons to be applied between segments.")
    resolved = []
    for burn in timeline.burns:
        drivers = burn.bridge.resolved_parameters(overrides.get(burn.burn_id))
        _scale(drivers['mass_flow_scale'])
        bridge = burn.bridge
        if initial is not None and (bridge.performance_surface is not None
                                    or bridge.stream_performance_surface is not None):
            raise DomainError("Native mass prediction needs time-only flow, without a performance surface.")
        offset = _fraction(burn.ignition_time_s)+_fraction(drivers["ignition_time_bias_s"])
        resolved.append((offset, burn, drivers))
    resolved.sort(key=lambda row: row[0])
    contexts = []

    def consumed_at(context, date):
        _offset, burn, drivers, horizon, _, _, _, _ = context
        point = _epoch_profile_time(burn.bridge.artifact.segments, date, burn.ignition_time_s,
            drivers["ignition_time_bias_s"], horizon=horizon, allow_aliases=False)
        return _consumption(burn.bridge.artifact.segments, _fraction(drivers["mass_flow_scale"]),
                            min(max(point, Fraction(0)), horizon))

    def mass_at(date, active):
        if initial is None:
            return None
        return initial-sum((consumed_at(row,date) for row in active),Fraction(0))

    for offset, burn, drivers in resolved:
        bridge = burn.bridge
        ignition = _finite_rational(offset, "Native ignition time")
        if not timeline.allow_overlapping_burns and contexts and ignition < contexts[-1][5][-1]:
            raise DomainError("Biased native burns overlap on the represented clock.")
        mass = mass_at(ignition, contexts)
        if mass is not None and mass <= 0:
            raise DomainError("Native mass history reaches nonpositive vehicle mass.")
        horizon = _fraction(bridge.artifact.duration_s)
        depletion = None
        deficit = None
        if mass is not None and bridge.protected_dry_mass_kg is not None:
            target = mass-_fraction(bridge.protected_dry_mass_kg)-_fraction(INVENTORY_LIMIT_TOLERANCE)
            cut = _cutoff(bridge.artifact.segments, _fraction(drivers["mass_flow_scale"]), target, offset)
            reached = None if cut is None else cut.profile_time
            if reached is not None and reached < bridge.artifact.duration_s:
                horizon = reached
                depletion = cut.date
                deficit = cut.shortfall
        boundaries = _profile_boundaries(bridge.artifact.segments, horizon)
        dates, _ = _clock_table(boundaries,burn.ignition_time_s,drivers["ignition_time_bias_s"])
        if depletion is not None and dates[-1] != depletion:
            raise DomainError("Native cutoff and active prefix clock do not agree.")
        contexts.append((offset,burn,drivers,horizon,ignition,dates,depletion,deficit))
    if initial is not None:
        final = initial-sum((_consumption(row[1].bridge.artifact.segments,
                            _fraction(row[2]["mass_flow_scale"]),row[3]) for row in contexts),Fraction(0))
        if final <= 0:
            raise DomainError("Native mass history reaches nonpositive vehicle mass.")
    windows=[]
    for offset,burn,drivers,horizon,ignition,dates,depletion,deficit in contexts:
        before=mass_at(ignition,contexts)
        after=mass_at(dates[-1],contexts)
        windows.append(LunarisBurnWindow(
            burn_id=burn.burn_id,ignition_time_s=ignition,cutoff_time_s=dates[-1],
            mass_at_ignition_kg=None if before is None else _finite_rational(before,"Native ignition mass"),
            depletion_time_s=depletion,time_events_s=dates,
            cutoff_numerics=_NATIVE_MASS_POLICY,
            cutoff_profile_time_ratio=(horizon.numerator,horizon.denominator),
            mass_at_cutoff_kg=None if after is None else _finite_rational(after,"Native cutoff mass"),
            depletion_mass_shortfall_bound_kg=None if deficit is None else _upper_bound(deficit)))
    return tuple(windows)


def lunaris_segmentation_times(
    timeline: PropulsionTimeline,
    *,
    parameter_overrides: Mapping[str, Mapping[str, float]] | None = None,
    initial_mass_kg: float | None = None,
) -> tuple[float, ...]:
    """Return every burn discontinuity time, strictly increasing.

    Times include each burn's resolved ``ignition_time_bias_s`` and, with
    ``initial_mass_kg``, the guarded dry-mass cut dates (see
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
    """Deduplicate identical dates; never merge distinct representable events."""
    return tuple(sorted(set(values)))


def _window(start_s: float | None, end_s: float | None) -> tuple[float, float]:
    start = float("-inf") if start_s is None else real(start_s, "Window start")
    end = float("inf") if end_s is None else real(end_s, "Window end")
    if end < start:
        raise DomainError("Window end must not precede its start.")
    return start, end


def _segmentation_plan(
    times: tuple[float, ...], start_s: float | None = None, end_s: float | None = None
) -> Any:
    start, end = _window(start_s, end_s)
    plans = _lunaris_plans()
    return plans.ManeuverPlan(
        tuple(
            plans.ImpulsiveManeuver(t_burn_s=value, dv_mps=(0.0, 0.0, 0.0))
            for value in times
            if start <= value <= end
        )
    )


def lunaris_segmentation_plan(
    timeline: PropulsionTimeline,
    *,
    parameter_overrides: Mapping[str, Mapping[str, float]] | None = None,
    initial_mass_kg: float | None = None,
    start_s: float | None = None,
    end_s: float | None = None,
) -> Any:
    """Return a Lunaris ``ManeuverPlan`` of zero-delta-v markers at discontinuities.

    Lunaris restarts integration at every maneuver time; zero impulses with no
    ``isp_s`` change neither velocity nor mass, so the finite-burn right-hand
    side is integrated piecewise-smoothly.  Lunaris refuses maneuvers outside
    the propagation window ``[t0, t0 + duration]``; pass that window as
    ``start_s``/``end_s`` to propagate only part of a timeline.

    References
    ----------
    Orekit EventDetector 13.1.5 API:
    https://www.orekit.org/static/apidocs/org/orekit/propagation/events/EventDetector.html
    """

    return _segmentation_plan(
        lunaris_segmentation_times(
            timeline, parameter_overrides=parameter_overrides, initial_mass_kg=initial_mass_kg
        ),
        start_s,
        end_s,
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
    (bias-shifted profile events and guarded dry-mass cuts) is a segment
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
            if (bridge.performance_surface is not None or bridge.stream_performance_surface is not None
                    or bridge.nozzle_exit_area_m2 is not None):
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
        self._overrides = {burn.burn_id:burn.bridge.resolved_parameters(
            (parameter_overrides or {}).get(burn.burn_id)) for burn in timeline.burns}
        # With a known initial mass the dry-mass cut is an exact boundary, and
        # the bridge's internal clamp would only bias the integrator's stages.
        self._limits = (
            initial_mass_kg is None if apply_inventory_limits is None else apply_inventory_limits
        )
        self._ids = tuple(burn.burn_id for burn in timeline.burns)
        self._windows = windows
        self._horizons = {window.burn_id:Fraction(*window.cutoff_profile_time_ratio) for window in windows}
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

    def segmentation_plan(self, start_s: float | None = None, end_s: float | None = None) -> Any:
        """Lunaris ``ManeuverPlan`` of zero-delta-v markers at :attr:`segmentation_times`.

        ``start_s``/``end_s`` keep only markers inside a propagation window
        (Lunaris refuses maneuvers outside ``[t0, t0 + duration]``).

        References
        ----------
        Orekit EventDetector 13.1.5 API:
        https://www.orekit.org/static/apidocs/org/orekit/propagation/events/EventDetector.html
        """

        return _segmentation_plan(self._boundaries, start_s, end_s)

    def predicted_mass_kg(self, time_s: float) -> float:
        """Vehicle mass implied by the timeline at ``time_s`` (needs ``initial_mass_kg``)."""

        return _finite_rational(self._predicted_mass(time_s), 'Native predicted vehicle mass')

    def _predicted_mass(self, time_s):
        """Keep the shared mass ledger exact through initial-state consistency checks."""

        if self._initial_mass is None:
            raise DomainError("Predicted mass needs initial_mass_kg.")
        time_s = real(time_s, "Time")
        mass = _fraction(self._initial_mass)
        for window in self._windows:
            if time_s <= window.ignition_time_s:
                break
            burn = self._burns[window.burn_id]
            drivers = burn.bridge.resolved_parameters(
                (self._overrides or {}).get(window.burn_id)
            )
            relative_exact = _epoch_profile_time(burn.bridge.artifact.segments,
                min(time_s, window.cutoff_time_s), burn.ignition_time_s,
                drivers['ignition_time_bias_s'],horizon=self._horizons[window.burn_id],allow_aliases=False)
            mass -= _consumption(burn.bridge.artifact.segments,_fraction(drivers["mass_flow_scale"]),
                                 min(max(relative_exact,Fraction(0)),self._horizons[window.burn_id]))
        return mass

    def _stopped(self, time_s: float, side: str) -> frozenset[str]:
        stopped = set()
        for window in self._windows:
            cut = window.depletion_time_s
            if cut is None:
                continue
            if time_s > cut or (
                time_s == cut and side == "right"
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
            time_s = real(t,'Native RHS time')
            mass_value = real(y[6],'Native RHS vehicle mass')
            if mass_value <= 0:
                raise DomainError('Native RHS vehicle mass must be positive.')
            if not first_time and self._initial_mass is not None:
                exact_expected = self._predicted_mass(time_s)
                expected = _finite_rational(exact_expected, 'Native predicted vehicle mass')
                if abs(_fraction(mass_value)-exact_expected) > (
                        _fraction(self._mass_tolerance)*_fraction(self._initial_mass)):
                    raise DomainError(
                        f"State mass {float(y[6])!r} kg at t={time_s!r} s differs from the "
                        f"timeline prediction {expected!r} kg; the dry-mass cut would be "
                        "wrong. Start from initial_mass_kg and add no other mass changes."
                    )
            on_boundary = time_s in boundaries
            side = (
                "left"
                if on_boundary
                and first_time
                and time_s != first_time[0]
                else "right"
            )
            value = timeline._evaluate(
                time_s=time_s,
                vehicle_mass_kg=mass_value,
                directions=_LazyDirections(self._law, self._ids, time_s, y),
                parameter_overrides=self._overrides,
                side=side,
                stopped_burns=self._stopped(time_s, side),
                apply_inventory_limits=self._limits,
                _clock_horizons=self._horizons,
                _canonical_epoch=True,
            )
            out = base_rhs(time_s, y)
            out[3] += value.acceleration_m_s2[0]
            out[4] += value.acceleration_m_s2[1]
            out[5] += value.acceleration_m_s2[2]
            out[6] += value.mass_derivative_kg_s
            # A rejected/failed evaluation must not initialize the segment's
            # clock or suppress the next initial-mass consistency check.
            if not first_time:
                first_time.append(time_s)
            return out

        return rhs

    def build_rhs(self) -> Callable[[float, Any], Any]:
        """Return the wrapped right-hand side (validated base builder)."""

        return self._wrap(self._base.build_rhs())

    def _build_solver_rhs(self) -> Callable[[float, Any], Any]:
        """Return the wrapped trusted right-hand side used by Lunaris' propagator."""

        trusted = getattr(self._base, "_build_solver_rhs", None)
        return self._wrap(trusted() if callable(trusted) else self._base.build_rhs())
