"""Multi-burn propulsion timelines with discrete mass jettison.

A single :class:`~.propagator.PropagatorPropulsionBridge` represents one burn.
A mission arc usually needs several: an ascent stage, a coast, a separation,
an upper-stage burn, a later restart.  :class:`PropulsionTimeline` places
burns on one timeline clock, checks that they do not overlap unless that is
declared, adds stage-separation / fairing jettison events as discrete mass
drops, and evaluates all active burns in one call.

The timeline computes no orbit.  Directions, ambient pressure and operating
conditions are supplied per burn by the consumer, which also owns event
location and the mode state that records which burns were stopped at an
inventory root (passed back as ``stopped_burns``).  The event and force
contract follows an event-driven propulsion model [1]_.  Jettisoned masses
are discrete vehicle-mass decrements, as in staging analysis [2]_.

Absolute timeline times are binary64 readings; each burn maps them to exact
profile-relative times with rational arithmetic [3]_, so a reading equal to a
canonical boundary date resolves to that boundary exactly.

References
----------
.. [1] OREKIT-13.1.5-PROPULSION-MODEL. Orekit 13.1.5, "PropulsionModel /
       EventDetector" API.
.. [2] SUTTON-BIBLARZ-2017. G. P. Sutton and O. Biblarz, "Rocket Propulsion
       Elements", 9th ed., Wiley, 2017.
.. [3] PYTHON-3.12-FRACTION. Python 3.12, "fractions": float inputs retain
       their exact binary value as a rational,
       https://docs.python.org/3.12/library/fractions.html.
"""

from __future__ import annotations

import json
from collections.abc import Mapping
from dataclasses import dataclass
from enum import Enum

from rocket_propulsion.core.errors import DomainError, InputError
from rocket_propulsion.core.validation import read_json_document, real, real_mapping, records, text

from .exchange_documents import _timeline_document
from .numeric_outputs import _finite_rational, _fraction
from .propagator import (
    BurnEventSurface,
    PropagatorPropulsionBridge,
    PropulsionDynamicsEvaluation,
)
from .serialization import _require_numerical_policies, canonical_sha256
from .time_coordinates import _TIME_POLICY, _epoch_profile_time, _epoch_table

PROPULSION_TIMELINE_SCHEMA = "rocket_propulsion_propulsion_timeline_v1"




@dataclass(frozen=True, slots=True)
class ScheduledBurn:
    """One propulsion bridge placed on the timeline clock.

    Parameters
    ----------
    burn_id : str
        Non-empty identifier; must not contain ``"/"`` or equal the reserved
        event namespace ``"jettison"``.
    bridge : PropagatorPropulsionBridge
        Burn model evaluated on its own profile-relative clock.
    ignition_time_s : float
        Timeline time of profile time zero, in s.  The bridge's
        ``ignition_time_bias_s`` driver shifts the burn relative to this
        nominal time.

    Raises
    ------
    DomainError
        If the identifier is blank or reserved, the bridge has the wrong type,
        or the ignition time is not a finite real number.

    References
    ----------
    .. [1] OREKIT-13.1.5-PROPULSION-MODEL. Orekit 13.1.5, "PropulsionModel /
           EventDetector" API.
    """

    burn_id: str
    bridge: PropagatorPropulsionBridge
    ignition_time_s: float

    def __post_init__(self) -> None:
        text(self.burn_id, "Burn identifier")
        if "/" in self.burn_id or self.burn_id == "jettison":
            raise DomainError(
                "Burn identifiers must not contain '/' or be 'jettison' (reserved event namespace)."
            )
        if not isinstance(self.bridge, PropagatorPropulsionBridge):
            raise DomainError("A scheduled burn requires a PropagatorPropulsionBridge.")
        object.__setattr__(
            self, "ignition_time_s", real(self.ignition_time_s, "Burn ignition time")
        )

    @property
    def cutoff_time_s(self) -> float:
        """Nominal timeline time of the profile end, in s.

        The sum ``ignition_time_s + duration`` is formed exactly and rounded
        once [1]_.

        Raises
        ------
        DomainError
            If the rounded sum is not finite.

        References
        ----------
        .. [1] PYTHON-3.12-FRACTION. Python 3.12, "fractions": float inputs retain
               their exact binary value as a rational,
               https://docs.python.org/3.12/library/fractions.html.
        """

        return _finite_rational(_fraction(self.ignition_time_s)+_fraction(self.bridge.artifact.duration_s),
                                'Scheduled cutoff time')


@dataclass(frozen=True, slots=True)
class MassJettison:
    """A discrete mass drop (stage, fairing, adapter) at a timeline time.

    Parameters
    ----------
    event_id : str
        Non-empty identifier without ``"/"``.
    time_s : float
        Timeline time of the drop, in s.
    mass_kg : float
        Dropped mass, in kg.  Finite and positive.
    label : str, optional
        Non-empty descriptive label.

    Raises
    ------
    DomainError
        If a field is blank, contains ``"/"``, or is outside its domain.

    Notes
    -----
    The vehicle mass decreases instantaneously by ``mass_kg``, the staging
    idealization of [1]_; separation dynamics and impulse are not modelled.

    References
    ----------
    .. [1] SUTTON-BIBLARZ-2017. G. P. Sutton and O. Biblarz, "Rocket Propulsion
           Elements", 9th ed., Wiley, 2017.
    """

    event_id: str
    time_s: float
    mass_kg: float
    label: str = "jettison"

    def __post_init__(self) -> None:
        text(self.event_id, "Jettison event identifier")
        if "/" in self.event_id:
            raise DomainError("Jettison identifiers must not contain '/'.")
        text(self.label, "Jettison label")
        object.__setattr__(self, "time_s", real(self.time_s, "Jettison time"))
        object.__setattr__(self, "mass_kg", real(self.mass_kg, "Jettison mass"))
        if self.mass_kg <= 0.0:
            raise DomainError("Jettison mass must be greater than zero.")


class TimelineEventKind(str, Enum):
    """Kinds of timeline roots.

    ``BURN_TIME`` is a burn boundary at a known time, ``BURN_STATE`` a burn
    root on a state variable (dry mass, tank reserve), and ``JETTISON`` a
    discrete mass drop.  The classification follows an event-detector
    contract [1]_.

    References
    ----------
    .. [1] OREKIT-13.1.5-PROPULSION-MODEL. Orekit 13.1.5, "PropulsionModel /
           EventDetector" API.
    """

    BURN_TIME = "burn_time"
    BURN_STATE = "burn_state"
    JETTISON = "jettison"


@dataclass(frozen=True, slots=True)
class TimelineEvent:
    """One timeline root with its owner and meaning.

    Parameters
    ----------
    event_id : str
        Namespaced identifier, ``"<burn_id>/<surface>"`` or
        ``"jettison/<event_id>"``.
    kind : TimelineEventKind
        Event classification.
    owner_id : str
        Identifier of the owning burn or jettison.
    time_s : float or None, optional
        Timeline time for time events (burn start, boundary or stop at
        nominal ignition; jettison); ``None`` for state events.
    surface : BurnEventSurface or None, optional
        The burn's own root surface.  When a state event fires, the consumer
        adds the owner to ``stopped_burns``.
    stops_burn : bool, optional
        Whether the event ends the owning burn.
    jettison_mass_kg : float, optional
        Dropped mass for jettison events, in kg.

    References
    ----------
    .. [1] OREKIT-13.1.5-PROPULSION-MODEL. Orekit 13.1.5, "PropulsionModel /
           EventDetector" API.
    """

    event_id: str
    kind: TimelineEventKind
    owner_id: str
    time_s: float | None = None
    surface: BurnEventSurface | None = None
    stops_burn: bool = False
    jettison_mass_kg: float = 0.0


@dataclass(frozen=True, slots=True)
class TimelineEvaluation:
    """Summed dynamics of every active burn at one timeline time.

    Parameters
    ----------
    time_s : float
        Timeline time of the evaluation, in s.
    active_burns : tuple of str
        Burns that produced thrust at ``time_s``.
    acceleration_m_s2 : tuple of float
        Sum of the burn accelerations, in m/s**2.
    mass_derivative_kg_s : float
        Sum of the vehicle mass derivatives, in kg/s (negative when firing).
    additional_state_derivatives : tuple of (str, float)
        Summed tank-state derivatives, sorted by state name.
    acceleration_wrt_mass : tuple of float
        Sum of the acceleration partials with respect to vehicle mass.
    burn_evaluations : tuple of (str, PropulsionDynamicsEvaluation)
        Per-burn evaluations, including inactive or inhibited burns inside
        their profile window.

    Notes
    -----
    Forces and mass flows of concurrent burns add linearly [1]_.

    References
    ----------
    .. [1] OREKIT-13.1.5-PROPULSION-MODEL. Orekit 13.1.5, "PropulsionModel /
           EventDetector" API.
    """

    time_s: float
    active_burns: tuple[str, ...]
    acceleration_m_s2: tuple[float, float, float]
    mass_derivative_kg_s: float
    additional_state_derivatives: tuple[tuple[str, float], ...]
    acceleration_wrt_mass: tuple[float, float, float]
    burn_evaluations: tuple[tuple[str, PropulsionDynamicsEvaluation], ...]


@dataclass(frozen=True, slots=True)
class PropulsionTimeline:
    """Ordered burns and jettisons on one clock.

    Parameters
    ----------
    burns : tuple of ScheduledBurn
        At least one burn; stored sorted by ignition time.
    jettisons : tuple of MassJettison, optional
        Discrete mass drops; stored sorted by time.
    allow_overlapping_burns : bool, optional
        Permit concurrent burns, whose forces and drains are then summed
        (for example attitude thrusters during a main burn).
    timeline_id : str, optional
        Non-empty timeline identifier.
    schema : str, optional
        Must equal ``PROPULSION_TIMELINE_SCHEMA``.

    Raises
    ------
    DomainError
        If the schema is unsupported, there is no burn, burn and jettison
        identifiers collide, or burns overlap without
        ``allow_overlapping_burns``.

    Notes
    -----
    Without ``allow_overlapping_burns``, burns may touch at an instant but not
    overlap; the comparison uses exact rational sums.  Burns may share tank
    states, so a restartable stage drains one inventory.  Burn events follow
    an event-detector contract [1]_; jettisons are instantaneous staging mass
    drops [2]_.

    References
    ----------
    .. [1] OREKIT-13.1.5-PROPULSION-MODEL. Orekit 13.1.5, "PropulsionModel /
           EventDetector" API.
    .. [2] SUTTON-BIBLARZ-2017. G. P. Sutton and O. Biblarz, "Rocket Propulsion
           Elements", 9th ed., Wiley, 2017.
    """

    burns: tuple[ScheduledBurn, ...]
    jettisons: tuple[MassJettison, ...] = ()
    allow_overlapping_burns: bool = False
    timeline_id: str = "propulsion-timeline"
    schema: str = PROPULSION_TIMELINE_SCHEMA

    def __post_init__(self) -> None:
        if self.schema != PROPULSION_TIMELINE_SCHEMA:
            raise DomainError("Unsupported propulsion timeline schema.")
        text(self.timeline_id, "Timeline identifier")
        if not isinstance(self.allow_overlapping_burns, bool):
            raise DomainError("allow_overlapping_burns must be a boolean.")
        for field_name, kind in (("burns", ScheduledBurn), ("jettisons", MassJettison)):
            object.__setattr__(
                self, field_name, records(getattr(self, field_name), kind, f"Timeline {field_name}")
            )
        if not self.burns:
            raise DomainError("A propulsion timeline requires at least one burn.")
        identifiers = [burn.burn_id for burn in self.burns] + [
            item.event_id for item in self.jettisons
        ]
        if len(identifiers) != len(set(identifiers)):
            raise DomainError("Burn and jettison identifiers must be unique.")
        object.__setattr__(
            self, "burns", tuple(sorted(self.burns, key=lambda burn: burn.ignition_time_s))
        )
        object.__setattr__(
            self, "jettisons", tuple(sorted(self.jettisons, key=lambda item: item.time_s))
        )
        if not self.allow_overlapping_burns:
            for earlier, later in zip(self.burns, self.burns[1:], strict=False):
                if _fraction(later.ignition_time_s) < (_fraction(earlier.ignition_time_s)
                                                     + _fraction(earlier.bridge.artifact.duration_s)):
                    raise DomainError(
                        f"Burns {earlier.burn_id!r} and {later.burn_id!r} overlap; "
                        "set allow_overlapping_burns to sum concurrent burns."
                    )

    @property
    def bound_states(self) -> tuple[str, ...]:
        """Sorted union of the tank states drained by any burn."""

        return tuple(
            sorted(
                {
                    binding.state_name
                    for burn in self.burns
                    for binding in burn.bridge.stream_bindings
                }
            )
        )

    def burn(self, burn_id: str) -> ScheduledBurn:
        """Return the scheduled burn with identifier ``burn_id``.

        Parameters
        ----------
        burn_id : str
            Identifier of a scheduled burn.

        Returns
        -------
        ScheduledBurn
            The matching burn.

        Raises
        ------
        DomainError
            If no burn has that identifier.
        """

        for burn in self.burns:
            if burn.burn_id == burn_id:
                return burn
        raise DomainError(f"Unknown burn {burn_id!r}.")

    def events(self) -> tuple[TimelineEvent, ...]:
        """Return every timeline root, time events first and sorted by time.

        Returns
        -------
        tuple of TimelineEvent
            Burn-time and jettison events sorted by ``(time_s, event_id)``,
            followed by burn state events in burn order.

        Raises
        ------
        DomainError
            If a burn's profile boundaries cannot be distinguished on the
            absolute clock.

        Notes
        -----
        Burn boundary times are the canonical dates
        ``round(ignition_time_s + t_k)``, each rounded once from the exact
        sum [2]_.  Time events use nominal ignition times; if a burn's
        ``ignition_time_bias_s`` driver is non-zero the consumer shifts that
        burn's time events by the same amount.  Event registration follows
        [1]_.

        References
        ----------
        .. [1] OREKIT-13.1.5-PROPULSION-MODEL. Orekit 13.1.5, "PropulsionModel /
               EventDetector" API.
        .. [2] PYTHON-3.12-FRACTION. Python 3.12, "fractions": float inputs retain
               their exact binary value as a rational,
               https://docs.python.org/3.12/library/fractions.html.
        """

        timed: list[TimelineEvent] = []
        state: list[TimelineEvent] = []
        for burn in self.burns:
            canonical, _ = _epoch_table(burn.bridge.artifact.segments, burn.ignition_time_s, 0.0)
            boundaries = (0.0,) + tuple(s.end_time_s for s in burn.bridge.artifact.segments)
            dates = dict(zip(boundaries, canonical, strict=True))
            for surface in burn.bridge.event_surfaces():
                if surface.relative_time_s is not None:
                    timed.append(
                        TimelineEvent(
                            f"{burn.burn_id}/{surface.event_id}",
                            TimelineEventKind.BURN_TIME,
                            burn.burn_id,
                            time_s=dates[surface.relative_time_s],
                            surface=surface,
                            stops_burn=surface.stops_burn,
                        )
                    )
                else:
                    state.append(
                        TimelineEvent(
                            f"{burn.burn_id}/{surface.event_id}",
                            TimelineEventKind.BURN_STATE,
                            burn.burn_id,
                            surface=surface,
                            stops_burn=surface.stops_burn,
                        )
                    )
        timed.extend(
            TimelineEvent(
                f"jettison/{item.event_id}",
                TimelineEventKind.JETTISON,
                item.event_id,
                time_s=item.time_s,
                jettison_mass_kg=item.mass_kg,
            )
            for item in self.jettisons
        )
        timed.sort(key=lambda event: (event.time_s, event.event_id))
        return tuple(timed) + tuple(state)

    def jettisoned_mass_kg(self, time_s: float, *, include_at_time: bool = True) -> float:
        """Return the total mass dropped up to ``time_s``, in kg.

        Parameters
        ----------
        time_s : float
            Timeline time, in s.
        include_at_time : bool, optional
            Include jettisons scheduled exactly at ``time_s``.

        Returns
        -------
        float
            Sum of the dropped masses.

        Raises
        ------
        DomainError
            If ``time_s`` is not a finite real number.

        References
        ----------
        .. [1] SUTTON-BIBLARZ-2017. G. P. Sutton and O. Biblarz, "Rocket Propulsion
               Elements", 9th ed., Wiley, 2017.
        """

        time_s = real(time_s, "Timeline time")
        return sum(
            item.mass_kg
            for item in self.jettisons
            if item.time_s < time_s or (include_at_time and item.time_s == time_s)
        )

    def apply_jettison(self, event_id: str, vehicle_mass_kg: float) -> float:
        """Return the vehicle mass after one jettison event.

        Parameters
        ----------
        event_id : str
            Identifier of a jettison on this timeline.
        vehicle_mass_kg : float
            Vehicle mass just before the drop, in kg.

        Returns
        -------
        float
            ``vehicle_mass_kg - mass_kg`` of the event, in kg.

        Raises
        ------
        DomainError
            If the event is unknown, the mass is not a finite real number, or
            the remaining mass would be non-positive.

        References
        ----------
        .. [1] SUTTON-BIBLARZ-2017. G. P. Sutton and O. Biblarz, "Rocket Propulsion
               Elements", 9th ed., Wiley, 2017.
        """

        vehicle_mass_kg = real(vehicle_mass_kg, "Vehicle mass")
        for item in self.jettisons:
            if item.event_id == event_id:
                remaining = vehicle_mass_kg - item.mass_kg
                if remaining <= 0.0:
                    raise DomainError(
                        f"Jettison {event_id!r} would leave a non-positive vehicle mass."
                    )
                return remaining
        raise DomainError(f"Unknown jettison event {event_id!r}.")

    def evaluate(self, *, time_s: float, vehicle_mass_kg: float,
                 directions: Mapping[str, tuple[float,float,float]],
                 additional_states: Mapping[str,float] | None = None,
                 stopped_burns: frozenset[str] | set[str] = frozenset(),
                 ambient_pressure_pa: float | None = None,
                 operating_conditions: Mapping[str,Mapping[str,float]] | None = None,
                 parameter_overrides: Mapping[str,Mapping[str,float]] | None = None,
                 side: str = 'right', apply_inventory_limits: bool = True) -> TimelineEvaluation:
        """Sum the dynamics of every burn that is firing at ``time_s``.

        Parameters
        ----------
        time_s : float
            Timeline time, in s.
        vehicle_mass_kg : float
            Current vehicle mass, in kg.
        directions : Mapping of str to tuple of float
            Thrust direction per burn identifier; required only for burns
            that are active at ``time_s``.
        additional_states : Mapping of str to float or None, optional
            Tank-state values, in kg.
        stopped_burns : set of str, optional
            Burns the consumer stopped at a state root; they are skipped.
        ambient_pressure_pa : float or None, optional
            Ambient pressure, in Pa, forwarded only to burns whose bridge has
            a nozzle exit area.
        operating_conditions : Mapping of str to Mapping or None, optional
            Operating conditions keyed by burn identifier.
        parameter_overrides : Mapping of str to Mapping or None, optional
            Parameter-driver overrides keyed by burn identifier.
        side : {"right", "left"}, optional
            Limit taken at a profile discontinuity.
        apply_inventory_limits : bool, optional
            Inhibit burns at dry-mass and reserve limits; pass ``False`` when
            the consumer locates those roots itself.

        Returns
        -------
        TimelineEvaluation
            Summed accelerations, mass rates, tank drains and per-burn
            evaluations.

        Raises
        ------
        DomainError
            If an input is malformed, a stopped burn is unknown, an active
            burn has no direction, or more than one burn is active on a
            non-overlapping timeline.

        Notes
        -----
        Each burn is evaluated on the complete-profile clock: a timeline
        reading equal to a canonical boundary date (or a supported float
        alias of it) maps to that profile boundary exactly, and any other
        reading maps to ``time_s - ignition_time_s - bias`` without rounding
        [2]_.  The force and mass-derivative contract follows [1]_.

        References
        ----------
        .. [1] OREKIT-13.1.5-PROPULSION-MODEL. Orekit 13.1.5, "PropulsionModel /
               EventDetector" API.
        .. [2] PYTHON-3.12-FRACTION. Python 3.12, "fractions": float inputs retain
               their exact binary value as a rational,
               https://docs.python.org/3.12/library/fractions.html.
        """
        arguments=locals().copy()
        arguments.pop('self')
        return self._evaluate(**arguments)

    def _evaluate(
        self,
        *,
        time_s: float,
        vehicle_mass_kg: float,
        directions: Mapping[str, tuple[float, float, float]],
        additional_states: Mapping[str, float] | None = None,
        stopped_burns: frozenset[str] | set[str] = frozenset(),
        ambient_pressure_pa: float | None = None,
        operating_conditions: Mapping[str, Mapping[str, float]] | None = None,
        parameter_overrides: Mapping[str, Mapping[str, float]] | None = None,
        side: str = "right",
        apply_inventory_limits: bool = True,
        _clock_horizons: Mapping[str, object] | None = None,
        _canonical_epoch: bool = False,
    ) -> TimelineEvaluation:
        """Implement :meth:`evaluate` with optional clock horizons.

        ``_clock_horizons`` maps burn identifiers to exact cutoff horizons of
        their profile clocks; ``_canonical_epoch`` disables float aliases so
        only canonical boundary dates map exactly to boundaries [2]_.  The
        other arguments are those of :meth:`evaluate` and the contract
        follows [1]_.

        References
        ----------
        .. [1] OREKIT-13.1.5-PROPULSION-MODEL. Orekit 13.1.5, "PropulsionModel /
               EventDetector" API.
        .. [2] PYTHON-3.12-FRACTION. Python 3.12, "fractions": float inputs retain
               their exact binary value as a rational,
               https://docs.python.org/3.12/library/fractions.html.
        """

        time_s = real(time_s, "Timeline time")
        vehicle_mass_kg = real(vehicle_mass_kg, "Vehicle mass")
        if not isinstance(directions, Mapping):
            raise DomainError("Directions must be a mapping keyed by burn identifier.")
        states = real_mapping(additional_states, "Additional states")
        if isinstance(stopped_burns, str) or not hasattr(stopped_burns, "__contains__"):
            raise DomainError("stopped_burns must be a set of burn identifiers.")
        unknown = set(stopped_burns) - {burn.burn_id for burn in self.burns}
        if unknown:
            raise DomainError("Unknown stopped burns: " + ", ".join(sorted(map(str, unknown))))
        for name, value in (
            ("operating_conditions", operating_conditions),
            ("parameter_overrides", parameter_overrides),
        ):
            if value is not None and not isinstance(value, Mapping):
                raise DomainError(f"{name} must be a mapping keyed by burn identifier.")
        conditions_by_burn = dict(operating_conditions or {})
        overrides_by_burn = dict(parameter_overrides or {})

        acceleration = [0.0, 0.0, 0.0]
        acceleration_wrt_mass = [0.0, 0.0, 0.0]
        mass_rate = 0.0
        drains = {name: 0.0 for name in self.bound_states}
        evaluations: list[tuple[str, PropulsionDynamicsEvaluation]] = []
        active: list[str] = []
        for burn in self.burns:
            if burn.burn_id in stopped_burns:
                continue
            burn_overrides = real_mapping(
                overrides_by_burn.get(burn.burn_id), f"Parameter overrides for {burn.burn_id!r}"
            )
            # Only the timing driver decides activity; full validation happens
            # once, inside the bridge evaluation below.
            timing = next(
                parameter
                for parameter in burn.bridge.parameters
                if parameter.name == "ignition_time_bias_s"
            )
            bias = timing.validate(
                burn_overrides.get("ignition_time_bias_s", timing.reference_value)
            )
            effective = _epoch_profile_time(burn.bridge.artifact.segments, time_s,
                burn.ignition_time_s, bias,
                horizon=None if _clock_horizons is None else _clock_horizons[burn.burn_id],
                allow_aliases=not _canonical_epoch)
            if effective < 0 or effective > _fraction(burn.bridge.artifact.duration_s):
                continue
            if burn.burn_id not in directions:
                raise DomainError(f"Burn {burn.burn_id!r} is active but has no direction.")
            value = burn.bridge._evaluate(
                relative_time_s=_finite_rational(_fraction(time_s)-_fraction(burn.ignition_time_s),
                                                 'Timeline relative time'),
                _resolved_profile_time=effective,
                vehicle_mass_kg=vehicle_mass_kg,
                direction=directions[burn.burn_id],
                additional_states=states,
                operating_conditions=conditions_by_burn.get(burn.burn_id),
                parameter_overrides=burn_overrides,
                side=side,
                apply_inventory_limits=apply_inventory_limits,
                ambient_pressure_pa=(
                    ambient_pressure_pa if burn.bridge.nozzle_exit_area_m2 is not None else None
                ),
            )
            evaluations.append((burn.burn_id, value))
            if not value.active:
                continue
            active.append(burn.burn_id)
            for index in range(3):
                acceleration[index] += value.acceleration_m_s2[index]
                acceleration_wrt_mass[index] += value.partials.acceleration_wrt_mass[index]
            mass_rate += value.mass_derivative_kg_s
            for name, rate in value.additional_state_derivatives:
                drains[name] = drains.get(name, 0.0) + rate
        if not self.allow_overlapping_burns and len(active) > 1:
            raise DomainError("More than one burn is active on a non-overlapping timeline.")
        return TimelineEvaluation(
            time_s=time_s,
            active_burns=tuple(active),
            acceleration_m_s2=(acceleration[0], acceleration[1], acceleration[2]),
            mass_derivative_kg_s=mass_rate,
            additional_state_derivatives=tuple(sorted(drains.items())),
            acceleration_wrt_mass=(
                acceleration_wrt_mass[0],
                acceleration_wrt_mass[1],
                acceleration_wrt_mass[2],
            ),
            burn_evaluations=tuple(evaluations),
        )

    def manifest(self, *, include_hash: bool = True) -> dict[str, object]:
        """Return a deterministic manifest binding every burn's model hash.

        Parameters
        ----------
        include_hash : bool, optional
            Add ``timeline_hash`` computed over the other members.

        Returns
        -------
        dict
            Schema, time-numerics policy, burns with bridge model hashes,
            and jettisons.
        """

        payload: dict[str, object] = {
            "schema": self.schema,
            "time_numerics": _TIME_POLICY,
            "timeline_id": self.timeline_id,
            "allow_overlapping_burns": self.allow_overlapping_burns,
            "burns": [
                {
                    "burn_id": burn.burn_id,
                    "ignition_time_s": burn.ignition_time_s,
                    "bridge_model_hash": burn.bridge.model_hash,
                }
                for burn in self.burns
            ],
            "jettisons": [
                {
                    "event_id": item.event_id,
                    "time_s": item.time_s,
                    "mass_kg": item.mass_kg,
                    "label": item.label,
                }
                for item in self.jettisons
            ],
        }
        if include_hash:
            payload["timeline_hash"] = canonical_sha256(payload)
        return payload

    @property
    def timeline_hash(self) -> str:
        """SHA-256 content identity of the schedule, jettisons and every bridge."""

        return canonical_sha256(self.manifest(include_hash=False))

    def to_dict(self) -> dict[str, object]:
        """Return the complete document: the manifest plus every bridge document.

        Returns
        -------
        dict
            Manifest with a ``bridges`` member keyed by burn identifier.
        """

        document = self.manifest()
        document["bridges"] = {burn.burn_id: burn.bridge.to_dict() for burn in self.burns}
        return document

    def to_json(self, *, indent: int | None = 2) -> str:
        """Serialize the complete timeline document as JSON text.

        Parameters
        ----------
        indent : int or None, optional
            Indentation passed to :func:`json.dumps`; ``None`` is compact.

        Returns
        -------
        str
            Key-sorted JSON without ``NaN`` or ``Infinity`` tokens [1]_.

        References
        ----------
        .. [1] IETF RFC 8259, "The JavaScript Object Notation (JSON) Data
               Interchange Format", https://www.rfc-editor.org/rfc/rfc8259.
        """

        return json.dumps(
            self.to_dict(), allow_nan=False, ensure_ascii=False, indent=indent, sort_keys=True
        )

    @classmethod
    def from_dict(cls, payload: Mapping[str, object]) -> PropulsionTimeline:
        """Rebuild a timeline, verifying every bridge and the timeline hash.

        Parameters
        ----------
        payload : Mapping
            Document produced by :meth:`to_dict`.

        Returns
        -------
        PropulsionTimeline
            The reconstructed timeline.

        Raises
        ------
        InputError
            If the document is malformed, declares a different time-numerics
            policy, or a timeline or bridge hash does not match its contents.

        Notes
        -----
        Hash-verified reopening supports the reproducibility and traceability
        of model inputs required by [1]_.

        References
        ----------
        .. [1] NASA-STD-7009B. NASA, "Standard for Models and Simulations",
               https://standards.nasa.gov/standard/NASA/NASA-STD-7009.
        """

        _timeline_document(payload)
        _require_numerical_policies(payload, {'time_numerics':_TIME_POLICY})
        try:
            bridges = payload["bridges"]
            timeline = cls(
                burns=tuple(
                    ScheduledBurn(
                        item["burn_id"],
                        PropagatorPropulsionBridge.from_dict(bridges[item["burn_id"]]),
                        item["ignition_time_s"],
                    )
                    for item in payload["burns"]
                ),
                jettisons=tuple(
                    MassJettison(item["event_id"], item["time_s"], item["mass_kg"], item["label"])
                    for item in payload["jettisons"]
                ),
                allow_overlapping_burns=payload["allow_overlapping_burns"],
                timeline_id=payload["timeline_id"],
                schema=payload["schema"],
            )
        except InputError:
            raise
        except (KeyError, TypeError, ValueError, AttributeError) as error:
            raise InputError("Malformed propulsion timeline document.") from error
        if payload.get("timeline_hash") != timeline.timeline_hash:
            raise InputError("Propulsion timeline hash does not match its contents.")
        for item in payload['burns']:
            if item['bridge_model_hash'] != timeline.burn(item['burn_id']).bridge.model_hash:
                raise InputError("Scheduled burn bridge hash does not match its bridge.",
                                 field='bridge_model_hash')
        return timeline

    @classmethod
    def from_json(cls, text_value: str) -> PropulsionTimeline:
        """Parse JSON text and reopen it with :meth:`from_dict`.

        Parameters
        ----------
        text_value : str
            RFC 8259 JSON document [1]_.

        Returns
        -------
        PropulsionTimeline
            The reconstructed, hash-verified timeline.

        Raises
        ------
        InputError
            If the text is not a valid document or fails a hash check.

        References
        ----------
        .. [1] IETF RFC 8259, "The JavaScript Object Notation (JSON) Data
               Interchange Format", https://www.rfc-editor.org/rfc/rfc8259.
        """

        return read_json_document(text_value, "Propulsion timeline document", cls.from_dict)
