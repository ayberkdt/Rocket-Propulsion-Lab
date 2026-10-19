"""Multi-burn propulsion timelines with discrete mass jettison.

A single :class:`~.propagator.PropagatorPropulsionBridge` represents one burn.
A mission arc usually needs several: an ascent stage, a coast, a separation,
an upper-stage burn, a later restart.  :class:`PropulsionTimeline` places
burns on one timeline clock, checks that they do not overlap unless that is
declared, adds stage-separation / fairing jettison events as discrete mass
drops, and evaluates all active burns in one call.

The timeline still computes no orbit.  Directions, ambient pressure and
operating conditions are supplied per burn by the consumer, which also owns
event location and the mode state that records which burns were stopped at an
inventory root (passed back as ``stopped_burns``).
"""

from __future__ import annotations

import json
from collections.abc import Mapping
from dataclasses import dataclass
from enum import Enum

from rocket_propulsion.core.errors import DomainError, InputError
from rocket_propulsion.core.validation import read_json_document, real, real_mapping, records, text

from .propagator import (
    BurnEventSurface,
    PropagatorPropulsionBridge,
    PropulsionDynamicsEvaluation,
)
from .serialization import canonical_sha256, _require_numerical_policies
from .numeric_outputs import _fraction, _finite_rational
from .time_coordinates import _TIME_POLICY, _epoch_profile_time, _epoch_table

PROPULSION_TIMELINE_SCHEMA = "rocket_propulsion_propulsion_timeline_v1"




@dataclass(frozen=True, slots=True)
class ScheduledBurn:
    """One bridge placed on the timeline clock.

    ``ignition_time_s`` is the timeline time of profile time zero.  The
    bridge's own ``ignition_time_bias_s`` driver still shifts the burn
    relative to that nominal time.

    References
    ----------
    Orekit Maneuver 13.1.5 API:
    https://www.orekit.org/static/apidocs/org/orekit/forces/maneuvers/Maneuver.html
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
        """Return the nominal timeline time of the profile end."""

        return _finite_rational(_fraction(self.ignition_time_s)+_fraction(self.bridge.artifact.duration_s),
                                'Scheduled cutoff time')


@dataclass(frozen=True, slots=True)
class MassJettison:
    """A discrete mass drop (stage, fairing, adapter) at a timeline time.

    References
    ----------
    NASA-SP-125: https://ntrs.nasa.gov/citations/19710019929
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

    References
    ----------
    Orekit EventDetector 13.1.5 API:
    https://www.orekit.org/static/apidocs/org/orekit/propagation/events/EventDetector.html
    """

    BURN_TIME = "burn_time"
    BURN_STATE = "burn_state"
    JETTISON = "jettison"


@dataclass(frozen=True, slots=True)
class TimelineEvent:
    """One timeline root with its owner and meaning.

    ``time_s`` is set for time events (burn start/boundary/stop at nominal
    ignition, jettison).  State events carry the burn's own
    :class:`BurnEventSurface`; when one fires, the consumer adds ``burn_id``
    to ``stopped_burns``.

    References
    ----------
    Orekit EventDetector 13.1.5 API:
    https://www.orekit.org/static/apidocs/org/orekit/propagation/events/EventDetector.html
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

    References
    ----------
    Orekit PropulsionModel 13.1.5 API:
    https://www.orekit.org/static/apidocs/org/orekit/forces/maneuvers/propulsion/PropulsionModel.html
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

    Burns may not overlap (touching at an instant is allowed) unless
    ``allow_overlapping_burns`` is set, in which case concurrent forces and
    drains are summed (for example attitude thrusters during a main burn).
    Burns may share tank states, so a restartable stage drains one inventory.

    References
    ----------
    Orekit Maneuver 13.1.5 API:
    https://www.orekit.org/static/apidocs/org/orekit/forces/maneuvers/Maneuver.html
    NASA-SP-125: https://ntrs.nasa.gov/citations/19710019929
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
        """Return the union of tank states drained by any burn."""

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
        """Return one scheduled burn by identifier."""

        for burn in self.burns:
            if burn.burn_id == burn_id:
                return burn
        raise DomainError(f"Unknown burn {burn_id!r}.")

    def events(self) -> tuple[TimelineEvent, ...]:
        """Return every timeline root, time events sorted by time.

        Time events use nominal ignition times; if a burn's
        ``ignition_time_bias_s`` driver is non-zero the consumer shifts that
        burn's time events by the same amount.

        References
        ----------
        Orekit EventDetectorsProvider 13.1.5 API:
        https://www.orekit.org/static/apidocs/org/orekit/propagation/events/EventDetectorsProvider.html
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
        """Return the total mass dropped up to ``time_s``."""

        time_s = real(time_s, "Timeline time")
        return sum(
            item.mass_kg
            for item in self.jettisons
            if item.time_s < time_s or (include_at_time and item.time_s == time_s)
        )

    def apply_jettison(self, event_id: str, vehicle_mass_kg: float) -> float:
        """Return the vehicle mass after one jettison event, failing closed.

        References
        ----------
        NASA-SP-125: https://ntrs.nasa.gov/citations/19710019929
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

    def evaluate(
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
    ) -> TimelineEvaluation:
        """Sum the dynamics of every burn that is firing at ``time_s``.

        ``directions`` needs an entry only for burns that are active.
        ``operating_conditions`` and ``parameter_overrides`` are keyed by
        burn identifier.  ``ambient_pressure_pa`` is forwarded only to burns
        whose bridge has a nozzle exit area.

        References
        ----------
        Orekit PropulsionModel 13.1.5 API:
        https://www.orekit.org/static/apidocs/org/orekit/forces/maneuvers/propulsion/PropulsionModel.html
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
                                            burn.ignition_time_s, bias)
            if effective < 0 or effective > _fraction(burn.bridge.artifact.duration_s):
                continue
            if burn.burn_id not in directions:
                raise DomainError(f"Burn {burn.burn_id!r} is active but has no direction.")
            value = burn.bridge.evaluate_at_epoch(
                time_s=time_s,
                ignition_time_s=burn.ignition_time_s,
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
        """Return a deterministic manifest binding every burn's model hash."""

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
        """Content identity of the schedule, jettisons and every bridge."""

        return canonical_sha256(self.manifest(include_hash=False))

    def to_dict(self) -> dict[str, object]:
        """Return the complete document including every bridge document."""

        document = self.manifest()
        document["bridges"] = {burn.burn_id: burn.bridge.to_dict() for burn in self.burns}
        return document

    def to_json(self, *, indent: int | None = 2) -> str:
        """Serialize the complete timeline document."""

        return json.dumps(
            self.to_dict(), allow_nan=False, ensure_ascii=False, indent=indent, sort_keys=True
        )

    @classmethod
    def from_dict(cls, payload: Mapping[str, object]) -> PropulsionTimeline:
        """Rebuild a timeline, verifying every bridge and the timeline hash.

        References
        ----------
        NASA-STD-7009B: https://standards.nasa.gov/standard/NASA/NASA-STD-7009
        """

        if not isinstance(payload, Mapping):
            raise InputError("Propulsion timeline document root must be an object.")
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
        return timeline

    @classmethod
    def from_json(cls, text_value: str) -> PropulsionTimeline:
        """Deserialize and integrity-check a timeline document."""

        return read_json_document(text_value, "Propulsion timeline document", cls.from_dict)
