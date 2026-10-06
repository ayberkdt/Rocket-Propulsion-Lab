"""Restart, cooldown, and minimum-pulse validation for engine commands."""

from __future__ import annotations

from dataclasses import dataclass
from itertools import pairwise
from math import isfinite

from rocket_propulsion.core.errors import DomainError

from .models import ThrottleSchedule, ThrottleSegment


def _finite_nonnegative(value: float, label: str) -> None:
    if not isfinite(value) or value < 0.0:
        raise DomainError(f"{label} must be finite and non-negative.")


@dataclass(frozen=True, slots=True)
class EngineDutyCycleLimits:
    """Operational timing limits for one independently commanded engine.

    References
    ----------
    NASA-CR-122347: https://ntrs.nasa.gov/citations/19720011120
    """

    engine_id: str
    minimum_on_time_s: float = 0.0
    minimum_off_time_s: float = 0.0
    maximum_starts: int | None = None
    ignition_delay_s: float = 0.0
    response_time_constant_s: float = 0.0

    def __post_init__(self) -> None:
        if not self.engine_id.strip():
            raise DomainError("Engine identifier must not be empty.")
        for value, label in (
            (self.minimum_on_time_s, "Minimum on-time"),
            (self.minimum_off_time_s, "Minimum off-time"),
            (self.ignition_delay_s, "Ignition delay"),
            (self.response_time_constant_s, "Response time constant"),
        ):
            _finite_nonnegative(value, label)
        if self.maximum_starts is not None and self.maximum_starts < 1:
            raise DomainError("Maximum starts must be at least one when declared.")


@dataclass(frozen=True, slots=True)
class EnginePulse:
    """One commanded on-interval for a named engine.

    References
    ----------
    NASA-CR-122347: https://ntrs.nasa.gov/citations/19720011120
    """

    engine_id: str
    start_time_s: float
    duration_s: float
    throttle: float = 1.0

    def __post_init__(self) -> None:
        if not self.engine_id.strip():
            raise DomainError("Engine pulse requires an engine identifier.")
        _finite_nonnegative(self.start_time_s, "Pulse start time")
        if not isfinite(self.duration_s) or self.duration_s <= 0.0:
            raise DomainError("Pulse duration must be finite and greater than zero.")
        if not isfinite(self.throttle) or not 0.0 < self.throttle <= 1.0:
            raise DomainError("Pulse throttle must be finite and in (0, 1].")

    @property
    def end_time_s(self) -> float:
        return self.start_time_s + self.duration_s


@dataclass(frozen=True, slots=True)
class EngineSequenceValidation:
    """Normalized, successfully validated multi-engine command sequence.

    References
    ----------
    NASA-CR-122347: https://ntrs.nasa.gov/citations/19720011120
    """

    pulses: tuple[EnginePulse, ...]
    starts_by_engine: tuple[tuple[str, int], ...]
    duration_s: float

    def schedule_for_engine(
        self, engine_id: str, *, total_duration_s: float | None = None
    ) -> ThrottleSchedule:
        """Convert one engine's pulse train into a discontinuity-safe schedule."""

        selected = tuple(pulse for pulse in self.pulses if pulse.engine_id == engine_id)
        if not selected:
            raise DomainError(f"No validated pulses exist for engine {engine_id!r}.")
        stop = self.duration_s if total_duration_s is None else total_duration_s
        if not isfinite(stop) or stop < selected[-1].end_time_s:
            raise DomainError("Schedule duration must include every selected pulse.")
        segments: list[ThrottleSegment] = []
        cursor = 0.0
        for index, pulse in enumerate(selected, start=1):
            if pulse.start_time_s > cursor:
                segments.append(
                    ThrottleSegment(
                        pulse.start_time_s - cursor, 0.0, 0.0, f"off-before-start-{index}"
                    )
                )
            segments.append(
                ThrottleSegment(
                    pulse.duration_s,
                    pulse.throttle,
                    pulse.throttle,
                    f"firing-{index}",
                )
            )
            cursor = pulse.end_time_s
        if stop > cursor:
            segments.append(ThrottleSegment(stop - cursor, 0.0, 0.0, "off-after-final"))
        return ThrottleSchedule(tuple(segments), name=f"{engine_id} validated pulse sequence")


def validate_engine_sequence(
    pulses: tuple[EnginePulse, ...],
    limits: tuple[EngineDutyCycleLimits, ...],
) -> EngineSequenceValidation:
    """Validate timing rules and return a deterministic normalized sequence.

    Limits are supplied as hardware-specific constraints; this function only
    enforces them and does not claim generic minimum pulse or cooldown values.

    References
    ----------
    NASA-CR-122347: https://ntrs.nasa.gov/citations/19720011120
    """

    if not pulses:
        raise DomainError("An engine command sequence requires at least one pulse.")
    if not limits:
        raise DomainError("Engine duty-cycle limits are required.")
    limits_by_engine = {item.engine_id: item for item in limits}
    if len(limits_by_engine) != len(limits):
        raise DomainError("Engine duty-cycle limit identifiers must be unique.")
    unknown = sorted({pulse.engine_id for pulse in pulses} - limits_by_engine.keys())
    if unknown:
        raise DomainError(f"No duty-cycle limits declared for engines: {', '.join(unknown)}.")

    normalized = tuple(sorted(pulses, key=lambda item: (item.engine_id, item.start_time_s)))
    starts: list[tuple[str, int]] = []
    for engine_id in sorted({pulse.engine_id for pulse in normalized}):
        engine_pulses = tuple(pulse for pulse in normalized if pulse.engine_id == engine_id)
        engine_limits = limits_by_engine[engine_id]
        if (
            engine_limits.maximum_starts is not None
            and len(engine_pulses) > engine_limits.maximum_starts
        ):
            raise DomainError(
                f"Engine {engine_id!r} exceeds its maximum start count "
                f"({len(engine_pulses)} > {engine_limits.maximum_starts})."
            )
        for pulse in engine_pulses:
            if pulse.duration_s + 1e-12 < engine_limits.minimum_on_time_s:
                raise DomainError(
                    f"Engine {engine_id!r} pulse at {pulse.start_time_s:g} s is shorter "
                    "than its minimum on-time."
                )
        for previous, current in pairwise(engine_pulses):
            gap = current.start_time_s - previous.end_time_s
            if gap < -1e-12:
                raise DomainError(f"Engine {engine_id!r} pulses overlap.")
            if gap + 1e-12 < engine_limits.minimum_off_time_s:
                raise DomainError(
                    f"Engine {engine_id!r} cooldown is {gap:g} s; at least "
                    f"{engine_limits.minimum_off_time_s:g} s is required."
                )
        starts.append((engine_id, len(engine_pulses)))

    return EngineSequenceValidation(
        pulses=normalized,
        starts_by_engine=tuple(starts),
        duration_s=max(pulse.end_time_s for pulse in normalized),
    )
