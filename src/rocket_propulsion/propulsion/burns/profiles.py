"""Piecewise-linear throttle geometry with exact supplied-binary integrals."""

from __future__ import annotations

from dataclasses import dataclass
from decimal import Context, Decimal, localcontext
from fractions import Fraction
from math import isqrt

from rocket_propulsion.core.errors import DomainError
from rocket_propulsion.core.validation import real, text

from .models import BurnDefinition, EngineOperatingPoint, ThrottleSchedule
from .numeric_outputs import _finite_rational, _fraction


@dataclass(frozen=True, slots=True)
class ScheduleInterval:
    """One schedule interval with public binary64 endpoint coordinates."""

    start_time_s: float
    end_time_s: float
    start_throttle: float
    end_throttle: float
    phase: str

    def __post_init__(self):
        for name in ("start_time_s", "end_time_s", "start_throttle", "end_throttle"):
            object.__setattr__(self, name, real(getattr(self, name), f"Schedule interval {name}"))
        if self.start_time_s < 0 or self.end_time_s < self.start_time_s:
            raise DomainError("Schedule interval times are invalid.")
        if any(not 0 <= value <= 1 for value in (self.start_throttle, self.end_throttle)):
            raise DomainError("Schedule interval throttle is outside [0, 1].")
        text(self.phase, "Schedule interval phase")

    @property
    def duration_s(self) -> float:
        return _finite_rational(
            _fraction(self.end_time_s) - _fraction(self.start_time_s), "Interval duration"
        )

    @property
    def exposure_s(self) -> float:
        return _finite_rational(
            (_fraction(self.start_throttle) + _fraction(self.end_throttle))
            * (_fraction(self.end_time_s) - _fraction(self.start_time_s))
            / 2,
            "Interval exposure",
        )

    @property
    def throttle_first_moment_s2(self) -> float:
        duration = _fraction(self.end_time_s) - _fraction(self.start_time_s)
        a, b = _fraction(self.start_throttle), _fraction(self.end_throttle)
        value = _fraction(self.start_time_s) * (a + b) * duration / 2 + duration**2 * (
            a / 2 + (b - a) / 3
        )
        return _finite_rational(value, "Interval first moment")


@dataclass(frozen=True, slots=True)
class _ScheduleSlice:
    start: Fraction
    duration: Fraction
    a: Fraction
    b: Fraction
    phase: str

    @property
    def exposure(self):
        return (self.a + self.b) * self.duration / 2

    @property
    def moment(self):
        return self.start * self.exposure + self.duration**2 * (self.a / 2 + (self.b - self.a) / 3)


def _schedule_totals(schedule):
    """Do not require whole-schedule scalar outputs for an early finite query."""
    if not isinstance(schedule, ThrottleSchedule):
        raise DomainError("A throttle schedule is required.")
    duration = sum((_fraction(s.duration_s) for s in schedule.segments), Fraction(0))
    exposure = sum(
        (
            (_fraction(s.start_throttle) + _fraction(s.end_throttle)) * _fraction(s.duration_s) / 2
            for s in schedule.segments
        ),
        Fraction(0),
    )
    return duration, exposure


def _schedule_slices(schedule, stop_time=None):
    total, _ = _schedule_totals(schedule)
    stop = total if stop_time is None else _fraction(stop_time, "Schedule stop time")
    if stop < 0:
        raise DomainError("Schedule stop time must be non-negative.")
    stop = min(stop, total)
    cursor = Fraction(0)
    slices = []
    for segment in schedule.segments:
        if cursor >= stop:
            break
        duration = _fraction(segment.duration_s)
        used = min(duration, stop - cursor)
        a = _fraction(segment.start_throttle)
        b = a + (_fraction(segment.end_throttle) - a) * used / duration
        slices.append(_ScheduleSlice(cursor, used, a, b, segment.phase))
        cursor += duration
    return tuple(slices)


def _schedule_integrals(schedule, time):
    slices = _schedule_slices(schedule, time)
    return (
        sum((s.exposure for s in slices), Fraction(0)),
        sum((s.moment for s in slices), Fraction(0)),
    )


def _sqrt_fraction(value):
    """Preserve rational roots exactly; otherwise use an isolated Decimal1100 root."""
    numerator, denominator = isqrt(value.numerator), isqrt(value.denominator)
    if numerator * numerator == value.numerator and denominator * denominator == value.denominator:
        return Fraction(numerator, denominator)
    with localcontext(Context(prec=1100)):
        return Fraction((Decimal(value.numerator) / Decimal(value.denominator)).sqrt())


def _inverse_exposure(schedule, exposure):
    """Earliest inverse, without an absolute tolerance or a dimensional slope cutoff."""
    requested = _fraction(exposure, "Requested throttle exposure")
    if requested < 0:
        raise DomainError("Requested throttle exposure must be non-negative.")
    _, total = _schedule_totals(schedule)
    if requested > total:
        return None
    if requested == 0:
        return Fraction(0)
    accumulated = cursor = Fraction(0)
    for segment in schedule.segments:
        duration, a, b = map(
            _fraction, (segment.duration_s, segment.start_throttle, segment.end_throttle)
        )
        segment_exposure = (a + b) * duration / 2
        if accumulated + segment_exposure < requested:
            accumulated += segment_exposure
            cursor += duration
            continue
        remaining = requested - accumulated
        if remaining == 0:
            return cursor
        if remaining == segment_exposure:
            return cursor + duration
        if a == b:
            local_time = remaining / a
        else:
            normalized = remaining / duration
            discriminant = a * a + 2 * (b - a) * normalized
            root = _sqrt_fraction(discriminant)
            local_time = duration * 2 * normalized / (a + root)
        return cursor + min(max(local_time, Fraction(0)), duration)
    raise DomainError("Reachable exposure did not resolve to a schedule interval.")


def resolved_schedule_intervals(
    schedule: ThrottleSchedule, *, stop_time_s: float | None = None
) -> tuple[ScheduleInterval, ...]:
    """Export slices once; fail explicitly if two distinct endpoints collapse."""
    result = []
    for s in _schedule_slices(schedule, stop_time_s):
        start = _finite_rational(s.start, "Interval start time")
        end = _finite_rational(s.start + s.duration, "Interval end time")
        if end <= start:
            raise DomainError(
                "Interval endpoints cannot be distinguished in binary64.",
                field="Interval time resolution",
            )
        result.append(
            ScheduleInterval(
                start,
                end,
                _finite_rational(s.a, "Interval start throttle"),
                _finite_rational(s.b, "Interval end throttle"),
                s.phase,
            )
        )
    return tuple(result)


def integrated_throttle(schedule: ThrottleSchedule, time_s: float) -> float:
    """Export the integral of the supplied-binary piecewise-linear command."""
    return _finite_rational(
        _schedule_integrals(schedule, _fraction(time_s, "Schedule query time"))[0],
        "Integrated throttle exposure",
    )


def throttle_first_moment(schedule: ThrottleSchedule, time_s: float) -> float:
    """Export the first time moment, with a typed error for genuine range loss."""
    return _finite_rational(
        _schedule_integrals(schedule, _fraction(time_s, "Schedule query time"))[1],
        "Throttle first moment",
    )


def time_for_integrated_throttle(schedule: ThrottleSchedule, exposure_s: float) -> float | None:
    """Return the earliest reachable command time, or None beyond total exposure."""
    value = _inverse_exposure(schedule, exposure_s)
    return None if value is None else _finite_rational(value, "Schedule inverse time")


@dataclass(frozen=True, slots=True)
class ProfiledPerformanceProvider:
    """L1 provider applying one throttle schedule to a resolved rated point."""

    rated_point: EngineOperatingPoint
    schedule: ThrottleSchedule
    provider_id: str = "profiled-performance"
    provider_version: str = "1.0"

    def operating_point(self) -> EngineOperatingPoint:
        return self.rated_point

    def discontinuities_s(self, definition: BurnDefinition) -> tuple[float, ...]:
        del definition
        cursor = Fraction(0)
        boundaries = []
        for segment in self.schedule.segments[:-1]:
            cursor += _fraction(segment.duration_s)
            value = _finite_rational(cursor, "Schedule discontinuity time")
            if boundaries and value <= boundaries[-1]:
                raise DomainError(
                    "Schedule discontinuities cannot be distinguished in binary64.",
                    field="Discontinuity time resolution",
                )
            boundaries.append(value)
        return tuple(boundaries)
