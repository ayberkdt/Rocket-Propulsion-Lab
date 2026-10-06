"""Exact piecewise-linear throttle schedules and the L1 provider boundary."""

from __future__ import annotations

from dataclasses import dataclass
from math import isfinite, sqrt

from rocket_propulsion.core.errors import DomainError

from .models import BurnDefinition, EngineOperatingPoint, ThrottleSchedule


@dataclass(frozen=True, slots=True)
class ScheduleInterval:
    """One schedule segment resolved onto absolute command time."""

    start_time_s: float
    end_time_s: float
    start_throttle: float
    end_throttle: float
    phase: str

    @property
    def duration_s(self) -> float:
        return self.end_time_s - self.start_time_s

    @property
    def exposure_s(self) -> float:
        return 0.5 * (self.start_throttle + self.end_throttle) * self.duration_s

    @property
    def throttle_first_moment_s2(self) -> float:
        """Return ``integral(t * throttle(t) dt)`` for this interval."""

        duration = self.duration_s
        local_moment = duration * duration * (
            self.start_throttle / 2.0
            + (self.end_throttle - self.start_throttle) / 3.0
        )
        return self.start_time_s * self.exposure_s + local_moment


def resolved_schedule_intervals(
    schedule: ThrottleSchedule, *, stop_time_s: float | None = None
) -> tuple[ScheduleInterval, ...]:
    """Resolve ordered segments, optionally clipping exactly at ``stop_time_s``."""

    stop = schedule.total_duration_s if stop_time_s is None else float(stop_time_s)
    if not isfinite(stop) or stop < 0.0:
        raise DomainError("Schedule stop time must be finite and non-negative.")
    stop = min(stop, schedule.total_duration_s)
    intervals: list[ScheduleInterval] = []
    cursor = 0.0
    for segment in schedule.segments:
        if cursor >= stop:
            break
        used_duration = min(segment.duration_s, stop - cursor)
        fraction = used_duration / segment.duration_s
        end_throttle = segment.start_throttle + fraction * (
            segment.end_throttle - segment.start_throttle
        )
        intervals.append(
            ScheduleInterval(
                start_time_s=cursor,
                end_time_s=cursor + used_duration,
                start_throttle=segment.start_throttle,
                end_throttle=end_throttle,
                phase=segment.phase,
            )
        )
        cursor += segment.duration_s
    return tuple(intervals)


def integrated_throttle(schedule: ThrottleSchedule, time_s: float) -> float:
    """Return exact full-throttle-equivalent exposure through ``time_s``."""

    return sum(
        interval.exposure_s
        for interval in resolved_schedule_intervals(schedule, stop_time_s=time_s)
    )


def throttle_first_moment(schedule: ThrottleSchedule, time_s: float) -> float:
    """Return exact first time moment of the throttle command."""

    return sum(
        interval.throttle_first_moment_s2
        for interval in resolved_schedule_intervals(schedule, stop_time_s=time_s)
    )


def time_for_integrated_throttle(
    schedule: ThrottleSchedule, exposure_s: float
) -> float | None:
    """Invert ``integral(throttle dt)`` at the earliest reachable command time."""

    if not isfinite(exposure_s) or exposure_s < 0.0:
        raise DomainError("Requested throttle exposure must be finite and non-negative.")
    if exposure_s == 0.0:
        return 0.0
    tolerance = max(1e-13, schedule.total_exposure_s * 1e-12)
    if exposure_s > schedule.total_exposure_s + tolerance:
        return None

    accumulated = 0.0
    cursor = 0.0
    for segment in schedule.segments:
        segment_exposure = segment.exposure_s
        if accumulated + segment_exposure + tolerance < exposure_s:
            accumulated += segment_exposure
            cursor += segment.duration_s
            continue
        remaining = max(0.0, exposure_s - accumulated)
        if remaining <= tolerance:
            return cursor
        slope = (segment.end_throttle - segment.start_throttle) / segment.duration_s
        if abs(slope) <= 1e-15:
            if segment.start_throttle <= 0.0:
                accumulated += segment_exposure
                cursor += segment.duration_s
                continue
            local_time = remaining / segment.start_throttle
        else:
            discriminant = segment.start_throttle**2 + 2.0 * slope * remaining
            if discriminant < -tolerance:
                raise DomainError("Throttle exposure inversion left the segment domain.")
            root = sqrt(max(0.0, discriminant))
            denominator = segment.start_throttle + root
            if denominator <= 0.0:
                raise DomainError("Throttle exposure inversion is singular.")
            local_time = 2.0 * remaining / denominator
        return cursor + min(max(local_time, 0.0), segment.duration_s)
    return schedule.total_duration_s


@dataclass(frozen=True, slots=True)
class ProfiledPerformanceProvider:
    """L1 provider applying one throttle schedule to a resolved rated point."""

    rated_point: EngineOperatingPoint
    schedule: ThrottleSchedule
    provider_id: str = "profiled-performance"
    provider_version: str = "1.0"

    def operating_point(self) -> EngineOperatingPoint:
        """Return the full-throttle point scaled by the schedule in the solver."""

        return self.rated_point

    def discontinuities_s(self, definition: BurnDefinition) -> tuple[float, ...]:
        """Return internal segment boundaries for exact solver splitting."""

        del definition
        cursor = 0.0
        boundaries: list[float] = []
        for segment in self.schedule.segments[:-1]:
            cursor += segment.duration_s
            boundaries.append(cursor)
        return tuple(boundaries)
