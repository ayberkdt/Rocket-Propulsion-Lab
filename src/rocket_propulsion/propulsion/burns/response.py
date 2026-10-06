"""Analytic command-to-throttle response models for transient burns.

This module owns actuator dynamics only.  It deliberately contains no orbit,
reference-frame, attitude, or force-integration concepts, which keeps the
model reusable by Rocket Propulsion Lab and a future Sidera maneuver adapter.
"""

from __future__ import annotations

from dataclasses import dataclass
from math import ceil, exp, isfinite

from rocket_propulsion.core.errors import DomainError

from .models import ThrottleSchedule
from .profiles import resolved_schedule_intervals


def _unit_interval(value: float, label: str) -> None:
    if not isfinite(value) or not 0.0 <= value <= 1.0:
        raise DomainError(f"{label} must be finite and in [0, 1].")


@dataclass(frozen=True, slots=True)
class ResponseSegment:
    """Exact first-order response to one linear command segment.

    References
    ----------
    NASA-TM-107318: https://ntrs.nasa.gov/citations/19970010379
    """

    start_time_s: float
    end_time_s: float
    command_start: float
    command_end: float
    realized_start: float
    realized_end: float
    time_constant_s: float
    phase: str

    def __post_init__(self) -> None:
        if not isfinite(self.start_time_s) or self.start_time_s < 0.0:
            raise DomainError("Response-segment start time must be finite and non-negative.")
        if not isfinite(self.end_time_s) or self.end_time_s <= self.start_time_s:
            raise DomainError("Response-segment end time must be after its start time.")
        for value, label in (
            (self.command_start, "Command start throttle"),
            (self.command_end, "Command end throttle"),
            (self.realized_start, "Realized start throttle"),
            (self.realized_end, "Realized end throttle"),
        ):
            _unit_interval(value, label)
        if not isfinite(self.time_constant_s) or self.time_constant_s < 0.0:
            raise DomainError("Response time constant must be finite and non-negative.")
        if not self.phase.strip():
            raise DomainError("Response-segment phase must not be empty.")

    @property
    def duration_s(self) -> float:
        """Return segment duration."""

        return self.end_time_s - self.start_time_s

    @property
    def command_slope_s_inv(self) -> float:
        """Return the linear command slope."""

        return (self.command_end - self.command_start) / self.duration_s

    def realized_throttle(self, time_s: float) -> float:
        """Evaluate the exact realized throttle inside this segment."""

        if not isfinite(time_s) or not self.start_time_s <= time_s <= self.end_time_s:
            raise DomainError("Response evaluation time lies outside the segment.")
        local_time = time_s - self.start_time_s
        slope = self.command_slope_s_inv
        if self.time_constant_s == 0.0:
            value = self.command_start + slope * local_time
        else:
            tau = self.time_constant_s
            transient = self.realized_start - self.command_start + slope * tau
            value = (
                self.command_start
                + slope * (local_time - tau)
                + transient * exp(-local_time / tau)
            )
        return min(1.0, max(0.0, value))

    @property
    def exposure_s(self) -> float:
        """Return the exact integral of realized throttle over the segment."""

        duration = self.duration_s
        slope = self.command_slope_s_inv
        if self.time_constant_s == 0.0:
            return 0.5 * (self.command_start + self.command_end) * duration
        tau = self.time_constant_s
        transient = self.realized_start - self.command_start + slope * tau
        return (
            (self.command_start - slope * tau) * duration
            + 0.5 * slope * duration * duration
            + transient * tau * (1.0 - exp(-duration / tau))
        )

    @property
    def throttle_first_moment_s2(self) -> float:
        """Return the exact absolute-time moment ``integral(t*r(t) dt)``."""

        duration = self.duration_s
        slope = self.command_slope_s_inv
        if self.time_constant_s == 0.0:
            local_moment = duration * duration * (
                self.command_start / 2.0
                + (self.command_end - self.command_start) / 3.0
            )
        else:
            tau = self.time_constant_s
            exponential = exp(-duration / tau)
            transient = self.realized_start - self.command_start + slope * tau
            exponential_moment = tau * tau * (1.0 - exponential) - tau * duration * exponential
            local_moment = (
                (self.command_start - slope * tau) * duration * duration / 2.0
                + slope * duration**3 / 3.0
                + transient * exponential_moment
            )
        return self.start_time_s * self.exposure_s + local_moment


@dataclass(frozen=True, slots=True)
class RealizedThrottleHistory:
    """Contiguous exact response history for one engine command.

    References
    ----------
    NASA-TM-107318: https://ntrs.nasa.gov/citations/19970010379
    NASA transient-modeling overview:
    https://ntrs.nasa.gov/citations/20040000363
    """

    segments: tuple[ResponseSegment, ...]
    model_id: str = "first_order_linear_command_v1"

    def __post_init__(self) -> None:
        if not self.segments:
            raise DomainError("A realized throttle history requires at least one segment.")
        if abs(self.segments[0].start_time_s) > 1e-12:
            raise DomainError("A realized throttle history must begin at zero command time.")
        for previous, current in zip(self.segments, self.segments[1:], strict=False):
            if abs(previous.end_time_s - current.start_time_s) > 1e-12:
                raise DomainError("Response segments must be contiguous in command time.")
            if previous.time_constant_s > 0.0 and abs(
                previous.realized_end - current.realized_start
            ) > 1e-10:
                raise DomainError("Finite-lag response must remain continuous at boundaries.")

    @property
    def duration_s(self) -> float:
        return self.segments[-1].end_time_s

    @property
    def exposure_s(self) -> float:
        return sum(segment.exposure_s for segment in self.segments)

    @property
    def throttle_first_moment_s2(self) -> float:
        return sum(segment.throttle_first_moment_s2 for segment in self.segments)

    @property
    def centroid_time_s(self) -> float:
        return (
            self.throttle_first_moment_s2 / self.exposure_s
            if self.exposure_s > 0.0
            else 0.0
        )

    def realized_throttle(self, time_s: float, *, side: str = "right") -> float:
        """Evaluate throttle with an explicit side at command discontinuities."""

        if side not in {"left", "right"}:
            raise DomainError("Response boundary side must be 'left' or 'right'.")
        if not isfinite(time_s) or time_s < 0.0:
            raise DomainError("Response time must be finite and non-negative.")
        tolerance = max(1e-12, self.duration_s * 1e-12)
        if time_s > self.duration_s + tolerance:
            return 0.0
        if time_s >= self.duration_s:
            return self.segments[-1].realized_end if side == "left" else 0.0
        for index, segment in enumerate(self.segments):
            if segment.start_time_s <= time_s < segment.end_time_s:
                if (
                    side == "left"
                    and index > 0
                    and abs(time_s - segment.start_time_s) <= tolerance
                ):
                    return self.segments[index - 1].realized_end
                return segment.realized_throttle(time_s)
        return 0.0

    def sample_times(self, maximum_step_s: float) -> tuple[float, ...]:
        """Return a boundary-preserving sampling grid for exchange artifacts."""

        if not isfinite(maximum_step_s) or maximum_step_s <= 0.0:
            raise DomainError("Maximum response sample step must be positive and finite.")
        times = {0.0, self.duration_s}
        for segment in self.segments:
            divisions = max(1, ceil(segment.duration_s / maximum_step_s))
            for index in range(divisions + 1):
                times.add(segment.start_time_s + segment.duration_s * index / divisions)
        return tuple(sorted(times))


def _response_end(
    *, command_start: float, command_end: float, realized_start: float, duration_s: float,
    time_constant_s: float,
) -> float:
    if time_constant_s == 0.0:
        return command_end
    slope = (command_end - command_start) / duration_s
    transient = realized_start - command_start + slope * time_constant_s
    value = (
        command_start
        + slope * (duration_s - time_constant_s)
        + transient * exp(-duration_s / time_constant_s)
    )
    return min(1.0, max(0.0, value))


def realize_first_order_schedule(
    schedule: ThrottleSchedule,
    *,
    time_constant_s: float,
    ignition_delay_s: float = 0.0,
    settling_duration_s: float = 0.0,
    initial_throttle: float = 0.0,
) -> RealizedThrottleHistory:
    """Resolve a linear command through ``tau*dr/dt + r = command`` exactly.

    The equation is an explicitly declared engineering approximation, not a
    universal engine law. Its time constant must be identified from hardware
    data or a higher-fidelity component model.

    References
    ----------
    NASA-TM-107318: https://ntrs.nasa.gov/citations/19970010379
    NASA transient-modeling overview:
    https://ntrs.nasa.gov/citations/20040000363
    """

    for value, label in (
        (time_constant_s, "Response time constant"),
        (ignition_delay_s, "Ignition delay"),
        (settling_duration_s, "Settling duration"),
    ):
        if not isfinite(value) or value < 0.0:
            raise DomainError(f"{label} must be finite and non-negative.")
    _unit_interval(initial_throttle, "Initial realized throttle")

    specifications: list[tuple[float, float, float, str]] = []
    if ignition_delay_s > 0.0:
        specifications.append((ignition_delay_s, 0.0, 0.0, "ignition-delay"))
    specifications.extend(
        (interval.duration_s, interval.start_throttle, interval.end_throttle, interval.phase)
        for interval in resolved_schedule_intervals(schedule)
    )
    if settling_duration_s > 0.0:
        specifications.append((settling_duration_s, 0.0, 0.0, "shutdown-settling"))

    cursor = 0.0
    realized = initial_throttle
    segments: list[ResponseSegment] = []
    for duration, command_start, command_end, phase in specifications:
        realized_start = command_start if time_constant_s == 0.0 else realized
        realized_end = _response_end(
            command_start=command_start,
            command_end=command_end,
            realized_start=realized_start,
            duration_s=duration,
            time_constant_s=time_constant_s,
        )
        segments.append(
            ResponseSegment(
                start_time_s=cursor,
                end_time_s=cursor + duration,
                command_start=command_start,
                command_end=command_end,
                realized_start=realized_start,
                realized_end=realized_end,
                time_constant_s=time_constant_s,
                phase=phase,
            )
        )
        cursor += duration
        realized = realized_end
    return RealizedThrottleHistory(tuple(segments))
