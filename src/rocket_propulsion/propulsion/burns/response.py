"""Analytic command-to-throttle response models for transient burns.

This module covers actuator dynamics only.  It contains no orbit,
reference-frame, attitude, or force-integration concepts, so the same model
serves Rocket Propulsion Lab analyses and trajectory consumers such as a
Sidera maneuver adapter.

The realized throttle ``r(t)`` follows a first-order lag of a piecewise-linear
command ``c(t)``::

    tau dr/dt + r = c(t)

solved in closed form on each command segment.  This is a declared
engineering approximation of start, throttle and shutdown transients [1]_
[2]_ [3]_, not a universal engine law; the time constant must be identified
from hardware data or a higher-fidelity component model.

References
----------
.. [1] NASA-TM-107318. NASA, "RL10A-3-3A Rocket Engine Modeling Project",
       https://ntrs.nasa.gov/citations/19970010379.
.. [2] NASA-20040000363. NASA, "Transient Mathematical Modeling for Liquid
       Rocket Engine Systems", https://ntrs.nasa.gov/citations/20040000363.
.. [3] NASA CR-184099. NASA, transient engine simulation covering steady
       balance, throttle transients, starts and shutdowns,
       https://ntrs.nasa.gov/citations/19910011919.
"""

from __future__ import annotations

from dataclasses import dataclass
from math import ceil, exp, isfinite

from rocket_propulsion.core.errors import DomainError

from .models import ThrottleSchedule
from .profiles import resolved_schedule_intervals

#: Largest sampling grid a response history may request (library choice):
#: one million samples bound memory to tens of MB and run time to seconds.
MAXIMUM_RESPONSE_SAMPLES = 1_000_000


def _unit_interval(value: float, label: str) -> None:
    """Raise ``DomainError`` unless ``value`` is finite and in ``[0, 1]``."""
    if not isfinite(value) or not 0.0 <= value <= 1.0:
        raise DomainError(f"{label} must be finite and in [0, 1].")


@dataclass(frozen=True, slots=True)
class ResponseSegment:
    """Closed-form first-order response to one linear command segment.

    Parameters
    ----------
    start_time_s, end_time_s : float
        Segment bounds in command time, in s; ``0 <= start < end``.
    command_start, command_end : float
        Commanded throttle at the segment bounds, in ``[0, 1]``.
    realized_start, realized_end : float
        Realized throttle at the segment bounds, in ``[0, 1]``.
    time_constant_s : float
        First-order time constant ``tau``, in s; zero means the realized
        throttle equals the command.
    phase : str
        Non-empty phase label.

    Raises
    ------
    DomainError
        If a field is non-finite, outside its domain, or the phase is blank.

    Notes
    -----
    With local time ``s = t - start_time_s``, command slope ``k`` and
    ``A = realized_start - command_start + k tau``, the solution of
    ``tau dr/dt + r = command_start + k s`` is::

        r(s) = command_start + k (s - tau) + A exp(-s / tau)

    The first-order lag is an engineering approximation of engine transient
    behaviour [1]_; ``tau`` is a caller-identified input.

    References
    ----------
    .. [1] NASA-TM-107318. NASA, "RL10A-3-3A Rocket Engine Modeling Project",
           https://ntrs.nasa.gov/citations/19970010379.
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
        """Segment duration ``end_time_s - start_time_s``, in s."""

        return self.end_time_s - self.start_time_s

    @property
    def command_slope_s_inv(self) -> float:
        """Linear command slope ``(command_end - command_start) / duration``, in 1/s."""

        return (self.command_end - self.command_start) / self.duration_s

    def realized_throttle(self, time_s: float) -> float:
        """Evaluate the closed-form realized throttle inside this segment.

        Parameters
        ----------
        time_s : float
            Command time within ``[start_time_s, end_time_s]``, in s.

        Returns
        -------
        float
            Realized throttle, clipped to ``[0, 1]``.

        Raises
        ------
        DomainError
            If ``time_s`` is non-finite or outside the segment.

        Notes
        -----
        Evaluates the first-order lag solution given in the class Notes [1]_.

        References
        ----------
        .. [1] NASA-TM-107318. NASA, "RL10A-3-3A Rocket Engine Modeling
               Project", https://ntrs.nasa.gov/citations/19970010379.
        """

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
        """Throttle exposure ``integral(r(t) dt)`` over the segment, in s.

        Notes
        -----
        Closed-form integral of the lag solution::

            (c0 - k tau) D + k D**2 / 2 + A tau (1 - exp(-D / tau))

        with duration ``D``; for ``tau = 0`` it is the trapezoid of the
        linear command.  The integral is of the unclipped solution, which
        stays in ``[0, 1]`` for commands and initial values in ``[0, 1]``.
        Multiplied by thrust it gives the impulse contribution [1]_.

        References
        ----------
        .. [1] NASA-ISP. NASA Glenn Research Center, "Specific Impulse":
               total impulse I = integral F dt,
               https://www1.grc.nasa.gov/beginners-guide-to-aeronautics/specific-impulse/.
        """

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
        """Absolute-time first moment ``integral(t r(t) dt)`` over the segment, in s**2.

        Notes
        -----
        Computed as ``start_time_s * exposure_s`` plus the closed-form local
        moment of the lag solution; the exponential term contributes
        ``A (tau**2 (1 - exp(-D / tau)) - tau D exp(-D / tau))``.
        """

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
    """Contiguous closed-form response history for one engine command.

    Parameters
    ----------
    segments : tuple of ResponseSegment
        At least one segment; the first starts at command time zero and
        consecutive segments share their boundary times within ``1e-12`` s.
    model_id : str, optional
        Identifier of the response model.

    Raises
    ------
    DomainError
        If the history is empty, does not start at zero, has gaps, or a
        finite-lag history jumps in realized throttle at a boundary (beyond
        ``1e-10``).

    Notes
    -----
    A finite time constant makes the realized throttle continuous across
    command discontinuities; a zero time constant follows the command and may
    jump.  The model is a transient-simulation approximation [1]_ [2]_.

    References
    ----------
    .. [1] NASA-TM-107318. NASA, "RL10A-3-3A Rocket Engine Modeling Project",
           https://ntrs.nasa.gov/citations/19970010379.
    .. [2] NASA-20040000363. NASA, "Transient Mathematical Modeling for Liquid
           Rocket Engine Systems", https://ntrs.nasa.gov/citations/20040000363.
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
        """End time of the last segment, in s."""
        return self.segments[-1].end_time_s

    @property
    def exposure_s(self) -> float:
        """Total throttle exposure ``integral(r(t) dt)``, in s."""
        return sum(segment.exposure_s for segment in self.segments)

    @property
    def throttle_first_moment_s2(self) -> float:
        """Total absolute-time first moment ``integral(t r(t) dt)``, in s**2."""
        return sum(segment.throttle_first_moment_s2 for segment in self.segments)

    @property
    def centroid_time_s(self) -> float:
        """Exposure-weighted centroid time, in s; zero when the exposure is zero."""
        return (
            self.throttle_first_moment_s2 / self.exposure_s
            if self.exposure_s > 0.0
            else 0.0
        )

    def realized_throttle(self, time_s: float, *, side: str = "right") -> float:
        """Evaluate the realized throttle with an explicit side at boundaries.

        Parameters
        ----------
        time_s : float
            Command time, in s.  Finite and non-negative.
        side : {"right", "left"}, optional
            Limit taken at a segment boundary (within a relative tolerance of
            ``1e-12`` of the duration, at least ``1e-12`` s).  At the end of the
            history ``"left"`` returns the last realized value and
            ``"right"`` returns zero.

        Returns
        -------
        float
            Realized throttle; zero after the history ends.

        Raises
        ------
        DomainError
            If ``side`` is invalid or ``time_s`` is negative or non-finite.
        """

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
        """Return a sampling grid that contains every segment boundary.

        Parameters
        ----------
        maximum_step_s : float
            Largest spacing between samples inside a segment, in s.

        Returns
        -------
        tuple of float
            Sorted, unique sample times from zero to the history duration.

        Raises
        ------
        DomainError
            If ``maximum_step_s`` is not finite and positive, or if the grid
            would exceed ``MAXIMUM_RESPONSE_SAMPLES`` samples.
        """

        if not isfinite(maximum_step_s) or maximum_step_s <= 0.0:
            raise DomainError("Maximum response sample step must be positive and finite.")
        # Count before building: a 1e9 s burn at a 1e-300 s step must be refused, not looped over.
        requested = sum(segment.duration_s / maximum_step_s for segment in self.segments)
        if not requested <= MAXIMUM_RESPONSE_SAMPLES:
            raise DomainError(
                f"Response grid needs about {requested:.3g} samples; the limit is "
                f"{MAXIMUM_RESPONSE_SAMPLES}. Increase maximum_step_s or shorten the burn.",
                code="sample_budget_exceeded", field="maximum_step_s")
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
    """Return the first-order lag value at the end of one linear command segment.

    Evaluates ``c0 + k (D - tau) + A exp(-D / tau)`` with
    ``A = r0 - c0 + k tau`` and clips it to ``[0, 1]``; with ``tau = 0`` the
    command end value is returned.  The lag is the transient approximation of
    [1]_.

    References
    ----------
    .. [1] NASA-TM-107318. NASA, "RL10A-3-3A Rocket Engine Modeling Project",
           https://ntrs.nasa.gov/citations/19970010379.
    """
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
    """Resolve a throttle schedule through ``tau dr/dt + r = command`` in closed form.

    Parameters
    ----------
    schedule : ThrottleSchedule
        Piecewise-linear command schedule.
    time_constant_s : float
        First-order time constant ``tau``, in s.  Finite and non-negative.
    ignition_delay_s : float, optional
        Zero-command interval inserted before the schedule, in s.
    settling_duration_s : float, optional
        Zero-command interval appended after the schedule, in s, during
        which the lagged throttle decays.
    initial_throttle : float, optional
        Realized throttle at time zero, in ``[0, 1]``.

    Returns
    -------
    RealizedThrottleHistory
        Segments for the ignition delay (if any), every resolved schedule
        interval, and the shutdown settling interval (if any).

    Raises
    ------
    DomainError
        If a time input is negative or non-finite, or ``initial_throttle`` is
        outside ``[0, 1]``.

    Notes
    -----
    The equation is a declared engineering approximation of start, throttle
    and shutdown transients [1]_ [2]_ [3]_, not a universal engine law.  Its
    time constant must be identified from hardware data or a higher-fidelity
    component model.  With ``tau = 0`` each segment starts at its command
    value.

    References
    ----------
    .. [1] NASA-TM-107318. NASA, "RL10A-3-3A Rocket Engine Modeling Project",
           https://ntrs.nasa.gov/citations/19970010379.
    .. [2] NASA-20040000363. NASA, "Transient Mathematical Modeling for Liquid
           Rocket Engine Systems", https://ntrs.nasa.gov/citations/20040000363.
    .. [3] NASA CR-184099. NASA, transient engine simulation covering steady
           balance, throttle transients, starts and shutdowns,
           https://ntrs.nasa.gov/citations/19910011919.
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
