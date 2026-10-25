"""Piecewise-linear throttle schedules with exact rational integrals.

A throttle schedule is a sequence of segments on each of which the command
``u`` varies linearly in time.  This module evaluates the throttle exposure
``E(t) = integral_0^t u dtau``, its first time moment
``M(t) = integral_0^t tau u dtau`` and the earliest time at which ``E``
reaches a requested value.  All arithmetic is exact on the rational values of
the supplied binary64 inputs [1]_; each exported field is rounded once.  Under
the proportional-throttle assumption the delivered impulse and propellant use
are ``F_rated E(t)`` and ``mdot_rated E(t)`` [2]_.

References
----------
.. [1] PYTHON-3.12-FRACTION. Python 3.12, "fractions" module.
       https://docs.python.org/3.12/library/fractions.html
.. [2] NASA Glenn Research Center, "Specific Impulse".
       https://www1.grc.nasa.gov/beginners-guide-to-aeronautics/specific-impulse/
"""

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
    """One linear throttle interval with binary64 endpoint coordinates.

    Parameters
    ----------
    start_time_s, end_time_s : float
        Interval endpoints with ``0 <= start_time_s <= end_time_s``, s.
    start_throttle, end_throttle : float
        Throttle command at the endpoints, each in ``[0, 1]``.
    phase : str
        Non-blank phase label.

    Raises
    ------
    DomainError
        If a coordinate is not a finite real number, the times are out of
        order or negative, a throttle lies outside ``[0, 1]``, or the phase is
        blank.

    Notes
    -----
    The derived integrals are evaluated exactly on the rational endpoint values
    and rounded once.  Multiplied by the rated thrust they give the interval
    impulse and its first time moment [1]_.

    References
    ----------
    .. [1] NASA Glenn Research Center, "Specific Impulse".
           https://www1.grc.nasa.gov/beginners-guide-to-aeronautics/specific-impulse/
    """

    start_time_s: float
    end_time_s: float
    start_throttle: float
    end_throttle: float
    phase: str

    def __post_init__(self):
        """Coerce endpoint fields to finite reals and validate the interval."""
        for name in ("start_time_s", "end_time_s", "start_throttle", "end_throttle"):
            object.__setattr__(self, name, real(getattr(self, name), f"Schedule interval {name}"))
        if self.start_time_s < 0 or self.end_time_s < self.start_time_s:
            raise DomainError("Schedule interval times are invalid.")
        if any(not 0 <= value <= 1 for value in (self.start_throttle, self.end_throttle)):
            raise DomainError("Schedule interval throttle is outside [0, 1].")
        text(self.phase, "Schedule interval phase")

    @property
    def duration_s(self) -> float:
        """Interval duration ``t1 - t0``, s, rounded once from the exact difference."""
        return _finite_rational(
            _fraction(self.end_time_s) - _fraction(self.start_time_s), "Interval duration"
        )

    @property
    def exposure_s(self) -> float:
        """Throttle exposure ``(u0 + u1) (t1 - t0) / 2`` of the interval, s.

        The exposure is the exact integral of the linear command; times the
        rated thrust it is the interval impulse [1]_.

        References
        ----------
        .. [1] NASA Glenn Research Center, "Specific Impulse".
               https://www1.grc.nasa.gov/beginners-guide-to-aeronautics/specific-impulse/
        """
        return _finite_rational(
            (_fraction(self.start_throttle) + _fraction(self.end_throttle))
            * (_fraction(self.end_time_s) - _fraction(self.start_time_s))
            / 2,
            "Interval exposure",
        )

    @property
    def throttle_first_moment_s2(self) -> float:
        """First time moment ``integral_t0^t1 t u dt`` of the interval, s^2.

        With ``h = t1 - t0`` the exact value is
        ``t0 (u0 + u1) h / 2 + h**2 (u0 / 2 + (u1 - u0) / 3)``.  Times the rated
        thrust it is the impulse moment from which a thrust-centroid time is
        formed [1]_.

        References
        ----------
        .. [1] NASA Glenn Research Center, "Specific Impulse".
               https://www1.grc.nasa.gov/beginners-guide-to-aeronautics/specific-impulse/
        """
        duration = _fraction(self.end_time_s) - _fraction(self.start_time_s)
        a, b = _fraction(self.start_throttle), _fraction(self.end_throttle)
        value = _fraction(self.start_time_s) * (a + b) * duration / 2 + duration**2 * (
            a / 2 + (b - a) / 3
        )
        return _finite_rational(value, "Interval first moment")


@dataclass(frozen=True, slots=True)
class _ScheduleSlice:
    """Exact rational slice of one schedule segment, possibly truncated.

    ``start`` and ``duration`` are in seconds; ``a`` and ``b`` are the command
    at the slice ends.
    """

    start: Fraction
    duration: Fraction
    a: Fraction
    b: Fraction
    phase: str

    @property
    def exposure(self):
        """Exact exposure ``(a + b) duration / 2`` of the slice, s."""
        return (self.a + self.b) * self.duration / 2

    @property
    def moment(self):
        """Exact first time moment of the slice about time zero, s^2."""
        return self.start * self.exposure + self.duration**2 * (self.a / 2 + (self.b - self.a) / 3)


def _schedule_totals(schedule):
    """Return the exact total duration and total exposure of a schedule.

    The totals are summed from the segment fields, so a query does not depend
    on whole-schedule scalar outputs being representable in binary64.

    Raises
    ------
    DomainError
        If ``schedule`` is not a `ThrottleSchedule`.
    """
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
    """Return the exact schedule slices up to ``stop_time``.

    ``stop_time`` defaults to the schedule end and is clipped to it.  A
    segment cut by the stop time is truncated, with its end command
    interpolated linearly.

    Raises
    ------
    DomainError
        If the schedule is invalid or ``stop_time`` is negative or not a finite
        real number.
    """
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
    """Return the exact exposure ``E(time)`` and first moment ``M(time)``."""
    slices = _schedule_slices(schedule, time)
    return (
        sum((s.exposure for s in slices), Fraction(0)),
        sum((s.moment for s in slices), Fraction(0)),
    )


def _sqrt_fraction(value):
    """Return the square root of a non-negative fraction.

    When numerator and denominator are both perfect squares (tested with
    ``math.isqrt`` [1]_) the root is returned exactly.  Otherwise the root is
    computed in a local 1100-digit decimal context [2]_ and converted exactly
    to a fraction [3]_; the global decimal context is not modified.

    References
    ----------
    .. [1] PYTHON-3.12-MATH. Python 3.12, "math" module, "math.isqrt".
           https://docs.python.org/3.12/library/math.html
    .. [2] PYTHON-3.12-DECIMAL. Python 3.12, "decimal" module.
           https://docs.python.org/3.12/library/decimal.html
    .. [3] PYTHON-3.12-FRACTION. Python 3.12, "fractions" module.
           https://docs.python.org/3.12/library/fractions.html
    """
    numerator, denominator = isqrt(value.numerator), isqrt(value.denominator)
    if numerator * numerator == value.numerator and denominator * denominator == value.denominator:
        return Fraction(numerator, denominator)
    with localcontext(Context(prec=1100)):
        return Fraction((Decimal(value.numerator) / Decimal(value.denominator)).sqrt())


def _inverse_exposure(schedule, exposure):
    """Return the earliest time at which the exposure reaches a request.

    Parameters
    ----------
    schedule : ThrottleSchedule
        Piecewise-linear throttle schedule.
    exposure : float or Fraction
        Requested non-negative exposure ``r``, s.

    Returns
    -------
    Fraction or None
        Earliest exact time ``t`` with ``E(t) = r``, or ``None`` when ``r``
        exceeds the total schedule exposure.

    Raises
    ------
    DomainError
        If the request is negative or not a finite real number, or no segment
        resolves a reachable request.

    Notes
    -----
    Exposures are compared exactly, without an absolute tolerance or a
    dimensional slope cutoff.  On the resolving segment, with start command
    ``a``, end command ``b``, duration ``h`` and remaining exposure ``r``, the
    local time solves ``a s + (b - a) s**2 / (2 h) = r``.  The root is taken in
    the cancellation-free form [1]_
    ``s = 2 r / (a + sqrt(a**2 + 2 (b - a) r / h))``, which also covers a
    zero start command; a constant command gives ``s = r / a``.  The local
    time is clipped to ``[0, h]``.

    References
    ----------
    .. [1] D. Goldberg, "What Every Computer Scientist Should Know About
           Floating-Point Arithmetic", ACM Computing Surveys 23(1), 1991.
    """
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
    """Return the schedule as binary64 intervals, optionally truncated.

    Parameters
    ----------
    schedule : ThrottleSchedule
        Piecewise-linear throttle schedule.
    stop_time_s : float or None, optional
        Truncation time, s; ``None`` keeps the full schedule.  Larger values
        are clipped to the schedule end.

    Returns
    -------
    tuple of ScheduleInterval
        One interval per (possibly truncated) segment.  Every endpoint is
        rounded once from its exact rational value.

    Raises
    ------
    DomainError
        If the schedule or stop time is invalid, or two distinct exact
        endpoints round to the same binary64 value (field
        ``"Interval time resolution"``).
    """
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
    """Return the throttle exposure ``E(t) = integral_0^t u dtau``.

    Parameters
    ----------
    schedule : ThrottleSchedule
        Piecewise-linear throttle schedule.
    time_s : float
        Query time, s; times beyond the schedule end give the total exposure.

    Returns
    -------
    float
        Exposure, s, rounded once from the exact value.  Under proportional
        throttling it is the delivered impulse divided by the rated thrust
        [1]_.

    Raises
    ------
    DomainError
        If the schedule or time is invalid or the result is not representable
        as a finite binary64 value.

    References
    ----------
    .. [1] NASA Glenn Research Center, "Specific Impulse".
           https://www1.grc.nasa.gov/beginners-guide-to-aeronautics/specific-impulse/
    """
    return _finite_rational(
        _schedule_integrals(schedule, _fraction(time_s, "Schedule query time"))[0],
        "Integrated throttle exposure",
    )


def throttle_first_moment(schedule: ThrottleSchedule, time_s: float) -> float:
    """Return the first time moment ``M(t) = integral_0^t tau u dtau``.

    Parameters
    ----------
    schedule : ThrottleSchedule
        Piecewise-linear throttle schedule.
    time_s : float
        Query time, s; times beyond the schedule end give the full moment.

    Returns
    -------
    float
        First moment, s^2, rounded once from the exact value.  ``M / E`` is the
        exposure-weighted (thrust-centroid) time.

    Raises
    ------
    DomainError
        If the schedule or time is invalid, or the exact moment overflows or
        underflows the finite binary64 range.
    """
    return _finite_rational(
        _schedule_integrals(schedule, _fraction(time_s, "Schedule query time"))[1],
        "Throttle first moment",
    )


def time_for_integrated_throttle(schedule: ThrottleSchedule, exposure_s: float) -> float | None:
    """Return the earliest time at which the exposure reaches a request.

    Parameters
    ----------
    schedule : ThrottleSchedule
        Piecewise-linear throttle schedule.
    exposure_s : float
        Requested non-negative exposure, s.

    Returns
    -------
    float or None
        Earliest time, s, rounded once from the exact root, or ``None`` when
        the request exceeds the total schedule exposure.

    Raises
    ------
    DomainError
        If the request is negative or not a finite real number.

    Notes
    -----
    Each segment is inverted in closed form with the cancellation-free
    quadratic root [1]_; see `_inverse_exposure`.

    References
    ----------
    .. [1] D. Goldberg, "What Every Computer Scientist Should Know About
           Floating-Point Arithmetic", ACM Computing Surveys 23(1), 1991.
    """
    value = _inverse_exposure(schedule, exposure_s)
    return None if value is None else _finite_rational(value, "Schedule inverse time")


@dataclass(frozen=True, slots=True)
class ProfiledPerformanceProvider:
    """L1 performance provider applying one throttle schedule to a rated point.

    Parameters
    ----------
    rated_point : EngineOperatingPoint
        Resolved operating point at full command.
    schedule : ThrottleSchedule
        Piecewise-linear throttle command.
    provider_id : str, default "profiled-performance"
        Provider identifier.
    provider_version : str, default "1.0"
        Provider version.

    Notes
    -----
    Thrust and every propellant stream scale in proportion to the command and
    system specific impulse is constant, so ``mdot = F / (g0 Isp)`` [1]_ holds
    at every command with ``g0`` the registered ``STANDARD_GRAVITY_M_S2``.
    This proportional scaling is a declared modelling assumption; the exact
    numerical accounting does not validate it.

    References
    ----------
    .. [1] NASA Glenn Research Center, "Specific Impulse".
           https://www1.grc.nasa.gov/beginners-guide-to-aeronautics/specific-impulse/
    """

    rated_point: EngineOperatingPoint
    schedule: ThrottleSchedule
    provider_id: str = "profiled-performance"
    provider_version: str = "1.0"

    def operating_point(self) -> EngineOperatingPoint:
        """Return the rated (full-command) operating point."""
        return self.rated_point

    def discontinuities_s(self, definition: BurnDefinition) -> tuple[float, ...]:
        """Return the interior segment-boundary times of the schedule.

        Parameters
        ----------
        definition : BurnDefinition
            Burn definition; accepted for interface compatibility and unused.

        Returns
        -------
        tuple of float
            Cumulative segment end times, s, excluding the final end, each
            rounded once from its exact sum.

        Raises
        ------
        DomainError
            If two consecutive boundaries round to the same binary64 value
            (field ``"Discontinuity time resolution"``).
        """
        del definition
        cursor = Fraction(0)
        boundaries: list[float] = []
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
