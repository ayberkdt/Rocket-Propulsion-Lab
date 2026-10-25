"""Restart, cooldown, and minimum-pulse validation for engine commands.

The module enforces caller-declared, hardware-specific timing limits on a
multi-engine pulse train and converts each engine's validated pulses into a
:class:`~rocket_propulsion.propulsion.burns.models.ThrottleSchedule`.  It
ships no generic minimum-pulse, cooldown or restart values: pulse rise,
tailoff and minimum impulse bit depend on the valve and thruster [1]_ [2]_,
and start count is a hardware reliability driver [3]_.

References
----------
.. [1] NASA-CR-122347. NASA, "Monopropellant Hydrazine Resisto Jet",
       https://ntrs.nasa.gov/citations/19720011120.
.. [2] NASA-20070032005. J. M. Parker and M. J. Wilson, "The Minimum Impulse
       Thruster", 2005, https://ntrs.nasa.gov/citations/20070032005.
.. [3] NASA-20050207429. NASA, "Key Reliability Drivers of Liquid Propulsion
       Engines and a Reliability Model for Sensitivity Analysis",
       https://ntrs.nasa.gov/citations/20050207429.
"""

from __future__ import annotations

from dataclasses import dataclass
from itertools import pairwise
from math import isfinite

from rocket_propulsion.core.errors import DomainError

from .models import ThrottleSchedule, ThrottleSegment


def _finite_nonnegative(value: float, label: str) -> None:
    """Raise ``DomainError`` unless ``value`` is finite and non-negative."""
    if not isfinite(value) or value < 0.0:
        raise DomainError(f"{label} must be finite and non-negative.")


@dataclass(frozen=True, slots=True)
class EngineDutyCycleLimits:
    """Operational timing limits for one independently commanded engine.

    Parameters
    ----------
    engine_id : str
        Non-empty engine identifier.
    minimum_on_time_s : float, optional
        Shortest permitted pulse, in s.  Finite and non-negative.
    minimum_off_time_s : float, optional
        Shortest permitted cooldown between consecutive pulses, in s.
    maximum_starts : int or None, optional
        Largest permitted number of pulses; ``None`` means unlimited.
    ignition_delay_s : float, optional
        Command-to-ignition delay carried for downstream response models, s.
    response_time_constant_s : float, optional
        First-order response time constant carried for downstream response
        models, in s.

    Raises
    ------
    DomainError
        If the identifier is blank, a time is negative or non-finite, or
        ``maximum_starts`` is declared below one.

    Notes
    -----
    All limits are hardware-specific inputs.  Minimum on-time reflects the
    valve response and minimum impulse bit of the selected thruster [1]_ [2]_;
    the start count is a reliability driver of the selected engine [3]_.

    References
    ----------
    .. [1] NASA-CR-122347. NASA, "Monopropellant Hydrazine Resisto Jet",
           https://ntrs.nasa.gov/citations/19720011120.
    .. [2] NASA-20070032005. J. M. Parker and M. J. Wilson, "The Minimum
           Impulse Thruster", 2005, https://ntrs.nasa.gov/citations/20070032005.
    .. [3] NASA-20050207429. NASA, "Key Reliability Drivers of Liquid
           Propulsion Engines and a Reliability Model for Sensitivity
           Analysis", https://ntrs.nasa.gov/citations/20050207429.
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

    Parameters
    ----------
    engine_id : str
        Non-empty identifier of the commanded engine.
    start_time_s : float
        Command start time, in s.  Finite and non-negative.
    duration_s : float
        Commanded on-time, in s.  Finite and positive.
    throttle : float, optional
        Constant commanded throttle in ``(0, 1]``.

    Raises
    ------
    DomainError
        If the identifier is blank or a field is outside its domain.

    Notes
    -----
    The pulse is a command interval.  The realized thrust history, including
    rise and tailoff [1]_, is produced by a response model.

    References
    ----------
    .. [1] NASA-CR-122347. NASA, "Monopropellant Hydrazine Resisto Jet",
           https://ntrs.nasa.gov/citations/19720011120.
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
        """Command end time ``start_time_s + duration_s``, in s."""
        return self.start_time_s + self.duration_s


@dataclass(frozen=True, slots=True)
class EngineSequenceValidation:
    """Normalized, validated multi-engine command sequence.

    Instances are returned by :func:`validate_engine_sequence`.

    Parameters
    ----------
    pulses : tuple of EnginePulse
        Pulses sorted by engine identifier, then start time.
    starts_by_engine : tuple of (str, int)
        Pulse count of each engine, sorted by identifier.
    duration_s : float
        Latest pulse end time, in s.

    References
    ----------
    .. [1] NASA-CR-122347. NASA, "Monopropellant Hydrazine Resisto Jet",
           https://ntrs.nasa.gov/citations/19720011120.
    """

    pulses: tuple[EnginePulse, ...]
    starts_by_engine: tuple[tuple[str, int], ...]
    duration_s: float

    def schedule_for_engine(
        self, engine_id: str, *, total_duration_s: float | None = None
    ) -> ThrottleSchedule:
        """Convert one engine's pulse train into a throttle schedule.

        Gaps before, between and after pulses become explicit zero-throttle
        segments, so each on/off discontinuity is a segment boundary.

        Parameters
        ----------
        engine_id : str
            Engine whose pulses are converted.
        total_duration_s : float or None, optional
            Schedule length, in s; defaults to the sequence duration.  Must
            cover every selected pulse.

        Returns
        -------
        ThrottleSchedule
            Piecewise-constant schedule starting at time zero.

        Raises
        ------
        DomainError
            If the engine has no validated pulses or the duration is
            non-finite or shorter than its last pulse.
        """

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
    """Validate pulse timing rules and return a deterministic sequence.

    Parameters
    ----------
    pulses : tuple of EnginePulse
        At least one commanded pulse, in any order.
    limits : tuple of EngineDutyCycleLimits
        One limit record per engine; identifiers must be unique and cover
        every engine that appears in ``pulses``.

    Returns
    -------
    EngineSequenceValidation
        Pulses sorted by engine and start time, with start counts.

    Raises
    ------
    DomainError
        If either input is empty, limits are duplicated or missing for an
        engine, an engine exceeds its start count, a pulse is shorter than the
        minimum on-time, pulses of one engine overlap, or a cooldown gap is
        shorter than the minimum off-time.

    Notes
    -----
    Comparisons use an absolute tolerance of ``1e-12`` s.  Limits are
    hardware-specific constraints [1]_ [2]_; this function only enforces them
    and does not supply generic minimum-pulse or cooldown values.

    References
    ----------
    .. [1] NASA-CR-122347. NASA, "Monopropellant Hydrazine Resisto Jet",
           https://ntrs.nasa.gov/citations/19720011120.
    .. [2] NASA-20070032005. J. M. Parker and M. J. Wilson, "The Minimum
           Impulse Thruster", 2005, https://ntrs.nasa.gov/citations/20070032005.
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
