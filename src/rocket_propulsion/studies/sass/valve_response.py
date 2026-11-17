"""Compare issued valve commands and observed transitions on a common clock."""

from __future__ import annotations

from bisect import bisect_left, bisect_right
from dataclasses import dataclass
from typing import Literal

from rocket_propulsion.core.errors import DomainError
from rocket_propulsion.core.validation import integer, strings, text

from .command_life import ThrusterCommandHistory, ValveCommand, _command_usage
from .evidence import CitedScalar


@dataclass(frozen=True, slots=True)
class ObservedValveTrace:
    """Observed valve state and review records for clock alignment and coverage.

    The history stores transitions of observed states. A named head has
    coverage only if the reviewer established continuous observation and
    reliable detection over the entire declared window. Timing error is
    an absolute relative clock and timestamp error bound in common ticks.

    References
    ----------
    .. [1] NASA-7009B, 4.2.5 and 4.2.7, empirical referents and uncertainty.
    """

    history: ThrusterCommandHistory
    verified_head_ids: tuple[str, ...]
    coverage_review_record: str
    clock_alignment_record: str
    max_relative_timing_error_ticks: int

    def __post_init__(self) -> None:
        """Freeze the claimed coverage and integer timing error.

        References
        ----------
        .. [1] NASA-7009B, 4.2.5 and 4.2.7, documented referents and uncertainty.
        """
        if not isinstance(self.history, ThrusterCommandHistory):
            raise DomainError("Observed trace requires a valve state history.")
        heads = strings(self.verified_head_ids, "verified observation heads")
        if len(set(heads)) != len(heads) or not set(heads).issubset(self.history.head_ids):
            raise DomainError("Verified observation heads must be unique observed head IDs.")
        error = integer(self.max_relative_timing_error_ticks, "maximum relative timing error ticks")
        if error < 0:
            raise DomainError("The relative timing error bound must be nonnegative.")
        object.__setattr__(self, "verified_head_ids", heads)
        object.__setattr__(self, "max_relative_timing_error_ticks", error)
        object.__setattr__(self, "coverage_review_record", text(self.coverage_review_record, "coverage review"))
        object.__setattr__(self, "clock_alignment_record", text(self.clock_alignment_record, "clock alignment"))


@dataclass(frozen=True, slots=True)
class ValveLatencyLimits:
    """Approved upper opening and closing delays on the common tick scale.

    References
    ----------
    .. [1] NASA-7009B, 4.1.1.5, documented acceptance criteria.
    """

    opening_ticks: int
    closing_ticks: int
    approval_record: str

    def __post_init__(self) -> None:
        """Require finite window-independent integer limits and an approval reference.

        References
        ----------
        .. [1] NASA-7009B, 4.1.1.5, documented acceptance criteria.
        """
        for field in ("opening_ticks", "closing_ticks"):
            delay = integer(getattr(self, field), field)
            if delay < 0:
                raise DomainError("Valve latency limits must be nonnegative integer ticks.")
            object.__setattr__(self, field, delay)
        object.__setattr__(self, "approval_record", text(self.approval_record, "latency approval"))


@dataclass(frozen=True, slots=True)
class ValveCommandResponse:
    """Result for one issued transition, with a latency enclosure when matched.

    References
    ----------
    .. [1] NASA-7009B, 4.2.7 and 4.3.8, uncertainty and results reporting.
    """

    command: ValveCommand
    limit_ticks: int
    status: str
    observation_tick: int | None = None
    latency_bounds_ticks: tuple[int, int] | None = None


@dataclass(frozen=True, slots=True)
class ValveObservationIssue:
    """An observed transition without a unique matching command.

    References
    ----------
    .. [1] NASA-7009B, 4.2.5, retaining empirical differences.
    """

    observation: ValveCommand
    status: Literal["unexpected_observation", "ambiguous_observation"]


@dataclass(frozen=True, slots=True)
class ValveResponseReview:
    """Response evidence for the declared observation window and limits.

    References
    ----------
    .. [1] NASA-7009B, 4.3.8, reporting the scope of results.
    """

    status: Literal["pass", "fail", "unassessed"]
    command_source_record: str
    observation_source_record: str
    coverage_review_record: str
    clock_alignment_record: str
    max_relative_timing_error_ticks: int
    limits: ValveLatencyLimits
    duration_ticks: int
    tick_period: CitedScalar
    unverified_head_ids: tuple[str, ...]
    initial_state_mismatch_head_ids: tuple[str, ...]
    responses: tuple[ValveCommandResponse, ...]
    observation_issues: tuple[ValveObservationIssue, ...]


def _group_by_head(history: ThrusterCommandHistory) -> dict[str, list[ValveCommand]]:
    """Validate the complete state ledger and retain each head's transition order.

    References
    ----------
    .. [1] NASA-SE-6. NASA Systems Engineering Handbook, 6.6, controlled data.
    """
    _command_usage(history)
    result: dict[str, list[ValveCommand]] = {head: [] for head in history.head_ids}
    for event in history.commands:
        result[event.head_id].append(event)
    return result


def _review_head_responses(
    commands: list[ValveCommand], observations: list[ValveCommand],
    initially_open: bool, initial_mismatch: bool, duration: int, timing_error: int, limits: ValveLatencyLimits,
) -> tuple[list[ValveCommandResponse], list[ValveObservationIssue]]:
    """Find responses that belong unambiguously to one command interval.

    References
    ----------
    .. [1] NASA-7009B, 4.2.5 and 4.2.7, comparisons and uncertainty.
    """
    ticks = [event.tick for event in commands]
    observation_ticks = [event.tick for event in observations]
    matches: list[list[ValveCommand]] = [[] for _ in commands]
    ambiguity = [0] * (len(commands) + 1)
    issues = []
    for observation in observations:
        lower, upper = observation.tick - timing_error, observation.tick + timing_error
        first = bisect_right(ticks, lower) - 1
        last = bisect_right(ticks, upper) - 1
        if first != last or lower < 0 or upper > duration:
            # Difference counts cover all affected bins in linear total work.
            if last >= 0:
                ambiguity[max(first, 0)] += 1
                ambiguity[last + 1] -= 1
            issues.append(ValveObservationIssue(observation, "ambiguous_observation"))
        elif first < 0 and initial_mismatch:
            issues.append(ValveObservationIssue(observation, "ambiguous_observation"))
        elif first < 0 or commands[first].opened != observation.opened:
            issues.append(ValveObservationIssue(observation, "unexpected_observation"))
        else:
            matches[first].append(observation)
    responses = []
    active_ambiguity = 0
    for index, command in enumerate(commands):
        active_ambiguity += ambiguity[index]
        limit = limits.opening_ticks if command.opened else limits.closing_ticks
        matching = matches[index]
        if len(matching) > 1:
            responses.append(ValveCommandResponse(command, limit, "multiple_responses"))
        elif active_ambiguity:
            responses.append(ValveCommandResponse(command, limit, "ambiguous_timing"))
        elif matching:
            tick = matching[0].tick
            delay = tick - command.tick
            bounds = (delay - timing_error, delay + timing_error)
            status = ("late_response" if bounds[0] > limit else
                      "within_limit" if bounds[1] <= limit else "timing_indeterminate")
            responses.append(ValveCommandResponse(command, limit, status, tick, bounds))
        else:
            last_before = bisect_left(observation_ticks, command.tick - timing_error) - 1
            prior_open = initially_open if last_before < 0 else observations[last_before].opened
            deadline = command.tick + limit + timing_error
            if command.tick < timing_error:
                status = "left_censored"
            elif prior_open == command.opened:
                status = "already_in_target_state"
            elif index + 1 < len(ticks) and deadline >= ticks[index + 1]:
                status = "superseded"
            elif deadline > duration:
                status = "right_censored"
            else:
                status = "missing_response"
            responses.append(ValveCommandResponse(command, limit, status))
    return responses, issues


def review_valve_responses(
    commanded: ThrusterCommandHistory, observed: ObservedValveTrace, limits: ValveLatencyLimits,
) -> ValveResponseReview:
    """Assess per-head timing against approved opening and closing limits.

    Both traces must describe the same time window, installed head set and
    tick scale. The caller provides reviewed clock alignment and observation
    coverage. The review assigns observed transitions only to the command
    interval they provably occupy. Repeated, unexpected or definitely late
    responses fail. A missing response fails only when the complete deadline
    lies inside the observed interval before any superseding command.

    References
    ----------
    .. [1] NASA-7009B, 4.1.1.5 and 4.2.5–4.2.7.
    .. [2] NASA NESC, NASA/TM-20210022275, transient pressure event verification.
    """
    if not (isinstance(commanded, ThrusterCommandHistory) and isinstance(observed, ObservedValveTrace)
            and isinstance(limits, ValveLatencyLimits)):
        raise DomainError("Valve response review requires commands, an observed trace and approved limits.")
    measured = observed.history
    if (commanded.duration_ticks != measured.duration_ticks
            or commanded.tick_period.value != measured.tick_period.value
            or set(commanded.head_ids) != set(measured.head_ids)):
        raise DomainError("Command and observation histories must share head IDs, clock scale and window.")
    command_by_head = _group_by_head(commanded)
    measured_by_head = _group_by_head(measured)
    verified = set(observed.verified_head_ids)
    unverified = tuple(head for head in commanded.head_ids if head not in verified)
    command_initial = set(commanded.initially_open)
    observed_initial = set(measured.initially_open)
    initial_mismatch = tuple(head for head in commanded.head_ids if head in verified
                             and ((head in command_initial) != (head in observed_initial)))
    responses: list[ValveCommandResponse] = []
    issues: list[ValveObservationIssue] = []
    for head in commanded.head_ids:
        if head not in verified:
            responses.extend(ValveCommandResponse(
                cmd, limits.opening_ticks if cmd.opened else limits.closing_ticks, "coverage_unverified",
            ) for cmd in command_by_head[head])
            continue
        head_responses, head_issues = _review_head_responses(
            command_by_head[head], measured_by_head[head], head in observed_initial,
            head in initial_mismatch, commanded.duration_ticks,
            observed.max_relative_timing_error_ticks, limits,
        )
        responses.extend(head_responses)
        issues.extend(head_issues)
    failed = any(row.status in {"late_response", "missing_response", "multiple_responses"}
                 for row in responses) or any(row.status == "unexpected_observation" for row in issues)
    passed = bool(responses) and not (unverified or initial_mismatch or issues) and all(
        row.status == "within_limit" for row in responses)
    return ValveResponseReview(
        "fail" if failed else "pass" if passed else "unassessed",
        commanded.source_record, measured.source_record, observed.coverage_review_record,
        observed.clock_alignment_record, observed.max_relative_timing_error_ticks, limits,
        commanded.duration_ticks, commanded.tick_period, unverified, initial_mismatch,
        tuple(responses), tuple(issues),
    )
