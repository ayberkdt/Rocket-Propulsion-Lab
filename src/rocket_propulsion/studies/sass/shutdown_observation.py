"""Tie residual forward impulse to independently observed valve closures."""

from __future__ import annotations

from bisect import bisect_right
from dataclasses import dataclass

from rocket_propulsion.core.errors import DomainError
from rocket_propulsion.core.validation import integer, records, text

from .command_life import ThrusterCommandHistory, ValveCommand
from .transient_observation import TransientReview, TransientTrace, review_transient
from .valve_response import (
    ObservedValveTrace,
    ValveCommandResponse,
    ValveLatencyLimits,
    ValveResponseReview,
    review_valve_responses,
)


@dataclass(frozen=True, slots=True)
class ShutdownWindow:
    """The reviewed end and force-to-command clock error for one shutdown.

    The force and command records share a nominal time origin and tick scale.
    The bound covers residual clock offset, drift and timestamp quantization
    over this window. The valve-state trace supplies its own alignment bound.

    References
    ----------
    .. [1] NASA-7009B, 4.2.5 and 4.2.7, trace alignment and uncertainty.
    """

    tail_end_tick: int
    force_to_command_error_ticks: int
    clock_alignment_record: str
    tail_window_review_record: str

    def __post_init__(self) -> None:
        """Require explicit integer bounds and review records.

        References
        ----------
        .. [1] NASA-7009B, 4.2.1 and 4.2.7, recorded conditions and uncertainty.
        """
        end = integer(self.tail_end_tick, "shutdown tail end tick")
        error = integer(self.force_to_command_error_ticks, "force-to-command timing error")
        if end <= 0 or error < 0:
            raise DomainError("Shutdown end must be positive and timing error nonnegative.")
        object.__setattr__(self, "tail_end_tick", end)
        object.__setattr__(self, "force_to_command_error_ticks", error)
        text(self.clock_alignment_record, "force-to-command clock alignment")
        text(self.tail_window_review_record, "shutdown tail window review")


@dataclass(frozen=True, slots=True)
class ShutdownImpulseReview:
    """Bounded positive impulse or the reason closure-based reduction was refused.

    The reported valve rows retain their timing assessment. An identifiable
    late closure can anchor an impulse measurement even when timing fails.

    References
    ----------
    .. [1] NASA-7009B, 4.3.8, scope and limitations of results.
    """

    status: str
    window: ShutdownWindow
    valve_review: ValveResponseReview
    closures: tuple[ValveCommandResponse, ...]
    final_closure_ticks: tuple[int, int] | None = None
    thrust_review: TransientReview | None = None


def review_shutdown_impulse(
    commanded: ThrusterCommandHistory, observed: ObservedValveTrace, limits: ValveLatencyLimits,
    thrust: TransientTrace, close_commands: tuple[ValveCommand, ...], window: ShutdownWindow,
) -> ShutdownImpulseReview:
    """Measure forward impulse after the last observed closure in a shutdown.

    Select one issued close for each valve in the declared shutdown. The
    reviewer establishes the relationship between the force channel, the
    selected valves and any other active hardware. The tail-end review
    establishes the analysis endpoint. The result includes its timing scope.

    References
    ----------
    .. [1] NASA-7009B, 4.2.5 and 4.2.7, empirical channels and uncertainty.
    .. [2] NASA NESC, NASA/TM-20210022275, transient event measurement.
    """
    if not isinstance(window, ShutdownWindow) or not isinstance(thrust, TransientTrace) or thrust.unit != "N":
        raise DomainError("Shutdown impulse requires a reviewed window and a thrust trace in N.")
    commands: tuple[ValveCommand, ...] = records(close_commands, ValveCommand, "shutdown closes")
    heads = {command.head_id for command in commands}
    if not commands or len(heads) != len(commands) or any(command.opened for command in commands):
        raise DomainError("Select one issued closing command per shutdown valve.")
    valves = review_valve_responses(commanded, observed, limits)
    if thrust.tick_period.value != valves.tick_period.value:
        raise DomainError("Shutdown thrust and valve histories must use the same tick scale.")
    responses = {row.command: row for row in valves.responses}
    if not all(command in responses for command in commands):
        raise DomainError("Selected shutdown closes must be present in the command history.")
    closures = tuple(responses[command] for command in commands)
    context = (window, valves, closures)
    ticks = {row.command.head_id: row.observation_tick for row in closures
             if row.observation_tick is not None}
    if len(ticks) != len(closures) or heads.intersection(valves.initial_state_mismatch_head_ids):
        return ShutdownImpulseReview("closure_unresolved", *context)

    sync_error = window.force_to_command_error_ticks
    total_error = valves.max_relative_timing_error_ticks + sync_error
    last = max(ticks.values())
    final = (last - total_error, last + total_error)
    end = window.tail_end_tick
    if final[0] < 0 or end + total_error > valves.duration_ticks or final[1] >= end:
        return ShutdownImpulseReview("tail_window_not_covered", *context, final)

    selected = {command.head_id: command.tick for command in commands}
    if any(event.opened and event.head_id in heads and selected[event.head_id] < event.tick
           and event.tick - sync_error <= end for event in commanded.commands):
        return ShutdownImpulseReview("subsequent_opening_command", *context, final)
    if any(event.opened and event.head_id in heads and ticks[event.head_id] < event.tick
           and event.tick - total_error <= end for event in observed.history.commands):
        return ShutdownImpulseReview("possible_valve_reopening", *context, final)
    if any(issue.observation.head_id in heads
           and final[0] - total_error <= issue.observation.tick <= end + total_error
           for issue in valves.observation_issues):
        return ShutdownImpulseReview("closure_unresolved", *context, final)

    sample_ticks = [sample.tick for sample in thrust.samples]
    left = bisect_right(sample_ticks, final[0]) - 1
    if left < 0:
        return ShutdownImpulseReview("missing_window_boundary", *context, final)
    result = review_transient(thrust, sample_ticks[left], end,
                              window_review_record=window.tail_window_review_record,
                              positive_impulse_start_ticks=final)
    return ShutdownImpulseReview(result.status, *context, final, result)

