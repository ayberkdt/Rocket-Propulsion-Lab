"""Adversarial pairing of commanded and observed valve transitions."""

from dataclasses import replace
from itertools import product

import pytest

from rocket_propulsion.core.errors import DomainError
from rocket_propulsion.studies.sass import (
    ObservedValveTrace,
    ValveCommand,
    ValveLatencyLimits,
    review_valve_responses,
)
from tests.test_sass_command_life import history


def observed(case, *, error=0, verified=None):
    case = replace(case, source_record="test-only observed state ledger")
    return ObservedValveTrace(case, case.head_ids if verified is None else verified,
                             "test-only observation review", "test-only clock review", error)


def limits(opening=5, closing=3):
    return ValveLatencyLimits(opening, closing, "test-only delay requirements")


def test_time_order_and_delay_classification_against_all_possible_clock_offsets():
    # Enumerate every integer time in the stated uncertainty interval.
    # This oracle does not perform the production command-bin search.
    for command_tick, limit, error, observed_tick in product((3, 8), (0, 2, 4), (0, 1, 2), range(1, 12)):
        commanded = history(ValveCommand(command_tick, "A", True), duration=12, heads=("A",))
        actual = history(ValveCommand(observed_tick, "A", True), duration=12, heads=("A",))
        review = review_valve_responses(commanded, observed(actual, error=error), limits(limit))
        possible = list(range(observed_tick - error, observed_tick + error + 1))
        if any(t < 0 or t > 12 for t in possible) or (
            any(t < command_tick for t in possible) and any(t >= command_tick for t in possible)
        ):
            expected = "unassessed"
        elif all(t < command_tick for t in possible):
            expected = "fail"
            assert review.observation_issues[0].status == "unexpected_observation"
        else:
            delays = [t - command_tick for t in possible]
            expected = ("pass" if all(delay <= limit for delay in delays) else
                        "fail" if all(delay > limit for delay in delays) else "unassessed")
            assert review.responses[0].latency_bounds_ticks == (min(delays), max(delays))
        assert review.status == expected, (command_tick, limit, error, observed_tick)


def test_opening_and_closing_use_their_own_approved_limits():
    cmd = history(ValveCommand(5, "A", True), ValveCommand(12, "A", False), duration=20, heads=("A",))
    actual = history(ValveCommand(7, "A", True), ValveCommand(15, "A", False), duration=20, heads=("A",))
    review = review_valve_responses(cmd, observed(actual), limits(opening=4, closing=2))
    assert review.status == "fail"
    assert [row.status for row in review.responses] == ["within_limit", "late_response"]
    assert [row.latency_bounds_ticks for row in review.responses] == [(2, 2), (3, 3)]


def test_late_open_cannot_be_reassigned_across_a_superseding_close():
    cmd = history(ValveCommand(10, "A", True), ValveCommand(12, "A", False),
                  ValveCommand(20, "A", True), duration=30, heads=("A",))
    actual = history(ValveCommand(13, "A", True), ValveCommand(15, "A", False),
                     ValveCommand(22, "A", True), duration=30, heads=("A",))
    review = review_valve_responses(cmd, observed(actual), limits())
    assert review.status == "fail"
    assert [row.status for row in review.responses] == ["superseded", "within_limit", "within_limit"]
    assert [row.observation_tick for row in review.responses] == [None, 15, 22]
    assert review.observation_issues[0].observation.tick == 13
    assert review.observation_issues[0].status == "unexpected_observation"


def test_uncertain_transition_ownership_keeps_both_neighboring_commands_open():
    cmd = history(ValveCommand(5, "A", True), ValveCommand(7, "A", False),
                  ValveCommand(9, "A", True), duration=18, heads=("A",))
    actual = history(ValveCommand(8, "A", True), duration=18, heads=("A",))
    review = review_valve_responses(cmd, observed(actual, error=1), limits())
    assert review.status == "unassessed"
    assert [row.status for row in review.responses] == ["superseded", "ambiguous_timing", "ambiguous_timing"]
    assert review.observation_issues[0].status == "ambiguous_observation"


def test_multiple_matching_transitions_cannot_hide_oscillation_after_a_command():
    cmd = history(ValveCommand(5, "A", True), duration=12, heads=("A",))
    actual = history(ValveCommand(6, "A", True), ValveCommand(7, "A", False),
                     ValveCommand(8, "A", True), duration=12, heads=("A",))
    review = review_valve_responses(cmd, observed(actual), limits())
    assert review.status == "fail"
    assert review.responses[0].status == "multiple_responses"
    assert [issue.observation.tick for issue in review.observation_issues] == [7]


def test_absent_response_requires_verified_coverage_and_does_not_hide_other_failure():
    cmd = history(ValveCommand(10, "A", True), ValveCommand(10, "B", True), duration=20)
    actual = history(duration=20)
    review = review_valve_responses(cmd, observed(actual, verified=("A",)), limits())
    assert review.status == "fail"
    assert review.unverified_head_ids == ("B",)
    assert [row.status for row in review.responses] == ["missing_response", "coverage_unverified"]
    no_coverage = review_valve_responses(cmd, observed(actual, verified=()), limits())
    assert no_coverage.status == "unassessed"
    assert all(row.status == "coverage_unverified" for row in no_coverage.responses)


def test_observation_on_other_head_cannot_satisfy_missing_response():
    cmd = history(ValveCommand(10, "A", True), duration=20)
    actual = history(ValveCommand(11, "B", True), duration=20)
    review = review_valve_responses(cmd, observed(actual), limits())
    assert review.status == "fail"
    assert review.responses[0].status == "missing_response"
    assert review.observation_issues[0].observation.head_id == "B"


def test_unverified_idle_head_cannot_be_inferred_closed_from_silence():
    cmd = history(ValveCommand(3, "A", True), duration=20)
    actual = history(ValveCommand(4, "A", True), duration=20)
    complete = review_valve_responses(cmd, observed(actual), limits())
    assert complete.status == "pass"
    incomplete = review_valve_responses(cmd, observed(actual, verified=("A",)), limits())
    assert incomplete.responses[0].status == "within_limit"
    assert incomplete.status == "unassessed"
    assert incomplete.unverified_head_ids == ("B",)


@pytest.mark.parametrize(("tick", "duration", "error", "status"), [
    (15, 20, 0, "missing_response"),
    (15, 19, 0, "right_censored"),
    (15, 20, 1, "right_censored"),
    (0, 20, 1, "left_censored"),
])
def test_window_boundary_distinguishes_missing_from_unobservable(tick, duration, error, status):
    cmd = history(ValveCommand(tick, "A", True), duration=duration, heads=("A",))
    review = review_valve_responses(cmd, observed(history(duration=duration, heads=("A",)), error=error),
                                    limits())
    assert review.responses[0].status == status
    assert review.status == ("fail" if status == "missing_response" else "unassessed")


def test_left_boundary_state_disagreement_cannot_be_promoted_to_hardware_failure_or_pass():
    cmd = history(ValveCommand(10, "A", False), duration=20, heads=("A",), initially_open=("A",))
    actual = history(ValveCommand(2, "A", True), ValveCommand(11, "A", False), duration=20, heads=("A",))
    review = review_valve_responses(cmd, observed(actual), limits())
    assert review.status == "unassessed"
    assert review.initial_state_mismatch_head_ids == ("A",)
    assert review.observation_issues[0].status == "ambiguous_observation"
    assert review.responses[0].status == "within_limit"


def test_mismatched_clocks_and_invalid_observed_transitions_are_rejected():
    cmd = history(ValveCommand(3, "A", True), heads=("A",))
    for actual in (
        history(heads=("A",), tick_seconds=2.0),
        history(heads=("A",), duration=9),
        history(heads=("B",)),
        history(ValveCommand(3, "A", True), ValveCommand(2, "A", False), heads=("A",)),
    ):
        with pytest.raises(DomainError):
            review_valve_responses(cmd, observed(actual), limits())
