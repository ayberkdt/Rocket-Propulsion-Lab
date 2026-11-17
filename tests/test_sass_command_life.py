"""Adversarial per-head timing, censoring and project admission checks."""

import json
import random
from dataclasses import asdict, replace

import pytest

from rocket_propulsion.core.errors import DomainError
from rocket_propulsion.studies.sass import (
    ThrusterCommandHistory,
    ValveCommand,
    build_project_html,
    project_from_json,
    review_thruster_command_life,
    thruster,
)
from tests.test_sass_studies import coherent_life_catalogue, life_study, q


def history(*commands, duration=8, heads=("A", "B"), initially_open=(), tick_seconds=1.0):
    return ThrusterCommandHistory("test-only explicit command schedule", q(tick_seconds, "s"),
                                  duration, heads, initially_open, commands)


def review(case):
    return review_thruster_command_life(thruster("ariane-1n"), case)


def test_atomic_head_swap_and_boundary_censoring():
    case = history(ValveCommand(3, "B", True), ValveCommand(3, "A", False),
                   ValveCommand(7, "B", False), ValveCommand(7, "A", True),
                   duration=8, initially_open=("A",))
    result = review(case)
    assert result.peak_simultaneous_open == 1
    a, b = result.usage
    assert (a.observed_open_ticks, a.longest_open_ticks, a.completed_cycles) == (4, 3, 0)
    assert a.left_censored and a.right_censored
    assert (b.observed_open_ticks, b.longest_open_ticks, b.completed_cycles) == (4, 4, 1)
    assert not b.left_censored and not b.right_censored
    assert [row.status for row in result.comparisons if row.head_id == "A"] == ["unassessed"] * 3
    reordered = (case.commands[1], case.commands[0], case.commands[3], case.commands[2])
    assert review(replace(case, commands=reordered)) == result


def test_closing_at_final_tick_completes_cycle_without_instantaneous_peak():
    case = history(ValveCommand(0, "A", True), ValveCommand(8, "B", True),
                   ValveCommand(8, "A", False))
    result = review(case)
    assert result.peak_simultaneous_open == 1
    a, b = result.usage
    assert (a.observed_open_ticks, a.completed_cycles, a.right_censored) == (8, 1, False)
    assert (b.observed_open_ticks, b.completed_cycles, b.right_censored) == (0, 0, True)


def test_tick_ledger_against_independent_discrete_occupancy_oracle():
    random_state = random.Random(26052926)
    heads = ("A", "B", "C", "D")
    for _ in range(80):
        duration = random_state.randrange(1, 40)
        before = {head: bool(random_state.getrandbits(1)) for head in heads}
        # Each tick's occupancy and the terminal state are generated first.
        # The command stream is a derived encoding of this independent oracle.
        states = {head: [bool(random_state.getrandbits(1)) for _ in range(duration + 1)]
                  for head in heads}
        events = []
        for tick in range(duration + 1):
            for head in heads:
                old = before[head] if tick == 0 else states[head][tick - 1]
                if old != states[head][tick]:
                    events.append(ValveCommand(tick, head, states[head][tick]))
        result = review(history(*events, duration=duration, heads=heads,
                                initially_open=tuple(head for head in heads if before[head])))
        assert result.peak_simultaneous_open == max(
            sum(states[head][tick] for head in heads) for tick in range(duration)
        )
        for row in result.usage:
            vector = states[row.head_id]
            boundaries = [index for index in range(duration + 1) if not vector[index]]
            start = -1
            longest = cycles = 0
            for end in [*boundaries, duration]:
                length = end - start - 1
                longest = max(longest, length)
                if length > 0 and not vector[end] and (start >= 0 or not before[row.head_id]):
                    cycles += 1
                start = end
            assert row.observed_open_ticks == sum(vector[:duration])
            assert row.longest_open_ticks == longest
            assert row.completed_cycles == cycles
            assert (row.left_censored, row.right_censored) == (before[row.head_id], vector[-1])


@pytest.mark.parametrize("commands", [
    (ValveCommand(1, "A", True), ValveCommand(1, "A", False)),
    (ValveCommand(3, "A", True), ValveCommand(2, "A", False)),
    (ValveCommand(1, "A", True), ValveCommand(2, "A", True)),
    (ValveCommand(0, "A", False),),
    (ValveCommand(9, "A", True),),
    (ValveCommand(0, "uninstalled", True),),
])
def test_ambiguous_or_inconsistent_transitions_cannot_produce_evidence(commands):
    try:
        review(history(*commands))
    except DomainError:
        return
    raise AssertionError("An inconsistent valve history was accepted.")


def test_integer_ticks_preserve_short_window_at_large_epoch_and_strict_boundary():
    epoch = 10**20
    case = history(ValveCommand(epoch, "A", True), ValveCommand(epoch + 1, "A", False),
                   duration=epoch + 1, tick_seconds=0.01)
    usage = review(case).usage[0]
    assert (usage.observed_open_ticks, usage.longest_open_ticks, usage.completed_cycles) == (1, 1, 1)
    with pytest.raises(DomainError):
        ValveCommand(True, "A", True)


def test_command_cycles_cross_the_actual_59000_reference():
    commands = tuple(ValveCommand(i, "A", i % 2 == 0) for i in range(2 * 59_001))
    case = history(*commands, duration=2 * 59_001, heads=("A",), tick_seconds=0.001)
    at_limit = review(replace(case, commands=commands[:-2]))
    beyond = review(case)
    first = next(row for row in at_limit.comparisons if row.metric == "cycle_life")
    second = next(row for row in beyond.comparisons if row.metric == "cycle_life")
    assert (first.demand, first.published_limit.value, first.status) == (59_000, 59_000, "no_exceedance_identified")
    assert (second.demand, second.status) == (59_001, "exceeds_reference")


def test_censored_duration_can_still_prove_reference_exceedance():
    result = review(history(duration=13, heads=("A",), initially_open=("A",), tick_seconds=3600.0))
    checks = {row.metric: row for row in result.comparisons}
    assert checks["single_burn_life"].demand == 46_800
    assert checks["single_burn_life"].status == "exceeds_reference"
    assert checks["cycle_life"].demand == 0
    assert checks["cycle_life"].status == "unassessed"


def test_decimal_clock_representation_near_duration_limit_requires_review():
    for stop, expected in ((431_999, "no_exceedance_identified"), (432_000, "unassessed"),
                           (432_001, "exceeds_reference")):
        result = review(history(ValveCommand(0, "A", True), ValveCommand(stop, "A", False),
                                duration=stop, heads=("A",), tick_seconds=0.1))
        check = next(row for row in result.comparisons if row.metric == "single_burn_life")
        assert check.status == expected


def test_project_archives_history_and_blocks_a_head_exceedance_hidden_by_mean():
    request = life_study(100, installed=2, firing=1).input_data
    case = history(ValveCommand(0, "A", True), ValveCommand(13, "A", False),
                   duration=13, tick_seconds=3600.0)
    request["candidates"][0]["command_history"] = asdict(case)
    with coherent_life_catalogue():
        study = project_from_json(json.dumps(request))
    candidate = study.candidates[0]
    assert candidate.command_life_review.peak_simultaneous_open == 1
    assert candidate.catalogue_flow_check.status == "within_published_range"
    assert study.scores[0].eligibility == "incomplete"
    assert study.scores[0].score is None
    assert "Published command-life references exceeded: A: single_burn_life" in build_project_html(study)
    saved = study.to_dict()
    with coherent_life_catalogue():
        assert project_from_json(json.dumps(saved["inputs"])).to_dict() == saved


def test_project_flags_peak_concurrency_and_requires_every_installed_head():
    request = life_study(100, installed=2, firing=1).input_data
    case = history(ValveCommand(0, "A", True), ValveCommand(1, "B", True),
                   ValveCommand(2, "A", False), ValveCommand(3, "B", False))
    request["candidates"][0]["command_history"] = asdict(case)
    with coherent_life_catalogue():
        study = project_from_json(json.dumps(request))
    assert study.candidates[0].command_life_review.peak_simultaneous_open == 2
    assert study.scores[0].score is None
    assert any("simultaneous heads" in note for note in study.candidates[0].unresolved)
    request["candidates"][0]["command_history"] = asdict(history(heads=("A",)))
    with pytest.raises(DomainError, match="every installed head"):
        project_from_json(json.dumps(request))
