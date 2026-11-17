"""Adversarial closure anchoring and bounded impulse on uncertain time windows."""

from dataclasses import replace
from fractions import Fraction
from itertools import product

import pytest

from rocket_propulsion.core.errors import DomainError
from rocket_propulsion.studies.sass import ValveCommand, review_transient
from rocket_propulsion.studies.sass.shutdown_observation import (
    ShutdownWindow,
    review_shutdown_impulse,
)
from tests.test_sass_command_life import history
from tests.test_sass_transient_observation import trace
from tests.test_sass_valve_response import limits, observed


def review(commanded, measured, force, *selection, end=15, sync=1, error=1, verified=None,
           approved=None):
    return review_shutdown_impulse(
        commanded, observed(measured, error=error, verified=verified),
        limits(5, 5) if approved is None else approved, force, selection,
        ShutdownWindow(end, sync, "synthetic force clock alignment", "synthetic endpoint review"),
    )


def histories(*, duration=18):
    close_a, close_b = ValveCommand(2, "A", False), ValveCommand(5, "B", False)
    commanded = history(close_a, close_b, initially_open=("A", "B"), duration=duration)
    measured = history(ValveCommand(4, "A", False), ValveCommand(8, "B", False),
                       initially_open=("A", "B"), duration=duration)
    return commanded, measured, (close_a, close_b)


def constant_force():
    return trace((3, 2, 0), (7, 2, 0), (11, 2, 0), (15, 2, 0), unit="N", gap=4, slew=0)


def test_last_physical_closure_and_both_clock_errors_anchor_impulse():
    commanded, measured, selection = histories()
    result = review(commanded, measured, constant_force(), *selection)
    assert result.status == "bounded"
    assert result.final_closure_ticks == (6, 10)
    assert [row.observation_tick for row in result.closures] == [4, 8]
    assert result.thrust_review.window_ticks == (3, 15)
    assert result.thrust_review.positive_impulse_ns == (10, 18)
    assert result.thrust_review.signed_integral == (24, 24)


def test_analytic_linear_path_for_all_uncertain_closure_locations():
    base = 10**24
    # F(t) = t - 4 N and t is seconds since the first sample. The slope
    # saturates its physical bound, so the force path is uniquely determined.
    force = trace((base, -4, 0), (base + 8, 4, 0), unit="N", slew=1)
    for earliest, latest in product(range(9), repeat=2):
        if earliest > latest:
            continue
        result = review_transient(
            force, base, base + 8, window_review_record="synthetic",
            positive_impulse_start_ticks=(base + earliest, base + latest),
        )
        lower = Fraction(16 - max(latest - 4, 0)**2, 2)
        upper = Fraction(16 - max(earliest - 4, 0)**2, 2)
        assert result.positive_impulse_ns == (lower, upper), (earliest, latest)
        assert result.signed_integral == (0, 0)


def test_subsample_closure_does_not_erase_between_sample_force():
    force = trace((0, 0, 0), (4, 0, 0), unit="N", slew=2)
    result = review_transient(force, 0, 4, window_review_record="synthetic",
                              positive_impulse_start_ticks=(1, 3))
    # Upper force tent peaks at 4 N at t=2. Clipping t<1 removes just 1 N*s.
    assert result.positive_impulse_ns == (0, 7)
    assert result.positive_impulse_start_ticks == (1, 3)


def test_no_command_based_fallback_when_physical_closure_is_unresolved():
    commanded, measured, selection = histories()
    for uncertain in (replace(measured, commands=measured.commands[:1]),
                      replace(measured, commands=(ValveCommand(2, "A", False), measured.commands[1]))):
        result = review(commanded, uncertain, constant_force(), *selection)
        assert result.status == "closure_unresolved"
        assert result.final_closure_ticks is result.thrust_review is None
    uncovered = review(commanded, measured, constant_force(), *selection, verified=("A",))
    assert uncovered.status == "closure_unresolved"


def test_late_but_identified_closure_still_bounds_tail_and_retains_timing_failure():
    commanded, measured, selection = histories()
    result = review(commanded, measured, constant_force(), *selection, approved=limits(5, 1))
    assert result.status == "bounded"
    assert result.valve_review.status == "fail"
    assert [row.status for row in result.closures] == ["timing_indeterminate", "late_response"]
    assert result.thrust_review.positive_impulse_ns == (10, 18)


def test_reopening_before_the_last_closure_invalidates_shared_shutdown():
    commanded, measured, selection = histories()
    reopens = (measured.commands[0], ValveCommand(6, "A", True), measured.commands[1])
    altered = replace(measured, commands=reopens)
    result = review(commanded, altered, constant_force(), *selection)
    assert result.status == "possible_valve_reopening"
    assert result.thrust_review is None


def test_reopening_within_clock_error_at_tail_end_is_not_outside():
    commanded, measured, selection = histories(duration=20)
    for tick, expected in ((17, "possible_valve_reopening"), (18, "bounded")):
        with_open = replace(measured, commands=(*measured.commands, ValveCommand(tick, "B", True)))
        result = review(commanded, with_open, constant_force(), *selection)
        assert result.status == expected


def test_issued_reopening_at_tail_boundary_is_a_new_condition():
    commanded, measured, selection = histories(duration=20)
    for tick, expected in ((16, "subsequent_opening_command"), (17, "bounded")):
        with_open = replace(commanded, commands=(*commanded.commands, ValveCommand(tick, "B", True)))
        result = review(with_open, measured, constant_force(), *selection)
        assert result.status == expected


def test_valve_and_force_records_both_cover_the_tail():
    commanded, measured, selection = histories(duration=16)
    force = constant_force()
    result = review(commanded, measured, force, *selection)
    assert result.status == "tail_window_not_covered"
    commanded, measured, selection = histories()
    for samples in (force.samples[1:], force.samples[:-1]):
        result = review(commanded, measured, replace(force, samples=samples), *selection)
        assert result.status == "missing_window_boundary"
        assert result.thrust_review is None or result.thrust_review.positive_impulse_ns is None
    result = review(commanded, measured, replace(force, max_gap_ticks=3), *selection)
    assert result.status == "sample_gap"
    assert result.thrust_review.positive_impulse_ns is None
    result = review(commanded, measured, force, *selection, end=10)
    assert result.status == "tail_window_not_covered"


def test_invalid_units_clocks_and_selection_fail_closed():
    commanded, measured, selection = histories()
    force = constant_force()
    for invalid in ((), (selection[0], selection[0]), (ValveCommand(2, "A", True),),
                    (ValveCommand(3, "A", False),)):
        with pytest.raises(DomainError):
            review(commanded, measured, force, *invalid)
    for invalid_force in (trace((0, 1, 0), (15, 1, 0)),
                          trace((0, 1, 0), (15, 1, 0), unit="N", period=.5)):
        with pytest.raises(DomainError):
            review(commanded, measured, invalid_force, *selection)
    for malformed in ((5, 4), (0, True), (1, 5), (-1, 4), (0,), [0, 4]):
        with pytest.raises(DomainError):
            review_transient(force, 3, 7, window_review_record="synthetic",
                             positive_impulse_start_ticks=malformed)

