"""Hostile checks for incomplete, inconsistent and undersampled transients."""

from dataclasses import replace
from fractions import Fraction
from itertools import pairwise, product
from math import nextafter
from sys import float_info

import pytest

from rocket_propulsion.core.errors import DomainError
from rocket_propulsion.studies.sass import (
    TransientSample,
    TransientTrace,
    assumption,
    review_transient,
)


def trace(*rows, unit="Pa", period=1, slew=1, gap=8):
    return TransientTrace(
        "synthetic channel", unit, assumption(period, "s", "synthetic exact clock"),
        tuple(TransientSample(*row) for row in rows), "synthetic source", "synthetic acquisition",
        gap, None if slew is None else assumption(slew, f"{unit}/s", "synthetic physical slew"),
    )


def reduce(case, start=0, end=4):
    return review_transient(case, start, end, window_review_record="synthetic reviewed window")


def test_identical_samples_hide_the_largest_pressure_peak():
    # With 2 Pa/s and a four-second gap the reachable pressure tent peaks at 14 Pa.
    # The tent area above the 10 Pa baseline is 8 Pa s.
    result = reduce(trace((0, 10, 0), (4, 10, 0), slew=2))
    assert result.status == "bounded"
    assert result.peak == (10, 14)
    assert result.signed_integral == (32, 48)
    assert result.positive_impulse_ns is None


def test_asymmetric_cone_crossing_at_a_subsample_time():
    # 0 N to 2 N at 3 N/s over 2 s. The upper crossing is at 4/3 s and
    # the lower crossing at 2/3 s. Each side contributes a signed triangle.
    result = reduce(trace((0, 0, 0), (2, 2, 0), unit="N", slew=3), end=2)
    assert result.peak == (2, 4)
    lo, hi = result.signed_integral
    assert Fraction(lo) <= Fraction(-2, 3) <= Fraction(nextafter(lo, float("inf")))
    assert Fraction(nextafter(hi, -float("inf"))) <= Fraction(14, 3) <= Fraction(hi)


def test_reverse_thrust_cannot_erase_forward_impulse():
    # Endpoint separation saturates the slope bound, so there is only one
    # possible path. The two signed triangles cancel, the positive one remains.
    result = reduce(trace((0, -2, 0), (4, 2, 0), unit="N", slew=1))
    assert result.peak == (2, 2)
    assert result.signed_integral == (0, 0)
    assert result.positive_impulse_ns == (2, 2)


def test_negative_thrust_does_not_become_positive_impulse():
    result = reduce(trace((0, -4, 0), (4, -4, 0), unit="N", slew=1))
    assert result.signed_integral == (-20, -12)
    assert result.positive_impulse_ns == (0, 0)


def test_missing_edges_gaps_and_missing_slew_never_publish_bounds():
    cases = [
        (trace((1, 10, 0), (4, 10, 0)), "missing_window_boundary"),
        (trace((0, 10, 0), (3, 10, 0)), "missing_window_boundary"),
        (trace((-1, 10, 0), (5, 10, 0)), "missing_window_boundary"),
        (trace(), "missing_window_boundary"),
        (trace((0, 10, 0), (3, 10, 0), (4, 10, 0), gap=2), "sample_gap"),
        (trace((0, 10, 0), (4, 10, 0), slew=None), "missing_slew"),
    ]
    for case, expected in cases:
        result = reduce(case)
        assert result.status == expected
        assert result.peak is result.signed_integral is result.positive_impulse_ns is None
        assert result.trace is case and result.window_ticks == (0, 4)
    assert reduce(trace((0, 10, 0), (4, 10, 0), gap=4)).status == "bounded"


def test_pairwise_feasibility_does_not_establish_a_common_middle_value():
    # At most 1 Pa can change per second. Either adjacent pair is feasible
    # with a suitable middle value, but no middle value satisfies both.
    case = trace((0, 0, 0), (1, 50, 50), (2, 100, 0), gap=1)
    result = reduce(case, end=2)
    assert result.status == "inconsistent_samples"
    assert result.peak is result.signed_integral is result.positive_impulse_ns is None


def test_late_constraint_tightens_earlier_samples_even_with_zero_slew():
    case = trace((0, 9, 9), (1, 8, 8), (4, 3, 0), slew=0)
    result = reduce(case)
    assert result.peak == (3, 3)
    assert result.signed_integral == (12, 12)


def test_clock_uses_exact_tick_differences_and_requested_window():
    tick = 10**24
    case = trace(
        (tick - 10, 800, 0), (tick, 0, 0), (tick + 1, 0, 0),
        (tick + 3, 0, 0), (tick + 20, 900, 0),
        period=0.25, slew=8, gap=2,
    )
    result = reduce(case, tick, tick + 3)
    assert result.peak == (0, 2)
    assert result.signed_integral == (-0.625, 0.625)


def test_against_enumerated_piecewise_linear_truths():
    # Independent oracle enumerates every integer-valued truth at t=0..4
    # within this grid, interpolates linearly, and sums the analytic segment
    # areas. The observation keeps t=0, 2, 4 and broadens one measurement.
    accepted = 0
    for actual in product(range(-2, 3), repeat=5):
        if any(abs(right - left) > 1 for left, right in pairwise(actual)):
            continue
        accepted += 1
        result = reduce(trace(
            (0, actual[0], 0), (2, actual[2] + 1, 1), (4, actual[4], 0), unit="N"
        ))
        assert result.status == "bounded"
        assert result.peak[0] <= max(actual) <= result.peak[1], actual
        area = Fraction(sum(left + right for left, right in pairwise(actual)), 2)
        positive_area = Fraction(sum(max(left, 0) + max(right, 0)
                                     for left, right in pairwise(actual)), 2)
        assert result.signed_integral[0] <= area <= result.signed_integral[1], actual
        assert result.positive_impulse_ns[0] <= positive_area <= result.positive_impulse_ns[1], actual
    assert accepted == 259


def test_outward_rounding_keeps_nonzero_subnormal_area_enclosed():
    smallest = nextafter(0.0, 1.0)
    case = trace((0, smallest, 0), (1, smallest, 0), unit="N", slew=0, period=0.5)
    result = reduce(case, end=1)
    assert result.signed_integral == result.positive_impulse_ns == (0, smallest)


def test_nonrepresentable_interval_fails_closed():
    try:
        reduce(trace((0, float_info.max, 0), (4, float_info.max, 0), slew=0))
    except DomainError:
        return
    raise AssertionError("overflowed physical interval was published")


def test_units_error_and_sampling_order_fail_closed():
    baseline = trace((0, 0, 0), (4, 0, 0))
    for kwargs in (
        {"unit": "bar"},
        {"tick_period": assumption(1, "ms", "wrong units")},
        {"max_absolute_slew": assumption(1, "N/s", "wrong channel")},
        {"max_absolute_slew": assumption(-1, "Pa/s", "negative slew")},
        {"max_gap_ticks": True},
        {"acquisition_review_record": " "},
        {"samples": (TransientSample(0, 0, 0), TransientSample(0, 1, 0))},
        {"samples": tuple(reversed(baseline.samples))},
    ):
        with pytest.raises(DomainError):
            replace(baseline, **kwargs)
    for row in ((True, 0, 0), (0, float("nan"), 0), (0, 0, -1), (0, 0, float("inf"))):
        with pytest.raises(DomainError):
            TransientSample(*row)
    for start, end, review in ((0, 0, "review"), (1, 0, "review"), (0.0, 4, "review"), (0, 4, "")):
        with pytest.raises(DomainError):
            review_transient(baseline, start, end, window_review_record=review)

