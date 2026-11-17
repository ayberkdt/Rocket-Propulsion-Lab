"""Bound measured transients between discrete samples on a reviewed time base."""

from __future__ import annotations

from bisect import bisect_left, bisect_right
from dataclasses import dataclass
from fractions import Fraction
from itertools import pairwise
from typing import Literal

from rocket_propulsion.core.errors import DomainError
from rocket_propulsion.core.validation import integer, real, records, text

from .evidence import CitedScalar, quantity
from .rational_numerics import outward_endpoint


@dataclass(frozen=True, slots=True)
class TransientSample:
    """Calibrated sample with an absolute, simultaneous error bound.

    References
    ----------
    .. [1] NASA-7009B, 4.2.5 and 4.2.7, referents and uncertainty.
    """

    tick: int
    value: float
    absolute_error: float

    def __post_init__(self) -> None:
        """Normalize and validate one sample before combining bounds.

        References
        ----------
        .. [1] NASA-7009B, 4.2.7, uncertainty characterization.
        """
        tick = integer(self.tick, "transient sample tick")
        error = real(self.absolute_error, "transient sample absolute error")
        if error < 0:
            raise DomainError("Transient sample absolute error must be nonnegative.")
        object.__setattr__(self, "tick", tick)
        object.__setattr__(self, "value", real(self.value, "transient sample value"))
        object.__setattr__(self, "absolute_error", error)


@dataclass(frozen=True, slots=True)
class TransientTrace:
    """One physical measurement channel on a reviewed common clock.

    The acquisition record must cover sensor position and orientation, the
    usable bandwidth and sample quality, calibration, and any time uncertainty
    accounted for in the value bounds. Clock period is exact for this reduction.
    A slew bound applies to the physical signal throughout the analysis window.

    References
    ----------
    .. [1] NASA-7009B, 4.2.5 and 4.2.7, referents and uncertainty.
    """

    channel_id: str
    unit: Literal["Pa", "N"]
    tick_period: CitedScalar
    samples: tuple[TransientSample, ...]
    source_record: str
    acquisition_review_record: str
    max_gap_ticks: int
    max_absolute_slew: CitedScalar | None = None

    def __post_init__(self) -> None:
        """Validate the channel, time base, evidence identifiers and sample order.

        References
        ----------
        .. [1] NASA-7009B, 4.2.5, documented empirical referents.
        """
        if self.unit not in ("Pa", "N"):
            raise DomainError("Transient channel unit must be Pa or N.")
        quantity(self.tick_period, "s", strictly_positive=True)
        if self.max_absolute_slew is not None:
            quantity(self.max_absolute_slew, f"{self.unit}/s", minimum=0.0)
        max_gap = integer(self.max_gap_ticks, "maximum sample gap ticks")
        if max_gap <= 0:
            raise DomainError("Maximum sample gap must be positive.")
        samples: tuple[TransientSample, ...] = records(
            self.samples, TransientSample, "transient samples"
        )
        if any(left.tick >= right.tick for left, right in pairwise(samples)):
            raise DomainError("Transient sample ticks must be strictly increasing.")
        object.__setattr__(self, "samples", samples)
        object.__setattr__(self, "max_gap_ticks", max_gap)
        for field in ("channel_id", "source_record", "acquisition_review_record"):
            object.__setattr__(self, field, text(getattr(self, field), field))


@dataclass(frozen=True, slots=True)
class TransientReview:
    """Conservative bounds or a specific reason they could not be established.

    Peak uses the channel unit. Integrals use channel unit times seconds.
    For thrust, positive impulse is the integral of max(thrust, 0) along the
    declared measurement direction. Signed impulse can cancel across reversals.

    References
    ----------
    .. [1] NASA-7009B, 4.3.8, scope and limitations of results.
    """

    status: Literal[
        "bounded", "missing_window_boundary", "sample_gap", "missing_slew", "inconsistent_samples"
    ]
    trace: TransientTrace
    window_ticks: tuple[int, int]
    window_review_record: str
    peak: tuple[float, float] | None = None
    signed_integral: tuple[float, float] | None = None
    positive_impulse_ns: tuple[float, float] | None = None
    positive_impulse_start_ticks: tuple[int, int] | None = None


def _tighten(
    ranges: list[tuple[Fraction, Fraction]], steps: list[Fraction], slew: Fraction
) -> bool:
    """Propagate each sample constraint in both directions along the time chain.

    References
    ----------
    .. [1] NASA-7009B, 4.2.7, propagate bounded measurement uncertainty.
    """
    for indices in (range(1, len(ranges)), range(len(ranges) - 2, -1, -1)):
        for i in indices:
            j = i - indices.step
            radius = slew * steps[min(i, j)]
            lo, hi = ranges[i]
            other_lo, other_hi = ranges[j]
            lo, hi = max(lo, other_lo - radius), min(hi, other_hi + radius)
            if lo > hi:
                return False
            ranges[i] = lo, hi
    return True


def _segment_envelope(
    left: Fraction, right: Fraction, elapsed: Fraction, slew: Fraction, *, upper: bool
) -> tuple[tuple[Fraction, Fraction], ...]:
    """Return the vertices of the min of upper cones or max of lower cones.

    References
    ----------
    .. [1] NASA-7009B, 4.2.7, explicitly bounded interpolation uncertainty.
    """
    if slew == 0:
        return ((Fraction(0), left), (elapsed, right))
    gradient = slew if upper else -slew
    crossing = (right - left + gradient * elapsed) / (2 * gradient)
    return ((Fraction(0), left), (crossing, left + gradient * crossing), (elapsed, right))


def _area(
    vertices: tuple[tuple[Fraction, Fraction], ...], *, positive_only: bool,
    start_at: Fraction | None = None,
) -> Fraction:
    """Integrate each linear span exactly, splitting at a zero crossing if needed.

    References
    ----------
    .. [1] NASA-7009B, 4.2.4, numerical verification.
    """
    total = Fraction(0)
    for (left, y0), (right, y1) in pairwise(vertices):
        if start_at is not None:
            if right <= start_at:
                continue
            if left < start_at:
                y0 += (y1 - y0) * (start_at - left) / (right - left)
                left = start_at
        width = right - left
        if not positive_only or min(y0, y1) >= 0:
            total += width * (y0 + y1) / 2
        elif max(y0, y1) > 0:
            height = max(y0, y1)
            total += width * height**2 / (2 * abs(y1 - y0))
    return total


def _interval(lo: Fraction, hi: Fraction) -> tuple[float, float]:
    """Round the two published endpoints away from the exact interval.

    References
    ----------
    .. [1] NASA-7009B, 4.2.4, numerical representation error.
    """
    return outward_endpoint(lo, upward=False), outward_endpoint(hi, upward=True)


def review_transient(
    trace: TransientTrace, start_tick: int, end_tick: int, *, window_review_record: str,
    positive_impulse_start_ticks: tuple[int, int] | None = None,
) -> TransientReview:
    """Reduce a fully sampled reviewed window using a physical slew bound.

    For closure- or shutdown-relative measurements, the window reviewer must
    cover the event timing uncertainty and establish where the signal tail ends.
    The review uses only the requested window and requires both boundary samples.
    A bounded result represents channel data and its declared uncertainty.
    Optionally bound positive impulse from any start between the two supplied
    ticks to the common end. The peak and signed integral still use the whole
    sampled window. Start ticks may lie between samples.

    References
    ----------
    .. [1] NASA-7009B, 4.2.5, 4.2.7 and 4.3.8, scope and uncertainty.
    """
    if not isinstance(trace, TransientTrace):
        raise DomainError("Transient reduction requires a reviewed measurement trace.")
    start = integer(start_tick, "transient window start")
    end = integer(end_tick, "transient window end")
    if end <= start:
        raise DomainError("Transient window end must be after its start.")
    review_record = text(window_review_record, "transient window review")
    window = (start, end)
    starts = None
    if positive_impulse_start_ticks is not None:
        if (trace.unit != "N" or not isinstance(positive_impulse_start_ticks, tuple)
                or len(positive_impulse_start_ticks) != 2):
            raise DomainError("An uncertain positive impulse start requires two thrust clock ticks.")
        earliest = integer(positive_impulse_start_ticks[0], "earliest impulse start")
        latest = integer(positive_impulse_start_ticks[1], "latest impulse start")
        if not start <= earliest <= latest <= end:
            raise DomainError("Both positive impulse start ticks must be inside the sampled window.")
        starts = (earliest, latest)
    ticks = [sample.tick for sample in trace.samples]
    first, stop = bisect_left(ticks, start), bisect_right(ticks, end)
    if first == stop or ticks[first] != start or ticks[stop - 1] != end:
        return TransientReview("missing_window_boundary", trace, window, review_record,
                               positive_impulse_start_ticks=starts)
    gaps = [right - left for left, right in pairwise(ticks[first:stop])]
    if max(gaps) > trace.max_gap_ticks:
        return TransientReview("sample_gap", trace, window, review_record,
                               positive_impulse_start_ticks=starts)
    if trace.max_absolute_slew is None:
        return TransientReview("missing_slew", trace, window, review_record,
                               positive_impulse_start_ticks=starts)

    period = Fraction(trace.tick_period.value)
    slew = Fraction(trace.max_absolute_slew.value)
    steps = [gap * period for gap in gaps]
    values = [(Fraction(x.value), Fraction(x.absolute_error)) for x in trace.samples[first:stop]]
    ranges = [(value - error, value + error) for value, error in values]
    if not _tighten(ranges, steps, slew):
        return TransientReview("inconsistent_samples", trace, window, review_record,
                               positive_impulse_start_ticks=starts)

    peak_lo, peak_hi = ranges[0]
    signed_lo = signed_hi = positive_lo = positive_hi = Fraction(0)
    for offset, (((lo0, hi0), (lo1, hi1)), elapsed) in enumerate(zip(pairwise(ranges), steps)):
        lower = _segment_envelope(lo0, lo1, elapsed, slew, upper=False)
        upper = _segment_envelope(hi0, hi1, elapsed, slew, upper=True)
        peak_lo = max(peak_lo, lo1)
        peak_hi = max(peak_hi, *(value for _, value in upper))
        signed_lo += _area(lower, positive_only=False)
        signed_hi += _area(upper, positive_only=False)
        if trace.unit == "N":
            earliest_offset = latest_offset = None
            if starts is not None:
                left_tick = ticks[first + offset]
                earliest_offset = (starts[0] - left_tick) * period
                latest_offset = (starts[1] - left_tick) * period
            positive_lo += _area(lower, positive_only=True, start_at=latest_offset)
            positive_hi += _area(upper, positive_only=True, start_at=earliest_offset)

    return TransientReview(
        "bounded", trace, window, review_record, _interval(peak_lo, peak_hi),
        _interval(signed_lo, signed_hi),
        _interval(positive_lo, positive_hi) if trace.unit == "N" else None, starts,
    )

