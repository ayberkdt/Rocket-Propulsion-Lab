"""Binary clock coordinates; canonical events are rounded once from exact sums.

An epoch clock that maps distinct profile boundaries to the same binary64 time
cannot represent that profile. No minimum relative duration is imposed.
"""

from bisect import bisect_left
from fractions import Fraction
from functools import lru_cache
from math import isfinite

from rocket_propulsion.core.errors import DomainError
from .numeric_outputs import _finite_rational

_TIME_POLICY = 'exact_relative_clock_and_once_rounded_epoch_events_v1'


@lru_cache(maxsize=128)
def _epoch_table(segments, ignition, bias):
    """Return canonical dates and unambiguous supported formation aliases.

    Canonical dates have priority. Aliases from either binary addition order
    are accepted only when they name one boundary and no other canonical date.
    Losing information in an alternative addition order cannot invalidate an
    otherwise representable canonical clock.
    """
    boundaries = (0.0,) + tuple(segment.end_time_s for segment in segments)
    offset = Fraction(ignition) + Fraction(bias)
    canonical = tuple(_finite_rational(offset + Fraction(t), 'Absolute profile boundary')
                      for t in boundaries)
    if any(right <= left for left, right in zip(canonical, canonical[1:])):
        raise DomainError('Absolute clock cannot distinguish successive profile boundaries.',
                          field='Epoch resolution')
    primary = dict(zip(canonical, boundaries, strict=True))
    proposed = {}
    for boundary in boundaries:
        for date in ((ignition + bias) + boundary, ignition + (bias + boundary)):
            if isfinite(date):
                proposed.setdefault(date, set()).add(boundary)
    images = dict(primary)
    for date, owners in proposed.items():
        if date not in primary and len(owners) == 1:
            images[date] = next(iter(owners))
    return canonical, tuple(sorted(images.items()))


def _epoch_profile_time(segments, time, ignition, bias):
    """Preserve exact cancellation; map only an explicitly represented event."""
    _, images = _epoch_table(segments, ignition, bias)
    index = bisect_left(images, time, key=lambda row: row[0])
    if index < len(images) and images[index][0] == time:
        return Fraction(images[index][1])
    return Fraction(time) - Fraction(ignition) - Fraction(bias)
