"""Mapping between absolute binary64 epochs and exact profile-relative times.

A profile boundary ``t_k`` is placed on the absolute clock as the canonical
date ``round(ignition + bias + t_k)``, rounded once from the exact rational
sum [1]_.  An epoch clock that maps distinct profile boundaries to the same
binary64 time cannot represent that profile and is rejected.  No minimum
relative duration is imposed.

References
----------
.. [1] PYTHON-3.12-FRACTION. Python 3.12, "fractions": float inputs retain
       their exact binary value as a rational,
       https://docs.python.org/3.12/library/fractions.html.
"""

import itertools
from bisect import bisect_left
from fractions import Fraction
from functools import lru_cache
from math import isfinite

from rocket_propulsion.core.errors import DomainError

from .numeric_outputs import _finite_rational

_TIME_POLICY = 'exact_relative_clock_and_once_rounded_epoch_events_v1'


@lru_cache(maxsize=128)
def _epoch_table(segments, ignition, bias):
    """Return the clock table of a complete profile (no cutoff horizon)."""
    return _clock_table(_profile_boundaries(segments, None), ignition, bias)


@lru_cache(maxsize=256)
def _profile_boundaries(segments, horizon):
    """Return the profile boundaries up to an optional exact cutoff ``horizon``.

    Raises
    ------
    DomainError
        If ``horizon`` lies outside ``[0, end of history]``.
    """
    if horizon is None:
        return (0.0,) + tuple(segment.end_time_s for segment in segments)
    horizon = Fraction(horizon)
    if not 0 <= horizon <= segments[-1].end_time_s:
        raise DomainError('A profile clock horizon must be inside the history.')
    if horizon == 0:
        return (0.0,)
    return (0.0,) + tuple(s.end_time_s for s in segments if s.end_time_s < horizon) + (horizon,)


@lru_cache(maxsize=128)
def _clock_table(boundaries, ignition, bias):
    """Return canonical boundary dates and their unambiguous float aliases.

    Each canonical date is the exact sum ``ignition + bias + t_k`` rounded
    once [1]_.  Canonical dates have priority.  Dates produced by either
    binary addition order, ``(ignition + bias) + t_k`` or
    ``ignition + (bias + t_k)``, are accepted as aliases only when they name
    one boundary and no other canonical date.  A loss of information in an
    alternative addition order does not invalidate an otherwise representable
    canonical clock.

    Returns
    -------
    canonical : tuple of float
        Strictly increasing canonical dates.
    images : tuple of (float, float or Fraction)
        Sorted ``(date, boundary)`` pairs, canonical dates and aliases.

    Raises
    ------
    DomainError
        If two successive boundaries share one canonical date.

    References
    ----------
    .. [1] PYTHON-3.12-FRACTION. Python 3.12, "fractions": float inputs retain
           their exact binary value as a rational,
           https://docs.python.org/3.12/library/fractions.html.
    """
    offset = Fraction(ignition) + Fraction(bias)
    canonical = tuple(_finite_rational(offset + Fraction(t), 'Absolute profile boundary')
                      for t in boundaries)
    if any(right <= left for left, right in itertools.pairwise(canonical)):
        raise DomainError('Absolute clock cannot distinguish successive profile boundaries.',
                          field='Epoch resolution')
    primary = dict(zip(canonical, boundaries, strict=True))
    proposed = {}
    for boundary in boundaries:
        if not isinstance(boundary, float):
            continue  # Internal exact cutoffs have only their canonical date.
        for date in ((ignition + bias) + boundary, ignition + (bias + boundary)):
            if isfinite(date):
                proposed.setdefault(date, set()).add(boundary)
    images = dict(primary)
    for date, owners in proposed.items():
        if date not in primary and len(owners) == 1:
            images[date] = next(iter(owners))
    return canonical, tuple(sorted(images.items()))


def _epoch_profile_time(segments, time, ignition, bias, *, horizon=None, allow_aliases=True):
    """Return the exact profile-relative time of an absolute clock reading.

    A reading equal to a canonical boundary date (or, with
    ``allow_aliases``, a supported alias) maps to that boundary exactly;
    any other reading maps to ``time - ignition - bias`` evaluated without
    rounding [1]_.

    References
    ----------
    .. [1] PYTHON-3.12-FRACTION. Python 3.12, "fractions": float inputs retain
           their exact binary value as a rational,
           https://docs.python.org/3.12/library/fractions.html.
    """
    boundaries = _profile_boundaries(segments, horizon)
    dates, images = _clock_table(boundaries, ignition, bias)
    if not allow_aliases:
        index = bisect_left(dates, time)
        if index < len(dates) and dates[index] == time:
            return Fraction(boundaries[index])
    else:
        index = bisect_left(images, time, key=lambda row: row[0])
        if index < len(images) and images[index][0] == time:
            return Fraction(images[index][1])
    return Fraction(time) - Fraction(ignition) - Fraction(bias)
