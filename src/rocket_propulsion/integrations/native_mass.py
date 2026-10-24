"""Exact supplied-binary native mass ledgers and guarded clock-grid inversion.

These are numerical contracts, not physical uncertainty or hardware acceptance
rules. The existing caller's inventory margin is supplied explicitly.
"""

from dataclasses import dataclass
from fractions import Fraction
from math import inf, isfinite, nextafter
from struct import pack, unpack
from sys import float_info

from rocket_propulsion.core.errors import DomainError
from rocket_propulsion.propulsion.burns.numeric_outputs import _finite_rational, _fraction

_NATIVE_MASS_POLICY = 'exact_mass_ledger_guarded_epoch_cutoff_v1'


def _scale(value):
    result = _fraction(value, 'Native mass-flow scale')
    if result < 0:
        raise DomainError('Native mass-flow scale must be nonnegative.')
    return result


def _consumption(segments, scale, time):
    """Integrate positive endpoint weights before any scalar conversion."""
    total = Fraction(0)
    for segment in segments:
        start, end = _fraction(segment.start_time_s), _fraction(segment.end_time_s)
        if time <= start:
            break
        h = end-start
        x = min((time-start)/h, Fraction(1))
        first = _fraction(segment.start_total_mass_flow_kg_s)
        last = _fraction(segment.end_total_mass_flow_kg_s)
        total += h * (first*x*(1-x/2) + last*x*x/2)
    return total*scale


def _ordered(value):
    """Monotone IEEE binary64 key, including negative dates and both zeros."""
    bits = unpack('>Q', pack('>d', value))[0]
    return (~bits & ((1 << 64)-1)) if bits >> 63 else bits | (1 << 63)


def _from_ordered(key):
    bits = key & ((1 << 63)-1) if key >> 63 else ~key & ((1 << 64)-1)
    return unpack('>d', pack('>Q', bits))[0]


def _greatest_safe(lower, upper, safe):
    """Use at most 65 predicate calls across the finite binary64 lattice."""
    if safe(upper):
        return upper
    lo, hi = _ordered(lower), _ordered(upper)
    while hi-lo > 1:
        middle = (lo+hi)//2
        if safe(_from_ordered(middle)):
            lo = middle
        else:
            hi = middle
    return _from_ordered(lo)


def _upper_bound(value):
    """Round a nonnegative numerical deficit outwards, including subnormals."""
    if value < 0:
        raise DomainError('A native mass deficit bound cannot be negative.')
    try:
        result = float(value)
    except OverflowError as error:
        raise DomainError('Native mass deficit bound is outside finite range.') from error
    if not isfinite(result):
        raise DomainError('Native mass deficit bound is outside finite range.')
    if _fraction(result) < value:
        result = nextafter(result, inf)
    if not isfinite(result):
        raise DomainError('Native mass deficit bound is outside finite range.')
    return result


@dataclass(frozen=True)
class _Cut:
    date: float
    profile_time: Fraction
    shortfall: Fraction


def _cutoff(segments, scale, target, offset):
    """Find a safe represented cutoff at the first binding flow interval.

    Canonical profile dates name their earliest boundary if a later, unused
    boundary projects onto the same date. Only the selected prefix is finally
    admitted by the consumer clock. An actual zero-duration positive cut fails.
    """
    ignition = _finite_rational(offset, 'Native ignition time')
    if target <= 0:
        return _Cut(ignition, Fraction(0), Fraction(0))
    prefix = Fraction(0)
    binding = None
    for index, segment in enumerate(segments):
        h = _fraction(segment.end_time_s)-_fraction(segment.start_time_s)
        available = h*(_fraction(segment.start_total_mass_flow_kg_s)
                       + _fraction(segment.end_total_mass_flow_kg_s))*scale/2
        if prefix+available >= target:
            binding = index
            break
        prefix += available
    if binding is None:
        return None
    end = _fraction(segments[binding].end_time_s)
    primary = {}
    for boundary in (0.0,) + tuple(s.end_time_s for s in segments[:binding+1]):
        try:
            date = float(offset+_fraction(boundary))
        except OverflowError:
            continue
        if isfinite(date):
            primary.setdefault(date, _fraction(boundary))
    exact_end = offset+end
    if exact_end > _fraction(float_info.max):
        high = float_info.max
        if _consumption(segments, scale, _fraction(high)-offset) < target:
            raise DomainError('The native depletion date lies beyond the finite clock range.')
    else:
        high = _finite_rational(exact_end, 'Native depletion bracket')

    def coordinate(date):
        return primary.get(date, min(max(_fraction(date)-offset, Fraction(0)), end))

    def safe(date):
        return _consumption(segments, scale, coordinate(date)) <= target

    date = _greatest_safe(ignition, high, safe)
    point = coordinate(date)
    if date <= ignition or point <= 0:
        raise DomainError('A positive native depletion interval is not represented by this clock.',
                          field='Depletion clock resolution')
    return _Cut(date, point, target-_consumption(segments, scale, point))
