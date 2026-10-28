"""Exact native mass ledgers and guarded inversion on the binary64 clock.

Consumed propellant is integrated exactly, as rational arithmetic over the
supplied binary64 profile values [1]_, and depletion cutoffs are chosen as
the latest represented clock date whose exact consumption does not exceed the
target [2]_. These are numerical contracts, not physical uncertainty or
hardware acceptance rules. Callers fold their inventory margin into the
target mass explicitly.

References
----------
.. [1] PYTHON-3.12-FRACTION. Python 3.12, "fractions - Rational numbers",
       https://docs.python.org/3.12/library/fractions.html
.. [2] PYTHON-3.12-NEXTAFTER. Python 3.12, ``math.nextafter``,
       https://docs.python.org/3.12/library/math.html#math.nextafter
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
    """Return a nonnegative mass-flow scale driver as an exact rational."""
    result = _fraction(value, 'Native mass-flow scale')
    if result < 0:
        raise DomainError('Native mass-flow scale must be nonnegative.')
    return result


def _consumption(segments, scale, time):
    """Return the exact propellant mass consumed from ignition to profile time ``time``.

    Each segment's total mass flow varies linearly between its endpoint
    values. For a segment of duration ``h`` reached up to the fraction
    ``x = min((time - start) / h, 1)`` the exact contribution is
    ``h (m0 x (1 - x/2) + m1 x^2 / 2)``. Contributions are summed as exact
    rationals [1]_ and multiplied by ``scale``; nothing is rounded.

    References
    ----------
    .. [1] PYTHON-3.12-FRACTION. Python 3.12, "fractions - Rational numbers",
           https://docs.python.org/3.12/library/fractions.html
    """
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
    """Return an unsigned integer key whose order matches binary64 float order.

    The map is strictly monotone over finite floats, including negative
    values, and orders ``-0.0`` immediately before ``+0.0``.
    """
    bits = unpack('>Q', pack('>d', value))[0]
    return (~bits & ((1 << 64)-1)) if bits >> 63 else bits | (1 << 63)


def _from_ordered(key):
    """Invert :func:`_ordered`."""
    bits = key & ((1 << 63)-1) if key >> 63 else ~key & ((1 << 64)-1)
    return unpack('>d', pack('>Q', bits))[0]


def _greatest_safe(lower, upper, safe):
    """Return the greatest float in ``[lower, upper]`` for which ``safe`` holds.

    ``safe`` must be monotone (true up to some date, false afterwards) and
    true at ``lower``. Bisection runs over the ordered binary64 lattice
    rather than over real numbers, so the result is an adjacent-float
    boundary [1]_ and at most 65 predicate calls are made.

    References
    ----------
    .. [1] PYTHON-3.12-NEXTAFTER. Python 3.12, ``math.nextafter``,
           https://docs.python.org/3.12/library/math.html#math.nextafter
    """
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
    """Round a nonnegative exact deficit up to the next representable float.

    The returned float is never below ``value``; subnormal results are
    supported. A negative or unrepresentably large value raises
    :class:`DomainError`.

    References
    ----------
    .. [1] PYTHON-3.12-NEXTAFTER. Python 3.12, ``math.nextafter``,
           https://docs.python.org/3.12/library/math.html#math.nextafter
    """
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
    """Represented cutoff: clock date, exact profile time and exact mass shortfall."""

    date: float
    profile_time: Fraction
    shortfall: Fraction


def _cutoff(segments, scale, target, offset):
    """Return the represented cutoff at which consumption first reaches ``target``.

    Parameters
    ----------
    segments : sequence
        Profile segments with start/end times and endpoint total mass flows.
    scale : Fraction
        Mass-flow scale driver.
    target : Fraction
        Mass that may be consumed, already reduced by the caller's margin.
    offset : Fraction
        Exact ignition epoch on the consumer clock.

    Returns
    -------
    _Cut or None
        ``None`` when the whole profile consumes less than ``target``;
        otherwise the greatest represented date whose exact consumption does
        not exceed ``target``, its exact profile time and the remaining
        shortfall. A nonpositive target cuts at ignition.

    Raises
    ------
    DomainError
        If the depletion date lies beyond the finite clock range, or the
        cut would be a zero-duration interval although a positive one is
        required.

    Notes
    -----
    The search is confined to the first segment in which cumulative
    consumption reaches the target. When several profile boundaries map to
    the same clock date, that date names the earliest of them. The consumer
    clock admits only the selected prefix.
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
    primary: dict[float, Fraction] = {}
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
