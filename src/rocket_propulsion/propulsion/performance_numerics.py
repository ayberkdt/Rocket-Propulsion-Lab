"""Checked positive log/exp expressions with directed decimal enclosures.

The finite precision budget can fail unresolved. It is not a proof of physical
model validity or universal successful inversion for every accepted input.
"""
from decimal import Decimal
from fractions import Fraction
from functools import lru_cache
from math import isfinite
from sys import float_info

from rocket_propulsion.core.errors import ConvergenceError, DomainError

from .burns.marginal_numerics import _PRECISIONS, _Intervals, _log_interval

_PERFORMANCE_POLICY = 'exact_binary_rational_ledgers_and_enclosed_log_exp_v1'
_OVERFLOW_BOUNDARY = Fraction(float_info.max) + 2**970
_ZERO_BOUNDARY = Fraction(1, 2**1075)


def _positive_enclosed_output(evaluate, label):
    """Export only after both exact decimal bounds round to the same float."""
    for precision in _PRECISIONS:
        lower, upper = evaluate(precision)
        if not lower.is_finite() or not upper.is_finite() or lower > upper or upper <= 0:
            raise ConvergenceError('Invalid positive expression enclosure.', field=label)
        lower, upper = max(Fraction(0), Fraction(lower)), Fraction(upper)
        if lower >= _OVERFLOW_BOUNDARY or upper <= _ZERO_BOUNDARY:
            raise DomainError(f'{label} is outside the nonzero finite output range.', field=label)
        try:
            left, right = float(lower), float(upper)
        except OverflowError:
            continue
        if left == right and left > 0 and isfinite(left):
            return left
    raise ConvergenceError('Performance expression rounding is unresolved within the precision budget.',
                           field=label)


def _scaled_log(scale, ratio, label):
    def evaluate(precision):
        math = _Intervals(precision)
        return math.mul(math.rational(scale), _log_interval(ratio, precision))
    return _positive_enclosed_output(evaluate, label)


@lru_cache(maxsize=2048)
def _exp_interval(value, precision):
    if value == 0:
        return Decimal(1), Decimal(1)
    math = _Intervals(precision)
    lower, upper = math.rational(value)
    result = (math.near.next_minus(math.near.exp(lower)),
              math.near.next_plus(math.near.exp(upper)))
    # Monotonicity gives these exact bounds for nonzero signed exponents.
    return ((max(Decimal(1), result[0]), result[1]) if value > 0
            else (max(Decimal(0), result[0]), min(Decimal(1), result[1])))


def _scaled_exp(scale, exponent, label):
    def evaluate(precision):
        math = _Intervals(precision)
        return math.mul(math.rational(scale), _exp_interval(exponent, precision))
    return _positive_enclosed_output(evaluate, label)


def _scaled_exp_complement(scale, exponent, label):
    """Export scale*(1-exp(-exponent)) without exporting survival or its complement."""
    def evaluate(precision):
        math = _Intervals(precision)
        complement = math.sub((Decimal(1), Decimal(1)), _exp_interval(-exponent, precision))
        return math.mul(math.rational(scale), complement)
    return _positive_enclosed_output(evaluate, label)


def _inverse_mass_outputs(initial, exponent):
    # e > 2, so exp(x) > 2**1024 for x > 1024, beyond the binary64
    # overflow midpoint. This is an output impossibility proof, not a clamp.
    if exponent > 1024:
        raise DomainError('Mass ratio is outside the finite output range.', field='Mass ratio')
    ratio = _scaled_exp(Fraction(1), exponent, 'Mass ratio')
    final = _scaled_exp(initial, -exponent, 'Final mass')
    consumed = _scaled_exp_complement(initial, exponent, 'Propellant mass')
    fraction = _scaled_exp_complement(Fraction(1), exponent, 'Propellant mass fraction')
    return ratio, final, consumed, fraction
