"""Positive logarithm/exponential expressions evaluated with directed decimal enclosures.

Each expression is enclosed between lower and upper decimal bounds at increasing
working precision and exported only when both bounds round to the same positive
finite binary64 value. A finite precision budget can end unresolved; that outcome
is reported as :class:`ConvergenceError`. Successful rounding is a numerical
statement only, not a proof of physical model validity or of successful inversion
for every accepted input.

References
----------
.. [1] PYTHON-3.12-DECIMAL. Python 3.12 documentation, decimal module: precision
       contexts and exact conversion of binary floats.
.. [2] PYTHON-3.12-FRACTION. Python 3.12 documentation, fractions module: a float
       converts to the rational equal to its exact binary value.
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
    """Round an enclosed positive expression to binary64 once both bounds agree.

    Parameters
    ----------
    evaluate : callable
        ``evaluate(precision) -> (lower, upper)`` returning finite ``Decimal`` bounds
        of the exact value at the given working precision.
    label : str
        Output label used in errors.

    Returns
    -------
    float
        The unique positive finite binary64 value both bounds round to.

    Raises
    ------
    ConvergenceError
        If an enclosure is invalid or the bounds do not agree within the precision
        budget.
    DomainError
        If the enclosed value lies outside the nonzero finite binary64 range.

    References
    ----------
    .. [1] PYTHON-3.12-DECIMAL. Python 3.12 documentation, decimal module: precision
           contexts and exact conversion of binary floats.
    .. [2] PYTHON-3.12-FRACTION. Python 3.12 documentation, fractions module: a float
           converts to the rational equal to its exact binary value.
    """
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
    """Export ``scale * ln(ratio)`` for a positive result.

    Parameters
    ----------
    scale : Fraction
        Exact positive multiplier.
    ratio : Fraction
        Exact argument of the logarithm, greater than one.
    label : str
        Output label used in errors.

    Returns
    -------
    float
        The correctly enclosed binary64 value.

    References
    ----------
    .. [1] PYTHON-3.12-DECIMAL. Python 3.12 documentation, decimal module: precision
           contexts and exact conversion of binary floats.
    """
    def evaluate(precision):
        math = _Intervals(precision)
        return math.mul(math.rational(scale), _log_interval(ratio, precision))
    return _positive_enclosed_output(evaluate, label)


@lru_cache(maxsize=2048)
def _exp_interval(value, precision):
    """Return directed decimal bounds of ``exp(value)`` at the given precision.

    Parameters
    ----------
    value : Fraction
        Exact exponent.
    precision : int
        Decimal working precision in significant digits.

    Returns
    -------
    tuple of Decimal
        ``(lower, upper)`` with ``lower <= exp(value) <= upper``. Monotonicity
        clamps the bounds to ``[1, inf)`` for positive and ``[0, 1]`` for negative
        exponents.

    References
    ----------
    .. [1] PYTHON-3.12-DECIMAL. Python 3.12 documentation, decimal module: precision
           contexts and exact conversion of binary floats.
    """
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
    """Export ``scale * exp(exponent)`` for a positive result.

    Parameters
    ----------
    scale : Fraction
        Exact positive multiplier.
    exponent : Fraction
        Exact exponent.
    label : str
        Output label used in errors.

    Returns
    -------
    float
        The correctly enclosed binary64 value.

    References
    ----------
    .. [1] PYTHON-3.12-DECIMAL. Python 3.12 documentation, decimal module: precision
           contexts and exact conversion of binary floats.
    """
    def evaluate(precision):
        math = _Intervals(precision)
        return math.mul(math.rational(scale), _exp_interval(exponent, precision))
    return _positive_enclosed_output(evaluate, label)


def _scaled_exp_complement(scale, exponent, label):
    """Export ``scale * (1 - exp(-exponent))`` without exporting either factor.

    Neither the survival term ``exp(-exponent)`` nor its complement is rounded
    before the product, so small exponents keep their relative accuracy.

    Parameters
    ----------
    scale : Fraction
        Exact positive multiplier.
    exponent : Fraction
        Exact positive exponent.
    label : str
        Output label used in errors.

    Returns
    -------
    float
        The correctly enclosed binary64 value.

    References
    ----------
    .. [1] PYTHON-3.12-DECIMAL. Python 3.12 documentation, decimal module: precision
           contexts and exact conversion of binary floats.
    """
    def evaluate(precision):
        math = _Intervals(precision)
        complement = math.sub((Decimal(1), Decimal(1)), _exp_interval(-exponent, precision))
        return math.mul(math.rational(scale), complement)
    return _positive_enclosed_output(evaluate, label)


def _inverse_mass_outputs(initial, exponent):
    # e > 2, so exp(x) > 2**1024 for x > 1024, beyond the binary64
    # overflow midpoint. This is an output impossibility proof, not a clamp.
    """Export the inverse rocket-equation mass ledger for a given exponent.

    Parameters
    ----------
    initial : Fraction
        Exact initial mass m0, kg.
    exponent : Fraction
        Exact ratio x = delta_v / c.

    Returns
    -------
    tuple of float
        ``(mass_ratio, final_mass, propellant_mass, propellant_fraction)``.

    Raises
    ------
    DomainError
        If the mass ratio exp(x) cannot be represented as a finite binary64 value.

    Notes
    -----
    MR = exp(x), mf = m0 * exp(-x), mp = m0 * (1 - exp(-x)) and mp / m0 =
    1 - exp(-x) [1]_. Each output is enclosed and rounded independently.

    References
    ----------
    .. [1] NASA Glenn Research Center, "Ideal Rocket Equation", Beginner's Guide to
           Aeronautics (Tsiolkovsky equation, mass ratio, effective exhaust velocity).
    .. [2] PYTHON-3.12-DECIMAL. Python 3.12 documentation, decimal module: precision
           contexts and exact conversion of binary floats.
    """
    if exponent > 1024:
        raise DomainError('Mass ratio is outside the finite output range.', field='Mass ratio')
    ratio = _scaled_exp(Fraction(1), exponent, 'Mass ratio')
    final = _scaled_exp(initial, -exponent, 'Final mass')
    consumed = _scaled_exp_complement(initial, exponent, 'Propellant mass')
    fraction = _scaled_exp_complement(Fraction(1), exponent, 'Propellant mass fraction')
    return ratio, final, consumed, fraction
