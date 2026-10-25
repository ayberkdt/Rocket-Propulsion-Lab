"""Exact rational intermediates and checked binary64 outputs for propulsion ledgers.

Every supplied float is converted to the rational number it represents exactly
[1]_, so intermediate sums and products carry no rounding error. Exactness is
relative to the supplied binary input values, not to measured physics. Groups of
terms are summed before the single output conversion so that small terms can
combine or cancel. A nonzero value that would overflow or round to zero on export
is reported as a domain failure instead of being returned silently.

References
----------
.. [1] PYTHON-3.12-FRACTION. Python 3.12, ``fractions``: float inputs retain
       their exact binary value as a rational.
       https://docs.python.org/3.12/library/fractions.html
"""

from fractions import Fraction
from math import isfinite

from rocket_propulsion.core.errors import DomainError
from rocket_propulsion.core.validation import real


def _fraction(value, label: str = "Propulsion arithmetic input") -> Fraction:
    """Return ``value`` as an exact rational.

    Parameters
    ----------
    value : Fraction or real
        An existing exact intermediate is returned unchanged; any other value is
        validated as a finite real scalar and converted exactly [1]_.
    label : str, optional
        Field name used in validation errors.

    Returns
    -------
    Fraction
        The exact rational value.

    Raises
    ------
    DomainError
        If a non-``Fraction`` value is not a finite real scalar.

    References
    ----------
    .. [1] PYTHON-3.12-FRACTION. Python 3.12, ``fractions``: float inputs retain
           their exact binary value as a rational.
           https://docs.python.org/3.12/library/fractions.html
    """
    return value if isinstance(value, Fraction) else Fraction(real(value, label))


def _finite_rational(value: Fraction, label: str) -> float:
    """Round an exact rational once to a finite binary64 float.

    Parameters
    ----------
    value : Fraction
        Exact value to export.
    label : str
        Output name used in error messages and the error ``field``.

    Returns
    -------
    float
        ``float(value)``.

    Raises
    ------
    DomainError
        If the value overflows the finite binary64 range, or if a nonzero value
        rounds to zero (underflow erases it).
    """
    try:
        result = float(value)
    except OverflowError as error:
        raise DomainError(f"{label} is outside the finite output range.", field=label) from error
    if not isfinite(result) or (value != 0 and result == 0.0):
        raise DomainError(f"{label} is outside the nonzero finite output range.", field=label)
    return result


def _finite_ratio(numerator: int, denominator: int, label: str) -> float:
    """Export ``numerator / denominator`` without building a reduced ``Fraction``.

    Parameters
    ----------
    numerator : int
        Integer numerator.
    denominator : int
        Positive integer denominator.
    label : str
        Output name used in error messages and the error ``field``.

    Returns
    -------
    float
        Python true division of the two integers.

    Raises
    ------
    DomainError
        If the denominator is not positive, the quotient overflows the finite
        binary64 range, or a nonzero quotient rounds to zero.
    """
    if denominator <= 0:
        raise DomainError(f"{label} requires a positive denominator.", field=label)
    try:
        result = numerator / denominator
    except OverflowError as error:
        raise DomainError(f"{label} is outside the finite output range.", field=label) from error
    if not isfinite(result) or (numerator != 0 and result == 0.0):
        raise DomainError(f"{label} is outside the nonzero finite output range.", field=label)
    return result
