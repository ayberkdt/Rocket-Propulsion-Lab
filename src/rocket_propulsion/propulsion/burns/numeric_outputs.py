"""Internal exact arithmetic and checked binary outputs for propulsion ledgers.

Exactness is relative to supplied binary input values, not measured physics.
Groups must be summed before output conversion so small terms can combine or
cancel. Unrepresentable nonzero exported values are explicit domain failures.
"""

from fractions import Fraction
from math import isfinite

from rocket_propulsion.core.errors import DomainError
from rocket_propulsion.core.validation import real


def _fraction(value, label: str = "Propulsion arithmetic input") -> Fraction:
    """Preserve existing exact intermediates and validate external scalar factors."""
    return value if isinstance(value, Fraction) else Fraction(real(value, label))


def _finite_rational(value: Fraction, label: str) -> float:
    """Round once to binary float; reject overflow and erased nonzero outputs."""
    try:
        result = float(value)
    except OverflowError as error:
        raise DomainError(f"{label} is outside the finite output range.", field=label) from error
    if not isfinite(result) or (value != 0 and result == 0.0):
        raise DomainError(f"{label} is outside the nonzero finite output range.", field=label)
    return result


def _finite_ratio(numerator: int, denominator: int, label: str) -> float:
    """Export a compiled integer ratio without constructing/reducing a Fraction."""
    if denominator <= 0:
        raise DomainError(f"{label} requires a positive denominator.", field=label)
    try:
        result = numerator / denominator
    except OverflowError as error:
        raise DomainError(f"{label} is outside the finite output range.", field=label) from error
    if not isfinite(result) or (numerator != 0 and result == 0.0):
        raise DomainError(f"{label} is outside the nonzero finite output range.", field=label)
    return result
