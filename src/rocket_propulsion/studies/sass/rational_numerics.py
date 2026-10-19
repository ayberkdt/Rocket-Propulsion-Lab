"""Final binary conversion for exact rational algebraic ledgers.

References
----------
NASA-7009B: numerical representation limits, distinct from physical uncertainty.
BIPM-SI: retain the declared units of each physical result.
"""

from fractions import Fraction
from math import isfinite

from rocket_propulsion.core.errors import DomainError


def binary_ledger(rows: tuple[tuple[str, Fraction | int, str], ...], label: str) -> tuple:
    """Round final channels once; reject nonzero values rounded to zero/infinity.

    Inputs to this internal converter are already validated and computed in
    rational arithmetic. It checks output representability, not model validity.

    References
    ----------
    NASA-7009B: declared numerical domain and output conversion boundaries.
    BIPM-SI: preserve physical units through numerical evaluation.
    """
    converted = []
    for name, value, unit in rows:
        try:
            binary = float(value)
        except OverflowError as error:
            raise DomainError(f"{label} {name} exceeds the finite binary output domain.") from error
        if not isfinite(binary) or (value != 0 and binary == 0):
            raise DomainError(f"{label} {name} is outside the representable binary output domain.")
        converted.append((name, binary, unit))
    return tuple(converted)
