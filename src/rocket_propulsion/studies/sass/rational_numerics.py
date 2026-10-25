"""Final binary conversion for exact rational algebraic ledgers.

References
----------
.. [1] NASA-7009B. NASA, "Standard for Models and Simulations, NASA-STD-7009B";
       numerical representation limits, distinct from physical uncertainty.
.. [2] BIPM-SI. BIPM, "SI Brochure, 9th edition"; retain the declared units of each
       physical result.
"""

from fractions import Fraction
from math import isfinite

from rocket_propulsion.core.errors import DomainError


def binary_ledger(rows: tuple[tuple[str, Fraction | int, str], ...], label: str) -> tuple:
    """Round exact ledger channels to binary64 once, refusing lost nonzero values.

    Parameters
    ----------
    rows : tuple of (str, Fraction or int, str)
        ``(name, exact value, unit)`` channels computed in rational arithmetic.
    label : str
        Ledger label used in errors.

    Returns
    -------
    tuple of (str, float, str)

    Raises
    ------
    DomainError
        If a channel overflows, is not finite, or is nonzero but rounds to zero.

    Notes
    -----
    This checks output representability [1]_, not model validity; units are carried
    unchanged [2]_.

    References
    ----------
    .. [1] NASA-7009B. NASA, "Standard for Models and Simulations, NASA-STD-7009B";
           declared numerical domain and output conversion boundaries.
    .. [2] BIPM-SI. BIPM, "SI Brochure, 9th edition"; preserve physical units through
           numerical evaluation.
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
