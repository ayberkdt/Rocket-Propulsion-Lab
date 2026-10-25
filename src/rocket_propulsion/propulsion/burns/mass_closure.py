"""Exact mass-closure residuals for vehicle mass versus tracked propellant.

The closure residual is ``m_vehicle - m_other - sum(m_propellant)``; it is
formed from the exact rational values of the supplied binary64 masses [1]_ and
exported once, together with the unrounded residual as an integer ratio so a
consumer can make the same sign decision portably.

References
----------
.. [1] PYTHON-3.12-FRACTION. Python 3.12, "fractions": float inputs retain
       their exact binary value as a rational,
       https://docs.python.org/3.12/library/fractions.html.
"""

from fractions import Fraction

from rocket_propulsion.core.errors import DomainError

from .numeric_outputs import _finite_rational, _fraction

_MASS_CLOSURE_POLICY = "exact_binary_inventory_residual_and_ratio_decision_v2"


def _inventory_ledger(vehicle_mass, other_mass, inventories):
    """Return tracked mass, residual, relative residual and the exact residual ratio.

    The tracked inventories are summed exactly before export and the residual
    ``vehicle - other - tracked`` keeps its exact rational value [1]_ as
    ``(numerator, denominator)`` strings.

    Raises
    ------
    DomainError
        If the vehicle mass is not positive, the other mass is negative, or an
        exported value is not a finite float.

    References
    ----------
    .. [1] PYTHON-3.12-FRACTION. Python 3.12, "fractions": float inputs retain
           their exact binary value as a rational,
           https://docs.python.org/3.12/library/fractions.html.
    """
    vehicle, other = _fraction(vehicle_mass), _fraction(other_mass)
    if vehicle <= 0 or other < 0:
        raise DomainError("Vehicle mass must be positive and other mass non-negative.")
    tracked = sum((_fraction(value) for value in inventories), Fraction(0))
    residual = vehicle - other - tracked
    return (_finite_rational(tracked, "Tracked propellant mass"),
            _finite_rational(residual, "Mass closure residual"),
            _finite_rational(residual / vehicle, "Relative mass closure residual"),
            (str(residual.numerator), str(residual.denominator)))


def _closure_ratio(value, rounded_residual):
    """Return the exact residual ratio, or the reported residual when none is supplied.

    Raises
    ------
    DomainError
        If the ratio is not two integers with a positive denominator, or its
        nearest float differs from ``rounded_residual``.
    """
    label = "Mass closure residual ratio"
    if value is None:
        return _fraction(rounded_residual)
    if not isinstance(value, (tuple, list)) or len(value) != 2:
        raise DomainError("Residual ratio must contain two integers.", field=label)
    if any(type(part) not in (str, int) for part in value):
        raise DomainError("Residual ratio must contain integer text or integers.", field=label)
    try:
        numerator, denominator = map(int, value)
        if denominator <= 0:
            raise ValueError("Nonpositive denominator")
        ratio = Fraction(numerator, denominator)
    except (ValueError, TypeError, ZeroDivisionError) as error:
        raise DomainError("Residual ratio must have a positive integer denominator.", field=label) from error
    if _finite_rational(ratio, label) != rounded_residual:
        raise DomainError("Residual ratio contradicts the reported residual.", field=label)
    return ratio
