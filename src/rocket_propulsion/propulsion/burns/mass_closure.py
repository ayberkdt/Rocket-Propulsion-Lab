"""Exact supplied-binary inventory arithmetic and portable decision evidence."""

from fractions import Fraction

from rocket_propulsion.core.errors import DomainError

from .numeric_outputs import _finite_rational, _fraction

_MASS_CLOSURE_POLICY = "exact_binary_inventory_residual_and_ratio_decision_v2"


def _inventory_ledger(vehicle_mass, other_mass, inventories):
    """Sum selected states before exporting; retain the unrounded residual ratio."""
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
    """Validate a supplied numerical ratio, or use a manual record's point value."""
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
