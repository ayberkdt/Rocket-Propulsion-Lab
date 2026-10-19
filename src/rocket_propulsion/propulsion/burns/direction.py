"""Internal scale-safe direction projections and signed acceleration tangents.

Rational intermediates preserve the supplied binary coordinates. The scaled
hypot norm remains a floating-point estimate, not an exact square root.
"""

from dataclasses import dataclass
from fractions import Fraction
from math import hypot

from rocket_propulsion.core.errors import DomainError
from rocket_propulsion.core.validation import vector3

from .numeric_outputs import _finite_rational, _fraction


@dataclass(frozen=True, slots=True)
class _DirectionGeometry:
    coordinates: tuple[Fraction, Fraction, Fraction]
    squared_norm: Fraction
    norm_estimate: Fraction

    def project(self, coefficient, divisor=1.0, *, label: str):
        """Convert only the final vector components, preserving rescued terms."""
        denominator = _fraction(divisor)
        if denominator <= 0:
            raise DomainError("Direction projection divisor must be positive.")
        gain = _fraction(coefficient) / (denominator * self.norm_estimate)
        return tuple(_finite_rational(gain * value, f"{label} ({index})")
                     for index, value in enumerate(self.coordinates))

    def jacobian(self, thrust, mass):
        """Return d(T d/(m |d|))/dd, retaining the sign of axial thrust T."""
        vehicle_mass = _fraction(mass)
        if vehicle_mass <= 0:
            raise DomainError("Direction Jacobian mass must be positive.")
        gain = _fraction(thrust) / (vehicle_mass * self.squared_norm * self.norm_estimate)
        return tuple(tuple(_finite_rational(
            gain * ((self.squared_norm if row == column else 0)
                    - self.coordinates[row] * self.coordinates[column]),
            f"Direction Jacobian ({row}, {column})")
            for column in range(3)) for row in range(3))


def _direction_geometry(direction):
    """Keep original coordinates even when a unit component cannot be stored."""
    vector = vector3(direction, "Thrust direction")
    scale = max(abs(value) for value in vector)
    if scale == 0:
        raise DomainError("Thrust direction must be non-zero.")
    scaled_norm = hypot(*(value / scale for value in vector))
    coordinates = tuple(_fraction(value) for value in vector)
    return _DirectionGeometry(coordinates, sum((value * value for value in coordinates), Fraction(0)),
                              _fraction(scale) * _fraction(scaled_norm))
