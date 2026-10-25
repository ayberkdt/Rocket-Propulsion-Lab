"""Internal scale-safe thrust-direction projections and acceleration Jacobians.

A thrust direction is kept as exact rationals of the supplied binary
coordinates [2]_; only final vector components are rounded to floats.  The
norm is estimated as ``s * hypot(d / s)`` with ``s = max |d_i|`` [1]_, which
avoids overflow and underflow but remains a floating-point estimate, not an
exact square root.

References
----------
.. [1] PYTHON-3.12-HYPOT. Python 3.12, "math.hypot", Euclidean norm and
       floating-point accuracy.
       https://docs.python.org/3.12/library/math.html#math.hypot
.. [2] PYTHON-3.12-FRACTION. Python 3.12, "fractions"; float inputs retain
       their exact binary value as a rational.
       https://docs.python.org/3.12/library/fractions.html
"""

from dataclasses import dataclass
from fractions import Fraction
from math import hypot

from rocket_propulsion.core.errors import DomainError
from rocket_propulsion.core.validation import vector3

from .numeric_outputs import _finite_rational, _fraction


@dataclass(frozen=True, slots=True)
class _DirectionGeometry:
    """Exact direction coordinates with their squared norm and norm estimate.

    Attributes
    ----------
    coordinates : tuple of Fraction
        Exact rationals of the supplied coordinates.
    squared_norm : Fraction
        Exact ``|d|**2``.
    norm_estimate : Fraction
        Rational value of the scaled ``hypot`` estimate of ``|d|``.
    """

    coordinates: tuple[Fraction, Fraction, Fraction]
    squared_norm: Fraction
    norm_estimate: Fraction

    def project(self, coefficient, divisor=1.0, *, label: str):
        """Return ``coefficient * d / (divisor * |d|)`` as finite floats.

        The gain and products stay exact rationals, so components that would
        underflow or overflow in an intermediate float are preserved; only the
        three final components are rounded.  With ``coefficient`` a thrust and
        ``divisor`` a mass this is the thrust acceleration vector [1]_.

        Raises
        ------
        DomainError
            If ``divisor`` is not positive or a component is not a finite
            float.

        References
        ----------
        .. [1] MIT-522-L5. MIT OCW, "16.522 Session 5: Review of Classical
               Astrodynamics", equation (19), thrust acceleration F/m.
        """
        denominator = _fraction(divisor)
        if denominator <= 0:
            raise DomainError("Direction projection divisor must be positive.")
        gain = _fraction(coefficient) / (denominator * self.norm_estimate)
        return tuple(_finite_rational(gain * value, f"{label} ({index})")
                     for index, value in enumerate(self.coordinates))

    def jacobian(self, thrust, mass):
        """Return the Jacobian of the thrust acceleration by the direction vector.

        For ``a = T d / (m |d|)`` [1]_,
        ``da_i/dd_j = T (|d|**2 delta_ij - d_i d_j) / (m |d|**3)``; the sign of
        the axial thrust ``T`` is retained.  The matrix is evaluated in exact
        rationals and rounded once per entry.

        Raises
        ------
        DomainError
            If ``mass`` is not positive or an entry is not a finite float.

        References
        ----------
        .. [1] MIT-522-L5. MIT OCW, "16.522 Session 5: Review of Classical
               Astrodynamics", equation (19), thrust acceleration F/m.
        """
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
    """Validate a thrust direction and build its exact geometry.

    The original coordinates are retained even when a unit-vector component
    would not be representable as a float.

    Raises
    ------
    DomainError
        If the direction is not three finite reals or is the zero vector.
    """
    vector = vector3(direction, "Thrust direction")
    scale = max(abs(value) for value in vector)
    if scale == 0:
        raise DomainError("Thrust direction must be non-zero.")
    scaled_norm = hypot(*(value / scale for value in vector))
    coordinates = tuple(_fraction(value) for value in vector)
    return _DirectionGeometry(coordinates, sum((value * value for value in coordinates), Fraction(0)),
                              _fraction(scale) * _fraction(scaled_norm))
