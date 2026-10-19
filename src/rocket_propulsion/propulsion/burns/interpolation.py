"""Shared shape-preserving interpolation kernels.

References
----------
F. N. Fritsch and R. E. Carlson, Monotone Piecewise Cubic Interpolation,
SIAM Journal on Numerical Analysis 17(2), 1980.
F. N. Fritsch and J. Butland, A Method for Constructing Local Monotone
Piecewise Cubic Interpolants, SIAM J. Sci. Stat. Comput. 5(2), 1984.
"""

from __future__ import annotations

from itertools import pairwise, product
from math import lcm, prod

from rocket_propulsion.core.errors import DomainError
from rocket_propulsion.core.validation import real_tuple

from .numeric_outputs import _finite_ratio, _finite_rational, _fraction

_MULTILINEAR_POLICY = "rational_multilinear_weights_and_exports_v1"


def _multilinear_terms(brackets):
    """Yield exact corner and derivative weights; convert only final fields."""
    fractions = tuple(_fraction(item[2]) for item in brackets)
    widths = tuple(_fraction(item[3]) for item in brackets)
    for corner in product((0, 1), repeat=len(brackets)):
        factors = tuple(fraction if upper else 1-fraction
                        for fraction, upper in zip(fractions, corner, strict=True))
        weight = prod(factors)
        gradients = tuple((1 if corner[axis] else -1)
                          * prod(value for index, value in enumerate(factors) if index != axis) / width
                          for axis, width in enumerate(widths))
        yield corner, weight, gradients


def pchip_slopes(x: tuple[float, ...], y: tuple[float, ...]) -> tuple[float, ...]:
    """Return Fritsch-Carlson / Fritsch-Butland monotone nodal slopes.

    Interior slopes are the weighted harmonic mean of adjacent secants (zero
    at a local extremum); end slopes use the shape-preserving three-point
    formula.  The resulting cubic Hermite interpolant is C1 and monotone
    wherever the data are.

    References
    ----------
    F. N. Fritsch and R. E. Carlson, Monotone Piecewise Cubic Interpolation,
    SIAM Journal on Numerical Analysis 17(2), 1980.
    """

    exact = _pchip_rational_slopes(real_tuple(x, "PCHIP coordinates"),
                                   real_tuple(y, "PCHIP ordinates"))
    return tuple(_finite_rational(value, f"PCHIP nodal slope ({index})")
                 for index, value in enumerate(exact))


def _pchip_rational_slopes(x, y):
    """Keep nodal and mixed derivatives exact until a consumer exports a field."""
    x = tuple(_fraction(value) for value in x)
    y = tuple(_fraction(value) for value in y)
    count = len(x)
    if count < 2 or count != len(y):
        raise DomainError("PCHIP needs at least two matching coordinate/ordinate nodes.")
    if any(b <= a for a, b in pairwise(x)):
        raise DomainError("PCHIP coordinates must be strictly increasing.")
    widths = [b - a for a, b in pairwise(x)]
    secants = [(y[i + 1] - y[i]) / widths[i] for i in range(count - 1)]
    if count == 2:
        return (secants[0], secants[0])
    slopes = [_fraction(0)] * count
    for i in range(1, count - 1):
        left, right = secants[i - 1], secants[i]
        if left == 0 or right == 0 or (left > 0) != (right > 0):
            slopes[i] = _fraction(0)
        else:
            w1 = 2 * widths[i] + widths[i - 1]
            w2 = widths[i] + 2 * widths[i - 1]
            slopes[i] = (w1 + w2) / (w1 / left + w2 / right)

    def end_slope(h0, h1, d0, d1):
        slope = ((2 * h0 + h1) * d0 - h0 * d1) / (h0 + h1)
        if d0 == 0 or slope == 0 or (slope > 0) != (d0 > 0):
            return _fraction(0)
        if (d0 == 0 or d1 == 0 or (d0 > 0) != (d1 > 0)) and abs(slope) > abs(3 * d0):
            return 3 * d0
        return slope

    slopes[0] = end_slope(widths[0], widths[1], secants[0], secants[1])
    slopes[-1] = end_slope(widths[-1], widths[-2], secants[-1], secants[-2])
    return tuple(slopes)


def _hermite_coefficients(x0, x1, y0, y1, m0, m1):
    """Compile exact power coefficients to integers with one common denominator."""
    x0, x1, y0, y1 = map(_fraction, (x0, x1, y0, y1))
    m0, m1 = _fraction(m0), _fraction(m1)
    h = x1 - x0
    delta = y1 - y0
    b = m0
    c = (3 * delta - h * (2 * m0 + m1)) / h**2
    d = (-2 * delta + h * (m0 + m1)) / h**3
    coefficients = (y0, b, c, d)
    denominator = lcm(*(value.denominator for value in coefficients))
    integers = tuple(value.numerator * (denominator // value.denominator) for value in coefficients)
    return (x0.numerator, x0.denominator, denominator, *integers)


def _hermite_segment(coefficients, query, *, derivative=True):
    """Evaluate the same rational polynomial using integer arithmetic per call."""
    base_numerator, base_denominator, denominator, a, b, c, d = coefficients
    query_numerator, query_denominator = query.as_integer_ratio()
    scale = max(base_denominator, query_denominator)  # binary input denominators are powers of two
    coordinate = (query_numerator * (scale // query_denominator)
                  - base_numerator * (scale // base_denominator))
    square, cube = scale**2, scale**3
    coordinate_square = coordinate**2
    value_numerator = (a * cube + b * coordinate * square
                       + c * coordinate_square * scale + d * coordinate_square * coordinate)
    value = _finite_ratio(value_numerator, denominator * cube, "Curve value")
    if not derivative:
        return value
    derivative_numerator = b * square + 2 * c * coordinate * scale + 3 * d * coordinate_square
    slope = _finite_ratio(derivative_numerator, denominator * square, "Curve derivative")
    return value, slope

