"""Shared multilinear and shape-preserving cubic interpolation kernels.

The kernels evaluate tensor-product multilinear weights on rectilinear grids
[4]_ and monotone piecewise cubic Hermite (PCHIP) slopes [1]_, [2]_, [3]_ in
exact rational arithmetic [5]_.  Binary64 inputs are converted to their exact
rational values; only exported fields are rounded, once, to binary64.

References
----------
.. [1] F. N. Fritsch and R. E. Carlson, "Monotone Piecewise Cubic
       Interpolation", SIAM J. Numer. Anal. 17(2), 1980.
.. [2] F. N. Fritsch and J. Butland, "A Method for Constructing Local
       Monotone Piecewise Cubic Interpolants", SIAM J. Sci. Stat. Comput.
       5(2), 1984.
.. [3] SCIPY-1.14.1-PCHIP. SciPy 1.14.1, documentation and implementation
       of "scipy.interpolate.PchipInterpolator".
.. [4] SCIPY-1.14.1-REGULAR-GRID. SciPy 1.14.1, documentation of
       "scipy.interpolate.RegularGridInterpolator".
.. [5] PYTHON-3.12-FRACTION. Python 3.12, "fractions" module.
       https://docs.python.org/3.12/library/fractions.html
"""

from __future__ import annotations

from itertools import pairwise, product
from math import lcm, prod

from rocket_propulsion.core.errors import DomainError
from rocket_propulsion.core.validation import real_tuple

from .numeric_outputs import _finite_ratio, _finite_rational, _fraction

_MULTILINEAR_POLICY = "rational_multilinear_weights_and_exports_v1"


def _multilinear_terms(brackets):
    """Yield exact multilinear corner weights and their coordinate gradients.

    Parameters
    ----------
    brackets : sequence
        One entry per axis; item ``[2]`` is the local fraction ``t`` in
        ``[0, 1]`` and item ``[3]`` the bracket width ``h``.

    Yields
    ------
    tuple
        ``(corner, weight, gradients)``.  ``corner`` is a tuple of 0/1 flags;
        ``weight`` is the tensor product of ``1 - t`` (lower) or ``t`` (upper)
        factors [1]_; ``gradients[k]`` is the derivative of ``weight`` with
        respect to axis ``k``, ``+-prod(other factors) / h_k``.  Values are
        exact `fractions.Fraction` objects; rounding is left to the consumer.

    References
    ----------
    .. [1] SCIPY-1.14.1-REGULAR-GRID. SciPy 1.14.1, documentation of
           "scipy.interpolate.RegularGridInterpolator".
    """
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

    Parameters
    ----------
    x : tuple of float
        Strictly increasing node coordinates, at least two.
    y : tuple of float
        Node ordinates, same length as ``x``.

    Returns
    -------
    tuple of float
        Nodal slopes of the monotone cubic Hermite interpolant, each rounded
        once from its exact rational value.

    Raises
    ------
    DomainError
        If fewer than two nodes are given, the lengths differ, the coordinates
        are not strictly increasing, or a slope is not representable as a
        finite binary64 value.

    Notes
    -----
    The slopes follow the monotone PCHIP construction [1]_, [2]_, [3]_.
    Interior slopes are the weighted harmonic mean of the adjacent secants,
    with weights ``w1 = 2 h_k + h_(k-1)`` and ``w2 = h_k + 2 h_(k-1)``; the
    slope is zero where the secants change sign or either is zero.  End slopes
    use the shape-preserving three-point formula, zeroed when its sign differs
    from the end secant and limited to three times the end secant when the two
    end secants do not share a sign.  Two nodes give the common secant at both
    ends.  The resulting cubic Hermite interpolant is C1 and monotone wherever
    the data are.

    References
    ----------
    .. [1] F. N. Fritsch and R. E. Carlson, "Monotone Piecewise Cubic
           Interpolation", SIAM J. Numer. Anal. 17(2), 1980.
    .. [2] F. N. Fritsch and J. Butland, "A Method for Constructing Local
           Monotone Piecewise Cubic Interpolants", SIAM J. Sci. Stat. Comput.
           5(2), 1984.
    .. [3] SCIPY-1.14.1-PCHIP. SciPy 1.14.1, documentation and implementation
           of "scipy.interpolate.PchipInterpolator".
    """

    exact = _pchip_rational_slopes(real_tuple(x, "PCHIP coordinates"),
                                   real_tuple(y, "PCHIP ordinates"))
    return tuple(_finite_rational(value, f"PCHIP nodal slope ({index})")
                 for index, value in enumerate(exact))


def _pchip_rational_slopes(x, y):
    """Return exact rational PCHIP nodal slopes.

    Same algorithm as `pchip_slopes` [1]_, [2]_, [3]_, but the inputs are
    converted to `fractions.Fraction` and the slopes are returned unrounded,
    so that consumers building nodal or mixed derivatives round only their
    final exported fields.

    Raises
    ------
    DomainError
        If fewer than two nodes are given, the lengths differ, or the
        coordinates are not strictly increasing.

    References
    ----------
    .. [1] F. N. Fritsch and R. E. Carlson, "Monotone Piecewise Cubic
           Interpolation", SIAM J. Numer. Anal. 17(2), 1980.
    .. [2] F. N. Fritsch and J. Butland, "A Method for Constructing Local
           Monotone Piecewise Cubic Interpolants", SIAM J. Sci. Stat. Comput.
           5(2), 1984.
    .. [3] SCIPY-1.14.1-PCHIP. SciPy 1.14.1, documentation and implementation
           of "scipy.interpolate.PchipInterpolator".
    """
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
    """Compile one cubic Hermite segment to integer power-basis coefficients.

    With ``h = x1 - x0`` and ``delta = y1 - y0`` the segment is
    ``p(s) = y0 + m0 s + c s**2 + d s**3`` for ``s = x - x0``, where
    ``c = (3 delta - h (2 m0 + m1)) / h**2`` and
    ``d = (-2 delta + h (m0 + m1)) / h**3`` [1]_.  The four exact rational
    coefficients are scaled to integers over their least common denominator.

    Returns
    -------
    tuple of int
        ``(x0_numerator, x0_denominator, common_denominator, a, b, c, d)``.

    References
    ----------
    .. [1] F. N. Fritsch and R. E. Carlson, "Monotone Piecewise Cubic
           Interpolation", SIAM J. Numer. Anal. 17(2), 1980.
    """
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
    """Evaluate a compiled cubic Hermite segment exactly in integer arithmetic.

    Parameters
    ----------
    coefficients : tuple of int
        Output of `_hermite_coefficients`.
    query : float
        Binary64 abscissa; its exact rational value is used.
    derivative : bool, default True
        Also return the first derivative.

    Returns
    -------
    float or tuple of float
        The segment value, or ``(value, slope)``, each rounded once to binary64.

    Raises
    ------
    DomainError
        If the value or slope is not representable as a finite binary64 value.

    Notes
    -----
    Binary64 denominators are powers of two, so the larger of the query and
    base denominators is a common denominator of both; the polynomial and its
    derivative are each evaluated as one integer ratio before rounding.
    """
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

