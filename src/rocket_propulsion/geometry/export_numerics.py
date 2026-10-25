"""Exact vertex arithmetic for text mesh exports, not geometric error bounds.

The helpers check that redundant station quantities agree within their binary64
representation cells [1]_, [2]_ and build STL facet vertices and unit normals
[3]_ with exact rational arithmetic on the serialized coordinates.

References
----------
.. [1] PYTHON-3.12-FRACTION. Python, ``fractions``; float inputs retain their
       exact binary value as a rational.
       https://docs.python.org/3.12/library/fractions.html
.. [2] PYTHON-3.12-NEXTAFTER. Python, ``math.nextafter``, adjacent representable
       floating-point values.
       https://docs.python.org/3.12/library/math.html#math.nextafter
.. [3] LOC-STL-FAMILY. Library of Congress, STL format family: oriented
       triangular surface and format limitations.
"""
from fractions import Fraction
from math import cos, inf, isfinite, nextafter, pi, sin

from rocket_propulsion.core.errors import DomainError
from rocket_propulsion.propulsion.burns.numeric_outputs import _finite_rational
from rocket_propulsion.propulsion.burns.probability_numerics import _sqrt_output


def _rounding_cell(value):
    """Return the closed round-to-nearest-even cell of a positive binary64 value.

    Parameters
    ----------
    value : float
        Finite positive float.

    Returns
    -------
    tuple of Fraction
        ``(lower, upper)``: the midpoints to the adjacent representable values
        [1]_, [2]_. When the upper neighbour overflows, the upper bound mirrors
        the lower half-spacing. The cell is a representation interval, not a
        measurement uncertainty interval.

    References
    ----------
    .. [1] PYTHON-3.12-FRACTION. Python, ``fractions``; float inputs retain their
           exact binary value as a rational.
           https://docs.python.org/3.12/library/fractions.html
    .. [2] PYTHON-3.12-NEXTAFTER. Python, ``math.nextafter``, adjacent representable
           floating-point values.
           https://docs.python.org/3.12/library/math.html#math.nextafter
    """
    center = Fraction(value)
    before = Fraction(nextafter(value, -inf))
    after = nextafter(value, inf)
    return ((center + before) / 2,
            (center + Fraction(after)) / 2 if isfinite(after) else center + (center - before) / 2)


def _check_geometry(stations, throat_area, exit_ratio):
    """Check that station radius, area and area ratio agree within rounding cells.

    For each station, ``pi r^2`` (with ``pi`` the binary64 constant) evaluated
    over the radius cell must overlap the area cell and not fall below the
    throat-area cell, and the area-ratio cell must overlap the implied
    ``pi r^2 / A_t`` interval. The last station's area ratio must overlap the
    nominal exit area ratio.

    This checks redundant station geometry, not flow validity or the accuracy
    of a continuous contour between stations. Cell endpoints are conservatively
    closed; no physical measurement tolerance is implied.

    Raises
    ------
    DomainError
        If any of the redundant quantities disagree beyond representation
        cells.
    """
    at_lo, at_hi = _rounding_cell(throat_area)
    fixed_pi = Fraction(pi)
    for station in stations:
        r_lo, r_hi = _rounding_cell(station.radius_m)
        a_lo, a_hi = _rounding_cell(station.area_m2)
        ar_lo, ar_hi = _rounding_cell(station.area_ratio)
        geom_lo, geom_hi = fixed_pi * r_lo**2, fixed_pi * r_hi**2
        if a_hi < geom_lo or a_lo > geom_hi or geom_hi < at_lo:
            raise DomainError('Area, radius and throat area disagree beyond representation cells.', field='Export geometry')
        if ar_hi < geom_lo / at_hi or ar_lo > geom_hi / at_lo:
            raise DomainError('Area ratio and radius disagree beyond representation cells.', field='Export geometry')
    ratio_lo, ratio_hi = _rounding_cell(exit_ratio)
    last_lo, last_hi = _rounding_cell(stations[-1].area_ratio)
    if ratio_hi < last_lo or ratio_lo > last_hi:
        raise DomainError('Exit station and nominal exit area ratio disagree.', field='Export geometry')


def _ring(station, count):
    """Return ``count`` vertices of the station circle in the y-z plane.

    Angles that are multiples of 90 degrees use exact unit cosines and sines.
    Other angles use the platform binary64 sine and cosine; their geometric
    sampling error is not bounded by this export policy. Each coordinate is the
    exact product ``r cos(phi)`` or ``r sin(phi)`` rounded once to binary64.
    """
    points = []
    for index in range(count):
        if (4 * index) % count == 0:
            cosine, sine = ((1., 0.), (0., 1.), (-1., 0.), (0., -1.))[(4 * index) // count]
        else:
            angle = 2.0 * pi * index / count
            cosine, sine = cos(angle), sin(angle)
        y = _finite_rational(Fraction(station.radius_m) * Fraction(cosine), 'STL y coordinate')
        z = _finite_rational(Fraction(station.radius_m) * Fraction(sine), 'STL z coordinate')
        points.append((station.x_m, y, z))
    return points


def _unit_normal(a, b, c):
    """Return the unit facet normal from the exact cross product of the vertices.

    The cross product ``(b - a) x (c - a)`` is formed exactly from the
    serialized binary64 coordinates; each component of the normalized vector
    is rounded once to nearest-even. The right-hand vertex order defines the
    outward orientation of an STL facet [1]_.

    Raises
    ------
    DomainError
        If the three vertices are collinear or coincident.

    References
    ----------
    .. [1] LOC-STL-FAMILY. Library of Congress, STL format family: oriented
           triangular surface and format limitations.
    """
    ab = tuple(Fraction(right) - Fraction(left) for left, right in zip(a, b))
    ac = tuple(Fraction(right) - Fraction(left) for left, right in zip(a, c))
    normal = (ab[1] * ac[2] - ab[2] * ac[1],
              ab[2] * ac[0] - ab[0] * ac[2],
              ab[0] * ac[1] - ab[1] * ac[0])
    squared = sum((value * value for value in normal), Fraction(0))
    if squared == 0:
        raise DomainError('Serialized STL facet is degenerate.', field='STL facet')
    return tuple((1 if value > 0 else -1) * _sqrt_output(value * value / squared, 'STL normal component')
                 if value else 0.0 for value in normal)
