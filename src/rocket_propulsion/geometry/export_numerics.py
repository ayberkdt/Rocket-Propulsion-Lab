"""Exact vertex arithmetic for text mesh exports, not geometric error bounds."""
from fractions import Fraction
from math import cos, inf, isfinite, nextafter, pi, sin

from rocket_propulsion.core.errors import DomainError
from rocket_propulsion.propulsion.burns.numeric_outputs import _finite_rational
from rocket_propulsion.propulsion.burns.probability_numerics import _sqrt_output


def _rounding_cell(value):
    """Positive binary64 RNE cell closure; not a measurement uncertainty interval."""
    center = Fraction(value)
    before = Fraction(nextafter(value, -inf))
    after = nextafter(value, inf)
    return ((center + before) / 2,
            (center + Fraction(after)) / 2 if isfinite(after) else center + (center - before) / 2)


def _check_geometry(stations, throat_area, exit_ratio):
    """Necessary overlap of representation cells under the fixed binary64 pi basis.

    This checks redundant station geometry, not flow validity or the accuracy
    of a continuous contour between stations. Cell endpoints are conservatively
    closed; no physical measurement tolerance is implied.
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
    """Sample a circle, preserving exact cardinal axes and checked products.

    Other angles use platform binary64 sine/cosine. Their geometric sampling
    error is not bounded by this export policy.
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
    """RNE components from the exact cross product of the serialized vertices."""
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
