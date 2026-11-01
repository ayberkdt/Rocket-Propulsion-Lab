"""Complete dimensional expressions on sampled binary64 nozzle coordinates.

Exact rational products/sums and the shared correctly rounded rational square
root avoid erasing small areas or overflowing private segment squares. Pi is
the binary64 math.pi basis also used by geometry exchange; these arithmetic
checks are not physical uncertainty or continuous-wall error bounds.

References
----------
.. [1] PYTHON-3.12-FRACTION. https://docs.python.org/3.12/library/fractions.html
.. [2] NACA-TR-1135. https://ntrs.nasa.gov/citations/19930091059
"""
from fractions import Fraction as F
from itertools import pairwise
from math import pi, sqrt

from rocket_propulsion.core.errors import DomainError
from rocket_propulsion.propulsion.burns.numeric_outputs import _finite_rational
from rocket_propulsion.propulsion.burns.probability_numerics import _sqrt_output


def dimensional(value, label):
    """Export a complete rational expression without nonfinite/false-zero output."""
    try:
        return _finite_rational(value, label)
    except DomainError as error:
        raise DomainError(error.message, code='numeric_range', field=label) from error


def throat_radius(area):
    """Round sqrt(A/pi) once, without rounding the private A/pi first."""
    return _sqrt_output(F(area)/F(pi), 'throat_radius_m')


def radius_exit(area, ratio):
    """Round the nominal exit radius from complete At*epsilon/pi."""
    return _sqrt_output(F(area)*F(ratio)/F(pi), 'exit_radius_m')


def radial_growth(radius, area_ratio):
    """Evaluate rt*(epsilon-1)/(sqrt(epsilon)+1) before rounding the increment."""
    return dimensional(F(radius)*(F(area_ratio)-1)/(F(sqrt(area_ratio))+1), 'radial_growth_m')


def station_geometry(radius, throat, throat_area, prescribed_ratio=None):
    """Complete local geometry with exact design throat/exit boundary ratios.

    Interior ratios use the nominal At, not a second rounded throat inferred
    from rt. The constructed throat and exit retain their defining ratios;
    those are boundary data, not clamps of arbitrary caller geometry.
    """
    area = F(pi)*F(radius)**2
    ratio = (prescribed_ratio if prescribed_ratio is not None else 1.0 if radius == throat
             else dimensional(area/F(throat_area), 'area_ratio'))
    return dimensional(area, 'area_m2'), ratio


def wall_area(stations):
    """Sum all frusta before exporting; private positive segment areas may underflow."""
    area = F(0)
    for left, right in pairwise(stations):
        dx, dr = F(right.x_m)-F(left.x_m), F(right.radius_m)-F(left.radius_m)
        length = _sqrt_output(dx*dx+dr*dr, 'segment_length_m')
        area += F(pi)*(F(left.radius_m)+F(right.radius_m))*F(length)
    return dimensional(area, 'surface_area_m2')
