"""Rao-style thrust-optimized parabolic bell nozzle construction.

A quadratic approximation of the thrust-optimized contour [1]_, [2]_ between
prescribed throat and exit wall angles.

References
----------
.. [1] RAO-1958. G. V. R. Rao, "Exhaust Nozzle Contour for Optimum Thrust",
       Jet Propulsion 28(6), 377-382, 1958. doi:10.2514/8.7324
.. [2] SUTTON-BIBLARZ-2017. G. P. Sutton and O. Biblarz, "Rocket Propulsion
       Elements", 9th ed., Wiley, 2017.
"""

from __future__ import annotations

from fractions import Fraction as F
from math import cos, pi, radians, tan

from rocket_propulsion.compressible import calculate_isentropic, mach_from_isentropic
from rocket_propulsion.core.errors import DomainError
from rocket_propulsion.core.validation import integer, numeric_range_guard, real

from .dimensional_numerics import (
    dimensional,
    radial_growth,
    radius_exit,
    station_geometry,
    throat_radius,
    wall_area,
)
from .nozzle_contour import NozzleContour, NozzleStation, _checked_contour, _positive


@numeric_range_guard('geometry')
def generate_rao_contour(
    *,
    throat_area_m2: float,
    area_ratio: float,
    gamma: float = 1.22,
    length_fraction: float = 0.8,
    throat_angle_deg: float = 30.0,
    exit_angle_deg: float = 8.0,
    chamber_radius_ratio: float = 2.5,
    station_count: int = 121,
) -> NozzleContour:
    """Generate a quadratic Rao/TOP approximation with explicit end angles.

    Parameters
    ----------
    throat_area_m2 : float
        Throat area in m^2, finite and positive.
    area_ratio : float
        Exit-to-throat area ratio, greater than one.
    gamma : float, default 1.22
        Ratio of specific heats for the station solution.
    length_fraction : float, default 0.8
        Diverging length as a fraction of a 15 degree cone with the same area
        ratio, in ``[0.6, 1.0]``.
    throat_angle_deg : float, default 30.0
        Initial wall angle ``theta_n`` at the throat, in ``[15, 40]`` degrees.
    exit_angle_deg : float, default 8.0
        Exit wall angle ``theta_e``, non-negative and below ``theta_n``.
    chamber_radius_ratio : float, default 2.5
        Chamber-to-throat radius ratio, greater than one.
    station_count : int, default 121
        Total number of stations, between 41 and 401.

    Returns
    -------
    NozzleContour
        Contour with ``contour="rao"``, sampled isentropic stations, wall area
        and a divergence estimate.

    Raises
    ------
    DomainError
        If an input is outside its domain, or the throat and exit tangent lines
        do not intersect inside the diverging length
        (code ``invalid_rao_control_point``).

    Notes
    -----
    The diverging wall is a quadratic Bezier curve from the throat
    ``(0, r_t)`` to the exit ``(L_n, r_e)``. Its control point is the
    intersection of the throat tangent (slope ``tan(theta_n)``) and the exit
    tangent (slope ``tan(theta_e)``), so both end angles are met exactly. This
    is the standard thrust-optimized-parabola approximation of the
    thrust-optimized contour [1]_, [2]_, not a digitized reproduction of a
    particular Rao optimization chart and not a method-of-characteristics
    solution. The converging side is a cosine contraction of length
    ``1.5 r_t``. The divergence efficiency uses ``(1 + cos(alpha))/2`` [2]_
    with ``alpha`` the mean of the throat and exit angles, a screening value
    only.

    References
    ----------
    .. [1] RAO-1958. G. V. R. Rao, "Exhaust Nozzle Contour for Optimum Thrust",
           Jet Propulsion 28(6), 377-382, 1958. doi:10.2514/8.7324
    .. [2] SUTTON-BIBLARZ-2017. G. P. Sutton and O. Biblarz, "Rocket Propulsion
           Elements", 9th ed., Wiley, 2017.
    """

    throat_area_m2 = _positive(throat_area_m2, 'throat_area_m2')
    area_ratio = _positive(area_ratio, 'area_ratio')
    gamma = _positive(gamma, 'gamma')
    length_fraction = _positive(length_fraction, 'length_fraction')
    chamber_radius_ratio = _positive(chamber_radius_ratio, 'chamber_radius_ratio')
    throat_angle_deg = real(throat_angle_deg, 'throat_angle_deg')
    exit_angle_deg = real(exit_angle_deg, 'exit_angle_deg')
    station_count = integer(station_count, 'station_count')
    if gamma <= 1.0:
        raise DomainError('Gamma must exceed one.', field='gamma')
    if area_ratio <= 1.0:
        raise DomainError("Rao area ratio must exceed one.", field="area_ratio")
    if not 0.6 <= length_fraction <= 1.0:
        raise DomainError("Rao length fraction must be in [0.6, 1.0].", field='length_fraction')
    if not 15.0 <= throat_angle_deg <= 40.0:
        raise DomainError("Rao throat angle must be in [15, 40] degrees.", field='throat_angle_deg')
    if not 0.0 <= exit_angle_deg < throat_angle_deg:
        raise DomainError("Rao exit angle must be non-negative and below throat angle.", field='exit_angle_deg')
    if chamber_radius_ratio <= 1.0:
        raise DomainError("Chamber radius ratio must exceed one.", field='chamber_radius_ratio')
    if not 41 <= station_count <= 401:
        raise DomainError("Station count must be between 41 and 401.", field='station_count')

    radius_throat = throat_radius(throat_area_m2)
    radius_increment = radial_growth(radius_throat, area_ratio)
    exit_radius = radius_exit(throat_area_m2, area_ratio)
    chamber_radius = chamber_radius_ratio * radius_throat
    converging_length = 1.5 * radius_throat
    cone_length = dimensional(F(radius_increment)/F(tan(radians(15.0))), 'diverging_length_m')
    diverging_length = length_fraction * cone_length
    throat_slope = tan(radians(throat_angle_deg))
    exit_slope = tan(radians(exit_angle_deg))
    control_x = dimensional((F(radius_increment)-F(exit_slope)*F(diverging_length))
                            /(F(throat_slope)-F(exit_slope)), 'control_x_m')
    control_increment = dimensional(F(throat_slope)*F(control_x), 'control_radius_increment_m')
    if not 0.0 < control_x < diverging_length:
        raise DomainError(
            "Rao tangent lines do not intersect inside the nozzle length.",
            code="invalid_rao_control_point",
            details={"control_x_m": control_x, "diverging_length_m": diverging_length},
        )

    converging_count = max(15, round(station_count * 0.3))
    diverging_count = station_count - converging_count + 1
    stations: list[NozzleStation] = []

    def append(x_m: float, radius_m: float, branch: str, prescribed_ratio=None) -> None:
        area, local_ratio = station_geometry(radius_m, radius_throat, throat_area_m2,prescribed_ratio)
        mach = (
            1.0
            if local_ratio == 1.0
            else mach_from_isentropic(f"area_ratio_{branch}", local_ratio, gamma)
        )
        flow = calculate_isentropic(mach, gamma)
        stations.append(
            NozzleStation(
                x_m=x_m,
                radius_m=radius_m,
                area_m2=area,
                area_ratio=local_ratio,
                mach=mach,
                pressure_total_ratio=flow.pressure_total_ratio,
                temperature_total_ratio=flow.temperature_total_ratio,
            )
        )

    for index in range(converging_count):
        fraction = index / (converging_count - 1)
        radius = radius_throat + (chamber_radius-radius_throat) * (
            0.5 + 0.5 * cos(pi * fraction)
        )
        append(-converging_length * (1.0 - fraction), radius, "subsonic")
    for index in range(1, diverging_count):
        fraction = index / (diverging_count - 1)
        one_minus = 1.0 - fraction
        x_m = 2.0 * one_minus * fraction * control_x + fraction**2 * diverging_length
        radius = dimensional(F(radius_throat)+2*F(one_minus)*F(fraction)*F(control_increment)
                             + F(fraction)**2*F(radius_increment), 'radius_m')
        last = index == diverging_count-1
        append(x_m, exit_radius if last else radius, "supersonic", area_ratio if last else None)

    surface_area = wall_area(stations)
    effective_angle = 0.5 * (throat_angle_deg + exit_angle_deg)
    return _checked_contour(NozzleContour(
        contour="rao",
        throat_radius_m=radius_throat,
        chamber_radius_m=chamber_radius,
        exit_radius_m=exit_radius,
        converging_length_m=converging_length,
        diverging_length_m=diverging_length,
        total_length_m=converging_length + diverging_length,
        throat_area_m2=throat_area_m2,
        exit_area_m2=throat_area_m2 * area_ratio,
        area_ratio=area_ratio,
        half_angle_deg=effective_angle,
        stations=tuple(stations),
        model_fidelity="Rao thrust-optimized parabola approximation",
        wall_tangency_residual_deg=0.0,
        estimated_divergence_efficiency=0.5 * (1.0 + cos(radians(effective_angle))),
        surface_area_m2=surface_area,
        source="Rao-style quadratic tangent construction; NASA SP-8120 context",
    ))

