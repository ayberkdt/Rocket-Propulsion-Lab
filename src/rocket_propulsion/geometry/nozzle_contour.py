"""Parametric axisymmetric converging-diverging nozzle contours.

Contours carry quasi-one-dimensional isentropic station solutions [2]_ and a
conical divergence-factor estimate [1]_.

References
----------
.. [1] SUTTON-BIBLARZ-2017. G. P. Sutton and O. Biblarz, "Rocket Propulsion
       Elements", 9th ed., Wiley, 2017.
.. [2] NACA-TR-1135. Ames Research Staff, "Equations, Tables, and Charts for
       Compressible Flow", NACA Report 1135, 1953.
"""

from dataclasses import dataclass, fields
from fractions import Fraction as F
from math import cos, isfinite, pi, radians, tan

from rocket_propulsion.compressible.isentropic import (
    calculate_isentropic,
    mach_from_isentropic,
)
from rocket_propulsion.core.errors import DomainError
from rocket_propulsion.core.validation import integer, numeric_range_guard, real, text

from .dimensional_numerics import (
    dimensional,
    radial_growth,
    radius_exit,
    station_geometry,
    throat_radius,
    wall_area,
)


@dataclass(frozen=True, slots=True)
class NozzleStation:
    """Geometry and ideal-flow properties at one axial location.

    Attributes
    ----------
    x_m : float
        Axial position in m; the throat is at ``x = 0`` for generated contours.
    radius_m, area_m2 : float
        Local wall radius and cross-sectional area.
    area_ratio : float
        Local area divided by the throat area.
    mach : float
        Quasi-one-dimensional isentropic Mach number for ``area_ratio`` on the
        subsonic (converging) or supersonic (diverging) branch [1]_.
    pressure_total_ratio, temperature_total_ratio : float
        Static-to-stagnation ratios ``p/p0`` and ``T/T0`` at ``mach``.

    References
    ----------
    .. [1] NACA-TR-1135. Ames Research Staff, "Equations, Tables, and Charts for
           Compressible Flow", NACA Report 1135, 1953.
    """

    x_m: float
    radius_m: float
    area_m2: float
    area_ratio: float
    mach: float
    pressure_total_ratio: float
    temperature_total_ratio: float


@dataclass(frozen=True, slots=True)
class NozzleContour:
    """Dimensions and sampled stations for an axisymmetric C-D nozzle.

    The divergence-efficiency estimate, when present, uses the conical
    divergence factor ``(1 + cos(alpha))/2`` [1]_.

    Attributes
    ----------
    contour : str
        Contour family, for example ``"conical"``, ``"preliminary"``,
        ``"rao"``, ``"moc-planar"`` or ``"moc-axisymmetric-corrected"``.
    throat_radius_m, chamber_radius_m, exit_radius_m : float
        Characteristic radii in m.
    converging_length_m, diverging_length_m, total_length_m : float
        Section lengths in m.
    throat_area_m2, exit_area_m2 : float
        Throat and exit areas in m^2.
    area_ratio : float
        Exit-to-throat area ratio.
    half_angle_deg : float
        Divergence half-angle (conical), or the representative angle reported
        by the generator for other families.
    stations : tuple of NozzleStation
        Sampled stations ordered by increasing ``x_m``.
    model_fidelity : str, default "preliminary"
        Description of the construction fidelity.
    wall_tangency_residual_deg, centerline_symmetry_residual_deg : float or None
        Characteristic-mesh residuals, when the generator computes them.
    convergence_delta_exit_mach : float or None
        Exit-Mach discrepancy reported by the generator, if any.
    estimated_divergence_efficiency : float or None
        Divergence factor estimate, if available.
    surface_area_m2 : float or None
        Wetted wall area from the sampled stations, if available.
    source : str or None
        Short statement of the construction basis.

    References
    ----------
    .. [1] SUTTON-BIBLARZ-2017. G. P. Sutton and O. Biblarz, "Rocket Propulsion
           Elements", 9th ed., Wiley, 2017.
    """

    contour: str
    throat_radius_m: float
    chamber_radius_m: float
    exit_radius_m: float
    converging_length_m: float
    diverging_length_m: float
    total_length_m: float
    throat_area_m2: float
    exit_area_m2: float
    area_ratio: float
    half_angle_deg: float
    stations: tuple[NozzleStation, ...]
    model_fidelity: str = "preliminary"
    wall_tangency_residual_deg: float | None = None
    centerline_symmetry_residual_deg: float | None = None
    convergence_delta_exit_mach: float | None = None
    estimated_divergence_efficiency: float | None = None
    surface_area_m2: float | None = None
    source: str | None = None


def _positive(value: float, label: str) -> float:
    """Own a finite positive scalar before geometry or flow work."""
    value = real(value, label)
    if value <= 0.0:
        raise DomainError(f"{label} must be finite and greater than zero.", field=label)
    return value


def _checked_contour(contour: NozzleContour) -> NozzleContour:
    """Refuse nonfinite, erased positive or unordered generated output.

    This postflight is for producer-owned records. It does not change the
    separate admission contracts of caller-created NozzleContour records or
    certify a continuous wall, a flow field, or a manufactured nozzle.
    """
    positive = {'throat_radius_m', 'chamber_radius_m', 'exit_radius_m',
                'converging_length_m', 'total_length_m', 'throat_area_m2',
                'exit_area_m2', 'area_ratio', 'surface_area_m2'}
    for item in fields(contour):
        value = getattr(contour, item.name)
        if isinstance(value, (float, int)) and (not isfinite(value) or value < 0.
                                               or (item.name in positive and value == 0.)):
            raise DomainError('Generated contour is outside its finite output range.',
                              code='numeric_range', field=item.name)
    for index, station in enumerate(contour.stations):
        for item in fields(station):
            value = getattr(station, item.name)
            if (not isfinite(value) or (item.name not in ('x_m', 'mach') and value <= 0.)
                    or (item.name == 'mach' and value < 0.)):
                raise DomainError('Generated station is outside its finite output range.',
                                  code='numeric_range', field=item.name)
        if index and station.x_m <= contour.stations[index-1].x_m:
            raise DomainError('Generated axial stations cannot be resolved as increasing floats.',
                              code='numeric_range', field='x_m')
    return contour


def _hermite_radius(
    fraction: float,
    length: float,
    throat_radius: float,
    radius_increment: float,
) -> float:
    """Return a smooth bell-like radius using endpoint slope constraints.

    Parameters
    ----------
    fraction : float
        Normalized position ``s`` in ``[0, 1]`` along the diverging section.
    length : float
        Diverging-section length ``L`` in m.
    throat_radius, radius_increment : float
        Throat radius and exit-minus-throat radius in m, retained separately.

    Returns
    -------
    float
        Radius from the cubic Hermite interpolant with end values
        ``throat_radius`` and ``exit_radius`` and end slopes ``tan(30 deg)``
        and ``tan(8 deg)``, scaled by ``L``. This is a geometric interpolant,
        not an optimized or characteristic-based wall.
    """

    h10 = fraction**3 - 2.0 * fraction**2 + fraction
    h01 = -2.0 * fraction**3 + 3.0 * fraction**2
    h11 = fraction**3 - fraction**2
    inlet_slope = tan(radians(30.0))
    exit_slope = tan(radians(8.0))
    return dimensional(F(throat_radius) + F(h10)*F(length)*F(inlet_slope)
                       + F(h01)*F(radius_increment) + F(h11)*F(length)*F(exit_slope),
                       'radius_m')


@numeric_range_guard('geometry')
def generate_nozzle_contour(
    *,
    throat_area_m2: float,
    area_ratio: float,
    gamma: float = 1.22,
    contour: str = "preliminary",
    half_angle_deg: float = 15.0,
    chamber_radius_ratio: float = 2.5,
    bell_length_fraction: float = 0.8,
    station_count: int = 121,
) -> NozzleContour:
    """Generate an axisymmetric nozzle contour and ideal station solution.

    Parameters
    ----------
    throat_area_m2 : float
        Throat area in m^2, finite and positive.
    area_ratio : float
        Exit-to-throat area ratio, at least one. At one, all stations sample
        the converging section ending at the sonic throat; no diverging span.
    gamma : float, default 1.22
        Ratio of specific heats for the station solution, greater than one.
    contour : str, default "preliminary"
        ``"conical"`` or ``"preliminary"`` (``"bell"`` is accepted as an alias
        of ``"preliminary"``); case-insensitive.
    half_angle_deg : float, default 15.0
        Conical divergence half-angle in degrees, in ``[2, 45)``.
    chamber_radius_ratio : float, default 2.5
        Chamber-to-throat radius ratio, greater than one.
    bell_length_fraction : float, default 0.8
        Preliminary diverging length as a fraction of a 15 degree cone with the
        same area ratio, in ``[0.5, 1.2]``.
    station_count : int, default 121
        Total number of stations, between 41 and 401.

    Returns
    -------
    NozzleContour
        Dimensions, sampled stations, surface area and divergence estimate.

    Raises
    ------
    DomainError
        If an input is outside its domain.

    Notes
    -----
    The converging side is a cosine contraction of length ``1.5 r_t`` from the
    chamber radius to the throat. A conical diverging side uses the selected
    half-angle; the preliminary option uses a cubic-Hermite curve with fixed
    30 degree throat and 8 degree exit slopes. Station Mach numbers follow the
    isentropic area relation [2]_. The divergence efficiency is
    ``(1 + cos(alpha))/2`` [1]_ with ``alpha`` the conical half-angle, or a
    representative 19 degrees for the preliminary contour. The wall area sums
    conical frusta between stations. Use :func:`generate_rao_contour` when an
    explicit Rao/TOP construction is needed. At unit area ratio the entire
    sample budget resolves the cosine contraction ending at the sonic exit;
    reported divergence angle and length are zero and divergence factor one.

    References
    ----------
    .. [1] SUTTON-BIBLARZ-2017. G. P. Sutton and O. Biblarz, "Rocket Propulsion
           Elements", 9th ed., Wiley, 2017.
    .. [2] NACA-TR-1135. Ames Research Staff, "Equations, Tables, and Charts for
           Compressible Flow", NACA Report 1135, 1953.
    """

    throat_area_m2 = _positive(throat_area_m2, 'throat_area_m2')
    area_ratio = _positive(area_ratio, 'area_ratio')
    gamma = _positive(gamma, 'gamma')
    half_angle_deg = _positive(half_angle_deg, 'half_angle_deg')
    chamber_radius_ratio = _positive(chamber_radius_ratio, 'chamber_radius_ratio')
    bell_length_fraction = _positive(bell_length_fraction, 'bell_length_fraction')
    station_count = integer(station_count, 'station_count')
    normalized_contour = text(contour, 'contour').strip().lower()
    if area_ratio < 1.0:
        raise DomainError("Exit-to-throat area ratio must be at least one.", field='area_ratio')
    if not 1.0 < gamma:
        raise DomainError("Gamma must be greater than one.", field='gamma')
    if normalized_contour == "bell":
        normalized_contour = "preliminary"
    if normalized_contour not in {"conical", "preliminary"}:
        raise DomainError("Nozzle contour must be conical or preliminary.", field='contour')
    if not 2.0 <= half_angle_deg < 45.0:
        raise DomainError("Conical half-angle must be in [2, 45) degrees.", field='half_angle_deg')
    if chamber_radius_ratio <= 1.0:
        raise DomainError("Chamber radius ratio must be greater than one.", field='chamber_radius_ratio')
    if not 0.5 <= bell_length_fraction <= 1.2:
        raise DomainError("Bell length fraction must be in [0.5, 1.2].", field='bell_length_fraction')
    if not 41 <= station_count <= 401:
        raise DomainError("Station count must be between 41 and 401.", field='station_count')

    radius_throat = throat_radius(throat_area_m2)
    radius_increment = radial_growth(radius_throat, area_ratio)
    exit_radius = radius_exit(throat_area_m2, area_ratio)
    chamber_radius = chamber_radius_ratio * radius_throat
    converging_length = 1.5 * radius_throat
    reference_cone_length = dimensional(F(radius_increment)/F(tan(radians(15.0))), 'diverging_length_m')
    if normalized_contour == "conical":
        diverging_length = dimensional(F(radius_increment)/F(tan(radians(half_angle_deg))), 'diverging_length_m')
    else:
        diverging_length = bell_length_fraction * reference_cone_length

    converging_count = station_count if area_ratio == 1.0 else max(15, round(station_count * 0.3))
    diverging_count = station_count - converging_count + 1
    stations: list[NozzleStation] = []

    def append_station(x_position: float, radius: float, branch: str, prescribed_ratio=None) -> None:
        station_area, local_area_ratio = station_geometry(radius, radius_throat, throat_area_m2,
                                                         prescribed_ratio)
        if local_area_ratio == 1.0:
            mach = 1.0
        else:
            mach = mach_from_isentropic(f"area_ratio_{branch}", local_area_ratio, gamma)
        ratios = calculate_isentropic(mach, gamma)
        stations.append(
            NozzleStation(
                x_m=x_position,
                radius_m=radius,
                area_m2=station_area,
                area_ratio=local_area_ratio,
                mach=mach,
                pressure_total_ratio=ratios.pressure_total_ratio,
                temperature_total_ratio=ratios.temperature_total_ratio,
            )
        )

    for index in range(converging_count):
        fraction = index / (converging_count - 1)
        x_position = -converging_length * (1.0 - fraction)
        radius = radius_throat + 0.5 * (chamber_radius - radius_throat) * (
            1.0 + cos(pi * fraction)
        )
        append_station(x_position, radius, "subsonic")

    for index in range(1, diverging_count):
        fraction = index / (diverging_count - 1)
        x_position = fraction * diverging_length
        if normalized_contour == "conical":
            radius = dimensional(F(radius_throat)+F(fraction)*F(radius_increment), 'radius_m')
        else:
            radius = _hermite_radius(fraction, diverging_length, radius_throat, radius_increment)
        last = index == diverging_count-1
        append_station(x_position, exit_radius if last else radius, "supersonic",
                       area_ratio if last else None)

    surface_area = wall_area(stations)
    effective_angle = (0.0 if area_ratio == 1.0 else
                       half_angle_deg if normalized_contour == "conical" else 19.0)
    return _checked_contour(NozzleContour(
        contour=normalized_contour,
        throat_radius_m=radius_throat,
        chamber_radius_m=chamber_radius,
        exit_radius_m=exit_radius,
        converging_length_m=converging_length,
        diverging_length_m=diverging_length,
        total_length_m=converging_length + diverging_length,
        throat_area_m2=throat_area_m2,
        exit_area_m2=area_ratio * throat_area_m2,
        area_ratio=area_ratio,
        half_angle_deg=0.0 if area_ratio == 1.0 else half_angle_deg,
        stations=tuple(stations),
        model_fidelity=(
            "analytical cosine contraction; sonic exit" if area_ratio == 1.0 else
            "analytical conical" if normalized_contour == "conical" else "preliminary Hermite"
        ),
        estimated_divergence_efficiency=0.5 * (1.0 + cos(radians(effective_angle))),
        surface_area_m2=surface_area,
        source=(
            "cosine contraction ending at sonic throat; no diverging span"
            if area_ratio == 1.0 else "analytical cone"
            if normalized_contour == "conical"
            else "cubic Hermite preliminary contour; not Rao optimized"
        ),
    ))

