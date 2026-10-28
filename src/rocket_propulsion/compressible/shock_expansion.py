"""Shock-expansion theory for two-dimensional polygonal supersonic sections.

A body made of straight segments in supersonic flow, with deflections small
enough that no shock detaches, sees a sequence of discrete oblique shocks and
Prandtl-Meyer expansions; the surface pressures, and hence lift, drag and
moment, then follow exactly from the oblique-shock and expansion relations
(Anderson section 9.7) [1]_. Wave interactions downstream of the trailing
edge do not affect the surface and are ignored, as in the source.

References
----------
.. [1] ANDERSON-FOA-2017. J. D. Anderson Jr., "Fundamentals of Aerodynamics",
       6th ed., McGraw-Hill, 2017, section 9.7, eqs. (9.43)-(9.49),
       pp. 660-664; Examples 9.6, 9.12, 9.13, 14.1.
"""

from __future__ import annotations

import itertools
from dataclasses import dataclass
from math import atan2, cos, degrees, isfinite, pi, radians, sin, sqrt, tan

from rocket_propulsion.core.errors import DomainError
from rocket_propulsion.core.metadata import WarningMessage

from .isentropic import calculate_isentropic, mach_from_isentropic
from .oblique_shock import calculate_oblique_shock

Point = tuple[float, float]


@dataclass(frozen=True, slots=True)
class AirfoilSection:
    """Polygonal two-dimensional section in body axes [1]_.

    Attributes
    ----------
    upper, lower : tuple of (float, float)
        Upper- and lower-surface vertices ``(x, y)`` from the leading edge to
        the trailing edge; both start at the same leading-edge point and ``x``
        increases strictly along each surface.
    name : str
        Label.

    Raises
    ------
    DomainError
        If a surface has fewer than two points, a coordinate is not finite,
        ``x`` does not increase, the surfaces do not share the leading edge or
        the upper surface lies below the lower one at the trailing edge.

    Notes
    -----
    ``x`` points from the leading edge toward the trailing edge and ``y`` is
    normal to it. The reference chord is the ``x`` extent of the section. If
    the two trailing-edge points differ, the segment joining them is a base
    whose pressure is set by the analysis.

    References
    ----------
    .. [1] ANDERSON-FOA-2017. J. D. Anderson Jr., "Fundamentals of
           Aerodynamics", 6th ed., McGraw-Hill, 2017, Figs. 9.36-9.39.
    """

    upper: tuple[Point, ...]
    lower: tuple[Point, ...]
    name: str = "section"

    def __post_init__(self) -> None:
        upper = tuple((float(x), float(y)) for x, y in self.upper)
        lower = tuple((float(x), float(y)) for x, y in self.lower)
        object.__setattr__(self, "upper", upper)
        object.__setattr__(self, "lower", lower)
        for label, surface in (("upper", upper), ("lower", lower)):
            if len(surface) < 2:
                raise DomainError(f"The {label} surface needs at least two points.", field=label)
            for x, y in surface:
                if not (isfinite(x) and isfinite(y)):
                    raise DomainError(f"The {label} surface has a non-finite coordinate.",
                                      field=label)
            for (x0, _), (x1, _) in itertools.pairwise(surface):
                if x1 <= x0:
                    raise DomainError(f"x must increase strictly along the {label} surface.",
                                      field=label)
        if upper[0] != lower[0]:
            raise DomainError("Upper and lower surfaces must share the leading edge.")
        if upper[-1][1] < lower[-1][1] or upper[-1][0] != lower[-1][0]:
            raise DomainError("Trailing-edge points must share x with upper above lower.")

    @property
    def chord(self) -> float:
        """Reference chord, the ``x`` extent of the section [1]_.

        References
        ----------
        .. [1] ANDERSON-FOA-2017. J. D. Anderson Jr., "Fundamentals of
               Aerodynamics", 6th ed., McGraw-Hill, 2017, Fig. 9.39.
        """
        return self.upper[-1][0] - self.upper[0][0]


@dataclass(frozen=True, slots=True)
class PanelSolution:
    """Flow on one straight panel of a section [1]_.

    Attributes
    ----------
    surface : str
        ``"upper"`` or ``"lower"``.
    index : int
        Panel number from the leading edge.
    start, end : (float, float)
        Panel end points in body axes.
    inclination_deg : float
        Angle by which the panel turns the free stream into itself (positive
        compression, negative expansion), degrees.
    wave : str
        Wave at the panel's leading corner: ``"shock"``, ``"expansion"`` or
        ``"none"``.
    mach : float
        Panel Mach number.
    pressure_ratio : float
        ``p/p_inf``.
    pressure_coefficient : float
        ``(p - p_inf)/q_inf``.

    References
    ----------
    .. [1] ANDERSON-FOA-2017. J. D. Anderson Jr., "Fundamentals of
           Aerodynamics", 6th ed., McGraw-Hill, 2017, section 9.7.
    """

    surface: str
    index: int
    start: Point
    end: Point
    inclination_deg: float
    wave: str
    mach: float
    pressure_ratio: float
    pressure_coefficient: float


@dataclass(frozen=True, slots=True)
class SectionCoefficients:
    """Force and moment coefficients of a section [1]_.

    Attributes
    ----------
    method : str
        Theory used (``"shock-expansion"``, ``"ackeret"`` or ``"busemann"``).
    mach, angle_of_attack_deg, gamma : float
        Free-stream conditions.
    lift_coefficient, drag_coefficient : float
        ``c_l`` and ``c_d`` per unit span, referred to ``q_inf c``.
    normal_force_coefficient, axial_force_coefficient : float
        Body-axis coefficients ``c_n`` and ``c_a``.
    moment_coefficient : float
        Pitching moment about ``moment_reference`` per ``q_inf c^2``,
        positive nose up.
    moment_reference : (float, float)
        Moment reference point in body axes.
    base_pressure_ratio : float
        ``p_base/p_inf`` applied to the base segment, if any.
    panels : tuple of PanelSolution
        Panel-by-panel flow.
    warnings : tuple of WarningMessage
        Validity notes.

    References
    ----------
    .. [1] ANDERSON-FOA-2017. J. D. Anderson Jr., "Fundamentals of
           Aerodynamics", 6th ed., McGraw-Hill, 2017, eqs. (1.15)-(1.19),
           (9.46)-(9.49).
    """

    method: str
    mach: float
    angle_of_attack_deg: float
    gamma: float
    lift_coefficient: float
    drag_coefficient: float
    normal_force_coefficient: float
    axial_force_coefficient: float
    moment_coefficient: float
    moment_reference: Point
    base_pressure_ratio: float
    panels: tuple[PanelSolution, ...]
    warnings: tuple[WarningMessage, ...] = ()


def flat_plate() -> AirfoilSection:
    """Return a unit-chord flat plate (Anderson Fig. 9.36) [1]_.

    References
    ----------
    .. [1] ANDERSON-FOA-2017. J. D. Anderson Jr., "Fundamentals of
           Aerodynamics", 6th ed., McGraw-Hill, 2017, Fig. 9.36.
    """
    return AirfoilSection(((0.0, 0.0), (1.0, 0.0)), ((0.0, 0.0), (1.0, 0.0)), "flat plate")


def diamond_airfoil(half_angle_deg: float, crest_fraction: float = 0.5) -> AirfoilSection:
    """Return a unit-chord symmetric diamond-wedge section (Anderson Fig. 9.37) [1]_.

    Parameters
    ----------
    half_angle_deg : float
        Front-face half-angle ``epsilon`` in degrees, in ``(0, 45)``.
    crest_fraction : float, default 0.5
        Chordwise position of maximum thickness, in ``(0, 1)``.

    Returns
    -------
    AirfoilSection
        Section with thickness ratio ``2 crest_fraction tan(epsilon)``.

    References
    ----------
    .. [1] ANDERSON-FOA-2017. J. D. Anderson Jr., "Fundamentals of
           Aerodynamics", 6th ed., McGraw-Hill, 2017, Fig. 9.37.
    """
    if not (isfinite(half_angle_deg) and 0.0 < half_angle_deg < 45.0):
        raise DomainError("Half-angle must be in (0, 45) degrees.", field="half_angle_deg")
    if not (isfinite(crest_fraction) and 0.0 < crest_fraction < 1.0):
        raise DomainError("Crest fraction must be in (0, 1).", field="crest_fraction")
    half = crest_fraction * tan(radians(half_angle_deg))
    return AirfoilSection(
        ((0.0, 0.0), (crest_fraction, half), (1.0, 0.0)),
        ((0.0, 0.0), (crest_fraction, -half), (1.0, 0.0)),
        "diamond",
    )


def wedge_section(half_angle_deg: float) -> AirfoilSection:
    """Return a unit-chord symmetric wedge with a flat base (Anderson Fig. 9.17) [1]_.

    References
    ----------
    .. [1] ANDERSON-FOA-2017. J. D. Anderson Jr., "Fundamentals of
           Aerodynamics", 6th ed., McGraw-Hill, 2017, Figs. 9.17, 9.39b.
    """
    if not (isfinite(half_angle_deg) and 0.0 < half_angle_deg < 45.0):
        raise DomainError("Half-angle must be in (0, 45) degrees.", field="half_angle_deg")
    half = tan(radians(half_angle_deg))
    return AirfoilSection(((0.0, 0.0), (1.0, half)), ((0.0, 0.0), (1.0, -half)), "wedge")


def parabolic_biconvex(thickness_ratio: float, panels: int = 64) -> AirfoilSection:
    """Return a symmetric parabolic-arc biconvex section, ``y = +-2 t x (1 - x)`` [1]_.

    Parameters
    ----------
    thickness_ratio : float
        Maximum thickness over chord, in ``(0, 0.5)``.
    panels : int, default 64
        Straight panels per surface (the arcs are inscribed polygons).

    References
    ----------
    .. [1] ANDERSON-FOA-2017. J. D. Anderson Jr., "Fundamentals of
           Aerodynamics", 6th ed., McGraw-Hill, 2017, Fig. 12.3 (bi-convex
           section); the parabolic arc is this library's geometric choice.
    """
    if not (isfinite(thickness_ratio) and 0.0 < thickness_ratio < 0.5):
        raise DomainError("Thickness ratio must be in (0, 0.5).", field="thickness_ratio")
    if panels < 2:
        raise DomainError("At least two panels are needed.", field="panels")
    xs = [k / panels for k in range(panels + 1)]
    upper = tuple((x, 2.0 * thickness_ratio * x * (1.0 - x)) for x in xs)
    lower = tuple((x, -y) for x, y in upper)
    return AirfoilSection(upper, lower, "parabolic biconvex")


def _panels(section: AirfoilSection, surface: str) -> list[tuple[Point, Point, float]]:
    """Return ``(start, end, slope angle)`` per panel."""
    points = section.upper if surface == "upper" else section.lower
    return [
        (a, b, atan2(b[1] - a[1], b[0] - a[0]))
        for a, b in itertools.pairwise(points)
    ]


def _inclination(surface: str, slope: float, alpha: float) -> float:
    """Compression angle of a panel relative to the free stream, radians.

    The free stream is ``(cos(alpha), sin(alpha))`` in body axes, so an upper
    panel of slope ``phi`` turns it by ``phi - alpha`` into itself and a lower
    panel by ``alpha - phi``.
    """
    return slope - alpha if surface == "upper" else alpha - slope


def integrate_section_forces(
    section: AirfoilSection,
    panel_pressure_coefficients: dict[str, list[float]],
    alpha_rad: float,
    base_pressure_coefficient: float,
    moment_reference: Point,
) -> tuple[float, float, float, float, float]:
    """Integrate panel pressures into ``(c_l, c_d, c_n, c_a, c_m)`` [1]_.

    Parameters
    ----------
    section : AirfoilSection
        Geometry.
    panel_pressure_coefficients : dict
        ``{"upper": [...], "lower": [...]}`` uniform ``Cp`` per panel.
    alpha_rad : float
        Angle of attack, radians.
    base_pressure_coefficient : float
        ``Cp`` on the base segment (ignored when there is no base).
    moment_reference : (float, float)
        Moment reference point.

    Returns
    -------
    tuple of float
        Lift, drag, normal-force, axial-force and pitching-moment
        coefficients (moment positive nose up).

    Notes
    -----
    Each panel carries the force ``-Cp q n L`` with ``n`` the outward normal,
    acting at the panel midpoint (uniform pressure). Lift and drag are the
    components normal and parallel to the free stream, Anderson eqs.
    (1.18)-(1.19) [1]_; the closed contour makes the free-stream pressure
    integrate to zero, so gauge pressure is used throughout.

    References
    ----------
    .. [1] ANDERSON-FOA-2017. J. D. Anderson Jr., "Fundamentals of
           Aerodynamics", 6th ed., McGraw-Hill, 2017, eqs. (1.15)-(1.19).
    """

    fx = fy = moment_z = 0.0
    x0, y0 = moment_reference
    segments: list[tuple[Point, Point, float, float, float]] = []
    for surface in ("upper", "lower"):
        for (a, b, _), cp in zip(_panels(section, surface),
                                 panel_pressure_coefficients[surface], strict=True):
            dx, dy = b[0] - a[0], b[1] - a[1]
            if surface == "upper":
                nx, ny = -dy, dx
            else:
                nx, ny = dy, -dx
            segments.append((a, b, nx, ny, cp))
    upper_te, lower_te = section.upper[-1], section.lower[-1]
    if upper_te != lower_te:
        dx, dy = lower_te[0] - upper_te[0], lower_te[1] - upper_te[1]
        segments.append((upper_te, lower_te, -dy, dx, base_pressure_coefficient))
    for a, b, nx, ny, cp in segments:
        force_x, force_y = -cp * nx, -cp * ny
        mx, my = 0.5 * (a[0] + b[0]) - x0, 0.5 * (a[1] + b[1]) - y0
        fx += force_x
        fy += force_y
        moment_z += mx * force_y - my * force_x
    chord = section.chord
    ca, cn = fx / chord, fy / chord
    lift = -ca * sin(alpha_rad) + cn * cos(alpha_rad)
    drag = ca * cos(alpha_rad) + cn * sin(alpha_rad)
    return lift, drag, cn, ca, -moment_z / (chord * chord)


def _expand(mach: float, turn: float, gamma: float) -> float:
    """Mach number after a Prandtl-Meyer expansion through ``turn`` radians."""
    nu = calculate_isentropic(mach, gamma).prandtl_meyer_deg or 0.0
    nu_max = degrees(0.5 * pi * (sqrt((gamma + 1.0) / (gamma - 1.0)) - 1.0))
    target = nu + degrees(turn)
    if target >= nu_max:
        raise DomainError(
            "Expansion exceeds the Prandtl-Meyer limit; the flow would separate into vacuum.",
            code="expansion_limit",
        )
    return mach_from_isentropic("prandtl_meyer_deg", target, gamma)


def shock_expansion_coefficients(
    section: AirfoilSection,
    mach: float,
    angle_of_attack_deg: float,
    *,
    gamma: float = 1.4,
    base_pressure_ratio: float = 1.0,
    moment_reference: Point | None = None,
) -> SectionCoefficients:
    """Exact inviscid section coefficients by shock-expansion theory [1]_.

    Parameters
    ----------
    section : AirfoilSection
        Polygonal geometry.
    mach : float
        Free-stream Mach number, greater than 1.
    angle_of_attack_deg : float
        Angle between the free stream and the body ``x`` axis, positive nose
        up, degrees.
    gamma : float, default 1.4
        Ratio of specific heats.
    base_pressure_ratio : float, default 1.0
        ``p_base/p_inf`` on the base segment of a blunt-based section
        (Anderson Example 9.6 assumes 1). Anderson Example 9.13b neglects
        the base force altogether, which is ``0.0`` here.
    moment_reference : (float, float), optional
        Moment reference point; defaults to the leading edge.

    Returns
    -------
    SectionCoefficients
        Lift, drag and moment coefficients plus the panel solution.

    Raises
    ------
    DomainError
        If an input is invalid, a compression corner detaches the shock, or
        an expansion exceeds the Prandtl-Meyer limit.

    Notes
    -----
    Marching along each surface from the leading edge, the turn of every
    corner relative to the upstream panel is computed; a compression passes
    through an attached weak oblique shock (theta-beta-M relation and the
    normal-shock jump of the normal Mach component) and an expansion through
    a Prandtl-Meyer fan, ``nu_2 = nu_1 + theta`` (eq. (9.43)), with
    ``p0`` constant [1]_. The surface pressures are integrated by
    :func:`integrate_section_forces`.

    References
    ----------
    .. [1] ANDERSON-FOA-2017. J. D. Anderson Jr., "Fundamentals of
           Aerodynamics", 6th ed., McGraw-Hill, 2017, section 9.7, pp. 660-668.
    """

    if not isfinite(mach) or mach <= 1.0:
        raise DomainError("Free-stream Mach number must be greater than 1.", field="mach")
    if not isfinite(gamma) or gamma <= 1.0:
        raise DomainError("Gamma must be finite and greater than 1.", field="gamma")
    if not isfinite(angle_of_attack_deg) or abs(angle_of_attack_deg) >= 90.0:
        raise DomainError("Angle of attack must be finite and below 90 degrees.",
                          field="angle_of_attack_deg")
    if not isfinite(base_pressure_ratio) or base_pressure_ratio < 0.0:
        raise DomainError("Base pressure ratio must be finite and non-negative.",
                          field="base_pressure_ratio")
    alpha = radians(angle_of_attack_deg)
    reference = section.upper[0] if moment_reference is None else moment_reference
    dynamic = 0.5 * gamma * mach * mach
    panels: list[PanelSolution] = []
    coefficients: dict[str, list[float]] = {"upper": [], "lower": []}
    for surface in ("upper", "lower"):
        local_mach = mach
        pressure_ratio = 1.0
        previous = 0.0
        for index, (a, b, slope) in enumerate(_panels(section, surface)):
            inclination = _inclination(surface, slope, alpha)
            turn = inclination - previous
            if turn > 1.0e-12:
                shock = calculate_oblique_shock(
                    upstream_mach=local_mach, turning_angle_deg=degrees(turn), gamma=gamma
                )
                pressure_ratio *= shock.pressure_ratio
                local_mach = shock.downstream_mach
                wave = "shock"
            elif turn < -1.0e-12:
                total_before = calculate_isentropic(local_mach, gamma).pressure_total_ratio
                local_mach = _expand(local_mach, -turn, gamma)
                total_after = calculate_isentropic(local_mach, gamma).pressure_total_ratio
                pressure_ratio *= total_after / total_before
                wave = "expansion"
            else:
                wave = "none"
            cp = (pressure_ratio - 1.0) / dynamic
            coefficients[surface].append(cp)
            panels.append(PanelSolution(surface, index, a, b, degrees(inclination), wave,
                                        local_mach, pressure_ratio, cp))
            previous = inclination
    lift, drag, cn, ca, cm = integrate_section_forces(
        section, coefficients, alpha, (base_pressure_ratio - 1.0) / dynamic, reference
    )
    return SectionCoefficients(
        method="shock-expansion",
        mach=mach,
        angle_of_attack_deg=angle_of_attack_deg,
        gamma=gamma,
        lift_coefficient=lift,
        drag_coefficient=drag,
        normal_force_coefficient=cn,
        axial_force_coefficient=ca,
        moment_coefficient=cm,
        moment_reference=reference,
        base_pressure_ratio=base_pressure_ratio,
        panels=tuple(panels),
    )
