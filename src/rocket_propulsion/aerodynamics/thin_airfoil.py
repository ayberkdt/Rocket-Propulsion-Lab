"""Classical thin airfoil theory for symmetric and cambered sections.

The vortex sheet on the chord line is expanded in Glauert's series, Anderson
Eq. (4.43); its coefficients follow from the camber-line slope by the
Fourier integrals (4.50)-(4.51) with ``x = c (1 - cos theta) / 2`` (4.21).
The aerodynamic coefficients are

* ``c_l = pi (2 A0 + A1)`` (4.56), lift slope ``2 pi`` per radian (4.58);
* ``alpha_L=0 = -(1/pi) int_0^pi (dz/dx)(cos theta - 1) d theta`` (4.61);
* ``c_m,le = -(pi/2)(A0 + A1 - A2/2)`` (4.62);
* ``c_m,c/4 = (pi/4)(A2 - A1)`` (4.64), independent of incidence, so the
  quarter chord is the theoretical aerodynamic centre;
* ``x_cp / c = (1/4)(1 + pi (A1 - A2) / c_l)`` (4.66).

The symmetric section is the special case ``dz/dx = 0`` (Eqs. 4.33, 4.39,
4.41). The ideal angle of attack, at which the leading-edge term ``A0`` of
(4.43) vanishes, is ``alpha_i = (1/pi) int_0^pi (dz/dx) d theta`` from (4.50),
with ideal lift coefficient ``pi A1`` from (4.56); NACA Report 824 tabulates
the same pair for its mean lines and relates them to the zero-lift angle by
its Eq. (16).

Angles are in radians; coefficients are per unit span on the chord.

References
----------
.. [1] ANDERSON-FOA-2017. J. D. Anderson Jr., Fundamentals of Aerodynamics,
       6th ed., McGraw-Hill, 2017, Sections 4.7-4.9, Eqs. (4.21)-(4.66).
.. [2] NACA-TR-824. I. H. Abbott, A. E. von Doenhoff, L. S. Stivers Jr.,
       "Summary of Airfoil Data", NACA Report 824, 1945, p. 14, Eq. (16).
"""

from __future__ import annotations

from collections.abc import Callable, Sequence
from dataclasses import dataclass
from math import acos, cos, isfinite, pi

from rocket_propulsion.core.errors import DomainError
from rocket_propulsion.core.metadata import WarningMessage

from ._subsonic_numerics import integrate
from .airfoil_geometry import NacaSection

# Anderson (Section 4.10 introduction) states that thin airfoil results compare
# favourably with experiment for sections of about 12 percent thickness or less.
_THIN_SECTION_LIMIT = 0.12
# Library threshold for the small-angle assumption of Eq. (4.14); not a source value.
_SMALL_ANGLE_LIMIT_RAD = 0.15


@dataclass(frozen=True, slots=True)
class ThinAirfoilResult:
    """Thin-airfoil aerodynamic coefficients of one section at one incidence [1]_ [2]_.

    Attributes
    ----------
    angle_of_attack_rad : float
        Geometric incidence of the chord line.
    fourier_a0 : float
        ``A0`` of Eq. (4.50).
    fourier_coefficients : tuple of float
        ``A1 ... AN`` of Eq. (4.51).
    lift_coefficient : float
        ``c_l`` (4.56).
    lift_slope_per_rad : float
        ``d c_l / d alpha = 2 pi`` (4.58).
    zero_lift_angle_rad : float
        ``alpha_L=0`` (4.61).
    moment_coefficient_leading_edge : float
        ``c_m,le`` (4.62), nose-up positive.
    moment_coefficient_quarter_chord : float
        ``c_m,c/4`` (4.64).
    center_of_pressure : float or None
        ``x_cp / c`` (4.66); ``None`` at zero lift, where it is undefined.
    ideal_angle_of_attack_rad : float
        Incidence with ``A0 = 0``.
    ideal_lift_coefficient : float
        ``pi A1``.
    warnings : tuple of WarningMessage
        Validity warnings (thickness, incidence).

    References
    ----------
    .. [1] ANDERSON-FOA-2017. J. D. Anderson Jr., Fundamentals of
           Aerodynamics, 6th ed., McGraw-Hill, 2017, Sections 4.7-4.8.
    .. [2] NACA-TR-824. I. H. Abbott, A. E. von Doenhoff, L. S. Stivers Jr.,
           "Summary of Airfoil Data", NACA Report 824, 1945, p. 14.
    """

    angle_of_attack_rad: float
    fourier_a0: float
    fourier_coefficients: tuple[float, ...]
    lift_coefficient: float
    lift_slope_per_rad: float
    zero_lift_angle_rad: float
    moment_coefficient_leading_edge: float
    moment_coefficient_quarter_chord: float
    center_of_pressure: float | None
    ideal_angle_of_attack_rad: float
    ideal_lift_coefficient: float
    warnings: tuple[WarningMessage, ...] = ()


def thin_airfoil(
    camber_slope: Callable[[float], float],
    angle_of_attack_rad: float,
    *,
    slope_breakpoints: Sequence[float] = (),
    coefficient_count: int = 4,
) -> ThinAirfoilResult:
    """Apply thin airfoil theory to a mean line given by its slope ``dz/dx``.

    Parameters
    ----------
    camber_slope : callable
        ``dz/dx`` as a function of ``x/c`` on ``[0, 1]``; ``lambda x: 0.0``
        for a symmetric section or flat plate.
    angle_of_attack_rad : float
        Incidence of the chord line, finite.
    slope_breakpoints : sequence of float, optional
        Stations ``x/c`` in ``(0, 1)`` where the slope or its derivatives are
        discontinuous (e.g. the four-digit ``p``); the quadrature is split
        there.
    coefficient_count : int, default 4
        Number of coefficients ``A1 ... AN`` returned, at least 2.

    Returns
    -------
    ThinAirfoilResult
        Coefficients of Eqs. (4.50)-(4.66).

    Raises
    ------
    DomainError
        If an input is not finite, a breakpoint is outside ``(0, 1)``, the
        slope returns a non-finite value, or fewer than two coefficients are
        requested.

    Notes
    -----
    The integrals over ``theta`` in ``[0, pi]`` are evaluated with composite
    Gauss-Legendre quadrature split at the breakpoints, so the camber-slope
    function is never evaluated exactly at a leading-edge or breakpoint
    singularity. A warning is issued for ``|alpha| > 0.15 rad``, a library
    threshold for the small-angle assumption of Eq. (4.14) [1]_.

    References
    ----------
    .. [1] ANDERSON-FOA-2017. J. D. Anderson Jr., Fundamentals of
           Aerodynamics, 6th ed., McGraw-Hill, 2017, Sections 4.7-4.8,
           Eqs. (4.14), (4.21), (4.43), (4.50), (4.51), (4.56), (4.58),
           (4.61), (4.62), (4.64), (4.66).
    .. [2] NACA-TR-824. I. H. Abbott, A. E. von Doenhoff, L. S. Stivers Jr.,
           "Summary of Airfoil Data", NACA Report 824, 1945, p. 14, Eq. (16).
    """

    alpha = float(angle_of_attack_rad)
    if not isfinite(alpha):
        raise DomainError("Angle of attack must be finite.", field="angle_of_attack_rad")
    if coefficient_count < 2:
        raise DomainError("At least two Fourier coefficients are required.", field="coefficient_count")
    stations = sorted(float(value) for value in slope_breakpoints)
    if any(not isfinite(value) or value <= 0.0 or value >= 1.0 for value in stations):
        raise DomainError("Slope breakpoints must lie strictly between 0 and 1.", field="slope_breakpoints")
    limits = [0.0, *(acos(1.0 - 2.0 * value) for value in stations), pi]

    def slope(theta: float) -> float:
        value = float(camber_slope(0.5 * (1.0 - cos(theta))))
        if not isfinite(value):
            raise DomainError("Camber slope must be finite inside the chord.", field="camber_slope")
        return value

    mean_slope = integrate(slope, limits) / pi
    a0 = alpha - mean_slope
    coefficients = tuple(
        2.0 / pi * integrate(lambda t, n=n: slope(t) * cos(n * t), limits)
        for n in range(1, max(coefficient_count, 2) + 1)
    )
    a1, a2 = coefficients[0], coefficients[1]
    zero_lift = -integrate(lambda t: slope(t) * (cos(t) - 1.0), limits) / pi
    lift = pi * (2.0 * a0 + a1)
    moment_le = -0.5 * pi * (a0 + a1 - 0.5 * a2)
    moment_quarter = 0.25 * pi * (a2 - a1)
    center = None if lift == 0.0 else 0.25 * (1.0 + pi * (a1 - a2) / lift)
    warnings: list[WarningMessage] = []
    if abs(alpha) > _SMALL_ANGLE_LIMIT_RAD:
        warnings.append(WarningMessage(
            "thin_airfoil_large_incidence",
            "Thin airfoil theory assumes small angles of attack (Anderson Eq. 4.14).",
            "warning",
            "angle_of_attack_rad",
        ))
    return ThinAirfoilResult(
        angle_of_attack_rad=alpha,
        fourier_a0=a0,
        fourier_coefficients=coefficients[:coefficient_count],
        lift_coefficient=lift,
        lift_slope_per_rad=2.0 * pi,
        zero_lift_angle_rad=zero_lift,
        moment_coefficient_leading_edge=moment_le,
        moment_coefficient_quarter_chord=moment_quarter,
        center_of_pressure=center,
        ideal_angle_of_attack_rad=mean_slope,
        ideal_lift_coefficient=pi * a1,
        warnings=tuple(warnings),
    )


def thin_airfoil_naca(section: NacaSection, angle_of_attack_rad: float) -> ThinAirfoilResult:
    """Thin airfoil theory for a NACA four- or five-digit mean line.

    Parameters
    ----------
    section : NacaSection
        Section from :func:`~rocket_propulsion.aerodynamics.airfoil_geometry.naca_section`.
    angle_of_attack_rad : float
        Incidence of the chord line.

    Returns
    -------
    ThinAirfoilResult
        As :func:`thin_airfoil`; a ``thin_airfoil_thick_section`` warning is
        added when the thickness exceeds 12 percent of chord, the range in
        which Anderson reports favourable agreement with experiment [1]_.

    Notes
    -----
    The mean-line slope is that of the NACA definitions (TR-460, TR-537)
    [2]_ [3]_; the quadrature is split at the arc junction.

    References
    ----------
    .. [1] ANDERSON-FOA-2017. J. D. Anderson Jr., Fundamentals of
           Aerodynamics, 6th ed., McGraw-Hill, 2017, Sections 4.8 and 4.10.
    .. [2] NACA-TR-460. E. N. Jacobs, K. E. Ward, R. M. Pinkerton, "The
           Characteristics of 78 Related Airfoil Sections from Tests in the
           Variable-Density Wind Tunnel", NACA Report 460, 1933, p. 5.
    .. [3] NACA-TR-537. E. N. Jacobs, R. M. Pinkerton, "Tests in the
           Variable-Density Wind Tunnel of Related Airfoils Having the Maximum
           Camber Unusually Far Forward", NACA Report 537, 1936, pp. 521-522.
    """

    if not isinstance(section, NacaSection):
        raise DomainError("A NacaSection is required.", field="section")
    breakpoints = () if section.is_symmetric else (section.junction,)
    result = thin_airfoil(
        lambda x: section.mean_line(x)[1],
        angle_of_attack_rad,
        slope_breakpoints=breakpoints,
    )
    if section.thickness_ratio > _THIN_SECTION_LIMIT:
        warning = WarningMessage(
            "thin_airfoil_thick_section",
            "Thin airfoil theory compares favourably with experiment for sections of about "
            "12 percent thickness or less (Anderson Section 4.10).",
            "warning",
            "section",
        )
        result = ThinAirfoilResult(
            **{name: getattr(result, name) for name in ThinAirfoilResult.__slots__ if name != "warnings"},
            warnings=(*result.warnings, warning),
        )
    return result
