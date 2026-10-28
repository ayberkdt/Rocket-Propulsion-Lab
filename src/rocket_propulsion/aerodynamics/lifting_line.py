"""Prandtl lifting-line theory solved with Glauert's Fourier series.

The spanwise circulation is expanded as
``Gamma(theta) = 2 b V sum_n A_n sin(n theta)`` with
``y = -(b/2) cos(theta)`` (Anderson Eqs. 5.46, 5.48). Substituting in the
fundamental equation (5.23) gives at each spanwise station ``theta_0``

    alpha(theta_0) = 4 b / (a_0 c(theta_0)) sum A_n sin(n theta_0)
                     + alpha_L=0(theta_0) + sum n A_n sin(n theta_0) / sin(theta_0)

which is Anderson Eq. (5.51) for the thin-airfoil section slope
``a_0 = 2 pi`` and NACA Report 572 Eq. (5) for a general section slope
``m_0`` (from ``c_l = a_0 (alpha_eff - alpha_L=0)``, Eq. 5.19, and
``c_l = 2 Gamma / (V c)``, Eq. 5.21). Collocating at ``N`` stations gives a
linear system for ``A_1 ... A_N``. Then

* ``C_L = pi AR A_1`` (5.53);
* ``C_D,i = pi AR sum n A_n^2 = C_L^2 (1 + delta) / (pi AR)`` with
  ``delta = sum_{n>=2} n (A_n / A_1)^2`` (5.60)-(5.61);
* span efficiency ``e = 1 / (1 + delta)`` (5.62);
* induced angle ``alpha_i(theta) = sum n A_n sin(n theta) / sin(theta)``
  (5.57);
* lift-slope factor ``tau`` from ``a = a_0 / (1 + (a_0 / (pi AR))(1 + tau))``
  (5.70); for the elliptic wing ``delta = tau = 0`` and
  ``a = a_0 / (1 + a_0 / (pi AR))`` (5.69), ``C_D,i = C_L^2 / (pi AR)``
  (5.43).

Angles in radians, lengths in metres.

References
----------
.. [1] ANDERSON-FOA-2017. J. D. Anderson Jr., Fundamentals of Aerodynamics,
       6th ed., McGraw-Hill, 2017, Section 5.3, Eqs. (5.19)-(5.70).
.. [2] NACA-TR-572. R. F. Anderson, "Determination of the Characteristics of
       Tapered Wings", NACA Report 572, 1937, pp. 3-5, Eq. (5).
"""

from __future__ import annotations

from collections.abc import Callable
from dataclasses import dataclass
from math import acos, cos, isfinite, pi, sin, sqrt

from rocket_propulsion.core.errors import DomainError
from rocket_propulsion.core.validation import numeric_range_guard

from ._subsonic_numerics import integrate, solve_dense

ChordFunction = Callable[[float], float]


@dataclass(frozen=True, slots=True)
class Wing:
    """Straight (unswept) finite wing described by spanwise distributions [1]_ [2]_.

    Attributes
    ----------
    span : float
        Tip-to-tip span ``b`` in m.
    chord : callable
        ``c(y)`` in m for ``-b/2 <= y <= b/2``; positive inside the span.
    twist : callable
        Geometric incidence of the local chord relative to the root chord,
        ``alpha(y) - alpha_root`` in rad (negative for washout).
    zero_lift_angle : callable
        Section ``alpha_L=0(y)`` in rad.
    section_lift_slope : callable
        Section lift slope ``a_0(y)`` per rad (``2 pi`` from thin airfoil
        theory).
    area : float
        Planform area ``S = int c dy`` in m^2.
    breakpoints : tuple of float
        Interior spanwise stations ``y`` where ``c(y)`` has a kink (for the
        area quadrature).

    Notes
    -----
    Build with :func:`trapezoidal_wing`, :func:`elliptic_wing` or
    :func:`custom_wing`.

    References
    ----------
    .. [1] ANDERSON-FOA-2017. J. D. Anderson Jr., Fundamentals of
           Aerodynamics, 6th ed., McGraw-Hill, 2017, Sections 5.3.1-5.3.2.
    .. [2] NACA-TR-572. R. F. Anderson, "Determination of the
           Characteristics of Tapered Wings", NACA Report 572, 1937, p. 5.
    """

    span: float
    chord: ChordFunction
    twist: ChordFunction
    zero_lift_angle: ChordFunction
    section_lift_slope: ChordFunction
    area: float
    breakpoints: tuple[float, ...] = ()

    @property
    @numeric_range_guard("span")
    def aspect_ratio(self) -> float:
        """``AR = b^2 / S``."""

        return self.span**2 / self.area


def _positive(value: float, label: str) -> float:
    number = float(value)
    if not isfinite(number) or number <= 0.0:
        raise DomainError(f"{label} must be finite and positive.", field=label)
    return number


def _finite(value: float, label: str) -> float:
    number = float(value)
    if not isfinite(number):
        raise DomainError(f"{label} must be finite.", field=label)
    return number


def trapezoidal_wing(
    span: float,
    root_chord: float,
    taper_ratio: float,
    *,
    tip_twist_rad: float = 0.0,
    zero_lift_angle_rad: float = 0.0,
    section_lift_slope: float = 2.0 * pi,
) -> Wing:
    """Straight-tapered wing with linear twist and uniform section properties.

    Parameters
    ----------
    span : float
        ``b`` in m.
    root_chord : float
        ``c_r`` in m.
    taper_ratio : float
        ``c_t / c_r`` in ``[0, 1]`` (Anderson Fig. 5.20 abscissa).
    tip_twist_rad : float, default 0
        Geometric twist at the tip relative to the root, varying linearly
        with ``|y|`` (negative = washout).
    zero_lift_angle_rad : float, default 0
        Section zero-lift angle, uniform.
    section_lift_slope : float, default ``2 pi``
        Section lift slope per radian, uniform.

    Returns
    -------
    Wing
        ``S = b c_r (1 + lambda) / 2``.

    References
    ----------
    .. [1] ANDERSON-FOA-2017. J. D. Anderson Jr., Fundamentals of
           Aerodynamics, 6th ed., McGraw-Hill, 2017, Section 5.3.2,
           Fig. 5.18.
    """

    b = _positive(span, "span")
    c_r = _positive(root_chord, "root_chord")
    taper = _finite(taper_ratio, "taper_ratio")
    if not 0.0 <= taper <= 1.0:
        raise DomainError("Taper ratio must lie in [0, 1].", field="taper_ratio")
    twist = _finite(tip_twist_rad, "tip_twist_rad")
    zero_lift = _finite(zero_lift_angle_rad, "zero_lift_angle_rad")
    slope = _positive(section_lift_slope, "section_lift_slope")
    half = 0.5 * b
    return Wing(
        span=b,
        chord=lambda y: c_r * (1.0 - (1.0 - taper) * abs(y) / half),
        twist=lambda y: twist * abs(y) / half,
        zero_lift_angle=lambda y: zero_lift,
        section_lift_slope=lambda y: slope,
        area=half * c_r * (1.0 + taper),
        breakpoints=(0.0,),
    )


def elliptic_wing(
    span: float,
    root_chord: float,
    *,
    zero_lift_angle_rad: float = 0.0,
    section_lift_slope: float = 2.0 * pi,
) -> Wing:
    """Untwisted elliptic-planform wing ``c(y) = c_r sqrt(1 - (2y/b)^2)``.

    Parameters
    ----------
    span : float
        ``b`` in m.
    root_chord : float
        ``c_r`` in m.
    zero_lift_angle_rad : float, default 0
        Uniform section zero-lift angle.
    section_lift_slope : float, default ``2 pi``
        Uniform section lift slope per radian.

    Returns
    -------
    Wing
        ``S = pi b c_r / 4``.

    Notes
    -----
    Without twist this planform carries the elliptic lift distribution
    (Anderson Section 5.3.1, Eq. 5.45) [1]_.

    References
    ----------
    .. [1] ANDERSON-FOA-2017. J. D. Anderson Jr., Fundamentals of
           Aerodynamics, 6th ed., McGraw-Hill, 2017, Section 5.3.1,
           Eqs. (5.31)-(5.45).
    """

    b = _positive(span, "span")
    c_r = _positive(root_chord, "root_chord")
    zero_lift = _finite(zero_lift_angle_rad, "zero_lift_angle_rad")
    slope = _positive(section_lift_slope, "section_lift_slope")
    half = 0.5 * b
    return Wing(
        span=b,
        chord=lambda y: c_r * sqrt(max(0.0, 1.0 - (y / half) ** 2)),
        twist=lambda y: 0.0,
        zero_lift_angle=lambda y: zero_lift,
        section_lift_slope=lambda y: slope,
        area=0.25 * pi * b * c_r,
    )


def custom_wing(
    span: float,
    chord: ChordFunction,
    *,
    twist: ChordFunction | None = None,
    zero_lift_angle: ChordFunction | None = None,
    section_lift_slope: ChordFunction | None = None,
    breakpoints: tuple[float, ...] = (),
) -> Wing:
    """Wing from user spanwise distributions; the area is integrated numerically.

    Parameters
    ----------
    span : float
        ``b`` in m.
    chord : callable
        ``c(y)`` in m, finite and positive for ``|y| < b/2``.
    twist, zero_lift_angle : callable, optional
        Spanwise twist and section zero-lift angle in rad (default zero).
    section_lift_slope : callable, optional
        Section lift slope per rad (default ``2 pi``).
    breakpoints : tuple of float, optional
        Interior stations where ``c(y)`` has a kink.

    Returns
    -------
    Wing
        Wing with ``S`` from composite Gauss-Legendre quadrature in
        ``theta`` (``y = -(b/2) cos(theta)``, Eq. 5.46).

    References
    ----------
    .. [1] ANDERSON-FOA-2017. J. D. Anderson Jr., Fundamentals of
           Aerodynamics, 6th ed., McGraw-Hill, 2017, Section 5.3.2.
    """

    b = _positive(span, "span")
    half = 0.5 * b
    stations = tuple(sorted(float(y) for y in breakpoints))
    if any(not isfinite(y) or abs(y) >= half for y in stations):
        raise DomainError("Breakpoints must lie strictly inside the span.", field="breakpoints")
    limits = sorted([0.0, *(acos(-y / half) for y in stations), pi])
    area = integrate(lambda t: float(chord(-half * cos(t))) * half * sin(t), limits)
    return Wing(
        span=b,
        chord=chord,
        twist=twist or (lambda y: 0.0),
        zero_lift_angle=zero_lift_angle or (lambda y: 0.0),
        section_lift_slope=section_lift_slope or (lambda y: 2.0 * pi),
        area=_positive(area, "area"),
        breakpoints=stations,
    )


@dataclass(frozen=True, slots=True)
class LiftingLineResult:
    """Lifting-line solution for one wing at one root incidence [1]_ [2]_.

    Attributes
    ----------
    root_angle_of_attack_rad : float
        Geometric incidence of the root chord.
    coefficients : tuple of float
        Glauert coefficients ``A_1 ... A_N`` (Eq. 5.48).
    lift_coefficient : float
        ``C_L = pi AR A_1`` (5.53).
    induced_drag_coefficient : float
        ``C_D,i = pi AR sum n A_n^2`` (5.60).
    induced_drag_factor : float
        ``delta = sum_{n>=2} n (A_n / A_1)^2`` (5.61); ``nan`` at zero lift.
    span_efficiency : float
        ``e = 1 / (1 + delta)`` (5.62); ``nan`` at zero lift.
    aspect_ratio : float
        ``b^2 / S``.
    stations_y : tuple of float
        Collocation stations in m.
    circulation : tuple of float
        ``Gamma(y) / V_inf`` at the stations in m (Eq. 5.48).
    section_lift_coefficients : tuple of float
        ``c_l(y) = 2 Gamma / (V c)`` (Eq. 5.21).
    induced_angles_rad : tuple of float
        ``alpha_i(y)`` (Eq. 5.57).

    References
    ----------
    .. [1] ANDERSON-FOA-2017. J. D. Anderson Jr., Fundamentals of
           Aerodynamics, 6th ed., McGraw-Hill, 2017, Section 5.3.2.
    .. [2] NACA-TR-572. R. F. Anderson, "Determination of the
           Characteristics of Tapered Wings", NACA Report 572, 1937, pp. 3-5.
    """

    root_angle_of_attack_rad: float
    coefficients: tuple[float, ...]
    lift_coefficient: float
    induced_drag_coefficient: float
    induced_drag_factor: float
    span_efficiency: float
    aspect_ratio: float
    stations_y: tuple[float, ...]
    circulation: tuple[float, ...]
    section_lift_coefficients: tuple[float, ...]
    induced_angles_rad: tuple[float, ...]


def solve_lifting_line(wing: Wing, root_angle_of_attack_rad: float, *, terms: int = 40) -> LiftingLineResult:
    """Solve Prandtl's lifting-line equation for a wing with Glauert's series.

    Parameters
    ----------
    wing : Wing
        Planform and section data.
    root_angle_of_attack_rad : float
        Geometric angle of attack of the root chord.
    terms : int, default 40
        Number ``N`` of series terms and collocation stations
        ``theta_k = k pi / (N + 1)``, ``k = 1 ... N``; at least 3.

    Returns
    -------
    LiftingLineResult
        Coefficients and spanwise distributions.

    Raises
    ------
    DomainError
        If the chord or slope is not finite and positive at a station, an
        input is not finite, or ``terms < 3``.

    Notes
    -----
    Collocation of Eq. (5.51) in its general-slope form (NACA TR-572 Eq. 5);
    all ``N`` sine terms are kept so asymmetric twist is admissible, and the
    even coefficients vanish for symmetric wings.

    References
    ----------
    .. [1] ANDERSON-FOA-2017. J. D. Anderson Jr., Fundamentals of
           Aerodynamics, 6th ed., McGraw-Hill, 2017, Section 5.3.2,
           Eqs. (5.46)-(5.62).
    .. [2] NACA-TR-572. R. F. Anderson, "Determination of the
           Characteristics of Tapered Wings", NACA Report 572, 1937, p. 5,
           Eq. (5).
    """

    if not isinstance(wing, Wing):
        raise DomainError("A Wing is required.", field="wing")
    alpha = _finite(root_angle_of_attack_rad, "root_angle_of_attack_rad")
    if terms < 3:
        raise DomainError("At least three series terms are required.", field="terms")
    half = 0.5 * wing.span
    thetas = [k * pi / (terms + 1) for k in range(1, terms + 1)]
    ys = [-half * cos(theta) for theta in thetas]
    matrix, rhs, chords = [], [], []
    for theta, y in zip(thetas, ys, strict=True):
        chord = _positive(wing.chord(y), "chord")
        slope = _positive(wing.section_lift_slope(y), "section_lift_slope")
        local = 4.0 * wing.span / (slope * chord)
        matrix.append([
            sin(n * theta) * (local + n / sin(theta)) for n in range(1, terms + 1)
        ])
        rhs.append(alpha + _finite(wing.twist(y), "twist") - _finite(wing.zero_lift_angle(y), "zero_lift_angle"))
        chords.append(chord)
    coefficients = solve_dense(matrix, rhs)
    aspect_ratio = wing.aspect_ratio
    a1 = coefficients[0]
    lift = pi * aspect_ratio * a1
    drag = pi * aspect_ratio * sum(n * a * a for n, a in enumerate(coefficients, start=1))
    if a1 != 0.0:
        delta = sum(n * (a / a1) ** 2 for n, a in enumerate(coefficients, start=1) if n >= 2)
        efficiency = 1.0 / (1.0 + delta)
    else:
        delta = efficiency = float("nan")
    circulation = tuple(
        2.0 * wing.span * sum(a * sin(n * theta) for n, a in enumerate(coefficients, start=1))
        for theta in thetas
    )
    return LiftingLineResult(
        root_angle_of_attack_rad=alpha,
        coefficients=tuple(coefficients),
        lift_coefficient=lift,
        induced_drag_coefficient=drag,
        induced_drag_factor=delta,
        span_efficiency=efficiency,
        aspect_ratio=aspect_ratio,
        stations_y=tuple(ys),
        circulation=circulation,
        section_lift_coefficients=tuple(2.0 * g / c for g, c in zip(circulation, chords, strict=True)),
        induced_angles_rad=tuple(
            sum(n * a * sin(n * theta) for n, a in enumerate(coefficients, start=1)) / sin(theta)
            for theta in thetas
        ),
    )


@dataclass(frozen=True, slots=True)
class LiftSlopeResult:
    """Finite-wing lift slope and Glauert factors of an untwisted wing [1]_.

    Attributes
    ----------
    lift_slope_per_rad : float
        ``a = d C_L / d alpha``.
    section_lift_slope_per_rad : float
        ``a_0`` used.
    lift_slope_factor : float
        ``tau`` of Eq. (5.70).
    induced_drag_factor : float
        ``delta`` of Eq. (5.61).
    span_efficiency : float
        ``1 / (1 + delta)``.
    aspect_ratio : float
        ``AR``.

    References
    ----------
    .. [1] ANDERSON-FOA-2017. J. D. Anderson Jr., Fundamentals of
           Aerodynamics, 6th ed., McGraw-Hill, 2017, Section 5.3.3,
           Eqs. (5.61), (5.69), (5.70).
    """

    lift_slope_per_rad: float
    section_lift_slope_per_rad: float
    lift_slope_factor: float
    induced_drag_factor: float
    span_efficiency: float
    aspect_ratio: float


def lift_slope(wing: Wing, *, terms: int = 40) -> LiftSlopeResult:
    """Lift slope, ``tau`` and ``delta`` of a wing with uniform section slope.

    Parameters
    ----------
    wing : Wing
        Untwisted wing with uniform ``a_0`` and ``alpha_L=0``.
    terms : int, default 40
        Series terms.

    Returns
    -------
    LiftSlopeResult
        ``a``, ``tau`` (Eq. 5.70 solved for ``tau``) and ``delta``.

    Raises
    ------
    DomainError
        If the wing is twisted or its section slope or zero-lift angle varies
        along the span (``tau`` is then not defined by Eq. 5.70).

    Notes
    -----
    ``C_L`` is linear in incidence, so ``a`` is obtained from one solution at
    ``alpha - alpha_L=0 = 1 rad`` and
    ``1 + tau = (a_0 / a - 1) pi AR / a_0`` [1]_.

    References
    ----------
    .. [1] ANDERSON-FOA-2017. J. D. Anderson Jr., Fundamentals of
           Aerodynamics, 6th ed., McGraw-Hill, 2017, Section 5.3.3,
           Eqs. (5.69)-(5.70).
    """

    half = 0.5 * wing.span
    probes = [half * (k / 10.0 - 0.95) for k in range(20)]
    slope0 = wing.section_lift_slope(0.0)
    zero0 = wing.zero_lift_angle(0.0)
    for y in probes:
        if wing.twist(y) != 0.0 or wing.section_lift_slope(y) != slope0 or wing.zero_lift_angle(y) != zero0:
            raise DomainError("lift_slope requires an untwisted wing with uniform section data.", field="wing")
    result = solve_lifting_line(wing, zero0 + 1.0, terms=terms)
    a = result.lift_coefficient
    tau = (slope0 / a - 1.0) * pi * wing.aspect_ratio / slope0 - 1.0
    return LiftSlopeResult(a, slope0, tau, result.induced_drag_factor, result.span_efficiency,
                           wing.aspect_ratio)


def finite_wing_lift_slope(section_lift_slope: float, aspect_ratio: float, lift_slope_factor: float = 0.0) -> float:
    """Finite-wing lift slope ``a = a_0 / (1 + (a_0 / (pi AR))(1 + tau))``.

    Parameters
    ----------
    section_lift_slope : float
        ``a_0`` per rad, positive.
    aspect_ratio : float
        ``AR``, positive.
    lift_slope_factor : float, default 0
        Glauert's ``tau`` (zero for the elliptic wing, Eq. 5.69).

    Returns
    -------
    float
        ``a`` per rad (Eq. 5.70).

    References
    ----------
    .. [1] ANDERSON-FOA-2017. J. D. Anderson Jr., Fundamentals of
           Aerodynamics, 6th ed., McGraw-Hill, 2017, Section 5.3.3,
           Eqs. (5.69)-(5.70).
    """

    a0 = _positive(section_lift_slope, "section_lift_slope")
    ar = _positive(aspect_ratio, "aspect_ratio")
    tau = _finite(lift_slope_factor, "lift_slope_factor")
    if tau <= -1.0:
        raise DomainError("Lift-slope factor must exceed -1.", field="lift_slope_factor")
    return a0 / (1.0 + a0 / (pi * ar) * (1.0 + tau))


def induced_drag_coefficient(lift_coefficient: float, aspect_ratio: float, induced_drag_factor: float = 0.0) -> float:
    """Induced drag coefficient ``C_D,i = C_L^2 (1 + delta) / (pi AR)``.

    Parameters
    ----------
    lift_coefficient : float
        ``C_L``.
    aspect_ratio : float
        ``AR``, positive.
    induced_drag_factor : float, default 0
        ``delta >= 0`` (zero for the elliptic lift distribution, Eq. 5.43).

    Returns
    -------
    float
        ``C_D,i`` (Eq. 5.61).

    References
    ----------
    .. [1] ANDERSON-FOA-2017. J. D. Anderson Jr., Fundamentals of
           Aerodynamics, 6th ed., McGraw-Hill, 2017, Section 5.3, Eqs. (5.43),
           (5.61).
    """

    cl = _finite(lift_coefficient, "lift_coefficient")
    ar = _positive(aspect_ratio, "aspect_ratio")
    delta = _finite(induced_drag_factor, "induced_drag_factor")
    if delta < 0.0:
        raise DomainError("Induced drag factor must be non-negative.", field="induced_drag_factor")
    return cl * cl * (1.0 + delta) / (pi * ar)
