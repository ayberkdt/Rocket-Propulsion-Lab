"""Inviscid hypersonic approximations: Newtonian flow, shock limits, similarity.

Newtonian and modified Newtonian surface pressure, the Newtonian flat plate,
the hypersonic (``M -> infinity``) limits of the oblique-shock relations and
the Mach-number-independence principle follow Anderson chapter 14 [1]_.
Hypersonic small-disturbance similarity and its unified supersonic-hypersonic
form, with the thin-wedge solution of Linnell, follow Van Dyke [2]_.
All relations assume a calorically perfect gas.

References
----------
.. [1] ANDERSON-FOA-2017. J. D. Anderson Jr., "Fundamentals of Aerodynamics",
       6th ed., McGraw-Hill, 2017, sections 14.3-14.6, pp. 877-898.
.. [2] NACA-TR-1194. M. D. Van Dyke, "A Study of Hypersonic
       Small-Disturbance Theory", NACA Report 1194, 1954, eqs. (1), (19),
       (33)-(35).
"""

from __future__ import annotations

from dataclasses import dataclass
from math import asin, atan, cos, degrees, isfinite, pi, radians, sin, sqrt, tan

from rocket_propulsion.core.errors import DomainError

from ._advanced_common import bracketed_root, maximize_unimodal
from .isentropic import calculate_isentropic
from .normal_shock import calculate_normal_shock


@dataclass(frozen=True, slots=True)
class NewtonianFlatPlateResult:
    """Newtonian flat-plate coefficients [1]_.

    Attributes
    ----------
    angle_of_attack_deg : float
        Angle of attack, degrees.
    maximum_pressure_coefficient : float
        ``Cp,max`` used on the windward side (2 for classical Newtonian).
    lower_pressure_coefficient, upper_pressure_coefficient : float
        ``Cp,max sin^2(alpha)`` and 0 (shadow).
    normal_force_coefficient, lift_coefficient, drag_coefficient : float
        ``c_n``, ``c_l = c_n cos(alpha)``, ``c_d = c_n sin(alpha)``.
    lift_to_drag_ratio : float
        ``cot(alpha)``.

    References
    ----------
    .. [1] ANDERSON-FOA-2017. J. D. Anderson Jr., "Fundamentals of
           Aerodynamics", 6th ed., McGraw-Hill, 2017, eqs. (14.8)-(14.16).
    """

    angle_of_attack_deg: float
    maximum_pressure_coefficient: float
    lower_pressure_coefficient: float
    upper_pressure_coefficient: float
    normal_force_coefficient: float
    lift_coefficient: float
    drag_coefficient: float
    lift_to_drag_ratio: float


@dataclass(frozen=True, slots=True)
class HypersonicShockLimits:
    """Exact oblique-shock ratios and their ``M -> infinity`` limits [1]_.

    Attributes
    ----------
    mach, wave_angle_deg, gamma : float
        Upstream Mach number, wave angle and ratio of specific heats.
    pressure_ratio, pressure_ratio_limit : float
        Eq. (14.28) and its limit eq. (14.29).
    density_ratio, density_ratio_limit : float
        Eq. (14.30) and eq. (14.31).
    temperature_ratio, temperature_ratio_limit : float
        ``(p2/p1)/(rho2/rho1)`` and eq. (14.32).
    pressure_coefficient, pressure_coefficient_limit : float
        Eq. (14.38) and eq. (14.39).

    References
    ----------
    .. [1] ANDERSON-FOA-2017. J. D. Anderson Jr., "Fundamentals of
           Aerodynamics", 6th ed., McGraw-Hill, 2017, section 14.5.
    """

    mach: float
    wave_angle_deg: float
    gamma: float
    pressure_ratio: float
    pressure_ratio_limit: float
    density_ratio: float
    density_ratio_limit: float
    temperature_ratio: float
    temperature_ratio_limit: float
    pressure_coefficient: float
    pressure_coefficient_limit: float


@dataclass(frozen=True, slots=True)
class SmallDisturbanceWedgeResult:
    """Hypersonic small-disturbance solution for a thin wedge [1]_.

    Attributes
    ----------
    mach, wedge_angle_deg, gamma : float
        Free-stream Mach number, wedge semivertex angle and gamma.
    unified : bool
        True when ``M`` was replaced by ``beta = sqrt(M^2 - 1)`` (eq. (35)).
    thickness_parameter : float
        ``delta = tan(wedge angle)`` (Van Dyke's choice for Fig. 5).
    similarity_parameter : float
        ``M delta`` or ``beta delta``.
    pressure_coefficient : float
        Eq. (34) or (35).
    shock_slope : float
        Shock-wave slope ``tau`` from eq. (33).
    shock_angle_deg : float
        ``atan(tau)`` in degrees.

    References
    ----------
    .. [1] NACA-TR-1194. M. D. Van Dyke, "A Study of Hypersonic
           Small-Disturbance Theory", NACA Report 1194, 1954, p. 8.
    """

    mach: float
    wedge_angle_deg: float
    gamma: float
    unified: bool
    thickness_parameter: float
    similarity_parameter: float
    pressure_coefficient: float
    shock_slope: float
    shock_angle_deg: float


def _check_gamma(gamma: float) -> None:
    if not isfinite(gamma) or gamma <= 1.0:
        raise DomainError("Gamma must be finite and greater than 1.", field="gamma")


def _check_angle(angle_deg: float, field: str) -> float:
    if not isfinite(angle_deg) or not -90.0 <= angle_deg <= 90.0:
        raise DomainError("Angle must be finite and within [-90, 90] degrees.", field=field)
    return radians(angle_deg)


def newtonian_pressure_coefficient(incidence_deg: float) -> float:
    """Newton's sine-squared law ``Cp = 2 sin^2(theta)`` [1]_.

    Parameters
    ----------
    incidence_deg : float
        Angle between the local surface tangent and the free stream, degrees;
        non-positive values are in the shadow (``Cp = 0``, eq. (14.9)).

    Returns
    -------
    float
        Pressure coefficient, Anderson eq. (14.4).

    References
    ----------
    .. [1] ANDERSON-FOA-2017. J. D. Anderson Jr., "Fundamentals of
           Aerodynamics", 6th ed., McGraw-Hill, 2017, pp. 879, 883.
    """
    theta = _check_angle(incidence_deg, "incidence_deg")
    return 0.0 if theta <= 0.0 else 2.0 * sin(theta) ** 2


def stagnation_pressure_coefficient(mach: float, gamma: float = 1.4) -> float:
    """Stagnation-point pressure coefficient ``Cp,max`` [1]_.

    Parameters
    ----------
    mach : float
        Free-stream Mach number, finite and positive.
    gamma : float, default 1.4
        Ratio of specific heats.

    Returns
    -------
    float
        ``(2/(gamma M^2)) (p0,2/p_inf - 1)``; for ``M > 1`` the stagnation
        pressure behind a normal shock, otherwise the isentropic ``p0``.

    Notes
    -----
    Anderson p. 880 [1]_: ``p0,2/p_inf`` from normal-shock theory and
    ``Cp,max = (2/(gamma M^2))(p0,2/p_inf - 1)``.

    References
    ----------
    .. [1] ANDERSON-FOA-2017. J. D. Anderson Jr., "Fundamentals of
           Aerodynamics", 6th ed., McGraw-Hill, 2017, pp. 879-880.
    """
    _check_gamma(gamma)
    if not isfinite(mach) or mach <= 0.0:
        raise DomainError("Mach number must be finite and positive.", field="mach")
    total_over_static = 1.0 / calculate_isentropic(mach, gamma).pressure_total_ratio
    if mach > 1.0:
        total_over_static *= calculate_normal_shock(mach, gamma).stagnation_pressure_ratio
    return 2.0 / (gamma * mach * mach) * (total_over_static - 1.0)


def modified_newtonian_pressure_coefficient(
    incidence_deg: float, mach: float, gamma: float = 1.4
) -> float:
    """Modified Newtonian law ``Cp = Cp,max sin^2(theta)``, Anderson eq. (14.7) [1]_.

    Parameters
    ----------
    incidence_deg : float
        Surface incidence, degrees (shadow for non-positive values).
    mach : float
        Free-stream Mach number.
    gamma : float, default 1.4
        Ratio of specific heats.

    Returns
    -------
    float
        Pressure coefficient.

    References
    ----------
    .. [1] ANDERSON-FOA-2017. J. D. Anderson Jr., "Fundamentals of
           Aerodynamics", 6th ed., McGraw-Hill, 2017, p. 880.
    """
    theta = _check_angle(incidence_deg, "incidence_deg")
    if theta <= 0.0:
        return 0.0
    return stagnation_pressure_coefficient(mach, gamma) * sin(theta) ** 2


def newtonian_flat_plate(
    angle_of_attack_deg: float, maximum_pressure_coefficient: float = 2.0
) -> NewtonianFlatPlateResult:
    """Newtonian flat-plate lift, drag and L/D [1]_.

    Parameters
    ----------
    angle_of_attack_deg : float
        Angle of attack in ``(0, 90]`` degrees.
    maximum_pressure_coefficient : float, default 2.0
        2 for Newton's law; pass :func:`stagnation_pressure_coefficient` for
        the modified law.

    Returns
    -------
    NewtonianFlatPlateResult
        Coefficients from Anderson eqs. (14.8)-(14.16).

    References
    ----------
    .. [1] ANDERSON-FOA-2017. J. D. Anderson Jr., "Fundamentals of
           Aerodynamics", 6th ed., McGraw-Hill, 2017, pp. 883-884.
    """
    alpha = _check_angle(angle_of_attack_deg, "angle_of_attack_deg")
    if alpha <= 0.0:
        raise DomainError("Angle of attack must be positive.", field="angle_of_attack_deg")
    if not isfinite(maximum_pressure_coefficient) or maximum_pressure_coefficient <= 0.0:
        raise DomainError("Cp,max must be finite and positive.",
                          field="maximum_pressure_coefficient")
    lower = maximum_pressure_coefficient * sin(alpha) ** 2
    normal = lower
    return NewtonianFlatPlateResult(
        angle_of_attack_deg=angle_of_attack_deg,
        maximum_pressure_coefficient=maximum_pressure_coefficient,
        lower_pressure_coefficient=lower,
        upper_pressure_coefficient=0.0,
        normal_force_coefficient=normal,
        lift_coefficient=normal * cos(alpha),
        drag_coefficient=normal * sin(alpha),
        lift_to_drag_ratio=cos(alpha) / sin(alpha),
    )


def newtonian_maximum_lift() -> tuple[float, float]:
    """Angle and value of the Newtonian flat-plate maximum lift [1]_.

    Returns
    -------
    tuple of float
        ``(alpha_deg, c_l,max)`` with ``sin^2(alpha) = 2/3``.

    References
    ----------
    .. [1] ANDERSON-FOA-2017. J. D. Anderson Jr., "Fundamentals of
           Aerodynamics", 6th ed., McGraw-Hill, 2017, p. 886.
    """
    alpha = asin(sqrt(2.0 / 3.0))
    return degrees(alpha), 2.0 * sin(alpha) ** 2 * cos(alpha)


def newtonian_maximum_lift_to_drag(zero_lift_drag_coefficient: float) -> tuple[float, float]:
    """Small-angle Newtonian ``(L/D)max`` with skin friction [1]_.

    Parameters
    ----------
    zero_lift_drag_coefficient : float
        Friction drag coefficient ``c_d,0``, positive.

    Returns
    -------
    tuple of float
        ``(alpha_deg, (L/D)max)`` from ``alpha = c_d,0^(1/3)`` (eq. (14.23))
        and ``(L/D)max = (2/3)/c_d,0^(1/3)`` (eq. (14.24), printed as 0.67).

    References
    ----------
    .. [1] ANDERSON-FOA-2017. J. D. Anderson Jr., "Fundamentals of
           Aerodynamics", 6th ed., McGraw-Hill, 2017, pp. 887-888.
    """
    if not isfinite(zero_lift_drag_coefficient) or zero_lift_drag_coefficient <= 0.0:
        raise DomainError("c_d,0 must be finite and positive.",
                          field="zero_lift_drag_coefficient")
    root = zero_lift_drag_coefficient ** (1.0 / 3.0)
    return degrees(root), (2.0 / 3.0) / root


def hypersonic_shock_limits(
    mach: float, wave_angle_deg: float, gamma: float = 1.4
) -> HypersonicShockLimits:
    """Exact oblique-shock ratios beside their hypersonic limits [1]_.

    Parameters
    ----------
    mach : float
        Upstream Mach number with ``M sin(beta) > 1``.
    wave_angle_deg : float
        Shock-wave angle, degrees.
    gamma : float, default 1.4
        Ratio of specific heats.

    Returns
    -------
    HypersonicShockLimits
        Anderson eqs. (14.28)-(14.32), (14.38), (14.39).

    References
    ----------
    .. [1] ANDERSON-FOA-2017. J. D. Anderson Jr., "Fundamentals of
           Aerodynamics", 6th ed., McGraw-Hill, 2017, pp. 892-894.
    """
    _check_gamma(gamma)
    beta = _check_angle(wave_angle_deg, "wave_angle_deg")
    if not isfinite(mach) or mach * sin(beta) <= 1.0:
        raise DomainError("The normal Mach number M sin(beta) must exceed 1.", field="mach")
    s2 = sin(beta) ** 2
    mn2 = mach * mach * s2
    pressure = 1.0 + 2.0 * gamma / (gamma + 1.0) * (mn2 - 1.0)
    pressure_limit = 2.0 * gamma / (gamma + 1.0) * mn2
    density = (gamma + 1.0) * mn2 / ((gamma - 1.0) * mn2 + 2.0)
    density_limit = (gamma + 1.0) / (gamma - 1.0)
    return HypersonicShockLimits(
        mach=mach,
        wave_angle_deg=wave_angle_deg,
        gamma=gamma,
        pressure_ratio=pressure,
        pressure_ratio_limit=pressure_limit,
        density_ratio=density,
        density_ratio_limit=density_limit,
        temperature_ratio=pressure / density,
        temperature_ratio_limit=2.0 * gamma * (gamma - 1.0) / (gamma + 1.0) ** 2 * mn2,
        pressure_coefficient=4.0 / (gamma + 1.0) * (s2 - 1.0 / (mach * mach)),
        pressure_coefficient_limit=4.0 / (gamma + 1.0) * s2,
    )


def limiting_wave_angle_ratio(gamma: float = 1.4) -> float:
    """Thin-wedge hypersonic limit ``beta/theta = (gamma + 1)/2``, eq. (14.35) [1]_.

    References
    ----------
    .. [1] ANDERSON-FOA-2017. J. D. Anderson Jr., "Fundamentals of
           Aerodynamics", 6th ed., McGraw-Hill, 2017, p. 893.
    """
    _check_gamma(gamma)
    return 0.5 * (gamma + 1.0)


def mach_independent_wedge(deflection_deg: float, gamma: float = 1.4) -> tuple[float, float]:
    """Wedge wave angle and ``Cp`` in the limit ``M -> infinity`` [1]_.

    Parameters
    ----------
    deflection_deg : float
        Wedge semivertex angle, degrees, below the limiting detachment angle.
    gamma : float, default 1.4
        Ratio of specific heats.

    Returns
    -------
    tuple of float
        ``(wave_angle_deg, Cp)``.

    Raises
    ------
    DomainError
        If the deflection exceeds the ``M -> infinity`` attachment limit.

    Notes
    -----
    Dividing the theta-beta-M relation (Anderson eq. (9.23)) by ``M^2`` and
    letting ``M -> infinity`` leaves
    ``tan(theta) = 2 cot(beta) sin^2(beta)/(gamma + cos(2 beta))``, which is
    solved for the weak ``beta``; ``Cp`` is eq. (14.39). Both are free of
    ``M`` - the Mach-number-independence principle of section 14.6 [1]_.

    References
    ----------
    .. [1] ANDERSON-FOA-2017. J. D. Anderson Jr., "Fundamentals of
           Aerodynamics", 6th ed., McGraw-Hill, 2017, pp. 892-897.
    """
    _check_gamma(gamma)
    theta = _check_angle(deflection_deg, "deflection_deg")
    if theta <= 0.0:
        raise DomainError("Deflection must be positive.", field="deflection_deg")

    def limit_deflection(beta: float) -> float:
        return atan(2.0 / tan(beta) * sin(beta) ** 2 / (gamma + cos(2.0 * beta)))

    peak_beta, peak_theta = maximize_unimodal(limit_deflection, 1.0e-9, 0.5 * pi - 1.0e-9)
    if theta > peak_theta:
        raise DomainError(
            f"Detached shock: limiting maximum deflection is {degrees(peak_theta):.6g} degrees.",
            code="detached_shock",
        )
    beta = bracketed_root(lambda b: limit_deflection(b) - theta, 1.0e-12, peak_beta,
                          label="limiting wave angle")
    return degrees(beta), 4.0 / (gamma + 1.0) * sin(beta) ** 2


def hypersonic_similarity_parameter(mach: float, thickness: float, *, unified: bool = False) -> float:
    """Hypersonic similarity parameter ``K = M tau`` (or ``beta tau``) [1]_.

    Parameters
    ----------
    mach : float
        Free-stream Mach number (greater than 1 when ``unified``).
    thickness : float
        Thickness ratio or slope ``tau``, positive.
    unified : bool, default False
        Use ``beta = sqrt(M^2 - 1)`` (unified supersonic-hypersonic rule).

    Returns
    -------
    float
        Similarity parameter. Bodies of affinely related shape with equal
        ``K`` have equal ``Cp/tau^2`` (Van Dyke eq. (19b) and the similarity
        rule on p. 6) [1]_.

    References
    ----------
    .. [1] NACA-TR-1194. M. D. Van Dyke, "A Study of Hypersonic
           Small-Disturbance Theory", NACA Report 1194, 1954, pp. 2, 6.
    """
    if not isfinite(thickness) or thickness <= 0.0:
        raise DomainError("Thickness parameter must be finite and positive.", field="thickness")
    if not isfinite(mach) or mach <= (1.0 if unified else 0.0):
        raise DomainError("Mach number out of range.", field="mach")
    speed = sqrt(mach * mach - 1.0) if unified else mach
    return speed * thickness


def small_disturbance_wedge(
    mach: float, wedge_angle_deg: float, gamma: float = 1.4, *, unified: bool = True
) -> SmallDisturbanceWedgeResult:
    """Linnell's thin-wedge result in hypersonic small-disturbance theory [1]_.

    Parameters
    ----------
    mach : float
        Free-stream Mach number (greater than 1).
    wedge_angle_deg : float
        Wedge semivertex angle, degrees, in ``(0, 45)``.
    gamma : float, default 1.4
        Ratio of specific heats.
    unified : bool, default True
        Replace ``M`` by ``beta = sqrt(M^2 - 1)`` (Van Dyke eq. (35)),
        which reduces to linearized theory as ``beta delta -> 0``.

    Returns
    -------
    SmallDisturbanceWedgeResult
        ``Cp = delta^2 [(gamma+1)/2 + sqrt(((gamma+1)/2)^2 + 4/K^2)]``
        (eq. (34)/(35)) and the shock slope
        ``tau = delta [(gamma+1)/4 + sqrt(((gamma+1)/4)^2 + 1/K^2)]`` (eq. (33)).

    Notes
    -----
    ``delta = tan(wedge angle)`` as in Van Dyke's Fig. 5. Applying the
    unified replacement to eq. (33) follows the general rule of eq. (19a)
    (replace ``M tau`` by ``beta tau``) [1]_.

    References
    ----------
    .. [1] NACA-TR-1194. M. D. Van Dyke, "A Study of Hypersonic
           Small-Disturbance Theory", NACA Report 1194, 1954, pp. 6-9.
    """
    _check_gamma(gamma)
    if not isfinite(mach) or mach <= 1.0:
        raise DomainError("Mach number must be greater than 1.", field="mach")
    if not isfinite(wedge_angle_deg) or not 0.0 < wedge_angle_deg < 45.0:
        raise DomainError("Wedge angle must be in (0, 45) degrees.", field="wedge_angle_deg")
    delta = tan(radians(wedge_angle_deg))
    k = hypersonic_similarity_parameter(mach, delta, unified=unified)
    half = 0.5 * (gamma + 1.0)
    quarter = 0.25 * (gamma + 1.0)
    cp = delta * delta * (half + sqrt(half * half + 4.0 / (k * k)))
    slope = delta * (quarter + sqrt(quarter * quarter + 1.0 / (k * k)))
    return SmallDisturbanceWedgeResult(
        mach=mach,
        wedge_angle_deg=wedge_angle_deg,
        gamma=gamma,
        unified=unified,
        thickness_parameter=delta,
        similarity_parameter=k,
        pressure_coefficient=cp,
        shock_slope=slope,
        shock_angle_deg=degrees(atan(slope)),
    )
