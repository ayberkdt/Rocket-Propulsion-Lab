"""Linearized compressible-flow corrections and thin-airfoil supersonic theory.

Subsonic: the Prandtl-Glauert, Karman-Tsien and Laitone compressibility
corrections and the critical-Mach-number construction of Anderson section
11.6 [1]_. Supersonic: Ackeret's linearized pressure coefficient and the
flat-plate lift and wave drag (Anderson chapter 12) [1]_, and Busemann's
second-order pressure coefficient, whose coefficients are the first two
terms of the weak oblique-shock series of NACA Report 1135, eq. (152), which
the report shows is shared by the Prandtl-Meyer expansion series, eq. (174),
to second order [2]_. Validity notes follow Anderson (Prandtl-Glauert above
``M = 0.7``, transonic band ``0.8 <= M <= 1.2``) [1]_ and the small-disturbance
criterion of Van Dyke (``beta tau << 1`` for linear theory) [3]_.

References
----------
.. [1] ANDERSON-FOA-2017. J. D. Anderson Jr., "Fundamentals of Aerodynamics",
       6th ed., McGraw-Hill, 2017, sections 11.4-11.6 (pp. 750-764) and
       12.2-12.3 (pp. 812-819).
.. [2] NACA-TR-1135. Ames Research Staff, "Equations, Tables, and Charts for
       Compressible Flow", NACA Report 1135, 1953, eqs. (151)-(152), (174).
.. [3] NACA-TR-1194. M. D. Van Dyke, "A Study of Hypersonic Small-Disturbance
       Theory", NACA Report 1194, 1954, eq. (1).
"""

from __future__ import annotations

from dataclasses import dataclass
from math import degrees, isfinite, radians, sqrt

from rocket_propulsion.core.errors import DomainError
from rocket_propulsion.core.metadata import WarningMessage

from ._advanced_common import bracketed_root
from .shock_expansion import (
    AirfoilSection,
    PanelSolution,
    SectionCoefficients,
    _inclination,
    _panels,
    integrate_section_forces,
)

_RULES = ("prandtl-glauert", "karman-tsien", "laitone")


@dataclass(frozen=True, slots=True)
class CriticalMachResult:
    """Critical Mach number from a low-speed minimum pressure coefficient [1]_.

    Attributes
    ----------
    rule : str
        Compressibility correction used.
    incompressible_minimum_cp : float
        Low-speed minimum ``Cp,0`` on the section.
    critical_mach : float
        Free-stream Mach number at which the minimum-pressure point is sonic.
    critical_pressure_coefficient : float
        ``Cp,cr`` at that Mach number.

    References
    ----------
    .. [1] ANDERSON-FOA-2017. J. D. Anderson Jr., "Fundamentals of
           Aerodynamics", 6th ed., McGraw-Hill, 2017, section 11.6.
    """

    rule: str
    incompressible_minimum_cp: float
    critical_mach: float
    critical_pressure_coefficient: float


@dataclass(frozen=True, slots=True)
class FlatPlateLinearResult:
    """Ackeret flat-plate lift and wave drag [1]_.

    Attributes
    ----------
    mach, angle_of_attack_deg : float
        Free-stream conditions.
    lift_coefficient, drag_coefficient : float
        ``4 alpha/beta`` and ``4 alpha^2/beta`` (alpha in radians).
    warnings : tuple of WarningMessage
        Validity notes.

    References
    ----------
    .. [1] ANDERSON-FOA-2017. J. D. Anderson Jr., "Fundamentals of
           Aerodynamics", 6th ed., McGraw-Hill, 2017, eqs. (12.23)-(12.24).
    """

    mach: float
    angle_of_attack_deg: float
    lift_coefficient: float
    drag_coefficient: float
    warnings: tuple[WarningMessage, ...] = ()


def _subsonic_mach(mach: float) -> None:
    if not isfinite(mach) or not 0.0 <= mach < 1.0:
        raise DomainError("Free-stream Mach number must be in [0, 1).", field="mach")


def _supersonic_mach(mach: float) -> None:
    if not isfinite(mach) or mach <= 1.0:
        raise DomainError("Free-stream Mach number must be greater than 1.", field="mach")


def _gamma(gamma: float) -> None:
    if not isfinite(gamma) or gamma <= 1.0:
        raise DomainError("Gamma must be finite and greater than 1.", field="gamma")


def subsonic_validity_warnings(mach: float) -> tuple[WarningMessage, ...]:
    """Validity notes for the subsonic linearized corrections [1]_.

    Parameters
    ----------
    mach : float
        Free-stream Mach number.

    Returns
    -------
    tuple of WarningMessage
        ``prandtl_glauert_high_mach`` at ``M >= 0.7`` and ``transonic_regime``
        at ``M >= 0.8``.

    References
    ----------
    .. [1] ANDERSON-FOA-2017. J. D. Anderson Jr., "Fundamentals of
           Aerodynamics", 6th ed., McGraw-Hill, 2017, pp. 751, 757.
    """
    warnings = []
    if mach >= 0.7:
        warnings.append(WarningMessage(
            "prandtl_glauert_high_mach",
            "Anderson notes that Prandtl-Glauert begins to give inappropriate results at "
            "M = 0.7 and above.", "warning", "mach"))
    if mach >= 0.8:
        warnings.append(WarningMessage(
            "transonic_regime",
            "Linearized theory does not apply in the transonic band 0.8 <= M <= 1.2.",
            "critical", "mach"))
    return tuple(warnings)


def prandtl_glauert(incompressible_value: float, mach: float) -> float:
    """Prandtl-Glauert rule, ``C = C_0/sqrt(1 - M^2)`` [1]_.

    Parameters
    ----------
    incompressible_value : float
        Low-speed ``Cp``, ``c_l`` or ``c_m``.
    mach : float
        Free-stream Mach number in ``[0, 1)``.

    Returns
    -------
    float
        Compressible value.

    Notes
    -----
    Anderson eqs. (11.51)-(11.53) [1]_.

    References
    ----------
    .. [1] ANDERSON-FOA-2017. J. D. Anderson Jr., "Fundamentals of
           Aerodynamics", 6th ed., McGraw-Hill, 2017, p. 754.
    """
    _subsonic_mach(mach)
    if not isfinite(incompressible_value):
        raise DomainError("Incompressible value must be finite.", field="incompressible_value")
    return incompressible_value / sqrt(1.0 - mach * mach)


def karman_tsien(incompressible_cp: float, mach: float) -> float:
    """Karman-Tsien rule, Anderson eq. (11.54) [1]_.

    Parameters
    ----------
    incompressible_cp : float
        Low-speed pressure coefficient ``Cp,0``.
    mach : float
        Free-stream Mach number in ``[0, 1)``.

    Returns
    -------
    float
        ``Cp = Cp,0 / (sqrt(1 - M^2) + [M^2/(1 + sqrt(1 - M^2))] Cp,0/2)``.

    References
    ----------
    .. [1] ANDERSON-FOA-2017. J. D. Anderson Jr., "Fundamentals of
           Aerodynamics", 6th ed., McGraw-Hill, 2017, p. 755.
    """
    _subsonic_mach(mach)
    if not isfinite(incompressible_cp):
        raise DomainError("Cp,0 must be finite.", field="incompressible_cp")
    beta = sqrt(1.0 - mach * mach)
    return incompressible_cp / (beta + mach * mach / (1.0 + beta) * incompressible_cp / 2.0)


def laitone(incompressible_cp: float, mach: float, gamma: float = 1.4) -> float:
    """Laitone rule, Anderson eq. (11.55) [1]_.

    Parameters
    ----------
    incompressible_cp : float
        Low-speed pressure coefficient ``Cp,0``.
    mach : float
        Free-stream Mach number in ``[0, 1)``.
    gamma : float, default 1.4
        Ratio of specific heats.

    Returns
    -------
    float
        ``Cp = Cp,0 / (sqrt(1 - M^2) + M^2 (1 + (gamma-1) M^2/2) Cp,0
        / (2 sqrt(1 - M^2)))``.

    References
    ----------
    .. [1] ANDERSON-FOA-2017. J. D. Anderson Jr., "Fundamentals of
           Aerodynamics", 6th ed., McGraw-Hill, 2017, p. 755.
    """
    _subsonic_mach(mach)
    _gamma(gamma)
    if not isfinite(incompressible_cp):
        raise DomainError("Cp,0 must be finite.", field="incompressible_cp")
    beta = sqrt(1.0 - mach * mach)
    nonlinear = mach * mach * (1.0 + 0.5 * (gamma - 1.0) * mach * mach) / (2.0 * beta)
    return incompressible_cp / (beta + nonlinear * incompressible_cp)


def corrected_pressure_coefficient(
    incompressible_cp: float, mach: float, rule: str = "prandtl-glauert", gamma: float = 1.4
) -> float:
    """Apply a named subsonic compressibility correction [1]_.

    Parameters
    ----------
    incompressible_cp : float
        Low-speed ``Cp,0``.
    mach : float
        Free-stream Mach number in ``[0, 1)``.
    rule : {"prandtl-glauert", "karman-tsien", "laitone"}
        Correction (case-insensitive).
    gamma : float, default 1.4
        Ratio of specific heats (Laitone only).

    Returns
    -------
    float
        Corrected pressure coefficient.

    References
    ----------
    .. [1] ANDERSON-FOA-2017. J. D. Anderson Jr., "Fundamentals of
           Aerodynamics", 6th ed., McGraw-Hill, 2017, eqs. (11.51), (11.54),
           (11.55).
    """
    name = rule.strip().lower()
    if name == "prandtl-glauert":
        return prandtl_glauert(incompressible_cp, mach)
    if name == "karman-tsien":
        return karman_tsien(incompressible_cp, mach)
    if name == "laitone":
        return laitone(incompressible_cp, mach, gamma)
    raise DomainError(f"Unknown compressibility rule {rule!r}; use one of {_RULES}.",
                      field="rule")


def critical_pressure_coefficient(mach: float, gamma: float = 1.4) -> float:
    """Pressure coefficient at a locally sonic point, Anderson eq. (11.59) [1]_.

    Parameters
    ----------
    mach : float
        Free-stream Mach number, finite and positive.
    gamma : float, default 1.4
        Ratio of specific heats.

    Returns
    -------
    float
        ``Cp,cr = 2/(gamma M^2) [((1 + (gamma-1)/2 M^2)/(1 + (gamma-1)/2))^(gamma/(gamma-1)) - 1]``.

    References
    ----------
    .. [1] ANDERSON-FOA-2017. J. D. Anderson Jr., "Fundamentals of
           Aerodynamics", 6th ed., McGraw-Hill, 2017, eqs. (11.59)-(11.60).
    """
    if not isfinite(mach) or mach <= 0.0:
        raise DomainError("Mach number must be finite and positive.", field="mach")
    _gamma(gamma)
    ratio = (1.0 + 0.5 * (gamma - 1.0) * mach * mach) / (1.0 + 0.5 * (gamma - 1.0))
    return 2.0 / (gamma * mach * mach) * (ratio ** (gamma / (gamma - 1.0)) - 1.0)


def critical_mach_number(
    incompressible_minimum_cp: float, rule: str = "prandtl-glauert", gamma: float = 1.4
) -> CriticalMachResult:
    """Critical Mach number of a section from its low-speed minimum ``Cp`` [1]_.

    Parameters
    ----------
    incompressible_minimum_cp : float
        Low-speed minimum pressure coefficient, negative.
    rule : {"prandtl-glauert", "karman-tsien", "laitone"}
        Compressibility correction for the minimum-pressure point.
    gamma : float, default 1.4
        Ratio of specific heats.

    Returns
    -------
    CriticalMachResult
        Critical Mach number and pressure coefficient.

    Raises
    ------
    DomainError
        If ``Cp,0`` is not negative or the intersection is not bracketed.

    Notes
    -----
    Solves the analytical form of the Anderson construction, eq. (11.63):
    the corrected minimum ``Cp`` equals ``Cp,cr`` of eq. (11.60), by bisection
    on ``0 < M < 1`` [1]_. For the Karman-Tsien and Laitone rules the
    correction becomes singular below ``M = 1`` when ``Cp,0 < 0``; beyond the
    singularity the mismatch is treated as negative so that the first
    (physical) intersection is returned.

    References
    ----------
    .. [1] ANDERSON-FOA-2017. J. D. Anderson Jr., "Fundamentals of
           Aerodynamics", 6th ed., McGraw-Hill, 2017, pp. 758-763.
    """
    if not isfinite(incompressible_minimum_cp) or incompressible_minimum_cp >= 0.0:
        raise DomainError("The low-speed minimum Cp must be negative.",
                          field="incompressible_minimum_cp")
    name = rule.strip().lower()

    corrected_pressure_coefficient(incompressible_minimum_cp, 0.5, name, gamma)

    def mismatch(mach: float) -> float:
        corrected = corrected_pressure_coefficient(incompressible_minimum_cp, mach, name, gamma)
        if corrected > 0.0:
            # Karman-Tsien and Laitone denominators vanish below M = 1 for a
            # negative Cp,0; past that point the corrected suction is unbounded.
            return -1.0
        return corrected - critical_pressure_coefficient(mach, gamma)

    mach = bracketed_root(mismatch, 1.0e-3, 1.0 - 1.0e-9, x_tolerance=1.0e-14,
                          label="critical Mach number")
    return CriticalMachResult(name, incompressible_minimum_cp, mach,
                              critical_pressure_coefficient(mach, gamma))


def supersonic_validity_warnings(
    mach: float, max_inclination_rad: float
) -> tuple[WarningMessage, ...]:
    """Validity notes for linearized supersonic theory [1]_ [2]_.

    Parameters
    ----------
    mach : float
        Free-stream Mach number.
    max_inclination_rad : float
        Largest absolute surface inclination to the free stream, radians.

    Returns
    -------
    tuple of WarningMessage
        ``transonic_regime`` for ``M < 1.2`` and ``hypersonic_similarity``
        when ``beta tau >= 0.3``.

    Notes
    -----
    Van Dyke states linear theory is adequate when ``beta tau << 1`` [2]_;
    the 0.3 threshold is this library's flag level, not a sourced value.

    References
    ----------
    .. [1] ANDERSON-FOA-2017. J. D. Anderson Jr., "Fundamentals of
           Aerodynamics", 6th ed., McGraw-Hill, 2017, p. 757.
    .. [2] NACA-TR-1194. M. D. Van Dyke, "A Study of Hypersonic
           Small-Disturbance Theory", NACA Report 1194, 1954, p. 2.
    """
    warnings = []
    if mach < 1.2:
        warnings.append(WarningMessage(
            "transonic_regime",
            "Linearized theory does not apply in the transonic band 0.8 <= M <= 1.2.",
            "critical", "mach"))
    parameter = sqrt(max(mach * mach - 1.0, 0.0)) * max_inclination_rad
    if parameter >= 0.3:
        warnings.append(WarningMessage(
            "hypersonic_similarity",
            f"beta*tau = {parameter:.3g}; linear theory needs beta*tau << 1 (Van Dyke).",
            "warning", "mach"))
    return tuple(warnings)


def ackeret_pressure_coefficient(inclination_rad: float, mach: float) -> float:
    """Linearized supersonic ``Cp = 2 theta/sqrt(M^2 - 1)``, Anderson eq. (12.15) [1]_.

    Parameters
    ----------
    inclination_rad : float
        Surface inclination into the free stream (positive compression), radians.
    mach : float
        Free-stream Mach number, greater than 1.

    Returns
    -------
    float
        Pressure coefficient.

    References
    ----------
    .. [1] ANDERSON-FOA-2017. J. D. Anderson Jr., "Fundamentals of
           Aerodynamics", 6th ed., McGraw-Hill, 2017, p. 814.
    """
    _supersonic_mach(mach)
    return 2.0 * inclination_rad / sqrt(mach * mach - 1.0)


def busemann_coefficients(mach: float, gamma: float = 1.4) -> tuple[float, float]:
    """Busemann coefficients ``(C1, C2)`` of ``Cp = C1 theta + C2 theta^2`` [1]_.

    Parameters
    ----------
    mach : float
        Free-stream Mach number, greater than 1.
    gamma : float, default 1.4
        Ratio of specific heats.

    Returns
    -------
    tuple of float
        ``C1 = 2/sqrt(M^2 - 1)`` and
        ``C2 = ((gamma+1) M^4 - 4 (M^2 - 1)) / (2 (M^2 - 1)^2)``.

    Notes
    -----
    The first two terms of NACA 1135 eq. (152) for ``(p2 - p1)/q1`` across a
    weak oblique shock of deflection ``delta``; eq. (174) shows the
    Prandtl-Meyer expansion series agrees with eq. (151) up to and including
    the square term, so the same coefficients hold for ``theta < 0`` [1]_.

    References
    ----------
    .. [1] NACA-TR-1135. Ames Research Staff, "Equations, Tables, and Charts
           for Compressible Flow", NACA Report 1135, 1953, pp. 622, 626.
    """
    _supersonic_mach(mach)
    _gamma(gamma)
    m2 = mach * mach
    return 2.0 / sqrt(m2 - 1.0), ((gamma + 1.0) * m2 * m2 - 4.0 * (m2 - 1.0)) / (
        2.0 * (m2 - 1.0) ** 2
    )


def busemann_pressure_coefficient(inclination_rad: float, mach: float, gamma: float = 1.4) -> float:
    """Second-order (Busemann) pressure coefficient [1]_.

    Parameters
    ----------
    inclination_rad : float
        Surface inclination into the free stream, radians.
    mach : float
        Free-stream Mach number, greater than 1.
    gamma : float, default 1.4
        Ratio of specific heats.

    Returns
    -------
    float
        ``C1 theta + C2 theta^2`` with :func:`busemann_coefficients`.

    References
    ----------
    .. [1] NACA-TR-1135. Ames Research Staff, "Equations, Tables, and Charts
           for Compressible Flow", NACA Report 1135, 1953, eq. (152).
    """
    c1, c2 = busemann_coefficients(mach, gamma)
    return c1 * inclination_rad + c2 * inclination_rad * inclination_rad


def ackeret_flat_plate(mach: float, angle_of_attack_deg: float) -> FlatPlateLinearResult:
    """Linearized flat-plate lift and wave drag, Anderson eqs. (12.23)-(12.24) [1]_.

    Parameters
    ----------
    mach : float
        Free-stream Mach number, greater than 1.
    angle_of_attack_deg : float
        Angle of attack in degrees.

    Returns
    -------
    FlatPlateLinearResult
        ``c_l = 4 alpha/beta`` and ``c_d = 4 alpha^2/beta``.

    References
    ----------
    .. [1] ANDERSON-FOA-2017. J. D. Anderson Jr., "Fundamentals of
           Aerodynamics", 6th ed., McGraw-Hill, 2017, p. 818.
    """
    _supersonic_mach(mach)
    if not isfinite(angle_of_attack_deg):
        raise DomainError("Angle of attack must be finite.", field="angle_of_attack_deg")
    alpha = radians(angle_of_attack_deg)
    beta = sqrt(mach * mach - 1.0)
    return FlatPlateLinearResult(
        mach, angle_of_attack_deg, 4.0 * alpha / beta, 4.0 * alpha * alpha / beta,
        supersonic_validity_warnings(mach, abs(alpha)),
    )


def linearized_section_coefficients(
    section: AirfoilSection,
    mach: float,
    angle_of_attack_deg: float,
    *,
    order: int = 1,
    gamma: float = 1.4,
    base_pressure_ratio: float = 1.0,
    moment_reference: tuple[float, float] | None = None,
) -> SectionCoefficients:
    """Ackeret (first-order) or Busemann (second-order) section coefficients [1]_ [2]_.

    Parameters
    ----------
    section : AirfoilSection
        Polygonal geometry (see :mod:`.shock_expansion`).
    mach : float
        Free-stream Mach number, greater than 1.
    angle_of_attack_deg : float
        Angle of attack, degrees.
    order : {1, 2}, default 1
        1 for Ackeret, eq. (12.15) [1]_; 2 for Busemann, NACA 1135 eq. (152) [2]_.
    gamma : float, default 1.4
        Ratio of specific heats (second order only).
    base_pressure_ratio : float, default 1.0
        ``p_base/p_inf`` on a blunt base.
    moment_reference : (float, float), optional
        Moment reference; defaults to the leading edge.

    Returns
    -------
    SectionCoefficients
        Coefficients, panel ``Cp`` values and validity warnings.

    Notes
    -----
    Each panel's ``Cp`` depends only on its own inclination to the free
    stream (no wave interaction), and the panel forces are integrated
    exactly with :func:`integrate_section_forces` rather than with the
    small-angle projections ``c_l = c_n``, ``c_d = c_n alpha`` of Anderson
    eqs. (12.21)-(12.22); the difference is of higher order in the angle.

    References
    ----------
    .. [1] ANDERSON-FOA-2017. J. D. Anderson Jr., "Fundamentals of
           Aerodynamics", 6th ed., McGraw-Hill, 2017, section 12.3.
    .. [2] NACA-TR-1135. Ames Research Staff, "Equations, Tables, and Charts
           for Compressible Flow", NACA Report 1135, 1953, eq. (152).
    """
    _supersonic_mach(mach)
    _gamma(gamma)
    if order not in (1, 2):
        raise DomainError("Order must be 1 (Ackeret) or 2 (Busemann).", field="order")
    if not isfinite(angle_of_attack_deg):
        raise DomainError("Angle of attack must be finite.", field="angle_of_attack_deg")
    alpha = radians(angle_of_attack_deg)
    dynamic = 0.5 * gamma * mach * mach
    coefficients: dict[str, list[float]] = {"upper": [], "lower": []}
    panels = []
    largest = 0.0
    for surface in ("upper", "lower"):
        for index, (a, b, slope) in enumerate(_panels(section, surface)):
            theta = _inclination(surface, slope, alpha)
            largest = max(largest, abs(theta))
            if order == 1:
                cp = ackeret_pressure_coefficient(theta, mach)
            else:
                cp = busemann_pressure_coefficient(theta, mach, gamma)
            coefficients[surface].append(cp)
            panels.append(PanelSolution(surface, index, a, b, degrees(theta),
                                        "linear", mach, 1.0 + dynamic * cp, cp))
    reference = section.upper[0] if moment_reference is None else moment_reference
    lift, drag, cn, ca, cm = integrate_section_forces(
        section, coefficients, alpha, (base_pressure_ratio - 1.0) / dynamic, reference
    )
    return SectionCoefficients(
        method="ackeret" if order == 1 else "busemann",
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
        warnings=supersonic_validity_warnings(mach, largest),
    )
