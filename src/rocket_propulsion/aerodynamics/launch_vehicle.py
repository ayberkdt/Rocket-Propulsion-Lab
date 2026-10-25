"""Launch-vehicle ascent aerodynamics: dynamic pressure, drag build-up and loads.

Contents
--------
* Dynamic pressure from Mach number and static pressure, ``q = gamma p M^2 / 2``
  (from ``q = rho u^2 / 2`` [1]_, ``a = sqrt(gamma R T)`` [2]_ and
  ``p = rho R T`` [3]_), on the U.S. Standard Atmosphere 1976 pressure [4]_
  already implemented in :mod:`rocket_propulsion.compressible.off_design_nozzle`,
  and the maximum-q point of a sampled trajectory.
* Barrowman's component drag build-up for a slender body [5]_ (chapter 4):
  skin friction with roughness and compressibility, the fineness-ratio
  body-friction correction, subsonic forebody pressure drag and base drag
  (subsonic and supersonic).
* Aerodynamic loads: dimensional forces from coefficients, shear and bending
  moment from a sectional normal-force (line-load) distribution, the NASA
  SP-8035 discrete-gust profile and its rule for combining random load
  contributions [6]_, and the wind-induced angle of attack.

What is not here
----------------
Supersonic forebody wave drag (Barrowman uses second-order shock-expansion
theory, equation 4-47), supersonic fin aerodynamics (Barrowman appendix A),
fin drag, buffet, plume and multi-body (booster separation) effects. Those
need methods or data that were not implemented; see
``docs/LAUNCH_VEHICLE_AERODYNAMICS.md``.

References
----------
.. [1] NASA-GRC-DYNAMIC-PRESSURE. NASA Glenn Research Center, "Dynamic
       Pressure", Beginner's Guide to Aeronautics.
.. [2] NASA-GRC-SPEED-OF-SOUND. NASA Glenn Research Center, "Speed of Sound",
       Beginner's Guide to Aeronautics.
.. [3] NASA-GAS. NASA Glenn Research Center, "Equation of State (Ideal Gas)".
.. [4] NOAA/NASA/USAF, "U.S. Standard Atmosphere, 1976", NOAA-S/T 76-1562 /
       NASA-TM-X-74335. https://ntrs.nasa.gov/citations/19770009539
.. [5] BARROWMAN-1967. J. S. Barrowman, "The Practical Calculation of the
       Aerodynamic Characteristics of Slender Finned Vehicles", M.S. thesis,
       1967; NTRS 20010047838.
.. [6] NASA-SP-8035. NASA, "Wind Loads During Ascent", NASA SP-8035, 1970;
       NTRS 19700027637.
"""

from __future__ import annotations

import itertools
from collections.abc import Sequence
from dataclasses import dataclass
from math import atan2, cos, degrees, exp, isfinite, log10, pi, radians, sin, sqrt

from rocket_propulsion.compressible.off_design_nozzle import standard_atmosphere_pressure_pa
from rocket_propulsion.core.errors import DomainError, InputError
from rocket_propulsion.core.metadata import WarningMessage

# Barrowman (1967), section 4.11: transition Reynolds number "about 5 x 10^5".
BARROWMAN_TRANSITION_REYNOLDS = 5.0e5

# NASA SP-8035 (1970), Appendix A.2: gust geometry constants.
SP8035_GUST_AMPLITUDE_M_S = 7.65
SP8035_GUST_HALF_WAVELENGTH_M = 30.0
SP8035_A1_M = 0.9215
SP8035_A2_M = 0.9137
SP8035_THICKNESS_RANGE_M = (30.0, 275.0)


def _number(value: float, field: str) -> float:
    """Return ``value`` as a finite float.

    Raises
    ------
    InputError
        If the value is not a finite real number.

    References
    ----------
    .. [1] BARROWMAN-1967. J. S. Barrowman, thesis, 1967, NTRS 20010047838.
    """

    if not isinstance(value, (int, float)) or isinstance(value, bool):
        raise InputError(f"{field} must be a number.", field=field)
    number = float(value)
    if not isfinite(number):
        raise InputError(f"{field} must be finite.", field=field)
    return number


def _positive(value: float, field: str) -> float:
    """Return ``value`` if finite and strictly positive.

    References
    ----------
    .. [1] BARROWMAN-1967. J. S. Barrowman, thesis, 1967, NTRS 20010047838.
    """

    number = _number(value, field)
    if number <= 0.0:
        raise InputError(f"{field} must be positive.", field=field)
    return number


def _non_negative(value: float, field: str) -> float:
    """Return ``value`` if finite and not negative.

    References
    ----------
    .. [1] BARROWMAN-1967. J. S. Barrowman, thesis, 1967, NTRS 20010047838.
    """

    number = _number(value, field)
    if number < 0.0:
        raise InputError(f"{field} must be non-negative.", field=field)
    return number


# ---------------------------------------------------------------------------
# Dynamic pressure
# ---------------------------------------------------------------------------


def dynamic_pressure_from_mach(*, static_pressure_pa: float, mach: float, gamma: float = 1.4) -> float:
    """Return the dynamic pressure ``q = gamma p M^2 / 2`` of a perfect gas, Pa.

    Parameters
    ----------
    static_pressure_pa : float
        Free-stream static pressure, Pa (non-negative).
    mach : float
        Free-stream Mach number (non-negative).
    gamma : float, default 1.4
        Ratio of specific heats, greater than one.

    Returns
    -------
    float
        Dynamic pressure, Pa.

    Notes
    -----
    ``q = rho u^2 / 2`` [1]_; with ``u = M a``, ``a^2 = gamma R T`` [2]_ and
    ``p = rho R T`` [3]_ this is ``q = gamma p M^2 / 2``.

    References
    ----------
    .. [1] NASA-GRC-DYNAMIC-PRESSURE. NASA Glenn Research Center, "Dynamic
           Pressure": ``q = rho u^2 / 2``.
    .. [2] NASA-GRC-SPEED-OF-SOUND. NASA Glenn Research Center, "Speed of
           Sound": ``a = sqrt(gamma R T)``.
    .. [3] NASA-GAS. NASA Glenn Research Center, "Equation of State (Ideal
           Gas)": ``p = rho R T``.
    """

    pressure = _non_negative(static_pressure_pa, "static_pressure_pa")
    mach_number = _non_negative(mach, "mach")
    ratio = _number(gamma, "gamma")
    if ratio <= 1.0:
        raise InputError("gamma must exceed one.", field="gamma")
    return 0.5 * ratio * pressure * mach_number * mach_number


def dynamic_pressure_at_altitude(*, geopotential_altitude_m: float, mach: float, gamma: float = 1.4) -> float:
    """Return ``q`` at a U.S. Standard Atmosphere 1976 altitude and Mach number.

    Parameters
    ----------
    geopotential_altitude_m : float
        Geopotential altitude, 0 to 84 852 m'.
    mach : float
        Mach number relative to the air.
    gamma : float, default 1.4
        Ratio of specific heats.

    Returns
    -------
    float
        Dynamic pressure, Pa.

    References
    ----------
    .. [1] NOAA/NASA/USAF, "U.S. Standard Atmosphere, 1976", NOAA-S/T 76-1562 /
           NASA-TM-X-74335. https://ntrs.nasa.gov/citations/19770009539
    .. [2] NASA-GRC-DYNAMIC-PRESSURE. NASA Glenn Research Center, "Dynamic
           Pressure".
    """

    pressure = standard_atmosphere_pressure_pa(_number(geopotential_altitude_m, "geopotential_altitude_m"))
    return dynamic_pressure_from_mach(static_pressure_pa=pressure, mach=mach, gamma=gamma)


@dataclass(frozen=True, slots=True)
class MaxQPoint:
    """The sampled trajectory point of maximum dynamic pressure.

    Attributes
    ----------
    index : int
        Index of the point in the input sequence.
    time_s : float
        Time of the point, s.
    dynamic_pressure_pa : float
        Maximum sampled ``q``, Pa.
    dynamic_pressures_pa : tuple of float
        ``q`` at every input point, Pa.

    References
    ----------
    .. [1] NASA-SP-8035. NASA, "Wind Loads During Ascent", 1970, section 1
           (maximum loads in the high-dynamic-pressure region).
    """

    index: int
    time_s: float
    dynamic_pressure_pa: float
    dynamic_pressures_pa: tuple[float, ...]


def maximum_dynamic_pressure(
    *,
    times_s: Sequence[float],
    geopotential_altitudes_m: Sequence[float],
    machs: Sequence[float],
    gamma: float = 1.4,
) -> MaxQPoint:
    """Return the maximum of ``q`` over sampled trajectory points.

    Parameters
    ----------
    times_s, geopotential_altitudes_m, machs : sequence of float
        Equal-length samples of a trajectory.
    gamma : float, default 1.4
        Ratio of specific heats.

    Returns
    -------
    MaxQPoint
        The largest sampled ``q``. No interpolation between samples is made.

    Notes
    -----
    SP-8035 [1]_ notes that ascent wind loads are usually largest in the
    high-dynamic-pressure region (6-15 km); this helper only locates the
    sampled maximum.

    References
    ----------
    .. [1] NASA-SP-8035. NASA, "Wind Loads During Ascent", 1970, section 1.
    .. [2] NOAA/NASA/USAF, "U.S. Standard Atmosphere, 1976", NOAA-S/T 76-1562 /
           NASA-TM-X-74335. https://ntrs.nasa.gov/citations/19770009539
    """

    times = tuple(_number(value, "times_s") for value in times_s)
    altitudes = tuple(geopotential_altitudes_m)
    mach_values = tuple(machs)
    if not times or len(times) != len(altitudes) or len(times) != len(mach_values):
        raise InputError("Trajectory samples must be non-empty and of equal length.", field="times_s")
    q_values = tuple(
        dynamic_pressure_at_altitude(geopotential_altitude_m=h, mach=m, gamma=gamma)
        for h, m in zip(altitudes, mach_values, strict=True)
    )
    index = max(range(len(q_values)), key=q_values.__getitem__)
    return MaxQPoint(index, times[index], q_values[index], q_values)


# ---------------------------------------------------------------------------
# Barrowman drag build-up (chapter 4)
# ---------------------------------------------------------------------------


@dataclass(frozen=True, slots=True)
class SkinFrictionResult:
    """Barrowman flat-plate skin-friction coefficient [1]_.

    Attributes
    ----------
    reynolds_number : float
        Reynolds number based on the surface length.
    mach : float
        Free-stream Mach number.
    regime : str
        ``"laminar"``, ``"transitional"``, ``"turbulent_smooth"`` or
        ``"turbulent_rough"`` (incompressible branch selected).
    incompressible_coefficient : float
        ``C_f`` before the compressibility correction.
    compressible_coefficient : float
        ``C_fc`` after the correction of section 4.12 or 4.13.
    warnings : tuple of WarningMessage
        Validity notes.

    References
    ----------
    .. [1] BARROWMAN-1967. J. S. Barrowman, thesis, 1967, NTRS 20010047838,
           sections 4.10-4.13.
    """

    reynolds_number: float
    mach: float
    regime: str
    incompressible_coefficient: float
    compressible_coefficient: float
    warnings: tuple[WarningMessage, ...]


def hama_turbulent_skin_friction(reynolds_number: float) -> float:
    """Return Hama's smooth turbulent ``C_f = 1 / (3.46 log10 R - 5.6)^2``.

    Barrowman [1]_ (equation 4-4) states that it agrees with the
    Schoenherr relation (equation 4-3) within +-2%.

    References
    ----------
    .. [1] BARROWMAN-1967. J. S. Barrowman, thesis, 1967, NTRS 20010047838,
           equation 4-4.
    """

    reynolds = _positive(reynolds_number, "reynolds_number")
    denominator = 3.46 * log10(reynolds) - 5.6
    if denominator <= 0.0:
        raise DomainError("Reynolds number too low for the turbulent relation.", field="reynolds_number")
    return 1.0 / (denominator * denominator)


def schoenherr_turbulent_skin_friction(reynolds_number: float) -> float:
    """Solve Schoenherr's ``log10(R C_f) = 0.242 / sqrt(C_f)`` for ``C_f``.

    Parameters
    ----------
    reynolds_number : float
        Reynolds number, at least ``1e4``.

    Returns
    -------
    float
        Turbulent flat-plate ``C_f``.

    Notes
    -----
    Solved by bisection on ``C_f`` in ``(1e-5, 0.05)``; the residual
    ``log10(R C_f) - 0.242 / sqrt(C_f)`` increases monotonically with
    ``C_f``. Barrowman [1]_ quotes the relation as equation 4-3.

    References
    ----------
    .. [1] BARROWMAN-1967. J. S. Barrowman, thesis, 1967, NTRS 20010047838,
           equation 4-3.
    """

    reynolds = _positive(reynolds_number, "reynolds_number")
    if reynolds < 1.0e4:
        raise DomainError("Schoenherr relation used here for R >= 1e4.", field="reynolds_number")
    lower, upper = 1.0e-5, 0.05

    def residual(cf: float) -> float:
        return log10(reynolds * cf) - 0.242 / sqrt(cf)

    if residual(lower) > 0.0 or residual(upper) < 0.0:
        raise DomainError("Schoenherr relation has no root in the search bracket.")
    for _ in range(200):
        middle = 0.5 * (lower + upper)
        if residual(middle) > 0.0:
            upper = middle
        else:
            lower = middle
    return 0.5 * (lower + upper)


def barrowman_skin_friction(
    *,
    reynolds_number: float,
    mach: float,
    roughness_to_length: float = 0.0,
    wall: str = "adiabatic",
    fully_turbulent: bool = False,
) -> SkinFrictionResult:
    """Return Barrowman's compressible flat-plate skin-friction coefficient.

    Parameters
    ----------
    reynolds_number : float
        Reynolds number based on the surface length (positive).
    mach : float
        Free-stream Mach number (non-negative).
    roughness_to_length : float, default 0.0
        Equivalent sand-roughness height over length, ``R_s / L_r``. Zero
        selects the smooth-surface relations.
    wall : {"adiabatic", "cooled"}, default "adiabatic"
        Supersonic smooth turbulent wall condition of equation 4-12
        (``k = 0.15`` without heat transfer, ``k = 0.0512`` cooled wall).
    fully_turbulent : bool, default False
        Use Hama's turbulent relation (4-4) without the transitional
        correction for ``R >= 5e5``.

    Returns
    -------
    SkinFrictionResult
        Incompressible and compressible ``C_f``.

    Notes
    -----
    Incompressible (section 4.11): laminar ``1.328 / sqrt(R)`` (4-2) below
    ``R = 5e5``; above it, the transitional ``1/(3.46 log10 R - 5.6)^2 -
    1700/R`` (4-6), or Hama (4-4) if ``fully_turbulent``. With roughness the
    flow is taken as turbulent (Hama) up to the critical Reynolds number
    ``51 (R_s/L)^-1.039`` (4-7) and ``0.032 (R_s/L)^0.2`` (4-8) above it.
    Compressibility: laminar unchanged subsonic and ``/(1 + 0.045 M^2)^0.25``
    supersonic (4-11); smooth turbulent ``(1 - 0.09 M^2)`` subsonic (4-9) and
    ``/(1 + k M^2)^0.58`` supersonic (4-12); rough ``(1 - 0.12 M^2)``
    subsonic (4-10) and ``/(1 + 0.18 M^2)`` supersonic (4-13), never below
    the smooth value at the same Mach number, as Barrowman requires. The
    subsonic forms are used for ``M < 1``; the branches differ by about 1%
    at ``M = 1``. The transitional branch uses the smooth turbulent
    correction (Barrowman states corrections only for laminar and turbulent
    flow).

    References
    ----------
    .. [1] BARROWMAN-1967. J. S. Barrowman, thesis, 1967, NTRS 20010047838,
           equations 4-2 to 4-13, printed pp. 43-48.
    """

    reynolds = _positive(reynolds_number, "reynolds_number")
    mach_number = _non_negative(mach, "mach")
    roughness = _non_negative(roughness_to_length, "roughness_to_length")
    if wall not in ("adiabatic", "cooled"):
        raise InputError("wall must be 'adiabatic' or 'cooled'.", field="wall")
    k_wall = 0.15 if wall == "adiabatic" else 0.0512
    warnings: list[WarningMessage] = []
    supersonic = mach_number >= 1.0

    def smooth_turbulent_correction(cf: float) -> float:
        if supersonic:
            return cf / (1.0 + k_wall * mach_number * mach_number) ** 0.58
        return cf * (1.0 - 0.09 * mach_number * mach_number)

    if roughness > 0.0:
        critical = 51.0 * roughness ** (-1.039)
        smooth = hama_turbulent_skin_friction(reynolds)
        smooth_c = smooth_turbulent_correction(smooth)
        if reynolds < critical:
            regime = "turbulent_smooth"
            cf = smooth
            cfc = smooth_c
        else:
            regime = "turbulent_rough"
            cf = 0.032 * roughness**0.2
            if supersonic:
                cfc = cf / (1.0 + 0.18 * mach_number * mach_number)
            else:
                cfc = cf * (1.0 - 0.12 * mach_number * mach_number)
            cfc = max(cfc, smooth_c)
    elif reynolds < BARROWMAN_TRANSITION_REYNOLDS:
        regime = "laminar"
        cf = 1.328 / sqrt(reynolds)
        cfc = cf / (1.0 + 0.045 * mach_number * mach_number) ** 0.25 if supersonic else cf
    else:
        cf = hama_turbulent_skin_friction(reynolds)
        if fully_turbulent:
            regime = "turbulent_smooth"
        else:
            regime = "transitional"
            cf -= 1700.0 / reynolds
        cfc = smooth_turbulent_correction(cf)
    if cfc <= 0.0:
        raise DomainError("Skin-friction coefficient is not positive at this condition.")
    if not supersonic and mach_number > 0.9:
        warnings.append(
            WarningMessage(
                "skin_friction_transonic",
                "Barrowman's subsonic compressibility factors are used up to M < 1.",
                field="mach",
            )
        )
    return SkinFrictionResult(reynolds, mach_number, regime, cf, cfc, tuple(warnings))


def body_friction_drag_coefficient(
    *,
    compressible_skin_friction: float,
    body_fineness_ratio: float,
    wetted_area_m2: float,
    reference_area_m2: float,
) -> float:
    """Return Barrowman's body skin-friction drag coefficient.

    ``C_Df_B = [1 + 0.5 / f_B] C_fc A_wB / A_r`` (equation 4-16), the
    flat-plate value corrected for the body not being a flat plate [1]_.

    References
    ----------
    .. [1] BARROWMAN-1967. J. S. Barrowman, thesis, 1967, NTRS 20010047838,
           equation 4-16, printed p. 48.
    """

    cfc = _positive(compressible_skin_friction, "compressible_skin_friction")
    fineness = _positive(body_fineness_ratio, "body_fineness_ratio")
    wetted = _positive(wetted_area_m2, "wetted_area_m2")
    reference = _positive(reference_area_m2, "reference_area_m2")
    return (1.0 + 0.5 / fineness) * cfc * wetted / reference


def subsonic_forebody_pressure_drag_coefficient(
    *,
    mach: float,
    nose_shape: str,
    nose_fineness_ratio: float,
    body_fineness_ratio: float,
    wetted_area_m2: float,
    reference_area_m2: float,
    compressible_skin_friction: float,
    sonic_skin_friction: float,
) -> float:
    """Return Barrowman's subsonic compressible forebody pressure drag.

    Parameters
    ----------
    mach : float
        Free-stream Mach number, ``0 <= M < 1``.
    nose_shape : {"cone", "tangent_ogive"}
        Nose shape for the sonic drag fit.
    nose_fineness_ratio : float
        ``f_N`` = nose length / nose base diameter, within ``1 <= f_N <= 10``
        (the range of the fits 4-42 and 4-43).
    body_fineness_ratio : float
        ``f_B`` = body length / diameter.
    wetted_area_m2 : float
        Body wetted area ``A_wB``, m^2.
    reference_area_m2 : float
        Reference area ``A_r``, m^2.
    compressible_skin_friction : float
        ``C_fc`` at the flight condition.
    sonic_skin_friction : float
        ``C*_fc``, the compressible skin friction at Mach one.

    Returns
    -------
    float
        ``C_Dp = 6 A_wB C_fc / (f_B^3 A_r (K - M^2)^0.6)`` (equation 4-41)
        with ``K = 1 + [6.82 A_wB C*_fc (f_N + 0.7)^1.29 / (f_B^3 A_r)]^(5/3)``
        for a cone (4-45) or ``(f_N + 1)^2.22`` for an ogive (4-46); ``K`` is
        chosen so that ``C_Dp`` equals the Mach-one value
        ``0.88 / (f_N + 0.7)^1.29`` (cone, 4-42) or ``0.88 / (f_N + 1)^2.22``
        (ogive, 4-43); the printed 6.82 is ``6 / 0.88 = 6.818`` rounded, so
        the match at Mach one is to about 0.03%.

    References
    ----------
    .. [1] BARROWMAN-1967. J. S. Barrowman, thesis, 1967, NTRS 20010047838,
           equations 4-40 to 4-46, printed pp. 57-59.
    """

    mach_number = _non_negative(mach, "mach")
    if mach_number >= 1.0:
        raise DomainError("Equation 4-41 is subsonic; supersonic wave drag is not implemented.", field="mach")
    if nose_shape not in ("cone", "tangent_ogive"):
        raise InputError("nose_shape must be 'cone' or 'tangent_ogive'.", field="nose_shape")
    f_n = _positive(nose_fineness_ratio, "nose_fineness_ratio")
    if not 1.0 <= f_n <= 10.0:
        raise DomainError("Barrowman's sonic forebody-drag fits cover 1 <= f_N <= 10.", field="nose_fineness_ratio")
    f_b = _positive(body_fineness_ratio, "body_fineness_ratio")
    wetted = _positive(wetted_area_m2, "wetted_area_m2")
    reference = _positive(reference_area_m2, "reference_area_m2")
    cfc = _positive(compressible_skin_friction, "compressible_skin_friction")
    cfc_star = _positive(sonic_skin_friction, "sonic_skin_friction")
    shape_term = (f_n + 0.7) ** 1.29 if nose_shape == "cone" else (f_n + 1.0) ** 2.22
    k_constant = 1.0 + (6.82 * wetted * cfc_star * shape_term / (f_b**3 * reference)) ** (5.0 / 3.0)
    return 6.0 * wetted * cfc / (f_b**3 * reference * (k_constant - mach_number * mach_number) ** 0.6)


@dataclass(frozen=True, slots=True)
class BaseDragResult:
    """Barrowman body base-drag coefficient [1]_.

    Attributes
    ----------
    mach : float
        Free-stream Mach number.
    base_drag_coefficient : float
        ``C_DB`` referred to the reference area.
    sonic_base_drag_coefficient : float
        ``C*_DB`` of equation 4-51.
    critical_mach : float
        ``M_cr = 0.892 / sqrt(C*_DB)`` (equation 4-56).
    branch : str
        ``"subsonic"``, ``"supersonic_exponential"`` or
        ``"supersonic_inverse_square"``.
    warnings : tuple of WarningMessage
        Validity notes.

    References
    ----------
    .. [1] BARROWMAN-1967. J. S. Barrowman, thesis, 1967, NTRS 20010047838,
           equations 4-48 to 4-56.
    """

    mach: float
    base_drag_coefficient: float
    sonic_base_drag_coefficient: float
    critical_mach: float
    branch: str
    warnings: tuple[WarningMessage, ...]


def barrowman_base_drag(
    *,
    mach: float,
    base_area_m2: float,
    reference_area_m2: float,
    wetted_area_m2: float,
    compressible_skin_friction: float,
    sonic_skin_friction: float,
    fin_trailing_edge_thickness_to_root_chord: float = 0.0,
) -> BaseDragResult:
    """Return Barrowman's body base-drag coefficient at any Mach number.

    Parameters
    ----------
    mach : float
        Free-stream Mach number (non-negative).
    base_area_m2, reference_area_m2, wetted_area_m2 : float
        ``A_B``, ``A_r`` and the body wetted area ``A_w``, m^2.
    compressible_skin_friction : float
        ``C_fc`` of the body at the flight condition (subsonic branch).
    sonic_skin_friction : float
        ``C*_fc`` at Mach one (subsonic branch constant ``K``).
    fin_trailing_edge_thickness_to_root_chord : float, default 0.0
        ``h_r / c_r`` of fins near the base (0 for none or zero thickness).

    Returns
    -------
    BaseDragResult
        Coefficient and branch.

    Notes
    -----
    Sonic value (4-51): ``C*_DB = (A_B/A_r)(0.185 + 1.15 h_r/c_r)``.
    Subsonic (4-50 with the 0.029 of equations 4-48, 4-49 and 4-52)::

        C_DB = 0.029 (A_B/A_r) / sqrt(C_fc (A_w/A_B) (K - M^2))
        K = 1 + (0.029 A_B/A_r)^2 / (C*_fc (A_w/A_B) (C*_DB)^2)   (4-52)

    Supersonic: ``C_DB = C*_DB [0.88 + 0.12 exp(-3.58 (M - 1))]`` for
    ``1 <= M <= M_cr`` (4-54) and ``0.7 (A_B/A_r) / M^2`` above (4-55), with
    ``M_cr = 0.892 / sqrt(C*_DB)`` (4-56).

    The thesis prints 0.29 in equation 4-50 but 0.029 in 4-48, 4-49 and 4-52,
    and equation 4-53 as printed lacks the square that its derivation from
    4-52 requires (``0.185/0.029 = 6.38``, ``1.15/0.029 = 39.7``); this
    function uses 0.029 and equation 4-52 directly. Equation 4-56 follows
    from equating 0.88 C*_DB and 0.7/M^2 with ``A_B = A_r``; the two
    supersonic branches therefore do not meet exactly, and a warning reports
    the jump at ``M_cr``.

    References
    ----------
    .. [1] BARROWMAN-1967. J. S. Barrowman, thesis, 1967, NTRS 20010047838,
           equations 4-48 to 4-56, printed pp. 59-62.
    """

    mach_number = _non_negative(mach, "mach")
    base = _positive(base_area_m2, "base_area_m2")
    reference = _positive(reference_area_m2, "reference_area_m2")
    wetted = _positive(wetted_area_m2, "wetted_area_m2")
    cfc = _positive(compressible_skin_friction, "compressible_skin_friction")
    cfc_star = _positive(sonic_skin_friction, "sonic_skin_friction")
    thickness = _non_negative(fin_trailing_edge_thickness_to_root_chord, "fin_trailing_edge_thickness_to_root_chord")
    sonic = (base / reference) * (0.185 + 1.15 * thickness)
    critical = 0.892 / sqrt(sonic)
    warnings: list[WarningMessage] = []
    if mach_number < 1.0:
        k_constant = 1.0 + (0.029 * base / reference) ** 2 / (cfc_star * (wetted / base) * sonic * sonic)
        value = 0.029 * (base / reference) / sqrt(cfc * (wetted / base) * (k_constant - mach_number**2))
        branch = "subsonic"
    else:
        exponential_at_critical = sonic * (0.88 + 0.12 * exp(-3.58 * (critical - 1.0)))
        inverse_at_critical = 0.7 * (base / reference) / (critical * critical)
        if mach_number <= critical:
            value = sonic * (0.88 + 0.12 * exp(-3.58 * (mach_number - 1.0)))
            branch = "supersonic_exponential"
        else:
            value = 0.7 * (base / reference) / (mach_number * mach_number)
            branch = "supersonic_inverse_square"
        jump = abs(exponential_at_critical - inverse_at_critical)
        if jump > 1.0e-3 * exponential_at_critical:
            warnings.append(
                WarningMessage(
                    "base_drag_branch_jump",
                    f"Barrowman's supersonic base-drag branches differ by {jump:.3g} at M_cr = {critical:.3g}.",
                    severity="info",
                    field="mach",
                )
            )
        if critical < 1.0:
            warnings.append(
                WarningMessage(
                    "base_drag_critical_below_one",
                    "M_cr < 1: the inverse-square branch applies from Mach one.",
                    field="fin_trailing_edge_thickness_to_root_chord",
                )
            )
    return BaseDragResult(mach_number, value, sonic, critical, branch, tuple(warnings))


# ---------------------------------------------------------------------------
# Loads
# ---------------------------------------------------------------------------


@dataclass(frozen=True, slots=True)
class AerodynamicLoads:
    """Body-axis aerodynamic forces and moments from coefficients.

    Attributes
    ----------
    axial_force_n, side_force_n, normal_force_n : float
        ``q S C_A``, ``q S C_Y``, ``q S C_N``, N.
    rolling_moment_n_m, pitching_moment_n_m, yawing_moment_n_m : float
        ``q S l C_l``, ``q S l C_m``, ``q S l C_n``, N m.

    References
    ----------
    .. [1] BARROWMAN-1967. J. S. Barrowman, thesis, 1967, NTRS 20010047838,
           equations 3-60 and 3-78 (coefficient definitions).
    .. [2] SHEA-2020-SLS-ASCENT-DB. P. R. Shea et al., "Ascent Aerodynamic
           Force and Moment Database Development for the Space Launch System",
           NTRS 20200002382 (six body-axis coefficients CA, CY, CN, CLL, CLM,
           CLN).
    """

    axial_force_n: float
    side_force_n: float
    normal_force_n: float
    rolling_moment_n_m: float
    pitching_moment_n_m: float
    yawing_moment_n_m: float


def loads_from_coefficients(
    *,
    dynamic_pressure_pa: float,
    reference_area_m2: float,
    reference_length_m: float,
    axial: float = 0.0,
    side: float = 0.0,
    normal: float = 0.0,
    rolling: float = 0.0,
    pitching: float = 0.0,
    yawing: float = 0.0,
) -> AerodynamicLoads:
    """Dimensionalise six body-axis coefficients with ``q``, ``S`` and ``l``.

    Forces are ``q S C`` and moments ``q S l C``, the definitions
    ``C_N = n / (q A_r)`` and ``C_m = m / (q A_r L_r)`` of Barrowman [1]_
    (equations 3-60 and 3-78), applied to the six body-axis coefficients that
    the SLS ascent database tabulates [2]_. Sign conventions are those of the
    coefficients supplied.

    References
    ----------
    .. [1] BARROWMAN-1967. J. S. Barrowman, thesis, 1967, NTRS 20010047838,
           equations 3-60 and 3-78.
    .. [2] SHEA-2020-SLS-ASCENT-DB. P. R. Shea et al., NTRS 20200002382,
           nomenclature.
    """

    q = _non_negative(dynamic_pressure_pa, "dynamic_pressure_pa")
    area = _positive(reference_area_m2, "reference_area_m2")
    length = _positive(reference_length_m, "reference_length_m")
    force = q * area
    moment = force * length
    return AerodynamicLoads(
        axial_force_n=force * _number(axial, "axial"),
        side_force_n=force * _number(side, "side"),
        normal_force_n=force * _number(normal, "normal"),
        rolling_moment_n_m=moment * _number(rolling, "rolling"),
        pitching_moment_n_m=moment * _number(pitching, "pitching"),
        yawing_moment_n_m=moment * _number(yawing, "yawing"),
    )


@dataclass(frozen=True, slots=True)
class ShearBendingDistribution:
    """Running shear and bending from a sectional normal-force distribution.

    Attributes
    ----------
    stations_m : tuple of float
        Axial stations from the nose, m.
    line_load_n_m : tuple of float
        Normal force per unit length at the stations, N/m.
    shear_n : tuple of float
        ``V(x) = integral_0^x n(xi) dxi``, N.
    bending_moment_n_m : tuple of float
        ``M(x) = integral_0^x (x - xi) n(xi) dxi``, N m.
    total_normal_force_n : float
        ``V`` at the last station, N.
    center_of_pressure_m : float or None
        ``x_last - M(x_last) / V(x_last)``, m; ``None`` if ``V`` is zero.

    References
    ----------
    .. [1] HEMSCH-2016-SLS-CFD-WT. M. J. Hemsch, "Global Comparison of CFD and
           Wind-Tunnel-Derived Force and Moment Databases for the Space Launch
           System", NTRS 20160007695 (line loads integrated to forces and
           moments).
    """

    stations_m: tuple[float, ...]
    line_load_n_m: tuple[float, ...]
    shear_n: tuple[float, ...]
    bending_moment_n_m: tuple[float, ...]
    total_normal_force_n: float
    center_of_pressure_m: float | None


def integrate_line_load(
    *, stations_m: Sequence[float], line_load_n_m: Sequence[float]
) -> ShearBendingDistribution:
    """Integrate a sectional normal-force distribution to shear and bending.

    Parameters
    ----------
    stations_m : sequence of float
        Strictly increasing stations from the nose tip, m (at least two).
    line_load_n_m : sequence of float
        Normal force per unit length at the stations, N/m (for example from
        :func:`rocket_propulsion.aerodynamics.slender_body.crossflow_line_load`
        or a CFD line-load database).

    Returns
    -------
    ShearBendingDistribution
        Free-free aerodynamic shear and bending, integrated from the nose with
        the trapezoidal rule for ``V`` and for ``M`` (``dM/dx = V``).

    Notes
    -----
    These are the aerodynamic contributions only. SP-8035 [1]_ lists the
    other ascent load sources (control forces, inertial forces from the
    vehicle's lateral and rotational motion, gust-excited elastic and slosh
    response); a structural load needs the mass distribution and a flight
    simulation and is outside this function. The SLS programme integrates CFD
    line loads to integrated forces and moments for comparison with its
    force-and-moment database [2]_.

    References
    ----------
    .. [1] NASA-SP-8035. NASA, "Wind Loads During Ascent", 1970, section 1.
    .. [2] HEMSCH-2016-SLS-CFD-WT. M. J. Hemsch, NTRS 20160007695, abstract.
    """

    stations = tuple(_number(value, "stations_m") for value in stations_m)
    loads = tuple(_number(value, "line_load_n_m") for value in line_load_n_m)
    if len(stations) < 2 or len(stations) != len(loads):
        raise InputError("Need at least two stations and one load per station.", field="stations_m")
    if any(b <= a for a, b in itertools.pairwise(stations)):
        raise InputError("Stations must be strictly increasing.", field="stations_m")
    shear = [0.0]
    bending = [0.0]
    for index in range(1, len(stations)):
        width = stations[index] - stations[index - 1]
        shear.append(shear[-1] + 0.5 * width * (loads[index] + loads[index - 1]))
        bending.append(bending[-1] + 0.5 * width * (shear[-1] + shear[-2]))
    total = shear[-1]
    center = None if total == 0.0 else stations[-1] - bending[-1] / total
    return ShearBendingDistribution(stations, loads, tuple(shear), tuple(bending), total, center)


def sp8035_discrete_gust_speed(*, delta_altitude_m: float, thickness_m: float) -> float:
    """Return the SP-8035 idealised one-minus-cosine gust speed, m/s.

    Parameters
    ----------
    delta_altitude_m : float
        ``Delta H``, altitude above point O where the shear build-up meets
        the wind-speed envelope, m (non-negative).
    thickness_m : float
        Gust thickness ``th`` between the inflection points of build-up and
        tail-off, 30 to 275 m.

    Returns
    -------
    float
        ``Delta W_G`` per Appendix A.2, equations (2)-(5)::

            a2 <= dH <= 30 - a1        3.825 [1 - cos(pi (dH + a1) / 30)]
            30 - a1 <= dH <= th - a1   7.65
            th - a1 <= dH <= th+30-a1  3.825 [1 - cos(pi (dH + 30 + a1 - th)/30)]
            th + 30 - a1 <= dH         0

        with ``a1 = 0.9215 m`` and ``a2 = 0.9137 m``. For ``0 <= dH <= a2``
        (equation 1, the linear extension of the shear build-up) the slope is
        the value that makes the line tangent to equation (2) at ``a2``,
        ``3.825 (pi/30) sin(pi (a2 + a1)/30) = 0.0765 (m/s)/m``. With the
        rounded ``a1`` and ``a2`` of the monograph the line and equation (2)
        then differ by about 5e-4 m/s at ``a2``.

    Notes
    -----
    The monograph prints equation (1) as ``0.765 dH``. That is ten times the
    slope of equation (2) at ``a2`` and does not meet equation (2) there,
    although the text defines ``a2`` as the tangent point; the tangent slope
    is therefore used. The 7.65 m/s amplitude is a 9 m/s embedded-gust
    maximum times 0.85 for less-than-unity gust/shear correlation.

    References
    ----------
    .. [1] NASA-SP-8035. NASA, "Wind Loads During Ascent", NASA SP-8035, 1970,
           Appendix A.2, printed p. 24, equations (1)-(5).
    """

    delta = _non_negative(delta_altitude_m, "delta_altitude_m")
    thickness = _number(thickness_m, "thickness_m")
    low, high = SP8035_THICKNESS_RANGE_M
    if not low <= thickness <= high:
        raise DomainError("SP-8035 gust thickness must lie between 30 and 275 m.", field="thickness_m")
    half_amplitude = 0.5 * SP8035_GUST_AMPLITUDE_M_S
    wavelength = SP8035_GUST_HALF_WAVELENGTH_M
    a1 = SP8035_A1_M
    a2 = SP8035_A2_M
    if delta <= a2:
        slope = half_amplitude * (pi / wavelength) * sin(pi * (a2 + a1) / wavelength)
        return slope * delta
    if delta <= wavelength - a1:
        return half_amplitude * (1.0 - cos(pi * (delta + a1) / wavelength))
    if delta <= thickness - a1:
        return SP8035_GUST_AMPLITUDE_M_S
    if delta <= thickness + wavelength - a1:
        return half_amplitude * (1.0 - cos(pi * (delta + wavelength + a1 - thickness) / wavelength))
    return 0.0


@dataclass(frozen=True, slots=True)
class CombinedLoad:
    """SP-8035 combination of quasi-steady and random load contributions.

    Attributes
    ----------
    quasi_steady : float
        Quasi-steady wind and wind-shear load.
    random_rss : float
        Root-sum-square of the independent random contributions.
    design_load : float
        ``quasi_steady + sign(quasi_steady) random_rss`` (magnitude-increasing
        superposition).
    lower, upper : float
        ``quasi_steady -+ random_rss``.

    References
    ----------
    .. [1] NASA-SP-8035. NASA, "Wind Loads During Ascent", 1970, section
           4.2.2 (intermediate and final design analyses), printed p. 17.
    """

    quasi_steady: float
    random_rss: float
    design_load: float
    lower: float
    upper: float


def combine_wind_loads(*, quasi_steady: float, random_components: Sequence[float]) -> CombinedLoad:
    """Combine loads as SP-8035 recommends for the final design analysis.

    Parameters
    ----------
    quasi_steady : float
        Quasi-steady-state wind and wind-shear load (any consistent unit).
    random_components : sequence of float
        Independent random contributions in the same unit: dispersion loads
        (at a consistent deviation level, e.g. 2 sigma), buffeting and gust
        loads.

    Returns
    -------
    CombinedLoad
        The random contributions are root-sum-squared and the total random
        load is superimposed on the quasi-steady load.

    Notes
    -----
    SP-8035 [1]_: the dispersion loads "are independent and random;
    therefore, the individual dispersion loads should be root sum squared
    along with the buffeting and gust loads to derive the total lateral loads
    from random effects", which are then superimposed on the quasi-steady
    wind and wind-shear loads. (Section 2.2.1.3 notes that many analyses
    instead superimpose gust loads directly, a unity-correlation assumption.)

    References
    ----------
    .. [1] NASA-SP-8035. NASA, "Wind Loads During Ascent", 1970, section
           2.2.1.3 (printed p. 9) and section 4.2.2 (printed p. 17).
    """

    base = _number(quasi_steady, "quasi_steady")
    components = tuple(_number(value, "random_components") for value in random_components)
    rss = sqrt(sum(value * value for value in components))
    sign = 1.0 if base >= 0.0 else -1.0
    return CombinedLoad(base, rss, base + sign * rss, base - rss, base + rss)


@dataclass(frozen=True, slots=True)
class WindIncidence:
    """Air-relative velocity and wind-induced angle of attack in a plane.

    Attributes
    ----------
    airspeed_m_s : float
        Magnitude of the air-relative velocity, m/s.
    angle_of_attack_deg : float
        Angle from the vehicle's ground-relative velocity vector to its
        air-relative velocity vector, degrees; positive when the air-relative
        velocity is steeper (a tail wind), negative for a head wind.
    flight_path_angle_air_deg : float
        Elevation of the air-relative velocity above the horizontal, degrees.

    References
    ----------
    .. [1] NASA-SP-8035. NASA, "Wind Loads During Ascent", 1970, section
           4.2.2, printed p. 17.
    """

    airspeed_m_s: float
    angle_of_attack_deg: float
    flight_path_angle_air_deg: float


def wind_induced_incidence(
    *, ground_speed_m_s: float, flight_path_angle_deg: float, horizontal_wind_m_s: float
) -> WindIncidence:
    """Combine the wind vector with the vehicle velocity vector (planar).

    Parameters
    ----------
    ground_speed_m_s : float
        Magnitude of the vehicle velocity relative to the ground, m/s.
    flight_path_angle_deg : float
        Elevation of that velocity above the local horizontal, degrees.
    horizontal_wind_m_s : float
        Horizontal wind component in the plane of the trajectory, positive
        when blowing in the downrange (tail-wind) direction, m/s.

    Returns
    -------
    WindIncidence
        With ``V_air = (V cos g - W, V sin g)``, the air-relative flight-path
        angle is ``g_air = atan2(V sin g, V cos g - W)`` and the incidence
        added by the wind, for a vehicle whose axis is aligned with its
        ground-relative velocity, is ``g_air - g``.

    Notes
    -----
    SP-8035 [1]_ directs that the wind-velocity vector "should be combined
    directly with the vehicle's velocity vector to determine the wind-induced
    angle of attack" at conditions where vehicle drift is negligible. The
    vector subtraction is elementary kinematics; whether the vehicle axis is
    aligned with the ground-relative velocity depends on the steering and is
    an assumption of this helper.

    References
    ----------
    .. [1] NASA-SP-8035. NASA, "Wind Loads During Ascent", 1970, section
           4.2.2, printed p. 17.
    """

    speed = _non_negative(ground_speed_m_s, "ground_speed_m_s")
    gamma = radians(_number(flight_path_angle_deg, "flight_path_angle_deg"))
    wind = _number(horizontal_wind_m_s, "horizontal_wind_m_s")
    horizontal = speed * cos(gamma) - wind
    vertical = speed * sin(gamma)
    airspeed = sqrt(horizontal * horizontal + vertical * vertical)
    if airspeed == 0.0:
        raise DomainError("Air-relative velocity is zero; incidence undefined.")
    gamma_air = atan2(vertical, horizontal)
    return WindIncidence(airspeed, degrees(gamma_air - gamma), degrees(gamma_air))
