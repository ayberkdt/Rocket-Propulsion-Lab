"""Engineering convective aeroheating: stagnation point, cones, plates, leading edges.

Included correlations, each used exactly as printed in the opened source:

* Boundary-layer stagnation-point correlations for a cylinder and a sphere
  with the Newtonian velocity gradient (Anderson Eqs. (18.65), (18.70),
  (18.82)) [1]_.
* Sutton-Graves: ``qw = K (ps/R)^0.5 (hs - hw)`` with ``K`` either from the
  computer-correlation Table II or from the general mixture relation Eq. (42)
  using Table I (``K`` in kg/(s m^1.5 atm^0.5), ``ps`` in atm) [2]_.
* Tauber's closed-form air relations: hemisphere stagnation point Eq. (40)
  (Anderson Eq. (14.62)), swept cylinder Eq. (41), laminar sharp cone
  Eq. (46) and wedge/flat plate Eq. (48), wing leading edge Eq. (49),
  turbulent flat plate Eqs. (56a)/(56b) and cone Eq. (59) [3]_.

Not included (no verified source opened in this work package): the
Fay-Riddell equation in its original form, Detra-Kemp-Riddell, and any
radiative-heating correlation (e.g. Tauber-Sutton 1991). All results are
cold-wall convective estimates; real-gas, catalytic, ablation, transition
and low-density effects are outside these formulas.

References
----------
.. [1] ANDERSON-FOA-2017. J. D. Anderson, Fundamentals of Aerodynamics, 6th
       ed., McGraw-Hill, 2017, Secs. 14.8 and 18.5, printed pp. 901-908,
       1037-1042.
.. [2] SUTTON-GRAVES-1971-TRR376. K. Sutton and R. A. Graves Jr., NASA
       TR R-376, 1971, Eqs. (33), (41)-(43), Tables I-II,
       https://ntrs.nasa.gov/citations/19720003329
.. [3] TAUBER-1989-TP2914. M. E. Tauber, NASA TP-2914, 1989, Eqs. (39)-(41),
       (45)-(49), (56)-(59), https://ntrs.nasa.gov/citations/19890017745
"""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass
from math import cos, isclose, radians, sin, sqrt
from types import MappingProxyType

from rocket_propulsion.core.constants import USSA1976_SEA_LEVEL_PRESSURE_PA
from rocket_propulsion.core.errors import DomainError, InputError
from rocket_propulsion.core.validation import real

__all__ = [
    "SUTTON_GRAVES_BASE_GASES",
    "SUTTON_GRAVES_TABLE_II",
    "SuttonGravesBaseGas",
    "newtonian_stagnation_velocity_gradient",
    "stagnation_point_heating_boundary_layer",
    "sutton_graves_coefficient",
    "sutton_graves_heating",
    "tauber_laminar_cone_heating",
    "tauber_laminar_flat_plate_heating",
    "tauber_leading_edge_heating",
    "tauber_stagnation_heating",
    "tauber_swept_cylinder_heating",
    "tauber_turbulent_cone_heating",
    "tauber_turbulent_flat_plate_heating",
]

_ATMOSPHERE_PA = USSA1976_SEA_LEVEL_PRESSURE_PA
"""Sutton-Graves state 1 atmosphere = 101.325 kN/m^2 (TR R-376, p. 2)."""


def _positive(value: object, label: str) -> float:
    number = real(value, label)
    if number <= 0.0:
        raise DomainError(f"{label} must be positive.", field=label)
    return number


def _nonnegative(value: object, label: str) -> float:
    number = real(value, label)
    if number < 0.0:
        raise DomainError(f"{label} must not be negative.", field=label)
    return number


def _angle(value: object, label: str, upper: float = 90.0) -> float:
    degrees = real(value, label)
    if not 0.0 <= degrees <= upper:
        raise DomainError(f"{label} must be between 0 and {upper} degrees.", field=label)
    return radians(degrees)


def _cold_wall_factor(wall_enthalpy: float, reference_enthalpy: float, label: str) -> float:
    hw = _nonnegative(wall_enthalpy, "wall_enthalpy_j_kg")
    href = _positive(reference_enthalpy, label)
    if hw >= href:
        raise DomainError("Wall enthalpy must be below the driving enthalpy.", field="wall_enthalpy_j_kg")
    return 1.0 - hw / href


# ---------------------------------------------------------------------------
# Boundary-layer stagnation-point correlations
# ---------------------------------------------------------------------------


def newtonian_stagnation_velocity_gradient(
    nose_radius_m: float, edge_pressure_pa: float, freestream_pressure_pa: float, edge_density_kg_m3: float
) -> float:
    """Return ``(du_e/dx)_s = (1/R) sqrt(2 (pe - p_inf)/rho_e)``.

    Parameters
    ----------
    nose_radius_m : float
        Local radius of curvature at the stagnation point.
    edge_pressure_pa : float
        Stagnation (boundary-layer edge) pressure, above the free-stream value.
    freestream_pressure_pa : float
        Free-stream static pressure.
    edge_density_kg_m3 : float
        Density at the stagnation-point boundary-layer edge.

    Returns
    -------
    float
        Velocity gradient in 1/s.

    Notes
    -----
    Newtonian theory; Tauber reports it accurate for circular cylinders but
    7.7 percent low for spheres, and inaccurate for blunter faces (Fig. 4).

    References
    ----------
    .. [1] ANDERSON-FOA-2017. Eq. (18.82), printed p. 1042.
    .. [2] TAUBER-1989-TP2914. Eq. (39a) and p. 7, printed pp. 6-7.
    """

    radius = _positive(nose_radius_m, "nose_radius_m")
    pe = _positive(edge_pressure_pa, "edge_pressure_pa")
    pinf = _nonnegative(freestream_pressure_pa, "freestream_pressure_pa")
    rho = _positive(edge_density_kg_m3, "edge_density_kg_m3")
    if pe <= pinf:
        raise DomainError("edge_pressure_pa must exceed freestream_pressure_pa.", field="edge_pressure_pa")
    return sqrt(2.0 * (pe - pinf) / rho) / radius


def stagnation_point_heating_boundary_layer(
    *,
    geometry: str,
    prandtl: float,
    edge_density_kg_m3: float,
    edge_viscosity_pa_s: float,
    velocity_gradient_1_s: float,
    adiabatic_wall_enthalpy_j_kg: float,
    wall_enthalpy_j_kg: float,
) -> float:
    """Return laminar stagnation heating ``qw`` from the self-similar solution fit.

    Parameters
    ----------
    geometry : {"cylinder", "sphere"}
        Two-dimensional stagnation line or axisymmetric stagnation point.
    prandtl : float
        Prandtl number.
    edge_density_kg_m3, edge_viscosity_pa_s : float
        ``rho_e`` and ``mu_e`` at the stagnation-point edge.
    velocity_gradient_1_s : float
        ``du_e/dx`` (see :func:`newtonian_stagnation_velocity_gradient`).
    adiabatic_wall_enthalpy_j_kg, wall_enthalpy_j_kg : float
        ``h_aw`` and ``h_w``.

    Returns
    -------
    float
        Heat flux in W/m^2, positive into the wall.

    Notes
    -----
    ``qw = C Pr^-0.6 (rho_e mu_e)^0.5 (du_e/dx)^0.5 (haw - hw)`` with
    ``C = 0.57`` (cylinder, Eq. (18.65)) or ``C = 0.763`` (sphere,
    Eq. (18.70)) [1]_. Eq. (18.70) as printed shows ``Pr^-0.65`` and no
    exponent on ``rho_e mu_e``; the text on the same page states that the two
    expressions differ only in the leading coefficient, and the 1/2 power is
    required for dimensional consistency, so the cylinder exponents are used
    for both. Anderson attributes both fits to Van Driest (1956).

    References
    ----------
    .. [1] ANDERSON-FOA-2017. Eqs. (18.65), (18.70), printed p. 1039.
    """

    coefficients = {"cylinder": 0.57, "sphere": 0.763}
    if geometry not in coefficients:
        raise InputError("geometry must be 'cylinder' or 'sphere'.", field="geometry")
    pr = _positive(prandtl, "prandtl")
    rho = _positive(edge_density_kg_m3, "edge_density_kg_m3")
    mu = _positive(edge_viscosity_pa_s, "edge_viscosity_pa_s")
    gradient = _positive(velocity_gradient_1_s, "velocity_gradient_1_s")
    haw = _positive(adiabatic_wall_enthalpy_j_kg, "adiabatic_wall_enthalpy_j_kg")
    hw = _nonnegative(wall_enthalpy_j_kg, "wall_enthalpy_j_kg")
    return coefficients[geometry] * pr**-0.6 * sqrt(rho * mu * gradient) * (haw - hw)


# ---------------------------------------------------------------------------
# Sutton-Graves
# ---------------------------------------------------------------------------


@dataclass(frozen=True, slots=True)
class SuttonGravesBaseGas:
    """Base-gas transport data of Sutton-Graves Table I.

    Attributes
    ----------
    molecular_weight : float
        ``M`` in kg/kmol.
    transport_parameter : float
        ``gamma`` of their Eq. (37).
    m_gamma : float
        Tabulated product ``M gamma`` used in Eq. (42).

    References
    ----------
    .. [1] SUTTON-GRAVES-1971-TRR376. Table I, printed p. 38.
    """

    molecular_weight: float
    transport_parameter: float
    m_gamma: float


SUTTON_GRAVES_BASE_GASES: Mapping[str, SuttonGravesBaseGas] = MappingProxyType({
    "N2": SuttonGravesBaseGas(28.014, 0.03654, 1.024),
    "O2": SuttonGravesBaseGas(32.000, 0.04129, 1.321),
    "H2": SuttonGravesBaseGas(2.016, 0.06775, 0.137),
    "He": SuttonGravesBaseGas(4.003, 0.10845, 0.434),
    "Ne": SuttonGravesBaseGas(20.183, 0.07449, 1.504),
    "Ar": SuttonGravesBaseGas(39.948, 0.04036, 1.612),
    "CO2": SuttonGravesBaseGas(44.011, 0.02919, 1.285),
    "NH3": SuttonGravesBaseGas(17.031, 0.04605, 0.784),
    "CH4": SuttonGravesBaseGas(16.043, 0.03345, 0.537),
})

SUTTON_GRAVES_TABLE_II: Mapping[str, float] = MappingProxyType({
    "N2": 0.1112,
    "O2": 0.1201,
    "H2": 0.0395,
    "He": 0.0797,
    "Ne": 0.1474,
    "Ar": 0.1495,
    "CO2": 0.1210,
    "NH3": 0.0990,
    "CH4": 0.0807,
    "air (0.2320 O2-0.7680 N2)": 0.1113,
    "0.1345 CO2-0.8655 N2": 0.1113,
    "0.5000 CO2-0.5000 N2": 0.1131,
    "0.8500 CO2-0.1500 N2": 0.1173,
    "0.6714 CO2-0.3286 Ar": 0.1262,
    "0.3626 CO2-0.3077 N2-0.3297 Ar": 0.1218,
    "0.1339 CO2-0.8525 N2-0.0136 Ar": 0.1118,
})
"""Least-squares heat-transfer coefficients ``K`` of TR R-376 Table II (p. 39).

Units kg/(s m^1.5 atm^0.5). Table II also lists the average and maximum
correlation errors; for air they are 3.3 and 9.8 percent.
"""


def sutton_graves_coefficient(mass_fractions: Mapping[str, float], wall_prandtl: float | None = None) -> float:
    """Return the Sutton-Graves ``K`` of a gas mixture from its composition.

    Parameters
    ----------
    mass_fractions : mapping of str to float
        Mass fractions of the base gases in :data:`SUTTON_GRAVES_BASE_GASES`;
        must sum to one within 1e-6.
    wall_prandtl : float, optional
        ``None`` uses Eq. (42) (``Pr_w = 0.69``, coefficient 0.1106);
        otherwise Eq. (41) with coefficient ``0.0885 / Pr_w^0.6``.

    Returns
    -------
    float
        ``K`` in kg/(s m^1.5 atm^0.5).

    Notes
    -----
    ``K = c (sum_i c_i / (M_i gamma_i))^-1/2`` with the tabulated ``M gamma``.
    For air this gives 0.1150, 3.3 percent above the Table II correlation of
    the authors' boundary-layer computations (0.1113).

    References
    ----------
    .. [1] SUTTON-GRAVES-1971-TRR376. Eqs. (41)-(42), Table I, printed
           pp. 15, 38.
    """

    if not isinstance(mass_fractions, Mapping) or not mass_fractions:
        raise InputError("mass_fractions must be a non-empty mapping.", field="mass_fractions")
    total = 0.0
    weighted = 0.0
    for name, fraction in mass_fractions.items():
        if name not in SUTTON_GRAVES_BASE_GASES:
            raise InputError(f"Unknown Sutton-Graves base gas {name!r}.", field="mass_fractions")
        value = _nonnegative(fraction, f"mass_fractions[{name}]")
        total += value
        weighted += value / SUTTON_GRAVES_BASE_GASES[name].m_gamma
    if not isclose(total, 1.0, abs_tol=1.0e-6):
        raise DomainError("mass_fractions must sum to 1.", field="mass_fractions")
    if wall_prandtl is None:
        return 0.1106 * weighted**-0.5
    return 0.0885 / _positive(wall_prandtl, "wall_prandtl") ** 0.6 * weighted**-0.5


def sutton_graves_heating(
    *,
    stagnation_pressure_pa: float,
    nose_radius_m: float,
    stagnation_enthalpy_j_kg: float,
    wall_enthalpy_j_kg: float,
    coefficient: float = SUTTON_GRAVES_TABLE_II["air (0.2320 O2-0.7680 N2)"],
) -> float:
    """Return Sutton-Graves stagnation-point convective heating in W/m^2.

    Parameters
    ----------
    stagnation_pressure_pa : float
        Stagnation pressure at the boundary-layer edge, Pa (converted to atm
        with 1 atm = 101.325 kN/m^2 as stated in the report).
    nose_radius_m : float
        Equivalent body (nose) radius ``R``.
    stagnation_enthalpy_j_kg, wall_enthalpy_j_kg : float
        ``h_s`` and ``h_w`` (``h_s > h_w``).
    coefficient : float, optional
        ``K`` in kg/(s m^1.5 atm^0.5); default is the Table II value for air.

    Returns
    -------
    float
        ``qw = K sqrt(ps/R) (hs - hw)``, Eq. (33).

    Notes
    -----
    The correlation database spans 2.3-116.2 MJ/kg, 0.001-100 atm and wall
    temperatures of 300 and 1111 K, chemical equilibrium, axisymmetric blunt
    body; radiative heating was neglected [1]_.

    References
    ----------
    .. [1] SUTTON-GRAVES-1971-TRR376. Summary, Eq. (33), Table II, printed
           pp. 1, 13, 39.
    """

    ps_atm = _positive(stagnation_pressure_pa, "stagnation_pressure_pa") / _ATMOSPHERE_PA
    radius = _positive(nose_radius_m, "nose_radius_m")
    hs = _positive(stagnation_enthalpy_j_kg, "stagnation_enthalpy_j_kg")
    hw = _nonnegative(wall_enthalpy_j_kg, "wall_enthalpy_j_kg")
    k = _positive(coefficient, "coefficient")
    if hw >= hs:
        raise DomainError("wall_enthalpy_j_kg must be below the stagnation enthalpy.", field="wall_enthalpy_j_kg")
    return k * sqrt(ps_atm / radius) * (hs - hw)


# ---------------------------------------------------------------------------
# Tauber TP-2914 closed forms (air, SI units as printed)
# ---------------------------------------------------------------------------


def tauber_stagnation_heating(
    *,
    freestream_density_kg_m3: float,
    velocity_m_s: float,
    nose_radius_m: float,
    wall_enthalpy_j_kg: float,
    total_enthalpy_j_kg: float,
) -> float:
    """Return hemisphere stagnation heating ``1.83e-4 sqrt(rho/rn) V^3 (1 - hw/Hs)``.

    Parameters
    ----------
    freestream_density_kg_m3, velocity_m_s, nose_radius_m : float
        ``rho_inf``, ``V_inf`` and ``r_n``.
    wall_enthalpy_j_kg, total_enthalpy_j_kg : float
        ``h_w`` and ``H_s`` (Anderson's Example 14.2 uses ``H_s ~ V^2/2``).

    Returns
    -------
    float
        Heat flux in W/m^2.

    Notes
    -----
    Tauber Eq. (40), derived for air with ``h_inf << V^2/2`` [2]_; identical
    to Anderson Eq. (14.62) written in W/cm^2 [1]_. Tauber notes the constant
    rises to 1.90 if the 7.7 percent Newtonian velocity-gradient error for a
    sphere is removed; the printed 1.83 is used.

    References
    ----------
    .. [1] ANDERSON-FOA-2017. Eq. (14.62), Example 14.2, printed pp. 906-908.
    .. [2] TAUBER-1989-TP2914. Eq. (40), printed pp. 6-7.
    """

    rho = _positive(freestream_density_kg_m3, "freestream_density_kg_m3")
    v = _positive(velocity_m_s, "velocity_m_s")
    radius = _positive(nose_radius_m, "nose_radius_m")
    factor = _cold_wall_factor(wall_enthalpy_j_kg, total_enthalpy_j_kg, "total_enthalpy_j_kg")
    return 1.83e-4 * sqrt(rho / radius) * v**3 * factor


def tauber_swept_cylinder_heating(
    *,
    freestream_density_kg_m3: float,
    velocity_m_s: float,
    cylinder_radius_m: float,
    sweep_deg: float,
    wall_enthalpy_j_kg: float,
    freestream_enthalpy_j_kg: float,
) -> float:
    """Return stagnation-line heating of a swept infinite cylinder, Tauber Eq. (41).

    Parameters
    ----------
    freestream_density_kg_m3, velocity_m_s, cylinder_radius_m : float
        ``rho_inf``, ``V_inf`` and ``r_cyl``.
    sweep_deg : float
        Sweep angle ``Lambda`` in degrees, ``[0, 90)``.
    wall_enthalpy_j_kg, freestream_enthalpy_j_kg : float
        ``h_w`` and ``h_inf``; ``h_aw = h_inf + 0.5 V^2 (1 - 0.18 sin^2 Lambda)``.

    Returns
    -------
    float
        ``1.29e-4 sqrt(rho/r) (1 - 0.18 sin^2 L) V^3 (1 - hw/haw) cos L`` in W/m^2.

    References
    ----------
    .. [1] TAUBER-1989-TP2914. Eq. (41), printed p. 7.
    """

    rho = _positive(freestream_density_kg_m3, "freestream_density_kg_m3")
    v = _positive(velocity_m_s, "velocity_m_s")
    radius = _positive(cylinder_radius_m, "cylinder_radius_m")
    sweep = _angle(sweep_deg, "sweep_deg")
    if sweep >= radians(90.0):
        raise DomainError("sweep_deg must be below 90 degrees.", field="sweep_deg")
    shape = 1.0 - 0.18 * sin(sweep) ** 2
    haw = _nonnegative(freestream_enthalpy_j_kg, "freestream_enthalpy_j_kg") + 0.5 * v * v * shape
    factor = _cold_wall_factor(wall_enthalpy_j_kg, haw, "freestream_enthalpy_j_kg")
    return 1.29e-4 * sqrt(rho / radius) * shape * v**3 * factor * cos(sweep)


def _laminar_recovery_enthalpy(freestream_enthalpy: float, velocity: float) -> float:
    return _nonnegative(freestream_enthalpy, "freestream_enthalpy_j_kg") + 0.40 * velocity * velocity


def tauber_laminar_cone_heating(
    *,
    freestream_density_kg_m3: float,
    velocity_m_s: float,
    distance_m: float,
    cone_half_angle_deg: float,
    wall_enthalpy_j_kg: float,
    freestream_enthalpy_j_kg: float,
) -> float:
    """Return laminar heating on a sharp cone in air, Tauber Eq. (46).

    Parameters
    ----------
    freestream_density_kg_m3, velocity_m_s : float
        ``rho_inf`` and ``V_inf``.
    distance_m : float
        Surface distance ``x`` from the apex.
    cone_half_angle_deg : float
        ``delta_c`` in degrees.
    wall_enthalpy_j_kg, freestream_enthalpy_j_kg : float
        ``h_w`` and ``h_inf``; ``h_aw = h_inf + 0.40 V^2`` as printed.

    Returns
    -------
    float
        ``4.03e-5 (rho cos d / x)^0.5 V^3.2 sin d (1 - hw/haw)`` in W/m^2.

    Notes
    -----
    Newtonian edge flow and high-speed approximation (``h_inf << V^2/2``).
    Tauber notes nose blunting usually lowers the heating aft of the nose.

    References
    ----------
    .. [1] TAUBER-1989-TP2914. Eqs. (45)-(46), printed p. 9.
    """

    rho = _positive(freestream_density_kg_m3, "freestream_density_kg_m3")
    v = _positive(velocity_m_s, "velocity_m_s")
    x = _positive(distance_m, "distance_m")
    delta = _angle(cone_half_angle_deg, "cone_half_angle_deg")
    factor = _cold_wall_factor(wall_enthalpy_j_kg, _laminar_recovery_enthalpy(freestream_enthalpy_j_kg, v),
                               "freestream_enthalpy_j_kg")
    return 4.03e-5 * sqrt(rho * cos(delta) / x) * v**3.2 * sin(delta) * factor


def tauber_laminar_flat_plate_heating(
    *,
    freestream_density_kg_m3: float,
    velocity_m_s: float,
    distance_m: float,
    plate_angle_deg: float,
    wall_enthalpy_j_kg: float,
    freestream_enthalpy_j_kg: float,
) -> float:
    """Return laminar heating on an inclined flat plate or wedge, Tauber Eq. (48).

    Parameters
    ----------
    freestream_density_kg_m3, velocity_m_s : float
        ``rho_inf`` and ``V_inf``.
    distance_m : float
        Distance ``x`` from the leading edge.
    plate_angle_deg : float
        ``delta_FP``: wedge angle or angle of attack, degrees.
    wall_enthalpy_j_kg, freestream_enthalpy_j_kg : float
        ``h_w`` and ``h_inf``; ``h_aw = h_inf + 0.40 V^2`` (stated with Eq. (46),
        to which Eq. (48) is equivalent).

    Returns
    -------
    float
        ``2.42e-5 (rho cos d / x)^0.5 V^3.2 sin d (1 - hw/haw)`` in W/m^2.

    References
    ----------
    .. [1] TAUBER-1989-TP2914. Eqs. (47)-(48), printed p. 9.
    """

    rho = _positive(freestream_density_kg_m3, "freestream_density_kg_m3")
    v = _positive(velocity_m_s, "velocity_m_s")
    x = _positive(distance_m, "distance_m")
    delta = _angle(plate_angle_deg, "plate_angle_deg")
    factor = _cold_wall_factor(wall_enthalpy_j_kg, _laminar_recovery_enthalpy(freestream_enthalpy_j_kg, v),
                               "freestream_enthalpy_j_kg")
    return 2.42e-5 * sqrt(rho * cos(delta) / x) * v**3.2 * sin(delta) * factor


def tauber_leading_edge_heating(cylinder_heating_w_m2: float, flat_plate_heating_w_m2: float,
                                sweep_deg: float, angle_of_attack_deg: float) -> float:
    """Combine swept-cylinder and flat-plate heating on a wing leading edge.

    Parameters
    ----------
    cylinder_heating_w_m2 : float
        ``q_cyl`` from :func:`tauber_swept_cylinder_heating`.
    flat_plate_heating_w_m2 : float
        ``q_FP`` (laminar Eq. (48) or turbulent Eq. (56)).
    sweep_deg, angle_of_attack_deg : float
        ``Lambda`` and ``alpha`` in degrees.

    Returns
    -------
    float
        ``(q_cyl^2 + q_FP^2 sin^2 Lambda)^0.5 cos alpha`` (Rubesin, as quoted).

    Notes
    -----
    Tauber limits this to small angles of attack for highly swept edges.

    References
    ----------
    .. [1] TAUBER-1989-TP2914. Eq. (49), printed p. 9.
    """

    q_cyl = _nonnegative(cylinder_heating_w_m2, "cylinder_heating_w_m2")
    q_fp = _nonnegative(flat_plate_heating_w_m2, "flat_plate_heating_w_m2")
    sweep = _angle(sweep_deg, "sweep_deg")
    alpha = _angle(angle_of_attack_deg, "angle_of_attack_deg")
    return sqrt(q_cyl * q_cyl + (q_fp * sin(sweep)) ** 2) * cos(alpha)


def tauber_turbulent_flat_plate_heating(
    *,
    freestream_density_kg_m3: float,
    velocity_m_s: float,
    distance_from_transition_m: float,
    plate_angle_deg: float,
    wall_temperature_k: float,
    wall_enthalpy_j_kg: float,
    total_enthalpy_j_kg: float,
) -> float:
    """Return turbulent heating on an inclined flat plate, Tauber Eqs. (56a)/(56b).

    Parameters
    ----------
    freestream_density_kg_m3, velocity_m_s : float
        ``rho_inf`` and ``V_inf``; Eq. (56a) applies for
        ``1500 < V <= 3960`` m/s and Eq. (56b) above 3960 m/s.
    distance_from_transition_m : float
        ``x - x_bt``, distance from the beginning of transition (Tauber's
        effective turbulent origin), positive.
    plate_angle_deg : float
        ``delta``: plate angle to the free stream, degrees.
    wall_temperature_k : float
        ``T_w`` (used by Eq. (56a) only).
    wall_enthalpy_j_kg, total_enthalpy_j_kg : float
        ``h_w`` and the total enthalpy ``H`` of the inviscid flow.

    Returns
    -------
    float
        Heat flux in W/m^2.

    Notes
    -----
    Correlation of reference-enthalpy calculations (Eckert ``h'``, Blasius
    turbulent law below ``Re = 1e7`` and Schultz-Grunow above), stated valid
    to +/-15 percent for altitudes 3-60 km and speeds of about 1-11 km/s.

    References
    ----------
    .. [1] TAUBER-1989-TP2914. Eqs. (53)-(58), printed pp. 13-14.
    """

    rho = _positive(freestream_density_kg_m3, "freestream_density_kg_m3")
    v = _positive(velocity_m_s, "velocity_m_s")
    x = _positive(distance_from_transition_m, "distance_from_transition_m")
    delta = _angle(plate_angle_deg, "plate_angle_deg")
    tw = _positive(wall_temperature_k, "wall_temperature_k")
    hw = _nonnegative(wall_enthalpy_j_kg, "wall_enthalpy_j_kg")
    total = _positive(total_enthalpy_j_kg, "total_enthalpy_j_kg")
    driving = 0.9 - hw / total
    if driving <= 0.0:
        raise DomainError("0.9 - hw/H must be positive.", field="wall_enthalpy_j_kg")
    if v <= 1500.0:
        raise DomainError("Tauber Eq. (56) is stated for V > 1500 m/s.", field="velocity_m_s")
    s2 = sin(delta) ** 2
    if v <= 3960.0:
        return 3.72e-4 * (rho * s2 * cos(delta) ** 2.22) ** 0.8 / (x**0.2 * (tw / 555.0) ** 0.25) * v**3.37 * driving
    return 2.45e-5 * (rho * s2 * cos(delta) ** 2.62) ** 0.8 / x**0.2 * v**3.7 * driving


def tauber_turbulent_cone_heating(flat_plate_heating_w_m2: float) -> float:
    """Return approximate turbulent sharp-cone heating ``1.15 q_FP``, Tauber Eq. (59).

    Parameters
    ----------
    flat_plate_heating_w_m2 : float
        Turbulent flat-plate value evaluated with the cone half angle.

    Returns
    -------
    float
        Heat flux in W/m^2.

    References
    ----------
    .. [1] TAUBER-1989-TP2914. Eq. (59), printed p. 14.
    """

    return 1.15 * _nonnegative(flat_plate_heating_w_m2, "flat_plate_heating_w_m2")
