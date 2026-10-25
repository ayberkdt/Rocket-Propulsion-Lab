"""Propellant vaporization in the chamber: Priem-Heidmann generalized length.

Priem and Heidmann computed how much of an injected spray vaporizes in a
given chamber length from a single-drop model (Ranz-Marshall heat and mass
transfer, drag ``C_D = 27 Re^-0.84``, gas velocity proportional to the mass
already vaporized) applied to a log-normal spray, and correlated the result
for several propellants against one "generalized length" [1]_ [2]_. This
module provides

* the generalized length ``L_gen`` in the form NASA CR-168066 implements from
  Priem and Heidmann's NASA TR R-67 (1960) [3]_,
* the percent of mass vaporized as a function of ``L_gen``, the 35-point
  average curve tabulated in CR-168066 and interpolated linearly as its
  program does [3]_,
* the effective length of NACA TN 3883, the propellant-independent
  correlation that preceded R-67 [2]_,
* the c* efficiency of a partly vaporized propellant pair, through Heidmann
  and Priem's NACA TN 4219 Eq. (5) [4]_.

``L_gen`` (inches) =
``L (p_c/300)^0.66 / CR^0.44 / [(1 - T/T_c)^0.4 (H_v/140)^0.8 (M/100)^0.35
(r/0.003)^1.45 (V/100)^0.75]`` with ``L`` in inches, ``p_c`` in psia,
``H_v`` in Btu/lb, ``r`` (mass-median drop radius) in inches and ``V``
(injection velocity) in ft/s; this module takes and returns SI and converts.

Limits: TR R-67 itself is not on NTRS and was not read; the generalized
length and the curve are quoted from CR-168066, a secondary source, whose
printed SI form (eq. 1, p. 38) is dimensionally inconsistent with its own
program (``V/328`` with ``V`` in m/s, ``H_v`` "in J/kg" over 325.64), so the
program's English-unit form is used. Only the cylindrical chamber length is
included (the program omits the convergent-section term of eq. 1). The curve
is one average curve; the spray spread it belongs to is not stated
(Priem's sprays use ``sigma_g = 2.3``, TN 3883 p. 5). Combustion is assumed
instantaneous once vaporized.

References
----------
.. [1] NACA-TN-3985-PRIEM-1957. R. J. Priem, NACA TN 3985, 1957,
       Eqs. (1)-(12), pp. 2-4, https://ntrs.nasa.gov/citations/19930084710
.. [2] NACA-TN-3883-PRIEM-1958. R. J. Priem, NACA TN 3883, 1958, effective
       length, p. 5, and Table I, https://ntrs.nasa.gov/citations/19930084722
.. [3] NASA-CR-168066-AEROJET-1983. Aerojet, NASA CR-168066, 1983, eq. (1),
       p. 38, and Appendix C program listing (``DATA LG/EV``, statement
       ``CVO``, subroutine ``SINTP``), pp. C-6 to C-8,
       https://ntrs.nasa.gov/citations/19830023615
.. [4] NACA-TN4219. Heidmann and Priem, NACA TN 4219 (1958), equations
       (1)-(7), printed pp. 3-4. https://ntrs.nasa.gov/citations/19930085136
"""

from __future__ import annotations

from bisect import bisect_right
from collections.abc import Callable
from dataclasses import dataclass

from rocket_propulsion.core.errors import DomainError
from rocket_propulsion.core.validation import numeric_range_guard, real

from .combustion_efficiency import vaporization_limited_cstar_efficiency

__all__ = [
    "PRIEM_HEIDMANN_CURVE",
    "VaporizationEfficiency",
    "generalized_vaporization_length",
    "percent_vaporized",
    "priem_effective_length",
    "vaporization_cstar_efficiency",
]

_INCH_M = 0.0254
_FOOT_M = 0.3048
_PSI_PA = 6894.757293168361
#: International-table Btu per pound in J/kg (1055.05585262 J / 0.45359237 kg).
_BTU_PER_LB_J_KG = 2326.0

#: CR-168066 Appendix C ``DATA LG/EV`` (p. C-6): generalized length in inches
#: and percent of mass vaporized, the average curve of NASA TR R-67.
PRIEM_HEIDMANN_CURVE: tuple[tuple[float, float], ...] = (
    (0.01, 0.083), (0.02, 0.31), (0.04, 0.91), (0.06, 1.61), (0.1, 3.10),
    (0.2, 6.5), (0.4, 12.2), (0.6, 17.3), (0.8, 21.8), (1.0, 25.9),
    (1.5, 34.2), (2.0, 41.5), (3.0, 52.2), (4.0, 60.0), (5.0, 66.0),
    (6.0, 70.5), (7.0, 74.5), (8.0, 77.8), (9.0, 80.5), (10.0, 82.7),
    (13.0, 87.8), (15.0, 90.2), (17.0, 92.0), (20.0, 94.0), (25.0, 96.2),
    (30.0, 97.6), (35.0, 98.5), (40.0, 99.1), (45.0, 99.4), (50.0, 99.6),
    (60.0, 99.82), (70.0, 99.92), (80.0, 99.96), (90.0, 99.98), (100.0, 99.99),
)
_LG = tuple(p[0] for p in PRIEM_HEIDMANN_CURVE)
_EV = tuple(p[1] for p in PRIEM_HEIDMANN_CURVE)


def _positive(value: object, label: str) -> float:
    number = real(value, label)
    if number <= 0.0:
        raise DomainError(f"{label} must be positive.", field=label)
    return number


def _reduced_temperature_term(liquid_temperature_k: float, critical_temperature_k: float) -> float:
    t = _positive(liquid_temperature_k, "liquid_temperature_k")
    tc = _positive(critical_temperature_k, "critical_temperature_k")
    if t >= tc:
        raise DomainError("The model does not apply at or above the critical temperature "
                          "(TN 3883 p. 5: 1 - T_r must be positive).",
                          field="liquid_temperature_k")
    return 1.0 - t / tc


@numeric_range_guard("chamber_pressure_pa")
def generalized_vaporization_length(*, chamber_length_m: float, chamber_pressure_pa: float,
                                    contraction_ratio: float, liquid_temperature_k: float,
                                    critical_temperature_k: float,
                                    heat_of_vaporization_j_kg: float, molecular_weight: float,
                                    drop_radius_m: float, injection_velocity_m_s: float) -> float:
    """Priem-Heidmann generalized length ``L_gen`` of one propellant.

    Parameters
    ----------
    chamber_length_m : float
        Cylindrical chamber length from the injector face, ``L``.
    chamber_pressure_pa : float
        Chamber pressure.
    contraction_ratio : float
        Chamber-to-throat area ratio ``CR >= 1``.
    liquid_temperature_k : float
        Injected liquid temperature ``T``.
    critical_temperature_k : float
        Critical temperature ``T_c`` of the propellant.
    heat_of_vaporization_j_kg : float
        Heat of vaporization ``H_v``.
    molecular_weight : float
        Molecular weight ``M`` in g/mol.
    drop_radius_m : float
        Mass-median drop radius ``r``.
    injection_velocity_m_s : float
        Drop injection velocity ``V``.

    Returns
    -------
    float
        ``L_gen`` in metres (the source's inches times 0.0254).

    References
    ----------
    .. [1] NASA-CR-168066-AEROJET-1983. Aerojet, NASA CR-168066, 1983, eq.
           (1), p. 38, and program statement ``CVO`` with ``LG01 = CVO *
           DLCO``, pp. C-7 to C-8, https://ntrs.nasa.gov/citations/19830023615
    """
    length_in = _positive(chamber_length_m, "chamber_length_m") / _INCH_M
    pc_psia = _positive(chamber_pressure_pa, "chamber_pressure_pa") / _PSI_PA
    cr = _positive(contraction_ratio, "contraction_ratio")
    if cr < 1.0:
        raise DomainError("contraction_ratio must be at least 1.", field="contraction_ratio")
    reduced = _reduced_temperature_term(liquid_temperature_k, critical_temperature_k)
    hv = _positive(heat_of_vaporization_j_kg, "heat_of_vaporization_j_kg") / _BTU_PER_LB_J_KG
    mw = _positive(molecular_weight, "molecular_weight")
    r_in = _positive(drop_radius_m, "drop_radius_m") / _INCH_M
    v_fps = _positive(injection_velocity_m_s, "injection_velocity_m_s") / _FOOT_M
    factor = ((pc_psia / 300.0) ** 0.66 / cr ** 0.44) / (
        reduced ** 0.4 * (hv / 140.0) ** 0.8 * (mw / 100.0) ** 0.35
        * (r_in / 0.003) ** 1.45 * (v_fps / 100.0) ** 0.75)
    return length_in * factor * _INCH_M


def percent_vaporized(generalized_length_m: float) -> float:
    """Percent of mass vaporized at a generalized length (TR R-67 average curve).

    Linear interpolation in the 35-point table of CR-168066, as its
    subroutine ``SINTP`` does; outside 0.01-100 in. the program stops with
    "TABLE TOO SMALL" and this function raises.

    Parameters
    ----------
    generalized_length_m : float
        ``L_gen`` from :func:`generalized_vaporization_length`.

    Returns
    -------
    float
        Percent vaporized, 0.083 to 99.99.

    References
    ----------
    .. [1] NASA-CR-168066-AEROJET-1983. Aerojet, NASA CR-168066, 1983,
           ``DATA LG/EV`` and ``SUBROUTINE SINTP``, pp. C-6 and C-8, and
           Fig. 26, p. 39, https://ntrs.nasa.gov/citations/19830023615
    """
    x = _positive(generalized_length_m, "generalized_length_m") / _INCH_M
    if x < _LG[0] or x > _LG[-1]:
        raise DomainError(
            f"L_gen = {x:.4g} in lies outside the tabulated 0.01-100 in.",
            code="outside_correlation_range", field="generalized_length_m",
            details={"value_in": x, "minimum_in": _LG[0], "maximum_in": _LG[-1]})
    i = min(max(bisect_right(_LG, x), 1), len(_LG) - 1)
    x0, x1, y0, y1 = _LG[i - 1], _LG[i], _EV[i - 1], _EV[i]
    return y0 + (y1 - y0) * (x - x0) / (x1 - x0)


def priem_effective_length(*, chamber_length_m: float, chamber_pressure_pa: float,
                           final_gas_velocity_m_s: float, injection_velocity_m_s: float,
                           mass_median_radius_m: float, reduced_temperature: float) -> float:
    """Effective length of NACA TN 3883, in the source's mixed English units.

    ``L_eff = L P^0.66 u_fin^0.40 (1.9e-5) / [(1 - T_r)^0.40 M^1.45 v_0^0.75]``
    with ``L`` and ``M`` (mass-median radius) in inches, ``P`` in psia and
    velocities in in./s (units inferred from the axis dimension of TN 3883
    Figs. 2-3, (lb)^0.66 (sec)^0.35 / (in.)^2.12, which the text does not
    spell out); the value is on the abscissa of those figures,
    which hold for all five propellants of that report. It shares the
    pressure, radius, velocity and ``(1 - T_r)`` exponents with
    :func:`generalized_vaporization_length`; the final-gas-velocity term is
    replaced there by the property and contraction-ratio factors of R-67.

    Parameters
    ----------
    chamber_length_m, chamber_pressure_pa : float
        ``L`` and ``P``.
    final_gas_velocity_m_s : float
        Gas velocity when all propellant has burned, ``u_fin``.
    injection_velocity_m_s : float
        Drop injection velocity ``v_0``.
    mass_median_radius_m : float
        ``M_g,M``.
    reduced_temperature : float
        ``T_r = T_l,0 / T_c`` in (0, 1).

    Returns
    -------
    float
        Effective length on the TN 3883 scale.

    References
    ----------
    .. [1] NACA-TN-3883-PRIEM-1958. R. J. Priem, NACA TN 3883, 1958, p. 5,
           https://ntrs.nasa.gov/citations/19930084722
    """
    length_in = _positive(chamber_length_m, "chamber_length_m") / _INCH_M
    p_psia = _positive(chamber_pressure_pa, "chamber_pressure_pa") / _PSI_PA
    u_fin = _positive(final_gas_velocity_m_s, "final_gas_velocity_m_s") / _INCH_M
    v0 = _positive(injection_velocity_m_s, "injection_velocity_m_s") / _INCH_M
    radius_in = _positive(mass_median_radius_m, "mass_median_radius_m") / _INCH_M
    tr = _positive(reduced_temperature, "reduced_temperature")
    if tr >= 1.0:
        raise DomainError("reduced_temperature must be below 1.", field="reduced_temperature")
    return (length_in * p_psia ** 0.66 * u_fin ** 0.40 * 1.9e-5
            / ((1.0 - tr) ** 0.40 * radius_in ** 1.45 * v0 ** 0.75))


@dataclass(frozen=True, slots=True)
class VaporizationEfficiency:
    """Vaporized fractions and the resulting c* efficiency.

    Attributes
    ----------
    oxidizer_percent_vaporized, fuel_percent_vaporized : float
        From :func:`percent_vaporized`.
    cstar_efficiency : float
        NACA TN 4219 Eq. (5).

    References
    ----------
    .. [1] NACA-TN4219. Heidmann and Priem, NACA TN 4219 (1958), Eq. (5),
           https://ntrs.nasa.gov/citations/19930085136
    """

    oxidizer_percent_vaporized: float
    fuel_percent_vaporized: float
    cstar_efficiency: float


def vaporization_cstar_efficiency(*, characteristic_velocity: Callable[[float], float],
                                  oxidizer_fuel_ratio: float,
                                  oxidizer_generalized_length_m: float,
                                  fuel_generalized_length_m: float) -> VaporizationEfficiency:
    """c* efficiency of a pair whose vaporization is set by the generalized lengths.

    A propellant injected as a gas is represented by a generalized length at
    the top of the table (100 in., 99.99 %).

    Parameters
    ----------
    characteristic_velocity : callable
        Theoretical ``c*(O/F)`` in m/s.
    oxidizer_fuel_ratio : float
        Injected O/F.
    oxidizer_generalized_length_m, fuel_generalized_length_m : float
        ``L_gen`` of each propellant.

    Returns
    -------
    VaporizationEfficiency

    References
    ----------
    .. [1] NACA-TN4219. Heidmann and Priem, NACA TN 4219 (1958), Eq. (5),
           https://ntrs.nasa.gov/citations/19930085136
    .. [2] NASA-CR-168066-AEROJET-1983. Aerojet, NASA CR-168066, 1983,
           https://ntrs.nasa.gov/citations/19830023615
    """
    ox = percent_vaporized(oxidizer_generalized_length_m)
    fu = percent_vaporized(fuel_generalized_length_m)
    eta = vaporization_limited_cstar_efficiency(
        characteristic_velocity=characteristic_velocity, oxidizer_fuel_ratio=oxidizer_fuel_ratio,
        oxidizer_vaporized_fraction=ox / 100.0, fuel_vaporized_fraction=fu / 100.0)
    return VaporizationEfficiency(ox, fu, eta)
