"""Thermal radiation for rocket engines and vehicles.

Contents, each implemented as printed in the source opened for this module:

* Black-body physics: Stefan-Boltzmann law and Planck's law (Lienhard
  Eqs. (1.28), (1.30)) [1]_ with the SI defining constants as listed by
  NIST CODATA [2]_; black-body band fractions from a series derived here
  from NIST DLMF identities [3]_ (see :func:`blackbody_fraction`).
* View factors: Lienhard Tables 10.2 and 10.3 [1]_ and Howell's catalog
  factors C-11, C-14, C-40, C-41, C-78 and C-111 [4]_; summation and
  reciprocity (Lienhard Eqs. (10.12), (10.15)); a frustum (nozzle-extension)
  enclosure built from C-41 with that algebra.
* Gray diffuse surfaces: two-surface circuit (Eqs. (10.23)-(10.27)),
  radiation shields (Eq. (10.31)) and the N-surface net-radiation method
  (Eqs. (10.39)-(10.44)) [1]_.
* Gas radiation: mean beam length (Eq. (10.58), Table 10.4), the total
  absorptivity scaling rule Eq. (10.57) and the gas-to-black-wall exchange of
  Example 10.11 [1]_; weighted-sum-of-gray-gases (WSGG) coefficients for
  pure H2O and pure CO2 at 1 atm (Bordbar, Fraga and Hostikka 2020,
  Tables 2-3) [5]_ and for soot (Sadeghi et al. 2021, Table B9), combined
  with the superposition method (Sadeghi et al. Eqs. (12)-(14)) [6]_.
* Radiative-equilibrium wall temperature (energy balance with the
  transfer-factor form of Lienhard Eq. (10.2)) and the gray interior of a
  conical nozzle extension radiating through its two openings.

Not included because no implementable relation was opened and verified in
this work package: alumina (Al2O3) particle radiation, gray-wall gas
enclosures beyond the black-wall exchange of Example 10.11, pressure-
dependent WSGG models (all WSGG data here are 1 atm fits), and the CO WSGG
coefficients of Sadeghi et al. Table B8 (as printed, weight ``a_2`` becomes
negative above about 1300 K, so the table could not be confirmed).

References
----------
.. [1] LIENHARD-AHTT-6.00. J. H. Lienhard IV and J. H. Lienhard V, A Heat
       Transfer Textbook, 6th ed., version 6.00 (19 April 2024), Chapters 1
       and 10, https://ahtt.mit.edu/
.. [2] NIST-CODATA-2022. NIST, CODATA recommended values of the
       fundamental physical constants (2022): h, c, k, sigma, b and c2,
       https://physics.nist.gov/cuu/Constants/index.html
.. [3] DLMF-25.5. NIST DLMF Eqs. (25.5.1), (25.6.1), (8.2.2), (8.4.8),
       (24.2.1) and Table 24.2.1, https://dlmf.nist.gov/25.5
.. [4] HOWELL-CATALOG-C. J. R. Howell, A Catalog of Radiation Heat Transfer
       Configuration Factors, Section C,
       https://www.thermalradiation.net/indexCat.html
.. [5] BORDBAR-FRAGA-HOSTIKKA-2020. H. Bordbar, G. C. Fraga and S. Hostikka,
       Int. Commun. Heat Mass Transfer 110 (2020) 104400, accepted
       manuscript, Eqs. (2)-(3) and Tables 2-3,
       https://doi.org/10.1016/j.icheatmasstransfer.2019.104400
.. [6] SADEGHI-2021-FSJ. H. Sadeghi, S. Hostikka, G. C. Fraga and
       H. Bordbar, Fire Safety Journal 125 (2021) 103420, Eqs. (9)-(14),
       (24) and Table B9, https://doi.org/10.1016/j.firesaf.2021.103420
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from math import asin, atan, exp, expm1, factorial, inf, log, pi, sin, sqrt
from types import MappingProxyType

from rocket_propulsion.core.constants import (
    BOLTZMANN_CONSTANT_J_K,
    PLANCK_CONSTANT_J_S,
    SPEED_OF_LIGHT_M_S,
    STEFAN_BOLTZMANN_W_M2_K4,
    USSA1976_SEA_LEVEL_PRESSURE_PA,
)
from rocket_propulsion.core.errors import ConvergenceError, DomainError, InputError
from rocket_propulsion.core.metadata import WarningMessage
from rocket_propulsion.core.validation import real

__all__ = [
    "BOLTZMANN_CONSTANT_J_K",
    "BORDBAR_2020_CO2",
    "BORDBAR_2020_H2O",
    "FIRST_RADIATION_CONSTANT_W_M2",
    "MEAN_BEAM_LENGTH_FACTORS",
    "PLANCK_CONSTANT_J_S",
    "SADEGHI_2021_SOOT",
    "SECOND_RADIATION_CONSTANT_M_K",
    "SPEED_OF_LIGHT_M_S",
    "STEFAN_BOLTZMANN_W_M2_K4",
    "WIEN_DISPLACEMENT_CONSTANT_M_K",
    "EnclosureSolution",
    "FrustumEnclosure",
    "GasRadiation",
    "GraySurface",
    "RadiationShield",
    "RadiativeEquilibrium",
    "RadiativeSink",
    "ShieldedExchange",
    "WsggSpecies",
    "blackbody_band_fraction",
    "blackbody_emissive_power",
    "blackbody_fraction",
    "enclosed_body_transfer_factor",
    "frustum_enclosure_view_factors",
    "frustum_interior_transfer_factors",
    "gas_to_black_wall_heat_flux",
    "mean_beam_length",
    "mean_beam_length_for_geometry",
    "parallel_plate_transfer_factor",
    "planck_spectral_emissive_power",
    "radiative_equilibrium_wall_temperature",
    "reciprocal_view_factor",
    "shielded_heat_rate",
    "solve_gray_enclosure",
    "soot_planck_mean_absorption_coefficient",
    "station_gas_emissivity",
    "two_surface_heat_rate",
    "view_factor_coaxial_disks",
    "view_factor_concentric_cylinders",
    "view_factor_concentric_spheres",
    "view_factor_cylinder_interior_to_itself",
    "view_factor_frustum_to_base",
    "view_factor_hinged_strips",
    "view_factor_matrix_residuals",
    "view_factor_parallel_cylinders",
    "view_factor_parallel_rectangles",
    "view_factor_parallel_strips",
    "view_factor_perpendicular_rectangles",
    "view_factor_perpendicular_strips",
    "view_factor_strip_to_cylinder",
    "view_factor_triangular_duct",
    "wien_peak_wavelength",
    "wsgg_mixture_absorptivity",
    "wsgg_mixture_emissivity",
    "wsgg_species_emissivity",
    "wsgg_weights",
]

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

# The SI defining constants and sigma are registered in rocket_propulsion.core.constants
# (BIPM SI Brochure 2019); sigma = 2 pi^5 k^4 / (15 h^3 c^2) agrees with NIST CODATA 2022.


FIRST_RADIATION_CONSTANT_W_M2 = 2.0 * pi * PLANCK_CONSTANT_J_S * SPEED_OF_LIGHT_M_S**2
"""``2 pi h c^2``, the numerator of Planck's law, Lienhard Eq. (1.30)."""

SECOND_RADIATION_CONSTANT_M_K = PLANCK_CONSTANT_J_S * SPEED_OF_LIGHT_M_S / BOLTZMANN_CONSTANT_J_K
"""``h c / k_B``, the exponent constant of Lienhard Eq. (1.30) (NIST ``c2``)."""


def _wien_root() -> float:
    # d/d(lambda) of Eq. (1.30) vanishes where x = 5 (1 - exp(-x)), x = c2/(lambda T).
    x = 5.0
    for _ in range(50):
        step = (x - 5.0 * (1.0 - exp(-x))) / (1.0 - 5.0 * exp(-x))
        x -= step
        if abs(step) < 1e-16 * x:
            break
    return x


WIEN_DISPLACEMENT_CONSTANT_M_K = SECOND_RADIATION_CONSTANT_M_K / _wien_root()
"""``(lambda T)`` at the maximum of Eq. (1.30); Lienhard Eq. (1.29) prints 2897.77 um K."""

_ATM_PA = USSA1976_SEA_LEVEL_PRESSURE_PA
"""1 atm = 101325 Pa, the total pressure of the WSGG fits."""

_TOL = 1e-12


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


def _emissivity(value: object, label: str) -> float:
    number = real(value, label)
    if not 0.0 < number <= 1.0:
        raise DomainError(f"{label} must be in (0, 1].", field=label)
    return number


def _fraction(value: object, label: str) -> float:
    number = real(value, label)
    if not 0.0 <= number <= 1.0:
        raise DomainError(f"{label} must be in [0, 1].", field=label)
    return number


# ---------------------------------------------------------------------------
# Black body
# ---------------------------------------------------------------------------


def blackbody_emissive_power(temperature_k: float) -> float:
    """Return the black-body emissive power ``e_b = sigma T^4`` in W/m^2.

    Parameters
    ----------
    temperature_k : float
        Absolute temperature, nonnegative.

    Returns
    -------
    float
        Emissive power in W/m^2.

    References
    ----------
    .. [1] LIENHARD-AHTT-6.00. Stefan-Boltzmann law, Eq. (1.28), p. 30.
    """

    temperature = _nonnegative(temperature_k, "temperature_k")
    return STEFAN_BOLTZMANN_W_M2_K4 * temperature**4


def planck_spectral_emissive_power(wavelength_m: float, temperature_k: float) -> float:
    """Return the monochromatic black-body emissive power ``e_lambda,b`` in W/m^3.

    Parameters
    ----------
    wavelength_m : float
        Vacuum wavelength, positive.
    temperature_k : float
        Absolute temperature, positive.

    Returns
    -------
    float
        ``2 pi h c^2 / (lambda^5 [exp(h c / (k T lambda)) - 1])`` per metre
        of wavelength (multiply by 1e-6 for W/(m^2 um)).

    Notes
    -----
    Evaluated as ``c1 lambda^-5 exp(-x) / (1 - exp(-x))`` with ``expm1`` so
    that neither the Wien tail nor the Rayleigh-Jeans limit overflows.

    References
    ----------
    .. [1] LIENHARD-AHTT-6.00. Planck's law, Eq. (1.30), p. 30.
    """

    wavelength = _positive(wavelength_m, "wavelength_m")
    temperature = _positive(temperature_k, "temperature_k")
    x = SECOND_RADIATION_CONSTANT_M_K / (wavelength * temperature)
    return FIRST_RADIATION_CONSTANT_W_M2 / wavelength**5 * exp(-x) / -expm1(-x)


# Bernoulli numbers B_0..B_16 (DLMF Table 24.2.1); odd ones beyond B_1 vanish.
_BERNOULLI = {0: 1.0, 1: -0.5, 2: 1.0 / 6.0, 4: -1.0 / 30.0, 6: 1.0 / 42.0, 8: -1.0 / 30.0,
              10: 5.0 / 66.0, 12: -691.0 / 2730.0, 14: 7.0 / 6.0, 16: -3617.0 / 510.0}
_PLANCK_INTEGRAL = pi**4 / 15.0  # Gamma(4) zeta(4) = 6 pi^4/90, DLMF (25.5.1), (25.6.1)
_SERIES_SWITCH = 1.0


def _upper_planck_integral(z: float) -> float:
    # int_z^inf x^3/(e^x - 1) dx = sum_n Gamma(4, n z)/n^4 (geometric series + DLMF 8.2.2),
    # Gamma(4, y) = 6 e^-y (1 + y + y^2/2 + y^3/6) (DLMF 8.4.8, 8.4.11).
    total = 0.0
    n = 1
    while True:
        y = n * z
        term = 6.0 * exp(-y) * (1.0 + y + y * y / 2.0 + y**3 / 6.0) / n**4
        total += term
        if term <= 1e-18 * total or n > 100_000:
            return total
        n += 1


def _lower_planck_integral(z: float) -> float:
    # int_0^z x^2 [x/(e^x - 1)] dx with x/(e^x - 1) = sum B_k x^k / k! (DLMF 24.2.1, |x| < 2 pi).
    return sum(b * z ** (k + 3) / ((k + 3) * factorial(k)) for k, b in _BERNOULLI.items())


def blackbody_fraction(wavelength_temperature_m_k: float) -> float:
    """Return ``F(0 -> lambda T)``, the black-body energy fraction below ``lambda``.

    Parameters
    ----------
    wavelength_temperature_m_k : float
        Product ``lambda T`` in m K (2897.77 um K = 2.89777e-3 m K); ``inf``
        returns 1.

    Returns
    -------
    float
        ``int_0^lambda e_lambda,b dlambda / (sigma T^4)``, in [0, 1].

    Notes
    -----
    With ``x = c2/(lambda T)`` Lienhard Eq. (1.30) gives
    ``F = (15/pi^4) int_z^inf x^3/(e^x - 1) dx``, ``z = c2/(lambda T)``.
    The normalisation ``pi^4/15 = Gamma(4) zeta(4)`` follows from DLMF
    (25.5.1) and (25.6.1). For ``z >= 1`` the integral is summed as
    ``sum_n Gamma(4, n z)/n^4`` (expanding ``1/(e^x-1) = sum e^{-n x}`` and
    using DLMF (8.2.2), (8.4.8), (8.4.11)), i.e.
    ``e^{-n z}(z^3/n + 3 z^2/n^2 + 6 z/n^3 + 6/n^4)``. For ``z < 1`` the
    complement ``int_0^z`` is integrated term by term from the Bernoulli
    generating function, DLMF (24.2.1), with ``B_0..B_16`` of DLMF
    Table 24.2.1; the first omitted term is below 3e-15 of the integral.
    The derivation is this module's; it is checked against direct
    quadrature of Eq. (1.30) in the tests.

    References
    ----------
    .. [1] LIENHARD-AHTT-6.00. Eq. (1.30), p. 30.
    .. [2] DLMF-25.5. Eqs. (25.5.1), (25.6.1), (8.2.2), (8.4.8), (8.4.11),
           (24.2.1), Table 24.2.1.
    """

    product = real(wavelength_temperature_m_k, "wavelength_temperature_m_k") if \
        wavelength_temperature_m_k != inf else inf
    if product == inf:
        return 1.0
    if product <= 0.0:
        raise DomainError("wavelength_temperature_m_k must be positive.",
                          field="wavelength_temperature_m_k")
    z = SECOND_RADIATION_CONSTANT_M_K / product
    if z >= _SERIES_SWITCH:
        value = _upper_planck_integral(z) / _PLANCK_INTEGRAL
    else:
        value = 1.0 - _lower_planck_integral(z) / _PLANCK_INTEGRAL
    return min(1.0, max(0.0, value))


def blackbody_band_fraction(lower_wavelength_m: float, upper_wavelength_m: float,
                            temperature_k: float) -> float:
    """Return the black-body energy fraction between two wavelengths.

    Parameters
    ----------
    lower_wavelength_m, upper_wavelength_m : float
        Band limits; ``0 <= lower < upper``, ``upper`` may be ``inf``.
    temperature_k : float
        Absolute temperature, positive.

    Returns
    -------
    float
        ``F(0 -> upper T) - F(0 -> lower T)``.

    References
    ----------
    .. [1] LIENHARD-AHTT-6.00. Eq. (1.30), p. 30.
    """

    temperature = _positive(temperature_k, "temperature_k")
    lower = _nonnegative(lower_wavelength_m, "lower_wavelength_m")
    upper = inf if upper_wavelength_m == inf else _positive(upper_wavelength_m, "upper_wavelength_m")
    if not upper > lower:
        raise DomainError("upper_wavelength_m must exceed lower_wavelength_m.", field="upper_wavelength_m")
    low = 0.0 if lower == 0.0 else blackbody_fraction(lower * temperature)
    high = 1.0 if upper == inf else blackbody_fraction(upper * temperature)
    return high - low


def wien_peak_wavelength(temperature_k: float) -> float:
    """Return the wavelength (m) of maximum ``e_lambda,b`` at ``temperature_k``.

    References
    ----------
    .. [1] LIENHARD-AHTT-6.00. Wien's displacement law, Eq. (1.29), p. 30.
    """

    return WIEN_DISPLACEMENT_CONSTANT_M_K / _positive(temperature_k, "temperature_k")


# ---------------------------------------------------------------------------
# View factors
# ---------------------------------------------------------------------------


def view_factor_parallel_strips(width_m: float, separation_m: float) -> float:
    """Return ``F1-2 = F2-1`` between two directly opposed infinite strips.

    Parameters
    ----------
    width_m : float
        Common strip width ``w``.
    separation_m : float
        Distance ``h`` between the strips.

    References
    ----------
    .. [1] LIENHARD-AHTT-6.00. Table 10.2, configuration 1, p. 561.
    """

    ratio = _positive(separation_m, "separation_m") / _positive(width_m, "width_m")
    return sqrt(1.0 + ratio * ratio) - ratio


def view_factor_hinged_strips(included_angle_rad: float) -> float:
    """Return ``F1-2 = F2-1`` between equal infinite strips sharing an edge.

    Parameters
    ----------
    included_angle_rad : float
        Angle ``alpha`` between the strips, in (0, pi].

    References
    ----------
    .. [1] LIENHARD-AHTT-6.00. Table 10.2, configuration 2, p. 561.
    """

    alpha = real(included_angle_rad, "included_angle_rad")
    if not 0.0 < alpha <= pi:
        raise DomainError("included_angle_rad must be in (0, pi].", field="included_angle_rad")
    return 1.0 - sin(alpha / 2.0)


def view_factor_perpendicular_strips(width_1_m: float, width_2_m: float) -> float:
    """Return ``F1-2`` from strip 1 (width ``w``) to a perpendicular strip 2 (``h``).

    The infinite strips share one edge at 90 degrees.

    References
    ----------
    .. [1] LIENHARD-AHTT-6.00. Table 10.2, configuration 3, p. 561.
    """

    ratio = _positive(width_2_m, "width_2_m") / _positive(width_1_m, "width_1_m")
    return 0.5 * (1.0 + ratio - sqrt(1.0 + ratio * ratio))


def view_factor_triangular_duct(side_1_m: float, side_2_m: float, side_3_m: float) -> float:
    """Return ``F1-2`` inside a long duct of triangular cross-section.

    Parameters
    ----------
    side_1_m, side_2_m, side_3_m : float
        Widths of the three flat walls; they must form a triangle.

    Returns
    -------
    float
        ``(A1 + A2 - A3) / (2 A1)``.

    References
    ----------
    .. [1] LIENHARD-AHTT-6.00. Table 10.2, configuration 4, p. 561; Example
           10.10, p. 579.
    """

    a1 = _positive(side_1_m, "side_1_m")
    a2 = _positive(side_2_m, "side_2_m")
    a3 = _positive(side_3_m, "side_3_m")
    if not (a1 < a2 + a3 and a2 < a1 + a3 and a3 < a1 + a2):
        raise DomainError("The three sides do not form a triangle.", field="side_3_m")
    return (a1 + a2 - a3) / (2.0 * a1)


def view_factor_strip_to_cylinder(radius_m: float, near_offset_m: float, far_offset_m: float,
                                  height_m: float) -> float:
    """Return ``F1-2`` from an infinite strip to a parallel infinite cylinder.

    Parameters
    ----------
    radius_m : float
        Cylinder radius ``r``.
    near_offset_m, far_offset_m : float
        Strip edges ``a < b`` measured along the strip from the foot of the
        perpendicular through the cylinder axis.
    height_m : float
        Axis-to-strip-plane distance ``c``; must exceed ``r``.

    Returns
    -------
    float
        ``r/(b - a) [atan(b/c) - atan(a/c)]``.

    References
    ----------
    .. [1] LIENHARD-AHTT-6.00. Table 10.2, configuration 5, p. 561.
    """

    r = _positive(radius_m, "radius_m")
    a = real(near_offset_m, "near_offset_m")
    b = real(far_offset_m, "far_offset_m")
    c = _positive(height_m, "height_m")
    if not b > a:
        raise DomainError("far_offset_m must exceed near_offset_m.", field="far_offset_m")
    if not c > r:
        raise DomainError("The cylinder must not reach the strip plane (height_m > radius_m).",
                          field="height_m")
    return r / (b - a) * (atan(b / c) - atan(a / c))


def view_factor_parallel_cylinders(diameter_m: float, gap_m: float) -> float:
    """Return ``F1-2 = F2-1`` between two parallel infinite cylinders of equal diameter.

    Parameters
    ----------
    diameter_m : float
        Cylinder diameter ``D``.
    gap_m : float
        Surface-to-surface gap ``s``, nonnegative.

    References
    ----------
    .. [1] LIENHARD-AHTT-6.00. Table 10.2, configuration 6, p. 561.
    """

    x = 1.0 + _nonnegative(gap_m, "gap_m") / _positive(diameter_m, "diameter_m")
    return (sqrt(x * x - 1.0) + asin(1.0 / x) - x) / pi


def view_factor_concentric_cylinders(inner_radius_m: float, outer_radius_m: float) -> tuple[
        tuple[float, float], tuple[float, float]]:
    """Return ``((F11, F12), (F21, F22))`` for long concentric cylinders.

    Surface 1 is the outside of the inner cylinder, 2 the inside of the outer.

    References
    ----------
    .. [1] LIENHARD-AHTT-6.00. Table 10.2, configuration 7, p. 561.
    """

    r1, r2 = _nested_radii(inner_radius_m, outer_radius_m)
    f21 = r1 / r2
    return ((0.0, 1.0), (f21, 1.0 - f21))


def view_factor_concentric_spheres(inner_radius_m: float, outer_radius_m: float) -> tuple[
        tuple[float, float], tuple[float, float]]:
    """Return ``((F11, F12), (F21, F22))`` for concentric spheres.

    References
    ----------
    .. [1] LIENHARD-AHTT-6.00. Table 10.3, configuration 4, p. 562.
    """

    r1, r2 = _nested_radii(inner_radius_m, outer_radius_m)
    f21 = (r1 / r2) ** 2
    return ((0.0, 1.0), (f21, 1.0 - f21))


def _nested_radii(inner: object, outer: object) -> tuple[float, float]:
    r1 = _positive(inner, "inner_radius_m")
    r2 = _positive(outer, "outer_radius_m")
    if not r2 >= r1:
        raise DomainError("outer_radius_m must not be smaller than inner_radius_m.", field="outer_radius_m")
    return r1, r2


def view_factor_parallel_rectangles(side_a_m: float, side_b_m: float, separation_m: float) -> float:
    """Return ``F1-2`` between identical, parallel, directly opposed rectangles.

    Parameters
    ----------
    side_a_m, side_b_m : float
        Rectangle sides ``a`` and ``b``.
    separation_m : float
        Plane spacing ``c``.

    References
    ----------
    .. [1] LIENHARD-AHTT-6.00. Table 10.3, configuration 1, p. 562.
    .. [2] HOWELL-CATALOG-C. Factor C-11.
    """

    x = _positive(side_a_m, "side_a_m") / _positive(separation_m, "separation_m")
    y = _positive(side_b_m, "side_b_m") / separation_m
    x2, y2 = x * x, y * y
    braces = (0.5 * log((1.0 + x2) * (1.0 + y2) / (1.0 + x2 + y2))
              + x * sqrt(1.0 + y2) * atan(x / sqrt(1.0 + y2))
              + y * sqrt(1.0 + x2) * atan(y / sqrt(1.0 + x2))
              - x * atan(x) - y * atan(y))
    return 2.0 / (pi * x * y) * braces


def view_factor_perpendicular_rectangles(width_1_m: float, height_2_m: float, common_edge_m: float) -> float:
    """Return ``F1-2`` between perpendicular rectangles sharing a full edge.

    Parameters
    ----------
    width_1_m : float
        Width ``w`` of rectangle 1 (normal to the common edge).
    height_2_m : float
        Height ``h`` of rectangle 2 (normal to the common edge).
    common_edge_m : float
        Length ``l`` of the shared edge.

    References
    ----------
    .. [1] LIENHARD-AHTT-6.00. Table 10.3, configuration 2, p. 562.
    .. [2] HOWELL-CATALOG-C. Factor C-14.
    """

    length = _positive(common_edge_m, "common_edge_m")
    w = _positive(width_1_m, "width_1_m") / length
    h = _positive(height_2_m, "height_2_m") / length
    w2, h2 = w * w, h * h
    s = h2 + w2
    log_term = (log((1.0 + w2) * (1.0 + h2) / (1.0 + s))
                + w2 * log(w2 * (1.0 + s) / ((1.0 + w2) * s))
                + h2 * log(h2 * (1.0 + s) / ((1.0 + h2) * s)))
    return (w * atan(1.0 / w) + h * atan(1.0 / h) - sqrt(s) * atan(1.0 / sqrt(s))
            + 0.25 * log_term) / (pi * w)


def view_factor_coaxial_disks(radius_1_m: float, radius_2_m: float, separation_m: float) -> float:
    """Return ``F1-2`` from disk 1 to a parallel coaxial disk 2.

    Parameters
    ----------
    radius_1_m, radius_2_m : float
        Disk radii ``r1`` (emitter) and ``r2``.
    separation_m : float
        Distance between the disk planes.

    Notes
    -----
    ``R_i = r_i/h``, ``X = 1 + (1 + R2^2)/R1^2``,
    ``F1-2 = (X - sqrt(X^2 - 4 (R2/R1)^2)) / 2``.

    References
    ----------
    .. [1] LIENHARD-AHTT-6.00. Table 10.3, configuration 3, p. 562.
    .. [2] HOWELL-CATALOG-C. Factor C-41 (C-40 for equal radii).
    """

    distance = _positive(separation_m, "separation_m")
    r1 = _positive(radius_1_m, "radius_1_m") / distance
    r2 = _positive(radius_2_m, "radius_2_m") / distance
    x = 1.0 + (1.0 + r2 * r2) / (r1 * r1)
    ratio = r2 / r1
    return 0.5 * (x - sqrt(max(0.0, x * x - 4.0 * ratio * ratio)))


def view_factor_cylinder_interior_to_itself(radius_m: float, length_m: float) -> float:
    """Return ``F2-2`` of the inside of a finite right circular cylinder to itself.

    Notes
    -----
    ``H = h/(2 r)``, ``F2-2 = 1 + H - sqrt(1 + H^2)``.

    References
    ----------
    .. [1] HOWELL-CATALOG-C. Factor C-78.
    """

    h = _positive(length_m, "length_m") / (2.0 * _positive(radius_m, "radius_m"))
    return 1.0 + h - sqrt(1.0 + h * h)


def view_factor_frustum_to_base(top_radius_m: float, base_radius_m: float, height_m: float) -> float:
    """Return ``F1-2`` from the interior of a conical frustum to its base disk.

    Parameters
    ----------
    top_radius_m : float
        Radius ``r3`` of the other (open) end.
    base_radius_m : float
        Radius ``r`` of the base disk ``A2``.
    height_m : float
        Axial height ``h``.

    Notes
    -----
    ``H = h/r3``, ``R = r/r3``, ``X = 1 + R^2 + H^2``,
    ``F1-2 = (2 R^2 - X + sqrt(X^2 - 4 R^2)) / (2 sqrt(X - 2R) (1 + R))``.

    References
    ----------
    .. [1] HOWELL-CATALOG-C. Factor C-111.
    """

    r3 = _positive(top_radius_m, "top_radius_m")
    big_r = _positive(base_radius_m, "base_radius_m") / r3
    big_h = _positive(height_m, "height_m") / r3
    x = 1.0 + big_r * big_r + big_h * big_h
    return ((2.0 * big_r * big_r - x + sqrt(x * x - 4.0 * big_r * big_r))
            / (2.0 * sqrt(x - 2.0 * big_r) * (1.0 + big_r)))


def reciprocal_view_factor(area_from_m2: float, view_factor: float, area_to_m2: float) -> float:
    """Return ``F2-1 = A1 F1-2 / A2`` (view-factor reciprocity).

    References
    ----------
    .. [1] LIENHARD-AHTT-6.00. Eq. (10.15), p. 557.
    """

    factor = _fraction(view_factor, "view_factor")
    result = _positive(area_from_m2, "area_from_m2") * factor / _positive(area_to_m2, "area_to_m2")
    if result > 1.0 + 1e-12:
        raise DomainError("Reciprocity gives a view factor above one; check the areas.", field="view_factor")
    return min(result, 1.0)


def view_factor_matrix_residuals(areas_m2: Sequence[float],
                                 view_factors: Sequence[Sequence[float]]) -> tuple[float, float]:
    """Return the largest summation and reciprocity residuals of an enclosure.

    Parameters
    ----------
    areas_m2 : sequence of float
        Surface areas ``A_i``.
    view_factors : sequence of sequence of float
        Square matrix ``F[i][j] = F_i-j``.

    Returns
    -------
    tuple of float
        ``max_i |sum_j F_ij - 1|`` and
        ``max_ij |A_i F_ij - A_j F_ji| / max(A_i, A_j)``.

    References
    ----------
    .. [1] LIENHARD-AHTT-6.00. Eqs. (10.12) and (10.15), pp. 555-557.
    """

    areas, matrix = _enclosure_geometry(areas_m2, view_factors)
    return _residuals(areas, matrix)


def _enclosure_geometry(areas_m2: Sequence[float], view_factors: Sequence[Sequence[float]]) -> tuple[
        tuple[float, ...], tuple[tuple[float, ...], ...]]:
    areas = tuple(_positive(a, f"areas_m2[{i}]") for i, a in enumerate(areas_m2))
    n = len(areas)
    if n < 2:
        raise InputError("An enclosure needs at least two surfaces.", field="areas_m2")
    rows = tuple(view_factors)
    if len(rows) != n or any(len(row) != n for row in rows):
        raise InputError("view_factors must be an n x n matrix matching areas_m2.", field="view_factors")
    matrix = tuple(tuple(real(v, f"view_factors[{i}][{j}]") for j, v in enumerate(row))
                   for i, row in enumerate(rows))
    for i, row in enumerate(matrix):
        for j, value in enumerate(row):
            if not -1e-12 <= value <= 1.0 + 1e-12:
                raise DomainError("View factors must lie in [0, 1].", field=f"view_factors[{i}][{j}]")
    return areas, matrix


def _residuals(areas: tuple[float, ...], matrix: tuple[tuple[float, ...], ...]) -> tuple[float, float]:
    n = len(areas)
    summation = max(abs(sum(row) - 1.0) for row in matrix)
    reciprocity = max((abs(areas[i] * matrix[i][j] - areas[j] * matrix[j][i]) / max(areas[i], areas[j])
                       for i in range(n) for j in range(n)), default=0.0)
    return summation, reciprocity


@dataclass(frozen=True, slots=True)
class FrustumEnclosure:
    """Three-surface enclosure of a conical frustum (e.g. a nozzle extension).

    Attributes
    ----------
    areas_m2 : tuple of float
        ``(A1, A2, A3)``: small opening disk, large opening disk, conical wall.
    view_factors : tuple of tuple of float
        ``F[i][j]`` for ``i, j`` in (small opening, large opening, wall).
    """

    areas_m2: tuple[float, float, float]
    view_factors: tuple[tuple[float, float, float], ...]


def frustum_enclosure_view_factors(small_radius_m: float, large_radius_m: float,
                                   length_m: float) -> FrustumEnclosure:
    """Return the view factors among the two openings and the wall of a frustum.

    Parameters
    ----------
    small_radius_m, large_radius_m : float
        End radii ``r1 <= r2`` (equal radii give a cylinder).
    length_m : float
        Axial length.

    Returns
    -------
    FrustumEnclosure
        Areas and the full 3 x 3 view-factor matrix.

    Notes
    -----
    ``F1-2`` is the coaxial-disk factor (Lienhard Table 10.3 case 3, Howell
    C-41); the remaining factors follow from summation (Eq. (10.12)) and
    reciprocity (Eq. (10.15)) as in Lienhard Examples 10.3 and 10.4, where a
    plane is laid across the open end of a conical shield. The flat ends do
    not see themselves. Wall-to-base values agree with Howell C-111.

    References
    ----------
    .. [1] LIENHARD-AHTT-6.00. Eqs. (10.12), (10.15), Table 10.3 case 3,
           Examples 10.3-10.4, pp. 555-567.
    .. [2] HOWELL-CATALOG-C. Factors C-41 and C-111.
    """

    r1 = _positive(small_radius_m, "small_radius_m")
    r2 = _positive(large_radius_m, "large_radius_m")
    length = _positive(length_m, "length_m")
    if r2 < r1:
        raise DomainError("large_radius_m must not be smaller than small_radius_m.", field="large_radius_m")
    a1, a2 = pi * r1 * r1, pi * r2 * r2
    a3 = pi * (r1 + r2) * sqrt(length * length + (r2 - r1) ** 2)
    f12 = view_factor_coaxial_disks(r1, r2, length)
    f13 = 1.0 - f12
    f21 = a1 * f12 / a2
    f23 = 1.0 - f21
    f31 = a1 * f13 / a3
    f32 = a2 * f23 / a3
    f33 = 1.0 - f31 - f32
    return FrustumEnclosure((a1, a2, a3), ((0.0, f12, f13), (f21, 0.0, f23), (f31, f32, f33)))


# ---------------------------------------------------------------------------
# Gray diffuse surfaces
# ---------------------------------------------------------------------------


def two_surface_heat_rate(*, temperature_1_k: float, temperature_2_k: float, area_1_m2: float,
                          area_2_m2: float, emissivity_1: float, emissivity_2: float,
                          view_factor_12: float) -> float:
    """Return ``Q_net,1-2`` (W) for two gray diffuse surfaces that see only each other.

    Notes
    -----
    ``sigma (T1^4 - T2^4) / [(1-e1)/(e1 A1) + 1/(A1 F12) + (1-e2)/(e2 A2)]``.

    References
    ----------
    .. [1] LIENHARD-AHTT-6.00. Eq. (10.23), p. 568.
    """

    t1 = _nonnegative(temperature_1_k, "temperature_1_k")
    t2 = _nonnegative(temperature_2_k, "temperature_2_k")
    a1 = _positive(area_1_m2, "area_1_m2")
    a2 = _positive(area_2_m2, "area_2_m2")
    e1 = _emissivity(emissivity_1, "emissivity_1")
    e2 = _emissivity(emissivity_2, "emissivity_2")
    f12 = _fraction(view_factor_12, "view_factor_12")
    if f12 == 0.0:
        return 0.0
    resistance = (1.0 - e1) / (e1 * a1) + 1.0 / (a1 * f12) + (1.0 - e2) / (e2 * a2)
    return STEFAN_BOLTZMANN_W_M2_K4 * (t1**4 - t2**4) / resistance


def parallel_plate_transfer_factor(emissivity_1: float, emissivity_2: float) -> float:
    """Return the transfer factor ``1/(1/e1 + 1/e2 - 1)`` of infinite parallel plates.

    Also valid for specularly reflecting surfaces (Lienhard Eq. (10.32)).

    References
    ----------
    .. [1] LIENHARD-AHTT-6.00. Eqs. (10.25), (10.32), pp. 569, 571.
    """

    e1 = _emissivity(emissivity_1, "emissivity_1")
    e2 = _emissivity(emissivity_2, "emissivity_2")
    return 1.0 / (1.0 / e1 + 1.0 / e2 - 1.0)


def enclosed_body_transfer_factor(emissivity_1: float, emissivity_2: float, area_ratio_12: float) -> float:
    """Return ``F1-2`` for a convex gray body 1 inside a gray enclosure 2.

    Parameters
    ----------
    emissivity_1, emissivity_2 : float
        Emittances of the inner body and the enclosure.
    area_ratio_12 : float
        ``A1/A2`` in (0, 1]; ``A1 << A2`` gives ``F1-2 -> e1`` (Eq. (10.30)).

    Notes
    -----
    ``1 / [1/e1 + (A1/A2)(1/e2 - 1)]``; valid only when body 1 does not
    view itself.

    References
    ----------
    .. [1] LIENHARD-AHTT-6.00. Eqs. (10.27), (10.30), p. 570.
    """

    e1 = _emissivity(emissivity_1, "emissivity_1")
    e2 = _emissivity(emissivity_2, "emissivity_2")
    ratio = real(area_ratio_12, "area_ratio_12")
    if not 0.0 <= ratio <= 1.0:
        raise DomainError("area_ratio_12 must be in [0, 1].", field="area_ratio_12")
    return 1.0 / (1.0 / e1 + ratio * (1.0 / e2 - 1.0))


@dataclass(frozen=True, slots=True)
class RadiationShield:
    """A thin, isothermal radiation shield between two nested surfaces.

    Attributes
    ----------
    area_m2 : float
        Shield area (both faces taken equal for a thin shield).
    inner_emissivity, outer_emissivity : float
        Emittances of the face toward the hot (inner) and the outer surface.
    """

    area_m2: float
    inner_emissivity: float
    outer_emissivity: float


@dataclass(frozen=True, slots=True)
class ShieldedExchange:
    """Net exchange through nested gray surfaces and shields.

    Attributes
    ----------
    heat_rate_w : float
        Net heat rate from the inner to the outer surface (negative: inward).
    shield_temperatures_k : tuple of float
        Shield temperatures from inner to outer.
    total_resistance_1_m2 : float
        Sum of the surface and space resistances, in 1/m^2.
    """

    heat_rate_w: float
    shield_temperatures_k: tuple[float, ...]
    total_resistance_1_m2: float


def shielded_heat_rate(*, inner_temperature_k: float, outer_temperature_k: float, inner_area_m2: float,
                       inner_emissivity: float, outer_area_m2: float | None, outer_emissivity: float,
                       shields: Sequence[RadiationShield] = ()) -> ShieldedExchange:
    """Return the net radiation between nested gray surfaces with thin shields.

    Parameters
    ----------
    inner_temperature_k, outer_temperature_k : float
        Temperatures of the innermost body and the outer enclosure.
    inner_area_m2, inner_emissivity : float
        Inner body area and emittance.
    outer_area_m2 : float or None
        Outer enclosure area; ``None`` treats it as much larger than every
        other area, so its surface resistance is neglected (Example 10.8).
    outer_emissivity : float
        Outer enclosure emittance (unused when ``outer_area_m2`` is None).
    shields : sequence of RadiationShield
        Shields ordered from the inner body outward; areas must not
        decrease.

    Returns
    -------
    ShieldedExchange

    Notes
    -----
    Each surface sees the next one completely (``F = 1``) and each gap adds
    the space resistance ``1/A`` of its inner surface. For one shield with
    equal face emittances this is Lienhard Eq. (10.31); with no shields it is
    Eq. (10.26). Shield temperatures are the circuit node potentials.

    References
    ----------
    .. [1] LIENHARD-AHTT-6.00. Eqs. (10.26), (10.31), Example 10.8,
           pp. 569-572.
    """

    t1 = _nonnegative(inner_temperature_k, "inner_temperature_k")
    t2 = _nonnegative(outer_temperature_k, "outer_temperature_k")
    a1 = _positive(inner_area_m2, "inner_area_m2")
    e1 = _emissivity(inner_emissivity, "inner_emissivity")
    areas = [a1]
    segments: list[tuple[str, float]] = [("surface", (1.0 - e1) / (e1 * a1))]
    previous = a1
    for index, shield in enumerate(shields):
        area = _positive(shield.area_m2, f"shields[{index}].area_m2")
        if area < previous * (1.0 - 1e-12):
            raise DomainError("Shield areas must not decrease outward.", field=f"shields[{index}].area_m2")
        ein = _emissivity(shield.inner_emissivity, f"shields[{index}].inner_emissivity")
        eout = _emissivity(shield.outer_emissivity, f"shields[{index}].outer_emissivity")
        segments.append(("space", 1.0 / previous))
        segments.append(("shield_in", (1.0 - ein) / (ein * area)))
        segments.append(("shield_out", (1.0 - eout) / (eout * area)))
        areas.append(area)
        previous = area
    segments.append(("space", 1.0 / previous))
    if outer_area_m2 is not None:
        a2 = _positive(outer_area_m2, "outer_area_m2")
        if a2 < previous * (1.0 - 1e-12):
            raise DomainError("outer_area_m2 must not be smaller than the surfaces it encloses.",
                              field="outer_area_m2")
        e2 = _emissivity(outer_emissivity, "outer_emissivity")
        segments.append(("surface", (1.0 - e2) / (e2 * a2)))
    else:
        _emissivity(outer_emissivity, "outer_emissivity")
    total = sum(value for _, value in segments)
    q = STEFAN_BOLTZMANN_W_M2_K4 * (t1**4 - t2**4) / total
    potential = STEFAN_BOLTZMANN_W_M2_K4 * t1**4
    temperatures = []
    for kind, value in segments:
        potential -= q * value
        if kind == "shield_in":
            if potential < -1e-9 * STEFAN_BOLTZMANN_W_M2_K4 * max(t1, t2) ** 4:
                raise ConvergenceError("Negative shield emissive power.")
            temperatures.append((max(potential, 0.0) / STEFAN_BOLTZMANN_W_M2_K4) ** 0.25)
    return ShieldedExchange(q, tuple(temperatures), total)


@dataclass(frozen=True, slots=True)
class GraySurface:
    """One isothermal, diffuse, gray, opaque surface of an enclosure.

    Exactly one of ``temperature_k`` and ``net_heat_flux_w_m2`` is given.

    Attributes
    ----------
    name : str
        Label.
    area_m2 : float
        Surface area (or area per unit length for 2-D enclosures).
    emissivity : float
        Total hemispherical emittance in (0, 1].
    temperature_k : float or None
        Prescribed temperature.
    net_heat_flux_w_m2 : float or None
        Prescribed net flux leaving the surface, ``q = B - H``; zero is an
        adiabatic (refractory) wall.
    """

    name: str
    area_m2: float
    emissivity: float
    temperature_k: float | None = None
    net_heat_flux_w_m2: float | None = None


@dataclass(frozen=True, slots=True)
class EnclosureSolution:
    """Result of the net-radiation method.

    Attributes
    ----------
    names : tuple of str
        Surface labels in input order.
    radiosities_w_m2, irradiations_w_m2 : tuple of float
        ``B_i`` and ``H_i = sum_j F_ij B_j``.
    net_heat_fluxes_w_m2, net_heat_rates_w : tuple of float
        ``q_i = B_i - H_i`` (leaving positive) and ``Q_i = A_i q_i``.
    temperatures_k : tuple of float
        Prescribed or solved surface temperatures.
    energy_residual_w : float
        ``sum_i Q_i``; zero for a closed enclosure up to round-off.
    summation_residual, reciprocity_residual : float
        View-factor checks of the input (see view_factor_matrix_residuals).
    """

    names: tuple[str, ...]
    radiosities_w_m2: tuple[float, ...]
    irradiations_w_m2: tuple[float, ...]
    net_heat_fluxes_w_m2: tuple[float, ...]
    net_heat_rates_w: tuple[float, ...]
    temperatures_k: tuple[float, ...]
    energy_residual_w: float
    summation_residual: float
    reciprocity_residual: float


def _solve_linear(matrix: list[list[float]], rhs: list[float]) -> list[float]:
    n = len(rhs)
    a = [row[:] + [rhs[i]] for i, row in enumerate(matrix)]
    scale = max(abs(v) for row in matrix for v in row)
    for col in range(n):
        pivot = max(range(col, n), key=lambda r: abs(a[r][col]))
        if abs(a[pivot][col]) <= 1e-13 * scale:
            raise DomainError("The enclosure equations are singular; prescribe at least one temperature "
                              "and check the view factors.", field="surfaces")
        a[col], a[pivot] = a[pivot], a[col]
        for row in range(col + 1, n):
            factor = a[row][col] / a[col][col]
            if factor:
                for k in range(col, n + 1):
                    a[row][k] -= factor * a[col][k]
    x = [0.0] * n
    for row in range(n - 1, -1, -1):
        x[row] = (a[row][n] - sum(a[row][k] * x[k] for k in range(row + 1, n))) / a[row][row]
    return x


def solve_gray_enclosure(surfaces: Sequence[GraySurface], view_factors: Sequence[Sequence[float]], *,
                         view_factor_tolerance: float = 1e-6) -> EnclosureSolution:
    """Solve radiation exchange in an enclosure of N gray diffuse surfaces.

    Parameters
    ----------
    surfaces : sequence of GraySurface
        Every surface of the enclosure, with a temperature or a net flux.
    view_factors : sequence of sequence of float
        ``F[i][j] = F_i-j`` including self-viewing terms ``F_ii``.
    view_factor_tolerance : float, default 1e-6
        Largest accepted summation and reciprocity residuals; larger ones
        raise :class:`DomainError` (fail-closed).

    Returns
    -------
    EnclosureSolution

    Notes
    -----
    For a prescribed temperature, Eq. (10.44):
    ``B_i - (1 - e_i) sum_j F_ij B_j = e_i sigma T_i^4``. For a prescribed
    flux, Eq. (10.39) with Eq. (10.41) gives
    ``B_i - sum_j F_ij B_j = q_i``, and the temperature follows from
    ``sigma T_i^4 = B_i + q_i (1 - e_i)/e_i``. The medium is transparent
    and every property is uniform over each surface (Lienhard p. 577).

    References
    ----------
    .. [1] LIENHARD-AHTT-6.00. Eqs. (10.39)-(10.44), pp. 577-579.
    """

    items = tuple(surfaces)
    areas, matrix = _enclosure_geometry([s.area_m2 for s in items], view_factors)
    tolerance = _positive(view_factor_tolerance, "view_factor_tolerance")
    summation, reciprocity = _residuals(areas, matrix)
    if summation > tolerance or reciprocity > tolerance:
        raise DomainError(
            f"View factors violate summation ({summation:.3g}) or reciprocity ({reciprocity:.3g}) "
            f"beyond view_factor_tolerance={tolerance:g}.", field="view_factors")
    n = len(items)
    system: list[list[float]] = []
    rhs: list[float] = []
    emissivities = []
    any_temperature = False
    for i, surface in enumerate(items):
        e = _emissivity(surface.emissivity, f"surfaces[{i}].emissivity")
        emissivities.append(e)
        has_t = surface.temperature_k is not None
        has_q = surface.net_heat_flux_w_m2 is not None
        if has_t == has_q:
            raise InputError("Give exactly one of temperature_k and net_heat_flux_w_m2.",
                             field=f"surfaces[{i}]")
        if has_t:
            any_temperature = True
            t = _nonnegative(surface.temperature_k, f"surfaces[{i}].temperature_k")
            system.append([(1.0 if i == j else 0.0) - (1.0 - e) * matrix[i][j] for j in range(n)])
            rhs.append(e * STEFAN_BOLTZMANN_W_M2_K4 * t**4)
        else:
            q = real(surface.net_heat_flux_w_m2, f"surfaces[{i}].net_heat_flux_w_m2")
            system.append([(1.0 if i == j else 0.0) - matrix[i][j] for j in range(n)])
            rhs.append(q)
    if not any_temperature:
        raise InputError("Prescribe the temperature of at least one surface.", field="surfaces")
    radiosity = _solve_linear(system, rhs)
    irradiation = [sum(matrix[i][j] * radiosity[j] for j in range(n)) for i in range(n)]
    flux = [radiosity[i] - irradiation[i] for i in range(n)]
    temperatures = []
    for i, surface in enumerate(items):
        if surface.temperature_k is not None:
            temperatures.append(float(surface.temperature_k))
            continue
        e = emissivities[i]
        q = float(surface.net_heat_flux_w_m2)  # type: ignore[arg-type]
        emissive = radiosity[i] + q * (1.0 - e) / e
        if emissive < 0.0:
            raise DomainError("The prescribed flux needs a negative emissive power.",
                              field=f"surfaces[{i}].net_heat_flux_w_m2")
        temperatures.append((emissive / STEFAN_BOLTZMANN_W_M2_K4) ** 0.25)
        flux[i] = q
    rates = [areas[i] * flux[i] for i in range(n)]
    return EnclosureSolution(tuple(s.name for s in items), tuple(radiosity), tuple(irradiation), tuple(flux),
                             tuple(rates), tuple(temperatures), sum(rates), summation, reciprocity)


def frustum_interior_transfer_factors(*, small_radius_m: float, large_radius_m: float, length_m: float,
                                      wall_emissivity: float) -> tuple[float, float, float]:
    """Return the gray interior wall's transfer factors to the two openings.

    Parameters
    ----------
    small_radius_m, large_radius_m, length_m : float
        Frustum geometry (e.g. the attachment and exit radii of a conical
        nozzle extension and its axial length).
    wall_emissivity : float
        Gray, diffuse emittance of the isothermal interior wall.

    Returns
    -------
    tuple of float
        ``(T_small, T_large, F_wall_wall)`` such that the net flux leaving
        the wall is ``sigma [T_small (Tw^4 - Ts^4) + T_large (Tw^4 - Tl^4)]``
        when each opening behaves as a black surface at ``Ts`` or ``Tl``.

    Notes
    -----
    Applying Eq. (10.42) to the wall (the only gray surface) with the two
    openings black gives ``T_k = e F_w-k / (1 - (1 - e) F_w-w)``; this
    module derives the closed form, and the tests check it against
    :func:`solve_gray_enclosure`. An isothermal wall is assumed; a real
    extension with an axial temperature gradient should be divided into
    several isothermal bands (Lienhard p. 581).

    References
    ----------
    .. [1] LIENHARD-AHTT-6.00. Eqs. (10.39), (10.42), pp. 577-578.
    .. [2] HOWELL-CATALOG-C. Factor C-41.
    """

    e = _emissivity(wall_emissivity, "wall_emissivity")
    enclosure = frustum_enclosure_view_factors(small_radius_m, large_radius_m, length_m)
    f31, f32, f33 = enclosure.view_factors[2]
    denominator = 1.0 - (1.0 - e) * f33
    return e * f31 / denominator, e * f32 / denominator, f33


# ---------------------------------------------------------------------------
# Radiative-equilibrium wall temperature
# ---------------------------------------------------------------------------


@dataclass(frozen=True, slots=True)
class RadiativeSink:
    """A radiation sink seen by a wall: net flux ``sigma F (Tw^4 - T^4)``.

    Attributes
    ----------
    transfer_factor : float
        Transfer factor ``F`` per unit wall area: the emittance for a small
        or convex surface in a large environment (Eq. (10.30)), or a value
        from :func:`frustum_interior_transfer_factors`.
    temperature_k : float
        Sink temperature.
    """

    transfer_factor: float
    temperature_k: float


@dataclass(frozen=True, slots=True)
class RadiativeEquilibrium:
    """Wall temperature where absorbed and convected heat equal radiated heat.

    Attributes
    ----------
    wall_temperature_k : float
        Solution of the energy balance.
    convective_heat_flux_w_m2 : float
        ``h (T_r - T_w)`` at the solution.
    radiated_heat_flux_w_m2 : float
        ``sigma sum F_k (T_w^4 - T_k^4)`` at the solution.
    residual_w_m2 : float
        Energy-balance residual.
    """

    wall_temperature_k: float
    convective_heat_flux_w_m2: float
    radiated_heat_flux_w_m2: float
    residual_w_m2: float


def radiative_equilibrium_wall_temperature(*, sinks: Sequence[RadiativeSink],
                                           heat_transfer_coefficient_w_m2_k: float = 0.0,
                                           recovery_temperature_k: float | None = None,
                                           imposed_heat_flux_w_m2: float = 0.0) -> RadiativeEquilibrium:
    """Solve ``q_imp + h (T_r - T_w) = sigma sum_k F_k (T_w^4 - T_k^4)`` for ``T_w``.

    Parameters
    ----------
    sinks : sequence of RadiativeSink
        Radiation sinks; at least one with a positive transfer factor. A
        thin nozzle extension cooled on both faces has one sink per face.
    heat_transfer_coefficient_w_m2_k : float, default 0
        Convective coefficient ``h`` from the hot gas to the wall.
    recovery_temperature_k : float, optional
        Gas recovery (adiabatic-wall) temperature; required when ``h > 0``.
    imposed_heat_flux_w_m2 : float, default 0
        Any other absorbed flux (e.g. gas radiation), nonnegative.

    Returns
    -------
    RadiativeEquilibrium

    Notes
    -----
    The left side decreases and the right side increases with ``T_w``, so
    the root is unique; it is bracketed and found by bisection to round-off.
    Lienhard Examples 1.5-1.6 solve the same balance for a thermocouple
    (``F = e`` in a large enclosure, Eq. (1.35)); the transfer-factor form of
    the net exchange is Eq. (10.2). Conduction along and through the wall is
    neglected.

    References
    ----------
    .. [1] LIENHARD-AHTT-6.00. Eqs. (1.34)-(1.35), Examples 1.5-1.6,
           pp. 32-33; Eqs. (10.2), (10.30), pp. 545, 570.
    """

    items = tuple(sinks)
    if not items:
        raise InputError("Give at least one radiative sink.", field="sinks")
    factors = []
    for index, sink in enumerate(items):
        factor = _fraction(sink.transfer_factor, f"sinks[{index}].transfer_factor")
        temperature = _nonnegative(sink.temperature_k, f"sinks[{index}].temperature_k")
        factors.append((factor, temperature))
    total = sum(f for f, _ in factors)
    if total <= 0.0:
        raise DomainError("At least one sink needs a positive transfer factor.", field="sinks")
    h = _nonnegative(heat_transfer_coefficient_w_m2_k, "heat_transfer_coefficient_w_m2_k")
    q_imp = _nonnegative(imposed_heat_flux_w_m2, "imposed_heat_flux_w_m2")
    if h > 0.0 and recovery_temperature_k is None:
        raise InputError("recovery_temperature_k is required when h > 0.", field="recovery_temperature_k")
    t_r = 0.0 if recovery_temperature_k is None else _nonnegative(recovery_temperature_k,
                                                                  "recovery_temperature_k")

    def radiated(t: float) -> float:
        return STEFAN_BOLTZMANN_W_M2_K4 * sum(f * (t**4 - ts**4) for f, ts in factors)

    def balance(t: float) -> float:
        return q_imp + h * (t_r - t) - radiated(t)

    t_sink = max(ts for _, ts in factors)
    low = 0.0
    high = max(t_r, (q_imp / (STEFAN_BOLTZMANN_W_M2_K4 * total) + t_sink**4) ** 0.25)
    if balance(high) > 0.0:
        high *= 1.0 + 1e-9
        if balance(high) > 0.0:
            raise ConvergenceError("Could not bracket the radiative-equilibrium temperature.")
    for _ in range(400):
        mid = 0.5 * (low + high)
        if balance(mid) > 0.0:
            low = mid
        else:
            high = mid
        if high - low <= 4e-16 * high:
            break
    t_w = 0.5 * (low + high)
    return RadiativeEquilibrium(t_w, h * (t_r - t_w), radiated(t_w), balance(t_w))


# ---------------------------------------------------------------------------
# Gas radiation: mean beam length, WSGG, gas-to-wall exchange
# ---------------------------------------------------------------------------

MEAN_BEAM_LENGTH_FACTORS: Mapping[str, float] = MappingProxyType({
    "cube": 2.0 / 3.0,
    "sphere": 2.0 / 3.0,
    "infinite_cylinder": 1.0,
    "cylinder_length_equals_diameter": 2.0 / 3.0,
    "infinite_parallel_plates": 2.0,
})
"""``L0 / dimension`` for Lienhard Table 10.4 (side, diameter, diameter,
diameter) and the parallel-plate result ``L0 = 2 l`` (spacing) printed
under Eq. (10.58)."""


def mean_beam_length(volume_m3: float, boundary_area_m2: float, *, optically_thick: bool = False) -> float:
    """Return the geometric mean beam length ``L0 = 4 V / A``.

    Parameters
    ----------
    volume_m3 : float
        Gas volume.
    boundary_area_m2 : float
        Irradiated boundary area.
    optically_thick : bool, default False
        Use 3.6 in place of 4 (10 % shorter) for strongly absorbing gas.

    References
    ----------
    .. [1] LIENHARD-AHTT-6.00. Eq. (10.58) and the text below it, p. 589.
    """

    factor = 3.6 if optically_thick else 4.0
    return factor * _positive(volume_m3, "volume_m3") / _positive(boundary_area_m2, "boundary_area_m2")


def mean_beam_length_for_geometry(geometry: str, dimension_m: float) -> float:
    """Return the tabulated mean beam length of a standard gas volume.

    Parameters
    ----------
    geometry : str
        A key of :data:`MEAN_BEAM_LENGTH_FACTORS`.
    dimension_m : float
        Cube side, sphere or cylinder diameter, or plate spacing.

    References
    ----------
    .. [1] LIENHARD-AHTT-6.00. Table 10.4 and Eq. (10.58), pp. 589, 593.
    """

    if geometry not in MEAN_BEAM_LENGTH_FACTORS:
        raise InputError(f"geometry must be one of {sorted(MEAN_BEAM_LENGTH_FACTORS)}.", field="geometry")
    return MEAN_BEAM_LENGTH_FACTORS[geometry] * _positive(dimension_m, "dimension_m")


@dataclass(frozen=True, slots=True)
class WsggSpecies:
    """Published weighted-sum-of-gray-gases coefficients of one absorber.

    Attributes
    ----------
    name : str
        Species label.
    absorption_coefficients : tuple of float
        ``kappa_j`` of the gray gases: atm^-1 m^-1 (pressure-based) for
        gases, m^-1 per unit soot volume fraction for soot.
    weight_polynomials : tuple of tuple of float
        ``b_j,n`` with ``a_j = sum_n b_j,n (T/T_ref)^n`` (ascending powers).
    reference_temperature_k : float
        ``T_ref`` of the fit.
    minimum_temperature_k, maximum_temperature_k : float
        Temperature range of the fitted database.
    minimum_path, maximum_path : float or None
        Range of the path variable of the fitted database (atm m for
        gases), ``None`` when the source does not state it.
    transparent_window : bool
        Whether a clear gas ``a_0 = 1 - sum a_j`` (``kappa_0 = 0``) exists.
    clip_negative_weights : bool
        Whether the source defines ``a_j = max(polynomial, 0)``.
    total_pressure_atm : float
        Total pressure of the fitted database.
    source : str
        Source key and locator.
    """

    name: str
    absorption_coefficients: tuple[float, ...]
    weight_polynomials: tuple[tuple[float, ...], ...]
    reference_temperature_k: float
    minimum_temperature_k: float
    maximum_temperature_k: float
    minimum_path: float | None
    maximum_path: float | None
    transparent_window: bool
    clip_negative_weights: bool
    total_pressure_atm: float
    source: str


BORDBAR_2020_CO2 = WsggSpecies(
    "CO2",
    (3.272772e-2, 4.229655e-1, 4.905367, 1.085440e2),
    ((8.495135e-1, -1.496812, 1.361406, -5.551699e-1, 8.076589e-2),
     (-1.103102e-1, 9.363958e-1, -1.250799, 6.527827e-1, -1.206959e-1),
     (1.731716e-1, -5.174223e-1, 8.256840e-1, -4.998864e-1, 1.008743e-1),
     (3.995426e-2, 1.423006e-1, -1.649481e-1, 5.140768e-2, -3.497246e-3)),
    1200.0, 300.0, 2400.0, 0.01, 60.0, True, False, 1.0,
    "BORDBAR-FRAGA-HOSTIKKA-2020, Table 2 and Eqs. (2)-(3) (pure CO2, pt = 1 atm, Y = 1, "
    "L = 0.01-60 m, T = 300-2400 K)",
)
"""Pure-CO2 WSGG coefficients at 1 atm (Bordbar, Fraga and Hostikka 2020, Table 2)."""

BORDBAR_2020_H2O = WsggSpecies(
    "H2O",
    (8.047859e-2, 9.557208e-1, 8.005283, 7.613186e1),
    ((6.670204e-1, -1.228413, 1.428908, -6.267906e-1, 9.628539e-2),
     (2.343433e-1, -3.192256e-1, 8.867348e-1, -5.927787e-1, 1.185824e-1),
     (-1.793041e-1, 1.683454, -2.136989, 1.020422, -1.723960e-1),
     (3.455969e-1, -7.510442e-1, 6.313180e-1, -2.416500e-1, 3.530972e-2)),
    1200.0, 300.0, 2400.0, 0.01, 60.0, True, False, 1.0,
    "BORDBAR-FRAGA-HOSTIKKA-2020, Table 3 and Eqs. (2)-(3) (pure H2O, pt = 1 atm, Y = 1, "
    "L = 0.01-60 m, T = 300-2400 K)",
)
"""Pure-H2O WSGG coefficients at 1 atm (Bordbar, Fraga and Hostikka 2020, Table 3)."""

SADEGHI_2021_SOOT = WsggSpecies(
    "soot",
    (1.4530e5, 2.0836e6, 7.5475e5, 4.2113e6),
    ((1.8613, -7.7857, 1.2809e1, -1.0158e1, 3.8717, -5.6880e-1),
     (2.5975e-1, -2.9708, 9.6830, -9.7681, 4.1073, -6.3057e-1),
     (-1.1374, 1.0625e1, -2.1665e1, 1.8750e1, -7.4578, 1.1219),
     (-3.8367e-2, 4.2159e-1, -1.4101, 1.7225, -7.5889e-1, 1.1454e-1)),
    1400.0, 300.0, 3000.0, None, None, False, True, 1.0,
    "SADEGHI-2021-FSJ, Table B9 and Eqs. (9)-(10) (Rayleigh soot with the Chang-Charalampopoulos "
    "refractive index, T = 300-3000 K, no transparent window)",
)
"""Soot WSGG coefficients (Sadeghi et al. 2021, Table B9); kappa in m^-1 per unit f_v."""


def wsgg_weights(species: WsggSpecies, temperature_k: float) -> tuple[float, ...]:
    """Return the gray-gas weights ``(a_0, a_1, ..., a_N)`` at ``temperature_k``.

    ``a_0`` (clear gas) is present only when ``species.transparent_window``.

    Raises
    ------
    DomainError
        Outside the fitted temperature range, or when a weight of a model
        without clipping is negative.

    References
    ----------
    .. [1] BORDBAR-FRAGA-HOSTIKKA-2020. Eq. (3) and ``a_0 = 1 - sum a_i``.
    .. [2] SADEGHI-2021-FSJ. Eqs. (10)-(11).
    """

    t = _positive(temperature_k, "temperature_k")
    if not species.minimum_temperature_k <= t <= species.maximum_temperature_k:
        raise DomainError(
            f"{species.name} WSGG data cover {species.minimum_temperature_k:g}-"
            f"{species.maximum_temperature_k:g} K; got {t:g} K.", field="temperature_k")
    tr = t / species.reference_temperature_k
    weights = []
    for row in species.weight_polynomials:
        value = sum(b * tr**n for n, b in enumerate(row))
        if species.clip_negative_weights:
            value = max(value, 0.0)
        elif value < 0.0:
            raise DomainError(f"{species.name} WSGG weight is negative at {t:g} K.", field="temperature_k")
        weights.append(value)
    if species.transparent_window:
        clear = 1.0 - sum(weights)
        if clear < -1e-12:
            raise DomainError(f"{species.name} WSGG weights exceed one at {t:g} K.", field="temperature_k")
        return (max(clear, 0.0), *weights)
    return tuple(weights)


def _path_check(species: WsggSpecies, path: float, allow: bool, warnings: list[WarningMessage]) -> None:
    low, high = species.minimum_path, species.maximum_path
    if low is None or high is None or path == 0.0:
        return
    if low <= path <= high:
        return
    message = (f"{species.name} pressure path length {path:.4g} atm m is outside the fitted "
               f"{low:g}-{high:g} atm m.")
    if not allow:
        raise DomainError(message + " Pass allow_path_extrapolation=True to accept.", field="path_length_m")
    warnings.append(WarningMessage("wsgg_path_extrapolated", message, "warning", "path_length_m"))


def wsgg_species_emissivity(species: WsggSpecies, temperature_k: float, path: float, *,
                            allow_path_extrapolation: bool = False) -> GasRadiation:
    """Return the total emissivity of one absorber, ``sum_j a_j (1 - exp(-kappa_j path))``.

    Parameters
    ----------
    species : WsggSpecies
        Coefficient set.
    temperature_k : float
        Gas temperature.
    path : float
        ``p_a L`` in atm m for gases, ``f_v L`` in m for soot.
    allow_path_extrapolation : bool, default False
        Accept a path outside the fitted range (warning recorded).

    Returns
    -------
    GasRadiation
        ``gray_gas_absorption_coefficients_1_m`` holds ``kappa_j`` per unit
        path variable (atm^-1 m^-1, or m^-1 per unit ``f_v``).

    References
    ----------
    .. [1] BORDBAR-FRAGA-HOSTIKKA-2020. Eq. (2).
    .. [2] SADEGHI-2021-FSJ. Eq. (9).
    """

    value = _nonnegative(path, "path")
    warnings: list[WarningMessage] = []
    _path_check(species, value, allow_path_extrapolation, warnings)
    weights = wsgg_weights(species, temperature_k)
    kappas = ((0.0,) if species.transparent_window else ()) + species.absorption_coefficients
    emissivity = sum(a * -expm1(-k * value) for a, k in zip(weights, kappas, strict=True))
    return GasRadiation(emissivity, kappas, weights, tuple(warnings))


@dataclass(frozen=True, slots=True)
class GasRadiation:
    """Total emissivity or absorptivity of an absorbing mixture.

    Attributes
    ----------
    value : float
        Total emissivity (or absorptivity).
    gray_gas_absorption_coefficients_1_m : tuple of float
        Mixture gray-gas absorption coefficients (m^-1), Sadeghi Eq. (12).
    gray_gas_weights : tuple of float
        Mixture weights, Sadeghi Eq. (13), evaluated at the emitting
        temperature.
    warnings : tuple of WarningMessage
        Extrapolation notices.
    """

    value: float
    gray_gas_absorption_coefficients_1_m: tuple[float, ...]
    gray_gas_weights: tuple[float, ...]
    warnings: tuple[WarningMessage, ...]


def _mixture(temperature_k: float, path_length_m: float, total_pressure_pa: float, h2o: float, co2: float,
             soot: float, allow_pressure: bool, allow_path: bool) -> GasRadiation:
    t = _positive(temperature_k, "temperature_k")
    length = _positive(path_length_m, "path_length_m")
    pressure = _positive(total_pressure_pa, "total_pressure_pa")
    x_w = _fraction(h2o, "h2o_mole_fraction")
    x_c = _fraction(co2, "co2_mole_fraction")
    f_v = _fraction(soot, "soot_volume_fraction")
    if x_w + x_c > 1.0 + 1e-12:
        raise DomainError("Mole fractions sum above one.", field="co2_mole_fraction")
    if x_w == 0.0 and x_c == 0.0 and f_v == 0.0:
        raise InputError("Give at least one absorbing component.", field="h2o_mole_fraction")
    warnings: list[WarningMessage] = []
    if abs(pressure - _ATM_PA) > 1e-9 * _ATM_PA:
        message = (f"WSGG coefficients are fitted at 1 atm total pressure; {pressure:.6g} Pa is used only "
                   "through the partial pressures, without any pressure-broadening correction.")
        if not allow_pressure:
            raise DomainError(message + " Pass allow_pressure_extrapolation=True to accept.",
                              field="total_pressure_pa")
        warnings.append(WarningMessage("wsgg_pressure_extrapolated", message, "critical", "total_pressure_pa"))
    components: list[tuple[tuple[float, ...], tuple[float, ...]]] = []
    for species, fraction in ((BORDBAR_2020_H2O, x_w), (BORDBAR_2020_CO2, x_c)):
        if fraction == 0.0:
            continue
        partial_atm = fraction * pressure / _ATM_PA
        _path_check(species, partial_atm * length, allow_path, warnings)
        weights = wsgg_weights(species, t)
        kappas = (0.0,) + tuple(k * partial_atm for k in species.absorption_coefficients)
        components.append((kappas, weights))
    if f_v > 0.0:
        weights = wsgg_weights(SADEGHI_2021_SOOT, t)
        kappas = tuple(k * f_v for k in SADEGHI_2021_SOOT.absorption_coefficients)
        components.append((kappas, weights))
        warnings.append(WarningMessage(
            "wsgg_soot_path_unstated",
            "The f_v L range of the soot WSGG fit is not stated in the source; its weights sum to "
            f"{sum(weights):.4f} (not exactly one) at {t:g} K.", "info", "soot_volume_fraction"))
    kappa_mix: list[float] = [0.0]
    weight_mix: list[float] = [1.0]
    for kappas, weights in components:
        kappa_mix = [km + k for km in kappa_mix for k in kappas]
        weight_mix = [wm * w for wm in weight_mix for w in weights]
    value = sum(w * -expm1(-k * length) for k, w in zip(kappa_mix, weight_mix, strict=True))
    return GasRadiation(value, tuple(kappa_mix), tuple(weight_mix), tuple(warnings))


def wsgg_mixture_emissivity(*, temperature_k: float, path_length_m: float, total_pressure_pa: float,
                            h2o_mole_fraction: float = 0.0, co2_mole_fraction: float = 0.0,
                            soot_volume_fraction: float = 0.0, allow_pressure_extrapolation: bool = False,
                            allow_path_extrapolation: bool = False) -> GasRadiation:
    """Return the total emissivity of an H2O/CO2/soot mixture along a path.

    Parameters
    ----------
    temperature_k : float
        Isothermal gas temperature; 300-2400 K for H2O/CO2, 300-3000 K for
        soot.
    path_length_m : float
        Path length, e.g. a mean beam length.
    total_pressure_pa : float
        Mixture pressure. The fits are at 1 atm; other values raise unless
        ``allow_pressure_extrapolation`` (critical warning).
    h2o_mole_fraction, co2_mole_fraction : float
        Absorber mole fractions (the rest is taken as non-absorbing).
    soot_volume_fraction : float
        Soot volume fraction ``f_v``.
    allow_pressure_extrapolation, allow_path_extrapolation : bool
        Accept totals other than 1 atm, or ``p_a L`` outside 0.01-60 atm m.

    Returns
    -------
    GasRadiation

    Notes
    -----
    Pure-species coefficients (Bordbar et al. Tables 2-3, Table B9 for
    soot) are combined by superposition: ``kappa_m = sum_i kappa_i,j p_i``
    and ``a_m = prod_i a_i,j`` over all combinations including each
    species' clear gas (Sadeghi et al. Eqs. (12)-(14)), assuming
    uncorrelated spectra. Bordbar et al. evaluated this double-integration
    use of their Tables 2-3 against line-by-line solutions in four 1-D slab
    benchmarks (their Table 4): domain-average errors 2.7-7.3 %, maxima up
    to 21.5 % (radiative flux, case 1). The pure-species tables were fitted
    with self-broadened gas (``Y = 1``) at 1 atm; applying them to diluted
    absorbers through ``p_a L`` neglects the broadening difference.

    References
    ----------
    .. [1] BORDBAR-FRAGA-HOSTIKKA-2020. Eqs. (2)-(3), Tables 2-4.
    .. [2] SADEGHI-2021-FSJ. Eqs. (9)-(14), Table B9.
    """

    return _mixture(temperature_k, path_length_m, total_pressure_pa, h2o_mole_fraction, co2_mole_fraction,
                    soot_volume_fraction, allow_pressure_extrapolation, allow_path_extrapolation)


def wsgg_mixture_absorptivity(*, gas_temperature_k: float, wall_temperature_k: float, path_length_m: float,
                              total_pressure_pa: float, h2o_mole_fraction: float = 0.0,
                              co2_mole_fraction: float = 0.0, allow_pressure_extrapolation: bool = False,
                              allow_path_extrapolation: bool = False) -> GasRadiation:
    """Return the total absorptivity of an H2O/CO2 gas for radiation from a wall.

    Notes
    -----
    Hottel's scaling rule as given by Lienhard Eq. (10.57):
    ``alpha_g = eps_g(p_a L T_w/T_g, p, T_w) (T_g/T_w)^(1/2)``, with
    ``eps_g`` from :func:`wsgg_mixture_emissivity` (every partial pressure
    scaled by ``T_w/T_g``). Lienhard's footnote 6 reports that both theory
    and later scaling work support the exponent 1/2. Soot is not covered by
    this rule in the opened source and is therefore not accepted.

    References
    ----------
    .. [1] LIENHARD-AHTT-6.00. Eq. (10.57) and footnote 6, p. 588.
    .. [2] BORDBAR-FRAGA-HOSTIKKA-2020. Tables 2-3.
    """

    t_g = _positive(gas_temperature_k, "gas_temperature_k")
    t_w = _positive(wall_temperature_k, "wall_temperature_k")
    length = _positive(path_length_m, "path_length_m") * t_w / t_g
    result = _mixture(t_w, length, total_pressure_pa, h2o_mole_fraction, co2_mole_fraction, 0.0,
                      allow_pressure_extrapolation, allow_path_extrapolation)
    return GasRadiation(result.value * sqrt(t_g / t_w), result.gray_gas_absorption_coefficients_1_m,
                        result.gray_gas_weights, result.warnings)


def gas_to_black_wall_heat_flux(*, gas_emissivity: float, gas_absorptivity: float, gas_temperature_k: float,
                                wall_temperature_k: float) -> float:
    """Return ``q = eps_g sigma T_g^4 - alpha_g sigma T_w^4`` (W/m^2) into a black wall.

    For one isothermal gas filling an enclosure with black walls at one
    temperature, the mean beam length carries all geometry (no view factor).

    References
    ----------
    .. [1] LIENHARD-AHTT-6.00. Example 10.11, p. 593.
    """

    eps = _fraction(gas_emissivity, "gas_emissivity")
    alpha = _fraction(gas_absorptivity, "gas_absorptivity")
    t_g = _nonnegative(gas_temperature_k, "gas_temperature_k")
    t_w = _nonnegative(wall_temperature_k, "wall_temperature_k")
    return STEFAN_BOLTZMANN_W_M2_K4 * (eps * t_g**4 - alpha * t_w**4)


def soot_planck_mean_absorption_coefficient(temperature_k: float, soot_volume_fraction: float) -> float:
    """Return the Planck-mean soot absorption coefficient (m^-1) for gray calculations.

    Notes
    -----
    ``kappa_s = f_v (a T^3 + b T^2 + c T + d)`` with ``a = 2.156e-5``,
    ``b = -0.2889``, ``c = 1804``, ``d = -2.012e5`` (fit at ``f_v = 1`` of
    the Rayleigh/Chang-Charalampopoulos spectrum). The source plots it in
    its Fig. 9; this module limits it to the 300-3000 K range of the soot
    WSGG data of the same paper.

    References
    ----------
    .. [1] SADEGHI-2021-FSJ. Eq. (24), Fig. 9.
    """

    t = _positive(temperature_k, "temperature_k")
    if not 300.0 <= t <= 3000.0:
        raise DomainError("The soot Planck-mean fit is limited to 300-3000 K.", field="temperature_k")
    f_v = _fraction(soot_volume_fraction, "soot_volume_fraction")
    return f_v * (2.156e-5 * t**3 - 0.2889 * t**2 + 1804.0 * t - 2.012e5)


def station_gas_emissivity(station: object, path_length_m: float, *, allow_pressure_extrapolation: bool = False,
                           allow_path_extrapolation: bool = False) -> GasRadiation:
    """Return the WSGG emissivity of a rocket nozzle station's H2O/CO2 content.

    Parameters
    ----------
    station : RocketStation
        A station from
        :func:`rocket_propulsion.thermochemistry.rocket_performance.solve_rocket_performance`;
        its ``temperature_k``, ``pressure_pa`` and ``mole_fractions`` are
        used.
    path_length_m : float
        Path length, e.g. :func:`mean_beam_length_for_geometry` with
        ``"infinite_cylinder"`` and the local diameter.
    allow_pressure_extrapolation, allow_path_extrapolation : bool
        See :func:`wsgg_mixture_emissivity`. Rocket chambers (several MPa,
        ~3000-3600 K) lie outside the 1 atm, 300-2400 K fits; temperatures
        outside the fit always raise.

    Returns
    -------
    GasRadiation

    References
    ----------
    .. [1] BORDBAR-FRAGA-HOSTIKKA-2020. Tables 2-3.
    .. [2] SADEGHI-2021-FSJ. Eqs. (12)-(14).
    """

    try:
        fractions = dict(station.mole_fractions)  # type: ignore[attr-defined]
        temperature = station.temperature_k  # type: ignore[attr-defined]
        pressure = station.pressure_pa  # type: ignore[attr-defined]
    except (AttributeError, TypeError, ValueError) as error:
        raise InputError("station must provide temperature_k, pressure_pa and mole_fractions.",
                         field="station") from error
    return wsgg_mixture_emissivity(temperature_k=temperature, path_length_m=path_length_m,
                                   total_pressure_pa=pressure, h2o_mole_fraction=fractions.get("H2O", 0.0),
                                   co2_mole_fraction=fractions.get("CO2", 0.0),
                                   allow_pressure_extrapolation=allow_pressure_extrapolation,
                                   allow_path_extrapolation=allow_path_extrapolation)

