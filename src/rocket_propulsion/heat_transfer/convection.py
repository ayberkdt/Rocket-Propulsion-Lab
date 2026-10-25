"""Convective heat-transfer correlations, with rocket thrust-chamber cooling in focus.

The module collects correlations that were read, equation by equation, from
the sources listed below; every function names its equation and page. It
covers

* dimensionless groups, hydraulic diameter and smooth/rough friction factors;
* internal laminar flow (fully developed values, Graetz entry fits), internal
  turbulent flow (Gnielinski, the Petukhov-type analogy, the historical
  power laws, Sieder-Tate, liquid metals, fully rough pipes), turbulent
  entrance ratios and property-ratio corrections [1]_;
* coolant-side correlations used for regenerative cooling: Taylor's hydrogen
  correlation with its x/D-dependent temperature-ratio exponent [2]_, the
  Hess-Kunz film-temperature equation and its tabulated coefficient ``C_L``
  [2]_ [3]_, the NASA SP-125 Sieder-Tate form and supercritical vapour-film
  equation [4]_, the Ito/Niino curvature multiplier as used by NASA's RTE
  code [3]_ [5]_, and wrappers that evaluate these with the package's
  Helmholtz equations of state for density and heat capacity;
* external flat-plate and cylinder correlations, natural convection on a
  vertical plate, pool boiling (Rohsenow, peak heat flux), laminar film
  condensation [1]_, a burnout-margin screen [6]_, subcooled flow-boiling
  burnout fluxes for ethanol [8]_ and MMH [9]_, and the Hatch-Papell slot
  film-cooling effectiveness [7]_.

Conventions: SI units throughout; Darcy friction factors; temperature ratios
in kelvin. Correlations are empirical. Outside the validity range stated by
the source a function raises :class:`~rocket_propulsion.core.errors.DomainError`
with code ``"outside_correlation_range"``; with ``extrapolate=True`` it returns
the value and a ``WarningMessage`` instead. Where a source states no range,
the result carries an ``"range_not_stated"`` info warning. No transport-property
(viscosity, conductivity) model for cryogenic coolants is bundled: those values
are always caller inputs. Results are not validated against hardware.

References
----------
.. [1] LIENHARD-AHTT-6.00. J. H. Lienhard V and J. H. Lienhard IV, A Heat
       Transfer Textbook, 6th ed., Phlogiston Press, 2024, Chs. 6-9,
       https://ahtt.mit.edu/
.. [2] NASA-TN-D-4332-TAYLOR-1968. M. F. Taylor, NASA TN D-4332, 1968,
       Eqs. (14), (19), (20), Table I, pp. 2-14,
       https://ntrs.nasa.gov/citations/19680004771
.. [3] NASA-CR-678-AEROJET-1967. Aerojet-General Corp., NASA CR-678, 1967,
       Eqs. (1)-(4), pp. 2-4, https://ntrs.nasa.gov/citations/19670008827
.. [4] NASA-SP125-COOLANT-SIDE. D. K. Huzel and D. H. Huang, NASA SP-125,
       2nd ed., 1971, printed p. 106, Eqs. (4-23) and (4-25),
       https://ntrs.nasa.gov/api/citations/19710019929/downloads/19710019929.pdf
.. [5] NASA-TM-1998-206313-WADEL. M. F. Wadel, NASA/TM-1998-206313, 1998,
       Eqs. (6)-(13), pp. 15-17, https://ntrs.nasa.gov/citations/19980017619
.. [6] NASA-SP-8087-1972. NASA SP-8087, Liquid Rocket Engine Fluid-Cooled
       Combustion Chambers, 1972, pp. 15-16 and 61,
       https://ntrs.nasa.gov/citations/19730022965
.. [7] NASA-TN-D-130-HATCH-PAPELL-1959. J. E. Hatch and S. S. Papell, NASA
       TN D-130, 1959, Eqs. (10)-(13), pp. 1-12,
       https://ntrs.nasa.gov/citations/19890068390
.. [8] NASA-TM-1998-206612-MEYER. M. L. Meyer, D. L. Linne and D. C. Rousar,
       NASA/TM-1998-206612, 1998, Eqs. (1)-(3) and Table III, pp. 5-8,
       https://ntrs.nasa.gov/citations/19980017166
.. [9] NASA-CR-134282-TOBIN-1974. R. D. Tobin, NASA CR-134282, 1974, Eq.
       (A-1) and Fig. A1, pp. A1-A3,
       https://ntrs.nasa.gov/citations/19740018134
"""

from __future__ import annotations

import itertools
from dataclasses import dataclass
from math import atan, exp, isfinite, log, log10, pi, sqrt, tanh
from typing import Literal, TypeVar

from rocket_propulsion.core.constants import STANDARD_GRAVITY_M_S2
from rocket_propulsion.core.errors import DomainError, InputError
from rocket_propulsion.core.metadata import WarningMessage
from rocket_propulsion.core.validation import numeric_range_guard, real
from rocket_propulsion.thermochemistry.transport import prandtl_number
from rocket_propulsion.thermodynamics.real_fluids import (
    FluidState,
    reference_fluid,
    reference_fluid_model,
)

__all__ = [
    "CoolantHeatTransfer",
    "CorrelationResult",
    "NearCriticalScreen",
    "burnout_margin_screen",
    "chen_friction_factor",
    "churchill_bernstein_cylinder_nusselt",
    "churchill_chu_vertical_plate_nusselt",
    "corrected_latent_heat",
    "curvature_heat_transfer_factor",
    "ethanol_burnout_heat_flux",
    "film_condensation_vertical_plate_nusselt",
    "film_cooled_wall_temperature",
    "filonenko_friction_factor",
    "flat_plate_laminar_nusselt",
    "flat_plate_turbulent_nusselt",
    "flat_plate_turbulent_power_law_nusselt",
    "gnielinski_nusselt",
    "gnielinski_power_law_nusselt",
    "graetz_number",
    "haaland_friction_factor",
    "hatch_papell_effectiveness",
    "hatch_papell_required_coolant_flow",
    "hatch_papell_velocity_ratio_function",
    "heat_transfer_coefficient",
    "hess_kunz_cl_factor",
    "hess_kunz_coefficient",
    "hess_kunz_nusselt",
    "hydraulic_diameter",
    "ito_curvature_friction_factor",
    "laminar_darcy_friction_factor",
    "laminar_entry_nusselt",
    "laminar_fully_developed_nusselt",
    "laminar_thermal_entry_length",
    "liquid_metal_nusselt",
    "mass_flux_reynolds_number",
    "mmh_burnout_heat_flux",
    "nusselt_number",
    "peak_pool_boiling_heat_flux",
    "peclet_number",
    "power_law_nusselt",
    "prandtl_number",
    "property_ratio_friction_factor",
    "property_ratio_nusselt_factor",
    "pseudocritical_temperature",
    "rayleigh_number",
    "rectangular_duct_hydraulic_diameter",
    "reynolds_number",
    "rohsenow_heat_flux",
    "rohsenow_wall_superheat",
    "rough_pipe_nusselt",
    "roughness_regime",
    "roughness_reynolds_number",
    "rte_coolant_nusselt",
    "sieder_tate_nusselt",
    "smooth_pipe_analogy_nusselt",
    "sp125_sieder_tate_nusselt",
    "sp125_supercritical_film_coefficient",
    "sp125_supercritical_film_heat_transfer",
    "stanton_number",
    "taylor_hydrogen_coefficient",
    "taylor_hydrogen_nusselt",
    "taylor_near_critical_region",
    "taylor_temperature_ratio_exponent",
    "taylor_wavelength",
    "turbulent_entry_nusselt_ratio",
]

#: Transition Reynolds number used by the source for circular and noncircular
#: ducts (Lienhard 6th ed., Example 7.7, p. 384; Fig. 7.6, p. 370).
PIPE_TRANSITION_REYNOLDS = 2300.0

_LIENHARD = "Lienhard & Lienhard, A Heat Transfer Textbook, 6th ed. (2024)"
_TAYLOR = "Taylor, NASA TN D-4332 (1968)"
_CR678 = "Aerojet-General, NASA CR-678 (1967)"
_SP125 = "Huzel & Huang, NASA SP-125 (1971)"
_WADEL = "Wadel, NASA/TM-1998-206313"
_SP8087 = "NASA SP-8087 (1972)"
_HATCH = "Hatch & Papell, NASA TN D-130 (1959)"
_MEYER = "Meyer, Linne & Rousar, NASA/TM-1998-206612"
_TOBIN = "Tobin, NASA CR-134282 (1974)"
#: One Btu/(in^2 s) in W/m^2 (IT Btu 1055.05585262 J over (0.0254 m)^2).
_BTU_IN2_S = 1055.05585262 / 0.0254 ** 2

# Exact unit conversions used to restate the source ranges in SI.
_RANKINE_TO_KELVIN = 5.0 / 9.0
_FOOT_M = 0.3048
_INCH_M = 0.0254
_PSI_PA = 6894.757293168361

_HYDROGEN_KEYS = ("parahydrogen", "normal_hydrogen")


@dataclass(frozen=True, slots=True)
class CorrelationResult:
    """Value of one correlation with its provenance and range warnings.

    Attributes
    ----------
    value : float
        Correlation output (Nusselt number, multiplier, heat flux, ...).
    quantity : str
        Name of the output quantity, with its unit when dimensional.
    correlation : str
        Short correlation name.
    source : str
        Source and equation locator.
    warnings : tuple of WarningMessage
        Extrapolation and applicability notes; empty inside the stated range.

    References
    ----------
    .. [1] LIENHARD-AHTT-6.00. J. H. Lienhard V and J. H. Lienhard IV, A
           Heat Transfer Textbook, 6th ed., 2024, https://ahtt.mit.edu/
    """

    value: float
    quantity: str
    correlation: str
    source: str
    warnings: tuple[WarningMessage, ...] = ()


@dataclass(frozen=True, slots=True)
class CoolantHeatTransfer:
    """Coolant-side coefficient evaluated with real-fluid density and heat capacity.

    Attributes
    ----------
    heat_transfer_coefficient_w_m2_k : float
        Local coefficient ``h = Nu k / D``.
    nusselt : float
        Nusselt number of the correlation, on its own reference state.
    reynolds, prandtl : float
        Reynolds and Prandtl numbers entering the correlation.
    reference_temperature_k : float
        Temperature at which the correlation's properties were evaluated.
    density_kg_m3, specific_heat_j_kg_k : float
        Equation-of-state density and isobaric heat capacity at the reference
        state.
    correlation, source : str
        Correlation name and equation locator.
    warnings : tuple of WarningMessage
        Correlation-range and equation-of-state warnings.

    References
    ----------
    .. [1] NASA-TN-D-4332-TAYLOR-1968. M. F. Taylor, NASA TN D-4332, 1968,
           https://ntrs.nasa.gov/citations/19680004771
    """

    heat_transfer_coefficient_w_m2_k: float
    nusselt: float
    reynolds: float
    prandtl: float
    reference_temperature_k: float
    density_kg_m3: float
    specific_heat_j_kg_k: float
    correlation: str
    source: str
    warnings: tuple[WarningMessage, ...] = ()


@dataclass(frozen=True, slots=True)
class NearCriticalScreen:
    """Result of Taylor's near-critical ("region 2") screen for hydrogen.

    Attributes
    ----------
    in_region : bool
        True when ``25 K < T_inlet < T*`` and ``p_c < p_inlet < 3.65 MPa``.
    pseudocritical_temperature_k : float or None
        ``T*``, the temperature of maximum ``c_p`` at the inlet pressure; None
        when the inlet pressure does not exceed the critical pressure.
    critical_pressure_pa : float
        Critical pressure of the equation of state used.
    warnings : tuple of WarningMessage
        A warning quoting Taylor's region-2 accuracy when ``in_region``.

    References
    ----------
    .. [1] NASA-TN-D-4332-TAYLOR-1968. M. F. Taylor, NASA TN D-4332, 1968,
           p. 14 and Fig. 3, https://ntrs.nasa.gov/citations/19680004771
    """

    in_region: bool
    pseudocritical_temperature_k: float | None
    critical_pressure_pa: float
    warnings: tuple[WarningMessage, ...] = ()


# ---------------------------------------------------------------------------
# Validation helpers
# ---------------------------------------------------------------------------


def _positive(value: object, label: str) -> float:
    """Return ``value`` as a finite positive float or raise :class:`DomainError`."""
    number = real(value, label)
    if number <= 0.0:
        raise DomainError(f"{label} must be positive.", field=label)
    return number


_Option = TypeVar("_Option", bound=str)


def _nonnegative(value: object, label: str) -> float:
    """Return ``value`` as a finite nonnegative float or raise :class:`DomainError`."""
    number = real(value, label)
    if number < 0.0:
        raise DomainError(f"{label} must be nonnegative.", field=label)
    return number


def _choice(value: object, options: tuple[_Option, ...], label: str) -> _Option:
    """Return ``value`` if it is one of ``options``; raise :class:`InputError` otherwise."""
    if not isinstance(value, str) or value not in options:
        raise InputError(f"{label} must be one of {', '.join(options)}.", field=label)
    return options[options.index(value)]


def _check_range(value: float, lo: float | None, hi: float | None, label: str,
                 correlation: str, extrapolate: bool, warnings: list[WarningMessage],
                 *, lo_open: bool = False, hi_open: bool = False) -> None:
    """Fail closed (or warn when ``extrapolate``) if ``value`` leaves ``[lo, hi]``.

    Open ends are used where the source states a strict inequality.
    """
    below = lo is not None and (value <= lo if lo_open else value < lo)
    above = hi is not None and (value >= hi if hi_open else value > hi)
    if not (below or above):
        return
    message = (f"{label} = {value:.6g} lies outside the stated range of {correlation} "
               f"[{lo if lo is not None else '-inf'}, {hi if hi is not None else 'inf'}].")
    if not extrapolate:
        raise DomainError(message, code="outside_correlation_range", field=label,
                          details={"value": value, "minimum": lo, "maximum": hi,
                                   "correlation": correlation})
    warnings.append(WarningMessage("correlation_extrapolated", message, "warning", label))


def _range_not_stated(correlation: str, detail: str) -> WarningMessage:
    """Info warning for a correlation whose source states no validity range."""
    return WarningMessage(
        "range_not_stated",
        f"{correlation}: the sources read state no numerical validity range ({detail}).",
        "info")


def _result(value: float, quantity: str, correlation: str, source: str,
            warnings: list[WarningMessage] | tuple[WarningMessage, ...] = ()) -> CorrelationResult:
    """Build a :class:`CorrelationResult` after a final finiteness check."""
    if not isfinite(value):
        raise DomainError(f"{correlation} produced a non-finite value.",
                          code="non_finite_result", details={"correlation": correlation})
    return CorrelationResult(value, quantity, correlation, source, tuple(warnings))


# ---------------------------------------------------------------------------
# Dimensionless groups and geometry
# ---------------------------------------------------------------------------


def reynolds_number(density_kg_m3: float, velocity_m_s: float, length_m: float,
                    viscosity_pa_s: float) -> float:
    """Reynolds number ``Re = rho u L / mu``.

    Parameters
    ----------
    density_kg_m3, velocity_m_s, length_m, viscosity_pa_s : float
        Density, (bulk) velocity, characteristic length and dynamic viscosity;
        all positive.

    Returns
    -------
    float
        Reynolds number.

    Raises
    ------
    DomainError
        If any input is not finite and positive.

    References
    ----------
    .. [1] LIENHARD-AHTT-6.00. Lienhard and Lienhard, A Heat Transfer
           Textbook, 6th ed., 2024, Sec. 7.1 (``Re_D = u_av D / nu``),
           https://ahtt.mit.edu/
    """
    rho = _positive(density_kg_m3, "density_kg_m3")
    u = _positive(velocity_m_s, "velocity_m_s")
    length = _positive(length_m, "length_m")
    mu = _positive(viscosity_pa_s, "viscosity_pa_s")
    return rho * u * length / mu


def mass_flux_reynolds_number(mass_flux_kg_m2_s: float, diameter_m: float,
                              viscosity_pa_s: float) -> float:
    """Bulk Reynolds number from the mass flux, ``Re_b = G D / mu_b``.

    Parameters
    ----------
    mass_flux_kg_m2_s : float
        Mass flow per unit cross-sectional area ``G``.
    diameter_m : float
        Inside (or hydraulic) diameter.
    viscosity_pa_s : float
        Dynamic viscosity at the reference temperature.

    Returns
    -------
    float
        Reynolds number.

    References
    ----------
    .. [1] NASA-TN-D-4332-TAYLOR-1968. M. F. Taylor, NASA TN D-4332, 1968,
           Symbols, p. 3 (``Re_b = G D / mu_b``),
           https://ntrs.nasa.gov/citations/19680004771
    """
    g = _positive(mass_flux_kg_m2_s, "mass_flux_kg_m2_s")
    d = _positive(diameter_m, "diameter_m")
    mu = _positive(viscosity_pa_s, "viscosity_pa_s")
    return g * d / mu


def peclet_number(reynolds: float, prandtl: float) -> float:
    """Peclet number ``Pe = Re Pr = u D / alpha``.

    Parameters
    ----------
    reynolds, prandtl : float
        Positive Reynolds and Prandtl numbers.

    Returns
    -------
    float
        Peclet number.

    References
    ----------
    .. [1] LIENHARD-AHTT-6.00. Lienhard and Lienhard, A Heat Transfer
           Textbook, 6th ed., 2024, Sec. 7.3, p. 376 (``Pe_D = u_av D / alpha``),
           https://ahtt.mit.edu/
    """
    return _positive(reynolds, "reynolds") * _positive(prandtl, "prandtl")


def graetz_number(reynolds: float, prandtl: float, diameter_m: float, distance_m: float) -> float:
    """Graetz number ``Gz = Re_D Pr D / x``, Lienhard Eq. (7.26).

    Parameters
    ----------
    reynolds, prandtl : float
        Reynolds and Prandtl numbers.
    diameter_m : float
        Pipe diameter ``D``.
    distance_m : float
        Distance ``x`` from the start of heating.

    Returns
    -------
    float
        Graetz number.

    References
    ----------
    .. [1] LIENHARD-AHTT-6.00. Lienhard and Lienhard, A Heat Transfer
           Textbook, 6th ed., 2024, Eq. (7.26), p. 362, https://ahtt.mit.edu/
    """
    return (_positive(reynolds, "reynolds") * _positive(prandtl, "prandtl")
            * _positive(diameter_m, "diameter_m") / _positive(distance_m, "distance_m"))


def rayleigh_number(expansion_coefficient_1_k: float, temperature_difference_k: float,
                    length_m: float, kinematic_viscosity_m2_s: float,
                    thermal_diffusivity_m2_s: float,
                    gravity_m_s2: float = STANDARD_GRAVITY_M_S2) -> float:
    """Rayleigh number ``Ra_L = g beta dT L^3 / (alpha nu)``.

    Parameters
    ----------
    expansion_coefficient_1_k : float
        Volumetric expansion coefficient ``beta`` (``1/T`` for an ideal gas).
    temperature_difference_k : float
        Magnitude of the wall-to-fluid temperature difference.
    length_m : float
        Characteristic length (plate height).
    kinematic_viscosity_m2_s, thermal_diffusivity_m2_s : float
        ``nu`` and ``alpha``.
    gravity_m_s2 : float, optional
        Gravitational acceleration; defaults to standard gravity.

    Returns
    -------
    float
        Rayleigh number.

    References
    ----------
    .. [1] LIENHARD-AHTT-6.00. Lienhard and Lienhard, A Heat Transfer
           Textbook, 6th ed., 2024, Example 8.1, p. 421, https://ahtt.mit.edu/
    """
    return (_positive(gravity_m_s2, "gravity_m_s2")
            * _positive(expansion_coefficient_1_k, "expansion_coefficient_1_k")
            * _positive(temperature_difference_k, "temperature_difference_k")
            * _positive(length_m, "length_m") ** 3
            / (_positive(kinematic_viscosity_m2_s, "kinematic_viscosity_m2_s")
               * _positive(thermal_diffusivity_m2_s, "thermal_diffusivity_m2_s")))


def nusselt_number(heat_transfer_coefficient_w_m2_k: float, length_m: float,
                   conductivity_w_m_k: float) -> float:
    """Nusselt number ``Nu = h L / k``.

    Parameters
    ----------
    heat_transfer_coefficient_w_m2_k, length_m, conductivity_w_m_k : float
        Coefficient, characteristic length and fluid conductivity.

    Returns
    -------
    float
        Nusselt number.

    References
    ----------
    .. [1] NASA-TN-D-4332-TAYLOR-1968. M. F. Taylor, NASA TN D-4332, 1968,
           Symbols, p. 3 (``Nu = h D / k``),
           https://ntrs.nasa.gov/citations/19680004771
    """
    return (_positive(heat_transfer_coefficient_w_m2_k, "heat_transfer_coefficient_w_m2_k")
            * _positive(length_m, "length_m") / _positive(conductivity_w_m_k, "conductivity_w_m_k"))


def heat_transfer_coefficient(nusselt: float, conductivity_w_m_k: float, length_m: float) -> float:
    """Heat-transfer coefficient ``h = Nu k / L``.

    Parameters
    ----------
    nusselt, conductivity_w_m_k, length_m : float
        Nusselt number, fluid conductivity and characteristic length.

    Returns
    -------
    float
        Coefficient in W/(m^2 K).

    References
    ----------
    .. [1] LIENHARD-AHTT-6.00. Lienhard and Lienhard, A Heat Transfer
           Textbook, 6th ed., 2024, Example 7.3, p. 372, https://ahtt.mit.edu/
    """
    return (_positive(nusselt, "nusselt") * _positive(conductivity_w_m_k, "conductivity_w_m_k")
            / _positive(length_m, "length_m"))


def stanton_number(nusselt: float, reynolds: float, prandtl: float) -> float:
    """Stanton number ``St = Nu / (Re Pr) = h / (rho c_p u)``.

    Parameters
    ----------
    nusselt, reynolds, prandtl : float
        Positive dimensionless groups.

    Returns
    -------
    float
        Stanton number.

    References
    ----------
    .. [1] LIENHARD-AHTT-6.00. Lienhard and Lienhard, A Heat Transfer
           Textbook, 6th ed., 2024, Eq. (6.111), p. 325, https://ahtt.mit.edu/
    """
    return _positive(nusselt, "nusselt") / (_positive(reynolds, "reynolds")
                                            * _positive(prandtl, "prandtl"))


def hydraulic_diameter(flow_area_m2: float, wetted_perimeter_m: float) -> float:
    """Hydraulic diameter ``D_h = 4 A_c / P``, Lienhard Eq. (7.60).

    Parameters
    ----------
    flow_area_m2 : float
        Cross-sectional flow area.
    wetted_perimeter_m : float
        Whole wetted perimeter (also when only part of it is heated).

    Returns
    -------
    float
        Hydraulic diameter, m.

    Notes
    -----
    For turbulent flow the circular-tube correlations with ``D -> D_h`` are
    typically within about +/-20 % (often +/-10 %), worse for sharp-cornered
    sections; laminar values depend on the shape (see
    :func:`laminar_fully_developed_nusselt`) [1]_.

    References
    ----------
    .. [1] LIENHARD-AHTT-6.00. Lienhard and Lienhard, A Heat Transfer
           Textbook, 6th ed., 2024, Eqs. (7.60)-(7.61), pp. 382-385,
           https://ahtt.mit.edu/
    """
    return 4.0 * _positive(flow_area_m2, "flow_area_m2") / _positive(wetted_perimeter_m,
                                                                     "wetted_perimeter_m")


def rectangular_duct_hydraulic_diameter(width_m: float, height_m: float) -> float:
    """Hydraulic diameter of a rectangular duct, ``2 a b / (a + b)``, Eq. (7.61a).

    Parameters
    ----------
    width_m, height_m : float
        Side lengths ``a`` and ``b``.

    Returns
    -------
    float
        Hydraulic diameter, m.

    References
    ----------
    .. [1] LIENHARD-AHTT-6.00. Lienhard and Lienhard, A Heat Transfer
           Textbook, 6th ed., 2024, Eq. (7.61a), p. 383, https://ahtt.mit.edu/
    """
    a = _positive(width_m, "width_m")
    b = _positive(height_m, "height_m")
    return 2.0 * a * b / (a + b)


def roughness_reynolds_number(reynolds: float, relative_roughness: float,
                              darcy_friction_factor: float) -> float:
    """Roughness Reynolds number ``Re_eps = Re_D (eps/D) sqrt(f/8)``, Eq. (7.48a).

    Parameters
    ----------
    reynolds : float
        Pipe Reynolds number.
    relative_roughness : float
        ``eps / D``.
    darcy_friction_factor : float
        Darcy friction factor ``f``.

    Returns
    -------
    float
        ``Re_eps``.

    References
    ----------
    .. [1] LIENHARD-AHTT-6.00. Lienhard and Lienhard, A Heat Transfer
           Textbook, 6th ed., 2024, Eq. (7.48a), p. 373, https://ahtt.mit.edu/
    """
    return (_positive(reynolds, "reynolds") * _nonnegative(relative_roughness, "relative_roughness")
            * sqrt(_positive(darcy_friction_factor, "darcy_friction_factor") / 8.0))


def roughness_regime(roughness_reynolds: float) -> str:
    """Classify a wall by ``Re_eps``, Lienhard Eqs. (7.48b)-(7.48d).

    Parameters
    ----------
    roughness_reynolds : float
        Roughness Reynolds number.

    Returns
    -------
    {"hydraulically_smooth", "transitionally_rough", "fully_rough"}
        ``Re_eps < 5``, ``5 <= Re_eps <= 70`` and ``Re_eps > 70``.

    References
    ----------
    .. [1] LIENHARD-AHTT-6.00. Lienhard and Lienhard, A Heat Transfer
           Textbook, 6th ed., 2024, Eqs. (7.48b)-(7.48d), p. 373,
           https://ahtt.mit.edu/
    """
    value = _nonnegative(roughness_reynolds, "roughness_reynolds")
    if value < 5.0:
        return "hydraulically_smooth"
    if value <= 70.0:
        return "transitionally_rough"
    return "fully_rough"


# ---------------------------------------------------------------------------
# Friction factors (Darcy)
# ---------------------------------------------------------------------------


def laminar_darcy_friction_factor(reynolds: float, *, extrapolate: bool = False) -> CorrelationResult:
    """Fully developed laminar circular-pipe friction factor ``f = 64 / Re_D``.

    Parameters
    ----------
    reynolds : float
        Reynolds number; the source places transition at about 2300.
    extrapolate : bool, optional
        Return a warned value above 2300 instead of raising.

    Returns
    -------
    CorrelationResult
        Darcy friction factor.

    References
    ----------
    .. [1] LIENHARD-AHTT-6.00. Lienhard and Lienhard, A Heat Transfer
           Textbook, 6th ed., 2024, Fig. 7.6, p. 370, and Example 7.7, p. 384,
           https://ahtt.mit.edu/
    """
    re = _positive(reynolds, "reynolds")
    warnings: list[WarningMessage] = []
    _check_range(re, None, PIPE_TRANSITION_REYNOLDS, "reynolds", "laminar f = 64/Re",
                 extrapolate, warnings)
    return _result(64.0 / re, "darcy_friction_factor", "laminar f = 64/Re",
                   f"{_LIENHARD}, Fig. 7.6", warnings)


def filonenko_friction_factor(reynolds: float, *, extrapolate: bool = False) -> CorrelationResult:
    """Smooth-pipe friction factor of Filonenko, Lienhard Eq. (7.42).

    ``f = (1.82 log10 Re_D - 1.64)^-2``.

    Parameters
    ----------
    reynolds : float
        Reynolds number. The source gives no separate range; the range of
        the Gnielinski equation it serves, ``2300 <= Re_D <= 5e6``, is enforced.
    extrapolate : bool, optional
        Return a warned value outside that range instead of raising.

    Returns
    -------
    CorrelationResult
        Darcy friction factor.

    References
    ----------
    .. [1] LIENHARD-AHTT-6.00. Lienhard and Lienhard, A Heat Transfer
           Textbook, 6th ed., 2024, Eq. (7.42), p. 369, https://ahtt.mit.edu/
    """
    re = _positive(reynolds, "reynolds")
    warnings: list[WarningMessage] = []
    _check_range(re, PIPE_TRANSITION_REYNOLDS, 5.0e6, "reynolds", "Filonenko (with Gnielinski)",
                 extrapolate, warnings)
    return _result((1.82 * log10(re) - 1.64) ** -2, "darcy_friction_factor", "Filonenko",
                   f"{_LIENHARD}, Eq. (7.42)", warnings)


def haaland_friction_factor(reynolds: float, relative_roughness: float, *,
                            extrapolate: bool = False) -> CorrelationResult:
    """Haaland's explicit rough-pipe friction factor, Lienhard Eq. (7.50).

    ``f = {1.8 log10[6.9/Re_D + ((eps/D)/3.7)^1.11]}^-2``.

    Parameters
    ----------
    reynolds : float
        Reynolds number, ``4000 <= Re_D <= 1e8``.
    relative_roughness : float
        ``eps / D``, ``0 <= eps/D <= 0.05``.
    extrapolate : bool, optional
        Return a warned value outside the range instead of raising.

    Returns
    -------
    CorrelationResult
        Darcy friction factor.

    References
    ----------
    .. [1] LIENHARD-AHTT-6.00. Lienhard and Lienhard, A Heat Transfer
           Textbook, 6th ed., 2024, Eq. (7.50), p. 373, https://ahtt.mit.edu/
    """
    re = _positive(reynolds, "reynolds")
    rel = _nonnegative(relative_roughness, "relative_roughness")
    warnings: list[WarningMessage] = []
    _check_range(re, 4000.0, 1.0e8, "reynolds", "Haaland", extrapolate, warnings)
    _check_range(rel, 0.0, 0.05, "relative_roughness", "Haaland", extrapolate, warnings)
    value = (1.8 * log10(6.9 / re + (rel / 3.7) ** 1.11)) ** -2
    return _result(value, "darcy_friction_factor", "Haaland", f"{_LIENHARD}, Eq. (7.50)",
                   warnings)


@numeric_range_guard("reynolds")
def chen_friction_factor(reynolds: float, relative_roughness: float) -> CorrelationResult:
    """Explicit Colebrook form (Chen) used by NASA's RTE code, Wadel Eq. (12).

    ``1/sqrt(f) = -2.0 log10[e/(3.7065 d) - (5.0452/Re) log10((e/d)^1.1098/2.8257
    + 5.8506/Re^0.8981)]``.

    Parameters
    ----------
    reynolds : float
        Coolant Reynolds number (turbulent per the source).
    relative_roughness : float
        ``e / d``, nonnegative.

    Returns
    -------
    CorrelationResult
        Darcy friction factor, with a ``range_not_stated`` info warning.

    Raises
    ------
    DomainError
        If the logarithm arguments become nonpositive (very low ``Re``).

    References
    ----------
    .. [1] NASA-TM-1998-206313-WADEL. M. F. Wadel, NASA/TM-1998-206313,
           Eq. (12), p. 17, https://ntrs.nasa.gov/citations/19980017619
    """
    re = _positive(reynolds, "reynolds")
    rel = _nonnegative(relative_roughness, "relative_roughness")
    inner = rel ** 1.1098 / 2.8257 + 5.8506 / re ** 0.8981
    outer = rel / 3.7065 - 5.0452 / re * log10(inner)
    if outer <= 0.0:
        raise DomainError("Chen friction-factor argument is not positive.",
                          code="outside_correlation_range", field="reynolds")
    value = (-2.0 * log10(outer)) ** -2
    return _result(value, "darcy_friction_factor", "Chen explicit Colebrook (RTE)",
                   f"{_WADEL}, Eq. (12)",
                   [_range_not_stated("Chen (RTE)", "turbulent coolant flow is assumed")])


def ito_curvature_friction_factor(reynolds: float, tube_radius_m: float,
                                  radius_of_curvature_m: float) -> CorrelationResult:
    """Ito's curved-pipe friction ratio ``f_c/f_s = [Re (r/R)^2]^0.05``.

    Parameters
    ----------
    reynolds : float
        Reynolds number.
    tube_radius_m : float
        Tube radius ``r`` (hydraulic radius ``d_h / 2`` in Wadel's use).
    radius_of_curvature_m : float
        Radius of curvature ``R`` of the passage centreline.

    Returns
    -------
    CorrelationResult
        Ratio of curved- to straight-tube friction, with a
        ``range_not_stated`` info warning.

    References
    ----------
    .. [1] NASA-CR-678-AEROJET-1967. Aerojet-General Corp., NASA CR-678,
           1967, Eq. (4), p. 4, https://ntrs.nasa.gov/citations/19670008827
    .. [2] NASA-TM-1998-206313-WADEL. M. F. Wadel, NASA/TM-1998-206313,
           Eq. (13), p. 17, https://ntrs.nasa.gov/citations/19980017619
    """
    re = _positive(reynolds, "reynolds")
    ratio = _positive(tube_radius_m, "tube_radius_m") / _positive(radius_of_curvature_m,
                                                                  "radius_of_curvature_m")
    return _result((re * ratio ** 2) ** 0.05, "friction_ratio", "Ito curved-pipe friction",
                   f"{_CR678}, Eq. (4); {_WADEL}, Eq. (13)",
                   [_range_not_stated("Ito friction ratio", "Re and r/R limits are not given")])


# ---------------------------------------------------------------------------
# Internal laminar flow
# ---------------------------------------------------------------------------

_LAMINAR_TABLE = {
    "circular": (3.657, 4.364),
    "square": (2.976, 3.608),
    "rectangle_2to1": (3.391, 4.123),
    "rectangle_4to1": (4.439, 5.331),
    "rectangle_8to1": (5.597, 6.490),
    "parallel_plates": (7.541, 8.235),
}

WallCondition = Literal["uniform_temperature", "uniform_heat_flux"]


def laminar_fully_developed_nusselt(cross_section: str, wall: WallCondition) -> CorrelationResult:
    """Fully developed laminar Nusselt number on the hydraulic diameter.

    Parameters
    ----------
    cross_section : str
        ``"circular"``, ``"square"``, ``"rectangle_2to1"``, ``"rectangle_4to1"``,
        ``"rectangle_8to1"``, ``"parallel_plates"`` (Table 7.4) or
        ``"parallel_plates_one_side_insulated"`` (Eq. 7.63, uniform flux on
        one plate, the other adiabatic; only with ``wall="uniform_heat_flux"``).
    wall : {"uniform_temperature", "uniform_heat_flux"}
        Thermal boundary condition.

    Returns
    -------
    CorrelationResult
        ``Nu_Dh``.

    Notes
    -----
    The rectangular uniform-flux values assume a peripherally uniform wall
    temperature (no wall conduction resistance around the perimeter) [1]_.

    References
    ----------
    .. [1] LIENHARD-AHTT-6.00. Lienhard and Lienhard, A Heat Transfer
           Textbook, 6th ed., 2024, Eqs. (7.22)-(7.23), (7.63) and Table 7.4,
           pp. 359-361, 386-387, https://ahtt.mit.edu/
    """
    wall = _choice(wall, ("uniform_temperature", "uniform_heat_flux"), "wall")
    options = (*_LAMINAR_TABLE, "parallel_plates_one_side_insulated")
    section = _choice(cross_section, options, "cross_section")
    if section == "parallel_plates_one_side_insulated":
        if wall != "uniform_heat_flux":
            raise InputError("The one-side-insulated plate value is for uniform heat flux.",
                             field="wall")
        return _result(5.385, "nusselt", "laminar parallel plates, one side adiabatic",
                       f"{_LIENHARD}, Eq. (7.63)")
    value = _LAMINAR_TABLE[section][0 if wall == "uniform_temperature" else 1]
    return _result(value, "nusselt", f"laminar fully developed ({section})",
                   f"{_LIENHARD}, Table 7.4")


def laminar_thermal_entry_length(reynolds: float, prandtl: float, diameter_m: float,
                                 wall: WallCondition) -> CorrelationResult:
    """Laminar thermal entry length, Lienhard Eq. (7.24).

    ``x_et / D ~= 0.034 Re_D Pr`` (uniform wall temperature) or
    ``0.043 Re_D Pr`` (uniform heat flux), for a hydrodynamically developed
    flow; beyond ``x_et`` the Nusselt number is within 5 % of the fully
    developed value.

    Parameters
    ----------
    reynolds, prandtl : float
        Reynolds (laminar, ``<= 2300``) and Prandtl numbers.
    diameter_m : float
        Pipe diameter.
    wall : {"uniform_temperature", "uniform_heat_flux"}
        Wall condition.

    Returns
    -------
    CorrelationResult
        Entry length in m.

    References
    ----------
    .. [1] LIENHARD-AHTT-6.00. Lienhard and Lienhard, A Heat Transfer
           Textbook, 6th ed., 2024, Eq. (7.24) and footnote 4, p. 361,
           https://ahtt.mit.edu/
    """
    wall = _choice(wall, ("uniform_temperature", "uniform_heat_flux"), "wall")
    re = _positive(reynolds, "reynolds")
    warnings: list[WarningMessage] = []
    _check_range(re, None, PIPE_TRANSITION_REYNOLDS, "reynolds", "laminar entry length",
                 False, warnings)
    coefficient = 0.034 if wall == "uniform_temperature" else 0.043
    value = coefficient * re * _positive(prandtl, "prandtl") * _positive(diameter_m, "diameter_m")
    return _result(value, "length_m", "laminar thermal entry length", f"{_LIENHARD}, Eq. (7.24)")


def laminar_entry_nusselt(graetz: float, wall: WallCondition, *,
                          average: bool = False) -> CorrelationResult:
    """Thermal-entry Nusselt number for hydrodynamically developed laminar pipe flow.

    Uniform wall temperature: local Eq. (7.28), average Eq. (7.29); uniform
    heat flux: local Eq. (7.30). The fits reproduce the exact solutions to
    within 1 %.

    Parameters
    ----------
    graetz : float
        ``Gz = Re_D Pr D / x`` (:func:`graetz_number`).
    wall : {"uniform_temperature", "uniform_heat_flux"}
        Wall condition.
    average : bool, optional
        Average Nusselt number over ``0..x`` (uniform temperature only).

    Returns
    -------
    CorrelationResult
        ``Nu_D`` (local or average).

    Raises
    ------
    InputError
        If an average is requested for the uniform-flux case (the source gives
        none).

    References
    ----------
    .. [1] LIENHARD-AHTT-6.00. Lienhard and Lienhard, A Heat Transfer
           Textbook, 6th ed., 2024, Eqs. (7.28)-(7.30), pp. 362-363,
           https://ahtt.mit.edu/
    """
    wall = _choice(wall, ("uniform_temperature", "uniform_heat_flux"), "wall")
    gz = _positive(graetz, "graetz")
    if wall == "uniform_temperature":
        if average:
            value = 3.657 / tanh(2.264 * gz ** (-1.0 / 3.0) + 1.7 * gz ** (-2.0 / 3.0)) \
                + 0.0499 * gz * tanh(1.0 / gz)
            return _result(value, "nusselt", "laminar entry, isothermal wall, average",
                           f"{_LIENHARD}, Eq. (7.29)")
        if gz <= 1000.0:
            value = 3.657 + 0.2362 * gz ** 0.488 * exp(-57.2 / gz)
        else:
            value = 1.077 * gz ** (1.0 / 3.0) - 0.7
        return _result(value, "nusselt", "laminar entry, isothermal wall, local",
                       f"{_LIENHARD}, Eq. (7.28)")
    if average:
        raise InputError("The source gives no average Nusselt fit for uniform heat flux.",
                         field="average")
    if gz <= 667.0:
        value = 4.364 + 0.263 * gz ** 0.506 * exp(-41.0 / gz)
    elif gz <= 2.0e4:
        value = 1.302 * gz ** (1.0 / 3.0) - 0.5
    else:
        value = 1.302 * gz ** (1.0 / 3.0) - 1.0
    return _result(value, "nusselt", "laminar entry, uniform heat flux, local",
                   f"{_LIENHARD}, Eq. (7.30)")


# ---------------------------------------------------------------------------
# Internal turbulent flow
# ---------------------------------------------------------------------------


def _darcy_or_filonenko(reynolds: float, darcy_friction_factor: float | None,
                        extrapolate: bool, warnings: list[WarningMessage]) -> float:
    """Return the supplied friction factor or Filonenko's smooth-pipe value."""
    if darcy_friction_factor is not None:
        return _positive(darcy_friction_factor, "darcy_friction_factor")
    friction = filonenko_friction_factor(reynolds, extrapolate=extrapolate)
    warnings.extend(friction.warnings)
    return friction.value


def gnielinski_nusselt(reynolds: float, prandtl: float, *,
                       darcy_friction_factor: float | None = None,
                       extrapolate: bool = False) -> CorrelationResult:
    """Gnielinski's turbulent/transitional pipe correlation, Lienhard Eq. (7.41).

    ``Nu_D = (f/8)(Re_D - 1000) Pr / [1 + 12.7 sqrt(f/8)(Pr^(2/3) - 1)]``,
    for ``2300 <= Re_D <= 5e6`` and ``0.6 <= Pr <= 1e5``; uniform wall
    temperature or heat flux. More than 90 % of the liquid data examined lie
    within +/-20 %, gases closer.

    Parameters
    ----------
    reynolds, prandtl : float
        Bulk Reynolds and Prandtl numbers.
    darcy_friction_factor : float, optional
        Darcy ``f``; defaults to Filonenko's smooth-pipe Eq. (7.42).
    extrapolate : bool, optional
        Return a warned value outside the range instead of raising.

    Returns
    -------
    CorrelationResult
        ``Nu_D`` with properties at the bulk temperature.

    References
    ----------
    .. [1] LIENHARD-AHTT-6.00. Lienhard and Lienhard, A Heat Transfer
           Textbook, 6th ed., 2024, Eqs. (7.41)-(7.42), p. 369,
           https://ahtt.mit.edu/
    """
    re = _positive(reynolds, "reynolds")
    pr = _positive(prandtl, "prandtl")
    warnings: list[WarningMessage] = []
    _check_range(re, PIPE_TRANSITION_REYNOLDS, 5.0e6, "reynolds", "Gnielinski", extrapolate,
                 warnings)
    _check_range(pr, 0.6, 1.0e5, "prandtl", "Gnielinski", extrapolate, warnings)
    f = _darcy_or_filonenko(re, darcy_friction_factor, True, warnings)
    value = (f / 8.0) * (re - 1000.0) * pr / (1.0 + 12.7 * sqrt(f / 8.0) * (pr ** (2.0 / 3.0) - 1.0))
    return _result(value, "nusselt", "Gnielinski", f"{_LIENHARD}, Eq. (7.41)",
                   _dedupe(warnings))


def smooth_pipe_analogy_nusselt(reynolds: float, prandtl: float, *,
                                darcy_friction_factor: float | None = None,
                                extrapolate: bool = False) -> CorrelationResult:
    """Petukhov-type Reynolds-Colburn analogy for smooth pipes, Lienhard Eq. (7.35).

    ``Nu_D = (f/8) Re_D Pr / [1 + 12.7 sqrt(f/8)(Pr^(2/3) - 1)]`` for
    ``Pr >= 0.7`` (Eq. 7.32 condition), fully developed turbulent flow.

    Parameters
    ----------
    reynolds, prandtl : float
        Bulk Reynolds (above the 2300 transition) and Prandtl numbers.
    darcy_friction_factor : float, optional
        Darcy ``f``; defaults to Filonenko's Eq. (7.42).
    extrapolate : bool, optional
        Return a warned value outside the range instead of raising.

    Returns
    -------
    CorrelationResult
        ``Nu_D``.

    References
    ----------
    .. [1] LIENHARD-AHTT-6.00. Lienhard and Lienhard, A Heat Transfer
           Textbook, 6th ed., 2024, Eqs. (7.32)-(7.35), pp. 366-367,
           https://ahtt.mit.edu/
    """
    re = _positive(reynolds, "reynolds")
    pr = _positive(prandtl, "prandtl")
    warnings: list[WarningMessage] = []
    _check_range(re, PIPE_TRANSITION_REYNOLDS, None, "reynolds", "smooth-pipe analogy (7.35)",
                 extrapolate, warnings, lo_open=True)
    _check_range(pr, 0.7, None, "prandtl", "smooth-pipe analogy (7.35)", extrapolate, warnings)
    f = _darcy_or_filonenko(re, darcy_friction_factor, True, warnings)
    value = (f / 8.0) * re * pr / (1.0 + 12.7 * sqrt(f / 8.0) * (pr ** (2.0 / 3.0) - 1.0))
    return _result(value, "nusselt", "smooth-pipe Reynolds-Colburn analogy",
                   f"{_LIENHARD}, Eq. (7.35)", _dedupe(warnings))


_POWER_LAWS = {
    "colburn": (0.023, 1.0 / 3.0, "Eq. (7.39a)"),
    "dittus_boelter": (0.0243, 0.4, "Eq. (7.39b)"),
    "mcadams": (0.0225, 0.4, "Eq. (7.39c)"),
}


def power_law_nusselt(reynolds: float, prandtl: float, form: str = "colburn", *,
                      extrapolate: bool = False) -> CorrelationResult:
    """Historical smooth-pipe power laws ``Nu_D = C Re_D^0.8 Pr^n``.

    ``form="colburn"``: ``C = 0.023, n = 1/3`` (Eq. 7.39a);
    ``"dittus_boelter"``: ``0.0243, 0.4`` (Eq. 7.39b, heating, as given by
    Lienhard); ``"mcadams"``: ``0.0225, 0.4`` (Eq. 7.39c). For ``Re_D >= 1e4``
    (Eq. 7.38 condition) and ``0.67 <= Pr <= 100``, small property variation;
    maximum errors +25 %/-40 %.

    Parameters
    ----------
    reynolds, prandtl : float
        Bulk Reynolds and Prandtl numbers.
    form : {"colburn", "dittus_boelter", "mcadams"}, optional
        Which power law.
    extrapolate : bool, optional
        Return a warned value outside the range instead of raising.

    Returns
    -------
    CorrelationResult
        ``Nu_D``.

    Notes
    -----
    The widely quoted ``0.023 Re^0.8 Pr^0.4`` is not one of the three printed
    forms; it is the leading factor of Taylor's correlation
    (:func:`taylor_hydrogen_nusselt`) and the RTE coolant default
    (:func:`rte_coolant_nusselt`).

    References
    ----------
    .. [1] LIENHARD-AHTT-6.00. Lienhard and Lienhard, A Heat Transfer
           Textbook, 6th ed., 2024, Eqs. (7.38)-(7.39c), pp. 368-369,
           https://ahtt.mit.edu/
    """
    form = _choice(form, tuple(_POWER_LAWS), "form")
    re = _positive(reynolds, "reynolds")
    pr = _positive(prandtl, "prandtl")
    warnings: list[WarningMessage] = []
    _check_range(re, 1.0e4, None, "reynolds", f"{form} power law", extrapolate, warnings)
    _check_range(pr, 0.67, 100.0, "prandtl", f"{form} power law", extrapolate, warnings)
    c, n, eq = _POWER_LAWS[form]
    return _result(c * re ** 0.8 * pr ** n, "nusselt", f"{form} power law", f"{_LIENHARD}, {eq}",
                   warnings)


def sieder_tate_nusselt(reynolds: float, prandtl: float, bulk_to_wall_viscosity_ratio: float, *,
                        extrapolate: bool = False) -> CorrelationResult:
    """Sieder-Tate liquid correlation, Lienhard Eq. (7.40).

    ``Nu_D = 0.023 Re_D^0.8 Pr^(1/3) (mu_b/mu_w)^0.14``, properties at the
    bulk temperature except ``mu_w``; same range as the power laws
    (``Re_D >= 1e4``, ``0.67 <= Pr <= 100``).

    Parameters
    ----------
    reynolds, prandtl : float
        Bulk Reynolds and Prandtl numbers.
    bulk_to_wall_viscosity_ratio : float
        ``mu_b / mu_w``.
    extrapolate : bool, optional
        Return a warned value outside the range instead of raising.

    Returns
    -------
    CorrelationResult
        ``Nu_D``.

    References
    ----------
    .. [1] LIENHARD-AHTT-6.00. Lienhard and Lienhard, A Heat Transfer
           Textbook, 6th ed., 2024, Eq. (7.40), pp. 368-369,
           https://ahtt.mit.edu/
    """
    base = power_law_nusselt(reynolds, prandtl, "colburn", extrapolate=extrapolate)
    ratio = _positive(bulk_to_wall_viscosity_ratio, "bulk_to_wall_viscosity_ratio")
    return _result(base.value * ratio ** 0.14, "nusselt", "Sieder-Tate", f"{_LIENHARD}, Eq. (7.40)",
                   base.warnings)


def gnielinski_power_law_nusselt(reynolds: float, prandtl: float, *,
                                 extrapolate: bool = False) -> CorrelationResult:
    """Gnielinski's smooth-pipe power-law fits, Lienhard Eqs. (7.43a)-(7.43b).

    ``0.0214 (Re^0.8 - 100) Pr^0.4`` for ``0.6 <= Pr <= 1.5`` and
    ``0.012 (Re^0.87 - 280) Pr^0.4`` for ``1.5 <= Pr <= 500``;
    ``2300 <= Re_D <= 5e6``.

    Parameters
    ----------
    reynolds, prandtl : float
        Bulk Reynolds and Prandtl numbers.
    extrapolate : bool, optional
        Return a warned value outside the range instead of raising.

    Returns
    -------
    CorrelationResult
        ``Nu_D``.

    References
    ----------
    .. [1] LIENHARD-AHTT-6.00. Lienhard and Lienhard, A Heat Transfer
           Textbook, 6th ed., 2024, Eqs. (7.43a)-(7.43b), p. 371,
           https://ahtt.mit.edu/
    """
    re = _positive(reynolds, "reynolds")
    pr = _positive(prandtl, "prandtl")
    warnings: list[WarningMessage] = []
    _check_range(re, PIPE_TRANSITION_REYNOLDS, 5.0e6, "reynolds", "Gnielinski power law",
                 extrapolate, warnings)
    _check_range(pr, 0.6, 500.0, "prandtl", "Gnielinski power law", extrapolate, warnings)
    if pr <= 1.5:
        return _result(0.0214 * (re ** 0.8 - 100.0) * pr ** 0.4, "nusselt",
                       "Gnielinski power law (7.43a)", f"{_LIENHARD}, Eq. (7.43a)", warnings)
    return _result(0.012 * (re ** 0.87 - 280.0) * pr ** 0.4, "nusselt",
                   "Gnielinski power law (7.43b)", f"{_LIENHARD}, Eq. (7.43b)", warnings)


def property_ratio_nusselt_factor(fluid: Literal["liquid", "gas"], ratio: float, *,
                                  heating: bool, extrapolate: bool = False) -> CorrelationResult:
    """Variable-property multiplier on a bulk-property Nusselt number.

    Liquids, Eq. (7.44): ``(mu_b/mu_w)^n``, ``n = 0.11`` (heating, ``T_w > T_b``)
    or ``0.25`` (cooling), for ``0.025 <= mu_b/mu_w <= 12.5``.
    Gases, Eq. (7.45): ``(T_b/T_w)^n``, ``n = 0.47`` (heating) or ``0``
    (cooling), for ``0.27 <= T_b/T_w <= 2.7`` (kelvin). The source applies
    these to Gnielinski's equation for smooth pipes.

    Parameters
    ----------
    fluid : {"liquid", "gas"}
        Fluid class.
    ratio : float
        ``mu_b / mu_w`` (liquid) or ``T_b / T_w`` (gas).
    heating : bool
        True when the wall is hotter than the bulk fluid.
    extrapolate : bool, optional
        Return a warned value outside the range instead of raising.

    Returns
    -------
    CorrelationResult
        Multiplier.

    References
    ----------
    .. [1] LIENHARD-AHTT-6.00. Lienhard and Lienhard, A Heat Transfer
           Textbook, 6th ed., 2024, Eqs. (7.44)-(7.45), p. 371,
           https://ahtt.mit.edu/
    """
    fluid = _choice(fluid, ("liquid", "gas"), "fluid")
    value = _positive(ratio, "ratio")
    warnings: list[WarningMessage] = []
    if fluid == "liquid":
        _check_range(value, 0.025, 12.5, "ratio", "liquid viscosity-ratio correction",
                     extrapolate, warnings)
        n, eq = (0.11 if heating else 0.25), "Eq. (7.44)"
    else:
        _check_range(value, 0.27, 2.7, "ratio", "gas temperature-ratio correction",
                     extrapolate, warnings)
        n, eq = (0.47 if heating else 0.0), "Eq. (7.45)"
    return _result(value ** n, "multiplier", f"{fluid} property-ratio correction",
                   f"{_LIENHARD}, {eq}", warnings)


def property_ratio_friction_factor(fluid: Literal["liquid", "gas"], ratio: float, *,
                                   heating: bool, extrapolate: bool = False) -> CorrelationResult:
    """Variable-property multiplier on the bulk friction factor.

    Liquids, Eq. (7.46): ``K = (7 - mu_b/mu_w)/6`` (heating) or
    ``(mu_b/mu_w)^-0.24`` (cooling), ``0.5 <= mu_b/mu_w <= 3``.
    Gases, Eq. (7.47): ``(T_b/T_w)^m`` with ``m ~= 0.23`` for both directions,
    ``0.14 <= T_b/T_w <= 3.3``.

    Parameters
    ----------
    fluid : {"liquid", "gas"}
        Fluid class.
    ratio : float
        ``mu_b / mu_w`` (liquid) or ``T_b / T_w`` (gas).
    heating : bool
        True when the wall is hotter than the bulk fluid.
    extrapolate : bool, optional
        Return a warned value outside the range instead of raising.

    Returns
    -------
    CorrelationResult
        Multiplier on ``f``.

    References
    ----------
    .. [1] LIENHARD-AHTT-6.00. Lienhard and Lienhard, A Heat Transfer
           Textbook, 6th ed., 2024, Eqs. (7.46)-(7.47), p. 371,
           https://ahtt.mit.edu/
    """
    fluid = _choice(fluid, ("liquid", "gas"), "fluid")
    value = _positive(ratio, "ratio")
    warnings: list[WarningMessage] = []
    if fluid == "liquid":
        _check_range(value, 0.5, 3.0, "ratio", "liquid friction correction", extrapolate,
                     warnings)
        factor = (7.0 - value) / 6.0 if heating else value ** -0.24
        eq = "Eq. (7.46)"
    else:
        _check_range(value, 0.14, 3.3, "ratio", "gas friction correction", extrapolate, warnings)
        factor, eq = value ** 0.23, "Eq. (7.47)"
    return _result(factor, "multiplier", f"{fluid} friction property correction",
                   f"{_LIENHARD}, {eq}", warnings)


_ENTRY_INLETS = {
    "long_straight_pipe": (0.9756, 0.760),
    "square_edged_inlet": (2.4254, 0.676),
    "bend_180_deg": (0.9759, 0.700),
    "bend_90_deg": (1.0517, 0.629),
    "sharp_elbow_90_deg": (2.0152, 0.614),
}


def turbulent_entry_nusselt_ratio(length_to_diameter: float, inlet: str, *,
                                  extrapolate: bool = False) -> CorrelationResult:
    """Average-to-fully-developed Nusselt ratio for turbulent gas entry flow.

    Bhatti-Shah as given by Lienhard Eq. (7.31): ``Nu_avg/Nu_inf = 1 + C/(L/D)^n``
    with ``C, n`` from Table 7.2, for ``Pr ~= 0.7`` and ``L/D > 3``.

    Parameters
    ----------
    length_to_diameter : float
        Heated length over diameter, ``> 3``.
    inlet : str
        ``"long_straight_pipe"``, ``"square_edged_inlet"``, ``"bend_180_deg"``,
        ``"bend_90_deg"`` or ``"sharp_elbow_90_deg"``.
    extrapolate : bool, optional
        Return a warned value for ``L/D <= 3`` instead of raising.

    Returns
    -------
    CorrelationResult
        Multiplier on the fully developed Nusselt number (average over
        ``0..L``), with an info warning restating the ``Pr ~= 0.7`` limit.

    References
    ----------
    .. [1] LIENHARD-AHTT-6.00. Lienhard and Lienhard, A Heat Transfer
           Textbook, 6th ed., 2024, Eq. (7.31) and Table 7.2, pp. 365-366,
           https://ahtt.mit.edu/
    """
    inlet = _choice(inlet, tuple(_ENTRY_INLETS), "inlet")
    ld = _positive(length_to_diameter, "length_to_diameter")
    warnings: list[WarningMessage] = []
    _check_range(ld, 3.0, None, "length_to_diameter", "turbulent entry (7.31)", extrapolate,
                 warnings, lo_open=True)
    warnings.append(WarningMessage("gas_prandtl_only",
                                   "Eq. (7.31) is for air or other fluids with Pr ~= 0.7.",
                                   "info"))
    c, n = _ENTRY_INLETS[inlet]
    return _result(1.0 + c / ld ** n, "multiplier", f"turbulent entry ({inlet})",
                   f"{_LIENHARD}, Eq. (7.31), Table 7.2", warnings)


def rough_pipe_nusselt(reynolds: float, prandtl: float, relative_roughness: float, *,
                       darcy_friction_factor: float | None = None,
                       extrapolate: bool = False) -> CorrelationResult:
    """Fully rough pipe Nusselt number (Bhatti-Shah), Lienhard Eq. (7.49).

    ``Nu_D = (f/8) Re_D Pr / [1 + sqrt(f/8)(4.5 Re_eps^0.2 Pr^0.5 - 8.48)]``
    for ``Re_D >= 1e4``, ``0.5 <= Pr <= 10``, ``0.002 <= eps/D <= 0.05`` in the
    fully rough regime ``Re_eps > 70``.

    Parameters
    ----------
    reynolds, prandtl : float
        Bulk Reynolds and Prandtl numbers.
    relative_roughness : float
        ``eps / D``.
    darcy_friction_factor : float, optional
        Darcy ``f``; defaults to Haaland's Eq. (7.50).
    extrapolate : bool, optional
        Return a warned value outside the range or regime instead of raising.

    Returns
    -------
    CorrelationResult
        Local ``Nu_D``. The source omits variable-property corrections for
        rough walls.

    References
    ----------
    .. [1] LIENHARD-AHTT-6.00. Lienhard and Lienhard, A Heat Transfer
           Textbook, 6th ed., 2024, Eqs. (7.48a)-(7.50) and Example 7.4,
           pp. 373-374, https://ahtt.mit.edu/
    """
    re = _positive(reynolds, "reynolds")
    pr = _positive(prandtl, "prandtl")
    rel = _positive(relative_roughness, "relative_roughness")
    warnings: list[WarningMessage] = []
    _check_range(re, 1.0e4, None, "reynolds", "Bhatti-Shah rough pipe", extrapolate, warnings)
    _check_range(pr, 0.5, 10.0, "prandtl", "Bhatti-Shah rough pipe", extrapolate, warnings)
    _check_range(rel, 0.002, 0.05, "relative_roughness", "Bhatti-Shah rough pipe", extrapolate,
                 warnings)
    if darcy_friction_factor is None:
        friction = haaland_friction_factor(re, rel, extrapolate=True)
        warnings.extend(friction.warnings)
        f = friction.value
    else:
        f = _positive(darcy_friction_factor, "darcy_friction_factor")
    re_eps = roughness_reynolds_number(re, rel, f)
    _check_range(re_eps, 70.0, None, "roughness_reynolds", "fully rough regime (7.48d)",
                 extrapolate, warnings, lo_open=True)
    value = (f / 8.0) * re * pr / (1.0 + sqrt(f / 8.0) * (4.5 * re_eps ** 0.2 * pr ** 0.5 - 8.48))
    return _result(value, "nusselt", "Bhatti-Shah fully rough pipe", f"{_LIENHARD}, Eq. (7.49)",
                   _dedupe(warnings))


def liquid_metal_nusselt(peclet: float, wall: WallCondition, *,
                         extrapolate: bool = False) -> CorrelationResult:
    """Fully developed turbulent liquid-metal pipe flow.

    Uniform heat flux: Skupinski et al., Eq. (7.55),
    ``Nu_D = 4.82 + 0.0185 Pe_D^0.827`` (fully developed, ``L/D >= 30``;
    68 % of accepted data within +/-20 %, ``Pe_D > 100`` per Taler); uniform
    wall temperature: Seban-Shimazaki, Eq. (7.54), ``5.0 + 0.025 Pe_D^0.8``
    (theoretical; no range stated).

    Parameters
    ----------
    peclet : float
        ``Pe_D = u_av D / alpha``.
    wall : {"uniform_temperature", "uniform_heat_flux"}
        Wall condition.
    extrapolate : bool, optional
        Return a warned value for ``Pe_D <= 100`` (uniform flux) instead of
        raising.

    Returns
    -------
    CorrelationResult
        ``Nu_D``.

    References
    ----------
    .. [1] LIENHARD-AHTT-6.00. Lienhard and Lienhard, A Heat Transfer
           Textbook, 6th ed., 2024, Eqs. (7.54)-(7.55) and Example 7.5,
           pp. 377-379, https://ahtt.mit.edu/
    """
    wall = _choice(wall, ("uniform_temperature", "uniform_heat_flux"), "wall")
    pe = _positive(peclet, "peclet")
    warnings: list[WarningMessage] = []
    if wall == "uniform_heat_flux":
        _check_range(pe, 100.0, None, "peclet", "Skupinski liquid metal", extrapolate, warnings,
                     lo_open=True)
        return _result(4.82 + 0.0185 * pe ** 0.827, "nusselt", "Skupinski liquid metal",
                       f"{_LIENHARD}, Eq. (7.55)", warnings)
    warnings.append(_range_not_stated("Seban-Shimazaki", "theoretical uniform-wall-temperature result"))
    return _result(5.0 + 0.025 * pe ** 0.8, "nusselt", "Seban-Shimazaki liquid metal",
                   f"{_LIENHARD}, Eq. (7.54)", warnings)


# ---------------------------------------------------------------------------
# Rocket coolant-side correlations
# ---------------------------------------------------------------------------


def taylor_temperature_ratio_exponent(distance_to_diameter: float) -> float:
    """Taylor's x/D-dependent exponent ``C2 = 0.57 - 1.59/(x/D)``, Eq. (19).

    Parameters
    ----------
    distance_to_diameter : float
        ``x / D`` from the start of the heated section.

    Returns
    -------
    float
        Exponent; the correlation multiplies by ``(T_s/T_b)^-C2``.

    References
    ----------
    .. [1] NASA-TN-D-4332-TAYLOR-1968. M. F. Taylor, NASA TN D-4332, 1968,
           Eq. (19) and Fig. 2, p. 12,
           https://ntrs.nasa.gov/citations/19680004771
    """
    return 0.57 - 1.59 / _positive(distance_to_diameter, "distance_to_diameter")


def taylor_hydrogen_nusselt(reynolds_bulk: float, prandtl_bulk: float,
                            wall_to_bulk_temperature_ratio: float, distance_to_diameter: float,
                            *, extrapolate: bool = False) -> CorrelationResult:
    """Taylor's single-phase hydrogen correlation, NASA TN D-4332 Eq. (20).

    ``Nu_b = 0.023 Re_b^0.8 Pr_b^0.4 (T_s/T_b)^-(0.57 - 1.59/(x/D))``,
    properties at the bulk temperature, ``Re_b = G D / mu_b``. Fitted to 3674
    hydrogen points in symmetrically heated straight tubes; range (Table I):
    ``2.0 <= x/D <= 252``, ``1.1 <= T_s/T_b <= 27.6``,
    ``7500 <= Re_b <= 1.38e7``. In region 1, 87 % of points lie within
    +/-25 %; in the near-critical region 2 only 40 % do (see
    :func:`taylor_near_critical_region`). Also predicted helium and nitrogen
    gas data (98 % within +/-25 %).

    Parameters
    ----------
    reynolds_bulk, prandtl_bulk : float
        Bulk Reynolds and Prandtl numbers.
    wall_to_bulk_temperature_ratio : float
        ``T_s / T_b`` in kelvin.
    distance_to_diameter : float
        ``x / D`` from the entrance of the heated section.
    extrapolate : bool, optional
        Return a warned value outside Table I instead of raising.

    Returns
    -------
    CorrelationResult
        ``Nu_b``.

    References
    ----------
    .. [1] NASA-TN-D-4332-TAYLOR-1968. M. F. Taylor, NASA TN D-4332, 1968,
           Summary p. 1, Table I p. 2, Eq. (20) p. 12, Table III p. 13, p. 14,
           https://ntrs.nasa.gov/citations/19680004771
    """
    re = _positive(reynolds_bulk, "reynolds_bulk")
    pr = _positive(prandtl_bulk, "prandtl_bulk")
    ratio = _positive(wall_to_bulk_temperature_ratio, "wall_to_bulk_temperature_ratio")
    xd = _positive(distance_to_diameter, "distance_to_diameter")
    warnings: list[WarningMessage] = []
    name = "Taylor hydrogen (TN D-4332)"
    _check_range(xd, 2.0, 252.0, "distance_to_diameter", name, extrapolate, warnings)
    _check_range(ratio, 1.1, 27.6, "wall_to_bulk_temperature_ratio", name, extrapolate, warnings)
    _check_range(re, 7500.0, 1.38e7, "reynolds_bulk", name, extrapolate, warnings)
    c2 = taylor_temperature_ratio_exponent(xd)
    value = 0.023 * re ** 0.8 * pr ** 0.4 * ratio ** (-c2)
    return _result(value, "nusselt", name, f"{_TAYLOR}, Eq. (20)", warnings)


def hess_kunz_nusselt(reynolds_film: float, prandtl_film: float,
                      wall_to_bulk_kinematic_viscosity_ratio: float, *,
                      cl_factor: float = 1.0) -> CorrelationResult:
    """Hess-Kunz supercritical-hydrogen film-temperature correlation.

    ``Nu_f = 0.0208 C_L Re_f^0.8 Pr_f^0.4 (1 + 0.01457 nu_w/nu_b)`` with
    ``Re_f = rho_f V_b D / mu_f`` and film properties at
    ``T_f = (T_w + T_b)/2``. ``C_L = 1`` is the original equation as reported
    by Taylor (Eq. 14); ``C_L != 1`` is the CR-678 form (Eq. 3) recommended
    for 600-1500 psia with ``C_L`` from :func:`hess_kunz_cl_factor`.

    Parameters
    ----------
    reynolds_film, prandtl_film : float
        Film Reynolds and Prandtl numbers.
    wall_to_bulk_kinematic_viscosity_ratio : float
        ``nu_w / nu_b``.
    cl_factor : float, optional
        Multiplying coefficient ``C_L``.

    Returns
    -------
    CorrelationResult
        ``Nu_f`` with a ``range_not_stated`` info warning (the original Hess
        and Kunz range is not given in the reports read).

    References
    ----------
    .. [1] NASA-TN-D-4332-TAYLOR-1968. M. F. Taylor, NASA TN D-4332, 1968,
           Eq. (14) and Symbols, pp. 3, 7,
           https://ntrs.nasa.gov/citations/19680004771
    .. [2] NASA-CR-678-AEROJET-1967. Aerojet-General Corp., NASA CR-678,
           1967, Eqs. (1) and (3), pp. 2-3,
           https://ntrs.nasa.gov/citations/19670008827
    """
    re = _positive(reynolds_film, "reynolds_film")
    pr = _positive(prandtl_film, "prandtl_film")
    ratio = _positive(wall_to_bulk_kinematic_viscosity_ratio,
                      "wall_to_bulk_kinematic_viscosity_ratio")
    cl = _positive(cl_factor, "cl_factor")
    value = 0.0208 * cl * re ** 0.8 * pr ** 0.4 * (1.0 + 0.01457 * ratio)
    return _result(value, "nusselt", "Hess-Kunz", f"{_TAYLOR}, Eq. (14); {_CR678}, Eq. (3)",
                   [_range_not_stated("Hess-Kunz", "reported without the original data range")])


_CL_RANKINE = (50.0, 55.0, 60.0, 65.0, 70.0, 75.0, 80.0, 85.0)
_CL_STRAIGHT = (2.0, 1.73, 1.48, 1.26, 1.07, 0.93, 0.87, 0.85)
_CL_CURVED = (None, None, 1.95, 1.64, 1.37, 1.21, 1.11, 1.10)


def hess_kunz_cl_factor(bulk_temperature_k: float, geometry: Literal["straight", "curved"],
                        *, extrapolate: bool = False) -> CorrelationResult:
    """Coefficient ``C_L`` of the modified Hess-Kunz equation, CR-678 p. 3.

    Linear interpolation (this module's choice) in the CR-678 table of
    ``C_L`` against coolant temperature: 50-85 degR (27.8-47.2 K) for
    straight tubes, 60-85 degR (33.3-47.2 K) for curved tubes.

    Parameters
    ----------
    bulk_temperature_k : float
        Coolant bulk temperature.
    geometry : {"straight", "curved"}
        Tube geometry column.
    extrapolate : bool, optional
        Return the end value with a warning outside the table instead of
        raising.

    Returns
    -------
    CorrelationResult
        ``C_L``.

    References
    ----------
    .. [1] NASA-CR-678-AEROJET-1967. Aerojet-General Corp., NASA CR-678,
           1967, Eq. (3) and table, p. 3,
           https://ntrs.nasa.gov/citations/19670008827
    """
    geometry = _choice(geometry, ("straight", "curved"), "geometry")
    t_r = _positive(bulk_temperature_k, "bulk_temperature_k") / _RANKINE_TO_KELVIN
    column = _CL_STRAIGHT if geometry == "straight" else _CL_CURVED
    points = [(t, c) for t, c in zip(_CL_RANKINE, column, strict=True) if c is not None]
    warnings: list[WarningMessage] = []
    lo, hi = points[0][0], points[-1][0]
    _check_range(t_r * _RANKINE_TO_KELVIN, lo * _RANKINE_TO_KELVIN, hi * _RANKINE_TO_KELVIN,
                 "bulk_temperature_k", f"CR-678 C_L table ({geometry})", extrapolate, warnings)
    t_r = min(max(t_r, lo), hi)
    for (t0, c0), (t1, c1) in itertools.pairwise(points):
        if t0 <= t_r <= t1:
            value = c0 + (c1 - c0) * (t_r - t0) / (t1 - t0)
            break
    return _result(value, "multiplier", f"Hess-Kunz C_L ({geometry})", f"{_CR678}, p. 3",
                   warnings)


def sp125_sieder_tate_nusselt(reynolds: float, prandtl: float, bulk_to_wall_viscosity_ratio: float,
                              coefficient: float) -> CorrelationResult:
    """NASA SP-125 coolant-side Sieder-Tate form, Eq. (4-23).

    ``Nu = C1 Re^0.8 Pr^0.4 (mu/mu_w)^0.14`` for non-boiling, subcritical
    coolant temperatures (regions A1-A2 and B1-B2 of SP-125 Fig. 4-27);
    ``C1`` is coolant specific and must be supplied.

    Parameters
    ----------
    reynolds, prandtl : float
        Bulk Reynolds and Prandtl numbers.
    bulk_to_wall_viscosity_ratio : float
        ``mu / mu_w``.
    coefficient : float
        ``C1`` (not tabulated in SP-125).

    Returns
    -------
    CorrelationResult
        ``Nu`` with a ``range_not_stated`` info warning.

    References
    ----------
    .. [1] NASA-SP125-COOLANT-SIDE. Huzel and Huang, NASA SP-125, 1971,
           printed pp. 105-106, Eq. (4-23),
           https://ntrs.nasa.gov/api/citations/19710019929/downloads/19710019929.pdf
    """
    re = _positive(reynolds, "reynolds")
    pr = _positive(prandtl, "prandtl")
    ratio = _positive(bulk_to_wall_viscosity_ratio, "bulk_to_wall_viscosity_ratio")
    c1 = _positive(coefficient, "coefficient")
    return _result(c1 * re ** 0.8 * pr ** 0.4 * ratio ** 0.14, "nusselt",
                   "SP-125 Sieder-Tate", f"{_SP125}, Eq. (4-23)",
                   [_range_not_stated("SP-125 Eq. (4-23)", "non-boiling subcritical regions only")])


def sp125_supercritical_film_coefficient(specific_heat_j_kg_k: float, viscosity_pa_s: float,
                                         prandtl: float, mass_flux_kg_m2_s: float,
                                         hydraulic_diameter_m: float, bulk_temperature_k: float,
                                         wall_temperature_k: float) -> CorrelationResult:
    """NASA SP-125 vapour-film coefficient, Eq. (4-25), as printed.

    ``h_c = 0.029 c_p mu^0.2 / Pr^(2/3) (G^0.8 / d^0.2) (T_co/T_wc)^0.55`` for
    coolant at supercritical pressure and temperature (region B2-B3 of
    SP-125 Fig. 4-27). The equation is dimensionally homogeneous, so SI
    inputs give W/(m^2 K); in Nusselt form it is
    ``Nu = 0.029 Re^0.8 Pr^(1/3) (T_co/T_wc)^0.55``.

    Parameters
    ----------
    specific_heat_j_kg_k, viscosity_pa_s, prandtl : float
        Coolant ``c_p``, ``mu`` and ``Pr`` (SP-125 does not state the
        reference temperature; :func:`sp125_supercritical_film_heat_transfer`
        uses bulk values).
    mass_flux_kg_m2_s : float
        Coolant mass flow per unit area ``G``.
    hydraulic_diameter_m : float
        Passage hydraulic diameter ``d``.
    bulk_temperature_k, wall_temperature_k : float
        ``T_co`` and ``T_wc``.

    Returns
    -------
    CorrelationResult
        Coefficient in W/(m^2 K), with a ``range_not_stated`` info warning.

    References
    ----------
    .. [1] NASA-SP125-COOLANT-SIDE. Huzel and Huang, NASA SP-125, 1971,
           printed p. 106, Eq. (4-25),
           https://ntrs.nasa.gov/api/citations/19710019929/downloads/19710019929.pdf
    """
    cp = _positive(specific_heat_j_kg_k, "specific_heat_j_kg_k")
    mu = _positive(viscosity_pa_s, "viscosity_pa_s")
    pr = _positive(prandtl, "prandtl")
    g = _positive(mass_flux_kg_m2_s, "mass_flux_kg_m2_s")
    d = _positive(hydraulic_diameter_m, "hydraulic_diameter_m")
    tco = _positive(bulk_temperature_k, "bulk_temperature_k")
    twc = _positive(wall_temperature_k, "wall_temperature_k")
    value = 0.029 * cp * mu ** 0.2 / pr ** (2.0 / 3.0) * (g ** 0.8 / d ** 0.2) * (tco / twc) ** 0.55
    return _result(value, "heat_transfer_coefficient_w_m2_k", "SP-125 supercritical vapour film",
                   f"{_SP125}, Eq. (4-25)",
                   [_range_not_stated("SP-125 Eq. (4-25)", "supercritical region B2-B3 only")])


def curvature_heat_transfer_factor(reynolds: float, hydraulic_radius_m: float,
                                   radius_of_curvature_m: float,
                                   side: Literal["concave", "convex"],
                                   coolant_state: Literal["gas", "liquid", "supercritical"]
                                   ) -> CorrelationResult:
    """Coolant-side curvature multiplier used by NASA's RTE code.

    ``phi_cur = [Re (r_c/R_cur)^2]^(+/-0.05)``, ``+`` on a concave (outer)
    wall, ``-`` on a convex wall, Wadel Eq. (11) (refs. Niino et al. and
    Owhadi et al. therein; the same group as Ito's friction ratio).

    Parameters
    ----------
    reynolds : float
        Coolant Reynolds number at the reference state.
    hydraulic_radius_m : float
        ``r_c``, half the coolant-channel hydraulic diameter.
    radius_of_curvature_m : float
        ``R_cur``.
    side : {"concave", "convex"}
        Heated wall relative to the bend.
    coolant_state : {"gas", "liquid", "supercritical"}
        Used for the SP-8087 guidance: for liquids an enhancement factor is
        not recommended (warning on a concave wall).

    Returns
    -------
    CorrelationResult
        Multiplier on the straight-channel Nusselt number.

    Notes
    -----
    NASA SP-8087 recommends an experimental enhancement factor for gases
    when curvature forces coolant against the heated wall, no factor for
    liquids, and a reduction when the bulk is forced away from the heated
    wall [2]_. CR-678 measured concave-side enhancement up to 2:1 in curved
    hydrogen tubes [3]_.

    References
    ----------
    .. [1] NASA-TM-1998-206313-WADEL. M. F. Wadel, NASA/TM-1998-206313,
           Eq. (11) and nomenclature, pp. vii, 16,
           https://ntrs.nasa.gov/citations/19980017619
    .. [2] NASA-SP-8087-1972. NASA SP-8087, 1972, sec. 3.1.1.5.2, p. 61,
           https://ntrs.nasa.gov/citations/19730022965
    .. [3] NASA-CR-678-AEROJET-1967. Aerojet-General Corp., NASA CR-678,
           1967, pp. 3-4 and 53-54,
           https://ntrs.nasa.gov/citations/19670008827
    """
    side = _choice(side, ("concave", "convex"), "side")
    state = _choice(coolant_state, ("gas", "liquid", "supercritical"), "coolant_state")
    re = _positive(reynolds, "reynolds")
    ratio = _positive(hydraulic_radius_m, "hydraulic_radius_m") / _positive(
        radius_of_curvature_m, "radius_of_curvature_m")
    group = re * ratio ** 2
    warnings = [_range_not_stated("RTE curvature factor", "no Re or r/R limits are given")]
    if side == "concave" and state == "liquid":
        warnings.append(WarningMessage(
            "curvature_factor_not_recommended",
            "NASA SP-8087 p. 61: use of a curvature enhancement factor is not recommended "
            "for liquids.", "warning", "coolant_state"))
    if group < 1.0:
        warnings.append(WarningMessage(
            "curvature_group_below_one",
            "Re (r/R)^2 < 1 inverts the sign of the curvature effect.", "warning", "reynolds"))
    exponent = 0.05 if side == "concave" else -0.05
    return _result(group ** exponent, "multiplier", f"RTE curvature ({side})",
                   f"{_WADEL}, Eq. (11)", warnings)


def rte_coolant_nusselt(reynolds: float, prandtl: float, *, coefficient: float = 0.023,
                        curvature_factor: float = 1.0,
                        entrance_factor: float = 1.0) -> CorrelationResult:
    """Coolant Nusselt number of NASA's RTE code, Wadel Eq. (6).

    ``Nu = C_c Re^0.8 Pr^0.4 phi_ent phi_cur`` with reference-state properties
    and ``Re`` corrected by Wadel Eq. (7); ``C_c = 0.023`` was used in the
    TM. ``phi_cur`` comes from :func:`curvature_heat_transfer_factor`; the
    TM's entrance factor (Eq. 10) is printed with a ``+1`` in the numerator
    whose meaning is not defined, so it is not implemented and is taken as
    an input here.

    Parameters
    ----------
    reynolds, prandtl : float
        Reference-state Reynolds and Prandtl numbers.
    coefficient : float, optional
        ``C_c`` (user input from data in RTE).
    curvature_factor, entrance_factor : float, optional
        Multipliers ``phi_cur`` and ``phi_ent``.

    Returns
    -------
    CorrelationResult
        ``Nu``, with a ``range_not_stated`` info warning.

    References
    ----------
    .. [1] NASA-TM-1998-206313-WADEL. M. F. Wadel, NASA/TM-1998-206313,
           Eqs. (6)-(10), pp. 15-16, https://ntrs.nasa.gov/citations/19980017619
    """
    re = _positive(reynolds, "reynolds")
    pr = _positive(prandtl, "prandtl")
    value = (_positive(coefficient, "coefficient") * re ** 0.8 * pr ** 0.4
             * _positive(curvature_factor, "curvature_factor")
             * _positive(entrance_factor, "entrance_factor"))
    return _result(value, "nusselt", "RTE coolant Nusselt", f"{_WADEL}, Eq. (6)",
                   [_range_not_stated("RTE coolant correlation", "C_c is user input from data")])


# ---------------------------------------------------------------------------
# Real-fluid evaluation
# ---------------------------------------------------------------------------


def _fluid_key(fluid: object, allowed: tuple[str, ...]) -> str:
    """Return a reference-fluid key accepted by the correlation."""
    return _choice(fluid, allowed, "fluid")


def _state(fluid: str, temperature_k: float, pressure_pa: float) -> FluidState:
    """Equation-of-state state at ``(T, p)`` from the package's Helmholtz models."""
    return reference_fluid(fluid).state_tp(temperature_k, pressure_pa)


def pseudocritical_temperature(fluid: str, pressure_pa: float) -> float:
    """Transposed (pseudo-)critical temperature ``T*``: the maximum of ``c_p(T)``.

    Taylor defines ``T*`` as the temperature at which ``c_p`` of the fluid at
    constant pressure reaches a maximum. It is located with the package's
    reference equation of state by a scan from ``T_c`` to ``3 T_c`` followed
    by golden-section refinement.

    Parameters
    ----------
    fluid : str
        Reference-fluid key (see
        :data:`~rocket_propulsion.thermodynamics.real_fluids.REFERENCE_FLUID_KEYS`).
    pressure_pa : float
        Pressure; must exceed the critical pressure.

    Returns
    -------
    float
        ``T*`` in K.

    Raises
    ------
    DomainError
        If the pressure is not supercritical or no interior maximum is found.

    References
    ----------
    .. [1] NASA-TN-D-4332-TAYLOR-1968. M. F. Taylor, NASA TN D-4332, 1968,
           Symbols, p. 3, https://ntrs.nasa.gov/citations/19680004771
    .. [2] LEACHMAN-2009-H2. Leachman et al., J. Phys. Chem. Ref. Data 38,
           721 (2009), https://srd.nist.gov/jpcrdreprint/1.3160306.pdf
    """
    model = reference_fluid_model(fluid)
    p = _positive(pressure_pa, "pressure_pa")
    if p <= model.p_crit_pa:
        raise DomainError("T* is defined only above the critical pressure.", field="pressure_pa")
    calculator = reference_fluid(fluid)

    def cp(t: float) -> float:
        return calculator.state_tp(t, p).cp_j_kg_k

    lo, hi = model.t_crit_k, min(3.0 * model.t_crit_k, model.t_max_k)
    count = 120
    grid = [lo + (hi - lo) * i / count for i in range(count + 1)]
    values = [cp(t) for t in grid]
    k = max(range(len(values)), key=values.__getitem__)
    if k in (0, count):
        raise DomainError("No interior c_p maximum between T_c and 3 T_c.",
                          code="non_convergence", field="pressure_pa")
    a, b = grid[k - 1], grid[k + 1]
    ratio = (sqrt(5.0) - 1.0) / 2.0
    c, d = b - ratio * (b - a), a + ratio * (b - a)
    fc, fd = cp(c), cp(d)
    while b - a > 1.0e-6 * model.t_crit_k:
        if fc > fd:
            b, d, fd = d, c, fc
            c = b - ratio * (b - a)
            fc = cp(c)
        else:
            a, c, fc = c, d, fd
            d = a + ratio * (b - a)
            fd = cp(d)
    return 0.5 * (a + b)


def taylor_near_critical_region(inlet_temperature_k: float, inlet_pressure_pa: float,
                                fluid: str = "parahydrogen") -> NearCriticalScreen:
    """Screen a hydrogen inlet state against Taylor's near-critical region 2.

    Region 2 is ``45 degR (25 K) < T_i < T*`` and ``p_c < p_i < 530 psia
    (3.65 MN/m^2)``; there only 40 % of 948 points were within +/-25 % of
    Eq. (20), against 87 % elsewhere.

    Parameters
    ----------
    inlet_temperature_k, inlet_pressure_pa : float
        Inlet state of the heated passage.
    fluid : {"parahydrogen", "normal_hydrogen"}, optional
        Equation of state used for ``p_c`` and ``T*``.

    Returns
    -------
    NearCriticalScreen
        Screen result.

    References
    ----------
    .. [1] NASA-TN-D-4332-TAYLOR-1968. M. F. Taylor, NASA TN D-4332, 1968,
           pp. 1, 13-14, https://ntrs.nasa.gov/citations/19680004771
    """
    key = _fluid_key(fluid, _HYDROGEN_KEYS)
    t_i = _positive(inlet_temperature_k, "inlet_temperature_k")
    p_i = _positive(inlet_pressure_pa, "inlet_pressure_pa")
    p_c = reference_fluid_model(key).p_crit_pa
    t_lower = 45.0 * _RANKINE_TO_KELVIN
    p_upper = 530.0 * _PSI_PA
    if p_i <= p_c:
        return NearCriticalScreen(False, None, p_c)
    t_star = pseudocritical_temperature(key, p_i)
    inside = t_lower < t_i < t_star and p_i < p_upper
    warnings: tuple[WarningMessage, ...] = ()
    if inside:
        warnings = (WarningMessage(
            "taylor_near_critical_region",
            "Inlet state is in Taylor's near-critical region 2, where only 40 % of 948 "
            "points were within +/-25 % of TN D-4332 Eq. (20).", "warning"),)
    return NearCriticalScreen(inside, t_star, p_c, warnings)


def _state_warnings(*states: FluidState) -> list[WarningMessage]:
    """Collect the unique equation-of-state warnings of several states."""
    out: list[WarningMessage] = []
    for state in states:
        out.extend(state.warnings)
    return _dedupe(out)


def _dedupe(warnings: list[WarningMessage]) -> list[WarningMessage]:
    """Drop repeated warnings while keeping order."""
    seen: set[tuple[str, str]] = set()
    out = []
    for warning in warnings:
        key = (warning.code, warning.message)
        if key not in seen:
            seen.add(key)
            out.append(warning)
    return out


def taylor_hydrogen_coefficient(*, fluid: str, bulk_temperature_k: float,
                                wall_temperature_k: float, pressure_pa: float,
                                mass_flux_kg_m2_s: float, diameter_m: float,
                                distance_from_entrance_m: float, bulk_viscosity_pa_s: float,
                                bulk_conductivity_w_m_k: float,
                                extrapolate: bool = False) -> CoolantHeatTransfer:
    """Taylor's Eq. (20) evaluated with real-fluid bulk heat capacity.

    ``c_p`` and density at ``(T_b, p)`` come from the reference equation of
    state; viscosity and conductivity are caller inputs. Taylor reports the
    correlation for hydrogen and also for helium and nitrogen gas.

    Parameters
    ----------
    fluid : {"parahydrogen", "normal_hydrogen", "helium", "nitrogen"}
        Reference fluid.
    bulk_temperature_k, wall_temperature_k : float
        ``T_b`` and ``T_s``.
    pressure_pa : float
        Static pressure.
    mass_flux_kg_m2_s, diameter_m : float
        ``G`` and tube inside diameter ``D``.
    distance_from_entrance_m : float
        ``x`` from the entrance of the heated section.
    bulk_viscosity_pa_s, bulk_conductivity_w_m_k : float
        Transport properties at ``T_b``.
    extrapolate : bool, optional
        Passed to :func:`taylor_hydrogen_nusselt`.

    Returns
    -------
    CoolantHeatTransfer
        Coefficient and the groups used.

    References
    ----------
    .. [1] NASA-TN-D-4332-TAYLOR-1968. M. F. Taylor, NASA TN D-4332, 1968,
           Eq. (20), p. 12, https://ntrs.nasa.gov/citations/19680004771
    .. [2] LEACHMAN-2009-H2. Leachman et al., J. Phys. Chem. Ref. Data 38,
           721 (2009), https://srd.nist.gov/jpcrdreprint/1.3160306.pdf
    """
    key = _fluid_key(fluid, (*_HYDROGEN_KEYS, "helium", "nitrogen"))
    tb = _positive(bulk_temperature_k, "bulk_temperature_k")
    tw = _positive(wall_temperature_k, "wall_temperature_k")
    d = _positive(diameter_m, "diameter_m")
    k_b = _positive(bulk_conductivity_w_m_k, "bulk_conductivity_w_m_k")
    mu_b = _positive(bulk_viscosity_pa_s, "bulk_viscosity_pa_s")
    state = _state(key, tb, _positive(pressure_pa, "pressure_pa"))
    re = mass_flux_reynolds_number(mass_flux_kg_m2_s, d, mu_b)
    pr = prandtl_number(state.cp_j_kg_k, mu_b, k_b)
    nu = taylor_hydrogen_nusselt(re, pr, tw / tb,
                                 _positive(distance_from_entrance_m, "distance_from_entrance_m") / d,
                                 extrapolate=extrapolate)
    return CoolantHeatTransfer(nu.value * k_b / d, nu.value, re, pr, tb, state.density_kg_m3,
                               state.cp_j_kg_k, nu.correlation, nu.source,
                               (*nu.warnings, *_state_warnings(state)))


def hess_kunz_coefficient(*, fluid: str, bulk_temperature_k: float, wall_temperature_k: float,
                          pressure_pa: float, mass_flux_kg_m2_s: float, diameter_m: float,
                          film_viscosity_pa_s: float, film_conductivity_w_m_k: float,
                          bulk_viscosity_pa_s: float, wall_viscosity_pa_s: float,
                          cl_factor: float = 1.0) -> CoolantHeatTransfer:
    """Hess-Kunz equation evaluated with real-fluid densities and film ``c_p``.

    Film temperature ``T_f = (T_w + T_b)/2`` (CR-678 p. 2);
    ``V_b = G / rho_b``; ``Re_f = rho_f V_b D / mu_f`` (Taylor symbols);
    ``nu_w/nu_b = (mu_w/rho_w)/(mu_b/rho_b)``. Densities at ``T_b``, ``T_f``,
    ``T_w`` and ``c_p`` at ``T_f`` come from the equation of state; the four
    transport properties are caller inputs.

    Parameters
    ----------
    fluid : {"parahydrogen", "normal_hydrogen"}
        Hydrogen equation of state.
    bulk_temperature_k, wall_temperature_k, pressure_pa : float
        State.
    mass_flux_kg_m2_s, diameter_m : float
        ``G`` and ``D``.
    film_viscosity_pa_s, film_conductivity_w_m_k : float
        ``mu_f`` and ``k_f`` at ``T_f``.
    bulk_viscosity_pa_s, wall_viscosity_pa_s : float
        ``mu_b`` and ``mu_w``.
    cl_factor : float, optional
        ``C_L`` (see :func:`hess_kunz_cl_factor`).

    Returns
    -------
    CoolantHeatTransfer
        Coefficient ``h = Nu_f k_f / D``.

    References
    ----------
    .. [1] NASA-TN-D-4332-TAYLOR-1968. M. F. Taylor, NASA TN D-4332, 1968,
           Eq. (14), Symbols, pp. 3, 7,
           https://ntrs.nasa.gov/citations/19680004771
    .. [2] NASA-CR-678-AEROJET-1967. Aerojet-General Corp., NASA CR-678,
           1967, Eqs. (1) and (3), pp. 2-3,
           https://ntrs.nasa.gov/citations/19670008827
    .. [3] LEACHMAN-2009-H2. Leachman et al., J. Phys. Chem. Ref. Data 38,
           721 (2009), https://srd.nist.gov/jpcrdreprint/1.3160306.pdf
    """
    key = _fluid_key(fluid, _HYDROGEN_KEYS)
    tb = _positive(bulk_temperature_k, "bulk_temperature_k")
    tw = _positive(wall_temperature_k, "wall_temperature_k")
    p = _positive(pressure_pa, "pressure_pa")
    g = _positive(mass_flux_kg_m2_s, "mass_flux_kg_m2_s")
    d = _positive(diameter_m, "diameter_m")
    mu_f = _positive(film_viscosity_pa_s, "film_viscosity_pa_s")
    k_f = _positive(film_conductivity_w_m_k, "film_conductivity_w_m_k")
    mu_b = _positive(bulk_viscosity_pa_s, "bulk_viscosity_pa_s")
    mu_w = _positive(wall_viscosity_pa_s, "wall_viscosity_pa_s")
    tf = 0.5 * (tw + tb)
    bulk, film, wall = _state(key, tb, p), _state(key, tf, p), _state(key, tw, p)
    velocity = g / bulk.density_kg_m3
    re_f = reynolds_number(film.density_kg_m3, velocity, d, mu_f)
    pr_f = prandtl_number(film.cp_j_kg_k, mu_f, k_f)
    ratio = (mu_w / wall.density_kg_m3) / (mu_b / bulk.density_kg_m3)
    nu = hess_kunz_nusselt(re_f, pr_f, ratio, cl_factor=cl_factor)
    return CoolantHeatTransfer(nu.value * k_f / d, nu.value, re_f, pr_f, tf, film.density_kg_m3,
                               film.cp_j_kg_k, nu.correlation, nu.source,
                               (*nu.warnings, *_state_warnings(bulk, film, wall)))


def sp125_supercritical_film_heat_transfer(*, fluid: str, bulk_temperature_k: float,
                                           wall_temperature_k: float, pressure_pa: float,
                                           mass_flux_kg_m2_s: float, hydraulic_diameter_m: float,
                                           viscosity_pa_s: float, conductivity_w_m_k: float
                                           ) -> CoolantHeatTransfer:
    """SP-125 Eq. (4-25) with the bulk ``c_p`` from the equation of state.

    SP-125 does not state where ``c_p``, ``mu`` and ``Pr`` are evaluated; this
    wrapper uses the bulk state (an assumption, reported as an info warning).

    Parameters
    ----------
    fluid : str
        Reference-fluid key.
    bulk_temperature_k, wall_temperature_k, pressure_pa : float
        State.
    mass_flux_kg_m2_s, hydraulic_diameter_m : float
        ``G`` and ``d``.
    viscosity_pa_s, conductivity_w_m_k : float
        Bulk transport properties (caller inputs).

    Returns
    -------
    CoolantHeatTransfer
        Coefficient and the groups used.

    References
    ----------
    .. [1] NASA-SP125-COOLANT-SIDE. Huzel and Huang, NASA SP-125, 1971,
           printed p. 106, Eq. (4-25),
           https://ntrs.nasa.gov/api/citations/19710019929/downloads/19710019929.pdf
    """
    tb = _positive(bulk_temperature_k, "bulk_temperature_k")
    model = reference_fluid_model(fluid)
    p = _positive(pressure_pa, "pressure_pa")
    state = _state(fluid, tb, p)
    mu = _positive(viscosity_pa_s, "viscosity_pa_s")
    k = _positive(conductivity_w_m_k, "conductivity_w_m_k")
    d = _positive(hydraulic_diameter_m, "hydraulic_diameter_m")
    pr = prandtl_number(state.cp_j_kg_k, mu, k)
    h = sp125_supercritical_film_coefficient(state.cp_j_kg_k, mu, pr, mass_flux_kg_m2_s, d, tb,
                                             wall_temperature_k)
    warnings = list(h.warnings)
    warnings.append(WarningMessage("bulk_reference_state_assumed",
                                   "SP-125 Eq. (4-25) does not state the property reference "
                                   "temperature; bulk properties were used.", "info"))
    if p <= model.p_crit_pa or tb <= model.t_crit_k:
        warnings.append(WarningMessage(
            "not_supercritical",
            "SP-125 Eq. (4-25) is for coolant above critical pressure and temperature.",
            "warning", "pressure_pa"))
    re = mass_flux_reynolds_number(mass_flux_kg_m2_s, d, mu)
    return CoolantHeatTransfer(h.value, h.value * d / k, re, pr, tb, state.density_kg_m3,
                               state.cp_j_kg_k, h.correlation, h.source,
                               (*warnings, *_state_warnings(state)))


# ---------------------------------------------------------------------------
# External forced convection and natural convection
# ---------------------------------------------------------------------------


def flat_plate_laminar_nusselt(reynolds_x: float, prandtl: float, *, average: bool = False,
                               extrapolate: bool = False) -> CorrelationResult:
    """Laminar flat plate, uniform wall temperature, Lienhard Eqs. (6.58)/(6.68).

    Local ``Nu_x = 0.332 Re_x^(1/2) Pr^(1/3)``; average over ``0..L``
    ``Nu_L = 0.664 Re_L^(1/2) Pr^(1/3)``; ``Pr >= 0.6`` (within 2 % of the
    exact solution). Properties at the film temperature. The caller must
    ensure the boundary layer is laminar.

    Parameters
    ----------
    reynolds_x, prandtl : float
        Reynolds number on ``x`` (or ``L``) and Prandtl number.
    average : bool, optional
        Return the plate-average value.
    extrapolate : bool, optional
        Return a warned value for ``Pr < 0.6`` instead of raising.

    Returns
    -------
    CorrelationResult
        Nusselt number.

    References
    ----------
    .. [1] LIENHARD-AHTT-6.00. Lienhard and Lienhard, A Heat Transfer
           Textbook, 6th ed., 2024, Eqs. (6.58) and (6.68), pp. 308-311,
           https://ahtt.mit.edu/
    """
    re = _positive(reynolds_x, "reynolds_x")
    pr = _positive(prandtl, "prandtl")
    warnings: list[WarningMessage] = []
    _check_range(pr, 0.6, None, "prandtl", "laminar flat plate", extrapolate, warnings)
    coefficient, eq = (0.664, "Eq. (6.68)") if average else (0.332, "Eq. (6.58)")
    return _result(coefficient * sqrt(re) * pr ** (1.0 / 3.0), "nusselt", "laminar flat plate",
                   f"{_LIENHARD}, {eq}", warnings)


def flat_plate_turbulent_nusselt(reynolds_x: float, prandtl: float, *,
                                 extrapolate: bool = False) -> CorrelationResult:
    """Local turbulent flat-plate Nusselt number, Lienhard Eqs. (6.102) and (6.111).

    ``C_f = 0.455/[ln(0.06 Re_x)]^2`` (White) and
    ``Nu_x = Re_x Pr (C_f/2) / [1 + 12.7 (Pr^(2/3) - 1) sqrt(C_f/2)]`` for
    ``Pr >= 0.5``, smooth walls, zero or mild pressure gradient; air data
    within +/-15 % (standard deviation 5.5 %).

    Parameters
    ----------
    reynolds_x, prandtl : float
        Local Reynolds and Prandtl numbers.
    extrapolate : bool, optional
        Return a warned value for ``Pr < 0.5`` instead of raising.

    Returns
    -------
    CorrelationResult
        Local Nusselt number.

    References
    ----------
    .. [1] LIENHARD-AHTT-6.00. Lienhard and Lienhard, A Heat Transfer
           Textbook, 6th ed., 2024, Eqs. (6.102) and (6.111), pp. 322 and 325,
           https://ahtt.mit.edu/
    """
    re = _positive(reynolds_x, "reynolds_x")
    pr = _positive(prandtl, "prandtl")
    warnings: list[WarningMessage] = []
    _check_range(pr, 0.5, None, "prandtl", "turbulent flat plate (6.111)", extrapolate, warnings)
    argument = 0.06 * re
    if argument <= 1.0:
        raise DomainError("White's C_f needs 0.06 Re_x > 1.", code="outside_correlation_range",
                          field="reynolds_x")
    cf = 0.455 / log(argument) ** 2
    value = re * pr * (cf / 2.0) / (1.0 + 12.7 * (pr ** (2.0 / 3.0) - 1.0) * sqrt(cf / 2.0))
    return _result(value, "nusselt", "turbulent flat plate (Reynolds-Colburn analogy)",
                   f"{_LIENHARD}, Eqs. (6.102), (6.111)", warnings)


def flat_plate_turbulent_power_law_nusselt(reynolds_x: float, prandtl: float,
                                           fluid: Literal["gas", "nonmetallic_liquid"]
                                           ) -> CorrelationResult:
    """Turbulent flat-plate power laws, Lienhard Eqs. (6.112)-(6.113).

    Gases (Reynolds et al.): ``Nu_x = 0.0296 Re_x^0.8 Pr^0.6``; nonmetallic
    liquids (Zukauskas, ``Pr`` up to 85, about +/-15 %):
    ``Nu_x = 0.032 Re_x^0.8 Pr^0.43``.

    Parameters
    ----------
    reynolds_x, prandtl : float
        Local Reynolds and Prandtl numbers.
    fluid : {"gas", "nonmetallic_liquid"}
        Fluid class.

    Returns
    -------
    CorrelationResult
        Local Nusselt number (info warning: no numerical range stated for gases).

    References
    ----------
    .. [1] LIENHARD-AHTT-6.00. Lienhard and Lienhard, A Heat Transfer
           Textbook, 6th ed., 2024, Eqs. (6.112)-(6.113), pp. 326-327,
           https://ahtt.mit.edu/
    """
    fluid = _choice(fluid, ("gas", "nonmetallic_liquid"), "fluid")
    re = _positive(reynolds_x, "reynolds_x")
    pr = _positive(prandtl, "prandtl")
    if fluid == "gas":
        return _result(0.0296 * re ** 0.8 * pr ** 0.6, "nusselt", "turbulent flat plate, gases",
                       f"{_LIENHARD}, Eq. (6.112)",
                       [_range_not_stated("Eq. (6.112)", "gases with Pr near that of air")])
    warnings: list[WarningMessage] = []
    _check_range(pr, None, 85.0, "prandtl", "Zukauskas liquids", False, warnings)
    return _result(0.032 * re ** 0.8 * pr ** 0.43, "nusselt", "turbulent flat plate, liquids",
                   f"{_LIENHARD}, Eq. (6.113)")


def churchill_bernstein_cylinder_nusselt(reynolds: float, prandtl: float,
                                         form: Literal["full", "low_re", "mid_re"] = "full", *,
                                         extrapolate: bool = False) -> CorrelationResult:
    """Cylinder in cross flow (Churchill-Bernstein), Lienhard Eqs. (7.65)-(7.68).

    ``Nu_D = 0.3 + 0.62 Re^(1/2) Pr^(1/3)/[1 + (0.4/Pr)^(2/3)]^(1/4) * B`` with
    ``B = [1 + (Re/282000)^(5/8)]^(4/5)`` (``"full"``, Eq. 7.65, ``Re Pr >= 0.2``),
    ``B = 1`` (``"low_re"``, Eq. 7.66, ``Re < 4000``) or
    ``B = 1 + (Re/282000)^(1/2)`` (``"mid_re"``, Eq. 7.68,
    ``40000 < Re < 400000``). Properties at the film temperature.

    Parameters
    ----------
    reynolds, prandtl : float
        ``Re_D`` and ``Pr``.
    form : {"full", "low_re", "mid_re"}, optional
        Which form.
    extrapolate : bool, optional
        Return a warned value outside the form's range instead of raising.

    Returns
    -------
    CorrelationResult
        Average ``Nu_D``.

    References
    ----------
    .. [1] LIENHARD-AHTT-6.00. Lienhard and Lienhard, A Heat Transfer
           Textbook, 6th ed., 2024, Eqs. (7.65)-(7.68) and Example 7.8,
           pp. 390-392, https://ahtt.mit.edu/
    """
    form = _choice(form, ("full", "low_re", "mid_re"), "form")
    re = _positive(reynolds, "reynolds")
    pr = _positive(prandtl, "prandtl")
    warnings: list[WarningMessage] = []
    _check_range(re * pr, 0.2, None, "peclet", "Churchill-Bernstein", extrapolate, warnings)
    base = 0.62 * sqrt(re) * pr ** (1.0 / 3.0) / (1.0 + (0.4 / pr) ** (2.0 / 3.0)) ** 0.25
    if form == "full":
        factor, eq = (1.0 + (re / 282000.0) ** 0.625) ** 0.8, "Eq. (7.65)"
    elif form == "low_re":
        _check_range(re, None, 4000.0, "reynolds", "Churchill-Bernstein (7.66)", extrapolate,
                     warnings, hi_open=True)
        factor, eq = 1.0, "Eq. (7.66)"
    else:
        _check_range(re, 40000.0, 400000.0, "reynolds", "Churchill-Bernstein (7.68)",
                     extrapolate, warnings, lo_open=True, hi_open=True)
        factor, eq = 1.0 + sqrt(re / 282000.0), "Eq. (7.68)"
    return _result(0.3 + base * factor, "nusselt", f"Churchill-Bernstein ({form})",
                   f"{_LIENHARD}, {eq}", warnings)


def churchill_chu_vertical_plate_nusselt(rayleigh: float, prandtl: float,
                                         form: Literal["laminar", "all_range"] = "laminar", *,
                                         extrapolate: bool = False) -> CorrelationResult:
    """Natural convection on a vertical plate (Churchill-Chu), Eqs. (8.13a)-(8.13b).

    Laminar (``Ra_L <= 1e9``): ``Nu_L = 0.68 + 0.670 Ra_L^(1/4)
    [1 + (0.492/Pr)^(9/16)]^(-4/9)``; all ``Ra_L`` and ``Pr``:
    ``Nu_L = {0.825 + 0.387 Ra_L^(1/6)/[1 + (0.492/Pr)^(9/16)]^(8/27)}^2``.
    Uniform wall temperature or heat flux; properties at
    ``T_f = (T_inf + T_w)/2``.

    Parameters
    ----------
    rayleigh, prandtl : float
        ``Ra_L`` (:func:`rayleigh_number`) and ``Pr``.
    form : {"laminar", "all_range"}, optional
        Eq. (8.13a) or Eq. (8.13b).
    extrapolate : bool, optional
        Return a warned laminar value for ``Ra_L > 1e9`` instead of raising.

    Returns
    -------
    CorrelationResult
        Average ``Nu_L``.

    References
    ----------
    .. [1] LIENHARD-AHTT-6.00. Lienhard and Lienhard, A Heat Transfer
           Textbook, 6th ed., 2024, Eqs. (8.13a)-(8.13b) and Example 8.1,
           pp. 419-421, https://ahtt.mit.edu/
    """
    form = _choice(form, ("laminar", "all_range"), "form")
    ra = _positive(rayleigh, "rayleigh")
    pr = _positive(prandtl, "prandtl")
    psi = 1.0 + (0.492 / pr) ** (9.0 / 16.0)
    if form == "all_range":
        return _result((0.825 + 0.387 * ra ** (1.0 / 6.0) / psi ** (8.0 / 27.0)) ** 2, "nusselt",
                       "Churchill-Chu all-range", f"{_LIENHARD}, Eq. (8.13b)")
    warnings: list[WarningMessage] = []
    _check_range(ra, None, 1.0e9, "rayleigh", "Churchill-Chu laminar", extrapolate, warnings)
    return _result(0.68 + 0.670 * ra ** 0.25 * psi ** (-4.0 / 9.0), "nusselt",
                   "Churchill-Chu laminar", f"{_LIENHARD}, Eq. (8.13a)", warnings)


# ---------------------------------------------------------------------------
# Boiling, burnout and condensation
# ---------------------------------------------------------------------------


def _rohsenow_scale(viscosity_pa_s: float, latent_heat_j_kg: float, liquid_density_kg_m3: float,
                    vapor_density_kg_m3: float, surface_tension_n_m: float,
                    gravity_m_s2: float) -> float:
    """Return ``mu h_fg sqrt(g (rho_f - rho_g) / sigma)``, the heat-flux scale of Eq. (9.4)."""
    mu = _positive(viscosity_pa_s, "viscosity_pa_s")
    hfg = _positive(latent_heat_j_kg, "latent_heat_j_kg")
    drho = _positive(liquid_density_kg_m3, "liquid_density_kg_m3") - _positive(
        vapor_density_kg_m3, "vapor_density_kg_m3")
    if drho <= 0.0:
        raise DomainError("Liquid density must exceed vapour density.",
                          field="liquid_density_kg_m3")
    sigma = _positive(surface_tension_n_m, "surface_tension_n_m")
    return mu * hfg * sqrt(_positive(gravity_m_s2, "gravity_m_s2") * drho / sigma)


def rohsenow_heat_flux(wall_superheat_k: float, *, specific_heat_j_kg_k: float,
                       latent_heat_j_kg: float, prandtl: float, viscosity_pa_s: float,
                       liquid_density_kg_m3: float, vapor_density_kg_m3: float,
                       surface_tension_n_m: float, surface_constant: float,
                       prandtl_exponent: float,
                       gravity_m_s2: float = STANDARD_GRAVITY_M_S2) -> CorrelationResult:
    """Nucleate pool-boiling heat flux from Rohsenow's Eq. (9.4).

    ``c_p (T_w - T_sat)/(h_fg Pr^s) = C_sf [q/(mu h_fg) sqrt(sigma/(g(rho_f -
    rho_g)))]^0.33`` solved for ``q``; liquid properties at ``T_sat``. Typical
    errors are 100 % in ``q`` and 25 % in ``T_w - T_sat``.

    Parameters
    ----------
    wall_superheat_k : float
        ``T_w - T_sat``.
    specific_heat_j_kg_k, latent_heat_j_kg, prandtl, viscosity_pa_s : float
        Liquid ``c_p``, ``h_fg``, ``Pr`` and ``mu``.
    liquid_density_kg_m3, vapor_density_kg_m3, surface_tension_n_m : float
        ``rho_f``, ``rho_g`` and ``sigma``.
    surface_constant, prandtl_exponent : float
        ``C_sf`` and ``s`` (Lienhard Table 9.2 or measured).
    gravity_m_s2 : float, optional
        Gravitational acceleration.

    Returns
    -------
    CorrelationResult
        Heat flux in W/m^2, with an accuracy info warning.

    References
    ----------
    .. [1] LIENHARD-AHTT-6.00. Lienhard and Lienhard, A Heat Transfer
           Textbook, 6th ed., 2024, Eq. (9.4) and Example 9.2, pp. 482-484,
           https://ahtt.mit.edu/
    """
    dt = _positive(wall_superheat_k, "wall_superheat_k")
    scale = _rohsenow_scale(viscosity_pa_s, latent_heat_j_kg, liquid_density_kg_m3,
                            vapor_density_kg_m3, surface_tension_n_m, gravity_m_s2)
    group = (_positive(specific_heat_j_kg_k, "specific_heat_j_kg_k") * dt
             / (latent_heat_j_kg * _positive(prandtl, "prandtl")
                ** real(prandtl_exponent, "prandtl_exponent")
                * _positive(surface_constant, "surface_constant")))
    return _result(scale * group ** (1.0 / 0.33), "heat_flux_w_m2", "Rohsenow nucleate boiling",
                   f"{_LIENHARD}, Eq. (9.4)",
                   [WarningMessage("low_accuracy_correlation",
                                   "Rohsenow: typical errors 100 % in q, 25 % in superheat.",
                                   "info")])


def rohsenow_wall_superheat(heat_flux_w_m2: float, *, specific_heat_j_kg_k: float,
                            latent_heat_j_kg: float, prandtl: float, viscosity_pa_s: float,
                            liquid_density_kg_m3: float, vapor_density_kg_m3: float,
                            surface_tension_n_m: float, surface_constant: float,
                            prandtl_exponent: float,
                            gravity_m_s2: float = STANDARD_GRAVITY_M_S2) -> CorrelationResult:
    """Wall superheat ``T_w - T_sat`` for a given heat flux, Rohsenow Eq. (9.4).

    Parameters
    ----------
    heat_flux_w_m2 : float
        Wall heat flux ``q``.
    specific_heat_j_kg_k, latent_heat_j_kg, prandtl, viscosity_pa_s : float
        Liquid properties at ``T_sat``.
    liquid_density_kg_m3, vapor_density_kg_m3, surface_tension_n_m : float
        ``rho_f``, ``rho_g``, ``sigma``.
    surface_constant, prandtl_exponent : float
        ``C_sf`` and ``s``.
    gravity_m_s2 : float, optional
        Gravitational acceleration.

    Returns
    -------
    CorrelationResult
        Superheat in K.

    References
    ----------
    .. [1] LIENHARD-AHTT-6.00. Lienhard and Lienhard, A Heat Transfer
           Textbook, 6th ed., 2024, Eq. (9.4), p. 482, https://ahtt.mit.edu/
    """
    q = _positive(heat_flux_w_m2, "heat_flux_w_m2")
    scale = _rohsenow_scale(viscosity_pa_s, latent_heat_j_kg, liquid_density_kg_m3,
                            vapor_density_kg_m3, surface_tension_n_m, gravity_m_s2)
    value = (_positive(surface_constant, "surface_constant") * (q / scale) ** 0.33
             * latent_heat_j_kg * _positive(prandtl, "prandtl")
             ** real(prandtl_exponent, "prandtl_exponent")
             / _positive(specific_heat_j_kg_k, "specific_heat_j_kg_k"))
    return _result(value, "temperature_difference_k", "Rohsenow nucleate boiling (inverse)",
                   f"{_LIENHARD}, Eq. (9.4)")


def taylor_wavelength(liquid_density_kg_m3: float, vapor_density_kg_m3: float,
                      surface_tension_n_m: float,
                      gravity_m_s2: float = STANDARD_GRAVITY_M_S2) -> float:
    """One-dimensional Taylor wavelength ``lambda_d1 = 2 pi sqrt(3) sqrt(sigma/(g drho))``.

    Parameters
    ----------
    liquid_density_kg_m3, vapor_density_kg_m3, surface_tension_n_m : float
        ``rho_f``, ``rho_g`` and ``sigma``.
    gravity_m_s2 : float, optional
        Gravitational acceleration.

    Returns
    -------
    float
        Wavelength in m.

    References
    ----------
    .. [1] LIENHARD-AHTT-6.00. Lienhard and Lienhard, A Heat Transfer
           Textbook, 6th ed., 2024, Eq. (9.6a), p. 486, and p. 488,
           https://ahtt.mit.edu/
    """
    drho = _positive(liquid_density_kg_m3, "liquid_density_kg_m3") - _positive(
        vapor_density_kg_m3, "vapor_density_kg_m3")
    if drho <= 0.0:
        raise DomainError("Liquid density must exceed vapour density.",
                          field="liquid_density_kg_m3")
    return 2.0 * pi * sqrt(3.0) * sqrt(_positive(surface_tension_n_m, "surface_tension_n_m")
                                       / (_positive(gravity_m_s2, "gravity_m_s2") * drho))


def peak_pool_boiling_heat_flux(*, vapor_density_kg_m3: float, latent_heat_j_kg: float,
                                surface_tension_n_m: float, liquid_density_kg_m3: float,
                                form: Literal["flat_plate", "zuber"] = "flat_plate",
                                heater_width_m: float | None = None,
                                gravity_m_s2: float = STANDARD_GRAVITY_M_S2,
                                extrapolate: bool = False) -> CorrelationResult:
    """Peak (critical) pool-boiling heat flux on a large horizontal flat heater.

    ``q_max = C rho_g^(1/2) h_fg [sigma g (rho_f - rho_g)]^(1/4)`` with
    ``C = 0.149`` (``"flat_plate"``, Lienhard-Dhir Eq. 9.11) or ``0.131``
    (``"zuber"``, Zuber-Kutateladze Eq. 9.12, the normalising value
    ``q_maxz``). Eq. (9.11) needs saturated liquid, a well-wetted surface,
    side walls (no entrained side flow) and a heater wider than about
    ``3 lambda_d1``; when ``heater_width_m`` is given this is checked.

    Parameters
    ----------
    vapor_density_kg_m3, latent_heat_j_kg, surface_tension_n_m,
    liquid_density_kg_m3 : float
        Saturation properties.
    form : {"flat_plate", "zuber"}, optional
        Coefficient choice.
    heater_width_m : float, optional
        Heater width or diameter for the ``3 lambda_d1`` check (flat plate).
    gravity_m_s2 : float, optional
        Gravitational acceleration.
    extrapolate : bool, optional
        Return a warned value for a too-small heater instead of raising.

    Returns
    -------
    CorrelationResult
        ``q_max`` in W/m^2.

    References
    ----------
    .. [1] LIENHARD-AHTT-6.00. Lienhard and Lienhard, A Heat Transfer
           Textbook, 6th ed., 2024, Eqs. (9.11)-(9.12), Examples 9.5-9.6,
           pp. 493-496, https://ahtt.mit.edu/
    """
    form = _choice(form, ("flat_plate", "zuber"), "form")
    rho_g = _positive(vapor_density_kg_m3, "vapor_density_kg_m3")
    hfg = _positive(latent_heat_j_kg, "latent_heat_j_kg")
    sigma = _positive(surface_tension_n_m, "surface_tension_n_m")
    rho_f = _positive(liquid_density_kg_m3, "liquid_density_kg_m3")
    g = _positive(gravity_m_s2, "gravity_m_s2")
    if rho_f <= rho_g:
        raise DomainError("Liquid density must exceed vapour density.",
                          field="liquid_density_kg_m3")
    warnings: list[WarningMessage] = []
    if form == "flat_plate" and heater_width_m is not None:
        width = _positive(heater_width_m, "heater_width_m")
        wavelengths = width / taylor_wavelength(rho_f, rho_g, sigma, g)
        _check_range(wavelengths, 3.0, None, "heater_width_m",
                     "Eq. (9.11) large-heater condition (about 3 Taylor wavelengths)",
                     extrapolate, warnings)
    coefficient, eq = (0.149, "Eq. (9.11)") if form == "flat_plate" else (0.131, "Eq. (9.12)")
    value = coefficient * sqrt(rho_g) * hfg * (sigma * g * (rho_f - rho_g)) ** 0.25
    return _result(value, "heat_flux_w_m2", f"peak pool boiling ({form})", f"{_LIENHARD}, {eq}",
                   warnings)


def burnout_margin_screen(heat_flux_w_m2: float, burnout_heat_flux_w_m2: float, *,
                          allowable_fraction: float = 0.8) -> CorrelationResult:
    """Screen a coolant heat flux against a mean burnout heat flux (NASA SP-8087).

    SP-8087 recommends operating heat-flux-limited systems below 80 % of the
    mean burnout heat flux (burnout data have a fractional standard
    deviation of about 0.15; Titan chambers were designed 15 %, newer designs
    at least 18 %, below burnout).

    Parameters
    ----------
    heat_flux_w_m2 : float
        Predicted coolant-side heat flux.
    burnout_heat_flux_w_m2 : float
        Mean burnout (critical) heat flux for the coolant conditions,
        from data or a correlation.
    allowable_fraction : float, optional
        Allowed ``q / q_burnout``; 0.8 per SP-8087.

    Returns
    -------
    CorrelationResult
        ``q / q_burnout``; a ``"burnout_margin_exceeded"`` critical warning
        when above ``allowable_fraction``.

    References
    ----------
    .. [1] NASA-SP-8087-1972. NASA SP-8087, 1972, sec. 2.1.1.5, pp. 15-16,
           and sec. 3.1.1.5.2, p. 61,
           https://ntrs.nasa.gov/citations/19730022965
    """
    q = _nonnegative(heat_flux_w_m2, "heat_flux_w_m2")
    q_bo = _positive(burnout_heat_flux_w_m2, "burnout_heat_flux_w_m2")
    limit = _positive(allowable_fraction, "allowable_fraction")
    if limit > 1.0:
        raise DomainError("allowable_fraction must not exceed 1.", field="allowable_fraction")
    ratio = q / q_bo
    warnings = []
    if ratio > limit:
        warnings.append(WarningMessage(
            "burnout_margin_exceeded",
            f"q/q_burnout = {ratio:.3g} exceeds the allowed fraction {limit:.3g} "
            "(NASA SP-8087 p. 61 recommends below 0.8).", "critical", "heat_flux_w_m2"))
    return _result(ratio, "burnout_heat_flux_ratio", "SP-8087 burnout margin",
                   f"{_SP8087}, p. 61", warnings)


def _velocity_subcooling(velocity_m_s: float, subcooling_k: float) -> float:
    """``V dT_sub`` in the sources' unit, F ft/s."""
    return (_positive(velocity_m_s, "velocity_m_s") / _FOOT_M) * (
        _positive(subcooling_k, "subcooling_k") / _RANKINE_TO_KELVIN)


def ethanol_burnout_heat_flux(velocity_m_s: float, subcooling_k: float, *,
                              pressure_pa: float | None = None,
                              extrapolate: bool = False) -> CorrelationResult:
    """Critical (burnout) heat flux of subcooled ethanol in forced flow (NASA/TM-1998-206612).

    ``CHF(200) = 0.1003 + 0.05264 sqrt(V dT_sub)`` (Eq. 1, Btu/(in^2 s) with
    ``V dT_sub`` in F ft/s), fitted near 200 psia and quoted within +/-20 %;
    with ``pressure_pa`` the factor ``F_P = 1.17 - 8.56e-4 P`` (P in psia,
    Eq. 2) is applied, ``CHF = CHF(200) F_P`` (Eq. 3). ``dT_sub = T_sat -
    T_bulk`` and ``V`` are local values at the transition.

    The tests used SDA-3C alcohol (95 % ethyl, 5 % isopropyl) in 0.125 in.
    OD, 0.016 in. wall 304 stainless tubes; the stated range is checked from
    Table III: ``V dT_sub`` 667-24 230 F ft/s and 144-703 psia. Above 250
    psia, or below it at more than 60 ft/s, the authors observed a small
    50-200 F excursion rather than film boiling at the "CHF" (p. 5); the
    result carries a note in that regime.

    Parameters
    ----------
    velocity_m_s : float
        Local coolant velocity.
    subcooling_k : float
        Local ``T_sat - T_bulk``.
    pressure_pa : float, optional
        Local pressure for the Eq. (2) correction; None returns Eq. (1).
    extrapolate : bool, optional
        Warn instead of raising outside the tested range.

    Returns
    -------
    CorrelationResult
        Burnout heat flux in W/m^2.

    References
    ----------
    .. [1] NASA-TM-1998-206612-MEYER. Meyer, Linne and Rousar, NASA/TM-1998-
           206612, 1998, Eqs. (1)-(3), Table III and Fig. 9, pp. 5-8,
           https://ntrs.nasa.gov/citations/19980017166
    """
    product = _velocity_subcooling(velocity_m_s, subcooling_k)
    warnings: list[WarningMessage] = []
    _check_range(product, 667.0 * (1 - 1e-9), 24230.0 * (1 + 1e-9), "velocity_subcooling",
                 "Meyer-Linne-Rousar ethanol CHF (V dT_sub, F ft/s)", extrapolate, warnings)
    value = 0.1003 + 0.05264 * sqrt(product)
    equation = "Eq. (1)"
    if pressure_pa is not None:
        psia = _positive(pressure_pa, "pressure_pa") / _PSI_PA
        _check_range(psia, 144.0 * (1 - 1e-9), 703.0 * (1 + 1e-9), "pressure_pa", "Meyer-Linne-Rousar pressure factor (psia)",
                     extrapolate, warnings)
        value *= 1.17 - 8.56e-4 * psia
        equation = "Eqs. (1)-(3)"
        velocity_fps = velocity_m_s / _FOOT_M
        if psia > 250.0 or velocity_fps > 60.0:
            warnings.append(WarningMessage(
                "mild_transition_regime",
                "Above 250 psia or 60 ft/s the tests showed a small 50-200 F excursion at this "
                "flux rather than film boiling (NASA/TM-1998-206612 p. 5).", "info", "pressure_pa"))
    return _result(value * _BTU_IN2_S, "heat_flux_w_m2", "ethanol subcooled burnout",
                   f"{_MEYER}, {equation}", warnings)


def mmh_burnout_heat_flux(velocity_m_s: float, subcooling_k: float, *,
                          extrapolate: bool = False) -> CorrelationResult:
    """Burnout heat flux of subcooled monomethylhydrazine in a heated tube (NASA CR-134282).

    ``(Q/A)_BO = 0.025 (V dT_sub)^0.633`` (Eq. A-1; Btu/(in^2 s), V in ft/s,
    ``dT_sub`` in F), from uniformly electrically heated circular tubes; the
    source calls it somewhat conservative (zero flux at zero velocity or
    subcooling). The data of Fig. A1 cover 10-56 ft/s and 31-226 F of
    subcooling; both are checked. At ``V dT_sub = 11 000`` the equation gives
    9.0 Btu/(in^2 s) while the drawn design curve of Fig. A1 ends near 8.2,
    so the printed equation and curve differ at the top of the range.

    Parameters
    ----------
    velocity_m_s : float
        Local coolant velocity.
    subcooling_k : float
        ``T_sat - T_local``.
    extrapolate : bool, optional
        Warn instead of raising outside the tested range.

    Returns
    -------
    CorrelationResult
        Burnout heat flux in W/m^2.

    References
    ----------
    .. [1] NASA-CR-134282-TOBIN-1974. R. D. Tobin, NASA CR-134282, 1974,
           Eq. (A-1), p. A1; Fig. A1, p. A2; p. A3,
           https://ntrs.nasa.gov/citations/19740018134
    """
    # Ranges are the tested ones; the 1e-9 slack absorbs unit round-off at the ends.
    velocity_fps = _positive(velocity_m_s, "velocity_m_s") / _FOOT_M
    subcooling_f = _positive(subcooling_k, "subcooling_k") / _RANKINE_TO_KELVIN
    warnings: list[WarningMessage] = []
    _check_range(velocity_fps, 10.0 * (1 - 1e-9), 56.0 * (1 + 1e-9), "velocity_m_s", "Tobin MMH burnout (ft/s)", extrapolate, warnings)
    _check_range(subcooling_f, 31.0 * (1 - 1e-9), 226.0 * (1 + 1e-9), "subcooling_k", "Tobin MMH burnout (F)", extrapolate, warnings)
    value = 0.025 * (velocity_fps * subcooling_f) ** 0.633
    return _result(value * _BTU_IN2_S, "heat_flux_w_m2", "MMH subcooled burnout", f"{_TOBIN}, Eq. (A-1)",
                   warnings)


def corrected_latent_heat(latent_heat_j_kg: float, liquid_specific_heat_j_kg_k: float,
                          saturation_minus_wall_k: float, liquid_prandtl: float, *,
                          extrapolate: bool = False) -> CorrelationResult:
    """Sadasivan-Lienhard corrected latent heat, Lienhard Eq. (8.61).

    ``h'_fg = h_fg [1 + (0.683 - 0.228/Pr) Ja]`` with
    ``Ja = c_p (T_sat - T_w)/h_fg``, for ``Pr >= 0.6``.

    Parameters
    ----------
    latent_heat_j_kg, liquid_specific_heat_j_kg_k : float
        ``h_fg`` and liquid ``c_p``.
    saturation_minus_wall_k : float
        ``T_sat - T_w`` (nonnegative).
    liquid_prandtl : float
        Liquid ``Pr``.
    extrapolate : bool, optional
        Return a warned value for ``Pr < 0.6`` instead of raising.

    Returns
    -------
    CorrelationResult
        ``h'_fg`` in J/kg.

    References
    ----------
    .. [1] LIENHARD-AHTT-6.00. Lienhard and Lienhard, A Heat Transfer
           Textbook, 6th ed., 2024, Eq. (8.61) and Example 8.6, pp. 446-448,
           https://ahtt.mit.edu/
    """
    hfg = _positive(latent_heat_j_kg, "latent_heat_j_kg")
    pr = _positive(liquid_prandtl, "liquid_prandtl")
    warnings: list[WarningMessage] = []
    _check_range(pr, 0.6, None, "liquid_prandtl", "Sadasivan-Lienhard h'_fg", extrapolate,
                 warnings)
    ja = (_positive(liquid_specific_heat_j_kg_k, "liquid_specific_heat_j_kg_k")
          * _nonnegative(saturation_minus_wall_k, "saturation_minus_wall_k") / hfg)
    return _result(hfg * (1.0 + (0.683 - 0.228 / pr) * ja), "latent_heat_j_kg",
                   "corrected latent heat", f"{_LIENHARD}, Eq. (8.61)", warnings)


def film_condensation_vertical_plate_nusselt(*, length_m: float, liquid_density_kg_m3: float,
                                             vapor_density_kg_m3: float,
                                             corrected_latent_heat_j_kg: float,
                                             liquid_kinematic_viscosity_m2_s: float,
                                             liquid_conductivity_w_m_k: float,
                                             saturation_minus_wall_k: float,
                                             gravity_m_s2: float = STANDARD_GRAVITY_M_S2
                                             ) -> CorrelationResult:
    """Average Nusselt number of laminar film condensation on a vertical plate.

    Nusselt's result, Lienhard Eq. (8.62b):
    ``Nu_L = 0.9428 [g (rho_f - rho_g) h'_fg L^3 / (nu k (T_sat - T_w))]^(1/4)``
    for an isothermal wall and a laminar film; liquid properties at the mean
    film temperature.

    Parameters
    ----------
    length_m : float
        Plate height ``L``.
    liquid_density_kg_m3, vapor_density_kg_m3 : float
        ``rho_f`` and ``rho_g``.
    corrected_latent_heat_j_kg : float
        ``h'_fg`` (:func:`corrected_latent_heat`).
    liquid_kinematic_viscosity_m2_s, liquid_conductivity_w_m_k : float
        Liquid ``nu`` and ``k``.
    saturation_minus_wall_k : float
        ``T_sat - T_w``, positive.
    gravity_m_s2 : float, optional
        Gravitational acceleration.

    Returns
    -------
    CorrelationResult
        ``Nu_L = h L / k`` with an info warning that film laminarity is not
        checked.

    References
    ----------
    .. [1] LIENHARD-AHTT-6.00. Lienhard and Lienhard, A Heat Transfer
           Textbook, 6th ed., 2024, Eqs. (8.56)-(8.62b) and Example 8.6,
           pp. 444-449, https://ahtt.mit.edu/
    """
    drho = _positive(liquid_density_kg_m3, "liquid_density_kg_m3") - _positive(
        vapor_density_kg_m3, "vapor_density_kg_m3")
    if drho <= 0.0:
        raise DomainError("Liquid density must exceed vapour density.",
                          field="liquid_density_kg_m3")
    group = (_positive(gravity_m_s2, "gravity_m_s2") * drho
             * _positive(corrected_latent_heat_j_kg, "corrected_latent_heat_j_kg")
             * _positive(length_m, "length_m") ** 3
             / (_positive(liquid_kinematic_viscosity_m2_s, "liquid_kinematic_viscosity_m2_s")
                * _positive(liquid_conductivity_w_m_k, "liquid_conductivity_w_m_k")
                * _positive(saturation_minus_wall_k, "saturation_minus_wall_k")))
    return _result(0.9428 * group ** 0.25, "nusselt", "Nusselt laminar film condensation",
                   f"{_LIENHARD}, Eq. (8.62b)",
                   [WarningMessage("laminar_film_assumed",
                                   "Laminar, smooth film assumed; transition is not checked.",
                                   "info")])


# ---------------------------------------------------------------------------
# Film cooling (Hatch-Papell)
# ---------------------------------------------------------------------------

_HP_SLOT = (0.063 * _INCH_M, 0.50 * _INCH_M)
_HP_GAS_VELOCITY = (104.0 * _FOOT_M, 1040.0 * _FOOT_M)
_HP_VELOCITY_RATIO = (0.45, 33.3)
_HP_GAS_T = (502.0 * _RANKINE_TO_KELVIN, 1965.0 * _RANKINE_TO_KELVIN)
_HP_COOLANT_T = (530.0 * _RANKINE_TO_KELVIN, 980.0 * _RANKINE_TO_KELVIN)


def hatch_papell_velocity_ratio_function(gas_to_coolant_velocity_ratio: float) -> float:
    """Hatch-Papell velocity-ratio function ``f(V_g/V_c)``, Eqs. (10)-(11).

    ``1 + 0.4 arctan(V_g/V_c - 1)`` (radians) for ``V_g/V_c >= 1`` and
    ``(V_c/V_g)^(1.5 (V_c/V_g - 1))`` for ``V_g/V_c <= 1``.

    Parameters
    ----------
    gas_to_coolant_velocity_ratio : float
        ``V_g / V_c``.

    Returns
    -------
    float
        ``f``; equal to 1 at a velocity ratio of 1.

    References
    ----------
    .. [1] NASA-TN-D-130-HATCH-PAPELL-1959. Hatch and Papell, NASA TN D-130,
           1959, Eqs. (10)-(11), p. 5,
           https://ntrs.nasa.gov/citations/19890068390
    """
    r = _positive(gas_to_coolant_velocity_ratio, "gas_to_coolant_velocity_ratio")
    if r >= 1.0:
        return 1.0 + 0.4 * atan(r - 1.0)
    inverse = 1.0 / r
    return inverse ** (1.5 * (inverse - 1.0))


def _hatch_papell_checks(slot_height_m: float, gas_velocity_m_s: float,
                         coolant_velocity_m_s: float, gas_temperature_k: float | None,
                         coolant_temperature_k: float | None, extrapolate: bool,
                         warnings: list[WarningMessage]) -> float:
    """Apply the Hatch-Papell data ranges (TN D-130 summary) and return ``V_g/V_c``."""
    name = "Hatch-Papell data range"
    s = _positive(slot_height_m, "slot_height_m")
    vg = _positive(gas_velocity_m_s, "gas_velocity_m_s")
    ratio = vg / _positive(coolant_velocity_m_s, "coolant_velocity_m_s")
    _check_range(s, *_HP_SLOT, "slot_height_m", name, extrapolate, warnings)
    _check_range(vg, *_HP_GAS_VELOCITY, "gas_velocity_m_s", name, extrapolate, warnings)
    _check_range(ratio, *_HP_VELOCITY_RATIO, "velocity_ratio", name, extrapolate, warnings)
    if gas_temperature_k is not None:
        _check_range(_positive(gas_temperature_k, "gas_temperature_k"), *_HP_GAS_T,
                     "gas_temperature_k", name, extrapolate, warnings)
    if coolant_temperature_k is not None:
        _check_range(_positive(coolant_temperature_k, "coolant_temperature_k"), *_HP_COOLANT_T,
                     "coolant_temperature_k", name, extrapolate, warnings)
    return ratio


def hatch_papell_effectiveness(*, heat_transfer_coefficient_w_m2_k: float, slot_length_m: float,
                               distance_m: float, coolant_mass_flow_kg_s: float,
                               coolant_specific_heat_j_kg_k: float, slot_height_m: float,
                               gas_velocity_m_s: float, coolant_velocity_m_s: float,
                               coolant_thermal_diffusivity_m2_s: float,
                               gas_temperature_k: float | None = None,
                               coolant_temperature_k: float | None = None,
                               extrapolate: bool = False) -> CorrelationResult:
    """Adiabatic-wall film-cooling effectiveness for tangential slot injection of a gas.

    Hatch-Papell Eq. (12): ``ln eta = -[h L x/(w c_p)_c - 0.04]
    (S V_g/alpha_c)^0.125 f(V_g/V_c)`` for ``h L x/(w c_p)_c >= 0.04``, and
    ``eta = 1`` below (Eq. 12a), with
    ``eta = (T_ad - t_w)/(T_ad - t_c)``. ``h`` is the gas-side coefficient
    without film (TN D-130 uses ``0.0265 (k/D_h) Re^0.8 Pr^0.3`` with film
    properties at ``(t_g + t_c)/2``), ``L`` the slot length (cooled
    perimeter), ``alpha_c`` the coolant diffusivity at slot-exit static
    temperature. Wall temperature predicted within about 5 % for
    ``0.2 <= eta <= 1``; values below about 0.3 are pessimistic. Data range:
    slot height 0.063-0.50 in, main-gas velocity 104-1040 ft/s, velocity
    ratio 0.45-33.3, main-gas temperature 502-1965 degR, coolant
    temperature 530-980 degR (checked when given).

    Parameters
    ----------
    heat_transfer_coefficient_w_m2_k : float
        Gas-to-wall coefficient ``h``.
    slot_length_m : float
        Slot length ``L`` (perpendicular to the flow).
    distance_m : float
        Distance ``x`` downstream of the slot.
    coolant_mass_flow_kg_s, coolant_specific_heat_j_kg_k : float
        ``w_c`` and ``c_p,c``.
    slot_height_m : float
        Slot height ``S``.
    gas_velocity_m_s, coolant_velocity_m_s : float
        ``V_g`` and ``V_c``.
    coolant_thermal_diffusivity_m2_s : float
        ``alpha_c``.
    gas_temperature_k, coolant_temperature_k : float, optional
        Main-gas and coolant slot-exit static temperatures for the range check.
    extrapolate : bool, optional
        Return a warned value outside the data range instead of raising.

    Returns
    -------
    CorrelationResult
        Effectiveness ``eta``.

    References
    ----------
    .. [1] NASA-TN-D-130-HATCH-PAPELL-1959. Hatch and Papell, NASA TN D-130,
           1959, Summary p. 1, assumption (6) p. 2, Eqs. (10)-(12a) pp. 5-6,
           accuracy pp. 6-8 and 10, symbols pp. 11-12,
           https://ntrs.nasa.gov/citations/19890068390
    .. [2] NASA-SP-8087-1972. NASA SP-8087, 1972, p. 15,
           https://ntrs.nasa.gov/citations/19730022965
    """
    warnings: list[WarningMessage] = []
    ratio = _hatch_papell_checks(slot_height_m, gas_velocity_m_s, coolant_velocity_m_s,
                                 gas_temperature_k, coolant_temperature_k, extrapolate, warnings)
    group = (_positive(heat_transfer_coefficient_w_m2_k, "heat_transfer_coefficient_w_m2_k")
             * _positive(slot_length_m, "slot_length_m") * _nonnegative(distance_m, "distance_m")
             / (_positive(coolant_mass_flow_kg_s, "coolant_mass_flow_kg_s")
                * _positive(coolant_specific_heat_j_kg_k, "coolant_specific_heat_j_kg_k")))
    if group <= 0.04:
        return _result(1.0, "effectiveness", "Hatch-Papell", f"{_HATCH}, Eq. (12a)", warnings)
    peclet_slot = (slot_height_m * gas_velocity_m_s
                   / _positive(coolant_thermal_diffusivity_m2_s,
                               "coolant_thermal_diffusivity_m2_s"))
    eta = exp(-(group - 0.04) * peclet_slot ** 0.125 * hatch_papell_velocity_ratio_function(ratio))
    if eta < 0.2:
        warnings.append(WarningMessage(
            "below_validated_effectiveness",
            f"eta = {eta:.3g} is below the 0.2-1.0 range of the 5 % wall-temperature "
            "agreement; Eq. (12) is pessimistic there.", "warning"))
    return _result(eta, "effectiveness", "Hatch-Papell", f"{_HATCH}, Eq. (12)", warnings)


def hatch_papell_required_coolant_flow(*, effectiveness: float,
                                       heat_transfer_coefficient_w_m2_k: float,
                                       slot_length_m: float, distance_m: float,
                                       coolant_specific_heat_j_kg_k: float, slot_height_m: float,
                                       gas_velocity_m_s: float, coolant_velocity_m_s: float,
                                       coolant_thermal_diffusivity_m2_s: float,
                                       gas_temperature_k: float | None = None,
                                       coolant_temperature_k: float | None = None,
                                       extrapolate: bool = False) -> CorrelationResult:
    """Coolant mass flow giving a target effectiveness, Hatch-Papell Eq. (13).

    ``w_c = h L x / c_p,c * 1/[-ln(eta) / ((S V_g/alpha_c)^0.125 f(V_g/V_c)) + 0.04]``.
    The source iterates on slot height through continuity; this function
    evaluates one step with the given geometry and velocities.

    Parameters
    ----------
    effectiveness : float
        Target ``eta``, ``0 < eta < 1``.
    heat_transfer_coefficient_w_m2_k, slot_length_m, distance_m : float
        ``h``, ``L`` and ``x``.
    coolant_specific_heat_j_kg_k, slot_height_m : float
        ``c_p,c`` and ``S``.
    gas_velocity_m_s, coolant_velocity_m_s, coolant_thermal_diffusivity_m2_s : float
        ``V_g``, ``V_c`` and ``alpha_c``.
    gas_temperature_k, coolant_temperature_k : float, optional
        For the data-range check.
    extrapolate : bool, optional
        Return a warned value outside the data range instead of raising.

    Returns
    -------
    CorrelationResult
        ``w_c`` in kg/s.

    References
    ----------
    .. [1] NASA-TN-D-130-HATCH-PAPELL-1959. Hatch and Papell, NASA TN D-130,
           1959, Eq. (13) and example, pp. 8-10,
           https://ntrs.nasa.gov/citations/19890068390
    """
    eta = _positive(effectiveness, "effectiveness")
    if eta >= 1.0:
        raise DomainError("Target effectiveness must be below 1.", field="effectiveness")
    warnings: list[WarningMessage] = []
    ratio = _hatch_papell_checks(slot_height_m, gas_velocity_m_s, coolant_velocity_m_s,
                                 gas_temperature_k, coolant_temperature_k, extrapolate, warnings)
    peclet_slot = (slot_height_m * gas_velocity_m_s
                   / _positive(coolant_thermal_diffusivity_m2_s,
                               "coolant_thermal_diffusivity_m2_s"))
    denominator = -log(eta) / (peclet_slot ** 0.125
                               * hatch_papell_velocity_ratio_function(ratio)) + 0.04
    value = (_positive(heat_transfer_coefficient_w_m2_k, "heat_transfer_coefficient_w_m2_k")
             * _positive(slot_length_m, "slot_length_m") * _positive(distance_m, "distance_m")
             / _positive(coolant_specific_heat_j_kg_k, "coolant_specific_heat_j_kg_k")
             / denominator)
    return _result(value, "mass_flow_kg_s", "Hatch-Papell required coolant flow",
                   f"{_HATCH}, Eq. (13)", warnings)


def film_cooled_wall_temperature(effectiveness: float, adiabatic_wall_temperature_k: float,
                                 coolant_temperature_k: float) -> float:
    """Wall temperature from the Hatch-Papell effectiveness definition.

    ``t_w = T_ad - eta (T_ad - t_c)``.

    Parameters
    ----------
    effectiveness : float
        ``eta`` in ``[0, 1]``.
    adiabatic_wall_temperature_k : float
        Main-stream recovery temperature ``T_ad``.
    coolant_temperature_k : float
        Coolant slot-exit static temperature ``t_c``.

    Returns
    -------
    float
        Adiabatic film-cooled wall temperature, K.

    References
    ----------
    .. [1] NASA-TN-D-130-HATCH-PAPELL-1959. Hatch and Papell, NASA TN D-130,
           1959, symbols, p. 12, https://ntrs.nasa.gov/citations/19890068390
    """
    eta = _nonnegative(effectiveness, "effectiveness")
    if eta > 1.0:
        raise DomainError("Effectiveness must not exceed 1.", field="effectiveness")
    t_ad = _positive(adiabatic_wall_temperature_k, "adiabatic_wall_temperature_k")
    return t_ad - eta * (t_ad - _positive(coolant_temperature_k, "coolant_temperature_k"))
