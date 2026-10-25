"""Transpiration (porous-wall) cooling of rocket chamber and nozzle walls.

A coolant forced through a porous wall does two things: it absorbs heat as
it passes through the wall, and it thickens the gas-side boundary layer as it
leaves the surface, which lowers the gas-side heat-transfer coefficient. This
module provides

* the film-theory correction of a heat-transfer coefficient for blowing,
  ``h/h* = phi/(e^phi - 1)``, and its resistance-factor form
  ``phi = ln(1 + R)`` (Mickley et al., NACA TN 3208) [1]_;
* the coolant mass flux that holds a wall at a target temperature, and the
  wall temperature reached for a given mass flux, from that correction and
  the coolant energy balance with an optional porous-wall efficiency (Pratt
  & Whitney, NASA CR-107268) [2]_ and an incident radiative heat flux
  (Eckert & Livingood, NACA Report 1182) [3]_;
* the Eckert-Livingood turbulent flat-plate temperature ratios, with the
  Rannie-Friedman sublayer factor ``r = 2.11 Re^-0.1`` or the film-theory
  limit [3]_;
* Rannie's turbulent pipe relation as printed by Kelley & L'Ecuyer [4]_;
* the two-temperature (solid and coolant) solution through a porous wall
  and its heat-exchange efficiency [2]_.

Conventions: SI units; temperatures in kelvin. The gas-side driving
temperature is the adiabatic-wall (recovery) temperature; for low-speed
flow it equals the free-stream temperature. Blowing rates and rate factors
are positive for injection into the gas (``v_w > 0``) and negative for
suction.

Limits: no foreign-gas (molecular-weight) diffusion effect beyond the
coolant heat capacity entering the rate factor; Rubesin & Pappas (NACA TN
4149) show hydrogen to be about five times as effective as air on a mass
basis, which film theory with the coolant's heat capacity only partly
captures. Mickley et al. found film theory to predict turbulent heat transfer
with roughly uniform blowing within the precision of their data but
recommended no correction for other blowing distributions. No transpiration
model here is validated against a rocket firing.

References
----------
.. [1] NACA-TN-3208-MICKLEY-1954. H. S. Mickley, R. C. Ross, A. L. Squyers
       and W. E. Stewart, NACA TN 3208, 1954, Eqs. (21)-(30), pp. 14-18,
       Table VII, pp. 78-89, https://ntrs.nasa.gov/citations/19930084056
.. [2] NASA-CR-107268-PWA-1969. Pratt & Whitney Aircraft, PWA FR-3390 /
       NASA CR-107268, 1969, Eqs. (III-9)-(III-29), pp. III-11 to III-16,
       https://ntrs.nasa.gov/citations/19700004126
.. [3] NACA-REPORT-1182-ECKERT-LIVINGOOD-1954. E. R. G. Eckert and
       J. N. B. Livingood, NACA Report 1182, 1954, Eqs. (20)-(31),
       pp. 8-11, https://ntrs.nasa.gov/citations/19930092205
.. [4] NASA-CR-76367-KELLEY-LECUYER-1966. J. B. Kelley and M. R. L'Ecuyer,
       Purdue JPC TM-66-5 / NASA CR-76367, 1966, Eq. (2.17), p. 17, and
       Appendix B, pp. 147-152, https://ntrs.nasa.gov/citations/19660021566
"""

from __future__ import annotations

from dataclasses import dataclass
from math import exp, expm1, isfinite, log1p, sqrt
from typing import Literal

from rocket_propulsion.core.errors import DomainError
from rocket_propulsion.core.metadata import WarningMessage
from rocket_propulsion.core.validation import real
from rocket_propulsion.heat_transfer.convection import CorrelationResult

__all__ = [
    "PorousWallSolution",
    "TranspirationState",
    "eckert_livingood_temperature_ratio",
    "film_theory_correction",
    "film_theory_rate_factor",
    "porous_wall_solution",
    "rannie_pipe_blowing_rate",
    "rannie_pipe_temperature_ratio",
    "transpiration_coolant_flux",
    "transpiration_wall_temperature",
]

_MICKLEY = "Mickley et al., NACA TN 3208 (1954)"
_PWA = "Pratt & Whitney, NASA CR-107268 (1969)"
_ECKERT = "Eckert & Livingood, NACA Report 1182 (1954)"
_KELLEY = "Kelley & L'Ecuyer, NASA CR-76367 (1966)"

# Below this |phi| the series 1 - phi/2 + phi^2/12 - phi^4/720 is used; its
# truncation error (phi^6/30240) is below 1e-16 there.
_SERIES_LIMIT = 1.0e-3


def _positive(value: object, label: str) -> float:
    number = real(value, label)
    if number <= 0.0:
        raise DomainError(f"{label} must be positive.", field=label)
    return number


def _nonnegative(value: object, label: str) -> float:
    number = real(value, label)
    if number < 0.0:
        raise DomainError(f"{label} must be nonnegative.", field=label)
    return number


def _blowing_ratio(phi: float) -> float:
    """``phi/(e^phi - 1)``, continuous through ``phi = 0``."""
    if abs(phi) < _SERIES_LIMIT:
        return 1.0 - phi / 2.0 + phi * phi / 12.0 - phi ** 4 / 720.0
    if phi > 700.0:
        return 0.0
    return phi / expm1(phi)


# ---------------------------------------------------------------------------
# Film theory
# ---------------------------------------------------------------------------


def film_theory_correction(rate_factor: float) -> CorrelationResult:
    """Film-theory ratio of the blown to the unblown transfer coefficient.

    ``theta = phi/(e^phi - 1)`` with ``phi = (rho v)_w c_p,c / h*`` for heat
    (``phi_H``) or ``2 (rho v)_w/(rho_1 u_1 c_f*)`` for friction (``phi_F``),
    where the asterisk marks the value without mass transfer and ``c_p,c``
    is the heat capacity of the injected coolant.

    Parameters
    ----------
    rate_factor : float
        ``phi``; positive for blowing, negative for suction.

    Returns
    -------
    CorrelationResult
        ``theta = h/h*`` (or ``c_f/c_f*``). ``theta -> 1`` as ``phi -> 0``,
        ``theta -> 0`` as ``phi -> +inf`` and ``theta -> -phi`` as
        ``phi -> -inf`` (TN 3208 footnote 2).

    References
    ----------
    .. [1] NACA-TN-3208-MICKLEY-1954. Mickley et al., NACA TN 3208, 1954,
           Eqs. (23b), (24) and (27), pp. 15-16, footnote 2, p. 17,
           https://ntrs.nasa.gov/citations/19930084056
    """
    phi = real(rate_factor, "rate_factor")
    return CorrelationResult(_blowing_ratio(phi), "transfer_coefficient_ratio",
                             "film theory (blowing)", f"{_MICKLEY}, Eq. (27)")


def film_theory_rate_factor(resistance_factor: float) -> float:
    """Rate factor ``phi = ln(R + 1)`` from the resistance factor ``R = phi/theta``.

    For heat transfer ``R_H = (rho v)_w c_p,c / h = (T_0 - T_1)/(T_s - T_0)``
    (wall 0, main stream 1, coolant supply ``s``), so ``R`` follows from the
    temperatures alone and the coolant flow needs no iteration.

    Parameters
    ----------
    resistance_factor : float
        ``R > -1``.

    Returns
    -------
    float
        ``phi``.

    References
    ----------
    .. [1] NACA-TN-3208-MICKLEY-1954. Mickley et al., NACA TN 3208, 1954,
           Eqs. (28)-(30), pp. 17-18,
           https://ntrs.nasa.gov/citations/19930084056
    """
    r = real(resistance_factor, "resistance_factor")
    if r <= -1.0:
        raise DomainError("resistance_factor must exceed -1.", field="resistance_factor")
    return log1p(r)


@dataclass(frozen=True, slots=True)
class TranspirationState:
    """Transpiration-cooled wall station solved with film theory.

    Attributes
    ----------
    coolant_mass_flux_kg_m2_s : float
        Coolant mass flux through the wall, ``(rho v)_w``.
    wall_temperature_k : float
        Gas-side wall temperature.
    rate_factor : float
        ``phi = (rho v)_w c_p,c / h*``.
    blowing_correction : float
        ``h/h* = phi/(e^phi - 1)``.
    blown_heat_transfer_coefficient_w_m2_k : float
        ``h``, the gas-side coefficient with blowing.
    convective_heat_flux_w_m2 : float
        ``h (T_aw - T_w)`` reaching the wall.
    absorbed_heat_flux_w_m2 : float
        ``(rho v)_w c_p,c eta (T_w - T_c)`` taken up by the coolant.
    unblown_heat_flux_w_m2 : float
        ``h* (T_aw - T_w)``, the flux the same wall would see without blowing.
    warnings : tuple of WarningMessage
        Applicability notes.

    References
    ----------
    .. [1] NACA-TN-3208-MICKLEY-1954. Mickley et al., NACA TN 3208, 1954,
           https://ntrs.nasa.gov/citations/19930084056
    """

    coolant_mass_flux_kg_m2_s: float
    wall_temperature_k: float
    rate_factor: float
    blowing_correction: float
    blown_heat_transfer_coefficient_w_m2_k: float
    convective_heat_flux_w_m2: float
    absorbed_heat_flux_w_m2: float
    unblown_heat_flux_w_m2: float
    warnings: tuple[WarningMessage, ...] = ()


def _check_temperatures(recovery: float, wall: float, coolant: float) -> None:
    if not coolant < wall < recovery:
        raise DomainError(
            "Transpiration cooling needs coolant temperature < wall temperature < "
            "recovery temperature.", field="wall_temperature_k",
            details={"coolant_temperature_k": coolant, "wall_temperature_k": wall,
                     "recovery_temperature_k": recovery})


def _state(flux: float, wall: float, h_star: float, cp: float, recovery: float,
           coolant: float, efficiency: float, radiation: float) -> TranspirationState:
    phi = flux * cp / h_star
    theta = _blowing_ratio(phi)
    h = h_star * theta
    warnings = []
    if radiation > 0.0:
        warnings.append(WarningMessage(
            "radiation_not_blocked",
            "The incident radiative flux is added to the wall load unchanged; film "
            "theory does not reduce radiation.", "info", "radiative_heat_flux_w_m2"))
    return TranspirationState(
        flux, wall, phi, theta, h, h * (recovery - wall),
        flux * cp * efficiency * (wall - coolant), h_star * (recovery - wall),
        tuple(warnings))


def transpiration_coolant_flux(*, unblown_heat_transfer_coefficient_w_m2_k: float,
                               recovery_temperature_k: float, wall_temperature_k: float,
                               coolant_temperature_k: float,
                               coolant_specific_heat_j_kg_k: float,
                               wall_efficiency: float = 1.0,
                               radiative_heat_flux_w_m2: float = 0.0) -> TranspirationState:
    """Coolant mass flux that holds a porous wall at ``wall_temperature_k``.

    Steady energy balance on the wall (Eckert & Livingood Eq. 20, Pratt &
    Whitney Eq. III-9): the coolant, entering at ``T_c`` and leaving the hot
    face at ``T_c + eta (T_w - T_c)``, absorbs the blown convective flux and
    any radiation,

        ``(rho v)_w c_p,c eta (T_w - T_c) = h* theta(phi) (T_aw - T_w) + q_r``,

    with ``phi = (rho v)_w c_p,c / h*`` and ``theta`` from film theory.
    Without radiation the solution is closed form,
    ``phi = ln(1 + R)`` with ``R = (T_aw - T_w)/(eta (T_w - T_c))``
    (TN 3208 Eq. 30); with radiation the balance is solved by bisection.

    Parameters
    ----------
    unblown_heat_transfer_coefficient_w_m2_k : float
        ``h*``, the gas-side coefficient without blowing (e.g. from Bartz).
    recovery_temperature_k : float
        Adiabatic-wall (recovery) temperature of the gas.
    wall_temperature_k : float
        Target gas-side wall temperature.
    coolant_temperature_k : float
        Coolant supply temperature behind the wall.
    coolant_specific_heat_j_kg_k : float
        Mean coolant heat capacity between ``T_c`` and ``T_w``.
    wall_efficiency : float, optional
        ``eta = (T_co - T_c)/(T_w - T_c)`` in (0, 1]; 1 when the coolant leaves
        at the wall temperature (see :func:`porous_wall_solution`).
    radiative_heat_flux_w_m2 : float, optional
        Radiation absorbed by the wall, not reduced by blowing.

    Returns
    -------
    TranspirationState

    References
    ----------
    .. [1] NACA-TN-3208-MICKLEY-1954. Mickley et al., NACA TN 3208, 1954,
           Eqs. (27)-(30), pp. 16-18,
           https://ntrs.nasa.gov/citations/19930084056
    .. [2] NASA-CR-107268-PWA-1969. Pratt & Whitney, NASA CR-107268, 1969,
           Eq. (III-9), p. III-12,
           https://ntrs.nasa.gov/citations/19700004126
    .. [3] NACA-REPORT-1182-ECKERT-LIVINGOOD-1954. Eckert & Livingood, NACA
           Report 1182, 1954, Eq. (20), p. 8,
           https://ntrs.nasa.gov/citations/19930092205
    """
    h_star = _positive(unblown_heat_transfer_coefficient_w_m2_k,
                       "unblown_heat_transfer_coefficient_w_m2_k")
    recovery = _positive(recovery_temperature_k, "recovery_temperature_k")
    wall = _positive(wall_temperature_k, "wall_temperature_k")
    coolant = _positive(coolant_temperature_k, "coolant_temperature_k")
    cp = _positive(coolant_specific_heat_j_kg_k, "coolant_specific_heat_j_kg_k")
    eta = _positive(wall_efficiency, "wall_efficiency")
    if eta > 1.0:
        raise DomainError("wall_efficiency must not exceed 1.", field="wall_efficiency")
    q_r = _nonnegative(radiative_heat_flux_w_m2, "radiative_heat_flux_w_m2")
    _check_temperatures(recovery, wall, coolant)
    if q_r == 0.0:
        phi = log1p((recovery - wall) / (eta * (wall - coolant)))
    else:
        def residual(p: float) -> float:
            return h_star * p * eta * (wall - coolant) - h_star * _blowing_ratio(p) * (
                recovery - wall) - q_r
        lo, hi = 0.0, 1.0
        while residual(hi) < 0.0:
            hi *= 2.0
            if hi > 1.0e4:
                raise DomainError("No coolant flow balances the radiative load.",
                                  field="radiative_heat_flux_w_m2")
        for _ in range(200):
            mid = 0.5 * (lo + hi)
            if residual(mid) < 0.0:
                lo = mid
            else:
                hi = mid
            if hi - lo <= 1e-15 * max(1.0, hi):
                break
        phi = 0.5 * (lo + hi)
    return _state(phi * h_star / cp, wall, h_star, cp, recovery, coolant, eta, q_r)


def transpiration_wall_temperature(*, coolant_mass_flux_kg_m2_s: float,
                                   unblown_heat_transfer_coefficient_w_m2_k: float,
                                   recovery_temperature_k: float,
                                   coolant_temperature_k: float,
                                   coolant_specific_heat_j_kg_k: float,
                                   wall_efficiency: float = 1.0,
                                   radiative_heat_flux_w_m2: float = 0.0) -> TranspirationState:
    """Wall temperature of a porous wall for a given coolant mass flux.

    Inverse of :func:`transpiration_coolant_flux`. Without radiation,
    ``T_w = T_c + (T_aw - T_c)/(1 + eta (e^phi - 1))``; for ``eta = 1`` this
    is ``(T_w - T_c)/(T_aw - T_c) = e^-phi``, the film-theory limit of
    Eckert & Livingood Eq. (30). With radiation,
    ``T_w (h* phi eta + h* theta) = h* phi eta T_c + h* theta T_aw + q_r``
    is linear in ``T_w`` because ``phi`` depends only on the flux.

    Parameters
    ----------
    coolant_mass_flux_kg_m2_s : float
        ``(rho v)_w``, nonnegative.
    unblown_heat_transfer_coefficient_w_m2_k, recovery_temperature_k,
    coolant_temperature_k, coolant_specific_heat_j_kg_k, wall_efficiency,
    radiative_heat_flux_w_m2 : float
        As in :func:`transpiration_coolant_flux`.

    Returns
    -------
    TranspirationState

    References
    ----------
    .. [1] NACA-TN-3208-MICKLEY-1954. Mickley et al., NACA TN 3208, 1954,
           Eqs. (27)-(30), pp. 16-18,
           https://ntrs.nasa.gov/citations/19930084056
    .. [2] NACA-REPORT-1182-ECKERT-LIVINGOOD-1954. Eckert & Livingood, NACA
           Report 1182, 1954, Eqs. (20) and (30), pp. 8 and 11,
           https://ntrs.nasa.gov/citations/19930092205
    """
    flux = _nonnegative(coolant_mass_flux_kg_m2_s, "coolant_mass_flux_kg_m2_s")
    h_star = _positive(unblown_heat_transfer_coefficient_w_m2_k,
                       "unblown_heat_transfer_coefficient_w_m2_k")
    recovery = _positive(recovery_temperature_k, "recovery_temperature_k")
    coolant = _positive(coolant_temperature_k, "coolant_temperature_k")
    cp = _positive(coolant_specific_heat_j_kg_k, "coolant_specific_heat_j_kg_k")
    eta = _positive(wall_efficiency, "wall_efficiency")
    if eta > 1.0:
        raise DomainError("wall_efficiency must not exceed 1.", field="wall_efficiency")
    q_r = _nonnegative(radiative_heat_flux_w_m2, "radiative_heat_flux_w_m2")
    if not coolant < recovery:
        raise DomainError("coolant_temperature_k must be below recovery_temperature_k.",
                          field="coolant_temperature_k")
    phi = flux * cp / h_star
    theta = _blowing_ratio(phi)
    wall = (h_star * phi * eta * coolant + h_star * theta * recovery + q_r) / (
        h_star * phi * eta + h_star * theta)
    if not isfinite(wall):
        raise DomainError("Wall temperature is not finite.", code="non_finite_result")
    return _state(flux, wall, h_star, cp, recovery, coolant, eta, q_r)


# ---------------------------------------------------------------------------
# Closed-form temperature ratios
# ---------------------------------------------------------------------------


def eckert_livingood_temperature_ratio(coolant_flow_ratio: float, reynolds: float,
                                       prandtl: float, *,
                                       radiation_to_convection_ratio: float = 0.0,
                                       model: Literal["rannie_friedman", "film"] = "rannie_friedman",
                                       extrapolate: bool = False) -> CorrelationResult:
    """Turbulent flat-plate transpiration temperature ratio, NACA Report 1182.

    Eq. (31), which reduces to Eq. (30) without radiation:

        ``(T_w - T_a)/(T_g - T_a) = 1 / [1 + (rho_a v_a / rho_g V_g)
        (Re^0.2 Pr^(2/3) / 0.037) / (h_r/h_cv + r phi/(e^(r phi) - 1))]``

    with ``r phi = (2.11/0.037) Re^0.1 Pr^(2/3) rho_a v_a/(rho_g V_g)``
    (Rannie-Friedman, ``r = 2.11 Re^-0.1``). ``model="film"`` replaces ``r``
    by 1, the source's film-theory variant. ``T_a`` is the coolant supply
    temperature; the coolant leaves at the wall temperature; ``Re`` is the
    plate Reynolds number of the gas with the turbulent coefficient
    ``Nu = 0.037 Re^0.8 Pr^(1/3)`` (Eq. 8).

    Parameters
    ----------
    coolant_flow_ratio : float
        ``rho_a v_a / (rho_g V_g)``, nonnegative.
    reynolds : float
        Gas Reynolds number on plate length; the source computes 1e5 to 1e9.
    prandtl : float
        Gas Prandtl number (the source's curves use 0.7).
    radiation_to_convection_ratio : float, optional
        ``h_r / h_cv``.
    model : {"rannie_friedman", "film"}, optional
        Sublayer factor ``r``.
    extrapolate : bool, optional
        Warn instead of raising outside the computed coolant-flow range
        (0 to 0.012) and Reynolds range (1e5 to 1e9).

    Returns
    -------
    CorrelationResult
        Temperature-difference ratio ``(T_w - T_a)/(T_g - T_a)``.

    References
    ----------
    .. [1] NACA-REPORT-1182-ECKERT-LIVINGOOD-1954. Eckert & Livingood, NACA
           Report 1182, 1954, Eqs. (8), (24)-(27), (30), (31) and Fig. 7,
           pp. 8-11, https://ntrs.nasa.gov/citations/19930092205
    """
    ratio = _nonnegative(coolant_flow_ratio, "coolant_flow_ratio")
    re = _positive(reynolds, "reynolds")
    pr = _positive(prandtl, "prandtl")
    hr = _nonnegative(radiation_to_convection_ratio, "radiation_to_convection_ratio")
    if model not in ("rannie_friedman", "film"):
        raise DomainError("model must be 'rannie_friedman' or 'film'.", field="model")
    warnings: list[WarningMessage] = []
    for value, lo, hi, label in ((ratio, 0.0, 0.012, "coolant_flow_ratio"),
                                 (re, 1.0e5, 1.0e9, "reynolds")):
        if value < lo or value > hi:
            message = (f"{label} = {value:.6g} lies outside the range computed by "
                       f"Eckert & Livingood [{lo:g}, {hi:g}].")
            if not extrapolate:
                raise DomainError(message, code="outside_correlation_range", field=label)
            warnings.append(WarningMessage("correlation_extrapolated", message,
                                           "warning", label))
    r = 2.11 * re ** -0.1 if model == "rannie_friedman" else 1.0
    group = re ** 0.2 * pr ** (2.0 / 3.0) / 0.037
    r_phi = r * group * ratio
    value = 1.0 / (1.0 + ratio * group / (hr + _blowing_ratio(r_phi)))
    eq = "Eq. (31)" if hr > 0.0 else "Eq. (30)"
    return CorrelationResult(value, "temperature_difference_ratio",
                             f"Eckert-Livingood transpiration ({model})",
                             f"{_ECKERT}, {eq}", tuple(warnings))


def rannie_pipe_temperature_ratio(blowing_rate: float, reynolds: float,
                                  prandtl: float) -> CorrelationResult:
    """Rannie's porous-pipe relation as printed by Kelley & L'Ecuyer, Eq. (2.17).

    ``(T_inf - T_w)/(T_inf - T_0) = 1 - exp(-37 Pr F Re^0.1) /
    [1 + (1.18 Re^0.1 - 1)(1 - exp(-37 F Re^0.1))]``, from a laminar
    sublayer of ``y_a+ = 5.6`` under a turbulent core with
    ``c_f = 0.046 Re^-0.2``, Reynolds analogy and coolant properties in the
    sublayer. ``T_0`` is the coolant reservoir temperature and
    ``F = (rho v)_w/(rho u)_inf``.

    Parameters
    ----------
    blowing_rate : float
        ``F``, nonnegative.
    reynolds : float
        Pipe Reynolds number (Kelley's Appendix B example uses the plate
        ``Re_x``).
    prandtl : float
        Prandtl number.

    Returns
    -------
    CorrelationResult
        ``(T_inf - T_w)/(T_inf - T_0)``.

    Notes
    -----
    Rannie's original (JPL Progress Report 4-50, 1947) was not read; the
    equation is quoted from Kelley & L'Ecuyer. They report differences from
    Duwez-Wheeler data of about 10 % (mullite) and 30 % (copper) at
    ``F = 0.016``.

    References
    ----------
    .. [1] NASA-CR-76367-KELLEY-LECUYER-1966. Kelley & L'Ecuyer, NASA
           CR-76367, 1966, Eqs. (2.14)-(2.17), pp. 16-17, and Appendix B,
           pp. 147-148, https://ntrs.nasa.gov/citations/19660021566
    """
    f = _nonnegative(blowing_rate, "blowing_rate")
    re = _positive(reynolds, "reynolds")
    pr = _positive(prandtl, "prandtl")
    s = re ** 0.1
    value = 1.0 - exp(-37.0 * pr * f * s) / (1.0 + (1.18 * s - 1.0) * (-expm1(-37.0 * f * s)))
    return CorrelationResult(value, "temperature_difference_ratio", "Rannie porous pipe",
                             f"{_KELLEY}, Eq. (2.17)",
                             (WarningMessage("secondary_source",
                                             "Rannie's relation is quoted from Kelley & "
                                             "L'Ecuyer; the 1947 original was not read.",
                                             "info"),))


def rannie_pipe_blowing_rate(temperature_ratio: float, reynolds: float, prandtl: float) -> float:
    """Blowing rate ``F`` that gives a Rannie temperature ratio (inverse of Eq. 2.17).

    Parameters
    ----------
    temperature_ratio : float
        Target ``(T_inf - T_w)/(T_inf - T_0)`` in [0, 1).
    reynolds, prandtl : float
        As in :func:`rannie_pipe_temperature_ratio`.

    Returns
    -------
    float
        ``F``, found by bisection to 1e-14 relative (the ratio is monotone in
        ``F``).

    References
    ----------
    .. [1] NASA-CR-76367-KELLEY-LECUYER-1966. Kelley & L'Ecuyer, NASA
           CR-76367, 1966, Appendix B, pp. 147-148 (trial-and-error solution),
           https://ntrs.nasa.gov/citations/19660021566
    """
    target = real(temperature_ratio, "temperature_ratio")
    if not 0.0 <= target < 1.0:
        raise DomainError("temperature_ratio must lie in [0, 1).", field="temperature_ratio")
    lo, hi = 0.0, 1.0e-3
    while rannie_pipe_temperature_ratio(hi, reynolds, prandtl).value < target:
        hi *= 2.0
        if hi > 10.0:
            raise DomainError("No blowing rate reaches the temperature ratio.",
                              field="temperature_ratio")
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if rannie_pipe_temperature_ratio(mid, reynolds, prandtl).value < target:
            lo = mid
        else:
            hi = mid
        if hi - lo <= 1e-14 * hi:
            break
    return 0.5 * (lo + hi)


# ---------------------------------------------------------------------------
# Inside the porous wall
# ---------------------------------------------------------------------------


@dataclass(frozen=True, slots=True)
class PorousWallSolution:
    """Two-temperature solution through a porous wall (Pratt & Whitney III-24 to III-29).

    ``x = 0`` is the coolant (plenum) face, ``x = delta`` the heated face.

    Attributes
    ----------
    efficiency : float
        ``eta = (T_co - T_ci)/(T_w - T_ci)``, Eq. (III-29).
    coolant_exit_temperature_k : float
        ``T_co``, coolant temperature leaving the hot face.
    thickness_m : float
        ``delta``.
    a_per_m2, b_per_m : float
        ``a = h_i A_i / k_eff`` and ``b = h_i A_i / (m_c c_p,c)``.
    roots_per_m : tuple of float
        ``phi_1, phi_2 = -b/2 +/- sqrt(b^2 + 4a)/2``.
    coefficients_k : tuple of float
        ``C_1, C_2, C_3`` of Eqs. (III-26)-(III-28).

    References
    ----------
    .. [1] NASA-CR-107268-PWA-1969. Pratt & Whitney, NASA CR-107268, 1969,
           Eqs. (III-11)-(III-29), pp. III-13 to III-16,
           https://ntrs.nasa.gov/citations/19700004126
    """

    efficiency: float
    coolant_exit_temperature_k: float
    thickness_m: float
    a_per_m2: float
    b_per_m: float
    roots_per_m: tuple[float, float]
    coefficients_k: tuple[float, float, float]

    def solid_temperature_k(self, x_m: float) -> float:
        """Matrix temperature ``T_m(x)``, Eq. (III-24)."""
        x = real(x_m, "x_m")
        c1, c2, c3 = self.coefficients_k
        p1, p2 = self.roots_per_m
        return c1 + c2 * exp(p1 * x) + c3 * exp(p2 * x)

    def coolant_temperature_k(self, x_m: float) -> float:
        """Coolant temperature ``T_c(x)``, Eq. (III-25)."""
        x = real(x_m, "x_m")
        c1, c2, c3 = self.coefficients_k
        p1, p2 = self.roots_per_m
        k = self.b_per_m / self.a_per_m2
        return c1 + c2 * k * p1 * exp(p1 * x) + c3 * k * p2 * exp(p2 * x)


def porous_wall_solution(*, thickness_m: float, volumetric_heat_transfer_w_m3_k: float,
                         effective_conductivity_w_m_k: float, coolant_mass_flux_kg_m2_s: float,
                         coolant_specific_heat_j_kg_k: float, wall_temperature_k: float,
                         coolant_inlet_temperature_k: float) -> PorousWallSolution:
    """Solid and coolant temperatures through a porous wall without local equilibrium.

    Solves Eqs. (III-11) and (III-14),
    ``T_m'' - a (T_m - T_c) = 0`` and ``T_c' - b (T_m - T_c) = 0``, with
    ``T_m(delta) = T_w``, ``T_m'(0) = 0`` (no conduction into the plenum) and
    ``T_c(0) = T_ci``.

    Parameters
    ----------
    thickness_m : float
        Wall thickness ``delta``.
    volumetric_heat_transfer_w_m3_k : float
        ``h_i A_i``: internal coefficient times pore surface per unit volume.
    effective_conductivity_w_m_k : float
        Effective conductivity of the porous matrix, ``k_eff``.
    coolant_mass_flux_kg_m2_s : float
        Coolant mass flux through the wall, ``m_c``.
    coolant_specific_heat_j_kg_k : float
        ``c_p,c``.
    wall_temperature_k : float
        Hot-face temperature ``T_w``.
    coolant_inlet_temperature_k : float
        Plenum coolant temperature ``T_ci``.

    Returns
    -------
    PorousWallSolution

    References
    ----------
    .. [1] NASA-CR-107268-PWA-1969. Pratt & Whitney, NASA CR-107268, 1969,
           Eqs. (III-11)-(III-29), pp. III-13 to III-16,
           https://ntrs.nasa.gov/citations/19700004126
    """
    delta = _positive(thickness_m, "thickness_m")
    hia = _positive(volumetric_heat_transfer_w_m3_k, "volumetric_heat_transfer_w_m3_k")
    k_eff = _positive(effective_conductivity_w_m_k, "effective_conductivity_w_m_k")
    m_c = _positive(coolant_mass_flux_kg_m2_s, "coolant_mass_flux_kg_m2_s")
    cp = _positive(coolant_specific_heat_j_kg_k, "coolant_specific_heat_j_kg_k")
    t_w = _positive(wall_temperature_k, "wall_temperature_k")
    t_ci = _positive(coolant_inlet_temperature_k, "coolant_inlet_temperature_k")
    if t_w <= t_ci:
        raise DomainError("wall_temperature_k must exceed coolant_inlet_temperature_k.",
                          field="wall_temperature_k")
    a = hia / k_eff
    b = hia / (m_c * cp)
    beta = 0.5 * sqrt(b * b + 4.0 * a)
    p1, p2 = -0.5 * b + beta, -0.5 * b - beta
    if p1 * delta > 700.0 or -p2 * delta > 700.0:
        raise DomainError("Porous wall is too thick for the closed form (overflow); "
                          "the coolant leaves at the wall temperature.",
                          code="non_finite_result", field="thickness_m")
    dt = t_w - t_ci
    c2 = dt / (exp(p1 * delta) - (p1 / p2) * exp(p2 * delta))
    c3 = dt / (exp(p2 * delta) - (p2 / p1) * exp(p1 * delta))
    eta = (exp(beta * delta) - exp(-beta * delta)) / (
        0.5 * (2.0 * beta / b + 1.0) * exp(beta * delta)
        + 0.5 * (2.0 * beta / b - 1.0) * exp(-beta * delta))
    return PorousWallSolution(eta, t_ci + eta * dt, delta, a, b, (p1, p2), (t_ci, c2, c3))
