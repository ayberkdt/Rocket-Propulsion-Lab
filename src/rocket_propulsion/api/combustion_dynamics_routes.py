"""HTTP routes for transpiration cooling, combustion stability and propellant vaporization.

Handlers translate JSON to library calls; the physics and its sources are in
:mod:`rocket_propulsion.heat_transfer.transpiration`,
:mod:`rocket_propulsion.thermochemistry.combustion_instability`,
:mod:`rocket_propulsion.thermochemistry.feed_coupled_instability` and
:mod:`rocket_propulsion.thermochemistry.vaporization`. Responses are strict
JSON.

References
----------
.. [1] NACA-TN-3208-MICKLEY-1954. Mickley et al., NACA TN 3208, 1954,
       https://ntrs.nasa.gov/citations/19930084056
.. [2] NASA-SP-194-1972. Harrje and Reardon (eds.), NASA SP-194, 1972,
       https://ntrs.nasa.gov/citations/19720026079
.. [3] NASA-CR-168066-AEROJET-1983. J. Fang, NASA CR-168066, 1983,
       https://ntrs.nasa.gov/citations/19830023615
"""

from __future__ import annotations

import cmath
from math import pi, sqrt
from typing import Any

from rocket_propulsion.api.aero_routes import _clean, _linspace, _number, _positive, _text
from rocket_propulsion.compressible.isentropic import mach_from_isentropic
from rocket_propulsion.core.constants import USSA1976_SEA_LEVEL_PRESSURE
from rocket_propulsion.core.errors import InputError
from rocket_propulsion.core.metadata import ModelMetadata
from rocket_propulsion.core.validation import numeric_range_guard
from rocket_propulsion.heat_transfer import transpiration
from rocket_propulsion.thermochemistry import (
    EquilibriumProblem,
    LocalEquilibriumProvider,
    get_propellant_pair,
    reactants_for_mass_ratio,
    vaporization,
)
from rocket_propulsion.thermochemistry import combustion_instability as instability
from rocket_propulsion.thermochemistry import feed_coupled_instability as feed
from rocket_propulsion.thermochemistry.combustion_efficiency import (
    cea_characteristic_velocity_function,
)
from rocket_propulsion.thermodynamics.real_fluids import reference_fluid, reference_fluid_model


def transpiration_route(payload: dict[str, Any]) -> dict[str, Any]:
    """Coolant flux for a target porous-wall temperature and the flux sweep [1]_.

    Parameters
    ----------
    payload : dict
        unblown_heat_transfer_coefficient_w_m2_k : float, default 2.0e4
        recovery_temperature_k : float, default 3300
        target_wall_temperature_k : float, default 900
        coolant_temperature_k : float, default 100
        coolant_specific_heat_j_kg_k : float, default 14300 (hydrogen)
        wall_efficiency : float, default 1
        radiative_heat_flux_w_m2 : float, default 0

    Returns
    -------
    dict
        The required state and a sweep of wall temperature and blown
        coefficient over 0 to 2.5 times the required coolant flux.

    References
    ----------
    .. [1] NACA-TN-3208-MICKLEY-1954. Mickley et al., NACA TN 3208, 1954,
           Eqs. (27)-(30), https://ntrs.nasa.gov/citations/19930084056
    """
    common = {
        "unblown_heat_transfer_coefficient_w_m2_k": _number(payload, "unblown_heat_transfer_coefficient_w_m2_k", 2.0e4),
        "recovery_temperature_k": _number(payload, "recovery_temperature_k", 3300.0),
        "coolant_temperature_k": _number(payload, "coolant_temperature_k", 100.0),
        "coolant_specific_heat_j_kg_k": _number(payload, "coolant_specific_heat_j_kg_k", 14300.0),
        "wall_efficiency": _number(payload, "wall_efficiency", 1.0),
        "radiative_heat_flux_w_m2": _number(payload, "radiative_heat_flux_w_m2", 0.0),
    }
    target = transpiration.transpiration_coolant_flux(
        wall_temperature_k=_number(payload, "target_wall_temperature_k", 900.0), **common)
    sweep = []
    for flux in _linspace(0.0, 2.5 * target.coolant_mass_flux_kg_m2_s, 26):
        state = transpiration.transpiration_wall_temperature(coolant_mass_flux_kg_m2_s=flux, **common)
        sweep.append({"coolant_mass_flux_kg_m2_s": flux, "wall_temperature_k": state.wall_temperature_k,
                      "blowing_correction": state.blowing_correction,
                      "convective_heat_flux_w_m2": state.convective_heat_flux_w_m2})
    return _clean({"required": target, "sweep": sweep, "warnings": list(target.warnings)})


def _chamber_state(pair_key: str, pressure: float, ratio: float) -> dict[str, float]:
    pair = get_propellant_pair(pair_key)
    chamber = LocalEquilibriumProvider().solve(
        EquilibriumProblem.hp(reactants_for_mass_ratio(pair, ratio), pressure_pa=pressure))
    gas_constant = chamber.gas_constant_j_kg_k
    return {"temperature_k": chamber.temperature_k, "gamma": chamber.gamma,
            "sound_speed_m_s": sqrt(chamber.gamma * gas_constant * chamber.temperature_k),
            "density_kg_m3": pressure / (gas_constant * chamber.temperature_k)}


@numeric_range_guard("chamber_radius_m")
def combustion_stability_route(payload: dict[str, Any]) -> dict[str, Any]:
    """Chamber acoustic modes and chugging stability for a propellant pair [1]_.

    The chamber state (temperature, frozen gamma, density, sound speed
    ``sqrt(gamma R T)``) and ``c*(O/F)`` come from the CEA-gated equilibrium
    solver; ``dc*/dr`` is a central difference (+/- 1 % in O/F). The gas
    residence time is ``theta_c = rho V_c / m_dot`` with ``V_c = pi R^2 L`` and
    ``m_dot = p_c A* / c*``, ``A* = pi R^2 / CR``. Feeds are rigid (purely
    resistive injectors), as in SP-194 Example 5.5.1.1.

    Parameters
    ----------
    payload : dict
        propellant_pair : str, default "LOX_CH4"
        oxidizer_fuel_ratio : float, default the pair's reference
        chamber_pressure_pa : float, default 7e6
        chamber_radius_m, chamber_length_m : float, default 0.10, 0.30
        contraction_ratio : float, default 3
        oxidizer_drop_fraction, fuel_drop_fraction : float, default 0.2, 0.2
            Injector drop over chamber pressure.
        combustion_time_lag_s : float, default 0.8e-3

    Returns
    -------
    dict
        Chamber state, mode table up to 3 x the 1T frequency, chugging gain
        and limit, the Nyquist verdict and the ``f(i omega)`` curve.

    References
    ----------
    .. [1] NASA-SP-194-1972. NASA SP-194, 1972, Table 3.5.1, Eqs.
           (3.5.1-19), (5.2.2-6), (5.5.1-6), Sect. 5.5.1,
           https://ntrs.nasa.gov/citations/19720026079
    """
    pair_key = _text(payload, "propellant_pair", "LOX_CH4")
    pair = get_propellant_pair(pair_key)
    ratio = _number(payload, "oxidizer_fuel_ratio", pair.reference_oxidizer_fuel_ratio)
    pressure = _number(payload, "chamber_pressure_pa", 7.0e6)
    radius = _positive(payload, "chamber_radius_m", 0.10)
    length = _positive(payload, "chamber_length_m", 0.30)
    contraction = _number(payload, "contraction_ratio", 3.0)
    if contraction <= 1.0:
        raise InputError("contraction_ratio must exceed 1.", field="contraction_ratio")
    ox_drop = _number(payload, "oxidizer_drop_fraction", 0.2)
    fuel_drop = _number(payload, "fuel_drop_fraction", 0.2)
    lag = _number(payload, "combustion_time_lag_s", 0.8e-3)
    state = _chamber_state(pair_key, pressure, ratio)
    cstar, _ = cea_characteristic_velocity_function(propellant_pair=pair_key, chamber_pressure_pa=pressure)
    c0 = cstar(ratio)
    slope = (cstar(1.01 * ratio) - cstar(0.99 * ratio)) / (0.02 * ratio)
    throat_area = pi * radius * radius / contraction
    mass_flow = pressure * throat_area / c0
    theta = state["density_kg_m3"] * pi * radius * radius * length / mass_flow
    mach = mach_from_isentropic("area_ratio_subsonic", contraction, state["gamma"])
    one_t = instability.acoustic_mode_frequency(sound_speed_m_s=state["sound_speed_m_s"],
                                                chamber_radius_m=radius, chamber_length_m=length,
                                                tangential_order=1)
    modes = instability.chamber_acoustic_modes(
        sound_speed_m_s=state["sound_speed_m_s"], chamber_radius_m=radius, chamber_length_m=length,
        maximum_frequency_hz=3.0 * one_t.frequency_hz, mean_mach_number=mach, max_order=4)
    chug = instability.chugging_stability(
        mixture_ratio=ratio, characteristic_velocity_m_s=c0, cstar_mixture_ratio_derivative_m_s=slope,
        chamber_pressure_pa=pressure, oxidizer_injector_drop_pa=ox_drop * pressure,
        fuel_injector_drop_pa=fuel_drop * pressure, gas_residence_time_s=theta, time_lag_s=lag)
    g_ox = (ratio / (1.0 + ratio)) / (2.0 * ox_drop)
    g_f = (1.0 / (1.0 + ratio)) / (2.0 * fuel_drop)
    nyquist = feed.feed_coupled_stability(
        mixture_ratio=ratio, characteristic_velocity_m_s=c0, cstar_mixture_ratio_derivative_m_s=slope,
        chamber_pressure_pa=pressure, gas_residence_time_s=theta, combustion_time_lag_s=lag,
        oxidizer_admittance=lambda s: g_ox, fuel_admittance=lambda s: g_f,
        oxidizer_admittance_bounds=feed.AdmittanceBounds(g_ox, 0.0),
        fuel_admittance_bounds=feed.AdmittanceBounds(g_f, 0.0), points=4000)
    curve = []
    k_ox = 1.0 + (1.0 + ratio) / c0 * slope
    k_f = 1.0 - ratio * (1.0 + ratio) / c0 * slope
    for f in _linspace(0.0, 4.0 * max(chug.neutral_frequency_hz or 0.0, 1.0 / (2.0 * pi * theta)), 241):
        s = 2j * pi * f
        value = cmath.exp(-s * lag) / (1.0 + s * theta) * (k_ox * g_ox + k_f * g_f)
        curve.append({"frequency_hz": f, "real": value.real, "imaginary": value.imag})
    return _clean({
        "chamber": {**state, "characteristic_velocity_m_s": c0, "cstar_mixture_ratio_derivative_m_s": slope,
                    "mass_flow_kg_s": mass_flow, "residence_time_s": theta, "mean_mach_number": mach},
        "modes": [{"name": m.name, "frequency_hz": m.frequency_hz, "eigenvalue": m.eigenvalue} for m in modes],
        "chugging": chug,
        "nyquist": {"stable": nyquist.stable, "right_half_plane_zeros": nyquist.right_half_plane_zeros,
                    "gain_crossover_hz": nyquist.gain_crossover_hz,
                    "phase_margins_deg": nyquist.phase_margins_deg, "curve": curve},
    })


_PROPELLANT_FLUIDS = {"LOX_CH4": ("oxygen", "methane"), "LOX_LH2": ("oxygen", None)}


def _liquid_properties(fluid: str) -> dict[str, float]:
    model = reference_fluid_model(fluid)
    saturation = reference_fluid(fluid).saturation_p(USSA1976_SEA_LEVEL_PRESSURE.value)
    return {"injection_temperature_k": saturation.temperature_k, "critical_temperature_k": model.t_crit_k,
            "heat_of_vaporization_j_kg": saturation.enthalpy_of_vaporization_j_kg,
            "molecular_weight": model.molar_mass_kg_mol * 1000.0}


def propellant_vaporization_route(payload: dict[str, Any]) -> dict[str, Any]:
    """Priem-Heidmann vaporization and the vaporization-limited c* efficiency [1]_.

    Propellants are injected as saturated liquid at 1 atm; ``T_c``, ``M`` and
    ``H_v`` at that temperature come from the package's reference equations
    of state. Hydrogen in LOX/LH2 is treated as injected gas (fully
    vaporized).

    Parameters
    ----------
    payload : dict
        propellant_pair : {"LOX_CH4", "LOX_LH2"}, default "LOX_CH4"
        chamber_pressure_pa : float, default 7e6
        chamber_length_m : float, default 0.25 (cylindrical length)
        contraction_ratio : float, default 3
        drop_radius_m : float, default 7.62e-5 (0.003 in.)
        injection_velocity_m_s : float, default 30.48 (100 ft/s)

    Returns
    -------
    dict
        Per-propellant properties, generalized length and percent vaporized,
        the c* efficiency, and a sweep over chamber length.

    References
    ----------
    .. [1] NASA-CR-168066-AEROJET-1983. J. Fang, NASA CR-168066, 1983,
           https://ntrs.nasa.gov/citations/19830023615
    """
    pair_key = _text(payload, "propellant_pair", "LOX_CH4")
    if pair_key not in _PROPELLANT_FLUIDS:
        raise InputError("propellant_pair must be LOX_CH4 or LOX_LH2 (reference equations of state "
                         "exist for oxygen and methane).", field="propellant_pair")
    pressure = _number(payload, "chamber_pressure_pa", 7.0e6)
    length = _number(payload, "chamber_length_m", 0.25)
    contraction = _number(payload, "contraction_ratio", 3.0)
    radius = _number(payload, "drop_radius_m", 7.62e-5)
    velocity = _number(payload, "injection_velocity_m_s", 30.48)
    pair = get_propellant_pair(pair_key)
    ratio = pair.reference_oxidizer_fuel_ratio
    ox_fluid, fuel_fluid = _PROPELLANT_FLUIDS[pair_key]
    props = {"oxidizer": _liquid_properties(ox_fluid),
             "fuel": None if fuel_fluid is None else _liquid_properties(fuel_fluid)}
    gas_length = vaporization.PRIEM_HEIDMANN_CURVE[-1][0] * 0.0254

    def lengths(chamber_length: float) -> tuple[float, float]:
        out = []
        for key in ("oxidizer", "fuel"):
            p = props[key]
            if p is None:
                out.append(gas_length)
                continue
            value = vaporization.generalized_vaporization_length(
                chamber_length_m=chamber_length, chamber_pressure_pa=pressure, contraction_ratio=contraction,
                liquid_temperature_k=p["injection_temperature_k"], critical_temperature_k=p["critical_temperature_k"],
                heat_of_vaporization_j_kg=p["heat_of_vaporization_j_kg"], molecular_weight=p["molecular_weight"],
                drop_radius_m=radius, injection_velocity_m_s=velocity)
            out.append(min(max(value, 0.01 * 0.0254), gas_length))
        return out[0], out[1]

    cstar, _ = cea_characteristic_velocity_function(propellant_pair=pair_key, chamber_pressure_pa=pressure)
    ox_len, fuel_len = lengths(length)
    result = vaporization.vaporization_cstar_efficiency(
        characteristic_velocity=cstar, oxidizer_fuel_ratio=ratio,
        oxidizer_generalized_length_m=ox_len, fuel_generalized_length_m=fuel_len)
    sweep = []
    for chamber_length in _linspace(0.05, 2.0 * length, 16):
        a, b = lengths(chamber_length)
        sweep.append({"chamber_length_m": chamber_length,
                      "oxidizer_percent_vaporized": vaporization.percent_vaporized(a),
                      "fuel_percent_vaporized": vaporization.percent_vaporized(b)})
    return _clean({"properties": props, "oxidizer_fuel_ratio": ratio,
                   "oxidizer_generalized_length_m": ox_len, "fuel_generalized_length_m": fuel_len,
                   "result": result, "sweep": sweep})


COMBUSTION_DYNAMICS_ROUTES = {
    "/api/v1/transpiration-cooling": transpiration_route,
    "/api/v1/combustion-stability": combustion_stability_route,
    "/api/v1/propellant-vaporization": propellant_vaporization_route,
}

COMBUSTION_DYNAMICS_ROUTE_METADATA = {
    "/api/v1/transpiration-cooling": ModelMetadata(
        model="film-theory transpiration cooling with wall energy balance", version="1.0",
        assumptions=("Film theory h/h* = phi/(e^phi - 1)", "Coolant leaves at T_c + eta (T_w - T_c)",
                     "Radiation not blocked by blowing"),
        units={"temperature": "K", "mass_flux": "kg/(m^2 s)"},
        references=("NACA TN 3208", "NASA CR-107268", "NACA Report 1182")),
    "/api/v1/combustion-stability": ModelMetadata(
        model="chamber acoustics and single-time-lag chugging (SP-194)", version="1.0",
        assumptions=("Rigid closed cylinder, frozen sound speed", "Rigid resistive injectors",
                     "Constant combustion time lag"),
        units={"frequency": "Hz", "time": "s"}, references=("NASA SP-194",)),
    "/api/v1/propellant-vaporization": ModelMetadata(
        model="Priem-Heidmann generalized length and TN 4219 c* efficiency", version="1.0",
        assumptions=("TR R-67 curve as tabulated in NASA CR-168066", "Saturated liquid at 1 atm",
                     "Vaporized propellant burns at once"),
        units={"length": "m", "percent": "%"}, references=("NASA CR-168066", "NACA TN 3883", "NACA TN 4219")),
}
