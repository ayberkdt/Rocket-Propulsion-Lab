"""Coupled gas-wall-coolant thermal analysis of a regeneratively cooled thrust chamber.

The solver follows the station-marching scheme of NASA's Rocket Thermal
Evaluation code as documented by Wadel [1]_: the coolant is marched along its
flow direction and, at every wall station, the gas-side film coefficient, the
conduction through the hot wall and the coolant-side film coefficient (with
the channel lands treated as fins) are solved simultaneously for the hot-wall
temperature. Each piece reuses a verified module of this package:

* gas side: Bartz's correlation with the local wall temperature [2]_ — the same
  expression as :func:`rocket_propulsion.thermodynamics.calculate_bartz_heat_transfer`,
  evaluated per station — and an adiabatic-wall (recovery) temperature with
  recovery factor ``Pr^(1/3)`` [3]_;
* optional gas radiation from a weighted-sum-of-gray-gases model, only where
  that model's fit range holds (:mod:`.radiation`);
* hot-wall conduction ``q = k (T_wg - T_wc) / t`` and straight-fin lands
  (:mod:`.conduction`) [4]_;
* coolant side: Gnielinski with the Chen friction factor, or Taylor's hydrogen
  correlation (:mod:`.convection`), with density, heat capacity, viscosity and
  conductivity from the NIST-verified real-fluid equations of state and
  transport correlations (:mod:`rocket_propulsion.thermodynamics.real_fluids`);
* coolant energy ``dh = q dA / mdot`` and Darcy pressure loss along the channel.

The march is explicit in the coolant direction (first order in station
spacing); refinement studies are in ``tests/test_regen_chamber.py``.

References
----------
.. [1] NASA-TM-1998-206313-WADEL. M. F. Wadel, "Comparison of High Aspect Ratio Cooling
       Channel Designs for a Rocket Combustion Chamber With Development of an
       Optimized Design", NASA/TM-1998-206313.
.. [2] NASA-SP125-BARTZ-PAGES. NASA SP-125, printed pp. 100-101, equations
       4-10, 4-13 and 4-14,
       https://ntrs.nasa.gov/api/citations/19710019929/downloads/19710019929.pdf
.. [3] NASA-HOPPER-2012-THERMAL. Linne et al., AIAA 2012-0806, printed p. 8,
       equations 1, 2 and 4,
       https://ntrs.nasa.gov/api/citations/20120004103/downloads/20120004103.pdf
.. [4] LIENHARD-AHTT-6.00. J. H. Lienhard IV and J. H. Lienhard V, "A Heat
       Transfer Textbook", v6.00, 2024, Ch. 2 (resistances) and Ch. 4 (fins).
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from math import isfinite, pi, sqrt
from typing import Literal

from rocket_propulsion.core.errors import ConvergenceError, DomainError, InputError
from rocket_propulsion.core.metadata import WarningMessage
from rocket_propulsion.core.validation import numeric_range_guard
from rocket_propulsion.geometry import NozzleContour
from rocket_propulsion.heat_transfer import conduction, convection, radiation
from rocket_propulsion.thermodynamics.materials import get_wall_material
from rocket_propulsion.thermodynamics.real_fluids import reference_fluid
from rocket_propulsion.thermodynamics.real_fluids.transport import transport_tp

CoolantCorrelation = Literal["auto", "gnielinski", "taylor"]
GasRadiation = Literal["none", "wsgg"]
_HYDROGEN = {"parahydrogen", "normal_hydrogen"}


def _positive(value: float, label: str) -> float:
    number = float(value)
    if not isfinite(number) or number <= 0.0:
        raise DomainError(f"{label} must be finite and positive.", field=label)
    return number


@dataclass(frozen=True, slots=True)
class HotGas:
    """Chamber stagnation state and frozen gas properties for the gas side [1]_.

    Attributes
    ----------
    chamber_pressure_pa, chamber_temperature_k : float
        Stagnation pressure and temperature.
    characteristic_velocity_m_s : float
        ``c*`` of the propellant.
    specific_heat_j_kg_k, viscosity_pa_s, prandtl_number, gamma : float
        Frozen chamber properties used by Bartz's correlation.
    throat_radius_of_curvature_m : float
        Throat wall radius of curvature.
    h2o_mole_fraction, co2_mole_fraction : float
        Frozen chamber mole fractions used by the optional gas radiation.

    References
    ----------
    .. [1] NASA-SP125-BARTZ-PAGES. NASA SP-125, printed pp. 100-101,
           https://ntrs.nasa.gov/api/citations/19710019929/downloads/19710019929.pdf
    """

    chamber_pressure_pa: float
    chamber_temperature_k: float
    characteristic_velocity_m_s: float
    specific_heat_j_kg_k: float
    viscosity_pa_s: float
    prandtl_number: float
    gamma: float
    throat_radius_of_curvature_m: float
    h2o_mole_fraction: float = 0.0
    co2_mole_fraction: float = 0.0


@dataclass(frozen=True, slots=True)
class CoolingChannels:
    """Rectangular milled channels around the circumference [1]_.

    Attributes
    ----------
    count : int
        Number of channels.
    width_m, height_m : float
        Channel width (circumferential) and depth (radial).
    hot_wall_thickness_m : float
        Wall between the hot gas and the channel bottom.
    roughness_m : float
        Absolute surface roughness of the channel.

    References
    ----------
    .. [1] NASA-TM-1998-206313-WADEL. M. F. Wadel, NASA/TM-1998-206313 (channel geometry).
    """

    count: int
    width_m: float
    height_m: float
    hot_wall_thickness_m: float
    roughness_m: float = 0.0


@dataclass(frozen=True, slots=True)
class Coolant:
    """Coolant inlet state and flow arrangement.

    Attributes
    ----------
    fluid : str
        Real-fluid key (``"parahydrogen"``, ``"methane"``, ``"oxygen"``, ...).
    inlet_temperature_k, inlet_pressure_pa : float
        Channel inlet state.
    mass_flow_kg_s : float
        Total coolant mass flow over all channels.
    counterflow : bool, default True
        Enter at the nozzle exit and flow toward the injector.
    """

    fluid: str
    inlet_temperature_k: float
    inlet_pressure_pa: float
    mass_flow_kg_s: float
    counterflow: bool = True


@dataclass(frozen=True, slots=True)
class ThermalStation:
    """Solved temperatures, fluxes and coolant state at one wall station."""

    x_m: float
    radius_m: float
    area_ratio: float
    mach: float
    adiabatic_wall_temperature_k: float
    gas_side_coefficient_w_m2_k: float
    radiative_heat_flux_w_m2: float
    heat_flux_w_m2: float
    hot_wall_temperature_k: float
    cold_wall_temperature_k: float
    coolant_temperature_k: float
    coolant_pressure_pa: float
    coolant_side_coefficient_w_m2_k: float
    coolant_reynolds: float
    land_fin_efficiency: float
    wall_limit_margin_k: float | None


@dataclass(frozen=True, slots=True)
class RegenChamberResult:
    """Station results and totals of a regenerative cooling solution.

    Attributes
    ----------
    stations : tuple of ThermalStation
        In contour order (injector to exit).
    total_heat_load_w : float
        Heat absorbed by the coolant.
    coolant_outlet_temperature_k, coolant_outlet_pressure_pa : float
        Channel outlet state.
    coolant_pressure_drop_pa : float
        Inlet minus outlet pressure.
    maximum_hot_wall_temperature_k : float
        Largest gas-side wall temperature.
    wall_limit_temperature_k : float or None
        Material limit when a catalog material is used.
    energy_balance_residual_w : float
        ``mdot (h_out - h_in) - Q`` with ``h`` re-evaluated from the EOS at the
        outlet state; measures the enthalpy inversion, not the march.
    warnings : tuple of WarningMessage
        Model statements and range notes.
    """

    stations: tuple[ThermalStation, ...]
    total_heat_load_w: float
    coolant_outlet_temperature_k: float
    coolant_outlet_pressure_pa: float
    coolant_pressure_drop_pa: float
    maximum_hot_wall_temperature_k: float
    wall_limit_temperature_k: float | None
    energy_balance_residual_w: float
    warnings: tuple[WarningMessage, ...]


def bartz_coefficient(gas: HotGas, *, throat_radius_m: float, area_ratio: float, mach: float,
                      wall_temperature_k: float) -> tuple[float, float]:
    """Return ``(h_g, T_aw)`` at one station with its own wall temperature [1]_, [2]_.

    With ``C = 1 + (gamma - 1) M^2 / 2`` and ``B = (T_w/T_c C + 1) / 2``,
    ``h = 0.026 cp mu^0.2 Pr^-0.6 D^-0.2 (p_c/c*)^0.8 (D/R_c)^0.1 (A/A*)^-0.9 B^-0.68 C^-0.12``
    and ``T_aw = T_c (1 + Pr^(1/3) (C - 1)) / C`` — the expression used by
    :func:`rocket_propulsion.thermodynamics.calculate_bartz_heat_transfer`.

    References
    ----------
    .. [1] NASA-SP125-BARTZ-PAGES. NASA SP-125, printed pp. 100-101, equations
           4-10, 4-13 and 4-14,
           https://ntrs.nasa.gov/api/citations/19710019929/downloads/19710019929.pdf
    .. [2] NASA-HOPPER-2012-THERMAL. Linne et al., AIAA 2012-0806, printed p. 8,
           equations 1, 2 and 4,
           https://ntrs.nasa.gov/api/citations/20120004103/downloads/20120004103.pdf
    """

    diameter = 2.0 * throat_radius_m
    compression = 1.0 + 0.5 * (gas.gamma - 1.0) * mach * mach
    sigma_base = 0.5 * (wall_temperature_k / gas.chamber_temperature_k * compression + 1.0)
    coefficient = (
        0.026 * gas.specific_heat_j_kg_k * gas.viscosity_pa_s**0.2 * gas.prandtl_number**-0.6
        * diameter**-0.2 * (gas.chamber_pressure_pa / gas.characteristic_velocity_m_s) ** 0.8
        * (diameter / gas.throat_radius_of_curvature_m) ** 0.1 * area_ratio**-0.9
        * sigma_base**-0.68 * compression**-0.12
    )
    recovery = gas.chamber_temperature_k * (1.0 + gas.prandtl_number ** (1.0 / 3.0) * (compression - 1.0)) / compression
    return coefficient, recovery


def _temperature_from_enthalpy(fluid: str, enthalpy: float, pressure: float, guess: float) -> float:
    """Invert ``h(T, p) = enthalpy``: Newton on ``dh/dT = cp``, bisection as fallback."""
    model = reference_fluid(fluid)
    temperature = guess
    for _ in range(12):
        try:
            state = model.state_tp(temperature, pressure)
        except DomainError:
            break
        step = (state.enthalpy_j_kg - enthalpy) / state.cp_j_kg_k
        temperature -= step
        if abs(step) < 1.0e-9 * temperature:
            return temperature

    def residual(temperature: float) -> float:
        return model.state_tp(temperature, pressure).enthalpy_j_kg - enthalpy

    low, high = max(model.model.t_min_k * 1.0001, guess * 0.8), min(model.model.t_max_k, guess * 1.5 + 5.0)
    f_low, f_high = residual(low), residual(high)
    for _ in range(40):
        if f_low <= 0.0 <= f_high:
            break
        if f_low > 0.0:
            low = max(model.model.t_min_k * 1.0001, low * 0.8)
            f_low = residual(low)
        if f_high < 0.0:
            high = min(model.model.t_max_k, high * 1.4)
            f_high = residual(high)
    else:
        raise ConvergenceError("Coolant temperature could not be bracketed from enthalpy.",
                               code="coolant_enthalpy_inversion")
    for _ in range(200):
        middle = 0.5 * (low + high)
        value = residual(middle)
        if value > 0.0:
            high = middle
        else:
            low = middle
        if high - low < 1.0e-9 * middle:
            break
    return 0.5 * (low + high)


@numeric_range_guard("contour")
def solve_regen_chamber(
    *,
    contour: NozzleContour,
    gas: HotGas,
    channels: CoolingChannels,
    coolant: Coolant,
    wall_material: str | None = "grcop-84",
    wall_thermal_conductivity_w_m_k: float | None = None,
    coolant_correlation: CoolantCorrelation = "auto",
    gas_radiation: GasRadiation = "none",
) -> RegenChamberResult:
    """March the coupled gas-wall-coolant problem along the cooling channels [1]_, [2]_.

    Parameters
    ----------
    contour : NozzleContour
        Wall stations with radius, area ratio and Mach number.
    gas : HotGas
        Chamber state and frozen gas properties.
    channels : CoolingChannels
        Channel geometry; the land width follows from the local circumference.
    coolant : Coolant
        Coolant fluid, inlet state, flow and direction.
    wall_material : str or None, default "grcop-84"
        Catalog material for conductivity and the temperature limit.
    wall_thermal_conductivity_w_m_k : float, optional
        Overrides the catalog conductivity.
    coolant_correlation : {"auto", "gnielinski", "taylor"}, default "auto"
        Coolant-side correlation. ``"auto"`` uses Taylor's correlation for
        hydrogen, whose wall-to-bulk temperature-ratio term matters for
        supercritical hydrogen (NASA TN D-4332), and Gnielinski otherwise;
        Taylor applies to hydrogen only.
    gas_radiation : {"none", "wsgg"}, default "none"
        Add H2O/CO2 gas radiation where the WSGG fit range holds.

    Returns
    -------
    RegenChamberResult
        Station temperatures and fluxes, coolant state and totals.

    Raises
    ------
    InputError
        For unknown options or a hydrogen-only correlation on another fluid.
    DomainError
        For non-physical geometry (no land between channels), a coolant
        state outside its equation range, or a two-phase coolant.
    ConvergenceError
        If a station's wall temperature cannot be bracketed.

    Notes
    -----
    At each station the hot-wall temperature solves
    ``q_g(T_wg) = (T_wg - T_b) / (t/k + W / (h_c (w + 2 eta h)))``, where
    ``q_g = h_g(T_wg) (T_aw - T_wg) + q_rad``, ``W`` is the circumferential
    width per channel at the hot wall and ``eta`` the efficiency of a land fin
    of thickness ``W - w`` and height ``h`` [3]_. The coolant then absorbs
    ``q_g 2 pi r ds`` over the next segment and loses ``f (ds/D_h) G^2/(2 rho)``
    in pressure (Chen friction factor, as in RTE [1]_).

    References
    ----------
    .. [1] NASA-TM-1998-206313-WADEL. M. F. Wadel, NASA/TM-1998-206313.
    .. [2] NASA-SP125-BARTZ-PAGES. NASA SP-125, printed pp. 100-101,
           https://ntrs.nasa.gov/api/citations/19710019929/downloads/19710019929.pdf
    .. [3] LIENHARD-AHTT-6.00. J. H. Lienhard IV and J. H. Lienhard V, "A Heat
           Transfer Textbook", v6.00, 2024, Ch. 4.
    """

    if coolant_correlation not in ("auto", "gnielinski", "taylor"):
        raise InputError("Coolant correlation must be auto, gnielinski or taylor.", field="coolant_correlation")
    if coolant_correlation == "auto":
        coolant_correlation = "taylor" if coolant.fluid in _HYDROGEN else "gnielinski"
    if coolant_correlation == "taylor" and coolant.fluid not in _HYDROGEN:
        raise InputError("Taylor's correlation applies to hydrogen coolant only.", field="coolant_correlation")
    if gas_radiation not in ("none", "wsgg"):
        raise InputError("Gas radiation must be none or wsgg.", field="gas_radiation")
    if channels.count < 1:
        raise DomainError("At least one channel is required.", field="count")
    width = _positive(channels.width_m, "width_m")
    height = _positive(channels.height_m, "height_m")
    thickness = _positive(channels.hot_wall_thickness_m, "hot_wall_thickness_m")
    mass_flow = _positive(coolant.mass_flow_kg_s, "mass_flow_kg_s")
    limit = None
    if wall_thermal_conductivity_w_m_k is not None:
        conductivity = _positive(wall_thermal_conductivity_w_m_k, "wall_thermal_conductivity_w_m_k")
    elif wall_material is not None:
        catalog = get_wall_material(wall_material).thermal_conductivity
        if catalog is None:
            raise InputError(f"{wall_material} has no verified catalog conductivity; give "
                             "wall_thermal_conductivity_w_m_k.",
                             code="wall_conductivity_unavailable",
                             field="wall_thermal_conductivity_w_m_k")
        conductivity = catalog.value
    else:
        raise InputError("Give a wall material or a wall conductivity.", field="wall_material")
    if wall_material is not None:
        limit = get_wall_material(wall_material).temperature_limit()[0].value

    stations = list(contour.stations)
    order = list(reversed(range(len(stations)))) if coolant.counterflow else list(range(len(stations)))
    fluid = reference_fluid(coolant.fluid)
    flow_area = width * height
    hydraulic = 2.0 * flow_area / (width + height)
    mass_flux = mass_flow / (channels.count * flow_area)
    temperature, pressure = float(coolant.inlet_temperature_k), float(coolant.inlet_pressure_pa)
    inlet_state = fluid.state_tp(temperature, pressure)
    enthalpy = inlet_state.enthalpy_j_kg
    total_heat = 0.0
    radiation_skipped = 0
    distance = 0.0
    solved: dict[int, ThermalStation] = {}
    warnings: list[WarningMessage] = []
    # Coolant correlations run with extrapolate=True so the march can finish. Their range
    # warnings are kept only for the converged wall temperature (not for bisection probes)
    # and are reported once per out-of-range parameter at the end.
    extrapolated: dict[str, tuple[str, set[int]]] = {}
    recording = False

    def note(result_warnings: Sequence[WarningMessage], station_index: int) -> None:
        if not recording:
            return
        for item in result_warnings:
            if item.code == "correlation_extrapolated":
                _, _, stated = item.message.partition(" lies outside the stated range of ")
                extrapolated.setdefault(str(item.field), (stated.rstrip("."), set()))[1].add(station_index)

    for step, index in enumerate(order):
        station = stations[index]
        circumference_per_channel = 2.0 * pi * (station.radius_m + thickness) / channels.count
        land = circumference_per_channel - width
        if land <= 0.0:
            raise DomainError("Channels do not fit the circumference at this station.",
                              field="width_m", details={"x_m": station.x_m})
        properties = transport_tp(coolant.fluid, temperature, pressure)
        if properties.phase == "liquid" and fluid.model.p_crit_pa > pressure:
            saturation = fluid.saturation_p(pressure)
            if temperature >= saturation.temperature_k * 0.999:
                raise DomainError("Coolant reaches saturation; two-phase flow is not modelled.",
                                  code="coolant_boiling", field="coolant")
        reynolds = mass_flux * hydraulic / properties.viscosity_pa_s
        friction = convection.chen_friction_factor(reynolds, channels.roughness_m / hydraulic).value

        def coolant_coefficient(wall_temperature: float, *, temperature=temperature,
                                pressure=pressure, distance=distance, properties=properties,
                                index=index, reynolds=reynolds, friction=friction) -> float:
            if coolant_correlation == "taylor":
                taylor = convection.taylor_hydrogen_coefficient(
                    fluid=coolant.fluid, bulk_temperature_k=temperature, wall_temperature_k=wall_temperature,
                    pressure_pa=pressure, mass_flux_kg_m2_s=mass_flux, diameter_m=hydraulic,
                    distance_from_entrance_m=max(distance, 50.0 * hydraulic),
                    bulk_viscosity_pa_s=properties.viscosity_pa_s,
                    bulk_conductivity_w_m_k=properties.thermal_conductivity_w_m_k,
                    extrapolate=True)
                note(taylor.warnings, index)
                return taylor.heat_transfer_coefficient_w_m2_k
            gnielinski = convection.gnielinski_nusselt(reynolds, properties.prandtl_number,
                                                       darcy_friction_factor=friction, extrapolate=True)
            note(gnielinski.warnings, index)
            return gnielinski.value * properties.thermal_conductivity_w_m_k / hydraulic

        def radiative(wall_temperature: float, *, station=station) -> float:
            nonlocal radiation_skipped
            if gas_radiation == "none" or (gas.h2o_mole_fraction + gas.co2_mole_fraction) <= 0.0:
                return 0.0
            compression = 1.0 + 0.5 * (gas.gamma - 1.0) * station.mach**2
            gas_temperature = gas.chamber_temperature_k / compression
            static_pressure = gas.chamber_pressure_pa * compression ** (-gas.gamma / (gas.gamma - 1.0))
            beam = radiation.mean_beam_length_for_geometry("infinite_cylinder", 2.0 * station.radius_m)
            try:
                emission = radiation.wsgg_mixture_emissivity(
                    temperature_k=gas_temperature, path_length_m=beam, total_pressure_pa=static_pressure,
                    h2o_mole_fraction=gas.h2o_mole_fraction, co2_mole_fraction=gas.co2_mole_fraction,
                    allow_pressure_extrapolation=True)
                absorption = radiation.wsgg_mixture_absorptivity(
                    gas_temperature_k=gas_temperature, wall_temperature_k=wall_temperature, path_length_m=beam,
                    total_pressure_pa=static_pressure, h2o_mole_fraction=gas.h2o_mole_fraction,
                    co2_mole_fraction=gas.co2_mole_fraction, allow_pressure_extrapolation=True)
            except DomainError:
                radiation_skipped += 1
                return 0.0
            return radiation.gas_to_black_wall_heat_flux(
                gas_emissivity=emission.value, gas_absorptivity=absorption.value,
                gas_temperature_k=gas_temperature, wall_temperature_k=wall_temperature)

        def balance(wall_temperature: float, *, station=station, temperature=temperature,
                    land=land, circumference_per_channel=circumference_per_channel,
                    coolant_coefficient=coolant_coefficient,
                    radiative=radiative) -> tuple[float, float, float, float, float, float]:
            h_gas, recovery = bartz_coefficient(gas, throat_radius_m=contour.throat_radius_m,
                                                area_ratio=station.area_ratio, mach=station.mach,
                                                wall_temperature_k=wall_temperature)
            q_rad = radiative(wall_temperature)
            q_gas = h_gas * (recovery - wall_temperature) + q_rad
            # Taylor's T_w is the coolant-side wall, one conduction drop below the hot wall;
            # bisection probes far below the root are floored at the bulk temperature.
            h_cool = coolant_coefficient(max(wall_temperature - q_gas * thickness / conductivity, temperature))
            fin = conduction.straight_fin(
                thermal_conductivity_w_m_k=conductivity, heat_transfer_coefficient_w_m2_k=h_cool,
                perimeter_m=2.0, cross_section_area_m2=land, length_m=height,
                base_temperature_k=temperature + 1.0, fluid_temperature_k=temperature)
            wetted = width + 2.0 * fin.efficiency * height
            resistance = thickness / conductivity + circumference_per_channel / (h_cool * wetted)
            conducted = (wall_temperature - temperature) / resistance
            return q_gas - conducted, q_gas, h_gas, recovery, h_cool, fin.efficiency

        recovery_guess = bartz_coefficient(gas, throat_radius_m=contour.throat_radius_m,
                                           area_ratio=station.area_ratio, mach=station.mach,
                                           wall_temperature_k=temperature)[1]
        low, high = temperature, recovery_guess
        f_low, f_high = balance(low)[0], balance(high)[0]
        if f_low < 0.0 or f_high > 0.0:
            raise ConvergenceError("Hot-wall temperature is not bracketed.", code="wall_temperature_bracket",
                                   details={"x_m": station.x_m})
        for _ in range(200):
            middle = 0.5 * (low + high)
            if balance(middle)[0] > 0.0:
                low = middle
            else:
                high = middle
            if high - low < 1.0e-9 * middle:
                break
        hot_wall = 0.5 * (low + high)
        recording = True
        _, q_gas, h_gas, recovery, h_cool, efficiency = balance(hot_wall)
        recording = False
        cold_wall = hot_wall - q_gas * thickness / conductivity
        solved[index] = ThermalStation(
            x_m=station.x_m, radius_m=station.radius_m, area_ratio=station.area_ratio, mach=station.mach,
            adiabatic_wall_temperature_k=recovery, gas_side_coefficient_w_m2_k=h_gas,
            radiative_heat_flux_w_m2=q_gas - h_gas * (recovery - hot_wall), heat_flux_w_m2=q_gas,
            hot_wall_temperature_k=hot_wall, cold_wall_temperature_k=cold_wall,
            coolant_temperature_k=temperature, coolant_pressure_pa=pressure,
            coolant_side_coefficient_w_m2_k=h_cool, coolant_reynolds=reynolds,
            land_fin_efficiency=efficiency, wall_limit_margin_k=None if limit is None else limit - hot_wall)
        if step == len(order) - 1:
            break
        following = stations[order[step + 1]]
        segment = sqrt((following.x_m - station.x_m) ** 2 + (following.radius_m - station.radius_m) ** 2)
        heat = q_gas * 2.0 * pi * station.radius_m * segment
        total_heat += heat
        distance += segment
        enthalpy += heat / mass_flow
        pressure -= friction * segment / hydraulic * mass_flux**2 / (2.0 * properties.density_kg_m3)
        if pressure <= 0.0:
            raise DomainError("Coolant pressure drop exceeds the inlet pressure.", code="coolant_pressure_exhausted")
        temperature = _temperature_from_enthalpy(coolant.fluid, enthalpy, pressure, temperature)

    outlet_enthalpy = fluid.state_tp(temperature, pressure).enthalpy_j_kg
    residual = mass_flow * (outlet_enthalpy - inlet_state.enthalpy_j_kg) - total_heat
    ordered = tuple(solved[index] for index in range(len(stations)))
    hottest = max(item.hot_wall_temperature_k for item in ordered)
    warnings.append(WarningMessage(
        "rte_class_model",
        "Quasi-one-dimensional station march; no axial wall conduction, no curvature or entrance "
        "enhancement, frozen chamber gas properties; not validated against hardware data.",
        severity="info"))
    for parameter, (stated, indices) in extrapolated.items():
        warnings.append(WarningMessage(
            "coolant_correlation_extrapolated",
            f"{parameter} lies outside the stated range of {stated} at {len(indices)} of "
            f"{len(stations)} stations; the correlation was extrapolated there.",
            severity="warning", field="coolant"))
    if gas_radiation == "wsgg":
        warnings.append(WarningMessage(
            "wsgg_pressure_extrapolated",
            f"Gas radiation uses 1 atm WSGG fits extrapolated in pressure; {radiation_skipped} station "
            "evaluations were outside the 300-2400 K fit range and got no radiation.",
            severity="critical" if radiation_skipped else "warning"))
    if limit is not None and hottest > limit:
        warnings.append(WarningMessage(
            "wall_limit_exceeded",
            f"Hot-wall temperature {hottest:.0f} K exceeds the {wall_material} limit {limit:.0f} K.",
            severity="critical"))
    return RegenChamberResult(
        stations=ordered,
        total_heat_load_w=total_heat,
        coolant_outlet_temperature_k=temperature,
        coolant_outlet_pressure_pa=pressure,
        coolant_pressure_drop_pa=coolant.inlet_pressure_pa - pressure,
        maximum_hot_wall_temperature_k=hottest,
        wall_limit_temperature_k=limit,
        energy_balance_residual_w=residual,
        warnings=tuple(warnings),
    )


def hot_gas_from_reactants(*, reactants: Mapping[str, float], chamber_pressure_pa: float,
                           throat_radius_of_curvature_m: float) -> HotGas:
    """Build :class:`HotGas` for any reactant mixture from the CEA-gated solver [1]_, [2]_.

    ``reactants`` maps reactant keys (gas species such as ``"H2"`` at
    298.15 K, or condensed reactants such as ``"O2(L)"``) to moles. Frozen
    ``cp`` and ``gamma`` come from the equilibrium chamber composition;
    viscosity and Prandtl number from NASA CEA's transport data at the
    chamber state.

    References
    ----------
    .. [1] NASA RP-1311, "Computer Program for Calculation of Complex Chemical
           Equilibrium Compositions and Applications, Part I: Analysis".
           https://ntrs.nasa.gov/citations/19950013764
    .. [2] NASA-CEA-TRANS. NASA Glenn Research Center, CEA transport database
           ``data/trans.inp``. https://github.com/nasa/cea
    """

    from rocket_propulsion.thermochemistry import (
        EquilibriumProblem,
        LocalEquilibriumProvider,
        solve_rocket_performance,
    )
    from rocket_propulsion.thermochemistry.transport import equilibrium_transport_properties

    mixture = dict(reactants)
    performance = solve_rocket_performance(reactants=mixture, chamber_pressure_pa=chamber_pressure_pa,
                                           area_ratios=(2.0,))
    chamber = LocalEquilibriumProvider().solve(EquilibriumProblem.hp(mixture, pressure_pa=chamber_pressure_pa))
    transport = equilibrium_transport_properties(chamber)
    fractions = dict(performance.chamber.mole_fractions)
    return HotGas(
        chamber_pressure_pa=chamber_pressure_pa,
        chamber_temperature_k=chamber.temperature_k,
        characteristic_velocity_m_s=performance.characteristic_velocity_m_s,
        specific_heat_j_kg_k=chamber.cp_j_kg_k,
        viscosity_pa_s=transport.viscosity_pa_s,
        prandtl_number=transport.frozen_prandtl,
        gamma=chamber.gamma,
        throat_radius_of_curvature_m=throat_radius_of_curvature_m,
        h2o_mole_fraction=fractions.get("H2O", 0.0),
        co2_mole_fraction=fractions.get("CO2", 0.0),
    )


def hot_gas_from_propellant(*, propellant_pair: str, chamber_pressure_pa: float,
                            oxidizer_fuel_ratio: float | None = None,
                            throat_radius_of_curvature_m: float) -> HotGas:
    """Build :class:`HotGas` for a propellant preset; see :func:`hot_gas_from_reactants` [1]_, [2]_.

    References
    ----------
    .. [1] NASA RP-1311, "Computer Program for Calculation of Complex Chemical
           Equilibrium Compositions and Applications, Part I: Analysis".
           https://ntrs.nasa.gov/citations/19950013764
    .. [2] NASA-CEA-TRANS. NASA Glenn Research Center, CEA transport database
           ``data/trans.inp``. https://github.com/nasa/cea
    """

    from rocket_propulsion.thermochemistry import (
        get_propellant_pair,
        reactants_for_mass_ratio,
    )

    pair = get_propellant_pair(propellant_pair)
    ratio = pair.reference_oxidizer_fuel_ratio if oxidizer_fuel_ratio is None else oxidizer_fuel_ratio
    return hot_gas_from_reactants(reactants=reactants_for_mass_ratio(pair, ratio),
                                  chamber_pressure_pa=chamber_pressure_pa,
                                  throat_radius_of_curvature_m=throat_radius_of_curvature_m)


__all__: Sequence[str] = (
    "Coolant",
    "CoolingChannels",
    "HotGas",
    "RegenChamberResult",
    "ThermalStation",
    "bartz_coefficient",
    "hot_gas_from_propellant",
    "hot_gas_from_reactants",
    "solve_regen_chamber",
)
