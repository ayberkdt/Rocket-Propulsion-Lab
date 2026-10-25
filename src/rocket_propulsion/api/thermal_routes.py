"""HTTP routes for heat transfer, regenerative cooling and combustion products/losses.

Handlers translate JSON to library calls; physics lives in
:mod:`rocket_propulsion.heat_transfer` and :mod:`rocket_propulsion.thermochemistry`,
which document their equations and sources. Responses are strict JSON.

References
----------
.. [1] LIENHARD-AHTT-6.00. J. H. Lienhard IV and J. H. Lienhard V, "A Heat
       Transfer Textbook", v6.00, 2024.
.. [2] NASA RP-1311, "Computer Program for Calculation of Complex Chemical
       Equilibrium Compositions and Applications, Part I: Analysis".
       https://ntrs.nasa.gov/citations/19950013764
"""

from __future__ import annotations

from dataclasses import asdict
from math import pi, sqrt
from typing import Any

from rocket_propulsion.api.aero_routes import _clean, _integer, _linspace, _number, _text
from rocket_propulsion.core.errors import DomainError, InputError
from rocket_propulsion.core.metadata import ModelMetadata
from rocket_propulsion.geometry import generate_nozzle_contour
from rocket_propulsion.heat_transfer import conduction, convection, radiation
from rocket_propulsion.heat_transfer.regen_chamber import (
    Coolant,
    CoolingChannels,
    hot_gas_from_propellant,
    solve_regen_chamber,
)
from rocket_propulsion.thermochemistry import combustion_efficiency as efficiency
from rocket_propulsion.thermochemistry import (
    get_propellant_pair,
    reactants_for_mass_ratio,
    solve_rocket_performance,
)
from rocket_propulsion.thermochemistry.finite_area_combustor import solve_finite_area_combustor
from rocket_propulsion.thermochemistry.nozzle_kinetics import (
    HyperbolicNozzle,
    solve_nozzle_kinetics,
)
from rocket_propulsion.thermodynamics.real_fluids.transport import transport_tp


def heat_conduction_route(payload: dict[str, Any]) -> dict[str, Any]:
    """Transient two-layer plane wall heated by hot gas and cooled by a coolant [1]_.

    Parameters
    ----------
    payload : dict
        hot_wall_thickness_m, hot_wall_conductivity_w_m_k : float, default 1 mm, 305
        hot_wall_heat_capacity_j_m3_k : float, default 3.5e6
        gas_heat_transfer_coefficient_w_m2_k, gas_temperature_k : float, default 2e4, 3000
        coolant_heat_transfer_coefficient_w_m2_k, coolant_temperature_k : float, default 5e4, 150
        duration_s : float, default 0.2

    Returns
    -------
    dict
        Surface temperature histories, temperature profiles at five times and
        the energy audit.

    References
    ----------
    .. [1] LIENHARD-AHTT-6.00. J. H. Lienhard IV and J. H. Lienhard V, "A Heat
           Transfer Textbook", v6.00, 2024, Ch. 5 (transient conduction).
    """

    thickness = _number(payload, "hot_wall_thickness_m", 1.0e-3)
    layer = conduction.WallLayer(
        thickness_m=thickness, thermal_conductivity=_number(payload, "hot_wall_conductivity_w_m_k", 305.0),
        volumetric_heat_capacity=_number(payload, "hot_wall_heat_capacity_j_m3_k", 3.5e6))
    duration = _number(payload, "duration_s", 0.2)
    if duration <= 0.0:
        raise DomainError("Duration must be positive.", field="duration_s")
    steps = 200
    cells = 40
    result = conduction.solve_transient_wall(
        [layer], conduction.WallGeometry(kind="plane"),
        conduction.SurfaceConvection(heat_transfer_coefficient_w_m2_k=_number(payload, "coolant_heat_transfer_coefficient_w_m2_k", 5.0e4),
                                     fluid_temperature_k=_number(payload, "coolant_temperature_k", 150.0)),
        conduction.SurfaceConvection(heat_transfer_coefficient_w_m2_k=_number(payload, "gas_heat_transfer_coefficient_w_m2_k", 2.0e4),
                                     fluid_temperature_k=_number(payload, "gas_temperature_k", 3000.0)),
        initial_temperature_k=_number(payload, "initial_temperature_k", 150.0),
        time_step_s=duration / steps, step_count=steps, cells_per_layer=cells, snapshot_every=steps // 4)
    history = [{"time_s": t, "hot_surface_k": hot, "cold_surface_k": cold}
               for t, hot, cold in zip(result.times_s, result.inner_surface_temperatures_k,
                                       result.outer_surface_temperatures_k)]
    profiles = [{"time_s": t, "points": [{"x_mm": 1000.0 * x, "temperature_k": value}
                                         for x, value in zip(result.cell_centers_m, temperatures)]}
                for t, temperatures in zip(result.snapshot_times_s, result.snapshot_temperatures_k)]
    return _clean({
        "history": history, "profiles": profiles,
        "final_hot_surface_k": result.inner_surface_temperatures_k[-1],
        "final_cold_surface_k": result.outer_surface_temperatures_k[-1],
        "stored_energy_change_j_m2": result.stored_energy_change_j,
        "energy_residual_j_m2": result.energy_residual_j, "warnings": result.warnings,
    })


def heat_convection_route(payload: dict[str, Any]) -> dict[str, Any]:
    """Coolant-side film coefficients from several correlations over a flow sweep [1]_.

    Parameters
    ----------
    payload : dict
        fluid : str, default "methane"
        temperature_k, pressure_pa : float, default 200 K, 15 MPa
        hydraulic_diameter_m : float, default 2 mm
        mass_flux_kg_m2_s : float, default 2e4 (centre of the sweep)

    Returns
    -------
    dict
        Coolant transport state and a mass-flux sweep of h from Gnielinski,
        Dittus-Boelter and Colburn.

    References
    ----------
    .. [1] LIENHARD-AHTT-6.00. J. H. Lienhard IV and J. H. Lienhard V, "A Heat
           Transfer Textbook", v6.00, 2024, Ch. 7 (internal forced convection).
    """

    fluid = _text(payload, "fluid", "methane")
    diameter = _number(payload, "hydraulic_diameter_m", 2.0e-3)
    centre = _number(payload, "mass_flux_kg_m2_s", 2.0e4)
    properties = transport_tp(fluid, _number(payload, "temperature_k", 200.0), _number(payload, "pressure_pa", 1.5e7))
    sweep = []
    for flux in _linspace(0.25 * centre, 2.0 * centre, 15):
        reynolds = flux * diameter / properties.viscosity_pa_s
        row: dict[str, Any] = {"mass_flux_kg_m2_s": flux, "reynolds": reynolds}
        scale = properties.thermal_conductivity_w_m_k / diameter
        for key, call in (("gnielinski", lambda reynolds=reynolds: convection.gnielinski_nusselt(reynolds, properties.prandtl_number, extrapolate=True)),
                          ("dittus_boelter", lambda reynolds=reynolds: convection.power_law_nusselt(reynolds, properties.prandtl_number, "dittus_boelter", extrapolate=True)),
                          ("colburn", lambda reynolds=reynolds: convection.power_law_nusselt(reynolds, properties.prandtl_number, "colburn", extrapolate=True))):
            try:
                row[key] = call().value * scale
            except (DomainError, InputError):
                row[key] = None
        sweep.append(row)
    return _clean({"transport": properties, "sweep": sweep})


def heat_radiation_route(payload: dict[str, Any]) -> dict[str, Any]:
    """Radiative-equilibrium temperature of a radiation-cooled conical nozzle extension [1]_.

    Parameters
    ----------
    payload : dict
        inlet_radius_m, exit_radius_m, length_m : float, default 0.2, 0.6, 1.0
        heat_transfer_coefficient_w_m2_k, recovery_temperature_k : float, default 200, 2000
        emissivity : float, default 0.8

    Returns
    -------
    dict
        Frustum view factors, the wall temperature at the chosen emissivity and an
        emissivity sweep.

    References
    ----------
    .. [1] LIENHARD-AHTT-6.00. J. H. Lienhard IV and J. H. Lienhard V, "A Heat
           Transfer Textbook", v6.00, 2024, Ch. 10 (radiation exchange).
    """

    inlet = _number(payload, "inlet_radius_m", 0.2)
    exit_radius = _number(payload, "exit_radius_m", 0.6)
    length = _number(payload, "length_m", 1.0)
    h = _number(payload, "heat_transfer_coefficient_w_m2_k", 200.0)
    recovery = _number(payload, "recovery_temperature_k", 2000.0)
    emissivity = _number(payload, "emissivity", 0.8)
    enclosure = radiation.frustum_enclosure_view_factors(small_radius_m=inlet, large_radius_m=exit_radius, length_m=length)

    def wall_temperature(eps: float) -> float:
        inner = radiation.frustum_interior_transfer_factors(small_radius_m=inlet, large_radius_m=exit_radius,
                                                            length_m=length, wall_emissivity=eps)
        sinks = [radiation.RadiativeSink(transfer_factor=factor, temperature_k=3.0) for factor in inner]
        sinks.append(radiation.RadiativeSink(transfer_factor=eps, temperature_k=3.0))
        return radiation.radiative_equilibrium_wall_temperature(
            sinks=sinks, heat_transfer_coefficient_w_m2_k=h, recovery_temperature_k=recovery).wall_temperature_k

    sweep = [{"emissivity": eps, "wall_temperature_k": wall_temperature(eps)} for eps in _linspace(0.1, 1.0, 10)]
    return _clean({"view_factors": enclosure, "emissivity": emissivity,
                   "wall_temperature_k": wall_temperature(emissivity), "sweep": sweep})


def regen_chamber_route(payload: dict[str, Any]) -> dict[str, Any]:
    """Coupled regenerative cooling of a conical-nozzle thrust chamber [1]_.

    Parameters
    ----------
    payload : dict
        propellant_pair : str, default "LOX_CH4"
        chamber_pressure_pa, throat_area_m2, area_ratio : float, default 7 MPa, 0.01, 10
        coolant : str, default "methane"; coolant_inlet_temperature_k, coolant_inlet_pressure_pa
        coolant_mass_flow_kg_s : float, optional (default: the fuel flow)
        channel_count, channel_width_m, channel_height_m, hot_wall_thickness_m
        wall_material : str, default "grcop-84"

    Returns
    -------
    dict
        Station profiles (temperatures, heat flux, coolant state) and totals.

    References
    ----------
    .. [1] NASA-TM-1998-206313-WADEL. M. F. Wadel, NASA/TM-1998-206313
           (RTE station-marching coolant model).
    """

    pair = get_propellant_pair(_text(payload, "propellant_pair", "LOX_CH4"))
    pressure = _number(payload, "chamber_pressure_pa", 7.0e6)
    throat_area = _number(payload, "throat_area_m2", 0.01)
    area_ratio = _number(payload, "area_ratio", 10.0)
    gas = hot_gas_from_propellant(propellant_pair=pair.key, chamber_pressure_pa=pressure,
                                  throat_radius_of_curvature_m=sqrt(throat_area / pi))
    contour = generate_nozzle_contour(contour="conical", throat_area_m2=throat_area, area_ratio=area_ratio,
                                      gamma=gas.gamma, chamber_radius_ratio=2.0, station_count=41)
    total_flow = pressure * throat_area / gas.characteristic_velocity_m_s
    fuel_flow = total_flow / (1.0 + pair.reference_oxidizer_fuel_ratio)
    coolant_flow = payload.get("coolant_mass_flow_kg_s")
    flow = fuel_flow if coolant_flow in (None, "") else _number(payload, "coolant_mass_flow_kg_s", fuel_flow)
    result = solve_regen_chamber(
        contour=contour, gas=gas,
        channels=CoolingChannels(count=_integer(payload, "channel_count", 100, 4, 1000),
                                 width_m=_number(payload, "channel_width_m", 1.5e-3),
                                 height_m=_number(payload, "channel_height_m", 3.0e-3),
                                 hot_wall_thickness_m=_number(payload, "hot_wall_thickness_m", 1.0e-3),
                                 roughness_m=_number(payload, "roughness_m", 1.0e-6)),
        coolant=Coolant(fluid=_text(payload, "coolant", "methane"),
                        inlet_temperature_k=_number(payload, "coolant_inlet_temperature_k", 120.0),
                        inlet_pressure_pa=_number(payload, "coolant_inlet_pressure_pa", 1.5e7),
                        mass_flow_kg_s=flow),
        wall_material=_text(payload, "wall_material", "grcop-84"))
    return _clean({"gas": gas, "fuel_flow_kg_s": fuel_flow, "coolant_flow_kg_s": flow, **asdict(result)})


def combustion_products_route(payload: dict[str, Any]) -> dict[str, Any]:
    """Equilibrium product composition, radicals and dissociation along the nozzle [1]_.

    Parameters
    ----------
    payload : dict
        propellant_pair : str, default "LOX_LH2"
        oxidizer_fuel_ratio : float, optional (preset reference)
        chamber_pressure_pa, area_ratio : float, default 7 MPa, 40

    Returns
    -------
    dict
        Stoichiometry, per-station leading species and radical fractions, and
        the dissociation summary.

    References
    ----------
    .. [1] NASA RP-1311, "Computer Program for Calculation of Complex Chemical
           Equilibrium Compositions and Applications, Part I: Analysis".
           https://ntrs.nasa.gov/citations/19950013764
    """

    pair = get_propellant_pair(_text(payload, "propellant_pair", "LOX_LH2"))
    ratio = _number(payload, "oxidizer_fuel_ratio", pair.reference_oxidizer_fuel_ratio)
    report = efficiency.combustion_product_report(
        chamber_pressure_pa=_number(payload, "chamber_pressure_pa", 7.0e6),
        area_ratio=_number(payload, "area_ratio", 40.0), propellant_pair=pair.key, oxidizer_fuel_ratio=ratio)
    stations = []
    for station in report.stations:
        fractions = sorted(station.mole_fractions, key=lambda item: -item[1])
        stations.append({"name": station.name, "temperature_k": station.temperature_k,
                         "pressure_pa": station.pressure_pa,
                         "molecular_mass_g_mol": station.molecular_mass_kg_mol * 1000.0,
                         "radical_mole_fraction": station.radical_mole_fraction,
                         "leading_species": [{"species": key, "mole_fraction": value}
                                             for key, value in fractions if value >= 1.0e-4][:10],
                         "radicals": [{"species": key, "mole_fraction": value}
                                      for key, value in station.radicals if value > 0.0]})
    return _clean({"pair": pair.key, "oxidizer_fuel_ratio": ratio, "stoichiometry": report.stoichiometry,
                   "characteristic_velocity_m_s": report.characteristic_velocity_m_s,
                   "vacuum_specific_impulse_s": report.vacuum_specific_impulse_s,
                   "dissociation": report.dissociation, "stations": stations, "warnings": report.warnings})


def combustion_losses_route(payload: dict[str, Any]) -> dict[str, Any]:
    """Performance loss chain from ideal equilibrium to delivered impulse [1]_, [2]_.

    Parameters
    ----------
    payload : dict
        propellant_pair : str, default "LOX_LH2"
        chamber_pressure_pa, area_ratio : float, default 6.894757 MPa, 40
        contraction_ratio : float, default 2.5 (finite-area combustor)
        divergence_efficiency, characteristic_velocity_efficiency : float, default 0.985, 0.98
        boundary_layer_decrement_s : float, default 3
        compute_kinetics : bool, default true for LOX_LH2 (H2-O2 mechanism only)

    Returns
    -------
    dict
        Frozen/shifting limits, the finite-area combustor, the finite-rate
        result when available and the JANNAF-style loss chain.

    References
    ----------
    .. [1] NASA RP-1311, "Computer Program for Calculation of Complex Chemical
           Equilibrium Compositions and Applications, Part I: Analysis".
           https://ntrs.nasa.gov/citations/19950013764
    .. [2] NASA-RP-1328. K. Radhakrishnan, NASA RP-1328, 1994 (one-dimensional
           kinetics equations). https://ntrs.nasa.gov/citations/19950008642
    """

    pair = get_propellant_pair(_text(payload, "propellant_pair", "LOX_LH2"))
    ratio = _number(payload, "oxidizer_fuel_ratio", pair.reference_oxidizer_fuel_ratio)
    pressure = _number(payload, "chamber_pressure_pa", 6.894757e6)
    area_ratio = _number(payload, "area_ratio", 40.0)
    reactants = reactants_for_mass_ratio(pair, ratio)
    limits = efficiency.frozen_equilibrium_limits(propellant_pair=pair.key, oxidizer_fuel_ratio=ratio,
                                                  chamber_pressure_pa=pressure, area_ratio=area_ratio)
    combustor = solve_finite_area_combustor(reactants=reactants, injector_pressure_pa=pressure,
                                            contraction_ratio=_number(payload, "contraction_ratio", 2.5),
                                            area_ratios=(area_ratio,))
    kinetic_efficiency = 1.0
    kinetics = None
    raw = payload.get("compute_kinetics", pair.key == "LOX_LH2")
    compute = raw if isinstance(raw, bool) else str(raw).strip().lower() in {"true", "1", "on", "yes"}
    if compute:
        if pair.key != "LOX_LH2":
            raise InputError("Finite-rate kinetics is available for LOX_LH2 only (H2-O2 mechanism).",
                             field="compute_kinetics")
        throat_radius = _number(payload, "throat_radius_m", 0.05)
        solution = solve_nozzle_kinetics(reactants=reactants, chamber_pressure_pa=pressure,
                                         nozzle=HyperbolicNozzle(throat_radius, 3.0, area_ratio, 30.0, 15.0))
        kinetic_efficiency = solution.vacuum_specific_impulse_s / limits.equilibrium_specific_impulse_s
        kinetics = {"vacuum_specific_impulse_s": solution.vacuum_specific_impulse_s,
                    "stations": [{"area_ratio": s.area_ratio, "temperature_k": s.temperature_k,
                                  **{f"x_{k}": v for k, v in s.mole_fractions if k in {"H", "O", "OH", "H2O"}}}
                                 for s in solution.stations[:: max(1, len(solution.stations) // 120)]],
                    "warnings": solution.warnings}
    chain = efficiency.jannaf_isp_loss_chain(
        ode_specific_impulse_s=limits.equilibrium_specific_impulse_s, kinetic_efficiency=kinetic_efficiency,
        divergence_efficiency=_number(payload, "divergence_efficiency", 0.985),
        boundary_layer_decrement_s=_number(payload, "boundary_layer_decrement_s", 3.0),
        characteristic_velocity_efficiency=_number(payload, "characteristic_velocity_efficiency", 0.98))
    ideal = solve_rocket_performance(reactants=reactants, chamber_pressure_pa=pressure, area_ratios=(area_ratio,))
    return _clean({
        "limits": limits, "kinetics": kinetics,
        "finite_area_combustor": {
            "contraction_ratio": combustor.contraction_ratio,
            "stagnation_pressure_ratio": combustor.stagnation_pressure_ratio,
            "vacuum_specific_impulse_s": combustor.exits[0].vacuum_specific_impulse_s,
            "infinite_area_vacuum_specific_impulse_s": ideal.exits[0].vacuum_specific_impulse_s,
            "combustor_end_mach_number": combustor.combustor_end_mach_number},
        "chain": {key: value for key, value in asdict(chain).items() if key != "loss_budget"},
    })


THERMAL_ROUTES = {
    "/api/v1/heat-conduction": heat_conduction_route,
    "/api/v1/heat-convection": heat_convection_route,
    "/api/v1/heat-radiation": heat_radiation_route,
    "/api/v1/regen-chamber": regen_chamber_route,
    "/api/v1/combustion-products": combustion_products_route,
    "/api/v1/combustion-losses": combustion_losses_route,
}

_LIENHARD = ("Lienhard & Lienhard, A Heat Transfer Textbook, v6.00 (2024)",)
THERMAL_ROUTE_METADATA = {
    "/api/v1/heat-conduction": ModelMetadata(
        model="finite-volume transient conduction (theta method)", version="1.0",
        assumptions=("One-dimensional plane wall", "Constant boundary coefficients during the run"),
        units={"temperature": "K", "length": "mm", "time": "s"}, references=_LIENHARD),
    "/api/v1/heat-convection": ModelMetadata(
        model="internal forced-convection correlations with real-fluid properties", version="1.0",
        assumptions=("Fully developed turbulent pipe flow", "Bulk properties; no property-ratio correction"),
        units={"heat_transfer_coefficient": "W/(m^2 K)"}, references=_LIENHARD),
    "/api/v1/heat-radiation": ModelMetadata(
        model="gray diffuse enclosure and radiative equilibrium", version="1.0",
        assumptions=("Gray diffuse isothermal wall", "Openings radiate to 3 K space"),
        units={"temperature": "K"}, references=_LIENHARD + ("Howell configuration-factor catalog",)),
    "/api/v1/regen-chamber": ModelMetadata(
        model="RTE-class coupled regenerative cooling march", version="1.0",
        assumptions=("Quasi-1D station march", "Bartz gas side with local wall temperature",
                     "Gnielinski coolant side with NIST real-fluid properties"),
        units={"temperature": "K", "heat_flux": "W/m^2", "pressure": "Pa"},
        references=("NASA/TM-1998-206313", "NASA SP-125")),
    "/api/v1/combustion-products": ModelMetadata(
        model="equilibrium product composition and dissociation", version="1.0",
        assumptions=("Shifting equilibrium, gas phase", "Gated against NASA CEA 3.3.4"),
        units={"mole_fraction": "1", "temperature": "K"}, references=("NASA RP-1311",)),
    "/api/v1/combustion-losses": ModelMetadata(
        model="JANNAF-style impulse loss chain with FAC and finite-rate kinetics", version="1.0",
        assumptions=("Divergence, boundary-layer and c* efficiencies are inputs",
                     "Finite-rate kinetics: GRI-Mech 3.0 H2-O2 subset"),
        units={"specific_impulse": "s"}, references=("NASA RP-1311", "NASA RP-1328", "NASA TP-1999-208522")),
}
