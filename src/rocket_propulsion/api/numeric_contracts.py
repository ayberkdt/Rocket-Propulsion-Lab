"""Declared top-level numeric request fields, including shared API helpers.

This type-admission catalogue does not describe required/default/range
rules or nested object/array contracts. Runtime validation uses these
immutable declarations; tests independently audit reader sites and exercise
the dispatcher before model execution.
"""

from types import MappingProxyType

NUMBER_FIELDS = MappingProxyType(
    {
        "/api/v1/airfoil": frozenset({"angle_of_attack_deg", "panels_per_surface"}),
        "/api/v1/burn/export/sidera": frozenset({"t_start_s"}),
        "/api/v1/burn/report": frozenset(),
        "/api/v1/burn/simulate": frozenset(),
        "/api/v1/case-comparison": frozenset(),
        "/api/v1/chemistry-capabilities": frozenset(),
        "/api/v1/combustion-losses": frozenset(
            {
                "area_ratio",
                "boundary_layer_decrement_s",
                "chamber_pressure_pa",
                "characteristic_velocity_efficiency",
                "contraction_ratio",
                "divergence_efficiency",
                "oxidizer_fuel_ratio",
                "throat_radius_m",
            }
        ),
        "/api/v1/combustion-products": frozenset(
            {"area_ratio", "chamber_pressure_pa", "oxidizer_fuel_ratio"}
        ),
        "/api/v1/combustion-stability": frozenset(
            {
                "chamber_length_m",
                "chamber_pressure_pa",
                "chamber_radius_m",
                "combustion_time_lag_s",
                "contraction_ratio",
                "fuel_drop_fraction",
                "oxidizer_drop_fraction",
                "oxidizer_fuel_ratio",
            }
        ),
        "/api/v1/duct-flow": frozenset({"gamma", "mach"}),
        "/api/v1/engine-cycle": frozenset(
            {
                "chamber_mass_flow_kg_s",
                "chamber_pressure_pa",
                "injector_pressure_drop_fraction",
                "inlet_pressure_pa",
                "propellant_density_kg_m3",
                "pump_efficiency",
                "turbine_efficiency",
                "turbine_gamma",
                "turbine_inlet_temperature_k",
                "turbine_pressure_ratio",
                "turbine_specific_heat_j_kg_k",
            }
        ),
        "/api/v1/engine-operating-point": frozenset(
            {
                "active_engine_count",
                "ambient_pressure_pa",
                "chamber_pressure_pa",
                "chamber_specific_impulse_s",
                "commanded_throttle",
                "delivered_thrust_n",
                "dump_mass_flow_kg_s",
                "ideal_thrust_n",
                "oxidizer_fuel_ratio",
                "realized_throttle",
                "system_specific_impulse_s",
            }
        ),
        "/api/v1/entry": frozenset(
            {
                "drag_coefficient",
                "entry_altitude_m",
                "entry_flight_path_angle_deg",
                "entry_velocity_m_s",
                "lift_to_drag",
                "mass_kg",
                "nose_radius_m",
                "reference_area_m2",
            }
        ),
        "/api/v1/equilibrium": frozenset(
            {
                "pressure_pa",
                "reactant_temperature_k",
                "target_enthalpy_j",
                "target_internal_energy_j",
                "temperature_k",
                "volume_m3",
            }
        ),
        "/api/v1/flat-plate": frozenset({"altitude_m", "length_m", "mach", "wall_temperature_k"}),
        "/api/v1/heat-conduction": frozenset(
            {
                "coolant_heat_transfer_coefficient_w_m2_k",
                "coolant_temperature_k",
                "duration_s",
                "gas_heat_transfer_coefficient_w_m2_k",
                "gas_temperature_k",
                "hot_wall_conductivity_w_m_k",
                "hot_wall_heat_capacity_j_m3_k",
                "hot_wall_thickness_m",
                "initial_temperature_k",
            }
        ),
        "/api/v1/heat-convection": frozenset(
            {"hydraulic_diameter_m", "mass_flux_kg_m2_s", "pressure_pa", "temperature_k"}
        ),
        "/api/v1/heat-radiation": frozenset(
            {
                "emissivity",
                "exit_radius_m",
                "heat_transfer_coefficient_w_m2_k",
                "inlet_radius_m",
                "length_m",
                "recovery_temperature_k",
            }
        ),
        "/api/v1/ideal-gas": frozenset(
            {"density_kg_m3", "gamma", "gas_constant_j_kg_k", "pressure_pa", "temperature_k"}
        ),
        "/api/v1/integrations/sidera/capabilities": frozenset(),
        "/api/v1/isentropic": frozenset({"gamma", "value"}),
        "/api/v1/isentropic-sweep": frozenset(
            {"gamma", "maximum_mach", "minimum_mach", "point_count"}
        ),
        "/api/v1/launch-vehicle-aero": frozenset(
            {
                "altitude_m",
                "angle_of_attack_deg",
                "body_length_m",
                "crossflow_drag_coefficient",
                "diameter_m",
                "mach",
                "nose_length_m",
            }
        ),
        "/api/v1/loss-budget": frozenset(
            {
                "boundary_layer_efficiency",
                "combustion_efficiency",
                "discharge_coefficient",
                "divergence_efficiency",
                "ideal_specific_impulse_s",
                "ideal_thrust_n",
                "relative_uncertainty",
            }
        ),
        "/api/v1/normal-shock": frozenset({"gamma", "upstream_mach"}),
        "/api/v1/nozzle": frozenset(
            {
                "ambient_pressure_pa",
                "area_ratio",
                "chamber_pressure_pa",
                "chamber_temperature_k",
                "gamma",
                "gas_constant_j_kg_k",
                "throat_area_m2",
            }
        ),
        "/api/v1/nozzle-chemistry": frozenset(
            {
                "ambient_pressure_pa",
                "area_ratio",
                "chamber_pressure_pa",
                "maximum_altitude_m",
                "oxidizer_fuel_ratio",
                "point_count",
                "throat_area_m2",
            }
        ),
        "/api/v1/nozzle-contour": frozenset(
            {
                "angular_segments",
                "area_ratio",
                "bell_length_fraction",
                "chamber_radius_ratio",
                "characteristic_count",
                "exit_angle_deg",
                "gamma",
                "half_angle_deg",
                "station_count",
                "throat_angle_deg",
                "throat_area_m2",
            }
        ),
        "/api/v1/nozzle-off-design": frozenset(
            {
                "altitude_point_count",
                "ambient_pressure_pa",
                "area_ratio",
                "chamber_pressure_pa",
                "chamber_temperature_k",
                "gamma",
                "gas_constant_j_kg_k",
                "maximum_altitude_m",
                "station_count",
                "throat_area_m2",
            }
        ),
        "/api/v1/nozzle-thermal": frozenset(
            {
                "area_ratio",
                "bell_length_fraction",
                "chamber_pressure_pa",
                "chamber_temperature_k",
                "channel_count",
                "channel_flow_area_m2",
                "channel_length_m",
                "characteristic_count",
                "characteristic_velocity_m_s",
                "coolant_density_kg_m3",
                "coolant_dynamic_viscosity_pa_s",
                "coolant_inlet_temperature_k",
                "coolant_mass_flow_kg_s",
                "coolant_side_heat_transfer_coefficient_w_m2_k",
                "coolant_specific_heat_j_kg_k",
                "dynamic_viscosity_pa_s",
                "gamma",
                "half_angle_deg",
                "hydraulic_diameter_m",
                "laminar_darcy_poiseuille_number",
                "maximum_coolant_temperature_k",
                "prandtl_number",
                "roughness_m",
                "specific_heat_j_kg_k",
                "station_count",
                "throat_area_m2",
                "throat_radius_of_curvature_m",
                "wall_temperature_k",
                "wall_thickness_m",
            }
        ),
        "/api/v1/oblique-shock": frozenset({"gamma", "turning_angle_deg", "upstream_mach"}),
        "/api/v1/of-sweep": frozenset(
            {"area_ratio", "chamber_pressure_pa", "maximum_ratio", "minimum_ratio", "point_count"}
        ),
        "/api/v1/polytropic": frozenset(
            {
                "exponent",
                "gamma",
                "gas_constant_j_kg_k",
                "inlet_pressure_pa",
                "inlet_temperature_k",
                "outlet_pressure_pa",
                "point_count",
            }
        ),
        "/api/v1/propellant-vaporization": frozenset(
            {
                "chamber_length_m",
                "chamber_pressure_pa",
                "contraction_ratio",
                "drop_radius_m",
                "injection_velocity_m_s",
            }
        ),
        "/api/v1/propellants": frozenset(),
        "/api/v1/real-fluid": frozenset({"point_count", "pressure_pa", "temperature_k"}),
        "/api/v1/regen-chamber": frozenset(
            {
                "area_ratio",
                "chamber_pressure_pa",
                "channel_count",
                "channel_height_m",
                "channel_helix_angle_deg",
                "channel_width_m",
                "closeout_thickness_m",
                "coolant_inlet_pressure_pa",
                "coolant_inlet_temperature_k",
                "coolant_mass_flow_kg_s",
                "hot_wall_thickness_m",
                "roughness_m",
                "throat_area_m2",
            }
        ),
        "/api/v1/rocket-equation": frozenset(
            {
                "delta_v_m_s",
                "effective_exhaust_velocity_m_s",
                "final_mass_kg",
                "initial_mass_kg",
                "specific_impulse_s",
            }
        ),
        "/api/v1/rocket-performance": frozenset({"chamber_pressure_pa", "oxidizer_fuel_ratio"}),
        "/api/v1/sass/mission": frozenset(
            {"bus_dry_mass_kg", "fixed_reserve_kg", "residual_fraction"}
        ),
        "/api/v1/species": frozenset(),
        "/api/v1/staging": frozenset({"payload_mass_kg", "resolution", "target_delta_v_m_s"}),
        "/api/v1/stagnation": frozenset(
            {"gamma", "gas_constant_j_kg_k", "mach", "static_pressure_pa", "static_temperature_k"}
        ),
        "/api/v1/supersonic-body": frozenset(
            {"airfoil_half_angle_deg", "angle_of_attack_deg", "cone_half_angle_deg", "mach"}
        ),
        "/api/v1/thermal-properties": frozenset(
            {
                "maximum_temperature_k",
                "minimum_temperature_k",
                "point_count",
                "pressure_pa",
                "temperature_k",
            }
        ),
        "/api/v1/thermo-process": frozenset(
            {
                "efficiency",
                "gamma",
                "gas_constant_j_kg_k",
                "inlet_pressure_pa",
                "inlet_temperature_k",
                "outlet_pressure_pa",
            }
        ),
        "/api/v1/thrust": frozenset(
            {
                "ambient_pressure_pa",
                "burn_time_s",
                "exit_area_m2",
                "exit_pressure_pa",
                "exit_velocity_m_s",
                "mass_flow_kg_s",
            }
        ),
        "/api/v1/transient/blowdown": frozenset(
            {
                "cant_efficiency",
                "chamber_pressure_exponent",
                "chamber_pressure_pa",
                "delivered_thrust_n",
                "initial_pressure_pa",
                "initial_propellant_mass_kg",
                "initial_ullage_volume_m3",
                "max_points",
                "maximum_mass_step_kg",
                "maximum_refinements",
                "maximum_supply_pressure_pa",
                "minimum_supply_pressure_pa",
                "polytropic_exponent",
                "propellant_density_kg_m3",
                "reference_supply_pressure_pa",
                "relative_tolerance",
                "reserve_mass_kg",
                "specific_impulse_pressure_exponent",
                "system_specific_impulse_s",
                "thrust_pressure_exponent",
            }
        ),
        "/api/v1/transient/calibration": frozenset(
            {"coverage_factor", "maximum_iterations", "relative_tolerance"}
        ),
        "/api/v1/transient/electric": frozenset(
            {
                "available_power_w",
                "duty_cycle",
                "power_margin_w",
                "propellant_mass_kg",
                "propellant_reserve_kg",
                "sweep_maximum_power_w",
                "sweep_minimum_power_w",
                "sweep_points",
                "up_switch_hysteresis_w",
                "vehicle_mass_kg",
            }
        ),
        "/api/v1/transient/feed-line": frozenset(
            {
                "ambient_pressure_pa",
                "bottle_gas_temperature_k",
                "bottle_heat_transfer_w_k",
                "bottle_pressure_pa",
                "bottle_volume_m3",
                "bottle_wall_temperature_k",
                "feed_line_resistance_pa_s2_kg2",
                "line_compliance_kg_pa",
                "line_inertance_pa_s2_kg",
                "line_mass_flow_kg_s",
                "line_resistance_pa_s2_kg2",
                "manifold_pressure_pa",
                "maximum_absolute_flow_kg_s",
                "maximum_manifold_pressure_pa",
                "maximum_tank_pressure_pa",
                "minimum_bottle_pressure_margin_pa",
                "minimum_injector_pressure_pa",
                "minimum_manifold_pressure_pa",
                "minimum_propellant_mass_kg",
                "pressurant_gamma",
                "pressurant_gas_constant_j_kg_k",
                "propellant_density_kg_m3",
                "propellant_mass_kg",
                "regulator_discharge_coefficient",
                "regulator_full_open_error_pa",
                "regulator_maximum_area_m2",
                "regulator_opening",
                "regulator_set_pressure_pa",
                "relative_time_s",
                "sweep_points",
                "tank_gas_temperature_k",
                "tank_heat_transfer_w_k",
                "tank_internal_volume_m3",
                "tank_pressure_pa",
                "tank_wall_temperature_k",
                "valve_time_constant_s",
                "vehicle_mass_kg",
            }
        ),
        "/api/v1/transient/pulse-mode": frozenset(
            {
                "delta_v_m_s",
                "initial_mass_kg",
                "max_points",
                "maximum_on_time_s",
                "minimum_off_time_s",
                "propellant_budget_kg",
                "pulse_off_time_s",
                "steady_mass_flow_kg_s",
                "steady_thrust_n",
                "sweep_points",
            }
        ),
        "/api/v1/transient/regulated-feed": frozenset(
            {
                "ambient_pressure_pa",
                "bottle_gas_temperature_k",
                "bottle_heat_transfer_w_k",
                "bottle_pressure_pa",
                "bottle_volume_m3",
                "bottle_wall_temperature_k",
                "demand_sweep_maximum_flow_kg_s",
                "feed_line_resistance_pa_s2_kg2",
                "initial_mass_flow_kg_s",
                "maximum_tank_pressure_pa",
                "minimum_bottle_pressure_margin_pa",
                "minimum_injector_pressure_pa",
                "minimum_propellant_mass_kg",
                "pressurant_gamma",
                "pressurant_gas_constant_j_kg_k",
                "propellant_density_kg_m3",
                "propellant_mass_kg",
                "regulator_discharge_coefficient",
                "regulator_full_open_error_pa",
                "regulator_maximum_area_m2",
                "regulator_opening",
                "regulator_set_pressure_pa",
                "relative_time_s",
                "relaxation",
                "sweep_points",
                "tank_gas_temperature_k",
                "tank_heat_transfer_w_k",
                "tank_internal_volume_m3",
                "tank_pressure_pa",
                "tank_wall_temperature_k",
                "valve_time_constant_s",
                "vehicle_mass_kg",
            }
        ),
        "/api/v1/transient/scenarios": frozenset(
            {
                "burn_duration_s",
                "chamber_pressure_pa",
                "delivered_thrust_n",
                "max_points",
                "maximum_step_s",
                "system_specific_impulse_s",
                "throttle",
            }
        ),
        "/api/v1/transient/timeline": frozenset(
            {
                "end_time_s",
                "fuel_mass_kg",
                "initial_vehicle_mass_kg",
                "max_points",
                "oxidizer_mass_kg",
            }
        ),
        "/api/v1/transient/uncertainty": frozenset(
            {
                "cant_efficiency",
                "chamber_pressure_exponent",
                "chamber_pressure_pa",
                "convergence_tolerance",
                "delivered_thrust_n",
                "histogram_bins",
                "initial_pressure_pa",
                "initial_propellant_mass_kg",
                "initial_ullage_volume_m3",
                "integration_relative_tolerance",
                "maximum_mass_step_kg",
                "maximum_supply_pressure_pa",
                "minimum_supply_pressure_pa",
                "polytropic_exponent",
                "propellant_density_kg_m3",
                "reference_supply_pressure_pa",
                "reserve_mass_kg",
                "sample_count",
                "seed",
                "specific_impulse_pressure_exponent",
                "system_specific_impulse_s",
                "thrust_pressure_exponent",
            }
        ),
        "/api/v1/transpiration-cooling": frozenset(
            {
                "coolant_specific_heat_j_kg_k",
                "coolant_temperature_k",
                "radiative_heat_flux_w_m2",
                "recovery_temperature_k",
                "target_wall_temperature_k",
                "unblown_heat_transfer_coefficient_w_m2_k",
                "wall_efficiency",
            }
        ),
        "/api/v1/two-phase-band": frozenset(
            {
                "condensed_mass_fraction",
                "ideal_specific_impulse_s",
                "ideal_thrust_n",
                "maximum_particle_coupling_efficiency",
                "minimum_particle_coupling_efficiency",
                "reference_pressure_thrust_n",
            }
        ),
        "/api/v1/uncertainty": frozenset(),
        "/api/v1/wing": frozenset(
            {
                "angle_of_attack_deg",
                "aspect_ratio",
                "section_lift_slope_per_rad",
                "taper_ratio",
                "tip_twist_deg",
            }
        ),
        "/api/v1/workspace": frozenset(),
    }
)


EXACT_INTEGER_FIELDS = MappingProxyType(
    {
        "/api/v1/airfoil": frozenset(),
        "/api/v1/burn/export/sidera": frozenset(),
        "/api/v1/burn/report": frozenset(),
        "/api/v1/burn/simulate": frozenset(),
        "/api/v1/case-comparison": frozenset(),
        "/api/v1/chemistry-capabilities": frozenset(),
        "/api/v1/combustion-losses": frozenset(),
        "/api/v1/combustion-products": frozenset(),
        "/api/v1/combustion-stability": frozenset(),
        "/api/v1/duct-flow": frozenset(),
        "/api/v1/engine-cycle": frozenset(),
        "/api/v1/engine-operating-point": frozenset(),
        "/api/v1/entry": frozenset(),
        "/api/v1/equilibrium": frozenset(),
        "/api/v1/flat-plate": frozenset(),
        "/api/v1/heat-conduction": frozenset(),
        "/api/v1/heat-convection": frozenset(),
        "/api/v1/heat-radiation": frozenset(),
        "/api/v1/ideal-gas": frozenset(),
        "/api/v1/integrations/sidera/capabilities": frozenset(),
        "/api/v1/isentropic": frozenset(),
        "/api/v1/isentropic-sweep": frozenset(),
        "/api/v1/launch-vehicle-aero": frozenset(),
        "/api/v1/loss-budget": frozenset(),
        "/api/v1/normal-shock": frozenset(),
        "/api/v1/nozzle": frozenset(),
        "/api/v1/nozzle-chemistry": frozenset(),
        "/api/v1/nozzle-contour": frozenset(),
        "/api/v1/nozzle-off-design": frozenset(),
        "/api/v1/nozzle-thermal": frozenset(),
        "/api/v1/oblique-shock": frozenset(),
        "/api/v1/of-sweep": frozenset(),
        "/api/v1/polytropic": frozenset(),
        "/api/v1/propellant-vaporization": frozenset(),
        "/api/v1/propellants": frozenset(),
        "/api/v1/real-fluid": frozenset(),
        "/api/v1/regen-chamber": frozenset(),
        "/api/v1/rocket-equation": frozenset(),
        "/api/v1/rocket-performance": frozenset(),
        "/api/v1/sass/mission": frozenset(),
        "/api/v1/species": frozenset(),
        "/api/v1/staging": frozenset(),
        "/api/v1/stagnation": frozenset(),
        "/api/v1/supersonic-body": frozenset(),
        "/api/v1/thermal-properties": frozenset(),
        "/api/v1/thermo-process": frozenset(),
        "/api/v1/thrust": frozenset(),
        "/api/v1/transient/blowdown": frozenset(),
        "/api/v1/transient/calibration": frozenset(),
        "/api/v1/transient/electric": frozenset(),
        "/api/v1/transient/feed-line": frozenset(),
        "/api/v1/transient/pulse-mode": frozenset(),
        "/api/v1/transient/regulated-feed": frozenset(),
        "/api/v1/transient/scenarios": frozenset(),
        "/api/v1/transient/timeline": frozenset(),
        "/api/v1/transient/uncertainty": frozenset(),
        "/api/v1/transpiration-cooling": frozenset(),
        "/api/v1/two-phase-band": frozenset(),
        "/api/v1/uncertainty": frozenset({"sample_count", "seed", "sobol_base_samples"}),
        "/api/v1/wing": frozenset(),
        "/api/v1/workspace": frozenset(),
    }
)


NUMERIC_FIELDS = MappingProxyType(
    {path: NUMBER_FIELDS[path] | EXACT_INTEGER_FIELDS[path] for path in NUMBER_FIELDS}
)
