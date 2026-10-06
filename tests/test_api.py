"""Tests for JSON-shaped route adapters."""

import json
import unittest

from rocket_propulsion.api.routes import ROUTES, calculation_metadata, dispatch_calculation
from rocket_propulsion.core.errors import DomainError


class ApiRouteTests(unittest.TestCase):
    def test_every_route_has_model_metadata(self) -> None:
        for path in ROUTES:
            with self.subTest(path=path):
                metadata = calculation_metadata(path)
                self.assertTrue(metadata["model"])
                self.assertTrue(metadata["assumptions"])
                self.assertIsInstance(metadata["units"], dict)

    def test_isentropic_route_supports_string_form_values(self) -> None:
        result = dispatch_calculation(
            "/api/v1/isentropic",
            {"input_kind": "temperature_total_ratio", "value": "0.555555555555", "gamma": "1.4"},
        )
        self.assertAlmostEqual(result["mach"], 2.0, places=8)

    def test_normal_shock_route(self) -> None:
        result = dispatch_calculation(
            "/api/v1/normal-shock", {"upstream_mach": 2.0, "gamma": 1.4}
        )
        self.assertAlmostEqual(result["pressure_ratio"], 4.5)

    def test_ideal_gas_route_treats_empty_field_as_missing(self) -> None:
        result = dispatch_calculation(
            "/api/v1/ideal-gas",
            {
                "pressure_pa": "101325",
                "temperature_k": "288.15",
                "density_kg_m3": "",
                "gas_constant_j_kg_k": "287.05",
                "gamma": "1.4",
            },
        )
        self.assertGreater(result["density_kg_m3"], 1.2)

    def test_invalid_number_is_a_domain_error(self) -> None:
        with self.assertRaises(DomainError):
            dispatch_calculation(
                "/api/v1/normal-shock", {"upstream_mach": "not-a-number"}
            )

    def test_unknown_route_is_key_error(self) -> None:
        with self.assertRaises(KeyError):
            dispatch_calculation("/api/v1/unknown", {})

    def test_rocket_equation_route(self) -> None:
        result = dispatch_calculation(
            "/api/v1/rocket-equation",
            {
                "initial_mass_kg": "10000",
                "final_mass_kg": "4000",
                "delta_v_m_s": "",
                "specific_impulse_s": "320",
                "effective_exhaust_velocity_m_s": "",
            },
        )
        self.assertGreater(result["delta_v_m_s"], 2800.0)

    def test_duct_flow_and_propellant_routes(self) -> None:
        fanno = dispatch_calculation(
            "/api/v1/duct-flow", {"model": "fanno", "mach": 0.5, "gamma": 1.4}
        )
        catalog = dispatch_calculation("/api/v1/propellants", {"category": "hybrid"})
        self.assertGreater(fanno["friction_parameter_to_sonic"], 0.0)
        self.assertEqual(catalog["count"], 2)

    def test_thermodynamic_process_route(self) -> None:
        result = dispatch_calculation(
            "/api/v1/thermo-process",
            {
                "inlet_pressure_pa": 100_000.0,
                "outlet_pressure_pa": 400_000.0,
                "inlet_temperature_k": 300.0,
                "efficiency": 0.8,
                "gamma": 1.4,
                "gas_constant_j_kg_k": 287.05,
            },
        )
        self.assertEqual(result["process"], "compression")
        self.assertGreater(result["entropy_change_j_kg_k"], 0.0)

    def test_oblique_shock_geometry_and_curve_routes(self) -> None:
        shock = dispatch_calculation(
            "/api/v1/oblique-shock",
            {"upstream_mach": 2.0, "turning_angle_deg": 10.0, "branch": "weak"},
        )
        contour = dispatch_calculation(
            "/api/v1/nozzle-contour",
            {"throat_area_m2": 0.01, "area_ratio": 20.0, "contour": "bell"},
        )
        sweep = dispatch_calculation(
            "/api/v1/isentropic-sweep",
            {"minimum_mach": 0.1, "maximum_mach": 4.0, "point_count": 41},
        )
        self.assertGreater(shock["wave_angle_deg"], shock["turning_angle_deg"])
        self.assertEqual(len(contour["stations"]), 121)
        self.assertEqual(len(sweep["points"]), 41)

    def test_rao_and_moc_geometry_routes_include_evidence_and_exports(self) -> None:
        rao = dispatch_calculation(
            "/api/v1/nozzle-contour",
            {
                "throat_area_m2": 0.01,
                "area_ratio": 16.0,
                "contour": "rao",
                "include_exports": "true",
                "angular_segments": 8,
            },
        )
        moc = dispatch_calculation(
            "/api/v1/nozzle-contour",
            {
                "throat_area_m2": 0.01,
                "area_ratio": 16.0,
                "contour": "moc",
                "characteristic_count": 10,
            },
        )
        self.assertIn("solid rocket_propulsion_nozzle", rao["exports"]["stl"])
        self.assertTrue(moc["characteristic_lines"])
        self.assertLess(moc["residuals"]["exit_mach_absolute"], 1.0e-8)

    def test_loss_and_thermal_routes_close_balances(self) -> None:
        losses = dispatch_calculation(
            "/api/v1/loss-budget",
            {"ideal_thrust_n": 100_000.0, "ideal_specific_impulse_s": 450.0},
        )
        thermal = dispatch_calculation(
            "/api/v1/nozzle-thermal",
            {
                "contour": "rao",
                "throat_area_m2": 0.01,
                "area_ratio": 25.0,
                "gamma": 1.22,
                "chamber_pressure_pa": 5.0e6,
                "chamber_temperature_k": 3500.0,
                "wall_temperature_k": 800.0,
                "characteristic_velocity_m_s": 1700.0,
                "dynamic_viscosity_pa_s": 7.0e-5,
                "specific_heat_j_kg_k": 3000.0,
                "prandtl_number": 0.7,
                "throat_radius_of_curvature_m": 0.03,
                "include_cooling": True,
                "coolant_mass_flow_kg_s": 5.0,
                "coolant_specific_heat_j_kg_k": 10_000.0,
                "coolant_inlet_temperature_k": 100.0,
                "coolant_density_kg_m3": 70.0,
                "coolant_dynamic_viscosity_pa_s": 1.0e-5,
                "channel_count": 100,
                "channel_flow_area_m2": 2.0e-5,
                "hydraulic_diameter_m": 0.004,
                "channel_length_m": 1.0,
                "roughness_m": 1.0e-5,
                "maximum_coolant_temperature_k": 1000.0,
            },
        )
        self.assertAlmostEqual(losses["thrust_balance_residual_n"], 0.0)
        self.assertAlmostEqual(thermal["cooling"]["energy_balance_residual_w"], 0.0)
        self.assertGreater(thermal["peak_heat_flux_w_m2"], 0.0)

    def test_staging_cycle_uncertainty_and_workspace_routes(self) -> None:
        stages = dispatch_calculation(
            "/api/v1/staging",
            {
                "stages": [
                    {"name": "Booster", "specific_impulse_s": 300, "structural_fraction": 0.1},
                    {"name": "Upper", "specific_impulse_s": 450, "structural_fraction": 0.08},
                ],
                "target_delta_v_m_s": 9400,
                "payload_mass_kg": 1000,
                "resolution": 100,
            },
        )
        cycle = dispatch_calculation(
            "/api/v1/engine-cycle",
            {
                "cycle": "gas-generator",
                "chamber_pressure_pa": 10e6,
                "chamber_mass_flow_kg_s": 100,
                "propellant_density_kg_m3": 1000,
                "inlet_pressure_pa": 0.3e6,
                "injector_pressure_drop_fraction": 0.2,
                "pump_efficiency": 0.7,
                "turbine_efficiency": 0.6,
                "turbine_inlet_temperature_k": 900,
                "turbine_specific_heat_j_kg_k": 2500,
                "turbine_gamma": 1.25,
                "turbine_pressure_ratio": 8,
            },
        )
        uncertainty = dispatch_calculation(
            "/api/v1/uncertainty",
            {
                "parameters": [
                    {"name": "specific_impulse_s", "minimum": 430, "maximum": 460},
                    {"name": "mass_ratio", "minimum": 3.0, "maximum": 4.0},
                ],
                "sample_count": 40,
                "seed": 7,
            },
        )
        workspace_text = json.dumps(
            {
                "schema_version": 2,
                "name": "API demo",
                "cases": [
                    {
                        "case_id": "a",
                        "name": "A",
                        "model": "rocket",
                        "inputs": {},
                        "results": {"delta_v_m_s": 4000},
                    }
                ],
            }
        )
        workspace = dispatch_calculation(
            "/api/v1/workspace",
            {"workspace_json": workspace_text, "include_report": True},
        )
        self.assertAlmostEqual(stages["achieved_delta_v_m_s"], 9400.0)
        self.assertAlmostEqual(cycle["power_balance_residual_w"], 0.0)
        self.assertEqual(uncertainty["seed"], 7)
        self.assertIn("<!doctype html>", workspace["report"]["html"])

    def test_polytropic_route_includes_diagram_path(self) -> None:
        result = dispatch_calculation(
            "/api/v1/polytropic",
            {
                "inlet_pressure_pa": 100_000.0,
                "outlet_pressure_pa": 400_000.0,
                "inlet_temperature_k": 300.0,
                "exponent": 1.3,
            },
        )
        self.assertEqual(result["process"], "compression")
        self.assertEqual(len(result["path"]), 61)

    def test_species_catalog_and_variable_property_sweep(self) -> None:
        catalog = dispatch_calculation("/api/v1/species", {})
        result = dispatch_calculation(
            "/api/v1/thermal-properties",
            {
                "species": "H2O",
                "temperature_k": 3500.0,
                "minimum_temperature_k": 200.0,
                "maximum_temperature_k": 6000.0,
                "point_count": 41,
            },
        )
        self.assertEqual(catalog["count"], 10)
        self.assertEqual(result["model"]["kind"], "species")
        self.assertEqual(len(result["curve"]), 41)
        self.assertGreater(result["properties"]["cp_j_kg_k"], 2_000.0)

    def test_variable_property_mixture_route(self) -> None:
        result = dispatch_calculation(
            "/api/v1/thermal-properties",
            {
                "composition": {"N2": 0.79, "O2": 0.21},
                "basis": "mole",
                "temperature_k": 300.0,
                "point_count": 21,
            },
        )
        self.assertEqual(result["model"]["kind"], "mixture")
        self.assertAlmostEqual(result["properties"]["gamma"], 1.398476, places=5)

    def test_local_equilibrium_and_capability_routes(self) -> None:
        capabilities = dispatch_calculation("/api/v1/chemistry-capabilities", {})
        result = dispatch_calculation(
            "/api/v1/equilibrium",
            {
                "provider": "local",
                "problem_type": "TP",
                "reactants": {"H2": 2.0, "O2": 1.0},
                "temperature_k": 3_500.0,
                "pressure_pa": 5_000_000.0,
            },
        )
        self.assertEqual(capabilities["providers"][0]["key"], "local")
        self.assertEqual(result["provider"], "local-gibbs")
        self.assertLess(result["residuals"]["maximum_element_relative"], 1.0e-9)

    def test_oxidizer_fuel_sweep_route(self) -> None:
        result = dispatch_calculation(
            "/api/v1/of-sweep",
            {
                "propellant_pair": "LOX_LH2",
                "chamber_pressure_pa": 5_000_000.0,
                "area_ratio": 10.0,
                "minimum_ratio": 5.0,
                "maximum_ratio": 6.0,
                "point_count": 3,
            },
        )
        self.assertEqual(len(result["points"]), 3)
        self.assertGreater(result["points"][1]["adiabatic_flame_temperature_k"], 3_000.0)
        self.assertGreater(
            result["points"][1]["equilibrium_vacuum_specific_impulse_s"],
            result["points"][1]["frozen_vacuum_specific_impulse_s"],
        )

    def test_off_design_nozzle_route_includes_altitude_sweep(self) -> None:
        result = dispatch_calculation(
            "/api/v1/nozzle-off-design",
            {
                "chamber_pressure_pa": 7_000_000.0,
                "chamber_temperature_k": 3_500.0,
                "throat_area_m2": 0.01,
                "area_ratio": 20.0,
                "ambient_pressure_pa": 101_325.0,
                "gamma": 1.22,
                "gas_constant_j_kg_k": 355.0,
                "altitude_point_count": 7,
            },
        )
        self.assertEqual(result["regime"], "overexpanded")
        self.assertTrue(result["separation"]["predicted"])
        self.assertEqual(len(result["altitude_sweep"]["points"]), 7)

    def test_standalone_constant_burn_route(self) -> None:
        result = dispatch_calculation(
            "/api/v1/burn/simulate",
            {
                "definition": {
                    "name": "API burn",
                    "initial_mass_kg": 1000.0,
                    "protected_dry_mass_kg": 800.0,
                    "maximum_duration_s": 10.0,
                    "target": {"kind": "duration", "value": 1.0},
                    "tanks": [
                        {
                            "tank_id": "propellant",
                            "propellant_key": "hydrazine",
                            "role": "monopropellant",
                            "loaded_mass_kg": 100.0,
                            "reserve_mass_kg": 10.0,
                        }
                    ],
                },
                "operating_point": {
                    "name": "API point",
                    "ideal_thrust_n": 1000.0,
                    "delivered_thrust_n": 980.665,
                    "chamber_specific_impulse_s": 100.0,
                    "system_specific_impulse_s": 100.0,
                    "chamber_mass_flow_kg_s": 1.0,
                    "total_tank_flow_kg_s": 1.0,
                    "chamber_pressure_pa": 1_000_000.0,
                    "streams": [
                        {
                            "stream_id": "main",
                            "tank_id": "propellant",
                            "propellant_key": "hydrazine",
                            "role": "monopropellant",
                            "destination": "chamber",
                            "mass_flow_kg_s": 1.0,
                        }
                    ],
                },
            },
        )
        self.assertEqual(result["schema"], "rocket_propulsion_burn_v1")
        self.assertTrue(result["summary"]["target_achieved"])
        self.assertAlmostEqual(result["summary"]["final_mass_kg"], 999.0)
        self.assertEqual(len(result["input_hash"]), 64)
        self.assertEqual(len(result["result_hash"]), 64)

        export = dispatch_calculation(
            "/api/v1/burn/export/sidera",
            {
                "burn_result": result,
                "t_start_s": 60.0,
                "direction": [1.0, 0.0, 0.0],
                "frame": "inertial",
            },
        )
        capabilities = dispatch_calculation(
            "/api/v1/integrations/sidera/capabilities", {}
        )
        self.assertEqual(export["schema"], "sidera_maneuver_plan_v1")
        self.assertEqual(export["maneuvers"][0]["throttle"], 1.0)
        self.assertIn(capabilities["status"], {"unavailable", "incompatible", "supported"})

        report = dispatch_calculation(
            "/api/v1/burn/report", {"burn_result": result}
        )
        self.assertIn("<!doctype html>", report["html"])
        self.assertEqual(report["result_hash"], result["result_hash"])

    def test_engine_operating_point_provider_route(self) -> None:
        result = dispatch_calculation(
            "/api/v1/engine-operating-point",
            {
                "family": "bipropellant",
                "delivered_thrust_n": 190_000.0,
                "system_specific_impulse_s": 190_000.0 / (9.80665 * 51.0),
                "ideal_thrust_n": 200_000.0,
                "chamber_specific_impulse_s": 190_000.0 / (9.80665 * 50.0),
                "oxidizer_fuel_ratio": 40.0 / 11.0,
                "dump_mass_flow_kg_s": 1.0,
                "dump_role": "fuel",
                "chamber_pressure_pa": 7_000_000.0,
                "oxidizer_key": "lox",
                "fuel_key": "rp1",
            },
        )
        self.assertAlmostEqual(result["total_tank_flow_kg_s"], 51.0)
        self.assertAlmostEqual(result["chamber_mass_flow_kg_s"], 50.0)
        self.assertEqual(result["streams"][-1]["destination"], "gas_generator_dump")

    def test_profiled_burn_route_uses_piecewise_linear_schedule(self) -> None:
        result = dispatch_calculation(
            "/api/v1/burn/simulate",
            {
                "definition": {
                    "name": "API profiled burn",
                    "initial_mass_kg": 1000.0,
                    "protected_dry_mass_kg": 800.0,
                    "maximum_duration_s": 10.0,
                    "target": {"kind": "duration", "value": 3.0},
                    "tanks": [
                        {
                            "tank_id": "propellant",
                            "propellant_key": "hydrazine",
                            "role": "monopropellant",
                            "loaded_mass_kg": 100.0,
                            "reserve_mass_kg": 10.0,
                        }
                    ],
                },
                "operating_point": {
                    "name": "rated point",
                    "ideal_thrust_n": 1000.0,
                    "delivered_thrust_n": 980.665,
                    "chamber_specific_impulse_s": 100.0,
                    "system_specific_impulse_s": 100.0,
                    "chamber_mass_flow_kg_s": 1.0,
                    "total_tank_flow_kg_s": 1.0,
                    "chamber_pressure_pa": 1_000_000.0,
                    "streams": [
                        {
                            "stream_id": "main",
                            "tank_id": "propellant",
                            "propellant_key": "hydrazine",
                            "role": "monopropellant",
                            "destination": "chamber",
                            "mass_flow_kg_s": 1.0,
                        }
                    ],
                },
                "schedule": {
                    "name": "API ramp",
                    "segments": [
                        {
                            "duration_s": 2.0,
                            "start_throttle": 0.0,
                            "end_throttle": 1.0,
                            "phase": "ramp-up",
                        },
                        {
                            "duration_s": 1.0,
                            "start_throttle": 1.0,
                            "end_throttle": 1.0,
                            "phase": "steady",
                        },
                    ],
                },
            },
        )
        self.assertEqual(result["schema"], "rocket_propulsion_burn_v2")
        self.assertEqual(result["evidence"]["segment_count"], 2)
        self.assertAlmostEqual(result["summary"]["consumed_propellant_kg"], 2.0)
        self.assertFalse(result["summary"]["exact_constant_sidera_export_available"])
        report = dispatch_calculation(
            "/api/v1/burn/report", {"burn_result": result}
        )
        self.assertIn("Throttle schedule", report["html"])
        approximate = dispatch_calculation(
            "/api/v1/burn/export/sidera",
            {
                "burn_result": result,
                "t_start_s": 0.0,
                "direction": [1.0, 0.0, 0.0],
                "mode": "equivalent_constant",
            },
        )
        self.assertEqual(approximate["manifest"]["mapping"], "equivalent_constant")


if __name__ == "__main__":
    unittest.main()

