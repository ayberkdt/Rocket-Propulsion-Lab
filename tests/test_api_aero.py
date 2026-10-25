"""HTTP routes for real fluids and aerodynamics: defaults, strict JSON and references."""

from __future__ import annotations

import json
import unittest
from math import radians
from pathlib import Path

from rocket_propulsion.aerodynamics import airfoil_geometry, thin_airfoil
from rocket_propulsion.api.routes import ROUTE_METADATA, dispatch_calculation
from rocket_propulsion.compressible.cone_flow import solve_cone_flow
from rocket_propulsion.core.errors import DomainError, InputError

NIST = json.loads((Path(__file__).parent / "reference_data" / "nist_real_fluids.json").read_text(encoding="utf-8"))
ROUTES = ("real-fluid", "airfoil", "wing", "supersonic-body", "flat-plate", "entry", "launch-vehicle-aero")


class AeroRouteDefaultsTests(unittest.TestCase):
    def test_every_route_answers_an_empty_request_with_strict_json(self) -> None:
        for name in ROUTES:
            path = f"/api/v1/{name}"
            with self.subTest(route=name):
                response = dispatch_calculation(path, {})
                json.dumps(response, allow_nan=False)
                self.assertIn(path, ROUTE_METADATA)


class AeroRouteReferenceTests(unittest.TestCase):
    def test_real_fluid_state_reproduces_nist_rows(self) -> None:
        for row in NIST["fluids"]["oxygen"]["single_phase"][:4]:
            with self.subTest(temperature=row["T_K"], pressure=row["p_MPa"]):
                state = dispatch_calculation("/api/v1/real-fluid", {
                    "fluid": "oxygen", "temperature_k": row["T_K"], "pressure_pa": row["p_MPa"] * 1.0e6,
                    "point_count": 11})["state"]
                self.assertAlmostEqual(state["density_kg_m3"] / row["rho_kg_m3"], 1.0, delta=1.0e-5)
                self.assertAlmostEqual(state["speed_of_sound_m_s"] / row["w_m_s"], 1.0, delta=1.0e-5)

    def test_airfoil_route_matches_the_library_thin_airfoil_result(self) -> None:
        response = dispatch_calculation("/api/v1/airfoil", {"designation": "2412", "angle_of_attack_deg": 4})
        expected = thin_airfoil.thin_airfoil_naca(airfoil_geometry.naca_section("2412"), radians(4.0))
        self.assertAlmostEqual(response["thin"]["lift_coefficient"], expected.lift_coefficient, delta=1.0e-12)
        self.assertEqual([row["angle_of_attack_deg"] for row in response["sweep"]], [float(a) for a in range(-4, 13, 2)])

    def test_cone_route_matches_the_taylor_maccoll_solver(self) -> None:
        response = dispatch_calculation("/api/v1/supersonic-body", {"mach": 3, "cone_half_angle_deg": 15})
        expected = solve_cone_flow(3.0, 15.0)
        self.assertAlmostEqual(response["cone"]["shock_angle_deg"], expected.shock_angle_deg, delta=1.0e-12)

    def test_entry_route_terminates_and_reports_heating(self) -> None:
        response = dispatch_calculation("/api/v1/entry", {})
        self.assertIn(response["termination"], {"velocity", "altitude"})
        self.assertGreater(response["peak_heating"]["stagnation_heating_w_m2"], 0.0)
        self.assertLessEqual(len(response["trajectory"]), 302)

    def test_launch_vehicle_slender_body_slope_is_two_per_radian(self) -> None:
        response = dispatch_calculation("/api/v1/launch-vehicle-aero", {})
        # Slender-body theory: CN_alpha = 2 on the base area for a nose-cylinder body.
        self.assertAlmostEqual(response["normal_force"]["normal_force_slope_per_rad"], 2.0, delta=1.0e-12)


class AeroRouteInputTests(unittest.TestCase):
    def test_invalid_inputs_are_refused(self) -> None:
        with self.assertRaises(InputError):
            dispatch_calculation("/api/v1/real-fluid", {"fluid": "unobtainium"})
        with self.assertRaises(InputError):
            dispatch_calculation("/api/v1/airfoil", {"panels_per_surface": 5})
        with self.assertRaises(DomainError):
            dispatch_calculation("/api/v1/wing", {"taper_ratio": 0})
        with self.assertRaises(InputError):
            dispatch_calculation("/api/v1/launch-vehicle-aero", {"nose_shape": "blunt"})
        with self.assertRaises(InputError):
            dispatch_calculation("/api/v1/flat-plate", {"mach": True})


if __name__ == "__main__":
    unittest.main()
