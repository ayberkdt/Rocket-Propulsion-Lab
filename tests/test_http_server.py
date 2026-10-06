"""End-to-end tests for the versioned HTTP response contract."""

from __future__ import annotations

import json
import threading
import unittest
from http.server import ThreadingHTTPServer
from urllib.error import HTTPError
from urllib.request import Request, urlopen

from rocket_propulsion.api.server import RocketPropulsionHandler


class QuietHandler(RocketPropulsionHandler):
    """Avoid request-log noise in the test output."""

    def log_message(self, message_format: str, *args: object) -> None:
        pass


class HttpContractTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.server = ThreadingHTTPServer(("127.0.0.1", 0), QuietHandler)
        cls.thread = threading.Thread(target=cls.server.serve_forever, daemon=True)
        cls.thread.start()
        cls.base_url = f"http://127.0.0.1:{cls.server.server_port}"

    @classmethod
    def tearDownClass(cls) -> None:
        cls.server.shutdown()
        cls.server.server_close()
        cls.thread.join(timeout=2.0)

    def post(self, path: str, payload: object) -> tuple[int, dict[str, object]]:
        request = Request(
            self.base_url + path,
            data=json.dumps(payload).encode("utf-8"),
            headers={"Content-Type": "application/json"},
            method="POST",
        )
        try:
            with urlopen(request, timeout=2.0) as response:
                return response.status, json.load(response)
        except HTTPError as error:
            return error.code, json.load(error)

    def test_success_has_data_and_model_metadata(self) -> None:
        status, body = self.post(
            "/api/v1/isentropic", {"input_kind": "mach", "value": 2.0}
        )
        self.assertEqual(status, 200)
        self.assertAlmostEqual(body["data"]["mach"], 2.0)
        self.assertEqual(body["meta"]["model"], "perfect-gas isentropic relations")
        self.assertIn("assumptions", body["meta"])
        self.assertIn("units", body["meta"])

    def test_invalid_field_has_stable_error_contract(self) -> None:
        status, body = self.post(
            "/api/v1/normal-shock", {"upstream_mach": "not-a-number"}
        )
        self.assertEqual(status, 400)
        self.assertEqual(body["error"]["code"], "invalid_input")
        self.assertEqual(body["error"]["field"], "upstream_mach")
        self.assertIn("message", body["error"])

    def test_unknown_route_is_structured_not_found(self) -> None:
        status, body = self.post("/api/v1/not-real", {})
        self.assertEqual(status, 404)
        self.assertEqual(body["error"]["code"], "unknown_route")

    def test_sidera_capabilities_are_available_without_sidera(self) -> None:
        with urlopen(
            self.base_url + "/api/v1/integrations/sidera/capabilities", timeout=2.0
        ) as response:
            body = json.load(response)
        self.assertEqual(response.status, 200)
        self.assertIn(body["data"]["status"], {"unavailable", "incompatible", "supported"})
        self.assertEqual(body["meta"]["model"], "Sidera public finite-burn capability probe")

    def test_browser_contains_standalone_burn_workbench(self) -> None:
        with urlopen(self.base_url + "/", timeout=2.0) as response:
            html = response.read().decode("utf-8")
        self.assertIn('id="burn-form"', html)
        self.assertIn('id="burn-sidera"', html)
        self.assertIn('id="burn-mass-chart"', html)
        self.assertIn('name="profile_mode"', html)
        self.assertIn('id="burn-thrust-chart"', html)
        self.assertIn('id="sidera-mode"', html)

    def test_burn_provider_and_simulation_over_http(self) -> None:
        status, point_body = self.post(
            "/api/v1/engine-operating-point",
            {
                "family": "monopropellant",
                "delivered_thrust_n": 980.665,
                "system_specific_impulse_s": 100.0,
                "chamber_pressure_pa": 1_000_000.0,
                "propellant_key": "hydrazine",
            },
        )
        self.assertEqual(status, 200)
        status, burn_body = self.post(
            "/api/v1/burn/simulate",
            {
                "definition": {
                    "name": "HTTP burn",
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
                "operating_point": point_body["data"],
            },
        )
        self.assertEqual(status, 200)
        self.assertAlmostEqual(burn_body["data"]["summary"]["final_mass_kg"], 999.0)

        status, profile_body = self.post(
            "/api/v1/burn/simulate",
            {
                "definition": {
                    "name": "HTTP profile",
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
                "operating_point": point_body["data"],
                "schedule": {
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
                    ]
                },
            },
        )
        self.assertEqual(status, 200)
        self.assertEqual(profile_body["data"]["schema"], "rocket_propulsion_burn_v2")
        self.assertAlmostEqual(
            profile_body["data"]["summary"]["consumed_propellant_kg"], 2.0
        )


if __name__ == "__main__":
    unittest.main()
