"""Contract tests for the transient-core and SASS HTTP routes."""

from __future__ import annotations

import json
import math
import threading
import unittest
from http.server import ThreadingHTTPServer
from pathlib import Path
from unittest.mock import patch
from urllib.error import HTTPError
from urllib.request import Request, urlopen

from rocket_propulsion.api import transient_routes
from rocket_propulsion.api.routes import (
    ROUTE_METADATA,
    ROUTES,
    calculation_metadata,
    dispatch_calculation,
)
from rocket_propulsion.api.server import RocketPropulsionHandler
from rocket_propulsion.core.errors import ConvergenceError, DomainError, InputError

ROOT = Path(__file__).resolve().parents[1]
PATHS = tuple(transient_routes.TRANSIENT_ROUTES)
# Inexpensive mesh controls for the blowdown-based routes in non-default cases.
FAST_BLOWDOWN = {"relative_tolerance": 1e-4}
FAST_ENSEMBLE = {"sample_count": 8, "integration_relative_tolerance": 1e-3}

_DEFAULT_RESULTS: dict[str, dict] = {}


def default_result(path: str) -> dict:
    """Run each default request once per test process (blowdown routes take ~1 s)."""

    if path not in _DEFAULT_RESULTS:
        _DEFAULT_RESULTS[path] = dispatch_calculation(path, {})
    return _DEFAULT_RESULTS[path]


def assert_series_contract(test: unittest.TestCase, group: dict) -> None:
    lengths = {len(values) for values in group["columns"].values()}
    test.assertEqual(lengths, {group["point_count"]})
    test.assertIn(group["x"], group["columns"])
    test.assertLessEqual(group["point_count"], transient_routes.MAX_SERIES_POINTS)
    test.assertLessEqual(group["point_count"], group["source_point_count"])


class RegistryTests(unittest.TestCase):
    def test_routes_are_registered_with_metadata(self) -> None:
        self.assertEqual(len(PATHS), 10)
        for path in PATHS:
            self.assertIs(ROUTES[path], transient_routes.TRANSIENT_ROUTES[path])
            self.assertIn(path, ROUTE_METADATA)
            meta = calculation_metadata(path)
            self.assertTrue(meta["model"])
            self.assertTrue(meta["assumptions"])
            self.assertTrue(meta["references"])

    def test_decimation_is_deterministic_and_keeps_end_points(self) -> None:
        indices = transient_routes._decimated_indices(1001, 400)
        self.assertEqual(len(indices), 400)
        self.assertEqual((indices[0], indices[-1]), (0, 1000))
        self.assertTrue(all(a < b for a, b in zip(indices, indices[1:], strict=False)))
        self.assertEqual(indices, transient_routes._decimated_indices(1001, 400))
        self.assertEqual(transient_routes._decimated_indices(5, 400), [0, 1, 2, 3, 4])

    def test_default_sass_study_is_the_example_file(self) -> None:
        example = json.loads((ROOT / "examples" / "sass_mission.json").read_text(encoding="utf-8"))
        self.assertEqual(transient_routes.default_sass_study(), example)


class DefaultRequestTests(unittest.TestCase):
    def test_every_route_answers_an_empty_request(self) -> None:
        for path in PATHS:
            with self.subTest(path=path):
                result = default_result(path)
                json.dumps(result, allow_nan=False)
                for key in ("model", "inputs", "summary", "series", "reference_ids", "warnings"):
                    self.assertIn(key, result)
                self.assertTrue(result["series"])
                for group in result["series"].values():
                    assert_series_contract(self, group)
                self.assertTrue(all(isinstance(item, str) for item in result["reference_ids"]))

    def test_unknown_fields_are_refused_by_every_route(self) -> None:
        for path in PATHS:
            with self.subTest(path=path), self.assertRaises(InputError) as caught:
                dispatch_calculation(path, {"not_a_field": 1})
            self.assertEqual(caught.exception.code, "unknown_field")


class BlowdownRouteTests(unittest.TestCase):
    def test_default_reproduces_the_example_case(self) -> None:
        result = default_result("/api/v1/transient/blowdown")
        summary = result["summary"]
        self.assertEqual(summary["cutoff_reason"], "minimum_supply_pressure")
        self.assertAlmostEqual(summary["initial_supply_pressure_pa"], 2.0e6)
        self.assertAlmostEqual(summary["final_supply_pressure_pa"], 1.0e6, delta=1.0)
        history = result["series"]["history"]["columns"]
        self.assertEqual(history["time_s"][0], 0.0)
        self.assertAlmostEqual(history["time_s"][-1], summary["duration_s"])
        self.assertEqual(result["inputs"]["maximum_mass_step_kg"], 0.5)

    def test_reserve_cutoff_and_decimation(self) -> None:
        payload = {**FAST_BLOWDOWN, "reserve_mass_kg": 5.0, "max_points": 10}
        result = dispatch_calculation("/api/v1/transient/blowdown", payload)
        self.assertEqual(result["summary"]["cutoff_reason"], "propellant_reserve")
        self.assertAlmostEqual(result["evidence"]["consumed_propellant_kg"], 5.0, places=12)
        # The tabulated artifact closes mass within the requested relative tolerance.
        self.assertAlmostEqual(result["summary"]["consumed_propellant_kg"], 5.0, delta=5.0e-4)
        group = result["series"]["history"]
        assert_series_contract(self, group)
        self.assertEqual(group["point_count"], 10)
        self.assertGreater(group["source_point_count"], 10)
        self.assertEqual(result, dispatch_calculation("/api/v1/transient/blowdown", payload))

    def test_invalid_inputs(self) -> None:
        with self.assertRaises(InputError):
            dispatch_calculation("/api/v1/transient/blowdown", {"initial_pressure_pa": "high"})
        with self.assertRaises(InputError) as caught:
            dispatch_calculation("/api/v1/transient/blowdown", {"maximum_mass_step_kg": 1e-6})
        self.assertEqual(caught.exception.code, "work_limit")
        with self.assertRaises(InputError):
            dispatch_calculation("/api/v1/transient/blowdown", {"max_points": 1})
        with self.assertRaises(DomainError):
            dispatch_calculation(
                "/api/v1/transient/blowdown", {**FAST_BLOWDOWN, "initial_pressure_pa": 3.0e6}
            )


class FeedRouteTests(unittest.TestCase):
    def test_regulated_feed_default_and_sweeps(self) -> None:
        result = default_result("/api/v1/transient/regulated-feed")
        feed = result["summary"]["feed"]
        self.assertTrue(result["evidence"]["converged"])
        self.assertLess(feed["injector_pressure_pa"], feed["tank_pressure_pa"])
        command = result["series"]["regulator_sweep"]["columns"]["regulator_command"]
        self.assertGreater(command[0], command[-1])
        demand = result["series"]["demand_sweep"]["columns"]
        self.assertEqual(demand["propellant_mass_flow_kg_s"][0], 0.0)
        self.assertAlmostEqual(demand["propellant_mass_flow_kg_s"][-1], 0.5)

    def test_regulated_feed_resistance_increases_line_drop(self) -> None:
        base = default_result("/api/v1/transient/regulated-feed")["summary"]["feed"]
        result = dispatch_calculation(
            "/api/v1/transient/regulated-feed",
            {"feed_line_resistance_pa_s2_kg2": 1.0e6, "sweep_points": 7},
        )
        self.assertGreater(
            result["summary"]["feed"]["feed_line_pressure_drop_pa"],
            base["feed_line_pressure_drop_pa"],
        )
        self.assertEqual(result["series"]["demand_sweep"]["point_count"], 7)

    def test_regulated_feed_rejects_bad_state(self) -> None:
        with self.assertRaises(DomainError):
            dispatch_calculation("/api/v1/transient/regulated-feed", {"propellant_mass_kg": -1.0})
        with self.assertRaises(InputError):
            dispatch_calculation("/api/v1/transient/regulated-feed", {"sweep_points": 401})

    def test_feed_line_default_and_overfilled_line(self) -> None:
        line = default_result("/api/v1/transient/feed-line")["summary"]["line"]
        self.assertGreater(line["natural_angular_frequency_rad_s"], 0.0)
        result = dispatch_calculation(
            "/api/v1/transient/feed-line", {"line_mass_flow_kg_s": 0.3}
        )
        derivative = result["summary"]["line"]["derivative"]
        # Line flow above engine demand fills the compliant manifold.
        self.assertGreater(derivative["manifold_pressure_pa_s"], 0.0)
        with self.assertRaises(DomainError):
            dispatch_calculation("/api/v1/transient/feed-line", {"line_compliance_kg_pa": 0.0})


class ElectricAndPulseRouteTests(unittest.TestCase):
    def test_electric_default_selects_the_example_level(self) -> None:
        result = default_result("/api/v1/transient/electric")
        self.assertEqual(result["summary"]["level_id"], "mid")
        sweep = result["series"]["power_sweep"]["columns"]
        self.assertEqual(sweep["level_index"][0], -1)
        self.assertIsNone(sweep["specific_impulse_s"][0])
        self.assertEqual(sweep["level_index"][-1], 2)

    def test_electric_higher_power_and_invalid_table(self) -> None:
        result = dispatch_calculation("/api/v1/transient/electric", {"available_power_w": 3500})
        self.assertEqual(result["summary"]["level_id"], "high")
        with self.assertRaises(InputError):
            dispatch_calculation("/api/v1/transient/electric", {"levels": []})
        with self.assertRaises(DomainError):
            dispatch_calculation(
                "/api/v1/transient/electric",
                {"levels": [{"level_id": "x", "input_power_w": 1, "thrust_n": 1,
                             "mass_flow_kg_s": 1e-9}]},
            )

    def test_pulse_mode_default_and_request(self) -> None:
        result = default_result("/api/v1/transient/pulse-mode")
        self.assertAlmostEqual(result["summary"]["minimum_impulse_bit_n_s"], 0.006)
        self.assertEqual(len(result["quantization"]), 4)
        self.assertAlmostEqual(
            result["summary"]["pulse_train"]["free_space_delta_v_m_s"], 0.05, places=6
        )
        custom = dispatch_calculation(
            "/api/v1/transient/pulse-mode",
            {"delta_v_m_s": 0.01, "quantization_requests_m_s": [0.001], "sweep_points": 5},
        )
        self.assertAlmostEqual(custom["summary"]["plan"]["requested_m_s"], 0.01)
        self.assertEqual(len(custom["quantization"]), 1)
        self.assertEqual(custom["series"]["impulse_bit_curve"]["point_count"], 5)
        with self.assertRaises(DomainError):
            dispatch_calculation("/api/v1/transient/pulse-mode", {"delta_v_m_s": 1.0e6})
        with self.assertRaises(InputError):
            dispatch_calculation("/api/v1/transient/pulse-mode", {"on_times_s": "fast"})

    def test_convergence_failure_propagates_as_convergence_error(self) -> None:
        failure = ConvergenceError("metrics unresolved")
        with patch.object(transient_routes, "free_space_burn_metrics", side_effect=failure), \
                self.assertRaises(ConvergenceError):
            dispatch_calculation("/api/v1/transient/pulse-mode", {})


class EnsembleRouteTests(unittest.TestCase):
    def test_uncertainty_default_statistics_and_histograms(self) -> None:
        result = default_result("/api/v1/transient/uncertainty")
        self.assertEqual(result["summary"]["sample_count"], 16)
        self.assertEqual(result["summary"]["seed"], 42)
        metrics = {item["metric"] for item in result["statistics"]}
        self.assertIn("delivered_impulse_n_s", metrics)
        for histogram in result["histograms"].values():
            self.assertEqual(sum(histogram["counts"]), 16)
            self.assertEqual(len(histogram["edges"]), len(histogram["counts"]) + 1)
        self.assertEqual(result["series"]["realizations"]["point_count"], 16)

    def test_uncertainty_is_seed_deterministic(self) -> None:
        payload = {**FAST_ENSEMBLE, "seed": 7}
        first = dispatch_calculation("/api/v1/transient/uncertainty", payload)
        second = dispatch_calculation("/api/v1/transient/uncertainty", payload)
        self.assertEqual(first, second)
        other = dispatch_calculation("/api/v1/transient/uncertainty", {**payload, "seed": 8})
        self.assertNotEqual(
            first["series"]["realizations"]["columns"]["delivered_impulse_n_s"],
            other["series"]["realizations"]["columns"]["delivered_impulse_n_s"],
        )

    def test_uncertainty_rejects_bad_requests(self) -> None:
        with self.assertRaises(InputError):
            dispatch_calculation("/api/v1/transient/uncertainty", {"sample_count": 10_000})
        with self.assertRaises(InputError):
            dispatch_calculation(
                "/api/v1/transient/uncertainty", {"integration_relative_tolerance": 1e-9}
            )
        bad_parameter = {
            "name": "wall_thickness_m",
            "nominal": 1.0,
            "lower": 0.9,
            "upper": 1.1,
            "distribution": "uniform",
            "uncertainty_class": "epistemic",
            "source": "test",
            "rationale": "test",
        }
        with self.assertRaises(DomainError):
            dispatch_calculation(
                "/api/v1/transient/uncertainty", {**FAST_ENSEMBLE, "parameters": [bad_parameter]}
            )
        with self.assertRaises(InputError):
            dispatch_calculation(
                "/api/v1/transient/uncertainty",
                {"parameters": [{**bad_parameter, "distribution": "bimodal"}]},
            )

    def test_scenarios_default_and_custom_tree(self) -> None:
        result = default_result("/api/v1/transient/scenarios")
        self.assertAlmostEqual(result["summary"]["zero_thrust_probability"], 0.02)
        impulse = next(
            item for item in result["statistics"] if item["metric"] == "delivered_impulse_n_s"
        )
        self.assertAlmostEqual(impulse["expected_value"], 1880.0)
        self.assertEqual(set(result["scenario_histories"]), {row["scenario_id"]
                                                             for row in result["scenarios"]})
        derated = {
            "scenarios": [
                {"scenario_id": "nominal", "probability": 0.95, "events": [],
                 "source": "test", "rationale": "nominal"},
                {"scenario_id": "derated", "probability": 0.05, "source": "test",
                 "rationale": "derated engine",
                 "events": [{"event_type": "thrust-scale", "engine_ids": ["engine-a"],
                             "value": 0.9, "source": "test", "rationale": "low thrust"}]},
            ]
        }
        custom = dispatch_calculation("/api/v1/transient/scenarios", derated)
        impulse = next(
            item for item in custom["statistics"] if item["metric"] == "delivered_impulse_n_s"
        )
        self.assertAlmostEqual(impulse["expected_value"], 0.95 * 2000.0 + 0.05 * 1900.0)
        not_exhaustive = {"scenarios": [{**derated["scenarios"][0], "probability": 0.5}]}
        with self.assertRaises(DomainError):
            dispatch_calculation("/api/v1/transient/scenarios", not_exhaustive)
        with self.assertRaises(InputError):
            dispatch_calculation(
                "/api/v1/transient/scenarios",
                {"scenarios": [{**derated["scenarios"][0], "probability": 1.0, "source": ""}]},
            )


class TimelineCalibrationSassRouteTests(unittest.TestCase):
    def test_timeline_default_mass_ledger_and_events(self) -> None:
        result = default_result("/api/v1/transient/timeline")
        summary = result["summary"]
        self.assertAlmostEqual(summary["final_vehicle_mass_kg"], 6000.0 - 1200.0 - 150.0 - 400.0)
        kinds = {event["kind"] for event in result["events"]}
        self.assertIn("jettison", kinds)
        history = result["series"]["history"]["columns"]
        self.assertEqual(max(history["axial_thrust_n"]), 30_000.0)
        self.assertEqual(history["time_s"], sorted(history["time_s"]))

    def test_timeline_single_burn_and_overlap(self) -> None:
        burn = dict(transient_routes._TIMELINE_BURNS[0])
        result = dispatch_calculation(
            "/api/v1/transient/timeline", {"burns": [burn], "jettisons": [], "max_points": 50}
        )
        self.assertAlmostEqual(result["summary"]["final_vehicle_mass_kg"], 4800.0)
        self.assertEqual(result["series"]["history"]["point_count"], 50)
        overlapping = [burn, {**burn, "burn_id": "second", "ignition_time_s": 20.0}]
        with self.assertRaises(DomainError):
            dispatch_calculation("/api/v1/transient/timeline", {"burns": overlapping})

    def test_calibration_default_and_subset(self) -> None:
        result = default_result("/api/v1/transient/calibration")
        self.assertAlmostEqual(result["summary"]["consensus_scale"], 1.008, places=9)
        self.assertEqual(result["summary"]["engine_count"], 6)
        population = result["residual_parameters"]["population"]
        self.assertLess(population["lower"], 1.0)
        self.assertGreater(population["upper"], 1.0)
        subset = [dict(item) for item in transient_routes._CALIBRATION_MEASUREMENTS[:12]]
        custom = dispatch_calculation(
            "/api/v1/transient/calibration",
            {"measurements": subset, "parameter_name": "specific_impulse_scale"},
        )
        self.assertEqual(custom["summary"]["engine_count"], 4)
        self.assertEqual(custom["series"]["measurements"]["point_count"], 12)
        with self.assertRaises(DomainError):
            dispatch_calculation(
                "/api/v1/transient/calibration", {"measurements": subset[:1]}
            )

    def test_sass_default_override_and_document(self) -> None:
        result = default_result("/api/v1/sass/mission")
        self.assertEqual(result["summary"]["candidate_count"], 5)
        self.assertEqual(result["summary"]["eligible_count"], 3)
        self.assertNotIn("document", result)
        heavier = dispatch_calculation(
            "/api/v1/sass/mission", {"bus_dry_mass_kg": 60.0, "include_document": True}
        )
        base_mass = result["series"]["trade"]["columns"]["wet_mass_kg"]
        new_mass = heavier["series"]["trade"]["columns"]["wet_mass_kg"]
        for before, after in zip(base_mass, new_mass, strict=True):
            if before is not None:
                self.assertGreater(after, before)
        self.assertEqual(heavier["document"]["study_hash"], heavier["summary"]["study_hash"])
        self.assertNotEqual(heavier["summary"]["study_hash"], result["summary"]["study_hash"])
        with self.assertRaises(InputError):
            dispatch_calculation("/api/v1/sass/mission", {"study": []})
        with self.assertRaises(DomainError):
            dispatch_calculation("/api/v1/sass/mission", {"study": {"schema": "other"}})


class QuietHandler(RocketPropulsionHandler):
    def log_message(self, message_format: str, *args: object) -> None:
        pass


class TransientHttpRoundTripTests(unittest.TestCase):
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

    def post(self, path: str, payload: object) -> tuple[int, dict]:
        request = Request(
            self.base_url + path,
            data=json.dumps(payload).encode("utf-8"),
            headers={"Content-Type": "application/json"},
            method="POST",
        )
        try:
            with urlopen(request, timeout=10.0) as response:
                return response.status, json.load(response)
        except HTTPError as error:
            return error.code, json.load(error)

    def test_empty_request_round_trip(self) -> None:
        status, body = self.post("/api/v1/transient/electric", {})
        self.assertEqual(status, 200)
        self.assertEqual(body["data"], json.loads(json.dumps(default_result(
            "/api/v1/transient/electric"))))
        self.assertEqual(body["meta"]["model"], ROUTE_METADATA[
            "/api/v1/transient/electric"].model)
        sweep = body["data"]["series"]["power_sweep"]
        self.assertTrue(all(len(column) == sweep["point_count"]
                            for column in sweep["columns"].values()))
        self.assertTrue(all(value is None or math.isfinite(value)
                            for value in sweep["columns"]["thrust_n"]))

    def test_error_statuses(self) -> None:
        status, body = self.post("/api/v1/transient/calibration", {"maximum_iterations": "x"})
        self.assertEqual(status, 400)
        self.assertEqual(body["error"]["field"], "maximum_iterations")
        status, body = self.post("/api/v1/transient/scenarios", {"scenarios": [
            {"scenario_id": "only", "probability": 0.5, "events": [], "source": "t",
             "rationale": "t"}]})
        self.assertEqual(status, 400)
        failure = ConvergenceError("metrics unresolved")
        with patch.object(transient_routes, "free_space_burn_metrics", side_effect=failure):
            status, body = self.post("/api/v1/transient/pulse-mode", {})
        self.assertEqual(status, 422)
        self.assertEqual(body["error"]["code"], failure.code)
        status, body = self.post("/api/v1/transient/not-real", {})
        self.assertEqual(status, 404)


if __name__ == "__main__":
    unittest.main()
