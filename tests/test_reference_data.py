"""Regression against source-labelled external reference cases."""

from __future__ import annotations

import json
import unittest
from pathlib import Path

from rocket_propulsion.compressible import calculate_isentropic, calculate_normal_shock
from rocket_propulsion.geometry import generate_moc_nozzle
from rocket_propulsion.propulsion import calculate_loss_budget
from rocket_propulsion.studies import workspace_from_json

REFERENCE_FILE = Path(__file__).parent / "reference_data" / "ideal_cases.json"
VALIDATION_FILE = Path(__file__).parent / "reference_data" / "validation_matrix.json"
EXAMPLE_WORKSPACE = Path(__file__).parents[1] / "examples" / "two_case.rplab.json"


class ReferenceDataTests(unittest.TestCase):
    def test_reference_cases(self) -> None:
        document = json.loads(REFERENCE_FILE.read_text(encoding="utf-8"))
        self.assertTrue(document["source"])
        calculators = {
            "isentropic_mach_2_air": lambda inputs: calculate_isentropic(**inputs),
            "normal_shock_mach_2_air": lambda inputs: calculate_normal_shock(**inputs),
        }
        for case in document["cases"]:
            with self.subTest(case=case["name"]):
                result = calculators[case["name"]](case["inputs"])
                tolerance = case["relative_tolerance"]
                for field, expected in case["expected"].items():
                    self.assertAlmostEqual(
                        getattr(result, field),
                        expected,
                        delta=tolerance * max(1.0, abs(expected)),
                    )

    def test_validation_matrix_automated_invariants(self) -> None:
        document = json.loads(VALIDATION_FILE.read_text(encoding="utf-8"))
        cases = {case["id"]: case for case in document["cases"]}
        moc_case = cases["moc-boundary-compatibility"]
        design = generate_moc_nozzle(
            throat_area_m2=0.01,
            area_ratio=25.0,
            gamma=1.22,
            characteristic_count=20,
        )
        expected = moc_case["expected"]
        self.assertLess(design.residuals.wall_tangency_max_deg, expected["wall_tangency_max_deg"])
        self.assertLess(
            design.residuals.centerline_symmetry_max_deg,
            expected["centerline_symmetry_max_deg"],
        )
        self.assertLess(design.residuals.exit_mach_absolute, expected["exit_mach_absolute"])
        self.assertLess(
            design.refined_length_change_relative,
            expected["refined_length_change_relative"],
        )
        losses = calculate_loss_budget(
            ideal_thrust_n=100_000.0, ideal_specific_impulse_s=450.0
        )
        loss_expected = cases["loss-budget-conservation"]["expected"]
        self.assertLessEqual(
            losses.thrust_balance_residual_n,
            loss_expected["thrust_balance_residual_n"],
        )
        self.assertLessEqual(
            losses.specific_impulse_balance_residual_s,
            loss_expected["specific_impulse_balance_residual_s"],
        )

    def test_example_workspace_is_current_and_portable(self) -> None:
        workspace = workspace_from_json(EXAMPLE_WORKSPACE.read_text(encoding="utf-8"))
        self.assertEqual(workspace.schema_version, 2)
        self.assertEqual(len(workspace.cases), 2)

    def test_cea_release_fixture_remains_explicitly_external(self) -> None:
        document = json.loads(VALIDATION_FILE.read_text(encoding="utf-8"))
        cea_case = next(case for case in document["cases"] if case["id"].startswith("cea-"))
        self.assertEqual(cea_case["status"], "recorded-external-gate")
        self.assertGreater(cea_case["expected"]["chamber_temperature_k"], 3000.0)


if __name__ == "__main__":
    unittest.main()
