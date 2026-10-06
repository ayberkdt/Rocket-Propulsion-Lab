"""Tests for workspaces, uncertainty, staging, and cycle balances."""

import json
import unittest

from rocket_propulsion.propulsion import (
    StageDefinition,
    calculate_engine_cycle,
    optimize_staging,
)
from rocket_propulsion.propulsion.burns import (
    BurnDefinition,
    BurnTarget,
    BurnTargetKind,
    ConstantPerformanceProvider,
    PropellantRole,
    TankDefinition,
    burn_result_from_dict,
    burn_result_to_dict,
    simulate_burn,
)
from rocket_propulsion.studies import (
    ParameterRange,
    StudyCase,
    build_report_bundle,
    compare_cases,
    create_workspace,
    duplicate_case,
    propagate_uncertainty,
    workspace_from_json,
    workspace_to_json,
)


class WorkspaceTests(unittest.TestCase):
    def setUp(self) -> None:
        self.case = StudyCase(
            case_id="sea-level",
            name="Sea level",
            model="ideal-nozzle",
            inputs={"ambient_pressure_pa": 101325.0},
            results={"thrust_n": 95_000.0, "specific_impulse_s": 320.0},
            model_versions={"ideal-nozzle": "1.0"},
            warnings=("screening result",),
        )

    def test_round_trip_is_numerically_and_textually_stable(self) -> None:
        workspace = create_workspace("Demo", (self.case,))
        encoded = workspace_to_json(workspace)
        restored = workspace_from_json(encoded)
        self.assertEqual(restored, workspace)
        self.assertEqual(workspace_to_json(restored), encoded)

    def test_v1_migration_duplicate_and_case_comparison(self) -> None:
        legacy = json.dumps(
            {
                "version": 1,
                "name": "Legacy",
                "case": {
                    "case_id": "old",
                    "name": "Old",
                    "model": "nozzle",
                    "inputs": {},
                    "results": {"thrust_n": 100.0},
                },
            }
        )
        workspace = workspace_from_json(legacy)
        copied = duplicate_case(workspace.cases[0], case_id="new", name="New")
        changed = StudyCase(
            copied.case_id,
            copied.name,
            copied.model,
            copied.inputs,
            {"thrust_n": 110.0},
        )
        comparison = compare_cases(workspace.cases[0], changed)
        self.assertEqual(workspace.schema_version, 2)
        self.assertAlmostEqual(comparison.differences[0].relative_difference or 0.0, 0.1)

    def test_report_bundle_carries_versions_and_warnings(self) -> None:
        report = build_report_bundle(create_workspace("Demo", (self.case,)))
        self.assertIn("model_versions", report.json)
        self.assertIn("specific_impulse_s", report.csv)
        self.assertIn("ideal-nozzle: 1.0", report.html)
        self.assertIn("screening result", report.html)
        self.assertIn("<svg", report.svg)

    def test_burn_artifact_survives_workspace_round_trip(self) -> None:
        provider = ConstantPerformanceProvider.bipropellant(
            delivered_thrust_n=24_000.0,
            system_specific_impulse_s=328.0,
            oxidizer_fuel_ratio=3.4,
            ideal_thrust_n=25_000.0,
            chamber_pressure_pa=7_000_000.0,
        )
        definition = BurnDefinition(
            name="workspace burn",
            initial_mass_kg=2_400.0,
            protected_dry_mass_kg=1_400.0,
            tanks=(
                TankDefinition(
                    "oxidizer", "oxidizer", PropellantRole.OXIDIZER, 760.0, 20.0
                ),
                TankDefinition("fuel", "fuel", PropellantRole.FUEL, 240.0, 8.0),
            ),
            target=BurnTarget(BurnTargetKind.IDEAL_DELTA_V, 420.0),
            maximum_duration_s=300.0,
        )
        result = simulate_burn(definition, provider)
        case = StudyCase(
            case_id="burn-420",
            name="420 m/s burn",
            model="standalone-constant-burn",
            inputs={"target_kind": BurnTargetKind.IDEAL_DELTA_V.value},
            results={
                "artifact": burn_result_to_dict(result),
                "duration_s": result.summary.duration_s,
                "ideal_delta_v_m_s": result.summary.ideal_delta_v_m_s,
            },
            model_versions={"rocket_propulsion_burn_v1": "1"},
        )

        encoded = workspace_to_json(create_workspace("Burn study", (case,)))
        restored = workspace_from_json(encoded)
        restored_result = burn_result_from_dict(restored.cases[0].results["artifact"])

        self.assertEqual(restored_result, result)
        self.assertEqual(restored_result.result_hash, result.result_hash)
        self.assertEqual(workspace_to_json(restored), encoded)


class SamplingTests(unittest.TestCase):
    def test_seeded_lhs_is_reproducible_and_ranks_dominant_input(self) -> None:
        parameters = (
            ParameterRange("strong", 0.0, 1.0),
            ParameterRange("weak", 0.0, 1.0),
        )
        first = propagate_uncertainty(
            parameters,
            evaluator=lambda sample: 10.0 * sample["strong"] + sample["weak"],
            sample_count=100,
            seed=123,
        )
        second = propagate_uncertainty(
            parameters,
            evaluator=lambda sample: 10.0 * sample["strong"] + sample["weak"],
            sample_count=100,
            seed=123,
        )
        self.assertEqual(first, second)
        self.assertEqual(first.sensitivity[0].parameter, "strong")
        self.assertLess(first.output_percentiles["p05"], first.output_percentiles["p95"])


class SystemAnalysisTests(unittest.TestCase):
    def test_staging_optimizer_meets_target_and_refines(self) -> None:
        result = optimize_staging(
            stages=(
                StageDefinition("Booster", 300.0, 0.10),
                StageDefinition("Upper", 450.0, 0.08),
            ),
            target_delta_v_m_s=9400.0,
            payload_mass_kg=1000.0,
            resolution=120,
        )
        self.assertAlmostEqual(result.achieved_delta_v_m_s, 9400.0)
        self.assertLess(result.convergence_delta_initial_mass_relative, 0.01)
        self.assertGreater(result.initial_mass_kg, result.payload_mass_kg)
        self.assertTrue(result.feasible)

    def test_cycle_power_balance_and_pressure_fed_limit(self) -> None:
        common = {
            "chamber_pressure_pa": 10.0e6,
            "chamber_mass_flow_kg_s": 100.0,
            "propellant_density_kg_m3": 1000.0,
            "inlet_pressure_pa": 0.3e6,
            "injector_pressure_drop_fraction": 0.2,
            "pump_efficiency": 0.7,
            "turbine_efficiency": 0.6,
            "turbine_inlet_temperature_k": 900.0,
            "turbine_specific_heat_j_kg_k": 2500.0,
            "turbine_gamma": 1.25,
            "turbine_pressure_ratio": 8.0,
        }
        gas_generator = calculate_engine_cycle(cycle="gas-generator", **common)
        pressure_fed = calculate_engine_cycle(cycle="pressure-fed", **common)
        self.assertAlmostEqual(gas_generator.power_balance_residual_w, 0.0)
        self.assertGreater(gas_generator.turbine_mass_flow_kg_s, 0.0)
        self.assertLess(gas_generator.delivered_mass_fraction, 1.0)
        self.assertEqual(pressure_fed.pump_power_w, 0.0)
        self.assertEqual(pressure_fed.required_feed_pressure_pa, 12.0e6)


if __name__ == "__main__":
    unittest.main()
