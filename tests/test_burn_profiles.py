"""Analytic, persistence, and adapter tests for L1 throttle schedules."""

import json
import unittest
from contextlib import redirect_stdout
from io import StringIO
from math import log
from pathlib import Path
from tempfile import TemporaryDirectory

from rocket_propulsion.core.errors import DomainError, InputError
from rocket_propulsion.integrations import export_sidera_artifact
from rocket_propulsion.propulsion.burns import (
    BurnDefinition,
    BurnTarget,
    BurnTargetKind,
    ConstantPerformanceProvider,
    CutoffReason,
    ProfiledPerformanceProvider,
    PropellantRole,
    TankDefinition,
    ThrottleSchedule,
    ThrottleSegment,
    build_burn_html_report,
    burn_result_from_json,
    burn_result_to_dict,
    burn_result_to_json,
    integrated_throttle,
    reduce_to_constant,
    simulate_burn,
    time_for_integrated_throttle,
)
from rocket_propulsion.propulsion.burns.cli import main as burn_cli_main
from rocket_propulsion.studies import (
    StudyCase,
    create_workspace,
    workspace_from_json,
    workspace_to_json,
)


def monopropellant_point():
    return ConstantPerformanceProvider.monopropellant(
        delivered_thrust_n=980.665,
        system_specific_impulse_s=100.0,
        chamber_pressure_pa=1_000_000.0,
        propellant_key="hydrazine",
    ).operating_point()


def definition(target: BurnTarget, *, loaded_mass_kg: float = 20.0, reserve_kg: float = 0.0):
    return BurnDefinition(
        name="profiled monopropellant burn",
        initial_mass_kg=100.0,
        protected_dry_mass_kg=80.0,
        tanks=(
            TankDefinition(
                "propellant",
                "hydrazine",
                PropellantRole.MONOPROPELLANT,
                loaded_mass_kg,
                reserve_kg,
            ),
        ),
        target=target,
        maximum_duration_s=100.0,
    )


def triangular_schedule() -> ThrottleSchedule:
    return ThrottleSchedule(
        name="two-second ramps",
        segments=(
            ThrottleSegment(2.0, 0.0, 1.0, "ramp-up"),
            ThrottleSegment(3.0, 1.0, 1.0, "steady"),
            ThrottleSegment(2.0, 1.0, 0.0, "ramp-down"),
        ),
    )


class ProfileScheduleTests(unittest.TestCase):
    def test_constant_profile_is_proven_exact_for_sidera(self) -> None:
        schedule = ThrottleSchedule((ThrottleSegment(2.0, 0.5, 0.5, "half"),))
        result = simulate_burn(
            definition(BurnTarget(BurnTargetKind.DURATION, 2.0)),
            ProfiledPerformanceProvider(monopropellant_point(), schedule),
        )
        self.assertTrue(result.summary.exact_constant_sidera_export_available)
        artifact = export_sidera_artifact(
            result,
            t_start_s=0.0,
            direction=(1.0, 0.0, 0.0),
        )
        self.assertEqual(artifact["manifest"]["mapping"], "exact_constant")
        self.assertAlmostEqual(artifact["maneuvers"][0]["thrust_n"], 490.3325)
        self.assertNotIn("approximation_warning", artifact["manifest"])

    def test_triangle_integrals_mass_impulse_centroid_and_delta_v(self) -> None:
        schedule = triangular_schedule()
        result = simulate_burn(
            definition(BurnTarget(BurnTargetKind.DURATION, 7.0)),
            ProfiledPerformanceProvider(monopropellant_point(), schedule),
        )

        self.assertEqual(result.schema, "rocket_propulsion_burn_v2")
        self.assertAlmostEqual(schedule.total_exposure_s, 5.0)
        self.assertAlmostEqual(result.summary.consumed_propellant_kg, 5.0)
        self.assertAlmostEqual(result.summary.delivered_total_impulse_n_s, 4_903.325)
        self.assertAlmostEqual(result.summary.final_mass_kg, 95.0)
        self.assertAlmostEqual(result.summary.thrust_centroid_time_s, 3.5)
        self.assertAlmostEqual(
            result.summary.ideal_delta_v_m_s,
            980.665 * log(100.0 / 95.0),
        )
        self.assertEqual(result.evidence.segment_count, 3)
        self.assertFalse(result.summary.exact_constant_sidera_export_available)
        self.assertLess(result.summary.mass_closure_error_kg, 1e-12)
        self.assertLess(result.summary.impulse_closure_error_n_s, 1e-12)

    def test_exposure_inverse_crosses_a_throttle_step_exactly(self) -> None:
        schedule = ThrottleSchedule(
            segments=(
                ThrottleSegment(1.0, 0.5, 0.5, "half"),
                ThrottleSegment(2.0, 1.0, 1.0, "full"),
            )
        )
        self.assertAlmostEqual(integrated_throttle(schedule, 2.0), 1.5)
        self.assertAlmostEqual(time_for_integrated_throttle(schedule, 1.5) or 0.0, 2.0)
        result = simulate_burn(
            definition(BurnTarget(BurnTargetKind.TOTAL_IMPULSE, 1_470.9975)),
            ProfiledPerformanceProvider(monopropellant_point(), schedule),
        )
        self.assertAlmostEqual(result.summary.duration_s, 2.0)
        self.assertTrue(result.summary.target_achieved)

    def test_open_cycle_bipropellant_streams_scale_with_profile(self) -> None:
        provider = ConstantPerformanceProvider.bipropellant(
            delivered_thrust_n=190_000.0,
            system_specific_impulse_s=190_000.0 / (9.80665 * 51.0),
            oxidizer_fuel_ratio=40.0 / 11.0,
            chamber_specific_impulse_s=190_000.0 / (9.80665 * 50.0),
            dump_mass_flow_kg_s=1.0,
            dump_role=PropellantRole.FUEL,
            chamber_pressure_pa=7_000_000.0,
            oxidizer_key="lox",
            fuel_key="rp1",
        )
        burn_definition = BurnDefinition(
            name="profiled gas-generator burn",
            initial_mass_kg=1000.0,
            protected_dry_mass_kg=800.0,
            tanks=(
                TankDefinition("oxidizer", "lox", PropellantRole.OXIDIZER, 100.0),
                TankDefinition("fuel", "rp1", PropellantRole.FUEL, 30.0),
            ),
            target=BurnTarget(BurnTargetKind.DURATION, 2.0),
            maximum_duration_s=10.0,
        )
        schedule = ThrottleSchedule((ThrottleSegment(2.0, 0.5, 0.5, "half"),))
        result = simulate_burn(
            burn_definition,
            ProfiledPerformanceProvider(provider.operating_point(), schedule),
        )
        self.assertAlmostEqual(result.summary.consumed_propellant_kg, 51.0)
        self.assertAlmostEqual(dict(result.summary.consumed_by_tank_kg)["oxidizer"], 40.0)
        self.assertAlmostEqual(dict(result.summary.consumed_by_tank_kg)["fuel"], 11.0)

    def test_schedule_end_and_tank_reserve_are_typed_cutoffs(self) -> None:
        short = ThrottleSchedule((ThrottleSegment(1.0, 1.0, 1.0, "pulse"),))
        ended = simulate_burn(
            definition(BurnTarget(BurnTargetKind.DURATION, 5.0)),
            ProfiledPerformanceProvider(monopropellant_point(), short),
        )
        self.assertFalse(ended.summary.target_achieved)
        self.assertEqual(ended.summary.cutoff_reason, CutoffReason.SCHEDULE_END)

        ramp = ThrottleSchedule((ThrottleSegment(10.0, 0.0, 1.0, "ramp"),))
        limited = simulate_burn(
            definition(
                BurnTarget(BurnTargetKind.DURATION, 10.0),
                loaded_mass_kg=3.0,
                reserve_kg=1.0,
            ),
            ProfiledPerformanceProvider(monopropellant_point(), ramp),
        )
        self.assertEqual(limited.summary.cutoff_reason, CutoffReason.TANK_RESERVE)
        self.assertEqual(limited.summary.binding_constraint, "propellant")
        self.assertAlmostEqual(limited.summary.duration_s, (40.0) ** 0.5)
        self.assertAlmostEqual(dict(limited.final_state.tank_masses_kg)["propellant"], 1.0)

    def test_v2_round_trip_report_and_exact_sidera_refusal(self) -> None:
        result = simulate_burn(
            definition(BurnTarget(BurnTargetKind.DURATION, 7.0)),
            ProfiledPerformanceProvider(monopropellant_point(), triangular_schedule()),
        )
        encoded = burn_result_to_json(result)
        restored = burn_result_from_json(encoded)
        self.assertEqual(restored, result)
        self.assertIn("Throttle schedule", build_burn_html_report(restored))

        tampered = json.loads(encoded)
        tampered["schedule"]["segments"][0]["end_throttle"] = 0.9
        with self.assertRaises(InputError):
            burn_result_from_json(json.dumps(tampered))
        with self.assertRaises(DomainError):
            export_sidera_artifact(
                restored,
                t_start_s=0.0,
                direction=(1.0, 0.0, 0.0),
            )
        reduction = reduce_to_constant(restored)
        self.assertFalse(reduction.exact)
        self.assertAlmostEqual(reduction.equivalent_system_specific_impulse_s, 100.0)
        self.assertLess(reduction.impulse_residual_n_s, 1e-12)
        approximate = export_sidera_artifact(
            restored,
            t_start_s=0.0,
            direction=(1.0, 0.0, 0.0),
            mode="equivalent_constant",
        )
        self.assertEqual(approximate["manifest"]["mapping"], "equivalent_constant")
        self.assertIn("approximation_warning", approximate["manifest"])
        self.assertAlmostEqual(
            approximate["maneuvers"][0]["thrust_n"] * restored.summary.duration_s,
            restored.summary.delivered_total_impulse_n_s,
        )

    def test_l0_artifact_shape_remains_schema_v1(self) -> None:
        result = simulate_burn(
            definition(BurnTarget(BurnTargetKind.DURATION, 1.0)),
            ConstantPerformanceProvider(monopropellant_point()),
        )
        payload = burn_result_to_dict(result)
        self.assertEqual(payload["schema"], "rocket_propulsion_burn_v1")
        self.assertNotIn("schedule", payload)

    def test_profiled_artifact_survives_workspace_round_trip(self) -> None:
        result = simulate_burn(
            definition(BurnTarget(BurnTargetKind.DURATION, 7.0)),
            ProfiledPerformanceProvider(monopropellant_point(), triangular_schedule()),
        )
        case = StudyCase(
            case_id="l1-profile",
            name="L1 profile",
            model="profiled-propulsion-burn",
            inputs={"schedule": "ramp-steady-ramp"},
            results={"artifact": burn_result_to_dict(result)},
            model_versions={"rocket_propulsion_burn_v2": "2"},
        )
        encoded = workspace_to_json(create_workspace("Profile study", (case,)))
        restored_workspace = workspace_from_json(encoded)
        restored = burn_result_from_json(
            json.dumps(restored_workspace.cases[0].results["artifact"])
        )
        self.assertEqual(restored, result)
        self.assertEqual(workspace_to_json(restored_workspace), encoded)

    def test_profiled_cli_simulate_validate_report_and_reduced_export(self) -> None:
        example = (
            Path(__file__).resolve().parents[1] / "examples" / "profiled_burn_input.json"
        )
        with TemporaryDirectory() as directory:
            result_path = Path(directory) / "profile-result.json"
            report_path = Path(directory) / "profile-report.html"
            sidera_path = Path(directory) / "sidera-equivalent.json"
            with redirect_stdout(StringIO()):
                simulate_status = burn_cli_main(
                    [
                        "simulate",
                        str(example),
                        "--output",
                        str(result_path),
                        "--report",
                        str(report_path),
                    ]
                )
                validate_status = burn_cli_main(["validate", str(result_path)])
                export_status = burn_cli_main(
                    [
                        "export",
                        str(result_path),
                        "--target",
                        "sidera",
                        "--output",
                        str(sidera_path),
                        "--start",
                        "100",
                        "--direction",
                        "1",
                        "0",
                        "0",
                        "--mode",
                        "equivalent_constant",
                    ]
                )
            self.assertEqual(simulate_status, 0)
            self.assertEqual(validate_status, 0)
            self.assertEqual(export_status, 0)
            self.assertIn("Throttle schedule", report_path.read_text(encoding="utf-8"))
            artifact = json.loads(sidera_path.read_text(encoding="utf-8"))
            self.assertEqual(artifact["manifest"]["mapping"], "equivalent_constant")


if __name__ == "__main__":
    unittest.main()
