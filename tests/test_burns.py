"""Analytic, conservation, cutoff, and persistence tests for L0 burns."""

import json
import unittest
from contextlib import redirect_stdout
from io import StringIO
from math import exp, log
from pathlib import Path
from tempfile import TemporaryDirectory

from rocket_propulsion.core.errors import DomainError, InputError
from rocket_propulsion.propulsion import STANDARD_GRAVITY_M_S2
from rocket_propulsion.propulsion.burns import (
    BurnDefinition,
    BurnTarget,
    BurnTargetKind,
    ConstantPerformanceProvider,
    CutoffReason,
    EngineOperatingPoint,
    FlowDestination,
    PropellantFlow,
    PropellantRole,
    TankDefinition,
    build_burn_html_report,
    burn_result_from_json,
    burn_result_to_json,
    simulate_burn,
    simulate_constant_burn,
)
from rocket_propulsion.propulsion.burns.cli import main as burn_cli_main


def example_operating_point() -> EngineOperatingPoint:
    """Return a bipropellant point with an explicit gas-generator dump."""

    delivered_thrust_n = 190_000.0
    total_flow_kg_s = 51.0
    return EngineOperatingPoint(
        name="reference gas-generator point",
        commanded_throttle=1.0,
        realized_throttle=1.0,
        active_engine_count=1,
        ideal_thrust_n=200_000.0,
        delivered_thrust_n=delivered_thrust_n,
        chamber_specific_impulse_s=400.0,
        system_specific_impulse_s=(
            delivered_thrust_n / (STANDARD_GRAVITY_M_S2 * total_flow_kg_s)
        ),
        chamber_mass_flow_kg_s=50.0,
        total_tank_flow_kg_s=total_flow_kg_s,
        oxidizer_mass_flow_kg_s=40.0,
        fuel_mass_flow_kg_s=11.0,
        mixture_ratio=40.0 / 11.0,
        chamber_pressure_pa=7_000_000.0,
        ambient_pressure_pa=0.0,
        streams=(
            PropellantFlow(
                "main-oxidizer",
                "oxidizer",
                "lox",
                PropellantRole.OXIDIZER,
                FlowDestination.CHAMBER,
                40.0,
            ),
            PropellantFlow(
                "main-fuel",
                "fuel",
                "rp1",
                PropellantRole.FUEL,
                FlowDestination.CHAMBER,
                10.0,
            ),
            PropellantFlow(
                "gas-generator-fuel",
                "fuel",
                "rp1",
                PropellantRole.FUEL,
                FlowDestination.GAS_GENERATOR_DUMP,
                1.0,
                contributes_to_delivered_thrust=False,
            ),
        ),
        model_id="test_point_v1",
        source="analytic fixture",
    )


def example_definition(target: BurnTarget, *, cant_efficiency: float = 1.0) -> BurnDefinition:
    return BurnDefinition(
        name="reference burn",
        initial_mass_kg=1_000.0,
        protected_dry_mass_kg=800.0,
        tanks=(
            TankDefinition("oxidizer", "lox", PropellantRole.OXIDIZER, 100.0, 5.0),
            TankDefinition("fuel", "rp1", PropellantRole.FUEL, 30.0, 5.0),
        ),
        target=target,
        maximum_duration_s=100.0,
        cant_efficiency=cant_efficiency,
    )


class ConstantBurnTests(unittest.TestCase):
    def test_bipropellant_provider_closes_open_cycle_streams(self) -> None:
        provider = ConstantPerformanceProvider.bipropellant(
            delivered_thrust_n=190_000.0,
            system_specific_impulse_s=190_000.0 / (STANDARD_GRAVITY_M_S2 * 51.0),
            oxidizer_fuel_ratio=40.0 / 11.0,
            ideal_thrust_n=200_000.0,
            chamber_specific_impulse_s=190_000.0 / (STANDARD_GRAVITY_M_S2 * 50.0),
            dump_mass_flow_kg_s=1.0,
            dump_role=PropellantRole.FUEL,
            chamber_pressure_pa=7_000_000.0,
            oxidizer_key="lox",
            fuel_key="rp1",
        )
        point = provider.operating_point()
        self.assertAlmostEqual(point.total_tank_flow_kg_s, 51.0)
        self.assertAlmostEqual(point.chamber_mass_flow_kg_s, 50.0)
        self.assertAlmostEqual(point.oxidizer_mass_flow_kg_s, 40.0)
        self.assertAlmostEqual(point.fuel_mass_flow_kg_s, 11.0)
        self.assertEqual(point.streams[-1].destination, FlowDestination.GAS_GENERATOR_DUMP)
        result = simulate_burn(
            example_definition(BurnTarget(BurnTargetKind.DURATION, 0.5)), provider
        )
        self.assertAlmostEqual(result.summary.consumed_propellant_kg, 25.5)

    def test_monopropellant_provider_builds_single_closed_stream(self) -> None:
        provider = ConstantPerformanceProvider.monopropellant(
            delivered_thrust_n=980.665,
            system_specific_impulse_s=100.0,
            chamber_pressure_pa=1_000_000.0,
            propellant_key="hydrazine",
        )
        point = provider.operating_point()
        self.assertAlmostEqual(point.total_tank_flow_kg_s, 1.0)
        self.assertEqual(len(point.streams), 1)
        self.assertEqual(point.streams[0].role, PropellantRole.MONOPROPELLANT)

    def test_duration_target_closes_mass_impulse_and_delta_v(self) -> None:
        point = example_operating_point()
        definition = example_definition(BurnTarget(BurnTargetKind.DURATION, 1.0))
        result = simulate_constant_burn(definition, point)

        self.assertTrue(result.summary.target_achieved)
        self.assertEqual(result.summary.cutoff_reason, CutoffReason.TARGET_REACHED)
        self.assertAlmostEqual(result.summary.consumed_propellant_kg, 51.0)
        self.assertAlmostEqual(result.summary.final_mass_kg, 949.0)
        self.assertAlmostEqual(result.summary.delivered_total_impulse_n_s, 190_000.0)
        self.assertAlmostEqual(
            result.summary.ideal_delta_v_m_s,
            190_000.0 / 51.0 * log(1_000.0 / 949.0),
        )
        self.assertLess(result.summary.mass_closure_error_kg, 1e-12)
        self.assertLess(result.summary.impulse_closure_error_n_s, 1e-12)
        self.assertLess(result.summary.delta_v_closure_error_m_s, 1e-10)

    def test_system_isp_accounts_for_dump_flow(self) -> None:
        point = example_operating_point()
        self.assertLess(point.system_specific_impulse_s, point.chamber_specific_impulse_s)
        self.assertAlmostEqual(
            point.system_specific_impulse_s,
            point.delivered_thrust_n / (STANDARD_GRAVITY_M_S2 * 51.0),
        )

    def test_delta_v_inverse_recovers_exact_duration(self) -> None:
        point = example_operating_point()
        duration = 0.75
        final_mass = 1_000.0 - 51.0 * duration
        requested_delta_v = 190_000.0 / 51.0 * log(1_000.0 / final_mass)
        definition = example_definition(
            BurnTarget(BurnTargetKind.IDEAL_DELTA_V, requested_delta_v)
        )
        result = simulate_constant_burn(definition, point)
        self.assertAlmostEqual(result.summary.duration_s, duration, places=12)
        self.assertAlmostEqual(result.summary.ideal_delta_v_m_s, requested_delta_v, places=10)

    def test_cant_efficiency_reduces_delta_v_but_not_delivered_impulse(self) -> None:
        point = example_operating_point()
        straight = simulate_constant_burn(
            example_definition(BurnTarget(BurnTargetKind.DURATION, 0.5)), point
        )
        canted = simulate_constant_burn(
            example_definition(BurnTarget(BurnTargetKind.DURATION, 0.5), cant_efficiency=0.9),
            point,
        )
        self.assertEqual(
            straight.summary.delivered_total_impulse_n_s,
            canted.summary.delivered_total_impulse_n_s,
        )
        self.assertAlmostEqual(
            canted.summary.ideal_delta_v_m_s,
            straight.summary.ideal_delta_v_m_s * 0.9,
        )

    def test_unreachable_target_reports_binding_tank_reserve(self) -> None:
        point = example_operating_point()
        definition = example_definition(BurnTarget(BurnTargetKind.DURATION, 10.0))
        result = simulate_constant_burn(definition, point)

        self.assertFalse(result.summary.target_achieved)
        self.assertEqual(result.summary.cutoff_reason, CutoffReason.TANK_RESERVE)
        self.assertEqual(result.summary.binding_constraint, "fuel")
        self.assertAlmostEqual(result.summary.duration_s, 25.0 / 11.0)
        self.assertAlmostEqual(dict(result.final_state.tank_masses_kg)["fuel"], 5.0)

    def test_propellant_and_impulse_targets_use_closed_form(self) -> None:
        point = example_operating_point()
        mass_result = simulate_constant_burn(
            example_definition(BurnTarget(BurnTargetKind.PROPELLANT_MASS, 25.5)), point
        )
        impulse_result = simulate_constant_burn(
            example_definition(BurnTarget(BurnTargetKind.TOTAL_IMPULSE, 95_000.0)), point
        )
        self.assertAlmostEqual(mass_result.summary.duration_s, 0.5)
        self.assertAlmostEqual(impulse_result.summary.duration_s, 0.5)

    def test_json_round_trip_preserves_hashes_and_detects_tampering(self) -> None:
        result = simulate_constant_burn(
            example_definition(BurnTarget(BurnTargetKind.DURATION, 1.0)),
            example_operating_point(),
        )
        text = burn_result_to_json(result)
        restored = burn_result_from_json(text)
        self.assertEqual(restored, result)
        self.assertEqual(burn_result_to_json(restored), text)

        payload = json.loads(text)
        payload["summary"]["final_mass_kg"] = 900.0
        with self.assertRaises(InputError):
            burn_result_from_json(json.dumps(payload))

    def test_html_report_and_cli_use_the_canonical_result(self) -> None:
        result = simulate_constant_burn(
            example_definition(BurnTarget(BurnTargetKind.DURATION, 1.0)),
            example_operating_point(),
        )
        report = build_burn_html_report(result)
        self.assertIn("<!doctype html>", report)
        self.assertIn(result.result_hash, report)
        self.assertIn("Propulsion-only preliminary result", report)

        example = Path(__file__).resolve().parents[1] / "examples" / "constant_burn_input.json"
        with TemporaryDirectory() as directory:
            result_path = Path(directory) / "result.json"
            report_path = Path(directory) / "report.html"
            sidera_path = Path(directory) / "sidera.json"
            with redirect_stdout(StringIO()):
                status = burn_cli_main(
                    [
                        "simulate",
                        str(example),
                        "--output",
                        str(result_path),
                        "--report",
                        str(report_path),
                    ]
                )
                validation_status = burn_cli_main(["validate", str(result_path)])
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
                    ]
                )
            self.assertEqual(status, 0)
            self.assertEqual(validation_status, 0)
            self.assertEqual(export_status, 0)
            persisted = burn_result_from_json(result_path.read_text(encoding="utf-8"))
            self.assertTrue(persisted.summary.target_achieved)
            self.assertIn(persisted.result_hash, report_path.read_text(encoding="utf-8"))
            sidera = json.loads(sidera_path.read_text(encoding="utf-8"))
            self.assertEqual(sidera["manifest"]["source_result_hash"], persisted.result_hash)

    def test_operating_point_rejects_inconsistent_system_isp(self) -> None:
        point = example_operating_point()
        with self.assertRaises(DomainError):
            EngineOperatingPoint(
                name=point.name,
                commanded_throttle=point.commanded_throttle,
                realized_throttle=point.realized_throttle,
                active_engine_count=point.active_engine_count,
                ideal_thrust_n=point.ideal_thrust_n,
                delivered_thrust_n=point.delivered_thrust_n,
                chamber_specific_impulse_s=point.chamber_specific_impulse_s,
                system_specific_impulse_s=point.system_specific_impulse_s + 1.0,
                chamber_mass_flow_kg_s=point.chamber_mass_flow_kg_s,
                total_tank_flow_kg_s=point.total_tank_flow_kg_s,
                oxidizer_mass_flow_kg_s=point.oxidizer_mass_flow_kg_s,
                fuel_mass_flow_kg_s=point.fuel_mass_flow_kg_s,
                mixture_ratio=point.mixture_ratio,
                chamber_pressure_pa=point.chamber_pressure_pa,
                ambient_pressure_pa=point.ambient_pressure_pa,
                streams=point.streams,
            )

    def test_delta_v_target_equation_has_expected_mass_solution(self) -> None:
        point = example_operating_point()
        requested = 120.0
        result = simulate_constant_burn(
            example_definition(BurnTarget(BurnTargetKind.IDEAL_DELTA_V, requested)), point
        )
        expected_final = 1_000.0 * exp(-requested * 51.0 / 190_000.0)
        self.assertAlmostEqual(result.summary.final_mass_kg, expected_final, places=10)


if __name__ == "__main__":
    unittest.main()
