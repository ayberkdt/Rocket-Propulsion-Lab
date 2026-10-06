"""Isolated tests for future Sidera transient-burn capabilities."""

from __future__ import annotations

import math
import tempfile
import unittest
from pathlib import Path

from rocket_propulsion.core.errors import DomainError, InputError
from rocket_propulsion.propulsion.burns import (
    BurnDefinition,
    BurnTarget,
    BurnTargetKind,
    ConstantPerformanceProvider,
    EngineClusterMember,
    EngineDutyCycleLimits,
    EnginePulse,
    EngineTransientCommand,
    InterpolationPolicy,
    TabulatedBurnArtifact,
    TankDefinition,
    ThrottleSchedule,
    ThrottleSegment,
    build_engine_cluster_profile,
    realize_first_order_schedule,
    simulate_profiled_burn,
    tabulated_artifact_from_csv,
    tabulated_artifact_from_profiled_burn,
    validate_engine_sequence,
)
from rocket_propulsion.propulsion.burns.cli import main as burn_cli_main
from rocket_propulsion.propulsion.burns.models import PropellantRole


class ResponseModelTests(unittest.TestCase):
    def test_first_order_step_has_exact_exposure_and_centroid(self) -> None:
        schedule = ThrottleSchedule((ThrottleSegment(10.0, 1.0, 1.0, "on"),))
        history = realize_first_order_schedule(schedule, time_constant_s=2.0)

        exponential = math.exp(-5.0)
        expected_exposure = 10.0 - 2.0 * (1.0 - exponential)
        expected_moment = 50.0 - (
            4.0 * (1.0 - exponential) - 20.0 * exponential
        )
        self.assertAlmostEqual(history.exposure_s, expected_exposure, places=12)
        self.assertAlmostEqual(
            history.throttle_first_moment_s2, expected_moment, places=12
        )
        self.assertAlmostEqual(
            history.segments[0].realized_end, 1.0 - exponential, places=12
        )

    def test_zero_lag_reproduces_piecewise_linear_command_exactly(self) -> None:
        schedule = ThrottleSchedule(
            (
                ThrottleSegment(2.0, 0.0, 1.0, "start"),
                ThrottleSegment(3.0, 1.0, 1.0, "steady"),
                ThrottleSegment(2.0, 1.0, 0.0, "stop"),
            )
        )
        history = realize_first_order_schedule(schedule, time_constant_s=0.0)
        self.assertAlmostEqual(history.exposure_s, schedule.total_exposure_s)
        self.assertEqual(history.realized_throttle(2.0, side="left"), 1.0)
        self.assertEqual(history.realized_throttle(2.0, side="right"), 1.0)


class EngineSequenceTests(unittest.TestCase):
    def test_restart_cooldown_and_minimum_pulse_are_enforced(self) -> None:
        limits = (
            EngineDutyCycleLimits(
                "main", minimum_on_time_s=2.0, minimum_off_time_s=5.0, maximum_starts=2
            ),
        )
        pulses = (
            EnginePulse("main", 0.0, 3.0, 0.8),
            EnginePulse("main", 8.0, 2.0, 1.0),
        )
        validated = validate_engine_sequence(pulses, limits)
        schedule = validated.schedule_for_engine("main", total_duration_s=12.0)
        self.assertEqual(validated.starts_by_engine, (("main", 2),))
        self.assertAlmostEqual(schedule.total_exposure_s, 4.4)
        self.assertEqual(schedule.total_duration_s, 12.0)

        with self.assertRaisesRegex(DomainError, "cooldown"):
            validate_engine_sequence(
                (pulses[0], EnginePulse("main", 7.9, 2.0)), limits
            )
        with self.assertRaisesRegex(DomainError, "minimum on-time"):
            validate_engine_sequence((EnginePulse("main", 0.0, 1.0),), limits)
        with self.assertRaisesRegex(DomainError, "maximum start count"):
            validate_engine_sequence(
                pulses + (EnginePulse("main", 15.0, 2.0),), limits
            )


class TabulatedArtifactTests(unittest.TestCase):
    TRACE = """time_s,delivered_thrust_n,total_mass_flow_kg_s,stream:solid,phase
0,0,0,0,ignition
1,100,1,1,rise
3,0,0,0,tailoff
"""

    def test_imported_trace_preserves_samples_integrals_and_hash(self) -> None:
        artifact = tabulated_artifact_from_csv(
            self.TRACE, source_id="qualification-static-fire-001"
        )
        self.assertAlmostEqual(artifact.delivered_total_impulse_n_s, 150.0)
        self.assertAlmostEqual(artifact.consumed_propellant_kg, 1.5)
        self.assertAlmostEqual(artifact.thrust_centroid_time_s, 4.0 / 3.0)
        self.assertIn("without smoothing", artifact.warnings[0])
        self.assertNotIn("reference_ids", artifact.to_dict())

        reopened = TabulatedBurnArtifact.from_json(artifact.to_json())
        self.assertEqual(reopened.artifact_hash, artifact.artifact_hash)
        tampered = artifact.to_dict()
        tampered["segments"][0]["end_delivered_thrust_n"] = 101.0
        with self.assertRaisesRegex(InputError, "hash"):
            TabulatedBurnArtifact.from_dict(tampered)

    def test_hold_interpolation_is_explicit(self) -> None:
        artifact = tabulated_artifact_from_csv(
            self.TRACE,
            source_id="solid-hold",
            interpolation=InterpolationPolicy.HOLD,
        )
        self.assertAlmostEqual(artifact.delivered_total_impulse_n_s, 200.0)
        self.assertAlmostEqual(artifact.consumed_propellant_kg, 2.0)

    def test_cli_builds_reopenable_trace_artifact(self) -> None:
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            trace_path = root / "trace.csv"
            output_path = root / "artifact.json"
            trace_path.write_text(self.TRACE, encoding="utf-8")
            exit_code = burn_cli_main(
                (
                    "tabulate-trace",
                    str(trace_path),
                    "--source-id",
                    "cli-trace",
                    "--output",
                    str(output_path),
                )
            )
            self.assertEqual(exit_code, 0)
            artifact = TabulatedBurnArtifact.from_json(
                output_path.read_text(encoding="utf-8")
            )
            self.assertEqual(artifact.source_id, "cli-trace")
            self.assertAlmostEqual(artifact.delivered_total_impulse_n_s, 150.0)

    def test_l1_conversion_retains_a_command_jump_without_smearing(self) -> None:
        provider = ConstantPerformanceProvider.monopropellant(
            delivered_thrust_n=100.0,
            system_specific_impulse_s=200.0,
            chamber_pressure_pa=1.0e6,
        )
        definition = BurnDefinition(
            name="jump profile",
            initial_mass_kg=100.0,
            protected_dry_mass_kg=50.0,
            tanks=(
                TankDefinition(
                    "propellant", "monopropellant", PropellantRole.MONOPROPELLANT, 20.0
                ),
            ),
            target=BurnTarget(BurnTargetKind.DURATION, 3.0),
            maximum_duration_s=3.0,
        )
        schedule = ThrottleSchedule(
            (
                ThrottleSegment(1.0, 0.0, 1.0, "ramp"),
                ThrottleSegment(1.0, 1.0, 1.0, "steady"),
                ThrottleSegment(1.0, 0.0, 0.0, "cutoff"),
            )
        )
        result = simulate_profiled_burn(definition, provider.point, schedule)
        artifact = tabulated_artifact_from_profiled_burn(result)

        self.assertEqual(artifact.state_at(2.0, side="left")["delivered_thrust_n"], 100.0)
        self.assertEqual(artifact.state_at(2.0, side="right")["delivered_thrust_n"], 0.0)
        self.assertAlmostEqual(
            artifact.delivered_total_impulse_n_s,
            result.summary.delivered_total_impulse_n_s,
        )


class EngineClusterTests(unittest.TestCase):
    def test_cluster_sampling_converges_to_exact_engine_response(self) -> None:
        point = ConstantPerformanceProvider.monopropellant(
            delivered_thrust_n=500.0,
            system_specific_impulse_s=220.0,
            chamber_pressure_pa=1.0e6,
        ).point
        schedule = ThrottleSchedule(
            (
                ThrottleSegment(5.0, 1.0, 1.0, "on"),
                ThrottleSegment(1.0, 0.0, 0.0, "shutdown"),
            )
        )
        members = (
            EngineClusterMember("left", point, cant_efficiency=1.0),
            EngineClusterMember("right", point, cant_efficiency=0.9),
            EngineClusterMember("disabled-spare", point, enabled=False),
        )
        commands = (
            EngineTransientCommand(
                "left", schedule, response_time_constant_s=0.5, settling_duration_s=3.0
            ),
            EngineTransientCommand(
                "right",
                schedule,
                response_time_constant_s=0.5,
                ignition_delay_s=1.0,
                settling_duration_s=3.0,
            ),
        )
        result = build_engine_cluster_profile(
            members,
            commands,
            maximum_step_s=1.0,
            relative_tolerance=1e-6,
        )
        evidence = result.evidence
        self.assertLessEqual(
            evidence.impulse_error_n_s,
            evidence.exact_delivered_impulse_n_s * evidence.relative_tolerance,
        )
        self.assertLessEqual(
            evidence.propellant_error_kg,
            evidence.exact_propellant_mass_kg * evidence.relative_tolerance,
        )
        self.assertLess(
            evidence.exact_axial_impulse_n_s,
            evidence.exact_delivered_impulse_n_s,
        )
        self.assertGreater(evidence.refinement_count, 0)
        self.assertEqual(result.artifact.source_id, "engine-cluster-response")
        self.assertNotIn("disabled-spare", dict(result.histories))


if __name__ == "__main__":
    unittest.main()
