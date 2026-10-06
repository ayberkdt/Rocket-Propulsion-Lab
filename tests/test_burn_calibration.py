"""Reference and domain tests for hot-fire variance-component calibration."""

from __future__ import annotations

import math
import unittest

from rocket_propulsion.core.errors import DomainError
from rocket_propulsion.propulsion.burns import (
    CalibrationScope,
    DistributionKind,
    HotFireScaleMeasurement,
    UncertaintyClass,
    calibrate_hot_fire_scale,
)


def measurement(
    engine: str,
    run: int,
    value: float,
    *,
    standard_uncertainty: float = 0.0,
    condition: str = "rated",
) -> HotFireScaleMeasurement:
    return HotFireScaleMeasurement(
        engine_id=engine,
        run_id=f"run-{run}",
        condition_id=condition,
        value=value,
        standard_uncertainty=standard_uncertainty,
        source="calibrated hot-fire fixture",
    )


class HotFireCalibrationTests(unittest.TestCase):
    def test_nist_balanced_variance_component_reference_case(self) -> None:
        # NIST/SEMATECH e-Handbook section 7.4.4: five batches, three
        # observations each, error variance 1.8 and between variance 11.71.
        raw = {
            "B1": (74.0, 76.0, 75.0),
            "B2": (68.0, 71.0, 72.0),
            "B3": (75.0, 77.0, 77.0),
            "B4": (72.0, 74.0, 73.0),
            "B5": (79.0, 81.0, 79.0),
        }
        scale = 75.0
        observations = tuple(
            measurement(engine, index, value / scale)
            for engine, values in raw.items()
            for index, value in enumerate(values)
        )
        result = calibrate_hot_fire_scale(
            "thrust_scale",
            observations,
            source="NIST/SEMATECH section 7.4.4 reference table",
            rationale="reproduce the published one-way random-effects example",
        )
        self.assertAlmostEqual(
            result.observed_within_engine_variance * scale**2,
            1.8,
            places=10,
        )
        self.assertAlmostEqual(
            result.between_engine_process_variance * scale**2,
            (36.93333333333333 - 1.8) / 3.0,
            places=7,
        )
        self.assertAlmostEqual(result.consensus_scale * scale, 1123.0 / 15.0)
        self.assertTrue(result.evidence.converged)
        self.assertIn("NIST-SEMATECH-7.4.4", result.reference_ids)

    def test_unbalanced_measurements_produce_population_and_same_engine_scopes(self) -> None:
        offsets = (-0.03, -0.02, -0.01, 0.01, 0.02, 0.03)
        observations: list[HotFireScaleMeasurement] = []
        for engine_index, offset in enumerate(offsets):
            deviations = (-0.004, 0.0, 0.004) if engine_index % 2 == 0 else (-0.004, 0.004)
            observations.extend(
                measurement(
                    f"E{engine_index}",
                    run_index,
                    1.01 + offset + deviation,
                    standard_uncertainty=0.001,
                )
                for run_index, deviation in enumerate(deviations)
            )
        result = calibrate_hot_fire_scale(
            "specific_impulse_scale",
            tuple(observations),
            source="qualification campaign fixture",
            rationale="separate engine-to-engine and repeat-burn variability",
        )
        self.assertEqual(result.engine_count, 6)
        self.assertEqual(result.observation_count, 15)
        self.assertAlmostEqual(result.consensus_scale, 1.01, delta=5e-5)
        self.assertGreater(
            result.between_engine_process_variance,
            result.within_engine_process_variance,
        )
        self.assertAlmostEqual(
            result.within_engine_fraction + result.between_engine_fraction,
            1.0,
        )
        population = result.residual_uncertain_parameter(
            scope=CalibrationScope.POPULATION
        )
        same_engine = result.residual_uncertain_parameter(
            scope=CalibrationScope.SAME_ENGINE
        )
        self.assertEqual(population.nominal, 1.0)
        self.assertEqual(population.distribution, DistributionKind.TRUNCATED_NORMAL)
        self.assertEqual(population.uncertainty_class, UncertaintyClass.ALEATORY)
        self.assertGreater(
            population.standard_deviation,
            same_engine.standard_deviation,
        )
        self.assertIn(result.calibration_hash, population.source)

    def test_measurement_uncertainty_is_not_mislabelled_as_process_variance(self) -> None:
        standard_uncertainty = math.sqrt(0.0002)
        observations = tuple(
            measurement(
                engine,
                run,
                value,
                standard_uncertainty=standard_uncertainty,
            )
            for engine in ("A", "B")
            for run, value in enumerate((0.99, 1.01))
        )
        result = calibrate_hot_fire_scale(
            "thrust_scale",
            observations,
            source="measurement-dominated fixture",
            rationale="verify Type B contribution is removed from process scatter",
        )
        self.assertAlmostEqual(result.observed_within_engine_variance, 0.0002)
        self.assertAlmostEqual(
            result.measurement_within_variance_contribution,
            0.0002,
        )
        self.assertEqual(result.within_engine_process_variance, 0.0)
        self.assertEqual(result.between_engine_process_variance, 0.0)
        with self.assertRaisesRegex(DomainError, "zero"):
            result.residual_uncertain_parameter()

    def test_mixed_conditions_duplicate_runs_and_provenance_changes_fail_or_hash(self) -> None:
        observations = (
            measurement("A", 0, 0.99),
            measurement("A", 1, 1.01),
            measurement("B", 0, 1.02),
        )
        first = calibrate_hot_fire_scale(
            "time_scale",
            observations,
            source="campaign revision A",
            rationale="trace calibration provenance",
        )
        second = calibrate_hot_fire_scale(
            "time_scale",
            observations,
            source="campaign revision B",
            rationale="trace calibration provenance",
        )
        self.assertNotEqual(first.calibration_hash, second.calibration_hash)

        mixed = observations[:-1] + (
            measurement("B", 0, 1.02, condition="off-design"),
        )
        with self.assertRaisesRegex(DomainError, "different operating conditions"):
            calibrate_hot_fire_scale(
                "time_scale",
                mixed,
                source="invalid mixed fixture",
                rationale="prove pooling refusal",
            )
        duplicate = observations + (measurement("A", 0, 1.00),)
        with self.assertRaisesRegex(DomainError, "identities"):
            calibrate_hot_fire_scale(
                "time_scale",
                duplicate,
                source="invalid duplicate fixture",
                rationale="prove duplicate refusal",
            )


if __name__ == "__main__":
    unittest.main()

