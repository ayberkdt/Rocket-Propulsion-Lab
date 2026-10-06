"""Paired hot-fire covariance, persistence, and ensemble-contract tests."""

from __future__ import annotations

import copy
import unittest
from random import Random

from rocket_propulsion.core.errors import DomainError, InputError
from rocket_propulsion.propulsion.burns import (
    CalibrationScope,
    HotFireMultivariateCalibration,
    HotFireVectorMeasurement,
    calibrate_hot_fire_multivariate,
    sample_uncertain_parameters,
)

PARAMETERS = ("thrust_scale", "specific_impulse_scale", "time_scale")


def campaign() -> tuple[HotFireVectorMeasurement, ...]:
    engine_effects = (-0.030, -0.018, -0.006, 0.006, 0.018, 0.030)
    engine_secondary = (0.002, -0.001, 0.0015, -0.002, 0.0005, -0.001)
    run_patterns = (
        ((-0.004, 0.002, -0.001), (0.001, -0.003, 0.003), (0.003, 0.001, -0.002)),
        ((-0.003, 0.003, -0.002), (0.002, -0.002, 0.001), (0.001, -0.001, 0.001)),
    )
    covariance = (
        (1.0e-6, 0.20e-6, -0.10e-6),
        (0.20e-6, 0.64e-6, 0.08e-6),
        (-0.10e-6, 0.08e-6, 0.49e-6),
    )
    observations: list[HotFireVectorMeasurement] = []
    for engine_index, engine_effect in enumerate(engine_effects):
        for run_index, run_effect in enumerate(run_patterns[engine_index % 2]):
            observations.append(
                HotFireVectorMeasurement(
                    engine_id=f"E{engine_index + 1}",
                    run_id=f"R{run_index + 1}",
                    condition_id="rated-100pct",
                    parameter_names=PARAMETERS,
                    values=(
                        1.010 + engine_effect + run_effect[0],
                        0.997
                        + 0.55 * engine_effect
                        + engine_secondary[engine_index]
                        + run_effect[1],
                        1.004
                        - 0.30 * engine_effect
                        + 0.5 * engine_secondary[engine_index]
                        + run_effect[2],
                    ),
                    measurement_covariance=covariance,
                    source="paired qualification campaign fixture",
                )
            )
    return tuple(observations)


class HotFireMultivariateCalibrationTests(unittest.TestCase):
    def setUp(self) -> None:
        self.result = calibrate_hot_fire_multivariate(
            campaign(),
            source="qualification campaign Q2",
            rationale="retain paired thrust, Isp, and duration residual dependence",
        )

    def test_paired_fit_builds_ensemble_ready_population_model(self) -> None:
        self.assertEqual(self.result.parameter_names, PARAMETERS)
        self.assertEqual(len(self.result.marginals), 3)
        self.assertEqual(self.result.marginals[0].engine_count, 6)
        self.assertEqual(self.result.marginals[0].observation_count, 18)
        self.assertGreater(
            self.result.population_correlation.effective_matrix[1][0],
            0.5,
        )
        self.assertLess(
            self.result.population_correlation.effective_matrix[2][0],
            -0.4,
        )
        parameters, correlation = self.result.residual_model()
        self.assertEqual(tuple(item.name for item in parameters), PARAMETERS)
        self.assertEqual(correlation.parameter_names, PARAMETERS)
        self.assertIn(self.result.calibration_hash, correlation.source)
        draws = sample_uncertain_parameters(
            parameters,
            sample_count=128,
            seed=42,
            correlation=correlation,
        )
        self.assertEqual(len(draws), 128)
        for draw in draws:
            values = draw.as_dict()
            for parameter in parameters:
                self.assertGreaterEqual(values[parameter.name], parameter.lower)
                self.assertLessEqual(values[parameter.name], parameter.upper)

    def test_covariance_diagonals_reproduce_scalar_variance_components(self) -> None:
        for index, marginal in enumerate(self.result.marginals):
            self.assertAlmostEqual(
                self.result.within_engine_process_covariance[index][index],
                marginal.within_engine_process_variance,
            )
            self.assertAlmostEqual(
                self.result.between_engine_process_covariance[index][index],
                marginal.between_engine_process_variance,
            )
            self.assertAlmostEqual(
                self.result.population_process_covariance[index][index],
                marginal.population_process_variance,
            )
        if self.result.same_engine_correlation is not None:
            parameters, correlation = self.result.residual_model(scope=CalibrationScope.SAME_ENGINE)
            self.assertEqual(len(parameters), 3)
            self.assertEqual(
                correlation.matrix, self.result.same_engine_correlation.effective_matrix
            )

    def test_declared_cross_measurement_covariance_is_removed(self) -> None:
        covariance = ((0.0002, 0.0004), (0.0004, 0.0008))
        observations = tuple(
            HotFireVectorMeasurement(
                engine_id=engine,
                run_id=f"R{run}",
                condition_id="measurement-removal",
                parameter_names=("thrust_scale", "specific_impulse_scale"),
                values=values,
                measurement_covariance=covariance,
                source="analytic covariance-removal fixture",
            )
            for engine, pairs in (
                ("E1", ((0.98, 0.96), (1.00, 1.00))),
                ("E2", ((1.00, 1.00), (1.02, 1.04))),
            )
            for run, values in enumerate(pairs, start=1)
        )
        result = calibrate_hot_fire_multivariate(
            observations,
            source="analytic paired fixture",
            rationale="prove declared measurement cross-covariance subtraction",
        )
        self.assertAlmostEqual(result.within_engine_process_covariance[0][0], 0.0)
        self.assertAlmostEqual(result.within_engine_process_covariance[1][1], 0.0)
        self.assertAlmostEqual(result.within_engine_process_covariance[1][0], 0.0)
        self.assertIsNone(result.same_engine_correlation)

    def test_nonphysical_raw_pairwise_estimate_is_regularized_with_evidence(self) -> None:
        random = Random(0)
        covariance = tuple(tuple(0.0001 for _ in range(3)) for _ in range(3))
        observations = tuple(
            HotFireVectorMeasurement(
                engine_id=f"E{engine}",
                run_id=f"R{run}",
                condition_id="regularization-fixture",
                parameter_names=PARAMETERS,
                values=tuple(1.0 + random.uniform(-0.03, 0.03) for _ in PARAMETERS),
                measurement_covariance=covariance,
                source="deterministic regularization fixture",
            )
            for engine in range(6)
            for run in range(3)
        )
        result = calibrate_hot_fire_multivariate(
            observations,
            source="deterministic regularization fixture",
            rationale="exercise explicit positive-definite correlation repair",
        )
        evidence = result.population_correlation
        self.assertTrue(evidence.regularized)
        self.assertGreater(evidence.coefficient_clip_count, 0)
        self.assertGreater(evidence.identity_shrinkage, 0.0)
        _parameters, correlation = result.residual_model()
        self.assertEqual(correlation.matrix, evidence.effective_matrix)

    def test_versioned_round_trip_recomputes_and_refuses_semantic_tampering(self) -> None:
        text = self.result.to_json()
        reopened = HotFireMultivariateCalibration.from_json(text)
        self.assertEqual(reopened, self.result)
        self.assertEqual(reopened.calibration_hash, self.result.calibration_hash)

        changed_hash = copy.deepcopy(self.result.to_dict())
        changed_hash["calibration_hash"] = "0" * 64
        with self.assertRaisesRegex(InputError, "hash"):
            HotFireMultivariateCalibration.from_dict(changed_hash)

        changed_output = copy.deepcopy(self.result.to_dict())
        changed_output["population_process_covariance"][0][1] += 0.01
        with self.assertRaisesRegex(InputError, "semantics"):
            HotFireMultivariateCalibration.from_dict(changed_output)

    def test_invalid_measurement_covariance_and_parameter_order_fail_closed(self) -> None:
        with self.assertRaisesRegex(DomainError, "positive semidefinite"):
            HotFireVectorMeasurement(
                engine_id="E1",
                run_id="R1",
                condition_id="rated",
                parameter_names=("thrust_scale", "specific_impulse_scale"),
                values=(1.0, 1.0),
                measurement_covariance=((1.0, 2.0), (2.0, 1.0)),
                source="invalid fixture",
            )

        observations = list(campaign())
        item = observations[-1]
        observations[-1] = HotFireVectorMeasurement(
            engine_id=item.engine_id,
            run_id=item.run_id,
            condition_id=item.condition_id,
            parameter_names=tuple(reversed(item.parameter_names)),
            values=tuple(reversed(item.values)),
            measurement_covariance=tuple(
                tuple(item.measurement_covariance[row][column] for row in reversed(range(3)))
                for column in reversed(range(3))
            ),
            source=item.source,
        )
        with self.assertRaisesRegex(DomainError, "identical parameter order"):
            calibrate_hot_fire_multivariate(
                tuple(observations),
                source="invalid fixture",
                rationale="prove order refusal",
            )


if __name__ == "__main__":
    unittest.main()
