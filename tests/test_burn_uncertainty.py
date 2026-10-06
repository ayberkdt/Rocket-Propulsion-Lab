"""Deterministic and physical checks for propulsion uncertainty ensembles."""

from __future__ import annotations

import unittest

from rocket_propulsion.core.errors import DomainError
from rocket_propulsion.propulsion.burns import (
    BlowdownTankModel,
    ConstantPerformanceProvider,
    CorrelationModel,
    DistributionKind,
    PressurePerformanceLaw,
    UncertainParameter,
    UncertaintyClass,
    sample_uncertain_parameters,
    simulate_blowdown_uncertainty,
)


def parameter(
    name: str,
    nominal: float,
    lower: float,
    upper: float,
    distribution: DistributionKind = DistributionKind.UNIFORM,
    *,
    sigma: float | None = None,
    uncertainty_class: UncertaintyClass = UncertaintyClass.EPISTEMIC,
) -> UncertainParameter:
    """Build a sourced test uncertainty input."""

    return UncertainParameter(
        name=name,
        nominal=nominal,
        lower=lower,
        upper=upper,
        distribution=distribution,
        uncertainty_class=uncertainty_class,
        source="qualification-test synthetic fixture",
        rationale="exercise bounded and reproducible propagation",
        standard_deviation=sigma,
    )


class UncertaintySamplingTests(unittest.TestCase):
    def test_marginals_are_bounded_reproducible_and_stratified(self) -> None:
        parameters = (
            parameter("uniform", 1.0, 0.5, 1.5),
            parameter("triangle", 2.0, 1.0, 4.0, DistributionKind.TRIANGULAR),
            parameter(
                "normal",
                10.0,
                8.0,
                12.0,
                DistributionKind.TRUNCATED_NORMAL,
                sigma=0.5,
            ),
            parameter("log", 10.0, 1.0, 100.0, DistributionKind.LOG_UNIFORM),
        )
        first = sample_uncertain_parameters(parameters, sample_count=32, seed=17)
        second = sample_uncertain_parameters(parameters, sample_count=32, seed=17)
        self.assertEqual(first, second)
        for draw in first:
            values = draw.as_dict()
            for item in parameters:
                self.assertGreaterEqual(values[item.name], item.lower)
                self.assertLessEqual(values[item.name], item.upper)

    def test_gaussian_copula_applies_declared_positive_correlation(self) -> None:
        parameters = (
            parameter("a", 0.0, -1.0, 1.0),
            parameter("b", 0.0, -1.0, 1.0),
        )
        correlation = CorrelationModel(
            parameter_names=("a", "b"),
            matrix=((1.0, 0.8), (0.8, 1.0)),
            source="joint hot-fire covariance fixture",
            rationale="exercise dependent propulsion inputs",
        )
        draws = sample_uncertain_parameters(
            parameters,
            sample_count=512,
            seed=3,
            correlation=correlation,
        )
        a_values = [draw.as_dict()["a"] for draw in draws]
        b_values = [draw.as_dict()["b"] for draw in draws]
        a_mean = sum(a_values) / len(a_values)
        b_mean = sum(b_values) / len(b_values)
        covariance = sum(
            (a - a_mean) * (b - b_mean)
            for a, b in zip(a_values, b_values, strict=True)
        )
        a_variance = sum((value - a_mean) ** 2 for value in a_values)
        b_variance = sum((value - b_mean) ** 2 for value in b_values)
        empirical = covariance / (a_variance * b_variance) ** 0.5
        self.assertGreater(empirical, 0.7)

        with self.assertRaisesRegex(DomainError, "positive definite"):
            CorrelationModel(
                parameter_names=("a", "b"),
                matrix=((1.0, 1.0), (1.0, 1.0)),
                source="invalid fixture",
                rationale="prove singular matrices fail closed",
            )


class BlowdownUncertaintyTests(unittest.TestCase):
    def setUp(self) -> None:
        self.tank = BlowdownTankModel(
            initial_pressure_pa=2_000_000.0,
            initial_ullage_volume_m3=0.01,
            initial_propellant_mass_kg=10.0,
            propellant_density_kg_m3=1_000.0,
            polytropic_exponent=1.2,
        )
        self.point = ConstantPerformanceProvider.monopropellant(
            delivered_thrust_n=100.0,
            system_specific_impulse_s=220.0,
            chamber_pressure_pa=1_500_000.0,
            source="qualification-test rated point",
        ).operating_point()
        self.law = PressurePerformanceLaw(
            reference_supply_pressure_pa=2_000_000.0,
            minimum_supply_pressure_pa=1_000_000.0,
            maximum_supply_pressure_pa=2_200_000.0,
            thrust_pressure_exponent=1.0,
            specific_impulse_pressure_exponent=0.02,
            chamber_pressure_exponent=1.0,
        )
        self.parameters = (
            parameter(
                "initial_pressure_pa",
                2_000_000.0,
                1_850_000.0,
                2_150_000.0,
                DistributionKind.TRUNCATED_NORMAL,
                sigma=50_000.0,
                uncertainty_class=UncertaintyClass.ALEATORY,
            ),
            parameter(
                "rated_delivered_thrust_n",
                100.0,
                95.0,
                105.0,
                DistributionKind.TRIANGULAR,
            ),
            parameter(
                "rated_system_specific_impulse_s",
                220.0,
                216.0,
                224.0,
                DistributionKind.TRIANGULAR,
            ),
            parameter("polytropic_exponent", 1.2, 1.1, 1.3),
        )

    def test_seeded_ensemble_preserves_every_artifact_closure(self) -> None:
        first = simulate_blowdown_uncertainty(
            self.tank,
            self.point,
            self.law,
            self.parameters,
            sample_count=16,
            seed=91,
            maximum_mass_step_kg=0.5,
            integration_relative_tolerance=1e-6,
        )
        second = simulate_blowdown_uncertainty(
            self.tank,
            self.point,
            self.law,
            self.parameters,
            sample_count=16,
            seed=91,
            maximum_mass_step_kg=0.5,
            integration_relative_tolerance=1e-6,
        )
        self.assertEqual(
            tuple(sample.artifact.artifact_hash for sample in first.samples),
            tuple(sample.artifact.artifact_hash for sample in second.samples),
        )
        self.assertEqual(first.aleatory_parameter_count, 1)
        self.assertEqual(first.epistemic_parameter_count, 3)
        self.assertEqual(len(first.statistics), 5)
        self.assertEqual(len(first.sensitivity), 5 * len(self.parameters))
        self.assertEqual(first.convergence.full_sample_count, 16)
        for sample in first.samples:
            self.assertGreater(sample.duration_s, 0.0)
            self.assertGreater(sample.delivered_impulse_n_s, 0.0)
            self.assertAlmostEqual(
                sample.consumed_propellant_kg,
                sample.artifact.consumed_propellant_kg,
            )
            self.assertIn("NASA-STD-7009B", sample.artifact.reference_ids)

    def test_nominal_mismatch_and_pressure_extrapolation_fail_before_sampling(self) -> None:
        mismatch = (
            parameter("rated_delivered_thrust_n", 99.0, 95.0, 105.0),
        )
        with self.assertRaisesRegex(DomainError, "does not match"):
            simulate_blowdown_uncertainty(
                self.tank, self.point, self.law, mismatch, sample_count=8
            )
        invalid_pressure = (
            parameter("initial_pressure_pa", 2_000_000.0, 1_900_000.0, 2_300_000.0),
        )
        with self.assertRaisesRegex(DomainError, "validity"):
            simulate_blowdown_uncertainty(
                self.tank, self.point, self.law, invalid_pressure, sample_count=8
            )


if __name__ == "__main__":
    unittest.main()
