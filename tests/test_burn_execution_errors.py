"""Navigation-facing propulsion uncertainty and ambient-pressure thrust tests."""

import random
import unittest
from math import exp, log, sqrt

from test_burn_multivariate_calibration import campaign

from rocket_propulsion.core.errors import DomainError, InputError
from rocket_propulsion.propulsion.burns import (
    CalibrationScope,
    FirstOrderGaussMarkov,
    GatesExecutionErrorModel,
    PerformanceAxis,
    PerformanceScalePoint,
    PropagatorPropulsionBridge,
    PropulsionParameterPrior,
    RectilinearPerformanceSurface,
    TabulatedBurnArtifact,
    TabulatedBurnSegment,
    calibrate_hot_fire_multivariate,
    delta_v_log_sensitivities,
)


def calibration():
    return calibrate_hot_fire_multivariate(
        campaign(),
        source="qualification campaign Q2",
        rationale="paired thrust/Isp/time residuals",
    )


class ParameterPriorTests(unittest.TestCase):
    def test_from_standard_deviations_builds_correlated_covariance(self) -> None:
        prior = PropulsionParameterPrior.from_standard_deviations(
            {"thrust_scale": 0.02, "mass_flow_scale": 0.01, "ignition_time_bias_s": 0.05},
            correlation={("thrust_scale", "mass_flow_scale"): 0.5},
            source="test",
            rationale="unit test",
        )
        self.assertEqual(prior.mean, (1.0, 1.0, 0.0))
        self.assertAlmostEqual(prior.covariance[0][1], 0.5 * 0.02 * 0.01)
        self.assertAlmostEqual(prior.covariance[2][2], 0.0025)
        self.assertEqual(prior.variance_of("unknown"), 0.0)
        for kwargs in (
            {"correlation": {("thrust_scale", "thrust_scale"): 0.5}},
            {"correlation": {("thrust_scale", "mass_flow_scale"): 1.5}},
            {"mean": {"other": 1.0}},
        ):
            with self.assertRaises(DomainError):
                PropulsionParameterPrior.from_standard_deviations(
                    {"thrust_scale": 0.02, "mass_flow_scale": 0.01},
                    source="test",
                    rationale="unit test",
                    **kwargs,  # type: ignore[arg-type]
                )
        with self.assertRaises(DomainError):
            PropulsionParameterPrior.from_standard_deviations(
                {"custom": 0.1}, source="test", rationale="no default mean"
            )

    def test_invalid_covariances_fail_closed(self) -> None:
        with self.assertRaises(DomainError):
            PropulsionParameterPrior(("a", "b"), (1.0, 1.0), ((1.0, 2.0), (2.0, 1.0)), "s", "r")
        with self.assertRaises(DomainError):
            PropulsionParameterPrior(("a", "b"), (1.0, 1.0), ((1.0, 0.1), (0.0, 1.0)), "s", "r")
        with self.assertRaises(DomainError):
            PropulsionParameterPrior(("a",), (1.0, 2.0), ((1.0,),), "s", "r")

    def test_hash_round_trip_and_bridge_check(self) -> None:
        prior = PropulsionParameterPrior.from_standard_deviations(
            {"thrust_scale": 0.4, "ignition_time_bias_s": 0.01}, source="s", rationale="r"
        )
        self.assertEqual(PropulsionParameterPrior.from_dict(prior.to_dict()), prior)
        tampered = prior.to_dict()
        tampered["mean"] = [1.1, 0.0]
        with self.assertRaises(InputError):
            PropulsionParameterPrior.from_dict(tampered)
        bridge = PropagatorPropulsionBridge(steady_artifact())
        warnings = prior.check_against(bridge.parameters)
        self.assertEqual(len(warnings), 1)
        self.assertIn("thrust_scale", warnings[0])
        unknown = PropulsionParameterPrior.from_standard_deviations(
            {"custom": 0.1}, mean={"custom": 0.0}, source="s", rationale="r"
        )
        with self.assertRaises(DomainError):
            unknown.check_against(bridge.parameters)

    def test_calibration_mapping_matches_delta_method_and_monte_carlo(self) -> None:
        result = calibration()
        prior = PropulsionParameterPrior.from_multivariate_calibration(
            result, ignition_time_bias_sigma_s=0.02
        )
        self.assertEqual(
            prior.parameter_names, ("thrust_scale", "mass_flow_scale", "ignition_time_bias_s")
        )
        self.assertEqual(prior.mean, (1.0, 1.0, 0.0))
        self.assertIn(result.calibration_hash, prior.source)
        thrust, isp = result.marginals[0], result.marginals[1]
        t_mean, i_mean = thrust.consensus_scale, isp.consensus_scale
        s_t = sqrt(thrust.population_process_variance)
        s_i = sqrt(isp.population_process_variance)
        rho = result.population_correlation.effective_matrix[0][1]
        rel_t, rel_i = s_t / t_mean, s_i / i_mean
        self.assertAlmostEqual(prior.covariance[0][0], rel_t**2, places=14)
        self.assertAlmostEqual(
            prior.covariance[1][1], rel_t**2 + rel_i**2 - 2 * rho * rel_t * rel_i, places=14
        )
        self.assertAlmostEqual(prior.covariance[0][1], rel_t**2 - rho * rel_t * rel_i, places=14)
        self.assertAlmostEqual(prior.covariance[2][2], 0.0004)

        # Nonlinear Monte Carlo of T and T/I agrees with the linearization.
        generator = random.Random(7)
        samples = []
        for _ in range(40_000):
            z1, z2 = generator.gauss(0, 1), generator.gauss(0, 1)
            t = t_mean + s_t * z1
            i = i_mean + s_i * (rho * z1 + sqrt(1 - rho * rho) * z2)
            samples.append((t / t_mean, (t / i) / (t_mean / i_mean)))
        for a in range(2):
            for b in range(2):
                mean_a = sum(item[a] for item in samples) / len(samples)
                mean_b = sum(item[b] for item in samples) / len(samples)
                empirical = sum(
                    (item[a] - mean_a) * (item[b] - mean_b) for item in samples
                ) / (len(samples) - 1)
                self.assertAlmostEqual(
                    empirical / prior.covariance[a][b], 1.0, delta=0.05
                )

        absolute = PropulsionParameterPrior.from_multivariate_calibration(
            result, artifact_includes_consensus=False
        )
        self.assertAlmostEqual(absolute.mean[0], t_mean)
        self.assertAlmostEqual(absolute.mean[1], t_mean / i_mean)
        same_engine = PropulsionParameterPrior.from_multivariate_calibration(
            result, scope=CalibrationScope.SAME_ENGINE
        )
        self.assertLess(same_engine.covariance[0][0], prior.covariance[0][0])


class GatesModelTests(unittest.TestCase):
    model = GatesExecutionErrorModel(0.01, 0.02, 0.005, 0.01)

    def test_covariance_structure(self) -> None:
        delta_v = (3.0, 0.0, 4.0)
        covariance = self.model.covariance(delta_v)
        unit = (0.6, 0.0, 0.8)
        along = sum(unit[i] * covariance[i][j] * unit[j] for i in range(3) for j in range(3))
        self.assertAlmostEqual(along, 0.01**2 + 25 * 0.02**2)
        self.assertAlmostEqual(covariance[1][1], 0.005**2 + 25 * 0.01**2)
        trace = sum(covariance[i][i] for i in range(3))
        self.assertAlmostEqual(trace, along + 2 * (0.005**2 + 25 * 0.01**2))
        with self.assertRaises(DomainError):
            self.model.covariance((0.0, 0.0, 0.0))
        with self.assertRaises(DomainError):
            GatesExecutionErrorModel(-0.1, 0.0, 0.0, 0.0)

    def test_samples_reproduce_covariance(self) -> None:
        delta_v = (1.0, -2.0, 2.0)
        expected = self.model.covariance(delta_v)
        generator = random.Random(11)
        draws = [self.model.sample(delta_v, generator) for _ in range(40_000)]
        for i in range(3):
            for j in range(3):
                empirical = sum(
                    (d[i] - delta_v[i]) * (d[j] - delta_v[j]) for d in draws
                ) / len(draws)
                scale = sqrt(expected[i][i] * expected[j][j])
                self.assertAlmostEqual(empirical / scale, expected[i][j] / scale, delta=0.03)

    def test_delta_v_sensitivities_match_finite_differences(self) -> None:
        mu = 0.3

        def log_dv(s_t: float, s_m: float) -> float:
            return log(s_t / s_m * -log(1 - s_m * mu))

        g_t, g_m = delta_v_log_sensitivities(mu)
        h = 1e-6
        self.assertAlmostEqual(g_t, (log_dv(1 + h, 1) - log_dv(1 - h, 1)) / (2 * h), places=8)
        self.assertAlmostEqual(g_m, (log_dv(1, 1 + h) - log_dv(1, 1 - h)) / (2 * h), places=8)
        self.assertAlmostEqual(delta_v_log_sensitivities(1e-4)[1], 0.5e-4, places=7)
        with self.assertRaises(DomainError):
            delta_v_log_sensitivities(1.0)

    def test_proportional_term_from_prior(self) -> None:
        prior = PropulsionParameterPrior.from_standard_deviations(
            {"thrust_scale": 0.02, "mass_flow_scale": 0.01},
            correlation={("thrust_scale", "mass_flow_scale"): 0.3},
            source="s",
            rationale="r",
        )
        mu = 0.2
        g_t, g_m = delta_v_log_sensitivities(mu)
        model = GatesExecutionErrorModel.from_parameter_prior(
            prior,
            consumed_mass_fraction=mu,
            fixed_magnitude_sigma_m_s=0.003,
            fixed_pointing_sigma_m_s=0.002,
            proportional_pointing_sigma_rad=0.004,
        )
        expected = sqrt(
            g_t**2 * 0.02**2 + g_m**2 * 0.01**2 + 2 * g_t * g_m * 0.3 * 0.02 * 0.01
        )
        self.assertAlmostEqual(model.proportional_magnitude_sigma, expected)
        self.assertEqual(model.fixed_magnitude_sigma_m_s, 0.003)
        with_timing = PropulsionParameterPrior.from_standard_deviations(
            {"thrust_scale": 0.02, "mass_flow_scale": 0.01, "ignition_time_bias_s": 0.1},
            correlation={("thrust_scale", "mass_flow_scale"): 0.3},
            source="s",
            rationale="timing bias has zero mean and must not break the mapping",
        )
        self.assertAlmostEqual(
            GatesExecutionErrorModel.from_parameter_prior(
                with_timing,
                consumed_mass_fraction=mu,
                fixed_magnitude_sigma_m_s=0.003,
                fixed_pointing_sigma_m_s=0.002,
                proportional_pointing_sigma_rad=0.004,
            ).proportional_magnitude_sigma,
            expected,
        )
        timing_only = PropulsionParameterPrior.from_standard_deviations(
            {"ignition_time_bias_s": 0.1}, source="s", rationale="r"
        )
        self.assertEqual(
            GatesExecutionErrorModel.from_parameter_prior(
                timing_only,
                consumed_mass_fraction=mu,
                fixed_magnitude_sigma_m_s=0.0,
                fixed_pointing_sigma_m_s=0.0,
                proportional_pointing_sigma_rad=0.0,
            ).proportional_magnitude_sigma,
            0.0,
        )


class GaussMarkovTests(unittest.TestCase):
    process = FirstOrderGaussMarkov("thrust_scale", 0.01, 600.0)

    def test_discrete_moments_are_exact_and_stationary(self) -> None:
        phi = self.process.transition(300.0)
        self.assertAlmostEqual(phi, exp(-0.5))
        self.assertAlmostEqual(self.process.discrete_noise_variance(300.0), 1e-4 * (1 - phi**2))
        mean, variance = self.process.propagate(0.02, 1e-4, 300.0)
        self.assertAlmostEqual(mean, 0.02 * phi)
        self.assertAlmostEqual(variance, 1e-4)
        self.assertAlmostEqual(self.process.power_spectral_density, 2e-4 / 600.0)
        # Two half steps equal one full step.
        half = self.process.propagate(*self.process.propagate(0.02, 0.0, 150.0), 150.0)
        full = self.process.propagate(0.02, 0.0, 300.0)
        self.assertAlmostEqual(half[0], full[0])
        self.assertAlmostEqual(half[1], full[1])

    def test_sample_paths_have_the_declared_autocorrelation(self) -> None:
        generator = random.Random(3)
        lag = 300.0
        pairs = []
        for _ in range(20_000):
            first, second = self.process.sample_path((0.0, lag), generator)
            pairs.append((first, second))
        variance = sum(a * a for a, _ in pairs) / len(pairs)
        covariance = sum(a * b for a, b in pairs) / len(pairs)
        self.assertAlmostEqual(variance / 1e-4, 1.0, delta=0.04)
        self.assertAlmostEqual(covariance / variance, exp(-lag / 600.0), delta=0.03)
        with self.assertRaises(DomainError):
            self.process.sample_path((1.0, 1.0), generator)
        with self.assertRaises(DomainError):
            FirstOrderGaussMarkov("x", 0.1, 0.0)


def steady_artifact() -> TabulatedBurnArtifact:
    return TabulatedBurnArtifact(
        segments=(
            TabulatedBurnSegment(
                0.0, 10.0, 1000.0, 1000.0, 1000.0, 1000.0, 0.3, 0.3,
                phase="steady",
            ),
        ),
        source_id="ambient-test",
        source_sha256="9" * 64,
        reference_ids=("NASA-20040000363",),
    )


class AmbientPressureTests(unittest.TestCase):
    area = 0.004

    def bridge(self, **kwargs: object) -> PropagatorPropulsionBridge:
        return PropagatorPropulsionBridge(
            steady_artifact(), nozzle_exit_area_m2=self.area, **kwargs  # type: ignore[arg-type]
        )

    def evaluate(self, bridge: PropagatorPropulsionBridge, **kwargs: object):
        arguments = {
            "relative_time_s": 5.0,
            "vehicle_mass_kg": 50.0,
            "direction": (0.0, 0.0, 2.0),
        }
        arguments.update(kwargs)
        return bridge.evaluate(**arguments)  # type: ignore[arg-type]

    def test_pressure_term_mass_flow_and_partial(self) -> None:
        bridge = self.bridge()
        sea_level = self.evaluate(bridge, ambient_pressure_pa=101_325.0)
        self.assertAlmostEqual(sea_level.axial_thrust_n, 1000.0 - 101_325.0 * self.area)
        self.assertAlmostEqual(sea_level.pressure_thrust_loss_n, 101_325.0 * self.area)
        self.assertAlmostEqual(sea_level.net_thrust_margin_n, sea_level.axial_thrust_n)
        self.assertAlmostEqual(sea_level.mass_derivative_kg_s, -0.3)
        vacuum = self.evaluate(bridge, ambient_pressure_pa=0.0)
        self.assertAlmostEqual(vacuum.axial_thrust_n, 1000.0)
        step = 10.0
        upper = self.evaluate(bridge, ambient_pressure_pa=50_000.0 + step).acceleration_m_s2
        lower = self.evaluate(bridge, ambient_pressure_pa=50_000.0 - step).acceleration_m_s2
        partial = self.evaluate(
            bridge, ambient_pressure_pa=50_000.0
        ).partials.acceleration_wrt_ambient_pressure
        for index in range(3):
            self.assertAlmostEqual(
                partial[index], (upper[index] - lower[index]) / (2 * step), places=12
            )
        self.assertAlmostEqual(partial[2], -self.area / 50.0)
        scaled = self.evaluate(
            bridge, ambient_pressure_pa=101_325.0, parameter_overrides={"thrust_scale": 1.1}
        )
        self.assertAlmostEqual(scaled.axial_thrust_n, 1100.0 - 101_325.0 * self.area)

    def test_net_thrust_is_not_clamped_and_margin_goes_negative(self) -> None:
        value = self.evaluate(self.bridge(), ambient_pressure_pa=300_000.0)
        self.assertLess(value.net_thrust_margin_n, 0.0)
        self.assertLess(value.acceleration_m_s2[2], 0.0)

    def test_contract_refusals_and_identity(self) -> None:
        bridge = self.bridge()
        with self.assertRaises(DomainError):
            self.evaluate(bridge)
        self.assertFalse(self.evaluate(bridge, relative_time_s=20.0).active)
        with self.assertRaises(DomainError):
            self.evaluate(PropagatorPropulsionBridge(steady_artifact()), ambient_pressure_pa=1.0)
        with self.assertRaises(DomainError):
            self.evaluate(bridge, ambient_pressure_pa=-1.0)
        with self.assertRaises(DomainError):
            PropagatorPropulsionBridge(steady_artifact(), nozzle_exit_area_m2=0.0)
        surface = RectilinearPerformanceSurface(
            axes=(PerformanceAxis("ambient_pressure_pa", "Pa", (0.0, 1.0e5)),),
            points=(
                PerformanceScalePoint((0.0,), 1.0, 1.0),
                PerformanceScalePoint((1.0e5,), 0.9, 1.0),
            ),
            source_id="double-count",
            source_sha256="8" * 64,
        )
        with self.assertRaises(DomainError):
            PropagatorPropulsionBridge(
                steady_artifact(), performance_surface=surface, nozzle_exit_area_m2=0.01
            )
        vacuum_only = PropagatorPropulsionBridge(steady_artifact())
        self.assertNotIn("nozzle_exit_area_m2", vacuum_only.manifest())
        self.assertNotEqual(vacuum_only.model_hash, bridge.model_hash)
        self.assertEqual(bridge.manifest()["nozzle_exit_area_m2"], self.area)


if __name__ == "__main__":
    unittest.main()
