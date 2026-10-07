"""Free-space burn metrics: Gauss-Legendre rule, Dv and impulsive-equivalent time."""

import random
import unittest
from math import isclose, log

from test_burn_input_robustness import random_artifact
from test_burn_propagator_integration import STATES, EventAwareRk4

from rocket_propulsion.core.errors import DomainError
from rocket_propulsion.propulsion.burns import (
    PropagatorPropulsionBridge,
    StreamStateBinding,
    TabulatedBurnArtifact,
    TabulatedBurnSegment,
    free_space_burn_metrics,
    gauss_legendre,
)


def steady(thrust: float, flow: float, duration: float) -> TabulatedBurnArtifact:
    return TabulatedBurnArtifact(
        segments=(TabulatedBurnSegment(0.0, duration, thrust, thrust, thrust, thrust,
                                       flow, flow, phase="steady"),),
        source_id="metrics", source_sha256="6" * 64,
    )


class GaussLegendreTests(unittest.TestCase):
    def test_rule_integrates_polynomials_exactly(self) -> None:
        for points in (1, 2, 5, 16, 40):
            nodes, weights = gauss_legendre(points)
            self.assertAlmostEqual(sum(weights), 2.0, places=13)
            for degree in range(2 * points):
                exact = 0.0 if degree % 2 else 2.0 / (degree + 1)
                numeric = sum(w * x**degree for x, w in zip(nodes, weights, strict=True))
                self.assertAlmostEqual(numeric, exact, places=12)
        with self.assertRaises(DomainError):
            gauss_legendre(0)


class FreeSpaceMetricsTests(unittest.TestCase):
    def test_constant_thrust_matches_closed_forms(self) -> None:
        thrust, flow, duration, m0 = 1200.0, 0.4, 900.0, 1000.0
        metrics = free_space_burn_metrics(steady(thrust, flow, duration), m0)
        exhaust = thrust / flow
        mf = m0 - flow * duration
        self.assertAlmostEqual(metrics.final_mass_kg, mf)
        self.assertTrue(isclose(metrics.delta_v_m_s, exhaust * log(m0 / mf), rel_tol=1e-13))
        self.assertTrue(isclose(metrics.effective_exhaust_velocity_m_s, exhaust, rel_tol=1e-13))
        # t* = int t/(m0 - q t) dt / int 1/(m0 - q t) dt
        numerator = -duration / flow - m0 / flow**2 * log(mf / m0)
        denominator = log(m0 / mf) / flow
        self.assertTrue(isclose(metrics.impulsive_equivalent_time_s, numerator / denominator,
                                rel_tol=1e-12))
        self.assertGreater(metrics.impulsive_equivalent_time_s, duration / 2)  # mass falls
        self.assertLess(metrics.delta_v_error_estimate_m_s, 1e-9)

    def test_matches_event_aware_rk4_on_random_profiles(self) -> None:
        generator = random.Random(77)
        for _ in range(10):
            artifact = random_artifact(generator)
            bridge = PropagatorPropulsionBridge(
                artifact,
                stream_bindings=(StreamStateBinding("ox", STATES[0]),
                                 StreamStateBinding("fu", STATES[1])),
            )
            metrics = free_space_burn_metrics(artifact, 2000.0)
            final = EventAwareRk4(bridge, (1.0, 0.0, 0.0)).propagate(
                0.0, artifact.duration_s, [0, 0, 0, 0, 0, 0, 2000.0, 900.0, 300.0], 0.02
            )
            self.assertTrue(isclose(metrics.delta_v_m_s, final[3], rel_tol=1e-9, abs_tol=1e-9))
            self.assertAlmostEqual(metrics.final_mass_kg, final[6], places=9)

    def test_scales_and_refusals(self) -> None:
        base = free_space_burn_metrics(steady(100.0, 0.1, 10.0), 50.0)
        scaled = free_space_burn_metrics(steady(100.0, 0.1, 10.0), 50.0, thrust_scale=1.1)
        self.assertAlmostEqual(scaled.delta_v_m_s, 1.1 * base.delta_v_m_s)
        with self.assertRaises(DomainError):
            free_space_burn_metrics(steady(100.0, 0.1, 10.0), 1.0)
        with self.assertRaises(DomainError):
            free_space_burn_metrics(steady(100.0, 0.1, 10.0), 50.0, mass_flow_scale=0.0)
        with self.assertRaises(DomainError):
            free_space_burn_metrics("not an artifact", 50.0)  # type: ignore[arg-type]


if __name__ == "__main__":
    unittest.main()
