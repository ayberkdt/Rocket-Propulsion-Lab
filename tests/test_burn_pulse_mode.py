"""Pulse-mode thruster: impulse bits, quantization and pulse-train artifacts."""

import random
import unittest

from rocket_propulsion.core.errors import DomainError, InputError
from rocket_propulsion.propulsion.burns import (
    PropagatorPropulsionBridge,
    PropulsionTimeline,
    PulseModeThruster,
    PulsePlan,
    ScheduledBurn,
    free_space_burn_metrics,
)
from rocket_propulsion.propulsion.performance import STANDARD_GRAVITY_M_S2

THRUST = 1.0
FLOW = THRUST / (STANDARD_GRAVITY_M_S2 * 220.0)
ON_TIMES = (0.010, 0.020, 0.050, 0.100, 0.200)


def synthetic_thruster(**changes) -> PulseModeThruster:
    """Synthetic (not measured) 1 N thruster: 6 ms dead time, 2 mN·s tail-off."""

    values = {
        "on_times_s": ON_TIMES,
        "impulse_bits_n_s": tuple(THRUST * (t - 0.006) + 0.002 for t in ON_TIMES),
        "propellant_per_pulse_kg": tuple(FLOW * (t + 0.003) for t in ON_TIMES),
        "steady_thrust_n": THRUST,
        "steady_mass_flow_kg_s": FLOW,
        "minimum_off_time_s": 0.02,
        "source_id": "synthetic-1N",
        "source_sha256": "c" * 64,
        "maximum_on_time_s": 2.0,
    }
    values.update(changes)
    return PulseModeThruster(**values)


class PulseBitTests(unittest.TestCase):
    def test_nodes_extension_and_derivatives(self) -> None:
        thruster = synthetic_thruster()
        for t, impulse, mass in zip(ON_TIMES, thruster.impulse_bits_n_s,
                                    thruster.propellant_per_pulse_kg, strict=True):
            bit = thruster.pulse(t)
            self.assertAlmostEqual(bit.impulse_n_s, impulse, places=15)
            self.assertAlmostEqual(bit.propellant_kg, mass, places=15)
        # Linear data are reproduced exactly between nodes, derivative included.
        for t in (0.013, 0.037, 0.11, 0.19):
            bit = thruster.pulse(t)
            self.assertAlmostEqual(bit.impulse_n_s, THRUST * (t - 0.006) + 0.002, places=14)
            self.assertAlmostEqual(bit.impulse_wrt_on_time_n, THRUST, places=10)
        # Fully developed beyond the last node: steady thrust and flow.
        bit = thruster.pulse(1.2)
        self.assertAlmostEqual(bit.impulse_n_s, 0.196 + THRUST * 1.0, places=14)
        self.assertAlmostEqual(bit.propellant_kg, FLOW * 0.203 + FLOW * 1.0, places=15)
        self.assertEqual(thruster.minimum_impulse_bit_n_s, thruster.impulse_bits_n_s[0])
        # Short pulses are less efficient than steady state (220 s).
        self.assertAlmostEqual(thruster.pulse(0.01).specific_impulse_s, 220.0 * 0.006 / 0.013, places=10)
        self.assertLess(thruster.pulse(1.2).specific_impulse_s, 220.0)
        for bad in (0.005, 2.5):
            with self.assertRaises(DomainError):
                thruster.pulse(bad)

    def test_nonlinear_table_is_monotone_and_differentiated_exactly(self) -> None:
        impulses = (0.0011, 0.0042, 0.0170, 0.0410, 0.0920)
        masses = (1.6e-6, 4.1e-6, 1.1e-5, 2.4e-5, 4.9e-5)
        thruster = synthetic_thruster(impulse_bits_n_s=impulses, propellant_per_pulse_kg=masses)
        previous = 0.0
        h = 1e-7
        generator = random.Random(4)
        for k in range(2001):
            t = 0.01 + 0.19 * k / 2000
            value = thruster.pulse(t).impulse_n_s
            self.assertGreaterEqual(value, previous - 1e-15)
            previous = value
        for _ in range(100):
            t = generator.uniform(0.0101, 0.1999)
            bit = thruster.pulse(t)
            numeric = (thruster.pulse(t + h).impulse_n_s - thruster.pulse(t - h).impulse_n_s) / (2 * h)
            self.assertAlmostEqual(bit.impulse_wrt_on_time_n, numeric, places=6)
            back = thruster.on_time_for_impulse(bit.impulse_n_s)
            self.assertAlmostEqual(back, t, places=12)
        with self.assertRaises(DomainError):
            thruster.on_time_for_impulse(0.001)


class ImpulsePlanTests(unittest.TestCase):
    def test_exact_requests_use_the_fewest_equal_pulses(self) -> None:
        thruster = synthetic_thruster(maximum_on_time_s=None)
        plan = thruster.plan_impulse(0.5)
        self.assertEqual(plan.pulse_count, 3)  # ceil(0.5 / 0.196)
        self.assertFalse(plan.quantized)
        self.assertAlmostEqual(plan.delivered, 0.5, places=12)
        self.assertAlmostEqual(plan.propellant_kg, 3 * plan.pulse.propellant_kg, places=18)
        self.assertEqual(thruster.plan_impulse(0.0).pulse_count, 0)

    def test_requests_below_the_minimum_impulse_bit_are_quantized(self) -> None:
        thruster = synthetic_thruster()
        nearest = thruster.plan_impulse(0.0041)
        self.assertEqual(nearest.pulse_count, 1)
        self.assertTrue(nearest.quantized)
        self.assertAlmostEqual(nearest.residual, 0.0041 - 0.006, places=15)
        self.assertEqual(thruster.plan_impulse(0.0029).pulse_count, 0)
        self.assertEqual(thruster.plan_impulse(0.0041, quantization="no_overshoot").pulse_count, 0)

    def test_plans_are_optimal_against_brute_force(self) -> None:
        # A short-pulse thruster whose longest pulse is below two minimum bits
        # has gaps between n longest pulses and n + 1 minimum bits.
        thruster = synthetic_thruster(
            on_times_s=(0.010, 0.011, 0.012), impulse_bits_n_s=(0.006, 0.0068, 0.0075),
            propellant_per_pulse_kg=(3e-6, 3.2e-6, 3.4e-6), maximum_on_time_s=None,
        )
        grid = [0.010 + 0.002 * k / 400 for k in range(401)]
        bits = [thruster.pulse(t).impulse_n_s for t in grid]
        generator = random.Random(5)
        for _ in range(300):
            request = generator.uniform(0.0, 0.05)
            plan = thruster.plan_impulse(request)
            best = min([request] + [abs(request - n * b) for n in range(1, 12) for b in bits])
            self.assertLessEqual(abs(plan.residual), best + 1e-12)
            if not plan.quantized:
                exact_counts = [n for n in range(1, 12)
                                if n * bits[0] <= request + 1e-12 <= n * bits[-1] + 2e-12]
                self.assertEqual(plan.pulse_count, min(exact_counts))
            safe = thruster.plan_impulse(request, quantization="no_overshoot")
            self.assertGreaterEqual(safe.residual, -1e-15)
        self.assertEqual(thruster.plan_impulse(0.0085).pulse_count, 1)  # 0.0075 is nearer
        self.assertEqual(thruster.plan_impulse(0.0110).pulse_count, 2)  # 0.0120 is nearer
        with self.assertRaises(DomainError):
            thruster.plan_impulse(1.0, maximum_pulses=10)


class DeltaVPlanTests(unittest.TestCase):
    def test_delta_v_plan_matches_the_rocket_equation_and_the_artifact(self) -> None:
        thruster = synthetic_thruster()
        mass = 12.0
        for request in (0.004, 0.05, 0.7, 3.0):
            plan = thruster.plan_delta_v(request, mass)
            self.assertFalse(plan.quantized, request)
            self.assertAlmostEqual(plan.delivered, request, delta=1e-12 * max(1.0, request))
            expected = thruster.train_delta_v_m_s(plan.pulse_count, plan.pulse.on_time_s, mass)
            self.assertAlmostEqual(expected, request, delta=1e-12 * max(1.0, request))
            # Fewest pulses: one pulse fewer at the longest width falls short.
            if plan.pulse_count > 1:
                self.assertLess(
                    thruster.train_delta_v_m_s(plan.pulse_count - 1, 2.0, mass), request
                )
            if plan.pulse_count > 2000:
                continue
            artifact = thruster.pulse_train_artifact(plan, period_s=plan.pulse.on_time_s + 0.05)
            self.assertEqual(len(artifact.segments), 2 * plan.pulse_count - 1)
            self.assertAlmostEqual(artifact.delivered_total_impulse_n_s,
                                   plan.delivered_impulse_n_s, places=12)
            metrics = free_space_burn_metrics(artifact, mass)
            self.assertAlmostEqual(metrics.delta_v_m_s, request, delta=1e-10 * max(1.0, request))
            self.assertAlmostEqual(metrics.final_mass_kg, plan.final_mass_kg, places=12)

    def test_small_delta_v_is_quantized_and_mass_limits_fail_closed(self) -> None:
        thruster = synthetic_thruster()
        mib_dv = thruster.train_delta_v_m_s(1, 0.010, 12.0)
        plan = thruster.plan_delta_v(0.3 * mib_dv, 12.0)
        self.assertEqual(plan.pulse_count, 0)
        self.assertTrue(plan.quantized)
        plan = thruster.plan_delta_v(0.7 * mib_dv, 12.0)
        self.assertEqual(plan.pulse_count, 1)
        self.assertAlmostEqual(plan.delivered, mib_dv, places=15)
        self.assertEqual(
            thruster.plan_delta_v(0.7 * mib_dv, 12.0, quantization="no_overshoot").pulse_count, 0
        )
        with self.assertRaises(DomainError):
            thruster.plan_delta_v(1e6, 12.0)
        with self.assertRaises(DomainError):
            thruster.plan_delta_v(3.0, 12.0, maximum_pulses=5)

    def test_pulse_train_runs_through_the_bridge_and_lunaris_wrapper(self) -> None:
        from rocket_propulsion.integrations import LunarisFiniteBurnDynamics
        from tests.test_lunaris_integration import _FreeSpaceDynamics, integrate_like_lunaris

        thruster = synthetic_thruster()
        plan = thruster.plan_delta_v(0.05, 12.0)
        artifact = thruster.pulse_train_artifact(plan, period_s=plan.pulse.on_time_s + 0.05)
        bridge = PropagatorPropulsionBridge(artifact, require_complete_stream_binding=False)
        timeline = PropulsionTimeline(burns=(ScheduledBurn("trim", bridge, 5.0),))
        dynamics = LunarisFiniteBurnDynamics(
            _FreeSpaceDynamics(), timeline, lambda burn_id, t, y: (1.0, 0.0, 0.0)
        )
        final = integrate_like_lunaris(dynamics, timeline, [0, 0, 0, 0, 0, 0, 12.0], 0.0,
                                       5.0 + artifact.duration_s + 1.0, 0.01)
        self.assertAlmostEqual(final[3], 0.05, delta=1e-12)
        self.assertAlmostEqual(final[6], plan.final_mass_kg, places=12)


class PulsePersistenceAndValidationTests(unittest.TestCase):
    def test_round_trip_and_hash(self) -> None:
        thruster = synthetic_thruster()
        self.assertEqual(PulseModeThruster.from_json(thruster.to_json()), thruster)
        payload = thruster.manifest()
        payload["steady_thrust_n"] = 2.0
        with self.assertRaises(InputError):
            PulseModeThruster.from_dict(payload)
        with self.assertRaises(InputError):
            PulseModeThruster.from_json("{")

    def test_invalid_models_and_inputs_fail_closed(self) -> None:
        bad_models = (
            {"on_times_s": (0.02, 0.01, 0.05, 0.1, 0.2)},
            {"impulse_bits_n_s": (0.004, 0.004, 0.05, 0.1, 0.2)},
            {"propellant_per_pulse_kg": (2e-6, 1e-6, 3e-6, 4e-6, 5e-6)},
            {"on_times_s": (0.0, 0.02, 0.05, 0.1, 0.2)},
            {"steady_thrust_n": 0.0},
            {"maximum_on_time_s": 0.1},
            {"minimum_off_time_s": -1.0},
            {"source_sha256": "x"},
        )
        for change in bad_models:
            with self.assertRaises(DomainError, msg=str(change)):
                synthetic_thruster(**change)
        thruster = synthetic_thruster()
        plan = thruster.plan_impulse(0.3)
        with self.assertRaises(DomainError):
            thruster.pulse_train_artifact(plan, period_s=plan.pulse.on_time_s + 0.01)
        with self.assertRaises(DomainError):
            thruster.pulse_train_artifact(thruster.plan_impulse(0.0), period_s=1.0)
        accepted = (DomainError, InputError)
        for value in (None, "1", float("nan"), float("inf"), True, [], object(), -1.0):
            for call in (
                lambda v=value: thruster.pulse(v),
                lambda v=value: thruster.plan_impulse(v),
                lambda v=value: thruster.plan_impulse(0.1, quantization=v),
                lambda v=value: thruster.plan_impulse(0.1, maximum_pulses=v),
                lambda v=value: thruster.plan_delta_v(v, 10.0),
                lambda v=value: thruster.plan_delta_v(0.1, v),
                lambda v=value: thruster.pulse_train_artifact(plan, v),
                lambda v=value: thruster.pulse_train_artifact(v, 1.0),
                lambda v=value: thruster.on_time_for_impulse(v),
                lambda v=value: thruster.train_delta_v_m_s(v, 0.05, 10.0),
                lambda v=value: synthetic_thruster(steady_mass_flow_kg_s=v),
                lambda v=value: synthetic_thruster(on_times_s=v),
                lambda v=value: PulseModeThruster.from_dict(v),
            ):
                try:
                    call()
                except accepted:
                    pass
        self.assertIsInstance(plan, PulsePlan)


if __name__ == "__main__":
    unittest.main()
