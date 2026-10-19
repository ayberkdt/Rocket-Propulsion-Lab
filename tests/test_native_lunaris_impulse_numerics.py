"""Native impulse norm/mass-law stress with an independent exponential series."""

import importlib.util
import json
import math
import random
import unittest
from decimal import Decimal, localcontext
from fractions import Fraction
from sys import float_info

import numpy as np


def reference(vector, mass, isp, precision):
    with localcontext() as context:
        context.prec = precision
        squared = sum((Fraction(float(value)) ** 2 for value in vector), Fraction(0))
        length = (Decimal(squared.numerator) / Decimal(squared.denominator)).sqrt()
        if isp is None:
            return float(length), mass
        from lunaris.common.constants import STANDARD_GRAVITY

        exponent = length / (Decimal.from_float(STANDARD_GRAVITY) * Decimal.from_float(isp))
        reductions = 0
        reduced = exponent
        while reduced > Decimal(".125"):
            reduced /= 2
            reductions += 1
        term = total = Decimal(1)
        for index in range(1, 5000):
            term *= -reduced / index
            following = total + term
            if following == total:
                break
            total = following
        else:
            raise AssertionError("independent exponential series did not settle")
        for _ in range(reductions):
            total *= total
        return float(length), float(Decimal.from_float(mass) * total)


@unittest.skipUnless(
    importlib.util.find_spec("lunaris") is not None, "Explicit Lunaris source path required"
)
class NativeLunarisImpulseAttacks(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        from lunaris.core.propagation import plans

        cls.native = plans

    def state(self, mass=1000.0):
        values = [2e6, 0.0, 0.0, 0.0, 1.0, 0.0]
        return np.array(values if mass is None else [*values, mass])

    def close(self, value, expected):
        self.assertTrue(math.isfinite(value))
        self.assertLessEqual(abs(value - expected), 4 * math.ulp(expected), (value, expected))

    def oracle(self, vector, mass, isp):
        value = reference(vector, mass, isp, 300)
        self.assertEqual(value, reference(vector, mass, isp, 200))
        return value

    def apply(self, vector, mass=1000.0, isp=None, frame="inertial", state=None):
        state = self.state(mass) if state is None else state
        maneuver = self.native.ImpulsiveManeuver(10.0, tuple(vector), frame=frame, isp_s=isp)
        return self.native.apply_impulsive_maneuver(state, maneuver)

    def test_large_finite_norm_and_mass_match_independent_series(self):
        vector = (1e200, -2e200, 2e200)
        expected_norm, expected_mass = self.oracle(vector, 1000.0, 1e200)
        final, record = self.apply(vector, isp=1e200)
        self.close(record["dv_norm_mps"], expected_norm)
        self.close(final[6], expected_mass)

    def test_subnormal_vector_keeps_nonzero_norm_and_finite_mass_ratio(self):
        tiny = float.fromhex("0x0.0000000000001p-1022")
        vector = (tiny, -tiny, tiny)
        expected_norm, expected_mass = self.oracle(vector, 1000.0, tiny)
        final, record = self.apply(vector, isp=tiny)
        self.assertGreater(record["dv_norm_mps"], 0.0)
        self.close(record["dv_norm_mps"], expected_norm)
        self.close(final[6], expected_mass)

    def test_large_isp_product_is_not_rounded_to_infinity(self):
        vector = (1e308, 0.0, 0.0)
        expected_norm, expected_mass = self.oracle(vector, 1000.0, 1e308)
        final, record = self.apply(vector, isp=1e308)
        self.close(record["dv_norm_mps"], expected_norm)
        self.close(final[6], expected_mass)
        self.assertLess(final[6], 1000.0)

    def test_small_ratio_factors_do_not_erase_a_finite_exponent(self):
        vector = (3e-300, 4e-300, 0.0)
        expected_norm, expected_mass = self.oracle(vector, 1000.0, 1e-300)
        final, record = self.apply(vector, isp=1e-300)
        self.close(record["dv_norm_mps"], expected_norm)
        self.close(final[6], expected_mass)

    def test_exponential_underflow_can_recover_a_representable_final_mass(self):
        from lunaris.common.constants import STANDARD_GRAVITY

        for mass, exponent in ((1e308, 1000.0), (1e200, 1200.0)):
            vector = (exponent * STANDARD_GRAVITY, 0.0, 0.0)
            with self.subTest(mass=mass, exponent=exponent):
                expected_norm, expected_mass = self.oracle(vector, mass, 1.0)
                self.assertGreater(expected_mass, 0.0)
                final, record = self.apply(vector, mass=mass, isp=1.0)
                self.close(record["dv_norm_mps"], expected_norm)
                self.close(final[6], expected_mass)

    def test_genuine_final_mass_underflow_has_an_explicit_error(self):
        tiny = float.fromhex("0x0.0000000000001p-1022")
        with self.assertRaisesRegex(ValueError, "post-maneuver mass"):
            self.apply((10.0, 0.0, 0.0), mass=tiny, isp=1.0)

    def test_genuine_scalar_norm_overflow_has_an_explicit_error(self):
        with self.assertRaisesRegex(ValueError, "delta-v norm"):
            self.apply((float_info.max,) * 3, mass=None)

    def test_genuine_velocity_update_overflow_has_an_explicit_error(self):
        state = self.state(None)
        state[3] = float_info.max
        with self.assertRaisesRegex(ValueError, "velocity update"):
            self.apply((float_info.max, 0.0, 0.0), mass=None, state=state)

    def test_ordinary_ric_rotation_and_mass_law_are_preserved(self):
        state = np.array((0.0, 2e6, 0.0, -1500.0, 0.0, 0.0, 500.0))
        command = (2.0, 3.0, 4.0)
        expected_dv = (-3.0, 2.0, 4.0)
        expected_norm, expected_mass = self.oracle(expected_dv, 500.0, 320.0)
        final, record = self.apply(command, mass=500.0, isp=320.0, frame="ric", state=state)
        np.testing.assert_array_equal(final[3:6], state[3:6] + np.array(expected_dv))
        self.close(record["dv_norm_mps"], expected_norm)
        self.close(final[6], expected_mass)

    def test_zero_delta_v_and_massless_kinematics_remain_valid(self):
        for mass, isp in ((500.0, 320.0), (500.0, None), (None, None)):
            with self.subTest(mass=mass, isp=isp):
                state = self.state(mass)
                final, record = self.apply((0.0, 0.0, 0.0), mass=mass, isp=isp)
                np.testing.assert_array_equal(final, state)
                self.assertEqual(record["dv_norm_mps"], 0.0)
                self.assertEqual(record["mass_before_kg"], record["mass_after_kg"])

    def test_malformed_state_and_missing_mass_contract_still_fail(self):
        for state in (np.zeros(5), np.full(7, float("nan")), self.state(None), self.state(0.0)):
            with self.subTest(state=state), self.assertRaises(ValueError):
                self.apply((1.0, 0.0, 0.0), isp=320.0, state=state)

    def test_seeded_vectors_and_mass_scales_match_independent_paths(self):
        random_source = random.Random(5801)
        from lunaris.common.constants import STANDARD_GRAVITY

        for _ in range(24):
            scale = math.ldexp(1.0, random_source.randint(-600, 600))
            vector = tuple(scale * random_source.uniform(-1.0, 1.0) for _ in range(3))
            mass = math.ldexp(1.0, random_source.randint(-300, 300))
            target_exponent = random_source.choice((1e-4, 0.125, 1.0, 10.0))
            isp = math.hypot(*vector) / STANDARD_GRAVITY / target_exponent
            expected_norm, expected_mass = self.oracle(vector, mass, isp)
            with self.subTest(vector=vector, mass=mass):
                final, record = self.apply(vector, mass=mass, isp=isp)
                self.close(record["dv_norm_mps"], expected_norm)
                self.close(final[6], expected_mass)

    def test_records_match_state_and_remain_strict_json(self):
        final, record = self.apply((1e200, 0.0, 0.0), isp=1e200)
        reopened = json.loads(json.dumps(record, allow_nan=False))
        self.assertEqual(reopened["mass_after_kg"], float(final[6]))
        self.assertEqual(reopened["mass_before_kg"], 1000.0)
        self.assertEqual(reopened["dv_command_mps"], [1e200, 0.0, 0.0])
        self.assertEqual(reopened["output_state_semantics"], "post_burn_single_timestamp")

    def test_success_and_error_do_not_mutate_the_callers_state(self):
        for vector in ((1.0, 0.0, 0.0), (float_info.max,) * 3):
            state = self.state()
            before = state.copy()
            try:
                self.apply(vector, isp=320.0, state=state)
            except ValueError:
                pass
            np.testing.assert_array_equal(state, before)

    def test_rocketprop_native_reproduction_retains_declared_mass(self):
        from rocket_propulsion.integrations import lunaris_impulsive_equivalents
        from rocket_propulsion.propulsion.burns import PropulsionTimeline, ScheduledBurn
        from tests.test_lunaris_integration import mass_only_bridge

        bridge = mass_only_bridge(1e300, 0.1, 1.0)
        timeline = PropulsionTimeline((ScheduledBurn("burn", bridge, 0.0),))
        result = lunaris_impulsive_equivalents(
            timeline, initial_mass_kg=1000.0, directions={"burn": ((1.0, 0.0, 0.0), "inertial")}
        )[0]
        final, record = self.native.apply_impulsive_maneuver(self.state(), result.to_lunaris())
        self.assertLessEqual(
            abs(final[6] - result.final_mass_kg), 8 * math.ulp(result.final_mass_kg)
        )
        self.close(record["dv_norm_mps"], result.delta_v_m_s)
