"""The bounded normal-shock first stage must return what the 80-digit stage returns."""
import math
import random
import unittest
from unittest.mock import patch

from rocket_propulsion.compressible import normal_shock as ns
from rocket_propulsion.core.errors import RocketPropulsionError


def outcome(mach, gamma):
    try:
        return ('value', ns.calculate_normal_shock(mach, gamma))
    except (RocketPropulsionError, ArithmeticError) as error:
        return ('refusal', type(error).__name__, getattr(error, 'code', None), getattr(error, 'field', None))


def complete(mach, gamma):
    with patch.object(ns, '_FAST_PATH', False):
        return outcome(mach, gamma)


def first_stage(mach, gamma):
    token = ns._FIRST_STAGE.set(True)
    try:
        return ns._evaluate(mach, gamma, ns._FIRST_STAGE_DIGITS)
    finally:
        ns._FIRST_STAGE.reset(token)


class NormalShockFirstStageTests(unittest.TestCase):
    def test_outputs_and_refusals_match_the_complete_stage_bit_for_bit(self):
        rng = random.Random(20261117)
        gammas = (1.0001, 1.05, math.nextafter(1.05, 2.), 1.13, 1.22, 1.4, 5/3, 3., 99., 100., 101., 1e6)
        machs = [math.nextafter(1., 2.), 1+1e-12, 1+1e-6, 1.001, 1.5, 2., 7.5, 1e3, 1e6, math.nextafter(1e6, 2e6), 1e60]
        machs += [10**rng.uniform(0., 6.) for _ in range(40)]
        for gamma in gammas:
            for mach in machs:
                with self.subTest(gamma=gamma, mach=mach):
                    self.assertEqual(outcome(mach, gamma), complete(mach, gamma))

    def test_ordinary_states_are_proved_by_the_first_stage(self):
        # A silent fallback would preserve values but lose the speed-up.
        for gamma in (1.13, 1.22, 1.4, 5/3):
            for mach in (1.2, 2., 3.5, 8., 25.):
                with self.subTest(gamma=gamma, mach=mach):
                    self.assertEqual(first_stage(mach, gamma), complete(mach, gamma)[1])

    def test_wide_allowance_defers_every_output_to_the_complete_stage(self):
        # A 1e-3 relative interval spans many positive binary64 values.
        with patch.object(ns, '_ALLOWANCE_DIGITS', ns._FIRST_STAGE_DIGITS-3), self.assertRaises(ns._Unproved):
            first_stage(2., 1.4)
        with patch.object(ns, '_ALLOWANCE_DIGITS', ns._FIRST_STAGE_DIGITS), self.assertRaises(ns._Unproved):
            first_stage(2., 1.4)
        with patch.object(ns, '_ALLOWANCE_DIGITS', ns._FIRST_STAGE_DIGITS):
            self.assertEqual(outcome(2., 1.4), complete(2., 1.4))

    def test_stage_marker_is_restored_after_success_and_fallback(self):
        ns.calculate_normal_shock(2., 1.4)
        self.assertFalse(ns._FIRST_STAGE.get())
        with patch.object(ns, '_ALLOWANCE_DIGITS', ns._FIRST_STAGE_DIGITS):
            ns.calculate_normal_shock(2., 1.4)
        self.assertFalse(ns._FIRST_STAGE.get())
        with self.assertRaises(RocketPropulsionError):
            ns.calculate_normal_shock(1e155, 10.)
        self.assertFalse(ns._FIRST_STAGE.get())


if __name__ == '__main__':
    unittest.main()
