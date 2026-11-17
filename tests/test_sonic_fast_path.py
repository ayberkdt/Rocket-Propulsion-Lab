"""The bounded first stage must return exactly what the complete stage returns."""
import math
import random
import unittest
from decimal import Decimal
from unittest.mock import patch

from rocket_propulsion.compressible import sonic_numerics as sn
from rocket_propulsion.compressible.isentropic import mach_from_isentropic
from rocket_propulsion.core.errors import RocketPropulsionError

GAMMAS = (1.05, 1.13, 1.22, 1.4, 5/3, 3.0, 10.0, 1e3, 1e6)


def outcome(function, *args, **kwargs):
    try:
        return ('value', function(*args, **kwargs))
    except (RocketPropulsionError, ArithmeticError) as error:  # refusals must also agree
        return ('refusal', type(error).__name__, getattr(error, 'code', None), getattr(error, 'field', None))


def complete(function, *args, **kwargs):
    with patch.object(sn, '_FAST_PATH', False):
        return outcome(function, *args, **kwargs)


def machs(rng, count):
    fixed = [math.nextafter(1., 2.), 1+1e-8, 1.0001, 1.009, 1.01, 1.0100001, 1.1, 1.1000001, 2., 7.5, 99.9, 999.]
    return fixed+[math.exp(rng.uniform(0., math.log(1000.))) for _ in range(count)]


class BoundedFirstStageTests(unittest.TestCase):
    def test_forward_and_inverse_values_match_the_complete_stage_bit_for_bit(self):
        rng = random.Random(20261117)
        for gamma in GAMMAS:
            for mach in machs(rng, 12):
                angle = complete(sn.sonic_pm, mach, gamma)
                with self.subTest(kind='pm', gamma=gamma, mach=mach):
                    self.assertEqual(outcome(sn.sonic_pm, mach, gamma), angle)
                if angle[0] == 'value' and angle[1] > 0:
                    with self.subTest(kind='pm_inverse', gamma=gamma, angle=angle[1]):
                        self.assertEqual(outcome(mach_from_isentropic, 'prandtl_meyer_deg', angle[1], gamma),
                                         complete(mach_from_isentropic, 'prandtl_meyer_deg', angle[1], gamma))
            for ratio in (1+1e-8, 1.0001, 1.5, 4., 25., 1e4, rng.uniform(1., 50.), rng.uniform(50., 1e5)):
                for kind in ('area_ratio_subsonic', 'area_ratio_supersonic'):
                    with self.subTest(kind=kind, gamma=gamma, ratio=ratio):
                        self.assertEqual(outcome(mach_from_isentropic, kind, ratio, gamma),
                                         complete(mach_from_isentropic, kind, ratio, gamma))

    def test_ordinary_nozzle_states_are_answered_by_the_first_stage(self):
        # A silent fallback would preserve values but lose the speed-up.
        for gamma in (1.13, 1.22, 1.4):
            for mach in (1.2, 2., 3.5, 6.):
                angle = complete(sn.sonic_pm, mach, gamma)[1]
                with self.subTest(gamma=gamma, mach=mach):
                    self.assertIsNotNone(sn._fast_pm(mach, gamma))
                    self.assertIsNotNone(sn._fast_pm_inverse(angle, gamma, 1000.))
                    self.assertIsNotNone(sn._fast_area_inverse(mach*2, gamma, 1e-8, 1000., False))
                    self.assertIsNotNone(sn._fast_area_inverse(mach*2, gamma, 1e-8, 1000., True))

    def test_uncertain_signs_and_wide_intervals_defer_to_the_complete_stage(self):
        self.assertEqual(sn._sign(Decimal('1e-40'), Decimal('1e-30')), 0)
        self.assertEqual(sn._sign(Decimal('-1e-20'), Decimal('1e-30')), -1)
        self.assertEqual(sn._sign(Decimal('1e-20'), Decimal('1e-30')), 1)
        # An allowance this wide cannot prove one binary64 value.
        with patch.object(sn, '_BOUND_MARGIN', 26):
            self.assertIsNone(sn._fast_pm(2.3, 1.22))
            self.assertIsNone(sn._fast_pm_inverse(30., 1.22, 1000.))
            self.assertIsNone(sn._fast_area_inverse(4., 1.22, 1e-8, 1000., False))
            self.assertEqual(outcome(mach_from_isentropic, 'area_ratio_supersonic', 4., 1.22),
                             complete(mach_from_isentropic, 'area_ratio_supersonic', 4., 1.22))

    def test_table_logarithm_stays_within_a_few_units_of_the_working_precision(self):
        from decimal import Context, localcontext
        rng = random.Random(7)
        values = [Decimal(1), Decimal(2), Decimal('1.99999999999999999999'), Decimal(2)**-30, Decimal(2)**40,
                  Decimal(1)+Decimal(1)/256, Decimal(1)+Decimal(255)/256, Decimal('7.31'), Decimal('1e-300'),
                  Decimal('1e400'), Decimal('1e-400')]
        values += [Decimal(math.exp(rng.uniform(-40., 40.))) for _ in range(60)]
        for precision in (30, 36):
            for value in values:
                with self.subTest(precision=precision, value=value):
                    with localcontext(Context(prec=precision, Emin=-999999, Emax=999999)):
                        actual = sn._table_ln(value)
                    with localcontext(Context(prec=precision+30, Emin=-999999, Emax=999999)):
                        error = abs(actual-value.ln())
                        self.assertLessEqual(error, (abs(value.ln())+1)*Decimal(10)**(3-precision))

    def test_midpoints_are_exact_half_ulp_boundaries(self):
        for value in (1., 1.5, 999.0, math.nextafter(1., 2.)):
            upper = math.nextafter(value, math.inf)
            self.assertEqual(sn._midpoint(value, upper)*2, Decimal(value)+Decimal(upper))

    def test_wrong_native_estimates_are_corrected_or_deferred(self):
        expected = complete(mach_from_isentropic, 'prandtl_meyer_deg', 20., 1.4)
        for bad in (-1., 1.):
            with self.subTest(native_atan=bad), patch.object(sn, 'atan', return_value=bad):
                self.assertEqual(outcome(mach_from_isentropic, 'prandtl_meyer_deg', 20., 1.4), expected)
        expected = complete(mach_from_isentropic, 'area_ratio_supersonic', 7., 1.22)
        with patch.object(sn, '_native_root', return_value=1.5):
            self.assertEqual(outcome(mach_from_isentropic, 'area_ratio_supersonic', 7., 1.22), expected)


if __name__ == '__main__':
    unittest.main()
