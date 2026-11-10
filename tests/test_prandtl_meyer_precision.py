"""PM references use exact supplied floats; generating Mach is not an inverse oracle."""
import json
import math
import unittest
from decimal import ROUND_DOWN, Decimal, DefaultContext, FloatOperation, Inexact, localcontext
from functools import cache
from unittest.mock import patch

import mpmath as mp

from rocket_propulsion.api.routes import dispatch_calculation
from rocket_propulsion.compressible import sonic_numerics
from rocket_propulsion.compressible.isentropic import (
    _prandtl_meyer_deg,
    calculate_isentropic,
    mach_from_isentropic,
)
from rocket_propulsion.core.errors import ConvergenceError, DomainError


@cache
def angle_reference(mach, gamma):
    """Independent integration of d(nu)/dM, scaled to retain very large gamma."""
    with mp.workdps(120):
        m, g = mp.mpf(mach), mp.mpf(gamma)
        def differential(t):
            x = 1+t*t
            return 2*t*t*mp.sqrt(2+t*t)/x * ((g+1)/(1+(g-1)*x*x/2))
        return float(mp.quad(differential, [0, mp.sqrt(m-1)])/(g+1)*180/mp.pi)


@cache
def inverse_reference(angle, gamma):
    """450-digit library transcendental reference + independent 160-step bisection.

    Admission and roots are calculated from the actual binary64 angle, including
    its quantisation. No production callbacks, Decimal or Newton steps are used.
    """
    with mp.workdps(450):
        g, wanted = mp.mpf(gamma), mp.mpf(angle)
        factor = mp.sqrt((g+1)/(g-1))
        maximum = 180/((g-1)*(factor+1))
        if not 0 <= wanted < maximum:
            return 'physical_domain', None
        def value(m):
            root = mp.sqrt((m-1)*(m+1))
            return (factor*mp.atan(root/factor)-mp.atan(root))*180/mp.pi
        lower, upper = mp.mpf(1), mp.mpf(1000)
        if wanted > value(upper):
            return 'outside_inverse_range', None
        for _ in range(160):
            middle = (lower+upper)/2
            if value(middle) < wanted:
                lower = middle
            else:
                upper = middle
        return 'valid', float((lower+upper)/2)


class PrandtlMeyerPrecisionTests(unittest.TestCase):
    def compare(self, angle, gamma):
        status, expected = inverse_reference(angle, gamma)
        if status != 'valid':
            with self.assertRaises(DomainError) as caught:
                mach_from_isentropic('prandtl_meyer_deg', angle, gamma)
            self.assertEqual(caught.exception.code, status)
            self.assertEqual(caught.exception.field, 'value')
            return
        try:
            actual = mach_from_isentropic('prandtl_meyer_deg', angle, gamma)
        except (DomainError, ConvergenceError, ArithmeticError) as error:
            self.fail(f'Independently admitted PM inverse refused: {angle=}, {gamma=}: {error!r}')
        self.assertTrue(math.isfinite(actual) and 1 <= actual <= 1000)
        self.assertLessEqual(abs(actual-expected), 2*math.ulp(expected))

    def test_two_valid_upper_bracket_angles_are_not_falsely_refused(self):
        for gamma, angle in ((1000., .08993033878768658), (100000., .0008988585731864637)):
            with self.subTest(gamma=gamma):
                self.compare(angle, gamma)

    def test_inaccurate_native_bracket_proposals_cannot_exclude_the_true_root(self):
        for bad_atan in (-1., 1.):
            with (
                self.subTest(native_atan=bad_atan),
                patch.object(sonic_numerics, 'atan', return_value=bad_atan),
            ):
                self.compare(20., 1.4)

    def test_exact_supplied_angles_resolve_beyond_native_callback_plateaus(self):
        for gamma, angle in ((1.4, 130.16759847813037), (10., 8.228568481231619),
                             (10., 9.48601131827985), (1000., .001113615792028615),
                             (1000., .07859364131142432), (100000., .000351901895840347),
                             (100000., .000785603078059918)):
            with self.subTest(gamma=gamma, angle=angle):
                self.compare(angle, gamma)

    def test_upper_boundary_neighbors_classify_before_float_export(self):
        for gamma in (1.4, 10., 1000., 100000., 1e6, 1e308):
            boundary = angle_reference(1000., gamma)
            for angle in (math.nextafter(boundary, 0.), boundary,
                          math.nextafter(boundary, math.inf)):
                with self.subTest(gamma=gamma, angle=angle):
                    self.compare(angle, gamma)

    def test_forward_angles_match_independent_differential_quadrature(self):
        for gamma in (1.4, 10., 1000., 100000., 1e308):
            for mach in (1.02, 2., 10., 1000.):
                with self.subTest(gamma=gamma, mach=mach):
                    expected = angle_reference(mach, gamma)
                    actual = _prandtl_meyer_deg(mach, gamma)
                    self.assertLessEqual(abs(actual-expected), 2*math.ulp(expected))

    def test_physical_refusals_identify_the_supplied_value(self):
        for gamma, angle in ((1.4, -1.), (100000., 1.), (1e308, 1e-305)):
            with self.subTest(gamma=gamma, angle=angle):
                self.compare(angle, gamma)

    def test_exact_input_sweep_covers_old_dispatch_thresholds_and_extreme_gamma(self):
        for gamma in (math.nextafter(1., math.inf), 1.001, 1.4, 10., 1000.,
                      100000., 1e6, 1e100, 1e308):
            for mach in (math.nextafter(1., math.inf), 1.0001, 1.01,
                         math.nextafter(1.01, math.inf), 2., 10., 999., 1000.):
                angle = angle_reference(mach, gamma)
                with self.subTest(gamma=gamma, generating_mach=mach, supplied_angle=angle):
                    self.compare(angle, gamma)
        for gamma in (1.4, 1e100, 1e308):
            for angle in (0., math.ulp(0.), 1e-320, 1e-100,
                          math.nextafter(1e-4, 0.), 1e-4, math.nextafter(1e-4, math.inf)):
                with self.subTest(gamma=gamma, supplied_angle=angle):
                    self.compare(angle, gamma)

    def test_public_forward_keeps_the_independent_angle(self):
        for gamma, mach in ((1.4, 2.), (1000., 1.05), (100000., 10.), (1e308, 1000.)):
            with self.subTest(gamma=gamma, mach=mach):
                expected = angle_reference(mach, gamma)
                try:
                    actual = calculate_isentropic(mach, gamma).prandtl_meyer_deg
                except (DomainError, ConvergenceError, ArithmeticError) as error:
                    self.fail(f'Independently finite complete state refused: {error!r}')
                self.assertLessEqual(abs(actual-expected), 2*math.ulp(expected))

    def test_nonzero_unrepresentable_forward_angle_has_a_typed_refusal(self):
        mach = math.nextafter(1., math.inf)
        self.assertEqual(angle_reference(mach, 1e308), 0.)
        with self.assertRaises(DomainError) as caught:
            _prandtl_meyer_deg(mach, 1e308)
        self.assertEqual(caught.exception.code, 'numeric_range')
        self.assertEqual(caught.exception.field, 'prandtl_meyer_deg')
        self.assertEqual(_prandtl_meyer_deg(1., 1e308), 0.)

    def test_pm_inverse_owns_active_context_on_success_and_refusal(self):
        sonic_numerics._decimal_pi.cache_clear()
        try:
            with localcontext() as ctx:
                ctx.prec = 5
                ctx.rounding = ROUND_DOWN
                ctx.Emin, ctx.Emax = -9, 9
                ctx.traps[FloatOperation] = ctx.traps[Inexact] = True
                before = context_state(ctx)
                self.compare(.07859364131142432, 1000.)
                self.compare(1., 100000.)
                self.compare(9.49, 10.)
                self.assertEqual(context_state(ctx), before)
                expected = angle_reference(10., 100000.)
                try:
                    actual = _prandtl_meyer_deg(10., 100000.)
                except (DomainError, ConvergenceError, ArithmeticError) as error:
                    self.fail(f'Caller context controlled PM evaluation: {error!r}')
                self.assertLessEqual(abs(actual-expected), 2*math.ulp(expected))
                self.assertEqual(context_state(ctx), before)
        finally:
            sonic_numerics._decimal_pi.cache_clear()

    def test_pi_cache_seeded_by_low_precision_caller_cannot_corrupt_forward_or_inverse(self):
        gamma, mach = 100000., 2.
        supplied = angle_reference(mach, gamma)
        # Use the precision shared by these forward and inverse calls. Seed the
        # cache before either public entry point can install its own context.
        precision = sonic_numerics._precision(gamma)
        with mp.workdps(precision + 20):
            expected_pi = Decimal(mp.nstr(mp.pi, precision + 12))
        tolerance = Decimal(f'1e{2 - precision}')
        sonic_numerics._decimal_pi.cache_clear()
        try:
            with localcontext() as ctx:
                ctx.prec = 5
                ctx.rounding = ROUND_DOWN
                ctx.Emin, ctx.Emax = -9, 9
                before = context_state(ctx)
                seeded_pi = sonic_numerics._decimal_pi(precision)
                self.assertEqual(context_state(ctx), before)

            with localcontext() as ctx:
                ctx.prec = precision + 20
                self.assertLessEqual(abs(seeded_pi - expected_pi), tolerance)
            actual = calculate_isentropic(mach, gamma).prandtl_meyer_deg
            self.assertLessEqual(abs(actual - supplied), 2 * math.ulp(supplied))
            self.compare(supplied, gamma)
        finally:
            sonic_numerics._decimal_pi.cache_clear()

    def test_pi_cache_misses_own_mutable_default_context(self):
        sonic_numerics._decimal_pi.cache_clear()
        with localcontext():
            saved = DefaultContext.copy()
            try:
                DefaultContext.prec = 5
                DefaultContext.rounding = ROUND_DOWN
                DefaultContext.Emin, DefaultContext.Emax = -9, 9
                DefaultContext.traps[FloatOperation] = DefaultContext.traps[Inexact] = True
                before = context_state(DefaultContext)
                self.compare(.000351901895840347, 100000.)
                self.assertEqual(context_state(DefaultContext), before)
            finally:
                for name in ('prec', 'rounding', 'Emin', 'Emax', 'capitals', 'clamp'):
                    setattr(DefaultContext, name, getattr(saved, name))
                DefaultContext.traps = saved.traps.copy()
                DefaultContext.flags = saved.flags.copy()
                sonic_numerics._decimal_pi.cache_clear()

    def test_exhausted_inverse_does_not_export_an_unresolved_iterate(self):
        with (patch.object(sonic_numerics, '_PM_INVERSE_ITERATIONS', 1),
              self.assertRaises(ConvergenceError) as caught):
            mach_from_isentropic('prandtl_meyer_deg', 20., 1.4)
        self.assertEqual(caught.exception.code, 'prandtl_meyer_not_converged')
        self.assertEqual(caught.exception.field, 'value')

    def test_failed_angle_evaluation_preserves_the_inverse_failure_cause(self):
        original = ConvergenceError('Injected PM evaluation failure', field='prandtl_meyer_deg')
        with (patch.object(sonic_numerics, '_decimal_pi', side_effect=original),
              self.assertRaises(ConvergenceError) as caught):
            mach_from_isentropic('prandtl_meyer_deg', 20., 1.4)
        self.assertEqual(caught.exception.code, 'prandtl_meyer_not_converged')
        self.assertEqual(caught.exception.field, 'value')
        self.assertIs(caught.exception.__cause__, original)

    def test_forward_series_exhaustion_is_never_a_partial_angle(self):
        for budget, mach in (('_PM_SERIES_ITERATIONS', 1.0001), ('_ATAN_SERIES_ITERATIONS', 2.)):
            sonic_numerics._decimal_pi.cache_clear()
            try:
                with self.subTest(budget=budget):
                    with (patch.object(sonic_numerics, budget, 1),
                          self.assertRaises(ConvergenceError) as caught):
                        _prandtl_meyer_deg(mach, 1.4)
                    self.assertEqual(caught.exception.field, 'prandtl_meyer_deg')
            finally:
                sonic_numerics._decimal_pi.cache_clear()

    def test_real_route_keeps_valid_and_typed_refused_pm_inputs(self):
        for gamma, angle in ((1000., .08993033878768658), (100000., .0008988585731864637)):
            with self.subTest(gamma=gamma):
                try:
                    answer = dispatch_calculation('/api/v1/isentropic',
                        {'input_kind': 'prandtl_meyer_deg', 'value': angle, 'gamma': gamma})
                except (DomainError, ConvergenceError, ArithmeticError) as error:
                    self.fail(f'Valid PM route input refused: {error!r}')
                expected = inverse_reference(angle, gamma)[1]
                self.assertLessEqual(abs(answer['mach']-expected), 2*math.ulp(expected))
                json.dumps(answer, allow_nan=False)
        for angle, gamma, code in ((-1., 1.4, 'physical_domain'),
                                   (9.49, 10., 'outside_inverse_range')):
            with self.subTest(angle=angle):
                with self.assertRaises(DomainError) as caught:
                    dispatch_calculation('/api/v1/isentropic',
                        {'input_kind': 'prandtl_meyer_deg', 'value': angle, 'gamma': gamma})
                self.assertEqual(caught.exception.code, code)
                self.assertEqual(caught.exception.field, 'value')


def context_state(ctx):
    return (ctx.prec, ctx.rounding, ctx.Emin, ctx.Emax, ctx.capitals, ctx.clamp,
            dict(ctx.flags), dict(ctx.traps))


if __name__ == '__main__':
    unittest.main()
