"""Independent thermodynamic/mass-flux oracles for complete forward ratios.

Extremes check the stated mathematical domain, not a physical gamma envelope.
"""
import json
import math
import unittest
from dataclasses import asdict
from decimal import DefaultContext, FloatOperation, Inexact, localcontext
from unittest.mock import patch

import mpmath as mp

from rocket_propulsion.compressible.isentropic import _area_ratio, calculate_isentropic
from rocket_propulsion.core.errors import DomainError

RATIOS = ('pressure_total_ratio','density_total_ratio','temperature_total_ratio',
          'pressure_critical_ratio','density_critical_ratio','temperature_critical_ratio',
          'area_critical_ratio')


def oracle(mach, gamma):
    """Entropy ratios and continuity at the sonic state; independent area algebra."""
    with mp.workdps(1200):
        m,g=mp.mpf(mach),mp.mpf(gamma)
        temperature=1/(1+(g-1)*m*m/2)
        sonic_temperature=2/(g+1)
        density=mp.exp(mp.log(temperature)/(g-1))
        pressure=mp.exp(g*mp.log(temperature)/(g-1))
        tcrit=temperature/sonic_temperature
        rhocrit=mp.exp(mp.log(tcrit)/(g-1))
        pcrit=mp.exp(g*mp.log(tcrit)/(g-1))
        # rho*M*sqrt(T)*A is constant, compared with its M=1 station.
        area=1/(rhocrit*m*mp.sqrt(tcrit))
        return dict(zip(RATIOS,map(float,(pressure,density,temperature,pcrit,rhocrit,tcrit,area))))


class IsentropicForwardRangeTests(unittest.TestCase):
    def valid(self,mach,gamma):
        try:
            return calculate_isentropic(mach,gamma)
        except (DomainError, ArithmeticError) as error:
            self.fail(f'Finite complete case unexpectedly refused at M={mach}, gamma={gamma}: {error!r}')

    def valid_area(self,mach,gamma):
        try:
            return _area_ratio(mach,gamma)
        except (DomainError, ArithmeticError) as error:
            self.fail(f'Finite complete area unexpectedly refused at M={mach}, gamma={gamma}: {error!r}')

    def assert_ratios(self,mach,gamma,ulps=2):
        expected=oracle(mach,gamma)
        self.assertTrue(all(math.isfinite(v) and v>0 for v in expected.values()), expected)
        actual=asdict(self.valid(mach,gamma))
        for field,value in expected.items():
            with self.subTest(mach=mach,gamma=gamma,field=field):
                self.assertLessEqual(abs(actual[field]-value),ulps*math.ulp(value))
        json.dumps(actual,allow_nan=False)

    def test_balanced_high_gamma_outputs_are_not_refused(self):
        for mach in (2.,3.,10.,1e5,6e7):
            with self.subTest(mach=mach): self.assert_ratios(mach,1e308)

    def test_small_mach_high_gamma_area_does_not_become_inverse_mach(self):
        for mach in (math.ulp(0.),1e-200,1e-160,1e-154,1e-100):
            with self.subTest(mach=mach): self.assert_ratios(mach,1e308)

    def test_nearly_unit_gamma_retains_the_sonic_normalization(self):
        for gamma in (math.nextafter(1.,math.inf),1+1e-12,1.0001,1.001,1.01,1.05):
            for mach in (.25,1.,2.,10.):
                with self.subTest(mach=mach,gamma=gamma):self.assert_ratios(mach,gamma)

    def test_subnormal_rounding_is_not_confused_with_real_zero(self):
        r=self.valid(6e7,1e308)
        self.assertEqual(r.temperature_total_ratio,math.ulp(0.))
        self.assertEqual(r.pressure_total_ratio,math.ulp(0.))
        self.assertGreater(r.temperature_critical_ratio,0)

    def test_genuinely_unrepresentable_results_are_typed_refusals(self):
        for mach,gamma in ((1e8,1e308),(1e150,1e100),(math.ulp(0.),1.4),
                           (100.,math.nextafter(1.,math.inf)),(1e308,1.4)):
            expected=oracle(mach,gamma)
            bad={name for name,value in expected.items() if not math.isfinite(value) or value==0}
            self.assertTrue(bad)
            with self.subTest(mach=mach,gamma=gamma),self.assertRaises(DomainError) as caught:
                calculate_isentropic(mach,gamma)
            self.assertEqual(caught.exception.code,'numeric_range')
            self.assertIn(caught.exception.field,bad)

    def test_ordinary_grid_matches_independent_entropy_and_continuity(self):
        for gamma in (1.1,1.22,1.4,5/3,2.,10.,1e4,1e6):
            for mach in (.01,.25,1.,2.,10.,30.,100.,1000.):
                expected=oracle(mach,gamma)
                if any(not math.isfinite(v) or v<=0 for v in expected.values()):
                    with self.assertRaises(DomainError):calculate_isentropic(mach,gamma)
                else:
                    self.assert_ratios(mach,gamma,ulps=64)

    def test_area_only_solver_does_not_require_unrelated_stagnation_outputs(self):
        # At huge M the static pressure is unrepresentable while the area fits.
        expected=oracle(1e150,1e100)
        self.assertEqual(expected['pressure_total_ratio'],0)
        self.assertEqual(self.valid_area(1e150,1e100),expected['area_critical_ratio'])
        # Native normalized power overflows before division by Mach, although
        # the complete area fits. Area-only inversion must retain that result.
        expected=oracle(3e15,1.1)
        self.assertEqual(self.valid_area(3e15,1.1),expected['area_critical_ratio'])
        # A genuinely unrepresentable area cannot escape as
        # an untyped OverflowError or a successful infinity.
        with self.assertRaises(DomainError) as caught:_area_ratio(100.,1.001)
        self.assertEqual(caught.exception.field,'area_critical_ratio')

    def test_large_power_arguments_do_not_amplify_rounded_exponents(self):
        cases=[(1e3,1.4)]
        cases += [(m,1.1) for m in (1e4,1e6,1e8)]
        cases += [(m,g) for g in (1.4,4.,10.,1e4) for m in (1e4,1e8,1e20,1e40)]
        for mach,gamma in cases:
            with self.subTest(mach=mach,gamma=gamma):self.assert_ratios(mach,gamma)

    def test_active_decimal_environment_is_private_and_unchanged(self):
        baseline=asdict(self.valid(2.,1e308))
        with localcontext() as ctx:
            ctx.prec=5;ctx.Emin=-9;ctx.Emax=9
            ctx.traps[FloatOperation]=True;ctx.traps[Inexact]=True
            before=(ctx.prec,ctx.Emin,ctx.Emax,dict(ctx.flags),dict(ctx.traps))
            self.assertEqual(asdict(self.valid(2.,1e308)),baseline)
            self.assertEqual((ctx.prec,ctx.Emin,ctx.Emax,dict(ctx.flags),dict(ctx.traps)),before)

    def test_process_decimal_defaults_cannot_poison_ratio_context(self):
        baseline=asdict(self.valid(.25,math.nextafter(1.,math.inf)))
        saved=DefaultContext.copy()
        with localcontext():
            try:
                DefaultContext.prec=5;DefaultContext.Emax=9;DefaultContext.Emin=-9
                DefaultContext.traps[Inexact]=True;DefaultContext.traps[FloatOperation]=True
                self.assertEqual(asdict(self.valid(.25,math.nextafter(1.,math.inf))),baseline)
            finally:
                for key in ('prec','rounding','Emin','Emax','capitals','clamp'):
                    setattr(DefaultContext,key,getattr(saved,key))
                DefaultContext.flags=saved.flags;DefaultContext.traps=saved.traps

    def test_ordinary_case_keeps_native_fast_path(self):
        with patch('rocket_propulsion.compressible.isentropic.precise_ratios',
                   side_effect=AssertionError('unexpected precise dispatch')):
            result=calculate_isentropic(2.,1.4)
        self.assertAlmostEqual(result.area_critical_ratio,1.6875,places=13)


if __name__ == '__main__':unittest.main()
