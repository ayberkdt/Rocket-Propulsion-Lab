"""Invert supplied binary64 ratios against an independent energy/entropy oracle.

The oracle interprets the supplied float, not an unrounded generating Mach.
Synthetic gamma scales test arithmetic contracts, not real-gas applicability.
"""
import json
import math
import unittest
from decimal import ROUND_DOWN, DefaultContext, FloatOperation, Inexact, localcontext
from unittest.mock import patch

import mpmath as mp

from rocket_propulsion.api.routes import dispatch_calculation
from rocket_propulsion.compressible.isentropic import mach_from_isentropic
from rocket_propulsion.core.errors import DomainError

KINDS=tuple(f'{quantity}_{reference}_ratio'
            for reference in ('total','critical') for quantity in ('temperature','pressure','density'))
GAMMAS=(math.nextafter(1.,math.inf),1+1e-12,1.01,1.1,1.4,5/3,4.,10.,1e4,1e100,1e308)


def oracle(kind,value,gamma):
    """Entropy gives temperature; total enthalpy then gives velocity / sound speed."""
    with mp.workdps(1200):
        g,r=mp.mpf(gamma),mp.mpf(value)
        power=(g-1)/g if kind.startswith('pressure') else g-1 if kind.startswith('density') else 1
        temperature=r**power
        if '_critical_' in kind:temperature*=2/(g+1)
        squared=2*(1/temperature-1)/(g-1)
        if squared<0:return 'physical_domain',None
        mach=mp.sqrt(squared)
        if not mp.mpf(1e-8)<=mach<=mp.mpf(1000.):return 'outside_inverse_range',None
        return 'valid',float(mach)


def boundary_ratio(kind,mach,gamma):
    with mp.workdps(1200):
        g,m=mp.mpf(gamma),mp.mpf(mach)
        temperature=1/(1+(g-1)*m*m/2)
        if '_critical_' in kind:temperature/=2/(g+1)
        power=g/(g-1) if kind.startswith('pressure') else 1/(g-1) if kind.startswith('density') else 1
        return float(temperature**power)


class IsentropicRatioInverseTests(unittest.TestCase):
    def valid(self,kind,value,gamma):
        try:return mach_from_isentropic(kind,value,gamma)
        except (DomainError,ArithmeticError) as error:
            self.fail(f'Valid supplied ratio refused: {kind}, {value}, gamma={gamma}: {error!r}')

    def compare(self,kind,value,gamma):
        status,expected=oracle(kind,value,gamma)
        if status!='valid':
            try:
                with self.assertRaises(DomainError) as caught:mach_from_isentropic(kind,value,gamma)
            except ArithmeticError as error:
                self.fail(f'Expected a typed {status} refusal, got {error!r}')
            self.assertEqual(caught.exception.code,status)
            self.assertEqual(caught.exception.field,'value')
            return
        actual=self.valid(kind,value,gamma)
        self.assertIsInstance(actual,float)
        self.assertTrue(math.isfinite(actual) and 1e-8<=actual<=1000.)
        self.assertLessEqual(abs(actual-expected),2*math.ulp(expected))
        json.dumps({'mach':actual},allow_nan=False)

    def test_balanced_extreme_and_unit_gamma_false_refusals(self):
        for kind,value,gamma in [('temperature_total_ratio',5e-309,1e308),
                                 ('temperature_critical_ratio',.25,1e308),
                                 ('temperature_critical_ratio',1.,GAMMAS[0]),
                                 ('pressure_critical_ratio',1.,GAMMAS[0])]:
            with self.subTest(kind=kind,value=value,gamma=gamma):self.compare(kind,value,gamma)

    def test_sonic_critical_inputs_keep_the_complete_normalization(self):
        for gamma in GAMMAS:
            for kind in KINDS[3:]:
                with self.subTest(kind=kind,gamma=gamma):
                    self.compare(kind,1.,gamma)
                    self.assertEqual(mach_from_isentropic(kind,1.,gamma),1.)

    def test_raw_ratio_matrix_uses_exact_supplied_floats(self):
        values=(math.ulp(0.),1e-308,1e-100,1e-10,.1,.5,
                math.nextafter(1.,0.),1.,math.nextafter(1.,math.inf),1.5,2.,1e10,1e100,1e308)
        for gamma in GAMMAS:
            for kind in KINDS:
                for value in values:
                    with self.subTest(kind=kind,value=value,gamma=gamma):self.compare(kind,value,gamma)

    def test_neighboring_ratios_at_both_inverse_boundaries(self):
        for gamma in (GAMMAS[0],1.1,1.4,4.,1e4,1e308):
            for kind in KINDS:
                for mach in (1e-8,1000.):
                    value=boundary_ratio(kind,mach,gamma)
                    if not math.isfinite(value) or value<=0:continue  # No supplied binary64 ratio.
                    for neighbor in (math.nextafter(value,0.),value,math.nextafter(value,math.inf)):
                        if not math.isfinite(neighbor) or neighbor<=0:continue
                        with self.subTest(kind=kind,value=neighbor,gamma=gamma):self.compare(kind,neighbor,gamma)

    def test_unity_total_ratio_does_not_invent_a_small_positive_mach(self):
        for gamma in GAMMAS:
            for kind in KINDS[:3]:
                with self.subTest(kind=kind,gamma=gamma):self.compare(kind,1.,gamma)

    def test_pressure_and_density_exponents_are_distinct(self):
        for gamma in (1.01,1.4,4.,1e4):
            for kind in KINDS:
                for value in (.1,.5,math.nextafter(1.,0.),1.5):
                    with self.subTest(kind=kind,value=value,gamma=gamma):self.compare(kind,value,gamma)

    def test_physical_domain_is_not_the_inverse_bracket(self):
        for kind in KINDS:
            with self.subTest(kind=kind):self.compare(kind,100.,1.4)

    def test_huge_private_exponents_are_refused_before_power(self):
        for kind in ('density_total_ratio','density_critical_ratio'):
            with self.subTest(kind=kind):self.compare(kind,.5,1e308)

    def test_finite_roots_outside_the_inverse_bracket_are_refused(self):
        for kind,value,gamma in [('temperature_total_ratio',1e-250,1e308),
                                 ('pressure_total_ratio',1e-250,1e308),
                                 ('temperature_total_ratio',1e-30,1.4),
                                 ('density_total_ratio',1e-100,4.)]:
            with self.subTest(kind=kind,value=value,gamma=gamma):self.compare(kind,value,gamma)

    def test_direct_mach_identity_has_no_derived_inverse_bracket(self):
        for value in (math.ulp(0.),1e-12,1e4,1e308):
            with self.subTest(value=value):self.assertEqual(mach_from_isentropic('mach',value,1.4),value)

    def test_ratio_inverse_does_not_require_unrelated_forward_outputs(self):
        with patch('rocket_propulsion.compressible.isentropic.calculate_isentropic',
                   side_effect=AssertionError('An inverse ratio must stand alone')):
            self.compare('pressure_total_ratio',math.ulp(0.),GAMMAS[0])

    def test_active_decimal_environment_is_private_and_unchanged(self):
        baseline=self.valid('temperature_total_ratio',5e-309,1e308)
        with localcontext() as ctx:
            ctx.prec=5;ctx.rounding=ROUND_DOWN;ctx.Emin=-9;ctx.Emax=9
            ctx.traps[FloatOperation]=ctx.traps[Inexact]=True
            before=(ctx.prec,ctx.rounding,ctx.Emin,ctx.Emax,dict(ctx.flags),dict(ctx.traps))
            self.assertEqual(self.valid('temperature_total_ratio',5e-309,1e308),baseline)
            self.assertEqual((ctx.prec,ctx.rounding,ctx.Emin,ctx.Emax,dict(ctx.flags),dict(ctx.traps)),before)

    def test_default_decimal_environment_cannot_choose_inverse_behavior(self):
        baseline=self.valid('temperature_critical_ratio',1.,GAMMAS[0])
        saved=DefaultContext.copy()
        with localcontext():
            try:
                DefaultContext.prec=5;DefaultContext.rounding=ROUND_DOWN
                DefaultContext.Emax=9;DefaultContext.Emin=-9
                DefaultContext.traps[Inexact]=DefaultContext.traps[FloatOperation]=True
                self.assertEqual(self.valid('temperature_critical_ratio',1.,GAMMAS[0]),baseline)
            finally:
                for key in ('prec','rounding','Emin','Emax','capitals','clamp'):
                    setattr(DefaultContext,key,getattr(saved,key))
                DefaultContext.flags=saved.flags;DefaultContext.traps=saved.traps

    def test_real_http_route_retains_the_inverted_ratio(self):
        for kind,value,gamma in [('temperature_total_ratio',5e-309,1e308),
                                 ('temperature_critical_ratio',1.,GAMMAS[0])]:
            with self.subTest(kind=kind):
                answer=dispatch_calculation('/api/v1/isentropic',{'input_kind':kind,'value':value,'gamma':gamma})
                self.assertEqual(answer['mach'],oracle(kind,value,gamma)[1])
                json.dumps(answer,allow_nan=False)


if __name__=='__main__':unittest.main()
