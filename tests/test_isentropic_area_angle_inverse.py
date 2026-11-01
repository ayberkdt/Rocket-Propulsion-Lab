"""Independent exact-input area/angle inverses, not rounded forward round trips."""
import json
import math
import unittest
from decimal import ROUND_DOWN, DefaultContext, FloatOperation, Inexact, localcontext
from functools import cache
from unittest.mock import patch

import mpmath as mp

from rocket_propulsion.api.routes import dispatch_calculation
from rocket_propulsion.compressible.isentropic import mach_from_isentropic
from rocket_propulsion.core.errors import ConvergenceError, DomainError

NEXT_GAMMA=math.nextafter(1.,math.inf)


@cache
def area_oracle(ratio,gamma,subsonic):
    """Entropy and sonic mass-flux conservation, with independent bisection."""
    with mp.workdps(600):
        g,wanted=mp.mpf(gamma),mp.log(mp.mpf(ratio))
        lower,upper=(mp.mpf(1e-8),mp.mpf(1)) if subsonic else (mp.mpf(1),mp.mpf(1000))
        def log_area(m):
            temperature=(g+1)/(2+(g-1)*m*m)
            log_density=mp.log(temperature)/(g-1)
            return -log_density-mp.log(m)-mp.log(temperature)/2
        left,right=log_area(lower),log_area(upper)
        if not min(left,right)<=wanted<=max(left,right):return 'outside_inverse_range',None
        for _ in range(140):
            middle=(lower+upper)/2
            if (log_area(middle)>wanted)==(not subsonic):upper=middle
            else:lower=middle
        return 'valid',float((lower+upper)/2)


def angle_oracle(value):
    with mp.workdps(1200):
        root=1/mp.sin(mp.mpf(value)*mp.pi/180)
        if not mp.mpf(1e-8)<=root<=1000:return 'outside_inverse_range',None
        return 'valid',float(root)


def area_boundary_ratio(gamma,mach):
    with mp.workdps(600):
        g,m=mp.mpf(gamma),mp.mpf(mach)
        t=(g+1)/(2+(g-1)*m*m)
        return float(1/(mp.exp(mp.log(t)/(g-1))*m*mp.sqrt(t)))


class IsentropicAreaAngleInverseTests(unittest.TestCase):
    def valid(self,kind,value,gamma):
        try:return mach_from_isentropic(kind,value,gamma)
        except (DomainError,ConvergenceError,ArithmeticError) as error:
            self.fail(f'Admitted inverse refused: {kind}, {value}, gamma={gamma}: {error!r}')

    def compare(self,kind,value,gamma=1.4):
        status,expected=(angle_oracle(value) if kind=='mach_angle_deg'
                         else area_oracle(value,gamma,kind.endswith('subsonic')))
        if status!='valid':
            try:
                with self.assertRaises(DomainError) as caught:mach_from_isentropic(kind,value,gamma)
            except ArithmeticError as error:self.fail(f'Expected typed inverse refusal, got {error!r}')
            self.assertEqual(caught.exception.code,status)
            self.assertEqual(caught.exception.field,'value')
            return
        actual=self.valid(kind,value,gamma)
        self.assertTrue(math.isfinite(actual) and 1e-8<=actual<=1000.)
        self.assertLessEqual(abs(actual-expected),2*math.ulp(expected))
        json.dumps({'mach':actual},allow_nan=False)

    def test_valid_area_does_not_require_a_representable_endpoint_area(self):
        for gamma in (NEXT_GAMMA,1.001,1.01):
            with self.subTest(gamma=gamma):self.compare('area_ratio_supersonic',2.,gamma)

    def test_tiny_positive_mach_angles_have_typed_range_refusals(self):
        for value in (math.ulp(0.),math.nextafter(math.ulp(0.),math.inf),1e-320,1e-100,.01):
            with self.subTest(value=value):self.compare('mach_angle_deg',value)

    def test_both_area_branches_match_independent_entropy_and_continuity(self):
        for gamma in (NEXT_GAMMA,1.001,1.01,1.4,4.,10.,1e3,1e308):
            for ratio in (math.nextafter(1.,math.inf),1.001,1.01,2.,20.,1e308):
                for subsonic in (True,False):
                    with self.subTest(gamma=gamma,ratio=ratio,subsonic=subsonic):
                        self.compare('area_ratio_subsonic' if subsonic else 'area_ratio_supersonic',ratio,gamma)

    def test_area_boundary_neighbors_classify_before_float_export(self):
        for gamma in (NEXT_GAMMA,1.01,1.4,4.,1e3,1e308):
            for subsonic,mach in ((True,1e-8),(False,1000.)):
                ratio=area_boundary_ratio(gamma,mach)
                if not math.isfinite(ratio) or ratio<=1:continue  # No distinct admitted binary64 ratio.
                for neighbor in (math.nextafter(ratio,1.),ratio,math.nextafter(ratio,math.inf)):
                    if not math.isfinite(neighbor) or neighbor<1:continue
                    with self.subTest(gamma=gamma,value=neighbor,subsonic=subsonic):
                        self.compare('area_ratio_subsonic' if subsonic else 'area_ratio_supersonic',neighbor,gamma)

    def test_angle_degree_conversion_and_gamma_independence(self):
        for gamma in (NEXT_GAMMA,1.4,1e308):
            for value in (.1,.5,1.,10.,30.,45.,60.,89.,89.99,math.nextafter(90.,0.),90.):
                with self.subTest(gamma=gamma,value=value):self.compare('mach_angle_deg',value,gamma)

    def test_angle_inverse_boundary_neighbors(self):
        with mp.workdps(1200):boundary=float(180/mp.pi*mp.asin(mp.mpf(1)/1000))
        for value in (math.nextafter(boundary,0.),boundary,math.nextafter(boundary,math.inf)):
            with self.subTest(value=value):self.compare('mach_angle_deg',value)

    def test_area_inverse_does_not_materialize_binary64_endpoint_states(self):
        with patch('rocket_propulsion.compressible.isentropic._area_ratio',
                   side_effect=AssertionError('A private endpoint area is not an output requirement')):
            for subsonic in (False,True):
                self.compare('area_ratio_subsonic' if subsonic else 'area_ratio_supersonic',2.,1.01)

    def test_active_decimal_context_is_owned_and_unchanged(self):
        cases=[('area_ratio_supersonic',2.,1.01),('mach_angle_deg',30.,1.4)]
        baseline=[self.valid(*args) for args in cases]
        with localcontext() as ctx:
            ctx.prec=5;ctx.rounding=ROUND_DOWN;ctx.Emin=-9;ctx.Emax=9
            ctx.traps[FloatOperation]=ctx.traps[Inexact]=True
            before=(ctx.prec,ctx.rounding,ctx.Emin,ctx.Emax,dict(ctx.flags),dict(ctx.traps))
            self.assertEqual([self.valid(*args) for args in cases],baseline)
            self.assertEqual((ctx.prec,ctx.rounding,ctx.Emin,ctx.Emax,dict(ctx.flags),dict(ctx.traps)),before)

    def test_default_decimal_context_cannot_choose_solver_behavior(self):
        cases=[('area_ratio_supersonic',2.,1.01),('mach_angle_deg',30.,1.4)]
        baseline=[self.valid(*args) for args in cases]
        saved=DefaultContext.copy()
        with localcontext():
            try:
                DefaultContext.prec=5;DefaultContext.rounding=ROUND_DOWN
                DefaultContext.Emax=9;DefaultContext.Emin=-9
                DefaultContext.traps[FloatOperation]=DefaultContext.traps[Inexact]=True
                self.assertEqual([self.valid(*args) for args in cases],baseline)
            finally:
                for key in ('prec','rounding','Emin','Emax','capitals','clamp'):setattr(DefaultContext,key,getattr(saved,key))
                DefaultContext.flags=saved.flags;DefaultContext.traps=saved.traps

    def test_area_budget_exhaustion_is_not_a_last_iterate(self):
        with (patch('rocket_propulsion.compressible.sonic_numerics._AREA_INVERSE_ITERATIONS',0),
              self.assertRaises(ConvergenceError) as caught):
            mach_from_isentropic('area_ratio_supersonic',2.,1.4)
        self.assertEqual(caught.exception.code,'area_inverse_not_converged')
        self.assertEqual(caught.exception.field,'value')

    def test_angle_series_exhaustion_is_not_a_partial_sine(self):
        with (patch('rocket_propulsion.compressible.isentropic_inverse._ANGLE_SERIES_ITERATIONS',1),
              self.assertRaises(ConvergenceError) as caught):
            mach_from_isentropic('mach_angle_deg',30.,1.4)
        self.assertEqual(caught.exception.code,'mach_angle_not_converged')
        self.assertEqual(caught.exception.field,'value')

    def test_pi_failure_is_propagated_as_a_typed_inverse_failure(self):
        original=ConvergenceError('Injected arctangent exhaustion',field='prandtl_meyer_deg')
        with (patch('rocket_propulsion.compressible.isentropic_inverse._atan',side_effect=original),
              self.assertRaises(ConvergenceError) as caught):
            mach_from_isentropic('mach_angle_deg',30.,1.4)
        self.assertEqual(caught.exception.code,'mach_angle_not_converged')
        self.assertEqual(caught.exception.field,'value')
        self.assertIs(caught.exception.__cause__,original)

    def test_typed_refusal_restores_the_callers_decimal_state(self):
        with localcontext() as ctx:
            ctx.prec=5;ctx.rounding=ROUND_DOWN;ctx.Emin=-9;ctx.Emax=9
            ctx.traps[FloatOperation]=ctx.traps[Inexact]=True
            before=(ctx.prec,ctx.rounding,ctx.Emin,ctx.Emax,dict(ctx.flags),dict(ctx.traps))
            for args in [('area_ratio_supersonic',1e308,1.4),('mach_angle_deg',math.ulp(0.),1.4)]:
                try:
                    with self.assertRaises(DomainError):mach_from_isentropic(*args)
                except ArithmeticError as error:self.fail(f'Uncontrolled private arithmetic: {error!r}')
                self.assertEqual((ctx.prec,ctx.rounding,ctx.Emin,ctx.Emax,dict(ctx.flags),dict(ctx.traps)),before)

    def test_physical_inverse_admission_has_typed_value_errors(self):
        for kind,value in [('area_ratio_subsonic',math.nextafter(1.,0.)),
                           ('area_ratio_supersonic',0.),('mach_angle_deg',0.),
                           ('mach_angle_deg',math.nextafter(90.,math.inf))]:
            with self.subTest(kind=kind,value=value):
                try:
                    with self.assertRaises(DomainError) as caught:mach_from_isentropic(kind,value,1.4)
                except ArithmeticError as error:self.fail(f'Uncontrolled physical refusal: {error!r}')
                self.assertEqual(caught.exception.code,'physical_domain')
                self.assertEqual(caught.exception.field,'value')

    def test_real_route_keeps_typed_angle_refusal_and_valid_area(self):
        answer=dispatch_calculation('/api/v1/isentropic',{'input_kind':'area_ratio_supersonic','value':2.,'gamma':1.01})
        self.assertLessEqual(abs(answer['mach']-area_oracle(2.,1.01,False)[1]),2*math.ulp(answer['mach']))
        json.dumps(answer,allow_nan=False)
        try:
            with self.assertRaises(DomainError) as caught:
                dispatch_calculation('/api/v1/isentropic',{'input_kind':'mach_angle_deg','value':math.ulp(0.),'gamma':1.4})
        except ArithmeticError as error:self.fail(f'Uncontrolled route arithmetic: {error!r}')
        self.assertEqual(caught.exception.code,'outside_inverse_range')


if __name__=='__main__':unittest.main()
