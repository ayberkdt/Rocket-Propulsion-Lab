"""Adversarial shared-root contracts; sampled callback evidence is not certification."""
import json
import math
import unittest
from dataclasses import asdict, replace
from fractions import Fraction as F
from sys import float_info

import mpmath as mp
from scipy.optimize import brentq

from rocket_propulsion.core.errors import ConvergenceError, DomainError
from rocket_propulsion.core.numerics import RootEvaluation, bisect_monotonic, solve_bisect_monotonic


def decision_boundary(x):
    """Continuous monotonic manufactured function with a rounded-residual trap."""
    e=F(2)**-54;t=F(x)
    return float(-1+2*t*(1-e) if t<=F(1,2) else -e+2*(t-F(1,2))*(3+e))


class RootAttacks(unittest.TestCase):
    def valid(self, function, target, lower, upper, **kwargs):
        try:return solve_bisect_monotonic(function,target,lower,upper,**kwargs)
        except (DomainError,ArithmeticError) as error:self.fail(f'Valid representable root problem refused: {error}')

    def check(self,result):
        self.assertLessEqual(result.lower_bound,result.value)
        self.assertLessEqual(result.value,result.upper_bound)
        by_x={s.x:s.value for s in result.evaluations}
        self.assertEqual(F(result.residual_ratio[0])/F(result.residual_ratio[1]),abs(F(by_x[result.value])-F(result.target)))
        self.assertLessEqual(F(result.residual_ratio[0])/F(result.residual_ratio[1]),F(result.tolerance)*max(F(1),abs(F(result.target))))
        self.assertLessEqual(min(by_x[result.lower_bound],by_x[result.upper_bound]),result.target)
        self.assertLessEqual(result.target,max(by_x[result.lower_bound],by_x[result.upper_bound]))
        self.assertFalse(result.root_enclosure_qualified)
        self.assertFalse(result.physical_model_qualified)
        self.assertIsNone(result.coordinate_error_bound)
        json.dumps(asdict(result),allow_nan=False)

    def test_invalid_controls_fail_before_any_callback(self):
        def forbidden(x):raise AssertionError('Malformed controls reached callback')
        base={'target':.5,'lower':0,'upper':1,'tolerance':1e-12,'max_iterations':200}
        for name in ('target','lower','upper','tolerance'):
            for value in (True,None,'1',float('nan'),float('inf'),-float('inf')):
                args=dict(base);args[name]=value
                with self.subTest(field=name,value=value),self.assertRaises(DomainError):solve_bisect_monotonic(forbidden,**args)
        for value in (True,None,1.5,'200',0,-1,4097):
            args=dict(base);args['max_iterations']=value
            with self.subTest(budget=value),self.assertRaises(DomainError):solve_bisect_monotonic(forbidden,**args)
        for value in (0,-1):
            args=dict(base);args['tolerance']=value
            with self.subTest(tolerance=value),self.assertRaises(DomainError):solve_bisect_monotonic(forbidden,**args)

    def test_noncallable_and_reversed_bounds_are_typed(self):
        for function in (None,True,1,{}):
            with self.subTest(function=function),self.assertRaises(DomainError):solve_bisect_monotonic(function,.5,0,1)
        with self.assertRaises(DomainError):solve_bisect_monotonic(lambda x:x,.5,1,0)

    def test_invalid_callback_values_stop_at_first_observation(self):
        for value in (True,None,'1',float('nan'),float('inf'),-float('inf')):
            calls=[]
            def invalid(x, *, value=value,calls=calls):calls.append(x);return value
            with self.subTest(value=value),self.assertRaises(DomainError) as error:
                solve_bisect_monotonic(invalid,.5,0,1)
            self.assertEqual(calls,[0.0])
            json.dumps(error.exception.as_dict(),allow_nan=False)

    def test_invalid_interior_observation_cannot_be_skipped(self):
        calls=[]
        def function(x):
            calls.append(x)
            return float('nan') if x==.5 else x
        with self.assertRaises(DomainError):solve_bisect_monotonic(function,.5,0,1)
        self.assertEqual(calls,[0.0,1.0,.5])

    def test_callback_application_exceptions_propagate_unchanged(self):
        failure=RuntimeError('External evaluator failed')
        def broken(x):raise failure
        with self.assertRaises(RuntimeError) as error:solve_bisect_monotonic(broken,.5,0,1)
        self.assertIs(error.exception,failure)

    def test_target_outside_actual_bracket_is_not_admitted_by_tolerance(self):
        for target in (1+5e-13,-5e-13):
            for sign in (1,-1):
                with self.subTest(target=target,sign=sign):
                    with self.assertRaises(DomainError) as error:
                        solve_bisect_monotonic(lambda x, sign=sign:sign*x,sign*target,0,1)
                    self.assertEqual(error.exception.code,'target_not_bracketed')

    def test_increasing_and_decreasing_exact_affine_oracles(self):
        for sign in (1,-1):
            for target in (.1,.3,1,1.7):
                function=lambda x, sign=sign:sign*float(F(3)*F(x))
                result=self.valid(function,sign*target,0,1,tolerance=1e-14)
                exact=F(target)/3
                self.assertLessEqual(F(result.lower_bound),exact)
                self.assertLessEqual(exact,F(result.upper_bound))
                self.check(result)

    def test_huge_same_sign_interval_does_not_overflow_its_midpoint(self):
        result=self.valid(lambda x:x,1.05e308,1e308,1.1e308)
        self.assertTrue(math.isfinite(result.value));self.check(result)

    def test_opposite_finite_extremes_and_large_residual_scale_stay_json_safe(self):
        for target,tolerance in ((0,1e-12),(float_info.max/2,float_info.max)):
            result=self.valid(lambda x:x,target,-float_info.max,float_info.max,tolerance=tolerance)
            self.check(result)

    def test_subnormal_root_and_tolerance_are_not_erased(self):
        tiny=math.ulp(0.0)
        result=self.valid(lambda x:2*x,2*tiny,0,4*tiny,tolerance=tiny)
        self.assertEqual(result.value,tiny);self.assertEqual(result.residual,0);self.check(result)

    def test_exact_endpoint_and_singleton_need_no_interior_evaluation(self):
        for lower,upper,target in ((0,1,0),(0,1,1),(1,1,1)):
            calls=[]
            def function(x, *, calls=calls):calls.append(x);return x
            result=self.valid(function,target,lower,upper)
            self.assertEqual(result.value,target);self.assertEqual(result.iterations,0)
            self.assertEqual(len(calls),1 if lower==upper else 2);self.check(result)
        with self.assertRaises(DomainError):solve_bisect_monotonic(lambda x:x,2,1,1)

    def test_flat_matching_callback_does_not_assert_uniqueness_or_qualification(self):
        result=self.valid(lambda x:2,2,-1,1)
        self.assertEqual(result.residual,0);self.check(result)

    def test_exact_residual_decision_is_not_made_from_rounded_subtraction(self):
        result=self.valid(decision_boundary,1,0,1,tolerance=1)
        self.assertNotEqual(result.value,.5)
        self.assertLessEqual(abs(F(decision_boundary(result.value))-1),F(1));self.check(result)

    def test_budget_is_actual_integer_and_failure_counts_real_observations(self):
        calls=[]
        def function(x):calls.append(x);return x*x
        with self.assertRaises(ConvergenceError) as error:
            solve_bisect_monotonic(function,2,0,2,tolerance=1e-30,max_iterations=1)
        self.assertEqual(calls,[0.0,2.0,1.0]);self.assertEqual(error.exception.details['iterations'],1)
        self.assertEqual(error.exception.code,'non_convergence')
        self.assertEqual(error.exception.details['termination'],'non_convergence')
        self.assertIs(error.exception.details['acceptance_satisfied'],False)
        json.dumps(error.exception.as_dict(),allow_nan=False)

    def test_adjacent_float_resolution_failure_preserves_signs_without_relaxing_tolerance(self):
        lower=1.0;upper=math.nextafter(1,math.inf);target=math.nextafter(3,math.inf)
        with self.assertRaises(ConvergenceError) as error:
            solve_bisect_monotonic(lambda x:3*x,target,lower,upper,tolerance=math.ulp(0.0))
        self.assertEqual(error.exception.code,'root_float_resolution')
        self.assertEqual(error.exception.details['termination'],'root_float_resolution')
        self.assertIs(error.exception.details['acceptance_satisfied'],False)
        details=error.exception.details
        self.assertEqual(details['lower_bound'],lower);self.assertEqual(details['upper_bound'],upper)
        self.assertLess(details['lower_function_value'],target);self.assertGreater(details['upper_function_value'],target)
        self.assertGreater(details['residual'],math.ulp(0.0));json.dumps(error.exception.as_dict(),allow_nan=False)

    def test_monotonicity_violating_observation_is_refused(self):
        def function(x):return 2 if x==.5 else x
        with self.assertRaises(DomainError) as error:solve_bisect_monotonic(function,.5,0,1)
        self.assertEqual(error.exception.code,'callback_not_monotonic')

    def test_discontinuous_step_does_not_create_a_converged_root(self):
        with self.assertRaises(ConvergenceError) as error:
            solve_bisect_monotonic(lambda x:0 if x<0 else 1,.5,-1,1,tolerance=1e-12,max_iterations=4096)
        self.assertEqual(error.exception.code,'root_float_resolution')
        json.dumps(error.exception.as_dict(),allow_nan=False)

    def test_result_replay_does_not_call_the_external_function_again(self):
        calls=[];observations=[]
        def function(x):calls.append(x);observations.append((x,x*x));return x*x
        result=self.valid(function,2,0,2)
        count=len(calls);copy=replace(result)
        self.assertEqual(len(calls),count);self.assertEqual(copy,result)
        self.assertEqual(len(result.evaluations),result.iterations+2)
        self.assertEqual([(row.x,row.value) for row in result.evaluations],observations)
        self.assertIsNot(copy.evaluations[0],result.evaluations[0]);self.check(result)

    def test_derived_fields_and_trace_coordinates_cannot_be_forged(self):
        result=self.valid(lambda x:x*x,2,0,2)
        for name in ('value','residual','lower_bound','upper_bound','scaled_residual'):
            with self.subTest(field=name),self.assertRaises(DomainError):
                replace(result,**{name:math.nextafter(getattr(result,name),math.inf)})
        with self.assertRaises(DomainError):replace(result,iterations=result.iterations+1)
        altered=replace(result.evaluations[-1],x=math.nextafter(result.evaluations[-1].x,math.inf))
        with self.assertRaises(DomainError):replace(result,evaluations=(*result.evaluations[:-1],altered))
        with self.assertRaises(DomainError):replace(result,residual_ratio=('0','1'))

    def test_false_evidence_flags_and_coordinate_bound_cannot_be_promoted(self):
        result=self.valid(lambda x:x*x,2,0,2)
        for name in ('root_enclosure_qualified','physical_model_qualified'):
            for value in (True,0,None):
                with self.subTest(field=name,value=value),self.assertRaises(DomainError):replace(result,**{name:value})
        with self.assertRaises(DomainError):replace(result,coordinate_error_bound=0.0)
        with self.assertRaises(DomainError):replace(result,method='certified')
        with self.assertRaises(DomainError):replace(result,evidence_kind='certified')

    def test_small_slope_residual_acceptance_is_not_reported_as_x_accuracy(self):
        result=self.valid(lambda x:1e-300*x,1.25e-300,0,2)
        self.assertLessEqual(result.residual,result.tolerance)
        # Even x=0 satisfies this residual rule while being far from x=1.25.
        self.assertLessEqual(1.25e-300,result.tolerance)
        self.assertIsNone(result.coordinate_error_bound);self.assertFalse(result.root_enclosure_qualified)
        self.check(result)

    def test_nonlinear_roots_agree_with_high_precision_and_independent_brent_solver(self):
        cases=((lambda x:x*x,2,lambda x:x*x),(lambda x:x*x*x+2*x,7,lambda x:x*x*x+2*x),
               (math.exp,10,mp.exp))
        for function,target,oracle in cases:
            result=self.valid(function,target,0,4,tolerance=1e-13)
            with mp.workdps(180):exact=mp.findroot(lambda x,oracle=oracle,target=target:oracle(x)-target,1)
            independent=brentq(lambda x, function=function,target=target:function(x)-target,0,4,xtol=1e-14)
            self.assertAlmostEqual(result.value,float(exact),delta=2e-12)
            self.assertAlmostEqual(result.value,independent,delta=2e-12);self.check(result)

    def test_float_wrapper_keeps_the_same_verified_algorithm(self):
        result=self.valid(lambda x:x*x,2,0,2,tolerance=1e-13)
        self.assertEqual(bisect_monotonic(lambda x:x*x,2,0,2,tolerance=1e-13),result.value)
        with self.assertRaises(DomainError):bisect_monotonic(lambda x:x,True,0,1)

    def test_owned_trace_is_bounded_and_malformed_samples_are_rejected(self):
        result=self.valid(lambda x:x*x,2,0,2)
        for value in (None,[],{},'trace',set(result.evaluations),(object(),)):
            with self.subTest(trace=repr(value)[:30]),self.assertRaises(DomainError):replace(result,evaluations=value)
        with self.assertRaises(DomainError):replace(result,evaluations=(*result.evaluations,result.evaluations[-1]))

    def test_public_evaluation_records_refuse_nonfinite_and_boolean_fields(self):
        for field in ('x','value'):
            for bad in (True,None,'1',float('nan'),float('inf')):
                arguments={'x':0,'value':1};arguments[field]=bad
                with self.subTest(field=field,value=bad),self.assertRaises(DomainError):RootEvaluation(**arguments)
        result=RootEvaluation(0,1)
        self.assertIs(type(result.x),float);self.assertIs(type(result.value),float)

