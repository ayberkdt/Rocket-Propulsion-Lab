"""Independent radial-flow oracles and hostile axis-to-interior contracts."""
import itertools
import unittest
from dataclasses import asdict, replace
from math import inf, ldexp, nextafter, sqrt

import mpmath as mp

from rocket_propulsion.core.errors import ConvergenceError, DomainError
from rocket_propulsion.geometry import (
    MocAxisDepartureResult,
    MocFlowPoint,
    MocWallSegment,
    solve_moc_axis_departure,
    solve_moc_axis_point,
    solve_moc_wall_point,
)
from tests.test_moc_boundary import relative_error
from tests.test_moc_unit_process import exact_mp


def departure_fixture(*,axisymmetric=True,gamma=1.4,mach=2.,step=.05,compressing=False,precision=110):
    """Spherical/planar mass conservation and analytic characteristic angles.

    Input history and outer C- parent are independent analytic states. The
    axis parent itself must be computed by the production boundary process.
    """
    with mp.workdps(precision):
        g,m,h=map(exact_mp,(gamma,mach,step));d=2 if axisymmetric else 1;sign=-1 if compressing else 1
        def area(z):return (2/(g+1)*(1+(g-1)*z*z/2))**((g+1)/(2*(g-1)))/z
        def nu(z):return mp.sqrt((g+1)/(g-1))*mp.atan(mp.sqrt((g-1)/(g+1)*(z*z-1)))-mp.atan(mp.sqrt(z*z-1))
        def point(z,beta):
            radius=2*(area(z)/area(m))**(mp.mpf(1)/d);angle=beta*mp.pi/180
            x=radius*mp.cos(angle)
            return MocFlowPoint(float(x if sign==1 else 4-x),float(radius*mp.sin(angle)),float(z),float(sign*beta))
        parent=m-sign*h;incoming=parent-sign*h
        incoming_turn=abs(nu(parent)-nu(incoming))/d*180/mp.pi
        turn=abs(nu(m)-nu(parent))/d*180/mp.pi
        args={'incoming': point(incoming,incoming_turn),'axis_history': tuple(point(parent-sign*i*h,0) for i in (2,1)),
            'gamma': gamma,'axisymmetric': axisymmetric}
        target=point(m,turn);minus=point(parent,2*turn)
        a2=(g+1)/(2+(g-1)*parent*parent);radius=2*(area(parent)/area(m))**(mp.mpf(1)/d)
        source=float(sign*a2*parent*mp.sqrt(a2)/radius) if axisymmetric else 0.
        return args,minus,(target.x_m,target.radius_m,target.mach,target.flow_angle_deg),source


def coupled_axis_march(count,precision=110):
    """Computed axis -> interior -> axis with exact analytic outer inflow.

    This is a coupled marching stencil, not a complete nozzle calculation.
    All interior/axis states after the initial two history points are computed.
    """
    with mp.workdps(precision):
        g=exact_mp(1.4)
        def nu(z):return mp.sqrt((g+1)/(g-1))*mp.atan(mp.sqrt((g-1)/(g+1)*(z*z-1)))-mp.atan(mp.sqrt(z*z-1))
        def area(z):return (2/(g+1)*(1+(g-1)*z*z/2))**((g+1)/(2*(g-1)))/z
        start=nu(mp.mpf(2));delta=(nu(exact_mp(2.4))-start)/(2*count)
        def point(value,beta):
            m=mp.findroot(lambda z:nu(z)-value,(mp.mpf('1.5'),mp.mpf(3)))
            radius=2*mp.sqrt(area(m)/area(mp.mpf(2)));angle=beta*mp.pi/180
            return MocFlowPoint(float(radius*mp.cos(angle)),float(radius*mp.sin(angle)),float(m),float(beta))
        turn=delta/2*180/mp.pi
        history=[point(start-i*delta,0) for i in (4,2)]
        axis=solve_moc_axis_point(incoming=point(start-delta,turn),axis_history=tuple(history))
        history.append(axis.point)
        for index in range(count):
            minus=point(start+2*index*delta,2*turn)
            off=solve_moc_axis_departure(axis_result=axis,minus_parent=minus)
            axis=solve_moc_axis_point(incoming=off.point,axis_history=tuple(history[-2:]))
            history.append(axis.point)
        target=point(start+2*count*delta,0)
        return axis,tuple(history),max(abs(axis.point.x_m-target.x_m)/target.x_m,abs(axis.point.mach-target.mach)/target.mach)


def uniform_perturbed_march(count,seed=1e-10):
    """Fixed-length constant flow with a small numerical history perturbation."""
    step=1/count;mach=2.;slope=1/sqrt(3.)
    history=[MocFlowPoint(-2*step,0.,mach+seed,0.),MocFlowPoint(-step,0.,mach,0.)]
    axis=solve_moc_axis_point(incoming=MocFlowPoint(-step/2,slope*step/2,mach,0.),axis_history=tuple(history))
    history.append(axis.point)
    for index in range(count):
        off=solve_moc_axis_departure(axis_result=axis,minus_parent=MocFlowPoint(index*step,slope*step,mach,0.))
        axis=solve_moc_axis_point(incoming=off.point,axis_history=tuple(history[-2:]))
        history.append(axis.point)
    return tuple(history),max(abs(p.mach-mach) for p in history)


class MocAxisDepartureAttacks(unittest.TestCase):
    def solve(self,**inputs):
        try:return solve_moc_axis_departure(**inputs)
        except (DomainError,ArithmeticError) as error:self.fail(f'Valid axis-departure witness refused: {error}')

    def case(self,**kw):
        args,minus,target,source=departure_fixture(**kw)
        return {'axis_result': solve_moc_axis_point(**args),'minus_parent': minus},target,source

    def uniform(self,mach=2.):
        axis=solve_moc_axis_point(incoming=MocFlowPoint(0.,1.,mach,0.),
            axis_history=(MocFlowPoint(-2.,0.,mach,0.),MocFlowPoint(-1.,0.,mach,0.)))
        return {'axis_result': axis,'minus_parent': MocFlowPoint(axis.point.x_m,2.,mach,0.)}

    def convergence(self,axisymmetric,compressing):
        for gamma,mach in ((1.4,2.),(1.22,3.)):
            errors=[]
            for h in (.1,.05,.025,.0125,.00625):
                inputs,target,_=self.case(axisymmetric=axisymmetric,compressing=compressing,gamma=gamma,mach=mach,step=h)
                r=self.solve(**inputs);errors.append(relative_error(r.point,target))
                self.assertEqual(r.axis_source_term_per_m,inputs['axis_result'].source_term_per_m)
                self.assertLessEqual(r.rounding_adjusted_residual_max,64*2**-52)
            self.assertTrue(all(6<a/b<12 for a,b in itertools.pairwise(errors)),errors)
            self.assertLess(errors[-1],2e-7)

    def test_spherical_expansion_carries_nonzero_axis_source(self):self.convergence(True,False)
    def test_spherical_compression_carries_negative_axis_source(self):self.convergence(True,True)
    def test_planar_expansion_uses_planar_equations(self):self.convergence(False,False)
    def test_planar_compression_uses_planar_equations(self):self.convergence(False,True)

    def test_independent_fixtures_stable_at_110_and_180_digits(self):
        for mode in (False,True):
            for compress in (False,True):
                self.assertEqual(departure_fixture(axisymmetric=mode,compressing=compress,precision=110),
                                 departure_fixture(axisymmetric=mode,compressing=compress,precision=180))

    def test_uniform_branch_preserves_extreme_mach_and_analytic_intersection(self):
        for mach in (1.01,2.,8.,1e200):
            with self.subTest(mach=mach):
                inputs=self.uniform(mach);r=self.solve(**inputs);axis=inputs['axis_result'].point
                self.assertEqual(r.point.mach,mach);self.assertEqual(r.point.flow_angle_deg,0.)
                self.assertLess(abs(r.point.x_m/(2*axis.x_m)-1),1e-13);self.assertAlmostEqual(r.point.radius_m,1.,places=13)
                self.assertEqual(r.axis_source_term_per_m,0.);self.assertEqual(r.solution_kind,'constant_axis_departure')

    def test_finite_axis_source_has_correct_sign_and_independent_spherical_limit(self):
        for compress in (False,True):
            errors=[]
            for h in (.05,.025,.0125,.00625):
                inputs,_,expected=self.case(compressing=compress,step=h);r=self.solve(**inputs)
                self.assertGreater(r.axis_source_term_per_m*expected,0.)
                errors.append(abs(r.axis_source_term_per_m/expected-1))
            self.assertTrue(all(3<a/b<6 for a,b in itertools.pairwise(errors)),errors)
            self.assertLess(errors[-1],1e-3)

    def test_off_axis_parent_requires_positive_radius(self):
        inputs=self.uniform()
        for radius in (0.,-1.):
            with self.subTest(radius=radius),self.assertRaises(DomainError):
                solve_moc_axis_departure(**{**inputs,'minus_parent':replace(inputs['minus_parent'],radius_m=radius)})

    def test_plain_axis_point_and_wall_result_cannot_supply_missing_regularity_history(self):
        inputs=self.uniform();wall=solve_moc_wall_point(incoming=MocFlowPoint(0.,1.,2.,0.),wall_segment=MocWallSegment(0.,2.,4.,2.))
        for axis in (inputs['axis_result'].point,wall,None,{},asdict(inputs['axis_result'])|{'source_term_per_m':1.}):
            with self.subTest(axis=axis),self.assertRaises(DomainError):solve_moc_axis_departure(**{**inputs,'axis_result':axis})

    def test_private_noncausal_intersection_has_distinct_error(self):
        inputs=self.uniform();p=inputs['minus_parent'];inputs['minus_parent']=replace(p,x_m=p.x_m-10.)
        with self.assertRaises(ConvergenceError) as context:solve_moc_axis_departure(**inputs)
        self.assertEqual(context.exception.code,'moc_departure_noncausal')

    def test_stored_coordinate_cannot_collapse_a_downstream_step(self):
        mach=1.1;scale=1e293
        axis=solve_moc_axis_point(incoming=MocFlowPoint(1e308,scale,mach,0.),
            axis_history=(MocFlowPoint(0.,0.,mach,0.),MocFlowPoint(1.,0.,mach,0.)))
        with self.assertRaises(ConvergenceError) as context:
            solve_moc_axis_departure(axis_result=axis,minus_parent=MocFlowPoint(axis.point.x_m,1.,mach,0.))
        self.assertEqual(context.exception.code,'moc_departure_unresolved_topology')

    def test_iteration_budget_must_fail_instead_of_returning_an_unconverged_state(self):
        inputs,_,_=self.case();self.solve(**inputs)
        with self.assertRaises(ConvergenceError) as context:solve_moc_axis_departure(**inputs,max_iterations=1)
        self.assertEqual(context.exception.code,'moc_departure_not_converged')

    def test_resource_budget_rejects_boolean_float_and_unbounded_counts(self):
        for count in (True,1.,0,-1,257):
            with self.subTest(count=count),self.assertRaises(DomainError):solve_moc_axis_departure(**self.uniform(),max_iterations=count)

    def test_scale_covariance_uses_owned_axis_history_not_an_epsilon_radius(self):
        args,minus,_,_=departure_fixture();base=self.solve(axis_result=solve_moc_axis_point(**args),minus_parent=minus)
        for exponent in (-900,900):
            def point(p, *, exponent=exponent):return replace(p,x_m=ldexp(p.x_m,exponent),radius_m=ldexp(p.radius_m,exponent))
            scaled={**args,'incoming':point(args['incoming']),'axis_history':tuple(point(p) for p in args['axis_history'])}
            r=self.solve(axis_result=solve_moc_axis_point(**scaled),minus_parent=point(minus))
            self.assertEqual((r.point.mach,r.point.flow_angle_deg,r.iterations),(base.point.mach,base.point.flow_angle_deg,base.iterations))
            self.assertEqual((ldexp(r.point.x_m,-exponent),ldexp(r.point.radius_m,-exponent)),(base.point.x_m,base.point.radius_m))

    def test_nested_mapping_round_trip_owns_axis_history_and_off_axis_state(self):
        r=self.solve(**self.case()[0]);raw=asdict(r);restored=MocAxisDepartureResult(**raw);self.assertEqual(r,restored)
        raw['inputs']['axis_result']['inputs']['axis_history'][0]['mach']=100.
        raw['inputs']['minus_parent']['mach']=999.
        self.assertNotEqual(restored.inputs.axis_result.inputs.axis_history[0].mach,100.)
        self.assertNotEqual(restored.inputs.minus_parent.mach,999.)

    def test_result_state_diagnostics_source_and_labels_must_replay(self):
        r=self.solve(**self.case()[0])
        for change in ({'point': replace(r.point,mach=r.point.mach+.01)},{'iterations': r.iterations+1},{'equation_residual_max': 1.},
                       {'rounding_adjusted_residual_max': 1.},{'coordinate_rounding_allowance_max': 1.},{'representation_limited': not r.representation_limited},
                       {'axis_source_term_per_m': 0.},{'axis_source_term_status': 'unrepresentable'},{'solution_kind': 'constant_axis_departure'}):
            with self.subTest(change=change),self.assertRaises(DomainError):replace(r,**change)

    def test_hostile_policy_equality_does_not_bypass_plain_text_ownership(self):
        class EqualText(str):
            def __eq__(self,other):return True
            def __ne__(self,other):return False
        r=self.solve(**self.uniform())
        with self.assertRaises(DomainError):replace(r,numerical_policy=EqualText('forged'))

    def test_local_process_cannot_claim_global_error_or_hardware_qualification(self):
        r=self.solve(**self.uniform());self.assertIsNone(r.discretization_error_bound)
        for change in ({'physical_model_qualified': True},{'full_nozzle_solution': True},{'physical_model_qualified': 0},
                       {'full_nozzle_solution': 0},{'discretization_error_bound': 0.}):
            with self.subTest(change=change),self.assertRaises(DomainError):replace(r,**change)

    def test_computed_axis_interior_axis_march_controls_accumulated_error(self):
        errors=[]
        for count in (4,8,16):
            axis,history,error=coupled_axis_march(count);errors.append(error)
            self.assertTrue(all(b.x_m>a.x_m for a,b in itertools.pairwise(history)))
            self.assertGreater(axis.point.mach,2.)
        self.assertTrue(all(3<a/b<6 for a,b in itertools.pairwise(errors)),errors)
        self.assertLess(errors[-1],2e-4)

    def test_uniform_computed_march_does_not_amplify_small_history_error(self):
        for count in (8,16,32):
            with self.subTest(count=count):
                history,error=uniform_perturbed_march(count)
                self.assertLess(error,1e-9)
                self.assertTrue(all(b.x_m>a.x_m for a,b in itertools.pairwise(history)))

    def test_unenforced_plus_compatibility_defect_is_visible_and_replayed(self):
        r=self.solve(**self.case()[0])
        self.assertGreater(r.averaged_plus_compatibility_defect,64*2**-52)
        self.assertLessEqual(r.rounding_adjusted_residual_max,64*2**-52)
        with self.assertRaises(DomainError):replace(r,averaged_plus_compatibility_defect=0.)
        planar=self.solve(**self.case(axisymmetric=False)[0])
        self.assertLessEqual(planar.averaged_plus_compatibility_defect,64*2**-52)

    def test_unrepresentable_source_diagnostic_does_not_erase_private_regularity(self):
        args,minus,_,_=departure_fixture();base=self.solve(axis_result=solve_moc_axis_point(**args),minus_parent=minus)
        def point(p):return replace(p,x_m=ldexp(p.x_m,-1026),radius_m=ldexp(p.radius_m,-1026))
        args={**args,'incoming':point(args['incoming']),'axis_history':tuple(point(p) for p in args['axis_history'])}
        axis=solve_moc_axis_point(**args);self.assertIsNone(axis.source_term_per_m)
        r=self.solve(axis_result=axis,minus_parent=point(minus))
        self.assertIsNone(r.axis_source_term_per_m);self.assertEqual(r.axis_source_term_status,'unrepresentable')
        self.assertLess(abs(r.point.mach/base.point.mach-1),2e-13)

    def test_one_ulp_parent_difference_cannot_be_called_constant_flow(self):
        inputs=self.uniform();inputs['minus_parent']=replace(inputs['minus_parent'],mach=nextafter(2.,inf))
        r=self.solve(**inputs)
        self.assertEqual(r.solution_kind,'iterated_axis_departure')

    def test_zero_axis_gradient_with_curved_history_is_not_uniform(self):
        # Independently inverted energy states give represented velocities
        # u+4*2**-10 and u+2**-10 at x=-1,0: zero derivative at x=1,
        # but nonzero curvature. A zero source alone is not uniform flow.
        from math import asin, tan
        axis=solve_moc_axis_point(incoming=MocFlowPoint(0.,tan(asin(.5)),2.,0.),
            axis_history=(MocFlowPoint(-1.,0.,2.0086363055731224,0.),MocFlowPoint(0.,0.,2.002154418338075,0.)))
        self.assertEqual(axis.point,MocFlowPoint(1.,0.,2.,0.));self.assertEqual(axis.source_term_per_m,0.)
        r=self.solve(axis_result=axis,minus_parent=MocFlowPoint(1.,.01,2.,0.))
        self.assertEqual(r.solution_kind,'iterated_axis_departure')
        self.assertGreater(r.point.flow_angle_deg,0.)
        self.assertNotEqual(r.point.mach,axis.point.mach)


if __name__=='__main__':unittest.main()
