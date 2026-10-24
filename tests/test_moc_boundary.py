"""Analytic axis regularity, radial slip, marching histories and hostile records."""
import itertools
import unittest
from dataclasses import asdict, replace
from math import ldexp, sqrt

import mpmath as mp

from rocket_propulsion.core.errors import ConvergenceError, DomainError
from rocket_propulsion.geometry import (
    MocFlowPoint,
    MocWallSegment,
    solve_moc_axis_point,
    solve_moc_wall_point,
)
from tests.test_moc_unit_process import exact_mp


def radial_point(mach,angle,gamma,dimension,reference_mach=2.,precision=110):
    """Mass conservation rho*V*R**dimension and the analytic radial gas state."""
    with mp.workdps(precision):
        g,m,mref=map(exact_mp,(gamma,mach,reference_mach));theta=exact_mp(angle)*mp.pi/180
        def area(z):return (2/(g+1)*(1+(g-1)*z*z/2))**((g+1)/(2*(g-1)))/z
        radius=2*(area(m)/area(mref))**(mp.mpf(1)/dimension)
        return MocFlowPoint(float(radius*mp.cos(theta)),float(radius*mp.sin(theta)),mach,angle)


def fixture(kind,*,axisymmetric=True,step=.05,mach=2.,gamma=1.4,precision=110):
    with mp.workdps(precision):
        d=2 if axisymmetric else 1;g=exact_mp(gamma);m=exact_mp(mach);h=exact_mp(step)
        def nu(z):return mp.sqrt((g+1)/(g-1))*mp.atan(mp.sqrt((g-1)/(g+1)*(z*z-1)))-mp.atan(mp.sqrt(z*z-1))
        turn=(nu(m)-nu(m-h))/d*180/mp.pi
        a2=(g+1)/(2+(g-1)*m*m);source=float(a2*m*mp.sqrt(a2)/2) if axisymmetric else 0.
        if kind=='axis':
            history=tuple(radial_point(float(m-i*h),0.,gamma,d,mach,precision) for i in (2,1))
            incoming=radial_point(float(m-h),float(turn),gamma,d,mach,precision)
            return {'incoming': incoming,'axis_history': history,'gamma': gamma,'axisymmetric': axisymmetric},(2.,0.,mach,0.),source
        angle=20.;incoming=radial_point(float(m-h),float(exact_mp(angle)-turn),gamma,d,mach,precision)
        a=radial_point(float(m-2*h),angle,gamma,d,mach,precision);b=radial_point(float(m+h),angle,gamma,d,mach,precision)
        wall=MocWallSegment(a.x_m,a.radius_m,b.x_m,b.radius_m)
        target=radial_point(mach,angle,gamma,d,mach,precision)
        return {'incoming': incoming,'wall_segment': wall,'gamma': gamma,'axisymmetric': axisymmetric},(target.x_m,target.radius_m,mach,angle),source


def relative_error(point,target):
    return max(abs(point.x_m-target[0])/2,abs(point.radius_m-target[1])/2,
               abs(point.mach-target[2])/target[2],abs(point.flow_angle_deg-target[3])/180)


def sink_fixture(kind,step=.025):
    """The inward radial solution, reflected about x=4, has a negative source."""
    with mp.workdps(110):
        g=exact_mp(1.4);m=mp.mpf(2);h=exact_mp(step)
        def area(z):return (2/(g+1)*(1+(g-1)*z*z/2))**((g+1)/(2*(g-1)))/z
        def nu(z):return mp.sqrt((g+1)/(g-1))*mp.atan(mp.sqrt((g-1)/(g+1)*(z*z-1)))-mp.atan(mp.sqrt(z*z-1))
        def point(z,beta):
            radius=2*mp.sqrt(area(z)/area(m));angle=beta*mp.pi/180
            return MocFlowPoint(float(4-radius*mp.cos(angle)),float(radius*mp.sin(angle)),float(z),float(-beta))
        delta=(nu(m+h)-nu(m))/2*180/mp.pi
        if kind=='axis':
            return {'incoming': point(m+h,delta),'axis_history': (point(m+2*h,0),point(m+h,0))}, (2.,0.,2.,0.)
        angle=mp.mpf(20);a,b=point(m+2*h,angle),point(m-h,angle);target=point(m,angle)
        return ({'incoming': point(m+h,angle-delta),'wall_segment': MocWallSegment(a.x_m,a.radius_m,b.x_m,b.radius_m)},
                (target.x_m,target.radius_m,2.,-20.))


class MocBoundaryAttacks(unittest.TestCase):
    def solve(self,kind,**inputs):
        try:return (solve_moc_axis_point if kind=='axis' else solve_moc_wall_point)(**inputs)
        except (DomainError,ArithmeticError) as error:self.fail(f'Valid {kind} boundary witness refused: {error}')

    def uniform(self,kind,mach=2.):
        incoming=MocFlowPoint(0.,1.,mach,0.)
        if kind=='axis':return {'incoming': incoming,'axis_history': (MocFlowPoint(-2.,0.,mach,0.),MocFlowPoint(-1.,0.,mach,0.))}
        return {'incoming': incoming,'wall_segment': MocWallSegment(0.,2.,2*mach,2.)}

    def convergence(self,kind,axisymmetric):
        for gamma,mach in ((1.4,2.),(1.22,3.)):
            errors=[]
            for step in (.1,.05,.025,.0125,.00625,.003125):
                inputs,target,_=fixture(kind,axisymmetric=axisymmetric,step=step,mach=mach,gamma=gamma)
                result=self.solve(kind,**inputs);errors.append(relative_error(result.point,target))
                self.assertIs(result.physical_model_qualified,False);self.assertIs(result.full_nozzle_solution,False)
            self.assertTrue(all(6<a/b<12 for a,b in itertools.pairwise(errors)),errors)
            self.assertLess(errors[-1],2e-7)

    def test_axisymmetric_axis_reflection_has_third_order_local_error(self):self.convergence('axis',True)
    def test_planar_axis_reflection_has_third_order_local_error(self):self.convergence('axis',False)
    def test_axisymmetric_radial_wall_slip_has_third_order_local_error(self):self.convergence('wall',True)
    def test_planar_radial_wall_slip_has_third_order_local_error(self):self.convergence('wall',False)

    def test_fixture_precision_is_stable_at_110_and_180_digits(self):
        for kind in ('axis','wall'):
            for mode in (False,True):self.assertEqual(fixture(kind,axisymmetric=mode,precision=110),fixture(kind,axisymmetric=mode,precision=180))

    def sink_convergence(self,kind):
        errors=[]
        for step in (.1,.05,.025,.0125,.00625):
            inputs,target=sink_fixture(kind,step);r=self.solve(kind,**inputs)
            self.assertLess(r.source_term_per_m,0.)
            errors.append(relative_error(r.point,target))
        self.assertTrue(all(6<a/b<12 for a,b in itertools.pairwise(errors)),errors)
        self.assertLess(errors[-1],2e-7)

    def test_axis_compression_retains_negative_source_and_refines(self):self.sink_convergence('axis')
    def test_converging_wall_retains_negative_source_and_refines(self):self.sink_convergence('wall')

    def test_axis_source_limit_is_nonzero_and_converges_to_spherical_value(self):
        errors=[]
        for step in (.05,.025,.0125,.00625):
            inputs,_,source=fixture('axis',step=step);r=self.solve('axis',**inputs)
            self.assertGreater(r.source_term_per_m,0.)
            errors.append(abs(r.source_term_per_m/source-1))
        self.assertTrue(all(3<a/b<6 for a,b in itertools.pairwise(errors)),errors)
        self.assertLess(errors[-1],1e-3)

    def test_uniform_boundary_states_preserve_mach_and_closed_form_position(self):
        for kind in ('axis','wall'):
            for mach in (1.01,2.,8.,1e200):
                with self.subTest(kind=kind,mach=mach):
                    r=self.solve(kind,**self.uniform(kind,mach))
                    expected=mach if mach==1e200 else sqrt(mach*mach-1)
                    self.assertLess(abs(r.point.x_m/expected-1),1e-13)
                    self.assertEqual(r.point.mach,mach)
                    self.assertEqual(r.point.radius_m,0. if kind=='axis' else 2.)
                    self.assertEqual(r.source_term_per_m,0.)

    def test_axis_history_must_be_two_ordered_exact_symmetry_states(self):
        x=self.uniform('axis');a,b=x['axis_history']
        for history in ((),(a,),[b,a],(a,b,b),(a,replace(b,radius_m=1e-30)),(a,replace(b,flow_angle_deg=1e-30)),None,set()):
            with self.subTest(history=history),self.assertRaises(DomainError):solve_moc_axis_point(**{**x,'axis_history':history})

    def test_incoming_axis_or_wall_point_cannot_be_on_or_below_axis(self):
        for kind in ('axis','wall'):
            x=self.uniform(kind)
            for radius in (0.,-1.):
                with self.subTest(kind=kind,radius=radius),self.assertRaises(DomainError):
                    (solve_moc_axis_point if kind=='axis' else solve_moc_wall_point)(**{**x,'incoming':replace(x['incoming'],radius_m=radius)})

    def test_wall_segment_rejects_reversed_x_nonpositive_radius_and_boolean_scalars(self):
        w=self.uniform('wall')['wall_segment']
        for change in ({'x1_m': 0.},{'x0_m': 9.},{'radius0_m': 0.},{'radius1_m': -1.},{'x0_m': True},{'radius0_m': '2'}):
            with self.subTest(change=change),self.assertRaises(DomainError):replace(w,**change)

    def test_finite_wall_segment_is_not_extrapolated(self):
        x=self.uniform('wall');self.solve('wall',**x)
        x['wall_segment']=replace(x['wall_segment'],x1_m=1.)
        with self.assertRaises(ConvergenceError) as context:solve_moc_wall_point(**x)
        self.assertEqual(context.exception.code,'moc_wall_outside_segment')

    def test_causality_is_preserved_after_coordinate_rounding(self):
        for kind in ('axis','wall'):
            x=self.uniform(kind,1.1);self.solve(kind,**x);x['incoming']=replace(x['incoming'],x_m=1e308)
            if kind=='wall':x['wall_segment']=MocWallSegment(0.,2.,1.5e308,2.)
            with self.subTest(kind=kind):
                with self.assertRaises(ConvergenceError) as context:
                    (solve_moc_axis_point if kind=='axis' else solve_moc_wall_point)(**x)
                self.assertEqual(context.exception.code,'moc_'+kind+'_unresolved_topology')

    def test_iteration_budget_is_not_silently_ignored(self):
        for kind in ('axis','wall'):
            inputs,_,_=fixture(kind);self.solve(kind,**inputs)
            with self.subTest(kind=kind),self.assertRaises(ConvergenceError):
                (solve_moc_axis_point if kind=='axis' else solve_moc_wall_point)(**inputs,max_iterations=1)

    def test_scaling_preserves_flow_and_iteration_count(self):
        for kind in ('axis','wall'):
            inputs,_,_=fixture(kind);base=self.solve(kind,**inputs)
            for exponent in (-900,900):
                x=dict(inputs);p=x['incoming'];x['incoming']=replace(p,x_m=ldexp(p.x_m,exponent),radius_m=ldexp(p.radius_m,exponent))
                if kind=='axis':x['axis_history']=tuple(replace(p,x_m=ldexp(p.x_m,exponent)) for p in inputs['axis_history'])
                else:
                    w=x['wall_segment'];x['wall_segment']=MocWallSegment(*(ldexp(v,exponent) for v in (w.x0_m,w.radius0_m,w.x1_m,w.radius1_m)))
                r=self.solve(kind,**x)
                self.assertEqual((r.point.mach,r.point.flow_angle_deg,r.iterations),(base.point.mach,base.point.flow_angle_deg,base.iterations))
                self.assertEqual(ldexp(r.point.x_m,-exponent),base.point.x_m)
                self.assertEqual(ldexp(r.point.radius_m,-exponent),base.point.radius_m)

    def test_inputs_and_nested_result_records_reopen_without_aliases(self):
        for kind in ('axis','wall'):
            r=self.solve(kind,**fixture(kind)[0]);raw=asdict(r);restored=type(r)(**raw)
            self.assertEqual(restored,r);raw['inputs']['incoming']['mach']=999.
            self.assertNotEqual(restored.inputs.incoming.mach,999.)

    def test_retained_solution_source_diagnostics_and_policy_replay(self):
        for kind in ('axis','wall'):
            r=self.solve(kind,**fixture(kind)[0])
            for change in ({'point': replace(r.point,mach=r.point.mach+.1)},{'iterations': r.iterations+1},{'source_term_per_m': 0.},
                           {'source_term_status': 'unrepresentable'},{'solution_kind': 'forged'},{'numerical_policy': 'forged'}):
                with self.subTest(kind=kind,change=change),self.assertRaises(DomainError):replace(r,**change)

    def test_boundary_process_cannot_claim_a_nozzle_error_bound_or_hardware_qualification(self):
        for kind in ('axis','wall'):
            r=self.solve(kind,**self.uniform(kind));self.assertIsNone(r.discretization_error_bound)
            for change in ({'physical_model_qualified': True},{'full_nozzle_solution': True},{'physical_model_qualified': 0},{'discretization_error_bound': 0.}):
                with self.subTest(kind=kind,change=change),self.assertRaises(DomainError):replace(r,**change)

    def test_hostile_string_equality_cannot_hide_forged_boundary_policy(self):
        class EqualText(str):
            def __eq__(self,other):return True
            def __ne__(self,other):return False
        for kind in ('axis','wall'):
            r=self.solve(kind,**self.uniform(kind))
            with self.subTest(kind=kind),self.assertRaises(DomainError):replace(r,numerical_policy=EqualText('forged'))

    def test_strict_common_gamma_selector_and_iteration_budget(self):
        for kind in ('axis','wall'):
            for name,value in (('gamma',True),('gamma',1.),('axisymmetric',1),('axisymmetric','true'),('max_iterations',1.),('max_iterations',True),('max_iterations',257)):
                with self.subTest(kind=kind,field=name,value=value),self.assertRaises(DomainError):
                    (solve_moc_axis_point if kind=='axis' else solve_moc_wall_point)(**{**self.uniform(kind),name:value})

    def test_axis_marching_history_remains_consistent_across_successive_reflections(self):
        errors=[]
        for step,count in ((.1,6),(.05,12),(.025,24)):
            history=[radial_point(1.8+i*step,0.,1.4,2) for i in (0,1)]
            for index in range(2,count+1):
                mach=1.8+index*step
                args,_target,_=fixture('axis',step=step,mach=mach)
                # Rescale the independent fixture from its local Rf=2 to the
                # common radial-flow normalization with R(M=2)=2.
                desired=radial_point(mach,0.,1.4,2);factor=desired.x_m/2
                p=args['incoming'];incoming=replace(p,x_m=p.x_m*factor,radius_m=p.radius_m*factor)
                r=self.solve('axis',incoming=incoming,axis_history=tuple(history[-2:]))
                history.append(r.point)
            expected=radial_point(2.4,0.,1.4,2)
            errors.append(max(abs(history[-1].x_m-expected.x_m),abs(history[-1].mach-2.4)))
        self.assertTrue(errors[2]<errors[1]<errors[0],errors)
        self.assertLess(errors[-1],1e-3)


if __name__=='__main__':unittest.main()
