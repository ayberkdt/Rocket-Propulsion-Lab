"""Manufactured radial flow, independent Cartesian balances and hostile boundaries."""
from collections.abc import Mapping
from dataclasses import asdict, replace
from fractions import Fraction
from math import ldexp, nextafter, sqrt, ulp
import unittest
import mpmath as mp

from rocket_propulsion.core.errors import ConvergenceError, DomainError
from rocket_propulsion.geometry import MocFlowPoint, solve_moc_interior_point


def exact_mp(value):
    q=Fraction(value);return mp.mpf(q.numerator)/q.denominator


def manufactured(*,axisymmetric=True,step=.05,mach=2.,angle=20.,gamma=1.4,precision=110):
    """Exact cylindrical/spherical radial flow and its analytic characteristics.

    rho*V*R**d=constant, d=1 (planar) or 2 (axisymmetric), gives
    R/Rf=(A(M)/A(Mf))**(1/d). Along each Mach line the polar angle
    change is +/- (nu_f-nu)/d. No production MOC equations/roots are called.
    """
    with mp.workdps(precision):
        g=exact_mp(gamma);m=exact_mp(mach);p=m-exact_mp(step);theta=exact_mp(angle)*mp.pi/180
        d=2 if axisymmetric else 1
        def area(z):return (2/(g+1)*(1+(g-1)*z*z/2))**((g+1)/(2*(g-1)))/z
        def nu(z):return mp.sqrt((g+1)/(g-1))*mp.atan(mp.sqrt((g-1)/(g+1)*(z*z-1)))-mp.atan(mp.sqrt(z*z-1))
        radius=2*(area(p)/area(m))**(mp.mpf(1)/d);delta=(nu(m)-nu(p))/d
        parents=[]
        for sign in (-1,1):
            phi=theta+sign*delta
            parents.append(MocFlowPoint(float(radius*mp.cos(phi)),float(radius*mp.sin(phi)),float(p),float(phi*180/mp.pi)))
        expected=(float(2*mp.cos(theta)),float(2*mp.sin(theta)),mach,angle)
    return dict(plus_parent=parents[0],minus_parent=parents[1],gamma=gamma,axisymmetric=axisymmetric),expected


def error(point,expected):
    return max(abs(point.x_m-expected[0])/2,abs(point.radius_m-expected[1])/2,
               abs(point.mach-expected[2])/expected[2],abs(point.flow_angle_deg-expected[3])/expected[3])


class MocUnitProcessAttacks(unittest.TestCase):
    def valid(self,**inputs):
        try:return solve_moc_interior_point(**inputs)
        except (DomainError,ArithmeticError) as exc:self.fail(f'Representable characteristic witness rejected: {exc}')

    def uniform(self,mach=2.,axisymmetric=True):
        return dict(plus_parent=MocFlowPoint(0,1,mach,0),minus_parent=MocFlowPoint(0,3,mach,0),axisymmetric=axisymmetric)

    def test_parallel_uniform_flow_matches_closed_form_in_both_geometries(self):
        for axis in (False,True):
            for mach in (1.01,2.,8.):
                with self.subTest(axis=axis,mach=mach):
                    r=self.valid(**self.uniform(mach,axis))
                    self.assertLess(abs(r.point.x_m/sqrt(mach*mach-1)-1),8*ulp(1.))
                    self.assertEqual((r.point.radius_m,r.point.mach,r.point.flow_angle_deg),(2.,mach,0.))
                    self.assertEqual(r.solution_kind,'constant_state_intersection')

    def test_inclined_uniform_planar_flow_matches_independent_intersection(self):
        a=MocFlowPoint(0,1,2.3,7.);b=replace(a,radius_m=3.)
        r=self.valid(plus_parent=a,minus_parent=b,axisymmetric=False)
        with mp.workdps(110):
            theta=exact_mp(7.)*mp.pi/180;mu=mp.asin(1/exact_mp(2.3))
            plus=mp.tan(theta+mu);minus=mp.tan(theta-mu);x=2/(plus-minus);y=1+plus*x
        self.assertLess(abs(r.point.x_m/float(x)-1),1e-14)
        self.assertLess(abs(r.point.radius_m/float(y)-1),1e-14)
        self.assertEqual((r.point.mach,r.point.flow_angle_deg),(2.3,7.))

    def convergence(self,axis):
        for gamma,mach,angle in ((1.4,2.,20.),(1.22,3.,12.)):
            errors=[]
            for step in (.1,.05,.025,.0125,.00625):
                x,wanted=manufactured(axisymmetric=axis,step=step,gamma=gamma,mach=mach,angle=angle)
                r=self.valid(**x);errors.append(error(r.point,wanted))
                self.assertEqual(r.solution_kind,'iterated_midpoint_intersection')
            self.assertTrue(all(6<left/right<12 for left,right in zip(errors,errors[1:])),errors)
            self.assertLess(errors[-1],2e-7)

    def test_axisymmetric_spherical_flow_has_third_order_local_refinement(self):self.convergence(True)

    def test_planar_radial_flow_has_third_order_local_refinement(self):self.convergence(False)

    def test_independent_oracle_is_stable_at_110_and_180_digits(self):
        for axis in (False,True):
            self.assertEqual(manufactured(axisymmetric=axis,precision=110),manufactured(axisymmetric=axis,precision=180))

    def test_cylindrical_source_term_cannot_be_replaced_by_planar_compatibility(self):
        inputs,wanted=manufactured(step=.05)
        correct=self.valid(**inputs);planar=self.valid(**{**inputs,'axisymmetric':False})
        self.assertLess(error(correct.point,wanted),1e-5)
        self.assertGreater(error(planar.point,wanted),1e-3)

    def test_actual_unit_solution_satisfies_independent_cartesian_equations(self):
        for axis in (False,True):
            inputs,_=manufactured(axisymmetric=axis);r=self.valid(**inputs)
            with mp.workdps(110):
                g=exact_mp(inputs['gamma'])
                def state(p,sign):
                    m=exact_mp(p.mach);theta=exact_mp(p.flow_angle_deg)*mp.pi/180
                    a2=(g+1)/(2+(g-1)*m*m);v=m*mp.sqrt(a2)
                    u,w=v*mp.cos(theta),v*mp.sin(theta);lam=mp.tan(theta+sign*mp.asin(1/m))
                    aa=u*u-a2;bb=2*u*w-aa*lam
                    source=a2*w/exact_mp(p.radius_m) if axis else 0
                    return lam,aa,bb,source,u,w
                for parent,sign in ((inputs['plus_parent'],1),(inputs['minus_parent'],-1)):
                    left,right=state(parent,sign),state(r.point,sign)
                    lam,aa,bb,source=[(a+b)/2 for a,b in zip(left[:4],right[:4])]
                    dx=exact_mp(r.point.x_m)-exact_mp(parent.x_m);dy=exact_mp(r.point.radius_m)-exact_mp(parent.radius_m)
                    self.assertLess(float(abs(dy-lam*dx)),1e-12)
                    balance=aa*(right[4]-left[4])+bb*(right[5]-left[5])-source*dx
                    self.assertLess(float(abs(balance)),1e-12)

    def test_common_power_of_two_scaling_preserves_flow_and_iterations(self):
        inputs,_=manufactured();base=self.valid(**inputs)
        for exponent in (-1000,-450,450,1000):
            scaled=dict(inputs)
            for name in ('plus_parent','minus_parent'):
                p=inputs[name];scaled[name]=replace(p,x_m=ldexp(p.x_m,exponent),radius_m=ldexp(p.radius_m,exponent))
            r=self.valid(**scaled)
            self.assertEqual((r.point.mach,r.point.flow_angle_deg,r.iterations),(base.point.mach,base.point.flow_angle_deg,base.iterations))
            self.assertEqual(ldexp(r.point.x_m,-exponent),base.point.x_m)
            self.assertEqual(ldexp(r.point.radius_m,-exponent),base.point.radius_m)

    def test_uniform_high_mach_does_not_square_mach_in_binary64(self):
        r=self.valid(**self.uniform(1e200))
        self.assertLess(abs(r.point.x_m/1e200-1),1e-14)
        self.assertEqual(r.point.mach,1e200)

    def test_unsupported_axis_boundary_is_refused_but_planar_zero_y_is_valid(self):
        x=self.uniform();x['plus_parent']=replace(x['plus_parent'],radius_m=0.)
        with self.assertRaises(DomainError):solve_moc_interior_point(**x)
        r=self.valid(**{**x,'axisymmetric':False});self.assertGreater(r.point.x_m,0)

    def test_parent_family_order_cannot_return_an_upstream_branch(self):
        x=self.uniform()
        with self.assertRaises(ConvergenceError) as context:
            solve_moc_interior_point(plus_parent=x['minus_parent'],minus_parent=x['plus_parent'])
        self.assertEqual(context.exception.code,'moc_unit_noncausal')

    def test_distinct_parent_geometry_is_required(self):
        x=self.uniform();x['minus_parent']=replace(x['plus_parent'],mach=2.5)
        with self.assertRaises(DomainError):solve_moc_interior_point(**x)

    def test_iteration_budget_is_enforced_without_returning_unconverged_state(self):
        inputs,_=manufactured();self.valid(**inputs)
        with self.assertRaises(ConvergenceError):solve_moc_interior_point(**inputs,max_iterations=1)

    def test_tiny_nonconstant_parent_difference_is_not_labelled_exact(self):
        x=self.uniform();x['minus_parent']=replace(x['minus_parent'],mach=nextafter(2.,3.))
        self.assertEqual(self.valid(**x).solution_kind,'iterated_midpoint_intersection')

    def test_point_scalars_are_strict_and_sonic_or_transverse_states_are_refused(self):
        for name in ('x_m','radius_m','mach','flow_angle_deg'):
            for value in (True,None,'2',float('nan'),float('inf')):
                with self.subTest(field=name,value=value),self.assertRaises(DomainError):replace(self.uniform()['plus_parent'],**{name:value})
        for change in (dict(mach=1.),dict(mach=.8),dict(flow_angle_deg=90.),dict(flow_angle_deg=70.)):
            with self.subTest(change=change),self.assertRaises(DomainError):replace(self.uniform()['plus_parent'],**change)

    def test_selector_gamma_and_iteration_count_do_not_accept_python_truthiness(self):
        for name,values in (('axisymmetric',(1,0,'true',None)),('gamma',(True,1.,None,'1.4')),
                            ('max_iterations',(True,1.,'128',None,0,257))):
            for value in values:
                with self.subTest(field=name,value=value),self.assertRaises(DomainError):solve_moc_interior_point(**{**self.uniform(),name:value})

    def test_mapping_reopening_owns_nested_data(self):
        r=self.valid(**manufactured()[0]);raw=asdict(r);restored=type(r)(**raw)
        self.assertEqual(restored,r);raw['inputs']['plus_parent']['mach']=99.
        self.assertNotEqual(restored.inputs.plus_parent.mach,99.)

    def test_retained_state_residual_iteration_and_policy_cannot_be_forged(self):
        r=self.valid(**manufactured()[0])
        for change in (dict(point=replace(r.point,mach=r.point.mach+.01)),dict(iterations=r.iterations+1),
                       dict(equation_residual_max=0.),dict(solution_kind='constant_state_intersection'),dict(numerical_policy='forged')):
            with self.subTest(change=change),self.assertRaises(DomainError):replace(r,**change)

    def test_partial_unit_process_cannot_claim_error_bound_or_nozzle_qualification(self):
        r=self.valid(**manufactured()[0])
        self.assertIsNone(r.discretization_error_bound)
        for change in (dict(physical_model_qualified=True),dict(physical_model_qualified=0),
                       dict(full_nozzle_solution=True),dict(discretization_error_bound=0.)):
            with self.subTest(change=change),self.assertRaises(DomainError):replace(r,**change)

    def test_malformed_mappings_are_typed_and_unrelated_callbacks_propagate(self):
        class Broken(Mapping):
            def __init__(self,error):self.error=error
            def __iter__(self):raise self.error('callback')
            def __getitem__(self,key):raise self.error('callback')
            def __len__(self):return 1
        for kind in (KeyError,TypeError):
            with self.subTest(error=kind),self.assertRaises(DomainError):solve_moc_interior_point(**{**self.uniform(),'plus_parent':Broken(kind)})
        with self.assertRaises(RuntimeError):solve_moc_interior_point(**{**self.uniform(),'plus_parent':Broken(RuntimeError)})

    def test_hostile_string_equality_cannot_authenticate_a_forged_policy(self):
        class EqualText(str):
            def __eq__(self,other):return True
            def __ne__(self,other):return False
        r=self.valid(**self.uniform())
        with self.assertRaises(DomainError):replace(r,numerical_policy=EqualText('forged'))

    def test_coordinate_rounding_is_visible_without_relaxing_compatibility(self):
        r=self.valid(**manufactured(step=.00625)[0])
        self.assertIs(r.representation_limited,True)
        self.assertGreater(r.equation_residual_max,1.5e-14)
        self.assertLess(r.rounding_adjusted_residual_max,2e-15)
        self.assertGreater(r.coordinate_rounding_allowance_max,r.equation_residual_max-1.4e-14)
        for change in (dict(coordinate_rounding_allowance_max=0.),dict(representation_limited=False),dict(representation_limited=1)):
            with self.subTest(change=change),self.assertRaises(DomainError):replace(r,**change)

    def test_large_translated_planar_intersection_keeps_private_products_finite(self):
        a=MocFlowPoint(-1e308,-1e308,1.1,0.);b=replace(a,radius_m=1e308)
        r=self.valid(plus_parent=a,minus_parent=b,axisymmetric=False)
        with mp.workdps(110):
            expected=exact_mp(-1e308)+exact_mp(1e308)*mp.sqrt(exact_mp(1.1)**2-1)
        self.assertLess(abs(r.point.x_m/float(expected)-1),1e-14)
        self.assertEqual(r.point.radius_m,0.)
        self.assertEqual(r.point.mach,1.1)

    def test_coordinate_allowance_matches_independent_high_precision_cell_geometry(self):
        inputs,_=manufactured(step=.00625);r=self.valid(**inputs);p=r.point
        with mp.workdps(110):
            def half_cell(v):return max(abs(exact_mp(nextafter(v,float('inf')))-exact_mp(v)),
                                        abs(exact_mp(v)-exact_mp(nextafter(v,float('-inf')))))/2
            a,b=inputs['plus_parent'],inputs['minus_parent']
            scale=max(abs(exact_mp(a.x_m)-exact_mp(b.x_m)),abs(exact_mp(a.radius_m)-exact_mp(b.radius_m)))
            bounds=[]
            for parent,sign in ((a,1),(b,-1)):
                def slope(point):return mp.tan(exact_mp(point.flow_angle_deg)*mp.pi/180+sign*mp.asin(1/exact_mp(point.mach)))
                average=(slope(parent)+slope(p))/2
                bounds.append((half_cell(p.radius_m)+abs(average)*half_cell(p.x_m))/scale)
            expected=float(max(bounds))
        self.assertLess(abs(r.coordinate_rounding_allowance_max/expected-1),5e-14)

    def test_coordinate_cells_cannot_legitimize_collapsed_downstream_topology(self):
        self.valid(**self.uniform(1.1))
        a=MocFlowPoint(1e308,1.,1.1,0.);b=replace(a,radius_m=3.)
        with self.assertRaises(ConvergenceError) as context:
            solve_moc_interior_point(plus_parent=a,minus_parent=b)
        self.assertEqual(context.exception.code,'moc_unit_unresolved_topology')


if __name__=='__main__':unittest.main()
