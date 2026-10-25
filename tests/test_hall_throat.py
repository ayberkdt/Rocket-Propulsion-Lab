"""Primary-equation, independent derivative/root and hostile Hall-model checks."""
import itertools
import json
import unittest
from dataclasses import asdict, replace
from fractions import Fraction
from math import ldexp, nextafter, ulp
from pathlib import Path

import mpmath as mp

from rocket_propulsion.core.errors import ConvergenceError, DomainError
from rocket_propulsion.geometry import (
    HallSonicInputs,
    HallSonicResult,
    HallThroatInputs,
    HallThroatModel,
    HallThroatResult,
    evaluate_hall_throat,
    locate_hall_sonic_point,
)
from rocket_propulsion.geometry import hall_throat as kernel
from tests.test_moc_unit_process import exact_mp


def primary_velocity(w, radius, gamma, ratio, axisymmetric, order=2):
    """Independent printed equations, not the production generic coefficients.

    Planar: ARC R&M3440 A.5,A.6,A.8,A.9, retaining Cartesian terms to order 2.
    Axisymmetric: TRW02874 eqs4-1..4-4 in that report's different axial coordinate.
    All arguments must already be mpmath values under the caller's precision.
    """
    g, r, R = gamma, radius, ratio
    if not axisymmetric:
        q1 = r*r/2 - mp.mpf(1)/6 + w
        q2 = (g+6)*r**4/18 - (2*g+9)*r*r/18 + (g+30)/270 + (r*r-mp.mpf(1)/2)*w - (2*g-3)*w*w/6
        theta1 = r**3/6-r/6+r*w
        theta2 = (22*g+45)*r**5/360-(5*g+15)*r**3/54+(34*g+165)*r/1080+w*(2*g+6)*(r**3-r)/9
        u = 1+q1/R + (q2/R**2 if order==2 else 0)
        v = mp.sqrt(g+1)/R**mp.mpf('1.5')*(theta1 + ((theta2+q1*theta1)/R if order==2 else 0))
        return u,v
    A = mp.sqrt(2/(g+1))
    z = w/(A*R)
    u = 1+A*z+(-mp.mpf(1)/4+r*r/2)/R
    v = mp.sqrt(1/R)*(z*r+(-r+r**3)/(4*A*R))
    if order == 2:
        u += (3-2*g)*z*z/(3*(g+1))+A*z*(-mp.mpf(5)/8+r*r)/R
        u += ((10*g+57)/288-(4*g+15)*r*r/24+(2*g+9)*r**4/24)/R**2
        v += mp.sqrt(1/R)*(A*z*z*r + z*(-(4*g+15)*r/12+(2*g+9)*r**3/6)/R
             +((28*g+93)*r/288-(20*g+63)*r**3/96+(g+3)*r**5/9)/(A*R**2))
    return u,v


def reference_point(args, precision=110):
    with mp.workdps(precision):
        model=args['model'];g=exact_mp(model.gamma);rt=exact_mp(model.throat_radius_m)
        R=exact_mp(model.throat_curvature_radius_m)/rt;w=exact_mp(args['axial_parameter'])
        sigma=int(model.axisymmetric);k=(sigma+1)*R/(g+1)
        eta=exact_mp(args['radius_fraction']);wall2=1+w*w/(k*R)
        radius=eta*mp.sqrt(wall2) if args.get('radial_basis','throat')=='wall' else eta
        order=args.get('order',2)
        u,v=primary_velocity(w,radius,g,R,model.axisymmetric,order)
        sound2=(g+1-(g-1)*(u*u+v*v))/2
        return tuple(float(value) for value in (rt*w/mp.sqrt(k),rt*radius,u,v,mp.sqrt((u*u+v*v)/sound2)))


class HallThroatAttacks(unittest.TestCase):
    def model(self, **changes):
        return HallThroatModel(**{'throat_radius_m':.01,'throat_curvature_radius_m':.05,**changes})

    def evaluate(self, **args):
        try:
            return evaluate_hall_throat(**args)
        except (DomainError,ArithmeticError) as error:
            self.fail(f'Valid Hall witness refused: {error}')

    def sonic(self, **args):
        try:
            return locate_hall_sonic_point(**args)
        except (DomainError,ArithmeticError) as error:
            self.fail(f'Valid retained-model sonic witness refused: {error}')

    def test_both_geometries_match_independent_primary_cartesian_equations(self):
        for mode,gamma,order in itertools.product((False,True),(1.22,1.4,1.67),(1,2)):
            for w,r,basis in ((-.4,.2,'throat'),(.2,.7,'throat'),(.6,1.,'wall')):
                args={'model':self.model(axisymmetric=mode,gamma=gamma),'axial_parameter':w,'radius_fraction':r,'radial_basis':basis,'order':order}
                result=self.evaluate(**args);expected=reference_point(args)
                with self.subTest(mode=mode,gamma=gamma,w=w,r=r,basis=basis,order=order):
                    for actual,target in zip((result.x_m,result.radius_m,result.axial_velocity_over_sonic,result.radial_velocity_over_sonic,result.mach),expected):
                        self.assertLess(abs(actual/target-1),3e-14)
                    self.assertEqual(result.irrotational_residual,0.)

    def test_primary_reference_values_are_stable_at_110_and_180_digits(self):
        for mode in (False,True):
            args={'model':self.model(axisymmetric=mode),'axial_parameter':.6,'radius_fraction':.7,'radial_basis':'wall'}
            self.assertEqual(reference_point(args,110),reference_point(args,180))

    def test_nonlinear_potential_diagnostic_matches_independent_differentiation(self):
        for mode,order in itertools.product((False,True),(1,2)):
            model=self.model(axisymmetric=mode,gamma=1.22)
            result=self.evaluate(model=model,axial_parameter=.5,radius_fraction=.6,order=order)
            with mp.workdps(110):
                g=exact_mp(model.gamma);R=exact_mp(model.throat_curvature_radius_m)/exact_mp(model.throat_radius_m)
                k=(int(mode)+1)*R/(g+1);w=exact_mp(.5);r=exact_mp(.6);X=w/mp.sqrt(k)
                def velocity(x,y,*,k=k,g=g,R=R,mode=mode,order=order):return primary_velocity(mp.sqrt(k)*x,y,g,R,mode,order)
                u,v=velocity(X,r)
                ux=mp.diff(lambda x,velocity=velocity,r=r:velocity(x,r)[0],X)
                ur=mp.diff(lambda y,velocity=velocity,X=X:velocity(X,y)[0],r)
                vx=mp.diff(lambda x,velocity=velocity,r=r:velocity(x,r)[1],X)
                vr=mp.diff(lambda y,velocity=velocity,X=X:velocity(X,y)[1],r)
                a2=(g+1-(g-1)*(u*u+v*v))/2
                terms=((a2-u*u)*ux,-2*u*v*ur,(a2-v*v)*vr,int(mode)*a2*v/r)
                expected=abs(sum(terms))/(mp.sqrt(k)/R)
                relative=abs(sum(terms))/sum(map(abs,terms))
                self.assertLess(float(abs(ur-vx)),1e-100)
            self.assertLess(abs(result.potential_residual-float(expected)),3e-15)
            self.assertLess(abs(result.relative_term_balance-float(relative)),3e-15)
            self.assertGreater(result.potential_residual,0.)

    def test_second_order_defect_has_expected_asymptotic_improvement_in_both_modes(self):
        for mode in (False,True):
            defects=[]
            for R in (5.,10.,20.,40.):
                model=HallThroatModel(1.,R,axisymmetric=mode)
                pair=[self.evaluate(model=model,axial_parameter=.5,radius_fraction=.6,order=j).potential_residual for j in (1,2)]
                self.assertLess(pair[1],pair[0]/8)
                defects.append(pair)
            self.assertTrue(all(3<a[0]/b[0]<5 for a,b in itertools.pairwise(defects)),defects)
            self.assertTrue(all(7<a[1]/b[1]<10 for a,b in itertools.pairwise(defects)),defects)

    def test_throat_plane_is_subsonic_on_axis_and_supersonic_on_wall(self):
        for mode in (False,True):
            model=self.model(axisymmetric=mode)
            axis=self.evaluate(model=model,axial_parameter=0.,radius_fraction=0.)
            wall=self.evaluate(model=model,axial_parameter=0.,radius_fraction=1.,radial_basis='wall')
            self.assertLess(axis.mach,1.);self.assertEqual(axis.flow_regime,'subsonic')
            self.assertGreater(wall.mach,1.);self.assertEqual(wall.flow_regime,'supersonic')
            self.assertEqual(wall.wall_slip_velocity_over_sonic,0.)
            self.assertEqual(wall.wall_slip_status,'finite')

    def test_bracketed_sonic_axis_does_not_hide_its_nonzero_full_pde_defect(self):
        result=self.sonic(model=self.model(),wall_fraction=0.)
        self.assertLess(abs(result.point.mach-1.),2*ulp(1.))
        self.assertGreater(result.point.potential_residual,0.)
        self.assertGreater(result.point.relative_term_balance,.99)
        self.assertFalse(result.sonic_initialization_complete)
        self.assertFalse(result.point.sonic_initialization_complete)

    def test_wall_slip_away_from_throat_remains_visible(self):
        for mode in (False,True):
            result=self.sonic(model=self.model(axisymmetric=mode),wall_fraction=1.)
            self.assertGreater(result.point.wall_slip_velocity_over_sonic,1e-5)
            self.assertEqual(result.point.wall_slip_status,'finite')
            self.assertFalse(result.point.physical_model_qualified)

    def test_unrepresentable_wall_defect_sets_representation_flag_even_when_state_is_finite(self):
        model=HallThroatModel(1.,5.);w=ldexp(1.,-1070)
        result=self.evaluate(model=model,axial_parameter=w,radius_fraction=1.,radial_basis='wall')
        with mp.workdps(700):
            g=exact_mp(model.gamma);W=exact_mp(w);R=mp.mpf(5);k=2*R/(g+1)
            x=W/mp.sqrt(k);radius=mp.sqrt(1+x*x/R)
            u,v=primary_velocity(W,radius,g,R,True)
            slip=abs(v-x*u/(R*radius))
            self.assertGreater(slip,0.)
            self.assertLess(slip,exact_mp(ulp(0.))/2)
        self.assertEqual(result.radial_velocity_status,'finite')
        self.assertGreater(result.radial_velocity_over_sonic,0.)
        self.assertEqual(result.wall_slip_status,'unrepresentable')
        self.assertIsNone(result.wall_slip_velocity_over_sonic)
        self.assertTrue(result.representation_limited)

    def test_whole_sonic_curve_has_independently_enclosed_roots_and_real_wall_endpoint(self):
        for mode in (False,True):
            model=self.model(axisymmetric=mode)
            points=[]
            for eta in (0.,.25,.5,.75,1.):
                result=self.sonic(model=model,wall_fraction=eta);points.append(result.point)
                with mp.workdps(180):
                    g=exact_mp(model.gamma);rt=exact_mp(model.throat_radius_m);R=exact_mp(model.throat_curvature_radius_m)/rt
                    k=(int(mode)+1)*R/(g+1);fraction=exact_mp(eta)
                    def equation(w,*,fraction=fraction,k=k,R=R,g=g,mode=mode):
                        r=fraction*mp.sqrt(1+w*w/(R*k));u,v=primary_velocity(w,r,g,R,mode)
                        return u*u+v*v-1
                    root=mp.findroot(equation,(-mp.mpf('.5'),mp.mpf('.5')))
                    x=rt*root/mp.sqrt(k)
                    radius=rt*fraction*mp.sqrt(1+root*root/(R*k))
                    self.assertLessEqual(exact_mp(result.parameter_lower),root)
                    self.assertGreaterEqual(exact_mp(result.parameter_upper),root)
                    self.assertLessEqual(exact_mp(result.x_lower_m),x)
                    self.assertGreaterEqual(exact_mp(result.x_upper_m),x)
                    self.assertLessEqual(exact_mp(result.radius_lower_m),radius)
                    self.assertGreaterEqual(exact_mp(result.radius_upper_m),radius)
                    self.assertLessEqual(exact_mp(result.sonic_function_lower),equation(exact_mp(result.parameter_lower)))
                    self.assertGreaterEqual(exact_mp(result.sonic_function_upper),equation(exact_mp(result.parameter_upper)))
                self.assertLessEqual(result.sonic_function_lower,0.)
                self.assertGreaterEqual(result.sonic_function_upper,0.)
                self.assertLessEqual(result.x_lower_m,result.point.x_m)
                self.assertLessEqual(result.point.x_m,result.x_upper_m)
                self.assertLessEqual(result.radius_lower_m,result.point.radius_m)
                self.assertLessEqual(result.point.radius_m,result.radius_upper_m)
            self.assertGreater(points[0].x_m,0.);self.assertLess(points[-1].x_m,0.)
            self.assertTrue(all(b.radius_m>a.radius_m for a,b in itertools.pairwise(points)))
            self.assertEqual(points[-1].radius_m,points[-1].wall_radius_m)
            self.assertGreater(points[-1].radius_m,model.throat_radius_m)

    def test_published_theory_cells_match_original_printed_rounding_width(self):
        data=json.loads((Path(__file__).resolve().parents[1]/'docs/assurance/hall_throat_theory_reference.json').read_text(encoding='utf-8'))
        compared=0
        for row in data['rows']:
            if row['status']!='equation_comparison':continue
            result=self.sonic(model=HallThroatModel(1.,1/row['inverse_curvature_ratio'],gamma=row['gamma']),wall_fraction=0.,order=row['order'])
            self.assertLessEqual(abs(result.point.x_m-row['axis_sonic_x_over_rt']),data['rounding_half_width'])
            compared+=1
        self.assertEqual(compared,11)

    def test_inconsistent_source_cell_is_flagged_instead_of_fitting_the_coefficients(self):
        data=json.loads((Path(__file__).resolve().parents[1]/'docs/assurance/hall_throat_theory_reference.json').read_text(encoding='utf-8'))
        flagged=[r for r in data['rows'] if r['status']!='equation_comparison']
        self.assertEqual(len(flagged),1)
        row=flagged[0];self.assertEqual(row['status'],'inconsistent_with_equation_4_1')
        self.assertEqual(row['axis_sonic_x_over_rt'],.0856)
        result=self.sonic(model=HallThroatModel(1.,10.),wall_fraction=0.,order=1)
        with mp.workdps(110):expected=float(mp.sqrt(exact_mp(1.4)+1)/(4*mp.sqrt(20)))
        self.assertLess(abs(result.point.x_m-expected),2e-16)
        self.assertGreater(abs(result.point.x_m-row['axis_sonic_x_over_rt']),20*data['rounding_half_width'])

    def test_coordinate_scaling_changes_geometry_but_not_state_or_root_certificates(self):
        model=self.model();base=self.sonic(model=model,wall_fraction=.7)
        for exponent in (-900,900):
            scaled=replace(model,throat_radius_m=ldexp(model.throat_radius_m,exponent),throat_curvature_radius_m=ldexp(model.throat_curvature_radius_m,exponent))
            result=self.sonic(model=scaled,wall_fraction=.7)
            for name in ('parameter_lower','parameter_upper','sonic_function_lower','sonic_function_upper','iterations'):
                self.assertEqual(getattr(result,name),getattr(base,name))
            for name in ('x_m','radius_m','mach','axial_velocity_over_sonic','potential_residual'):
                expected=ldexp(getattr(base.point,name),exponent) if name in ('x_m','radius_m') else getattr(base.point,name)
                self.assertEqual(getattr(result.point,name),expected)
            for name in ('x_lower_m','x_upper_m','radius_lower_m','radius_upper_m'):
                self.assertEqual(getattr(result,name),ldexp(getattr(base,name),exponent))

    def test_radial_enclosure_accounts_for_the_interior_minimum_across_zero(self):
        model=HallThroatModel(1.,5.)
        lo,hi=kernel._radial_bounds(model,Fraction(1,2),Fraction(-1,2),Fraction(1,2))
        self.assertEqual(lo,.5)
        with mp.workdps(110):
            expected=mp.sqrt(1+mp.mpf('.2')*mp.mpf('.25')/(10/exact_mp(1.4+1)))/2
        self.assertGreaterEqual(exact_mp(hi),expected)
        self.assertLessEqual(exact_mp(nextafter(hi,-float('inf'))),expected)
        for left,right in ((Fraction(0),Fraction(1,2)),(Fraction(1,2),Fraction(1)),(Fraction(-1),Fraction(-1,2))):
            lo,hi=kernel._radial_bounds(model,Fraction(1,2),left,right)
            with mp.workdps(110):
                k=10/(exact_mp(1.4)+1)
                for w in (left,(left+right)/2,right):
                    radius=mp.sqrt(1+exact_mp(w)**2/(5*k))/2
                    self.assertLessEqual(exact_mp(lo),radius)
                    self.assertGreaterEqual(exact_mp(hi),radius)

    def test_radical_coordinate_enclosures_are_outward_for_both_signs_and_extreme_scales(self):
        for exponent in (-1000,0,1000):
            square=Fraction(2)*Fraction(2)**(2*exponent)
            for sign in (-1,1):
                lo,hi=kernel._root_bounds(square,sign)
                self.assertEqual(nextafter(lo,float('inf')),hi)
                if sign>0:
                    self.assertLessEqual(Fraction(lo)**2,square);self.assertGreaterEqual(Fraction(hi)**2,square)
                else:
                    self.assertGreaterEqual(Fraction(lo)**2,square);self.assertLessEqual(Fraction(hi)**2,square)

    def test_representable_sonic_position_retains_unrepresentable_field_perturbation(self):
        model=HallThroatModel(1.,1e308)
        result=self.sonic(model=model,wall_fraction=0.)
        self.assertGreater(result.point.x_m,0.)
        self.assertEqual(result.point.mach,1.)
        self.assertEqual(result.point.flow_regime,'supersonic')
        self.assertIsNone(result.point.axial_perturbation_over_sonic)
        self.assertTrue(result.point.representation_limited)
        with mp.workdps(700):
            g=exact_mp(model.gamma);eps=1/exact_mp(model.throat_curvature_radius_m)
            root=mp.findroot(lambda w:w-mp.mpf('.25')+eps*((3-2*g)*w*w/6-mp.mpf('.625')*w+(10*g+57)/288),mp.mpf('.25'))
            x=root*mp.sqrt((g+1)*eps/2)
            self.assertLessEqual(exact_mp(result.parameter_lower),root)
            self.assertGreaterEqual(exact_mp(result.parameter_upper),root)
            self.assertLessEqual(exact_mp(result.x_lower_m),x)
            self.assertGreaterEqual(exact_mp(result.x_upper_m),x)

    def test_exact_retained_model_sonic_root_at_zero_needs_no_mach_offset(self):
        result=self.sonic(model=HallThroatModel(1.,1.,gamma=1.5),wall_fraction=0.,max_iterations=1)
        self.assertEqual(result.point.inputs.axial_parameter,0.)
        self.assertEqual((result.point.x_m,result.x_lower_m,result.x_upper_m),(0.,0.,0.))
        self.assertEqual(result.point.mach,1.);self.assertEqual(result.point.flow_regime,'sonic')
        self.assertEqual(result.iterations,1)
        self.assertGreater(result.point.potential_residual,.1)
        self.assertFalse(result.sonic_initialization_complete)

    def test_private_perturbations_survive_extreme_curvature_ratio_and_erased_displays(self):
        model=HallThroatModel(ldexp(1.,-1000),ldexp(1.,1000))
        result=self.evaluate(model=model,axial_parameter=0.,radius_fraction=.5)
        self.assertEqual(result.mach,1.);self.assertEqual(result.flow_regime,'subsonic')
        self.assertTrue(result.representation_limited)
        self.assertIsNone(result.radial_velocity_over_sonic);self.assertIsNone(result.axial_perturbation_over_sonic)
        self.assertEqual(result.radial_velocity_status,'unrepresentable');self.assertEqual(result.axial_perturbation_status,'unrepresentable')
        self.assertGreater(result.potential_residual,0.)

    def test_unrepresentable_physical_sonic_position_is_not_exported_as_zero(self):
        with self.assertRaises(DomainError):
            locate_hall_sonic_point(model=HallThroatModel(ldexp(1.,-1000),ldexp(1.,1000)),wall_fraction=0.)

    def test_small_curvature_unbracketed_or_uncertified_root_is_not_blessed(self):
        for model,eta,code in ((HallThroatModel(1.,.1),0.,'hall_sonic_not_bracketed'),
                               (HallThroatModel(1.,2.,axisymmetric=False),1.,'hall_sonic_branch_unresolved')):
            with self.subTest(code=code):
                with self.assertRaises(ConvergenceError) as context:locate_hall_sonic_point(model=model,wall_fraction=eta)
                self.assertEqual(context.exception.code,code)

    def test_exact_bernstein_certificate_cannot_accept_hidden_negative_regions(self):
        self.assertTrue(kernel._positive_on_interval((Fraction(1),),Fraction(-1),Fraction(1)))
        self.assertFalse(kernel._positive_on_interval((),Fraction(-1),Fraction(1)))
        self.assertFalse(kernel._positive_on_interval((Fraction(0),),Fraction(-1),Fraction(1)))
        self.assertFalse(kernel._positive_on_interval((Fraction(-1,1000),Fraction(0),Fraction(1)),Fraction(-1),Fraction(1)))

    def test_root_iteration_budget_must_fail_instead_of_accepting_an_unresolved_interval(self):
        with self.assertRaises(ConvergenceError) as context:
            locate_hall_sonic_point(model=self.model(),wall_fraction=.7,max_iterations=1)
        self.assertEqual(context.exception.code,'hall_sonic_not_converged')

    def test_model_inputs_refuse_nonphysical_and_malformed_scalars(self):
        for field in ('throat_radius_m','throat_curvature_radius_m','gamma'):
            for value in (True,None,'1',float('nan'),float('inf'),0.,-1.):
                with self.subTest(field=field,value=value),self.assertRaises(DomainError):self.model(**{field:value})
        with self.assertRaises(DomainError):self.model(gamma=1.)
        with self.assertRaises(DomainError):self.model(axisymmetric=1)

    def test_coordinate_order_and_root_resources_have_strict_types(self):
        model=self.model()
        for changes in ({'axial_parameter':True},{'radius_fraction':-.1},{'radius_fraction':2.},{'radial_basis':'flat'},
                        {'order':True},{'order':1.},{'order':0},{'order':3}):
            with self.subTest(changes=changes),self.assertRaises(DomainError):
                HallThroatInputs(**{'model':model,'axial_parameter':0.,'radius_fraction':.5,**changes})
        for changes in ({'wall_fraction':-.1},{'wall_fraction':1.1},{'wall_fraction':True},{'order':3},
                        {'max_iterations':True},{'max_iterations':1.},{'max_iterations':0},{'max_iterations':257}):
            with self.subTest(changes=changes),self.assertRaises(DomainError):HallSonicInputs(**{'model':model,'wall_fraction':.5,**changes})

    def test_polynomial_branch_leaving_positive_energy_is_refused(self):
        with self.assertRaises(DomainError):
            evaluate_hall_throat(model=self.model(),axial_parameter=100.,radius_fraction=0.)

    def test_state_and_sonic_mapping_round_trips_own_all_inputs(self):
        result=self.sonic(model=self.model(),wall_fraction=.7);raw=asdict(result)
        copied=HallSonicResult(**raw);self.assertEqual(copied,result)
        raw['inputs']['model']['gamma']=2.;raw['point']['inputs']['model']['gamma']=3.
        self.assertNotEqual(copied.inputs.model.gamma,2.);self.assertNotEqual(copied.point.inputs.model.gamma,3.)
        self.assertEqual(HallThroatResult(**asdict(result.point)),result.point)

    def test_every_state_diagnostic_and_optional_display_must_replay(self):
        result=self.evaluate(model=self.model(),axial_parameter=.5,radius_fraction=1.,radial_basis='wall')
        changes=({'mach':2.},{'x_m':0.},{'radius_m':0.},{'wall_radius_m':1.},{'axial_velocity_over_sonic':1.},
            {'radial_velocity_over_sonic':0.},{'radial_velocity_status':'unrepresentable'},
            {'axial_perturbation_over_sonic':None},{'axial_perturbation_status':'unrepresentable'},
            {'potential_residual':0.},{'relative_term_balance':0.},{'irrotational_residual':1.},
            {'wall_slip_velocity_over_sonic':0.},{'wall_slip_status':'not_wall'},
            {'flow_regime':'sonic'},{'representation_limited':True},{'flow_angle_deg':0.})
        for change in changes:
            with self.subTest(change=change),self.assertRaises(DomainError):replace(result,**change)

    def test_root_state_brackets_iterations_and_signs_must_replay(self):
        result=self.sonic(model=self.model(),wall_fraction=.7)
        point=self.evaluate(model=self.model(),axial_parameter=result.point.inputs.axial_parameter+.05,
                            radius_fraction=.7,radial_basis='wall')
        with self.assertRaises(DomainError):replace(result,point=point)
        other=self.evaluate(model=self.model(),axial_parameter=result.point.inputs.axial_parameter,
                            radius_fraction=.7,radial_basis='wall',order=1)
        with self.assertRaises(DomainError):replace(result,point=other)

    def test_forged_root_brackets_and_evidence_flags_are_refused(self):
        result=self.sonic(model=self.model(),wall_fraction=.7)
        for change in ({'iterations':result.iterations+1},{'parameter_lower':result.parameter_lower-.01},
                       {'parameter_lower':nextafter(result.parameter_lower,-float('inf'))},
                       {'parameter_upper':result.parameter_upper+.01},{'x_lower_m':0.},{'x_upper_m':0.},
                       {'radius_lower_m':0.},{'radius_upper_m':0.},
                       {'sonic_function_lower':-1.},{'sonic_function_upper':1.}):
            with self.subTest(change=change),self.assertRaises(DomainError):replace(result,**change)
        for record in (result,result.point):
            for field in ('physical_model_qualified','sonic_initialization_complete'):
                for value in (True,0,1):
                    with self.subTest(field=field,value=value),self.assertRaises(DomainError):replace(record,**{field:value})
            for change in ({'approximation_error_bound':0.},{'numerical_policy':'certified'}):
                with self.subTest(change=change),self.assertRaises(DomainError):replace(record,**change)

    def test_hostile_policy_equality_cannot_authenticate_a_forged_record(self):
        class EqualText(str):
            def __eq__(self,other):return True
        result=self.sonic(model=self.model(),wall_fraction=.7)
        for record in (result,result.point):
            with self.assertRaises(DomainError):replace(record,numerical_policy=EqualText('forged'))

    def test_json_input_ownership_does_not_drop_radial_basis_or_order(self):
        inputs=HallThroatInputs(self.model(),.5,.7,'wall',1)
        restored=HallThroatInputs(**json.loads(json.dumps(asdict(inputs))))
        self.assertEqual(restored,inputs)
        self.assertEqual(evaluate_hall_throat(**asdict(restored)).inputs,inputs)
