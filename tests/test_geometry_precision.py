"""Independent high-precision witnesses for near-sonic and tiny geometry.

Arithmetic range is tested separately from empirical flow/hardware validity.
No test encodes the old refusal of a balanced finite input as expected behavior.
"""
import json
import math
import random
import unittest
from dataclasses import asdict
from decimal import ROUND_DOWN, DefaultContext, FloatOperation, Inexact, localcontext
from fractions import Fraction
from functools import cache

import mpmath as mp

from rocket_propulsion.compressible.isentropic import calculate_isentropic, mach_from_isentropic
from rocket_propulsion.core.errors import ConvergenceError, DomainError
from rocket_propulsion.geometry import (
    export_nozzle_contour,
    generate_moc_nozzle,
    generate_nozzle_contour,
    generate_rao_contour,
)
from rocket_propulsion.thermodynamics import snapshot_thermal_contour

GENERATORS = (generate_nozzle_contour,generate_rao_contour,generate_moc_nozzle)
NEXT_AREA = math.nextafter(1.,math.inf)
TINY_AREA = math.nextafter(0.,1.)


@cache
def geometry(generator, area=.01, ratio=10., count=5):
    args={'throat_area_m2':area,'area_ratio':ratio}
    args.update({'characteristic_count':count} if generator is generate_moc_nozzle else {'station_count':41})
    return generator(**args)


def contour(value):
    return value.contour if hasattr(value,'characteristic_lines') else value


def pm_integral(mach,gamma):
    """Quadrature of d(nu)/dM after t=sqrt(M-1); scaled before integration."""
    g,m=mp.mpf(gamma),mp.mpf(mach)
    return mp.quad(lambda t: (g+1)*2*t*t*mp.sqrt(2+t*t)
        /((1+t*t)*(1+(g-1)*(1+t*t)**2/2)),[0,mp.sqrt(m-1)])/(g+1)*180/mp.pi


def inverse_area(ratio,gamma,subsonic):
    g,a=mp.mpf(gamma),mp.mpf(ratio)
    first=mp.sqrt((g+1)*mp.log(a)/2)*(-1 if subsonic else 1)
    def residual(delta):
        m=1+delta
        return ((g+1)/(2*(g-1)))*mp.log((2+(g-1)*m*m)/(g+1))-mp.log(m)-mp.log(a)
    return 1+mp.findroot(residual,(first*mp.mpf('.8'),first*mp.mpf('1.2')))


def moc_reference_length(area,ratio,gamma,count):
    """Independent 100-digit compatibility net with normal-vector matrix intersections.

    Same planar fan discretization, independent roots/coordinates and dimensional
    mapping. This checks numerical realization, not the physical MOC model.
    """
    g=mp.mpf(gamma)
    exit_mach=inverse_area(ratio,gamma,False)
    def pm(m):
        root=mp.sqrt(m*m-1)
        factor=mp.sqrt((g+1)/(g-1))
        return factor*mp.atan(root/factor)-mp.atan(root)
    turn=pm(exit_mach)/2

    @cache
    def state(theta,nu):
        if nu==0:
            mach=mp.mpf(1)
        else:
            first=(3*(g+1)*nu/(4*mp.sqrt(2)))**(mp.mpf(2)/3)
            mach=1+mp.findroot(lambda d:pm(1+d)-nu,(first*mp.mpf('.5'),first*2))
        return theta,nu,mp.asin(1/mach)

    def intersection(first,angle1,second,angle2):
        normals=mp.matrix([[-mp.sin(angle1),mp.cos(angle1)],
                           [-mp.sin(angle2),mp.cos(angle2)]])
        rhs=mp.matrix([normals[0,0]*first[0]+normals[0,1]*first[1],
                       normals[1,0]*second[0]+normals[1,1]*second[1]])
        result=mp.lu_solve(normals,rhs)
        return result[0],result[1]

    rows=[];walls=[((mp.mpf(0),mp.mpf(1)),turn)]
    for row_index in range(count):
        fan=(row_index+1)*turn/count
        row=[]
        for column in range(count+1-row_index):
            if column==0:
                theta,nu=mp.mpf(0),2*fan
            elif column==count-row_index:
                theta=turn-fan
                nu=theta-(row[-1][1]-row[-1][2])
            else:
                plus=row[-1][1]-row[-1][2]
                minus=2*(row_index+1+column)*turn/count
                theta,nu=(plus+minus)/2,(minus-plus)/2
            theta,nu,mu=state(theta,nu)
            if column==0:
                origin,stheta,_,smu=(((mp.mpf(0),mp.mpf(1)),*state(fan,fan))
                                    if row_index==0 else rows[-1][1])
                xy=intersection(origin,(stheta-smu+theta-mu)/2,(0,0),mp.mpf(0))
            elif column==count-row_index:
                lower,ltheta,_,lmu=row[-1]
                wall,wtheta=walls[-1]
                xy=intersection(wall,(wtheta+theta)/2,lower,(ltheta+lmu+theta+mu)/2)
                walls.append((xy,theta))
            else:
                lower,ltheta,_,lmu=row[-1]
                upper,utheta,_,umu=(((mp.mpf(0),mp.mpf(1)),
                    *state((row_index+1+column)*turn/count,(row_index+1+column)*turn/count))
                    if row_index==0 else rows[-1][column+1])
                xy=intersection(lower,(ltheta+lmu+theta+mu)/2,upper,(utheta-umu+theta-mu)/2)
            row.append((xy,theta,nu,mu))
        rows.append(row)
    xy=walls[-1][0]
    radius=mp.sqrt(mp.mpf(area)/mp.mpf(math.pi))
    return xy[0]*radius*(mp.sqrt(mp.mpf(ratio))-1)/(xy[1]-1)


class GeometryPrecisionTests(unittest.TestCase):
    def valid(self,function,*args,**kwargs):
        try:
            return function(*args,**kwargs)
        except (DomainError,ArithmeticError) as error:
            self.fail(f'Balanced finite input failed: {type(error).__name__}: {error}')

    def test_area_inverse_distinguishes_the_two_roots_at_one_ulp_above_sonic(self):
        with mp.workdps(100):
            for gamma in (1.1,1.22,1.4,5/3):
                for ratio in (NEXT_AREA,1+8*math.ulp(1.),1+1e-10,1.0001,1.001):
                    for subsonic in (False,True):
                        expected=inverse_area(ratio,gamma,subsonic)
                        actual=self.valid(mach_from_isentropic,'area_ratio_subsonic' if subsonic
                                          else 'area_ratio_supersonic',ratio,gamma)
                        with self.subTest(gamma=gamma,ratio=ratio,subsonic=subsonic):
                            self.assertLessEqual(abs(mp.mpf(actual)-expected),mp.mpf(math.ulp(actual)))
                            self.assertLess(actual,1.) if subsonic else self.assertGreater(actual,1.)

    def test_near_sonic_area_evaluation_matches_an_independent_mass_flux_ratio(self):
        with mp.workdps(150):
            for gamma in (1.1,1.4,1e20,1e100):
                for mach in (math.nextafter(1.,0.),math.nextafter(1.,math.inf),.9999,1.0001):
                    g,m=mp.mpf(gamma),mp.mpf(mach)
                    # Area ratio = sonic mass flux / local mass flux, using
                    # density*Mach*sqrt(T), not the production log expression.
                    t=1/(1+(g-1)*m*m/2);ts=2/(g+1)
                    expected=(ts/t)**(1/(g-1))*mp.sqrt(ts/t)/m
                    actual=calculate_isentropic(mach,gamma).area_critical_ratio
                    with self.subTest(gamma=gamma,mach=mach):
                        self.assertEqual(actual,float(expected))

    def test_pm_evaluation_matches_independent_scaled_differential_quadrature(self):
        with mp.workdps(100):
            for gamma in (math.nextafter(1.,math.inf),1.22,1.4,1e20,1e100,1e200):
                for mach in (math.nextafter(1.,math.inf),1+1e-8,1.0001,1.01):
                    expected=pm_integral(mach,gamma)
                    actual=calculate_isentropic(mach,gamma).prandtl_meyer_deg
                    with self.subTest(gamma=gamma,mach=mach):
                        self.assertGreater(actual,0.)
                        self.assertLess(abs(mp.mpf(actual)/expected-1),mp.mpf('5e-16'))

    def test_small_pm_inverse_keeps_representable_mach_offsets(self):
        with mp.workdps(100):
            for gamma in (1.1,1.22,1.4,1e20,1e100):
                for mach in (math.nextafter(1.,math.inf),1+1e-8,1.0001,2.):
                    target=float(pm_integral(mach,gamma))
                    actual=self.valid(mach_from_isentropic,'prandtl_meyer_deg',target,gamma)
                    with self.subTest(gamma=gamma,mach=mach):
                        self.assertGreater(actual,1.)
                        self.assertLessEqual(abs(actual-mach),math.ulp(mach))

    def test_private_decimal_arithmetic_ignores_hostile_caller_context(self):
        reference=mach_from_isentropic('area_ratio_supersonic',NEXT_AREA,1.22)
        angle=calculate_isentropic(reference,1.22).prandtl_meyer_deg
        with localcontext() as context:
            context.prec=5;context.Emax=9;context.Emin=-9;context.traps[Inexact]=True
            context.traps[FloatOperation]=True
            original=context.copy()
            self.assertEqual(self.valid(mach_from_isentropic,'area_ratio_supersonic',NEXT_AREA,1.22),reference)
            self.assertEqual(self.valid(calculate_isentropic,reference,1.22).prandtl_meyer_deg,angle)
            self.assertEqual(context.flags,original.flags)
            self.assertEqual(context.traps,original.traps)
            self.assertEqual(context.prec,original.prec)

    def test_mutable_decimal_defaults_cannot_configure_the_private_solver(self):
        reference=mach_from_isentropic('area_ratio_supersonic',NEXT_AREA,1.22)
        angle=calculate_isentropic(reference,1.22).prandtl_meyer_deg
        # Initialize/preserve the caller context before changing the process
        # defaults, and restore every global setting even if an assertion fails.
        with localcontext():
            saved=DefaultContext.copy()
            try:
                DefaultContext.prec=5;DefaultContext.rounding=ROUND_DOWN
                DefaultContext.Emin=-9;DefaultContext.Emax=9
                DefaultContext.traps[FloatOperation]=True
                DefaultContext.traps[Inexact]=True
                self.assertEqual(self.valid(mach_from_isentropic,'area_ratio_supersonic',NEXT_AREA,1.22),reference)
                self.assertEqual(self.valid(calculate_isentropic,reference,1.22).prandtl_meyer_deg,angle)
            finally:
                for name in ('prec','rounding','Emin','Emax','capitals','clamp'):
                    setattr(DefaultContext,name,getattr(saved,name))
                DefaultContext.traps=saved.traps.copy()
                DefaultContext.flags=saved.flags.copy()

    def test_public_flow_boundary_owns_reals_and_refuses_malformed_values(self):
        for field in ('mach','gamma'):
            for value in (True,'1',None,math.nan,math.inf,Fraction(1,10**500)):
                with self.subTest(field=field,value=repr(value)[:20]), self.assertRaises(DomainError) as error:
                    calculate_isentropic(**({'mach':2.,'gamma':1.4}|{field:value}))
                self.assertEqual(error.exception.field,field)
        self.assertEqual(calculate_isentropic(Fraction(2),Fraction(7,5)),calculate_isentropic(2.,1.4))

    def test_all_six_finite_counterexamples_produce_resolved_usable_contours(self):
        for generator in GENERATORS:
            for area,ratio in ((.01,NEXT_AREA),(TINY_AREA,10.)):
                with self.subTest(generator=generator.__name__,area=area,ratio=ratio):
                    c=contour(self.valid(geometry,generator,area,ratio))
                    self.assertGreater(c.throat_radius_m,0.)
                    self.assertGreater(c.diverging_length_m,0.)
                    self.assertGreater(c.surface_area_m2,0.)
                    self.assertTrue(all(b.x_m>a.x_m for a,b in zip(c.stations,c.stations[1:])))
                    self.assertTrue(all(p.area_m2>0. and p.area_ratio>=1. for p in c.stations))
                    self.assertGreater(c.stations[-1].mach,1.)
                    snapshot_thermal_contour(c)
                    export=export_nozzle_contour(c,angular_segments=8)
                    self.assertFalse(export.physical_model_qualified)
                    self.assertIn('facet normal',export.stl)
                    json.dumps(asdict(c),allow_nan=False)

    def test_minimum_float_throat_radius_and_area_do_not_erase_private_products(self):
        with mp.workdps(160):
            expected=mp.sqrt(mp.mpf(TINY_AREA)/mp.mpf(math.pi))
            for generator in GENERATORS:
                c=contour(geometry(generator,TINY_AREA,10.))
                with self.subTest(generator=generator.__name__):
                    self.assertEqual(c.throat_radius_m,float(expected))
                    self.assertEqual(min(p.area_m2 for p in c.stations),TINY_AREA)
                    self.assertEqual(min(p.area_ratio for p in c.stations),1.)

    def test_subnormal_interior_area_ratios_are_not_quantized_before_division(self):
        with mp.workdps(160):
            for generator in GENERATORS:
                c=contour(geometry(generator,TINY_AREA,10.))
                for p in c.stations[:-1]:
                    if p.radius_m==c.throat_radius_m:
                        continue  # defining critical boundary, not inferred from rounded rt
                    expected=mp.mpf(math.pi)*mp.mpf(p.radius_m)**2/mp.mpf(TINY_AREA)
                    with self.subTest(generator=generator.__name__,x=p.x_m):
                        self.assertEqual(p.area_ratio,float(expected))

    def test_nominal_boundary_data_and_export_redundancies_remain_coherent_across_scale(self):
        from rocket_propulsion.geometry.export import ContourExportInputs

        rng=random.Random(20261031)
        for _ in range(30):
            area,ratio=10**rng.uniform(-300,300),10**rng.uniform(.01,2)
            for generator in GENERATORS:
                with self.subTest(generator=generator.__name__,area=area,ratio=ratio):
                    c=contour(self.valid(geometry,generator,area,ratio))
                    self.assertEqual(c.stations[-1].area_ratio,ratio)
                    self.valid(ContourExportInputs,c.contour,c.throat_area_m2,c.area_ratio,c.stations)

    def test_tiny_sampled_wall_area_matches_independent_complete_quadrature(self):
        with mp.workdps(160):
            for generator in GENERATORS:
                c=contour(geometry(generator,TINY_AREA,10.))
                expected=mp.fsum(mp.mpf(math.pi)*(mp.mpf(a.radius_m)+mp.mpf(b.radius_m))
                    *mp.sqrt((mp.mpf(b.x_m)-mp.mpf(a.x_m))**2
                             +(mp.mpf(b.radius_m)-mp.mpf(a.radius_m))**2)
                    for a,b in zip(c.stations,c.stations[1:]))
                with self.subTest(generator=generator.__name__):
                    self.assertLessEqual(abs(mp.mpf(c.surface_area_m2)-expected),mp.mpf(math.ulp(c.surface_area_m2)))

    def test_near_unit_conical_and_rao_spans_match_uncancelled_geometry(self):
        with mp.workdps(100):
            for generator in (generate_nozzle_contour,generate_rao_contour):
                c=contour(geometry(generator,.01,NEXT_AREA))
                growth=mp.mpf(c.throat_radius_m)*(mp.sqrt(mp.mpf(NEXT_AREA))-1)
                expected=mp.mpf(.8)*growth/mp.tan(mp.pi/12)
                with self.subTest(generator=generator.__name__):
                    self.assertLess(abs(mp.mpf(c.diverging_length_m)/expected-1),mp.mpf('8e-16'))

    def test_rao_samples_keep_the_control_point_increment_before_adding_the_radius(self):
        with mp.workdps(120):
            for area,ratio in ((.01,NEXT_AREA),(TINY_AREA,10.)):
                c=contour(self.valid(geometry,generate_rao_contour,area,ratio))
                radius=mp.mpf(c.throat_radius_m)
                growth=radius*(mp.sqrt(mp.mpf(ratio))-1)
                length=mp.mpf(.8)*growth/mp.tan(mp.pi/12)
                inlet,exit_=mp.tan(mp.pi/6),mp.tan(mp.pi*8/180)
                control=(growth-exit_*length)/(inlet-exit_)
                for n,p in enumerate(c.stations[15:],1):
                    fraction=mp.mpf(n/26)
                    expected_x=2*(1-fraction)*fraction*control+fraction*fraction*length
                    expected_r=radius+2*(1-fraction)*fraction*inlet*control+fraction*fraction*growth
                    with self.subTest(area=area,sample=n):
                        self.assertLessEqual(abs(mp.mpf(p.x_m)-expected_x),8*mp.mpf(math.ulp(p.x_m)))
                        self.assertLessEqual(abs(mp.mpf(p.radius_m)-expected_r),mp.mpf(math.ulp(p.radius_m)))

    def test_moc_refinement_indicator_is_dimensionless_even_at_minimum_area(self):
        standard=geometry(generate_moc_nozzle,.01,10.)
        tiny=geometry(generate_moc_nozzle,TINY_AREA,10.)
        self.assertGreater(tiny.refined_length_change_relative,1e-5)
        self.assertTrue(math.isclose(tiny.refined_length_change_relative,
                                    standard.refined_length_change_relative,rel_tol=2e-12))

    def test_moc_lengths_match_independent_high_precision_characteristic_intersections(self):
        with mp.workdps(100):
            for area,ratio in ((.01,NEXT_AREA),(.01,1+1e-10),(.01,1.0001),
                               (.01,2.),(TINY_AREA,10.)):
                result=self.valid(geometry,generate_moc_nozzle,area,ratio)
                expected=moc_reference_length(area,ratio,1.22,5)
                # At sonic speed, relative angular/length sensitivity grows
                # with ulp(M)/(M-1). This is a screening budget, not a rigorous
                # interval enclosure; a separate mesh error remains in the DTO.
                sensitivity=math.ulp(result.exit_mach_target)/(result.exit_mach_target-1.)
                tolerance=max(2e-12,32*5*sensitivity)
                with self.subTest(area=area,ratio=ratio):
                    self.assertLess(abs(mp.mpf(result.contour.diverging_length_m)/expected-1),tolerance)

    def test_isentropic_inverse_bracket_remains_a_numerical_limit(self):
        with self.assertRaises(DomainError) as error:
            mach_from_isentropic('area_ratio_supersonic',NEXT_AREA,1e100)
        self.assertEqual(error.exception.code,'outside_inverse_range')

    def test_exhausted_newton_budget_is_not_exported_as_a_converged_inverse(self):
        from rocket_propulsion.compressible.isentropic import _newton_bracketed

        # A smooth, strictly increasing finite callback with the known root
        # x=1e-4. Its steep curvature makes Newton's ~1/600 steps exceed 200
        # iterations from the initial midpoint. Returning that iterate is false.
        with self.assertRaises(ConvergenceError) as error:
            _newton_bracketed(lambda x:math.exp(600*x),lambda x,y:600*y,
                              math.exp(.06),0.,1.,'bounded_exponential')
        self.assertEqual(error.exception.code,'isentropic_inverse_not_converged')


if __name__=='__main__':
    unittest.main()
