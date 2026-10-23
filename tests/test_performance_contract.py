"""Independent performance arithmetic, retained-record and API attacks."""
import json
import unittest
from dataclasses import asdict, replace
from decimal import Decimal, localcontext
from fractions import Fraction
from math import isfinite, log, nextafter, ulp
from random import Random
from sys import float_info
from unittest.mock import patch

from rocket_propulsion.api.routes import calculation_metadata, dispatch_calculation
from rocket_propulsion.core.constants import STANDARD_GRAVITY_M_S2
from rocket_propulsion.core.errors import ConvergenceError, DomainError, InputError
from rocket_propulsion.propulsion import calculate_rocket_equation, calculate_thrust

TINY = ulp(0.)
MAX = float_info.max


def rocket_reference(inputs, precision=1600):
    """High-precision independent equations; no production helper imports."""
    with localcontext() as context:
        context.prec = precision
        values = {name: Decimal.from_float(value) for name,value in inputs.items()}
        g = Decimal.from_float(STANDARD_GRAVITY_M_S2)
        m = values['initial_mass_kg']
        if 'specific_impulse_s' in values:
            isp = values['specific_impulse_s']; velocity = isp*g
        else:
            velocity = values['effective_exhaust_velocity_m_s']; isp = velocity/g
        if 'final_mass_kg' in values:
            final = values['final_mass_kg']; ratio = m/final
            consumed = m-final; fraction = consumed/m; delta = velocity*ratio.ln()
        else:
            delta = values['delta_v_m_s']; x = delta/velocity
            ratio = x.exp(); survival = (-x).exp()
            final = m*survival; fraction = 1-survival; consumed = m*fraction
        return {name: float(value) for name,value in dict(initial_mass_kg=m,
            final_mass_kg=final,propellant_mass_kg=consumed,mass_ratio=ratio,
            propellant_mass_fraction=fraction,specific_impulse_s=isp,
            effective_exhaust_velocity_m_s=velocity,delta_v_m_s=delta).items()}


def thrust_reference(inputs):
    x = {name: Fraction(value) for name,value in inputs.items()}
    momentum = x['mass_flow_kg_s']*x['exit_velocity_m_s']
    pressure = x['exit_area_m2']*(x['exit_pressure_pa']-x['ambient_pressure_pa'])
    total = momentum+pressure
    return dict(momentum_thrust_n=float(momentum),pressure_thrust_n=float(pressure),
                total_thrust_n=float(total),effective_exhaust_velocity_m_s=float(total/x['mass_flow_kg_s']),
                specific_impulse_s=float(total/(x['mass_flow_kg_s']*Fraction(STANDARD_GRAVITY_M_S2))),
                total_impulse_n_s=float(total*x['burn_time_s']) if 'burn_time_s' in x else None)


class PerformanceAttacks(unittest.TestCase):
    def valid(self, function, **inputs):
        try:
            return function(**inputs)
        except (DomainError, ConvergenceError, ArithmeticError) as error:
            self.fail(f'Representable valid performance case rejected: {error}')

    def check_rocket(self, inputs, precision=1600):
        expected=rocket_reference(inputs,precision)
        self.assertEqual(expected,rocket_reference(inputs,precision+400))
        self.assertTrue(all(isfinite(x) and x>0 for x in expected.values()))
        result=self.valid(calculate_rocket_equation,**inputs)
        self.assertEqual({name:getattr(result,name) for name in expected},expected)
        return result

    def check_thrust(self, inputs):
        expected=thrust_reference(inputs)
        result=self.valid(calculate_thrust,**inputs)
        self.assertEqual({name:getattr(result,name) for name in expected},expected)
        return result

    def thrust_inputs(self):
        return dict(mass_flow_kg_s=10.,exit_velocity_m_s=3000.,exit_area_m2=.5,
                    exit_pressure_pa=50000.,ambient_pressure_pa=10000.,burn_time_s=20.)

    def test_near_equal_large_and_small_masses_preserve_delta_v(self):
        for initial in (1e308,1.,1e-200):
            with self.subTest(initial=initial):
                self.check_rocket(dict(initial_mass_kg=initial,final_mass_kg=nextafter(initial,0.),
                                       specific_impulse_s=300.))

    def test_forward_outputs_use_unrounded_performance_and_mass_ledgers(self):
        for inputs in (dict(initial_mass_kg=500000.,final_mass_kg=100000.,specific_impulse_s=450.),
                       dict(initial_mass_kg=3.,final_mass_kg=.7,effective_exhaust_velocity_m_s=nextafter(10.,20.)),
                       dict(initial_mass_kg=1e-100,final_mass_kg=1e-110,specific_impulse_s=1e100),
                       dict(initial_mass_kg=1e200,final_mass_kg=1e199,specific_impulse_s=1e-200)):
            with self.subTest(inputs=inputs):self.check_rocket(inputs)

    def test_forward_power_of_two_mass_scaling_does_not_change_delta_v(self):
        results=[]
        for exponent in (-600,0,600):
            scale=2.**exponent
            results.append(self.check_rocket(dict(initial_mass_kg=7.*scale,final_mass_kg=3.*scale,
                                                  effective_exhaust_velocity_m_s=1234.5)))
        self.assertEqual(len({r.delta_v_m_s for r in results}),1)
        self.assertEqual(len({r.propellant_mass_fraction for r in results}),1)

    def test_seeded_forward_cases_match_stable_independent_logarithms(self):
        rng=Random(221)
        for _ in range(8):
            initial=rng.uniform(1.,10.)*2.**rng.randint(-300,300)
            inputs=dict(initial_mass_kg=initial,final_mass_kg=initial*rng.uniform(.1,.99),
                        specific_impulse_s=rng.uniform(50.,500.))
            self.check_rocket(inputs,600)

    def test_inverse_outputs_are_independent_exports_of_the_same_exponent(self):
        for inputs in (dict(initial_mass_kg=500000.,delta_v_m_s=7100.,specific_impulse_s=450.),
                       dict(initial_mass_kg=1e308,delta_v_m_s=TINY,effective_exhaust_velocity_m_s=1.),
                       dict(initial_mass_kg=1e100,delta_v_m_s=1e-200,effective_exhaust_velocity_m_s=1e100),
                       dict(initial_mass_kg=1.,delta_v_m_s=700.,effective_exhaust_velocity_m_s=1.),
                       dict(initial_mass_kg=1e300,delta_v_m_s=1e-20,specific_impulse_s=300.)):
            with self.subTest(inputs=inputs):self.check_rocket(inputs,2200)

    def test_inverse_rounded_mass_equality_does_not_erase_propellant(self):
        result=self.check_rocket(dict(initial_mass_kg=1e308,delta_v_m_s=TINY,
                                      effective_exhaust_velocity_m_s=1.),2200)
        self.assertEqual(result.final_mass_kg,result.initial_mass_kg)
        self.assertEqual(result.mass_ratio,1.)
        self.assertGreater(result.propellant_mass_kg,0.)
        self.assertEqual(result.propellant_mass_fraction,TINY)

    def test_inverse_mass_ratio_overflow_boundary_is_not_an_arbitrary_input_cap(self):
        delta=log(MAX)
        self.check_rocket(dict(initial_mass_kg=MAX,delta_v_m_s=delta,effective_exhaust_velocity_m_s=1.),2200)
        with self.assertRaises(DomainError) as caught:
            calculate_rocket_equation(initial_mass_kg=MAX,delta_v_m_s=nextafter(delta,float('inf')),
                                      effective_exhaust_velocity_m_s=1.)
        self.assertEqual(caught.exception.field,'Mass ratio')

    def test_inverse_final_mass_underflow_boundary_preserves_the_last_nonzero_cell(self):
        delta=log(6.)
        self.check_rocket(dict(initial_mass_kg=3*TINY,delta_v_m_s=nextafter(delta,0.),
                              effective_exhaust_velocity_m_s=1.),2200)
        with self.assertRaises(DomainError) as caught:
            calculate_rocket_equation(initial_mass_kg=3*TINY,delta_v_m_s=nextafter(delta,float('inf')),
                                      effective_exhaust_velocity_m_s=1.)
        self.assertEqual(caught.exception.field,'Final mass')

    def test_forward_and_inverse_keep_their_distinct_source_mode(self):
        forward=self.check_rocket(dict(initial_mass_kg=500000.,final_mass_kg=100000.,specific_impulse_s=450.))
        inverse=self.check_rocket(dict(initial_mass_kg=500000.,delta_v_m_s=forward.delta_v_m_s,specific_impulse_s=450.))
        self.assertEqual((forward.solve_mode,inverse.solve_mode),('forward','inverse'))
        self.assertEqual(forward.performance_input,'specific_impulse_s')
        self.assertNotEqual(forward.numerical_policy,'')

    def test_all_public_rocket_quantities_reject_malformed_values_typed(self):
        forward=dict(initial_mass_kg=2.,final_mass_kg=1.,specific_impulse_s=1.)
        inverse=dict(initial_mass_kg=2.,delta_v_m_s=1.,effective_exhaust_velocity_m_s=1.)
        for inputs in (forward,inverse):
            for field in inputs:
                for bad in (True,'1',[],float('nan'),float('inf'),10**400,0.,-1.):
                    with self.subTest(field=field,bad=repr(bad)),self.assertRaises(DomainError):
                        calculate_rocket_equation(**{**inputs,field:bad})

    def test_mutually_exclusive_or_missing_rocket_inputs_fail_closed(self):
        for inputs in ({},dict(initial_mass_kg=2.,final_mass_kg=1.),
                       dict(initial_mass_kg=2.,specific_impulse_s=300.),
                       dict(initial_mass_kg=2.,final_mass_kg=1.,delta_v_m_s=1.,specific_impulse_s=300.),
                       dict(initial_mass_kg=2.,final_mass_kg=1.,specific_impulse_s=300.,effective_exhaust_velocity_m_s=1.),
                       dict(initial_mass_kg=1.,final_mass_kg=1.,specific_impulse_s=300.),
                       dict(initial_mass_kg=1.,final_mass_kg=2.,specific_impulse_s=300.)):
            with self.subTest(inputs=inputs),self.assertRaises(DomainError):calculate_rocket_equation(**inputs)

    def test_unrepresentable_rocket_outputs_are_named_domain_errors(self):
        for inputs,field in ((dict(initial_mass_kg=MAX,final_mass_kg=TINY,specific_impulse_s=300.),'Mass ratio'),
                             (dict(initial_mass_kg=100.,delta_v_m_s=1e308,specific_impulse_s=1.),'Mass ratio'),
                             (dict(initial_mass_kg=2.,final_mass_kg=1.,specific_impulse_s=MAX),'Effective exhaust velocity'),
                             (dict(initial_mass_kg=2.,final_mass_kg=1.,effective_exhaust_velocity_m_s=TINY),'Specific impulse'),
                             (dict(initial_mass_kg=TINY,delta_v_m_s=1.,effective_exhaust_velocity_m_s=1.),'Final mass')):
            with self.subTest(inputs=inputs),self.assertRaises(DomainError) as caught:
                calculate_rocket_equation(**inputs)
            self.assertEqual(caught.exception.field,field)

    def test_thrust_recovers_positive_cancellation_of_rounded_products(self):
        inputs=dict(mass_flow_kg_s=nextafter(1.,2.),exit_velocity_m_s=nextafter(1.,2.),exit_area_m2=1.,
                    exit_pressure_pa=1.,ambient_pressure_pa=nextafter(2.,3.),burn_time_s=1e100)
        result=self.check_thrust(inputs)
        self.assertEqual(result.total_thrust_n,2.**-104)
        self.assertGreater(result.total_impulse_n_s,0.)

    def test_thrust_derived_velocity_isp_and_impulse_do_not_reuse_rounded_force(self):
        for inputs in (self.thrust_inputs(),dict(mass_flow_kg_s=.7,exit_velocity_m_s=1.3,exit_area_m2=.2,
                       exit_pressure_pa=3.1,ambient_pressure_pa=.9,burn_time_s=1.7),
                       dict(mass_flow_kg_s=1e-150,exit_velocity_m_s=1e150,exit_area_m2=1e-150,
                            exit_pressure_pa=1e150,ambient_pressure_pa=0.,burn_time_s=1e100),
                       dict(mass_flow_kg_s=1.9278895965686447,exit_velocity_m_s=6.842186066281001,
                            exit_area_m2=9.564215230275066,exit_pressure_pa=5.614348552051298,
                            ambient_pressure_pa=4.490804331401026,burn_time_s=9.48046974867356),
                       dict(mass_flow_kg_s=7.517394510277857,exit_velocity_m_s=6.549235479927716,
                            exit_area_m2=.6483561707626898,exit_pressure_pa=8.295362635590601,
                            ambient_pressure_pa=.5978115962349206,burn_time_s=4.911995580020825)):
            with self.subTest(inputs=inputs):self.check_thrust(inputs)

    def test_zero_pressure_thrust_and_optional_impulse_are_valid(self):
        inputs=self.thrust_inputs();inputs.pop('burn_time_s');inputs['ambient_pressure_pa']=inputs['exit_pressure_pa']
        result=self.check_thrust(inputs)
        self.assertEqual(result.pressure_thrust_n,0.)
        self.assertIsNone(result.total_impulse_n_s)
        self.assertIsNone(result.burn_time_s)

    def test_negative_pressure_term_is_allowed_only_with_positive_net_force(self):
        inputs=self.thrust_inputs();inputs['ambient_pressure_pa']=60000.
        result=self.check_thrust(inputs)
        self.assertLess(result.pressure_thrust_n,0.)
        for pressure in (110000.,120000.):
            with self.subTest(pressure=pressure),self.assertRaises(DomainError):
                calculate_thrust(**{**inputs,'ambient_pressure_pa':pressure})

    def test_all_thrust_inputs_have_strict_typed_boundaries(self):
        inputs=self.thrust_inputs()
        for field in inputs:
            for bad in (True,'1',None,[],float('nan'),float('inf'),10**400,-1.):
                if field=='burn_time_s' and bad is None:continue
                with self.subTest(field=field,bad=repr(bad)),self.assertRaises(DomainError):
                    calculate_thrust(**{**inputs,field:bad})

    def test_unrepresentable_thrust_outputs_are_not_clipped_or_zeroed(self):
        for inputs,field in ((dict(mass_flow_kg_s=MAX,exit_velocity_m_s=1.,exit_area_m2=1.,exit_pressure_pa=MAX,ambient_pressure_pa=0.),'Total thrust'),
                             (dict(mass_flow_kg_s=TINY,exit_velocity_m_s=.25,exit_area_m2=1.,exit_pressure_pa=2.,ambient_pressure_pa=1.),'Momentum thrust'),
                             (dict(mass_flow_kg_s=TINY,exit_velocity_m_s=1.,exit_area_m2=1.,exit_pressure_pa=1.,ambient_pressure_pa=1.,burn_time_s=.25),'Total impulse')):
            with self.subTest(field=field),self.assertRaises(DomainError) as caught:calculate_thrust(**inputs)
            self.assertEqual(caught.exception.field,field)

    def test_forward_result_rejects_every_forged_derived_field_and_policy(self):
        result=self.valid(calculate_rocket_equation,initial_mass_kg=500000.,final_mass_kg=100000.,specific_impulse_s=450.)
        for field in ('propellant_mass_kg','mass_ratio','propellant_mass_fraction','effective_exhaust_velocity_m_s','delta_v_m_s'):
            with self.subTest(field=field),self.assertRaises(DomainError):replace(result,**{field:nextafter(getattr(result,field),float('inf'))})
        for change in ({'numerical_policy':'unknown'},{'solve_mode':'unknown'},{'performance_input':'unknown'},
                       {'final_mass_kg':result.initial_mass_kg},{'specific_impulse_s':True},{'delta_v_m_s':float('inf')}):
            with self.subTest(change=change),self.assertRaises(DomainError):replace(result,**change)

    def test_inverse_result_checks_original_delta_v_instead_of_rounded_mass_loss(self):
        result=self.valid(calculate_rocket_equation,initial_mass_kg=1e308,delta_v_m_s=TINY,effective_exhaust_velocity_m_s=1.)
        for change in ({'propellant_mass_kg':0.},{'propellant_mass_fraction':0.},
                       {'final_mass_kg':nextafter(result.final_mass_kg,0.)},{'mass_ratio':2.},
                       {'specific_impulse_s':nextafter(result.specific_impulse_s,0.)}):
            with self.subTest(change=change),self.assertRaises(DomainError):replace(result,**change)

    def test_thrust_result_replays_every_derived_field_and_duration(self):
        result=self.check_thrust(self.thrust_inputs())
        for field in thrust_reference(self.thrust_inputs()):
            with self.subTest(field=field),self.assertRaises(DomainError):replace(result,**{field:nextafter(getattr(result,field),float('inf'))})
        for change in ({'burn_time_s':None},{'burn_time_s':True},{'numerical_policy':'unknown'},{'mass_flow_kg_s':True},
                       {'ambient_pressure_pa':float('inf')}):
            with self.subTest(change=change),self.assertRaises(DomainError):replace(result,**change)

    def test_results_round_trip_as_finite_json_and_canonical_real_scalars(self):
        rocket=self.valid(calculate_rocket_equation,initial_mass_kg=2,final_mass_kg=1,specific_impulse_s=300)
        thrust=self.check_thrust(self.thrust_inputs())
        for result in (rocket,thrust):
            document=json.loads(json.dumps(asdict(result),allow_nan=False))
            self.assertEqual(type(result)(**document),result)
            self.assertEqual(document['numerical_policy'],result.numerical_policy)
        self.assertIs(type(rocket.initial_mass_kg),float)

    def test_decimal_context_cannot_change_performance_outputs(self):
        inputs=dict(initial_mass_kg=3.,final_mass_kg=1.1,specific_impulse_s=321.)
        reference=rocket_reference(inputs)
        with localcontext() as context:
            context.prec=3
            result=self.valid(calculate_rocket_equation,**inputs)
        self.assertEqual(result.delta_v_m_s,reference['delta_v_m_s'])

    def test_unresolved_log_enclosure_is_not_reported_as_unrepresentability(self):
        from rocket_propulsion.propulsion import performance_numerics as numerics
        with patch.object(numerics,'_log_interval',return_value=(Decimal(0),Decimal(2))):
            with self.assertRaises(ConvergenceError):
                calculate_rocket_equation(initial_mass_kg=2.,final_mass_kg=1.,effective_exhaust_velocity_m_s=1.)

    def test_intervals_straddling_export_boundaries_remain_unresolved(self):
        from rocket_propulsion.propulsion import performance_numerics as numerics
        for bounds in ((Decimal(1),Decimal(2)),(Decimal.from_float(MAX),Decimal(2**1024)),
                       (Decimal(0),Decimal.from_float(TINY))):
            with self.subTest(bounds=bounds):
                try:numerics._positive_enclosed_output(lambda precision:bounds,'Probe output')
                except ConvergenceError:pass
                except DomainError as error:self.fail(f'Unresolved bounds mislabelled as impossible: {error}')
                else:self.fail('Unresolved bounds were accepted as a proven output.')

    def test_enclosure_export_keeps_last_finite_and_subnormal_rounding_cells(self):
        from rocket_propulsion.propulsion import performance_numerics as numerics
        maximum_bounds=(Decimal(int(MAX)+2**968),Decimal(int(MAX)+2**969))
        with localcontext() as context:
            context.prec=2200
            tiny_bounds=(Decimal(3)/Decimal(2**1076),Decimal(7)/Decimal(2**1077))
        for bounds,expected in ((maximum_bounds,MAX),(tiny_bounds,TINY)):
            with self.subTest(expected=expected):
                actual=self.valid(numerics._positive_enclosed_output,evaluate=lambda precision:bounds,label='Probe output')
                self.assertEqual(actual,expected)


class PerformanceApiAttacks(unittest.TestCase):
    def test_api_recovers_near_equal_masses_and_positive_thrust_cancellation(self):
        inputs=dict(initial_mass_kg=1e308,final_mass_kg=nextafter(1e308,0.),specific_impulse_s=300.)
        result=dispatch_calculation('/api/v1/rocket-equation',inputs)
        self.assertEqual(result['delta_v_m_s'],rocket_reference(inputs)['delta_v_m_s'])
        thrust=dict(mass_flow_kg_s=nextafter(1.,2.),exit_velocity_m_s=nextafter(1.,2.),exit_area_m2=1.,
                    exit_pressure_pa=1.,ambient_pressure_pa=nextafter(2.,3.))
        try:result=dispatch_calculation('/api/v1/thrust',thrust)
        except DomainError as error:self.fail(str(error))
        self.assertEqual(result['total_thrust_n'],thrust_reference(thrust)['total_thrust_n'])
        json.dumps(result,allow_nan=False)

    def test_api_unrepresentable_outputs_remain_typed(self):
        for path,inputs in (('/api/v1/rocket-equation',dict(initial_mass_kg=100.,delta_v_m_s=1e308,specific_impulse_s=1.)),
                            ('/api/v1/thrust',dict(mass_flow_kg_s=MAX,exit_velocity_m_s=MAX,exit_area_m2=1.,exit_pressure_pa=1.))):
            with self.subTest(path=path),self.assertRaises((DomainError,InputError)):dispatch_calculation(path,inputs)

    def test_api_metadata_identifies_revised_arithmetic_and_retained_inputs(self):
        for path,inputs in (('/api/v1/rocket-equation',dict(initial_mass_kg=2.,final_mass_kg=1.,specific_impulse_s=300.)),
                            ('/api/v1/thrust',dict(mass_flow_kg_s=1.,exit_velocity_m_s=1.,exit_area_m2=1.,exit_pressure_pa=1.,burn_time_s=2.))):
            with self.subTest(path=path):
                self.assertEqual(calculation_metadata(path)['version'],'2.0')
                result=dispatch_calculation(path,inputs)
                self.assertIn('numerical_policy',result)
                if path.endswith('thrust'):self.assertEqual(result['burn_time_s'],2.)
