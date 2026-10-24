"""Exact energy/laminar ledgers, independent turbulent references and hostile records."""
from collections.abc import Mapping
from dataclasses import asdict, replace
from decimal import Decimal, localcontext
from fractions import Fraction
from math import inf, nextafter
import unittest

from rocket_propulsion.core.errors import DomainError, ConvergenceError
from rocket_propulsion.thermodynamics import calculate_regenerative_cooling
from rocket_propulsion.api import routes as api_routes


def inputs(**changes):
    return dict(dict(absorbed_heat_w=1000.,coolant_mass_flow_kg_s=1.,coolant_specific_heat_j_kg_k=1000.,
        coolant_inlet_temperature_k=100.,coolant_density_kg_m3=1000.,coolant_dynamic_viscosity_pa_s=.001,
        channel_count=10,channel_flow_area_m2=.001,hydraulic_diameter_m=.01,channel_length_m=1.,
        roughness_m=0.,maximum_coolant_temperature_k=1000.),**changes)


def rational_state(x):
    f=Fraction
    capacity=f(x['coolant_mass_flow_kg_s'])*f(x['coolant_specific_heat_j_kg_k'])
    rise=f(x['absorbed_heat_w'])/capacity
    outlet=f(x['coolant_inlet_temperature_k'])+rise
    velocity=f(x['coolant_mass_flow_kg_s'])/(f(x['coolant_density_kg_m3'])*x['channel_count']*f(x['channel_flow_area_m2']))
    reynolds=f(x['coolant_mass_flow_kg_s'])*f(x['hydraulic_diameter_m'])/(x['channel_count']*f(x['channel_flow_area_m2'])*f(x['coolant_dynamic_viscosity_pa_s']))
    head=f(x['channel_length_m'])/f(x['hydraulic_diameter_m'])*f(x['coolant_density_kg_m3'])*velocity**2/2
    return capacity,rise,outlet,velocity,reynolds,head


def optional_float(value):
    try:r=float(value)
    except OverflowError:return None
    return None if value and r==0 else r


def reference(x,precision=180):
    capacity,rise,outlet,velocity,reynolds,head=rational_state(x)
    with localcontext() as c:
        c.prec=precision
        d=lambda f:Decimal(f.numerator)/Decimal(f.denominator)
        if reynolds<2300:
            friction=Fraction(x.get('laminar_darcy_poiseuille_number',64.))/reynolds
            exported_f=float(friction);pressure=float(friction*head)
        else:
            arg=d(Fraction(x['roughness_m'])/(Fraction(37,10)*Fraction(x['hydraulic_diameter_m'])))
            arg+=Decimal('5.74')/((d(reynolds).ln()*Decimal('.9')).exp())
            friction=Decimal('.25')/(arg.ln()/Decimal(10).ln())**2
            exported_f=float(friction);pressure=float(friction*d(head))
    margin=Fraction(x['maximum_coolant_temperature_k'])-outlet
    displayed=abs(Fraction(x['absorbed_heat_w'])-capacity*(Fraction(float(outlet))-Fraction(x['coolant_inlet_temperature_k'])))
    return dict(coolant_capacity_w_k=float(capacity),temperature_rise_k=float(rise),outlet_temperature_k=float(outlet),
        inlet_temperature_k=x['coolant_inlet_temperature_k'],absorbed_heat_w=x['absorbed_heat_w'],
        channel_velocity_m_s=float(velocity),reynolds_number=float(reynolds),darcy_friction_factor=exported_f,
        pressure_drop_pa=pressure,temperature_limit_exceeded=margin<0,temperature_limit_margin_k=optional_float(margin),
        displayed_temperature_energy_residual_w=optional_float(displayed),
        flow_regime='laminar' if reynolds<2300 else 'transition' if reynolds<4000 else 'turbulent')


def pole_inputs(multiplier=1,near=False):
    root=4;power=root**9;diameter=500*power
    roughness=37*(50*power*multiplier-287)
    flow=float(root**10)
    return inputs(channel_count=1,channel_flow_area_m2=float(diameter),hydraulic_diameter_m=float(diameter),
        roughness_m=float(roughness),coolant_dynamic_viscosity_pa_s=1.,coolant_density_kg_m3=1.,
        coolant_mass_flow_kg_s=nextafter(flow,inf) if near else flow)


class RegenerativeAttacks(unittest.TestCase):
    def valid(self,x):
        try:return calculate_regenerative_cooling(**x)
        except (DomainError,ConvergenceError,ArithmeticError) as error:
            self.fail(f'Representable prescribed coolant witness rejected: {error}')

    def check(self,x):
        expected=reference(x);self.assertEqual(expected,reference(x,360),'Independent reference did not stabilize')
        r=self.valid(x)
        for name,value in expected.items():self.assertEqual(getattr(r,name),value,name)
        self.assertFalse(r.physical_model_qualified)
        self.assertFalse(r.friction_correlation_domain_qualified)
        self.assertEqual(r.energy_balance_residual_w,0.)
        return r

    def test_boolean_heat_has_a_valid_numeric_positive_control(self):
        self.valid(inputs(absorbed_heat_w=1.))
        with self.assertRaises(DomainError):calculate_regenerative_cooling(**inputs(absorbed_heat_w=True))

    def test_quantity_fields_reject_bool_text_none_nonfinite_and_erased_values(self):
        for name in inputs():
            if name=='channel_count':continue
            for value in (True,None,'1',inf,float('nan'),Fraction(1,10**1000)):
                with self.subTest(field=name,value=repr(value)),self.assertRaises(DomainError):
                    calculate_regenerative_cooling(**inputs(**{name:value}))

    def test_channel_count_is_a_positive_actual_integer(self):
        for value in (True,1.,None,'1',0,-1):
            with self.subTest(value=value),self.assertRaises(DomainError):
                calculate_regenerative_cooling(**inputs(channel_count=value))

    def test_laminar_profile_matches_independent_exact_rationals(self):self.check(inputs())

    def test_turbulent_profile_matches_independent_decimal_reference(self):
        self.check(inputs(coolant_mass_flow_kg_s=100.,roughness_m=1e-5))

    def test_recoverable_private_density_area_overflow(self):
        r=self.check(inputs(absorbed_heat_w=1.,coolant_mass_flow_kg_s=1e200,coolant_specific_heat_j_kg_k=1e-200,
            coolant_density_kg_m3=1e200,coolant_dynamic_viscosity_pa_s=1.,channel_count=1,channel_flow_area_m2=1e200,hydraulic_diameter_m=1.))
        self.assertEqual(r.pressure_drop_pa,3.2e-199)

    def test_unrepresentable_mandatory_capacity_fails_typed(self):
        with self.assertRaises(DomainError):calculate_regenerative_cooling(**inputs(coolant_mass_flow_kg_s=10.,coolant_specific_heat_j_kg_k=1e308))

    def test_true_limit_exceedance_is_not_hidden_by_display_rounding(self):
        x=inputs(absorbed_heat_w=nextafter(2.**-52,inf),coolant_mass_flow_kg_s=1.,coolant_specific_heat_j_kg_k=1.,
                 coolant_inlet_temperature_k=1.,maximum_coolant_temperature_k=nextafter(1.,inf))
        r=self.check(x)
        self.assertEqual(r.outlet_temperature_k,x['maximum_coolant_temperature_k'])
        self.assertTrue(r.temperature_limit_exceeded)
        self.assertTrue(any(w.code=='coolant_temperature_limit_exceeded' and w.severity=='critical' for w in r.warnings))

    def test_exact_limit_equality_is_not_reported_as_exceeded(self):
        r=self.check(inputs(absorbed_heat_w=1.,coolant_mass_flow_kg_s=1.,coolant_specific_heat_j_kg_k=1.,
            coolant_inlet_temperature_k=1.,maximum_coolant_temperature_k=2.))
        self.assertFalse(r.temperature_limit_exceeded)
        self.assertEqual(r.temperature_limit_margin_k,0.)
        self.assertFalse(any(w.code=='coolant_temperature_limit_exceeded' for w in r.warnings))

    def test_tiny_optional_margin_and_display_residual_do_not_invalidate_representable_outputs(self):
        tiny=float.fromhex('0x0.0000000000001p-1022')
        r=self.check(inputs(absorbed_heat_w=tiny,coolant_mass_flow_kg_s=1.,coolant_specific_heat_j_kg_k=nextafter(1.,inf),
            coolant_inlet_temperature_k=tiny,maximum_coolant_temperature_k=2*tiny))
        self.assertIsNone(r.temperature_limit_margin_k)
        self.assertIsNone(r.displayed_temperature_energy_residual_w)
        self.assertFalse(r.temperature_limit_exceeded)

    def test_positive_subnormal_heat_is_not_erased_when_outlet_does_not_move(self):
        tiny=float.fromhex('0x0.0000000000001p-1022')
        r=self.check(inputs(absorbed_heat_w=tiny,coolant_mass_flow_kg_s=1.,coolant_specific_heat_j_kg_k=1.,coolant_inlet_temperature_k=1.))
        self.assertEqual(r.temperature_rise_k,tiny);self.assertEqual(r.outlet_temperature_k,1.)
        self.assertEqual(r.displayed_temperature_energy_residual_w,tiny)

    def test_exact_reynolds_below_threshold_does_not_branch_on_rounded_export(self):
        delta=2.**-52
        for threshold,regime in ((2300.,'laminar'),(4000.,'transition')):
            x=inputs(channel_count=1,coolant_mass_flow_kg_s=threshold,hydraulic_diameter_m=1+2*delta,
                channel_flow_area_m2=1+delta,coolant_dynamic_viscosity_pa_s=1+delta)
            self.assertLess(rational_state(x)[4],Fraction(threshold))
            r=self.check(x);self.assertEqual(r.reynolds_number,threshold);self.assertEqual(r.flow_regime,regime)

    def test_non_circular_laminar_factor_is_an_explicit_input(self):
        x=inputs(laminar_darcy_poiseuille_number=96.)
        r=self.check(x);base=self.check(inputs())
        self.assertAlmostEqual(r.darcy_friction_factor/base.darcy_friction_factor,1.5)
        self.assertAlmostEqual(r.pressure_drop_pa/base.pressure_drop_pa,1.5)

    def test_density_changes_velocity_and_head_but_not_reynolds(self):
        a=self.check(inputs());b=self.check(inputs(coolant_density_kg_m3=2000.))
        self.assertEqual(a.reynolds_number,b.reynolds_number)
        self.assertEqual(a.channel_velocity_m_s/2,b.channel_velocity_m_s)
        self.assertEqual(a.pressure_drop_pa/2,b.pressure_drop_pa)

    def test_parallel_channels_use_per_channel_flow(self):
        a=self.check(inputs());b=self.check(inputs(channel_count=20))
        self.assertEqual(a.reynolds_number/2,b.reynolds_number)
        self.assertEqual(a.pressure_drop_pa/2,b.pressure_drop_pa)

    def test_large_exact_channel_count_does_not_overflow_before_small_area_product(self):
        self.check(inputs(channel_count=2**1024,channel_flow_area_m2=float.fromhex('0x0.0000000000001p-1022')))

    def test_laminar_coefficient_is_validated_even_in_turbulent_branch(self):
        for value in (True,None,'64',0.,-1.,inf):
            with self.subTest(value=value),self.assertRaises(DomainError):
                calculate_regenerative_cooling(**inputs(coolant_mass_flow_kg_s=100.,laminar_darcy_poiseuille_number=value))

    def test_pressure_drop_does_not_reuse_rounded_friction_factor(self):
        x=inputs(coolant_mass_flow_kg_s=293.1494670180798,channel_length_m=1.6575164190592035,roughness_m=1.719533884578709e-05)
        r=self.check(x)
        self.assertEqual(r.pressure_drop_pa,1662380.5401107145)
        self.assertEqual(float(Fraction(r.darcy_friction_factor)*rational_state(x)[5]),1662380.5401107143)

    def test_temperature_rise_does_not_reuse_rounded_heat_capacity(self):
        x=inputs(absorbed_heat_w=9.030606519217852,coolant_mass_flow_kg_s=1.447258791541652,
                 coolant_specific_heat_j_kg_k=2.9972871869898663,coolant_inlet_temperature_k=.1)
        r=self.check(x)
        self.assertEqual(r.temperature_rise_k,2.08181611534103)
        self.assertNotEqual(r.temperature_rise_k,float(Fraction(x['absorbed_heat_w'])/Fraction(r.coolant_capacity_w_k)))

    def test_exact_swamee_jain_pole_fails_with_a_typed_singularity(self):
        with self.assertRaises(DomainError) as caught:calculate_regenerative_cooling(**pole_inputs())
        self.assertEqual(caught.exception.code,'singular_swamee_jain_argument')

    def test_representable_near_pole_is_not_collapsed_to_log_zero(self):
        self.check(pole_inputs(near=True))

    def test_positive_log_extension_is_flagged_instead_of_qualified(self):
        r=self.check(pole_inputs(multiplier=10))
        self.assertEqual(r.darcy_friction_factor,.25)
        self.assertTrue(any(w.code=='swamee_jain_nonphysical_inverse_root' and w.severity=='critical' for w in r.warnings))

    def test_temperature_display_discrepancy_is_not_measured_energy_balance(self):
        r=self.check(inputs(absorbed_heat_w=1.,coolant_inlet_temperature_k=1e200,maximum_coolant_temperature_k=1e201))
        self.assertEqual(r.outlet_temperature_k,1e200)
        self.assertEqual(r.displayed_temperature_energy_residual_w,1.)
        self.assertEqual(r.residual_kind,'bulk_energy_ledger_identity_not_measured_energy_closure')

    def test_complete_record_replay_rejects_forged_fields_and_qualifications(self):
        r=self.valid(inputs())
        for change in (dict(outlet_temperature_k=1e99),dict(pressure_drop_pa=1e99),dict(physical_model_qualified=True),
                       dict(physical_model_qualified=0),dict(temperature_limit_exceeded=True),dict(numerical_policy='forged'),
                       dict(flow_regime='turbulent'),dict(energy_balance_residual_w=1.)):
            with self.subTest(change=change),self.assertRaises(DomainError):replace(r,**change)

    def test_mapping_round_trip_owns_inputs_and_warning_text(self):
        r=self.valid(inputs());raw=asdict(r);restored=type(r)(**raw)
        self.assertEqual(restored,r)
        raw['inputs']['coolant_density_kg_m3']=1.
        self.assertEqual(restored.inputs.coolant_density_kg_m3,1000.)
        class EqualText(str):
            def __eq__(self,other):return True
            def __ne__(self,other):return False
        with self.assertRaises(DomainError):replace(r,numerical_policy=EqualText('forged'))
        warning=replace(r.warnings[0],message=EqualText('forged'))
        with self.assertRaises(DomainError):replace(r,warnings=(warning,)+r.warnings[1:])

    def test_known_mapping_failures_are_typed_without_swallowing_programmer_errors(self):
        class Broken(Mapping):
            def __init__(self,error):self.error=error
            def __len__(self):return 1
            def __getitem__(self,key):raise self.error('callback')
            def __iter__(self):raise self.error('callback')
        r=self.valid(inputs())
        with self.assertRaises(DomainError):replace(r,inputs=Broken(KeyError))
        with self.assertRaises(RuntimeError):replace(r,inputs=Broken(RuntimeError))

    def test_warning_iterable_cardinality_is_bounded(self):
        r=self.valid(inputs())
        def endless():
            for _ in range(7):yield r.warnings[0]
            raise AssertionError('Read beyond coolant warning cardinality rejection')
        with self.assertRaises(DomainError):replace(r,warnings=endless())


class RegenerativeApiAttacks(unittest.TestCase):
    def payload(self):
        x=inputs();x.pop('absorbed_heat_w')
        x.update(throat_area_m2=.01,area_ratio=25.,station_count=41,chamber_pressure_pa=5e6,chamber_temperature_k=3500.,
            wall_temperature_k=800.,characteristic_velocity_m_s=1700.,dynamic_viscosity_pa_s=7e-5,
            specific_heat_j_kg_k=3000.,prandtl_number=.7,gamma=1.22,throat_radius_of_curvature_m=.03,include_cooling=True)
        return x

    def test_api_forwards_explicit_laminar_shape_factor(self):
        x=self.payload();x['laminar_darcy_poiseuille_number']='96'
        for invoke in (lambda:api_routes.dispatch_calculation('/api/v1/nozzle-thermal',x),
                       lambda:api_routes.calculate_nozzle_thermal_route(x)):
            with self.subTest(path=invoke):
                r=invoke()['cooling']
                self.assertEqual(r['inputs']['laminar_darcy_poiseuille_number'],96.)
                self.assertEqual(r['model_version'],'2.0')
                self.assertIs(r['physical_model_qualified'],False)

    def test_api_rejects_boolean_shape_factor_instead_of_ignoring_it(self):
        x=self.payload();x['laminar_darcy_poiseuille_number']=True
        with self.assertRaises(DomainError):api_routes.calculate_nozzle_thermal_route(x)

    def test_api_preserves_critical_limit_warning_when_displayed_outlet_equals_limit(self):
        x=self.payload();first=api_routes.calculate_nozzle_thermal_route(x)
        heat=Fraction(first['total_heat_load_w']);inlet=Fraction(x['coolant_inlet_temperature_k'])
        for cp in range(1000,1016):
            actual=inlet+heat/(Fraction(x['coolant_mass_flow_kg_s'])*cp)
            limit=float(actual)
            if actual>Fraction(limit):break
        else:self.fail('No independently established rounded-limit witness')
        x.update(coolant_specific_heat_j_kg_k=float(cp),maximum_coolant_temperature_k=limit)
        r=api_routes.calculate_nozzle_thermal_route(x)['cooling']
        self.assertEqual(r['outlet_temperature_k'],limit)
        self.assertIs(r['temperature_limit_exceeded'],True)
        self.assertTrue(any(w['code']=='coolant_temperature_limit_exceeded' and w['severity']=='critical' for w in r['warnings']))


if __name__=='__main__':unittest.main()
