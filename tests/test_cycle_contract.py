"""Independent component-cycle balances, routing contracts and hostile records."""
import json
import unittest
from collections.abc import Mapping
from dataclasses import asdict, replace
from decimal import Decimal, localcontext
from math import nextafter
from sys import float_info

from rocket_propulsion.api.routes import dispatch_calculation, calculation_metadata
from rocket_propulsion.core.errors import DomainError, ConvergenceError, InputError
from rocket_propulsion.propulsion import calculate_engine_cycle


def base_inputs(**changes):
    values=dict(cycle='gas-generator',chamber_pressure_pa=7e6,chamber_mass_flow_kg_s=25.,
        propellant_density_kg_m3=1030.,inlet_pressure_pa=2.5e5,injector_pressure_drop_fraction=.2,
        pump_efficiency=.68,turbine_efficiency=.65,turbine_inlet_temperature_k=900.,
        turbine_specific_heat_j_kg_k=2200.,turbine_gamma=1.25,turbine_pressure_ratio=8.)
    return dict(values,**changes)


def reference(inputs,precision):
    with localcontext() as context:
        context.prec=precision
        x={k:Decimal.from_float(v) for k,v in inputs.items() if isinstance(v,float)}
        flow=x['chamber_mass_flow_kg_s'];required=x['chamber_pressure_pa']*(1+x['injector_pressure_drop_fraction'])
        margin=x['inlet_pressure_pa']-required
        rise=max(Decimal(0),-margin)
        specific=rise/(x['propellant_density_kg_m3']*x['pump_efficiency'])
        exponent=(x['turbine_gamma']-1)/x['turbine_gamma']
        work=x['turbine_efficiency']*x['turbine_specific_heat_j_kg_k']*x['turbine_inlet_temperature_k']*(1-(-exponent*x['turbine_pressure_ratio'].ln()).exp())
        if inputs.get('pump_feed_policy','shared-pump')=='shared-pump' and inputs['cycle']=='gas-generator':
            turbine=flow*specific/(work-specific)
            total=flow*work/(work-specific);pumped=total
            fraction=specific/work;delivered=1-fraction
        else:
            turbine=flow*specific/work;pumped=flow
            total=flow+turbine if inputs['cycle']=='gas-generator' else flow
            fraction=turbine/total;delivered=flow/total
        power=pumped*specific
        return {k:float(v) for k,v in dict(required_feed_pressure_pa=required,
            inlet_pressure_margin_pa=margin,pump_pressure_rise_pa=rise,pump_power_w=power,
            turbine_specific_work_j_kg=work,turbine_mass_flow_kg_s=turbine,turbine_flow_fraction=fraction,
            turbine_power_w=power,power_balance_residual_w=Decimal(0),total_propellant_flow_kg_s=total,
            delivered_mass_fraction=delivered,pump_mass_flow_kg_s=pumped).items()}


class CycleAttacks(unittest.TestCase):
    def valid(self,inputs):
        try:return calculate_engine_cycle(**inputs)
        except (DomainError,ConvergenceError,ArithmeticError) as error:self.fail(f'Valid cycle case rejected: {error}')

    def check(self,inputs,precision=1000):
        expected=reference(inputs,precision)
        self.assertEqual(expected,reference(inputs,precision+400))
        actual=self.valid(inputs)
        for name in ('required_feed_pressure_pa','pump_pressure_rise_pa','pump_power_w',
                     'turbine_specific_work_j_kg','turbine_mass_flow_kg_s','total_propellant_flow_kg_s'):
            self.assertEqual(getattr(actual,name),expected[name],name)
        self.assertEqual({name:getattr(actual,name) for name in expected},expected)
        return actual

    def test_shared_pump_includes_gas_generator_flow_in_its_power_duty(self):
        result=self.check(base_inputs())
        self.assertGreater(result.pump_mass_flow_kg_s,result.chamber_mass_flow_kg_s)
        self.assertEqual(result.pump_mass_flow_kg_s,result.total_propellant_flow_kg_s)
        self.assertEqual(result.inputs.pump_feed_policy,'shared-pump')

    def test_external_generator_feed_is_explicit_and_warns_about_excluded_supply(self):
        result=self.check(base_inputs(pump_feed_policy='external-gg-feed'))
        self.assertEqual(result.pump_mass_flow_kg_s,result.chamber_mass_flow_kg_s)
        self.assertGreater(result.total_propellant_flow_kg_s,result.pump_mass_flow_kg_s)
        self.assertIn('external_gg_feed_not_modelled',{warning.code for warning in result.warnings})
        self.assertIn('excludes_external_gg_supply',result.energy_balance_scope)

    def test_required_pressure_rounding_does_not_erase_positive_pump_head(self):
        result=self.check(base_inputs(chamber_pressure_pa=1e8,chamber_mass_flow_kg_s=1.,
            propellant_density_kg_m3=1000.,inlet_pressure_pa=1e8,injector_pressure_drop_fraction=1e-20,pump_efficiency=.5))
        self.assertEqual(result.required_feed_pressure_pa,result.inputs.inlet_pressure_pa)
        self.assertGreater(result.pump_pressure_rise_pa,0.)

    def test_recoverable_pump_product_overflow_is_rescued_before_export(self):
        self.check(base_inputs(chamber_pressure_pa=1e200,chamber_mass_flow_kg_s=1e200,
            propellant_density_kg_m3=1e200,inlet_pressure_pa=1.4e200,injector_pressure_drop_fraction=.5,
            pump_efficiency=1.,turbine_efficiency=1.,turbine_inlet_temperature_k=1e200,
            turbine_specific_heat_j_kg_k=1e100))

    def test_recoverable_pump_denominator_underflow_is_not_division_by_zero(self):
        self.check(base_inputs(chamber_pressure_pa=1e-200,chamber_mass_flow_kg_s=1e-250,
            propellant_density_kg_m3=1e-200,inlet_pressure_pa=0.,injector_pressure_drop_fraction=.5,
            pump_efficiency=1e-200,turbine_efficiency=1.,turbine_inlet_temperature_k=1e5,
            turbine_specific_heat_j_kg_k=1e200))

    def test_recoverable_turbine_product_underflow_remains_supported_with_external_supply(self):
        self.check(base_inputs(pump_feed_policy='external-gg-feed',chamber_pressure_pa=1e-97,
            chamber_mass_flow_kg_s=1.,propellant_density_kg_m3=1000.,inlet_pressure_pa=0.,pump_efficiency=.5,
            turbine_efficiency=1e-300,turbine_inlet_temperature_k=1e300,turbine_specific_heat_j_kg_k=1e-100))

    def test_near_unity_gamma_and_pressure_ratio_do_not_erase_turbine_work(self):
        self.check(base_inputs(chamber_pressure_pa=1.,chamber_mass_flow_kg_s=1.,inlet_pressure_pa=1.19,
            propellant_density_kg_m3=1000.,pump_efficiency=.5,turbine_efficiency=1.,
            turbine_inlet_temperature_k=1e10,turbine_specific_heat_j_kg_k=1e20,
            turbine_gamma=nextafter(1.,2.),turbine_pressure_ratio=nextafter(1.,2.)))

    def test_shared_pump_feedback_can_be_infeasible_when_chamber_only_duty_looks_closed(self):
        # gamma=2 and pressure ratio=4 give exact turbine work=.5.
        inputs=base_inputs(chamber_pressure_pa=1.,chamber_mass_flow_kg_s=1.,inlet_pressure_pa=0.,
            injector_pressure_drop_fraction=.5,propellant_density_kg_m3=1.,pump_efficiency=1.,
            turbine_efficiency=1.,turbine_inlet_temperature_k=1.,turbine_specific_heat_j_kg_k=1.,
            turbine_gamma=2.,turbine_pressure_ratio=4.)
        external=self.valid(dict(inputs,pump_feed_policy='external-gg-feed'))
        self.assertEqual(external.turbine_mass_flow_kg_s,3.)
        with self.assertRaises(DomainError) as caught:calculate_engine_cycle(**inputs)
        self.assertEqual(caught.exception.code,'infeasible_cycle_power_balance')

    def exact_boundary(self,cycle):
        return base_inputs(cycle=cycle,chamber_pressure_pa=1.,chamber_mass_flow_kg_s=1.,
            inlet_pressure_pa=1.,injector_pressure_drop_fraction=.5,propellant_density_kg_m3=1.,pump_efficiency=1.,
            turbine_efficiency=1.,turbine_inlet_temperature_k=1.,turbine_specific_heat_j_kg_k=1.,
            turbine_gamma=2.,turbine_pressure_ratio=4.)

    def test_exact_staged_full_flow_equality_is_valid_without_a_tolerance(self):
        result=self.valid(self.exact_boundary('staged-combustion'))
        self.assertEqual(result.turbine_flow_fraction,1.)
        self.assertEqual(result.turbine_mass_flow_kg_s,1.)
        self.assertEqual(result.pump_power_w,.5)

    def test_rational_turbine_roots_cover_noninteger_ratio_and_higher_degree(self):
        for changes in ({'turbine_gamma':1.5,'turbine_pressure_ratio':8.},
                        {'turbine_pressure_ratio':2.25,'propellant_density_kg_m3':1.5}):
            with self.subTest(changes=changes):
                result=self.valid(dict(self.exact_boundary('staged-combustion'),**changes))
                self.assertEqual(result.turbine_flow_fraction,1.)
                self.assertEqual(result.turbine_mass_flow_kg_s,1.)

    def test_shared_generator_remains_supported_just_below_capacity(self):
        inputs=self.exact_boundary('gas-generator');inputs['inlet_pressure_pa']=nextafter(1.,2.)
        result=self.check(inputs)
        self.assertGreater(result.turbine_flow_fraction,.99)
        self.assertGreater(result.delivered_mass_fraction,0.)
        self.assertIn('high_gas_generator_flow_fraction',{w.code for w in result.warnings})

    def test_unrounded_turbine_work_decides_a_close_generator_feedback_boundary(self):
        self.check(base_inputs(chamber_pressure_pa=1.,chamber_mass_flow_kg_s=1.,inlet_pressure_pa=1.,
            injector_pressure_drop_fraction=.12944943670387585,propellant_density_kg_m3=1.,pump_efficiency=1.,
            turbine_efficiency=1.,turbine_inlet_temperature_k=1.,turbine_specific_heat_j_kg_k=1.,
            turbine_gamma=1.25,turbine_pressure_ratio=2.))

    def reject_power(self,inputs):
        try:calculate_engine_cycle(**inputs)
        except DomainError as error:self.assertEqual(error.code,'infeasible_cycle_power_balance')
        except (ConvergenceError,ArithmeticError) as error:self.fail(f'Known infeasible capacity was not rejected typed: {error}')
        else:self.fail('Known infeasible capacity was accepted.')

    def test_exact_shared_generator_equality_is_infeasible(self):
        self.reject_power(self.exact_boundary('gas-generator'))

    def test_one_ulp_beyond_staged_turbine_capacity_is_rejected(self):
        inputs=self.exact_boundary('staged-combustion');inputs['pump_efficiency']=nextafter(1.,0.)
        self.reject_power(inputs)

    def test_staged_partial_flow_matches_an_independent_component_balance(self):
        result=self.check(base_inputs(cycle='staged-combustion'))
        self.assertEqual(result.total_propellant_flow_kg_s,result.chamber_mass_flow_kg_s)
        self.assertEqual(result.delivered_mass_fraction,1.)
        self.assertLess(result.turbine_flow_fraction,1.)

    def test_zero_pump_head_has_zero_turbine_demand_and_positive_available_work(self):
        result=self.check(base_inputs(inlet_pressure_pa=1e8))
        self.assertEqual(result.pump_power_w,0.)
        self.assertEqual(result.turbine_mass_flow_kg_s,0.)
        self.assertGreater(result.turbine_specific_work_j_kg,0.)

    def test_pressure_fed_sizing_exposes_insufficient_supply_and_requires_no_turbine(self):
        inputs=dict(cycle='pressure-fed',chamber_pressure_pa=1e8,chamber_mass_flow_kg_s=1.,
                    inlet_pressure_pa=1e8,injector_pressure_drop_fraction=1e-20)
        result=self.valid(inputs)
        self.assertEqual(result.required_feed_pressure_pa,1e8)
        self.assertLess(result.inlet_pressure_margin_pa,0.)
        self.assertIs(result.pressure_fed_supply_sufficient,False)
        self.assertIn('insufficient_pressure_fed_supply',{w.code for w in result.warnings})
        self.assertEqual(result.pump_power_w,0.)
        self.assertIsNone(result.inputs.turbine_gamma)

    def test_pressure_fed_exact_margin_zero_is_sufficient_for_declared_losses(self):
        result=self.valid(dict(cycle='pressure-fed',chamber_pressure_pa=2.,chamber_mass_flow_kg_s=1.,
                              inlet_pressure_pa=3.,injector_pressure_drop_fraction=.5))
        self.assertIs(result.pressure_fed_supply_sufficient,True)
        self.assertEqual(result.inlet_pressure_margin_pa,0.)

    def test_unrepresentable_optional_margin_does_not_restrict_pressure_sizing(self):
        tiny=float.fromhex('0x0.0000000000001p-1022')
        result=self.valid(dict(cycle='pressure-fed',chamber_pressure_pa=tiny,chamber_mass_flow_kg_s=1.,
                              inlet_pressure_pa=tiny,injector_pressure_drop_fraction=tiny))
        self.assertIsNone(result.inlet_pressure_margin_pa)
        self.assertIs(result.inlet_pressure_margin_representable,False)
        self.assertIs(result.pressure_fed_supply_sufficient,False)

    def test_bad_cycle_and_routing_names_fail_typed(self):
        for changes in ({'cycle':None},{'cycle':True},{'cycle':'unknown'},{'pump_feed_policy':True},
                        {'pump_feed_policy':'unknown'},{'cycle':'staged-combustion','pump_feed_policy':'external-gg-feed'}):
            with self.subTest(changes=changes),self.assertRaises(DomainError):calculate_engine_cycle(**base_inputs(**changes))

    def test_all_supplied_quantities_reject_bool_text_nonfinite_and_overflow(self):
        inputs=base_inputs()
        for name in inputs:
            if name=='cycle':continue
            for bad in (True,'1',[],float('nan'),float('inf'),10**400):
                with self.subTest(name=name,bad=repr(bad)),self.assertRaises(DomainError):
                    calculate_engine_cycle(**dict(inputs,**{name:bad}))

    def test_pumped_modes_require_their_active_inputs(self):
        for name in ('propellant_density_kg_m3','pump_efficiency','turbine_efficiency',
                     'turbine_inlet_temperature_k','turbine_specific_heat_j_kg_k','turbine_gamma','turbine_pressure_ratio'):
            with self.subTest(name=name),self.assertRaises(DomainError):
                calculate_engine_cycle(**base_inputs(**{name:None}))

    def test_unrepresentable_required_feed_pressure_fails_named(self):
        with self.assertRaises(DomainError) as caught:
            calculate_engine_cycle(**base_inputs(chamber_pressure_pa=float_info.max,injector_pressure_drop_fraction=1.))
        self.assertEqual(caught.exception.field,'Required feed pressure')

    def test_record_rejects_every_forged_balance_field_and_identity(self):
        result=self.valid(base_inputs())
        for name in reference(base_inputs(),1000):
            with self.subTest(name=name),self.assertRaises(DomainError):
                replace(result,**{name:nextafter(getattr(result,name),float('inf'))})
        for changes in ({'numerical_policy':'unknown'},{'inputs':None},{'warnings':()},
                        {'assumptions':()},{'energy_balance_scope':'unknown'},{'pressure_fed_supply_sufficient':True},
                        {'inlet_pressure_margin_representable':False},{'inlet_pressure_margin_representable':1}):
            with self.subTest(changes=changes),self.assertRaises(DomainError):replace(result,**changes)

    def test_record_rejects_changed_original_inputs_and_mutable_aliases(self):
        result=self.valid(base_inputs());raw=asdict(result.inputs)
        owned=replace(result,inputs=raw);raw['chamber_pressure_pa']=1.
        self.assertEqual(owned,result)
        with self.assertRaises(DomainError):
            replace(result,inputs=replace(result.inputs,inlet_pressure_pa=0.))
        restored=type(result)(**json.loads(json.dumps(asdict(result),allow_nan=False)))
        self.assertEqual(restored,result)

    def test_unclosed_hot_gas_and_hardware_limits_are_explicit(self):
        result=self.valid(base_inputs())
        self.assertIn('hot_gas_pressure_not_closed',{w.code for w in result.warnings})
        self.assertTrue(any('cavitation' in assumption for assumption in result.assumptions))
        self.assertIn('ideal',result.energy_balance_scope)

    def test_json_warning_records_are_owned_and_reject_forged_schema(self):
        result=self.valid(base_inputs());document=asdict(result)
        restored=type(result)(**document)
        document['warnings'][0]['message']='changed'
        self.assertEqual(restored,result)
        for warning in ({'code':'hot_gas_pressure_not_closed','message':'forged'},
                        {'code':'x','message':'x','unexpected':True},
                        {'code':True,'message':[]},{},dict(asdict(result.warnings[0]),unexpected=True)):
            with self.subTest(warning=warning),self.assertRaises(DomainError):
                replace(result,warnings=(warning,))

    def test_mapping_protocol_errors_are_typed_without_hiding_callback_failures(self):
        class BrokenMapping(Mapping):
            def __init__(self,error):self.error=error
            def __len__(self):return 1
            def __iter__(self):return iter(('bad',))
            def __getitem__(self,key):raise self.error('broken mapping')
        result=self.valid(base_inputs())
        for name in ('inputs','warnings'):
            for error,expected in ((KeyError,DomainError),(RuntimeError,RuntimeError)):
                value=BrokenMapping(error)
                if name=='warnings':value=(value,)
                with self.subTest(name=name,error=error),self.assertRaises(expected):
                    replace(result,**{name:value})


class CycleApiAttacks(unittest.TestCase):
    def test_api_pressure_fed_sizing_works_without_inactive_pump_fields(self):
        result=dispatch_calculation('/api/v1/engine-cycle',dict(cycle='pressure-fed',chamber_pressure_pa='2',
            chamber_mass_flow_kg_s='1',inlet_pressure_pa='1',injector_pressure_drop_fraction='.5'))
        self.assertIs(result['pressure_fed_supply_sufficient'],False)
        self.assertEqual(result['inlet_pressure_margin_pa'],-2.)

    def test_api_retains_explicit_generator_routing_and_current_model_version(self):
        result=dispatch_calculation('/api/v1/engine-cycle',base_inputs(pump_feed_policy='external-gg-feed'))
        self.assertEqual(result['inputs']['pump_feed_policy'],'external-gg-feed')
        self.assertEqual(calculation_metadata('/api/v1/engine-cycle')['version'],'2.0')
        json.dumps(result,allow_nan=False)

    def test_api_rejects_wrong_cycle_control_types(self):
        for changes in ({'cycle':None},{'cycle':True},{'pump_feed_policy':None},{'pump_feed_policy':True}):
            with self.subTest(changes=changes),self.assertRaises((DomainError,InputError)):
                dispatch_calculation('/api/v1/engine-cycle',base_inputs(**changes))
