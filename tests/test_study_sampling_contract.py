"""Adversarial public study-sampling contracts; no hardware validation claims."""

import copy
import json
import pickle
import unittest
from dataclasses import asdict, replace
from decimal import Decimal, localcontext
from fractions import Fraction
from math import isfinite, log, nextafter
from random import Random
from sys import float_info
from unittest.mock import patch

from rocket_propulsion.api.routes import dispatch_calculation
from rocket_propulsion.core.constants import STANDARD_GRAVITY_M_S2
from rocket_propulsion.core.errors import DomainError, InputError
from rocket_propulsion.studies import (
    ParameterRange,
    SensitivityRecord,
    SweepPoint,
    latin_hypercube_samples,
    monte_carlo_samples,
    parameter_sweep,
    propagate_uncertainty,
)


def decimal_moments(values):
    """Independent high-precision oracle, never importing production kernels."""
    with localcontext() as context:
        context.prec = 2200
        rows = [Decimal.from_float(x) for x in values]
        mean = sum(rows) / len(rows)
        deviation = (sum((x-mean)**2 for x in rows)/(len(rows)-1)).sqrt()
        return float(mean), float(deviation)


class StudySamplingAttacks(unittest.TestCase):
    def valid(self, function, *args, **kwargs):
        try:
            return function(*args, **kwargs)
        except (DomainError, ArithmeticError) as error:
            self.fail(f"Representable valid case rejected: {error}")

    def result(self, **kwargs):
        options = {"evaluator": lambda row: row['x'], "sample_count": 7, "seed": 19}
        options.update(kwargs)
        return self.valid(propagate_uncertainty, (ParameterRange('x', 1., 2.),), **options)

    def test_finite_opposite_bounds_have_exact_affine_monte_carlo_samples(self):
        bound = float_info.max
        rng = Random(1)
        expected = [float((1-Fraction(u))*Fraction(-bound)+Fraction(u)*Fraction(bound))
                    for u in (rng.random() for _ in range(8))]
        points = self.valid(monte_carlo_samples, (ParameterRange('x', -bound, bound),),
                            sample_count=8, seed=1)
        self.assertEqual([p['x'] for p in points], expected)
        self.assertTrue(all(isfinite(p['x']) for p in points))

    def test_latin_hypercube_preserves_exact_latent_strata_and_rng_order(self):
        bound = float_info.max
        rng = Random(22)
        q = [(Fraction(i)+Fraction(rng.random()))/7 for i in range(7)]
        expected = [float((1-p)*Fraction(-bound)+p*Fraction(bound)) for p in q]
        rng.shuffle(expected)
        points = self.valid(latin_hypercube_samples, (ParameterRange('x', -bound, bound),),
                            sample_count=7, seed=22)
        self.assertEqual([p['x'] for p in points], expected)
        self.assertEqual(sorted(int((Fraction(p['x'])+Fraction(bound))/(2*Fraction(bound))*7)
                                for p in points), list(range(7)))

    def test_sweep_is_inclusive_and_does_not_overflow_the_span(self):
        bound = float_info.max
        points = self.valid(parameter_sweep, ParameterRange('x', -bound, bound),
                            point_count=3, evaluator=lambda row: row['x'])
        self.assertEqual([p.output for p in points], [-bound, 0., bound])

    def test_smallest_jitter_is_not_erased_before_affine_transport(self):
        tiny = float.fromhex('0x0.0000000000001p-1022')
        # Exact q=tiny/2 would round to zero before transport, but MAX*q is normal.
        with patch('rocket_propulsion.studies.sampling.Random.random', return_value=tiny):
            points = self.valid(latin_hypercube_samples,
                                (ParameterRange('x', 0., float_info.max),), sample_count=2, seed=1)
        expected = [float(Fraction(float_info.max)*(Fraction(i)+Fraction(tiny))/2)
                    for i in range(2)]
        self.assertEqual(sorted(p['x'] for p in points), expected)
        self.assertGreater(expected[0], 0.)

    def test_binary_rounding_at_adjacent_subnormal_bounds_is_allowed(self):
        tiny = float.fromhex('0x0.0000000000001p-1022')
        points = self.valid(monte_carlo_samples, (ParameterRange('x', 0., tiny),), sample_count=4, seed=1)
        self.assertEqual([p['x'] for p in points], [0., tiny, tiny, 0.])

    def test_constant_huge_output_has_finite_mean_zero_dispersion_and_null_pearson(self):
        result = self.result(evaluator=lambda row: 1e308)
        self.assertEqual(result.output_mean, 1e308)
        self.assertEqual(result.output_standard_deviation, 0.)
        self.assertEqual(dict(result.output_percentiles), dict.fromkeys(('p05','p50','p95'), 1e308))
        record, = result.sensitivity
        self.assertIsNone(record.correlation)
        self.assertIsNone(record.absolute_correlation)
        self.assertIsNone(record.rank)
        self.assertEqual(record.undefined_reason, 'constant_input_or_output')

    def test_tiny_and_huge_outputs_match_independent_decimal_moments(self):
        for scale in (2e-200, 2e200):
            with self.subTest(scale=scale):
                result = self.result(evaluator=lambda row, *, scale=scale: scale*row['x'])
                mean, deviation = decimal_moments([p.output for p in result.samples])
                self.assertEqual(result.output_mean, mean)
                self.assertEqual(result.output_standard_deviation, deviation)
                self.assertGreater(deviation, 0.)
                self.assertEqual(result.sensitivity[0].correlation, 1.)

    def test_named_percentiles_use_rational_five_percent_positions(self):
        points = tuple({'x': float(i)} for i in range(7))
        with patch('rocket_propulsion.studies.sampling.latin_hypercube_samples', return_value=points):
            result = self.valid(propagate_uncertainty, (ParameterRange('x',0.,6.),),
                                evaluator=lambda row: row['x'], sample_count=7, seed=1)
        self.assertEqual(dict(result.output_percentiles), {'p05': .3, 'p50': 3., 'p95': 5.7})

    def test_opposite_huge_outputs_have_representable_sample_dispersion(self):
        points=({'x':-1.},{'x':1.})
        with patch('rocket_propulsion.studies.sampling.monte_carlo_samples',return_value=points):
            result=self.valid(propagate_uncertainty,(ParameterRange('x',-1.,1.),),
                              evaluator=lambda row: row['x']*(float_info.max/2),
                              sample_count=2,method='monte-carlo')
        mean,deviation=decimal_moments([point.output for point in result.samples])
        self.assertEqual((result.output_mean,result.output_standard_deviation),(mean,deviation))
        self.assertEqual(mean,0.)

    def test_unrepresentable_required_summary_is_a_named_domain_failure(self):
        points=({'x':0.},{'x':1.})
        tiny=float.fromhex('0x0.0000000000001p-1022')
        with (
            patch('rocket_propulsion.studies.sampling.monte_carlo_samples', return_value=points),
            self.assertRaises(DomainError) as caught,
        ):
            propagate_uncertainty((ParameterRange('x',0.,1.),),
                                  evaluator=lambda row: row['x']*2*tiny,
                                  sample_count=2,method='monte-carlo')
        self.assertEqual(caught.exception.field,'Sample quantile')

    def test_undefined_constant_input_is_last_and_defined_zero_remains_ranked(self):
        points = ({'constant':1.,'zero':-1.,'linear':-1.},
                  {'constant':1.,'zero':1.,'linear':0.},
                  {'constant':1.,'zero':-1.,'linear':1.})
        parameters = tuple(ParameterRange(name,-2.,2.) for name in points[0])
        with patch('rocket_propulsion.studies.sampling.monte_carlo_samples', return_value=points):
            result = self.valid(propagate_uncertainty, parameters, evaluator=lambda row: -row['linear'],
                                sample_count=3, seed=1, method='monte-carlo')
        self.assertEqual([(r.parameter,r.correlation,r.rank) for r in result.sensitivity],
                         [('linear',-1.,1),('zero',0.,2),('constant',None,None)])

    def test_parameter_ranges_reject_malformed_names_and_bounds_at_construction(self):
        for values in ((True,0.,1.), ('',0.,1.), ('\ud800',0.,1.), ('x',True,2.),
                       ('x','0',1.), ('x',None,1.), ('x',0.,float('inf')),
                       ('x',0.,10**400), ('x',1.,1.), ('x',2.,1.)):
            with self.subTest(values=repr(values)), self.assertRaises(DomainError):
                ParameterRange(*values)

    def test_sample_counts_seeds_and_methods_fail_typed_before_evaluation(self):
        calls = []
        for options in ({'sample_count':True}, {'sample_count':3.5}, {'sample_count':'3'},
                        {'sample_count':1}, {'sample_count':100001}, {'seed':True},
                        {'seed':1.5}, {'seed':None}, {'seed':'1'}, {'method':None}, {'method':3}):
            with self.subTest(options=options), self.assertRaises(DomainError):
                propagate_uncertainty((ParameterRange('x',0.,1.),),
                                      evaluator=lambda row: calls.append(row) or 1., **options)
        self.assertEqual(calls, [])

    def test_parameter_collection_is_snapshotted_and_duplicates_rejected(self):
        parameters = [ParameterRange('x',1.,2.)]
        def evaluate(row):
            parameters.clear()
            return row['x']
        result = self.valid(propagate_uncertainty, parameters, evaluator=evaluate, sample_count=4)
        self.assertEqual(result.sensitivity[0].parameter, 'x')
        self.assertEqual(result.parameters, (ParameterRange('x',1.,2.),))
        for bad in (None, 'x', {'x':1}, [object()], [ParameterRange('x',1.,2.)]*2):
            with self.subTest(bad=repr(bad)), self.assertRaises(DomainError):
                monte_carlo_samples(bad, sample_count=3, seed=1)

    def test_evaluator_cannot_relabel_recorded_inputs(self):
        attempts = []
        def evaluator(row):
            try:
                row['x'] = -999.
            except TypeError:
                attempts.append(True)
            return row['x']
        result = self.result(evaluator=evaluator)
        self.assertEqual(len(attempts), 7)
        self.assertTrue(all(1. <= point.inputs['x'] <= 2. and point.output == point.inputs['x']
                            for point in result.samples))

    def test_evaluator_rejects_bool_text_nonfinite_and_unrepresentable_output(self):
        for value in (True, '1.0', None, float('nan'), float('inf'), 10**400):
            with self.subTest(value=repr(value)), self.assertRaises(DomainError):
                self.result_without_valid(evaluator=lambda row, *, value=value: value)
        with self.assertRaises(DomainError):
            propagate_uncertainty((ParameterRange('x',1.,2.),), evaluator=None, sample_count=4)

    def result_without_valid(self, **options):
        return propagate_uncertainty((ParameterRange('x',1.,2.),), sample_count=4, **options)

    def test_unrelated_callback_failure_is_preserved(self):
        class CallbackFailure(RuntimeError):
            pass
        def evaluator(row):
            raise CallbackFailure('application failure')
        with self.assertRaisesRegex(CallbackFailure, 'application failure'):
            self.result_without_valid(evaluator=evaluator)

    def test_records_snapshot_mappings_and_remain_asdict_json_pickle_compatible(self):
        source = {'x': 2.}
        point = SweepPoint(source, 3.)
        source['x'] = 99.
        self.assertEqual(point.inputs['x'], 2.)
        for record in (point, copy.deepcopy(point), pickle.loads(pickle.dumps(point))):
            for operation in (lambda*, record=record: record.inputs.update(x=99.), lambda*, record=record: record.inputs.clear(),
                              lambda*, record=record: record.inputs.pop('x'), lambda*, record=record: record.inputs.setdefault('y', 1.),
                              lambda*, record=record: record.inputs.__ior__({'x':99.}),
                              lambda*, record=record: record.inputs.__delitem__('x')):
                with self.assertRaises(TypeError):
                    operation()
            self.assertEqual(json.loads(json.dumps(asdict(record))), {'inputs':{'x':2.},'output':3.})
            record.inputs.__init__({'x':99.})
            self.assertEqual(record.inputs['x'],2.)
        result = self.result()
        with self.assertRaises(TypeError):
            result.output_percentiles['p05'] = 0.
        self.assertEqual(json.loads(json.dumps(asdict(result)))['sample_count'], 7)

    def test_sweep_validates_fixed_inputs_counts_and_callback_outputs(self):
        for options in ({'point_count':True}, {'point_count':3.5}, {'point_count':10001},
                        {'fixed_inputs':{'z':float('inf')}}, {'fixed_inputs':{'z':True}},
                        {'fixed_inputs':[]}, {'evaluator':lambda row: '3'}):
            args = {"point_count": 3, "evaluator": lambda row: row['x']}
            args.update(options)
            with self.subTest(options=repr(options)), self.assertRaises(DomainError):
                parameter_sweep(ParameterRange('x',0.,1.), **args)

    def test_summary_and_record_forgery_is_rejected_by_owning_result(self):
        result = self.result()
        forged = replace(result.samples[0], output=result.samples[0].output+1.)
        changes = ({'output_mean':nextafter(result.output_mean, float('inf'))},
                   {'output_standard_deviation':0.}, {'output_percentiles':{'p05':0.,'p50':0.,'p95':0.}},
                   {'samples':(forged,)+result.samples[1:]}, {'samples':result.samples[:-1]},
                   {'sample_count':True}, {'sample_count':result.sample_count+1}, {'seed':True}, {'sensitivity':()},
                   {'parameters':()}, {'sampling_numerics':'unknown'}, {'summary_numerics':'unknown'},
                   {'sensitivity':(replace(result.sensitivity[0],correlation=.5,absolute_correlation=.5),)},
                   {'sensitivity':(replace(result.sensitivity[0],rank=2),)},
                   {'samples':(replace(result.samples[0],inputs={'x':3.}),)+result.samples[1:]},
                   {'samples':(replace(result.samples[0],inputs={'other':1.}),)+result.samples[1:]})
        for change in changes:
            with self.subTest(change=repr(change)), self.assertRaises(DomainError):
                replace(result, **change)

    def test_local_records_reject_inconsistent_defined_and_undefined_sensitivity(self):
        for args in (('x',2.,2.,1), ('x',.5,.4,1), ('x',None,None,1),
                     ('x',0.,0.,True), ('x',float('nan'),0.,1)):
            with self.subTest(args=args), self.assertRaises(DomainError):
                SensitivityRecord(*args)
        for inputs, output in (({'x':True},1.), ({'x':1.},True), ({'x':1.},'2'), ([],1.)):
            with self.subTest(inputs=inputs, output=output), self.assertRaises(DomainError):
                SweepPoint(inputs,output)

    def test_full_declared_sample_and_sweep_limits_are_preserved(self):
        parameters = (ParameterRange('x',0.,1.),)
        for sampler in (monte_carlo_samples, latin_hypercube_samples):
            points = self.valid(sampler, parameters, sample_count=100000, seed=1)
            self.assertEqual(len(points),100000)
            self.assertTrue(all(0. <= point['x'] <= 1. for point in points))
        points = self.valid(parameter_sweep, parameters[0], point_count=10000, evaluator=lambda row: row['x'])
        self.assertEqual((len(points),points[0].output,points[-1].output),(10000,0.,1.))


class StudyUncertaintyApiAttacks(unittest.TestCase):
    def payload(self):
        return {"parameters": [{"name": 'specific_impulse_s',"minimum": 430.,"maximum": 460.},
                                {"name": 'mass_ratio',"minimum": 3.,"maximum": 4.}], "sample_count": 4, "seed": 7}

    def call(self, payload):
        return dispatch_calculation('/api/v1/uncertainty', payload)

    def test_api_preserves_exact_large_integer_seed_and_exposes_policy(self):
        seed = 2**60+3
        for raw_seed in (seed,str(seed)):
            payload = self.payload(); payload['seed'] = raw_seed
            result = self.call(payload)
            self.assertEqual(result['seed'], seed)
            self.assertIn('sampling_numerics',result)
            self.assertIn('summary_numerics',result)
            self.assertEqual(len(result['parameters']),2)
            self.assertTrue(any('empirical' in warning.lower() for warning in result['warnings']))
            json.dumps(result,allow_nan=False)

    def test_route_recovers_representable_delta_v_after_overflowing_intermediate(self):
        payload = self.payload()
        payload['parameters'] = [{"name": 'specific_impulse_s',"minimum": 5e307,"maximum": 6e307},
                                 {"name": 'mass_ratio',"minimum": nextafter(1.,2.),"maximum": nextafter(nextafter(1.,2.),2.)}]
        try:
            result = self.call(payload)
        except (DomainError, ArithmeticError) as error:
            self.fail(f'Representable route rejected: {error}')
        for point in result['samples']:
            inputs = point['inputs']
            with localcontext() as context:
                context.prec=2200
                expected=float(Decimal.from_float(STANDARD_GRAVITY_M_S2)*Decimal.from_float(inputs['specific_impulse_s'])
                               * Decimal.from_float(log(inputs['mass_ratio'])))
            self.assertEqual(point['output'],expected)

    def test_route_validates_entire_physical_range_before_any_draw(self):
        for name,lo,hi in (('specific_impulse_s',-1.,460.), ('mass_ratio',0.,4.), ('mass_ratio',.5,4.)):
            payload=self.payload()
            next(p for p in payload['parameters'] if p['name']==name).update(minimum=lo,maximum=hi)
            with self.subTest(name=name,lo=lo), patch('rocket_propulsion.api.routes.propagate_uncertainty') as evaluate:
                evaluate.side_effect=lambda *args,**kwargs: self.fail('Invalid physical range reached evaluation.')
                with self.assertRaises((DomainError,InputError)):
                    self.call(payload)
                evaluate.assert_not_called()

    def test_route_rejects_malformed_objects_duplicate_names_and_control_types(self):
        payloads=[]
        for key,value in (('method',None),('method',True),('sample_count',True),('sample_count','3.2'),('seed',True)):
            payload=self.payload();payload[key]=value;payloads.append(payload)
        payload=self.payload();payload['parameters'].append(dict(payload['parameters'][0]));payloads.append(payload)
        payload=self.payload();payload['parameters'][0]['minimum']=10**400;payloads.append(payload)
        payload=self.payload();payload['parameters'][0]['name']=True;payloads.append(payload)
        payload=self.payload();payload['parameters'].append([]);payloads.append(payload)
        for payload in payloads:
            with self.subTest(payload=repr(payload)), self.assertRaises((DomainError,InputError)):
                self.call(payload)

    def test_form_payload_is_supported_without_changing_its_seed(self):
        payload={"minimum_specific_impulse_s": '430',"maximum_specific_impulse_s": '460',
                     "minimum_mass_ratio": '3',"maximum_mass_ratio": '4',"sample_count": '4',"seed": '42.0',"method": ' Monte-Carlo '}
        result=self.call(payload)
        self.assertEqual((result['seed'],result['sample_count'],result['method']),(42,4,'monte-carlo'))
