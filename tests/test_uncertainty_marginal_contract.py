"""Bounded inverse-CDF references and adversarial sampler/correlation contracts."""

import json
import sys
import unittest
from dataclasses import replace
from decimal import Decimal, localcontext
from fractions import Fraction
from math import isfinite, nextafter
from pathlib import Path
from unittest.mock import patch

from rocket_propulsion.core.errors import ConvergenceError, DomainError
from rocket_propulsion.propulsion.burns import (
    CorrelationModel,
    DistributionKind,
    UncertainParameter,
    UncertaintyClass,
    sample_uncertain_parameters,
)
from rocket_propulsion.propulsion.burns import uncertainty as model

MAX=sys.float_info.max
TINY=float.fromhex('0x0.0000000000001p-1022')


def parameter(kind=DistributionKind.UNIFORM, lower=-1., upper=1., nominal=0., sigma=None, name='x'):
    return UncertainParameter(name, nominal, lower, upper, kind, UncertaintyClass.EPISTEMIC,
                              'synthetic marginal oracle', 'numerical contract only', sigma)


class MarginalNumericalAttacks(unittest.TestCase):
    def evaluate(self, item, probability):
        try:
            result=item.value_at_probability(probability)
        except (DomainError,ArithmeticError) as error:
            self.fail(f'Representable bounded quantile was rejected: {error}')
        self.assertTrue(isfinite(result))
        self.assertGreaterEqual(result,item.lower)
        self.assertLessEqual(result,item.upper)
        return result

    def test_uniform_finite_opposite_extremes_recover_exact_affine_values(self):
        item=parameter(lower=-MAX,upper=MAX)
        for p in (.25,.5,.75):
            expected=float((1-Fraction(p))*Fraction(-MAX)+Fraction(p)*Fraction(MAX))
            self.assertEqual(self.evaluate(item,p),expected)

    def test_uniform_subnormal_grid_distinguishes_exact_zero_from_erased_positive(self):
        item=parameter(lower=0.,upper=TINY,nominal=0.)
        with self.assertRaises(DomainError) as caught:item.value_at_probability(.25)
        self.assertIs(type(caught.exception),DomainError)
        self.assertEqual(caught.exception.field,'Marginal value')
        self.assertEqual(self.evaluate(item,.75),TINY)
        self.assertEqual(self.evaluate(parameter(lower=-TINY,upper=TINY),.5),0.)

    def test_triangular_wide_and_tiny_spans_match_decimal_root_ledger(self):
        for a,b,m in ((0.,MAX,MAX/2),(-MAX,MAX,0.),(0.,1e-300,5e-301),
                      (-1e-300,1e-300,0.),(0.,1.,0.),(0.,1.,1.)):
            item=parameter(DistributionKind.TRIANGULAR,a,b,m)
            for p in (.125,.25,.5,.75,.875):
                with localcontext() as context:
                    context.prec=2200
                    lo,hi,mode,u=map(Decimal.from_float,(a,b,m,p))
                    width=hi-lo
                    expected=lo+(u*width*(mode-lo)).sqrt() if u<=(mode-lo)/width else hi-((1-u)*width*(hi-mode)).sqrt()
                with self.subTest(a=a,b=b,p=p):
                    self.assertEqual(self.evaluate(item,p),float(expected))

    def test_triangular_subnormal_midpoint_ties_choose_the_even_destination(self):
        for upper,expected in ((3*TINY,2*TINY),(5*TINY,2*TINY),(9*TINY,4*TINY),
                               (21*TINY,10*TINY),(23*TINY,12*TINY)):
            item=parameter(DistributionKind.TRIANGULAR,0.,upper,upper)
            self.assertEqual(self.evaluate(item,.25),expected)
        with self.assertRaises(DomainError) as caught:
            parameter(DistributionKind.TRIANGULAR,0.,TINY,TINY).value_at_probability(.25)
        self.assertIs(type(caught.exception),DomainError)
        self.assertEqual(caught.exception.field,'Marginal value')

    def test_log_and_truncated_normal_match_independent_mpmath_reference_export(self):
        fixture=Path(__file__).with_name('fixtures')/'marginal_quantile_reference.json'
        rows=json.loads(fixture.read_text())['cases']
        self.assertGreaterEqual(len(rows),20)
        for row in rows:
            if not row['case'].startswith(('normal_bounded_','normal_shifted_')):
                self.assert_reference(row)

    def assert_reference(self,row):
        item=parameter(DistributionKind(row['distribution']),
            float.fromhex(row['lower_hex']),float.fromhex(row['upper_hex']),float.fromhex(row['nominal_hex']),
            None if row['sigma_hex'] is None else float.fromhex(row['sigma_hex']))
        p=float.fromhex(row['probability_hex'])
        with self.subTest(case=row['case']):
            if row['positive_value_rounds_to_zero']:
                with self.assertRaises(DomainError) as caught:item.value_at_probability(p)
                self.assertIs(type(caught.exception),DomainError)
                self.assertEqual(caught.exception.field,'Marginal value')
            else:self.assertEqual(self.evaluate(item,p).hex(),row['expected_float_hex'])

    def test_bounded_normal_quantiles_match_independent_mpmath_reference(self):
        fixture=Path(__file__).with_name('fixtures')/'marginal_quantile_reference.json'
        rows=json.loads(fixture.read_text())['cases']
        selected=[r for r in rows if r['case'].startswith(('normal_bounded_','normal_shifted_'))]
        self.assertEqual(len(selected),10)
        for row in selected:self.assert_reference(row)

    def test_external_decimal_context_cannot_change_quantiles(self):
        items=(parameter(DistributionKind.LOG_UNIFORM,1e-300,1e300,1.),
               parameter(DistributionKind.TRUNCATED_NORMAL,-1.,1.,0.,1.),
               parameter(DistributionKind.TRIANGULAR,-MAX,MAX,0.))
        for item in items:
            expected=self.evaluate(item,.125)
            for precision in (2,7):
                with localcontext() as context:
                    context.prec=precision
                    self.assertEqual(self.evaluate(item,.125),expected)

    def test_inaccurate_normal_hint_still_requires_a_certified_output_cell(self):
        from rocket_propulsion.propulsion.burns import marginal_numerics as numerics
        fixture=json.loads((Path(__file__).with_name('fixtures')/'marginal_quantile_reference.json').read_text())
        selected=[r for r in fixture['cases'] if r['case'] in ('normal_full_'+(.25).hex(),'normal_full_'+(.75).hex())]
        self.assertEqual(len(selected),2)
        for row in selected:
            p=float.fromhex(row['probability_hex']);expected=float.fromhex(row['expected_float_hex'])
            for hint in (nextafter(expected,0.),64. if p<.5 else -64.):
                class HintedNormal:
                    def inv_cdf(self,probability, *, hint=hint):return hint
                with patch.object(numerics,'NormalDist',HintedNormal):
                    self.assertEqual(numerics._normal_quantile(p),expected)

    def test_log_uniform_geometric_median_has_no_log_roundtrip_drift(self):
        for a,b,expected in ((1.,2.,None),(1.,4.,2.),(1e-300,1e300,1.),(TINY,4*TINY,2*TINY)):
            item=parameter(DistributionKind.LOG_UNIFORM,a,b,a)
            # The decimal 1e-300/1e300 inputs are not exact reciprocals.
            with localcontext() as context:
                context.prec=2200
                expected=float((Decimal.from_float(a)*Decimal.from_float(b)).sqrt())
            self.assertEqual(self.evaluate(item,.5),expected)

    def test_probability_and_parameter_numeric_controls_are_typed(self):
        item=parameter()
        for p in (None,True,'0.5',0.,1.,float('nan'),float('inf')):
            with self.subTest(p=repr(p)),self.assertRaises(DomainError):item.value_at_probability(p)
        for name in ('nominal','lower','upper','standard_deviation'):
            values={'distribution':DistributionKind.TRUNCATED_NORMAL,'standard_deviation':.5,name:True}
            with self.subTest(name=name),self.assertRaises(DomainError):replace(item,**values)

    def test_distribution_and_uncertainty_classes_cannot_fall_through_to_assertions(self):
        item=parameter()
        for changes in ({'distribution':'uniform'},{'distribution':None},{'uncertainty_class':'epistemic'}):
            with self.subTest(changes=changes),self.assertRaises(DomainError):replace(item,**changes)

    def test_owned_text_is_literal_and_not_a_subclass_callback(self):
        class Label(str):pass
        item=replace(parameter(),name=Label('x'),source=Label('source'))
        self.assertIs(type(item.name),str)
        self.assertIs(type(item.source),str)


class CorrelationAndSamplerAttacks(unittest.TestCase):
    def test_true_near_singular_positive_definite_matrix_has_no_absolute_pivot_floor(self):
        rho=nextafter(1.,0.)
        try:correlation=CorrelationModel(('a','b'),((1.,rho),(rho,1.)),'synthetic','exact SPD input')
        except DomainError as error:self.fail(f'A true positive definite matrix was rejected: {error}')
        lower=model._cholesky(correlation.matrix)
        with localcontext() as context:
            context.prec=100
            expected=float((1-Decimal.from_float(rho)**2).sqrt())
        self.assertEqual(lower[1][1],expected)
        self.assertGreater(lower[1][1],0.)

    def test_unit_diagonal_and_symmetry_are_not_silently_ignored_within_tolerance(self):
        for matrix in (((1.,.2),(.2+1e-13,1.)),((1.+1e-13,.2),(.2,1.)),
                       ((1.-1e-13,.2),(.2,1.)),((True,.2),(.2,1.))):
            with self.subTest(matrix=matrix),self.assertRaises(DomainError):
                CorrelationModel(('a','b'),matrix,'synthetic','declared correlation contract')

    def test_matrix_names_and_rows_are_copied_before_caller_mutation(self):
        names=['a','b'];matrix=[[1.,.2],[.2,1.]]
        value=CorrelationModel(names,matrix,'synthetic','snapshot contract')
        names[0]='changed';matrix[0][1]=.8
        self.assertEqual(value.parameter_names,('a','b'))
        self.assertEqual(value.matrix,((1.,.2),(.2,1.)))

    def test_sample_controls_refuse_boolean_fractional_and_malformed_arguments(self):
        items=(parameter(),)
        for kwargs in ({'sample_count':4.5,'seed':1},{'sample_count':True,'seed':1},
                       {'sample_count':4,'seed':False},{'sample_count':4,'seed':'1'},
                       {'sample_count':4,'seed':1,'method':None},
                       {'sample_count':4,'seed':1,'correlation':{} }):
            with self.subTest(kwargs=kwargs),self.assertRaises(DomainError):sample_uncertain_parameters(items,**kwargs)
        for items in ((None,),{'x':parameter()},None):
            with self.subTest(items=items),self.assertRaises(DomainError):sample_uncertain_parameters(items,sample_count=4,seed=1)

    def test_zero_random_endpoint_is_resampled_without_a_probability_floor(self):
        class RandomSource:
            def __init__(self):self.values=iter((0.,.125,.25,.5,.75))
            def random(self):return next(self.values)
            def shuffle(self,items):pass
        with patch.object(model,'Random',lambda seed:RandomSource()):
            try:draws=sample_uncertain_parameters((parameter(lower=0.,upper=1.,nominal=.5),),sample_count=4,seed=1,method='monte-carlo')
            except DomainError as error:self.fail(f'Zero RNG endpoint must be resampled: {error}')
        self.assertEqual([row.as_dict()['x'] for row in draws],[.125,.25,.5,.75])

    def test_failed_random_source_terminates_with_a_distinct_convergence_error(self):
        class RandomSource:
            def __init__(self):self.calls=0
            def random(self):self.calls+=1;return 0.
        source=RandomSource()
        with (
            patch.object(model, 'Random', lambda seed: source),
            self.assertRaises(ConvergenceError) as caught,
        ):
            sample_uncertain_parameters((parameter(),),sample_count=4,seed=1)
        self.assertEqual(caught.exception.field,'Unit random variate')
        self.assertGreater(source.calls,0)
        self.assertLessEqual(source.calls,1024)

    def test_smallest_positive_jitter_stays_in_the_first_open_stratum(self):
        class RandomSource:
            def random(self):return TINY
            def shuffle(self,items):pass
        for count in (4,7,10_000):
            with patch.object(model,'Random',lambda seed:RandomSource()):
                try:rows=model._unit_samples(1,sample_count=count,seed=1,method='latin-hypercube')
                except DomainError as error:self.fail(f'First stratum has a representable positive point: {error}')
            self.assertEqual(rows[0][0],TINY)
            for i,(p,) in enumerate(rows):
                self.assertTrue(0<p<1)
                self.assertEqual(int(Fraction(p)*count),i)

    def test_last_latin_hypercube_cells_cannot_round_into_the_next_stratum(self):
        class RandomSource:
            def random(self):return nextafter(1.,0.)
            def shuffle(self,items):pass
        with patch.object(model,'Random',lambda seed:RandomSource()):
            try:draws=sample_uncertain_parameters((parameter(lower=0.,upper=1.,nominal=.5),),sample_count=4,seed=1)
            except DomainError as error:self.fail(f'Open exact LHS cells need directed quantization: {error}')
        for i,row in enumerate(draws):
            p=row.as_dict()['x']
            self.assertTrue(0<p<1)
            self.assertEqual(int(Fraction(p)*4),i)

    def test_copula_does_not_clip_independent_tail_probabilities(self):
        items=(parameter(DistributionKind.TRUNCATED_NORMAL,-MAX,MAX,0.,1.,name='a'),parameter(name='b'))
        correlation=CorrelationModel(('a','b'),((1.,0.),(0.,1.)),'synthetic','identity copula')
        with patch.object(model,'_unit_samples',return_value=[[TINY,.5] for _ in range(4)]):
            try:draws=sample_uncertain_parameters(items,sample_count=4,seed=1,correlation=correlation)
            except DomainError as error:self.fail(f'A finite Gaussian tail quantile must survive: {error}')
        fixture=json.loads((Path(__file__).with_name('fixtures')/'marginal_quantile_reference.json').read_text())
        reference=next(r for r in fixture['cases'] if r['case']=='normal_full_'+TINY.hex())
        self.assertEqual(draws[0].as_dict()['a'].hex(),reference['expected_float_hex'])

    def test_independent_copula_inverse_uses_reference_rounding_at_both_tails(self):
        fixture=json.loads((Path(__file__).with_name('fixtures')/'marginal_quantile_reference.json').read_text())
        rows=[r for r in fixture['cases'] if r['case'].startswith('normal_full_')]
        self.assertEqual(len(rows),4)
        items=(parameter(DistributionKind.TRUNCATED_NORMAL,-MAX,MAX,0.,1.,name='a'),parameter(name='b'))
        correlation=CorrelationModel(('a','b'),((1.,0.),(0.,1.)),'synthetic','independent inverse reference')
        units=[[float.fromhex(r['probability_hex']),.5] for r in rows]
        with patch.object(model,'_unit_samples',return_value=units):
            try:draws=sample_uncertain_parameters(items,sample_count=4,seed=1,correlation=correlation)
            except DomainError as error:self.fail(f'Independent inverse reference must be representable: {error}')
        self.assertEqual([d.as_dict()['a'].hex() for d in draws],[r['expected_float_hex'] for r in rows])

    def test_correlated_gaussian_tail_is_not_exported_as_an_endpoint_probability(self):
        items=tuple(parameter(DistributionKind.TRUNCATED_NORMAL,-MAX,MAX,0.,1.,name=x) for x in ('a','b'))
        rho=.9
        correlation=CorrelationModel(('a','b'),((1.,rho),(rho,1.)),'synthetic','latent tail contract')
        unit=nextafter(1.,0.)
        with patch.object(model,'_unit_samples',return_value=[[unit,unit] for _ in range(4)]):
            try:draws=sample_uncertain_parameters(items,sample_count=4,seed=1,correlation=correlation)
            except DomainError as error:self.fail(f'Finite correlated latent quantile must not become CDF=1: {error}')
        fixture=json.loads((Path(__file__).with_name('fixtures')/'marginal_quantile_reference.json').read_text())
        reference=next(r for r in fixture['cases'] if r['case']=='normal_full_'+unit.hex())
        z=float.fromhex(reference['expected_float_hex'])
        with localcontext() as context:
            context.prec=100
            factor=float((1-Decimal.from_float(rho)**2).sqrt())
            expected=float(Decimal.from_float(rho)*Decimal.from_float(z)+Decimal.from_float(factor)*Decimal.from_float(z))
        self.assertEqual(draws[0].as_dict()['b'],expected)

    def test_draw_records_reject_duplicate_names_and_copy_pairs(self):
        pairs=[['x',1.0]]
        draw=model.UncertaintyDraw(0,pairs);pairs[0][1]=2.
        self.assertEqual(draw.values,(('x',1.),))
        with self.assertRaises(DomainError):model.UncertaintyDraw(True,(('x',1.),))
        with self.assertRaises(DomainError):model.UncertaintyDraw(0,(('x',1.),('x',2.)))

    def test_custom_marginal_outputs_are_real_finite_and_within_declared_bounds(self):
        class Custom(UncertainParameter):
            output=.25
            def value_at_probability(self,p):return self.output
        custom=Custom('x',0.,-1.,1.,DistributionKind.UNIFORM,UncertaintyClass.EPISTEMIC,
                      'synthetic callback','exercise checked extension boundary')
        valid=sample_uncertain_parameters((custom,),sample_count=4,seed=1)
        self.assertEqual([r.as_dict()['x'] for r in valid],[.25]*4)
        for output in (True,'0.25',float('nan'),float('inf'),2.):
            Custom.output=output
            with self.subTest(output=repr(output)),self.assertRaises(DomainError):
                sample_uncertain_parameters((custom,),sample_count=4,seed=1)

    def test_custom_float_callback_cannot_receive_a_silently_clipped_copula_endpoint(self):
        class Custom(UncertainParameter):
            calls=0
            def value_at_probability(self,p):
                type(self).calls+=1
                return .25
        custom=Custom('b',0.,-1.,1.,DistributionKind.UNIFORM,UncertaintyClass.EPISTEMIC,
                      'synthetic callback','finite open probability contract')
        self.assertEqual(sample_uncertain_parameters((custom,),sample_count=4,seed=1)[0].as_dict()['b'],.25)
        Custom.calls=0
        items=(parameter(name='a'),custom)
        correlation=CorrelationModel(('a','b'),((1.,.9),(.9,1.)),'synthetic','latent tail endpoint')
        unit=nextafter(1.,0.)
        with (
            patch.object(model, '_unit_samples', return_value=[[unit, unit] for _ in range(4)]),
            self.assertRaises(DomainError) as caught,
        ):
            sample_uncertain_parameters(items,sample_count=4,seed=1,correlation=correlation)
        self.assertEqual(caught.exception.field,'Sampling probability')
        self.assertEqual(Custom.calls,0)


class SamplingIdentityAttacks(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        from rocket_propulsion.propulsion.burns import build_hierarchical_propulsion_ensemble
        from tests.test_burn_hierarchical_ensemble import HierarchicalPropulsionEnsembleTests
        fixture=HierarchicalPropulsionEnsembleTests();fixture.setUp()
        cls.discrete,cls.definitions=fixture.discrete,fixture.definitions
        cls.ensemble=build_hierarchical_propulsion_ensemble(cls.discrete,cls.definitions,seed=71)

    def test_missing_or_rehashed_unsupported_sampling_policy_cannot_claim_current_identity(self):
        from rocket_propulsion.core.errors import InputError
        from rocket_propulsion.propulsion.burns import (
            hierarchical_ensemble_from_dict,
            hierarchical_ensemble_to_dict,
        )
        from tests.test_ensemble_document_contract import rehash
        raw=hierarchical_ensemble_to_dict(self.ensemble)
        self.assertEqual(hierarchical_ensemble_from_dict(raw),self.ensemble)
        self.assertEqual(raw.get('sampling_numerics'),'open_binary_strata_latent_gaussian_and_certified_bounded_quantiles_v1')
        for value in (None,'legacy_clipped_copula'):
            changed=dict(raw)
            if value is None:changed.pop('sampling_numerics',None)
            else:changed['sampling_numerics']=value
            with self.subTest(value=value):
                with self.assertRaises(InputError) as caught:
                    hierarchical_ensemble_from_dict(rehash(changed))
                self.assertEqual(caught.exception.field,'sampling_numerics')

    def test_sampling_policy_changes_conditional_artifact_identity(self):
        from rocket_propulsion.propulsion.burns import build_hierarchical_propulsion_ensemble
        from rocket_propulsion.propulsion.burns import ensemble as hierarchy
        with patch.object(hierarchy,'_SAMPLING_POLICY','future_test_policy'):
            other=build_hierarchical_propulsion_ensemble(self.discrete,self.definitions,seed=71)
        self.assertEqual([r.inputs for r in other.realizations],[r.inputs for r in self.ensemble.realizations])
        continuous=[i for i,r in enumerate(self.ensemble.realizations) if r.inputs]
        self.assertTrue(continuous)
        for i in continuous:
            self.assertNotEqual(other.realizations[i].artifact.artifact_hash,self.ensemble.realizations[i].artifact.artifact_hash)
