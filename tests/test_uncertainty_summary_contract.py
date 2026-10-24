"""Independent sample-summary ledgers and adversarial retained-result contracts.

Synthetic mocked solver outputs test assembly only, never physical validation.
"""
import sys
import unittest
from dataclasses import replace
from decimal import Decimal, localcontext
from fractions import Fraction
from math import nextafter
from types import SimpleNamespace
from unittest.mock import patch

from rocket_propulsion.core.errors import DomainError
from rocket_propulsion.propulsion.burns import uncertainty as model

MAX=sys.float_info.max
TINY=float.fromhex('0x0.0000000000001p-1022')


def independent_moments(values):
    with localcontext() as context:
        context.prec=2200
        rows=list(map(Decimal.from_float,values));mean=sum(rows)/len(rows)
        variance=sum((x-mean)**2 for x in rows)/(len(rows)-1)
        sd=variance.sqrt();sem=(variance/len(rows)).sqrt()
        width=Decimal(49)/25*sem
        return tuple(float(v) for v in (mean,sd,sem,mean-width,mean+width))


def independent_correlation(left,right):
    with localcontext() as context:
        context.prec=2200
        x=list(map(Decimal.from_float,left));y=list(map(Decimal.from_float,right))
        a=sum(x)/len(x);b=sum(y)/len(y)
        xx=sum((v-a)**2 for v in x);yy=sum((v-b)**2 for v in y)
        if not xx or not yy:return None
        return float(sum((u-a)*(v-b) for u,v in zip(x,y))/(xx*yy).sqrt())


def synthetic_assembly(method='latin-hypercube'):
    from rocket_propulsion.propulsion.burns import scale_tabulated_artifact
    from tests import test_burn_hierarchical_ensemble as artifact_fixture
    from tests import test_burn_uncertainty as physical_fixture
    fixture=physical_fixture.BlowdownUncertaintyTests();fixture.setUp()
    other=artifact_fixture.HierarchicalPropulsionEnsembleTests();other.setUp()
    base=other.discrete.realizations[0].artifact
    artifacts=[scale_tabulated_artifact(base,{'thrust_scale':scale},source_id=f'synthetic-summary:{i}',
                                       provenance={'scope':'synthetic interface assembly only'})
               for i,scale in enumerate((1.,1.125,1.25,1.375))]
    outputs=[SimpleNamespace(artifact=a,samples=(SimpleNamespace(supply_pressure_pa=1e6),),
                             cutoff_reason='synthetic-interface-fixture') for a in artifacts]
    with patch.object(model,'simulate_pressure_fed_blowdown',side_effect=outputs) as solve:
        result=model.simulate_blowdown_uncertainty(fixture.tank,fixture.point,fixture.law,
            fixture.parameters,sample_count=4,seed=91,method=method,maximum_mass_step_kg=.5)
    assert solve.call_count==4
    return result


class SummaryNumericalAttacks(unittest.TestCase):
    def statistic(self,values):
        try:return model._statistics('synthetic-kernel',values)
        except (DomainError,ArithmeticError) as error:
            self.fail(f'Representable sample summary was rejected: {error}')

    def convergence(self,values,tolerance):
        helper=getattr(model,'_sampling_convergence',None)
        self.assertTrue(callable(helper),'A shared exact-moment stability path is required')
        try:return helper(values,tolerance)
        except (DomainError,ArithmeticError) as error:
            self.fail(f'Representable stability diagnostic was rejected: {error}')

    def test_tiny_and_huge_dispersion_match_independent_decimal_moments(self):
        for values in ([0.,2e-200],[-1e200,1e200],[1e308,1e308],
                       [1.,nextafter(1.,2.)],[0.,1.,2.],[2*TINY,4*TINY],
                       [12.,37.],[12*TINY,37*TINY]):
            with self.subTest(values=values):
                expected=independent_moments(values);actual=self.statistic(values)
                self.assertEqual((actual.mean,actual.standard_deviation,actual.mean_standard_error,
                                  actual.mean_ci95_lower,actual.mean_ci95_upper),expected)

    def test_standard_error_is_not_computed_from_an_already_rounded_deviation(self):
        result=self.statistic([0.,1.,2.])
        self.assertEqual(result.mean_standard_error,independent_moments([0.,1.,2.])[2])

    def test_true_zero_is_retained_and_unrepresentable_mandatory_outputs_fail_typed(self):
        zero=self.statistic([0.,0.]);self.assertEqual(zero.mean,0.);self.assertEqual(zero.standard_deviation,0.)
        for values in ([0.,TINY],[-MAX,MAX],[0.,MAX]):
            with self.subTest(values=values),self.assertRaises(DomainError):
                model._statistics('mandatory-output',values)
        with self.assertRaises(DomainError) as caught:model._statistics('mandatory-percentile',[0.,2*TINY])
        self.assertEqual(caught.exception.field,'Sample quantile')
        with self.assertRaises(DomainError) as caught:model._statistics('mandatory-endpoint',[TINY,3*TINY])
        self.assertEqual(caught.exception.field,'Mean proxy lower endpoint')

    def test_linear_percentiles_recover_opposite_extremes_and_subnormal_endpoints(self):
        for values in ([-MAX,MAX],[0.,2*TINY],[-1e-300,1e-300]):
            for p in (Fraction(0),Fraction(1,20),Fraction(1,2),Fraction(19,20),Fraction(1)):
                exact=(1-p)*Fraction(values[0])+p*Fraction(values[1]);expected=float(exact)
                with self.subTest(values=values,p=p):
                    if exact and not expected:
                        with self.assertRaises(DomainError):model._percentile(values,p)
                    else:
                        try:actual=model._percentile(values,p)
                        except (DomainError,ArithmeticError) as error:self.fail(f'Finite percentile was lost: {error}')
                        self.assertEqual(actual,expected)

    def test_named_percentiles_use_rational_five_and_ninetyfive_percent_positions(self):
        values=[float(i) for i in range(7)]
        result=self.statistic(values)
        self.assertEqual(result.p05,float(Fraction(3,10)))
        self.assertEqual(result.p95,float(Fraction(57,10)))

    def test_summary_helpers_have_typed_cardinality_numeric_and_probability_boundaries(self):
        for values in ([],[1.],None,{'x':1.},[False,1.],['0',1.],[float('nan'),1.]):
            with self.subTest(values=values),self.assertRaises(DomainError):model._statistics('invalid',values)
        for p in (None,True,'0.5',-1.,1.01,float('inf')):
            with self.subTest(p=p),self.assertRaises(DomainError):model._percentile([0.,1.],p)
        with self.assertRaises(DomainError):model._percentile([], .5)

    def test_pearson_recovers_scale_invariance_and_signed_unit_correlations(self):
        for left,right in (([1e-200,2e-200,3e-200,4e-200],[1e200,2e200,3e200,4e200]),
                           ([1e-200,2e-200,3e-200,4e-200],[-1e200,-2e200,-3e200,-4e200])):
            expected=independent_correlation(left,right)
            try:actual=model._correlation(left,right)
            except (DomainError,ArithmeticError) as error:self.fail(f'Finite Pearson coefficient was lost: {error}')
            self.assertEqual(actual,expected)

    def test_nontrivial_pearson_coefficients_match_decimal_centered_products(self):
        for scale in (2.**-600,1.,2.**600):
            left=[x*scale for x in (-3.,-1.,2.,4.)]
            right=[x/scale for x in (1.,4.,-2.,3.)]
            try:actual=model._correlation(left,right)
            except (DomainError,ArithmeticError) as error:self.fail(f'Finite Pearson coefficient was lost: {error}')
            self.assertEqual(actual,independent_correlation(left,right))

    def test_pearson_centering_preserves_scatter_near_large_offsets(self):
        base=2.**900;step=2.**848
        left=[base+x*step for x in (0.,1.,3.,4.)]
        right=[-base+x*step for x in (4.,1.,2.,0.)]
        try:actual=model._correlation(left,right)
        except (DomainError,ArithmeticError) as error:self.fail(f'Representable offset correlation was lost: {error}')
        self.assertEqual(actual,independent_correlation(left,right))

    def test_nonzero_unrepresentable_pearson_is_not_fabricated_zero(self):
        with self.assertRaises(DomainError) as caught:
            model._correlation([1.,-1.,TINY,-TINY],[1.,1.,TINY,-TINY])
        self.assertEqual(caught.exception.field,'Pearson correlation')

    def test_undefined_pearson_is_not_reported_as_zero_linear_relationship(self):
        self.assertIsNone(model._correlation([1.,1.,1.],[1.,2.,3.]))
        self.assertIsNone(model._correlation([0.,0.,0.],[0.,0.,0.]))
        self.assertEqual(model._correlation([-1.,0.,1.],[1.,-2.,1.]),0.)

    def test_pearson_refuses_mismatched_empty_and_malformed_vectors(self):
        for left,right in (([],[]),([1.],[2.]),([1.,2.],[2.]),([True,2.],[1.,2.]),
                           (None,[1.,2.]),([1.,float('nan')],[1.,2.])):
            with self.subTest(left=left,right=right),self.assertRaises(DomainError):model._correlation(left,right)

    def test_relative_diagnostics_are_symmetric_scale_invariant_and_finite(self):
        for left,right,expected in ((1e-200,2e-200,.5),(2e-200,1e-200,.5),
                                    (-MAX,MAX,2.),(0.,0.,0.),(0.,1.,1.)):
            self.assertEqual(model._relative_change(left,right),expected)

    def test_half_full_diagnostic_retains_changes_hidden_by_rounded_means(self):
        values=[1.,1.,1.,nextafter(1.,2.)]
        result=self.convergence(values,.1)
        exact_mean=sum(map(Fraction,values))/4
        expected=float((exact_mean-1)/exact_mean)
        self.assertGreater(expected,0.)
        self.assertEqual(result.mean_relative_change,expected)
        self.assertEqual(result.standard_deviation_relative_change,1.)
        self.assertFalse(result.converged)

    def test_stability_threshold_uses_exact_values_even_at_a_display_tie(self):
        with localcontext() as context:
            context.prec=200
            exact=1-(Decimal(2)/3).sqrt()
            tolerance=float(exact)
            expected=exact<=Decimal.from_float(tolerance)
        self.assertFalse(expected,'This reference deliberately rounds the exact change down to the threshold')
        result=self.convergence([0.,2.,0.,2.],tolerance)
        self.assertEqual(result.standard_deviation_relative_change,tolerance)
        self.assertEqual(result.converged,expected)
        constant=self.convergence([1.,1.,1.,1.],.1)
        self.assertEqual((constant.mean_relative_change,constant.standard_deviation_relative_change),(0.,0.))
        self.assertTrue(constant.converged)

    def test_interval_is_explicitly_an_unqualified_normal_formula_proxy(self):
        result=self.statistic([0.,1.,2.])
        self.assertEqual(getattr(result,'mean_interval_kind',None),'iid_normal_formula_proxy_unqualified_coverage')
        self.assertIs(getattr(result,'coverage_established',None),False)


class SummaryRecordAttacks(unittest.TestCase):
    @classmethod
    def setUpClass(cls):cls.result=synthetic_assembly()

    def test_result_retains_normalized_sampler_and_summary_policy(self):
        self.assertEqual(getattr(self.result,'summary_numerics',None),
            'exact_binary_sample_moments_linear_rational_quantiles_and_pearson_v1')
        self.assertTrue(any('dependent' in w.lower() and 'interval' in w.lower() for w in self.result.warnings))

    def test_metric_records_refuse_booleans_negative_dispersion_and_inverted_quantiles(self):
        statistic=self.result.statistics[0]
        for changes in ({'mean':True},{'standard_deviation':True},{'standard_deviation':-1e-300},
                        {'coverage_established':True},{'mean_interval_kind':'qualified_95_percent_interval'},
                        {'p50':statistic.maximum+1.},{'mean_ci95_lower':statistic.mean+1.}):
            with self.subTest(changes=changes),self.assertRaises(DomainError):replace(statistic,**changes)

    def test_sensitivity_records_distinguish_undefined_constant_output_and_rank_it_last(self):
        records=[r for r in self.result.sensitivity if r.metric=='final_supply_pressure_pa']
        self.assertEqual(len(records),len(self.result.parameters))
        for row in records:
            self.assertIsNone(row.correlation);self.assertIsNone(row.absolute_correlation)
            self.assertIsNone(row.rank)
            self.assertEqual(getattr(row,'undefined_reason',None),'constant_input_or_output')

    def test_sensitivity_records_cannot_forge_absolute_value_or_rank(self):
        row=next(x for x in self.result.sensitivity if x.correlation is not None)
        for changes in ({'correlation':True},{'absolute_correlation':.5},{'rank':0}):
            with self.subTest(changes=changes),self.assertRaises(DomainError):replace(row,**changes)

    def test_convergence_records_refuse_boolean_counts_and_false_stability_claims(self):
        row=self.result.convergence
        for changes in ({'half_sample_count':True},{'full_sample_count':row.full_sample_count+2},
                        {'mean_relative_change':-1e-300},{'converged':not row.converged}):
            with self.subTest(changes=changes),self.assertRaises(DomainError):replace(row,**changes)

    def test_retained_samples_snapshot_input_pairs_and_reproduce_artifact_metrics(self):
        row=self.result.samples[0];pairs=[list(p) for p in row.inputs]
        owned=replace(row,inputs=pairs);pairs[0][1]=999.
        self.assertEqual(owned.inputs,row.inputs)
        for changes in ({'index':True},{'duration_s':nextafter(row.duration_s,float('inf'))},
                        {'inputs':row.inputs+row.inputs[:1]}):
            with self.subTest(changes=changes),self.assertRaises(DomainError):replace(row,**changes)

    def test_result_snapshots_samples_and_refuses_incomplete_or_duplicate_records(self):
        rows=list(self.result.samples);owned=replace(self.result,samples=rows);rows.pop()
        self.assertEqual(owned.samples,self.result.samples)
        parameters=list(self.result.parameters);owned=replace(self.result,parameters=parameters);parameters.pop()
        self.assertEqual(owned.parameters,self.result.parameters)
        for changes in ({'sample_count':True},{'sample_count':self.result.sample_count+1},
                        {'samples':self.result.samples[:-1]},
                        {'statistics':self.result.statistics[:-1]},
                        {'sensitivity':self.result.sensitivity+self.result.sensitivity[:1]},
                        {'summary_numerics':'unsupported'}):
            with self.subTest(changes=changes),self.assertRaises(DomainError):replace(self.result,**changes)

    def test_result_reordering_preserves_semantic_indices_rank_order_and_stability(self):
        other=replace(self.result,samples=tuple(reversed(self.result.samples)),
            statistics=tuple(reversed(self.result.statistics)),sensitivity=tuple(reversed(self.result.sensitivity)))
        self.assertEqual(other,self.result)

    def test_result_recomputes_statistics_and_sensitivities_from_retained_samples(self):
        row=next(r for r in self.result.statistics if r.metric=='delivered_impulse_n_s')
        changed=replace(row,mean=nextafter(row.mean,float('inf')))
        rows=tuple(changed if r is row else r for r in self.result.statistics)
        with self.assertRaises(DomainError):replace(self.result,statistics=rows)
        sensitivity=next(r for r in self.result.sensitivity if r.correlation not in (None,0.,1.,-1.))
        r=nextafter(sensitivity.correlation,0.)
        changed=replace(sensitivity,correlation=r,absolute_correlation=abs(r))
        rows=tuple(changed if x is sensitivity else x for x in self.result.sensitivity)
        with self.assertRaises(DomainError):replace(self.result,sensitivity=rows)
        undefined=next(r for r in self.result.sensitivity if r.correlation is None)
        forged=replace(undefined,correlation=0.,absolute_correlation=0.,rank=1,undefined_reason=None)
        rows=tuple(forged if x is undefined else x for x in self.result.sensitivity)
        with self.assertRaises(DomainError):replace(self.result,sensitivity=rows)

    def test_one_ulp_nominal_configuration_mismatch_is_rejected_before_solving(self):
        from tests import test_burn_uncertainty as fixture_module
        fixture=fixture_module.BlowdownUncertaintyTests();fixture.setUp()
        values=tuple(replace(p,nominal=nextafter(p.nominal,float('inf')))
                     if p.name=='rated_delivered_thrust_n' else p for p in fixture.parameters)
        with patch.object(model,'simulate_pressure_fed_blowdown',side_effect=AssertionError('Mismatched nominal reached solver')) as solve:
            with self.assertRaises(DomainError):
                model.simulate_blowdown_uncertainty(fixture.tank,fixture.point,fixture.law,values,sample_count=4)
            solve.assert_not_called()

    def test_consistently_relabelled_unknown_blowdown_parameter_is_not_accepted(self):
        old=self.result.parameters[0].name;new='unsupported_physical_parameter'
        parameters=(replace(self.result.parameters[0],name=new),)+self.result.parameters[1:]
        samples=tuple(replace(row,inputs=tuple((new if key==old else key,value) for key,value in row.inputs))
                      for row in self.result.samples)
        sensitivity=tuple(replace(row,parameter=new) if row.parameter==old else row for row in self.result.sensitivity)
        with self.assertRaises(DomainError):replace(self.result,parameters=parameters,samples=samples,sensitivity=sensitivity)

    def test_public_controls_fail_typed_before_any_solver_call(self):
        from tests import test_burn_uncertainty as fixture_module
        fixture=fixture_module.BlowdownUncertaintyTests();fixture.setUp()
        for arguments in ({'convergence_tolerance':True},{'convergence_tolerance':None},
                          {'reserve_mass_kg':True},{'sample_count':4.5},
                          {'parameters':None},{'method':None}):
            kwargs={'parameters': fixture.parameters,'sample_count': 4,'seed': 1};kwargs.update(arguments)
            with patch.object(model,'simulate_pressure_fed_blowdown',side_effect=AssertionError('Invalid control reached solver')) as solve:
                with self.subTest(arguments=arguments),self.assertRaises(DomainError):
                    model.simulate_blowdown_uncertainty(fixture.tank,fixture.point,fixture.law,**kwargs)
                solve.assert_not_called()
