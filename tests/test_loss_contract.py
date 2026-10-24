"""Independent loss ledgers, declared operating bases and covariance attacks."""
import json
import unittest
from dataclasses import asdict, replace
from decimal import Decimal, localcontext
from fractions import Fraction
from math import inf, nextafter
from sys import float_info

from rocket_propulsion.core.errors import DomainError, ConvergenceError
from rocket_propulsion.propulsion import calculate_loss_budget
from rocket_propulsion.api.routes import dispatch_calculation, calculation_metadata

FLAGS = ('enable_combustion_loss', 'enable_discharge_loss',
         'enable_boundary_layer_loss', 'enable_divergence_loss')
EFFS = ('combustion_efficiency', 'discharge_coefficient',
        'boundary_layer_efficiency', 'divergence_efficiency')
IDENTITY = [[int(i == j) for j in range(4)] for i in range(4)]


def inputs(**changes):
    return dict(dict(ideal_thrust_n=100000., ideal_specific_impulse_s=450.,
                     relative_uncertainty=0.), **changes)


def disabled(**changes):
    return inputs(**dict({name: False for name in FLAGS}, **changes))


def oracle(x, precision=1400):
    """Separate Decimal state transitions and analytic log-sensitivity covariance."""
    with localcontext() as c:
        c.prec = precision
        d = lambda v: Decimal.from_float(float(v))
        f0, i0 = d(x['ideal_thrust_n']), d(x['ideal_specific_impulse_s'])
        eff = [d(x.get(k, v)) for k, v in zip(EFFS, (.97, .98, .985, .98))]
        enabled = [x.get(k, True) for k in FLAGS]
        physical = x.get('operating_basis') == 'fixed-pressure-area'
        f, isp, flow = f0, i0, Decimal(1)
        terms = []
        for index, (eta, active) in enumerate(zip(eff, enabled)):
            before_f, before_i = f, isp
            if active:
                if physical:
                    if index == 0:
                        isp *= eta; flow /= eta
                    elif index == 1:
                        isp /= eta; flow *= eta
                    else:
                        f *= eta; isp *= eta
                else:
                    f *= eta; isp *= eta
            terms.append((float(before_f - f), float(before_i - isp)))
        matrix = x.get('correlation_matrix', IDENTITY)
        elastic = ([0, 0, 1, 1], [1, -1, 1, 1], [-1, 1, 0, 0]) if physical else ([1]*4, [1]*4, [0]*4)
        variances = []
        for row in elastic:
            e = [v if active else 0 for v, active in zip(row, enabled)]
            variances.append(d(x.get('relative_uncertainty', .02))**2 * sum(
                (d(matrix[a][b])*e[a]*e[b] for a in range(4) for b in range(4)), Decimal(0)))
        sf, si, sm = [v.sqrt() for v in variances]
        return dict(delivered_thrust_n=float(f), delivered_specific_impulse_s=float(isp),
            total_efficiency=float(f/f0), specific_impulse_ratio=float(isp/i0), mass_flow_ratio=float(flow),
            total_thrust_loss_n=float(f0-f), total_specific_impulse_loss_s=float(i0-isp),
            one_sigma_relative=float(sf), specific_impulse_one_sigma_relative=float(si),
            mass_flow_one_sigma_relative=float(sm),
            delivered_thrust_band_n=(float(f*(1-sf)), float(f*(1+sf))),
            delivered_specific_impulse_band_s=(float(isp*(1-si)), float(isp*(1+si))), terms=terms)


class LossAttacks(unittest.TestCase):
    def valid(self, x):
        try:
            return calculate_loss_budget(**x)
        except (DomainError, ConvergenceError, ArithmeticError) as error:
            self.fail(f'Representable loss-budget witness rejected: {error}')

    def check(self, x):
        a, b = oracle(x), oracle(x, 2200)
        self.assertEqual(a, b, 'Independent reference precision did not settle.')
        result = self.valid(x)
        for field, expected in a.items():
            if field == 'terms':
                self.assertEqual([(t.thrust_loss_n, t.specific_impulse_loss_s) for t in result.terms], expected)
            else:
                self.assertEqual(getattr(result, field, None), expected, field)
        return result

    def test_real_inputs_reject_booleans_text_nonfinite_and_erased_values(self):
        for field in ('ideal_thrust_n', 'ideal_specific_impulse_s', *EFFS, 'relative_uncertainty'):
            for bad in (True, '1', None, [], float('nan'), inf, 10**10000, Fraction(1, 10**1000)):
                with self.subTest(field=field, bad=type(bad).__name__):
                    with self.assertRaises(DomainError):
                        calculate_loss_budget(**inputs(**{field: bad}))

    def test_switches_are_actual_booleans(self):
        for field in FLAGS:
            for bad in (0, 1, 2, 'false', [], None):
                with self.subTest(field=field, bad=bad):
                    with self.assertRaises(DomainError):
                        calculate_loss_budget(**inputs(**{field: bad}))

    def test_disabled_factors_are_still_validated(self):
        for field in EFFS:
            for bad in (0., -1., nextafter(1., inf)):
                with self.subTest(field=field, bad=bad):
                    with self.assertRaises(DomainError):
                        calculate_loss_budget(**disabled(**{field: bad}))

    def test_representable_small_loss_is_not_rounded_thrust_subtraction(self):
        x = disabled(ideal_thrust_n=1e308, ideal_specific_impulse_s=1.,
                     enable_combustion_loss=True, combustion_efficiency=nextafter(1., 0.))
        result = self.valid(x)
        exact = float(Fraction(1e308) * (1-Fraction(x['combustion_efficiency'])))
        self.assertEqual(result.total_thrust_loss_n, exact)
        self.assertEqual(result.terms[0].thrust_loss_n, exact)
        self.assertEqual(result.thrust_balance_residual_n, 0.)

    def test_all_switch_combinations_follow_exact_sequential_ledgers(self):
        for basis in ('fixed-mass-flow-lumped', 'fixed-pressure-area'):
            for mask in range(16):
                with self.subTest(basis=basis, mask=mask):
                    self.check(inputs(operating_basis=basis, **dict(zip(FLAGS, [bool(mask & (1<<i)) for i in range(4)]))))

    def test_no_enabled_uncertainty_sources_means_zero_factor_uncertainty(self):
        result = self.valid(disabled(ideal_thrust_n=float_info.max, relative_uncertainty=.25))
        self.assertEqual(result.one_sigma_relative, 0.)
        self.assertEqual(result.delivered_thrust_band_n, (float_info.max, float_info.max))

    def test_mandatory_band_overflow_is_typed(self):
        with self.assertRaises(DomainError):
            calculate_loss_budget(**disabled(ideal_thrust_n=float_info.max, relative_uncertainty=.25,
                                             combustion_efficiency=1., enable_combustion_loss=True))

    def test_nonzero_unrepresentable_loss_is_not_silently_erased(self):
        with self.assertRaises(DomainError):
            calculate_loss_budget(**disabled(ideal_thrust_n=5e-324, enable_combustion_loss=True,
                                             combustion_efficiency=nextafter(1., 0.)))

    def test_subnormal_uncertainty_is_not_squared_to_zero(self):
        self.check(disabled(relative_uncertainty=5e-324, combustion_efficiency=1., enable_combustion_loss=True))

    def test_physical_discharge_changes_flow_and_isp_without_forced_thrust_loss(self):
        result = self.check(disabled(operating_basis='fixed-pressure-area', enable_discharge_loss=True,
                                     discharge_coefficient=.8))
        self.assertEqual(result.delivered_thrust_n, result.ideal_thrust_n)
        self.assertLess(result.total_specific_impulse_loss_s, 0.)
        self.assertTrue(any(w.code == 'coefficient_compatibility_not_validated' for w in result.warnings))

    def test_physical_combustion_changes_cstar_and_flow_at_fixed_pressure_area(self):
        result = self.check(disabled(operating_basis='fixed-pressure-area', enable_combustion_loss=True,
                                     combustion_efficiency=.8, relative_uncertainty=.1))
        self.assertEqual(result.one_sigma_relative, 0.)
        self.assertEqual(result.specific_impulse_one_sigma_relative, .1)
        self.assertEqual(result.mass_flow_one_sigma_relative, .1)

    def test_recoverable_intrinsic_cstar_discharge_ratio_uses_unexported_intermediates(self):
        self.check(disabled(operating_basis='fixed-pressure-area', enable_combustion_loss=True,
            enable_discharge_loss=True, combustion_efficiency=5e-324, discharge_coefficient=5e-324))

    def test_covariance_uses_signed_elasticities_and_allows_exact_cancellation(self):
        matrix = [[1]*4 for _ in range(4)]
        self.check(inputs(operating_basis='fixed-pressure-area', correlation_matrix=matrix,
                          relative_uncertainty=.25))
        result = self.check(disabled(operating_basis='fixed-pressure-area', correlation_matrix=matrix,
            enable_combustion_loss=True, enable_discharge_loss=True, relative_uncertainty=.25))
        self.assertEqual(result.specific_impulse_one_sigma_relative, 0.)
        self.assertEqual(result.mass_flow_one_sigma_relative, 0.)

    def test_perfect_correlation_can_produce_exact_zero_lower_band(self):
        result = self.check(inputs(correlation_matrix=[[1]*4 for _ in range(4)], relative_uncertainty=.25))
        self.assertEqual(result.one_sigma_relative, 1.)
        self.assertEqual(result.delivered_thrust_band_n[0], 0.)

    def test_partial_covariance_distinguishes_thrust_isp_and_flow(self):
        matrix = [[1., .5, 0., 0.], [.5, 1., 0., 0.], [0., 0., 1., -.25], [0., 0., -.25, 1.]]
        result = self.check(inputs(operating_basis='fixed-pressure-area', correlation_matrix=matrix,
                                   relative_uncertainty=.1))
        self.assertNotEqual(result.one_sigma_relative, result.specific_impulse_one_sigma_relative)
        self.assertNotEqual(result.one_sigma_relative, result.mass_flow_one_sigma_relative)

    def test_term_replay_rejects_forged_labels_sources_switches_and_policy(self):
        term = self.valid(inputs()).terms[0]
        for changes in ({'thrust_loss_n': nextafter(term.thrust_loss_n, inf)},
                        {'specific_impulse_loss_s': nextafter(term.specific_impulse_loss_s, inf)},
                        {'label': 'validated hardware'}, {'source': 'NASA certified'},
                        {'enabled': 1}, {'efficiency': nextafter(term.efficiency, inf)},
                        {'key': 'discharge'}, {'numerical_policy': 'unknown'}):
            with self.subTest(changes=changes):
                with self.assertRaises(DomainError):
                    replace(term, **changes)

    def test_bands_warnings_and_boolean_qualification_cannot_be_forged(self):
        result = self.valid(inputs(relative_uncertainty=.02))
        for changes in ({'delivered_thrust_band_n': (result.delivered_thrust_band_n[1], result.delivered_thrust_band_n[0])},
                        {'delivered_specific_impulse_band_s': (1.,)}, {'warnings': ()},
                        {'warnings': (replace(result.warnings[0], severity='info'), *result.warnings[1:])},
                        {'coverage_qualified': 0}, {'one_sigma_relative': True},
                        {'specific_impulse_ratio': nextafter(result.specific_impulse_ratio, inf)},
                        {'mass_flow_ratio': nextafter(result.mass_flow_ratio, inf)}):
            with self.subTest(changes=list(changes)):
                with self.assertRaises(DomainError):
                    replace(result, **changes)

    def test_retained_warning_text_drops_comparison_callbacks_before_replay(self):
        class EqualText(str):
            def __eq__(self, other):
                return True
            __hash__ = str.__hash__
        result = self.valid(inputs())
        for field in ('code', 'message', 'severity', 'field'):
            forged = replace(result.warnings[0], **{field: EqualText('fabricated qualification')})
            for warning in (forged, asdict(forged)):
                with self.subTest(field=field, mapping=isinstance(warning, dict)):
                    with self.assertRaises(DomainError):
                        replace(result, warnings=(warning, *result.warnings[1:]))
        original = result.warnings[0]
        valid = replace(original, code=EqualText(original.code), message=EqualText(original.message),
                        severity=EqualText(original.severity))
        reopened = replace(result, warnings=(valid, *result.warnings[1:]))
        for field in ('code', 'message', 'severity'):
            self.assertIs(type(getattr(reopened.warnings[0], field)), str)

    def test_legacy_default_is_numerically_identical_to_explicit_legacy_basis(self):
        self.assertEqual(self.valid(inputs()), self.valid(inputs(operating_basis='fixed-mass-flow-lumped')))

    def test_checked_exports_keep_extreme_but_representable_factor_products(self):
        self.check(disabled(ideal_thrust_n=float_info.max, ideal_specific_impulse_s=float_info.max,
                            combustion_efficiency=5e-324, enable_combustion_loss=True))

    def test_indefinite_and_one_ulp_outside_psd_are_rejected_without_tolerance(self):
        bad = [[1., .9, .9, 0.], [.9, 1., -.9, 0.], [.9, -.9, 1., 0.], [0., 0., 0., 1.]]
        r = nextafter(-.5, -inf)
        near = [[1., r, r, 0.], [r, 1., r, 0.], [r, r, 1., 0.], [0., 0., 0., 1.]]
        zero_pivot = [[1., 1., 0., 0.], [1., 1., .1, 0.], [0., .1, 1., 0.], [0., 0., 0., 1.]]
        for matrix in (bad, near, zero_pivot):
            with self.subTest(matrix=matrix):
                with self.assertRaises(DomainError):
                    calculate_loss_budget(**inputs(correlation_matrix=matrix))

    def test_correlation_shape_symmetry_diagonal_and_types_are_strict(self):
        cases = [[], [[1]], 'identity', None]
        for i, j, value in ((0, 1, .1), (0, 0, nextafter(1., 0.)), (0, 1, True), (0, 1, inf)):
            matrix = [row[:] for row in IDENTITY]; matrix[i][j] = value; cases.append(matrix)
        for matrix in cases:
            with self.subTest(matrix=matrix):
                with self.assertRaises(DomainError):
                    calculate_loss_budget(**inputs(correlation_matrix=matrix))

    def test_result_replay_rejects_scalar_term_policy_scope_and_input_forgeries(self):
        result = self.valid(inputs())
        for field in ('delivered_thrust_n', 'delivered_specific_impulse_s', 'total_efficiency',
                      'total_thrust_loss_n', 'one_sigma_relative'):
            with self.subTest(field=field):
                with self.assertRaises(DomainError):
                    replace(result, **{field: nextafter(getattr(result, field), inf)})
        self.assertTrue(hasattr(result, 'inputs'), 'Result must retain its operating inputs.')
        for changes in ({'numerical_policy': 'unknown'}, {'coverage_qualified': True},
                        {'uncertainty_kind': 'qualified'}, {'assumptions': ('validated hardware',)},
                        {'inputs': replace(result.inputs, combustion_efficiency=.9)},
                        {'terms': tuple(reversed(result.terms))}):
            with self.subTest(changes=list(changes)):
                with self.assertRaises(DomainError):
                    replace(result, **changes)

    def test_json_reopen_snapshots_matrix_inputs_terms_and_warnings(self):
        result = self.valid(inputs(correlation_matrix=[row[:] for row in IDENTITY]))
        document = json.loads(json.dumps(asdict(result), allow_nan=False))
        reopened = type(result)(**document)
        self.assertEqual(reopened, result)
        document['inputs']['correlation_matrix'][0][0] = 0.
        document['terms'][0]['label'] = 'changed'
        self.assertEqual(reopened, result)
        with self.assertRaises(DomainError):
            type(result)(**document)

    def test_basis_and_coverage_claims_are_explicit_and_algebra_is_not_validation(self):
        result = self.valid(inputs())
        self.assertEqual(getattr(result, 'operating_basis', None), 'fixed-mass-flow-lumped')
        self.assertIs(getattr(result, 'coverage_qualified', None), False)
        self.assertEqual(getattr(result, 'uncertainty_kind', None), 'first_order_relative_factor_covariance_proxy')
        self.assertTrue(any(w.code == 'legacy_discharge_is_performance_multiplier' for w in result.warnings))
        for value in ('unknown', True, 1, None):
            with self.subTest(value=value):
                with self.assertRaises(DomainError):
                    calculate_loss_budget(**inputs(operating_basis=value))


class LossApiAttacks(unittest.TestCase):
    def test_dispatch_accepts_explicit_basis_and_matrix_and_exports_version_two(self):
        x = disabled(operating_basis='fixed-pressure-area', enable_discharge_loss=True, discharge_coefficient=.8,
                     correlation_matrix=IDENTITY)
        result = dispatch_calculation('/api/v1/loss-budget', x)
        self.assertEqual(result.get('mass_flow_ratio'), .8)
        self.assertEqual(result.get('operating_basis'), 'fixed-pressure-area')
        self.assertEqual(calculation_metadata('/api/v1/loss-budget')['version'], '2.0')
        json.dumps(result, allow_nan=False)

    def test_dispatch_does_not_coerce_invalid_basis_or_hide_invalid_matrix(self):
        for change in ({'operating_basis': True}, {'correlation_matrix': [[1]]}):
            with self.subTest(change=change):
                with self.assertRaises(DomainError):
                    dispatch_calculation('/api/v1/loss-budget', inputs(**change))
