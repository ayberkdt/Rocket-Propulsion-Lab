"""Independent mass ledgers, exhaustive grid oracles and hostile staging records."""
import unittest
from dataclasses import asdict, replace
from decimal import Decimal, localcontext
from math import inf, nextafter
from unittest.mock import patch

from rocket_propulsion.api import routes as api_routes
from rocket_propulsion.api.routes import dispatch_calculation
from rocket_propulsion.core.constants import STANDARD_GRAVITY_M_S2
from rocket_propulsion.core.errors import ConvergenceError, DomainError
from rocket_propulsion.propulsion import staging_numerics
from rocket_propulsion.propulsion.staging import StageDefinition, optimize_staging


def case(**changes):
    return dict({'stages': (StageDefinition('S1', 300., .1),),
                     'target_delta_v_m_s': 100., 'payload_mass_kg': 1., 'resolution': 40}, **changes)


def reference(x, units, precision=1400):
    """Direct exponential rocket equation, independent of production helpers."""
    with localcontext() as c:
        c.prec = precision
        d = lambda v: Decimal.from_float(float(v))
        upper = d(x['payload_mass_kg']); result = []
        for stage, count in zip(reversed(x['stages']), reversed(units)):
            dv = d(x['target_delta_v_m_s'])*count/x['resolution']
            ratio = (dv/(d(STANDARD_GRAVITY_M_S2)*d(stage.specific_impulse_s))).exp()
            eps = d(stage.structural_fraction)
            if 1-eps*ratio <= 0:
                return None
            wet = upper*(ratio-1)/(1-eps*ratio)
            result.append({'delta_v_m_s': float(dv), 'mass_ratio': float(ratio),
                'initial_mass_kg': float(upper+wet), 'upper_stack_mass_kg': float(upper),
                'dry_mass_kg': float(eps*wet), 'propellant_mass_kg': float((1-eps)*wet)})
            upper += wet
        return upper, tuple(reversed(result))


def compositions(total, count):
    if count == 1:
        yield (total,)
    else:
        for first in range(total+1):
            for tail in compositions(total-first, count-1):
                yield (first,)+tail


def exhaustive(x):
    """Enumerate every grid allocation; compare unrounded 180-digit masses."""
    best = None
    for units in compositions(x['resolution'], len(x['stages'])):
        r = reference(x, units, 180)
        if r is None:
            continue
        # Smaller upper-stage allocations win an exactly equal objective.
        key = r[0], tuple(reversed(units))
        if best is None or key < best[0]:
            best = key, units
    return None if best is None else best[1]


class StagingAttacks(unittest.TestCase):
    def valid(self, x):
        try:
            return optimize_staging(**x)
        except (DomainError, ConvergenceError, ArithmeticError) as error:
            self.fail(f'Representable staging witness rejected: {error}')

    def check(self, x, units):
        r = self.valid(x)
        a, b = reference(x, units), reference(x, units, 2200)
        self.assertEqual((float(a[0]), a[1]), (float(b[0]), b[1]))
        self.assertEqual(r.initial_mass_kg, float(a[0]))
        with localcontext() as c:
            c.prec = 1400
            self.assertEqual(r.payload_fraction,
                             float(Decimal.from_float(float(x['payload_mass_kg']))/a[0]))
        for actual, expected in zip(r.stages, a[1]):
            for name, value in expected.items():
                self.assertEqual(getattr(actual, name), value, name)
        self.assertEqual(r.allocation_units, units)
        self.assertEqual(r.achieved_delta_v_m_s, x['target_delta_v_m_s'])
        return r

    def test_boolean_target_is_not_one_metre_per_second(self):
        with self.assertRaises(DomainError): optimize_staging(**case(target_delta_v_m_s=True))

    def test_boolean_payload_is_not_one_kilogram(self):
        with self.assertRaises(DomainError): optimize_staging(**case(payload_mass_kg=True))

    def test_stage_scalars_reject_booleans_with_numeric_positive_control(self):
        self.valid(case(stages=(StageDefinition('control', 1., .1),), target_delta_v_m_s=1.))
        for isp, eps in ((True, .1), (300., True)):
            with self.subTest(isp=isp, eps=eps), self.assertRaises(DomainError):
                optimize_staging(**case(stages=(StageDefinition('bad', isp, eps),), target_delta_v_m_s=1.))

    def test_resolution_is_an_actual_integer(self):
        for value in (40., True, '40', None, 39, 801):
            with self.subTest(value=value), self.assertRaises(DomainError):
                optimize_staging(**case(resolution=value))

    def test_malformed_scalars_and_names_are_typed(self):
        for field in ('target_delta_v_m_s', 'payload_mass_kg'):
            for value in (None, '1', inf, float('nan'), 0, -1):
                with self.subTest(field=field, value=value), self.assertRaises(DomainError):
                    optimize_staging(**case(**{field: value}))
        for name in (None, [], 42, '', ' ', '\ud800'):
            with self.subTest(name=repr(name)), self.assertRaises(DomainError):
                optimize_staging(**case(stages=(StageDefinition(name, 300., .1),)))

    def test_oversized_stage_iterables_have_a_bounded_read(self):
        stage = StageDefinition('same', 300., .1)
        def endless():
            for _ in range(6): yield stage
            raise AssertionError('Read beyond the stage-count rejection witness')
        with self.assertRaises(DomainError): optimize_staging(**case(stages=endless()))

    def test_large_isp_keeps_recoverable_positive_propellant(self):
        x = case(stages=(StageDefinition('S', 1e308, .1),), target_delta_v_m_s=1e6, payload_mass_kg=1e308)
        r = self.check(x, (40,))
        self.assertEqual(r.stages[0].propellant_mass_kg, 101971.62129779284)
        self.assertEqual(r.stages[0].mass_ratio, 1.)

    def test_last_float_below_capacity_is_not_rejected_by_rounded_ratio(self):
        x = case(stages=(StageDefinition('S', 1., .5),), target_delta_v_m_s=6.797451798238187,
                 payload_mass_kg=1e-200)
        r = self.check(x, (40,))
        self.assertEqual(r.initial_mass_kg, 3.1143148717373296e-183)
        self.assertEqual(r.stages[0].mass_ratio, 2.)

    def test_next_float_above_capacity_is_infeasible(self):
        with self.assertRaises(DomainError):
            optimize_staging(**case(stages=(StageDefinition('S', 1., .5),),
                target_delta_v_m_s=nextafter(6.797451798238187, inf), payload_mass_kg=1e-200))

    def test_unresolved_capacity_enclosure_fails_without_false_infeasibility(self):
        from fractions import Fraction
        stage = StageDefinition('S', 1., .5); target = 6.797451798238187
        # A conservative enclosure containing the independently computed true
        # survival straddles capacity; its upper bound rounds to epsilon.
        bounds = Decimal('.5'), Decimal('.50000000000000001')
        with localcontext() as c:
            c.prec = 2200
            actual = (-Decimal.from_float(target)/Decimal.from_float(STANDARD_GRAVITY_M_S2)).exp()
            self.assertLess(bounds[0], actual); self.assertLess(actual, bounds[1])
        staging_numerics._feasible.cache_clear()
        try:
            with (
                patch.object(staging_numerics, '_exp_interval', return_value=bounds),
                self.assertRaises(ConvergenceError),
            ):
                staging_numerics._feasible(stage, Fraction(target))
        finally:
            staging_numerics._feasible.cache_clear()

    def test_grid_infeasibility_does_not_claim_continuous_infeasibility(self):
        stages = (StageDefinition('a', 1., .5), StageDefinition('b', 1., .5))
        with self.assertRaises(DomainError) as caught:
            optimize_staging(**case(stages=stages, target_delta_v_m_s=13.5, resolution=41))
        self.assertEqual(caught.exception.code, 'infeasible_staging_grid')
        # Even grid has a feasible equal split of exactly the same total delta-v.
        self.check(case(stages=stages, target_delta_v_m_s=13.5), (20, 20))

    def test_unrepresentable_mandatory_outputs_fail_typed(self):
        with self.assertRaises(DomainError):
            optimize_staging(**case(target_delta_v_m_s=5000., payload_mass_kg=1e308))

    def test_two_stage_global_grid_optimum_matches_exhaustive_search(self):
        x = case(stages=(StageDefinition('lower', 300., .1), StageDefinition('upper', 450., .08)),
                 target_delta_v_m_s=9400., payload_mass_kg=1000.)
        self.check(x, exhaustive(x))

    def test_three_stage_global_grid_optimum_matches_exhaustive_search(self):
        x = case(stages=(StageDefinition('a', 260., .12), StageDefinition('b', 330., .09),
                         StageDefinition('c', 450., .08)), target_delta_v_m_s=8000., payload_mass_kg=23.)
        self.check(x, exhaustive(x))

    def test_identical_stage_ties_follow_lower_stage_first_policy(self):
        x = case(stages=(StageDefinition('a', 300., .1), StageDefinition('b', 300., .1)), resolution=41)
        self.check(x, (21, 20))

    def test_five_identical_stages_have_equal_allocation(self):
        x = case(stages=tuple(StageDefinition(str(i), 300., .1) for i in range(5)), target_delta_v_m_s=9000.)
        self.check(x, (8,)*5)

    def test_float_unit_mass_factors_do_not_erase_optimizer_costs(self):
        x = case(stages=(StageDefinition('less efficient', 1e307, .1),
                         StageDefinition('more efficient', 1e308, .1)),
                 target_delta_v_m_s=1e6, payload_mass_kg=1e308)
        self.check(x, (0, 40))

    def test_nearly_equal_isps_choose_the_cheaper_odd_unit(self):
        x = case(stages=(StageDefinition('a', 300., .1),
                         StageDefinition('b', nextafter(300., inf), .1)), resolution=41)
        self.check(x, exhaustive(x))
        self.assertEqual(self.valid(x).allocation_units, (20, 21))

    def test_unrounded_upper_stack_is_used_in_lower_stage_ledger(self):
        x = case(stages=(StageDefinition('a', 308.73502360600946, .052608970539099305),
                         StageDefinition('b', 308.73502360600946, .052608970539099305)),
                 target_delta_v_m_s=1123.3228889984566, payload_mass_kg=412.35628019652785)
        r = self.check(x, (20, 20))
        self.assertEqual(r.stages[0].propellant_mass_kg, 103.51419346520541)
        self.assertNotEqual(r.stages[0].propellant_mass_kg, 103.5141934652054,
                            'Rounded upper-stack arithmetic produces this neighbouring value.')

    def test_result_replay_rejects_forged_mass_and_qualification(self):
        r = self.valid(case())
        for changes in ({'initial_mass_kg': 1e99}, {'feasible': False}, {'feasible': 1},
                        {'achieved_delta_v_m_s': 101.}, {'payload_fraction': .5},
                        {'numerical_policy': 'other'}, {'continuous_optimum_qualified': True}):
            with self.subTest(changes=changes), self.assertRaises(DomainError): replace(r, **changes)

    def test_stage_record_replay_rejects_forged_propellant(self):
        r = self.valid(case())
        with self.assertRaises(DomainError): replace(r.stages[0], propellant_mass_kg=123.)

    def test_retained_literal_text_callbacks_cannot_disguise_policy_forgery(self):
        class EqualText(str):
            def __eq__(self, other): return True
            def __ne__(self, other): return False
        r = self.valid(case())
        for name in ('numerical_policy', 'optimization_scope', 'tie_break_policy'):
            with self.subTest(field=name), self.assertRaises(DomainError):
                replace(r, **{name: EqualText('forged')})

    def test_retained_allocation_iterables_have_a_bounded_read(self):
        r = self.valid(case())
        def endless():
            yield 40; yield 40
            raise AssertionError('Read beyond allocation cardinality rejection')
        with self.assertRaises(DomainError): replace(r, allocation_units=endless())

    def test_mapping_round_trip_snapshots_inputs_and_rejects_forged_allocations(self):
        r = self.valid(case()); raw = asdict(r)
        restored = type(r)(**raw)
        self.assertEqual(restored, r)
        raw['inputs']['stages'][0]['specific_impulse_s'] = 1.
        self.assertEqual(restored.inputs.stages[0].specific_impulse_s, 300.)
        with self.assertRaises(DomainError): replace(r, allocation_units=(39,))

    def test_maximum_grid_does_not_report_zero_as_refinement_evidence(self):
        r = self.valid(case(resolution=800))
        self.assertFalse(r.refinement_performed)
        self.assertIsNone(r.convergence_delta_initial_mass_relative)
        self.assertFalse(r.continuous_optimum_qualified)

    def test_non_nested_grid_difference_is_only_a_diagnostic(self):
        x = case(stages=(StageDefinition('a', 300., .1), StageDefinition('b', 450., .08)),
                 target_delta_v_m_s=9000., resolution=401)
        r = self.valid(x)
        coarse = reference(x, r.allocation_units)[0]
        refined = reference(dict(x, resolution=800), r.refined_allocation_units)[0]
        with localcontext() as c:
            c.prec = 1400
            expected = float(abs(coarse-refined)/refined)
        self.assertEqual(r.convergence_delta_initial_mass_relative, expected)
        self.assertFalse(r.refinement_grid_nested)
        self.assertFalse(r.continuous_optimum_qualified)

    def test_api_refuses_nontext_stage_names_instead_of_stringifying_them(self):
        payload = {'target_delta_v_m_s': 100., 'payload_mass_kg': 1., 'resolution': 40,
                       'stages': [{'name': None, 'specific_impulse_s': 300., 'structural_fraction': .1}]}
        for label, invoke in (('registered dispatch', lambda: dispatch_calculation('/api/v1/staging', payload)),
                              ('direct route boundary', lambda: api_routes.calculate_staging_route(payload))):
            with self.subTest(path=label), self.assertRaises(DomainError): invoke()


if __name__ == '__main__': unittest.main()
