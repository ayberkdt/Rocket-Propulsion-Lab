"""Independent conditional gas/particle ledgers and hostile retained records."""
import json
import unittest
from collections.abc import Mapping
from dataclasses import asdict, replace
from decimal import Decimal, localcontext
from fractions import Fraction
from math import inf, nextafter
from sys import float_info

from rocket_propulsion.core.errors import DomainError, ConvergenceError
from rocket_propulsion.propulsion import calculate_two_phase_band
from rocket_propulsion.api.routes import dispatch_calculation, calculation_metadata

VELOCITY = 'exit-velocity-ratio'
ENERGY = 'specific-kinetic-energy-ratio'


def inputs(**changes):
    return dict(dict(ideal_thrust_n=100000., ideal_specific_impulse_s=450., condensed_mass_fraction=.2,
                     minimum_particle_coupling_efficiency=.6, maximum_particle_coupling_efficiency=.9), **changes)


def reference(x, precision=1400):
    """Decimal momentum/pressure sums; no production arithmetic helpers."""
    with localcontext() as c:
        c.prec = precision
        d = lambda v: Decimal.from_float(float(v))
        force, isp, chi = [d(x[k]) for k in ('ideal_thrust_n', 'ideal_specific_impulse_s', 'condensed_mass_fraction')]
        pressure = d(x.get('reference_pressure_thrust_n', 0.))
        momentum = force-pressure
        velocity, mixture, forces, impulses, ratios, signs = [], [], [], [], [], []
        for name in ('minimum_particle_coupling_efficiency', 'maximum_particle_coupling_efficiency'):
            speed = d(x[name])
            if x.get('coupling_definition', VELOCITY) == ENERGY:
                speed = speed.sqrt()
            eta = 1-chi+chi*speed
            thrust = momentum*eta+pressure
            velocity.append(float(speed)); mixture.append(float(eta)); forces.append(float(thrust))
            ratios.append(float(thrust/force)); impulses.append(float(isp*thrust/force)); signs.append(thrust <= 0)
        return dict(particle_velocity_ratio_range=tuple(velocity), momentum_efficiency_range=tuple(mixture),
                    performance_efficiency_range=tuple(ratios), total_thrust_ratio_range=tuple(ratios),
                    thrust_range_n=tuple(forces), specific_impulse_range_s=tuple(impulses),
                    contains_nonpositive_thrust=any(signs))


class TwoPhaseAttacks(unittest.TestCase):
    def valid(self, x):
        try:
            return calculate_two_phase_band(**x)
        except (DomainError, ConvergenceError, ArithmeticError) as error:
            self.fail(f'Representable conditional mixture witness rejected: {error}')

    def check(self, x):
        a, b = reference(x), reference(x, 2200)
        self.assertEqual(a, b, 'Independent reference precision did not settle.')
        result = self.valid(x)
        for name, expected in a.items():
            self.assertEqual(getattr(result, name, None), expected, name)
        return result

    def test_all_quantities_reject_bool_text_nonfinite_overflow_and_erased_inputs(self):
        for field in (*inputs(), 'reference_pressure_thrust_n'):
            for bad in (True, '1', None, [], float('nan'), inf, 10**10000, Fraction(1, 10**1000)):
                with self.subTest(field=field, type=type(bad).__name__):
                    with self.assertRaises(DomainError):
                        calculate_two_phase_band(**inputs(**{field: bad}))

    def test_reference_and_fraction_domains_are_strict(self):
        for field, bad in (('ideal_thrust_n', 0.), ('ideal_specific_impulse_s', -1.),
                           ('condensed_mass_fraction', nextafter(1., inf)), ('condensed_mass_fraction', -5e-324),
                           ('minimum_particle_coupling_efficiency', -5e-324),
                           ('maximum_particle_coupling_efficiency', nextafter(1., inf))):
            with self.subTest(field=field, bad=bad):
                with self.assertRaises(DomainError):
                    calculate_two_phase_band(**inputs(**{field: bad}))
        with self.assertRaises(DomainError):
            calculate_two_phase_band(**inputs(minimum_particle_coupling_efficiency=.9,
                                              maximum_particle_coupling_efficiency=.6))

    def test_reference_momentum_must_be_positive_before_output_rounding(self):
        for pressure in (100000., nextafter(100000., inf)):
            with self.subTest(pressure=pressure):
                with self.assertRaises(DomainError):
                    calculate_two_phase_band(**inputs(reference_pressure_thrust_n=pressure))

    def test_positive_subnormal_particle_velocity_and_force_are_not_erased(self):
        x = inputs(ideal_thrust_n=float_info.max, ideal_specific_impulse_s=1., condensed_mass_fraction=1.,
                   minimum_particle_coupling_efficiency=5e-324, maximum_particle_coupling_efficiency=5e-324)
        result = self.valid(x)
        expected = float(Fraction(float_info.max)*Fraction(5e-324))
        self.assertEqual(result.performance_efficiency_range, (5e-324, 5e-324))
        self.assertEqual(result.thrust_range_n, (expected, expected))
        self.assertEqual(result.specific_impulse_range_s, (5e-324, 5e-324))

    def test_default_velocity_ledgers_cover_zero_and_full_phase_limits(self):
        for chi in (0., .2, .5, 1.):
            for low, high in ((0., 0.), (.25, .25), (.6, .9), (0., 1.), (1., 1.)):
                with self.subTest(chi=chi, range=(low, high)):
                    self.check(inputs(condensed_mass_fraction=chi,
                        minimum_particle_coupling_efficiency=low, maximum_particle_coupling_efficiency=high))

    def test_pressure_thrust_is_not_multiplied_by_particle_momentum_loss(self):
        result = self.check(inputs(ideal_thrust_n=1100., reference_pressure_thrust_n=100.,
            condensed_mass_fraction=.25, minimum_particle_coupling_efficiency=.5,
            maximum_particle_coupling_efficiency=.5))
        self.assertEqual(result.thrust_range_n, (975., 975.))
        self.assertNotEqual(result.thrust_range_n[0], 1100.*.875)

    def test_specific_energy_ratio_is_converted_to_velocity_before_mixing(self):
        common = inputs(ideal_thrust_n=1000., condensed_mass_fraction=.25,
            minimum_particle_coupling_efficiency=.25, maximum_particle_coupling_efficiency=.25)
        velocity = self.check(dict(common, coupling_definition=VELOCITY))
        energy = self.check(dict(common, coupling_definition=ENERGY))
        self.assertEqual(velocity.thrust_range_n[0], 812.5)
        self.assertEqual(energy.thrust_range_n[0], 875.)
        self.assertEqual(energy.particle_velocity_ratio_range, (.5, .5))

    def test_irrational_energy_endpoints_use_independent_unrounded_roots(self):
        self.check(inputs(coupling_definition=ENERGY, reference_pressure_thrust_n=-20000.,
            minimum_particle_coupling_efficiency=.1, maximum_particle_coupling_efficiency=.7))

    def test_exact_zero_energy_cancellation_is_allowed(self):
        result = self.check(inputs(ideal_thrust_n=1., reference_pressure_thrust_n=-1.,
            condensed_mass_fraction=1., coupling_definition=ENERGY,
            minimum_particle_coupling_efficiency=.25, maximum_particle_coupling_efficiency=.25))
        self.assertEqual(result.thrust_range_n, (0., 0.))
        self.assertTrue(result.contains_nonpositive_thrust)

    def test_close_energy_cancellation_and_sign_do_not_use_displayed_velocity(self):
        boundary = -.46247529557426437
        for pressure in (nextafter(boundary, -inf), boundary, nextafter(boundary, inf)):
            with self.subTest(pressure=pressure):
                self.check(inputs(ideal_thrust_n=1., ideal_specific_impulse_s=1., reference_pressure_thrust_n=pressure,
                    condensed_mass_fraction=1., coupling_definition=ENERGY,
                    minimum_particle_coupling_efficiency=.1, maximum_particle_coupling_efficiency=.1))

    def test_conditional_negative_force_is_retained_and_warned(self):
        result = self.check(inputs(ideal_thrust_n=1., reference_pressure_thrust_n=-1.,
            condensed_mass_fraction=1., minimum_particle_coupling_efficiency=0.,
            maximum_particle_coupling_efficiency=1.))
        self.assertEqual(result.thrust_range_n, (-1., 1.))
        self.assertTrue(any(w.code == 'nonpositive_conditional_thrust' for w in result.warnings))

    def test_private_momentum_overflow_can_cancel_to_finite_public_results(self):
        for definition, low in ((VELOCITY, .5), (ENERGY, .25)):
            with self.subTest(definition=definition):
                self.check(inputs(ideal_thrust_n=float_info.max, reference_pressure_thrust_n=-float_info.max,
                    condensed_mass_fraction=.5, minimum_particle_coupling_efficiency=low,
                    maximum_particle_coupling_efficiency=low, coupling_definition=definition))

    def test_underflow_midpoint_rejected_but_next_representable_output_kept(self):
        with self.assertRaises(DomainError):
            calculate_two_phase_band(**inputs(ideal_thrust_n=5e-324, ideal_specific_impulse_s=1.,
                condensed_mass_fraction=1., minimum_particle_coupling_efficiency=.5,
                maximum_particle_coupling_efficiency=.5))
        result = self.check(inputs(ideal_thrust_n=5e-324, ideal_specific_impulse_s=1., condensed_mass_fraction=1.,
            minimum_particle_coupling_efficiency=nextafter(.5, inf), maximum_particle_coupling_efficiency=1.))
        self.assertEqual(result.thrust_range_n[0], 5e-324)

    def test_specific_energy_of_smallest_subnormal_is_not_squared_or_erased(self):
        self.check(inputs(coupling_definition=ENERGY, condensed_mass_fraction=1., ideal_thrust_n=1.,
            ideal_specific_impulse_s=1., minimum_particle_coupling_efficiency=5e-324,
            maximum_particle_coupling_efficiency=5e-324))

    def test_unrepresentable_mandatory_ratio_is_typed(self):
        with self.assertRaises(DomainError):
            calculate_two_phase_band(**inputs(ideal_thrust_n=5e-324, reference_pressure_thrust_n=-float_info.max,
                condensed_mass_fraction=1., minimum_particle_coupling_efficiency=0.,
                maximum_particle_coupling_efficiency=0.))

    def test_perfect_coupling_keeps_maximum_outputs_despite_private_overflow(self):
        for definition in (VELOCITY, ENERGY):
            with self.subTest(definition=definition):
                self.check(inputs(ideal_thrust_n=float_info.max, ideal_specific_impulse_s=float_info.max,
                    reference_pressure_thrust_n=-float_info.max, condensed_mass_fraction=1.,
                    minimum_particle_coupling_efficiency=1., maximum_particle_coupling_efficiency=1.,
                    coupling_definition=definition))

    def test_isp_does_not_reuse_the_rounded_exported_force(self):
        result = self.check(inputs(ideal_thrust_n=1., ideal_specific_impulse_s=7.))
        wrong = float(Fraction(result.thrust_range_n[0])*7)
        self.assertEqual(result.specific_impulse_range_s[0], 6.4399999999999995)
        self.assertNotEqual(result.specific_impulse_range_s[0], wrong)

    def test_fixed_pairs_and_warning_collections_stop_at_schema_cardinality(self):
        result = self.valid(inputs())
        def pair():
            yield 1.; yield 2.; yield 3.
            raise AssertionError('Pair constructor read beyond its rejection boundary.')
        for field in ('thrust_range_n', 'particle_velocity_ratio_range', 'specific_impulse_range_s'):
            with self.subTest(field=field):
                with self.assertRaises(DomainError):
                    replace(result, **{field: pair()})
        def warnings():
            yield from (result.warnings[0] for _ in range(5))
            raise AssertionError('Warning constructor read beyond its rejection boundary.')
        with self.assertRaises(DomainError):
            replace(result, warnings=warnings())

    def test_unknown_or_nontext_coupling_definitions_are_rejected(self):
        for definition in (True, None, 1, 'total-energy-ratio', 'unknown'):
            with self.subTest(definition=definition):
                with self.assertRaises(DomainError):
                    calculate_two_phase_band(**inputs(coupling_definition=definition))

    def test_result_replay_rejects_changed_original_fields_and_all_ranges(self):
        result = self.valid(inputs())
        with self.assertRaises(DomainError):
            replace(result, condensed_mass_fraction=.3)
        for field in ('thrust_range_n', 'specific_impulse_range_s', 'performance_efficiency_range',
                      'particle_coupling_efficiency_range'):
            value = getattr(result, field)
            with self.subTest(field=field):
                with self.assertRaises(DomainError):
                    replace(result, **{field: (nextafter(value[0], inf), value[1])})
        self.assertTrue(hasattr(result, 'inputs'))
        for field in ('particle_velocity_ratio_range', 'momentum_efficiency_range', 'total_thrust_ratio_range'):
            value = getattr(result, field)
            with self.subTest(field=field):
                with self.assertRaises(DomainError):
                    replace(result, **{field: (nextafter(value[0], inf), value[1])})

    def test_retained_inputs_policies_flags_and_scope_cannot_be_forged(self):
        result = self.valid(inputs())
        self.assertTrue(hasattr(result, 'inputs'))
        for changes in ({'inputs': replace(result.inputs, condensed_mass_fraction=.3)},
                        {'numerical_policy': 'unknown'}, {'band_kind': 'confidence_interval'},
                        {'physical_bound_qualified': True}, {'statistical_coverage_qualified': 0},
                        {'contains_nonpositive_thrust': 0}, {'model': 'certified hardware'},
                        {'assumptions': ('validated gas-particle feedback',)}, {'warnings': ()}):
            with self.subTest(fields=list(changes)):
                with self.assertRaises(DomainError):
                    replace(result, **changes)

    def test_json_reopening_owns_nested_mappings_and_ordered_pairs(self):
        result = self.valid(inputs())
        document = json.loads(json.dumps(asdict(result), allow_nan=False))
        reopened = type(result)(**document)
        self.assertEqual(reopened, result)
        self.assertTrue(hasattr(reopened, 'inputs'))
        document['inputs']['condensed_mass_fraction'] = .9
        document['warnings'][0]['message'] = 'changed'
        document['thrust_range_n'][0] = 0.
        self.assertEqual(reopened, result)
        with self.assertRaises(DomainError):
            type(result)(**document)

    def test_warning_text_is_normalized_before_comparison(self):
        class EqualText(str):
            def __eq__(self, other): return True
            __hash__ = str.__hash__
        result = self.valid(inputs())
        for field in ('code', 'message', 'severity', 'field'):
            warning = replace(result.warnings[0], **{field: EqualText('fabricated certification')})
            for item in (warning, asdict(warning)):
                with self.subTest(field=field, mapping=isinstance(item, dict)):
                    with self.assertRaises(DomainError):
                        replace(result, warnings=(item, *result.warnings[1:]))

    def test_mapping_protocol_errors_are_typed_without_hiding_programmer_errors(self):
        class Broken(Mapping):
            def __init__(self, error): self.error = error
            def __iter__(self): return iter(('ideal_thrust_n',))
            def __len__(self): return 1
            def __getitem__(self, key): raise self.error('broken mapping')
        result = self.valid(inputs())
        self.assertTrue(hasattr(result, 'inputs'))
        with self.assertRaises(DomainError):
            replace(result, inputs=Broken(KeyError))
        with self.assertRaises(RuntimeError):
            replace(result, inputs=Broken(RuntimeError))

    def test_scenario_ranges_make_no_physical_or_statistical_qualification_claim(self):
        result = self.valid(inputs())
        self.assertIs(getattr(result, 'physical_bound_qualified', None), False)
        self.assertIs(getattr(result, 'statistical_coverage_qualified', None), False)
        self.assertEqual(getattr(result, 'band_kind', None), 'conditional_exit_state_scenario_range')
        self.assertEqual(getattr(result, 'coupling_definition', None), VELOCITY)
        codes = {w.code for w in result.warnings}
        self.assertIn('two_phase_band_only', codes)
        self.assertIn('zero_pressure_term_assumed', codes)
        all_particle = self.valid(inputs(condensed_mass_fraction=1.))
        self.assertTrue(any(w.code == 'no_gas_phase_limit' for w in all_particle.warnings))


class TwoPhaseApiAttacks(unittest.TestCase):
    def test_dispatch_retains_pressure_force_and_specific_energy_definition(self):
        x = inputs(ideal_thrust_n=1100., reference_pressure_thrust_n=100., condensed_mass_fraction=.25,
            minimum_particle_coupling_efficiency=.25, maximum_particle_coupling_efficiency=.25,
            coupling_definition=ENERGY)
        result = dispatch_calculation('/api/v1/two-phase-band', x)
        self.assertEqual(result['thrust_range_n'], (975., 975.))
        self.assertEqual(result.get('coupling_definition'), ENERGY)
        self.assertEqual(calculation_metadata('/api/v1/two-phase-band')['version'], '2.0')
        json.dumps(result, allow_nan=False)

    def test_dispatch_does_not_stringify_invalid_definition(self):
        for value in (True, None, 1):
            with self.subTest(value=value):
                with self.assertRaises(DomainError):
                    dispatch_calculation('/api/v1/two-phase-band', inputs(coupling_definition=value))
