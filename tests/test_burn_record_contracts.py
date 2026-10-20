"""Hostile public-record and artifact boundaries, with genuine solver controls."""

import copy
import json
import unittest
from dataclasses import fields, replace
from decimal import Decimal, localcontext
from sys import float_info

from rocket_propulsion.core.errors import DomainError, InputError
from rocket_propulsion.propulsion.burns import models as m
from rocket_propulsion.propulsion.burns import serialization as s
from rocket_propulsion.propulsion.burns.integration import simulate_constant_burn
from rocket_propulsion.propulsion.burns.profile_integration import simulate_profiled_burn
from tests.test_burn_constant_numerics import mono
from tests.test_burns import example_operating_point


class BurnRecordContractAttacks(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.definition, cls.point = mono()
        cls.constant = simulate_constant_burn(cls.definition, cls.point)
        cls.profiled = simulate_profiled_burn(cls.definition, cls.point,
            m.ThrottleSchedule((m.ThrottleSegment(2., 1., 1.),)))

    def rejects(self, build):
        with self.assertRaises(DomainError) as caught:
            build()
        self.assertIsNotNone(caught.exception.field)
        return caught.exception

    def test_boolean_and_text_physical_scalars_are_not_numbers(self):
        for record, field in ((self.point.streams[0], 'mass_flow_kg_s'),
                              (self.definition.tanks[0], 'loaded_mass_kg'),
                              (self.definition.tanks[0], 'reserve_mass_kg'),
                              (self.definition.target, 'value'),
                              (self.definition, 'maximum_duration_s'),
                              (self.definition, 'cant_efficiency'),
                              (self.point, 'commanded_throttle')):
            for bad in (True, False, '1', None, float('nan'), float('inf')):
                with self.subTest(record=type(record).__name__, field=field, bad=bad):
                    self.rejects(lambda: replace(record, **{field: bad}))

    def test_identifiers_roles_and_destinations_fail_with_field_context(self):
        for record, field in ((self.point.streams[0], 'stream_id'),
                              (self.point.streams[0], 'propellant_key'),
                              (self.point.streams[0], 'tank_id'),
                              (self.point.streams[0], 'role'),
                              (self.point.streams[0], 'destination'),
                              (self.definition.tanks[0], 'role'),
                              (self.definition.target, 'kind'),
                              (self.definition, 'name'), (self.point, 'source')):
            for bad in (None, True, 1, [], '', ' '):
                with self.subTest(record=type(record).__name__, field=field, bad=bad):
                    self.rejects(lambda: replace(record, **{field: bad}))

    def test_boolean_flags_and_fractional_engine_counts_are_rejected(self):
        for field in ('tank_depleting', 'contributes_to_delivered_thrust'):
            for bad in ('false', 'true', 0, 1, None, []):
                with self.subTest(field=field, bad=bad):
                    self.rejects(lambda: replace(self.point.streams[0], **{field: bad}))
        for bad in (True, False, 1.5, '1', None):
            with self.subTest(count=bad):
                self.rejects(lambda: replace(self.point, active_engine_count=bad))

    def test_recognized_enum_spellings_preserve_physical_dispatch_and_hash(self):
        point = example_operating_point()
        changed = replace(point, streams=tuple(replace(flow, role=flow.role.value,
                           destination=flow.destination.value) for flow in point.streams))
        self.assertEqual(s.canonical_sha256(point), s.canonical_sha256(changed))
        self.assertTrue(all(isinstance(flow.role, m.PropellantRole) for flow in changed.streams))
        self.assertTrue(all(isinstance(flow.destination, m.FlowDestination) for flow in changed.streams))
        definition = replace(self.definition, target=m.BurnTarget('duration', 1.))
        self.assertIs(definition.target.kind, m.BurnTargetKind.DURATION)
        self.assertEqual(simulate_constant_burn(definition, self.point), self.constant)

    def test_mutable_input_collections_are_snapshotted_before_hashing(self):
        streams, warnings, tanks = list(self.point.streams), ['note'], list(self.definition.tanks)
        point = replace(self.point, streams=streams, warnings=warnings)
        definition = replace(self.definition, tanks=tanks)
        before = s.burn_input_hash(definition, point)
        streams.clear()
        warnings[0] = 'changed'
        tanks.clear()
        self.assertEqual(s.burn_input_hash(definition, point), before)
        self.assertEqual(point.warnings, ('note',))
        self.assertEqual(len(point.streams), 1)
        self.assertEqual(len(definition.tanks), 1)
        self.assertEqual(simulate_constant_burn(definition, point).summary, self.constant.summary)

    def test_integer_and_float_physical_inputs_have_one_canonical_identity(self):
        def integral_copy(record):
            return replace(record, **{f.name:int(value) for f in fields(record)
                           if f.type == 'float' and isinstance(value := getattr(record, f.name), float)
                           and value.is_integer()})
        definition = integral_copy(self.definition)
        definition = replace(definition, tanks=tuple(integral_copy(t) for t in definition.tanks),
                             target=integral_copy(definition.target))
        point = integral_copy(self.point)
        point = replace(point, streams=tuple(integral_copy(flow) for flow in point.streams))
        self.assertEqual(s.burn_input_hash(definition, point), self.constant.input_hash)
        result = simulate_constant_burn(definition, point)
        self.assertEqual(result.result_hash, self.constant.result_hash)
        self.assertEqual(s.burn_result_from_json(s.burn_result_to_json(result)), result)

    def test_wrong_nested_record_types_are_typed_failures(self):
        for record, field, bad in ((self.point, 'streams', 'stream'),
                                   (self.point, 'streams', [None]),
                                   (self.point, 'warnings', 'warning'),
                                   (self.definition, 'tanks', {}),
                                   (self.definition, 'tanks', [None]),
                                   (self.definition, 'target', {})):
            with self.subTest(field=field, bad=bad):
                self.rejects(lambda: replace(record, **{field: bad}))

    def test_recoverable_isp_product_overflow_keeps_valid_full_input_scale(self):
        for precision in (100, 180):
            with localcontext() as context:
                context.prec = precision
                expected = float(Decimal.from_float(1e308) /
                    (Decimal.from_float(9.80665)*Decimal.from_float(1e308)))
            stream = replace(self.point.streams[0], mass_flow_kg_s=1e308)
            try:
                point = replace(self.point, ideal_thrust_n=1e308, delivered_thrust_n=1e308,
                                chamber_mass_flow_kg_s=1e308, total_tank_flow_kg_s=1e308,
                                system_specific_impulse_s=expected, streams=(stream,))
            except DomainError as error:
                self.fail(f'Independently finite Isp fixture was rejected: {error}')
            self.assertEqual(point.system_specific_impulse_s, expected)

    def test_unrepresentable_mixture_ratio_cannot_pass_inf_comparison(self):
        tiny = float.fromhex('0x0.0000000000001p-1022')
        ox = replace(self.point.streams[0], role=m.PropellantRole.OXIDIZER)
        fuel = replace(ox, stream_id='fuel', tank_id='fuel', role=m.PropellantRole.FUEL,
                       mass_flow_kg_s=tiny)
        error = self.rejects(lambda: replace(self.point, streams=(ox, fuel),
                                    oxidizer_mass_flow_kg_s=1., fuel_mass_flow_kg_s=tiny,
                                    mixture_ratio=float_info.max))
        self.assertEqual(error.field, 'mixture_ratio')

    def test_budget_comparison_does_not_turn_double_overflow_into_acceptance(self):
        tank = replace(self.definition.tanks[0], loaded_mass_kg=float_info.max)
        with localcontext() as context:
            context.prec = 500
            self.assertGreater(Decimal.from_float(tank.loaded_mass_kg)*2,
                Decimal.from_float(float_info.max)*Decimal.from_float(1.+1e-12))
        error = self.rejects(lambda: replace(self.definition, initial_mass_kg=float_info.max,
                                    protected_dry_mass_kg=float_info.max, tanks=(tank,)))
        self.assertEqual(error.field, 'initial_mass_kg')
        legal = replace(self.definition, initial_mass_kg=float_info.max,
                        protected_dry_mass_kg=0., tanks=(tank,))
        self.assertEqual(legal.tanks, (tank,))

    def test_untracked_paths_keep_blank_tank_and_zero_flow_controls(self):
        flow = replace(self.point.streams[0], tank_id='', mass_flow_kg_s=0.,
                       tank_depleting=False, contributes_to_delivered_thrust=False)
        self.assertEqual(flow.tank_id, '')
        self.assertFalse(flow.tank_depleting)
        self.assertEqual(flow.mass_flow_kg_s, 0.)
        self.rejects(lambda: replace(flow, tank_id=None))

    def test_state_pairs_are_strict_unique_nonnegative_and_immutable(self):
        for pairs in ((('tank', 1., 'ignored'),), (('tank', True),), (('tank', -1.),),
                      (('tank', 1.), ('tank', 2.)), (('', 1.),), (('\ud800', 1.),), None, {'tank':1.}):
            with self.subTest(pairs=pairs):
                self.rejects(lambda: m.PropulsionState(0., 10., pairs))
        pairs = [['tank', 1]]
        state = m.PropulsionState(0, 10, pairs)
        pairs[0][1] = 7.
        self.assertEqual(state.tank_masses_kg, (('tank', 1.),))
        for field, bad in (('time_s', True), ('time_s', -1.), ('vehicle_mass_kg', 0.),
                            ('vehicle_mass_kg', float('nan'))):
            with self.subTest(field=field):
                self.rejects(lambda: replace(state, **{field:bad}))

    def test_profile_point_scalars_and_phases_reject_malformed_data(self):
        for point in (self.constant.profile[0], self.profiled.profile[0]):
            for field, bad in (('time_s', True), ('phase', None), ('vehicle_mass_kg', 0.),
                                ('delivered_thrust_n', -1.), ('total_mass_flow_kg_s', float('nan'))):
                with self.subTest(point=type(point).__name__, field=field):
                    self.rejects(lambda: replace(point, **{field:bad}))
        for bad in (True, -1., 1.0001):
            self.rejects(lambda: replace(self.profiled.profile[0], commanded_throttle=bad))
        idle = replace(self.profiled.profile[0], commanded_throttle=0., realized_throttle=0.,
                       delivered_thrust_n=0., axial_thrust_n=0., total_mass_flow_kg_s=0.)
        self.assertEqual(idle.commanded_throttle, 0.)

    def test_summary_and_evidence_flags_counts_and_domains_are_checked(self):
        for record, field, bad in ((self.constant.summary,'target_achieved','false'),
                                  (self.constant.summary,'exact_constant_sidera_export_available',1),
                                  (self.constant.summary,'consumed_propellant_kg',-1.),
                                  (self.constant.summary,'duration_s',True),
                                  (self.constant.summary,'cutoff_reason','unknown'),
                                  (self.constant.evidence,'accepted_steps',1.5),
                                  (self.constant.evidence,'rejected_steps',True),
                                  (self.constant.evidence,'absolute_tolerance',-1.)):
            with self.subTest(record=type(record).__name__, field=field):
                self.rejects(lambda: replace(record, **{field:bad}))
        self.assertEqual(replace(self.constant.evidence, segment_count=0,
                                accepted_steps=0, absolute_tolerance=0.).accepted_steps, 0)
        self.assertEqual(replace(self.constant.evidence, accepted_steps=10**400).accepted_steps, 10**400)

    def test_schedule_lists_are_copied_and_malformed_commands_fail_typed(self):
        segments = [m.ThrottleSegment(1., 0., 1.)]
        schedule = m.ThrottleSchedule(segments)
        before = s.canonical_sha256(schedule)
        segments.clear()
        self.assertEqual(s.canonical_sha256(schedule), before)
        self.assertEqual(len(schedule.segments), 1)
        for bad in (None, {}, 'segments', [None], []):
            with self.subTest(segments=bad):
                self.rejects(lambda: m.ThrottleSchedule(bad))
        for bad in (-1., 1.01, True, '1'):
            with self.subTest(throttle=bad):
                self.rejects(lambda: m.ThrottleSegment(1., bad, 1.))

    def test_result_snapshots_and_top_level_record_types_are_checked(self):
        for result in (self.constant, self.profiled):
            profile, warnings = list(result.profile), list(result.warnings)
            changed = replace(result, profile=profile, warnings=warnings)
            profile.clear()
            warnings.append('injected')
            self.assertEqual(s.burn_result_hash(changed), result.result_hash)
            for field, bad in (('definition', {}), ('summary', None), ('evidence', {}),
                                ('profile', [None]), ('schema', 'unknown'), ('input_hash', 'bad')):
                with self.subTest(result=type(result).__name__, field=field):
                    self.rejects(lambda: replace(result, **{field:bad}))

    def test_json_reader_cannot_whiten_raw_scalar_types_using_old_valid_hashes(self):
        for result in (self.constant, self.profiled):
            baseline = json.loads(s.burn_result_to_json(result))
            attacks = (('summary','target_achieved','false'),
                       ('summary','exact_constant_sidera_export_available',1),
                       ('operating_point','commanded_throttle',True),
                       ('operating_point','active_engine_count',1.75),
                       ('evidence','accepted_steps',1.75),
                       ('definition','maximum_duration_s','100.0'),
                       ('operating_point','warnings',''))
            for parent, field, value in attacks:
                payload = copy.deepcopy(baseline)
                payload[parent][field] = value
                with self.subTest(schema=result.schema, parent=parent, field=field), self.assertRaises(InputError):
                    s.burn_result_from_json(json.dumps(payload))

    def test_json_pair_shape_cannot_discard_unhashed_trailing_content(self):
        for result in (self.constant, self.profiled):
            payload = json.loads(s.burn_result_to_json(result))
            for parent, field in (('final_state','tank_masses_kg'), ('summary','consumed_by_tank_kg')):
                bad = copy.deepcopy(payload)
                bad[parent][field][0].append('unhashed')
                with self.subTest(schema=result.schema, field=field), self.assertRaises(InputError):
                    s.burn_result_from_dict(bad)

    def test_unknown_fields_are_rejected_at_every_record_level(self):
        for result in (self.constant, self.profiled):
            baseline = json.loads(s.burn_result_to_json(result))
            paths = [(), ('definition',), ('definition','target'), ('definition','tanks',0),
                     ('operating_point',), ('operating_point','streams',0), ('summary',),
                     ('evidence',), ('initial_state',), ('profile',0)]
            if result is self.profiled:
                paths += [('schedule',), ('schedule','segments',0)]
            for path in paths:
                payload = copy.deepcopy(baseline)
                node = payload
                for part in path:
                    node = node[part]
                node['unhashed_extension'] = {'claim':'accepted'}
                with self.subTest(schema=result.schema, path=path), self.assertRaises(InputError):
                    s.burn_result_from_dict(payload)

    def test_duplicate_json_names_and_nonstandard_constants_are_rejected(self):
        for result in (self.constant, self.profiled):
            encoded = s.burn_result_to_json(result, indent=None)
            duplicated = encoded.replace('"target_achieved": true',
                                         '"target_achieved": false, "target_achieved": true')
            self.assertNotEqual(encoded, duplicated)
            with self.assertRaises(InputError):
                s.burn_result_from_json(duplicated)
            # Even discarded extension fields must not admit NaN/Infinity.
            for bad in ('NaN','Infinity','-Infinity'):
                raw = encoded[:-1]+', "discarded": '+bad+'}'
                with self.subTest(constant=bad), self.assertRaises(InputError):
                    s.burn_result_from_json(raw)
                raw = encoded.replace('"target_achieved": true', '"target_achieved": '+bad)
                with self.subTest(known_flag_constant=bad), self.assertRaises(InputError):
                    s.burn_result_from_json(raw)

    def test_direct_dict_reader_returns_typed_errors_for_malformed_nesting(self):
        baseline = json.loads(s.burn_result_to_json(self.constant))
        for field, value in (('schema', []), ('definition', None), ('profile', 'profile'),
                              ('operating_point', []), ('summary', None)):
            payload = copy.deepcopy(baseline)
            payload[field] = value
            with self.subTest(field=field), self.assertRaises(InputError):
                s.burn_result_from_dict(payload)

    def test_writer_rejects_provisional_stale_and_self_consistent_forged_results(self):
        for result in (self.constant, self.profiled):
            stale = replace(result, result_hash='')
            altered = replace(result, summary=replace(result.summary,
                              ideal_delta_v_m_s=result.summary.ideal_delta_v_m_s+1.))
            forged = replace(altered, result_hash=s.burn_result_hash(altered))
            for bad in (stale, altered, forged, {}, None):
                with self.subTest(schema=result.schema, bad=type(bad).__name__), self.assertRaises(InputError):
                    s.burn_result_to_dict(bad)

    def test_ordinary_v1_v2_artifacts_and_result_pair_snapshots_still_roundtrip(self):
        for result in (self.constant, self.profiled):
            for field in ('initial_state','final_state'):
                original = getattr(result, field)
                pairs = [list(item) for item in original.tank_masses_kg]
                changed = replace(original, tank_masses_kg=pairs)
                pairs[0][1] = 999.
                self.assertEqual(changed, original)
            payload = s.burn_result_to_dict(result)
            self.assertEqual(s.burn_result_from_dict(payload), result)
            self.assertEqual(s.burn_result_from_json(s.burn_result_to_json(result)), result)
            self.assertEqual(s.burn_result_hash(result), result.result_hash)
            self.assertEqual(s.burn_result_to_dict(result), payload)

    def test_unmodified_float_v1_v2_artifacts_reproduce_without_normalization_drift(self):
        for result in (self.constant, self.profiled):
            encoded = s.burn_result_to_json(result)
            self.assertEqual(s.burn_result_from_json(encoded), result)
            self.assertEqual(s.burn_result_from_dict(json.loads(encoded)), result)
            self.assertEqual(s.burn_result_to_json(result), encoded)
            self.assertEqual(s.burn_result_hash(result), result.result_hash)

    def test_unicode_scalar_text_roundtrips_but_unpaired_surrogates_fail_typed(self):
        point = replace(self.point, name='İtici ölçümü 🚀', warnings=['Δv denetimi'])
        result = simulate_constant_burn(self.definition, point)
        encoded = s.burn_result_to_json(result)
        self.assertIn('🚀', encoded)
        self.assertEqual(s.burn_result_from_json(encoded), result)
        for record, field, bad in ((self.point, 'name', '\ud800'),
                                  (self.point, 'warnings', ['\udfff']),
                                  (self.point.streams[0], 'stream_id', '\ud800'),
                                  (self.definition.tanks[0], 'tank_id', '\udfff'),
                                  (self.constant.summary, 'binding_constraint', '\ud800'),
                                  (self.profiled.profile[0], 'phase', '\udfff')):
            with self.subTest(field=field):
                self.rejects(lambda: replace(record, **{field:bad}))
        # Non-depleting streams permit an empty tank identifier, but its text
        # must still be encodable; this branch has no nonblank-name shortcut.
        self.rejects(lambda: replace(self.point.streams[0], tank_depleting=False,
                                     tank_id='\ud800'))
        payload = json.loads(s.burn_result_to_json(self.constant))
        payload['operating_point']['name'] = '\ud800'
        with self.assertRaises(InputError):
            s.burn_result_from_json(json.dumps(payload))
