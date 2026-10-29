"""Numeric request admission must reject the supplied field before model work.

Stub bodies isolate this dispatcher contract; they are not valid-physics
fixtures. Real route integration tests remain necessary. Source discovery
covers scalar top-level readers and API helpers, not nested document schemas.
"""

import unittest
from unittest.mock import Mock, patch

from rocket_propulsion.api import routes as api
from rocket_propulsion.api.numeric_contracts import (
    EXACT_INTEGER_FIELDS,
    NUMBER_FIELDS,
    NUMERIC_FIELDS,
)
from rocket_propulsion.api.routes import _number
from rocket_propulsion.core.errors import InputError
from tests.api_numeric_contract_tools import numeric_payload_fields

BAD_NUMBERS = {'nan': float('nan'), 'inf': float('inf'), '-inf': -float('inf'),
               'bool': True, 'false': False, 'text': 'not-a-number', 'list': [1, 2],
               'object': {'a': 1}, 'nan_text': 'NaN', 'inf_text': 'Infinity',
               'overflow_text': '1e309', 'whitespace': ' '}


def _keyword_payload_helper(data):
    return _number(payload=data, key='gamma')


def _keyword_payload_route(payload):
    return _keyword_payload_helper(data=payload)


def _computed_key_route(payload):
    return _number(payload, 'prefix_' + 'gamma')


class NumericRequestAdmission(unittest.TestCase):
    def test_every_declared_numeric_field_is_refused_before_the_route_body(self):
        for path, fields in NUMERIC_FIELDS.items():
            body = Mock(return_value={'accepted': True})
            with patch.dict(api.ROUTES, {path: body}):
                for field in sorted(fields):
                    # 1e309 overflows binary64 but is a valid exact-integer
                    # form. Range limits, if any, belong to the real handler.
                    bad_values = {label: value for label, value in BAD_NUMBERS.items()
                                  if label != 'overflow_text' or field in NUMBER_FIELDS[path]}
                    for label, bad in bad_values.items():
                        with self.subTest(route=path, field=field, kind=label):
                            body.reset_mock()
                            with self.assertRaises(InputError) as raised:
                                api.dispatch_calculation(path, {field: bad})
                            self.assertEqual(raised.exception.field, field)
                            body.assert_not_called()

    def test_non_object_requests_are_refused_before_the_route_body(self):
        body = Mock(return_value={'accepted': True})
        with patch.dict(api.ROUTES, {'/api/v1/isentropic': body}):
            for bad in (None, [], True, 1, 'payload'):
                with self.subTest(value=bad):
                    body.reset_mock()
                    with self.assertRaises(InputError) as raised:
                        api.dispatch_calculation('/api/v1/isentropic', bad)
                    self.assertEqual(raised.exception.field, 'payload')
                    body.assert_not_called()

    def test_numeric_forms_and_optional_placeholders_are_not_rewritten(self):
        # This is type-only admission. Ranges, integrality and placeholder
        # semantics remain with the real route, not this stub.
        for path, fields in NUMBER_FIELDS.items():
            body = Mock(return_value={'accepted': True})
            with patch.dict(api.ROUTES, {path: body}):
                for field in sorted(fields):
                    for value in (1, 1.25, '1.25', -0.0, None, ''):
                        with self.subTest(route=path, field=field, value=value):
                            payload = {field: value}
                            try:
                                api.dispatch_calculation(path, payload)
                            except InputError as error:
                                self.fail(f'Type-valid numeric form refused: {error}')
                            self.assertIs(body.call_args.args[0], payload)
                            self.assertIs(payload[field], value)

    def test_exact_integer_fields_keep_large_integer_and_form_text_semantics(self):
        large = 1 << 4096  # Reproducible seeds are not restricted to binary64.
        for path, fields in EXACT_INTEGER_FIELDS.items():
            body = Mock(return_value={'accepted': True})
            with patch.dict(api.ROUTES, {path: body}):
                for field in sorted(fields):
                    for value in (large, str(large), 1.0, '1.0', '1e309'):
                        with self.subTest(route=path, field=field, type=type(value).__name__):
                            payload = {field: value}
                            try:
                                api.dispatch_calculation(path, payload)
                            except InputError as error:
                                self.fail(f'Type-valid exact integer refused: {error}')
                            self.assertIs(body.call_args.args[0], payload)
                            self.assertIs(payload[field], value)

    def test_numeric_overflow_and_fractional_integer_controls_never_run_the_body(self):
        for path, fields in NUMBER_FIELDS.items():
            body = Mock(return_value={'accepted': True})
            with patch.dict(api.ROUTES, {path: body}):
                for field in sorted(fields):
                    with self.subTest(route=path, field=field):
                        with self.assertRaises(InputError) as raised:
                            api.dispatch_calculation(path, {field: 1 << 4096})
                        self.assertEqual(raised.exception.field, field)
                        body.assert_not_called()
        for path, fields in EXACT_INTEGER_FIELDS.items():
            body = Mock(return_value={'accepted': True})
            with patch.dict(api.ROUTES, {path: body}):
                for field in sorted(fields):
                    for value in (1.25, '1.25', None, ''):
                        with self.subTest(route=path, field=field, value=value):
                            with self.assertRaises(InputError) as raised:
                                api.dispatch_calculation(path, {field: value})
                            self.assertEqual(raised.exception.field, field)
                            body.assert_not_called()

    def test_source_audit_follows_keyword_aliases_and_reports_computed_keys(self):
        numbers, integers, unresolved = numeric_payload_fields('fixture', _keyword_payload_route)
        self.assertEqual(numbers, {'gamma'})
        self.assertEqual(integers, set())
        self.assertEqual(unresolved, [])
        numbers, integers, unresolved = numeric_payload_fields('fixture', _computed_key_route)
        self.assertEqual(numbers | integers, set())
        self.assertEqual(len(unresolved), 1)
        self.assertEqual(unresolved[0][0], '_computed_key_route')

    def test_contract_catalogue_matches_scalar_readers_and_shared_helpers(self):
        self.assertEqual(set(NUMERIC_FIELDS), set(api.ROUTES))
        self.assertEqual(set(NUMBER_FIELDS), set(api.ROUTES))
        self.assertEqual(set(EXACT_INTEGER_FIELDS), set(api.ROUTES))
        for path, route in api.ROUTES.items():
            with self.subTest(route=path):
                numbers, integers, unresolved = numeric_payload_fields(path, route)
                self.assertEqual(unresolved, [])
                self.assertEqual(NUMBER_FIELDS[path], numbers)
                self.assertEqual(EXACT_INTEGER_FIELDS[path], integers)
                self.assertFalse(NUMBER_FIELDS[path] & EXACT_INTEGER_FIELDS[path])
                self.assertEqual(NUMERIC_FIELDS[path], numbers | integers)
        self.assertGreater(sum(map(len, NUMERIC_FIELDS.values())), 400)

    def test_contract_catalogue_and_field_sets_are_immutable(self):
        for catalogue in (NUMBER_FIELDS, EXACT_INTEGER_FIELDS, NUMERIC_FIELDS):
            with self.assertRaises(TypeError):
                catalogue['/api/v1/isentropic'] = frozenset()
            for fields in catalogue.values():
                self.assertIsInstance(fields, frozenset)

    def test_inactive_but_numeric_geometry_control_still_refuses_a_boolean(self):
        valid = {'contour': 'bell', 'throat_area_m2': 0.01,
                 'area_ratio': 10.0, 'characteristic_count': 12}
        self.assertTrue(api.dispatch_calculation('/api/v1/nozzle-contour', valid)['stations'])
        with self.assertRaises(InputError) as raised:
            api.dispatch_calculation('/api/v1/nozzle-contour', {**valid, 'characteristic_count': True})
        self.assertEqual(raised.exception.field, 'characteristic_count')
