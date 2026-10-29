"""Reject compact integer-expansion attacks before allocation or model work.

An independent Decimal.__int__ tripwire makes these attacks safe even when
the production resource guard is absent. No billion-digit integer is built.
"""

from __future__ import annotations

import json
import threading
import unittest
from decimal import Decimal
from http.server import ThreadingHTTPServer
from unittest.mock import Mock, patch
from urllib.error import HTTPError
from urllib.request import Request, urlopen

from rocket_propulsion.api import routes as api
from rocket_propulsion.core.errors import InputError
from tests.test_http_server import QuietHandler

FIELDS = ('seed', 'sample_count', 'sobol_base_samples')
API_PATH = '/api/v1/uncertainty'


def tripwire_decimal(attempts):
    """Instrument conversion, allowing small controls but stopping huge allocation."""
    class GuardedDecimal(Decimal):
        def __int__(self):
            digits = 1 if self.is_zero() else self.adjusted() + 1
            attempts.append(digits)
            if digits > 10000:
                raise AssertionError('Unsafe integer expansion reached the allocation boundary.')
            return super().__int__()
    return GuardedDecimal


class IntegerExpansionBudgetTests(unittest.TestCase):
    def test_compact_exponents_refuse_before_expansion_and_model(self):
        for field in FIELDS:
            for text in ('1e1000000000', '-1e1000000000', '1e4300', '-1e4300'):
                attempts, body = [], Mock(return_value={'accepted': True})
                with (self.subTest(field=field, text=text),
                      patch.object(api, 'Decimal', tripwire_decimal(attempts)),
                      patch.dict(api.ROUTES, {API_PATH: body})):
                    with self.assertRaises(InputError) as raised:
                        api.dispatch_calculation(API_PATH, {field: text})
                    self.assertEqual(raised.exception.field, field)
                    self.assertEqual(raised.exception.code, 'integer_budget_exceeded')
                    self.assertEqual(raised.exception.details['maximum_decimal_digits'], 4300)
                    self.assertEqual(attempts, [])
                    body.assert_not_called()

    def test_allocated_oversized_integers_refuse_without_stringification(self):
        class NoStringInteger(int):
            def __str__(self):
                raise AssertionError('Resource guard tried decimal stringification.')
            __repr__ = __str__

        for field in FIELDS:
            for sign in (-1, 1):
                with self.subTest(field=field, sign=sign):
                    value = NoStringInteger(sign * 10**4300)
                    with self.assertRaises(InputError) as raised:
                        api._study_request_integer({field: value}, field, 0)
                    self.assertEqual(raised.exception.code, 'integer_budget_exceeded')
                    self.assertEqual(raised.exception.field, field)

    def test_boundary_and_supported_large_forms_stay_exact_and_json_encodable(self):
        largest = 10**4300 - 1
        large_seed = 1 << 4096
        forms = [(largest, largest), (-largest, -largest), (str(largest), largest),
                 ('1e4299', 10**4299), ('1e309', 10**309),
                 (large_seed, large_seed), (str(large_seed), large_seed),
                 (1.0, 1), ('1.0', 1)]
        for value, expected in forms:
            with self.subTest(input_type=type(value).__name__, expected_bits=expected.bit_length()):
                result = api._study_request_integer({'seed': value}, 'seed', 0)
                self.assertEqual(result, expected)
                self.assertIs(type(result), int)
                self.assertEqual(json.loads(json.dumps({'seed': result}))['seed'], expected)

    def test_lower_host_limit_is_honored_and_disabled_host_limit_keeps_api_cap(self):
        for configured, accepted, refused in [(640, '1e639', '1e640'),
                                               (0, '1e4299', '1e4300')]:
            with (self.subTest(configured=configured),
                  patch.object(api, 'get_int_max_str_digits', return_value=configured, create=True),
                  patch('sys.set_int_max_str_digits', side_effect=AssertionError('Global limit changed.'))):
                self.assertGreater(api._study_request_integer({'seed': accepted}, 'seed', 0), 0)
                with self.assertRaises(InputError) as raised:
                    api._study_request_integer({'seed': refused}, 'seed', 0)
                self.assertEqual(raised.exception.code, 'integer_budget_exceeded')
                self.assertEqual(raised.exception.details['maximum_decimal_digits'], configured or 4300)

    def test_zero_exponents_and_fractional_refusal_do_not_expand(self):
        attempts = []
        with patch.object(api, 'Decimal', tripwire_decimal(attempts)):
            for text in ('0e1000000000', '-0e1000000000', '0e-1000000000'):
                self.assertEqual(api._study_request_integer({'seed': text}, 'seed', 0), 0)
            attempts.clear()
            with self.assertRaises(InputError) as raised:
                api._study_request_integer({'seed': '1e-1000000000'}, 'seed', 0)
            self.assertEqual(raised.exception.field, 'seed')
            self.assertNotEqual(raised.exception.code, 'integer_budget_exceeded')
            self.assertEqual(attempts, [])


class IntegerBudgetHttpTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.server = ThreadingHTTPServer(('127.0.0.1', 0), QuietHandler)
        cls.thread = threading.Thread(target=cls.server.serve_forever, daemon=True)
        cls.thread.start()
        cls.url = f'http://127.0.0.1:{cls.server.server_port}{API_PATH}'

    @classmethod
    def tearDownClass(cls):
        cls.server.shutdown()
        cls.server.server_close()
        cls.thread.join(timeout=2.)

    def post(self, payload):
        request = Request(self.url, data=json.dumps(payload).encode(),
                          headers={'Content-Type': 'application/json'}, method='POST')
        try:
            with urlopen(request, timeout=5.) as response:
                return response.status, json.load(response)
        except HTTPError as error:
            return error.code, json.load(error)

    def test_compact_attack_returns_named_400_instead_of_internal_error_or_disconnect(self):
        for field in FIELDS:
            attempts = []
            with (self.subTest(field=field),
                  patch.object(api, 'Decimal', tripwire_decimal(attempts))):
                status, body = self.post({field: '1e1000000000'})
                self.assertEqual(status, 400)
                self.assertEqual(body['error']['code'], 'integer_budget_exceeded')
                self.assertEqual(body['error']['field'], field)
                self.assertEqual(attempts, [])

    def test_supported_large_seed_survives_real_uncertainty_and_http_serialization(self):
        status, body = self.post({
            'parameters':[
                {'name':'specific_impulse_s','minimum':430,'maximum':460},
                {'name':'mass_ratio','minimum':3.,'maximum':4.}],
            'sample_count':8,'sobol_base_samples':8,'seed':'1e309',
        })
        self.assertEqual(status, 200)
        self.assertEqual(body['data']['seed'], 10**309)
        self.assertEqual(len(body['data']['samples']), 8)


if __name__ == '__main__':
    unittest.main()
