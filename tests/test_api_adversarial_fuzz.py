"""Hostile inputs on declared top-level API fields.

The route contract (``api.routes.dispatch_calculation``) allows only
:class:`RocketPropulsionError` to escape on known routes. This module sends
discovered fields, one at a time, two families of values:

* **type garbage** - NaN, +/-inf, a boolean, text, a list and an object. A
  declared numeric fields must be refused with an error naming that field.
  Other fields receive a typed-error/finite-JSON smoke check; it does not
  establish that malformed flags, lists or nested objects are refused;
* **semantic extremes** - 0, -0.0, -1, 1e-300, +/-1e300 and 10^9. A route
  may compute or refuse, but it must never raise a bare ``OverflowError``,
  ``ZeroDivisionError``, ``ValueError`` or ``TypeError``, and anything it
  returns must serialise with ``allow_nan=False``.

Direct getters are discovered by AST and supplemented with the audited
scalar-helper catalogue. Nested schemas and branch-dependent semantic ranges
are not comprehensively covered. These single-field requests isolate numeric
admission; they are not baseline-valid physics fixtures. Each route gets its
own test method, named after the route.
"""

from __future__ import annotations

import ast
import inspect
import json
import textwrap
import unittest
from unittest.mock import patch

from rocket_propulsion.api.numeric_contracts import NUMERIC_FIELDS
from rocket_propulsion.api.routes import ROUTES, dispatch_calculation
from rocket_propulsion.core.errors import RocketPropulsionError

GARBAGE = {"nan": float("nan"), "inf": float("inf"), "-inf": float("-inf"), "bool": True,
           "text": "not-a-number", "list": [1, 2], "object": {"a": 1}}
EXTREMES = {"zero": 0.0, "negative_zero": -0.0, "minus_one": -1.0, "tiny": 1e-300,
            "huge": 1e300, "minus_huge": -1e300, "billion": 10 ** 9}


def payload_fields(route) -> list[str]:
    """String keys a route reads from its payload, found in its source."""
    try:
        tree = ast.parse(textwrap.dedent(inspect.getsource(route)))
    except (OSError, TypeError):
        return []
    keys: set[str] = set()
    for node in ast.walk(tree):
        if isinstance(node, ast.Call):
            args = node.args
            if (isinstance(node.func, ast.Attribute) and node.func.attr == "get" and args
                    and isinstance(args[0], ast.Constant) and isinstance(args[0].value, str)
                    and isinstance(node.func.value, ast.Name) and node.func.value.id == "payload"):
                keys.add(args[0].value)
            elif (len(args) >= 2 and isinstance(args[0], ast.Name) and args[0].id == "payload"
                  and isinstance(args[1], ast.Constant) and isinstance(args[1].value, str)):
                keys.add(args[1].value)
        elif (isinstance(node, ast.Subscript) and isinstance(node.value, ast.Name)
              and node.value.id == "payload" and isinstance(node.slice, ast.Constant)
              and isinstance(node.slice.value, str)):
            keys.add(node.slice.value)
    return sorted(keys)


def probe(path: str, field: str, value: object, *, require_refusal: bool = False) -> str | None:
    """Return a failure description, or None when the route kept its contract."""
    try:
        result = dispatch_calculation(path, {field: value})
    except RocketPropulsionError as exc:
        if require_refusal and exc.field != field:
            return f"{field}={value!r}: refusal names unrelated field {exc.field!r}"
        return None
    except Exception as exc:  # noqa: BLE001 - any other exception is the defect under test
        return f"{field}={value!r}: {type(exc).__name__}: {exc}"
    if require_refusal:
        return f"{field}={value!r}: malformed number returned success"
    try:
        json.dumps(result, allow_nan=False)
    except (TypeError, ValueError) as exc:
        return f"{field}={value!r}: response is not strict JSON ({exc})"
    return None


class TypeGarbage(unittest.TestCase):
    """Malformed values must be refused, never coerced or leaked."""


class SemanticExtremes(unittest.TestCase):
    """Extreme but well-typed values must be refused or answered in finite JSON."""


def _make(path: str, values: dict[str, object], *, numeric_refusal: bool = False):
    def test(self: unittest.TestCase) -> None:
        fields = sorted(set(payload_fields(ROUTES[path])) | NUMERIC_FIELDS[path])
        failures = [message for field in fields
                    for value in values.values()
                    if (message := probe(path, field, value,
                                         require_refusal=numeric_refusal and field in NUMERIC_FIELDS[path])) is not None]
        self.assertEqual(failures, [], f"{path}: {len(failures)} contract violations")
    return test


for _path in sorted(ROUTES):
    _name = "test_" + _path.removeprefix("/api/v1/").replace("/", "_").replace("-", "_")
    setattr(TypeGarbage, _name, _make(_path, GARBAGE, numeric_refusal=True))
    setattr(SemanticExtremes, _name, _make(_path, EXTREMES))


class FieldDiscovery(unittest.TestCase):
    def test_every_computing_route_exposes_fields(self) -> None:
        """Guards the fuzzer itself: an AST change must not silently test nothing."""
        empty = [p for p, r in ROUTES.items() if not payload_fields(r)]
        allowed = {"/api/v1/species", "/api/v1/chemistry-capabilities", "/api/v1/integrations/sidera/capabilities",
                   "/api/v1/transient/feed-line"}
        self.assertEqual(sorted(set(empty) - allowed), [])
        self.assertGreater(sum(len(payload_fields(r)) for r in ROUTES.values()), 350)


class RefusalOracle(unittest.TestCase):
    def test_finite_success_is_a_failure_when_refusal_is_required(self):
        with patch(__name__ + '.dispatch_calculation', return_value={'accepted': True}):
            self.assertIsNotNone(probe('/api/v1/isentropic', 'value', True, require_refusal=True))

    def test_an_unrelated_field_error_cannot_prove_refusal(self):
        error = RocketPropulsionError('Missing another field.', field='unrelated')
        with patch(__name__ + '.dispatch_calculation', side_effect=error):
            self.assertIsNotNone(probe('/api/v1/isentropic', 'value', True, require_refusal=True))

    def test_the_correct_field_error_proves_refusal(self):
        error = RocketPropulsionError('Malformed number.', field='value')
        with patch(__name__ + '.dispatch_calculation', side_effect=error):
            self.assertIsNone(probe('/api/v1/isentropic', 'value', True, require_refusal=True))

    def test_semantic_extremes_can_still_produce_finite_success(self):
        with patch(__name__ + '.dispatch_calculation', return_value={'accepted': True}):
            self.assertIsNone(probe('/api/v1/isentropic', 'value', 0.0))

    def test_generated_numeric_family_requires_refusal(self):
        generated = _make('/api/v1/isentropic', {'bool': True}, numeric_refusal=True)
        with (patch(__name__ + '.dispatch_calculation', return_value={'accepted': True}),
              self.assertRaises(AssertionError)):
            generated(self)

    def test_generated_family_covers_numbers_read_only_in_helpers(self):
        path = '/api/v1/transient/feed-line'
        self.assertNotIn('bottle_pressure_pa', payload_fields(ROUTES[path]))
        self.assertIn('bottle_pressure_pa', NUMERIC_FIELDS[path])
        generated = _make(path, {'bool': True}, numeric_refusal=True)
        with (patch(__name__ + '.dispatch_calculation', return_value={'accepted': True}) as call,
              self.assertRaises(AssertionError)):
            generated(self)
        self.assertTrue(any(payload == {'bottle_pressure_pa': True}
                            for _, payload in (item.args for item in call.call_args_list)))


if __name__ == "__main__":
    unittest.main()
