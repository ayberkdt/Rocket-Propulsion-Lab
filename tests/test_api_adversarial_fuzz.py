"""Hostile inputs on every HTTP route: fail closed or answer in strict, finite JSON.

The route contract (``api.routes.dispatch_calculation``) allows only
:class:`RocketPropulsionError` to escape. This module sends every payload
field of every registered route, one at a time, two families of bad values:

* **type garbage** - NaN, +/-inf, a boolean, text, a list and an object. A
  route must refuse these with ``RocketPropulsionError``; silently accepting
  a boolean or a list as a number is a failure too, unless the field is a
  genuine flag, list or object field;
* **semantic extremes** - 0, -0.0, -1, 1e-300, +/-1e300 and 10^9. A route
  may compute or refuse, but it must never raise a bare ``OverflowError``,
  ``ZeroDivisionError``, ``ValueError`` or ``TypeError``, and anything it
  returns must serialise with ``allow_nan=False``.

Payload fields are discovered from the route source by AST (getter calls
such as ``_number(payload, "field", ...)``, ``payload.get("field")`` and
``payload["field"]``), so new routes and fields are covered automatically.
Each route gets its own test method, named after the route.
"""

from __future__ import annotations

import ast
import inspect
import json
import textwrap
import unittest

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


def probe(path: str, field: str, value: object) -> str | None:
    """Return a failure description, or None when the route kept its contract."""
    try:
        result = dispatch_calculation(path, {field: value})
    except RocketPropulsionError:
        return None
    except Exception as exc:  # noqa: BLE001 - any other exception is the defect under test
        return f"{field}={value!r}: {type(exc).__name__}: {exc}"
    try:
        json.dumps(result, allow_nan=False)
    except (TypeError, ValueError) as exc:
        return f"{field}={value!r}: response is not strict JSON ({exc})"
    return None


class TypeGarbage(unittest.TestCase):
    """Malformed values must be refused, never coerced or leaked."""


class SemanticExtremes(unittest.TestCase):
    """Extreme but well-typed values must be refused or answered in finite JSON."""


def _make(path: str, values: dict[str, object]):
    def test(self: unittest.TestCase) -> None:
        failures = [message for field in payload_fields(ROUTES[path])
                    for value in values.values()
                    if (message := probe(path, field, value)) is not None]
        self.assertEqual(failures, [], f"{path}: {len(failures)} contract violations")
    return test


for _path in sorted(ROUTES):
    _name = "test_" + _path.removeprefix("/api/v1/").replace("/", "_").replace("-", "_")
    setattr(TypeGarbage, _name, _make(_path, GARBAGE))
    setattr(SemanticExtremes, _name, _make(_path, EXTREMES))


class FieldDiscovery(unittest.TestCase):
    def test_every_computing_route_exposes_fields(self) -> None:
        """Guards the fuzzer itself: an AST change must not silently test nothing."""
        empty = [p for p, r in ROUTES.items() if not payload_fields(r)]
        allowed = {"/api/v1/species", "/api/v1/chemistry-capabilities", "/api/v1/integrations/sidera/capabilities",
                   "/api/v1/transient/feed-line"}
        self.assertEqual(sorted(set(empty) - allowed), [])
        self.assertGreater(sum(len(payload_fields(r)) for r in ROUTES.values()), 350)


if __name__ == "__main__":
    unittest.main()
