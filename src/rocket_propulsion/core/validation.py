"""Strict input coercion for public numerical boundaries.

These helpers validate external values at the constructors and evaluations
that call them. Their numerical coercion provides two properties:

* malformed input (``None``, text, ``bool``, containers, NaN, infinity) fails
  closed with :class:`DomainError` naming the field, never with a bare
  ``TypeError`` from deep inside an equation; and
* accepted real quantities are normalized to ``float`` without erasing nonzero
  inputs; integer fields remain exact integers. Hashed exchange artifacts are
  canonical JSON, where ``1`` and ``1.0`` serialize differently; normalizing at
  construction makes an object built with integers hash exactly like the same
  object reopened from JSON, so it can be written and read back.

``bool`` is rejected on purpose: ``True`` is an ``int`` in Python, and a
boolean passed where a physical quantity is expected is always a caller bug.
Ordered containers preserve component order. Known malformed-value protocol
failures carry a field label; unrelated application callback failures propagate.
Programmer misuse of helper configuration (record kinds or dataclass classes)
raises TypeError. Variable-size Python iterables must be finite; fixed vectors
and pairs instead perform bounded cardinality checks.
"""

from __future__ import annotations

import json
from collections.abc import Callable, Mapping
from collections.abc import Set as AbstractSet
from dataclasses import fields, is_dataclass
from itertools import islice
from math import isfinite
from numbers import Integral, Real
from types import UnionType
from typing import Union, get_args, get_origin

from .errors import DomainError, InputError


def real(value: object, label: str) -> float:
    """Normalize a finite Real, rejecting erased nonzero values and booleans."""

    if isinstance(value, bool) or not isinstance(value, Real):
        raise DomainError(f"{label} must be a real number.", field=label)
    try:
        number = float(value)
        erased_nonzero = number == 0.0 and value != 0
    except DomainError:
        raise
    except (TypeError, ValueError, OverflowError) as error:
        raise DomainError(f"{label} cannot be normalized to a finite float.", field=label) from error
    if not isfinite(number) or erased_nonzero:
        raise DomainError(f"{label} is outside the nonzero finite float range.", field=label)
    return number


def integer(value: object, label: str) -> int:
    """Return ``value`` as an ``int`` (``bool`` and floats are rejected)."""

    if isinstance(value, bool) or not isinstance(value, Integral):
        raise DomainError(f"{label} must be an integer.", field=label)
    try:
        return int(value)
    except DomainError:
        raise
    except (TypeError, ValueError, OverflowError) as error:
        raise DomainError(f"{label} cannot be normalized to an integer.", field=label) from error


def _text_value(value, label, *, allow_blank=False):
    if not isinstance(value, str) or (not allow_blank and not str.strip(value)):
        kind = 'string' if allow_blank else 'non-empty string'
        raise DomainError(f"{label} must be a {kind}.", field=label)
    try:
        str.encode(value, 'utf-8')
    except UnicodeError as error:
        raise DomainError(f"{label} must be UTF-8 encodable.", field=label) from error
    return str.__str__(value)


def text(value: object, label: str) -> str:
    """Return plain UTF-8 text without changing its literal character content."""

    return _text_value(value, label)


def _ordered_items(value, label, *, size=None):
    """Snapshot an ordered finite iterable; fixed cardinalities use bounded reads."""
    if isinstance(value, (str, bytes, bytearray, memoryview, Mapping, AbstractSet)):
        raise DomainError(f"{label} must be an ordered, non-binary iterable.", field=label)
    try:
        iterator = iter(value)
        items = tuple(iterator) if size is None else tuple(islice(iterator, size+1))
    except DomainError:
        raise
    except (TypeError, ValueError, KeyError, IndexError, OverflowError) as error:
        raise DomainError(f"{label} has an invalid iteration protocol.", field=label) from error
    if size is not None and len(items) != size:
        raise DomainError(f"{label} requires exactly {size} entries.", field=label)
    return items


def vector3(value: object, label: str) -> tuple[float, float, float]:
    """Return a three-component finite ``float`` tuple."""

    items = _ordered_items(value, label, size=3)
    return (real(items[0], label), real(items[1], label), real(items[2], label))


def real_tuple(value: object, label: str) -> tuple[float, ...]:
    """Return a tuple of finite floats from any non-string sequence."""

    return tuple(real(item, label) for item in _ordered_items(value, label))


def real_mapping(value: object, label: str) -> dict[str, float]:
    """Return a ``{non-blank str: finite float}`` copy of a mapping (``None`` -> empty)."""

    if value is None:
        return {}
    if not isinstance(value, Mapping):
        raise DomainError(f"{label} must be a mapping of names to real numbers.", field=label)
    result: dict[str, float] = {}
    try:
        for key, item in value.items():
            name = text(key, f"{label} key")
            if name in result:
                raise DomainError(f"{label} exposes duplicate mapping keys.", field=label)
            result[name] = real(item, f"{label}[{name!r}]")
    except DomainError:
        raise
    except (TypeError, ValueError, KeyError, IndexError, OverflowError) as error:
        raise DomainError(f"{label} has an invalid mapping protocol.", field=label) from error
    return result


def records(value: object, kind: type | tuple[type, ...], label: str) -> tuple:
    """Return a tuple of records of ``kind`` from any non-string iterable."""

    kinds = kind if isinstance(kind, tuple) else (kind,)
    if not kinds or any(not isinstance(item, type) for item in kinds):
        raise TypeError('records requires a type or a nonempty tuple of types.')
    items = _ordered_items(value, label)
    if any(not isinstance(item, kind) for item in items):
        raise DomainError(f"{label} contains a record of the wrong type.", field=label)
    return items


def strings(value: object, label: str, *, allow_blank: bool = False) -> tuple[str, ...]:
    """Return a tuple of strings (non-blank unless ``allow_blank``)."""

    if not isinstance(allow_blank, bool):
        raise DomainError(f"{label} allow_blank must be a boolean.", field=label)
    items = _ordered_items(value, label)
    return tuple(_text_value(item, label, allow_blank=allow_blank) for item in items)


def name_value_pairs(value: object, label: str) -> tuple[tuple[str, float], ...]:
    """Return ``((name, float), ...)`` from a sequence of two-item pairs."""

    result: list[tuple[str, float]] = []
    for item in _ordered_items(value, label):
        name, number = _ordered_items(item, label, size=2)
        result.append((text(name, f"{label} name"), real(number, label)))
    return tuple(result)


def real_matrix(value: object, label: str) -> tuple[tuple[float, ...], ...]:
    """Return a rectangular finite float matrix; empty dimensions remain valid."""

    rows = tuple(real_tuple(row, f'{label}[{index}]')
                 for index, row in enumerate(_ordered_items(value, label)))
    if rows and any(len(row) != len(rows[0]) for row in rows):
        raise DomainError(f"{label} must be rectangular.", field=label)
    return rows


def _scalar_kind(annotation):
    """Recognize supported scalar spellings without evaluating annotation code."""
    if annotation is float:
        return 'float'
    if annotation is int:
        return 'int'
    if annotation is str:
        return 'str'
    if isinstance(annotation, str):
        name = annotation.replace(' ', '')
        while len(name) >= 2 and name[0] in ('\"', "'") and name[-1] == name[0]:
            name = name[1:-1]
        name = name.replace('typing.', '').replace('builtins.', '')
        if name in ('float', 'int', 'str'):
            return name
        if name in ('float|None', 'None|float', 'Optional[float]',
                    'Union[float,None]', 'Union[None,float]',
                    'Union[float,NoneType]', 'Union[NoneType,float]'):
            return 'optional_float'
    elif get_origin(annotation) in (Union, UnionType):
        args = get_args(annotation)
        if len(args) == 2 and any(item is float for item in args) and any(item is type(None) for item in args):
            return 'optional_float'
    return None


def coerce_fields(instance: object, label: str) -> None:
    """Validate and normalize a frozen dataclass's annotated scalar fields in place.

    Fields annotated ``float`` become finite floats, ``float | None`` likewise
    unless ``None``, ``int`` become integers and ``str`` must be non-blank
    strings.  Other annotations are left to the class's own checks.  Call it
    after supported fields are initialized and before model-domain checks. Recognized
    validation completes before field writes; unrelated annotations are not
    evaluated or accessed. Optional[float] and quoted scalar spellings are
    recognized as well as the corresponding actual types.
    """

    if not is_dataclass(instance) or isinstance(instance, type):
        raise TypeError("coerce_fields requires a dataclass instance.")
    prepared = []
    for item in fields(instance):
        annotation = _scalar_kind(item.type)
        if annotation is None:
            continue
        value = getattr(instance, item.name)
        name = f"{label} {item.name}"
        if annotation == "float":
            value = real(value, name)
        elif annotation == 'optional_float':
            value = None if value is None else real(value, name)
        elif annotation == "int":
            value = integer(value, name)
        elif annotation == "str":
            value = text(value, name)
        prepared.append((item.name, value))
    for name, value in prepared:
        object.__setattr__(instance, name, value)


def json_text(value: object, label: str) -> str:
    """Return ``value`` if it is a string suitable for ``json.loads``."""

    if not isinstance(value, str):
        raise InputError(f"{label} must be JSON text.", field=label)
    return str.__str__(value)


# Every way a decoded document can fail to build a model: missing keys, wrong
# types, bad values (DomainError is a ValueError), short lists, integers beyond
# the double range, and pathological nesting.
DOCUMENT_ERRORS = (
    KeyError, TypeError, ValueError, IndexError, AttributeError, OverflowError, RecursionError,
)


def finite_float_text(value: object, label: str) -> float:
    """Parse numeric text without inventing infinity or erasing nonzero input."""
    if not isinstance(value, str):
        raise InputError(f'{label} must be numeric text.', field=label)
    literal = str.__str__(value)
    try:
        result = float(literal)
    except (ValueError, OverflowError) as error:
        raise InputError(f'{label} must contain a finite number.', field=label) from error
    erased_nonzero = result == 0.0 and any(
        c.isdecimal() and int(c) != 0 for c in literal.lower().split('e', 1)[0])
    if not isfinite(result) or erased_nonzero:
        raise InputError(f'{label} is outside the nonzero finite binary range.', field=label)
    return result


def _document_float(token):
    return finite_float_text(token, 'JSON number')


def _unique_document_object(pairs):
    """Reject duplicate names before last-wins information loss."""
    result = {}
    for key, value in pairs:
        if key in result:
            raise ValueError('Duplicate JSON member name.')
        result[key] = value
    return result


def _reject_document_constant(token):
    raise ValueError('Non-finite JSON constants are unsupported.')


def _check_document_values(payload):
    """Check decoded text and overflowed numbers, including unused fields."""
    pending = [payload]
    while pending:
        item = pending.pop()
        if isinstance(item, str):
            item.encode('utf-8')
        elif isinstance(item, float) and not isfinite(item):
            raise ValueError('JSON number is outside the finite binary range.')
        elif isinstance(item, dict):
            pending.extend(item)
            pending.extend(item.values())
        elif isinstance(item, list):
            pending.extend(item)


def read_json_document(value: object, label: str, build: Callable[[object], object]) -> object:
    """Parse supported JSON and type recognized document failures.

    ``value`` must be ``str`` (bytes are refused so the decoding is explicit).
    The root must be a JSON object.  Errors raised while building, including
    those raised by an integrity-hash check, are reported as
    :class:`InputError` (a :class:`DomainError` from a validator passes
    through unchanged, it already names the bad field).
    Duplicate member names, invalid Unicode, non-finite numbers and nonzero
    decimal literals that would round to zero are rejected before the builder.
    These are repository interoperability rules, not authentication or physical
    validation. Other readers that use different parsers require separate review.

    References
    ----------
    RFC 8259, The JavaScript Object Notation (JSON) Data Interchange Format:
    https://www.rfc-editor.org/rfc/rfc8259
    """

    source = json_text(value, label)
    try:
        source.encode('utf-8')
        payload = json.loads(source, object_pairs_hook=_unique_document_object,
                             parse_constant=_reject_document_constant, parse_float=_document_float)
        _check_document_values(payload)
    except InputError:
        raise
    except (ValueError, RecursionError) as error:  # malformed, oversized or too deep
        raise InputError(f"{label} violates the supported JSON contract.", field=label) from error
    if not isinstance(payload, dict):
        raise InputError(f"{label} root must be a JSON object.", field=label)
    try:
        return build(payload)
    except (InputError, DomainError):
        raise
    except DOCUMENT_ERRORS as error:
        raise InputError(f"{label} is malformed ({type(error).__name__}).", field=label) from error
