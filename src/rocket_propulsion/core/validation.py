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
boolean passed where a physical quantity is expected is always a caller error.
Ordered containers preserve component order. Known malformed-value protocol
failures carry a field label; unrelated application callback failures propagate.
Programmer misuse of helper configuration (record kinds or dataclass classes)
raises TypeError. Variable-size Python iterables must be finite; fixed vectors
and pairs instead perform bounded cardinality checks.

References
----------
.. [1] PYTHON-3.12-NUMBERS. Python 3.12, "numbers - Numeric abstract base
       classes", https://docs.python.org/3.12/library/numbers.html
.. [2] IETF RFC 8259, "The JavaScript Object Notation (JSON) Data Interchange
       Format", https://www.rfc-editor.org/rfc/rfc8259
.. [3] PYTHON-3.12-JSON. Python 3.12, "json - JSON encoder and decoder",
       https://docs.python.org/3.12/library/json.html
"""

from __future__ import annotations

import functools
import json
from collections.abc import Callable, Mapping
from collections.abc import Set as AbstractSet
from dataclasses import fields, is_dataclass
from itertools import islice
from math import isfinite
from numbers import Integral, Real
from types import UnionType
from typing import ParamSpec, TypeVar, Union, get_args, get_origin

from .errors import DomainError, InputError

_P = ParamSpec("_P")
_R = TypeVar("_R")


def real(value: object, label: str) -> float:
    """Normalize a finite real number to ``float``.

    Parameters
    ----------
    value : object
        Candidate value; any :class:`numbers.Real` except ``bool`` [1]_.
    label : str
        Field name used in error messages and as the error ``field``.

    Returns
    -------
    float
        ``float(value)``.

    Raises
    ------
    DomainError
        If ``value`` is a ``bool`` or not a ``Real``, cannot be converted, is
        not finite after conversion, or is a nonzero value that converts to
        ``0.0`` (underflow erasing a nonzero quantity).

    References
    ----------
    .. [1] PYTHON-3.12-NUMBERS. Python 3.12, "numbers - Numeric abstract base
           classes", https://docs.python.org/3.12/library/numbers.html
    """

    # A native float is already normalized: no conversion can erase a
    # nonzero value. Keep checking finiteness on *every* call. Subclasses
    # must take the protocol path below so their declared errors/conversion
    # behavior are preserved instead of being bypassed by this fast path.
    if type(value) is float:
        if not isfinite(value):
            raise DomainError(f"{label} is outside the nonzero finite float range.", field=label)
        return value
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
    """Return an integral value as an exact ``int``.

    Parameters
    ----------
    value : object
        Candidate value; any :class:`numbers.Integral` except ``bool`` [1]_.
        Floats, including integral-valued floats, are rejected.
    label : str
        Field name used in error messages.

    Returns
    -------
    int
        ``int(value)``.

    Raises
    ------
    DomainError
        If ``value`` is a ``bool``, not ``Integral``, or cannot be converted.

    References
    ----------
    .. [1] PYTHON-3.12-NUMBERS. Python 3.12, "numbers - Numeric abstract base
           classes", https://docs.python.org/3.12/library/numbers.html
    """

    if isinstance(value, bool) or not isinstance(value, Integral):
        raise DomainError(f"{label} must be an integer.", field=label)
    try:
        return int(value)
    except DomainError:
        raise
    except (TypeError, ValueError, OverflowError) as error:
        raise DomainError(f"{label} cannot be normalized to an integer.", field=label) from error


def _text_value(value, label, *, allow_blank=False):
    """Return ``value`` as a plain ``str`` after type, blank and UTF-8 checks.

    A ``str`` subclass is reduced to a plain ``str`` with identical
    characters. Unless ``allow_blank`` is true, whitespace-only text is
    refused. Text that cannot be encoded as UTF-8 (for example a lone
    surrogate) raises :class:`DomainError`.
    """
    if not isinstance(value, str) or (not allow_blank and not str.strip(value)):
        kind = 'string' if allow_blank else 'non-empty string'
        raise DomainError(f"{label} must be a {kind}.", field=label)
    try:
        str.encode(value, 'utf-8')
    except UnicodeError as error:
        raise DomainError(f"{label} must be UTF-8 encodable.", field=label) from error
    return str.__str__(value)


def text(value: object, label: str) -> str:
    """Return non-blank, UTF-8 encodable text as a plain ``str``.

    Parameters
    ----------
    value : object
        Candidate text. ``str`` subclasses are accepted and reduced to a
        plain ``str``; the characters are not stripped or normalized.
    label : str
        Field name used in error messages.

    Returns
    -------
    str
        The same characters as ``value``.

    Raises
    ------
    DomainError
        If ``value`` is not a ``str``, is blank, or cannot be encoded as
        UTF-8.
    """

    return _text_value(value, label)


def _ordered_items(value, label, *, size=None):
    """Snapshot an ordered iterable into a tuple.

    Strings, binary buffers, mappings and sets are refused because they are
    either scalars in disguise or unordered. With ``size`` given, at most
    ``size + 1`` items are read, so an endless iterator cannot block a
    fixed-cardinality check; without it the iterable must be finite. Known
    iteration-protocol failures become :class:`DomainError` naming ``label``.
    """
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
    """Return a three-component tuple of finite floats.

    Parameters
    ----------
    value : object
        Ordered iterable with exactly three real components.
    label : str
        Field name used in error messages.

    Returns
    -------
    tuple of float
        The components, normalized with :func:`real`.

    Raises
    ------
    DomainError
        If ``value`` is unordered, string-like, has another length, or any
        component is not a finite real.
    """

    items = _ordered_items(value, label, size=3)
    return (real(items[0], label), real(items[1], label), real(items[2], label))


def real_tuple(value: object, label: str) -> tuple[float, ...]:
    """Return a tuple of finite floats from an ordered, finite iterable.

    Parameters
    ----------
    value : object
        Ordered iterable of real numbers (strings, mappings and sets are
        refused).
    label : str
        Field name used in error messages.

    Returns
    -------
    tuple of float
        Components in their original order; may be empty.

    Raises
    ------
    DomainError
        If ``value`` is not an accepted iterable or an item is not a finite
        real.
    """

    return tuple(real(item, label) for item in _ordered_items(value, label))


def real_mapping(value: object, label: str) -> dict[str, float]:
    """Return a ``{non-blank str: finite float}`` copy of a mapping.

    Parameters
    ----------
    value : object
        A :class:`~collections.abc.Mapping` of names to real numbers, or
        ``None``.
    label : str
        Field name used in error messages.

    Returns
    -------
    dict of str to float
        A new dictionary; empty when ``value`` is ``None``.

    Raises
    ------
    DomainError
        If ``value`` is not a mapping, a key is not non-blank text, a value
        is not a finite real, two keys normalize to the same text, or the
        mapping's iteration protocol fails.
    """

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
    """Return a tuple of records that are all instances of ``kind``.

    Parameters
    ----------
    value : object
        Ordered, finite iterable of records.
    kind : type or tuple of type
        Accepted record type(s).
    label : str
        Field name used in error messages.

    Returns
    -------
    tuple
        The records in their original order.

    Raises
    ------
    TypeError
        If ``kind`` is not a type or a non-empty tuple of types (a
        programming error, not an input error).
    DomainError
        If ``value`` is not an accepted iterable or contains a record of
        another type.
    """

    kinds = kind if isinstance(kind, tuple) else (kind,)
    if not kinds or any(not isinstance(item, type) for item in kinds):
        raise TypeError('records requires a type or a nonempty tuple of types.')
    items = _ordered_items(value, label)
    if any(not isinstance(item, kind) for item in items):
        raise DomainError(f"{label} contains a record of the wrong type.", field=label)
    return items


def strings(value: object, label: str, *, allow_blank: bool = False) -> tuple[str, ...]:
    """Return a tuple of UTF-8 encodable strings.

    Parameters
    ----------
    value : object
        Ordered, finite iterable of strings.
    label : str
        Field name used in error messages.
    allow_blank : bool, optional
        Accept empty or whitespace-only strings when true.

    Returns
    -------
    tuple of str
        Plain ``str`` copies in their original order.

    Raises
    ------
    DomainError
        If ``allow_blank`` is not a ``bool``, ``value`` is not an accepted
        iterable, or an item fails the text checks.
    """

    if not isinstance(allow_blank, bool):
        raise DomainError(f"{label} allow_blank must be a boolean.", field=label)
    items = _ordered_items(value, label)
    return tuple(_text_value(item, label, allow_blank=allow_blank) for item in items)


def name_value_pairs(value: object, label: str) -> tuple[tuple[str, float], ...]:
    """Return ``((name, float), ...)`` from an iterable of two-item pairs.

    Parameters
    ----------
    value : object
        Ordered, finite iterable whose items are ordered pairs
        ``(name, number)``.
    label : str
        Field name used in error messages.

    Returns
    -------
    tuple of (str, float)
        Validated pairs in their original order; names may repeat.

    Raises
    ------
    DomainError
        If an item is not a pair, a name is not non-blank text, or a number
        is not a finite real.
    """

    result: list[tuple[str, float]] = []
    for item in _ordered_items(value, label):
        name, number = _ordered_items(item, label, size=2)
        result.append((text(name, f"{label} name"), real(number, label)))
    return tuple(result)


def real_matrix(value: object, label: str) -> tuple[tuple[float, ...], ...]:
    """Return a rectangular matrix of finite floats.

    Parameters
    ----------
    value : object
        Ordered iterable of ordered rows of real numbers.
    label : str
        Field name used in error messages; rows are labelled
        ``label[index]``.

    Returns
    -------
    tuple of tuple of float
        The rows. A matrix with no rows, or with rows of length zero, is
        valid.

    Raises
    ------
    DomainError
        If a row is not an accepted iterable, an entry is not a finite real,
        or the rows differ in length.
    """

    rows = tuple(real_tuple(row, f'{label}[{index}]')
                 for index, row in enumerate(_ordered_items(value, label)))
    if rows and any(len(row) != len(rows[0]) for row in rows):
        raise DomainError(f"{label} must be rectangular.", field=label)
    return rows


def _scalar_kind(annotation):
    """Classify a dataclass field annotation without evaluating it.

    Returns ``"float"``, ``"int"``, ``"str"``, ``"optional_float"`` or
    ``None``. Both actual types and string annotations (as produced by
    ``from __future__ import annotations``) are recognized, including quoted,
    ``typing.``-qualified, ``Optional[float]`` and ``Union`` spellings.
    """
    if annotation is float:
        return 'float'
    if annotation is int:
        return 'int'
    if annotation is str:
        return 'str'
    if type(annotation) is str and annotation in ('float', 'int', 'str'):
        return annotation
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

    Parameters
    ----------
    instance : object
        A dataclass instance, typically ``self`` inside ``__post_init__``.
    label : str
        Prefix for error labels; each field is reported as
        ``"<label> <field name>"``.

    Raises
    ------
    TypeError
        If ``instance`` is not a dataclass instance (a programming error).
    DomainError
        If a recognized field fails :func:`real`, :func:`integer` or
        :func:`text`.

    Notes
    -----
    Fields annotated ``float`` become finite floats, ``float | None``
    likewise unless ``None``, ``int`` become exact integers and ``str`` must
    be non-blank strings. Other annotations are left to the class's own
    checks and are neither evaluated nor read. ``Optional[float]`` and quoted
    scalar spellings are recognized as well as the corresponding actual
    types. Every recognized field is validated before any field is written,
    so a failure leaves the instance unchanged. Call it after the fields are
    initialized and before model-domain checks.
    """

    if not is_dataclass(instance) or isinstance(instance, type):
        raise TypeError("coerce_fields requires a dataclass instance.")
    prepared = []
    for item in fields(instance):
        annotation = _scalar_kind(item.type)
        if annotation is None:
            continue
        value = getattr(instance, item.name)
        if annotation == 'float' and type(value) is float and isfinite(value):
            # Already normalized; finiteness is checked on every visit.
            # Subclasses and nonfinite values retain real()'s protocol.
            continue
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
    """Return ``value`` as a plain ``str`` suitable for ``json.loads``.

    Parameters
    ----------
    value : object
        Candidate JSON document text. ``bytes`` are refused so that the
        text decoding is always explicit.
    label : str
        Field name used in error messages.

    Returns
    -------
    str
        The same characters as ``value``.

    Raises
    ------
    InputError
        If ``value`` is not a ``str``.
    """

    if not isinstance(value, str):
        raise InputError(f"{label} must be JSON text.", field=label)
    return str.__str__(value)


# Exceptions that signal a decoded document cannot build a model: missing keys,
# wrong types, bad values (DomainError is a ValueError), short lists, integers
# beyond the binary64 range, and pathological nesting.
DOCUMENT_ERRORS = (
    KeyError, TypeError, ValueError, IndexError, AttributeError, OverflowError, RecursionError,
)


def finite_float_text(value: object, label: str) -> float:
    """Parse numeric text to a finite float without erasing nonzero input.

    Parameters
    ----------
    value : object
        Text accepted by :class:`float` (decimal or exponent notation).
    label : str
        Field name used in error messages.

    Returns
    -------
    float
        The parsed value.

    Raises
    ------
    InputError
        If ``value`` is not a ``str``, does not parse, overflows to infinity,
        names a non-finite value, or has a nonzero mantissa digit but rounds
        to zero.

    Notes
    -----
    The underflow test inspects the decimal digits of the mantissa (the
    part before any exponent marker), so ``"1e-400"`` is rejected while
    ``"0e5"`` is accepted as zero.
    """
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
    """``parse_float`` hook: decode a JSON number token with :func:`finite_float_text`."""
    return finite_float_text(token, 'JSON number')


def _unique_document_object(pairs):
    """``object_pairs_hook`` that refuses duplicate member names.

    RFC 8259 leaves the behaviour for duplicate names unspecified [1]_, and
    the standard-library decoder keeps the last value by default [2]_;
    refusing them prevents silent information loss.

    References
    ----------
    .. [1] IETF RFC 8259, "The JavaScript Object Notation (JSON) Data
           Interchange Format", https://www.rfc-editor.org/rfc/rfc8259
    .. [2] PYTHON-3.12-JSON. Python 3.12, "json - JSON encoder and decoder",
           https://docs.python.org/3.12/library/json.html
    """
    result = {}
    for key, value in pairs:
        if key in result:
            raise ValueError('Duplicate JSON member name.')
        result[key] = value
    return result


def _reject_document_constant(token):
    """``parse_constant`` hook: refuse ``NaN``, ``Infinity`` and ``-Infinity``.

    These tokens are outside the RFC 8259 grammar but are accepted by the
    standard-library decoder by default.
    """
    raise ValueError('Non-finite JSON constants are unsupported.')


def _check_document_values(payload):
    """Refuse unencodable text and non-finite numbers anywhere in a decoded document.

    Every string, keys included, must be UTF-8 encodable (which rejects lone
    surrogates produced by escape sequences) and every float must be finite.
    The walk is iterative and also covers fields the builder never reads.
    """
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
    """Parse a JSON document under the library's strict reading rules and build it.

    Parameters
    ----------
    value : object
        Document text; must be a ``str`` (``bytes`` are refused so that the
        decoding is explicit).
    label : str
        Document name used in error messages.
    build : callable
        Function that turns the decoded root object into a model.

    Returns
    -------
    object
        Whatever ``build`` returns.

    Raises
    ------
    InputError
        If the text is not a ``str``, violates the reading rules below, has a
        root that is not a JSON object, or ``build`` raises one of
        ``DOCUMENT_ERRORS`` (including integrity-hash failures).
    DomainError
        Passed through unchanged from ``build``; it already names the
        offending field.

    Notes
    -----
    Before ``build`` runs, the reader rejects duplicate member names, text
    that is not UTF-8 encodable, the tokens ``NaN`` and ``Infinity`` (outside
    the RFC 8259 grammar [1]_ but accepted by the Python decoder by
    default [2]_), numbers that overflow to infinity, and nonzero decimal
    literals that would round to zero. These are interoperability rules of
    this library, not authentication or physical validation; readers that
    use other parsers are not covered by them.

    References
    ----------
    .. [1] IETF RFC 8259, "The JavaScript Object Notation (JSON) Data
           Interchange Format", https://www.rfc-editor.org/rfc/rfc8259
    .. [2] PYTHON-3.12-JSON. Python 3.12, "json - JSON encoder and decoder",
           https://docs.python.org/3.12/library/json.html
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


def numeric_range_guard(field: str) -> Callable[[Callable[_P, _R]], Callable[_P, _R]]:
    """Make a numerical function fail closed on floating-point range errors.

    Python's float arithmetic raises ``OverflowError`` from ``**`` and
    ``math`` functions and ``ZeroDivisionError`` when a denominator underflows
    to zero; for inputs that are finite and well typed but extreme (``1e300``,
    ``1e-300``) those bare exceptions would escape a public function. The
    decorator converts exactly these two into :class:`DomainError` with code
    ``"numeric_range"``, naming ``field`` and the function; every other
    exception, including the package's own errors, passes through unchanged.

    Parameters
    ----------
    field : str
        Input the error is reported against (usually the first argument).

    Returns
    -------
    callable
        Decorator preserving the wrapped function's signature.

    References
    ----------
    .. [1] PYTHON-3.12-EXCEPTIONS. Python 3.12, "Built-in Exceptions",
           OverflowError and ZeroDivisionError,
           https://docs.python.org/3.12/library/exceptions.html
    """

    def decorate(function: Callable[_P, _R]) -> Callable[_P, _R]:
        @functools.wraps(function)
        def guarded(*args: _P.args, **kwargs: _P.kwargs) -> _R:
            try:
                return function(*args, **kwargs)
            except (OverflowError, ZeroDivisionError) as exc:
                raise DomainError(
                    f"{function.__qualname__}: inputs are outside the floating-point range of the "
                    f"model ({type(exc).__name__}: {exc}).",
                    code="numeric_range", field=field) from exc
        return guarded
    return decorate
