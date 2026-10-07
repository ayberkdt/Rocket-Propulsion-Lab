"""Strict input coercion for public numerical boundaries.

Every public constructor and evaluation that receives numbers from outside
the package passes them through these helpers.  They guarantee two things:

* malformed input (``None``, text, ``bool``, containers, NaN, infinity) fails
  closed with :class:`DomainError` naming the field, never with a bare
  ``TypeError`` from deep inside an equation; and
* accepted numbers are normalized to ``float``.  Hashed exchange artifacts are
  canonical JSON, where ``1`` and ``1.0`` serialize differently; normalizing at
  construction makes an object built with integers hash exactly like the same
  object reopened from JSON, so it can be written and read back.

``bool`` is rejected on purpose: ``True`` is an ``int`` in Python, and a
boolean passed where a physical quantity is expected is always a caller bug.
"""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import fields, is_dataclass
from math import isfinite
from numbers import Integral, Real

from .errors import DomainError, InputError


def real(value: object, label: str) -> float:
    """Return ``value`` as a finite ``float`` or raise :class:`DomainError`."""

    if isinstance(value, bool) or not isinstance(value, Real):
        raise DomainError(f"{label} must be a real number.", field=label)
    number = float(value)
    if not isfinite(number):
        raise DomainError(f"{label} must be finite.", field=label)
    return number


def integer(value: object, label: str) -> int:
    """Return ``value`` as an ``int`` (``bool`` and floats are rejected)."""

    if isinstance(value, bool) or not isinstance(value, Integral):
        raise DomainError(f"{label} must be an integer.", field=label)
    return int(value)


def text(value: object, label: str) -> str:
    """Return a non-blank string or raise :class:`DomainError`."""

    if not isinstance(value, str) or not value.strip():
        raise DomainError(f"{label} must be a non-empty string.", field=label)
    return value


def vector3(value: object, label: str) -> tuple[float, float, float]:
    """Return a three-component finite ``float`` tuple."""

    if isinstance(value, (str, bytes, Mapping)) or not hasattr(value, "__len__"):
        raise DomainError(f"{label} must contain three finite components.", field=label)
    items = tuple(value)  # type: ignore[call-overload]
    if len(items) != 3:
        raise DomainError(f"{label} must contain three finite components.", field=label)
    return (real(items[0], label), real(items[1], label), real(items[2], label))


def real_tuple(value: object, label: str) -> tuple[float, ...]:
    """Return a tuple of finite floats from any non-string sequence."""

    if isinstance(value, (str, bytes, Mapping)) or not hasattr(value, "__iter__"):
        raise DomainError(f"{label} must be a sequence of real numbers.", field=label)
    return tuple(real(item, label) for item in value)  # type: ignore[union-attr]


def real_mapping(value: object, label: str) -> dict[str, float]:
    """Return a ``{non-blank str: finite float}`` copy of a mapping (``None`` -> empty)."""

    if value is None:
        return {}
    if not isinstance(value, Mapping):
        raise DomainError(f"{label} must be a mapping of names to real numbers.", field=label)
    result: dict[str, float] = {}
    for key, item in value.items():
        name = text(key, f"{label} key")
        result[name] = real(item, f"{label}[{name!r}]")
    return result


def records(value: object, kind: type | tuple[type, ...], label: str) -> tuple:
    """Return a tuple of records of ``kind`` from any non-string iterable."""

    if isinstance(value, (str, bytes, Mapping)) or not hasattr(value, "__iter__"):
        raise DomainError(f"{label} must be a sequence.", field=label)
    items = tuple(value)  # type: ignore[call-overload]
    if any(not isinstance(item, kind) for item in items):
        raise DomainError(f"{label} contains a record of the wrong type.", field=label)
    return items


def strings(value: object, label: str, *, allow_blank: bool = False) -> tuple[str, ...]:
    """Return a tuple of strings (non-blank unless ``allow_blank``)."""

    if isinstance(value, (str, bytes, Mapping)) or not hasattr(value, "__iter__"):
        raise DomainError(f"{label} must be a sequence of strings.", field=label)
    items = tuple(value)  # type: ignore[call-overload]
    for item in items:
        if not isinstance(item, str) or (not allow_blank and not item.strip()):
            raise DomainError(f"{label} must contain only non-empty strings.", field=label)
    return items


def name_value_pairs(value: object, label: str) -> tuple[tuple[str, float], ...]:
    """Return ``((name, float), ...)`` from a sequence of two-item pairs."""

    if isinstance(value, (str, bytes, Mapping)) or not hasattr(value, "__iter__"):
        raise DomainError(f"{label} must be a sequence of (name, value) pairs.", field=label)
    result: list[tuple[str, float]] = []
    for item in value:  # type: ignore[union-attr]
        if isinstance(item, (str, bytes)) or not hasattr(item, "__len__") or len(item) != 2:
            raise DomainError(f"{label} must be a sequence of (name, value) pairs.", field=label)
        result.append((text(item[0], f"{label} name"), real(item[1], label)))
    return tuple(result)


def real_matrix(value: object, label: str) -> tuple[tuple[float, ...], ...]:
    """Return a tuple-of-tuples of finite floats (rows may be any sequence)."""

    if isinstance(value, (str, bytes, Mapping)) or not hasattr(value, "__iter__"):
        raise DomainError(f"{label} must be a matrix.", field=label)
    return tuple(real_tuple(row, label) for row in value)  # type: ignore[union-attr]


def coerce_fields(instance: object, label: str) -> None:
    """Validate and normalize a frozen dataclass's annotated scalar fields in place.

    Fields annotated ``float`` become finite floats, ``float | None`` likewise
    unless ``None``, ``int`` become integers and ``str`` must be non-blank
    strings.  Other annotations are left to the class's own checks.  Call it
    first in ``__post_init__`` so later checks only see clean values.
    """

    if not is_dataclass(instance):
        raise TypeError("coerce_fields requires a dataclass instance.")
    for item in fields(instance):
        annotation = item.type if isinstance(item.type, str) else getattr(
            item.type, "__name__", str(item.type)
        )
        annotation = annotation.replace(" ", "")
        value = getattr(instance, item.name)
        name = f"{label} {item.name}"
        if annotation == "float":
            object.__setattr__(instance, item.name, real(value, name))
        elif annotation in {"float|None", "None|float"}:
            if value is not None:
                object.__setattr__(instance, item.name, real(value, name))
        elif annotation == "int":
            object.__setattr__(instance, item.name, integer(value, name))
        elif annotation == "str":
            text(value, name)


def json_text(value: object, label: str) -> str:
    """Return ``value`` if it is a string suitable for ``json.loads``."""

    if not isinstance(value, str):
        raise InputError(f"{label} must be JSON text.", field=label)
    return value
