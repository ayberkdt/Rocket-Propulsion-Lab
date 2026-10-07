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


def json_text(value: object, label: str) -> str:
    """Return ``value`` if it is a string suitable for ``json.loads``."""

    if not isinstance(value, str):
        raise InputError(f"{label} must be JSON text.", field=label)
    return value
