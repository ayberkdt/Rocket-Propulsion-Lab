"""Owned nested data and caloric identities for thermodynamic records.

These checks establish record integrity, not the authenticity of a source
or the physical validity of a polynomial over its whole fitted range.
"""

from collections.abc import Mapping
from math import isclose, isfinite, ulp
from types import MappingProxyType
from typing import TypeVar

from rocket_propulsion.core.errors import DomainError
from rocket_propulsion.core.validation import _ordered_items, real, text

_T = TypeVar('_T')


def _owned_formula(value: object) -> tuple[tuple[str, float], ...]:
    formula = []
    seen = set()
    for entry in _ordered_items(value, 'formula'):
        symbol, amount = _ordered_items(entry, 'formula entry', size=2)
        symbol = text(symbol, 'element').strip().capitalize()
        amount = real(amount, 'atom count')
        # E is the signed electron inventory used for ionic formulae.
        if amount == 0 or (symbol != 'E' and amount < 0) or symbol in seen:
            raise DomainError('Formula entries must be unique with positive atom counts (signed E allowed).',
                              field='formula')
        seen.add(symbol)
        formula.append((symbol, amount))
    if not formula:
        raise DomainError('A formula must contain at least one element.', field='formula')
    return tuple(formula)


def _owned_catalog(value: Mapping[str, _T], kind: type[_T], label: str) -> Mapping[str, _T]:
    if not isinstance(value, Mapping):
        raise DomainError(f'{label} must be a mapping.', field=label)
    result = {}
    for key, record in value.items():
        key = text(key, label).strip().upper()
        record_key = getattr(record, 'key', None)
        if (not isinstance(record, kind) or not isinstance(record_key, str)
                or record_key.strip().upper() != key or key in result):
            raise DomainError(f'{label} needs unique canonical keys matching its records.', field=label)
        result[key] = record
    return MappingProxyType(result)


def _property_identity(actual: float, expected: float, label: str) -> None:
    if not isfinite(expected) or not isclose(actual, expected, rel_tol=5e-13, abs_tol=0):
        raise DomainError(f'Inconsistent thermodynamic property: {label}.',
                          code='invalid_thermodynamic_property', field=label)


def _caloric_identity(cp: float, cv: float, gas_constant: float, gamma: float) -> None:
    if not (gas_constant > 0 and cp >= cv > 0 and gamma >= 1):
        raise DomainError('Rounded ideal-gas states require cp >= cv > 0, R > 0 and gamma >= 1.',
                          code='invalid_thermodynamic_property')
    # Separately rounded cp and cv can coincide even when their retained
    # positive difference is representable. Check the gap against its
    # operand resolution, not just against the much smaller value of R.
    rounding = ulp(cp) + ulp(cv) + ulp(gas_constant)
    if not isclose(cp - cv, gas_constant, rel_tol=5e-13, abs_tol=rounding):
        raise DomainError('Inconsistent thermodynamic property: cp-cv=R.',
                          code='invalid_thermodynamic_property', field='cp-cv=R')
    _property_identity(cp / cv, gamma, 'gamma=cp/cv')
