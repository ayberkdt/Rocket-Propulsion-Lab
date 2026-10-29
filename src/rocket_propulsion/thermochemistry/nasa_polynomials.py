"""NASA Glenn seven-term heat-capacity polynomial evaluation.

The representation is often called NASA9 because seven heat-capacity
coefficients are accompanied by two integration constants. Equations follow
NASA/TP-2002-211556 [1]_. The molar gas constant is that model's registered
``NASA_GLENN_2002_GAS_CONSTANT_J_MOL_K`` from Appendix A. The general SI
constant in :mod:`rocket_propulsion.core.constants` retains its SI value.

Rational terms retain the normalized binary inputs exactly. A bounded adaptive
Decimal log enclosure resolves each exported binary64 rounding; it is numerical
evidence for the retained fit, not a thermodynamic accuracy bound.

References
----------
.. [1] NASA/TP-2002-211556. NASA Glenn, "NASA Glenn Coefficients for
       Calculating Thermodynamic Properties of Individual Species".
       https://ntrs.nasa.gov/citations/20020085330
"""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass
from decimal import MAX_EMAX, MIN_EMIN, ROUND_HALF_EVEN, Context, Decimal
from fractions import Fraction
from functools import lru_cache
from itertools import pairwise
from math import inf, log
from typing import TypedDict

from rocket_propulsion.core.constants import NASA_GLENN_2002_GAS_CONSTANT_J_MOL_K
from rocket_propulsion.core.errors import ConvergenceError, DomainError
from rocket_propulsion.core.validation import (
    _ordered_items,
    coerce_fields,
    numeric_range_guard,
    real,
)

UNIVERSAL_GAS_CONSTANT_J_MOL_K = NASA_GLENN_2002_GAS_CONSTANT_J_MOL_K


def _rounded(value):
    """Allow temporary infinities while checking an enclosing interval."""
    try:
        return float(value)
    except OverflowError:
        return inf if value > 0 else -inf


@lru_cache(maxsize=4096)
def _log_enclosure(temperature, precision):
    """Return ln(T) at a precision and its conservative one-ULP radius, exactly.

    Decimal.ln is correctly rounded, ROUND_HALF_EVEN, in an isolated context
    whose settings, traps and flags are all explicit. One decimal ULP is
    conservative and is retained as an exact Fraction. The enclosure depends
    only on the binary temperature and precision, so it is shared.
    """
    context = Context(prec=precision, Emax=MAX_EMAX, Emin=MIN_EMIN,
                      rounding=ROUND_HALF_EVEN, traps=[], flags=[], capitals=1, clamp=0)
    logarithm = context.ln(Decimal.from_float(temperature))
    unit = Fraction(Decimal((0, (1,), logarithm.adjusted()-precision+1)))
    return Fraction(logarithm), unit


def _export(base, log_coefficient, temperature, label):
    """Resolve one rounding of base + coefficient*ln(T), or fail explicitly."""
    if not log_coefficient or temperature == 1.0:
        result = _rounded(base)
        nonzero = bool(base)
    else:
        for precision in (48, 96, 192, 384, 768, 1536, 3072):
            logarithm, unit = _log_enclosure(temperature, precision)
            center = base + log_coefficient*logarithm
            radius = abs(log_coefficient)*unit
            lower, upper = _rounded(center-radius), _rounded(center+radius)
            if lower == upper:
                result = lower
                break
        else:
            raise ConvergenceError('NASA9 logarithm rounding is unresolved.',
                                   code='nasa9_rounding_unresolved', field=label)
        # For rational positive T != 1, ln(T) is not a nonzero rational;
        # a nonzero rational log coefficient cannot cancel a rational base.
        nonzero = True
    if result in (inf, -inf) or (nonzero and result == 0.0):
        raise DomainError('NASA9 property is outside the nonzero finite output range.',
                          code='nasa9_output_range', field=label)
    return result


@lru_cache(maxsize=8192)
def _exact_terms(coefficients, temperature):
    """Return molar (rational, log coefficient) pairs for cp, H and S.

    Cached on the owned coefficient tuple and binary temperature: the result
    is an immutable exact function of those inputs, so reuse cannot change it.
    """
    t = Fraction(temperature)
    a, b, c, d, e, f, g, h, s = map(Fraction, coefficients)
    r = Fraction(UNIVERSAL_GAS_CONSTANT_J_MOL_K)
    cp = r*(a/t**2+b/t+c+d*t+e*t**2+f*t**3+g*t**4)
    enthalpy = r*(-a/t+c*t+d*t**2/2+e*t**3/3+f*t**4/4+g*t**5/5+h)
    entropy = r*(-a/(2*t**2)-b/t+d*t+e*t**2/2+f*t**3/3+g*t**4/4+s)
    return (cp, Fraction(0)), (enthalpy, r*b), (entropy, r*c)


@dataclass(frozen=True, slots=True)
class NasaPolynomialInterval:
    """One temperature interval of a NASA Glenn polynomial fit.

    Attributes
    ----------
    minimum_temperature_k, maximum_temperature_k : float
        Interval bounds in K; the minimum is positive and below the maximum.
    coefficients : tuple of float
        ``(a1, ..., a7, b1, b2)``: seven heat-capacity coefficients and the
        enthalpy and entropy integration constants [1]_.

    Raises
    ------
    DomainError
        If bounds are not finite, positive and increasing, or there are not
        exactly nine finite real coefficients. Nested input is snapshotted.

    References
    ----------
    .. [1] NASA/TP-2002-211556. NASA Glenn, "NASA Glenn Coefficients for
           Calculating Thermodynamic Properties of Individual Species".
           https://ntrs.nasa.gov/citations/20020085330
    """

    minimum_temperature_k: float
    maximum_temperature_k: float
    coefficients: tuple[float, float, float, float, float, float, float, float, float]

    def __post_init__(self) -> None:
        coerce_fields(self, 'NASA interval')
        if self.minimum_temperature_k <= 0.0:
            raise DomainError("NASA polynomial minimum temperature must be positive.")
        if self.maximum_temperature_k <= self.minimum_temperature_k:
            raise DomainError("NASA polynomial temperature interval must be increasing.")
        values = _ordered_items(self.coefficients, 'coefficients', size=9)
        object.__setattr__(self, 'coefficients', tuple(real(value, 'coefficients') for value in values))

    def _temperature(self, value: float) -> float:
        temperature = real(value, 'temperature_k')
        if temperature <= 0:
            raise DomainError('Polynomial evaluation needs a positive temperature.', field='temperature_k')
        return temperature

    def contains(self, temperature_k: float, *, include_upper: bool = True) -> bool:
        """Return whether this interval owns a temperature.

        Parameters
        ----------
        temperature_k : float
            Temperature in K.
        include_upper : bool, default True
            Treat the upper bound as part of the interval. Adjacent intervals
            share bounds, so all but the last interval of a species exclude
            it.

        Returns
        -------
        bool
            ``True`` if the temperature lies in the interval.
        """

        temperature_k = real(temperature_k, 'temperature_k')
        if not isinstance(include_upper, bool):
            raise DomainError('include_upper must be boolean.', field='include_upper')
        upper_ok = (
            temperature_k <= self.maximum_temperature_k
            if include_upper
            else temperature_k < self.maximum_temperature_k
        )
        return self.minimum_temperature_k <= temperature_k and upper_ok

    def _terms(self, temperature):
        """Return exact molar (rational, log coefficient) pairs for cp, H, S."""
        return _exact_terms(self.coefficients, self._temperature(temperature))

    def _binary64_molar(self, temperature_k: float) -> tuple[float, float, float]:
        """Return molar (cp, H, S) by binary64 evaluation, for solver internals.

        These are main's original floating expressions. Iterative solvers use
        them for speed; exported records use the exactly rounded methods.
        """
        a1, a2, a3, a4, a5, a6, a7, b1, b2 = self.coefficients
        t = self._temperature(temperature_k)
        cp_over_r = a1 / t**2 + a2 / t + a3 + a4 * t + a5 * t**2 + a6 * t**3 + a7 * t**4
        h_over_rt = (-a1 / t**2 + a2 * log(t) / t + a3 + a4 * t / 2.0 + a5 * t**2 / 3.0
                     + a6 * t**3 / 4.0 + a7 * t**4 / 5.0 + b1 / t)
        s_over_r = (-a1 / (2.0 * t**2) - a2 / t + a3 * log(t) + a4 * t + a5 * t**2 / 2.0
                    + a6 * t**3 / 3.0 + a7 * t**4 / 4.0 + b2)
        return (real(cp_over_r * UNIVERSAL_GAS_CONSTANT_J_MOL_K, 'cp_molar_j_mol_k'),
                real(h_over_rt * UNIVERSAL_GAS_CONSTANT_J_MOL_K * t, 'enthalpy_molar_j_mol'),
                real(s_over_r * UNIVERSAL_GAS_CONSTANT_J_MOL_K, 'entropy_molar_j_mol_k'))

    @numeric_range_guard('temperature_k')
    def cp_molar_j_mol_k(self, temperature_k: float) -> float:
        """Return standard-state molar heat capacity at constant pressure.

        Parameters
        ----------
        temperature_k : float
            Positive finite real temperature in K. This polynomial primitive
            does not enforce the fit range; Species.properties does.

        Returns
        -------
        float
            ``cp`` in J/(mol K).

        Notes
        -----
        ``cp/R = a1 T^-2 + a2 T^-1 + a3 + a4 T + a5 T^2 + a6 T^3 + a7 T^4``
        [1]_, evaluated exactly on the binary inputs and rounded once.

        References
        ----------
        .. [1] NASA/TP-2002-211556. NASA Glenn, "NASA Glenn Coefficients for
               Calculating Thermodynamic Properties of Individual Species".
               https://ntrs.nasa.gov/citations/20020085330
        """

        t = self._temperature(temperature_k)
        base, log_coefficient = self._terms(t)[0]
        return _export(base, log_coefficient, t, 'cp_molar_j_mol_k')

    @numeric_range_guard('temperature_k')
    def enthalpy_molar_j_mol(self, temperature_k: float) -> float:
        """Return standard-state molar enthalpy including heat of formation.

        Parameters
        ----------
        temperature_k : float
            Positive finite real temperature in K. This polynomial primitive
            also supports reference-enthalpy checks outside the fitted range;
            Species.properties enforces the physical-state range.

        Returns
        -------
        float
            ``H`` in J/mol on the NASA Glenn assigned-enthalpy basis.

        Notes
        -----
        ``H/(R T) = -a1 T^-2 + a2 ln(T)/T + a3 + a4 T/2 + a5 T^2/3 + a6 T^3/4
        + a7 T^4/5 + b1/T`` [1]_. The ``ln(T)`` term is enclosed and the
        result is rounded once.

        References
        ----------
        .. [1] NASA/TP-2002-211556. NASA Glenn, "NASA Glenn Coefficients for
               Calculating Thermodynamic Properties of Individual Species".
               https://ntrs.nasa.gov/citations/20020085330
        """

        t = self._temperature(temperature_k)
        base, log_coefficient = self._terms(t)[1]
        return _export(base, log_coefficient, t, 'enthalpy_molar_j_mol')

    @numeric_range_guard('temperature_k')
    def entropy_molar_j_mol_k(self, temperature_k: float) -> float:
        """Return standard-state molar entropy at one reference bar.

        Parameters
        ----------
        temperature_k : float
            Positive finite real temperature in K. This polynomial primitive
            does not enforce the fit range; Species.properties does.

        Returns
        -------
        float
            ``S`` in J/(mol K) at the 1 bar standard-state pressure.

        Notes
        -----
        ``S/R = -a1 T^-2/2 - a2/T + a3 ln(T) + a4 T + a5 T^2/2 + a6 T^3/3
        + a7 T^4/4 + b2`` [1]_. The ``ln(T)`` term is enclosed and the
        result is rounded once.

        References
        ----------
        .. [1] NASA/TP-2002-211556. NASA Glenn, "NASA Glenn Coefficients for
               Calculating Thermodynamic Properties of Individual Species".
               https://ntrs.nasa.gov/citations/20020085330
        """

        t = self._temperature(temperature_k)
        base, log_coefficient = self._terms(t)[2]
        return _export(base, log_coefficient, t, 'entropy_molar_j_mol_k')



class NasaIntervalRecord(TypedDict):
    """One NASA-9 interval as stored in the packaged JSON data files.

    Attributes
    ----------
    minimum_temperature_k, maximum_temperature_k : float
        Interval bounds in K.
    coefficients : list of float
        ``[a1, ..., a7, b1, b2]`` as printed in the source [1]_.

    References
    ----------
    .. [1] NASA/TP-2002-211556. NASA Glenn, "NASA Glenn Coefficients for
           Calculating Thermodynamic Properties of Individual Species".
           https://ntrs.nasa.gov/citations/20020085330
    """

    minimum_temperature_k: float
    maximum_temperature_k: float
    coefficients: list[float]


def interval_from_record(record: NasaIntervalRecord) -> NasaPolynomialInterval:
    """Build a :class:`NasaPolynomialInterval` from one packaged JSON record.

    Parameters
    ----------
    record : NasaIntervalRecord
        Interval bounds and the nine coefficients ``a1..a7, b1, b2``.

    Returns
    -------
    NasaPolynomialInterval

    Raises
    ------
    DomainError
        If the record does not hold exactly nine coefficients, or the bounds
        are not positive and increasing.

    References
    ----------
    .. [1] NASA/TP-2002-211556. NASA Glenn, "NASA Glenn Coefficients for
           Calculating Thermodynamic Properties of Individual Species".
           https://ntrs.nasa.gov/citations/20020085330
    """
    if not isinstance(record, Mapping):
        raise DomainError('A NASA interval record must be a mapping.', field='record')
    try:
        minimum = record["minimum_temperature_k"]
        maximum = record["maximum_temperature_k"]
        values = _ordered_items(record['coefficients'], 'coefficients', size=9)
    except KeyError as exc:
        raise DomainError('A NASA interval record is missing a required field.', field='record') from exc
    a1, a2, a3, a4, a5, a6, a7, b1, b2 = values
    return NasaPolynomialInterval(
        minimum_temperature_k=minimum,
        maximum_temperature_k=maximum,
        coefficients=(a1, a2, a3, a4, a5, a6, a7, b1, b2),
    )


def _owned_intervals(value: object) -> tuple[NasaPolynomialInterval, ...]:
    """Own ordered, nonempty, contiguous intervals; never sort or close gaps."""
    items = _ordered_items(value, 'intervals')
    if not items or any(not isinstance(item, NasaPolynomialInterval) for item in items):
        raise DomainError('Expected nonempty NASA polynomial intervals.', field='intervals')
    for previous, following in pairwise(items):
        if previous.maximum_temperature_k != following.minimum_temperature_k:
            raise DomainError('NASA intervals must be ordered and contiguous, without gaps or overlap.',
                              field='intervals')
    return items
