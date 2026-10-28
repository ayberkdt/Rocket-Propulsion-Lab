"""NASA Glenn seven-term heat-capacity polynomial evaluation.

The representation is often called NASA9 because seven heat-capacity
coefficients are accompanied by two integration constants. Equations follow
NASA/TP-2002-211556 [1]_. The molar gas constant is that model's registered
``NASA_GLENN_2002_GAS_CONSTANT_J_MOL_K`` from Appendix A. The general SI
constant in :mod:`rocket_propulsion.core.constants` retains its SI value.

References
----------
.. [1] NASA/TP-2002-211556. NASA Glenn, "NASA Glenn Coefficients for
       Calculating Thermodynamic Properties of Individual Species".
       https://ntrs.nasa.gov/citations/20020085330
"""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass
from itertools import pairwise
from math import log
from typing import TypedDict

from rocket_propulsion.core.constants import NASA_GLENN_2002_GAS_CONSTANT_J_MOL_K
from rocket_propulsion.core.errors import DomainError
from rocket_propulsion.core.validation import (
    _ordered_items,
    coerce_fields,
    numeric_range_guard,
    real,
)

UNIVERSAL_GAS_CONSTANT_J_MOL_K = NASA_GLENN_2002_GAS_CONSTANT_J_MOL_K


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
        [1]_.

        References
        ----------
        .. [1] NASA/TP-2002-211556. NASA Glenn, "NASA Glenn Coefficients for
               Calculating Thermodynamic Properties of Individual Species".
               https://ntrs.nasa.gov/citations/20020085330
        """

        a1, a2, a3, a4, a5, a6, a7, _, _ = self.coefficients
        t = self._temperature(temperature_k)
        cp_over_r = (
            a1 / t**2
            + a2 / t
            + a3
            + a4 * t
            + a5 * t**2
            + a6 * t**3
            + a7 * t**4
        )
        return real(cp_over_r * UNIVERSAL_GAS_CONSTANT_J_MOL_K, 'cp_molar_j_mol_k')

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
        + a7 T^4/5 + b1/T`` [1]_.

        References
        ----------
        .. [1] NASA/TP-2002-211556. NASA Glenn, "NASA Glenn Coefficients for
               Calculating Thermodynamic Properties of Individual Species".
               https://ntrs.nasa.gov/citations/20020085330
        """

        a1, a2, a3, a4, a5, a6, a7, b1, _ = self.coefficients
        t = self._temperature(temperature_k)
        h_over_rt = (
            -a1 / t**2
            + a2 * log(t) / t
            + a3
            + a4 * t / 2.0
            + a5 * t**2 / 3.0
            + a6 * t**3 / 4.0
            + a7 * t**4 / 5.0
            + b1 / t
        )
        return real(h_over_rt * UNIVERSAL_GAS_CONSTANT_J_MOL_K * t, 'enthalpy_molar_j_mol')

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
        + a7 T^4/4 + b2`` [1]_.

        References
        ----------
        .. [1] NASA/TP-2002-211556. NASA Glenn, "NASA Glenn Coefficients for
               Calculating Thermodynamic Properties of Individual Species".
               https://ntrs.nasa.gov/citations/20020085330
        """

        a1, a2, a3, a4, a5, a6, a7, _, b2 = self.coefficients
        t = self._temperature(temperature_k)
        s_over_r = (
            -a1 / (2.0 * t**2)
            - a2 / t
            + a3 * log(t)
            + a4 * t
            + a5 * t**2 / 2.0
            + a6 * t**3 / 3.0
            + a7 * t**4 / 4.0
            + b2
        )
        return real(s_over_r * UNIVERSAL_GAS_CONSTANT_J_MOL_K, 'entropy_molar_j_mol_k')



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
