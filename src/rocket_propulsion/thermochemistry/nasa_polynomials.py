"""NASA Glenn seven-term heat-capacity polynomial evaluation.

The representation is often called NASA9 because seven heat-capacity
coefficients are accompanied by two integration constants. Equations follow
NASA/TP-2002-211556 exactly [1]_. The molar gas constant is the registered
``MOLAR_GAS_CONSTANT_J_MOL_K`` of :mod:`rocket_propulsion.core.constants`.

References
----------
.. [1] NASA/TP-2002-211556. NASA Glenn, "NASA Glenn Coefficients for
       Calculating Thermodynamic Properties of Individual Species".
       https://ntrs.nasa.gov/citations/20020085330
"""

from __future__ import annotations

from dataclasses import dataclass
from math import log
from typing import TypedDict

from rocket_propulsion.core.constants import MOLAR_GAS_CONSTANT_J_MOL_K
from rocket_propulsion.core.errors import DomainError

UNIVERSAL_GAS_CONSTANT_J_MOL_K = MOLAR_GAS_CONSTANT_J_MOL_K


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
        If the interval bounds are not positive and increasing.

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
        if self.minimum_temperature_k <= 0.0:
            raise DomainError("NASA polynomial minimum temperature must be positive.")
        if self.maximum_temperature_k <= self.minimum_temperature_k:
            raise DomainError("NASA polynomial temperature interval must be increasing.")

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

        upper_ok = (
            temperature_k <= self.maximum_temperature_k
            if include_upper
            else temperature_k < self.maximum_temperature_k
        )
        return self.minimum_temperature_k <= temperature_k and upper_ok

    def cp_molar_j_mol_k(self, temperature_k: float) -> float:
        """Return standard-state molar heat capacity at constant pressure.

        Parameters
        ----------
        temperature_k : float
            Temperature in K inside the interval (not checked here).

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
        t = temperature_k
        cp_over_r = (
            a1 / t**2
            + a2 / t
            + a3
            + a4 * t
            + a5 * t**2
            + a6 * t**3
            + a7 * t**4
        )
        return cp_over_r * UNIVERSAL_GAS_CONSTANT_J_MOL_K

    def enthalpy_molar_j_mol(self, temperature_k: float) -> float:
        """Return standard-state molar enthalpy including heat of formation.

        Parameters
        ----------
        temperature_k : float
            Temperature in K inside the interval (not checked here).

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
        t = temperature_k
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
        return h_over_rt * UNIVERSAL_GAS_CONSTANT_J_MOL_K * t

    def entropy_molar_j_mol_k(self, temperature_k: float) -> float:
        """Return standard-state molar entropy at one reference bar.

        Parameters
        ----------
        temperature_k : float
            Temperature in K inside the interval (not checked here).

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
        t = temperature_k
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
        return s_over_r * UNIVERSAL_GAS_CONSTANT_J_MOL_K



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
    minimum = float(record["minimum_temperature_k"])
    maximum = float(record["maximum_temperature_k"])
    values = tuple(float(value) for value in record["coefficients"])
    if len(values) != 9:
        raise DomainError("A NASA-9 interval needs nine coefficients (a1..a7, b1, b2).",
                          field="coefficients")
    a1, a2, a3, a4, a5, a6, a7, b1, b2 = values
    return NasaPolynomialInterval(
        minimum_temperature_k=minimum,
        maximum_temperature_k=maximum,
        coefficients=(a1, a2, a3, a4, a5, a6, a7, b1, b2),
    )
