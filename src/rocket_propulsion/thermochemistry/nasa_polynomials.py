"""NASA Glenn seven-term heat-capacity polynomial evaluation.

The representation is often called NASA9 because seven heat-capacity
coefficients are accompanied by two integration constants. Equations follow
NASA/TP-2002-211556 exactly.
"""

from __future__ import annotations

from dataclasses import dataclass
from math import log

from rocket_propulsion.core.errors import DomainError

UNIVERSAL_GAS_CONSTANT_J_MOL_K = 8.31446261815324


@dataclass(frozen=True, slots=True)
class NasaPolynomialInterval:
    """One temperature interval of a NASA Glenn polynomial fit."""

    minimum_temperature_k: float
    maximum_temperature_k: float
    coefficients: tuple[float, float, float, float, float, float, float, float, float]

    def __post_init__(self) -> None:
        if self.minimum_temperature_k <= 0.0:
            raise DomainError("NASA polynomial minimum temperature must be positive.")
        if self.maximum_temperature_k <= self.minimum_temperature_k:
            raise DomainError("NASA polynomial temperature interval must be increasing.")

    def contains(self, temperature_k: float, *, include_upper: bool = True) -> bool:
        """Return whether this interval owns a temperature."""

        upper_ok = (
            temperature_k <= self.maximum_temperature_k
            if include_upper
            else temperature_k < self.maximum_temperature_k
        )
        return self.minimum_temperature_k <= temperature_k and upper_ok

    def cp_molar_j_mol_k(self, temperature_k: float) -> float:
        """Return standard-state molar heat capacity at constant pressure."""

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
        """Return standard-state molar enthalpy including heat of formation."""

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
        """Return standard-state molar entropy at one reference bar."""

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

