"""Elementary gas-phase reaction kinetics and a stiff Rosenbrock integrator.

A mechanism is a set of elementary reactions with modified-Arrhenius forward
rate coefficients, optional third-body efficiencies and optional Troe falloff.
Rates follow the mass-action law; reverse rate coefficients come from the
concentration equilibrium constant computed from the NASA Glenn polynomials of
:mod:`.species`, so forward and reverse rates are thermodynamically consistent
with the repository's equilibrium solver [1]_, [2]_, [3]_.

Equations (NASA RP-1328 [1]_ and the MFiX/CHEMKIN summary [2]_):

* ``k_f = A T^b exp(-E/(R T))`` (RP-1328 eq. (2.4); MFiX eq. (4.8));
* ``k_r = k_f / K_c`` and ``K_c = (p0/(R T))^dn exp(-dG0/(R T))`` with
  ``dn`` the net change of moles and ``dG0`` the standard-state Gibbs change
  (RP-1328 eqs. (2.6)--(2.10); MFiX eqs. (4.9)--(4.13)). ``p0`` is the 1 bar
  standard state of the NASA Glenn polynomials [3]_;
* third bodies ``[M] = sum_l eps_l C_l`` multiply both directions
  (RP-1328 eqs. (2.29)--(2.32); MFiX eqs. (4.14)--(4.17));
* falloff ``k = k_inf (P_r/(1+P_r)) F`` with ``P_r = k_0 [M]/k_inf`` and the
  Troe form of ``F`` (MFiX eq. (4.18)--(4.20) [2]_; Cantera [4]_, which
  states the base-10 logarithms);
* ``omega_i = sum_j (nu''_ij - nu'_ij)(R_j - R_-j)`` (RP-1328 eq. (2.33)).

The stiff integrator is the RODAS embedded Rosenbrock method of order 4(3)
with the coefficient set ``METH=1`` and the step-size controller of Hairer and
Wanner's ``rodas.f`` [5]_.

References
----------
.. [1] NASA-RP-1328. K. Radhakrishnan, "LSENS -- A General Chemical Kinetics
       and Sensitivity Analysis Code for Homogeneous Gas-Phase Reactions.
       I. Theory and Numerical Solution Procedures", NASA RP-1328, 1994,
       chapter 2. https://ntrs.nasa.gov/citations/19950008642
.. [2] MFIX-25.2-CHEMKIN. NETL, MFiX 25.2 documentation, section 4.14.4
       "Implementation of CHEMKIN mechanism in MFiX".
       https://mfix.netl.doe.gov/doc/mfix/25.2/html/model_setup/reactions/chemkin.html
.. [3] NASA/TP-2002-211556. NASA Glenn, "NASA Glenn Coefficients for
       Calculating Thermodynamic Properties of Individual Species".
       https://ntrs.nasa.gov/citations/20020085330
.. [4] CANTERA-3.2-RATES. Cantera 3.2.0 documentation, "Rate Constant
       Parameterizations" and "Reaction Rates".
       https://cantera.org/stable/reference/kinetics/rate-constants.html
.. [5] HAIRER-RODAS. E. Hairer and G. Wanner, ``rodas.f`` (version of
       October 28, 1996) and ``dc_decsol.f``.
       https://www.unige.ch/~hairer/prog/stiff/rodas.f
"""

from __future__ import annotations

import itertools
from collections.abc import Callable, Mapping, Sequence
from dataclasses import dataclass, field
from functools import lru_cache
from math import exp, inf, isfinite, log, log10, sqrt

from rocket_propulsion.core.constants import MOLAR_GAS_CONSTANT_J_MOL_K
from rocket_propulsion.core.errors import ConvergenceError, DomainError, InputError

from .mixtures import REFERENCE_PRESSURE_PA
from .nasa_polynomials import UNIVERSAL_GAS_CONSTANT_J_MOL_K as _THERMO_R
from .species import Species, get_species

CALORIE_TH_J = 4.184
"""Thermochemical calorie in J, exact (NIST SP 811, Appendix B.8 [1]_).

References
----------
.. [1] NIST-SP811-B8. NIST Guide to the SI, Appendix B.8, "Factors for Units
       Listed Alphabetically", calorie_th to joule.
       https://www.nist.gov/pml/special-publication-811/nist-guide-si-appendix-b-conversion-factors/nist-guide-si-appendix-b8
"""

_UNIT_ROUNDOFF = 2.220446049250313e-16


# ---------------------------------------------------------------------------
# Mechanism representation
# ---------------------------------------------------------------------------


@dataclass(frozen=True, slots=True)
class ArrheniusRate:
    """Modified Arrhenius rate coefficient ``k = A T^b exp(-E/(R T))`` [1]_.

    Attributes
    ----------
    pre_exponential_si : float
        ``A`` in SI units ``(m^3/mol)^(order-1) / s`` (with ``K^-b``).
    temperature_exponent : float
        ``b``.
    activation_energy_j_mol : float
        ``E`` in J/mol (may be negative, as in GRI-Mech 3.0).

    References
    ----------
    .. [1] NASA-RP-1328. K. Radhakrishnan, NASA RP-1328, 1994, eq. (2.4).
           https://ntrs.nasa.gov/citations/19950008642
    """

    pre_exponential_si: float
    temperature_exponent: float
    activation_energy_j_mol: float

    def __post_init__(self) -> None:
        for name in ("pre_exponential_si", "temperature_exponent", "activation_energy_j_mol"):
            if not isfinite(getattr(self, name)):
                raise InputError(f"Arrhenius {name} must be finite.", field=name)
        if self.pre_exponential_si < 0.0:
            raise DomainError("Arrhenius pre-exponential factor must be non-negative.",
                              field="pre_exponential_si")

    @classmethod
    def from_chemkin(cls, pre_exponential_cm_mol_s: float, temperature_exponent: float,
                     activation_energy_cal_mol: float, *, order: int) -> ArrheniusRate:
        """Convert CHEMKIN default units (cm, mol, s, cal/mol) to SI [1]_, [2]_.

        Parameters
        ----------
        pre_exponential_cm_mol_s : float
            ``A`` in ``(cm^3/mol)^(order-1) / s``.
        temperature_exponent : float
            ``b``.
        activation_energy_cal_mol : float
            ``E`` in thermochemical cal/mol.
        order : int
            Overall molecularity including a third body ``M``, for example
            2 for ``O+H2``, 3 for ``2O+M`` and ``H+2O2``.

        Returns
        -------
        ArrheniusRate
            ``A_SI = A (1e-6)^(order-1)``; ``E_SI = 4.184 E``.

        Raises
        ------
        InputError
            If ``order`` is not 1, 2 or 3.

        References
        ----------
        .. [1] MFIX-25.2-CHEMKIN. NETL, MFiX 25.2 documentation, section
               4.14.4.3.2 "Units in the Mechanism" (default cm-mole-sec-K and
               cal/mole).
               https://mfix.netl.doe.gov/doc/mfix/25.2/html/model_setup/reactions/chemkin.html
        .. [2] NIST-SP811-B8. NIST Guide to the SI, Appendix B.8.
               https://www.nist.gov/pml/special-publication-811/nist-guide-si-appendix-b-conversion-factors/nist-guide-si-appendix-b8
        """

        if order not in (1, 2, 3):
            raise InputError("Reaction order must be 1, 2 or 3.", field="order")
        return cls(pre_exponential_cm_mol_s * 1.0e-6 ** (order - 1), temperature_exponent,
                   activation_energy_cal_mol * CALORIE_TH_J)

    def __call__(self, temperature_k: float) -> float:
        """Return ``k(T)`` with SI R for activation energies converted to J/mol.

        This rate-parameter convention is separate from the NASA/CEA R
        convention used for thermodynamic equilibrium and concentrations.
        """
        return (self.pre_exponential_si * temperature_k ** self.temperature_exponent
                * exp(-self.activation_energy_j_mol / (MOLAR_GAS_CONSTANT_J_MOL_K * temperature_k)))


@dataclass(frozen=True, slots=True)
class TroeFalloff:
    """Troe falloff blending of low- and high-pressure limits [1]_, [2]_.

    Attributes
    ----------
    low_pressure : ArrheniusRate
        ``k_0`` (one order higher than the high-pressure limit).
    alpha, t3_k, t1_k : float
        Troe parameters ``alpha``, ``T***`` and ``T*``.
    t2_k : float or None
        Optional ``T**``; the third term of ``F_cent`` is omitted when ``None``.

    References
    ----------
    .. [1] MFIX-25.2-CHEMKIN. NETL, MFiX 25.2 documentation, eqs. (4.18)--(4.20)
           and the Troe form.
           https://mfix.netl.doe.gov/doc/mfix/25.2/html/model_setup/reactions/chemkin.html
    .. [2] CANTERA-3.2-RATES. Cantera 3.2.0, "The Troe Falloff Function".
           https://cantera.org/stable/reference/kinetics/rate-constants.html
    """

    low_pressure: ArrheniusRate
    alpha: float
    t3_k: float
    t1_k: float
    t2_k: float | None = None

    def __post_init__(self) -> None:
        if not (isfinite(self.alpha) and self.t3_k > 0.0 and self.t1_k > 0.0):
            raise DomainError("Troe parameters must be finite with positive T*** and T*.")
        if self.t2_k is not None and not isfinite(self.t2_k):
            raise DomainError("Troe T** must be finite.")

    def center(self, temperature_k: float) -> float:
        """Return ``F_cent(T)``."""
        value = ((1.0 - self.alpha) * exp(-temperature_k / self.t3_k)
                 + self.alpha * exp(-temperature_k / self.t1_k))
        if self.t2_k is not None:
            value += exp(-self.t2_k / temperature_k)
        return value

    def broadening(self, temperature_k: float, reduced_pressure: float) -> tuple[float, float]:
        """Return ``F`` and ``dF/dP_r`` at one temperature and reduced pressure.

        Notes
        -----
        ``log10 F = log10 F_cent / (1 + f1^2)``, ``f1 = (log10 P_r + C) /
        (N - 0.14 (log10 P_r + C))``, ``C = -0.4 - 0.67 log10 F_cent``,
        ``N = 0.75 - 1.27 log10 F_cent`` [1]_. The derivative follows by the
        chain rule.

        References
        ----------
        .. [1] CANTERA-3.2-RATES. Cantera 3.2.0, "The Troe Falloff Function".
               https://cantera.org/stable/reference/kinetics/rate-constants.html
        """

        if reduced_pressure <= 0.0:
            return 1.0, 0.0
        log_center = log10(max(self.center(temperature_k), 1.0e-300))
        c = -0.4 - 0.67 * log_center
        n = 0.75 - 1.27 * log_center
        u = log10(reduced_pressure) + c
        denominator = n - 0.14 * u
        f1 = u / denominator
        log_f = log_center / (1.0 + f1 * f1)
        broadening = 10.0 ** log_f
        df1_du = n / (denominator * denominator)
        dlogf_dlogpr = -log_center * 2.0 * f1 * df1_du / (1.0 + f1 * f1) ** 2
        return broadening, broadening * dlogf_dlogpr / reduced_pressure


@dataclass(frozen=True, slots=True)
class ThirdBody:
    """Third-body collision efficiencies ``[M] = sum_l eps_l C_l`` [1]_.

    Attributes
    ----------
    efficiencies : tuple of (str, float)
        Species with an efficiency different from the default.
    default_efficiency : float, default 1.0
        Efficiency of every other species.

    References
    ----------
    .. [1] NASA-RP-1328. K. Radhakrishnan, NASA RP-1328, 1994, eqs.
           (2.29)--(2.32). https://ntrs.nasa.gov/citations/19950008642
    """

    efficiencies: tuple[tuple[str, float], ...] = ()
    default_efficiency: float = 1.0

    def __post_init__(self) -> None:
        if any(not isfinite(value) or value < 0.0 for _, value in self.efficiencies):
            raise DomainError("Third-body efficiencies must be finite and non-negative.")


@dataclass(frozen=True, slots=True)
class Reaction:
    """One reversible elementary reaction with its source locator.

    Attributes
    ----------
    equation : str
        Reaction as written in the source.
    reactants, products : tuple of (str, int)
        Species and stoichiometric coefficients as written (a collider named
        explicitly, such as ``O2`` in ``H+2O2<=>HO2+O2``, appears on both
        sides and so cancels in the net stoichiometry but not in the rate).
    rate : ArrheniusRate
        Forward rate (the high-pressure limit when ``falloff`` is set).
    third_body : ThirdBody or None
        Generic third body ``+M`` or ``(+M)``.
    falloff : TroeFalloff or None
        Pressure-dependent ``(+M)`` reaction; requires ``third_body``.
    duplicate : bool
        Marked ``DUPLICATE`` in the source; the rates add.
    source : str
        Exact source locator of the rate parameters.
    """

    equation: str
    reactants: tuple[tuple[str, int], ...]
    products: tuple[tuple[str, int], ...]
    rate: ArrheniusRate
    third_body: ThirdBody | None = None
    falloff: TroeFalloff | None = None
    duplicate: bool = False
    source: str = ""

    def __post_init__(self) -> None:
        if not self.reactants or not self.products:
            raise InputError(f"Reaction {self.equation!r} needs reactants and products.")
        if any(count <= 0 for _, count in self.reactants + self.products):
            raise InputError(f"Reaction {self.equation!r} has a non-positive coefficient.")
        if self.falloff is not None and self.third_body is None:
            raise InputError(f"Falloff reaction {self.equation!r} needs a third body.")

    @property
    def net_stoichiometry(self) -> dict[str, int]:
        """Return ``nu''_i - nu'_i`` for every species with a non-zero change."""
        result: dict[str, int] = {}
        for key, count in self.products:
            result[key] = result.get(key, 0) + count
        for key, count in self.reactants:
            result[key] = result.get(key, 0) - count
        return {key: value for key, value in result.items() if value != 0}


@dataclass(frozen=True, slots=True)
class Mechanism:
    """A validated set of elementary reactions over a fixed species list.

    Attributes
    ----------
    name : str
        Mechanism name.
    species : tuple of str
        Species keys of the repository catalog [1]_, in state-vector order.
    reactions : tuple of Reaction
        Reactions; each must conserve every element.
    source, source_url : str
        Mechanism source and URL.
    notes : tuple of str
        Recorded omissions and adaptations.

    Raises
    ------
    InputError
        If a species is unknown, duplicated, or a reaction references a
        species outside the list.
    DomainError
        If a reaction does not conserve elements.

    References
    ----------
    .. [1] NASA/TP-2002-211556. NASA Glenn, "NASA Glenn Coefficients for
           Calculating Thermodynamic Properties of Individual Species".
           https://ntrs.nasa.gov/citations/20020085330
    """

    name: str
    species: tuple[str, ...]
    reactions: tuple[Reaction, ...]
    source: str = ""
    source_url: str = ""
    notes: tuple[str, ...] = field(default=())

    def __post_init__(self) -> None:
        if len(set(self.species)) != len(self.species) or not self.species:
            raise InputError("Mechanism species must be non-empty and unique.", field="species")
        catalog = {key: get_species(key) for key in self.species}
        known = set(self.species)
        for reaction in self.reactions:
            for key, _ in reaction.reactants + reaction.products:
                if key not in known:
                    raise InputError(f"Reaction {reaction.equation!r} uses unknown species {key}.",
                                     field="reactions")
            balance: dict[str, float] = {}
            for key, change in reaction.net_stoichiometry.items():
                for element, atoms in catalog[key].formula:
                    balance[element] = balance.get(element, 0.0) + change * atoms
            if any(abs(value) > 1.0e-12 for value in balance.values()):
                raise DomainError(f"Reaction {reaction.equation!r} does not conserve elements.",
                                  code="reaction_not_balanced", details={"balance": balance})

    @property
    def elements(self) -> tuple[str, ...]:
        """Sorted element symbols of the species."""
        return tuple(sorted({element for key in self.species
                             for element, _ in get_species(key).formula}))

    def element_matrix(self) -> tuple[tuple[float, ...], ...]:
        """Return ``a[e][k]``, atoms of element ``e`` in species ``k``."""
        formulas = [dict(get_species(key).formula) for key in self.species]
        return tuple(tuple(formula.get(element, 0.0) for formula in formulas)
                     for element in self.elements)


@dataclass(frozen=True, slots=True)
class RateEvaluation:
    """Forward and reverse rates and production rates at one state.

    Attributes
    ----------
    temperature_k : float
        Temperature.
    forward_rate_constants, equilibrium_constants : tuple of float
        ``k_f`` (falloff-blended where applicable, multiplier applied) and
        ``K_c`` in SI concentration units.
    forward_rates, reverse_rates : tuple of float
        ``R_j`` and ``R_-j`` in mol/(m^3 s), third body included.
    net_production_rates : tuple of float
        ``omega_i`` in mol/(m^3 s), species order of the mechanism.
    """

    temperature_k: float
    forward_rate_constants: tuple[float, ...]
    equilibrium_constants: tuple[float, ...]
    forward_rates: tuple[float, ...]
    reverse_rates: tuple[float, ...]
    net_production_rates: tuple[float, ...]


class _Compiled:
    """Index tables of a mechanism for fast evaluation (internal)."""

    def __init__(self, mechanism: Mechanism) -> None:
        index = {key: position for position, key in enumerate(mechanism.species)}
        self.species: tuple[Species, ...] = tuple(get_species(key) for key in mechanism.species)
        self.size = len(mechanism.species)
        self.reactions = mechanism.reactions
        self.reactant_orders: list[list[tuple[int, int]]] = []
        self.product_orders: list[list[tuple[int, int]]] = []
        self.net: list[list[tuple[int, int]]] = []
        self.delta_moles: list[int] = []
        self.efficiencies: list[list[float] | None] = []
        for reaction in mechanism.reactions:
            self.reactant_orders.append([(index[key], count) for key, count in reaction.reactants])
            self.product_orders.append([(index[key], count) for key, count in reaction.products])
            net = reaction.net_stoichiometry
            self.net.append([(index[key], change) for key, change in net.items()])
            self.delta_moles.append(sum(net.values()))
            if reaction.third_body is None:
                self.efficiencies.append(None)
            else:
                values = [reaction.third_body.default_efficiency] * self.size
                for key, value in reaction.third_body.efficiencies:
                    if key in index:
                        values[index[key]] = value
                self.efficiencies.append(values)
        self.intervals = tuple(
            tuple((interval.minimum_temperature_k, interval.maximum_temperature_k,
                   interval.coefficients) for interval in species.intervals)
            for species in self.species)

    def thermo(self, temperature_k: float) -> tuple[list[float], list[float], list[float]]:
        """Return ``cp/R``, ``h/(R T)`` and ``s0/R`` of every species.

        The NASA Glenn seven-term expressions of NASA/TP-2002-211556, as in
        :mod:`.nasa_polynomials`, evaluated with shared powers of ``T``; no
        extrapolation outside the polynomial ranges.
        """
        t = temperature_k
        if not isfinite(t) or t <= 0.0:
            raise DomainError("Temperature must be finite and positive.", field="temperature_k")
        t2, t3, t4 = t * t, t * t * t, t * t * t * t
        inverse, inverse2, log_t = 1.0 / t, 1.0 / t2, log(t)
        cp_r, h_rt, s_r = [], [], []
        for species, intervals in zip(self.species, self.intervals, strict=True):
            last = len(intervals) - 1
            for position, (low, high, coefficients) in enumerate(intervals):
                if low <= t and (t < high or (position == last and t <= high)):
                    break
            else:
                raise DomainError(
                    f"Temperature for {species.key} must be in "
                    f"[{intervals[0][0]:g}, {intervals[-1][1]:g}] K.",
                    code="property_temperature_out_of_range", field="temperature_k",
                    details={"species": species.key, "temperature_k": t})
            a1, a2, a3, a4, a5, a6, a7, b1, b2 = coefficients
            cp_r.append(a1 * inverse2 + a2 * inverse + a3 + a4 * t + a5 * t2 + a6 * t3 + a7 * t4)
            h_rt.append(-a1 * inverse2 + a2 * log_t * inverse + a3 + a4 * t / 2.0
                        + a5 * t2 / 3.0 + a6 * t3 / 4.0 + a7 * t4 / 5.0 + b1 * inverse)
            s_r.append(-a1 * inverse2 / 2.0 - a2 * inverse + a3 * log_t + a4 * t
                       + a5 * t2 / 2.0 + a6 * t3 / 3.0 + a7 * t4 / 4.0 + b2)
        return cp_r, h_rt, s_r

    def equilibrium_constants(self, temperature_k: float, h_rt: Sequence[float],
                              s_r: Sequence[float]) -> list[float]:
        """Return ``K_c`` on the NASA thermodynamic ``h/(R T)``, ``s0/R`` basis."""
        log_scale = log(REFERENCE_PRESSURE_PA / (_THERMO_R * temperature_k))
        values = []
        for net, delta in zip(self.net, self.delta_moles, strict=True):
            exponent = -sum(change * (h_rt[i] - s_r[i]) for i, change in net) + delta * log_scale
            values.append(exp(min(exponent, 700.0)))
        return values

    def rates(self, temperature_k: float, concentrations: Sequence[float],
              multiplier: float, constants: Sequence[float]
              ) -> tuple[list[float], list[float], list[float], list[float]]:
        """Return ``k_f``, ``R_j``, ``R_-j`` and ``omega_i`` (see :func:`evaluate_rates`)."""
        forward_constants, forward, reverse = [], [], []
        production = [0.0] * self.size
        for j, reaction in enumerate(self.reactions):
            efficiencies = self.efficiencies[j]
            third = (None if efficiencies is None
                     else sum(e * c for e, c in zip(efficiencies, concentrations, strict=True)))
            k_forward = _forward_constant(reaction, temperature_k, third)[0] * multiplier
            rate_forward = k_forward
            rate_reverse = k_forward / constants[j]
            for i, order in self.reactant_orders[j]:
                rate_forward *= concentrations[i] if order == 1 else concentrations[i] ** order
            for i, order in self.product_orders[j]:
                rate_reverse *= concentrations[i] if order == 1 else concentrations[i] ** order
            if third is not None and reaction.falloff is None:
                rate_forward *= third
                rate_reverse *= third
            forward_constants.append(k_forward)
            forward.append(rate_forward)
            reverse.append(rate_reverse)
            net_rate = rate_forward - rate_reverse
            for i, change in self.net[j]:
                production[i] += change * net_rate
        return forward_constants, forward, reverse, production

    def rates_jacobian(self, temperature_k: float, concentrations: Sequence[float],
                       multiplier: float, constants: Sequence[float]) -> list[list[float]]:
        """Return ``d omega_i/d C_k`` (see :func:`production_rate_jacobian`)."""
        size = self.size
        jacobian = [[0.0] * size for _ in range(size)]
        for j, reaction in enumerate(self.reactions):
            efficiencies = self.efficiencies[j]
            third = (None if efficiencies is None
                     else sum(e * c for e, c in zip(efficiencies, concentrations, strict=True)))
            k_forward, dk_dm = _forward_constant(reaction, temperature_k, third)
            k_forward *= multiplier
            dk_dm *= multiplier
            inverse_kc = 1.0 / constants[j]
            dq = [0.0] * size
            uses_third = third is not None and reaction.falloff is None
            factor = third if third is not None and uses_third else 1.0
            for orders, sign in ((self.reactant_orders[j], 1.0),
                                 (self.product_orders[j], -inverse_kc)):
                derivative_keys = {k for k, _ in orders}
                for k in derivative_keys:
                    dq[k] += sign * k_forward * factor * _mass_action_derivative(
                        orders, concentrations, k)
            if efficiencies is not None and (uses_third or reaction.falloff is not None):
                difference = (_mass_action(self.reactant_orders[j], concentrations)
                              - inverse_kc * _mass_action(self.product_orders[j], concentrations))
                scale = k_forward if uses_third else dk_dm
                for k in range(size):
                    dq[k] += scale * efficiencies[k] * difference
            for i, change in self.net[j]:
                row = jacobian[i]
                for k in range(size):
                    if dq[k] != 0.0:
                        row[k] += change * dq[k]
        return jacobian


@lru_cache(maxsize=16)
def _compile(mechanism: Mechanism) -> _Compiled:
    return _Compiled(mechanism)


def _check_state(temperature_k: float, concentrations: Sequence[float], size: int) -> None:
    if not isfinite(temperature_k) or temperature_k <= 0.0:
        raise DomainError("Temperature must be finite and positive.", field="temperature_k")
    if len(concentrations) != size:
        raise InputError(f"Expected {size} concentrations.", field="concentrations")
    if any(not isfinite(value) for value in concentrations):
        raise InputError("Concentrations must be finite.", field="concentrations")


def standard_gibbs_over_rt(mechanism: Mechanism, temperature_k: float) -> tuple[float, ...]:
    """Return ``g0_i/(R T) = h_i/(R T) - s0_i/R`` for every species [1]_.

    Parameters
    ----------
    mechanism : Mechanism
        Mechanism whose species order is used.
    temperature_k : float
        Temperature in K inside every species' polynomial range.

    Returns
    -------
    tuple of float
        Dimensionless standard-state (1 bar) Gibbs energies.

    Raises
    ------
    DomainError
        If the temperature is outside a species' polynomial range.

    References
    ----------
    .. [1] NASA/TP-2002-211556. NASA Glenn, "NASA Glenn Coefficients for
           Calculating Thermodynamic Properties of Individual Species".
           https://ntrs.nasa.gov/citations/20020085330
    """

    _, h_rt, s_r = _compile(mechanism).thermo(temperature_k)
    return tuple(h - s for h, s in zip(h_rt, s_r, strict=True))


def equilibrium_constants(mechanism: Mechanism, temperature_k: float) -> tuple[float, ...]:
    """Return the concentration equilibrium constants ``K_c`` of every reaction [1]_, [2]_.

    Parameters
    ----------
    mechanism : Mechanism
        Mechanism.
    temperature_k : float
        Temperature in K.

    Returns
    -------
    tuple of float
        ``K_c = (p0/(R T))^dn exp(-sum_i nu_i g0_i/(R T))`` in ``(mol/m^3)^dn``
        with ``p0 = 1 bar``, the standard state of the NASA Glenn polynomials.

    References
    ----------
    .. [1] NASA-RP-1328. K. Radhakrishnan, NASA RP-1328, 1994, eqs.
           (2.7)--(2.10). https://ntrs.nasa.gov/citations/19950008642
    .. [2] MFIX-25.2-CHEMKIN. NETL, MFiX 25.2 documentation, eqs. (4.10)--(4.13).
           https://mfix.netl.doe.gov/doc/mfix/25.2/html/model_setup/reactions/chemkin.html
    """

    compiled = _compile(mechanism)
    _, h_rt, s_r = compiled.thermo(temperature_k)
    return tuple(compiled.equilibrium_constants(temperature_k, h_rt, s_r))


def _forward_constant(reaction: Reaction, temperature_k: float,
                      third_body: float | None) -> tuple[float, float]:
    """Return ``k_f`` and ``dk_f/d[M]`` (zero unless falloff)."""
    high = reaction.rate(temperature_k)
    if reaction.falloff is None:
        return high, 0.0
    low = reaction.falloff.low_pressure(temperature_k)
    concentration = max(third_body or 0.0, 0.0)
    if high <= 0.0:
        return 0.0, 0.0
    reduced = low * concentration / high
    broadening, dbroadening = reaction.falloff.broadening(temperature_k, reduced)
    blend = reduced / (1.0 + reduced)
    value = high * blend * broadening
    dblend = 1.0 / (1.0 + reduced) ** 2
    dvalue_dreduced = high * (dblend * broadening + blend * dbroadening)
    return value, dvalue_dreduced * low / high


def evaluate_rates(mechanism: Mechanism, temperature_k: float,
                   concentrations: Sequence[float], *, rate_multiplier: float = 1.0
                   ) -> RateEvaluation:
    """Evaluate forward, reverse and net species production rates [1]_.

    Parameters
    ----------
    mechanism : Mechanism
        Mechanism.
    temperature_k : float
        Temperature in K.
    concentrations : sequence of float
        Molar concentrations ``C_i`` in mol/m^3, mechanism species order.
    rate_multiplier : float, default 1.0
        Common factor on every forward rate coefficient; reverse rates follow
        through ``K_c`` so equilibrium is unchanged. ``0`` freezes the
        chemistry; large values approach local equilibrium.

    Returns
    -------
    RateEvaluation
        Rate constants, directional rates and ``omega_i`` in mol/(m^3 s).

    Raises
    ------
    InputError
        If the state vector has the wrong length or is not finite.
    DomainError
        If the temperature is outside the thermodynamic data range or the
        multiplier is negative.

    Notes
    -----
    ``R_j = k_f [M] prod C^nu'``, ``R_-j = (k_f/K_c) [M] prod C^nu''`` with
    ``[M]`` present only for ``+M`` reactions; for ``(+M)`` falloff reactions
    ``[M]`` enters ``k_f`` through ``P_r`` instead [1]_, [2]_.

    References
    ----------
    .. [1] NASA-RP-1328. K. Radhakrishnan, NASA RP-1328, 1994, eqs.
           (2.16)--(2.33). https://ntrs.nasa.gov/citations/19950008642
    .. [2] MFIX-25.2-CHEMKIN. NETL, MFiX 25.2 documentation, eqs. (4.4)--(4.20).
           https://mfix.netl.doe.gov/doc/mfix/25.2/html/model_setup/reactions/chemkin.html
    """

    compiled = _compile(mechanism)
    _check_state(temperature_k, concentrations, compiled.size)
    if not isfinite(rate_multiplier) or rate_multiplier < 0.0:
        raise DomainError("Rate multiplier must be finite and non-negative.",
                          field="rate_multiplier")
    constants = equilibrium_constants(mechanism, temperature_k)
    forward_constants, forward, reverse, production = compiled.rates(
        temperature_k, concentrations, rate_multiplier, constants)
    return RateEvaluation(temperature_k, tuple(forward_constants), constants, tuple(forward),
                          tuple(reverse), tuple(production))


def production_rate_jacobian(mechanism: Mechanism, temperature_k: float,
                             concentrations: Sequence[float], *, rate_multiplier: float = 1.0
                             ) -> tuple[tuple[float, ...], ...]:
    """Return the analytic Jacobian ``d omega_i / d C_k`` at fixed temperature [1]_.

    Parameters
    ----------
    mechanism : Mechanism
        Mechanism.
    temperature_k : float
        Temperature in K.
    concentrations : sequence of float
        Molar concentrations in mol/m^3.
    rate_multiplier : float, default 1.0
        Common forward-rate factor, as in :func:`evaluate_rates`.

    Returns
    -------
    tuple of tuple of float
        Square matrix in mechanism species order, units 1/s.

    Notes
    -----
    Differentiates the mass-action products, the third-body factor
    ``[M] = sum eps_l C_l`` and, for falloff reactions, ``k_f(P_r)`` through
    ``dP_r/dC_l = eps_l k_0/k_inf`` and the Troe ``dF/dP_r``. The equations
    differentiated are those of :func:`evaluate_rates` [1]_. Concentrations
    are used as given; a zero concentration raised to the power one has
    derivative one.

    References
    ----------
    .. [1] NASA-RP-1328. K. Radhakrishnan, NASA RP-1328, 1994, eqs.
           (2.29)--(2.33). https://ntrs.nasa.gov/citations/19950008642
    """

    compiled = _compile(mechanism)
    _check_state(temperature_k, concentrations, compiled.size)
    if not isfinite(rate_multiplier) or rate_multiplier < 0.0:
        raise DomainError("Rate multiplier must be finite and non-negative.",
                          field="rate_multiplier")
    constants = equilibrium_constants(mechanism, temperature_k)
    jacobian = compiled.rates_jacobian(temperature_k, concentrations, rate_multiplier, constants)
    return tuple(tuple(row) for row in jacobian)


def _mass_action(orders: list[tuple[int, int]], concentrations: Sequence[float]) -> float:
    value = 1.0
    for i, order in orders:
        value *= concentrations[i] ** order
    return value


def _mass_action_derivative(orders: list[tuple[int, int]], concentrations: Sequence[float],
                            k: int) -> float:
    total_order = sum(order for i, order in orders if i == k)
    if total_order == 0:
        return 0.0
    value = total_order * concentrations[k] ** (total_order - 1)
    for i, order in orders:
        if i != k:
            value *= concentrations[i] ** order
    return value


def concentrations_from_mole_fractions(mechanism: Mechanism, temperature_k: float,
                                       pressure_pa: float,
                                       mole_fractions: Mapping[str, float]) -> tuple[float, ...]:
    """Return ``C_i = x_i p/(R T)`` on the NASA thermodynamic basis, mol/m^3.

    Parameters
    ----------
    mechanism : Mechanism
        Mechanism defining the species order.
    temperature_k, pressure_pa : float
        State.
    mole_fractions : mapping of str to float
        Mole fractions; missing species are zero, unknown keys are rejected.

    Returns
    -------
    tuple of float
        Concentrations in mechanism order.

    Raises
    ------
    InputError
        If a key is not a mechanism species.
    DomainError
        If the state is not positive.
    """

    if temperature_k <= 0.0 or pressure_pa <= 0.0:
        raise DomainError("Temperature and pressure must be positive.")
    unknown = set(mole_fractions) - set(mechanism.species)
    if unknown:
        raise InputError(f"Species not in mechanism: {sorted(unknown)}.", field="mole_fractions")
    total = pressure_pa / (_THERMO_R * temperature_k)
    return tuple(total * float(mole_fractions.get(key, 0.0)) for key in mechanism.species)


# ---------------------------------------------------------------------------
# Dense linear algebra (LU with partial pivoting)
# ---------------------------------------------------------------------------


def _lu_factor(matrix: list[list[float]]) -> tuple[list[list[float]], list[int]] | None:
    """Return the in-place LU factors and pivots, or ``None`` if singular."""
    size = len(matrix)
    lu = [row[:] for row in matrix]
    pivots = list(range(size))
    for column in range(size):
        pivot = max(range(column, size), key=lambda row: abs(lu[row][column]))
        if lu[pivot][column] == 0.0 or not isfinite(lu[pivot][column]):
            return None
        if pivot != column:
            lu[column], lu[pivot] = lu[pivot], lu[column]
            pivots[column], pivots[pivot] = pivots[pivot], pivots[column]
        diagonal = lu[column][column]
        row_c = lu[column]
        for row in range(column + 1, size):
            row_r = lu[row]
            factor = row_r[column] / diagonal
            if factor != 0.0:
                row_r[column] = factor
                for k in range(column + 1, size):
                    row_r[k] -= factor * row_c[k]
            else:
                row_r[column] = 0.0
    return lu, pivots


def _lu_solve(factors: tuple[list[list[float]], list[int]], rhs: Sequence[float]) -> list[float]:
    lu, pivots = factors
    size = len(lu)
    y = [rhs[pivots[i]] for i in range(size)]
    for i in range(size):
        row = lu[i]
        y[i] -= sum(row[k] * y[k] for k in range(i))
    for i in range(size - 1, -1, -1):
        row = lu[i]
        y[i] = (y[i] - sum(row[k] * y[k] for k in range(i + 1, size))) / row[i]
    return y


# ---------------------------------------------------------------------------
# RODAS stiff integrator
# ---------------------------------------------------------------------------

# Coefficients of rodas.f subroutine ROCOE, METH=1 (lines 895-930 of the
# version of October 28, 1996) [HAIRER-RODAS].
_C2, _C3, _C4 = 0.386, 0.21, 0.63
_D1, _D2, _D3, _D4 = 0.25, -0.1043, 0.1035, -0.3620000000000023e-01
_A21 = 0.1544e+01
_A31, _A32 = 0.9466785280815826, 0.2557011698983284
_A41, _A42, _A43 = 0.3314825187068521e+01, 0.2896124015972201e+01, 0.9986419139977817
_A51, _A52, _A53, _A54 = (0.1221224509226641e+01, 0.6019134481288629e+01,
                          0.1253708332932087e+02, -0.6878860361058950)
_C21 = -0.56688e+01
_C31, _C32 = -0.2430093356833875e+01, -0.2063599157091915
_C41, _C42, _C43 = -0.1073529058151375, -0.9594562251023355e+01, -0.2047028614809616e+02
_C51, _C52, _C53, _C54 = (0.7496443313967647e+01, -0.1024680431464352e+02,
                          -0.3399990352819905e+02, 0.1170890893206160e+02)
_C61, _C62, _C63, _C64, _C65 = (0.8083246795921522e+01, -0.7981132988064893e+01,
                                -0.3152159432874371e+02, 0.1631930543123136e+02,
                                -0.6058818238834054e+01)
_GAMMA = 0.25


@dataclass(frozen=True, slots=True)
class StiffSolution:
    """Output of :func:`solve_stiff`.

    Attributes
    ----------
    t : tuple of float
        Output abscissae (``t_eval`` or every accepted step).
    y : tuple of tuple of float
        Solution at each output abscissa.
    accepted_steps, rejected_steps : int
        Step counts (rejections after the first accepted step, as in
        ``rodas.f``).
    function_evaluations, jacobian_evaluations, decompositions : int
        Work counters.
    last_step : float
        Step size proposed for continuation.
    """

    t: tuple[float, ...]
    y: tuple[tuple[float, ...], ...]
    accepted_steps: int
    rejected_steps: int
    function_evaluations: int
    jacobian_evaluations: int
    decompositions: int
    last_step: float


class StepRejected(Exception):
    """Raised by a right-hand side to reject the current trial step.

    The integrator halves the step, like the singular-matrix branch of
    ``rodas.f``; five consecutive occurrences end the integration.
    """


def solve_stiff(
    function: Callable[[float, Sequence[float]], Sequence[float]],
    t_span: tuple[float, float],
    y0: Sequence[float],
    *,
    rtol: float = 1.0e-6,
    atol: float | Sequence[float] = 1.0e-10,
    jacobian: Callable[[float, Sequence[float]], Sequence[Sequence[float]]] | None = None,
    autonomous: bool = False,
    t_eval: Sequence[float] | None = None,
    first_step: float = 1.0e-6,
    max_step: float = inf,
    max_steps: int = 100_000,
    step_callback: Callable[[float, tuple[float, ...]], None] | None = None,
) -> StiffSolution:
    """Integrate ``y' = f(t, y)`` with the RODAS Rosenbrock method of order 4(3) [1]_.

    Parameters
    ----------
    function : callable
        ``f(t, y)`` returning a sequence of the same length as ``y``. It may
        raise :class:`StepRejected` to reject a trial step.
    t_span : tuple of float
        ``(t0, t_end)``; integration may run backwards.
    y0 : sequence of float
        Initial state.
    rtol : float, default 1e-6
        Relative tolerance.
    atol : float or sequence of float, default 1e-10
        Absolute tolerance, scalar or per component.
    jacobian : callable, optional
        ``J(t, y)`` as a full matrix; finite differences with the increment
        ``sqrt(uround max(1e-5, |y_i|))`` of ``rodas.f`` otherwise.
    autonomous : bool, default False
        Skip the ``df/dt`` term (estimated by forward differences otherwise).
    t_eval : sequence of float, optional
        Monotonic output abscissae inside the span; every accepted step is
        returned when omitted.
    first_step : float, default 1e-6
        Initial step magnitude (the ``rodas.f`` default).
    max_step : float, default inf
        Largest step magnitude.
    max_steps : int, default 100000
        Step budget.
    step_callback : callable, optional
        Called as ``callback(t, y)`` after every accepted step, with ``y`` a
        tuple.

    Returns
    -------
    StiffSolution
        Output and work statistics.

    Raises
    ------
    InputError
        If the inputs are inconsistent.
    ConvergenceError
        If the step budget is exhausted, the step becomes too small, the
        iteration matrix is singular or the step is rejected five consecutive
        times by :class:`StepRejected` (codes ``stiff_max_steps``,
        ``stiff_step_too_small``, ``stiff_singular_matrix``).

    Notes
    -----
    One step solves ``(1/(h gamma) I - J) k_i = f(t + c_i h, y + sum a_ij k_j)
    + sum (c_ij/h) k_j + h d_i f_t`` for six stages; the fifth stage gives the
    embedded third-order solution and the sixth the fourth-order solution
    whose difference is the error estimate. The error norm is
    ``sqrt(mean((k_6/sk)^2))``, ``sk = atol + rtol max(|y|, |y_new|)``. The new
    step is ``h / max(1/6, min(5, err^0.25/0.9))`` with the predictive
    (Gustafsson) controller; after a rejection the step may not grow. All
    constants and the stage algebra are those of ``rodas.f`` subroutines
    ``ROSCOR`` and ``ROCOE`` (``METH=1``) with the identity mass matrix and a
    full Jacobian (``dc_decsol.f`` ``DECOMR``/``SLVROD`` case 1) [1]_. One
    deviation: the counter of singular (or :class:`StepRejected`) stages is
    reset after every accepted step, so five *consecutive* failures end the
    integration, whereas ``rodas.f`` counts them over the whole call.

    References
    ----------
    .. [1] HAIRER-RODAS. E. Hairer and G. Wanner, ``rodas.f`` (version of
           October 28, 1996), subroutines ROSCOR and ROCOE, and
           ``dc_decsol.f``. https://www.unige.ch/~hairer/prog/stiff/rodas.f
    """

    t0, t_end = float(t_span[0]), float(t_span[1])
    if not (isfinite(t0) and isfinite(t_end)) or t0 == t_end:
        raise InputError("t_span must be two distinct finite values.", field="t_span")
    y = [float(value) for value in y0]
    size = len(y)
    if size == 0 or any(not isfinite(value) for value in y):
        raise InputError("Initial state must be non-empty and finite.", field="y0")
    if not isfinite(rtol) or rtol <= 0.0:
        raise InputError("rtol must be positive.", field="rtol")
    atol_vector = ([float(atol)] * size if isinstance(atol, (int, float))
                   else [float(value) for value in atol])
    if len(atol_vector) != size or any(not isfinite(v) or v < 0.0 for v in atol_vector):
        raise InputError("atol must be non-negative with one value per component.", field="atol")
    direction = 1.0 if t_end > t0 else -1.0
    targets = [t_end] if t_eval is None else [float(value) for value in t_eval]
    if any((target - t0) * direction <= 0.0 for target in targets[:1]) and t_eval is not None:
        raise InputError("t_eval must lie after t0 in the integration direction.", field="t_eval")
    for previous, current in itertools.pairwise(targets):
        if (current - previous) * direction <= 0.0:
            raise InputError("t_eval must be strictly monotonic.", field="t_eval")
    if (targets[-1] - t_end) * direction > 0.0:
        raise InputError("t_eval must not exceed t_end.", field="t_eval")

    counters = {"fcn": 0, "jac": 0, "dec": 0, "acc": 0, "rej": 0, "steps": 0}

    def call(t: float, state: Sequence[float]) -> list[float]:
        counters["fcn"] += 1
        values = [float(value) for value in function(t, state)]
        if len(values) != size:
            raise InputError("Right-hand side returned the wrong length.", field="function")
        return values

    output_t: list[float] = []
    output_y: list[tuple[float, ...]] = []
    t = t0
    h = abs(first_step) if first_step and abs(first_step) > 10.0 * _UNIT_ROUNDOFF else 1.0e-6
    for target in targets:
        t, y, h = _rodas_core(call, t, y, target, h, rtol, atol_vector, jacobian, autonomous,
                              max_step, max_steps, counters, step_callback,
                              output_t if t_eval is None else None,
                              output_y if t_eval is None else None)
        if t_eval is not None:
            output_t.append(t)
            output_y.append(tuple(y))
    return StiffSolution(tuple(output_t), tuple(output_y), counters["acc"], counters["rej"],
                         counters["fcn"], counters["jac"], counters["dec"], h)


def _rodas_core(call, t, y, t_end, h, rtol, atol, jacobian, autonomous, max_step, max_steps,
                counters, step_callback, record_t, record_y):
    """Advance from ``t`` to ``t_end``; returns ``(t, y, h_proposed)`` (rodas.f ROSCOR)."""
    size = len(y)
    fac1, fac2, safe = 5.0, 1.0 / 6.0, 0.9
    direction = 1.0 if t_end > t else -1.0
    h_max = min(abs(max_step), abs(t_end - t))
    h = direction * min(abs(h), h_max)
    reject = False
    last = False
    singular = 0
    h_accepted = 0.0
    error_accepted = 0.0
    accepted_here = 0
    h_optimal = h
    while True:
        if counters["steps"] > max_steps:
            raise ConvergenceError("Stiff integrator exceeded its step budget.",
                                   code="stiff_max_steps", details={"t": t})
        if 0.1 * abs(h) <= abs(t) * _UNIT_ROUNDOFF:
            raise ConvergenceError("Stiff integrator step became too small.",
                                   code="stiff_step_too_small", details={"t": t, "h": h})
        if last:
            return t, y, h_optimal
        h_optimal = h
        if (t + h * 1.0001 - t_end) * direction >= 0.0:
            h = t_end - t
            last = True
        f0 = call(t, y)
        counters["jac"] += 1
        if jacobian is not None:
            jac = [[float(v) for v in row] for row in jacobian(t, y)]
        else:
            jac = [[0.0] * size for _ in range(size)]
            for i in range(size):
                saved = y[i]
                delta = sqrt(_UNIT_ROUNDOFF * max(1.0e-5, abs(saved)))
                y[i] = saved + delta
                try:
                    shifted = call(t, y)
                finally:
                    y[i] = saved
                for row in range(size):
                    jac[row][i] = (shifted[row] - f0[row]) / delta
        if autonomous:
            ft = None
        else:
            delta = sqrt(_UNIT_ROUNDOFF * max(1.0e-5, abs(t)))
            shifted = call(t + delta, y)
            ft = [(shifted[i] - f0[i]) / delta for i in range(size)]
        while True:  # label 2 of rodas.f: compute the stages with the current h
            counters["steps"] += 1
            outcome = _rodas_step(call, t, y, h, f0, jac, ft, counters)
            if outcome is None:
                singular += 1
                if singular >= 5:
                    raise ConvergenceError("Stiff integrator matrix repeatedly singular or "
                                           "step repeatedly rejected by the right-hand side.",
                                           code="stiff_singular_matrix",
                                           details={"t": t, "h": h})
                h *= 0.5
                reject = True
                last = False
                continue
            y_new, error_vector = outcome
            error = 0.0
            for i in range(size):
                scale = atol[i] + rtol * max(abs(y[i]), abs(y_new[i]))
                error += (error_vector[i] / scale) ** 2
            error = sqrt(error / size)
            if not isfinite(error):
                singular += 1
                if singular >= 5:
                    raise ConvergenceError("Stiff integrator produced a non-finite error.",
                                           code="stiff_singular_matrix", details={"t": t})
                h *= 0.5
                reject = True
                last = False
                continue
            factor = max(fac2, min(fac1, error ** 0.25 / safe))
            h_new = h / factor
            if error <= 1.0:
                singular = 0
                counters["acc"] += 1
                accepted_here += 1
                if accepted_here > 1:
                    gustafsson = (h_accepted / h) * (error ** 2 / error_accepted) ** 0.25 / safe
                    gustafsson = max(fac2, min(fac1, gustafsson))
                    factor = max(factor, gustafsson)
                    h_new = h / factor
                h_accepted = h
                error_accepted = max(1.0e-2, error)
                y = y_new
                t = t + h
                if step_callback is not None:
                    step_callback(t, tuple(y))
                if record_t is not None:
                    record_t.append(t)
                    record_y.append(tuple(y))
                if abs(h_new) > h_max:
                    h_new = direction * h_max
                if reject:
                    h_new = direction * min(abs(h_new), abs(h))
                reject = False
                h = h_new
                break
            reject = True
            last = False
            h = h_new
            if accepted_here >= 1:
                counters["rej"] += 1


def _rodas_step(call, t, y, h, f0, jac, ft, counters):
    """Six RODAS stages; returns ``(y_new, k6)`` or ``None`` on a singular matrix."""
    size = len(y)
    inverse = 1.0 / (h * _GAMMA)
    matrix = [[-value for value in row] for row in jac]
    for i in range(size):
        matrix[i][i] += inverse
    factors = _lu_factor(matrix)
    counters["dec"] += 1
    if factors is None:
        return None

    def rhs(values, hd, combination):
        out = list(values)
        if hd != 0.0 and ft is not None:
            out = [out[i] + hd * ft[i] for i in range(size)]
        if combination is not None:
            out = [out[i] + combination[i] for i in range(size)]
        return _lu_solve(factors, out)

    try:
        k1 = rhs(f0, h * _D1, None)
        stage = [y[i] + _A21 * k1[i] for i in range(size)]
        f = call(t + _C2 * h, stage)
        k2 = rhs(f, h * _D2, [_C21 / h * k1[i] for i in range(size)])
        stage = [y[i] + _A31 * k1[i] + _A32 * k2[i] for i in range(size)]
        f = call(t + _C3 * h, stage)
        k3 = rhs(f, h * _D3, [(_C31 * k1[i] + _C32 * k2[i]) / h for i in range(size)])
        stage = [y[i] + _A41 * k1[i] + _A42 * k2[i] + _A43 * k3[i] for i in range(size)]
        f = call(t + _C4 * h, stage)
        k4 = rhs(f, h * _D4, [(_C41 * k1[i] + _C42 * k2[i] + _C43 * k3[i]) / h
                              for i in range(size)])
        stage = [y[i] + _A51 * k1[i] + _A52 * k2[i] + _A53 * k3[i] + _A54 * k4[i]
                 for i in range(size)]
        f = call(t + h, stage)
        k5 = rhs(f, 0.0, [(_C51 * k1[i] + _C52 * k2[i] + _C53 * k3[i] + _C54 * k4[i]) / h
                          for i in range(size)])
        stage = [stage[i] + k5[i] for i in range(size)]
        f = call(t + h, stage)
        k6 = rhs(f, 0.0, [(_C61 * k1[i] + _C62 * k2[i] + _C63 * k3[i] + _C64 * k4[i]
                           + _C65 * k5[i]) / h for i in range(size)])
    except StepRejected:
        return None
    y_new = [stage[i] + k6[i] for i in range(size)]
    if any(not isfinite(value) for value in y_new):
        return None
    return y_new, k6
