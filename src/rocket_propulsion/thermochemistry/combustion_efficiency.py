"""Product composition, stoichiometry, combustion-efficiency models and O/F optimisation.

Everything between the injector and the throat that sets delivered performance:

* **Stoichiometry.** Oxidation-state (valence) bookkeeping of NASA CEA: the
  chemical equivalence ratio ``r`` and the weight equivalence ratio ``phi``
  of RP-1311 Part I equations (9.14)-(9.19) [1]_, with the valences of CEA's
  atomic data (H +1, C +4, N 0, O -2) [2]_; total-reactant density, equation
  (9.13) [1]_, and density impulse ``I_s d`` [3]_.
* **Product composition.** Full equilibrium mole and mass fractions at the
  chamber, throat and exit of the CEA-equivalent infinite-area rocket solve,
  a radical inventory, element and mass-fraction balances, and a dissociation
  energy split built on the assigned-enthalpy definition
  ``H(T) = Delta_f H(298.15) + [H(T) - H(298.15)]`` of RP-1311 equation (4.2)
  [4]_.
* **Combustion efficiency.** The c*, thrust-coefficient and vacuum-Isp
  efficiencies of NASA/TP-1999-208522 equations (6)-(11) [5]_; the
  vaporization-limited c* efficiency of Heidmann and Priem, NACA TN 4219
  equations (1)-(7) [6]_; the two-streamtube Rupe mixing model and the
  c*- and Isp-based energy-release efficiencies of the JANNAF simplified
  procedure as implemented in ROCCID, NASA CR-187109 equations (2.4)-(2.11)
  [7]_; the SP-8089 maximum-mixing orifice correlation [8]_; characteristic
  length and stay time, SP-125 equations (4-3)-(4-4) [9]_; and the
  momentum pressure-loss relation of TP-1999-208522 equation (4) [5]_.
* **Isp loss chain.** ODE -> ODK (kinetic) -> two-dimensional inviscid
  (divergence and shock) -> boundary layer -> energy release, in the order of
  the JANNAF methodology shown in TP-1999-208522 Figures 20 and 24 [5]_, with
  the boundary-layer thrust decrement converted to an impulse decrement and
  subtracted [10]_. The factors after kinetics are handed to
  :func:`rocket_propulsion.propulsion.losses.calculate_loss_budget`.
* **O/F optimisation.** Golden-section maximisation of vacuum Isp or density
  impulse evaluated with the CEA-equivalent shifting-equilibrium solver.

References
----------
.. [1] NASA-RP1311-P1-EQRATIO. Gordon and McBride, NASA RP-1311 Part I (1994),
       chapter 9, equations (9.10)-(9.19), printed pp. 42-43.
       https://ntrs.nasa.gov/citations/19950013764
.. [2] NASA-CEA-SRC-VALENCE. NASA Glenn, CEA source ``source/atomic_data.f90``
       and ``source/mixture.f90``, commit 2f79a647. https://github.com/nasa/cea
.. [3] NASA-SP125-DENSITY-IMPULSE. Huzel and Huang, NASA SP-125 (1967),
       equation (1-45), printed p. 20. https://ntrs.nasa.gov/citations/19710019929
.. [4] NASA-RP1311-P1-ENTHALPY. Gordon and McBride, NASA RP-1311 Part I (1994),
       section 4.1, equations (4.1)-(4.2), printed p. 19.
       https://ntrs.nasa.gov/citations/19950013764
.. [5] NASA-TP1999-208522. Jankovsky, Smith and Pavli, NASA/TP-1999-208522,
       equations (3)-(11), printed pp. 12-13; loss chain printed p. 25 and
       Figures 20 and 24. https://ntrs.nasa.gov/citations/19990051032
.. [6] NACA-TN4219. Heidmann and Priem, NACA TN 4219 (1958), equations
       (1)-(7), printed pp. 3-4. https://ntrs.nasa.gov/citations/19930085136
.. [7] NASA-CR187109-ROCCID. Muss, Nguyen and Johnson, NASA CR-187109 (1991),
       section 2.2, equations (2.4)-(2.11), printed pp. 16-18.
       https://ntrs.nasa.gov/citations/19910014917
.. [8] NASA-SP8089-MIXING. NASA SP-8089 (1976), equation (1) and Table IV,
       printed pp. 19-20. https://ntrs.nasa.gov/citations/19760023196
.. [9] NASA-SP125-LSTAR. Huzel and Huang, NASA SP-125 (1967), equations (4-3)
       and (4-4) and Table 4-1, printed pp. 86-87.
       https://ntrs.nasa.gov/citations/19710019929
.. [10] NASA-TM89814-JANNAF. Davidian, NASA TM-89814 (1987), printed pp. 1-4.
        https://ntrs.nasa.gov/citations/19870010946
"""

from __future__ import annotations

from collections.abc import Callable, Mapping, Sequence
from dataclasses import dataclass
from math import isfinite, sqrt
from typing import Literal

from rocket_propulsion.core.constants import STANDARD_GRAVITY_M_S2
from rocket_propulsion.core.errors import DomainError, InputError
from rocket_propulsion.core.metadata import WarningMessage

from .equilibrium import (
    EquilibriumProblem,
    EquilibriumProvider,
    EquilibriumResult,
    LocalEquilibriumProvider,
)
from .rocket_equilibrium import PropellantPair, get_propellant_pair, reactants_for_mass_ratio
from .rocket_performance import Chemistry, solve_rocket_performance
from .species import get_reactant, get_species

#: Oxidation states used by NASA CEA for the chemical equivalence ratio
#: (``source/atomic_data.f90``, commit 2f79a647, entries H, C, N, O).
CEA_OXIDATION_STATES: dict[str, float] = {"H": 1.0, "C": 4.0, "N": 0.0, "O": -2.0}

#: Open-shell species of the bundled catalog reported as radicals. This list is
#: a reporting choice of this module, not a quantity from a cited source.
DEFAULT_RADICALS: tuple[str, ...] = (
    "H", "O", "OH", "HO2", "N", "NH", "NH2", "NO", "NO2",
    "C", "CH", "CH2", "CH3", "HCO", "CN", "NCO", "COOH",
)

#: Rupe-mixing-factor ``M`` of SP-8089 Table IV (printed p. 20) by element type.
SP8089_MIXING_FACTOR: dict[str, float] = {
    "1-on-1": 1.0, "2-on-2": 1.0, "2-on-1": 1.6, "3-on-1": 3.5, "4-on-1": 9.4, "5-on-1": 27.5,
}

_INCH_M = 0.0254
#: Recommended chamber characteristic length ``L*`` ranges of SP-125 Table 4-1
#: (printed p. 87), converted from inches to metres.
SP125_RECOMMENDED_CHARACTERISTIC_LENGTH_M: dict[str, tuple[float, float]] = {
    "LOX_LH2_GH2_INJECTION": (22.0 * _INCH_M, 28.0 * _INCH_M),
    "LOX_LH2_LH2_INJECTION": (30.0 * _INCH_M, 40.0 * _INCH_M),
    "LOX_RP1": (40.0 * _INCH_M, 50.0 * _INCH_M),
    "NTO_HYDRAZINE_BASE": (30.0 * _INCH_M, 35.0 * _INCH_M),
}

Objective = Literal["vacuum_specific_impulse", "density_vacuum_impulse"]
_GOLDEN = (sqrt(5.0) - 1.0) / 2.0
#: Shifting then frozen chemistry, the two limits of :func:`frozen_equilibrium_limits`.
_CHEMISTRY_LIMITS: tuple[Chemistry, Chemistry] = ("equilibrium", "frozen")


def _positive(value: float, field: str) -> float:
    number = float(value)
    if not isfinite(number) or number <= 0.0:
        raise DomainError(f"{field} must be finite and positive.", field=field)
    return number


def _fraction(value: float, field: str, *, allow_zero: bool = True) -> float:
    number = float(value)
    low_ok = number >= 0.0 if allow_zero else number > 0.0
    if not isfinite(number) or not low_ok or number > 1.0:
        raise DomainError(f"{field} must lie in {'[0' if allow_zero else '(0'}, 1].", field=field)
    return number


# ---------------------------------------------------------------------------
# Stoichiometry
# ---------------------------------------------------------------------------


@dataclass(frozen=True, slots=True)
class ValenceTotals:
    """Positive and negative oxidation states per kilogram of one reactant [1]_.

    Attributes
    ----------
    reactant : str
        Reactant key.
    positive_per_kg, negative_per_kg : float
        ``V+ = sum_i V_i+ b_i`` and ``V- = sum_i V_i- b_i`` in kg-atom
        valence per kg (RP-1311 eqs. (9.14)-(9.15)).

    References
    ----------
    .. [1] NASA-RP1311-P1-EQRATIO. Gordon and McBride, NASA RP-1311 Part I
           (1994), chapter 9, equations (9.10)-(9.19), printed pp. 42-43.
           https://ntrs.nasa.gov/citations/19950013764
    """

    reactant: str
    positive_per_kg: float
    negative_per_kg: float


def reactant_valence(reactant: str) -> ValenceTotals:
    """Return the CEA oxidation-state totals of a reactant per kilogram [1]_, [2]_.

    Parameters
    ----------
    reactant : str
        Gas species or condensed reactant key.

    Returns
    -------
    ValenceTotals
        ``V+`` and ``V-`` per kg.

    Raises
    ------
    InputError
        If the reactant is unknown or contains an element without a CEA
        valence in :data:`CEA_OXIDATION_STATES`.

    References
    ----------
    .. [1] NASA-RP1311-P1-EQRATIO. Gordon and McBride, NASA RP-1311 Part I
           (1994), chapter 9, equations (9.10)-(9.19), printed pp. 42-43.
           https://ntrs.nasa.gov/citations/19950013764
    .. [2] NASA-CEA-SRC-VALENCE. NASA Glenn, CEA source
           ``source/atomic_data.f90`` and ``source/mixture.f90``, commit
           2f79a647. https://github.com/nasa/cea
    """

    record = get_reactant(reactant)
    positive = negative = 0.0
    for element, count in record.formula:
        if element not in CEA_OXIDATION_STATES:
            raise InputError(f"No CEA valence registered for element {element!r}.", field="reactant",
                             code="unknown_valence")
        valence = CEA_OXIDATION_STATES[element] * count / record.molecular_mass_kg_mol
        if valence > 0.0:
            positive += valence
        else:
            negative += valence
    return ValenceTotals(record.key, positive, negative)


def stoichiometric_oxidizer_fuel_ratio(fuel: str, oxidizer: str) -> float:
    """Return the O/F mass ratio at chemical equivalence ratio ``r = 1`` [1]_, [2]_.

    Parameters
    ----------
    fuel, oxidizer : str
        Reactant keys.

    Returns
    -------
    float
        ``-(V+_f + V-_f) / (V+_o + V-_o)``, RP-1311 eq. (9.18) with ``r = 1``.

    Raises
    ------
    DomainError
        If the oxidizer carries no net negative valence or the fuel no net
        positive valence.

    References
    ----------
    .. [1] NASA-RP1311-P1-EQRATIO. Gordon and McBride, NASA RP-1311 Part I
           (1994), chapter 9, equations (9.10)-(9.19), printed pp. 42-43.
           https://ntrs.nasa.gov/citations/19950013764
    .. [2] NASA-CEA-SRC-VALENCE. NASA Glenn, CEA source
           ``source/atomic_data.f90`` and ``source/mixture.f90``, commit
           2f79a647. https://github.com/nasa/cea
    """

    return oxidizer_fuel_ratio_for_chemical_equivalence_ratio(fuel, oxidizer, 1.0)


def oxidizer_fuel_ratio_for_chemical_equivalence_ratio(fuel: str, oxidizer: str, r: float) -> float:
    """Return O/F for an assigned chemical equivalence ratio ``r`` [1]_, [2]_.

    Parameters
    ----------
    fuel, oxidizer : str
        Reactant keys.
    r : float
        Chemical equivalence ratio, positive.

    Returns
    -------
    float
        ``-(V+_f + r V-_f) / (V+_o + r V-_o)`` (CEA
        ``mixture_chem_eq_ratio_to_of_ratio``).

    Raises
    ------
    DomainError
        If ``r`` is not positive or the result is not a positive ratio.

    References
    ----------
    .. [1] NASA-RP1311-P1-EQRATIO. Gordon and McBride, NASA RP-1311 Part I
           (1994), chapter 9, equations (9.10)-(9.19), printed pp. 42-43.
           https://ntrs.nasa.gov/citations/19950013764
    .. [2] NASA-CEA-SRC-VALENCE. NASA Glenn, CEA source
           ``source/atomic_data.f90`` and ``source/mixture.f90``, commit
           2f79a647. https://github.com/nasa/cea
    """

    r = _positive(r, "r")
    f, o = reactant_valence(fuel), reactant_valence(oxidizer)
    denominator = o.positive_per_kg + r * o.negative_per_kg
    ratio = -(f.positive_per_kg + r * f.negative_per_kg) / denominator if denominator else float("nan")
    if not isfinite(ratio) or ratio <= 0.0:
        raise DomainError("Reactants do not form a fuel/oxidizer pair at this equivalence ratio.",
                          code="invalid_equivalence_ratio", field="r")
    return ratio


def chemical_equivalence_ratio(fuel: str, oxidizer: str, oxidizer_fuel_ratio: float) -> float:
    """Return the chemical equivalence ratio ``r = -V+/V-`` [1]_, [2]_.

    Parameters
    ----------
    fuel, oxidizer : str
        Reactant keys.
    oxidizer_fuel_ratio : float
        O/F mass ratio, positive.

    Returns
    -------
    float
        ``-(V+_f + (o/f) V+_o) / (V-_f + (o/f) V-_o)``, RP-1311 eqs.
        (9.16)-(9.18); greater than one is fuel rich.

    Raises
    ------
    DomainError
        If the ratio is not positive or the mixture has no negative valence.

    References
    ----------
    .. [1] NASA-RP1311-P1-EQRATIO. Gordon and McBride, NASA RP-1311 Part I
           (1994), chapter 9, equations (9.10)-(9.19), printed pp. 42-43.
           https://ntrs.nasa.gov/citations/19950013764
    .. [2] NASA-CEA-SRC-VALENCE. NASA Glenn, CEA source
           ``source/atomic_data.f90`` and ``source/mixture.f90``, commit
           2f79a647. https://github.com/nasa/cea
    """

    ratio = _positive(oxidizer_fuel_ratio, "oxidizer_fuel_ratio")
    f, o = reactant_valence(fuel), reactant_valence(oxidizer)
    negative = f.negative_per_kg + ratio * o.negative_per_kg
    if negative >= 0.0:
        raise DomainError("Mixture carries no negative valence.", field="oxidizer")
    return -(f.positive_per_kg + ratio * o.positive_per_kg) / negative


def equivalence_ratio(fuel: str, oxidizer: str, oxidizer_fuel_ratio: float) -> float:
    """Return the weight equivalence ratio ``phi = (f/o) / (f/o)_stoich`` [1]_.

    Parameters
    ----------
    fuel, oxidizer : str
        Reactant keys.
    oxidizer_fuel_ratio : float
        O/F mass ratio, positive.

    Returns
    -------
    float
        ``(o/f)_stoich / (o/f)``, RP-1311 eq. (9.19); greater than one is
        fuel rich.

    References
    ----------
    .. [1] NASA-RP1311-P1-EQRATIO. Gordon and McBride, NASA RP-1311 Part I
           (1994), chapter 9, equations (9.10)-(9.19), printed pp. 42-43.
           https://ntrs.nasa.gov/citations/19950013764
    """

    ratio = _positive(oxidizer_fuel_ratio, "oxidizer_fuel_ratio")
    return stoichiometric_oxidizer_fuel_ratio(fuel, oxidizer) / ratio


def bulk_density(oxidizer_fuel_ratio: float, oxidizer_density_kg_m3: float,
                 fuel_density_kg_m3: float) -> float:
    """Return the total-reactant density of RP-1311 eq. (9.13) [1]_.

    Parameters
    ----------
    oxidizer_fuel_ratio : float
        O/F mass ratio.
    oxidizer_density_kg_m3, fuel_density_kg_m3 : float
        Densities of the total oxidant and total fuel.

    Returns
    -------
    float
        ``(o/f + 1) rho_o rho_f / (rho_o + (o/f) rho_f)``.

    Notes
    -----
    SP-125 equation (1-46) prints ``(r_w - 1)`` in the numerator; that
    disagrees with RP-1311 (9.13) and with mass conservation, so the RP-1311
    form is used.

    References
    ----------
    .. [1] NASA-RP1311-P1-EQRATIO. Gordon and McBride, NASA RP-1311 Part I
           (1994), chapter 9, equations (9.10)-(9.19), printed pp. 42-43.
           https://ntrs.nasa.gov/citations/19950013764
    """

    ratio = _positive(oxidizer_fuel_ratio, "oxidizer_fuel_ratio")
    rho_o = _positive(oxidizer_density_kg_m3, "oxidizer_density_kg_m3")
    rho_f = _positive(fuel_density_kg_m3, "fuel_density_kg_m3")
    return (ratio + 1.0) * rho_o * rho_f / (rho_o + ratio * rho_f)


def density_impulse(specific_impulse_s: float, density: float) -> float:
    """Return density impulse ``I_s d`` of SP-125 eq. (1-45) [1]_.

    Parameters
    ----------
    specific_impulse_s : float
        Specific impulse in s.
    density : float
        Bulk propellant density in any unit; the product carries that unit
        times seconds.

    Returns
    -------
    float
        ``I_s d``.

    References
    ----------
    .. [1] NASA-SP125-DENSITY-IMPULSE. Huzel and Huang, NASA SP-125 (1967),
           equation (1-45), printed p. 20.
           https://ntrs.nasa.gov/citations/19710019929
    """

    return _positive(specific_impulse_s, "specific_impulse_s") * _positive(density, "density")


@dataclass(frozen=True, slots=True)
class StoichiometryReport:
    """Mixture-ratio descriptors of one fuel/oxidizer pair at one O/F [1]_, [2]_.

    Attributes
    ----------
    fuel, oxidizer : str
        Reactant keys.
    oxidizer_fuel_ratio : float
        O/F mass ratio.
    stoichiometric_oxidizer_fuel_ratio : float
        O/F at ``r = 1``.
    chemical_equivalence_ratio : float
        ``r`` of RP-1311 eq. (9.18).
    equivalence_ratio : float
        ``phi`` of RP-1311 eq. (9.19).
    fuel_mass_fraction : float
        ``1 / (1 + o/f)``.

    References
    ----------
    .. [1] NASA-RP1311-P1-EQRATIO. Gordon and McBride, NASA RP-1311 Part I
           (1994), chapter 9, equations (9.10)-(9.19), printed pp. 42-43.
           https://ntrs.nasa.gov/citations/19950013764
    .. [2] NASA-CEA-SRC-VALENCE. NASA Glenn, CEA source
           ``source/atomic_data.f90`` and ``source/mixture.f90``, commit
           2f79a647. https://github.com/nasa/cea
    """

    fuel: str
    oxidizer: str
    oxidizer_fuel_ratio: float
    stoichiometric_oxidizer_fuel_ratio: float
    chemical_equivalence_ratio: float
    equivalence_ratio: float
    fuel_mass_fraction: float


def stoichiometry_report(fuel: str, oxidizer: str, oxidizer_fuel_ratio: float) -> StoichiometryReport:
    """Collect stoichiometric O/F, ``r`` and ``phi`` for one mixture [1]_, [2]_.

    Parameters
    ----------
    fuel, oxidizer : str
        Reactant keys.
    oxidizer_fuel_ratio : float
        O/F mass ratio.

    Returns
    -------
    StoichiometryReport
        Mixture descriptors.

    References
    ----------
    .. [1] NASA-RP1311-P1-EQRATIO. Gordon and McBride, NASA RP-1311 Part I
           (1994), chapter 9, equations (9.10)-(9.19), printed pp. 42-43.
           https://ntrs.nasa.gov/citations/19950013764
    .. [2] NASA-CEA-SRC-VALENCE. NASA Glenn, CEA source
           ``source/atomic_data.f90`` and ``source/mixture.f90``, commit
           2f79a647. https://github.com/nasa/cea
    """

    ratio = _positive(oxidizer_fuel_ratio, "oxidizer_fuel_ratio")
    stoich = stoichiometric_oxidizer_fuel_ratio(fuel, oxidizer)
    return StoichiometryReport(fuel, oxidizer, ratio, stoich, chemical_equivalence_ratio(fuel, oxidizer, ratio),
                               stoich / ratio, 1.0 / (1.0 + ratio))


# ---------------------------------------------------------------------------
# Product composition and dissociation
# ---------------------------------------------------------------------------


@dataclass(frozen=True, slots=True)
class StationComposition:
    """Full equilibrium composition and energy split at one station [1]_, [2]_.

    Attributes
    ----------
    name : str
        ``"chamber"``, ``"throat"``, ``"exit"`` or ``"reference"``.
    temperature_k, pressure_pa : float
        State.
    molecular_mass_kg_mol : float
        Mixture molar mass.
    mole_fractions, mass_fractions : tuple of (str, float)
        Every retained species (local provider keeps amounts above 1e-14 mol),
        largest first.
    radicals : tuple of (str, float)
        Mole fractions of the species in the report's radical list.
    radical_mole_fraction : float
        Sum of ``radicals``.
    enthalpy_j_kg : float
        Mixture specific enthalpy (RP-1311 assigned-enthalpy scale).
    formation_enthalpy_j_kg : float
        ``sum_j (n_j/m) Delta_f H_j(298.15)``, the chemical part of eq. (4.2).
    sensible_enthalpy_j_kg : float
        ``h - formation_enthalpy_j_kg``.
    element_relative_residual : float
        Largest relative difference between product and reactant elements.
    mass_fraction_sum_residual : float
        ``|sum(w_j) - 1|``.

    References
    ----------
    .. [1] NASA-RP1311-P1-ENTHALPY. Gordon and McBride, NASA RP-1311 Part I
           (1994), section 4.1, equations (4.1)-(4.2), printed p. 19.
           https://ntrs.nasa.gov/citations/19950013764
    .. [2] NASA-RP1311-P1-FAC. Gordon and McBride, NASA RP-1311 Part I (1994),
           chapter 6, sections 6.1-6.4, equations (6.1)-(6.33).
           https://ntrs.nasa.gov/citations/19950013764
    """

    name: str
    temperature_k: float
    pressure_pa: float
    molecular_mass_kg_mol: float
    mole_fractions: tuple[tuple[str, float], ...]
    mass_fractions: tuple[tuple[str, float], ...]
    radicals: tuple[tuple[str, float], ...]
    radical_mole_fraction: float
    enthalpy_j_kg: float
    formation_enthalpy_j_kg: float
    sensible_enthalpy_j_kg: float
    element_relative_residual: float
    mass_fraction_sum_residual: float


@dataclass(frozen=True, slots=True)
class DissociationSummary:
    """Chemical energy held in dissociation and returned by recombination [1]_.

    Attributes
    ----------
    reference_temperature_k : float
        Temperature of the recombined reference state (equilibrium at the
        chamber pressure).
    available_energy_j_kg : float
        ``h_c - h_ref``: energy released from reactants to recombined
        products at the reference temperature.
    dissociation_energy_j_kg : float
        ``h_f,c - h_f,ref``: chemical energy still held in the chamber
        composition.
    dissociation_fraction : float
        ``dissociation_energy / available_energy``.
    sensible_fraction : float
        ``(h_s,c - h_s,ref) / available_energy``; with
        ``dissociation_fraction`` it sums to one by construction.
    throat_recombination_j_kg, exit_recombination_j_kg : float
        ``h_f,c - h_f,station``: chemical energy returned to the flow by
        shifting-equilibrium recombination up to the station.
    exit_recovered_fraction : float
        ``exit_recombination / dissociation_energy``.

    Notes
    -----
    The split uses only the assigned-enthalpy decomposition of RP-1311
    eq. (4.2) [1]_; the energy names are definitions of this module. The
    reference composition is gas-phase equilibrium from the bundled species
    (water stays vapour; no graphite).

    References
    ----------
    .. [1] NASA-RP1311-P1-ENTHALPY. Gordon and McBride, NASA RP-1311 Part I
           (1994), section 4.1, equations (4.1)-(4.2), printed p. 19.
           https://ntrs.nasa.gov/citations/19950013764
    """

    reference_temperature_k: float
    available_energy_j_kg: float
    dissociation_energy_j_kg: float
    dissociation_fraction: float
    sensible_fraction: float
    throat_recombination_j_kg: float
    exit_recombination_j_kg: float
    exit_recovered_fraction: float


@dataclass(frozen=True, slots=True)
class CombustionProductReport:
    """Composition, stoichiometry and dissociation report of one operating point [1]_.

    Attributes
    ----------
    pair : str or None
        Propellant preset key when given.
    reactants : tuple of (str, float)
        Reactant amounts in mol.
    stoichiometry : StoichiometryReport or None
        Mixture descriptors when the pair is known.
    chamber_pressure_pa, area_ratio : float
        Operating point.
    characteristic_velocity_m_s, vacuum_specific_impulse_s : float
        Shifting-equilibrium IAC performance of the same solve.
    stations : tuple of StationComposition
        Chamber, throat and exit.
    reference : StationComposition
        Recombined reference state.
    dissociation : DissociationSummary
        Energy split.
    element_moles : tuple of (str, float)
        Reactant element inventory.
    warnings : tuple of WarningMessage
        Model statements.

    References
    ----------
    .. [1] NASA-RP1311-P1-FAC. Gordon and McBride, NASA RP-1311 Part I (1994),
           chapter 6, sections 6.1-6.4, equations (6.1)-(6.33).
           https://ntrs.nasa.gov/citations/19950013764
    """

    pair: str | None
    reactants: tuple[tuple[str, float], ...]
    stoichiometry: StoichiometryReport | None
    chamber_pressure_pa: float
    area_ratio: float
    characteristic_velocity_m_s: float
    vacuum_specific_impulse_s: float
    stations: tuple[StationComposition, ...]
    reference: StationComposition
    dissociation: DissociationSummary
    element_moles: tuple[tuple[str, float], ...]
    warnings: tuple[WarningMessage, ...]


def _station_composition(name: str, result: EquilibriumResult, radicals: Sequence[str]) -> StationComposition:
    mass = result.total_mass_kg
    formation = sum(item.moles * get_species(item.species).heat_of_formation_j_mol_298
                    for item in result.species) / mass
    enthalpy = result.enthalpy_j / mass
    inventory = dict(result.element_moles)
    products: dict[str, float] = {}
    for item in result.species:
        for element, count in get_species(item.species).formula:
            products[element] = products.get(element, 0.0) + item.moles * count
    element_residual = max(abs(products.get(key, 0.0) - value) / value for key, value in inventory.items())
    fractions = {item.species: item.mole_fraction for item in result.species}
    selected = tuple((key, fractions[key]) for key in radicals if key in fractions)
    return StationComposition(
        name=name,
        temperature_k=result.temperature_k,
        pressure_pa=result.pressure_pa,
        molecular_mass_kg_mol=result.molecular_mass_kg_mol,
        mole_fractions=tuple((item.species, item.mole_fraction) for item in result.species),
        mass_fractions=tuple(sorted(((item.species, item.mass_fraction) for item in result.species),
                                    key=lambda pair: -pair[1])),
        radicals=selected,
        radical_mole_fraction=sum(value for _, value in selected),
        enthalpy_j_kg=enthalpy,
        formation_enthalpy_j_kg=formation,
        sensible_enthalpy_j_kg=enthalpy - formation,
        element_relative_residual=element_residual,
        mass_fraction_sum_residual=abs(sum(item.mass_fraction for item in result.species) - 1.0),
    )


def combustion_product_report(
    *,
    chamber_pressure_pa: float,
    area_ratio: float,
    propellant_pair: str | None = None,
    oxidizer_fuel_ratio: float | None = None,
    reactants: Mapping[str, float] | None = None,
    radicals: Sequence[str] = DEFAULT_RADICALS,
    reference_temperature_k: float = 298.15,
    provider: EquilibriumProvider | None = None,
) -> CombustionProductReport:
    """Report equilibrium products, radicals and dissociation energy [1]_, [2]_, [3]_.

    Parameters
    ----------
    chamber_pressure_pa : float
        Chamber pressure in Pa.
    area_ratio : float
        Supersonic exit-to-throat area ratio, greater than one.
    propellant_pair : str, optional
        Preset key (with ``oxidizer_fuel_ratio``); alternatively give
        ``reactants``.
    oxidizer_fuel_ratio : float, optional
        O/F mass ratio for the preset.
    reactants : mapping of str to float, optional
        Reactant amounts in mol when no preset is used.
    radicals : sequence of str, default :data:`DEFAULT_RADICALS`
        Species reported as the radical inventory.
    reference_temperature_k : float, default 298.15
        Temperature of the recombined reference state, 200 to 1000 K.
    provider : EquilibriumProvider, optional
        Equilibrium provider; the local provider by default.

    Returns
    -------
    CombustionProductReport
        Station compositions, radical inventory, element and mass balances,
        stoichiometry and dissociation energy split.

    Raises
    ------
    InputError
        If neither or both of the preset and the reactants are given.
    DomainError
        If a state or the reference temperature is out of range.

    Notes
    -----
    Stations come from
    :func:`~rocket_propulsion.thermochemistry.rocket_performance.solve_rocket_performance`
    with shifting equilibrium (the CEA infinite-area-combustor problem [2]_);
    each station is re-solved at its temperature and pressure to keep the
    full species list. The enthalpy split follows RP-1311 eq. (4.2) [1]_.
    Stoichiometry uses CEA valences [3]_.

    References
    ----------
    .. [1] NASA-RP1311-P1-ENTHALPY. Gordon and McBride, NASA RP-1311 Part I
           (1994), section 4.1, equations (4.1)-(4.2), printed p. 19.
           https://ntrs.nasa.gov/citations/19950013764
    .. [2] NASA-RP1311-P1-FAC. Gordon and McBride, NASA RP-1311 Part I (1994),
           chapter 6, sections 6.1-6.4, equations (6.1)-(6.33).
           https://ntrs.nasa.gov/citations/19950013764
    .. [3] NASA-RP1311-P1-EQRATIO. Gordon and McBride, NASA RP-1311 Part I
           (1994), chapter 9, equations (9.10)-(9.19), printed pp. 42-43.
           https://ntrs.nasa.gov/citations/19950013764
    """

    pair: PropellantPair | None = None
    if (propellant_pair is None) == (reactants is None):
        raise InputError("Give either propellant_pair with oxidizer_fuel_ratio or reactants.",
                         field="propellant_pair")
    if propellant_pair is not None:
        if oxidizer_fuel_ratio is None:
            raise InputError("A propellant pair needs oxidizer_fuel_ratio.", field="oxidizer_fuel_ratio")
        pair = get_propellant_pair(propellant_pair)
        mixture = reactants_for_mass_ratio(pair, oxidizer_fuel_ratio)
    else:
        mixture = dict(reactants)  # type: ignore[arg-type]
    if not isfinite(reference_temperature_k) or not 200.0 <= reference_temperature_k <= 1000.0:
        raise DomainError("Reference temperature must be within 200-1000 K.", field="reference_temperature_k")
    active = provider or LocalEquilibriumProvider()
    performance = solve_rocket_performance(reactants=mixture, chamber_pressure_pa=chamber_pressure_pa,
                                           area_ratios=(area_ratio,), chemistry="equilibrium", provider=active)
    states = (("chamber", performance.chamber), ("throat", performance.throat),
              ("exit", performance.exits[0].station))
    stations: list[StationComposition] = []
    results: list[EquilibriumResult] = []
    for name, station in states:
        result = active.solve(EquilibriumProblem.tp(mixture, temperature_k=station.temperature_k,
                                                    pressure_pa=station.pressure_pa))
        results.append(result)
        stations.append(_station_composition(name, result, radicals))
    reference_result = active.solve(EquilibriumProblem.tp(mixture, temperature_k=reference_temperature_k,
                                                          pressure_pa=chamber_pressure_pa))
    reference = _station_composition("reference", reference_result, radicals)
    chamber, throat, exit_ = stations
    available = chamber.enthalpy_j_kg - reference.enthalpy_j_kg
    dissociation = chamber.formation_enthalpy_j_kg - reference.formation_enthalpy_j_kg
    summary = DissociationSummary(
        reference_temperature_k=reference_temperature_k,
        available_energy_j_kg=available,
        dissociation_energy_j_kg=dissociation,
        dissociation_fraction=dissociation / available,
        sensible_fraction=(chamber.sensible_enthalpy_j_kg - reference.sensible_enthalpy_j_kg) / available,
        throat_recombination_j_kg=chamber.formation_enthalpy_j_kg - throat.formation_enthalpy_j_kg,
        exit_recombination_j_kg=chamber.formation_enthalpy_j_kg - exit_.formation_enthalpy_j_kg,
        exit_recovered_fraction=(chamber.formation_enthalpy_j_kg - exit_.formation_enthalpy_j_kg) / dissociation,
    )
    stoich = None
    if pair is not None and oxidizer_fuel_ratio is not None:
        stoich = stoichiometry_report(pair.fuel_reactant, pair.oxidizer_reactant, float(oxidizer_fuel_ratio))
    return CombustionProductReport(
        pair=pair.key if pair is not None else None,
        reactants=performance.reactants,
        stoichiometry=stoich,
        chamber_pressure_pa=chamber_pressure_pa,
        area_ratio=float(area_ratio),
        characteristic_velocity_m_s=performance.characteristic_velocity_m_s,
        vacuum_specific_impulse_s=performance.exits[0].vacuum_specific_impulse_s,
        stations=tuple(stations),
        reference=reference,
        dissociation=summary,
        element_moles=results[0].element_moles,
        warnings=(WarningMessage(
            code="gas_phase_products",
            message=("Gas-phase equilibrium from the bundled species set; no condensed products "
                     "(the recombined reference keeps water as vapour)."),
            severity="info"),),
    )


# ---------------------------------------------------------------------------
# Efficiencies from measurements and loss models
# ---------------------------------------------------------------------------


def characteristic_velocity(chamber_pressure_pa: float, throat_area_m2: float, mass_flow_kg_s: float) -> float:
    """Return ``c* = p_c A_t / mdot`` [1]_, [2]_.

    Parameters
    ----------
    chamber_pressure_pa, throat_area_m2, mass_flow_kg_s : float
        Positive chamber (stagnation) pressure, throat area and mass flow.

    Returns
    -------
    float
        Characteristic velocity in m/s.

    References
    ----------
    .. [1] NASA-TP1999-208522. Jankovsky, Smith and Pavli, NASA/TP-1999-208522,
           equations (3)-(11), printed pp. 12-13; loss chain printed p. 25 and
           Figures 20 and 24. https://ntrs.nasa.gov/citations/19990051032
    .. [2] NASA-RP1311-P1-FAC. Gordon and McBride, NASA RP-1311 Part I (1994),
           chapter 6, sections 6.1-6.4, equations (6.1)-(6.33).
           https://ntrs.nasa.gov/citations/19950013764
    """

    return (_positive(chamber_pressure_pa, "chamber_pressure_pa") * _positive(throat_area_m2, "throat_area_m2")
            / _positive(mass_flow_kg_s, "mass_flow_kg_s"))


@dataclass(frozen=True, slots=True)
class MeasuredEfficiencies:
    """Measured-to-ODE efficiencies of NASA/TP-1999-208522 eqs. (6)-(11) [1]_.

    Attributes
    ----------
    characteristic_velocity_m_s : float
        ``p A_t / mdot``.
    vacuum_thrust_coefficient : float
        ``F_v / (p A_t)``.
    vacuum_specific_impulse_s : float
        ``F_v / (mdot g0)``.
    characteristic_velocity_efficiency, thrust_coefficient_efficiency : float
        Ratios to the ODE values.
    specific_impulse_efficiency : float
        Ratio to the ODE value; equals the product of the other two when the
        ODE values satisfy ``I g0 = c* C_F``.

    References
    ----------
    .. [1] NASA-TP1999-208522. Jankovsky, Smith and Pavli, NASA/TP-1999-208522,
           equations (3)-(11), printed pp. 12-13; loss chain printed p. 25 and
           Figures 20 and 24. https://ntrs.nasa.gov/citations/19990051032
    """

    characteristic_velocity_m_s: float
    vacuum_thrust_coefficient: float
    vacuum_specific_impulse_s: float
    characteristic_velocity_efficiency: float
    thrust_coefficient_efficiency: float
    specific_impulse_efficiency: float


def measured_efficiencies(*, chamber_pressure_pa: float, throat_area_m2: float, mass_flow_kg_s: float,
                          vacuum_thrust_n: float, ideal_characteristic_velocity_m_s: float,
                          ideal_vacuum_specific_impulse_s: float) -> MeasuredEfficiencies:
    """Reduce test-stand data to c*, C_F and Isp efficiencies [1]_.

    Parameters
    ----------
    chamber_pressure_pa : float
        Effective chamber total pressure ``P_c,e`` (TP-1999-208522 eq. (3)).
    throat_area_m2, mass_flow_kg_s, vacuum_thrust_n : float
        Measured throat area, mass flow and vacuum thrust.
    ideal_characteristic_velocity_m_s, ideal_vacuum_specific_impulse_s : float
        ODE values at the same O/F and pressure.

    Returns
    -------
    MeasuredEfficiencies
        Equations (6)-(11) of [1]_; the ODE thrust coefficient is
        ``I_vac g0 / c*``.

    References
    ----------
    .. [1] NASA-TP1999-208522. Jankovsky, Smith and Pavli, NASA/TP-1999-208522,
           equations (3)-(11), printed pp. 12-13; loss chain printed p. 25 and
           Figures 20 and 24. https://ntrs.nasa.gov/citations/19990051032
    """

    cstar = characteristic_velocity(chamber_pressure_pa, throat_area_m2, mass_flow_kg_s)
    thrust = _positive(vacuum_thrust_n, "vacuum_thrust_n")
    coefficient = thrust / (chamber_pressure_pa * throat_area_m2)
    isp = thrust / (mass_flow_kg_s * STANDARD_GRAVITY_M_S2)
    ideal_cstar = _positive(ideal_characteristic_velocity_m_s, "ideal_characteristic_velocity_m_s")
    ideal_isp = _positive(ideal_vacuum_specific_impulse_s, "ideal_vacuum_specific_impulse_s")
    ideal_coefficient = ideal_isp * STANDARD_GRAVITY_M_S2 / ideal_cstar
    return MeasuredEfficiencies(cstar, coefficient, isp, cstar / ideal_cstar, coefficient / ideal_coefficient,
                                isp / ideal_isp)


def vaporization_limited_cstar_efficiency(
    *,
    characteristic_velocity: Callable[[float], float],
    oxidizer_fuel_ratio: float,
    oxidizer_vaporized_fraction: float,
    fuel_vaporized_fraction: float,
) -> float:
    """Return the c* efficiency of incompletely vaporized propellants [1]_.

    Parameters
    ----------
    characteristic_velocity : callable
        Theoretical ``c*(O/F)`` in m/s, for example
        :func:`cea_characteristic_velocity_function`.
    oxidizer_fuel_ratio : float
        Injected O/F ``w_o / w_f``.
    oxidizer_vaporized_fraction, fuel_vaporized_fraction : float
        Fractions ``O`` and ``F`` vaporized and burned, each in (0, 1].

    Returns
    -------
    float
        ``eta_c* = [c*(O w_o / F w_f) / c*(w_o / w_f)] (O w_o + F w_f)/(w_o + w_f)``,
        NACA TN 4219 eq. (5); eqs. (6) and (7) are the cases ``F = 1`` and
        ``O = 1``.

    Notes
    -----
    The model of [1]_ treats the liquid drops as carrying no volume or kinetic
    energy and the vapour as burning to equilibrium at the vapour mixture
    ratio; only the vaporized mass contributes to ``p_c A_t``.

    References
    ----------
    .. [1] NACA-TN4219. Heidmann and Priem, NACA TN 4219 (1958), equations
           (1)-(7), printed pp. 3-4. https://ntrs.nasa.gov/citations/19930085136
    """

    ratio = _positive(oxidizer_fuel_ratio, "oxidizer_fuel_ratio")
    oxidant = _fraction(oxidizer_vaporized_fraction, "oxidizer_vaporized_fraction", allow_zero=False)
    fuel = _fraction(fuel_vaporized_fraction, "fuel_vaporized_fraction", allow_zero=False)
    vapour_ratio = oxidant * ratio / fuel
    burned = (oxidant * ratio + fuel) / (ratio + 1.0)
    return characteristic_velocity(vapour_ratio) / characteristic_velocity(ratio) * burned


@dataclass(frozen=True, slots=True)
class StreamtubeSplit:
    """Two-streamtube mixture-ratio distribution of one zone [1]_.

    Attributes
    ----------
    rupe_mixing_efficiency : float
        ``E_m`` as a fraction.
    mixing_parameter : float
        ``eta_m = 0.449 - 0.399 E_m + 0.944 E_m^2`` (eq. (2.5)).
    mean_mixture_ratio : float
        Zone mean O/F ``MR_M``.
    low_mixture_ratio, high_mixture_ratio : float
        ``MR_R = MR_M eta_m`` and ``MR_L = MR_M / eta_m``.
    low_mass_fraction, high_mass_fraction : float
        ``X_R = (1 + MR_R) / ((1 + eta_m)(1 + MR_M))`` and ``X_L = 1 - X_R``.

    Notes
    -----
    ROCCID quotes ``E_m`` as a percentage in its element-density correlation
    (eq. (2.6)); eq. (2.5) only yields a mixture-ratio multiplier between 0
    and 1 with ``E_m`` as a fraction, which is how it is evaluated here. At
    ``E_m = 1`` the fit gives ``eta_m = 0.994``, not exactly 1.

    References
    ----------
    .. [1] NASA-CR187109-ROCCID. Muss, Nguyen and Johnson, NASA CR-187109
           (1991), section 2.2, equations (2.4)-(2.11), printed pp. 16-18.
           https://ntrs.nasa.gov/citations/19910014917
    """

    rupe_mixing_efficiency: float
    mixing_parameter: float
    mean_mixture_ratio: float
    low_mixture_ratio: float
    high_mixture_ratio: float
    low_mass_fraction: float
    high_mass_fraction: float


def rupe_streamtube_split(rupe_mixing_efficiency: float, mean_mixture_ratio: float) -> StreamtubeSplit:
    """Split a zone into two streamtubes from its Rupe ``E_m`` [1]_.

    Parameters
    ----------
    rupe_mixing_efficiency : float
        ``E_m`` as a fraction in (0.25, 1]; below about 0.25 the fit of
        eq. (2.5) stops decreasing and the split is not monotone.
    mean_mixture_ratio : float
        Zone mean O/F.

    Returns
    -------
    StreamtubeSplit
        Equations (2.4a), (2.4b) and (2.5) of [1]_.

    Raises
    ------
    DomainError
        If ``E_m`` or the mixture ratio is out of range.

    References
    ----------
    .. [1] NASA-CR187109-ROCCID. Muss, Nguyen and Johnson, NASA CR-187109
           (1991), section 2.2, equations (2.4)-(2.11), printed pp. 16-18.
           https://ntrs.nasa.gov/citations/19910014917
    """

    em = float(rupe_mixing_efficiency)
    if not isfinite(em) or not 0.25 < em <= 1.0:
        raise DomainError("Rupe mixing efficiency must lie in (0.25, 1].", field="rupe_mixing_efficiency")
    mean = _positive(mean_mixture_ratio, "mean_mixture_ratio")
    eta = 0.449 - 0.399 * em + 0.944 * em * em
    low, high = mean * eta, mean / eta
    low_fraction = (1.0 + low) / ((1.0 + eta) * (1.0 + mean))
    return StreamtubeSplit(em, eta, mean, low, high, low_fraction, 1.0 - low_fraction)


@dataclass(frozen=True, slots=True)
class EnergyReleaseEfficiency:
    """Mixing- and vaporization-limited c* and Isp efficiencies [1]_.

    Attributes
    ----------
    split : StreamtubeSplit
        Streamtube distribution.
    vaporization_efficiency : float
        ``eta_vap = W_vapor / W_inj`` (eq. (2.7)).
    delivered_characteristic_velocity_m_s : float
        ``eta_vap [X_R c*(MR_R) + X_L c*(MR_L)]`` (eq. (2.8), one zone).
    characteristic_velocity_efficiency : float
        Eq. (2.9).
    delivered_specific_impulse_s : float
        ``eta_vap [X_R I(MR_R) + X_L I(MR_L)]`` (eq. (2.10), one zone).
    specific_impulse_energy_release_efficiency : float
        Eq. (2.11).

    References
    ----------
    .. [1] NASA-CR187109-ROCCID. Muss, Nguyen and Johnson, NASA CR-187109
           (1991), section 2.2, equations (2.4)-(2.11), printed pp. 16-18.
           https://ntrs.nasa.gov/citations/19910014917
    """

    split: StreamtubeSplit
    vaporization_efficiency: float
    delivered_characteristic_velocity_m_s: float
    characteristic_velocity_efficiency: float
    delivered_specific_impulse_s: float
    specific_impulse_energy_release_efficiency: float


def energy_release_efficiency(
    *,
    rupe_mixing_efficiency: float,
    oxidizer_fuel_ratio: float,
    characteristic_velocity: Callable[[float], float],
    specific_impulse: Callable[[float], float],
    vaporization_efficiency: float = 1.0,
) -> EnergyReleaseEfficiency:
    """JANNAF simplified energy-release efficiency of one zone (ROCCID) [1]_.

    Parameters
    ----------
    rupe_mixing_efficiency : float
        Zone ``E_m`` as a fraction (see :func:`rupe_streamtube_split`).
    oxidizer_fuel_ratio : float
        Injected (zone mean) O/F.
    characteristic_velocity, specific_impulse : callable
        Theoretical ``c*(O/F)`` and ``I(O/F)``, for example from
        :func:`cea_characteristic_velocity_function`.
    vaporization_efficiency : float, default 1
        ``eta_vap`` in (0, 1].

    Returns
    -------
    EnergyReleaseEfficiency
        Equations (2.7)-(2.11) of [1]_ with one zone (``xi_C = 1``).

    Notes
    -----
    ROCCID evaluates the streamtube mixture ratios with the vapour mixture
    ratio at the throat; here the injected ratio is used, which is the case
    of equal vaporization of both propellants.

    References
    ----------
    .. [1] NASA-CR187109-ROCCID. Muss, Nguyen and Johnson, NASA CR-187109
           (1991), section 2.2, equations (2.4)-(2.11), printed pp. 16-18.
           https://ntrs.nasa.gov/citations/19910014917
    """

    split = rupe_streamtube_split(rupe_mixing_efficiency, oxidizer_fuel_ratio)
    vaporized = _fraction(vaporization_efficiency, "vaporization_efficiency", allow_zero=False)
    cstar = vaporized * (split.low_mass_fraction * characteristic_velocity(split.low_mixture_ratio)
                         + split.high_mass_fraction * characteristic_velocity(split.high_mixture_ratio))
    isp = vaporized * (split.low_mass_fraction * specific_impulse(split.low_mixture_ratio)
                       + split.high_mass_fraction * specific_impulse(split.high_mixture_ratio))
    return EnergyReleaseEfficiency(split, vaporized, cstar, cstar / characteristic_velocity(oxidizer_fuel_ratio),
                                   isp, isp / specific_impulse(oxidizer_fuel_ratio))


def maximum_mixing_diameter_ratio(*, element_type: str, outside_density_kg_m3: float,
                                  center_density_kg_m3: float, center_mass_flow_kg_s: float,
                                  outside_mass_flow_kg_s: float) -> float:
    """Return ``d_c/d_ou`` for maximum mixing efficiency, SP-8089 eq. (1) [1]_.

    Parameters
    ----------
    element_type : str
        Key of :data:`SP8089_MIXING_FACTOR` (``"1-on-1"`` ... ``"5-on-1"``).
    outside_density_kg_m3, center_density_kg_m3 : float
        Liquid densities of the outside and center streams.
    center_mass_flow_kg_s, outside_mass_flow_kg_s : float
        Total mass flow through all center and all outside orifices.

    Returns
    -------
    float
        ``sqrt(M [(rho_ou/rho_c)(w_c/w_ou)^2]^0.7)``.

    Raises
    ------
    InputError
        If the element type has no tabulated ``M``.

    Notes
    -----
    SP-8089 states the correlation for unlike-impinging liquid/liquid elements
    designed for a 60-degree included impingement angle [1]_.

    References
    ----------
    .. [1] NASA-SP8089-MIXING. NASA SP-8089 (1976), equation (1) and Table IV,
           printed pp. 19-20. https://ntrs.nasa.gov/citations/19760023196
    """

    try:
        factor = SP8089_MIXING_FACTOR[element_type]
    except KeyError as error:
        raise InputError(f"Unknown SP-8089 element type {element_type!r}.", field="element_type",
                         details={"supported": sorted(SP8089_MIXING_FACTOR)}) from error
    group = (_positive(outside_density_kg_m3, "outside_density_kg_m3")
             / _positive(center_density_kg_m3, "center_density_kg_m3")
             * (_positive(center_mass_flow_kg_s, "center_mass_flow_kg_s")
                / _positive(outside_mass_flow_kg_s, "outside_mass_flow_kg_s")) ** 2)
    return sqrt(factor * group ** 0.7)


@dataclass(frozen=True, slots=True)
class ChamberResidence:
    """Characteristic length and propellant stay time [1]_.

    Attributes
    ----------
    characteristic_length_m : float
        ``L* = V_c / A_t`` (SP-125 eq. (4-4)).
    stay_time_s : float
        ``t_s = V_c / (mdot v)`` (SP-125 eq. (4-3)).
    mass_flow_kg_s : float
        Mass flow used.
    recommended_range_m : tuple of float or None
        SP-125 Table 4-1 range for the requested propellant key.
    warnings : tuple of WarningMessage
        Out-of-range statements.

    References
    ----------
    .. [1] NASA-SP125-LSTAR. Huzel and Huang, NASA SP-125 (1967), equations
           (4-3) and (4-4) and Table 4-1, printed pp. 86-87.
           https://ntrs.nasa.gov/citations/19710019929
    """

    characteristic_length_m: float
    stay_time_s: float
    mass_flow_kg_s: float
    recommended_range_m: tuple[float, float] | None
    warnings: tuple[WarningMessage, ...]


def chamber_residence(*, chamber_volume_m3: float, throat_area_m2: float, mass_flow_kg_s: float,
                      gas_density_kg_m3: float, propellant_key: str | None = None) -> ChamberResidence:
    """Return ``L*`` and the stay time of a combustion chamber [1]_.

    Parameters
    ----------
    chamber_volume_m3, throat_area_m2 : float
        Chamber volume up to the throat and throat area.
    mass_flow_kg_s : float
        Propellant mass flow.
    gas_density_kg_m3 : float
        Average gas density, the inverse of SP-125's average specific volume.
    propellant_key : str, optional
        Key of :data:`SP125_RECOMMENDED_CHARACTERISTIC_LENGTH_M` to compare
        with.

    Returns
    -------
    ChamberResidence
        ``L*``, stay time and range check.

    Raises
    ------
    InputError
        If the propellant key is unknown.

    References
    ----------
    .. [1] NASA-SP125-LSTAR. Huzel and Huang, NASA SP-125 (1967), equations
           (4-3) and (4-4) and Table 4-1, printed pp. 86-87.
           https://ntrs.nasa.gov/citations/19710019929
    """

    volume = _positive(chamber_volume_m3, "chamber_volume_m3")
    length = volume / _positive(throat_area_m2, "throat_area_m2")
    flow = _positive(mass_flow_kg_s, "mass_flow_kg_s")
    stay = volume * _positive(gas_density_kg_m3, "gas_density_kg_m3") / flow
    recommended = None
    warnings: list[WarningMessage] = []
    if propellant_key is not None:
        key = propellant_key.strip().upper()
        if key not in SP125_RECOMMENDED_CHARACTERISTIC_LENGTH_M:
            raise InputError(f"No SP-125 L* range for {propellant_key!r}.", field="propellant_key",
                             details={"supported": sorted(SP125_RECOMMENDED_CHARACTERISTIC_LENGTH_M)})
        recommended = SP125_RECOMMENDED_CHARACTERISTIC_LENGTH_M[key]
        if not recommended[0] <= length <= recommended[1]:
            warnings.append(WarningMessage(
                code="characteristic_length_outside_sp125_range",
                message=f"L* = {length:.3f} m lies outside SP-125 Table 4-1 range "
                        f"{recommended[0]:.3f}-{recommended[1]:.3f} m for {key}.",
                field="chamber_volume_m3"))
    return ChamberResidence(length, stay, flow, recommended, tuple(warnings))


def injector_to_throat_total_pressure_ratio(*, static_to_total_pressure_ratio: float,
                                            subsonic_specific_impulse_m_s: float,
                                            injection_velocity_m_s: float,
                                            characteristic_velocity_m_s: float,
                                            contraction_ratio: float) -> float:
    """Return ``P_c,T / P_c,a`` from the momentum balance of TP-1999-208522 eq. (4) [1]_.

    Parameters
    ----------
    static_to_total_pressure_ratio : float
        ``P_s / P_T`` at the chamber end.
    subsonic_specific_impulse_m_s : float
        Theoretical subsonic specific impulse at the chamber end in SI
        (``I g_c``, the flow velocity for ``p_e = p_a``, RP-1311 eq. (6.8)).
    injection_velocity_m_s : float
        Mass-averaged propellant injection velocity ``V_av``, non-negative.
    characteristic_velocity_m_s : float
        Theoretical ``c*``.
    contraction_ratio : float
        ``A_c / A_t``, greater than one.

    Returns
    -------
    float
        ``[P_s/P_T + (I g_c - V_av) / (c* eps_c)]^-1``.

    Notes
    -----
    With ``V_av = 0`` and ``c* = p_inf A_t / mdot`` the relation is the
    constant-area momentum balance ``p_inj = (p + rho u^2)_c`` of RP-1311
    eq. (6.28); the finite-area-combustor solution satisfies it identically.

    References
    ----------
    .. [1] NASA-TP1999-208522. Jankovsky, Smith and Pavli, NASA/TP-1999-208522,
           equations (3)-(11), printed pp. 12-13; loss chain printed p. 25 and
           Figures 20 and 24. https://ntrs.nasa.gov/citations/19990051032
    """

    ratio = _fraction(static_to_total_pressure_ratio, "static_to_total_pressure_ratio", allow_zero=False)
    impulse = _positive(subsonic_specific_impulse_m_s, "subsonic_specific_impulse_m_s")
    velocity = float(injection_velocity_m_s)
    if not isfinite(velocity) or velocity < 0.0:
        raise DomainError("Injection velocity must be non-negative.", field="injection_velocity_m_s")
    contraction = _positive(contraction_ratio, "contraction_ratio")
    if contraction <= 1.0:
        raise DomainError("Contraction ratio must exceed one.", field="contraction_ratio")
    denominator = ratio + (impulse - velocity) / (_positive(characteristic_velocity_m_s,
                                                            "characteristic_velocity_m_s") * contraction)
    if denominator <= 0.0:
        raise DomainError("Momentum balance gives a non-positive pressure ratio.",
                          field="injection_velocity_m_s")
    return 1.0 / denominator


# ---------------------------------------------------------------------------
# Isp loss chain
# ---------------------------------------------------------------------------


@dataclass(frozen=True, slots=True)
class IspLossChain:
    """JANNAF-ordered vacuum-Isp loss chain [1]_, [2]_.

    Attributes
    ----------
    ode_specific_impulse_s : float
        One-dimensional equilibrium (ideal) vacuum Isp.
    odk_specific_impulse_s : float
        After kinetic losses, ``I_ODE eta_kin``.
    two_dimensional_specific_impulse_s : float
        After divergence and shock losses, ``I_ODK eta_div``.
    boundary_layer_specific_impulse_s : float
        ``I_2D - Delta I_BL``.
    delivered_specific_impulse_s : float
        ``I_BL eta_c*`` (energy-release losses).
    kinetic_efficiency, divergence_efficiency : float
        Inputs.
    boundary_layer_efficiency : float
        ``I_BL / I_2D``.
    energy_release_efficiency : float
        The c* efficiency input.
    overall_efficiency : float
        ``I_delivered / I_ODE``, the product of the four factors.
    kinetic_loss_s, divergence_loss_s, boundary_layer_loss_s, energy_release_loss_s : float
        Isp decrements in chain order; they sum to ``I_ODE - I_delivered``.
    loss_budget : object
        :class:`rocket_propulsion.propulsion.losses.PerformanceLossBudget`
        of the post-kinetic factors with ``I_ODK`` as ideal and the discharge
        factor disabled.

    References
    ----------
    .. [1] NASA-TP1999-208522. Jankovsky, Smith and Pavli, NASA/TP-1999-208522,
           equations (3)-(11), printed pp. 12-13; loss chain printed p. 25 and
           Figures 20 and 24. https://ntrs.nasa.gov/citations/19990051032
    .. [2] NASA-TM89814-JANNAF. Davidian, NASA TM-89814 (1987), printed
           pp. 1-4. https://ntrs.nasa.gov/citations/19870010946
    """

    ode_specific_impulse_s: float
    odk_specific_impulse_s: float
    two_dimensional_specific_impulse_s: float
    boundary_layer_specific_impulse_s: float
    delivered_specific_impulse_s: float
    kinetic_efficiency: float
    divergence_efficiency: float
    boundary_layer_efficiency: float
    energy_release_efficiency: float
    overall_efficiency: float
    kinetic_loss_s: float
    divergence_loss_s: float
    boundary_layer_loss_s: float
    energy_release_loss_s: float
    loss_budget: object


def jannaf_isp_loss_chain(*, ode_specific_impulse_s: float, kinetic_efficiency: float,
                          divergence_efficiency: float, boundary_layer_decrement_s: float,
                          characteristic_velocity_efficiency: float,
                          reference_thrust_n: float = 1.0) -> IspLossChain:
    """Apply the JANNAF loss sequence to an ODE vacuum Isp [1]_, [2]_.

    Parameters
    ----------
    ode_specific_impulse_s : float
        One-dimensional equilibrium vacuum Isp.
    kinetic_efficiency : float
        ``I_ODK / I_ODE`` in (0, 1] (from a kinetics solver; frozen and
        shifting equilibrium are its two limits, see
        :func:`frozen_equilibrium_limits`).
    divergence_efficiency : float
        Two-dimensional inviscid (divergence and shock) factor in (0, 1].
    boundary_layer_decrement_s : float
        Boundary-layer impulse decrement, non-negative and smaller than the
        two-dimensional Isp.
    characteristic_velocity_efficiency : float
        Energy-release (c*) efficiency in (0, 1].
    reference_thrust_n : float, default 1
        Thrust attached to ``I_ODK`` for the hand-off to
        :func:`~rocket_propulsion.propulsion.losses.calculate_loss_budget`.

    Returns
    -------
    IspLossChain
        Stations, factors and decrements.

    Notes
    -----
    Order and treatment follow NASA/TP-1999-208522 printed p. 25 and
    Figures 20 and 24: ODE -> ODK (kinetic losses) -> MOC (divergence and
    shock) -> TDK/MABL (boundary layer) -> multiplication by ``eta_c*``
    (energy-release losses) [1]_; the boundary-layer thrust decrement is
    converted to an impulse decrement and subtracted [2]_.

    References
    ----------
    .. [1] NASA-TP1999-208522. Jankovsky, Smith and Pavli, NASA/TP-1999-208522,
           equations (3)-(11), printed pp. 12-13; loss chain printed p. 25 and
           Figures 20 and 24. https://ntrs.nasa.gov/citations/19990051032
    .. [2] NASA-TM89814-JANNAF. Davidian, NASA TM-89814 (1987), printed
           pp. 1-4. https://ntrs.nasa.gov/citations/19870010946
    """

    from rocket_propulsion.propulsion.losses import calculate_loss_budget

    ode = _positive(ode_specific_impulse_s, "ode_specific_impulse_s")
    kinetic = _fraction(kinetic_efficiency, "kinetic_efficiency", allow_zero=False)
    divergence = _fraction(divergence_efficiency, "divergence_efficiency", allow_zero=False)
    combustion = _fraction(characteristic_velocity_efficiency, "characteristic_velocity_efficiency",
                           allow_zero=False)
    decrement = float(boundary_layer_decrement_s)
    odk = ode * kinetic
    two_d = odk * divergence
    if not isfinite(decrement) or not 0.0 <= decrement < two_d:
        raise DomainError("Boundary-layer decrement must be non-negative and below the 2-D Isp.",
                          field="boundary_layer_decrement_s")
    boundary = two_d - decrement
    delivered = boundary * combustion
    boundary_efficiency = boundary / two_d
    budget = calculate_loss_budget(
        ideal_thrust_n=_positive(reference_thrust_n, "reference_thrust_n"), ideal_specific_impulse_s=odk,
        combustion_efficiency=combustion, boundary_layer_efficiency=boundary_efficiency,
        divergence_efficiency=divergence, enable_discharge_loss=False)
    return IspLossChain(
        ode_specific_impulse_s=ode,
        odk_specific_impulse_s=odk,
        two_dimensional_specific_impulse_s=two_d,
        boundary_layer_specific_impulse_s=boundary,
        delivered_specific_impulse_s=delivered,
        kinetic_efficiency=kinetic,
        divergence_efficiency=divergence,
        boundary_layer_efficiency=boundary_efficiency,
        energy_release_efficiency=combustion,
        overall_efficiency=delivered / ode,
        kinetic_loss_s=ode - odk,
        divergence_loss_s=odk - two_d,
        boundary_layer_loss_s=decrement,
        energy_release_loss_s=boundary - delivered,
        loss_budget=budget,
    )


@dataclass(frozen=True, slots=True)
class ChemistryLimits:
    """Shifting-equilibrium and frozen vacuum Isp at one point [1]_.

    Attributes
    ----------
    equilibrium_specific_impulse_s, frozen_specific_impulse_s : float
        Vacuum Isp with instantaneous and with infinitely slow reactions.
    frozen_to_equilibrium_ratio : float
        ``I_frozen / I_equilibrium``.

    References
    ----------
    .. [1] NASA-RP1311-P1-FAC. Gordon and McBride, NASA RP-1311 Part I (1994),
           chapter 6, sections 6.1-6.4, equations (6.1)-(6.33).
           https://ntrs.nasa.gov/citations/19950013764
    """

    equilibrium_specific_impulse_s: float
    frozen_specific_impulse_s: float
    frozen_to_equilibrium_ratio: float


def frozen_equilibrium_limits(*, propellant_pair: str, oxidizer_fuel_ratio: float, chamber_pressure_pa: float,
                              area_ratio: float,
                              provider: EquilibriumProvider | None = None) -> ChemistryLimits:
    """Return the two chemistry limits of the nozzle expansion [1]_.

    Parameters
    ----------
    propellant_pair : str
        Preset key.
    oxidizer_fuel_ratio, chamber_pressure_pa, area_ratio : float
        Operating point.
    provider : EquilibriumProvider, optional
        Equilibrium provider.

    Returns
    -------
    ChemistryLimits
        Vacuum Isp for instantaneous (shifting) and infinitely slow (frozen)
        reactions, RP-1311 sections 6.1 and 6.5 [1]_. A finite-rate (ODK)
        result is not computed here.

    References
    ----------
    .. [1] NASA-RP1311-P1-FAC. Gordon and McBride, NASA RP-1311 Part I (1994),
           chapter 6, sections 6.1-6.4, equations (6.1)-(6.33).
           https://ntrs.nasa.gov/citations/19950013764
    """

    pair = get_propellant_pair(propellant_pair)
    mixture = reactants_for_mass_ratio(pair, oxidizer_fuel_ratio)
    active = provider or LocalEquilibriumProvider()
    values = [solve_rocket_performance(reactants=mixture, chamber_pressure_pa=chamber_pressure_pa,
                                       area_ratios=(area_ratio,), chemistry=chemistry,
                                       provider=active).exits[0].vacuum_specific_impulse_s
              for chemistry in _CHEMISTRY_LIMITS]
    return ChemistryLimits(values[0], values[1], values[1] / values[0])


# ---------------------------------------------------------------------------
# CEA-consistent property callables and O/F optimisation
# ---------------------------------------------------------------------------


def cea_characteristic_velocity_function(
    *, propellant_pair: str, chamber_pressure_pa: float, area_ratio: float = 40.0,
    provider: EquilibriumProvider | None = None,
) -> tuple[Callable[[float], float], Callable[[float], float]]:
    """Return memoized ``c*(O/F)`` and vacuum ``I(O/F)`` callables [1]_.

    Parameters
    ----------
    propellant_pair : str
        Preset key.
    chamber_pressure_pa, area_ratio : float
        Operating point.
    provider : EquilibriumProvider, optional
        Equilibrium provider.

    Returns
    -------
    tuple of callable
        ``(c*, I_vac)`` in m/s and s from the shifting-equilibrium
        infinite-area solve, cached by O/F.

    References
    ----------
    .. [1] NASA-RP1311-P1-FAC. Gordon and McBride, NASA RP-1311 Part I (1994),
           chapter 6, sections 6.1-6.4, equations (6.1)-(6.33).
           https://ntrs.nasa.gov/citations/19950013764
    """

    pair = get_propellant_pair(propellant_pair)
    active = provider or LocalEquilibriumProvider()
    cache: dict[float, tuple[float, float]] = {}

    def evaluate(ratio: float) -> tuple[float, float]:
        key = round(float(ratio), 12)
        if key not in cache:
            result = solve_rocket_performance(reactants=reactants_for_mass_ratio(pair, key),
                                              chamber_pressure_pa=chamber_pressure_pa, area_ratios=(area_ratio,),
                                              provider=active)
            cache[key] = (result.characteristic_velocity_m_s, result.exits[0].vacuum_specific_impulse_s)
        return cache[key]

    return (lambda ratio: evaluate(ratio)[0]), (lambda ratio: evaluate(ratio)[1])


@dataclass(frozen=True, slots=True)
class MixtureRatioOptimum:
    """Best O/F for one objective [1]_, [2]_.

    Attributes
    ----------
    pair : str
        Preset key.
    objective : str
        ``"vacuum_specific_impulse"`` (s) or ``"density_vacuum_impulse"``
        (``rho_0 I_vac``, kg s/m^3).
    oxidizer_fuel_ratio : float
        Maximizing O/F.
    objective_value : float
        Objective at the optimum.
    vacuum_specific_impulse_s : float
        Vacuum Isp at the optimum.
    bulk_density_kg_m3 : float or None
        ``rho_0`` at the optimum for the density objective.
    bracket : tuple of float
        Search interval.
    evaluations : int
        Rocket solves performed.
    at_bracket_edge : bool
        Whether the optimum lies within the tolerance of a bracket end.
    warnings : tuple of WarningMessage
        Bracket-edge statement.

    References
    ----------
    .. [1] NASA-RP1311-P1-EQRATIO. Gordon and McBride, NASA RP-1311 Part I
           (1994), chapter 9, equations (9.10)-(9.19), printed pp. 42-43.
           https://ntrs.nasa.gov/citations/19950013764
    .. [2] NASA-SP125-DENSITY-IMPULSE. Huzel and Huang, NASA SP-125 (1967),
           equation (1-45), printed p. 20.
           https://ntrs.nasa.gov/citations/19710019929
    """

    pair: str
    objective: str
    oxidizer_fuel_ratio: float
    objective_value: float
    vacuum_specific_impulse_s: float
    bulk_density_kg_m3: float | None
    bracket: tuple[float, float]
    evaluations: int
    at_bracket_edge: bool
    warnings: tuple[WarningMessage, ...]


def optimize_mixture_ratio(
    *,
    propellant_pair: str,
    chamber_pressure_pa: float,
    area_ratio: float,
    objective: Objective = "vacuum_specific_impulse",
    bracket: tuple[float, float] | None = None,
    oxidizer_density_kg_m3: float | None = None,
    fuel_density_kg_m3: float | None = None,
    tolerance: float = 1.0e-3,
    provider: EquilibriumProvider | None = None,
) -> MixtureRatioOptimum:
    """Maximize vacuum Isp or density impulse over O/F [1]_, [2]_.

    Parameters
    ----------
    propellant_pair : str
        Preset key.
    chamber_pressure_pa, area_ratio : float
        Operating point (vacuum Isp at this area ratio).
    objective : {"vacuum_specific_impulse", "density_vacuum_impulse"}
        Quantity to maximize; the density objective is ``rho_0 I_vac`` with
        ``rho_0`` from RP-1311 eq. (9.13) and density impulse from SP-125
        eq. (1-45).
    bracket : tuple of float, optional
        O/F interval; the preset sweep range by default.
    oxidizer_density_kg_m3, fuel_density_kg_m3 : float, optional
        Required for the density objective; no propellant densities are
        bundled.
    tolerance : float, default 1e-3
        Final O/F interval width, in ``[1e-5, 0.5]``.
    provider : EquilibriumProvider, optional
        Equilibrium provider.

    Returns
    -------
    MixtureRatioOptimum
        Optimum O/F and objective.

    Raises
    ------
    InputError
        If the objective is unknown or densities are missing.
    DomainError
        If the bracket or tolerance is invalid.

    Notes
    -----
    Golden-section search assumes a unimodal objective in the bracket. Each
    evaluation is the shifting-equilibrium infinite-area rocket solve that
    reproduces NASA CEA 3.3.4 (``tests/reference_data/cea_gate_points.json``).

    References
    ----------
    .. [1] NASA-RP1311-P1-EQRATIO. Gordon and McBride, NASA RP-1311 Part I
           (1994), chapter 9, equations (9.10)-(9.19), printed pp. 42-43.
           https://ntrs.nasa.gov/citations/19950013764
    .. [2] NASA-SP125-DENSITY-IMPULSE. Huzel and Huang, NASA SP-125 (1967),
           equation (1-45), printed p. 20.
           https://ntrs.nasa.gov/citations/19710019929
    """

    if objective not in ("vacuum_specific_impulse", "density_vacuum_impulse"):
        raise InputError("Unknown mixture-ratio objective.", field="objective")
    pair = get_propellant_pair(propellant_pair)
    low, high = bracket if bracket is not None else (pair.minimum_oxidizer_fuel_ratio,
                                                     pair.maximum_oxidizer_fuel_ratio)
    low, high = float(low), float(high)
    if not (isfinite(low) and isfinite(high)) or low <= 0.0 or high <= low:
        raise DomainError("O/F bracket must be positive and increasing.", field="bracket")
    if not isfinite(tolerance) or not 1.0e-5 <= tolerance <= 0.5:
        raise DomainError("O/F tolerance must lie in [1e-5, 0.5].", field="tolerance")
    densities = None
    if objective == "density_vacuum_impulse":
        if oxidizer_density_kg_m3 is None or fuel_density_kg_m3 is None:
            raise InputError("The density objective needs oxidizer and fuel densities.",
                             field="oxidizer_density_kg_m3")
        densities = (_positive(oxidizer_density_kg_m3, "oxidizer_density_kg_m3"),
                     _positive(fuel_density_kg_m3, "fuel_density_kg_m3"))
    _, isp = cea_characteristic_velocity_function(propellant_pair=pair.key,
                                                  chamber_pressure_pa=chamber_pressure_pa,
                                                  area_ratio=area_ratio, provider=provider)
    count = 0

    def value(ratio: float) -> float:
        nonlocal count
        count += 1
        result = isp(ratio)
        if densities is not None:
            result = density_impulse(result, bulk_density(ratio, *densities))
        return result

    a, b = low, high
    first, second = b - _GOLDEN * (b - a), a + _GOLDEN * (b - a)
    f_first, f_second = value(first), value(second)
    while b - a > tolerance:
        if f_first > f_second:
            b, second, f_second = second, first, f_first
            first = b - _GOLDEN * (b - a)
            f_first = value(first)
        else:
            a, first, f_first = first, second, f_second
            second = a + _GOLDEN * (b - a)
            f_second = value(second)
    best = 0.5 * (a + b)
    objective_value = value(best)
    edge = best - low <= tolerance or high - best <= tolerance
    warnings = (WarningMessage(code="optimum_at_bracket_edge",
                               message="The optimum lies at the end of the O/F bracket; widen the bracket.",
                               field="bracket"),) if edge else ()
    return MixtureRatioOptimum(
        pair=pair.key,
        objective=objective,
        oxidizer_fuel_ratio=best,
        objective_value=objective_value,
        vacuum_specific_impulse_s=isp(best),
        bulk_density_kg_m3=bulk_density(best, *densities) if densities is not None else None,
        bracket=(low, high),
        evaluations=count,
        at_bracket_edge=edge,
        warnings=warnings,
    )
