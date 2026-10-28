"""Chemical equilibrium with pure condensed-phase products (Gordon-McBride method).

:class:`CondensedEquilibriumProvider` minimizes the Gibbs energy of an ideal-gas
mixture plus pure condensed species (graphite, ice, liquid water, aluminium,
alumina, ...) subject to element conservation, following NASA RP-1311 Part I
[1]_: the reduced Newton-Raphson iteration of section 2.3.2 (eqs. (2.24)-(2.28)
with eq. (2.18) for the gas corrections), the control factor and convergence
tests of section 3.3 (eqs. (3.1)-(3.6c)), the condensed-phase inclusion test of
section 3.4 (eq. (3.7)) with removal of species whose amount becomes negative,
and the phase-transition rules of section 3.5. Thermodynamic data are NASA's
own ``thermo.inp`` records (the bundled gas catalog plus
``data/nasa_cea_condensed.json``) [2]_.

:func:`solve_condensed_rocket_performance` solves the infinite-area-combustor
rocket problem of RP-1311 chapter 6 [1]_ with these two-phase equilibrium
states: zero velocity and thermal lag between condensed and gaseous species,
condensed species occupy no volume (section 2.2), isentropic expansion from
the assigned-enthalpy chamber. It reuses the station, exit and performance
records and the root helpers of
:mod:`~rocket_propulsion.thermochemistry.rocket_performance`. The provider also
implements the
:class:`~rocket_propulsion.thermochemistry.equilibrium.EquilibriumProvider`
protocol; :func:`~rocket_propulsion.thermochemistry.rocket_performance.solve_rocket_performance`
accepts it only when every chamber species above a mole fraction of 1e-10 is a
gas of the base catalog, because that function freezes the chamber
composition as a gas :class:`~rocket_propulsion.thermochemistry.mixtures.Mixture`.

References
----------
.. [1] NASA-RP1311-P1-CONDENSED. Gordon and McBride, NASA RP-1311 Part I
       (1994). https://ntrs.nasa.gov/citations/19950013764
.. [2] NASA-CEA-THERMO. NASA Glenn Research Center, CEA thermodynamic database
       ``data/thermo.inp``. https://github.com/nasa/cea
"""

from __future__ import annotations

import json
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from functools import lru_cache
from math import exp, isfinite, log, sqrt
from pathlib import Path
from typing import TypedDict

from rocket_propulsion.core.constants import STANDARD_GRAVITY_M_S2
from rocket_propulsion.core.errors import ConvergenceError, DomainError, InputError
from rocket_propulsion.core.metadata import WarningMessage
from rocket_propulsion.core.validation import coerce_fields, real, text

from ._records import _owned_catalog, _owned_formula
from .equilibrium import (
    EquilibriumConvergence,
    EquilibriumProblem,
    EquilibriumResiduals,
    EquilibriumResult,
    EquilibriumSpeciesAmount,
    ReactantAmount,
)
from .mixtures import REFERENCE_PRESSURE_PA
from .nasa_polynomials import (
    UNIVERSAL_GAS_CONSTANT_J_MOL_K,
    NasaIntervalRecord,
    NasaPolynomialInterval,
    _owned_intervals,
    interval_from_record,
)
from .rocket_performance import (
    RocketExit,
    RocketPerformance,
    RocketStation,
    _bracket,
    _illinois,
    _maximize_bounded,
)
from .species import CondensedReactant, Species, get_reactant, species_catalog

CONDENSED_DATA_FILE = Path(__file__).parent / "data" / "nasa_cea_condensed.json"

# RP-1311 Part I section 3.2 and 3.3 constants.
_SIZE = 18.420681  # -ln 1e-8
_LAMBDA2_LOG = 9.2103404  # -ln 1e-4 in eq. (3.2)
# Weight of the Tikhonov term used only when the Newton matrix is singular
# (see _regularized); relative to the unit diagonal of the scaled normal matrix.
_REGULARIZATION = 1.0e-12
# Gas moles below this fraction of all moles mean that no gas phase exists
# (e.g. stoichiometric H2/O2 below the boiling point), which the ideal-gas
# formulation of RP-1311 cannot represent; such states are rejected.
_MINIMUM_GAS_FRACTION = 1.0e-8
_INITIAL_TOTAL_MOL_KG = 100.0  # n = 0.1 kg-mol/kg of section 3.1
_INITIAL_TEMPERATURE_K = 3800.0  # section 3.1
_TRANSITION_BAND_K = 50.0  # section 3.5
_MINIMUM_TEMPERATURE_K = 200.0
_MAXIMUM_TEMPERATURE_K = 6000.0
# Species above this mole fraction must lie inside their data range at the
# converged temperature; smaller ones are evaluated with the nearest interval
# and reported in a warning.
_RANGE_CHECK_MOLE_FRACTION = 1.0e-6


@dataclass(frozen=True, slots=True)
class CondensedSpecies:
    """One pure condensed product phase from NASA ``thermo.inp`` [1]_.

    Attributes
    ----------
    key : str
        ``thermo.inp`` name, for example ``"AL2O3(L)"``.
    formula : tuple of (str, float)
        Element symbols and atom counts.
    molecular_mass_kg_mol : float
        Molar mass in kg/mol.
    heat_of_formation_j_mol_298 : float
        Header heat of formation of the substance at 298.15 K, J/mol.
    phase_group : str
        Substance shared by all phases, for example ``"AL2O3"``.
    phase_index : int
        Phase number of the ``thermo.inp`` header.
    intervals : tuple of NasaPolynomialInterval
        Contiguous polynomial intervals; their union is the phase's range.
    source, source_url : str
        Title and URL of the data source.

    References
    ----------
    .. [1] NASA-CEA-THERMO. NASA Glenn Research Center, CEA thermodynamic
           database ``data/thermo.inp``. https://github.com/nasa/cea
    """

    key: str
    formula: tuple[tuple[str, float], ...]
    molecular_mass_kg_mol: float
    heat_of_formation_j_mol_298: float
    phase_group: str
    phase_index: int
    intervals: tuple[NasaPolynomialInterval, ...]
    source: str
    source_url: str

    def __post_init__(self) -> None:
        coerce_fields(self, 'Condensed species')
        if self.molecular_mass_kg_mol <= 0 or self.phase_index < 1:
            raise DomainError('A condensed species needs positive molar mass and phase index.')
        object.__setattr__(self, 'key', self.key.strip())
        object.__setattr__(self, 'phase_group', self.phase_group.strip())
        object.__setattr__(self, 'formula', _owned_formula(self.formula))
        object.__setattr__(self, 'intervals', _owned_intervals(self.intervals))

    @property
    def minimum_temperature_k(self) -> float:
        """Lower bound of the phase's data range, K."""
        return self.intervals[0].minimum_temperature_k

    @property
    def maximum_temperature_k(self) -> float:
        """Upper bound of the phase's data range, K."""
        return self.intervals[-1].maximum_temperature_k

    def covers(self, temperature_k: float) -> bool:
        """Return whether ``temperature_k`` lies in the phase's data range."""
        temperature_k = real(temperature_k, 'temperature_k')
        return self.minimum_temperature_k <= temperature_k <= self.maximum_temperature_k


@dataclass(frozen=True, slots=True)
class CondensedCatalog:
    """Owned snapshot of products and reactants for the condensed-phase provider.

    Its dictionaries may be edited locally; no edit changes a subsequent
    catalogue or the provider's internal data. Species/phase records own
    immutable nested formulae and intervals.

    Attributes
    ----------
    gases : dict of str to Species
        Base gas catalog plus the gases of ``nasa_cea_condensed.json``, keyed
        by upper-case name.
    condensed : dict of str to CondensedSpecies
        Condensed product phases keyed by upper-case name.
    reactants : dict of str to CondensedReactant
        Reactants with an enthalpy assigned at 298.15 K (``AL(CR)``,
        ``NH4CLO4(I)``).
    """

    gases: dict[str, Species]
    condensed: dict[str, CondensedSpecies]
    reactants: dict[str, CondensedReactant]

    def __post_init__(self) -> None:
        gases = _owned_catalog(self.gases, Species, 'gases')
        condensed = _owned_catalog(self.condensed, CondensedSpecies, 'condensed')
        reactants = _owned_catalog(self.reactants, CondensedReactant, 'reactants')
        if gases.keys() & condensed.keys():
            raise DomainError('Gas and condensed product keys must be disjoint.')
        object.__setattr__(self, 'gases', dict(gases))
        object.__setattr__(self, 'condensed', dict(condensed))
        object.__setattr__(self, 'reactants', dict(reactants))


class _PhaseRecord(TypedDict):
    """Fields shared by every species record of ``nasa_cea_condensed.json``."""

    formula: dict[str, float]
    intervals: list[NasaIntervalRecord]


def _intervals(record: _PhaseRecord) -> tuple[NasaPolynomialInterval, ...]:
    return tuple(interval_from_record(item) for item in record["intervals"])


def _formula(record: _PhaseRecord) -> tuple[tuple[str, float], ...]:
    return _owned_formula(tuple(record["formula"].items()))


@lru_cache(maxsize=1)
def _condensed_catalog_tables() -> tuple[Mapping[str, Species], Mapping[str, CondensedSpecies], Mapping[str, CondensedReactant]]:
    """Load the gas, condensed-product and reactant records once [1]_.

    Returns
    -------
    tuple of mappings
        Private immutable gas, condensed-product and reactant tables;
        gas keys of the base catalog
        (:func:`~rocket_propulsion.thermochemistry.species.species_catalog`)
        and of ``nasa_cea_condensed.json`` are disjoint.

    Raises
    ------
    RuntimeError
        If a gas key occurs in both files or a key is duplicated.

    References
    ----------
    .. [1] NASA-CEA-THERMO. NASA Glenn Research Center, CEA thermodynamic
           database ``data/thermo.inp``. https://github.com/nasa/cea
    """

    document = json.loads(CONDENSED_DATA_FILE.read_text(encoding="utf-8"))
    source = document["source"]
    title, url = source["title"], source["url"]
    gases = dict(species_catalog())
    for record in document["gas_species"]:
        key = record["key"].upper()
        if key in gases:
            raise RuntimeError(f"Gas species {key} is defined twice.")
        gases[key] = Species(
            key=record["key"],
            name=record["name"],
            formula=_formula(record),
            molecular_mass_kg_mol=real(record["molecular_mass_g_mol"], 'molecular_mass_g_mol') / 1000.0,
            heat_of_formation_j_mol_298=record["heat_of_formation_j_mol_298"],
            intervals=_intervals(record),
            source=title,
            source_url=url,
        )
    condensed: dict[str, CondensedSpecies] = {}
    for record in document["condensed_species"]:
        key = record["key"].upper()
        if key in condensed or key in gases:
            raise RuntimeError(f"Condensed species {key} is defined twice.")
        condensed[key] = CondensedSpecies(
            key=record["key"],
            formula=_formula(record),
            molecular_mass_kg_mol=real(record["molecular_mass_g_mol"], 'molecular_mass_g_mol') / 1000.0,
            heat_of_formation_j_mol_298=record["heat_of_formation_j_mol_298"],
            phase_group=record["phase_group"],
            phase_index=record["phase_index"],
            intervals=_intervals(record),
            source=title,
            source_url=url,
        )
    reactants: dict[str, CondensedReactant] = {}
    for record in document['reactants']:
        parsed = CondensedReactant(
            key=record["key"],
            name=record["name"],
            formula=_formula(record),
            molecular_mass_kg_mol=real(record["molecular_mass_g_mol"], 'molecular_mass_g_mol') / 1000.0,
            enthalpy_j_mol=record["assigned_enthalpy_j_mol"],
            temperature_k=record["assigned_temperature_k"],
            source=title,
            source_url=url,
        )
        key = parsed.key.upper()
        if key in reactants:
            raise RuntimeError('Duplicate condensed-reactant keys.')
        reactants[key] = parsed
    return (_owned_catalog(gases, Species, 'gases'),
            _owned_catalog(condensed, CondensedSpecies, 'condensed'),
            _owned_catalog(reactants, CondensedReactant, 'reactants'))


def condensed_catalog() -> CondensedCatalog:
    """Return an independent, serializable snapshot of the cached data [1]_.

    All three dictionaries are owned copies. Edits cannot change later
    callers or thermodynamic calculations; source-labelled records are
    immutable. ``dataclasses.asdict`` remains supported.

    References
    ----------
    .. [1] NASA-CEA-THERMO. NASA Glenn Research Center, CEA thermodynamic
           database ``data/thermo.inp``. https://github.com/nasa/cea
    """
    gases, condensed, reactants = _condensed_catalog_tables()
    return CondensedCatalog(gases=dict(gases), condensed=dict(condensed), reactants=dict(reactants))


def get_condensed_reactant(key: str) -> Species | CondensedReactant:
    """Return a reactant of the condensed catalog or of the base catalog.

    Parameters
    ----------
    key : str
        Reactant key such as ``"AL(cr)"``, ``"NH4CLO4(I)"``, ``"O2(L)"`` or
        ``"O2"``; case and surrounding blanks are ignored.

    Returns
    -------
    Species or CondensedReactant
        Gas reactants are evaluated at the problem's reactant temperature;
        condensed reactants carry an assigned enthalpy.

    Raises
    ------
    InputError
        If the key is unknown (code ``unknown_species``).

    References
    ----------
    .. [1] NASA-CEA-THERMO. NASA Glenn Research Center, CEA thermodynamic
           database ``data/thermo.inp``. https://github.com/nasa/cea
    """

    reactant = condensed_catalog().reactants.get(text(key, 'species').strip().upper())
    return reactant if reactant is not None else get_reactant(key)


def _validated_reactants(values: Mapping[str, float]) -> tuple[ReactantAmount, ...]:
    if not values:
        raise InputError("At least one reactant is required.", field="reactants")
    result = []
    for key, raw in values.items():
        try:
            amount = float(raw)
        except (TypeError, ValueError) as error:
            raise InputError(f"Reactant amount for {key} must be numeric.", field="reactants") from error
        if not isfinite(amount) or amount <= 0.0:
            raise DomainError(f"Reactant amount for {key} must be finite and positive.",
                              field="reactants")
        get_condensed_reactant(str(key))
        result.append(ReactantAmount(str(key).strip().upper(), amount))
    return tuple(sorted(result, key=lambda item: item.species))


@dataclass(frozen=True, slots=True)
class CondensedEquilibriumProblem:
    """A TP, HP or SP equilibrium problem whose products may include condensed phases.

    Use the :meth:`tp`, :meth:`hp` and :meth:`sp` constructors, which accept
    the reactants of :func:`get_condensed_reactant`.

    Attributes
    ----------
    problem_type : {"TP", "HP", "SP"}
        Assigned temperature, enthalpy or entropy, each with pressure, as in
        RP-1311 Part I eqs. (2.12), (2.13) and (2.15) [1]_.
    reactants : tuple of ReactantAmount
        Reactant amounts in mol, sorted by key.
    pressure_pa : float
        Pressure in Pa.
    temperature_k : float or None
        Temperature of a TP problem.
    reactant_temperature_k : float
        Temperature of gas reactants for the default HP enthalpy target.
    target_enthalpy_j, target_entropy_j_k : float or None
        Total enthalpy (HP, default: reactant enthalpy) or total entropy (SP,
        required) of the products.
    candidate_species : tuple of str or None
        Product keys (gas or condensed) to consider; ``None`` means every
        bundled product made of the reactant elements.
    schema_version : str, default "1.0"
        Problem schema version.

    References
    ----------
    .. [1] NASA-RP1311-P1-CONDENSED. Gordon and McBride, NASA RP-1311 Part I
           (1994). https://ntrs.nasa.gov/citations/19950013764
    """

    problem_type: str
    reactants: tuple[ReactantAmount, ...]
    pressure_pa: float
    temperature_k: float | None = None
    reactant_temperature_k: float = 298.15
    target_enthalpy_j: float | None = None
    target_entropy_j_k: float | None = None
    candidate_species: tuple[str, ...] | None = None
    schema_version: str = "1.0"

    @classmethod
    def tp(cls, reactants: Mapping[str, float], *, temperature_k: float, pressure_pa: float,
           candidate_species: tuple[str, ...] | None = None) -> CondensedEquilibriumProblem:
        """Build an assigned temperature and pressure problem.

        Parameters
        ----------
        reactants : mapping of str to float
            Positive amount in mol per reactant key.
        temperature_k, pressure_pa : float
            Assigned state.
        candidate_species : tuple of str, optional
            Product keys to consider.

        Returns
        -------
        CondensedEquilibriumProblem
            Problem with ``problem_type="TP"``.

        References
        ----------
        .. [1] NASA-RP1311-P1-CONDENSED. Gordon and McBride, NASA RP-1311 Part I
               (1994). https://ntrs.nasa.gov/citations/19950013764
        """
        return cls("TP", _validated_reactants(reactants), pressure_pa, temperature_k=temperature_k,
                   candidate_species=candidate_species)

    @classmethod
    def hp(cls, reactants: Mapping[str, float], *, pressure_pa: float,
           reactant_temperature_k: float = 298.15, target_enthalpy_j: float | None = None,
           candidate_species: tuple[str, ...] | None = None) -> CondensedEquilibriumProblem:
        """Build an assigned enthalpy and pressure (adiabatic combustion) problem.

        Parameters
        ----------
        reactants : mapping of str to float
            Positive amount in mol per reactant key.
        pressure_pa : float
            Combustion pressure in Pa.
        reactant_temperature_k : float, default 298.15
            Temperature of gas reactants; condensed reactants use their
            assigned enthalpy.
        target_enthalpy_j : float, optional
            Explicit total enthalpy in J.
        candidate_species : tuple of str, optional
            Product keys to consider.

        Returns
        -------
        CondensedEquilibriumProblem
            Problem with ``problem_type="HP"``.

        References
        ----------
        .. [1] NASA-RP1311-P1-CONDENSED. Gordon and McBride, NASA RP-1311 Part I
               (1994). https://ntrs.nasa.gov/citations/19950013764
        """
        return cls("HP", _validated_reactants(reactants), pressure_pa,
                   reactant_temperature_k=reactant_temperature_k,
                   target_enthalpy_j=target_enthalpy_j, candidate_species=candidate_species)

    @classmethod
    def sp(cls, reactants: Mapping[str, float], *, pressure_pa: float, target_entropy_j_k: float,
           candidate_species: tuple[str, ...] | None = None) -> CondensedEquilibriumProblem:
        """Build an assigned entropy and pressure (isentropic) problem.

        Parameters
        ----------
        reactants : mapping of str to float
            Positive amount in mol per reactant key; fixes the element
            inventory and the mass.
        pressure_pa : float
            Pressure in Pa.
        target_entropy_j_k : float
            Total product entropy in J/K.
        candidate_species : tuple of str, optional
            Product keys to consider.

        Returns
        -------
        CondensedEquilibriumProblem
            Problem with ``problem_type="SP"``.

        References
        ----------
        .. [1] NASA-RP1311-P1-CONDENSED. Gordon and McBride, NASA RP-1311 Part I
               (1994). https://ntrs.nasa.gov/citations/19950013764
        """
        return cls("SP", _validated_reactants(reactants), pressure_pa,
                   target_entropy_j_k=target_entropy_j_k, candidate_species=candidate_species)


@dataclass(frozen=True, slots=True)
class CondensedPhaseAmount:
    """Equilibrium amount of one condensed phase.

    Attributes
    ----------
    species, phase_group : str
        Phase key (``"AL2O3(L)"``) and substance (``"AL2O3"``).
    moles : float
        Amount in mol.
    mass_fraction : float
        Mass of the phase over the total product mass.
    """

    species: str
    phase_group: str
    moles: float
    mass_fraction: float


@dataclass(frozen=True, slots=True)
class CondensedEquilibriumState:
    """Converged composition used to warm-start a neighbouring problem.

    Attributes
    ----------
    gas_keys : tuple of str
        Gas species in solver order.
    log_gas_amounts : tuple of float
        ``ln n_j`` in mol per kg of mixture.
    log_gas_total : float
        ``ln n`` in mol per kg.
    condensed : tuple of (str, float)
        Included condensed phases and their amounts in mol per kg.
    temperature_k : float
        Temperature of the state.
    """

    gas_keys: tuple[str, ...]
    log_gas_amounts: tuple[float, ...]
    log_gas_total: float
    condensed: tuple[tuple[str, float], ...]
    temperature_k: float


@dataclass(frozen=True, slots=True)
class CondensedEquilibriumSolution:
    """Equilibrium result plus the condensed-phase detail of one solve.

    Attributes
    ----------
    result : EquilibriumResult
        Provider-neutral result. ``total_moles`` counts gas and condensed
        species and mole fractions refer to that total; ``molecular_mass_kg_mol``
        is RP-1311's ``M = 1/n`` (total mass over gas moles, eq. (2.3)) so that
        ``p/(R/M T)`` is the mixture density of eq. (2.1).
    gas_moles : float
        Gas amount in mol.
    condensed : tuple of CondensedPhaseAmount
        Included condensed phases.
    condensed_mass_fraction : float
        Total condensed mass over product mass.
    mixture_molecular_mass_kg_mol : float
        ``MW`` of RP-1311 eq. (2.4a): mass over all moles.
    element_potentials : tuple of (str, float)
        Dimensionless ``pi_i = -lambda_i/(R T)`` of the final iteration.
    inclusion_tests : tuple of (str, float)
        ``(mu_c/RT) - sum_i pi_i a_ic`` (eq. (3.7)) for every condensed
        candidate that covers the final temperature but is not included; all
        are non-negative at a Gibbs minimum.
    events : tuple of str
        Phase changes made by the solver (``"include C(gr)"``, ``"remove ..."``,
        ``"switch ..."``, ``"add ... at transition"``).
    newton_iterations : int
        Newton iterations over all passes.
    state : CondensedEquilibriumState
        Converged state for warm starts.

    References
    ----------
    .. [1] NASA-RP1311-P1-CONDENSED. Gordon and McBride, NASA RP-1311 Part I
           (1994). https://ntrs.nasa.gov/citations/19950013764
    """

    result: EquilibriumResult
    gas_moles: float
    condensed: tuple[CondensedPhaseAmount, ...]
    condensed_mass_fraction: float
    mixture_molecular_mass_kg_mol: float
    element_potentials: tuple[tuple[str, float], ...]
    inclusion_tests: tuple[tuple[str, float], ...]
    events: tuple[str, ...]
    newton_iterations: int
    state: CondensedEquilibriumState


class _Singular(Exception):
    """Raised by :func:`_gauss` for a (numerically) singular matrix."""


def _gauss(matrix: list[list[float]], vector: list[float]) -> list[float]:
    """Solve a small dense system by Gaussian elimination with partial pivoting."""
    size = len(vector)
    rows = [row[:] + [value] for row, value in zip(matrix, vector, strict=True)]
    scale = max((abs(value) for row in matrix for value in row), default=0.0)
    if scale == 0.0:
        raise _Singular
    for column in range(size):
        pivot = max(range(column, size), key=lambda index: abs(rows[index][column]))
        if abs(rows[pivot][column]) <= 1.0e-13 * scale:
            raise _Singular
        rows[column], rows[pivot] = rows[pivot], rows[column]
        head = rows[column]
        inverse = 1.0 / head[column]
        for index in range(column + 1, size):
            factor = rows[index][column] * inverse
            if factor != 0.0:
                row = rows[index]
                for entry in range(column, size + 1):
                    row[entry] -= factor * head[entry]
    solution = [0.0] * size
    for index in range(size - 1, -1, -1):
        row = rows[index]
        total = row[size] - sum(row[entry] * solution[entry] for entry in range(index + 1, size))
        solution[index] = total / row[index]
    return solution


def _regularized(matrix: list[list[float]], vector: list[float], target: list[float]) -> list[float]:
    """Solve a singular Newton system by column-scaled Tikhonov least squares.

    Minimizes ``|A x - b|^2 + eps |D (x - target)|^2`` with ``D`` the column
    norms of ``A`` and ``eps = 1e-12``: directions the equations determine are
    unchanged to that order, and an undetermined direction keeps the value of
    the previous iterate (``target``).
    """
    size = len(vector)
    norms = [sqrt(sum(matrix[row][column] ** 2 for row in range(size))) or 1.0
             for column in range(size)]
    scaled = [[matrix[row][column] / norms[column] for column in range(size)] for row in range(size)]
    normal = [[sum(scaled[row][left] * scaled[row][right] for row in range(size))
               for right in range(size)] for left in range(size)]
    right_side = [sum(scaled[row][column] * vector[row] for row in range(size)) for column in range(size)]
    for column in range(size):
        normal[column][column] += _REGULARIZATION
        right_side[column] += _REGULARIZATION * target[column] * norms[column]
    solution = _gauss(normal, right_side)
    return [value / norm for value, norm in zip(solution, norms, strict=True)]


def _reduced(intervals: tuple[NasaPolynomialInterval, ...], temperature_k: float,
             log_temperature: float) -> tuple[float, float, float]:
    """Return ``(cp/R, H/RT, S/R)`` from the interval owning ``T`` (nearest end outside).

    The NASA Glenn forms are those of NASA/TP-2002-211556, as implemented in
    :class:`~rocket_propulsion.thermochemistry.nasa_polynomials.NasaPolynomialInterval`.
    """
    chosen = intervals[-1]
    if temperature_k < intervals[0].minimum_temperature_k:
        chosen = intervals[0]
    else:
        for interval in intervals:
            if temperature_k < interval.maximum_temperature_k:
                chosen = interval
                break
    a1, a2, a3, a4, a5, a6, a7, b1, b2 = chosen.coefficients
    t = temperature_k
    inverse = 1.0 / t
    inverse2 = inverse * inverse
    t2 = t * t
    t3 = t2 * t
    t4 = t3 * t
    cp = a1 * inverse2 + a2 * inverse + a3 + a4 * t + a5 * t2 + a6 * t3 + a7 * t4
    enthalpy = (-a1 * inverse2 + a2 * log_temperature * inverse + a3 + a4 * t / 2.0
                + a5 * t2 / 3.0 + a6 * t3 / 4.0 + a7 * t4 / 5.0 + b1 * inverse)
    entropy = (-a1 * inverse2 / 2.0 - a2 * inverse + a3 * log_temperature + a4 * t
               + a5 * t2 / 2.0 + a6 * t3 / 3.0 + a7 * t4 / 4.0 + b2)
    return cp, enthalpy, entropy


@dataclass(slots=True)
class _System:
    """Element inventory and candidate products of one problem, per kg of mixture."""

    elements: tuple[str, ...]
    b0: list[float]
    mass_kg: float
    gases: tuple[Species, ...]
    gas_atoms: tuple[tuple[tuple[int, float], ...], ...]
    gas_dense: tuple[tuple[float, ...], ...]
    condensed: tuple[CondensedSpecies, ...]
    condensed_dense: tuple[tuple[float, ...], ...]
    log_pressure: float
    mode: str
    temperature_k: float | None
    target_h_over_r: float | None
    target_s_over_r: float | None

    def fixed_temperature_k(self) -> float:
        """Prescribed temperature of a TP problem, K."""
        if self.temperature_k is None:
            raise InputError("TP requires a temperature.", field="temperature_k")
        return self.temperature_k

    def enthalpy_target(self) -> float:
        """Target ``H/R`` of an HP problem, K mol/kg."""
        if self.target_h_over_r is None:
            raise InputError("HP requires a target enthalpy.", field="target_enthalpy_j")
        return self.target_h_over_r

    def entropy_target(self) -> float:
        """Target ``S/R`` of an SP problem, mol/kg."""
        if self.target_s_over_r is None:
            raise InputError("SP requires a target entropy.", field="target_entropy_j_k")
        return self.target_s_over_r


@dataclass(slots=True)
class _Iterate:
    temperature_k: float
    log_total: float
    log_amounts: list[float]
    active: list[int]
    condensed_amounts: dict[int, float]


class CondensedEquilibriumProvider:
    """Gibbs-energy equilibrium of ideal gases and pure condensed phases.

    Attributes
    ----------
    name, version : str
        Provider identification copied into every result.
    composition_tolerance : float
        Bound on ``n_j |Delta ln n_j| / sum n``, ``|Delta n_c| / sum n`` and
        ``n |Delta ln n| / sum n`` (RP-1311 eq. (3.5) uses 0.5e-5).
    element_tolerance : float
        Bound on ``|b0_i - b_i| / max b0`` (eq. (3.6a) uses 1e-6).
    temperature_tolerance : float
        Bound on ``|Delta ln T|`` (eq. (3.6b) uses 1e-4).
    state_tolerance : float
        Bound on the relative closure of the assigned entropy (SP; eq. (3.6c)
        uses an absolute 0.5e-4 in kg-mol units) or enthalpy (HP).
    maximum_iterations : int
        Newton iterations allowed per pass (CEA uses ``50 + NS/2``).

    Notes
    -----
    Assumptions [1]_: ideal gases, pure condensed phases occupying no volume
    (section 2.2), no ions, NASA ``thermo.inp`` data between 200 and 6000 K.
    The algorithm is RP-1311 Part I chapters 2-3; the tolerances above are
    tighter than CEA's so that comparisons with CEA are dominated by CEA's
    own convergence. Selection among several condensed candidates uses the
    most negative eq. (3.7) value per mole. Gas polynomials are evaluated with
    the nearest interval when the temperature lies outside a species' range;
    a converged state in which such a species exceeds a mole fraction of 1e-6
    is rejected. A singular Newton matrix (section 3.6, eq. (3.11)) is solved
    by a column-scaled Tikhonov least-squares step that keeps the
    undetermined direction of the element potentials, instead of CEA's
    reset and component switch; states with gas below 1e-8 of all moles
    raise ``DomainError`` (code ``no_gas_phase``).

    References
    ----------
    .. [1] NASA-RP1311-P1-CONDENSED. Gordon and McBride, NASA RP-1311 Part I
           (1994). https://ntrs.nasa.gov/citations/19950013764
    """

    name = "local-gibbs-condensed"
    version = "0.2-nasa2002-experimental"
    composition_tolerance = 1.0e-10
    element_tolerance = 1.0e-10
    temperature_tolerance = 1.0e-10
    state_tolerance = 1.0e-10
    maximum_iterations = 300

    def __init__(self) -> None:
        self._warm: dict[tuple, CondensedEquilibriumState] = {}

    # ------------------------------------------------------------------ API
    def solve(self, problem: EquilibriumProblem | CondensedEquilibriumProblem) -> EquilibriumResult:
        """Solve a TP, HP or SP problem and return the provider-neutral result.

        Parameters
        ----------
        problem : EquilibriumProblem or CondensedEquilibriumProblem
            ``EquilibriumProblem`` instances of type TP or HP are accepted so
            that the provider plugs into
            :func:`~rocket_propulsion.thermochemistry.rocket_performance.solve_rocket_performance`.

        Returns
        -------
        EquilibriumResult
            See :class:`CondensedEquilibriumSolution` for the meaning of the
            fields when condensed phases are present.

        Raises
        ------
        InputError
            For an unsupported schema or problem type (UV) or missing state.
        DomainError
            For states outside the data range.
        ConvergenceError
            If the iteration fails.

        References
        ----------
        .. [1] NASA-RP1311-P1-CONDENSED. Gordon and McBride, NASA RP-1311 Part I
               (1994). https://ntrs.nasa.gov/citations/19950013764
        """
        return self.solve_detailed(problem).result

    def solve_detailed(
        self,
        problem: EquilibriumProblem | CondensedEquilibriumProblem,
        *,
        initial: CondensedEquilibriumState | None = None,
    ) -> CondensedEquilibriumSolution:
        """Solve a problem and return the condensed-phase detail.

        Parameters
        ----------
        problem : EquilibriumProblem or CondensedEquilibriumProblem
            TP, HP or SP problem.
        initial : CondensedEquilibriumState, optional
            Warm start (composition, included condensed phases, temperature).
            Without it the provider uses the last state it converged for the
            same reactant proportions, and otherwise the initial estimates of
            RP-1311 section 3.1 [1]_ (``n = 0.1`` kg-mol/kg spread evenly over
            the gases, no condensed species, 3800 K).

        Returns
        -------
        CondensedEquilibriumSolution
            Result, condensed amounts, element potentials, eq. (3.7) test
            values and the warm-start state.

        Raises
        ------
        InputError, DomainError, ConvergenceError
            As for :meth:`solve`.

        Notes
        -----
        A warm start that fails to converge is retried from the section 3.1
        estimates before an error is raised. Mixture properties follow
        RP-1311 eqs. (2.3), (2.14) and (2.16)-(2.17) [1]_: ``M`` is mass per
        gas mole, ``H = sum_j n_j H_j`` and ``S`` adds ``-R ln(n_j/n) - R ln P``
        for gases only. The reported ``cp`` is the frozen ``sum_j n_j cp_j``
        over all species per unit mass and ``cv = cp - n R`` (eq. (2.70) with
        the frozen derivatives ``(d ln V/d ln T)_p = 1`` and
        ``(d ln V/d ln p)_T = -1`` of eq. (2.1)).

        References
        ----------
        .. [1] NASA-RP1311-P1-CONDENSED. Gordon and McBride, NASA RP-1311 Part I
               (1994). https://ntrs.nasa.gov/citations/19950013764
        """
        if problem.schema_version != "1.0":
            raise InputError("Unsupported equilibrium problem schema.", code="unsupported_schema",
                             field="schema_version")
        system = self._system(problem)
        signature = (
            system.elements,
            tuple(item.species for item in problem.reactants),
            tuple(round(item.moles / sum(r.moles for r in problem.reactants), 12)
                  for item in problem.reactants),
            problem.candidate_species,
        )
        start = initial if initial is not None else self._warm.get(signature)
        try:
            iterate, events, iterations, potentials = self._solve_system(system, start)
        except DomainError:  # includes ConvergenceError
            if start is None:
                raise
            iterate, events, iterations, potentials = self._solve_system(system, None)
        solution = self._build(problem, system, iterate, events, iterations, potentials)
        self._warm[signature] = solution.state
        return solution

    # ------------------------------------------------------------ set-up
    def _system(self, problem: EquilibriumProblem | CondensedEquilibriumProblem) -> _System:
        kind = problem.problem_type.strip().upper()
        if kind not in ("TP", "HP", "SP"):
            raise InputError("Condensed equilibrium supports TP, HP and SP problems.",
                             code="unsupported_problem_type", field="problem_type")
        pressure = problem.pressure_pa
        if pressure is None or not isfinite(pressure) or pressure <= 0.0:
            raise DomainError("Equilibrium pressure must be finite and positive.", field="pressure_pa")
        catalog = condensed_catalog()
        mass = 0.0
        inventory: dict[str, float] = {}
        enthalpy = 0.0
        reactant_temperature = getattr(problem, "reactant_temperature_k", None)
        for reactant in problem.reactants:
            record = get_condensed_reactant(reactant.species)
            mass += reactant.moles * record.molecular_mass_kg_mol
            for element, count in record.formula:
                inventory[element] = inventory.get(element, 0.0) + reactant.moles * count
            if kind == "HP":
                if isinstance(record, CondensedReactant):
                    enthalpy += reactant.moles * record.enthalpy_j_mol
                else:
                    if reactant_temperature is None:
                        raise InputError("HP problems with gas reactants need reactant_temperature_k.",
                                         field="reactant_temperature_k")
                    enthalpy += reactant.moles * record.properties(reactant_temperature).enthalpy_molar_j_mol
        elements = tuple(sorted(inventory))
        allowed = set(elements)
        if problem.candidate_species is None:
            gases = tuple(item for item in catalog.gases.values()
                          if set(dict(item.formula)) <= allowed)
            condensed = tuple(item for item in catalog.condensed.values()
                              if set(dict(item.formula)) <= allowed)
        else:
            gas_list, condensed_list = [], []
            for key in problem.candidate_species:
                normalized = key.strip().upper()
                if normalized in catalog.condensed:
                    condensed_list.append(catalog.condensed[normalized])
                elif normalized in catalog.gases:
                    gas_list.append(catalog.gases[normalized])
                else:
                    raise InputError(f"Unknown product species: {key!r}.", code="unknown_species",
                                     field="candidate_species")
            gases = tuple(item for item in gas_list if set(dict(item.formula)) <= allowed)
            condensed = tuple(item for item in condensed_list if set(dict(item.formula)) <= allowed)
        represented = {element for item in gases for element, _ in item.formula}
        missing = sorted(allowed - represented)
        if missing:
            raise DomainError("Gas candidates do not represent every reactant element.",
                              code="incomplete_element_basis", details={"missing_elements": missing})
        index = {element: position for position, element in enumerate(elements)}

        def sparse(formula: tuple[tuple[str, float], ...]) -> tuple[tuple[int, float], ...]:
            return tuple((index[element], count) for element, count in formula if count != 0.0)

        def dense(formula: tuple[tuple[str, float], ...]) -> tuple[float, ...]:
            values = dict(formula)
            return tuple(values.get(element, 0.0) for element in elements)

        temperature = None
        target_h = None
        target_s = None
        if kind == "TP":
            temperature = problem.temperature_k
            if temperature is None or not isfinite(temperature):
                raise InputError("TP requires a temperature.", field="temperature_k")
            if not _MINIMUM_TEMPERATURE_K <= temperature <= _MAXIMUM_TEMPERATURE_K:
                raise DomainError("Temperature must lie between 200 and 6000 K.", field="temperature_k")
        elif kind == "HP":
            target = problem.target_enthalpy_j if problem.target_enthalpy_j is not None else enthalpy
            target_h = target / mass / UNIVERSAL_GAS_CONSTANT_J_MOL_K
        else:
            entropy = getattr(problem, "target_entropy_j_k", None)
            if entropy is None or not isfinite(entropy):
                raise InputError("SP requires a target entropy.", field="target_entropy_j_k")
            target_s = entropy / mass / UNIVERSAL_GAS_CONSTANT_J_MOL_K
        return _System(
            elements=elements,
            b0=[inventory[element] / mass for element in elements],
            mass_kg=mass,
            gases=gases,
            gas_atoms=tuple(sparse(item.formula) for item in gases),
            gas_dense=tuple(dense(item.formula) for item in gases),
            condensed=condensed,
            condensed_dense=tuple(dense(item.formula) for item in condensed),
            log_pressure=log(pressure / REFERENCE_PRESSURE_PA),
            mode=kind,
            temperature_k=temperature,
            target_h_over_r=target_h,
            target_s_over_r=target_s,
        )

    def _initial(self, system: _System, start: CondensedEquilibriumState | None) -> _Iterate:
        temperature = system.fixed_temperature_k() if system.mode == "TP" else _INITIAL_TEMPERATURE_K
        keys = tuple(item.key for item in system.gases)
        if start is not None and start.gas_keys == keys:
            condensed_index = {item.key: position for position, item in enumerate(system.condensed)}
            active = [condensed_index[key] for key, _ in start.condensed if key in condensed_index]
            amounts = {condensed_index[key]: value for key, value in start.condensed
                       if key in condensed_index}
            return _Iterate(
                temperature_k=temperature if system.mode == "TP" else start.temperature_k,
                log_total=start.log_gas_total,
                log_amounts=list(start.log_gas_amounts),
                active=active,
                condensed_amounts=amounts,
            )
        count = len(system.gases)
        return _Iterate(temperature, log(_INITIAL_TOTAL_MOL_KG),
                        [log(_INITIAL_TOTAL_MOL_KG / count)] * count, [], {})

    # ------------------------------------------------------------ Newton
    def _newton(self, system: _System, state: _Iterate) -> tuple[int, list[float]]:
        """Converge one fixed set of condensed species (RP-1311 eqs. (2.18)-(2.28), (3.1)-(3.6))."""
        elements = len(system.elements)
        variable_t = system.mode != "TP"
        previous_potentials = [0.0] * elements
        for iteration in range(1, self.maximum_iterations + 1):
            temperature = state.temperature_k
            log_t = log(temperature)
            log_p = system.log_pressure
            active = state.active
            ncond = len(active)
            gas_h, gas_s0, gas_cp = [], [], []
            for species in system.gases:
                cp, h, s = _reduced(species.intervals, temperature, log_t)
                gas_cp.append(cp)
                gas_h.append(h)
                gas_s0.append(s)
            cond_props = [_reduced(system.condensed[c].intervals, temperature, log_t) for c in active]
            amounts = [exp(value) for value in state.log_amounts]
            log_n = state.log_total
            total_gas = exp(log_n)
            mu = [gas_h[j] - gas_s0[j] + state.log_amounts[j] - log_n + log_p
                  for j in range(len(amounts))]
            size = elements + ncond + 1 + (1 if variable_t else 0)
            matrix = [[0.0] * size for _ in range(size)]
            rhs = [0.0] * size
            n_col = elements + ncond
            t_col = n_col + 1
            b = [0.0] * elements
            sum_n = 0.0
            sum_nh = 0.0
            sum_nmu = 0.0
            sum_ncp = 0.0
            sum_nh2 = 0.0
            sum_nhmu = 0.0
            sum_ns = 0.0
            sum_nhs = 0.0
            sum_nsmu = 0.0
            row_an = [0.0] * elements
            row_anh = [0.0] * elements
            row_anmu = [0.0] * elements
            row_ans = [0.0] * elements
            for j, atoms in enumerate(system.gas_atoms):
                nj = amounts[j]
                if nj == 0.0:
                    continue
                hj = gas_h[j]
                muj = mu[j]
                sum_n += nj
                sum_nh += nj * hj
                sum_nmu += nj * muj
                sum_ncp += nj * gas_cp[j]
                sum_nh2 += nj * hj * hj
                sum_nhmu += nj * hj * muj
                if system.mode == "SP":
                    sj = gas_s0[j] - state.log_amounts[j] + log_n - log_p
                    sum_ns += nj * sj
                    sum_nhs += nj * hj * sj
                    sum_nsmu += nj * sj * muj
                for k, akj in atoms:
                    weight = akj * nj
                    b[k] += weight
                    row_an[k] += weight
                    row_anh[k] += weight * hj
                    row_anmu[k] += weight * muj
                    if system.mode == "SP":
                        row_ans[k] += weight * (gas_s0[j] - state.log_amounts[j] + log_n - log_p)
                    for i, aij in atoms:
                        matrix[k][i] += weight * aij
            for position, c in enumerate(active):
                amount = state.condensed_amounts[c]
                for k, akc in enumerate(system.condensed_dense[c]):
                    b[k] += akc * amount
            # Element rows, eq. (2.24).
            for k in range(elements):
                for position, c in enumerate(active):
                    matrix[k][elements + position] = system.condensed_dense[c][k]
                matrix[k][n_col] = row_an[k]
                if variable_t:
                    matrix[k][t_col] = row_anh[k]
                rhs[k] = system.b0[k] - b[k] + row_anmu[k]
            # Condensed rows, eq. (2.25).
            for position, c in enumerate(active):
                row = elements + position
                _cp_c, h_c, s_c = cond_props[position]
                for i, aic in enumerate(system.condensed_dense[c]):
                    matrix[row][i] = aic
                if variable_t:
                    matrix[row][t_col] = h_c
                rhs[row] = h_c - s_c
            # Total-moles row, eq. (2.26).
            for i in range(elements):
                matrix[n_col][i] = row_an[i]
            matrix[n_col][n_col] = sum_n - total_gas
            if variable_t:
                matrix[n_col][t_col] = sum_nh
            rhs[n_col] = total_gas - sum_n + sum_nmu
            cond_cp = sum(state.condensed_amounts[c] * cond_props[p][0] for p, c in enumerate(active))
            if system.mode == "HP":
                h_total = sum_nh + sum(state.condensed_amounts[c] * cond_props[p][1]
                                       for p, c in enumerate(active))
                for i in range(elements):
                    matrix[t_col][i] = row_anh[i]
                for position, _ in enumerate(active):
                    matrix[t_col][elements + position] = cond_props[position][1]
                matrix[t_col][n_col] = sum_nh
                matrix[t_col][t_col] = sum_ncp + cond_cp + sum_nh2
                rhs[t_col] = system.enthalpy_target() / temperature - h_total + sum_nhmu
            elif system.mode == "SP":
                s_total = sum_ns + sum(state.condensed_amounts[c] * cond_props[p][2]
                                       for p, c in enumerate(active))
                for i in range(elements):
                    matrix[t_col][i] = row_ans[i]
                for position, _ in enumerate(active):
                    matrix[t_col][elements + position] = cond_props[position][2]
                matrix[t_col][n_col] = sum_ns
                matrix[t_col][t_col] = sum_ncp + cond_cp + sum_nhs
                rhs[t_col] = system.entropy_target() - s_total + total_gas - sum_n + sum_nsmu
            try:
                correction = _gauss(matrix, rhs)
            except _Singular:
                condensed_total = sum(state.condensed_amounts[c] for c in active)
                if sum_n < _MINIMUM_GAS_FRACTION * (sum_n + condensed_total):
                    raise DomainError(
                        "The equilibrium has no gas phase at this pressure; the ideal-gas "
                        "formulation does not apply.", code="no_gas_phase",
                        details={"gas_mole_fraction": sum_n / (sum_n + condensed_total)},
                    ) from None
                # Section 3.6 singularities (e.g. eq. (3.11), all gas species
                # in one element ratio): keep the undetermined direction of
                # the previous iterate instead of CEA's reset/component switch.
                target = previous_potentials + [0.0] * (size - elements)
                try:
                    correction = _regularized(matrix, rhs, target)
                except _Singular:
                    raise ConvergenceError("Condensed equilibrium matrix is singular.",
                                           code="equilibrium_singular_jacobian",
                                           details={"condensed": [system.condensed[c].key for c in active]}
                                           ) from None
            potentials = correction[:elements]
            previous_potentials = potentials
            delta_cond = correction[elements:elements + ncond]
            delta_log_n = correction[n_col]
            delta_log_t = correction[t_col] if variable_t else 0.0
            delta = []
            for j, atoms in enumerate(system.gas_atoms):
                value = -mu[j] + delta_log_n + gas_h[j] * delta_log_t
                for i, aij in atoms:
                    value += aij * potentials[i]
                delta.append(value)
            # Convergence tests, eqs. (3.5)-(3.6c), on the corrections of this iterate.
            total_all = sum_n + sum(state.condensed_amounts[c] for c in active)
            converged = (
                all(amounts[j] * abs(delta[j]) <= self.composition_tolerance * total_all
                    for j in range(len(delta)))
                and all(abs(value) <= self.composition_tolerance * total_all for value in delta_cond)
                and total_gas * abs(delta_log_n) <= self.composition_tolerance * total_all
                and all(abs(system.b0[k] - b[k]) <= self.element_tolerance * max(system.b0)
                        for k in range(elements))
                and abs(delta_log_t) <= self.temperature_tolerance
            )
            if converged and system.mode == "SP":
                target_s = system.entropy_target()
                converged = abs(target_s - s_total) <= self.state_tolerance * max(
                    abs(target_s), 1.0)
            if converged and system.mode == "HP":
                target_h = system.enthalpy_target()
                converged = abs(target_h / temperature - h_total) <= (
                    self.state_tolerance * max(abs(target_h / temperature), 1.0))
            # Control factor, eqs. (3.1)-(3.3).
            largest = max(5.0 * abs(delta_log_t), 5.0 * abs(delta_log_n))
            lambda2 = 1.0
            for j, value in enumerate(delta):
                relative = state.log_amounts[j] - log_n
                if relative > -_SIZE:
                    largest = max(largest, abs(value))
                elif value >= 0.0 and value - delta_log_n > 0.0:
                    lambda2 = min(lambda2, abs((-relative - _LAMBDA2_LOG) / (value - delta_log_n)))
            lambda1 = 2.0 / largest if largest > 2.0 else 1.0
            factor = min(1.0, lambda1, lambda2)
            # Updates, eq. (3.4).
            state.log_amounts = [value + factor * change
                                 for value, change in zip(state.log_amounts, delta, strict=True)]
            for position, c in enumerate(active):
                state.condensed_amounts[c] += factor * delta_cond[position]
            state.log_total = log_n + factor * delta_log_n
            if variable_t:
                state.temperature_k = min(max(temperature * exp(factor * delta_log_t), 100.0), 7000.0)
            if converged:
                return iteration, potentials
        raise ConvergenceError("Condensed equilibrium did not converge.",
                               code="equilibrium_not_converged",
                               details={"iterations": self.maximum_iterations,
                                        "condensed": [system.condensed[c].key for c in state.active]})

    def _test_values(self, system: _System, state: _Iterate, potentials: list[float]) -> dict[int, float]:
        """Return eq. (3.7) values for candidates covering T, outside the active groups."""
        temperature = state.temperature_k
        log_t = log(temperature)
        groups = {system.condensed[c].phase_group for c in state.active}
        values = {}
        for c, species in enumerate(system.condensed):
            if c in state.active or species.phase_group in groups or not species.covers(temperature):
                continue
            _, h, s = _reduced(species.intervals, temperature, log_t)
            values[c] = (h - s) - sum(a * p for a, p in zip(system.condensed_dense[c], potentials,
                                                            strict=True))
        return values

    def _solve_system(self, system: _System, start: CondensedEquilibriumState | None):
        """Alternate Newton convergence with the condensed-phase tests of RP-1311 sections 3.4-3.5."""
        state = self._initial(system, start)
        events: list[str] = []
        total_iterations = 0
        seen: set[tuple[tuple[int, ...], int]] = set()
        variable_t = system.mode != "TP"
        for _ in range(4 * len(system.condensed) + 12):
            if system.mode == "TP":
                # With T assigned only phases covering T are meaningful (section 3.5).
                for c in list(state.active):
                    if not system.condensed[c].covers(state.temperature_k):
                        state.active.remove(c)
                        state.condensed_amounts.pop(c, None)
                        events.append(f"remove {system.condensed[c].key} (outside range)")
            iterations, potentials = self._newton(system, state)
            total_iterations += iterations
            temperature = state.temperature_k
            gas_total = sum(exp(value) for value in state.log_amounts)
            condensed_total = sum(state.condensed_amounts[c] for c in state.active)
            if gas_total < _MINIMUM_GAS_FRACTION * (gas_total + max(condensed_total, 0.0)):
                raise DomainError(
                    "The equilibrium has no gas phase at this pressure; the ideal-gas "
                    "formulation does not apply.", code="no_gas_phase",
                    details={"gas_mole_fraction": gas_total / (gas_total + condensed_total)},
                )
            # Negative condensed amounts are removed (section 3.4).
            negative = [c for c in state.active if state.condensed_amounts[c] <= 0.0]
            if negative:
                worst = min(negative, key=lambda c: state.condensed_amounts[c])
                state.active.remove(worst)
                state.condensed_amounts.pop(worst)
                events.append(f"remove {system.condensed[worst].key}")
                continue
            # Phase transitions of included species with T variable (section 3.5).
            changed = False
            if variable_t:
                for c in list(state.active):
                    species = system.condensed[c]
                    if species.covers(temperature):
                        continue
                    above = temperature > species.maximum_temperature_k
                    boundary = species.maximum_temperature_k if above else species.minimum_temperature_k
                    neighbour = next(
                        (d for d, other in enumerate(system.condensed)
                         if other.phase_group == species.phase_group and d != c
                         and (other.minimum_temperature_k == boundary if above
                              else other.maximum_temperature_k == boundary)),
                        None,
                    )
                    if neighbour is None:
                        state.active.remove(c)
                        state.condensed_amounts.pop(c)
                        events.append(f"remove {species.key} (outside range)")
                    elif neighbour in state.active:
                        continue
                    elif abs(temperature - boundary) > _TRANSITION_BAND_K:
                        amount = state.condensed_amounts.pop(c)
                        state.active[state.active.index(c)] = neighbour
                        state.condensed_amounts[neighbour] = amount
                        events.append(f"switch {species.key} to {system.condensed[neighbour].key}")
                    else:
                        state.active.append(neighbour)
                        state.condensed_amounts[neighbour] = 0.0
                        events.append(f"add {system.condensed[neighbour].key} at transition")
                    changed = True
                    break
            if changed:
                continue
            # Inclusion test, eq. (3.7): at most one new species per convergence.
            tests = self._test_values(system, state, potentials)
            candidate = min(tests, key=tests.__getitem__, default=None)
            if candidate is not None and tests[candidate] < 0.0:
                key = (tuple(sorted(state.active + [candidate])), round(temperature))
                if key in seen:
                    raise ConvergenceError("Condensed-phase selection is cycling.",
                                           code="condensed_phase_cycle",
                                           details={"species": system.condensed[candidate].key})
                seen.add(key)
                state.active.append(candidate)
                state.condensed_amounts[candidate] = 0.0
                events.append(f"include {system.condensed[candidate].key}")
                continue
            return state, events, total_iterations, potentials
        raise ConvergenceError("Condensed-phase selection did not settle.",
                               code="condensed_phase_cycle")

    # ------------------------------------------------------------ result
    def _build(self, problem, system: _System, state: _Iterate, events: list[str], iterations: int,
               potentials: list[float]) -> CondensedEquilibriumSolution:
        temperature = state.temperature_k
        if not _MINIMUM_TEMPERATURE_K <= temperature <= _MAXIMUM_TEMPERATURE_K:
            raise DomainError("Equilibrium temperature left the 200-6000 K data range.",
                              code="equilibrium_temperature_out_of_range",
                              details={"temperature_k": temperature})
        mass = system.mass_kg
        r = UNIVERSAL_GAS_CONSTANT_J_MOL_K
        log_t = log(temperature)
        gas_amounts = {species.key: exp(value) * mass
                       for species, value in zip(system.gases, state.log_amounts, strict=True)}
        condensed_amounts = {system.condensed[c].key: state.condensed_amounts[c] * mass
                             for c in state.active}
        gas_moles = sum(gas_amounts.values())
        total_moles = gas_moles + sum(condensed_amounts.values())
        enthalpy = 0.0
        entropy = 0.0
        cp_total = 0.0
        out_of_range_trace = []
        log_p = system.log_pressure
        for species in system.gases:
            amount = gas_amounts[species.key]
            if amount <= 0.0:
                continue
            cp, h, s0 = _reduced(species.intervals, temperature, log_t)
            enthalpy += amount * h * r * temperature
            entropy += amount * r * (s0 - log(amount / gas_moles) - log_p)
            cp_total += amount * cp * r
            if not species.minimum_temperature_k <= temperature <= species.maximum_temperature_k:
                if amount / total_moles >= _RANGE_CHECK_MOLE_FRACTION:
                    raise DomainError(
                        f"{species.key} is significant outside its data range.",
                        code="species_temperature_out_of_range",
                        details={"species": species.key, "temperature_k": temperature},
                    )
                out_of_range_trace.append(species.key)
        by_key = {item.key: item for item in system.condensed}
        for key, amount in condensed_amounts.items():
            phase = by_key[key]
            cp, h, s0 = _reduced(phase.intervals, temperature, log_t)
            enthalpy += amount * h * r * temperature
            entropy += amount * s0 * r
            cp_total += amount * cp * r
            outside = max(phase.minimum_temperature_k - temperature,
                          temperature - phase.maximum_temperature_k)
            coexisting = any(by_key[other].phase_group == phase.phase_group
                             and by_key[other].covers(temperature)
                             for other in condensed_amounts if other != key)
            # Two phases of one substance pin T to their transition; the data
            # meet there to |dG/RT| < 3e-8 (sync_nasa_condensed.py), i.e. a few
            # 1e-5 K for alumina, so 0.01 K is accepted.
            if outside > 0.0 and not (coexisting and outside <= 0.01):
                raise DomainError(f"{key} converged outside its data range.",
                                  code="species_temperature_out_of_range",
                                  details={"species": key, "temperature_k": temperature})
        molar_masses = {species.key: species.molecular_mass_kg_mol for species in system.gases}
        molar_masses.update({key: by_key[key].molecular_mass_kg_mol for key in condensed_amounts})
        product_mass = sum(amount * molar_masses[key] for key, amount in gas_amounts.items())
        product_mass += sum(amount * molar_masses[key] for key, amount in condensed_amounts.items())
        everything = {**gas_amounts, **condensed_amounts}
        retained = {key: value for key, value in everything.items() if value > 1.0e-14 * total_moles}
        formulas = {species.key: dict(species.formula) for species in system.gases}
        formulas.update({key: dict(by_key[key].formula) for key in condensed_amounts})
        element_residual = 0.0
        for position, element in enumerate(system.elements):
            target = system.b0[position] * mass
            calculated = sum(amount * formulas[key].get(element, 0.0) for key, amount in everything.items())
            element_residual = max(element_residual, abs(calculated - target) / target)
        energy_target = None
        energy_residual = None
        if system.mode == "HP":
            energy_target = system.enthalpy_target() * r * mass
            energy_residual = abs(enthalpy - energy_target) / max(abs(energy_target), 1.0)
        elif system.mode == "SP":
            target_s = system.entropy_target() * r * mass
            energy_residual = abs(entropy - target_s) / max(abs(target_s), 1.0e-12)
        molecular_mass = mass / gas_moles
        gas_constant = r / molecular_mass
        cp_mass = cp_total / mass
        cv_mass = cp_mass - gas_moles * r / mass
        species_amounts = tuple(
            EquilibriumSpeciesAmount(
                species=key,
                moles=value,
                mole_fraction=value / total_moles,
                mass_fraction=value * molar_masses[key] / product_mass,
            )
            for key, value in sorted(retained.items(), key=lambda item: -item[1])
        )
        warnings = [WarningMessage(
            code="experimental_equilibrium_provider",
            message=("The local condensed-phase equilibrium solver is experimental; validate "
                     "design results against NASA CEA."),
        )]
        if out_of_range_trace:
            warnings.append(WarningMessage(
                code="trace_species_extrapolated",
                message=(f"{len(out_of_range_trace)} trace species (mole fraction below 1e-6) were "
                         "evaluated outside their polynomial range."),
                severity="info",
            ))
        mole_sum = sum(item.mole_fraction for item in species_amounts)
        residuals = EquilibriumResiduals(
            maximum_element_relative=element_residual,
            mole_fraction_sum=abs(mole_sum - 1.0),
            energy_relative=energy_residual,
        )
        result = EquilibriumResult(
            schema_version="1.0",
            provider=self.name,
            provider_version=self.version,
            problem_type=system.mode,
            temperature_k=temperature,
            pressure_pa=REFERENCE_PRESSURE_PA * exp(log_p),
            total_moles=total_moles,
            total_mass_kg=product_mass,
            molecular_mass_kg_mol=molecular_mass,
            gas_constant_j_kg_k=gas_constant,
            cp_j_kg_k=cp_mass,
            cv_j_kg_k=cv_mass,
            gamma=cp_mass / cv_mass,
            enthalpy_j=enthalpy,
            internal_energy_j=enthalpy - gas_moles * r * temperature,
            entropy_j_k=entropy,
            species=species_amounts,
            element_moles=tuple((element, system.b0[i] * mass) for i, element in enumerate(system.elements)),
            residuals=residuals,
            convergence=EquilibriumConvergence(
                converged=True,
                composition_iterations=iterations,
                temperature_iterations=iterations if system.mode != "TP" else 0,
                pressure_iterations=0,
                maximum_residual=max(element_residual, residuals.mole_fraction_sum, energy_residual or 0.0),
            ),
            assumptions=(
                "ideal-gas mixture plus pure condensed phases of negligible volume",
                "no velocity or thermal lag between phases",
                "no ionized species",
                "NASA CEA thermo.inp data",
                "experimental local Gibbs minimizer (RP-1311 method)",
            ),
            warnings=tuple(warnings),
            source="Gordon and McBride, NASA RP-1311 Part I; NASA CEA thermo.inp",
            source_url="https://ntrs.nasa.gov/citations/19950013764",
        )
        condensed = tuple(
            CondensedPhaseAmount(species=key, phase_group=by_key[key].phase_group, moles=value,
                                 mass_fraction=value * by_key[key].molecular_mass_kg_mol / product_mass)
            for key, value in sorted(condensed_amounts.items())
        )
        tests = self._test_values(system, state, potentials)
        return CondensedEquilibriumSolution(
            result=result,
            gas_moles=gas_moles,
            condensed=condensed,
            condensed_mass_fraction=sum(item.mass_fraction for item in condensed),
            mixture_molecular_mass_kg_mol=product_mass / total_moles,
            element_potentials=tuple(zip(system.elements, potentials, strict=True)),
            inclusion_tests=tuple((system.condensed[c].key, value) for c, value in sorted(tests.items())),
            events=tuple(events),
            newton_iterations=iterations,
            state=CondensedEquilibriumState(
                gas_keys=tuple(species.key for species in system.gases),
                log_gas_amounts=tuple(state.log_amounts),
                log_gas_total=state.log_total,
                condensed=tuple((system.condensed[c].key, state.condensed_amounts[c]) for c in state.active),
                temperature_k=temperature,
            ),
        )


# ---------------------------------------------------------------- rocket
@dataclass(frozen=True, slots=True)
class CondensedStationDetail:
    """Two-phase detail of one rocket station.

    Attributes
    ----------
    name : {"chamber", "throat", "exit"}
        Station kind.
    area_ratio : float or None
        Local area over throat area; ``None`` for the chamber.
    temperature_k, pressure_pa : float
        Station state.
    condensed_mass_fraction : float
        Condensed mass over total mass.
    condensed_species : tuple of (str, float)
        Condensed phases and their mass fractions.
    molecular_mass_kg_mol, mixture_molecular_mass_kg_mol : float
        ``M`` (mass per gas mole) and ``MW`` (mass per mole of all species) of
        RP-1311 eqs. (2.3) and (2.4a).
    events : tuple of str
        Phase changes the solver made at this station.

    References
    ----------
    .. [1] NASA-RP1311-P1-CONDENSED. Gordon and McBride, NASA RP-1311 Part I
           (1994). https://ntrs.nasa.gov/citations/19950013764
    """

    name: str
    area_ratio: float | None
    temperature_k: float
    pressure_pa: float
    condensed_mass_fraction: float
    condensed_species: tuple[tuple[str, float], ...]
    molecular_mass_kg_mol: float
    mixture_molecular_mass_kg_mol: float
    events: tuple[str, ...]


@dataclass(frozen=True, slots=True)
class CondensedRocketPerformance:
    """Two-phase infinite-area-combustor rocket result.

    Attributes
    ----------
    performance : RocketPerformance
        Chamber, throat, exits, ``c*``, vacuum Isp and ``C_F`` in the format of
        :func:`~rocket_propulsion.thermochemistry.rocket_performance.solve_rocket_performance`;
        station mole fractions include condensed species.
    stations : tuple of CondensedStationDetail
        Chamber, throat and exits in that order.
    """

    performance: RocketPerformance
    stations: tuple[CondensedStationDetail, ...]


def _detail(name: str, ratio: float | None, solution: CondensedEquilibriumSolution) -> CondensedStationDetail:
    result = solution.result
    return CondensedStationDetail(
        name=name,
        area_ratio=ratio,
        temperature_k=result.temperature_k,
        pressure_pa=result.pressure_pa,
        condensed_mass_fraction=solution.condensed_mass_fraction,
        condensed_species=tuple((item.species, item.mass_fraction) for item in solution.condensed),
        molecular_mass_kg_mol=result.molecular_mass_kg_mol,
        mixture_molecular_mass_kg_mol=solution.mixture_molecular_mass_kg_mol,
        events=solution.events,
    )


def solve_condensed_rocket_performance(
    *,
    reactants: Mapping[str, float],
    chamber_pressure_pa: float,
    area_ratios: Sequence[float] = (10.0, 25.0, 50.0),
    reactant_temperature_k: float = 298.15,
    provider: CondensedEquilibriumProvider | None = None,
) -> CondensedRocketPerformance:
    """Solve the infinite-area-combustor rocket problem with condensed products [1]_.

    Parameters
    ----------
    reactants : mapping of str to float
        Reactant amounts in mol; any key of :func:`get_condensed_reactant`
        (for example ``"AL(cr)"``, ``"NH4CLO4(I)"``, ``"O2(L)"``, ``"RP-1"``).
    chamber_pressure_pa : float
        Chamber pressure in Pa.
    area_ratios : sequence of float, default (10, 25, 50)
        Supersonic exit-to-throat area ratios, each greater than one.
    reactant_temperature_k : float, default 298.15
        Temperature of gas reactants.
    provider : CondensedEquilibriumProvider, optional
        Solver instance; a new one by default.

    Returns
    -------
    CondensedRocketPerformance
        Rocket result and per-station condensed detail.

    Raises
    ------
    DomainError
        For a non-positive pressure, an area ratio not above one, or a
        station that cannot be bracketed.
    ConvergenceError
        If an equilibrium or a station root fails to converge.

    Notes
    -----
    Assumptions of RP-1311 section 6.1 [1]_: one-dimensional flow, complete
    adiabatic combustion, isentropic expansion, homogeneous mixing, ideal-gas
    law, zero temperature and velocity lags between condensed and gaseous
    species, shifting equilibrium. Every station is an SP equilibrium at the
    chamber entropy, warm-started from the nearest station already solved.
    Velocity follows eq. (6.5), ``u = sqrt(2 (h_inf - h))``; density is
    ``p / (n R T)`` with ``n`` the gas moles per unit mass (eq. (2.1b)). The
    throat is the pressure of maximum ``rho u`` (equivalent to RP-1311's
    sonic-throat condition of section 6.3.3), found by golden-section search
    on ``ln p`` between ``p_c/2.5`` and ``p_c/1.2`` to 1e-5; exits solve the
    continuity relation (6.1) for the area ratio on the supersonic branch.
    ``c* = p_c / (rho u)_t`` and the vacuum impulse ``u_e + p_e A_e/mdot``
    (eq. (6.7) with zero ambient pressure). The root helpers are those of
    :mod:`~rocket_propulsion.thermochemistry.rocket_performance`.

    References
    ----------
    .. [1] NASA-RP1311-P1-CONDENSED. Gordon and McBride, NASA RP-1311 Part I
           (1994). https://ntrs.nasa.gov/citations/19950013764
    """

    if not isfinite(chamber_pressure_pa) or chamber_pressure_pa <= 0.0:
        raise DomainError("Chamber pressure must be positive.", field="chamber_pressure_pa")
    ratios = tuple(float(value) for value in area_ratios)
    if not ratios or any(not isfinite(value) or value <= 1.0 for value in ratios):
        raise DomainError("Rocket area ratios must all exceed one.", field="area_ratios")
    active = provider or CondensedEquilibriumProvider()
    pc = chamber_pressure_pa
    chamber = active.solve_detailed(CondensedEquilibriumProblem.hp(
        reactants, pressure_pa=pc, reactant_temperature_k=reactant_temperature_k))
    entropy = chamber.result.entropy_j_k
    mass = chamber.result.total_mass_kg
    h_chamber = chamber.result.enthalpy_j / mass
    solved: dict[float, CondensedEquilibriumSolution] = {log(pc): chamber}

    def station(log_pressure: float) -> CondensedEquilibriumSolution:
        if log_pressure in solved:
            return solved[log_pressure]
        nearest = min(solved, key=lambda value: abs(value - log_pressure))
        solution = active.solve_detailed(
            CondensedEquilibriumProblem.sp(reactants, pressure_pa=exp(log_pressure),
                                           target_entropy_j_k=entropy),
            initial=solved[nearest].state,
        )
        solved[log_pressure] = solution
        return solution

    def velocity(solution: CondensedEquilibriumSolution) -> float:
        drop = h_chamber - solution.result.enthalpy_j / mass
        if drop <= 0.0:
            raise DomainError("Isentropic expansion produced no enthalpy drop.",
                              code="invalid_rocket_expansion")
        return sqrt(2.0 * drop)

    def density(solution: CondensedEquilibriumSolution) -> float:
        result = solution.result
        return result.pressure_pa / (result.gas_constant_j_kg_k * result.temperature_k)

    def flux(log_pressure: float) -> float:
        solution = station(log_pressure)
        return density(solution) * velocity(solution)

    log_throat = _maximize_bounded(flux, log(pc / 2.5), log(pc / 1.2), 1.0e-5)
    throat = station(log_throat)
    throat_velocity = velocity(throat)
    throat_flux = density(throat) * throat_velocity

    def make_station(name: str, ratio: float | None, solution: CondensedEquilibriumSolution,
                     speed: float) -> RocketStation:
        result = solution.result
        return RocketStation(
            name=name, area_ratio=ratio, pressure_pa=result.pressure_pa,
            temperature_k=result.temperature_k, density_kg_m3=density(solution),
            molecular_mass_kg_mol=result.molecular_mass_kg_mol,
            enthalpy_j_kg=result.enthalpy_j / mass, entropy_j_kg_k=result.entropy_j_k / mass,
            velocity_m_s=speed,
            mole_fractions=tuple((item.species, item.mole_fraction) for item in result.species
                                 if item.mole_fraction >= 1.0e-5),
        )

    exits: list[RocketExit] = []
    details = [_detail("chamber", None, chamber), _detail("throat", 1.0, throat)]
    for ratio in ratios:
        def residual(log_pressure: float, ratio: float = ratio) -> float:
            return log(throat_flux / flux(log_pressure)) - log(ratio)

        guess = log_throat - log(ratio) * chamber.result.gamma
        lo, hi, f_lo, f_hi = _bracket(residual, guess, log_throat - 40.0, log_throat - 1.0e-3)
        root = _illinois(residual, lo, hi, f_lo, f_hi, 1.0e-12)
        solution = station(root)
        speed = velocity(solution)
        vacuum = speed + solution.result.pressure_pa * ratio / throat_flux
        exits.append(RocketExit(
            area_ratio=ratio,
            station=make_station("exit", ratio, solution, speed),
            effective_exhaust_velocity_m_s=speed,
            vacuum_effective_exhaust_velocity_m_s=vacuum,
            vacuum_specific_impulse_s=vacuum / STANDARD_GRAVITY_M_S2,
            vacuum_thrust_coefficient=(throat_flux * speed + solution.result.pressure_pa * ratio) / pc,
        ))
        details.append(_detail("exit", ratio, solution))
    performance = RocketPerformance(
        provider=active.name,
        provider_version=active.version,
        chemistry="equilibrium",
        reactants=tuple(sorted((key.strip().upper(), float(value)) for key, value in reactants.items())),
        chamber_pressure_pa=pc,
        chamber=make_station("chamber", None, chamber, 0.0),
        throat=make_station("throat", 1.0, throat, throat_velocity),
        exits=tuple(exits),
        characteristic_velocity_m_s=pc / throat_flux,
        warnings=(WarningMessage(
            code="ideal_two_phase_rocket_assumptions",
            message=("Infinite-area combustor, adiabatic one-dimensional isentropic flow, shifting "
                     "equilibrium with condensed products in velocity and thermal equilibrium with "
                     "the gas; no kinetics, friction, divergence or particle lag."),
        ),),
    )
    return CondensedRocketPerformance(performance=performance, stations=tuple(details))
