"""Chemical-equilibrium provider contract and an experimental local solver.

The local provider minimizes ideal-gas Gibbs energy subject to element
conservation [1]_ using NASA Glenn species data [2]_. It is limited to the
bundled gas species and 200--6000 K. Bounded solver comparisons against NASA
CEA are recorded in ``tests/test_cea_gate.py``; shared thermodynamic data do
not independently validate the catalogue or hardware. The provider protocol also lets the NASA CEA adapter of
:mod:`.cea` be selected without coupling the rest of the application to an
external executable.

References
----------
.. [1] NASA RP-1311, "Computer Program for Calculation of Complex Chemical
       Equilibrium Compositions and Applications, Part I: Analysis".
       https://ntrs.nasa.gov/citations/19950013764
.. [2] NASA/TP-2002-211556. NASA Glenn, "NASA Glenn Coefficients for
       Calculating Thermodynamic Properties of Individual Species".
       https://ntrs.nasa.gov/citations/20020085330
"""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass, replace
from itertools import combinations
from math import exp, isfinite, log
from operator import mul
from typing import Protocol, runtime_checkable

from rocket_propulsion.core.errors import ConvergenceError, DomainError, InputError
from rocket_propulsion.core.metadata import WarningMessage

from .mixtures import REFERENCE_PRESSURE_PA, Mixture
from .nasa_polynomials import UNIVERSAL_GAS_CONSTANT_J_MOL_K
from .species import CondensedReactant, Species, get_reactant, get_species, species_catalog

#: Library bound, not a source value: the gas phase is the ideal-gas mixture of
#: :mod:`.mixtures`, an assumption that fails long before 1 GPa (rocket
#: chambers run below about 30 MPa). NASA CEA itself imposes no bound; above
#: this pressure solves are refused instead of returning ideal-gas numbers for
#: states that cannot be ideal gases.
MAXIMUM_IDEAL_GAS_PRESSURE_PA = 1.0e9


@dataclass(frozen=True, slots=True)
class ReactantAmount:
    """A positive molar amount of one bundled species or condensed reactant.

    Attributes
    ----------
    species : str
        Upper-case gas species key (``"O2"``) or condensed reactant key
        (``"O2(L)"``, ``"RP-1"``).
    moles : float
        Finite, positive amount in mol.
    """

    species: str
    moles: float


@dataclass(frozen=True, slots=True)
class EquilibriumProblem:
    """A versioned TP, HP, or UV ideal-gas equilibrium problem.

    Use the :meth:`tp`, :meth:`hp` and :meth:`uv` constructors; they validate
    and sort the reactants.

    Attributes
    ----------
    problem_type : str
        ``"TP"`` (assigned temperature and pressure), ``"HP"`` (assigned
        enthalpy and pressure) or ``"UV"`` (assigned internal energy and
        volume) [1]_.
    reactants : tuple of ReactantAmount
        Reactant amounts sorted by species key.
    pressure_pa : float or None
        Pressure for TP and HP problems.
    temperature_k : float or None
        Temperature for TP problems.
    reactant_temperature_k : float or None
        Reactant temperature defining the default HP/UV energy target.
    target_enthalpy_j, target_internal_energy_j : float or None
        Explicit energy targets in J; when ``None`` the reactant enthalpy or
        internal energy at ``reactant_temperature_k`` is used.
    volume_m3 : float or None
        Volume for UV problems.
    candidate_species : tuple of str or None
        Product species to consider; ``None`` uses every bundled species whose
        elements are present in the reactants.
    schema_version : str, default "1.0"
        Problem schema version.

    References
    ----------
    .. [1] NASA RP-1311, "Computer Program for Calculation of Complex Chemical
           Equilibrium Compositions and Applications, Part I: Analysis".
           https://ntrs.nasa.gov/citations/19950013764
    """

    problem_type: str
    reactants: tuple[ReactantAmount, ...]
    pressure_pa: float | None = None
    temperature_k: float | None = None
    reactant_temperature_k: float | None = None
    target_enthalpy_j: float | None = None
    target_internal_energy_j: float | None = None
    volume_m3: float | None = None
    candidate_species: tuple[str, ...] | None = None
    schema_version: str = "1.0"

    @classmethod
    def tp(
        cls,
        reactants: Mapping[str, float],
        *,
        temperature_k: float,
        pressure_pa: float,
        candidate_species: tuple[str, ...] | None = None,
    ) -> EquilibriumProblem:
        """Build an assigned-temperature-and-pressure problem [1]_.

        Parameters
        ----------
        reactants : mapping of str to float
            Positive molar amount per species key.
        temperature_k : float
            Equilibrium temperature in K.
        pressure_pa : float
            Equilibrium pressure in Pa.
        candidate_species : tuple of str, optional
            Product species to consider.

        Returns
        -------
        EquilibriumProblem
            Problem with ``problem_type="TP"``.

        Raises
        ------
        InputError
            If no reactant is given, an amount is not numeric or a species is
            unknown.
        DomainError
            If an amount is not finite and positive.

        References
        ----------
        .. [1] NASA RP-1311, "Computer Program for Calculation of Complex Chemical
               Equilibrium Compositions and Applications, Part I: Analysis".
               https://ntrs.nasa.gov/citations/19950013764
        """

        return cls(
            problem_type="TP",
            reactants=_reactants(reactants),
            temperature_k=temperature_k,
            pressure_pa=pressure_pa,
            candidate_species=candidate_species,
        )

    @classmethod
    def hp(
        cls,
        reactants: Mapping[str, float],
        *,
        pressure_pa: float,
        reactant_temperature_k: float = 298.15,
        target_enthalpy_j: float | None = None,
        candidate_species: tuple[str, ...] | None = None,
    ) -> EquilibriumProblem:
        """Build an assigned-enthalpy-and-pressure combustion problem [1]_.

        Parameters
        ----------
        reactants : mapping of str to float
            Positive molar amount per species key.
        pressure_pa : float
            Combustion pressure in Pa.
        reactant_temperature_k : float, default 298.15
            Reactant temperature in K that sets the default enthalpy target.
        target_enthalpy_j : float, optional
            Explicit total enthalpy target in J.
        candidate_species : tuple of str, optional
            Product species to consider.

        Returns
        -------
        EquilibriumProblem
            Problem with ``problem_type="HP"``.

        Raises
        ------
        InputError
            If no reactant is given, an amount is not numeric or a species is
            unknown.
        DomainError
            If an amount is not finite and positive.

        References
        ----------
        .. [1] NASA RP-1311, "Computer Program for Calculation of Complex Chemical
               Equilibrium Compositions and Applications, Part I: Analysis".
               https://ntrs.nasa.gov/citations/19950013764
        """

        return cls(
            problem_type="HP",
            reactants=_reactants(reactants),
            pressure_pa=pressure_pa,
            reactant_temperature_k=reactant_temperature_k,
            target_enthalpy_j=target_enthalpy_j,
            candidate_species=candidate_species,
        )

    @classmethod
    def uv(
        cls,
        reactants: Mapping[str, float],
        *,
        volume_m3: float,
        reactant_temperature_k: float = 298.15,
        target_internal_energy_j: float | None = None,
        candidate_species: tuple[str, ...] | None = None,
    ) -> EquilibriumProblem:
        """Build an assigned-internal-energy-and-volume problem [1]_.

        Parameters
        ----------
        reactants : mapping of str to float
            Positive molar amount per species key.
        volume_m3 : float
            Closed volume in m^3.
        reactant_temperature_k : float, default 298.15
            Reactant temperature in K that sets the default energy target.
        target_internal_energy_j : float, optional
            Explicit total internal-energy target in J.
        candidate_species : tuple of str, optional
            Product species to consider.

        Returns
        -------
        EquilibriumProblem
            Problem with ``problem_type="UV"``.

        Raises
        ------
        InputError
            If no reactant is given, an amount is not numeric or a species is
            unknown.
        DomainError
            If an amount is not finite and positive.

        References
        ----------
        .. [1] NASA RP-1311, "Computer Program for Calculation of Complex Chemical
               Equilibrium Compositions and Applications, Part I: Analysis".
               https://ntrs.nasa.gov/citations/19950013764
        """

        return cls(
            problem_type="UV",
            reactants=_reactants(reactants),
            reactant_temperature_k=reactant_temperature_k,
            target_internal_energy_j=target_internal_energy_j,
            volume_m3=volume_m3,
            candidate_species=candidate_species,
        )


@dataclass(frozen=True, slots=True)
class EquilibriumSpeciesAmount:
    """Equilibrium amount and fractions for one retained species.

    Attributes
    ----------
    species : str
        Species key.
    moles : float
        Amount in mol.
    mole_fraction, mass_fraction : float
        Fractions in the product mixture.
    """

    species: str
    moles: float
    mole_fraction: float
    mass_fraction: float


@dataclass(frozen=True, slots=True)
class EquilibriumResiduals:
    """Normalized conservation and state-closure residuals.

    Attributes
    ----------
    maximum_element_relative : float
        Largest relative difference between product and reactant element
        amounts.
    mole_fraction_sum : float
        ``|sum(x_j) - 1|``.
    energy_relative : float or None
        Relative enthalpy (HP) or internal-energy (UV) closure residual,
        normalized by ``max(|target|, 1 J)``.
    volume_relative : float or None
        Relative volume closure residual (UV).
    """

    maximum_element_relative: float
    mole_fraction_sum: float
    energy_relative: float | None = None
    volume_relative: float | None = None


@dataclass(frozen=True, slots=True)
class EquilibriumConvergence:
    """Iteration record retained with every successful solution.

    Attributes
    ----------
    converged : bool
        Whether the solution was accepted.
    composition_iterations : int
        Iterations of the composition solve at the final state.
    temperature_iterations, pressure_iterations : int
        Iterations of the temperature and pressure closures.
    maximum_residual : float
        Largest of the composition, element, mole-fraction and closure
        residuals.
    """

    converged: bool
    composition_iterations: int
    temperature_iterations: int
    pressure_iterations: int
    maximum_residual: float


@dataclass(frozen=True, slots=True)
class EquilibriumResult:
    """Provider-neutral, versioned equilibrium result.

    The fields follow the equilibrium-state description of [1]_ and are filled
    by the local Gibbs solver or the NASA CEA adapter.

    Attributes
    ----------
    schema_version, provider, provider_version : str
        Result schema and producing provider.
    problem_type : str
        ``"TP"``, ``"HP"`` or ``"UV"``.
    temperature_k, pressure_pa : float
        Equilibrium state.
    total_moles, total_mass_kg : float
        Product amount and mass.
    molecular_mass_kg_mol, gas_constant_j_kg_k : float
        Mixture molar mass and specific gas constant.
    cp_j_kg_k, cv_j_kg_k, gamma : float
        Specific heats and their ratio as reported by the provider (frozen
        composition for the local provider).
    enthalpy_j, internal_energy_j, entropy_j_k : float
        Extensive state properties of the product mixture.
    species : tuple of EquilibriumSpeciesAmount
        Retained product species in decreasing amount.
    element_moles : tuple of (str, float)
        Reactant element inventory.
    residuals : EquilibriumResiduals
        Conservation and closure residuals.
    convergence : EquilibriumConvergence
        Iteration record.
    assumptions : tuple of str
        Model assumptions of the provider.
    warnings : tuple of WarningMessage
        Provider warnings.
    source, source_url : str
        Source of the model and data.

    References
    ----------
    .. [1] NASA RP-1311, "Computer Program for Calculation of Complex Chemical
           Equilibrium Compositions and Applications, Part I: Analysis".
           https://ntrs.nasa.gov/citations/19950013764
    """

    schema_version: str
    provider: str
    provider_version: str
    problem_type: str
    temperature_k: float
    pressure_pa: float
    total_moles: float
    total_mass_kg: float
    molecular_mass_kg_mol: float
    gas_constant_j_kg_k: float
    cp_j_kg_k: float
    cv_j_kg_k: float
    gamma: float
    enthalpy_j: float
    internal_energy_j: float
    entropy_j_k: float
    species: tuple[EquilibriumSpeciesAmount, ...]
    element_moles: tuple[tuple[str, float], ...]
    residuals: EquilibriumResiduals
    convergence: EquilibriumConvergence
    assumptions: tuple[str, ...]
    warnings: tuple[WarningMessage, ...]
    source: str
    source_url: str


@runtime_checkable
class EquilibriumProvider(Protocol):
    """Common interface implemented by local and external equilibrium engines.

    Attributes
    ----------
    name, version : str
        Provider identification copied into every result.
    """

    name: str
    version: str

    def solve(self, problem: EquilibriumProblem) -> EquilibriumResult:
        """Solve a validated equilibrium problem.

        Parameters
        ----------
        problem : EquilibriumProblem
            Problem to solve.

        Returns
        -------
        EquilibriumResult
            Provider-neutral result.
        """


def _reactants(values: Mapping[str, float]) -> tuple[ReactantAmount, ...]:
    """Validate reactant amounts and return them sorted by upper-case species key."""
    if not values:
        raise InputError("At least one reactant is required.", field="reactants")
    result: list[ReactantAmount] = []
    for key, raw_value in values.items():
        try:
            value = float(raw_value)
        except (TypeError, ValueError) as error:
            raise InputError(
                f"Reactant amount for {key} must be numeric.", field="reactants"
            ) from error
        if not isfinite(value) or value <= 0.0:
            raise DomainError(
                f"Reactant amount for {key} must be finite and positive.",
                field="reactants",
            )
        get_reactant(str(key))
        result.append(ReactantAmount(str(key).upper(), value))
    return tuple(sorted(result, key=lambda item: item.species))


def _solve_linear(matrix: list[list[float]], vector: list[float]) -> list[float]:
    """Solve a small dense system by Gauss-Jordan elimination with partial pivoting.

    Raises
    ------
    ConvergenceError
        If a pivot magnitude falls below 1e-14
        (code ``equilibrium_singular_jacobian``).
    """

    size = len(vector)
    augmented = [row[:] + [value] for row, value in zip(matrix, vector, strict=True)]
    for column in range(size):
        pivot = max(range(column, size), key=lambda row: abs(augmented[row][column]))
        if abs(augmented[pivot][column]) < 1.0e-14:
            raise ConvergenceError(
                "Equilibrium Jacobian became singular.",
                code="equilibrium_singular_jacobian",
                details={"column": column},
            )
        augmented[column], augmented[pivot] = augmented[pivot], augmented[column]
        divisor = augmented[column][column]
        for entry in range(column, size + 1):
            augmented[column][entry] /= divisor
        for row in range(size):
            if row == column:
                continue
            factor = augmented[row][column]
            for entry in range(column, size + 1):
                augmented[row][entry] -= factor * augmented[column][entry]
    return [augmented[row][-1] for row in range(size)]


def _logsumexp(values: list[float]) -> float:
    """Return ``log(sum(exp(values)))`` evaluated without overflow."""
    maximum = max(values)
    return maximum + log(sum(exp(value - maximum) for value in values))


# The ten NASA/TP-2002-211556 combustion species; enough to bound the major
# composition of H/C/N/O systems quickly when Newton needs a warm start.
_MAJOR_SPECIES = frozenset({"H2", "O2", "H2O", "OH", "H", "O", "N2", "CO", "CO2", "CH4"})


def _potentials_from_basis(
    amounts: Mapping[str, float],
    candidates: tuple[Species, ...],
    g_over_rt: list[float],
    elements: tuple[str, ...],
    pressure_term: float,
) -> list[float] | None:
    """Return element potentials whose weights reproduce a basis composition.

    For every basis species ``j``, ``sum_i a_ij pi_i = ln x_j + g_j/(RT) -
    pressure_term``; the minimum-norm solution ``pi = A^T (A A^T)^-1 b`` is
    used when the basis has fewer species than elements. Returns ``None``
    when the basis is empty or singular.
    """
    total = sum(value for value in amounts.values() if value > 0.0)
    index = {species.key: position for position, species in enumerate(candidates)}
    rows: list[list[float]] = []
    vector: list[float] = []
    for key, amount in amounts.items():
        if amount <= 0.0 or key not in index:
            continue
        formula = dict(candidates[index[key]].formula)
        rows.append([formula.get(element, 0.0) for element in elements])
        vector.append(log(amount / total) + g_over_rt[index[key]] - pressure_term)
    if not rows:
        return None
    gram = [[sum(a * b for a, b in zip(left, right, strict=True)) for right in rows] for left in rows]
    try:
        weights = _solve_linear(gram, vector)
    except ConvergenceError:
        return None
    return [
        sum(weight * row[column] for weight, row in zip(weights, rows, strict=True))
        for column in range(len(elements))
    ]


class LocalEquilibriumProvider:
    """Experimental ideal-gas equilibrium solver for the bundled species set.

    Attributes
    ----------
    name, version : str
        Provider identification.
    maximum_iterations : int
        Newton iteration limit of the composition solve.
    residual_tolerance : float
        Acceptance tolerance of the composition residuals.

    Notes
    -----
    Assumptions: ideal-gas mixture, gas-phase species only, bundled NASA Glenn
    thermodynamic data [2]_, temperatures between 200 and 6000 K. Results carry
    an ``experimental_equilibrium_provider`` warning; rocket design results
    should be checked against NASA CEA [1]_.

    References
    ----------
    .. [1] NASA RP-1311, "Computer Program for Calculation of Complex Chemical
           Equilibrium Compositions and Applications, Part I: Analysis".
           https://ntrs.nasa.gov/citations/19950013764
    .. [2] NASA/TP-2002-211556. NASA Glenn, "NASA Glenn Coefficients for
           Calculating Thermodynamic Properties of Individual Species".
           https://ntrs.nasa.gov/citations/20020085330
    """

    name = "local-gibbs"
    version = "0.2-experimental"
    maximum_iterations = 100
    residual_tolerance = 2.0e-9

    def __init__(self) -> None:
        # Last converged element potentials per element tuple (warm start only;
        # every returned state still satisfies the residual tolerance).
        self._last_potentials: dict[tuple[str, ...], list[float]] = {}

    def solve(self, problem: EquilibriumProblem) -> EquilibriumResult:
        """Solve TP directly and HP/UV through nested state closure.

        Parameters
        ----------
        problem : EquilibriumProblem
            Problem with ``schema_version="1.0"``.

        Returns
        -------
        EquilibriumResult
            Equilibrium composition, properties, residuals and warnings.

        Raises
        ------
        InputError
            If the schema or problem type is unsupported or a required state
            variable is missing.
        DomainError
            If the state is outside the model range or the energy/volume
            target is not bracketed.
        ConvergenceError
            If the composition or closure iterations fail.

        Notes
        -----
        TP problems minimize Gibbs energy at the assigned state. HP problems
        bisect temperature in ``[1800, 6000]`` K on the enthalpy residual. UV
        problems bisect temperature on the internal-energy residual, with an
        inner bisection on ``ln p`` for the volume at each temperature [1]_.

        References
        ----------
        .. [1] NASA RP-1311, "Computer Program for Calculation of Complex Chemical
               Equilibrium Compositions and Applications, Part I: Analysis".
               https://ntrs.nasa.gov/citations/19950013764
        """

        kind = problem.problem_type.strip().upper()
        if problem.schema_version != "1.0":
            raise InputError(
                "Unsupported equilibrium problem schema.",
                code="unsupported_schema",
                field="schema_version",
            )
        if kind == "TP":
            return self._solve_tp_problem(problem)
        if kind == "HP":
            return self._solve_hp_problem(problem)
        if kind == "UV":
            return self._solve_uv_problem(problem)
        raise InputError(
            "Equilibrium problem type must be TP, HP, or UV.",
            field="problem_type",
        )

    def _inventory(self, problem: EquilibriumProblem) -> dict[str, float]:
        """Return reactant element amounts in mol, sorted by element symbol."""
        inventory: dict[str, float] = {}
        for reactant in problem.reactants:
            for element, coefficient in get_reactant(reactant.species).formula:
                inventory[element] = inventory.get(element, 0.0) + reactant.moles * coefficient
        return dict(sorted(inventory.items()))

    def _candidates(
        self, problem: EquilibriumProblem, inventory: Mapping[str, float]
    ) -> tuple[Species, ...]:
        """Return candidate product species composed only of reactant elements.

        Raises
        ------
        DomainError
            If no candidate matches (code ``no_equilibrium_species``) or the
            candidates do not contain every reactant element
            (code ``incomplete_element_basis``).
        """
        active_elements = set(inventory)
        requested = (
            tuple(get_species(key) for key in problem.candidate_species)
            if problem.candidate_species is not None
            else tuple(species_catalog().values())
        )
        candidates = tuple(
            species
            for species in requested
            if set(dict(species.formula)).issubset(active_elements)
        )
        if not candidates:
            raise DomainError(
                "No bundled product species match the reactant elements.",
                code="no_equilibrium_species",
            )
        represented = {
            element
            for species in candidates
            for element, coefficient in species.formula
            if coefficient > 0.0
        }
        missing = sorted(active_elements - represented)
        if missing:
            raise DomainError(
                "Candidate species do not represent every reactant element.",
                code="incomplete_element_basis",
                details={"missing_elements": missing},
            )
        return candidates

    def _composition_at_tp(
        self,
        *,
        inventory: Mapping[str, float],
        candidates: tuple[Species, ...],
        temperature_k: float,
        pressure_pa: float,
    ) -> tuple[dict[str, float], int, float]:
        """Return the equilibrium species amounts at fixed temperature and pressure.

        Returns
        -------
        tuple
            ``(amounts, iterations, residual)``: moles per species key, Newton
            iterations and the final residual norm.

        Raises
        ------
        DomainError
            If the temperature is outside 200--6000 K or the pressure is not
            positive.
        ConvergenceError
            If the Newton solve fails and no feasible boundary composition
            exists.

        Notes
        -----
        The element-potential form of Gibbs minimization is used [1]_: with
        element potentials ``pi_i``, mole fractions are
        ``x_j proportional to exp(-g_j/(R T) + sum_i a_ij pi_i - ln(p/p_ref))``
        with standard-state ``g_j = h_j - T s_j`` from NASA Glenn data [2]_.
        The unknowns are found by damped Newton iteration (forward-difference
        Jacobian, step 1e-5, halving line search) so that the fractions sum to
        one and every element ratio matches the reactant inventory. If the
        Jacobian is singular or no step reduces the residual, the lowest-Gibbs
        feasible boundary composition is returned instead. Candidates whose
        polynomial range does not cover ``temperature_k`` are left out.

        References
        ----------
        .. [1] NASA RP-1311, "Computer Program for Calculation of Complex Chemical
               Equilibrium Compositions and Applications, Part I: Analysis".
               https://ntrs.nasa.gov/citations/19950013764
        .. [2] NASA/TP-2002-211556. NASA Glenn, "NASA Glenn Coefficients for
               Calculating Thermodynamic Properties of Individual Species".
               https://ntrs.nasa.gov/citations/20020085330
        """
        if not isfinite(temperature_k) or not 200.0 <= temperature_k <= 6000.0:
            raise DomainError(
                "Local equilibrium temperature must be between 200 and 6000 K.",
                field="temperature_k",
            )
        if not isfinite(pressure_pa) or pressure_pa <= 0.0:
            raise DomainError("Equilibrium pressure must be positive.", field="pressure_pa")
        if pressure_pa > MAXIMUM_IDEAL_GAS_PRESSURE_PA:
            raise DomainError(
                f"Equilibrium pressure {pressure_pa:.4g} Pa exceeds the {MAXIMUM_IDEAL_GAS_PRESSURE_PA:.0e} Pa "
                "limit of the ideal-gas mixture model.",
                code="pressure_outside_ideal_gas_range", field="pressure_pa")
        # Some NASA records start at 300 K; below that they are simply absent.
        candidates = tuple(
            species
            for species in candidates
            if species.minimum_temperature_k <= temperature_k <= species.maximum_temperature_k
        )

        elements = tuple(inventory)
        reference_element = elements[0]
        target_ratios = {
            element: inventory[element] / inventory[reference_element]
            for element in elements[1:]
        }
        atom_matrix = [dict(species.formula) for species in candidates]
        atom_rows = [tuple(record.get(element, 0.0) for element in elements)
                     for record in atom_matrix]
        atom_columns = {element: tuple(row[column] for row in atom_rows)
                        for column, element in enumerate(elements)}
        g_over_rt = []
        for species in candidates:
            properties = species.properties(temperature_k)
            gibbs = properties.enthalpy_molar_j_mol - (
                temperature_k * properties.entropy_molar_j_mol_k
            )
            g_over_rt.append(
                gibbs / (UNIVERSAL_GAS_CONSTANT_J_MOL_K * temperature_k)
            )
        pressure_term = log(REFERENCE_PRESSURE_PA / pressure_pa)

        def residuals(potentials: list[float]) -> tuple[list[float], list[float]]:
            log_weights = [
                pressure_term
                - gibbs
                + sum(map(mul, atom_rows[index], potentials))
                for index, gibbs in enumerate(g_over_rt)
            ]
            normalization = _logsumexp(log_weights)
            fractions = [exp(value - normalization) for value in log_weights]
            atoms_per_mole = {
                element: sum(map(mul, fractions, counts))
                for element, counts in atom_columns.items()
            }
            values = [normalization]
            reference_atoms = atoms_per_mole[reference_element]
            for element in elements[1:]:
                ratio = (
                    atoms_per_mole[element] / reference_atoms
                    if reference_atoms > 0.0
                    else 0.0
                )
                if ratio <= 0.0 or not isfinite(ratio):
                    values.append(1.0e6)
                else:
                    values.append(
                        log(ratio) - log(target_ratios[element])
                    )
            return values, fractions

        solved_potentials: list[list[float]] = [[]]

        def newton(start: list[float]) -> tuple[list[float], int, float] | None:
            potentials = start[:]
            for iteration in range(1, self.maximum_iterations + 1):
                values, fractions = residuals(potentials)
                norm = max(abs(value) for value in values)
                if norm <= self.residual_tolerance:
                    solved_potentials[0] = potentials
                    return fractions, iteration, norm
                step_size = 1.0e-5
                # One residual vector supplies every row of a difference column.
                # Re-evaluating that identical vector per row adds no information.
                shifted_columns: list[list[float]] = []
                for column in range(len(elements)):
                    shifted = potentials[:]
                    shifted[column] += step_size
                    shifted_values, _ = residuals(shifted)
                    shifted_columns.append(shifted_values)
                jacobian = [[(shifted_columns[column][row] - values[row]) / step_size
                             for column in range(len(elements))]
                            for row in range(len(elements))]
                try:
                    delta = _solve_linear(jacobian, [-value for value in values])
                except ConvergenceError:
                    return None
                damping = 1.0
                while damping >= 1.0e-6:
                    trial = [
                        value + damping * change
                        for value, change in zip(potentials, delta, strict=True)
                    ]
                    trial_values, _ = residuals(trial)
                    if max(abs(value) for value in trial_values) < norm:
                        potentials = trial
                        break
                    damping *= 0.5
                else:
                    return None
            return None

        # Continuation: successive states of one expansion are close, so the
        # last converged potentials for this element set are tried first.
        previous = self._last_potentials.get(elements)
        solved = newton(previous) if previous is not None else None
        if solved is None:
            solved = newton([0.0] * len(elements))
        if solved is None:
            # Warm start: element potentials that reproduce the zero-trace
            # boundary composition of the major species exactly.
            major = tuple(species for species in candidates if species.key in _MAJOR_SPECIES)
            try:
                boundary, _, _ = self._boundary_composition(
                    inventory=inventory,
                    candidates=major or candidates,
                    temperature_k=temperature_k,
                    pressure_pa=pressure_pa,
                    iterations=0,
                )
            except ConvergenceError:
                boundary = {}
            start = _potentials_from_basis(
                boundary, candidates, g_over_rt, elements, pressure_term
            )
            if start is not None:
                solved = newton(start)
        if solved is None:
            return self._boundary_composition(
                inventory=inventory,
                candidates=candidates,
                temperature_k=temperature_k,
                pressure_pa=pressure_pa,
                iterations=self.maximum_iterations,
            )
        fractions, iteration, norm = solved
        self._last_potentials[elements] = solved_potentials[0]

        atoms_per_mole = {
            element: sum(map(mul, fractions, counts))
            for element, counts in atom_columns.items()
        }
        total_moles = inventory[reference_element] / atoms_per_mole[reference_element]
        amounts = {
            species.key: total_moles * fraction
            for species, fraction in zip(candidates, fractions, strict=True)
        }
        return amounts, iteration, norm

    def _boundary_composition(
        self,
        *,
        inventory: Mapping[str, float],
        candidates: tuple[Species, ...],
        temperature_k: float,
        pressure_pa: float,
        iterations: int,
    ) -> tuple[dict[str, float], int, float]:
        """Select the lowest-G feasible boundary basis when traces underflow.

        At low combustion temperatures, dissociated fractions can be far below
        double-precision sensitivity and the interior atomic-potential Newton
        system becomes singular. Enumerating species subsets of at most the
        number of elements, solving the scaled element balance by least
        squares, and keeping non-negative solutions that conserve every
        element within 1e-8 retains conservation while treating the omitted
        trace species as zero. The candidate with the lowest ideal-mixture
        Gibbs energy [1]_ is returned.

        Raises
        ------
        ConvergenceError
            If no feasible basis exists (code ``equilibrium_boundary_failure``).

        References
        ----------
        .. [1] NASA RP-1311, "Computer Program for Calculation of Complex Chemical
               Equilibrium Compositions and Applications, Part I: Analysis".
               https://ntrs.nasa.gov/citations/19950013764
        """

        elements = tuple(inventory)
        maximum_basis_size = min(len(elements), len(candidates))
        best: tuple[float, dict[str, float], float] | None = None
        for basis_size in range(1, maximum_basis_size + 1):
            for basis in combinations(candidates, basis_size):
                scaled = [
                    [dict(species.formula).get(element, 0.0) / inventory[element]
                     for species in basis]
                    for element in elements
                ]
                normal = [
                    [
                        sum(row[left] * row[right] for row in scaled)
                        for right in range(basis_size)
                    ]
                    for left in range(basis_size)
                ]
                vector = [
                    sum(row[column] for row in scaled)
                    for column in range(basis_size)
                ]
                try:
                    solution = _solve_linear(normal, vector)
                except ConvergenceError:
                    continue
                if any(value < -1.0e-12 for value in solution):
                    continue
                solution = [max(0.0, value) for value in solution]
                calculated = {
                    element: sum(
                        amount * dict(species.formula).get(element, 0.0)
                        for species, amount in zip(basis, solution, strict=True)
                    )
                    for element in elements
                }
                residual = max(
                    abs(calculated[element] - target) / target
                    for element, target in inventory.items()
                )
                if residual > 1.0e-8:
                    continue
                total_moles = sum(solution)
                if total_moles <= 0.0:
                    continue
                gibbs = 0.0
                for species, amount in zip(basis, solution, strict=True):
                    if amount <= 0.0:
                        continue
                    properties = species.properties(temperature_k)
                    standard_gibbs = properties.enthalpy_molar_j_mol - (
                        temperature_k * properties.entropy_molar_j_mol_k
                    )
                    chemical_potential = standard_gibbs + (
                        UNIVERSAL_GAS_CONSTANT_J_MOL_K
                        * temperature_k
                        * log((amount / total_moles) * pressure_pa / REFERENCE_PRESSURE_PA)
                    )
                    gibbs += amount * chemical_potential
                amounts = {
                    species.key: amount
                    for species, amount in zip(basis, solution, strict=True)
                    if amount > 0.0
                }
                if best is None or gibbs < best[0]:
                    best = (gibbs, amounts, residual)
        if best is None:
            raise ConvergenceError(
                "Local equilibrium could not construct a feasible boundary composition.",
                code="equilibrium_boundary_failure",
                details={"iterations": iterations},
            )
        _, amounts, residual = best
        return amounts, iterations, residual

    def _build_result(
        self,
        *,
        problem: EquilibriumProblem,
        temperature_k: float,
        pressure_pa: float,
        amounts: Mapping[str, float],
        composition_iterations: int,
        composition_residual: float,
        temperature_iterations: int = 0,
        pressure_iterations: int = 0,
        energy_target_j: float | None = None,
        volume_target_m3: float | None = None,
    ) -> EquilibriumResult:
        """Assemble an :class:`EquilibriumResult` from solved species amounts.

        Species below 1e-14 mol are dropped. Mixture properties use frozen
        composition at the final state; ``U = H - n R T`` for the ideal-gas
        mixture. Element, mole-fraction, energy and volume residuals are
        recomputed from the retained amounts.
        """
        retained = {key: value for key, value in amounts.items() if value > 1.0e-14}
        total_moles = sum(retained.values())
        mole_fractions = {key: value / total_moles for key, value in retained.items()}
        mixture = Mixture.from_fractions(mole_fractions)
        properties = mixture.properties(temperature_k, pressure_pa=pressure_pa)
        total_mass = sum(
            amount * get_species(key).molecular_mass_kg_mol
            for key, amount in retained.items()
        )
        enthalpy = sum(
            amount * get_species(key).properties(temperature_k).enthalpy_molar_j_mol
            for key, amount in retained.items()
        )
        internal_energy = enthalpy - (
            total_moles * UNIVERSAL_GAS_CONSTANT_J_MOL_K * temperature_k
        )
        entropy = properties.entropy_j_kg_k * total_mass
        mass_fractions = {
            component.species.key: component.mass_fraction
            for component in mixture.components
        }
        inventory = self._inventory(problem)
        calculated_inventory = {
            element: sum(
                amount * dict(get_species(key).formula).get(element, 0.0)
                for key, amount in retained.items()
            )
            for element in inventory
        }
        element_residual = max(
            abs(calculated_inventory[element] - target) / target
            for element, target in inventory.items()
        )
        energy_residual = (
            None
            if energy_target_j is None
            else abs(internal_energy if problem.problem_type == "UV" else enthalpy - energy_target_j)
            / max(abs(energy_target_j), 1.0)
        )
        if energy_target_j is not None and problem.problem_type == "UV":
            energy_residual = abs(internal_energy - energy_target_j) / max(
                abs(energy_target_j), 1.0
            )
        actual_volume = (
            total_moles * UNIVERSAL_GAS_CONSTANT_J_MOL_K * temperature_k / pressure_pa
        )
        volume_residual = (
            None
            if volume_target_m3 is None
            else abs(actual_volume - volume_target_m3) / volume_target_m3
        )
        residuals = EquilibriumResiduals(
            maximum_element_relative=element_residual,
            mole_fraction_sum=abs(sum(mole_fractions.values()) - 1.0),
            energy_relative=energy_residual,
            volume_relative=volume_residual,
        )
        maximum_residual = max(
            composition_residual,
            element_residual,
            residuals.mole_fraction_sum,
            energy_residual or 0.0,
            volume_residual or 0.0,
        )
        species_amounts = tuple(
            EquilibriumSpeciesAmount(
                species=key,
                moles=amount,
                mole_fraction=mole_fractions[key],
                mass_fraction=mass_fractions[key],
            )
            for key, amount in sorted(retained.items(), key=lambda item: -item[1])
        )
        return EquilibriumResult(
            schema_version="1.0",
            provider=self.name,
            provider_version=self.version,
            problem_type=problem.problem_type,
            temperature_k=temperature_k,
            pressure_pa=pressure_pa,
            total_moles=total_moles,
            total_mass_kg=total_mass,
            molecular_mass_kg_mol=properties.molecular_mass_kg_mol,
            gas_constant_j_kg_k=properties.gas_constant_j_kg_k,
            cp_j_kg_k=properties.cp_j_kg_k,
            cv_j_kg_k=properties.cv_j_kg_k,
            gamma=properties.gamma,
            enthalpy_j=enthalpy,
            internal_energy_j=internal_energy,
            entropy_j_k=entropy,
            species=species_amounts,
            element_moles=tuple(inventory.items()),
            residuals=residuals,
            convergence=EquilibriumConvergence(
                converged=True,
                composition_iterations=composition_iterations,
                temperature_iterations=temperature_iterations,
                pressure_iterations=pressure_iterations,
                maximum_residual=maximum_residual,
            ),
            assumptions=(
                "ideal-gas mixture",
                "gas-phase species only",
                "bundled NASA Glenn thermodynamic data",
                "experimental local Gibbs minimizer",
            ),
            warnings=(
                WarningMessage(
                    code="experimental_equilibrium_provider",
                    message=(
                        "The local equilibrium solver is experimental; validate rocket "
                        "design results against NASA CEA."
                    ),
                ),
            ),
            source="Gordon and McBride, NASA RP-1311; NASA/TP-2002-211556",
            source_url="https://ntrs.nasa.gov/citations/19950013764",
        )

    def _solve_tp_problem(self, problem: EquilibriumProblem) -> EquilibriumResult:
        """Solve an assigned temperature and pressure problem."""
        if problem.temperature_k is None or problem.pressure_pa is None:
            raise InputError("TP requires temperature and pressure.", field="problem")
        inventory = self._inventory(problem)
        candidates = self._candidates(problem, inventory)
        amounts, iterations, residual = self._composition_at_tp(
            inventory=inventory,
            candidates=candidates,
            temperature_k=problem.temperature_k,
            pressure_pa=problem.pressure_pa,
        )
        return self._build_result(
            problem=problem,
            temperature_k=problem.temperature_k,
            pressure_pa=problem.pressure_pa,
            amounts=amounts,
            composition_iterations=iterations,
            composition_residual=residual,
        )

    def _reactant_energy(self, problem: EquilibriumProblem) -> tuple[float, float]:
        """Return reactant ``(H, U)`` in J.

        Gas reactants are evaluated at ``reactant_temperature_k``; condensed
        reactants contribute their assigned enthalpy at their own storage
        temperature, as in NASA CEA [1]_, and no ``pV`` term.

        Raises
        ------
        InputError
            If a gas reactant is present and ``reactant_temperature_k`` is
            missing.

        References
        ----------
        .. [1] NASA RP-1311, "Computer Program for Calculation of Complex Chemical
               Equilibrium Compositions and Applications, Part II: Users Manual
               and Program Description". https://ntrs.nasa.gov/citations/19960044559
        """
        temperature = problem.reactant_temperature_k
        enthalpy = 0.0
        gas_pv = 0.0
        for reactant in problem.reactants:
            record = get_reactant(reactant.species)
            if isinstance(record, CondensedReactant):
                enthalpy += reactant.moles * record.enthalpy_j_mol
                continue
            if temperature is None:
                raise InputError(
                    "Combustion problems with gas reactants require reactant_temperature_k.",
                    field="reactant_temperature_k",
                )
            enthalpy += reactant.moles * record.properties(temperature).enthalpy_molar_j_mol
            gas_pv += reactant.moles * UNIVERSAL_GAS_CONSTANT_J_MOL_K * temperature
        return enthalpy, enthalpy - gas_pv

    def _bisect_state(
        self,
        function,
        *,
        lower: float,
        upper: float,
        tolerance: float = 1.0e-8,
        maximum_iterations: int = 80,
    ) -> tuple[float, int]:
        """Return a bracketed root of ``function`` by bisection and the iteration count.

        Raises
        ------
        DomainError
            If the end values do not bracket a root
            (code ``equilibrium_state_not_bracketed``).
        ConvergenceError
            If the iteration limit is reached.
        """
        low_value = function(lower)
        high_value = function(upper)
        if low_value == 0.0:
            return lower, 0
        if high_value == 0.0:
            return upper, 0
        if low_value * high_value > 0.0:
            raise DomainError(
                "Equilibrium state target is not bracketed by model limits.",
                code="equilibrium_state_not_bracketed",
                details={"lower_residual": low_value, "upper_residual": high_value},
            )
        midpoint = lower
        for iteration in range(1, maximum_iterations + 1):
            midpoint = 0.5 * (lower + upper)
            value = function(midpoint)
            if abs(value) <= tolerance or (upper - lower) <= tolerance:
                return midpoint, iteration
            if low_value * value <= 0.0:
                upper = midpoint
            else:
                lower = midpoint
                low_value = value
        raise ConvergenceError(
            "Equilibrium state closure exceeded its iteration limit.",
            details={"iterations": maximum_iterations},
        )

    def _solve_hp_problem(self, problem: EquilibriumProblem) -> EquilibriumResult:
        """Solve an assigned enthalpy and pressure (adiabatic combustion) problem."""
        if problem.pressure_pa is None:
            raise InputError("HP requires pressure.", field="pressure_pa")
        pressure_pa = problem.pressure_pa
        reactant_enthalpy, _ = self._reactant_energy(problem)
        target = (
            problem.target_enthalpy_j
            if problem.target_enthalpy_j is not None
            else reactant_enthalpy
        )
        inventory = self._inventory(problem)
        candidates = self._candidates(problem, inventory)
        last: tuple[dict[str, float], int, float] | None = None

        def energy_residual(temperature: float) -> float:
            nonlocal last
            last = self._composition_at_tp(
                inventory=inventory,
                candidates=candidates,
                temperature_k=temperature,
                pressure_pa=pressure_pa,
            )
            amounts, _, _ = last
            enthalpy = sum(
                amount * get_species(key).properties(temperature).enthalpy_molar_j_mol
                for key, amount in amounts.items()
            )
            return (enthalpy - target) / max(abs(target), 1.0)

        temperature, temperature_iterations = self._bisect_state(
            # Below roughly 1800 K the dissociated fractions in the bundled
            # combustion systems can fall below double-precision sensitivity.
            # HP combustion roots are far above that numerical floor.
            energy_residual,
            lower=1800.0,
            upper=6000.0,
        )
        energy_residual(temperature)
        assert last is not None
        amounts, composition_iterations, composition_residual = last
        return self._build_result(
            problem=problem,
            temperature_k=temperature,
            pressure_pa=problem.pressure_pa,
            amounts=amounts,
            composition_iterations=composition_iterations,
            composition_residual=composition_residual,
            temperature_iterations=temperature_iterations,
            energy_target_j=target,
        )

    def _solve_uv_problem(self, problem: EquilibriumProblem) -> EquilibriumResult:
        """Solve an assigned internal energy and volume (constant-volume) problem."""
        if problem.volume_m3 is None or problem.volume_m3 <= 0.0:
            raise InputError("UV requires a positive volume.", field="volume_m3")
        volume_m3 = problem.volume_m3
        _, reactant_internal_energy = self._reactant_energy(problem)
        target = (
            problem.target_internal_energy_j
            if problem.target_internal_energy_j is not None
            else reactant_internal_energy
        )
        inventory = self._inventory(problem)
        candidates = self._candidates(problem, inventory)
        pressure_iteration_total = 0
        last: tuple[float, dict[str, float], int, float] | None = None

        def pressure_for_temperature(temperature: float) -> tuple[float, dict[str, float], int, float]:
            nonlocal pressure_iteration_total
            latest: tuple[dict[str, float], int, float] | None = None

            def volume_residual(log_pressure: float) -> float:
                nonlocal latest
                pressure = exp(log_pressure)
                latest = self._composition_at_tp(
                    inventory=inventory,
                    candidates=candidates,
                    temperature_k=temperature,
                    pressure_pa=pressure,
                )
                amounts, _, _ = latest
                actual_volume = (
                    sum(amounts.values())
                    * UNIVERSAL_GAS_CONSTANT_J_MOL_K
                    * temperature
                    / pressure
                )
                return log(actual_volume / volume_m3)

            log_pressure, iterations = self._bisect_state(
                volume_residual,
                lower=log(1.0),
                upper=log(MAXIMUM_IDEAL_GAS_PRESSURE_PA),
                tolerance=1.0e-10,
            )
            pressure_iteration_total += iterations
            volume_residual(log_pressure)
            assert latest is not None
            amounts, composition_iterations, composition_residual = latest
            return exp(log_pressure), amounts, composition_iterations, composition_residual

        def internal_energy_residual(temperature: float) -> float:
            nonlocal last
            last = pressure_for_temperature(temperature)
            _, amounts, _, _ = last
            enthalpy = sum(
                amount * get_species(key).properties(temperature).enthalpy_molar_j_mol
                for key, amount in amounts.items()
            )
            internal_energy = enthalpy - (
                sum(amounts.values()) * UNIVERSAL_GAS_CONSTANT_J_MOL_K * temperature
            )
            return (internal_energy - target) / max(abs(target), 1.0)

        temperature, temperature_iterations = self._bisect_state(
            internal_energy_residual, lower=1800.0, upper=6000.0
        )
        internal_energy_residual(temperature)
        assert last is not None
        pressure, amounts, composition_iterations, composition_residual = last
        result = self._build_result(
            problem=problem,
            temperature_k=temperature,
            pressure_pa=pressure,
            amounts=amounts,
            composition_iterations=composition_iterations,
            composition_residual=composition_residual,
            temperature_iterations=temperature_iterations,
            pressure_iterations=pressure_iteration_total,
            energy_target_j=target,
            volume_target_m3=problem.volume_m3,
        )
        return replace(result, problem_type="UV")

