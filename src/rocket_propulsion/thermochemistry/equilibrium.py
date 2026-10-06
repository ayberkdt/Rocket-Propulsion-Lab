"""Chemical-equilibrium provider contract and an experimental local solver.

The local provider minimizes ideal-gas Gibbs energy with elemental constraints.
It is intentionally limited to the bundled gas species and 200--6000 K.  The
provider boundary also permits the validated NASA CEA adapter to be selected
without coupling the rest of the application to an external executable.
"""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass, replace
from itertools import combinations
from math import exp, isfinite, log
from typing import Protocol, runtime_checkable

from rocket_propulsion.core.errors import ConvergenceError, DomainError, InputError
from rocket_propulsion.core.metadata import WarningMessage

from .mixtures import REFERENCE_PRESSURE_PA, Mixture
from .nasa_polynomials import UNIVERSAL_GAS_CONSTANT_J_MOL_K
from .species import Species, get_species, species_catalog


@dataclass(frozen=True, slots=True)
class ReactantAmount:
    """A positive molar amount of one bundled species."""

    species: str
    moles: float


@dataclass(frozen=True, slots=True)
class EquilibriumProblem:
    """A versioned TP, HP, or UV ideal-gas equilibrium problem."""

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
        """Build an assigned-temperature-and-pressure problem."""

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
        """Build an assigned-enthalpy-and-pressure combustion problem."""

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
        """Build an assigned-internal-energy-and-volume problem."""

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
    """Equilibrium amount and fractions for one retained species."""

    species: str
    moles: float
    mole_fraction: float
    mass_fraction: float


@dataclass(frozen=True, slots=True)
class EquilibriumResiduals:
    """Normalized conservation and state-closure residuals."""

    maximum_element_relative: float
    mole_fraction_sum: float
    energy_relative: float | None = None
    volume_relative: float | None = None


@dataclass(frozen=True, slots=True)
class EquilibriumConvergence:
    """Iteration evidence retained with every successful solution."""

    converged: bool
    composition_iterations: int
    temperature_iterations: int
    pressure_iterations: int
    maximum_residual: float


@dataclass(frozen=True, slots=True)
class EquilibriumResult:
    """Provider-neutral, versioned equilibrium result."""

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
    """Common interface implemented by local and external equilibrium engines."""

    name: str
    version: str

    def solve(self, problem: EquilibriumProblem) -> EquilibriumResult:
        """Solve a validated equilibrium problem."""


def _reactants(values: Mapping[str, float]) -> tuple[ReactantAmount, ...]:
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
        get_species(str(key))
        result.append(ReactantAmount(str(key).upper(), value))
    return tuple(sorted(result, key=lambda item: item.species))


def _solve_linear(matrix: list[list[float]], vector: list[float]) -> list[float]:
    """Solve a small dense system with scaled partial pivoting."""

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
    maximum = max(values)
    return maximum + log(sum(exp(value - maximum) for value in values))


class LocalEquilibriumProvider:
    """Experimental ideal-gas equilibrium solver for the bundled species set."""

    name = "local-gibbs"
    version = "0.1-experimental"
    maximum_iterations = 100
    residual_tolerance = 2.0e-9

    def solve(self, problem: EquilibriumProblem) -> EquilibriumResult:
        """Solve TP directly and HP/UV through nested state closure."""

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
        inventory: dict[str, float] = {}
        for reactant in problem.reactants:
            for element, coefficient in get_species(reactant.species).formula:
                inventory[element] = inventory.get(element, 0.0) + reactant.moles * coefficient
        return dict(sorted(inventory.items()))

    def _candidates(
        self, problem: EquilibriumProblem, inventory: Mapping[str, float]
    ) -> tuple[Species, ...]:
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
        if not isfinite(temperature_k) or not 200.0 <= temperature_k <= 6000.0:
            raise DomainError(
                "Local equilibrium temperature must be between 200 and 6000 K.",
                field="temperature_k",
            )
        if not isfinite(pressure_pa) or pressure_pa <= 0.0:
            raise DomainError("Equilibrium pressure must be positive.", field="pressure_pa")

        elements = tuple(inventory)
        reference_element = elements[0]
        target_ratios = {
            element: inventory[element] / inventory[reference_element]
            for element in elements[1:]
        }
        atom_matrix = [dict(species.formula) for species in candidates]
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
                + sum(
                    atom_matrix[index].get(element, 0.0) * potentials[element_index]
                    for element_index, element in enumerate(elements)
                )
                for index, gibbs in enumerate(g_over_rt)
            ]
            normalization = _logsumexp(log_weights)
            fractions = [exp(value - normalization) for value in log_weights]
            atoms_per_mole = {
                element: sum(
                    fraction * atom_matrix[index].get(element, 0.0)
                    for index, fraction in enumerate(fractions)
                )
                for element in elements
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

        potentials = [0.0] * len(elements)
        for iteration in range(1, self.maximum_iterations + 1):
            values, fractions = residuals(potentials)
            norm = max(abs(value) for value in values)
            if norm <= self.residual_tolerance:
                break
            step_size = 1.0e-5
            jacobian: list[list[float]] = []
            for row in range(len(elements)):
                jacobian.append([])
                for column in range(len(elements)):
                    shifted = potentials[:]
                    shifted[column] += step_size
                    shifted_values, _ = residuals(shifted)
                    jacobian[row].append((shifted_values[row] - values[row]) / step_size)
            try:
                delta = _solve_linear(jacobian, [-value for value in values])
            except ConvergenceError:
                return self._boundary_composition(
                    inventory=inventory,
                    candidates=candidates,
                    temperature_k=temperature_k,
                    pressure_pa=pressure_pa,
                    iterations=iteration,
                )
            damping = 1.0
            accepted = False
            while damping >= 1.0e-6:
                trial = [
                    value + damping * change
                    for value, change in zip(potentials, delta, strict=True)
                ]
                trial_values, _ = residuals(trial)
                trial_norm = max(abs(value) for value in trial_values)
                if trial_norm < norm:
                    potentials = trial
                    accepted = True
                    break
                damping *= 0.5
            if not accepted:
                return self._boundary_composition(
                    inventory=inventory,
                    candidates=candidates,
                    temperature_k=temperature_k,
                    pressure_pa=pressure_pa,
                    iterations=iteration,
                )
        else:
            return self._boundary_composition(
                inventory=inventory,
                candidates=candidates,
                temperature_k=temperature_k,
                pressure_pa=pressure_pa,
                iterations=self.maximum_iterations,
            )

        atoms_per_mole = {
            element: sum(
                fraction * atom_matrix[index].get(element, 0.0)
                for index, fraction in enumerate(fractions)
            )
            for element in elements
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
        system becomes singular.  Enumerating element-sized feasible bases
        retains exact conservation while explicitly treating those trace
        species as zero.
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
        temperature = problem.reactant_temperature_k
        if temperature is None:
            raise InputError(
                "Combustion problems require reactant_temperature_k.",
                field="reactant_temperature_k",
            )
        enthalpy = 0.0
        total_moles = 0.0
        for reactant in problem.reactants:
            enthalpy += (
                reactant.moles
                * get_species(reactant.species).properties(temperature).enthalpy_molar_j_mol
            )
            total_moles += reactant.moles
        internal_energy = enthalpy - (
            total_moles * UNIVERSAL_GAS_CONSTANT_J_MOL_K * temperature
        )
        return enthalpy, internal_energy

    def _bisect_state(
        self,
        function,
        *,
        lower: float,
        upper: float,
        tolerance: float = 1.0e-8,
        maximum_iterations: int = 80,
    ) -> tuple[float, int]:
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
        if problem.pressure_pa is None:
            raise InputError("HP requires pressure.", field="pressure_pa")
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
                pressure_pa=problem.pressure_pa,
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
        if problem.volume_m3 is None or problem.volume_m3 <= 0.0:
            raise InputError("UV requires a positive volume.", field="volume_m3")
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
                return log(actual_volume / problem.volume_m3)

            log_pressure, iterations = self._bisect_state(
                volume_residual,
                lower=log(1.0),
                upper=log(1.0e10),
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

