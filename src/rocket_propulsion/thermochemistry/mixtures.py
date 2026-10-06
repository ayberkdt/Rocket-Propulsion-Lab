"""Ideal-gas mixture properties with explicit composition bases."""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass
from math import isclose, log

from rocket_propulsion.core.errors import DomainError, InputError
from rocket_propulsion.core.numerics import solve_bisect_monotonic

from .nasa_polynomials import UNIVERSAL_GAS_CONSTANT_J_MOL_K
from .species import Species, get_species

REFERENCE_PRESSURE_PA = 100_000.0


@dataclass(frozen=True, slots=True)
class MixtureComponent:
    """A resolved species with consistent mole and mass fractions."""

    species: Species
    mole_fraction: float
    mass_fraction: float


@dataclass(frozen=True, slots=True)
class MixtureProperties:
    """Temperature-dependent properties of an ideal gas mixture."""

    temperature_k: float
    pressure_pa: float
    molecular_mass_kg_mol: float
    gas_constant_j_kg_k: float
    cp_j_kg_k: float
    cv_j_kg_k: float
    gamma: float
    enthalpy_j_kg: float
    entropy_j_kg_k: float


@dataclass(frozen=True, slots=True)
class Mixture:
    """An immutable ideal-gas mixture defined by resolved components."""

    components: tuple[MixtureComponent, ...]

    @classmethod
    def from_fractions(
        cls, fractions: Mapping[str, float], *, basis: str = "mole"
    ) -> Mixture:
        """Build a mixture from fractions that must sum to one."""

        normalized_basis = basis.strip().lower()
        if normalized_basis not in {"mole", "mass"}:
            raise InputError(
                "Composition basis must be mole or mass.",
                field="basis",
                details={"supported": ["mole", "mass"]},
            )
        if not fractions:
            raise InputError("Composition must contain at least one species.", field="fractions")
        entries = tuple((get_species(key), float(value)) for key, value in fractions.items())
        for species, fraction in entries:
            if fraction <= 0.0:
                raise DomainError(
                    f"Fraction for {species.key} must be greater than zero.",
                    field="fractions",
                )
        total = sum(fraction for _, fraction in entries)
        if not isclose(total, 1.0, rel_tol=0.0, abs_tol=1.0e-9):
            raise DomainError(
                "Composition fractions must sum to one.",
                code="composition_sum",
                field="fractions",
                details={"sum": total},
            )

        if normalized_basis == "mole":
            molecular_mass = sum(
                fraction * species.molecular_mass_kg_mol for species, fraction in entries
            )
            components = tuple(
                MixtureComponent(
                    species=species,
                    mole_fraction=fraction,
                    mass_fraction=fraction * species.molecular_mass_kg_mol / molecular_mass,
                )
                for species, fraction in entries
            )
        else:
            mole_denominator = sum(
                fraction / species.molecular_mass_kg_mol for species, fraction in entries
            )
            components = tuple(
                MixtureComponent(
                    species=species,
                    mole_fraction=(fraction / species.molecular_mass_kg_mol)
                    / mole_denominator,
                    mass_fraction=fraction,
                )
                for species, fraction in entries
            )
        return cls(components=components)

    @property
    def molecular_mass_kg_mol(self) -> float:
        return sum(
            component.mole_fraction * component.species.molecular_mass_kg_mol
            for component in self.components
        )

    @property
    def minimum_temperature_k(self) -> float:
        return max(component.species.minimum_temperature_k for component in self.components)

    @property
    def maximum_temperature_k(self) -> float:
        return min(component.species.maximum_temperature_k for component in self.components)

    def properties(
        self, temperature_k: float, *, pressure_pa: float = REFERENCE_PRESSURE_PA
    ) -> MixtureProperties:
        """Evaluate ideal-mixture properties including entropy of mixing."""

        if pressure_pa <= 0.0:
            raise DomainError("Mixture pressure must be greater than zero.", field="pressure_pa")
        species_properties = tuple(
            (component, component.species.properties(temperature_k))
            for component in self.components
        )
        cp_molar = sum(
            component.mole_fraction * properties.cp_molar_j_mol_k
            for component, properties in species_properties
        )
        enthalpy_molar = sum(
            component.mole_fraction * properties.enthalpy_molar_j_mol
            for component, properties in species_properties
        )
        standard_entropy_molar = sum(
            component.mole_fraction * properties.entropy_molar_j_mol_k
            for component, properties in species_properties
        )
        mixing_entropy = -UNIVERSAL_GAS_CONSTANT_J_MOL_K * sum(
            component.mole_fraction * log(component.mole_fraction)
            for component in self.components
        )
        pressure_entropy = -UNIVERSAL_GAS_CONSTANT_J_MOL_K * log(
            pressure_pa / REFERENCE_PRESSURE_PA
        )
        entropy_molar = standard_entropy_molar + mixing_entropy + pressure_entropy
        molecular_mass = self.molecular_mass_kg_mol
        gas_constant = UNIVERSAL_GAS_CONSTANT_J_MOL_K / molecular_mass
        cp_mass = cp_molar / molecular_mass
        cv_mass = cp_mass - gas_constant
        if cv_mass <= 0.0:
            raise DomainError(
                "Mixture produced non-positive cv.",
                code="invalid_thermodynamic_property",
                details={"temperature_k": temperature_k, "cv_j_kg_k": cv_mass},
            )
        return MixtureProperties(
            temperature_k=temperature_k,
            pressure_pa=pressure_pa,
            molecular_mass_kg_mol=molecular_mass,
            gas_constant_j_kg_k=gas_constant,
            cp_j_kg_k=cp_mass,
            cv_j_kg_k=cv_mass,
            gamma=cp_mass / cv_mass,
            enthalpy_j_kg=enthalpy_molar / molecular_mass,
            entropy_j_kg_k=entropy_molar / molecular_mass,
        )

    def temperature_from_enthalpy(
        self, enthalpy_j_kg: float, *, pressure_pa: float = REFERENCE_PRESSURE_PA
    ) -> float:
        """Invert mixture enthalpy over the common species range."""

        return solve_bisect_monotonic(
            lambda temperature: self.properties(
                temperature, pressure_pa=pressure_pa
            ).enthalpy_j_kg,
            enthalpy_j_kg,
            self.minimum_temperature_k,
            self.maximum_temperature_k,
            tolerance=1.0e-10,
        ).value

    def temperature_from_entropy(
        self, entropy_j_kg_k: float, *, pressure_pa: float = REFERENCE_PRESSURE_PA
    ) -> float:
        """Invert mixture entropy at constant pressure."""

        return solve_bisect_monotonic(
            lambda temperature: self.properties(
                temperature, pressure_pa=pressure_pa
            ).entropy_j_kg_k,
            entropy_j_kg_k,
            self.minimum_temperature_k,
            self.maximum_temperature_k,
            tolerance=1.0e-10,
        ).value

