"""Ideal-gas mixture properties with explicit composition bases.

Mixture properties are mole-fraction averages of NASA Glenn species
properties [2]_ plus the ideal entropy of mixing and the pressure correction
to the 1 bar standard state [1]_.

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
from dataclasses import dataclass
from math import isclose, log

from rocket_propulsion.core.errors import DomainError, InputError
from rocket_propulsion.core.numerics import solve_bisect_monotonic

from .nasa_polynomials import UNIVERSAL_GAS_CONSTANT_J_MOL_K
from .species import Species, get_species

REFERENCE_PRESSURE_PA = 100_000.0


@dataclass(frozen=True, slots=True)
class MixtureComponent:
    """A resolved species with consistent mole and mass fractions.

    Attributes
    ----------
    species : Species
        Catalog species.
    mole_fraction, mass_fraction : float
        Fractions of the component in the mixture.
    """

    species: Species
    mole_fraction: float
    mass_fraction: float


@dataclass(frozen=True, slots=True)
class MixtureProperties:
    """Temperature-dependent properties of an ideal gas mixture.

    Attributes
    ----------
    temperature_k, pressure_pa : float
        State at which the properties are evaluated.
    molecular_mass_kg_mol : float
        Mixture molar mass.
    gas_constant_j_kg_k : float
        Specific gas constant ``R/M``.
    cp_j_kg_k, cv_j_kg_k : float
        Frozen-composition specific heats.
    gamma : float
        Frozen ratio of specific heats ``cp/cv``.
    enthalpy_j_kg : float
        Mass-specific enthalpy.
    entropy_j_kg_k : float
        Mass-specific entropy including mixing and pressure terms [1]_.

    References
    ----------
    .. [1] NASA RP-1311, "Computer Program for Calculation of Complex Chemical
           Equilibrium Compositions and Applications, Part I: Analysis".
           https://ntrs.nasa.gov/citations/19950013764
    """

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
    """An immutable ideal-gas mixture defined by resolved components [1]_.

    Attributes
    ----------
    components : tuple of MixtureComponent
        Components with consistent mole and mass fractions.

    References
    ----------
    .. [1] NASA RP-1311, "Computer Program for Calculation of Complex Chemical
           Equilibrium Compositions and Applications, Part I: Analysis".
           https://ntrs.nasa.gov/citations/19950013764
    """

    components: tuple[MixtureComponent, ...]

    @classmethod
    def from_fractions(
        cls, fractions: Mapping[str, float], *, basis: str = "mole"
    ) -> Mixture:
        """Build a mixture from fractions that must sum to one.

        Parameters
        ----------
        fractions : mapping of str to float
            Positive fraction per species key.
        basis : str, default "mole"
            ``"mole"`` or ``"mass"`` (case-insensitive).

        Returns
        -------
        Mixture
            Mixture with both mole and mass fractions resolved.

        Raises
        ------
        InputError
            If the basis is unknown, the composition is empty or a species is
            unknown.
        DomainError
            If a fraction is not positive or the fractions do not sum to one
            within 1e-9 (code ``composition_sum``).
        """

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
        """Mole-fraction-weighted molar mass, kg/mol."""
        return sum(
            component.mole_fraction * component.species.molecular_mass_kg_mol
            for component in self.components
        )

    @property
    def minimum_temperature_k(self) -> float:
        """Lowest temperature valid for every component, K."""
        return max(component.species.minimum_temperature_k for component in self.components)

    @property
    def maximum_temperature_k(self) -> float:
        """Highest temperature valid for every component, K."""
        return min(component.species.maximum_temperature_k for component in self.components)

    def properties(
        self, temperature_k: float, *, pressure_pa: float = REFERENCE_PRESSURE_PA
    ) -> MixtureProperties:
        """Evaluate ideal-mixture properties including entropy of mixing.

        Parameters
        ----------
        temperature_k : float
            Temperature in K, within every component's range.
        pressure_pa : float, default REFERENCE_PRESSURE_PA
            Mixture pressure in Pa, positive.

        Returns
        -------
        MixtureProperties
            Frozen-composition properties at the state.

        Raises
        ------
        DomainError
            If the pressure is not positive, a component is evaluated outside
            its range, or the mixture ``cv`` is not positive.

        Notes
        -----
        Molar ``cp``, ``H`` and standard-state ``S`` are mole-fraction
        averages of the species values [2]_. The ideal-mixture entropy adds
        ``-R sum(x_j ln x_j) - R ln(p/p_ref)`` with ``p_ref`` the 1 bar
        standard state [1]_.

        References
        ----------
        .. [1] NASA RP-1311, "Computer Program for Calculation of Complex Chemical
               Equilibrium Compositions and Applications, Part I: Analysis".
               https://ntrs.nasa.gov/citations/19950013764
        .. [2] NASA/TP-2002-211556. NASA Glenn, "NASA Glenn Coefficients for
               Calculating Thermodynamic Properties of Individual Species".
               https://ntrs.nasa.gov/citations/20020085330
        """

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
        """Invert mixture enthalpy over the common species range.

        Parameters
        ----------
        enthalpy_j_kg : float
            Target mass-specific enthalpy in J/kg.
        pressure_pa : float, default REFERENCE_PRESSURE_PA
            Pressure in Pa (enthalpy of an ideal gas does not depend on it).

        Returns
        -------
        float
            Temperature in K by monotonic bisection (tolerance 1e-10) at fixed
            composition, using :meth:`properties` [1]_.

        Raises
        ------
        DomainError
            If the target is not bracketed by the common range.

        References
        ----------
        .. [1] NASA RP-1311, "Computer Program for Calculation of Complex Chemical
               Equilibrium Compositions and Applications, Part I: Analysis".
               https://ntrs.nasa.gov/citations/19950013764
        """

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
        """Invert mixture entropy at constant pressure.

        Parameters
        ----------
        entropy_j_kg_k : float
            Target mass-specific entropy in J/(kg K).
        pressure_pa : float, default REFERENCE_PRESSURE_PA
            Pressure in Pa.

        Returns
        -------
        float
            Temperature in K by monotonic bisection (tolerance 1e-10) at fixed
            composition, using :meth:`properties` [1]_.

        Raises
        ------
        DomainError
            If the target is not bracketed by the common range.

        References
        ----------
        .. [1] NASA RP-1311, "Computer Program for Calculation of Complex Chemical
               Equilibrium Compositions and Applications, Part I: Analysis".
               https://ntrs.nasa.gov/citations/19950013764
        """

        return solve_bisect_monotonic(
            lambda temperature: self.properties(
                temperature, pressure_pa=pressure_pa
            ).entropy_j_kg_k,
            entropy_j_kg_k,
            self.minimum_temperature_k,
            self.maximum_temperature_k,
            tolerance=1.0e-10,
        ).value

