"""Temperature-dependent ideal-gas species built from NASA Glenn data."""

from __future__ import annotations

import json
from dataclasses import dataclass
from functools import lru_cache
from math import isfinite
from pathlib import Path

from rocket_propulsion.core.errors import DomainError, InputError
from rocket_propulsion.core.numerics import solve_bisect_monotonic

from .nasa_polynomials import UNIVERSAL_GAS_CONSTANT_J_MOL_K, NasaPolynomialInterval

DATA_FILE = Path(__file__).parent / "data" / "nasa_glenn_2002.json"


@dataclass(frozen=True, slots=True)
class SpeciesProperties:
    """Molar and mass-specific ideal-gas properties at one temperature."""

    species: str
    temperature_k: float
    molecular_mass_kg_mol: float
    cp_molar_j_mol_k: float
    enthalpy_molar_j_mol: float
    entropy_molar_j_mol_k: float
    gas_constant_j_kg_k: float
    cp_j_kg_k: float
    cv_j_kg_k: float
    gamma: float
    enthalpy_j_kg: float
    entropy_j_kg_k: float


@dataclass(frozen=True, slots=True)
class Species:
    """One gas species and its piecewise NASA polynomial representation."""

    key: str
    name: str
    formula: tuple[tuple[str, float], ...]
    molecular_mass_kg_mol: float
    heat_of_formation_j_mol_298: float
    intervals: tuple[NasaPolynomialInterval, ...]
    source: str
    source_url: str

    @property
    def minimum_temperature_k(self) -> float:
        return self.intervals[0].minimum_temperature_k

    @property
    def maximum_temperature_k(self) -> float:
        return self.intervals[-1].maximum_temperature_k

    def _interval(self, temperature_k: float) -> NasaPolynomialInterval:
        if not isfinite(temperature_k):
            raise InputError("Temperature must be finite.", field="temperature_k")
        for index, interval in enumerate(self.intervals):
            if interval.contains(
                temperature_k, include_upper=index == len(self.intervals) - 1
            ):
                return interval
        raise DomainError(
            f"Temperature for {self.key} must be in "
            f"[{self.minimum_temperature_k:g}, {self.maximum_temperature_k:g}] K.",
            code="property_temperature_out_of_range",
            field="temperature_k",
            details={
                "species": self.key,
                "minimum_temperature_k": self.minimum_temperature_k,
                "maximum_temperature_k": self.maximum_temperature_k,
            },
        )

    def properties(self, temperature_k: float) -> SpeciesProperties:
        """Evaluate all ideal-gas properties without extrapolation."""

        interval = self._interval(temperature_k)
        cp_molar = interval.cp_molar_j_mol_k(temperature_k)
        enthalpy_molar = interval.enthalpy_molar_j_mol(temperature_k)
        entropy_molar = interval.entropy_molar_j_mol_k(temperature_k)
        gas_constant = UNIVERSAL_GAS_CONSTANT_J_MOL_K / self.molecular_mass_kg_mol
        cp_mass = cp_molar / self.molecular_mass_kg_mol
        cv_mass = cp_mass - gas_constant
        if cv_mass <= 0.0:
            raise DomainError(
                f"NASA polynomial for {self.key} produced non-positive cv.",
                code="invalid_thermodynamic_property",
                details={"temperature_k": temperature_k, "cv_j_kg_k": cv_mass},
            )
        return SpeciesProperties(
            species=self.key,
            temperature_k=temperature_k,
            molecular_mass_kg_mol=self.molecular_mass_kg_mol,
            cp_molar_j_mol_k=cp_molar,
            enthalpy_molar_j_mol=enthalpy_molar,
            entropy_molar_j_mol_k=entropy_molar,
            gas_constant_j_kg_k=gas_constant,
            cp_j_kg_k=cp_mass,
            cv_j_kg_k=cv_mass,
            gamma=cp_mass / cv_mass,
            enthalpy_j_kg=enthalpy_molar / self.molecular_mass_kg_mol,
            entropy_j_kg_k=entropy_molar / self.molecular_mass_kg_mol,
        )

    def temperature_from_enthalpy(self, enthalpy_j_kg: float) -> float:
        """Invert specific enthalpy over the species' documented range."""

        result = solve_bisect_monotonic(
            lambda temperature: self.properties(temperature).enthalpy_j_kg,
            enthalpy_j_kg,
            self.minimum_temperature_k,
            self.maximum_temperature_k,
            tolerance=1.0e-10,
        )
        return result.value

    def temperature_from_entropy(self, entropy_j_kg_k: float) -> float:
        """Invert standard-state specific entropy over the documented range."""

        result = solve_bisect_monotonic(
            lambda temperature: self.properties(temperature).entropy_j_kg_k,
            entropy_j_kg_k,
            self.minimum_temperature_k,
            self.maximum_temperature_k,
            tolerance=1.0e-10,
        )
        return result.value


def _parse_species(record: dict[str, object], source: dict[str, str]) -> Species:
    intervals = tuple(
        NasaPolynomialInterval(
            minimum_temperature_k=float(interval["minimum_temperature_k"]),
            maximum_temperature_k=float(interval["maximum_temperature_k"]),
            coefficients=tuple(float(value) for value in interval["coefficients"]),
        )
        for interval in record["intervals"]
    )
    return Species(
        key=str(record["key"]),
        name=str(record["name"]),
        formula=tuple((str(key), float(value)) for key, value in record["formula"].items()),
        molecular_mass_kg_mol=float(record["molecular_mass_g_mol"]) / 1000.0,
        heat_of_formation_j_mol_298=float(record["heat_of_formation_j_mol_298"]),
        intervals=intervals,
        source=source["title"],
        source_url=source["url"],
    )


@lru_cache(maxsize=1)
def species_catalog() -> dict[str, Species]:
    """Load and validate the bundled, source-labelled species catalog once."""

    document = json.loads(DATA_FILE.read_text(encoding="utf-8"))
    source = document["source"]
    result = {
        record["key"].upper(): _parse_species(record, source)
        for record in document["species"]
    }
    if len(result) != len(document["species"]):
        raise RuntimeError("Duplicate species keys in thermochemistry data file.")
    return result


def get_species(key: str) -> Species:
    """Return a species by case-insensitive key."""

    normalized = key.strip().upper()
    try:
        return species_catalog()[normalized]
    except KeyError as error:
        raise InputError(
            f"Unknown species: {key!r}.",
            code="unknown_species",
            field="species",
            details={"supported": sorted(species_catalog())},
        ) from error

