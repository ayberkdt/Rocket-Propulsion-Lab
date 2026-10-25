"""Temperature-dependent ideal-gas species built from NASA Glenn data.

Species records and their polynomial coefficients are loaded from the bundled
``data/nasa_glenn_2002.json`` file, whose source title and URL are retained
on every :class:`Species` [1]_.

References
----------
.. [1] NASA/TP-2002-211556. NASA Glenn, "NASA Glenn Coefficients for
       Calculating Thermodynamic Properties of Individual Species".
       https://ntrs.nasa.gov/citations/20020085330
"""

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
    """Molar and mass-specific ideal-gas properties at one temperature [1]_.

    Attributes
    ----------
    species : str
        Species key.
    temperature_k : float
        Temperature in K.
    molecular_mass_kg_mol : float
        Molar mass in kg/mol.
    cp_molar_j_mol_k, enthalpy_molar_j_mol, entropy_molar_j_mol_k : float
        Standard-state molar heat capacity, enthalpy and entropy.
    gas_constant_j_kg_k : float
        Specific gas constant ``R/M``.
    cp_j_kg_k, cv_j_kg_k : float
        Specific heats, with ``cv = cp - R/M``.
    gamma : float
        ``cp/cv``.
    enthalpy_j_kg, entropy_j_kg_k : float
        Mass-specific standard-state enthalpy and entropy.

    References
    ----------
    .. [1] NASA/TP-2002-211556. NASA Glenn, "NASA Glenn Coefficients for
           Calculating Thermodynamic Properties of Individual Species".
           https://ntrs.nasa.gov/citations/20020085330
    """

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
    """One gas species and its piecewise NASA polynomial representation.

    Attributes
    ----------
    key, name : str
        Catalog key and descriptive name.
    formula : tuple of (str, float)
        Element symbols and atom counts.
    molecular_mass_kg_mol : float
        Molar mass in kg/mol.
    heat_of_formation_j_mol_298 : float
        Heat of formation at the reference temperature, J/mol.
    intervals : tuple of NasaPolynomialInterval
        Contiguous temperature intervals in increasing order [1]_.
    source, source_url : str
        Title and URL of the data source.

    References
    ----------
    .. [1] NASA/TP-2002-211556. NASA Glenn, "NASA Glenn Coefficients for
           Calculating Thermodynamic Properties of Individual Species".
           https://ntrs.nasa.gov/citations/20020085330
    """

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
        """Lower bound of the first polynomial interval, K."""
        return self.intervals[0].minimum_temperature_k

    @property
    def maximum_temperature_k(self) -> float:
        """Upper bound of the last polynomial interval, K."""
        return self.intervals[-1].maximum_temperature_k

    def _interval(self, temperature_k: float) -> NasaPolynomialInterval:
        """Return the interval owning ``temperature_k``; never extrapolate.

        Raises
        ------
        InputError
            If the temperature is not finite.
        DomainError
            If it lies outside the species range
            (code ``property_temperature_out_of_range``).
        """
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
        """Evaluate all ideal-gas properties without extrapolation.

        Parameters
        ----------
        temperature_k : float
            Temperature in K within the species range.

        Returns
        -------
        SpeciesProperties
            Molar and mass-specific properties.

        Raises
        ------
        InputError
            If the temperature is not finite.
        DomainError
            If the temperature is outside the polynomial range, or the
            polynomial yields a non-positive ``cv``.

        Notes
        -----
        Molar properties come from the NASA Glenn polynomials [1]_; mass
        quantities divide by the molar mass, ``R/M`` uses the registered molar
        gas constant, ``cv = cp - R/M`` and ``gamma = cp/cv`` for an ideal
        gas.

        References
        ----------
        .. [1] NASA/TP-2002-211556. NASA Glenn, "NASA Glenn Coefficients for
               Calculating Thermodynamic Properties of Individual Species".
               https://ntrs.nasa.gov/citations/20020085330
        """

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
        """Invert specific enthalpy over the species' documented range.

        Parameters
        ----------
        enthalpy_j_kg : float
            Target mass-specific enthalpy in J/kg.

        Returns
        -------
        float
            Temperature in K, found by monotonic bisection with tolerance
            1e-10 on the NASA Glenn polynomial enthalpy [1]_.

        Raises
        ------
        DomainError
            If the target is not bracketed by the species range.

        References
        ----------
        .. [1] NASA/TP-2002-211556. NASA Glenn, "NASA Glenn Coefficients for
               Calculating Thermodynamic Properties of Individual Species".
               https://ntrs.nasa.gov/citations/20020085330
        """

        result = solve_bisect_monotonic(
            lambda temperature: self.properties(temperature).enthalpy_j_kg,
            enthalpy_j_kg,
            self.minimum_temperature_k,
            self.maximum_temperature_k,
            tolerance=1.0e-10,
        )
        return result.value

    def temperature_from_entropy(self, entropy_j_kg_k: float) -> float:
        """Invert standard-state specific entropy over the documented range.

        Parameters
        ----------
        entropy_j_kg_k : float
            Target mass-specific standard-state entropy in J/(kg K).

        Returns
        -------
        float
            Temperature in K, found by monotonic bisection with tolerance
            1e-10 on the NASA Glenn polynomial entropy [1]_.

        Raises
        ------
        DomainError
            If the target is not bracketed by the species range.

        References
        ----------
        .. [1] NASA/TP-2002-211556. NASA Glenn, "NASA Glenn Coefficients for
               Calculating Thermodynamic Properties of Individual Species".
               https://ntrs.nasa.gov/citations/20020085330
        """

        result = solve_bisect_monotonic(
            lambda temperature: self.properties(temperature).entropy_j_kg_k,
            entropy_j_kg_k,
            self.minimum_temperature_k,
            self.maximum_temperature_k,
            tolerance=1.0e-10,
        )
        return result.value


def _parse_species(record: dict[str, object], source: dict[str, str]) -> Species:
    """Build a :class:`Species` from one JSON record; molar mass is given in g/mol."""
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
    """Load and validate the bundled, source-labelled species catalog once [1]_.

    Returns
    -------
    dict of str to Species
        Species keyed by upper-case key; the result is cached.

    Raises
    ------
    RuntimeError
        If the data file contains duplicate species keys.

    References
    ----------
    .. [1] NASA/TP-2002-211556. NASA Glenn, "NASA Glenn Coefficients for
           Calculating Thermodynamic Properties of Individual Species".
           https://ntrs.nasa.gov/citations/20020085330
    """

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
    """Return a species by case-insensitive key.

    Parameters
    ----------
    key : str
        Species key; surrounding whitespace and case are ignored.

    Returns
    -------
    Species
        Catalog entry.

    Raises
    ------
    InputError
        If the key is unknown (code ``unknown_species``); the supported keys
        are listed in the error details.
    """

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

