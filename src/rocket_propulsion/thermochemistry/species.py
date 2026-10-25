"""Temperature-dependent ideal-gas species and condensed reactants from NASA data.

The ten combustion species of ``data/nasa_glenn_2002.json`` come from
NASA/TP-2002-211556 [1]_. Additional gas-phase products and the condensed
reactants with an assigned enthalpy (liquid oxygen, liquid hydrogen, liquid
methane, RP-1, NTO, MMH) come from ``data/nasa_cea_thermo.json``, extracted by
``scripts/sync_nasa_thermo.py`` from NASA's CEA ``thermo.inp`` at a pinned
commit [2]_. Every record keeps the title and URL of its own source.

References
----------
.. [1] NASA/TP-2002-211556. NASA Glenn, "NASA Glenn Coefficients for
       Calculating Thermodynamic Properties of Individual Species".
       https://ntrs.nasa.gov/citations/20020085330
.. [2] NASA-CEA-THERMO. NASA Glenn Research Center, CEA thermodynamic database
       ``data/thermo.inp``, commit 2f79a647e85737742fb08a127da00d1ea1d069da.
       https://github.com/nasa/cea
"""

from __future__ import annotations

import json
from dataclasses import dataclass
from functools import lru_cache
from math import isfinite
from pathlib import Path
from typing import TypedDict

from rocket_propulsion.core.errors import DomainError, InputError
from rocket_propulsion.core.numerics import solve_bisect_monotonic

from .nasa_polynomials import (
    UNIVERSAL_GAS_CONSTANT_J_MOL_K,
    NasaIntervalRecord,
    NasaPolynomialInterval,
    interval_from_record,
)

DATA_FILE = Path(__file__).parent / "data" / "nasa_glenn_2002.json"
CEA_DATA_FILE = Path(__file__).parent / "data" / "nasa_cea_thermo.json"


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


class _SpeciesRecord(TypedDict):
    """One gas-species record of the packaged JSON data files."""

    key: str
    name: str
    formula: dict[str, float]
    molecular_mass_g_mol: float
    heat_of_formation_j_mol_298: float
    intervals: list[NasaIntervalRecord]


def _parse_species(record: _SpeciesRecord, source: dict[str, str]) -> Species:
    """Build a :class:`Species` from one JSON record; molar mass is given in g/mol."""
    intervals = tuple(interval_from_record(interval) for interval in record["intervals"])
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
    """Load and validate the bundled, source-labelled species catalog once [1]_, [2]_.

    Returns
    -------
    dict of str to Species
        Species keyed by upper-case key; the result is cached.

    Raises
    ------
    RuntimeError
        If a species key occurs in more than one record.

    References
    ----------
    .. [1] NASA/TP-2002-211556. NASA Glenn, "NASA Glenn Coefficients for
           Calculating Thermodynamic Properties of Individual Species".
           https://ntrs.nasa.gov/citations/20020085330
    .. [2] NASA-CEA-THERMO. NASA Glenn Research Center, CEA thermodynamic
           database ``data/thermo.inp``. https://github.com/nasa/cea
    """

    result: dict[str, Species] = {}
    count = 0
    for path in (DATA_FILE, CEA_DATA_FILE):
        document = json.loads(path.read_text(encoding="utf-8"))
        source = document["source"]
        for record in document["species"]:
            result[record["key"].upper()] = _parse_species(record, source)
        count += len(document["species"])
    if len(result) != count:
        raise RuntimeError("Duplicate species keys in thermochemistry data files.")
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



@dataclass(frozen=True, slots=True)
class CondensedReactant:
    """A condensed reactant with an assigned enthalpy at one temperature [1]_.

    NASA CEA describes cryogenic and storable liquid propellants by their
    enthalpy at the storage temperature instead of a polynomial; the reactant
    enters only the element and energy balance, never the product set [2]_.

    Attributes
    ----------
    key, name : str
        Catalog key (the CEA name, for example ``"O2(L)"``) and description.
    formula : tuple of (str, float)
        Element symbols and atom counts; RP-1 is the CEA formula ``CH1.95``.
    molecular_mass_kg_mol : float
        Molar mass in kg/mol.
    enthalpy_j_mol : float
        Assigned molar enthalpy at ``temperature_k`` on the NASA scale, J/mol.
    temperature_k : float
        Temperature at which the enthalpy is assigned, K.
    source, source_url : str
        Title and URL of the data source.

    References
    ----------
    .. [1] NASA-CEA-THERMO. NASA Glenn Research Center, CEA thermodynamic
           database ``data/thermo.inp``. https://github.com/nasa/cea
    .. [2] NASA RP-1311, "Computer Program for Calculation of Complex Chemical
           Equilibrium Compositions and Applications, Part II: Users Manual and
           Program Description". https://ntrs.nasa.gov/citations/19960044559
    """

    key: str
    name: str
    formula: tuple[tuple[str, float], ...]
    molecular_mass_kg_mol: float
    enthalpy_j_mol: float
    temperature_k: float
    source: str
    source_url: str


@lru_cache(maxsize=1)
def condensed_reactant_catalog() -> dict[str, CondensedReactant]:
    """Load the bundled assigned-enthalpy reactants once [1]_.

    Returns
    -------
    dict of str to CondensedReactant
        Reactants keyed by upper-case CEA name; the result is cached.

    References
    ----------
    .. [1] NASA-CEA-THERMO. NASA Glenn Research Center, CEA thermodynamic
           database ``data/thermo.inp``. https://github.com/nasa/cea
    """

    document = json.loads(CEA_DATA_FILE.read_text(encoding="utf-8"))
    source = document["source"]
    return {
        record["key"].upper(): CondensedReactant(
            key=str(record["key"]),
            name=str(record["name"]),
            formula=tuple((str(key), float(value)) for key, value in record["formula"].items()),
            molecular_mass_kg_mol=float(record["molecular_mass_g_mol"]) / 1000.0,
            enthalpy_j_mol=float(record["assigned_enthalpy_j_mol"]),
            temperature_k=float(record["assigned_temperature_k"]),
            source=str(source["title"]),
            source_url=str(source["url"]),
        )
        for record in document["reactants"]
    }


def get_reactant(key: str) -> Species | CondensedReactant:
    """Return a gas species or a condensed reactant by case-insensitive key.

    Parameters
    ----------
    key : str
        Species key (``"O2"``) or condensed reactant key (``"O2(L)"``).

    Returns
    -------
    Species or CondensedReactant
        The condensed reactant with that key, otherwise the gas species; the
        two key sets are disjoint because condensed keys carry a phase or a
        mixture name.

    Raises
    ------
    InputError
        If the key is neither a species nor a condensed reactant
        (code ``unknown_species``).
    """

    normalized = key.strip().upper()
    reactant = condensed_reactant_catalog().get(normalized)
    if reactant is not None:
        return reactant
    return get_species(key)
