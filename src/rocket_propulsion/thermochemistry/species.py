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

from collections.abc import Mapping
from dataclasses import dataclass
from fractions import Fraction
from functools import lru_cache
from pathlib import Path
from typing import Any, TypedDict, cast

from rocket_propulsion.core.errors import DomainError, InputError
from rocket_propulsion.core.numerics import solve_bisect_monotonic
from rocket_propulsion.core.validation import (
    coerce_fields,
    numeric_range_guard,
    read_json_document,
    real,
    text,
)

from ._records import _caloric_identity, _owned_catalog, _owned_formula, _property_identity
from .nasa_polynomials import (
    UNIVERSAL_GAS_CONSTANT_J_MOL_K,
    NasaIntervalRecord,
    NasaPolynomialInterval,
    _export,
    _owned_intervals,
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

    def __post_init__(self) -> None:
        coerce_fields(self, 'Species properties')
        if self.temperature_k <= 0 or self.molecular_mass_kg_mol <= 0:
            raise DomainError('Species temperature and molar mass must be positive.')
        object.__setattr__(self, 'species', self.species.strip())
        _caloric_identity(self.cp_j_kg_k, self.cv_j_kg_k, self.gas_constant_j_kg_k, self.gamma)
        mass = self.molecular_mass_kg_mol
        _property_identity(self.gas_constant_j_kg_k * mass, UNIVERSAL_GAS_CONSTANT_J_MOL_K, 'R')
        _property_identity(self.cp_j_kg_k * mass, self.cp_molar_j_mol_k, 'molar cp')
        _property_identity(self.enthalpy_j_kg * mass, self.enthalpy_molar_j_mol, 'molar enthalpy')
        _property_identity(self.entropy_j_kg_k * mass, self.entropy_molar_j_mol_k, 'molar entropy')


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

    def __post_init__(self) -> None:
        coerce_fields(self, 'Species')
        if self.molecular_mass_kg_mol <= 0:
            raise DomainError('Species molar mass must be positive.', field='molecular_mass_kg_mol')
        object.__setattr__(self, 'key', self.key.strip())
        object.__setattr__(self, 'formula', _owned_formula(self.formula))
        object.__setattr__(self, 'intervals', _owned_intervals(self.intervals))

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
        temperature_k = real(temperature_k, 'temperature_k')
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

    @numeric_range_guard('temperature_k')
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
        gas. Every field is formed from the exact retained terms and rounded
        once, so a positive ``cv`` is decided before display rounding.

        References
        ----------
        .. [1] NASA/TP-2002-211556. NASA Glenn, "NASA Glenn Coefficients for
               Calculating Thermodynamic Properties of Individual Species".
               https://ntrs.nasa.gov/citations/20020085330
        """

        return _species_properties(self, real(temperature_k, 'temperature_k'))

    @numeric_range_guard('temperature_k')
    def _binary64_properties(self, temperature_k: float) -> SpeciesProperties:
        """Evaluate a species state in binary64 for iterative solver internals.

        Same range checks and record identities as :meth:`properties`, but
        with main's original floating evaluation instead of one exact rounding
        per field. Results agree to a few ULPs.
        """
        temperature_k = real(temperature_k, 'temperature_k')
        cp_molar, enthalpy_molar, entropy_molar = self._interval(temperature_k)._binary64_molar(temperature_k)
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
            real(enthalpy_j_kg, 'enthalpy_j_kg'),
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
            real(entropy_j_kg_k, 'entropy_j_kg_k'),
            self.minimum_temperature_k,
            self.maximum_temperature_k,
            tolerance=1.0e-10,
        )
        return result.value


@lru_cache(maxsize=4096)
def _species_properties(species: Species, temperature_k: float) -> SpeciesProperties:
    """Evaluate one owned species state; equal immutable inputs share a result."""
    interval = species._interval(temperature_k)
    (cp, _), (h, hlog), (s, slog) = interval._terms(temperature_k)
    mass = Fraction(species.molecular_mass_kg_mol)
    r = Fraction(UNIVERSAL_GAS_CONSTANT_J_MOL_K)
    cv = cp - r
    if cv <= 0:
        raise DomainError(
            f"NASA polynomial for {species.key} produced non-positive cv.",
            code="invalid_thermodynamic_property",
            details={"temperature_k": temperature_k},
        )
    return SpeciesProperties(
        species=species.key,
        temperature_k=temperature_k,
        molecular_mass_kg_mol=species.molecular_mass_kg_mol,
        cp_molar_j_mol_k=_export(cp, 0, temperature_k, 'cp_molar_j_mol_k'),
        enthalpy_molar_j_mol=_export(h, hlog, temperature_k, 'enthalpy_molar_j_mol'),
        entropy_molar_j_mol_k=_export(s, slog, temperature_k, 'entropy_molar_j_mol_k'),
        gas_constant_j_kg_k=_export(r/mass, 0, temperature_k, 'gas_constant_j_kg_k'),
        cp_j_kg_k=_export(cp/mass, 0, temperature_k, 'cp_j_kg_k'),
        cv_j_kg_k=_export(cv/mass, 0, temperature_k, 'cv_j_kg_k'),
        gamma=_export(cp/cv, 0, temperature_k, 'gamma'),
        enthalpy_j_kg=_export(h/mass, hlog/mass, temperature_k, 'enthalpy_j_kg'),
        entropy_j_kg_k=_export(s/mass, slog/mass, temperature_k, 'entropy_j_kg_k'),
    )


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
        key=record["key"],
        name=record["name"],
        formula=tuple(record["formula"].items()),
        molecular_mass_kg_mol=real(record["molecular_mass_g_mol"], 'molecular_mass_g_mol') / 1000.0,
        heat_of_formation_j_mol_298=record["heat_of_formation_j_mol_298"],
        intervals=intervals,
        source=source["title"],
        source_url=source["url"],
    )


_COEFFICIENT_ORDER = ['a1', 'a2', 'a3', 'a4', 'a5', 'a6', 'a7', 'b1', 'b2']


def _catalogue_document(path: Path, *, defining_constant: bool) -> dict[str, Any]:
    """Read one packaged catalogue under the strict JSON rules and check its model.

    Duplicate members, non-finite tokens, an unsupported schema, coefficient
    order, standard pressure or (when declared or required) defining gas
    constant are refused. This checks packaged data integrity only.
    """

    def build(document):
        if type(document.get('schema_version')) is not int or document['schema_version'] != 1:
            raise InputError('Unsupported thermochemistry catalogue schema.', field='schema_version')
        model = document['model']
        if (model['coefficient_order'] != _COEFFICIENT_ORDER
                or type(model['reference_pressure_pa']) not in (int, float)
                or model['reference_pressure_pa'] != 100_000):
            raise InputError('Unsupported NASA9 coefficient order or reference pressure.', field='model')
        declared = defining_constant or 'gas_constant_j_mol_k' in model
        if declared and real(model['gas_constant_j_mol_k'],
                             'model.gas_constant_j_mol_k') != UNIVERSAL_GAS_CONSTANT_J_MOL_K:
            raise InputError('Catalogue gas constant differs from the fit model.', field='model')
        text(document['source']['title'], 'source.title')
        text(document['source']['url'], 'source.url')
        if not isinstance(document['species'], list) or not document['species']:
            raise InputError('Catalogue species must be a nonempty array.', field='species')
        return document

    return cast(dict[str, Any], read_json_document(path.read_text(encoding='utf-8'), path.name, build))


@lru_cache(maxsize=1)
def _species_catalog() -> Mapping[str, Species]:
    """Load and validate the bundled, source-labelled species catalog once [1]_, [2]_.

    Returns
    -------
    mapping of str to Species
        Private immutable table keyed by canonical upper-case lookup keys;
        source identifiers retain their printed spelling.

    Raises
    ------
    InputError
        If a document violates the strict JSON or catalogue schema, or a
        species key occurs in more than one record.

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
    for path, defining in ((DATA_FILE, True), (CEA_DATA_FILE, False)):
        document = _catalogue_document(path, defining_constant=defining)
        source = document["source"]
        for record in document["species"]:
            if not isinstance(record, Mapping):
                raise InputError('Malformed catalogue species record.', field='species')
            try:
                parsed = _parse_species(cast(_SpeciesRecord, record), source)
            except (KeyError, TypeError, AttributeError) as error:
                raise InputError('Malformed catalogue species record.', field='species') from error
            result[parsed.key.upper()] = parsed
        count += len(document["species"])
    if len(result) != count:
        raise InputError("Duplicate species keys in thermochemistry data files.", field='species')
    return _owned_catalog(result, Species, 'species')


def species_catalog() -> dict[str, Species]:
    """Return an independent mapping of the bundled immutable species [1]_.

    The validated internal catalogue is cached; editing this dictionary
    cannot affect lookup or any later caller. Nested formulae and intervals
    are owned immutable tuples.

    References
    ----------
    .. [1] NASA/TP-2002-211556. NASA Glenn, "NASA Glenn Coefficients for
           Calculating Thermodynamic Properties of Individual Species".
           https://ntrs.nasa.gov/citations/20020085330
    """
    return dict(_species_catalog())


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

    normalized = text(key, 'species').strip().upper()
    try:
        return _species_catalog()[normalized]
    except KeyError as error:
        raise InputError(
            f"Unknown species: {key!r}.",
            code="unknown_species",
            field="species",
            details={"supported": sorted(_species_catalog())},
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

    def __post_init__(self) -> None:
        coerce_fields(self, 'Condensed reactant')
        if self.molecular_mass_kg_mol <= 0 or self.temperature_k <= 0:
            raise DomainError('Condensed-reactant molar mass and assigned temperature must be positive.')
        object.__setattr__(self, 'key', self.key.strip())
        object.__setattr__(self, 'formula', _owned_formula(self.formula))


@lru_cache(maxsize=1)
def _condensed_reactant_catalog() -> Mapping[str, CondensedReactant]:
    """Load the bundled assigned-enthalpy reactants once [1]_.

    Returns
    -------
    mapping of str to CondensedReactant
        Private immutable table keyed by upper-case CEA name.

    References
    ----------
    .. [1] NASA-CEA-THERMO. NASA Glenn Research Center, CEA thermodynamic
           database ``data/thermo.inp``. https://github.com/nasa/cea
    """

    document = _catalogue_document(CEA_DATA_FILE, defining_constant=False)
    source = document["source"]
    result: dict[str, CondensedReactant] = {}
    for record in document['reactants']:
        parsed = CondensedReactant(
            key=record["key"],
            name=record["name"],
            formula=tuple(record["formula"].items()),
            molecular_mass_kg_mol=real(record["molecular_mass_g_mol"], 'molecular_mass_g_mol') / 1000.0,
            enthalpy_j_mol=record["assigned_enthalpy_j_mol"],
            temperature_k=record["assigned_temperature_k"],
            source=source["title"],
            source_url=source["url"],
        )
        key = parsed.key.upper()
        if key in result:
            raise RuntimeError('Duplicate condensed-reactant keys.')
        result[key] = parsed
    return _owned_catalog(result, CondensedReactant, 'reactants')


def condensed_reactant_catalog() -> dict[str, CondensedReactant]:
    """Return an independent mapping of assigned-enthalpy reactants [1]_.

    Editing the returned dictionary cannot change the cached catalogue.

    References
    ----------
    .. [1] NASA-CEA-THERMO. NASA Glenn Research Center, CEA thermodynamic
           database ``data/thermo.inp``. https://github.com/nasa/cea
    """
    return dict(_condensed_reactant_catalog())


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

    normalized = text(key, 'species').strip().upper()
    reactant = _condensed_reactant_catalog().get(normalized)
    if reactant is not None:
        return reactant
    return get_species(key)
