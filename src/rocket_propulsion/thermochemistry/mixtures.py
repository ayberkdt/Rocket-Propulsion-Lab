"""Ideal-gas mixture properties with explicit composition bases.

Mixture properties are mole-fraction averages of NASA Glenn species
properties [2]_ plus the ideal entropy of mixing and the pressure correction
to the 1 bar standard state [1]_. Accepted raw fractions and their basis are
retained and normalized exactly once; forward properties are formed from
exact rational terms with enclosed logarithms (see ``mixture_numerics``).
These are rounding bounds on the retained model, not physical error bounds,
and they do not qualify real gases, equilibrium composition or inverse
accuracy [3]_.

References
----------
.. [1] NASA RP-1311, "Computer Program for Calculation of Complex Chemical
       Equilibrium Compositions and Applications, Part I: Analysis".
       https://ntrs.nasa.gov/citations/19950013764
.. [2] NASA/TP-2002-211556. NASA Glenn, "NASA Glenn Coefficients for
       Calculating Thermodynamic Properties of Individual Species".
       https://ntrs.nasa.gov/citations/20020085330
.. [3] CANTERA-32-IDEAL-MIXTURE. Cantera 3.2, IdealGasPhase standard and
       partial molar relations.
       https://cantera.org/3.2/cxx/d7/dfa/classCantera_1_1IdealGasPhase.html
"""

from __future__ import annotations

from collections.abc import Mapping
from collections.abc import Set as AbstractSet
from dataclasses import dataclass, field, fields
from itertools import islice
from math import fsum, isclose, log, log1p
from typing import NamedTuple

from rocket_propulsion.core.errors import DomainError, InputError
from rocket_propulsion.core.numerics import solve_bisect_monotonic
from rocket_propulsion.core.validation import (
    _ordered_items,
    coerce_fields,
    numeric_range_guard,
    real,
    text,
)

from ._records import _caloric_identity, _property_identity
from .mixture_numerics import _coalesced, _composition, _forward, _ledger, _rational
from .nasa_polynomials import UNIVERSAL_GAS_CONSTANT_J_MOL_K
from .species import Species, get_species

REFERENCE_PRESSURE_PA = 100_000.0
NORMALIZATION_POLICY = 'exact_normalization_of_retained_fraction_inputs_1e-9_v1'
MAXIMUM_COMPONENTS = 64


def _owned_species(value):
    if not isinstance(value,Species):
        raise DomainError('Mixture component requires a Species record.',field='species')
    return Species(**{f.name:getattr(value,f.name) for f in fields(Species)})


def _bounded(values, label):
    if isinstance(values,(str,bytes,bytearray,memoryview,Mapping,AbstractSet)):
        raise DomainError('Mixture records must be an ordered iterable.',field=label)
    try: rows = tuple(islice(iter(values),65))
    except (TypeError,ValueError,KeyError,IndexError,OverflowError) as error:
        raise DomainError('Invalid mixture record iterable.',field=label) from error
    if not 1 <= len(rows) <= MAXIMUM_COMPONENTS:
        raise DomainError('Mixture requires 1..64 components.',code='mixture_component_resource',field=label)
    return rows


def _basis(value, label='composition_basis'):
    result = text(value,label).strip().lower()
    if result not in {'mole','mass'}:
        raise InputError('Composition basis must be mole or mass.',field=label,
                         details={'supported': ['mole', 'mass']})
    return result


def _positive(value,label):
    number = real(value,label)
    if number <= 0:
        raise DomainError('Mixture thermodynamic control must be positive.',field=label)
    return number


def _canonical(species):
    return species.key.upper()


class _Binary64State(NamedTuple):
    """Frozen-mixture state from binary64 evaluation, for solver internals only."""

    molecular_mass_kg_mol: float
    gas_constant_j_kg_k: float
    cp_j_kg_k: float
    cv_j_kg_k: float
    gamma: float
    enthalpy_j_kg: float
    entropy_j_kg_k: float


@dataclass(frozen=True, slots=True)
class MixtureComponent:
    """A resolved species with consistent mole and mass fractions.

    Attributes
    ----------
    species : Species
        Owned copy of a catalog species.
    mole_fraction, mass_fraction : float
        Binary64 displays of the component's normalized fractions.
    """

    species: Species
    mole_fraction: float
    mass_fraction: float

    def __post_init__(self) -> None:
        coerce_fields(self, 'Mixture component')
        if not (0 < self.mole_fraction <= 1 and 0 < self.mass_fraction <= 1):
            raise DomainError('Component fractions must be in (0, 1].', field='fractions')
        object.__setattr__(self, 'species', _owned_species(self.species))

    @classmethod
    def _sharing(cls, species: Species, mole_fraction: float, mass_fraction: float) -> MixtureComponent:
        """Build a component around a Species this module already owns, without another copy."""
        component = object.__new__(cls)
        for name, value in (('mole_fraction', mole_fraction), ('mass_fraction', mass_fraction)):
            object.__setattr__(component, name, real(value, name))
        if not (0 < component.mole_fraction <= 1 and 0 < component.mass_fraction <= 1):
            raise DomainError('Component fractions must be in (0, 1].', field='fractions')
        object.__setattr__(component, 'species', species)
        return component


@dataclass(frozen=True, slots=True)
class Mixture:
    """An immutable ideal-gas mixture defined by resolved components [1]_.

    Attributes
    ----------
    components : tuple of MixtureComponent
        Owned components, one per canonical species key (at most 64).
    composition_basis : str, default "mole"
        Basis of ``input_fractions``: ``"mole"`` or ``"mass"``.
    input_fractions : tuple of (str, float)
        Accepted raw fractions by species key, in input order; aliases of one
        species may repeat a key and are combined exactly. When omitted, the
        declared basis' component displays are used after checking that both
        display sets sum to one within 1e-12 and agree with each other.
    normalization_policy : str
        Identifier of the single exact normalization applied to the inputs.

    References
    ----------
    .. [1] NASA RP-1311, "Computer Program for Calculation of Complex Chemical
           Equilibrium Compositions and Applications, Part I: Analysis".
           https://ntrs.nasa.gov/citations/19950013764
    """

    components: tuple[MixtureComponent, ...]
    composition_basis: str = 'mole'
    input_fractions: tuple[tuple[str, float], ...] = ()
    normalization_policy: str = NORMALIZATION_POLICY

    def __post_init__(self) -> None:
        basis = _basis(self.composition_basis)
        policy = text(self.normalization_policy, 'normalization_policy')
        if policy != NORMALIZATION_POLICY:
            raise DomainError('Unsupported mixture normalization policy.', field='normalization_policy')
        source = _bounded(self.components, 'components')
        if any(not isinstance(c, MixtureComponent) for c in source):
            raise DomainError('Mixture needs MixtureComponent records.', field='components')
        given = tuple(MixtureComponent(c.species, c.mole_fraction, c.mass_fraction) for c in source)
        keys = tuple(_canonical(c.species) for c in given)
        if len(set(keys)) != len(keys):
            raise DomainError('Direct mixtures must have one component per canonical species key.',
                              code='duplicate_species', field='components')
        species = tuple(c.species for c in given)
        if type(self.input_fractions) is tuple and not self.input_fractions:
            raw = self._display_inputs(given, basis)
        else:
            raw = []
            for item in _bounded(self.input_fractions, 'input_fractions'):
                key, value = _ordered_items(item, 'input_fractions', size=2)
                raw.append((text(key, 'species').strip(), real(value, 'fraction')))
            raw = tuple(raw)
        resolved = {_canonical(item): item for item in species}
        try:
            order, values = _coalesced(tuple(resolved[key.upper()] for key, _ in raw), [v for _, v in raw])
        except KeyError as error:
            raise DomainError('Raw fractions must identify the owned components.', field='input_fractions') from error
        if tuple(_canonical(item) for item in order) != keys:
            raise DomainError('Raw fractions must identify the owned components in order.', field='input_fractions')
        mole, mass, _ = _composition(species, values, basis)
        # The species in ``given`` are this mixture's own copies; share them.
        components = tuple(MixtureComponent._sharing(item, _rational(x, 'mole_fraction'), _rational(y, 'mass_fraction'))
                           for item, x, y in zip(species, mole, mass, strict=True))
        if self.input_fractions and any(
                (c.mole_fraction, c.mass_fraction) != (o.mole_fraction, o.mass_fraction)
                for c, o in zip(given, components, strict=True)):
            raise DomainError('Component displays contradict retained composition inputs.', field='components')
        if max(s.minimum_temperature_k for s in species) > min(s.maximum_temperature_k for s in species):
            raise DomainError('Mixture species have no common temperature range.', field='components')
        object.__setattr__(self, 'components', components)
        object.__setattr__(self, 'composition_basis', basis)
        object.__setattr__(self, 'input_fractions', raw)
        object.__setattr__(self, 'normalization_policy', policy)

    @staticmethod
    def _display_inputs(components, basis):
        """Check direct component displays and return the declared basis as raw inputs."""
        mole_total = fsum(c.mole_fraction for c in components)
        mass_total = fsum(c.mass_fraction for c in components)
        for total in (mole_total, mass_total):
            if not isclose(total, 1.0, rel_tol=0, abs_tol=1e-12):
                raise DomainError('Direct component fractions must sum to one.', code='composition_sum',
                                  field='components')
        species = tuple(c.species for c in components)
        declared = tuple(c.mole_fraction if basis == 'mole' else c.mass_fraction for c in components)
        mole, mass, _ = _composition(species, declared, basis)
        other = mass if basis == 'mole' else mole
        for component, expected in zip(components, other, strict=True):
            shown = component.mass_fraction / mass_total if basis == 'mole' else component.mole_fraction / mole_total
            _property_identity(shown, float(expected), 'component fraction')
        return tuple((c.species.key, value) for c, value in zip(components, declared, strict=True))

    @classmethod
    @numeric_range_guard('fractions')
    def from_fractions(
        cls, fractions: Mapping[str, float], *, basis: str = "mole"
    ) -> Mixture:
        """Build a mixture from fractions that must sum to one.

        Parameters
        ----------
        fractions : mapping of str to float
            Positive fraction per species key (at most 64 entries). Aliases
            of one species are combined exactly before any mixing calculation.
            A sum within the stated 1e-9 tolerance is normalized exactly once
            to remove input rounding only; compositions outside that tolerance
            are rejected. The accepted raw values are retained.
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
            If a fraction is not positive, a nonzero fraction cannot be
            represented, or the fractions do not sum to one within 1e-9 (code
            ``composition_sum``).
        """

        normalized_basis = _basis(basis, 'basis')
        if not isinstance(fractions, Mapping):
            raise InputError("Composition must contain at least one species.", field="fractions")
        try:
            items = tuple(islice(fractions.items(), MAXIMUM_COMPONENTS + 1))
        except (TypeError, ValueError, KeyError, IndexError, OverflowError) as error:
            raise InputError('Invalid composition mapping.', field='fractions') from error
        if not items:
            raise InputError("Composition must contain at least one species.", field="fractions")
        if len(items) > MAXIMUM_COMPONENTS:
            raise DomainError('Mixture requires 1..64 species.', code='mixture_component_resource',
                              field='fractions')
        raw = []
        for item in items:
            key, value = _ordered_items(item, 'fractions', size=2)
            species = get_species(key)
            fraction = real(value, 'fractions')
            if fraction <= 0.0:
                raise DomainError(f"Fraction for {species.key} must be greater than zero.", field="fractions")
            raw.append((species, fraction))
        order, values = _coalesced(tuple(s for s, _ in raw), [v for _, v in raw])
        mole, mass, _ = _composition(order, values, normalized_basis)
        # Catalogue records are shared here; the constructor takes the one owned copy.
        components = tuple(MixtureComponent._sharing(s, _rational(x, 'mole_fraction'), _rational(y, 'mass_fraction'))
                           for s, x, y in zip(order, mole, mass, strict=True))
        return cls(components, normalized_basis, tuple((s.key, v) for s, v in raw))

    @property
    def molecular_mass_kg_mol(self) -> float:
        """Mole-fraction-weighted molar mass from the exact composition, kg/mol."""
        return _rational(_ledger(self)[1][2], 'molecular_mass_kg_mol')

    @property
    def minimum_temperature_k(self) -> float:
        """Lowest temperature valid for every component, K."""
        return max(component.species.minimum_temperature_k for component in self.components)

    @property
    def maximum_temperature_k(self) -> float:
        """Highest temperature valid for every component, K."""
        return min(component.species.maximum_temperature_k for component in self.components)

    @numeric_range_guard('mixture_state')
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
            Frozen-composition properties at the state, owning this mixture.

        Raises
        ------
        DomainError
            If the pressure is not positive, a component is evaluated outside
            its range, any component ``cv`` is not positive, or a mandatory
            nonzero output is not representable.

        Notes
        -----
        Molar ``cp``, ``H`` and standard-state ``S`` are mole-fraction
        averages of the species values [2]_. The ideal-mixture entropy adds
        ``-R sum(x_j ln x_j) - R ln(p/p_ref)`` with ``p_ref`` the 1 bar
        standard state [1]_. Sums are exact; logarithms are reduced to
        pairwise-coprime integer bases and enclosed before one rounding.

        References
        ----------
        .. [1] NASA RP-1311, "Computer Program for Calculation of Complex Chemical
               Equilibrium Compositions and Applications, Part I: Analysis".
               https://ntrs.nasa.gov/citations/19950013764
        .. [2] NASA/TP-2002-211556. NASA Glenn, "NASA Glenn Coefficients for
               Calculating Thermodynamic Properties of Individual Species".
               https://ntrs.nasa.gov/citations/20020085330
        """

        t = _positive(temperature_k, 'temperature_k')
        p = _positive(pressure_pa, 'pressure_pa')
        return MixtureProperties(temperature_k=t, pressure_pa=p, mixture=self, **_forward(self, t, p))

    @numeric_range_guard('mixture_state')
    def _binary64_state(self, temperature_k: float, pressure_pa: float = REFERENCE_PRESSURE_PA) -> _Binary64State:
        """Evaluate the frozen state in binary64 for iterative solver internals.

        This is main's original mixing algorithm over the component displays
        and binary64 species states. Exported records use :meth:`properties`.
        """
        return _binary64_mixing(tuple((c.species, c.mole_fraction) for c in self.components),
                                _positive(temperature_k, 'temperature_k'), _positive(pressure_pa, 'pressure_pa'))

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
            composition, using :meth:`properties` [1]_. No inverse enclosure
            is asserted.

        Raises
        ------
        DomainError
            If the target or pressure is malformed (before any fit is
            evaluated) or the target is not bracketed by the common range.

        References
        ----------
        .. [1] NASA RP-1311, "Computer Program for Calculation of Complex Chemical
               Equilibrium Compositions and Applications, Part I: Analysis".
               https://ntrs.nasa.gov/citations/19950013764
        """

        target = real(enthalpy_j_kg, 'enthalpy_j_kg')
        pressure_pa = _positive(pressure_pa, 'pressure_pa')
        return solve_bisect_monotonic(
            lambda temperature: self.properties(
                temperature, pressure_pa=pressure_pa
            ).enthalpy_j_kg,
            target,
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
            composition, using :meth:`properties` [1]_. No inverse enclosure
            is asserted.

        Raises
        ------
        DomainError
            If the target or pressure is malformed (before any fit is
            evaluated) or the target is not bracketed by the common range.

        References
        ----------
        .. [1] NASA RP-1311, "Computer Program for Calculation of Complex Chemical
               Equilibrium Compositions and Applications, Part I: Analysis".
               https://ntrs.nasa.gov/citations/19950013764
        """

        target = real(entropy_j_kg_k, 'entropy_j_kg_k')
        pressure_pa = _positive(pressure_pa, 'pressure_pa')
        return solve_bisect_monotonic(
            lambda temperature: self.properties(
                temperature, pressure_pa=pressure_pa
            ).entropy_j_kg_k,
            target,
            self.minimum_temperature_k,
            self.maximum_temperature_k,
            tolerance=1.0e-10,
        ).value


def _binary64_mixing(entries, temperature_k, pressure_pa):
    """Main's original binary64 ideal-mixing evaluation over (species, mole fraction) pairs."""
    states = tuple((x, species._binary64_properties(temperature_k)) for species, x in entries)
    cp_molar = fsum(x * s.cp_molar_j_mol_k for x, s in states)
    enthalpy_molar = fsum(x * s.enthalpy_molar_j_mol for x, s in states)
    standard_entropy_molar = fsum(x * s.entropy_molar_j_mol_k for x, s in states)
    mixing_entropy = -UNIVERSAL_GAS_CONSTANT_J_MOL_K * fsum(x * log(x) for _, x in entries)
    # Preserve near-reference changes with log1p, and do not underflow
    # p/p_ref to zero for representable, positive subnormal pressures.
    log_pressure_ratio = (log1p((pressure_pa - REFERENCE_PRESSURE_PA) / REFERENCE_PRESSURE_PA)
                          if 0.5 * REFERENCE_PRESSURE_PA <= pressure_pa <= 1.5 * REFERENCE_PRESSURE_PA
                          else log(pressure_pa) - log(REFERENCE_PRESSURE_PA))
    entropy_molar = standard_entropy_molar + mixing_entropy - UNIVERSAL_GAS_CONSTANT_J_MOL_K * log_pressure_ratio
    molecular_mass = fsum(x * species.molecular_mass_kg_mol for species, x in entries)
    gas_constant = UNIVERSAL_GAS_CONSTANT_J_MOL_K / molecular_mass
    cp_mass = cp_molar / molecular_mass
    cv_mass = cp_mass - gas_constant
    if cv_mass <= 0.0:
        raise DomainError("Mixture produced non-positive cv.", code="invalid_thermodynamic_property",
                          details={"temperature_k": temperature_k, "cv_j_kg_k": cv_mass})
    return _Binary64State(molecular_mass, gas_constant, cp_mass, cv_mass, cp_mass / cv_mass,
                          enthalpy_molar / molecular_mass, entropy_molar / molecular_mass)


def _binary64_frozen(mole_fractions, temperature_k, pressure_pa):
    """Return (mass fractions, binary64 state) of a frozen mole composition for solvers.

    Reproduces main's original binary64 composition steps (alias combination,
    1e-9 sum admission and normalization, the direct constructor's second
    normalization) and mixing, without building owned records. Exported
    mixtures use :meth:`Mixture.from_fractions` and :meth:`Mixture.properties`.
    """
    amounts: dict[str, list[float]] = {}
    species_by_key: dict[str, Species] = {}
    for key, value in mole_fractions.items():
        species = get_species(key)
        fraction = real(value, 'fractions')
        if fraction <= 0.0:
            raise DomainError(f"Fraction for {species.key} must be greater than zero.", field="fractions")
        species_by_key[_canonical(species)] = species
        amounts.setdefault(_canonical(species), []).append(fraction)
    coalesced = tuple((species_by_key[key], fsum(values)) for key, values in amounts.items())
    total = fsum(fraction for _, fraction in coalesced)
    if not isclose(total, 1.0, rel_tol=0.0, abs_tol=1.0e-9):
        raise DomainError("Composition fractions must sum to one.", code="composition_sum",
                          field="fractions", details={"sum": total})
    entries = tuple((species, fraction / total) for species, fraction in coalesced)
    mole_total = fsum(fraction for _, fraction in entries)
    entries = tuple((species, fraction / mole_total) for species, fraction in entries)
    molar_mass = fsum(fraction * species.molecular_mass_kg_mol for species, fraction in entries)
    mass_fractions = {species.key: fraction * species.molecular_mass_kg_mol / molar_mass
                      for species, fraction in entries}
    return mass_fractions, _binary64_mixing(entries, real(temperature_k, 'temperature_k'),
                                            real(pressure_pa, 'pressure_pa'))


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
    mixture : Mixture
        Owned copy of the defining mixture. All seven derived fields are
        replayed from it; it does not take part in record equality.
    physical_model_qualified : bool, default False
        Always ``False``: these records do not qualify physical applicability.

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
    mixture: Mixture = field(compare=False)
    physical_model_qualified: bool = False

    def __post_init__(self) -> None:
        coerce_fields(self, 'Mixture properties')
        if type(self.physical_model_qualified) is not bool or self.physical_model_qualified:
            raise DomainError('Mixture properties do not establish physical qualification.',
                              field='physical_model_qualified')
        if self.temperature_k <= 0 or self.pressure_pa <= 0 or self.molecular_mass_kg_mol <= 0:
            raise DomainError('Mixture temperature, pressure and molar mass must be positive.')
        _caloric_identity(self.cp_j_kg_k, self.cv_j_kg_k, self.gas_constant_j_kg_k, self.gamma)
        _property_identity(self.gas_constant_j_kg_k * self.molecular_mass_kg_mol,
                           UNIVERSAL_GAS_CONSTANT_J_MOL_K, 'R')
        if not isinstance(self.mixture, Mixture):
            raise DomainError('Properties require an owned Mixture.', field='mixture')
        owned = Mixture(self.mixture.components, self.mixture.composition_basis,
                        self.mixture.input_fractions, self.mixture.normalization_policy)
        expected = _forward(owned, self.temperature_k, self.pressure_pa)
        if any(getattr(self, name) != value for name, value in expected.items()):
            raise DomainError('Mixture properties contradict their owned inputs.', field='properties')
        object.__setattr__(self, 'mixture', owned)
