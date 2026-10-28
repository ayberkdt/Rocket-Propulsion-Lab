"""Admitted-record attacks and independent ideal-mixture identities.

These tests exercise public constructors as well as factories. A frozen
dataclass must own its nested input data; a caller's later list/dict edits
must not alter a thermodynamic state or a process-wide catalogue. Entropy
checks use species-identity invariance and ideal mixing/pressure identities,
not another copy of the record-building code.
"""

import json
import unittest
from dataclasses import asdict, replace
from fractions import Fraction
from math import inf, isfinite, nan

import mpmath as mp

from rocket_propulsion.core.errors import DomainError
from rocket_propulsion.thermochemistry import Mixture, get_species, species_catalog
from rocket_propulsion.thermochemistry.condensed_equilibrium import condensed_catalog
from rocket_propulsion.thermochemistry.mixtures import MixtureComponent
from rocket_propulsion.thermochemistry.nasa_polynomials import (
    UNIVERSAL_GAS_CONSTANT_J_MOL_K as GAS_CONSTANT,
)
from rocket_propulsion.thermochemistry.nasa_polynomials import (
    NasaPolynomialInterval,
    interval_from_record,
)
from rocket_propulsion.thermochemistry.species import condensed_reactant_catalog, get_reactant

COEFFICIENTS = (0.0, 0.0, 3.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)


def interval(lower=200.0, upper=1000.0):
    return NasaPolynomialInterval(lower, upper, COEFFICIENTS)


class PolynomialAdmission(unittest.TestCase):
    def test_raw_polynomial_reference_evaluation_and_state_range_are_distinct(self):
        fit = NasaPolynomialInterval(300.0, 1000.0, COEFFICIENTS)
        with mp.workdps(60):
            t, r = mp.mpf(298.15), mp.mpf(GAS_CONSTANT)
            self.assertAlmostEqual(fit.cp_molar_j_mol_k(float(t)), float(mp.mpf('3.5') * r), places=12)
            self.assertAlmostEqual(fit.enthalpy_molar_j_mol(float(t)), float(mp.mpf('3.5') * r * t), places=10)
            self.assertAlmostEqual(fit.entropy_molar_j_mol_k(float(t)), float(mp.mpf('3.5') * r * mp.log(t)), places=11)
        gas = replace(get_species('H2'), intervals=(fit,))
        with self.assertRaises(DomainError):
            gas.properties(298.15)
        self.assertFalse(fit.contains(298.15))
        self.assertTrue(fit.contains(300.0))
        self.assertFalse(fit.contains(1000.0, include_upper=False))

    def test_malformed_temperature_and_membership_flags_are_refused(self):
        fit = interval()
        for bad in (True, '300', None, nan, inf, -1.0, 0.0):
            for name in ('cp_molar_j_mol_k', 'enthalpy_molar_j_mol', 'entropy_molar_j_mol_k'):
                with self.subTest(value=bad, method=name), self.assertRaises(DomainError):
                    getattr(fit, name)(bad)
        for flag in ('false', 0, 1, None):
            with self.subTest(flag=flag), self.assertRaises(DomainError):
                fit.contains(1000.0, include_upper=flag)

    def test_coefficient_list_is_owned(self):
        values = list(COEFFICIENTS)
        record = NasaPolynomialInterval(200.0, 1000.0, values)
        before = record.cp_molar_j_mol_k(300.0)
        values[2] *= 2
        values.clear()
        self.assertEqual(record.cp_molar_j_mol_k(300.0), before)
        self.assertIsInstance(record.coefficients, tuple)

    def test_bounds_and_coefficients_are_finite_real_and_fixed_size(self):
        for field in ('minimum_temperature_k', 'maximum_temperature_k'):
            for bad in (True, '300', None, nan, inf, -inf, 0.0, -1.0):
                p = {'minimum_temperature_k': 200.0, 'maximum_temperature_k': 1000.0, 'coefficients': COEFFICIENTS}
                p[field] = bad
                with self.subTest(field=field, bad=bad), self.assertRaises(DomainError):
                    NasaPolynomialInterval(**p)
        for values in (COEFFICIENTS[:-1], (*COEFFICIENTS, 0.0), '123456789',
                       {str(i): x for i, x in enumerate(COEFFICIENTS)}):
            with self.subTest(kind=type(values).__name__), self.assertRaises(DomainError):
                NasaPolynomialInterval(200.0, 1000.0, values)

    def test_fixed_arity_does_not_exhaust_an_unbounded_input(self):
        consumed = []

        def guarded_stream():
            for i in range(10):
                consumed.append(i)
                yield 0.0
            # A deliberately broken full-tuple conversion fails safely in
            # the test instead of consuming unbounded RAM on the host.
            raise AssertionError('The fixed-size constructor over-consumed its input.')

        with self.assertRaises(DomainError):
            NasaPolynomialInterval(200.0, 1000.0, guarded_stream())
        self.assertEqual(len(consumed), 10)
        for bad in (True, '0', None, nan, inf, 1j):
            values = (*COEFFICIENTS[:2], bad, *COEFFICIENTS[3:])
            with self.subTest(coefficient=bad), self.assertRaises(DomainError):
                NasaPolynomialInterval(200.0, 1000.0, values)

    def test_record_factory_cannot_bypass_admission(self):
        for record in ({}, None,
                       {'minimum_temperature_k': True, 'maximum_temperature_k': 1000.0, 'coefficients': COEFFICIENTS},
                       {'minimum_temperature_k': 200.0, 'maximum_temperature_k': 1000.0,
                        'coefficients': (*COEFFICIENTS[:-1], nan)}):
            with self.subTest(record=record), self.assertRaises(DomainError):
                interval_from_record(record)

    def test_evaluation_cannot_return_nonfinite_values(self):
        overflowing = NasaPolynomialInterval(200.0, 1000.0, (0.0, 0.0, 1e308, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0))
        for method in ('cp_molar_j_mol_k', 'enthalpy_molar_j_mol', 'entropy_molar_j_mol_k'):
            with self.subTest(method=method), self.assertRaises(DomainError):
                getattr(overflowing, method)(300.0)


class OwnedSpeciesAndPhases(unittest.TestCase):
    def test_positive_caloric_state_can_have_equal_separately_rounded_outputs(self):
        # An admitted constant-cp fit with a positive exact cp-cv gap. Each
        # required field is representable, but cp and cv round to the same
        # binary64 number and gamma rounds to 1. Rejecting those display
        # values would unnecessarily narrow the supported numerical domain.
        base = get_species('H2')
        fit = NasaPolynomialInterval(200.0, 1000.0, (0.0, 0.0, 1e30, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0))
        gas = replace(base, intervals=(fit,))
        with mp.workdps(100):
            r = mp.mpf(GAS_CONSTANT) / mp.mpf(base.molecular_mass_kg_mol)
            cp = mp.mpf(1e30) * r
            cv = cp - r
            self.assertGreater(cv, 0)
            self.assertGreater(cp / cv, 1)
            self.assertEqual(float(cp), float(cv))
            self.assertEqual(float(cp / cv), 1.0)
            state = gas.properties(300.0)
            self.assertEqual(state.cp_j_kg_k, state.cv_j_kg_k)
            self.assertEqual(state.gamma, 1.0)
            self.assertLess(abs(mp.mpf(state.cp_j_kg_k) / cp - 1), mp.mpf('1e-15'))
            self.assertLess(abs(mp.mpf(state.cv_j_kg_k) / cv - 1), mp.mpf('1e-15'))

    def test_source_names_are_preserved_while_lookup_is_case_insensitive(self):
        catalog = condensed_catalog()
        self.assertEqual(catalog.condensed['C(GR)'].key, 'C(gr)')
        self.assertEqual(catalog.condensed['AL2O3(A)'].key, 'AL2O3(a)')
        self.assertEqual(catalog.gases['C2H2,ACETYLENE'].key, 'C2H2,acetylene')
        self.assertIs(get_species(' h2 '), get_species('H2'))

    def test_formula_and_intervals_are_owned(self):
        for base in (get_species('H2'), condensed_catalog().condensed['H2O(L)']):
            formula = [[key, amount] for key, amount in base.formula]
            fits = list(base.intervals)
            owned = replace(base, formula=formula, intervals=fits)
            expected_formula = tuple(base.formula)
            expected_fits = tuple(base.intervals)
            formula[0][1] *= 2
            formula.clear()
            fits.clear()
            with self.subTest(species=base.key):
                self.assertEqual(owned.formula, expected_formula)
                self.assertEqual(owned.intervals, expected_fits)

    def test_invalid_species_records_are_refused(self):
        base = get_species('H2')
        for changes in ({'molecular_mass_kg_mol': 0.0}, {'molecular_mass_kg_mol': nan},
                        {'molecular_mass_kg_mol': True}, {'heat_of_formation_j_mol_298': nan},
                        {'key': ''}, {'source': ''}, {'formula': []}, {'formula': [('H', nan)]},
                        {'formula': [('H', -2.0)]}, {'formula': [('H', 1.0), ('H', 1.0)]},
                        {'intervals': []}, {'intervals': [object()]},
                        {'intervals': [interval(200.0, 500.0), interval(501.0, 1000.0)]},
                        {'intervals': [interval(200.0, 600.0), interval(500.0, 1000.0)]},
                        {'intervals': [interval(500.0, 1000.0), interval(200.0, 500.0)]}):
            with self.subTest(changes=changes), self.assertRaises(DomainError):
                replace(base, **changes)

    def test_condensed_records_share_the_admission_boundary(self):
        liquid = get_reactant('O2(L)')
        phase = condensed_catalog().condensed['H2O(L)']
        for base, changes in ((liquid, {'temperature_k': 0.0}), (liquid, {'enthalpy_j_mol': nan}),
                              (liquid, {'molecular_mass_kg_mol': -1.0}),
                              (liquid, {'formula': [('O', inf)]}),
                              (phase, {'intervals': []}), (phase, {'phase_index': True}),
                              (phase, {'molecular_mass_kg_mol': 0.0})):
            with self.subTest(key=base.key, changes=changes), self.assertRaises(DomainError):
                replace(base, **changes)

    def test_result_records_cannot_claim_inconsistent_or_nonfinite_states(self):
        for state in (get_species('H2').properties(300.0), Mixture.from_fractions({'H2': 1.0}).properties(300.0)):
            for changes in ({'cp_j_kg_k': nan}, {'gamma': inf}, {'gamma': 2.5},
                            {'cv_j_kg_k': -1.0}, {'molecular_mass_kg_mol': 0.0},
                            {'gas_constant_j_kg_k': state.gas_constant_j_kg_k * 2}):
                with self.subTest(kind=type(state).__name__, changes=changes), self.assertRaises(DomainError):
                    replace(state, **changes)


class CatalogueIsolation(unittest.TestCase):
    def test_owned_catalogue_snapshots_remain_json_serializable(self):
        catalog = condensed_catalog()
        encoded = json.dumps(asdict(catalog), allow_nan=False)
        decoded = json.loads(encoded)
        self.assertEqual(decoded['condensed']['C(GR)']['key'], 'C(gr)')
        self.assertEqual(len(decoded['gases']), len(catalog.gases))

    def test_catalogue_keys_must_identify_their_record(self):
        catalog = condensed_catalog()
        h2, n2 = get_species('H2'), get_species('N2')
        for gases in ({'H2': n2}, {'H2': h2, 'h2': h2}, {'H2': object()}, {True: h2}):
            with self.subTest(keys=list(gases)), self.assertRaises(DomainError):
                replace(catalog, gases=gases)

    def test_gas_and_liquid_catalogue_edits_do_not_pollute_lookup(self):
        for factory, lookup, key, other in ((species_catalog, get_species, 'H2', 'N2'),
                                           (condensed_reactant_catalog, get_reactant, 'O2(L)', 'H2(L)')):
            public = factory()
            original = dict(public)
            try:
                try:
                    public[key] = public[other]
                except TypeError:
                    pass
                observed = lookup(key)
                next_snapshot = factory()[key]
            finally:
                if isinstance(public, dict):
                    public.clear()
                    public.update(original)
            with self.subTest(factory=factory.__name__):
                self.assertIs(observed, original[key])
                self.assertIs(next_snapshot, original[key])

    def test_extended_catalogue_mappings_are_owned_and_protected(self):
        catalog = condensed_catalog()
        for name in ('gases', 'condensed', 'reactants'):
            mapping = getattr(catalog, name)
            original = dict(mapping)
            key, other = list(original)[:2]
            try:
                try:
                    mapping[key] = mapping[other]
                except TypeError:
                    pass
                observed = getattr(condensed_catalog(), name)[key]
            finally:
                if isinstance(mapping, dict):
                    mapping.clear()
                    mapping.update(original)
            with self.subTest(mapping=name):
                self.assertIs(observed, original[key])
            inputs = dict(original)
            owned = replace(catalog, **{name: inputs})
            inputs.clear()
            self.assertEqual(dict(getattr(owned, name)), original)


class MixtureIdentityAndOwnership(unittest.TestCase):
    def test_distinct_isomers_do_not_lose_real_mixing_entropy(self):
        catalog = condensed_catalog().gases
        first, second = catalog['C3H4,ALLENE'], catalog['C3H4,PROPYNE']
        self.assertEqual(first.formula, second.formula)
        self.assertEqual(first.molecular_mass_kg_mol, second.molecular_mass_kg_mol)
        mixture = Mixture((MixtureComponent(first, 0.5, 0.5), MixtureComponent(second, 0.5, 0.5)))
        state = mixture.properties(1000.0)
        unmixed = 0.5 * (first.properties(1000.0).entropy_molar_j_mol_k
                         + second.properties(1000.0).entropy_molar_j_mol_k)
        with mp.workdps(80):
            expected = mp.mpf(GAS_CONSTANT) * mp.log(2)
            actual = state.entropy_j_kg_k * state.molecular_mass_kg_mol - unmixed
            self.assertAlmostEqual(actual, float(expected), delta=1e-12)

    def test_direct_pure_component_roundoff_is_normalized_before_properties(self):
        species = get_species('H2')
        pure = Mixture.from_fractions({'H2': 1.0}).properties(300.0)
        for mole, mass in ((1 - 5e-13, 1.0), (1.0, 1 - 5e-13)):
            with self.subTest(mole=mole, mass=mass):
                direct = Mixture((MixtureComponent(species, mole, mass),))
                self.assertEqual(direct.properties(300.0), pure)

    def test_fraction_aliases_do_not_hide_negative_amounts_and_reals_remain_supported(self):
        with self.assertRaises(DomainError):
            Mixture.from_fractions({'H2': 2.0, ' h2 ': -1.0})
        with self.assertRaises(DomainError):
            Mixture.from_fractions({'H2': 1e308, 'N2': 1e308})
        rational = Mixture.from_fractions({'H2': Fraction(1, 4), 'N2': Fraction(3, 4)})
        binary = Mixture.from_fractions({'H2': 0.25, 'N2': 0.75})
        self.assertEqual(rational.properties(300.0), binary.properties(300.0))

    def test_aliases_are_one_physical_species_in_both_bases(self):
        for basis in ('mole', 'mass'):
            for temperature in (300.0, 2000.0, 6000.0):
                pure = Mixture.from_fractions({'H2': 1.0}, basis=basis)
                aliases = Mixture.from_fractions({'H2': 0.4, ' h2 ': 0.6}, basis=basis)
                with self.subTest(basis=basis, temperature=temperature):
                    self.assertEqual(len(aliases.components), 1)
                    self.assertEqual(aliases.properties(temperature), pure.properties(temperature))
        reference = Mixture.from_fractions({'H2': 0.5, 'N2': 0.5}, basis='mass')
        aliases = Mixture.from_fractions({'H2': 0.2, ' h2 ': 0.3, 'N2': 0.5}, basis='mass')
        self.assertEqual(aliases.properties(300.0), reference.properties(300.0))

    def test_accepted_sum_roundoff_does_not_make_a_pure_gas_mix_with_itself(self):
        pure = Mixture.from_fractions({'H2': 1.0}).properties(300.0)
        for value in (1 + 5e-10, 1 - 5e-10):
            with self.subTest(value=value):
                self.assertEqual(Mixture.from_fractions({'H2': value}).properties(300.0), pure)

    def test_component_list_is_owned_and_invalid_direct_mixtures_are_refused(self):
        base = Mixture.from_fractions({'H2': 1.0})
        parts = list(base.components)
        owned = Mixture(parts)
        expected = owned.properties(300.0)
        parts.clear()
        self.assertEqual(owned.properties(300.0), expected)
        h2, n2 = get_species('H2'), get_species('N2')
        for parts in ([], [object()],
                      [MixtureComponent(h2, 0.5, 0.5), MixtureComponent(h2, 0.5, 0.5)],
                      [MixtureComponent(h2, 0.5, 0.5), MixtureComponent(n2, 0.5, 0.5)]):
            with self.subTest(parts=parts), self.assertRaises(DomainError):
                Mixture(parts)

    def test_malformed_fraction_and_pressure_values_are_not_coerced(self):
        mixture = Mixture.from_fractions({'H2': 1.0})
        for bad in (True, '1', None, nan, inf, -1.0, 0.0):
            with self.subTest(fraction=bad), self.assertRaises(DomainError):
                Mixture.from_fractions({'H2': bad})
            with self.subTest(pressure=bad), self.assertRaises(DomainError):
                mixture.properties(300.0, pressure_pa=bad)

    def test_mixing_and_pressure_entropy_against_high_precision_identities(self):
        mixture = Mixture.from_fractions({'H2': 0.25, 'N2': 0.75})
        state = mixture.properties(300.0)
        unmixed = sum(c.mole_fraction * c.species.properties(300.0).entropy_molar_j_mol_k
                      for c in mixture.components)
        with mp.workdps(80):
            x = (mp.mpf('0.25'), mp.mpf('0.75'))
            expected = -mp.mpf(GAS_CONSTANT) * sum(v * mp.log(v) for v in x)
            actual = state.entropy_j_kg_k * state.molecular_mass_kg_mol - unmixed
            self.assertAlmostEqual(actual, float(expected), delta=1e-12)
            for pressure in (5e-324, 1e-100, 1e308):
                result = mixture.properties(300.0, pressure_pa=pressure)
                expected_delta = -mp.mpf(state.gas_constant_j_kg_k) * mp.log(mp.mpf(pressure) / 100000)
                self.assertTrue(isfinite(result.entropy_j_kg_k))
                self.assertAlmostEqual(result.entropy_j_kg_k - state.entropy_j_kg_k,
                                       float(expected_delta), delta=1e-8)


if __name__ == '__main__':
    unittest.main()
