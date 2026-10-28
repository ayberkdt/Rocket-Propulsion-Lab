"""Primary NASA records, dimensional scale and coupled thermochemical identities.

The raw reference records come from pinned NASA input, not the bundled JSON
or its production parser. Appendix A of NASA/TP-2002-211556 defines R; its
Appendix D prints the H digits. This is model/data verification, not hardware
validation or a uniform polynomial-error bound.
"""

import hashlib
import json
import os
import re
import unittest
from decimal import Decimal
from pathlib import Path

import mpmath as mp

from rocket_propulsion.core import constants
from rocket_propulsion.thermochemistry import (
    EquilibriumProblem,
    LocalEquilibriumProvider,
    Mixture,
    get_reactant,
    get_species,
)
from rocket_propulsion.thermochemistry.cea import render_cea_equilibrium_input
from rocket_propulsion.thermochemistry.equilibrium_derivatives import derivatives_from_composition
from rocket_propulsion.thermochemistry.kinetics import (
    ArrheniusRate,
    Mechanism,
    Reaction,
    _compile,
    concentrations_from_mole_fractions,
    equilibrium_constants,
)
from rocket_propulsion.thermochemistry.nasa_polynomials import UNIVERSAL_GAS_CONSTANT_J_MOL_K
from rocket_propulsion.thermochemistry.nozzle_kinetics import _Expansion
from rocket_propulsion.thermochemistry.species import DATA_FILE
from rocket_propulsion.thermochemistry.transport import modified_eucken_conductivity

REFERENCE = json.loads((Path(__file__).parent / 'reference_data/nasa_glenn_primary_records.json').read_text())
SOURCE_R = Decimal('8.314510')
PINNED_SHA = 'fa7746572952d74e249e818a82a35c113829742fb421a308e167185528884363'


def number(token):
    return Decimal(token.strip().replace('D', 'E'))


def source_intervals(key):
    lines = REFERENCE['records'][key]['raw_record'].splitlines()
    result = []
    for start in (2, 5):
        bounds, first, second = lines[start:start + 3]
        # Seven cp terms, blank field, and two integration constants.
        coefficients = [number(first[i:i + 16]) for i in range(0, 80, 16)]
        coefficients += [number(second[i:i + 16]) for i in (0, 16, 48, 64)]
        result.append((number(bounds[:11]), number(bounds[11:21]), coefficients))
    return result


def source_properties(key, temperature, interval=None):
    """90-digit evaluation from printed decimal coefficients, with term scales."""
    t = mp.mpf(temperature)
    intervals = source_intervals(key)
    if interval is None:
        interval = int(t >= float(intervals[0][1]))
    a, b, c, d, e, f, g, h, s = [mp.mpf(str(v)) for v in intervals[interval][2]]
    cp = [a / t**2, b / t, c, d*t, e*t**2, f*t**3, g*t**4]
    enthalpy = [-a/t, b*mp.log(t), c*t, d*t**2/2, e*t**3/3, f*t**4/4, g*t**5/5, h]
    entropy = [-a/(2*t**2), -b/t, c*mp.log(t), d*t, e*t**2/2, f*t**3/3, g*t**4/4, s]
    r = mp.mpf(str(SOURCE_R))
    groups = (cp, enthalpy, entropy)
    return tuple(r*mp.fsum(v) for v in groups), tuple(r*mp.fsum(abs(x) for x in v) for v in groups)


class PrimaryRecords(unittest.TestCase):
    def test_raw_reference_provenance_and_available_full_source(self):
        self.assertEqual(REFERENCE['source']['thermo_inp_sha256'], PINNED_SHA)
        self.assertEqual(REFERENCE['gas_constant_j_mol_k'], str(SOURCE_R))
        self.assertEqual(len(REFERENCE['records']), 10)
        raw_path = os.environ.get('RPLAB_THERMO_INP')
        raw = None if raw_path is None else Path(raw_path).read_bytes()
        if raw is not None:
            self.assertEqual(hashlib.sha256(raw).hexdigest(), PINNED_SHA)
        for record in (*REFERENCE['records'].values(), *REFERENCE['assigned_reactants'].values()):
            text = record['raw_record']
            self.assertEqual(hashlib.sha256(text.encode()).hexdigest(), record['raw_record_sha256'])
            if raw is not None:
                lines = raw.decode().splitlines()
                start = record['first_line_in_pinned_input'] - 1
                self.assertEqual(lines[start:start + len(text.splitlines())], text.splitlines())

    def test_every_bundled_literal_matches_the_primary_record(self):
        doc = json.loads(DATA_FILE.read_text(), parse_float=Decimal, parse_int=Decimal)
        self.assertEqual({r['key'] for r in doc['species']}, set(REFERENCE['records']))
        for record in doc['species']:
            expected = source_intervals(record['key'])
            self.assertEqual(len(record['intervals']), 2)
            for i, (fit, (low, high, coefficients)) in enumerate(zip(record['intervals'], expected, strict=True)):
                self.assertEqual((fit['minimum_temperature_k'], fit['maximum_temperature_k']), (low, high))
                for j, (actual, source) in enumerate(zip(fit['coefficients'], coefficients, strict=True)):
                    with self.subTest(species=record['key'], interval=i, coefficient=j):
                        self.assertEqual(actual, source)

    def test_molar_masses_and_assigned_enthalpies_keep_source_units(self):
        doc = json.loads(DATA_FILE.read_text(), parse_float=Decimal, parse_int=Decimal)
        for record in doc['species']:
            header = REFERENCE['records'][record['key']]['raw_record'].splitlines()[1]
            with self.subTest(species=record['key']):
                self.assertEqual(record['molecular_mass_g_mol'], number(header[52:65]))
                self.assertEqual(record['heat_of_formation_j_mol_298'], number(header[65:80]))
                formula = {header[i:i + 2].strip().capitalize(): number(header[i + 2:i + 8])
                           for i in range(10, 50, 8) if header[i:i + 2].strip()}
                self.assertEqual(record['formula'], formula)
        for key, reference in REFERENCE['assigned_reactants'].items():
            lines = reference['raw_record'].splitlines()
            reactant = get_reactant(key)
            self.assertEqual(reactant.temperature_k, float(number(lines[2][:11])))
            self.assertEqual(reactant.enthalpy_j_mol, float(number(lines[1][65:80])))

    def test_nasa_model_constant_is_distinct_from_the_general_si_constant(self):
        self.assertEqual(UNIVERSAL_GAS_CONSTANT_J_MOL_K, float(SOURCE_R))
        self.assertEqual(constants.MOLAR_GAS_CONSTANT_J_MOL_K, 8.31446261815324)
        self.assertNotEqual(UNIVERSAL_GAS_CONSTANT_J_MOL_K, constants.MOLAR_GAS_CONSTANT_J_MOL_K)

    def test_dimensional_properties_match_high_precision_primary_evaluation(self):
        with mp.workdps(90):
            for key in REFERENCE['records']:
                for index, (low, high, _) in enumerate(source_intervals(key)):
                    fit = get_species(key).intervals[index]
                    for t in (float(low), float((low + high)/2), float(high)):
                        reference, scales = source_properties(key, t, index)
                        actual = (fit.cp_molar_j_mol_k(t), fit.enthalpy_molar_j_mol(t), fit.entropy_molar_j_mol_k(t))
                        for field, value, expected, scale in zip(('cp', 'h', 's'), actual, reference, scales, strict=True):
                            # A sample-wise roundoff allowance proportional to
                            # the absolute term sum, not a near-zero result.
                            bound = 64 * mp.mpf(2)**-53 * scale
                            with self.subTest(species=key, interval=index, temperature=t, field=field):
                                self.assertLessEqual(abs(mp.mpf(value)-expected), bound)


class CoupledNasaBasis(unittest.TestCase):
    def test_mixing_pressure_entropy_and_caloric_gap_share_the_nasa_basis(self):
        gas = Mixture.from_fractions({'H2': 0.375, 'N2': 0.625})
        state = gas.properties(1400.0, pressure_pa=230000.0)
        with mp.workdps(90):
            r = mp.mpf(str(SOURCE_R));x = (mp.mpf('0.375'), mp.mpf('0.625'))
            standard = mp.fsum(a*source_properties(k, 1400.0)[0][2] for a, k in zip(x, ('H2', 'N2'), strict=True))
            expected = standard-r*mp.fsum(a*mp.log(a) for a in x)-r*mp.log(mp.mpf('2.3'))
            self.assertAlmostEqual(state.entropy_j_kg_k*state.molecular_mass_kg_mol,
                                   float(expected), delta=2e-10)
        self.assertAlmostEqual((state.cp_j_kg_k-state.cv_j_kg_k)*state.molecular_mass_kg_mol,
                               float(SOURCE_R), delta=2e-13)

    def test_pure_gas_uv_recovers_a_source_defined_energy_and_pressure(self):
        with mp.workdps(90):
            target = float(source_properties('N2', 2800.0)[0][1]-mp.mpf(str(SOURCE_R))*2800)
        result = LocalEquilibriumProvider().solve(EquilibriumProblem.uv(
            {'N2': 1.0}, volume_m3=1.0, target_internal_energy_j=target, candidate_species=('N2',)))
        # The provider's UV search is currently restricted to 1800..6000 K.
        # Isolate the dimensional basis here; broader reachability is open.
        self.assertAlmostEqual(result.temperature_k, 2800.0, delta=5e-5)
        self.assertAlmostEqual(result.pressure_pa, float(SOURCE_R)*2800.0, delta=6e-4)

    def test_assigned_liquid_hp_balance_closes_in_primary_dimensional_units(self):
        target = sum(weight*float(number(REFERENCE['assigned_reactants'][key]['raw_record'].splitlines()[1][65:80]))
                     for key, weight in (('H2(L)', 1.0), ('O2(L)', 0.5)))
        result = LocalEquilibriumProvider().solve(EquilibriumProblem.hp(
            {'H2(L)': 1.0, 'O2(L)': 0.5}, pressure_pa=5e6,
            candidate_species=('H2', 'O2', 'H2O', 'OH', 'H', 'O')))
        with mp.workdps(90):
            h = mp.fsum(mp.mpf(item.moles)*source_properties(item.species, result.temperature_k)[0][1]
                        for item in result.species)
            self.assertAlmostEqual(float(h), target, delta=2e-4)

    def test_pure_species_derivatives_use_the_same_caloric_and_eos_basis(self):
        d = derivatives_from_composition({'N2': 1.0}, 1300.0, 2e5)
        with mp.workdps(90):
            cp = float(source_properties('N2', 1300.0)[0][0])
        mass = get_species('N2').molecular_mass_kg_mol
        self.assertAlmostEqual(d.frozen_cp_j_kg_k*mass, cp, delta=2e-10)
        self.assertAlmostEqual((d.frozen_cp_j_kg_k-d.frozen_cv_j_kg_k)*mass, float(SOURCE_R), delta=2e-12)

    def test_reaction_equilibrium_and_concentration_use_the_nasa_basis(self):
        mechanism = Mechanism(name='primary hydrogen', species=('H2', 'H'), reactions=(Reaction(
            equation='H2 <=> 2 H', reactants=(('H2', 1),), products=(('H', 2),),
            rate=ArrheniusRate(1.0, 0.0, 0.0)),))
        concentrations = concentrations_from_mole_fractions(mechanism, 2800.0, 3e5, {'H2': 1.0})
        self.assertAlmostEqual(sum(concentrations)*float(SOURCE_R)*2800.0, 3e5, delta=1e-8)
        with mp.workdps(90):
            h2, _ = source_properties('H2', 2800.0);h, _ = source_properties('H', 2800.0)
            dg = 2*(h[1]-2800*h[2])-(h2[1]-2800*h2[2])
            expected = mp.exp(-dg/(mp.mpf(str(SOURCE_R))*2800))*1e5/(mp.mpf(str(SOURCE_R))*2800)
            self.assertAlmostEqual(equilibrium_constants(mechanism, 2800.0)[0]/float(expected), 1.0, delta=2e-12)

    def test_arrhenius_activation_energy_retains_its_si_rate_convention(self):
        with mp.workdps(90):
            si = mp.mpf('6.02214076e23')*mp.mpf('1.380649e-23')
            expected = 1000*mp.exp(-mp.mpf(12345)/(si*1500))
            self.assertAlmostEqual(ArrheniusRate(1000.0, 0.0, 12345.0)(1500.0)/float(expected), 1.0, delta=3e-15)

    def test_finite_rate_fixed_composition_kernel_keeps_source_energy_entropy_and_eos(self):
        mechanism = Mechanism(name='frozen hydrogen kernel', species=('H2', 'H'), reactions=(Reaction(
            equation='H2 <=> 2 H', reactants=(('H2', 1),), products=(('H', 2),),
            rate=ArrheniusRate(1.0, 0.0, 0.0)),))
        header = REFERENCE['records']['H2']['raw_record'].splitlines()[1]
        mass = float(number(header[52:65]))/1000
        # Isolate the dimensional kernel; this fixture does not exercise
        # chamber HP solution or the nozzle's coupled ODE integration.
        expansion = _Expansion.__new__(_Expansion)
        expansion.compiled = _compile(mechanism)
        with mp.workdps(90):
            r = mp.mpf(str(SOURCE_R));cp_h_s, _ = source_properties('H2', 2000.0)
            h0 = source_properties('H2', 2800.0)[0][1]
            expansion.h0 = float(h0/mass)
            entropy = float((cp_h_s[2]-r*mp.log(2))/mass)
            expected = (cp_h_s[0]/mass, mp.log(2), 2*(h0-cp_h_s[1])/mass,
                        cp_h_s[0]/(cp_h_s[0]-r), r/mass)
            actual = expansion._evaluate(2000.0, (1/mass, 0.0), 1/mass, 0.0, entropy)
            for field, value, reference in zip(('cp', 'log pressure ratio', 'velocity squared', 'gamma', 'R/M'),
                                               actual[:5], expected, strict=True):
                with self.subTest(field=field):
                    self.assertAlmostEqual(value/float(reference), 1.0, delta=3e-12)

    def test_cea_eucken_estimate_keeps_cea_dimensional_scale(self):
        r = float(SOURCE_R);expected = 2e-5*(r/0.04)*3.75
        self.assertAlmostEqual(modified_eucken_conductivity(2e-5, 2.5*r, 0.04)/expected, 1.0, delta=3e-15)

    def test_external_cea_enthalpy_input_uses_cea_r_units(self):
        problem = EquilibriumProblem.hp({'N2': 1.0}, pressure_pa=2e5, target_enthalpy_j=12345.0)
        deck = render_cea_equilibrium_input(problem)
        value = float(re.search(r'h/r=([^ ]+)', deck.splitlines()[0]).group(1))
        expected = 12345.0/(get_species('N2').molecular_mass_kg_mol*float(SOURCE_R)*1000)
        self.assertAlmostEqual(value/expected, 1.0, delta=3e-12)
