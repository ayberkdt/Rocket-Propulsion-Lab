"""Primary-data, arithmetic and catalogue corruption attacks; no hardware claim."""

import copy
import json
import math
import tempfile
import unittest
from dataclasses import replace
from decimal import (
    ROUND_UP,
    DecimalException,
    DefaultContext,
    FloatOperation,
    Inexact,
    localcontext,
)
from fractions import Fraction
from pathlib import Path
from unittest.mock import patch

import mpmath as mp

from rocket_propulsion.core.errors import DomainError
from rocket_propulsion.thermochemistry import NasaPolynomialInterval, get_species, species_catalog
from rocket_propulsion.thermochemistry import species as models

ROOT = Path(__file__).resolve().parents[1]
REFERENCE = json.loads((ROOT / 'docs/assurance/nasa9_catalogue_reference.json').read_text())
SOURCE_R = 8.314510  # NASA/TP-2002-211556 Appendix A, printed p.7/PDF p.15.


def independent(interval, temperature):
    """180-digit source equations, with binary input values preserved exactly."""
    with mp.workdps(180):
        t = mp.mpf(temperature)
        a, b, c, d, e, f, g, h, s = map(mp.mpf, interval.coefficients)
        r = mp.mpf(SOURCE_R)
        return tuple(map(float, (
            r * (a/t**2 + b/t + c + d*t + e*t**2 + f*t**3 + g*t**4),
            r * (-a/t + b*mp.log(t) + c*t + d*t**2/2 + e*t**3/3 + f*t**4/4 + g*t**5/5 + h),
            r * (-a/(2*t**2) - b/t + c*mp.log(t) + d*t + e*t**2/2 + f*t**3/3 + g*t**4/4 + s),
        )))


class Nasa9Attacks(unittest.TestCase):
    def test_all_180_coefficients_and_ten_species_metadata_match_primary_extraction(self):
        for row in REFERENCE['species']:
            species = get_species(row['key'])
            with self.subTest(species=species.key):
                self.assertEqual(dict(species.formula), row['formula'])
                self.assertEqual(species.molecular_mass_kg_mol, float(Fraction(float(row['molecular_mass_g_mol']))/1000))
                self.assertEqual(species.heat_of_formation_j_mol_298, float(row['heat_of_formation_j_mol_298']))
                self.assertEqual(len(species.intervals), 2)
                for actual, expected in zip(species.intervals, row['intervals'], strict=True):
                    self.assertEqual(actual.minimum_temperature_k, expected['minimum_temperature_k'])
                    self.assertEqual(actual.maximum_temperature_k, expected['maximum_temperature_k'])
                    self.assertEqual(actual.coefficients, tuple(map(float, expected['coefficients'])))

    def test_source_R_and_independent_equations_over_every_bundled_interval(self):
        for species in species_catalog().values():
            for iv in species.intervals:
                for t in (iv.minimum_temperature_k, (iv.minimum_temperature_k+iv.maximum_temperature_k)/2,
                          iv.maximum_temperature_k):
                    with self.subTest(species=species.key, temperature=t, lower=iv.minimum_temperature_k):
                        expected = independent(iv, t)
                        actual = (iv.cp_molar_j_mol_k(t), iv.enthalpy_molar_j_mol(t), iv.entropy_molar_j_mol_k(t))
                        self.assertEqual(actual, expected)

    def test_derivative_identities_use_independent_high_precision_not_inverse_round_trip(self):
        with mp.workdps(110):
            for species in species_catalog().values():
                for iv in species.intervals:
                    t = mp.mpf((iv.minimum_temperature_k+iv.maximum_temperature_k)/2)
                    a,b,c,d,e,f,g,h,s = map(mp.mpf, iv.coefficients)
                    def enthalpy(x, *, a=a,b=b,c=c,d=d,e=e,f=f,g=g,h=h):
                        return mp.mpf(SOURCE_R)*(-a/x+b*mp.log(x)+c*x+d*x*x/2+e*x**3/3+f*x**4/4+g*x**5/5+h)
                    def entropy(x, *, a=a,b=b,c=c,d=d,e=e,f=f,g=g,s=s):
                        return mp.mpf(SOURCE_R)*(-a/(2*x*x)-b/x+c*mp.log(x)+d*x+e*x*x/2+f*x**3/3+g*x**4/4+s)
                    with self.subTest(species=species.key, lower=iv.minimum_temperature_k):
                        self.assertEqual(iv.cp_molar_j_mol_k(float(t)), float(mp.diff(enthalpy,t)))
                        self.assertEqual(iv.cp_molar_j_mol_k(float(t)), float(t*mp.diff(entropy,t)))

    def test_catalogue_dictionary_edits_cannot_change_later_resolution(self):
        catalog = species_catalog()
        original = dict(catalog)
        try:
            catalog['H2'] = catalog['O2']
            catalog.pop('CO')
            catalog['FAKE'] = original['H2']
            self.assertEqual(get_species('H2').key, 'H2')
            self.assertEqual(set(species_catalog()), set(original))
            self.assertEqual(get_species('CO').key, 'CO')
        finally:
            catalog.clear()
            catalog.update(original)

    def test_intervals_snapshot_coefficients_instead_of_borrowing_mutable_lists(self):
        coefficients = [0,0,3.5,0,0,0,0,0,1]
        iv = NasaPolynomialInterval(200,1000,coefficients)
        before = iv.cp_molar_j_mol_k(300)
        coefficients[2] = 100
        self.assertEqual(iv.cp_molar_j_mol_k(300), before)
        self.assertIs(type(iv.coefficients), tuple)
        self.assertTrue(all(type(x) is float for x in iv.coefficients))

    def test_invalid_polynomial_scalars_and_cardinality_fail_at_construction(self):
        for value in (True,False,None,'200',float('nan'),float('inf'),-float('inf')):
            for field in ('minimum_temperature_k','maximum_temperature_k'):
                args = {'minimum_temperature_k':200,'maximum_temperature_k':1000,'coefficients':(0,)*9}
                args[field] = value
                with self.subTest(field=field,value=value), self.assertRaises(DomainError):
                    NasaPolynomialInterval(**args)
            with self.subTest(coefficient=value), self.assertRaises(DomainError):
                NasaPolynomialInterval(200,1000,(0,)*8+(value,))
        for values in (None,'123456789',{0,1,2,3,4,5,6,7,8},(0,)*8,(0,)*10):
            with self.subTest(coefficients=values), self.assertRaises(DomainError):
                NasaPolynomialInterval(200,1000,values)

    def test_fixed_coefficient_cardinality_reads_only_ten_elements(self):
        reads = []
        def values():
            for i in range(11):
                reads.append(i)
                if i == 10:
                    raise AssertionError('Unbounded coefficient read')
                yield 0
        with self.assertRaises(DomainError):
            NasaPolynomialInterval(200,1000,values())
        self.assertEqual(reads,list(range(10)))

    def test_direct_polynomial_calls_refuse_nonpositive_and_non_numeric_temperature(self):
        iv = get_species('H2').intervals[0]
        for t in (0,-1,True,None,'300',float('nan'),float('inf')):
            for method in (iv.cp_molar_j_mol_k,iv.enthalpy_molar_j_mol,iv.entropy_molar_j_mol_k):
                with self.subTest(t=t,method=method.__name__), self.assertRaises(DomainError):
                    method(t)
        with self.assertRaises(DomainError):
            iv.contains(300,include_upper=1)

    def test_piecewise_breakpoint_belongs_to_upper_interval_and_endpoints_are_allowed(self):
        low = NasaPolynomialInterval(200,1000,(0,0,3,0,0,0,0,0,1))
        high = NasaPolynomialInterval(1000,6000,(0,0,4,0,0,0,0,0,1))
        species = replace(get_species('H2'), intervals=(low,high))
        self.assertEqual(species.properties(1000).cp_molar_j_mol_k,float(Fraction(SOURCE_R)*4))
        self.assertEqual(species.properties(math.nextafter(1000,0)).cp_molar_j_mol_k,float(Fraction(SOURCE_R)*3))
        self.assertEqual(species.properties(6000).cp_molar_j_mol_k,float(Fraction(SOURCE_R)*4))

    def test_positive_finite_temperature_powers_cannot_overflow_or_underflow_first(self):
        for t in (2.0**900,2.0**-1000):
            iv = NasaPolynomialInterval(t/2,t*2,(0,0,2.5,0,0,0,0,0,1))
            with self.subTest(t=t):
                self.assertEqual((iv.cp_molar_j_mol_k(t),iv.enthalpy_molar_j_mol(t),iv.entropy_molar_j_mol_k(t)),independent(iv,t))

    def test_cancellation_preserves_a_finite_heat_capacity(self):
        iv = NasaPolynomialInterval(1,3,(2.0**1023,8,-2.0**1021,0,0,0,0,0,1))
        self.assertEqual(iv.cp_molar_j_mol_k(2),float(Fraction(SOURCE_R)*4))
        # Compensated floating sum can retain the previous dyadic witness.
        # Here the exact division residual is lost before summation starts.
        a = 2.0**1023
        c = -float(Fraction(a)/9)
        iv = NasaPolynomialInterval(1,4,(a,8,c,0,0,0,0,0,1))
        expected = float(Fraction(SOURCE_R)*(Fraction(a)/9+Fraction(8,3)+Fraction(c)))
        self.assertNotEqual(expected,float(Fraction(SOURCE_R)*Fraction(8,3)))
        self.assertEqual(iv.cp_molar_j_mol_k(3),expected)

    def test_one_ulp_outside_each_fit_is_a_raw_reference_but_never_a_species_state(self):
        # Main's contract: raw primitives also serve reference checks outside
        # their interval; Species.properties refuses any state outside range.
        iv = NasaPolynomialInterval(200,1000,get_species('H2').intervals[0].coefficients)
        species = replace(get_species('H2'),intervals=(iv,))
        for t in (math.nextafter(200,-math.inf),math.nextafter(1000,math.inf)):
            with self.subTest(t=t):
                self.assertEqual((iv.cp_molar_j_mol_k(t),iv.enthalpy_molar_j_mol(t),iv.entropy_molar_j_mol_k(t)),
                                 independent(iv,t))
                with self.assertRaises(DomainError):
                    species.properties(t)

    def test_nonfinite_interval_bounds_cannot_pass_order_comparisons(self):
        for lower,upper in ((float('nan'),1000),(200,float('nan')),(200,float('inf'))):
            with self.subTest(lower=lower,upper=upper), self.assertRaises(DomainError):
                NasaPolynomialInterval(lower,upper,(0,0,3.5,0,0,0,0,0,1))

    def test_unrepresentable_nonzero_and_overflow_outputs_fail_typed(self):
        for iv,t in ((NasaPolynomialInterval(1,3,(0,0,1e308,0,0,0,0,0,0)),2),
                     (NasaPolynomialInterval(1e200,1e201,(math.ulp(0.0),0,0,0,0,0,0,0,0)),1e200)):
            with self.subTest(t=t), self.assertRaises(DomainError):
                iv.cp_molar_j_mol_k(t)

    def test_log_evaluation_is_independent_of_callers_decimal_context(self):
        iv = get_species('H2O').intervals[1]
        expected = independent(iv,3500)
        with localcontext() as context:
            context.prec = 3
            context.Emax = 3
            context.traps[FloatOperation] = True
            self.assertEqual((iv.cp_molar_j_mol_k(3500),iv.enthalpy_molar_j_mol(3500),iv.entropy_molar_j_mol_k(3500)),expected)

    def test_mutable_default_context_cannot_inject_traps_into_a_valid_fit(self):
        iv = get_species('H2O').intervals[1]
        expected = independent(iv,3500)
        saved = DefaultContext.copy()
        try:
            DefaultContext.prec = 3
            DefaultContext.rounding = ROUND_UP
            DefaultContext.Emax = 3
            DefaultContext.Emin = -3
            DefaultContext.clamp = 1
            DefaultContext.traps[Inexact] = True
            try:
                actual = (iv.cp_molar_j_mol_k(3500),iv.enthalpy_molar_j_mol(3500),iv.entropy_molar_j_mol_k(3500))
            except DecimalException as error:
                self.fail(f'Valid fit inherited an external decimal trap: {error}')
            self.assertEqual(actual,expected)
        finally:
            for name in ('prec','rounding','Emax','Emin','capitals','clamp'):
                setattr(DefaultContext,name,getattr(saved,name))
            DefaultContext.traps = saved.traps.copy()
            DefaultContext.flags = saved.flags.copy()

    def test_species_constructor_rejects_gaps_overlaps_bad_formula_and_mass(self):
        species = get_species('H2')
        for intervals in ((),(species.intervals[1],species.intervals[0]),
                          (species.intervals[0],replace(species.intervals[1],minimum_temperature_k=1001))):
            with self.subTest(intervals=intervals), self.assertRaises(DomainError):
                replace(species,intervals=intervals)
        for formula in ((),(('H',0),),(('H',-1),),(('H',True),),(('H',1),('H',2))):
            with self.subTest(formula=formula), self.assertRaises(DomainError):
                replace(species,formula=formula)
        for mass in (0,-1,True,None,float('nan'),float('inf')):
            with self.subTest(mass=mass), self.assertRaises(DomainError):
                replace(species,molecular_mass_kg_mol=mass)

    def test_species_owns_nested_records_and_composition(self):
        original = get_species('H2')
        formula = [['H',2]]
        intervals = list(original.intervals)
        result = replace(original,formula=formula,intervals=intervals)
        formula[0][1] = 9
        intervals.clear()
        self.assertEqual(result.formula,(('H',2.0),))
        self.assertEqual(result.intervals,original.intervals)
        self.assertIs(type(result.intervals),tuple)

    def test_species_lookup_and_inverse_targets_reject_boolean_or_protocol_errors(self):
        species = get_species('H2')
        for value in (None,True,4,{},b'H2','\ud800'):
            with self.subTest(key=repr(value)), self.assertRaises(DomainError):
                get_species(value)
        for value in (None,True,'300',float('nan'),float('inf')):
            with self.subTest(t=value), self.assertRaises(DomainError):
                species.properties(value)
            for method in (species.temperature_from_enthalpy,species.temperature_from_entropy):
                with self.subTest(target=value,method=method.__name__), self.assertRaises(DomainError):
                    method(value)

    def test_species_cv_uses_retained_heat_capacity_before_display_rounding(self):
        iv = NasaPolynomialInterval(.5,2,(2.0**-54,0,1,0,0,0,0,0,1))
        species = replace(get_species('H2'),molecular_mass_kg_mol=1,intervals=(iv,))
        state = species.properties(1)
        self.assertEqual(state.cp_j_kg_k,state.gas_constant_j_kg_k)
        self.assertEqual(state.cv_j_kg_k,float(Fraction(SOURCE_R)*Fraction(2)**-54))
        self.assertEqual(state.gamma,float((1+Fraction(2)**-54)/(Fraction(2)**-54)))

    def test_species_mandatory_mass_properties_cannot_return_infinity(self):
        species = replace(get_species('H2'),molecular_mass_kg_mol=math.ulp(0.0))
        with self.assertRaises(DomainError):
            species.properties(300)


class CatalogueDocumentAttacks(unittest.TestCase):
    def test_malformed_json_duplicate_keys_and_nonfinite_constants_are_typed(self):
        valid = models.DATA_FILE.read_text()
        duplicate = valid.replace('"schema_version": 1,','"schema_version": 1, "schema_version": 1,',1)
        nonfinite = valid.replace('"schema_version": 1,','"schema_version": 1, "unused": NaN,',1)
        for content in ('{','{"schema_version":1,"schema_version":1}', '{"schema_version": NaN}',duplicate,nonfinite):
            with self.subTest(content=content):
                self.check_document(content)

    def check_document(self, document):
        content = document if isinstance(document,str) else json.dumps(document)
        cache = models._species_catalog
        with tempfile.TemporaryDirectory() as folder:
            path = Path(folder)/'catalog.json'
            path.write_text(content,encoding='utf8')
            cache.cache_clear()
            try:
                with patch.object(models,'DATA_FILE',path), self.assertRaises(DomainError):
                    species_catalog()
            finally:
                cache.cache_clear()

    def test_wrong_schema_pressure_coefficient_order_or_missing_fields_are_refused(self):
        original = json.loads(models.DATA_FILE.read_text())
        cases=[]
        for version in (True,2,'1'):
            doc=copy.deepcopy(original);doc['schema_version']=version;cases.append(doc)
        doc=copy.deepcopy(original);doc['model']['reference_pressure_pa']=101325;cases.append(doc)
        doc=copy.deepcopy(original);doc['model']['coefficient_order'].reverse();cases.append(doc)
        doc=copy.deepcopy(original);doc['model']['gas_constant_j_mol_k']=8.31446261815324;cases.append(doc)
        doc=copy.deepcopy(original);del doc['species'][0]['intervals'][0]['coefficients'];cases.append(doc)
        doc=copy.deepcopy(original);doc['species']=[];cases.append(doc)
        doc=copy.deepcopy(original);doc['species'].append(doc['species'][0]);cases.append(doc)
        doc=copy.deepcopy(original);doc['species'][0]['molecular_mass_g_mol']=True;cases.append(doc)
        for index,doc in enumerate(cases):
            with self.subTest(case=index):
                self.check_document(doc)

