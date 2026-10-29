"""Independent ideal-mixture equations and actual public-boundary attacks."""
import itertools
import json
import math
import unittest
from dataclasses import asdict, fields, replace
from decimal import DefaultContext, Inexact, localcontext
from fractions import Fraction as F
from pathlib import Path

import mpmath as mp

from rocket_propulsion.core.errors import DomainError
from rocket_propulsion.thermochemistry import Mixture, MixtureComponent, get_species

ROOT=Path(__file__).resolve().parents[1]
R=8.314510


def oracle(entries, temperature, pressure, *, basis='mole', digits=180):
    """No production property/normalization helpers; printed NASA9 + ideal mixing."""
    with mp.workdps(digits):
        source=[(s,mp.mpf(v)) for s,v in entries]
        raw=[v if basis=='mole' else v/mp.mpf(s.molecular_mass_kg_mol) for s,v in source]
        weights=[v/sum(raw) for v in raw]
        mass=sum(x*mp.mpf(s.molecular_mass_kg_mol) for (s,_),x in zip(source,weights,strict=True))
        cp=h=entropy=mp.mpf(0);t=mp.mpf(temperature);r=mp.mpf(R)
        for (s,_),x in zip(source,weights,strict=True):
            iv=next(iv for index,iv in enumerate(s.intervals) if iv.minimum_temperature_k<=temperature and
                    (temperature<iv.maximum_temperature_k or index==len(s.intervals)-1 and temperature<=iv.maximum_temperature_k))
            a,b,c,d,e,f,g,j,k=map(mp.mpf,iv.coefficients)
            cp+=x*r*(a/t**2+b/t+c+d*t+e*t**2+f*t**3+g*t**4)
            h+=x*r*(-a/t+b*mp.log(t)+c*t+d*t**2/2+e*t**3/3+f*t**4/4+g*t**5/5+j)
            entropy+=x*(r*(-a/(2*t**2)-b/t+c*mp.log(t)+d*t+e*t**2/2+f*t**3/3+g*t**4/4+k)-r*mp.log(x))
        entropy-=r*(mp.log(mp.mpf(pressure))-mp.log(100_000))
        return {'molecular_mass_kg_mol': float(mass),'gas_constant_j_kg_k': float(r/mass),
                    'cp_j_kg_k': float(cp/mass),'cv_j_kg_k': float((cp-r)/mass),'gamma': float(cp/(cp-r)),
                    'enthalpy_j_kg': float(h/mass),'entropy_j_kg_k': float(entropy/mass)}


def custom(key,mass,coefficients,*,lower=.5,upper=6000):
    from rocket_propulsion.thermochemistry import NasaPolynomialInterval
    return replace(get_species('H2'),key=key,molecular_mass_kg_mol=mass,
                   intervals=(NasaPolynomialInterval(lower,upper,coefficients),))


class MixtureAttacks(unittest.TestCase):
    def assert_oracle(self,mixture,entries,t,p,*,basis='mole',digits=180):
        expected=oracle(entries,t,p,basis=basis,digits=digits)
        try: actual=mixture.properties(t,pressure_pa=p)
        except (DomainError,ArithmeticError) as error: self.fail(f'Independently representable state was refused: {error}')
        for name,value in expected.items():
            with self.subTest(field=name,t=t,p=p): self.assertEqual(getattr(actual,name),value)
        return actual

    def test_independent_source_equations_for_mole_and_mass_compositions(self):
        for basis in ('mole','mass'):
            for values in ({'N2':.79,'O2':.21},{'H2O':.7,'CO2':.3},{'H2':.2,'O2':.3,'OH':.5}):
                mixture=Mixture.from_fractions(values,basis=basis);entries=[(get_species(k),v) for k,v in values.items()]
                with mp.workdps(180):
                    amounts=[mp.mpf(v) if basis=='mole' else mp.mpf(v)/mp.mpf(s.molecular_mass_kg_mol) for s,v in entries]
                    mole=[n/sum(amounts) for n in amounts]
                    mass=[x*mp.mpf(s.molecular_mass_kg_mol) for (s,_),x in zip(entries,mole,strict=True)]
                    for c,x,y in zip(mixture.components,mole,mass,strict=True):
                        self.assertEqual(c.mole_fraction,float(x));self.assertEqual(c.mass_fraction,float(y/sum(mass)))
                for t,p in ((200,100_000),(298.15,2e6),(1000,1e-6),(3500,1e8),(6000,100_000)):
                    with self.subTest(basis=basis,values=values):self.assert_oracle(mixture,entries,t,p,basis=basis)

    def test_minimum_positive_and_maximum_finite_pressures_have_finite_entropy(self):
        mixture=Mixture.from_fractions({'H2':1.0})
        for p in (math.ulp(0.0),float.fromhex('0x1.fffffffffffffp+1023')):
            self.assert_oracle(mixture,[(get_species('H2'),1.0)],300,p,digits=700)

    def test_aliases_of_same_species_are_combined_exactly_before_mixing_entropy(self):
        # Main's contract: aliases are one physical species. Their raw entries
        # are retained, and no mixing entropy is invented between them.
        for values in ({'H2':.5,'h2':.5},{'H2':.4,' H2 ':.6}):
            with self.subTest(values=values):
                mixture=Mixture.from_fractions(values)
                self.assertEqual(len(mixture.components),1)
                self.assertEqual(tuple(v for _,v in mixture.input_fractions),tuple(values.values()))
                self.assert_oracle(mixture,[(get_species('H2'),1.0)],300,100_000)

    def test_boolean_nonfinite_and_textual_fractions_are_not_cast(self):
        for value in (True,False,None,'1',float('nan'),float('inf'),-1,0):
            with self.subTest(value=value),self.assertRaises(DomainError):Mixture.from_fractions({'H2':value})
        for value in (None,[],(),1,'H2'):
            with self.subTest(mapping=value),self.assertRaises(DomainError):Mixture.from_fractions(value)

    def test_basis_and_thermodynamic_controls_are_strict(self):
        for basis in (None,True,4,'volume'):
            with self.subTest(basis=basis),self.assertRaises(DomainError):Mixture.from_fractions({'H2':1},basis=basis)
        mixture=Mixture.from_fractions({'H2':1})
        for p in (True,None,'100000',0,-1,float('inf'),float('nan')):
            with self.subTest(p=p),self.assertRaises(DomainError):mixture.properties(300,pressure_pa=p)
        for t in (True,None,'300',0,199,6001,float('inf'),float('nan')):
            with self.subTest(t=t),self.assertRaises(DomainError):mixture.properties(t)

    def test_unqualified_sum_is_refused_and_accepted_sum_is_explicitly_normalized(self):
        with self.assertRaises(DomainError) as error:Mixture.from_fractions({'N2':.8,'O2':.3})
        self.assertEqual(error.exception.code,'composition_sum')
        source={'N2':.8,'O2':.2000000005};mixture=Mixture.from_fractions(source)
        self.assertEqual(dict(mixture.input_fractions),source)
        self.assertEqual(mixture.composition_basis,'mole')
        total=sum(map(F,source.values()))
        for comp in mixture.components:self.assertEqual(comp.mole_fraction,float(F(source[comp.species.key])/total))
        self.assert_oracle(mixture,[(get_species(k),v) for k,v in source.items()],300,100_000)

    def test_fraction_conversion_does_not_drop_representable_mass_basis_trace(self):
        source={'H2':math.ulp(0.0),'O2':1.0};mixture=Mixture.from_fractions(source,basis='mass')
        x=[F(v)/F(get_species(k).molecular_mass_kg_mol) for k,v in source.items()];total=sum(x)
        for comp,n in zip(mixture.components,x,strict=True):
            self.assertEqual(comp.mole_fraction,float(n/total));self.assertGreater(comp.mass_fraction,0)
        self.assert_oracle(mixture,[(get_species(k),v) for k,v in source.items()],300,100_000,basis='mass',digits=700)

    def test_genuinely_unrepresentable_component_fraction_is_not_silently_erased(self):
        with self.assertRaises(DomainError):Mixture.from_fractions({'H2':math.ulp(0.0),'O2':1.0})

    def test_component_and_mixture_records_refuse_invalid_or_inconsistent_fractions(self):
        s=get_species('H2')
        for number in (0,-1,True,None,float('nan'),float('inf'),1.1):
            with self.subTest(number=number),self.assertRaises(DomainError):MixtureComponent(s,number,1)
        for components in ((),[],None,'gas',{MixtureComponent(s,1,1)},(object(),)):
            with self.subTest(components=components),self.assertRaises(DomainError):Mixture(components)
        with self.assertRaises(DomainError):Mixture((MixtureComponent(s,.5,.5),MixtureComponent(s,.5,.5)))
        with self.assertRaises(DomainError):Mixture((MixtureComponent(s,.5,.5),MixtureComponent(get_species('O2'),.5,.5)))

    def test_reopened_mixture_refuses_displays_that_contradict_its_raw_inputs(self):
        mixture=Mixture.from_fractions({'N2':.79,'O2':.21})
        self.assertEqual(Mixture(mixture.components,mixture.composition_basis,mixture.input_fractions),mixture)
        n2,o2=mixture.components
        forged=(MixtureComponent(n2.species,math.nextafter(n2.mole_fraction,1),n2.mass_fraction),o2)
        with self.assertRaises(DomainError):Mixture(forged,mixture.composition_basis,mixture.input_fractions)

    def test_records_own_nested_species_and_composition_inputs(self):
        s=get_species('H2');component=MixtureComponent(s,1,1);source=[component]
        mixture=Mixture(source);source.clear()
        self.assertEqual(len(mixture.components),1)
        self.assertIsNot(component.species,s)
        self.assertIsNot(mixture.components[0],component)
        self.assertIsNot(mixture.components[0].species,component.species)
        self.assertIs(type(mixture.components),tuple)

    def test_order_cannot_change_any_exported_property(self):
        entries=(('H2',.2),('O2',.3),('H2O',.5));reference=None
        for order in itertools.permutations(entries):
            result=Mixture.from_fractions(dict(order)).properties(3500,pressure_pa=2e6)
            values=tuple(getattr(result,k) for k in oracle([(get_species(k),v) for k,v in entries],3500,2e6))
            if reference is None:reference=values
            self.assertEqual(values,reference)

    def test_pure_gas_entropy_pressure_change_matches_independent_log_difference(self):
        mixture=Mixture.from_fractions({'H2':1});a=mixture.properties(3500,pressure_pa=100_000);b=mixture.properties(3500,pressure_pa=1e6)
        with mp.workdps(110):expected=float(-mp.mpf(R)/mp.mpf(get_species('H2').molecular_mass_kg_mol)*mp.log(10))
        self.assertAlmostEqual(b.entropy_j_kg_k-a.entropy_j_kg_k,expected,delta=8*math.ulp(expected))
        self.assertEqual(a.enthalpy_j_kg,b.enthalpy_j_kg);self.assertEqual(a.cp_j_kg_k,b.cp_j_kg_k)

    def test_custom_component_overflow_can_cancel_before_mixture_output(self):
        a=custom('A',1,(0,0,1e308,0,0,0,0,0,1));b=custom('B',1,(0,0,3,0,0,0,0,0,1))
        mixture=Mixture((MixtureComponent(a,1e-100,1e-100),MixtureComponent(b,1,1)))
        self.assert_oracle(mixture,[(a,1e-100),(b,1)],1,100_000,digits=700)

    def test_custom_tiny_positive_cv_survives_display_cancellation(self):
        a=custom('A',1,(2.0**-54,0,1,0,0,0,0,0,1),lower=.5,upper=2)
        mixture=Mixture((MixtureComponent(a,1,1),));state=self.assert_oracle(mixture,[(a,1)],1,100_000)
        self.assertEqual(state.cp_j_kg_k,state.gas_constant_j_kg_k);self.assertGreater(state.cv_j_kg_k,0)

    def test_exact_cancellation_between_dependent_logarithms_is_preserved(self):
        a=custom('A',1,(0,0,3,0,0,0,0,0,0))
        mixture=Mixture((MixtureComponent(a,1,1),))
        try:state=mixture.properties(2,pressure_pa=800_000)
        except DomainError as error:self.fail(f'Exact zero entropy must remain representable: {error}')
        self.assertEqual(state.entropy_j_kg_k,0.0)
        b=replace(a,key='B')
        mixed=Mixture((MixtureComponent(a,.5,.5),MixtureComponent(b,.5,.5)))
        try:state=mixed.properties(2,pressure_pa=1_600_000)
        except DomainError as error:self.fail(f'Mixing/temperature/pressure logs cancel exactly: {error}')
        self.assertEqual(state.entropy_j_kg_k,0.0)

    def test_trace_mixing_entropy_remains_finite_and_nonzero(self):
        a=custom('A',1,(0,0,3,0,0,0,0,0,0));b=replace(a,key='B')
        mixture=Mixture((MixtureComponent(a,math.ulp(0.0),math.ulp(0.0)),MixtureComponent(b,1,1)))
        state=self.assert_oracle(mixture,[(a,math.ulp(0.0)),(b,1)],1,100_000,digits=700)
        self.assertGreater(state.entropy_j_kg_k,0)

    def test_invalid_component_cv_cannot_be_hidden_by_positive_mixture_average(self):
        a=custom('A',1,(0,0,.5,0,0,0,0,0,1));b=custom('B',1,(0,0,10,0,0,0,0,0,1))
        mixture=Mixture((MixtureComponent(a,.01,.01),MixtureComponent(b,.99,.99)))
        with self.assertRaises(DomainError):mixture.properties(300)

    def test_no_common_temperature_interval_is_refused_at_construction(self):
        a=custom('A',1,(0,0,3,0,0,0,0,0,1),lower=200,upper=500)
        b=custom('B',1,(0,0,3,0,0,0,0,0,1),lower=600,upper=1000)
        with self.assertRaises(DomainError):Mixture((MixtureComponent(a,.5,.5),MixtureComponent(b,.5,.5)))

    def test_shared_endpoint_is_a_valid_forward_state_even_without_inverse_bracket(self):
        a=custom('A',1,(0,0,3,0,0,0,0,0,1),lower=200,upper=1000)
        b=custom('B',1,(0,0,3,0,0,0,0,0,1),lower=1000,upper=6000)
        try:mixture=Mixture((MixtureComponent(a,.5,.5),MixtureComponent(b,.5,.5)))
        except DomainError as error:self.fail(f'Both owned fits admit exactly 1000 K: {error}')
        self.assert_oracle(mixture,[(a,.5),(b,.5)],1000,100_000)
        for t in (math.nextafter(1000,0),math.nextafter(1000,math.inf)):
            with self.subTest(t=t),self.assertRaises(DomainError):mixture.properties(t)

    def test_environment_context_changes_cannot_alter_mixture_exports(self):
        mixture=Mixture.from_fractions({'H2':.2,'O2':.8});expected=oracle([(get_species('H2'),.2),(get_species('O2'),.8)],3500,2e6)
        saved=DefaultContext.copy()
        try:
            DefaultContext.traps[Inexact]=True
            with localcontext() as context:
                context.prec=3;context.Emax=3;context.Emin=-3
                try:actual=mixture.properties(3500,pressure_pa=2e6)
                except (DomainError,ArithmeticError) as error:self.fail(f'Valid mixture depends on external decimal context: {error}')
            for k,v in expected.items():self.assertEqual(getattr(actual,k),v)
        finally:DefaultContext.traps=saved.traps.copy();DefaultContext.flags=saved.flags.copy()

    def test_inverse_targets_fail_at_boundary_before_evaluating_a_fit(self):
        from unittest.mock import patch
        mixture=Mixture.from_fractions({'H2':1})
        for value in (True,None,'1',float('inf'),float('nan')):
            for method in (mixture.temperature_from_enthalpy,mixture.temperature_from_entropy):
                with (patch.object(Mixture,'properties',side_effect=AssertionError('Invalid target reached fit')),
                      self.subTest(value=value,method=method.__name__),self.assertRaises(DomainError)):
                    method(value)

    def test_result_owns_inputs_and_refuses_one_ulp_or_qualification_forgery(self):
        mixture=Mixture.from_fractions({'H2':.2,'O2':.8});state=mixture.properties(3500,pressure_pa=2e6)
        self.assertIsNot(state.mixture,mixture)
        for field in fields(state):
            value=getattr(state,field.name)
            if type(value) is float and field.name not in {'temperature_k','pressure_pa'}:
                with self.subTest(field=field.name),self.assertRaises(DomainError):replace(state,**{field.name:math.nextafter(value,math.inf)})
        for value in (True,0,None):
            with self.subTest(qualified=value),self.assertRaises(DomainError):replace(state,physical_model_qualified=value)
        json.dumps(asdict(state),allow_nan=False)

    def test_component_resource_boundary_reads_at_most_65_records(self):
        component=MixtureComponent(get_species('H2'),1,1);reads=[]
        def records():
            for index in range(66):
                if index==65:raise AssertionError('Component iterator was read beyond budget')
                reads.append(index);yield component
        with self.assertRaises(DomainError):Mixture(records())
        self.assertEqual(reads,list(range(65)))

    def test_malformed_mapping_items_are_typed_and_mapping_read_is_bounded(self):
        from collections.abc import Mapping
        class Values(Mapping):
            def __iter__(self):return iter(())
            def __len__(self):return 10**12
            def __getitem__(self,key):return 1
            def items(self):return iter((('H2',1,2),))
        with self.assertRaises(DomainError):Mixture.from_fractions(Values())
        reads=[]
        class Endless(Values):
            def items(self):
                for index in range(66):
                    if index==65:raise AssertionError('Mapping was read beyond budget')
                    reads.append(index);yield ('H2',1)
        with self.assertRaises(DomainError):Mixture.from_fractions(Endless())
        self.assertEqual(reads,list(range(65)))

    def test_sixty_four_component_limit_includes_the_boundary(self):
        s=custom('S',1,(0,0,3,0,0,0,0,0,1))
        components=tuple(MixtureComponent(replace(s,key=f'S{i}'),1/64,1/64) for i in range(64))
        mixture=Mixture(components);state=mixture.properties(1)
        self.assertEqual(len(state.mixture.components),64)
        self.assertEqual(state.gamma,1.5)
        with mp.workdps(110):expected=float(mp.mpf(R)*(1+mp.log(64)))
        self.assertEqual(state.entropy_j_kg_k,expected)

    def test_nonzero_formation_offset_cannot_disappear_in_mass_conversion(self):
        s=custom('S',1e308,(math.ulp(0.0),0,3,0,0,0,0,-3,1))
        mixture=Mixture((MixtureComponent(s,1,1),))
        with self.assertRaises(DomainError):mixture.properties(1)

    def test_frozen_cantera_320_cross_implementation_states(self):
        import hashlib
        reference=json.loads((ROOT/'docs/assurance/mixture_cantera_reference.json').read_text())
        self.assertEqual(reference['cantera_version'],'3.2.0')
        self.assertEqual(reference['tolerances'],{'relative':2e-12,'cp_over_R_absolute':5e-12,
                                                'enthalpy_over_R_k_absolute':3e-8,'entropy_over_R_absolute':5e-11})
        self.assertEqual(reference['source_fixture_sha256'],hashlib.sha256((ROOT/'docs/assurance/nasa9_catalogue_reference.json').read_bytes()).hexdigest())
        self.assertEqual(reference['generator_sha256'],hashlib.sha256((ROOT/'scripts/generate_mixture_cantera_reference.py').read_bytes()).hexdigest())
        self.assertEqual(len(reference['cases']),150)
        for row in reference['cases']:
            mixture=Mixture.from_fractions(row['input_fractions'],basis=row['basis'])
            state=mixture.properties(row['temperature_k'],pressure_pa=row['pressure_pa'])
            values={'cp_over_R':state.cp_j_kg_k*state.molecular_mass_kg_mol/R,
                    'enthalpy_over_R_k':state.enthalpy_j_kg*state.molecular_mass_kg_mol/R,
                    'entropy_over_R':state.entropy_j_kg_k*state.molecular_mass_kg_mol/R}
            for name,value in values.items():
                with self.subTest(composition=row['input_fractions'],basis=row['basis'],t=row['temperature_k'],p=row['pressure_pa'],field=name):
                    self.assertTrue(math.isclose(value,row[name],rel_tol=reference['tolerances']['relative'],abs_tol=reference['tolerances'][name+'_absolute']),
                                    f'{name}: {value!r} versus Cantera {row[name]!r}')

