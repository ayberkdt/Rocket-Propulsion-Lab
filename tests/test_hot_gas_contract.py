"""Adversarial gas-state ownership and independent Bartz/recovery controls."""
import math
import unittest
from dataclasses import fields, replace
from decimal import Inexact, localcontext
from fractions import Fraction
from unittest.mock import patch

import mpmath as mp

from rocket_propulsion.core.errors import DomainError
from rocket_propulsion.heat_transfer.regen_chamber import (
    HotGas,
    bartz_coefficient,
    hot_gas_from_reactants,
    solve_regen_chamber,
)

BASE = {'chamber_pressure_pa': 7e6, 'chamber_temperature_k': 3500.,
            'characteristic_velocity_m_s': 1800., 'specific_heat_j_kg_k': 3200.,
            'viscosity_pa_s': 1e-4, 'prandtl_number': .7, 'gamma': 1.2,
            'throat_radius_of_curvature_m': .05, 'h2o_mole_fraction': .2, 'co2_mole_fraction': .1}
QUERY = {'throat_radius_m': .05, 'area_ratio': 1., 'mach': 1., 'wall_temperature_k': 900.}


class MutableScalar(float):
    def __new__(cls, value):
        obj = super().__new__(cls, value)
        obj.factor = 1.
        return obj

    def __float__(self):
        return super().__float__() * self.factor

    def __mul__(self, other):
        return float(self) * other


def reference(gas, query):
    """Log-space transcription with independent arbitrary precision arithmetic."""
    with mp.workdps(100):
        def exact(value):
            return mp.mpf(float(value))
        p, temp, cstar, cp, mu, pr, gamma, rc = (
            exact(getattr(gas, name)) for name in list(BASE)[:8])
        radius, area, mach, wall = (exact(query[name]) for name in QUERY)
        diameter = 2 * radius
        kinetic = (gamma - 1) * mach**2 / 2
        compression = 1 + kinetic
        sigma_base = (1 + wall * compression / temp) / 2
        log_h = (mp.log(mp.mpf('0.026')) + mp.log(cp) + mp.log(mu) / 5
                 - mp.mpf('0.6') * mp.log(pr) - mp.log(diameter) / 5
                 + mp.mpf('0.8') * (mp.log(p) - mp.log(cstar))
                 + (mp.log(diameter) - mp.log(rc)) / 10
                 - mp.mpf('0.9') * mp.log(area) - mp.mpf('0.68') * mp.log(sigma_base)
                 - mp.mpf('0.12') * mp.log(compression))
        # Preserve the kinetic term separately from 1 + kinetic for small M.
        recovery = temp * (1 + mp.root(pr, 3) * kinetic) / compression
        return float(mp.exp(log_h)), float(recovery)


class HotGasContractTests(unittest.TestCase):
    def test_gas_scalars_refuse_malformed_values_before_any_equation(self):
        for name in BASE:
            for value in (True, '1', None, (), complex(1), float('nan'), float('inf'),
                          Fraction(1, 10**500)):
                with self.subTest(field=name, value=repr(value)[:30]):
                    with self.assertRaises(DomainError) as caught:
                        HotGas(**(BASE | {name: value}))
                    self.assertIn(name, caught.exception.field or '')

    def test_physical_domains_include_gamma_and_composition_closure(self):
        for name in list(BASE)[:8]:
            for bad in (0., -1.):
                with self.subTest(field=name, value=bad), self.assertRaises(DomainError):
                    HotGas(**(BASE | {name: bad}))
        for gamma in (.9, 1.):
            with self.assertRaises(DomainError):
                HotGas(**(BASE | {'gamma': gamma}))
        for h2o, co2 in ((-.01, .1), (.1, -.01), (1.01, 0), (0, 1.01), (.8, .8)):
            with self.subTest(h2o=h2o, co2=co2), self.assertRaises(DomainError):
                HotGas(**(BASE | {'h2o_mole_fraction': h2o, 'co2_mole_fraction': co2}))
        for h2o, co2 in ((0., 0.), (1., 0.), (0., 1.), (.1, .9)):
            gas = HotGas(**(BASE | {'h2o_mole_fraction': h2o, 'co2_mole_fraction': co2}))
            self.assertEqual(gas.h2o_mole_fraction, h2o)

    def test_owned_native_floats_do_not_retain_caller_arithmetic(self):
        scalars = {name: MutableScalar(value) for name, value in BASE.items()}
        gas = HotGas(**scalars)
        for value in scalars.values():
            value.factor = 10.
        for field in fields(gas):
            self.assertIs(type(getattr(gas, field.name)), float)
            self.assertEqual(getattr(gas, field.name), BASE[field.name])
        expected = bartz_coefficient(HotGas(**BASE), **QUERY)
        self.assertEqual(bartz_coefficient(gas, **QUERY), expected)

    def test_direct_query_refuses_invalid_scalars_and_physical_domains(self):
        gas = HotGas(**BASE)
        for name in QUERY:
            for value in (True, '1', None, float('nan'), float('inf'), Fraction(1, 10**500), -1.):
                with self.subTest(field=name, value=repr(value)[:20]):
                    with self.assertRaises(DomainError) as caught:
                        bartz_coefficient(gas, **(QUERY | {name: value}))
                    self.assertEqual(caught.exception.field, name)
        for name in ('throat_radius_m', 'wall_temperature_k'):
            with self.assertRaises(DomainError):
                bartz_coefficient(gas, **(QUERY | {name: 0.}))
        for area in (0., .9):
            with self.assertRaises(DomainError):
                bartz_coefficient(gas, **(QUERY | {'area_ratio': area}))

    def test_station_query_owns_mutable_conversion_once(self):
        values = {name: MutableScalar(value) for name, value in QUERY.items()}
        self.assertEqual(bartz_coefficient(HotGas(**BASE), **values),
                         bartz_coefficient(HotGas(**BASE), **QUERY))

    def test_normal_operating_points_agree_with_independent_log_oracle(self):
        gas = HotGas(**BASE)
        for mach, area, wall in ((0., 4., 300.), (1., 1., 900.), (3., 5., 1400.)):
            query = QUERY | {'mach': mach, 'area_ratio': area, 'wall_temperature_k': wall}
            actual, expected = bartz_coefficient(gas, **query), reference(gas, query)
            for value, oracle in zip(actual, expected, strict=True):
                self.assertTrue(math.isclose(value, oracle, rel_tol=3e-14), (value, oracle))
        self.assertEqual(bartz_coefficient(gas, **(QUERY | {'mach': 0.}))[1], 3500.)

    def test_finite_results_survive_intermediate_overflow_and_underflow(self):
        for gas, query in (
            (HotGas(**(BASE | {'chamber_pressure_pa': 1e300,
                                  'characteristic_velocity_m_s': 1e-300,
                                  'specific_heat_j_kg_k': 1e-300})), QUERY),
            (HotGas(**BASE), QUERY | {'throat_radius_m': 1e308}),
            (HotGas(**(BASE | {'chamber_pressure_pa': 1e-323,
                              'characteristic_velocity_m_s': 3.,
                              'specific_heat_j_kg_k': 1e200})), QUERY),
            (HotGas(**(BASE | {'specific_heat_j_kg_k': 1e-320,
                              'viscosity_pa_s': 1e100})), QUERY),
        ):
            for value, oracle in zip(bartz_coefficient(gas, **query), reference(gas, query), strict=True):
                self.assertTrue(math.isfinite(value) and value > 0)
                self.assertTrue(math.isclose(value, oracle, rel_tol=5e-13), (value, oracle))

    def test_large_mach_recovery_and_film_coefficient_are_not_nan_or_false_zero(self):
        gas, query = HotGas(**(BASE | {'specific_heat_j_kg_k': 1e300})), QUERY | {'mach': 1e300}
        actual, expected = bartz_coefficient(gas, **query), reference(gas, query)
        for value, oracle in zip(actual, expected, strict=True):
            self.assertTrue(math.isfinite(value) and value > 0)
            self.assertTrue(math.isclose(value, oracle, rel_tol=5e-13), (value, oracle))

    def test_small_kinetic_term_is_not_erased_before_large_recovery_factor(self):
        gas = HotGas(**(BASE | {'prandtl_number': 1e300}))
        query = QUERY | {'mach': 1e-10}
        actual, expected = bartz_coefficient(gas, **query), reference(gas, query)
        self.assertGreater(actual[1], 1e80)
        for value, oracle in zip(actual, expected, strict=True):
            self.assertTrue(math.isclose(value, oracle, rel_tol=5e-13), (value, oracle))

    def test_extended_arithmetic_is_independent_of_caller_decimal_context(self):
        gas = HotGas(**(BASE | {'specific_heat_j_kg_k': 1e300}))
        query = QUERY | {'mach': 1e300}
        expected = bartz_coefficient(gas, **query)
        with localcontext() as context:
            context.prec = 5
            context.Emax = 9
            context.Emin = -9
            context.traps[Inexact] = True
            self.assertEqual(bartz_coefficient(gas, **query), expected)

    def test_unrepresentable_positive_outputs_raise_typed_numeric_range(self):
        for values in ({'specific_heat_j_kg_k': 5e-324, 'chamber_pressure_pa': 5e-324},
                       {'specific_heat_j_kg_k': 1e308, 'chamber_pressure_pa': 1e308,
                            'characteristic_velocity_m_s': 5e-324}):
            with self.subTest(values=values), self.assertRaises(DomainError) as caught:
                bartz_coefficient(HotGas(**(BASE | values)), **QUERY)
            self.assertEqual(caught.exception.code, 'numeric_range')

    def test_invalid_factory_radius_precedes_expensive_chemistry(self):
        for bad in (0., -1., True, '1', float('nan')):
            with patch('rocket_propulsion.thermochemistry.solve_rocket_performance') as solve:
                with self.subTest(radius=bad), self.assertRaises(DomainError):
                    hot_gas_from_reactants(reactants={'H2': 1., 'O2': .5},
                                          chamber_pressure_pa=7e6,
                                          throat_radius_of_curvature_m=bad)
                solve.assert_not_called()

    def test_real_coupled_heat_load_is_invariant_to_caller_gas_mutation(self):
        from tests.test_regen_channel_geometry import CONSTANT, contour, coolant, gas

        owned = gas()
        raw_cp = MutableScalar(owned.specific_heat_j_kg_k)
        candidate = replace(owned, specific_heat_j_kg_k=raw_cp)
        def solve():
            return solve_regen_chamber(contour=contour(), gas=candidate, channels=CONSTANT,
                                       coolant=coolant())
        before = solve()
        raw_cp.factor = 1.1
        after = solve()
        self.assertGreater(before.total_heat_load_w, 0.)
        for name in ('total_heat_load_w', 'coolant_outlet_temperature_k',
                     'coolant_outlet_pressure_pa', 'maximum_hot_wall_temperature_k'):
            self.assertTrue(math.isclose(getattr(before, name), getattr(after, name), rel_tol=1e-12))


if __name__ == '__main__':
    unittest.main()
