"""Finite ideal-gas pressures and Jacobians across hostile volume scales.

Decimal oracles use p=mRT/V and its analytic derivatives (NASA SP-8112).
Synthetic range checks are numerical verification, not a real tank design.
"""

from __future__ import annotations

import math
import unittest
from dataclasses import replace
from decimal import Decimal, localcontext

from rocket_propulsion.propulsion.burns import RegulatedFeedSystem
from tests.test_burn_feed_system import example_feed_system


class FeedPressureRangeTests(unittest.TestCase):
    def test_representable_pressure_partials_survive_overflowing_ullage_square(self):
        for volume in (.02, 1e160, 1e300):
            system = RegulatedFeedSystem(replace(example_feed_system().configuration,
                tank_internal_volume_m3=volume))
            state = system.initial_state(propellant_mass_kg=10., tank_pressure_pa=2e6,
                tank_gas_temperature_k=300., bottle_pressure_pa=20e6,
                bottle_gas_temperature_k=300., regulator_opening=0.)
            with self.subTest(volume=volume):
                result = system.evaluate(state, .05)
                with localcontext() as context:
                    context.prec = 90
                    mass, gas_r, temp, rho, ullage = map(Decimal,
                        (state.tank_pressurant_mass_kg,
                         system.configuration.gas.specific_gas_constant_j_kg_k,
                         state.tank_gas_temperature_k,
                         system.configuration.propellant_density_kg_m3,
                         result.tank_ullage_volume_m3))
                    expected = {'propellant_mass_kg':mass*gas_r*temp/(rho*ullage**2),
                                'tank_pressurant_mass_kg':gas_r*temp/ullage,
                                'tank_gas_temperature_k':mass*gas_r/ullage}
                actual = dict(result.pressure_partials.injector_pressure_wrt_state)
                for field,value in expected.items():
                    self.assertGreater(actual[field], 0., field)
                    self.assertTrue(math.isclose(actual[field], float(value), rel_tol=5e-14),
                                    (field, actual[field], value))

    def test_density_compensation_cannot_be_lost_by_squaring_volume_first(self):
        system = RegulatedFeedSystem(replace(example_feed_system().configuration,
            tank_internal_volume_m3=1e200, propellant_density_kg_m3=1e-200,
            minimum_propellant_mass_kg=1e-12))
        state = system.initial_state(propellant_mass_kg=1e-10, tank_pressure_pa=2e6,
            tank_gas_temperature_k=300., bottle_pressure_pa=20e6,
            bottle_gas_temperature_k=300., regulator_opening=0.)
        result = system.evaluate(state, 0.)
        with localcontext() as context:
            context.prec = 90
            expected = (Decimal(state.tank_pressurant_mass_kg)
                * Decimal(system.configuration.gas.specific_gas_constant_j_kg_k)
                * Decimal.from_float(300.)
                / (Decimal.from_float(1e-200)*Decimal(result.tank_ullage_volume_m3)**2))
        actual = dict(result.pressure_partials.injector_pressure_wrt_state)['propellant_mass_kg']
        self.assertTrue(math.isclose(actual, float(expected), rel_tol=5e-14))

    def test_initial_mass_and_recovered_pressure_survive_overflowing_pv_product(self):
        system = RegulatedFeedSystem(replace(example_feed_system().configuration,
            tank_internal_volume_m3=1e307))
        state = system.initial_state(propellant_mass_kg=10., tank_pressure_pa=2e6,
            tank_gas_temperature_k=300., bottle_pressure_pa=20e6,
            bottle_gas_temperature_k=300., regulator_opening=0.)
        with localcontext() as context:
            context.prec = 90
            expected = Decimal.from_float(2e6)*Decimal.from_float(1e307)/(
                Decimal(system.configuration.gas.specific_gas_constant_j_kg_k)
                * Decimal.from_float(300.))
        self.assertTrue(math.isclose(state.tank_pressurant_mass_kg,float(expected),rel_tol=5e-14))
        _,pressure,_ = system._pressures(state)
        self.assertTrue(math.isclose(pressure,2e6,rel_tol=5e-14))

    def test_partial_sign_and_mass_pressure_elasticities_remain_physical(self):
        system = example_feed_system()
        state = system.initial_state(propellant_mass_kg=10., tank_pressure_pa=2e6,
            tank_gas_temperature_k=300., bottle_pressure_pa=20e6,
            bottle_gas_temperature_k=300., regulator_opening=0.)
        result = system.evaluate(state,.05)
        jac = dict(result.pressure_partials.injector_pressure_wrt_state)
        self.assertGreater(jac['propellant_mass_kg'],0.)
        self.assertTrue(math.isclose(jac['tank_pressurant_mass_kg']*state.tank_pressurant_mass_kg,
                                    result.tank_pressure_pa,rel_tol=5e-14))
        self.assertTrue(math.isclose(jac['tank_gas_temperature_k']*state.tank_gas_temperature_k,
                                    result.tank_pressure_pa,rel_tol=5e-14))
        self.assertLess(result.pressure_partials.injector_pressure_wrt_propellant_flow,0.)


if __name__ == '__main__':
    unittest.main()
