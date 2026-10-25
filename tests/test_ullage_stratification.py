"""Stratified ullage: limits, conservation, convergence, comparisons and errors.

Oracles and sources
-------------------
* Ideal (no heat transfer) pressurant mass ``m_i = p dV / (R T_g)`` and the
  mass ratio (collapse factor): NASA TN D-2797, printed p. 8, eqs. 23-25,
  https://ntrs.nasa.gov/citations/19650012768
* Natural convection ``h = K_H C (k/L) (Gr Pr)**n``; ullage-wall C = 0.54,
  ullage-liquid C = 0.27, n = 0.25, K_H = 1, L = tank diameter: Majumdar and
  Steadman, AIAA 99-0879 draft, pp. 4-5, https://ntrs.nasa.gov/citations/19980236009
  ``Gr = g beta dT L**3 rho**2 / mu**2``: NASA 20190030427, PDF p. 1,
  https://ntrs.nasa.gov/citations/20190030427
* SP-8112 eq. (14) first-cut ullage mass at the coldest temperature and
  p. 90 statement that stratification lowers heat transfer to the propellant:
  https://ntrs.nasa.gov/citations/19760015212
* Adiabatic ideal-gas limits (``p V**gamma`` constant for isentropic
  expansion; ``d(pV)/dt = gamma R T_in mdot`` for adiabatic filling) follow
  from the specific-heat and isentropic relations of NASA Glenn:
  https://www.grc.nasa.gov/www/k-12/airplane/specheat.html,
  https://www.grc.nasa.gov/www/k-12/airplane/isentrop.html
* Tank, gas transport and temperatures are synthetic, not hardware data.
"""

from __future__ import annotations

import itertools
import time
import unittest
from dataclasses import replace

from rocket_propulsion.core.errors import DomainError
from rocket_propulsion.propulsion.burns.feed_system import IdealPressurantGas
from rocket_propulsion.propulsion.burns.ullage_stratification import (
    ULLAGE_LIQUID_CORRELATION,
    ULLAGE_WALL_CORRELATION,
    CylindricalTankWall,
    NaturalConvectionCorrelation,
    UllageConfiguration,
    UllageGasTransport,
    compare_ullage_models,
    natural_convection_coefficient_w_m2_k,
    simulate_ullage,
)

GAS = IdealPressurantGas(2077.0, 1.667, "helium-like")
TRANSPORT = UllageGasTransport(0.06, 1.0e-5)
TANK = CylindricalTankWall(0.5, 3.0, 0.003, 2700.0, 900.0)


def configuration(**changes) -> UllageConfiguration:
    base = UllageConfiguration(
        gas=GAS, transport=TRANSPORT, tank=TANK, initial_ullage_height_m=0.3,
        initial_pressure_pa=3.0e5, initial_gas_temperature_k=100.0,
        initial_wall_temperature_k=90.0, liquid_temperature_k=90.0, inlet_temperature_k=300.0,
        liquid_outflow_m3_s=0.02, acceleration_m_s2=10.0, duration_s=60.0, step_count=120,
        set_pressure_pa=3.0e5)
    return replace(base, **changes)


def adiabatic(**changes) -> UllageConfiguration:
    return configuration(wall_correlation=None, interface_correlation=None, **changes)


class CorrelationTests(unittest.TestCase):
    def test_constants_are_the_cited_values(self) -> None:
        self.assertEqual((ULLAGE_WALL_CORRELATION.coefficient, ULLAGE_WALL_CORRELATION.exponent,
                          ULLAGE_WALL_CORRELATION.adjustment_factor), (0.54, 0.25, 1.0))
        self.assertEqual((ULLAGE_LIQUID_CORRELATION.coefficient,
                          ULLAGE_LIQUID_CORRELATION.exponent), (0.27, 0.25))

    def test_coefficient_matches_hand_evaluation(self) -> None:
        p, t_gas, t_wall, length, g = 3.0e5, 250.0, 90.0, 1.0, 10.0
        rho = p / (GAS.specific_gas_constant_j_kg_k * t_gas)
        grashof = g * (1.0 / t_gas) * (t_gas - t_wall) * length ** 3 * rho ** 2 / 1.0e-5 ** 2
        prandtl = 1.0e-5 * GAS.specific_heat_cp_j_kg_k / 0.06
        expected = 0.54 * 0.06 / length * (grashof * prandtl) ** 0.25
        value = natural_convection_coefficient_w_m2_k(
            ULLAGE_WALL_CORRELATION, gas=GAS, transport=TRANSPORT, pressure_pa=p,
            gas_temperature_k=t_gas, surface_temperature_k=t_wall, length_scale_m=length,
            acceleration_m_s2=g)
        self.assertAlmostEqual(value / expected, 1.0, places=12)
        halved = replace(ULLAGE_WALL_CORRELATION, adjustment_factor=0.5)
        self.assertAlmostEqual(natural_convection_coefficient_w_m2_k(
            halved, gas=GAS, transport=TRANSPORT, pressure_pa=p, gas_temperature_k=t_gas,
            surface_temperature_k=t_wall, length_scale_m=length, acceleration_m_s2=g) / value, 0.5)

    def test_no_buoyancy_or_no_temperature_difference_gives_zero(self) -> None:
        common = {"gas": GAS, "transport": TRANSPORT, "pressure_pa": 3e5, "length_scale_m": 1.0}
        self.assertEqual(natural_convection_coefficient_w_m2_k(
            ULLAGE_WALL_CORRELATION, gas_temperature_k=200.0, surface_temperature_k=90.0,
            acceleration_m_s2=0.0, **common), 0.0)
        self.assertEqual(natural_convection_coefficient_w_m2_k(
            ULLAGE_WALL_CORRELATION, gas_temperature_k=90.0, surface_temperature_k=90.0,
            acceleration_m_s2=10.0, **common), 0.0)


class AdiabaticLimitTests(unittest.TestCase):
    def test_regulated_adiabatic_requirement_equals_ideal_mass(self) -> None:
        comparison = compare_ullage_models(adiabatic())
        ideal = 3.0e5 * 0.02 * 60.0 / (GAS.specific_gas_constant_j_kg_k * 300.0)
        self.assertAlmostEqual(comparison.ideal_added_mass_kg / ideal, 1.0, places=12)
        self.assertAlmostEqual(comparison.stratified_collapse_factor, 1.0, places=10)
        self.assertAlmostEqual(comparison.well_mixed_collapse_factor, 1.0, places=10)
        for sample in comparison.stratified.samples:
            self.assertAlmostEqual(sample.pressure_pa, 3.0e5, delta=1e-6)

    def test_isentropic_expansion_without_inflow(self) -> None:
        config = adiabatic(set_pressure_pa=None, step_count=40)
        for mixed in (False, True):
            history = simulate_ullage(config, well_mixed=mixed)
            gamma = GAS.gamma
            first = history.samples[0]
            invariant = first.pressure_pa * first.ullage_volume_m3 ** gamma
            for sample in history.samples:
                self.assertAlmostEqual(
                    sample.pressure_pa * sample.ullage_volume_m3 ** gamma / invariant, 1.0,
                    places=12)
            self.assertEqual(history.added_mass_kg, 0.0)

    def test_closed_tank_adiabatic_filling_converges_to_linear_pressure_rise(self) -> None:
        mdot = 1.0e-3
        errors, residuals = [], []
        for steps in (25, 50, 100):
            config = adiabatic(set_pressure_pa=None, prescribed_inflow_kg_s=mdot,
                               liquid_outflow_m3_s=0.0, step_count=steps)
            history = simulate_ullage(config)
            volume = history.samples[0].ullage_volume_m3
            slope = GAS.gamma * GAS.specific_gas_constant_j_kg_k * 300.0 * mdot / volume
            errors.append(max(abs(s.pressure_pa / (3.0e5 + slope * s.time_s) - 1.0)
                              for s in history.samples))
            residuals.append(abs(history.energy_balance_residual_j) / history.energy_scale_j)
            self.assertAlmostEqual(history.added_mass_kg, mdot * 60.0, places=12)
        # First order in dt: halving the step halves both the pressure error and
        # the energy residual (the parcel enters at the end-of-step pressure).
        self.assertLess(errors[-1], 1e-3)
        for coarse, fine in itertools.pairwise(errors):
            self.assertAlmostEqual(fine / coarse, 0.5, delta=0.02)
        for coarse, fine in itertools.pairwise(residuals):
            self.assertAlmostEqual(fine / coarse, 0.5, delta=0.02)


class ConservationAndConvergenceTests(unittest.TestCase):
    def test_mass_energy_and_volume_close_with_heat_transfer(self) -> None:
        for mixed in (False, True):
            history = simulate_ullage(configuration(), well_mixed=mixed)
            initial = history.samples[0]
            for sample in history.samples:
                self.assertAlmostEqual(sample.gas_mass_kg,
                                       initial.gas_mass_kg + sample.added_mass_kg, places=12)
            self.assertLess(abs(history.energy_balance_residual_j),
                            1e-9 * history.energy_scale_j)
            last = history.samples[-1]
            self.assertGreater(last.heat_to_wall_j, 0.0)
            self.assertGreater(last.heat_to_liquid_j, 0.0)
            self.assertAlmostEqual(last.ullage_volume_m3,
                                   TANK.cross_section_m2 * 0.3 + 0.02 * 60.0, places=12)

    def test_added_mass_converges_with_time_step(self) -> None:
        masses = [simulate_ullage(configuration(step_count=n)).added_mass_kg for n in (40, 80, 160)]
        first, second = abs(masses[1] - masses[0]), abs(masses[2] - masses[1])
        self.assertLess(second, 0.7 * first)
        self.assertLess(second / masses[2], 2e-3)

    def test_stratified_profile_is_hot_on_top_and_bounded(self) -> None:
        history = simulate_ullage(configuration())
        temperatures = history.final_layer_temperatures_k
        self.assertGreater(temperatures[0], temperatures[-1])
        self.assertLessEqual(max(temperatures), 300.0 + 1e-9)
        self.assertGreaterEqual(min(temperatures), 90.0 - 1e-9)
        walls = history.final_wall_temperatures_k
        self.assertTrue(all(90.0 - 1e-9 <= t <= 300.0 for t in walls))
        self.assertEqual(walls[-1], 90.0)  # still wetted: held at liquid temperature

    def test_layer_cap_merges_but_keeps_mass(self) -> None:
        capped = simulate_ullage(configuration(maximum_layer_count=20))
        self.assertGreater(capped.layer_merges, 0)
        self.assertEqual(len(capped.final_layer_temperatures_k), 20)
        free = simulate_ullage(configuration())
        self.assertAlmostEqual(capped.added_mass_kg / free.added_mass_kg, 1.0, delta=0.02)

    def test_runs_fast(self) -> None:
        start = time.perf_counter()
        compare_ullage_models(configuration())
        self.assertLess(time.perf_counter() - start, 2.0)


class ComparisonTests(unittest.TestCase):
    def test_heat_loss_raises_requirement_above_ideal_and_below_sp8112_bound(self) -> None:
        comparison = compare_ullage_models(configuration())
        self.assertGreater(comparison.stratified_collapse_factor, 1.0)
        self.assertGreater(comparison.well_mixed_collapse_factor, 1.0)
        for history in (comparison.stratified, comparison.well_mixed):
            self.assertLess(history.samples[-1].gas_mass_kg, comparison.sp8112_ullage_mass_kg)

    def test_stratification_reduces_heat_to_liquid_and_pressurant_need(self) -> None:
        # Interface-dominated case (adiabatic wall): SP-8112 p. 90 expectation.
        comparison = compare_ullage_models(configuration(wall_correlation=None))
        stratified, mixed = comparison.stratified, comparison.well_mixed
        self.assertLess(stratified.samples[-1].heat_to_liquid_j,
                        mixed.samples[-1].heat_to_liquid_j)
        self.assertLess(stratified.added_mass_kg, mixed.added_mass_kg)
        self.assertGreater(stratified.final_layer_temperatures_k[0],
                           mixed.final_layer_temperatures_k[0])

    def test_constant_wall_coefficient_option(self) -> None:
        zero = simulate_ullage(configuration(wall_correlation=None,
                                             wall_heat_transfer_coefficient_w_m2_k=0.0))
        none = simulate_ullage(configuration(wall_correlation=None))
        self.assertAlmostEqual(zero.added_mass_kg, none.added_mass_kg, places=12)
        hot = simulate_ullage(configuration(wall_correlation=None,
                                            wall_heat_transfer_coefficient_w_m2_k=20.0))
        self.assertGreater(hot.added_mass_kg, none.added_mass_kg)

    def test_heating_above_set_pressure_locks_out_inflow(self) -> None:
        # Warm wall and liquid heat a cold ullage: pressure rises above the set
        # value and no pressurant is admitted while the outflow is too small.
        config = configuration(initial_gas_temperature_k=100.0, initial_wall_temperature_k=300.0,
                               liquid_temperature_k=300.0, liquid_outflow_m3_s=1.0e-5)
        history = simulate_ullage(config)
        self.assertGreater(history.samples[-1].pressure_pa, 3.0e5)
        self.assertEqual(history.added_mass_kg, 0.0)


class UllageErrorTests(unittest.TestCase):
    def test_invalid_configurations_fail_closed(self) -> None:
        for changes in (
            {"initial_ullage_height_m": 3.0}, {"initial_pressure_pa": 0.0},
            {"liquid_outflow_m3_s": -1.0}, {"step_count": 0}, {"duration_s": float("nan")},
            {"prescribed_inflow_kg_s": 1.0}, {"wall_heat_transfer_coefficient_w_m2_k": 5.0},
            {"liquid_outflow_m3_s": 1.0}, {"maximum_layer_count": 2}, {"schema": "x"},
            {"wall_correlation": "0.54"}, {"acceleration_m_s2": -1.0}, {"gas": None},
        ):
            with self.subTest(changes=changes), self.assertRaises(DomainError):
                configuration(**changes)
        with self.assertRaises(DomainError):
            NaturalConvectionCorrelation(0.54, 1.0)
        with self.assertRaises(DomainError):
            UllageGasTransport(0.0, 1e-5)
        with self.assertRaises(DomainError):
            CylindricalTankWall(0.5, 3.0, -0.003, 2700.0, 900.0)
        with self.assertRaises(DomainError):
            compare_ullage_models(configuration(set_pressure_pa=4.0e5))
        with self.assertRaises(DomainError):
            compare_ullage_models(configuration(liquid_outflow_m3_s=0.0))
        with self.assertRaises(DomainError):
            simulate_ullage("configuration")
        with self.assertRaises(DomainError):
            natural_convection_coefficient_w_m2_k(
                ULLAGE_WALL_CORRELATION, gas=GAS, transport=TRANSPORT, pressure_pa=3e5,
                gas_temperature_k=200.0, surface_temperature_k=90.0, length_scale_m=0.0,
                acceleration_m_s2=1.0)


if __name__ == "__main__":
    unittest.main()
