"""Direct-acting regulator characteristic: force balance, lockup, droop, hysteresis.

Oracles and sources
-------------------
* NASA SP-8080 (1973), printed p. 72 eq. (2): steady force balance
  ``P_A A_A + L - Y X - P_R A_A + P_R A_i - B - Z X - P_S A_i = 0`` and p. 73
  eq. (3): ``dP_R = -dP_S A_i/(A_A - A_i) - dX (Y + Z)/(A_A - A_i)``;
  p. 104 Glossary: lockup (all flow stops) and Coulomb damping (constant
  friction opposing motion); p. 76: Cd = 0.6 for preliminary sizing.
  https://ntrs.nasa.gov/citations/19740002611
* NASA Glenn "Mass Flow Choking", eq. 10 and the choked form, and "Isentropic
  Flow", eq. 6: https://www.grc.nasa.gov/www/k-12/BGP/mflchk.html,
  https://www.grc.nasa.gov/www/k-12/airplane/isentrop.html
* Regulator geometry and gas values below are synthetic, not hardware data.
"""

from __future__ import annotations

import unittest
from dataclasses import replace
from math import sqrt

from rocket_propulsion.core.errors import DomainError
from rocket_propulsion.propulsion.burns.feed_system import IdealPressurantGas
from rocket_propulsion.propulsion.burns.regulator_hysteresis import (
    SP8080_DEFAULT_DISCHARGE_COEFFICIENT,
    DirectActingRegulator,
    RegulatorBranch,
    compressible_valve_mass_flow_kg_s,
    feed_configuration_with_regulated_pressure,
    regulated_pressure_sensitivities,
    regulator_operating_point,
    sweep_regulator,
)
from tests.test_burn_feed_system import example_feed_system

GAS = IdealPressurantGas(2077.0, 1.667, "helium-like")
SUPPLY = 20.0e6
TEMPERATURE = 300.0


def regulator(friction: float = 10.0) -> DirectActingRegulator:
    return DirectActingRegulator(
        gas=GAS, actuator_area_m2=2.0e-3, inlet_valve_area_m2=2.0e-5,
        reference_spring_load_n=4400.0, reference_spring_rate_n_m=2.0e4,
        bias_spring_load_n=20.0, bias_spring_rate_n_m=1.0e3, reference_pressure_pa=0.0,
        maximum_stroke_m=1.0e-3, flow_area_per_stroke_m=3.0e-3, coulomb_friction_n=friction)


UNCHOKED_SUPPLY = 2.6e6


def unchoked_regulator() -> DirectActingRegulator:
    """Lockup 2.0 MPa at 2.6 MPa supply: the inlet valve stays subsonic."""

    return replace(regulator(friction=50.0), reference_spring_load_n=3982.0)


class ForceBalanceTests(unittest.TestCase):
    def test_branch_pressure_satisfies_sp8080_equation_2(self) -> None:
        reg = regulator()
        for branch, sign in ((RegulatorBranch.OPENING, 1.0), (RegulatorBranch.CLOSING, -1.0),
                             (RegulatorBranch.FRICTIONLESS, 0.0)):
            for stroke in (0.0, 2.0e-4, 9.0e-4):
                p_r = reg.branch_pressure_pa(stroke_m=stroke, supply_pressure_pa=SUPPLY,
                                             branch=branch)
                force = (0.0 * 2.0e-3 + 4400.0 - 2.0e4 * stroke - p_r * 2.0e-3 + p_r * 2.0e-5
                         - 20.0 - 1.0e3 * stroke - SUPPLY * 2.0e-5)
                self.assertAlmostEqual(force, sign * 10.0, delta=1e-6)

    def test_sensitivities_match_equation_3_by_finite_difference(self) -> None:
        reg = regulator()
        supply_slope, stroke_slope = regulated_pressure_sensitivities(reg)
        self.assertAlmostEqual(supply_slope, -2.0e-5 / (2.0e-3 - 2.0e-5))
        self.assertAlmostEqual(stroke_slope, -(2.0e4 + 1.0e3) / (2.0e-3 - 2.0e-5))
        p = reg.branch_pressure_pa
        h = 1.0e5
        fd_supply = (p(stroke_m=4e-4, supply_pressure_pa=SUPPLY + h, branch="opening")
                     - p(stroke_m=4e-4, supply_pressure_pa=SUPPLY - h, branch="opening")) / (2 * h)
        self.assertAlmostEqual(fd_supply, supply_slope, places=9)
        dx = 1.0e-5
        fd_stroke = (p(stroke_m=4e-4 + dx, supply_pressure_pa=SUPPLY, branch="closing")
                     - p(stroke_m=4e-4 - dx, supply_pressure_pa=SUPPLY, branch="closing")) / (2 * dx)
        self.assertAlmostEqual(fd_stroke / stroke_slope, 1.0, places=9)

    def test_default_discharge_coefficient_is_sp8080_value(self) -> None:
        self.assertEqual(SP8080_DEFAULT_DISCHARGE_COEFFICIENT, 0.6)
        self.assertEqual(regulator().discharge_coefficient, 0.6)


class FlowLawTests(unittest.TestCase):
    def test_choked_flow_matches_nasa_glenn_choked_equation(self) -> None:
        g, r = GAS.gamma, GAS.specific_gas_constant_j_kg_k
        area = 1.0e-6
        flow = compressible_valve_mass_flow_kg_s(gas=GAS, effective_area_m2=area,
                                                 upstream_pressure_pa=SUPPLY,
                                                 upstream_temperature_k=TEMPERATURE,
                                                 downstream_pressure_pa=1.0e6)
        expected = (area * SUPPLY / sqrt(TEMPERATURE) * sqrt(g / r)
                    * ((g + 1.0) / 2.0) ** (-(g + 1.0) / (g - 1.0) / 2.0))
        self.assertAlmostEqual(flow / expected, 1.0, places=12)

    def test_unchoked_flow_matches_mach_form_and_vanishes_at_equal_pressure(self) -> None:
        g, r = GAS.gamma, GAS.specific_gas_constant_j_kg_k
        downstream = 0.8 * SUPPLY
        mach = sqrt(2.0 / (g - 1.0) * ((SUPPLY / downstream) ** ((g - 1.0) / g) - 1.0))
        expected = (1.0e-6 * SUPPLY / sqrt(TEMPERATURE) * sqrt(g / r) * mach
                    * (1.0 + 0.5 * (g - 1.0) * mach * mach) ** (-(g + 1.0) / (g - 1.0) / 2.0))
        flow = compressible_valve_mass_flow_kg_s(gas=GAS, effective_area_m2=1.0e-6,
                                                 upstream_pressure_pa=SUPPLY,
                                                 upstream_temperature_k=TEMPERATURE,
                                                 downstream_pressure_pa=downstream)
        self.assertAlmostEqual(flow / expected, 1.0, places=12)
        self.assertEqual(compressible_valve_mass_flow_kg_s(
            gas=GAS, effective_area_m2=1.0e-6, upstream_pressure_pa=SUPPLY,
            upstream_temperature_k=TEMPERATURE, downstream_pressure_pa=SUPPLY), 0.0)


class CharacteristicTests(unittest.TestCase):
    def test_zero_demand_locks_up_at_the_closing_branch_zero_stroke_pressure(self) -> None:
        reg = regulator()
        point = regulator_operating_point(reg, mass_flow_kg_s=0.0, supply_pressure_pa=SUPPLY,
                                          supply_temperature_k=TEMPERATURE)
        expected = (4400.0 - 20.0 - SUPPLY * 2.0e-5 + 10.0) / (2.0e-3 - 2.0e-5)
        self.assertIs(point.branch, RegulatorBranch.LOCKUP)
        self.assertEqual(point.stroke_m, 0.0)
        self.assertAlmostEqual(point.regulated_pressure_pa, expected, places=6)
        self.assertEqual(point.droop_pa, 0.0)

    def test_regulated_pressure_droops_with_flow(self) -> None:
        reg = regulator(friction=0.0)
        points = sweep_regulator(reg, mass_flows_kg_s=(0.005, 0.01, 0.02, 0.03),
                                 supply_pressures_pa=SUPPLY, supply_temperature_k=TEMPERATURE)
        pressures = [p.regulated_pressure_pa for p in points]
        self.assertEqual(pressures, sorted(pressures, reverse=True))
        self.assertTrue(all(p.choked for p in points))
        # Choked flow is linear in stroke, so droop is linear in flow (eq. 3).
        slope_a = (pressures[1] - pressures[0]) / 0.005
        slope_b = (pressures[3] - pressures[2]) / 0.01
        self.assertAlmostEqual(slope_a / slope_b, 1.0, places=6)
        for point in points:
            flow = reg.mass_flow_kg_s(stroke_m=point.stroke_m, supply_pressure_pa=SUPPLY,
                                      supply_temperature_k=TEMPERATURE,
                                      regulated_pressure_pa=point.regulated_pressure_pa)
            self.assertAlmostEqual(flow / point.mass_flow_kg_s, 1.0, places=9)

    def test_frictionless_characteristic_is_path_independent(self) -> None:
        reg = regulator(friction=0.0)
        up = sweep_regulator(reg, mass_flows_kg_s=(0.0, 0.01, 0.02, 0.03),
                             supply_pressures_pa=SUPPLY, supply_temperature_k=TEMPERATURE)
        down = sweep_regulator(reg, mass_flows_kg_s=(0.03, 0.02, 0.01, 0.0),
                               supply_pressures_pa=SUPPLY, supply_temperature_k=TEMPERATURE)
        for a, b in zip(up, reversed(down), strict=True):
            self.assertAlmostEqual(a.regulated_pressure_pa, b.regulated_pressure_pa, delta=1e-3)

    def test_hysteresis_band_separates_opening_and_closing_branches(self) -> None:
        reg = regulator(friction=10.0)
        flows = (0.0, 0.01, 0.02, 0.03, 0.02, 0.01, 0.0)
        points = sweep_regulator(reg, mass_flows_kg_s=flows, supply_pressures_pa=SUPPLY,
                                 supply_temperature_k=TEMPERATURE)
        self.assertEqual([p.branch for p in points],
                         [RegulatorBranch.LOCKUP] + [RegulatorBranch.OPENING] * 3
                         + [RegulatorBranch.CLOSING] * 2 + [RegulatorBranch.LOCKUP])
        band = 2.0 * 10.0 / (2.0e-3 - 2.0e-5)
        self.assertAlmostEqual(reg.hysteresis_band_pa, band)
        for opening, closing in ((points[1], points[5]), (points[2], points[4])):
            self.assertAlmostEqual(opening.stroke_m, closing.stroke_m, delta=1e-12)
            self.assertAlmostEqual(closing.regulated_pressure_pa - opening.regulated_pressure_pa,
                                   band, delta=1e-3)

    def test_unchoked_small_reversal_stays_stuck_inside_the_band(self) -> None:
        reg = unchoked_regulator()
        first = regulator_operating_point(reg, mass_flow_kg_s=0.003,
                                          supply_pressure_pa=UNCHOKED_SUPPLY,
                                          supply_temperature_k=TEMPERATURE)
        self.assertFalse(first.choked)
        self.assertIs(first.branch, RegulatorBranch.OPENING)
        second = regulator_operating_point(reg, mass_flow_kg_s=0.00299,
                                           supply_pressure_pa=UNCHOKED_SUPPLY,
                                           supply_temperature_k=TEMPERATURE,
                                           previous_stroke_m=first.stroke_m)
        self.assertIs(second.branch, RegulatorBranch.STUCK)
        self.assertEqual(second.stroke_m, first.stroke_m)
        self.assertGreater(second.regulated_pressure_pa, first.regulated_pressure_pa)
        self.assertLess(second.regulated_pressure_pa - first.regulated_pressure_pa,
                        reg.hysteresis_band_pa)
        third = regulator_operating_point(reg, mass_flow_kg_s=0.002,
                                          supply_pressure_pa=UNCHOKED_SUPPLY,
                                          supply_temperature_k=TEMPERATURE,
                                          previous_stroke_m=second.stroke_m)
        self.assertIs(third.branch, RegulatorBranch.CLOSING)
        self.assertLess(third.stroke_m, second.stroke_m)

    def test_demand_beyond_full_open_capacity_fails_closed(self) -> None:
        with self.assertRaises(DomainError):
            regulator_operating_point(regulator(), mass_flow_kg_s=1.0, supply_pressure_pa=SUPPLY,
                                      supply_temperature_k=TEMPERATURE)

    def test_saturated_point_is_flagged(self) -> None:
        reg = replace(unchoked_regulator(), maximum_stroke_m=5.0e-5)
        point = regulator_operating_point(reg, mass_flow_kg_s=2.0e-4,
                                          supply_pressure_pa=UNCHOKED_SUPPLY,
                                          supply_temperature_k=TEMPERATURE)
        self.assertTrue(point.saturated)
        self.assertEqual(point.stroke_m, 5.0e-5)
        lower = reg.branch_pressure_pa(stroke_m=5.0e-5, supply_pressure_pa=UNCHOKED_SUPPLY,
                                       branch=RegulatorBranch.OPENING)
        self.assertLess(point.regulated_pressure_pa, lower)


class FeedCouplingAndErrorTests(unittest.TestCase):
    def test_operating_point_drives_the_regulated_feed_set_pressure(self) -> None:
        config = example_feed_system().configuration
        point = regulator_operating_point(regulator(), mass_flow_kg_s=0.01,
                                          supply_pressure_pa=SUPPLY,
                                          supply_temperature_k=TEMPERATURE)
        updated = feed_configuration_with_regulated_pressure(config, point)
        self.assertEqual(updated.regulator_set_pressure_pa, point.regulated_pressure_pa)
        self.assertEqual(replace(updated, regulator_set_pressure_pa=config.regulator_set_pressure_pa),
                         config)
        # A regulator delivering more than the tank overpressure limit is refused.
        high = regulator_operating_point(replace(regulator(), reference_spring_load_n=5000.0),
                                         mass_flow_kg_s=0.01, supply_pressure_pa=SUPPLY,
                                         supply_temperature_k=TEMPERATURE)
        self.assertGreater(high.regulated_pressure_pa, config.maximum_tank_pressure_pa)
        with self.assertRaises(DomainError):
            feed_configuration_with_regulated_pressure(config, high)
        with self.assertRaises(DomainError):
            feed_configuration_with_regulated_pressure(config, "point")

    def test_invalid_regulators_and_inputs_fail_closed(self) -> None:
        base = regulator()
        for kwargs in (dict(inlet_valve_area_m2=2.0e-3), dict(actuator_area_m2=-1.0),
                       dict(reference_spring_rate_n_m=0.0, bias_spring_rate_n_m=0.0),
                       dict(coulomb_friction_n=-1.0), dict(discharge_coefficient=1.2),
                       dict(maximum_flow_area_m2=0.0), dict(model_id=" "), dict(schema="x"),
                       dict(gas="helium"), dict(maximum_stroke_m=float("nan"))):
            with self.subTest(kwargs=kwargs), self.assertRaises(DomainError):
                replace(base, **kwargs)
        for kwargs in (dict(mass_flow_kg_s=-1.0), dict(supply_pressure_pa=0.0),
                       dict(supply_temperature_k=float("inf")), dict(previous_stroke_m=1.0)):
            arguments = dict(mass_flow_kg_s=0.01, supply_pressure_pa=SUPPLY,
                             supply_temperature_k=TEMPERATURE) | kwargs
            with self.subTest(kwargs=kwargs), self.assertRaises(DomainError):
                regulator_operating_point(base, **arguments)
        with self.assertRaises(DomainError):  # lockup above supply
            regulator_operating_point(base, mass_flow_kg_s=0.01, supply_pressure_pa=1.0e6,
                                      supply_temperature_k=TEMPERATURE)
        with self.assertRaises(DomainError):
            sweep_regulator(base, mass_flows_kg_s=(), supply_pressures_pa=SUPPLY,
                            supply_temperature_k=TEMPERATURE)
        with self.assertRaises(DomainError):
            sweep_regulator(base, mass_flows_kg_s=(0.0, 0.01), supply_pressures_pa=(SUPPLY,),
                            supply_temperature_k=TEMPERATURE)


if __name__ == "__main__":
    unittest.main()
