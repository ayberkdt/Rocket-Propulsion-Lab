"""Specific-speed definitions, similarity identities, NPSH margins and the homologous pump map."""

import unittest
from math import pi

from rocket_propulsion.core.constants import STANDARD_GRAVITY_M_S2
from rocket_propulsion.core.errors import DomainError, InputError
from rocket_propulsion.propulsion import (
    PumpCharacteristic,
    calculate_engine_cycle,
    evaluate_cycle_pump_point,
    evaluate_pump_map,
    evaluate_pump_suction,
    net_positive_suction_head,
    pump_specific_speed,
    pump_suction_specific_speed,
    scale_pump_affinity,
)

GPM = 6.309020e-05  # NIST SP 811 B.9
FOOT = 0.3048  # NIST SP 811 B.9
RPM = 2 * pi / 60


def characteristic(**changes):
    values = dict(design_flow_m3_s=0.1, design_head_m=1000.0, design_speed_rad_s=3000.0,
                  shutoff_coefficient=1.2, linear_coefficient=0.1, quadratic_coefficient=-0.3)
    values.update(changes)
    return PumpCharacteristic(**values)


class SpecificSpeed(unittest.TestCase):
    def test_us_customary_definition(self):
        value = pump_specific_speed(shaft_speed_rad_s=3600 * RPM, flow_m3_s=100 * GPM,
                                    head_m=100 * FOOT)
        self.assertAlmostEqual(value, 3600 * 100 ** 0.5 / 100 ** 0.75, delta=1e-9 * value)
        suction = pump_suction_specific_speed(shaft_speed_rad_s=10000 * RPM, flow_m3_s=400 * GPM,
                                              npsh_m=16 * FOOT)
        self.assertAlmostEqual(suction, 10000 * 20 / 8, delta=1e-9 * suction)

    def test_specific_speed_is_invariant_under_homologous_speed_change(self):
        base = scale_pump_affinity(flow_m3_s=0.05, head_m=800.0, shaft_power_w=5e5, speed_ratio=1.0)
        for ratio in (0.5, 0.8, 1.3):
            scaled = scale_pump_affinity(flow_m3_s=0.05, head_m=800.0, shaft_power_w=5e5,
                                         speed_ratio=ratio)
            ns0 = pump_specific_speed(shaft_speed_rad_s=2000.0, flow_m3_s=base.flow_m3_s,
                                      head_m=base.head_m)
            ns1 = pump_specific_speed(shaft_speed_rad_s=2000.0 * ratio, flow_m3_s=scaled.flow_m3_s,
                                      head_m=scaled.head_m)
            self.assertAlmostEqual(ns0, ns1, delta=1e-12 * ns0)

    def test_inputs_must_be_finite_and_positive(self):
        for name in ("shaft_speed_rad_s", "flow_m3_s", "head_m"):
            for value in (0.0, -1.0, float("inf"), True, None):
                args = dict(shaft_speed_rad_s=1000.0, flow_m3_s=0.1, head_m=100.0)
                args[name] = value
                with self.subTest(field=name, value=value), self.assertRaises(DomainError):
                    pump_specific_speed(**args)


class Affinity(unittest.TestCase):
    def test_speed_and_trim_exponents(self):
        r = scale_pump_affinity(flow_m3_s=2.0, head_m=10.0, shaft_power_w=100.0, speed_ratio=2.0,
                                diameter_ratio=0.9)
        self.assertAlmostEqual(r.flow_m3_s, 2.0 * 2.0 * 0.9)
        self.assertAlmostEqual(r.head_m, 10.0 * (2.0 * 0.9) ** 2)
        self.assertAlmostEqual(r.shaft_power_w, 100.0 * (2.0 * 0.9) ** 3)
        self.assertFalse(r.physical_model_qualified)

    def test_inverse_scaling_returns_the_original_point(self):
        forward = scale_pump_affinity(flow_m3_s=0.3, head_m=500.0, shaft_power_w=2e6,
                                      speed_ratio=1.25, diameter_ratio=0.95)
        back = scale_pump_affinity(flow_m3_s=forward.flow_m3_s, head_m=forward.head_m,
                                   shaft_power_w=forward.shaft_power_w, speed_ratio=1 / 1.25,
                                   diameter_ratio=1 / 0.95)
        for name, original in (("flow_m3_s", 0.3), ("head_m", 500.0), ("shaft_power_w", 2e6)):
            self.assertAlmostEqual(getattr(back, name), original, delta=1e-12 * original)

    def test_trim_scope_warnings(self):
        def codes(d):
            result = scale_pump_affinity(flow_m3_s=1.0, head_m=1.0, shaft_power_w=1.0,
                                         speed_ratio=1.0, diameter_ratio=d)
            return [w.code for w in result.warnings]
        self.assertIn("impeller_enlargement_outside_trim_relation", codes(1.1))
        self.assertIn("impeller_trim_below_doe_practice", codes(0.69))
        self.assertNotIn("impeller_trim_below_doe_practice", codes(0.7))
        with self.assertRaises(DomainError):
            scale_pump_affinity(flow_m3_s=1.0, head_m=1.0, shaft_power_w=1.0, speed_ratio=0.0)


class Suction(unittest.TestCase):
    def test_npsh_available_definition(self):
        value = net_positive_suction_head(inlet_total_pressure_pa=3e5, vapor_pressure_pa=1e5,
                                          density_kg_m3=1000.0)
        self.assertEqual(value, 2e5 / (1000.0 * STANDARD_GRAVITY_M_S2))
        below = net_positive_suction_head(inlet_total_pressure_pa=1e5, vapor_pressure_pa=2e5,
                                          density_kg_m3=1000.0)
        self.assertLess(below, 0)
        for args in ({"vapor_pressure_pa": -1.0}, {"density_kg_m3": 0.0},
                     {"inlet_total_pressure_pa": 0.0}):
            base = dict(inlet_total_pressure_pa=3e5, vapor_pressure_pa=1e5, density_kg_m3=1000.0)
            base.update(args)
            with self.subTest(args=args), self.assertRaises(DomainError):
                net_positive_suction_head(**base)

    def test_margin_decision_is_inclusive_and_ratio_scaled(self):
        density = 1000.0
        npsha = 2e5 / (density * STANDARD_GRAVITY_M_S2)
        exact = evaluate_pump_suction(inlet_total_pressure_pa=3e5, vapor_pressure_pa=1e5,
                                      density_kg_m3=density, npsh_required_m=npsha)
        self.assertTrue(exact.cavitation_margin_adequate)
        self.assertEqual(exact.npsh_margin_m, 0.0)
        strict = evaluate_pump_suction(inlet_total_pressure_pa=3e5, vapor_pressure_pa=1e5,
                                       density_kg_m3=density, npsh_required_m=npsha,
                                       required_margin_ratio=1.25)
        self.assertFalse(strict.cavitation_margin_adequate)
        self.assertIn("insufficient_npsh_margin", [w.code for w in strict.warnings])

    def test_suction_specific_speed_and_input_pairing(self):
        r = evaluate_pump_suction(inlet_total_pressure_pa=5e5, vapor_pressure_pa=1e5,
                                  density_kg_m3=800.0, npsh_required_m=16 * FOOT,
                                  shaft_speed_rad_s=10000 * RPM, flow_m3_s=400 * GPM)
        self.assertAlmostEqual(r.suction_specific_speed, 25000.0, delta=1e-6)
        with self.assertRaises(InputError):
            evaluate_pump_suction(inlet_total_pressure_pa=5e5, vapor_pressure_pa=1e5,
                                  density_kg_m3=800.0, npsh_required_m=5.0, shaft_speed_rad_s=100.0)
        with self.assertRaises(DomainError):
            evaluate_pump_suction(inlet_total_pressure_pa=5e5, vapor_pressure_pa=1e5,
                                  density_kg_m3=800.0, npsh_required_m=5.0,
                                  required_margin_ratio=0.9)


class PumpMap(unittest.TestCase):
    def test_design_point_and_homologous_identity(self):
        c = characteristic()
        self.assertAlmostEqual(c.head_m(0.1, 3000.0), 1000.0, delta=1e-9)
        for n in (0.6, 1.0, 1.4):
            for q in (0.0, 0.05, 0.12):
                self.assertAlmostEqual(c.head_m(q * n, 3000.0 * n), n * n * c.head_m(q, 3000.0),
                                       delta=1e-9)

    def test_slope_zero_slope_flow_and_sp8109_rule(self):
        c = characteristic()
        point = evaluate_pump_map(c, flow_m3_s=0.1, shaft_speed_rad_s=3000.0)
        self.assertAlmostEqual(point.head_slope_s_m2, 1e4 * (0.1 - 0.6))
        self.assertTrue(point.negative_slope)
        self.assertAlmostEqual(point.zero_slope_flow_m3_s, 0.1 * 0.1 / 0.6)
        self.assertTrue(point.meets_sp8109_zero_slope_rule)
        near = evaluate_pump_map(c, flow_m3_s=0.017, shaft_speed_rad_s=3000.0)
        self.assertFalse(near.meets_sp8109_zero_slope_rule)
        self.assertIn("pump_zero_slope_margin_below_sp8109", [w.code for w in near.warnings])
        rising = evaluate_pump_map(characteristic(linear_coefficient=0.4, quadratic_coefficient=-0.6),
                                   flow_m3_s=0.01, shaft_speed_rad_s=3000.0)
        self.assertFalse(rising.negative_slope)
        self.assertIn("pump_head_slope_not_negative", [w.code for w in rising.warnings])
        convex = evaluate_pump_map(characteristic(linear_coefficient=-0.5, quadratic_coefficient=0.3),
                                   flow_m3_s=0.1, shaft_speed_rad_s=3000.0)
        self.assertIsNone(convex.meets_sp8109_zero_slope_rule)

    def test_runout_beyond_positive_head_is_critical(self):
        point = evaluate_pump_map(characteristic(), flow_m3_s=0.3, shaft_speed_rad_s=3000.0)
        self.assertLessEqual(point.head_m, 0)
        self.assertIsNone(point.specific_speed)
        self.assertIn("pump_map_nonpositive_head", [w.code for w in point.warnings])

    def test_characteristic_validation(self):
        with self.assertRaises(DomainError) as caught:
            characteristic(shutoff_coefficient=1.0)
        self.assertEqual(caught.exception.code, "characteristic_not_normalized")
        for changes in ({"shutoff_coefficient": 0.0, "linear_coefficient": 1.3},
                        {"design_head_m": 0.0}, {"quadratic_coefficient": float("nan")}):
            with self.subTest(changes=changes), self.assertRaises(DomainError):
                characteristic(**changes)
        with self.assertRaises(DomainError):
            characteristic().head_m(-0.1, 3000.0)
        with self.assertRaises(DomainError):
            evaluate_pump_map(object(), flow_m3_s=0.1, shaft_speed_rad_s=3000.0)


class CyclePoint(unittest.TestCase):
    def cycle(self, **changes):
        values = dict(cycle="gas-generator", chamber_pressure_pa=7e6, chamber_mass_flow_kg_s=25.0,
                      propellant_density_kg_m3=1030.0, inlet_pressure_pa=2.5e5,
                      injector_pressure_drop_fraction=0.2, pump_efficiency=0.68,
                      turbine_efficiency=0.65, turbine_inlet_temperature_k=900.0,
                      turbine_specific_heat_j_kg_k=2200.0, turbine_gamma=1.25,
                      turbine_pressure_ratio=8.0)
        values.update(changes)
        return calculate_engine_cycle(**values)

    def test_cycle_duty_is_converted_and_compared(self):
        cycle = self.cycle()
        required = cycle.pump_pressure_rise_pa / (1030.0 * STANDARD_GRAVITY_M_S2)
        flow = cycle.pump_mass_flow_kg_s / 1030.0
        pump = PumpCharacteristic(design_flow_m3_s=flow, design_head_m=required * 1.1,
                                  design_speed_rad_s=3000.0, shutoff_coefficient=1.2,
                                  linear_coefficient=0.1, quadratic_coefficient=-0.3)
        ok = evaluate_cycle_pump_point(cycle, pump, shaft_speed_rad_s=3000.0)
        self.assertEqual((ok.required_head_m, ok.volumetric_flow_m3_s), (required, flow))
        self.assertTrue(ok.head_sufficient)
        self.assertAlmostEqual(ok.head_margin_m, 0.1 * required, delta=1e-9 * required)
        slow = evaluate_cycle_pump_point(cycle, pump, shaft_speed_rad_s=2500.0)
        self.assertFalse(slow.head_sufficient)
        self.assertIn("pump_head_below_cycle_requirement", [w.code for w in slow.warnings])
        self.assertEqual(self.cycle(), cycle)

    def test_pressure_fed_cycle_has_no_pump(self):
        fed = calculate_engine_cycle(cycle="pressure-fed", chamber_pressure_pa=1e6,
                                     chamber_mass_flow_kg_s=1.0, inlet_pressure_pa=2e6,
                                     injector_pressure_drop_fraction=0.2)
        with self.assertRaises(InputError) as caught:
            evaluate_cycle_pump_point(fed, characteristic(), shaft_speed_rad_s=3000.0)
        self.assertEqual(caught.exception.code, "cycle_has_no_pump")
        with self.assertRaises(InputError):
            evaluate_cycle_pump_point(object(), characteristic(), shaft_speed_rad_s=3000.0)


if __name__ == "__main__":
    unittest.main()
