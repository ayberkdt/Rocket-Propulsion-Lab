"""Regime-boundary, conservation, and altitude tests for C-D nozzles."""

import unittest
from itertools import pairwise

from rocket_propulsion.compressible import (
    calculate_off_design_nozzle,
    generate_nozzle_altitude_sweep,
    nozzle_regime_thresholds,
    standard_atmosphere_pressure_pa,
)


class OffDesignNozzleTests(unittest.TestCase):
    def setUp(self) -> None:
        self.inputs = {
            "chamber_pressure_pa": 7_000_000.0,
            "chamber_temperature_k": 3_500.0,
            "throat_area_m2": 0.01,
            "area_ratio": 20.0,
            "gamma": 1.22,
            "gas_constant_j_kg_k": 355.0,
        }
        self.thresholds = nozzle_regime_thresholds(
            chamber_pressure_pa=self.inputs["chamber_pressure_pa"],
            area_ratio=self.inputs["area_ratio"],
            gamma=self.inputs["gamma"],
        )

    def solve(self, ambient_pressure_pa: float):
        return calculate_off_design_nozzle(
            **self.inputs, ambient_pressure_pa=ambient_pressure_pa
        )

    def test_analytical_threshold_order(self) -> None:
        values = (
            self.thresholds.chamber_pressure_pa,
            self.thresholds.just_choked_back_pressure_pa,
            self.thresholds.shock_at_exit_back_pressure_pa,
            self.thresholds.design_exit_pressure_pa,
            0.0,
        )
        self.assertTrue(all(left > right for left, right in pairwise(values)))

    def test_regime_order_as_back_pressure_falls(self) -> None:
        just = self.thresholds.just_choked_back_pressure_pa
        shock = self.thresholds.shock_at_exit_back_pressure_pa
        design = self.thresholds.design_exit_pressure_pa
        pressures = (
            self.inputs["chamber_pressure_pa"],
            0.5 * (self.inputs["chamber_pressure_pa"] + just),
            just,
            0.5 * (just + shock),
            shock,
            0.5 * (shock + design),
            design,
            0.5 * design,
        )
        self.assertEqual(
            [self.solve(value).regime for value in pressures],
            [
                "no-flow",
                "unchoked",
                "just-choked",
                "internal-normal-shock",
                "shock-at-exit",
                "overexpanded",
                "design",
                "underexpanded",
            ],
        )

    def test_internal_shock_conserves_mass_and_moves_downstream(self) -> None:
        just = self.thresholds.just_choked_back_pressure_pa
        shock_exit = self.thresholds.shock_at_exit_back_pressure_pa
        high_back_pressure = self.solve(0.75 * just + 0.25 * shock_exit)
        low_back_pressure = self.solve(0.25 * just + 0.75 * shock_exit)
        self.assertLess(high_back_pressure.mass_flow_relative_residual, 1.0e-9)
        self.assertLess(low_back_pressure.mass_flow_relative_residual, 1.0e-9)
        self.assertLess(
            high_back_pressure.shock.x_fraction,
            low_back_pressure.shock.x_fraction,
        )
        self.assertLess(high_back_pressure.exit_total_pressure_pa, 7_000_000.0)

    def test_summerfield_warning_is_explicitly_empirical(self) -> None:
        result = self.solve(101_325.0)
        self.assertEqual(result.regime, "overexpanded")
        self.assertTrue(result.separation.predicted)
        self.assertIn("15 degree", result.separation.valid_for)
        self.assertEqual(result.warnings[0].code, "possible_flow_separation")

    def test_standard_atmosphere_and_altitude_sweep(self) -> None:
        self.assertAlmostEqual(standard_atmosphere_pressure_pa(0.0), 101_325.0, places=3)
        self.assertAlmostEqual(
            standard_atmosphere_pressure_pa(11_000.0), 22_632.06, delta=0.1
        )
        sweep = generate_nozzle_altitude_sweep(
            **self.inputs, maximum_altitude_m=50_000.0, point_count=11
        )
        self.assertEqual(len(sweep.points), 11)
        self.assertGreater(sweep.points[-1].thrust_n, sweep.points[0].thrust_n)
        self.assertEqual(sweep.points[-1].regime, "underexpanded")


if __name__ == "__main__":
    unittest.main()

