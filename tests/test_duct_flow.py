"""Tests for Fanno and Rayleigh line relations."""

import unittest

from rocket_propulsion.compressible.duct_flow import calculate_fanno, calculate_rayleigh


class DuctFlowTests(unittest.TestCase):
    def test_sonic_fanno_state_is_starred_state(self) -> None:
        result = calculate_fanno(1.0, 1.4)
        self.assertAlmostEqual(result.temperature_critical_ratio, 1.0)
        self.assertAlmostEqual(result.pressure_critical_ratio, 1.0)
        self.assertAlmostEqual(result.stagnation_pressure_critical_ratio, 1.0)
        self.assertAlmostEqual(result.friction_parameter_to_sonic, 0.0)
        self.assertAlmostEqual(result.entropy_difference_over_r, 0.0)

    def test_friction_drives_both_branches_toward_choking(self) -> None:
        self.assertGreater(calculate_fanno(0.5).friction_parameter_to_sonic, 0.0)
        self.assertGreater(calculate_fanno(2.0).friction_parameter_to_sonic, 0.0)

    def test_sonic_rayleigh_state_is_starred_state(self) -> None:
        result = calculate_rayleigh(1.0, 1.4)
        self.assertAlmostEqual(result.stagnation_temperature_critical_ratio, 1.0)
        self.assertAlmostEqual(result.temperature_critical_ratio, 1.0)
        self.assertAlmostEqual(result.pressure_critical_ratio, 1.0)
        self.assertAlmostEqual(result.stagnation_pressure_critical_ratio, 1.0)
        self.assertAlmostEqual(result.entropy_difference_over_r, 0.0)


if __name__ == "__main__":
    unittest.main()

