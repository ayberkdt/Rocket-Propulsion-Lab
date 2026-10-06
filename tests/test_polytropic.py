"""Tests for polytropic paths and first-law energy accounting."""

import unittest

from rocket_propulsion.thermodynamics.polytropic import calculate_polytropic_process


class PolytropicTests(unittest.TestCase):
    def test_gamma_exponent_is_reversible_adiabatic(self) -> None:
        result = calculate_polytropic_process(
            inlet_pressure_pa=100_000.0,
            outlet_pressure_pa=500_000.0,
            inlet_temperature_k=300.0,
            exponent=1.4,
            gamma=1.4,
        )
        self.assertAlmostEqual(result.heat_into_gas_j_kg, 0.0, places=7)
        self.assertAlmostEqual(result.entropy_change_j_kg_k, 0.0, places=9)
        self.assertLess(result.boundary_work_by_gas_j_kg, 0.0)
        self.assertEqual(len(result.path), 61)

    def test_isothermal_compression_rejects_heat(self) -> None:
        result = calculate_polytropic_process(
            inlet_pressure_pa=100_000.0,
            outlet_pressure_pa=200_000.0,
            inlet_temperature_k=300.0,
            exponent=1.0,
        )
        self.assertAlmostEqual(result.outlet_temperature_k, 300.0)
        self.assertAlmostEqual(result.internal_energy_change_j_kg, 0.0)
        self.assertAlmostEqual(result.heat_into_gas_j_kg, result.boundary_work_by_gas_j_kg)
        self.assertLess(result.entropy_change_j_kg_k, 0.0)


if __name__ == "__main__":
    unittest.main()

