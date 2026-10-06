"""Tests connecting static, total, and process thermodynamic states."""

import unittest

from rocket_propulsion.thermodynamics.flow_state import (
    calculate_isentropic_process,
    calculate_stagnation_state,
)


class FlowThermodynamicsTests(unittest.TestCase):
    def test_stagnation_state_at_mach_two(self) -> None:
        result = calculate_stagnation_state(
            mach=2.0,
            static_pressure_pa=100_000.0,
            static_temperature_k=300.0,
        )
        self.assertAlmostEqual(result.stagnation_temperature_k, 540.0)
        self.assertAlmostEqual(result.stagnation_pressure_pa, 782_444.9067, places=3)
        self.assertAlmostEqual(result.dynamic_pressure_pa, 280_000.0, places=6)
        self.assertAlmostEqual(
            result.stagnation_enthalpy_j_kg - result.static_enthalpy_j_kg,
            0.5 * result.velocity_m_s**2,
            places=6,
        )

    def test_inefficient_compression_increases_entropy(self) -> None:
        result = calculate_isentropic_process(
            inlet_pressure_pa=100_000.0,
            outlet_pressure_pa=400_000.0,
            inlet_temperature_k=300.0,
            efficiency=0.8,
        )
        self.assertEqual(result.process, "compression")
        self.assertGreater(result.actual_outlet_temperature_k, result.ideal_outlet_temperature_k)
        self.assertGreater(result.entropy_change_j_kg_k, 0.0)
        self.assertGreater(result.specific_work_j_kg, 0.0)

    def test_inefficient_expansion_reduces_work_output(self) -> None:
        result = calculate_isentropic_process(
            inlet_pressure_pa=400_000.0,
            outlet_pressure_pa=100_000.0,
            inlet_temperature_k=800.0,
            efficiency=0.85,
        )
        self.assertEqual(result.process, "expansion")
        self.assertLess(result.specific_work_j_kg, 0.0)
        self.assertGreater(result.entropy_change_j_kg_k, 0.0)


if __name__ == "__main__":
    unittest.main()

