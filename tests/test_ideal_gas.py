"""Tests for the calorically perfect ideal-gas state."""

import unittest

from rocket_propulsion.core.errors import DomainError
from rocket_propulsion.thermodynamics.ideal_gas import complete_ideal_gas_state


class IdealGasTests(unittest.TestCase):
    def test_completes_density_and_derived_properties(self) -> None:
        state = complete_ideal_gas_state(pressure_pa=101_325.0, temperature_k=288.15)

        self.assertAlmostEqual(state.density_kg_m3, 1.225012, places=5)
        self.assertAlmostEqual(state.cp_j_kg_k, 1004.675, places=6)
        self.assertAlmostEqual(state.cv_j_kg_k, 717.625, places=6)
        self.assertAlmostEqual(state.speed_of_sound_m_s, 340.292, places=3)

    def test_each_state_variable_can_be_missing(self) -> None:
        reference = complete_ideal_gas_state(pressure_pa=200_000.0, temperature_k=500.0)
        pressure = complete_ideal_gas_state(
            temperature_k=reference.temperature_k,
            density_kg_m3=reference.density_kg_m3,
        )
        temperature = complete_ideal_gas_state(
            pressure_pa=reference.pressure_pa,
            density_kg_m3=reference.density_kg_m3,
        )
        self.assertAlmostEqual(pressure.pressure_pa, reference.pressure_pa)
        self.assertAlmostEqual(temperature.temperature_k, reference.temperature_k)

    def test_requires_exactly_two_state_values(self) -> None:
        with self.assertRaises(DomainError):
            complete_ideal_gas_state(pressure_pa=101_325.0)
        with self.assertRaises(DomainError):
            complete_ideal_gas_state(
                pressure_pa=101_325.0,
                temperature_k=288.15,
                density_kg_m3=1.225,
            )


if __name__ == "__main__":
    unittest.main()

