"""Tests for ideal choked-nozzle performance."""

import unittest

from rocket_propulsion.compressible.nozzle import calculate_ideal_nozzle
from rocket_propulsion.core.errors import DomainError


class NozzleTests(unittest.TestCase):
    def test_nozzle_state_is_internally_consistent(self) -> None:
        result = calculate_ideal_nozzle(
            chamber_pressure_pa=7_000_000.0,
            chamber_temperature_k=3500.0,
            throat_area_m2=0.01,
            area_ratio=20.0,
            ambient_pressure_pa=101_325.0,
            gamma=1.22,
            gas_constant_j_kg_k=355.0,
        )

        self.assertGreater(result.exit_mach, 3.0)
        self.assertGreater(result.mass_flow_kg_s, 0.0)
        self.assertAlmostEqual(result.exit_area_m2, 0.2)
        self.assertAlmostEqual(
            result.thrust_n, result.momentum_thrust_n + result.pressure_thrust_n
        )
        self.assertAlmostEqual(
            result.characteristic_velocity_m_s,
            7_000_000.0 * 0.01 / result.mass_flow_kg_s,
        )

    def test_rejects_nonphysical_area_ratio(self) -> None:
        with self.assertRaises(DomainError):
            calculate_ideal_nozzle(
                chamber_pressure_pa=1_000_000.0,
                chamber_temperature_k=3000.0,
                throat_area_m2=0.01,
                area_ratio=0.9,
                ambient_pressure_pa=0.0,
            )


if __name__ == "__main__":
    unittest.main()

