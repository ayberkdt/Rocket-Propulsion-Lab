"""Tests for rocket performance and propellant reference records."""

import unittest
from math import log

from rocket_propulsion.core.errors import DomainError
from rocket_propulsion.propulsion import (
    STANDARD_GRAVITY_M_S2,
    calculate_rocket_equation,
    calculate_thrust,
    list_propellants,
)


class PropulsionTests(unittest.TestCase):
    def test_tsiolkovsky_forward_and_inverse(self) -> None:
        forward = calculate_rocket_equation(
            initial_mass_kg=500_000.0,
            final_mass_kg=100_000.0,
            specific_impulse_s=450.0,
        )
        self.assertAlmostEqual(
            forward.delta_v_m_s, 450.0 * STANDARD_GRAVITY_M_S2 * log(5.0)
        )
        inverse = calculate_rocket_equation(
            initial_mass_kg=500_000.0,
            delta_v_m_s=forward.delta_v_m_s,
            specific_impulse_s=450.0,
        )
        self.assertAlmostEqual(inverse.final_mass_kg, 100_000.0, places=6)

    def test_thrust_separates_momentum_and_pressure(self) -> None:
        result = calculate_thrust(
            mass_flow_kg_s=10.0,
            exit_velocity_m_s=3000.0,
            exit_area_m2=0.5,
            exit_pressure_pa=50_000.0,
            ambient_pressure_pa=10_000.0,
            burn_time_s=20.0,
        )
        self.assertAlmostEqual(result.momentum_thrust_n, 30_000.0)
        self.assertAlmostEqual(result.pressure_thrust_n, 20_000.0)
        self.assertAlmostEqual(result.total_thrust_n, 50_000.0)
        self.assertAlmostEqual(result.total_impulse_n_s or 0.0, 1_000_000.0)

    def test_catalog_covers_liquid_solid_and_hybrid(self) -> None:
        self.assertEqual(len(list_propellants("liquid")), 4)
        self.assertEqual(len(list_propellants("solid")), 1)
        self.assertEqual(len(list_propellants("hybrid")), 2)

    def test_invalid_catalog_category_is_rejected(self) -> None:
        with self.assertRaises(DomainError):
            list_propellants("plasma")


if __name__ == "__main__":
    unittest.main()

