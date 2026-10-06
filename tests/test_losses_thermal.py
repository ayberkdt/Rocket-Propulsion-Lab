"""Tests for performance loss accounting and preliminary thermal models."""

import unittest

from rocket_propulsion.core.errors import DomainError
from rocket_propulsion.geometry import generate_rao_contour
from rocket_propulsion.propulsion import calculate_loss_budget, calculate_two_phase_band
from rocket_propulsion.thermodynamics import (
    calculate_bartz_heat_transfer,
    calculate_regenerative_cooling,
)


class PerformanceLossTests(unittest.TestCase):
    def test_loss_terms_reconcile_and_can_all_be_disabled(self) -> None:
        result = calculate_loss_budget(
            ideal_thrust_n=100_000.0,
            ideal_specific_impulse_s=450.0,
        )
        expected = 0.97 * 0.98 * 0.985 * 0.98
        self.assertAlmostEqual(result.total_efficiency, expected)
        self.assertAlmostEqual(result.thrust_balance_residual_n, 0.0)
        self.assertAlmostEqual(result.specific_impulse_balance_residual_s, 0.0)
        ideal = calculate_loss_budget(
            ideal_thrust_n=100_000.0,
            ideal_specific_impulse_s=450.0,
            enable_discharge_loss=False,
            enable_divergence_loss=False,
            enable_boundary_layer_loss=False,
            enable_combustion_loss=False,
        )
        self.assertEqual(ideal.delivered_thrust_n, ideal.ideal_thrust_n)
        self.assertTrue(all(term.thrust_loss_n == 0.0 for term in ideal.terms))

    def test_two_phase_band_bounds_perfect_coupling(self) -> None:
        band = calculate_two_phase_band(
            ideal_thrust_n=100_000.0,
            ideal_specific_impulse_s=450.0,
            condensed_mass_fraction=0.2,
            minimum_particle_coupling_efficiency=0.6,
            maximum_particle_coupling_efficiency=0.9,
        )
        self.assertAlmostEqual(band.thrust_range_n[0], 92_000.0)
        self.assertAlmostEqual(band.thrust_range_n[1], 98_000.0)
        self.assertTrue(band.warnings)


class ThermalLoadTests(unittest.TestCase):
    def setUp(self) -> None:
        self.contour = generate_rao_contour(throat_area_m2=0.01, area_ratio=25.0)

    def _heat_result(self):
        return calculate_bartz_heat_transfer(
            contour=self.contour,
            chamber_pressure_pa=5.0e6,
            chamber_temperature_k=3500.0,
            wall_temperature_k=800.0,
            characteristic_velocity_m_s=1700.0,
            dynamic_viscosity_pa_s=7.0e-5,
            specific_heat_j_kg_k=3000.0,
            prandtl_number=0.7,
            gamma=1.22,
            throat_radius_of_curvature_m=0.03,
        )

    def test_bartz_profile_is_positive_and_energy_integral_closes(self) -> None:
        result = self._heat_result()
        self.assertGreater(result.peak_heat_flux_w_m2, 1.0e6)
        self.assertLess(abs(result.peak_location_m), 0.02)
        self.assertGreater(result.total_heat_load_w, 0.0)
        self.assertAlmostEqual(result.energy_balance_residual_w, 0.0)
        self.assertEqual(len(result.stations), len(self.contour.stations))

    def test_bartz_does_not_hide_invalid_wall_temperature(self) -> None:
        with self.assertRaises(DomainError):
            calculate_bartz_heat_transfer(
                contour=self.contour,
                chamber_pressure_pa=5.0e6,
                chamber_temperature_k=3500.0,
                wall_temperature_k=3600.0,
                characteristic_velocity_m_s=1700.0,
                dynamic_viscosity_pa_s=7.0e-5,
                specific_heat_j_kg_k=3000.0,
                prandtl_number=0.7,
                gamma=1.22,
                throat_radius_of_curvature_m=0.03,
            )

    def test_regenerative_energy_balance_and_pressure_loss(self) -> None:
        heat = self._heat_result()
        cooling = calculate_regenerative_cooling(
            absorbed_heat_w=heat.total_heat_load_w,
            coolant_mass_flow_kg_s=5.0,
            coolant_specific_heat_j_kg_k=10_000.0,
            coolant_inlet_temperature_k=100.0,
            coolant_density_kg_m3=70.0,
            coolant_dynamic_viscosity_pa_s=1.0e-5,
            channel_count=100,
            channel_flow_area_m2=2.0e-5,
            hydraulic_diameter_m=0.004,
            channel_length_m=1.0,
            roughness_m=1.0e-5,
            maximum_coolant_temperature_k=1000.0,
        )
        self.assertAlmostEqual(cooling.energy_balance_residual_w, 0.0)
        self.assertGreater(cooling.outlet_temperature_k, cooling.inlet_temperature_k)
        self.assertGreater(cooling.pressure_drop_pa, 0.0)


if __name__ == "__main__":
    unittest.main()
