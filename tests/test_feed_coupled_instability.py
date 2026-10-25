"""Feed-coupled chugging: SP-194 Sect. 5.4 admittances and the Nyquist count of Sect. 5.5.1.

Checks:

* the closed-form injection admittance (Eq. 5.4.2-1) equals the composition of
  the component ratios (tank capacitance, Eq. 5.4.1-2 line, Eq. 5.4.1-16
  manifold, Eq. 5.4.1-12 orifice) started from ``G_0 = 0``, and tends to Eq.
  (5.4.2-2) as the tank capacitance grows;
* at ``s = i omega``, Eq. (5.4.2-2) equals the printed real-frequency form
  Eq. (5.4.2-3); without a manifold the admittance vanishes at the line's
  quarter-wave resonance and returns to ``1/Z_i`` at its half-wave
  resonance, as SP-194 p. 244 states;
* rigid feeds reproduce Example 5.5.1.1: no right-half-plane zero at
  ``K = 0.894``; two (a complex pair, as the book says) at ``K = 1.372`` with
  a 2.75 ms lag; the count changes across the printed 1.78 ms limit, with the
  gain crossover at the neutral frequency of Eq. (5.5.1-6);
* a physical feed whose line is short and manifold empty reduces to the
  rigid result; Eq. (5.4.1-5) reduces to the unconfined sound speed for a
  rigid wall.
"""

from __future__ import annotations

import unittest
from math import atan, pi, sqrt, tan

from rocket_propulsion.core.errors import DomainError
from rocket_propulsion.thermochemistry import feed_coupled_instability as fc

FT = 0.3048
R = 5.2
CSTAR = 7920 * FT
SLOPE = -185 * FT
THETA = 0.7e-3
PC = 1.0e6


def rigid(ox_drop: float, fuel_drop: float = 0.8):
    return (lambda s: (R / (1 + R)) / (2 * ox_drop)), (lambda s: (1 / (1 + R)) / (2 * fuel_drop))


def gain(ox_drop: float) -> float:
    return ((1 + (1 + R) / CSTAR * SLOPE) * R / (1 + R) / (2 * ox_drop)
            + (1 - R * (1 + R) / CSTAR * SLOPE) / (1 + R) / (2 * 0.8))


def stability(tau: float, ox_drop: float) -> fc.FeedCoupledStability:
    g_ox, g_f = rigid(ox_drop)
    return fc.feed_coupled_stability(
        mixture_ratio=R, characteristic_velocity_m_s=CSTAR, cstar_mixture_ratio_derivative_m_s=SLOPE,
        chamber_pressure_pa=PC, gas_residence_time_s=THETA, combustion_time_lag_s=tau,
        oxidizer_admittance=g_ox, fuel_admittance=g_f)


class Admittances(unittest.TestCase):
    ZI, ZL, TL, CM, CT = 0.9, 0.35, 0.004, 0.0012, 0.05

    def chain(self, s: complex, tank: float | None) -> complex:
        g = None if tank is None else fc.manifold_admittance_ratio(0j, s, manifold_capacitance_s=tank)
        g = fc.line_admittance_ratio(g, s, line_impedance=self.ZL, line_time_s=self.TL)
        g = fc.manifold_admittance_ratio(g, s, manifold_capacitance_s=self.CM)
        return fc.orifice_admittance_ratio(g, s, orifice_impedance=self.ZI)

    def closed(self, s: complex, tank: float | None) -> complex:
        return fc.pressurized_tank_injection_admittance(
            s, injector_impedance=self.ZI, line_impedance=self.ZL, line_time_s=self.TL,
            manifold_capacitance_s=self.CM, tank_capacitance_s=tank)

    def test_closed_form_equals_component_chain(self) -> None:
        for s in (2j * pi * 37.0, 2j * pi * 137.0 + 30.0, 2j * pi * 900.0 - 12.0):
            for tank in (0.05, 3.0, None):
                with self.subTest(s=s, tank=tank):
                    self.assertLess(abs(self.chain(s, tank) - self.closed(s, tank)), 1e-12)

    def test_large_tank_tends_to_eq_5_4_2_2(self) -> None:
        s = 2j * pi * 250.0
        self.assertLess(abs(self.closed(s, 1e9) - self.closed(s, None)), 1e-9)

    def test_real_frequency_form_eq_5_4_2_3(self) -> None:
        for f in (20.0, 75.0, 310.0):
            w = 2 * pi * f
            t = tan(w * self.TL)
            printed = (1 / self.ZI) * (1 - w * self.CM * self.ZL * t) / (
                1 - w * self.CM * self.ZL * t + 1j * (self.ZL / self.ZI) * t)
            self.assertLess(abs(self.closed(1j * w, None) - printed), 1e-12)

    def test_quarter_and_half_wave_without_manifold(self) -> None:
        def g(f: float) -> complex:
            return fc.pressurized_tank_injection_admittance(
                2j * pi * f, injector_impedance=self.ZI, line_impedance=self.ZL,
                line_time_s=self.TL, manifold_capacitance_s=0.0)
        quarter = 1 / (4 * self.TL)
        self.assertLess(abs(g(quarter)), 1e-10)
        self.assertLess(abs(g(2 * quarter) - 1 / self.ZI), 1e-10)

    def test_wall_compliance(self) -> None:
        stiff = fc.wall_compliant_sound_speed(liquid_sound_speed_m_s=1000.0, liquid_density_kg_m3=1100.0,
                                              line_diameter_m=0.05, wall_thickness_m=0.002,
                                              wall_elastic_modulus_pa=1e30)
        self.assertLess(abs(stiff - 1000.0), 1e-9)
        soft = fc.wall_compliant_sound_speed(liquid_sound_speed_m_s=1000.0, liquid_density_kg_m3=1100.0,
                                             line_diameter_m=0.05, wall_thickness_m=0.002,
                                             wall_elastic_modulus_pa=2e11)
        self.assertLess(abs(soft - 1000.0 / sqrt(1 + 1100 * 1e6 * 0.05 / (2e11 * 0.002))), 1e-9)


class ExampleNyquist(unittest.TestCase):
    def test_case_a_stable(self) -> None:
        result = stability(2.75e-3, 0.5)
        self.assertTrue(result.stable)
        self.assertEqual(result.right_half_plane_zeros, 0)

    def test_case_b_two_zeros(self) -> None:
        result = stability(2.75e-3, 0.3)
        self.assertFalse(result.stable)
        self.assertEqual(result.right_half_plane_zeros, 2)  # "two zeros (complex conjugates)", p. 253
        k = gain(0.3)
        neutral = sqrt(k * k - 1) / (2 * pi * THETA)
        self.assertLess(abs(result.gain_crossover_hz[0] / neutral - 1.0), 2e-3)
        self.assertLess(result.phase_margins_deg[0], 0.0)

    def test_count_changes_at_the_printed_limit(self) -> None:
        k = gain(0.3)
        limit = THETA * (pi - atan(sqrt(k * k - 1))) / sqrt(k * k - 1)
        self.assertLess(abs(limit - 1.78e-3), 5e-6)
        self.assertTrue(stability(0.99 * limit, 0.3).stable)
        self.assertFalse(stability(1.01 * limit, 0.3).stable)


class PhysicalFeeds(unittest.TestCase):
    def feeds(self, line_length: float, manifold: float, ox_drop: float) -> tuple[fc.FeedSystem, fc.FeedSystem]:
        m_total = 10.0
        m_ox, m_f = m_total * R / (1 + R), m_total / (1 + R)
        common = {"line_area_m2": 2e-3, "liquid_sound_speed_m_s": 900.0, "liquid_density_kg_m3": 1000.0,
                  "liquid_bulk_modulus_pa": 8e8, "manifold_volume_m3": manifold, "line_length_m": line_length}
        return (fc.FeedSystem(propellant_mass_flow_kg_s=m_ox, injector_pressure_drop_pa=ox_drop * PC, **common),
                fc.FeedSystem(propellant_mass_flow_kg_s=m_f, injector_pressure_drop_pa=0.8 * PC, **common))

    def solve(self, tau: float, ox_drop: float, line_length: float, manifold: float) -> fc.FeedCoupledStability:
        ox, fu = self.feeds(line_length, manifold, ox_drop)
        return fc.feed_coupled_stability(
            mixture_ratio=R, characteristic_velocity_m_s=CSTAR, cstar_mixture_ratio_derivative_m_s=SLOPE,
            chamber_pressure_pa=PC, gas_residence_time_s=THETA, combustion_time_lag_s=tau,
            oxidizer_feed=ox, fuel_feed=fu)

    def test_short_line_reduces_to_rigid(self) -> None:
        ox, _ = self.feeds(1e-9, 0.0, 0.3)
        g = ox.admittance(2j * pi * 200.0, total_mass_flow_kg_s=10.0, chamber_pressure_pa=PC)
        self.assertLess(abs(g - (R / (1 + R)) / 0.6), 1e-6)
        self.assertEqual(self.solve(2.75e-3, 0.3, 1e-9, 0.0).right_half_plane_zeros,
                         stability(2.75e-3, 0.3).right_half_plane_zeros)
        self.assertTrue(self.solve(2.75e-3, 0.5, 1e-9, 0.0).stable)

    def test_feed_sweep_matches_admittance(self) -> None:
        ox, _ = self.feeds(2.0, 1e-3, 0.5)
        sweep = fc.feed_admittance_sweep(ox, total_mass_flow_kg_s=10.0, chamber_pressure_pa=PC,
                                         frequencies_hz=(50.0, 112.5))
        self.assertEqual(sweep[1].admittance,
                         ox.admittance(2j * pi * 112.5, total_mass_flow_kg_s=10.0, chamber_pressure_pa=PC))

    def test_inconsistent_inputs_are_refused(self) -> None:
        ox, fu = self.feeds(2.0, 1e-3, 0.5)
        with self.assertRaises(DomainError):
            fc.feed_coupled_stability(mixture_ratio=4.0, characteristic_velocity_m_s=CSTAR,
                                      cstar_mixture_ratio_derivative_m_s=SLOPE, chamber_pressure_pa=PC,
                                      gas_residence_time_s=THETA, combustion_time_lag_s=1e-3,
                                      oxidizer_feed=ox, fuel_feed=fu)
        with self.assertRaises(DomainError):
            fc.feed_coupled_stability(mixture_ratio=R, characteristic_velocity_m_s=CSTAR,
                                      cstar_mixture_ratio_derivative_m_s=SLOPE, chamber_pressure_pa=PC,
                                      gas_residence_time_s=THETA, combustion_time_lag_s=1e-3,
                                      oxidizer_feed=ox)


if __name__ == "__main__":
    unittest.main()
