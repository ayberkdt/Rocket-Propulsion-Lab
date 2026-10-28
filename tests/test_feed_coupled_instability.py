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
from math import atan, ceil, floor, pi, sqrt, tan

import mpmath as mp

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


def stability(tau: float, ox_drop: float, points: int = 40000) -> fc.FeedCoupledStability:
    g_ox, g_f = rigid(ox_drop)
    return fc.feed_coupled_stability(
        mixture_ratio=R, characteristic_velocity_m_s=CSTAR, cstar_mixture_ratio_derivative_m_s=SLOPE,
        chamber_pressure_pa=PC, gas_residence_time_s=THETA, combustion_time_lag_s=tau,
        oxidizer_admittance=g_ox, fuel_admittance=g_f,
        oxidizer_admittance_bounds=fc.AdmittanceBounds(g_ox(0), 0.0),
        fuel_admittance_bounds=fc.AdmittanceBounds(g_f(0), 0.0), points=points)


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


class NumericalResolution(unittest.TestCase):
    def test_unstable_pair_below_the_old_first_frequency(self) -> None:
        base_ox, base_fuel = rigid(0.3)
        scale = 1.0001 / gain(0.3)
        tau = 0.36
        # Independent characteristic root, not another phase-grid implementation.
        actual_gain = ((1 + (1 + R) / CSTAR * SLOPE) * base_ox(0) * scale
                       + (1 - R * (1 + R) / CSTAR * SLOPE) * base_fuel(0) * scale)
        with mp.workdps(80):
            k, delay, theta = mp.mpf(actual_gain), mp.mpf(tau), mp.mpf(THETA)
            root = mp.lambertw(-k * delay / theta * mp.exp(delay / theta), 0) / delay - 1 / theta
            self.assertGreater(mp.re(root), 0)
            self.assertLess(abs(1 + theta * root + k * mp.exp(-delay * root)), mp.mpf('1e-65'))
        result = fc.feed_coupled_stability(
            mixture_ratio=R, characteristic_velocity_m_s=CSTAR, cstar_mixture_ratio_derivative_m_s=SLOPE,
            chamber_pressure_pa=PC, gas_residence_time_s=THETA, combustion_time_lag_s=tau,
            oxidizer_admittance=lambda s: base_ox(s) * scale, fuel_admittance=lambda s: base_fuel(s) * scale,
            oxidizer_admittance_bounds=fc.AdmittanceBounds(base_ox(0) * scale, 0.0),
            fuel_admittance_bounds=fc.AdmittanceBounds(base_fuel(0) * scale, 0.0))
        self.assertFalse(result.stable)
        self.assertEqual(result.right_half_plane_zeros, 2)

    def test_long_lag_counts_against_all_analytic_crossings(self) -> None:
        k = gain(0.3)
        omega = sqrt(k * k - 1) / THETA
        for tau in (0.1, 0.5, 1.0, 5.0, 10.0):
            expected = 2 * max(0, ceil((tau * omega + atan(THETA * omega) - pi) / (2 * pi)))
            with self.subTest(tau=tau):
                self.assertEqual(stability(tau, 0.3, points=200000).right_half_plane_zeros, expected)

    def test_insufficient_budget_does_not_return_a_verdict(self) -> None:
        with self.assertRaises(DomainError) as raised:
            stability(10.0, 0.3, points=1000)
        self.assertEqual(raised.exception.code, "nyquist_resolution")

    def test_neutral_mode_does_not_get_rounded_to_stable(self) -> None:
        k = gain(0.3)
        omega = sqrt(k * k - 1) / THETA
        tau = (pi - atan(THETA * omega)) / omega
        with self.assertRaises(DomainError) as raised:
            stability(tau, 0.3)
        self.assertIn(raised.exception.code, {"nyquist_resolution", "nyquist_boundary_unresolved"})


class BoundContracts(unittest.TestCase):
    def solve(self, response, bounds, **changes):
        p = {"mixture_ratio": 1.0, "characteristic_velocity_m_s": 1.0,
             "cstar_mixture_ratio_derivative_m_s": 0.0, "chamber_pressure_pa": PC,
             "gas_residence_time_s": 0.001, "combustion_time_lag_s": 0.0,
             "oxidizer_admittance": response, "fuel_admittance": lambda s: 0.0,
             "oxidizer_admittance_bounds": bounds, "fuel_admittance_bounds": fc.AdmittanceBounds(0.0, 0.0)}
        return fc.feed_coupled_stability(**{**p, **changes})

    def test_callbacks_require_bounds_and_real_finite_response(self) -> None:
        for response, bounds, code in ((lambda s: 0.5, None, "admittance_bounds_required"),
                                       (lambda s: 2.0, fc.AdmittanceBounds(1.0, 0.0), "nyquist_invalid_bound"),
                                       (lambda s: complex('nan'), fc.AdmittanceBounds(1.0, 0.0), "nyquist_response"),
                                       (lambda s: 1j, fc.AdmittanceBounds(1.0, 0.0), "nyquist_response")):
            with self.subTest(code=code), self.assertRaises(DomainError) as raised:
                self.solve(response, bounds)
            self.assertEqual(raised.exception.code, code)

    def test_understated_slope_is_refused(self) -> None:
        with self.assertRaises(DomainError) as raised:
            self.solve(lambda s: 2 / (1 + 100 * s), fc.AdmittanceBounds(2.0, 0.0))
        self.assertEqual(raised.exception.code, "nyquist_invalid_bound")

    def test_small_terminal_value_is_not_a_tail_bound(self) -> None:
        omega, damping, amplitude = 1000.0, 0.01, 10.0

        def resonance(s):
            return amplitude * 2 * damping * omega * s / (s * s + 2 * damping * omega * s + omega * omega)

        self.assertLess(abs(resonance(2j * pi * 10) / (1 + 2j * pi * 10 * 0.001)), 0.05)
        self.assertGreater(abs(resonance(1j * omega) / (1 + 1j * omega * 0.001)), 1)
        with self.assertRaises(DomainError) as raised:
            self.solve(resonance, fc.AdmittanceBounds(amplitude, 1e6), maximum_frequency_hz=10.0)
        self.assertEqual(raised.exception.code, "nyquist_tail_unresolved")

    def test_negative_dc_gain_can_have_an_odd_real_root(self) -> None:
        # F(s)=1-1.5/(1+0.001s), so the only closed-loop root is +500/s.
        result = self.solve(lambda s: 0.5, fc.AdmittanceBounds(0.5, 0.0),
                            cstar_mixture_ratio_derivative_m_s=-2.0)
        self.assertFalse(result.stable)
        self.assertEqual(result.right_half_plane_zeros, 1)

    def test_small_gain_margin_is_a_bound_over_all_frequencies(self) -> None:
        result = self.solve(lambda s: 0.5, fc.AdmittanceBounds(0.5, 0.0))
        self.assertTrue(result.stable)
        self.assertEqual(result.minimum_distance_to_minus_one, 0.5)

    def test_signed_gains_and_time_scale_invariance(self) -> None:
        for k in (-2.0, -1.01, -0.5, 0.5, 1.0, 1.0001, 1.372, 3.0):
            for theta in (1e-6, 0.0007, 0.3):
                for lag_ratio in (0.2, 2.0, 50.0):
                    expected = 0
                    if abs(k) > 1:
                        frequency = sqrt(k * k - 1)
                        angle = lag_ratio * frequency + atan(frequency)
                        expected = (2 * max(0, ceil((angle - pi) / (2 * pi))) if k > 0
                                    else 1 + 2 * floor(angle / (2 * pi)))
                    with self.subTest(gain=k, theta=theta, lag_ratio=lag_ratio):
                        result = self.solve(lambda s: 1.0, fc.AdmittanceBounds(1.0, 0.0),
                                            cstar_mixture_ratio_derivative_m_s=(k - 1) / 2,
                                            gas_residence_time_s=theta, combustion_time_lag_s=lag_ratio * theta)
                        self.assertEqual(result.right_half_plane_zeros, expected)
                        self.assertEqual(result.stable, expected == 0)
                        self.assertGreater(result.minimum_distance_to_minus_one, 0)


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

    def test_finite_tank_zero_length_line_against_polynomial_roots(self) -> None:
        # G_ox=G_f=0.01*s/(1+0.01*s), including the finite tank DC limit.
        # With zero combustion lag the closed-loop equation is quadratic.
        common = {"line_area_m2": 1e-3, "liquid_sound_speed_m_s": 1000.0, "liquid_density_kg_m3": 1000.0,
                  "liquid_bulk_modulus_pa": 1e9, "manifold_volume_m3": 0.0, "line_length_m": 0.0,
                  "tank_volume_m3": 0.03}
        ox = fc.FeedSystem(propellant_mass_flow_kg_s=2.0, injector_pressure_drop_pa=PC / 3, **common)
        fu = fc.FeedSystem(propellant_mass_flow_kg_s=1.0, injector_pressure_drop_pa=PC / 6, **common)
        for slope in (1.0, 1.2):
            with self.subTest(slope=slope), mp.workdps(60):
                roots = mp.polyroots([mp.mpf('0.00001'), mp.mpf('0.011') + (2 - 3 * mp.mpf(slope)) * mp.mpf('0.01'), 1])
                expected = sum(mp.re(root) > 0 for root in roots)
                result = fc.feed_coupled_stability(
                    mixture_ratio=2.0, characteristic_velocity_m_s=1.0,
                    cstar_mixture_ratio_derivative_m_s=slope, chamber_pressure_pa=PC,
                    gas_residence_time_s=0.001, combustion_time_lag_s=0.0,
                    oxidizer_feed=ox, fuel_feed=fu)
                self.assertEqual(result.right_half_plane_zeros, expected)

    def test_derived_slope_envelopes_bound_high_precision_derivatives(self) -> None:
        # Independent high-precision differentiation of the original transfer
        # function, including intervals around acoustic resonances.
        p = {"injector_impedance": 0.9, "line_impedance": 0.35, "line_time_s": 0.004,
             "manifold_capacitance_s": 0.0012}
        for tank in (None, 0.05):
            normalized = {**p, "tank_capacitance_s": tank}
            envelope = fc._feed_admittance_bounds(normalized)
            for lower, upper in ((0.0, 1.0), (390.0, 396.0), (780.0, 790.0), (1e3, 1e4)):
                with self.subTest(tank=tank, lower=lower), mp.workdps(60):
                    def g(w, tank=tank):
                        s = mp.j * w
                        zi, zl, tl, cm = map(mp.mpf, (0.9, 0.35, 0.004, 0.0012))
                        th = mp.tanh(s * tl)
                        if tank is None:
                            return (1 + s * cm * zl * th) / (zi * (1 + (zl / zi + s * cm * zl) * th))
                        ct = mp.mpf(tank)
                        num = s * (ct + cm) + (1 / zl + s**2 * cm * ct * zl) * th
                        den = 1 + s * zi * (ct + cm) + (zi / zl + s * ct * zl + s**2 * cm * ct * zl * zi) * th
                        return num / den
                    bound = envelope.slope(lower, upper)
                    for w in (lower, (lower + upper) / 2, upper):
                        self.assertLessEqual(abs(mp.diff(g, mp.mpf(w))), bound)

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
