"""Combustion-instability relations, checked against NASA SP-194.

Reference values are read from
``tests/reference_data/combustion_instability.json``: Table 3.5.1 (transverse
eigenvalues), Example 5.5.1.1 (single-time-lag chugging limit), Fig. 4.2.2a
(neutral n-tau curve), the Q-factor example of p. 524 and the end-correction
constants of Eqs. (8.3.2-1)-(8.3.2-2).

Analytic checks: each eigenvalue is a root of ``J_nu'``; the chugging and
Crocco-Cheng limits satisfy their characteristic equations at the neutral
frequency; mode frequencies follow the dimensional form of Eq. (3.5.1-19).
"""

from __future__ import annotations

import cmath
import json
import unittest
from math import cos, pi, sqrt
from pathlib import Path
from typing import ClassVar

from rocket_propulsion.core.errors import DomainError
from rocket_propulsion.heat_transfer.conduction import bessel_jn
from rocket_propulsion.thermochemistry import combustion_instability as ci

DATA = json.loads((Path(__file__).parent / "reference_data" / "combustion_instability.json")
                  .read_text(encoding="utf-8"))
FT = 0.3048


def jn_prime(n: int, x: float) -> float:
    return -bessel_jn(1, x) if n == 0 else 0.5 * (bessel_jn(n - 1, x) - bessel_jn(n + 1, x))


class TransverseEigenvalues(unittest.TestCase):
    def test_table_3_5_1(self) -> None:
        table = DATA["table_3_5_1"]
        misprint = tuple(table["misprint"]["mode"])
        for nu, eta, printed in table["rows"]:
            s = ci.transverse_eigenvalue(nu, eta)
            with self.subTest(nu=nu, eta=eta):
                self.assertLess(abs(jn_prime(nu, s)), 1e-12)
                expected = table["misprint"]["root"] if (nu, eta) == misprint else printed
                self.assertLess(abs(s - expected), table["printed_abs_tol"])

    def test_misprint_is_a_misprint(self) -> None:
        """8.5263 is not a root of J_1'; 8.5363 is (and lies between 5.33 and 11.7)."""
        self.assertGreater(abs(jn_prime(1, 8.5263)), 1e-3)
        self.assertLess(abs(ci.transverse_eigenvalue(1, 3) - 8.5363), 1e-4)

    def test_zero_root_counts_for_axisymmetric_modes(self) -> None:
        self.assertEqual(ci.transverse_eigenvalue(0, 1), 0.0)
        with self.assertRaises(DomainError):
            ci.transverse_eigenvalue(-1, 1)


class ModeFrequencies(unittest.TestCase):
    KW: ClassVar[dict[str, float]] = {"sound_speed_m_s": 1200.0, "chamber_radius_m": 0.15,
                                      "chamber_length_m": 0.40}

    def test_dimensional_form_of_eq_3_5_1_19(self) -> None:
        mode = ci.acoustic_mode_frequency(tangential_order=1, radial_index=1, longitudinal_order=1, **self.KW)
        s = ci.transverse_eigenvalue(1, 1)
        expected = 1200.0 / (2 * pi) * sqrt((s / 0.15) ** 2 + (pi / 0.40) ** 2)
        self.assertLess(abs(mode.frequency_hz / expected - 1.0), 1e-14)
        self.assertEqual(mode.name, "1T1L")

    def test_pure_modes(self) -> None:
        one_l = ci.acoustic_mode_frequency(longitudinal_order=1, **self.KW)
        self.assertLess(abs(one_l.frequency_hz - 1200.0 / 0.8), 1e-9)  # a/2L
        one_t = ci.acoustic_mode_frequency(tangential_order=1, **self.KW)
        self.assertLess(abs(one_t.frequency_hz - ci.transverse_eigenvalue(1, 1) * 1200 / (2 * pi * 0.15)), 1e-9)
        one_r = ci.acoustic_mode_frequency(radial_index=2, **self.KW)
        self.assertEqual(one_r.name, "1R")

    def test_mach_factor_for_longitudinal_and_warning_otherwise(self) -> None:
        still = ci.acoustic_mode_frequency(longitudinal_order=1, **self.KW).frequency_hz
        moving = ci.acoustic_mode_frequency(longitudinal_order=1, mean_mach_number=0.2, **self.KW)
        self.assertLess(abs(moving.frequency_hz / still - 0.96), 1e-14)
        transverse = ci.acoustic_mode_frequency(tangential_order=1, mean_mach_number=0.2, **self.KW)
        self.assertEqual(transverse.warnings[0].code, "mach_correction_not_modelled")

    def test_mode_list_is_sorted_and_bounded(self) -> None:
        modes = ci.chamber_acoustic_modes(maximum_frequency_hz=6000.0, **self.KW)
        freqs = [m.frequency_hz for m in modes]
        self.assertEqual(freqs, sorted(freqs))
        self.assertLessEqual(freqs[-1], 6000.0)
        self.assertEqual(modes[0].name, "1L")
        self.assertIn("1T", [m.name for m in modes])


class ChuggingExample(unittest.TestCase):
    def solve(self, ox_drop: float) -> ci.ChuggingStability:
        e = DATA["example_5_5_1_1"]
        pc = 2.0e6
        return ci.chugging_stability(
            mixture_ratio=e["mixture_ratio"], characteristic_velocity_m_s=e["cstar_ft_s"] * FT,
            cstar_mixture_ratio_derivative_m_s=e["dcstar_dr_ft_s"] * FT, chamber_pressure_pa=pc,
            oxidizer_injector_drop_pa=ox_drop * pc, fuel_injector_drop_pa=e["fuel_drop_over_pc"] * pc,
            gas_residence_time_s=e["residence_time_ms"] * 1e-3, time_lag_s=e["time_lag_ms"] * 1e-3)

    def test_case_a_is_unconditionally_stable(self) -> None:
        e = DATA["example_5_5_1_1"]
        result = self.solve(e["oxidizer_drop_over_pc"])
        self.assertLess(abs(result.gain - e["gain_case_a"]), 5e-4)
        self.assertTrue(result.unconditionally_stable)
        self.assertTrue(result.stable)

    def test_case_b_limit(self) -> None:
        e = DATA["example_5_5_1_1"]
        result = self.solve(e["oxidizer_drop_over_pc_case_b"])
        self.assertLess(abs(result.gain - e["gain_case_b"]), 5e-4)
        self.assertLess(abs(result.critical_time_lag_s * 1e3 - e["critical_lag_ms_case_b"]), 5e-3)
        self.assertFalse(result.stable)  # 2.75 ms > 1.78 ms

    def test_limit_satisfies_the_characteristic_equation(self) -> None:
        e = DATA["example_5_5_1_1"]
        r = self.solve(e["oxidizer_drop_over_pc_case_b"])
        s = 2j * pi * r.neutral_frequency_hz
        theta = e["residence_time_ms"] * 1e-3
        f = r.gain * cmath.exp(-s * r.critical_time_lag_s) / (1 + s * theta)
        self.assertLess(abs(f + 1.0), 1e-12)


class IntrinsicAndNeutralCurves(unittest.TestCase):
    def test_crocco_cheng_limit_satisfies_eq_5_3_1_18(self) -> None:
        theta = 1.0e-3
        for n in (0.6, 0.75, 1.0, 2.0, 5.0):
            tau = ci.intrinsic_critical_time_lag(n, theta)
            omega = sqrt(2 * n - 1) / theta
            residual = n * (1 - cmath.exp(-1j * omega * tau)) - (1 + 1j * omega * theta)
            with self.subTest(n=n):
                self.assertLess(abs(residual), 1e-12)
        self.assertLess(abs(ci.intrinsic_critical_time_lag(1.0, theta) - pi / 2 * theta), 1e-18)

    def test_no_intrinsic_instability_below_one_half(self) -> None:
        self.assertIsNone(ci.intrinsic_critical_time_lag(0.5, 1e-3))
        self.assertIsNone(ci.intrinsic_critical_time_lag(0.3, 1e-3))

    def test_figure_4_2_2a_minimum(self) -> None:
        fig = DATA["figure_4_2_2a"]
        phase = 2 * pi * fig["at_tau_over_tw"]
        self.assertLess(abs(ci.neutral_interaction_index(phase, 1.0) - fig["minimum_n_over_gr"]), 1e-15)
        for tau_over_tw in (0.3, 0.45, 0.55, 0.7, 1.5):
            self.assertGreaterEqual(ci.neutral_interaction_index(2 * pi * tau_over_tw, 1.0), 0.5)
        with self.assertRaises(DomainError):
            ci.neutral_interaction_index(2 * pi, 1.0)


class ResonatorsAndDamping(unittest.TestCase):
    def test_end_correction(self) -> None:
        c = DATA["end_correction"]
        d, a = 0.004, 0.04
        delta = c["coefficient"] * d * (1 - c["open_area_factor"] * sqrt(a))
        self.assertLess(abs(ci.orifice_effective_length(orifice_length_m=0.006, orifice_diameter_m=d,
                                                        open_area_ratio=a) - (0.006 + delta)), 1e-15)
        high = ci.orifice_effective_length(orifice_length_m=0.006, orifice_diameter_m=d,
                                           open_area_ratio=a, regime="high_amplitude")
        self.assertLess(abs(high - (0.006 + c["high_amplitude_factor"] * delta)), 1e-15)

    def test_helmholtz_and_liner_forms_agree_for_a_flat_array(self) -> None:
        """Eq. (8.3.3-1): A_o/V_cav = a/L_cav for a large flat array."""
        a, depth, l_eff, area = 0.05, 0.02, 0.008, 1.0e-5
        volume = area * depth / a  # cell volume behind one orifice
        one = ci.helmholtz_resonator_frequency(sound_speed_m_s=900.0, orifice_area_m2=area,
                                               effective_length_m=l_eff, cavity_volume_m3=volume)
        array = ci.liner_resonant_frequency(sound_speed_m_s=900.0, open_area_ratio=a,
                                            cavity_depth_m=depth, effective_length_m=l_eff)
        self.assertLess(abs(one / array - 1.0), 1e-14)
        self.assertLess(abs(one - 900.0 / (2 * pi) * sqrt(area / (l_eff * volume))), 1e-9)

    def test_q_factor_example(self) -> None:
        e = DATA["q_factor_example"]
        bandwidth = e["frequency_hz"] / e["q"]
        f1, f2 = 550.0 - bandwidth / 2, 550.0 + bandwidth / 2
        self.assertLess(abs(ci.quality_factor(550.0, f1, f2) - 55.0), 1e-12)
        shifted = ci.quality_factor(550.0, f1, f2 + 1.0)
        self.assertLess(abs((55.0 - shifted) / 55.0 - 0.10), 0.01)  # "a 10 % error in Q"

    def test_decay_rate_and_response_factor(self) -> None:
        times = [k * 1e-3 for k in range(10)]
        alpha = 120.0
        amps = [2.0 * cmath.exp(-alpha * t).real for t in times]
        self.assertLess(abs(ci.damping_rate_from_decay(times, amps) - alpha), 1e-9)
        samples = [2 * pi * k / 64 for k in range(64)]
        p = [cos(x) for x in samples]
        m = [0.8 * cos(x - pi / 3) for x in samples]
        self.assertLess(abs(ci.response_factor(p, m) - 0.8 * cos(pi / 3)), 1e-12)  # Eq. (4.4.1-5)
        m_out = [0.8 * cos(x - 2 * pi / 3) for x in samples]
        self.assertLess(ci.response_factor(p, m_out), 0.0)


if __name__ == "__main__":
    unittest.main()
