"""Priem-Heidmann generalized length and the vaporization-limited c* efficiency.

The 35-point curve is checked against an independent transcription of NASA
CR-168066's ``DATA LG/EV`` statement; the generalized length against the
reference conditions of the CR-168066 program (300 psia, 140 Btu/lb,
M = 100, r = 0.003 in., V = 100 ft/s) and against the exponents it shares
with Priem's NACA TN 3883 effective length; the c* efficiency against NACA
TN 4219 Eq. (5).
"""

from __future__ import annotations

import unittest

from rocket_propulsion.core.errors import DomainError
from rocket_propulsion.thermochemistry import vaporization as vp

INCH = 0.0254
PSI = 6894.757293168361
FT = 0.3048

# Independent transcription of CR-168066 p. C-6 (DATA LG/EV), inches and percent.
LG = (0.01, 0.02, 0.04, 0.06, 0.1, 0.2, 0.4, 0.6, 0.8, 1.0, 1.5, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0,
      8.0, 9.0, 10.0, 13.0, 15.0, 17.0, 20.0, 25.0, 30.0, 35.0, 40.0, 45.0, 50.0, 60.0, 70.0,
      80.0, 90.0, 100.0)
EV = (0.083, 0.31, 0.91, 1.61, 3.10, 6.5, 12.2, 17.3, 21.8, 25.9, 34.2, 41.5, 52.2, 60.0, 66.0,
      70.5, 74.5, 77.8, 80.5, 82.7, 87.8, 90.2, 92.0, 94.0, 96.2, 97.6, 98.5, 99.1, 99.4, 99.6,
      99.82, 99.92, 99.96, 99.98, 99.99)

REFERENCE = {"chamber_length_m": 10 * INCH, "chamber_pressure_pa": 300 * PSI, "contraction_ratio": 1.0,
             "liquid_temperature_k": 1e-9, "critical_temperature_k": 500.0,
             "heat_of_vaporization_j_kg": 140 * 2326.0, "molecular_weight": 100.0,
             "drop_radius_m": 0.003 * INCH, "injection_velocity_m_s": 100 * FT}


class Curve(unittest.TestCase):
    def test_table_nodes(self) -> None:
        self.assertEqual(len(vp.PRIEM_HEIDMANN_CURVE), 35)
        for (x, y), lg, ev in zip(vp.PRIEM_HEIDMANN_CURVE, LG, EV):
            self.assertEqual((x, y), (lg, ev))
            self.assertLess(abs(vp.percent_vaporized(lg * INCH) - ev), 1e-12)

    def test_linear_interpolation_and_monotonicity(self) -> None:
        self.assertLess(abs(vp.percent_vaporized(2.5 * INCH) - (41.5 + 52.2) / 2), 1e-12)
        values = [vp.percent_vaporized(x * INCH) for x in (0.01, 0.3, 1.2, 7.7, 33.0, 99.0)]
        self.assertEqual(values, sorted(values))

    def test_outside_the_table(self) -> None:
        for x in (0.005, 150.0):
            with self.assertRaises(DomainError):
                vp.percent_vaporized(x * INCH)


class GeneralizedLength(unittest.TestCase):
    def test_reference_conditions_give_the_chamber_length(self) -> None:
        self.assertLess(abs(vp.generalized_vaporization_length(**REFERENCE) / (10 * INCH) - 1.0), 1e-8)

    def test_exponents(self) -> None:
        base = vp.generalized_vaporization_length(**REFERENCE)
        cases = {"chamber_pressure_pa": (2.0, 0.66), "contraction_ratio": (2.0, -0.44),
                 "heat_of_vaporization_j_kg": (2.0, -0.8), "molecular_weight": (2.0, -0.35),
                 "drop_radius_m": (2.0, -1.45), "injection_velocity_m_s": (2.0, -0.75),
                 "chamber_length_m": (2.0, 1.0)}
        for name, (factor, exponent) in cases.items():
            with self.subTest(name=name):
                value = vp.generalized_vaporization_length(**{**REFERENCE, name: REFERENCE[name] * factor})
                self.assertLess(abs(value / base - factor ** exponent), 1e-12)
        hot = vp.generalized_vaporization_length(**{**REFERENCE, "liquid_temperature_k": 250.0})
        self.assertLess(abs(hot / base - 0.5 ** -0.4), 1e-8)

    def test_shares_tn_3883_exponents(self) -> None:
        kw = {"chamber_length_m": 0.2, "chamber_pressure_pa": 2e6, "final_gas_velocity_m_s": 240.0,
              "injection_velocity_m_s": 30.0, "mass_median_radius_m": 75e-6, "reduced_temperature": 0.5}
        base = vp.priem_effective_length(**kw)
        for name, exponent in (("chamber_pressure_pa", 0.66), ("injection_velocity_m_s", -0.75),
                               ("mass_median_radius_m", -1.45)):
            value = vp.priem_effective_length(**{**kw, name: kw[name] * 2.0})
            self.assertLess(abs(value / base - 2.0 ** exponent), 1e-12)

    def test_supercritical_injection_is_refused(self) -> None:
        with self.assertRaises(DomainError):
            vp.generalized_vaporization_length(**{**REFERENCE, "liquid_temperature_k": 600.0})


class CstarEfficiency(unittest.TestCase):
    def test_full_vaporization_approaches_unity(self) -> None:
        result = vp.vaporization_cstar_efficiency(
            characteristic_velocity=lambda r: 2300.0 - 20.0 * (r - 3.0) ** 2, oxidizer_fuel_ratio=3.0,
            oxidizer_generalized_length_m=100 * INCH, fuel_generalized_length_m=100 * INCH)
        self.assertLess(abs(result.cstar_efficiency - 1.0), 2e-4)

    def test_tn_4219_eq_5(self) -> None:
        def cstar(r: float) -> float:
            return 2300.0 - 20.0 * (r - 3.0) ** 2
        result = vp.vaporization_cstar_efficiency(
            characteristic_velocity=cstar, oxidizer_fuel_ratio=3.0,
            oxidizer_generalized_length_m=2.0 * INCH, fuel_generalized_length_m=10.0 * INCH)
        o, f = 0.415, 0.827
        expected = cstar(o * 3.0 / f) / cstar(3.0) * (o * 3.0 + f) / 4.0
        self.assertLess(abs(result.cstar_efficiency - expected), 1e-12)
        self.assertEqual((result.oxidizer_percent_vaporized, result.fuel_percent_vaporized), (41.5, 82.7))


if __name__ == "__main__":
    unittest.main()
