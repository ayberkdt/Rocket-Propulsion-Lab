"""Subcooled flow-boiling burnout correlations against their own test data.

``tests/reference_data/burnout.json`` holds Table III of NASA/TM-1998-206612
(38 ethanol tests) and the constants of NASA CR-134282 Eq. (A-1). The
ethanol correlation is checked the way its authors state it: Eq. (1) within
+/-20 % of the tests near 200 psia (all below 250 psia), and Eqs. (1)-(3)
with the pressure factor across all pressures.
"""

from __future__ import annotations

import json
import unittest
from pathlib import Path
from statistics import mean, pstdev

from rocket_propulsion.core.errors import DomainError
from rocket_propulsion.heat_transfer import convection as cv

DATA = json.loads((Path(__file__).parent / "reference_data" / "burnout.json").read_text(encoding="utf-8"))
FT = 0.3048
PSI = 6894.757293168361
BTU_IN2_S = 1055.05585262 / 0.0254 ** 2


def ethanol(row: dict, *, pressure: bool) -> float:
    velocity = row["velocity_ft_s"]
    subcooling_k = row["velocity_subcooling_f_ft_s"] / velocity * 5.0 / 9.0
    result = cv.ethanol_burnout_heat_flux(velocity * FT, subcooling_k,
                                          pressure_pa=row["pressure_psia"] * PSI if pressure else None)
    return result.value / BTU_IN2_S


class EthanolBurnout(unittest.TestCase):
    TABLE = DATA["meyer_linne_rousar_ethanol"]["table_iii"]

    def test_low_pressure_data_within_the_stated_band(self) -> None:
        rows = [r for r in self.TABLE if r["pressure_psia"] < 250.0]
        self.assertEqual(len(rows), 23)
        errors = [r["chf_btu_in2_s"] / ethanol(r, pressure=False) - 1.0 for r in rows]
        # 22 of 23 inside +/-20 %; test 4 (long tube) is at -23.7 %.
        inside = [e for e in errors if abs(e) <= DATA["meyer_linne_rousar_ethanol"]["stated_band"]]
        self.assertGreaterEqual(len(inside), 22)
        self.assertLess(max(abs(e) for e in errors), 0.24)
        self.assertLess(abs(mean(errors)), 0.03)
        self.assertLess(pstdev(errors), 0.12)

    def test_pressure_corrected_over_all_tests(self) -> None:
        errors = [r["chf_btu_in2_s"] / ethanol(r, pressure=True) - 1.0 for r in self.TABLE]
        self.assertLess(abs(mean(errors)), 0.03)
        self.assertLess(pstdev(errors), 0.14)
        worst = max(self.TABLE, key=lambda r: abs(r["chf_btu_in2_s"] / ethanol(r, pressure=True) - 1.0))
        self.assertEqual(worst["test"], 15)  # 514 psia, the authors' caveat on Fig. 9

    def test_unit_conversion_on_test_1(self) -> None:
        row = self.TABLE[0]
        self.assertLess(abs(ethanol(row, pressure=False) - (0.1003 + 0.05264 * 825 ** 0.5)), 1e-12)

    def test_range_and_regime_notes(self) -> None:
        with self.assertRaises(DomainError):
            cv.ethanol_burnout_heat_flux(3.0 * FT, 10.0)  # V dT_sub = 54
        high = cv.ethanol_burnout_heat_flux(30 * FT, 100.0, pressure_pa=500 * PSI)
        self.assertEqual(high.warnings[-1].code, "mild_transition_regime")
        self.assertLess(high.value, cv.ethanol_burnout_heat_flux(30 * FT, 100.0).value)


class MmhBurnout(unittest.TestCase):
    def test_equation_a1(self) -> None:
        c = DATA["tobin_mmh"]
        for v, dt in ((10, 100), (40, 200), (56, 31)):
            value = cv.mmh_burnout_heat_flux(v * FT, dt * 5.0 / 9.0).value / BTU_IN2_S
            self.assertLess(abs(value / (c["coefficient"] * (v * dt) ** c["exponent"]) - 1.0), 1e-12)

    def test_range(self) -> None:
        with self.assertRaises(DomainError):
            cv.mmh_burnout_heat_flux(80 * FT, 100 * 5.0 / 9.0)
        warned = cv.mmh_burnout_heat_flux(80 * FT, 100 * 5.0 / 9.0, extrapolate=True)
        self.assertEqual(warned.warnings[0].code, "correlation_extrapolated")

    def test_feeds_the_sp8087_margin_screen(self) -> None:
        burnout = cv.mmh_burnout_heat_flux(40 * FT, 200 * 5.0 / 9.0).value
        screen = cv.burnout_margin_screen(0.85 * burnout, burnout)
        self.assertEqual(screen.warnings[0].code, "burnout_margin_exceeded")


if __name__ == "__main__":
    unittest.main()
