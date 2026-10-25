"""Gas-side heat transfer against a hot firing: NASA TN D-6396.

Schacht and Quentmeyer measured the hot-gas-side coefficient of a 5-in.
(12.70 cm) throat LOX/GH2 chamber at about 2.09 MPa with copper heat-sink
plugs at six circumferential positions in the cylindrical chamber
(A/A* = 4.64) and at the throat (TN D-6396 Table I, five readings,
``tests/reference_data/hotfire_tn_d6396.json``).

The package's station Bartz coefficient (SP-125 Eqs. 4-13 and 4-14 with the
local wall temperature in the sigma factor) is evaluated with the
CEA-equivalent chamber state of the measured O/F and pressure, gaseous
hydrogen at 298.15 K and liquid oxygen. Measured combustion efficiency was
0.9885 to 1.00, so the ideal c* is used.

What the comparison shows (regression bands, not accuracy targets):

* recovery temperatures agree with the authors' computed ``T_aw`` within
  3.5 %;
* Bartz over-predicts the measured chamber coefficient by 7-30 % and the
  throat coefficient by 13-61 % (circumferential means). The authors'
  companion report TN D-2832 found the same: their fitted throat constant
  ``C = 0.0151`` against Bartz's 0.026 (about 42 % lower), with the chamber
  constant 0.0257 near 0.026. The over-prediction at the throat is a
  property of the correlation, not of this implementation, and it is
  conservative for wall design.
"""

from __future__ import annotations

import json
import unittest
from pathlib import Path
from statistics import mean

from rocket_propulsion.compressible.isentropic import mach_from_isentropic
from rocket_propulsion.heat_transfer.regen_chamber import bartz_coefficient, hot_gas_from_reactants
from rocket_propulsion.thermochemistry.species import get_reactant

DATA = json.loads((Path(__file__).parent / "reference_data" / "hotfire_tn_d6396.json")
                  .read_text(encoding="utf-8"))


def _predictions() -> dict[str, dict[str, float]]:
    geometry = DATA["geometry"]
    throat_radius = geometry["throat_diameter_m"] / 2.0
    h2, o2 = get_reactant("H2"), get_reactant("O2(L)")
    out = {}
    for key, run in DATA["readings"].items():
        ratio = run["oxidizer_fuel_ratio_used"]
        gas = hot_gas_from_reactants(
            reactants={h2.key: 1.0, o2.key: ratio * h2.molecular_mass_kg_mol / o2.molecular_mass_kg_mol},
            chamber_pressure_pa=run["chamber_pressure_mpa"] * 1e6,
            throat_radius_of_curvature_m=geometry["throat_wall_radius_m"])
        mach = mach_from_isentropic("area_ratio_subsonic", geometry["chamber_area_ratio"], gas.gamma)
        h_ch, taw_ch = bartz_coefficient(gas, throat_radius_m=throat_radius,
                                         area_ratio=geometry["chamber_area_ratio"], mach=mach,
                                         wall_temperature_k=mean(run["chamber"]["wall_k"]))
        h_th, taw_th = bartz_coefficient(gas, throat_radius_m=throat_radius, area_ratio=1.0, mach=1.0,
                                         wall_temperature_k=mean(run["throat"]["wall_k"]))
        out[key] = {
            "chamber_ratio": h_ch / 1e3 / mean(run["chamber"]["h_kw_m2_k"]),
            "throat_ratio": h_th / 1e3 / mean(run["throat"]["h_kw_m2_k"]),
            "chamber_taw_ratio": taw_ch / run["chamber"]["adiabatic_wall_k"],
            "throat_taw_ratio": taw_th / run["throat"]["adiabatic_wall_k"],
        }
    return out


class TnD6396HotFire(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.results = _predictions()

    def test_reference_data_is_self_consistent(self) -> None:
        for key, run in DATA["readings"].items():
            with self.subTest(reading=key):
                fuel = run["percent_fuel"]
                self.assertLess(abs((100.0 - fuel) / fuel / run["oxidizer_fuel_ratio_used"] - 1.0), 2e-3)

    def test_recovery_temperature_matches_the_authors(self) -> None:
        for key, r in self.results.items():
            with self.subTest(reading=key):
                self.assertLess(abs(r["chamber_taw_ratio"] - 1.0), 0.035)
                self.assertLess(abs(r["throat_taw_ratio"] - 1.0), 0.035)

    def test_chamber_coefficient_band(self) -> None:
        ratios = [r["chamber_ratio"] for r in self.results.values()]
        for key, value in zip(self.results, ratios):
            with self.subTest(reading=key):
                self.assertGreater(value, 1.0)
                self.assertLess(value, 1.35)
        self.assertLess(abs(mean(ratios) - 1.17), 0.05)

    def test_throat_coefficient_band(self) -> None:
        ratios = [r["throat_ratio"] for r in self.results.values()]
        for key, value in zip(self.results, ratios):
            with self.subTest(reading=key):
                self.assertGreater(value, 1.05)
                self.assertLess(value, 1.70)
        self.assertLess(abs(mean(ratios) - 1.37), 0.05)

    def test_throat_overprediction_agrees_with_tn_d_2832(self) -> None:
        """TN D-2832's throat constant implies measured/0.026-form = 0.58 +/- 13 %."""
        c = DATA["tn_d_2832_correlation_constant"]
        implied = c["bartz_type_constant"] / c["by_station"]["3 bar (throat mean)"]
        sigma = c["throat_one_sigma_percent"] / 100.0
        observed = mean(r["throat_ratio"] for r in self.results.values())
        self.assertLess(abs(observed / implied - 1.0), 2.0 * sigma)


if __name__ == "__main__":
    unittest.main()
