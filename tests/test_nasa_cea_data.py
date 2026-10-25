"""Bundled NASA CEA thermodynamic subset, condensed reactants and catalog links."""

from __future__ import annotations

import json
import os
import subprocess
import sys
import unittest
from pathlib import Path

from rocket_propulsion.propulsion.propellants import PROPELLANTS
from rocket_propulsion.thermochemistry import (
    PROPELLANT_PAIRS,
    EquilibriumProblem,
    LocalEquilibriumProvider,
    condensed_reactant_catalog,
    get_reactant,
    reactants_for_mass_ratio,
    solve_rocket_performance,
    species_catalog,
)
from rocket_propulsion.thermochemistry.species import CEA_DATA_FILE, CondensedReactant

ROOT = Path(__file__).resolve().parents[1]
DOCUMENT = json.loads(CEA_DATA_FILE.read_text(encoding="utf-8"))
GATE = json.loads((ROOT / "tests" / "reference_data" / "cea_gate_points.json").read_text(encoding="utf-8"))


class BundledSubsetTests(unittest.TestCase):
    def test_provenance_pins_the_nasa_commit_and_hash(self) -> None:
        source = DOCUMENT["source"]
        self.assertEqual(source["repository"], "https://github.com/nasa/cea")
        self.assertEqual(source["commit"], "2f79a647e85737742fb08a127da00d1ea1d069da")
        self.assertEqual(source["sha256"],
                         "fa7746572952d74e249e818a82a35c113829742fb421a308e167185528884363")
        self.assertEqual(source["license"], "Apache-2.0")

    def test_polynomials_reproduce_the_recorded_heat_of_formation(self) -> None:
        catalog = species_catalog()
        for record in DOCUMENT["species"]:
            species = catalog[record["key"].upper()]
            with self.subTest(species=record["key"]):
                if species.minimum_temperature_k > 298.15:
                    continue
                self.assertAlmostEqual(
                    species.properties(298.15).enthalpy_molar_j_mol,
                    record["heat_of_formation_j_mol_298"],
                    delta=1.0e-3 * max(1.0, abs(record["heat_of_formation_j_mol_298"])),
                )

    def test_intervals_are_contiguous(self) -> None:
        for record in DOCUMENT["species"]:
            intervals = record["intervals"]
            for left, right in zip(intervals, intervals[1:], strict=False):
                with self.subTest(species=record["key"]):
                    self.assertEqual(left["maximum_temperature_k"], right["minimum_temperature_k"])
                    self.assertEqual(len(left["coefficients"]), 9)

    def test_condensed_reactants_keep_their_storage_state(self) -> None:
        reactants = condensed_reactant_catalog()
        self.assertEqual(set(reactants), {"H2(L)", "O2(L)", "CH4(L)", "RP-1", "N2O4(L)", "CH6N2(L)"})
        self.assertEqual(reactants["H2(L)"].temperature_k, 20.27)
        self.assertEqual(reactants["O2(L)"].temperature_k, 90.17)
        self.assertEqual(reactants["H2(L)"].enthalpy_j_mol, -9012.0)
        self.assertEqual(dict(reactants["RP-1"].formula), {"C": 1.0, "H": 1.95})

    def test_condensed_reactants_enter_only_the_balance_never_the_products(self) -> None:
        result = LocalEquilibriumProvider().solve(
            EquilibriumProblem.hp({"H2(L)": 1.0, "O2(L)": 0.5}, pressure_pa=5.0e6))
        self.assertTrue(all(item.species not in condensed_reactant_catalog() for item in result.species))
        self.assertLess(result.residuals.maximum_element_relative, 1.0e-9)

    def test_condensed_only_hp_problem_needs_no_reactant_temperature(self) -> None:
        problem = EquilibriumProblem.hp({"RP-1": 1.0, "O2(L)": 1.1}, pressure_pa=7.0e6)
        problem = EquilibriumProblem(**{**{slot: getattr(problem, slot) for slot in problem.__slots__},
                                        "reactant_temperature_k": None})
        result = LocalEquilibriumProvider().solve(problem)
        self.assertGreater(result.temperature_k, 3000.0)


class CatalogLinkTests(unittest.TestCase):
    def test_linked_records_resolve_to_computable_pairs(self) -> None:
        linked = [record for record in PROPELLANTS if record.thermochemistry_pair]
        self.assertEqual({record.thermochemistry_pair for record in linked},
                         {"LOX_LH2", "LOX_RP1", "LOX_CH4", "NTO_MMH"})
        for record in linked:
            pair = PROPELLANT_PAIRS[record.thermochemistry_pair]
            with self.subTest(record=record.key):
                self.assertIsInstance(get_reactant(pair.fuel_reactant), CondensedReactant)
                self.assertIsInstance(get_reactant(pair.oxidizer_reactant), CondensedReactant)
                self.assertTrue(reactants_for_mass_ratio(pair, pair.reference_oxidizer_fuel_ratio))

    def test_unlinked_records_state_why_and_published_values_carry_a_locator(self) -> None:
        for record in PROPELLANTS:
            with self.subTest(record=record.key):
                self.assertTrue(record.computability_note)
                if record.thermochemistry_pair is None:
                    self.assertTrue(record.computability_note.startswith("Not computed"))
                self.assertEqual(record.published_isp_s is None, not record.published_isp_source)


class PressureDomainTests(unittest.TestCase):
    def test_every_pair_matches_cea_at_both_ends_of_the_verified_range(self) -> None:
        points = GATE["pressure_domain_points"]
        self.assertEqual({point["pair"] for point in points}, set(PROPELLANT_PAIRS))
        pressures = sorted({point["chamber_pressure_pa"] for point in points})
        for pair in PROPELLANT_PAIRS.values():
            self.assertEqual(tuple(pressures), pair.verified_chamber_pressure_range_pa)
        for point in points:
            pair = PROPELLANT_PAIRS[point["pair"]]
            result = solve_rocket_performance(
                reactants=reactants_for_mass_ratio(pair, point["oxidizer_fuel_ratio"]),
                chamber_pressure_pa=point["chamber_pressure_pa"],
                area_ratios=(GATE["conditions"]["area_ratio"],),
            )
            measured = {
                "chamber_temperature_k": result.chamber.temperature_k,
                "characteristic_velocity_m_s": result.characteristic_velocity_m_s,
                "vacuum_specific_impulse_s": result.exits[0].vacuum_specific_impulse_s,
            }
            for name, expected in point["expected"].items():
                with self.subTest(pair=point["pair"], pressure=point["chamber_pressure_pa"], quantity=name):
                    self.assertLessEqual(abs(measured[name] / expected - 1.0), 2.0e-5)


@unittest.skipUnless(os.environ.get("RPLAB_THERMO_INP"), "set RPLAB_THERMO_INP to the pinned thermo.inp")
class SyncScriptTests(unittest.TestCase):
    def test_bundle_matches_the_pinned_thermo_inp(self) -> None:
        result = subprocess.run(
            [sys.executable, str(ROOT / "scripts" / "sync_nasa_thermo.py"),
             "--source", os.environ["RPLAB_THERMO_INP"], "--check"],
            capture_output=True, text=True, check=False,
        )
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)


if __name__ == "__main__":
    unittest.main()
