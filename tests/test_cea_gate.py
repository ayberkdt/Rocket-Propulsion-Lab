"""Phase 2 release gate: local equilibrium rocket performance against NASA CEA.

The reference values in ``reference_data/cea_gate_points.json`` were produced
by NASA CEA 3.3.4 (``scripts/generate_cea_gate.py``) with tolerances declared
before the local solver was compared. The tests here need no CEA install;
the ``cea``-dependent tests are skipped when the package is absent.
"""

from __future__ import annotations

import json
import unittest
from functools import cache
from pathlib import Path
from unittest.mock import patch

from rocket_propulsion.core.errors import DomainError, FeatureUnavailableError, InputError
from rocket_propulsion.thermochemistry.cea_library import cea_library_info, solve_cea_rocket
from rocket_propulsion.thermochemistry.rocket_performance import solve_rocket_performance
from rocket_propulsion.thermochemistry.species import get_reactant

DATA = Path(__file__).parent / "reference_data"
GATE = json.loads((DATA / "cea_gate_points.json").read_text(encoding="utf-8"))
MATRIX = json.loads((DATA / "validation_matrix.json").read_text(encoding="utf-8"))


def _reactants(fuel: str, oxidizer: str, ratio: float) -> dict[str, float]:
    fuel_mass = get_reactant(fuel).molecular_mass_kg_mol
    oxidizer_mass = get_reactant(oxidizer).molecular_mass_kg_mol
    return {fuel: 1.0, oxidizer: ratio * fuel_mass / oxidizer_mass}


@cache
def _local(fuel: str, oxidizer: str, ratio: float, pressure: float, areas: tuple[float, ...]):
    return solve_rocket_performance(
        reactants=_reactants(fuel, oxidizer, ratio), chamber_pressure_pa=pressure, area_ratios=areas
    )


class NinePointGateTests(unittest.TestCase):
    def test_fixture_declares_nine_points_over_three_pairs_and_its_provenance(self) -> None:
        self.assertEqual(len(GATE["points"]), 9)
        self.assertEqual({point["pair"] for point in GATE["points"]}, {"LOX_LH2", "LOX_CH4", "LOX_RP1"})
        self.assertEqual(GATE["reference"]["version"], "3.3.4")
        self.assertRegex(GATE["reference"]["thermo_lib_sha256"], r"^[0-9a-f]{64}$")
        self.assertEqual(set(GATE["declared_relative_tolerance"]),
                         {"chamber_temperature_k", "characteristic_velocity_m_s",
                          "vacuum_specific_impulse_s"})

    def test_local_solver_matches_cea_within_declared_tolerance(self) -> None:
        pressure = GATE["conditions"]["chamber_pressure_pa"]
        area = GATE["conditions"]["area_ratio"]
        for point in GATE["points"]:
            result = _local(point["fuel"], point["oxidizer"], point["oxidizer_fuel_ratio"], pressure, (area,))
            measured = {
                "chamber_temperature_k": result.chamber.temperature_k,
                "characteristic_velocity_m_s": result.characteristic_velocity_m_s,
                "vacuum_specific_impulse_s": result.exits[0].vacuum_specific_impulse_s,
            }
            for name, tolerance in GATE["declared_relative_tolerance"].items():
                with self.subTest(point=point["id"], quantity=name):
                    expected = point["expected"][name]
                    self.assertLessEqual(abs(measured[name] / expected - 1.0), tolerance)
                    # Same NASA data and the same problem: agreement is far tighter.
                    self.assertLessEqual(abs(measured[name] / expected - 1.0), 2.0e-5)

    def test_rp1311_example_8_is_reproduced_at_its_published_station(self) -> None:
        case = next(item for item in MATRIX["cases"] if item["id"] == "cea-rp1311-example-8-lox-lh2")
        conditions = case["conditions"]
        result = _local(conditions["fuel"], conditions["oxidizer"], conditions["oxidizer_fuel_ratio"],
                        conditions["chamber_pressure_pa"], (25.0, case["station"]["area_ratio"]))
        exit_75 = result.exits[1]
        measured = {
            "chamber_temperature_k": result.chamber.temperature_k,
            "characteristic_velocity_m_s": result.characteristic_velocity_m_s,
            "vacuum_effective_exhaust_velocity_m_s": exit_75.vacuum_effective_exhaust_velocity_m_s,
            "exit_effective_exhaust_velocity_m_s": exit_75.effective_exhaust_velocity_m_s,
        }
        for name, expected in case["expected"].items():
            with self.subTest(quantity=name):
                self.assertLessEqual(abs(measured[name] / expected - 1.0), case["relative_tolerance"])
        nasa = GATE["nasa_reference_output_example8"]["expected"]
        self.assertAlmostEqual(result.chamber.temperature_k, nasa["chamber_temperature_k"], delta=0.01)
        self.assertAlmostEqual(result.throat.temperature_k, nasa["throat_temperature_k"], delta=0.05)
        self.assertAlmostEqual(result.characteristic_velocity_m_s,
                               nasa["characteristic_velocity_m_s"], delta=0.05)
        self.assertAlmostEqual(result.exits[0].vacuum_effective_exhaust_velocity_m_s,
                               nasa["area_ratio_25_vacuum_effective_exhaust_velocity_m_s"], delta=0.05)

    def test_rocket_stations_are_physically_ordered(self) -> None:
        point = GATE["points"][1]
        result = _local(point["fuel"], point["oxidizer"], point["oxidizer_fuel_ratio"],
                        GATE["conditions"]["chamber_pressure_pa"], (GATE["conditions"]["area_ratio"],))
        chamber, throat, exit_station = result.chamber, result.throat, result.exits[0].station
        self.assertGreater(chamber.pressure_pa, throat.pressure_pa)
        self.assertGreater(throat.pressure_pa, exit_station.pressure_pa)
        self.assertGreater(chamber.temperature_k, throat.temperature_k)
        self.assertGreater(exit_station.velocity_m_s, throat.velocity_m_s)
        self.assertAlmostEqual(throat.entropy_j_kg_k, chamber.entropy_j_kg_k, delta=1.0e-6 * chamber.entropy_j_kg_k)
        self.assertAlmostEqual(exit_station.entropy_j_kg_k, chamber.entropy_j_kg_k,
                               delta=1.0e-6 * chamber.entropy_j_kg_k)

    def test_invalid_rocket_inputs_are_rejected(self) -> None:
        reactants = _reactants("H2(L)", "O2(L)", 5.5)
        with self.assertRaises(DomainError):
            solve_rocket_performance(reactants=reactants, chamber_pressure_pa=-1.0)
        with self.assertRaises(DomainError):
            solve_rocket_performance(reactants=reactants, chamber_pressure_pa=5.0e6, area_ratios=(0.5,))
        with self.assertRaises(InputError):
            solve_rocket_performance(reactants=reactants, chamber_pressure_pa=5.0e6, chemistry="kinetic")


@unittest.skipUnless(cea_library_info().available, "NASA cea package not installed")
class InstalledCeaTests(unittest.TestCase):
    def test_installed_cea_reproduces_the_recorded_points(self) -> None:
        conditions = GATE["conditions"]
        for point in GATE["points"]:
            result = solve_cea_rocket(
                fuel=point["fuel"], oxidizer=point["oxidizer"],
                oxidizer_fuel_ratio=point["oxidizer_fuel_ratio"],
                chamber_pressure_pa=conditions["chamber_pressure_pa"],
                area_ratios=(conditions["area_ratio"],),
            )
            with self.subTest(point=point["id"]):
                self.assertAlmostEqual(result.chamber.temperature_k,
                                       point["expected"]["chamber_temperature_k"], delta=1.0e-6)
                self.assertAlmostEqual(result.exits[0].vacuum_specific_impulse_s,
                                       point["expected"]["vacuum_specific_impulse_s"], delta=1.0e-6)


class MissingCeaTests(unittest.TestCase):
    def test_missing_package_raises_feature_unavailable(self) -> None:
        cea_library_info.cache_clear()
        try:
            with patch('rocket_propulsion.thermochemistry.cea_library._import_cea',
                       side_effect=ModuleNotFoundError('deliberately unavailable CEA')):
                self.assertFalse(cea_library_info().available)
                with self.assertRaises(FeatureUnavailableError):
                    solve_cea_rocket(fuel="H2(L)", oxidizer="O2(L)", oxidizer_fuel_ratio=5.5,
                                     chamber_pressure_pa=5.0e6)
        finally:
            cea_library_info.cache_clear()


if __name__ == "__main__":
    unittest.main()
