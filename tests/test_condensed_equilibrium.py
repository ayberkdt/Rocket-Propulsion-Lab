"""Equilibrium with condensed-phase products against NASA CEA 3.3.4.

Reference values: ``reference_data/cea_condensed_points.json``, produced by
NASA CEA 3.3.4 with its full product database for the reactant elements
(``scripts/generate_cea_condensed_points.py``). The provider's species set is
the same set, drawn from the same ``thermo.inp``; remaining differences come
from CEA's convergence tolerances (RP-1311 Part I eqs. (3.5)-(3.6c)), from
CEA's area-ratio iteration and from the gas constant (8.31451 in CEA,
8.314462618 J/(mol K) here, a 5.7e-6 relative difference in entropy).

Tolerances were declared after one exploratory comparison, not before: the
largest deviation then was 8.9e-6 (relative throat pressure, LOX/CH4 O/F 0.8;
relative for temperatures, pressures, molar masses, c* and Isp, absolute for
mass and mole fractions), and the tolerances are 3e-5, about 3.4 times that.
Measured values are listed in docs/CONDENSED_EQUILIBRIUM.md.
"""

from __future__ import annotations

import json
import unittest
from functools import cache
from math import log
from pathlib import Path

from rocket_propulsion.core.errors import ConvergenceError, DomainError, InputError
from rocket_propulsion.thermochemistry.condensed_equilibrium import (
    CondensedEquilibriumProblem,
    CondensedEquilibriumProvider,
    condensed_catalog,
    get_condensed_reactant,
    solve_condensed_rocket_performance,
)
from rocket_propulsion.thermochemistry.equilibrium import (
    EquilibriumProblem,
    EquilibriumProvider,
    ReactantAmount,
)
from rocket_propulsion.thermochemistry.nasa_polynomials import UNIVERSAL_GAS_CONSTANT_J_MOL_K
from rocket_propulsion.thermochemistry.rocket_performance import solve_rocket_performance

DATA = Path(__file__).parent / "reference_data"
REFERENCE = json.loads((DATA / "cea_condensed_points.json").read_text(encoding="utf-8"))
BUNDLE_PATH = (Path(__file__).parents[1] / "src" / "rocket_propulsion" / "thermochemistry"
               / "data" / "nasa_cea_condensed.json")
BUNDLE = json.loads(BUNDLE_PATH.read_text(encoding="utf-8"))
ROCKETS = {point["id"]: point for point in REFERENCE["rocket_points"]}
EQUILIBRIA = {point["id"]: point for point in REFERENCE["equilibrium_points"]}

RELATIVE = 3.0e-5  # temperature, pressure, M, MW, c*, vacuum Isp
FRACTION = 3.0e-5  # absolute: condensed mass fractions, gas mole fractions
ENTROPY_RELATIVE = 2.0e-5
ENTHALPY_ABSOLUTE_J_KG = 1.0
ENTHALPY_RELATIVE = 1.0e-5
R_CEA_J_MOL_K = 8.31451  # RP-1311 Part I section 2.1


def _moles(names: list[str], mass_fractions: list[float]) -> dict[str, float]:
    return {name: fraction / get_condensed_reactant(name).molecular_mass_kg_mol
            for name, fraction in zip(names, mass_fractions, strict=True)}


@cache
def _rocket(identifier: str):
    point = ROCKETS[identifier]
    ratios = tuple(station["area_ratio"] for station in point["stations"] if station["name"] == "exit")
    return solve_condensed_rocket_performance(
        reactants=_moles(point["reactants"], point["mass_fractions"]),
        chamber_pressure_pa=point["chamber_pressure_pa"],
        area_ratios=ratios,
    )


def _problem(point: dict, **overrides) -> CondensedEquilibriumProblem:
    reactants = _moles(point["reactants"], point["mass_fractions"])
    if point["problem_type"] == "TP":
        return CondensedEquilibriumProblem.tp(reactants, temperature_k=point["temperature_k"],
                                              pressure_pa=point["pressure_pa"], **overrides)
    return CondensedEquilibriumProblem.hp(reactants, pressure_pa=point["pressure_pa"], **overrides)


@cache
def _equilibrium(identifier: str):
    return CondensedEquilibriumProvider().solve_detailed(_problem(EQUILIBRIA[identifier]))


class BundledDataTests(unittest.TestCase):
    def test_provenance_and_selection_counts(self) -> None:
        source = BUNDLE["source"]
        self.assertEqual(source["commit"], "2f79a647e85737742fb08a127da00d1ea1d069da")
        self.assertEqual(source["sha256"],
                         "fa7746572952d74e249e818a82a35c113829742fb421a308e167185528884363")
        self.assertEqual(REFERENCE["local_data"]["thermo_inp_sha256"], source["sha256"])
        self.assertEqual((len(BUNDLE["gas_species"]), len(BUNDLE["condensed_species"]),
                          len(BUNDLE["reactants"])), (186, 16, 2))
        catalog = condensed_catalog()
        for key in ("C(GR)", "H2O(CR)", "H2O(L)", "AL(CR)", "AL(L)", "AL2O3(A)", "AL2O3(L)"):
            self.assertIn(key, catalog.condensed)
        for key in ("AL", "ALO", "AL2O", "ALCL", "ALCL2", "ALCL3", "ALOH", "HALO2", "HCL", "CL", "CL2"):
            self.assertIn(key, catalog.gases)

    def test_heat_of_formation_identity_from_bundled_coefficients(self) -> None:
        # RP-1311 Part I eqs. (4.1)-(4.2): H(298.15) = Delta_f H(298.15).
        checked = 0
        catalog = condensed_catalog()
        records = [catalog.gases[item["key"].upper()] for item in BUNDLE["gas_species"]]
        records += [catalog.condensed[item["key"].upper()] for item in BUNDLE["condensed_species"]]
        for record in records:
            first = record.intervals[0]
            if not (first.minimum_temperature_k <= 300.0 and first.maximum_temperature_k >= 298.15):
                continue
            reduced = first.enthalpy_molar_j_mol(298.15) / (UNIVERSAL_GAS_CONSTANT_J_MOL_K * 298.15)
            with self.subTest(species=record.key):
                self.assertLessEqual(abs(reduced * R_CEA_J_MOL_K * 298.15
                                         - record.heat_of_formation_j_mol_298), 0.5)
            checked += 1
        self.assertEqual(checked, len(BUNDLE["verification"]["heat_of_formation_identity"]))

    def test_condensed_phases_have_equal_gibbs_energy_at_transitions(self) -> None:
        catalog = condensed_catalog().condensed
        for transition in BUNDLE["verification"]["phase_transitions"]:
            lower = catalog[transition["lower"].upper()]
            upper = catalog[transition["upper"].upper()]
            temperature = transition["temperature_k"]
            self.assertEqual(lower.maximum_temperature_k, upper.minimum_temperature_k)

            def gibbs(species, *, upper=upper, temperature=temperature) -> float:
                interval = species.intervals[0] if species is upper else species.intervals[-1]
                return (interval.enthalpy_molar_j_mol(temperature)
                        - temperature * interval.entropy_molar_j_mol_k(temperature))

            difference = (gibbs(upper) - gibbs(lower)) / (UNIVERSAL_GAS_CONSTANT_J_MOL_K * temperature)
            with self.subTest(transition=transition["lower"]):
                self.assertLess(abs(difference), 1.0e-3)
        alumina = next(item for item in BUNDLE["verification"]["phase_transitions"]
                       if item["lower"] == "AL2O3(a)")
        self.assertEqual(alumina["temperature_k"], 2327.0)
        self.assertAlmostEqual(alumina["transition_enthalpy_j_mol"], 111400.0, delta=1.0)
        self.assertLess(abs(alumina["gibbs_difference_over_rt"]), 1.0e-7)

    def test_reactants_carry_their_298_k_heat_of_formation(self) -> None:
        aluminium = get_condensed_reactant("Al(cr)")
        perchlorate = get_condensed_reactant("NH4CLO4(I)")
        self.assertEqual((aluminium.enthalpy_j_mol, aluminium.temperature_k), (0.0, 298.15))
        self.assertEqual((perchlorate.enthalpy_j_mol, perchlorate.temperature_k), (-295767.0, 298.15))
        self.assertEqual(dict(perchlorate.formula), {"N": 1.0, "H": 4.0, "Cl": 1.0, "O": 4.0})
        for record in BUNDLE["reactants"]:
            self.assertLess(abs(record["polynomial_enthalpy_residual_j_mol"]), 0.5)


class CeaEquilibriumPointTests(unittest.TestCase):
    def test_reference_provenance(self) -> None:
        self.assertEqual(REFERENCE["reference"]["version"], "3.3.4")
        self.assertRegex(REFERENCE["reference"]["thermo_lib_sha256"], r"^[0-9a-f]{64}$")
        self.assertEqual(len(EQUILIBRIA), 7)
        self.assertEqual(len(ROCKETS), 4)

    def test_tp_and_hp_points_match_cea(self) -> None:
        for identifier, point in EQUILIBRIA.items():
            solution = _equilibrium(identifier)
            result = solution.result
            mass = result.total_mass_kg
            with self.subTest(point=identifier):
                self.assertLessEqual(abs(result.temperature_k / point["temperature_k"] - 1.0), RELATIVE)
                self.assertLessEqual(abs(result.molecular_mass_kg_mol * 1e3 / point["molecular_mass_g_mol"]
                                         - 1.0), RELATIVE)
                self.assertLessEqual(abs(solution.mixture_molecular_mass_kg_mol * 1e3
                                         / point["mixture_molecular_mass_g_mol"] - 1.0), RELATIVE)
                self.assertLessEqual(abs(result.enthalpy_j / mass - point["enthalpy_j_kg"]),
                                     ENTHALPY_ABSOLUTE_J_KG + ENTHALPY_RELATIVE * abs(point["enthalpy_j_kg"]))
                self.assertLessEqual(abs(result.entropy_j_k / mass / point["entropy_j_kg_k"] - 1.0),
                                     ENTROPY_RELATIVE)
                condensed = {item.species: item.mass_fraction for item in solution.condensed}
                self.assertEqual(set(condensed), set(point["condensed_mass_fractions"]))
                for key, value in point["condensed_mass_fractions"].items():
                    self.assertLessEqual(abs(condensed[key] - value), FRACTION)
                fractions = {item.species: item.mole_fraction for item in result.species}
                for key, value in point["gas_mole_fractions"].items():
                    self.assertLessEqual(abs(fractions[key] - value), FRACTION)

    def test_expected_phases_are_selected(self) -> None:
        expected = {"h2-o2-tp-250": "H2O(cr)", "h2-o2-tp-300": "H2O(L)", "h2-o2-tp-350": "H2O(L)",
                    "al-o2-tp-2000": "AL2O3(a)", "al-o2-tp-3000": "AL2O3(L)",
                    "al-o2-hp-20bar": "AL2O3(L)", "ch4-o2-of-1.0-tp-1000": "C(gr)"}
        for identifier, phase in expected.items():
            with self.subTest(point=identifier):
                self.assertEqual([item.species for item in _equilibrium(identifier).condensed], [phase])


class CeaRocketPointTests(unittest.TestCase):
    def test_rocket_stations_match_cea(self) -> None:
        for identifier, point in ROCKETS.items():
            output = _rocket(identifier)
            performance = output.performance
            stations = [performance.chamber, performance.throat] + [item.station for item in performance.exits]
            self.assertLessEqual(abs(performance.characteristic_velocity_m_s
                                     / point["characteristic_velocity_m_s"] - 1.0), RELATIVE)
            for station, detail, reference in zip(stations, output.stations, point["stations"], strict=True):
                with self.subTest(point=identifier, station=reference["name"], ratio=reference["area_ratio"]):
                    self.assertLessEqual(abs(station.temperature_k / reference["temperature_k"] - 1.0),
                                         RELATIVE)
                    self.assertLessEqual(abs(station.pressure_pa / reference["pressure_pa"] - 1.0), RELATIVE)
                    self.assertLessEqual(abs(station.molecular_mass_kg_mol * 1e3
                                             / reference["molecular_mass_g_mol"] - 1.0), RELATIVE)
                    self.assertLessEqual(abs(detail.mixture_molecular_mass_kg_mol * 1e3
                                             / reference["mixture_molecular_mass_g_mol"] - 1.0), RELATIVE)
                    condensed = dict(detail.condensed_species)
                    self.assertEqual(set(condensed), set(reference["condensed_mass_fractions"]))
                    for key, value in reference["condensed_mass_fractions"].items():
                        self.assertLessEqual(abs(condensed[key] - value), FRACTION)
                    fractions = dict(station.mole_fractions)
                    for key, value in reference["gas_mole_fractions"].items():
                        self.assertLessEqual(abs(fractions.get(key, 0.0) - value), FRACTION)
            for exit_, reference in zip(performance.exits, point["stations"][2:], strict=True):
                with self.subTest(point=identifier, exit=reference["area_ratio"]):
                    self.assertLessEqual(abs(exit_.vacuum_specific_impulse_s
                                             / reference["vacuum_specific_impulse_s"] - 1.0), RELATIVE)

    def test_graphite_appears_only_in_the_nozzle_at_lox_ch4_of_1_2(self) -> None:
        stations = _rocket("lox-ch4-of-1.2").stations
        self.assertEqual(stations[0].condensed_species, ())
        self.assertEqual(stations[1].condensed_species, ())
        self.assertEqual([key for key, _ in stations[2].condensed_species], ["C(gr)"])

    def test_alumina_freezes_at_its_melting_point_in_the_nozzle(self) -> None:
        stations = _rocket("al-ap-20-80").stations
        phases = [tuple(key for key, _ in item.condensed_species) for item in stations]
        self.assertEqual(phases, [("AL2O3(L)",)] * 4 + [("AL2O3(L)", "AL2O3(a)")] * 2 + [("AL2O3(a)",)])
        for detail in stations[4:6]:
            self.assertAlmostEqual(detail.temperature_k, 2327.0, delta=0.01)
        liquid = [dict(item.condensed_species)["AL2O3(L)"] for item in stations[4:6]]
        self.assertGreater(liquid[0], liquid[1])


class GibbsMinimumTests(unittest.TestCase):
    def _check_minimum(self, identifier: str) -> None:
        solution = _equilibrium(identifier)
        result = solution.result
        potentials = dict(solution.element_potentials)
        catalog = condensed_catalog()
        temperature = result.temperature_k
        log_p = log(result.pressure_pa / 1.0e5)
        # Element conservation, RP-1311 eq. (2.7).
        self.assertLess(result.residuals.maximum_element_relative, 1.0e-9)
        # Phase rule at assigned T and p: condensed phases <= elements.
        self.assertLessEqual(len(solution.condensed), len(potentials))
        # Included condensed species satisfy eq. (2.10) with eq. (2.11).
        for item in solution.condensed:
            species = catalog.condensed[item.species.upper()]
            interval = next(i for i in species.intervals
                            if i.minimum_temperature_k <= temperature <= i.maximum_temperature_k)
            mu = (interval.enthalpy_molar_j_mol(temperature)
                  - temperature * interval.entropy_molar_j_mol_k(temperature)) / (
                UNIVERSAL_GAS_CONSTANT_J_MOL_K * temperature)
            balance = sum(count * potentials[element] for element, count in species.formula)
            self.assertLess(abs(mu - balance), 1.0e-8)
        # Gases above 1e-8 satisfy eq. (2.10).
        condensed_keys = {item.species for item in solution.condensed}
        for item in result.species:
            if item.species in condensed_keys or item.mole_fraction < 1.0e-8:
                continue
            species = catalog.gases[item.species.upper()]
            properties = species.properties(temperature)
            mu = (properties.enthalpy_molar_j_mol - temperature * properties.entropy_molar_j_mol_k) / (
                UNIVERSAL_GAS_CONSTANT_J_MOL_K * temperature) + log(item.moles / solution.gas_moles) + log_p
            balance = sum(count * potentials[element] for element, count in species.formula)
            self.assertLess(abs(mu - balance), 1.0e-6, item.species)
        # No excluded condensed candidate lowers G, eq. (3.7).
        for key, value in solution.inclusion_tests:
            self.assertGreaterEqual(value, -1.0e-9, key)

    def test_converged_states_are_constrained_gibbs_minima(self) -> None:
        for identifier in ("h2-o2-tp-250", "h2-o2-tp-300", "al-o2-tp-2000", "al-o2-tp-3000",
                           "ch4-o2-of-1.0-tp-1000"):
            with self.subTest(point=identifier):
                self._check_minimum(identifier)

    def test_condensed_products_lower_the_gibbs_energy_of_the_gas_only_minimum(self) -> None:
        catalog = condensed_catalog()
        for identifier in ("h2-o2-tp-250", "al-o2-tp-2000", "ch4-o2-of-1.0-tp-1000"):
            point = EQUILIBRIA[identifier]
            full = _equilibrium(identifier).result
            elements = {element for element, _ in full.element_moles}
            gases = tuple(key for key, species in catalog.gases.items()
                          if set(dict(species.formula)) <= elements)
            gas_only = CondensedEquilibriumProvider().solve(_problem(point, candidate_species=gases))
            self.assertFalse({item.species.upper() for item in gas_only.species} & set(catalog.condensed))
            g_full = full.enthalpy_j - full.temperature_k * full.entropy_j_k
            g_gas = gas_only.enthalpy_j - gas_only.temperature_k * gas_only.entropy_j_k
            with self.subTest(point=identifier):
                self.assertLess(g_full, g_gas)


class PhaseSelectionTests(unittest.TestCase):
    def test_graphite_is_removed_when_a_warm_start_carries_it_into_a_carbon_free_state(self) -> None:
        point = EQUILIBRIA["ch4-o2-of-1.0-tp-1000"]
        reactants = _moles(point["reactants"], point["mass_fractions"])
        cold_provider = CondensedEquilibriumProvider()
        start = _equilibrium("ch4-o2-of-1.0-tp-1000")
        problem = CondensedEquilibriumProblem.tp(reactants, temperature_k=1400.0, pressure_pa=1.0e6)
        warm = CondensedEquilibriumProvider().solve_detailed(problem, initial=start.state)
        cold = cold_provider.solve_detailed(problem)
        self.assertIn("remove C(gr)", warm.events)
        self.assertEqual(warm.condensed, ())
        for left, right in zip(warm.result.species, cold.result.species, strict=True):
            self.assertEqual(left.species, right.species)
            self.assertAlmostEqual(left.mole_fraction, right.mole_fraction, delta=1.0e-10)

    def test_assigned_temperature_switches_alumina_phase(self) -> None:
        start = _equilibrium("al-o2-tp-2000")
        point = EQUILIBRIA["al-o2-tp-3000"]
        solution = CondensedEquilibriumProvider().solve_detailed(_problem(point), initial=start.state)
        self.assertEqual(solution.events, ("remove AL2O3(a) (outside range)", "include AL2O3(L)"))
        self.assertAlmostEqual(solution.condensed_mass_fraction,
                               point["condensed_mass_fractions"]["AL2O3(L)"], delta=FRACTION)

    def test_isentropic_state_adds_solid_alumina_at_the_melting_point(self) -> None:
        point = ROCKETS["al-ap-20-80"]
        reactants = _moles(point["reactants"], point["mass_fractions"])
        provider = CondensedEquilibriumProvider()
        chamber = provider.solve_detailed(CondensedEquilibriumProblem.hp(
            reactants, pressure_pa=point["chamber_pressure_pa"]))
        station = point["stations"][4]  # area ratio 60, two-phase alumina
        solution = CondensedEquilibriumProvider().solve_detailed(
            CondensedEquilibriumProblem.sp(reactants, pressure_pa=station["pressure_pa"],
                                           target_entropy_j_k=chamber.result.entropy_j_k),
            initial=chamber.state)
        self.assertEqual(solution.events, ("add AL2O3(a) at transition",))
        self.assertAlmostEqual(solution.result.temperature_k, 2327.0, delta=0.01)
        self.assertLess(solution.result.residuals.energy_relative, 1.0e-9)
        cold = CondensedEquilibriumProvider().solve_detailed(
            CondensedEquilibriumProblem.sp(reactants, pressure_pa=station["pressure_pa"],
                                           target_entropy_j_k=chamber.result.entropy_j_k))
        self.assertAlmostEqual(cold.result.temperature_k, solution.result.temperature_k, delta=1.0e-6)
        self.assertEqual(sorted(item.species for item in cold.condensed), ["AL2O3(L)", "AL2O3(a)"])


class DegenerateStateTests(unittest.TestCase):
    def test_gas_free_equilibrium_is_rejected(self) -> None:
        # Stoichiometric H2/O2 at 300 K and 1 bar is all liquid water; NASA CEA
        # 3.3.4 reports no convergence for this state as well.
        with self.assertRaises(DomainError) as caught:
            CondensedEquilibriumProvider().solve(
                CondensedEquilibriumProblem.tp({"H2": 2.0, "O2": 1.0}, temperature_k=300.0, pressure_pa=1e5))
        self.assertEqual(caught.exception.code, "no_gas_phase")

    def test_single_gas_ratio_singularity_is_resolved(self) -> None:
        # RP-1311 Part I eq. (3.11): stoichiometric H2/O2 at 400 K, only H2O(g).
        result = CondensedEquilibriumProvider().solve(
            CondensedEquilibriumProblem.tp({"H2": 2.0, "O2": 1.0}, temperature_k=400.0, pressure_pa=1e5))
        self.assertEqual(result.species[0].species, "H2O")
        self.assertGreater(result.species[0].mole_fraction, 1.0 - 1.0e-10)
        self.assertLess(result.residuals.maximum_element_relative, 1.0e-9)


class ProviderContractTests(unittest.TestCase):
    def test_provider_implements_the_equilibrium_protocol(self) -> None:
        provider = CondensedEquilibriumProvider()
        self.assertIsInstance(provider, EquilibriumProvider)
        reactants = {"CH4(L)": 1.0, "O2(L)": 0.4}
        base = provider.solve(EquilibriumProblem.hp(reactants, pressure_pa=50.0e5))
        own = CondensedEquilibriumProvider().solve(CondensedEquilibriumProblem.hp(reactants, pressure_pa=50.0e5))
        self.assertAlmostEqual(base.temperature_k, own.temperature_k, delta=1.0e-6)
        self.assertIn("C(gr)", [item.species for item in base.species])

    def test_rocket_performance_accepts_the_provider_for_gas_only_chambers(self) -> None:
        reactants = {"H2(L)": 1.0, "O2(L)": 5.5 * 2.01588 / 31.9988}
        through_base = solve_rocket_performance(reactants=reactants, chamber_pressure_pa=68.9476e5,
                                                area_ratios=(40.0,), provider=CondensedEquilibriumProvider())
        direct = solve_condensed_rocket_performance(reactants=reactants, chamber_pressure_pa=68.9476e5,
                                                    area_ratios=(40.0,)).performance
        self.assertAlmostEqual(through_base.chamber.temperature_k / direct.chamber.temperature_k, 1.0, delta=1e-9)
        self.assertAlmostEqual(through_base.characteristic_velocity_m_s
                               / direct.characteristic_velocity_m_s, 1.0, delta=1.0e-7)
        self.assertAlmostEqual(through_base.exits[0].vacuum_specific_impulse_s
                               / direct.exits[0].vacuum_specific_impulse_s, 1.0, delta=1.0e-7)

    def test_rocket_performance_cannot_freeze_a_chamber_holding_graphite(self) -> None:
        # solve_rocket_performance builds a frozen gas mixture from the chamber
        # species; condensed products are not gas species, hence the dedicated
        # two-phase solver.
        with self.assertRaises(InputError):
            solve_rocket_performance(reactants={"CH4(L)": 1.0, "O2(L)": 0.4}, chamber_pressure_pa=68.9476e5,
                                     area_ratios=(10.0,), provider=CondensedEquilibriumProvider())


class InputValidationTests(unittest.TestCase):
    def test_invalid_problems_fail_closed(self) -> None:
        provider = CondensedEquilibriumProvider()
        with self.assertRaises(InputError):
            provider.solve(EquilibriumProblem.uv({"H2": 2.0, "O2": 1.0}, volume_m3=1.0))
        with self.assertRaises(InputError):
            CondensedEquilibriumProblem.tp({"UNOBTAINIUM": 1.0}, temperature_k=1000.0, pressure_pa=1e5)
        with self.assertRaises(DomainError):
            CondensedEquilibriumProblem.tp({"H2": 0.0}, temperature_k=1000.0, pressure_pa=1e5)
        with self.assertRaises(DomainError):
            provider.solve(CondensedEquilibriumProblem.tp({"H2": 1.0}, temperature_k=1000.0, pressure_pa=0.0))
        with self.assertRaises(DomainError):
            provider.solve(CondensedEquilibriumProblem.tp({"H2": 1.0}, temperature_k=100.0, pressure_pa=1e5))
        with self.assertRaises(InputError):
            provider.solve(CondensedEquilibriumProblem("SP", (ReactantAmount("H2", 1.0),), 1.0e5))
        with self.assertRaises(InputError):
            provider.solve(CondensedEquilibriumProblem.tp({"H2": 1.0}, temperature_k=1000.0, pressure_pa=1e5,
                                                          candidate_species=("H2", "XYZ")))
        with self.assertRaises(DomainError):
            provider.solve(CondensedEquilibriumProblem.tp({"AL(cr)": 1.0, "O2": 1.0}, temperature_k=2000.0,
                                                          pressure_pa=1e5,
                                                          candidate_species=("O2", "O", "AL2O3(a)")))
        with self.assertRaises(DomainError):
            solve_condensed_rocket_performance(reactants={"H2": 2.0, "O2": 1.0}, chamber_pressure_pa=1e6,
                                               area_ratios=(1.0,))
        self.assertTrue(issubclass(ConvergenceError, DomainError))


if __name__ == "__main__":
    unittest.main()
