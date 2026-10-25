"""Combustion efficiency, product composition and finite-area combustor tests.

Reference values:

* ``reference_data/cea_fac_points.json`` -- NASA CEA 3.3.4 finite-area-combustor
  runs and the NASA CEA2 regression outputs of RP-1311 examples 9 and 10
  (``scripts/generate_cea_fac_points.py``);
* ``reference_data/combustion_efficiency.json`` -- NASA CEA 3.3.4 mole and
  mass fractions (with trace radicals), stoichiometric O/F, O/F optima and
  published values read from NASA SP-8089 and SP-125 with page locators
  (``scripts/generate_combustion_efficiency_reference.py``).

Tolerances were declared in the generators before any local value was
compared. The tests need no CEA install. ``RPL_SLOW_TESTS=1`` also runs the
LOX/CH4 and LOX/RP-1 O/F optima (about a minute).
"""

from __future__ import annotations

import json
import os
import unittest
from functools import lru_cache
from math import isclose
from pathlib import Path

from rocket_propulsion.core.errors import DomainError, InputError
from rocket_propulsion.thermochemistry import combustion_efficiency as ce
from rocket_propulsion.thermochemistry.finite_area_combustor import (
    rp1311_stagnation_pressure_estimate,
    solve_finite_area_combustor,
)
from rocket_propulsion.thermochemistry.rocket_equilibrium import PROPELLANT_PAIRS, reactants_for_mass_ratio
from rocket_propulsion.thermochemistry.species import get_reactant, get_species

DATA = Path(__file__).parent / "reference_data"
FAC = json.loads((DATA / "cea_fac_points.json").read_text(encoding="utf-8"))
REF = json.loads((DATA / "combustion_efficiency.json").read_text(encoding="utf-8"))
FAC_POINTS = {point["id"]: point for point in FAC["points"]}
CASES = {case["id"]: case for case in REF["compositions"]}


@lru_cache(maxsize=None)
def _fac(point_id: str):
    point = FAC_POINTS[point_id]
    option = ({"contraction_ratio": point["option_value"]} if point["option"] == "ac_at"
              else {"mass_flux_per_chamber_area_kg_m2_s": point["option_value"]})
    return solve_finite_area_combustor(
        reactants=reactants_for_mass_ratio(PROPELLANT_PAIRS[point["pair"]], point["oxidizer_fuel_ratio"]),
        injector_pressure_pa=point["injector_pressure_pa"],
        area_ratios=tuple(item["area_ratio"] for item in point["expected"]["exits"]),
        **option,
    )


@lru_cache(maxsize=None)
def _report(case_id: str):
    case = CASES[case_id]
    return ce.combustion_product_report(
        chamber_pressure_pa=case["chamber_pressure_pa"], area_ratio=case["area_ratio"],
        propellant_pair=case["pair"], oxidizer_fuel_ratio=case["oxidizer_fuel_ratio"])


def _relative(local: float, reference: float) -> float:
    return abs(local - reference) / abs(reference)


class FiniteAreaCombustorTests(unittest.TestCase):
    def test_fixture_records_cea_provenance(self) -> None:
        self.assertEqual(FAC["reference"]["version"], "3.3.4")
        self.assertEqual(len(FAC["reference"]["thermo_lib_sha256"]), 64)
        self.assertEqual({point["pair"] for point in FAC["points"]}, {"LOX_LH2", "LOX_CH4", "LOX_RP1"})

    def test_matches_cea_334_fac_points_within_declared_tolerances(self) -> None:
        tolerance = FAC["declared_relative_tolerance"]
        for point_id, point in FAC_POINTS.items():
            result, expected = _fac(point_id), point["expected"]
            local = {
                "stagnation_pressure_ratio": result.stagnation_pressure_ratio,
                "combustor_end_pressure_pa": result.combustor_end.pressure_pa,
                "combustor_end_temperature_k": result.combustor_end.temperature_k,
                "combustor_end_mach_number": result.combustor_end_mach_number,
                "throat_pressure_pa": result.throat.pressure_pa,
                "characteristic_velocity_m_s": result.characteristic_velocity_m_s,
                "contraction_ratio": result.contraction_ratio,
            }
            for key, value in local.items():
                with self.subTest(point=point_id, quantity=key):
                    self.assertLessEqual(_relative(value, expected[key]), tolerance[key])
            for local_exit, cea_exit in zip(result.exits, expected["exits"], strict=True):
                with self.subTest(point=point_id, area_ratio=cea_exit["area_ratio"]):
                    self.assertLessEqual(
                        _relative(local_exit.vacuum_effective_exhaust_velocity_m_s,
                                  cea_exit["vacuum_effective_exhaust_velocity_m_s"]),
                        tolerance["vacuum_effective_exhaust_velocity_m_s"])
            with self.subTest(point=point_id, quantity="injector_entropy"):
                self.assertLessEqual(_relative(result.injector.entropy_j_kg_k,
                                               expected["injector_entropy_j_kg_k"]), 1.0e-3)

    def test_matches_nasa_cea2_regression_outputs_of_rp1311_examples_9_and_10(self) -> None:
        outputs = FAC["nasa_reference_outputs"]
        nine, ten = _fac("rp1311-ex9"), _fac("rp1311-ex10")
        e9, e10 = outputs["example9"]["expected"], outputs["example10"]["expected"]
        self.assertLessEqual(_relative(1.0 / nine.stagnation_pressure_ratio,
                                       e9["injector_to_infinite_pressure_ratio"]), 2.0e-6 + 5.0e-7)
        self.assertLessEqual(_relative(nine.combustor_end.pressure_pa, e9["combustor_end_pressure_pa"]), 5.0e-4)
        self.assertLessEqual(_relative(nine.combustor_end.temperature_k, e9["combustor_end_temperature_k"]),
                             1.0e-3)
        self.assertLessEqual(_relative(nine.combustor_end_mach_number, e9["combustor_end_mach_number"]), 2.0e-3)
        self.assertLessEqual(_relative(nine.characteristic_velocity_m_s, e9["characteristic_velocity_m_s"]),
                             1.0e-3)
        for exit_, key in zip(nine.exits[:2], ("area_ratio_25_vacuum_effective_exhaust_velocity_m_s",
                                               "area_ratio_50_vacuum_effective_exhaust_velocity_m_s"), strict=True):
            self.assertLessEqual(_relative(exit_.vacuum_effective_exhaust_velocity_m_s, e9[key]), 1.0e-3)
        self.assertLessEqual(_relative(ten.contraction_ratio, e10["contraction_ratio"]), 2.0e-4)
        self.assertLessEqual(_relative(1.0 / ten.stagnation_pressure_ratio,
                                       e10["injector_to_infinite_pressure_ratio"]), 2.0e-4)
        self.assertLessEqual(_relative(ten.characteristic_velocity_m_s, e10["characteristic_velocity_m_s"]),
                             1.0e-3)

    def test_momentum_balance_closes_and_matches_tp1999_208522_equation_4(self) -> None:
        result = _fac("rp1311-ex9")
        self.assertLessEqual(result.momentum_relative_residual, 1.0e-7)
        p_inf = result.infinite_area.pressure_pa
        ratio = ce.injector_to_throat_total_pressure_ratio(
            static_to_total_pressure_ratio=result.combustor_end.pressure_pa / p_inf,
            subsonic_specific_impulse_m_s=result.combustor_end.velocity_m_s,
            injection_velocity_m_s=0.0,
            characteristic_velocity_m_s=result.characteristic_velocity_m_s,
            contraction_ratio=result.contraction_ratio)
        self.assertTrue(isclose(ratio, p_inf / result.injector_pressure_pa, rel_tol=1.0e-6))

    def test_combustion_in_the_duct_raises_entropy_and_lowers_stagnation_pressure(self) -> None:
        for point_id in FAC_POINTS:
            result = _fac(point_id)
            with self.subTest(point=point_id):
                self.assertGreater(result.entropy_rise_j_kg_k, 0.0)
                self.assertLess(result.stagnation_pressure_ratio, 1.0)
                self.assertLess(result.combustor_end_mach_number, 1.0)
                self.assertTrue(isclose(result.injector.enthalpy_j_kg, result.infinite_area.enthalpy_j_kg,
                                        rel_tol=1.0e-6, abs_tol=50.0))
        self.assertLess(_fac("lox_lh2-cr1.3").stagnation_pressure_ratio,
                        _fac("lox_lh2-cr3").stagnation_pressure_ratio)

    def test_rp1311_initial_estimate_and_input_validation(self) -> None:
        self.assertTrue(isclose(rp1311_stagnation_pressure_estimate(1.58),
                                (1.0257 - 1.2318 * 1.58) / (1.0 - 1.26505 * 1.58)))
        reactants = reactants_for_mass_ratio(PROPELLANT_PAIRS["LOX_LH2"], 6.0)
        with self.assertRaises(InputError):
            solve_finite_area_combustor(reactants=reactants, injector_pressure_pa=5.0e6)
        with self.assertRaises(InputError):
            solve_finite_area_combustor(reactants=reactants, injector_pressure_pa=5.0e6, contraction_ratio=2.0,
                                        mass_flux_per_chamber_area_kg_m2_s=1000.0)
        with self.assertRaises(DomainError):
            solve_finite_area_combustor(reactants=reactants, injector_pressure_pa=5.0e6, contraction_ratio=1.0)
        with self.assertRaises(DomainError):
            solve_finite_area_combustor(reactants=reactants, injector_pressure_pa=5.0e6,
                                        mass_flux_per_chamber_area_kg_m2_s=1.0e5)


class ProductCompositionTests(unittest.TestCase):
    def test_mole_and_mass_fractions_match_cea_by_declared_band(self) -> None:
        tolerance = REF["composition_tolerance"]
        bands = sorted(tolerance["bands"], key=lambda band: -band["minimum_mole_fraction"])
        floor = bands[-1]["minimum_mole_fraction"]
        for case_id, case in CASES.items():
            report = _report(case_id)
            for station in report.stations:
                expected = case["stations"][station.name]
                local_x = {key.upper(): value for key, value in station.mole_fractions}
                local_w = {key.upper(): value for key, value in station.mass_fractions}
                with self.subTest(case=case_id, station=station.name, quantity="T"):
                    self.assertLessEqual(_relative(station.temperature_k, expected["temperature_k"]),
                                         tolerance["temperature_relative"])
                    self.assertLessEqual(_relative(station.molecular_mass_kg_mol * 1000.0,
                                                   expected["molecular_mass_g_mol"]),
                                         tolerance["molecular_mass_relative"])
                for name, value in expected["mole_fractions"].items():
                    if value < floor:
                        continue
                    key = name.upper()
                    try:
                        get_species(key)
                    except InputError:
                        continue  # CEA species absent from the bundled catalog
                    allowed = next(band["relative"] for band in bands if value >= band["minimum_mole_fraction"])
                    with self.subTest(case=case_id, station=station.name, species=key):
                        self.assertIn(key, local_x)
                        self.assertLessEqual(_relative(local_x[key], value), allowed)
                        mass = expected["mass_fractions"][name]
                        self.assertLessEqual(_relative(local_w[key], mass), allowed)

    def test_trace_radicals_are_compared_and_recombine_along_the_nozzle(self) -> None:
        case = CASES["lox_lh2-of6-p68.9"]
        for name in ("H", "O", "OH", "HO2"):
            self.assertGreater(case["stations"]["chamber"]["mole_fractions"][name], 1.0e-5)
        report = _report("lox_lh2-of6-p68.9")
        chamber, throat, exit_ = report.stations
        self.assertGreater(chamber.radical_mole_fraction, throat.radical_mole_fraction)
        self.assertGreater(throat.radical_mole_fraction, exit_.radical_mole_fraction)
        cea_sum = [sum(case["stations"][station]["mole_fractions"].get(name, 0.0)
                       for name in ("H", "O", "OH", "HO2")) for station in ("chamber", "throat", "exit")]
        local_sum = [sum(dict(station.radicals).get(name, 0.0) for name in ("H", "O", "OH", "HO2"))
                     for station in report.stations]
        for local, reference in zip(local_sum, cea_sum, strict=True):
            self.assertLessEqual(_relative(local, reference), 1.0e-2)

    def test_element_and_mass_balances_and_energy_split_identities(self) -> None:
        for case_id in CASES:
            report = _report(case_id)
            for station in (*report.stations, report.reference):
                with self.subTest(case=case_id, station=station.name):
                    # The local Gibbs solver accepts composition residuals of 2e-9.
                    self.assertLessEqual(station.element_relative_residual, 5.0e-9)
                    self.assertLessEqual(station.mass_fraction_sum_residual, 1.0e-12)
                    self.assertLessEqual(abs(sum(value for _, value in station.mole_fractions) - 1.0), 1.0e-12)
                    self.assertTrue(isclose(station.formation_enthalpy_j_kg + station.sensible_enthalpy_j_kg,
                                            station.enthalpy_j_kg, rel_tol=1.0e-12, abs_tol=1.0e-6))
            summary = report.dissociation
            with self.subTest(case=case_id, quantity="dissociation"):
                self.assertTrue(isclose(summary.dissociation_fraction + summary.sensible_fraction, 1.0,
                                        rel_tol=1.0e-12))
                self.assertGreater(summary.dissociation_energy_j_kg, 0.0)
                self.assertGreater(summary.exit_recombination_j_kg, summary.throat_recombination_j_kg)
                self.assertGreater(summary.throat_recombination_j_kg, 0.0)

    def test_assigned_enthalpy_at_298_equals_heat_of_formation(self) -> None:
        # RP-1311 Part I eq. (4.2): H(298.15) = Delta_f H(298.15); polynomial fits agree to a few J/mol.
        for key in ("H2O", "OH", "H", "O", "CO2", "CO", "CH4", "H2", "O2"):
            species = get_species(key)
            with self.subTest(species=key):
                self.assertLessEqual(abs(species.properties(298.15).enthalpy_molar_j_mol
                                         - species.heat_of_formation_j_mol_298), 5.0)

    def test_report_input_validation(self) -> None:
        with self.assertRaises(InputError):
            ce.combustion_product_report(chamber_pressure_pa=5.0e6, area_ratio=10.0)
        with self.assertRaises(InputError):
            ce.combustion_product_report(chamber_pressure_pa=5.0e6, area_ratio=10.0, propellant_pair="LOX_LH2")
        with self.assertRaises(DomainError):
            ce.combustion_product_report(chamber_pressure_pa=5.0e6, area_ratio=10.0, propellant_pair="LOX_LH2",
                                         oxidizer_fuel_ratio=6.0, reference_temperature_k=100.0)


class StoichiometryTests(unittest.TestCase):
    def test_stoichiometric_and_equivalence_ratios_match_cea(self) -> None:
        for record in REF["stoichiometry"]:
            fuel, oxidizer = record["fuel"], record["oxidizer"]
            with self.subTest(pair=record["pair"]):
                self.assertTrue(isclose(ce.stoichiometric_oxidizer_fuel_ratio(fuel, oxidizer),
                                        record["stoichiometric_oxidizer_fuel_ratio"], rel_tol=1.0e-6))
                self.assertTrue(isclose(ce.oxidizer_fuel_ratio_for_chemical_equivalence_ratio(fuel, oxidizer, 0.9),
                                        record["oxidizer_fuel_ratio_at_r_0.9"], rel_tol=1.0e-6))
                for phi in (0.8, 1.25):
                    ratio = record[f"oxidizer_fuel_ratio_at_phi_{phi}"]
                    self.assertTrue(isclose(ce.equivalence_ratio(fuel, oxidizer, ratio), phi, rel_tol=1.0e-6))

    def test_hydrogen_oxygen_stoichiometry_from_formulas(self) -> None:
        # 2 H2 + O2: (o/f)_st = M_O2 / (2 M_H2).
        expected = get_reactant("O2(L)").molecular_mass_kg_mol / (2.0 * get_reactant("H2(L)").molecular_mass_kg_mol)
        self.assertTrue(isclose(ce.stoichiometric_oxidizer_fuel_ratio("H2(L)", "O2(L)"), expected, rel_tol=1.0e-12))
        # CH4 + 2 O2: (o/f)_st = 2 M_O2 / M_CH4.
        expected = 2.0 * get_reactant("O2(L)").molecular_mass_kg_mol / get_reactant("CH4(L)").molecular_mass_kg_mol
        self.assertTrue(isclose(ce.stoichiometric_oxidizer_fuel_ratio("CH4(L)", "O2(L)"), expected, rel_tol=1.0e-12))

    def test_r_equals_phi_when_valences_are_separated(self) -> None:
        # RP-1311 Part I p. 43: r and phi coincide when all positive valences are in the
        # fuel and all negative ones in the oxidizer (e.g. CH4 + O2).
        for fuel in ("H2(L)", "CH4(L)", "RP-1"):
            for ratio in (2.0, 3.0, 6.0):
                report = ce.stoichiometry_report(fuel, "O2(L)", ratio)
                with self.subTest(fuel=fuel, ratio=ratio):
                    self.assertTrue(isclose(report.chemical_equivalence_ratio, report.equivalence_ratio,
                                            rel_tol=1.0e-12))

    def test_density_impulse_reproduces_sp125_table_1_7_rows(self) -> None:
        for row in REF["published"]["sp125_table_1_7_density_impulse_rows"]["rows"]:
            value = ce.density_impulse(row["specific_impulse_s"], row["bulk_density_g_cc"])
            rounding = 0.5 * row["bulk_density_g_cc"] + 0.005 * row["specific_impulse_s"] + 0.5
            with self.subTest(row=row):
                self.assertLessEqual(abs(value - row["density_impulse_s_g_cc"]), rounding)

    def test_bulk_density_is_mass_over_volume(self) -> None:
        rho = ce.bulk_density(6.0, 1140.0, 70.0)
        self.assertTrue(isclose(rho, 7.0 / (6.0 / 1140.0 + 1.0 / 70.0), rel_tol=1.0e-14))
        self.assertTrue(isclose(ce.bulk_density(1.0, 800.0, 800.0), 800.0, rel_tol=1.0e-14))


class EfficiencyModelTests(unittest.TestCase):
    def test_measured_efficiencies_multiply(self) -> None:
        result = ce.measured_efficiencies(chamber_pressure_pa=6.9e6, throat_area_m2=0.01, mass_flow_kg_s=29.0,
                                          vacuum_thrust_n=1.25e5, ideal_characteristic_velocity_m_s=2400.0,
                                          ideal_vacuum_specific_impulse_s=455.0)
        self.assertTrue(isclose(result.characteristic_velocity_m_s, 6.9e6 * 0.01 / 29.0))
        self.assertTrue(isclose(result.specific_impulse_efficiency,
                                result.characteristic_velocity_efficiency * result.thrust_coefficient_efficiency,
                                rel_tol=1.0e-12))

    def test_vaporization_limited_cstar_reduces_to_tn4219_special_cases(self) -> None:
        def cstar(ratio: float) -> float:
            return 2000.0 + 150.0 * ratio - 15.0 * ratio * ratio

        ratio = 5.0
        full = ce.vaporization_limited_cstar_efficiency(characteristic_velocity=cstar, oxidizer_fuel_ratio=ratio,
                                                        oxidizer_vaporized_fraction=1.0, fuel_vaporized_fraction=1.0)
        self.assertEqual(full, 1.0)
        oxidant = 0.6  # eq. (6): fuel fully vaporized
        eq6 = cstar(oxidant * ratio) / cstar(ratio) * (oxidant * ratio + 1.0) / (ratio + 1.0)
        self.assertTrue(isclose(ce.vaporization_limited_cstar_efficiency(
            characteristic_velocity=cstar, oxidizer_fuel_ratio=ratio, oxidizer_vaporized_fraction=oxidant,
            fuel_vaporized_fraction=1.0), eq6, rel_tol=1.0e-14))
        fuel = 0.7  # eq. (7): oxidant fully vaporized
        eq7 = cstar(ratio / fuel) / cstar(ratio) * (ratio + fuel) / (ratio + 1.0)
        self.assertTrue(isclose(ce.vaporization_limited_cstar_efficiency(
            characteristic_velocity=cstar, oxidizer_fuel_ratio=ratio, oxidizer_vaporized_fraction=1.0,
            fuel_vaporized_fraction=fuel), eq7, rel_tol=1.0e-14))
        constant = ce.vaporization_limited_cstar_efficiency(
            characteristic_velocity=lambda _: 2000.0, oxidizer_fuel_ratio=ratio, oxidizer_vaporized_fraction=0.5,
            fuel_vaporized_fraction=0.8)
        self.assertTrue(isclose(constant, (0.5 * ratio + 0.8) / (ratio + 1.0), rel_tol=1.0e-14))

    def test_roccid_streamtubes_conserve_mass_and_oxidizer(self) -> None:
        for em in (0.5, 0.75, 0.9, 1.0):
            for mean in (1.6, 2.7, 6.0):
                split = ce.rupe_streamtube_split(em, mean)
                with self.subTest(em=em, mean=mean):
                    self.assertTrue(isclose(split.low_mass_fraction + split.high_mass_fraction, 1.0))
                    oxidizer = (split.low_mass_fraction * split.low_mixture_ratio / (1.0 + split.low_mixture_ratio)
                                + split.high_mass_fraction * split.high_mixture_ratio
                                / (1.0 + split.high_mixture_ratio))
                    self.assertTrue(isclose(oxidizer, mean / (1.0 + mean), rel_tol=1.0e-12))
        self.assertTrue(isclose(ce.rupe_streamtube_split(1.0, 2.0).mixing_parameter, 0.994, rel_tol=1.0e-12))
        with self.assertRaises(DomainError):
            ce.rupe_streamtube_split(0.2, 2.0)

    def test_energy_release_efficiency_with_cea_consistent_properties(self) -> None:
        cstar, isp = ce.cea_characteristic_velocity_function(propellant_pair="LOX_LH2",
                                                             chamber_pressure_pa=6894760.0, area_ratio=40.0)
        result = ce.energy_release_efficiency(rupe_mixing_efficiency=0.9, oxidizer_fuel_ratio=6.0,
                                              characteristic_velocity=cstar, specific_impulse=isp,
                                              vaporization_efficiency=0.98)
        self.assertLess(result.characteristic_velocity_efficiency, 0.98)
        self.assertGreater(result.characteristic_velocity_efficiency, 0.9)
        self.assertLess(result.specific_impulse_energy_release_efficiency, 0.98)
        perfect = ce.energy_release_efficiency(rupe_mixing_efficiency=0.9, oxidizer_fuel_ratio=6.0,
                                               characteristic_velocity=lambda _: 2300.0,
                                               specific_impulse=lambda _: 450.0)
        self.assertTrue(isclose(perfect.characteristic_velocity_efficiency, 1.0, rel_tol=1.0e-12))

    def test_sp8089_mixing_correlation_uses_table_iv(self) -> None:
        published = REF["published"]["sp8089_table_iv_mixing_factor"]["values"]
        self.assertEqual(ce.SP8089_MIXING_FACTOR, published)
        self.assertTrue(isclose(ce.maximum_mixing_diameter_ratio(
            element_type="1-on-1", outside_density_kg_m3=1000.0, center_density_kg_m3=1000.0,
            center_mass_flow_kg_s=1.0, outside_mass_flow_kg_s=1.0), 1.0))
        ratio = ce.maximum_mixing_diameter_ratio(element_type="2-on-1", outside_density_kg_m3=1140.0,
                                                 center_density_kg_m3=810.0, center_mass_flow_kg_s=1.0,
                                                 outside_mass_flow_kg_s=2.7)
        self.assertTrue(isclose(ratio ** 2, 1.6 * ((1140.0 / 810.0) * (1.0 / 2.7) ** 2) ** 0.7, rel_tol=1.0e-12))
        with self.assertRaises(InputError):
            ce.maximum_mixing_diameter_ratio(element_type="6-on-1", outside_density_kg_m3=1.0,
                                             center_density_kg_m3=1.0, center_mass_flow_kg_s=1.0,
                                             outside_mass_flow_kg_s=1.0)

    def test_characteristic_length_stay_time_and_sp125_ranges(self) -> None:
        published = REF["published"]["sp125_table_4_1_characteristic_length_in"]["values"]
        for key, (low, high) in published.items():
            self.assertEqual(ce.SP125_RECOMMENDED_CHARACTERISTIC_LENGTH_M[key.upper()],
                             (low * 0.0254, high * 0.0254))
        result = ce.chamber_residence(chamber_volume_m3=0.02, throat_area_m2=0.01, mass_flow_kg_s=30.0,
                                      gas_density_kg_m3=4.0, propellant_key="LOX_RP1")
        self.assertTrue(isclose(result.characteristic_length_m, 2.0))
        self.assertTrue(isclose(result.stay_time_s, 0.02 * 4.0 / 30.0))
        self.assertTrue(result.warnings)
        inside = ce.chamber_residence(chamber_volume_m3=0.0115, throat_area_m2=0.01, mass_flow_kg_s=30.0,
                                      gas_density_kg_m3=4.0, propellant_key="LOX_RP1")
        self.assertFalse(inside.warnings)

    def test_jannaf_loss_chain_is_multiplicative_and_matches_loss_budget(self) -> None:
        chain = ce.jannaf_isp_loss_chain(ode_specific_impulse_s=455.0, kinetic_efficiency=0.995,
                                         divergence_efficiency=0.99, boundary_layer_decrement_s=4.0,
                                         characteristic_velocity_efficiency=0.985)
        product = (chain.kinetic_efficiency * chain.divergence_efficiency * chain.boundary_layer_efficiency
                   * chain.energy_release_efficiency)
        self.assertTrue(isclose(chain.overall_efficiency, product, rel_tol=1.0e-12))
        self.assertTrue(isclose(chain.kinetic_loss_s + chain.divergence_loss_s + chain.boundary_layer_loss_s
                                + chain.energy_release_loss_s,
                                chain.ode_specific_impulse_s - chain.delivered_specific_impulse_s, rel_tol=1.0e-12))
        self.assertTrue(isclose(chain.loss_budget.delivered_specific_impulse_s, chain.delivered_specific_impulse_s,
                                rel_tol=1.0e-12))
        with self.assertRaises(DomainError):
            ce.jannaf_isp_loss_chain(ode_specific_impulse_s=455.0, kinetic_efficiency=1.2, divergence_efficiency=1.0,
                                     boundary_layer_decrement_s=0.0, characteristic_velocity_efficiency=1.0)


class MixtureRatioOptimumTests(unittest.TestCase):
    def _check(self, pair: str, objective: str, bracket_key: str, cea_key: str, tolerance: float) -> None:
        record = next(item for item in REF["optima"] if item["pair"] == pair)
        densities = REF["conditions"]["test_densities_kg_m3"]
        fuel, oxidizer = record["fuel"], record["oxidizer"]
        kwargs = {}
        if objective == "density_vacuum_impulse":
            kwargs = {"oxidizer_density_kg_m3": densities[oxidizer], "fuel_density_kg_m3": densities[fuel]}
        result = ce.optimize_mixture_ratio(propellant_pair=pair, chamber_pressure_pa=record["chamber_pressure_pa"],
                                           area_ratio=record["area_ratio"], objective=objective,
                                           bracket=tuple(record[bracket_key]), tolerance=tolerance, **kwargs)
        expected = record[cea_key]
        limits = REF["optimum_tolerance"]
        self.assertLessEqual(abs(result.oxidizer_fuel_ratio - expected["oxidizer_fuel_ratio"]),
                             limits["oxidizer_fuel_ratio_absolute"])
        self.assertLessEqual(_relative(result.objective_value, expected["objective"]), limits["objective_relative"])
        return result

    def test_lox_lh2_maximum_vacuum_isp_matches_cea(self) -> None:
        result = self._check("LOX_LH2", "vacuum_specific_impulse", "bracket", "maximum_vacuum_specific_impulse",
                             0.04)
        self.assertFalse(result.at_bracket_edge)

    @unittest.skipUnless(os.environ.get("RPL_SLOW_TESTS") == "1", "slow: set RPL_SLOW_TESTS=1")
    def test_lox_lh2_density_impulse_optimum_is_the_bracket_edge_as_in_cea(self) -> None:
        result = self._check("LOX_LH2", "density_vacuum_impulse", "density_bracket",
                             "maximum_density_vacuum_impulse", 0.04)
        self.assertTrue(result.at_bracket_edge)
        self.assertTrue(result.warnings)

    @unittest.skipUnless(os.environ.get("RPL_SLOW_TESTS") == "1", "slow: set RPL_SLOW_TESTS=1")
    def test_hydrocarbon_optima_match_cea(self) -> None:
        for pair in ("LOX_CH4", "LOX_RP1"):
            with self.subTest(pair=pair):
                self._check(pair, "vacuum_specific_impulse", "bracket", "maximum_vacuum_specific_impulse", 0.01)
                self._check(pair, "density_vacuum_impulse", "density_bracket", "maximum_density_vacuum_impulse",
                            0.01)

    def test_density_objective_requires_densities(self) -> None:
        with self.assertRaises(InputError):
            ce.optimize_mixture_ratio(propellant_pair="LOX_LH2", chamber_pressure_pa=7.0e6, area_ratio=40.0,
                                      objective="density_vacuum_impulse")
        with self.assertRaises(InputError):
            ce.optimize_mixture_ratio(propellant_pair="LOX_LH2", chamber_pressure_pa=7.0e6, area_ratio=40.0,
                                      objective="thrust")


if __name__ == "__main__":
    unittest.main()
