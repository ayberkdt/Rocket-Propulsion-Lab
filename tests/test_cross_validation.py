"""Cross-validation between independent models, data sources and physical laws.

Each package is verified against its own primary reference elsewhere (NASA CEA,
NIST WebBook, textbook tables). These tests check what no single reference can:
that independent data sources agree where their ranges overlap, that results
obey thermodynamic ordering and conservation laws across modules, and that
limits connect one model to another. Bands marked "measured 2026-10-25" are the
observed agreement between two independent published data sets and act as
regression guards; the others are laws that must hold exactly or by sign.
"""

from __future__ import annotations

import random
import unittest
from functools import lru_cache

from rocket_propulsion.heat_transfer import conduction, radiation
from rocket_propulsion.thermochemistry import (
    PROPELLANT_PAIRS,
    LocalEquilibriumProvider,
    get_species,
    reactants_for_mass_ratio,
    solve_rocket_performance,
)
from rocket_propulsion.thermochemistry import combustion_efficiency as ce
from rocket_propulsion.thermochemistry.condensed_equilibrium import CondensedEquilibriumProvider
from rocket_propulsion.thermochemistry.finite_area_combustor import solve_finite_area_combustor
from rocket_propulsion.thermodynamics.real_fluids import reference_fluid


@lru_cache(maxsize=None)
def _rocket(pair_key: str, pressure: float, chemistry: str):
    pair = PROPELLANT_PAIRS[pair_key]
    return solve_rocket_performance(
        reactants=reactants_for_mass_ratio(pair, pair.reference_oxidizer_fuel_ratio),
        chamber_pressure_pa=pressure, area_ratios=(10.0, 40.0), chemistry=chemistry)


class IndependentDataSourceTests(unittest.TestCase):
    """Reference Helmholtz EOS (ideal part) against NASA Glenn polynomials."""

    # Largest |cp_EOS(p -> 0) / cp_NASA - 1| measured 2026-10-25 per fluid.
    BANDS = {("nitrogen", "N2"): 1.0e-4, ("oxygen", "O2"): 1.0e-4,
             ("normal_hydrogen", "H2"): 1.5e-3, ("methane", "CH4"): 2.0e-3}

    def test_low_density_heat_capacity_agrees_between_eos_and_nasa_polynomials(self) -> None:
        for (fluid, species), band in self.BANDS.items():
            model = reference_fluid(fluid)
            for temperature in (300.0, 500.0, 900.0):
                if temperature > model.model.t_max_k:
                    continue
                with self.subTest(fluid=fluid, temperature=temperature):
                    eos = model.state_tp(temperature, 100.0).cp_j_kg_k
                    nasa = get_species(species).properties(temperature).cp_j_kg_k
                    self.assertLess(abs(eos / nasa - 1.0), band)

    def test_low_density_state_is_an_ideal_gas(self) -> None:
        for fluid in ("nitrogen", "oxygen", "methane", "normal_hydrogen", "helium"):
            with self.subTest(fluid=fluid):
                # At 1 Pa the second-virial departure B p/(R T) is below 1e-7 for
                # every fluid here (methane at 100 Pa already shows -2e-6).
                state = reference_fluid(fluid).state_tp(280.0, 1.0)
                self.assertAlmostEqual(state.compressibility_factor, 1.0, delta=1.0e-7)


class ThermochemicalOrderingTests(unittest.TestCase):
    """Second-law and expansion orderings that every rocket solution must obey."""

    def test_shifting_equilibrium_never_delivers_less_than_frozen(self) -> None:
        for pair in PROPELLANT_PAIRS:
            shifting = _rocket(pair, 7.0e6, "equilibrium")
            frozen = _rocket(pair, 7.0e6, "frozen")
            for a, b in zip(shifting.exits, frozen.exits, strict=True):
                with self.subTest(pair=pair, area_ratio=a.area_ratio):
                    self.assertGreater(a.vacuum_specific_impulse_s, b.vacuum_specific_impulse_s)

    def test_vacuum_impulse_rises_and_temperature_falls_with_area_ratio(self) -> None:
        for pair in PROPELLANT_PAIRS:
            result = _rocket(pair, 7.0e6, "equilibrium")
            low, high = result.exits
            with self.subTest(pair=pair):
                self.assertGreater(high.vacuum_specific_impulse_s, low.vacuum_specific_impulse_s)
                self.assertLess(high.station.temperature_k, low.station.temperature_k)
                self.assertLess(result.throat.temperature_k, result.chamber.temperature_k)

    def test_entropy_is_constant_along_the_isentrope(self) -> None:
        for pair in PROPELLANT_PAIRS:
            result = _rocket(pair, 7.0e6, "equilibrium")
            reference = result.chamber.entropy_j_kg_k
            for station in (result.throat, *[item.station for item in result.exits]):
                with self.subTest(pair=pair, station=station.name):
                    self.assertAlmostEqual(station.entropy_j_kg_k / reference, 1.0, delta=1.0e-8)

    def test_finite_area_combustor_loses_to_the_infinite_area_limit_and_recovers_it(self) -> None:
        pair = PROPELLANT_PAIRS["LOX_LH2"]
        reactants = reactants_for_mass_ratio(pair, pair.reference_oxidizer_fuel_ratio)
        iac = solve_rocket_performance(reactants=reactants, chamber_pressure_pa=7.0e6, area_ratios=(40.0,))
        narrow = solve_finite_area_combustor(reactants=reactants, injector_pressure_pa=7.0e6,
                                             contraction_ratio=1.5, area_ratios=(40.0,))
        wide = solve_finite_area_combustor(reactants=reactants, injector_pressure_pa=7.0e6,
                                           contraction_ratio=20.0, area_ratios=(40.0,))
        self.assertLess(narrow.stagnation_pressure_ratio, wide.stagnation_pressure_ratio)
        self.assertLess(wide.stagnation_pressure_ratio, 1.0)
        # Vacuum Isp at fixed area ratio depends on the stagnation pressure only weakly;
        # the large-contraction combustor must approach the infinite-area value.
        self.assertAlmostEqual(wide.exits[0].vacuum_specific_impulse_s / iac.exits[0].vacuum_specific_impulse_s,
                               1.0, delta=2.0e-3)

    def test_condensed_provider_matches_gas_only_provider_without_condensed_products(self) -> None:
        pair = PROPELLANT_PAIRS["LOX_LH2"]
        reactants = reactants_for_mass_ratio(pair, pair.reference_oxidizer_fuel_ratio)
        from rocket_propulsion.thermochemistry import EquilibriumProblem

        problem = EquilibriumProblem.hp(reactants, pressure_pa=7.0e6)
        gas = LocalEquilibriumProvider().solve(problem)
        condensed = CondensedEquilibriumProvider().solve(problem)
        self.assertAlmostEqual(condensed.temperature_k / gas.temperature_k, 1.0, delta=1.0e-9)


class CombustionChemistryLawTests(unittest.TestCase):
    """Le Chatelier and recombination trends of the product composition."""

    def test_higher_chamber_pressure_suppresses_dissociation(self) -> None:
        reports = [ce.combustion_product_report(chamber_pressure_pa=pressure, area_ratio=40.0,
                                                propellant_pair="LOX_LH2", oxidizer_fuel_ratio=6.0)
                   for pressure in (1.0e6, 7.0e6, 2.0e7)]
        radicals = [report.stations[0].radical_mole_fraction for report in reports]
        fractions = [report.dissociation.dissociation_fraction for report in reports]
        self.assertEqual(radicals, sorted(radicals, reverse=True))
        self.assertEqual(fractions, sorted(fractions, reverse=True))

    def test_radicals_recombine_along_the_expansion(self) -> None:
        for pair in ("LOX_LH2", "LOX_CH4", "LOX_RP1", "NTO_MMH"):
            report = ce.combustion_product_report(
                chamber_pressure_pa=7.0e6, area_ratio=40.0, propellant_pair=pair,
                oxidizer_fuel_ratio=PROPELLANT_PAIRS[pair].reference_oxidizer_fuel_ratio)
            radicals = [station.radical_mole_fraction for station in report.stations]
            with self.subTest(pair=pair):
                self.assertEqual(radicals, sorted(radicals, reverse=True))
                for station in report.stations:
                    self.assertLess(abs(station.mass_fraction_sum_residual), 1.0e-12)
                    self.assertLess(station.element_relative_residual, 1.0e-8)


class HeatTransferLawTests(unittest.TestCase):
    """Conservation and reciprocity on seeded random problems."""

    def test_gray_enclosure_conserves_energy_on_random_frustums(self) -> None:
        rng = random.Random(20261025)
        for _ in range(25):
            inlet, outlet = rng.uniform(0.1, 0.5), rng.uniform(0.5, 1.5)
            length = rng.uniform(0.2, 2.0)
            factors = radiation.frustum_enclosure_view_factors(
                small_radius_m=inlet, large_radius_m=outlet, length_m=length)
            summation, reciprocity = radiation.view_factor_matrix_residuals(
                factors.areas_m2, factors.view_factors)
            self.assertLess(summation, 1.0e-12)
            self.assertLess(reciprocity, 1.0e-12)

    def test_transient_wall_conserves_energy_on_random_layups(self) -> None:
        rng = random.Random(7123)
        for _ in range(6):
            layers = [conduction.WallLayer(thickness_m=rng.uniform(1e-3, 5e-3),
                                           thermal_conductivity=rng.uniform(10.0, 400.0),
                                           volumetric_heat_capacity=rng.uniform(1.0e6, 4.0e6))
                      for _ in range(rng.randint(1, 3))]
            result = conduction.solve_transient_wall(
                layers, conduction.WallGeometry(kind="plane"),
                conduction.SurfaceConvection(heat_transfer_coefficient_w_m2_k=rng.uniform(1e3, 3e4),
                                             fluid_temperature_k=rng.uniform(40.0, 300.0)),
                conduction.SurfaceConvection(heat_transfer_coefficient_w_m2_k=rng.uniform(1e3, 3e4),
                                             fluid_temperature_k=rng.uniform(1500.0, 3000.0)),
                initial_temperature_k=300.0, time_step_s=rng.uniform(0.01, 0.1), step_count=60,
                cells_per_layer=12)
            scale = max(abs(result.stored_energy_change_j), 1.0)
            self.assertLess(abs(result.energy_residual_j) / scale, 1.0e-9)


if __name__ == "__main__":
    unittest.main()
