"""Conservation, state closure, CEA rendering, and parser tests."""

import unittest

from rocket_propulsion.core.errors import ConvergenceError, FeatureUnavailableError
from rocket_propulsion.thermochemistry import (
    CeaReactant,
    CeaRocketProblem,
    CeaSubprocessAdapter,
    EquilibriumProblem,
    EquilibriumProvider,
    LocalEquilibriumProvider,
    cea_output_warnings,
    parse_cea_plot,
    parse_fortran_number,
    render_cea_equilibrium_input,
    render_cea_rocket_input,
)


class LocalEquilibriumTests(unittest.TestCase):
    def setUp(self) -> None:
        self.provider = LocalEquilibriumProvider()

    def test_provider_contract_and_tp_element_balance(self) -> None:
        self.assertIsInstance(self.provider, EquilibriumProvider)
        result = self.provider.solve(
            EquilibriumProblem.tp(
                {"H2": 2.0, "O2": 1.0},
                temperature_k=3_500.0,
                pressure_pa=5_000_000.0,
            )
        )
        fractions = {item.species: item.mole_fraction for item in result.species}
        self.assertGreater(fractions["H2O"], 0.6)
        self.assertGreater(fractions["OH"], 0.05)
        self.assertLess(result.residuals.maximum_element_relative, 1.0e-9)
        self.assertTrue(result.convergence.converged)

    def test_hp_combustion_closes_reactant_enthalpy(self) -> None:
        result = self.provider.solve(
            EquilibriumProblem.hp(
                {"H2": 2.0, "O2": 1.0},
                pressure_pa=5_000_000.0,
                reactant_temperature_k=298.15,
            )
        )
        self.assertGreater(result.temperature_k, 3_500.0)
        self.assertLess(result.temperature_k, 3_800.0)
        self.assertLess(result.residuals.energy_relative, 1.0e-5)
        self.assertGreater(result.convergence.temperature_iterations, 0)

    def test_uv_combustion_closes_energy_and_volume(self) -> None:
        problem = EquilibriumProblem.uv(
            {"H2": 2.0, "O2": 1.0},
            volume_m3=0.07436,
            reactant_temperature_k=298.15,
        )
        result = self.provider.solve(problem)
        self.assertGreater(result.pressure_pa, 100_000.0)
        self.assertLess(result.residuals.energy_relative, 1.0e-7)
        self.assertLess(result.residuals.volume_relative, 1.0e-8)


class CeaAdapterTests(unittest.TestCase):
    def test_fortran_number_and_plot_parser(self) -> None:
        self.assertAlmostEqual(parse_fortran_number("9.1864-05"), 9.1864e-5)
        self.assertAlmostEqual(parse_fortran_number("1.250D+03"), 1_250.0)
        table = parse_cea_plot(
            "3389.270 53.3172 12.723 1.1449\n3190.532 30.668 12.8 1.15\n",
            ("t", "p", "mw", "gam"),
        )
        self.assertEqual(len(table.rows), 2)
        self.assertAlmostEqual(table.rows[0]["t"], 3389.27)

    def test_generic_and_rocket_input_decks_are_deterministic(self) -> None:
        generic = render_cea_equilibrium_input(
            EquilibriumProblem.hp(
                {"H2": 2.0, "O2": 1.0},
                pressure_pa=5_331_720.0,
                reactant_temperature_k=298.15,
                candidate_species=("H2", "O2", "H2O", "OH", "H", "O"),
            )
        )
        self.assertIn("problem hp p,bar=53.3172", generic)
        self.assertIn("only H2 O2 H2O OH H O", generic)
        rocket = render_cea_rocket_input(
            CeaRocketProblem(
                fuel=CeaReactant("H2(L)", "fuel", 100.0, "wt%", 20.27),
                oxidizer=CeaReactant("O2(L)", "oxid", 100.0, "wt%", 90.17),
                oxidizer_fuel_ratio=5.55157,
                chamber_pressure_pa=5_331_720.0,
                area_ratios=(25.0, 50.0, 75.0),
            )
        )
        self.assertIn("problem rocket", rocket)
        self.assertIn("o/f=5.55157", rocket)
        self.assertIn("supersonic,ae/at=25,50,75", rocket)

    def test_missing_cea_installation_is_an_optional_feature_error(self) -> None:
        adapter = CeaSubprocessAdapter(None)
        with self.assertRaises(FeatureUnavailableError) as context:
            adapter.run("end\n", plot_columns=("t",))
        self.assertEqual(context.exception.code, "feature_unavailable")

    def test_cea_non_convergence_is_not_silently_accepted(self) -> None:
        with self.assertRaises(ConvergenceError) as context:
            cea_output_warnings("WARNING!! NO CONVERGENCE FOR POINT 2")
        self.assertEqual(context.exception.code, "cea_non_convergence")


if __name__ == "__main__":
    unittest.main()

