"""Reference, boundary, inverse, and mixture tests for thermochemistry."""

import unittest

from rocket_propulsion.core.errors import DomainError
from rocket_propulsion.thermochemistry import Mixture, get_species, species_catalog


class SpeciesPropertyTests(unittest.TestCase):
    def test_catalog_contains_initial_rocket_species(self) -> None:
        self.assertEqual(
            set(species_catalog()),
            {"H2", "O2", "H2O", "OH", "H", "O", "N2", "CO", "CO2", "CH4"},
        )

    def test_published_298_k_reference_values(self) -> None:
        references = {
            # cp [J/mol/K], h [J/mol], s [J/mol/K], NASA/TP-2002-211556 Appendix B.
            "H2": (28.836, 0.0, 130.681),
            "H2O": (33.588, -241_826.0, 188.829),
            "CH4": (35.691, -74_600.0, 186.371),
        }
        for key, (cp, enthalpy, entropy) in references.items():
            with self.subTest(species=key):
                result = get_species(key).properties(298.15)
                self.assertAlmostEqual(result.cp_molar_j_mol_k, cp, delta=0.002)
                self.assertAlmostEqual(result.enthalpy_molar_j_mol, enthalpy, delta=2.0)
                self.assertAlmostEqual(result.entropy_molar_j_mol_k, entropy, delta=0.002)

    def test_piecewise_fit_is_continuous_at_1000_k(self) -> None:
        for species in species_catalog().values():
            with self.subTest(species=species.key):
                below = species.properties(1000.0 - 1.0e-5)
                above = species.properties(1000.0)
                self.assertAlmostEqual(
                    below.cp_molar_j_mol_k, above.cp_molar_j_mol_k, delta=0.003
                )
                self.assertAlmostEqual(
                    below.enthalpy_molar_j_mol, above.enthalpy_molar_j_mol, delta=1.0
                )

    def test_no_silent_temperature_extrapolation(self) -> None:
        with self.assertRaises(DomainError) as context:
            get_species("H2O").properties(6_001.0)
        self.assertEqual(context.exception.code, "property_temperature_out_of_range")

    def test_species_enthalpy_and_entropy_inverse(self) -> None:
        species = get_species("H2O")
        state = species.properties(3_500.0)
        self.assertAlmostEqual(
            species.temperature_from_enthalpy(state.enthalpy_j_kg), 3_500.0, places=5
        )
        self.assertAlmostEqual(
            species.temperature_from_entropy(state.entropy_j_kg_k), 3_500.0, places=5
        )


class MixturePropertyTests(unittest.TestCase):
    def test_mole_and_mass_fraction_conversion_round_trip(self) -> None:
        by_mole = Mixture.from_fractions({"N2": 0.79, "O2": 0.21}, basis="mole")
        mass_fractions = {
            component.species.key: component.mass_fraction
            for component in by_mole.components
        }
        by_mass = Mixture.from_fractions(mass_fractions, basis="mass")
        for left, right in zip(by_mole.components, by_mass.components, strict=True):
            self.assertEqual(left.species.key, right.species.key)
            self.assertAlmostEqual(left.mole_fraction, right.mole_fraction, places=12)

    def test_air_like_mixture_at_300_k(self) -> None:
        state = Mixture.from_fractions({"N2": 0.79, "O2": 0.21}).properties(300.0)
        self.assertAlmostEqual(state.molecular_mass_kg_mol, 0.028850334, places=12)
        self.assertAlmostEqual(state.gas_constant_j_kg_k, 288.19294, places=4)
        self.assertAlmostEqual(state.gamma, 1.398476, places=5)

    def test_mixture_fraction_sum_is_not_silently_normalized(self) -> None:
        with self.assertRaises(DomainError) as context:
            Mixture.from_fractions({"N2": 0.8, "O2": 0.3})
        self.assertEqual(context.exception.code, "composition_sum")

    def test_mixture_enthalpy_and_entropy_inverse(self) -> None:
        mixture = Mixture.from_fractions({"H2O": 0.7, "CO2": 0.3})
        state = mixture.properties(2_750.0, pressure_pa=2_000_000.0)
        self.assertAlmostEqual(
            mixture.temperature_from_enthalpy(
                state.enthalpy_j_kg, pressure_pa=2_000_000.0
            ),
            2_750.0,
            places=5,
        )
        self.assertAlmostEqual(
            mixture.temperature_from_entropy(
                state.entropy_j_kg_k, pressure_pa=2_000_000.0
            ),
            2_750.0,
            places=5,
        )


if __name__ == "__main__":
    unittest.main()
