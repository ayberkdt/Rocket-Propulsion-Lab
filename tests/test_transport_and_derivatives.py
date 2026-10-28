"""Transport properties and equilibrium derivatives against NASA CEA and published values.

Reference values in ``reference_data/cea_transport_points.json`` were produced
by NASA CEA 3.3.4 with ``transport=True`` (``scripts/generate_cea_transport_points.py``);
the file records the CEA version, the SHA-256 of the ``thermo.lib`` and
``trans.lib`` it used, and CEA's gas constant. No CEA install is needed here.

Tolerances and why (measured deviations in ``docs/TRANSPORT_AND_DERIVATIVES.md``):

* species fits: 1e-12 relative. Same trans.inp coefficients and the same
  equation (5.1) of RP-1311 Part I; only rounding differs.
* fixed small mixtures whose species all have trans.inp data: viscosity and
  frozen conductivity 1e-7. Quantities proportional to the gas constant
  (specific heats, reaction conductivity, Prandtl numbers) retain the original
  2e-5 allowance. The earlier 5.7e-6 SI/CEA scale difference is removed:
  the NASA model now uses CEA's R = 8314.51 J/(kmol K).
* rocket chamber/throat states: CEA estimates properties of HO2, H2O2, O3,
  HCO, COOH, ... from an undocumented correlation; this library omits them
  (at most 3.3e-4 of the mixture here). Viscosity 1e-4, frozen conductivity
  3e-4, equilibrium conductivity and Prandtl numbers 1e-3.
* equilibrium derivatives evaluated on CEA's own composition: gamma_s 1e-7
  (dimensionless, independent of R); specific heats and sound speed 2e-6
  on the common NASA basis (residual from CEA trace species absent from the
  bundled thermodynamic catalog).
* derivatives of the local provider's own equilibrium at CEA's T and p:
  gamma_s 1e-6, sound speed and specific heat 1e-5. No fixture or numerical
  acceptance tolerance is widened for the basis repair.
"""

from __future__ import annotations

import json
import unittest
from math import sqrt
from pathlib import Path

from rocket_propulsion.core.constants import (
    USSA1976_SEA_LEVEL_TEMPERATURE,
)
from rocket_propulsion.core.errors import DomainError, InputError
from rocket_propulsion.thermochemistry.equilibrium import (
    EquilibriumProblem,
    LocalEquilibriumProvider,
)
from rocket_propulsion.thermochemistry.equilibrium_derivatives import (
    derivatives_from_composition,
    equilibrium_derivatives,
    finite_difference_derivatives,
)
from rocket_propulsion.thermochemistry.rocket_equilibrium import (
    get_propellant_pair,
    reactants_for_mass_ratio,
)
from rocket_propulsion.thermochemistry.species import species_catalog
from rocket_propulsion.thermochemistry.transport import (
    DATA_FILE,
    binary_interaction,
    chapman_enskog_viscosity,
    equilibrium_transport_properties,
    get_species_transport,
    independent_reactions,
    mixture_transport_properties,
    modified_eucken_conductivity,
    neufeld_collision_integral,
    prandtl_number,
    rigid_sphere_interaction_viscosity,
    sutherland_viscosity,
)

DATA = Path(__file__).parent / "reference_data"
POINTS = json.loads((DATA / "cea_transport_points.json").read_text(encoding="utf-8"))
# CEA output units recorded in the fixture: millipoise, mW/(cm K), kJ/(kg K).
MILLIPOISE = 1.0e-4
MILLIWATT_PER_CM_K = 0.1
KJ = 1.0e3
NASA_R = 8.314510  # NASA/TP-2002-211556 Appendix A, independent of production alias.
R_RATIO = NASA_R * 1.0e3 / POINTS["reference"]["universal_gas_constant_j_kmol_k"]


def _relative(actual: float, expected: float) -> float:
    return abs(actual / expected - 1.0)


def _cea_fractions(record: dict) -> dict[str, float]:
    return record["cea"]["mole_fractions"]


class FixtureProvenanceTests(unittest.TestCase):
    def test_reference_file_records_cea_version_libraries_and_gas_constant(self) -> None:
        reference = POINTS["reference"]
        self.assertEqual(reference["version"], "3.3.4")
        self.assertEqual(len(reference["thermo_lib_sha256"]), 64)
        self.assertEqual(len(reference["trans_lib_sha256"]), 64)
        self.assertEqual(reference["universal_gas_constant_j_kmol_k"], 8314.51)
        self.assertEqual(len(POINTS["rocket"]), 18)

    def test_bundled_transport_data_is_pinned_to_trans_inp(self) -> None:
        document = json.loads(DATA_FILE.read_text(encoding="utf-8"))
        source = document["source"]
        self.assertEqual(source["commit"], "2f79a647e85737742fb08a127da00d1ea1d069da")
        self.assertEqual(
            source["sha256"], "379c514a7f7638371d8d9254a1d653cd078084a9371061aa0b45c5a5c3334d41"
        )
        self.assertEqual(source["license"], "Apache-2.0")
        self.assertEqual(len(document["species"]), 66)
        self.assertEqual(len(document["binary_interactions"]), 41)
        intervals = sum(
            len(record["viscosity"]) + len(record["thermal_conductivity"])
            for record in document["species"] + document["binary_interactions"]
        )
        self.assertEqual(intervals, 379)  # number of " V"/" C" lines in trans.inp


class SpeciesFitTests(unittest.TestCase):
    def test_pure_species_viscosity_and_conductivity_match_cea(self) -> None:
        for record in POINTS["pure_species"]:
            with self.subTest(record["id"]):
                data = get_species_transport(record["species"])
                temperature = record["cea"]["T"]
                self.assertLess(_relative(data.viscosity_pa_s(temperature),
                                          record["cea"]["viscosity"] * MILLIPOISE), 1.0e-12)
                self.assertLess(_relative(data.thermal_conductivity_w_m_k(temperature),
                                          record["cea"]["conductivity_fr"] * MILLIWATT_PER_CM_K),
                                1.0e-12)

    def test_species_fit_fails_closed_outside_its_range(self) -> None:
        hydroxyl = get_species_transport("OH")  # fits start at 1000 K
        with self.assertRaises(DomainError) as caught:
            hydroxyl.viscosity_pa_s(900.0)
        self.assertEqual(caught.exception.code, "transport_temperature_out_of_range")
        self.assertGreater(hydroxyl.viscosity_pa_s(900.0, allow_extrapolation=True), 0.0)
        with self.assertRaises(DomainError):
            get_species_transport("HO2")
        with self.assertRaises(InputError):
            hydroxyl.viscosity_pa_s(float("nan"))

    def test_binary_lookup_is_order_independent(self) -> None:
        self.assertIs(binary_interaction("H2O", "H2"), binary_interaction("h2", "h2o"))
        self.assertIsNone(binary_interaction("OH", "CO2"))


class MixtureRuleTests(unittest.TestCase):
    def test_fixed_mixtures_match_cea(self) -> None:
        for record in POINTS["mixtures"]:
            with self.subTest(record["id"]):
                cea = record["cea"]
                result = mixture_transport_properties(_cea_fractions(record), cea["T"])
                self.assertEqual(result.omitted_species, ())
                self.assertLess(_relative(result.viscosity_pa_s, cea["viscosity"] * MILLIPOISE),
                                1.0e-7)
                self.assertLess(_relative(result.frozen_thermal_conductivity_w_m_k,
                                          cea["conductivity_fr"] * MILLIWATT_PER_CM_K), 1.0e-7)
                self.assertLess(_relative(result.equilibrium_thermal_conductivity_w_m_k,
                                          cea["conductivity_eq"] * MILLIWATT_PER_CM_K), 2.0e-5)
                self.assertLess(_relative(result.frozen_cp_j_kg_k, cea["cp_fr"] * KJ), 2.0e-5)
                self.assertLess(_relative(result.frozen_prandtl, cea["Pr_fr"]), 2.0e-5)
                self.assertLess(_relative(result.equilibrium_prandtl, cea["Pr_eq"]), 2.0e-5)

    def test_reaction_terms_exist_only_when_species_can_react(self) -> None:
        counts = {
            record["id"]: mixture_transport_properties(_cea_fractions(record),
                                                       record["cea"]["T"]).reaction_count
            for record in POINTS["mixtures"]
        }
        self.assertEqual(counts["air-like-n2-o2-1000k"], 0)
        self.assertEqual(counts["h2-h2o-2000k"], 0)
        self.assertEqual(counts["co-co2-o2-3000k"], 1)
        self.assertEqual(counts["h-o-system-3500k"], 4)
        frozen = mixture_transport_properties(_cea_fractions(POINTS["mixtures"][3]), 3500.0,
                                              reacting=False)
        self.assertEqual(frozen.reaction_thermal_conductivity_w_m_k, 0.0)
        self.assertEqual(frozen.equilibrium_prandtl, frozen.frozen_prandtl)

    def test_rocket_chamber_and_throat_states_match_cea(self) -> None:
        for record in POINTS["rocket"]:
            with self.subTest(record["id"]):
                cea = record["cea"]
                result = mixture_transport_properties(_cea_fractions(record), cea["T"])
                self.assertLess(sum(value for _, value in result.omitted_species), 5.0e-4)
                self.assertLess(_relative(result.viscosity_pa_s, cea["viscosity"] * MILLIPOISE),
                                1.0e-4)
                self.assertLess(_relative(result.frozen_thermal_conductivity_w_m_k,
                                          cea["conductivity_fr"] * MILLIWATT_PER_CM_K), 3.0e-4)
                self.assertLess(_relative(result.equilibrium_thermal_conductivity_w_m_k,
                                          cea["conductivity_eq"] * MILLIWATT_PER_CM_K), 1.0e-3)
                self.assertLess(_relative(result.frozen_prandtl, cea["Pr_fr"]), 1.0e-3)
                self.assertLess(_relative(result.equilibrium_prandtl, cea["Pr_eq"]), 1.0e-3)

    def test_too_much_mixture_without_data_fails_closed(self) -> None:
        with self.assertRaises(DomainError) as caught:
            mixture_transport_properties({"HO2": 0.1, "N2": 0.9}, 1500.0)
        self.assertEqual(caught.exception.code, "transport_data_missing")
        with self.assertRaises(InputError):
            mixture_transport_properties({"N2": -1.0}, 1000.0)
        with self.assertRaises(InputError):
            mixture_transport_properties({"N2": 1.0}, 0.0)

    def test_omission_and_extrapolation_are_reported(self) -> None:
        result = mixture_transport_properties({"HO2": 1.0e-4, "N2": 1.0}, 1500.0)
        self.assertEqual([key for key, _ in result.omitted_species], ["HO2"])
        self.assertIn("transport_species_omitted", [w.code for w in result.warnings])
        with self.assertRaises(DomainError):
            mixture_transport_properties({"H2O": 0.9, "OH": 0.1}, 900.0)
        extrapolated = mixture_transport_properties({"H2O": 0.9, "OH": 0.1}, 900.0,
                                                    allow_extrapolation=True)
        self.assertIn("transport_fit_extrapolated", [w.code for w in extrapolated.warnings])

    def test_equilibrium_result_wrapper(self) -> None:
        pair = get_propellant_pair("LOX_LH2")
        result = LocalEquilibriumProvider().solve(EquilibriumProblem.tp(
            reactants_for_mass_ratio(pair, 5.5), temperature_k=3398.368485222262,
            pressure_pa=6894760.0,
        ))
        record = next(r for r in POINTS["rocket"] if r["id"] == "cea-lox_lh2-of-5.5-chamber")
        transport = equilibrium_transport_properties(result)
        self.assertLess(_relative(transport.viscosity_pa_s,
                                  record["cea"]["viscosity"] * MILLIPOISE), 1.0e-4)
        self.assertLess(_relative(transport.equilibrium_thermal_conductivity_w_m_k,
                                  record["cea"]["conductivity_eq"] * MILLIWATT_PER_CM_K), 1.0e-3)


class MixtureIdentityTests(unittest.TestCase):
    def test_rigid_sphere_interaction_is_symmetric_and_reproduces_equation_5_5(self) -> None:
        eta_i, eta_j, m_i, m_j = 4.2e-5, 2.1e-5, 0.028, 0.002
        forward = rigid_sphere_interaction_viscosity(eta_i, eta_j, m_i, m_j)
        backward = rigid_sphere_interaction_viscosity(eta_j, eta_i, m_j, m_i)
        self.assertAlmostEqual(forward / backward, 1.0, places=13)
        phi_57 = eta_i / forward * 2.0 * m_j / (m_i + m_j)
        phi_55 = 0.25 * (1.0 + sqrt(eta_i / eta_j) * (m_j / m_i) ** 0.25) ** 2 * sqrt(
            2.0 * m_j / (m_i + m_j))
        self.assertAlmostEqual(phi_57 / phi_55, 1.0, places=13)

    def test_pure_gas_mixture_rule_returns_the_species_value(self) -> None:
        result = mixture_transport_properties({"N2": 1.0}, 1000.0)
        data = get_species_transport("N2")
        self.assertEqual(result.viscosity_pa_s, data.viscosity_pa_s(1000.0))
        self.assertEqual(result.frozen_thermal_conductivity_w_m_k,
                         data.thermal_conductivity_w_m_k(1000.0))

    def test_reaction_basis_conserves_elements(self) -> None:
        formulas = [{"H": 2.0, "O": 1.0}, {"H": 2.0}, {"O": 2.0}, {"O": 1.0, "H": 1.0},
                    {"H": 1.0}, {"O": 1.0}]
        reactions = independent_reactions(formulas)
        self.assertEqual(len(reactions), 4)
        for alpha in reactions:
            for element in ("H", "O"):
                self.assertAlmostEqual(
                    sum(a * f.get(element, 0.0) for a, f in zip(alpha, formulas, strict=True)),
                    0.0, places=12)

    def test_prandtl_and_eucken_helpers(self) -> None:
        self.assertAlmostEqual(prandtl_number(1000.0, 2.0e-5, 0.025), 0.8)
        # Monatomic gas, cp = 5/2 R: modified Eucken reduces to 15/4 (R/M) mu.
        value = modified_eucken_conductivity(2.0e-5, 2.5 * NASA_R, 0.04)
        self.assertAlmostEqual(value, 3.75 * NASA_R / 0.04 * 2.0e-5, places=15)
        with self.assertRaises(InputError):
            prandtl_number(1000.0, 0.0, 0.025)


class FallbackModelTests(unittest.TestCase):
    def test_sutherland_reproduces_ussa1976_sea_level_viscosity(self) -> None:
        # U.S. Standard Atmosphere 1976, Table III, printed p. 101, Z = 0 m:
        # mu = 1.7894e-5 N s/m^2 at T0 = 288.15 K (5 significant digits).
        mu = sutherland_viscosity(USSA1976_SEA_LEVEL_TEMPERATURE.value)
        self.assertLess(_relative(mu, 1.7894e-5), 3.0e-5)

    def test_chapman_enskog_against_svehla_nitrogen_table(self) -> None:
        # Svehla, NASA TR R-132 (1962), Table III, printed p. 90, N2 with
        # sigma = 3.798 A, eps/k = 71.4 K, M = 28.02: eta x 1e6 poise.
        # Svehla used tabulated collision integrals; the Neufeld fit differs
        # by up to 0.8 % at T* = 28, hence the 1 % tolerance.
        table = {100.0: 72.4, 300.0: 177.7, 1000.0: 397.1, 2000.0: 624.0}
        for temperature, micropoise in table.items():
            with self.subTest(temperature=temperature):
                value = chapman_enskog_viscosity(
                    temperature, molar_mass_kg_mol=0.02802, collision_diameter_m=3.798e-10,
                    well_depth_k=71.4,
                )
                self.assertLess(_relative(value, micropoise * 1.0e-7), 1.0e-2)

    def test_collision_integral_rejects_invalid_input(self) -> None:
        self.assertGreater(neufeld_collision_integral(1.0), neufeld_collision_integral(10.0))
        with self.assertRaises(InputError):
            neufeld_collision_integral(0.0)


class DerivativeTests(unittest.TestCase):
    catalog = species_catalog()

    def _cea_composition(self, record: dict) -> dict[str, float]:
        return {key.upper(): value for key, value in _cea_fractions(record).items()
                if key.upper() in self.catalog}

    def test_derivatives_on_cea_composition_match_cea(self) -> None:
        for record in POINTS["rocket"] + POINTS["mixtures"][2:]:
            with self.subTest(record["id"]):
                cea = record["cea"]
                result = derivatives_from_composition(self._cea_composition(record), cea["T"],
                                                      cea["P"] * 1.0e5)
                self.assertLess(_relative(result.gamma_s, cea["gamma_s"]), 1.0e-7)
                scale = 1.0 / R_RATIO
                self.assertLess(_relative(result.equilibrium_cp_j_kg_k * scale, cea["cp_eq"] * KJ),
                                2.0e-6)
                self.assertLess(_relative(result.equilibrium_cv_j_kg_k * scale, cea["cv_eq"] * KJ),
                                2.0e-6)
                self.assertLess(_relative(result.frozen_cp_j_kg_k * scale, cea["cp_fr"] * KJ),
                                2.0e-6)
                self.assertLess(_relative(result.frozen_cv_j_kg_k * scale, cea["cv_fr"] * KJ),
                                2.0e-6)
                if "sonic_velocity" in cea:
                    self.assertLess(_relative(result.equilibrium_sound_speed_m_s / sqrt(R_RATIO),
                                              cea["sonic_velocity"]), 2.0e-6)

    def test_local_provider_derivatives_at_cea_states(self) -> None:
        provider = LocalEquilibriumProvider()
        for record in POINTS["rocket"]:
            with self.subTest(record["id"]):
                cea = record["cea"]
                pair = get_propellant_pair(record["pair"])
                result = provider.solve(EquilibriumProblem.tp(
                    reactants_for_mass_ratio(pair, record["oxidizer_fuel_ratio"]),
                    temperature_k=cea["T"], pressure_pa=cea["P"] * 1.0e5,
                ))
                derivatives = equilibrium_derivatives(result)
                self.assertLess(_relative(derivatives.gamma_s, cea["gamma_s"]), 1.0e-6)
                self.assertLess(_relative(derivatives.equilibrium_sound_speed_m_s / sqrt(R_RATIO),
                                          cea["sonic_velocity"]), 1.0e-5)
                self.assertLess(_relative(derivatives.equilibrium_cp_j_kg_k / R_RATIO,
                                          cea["cp_eq"] * KJ), 1.0e-5)

    def test_frozen_limit_without_possible_reactions(self) -> None:
        record = next(r for r in POINTS["pure_species"] if r["id"] == "pure-N2-1000k")
        result = derivatives_from_composition({"N2": 1.0}, 1000.0, 1.0e5)
        self.assertAlmostEqual(result.dlnv_dlnt_p, 1.0, places=12)
        self.assertAlmostEqual(result.dlnv_dlnp_t, -1.0, places=12)
        self.assertAlmostEqual(result.gamma_s, result.frozen_gamma, places=12)
        self.assertAlmostEqual(result.equilibrium_cp_j_kg_k, result.frozen_cp_j_kg_k, places=9)
        self.assertLess(_relative(result.gamma_s, record["cea"]["gamma_s"]), 1.0e-7)

    def test_bridgman_identities(self) -> None:
        record = next(r for r in POINTS["rocket"] if r["id"] == "cea-lox_ch4-of-3.4-chamber")
        cea = record["cea"]
        d = derivatives_from_composition(self._cea_composition(record), cea["T"], cea["P"] * 1.0e5)
        gas_constant = NASA_R / d.molecular_mass_kg_mol
        # RP-1311 eq. (2.68): (dln p/dln V)_s = -gamma_s.
        dlnp_dlnv = d.equilibrium_cp_j_kg_k / (
            d.equilibrium_cp_j_kg_k * d.dlnv_dlnp_t + gas_constant * d.dlnv_dlnt_p ** 2)
        self.assertAlmostEqual(-dlnp_dlnv / d.gamma_s, 1.0, places=12)
        # Shifting equilibrium raises cp and lowers gamma_s below the frozen ratio.
        self.assertGreater(d.equilibrium_cp_j_kg_k, d.frozen_cp_j_kg_k)
        self.assertLess(d.gamma_s, d.frozen_gamma)
        self.assertGreater(d.dlnv_dlnt_p, 1.0)

    def test_finite_difference_check_of_linearized_derivatives(self) -> None:
        provider = LocalEquilibriumProvider()
        for pair_key, ratio, temperature in (("LOX_CH4", 3.4, 3549.196144498271),
                                             ("LOX_LH2", 5.5, 3398.368485222262)):
            with self.subTest(pair_key):
                result = provider.solve(EquilibriumProblem.tp(
                    reactants_for_mass_ratio(get_propellant_pair(pair_key), ratio),
                    temperature_k=temperature, pressure_pa=6894760.0,
                ))
                analytic = equilibrium_derivatives(result)
                numeric = finite_difference_derivatives(result, provider, relative_step=1.0e-3)
                self.assertLess(_relative(numeric.equilibrium_cp_j_kg_k,
                                          analytic.equilibrium_cp_j_kg_k), 2.0e-5)
                self.assertLess(_relative(numeric.dlnv_dlnt_p, analytic.dlnv_dlnt_p), 2.0e-5)
                self.assertLess(_relative(numeric.dlnv_dlnp_t, analytic.dlnv_dlnp_t), 2.0e-5)
                self.assertLess(_relative(numeric.equilibrium_sound_speed_m_s,
                                          analytic.equilibrium_sound_speed_m_s), 2.0e-5)

    def test_invalid_inputs(self) -> None:
        with self.assertRaises(InputError):
            derivatives_from_composition({}, 3000.0, 1.0e6)
        with self.assertRaises(InputError):
            derivatives_from_composition({"H2O": 1.0}, 3000.0, -1.0)


if __name__ == "__main__":
    unittest.main()
