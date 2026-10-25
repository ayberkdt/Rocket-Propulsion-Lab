"""Finite-rate kinetics: stiff integrator, mechanism consistency and nozzle limits.

Reference values live in ``reference_data/finite_rate_kinetics.json``: the
Robertson and HIRES solutions of Hairer's stiff test set, the GRI-Mech 3.0
transcription spot check and NASA CEA 3.3.4 frozen and shifting-equilibrium
bounds for the nozzle case (``scripts/generate_finite_rate_cea_bounds.py``).
"""

import json
import unittest
from math import exp, log10
from pathlib import Path

from rocket_propulsion.core.constants import MOLAR_GAS_CONSTANT_J_MOL_K
from rocket_propulsion.core.errors import ConvergenceError, DomainError, InputError
from rocket_propulsion.thermochemistry import (
    EquilibriumProblem,
    LocalEquilibriumProvider,
    get_propellant_pair,
    get_species,
    reactants_for_mass_ratio,
    solve_rocket_performance,
)
from rocket_propulsion.thermochemistry.kinetics import (
    ArrheniusRate,
    Mechanism,
    Reaction,
    concentrations_from_mole_fractions,
    equilibrium_constants,
    evaluate_rates,
    production_rate_jacobian,
    solve_stiff,
)
from rocket_propulsion.thermochemistry.kinetics_mechanisms import (
    H2O2_SPECIES,
    gri30_h2o2_mechanism,
)
from rocket_propulsion.thermochemistry.nozzle_kinetics import (
    HyperbolicNozzle,
    _Expansion,
    bray_sudden_freezing,
    solve_nozzle_kinetics,
)

REFERENCE = json.loads(
    (Path(__file__).parent / "reference_data" / "finite_rate_kinetics.json").read_text(
        encoding="utf-8"))
TOLERANCE = REFERENCE["declared_tolerances"]
CASE = REFERENCE["nozzle_case"]
CEA = REFERENCE["cea_bounds"]["results"]


def _robertson(t, y):
    f1 = -0.04 * y[0] + 1.0e4 * y[1] * y[2]
    f3 = 3.0e7 * y[1] * y[1]
    return [f1, -f1 - f3, f3]


def _robertson_jacobian(t, y):
    return [[-0.04, 1.0e4 * y[2], 1.0e4 * y[1]],
            [0.04, -1.0e4 * y[2] - 6.0e7 * y[1], -1.0e4 * y[1]],
            [0.0, 6.0e7 * y[1], 0.0]]


def _hires(t, y):
    d = [0.0] * 8
    d[0] = -1.71 * y[0] + 0.43 * y[1] + 8.32 * y[2] + 0.0007
    d[1] = 1.71 * y[0] - 8.75 * y[1]
    d[2] = -10.03 * y[2] + 0.43 * y[3] + 0.035 * y[4]
    d[3] = 8.32 * y[1] + 1.71 * y[2] - 1.12 * y[3]
    d[4] = -1.745 * y[4] + 0.43 * y[5] + 0.43 * y[6]
    d[5] = -280.0 * y[5] * y[7] + 0.69 * y[3] + 1.71 * y[4] - 0.43 * y[5] + 0.69 * y[6]
    d[6] = 280.0 * y[5] * y[7] - 1.81 * y[6]
    d[7] = -d[6]
    return d


class StiffIntegratorTests(unittest.TestCase):
    def test_robertson_matches_reference_solution(self) -> None:
        data = REFERENCE["stiff_benchmarks"]["robertson"]
        rtol = data["solver"]["rtol"]
        solution = solve_stiff(_robertson, (0.0, data["output_times"][-1]),
                               data["initial_state"], rtol=rtol,
                               atol=data["solver"]["atol_over_rtol"] * rtol,
                               jacobian=_robertson_jacobian, autonomous=True,
                               t_eval=data["output_times"])
        for t, computed, expected in zip(solution.t, solution.y, data["reference"], strict=True):
            for component, (value, reference) in enumerate(zip(computed, expected, strict=True)):
                with self.subTest(t=t, component=component):
                    self.assertLessEqual(abs(value - reference) / abs(reference),
                                         data["declared_relative_tolerance"])

    def test_robertson_error_decreases_with_tolerance(self) -> None:
        data = REFERENCE["stiff_benchmarks"]["robertson"]
        errors = []
        for rtol in (1.0e-4, 1.0e-6, 1.0e-8):
            solution = solve_stiff(_robertson, (0.0, 1.0e11), data["initial_state"], rtol=rtol,
                                   atol=1.0e-6 * rtol, jacobian=_robertson_jacobian,
                                   autonomous=True, t_eval=data["output_times"])
            errors.append(max(abs(v - r) / abs(r) for y, ref in
                              zip(solution.y, data["reference"], strict=True)
                              for v, r in zip(y, ref, strict=True)))
        self.assertGreater(errors[0], errors[1])
        self.assertGreater(errors[1], errors[2])

    def test_hires_matches_reference_solution(self) -> None:
        data = REFERENCE["stiff_benchmarks"]["hires"]
        rtol = data["solver"]["rtol"]
        solution = solve_stiff(_hires, (0.0, data["output_times"][-1]), data["initial_state"],
                               rtol=rtol, atol=data["solver"]["atol_over_rtol"] * rtol,
                               autonomous=True, t_eval=data["output_times"])
        for computed, expected in zip(solution.y, data["reference"], strict=True):
            for value, reference in zip(computed, expected, strict=True):
                self.assertLessEqual(abs(value - reference) / abs(reference),
                                     data["declared_relative_tolerance"])

    def test_step_budget_and_inputs_fail_closed(self) -> None:
        with self.assertRaises(ConvergenceError) as context:
            solve_stiff(_robertson, (0.0, 1.0e11), [1.0, 0.0, 0.0], rtol=1.0e-8,
                        atol=1.0e-14, jacobian=_robertson_jacobian, autonomous=True,
                        max_steps=10)
        self.assertEqual(context.exception.code, "stiff_max_steps")
        with self.assertRaises(InputError):
            solve_stiff(_robertson, (0.0, 0.0), [1.0, 0.0, 0.0])
        with self.assertRaises(InputError):
            solve_stiff(_robertson, (0.0, 1.0), [1.0, 0.0, 0.0], atol=[1.0e-8])


class MechanismTests(unittest.TestCase):
    def setUp(self) -> None:
        self.mechanism = gri30_h2o2_mechanism()

    def test_transcription_spot_check_and_units(self) -> None:
        data = REFERENCE["mechanism"]
        self.assertEqual(len(self.mechanism.reactions), data["reaction_count"])
        spot = data["spot_check"]
        reaction = next(r for r in self.mechanism.reactions
                        if r.source.endswith(f"line {spot['line']}"))
        self.assertEqual(reaction.equation, spot["equation"])
        self.assertAlmostEqual(reaction.rate.pre_exponential_si, spot["si"]["A_m3_mol_s"],
                               delta=1.0e-15)
        self.assertEqual(reaction.rate.temperature_exponent, spot["si"]["b"])
        self.assertAlmostEqual(reaction.rate.activation_energy_j_mol, spot["si"]["E_J_mol"],
                               delta=1.0e-9)
        self.assertTrue(all("grimech30.dat line" in r.source for r in self.mechanism.reactions))

    def test_unbalanced_reaction_is_rejected(self) -> None:
        with self.assertRaises(DomainError):
            Mechanism(name="bad", species=("H", "H2", "O"), reactions=(Reaction(
                equation="H+H<=>O", reactants=(("H", 2),), products=(("O", 1),),
                rate=ArrheniusRate(1.0, 0.0, 0.0)),))

    def test_forward_and_reverse_rates_balance_at_equilibrium(self) -> None:
        provider = LocalEquilibriumProvider()
        for temperature, pressure in ((3500.0, 7.0e6), (2500.0, 1.0e6), (1500.0, 1.0e5)):
            result = provider.solve(EquilibriumProblem.tp(
                {"H2": 2.0, "O2": 0.8}, temperature_k=temperature, pressure_pa=pressure,
                candidate_species=H2O2_SPECIES))
            fractions = {item.species: item.mole_fraction for item in result.species}
            concentrations = concentrations_from_mole_fractions(self.mechanism, temperature,
                                                                pressure, fractions)
            rates = evaluate_rates(self.mechanism, temperature, concentrations)
            for j, (forward, reverse) in enumerate(zip(rates.forward_rates, rates.reverse_rates,
                                                       strict=True)):
                with self.subTest(temperature=temperature, reaction=j):
                    self.assertLessEqual(abs(forward - reverse) / max(forward, reverse),
                                         TOLERANCE["detailed_balance_relative"])

    def test_equilibrium_constant_matches_species_gibbs_energies(self) -> None:
        temperature = 2800.0
        constants = equilibrium_constants(self.mechanism, temperature)
        r = MOLAR_GAS_CONSTANT_J_MOL_K
        for reaction, kc in zip(self.mechanism.reactions, constants, strict=True):
            gibbs = 0.0
            for key, change in reaction.net_stoichiometry.items():
                props = get_species(key).properties(temperature)
                gibbs += change * (props.enthalpy_molar_j_mol
                                   - temperature * props.entropy_molar_j_mol_k)
            delta = sum(reaction.net_stoichiometry.values())
            expected = (1.0e5 / (r * temperature)) ** delta * exp(-gibbs / (r * temperature))
            self.assertAlmostEqual(kc / expected, 1.0, delta=1.0e-10)

    def test_production_rates_conserve_elements_and_jacobian_matches(self) -> None:
        temperature, pressure = 2600.0, 2.0e6
        fractions = dict(zip(H2O2_SPECIES, (0.25, 0.02, 0.003, 0.004, 0.03, 0.69, 2.0e-5,
                                            1.0e-5), strict=True))
        concentrations = concentrations_from_mole_fractions(self.mechanism, temperature,
                                                            pressure, fractions)
        omega = evaluate_rates(self.mechanism, temperature, concentrations).net_production_rates
        scale = max(abs(value) for value in omega)
        for row in self.mechanism.element_matrix():
            self.assertLessEqual(abs(sum(a * w for a, w in zip(row, omega, strict=True))),
                                 TOLERANCE["element_conservation_relative"] * scale)
        jacobian = production_rate_jacobian(self.mechanism, temperature, concentrations)
        for k in range(len(concentrations)):
            step = 1.0e-6 * concentrations[k]
            plus = list(concentrations)
            minus = list(concentrations)
            plus[k] += step
            minus[k] -= step
            up = evaluate_rates(self.mechanism, temperature, plus).net_production_rates
            down = evaluate_rates(self.mechanism, temperature, minus).net_production_rates
            for i in range(len(concentrations)):
                row_scale = max(abs(v) for v in jacobian[i]) or 1.0
                with self.subTest(row=i, column=k):
                    self.assertLessEqual(abs((up[i] - down[i]) / (2.0 * step) - jacobian[i][k])
                                         / row_scale,
                                         TOLERANCE["rate_jacobian_vs_central_difference"])

    def test_troe_falloff_follows_documented_form(self) -> None:
        falloff = next(r for r in self.mechanism.reactions if r.falloff is not None)
        index = self.mechanism.reactions.index(falloff)
        troe = falloff.falloff
        for temperature, water in ((1200.0, 10.0), (2500.0, 300.0), (3500.0, 3000.0)):
            concentrations = [0.0] * len(H2O2_SPECIES)
            concentrations[H2O2_SPECIES.index("H2O")] = water
            k = evaluate_rates(self.mechanism, temperature,
                               concentrations).forward_rate_constants[index]
            # Independent evaluation of the Troe form (Cantera "Troe Falloff
            # Function"; MFiX eqs. (4.18)-(4.20)) with [M] = 6 [H2O] (line 155).
            k_inf = falloff.rate(temperature)
            reduced = troe.low_pressure(temperature) * 6.0 * water / k_inf
            center = ((1.0 - troe.alpha) * exp(-temperature / troe.t3_k)
                      + troe.alpha * exp(-temperature / troe.t1_k) + exp(-troe.t2_k / temperature))
            c = -0.4 - 0.67 * log10(center)
            n = 0.75 - 1.27 * log10(center)
            f1 = (log10(reduced) + c) / (n - 0.14 * (log10(reduced) + c))
            expected = k_inf * reduced / (1.0 + reduced) * 10.0 ** (log10(center) / (1.0 + f1 ** 2))
            with self.subTest(temperature=temperature):
                self.assertAlmostEqual(k / expected, 1.0, delta=1.0e-12)


def _case():
    reactants = reactants_for_mass_ratio(get_propellant_pair(CASE["propellant_pair"]),
                                         CASE["oxidizer_fuel_ratio"])
    nozzle = HyperbolicNozzle(**CASE["nozzle"])
    return reactants, nozzle, CASE["chamber_pressure_pa"]


class NozzleKineticsTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.reactants, cls.nozzle, cls.pc = _case()
        cls.repository = {
            chemistry: solve_rocket_performance(
                reactants=cls.reactants, chamber_pressure_pa=cls.pc,
                area_ratios=(cls.nozzle.exit_area_ratio,), chemistry=chemistry)
            .exits[0].vacuum_specific_impulse_s
            for chemistry in ("frozen", "equilibrium")}
        cls.baseline = solve_nozzle_kinetics(reactants=cls.reactants, chamber_pressure_pa=cls.pc,
                                             nozzle=cls.nozzle)

    def test_baseline_lies_between_frozen_and_shifting(self) -> None:
        isp = self.baseline.vacuum_specific_impulse_s
        self.assertLess(self.repository["frozen"], isp)
        self.assertLess(isp, self.repository["equilibrium"])
        self.assertLess(CEA["frozen_at_chamber"]["vacuum_specific_impulse_s"], isp)
        self.assertLess(isp, CEA["shifting_equilibrium"]["vacuum_specific_impulse_s"])
        self.assertLess(self.baseline.maximum_element_residual, 1.0e-8)

    def test_zero_rates_reproduce_frozen_chamber_expansion(self) -> None:
        frozen = solve_nozzle_kinetics(reactants=self.reactants, chamber_pressure_pa=self.pc,
                                       nozzle=self.nozzle, rate_multiplier=0.0, start="chamber")
        isp = frozen.vacuum_specific_impulse_s
        self.assertLessEqual(abs(isp / self.repository["frozen"] - 1.0),
                             TOLERANCE["zero_rate_vs_repository_frozen_isp_relative"])
        cea = CEA["frozen_at_chamber"]
        self.assertLessEqual(abs(isp / cea["vacuum_specific_impulse_s"] - 1.0),
                             TOLERANCE["zero_rate_vs_cea_frozen_isp_relative"])
        self.assertLessEqual(abs(frozen.characteristic_velocity_m_s
                                 / cea["characteristic_velocity_m_s"] - 1.0),
                             TOLERANCE["zero_rate_vs_cea_frozen_isp_relative"])
        self.assertEqual(frozen.maximum_element_residual, 0.0)

    def test_fast_rates_converge_to_shifting_equilibrium(self) -> None:
        equilibrium = self.repository["equilibrium"]
        gaps = [equilibrium - self.baseline.vacuum_specific_impulse_s]
        for multiplier in (1.0e2, 1.0e4):
            result = solve_nozzle_kinetics(reactants=self.reactants, chamber_pressure_pa=self.pc,
                                           nozzle=self.nozzle, rate_multiplier=multiplier)
            gaps.append(equilibrium - result.vacuum_specific_impulse_s)
        self.assertGreater(gaps[0], gaps[1])
        self.assertGreater(gaps[1], gaps[2])
        self.assertLessEqual(abs(gaps[2]), TOLERANCE["fast_rate_1e4_vs_repository_shifting_isp_s"])

    def test_zero_rates_after_equilibrium_start_freeze_downstream_of_throat(self) -> None:
        result = solve_nozzle_kinetics(reactants=self.reactants, chamber_pressure_pa=self.pc,
                                       nozzle=self.nozzle, rate_multiplier=0.0)
        isp = result.vacuum_specific_impulse_s
        self.assertLess(CEA["frozen_at_throat"]["vacuum_specific_impulse_s"], isp)
        self.assertLess(isp, self.baseline.vacuum_specific_impulse_s)

    def test_tolerance_and_start_station_convergence(self) -> None:
        coarse = solve_nozzle_kinetics(reactants=self.reactants, chamber_pressure_pa=self.pc,
                                       nozzle=self.nozzle, rtol=1.0e-5)
        early = solve_nozzle_kinetics(reactants=self.reactants, chamber_pressure_pa=self.pc,
                                      nozzle=self.nozzle, start_area_ratio=1.05)
        isp = self.baseline.vacuum_specific_impulse_s
        self.assertLessEqual(abs(coarse.vacuum_specific_impulse_s - isp),
                             TOLERANCE["rtol_1e-5_vs_1e-6_isp_s"])
        self.assertLessEqual(abs(early.vacuum_specific_impulse_s - isp),
                             TOLERANCE["start_area_ratio_1.05_vs_1.1_isp_s"])

    def test_radicals_recombine_along_the_nozzle(self) -> None:
        start = dict(self.baseline.stations[0].mole_fractions)
        exit_ = dict(self.baseline.exit.mole_fractions)
        for key in ("H", "O", "OH"):
            self.assertLess(exit_[key], start[key])
        self.assertGreater(exit_["H2O"], start["H2O"])
        temperatures = [station.temperature_k for station in self.baseline.stations]
        self.assertTrue(all(b < a for a, b in zip(temperatures, temperatures[1:], strict=False)))

    def test_bray_estimate_brackets(self) -> None:
        early = bray_sudden_freezing(reactants=self.reactants, chamber_pressure_pa=self.pc,
                                     nozzle=self.nozzle, freezing_constant=0.42)
        late = bray_sudden_freezing(reactants=self.reactants, chamber_pressure_pa=self.pc,
                                    nozzle=self.nozzle, freezing_constant=1.25)
        default = bray_sudden_freezing(reactants=self.reactants, chamber_pressure_pa=self.pc,
                                       nozzle=self.nozzle)
        self.assertLess(early.freezing_area_ratio, default.freezing_area_ratio)
        self.assertLess(default.freezing_area_ratio, late.freezing_area_ratio)
        self.assertLess(early.vacuum_specific_impulse_s, late.vacuum_specific_impulse_s)
        for estimate in (early, default, late):
            self.assertLess(self.repository["frozen"], estimate.vacuum_specific_impulse_s)
            self.assertLess(estimate.vacuum_specific_impulse_s,
                            estimate.equilibrium_vacuum_specific_impulse_s)
        self.assertLessEqual(abs(default.vacuum_specific_impulse_s
                                 - self.baseline.vacuum_specific_impulse_s),
                             TOLERANCE["bray_K1_vs_finite_rate_isp_s"])
        self.assertAlmostEqual(early.equilibrium_vacuum_specific_impulse_s,
                               CEA["shifting_equilibrium"]["vacuum_specific_impulse_s"],
                               delta=0.01)

    def test_analytic_nozzle_jacobian_matches_central_differences(self) -> None:
        expansion = _Expansion(self.reactants, self.pc, gri30_h2o2_mechanism(),
                               LocalEquilibriumProvider())
        size = expansion.size
        y = list(expansion.sigma0) + [expansion.s0]
        flux, branch = 2500.0, "sup"
        analytic = expansion.jacobian(y[:size], y[size], flux, branch)
        for k in range(size + 1):
            step = 1.0e-5 * (abs(y[k]) if k < size else 1.0)
            plus, minus = list(y), list(y)
            plus[k] += step
            minus[k] -= step
            expansion.last_temperature.clear()
            up = expansion.derivatives(plus[:size], plus[size], flux, branch)[0]
            expansion.last_temperature.clear()
            down = expansion.derivatives(minus[:size], minus[size], flux, branch)[0]
            for i in range(size + 1):
                row_scale = max(abs(v) for v in analytic[i]) or 1.0
                with self.subTest(row=i, column=k):
                    self.assertLessEqual(
                        abs((up[i] - down[i]) / (2.0 * step) - analytic[i][k]) / row_scale,
                        TOLERANCE["nozzle_jacobian_vs_central_difference"])

    def test_inputs_fail_closed(self) -> None:
        with self.assertRaises(InputError):
            solve_nozzle_kinetics(reactants={"CH4": 1.0, "O2": 2.0}, chamber_pressure_pa=self.pc,
                                  nozzle=self.nozzle)
        with self.assertRaises(InputError):
            solve_nozzle_kinetics(reactants=self.reactants, chamber_pressure_pa=self.pc,
                                  nozzle=self.nozzle, start="throat")
        with self.assertRaises(DomainError):
            HyperbolicNozzle(0.05, 3.0, 0.5, 30.0, 15.0)


if __name__ == "__main__":
    unittest.main()
