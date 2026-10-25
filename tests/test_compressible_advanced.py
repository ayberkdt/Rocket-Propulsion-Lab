"""Reference and consistency tests for the advanced compressible-flow modules.

Every reference number is loaded from ``reference_data/compressible_advanced.json``,
where each entry carries its source locator (report table/chart or textbook
example and page) and the tolerance basis.
"""

from __future__ import annotations

import json
import unittest
from math import asin, cos, degrees, pi, radians, sin, sqrt, tan
from pathlib import Path

from rocket_propulsion.compressible.cone_flow import (
    cone_flow_from_shock_angle,
    maximum_cone_half_angle,
    solve_cone_flow,
)
from rocket_propulsion.compressible.hypersonic_inviscid import (
    hypersonic_shock_limits,
    hypersonic_similarity_parameter,
    limiting_wave_angle_ratio,
    mach_independent_wedge,
    modified_newtonian_pressure_coefficient,
    newtonian_flat_plate,
    newtonian_maximum_lift,
    newtonian_maximum_lift_to_drag,
    newtonian_pressure_coefficient,
    small_disturbance_wedge,
    stagnation_pressure_coefficient,
)
from rocket_propulsion.compressible.isentropic import calculate_isentropic
from rocket_propulsion.compressible.linearized_flow import (
    ackeret_flat_plate,
    ackeret_pressure_coefficient,
    busemann_coefficients,
    busemann_pressure_coefficient,
    corrected_pressure_coefficient,
    critical_mach_number,
    critical_pressure_coefficient,
    karman_tsien,
    laitone,
    linearized_section_coefficients,
    prandtl_glauert,
    subsonic_validity_warnings,
)
from rocket_propulsion.compressible.normal_shock import calculate_normal_shock
from rocket_propulsion.compressible.oblique_shock import calculate_oblique_shock
from rocket_propulsion.compressible.shock_expansion import (
    AirfoilSection,
    diamond_airfoil,
    flat_plate,
    integrate_section_forces,
    parabolic_biconvex,
    shock_expansion_coefficients,
    wedge_section,
)
from rocket_propulsion.compressible.thermally_perfect import (
    HarmonicVibratorGas,
    PolynomialGas,
    thermally_perfect_isentropic,
    thermally_perfect_normal_shock,
    thermally_perfect_oblique_shock,
)
from rocket_propulsion.core.errors import DomainError
from rocket_propulsion.thermochemistry.species import get_species

REFERENCE = json.loads(
    (Path(__file__).parent / "reference_data" / "compressible_advanced.json").read_text(
        encoding="utf-8"
    )
)
RANKINE = 5.0 / 9.0  # 1 degR = 5/9 K (unit definition)
GAS = HarmonicVibratorGas()


def _section(entry: dict) -> AirfoilSection:
    if entry["section"] == "flat_plate":
        return flat_plate()
    return wedge_section(entry["half_angle_deg"])


class ThermallyPerfectGasTests(unittest.TestCase):
    def test_gamma_table_naca1135_p627(self) -> None:
        table = REFERENCE["thermally_perfect"]["gamma_table"]
        for temperature_r, printed in table["points_degR_gamma"]:
            with self.subTest(T_degR=temperature_r):
                self.assertAlmostEqual(GAS.gamma(temperature_r * RANKINE), printed,
                                       delta=table["tolerance_abs"])

    def test_isentropic_charts_11_12_14(self) -> None:
        for entry in REFERENCE["thermally_perfect"]["isentropic_charts"]:
            with self.subTest(chart=entry["chart"], Tt=entry["total_temperature_degR"]):
                result = thermally_perfect_isentropic(
                    entry["mach"], entry["total_temperature_degR"] * RANKINE, GAS
                )
                perfect = calculate_isentropic(entry["mach"], GAS.gamma_perfect)
                perfect_value = {
                    "density_ratio": perfect.density_total_ratio,
                    "pressure_ratio": perfect.pressure_total_ratio,
                    "area_ratio": perfect.area_critical_ratio,
                }[entry["quantity"]]
                ratio = getattr(result, entry["quantity"]) / perfect_value
                self.assertAlmostEqual(ratio, entry["value"], delta=entry["tolerance_abs"])

    def test_normal_shock_charts_15_to_20(self) -> None:
        for entry in REFERENCE["thermally_perfect"]["normal_shock_charts"]:
            with self.subTest(chart=entry["chart"], Tt=entry["total_temperature_degR"]):
                upstream = thermally_perfect_isentropic(
                    entry["mach"], entry["total_temperature_degR"] * RANKINE, GAS
                )
                result = thermally_perfect_normal_shock(entry["mach"], upstream.temperature_k, GAS)
                perfect = calculate_normal_shock(entry["mach"], GAS.gamma_perfect)
                ratio = getattr(result, entry["quantity"]) / getattr(perfect, entry["quantity"])
                self.assertAlmostEqual(ratio, entry["value"], delta=entry["tolerance_abs"])

    def test_low_temperature_limit_is_calorically_perfect(self) -> None:
        # Theta/T = 30: the vibrational mode is frozen to within 1e-9.
        total = GAS.vibrational_temperature_k / 30.0
        for mach in (0.5, 2.0, 5.0):
            with self.subTest(mach=mach):
                result = thermally_perfect_isentropic(mach, total, GAS)
                perfect = calculate_isentropic(mach, 1.4)
                self.assertAlmostEqual(result.pressure_ratio / perfect.pressure_total_ratio, 1.0,
                                       delta=1e-8)
                self.assertAlmostEqual(result.area_ratio / perfect.area_critical_ratio, 1.0,
                                       delta=1e-8)
        shock = thermally_perfect_normal_shock(3.0, total / 2.8, GAS)
        perfect = calculate_normal_shock(3.0, 1.4)
        self.assertAlmostEqual(shock.pressure_ratio / perfect.pressure_ratio, 1.0, delta=1e-8)
        self.assertAlmostEqual(shock.stagnation_pressure_ratio / perfect.stagnation_pressure_ratio,
                               1.0, delta=1e-8)

    def test_oblique_shock_reduces_to_perfect_gas_and_to_normal_component(self) -> None:
        frozen = HarmonicVibratorGas(vibrational_temperature_k=1.0e6)
        result = thermally_perfect_oblique_shock(3.0, 15.0, 300.0, frozen)
        perfect = calculate_oblique_shock(upstream_mach=3.0, turning_angle_deg=15.0)
        self.assertAlmostEqual(result.wave_angle_deg, perfect.wave_angle_deg, places=7)
        self.assertAlmostEqual(result.downstream_mach, perfect.downstream_mach, places=7)
        self.assertAlmostEqual(result.maximum_turning_angle_deg,
                               perfect.maximum_turning_angle_deg, places=5)
        hot = thermally_perfect_oblique_shock(6.0, 20.0, 600.0, GAS)
        normal = thermally_perfect_normal_shock(
            6.0 * sin(radians(hot.wave_angle_deg)), 600.0, GAS
        )
        self.assertAlmostEqual(hot.pressure_ratio, normal.pressure_ratio, places=10)
        # Vibrational excitation softens the gas: denser, weaker-angled shock.
        cold = calculate_oblique_shock(upstream_mach=6.0, turning_angle_deg=20.0)
        self.assertLess(hot.wave_angle_deg, cold.wave_angle_deg)
        self.assertGreater(hot.density_ratio, cold.density_ratio)

    def test_nasa_polynomial_nitrogen_near_room_temperature(self) -> None:
        gas = PolynomialGas(get_species("N2"))
        gamma = gas.gamma(300.0)
        self.assertAlmostEqual(gamma, 1.40, delta=0.002)
        shock = thermally_perfect_normal_shock(1.3, 300.0, gas)
        perfect = calculate_normal_shock(1.3, gamma)
        self.assertAlmostEqual(shock.pressure_ratio / perfect.pressure_ratio, 1.0, delta=1e-3)
        flow = thermally_perfect_isentropic(2.0, 600.0, gas)
        self.assertGreater(flow.temperature_k, gas.minimum_temperature_k)
        with self.assertRaises(DomainError):
            thermally_perfect_isentropic(10.0, 600.0, gas)  # static T below the fit range

    def test_inputs_fail_closed(self) -> None:
        with self.assertRaises(DomainError):
            thermally_perfect_normal_shock(0.9, 300.0)
        with self.assertRaises(DomainError):
            thermally_perfect_isentropic(float("nan"), 300.0)
        with self.assertRaises(DomainError):
            thermally_perfect_oblique_shock(2.0, 40.0, 300.0)
        with self.assertRaises(DomainError):
            HarmonicVibratorGas(gamma_perfect=1.0)


class ConeFlowTests(unittest.TestCase):
    def test_sims_sp3004_tables_7_9_11(self) -> None:
        data = REFERENCE["cone_flow"]
        tolerance = data["tolerance"]
        for point in data["points"]:
            with self.subTest(cone=point["cone_deg"], mach=point["mach"]):
                result = solve_cone_flow(point["mach"], point["cone_deg"])
                self.assertAlmostEqual(radians(result.shock_angle_deg), point["shock_angle_rad"],
                                       delta=tolerance["shock_angle_rad_abs"])
                self.assertAlmostEqual(result.surface_pressure_coefficient,
                                       point["pressure_coefficient"],
                                       delta=tolerance["pressure_coefficient_abs"])
                self.assertAlmostEqual(result.surface_mach / point["surface_mach"], 1.0,
                                       delta=tolerance["surface_mach_rel"])

    def test_anderson_figure_9_16_and_relieving_effect(self) -> None:
        entry = REFERENCE["cone_flow"]["anderson_figure_9_16"]
        cone = solve_cone_flow(entry["mach"], entry["cone_deg"])
        self.assertAlmostEqual(cone.shock_angle_deg, entry["shock_angle_deg"],
                               delta=entry["tolerance_deg"])
        wedge = calculate_oblique_shock(upstream_mach=2.0, turning_angle_deg=20.0)
        self.assertLess(cone.shock_angle_deg, wedge.wave_angle_deg)
        self.assertLess(cone.surface_pressure_ratio, wedge.pressure_ratio)
        self.assertGreater(cone.surface_mach, wedge.downstream_mach)

    def test_post_shock_state_matches_wedge_relation(self) -> None:
        result = cone_flow_from_shock_angle(4.0, 25.0)
        wedge = calculate_oblique_shock(upstream_mach=4.0,
                                        turning_angle_deg=result.post_shock_deflection_deg)
        self.assertAlmostEqual(wedge.wave_angle_deg, 25.0, places=8)
        self.assertAlmostEqual(wedge.downstream_mach, result.post_shock_mach, places=9)
        self.assertGreater(result.cone_half_angle_deg, result.post_shock_deflection_deg)

    def test_vanishing_cone_gives_mach_cone(self) -> None:
        result = solve_cone_flow(3.0, 0.05)
        self.assertAlmostEqual(result.shock_angle_deg, degrees(asin(1.0 / 3.0)), delta=0.05)
        self.assertAlmostEqual(result.surface_pressure_coefficient, 0.0, delta=1e-4)

    def test_detachment(self) -> None:
        maximum, _ = maximum_cone_half_angle(2.0)
        self.assertGreater(maximum, 40.0)
        with self.assertRaises(DomainError):
            solve_cone_flow(2.0, maximum + 0.5)


class ShockExpansionTests(unittest.TestCase):
    def test_anderson_examples_replay_printed_pressures(self) -> None:
        for entry in REFERENCE["shock_expansion"]["examples"]:
            with self.subTest(example=entry["id"]):
                section = _section(entry)
                dynamic = 0.5 * 1.4 * entry["mach"] ** 2
                cps = {surface: [(p - 1.0) / dynamic for p in values]
                       for surface, values in entry["replay_panel_pressure_ratios"].items()}
                lift, drag, cn, _, _ = integrate_section_forces(
                    section, cps, radians(entry["alpha_deg"]),
                    (entry["base_pressure_ratio"] - 1.0) / dynamic, section.upper[0])
                printed = entry["printed"]
                tol = entry["replay_tolerance_abs"]
                if "lift_coefficient" in printed:
                    self.assertAlmostEqual(lift, printed["lift_coefficient"], delta=tol)
                if "drag_coefficient" in printed:
                    self.assertAlmostEqual(drag, printed["drag_coefficient"], delta=tol)
                if "normal_force_coefficient" in printed:
                    self.assertAlmostEqual(cn, printed["normal_force_coefficient"], delta=tol)

    def test_anderson_examples_exact_theory(self) -> None:
        for entry in REFERENCE["shock_expansion"]["examples"]:
            if entry["exact_tolerance_rel"] is None:
                continue
            with self.subTest(example=entry["id"]):
                result = shock_expansion_coefficients(
                    _section(entry), entry["mach"], entry["alpha_deg"],
                    base_pressure_ratio=entry["base_pressure_ratio"])
                rel = entry["exact_tolerance_rel"]
                for key in ("lift_coefficient", "drag_coefficient", "normal_force_coefficient"):
                    if key in entry["printed"]:
                        value = getattr(result, key)
                        self.assertAlmostEqual(value / entry["printed"][key], 1.0, delta=rel,
                                               msg=key)

    def test_example_14_1_pressure_coefficients_and_lift_to_drag(self) -> None:
        printed = REFERENCE["shock_expansion"]["examples"][4]["printed"]
        result = shock_expansion_coefficients(flat_plate(), 8.0, 15.0)
        upper, lower = result.panels
        self.assertAlmostEqual(upper.pressure_coefficient,
                               printed["upper_pressure_coefficient"], delta=0.0005)
        self.assertAlmostEqual(lower.pressure_coefficient / printed["lower_pressure_coefficient"],
                               1.0, delta=0.02)
        self.assertAlmostEqual(result.lift_coefficient / result.drag_coefficient,
                               printed["lift_to_drag_ratio"], delta=0.005)

    def test_example_9_12_upper_surface_expansion(self) -> None:
        printed = REFERENCE["shock_expansion"]["examples"][1]["printed"]
        upper = shock_expansion_coefficients(flat_plate(), 3.0, 5.0).panels[0]
        self.assertAlmostEqual(upper.mach, printed["upper_mach"], delta=0.005)
        self.assertAlmostEqual(upper.pressure_ratio, printed["upper_pressure_ratio"], delta=0.0005)

    def test_ackeret_limit_for_small_angles(self) -> None:
        for mach in (1.5, 2.0, 3.0):
            beta = sqrt(mach * mach - 1.0)
            alpha = radians(0.2)
            exact = shock_expansion_coefficients(flat_plate(), mach, 0.2)
            self.assertAlmostEqual(exact.lift_coefficient / (4.0 * alpha / beta), 1.0, delta=2e-3)
            self.assertAlmostEqual(exact.drag_coefficient / (4.0 * alpha * alpha / beta), 1.0,
                                   delta=2e-3)
            epsilon = radians(0.2)
            diamond = shock_expansion_coefficients(diamond_airfoil(0.2), mach, 0.0)
            thickness = 2.0 * 0.5 * tan(epsilon)
            self.assertAlmostEqual(diamond.drag_coefficient / (4.0 * thickness**2 / beta), 1.0,
                                   delta=2e-3)

    def test_flat_plate_centre_of_pressure_and_geometry(self) -> None:
        result = shock_expansion_coefficients(flat_plate(), 4.0, 8.0)
        self.assertAlmostEqual(result.moment_coefficient, -0.5 * result.normal_force_coefficient,
                               places=12)
        self.assertAlmostEqual(result.axial_force_coefficient, 0.0, places=14)
        symmetric = shock_expansion_coefficients(parabolic_biconvex(0.06), 2.5, 0.0)
        self.assertAlmostEqual(symmetric.lift_coefficient, 0.0, places=12)
        with self.assertRaises(DomainError):
            shock_expansion_coefficients(wedge_section(30.0), 1.5, 0.0)  # detached
        with self.assertRaises(DomainError):
            AirfoilSection(((0.0, 0.0), (0.0, 1.0)), ((0.0, 0.0), (1.0, 0.0)))


class LinearizedFlowTests(unittest.TestCase):
    def test_anderson_chapter_11_examples(self) -> None:
        data = REFERENCE["linearized"]
        e1 = data["example_11_1"]
        self.assertAlmostEqual(prandtl_glauert(e1["cp0"], e1["mach"]), e1["cp"],
                               delta=e1["tolerance_abs"])
        e2 = data["example_11_2"]
        self.assertAlmostEqual(prandtl_glauert(2.0 * pi, e2["mach"]),
                               e2["lift_slope_per_rad"], delta=e2["tolerance_abs"])
        for mach, printed in data["critical_cp_table"]["points"]:
            digits = len(str(printed).split(".")[1]) if "." in str(printed) else 0
            with self.subTest(mach=mach):
                self.assertAlmostEqual(critical_pressure_coefficient(mach), printed,
                                       delta=0.5 * 10.0 ** -max(digits, 2) + 1e-12)
        for mach, printed in data["minimum_cp_table"]["points"]:
            with self.subTest(mach=mach):
                self.assertAlmostEqual(prandtl_glauert(-0.43, mach), printed, delta=0.0005 + 1e-9)
        table = data["critical_mach_iteration_table"]
        for mach, left, right in table["points"]:
            with self.subTest(mach=mach):
                self.assertAlmostEqual(prandtl_glauert(-0.43, mach), left,
                                       delta=table["tolerance_abs"])
                self.assertAlmostEqual(critical_pressure_coefficient(mach), right,
                                       delta=table["tolerance_abs"])
        naca = data["critical_mach_naca0012"]
        result = critical_mach_number(naca["cp0_min"])
        self.assertAlmostEqual(result.critical_mach, naca["critical_mach"],
                               delta=naca["tolerance_abs"])
        for mach, printed in data["critical_cp_naca0012_data"]["points"]:
            self.assertAlmostEqual(critical_pressure_coefficient(mach), printed, delta=0.0005)

    def test_karman_tsien_and_laitone_reduce_to_prandtl_glauert(self) -> None:
        for mach in (0.3, 0.6, 0.75):
            for cp0 in (-1e-4, 1e-4):
                pg = prandtl_glauert(cp0, mach)
                self.assertAlmostEqual(karman_tsien(cp0, mach) / pg, 1.0, delta=1e-4)
                self.assertAlmostEqual(laitone(cp0, mach) / pg, 1.0, delta=1e-4)
            self.assertAlmostEqual(karman_tsien(-0.5, 0.0), -0.5, places=14)
            self.assertAlmostEqual(laitone(-0.5, 0.0), -0.5, places=14)
        # The nonlinear rules predict more suction, giving a lower critical Mach.
        pg = critical_mach_number(-0.43).critical_mach
        kt = critical_mach_number(-0.43, "karman-tsien").critical_mach
        lt = critical_mach_number(-0.43, "laitone").critical_mach
        self.assertLess(kt, pg)
        self.assertLess(lt, kt)
        self.assertAlmostEqual(corrected_pressure_coefficient(-0.43, 0.6, "Karman-Tsien"),
                               karman_tsien(-0.43, 0.6), places=14)
        self.assertEqual({w.code for w in subsonic_validity_warnings(0.85)},
                         {"prandtl_glauert_high_mach", "transonic_regime"})

    def test_anderson_chapter_12_examples(self) -> None:
        data = REFERENCE["linearized"]
        e1 = data["example_12_1"]
        result = ackeret_flat_plate(e1["mach"], degrees(e1["alpha_rad_printed"]))
        self.assertAlmostEqual(result.lift_coefficient, e1["lift_coefficient"],
                               delta=e1["tolerance_abs"])
        self.assertAlmostEqual(result.drag_coefficient, e1["drag_coefficient"],
                               delta=e1["tolerance_abs"])
        e2 = data["example_12_2"]
        alpha = e2["lift_coefficient"] * sqrt(e2["mach"] ** 2 - 1.0) / 4.0
        self.assertAlmostEqual(alpha, e2["alpha_rad"], delta=e2["tolerance_abs"])
        self.assertAlmostEqual(
            ackeret_flat_plate(e2["mach"], degrees(alpha)).lift_coefficient,
            e2["lift_coefficient"], places=12)
        e3 = data["example_12_3"]
        self.assertAlmostEqual(ackeret_flat_plate(e3["mach"], degrees(e3["alpha_rad"]))
                               .drag_coefficient, e3["drag_coefficient"],
                               delta=e3["tolerance_abs"])
        exact = shock_expansion_coefficients(flat_plate(), 3.0, 5.0)
        # Anderson p. 819: linear theory within about 1.6 % of the exact result at 5 deg.
        self.assertAlmostEqual(result.lift_coefficient / exact.lift_coefficient, 1.0, delta=0.016)

    def test_busemann_is_second_order_accurate(self) -> None:
        c1, c2 = busemann_coefficients(2.0)
        self.assertAlmostEqual(c1, 2.0 / sqrt(3.0), places=14)
        self.assertAlmostEqual(c2, (2.4 * 16.0 - 12.0) / 18.0, places=14)
        errors_first, errors_second = [], []
        for theta_deg in (1.0, 2.0):
            theta = radians(theta_deg)
            exact = 2.0 / (1.4 * 4.0) * (
                calculate_oblique_shock(upstream_mach=2.0, turning_angle_deg=theta_deg)
                .pressure_ratio - 1.0)
            errors_first.append(abs(exact - ackeret_pressure_coefficient(theta, 2.0)))
            errors_second.append(abs(exact - busemann_pressure_coefficient(theta, 2.0)))
        self.assertAlmostEqual(errors_first[1] / errors_first[0], 4.0, delta=0.3)
        self.assertAlmostEqual(errors_second[1] / errors_second[0], 8.0, delta=0.8)
        self.assertLess(errors_second[0], 0.1 * errors_first[0])

    def test_linearized_section_matches_flat_plate_formulas(self) -> None:
        result = linearized_section_coefficients(flat_plate(), 2.0, 2.0)
        closed = ackeret_flat_plate(2.0, 2.0)
        alpha = radians(2.0)
        self.assertAlmostEqual(result.normal_force_coefficient, closed.lift_coefficient, places=12)
        self.assertAlmostEqual(result.lift_coefficient, closed.lift_coefficient * cos(alpha), places=12)
        thick = linearized_section_coefficients(diamond_airfoil(10.0), 1.1, 0.0)
        self.assertIn("transonic_regime", {w.code for w in thick.warnings})


class HypersonicTests(unittest.TestCase):
    def test_example_14_1_newtonian(self) -> None:
        data = REFERENCE["hypersonic"]["example_14_1b"]
        result = newtonian_flat_plate(data["alpha_deg"])
        self.assertAlmostEqual(result.lower_pressure_coefficient,
                               data["lower_pressure_coefficient"], delta=0.0005)
        self.assertAlmostEqual(result.lift_coefficient, data["lift_coefficient"], delta=0.00005)
        self.assertAlmostEqual(result.drag_coefficient, data["drag_coefficient"], delta=0.000005)
        self.assertAlmostEqual(result.lift_to_drag_ratio, data["lift_to_drag_ratio"],
                               delta=0.005)
        self.assertEqual(newtonian_pressure_coefficient(-5.0), 0.0)

    def test_maximum_lift_and_lift_to_drag(self) -> None:
        data = REFERENCE["hypersonic"]["maximum_lift"]
        alpha, cl = newtonian_maximum_lift()
        self.assertAlmostEqual(alpha, data["alpha_deg"], delta=0.05)
        self.assertAlmostEqual(cl, data["lift_coefficient"], delta=0.005)
        alpha_ld, ld = newtonian_maximum_lift_to_drag(0.008)
        a = radians(alpha_ld)
        # eq. (14.25): total drag 3 c_d,0 at (L/D)max, eq. (14.27): wave drag = 2 c_d,0.
        self.assertAlmostEqual(2.0 * a**3, 2.0 * 0.008, places=12)
        self.assertAlmostEqual(ld, 2.0 * a * a / (2.0 * a**3 + 0.008), places=12)

    def test_stagnation_pressure_coefficient(self) -> None:
        data = REFERENCE["hypersonic"]["stagnation_cp"]
        self.assertAlmostEqual(stagnation_pressure_coefficient(1.0), data["mach_1"], delta=0.005)
        self.assertAlmostEqual(stagnation_pressure_coefficient(1.0e4), data["mach_infinity"],
                               delta=0.025)
        self.assertAlmostEqual(stagnation_pressure_coefficient(1.0e-3), 1.0, delta=1e-6)
        self.assertAlmostEqual(modified_newtonian_pressure_coefficient(90.0, 6.0),
                               stagnation_pressure_coefficient(6.0), places=14)

    def test_strong_shock_limits_and_mach_independence(self) -> None:
        for gamma in (1.4, 1.2):
            near = hypersonic_shock_limits(1.0e4, 20.0, gamma)
            self.assertAlmostEqual(near.pressure_ratio / near.pressure_ratio_limit, 1.0, delta=1e-6)
            self.assertAlmostEqual(near.density_ratio / near.density_ratio_limit, 1.0, delta=1e-6)
            self.assertAlmostEqual(near.temperature_ratio / near.temperature_ratio_limit, 1.0,
                                   delta=1e-6)
            self.assertAlmostEqual(near.pressure_coefficient, near.pressure_coefficient_limit,
                                   delta=1e-7)
        self.assertAlmostEqual(limiting_wave_angle_ratio(1.4),
                               REFERENCE["hypersonic"]["limiting_wave_angle_ratio"]["ratio"],
                               places=14)
        beta_limit, cp_limit = mach_independent_wedge(15.0)
        previous = None
        for mach in (10.0, 20.0, 50.0, 200.0):
            shock = calculate_oblique_shock(upstream_mach=mach, turning_angle_deg=15.0)
            cp = 2.0 / (1.4 * mach * mach) * (shock.pressure_ratio - 1.0)
            gap = abs(cp - cp_limit)
            if previous is not None:
                self.assertLess(gap, previous)
            previous = gap
        self.assertLess(previous, 2e-4)
        self.assertAlmostEqual(shock.wave_angle_deg, beta_limit, delta=0.02)
        small_beta, _ = mach_independent_wedge(0.5)
        self.assertAlmostEqual(small_beta / 0.5, limiting_wave_angle_ratio(1.4), delta=2e-3)

    def test_small_disturbance_similarity(self) -> None:
        # Large K: Cp/delta^2 -> gamma + 1 and tau/delta -> (gamma + 1)/2 (eq. (14.35)).
        big = small_disturbance_wedge(1.0e4, 5.0, unified=False)
        self.assertAlmostEqual(big.pressure_coefficient / big.thickness_parameter**2, 2.4,
                               delta=1e-4)
        self.assertAlmostEqual(big.shock_slope / big.thickness_parameter, 1.2, delta=1e-4)
        # Small beta*delta: unified result -> Ackeret 2 delta/beta.
        thin = small_disturbance_wedge(2.0, 0.01)
        self.assertAlmostEqual(thin.pressure_coefficient /
                               ackeret_pressure_coefficient(thin.thickness_parameter, 2.0),
                               1.0, delta=1e-3)
        # Similarity: at fixed K = M delta the exact Cp/delta^2 approaches the
        # small-disturbance value as delta -> 0.
        k = 1.0
        gaps = []
        for wedge_deg in (8.0, 4.0, 2.0):
            delta = tan(radians(wedge_deg))
            mach = k / delta
            exact = calculate_oblique_shock(upstream_mach=mach, turning_angle_deg=wedge_deg)
            cp = 2.0 / (1.4 * mach * mach) * (exact.pressure_ratio - 1.0)
            approx = small_disturbance_wedge(mach, wedge_deg, unified=False)
            self.assertAlmostEqual(approx.similarity_parameter, k, places=12)
            gaps.append(abs(cp / delta**2 - approx.pressure_coefficient / delta**2))
        self.assertLess(gaps[2], gaps[1])
        self.assertLess(gaps[1], gaps[0])
        self.assertAlmostEqual(hypersonic_similarity_parameter(3.0, 0.1, unified=True),
                               0.1 * sqrt(8.0), places=14)


if __name__ == "__main__":
    unittest.main()
