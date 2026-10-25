"""Transpiration cooling, checked against the reports it comes from.

Reference values are read from ``tests/reference_data/transpiration.json``:

* Mickley et al., NACA TN 3208 (1954), Table VII: the authors' own
  film-theory coefficients ``h_T = h*_T phi/(e^phi - 1)`` for six
  compartments, and their measured coefficients for two downstream ones;
* Kelley & L'Ecuyer, NASA CR-76367 (1966), Appendix B: Rannie's worked
  example and Table B.1, four independent methods on one flat-plate case.

Analytic checks: the film-theory limits of TN 3208 footnote 2, the
closed-form wall temperature against Eckert & Livingood Eq. (30) with
``r = 1``, Eq. (31) reducing to Eq. (30) without radiation, the energy
balance of every solved state, and the Pratt & Whitney porous-wall solution
against its own boundary conditions and differential equations.
"""

from __future__ import annotations

import json
import unittest
from math import exp, log
from pathlib import Path
from typing import ClassVar

from rocket_propulsion.core.errors import DomainError
from rocket_propulsion.heat_transfer import transpiration as tp

DATA = json.loads((Path(__file__).parent / "reference_data" / "transpiration.json")
                  .read_text(encoding="utf-8"))
RANKINE = 5.0 / 9.0


def rel(a: float, b: float) -> float:
    return abs(a / b - 1.0)


class MickleyFilmTheory(unittest.TestCase):
    def test_reproduces_the_printed_film_theory_coefficients(self) -> None:
        table = DATA["mickley_table_vii"]
        for row in table["rows"]:
            with self.subTest(run=row["run"], compartment=row["compartment"]):
                h = row["h_star"] * tp.film_theory_correction(row["phi"]).value
                self.assertLess(rel(h, row["h_film"]), table["printed_rel_tol"])

    def test_limits_of_footnote_2(self) -> None:
        self.assertEqual(tp.film_theory_correction(0.0).value, 1.0)
        self.assertAlmostEqual(tp.film_theory_correction(1e-9).value, 1.0 - 5e-10, places=15)
        self.assertLess(tp.film_theory_correction(50.0).value, 1e-19)
        self.assertLess(rel(tp.film_theory_correction(-40.0).value, 40.0), 1e-15)

    def test_series_branch_is_continuous_with_the_closed_form(self) -> None:
        for phi in (-1.0001e-3, -0.9999e-3, 0.9999e-3, 1.0001e-3):
            with self.subTest(phi=phi):
                self.assertLess(rel(tp.film_theory_correction(phi).value, phi / (exp(phi) - 1.0)), 2e-13)

    def test_resistance_factor_identity(self) -> None:
        # Eq. (28) R = phi/theta with Eq. (30) phi = ln(1 + R).
        for r in (0.1, 1.0, 4.0, -0.5):
            phi = tp.film_theory_rate_factor(r)
            self.assertLess(rel(phi / tp.film_theory_correction(phi).value, r), 1e-12)
        self.assertLess(rel(tp.film_theory_rate_factor(1.0), log(2.0)), 1e-15)
        with self.assertRaises(DomainError):
            tp.film_theory_rate_factor(-1.0)

    def test_measured_downstream_coefficients_exceed_film_theory(self) -> None:
        """Characterisation, not accuracy: measured h is 25-28 % above film theory.

        TN 3208 pp. 53-54 report the porous plate's unblown turbulent friction
        15-24 % above smooth-plate values; the two downstream rows show the
        same direction for heat transfer.
        """
        rows = {f"{r['run']} {r['compartment']}": r for r in DATA["mickley_table_vii"]["rows"]}
        for key in DATA["mickley_table_vii"]["downstream_rows_for_measurement"]:
            row = rows[key]
            ratio = row["h_measured"] / (row["h_star"] * tp.film_theory_correction(row["phi"]).value)
            with self.subTest(row=key):
                self.assertGreater(ratio, 1.20)
                self.assertLess(ratio, 1.32)


class KelleyRannieExample(unittest.TestCase):
    def test_printed_worked_example(self) -> None:
        c = DATA["kelley_appendix_b"]["conditions"]
        printed = DATA["kelley_appendix_b"]["rannie_printed"]
        value = tp.rannie_pipe_temperature_ratio(printed["blowing_rate"], c["reynolds"], c["prandtl"]).value
        # The source rounds to 0.404; Eq. (2.17) evaluates to 0.4069 (0.7 %).
        self.assertLess(rel(value, printed["temperature_ratio"]), 0.01)
        target = (c["freestream_temperature_r"] - c["wall_temperature_r"]) / (
            c["freestream_temperature_r"] - c["coolant_reservoir_temperature_r"])
        # 1000/2470 = 0.4049, printed truncated as 0.404.
        self.assertLess(rel(target, printed["temperature_ratio"]), 0.003)
        f = tp.rannie_pipe_blowing_rate(target, c["reynolds"], c["prandtl"])
        self.assertLess(rel(f, printed["blowing_rate"]), 0.015)

    def test_rannie_6000r_column_of_table_b1(self) -> None:
        c = DATA["kelley_appendix_b"]["conditions"]
        target = (6000.0 - c["wall_temperature_r"]) / (6000.0 - c["coolant_reservoir_temperature_r"])
        f = tp.rannie_pipe_blowing_rate(target, c["reynolds"], c["prandtl"])
        self.assertLess(rel(f, DATA["kelley_appendix_b"]["table_b1_blowing_rate_6000r"]["Rannie"]), 0.03)

    def test_film_theory_lies_inside_the_four_published_methods(self) -> None:
        """Same case as Table B.1 (air into air, T_r = 3041 R, St0 = 1.5e-3).

        Film theory with the coolant energy balance gives F = St0 ln(1 + R),
        R = (T_r - T_w)/(T_w - T_0); it must fall inside the 0.63e-3 to
        1.07e-3 spread of Rannie, Rubesin, Spalding et al. and Bartle-Leadon,
        and lies within 5 % of Rubesin's boundary-layer solution, the one that
        uses the same St0.
        """
        c = DATA["kelley_appendix_b"]["conditions"]
        state = tp.transpiration_coolant_flux(
            unblown_heat_transfer_coefficient_w_m2_k=c["unblown_stanton"],
            recovery_temperature_k=c["recovery_temperature_r"] * RANKINE,
            wall_temperature_k=c["wall_temperature_r"] * RANKINE,
            coolant_temperature_k=c["coolant_reservoir_temperature_r"] * RANKINE,
            coolant_specific_heat_j_kg_k=1.0)
        f = state.coolant_mass_flux_kg_m2_s  # rho u c_p = 1, so the flux is F
        published = DATA["kelley_appendix_b"]["table_b1_blowing_rate_3000r"]
        self.assertGreater(f, min(published.values()))
        self.assertLess(f, max(published.values()))
        self.assertLess(rel(f, published["Rubesin"]), 0.05)


class WallEnergyBalance(unittest.TestCase):
    KW: ClassVar[dict[str, float]] = {
        "unblown_heat_transfer_coefficient_w_m2_k": 4000.0, "recovery_temperature_k": 3300.0,
        "coolant_temperature_k": 100.0, "coolant_specific_heat_j_kg_k": 14300.0}

    def test_coolant_flux_closes_the_balance(self) -> None:
        for eta in (1.0, 0.6):
            for q_r in (0.0, 2.0e6):
                with self.subTest(eta=eta, q_r=q_r):
                    s = tp.transpiration_coolant_flux(wall_temperature_k=900.0, wall_efficiency=eta,
                                                      radiative_heat_flux_w_m2=q_r, **self.KW)
                    load = s.convective_heat_flux_w_m2 + q_r
                    self.assertLess(rel(s.absorbed_heat_flux_w_m2, load), 1e-10)
                    self.assertLess(s.convective_heat_flux_w_m2, s.unblown_heat_flux_w_m2)

    def test_wall_temperature_inverts_coolant_flux(self) -> None:
        for eta, q_r in ((1.0, 0.0), (0.7, 0.0), (1.0, 3.0e6), (0.5, 1.0e6)):
            with self.subTest(eta=eta, q_r=q_r):
                s = tp.transpiration_coolant_flux(wall_temperature_k=1000.0, wall_efficiency=eta,
                                                  radiative_heat_flux_w_m2=q_r, **self.KW)
                back = tp.transpiration_wall_temperature(
                    coolant_mass_flux_kg_m2_s=s.coolant_mass_flux_kg_m2_s, wall_efficiency=eta,
                    radiative_heat_flux_w_m2=q_r, **self.KW)
                self.assertLess(abs(back.wall_temperature_k - 1000.0), 1e-8)

    def test_unit_efficiency_matches_eckert_livingood_film_limit(self) -> None:
        # Eq. (30) with r = 1 and Re^0.1/2.11 -> 1 is 1/(1 + (e^phi - 1)) = e^-phi.
        s = tp.transpiration_wall_temperature(coolant_mass_flux_kg_m2_s=0.3, **self.KW)
        ratio = (s.wall_temperature_k - 100.0) / (3300.0 - 100.0)
        self.assertLess(rel(ratio, exp(-s.rate_factor)), 1e-12)

    def test_zero_flow_gives_recovery_temperature(self) -> None:
        s = tp.transpiration_wall_temperature(coolant_mass_flux_kg_m2_s=0.0, **self.KW)
        self.assertAlmostEqual(s.wall_temperature_k, 3300.0, places=9)

    def test_rejects_impossible_temperatures(self) -> None:
        with self.assertRaises(DomainError):
            tp.transpiration_coolant_flux(wall_temperature_k=3400.0, **self.KW)
        with self.assertRaises(DomainError):
            tp.transpiration_coolant_flux(wall_temperature_k=50.0, **self.KW)


class EckertLivingood(unittest.TestCase):
    def test_eq31_without_radiation_is_eq30(self) -> None:
        for ratio in (0.0, 0.002, 0.006, 0.012):
            re, pr = 1.0e7, 0.7
            value = tp.eckert_livingood_temperature_ratio(ratio, re, pr).value
            r_phi = (2.11 / 0.037) * re ** 0.1 * pr ** (2.0 / 3.0) * ratio
            eq30 = 1.0 / (1.0 + re ** 0.1 / 2.11 * (exp(r_phi) - 1.0))
            with self.subTest(ratio=ratio):
                self.assertLess(abs(value - eq30), 1e-14)

    def test_rannie_friedman_is_more_conservative_than_film_theory(self) -> None:
        """Report 1182 p. 11: film-theory curves lie below Rannie-Friedman."""
        for ratio in (0.002, 0.005, 0.01):
            rf = tp.eckert_livingood_temperature_ratio(ratio, 1e7, 0.7).value
            film = tp.eckert_livingood_temperature_ratio(ratio, 1e7, 0.7, model="film").value
            self.assertGreater(rf, film)

    def test_radiation_raises_the_wall_temperature(self) -> None:
        base = tp.eckert_livingood_temperature_ratio(0.004, 1e7, 0.7).value
        hot = tp.eckert_livingood_temperature_ratio(0.004, 1e7, 0.7,
                                                    radiation_to_convection_ratio=1.0).value
        self.assertGreater(hot, base)

    def test_outside_computed_range(self) -> None:
        with self.assertRaises(DomainError):
            tp.eckert_livingood_temperature_ratio(0.02, 1e7, 0.7)
        warned = tp.eckert_livingood_temperature_ratio(0.02, 1e7, 0.7, extrapolate=True)
        self.assertEqual(warned.warnings[0].code, "correlation_extrapolated")


class PorousWall(unittest.TestCase):
    KW: ClassVar[dict[str, float]] = {
        "thickness_m": 0.003, "volumetric_heat_transfer_w_m3_k": 5.0e7,
        "effective_conductivity_w_m_k": 20.0, "coolant_mass_flux_kg_m2_s": 1.5,
        "coolant_specific_heat_j_kg_k": 14300.0, "wall_temperature_k": 900.0,
        "coolant_inlet_temperature_k": 100.0}

    def test_boundary_conditions_and_equations(self) -> None:
        s = tp.porous_wall_solution(**self.KW)
        d = s.thickness_m
        self.assertLess(abs(s.solid_temperature_k(d) - 900.0), 1e-9)
        self.assertLess(abs(s.coolant_temperature_k(0.0) - 100.0), 1e-9)
        self.assertLess(abs(s.coolant_temperature_k(d) - s.coolant_exit_temperature_k), 1e-9)
        _, c2, c3 = s.coefficients_k
        p1, p2 = s.roots_per_m
        self.assertLess(abs(c2 * p1 + c3 * p2), 1e-9 * abs(c2 * p1))  # T_m'(0) = 0
        for x in (0.0, 0.4 * d, d):
            tm, tc = s.solid_temperature_k(x), s.coolant_temperature_k(x)
            tm2 = c2 * p1 * p1 * exp(p1 * x) + c3 * p2 * p2 * exp(p2 * x)
            tc1 = (s.b_per_m / s.a_per_m2) * (c2 * p1 * p1 * exp(p1 * x) + c3 * p2 * p2 * exp(p2 * x))
            with self.subTest(x=x):
                self.assertLess(abs(tm2 - s.a_per_m2 * (tm - tc)), 1e-7 * abs(tm2) + 1e-6)
                self.assertLess(abs(tc1 - s.b_per_m * (tm - tc)), 1e-7 * abs(tc1) + 1e-6)

    def test_efficiency_rises_with_internal_exchange(self) -> None:
        low = tp.porous_wall_solution(**{**self.KW, "volumetric_heat_transfer_w_m3_k": 1e6}).efficiency
        high = tp.porous_wall_solution(**{**self.KW, "volumetric_heat_transfer_w_m3_k": 1e9}).efficiency
        self.assertLess(low, high)
        self.assertLess(high, 1.0)
        self.assertGreater(high, 0.95)


if __name__ == "__main__":
    unittest.main()
