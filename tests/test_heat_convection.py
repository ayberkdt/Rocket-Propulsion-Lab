"""Convective heat-transfer correlations, checked against the sources they come from.

Reference values are read from ``tests/reference_data/heat_convection.json``,
which records the source and locator of every number:

* Lienhard & Lienhard, A Heat Transfer Textbook, 6th ed. (2024): worked
  Examples 6.5, 7.1-7.8, 8.1, 8.6, 9.2, 9.5, 9.6, Tables 7.2 and 7.4;
* Taylor, NASA TN D-4332 (1968): Eq. (19)/(20), Table I, Fig. 2 (read-off);
* Aerojet-General, NASA CR-678 (1967): the Hess-Kunz ``C_L`` table, Ito's
  curvature ratio;
* Hatch & Papell, NASA TN D-130 (1959): the worked film-cooling example;
* NASA SP-8087 (1972): the 80 % burnout-margin rule.

Cross-checks between correlations in their overlap regions carry measured
agreement bands (regression guards, labelled as such in the JSON). Limits:
fully developed asymptotes, unit velocity ratio, zero temperature-ratio
exponent. Real-fluid wrappers are checked for consistency with the package's
Helmholtz equations of state.
"""

from __future__ import annotations

import json
import unittest
from math import exp, log
from pathlib import Path

from rocket_propulsion.core.errors import DomainError, InputError
from rocket_propulsion.heat_transfer import convection as cv
from rocket_propulsion.thermodynamics.real_fluids import reference_fluid

DATA = json.loads((Path(__file__).parent / "reference_data" / "heat_convection.json")
                  .read_text(encoding="utf-8"))
EX = DATA["lienhard_examples"]

# Exact English-to-SI factors used to restate the Hatch-Papell example.
FT = 0.3048
LB = 0.45359237
BTU = 1055.05585262
RANKINE = 5.0 / 9.0


def rel(a: float, b: float) -> float:
    return abs(a / b - 1.0)


class LienhardInternalFlowExamples(unittest.TestCase):
    def test_example_7_1_fully_developed_uniform_flux(self) -> None:
        e = EX["example_7_1"]
        nu = cv.laminar_fully_developed_nusselt("circular", "uniform_heat_flux").value
        self.assertEqual(nu, e["nusselt_uniform_flux_fully_developed"])
        h = cv.heat_transfer_coefficient(nu, e["conductivity_w_m_k"], e["diameter_m"])
        self.assertLess(rel(h, e["heat_transfer_coefficient_w_m2_k"]), e["printed_rel_tol"])

    def test_example_7_2_uniform_flux_entry(self) -> None:
        e = EX["example_7_2"]
        nu = cv.laminar_entry_nusselt(e["graetz"], "uniform_heat_flux").value
        self.assertLess(rel(nu, e["nusselt_local_uniform_flux"]), e["printed_rel_tol"])

    def test_example_7_3_gnielinski_with_property_corrections(self) -> None:
        e = EX["example_7_3"]
        f = cv.filonenko_friction_factor(e["reynolds"]).value
        self.assertLess(rel(f, e["darcy_friction_factor"]), e["printed_rel_tol"])
        nu = cv.gnielinski_nusselt(e["reynolds"], e["prandtl"]).value
        nu *= cv.property_ratio_nusselt_factor("liquid", e["bulk_to_wall_viscosity_ratio"],
                                               heating=True).value
        self.assertLess(rel(nu, e["nusselt"]), e["printed_rel_tol"])
        h = cv.heat_transfer_coefficient(nu, e["conductivity_w_m_k"], e["diameter_m"])
        self.assertLess(rel(h, e["heat_transfer_coefficient_w_m2_k"]), e["printed_rel_tol"])
        k = cv.property_ratio_friction_factor("liquid", e["friction_ratio_for_eq_7_46"],
                                              heating=True).value
        self.assertLess(rel(e["darcy_friction_factor"] * k, e["corrected_friction_factor"]),
                        e["printed_rel_tol"])

    def test_example_7_4_fully_rough_pipe(self) -> None:
        e, smooth = EX["example_7_4"], EX["example_7_3"]
        rel_rough = e["roughness_m"] / e["diameter_m"]
        f = cv.haaland_friction_factor(e["reynolds"], rel_rough).value
        self.assertLess(rel(f, e["darcy_friction_factor"]), e["printed_rel_tol"])
        re_eps = cv.roughness_reynolds_number(e["reynolds"], rel_rough, f)
        self.assertLess(rel(re_eps, e["roughness_reynolds"]), e["printed_rel_tol"])
        self.assertEqual(cv.roughness_regime(re_eps), "fully_rough")
        nu = cv.rough_pipe_nusselt(e["reynolds"], e["prandtl"], rel_rough).value
        self.assertLess(rel(nu, e["nusselt"]), e["printed_rel_tol"])
        h = cv.heat_transfer_coefficient(nu, smooth["conductivity_w_m_k"], e["diameter_m"])
        self.assertLess(rel(h, e["heat_transfer_coefficient_w_m2_k"]), e["h_printed_rel_tol"])
        h_smooth = smooth["heat_transfer_coefficient_w_m2_k"]
        self.assertAlmostEqual(h / h_smooth, e["h_ratio_rough_to_smooth"],
                               delta=e["ratio_printed_abs_tol"])

    def test_example_7_5_liquid_metal(self) -> None:
        e = EX["example_7_5"]
        nu = cv.liquid_metal_nusselt(e["peclet"], "uniform_heat_flux").value
        self.assertLess(rel(nu, e["nusselt"]), e["printed_rel_tol"])

    def test_example_7_6_isothermal_average_entry(self) -> None:
        e = EX["example_7_6"]
        nu = cv.laminar_entry_nusselt(e["graetz"], "uniform_temperature", average=True).value
        self.assertLess(rel(nu, e["average_nusselt_isothermal"]), e["printed_rel_tol"])

    def test_example_7_7_square_duct(self) -> None:
        e = EX["example_7_7"]
        dh = cv.rectangular_duct_hydraulic_diameter(e["side_m"], e["side_m"])
        self.assertAlmostEqual(dh, e["side_m"], places=15)
        f = cv.filonenko_friction_factor(e["reynolds"]).value
        self.assertLess(rel(f, e["darcy_friction_factor"]), e["printed_rel_tol"])
        nu = cv.gnielinski_nusselt(e["reynolds"], e["prandtl"]).value
        self.assertLess(rel(nu, e["nusselt"]), e["printed_rel_tol"])


class LienhardTables(unittest.TestCase):
    def test_table_7_4_laminar_values(self) -> None:
        t = DATA["lienhard_tables"]["table_7_4"]
        for wall in ("uniform_temperature", "uniform_heat_flux"):
            for section, value in t[wall].items():
                self.assertEqual(cv.laminar_fully_developed_nusselt(section, wall).value, value)
        self.assertEqual(cv.laminar_fully_developed_nusselt(
            "parallel_plates_one_side_insulated", "uniform_heat_flux").value,
            t["one_side_insulated_uniform_flux"])
        with self.assertRaises(InputError):
            cv.laminar_fully_developed_nusselt("parallel_plates_one_side_insulated",
                                               "uniform_temperature")

    def test_entry_fits_reach_fully_developed_limits(self) -> None:
        t = DATA["lienhard_tables"]["fig_7_4_limits"]
        gz = t["large_x_graetz"]
        flux = cv.laminar_entry_nusselt(gz, "uniform_heat_flux").value
        temp = cv.laminar_entry_nusselt(gz, "uniform_temperature").value
        self.assertLess(rel(flux, t["uniform_flux_limit"]), t["rel_tol"])
        self.assertLess(rel(temp, t["uniform_temperature_limit"]), t["rel_tol"])
        with self.assertRaises(InputError):
            cv.laminar_entry_nusselt(gz, "uniform_heat_flux", average=True)

    def test_turbulent_entry_ratio_table_7_2(self) -> None:
        t = DATA["lienhard_tables"]["table_7_2"]
        for inlet, pair in t.items():
            if inlet == "source":
                continue
            c, n = pair
            result = cv.turbulent_entry_nusselt_ratio(10.0, inlet)
            self.assertAlmostEqual(result.value, 1.0 + c / 10.0 ** n, places=14)
            self.assertIn("gas_prandtl_only", [w.code for w in result.warnings])
        with self.assertRaises(DomainError):
            cv.turbulent_entry_nusselt_ratio(3.0, "long_straight_pipe")


class LienhardExternalAndPhaseChange(unittest.TestCase):
    def test_example_6_5_laminar_plate(self) -> None:
        e = EX["example_6_5"]
        nu = cv.flat_plate_laminar_nusselt(e["reynolds"], e["prandtl"], average=True).value
        self.assertLess(rel(nu, e["average_nusselt"]), e["printed_rel_tol"])
        local = cv.flat_plate_laminar_nusselt(e["reynolds"], e["prandtl"]).value
        self.assertAlmostEqual(nu / local, 2.0, places=12)

    def test_example_7_8_cylinder(self) -> None:
        e = EX["example_7_8"]
        nu = cv.churchill_bernstein_cylinder_nusselt(e["reynolds"], e["prandtl"], "low_re").value
        self.assertLess(rel(nu, e["nusselt"]), e["printed_rel_tol"])
        full = cv.churchill_bernstein_cylinder_nusselt(e["reynolds"], e["prandtl"]).value
        self.assertGreater(full, nu)  # Eq. (7.65) bracket exceeds 1 for any Re > 0

    def test_example_8_1_churchill_chu(self) -> None:
        e = EX["example_8_1"]
        ra = cv.rayleigh_number(e["expansion_coefficient_1_k"], e["temperature_difference_k"],
                                e["length_m"], e["kinematic_viscosity_m2_s"],
                                e["thermal_diffusivity_m2_s"], gravity_m_s2=e["gravity_m_s2"])
        self.assertLess(rel(ra, e["rayleigh"]), e["rayleigh_printed_rel_tol"])
        nu = cv.churchill_chu_vertical_plate_nusselt(e["rayleigh"], e["prandtl"]).value
        self.assertLess(rel(nu, e["average_nusselt"]), e["printed_rel_tol"])
        with self.assertRaises(DomainError):
            cv.churchill_chu_vertical_plate_nusselt(1e10, 0.7)
        self.assertGreater(cv.churchill_chu_vertical_plate_nusselt(1e10, 0.7, "all_range").value,
                           0.0)

    def test_example_8_6_film_condensation(self) -> None:
        e = EX["example_8_6"]
        hfg = cv.corrected_latent_heat(e["latent_heat_j_kg"], e["liquid_specific_heat_j_kg_k"],
                                       e["saturation_minus_wall_k"], e["liquid_prandtl"]).value
        self.assertLess(rel(hfg, e["corrected_latent_heat_j_kg"]), e["printed_rel_tol"])
        nu = cv.film_condensation_vertical_plate_nusselt(
            length_m=e["length_m"], liquid_density_kg_m3=e["liquid_density_kg_m3"],
            vapor_density_kg_m3=e["vapor_density_kg_m3"],
            corrected_latent_heat_j_kg=e["corrected_latent_heat_j_kg"],
            liquid_kinematic_viscosity_m2_s=e["kinematic_viscosity_m2_s"],
            liquid_conductivity_w_m_k=e["conductivity_w_m_k"],
            saturation_minus_wall_k=e["saturation_minus_wall_k"],
            gravity_m_s2=e["gravity_m_s2"]).value
        self.assertLess(rel(nu, e["average_nusselt"]), e["nusselt_rel_tol_from_rounded_delta"])

    def test_example_9_2_rohsenow(self) -> None:
        e = EX["example_9_2"]
        props = dict(specific_heat_j_kg_k=e["specific_heat_j_kg_k"],
                     latent_heat_j_kg=e["latent_heat_j_kg"], prandtl=e["prandtl"],
                     viscosity_pa_s=e["viscosity_pa_s"],
                     liquid_density_kg_m3=e["density_difference_kg_m3"] + 1.0,
                     vapor_density_kg_m3=1.0, surface_tension_n_m=e["surface_tension_n_m"],
                     surface_constant=e["surface_constant"], prandtl_exponent=1.0,
                     gravity_m_s2=e["gravity_m_s2"])
        q = cv.rohsenow_heat_flux(e["wall_superheat_k"], **props).value
        self.assertLess(rel(q, e["heat_flux_w_m2"]), e["rel_tol"])
        dt = cv.rohsenow_wall_superheat(q, **props).value
        self.assertAlmostEqual(dt, e["wall_superheat_k"], places=9)

    def test_examples_9_5_and_9_6_peak_heat_flux(self) -> None:
        for key in ("example_9_5", "example_9_6"):
            e = EX[key]
            q = cv.peak_pool_boiling_heat_flux(
                vapor_density_kg_m3=e["vapor_density_kg_m3"], latent_heat_j_kg=e["latent_heat_j_kg"],
                surface_tension_n_m=e["surface_tension_n_m"],
                liquid_density_kg_m3=e["liquid_density_kg_m3"], gravity_m_s2=e["gravity_m_s2"])
            self.assertLess(rel(q.value, e["peak_heat_flux_w_m2"]), e["printed_rel_tol"], key)
        e = EX["example_9_5"]
        zuber = cv.peak_pool_boiling_heat_flux(
            vapor_density_kg_m3=e["vapor_density_kg_m3"], latent_heat_j_kg=e["latent_heat_j_kg"],
            surface_tension_n_m=e["surface_tension_n_m"],
            liquid_density_kg_m3=e["liquid_density_kg_m3"], gravity_m_s2=e["gravity_m_s2"],
            form="zuber").value
        self.assertAlmostEqual(zuber / (e["peak_heat_flux_w_m2"]), 0.131 / 0.149, delta=1e-3)

    def test_taylor_wavelength_and_small_heater_screen(self) -> None:
        w = EX["taylor_wavelength_water"]
        lam = cv.taylor_wavelength(w["liquid_density_kg_m3"], w["vapor_density_kg_m3"],
                                   w["surface_tension_n_m"], w["gravity_m_s2"])
        self.assertLess(rel(lam, w["wavelength_m"]), w["printed_rel_tol"])
        kwargs = dict(vapor_density_kg_m3=w["vapor_density_kg_m3"], latent_heat_j_kg=2.257e6,
                      surface_tension_n_m=w["surface_tension_n_m"],
                      liquid_density_kg_m3=w["liquid_density_kg_m3"],
                      gravity_m_s2=w["gravity_m_s2"])
        with self.assertRaises(DomainError):
            cv.peak_pool_boiling_heat_flux(heater_width_m=2.0 * lam, **kwargs)
        ok = cv.peak_pool_boiling_heat_flux(heater_width_m=4.0 * lam, **kwargs)
        self.assertEqual(ok.warnings, ())


class TaylorHydrogen(unittest.TestCase):
    def test_exponent_matches_figure_2_readoff(self) -> None:
        t = DATA["taylor_tn_d4332"]
        tol = t["readoff_abs_tol"]
        for xd, c2 in t["fig2_present_report_readoff"]:
            self.assertAlmostEqual(cv.taylor_temperature_ratio_exponent(xd), c2, delta=tol)
        # The printed comparison curves of Eqs. (9) and (13) read off the same figure.
        for xd, c2 in t["fig2_ref12_13_eq13_readoff"]:
            self.assertAlmostEqual(0.29 + 0.0019 * xd, c2, delta=tol)
        for xd, c2 in t["fig2_ref9_eq9_readoff"]:
            self.assertAlmostEqual(0.29 + 0.0056 * xd, c2, delta=tol)

    def test_zero_exponent_reduces_to_rte_default(self) -> None:
        xd = 1.59 / 0.57
        a = cv.taylor_hydrogen_nusselt(2e5, 0.8, 5.0, xd).value
        b = cv.rte_coolant_nusselt(2e5, 0.8).value
        self.assertAlmostEqual(a, b, delta=1e-9 * b)

    def test_table_i_range_fails_closed(self) -> None:
        r = DATA["taylor_tn_d4332"]["table_i"]
        cv.taylor_hydrogen_nusselt(r["reynolds_bulk"][0], 0.7, r["wall_to_bulk"][1],
                                   r["x_over_d"][0])
        for args in ((1e5, 0.7, 2.0, 1.0), (1e5, 0.7, 30.0, 50.0), (5000.0, 0.7, 2.0, 50.0)):
            with self.assertRaises(DomainError) as ctx:
                cv.taylor_hydrogen_nusselt(*args)
            self.assertEqual(ctx.exception.code, "outside_correlation_range")
        warned = cv.taylor_hydrogen_nusselt(1e5, 0.7, 1.0, 50.0, extrapolate=True)
        self.assertEqual(warned.warnings[0].code, "correlation_extrapolated")

    def test_real_fluid_coefficient_is_consistent_with_eos(self) -> None:
        kw = dict(bulk_temperature_k=60.0, wall_temperature_k=400.0, pressure_pa=6.0e6,
                  mass_flux_kg_m2_s=4000.0, diameter_m=0.004, distance_from_entrance_m=0.2,
                  bulk_viscosity_pa_s=4.0e-6, bulk_conductivity_w_m_k=0.09)
        result = cv.taylor_hydrogen_coefficient(fluid="parahydrogen", **kw)
        state = reference_fluid("parahydrogen").state_tp(60.0, 6.0e6)
        re = 4000.0 * 0.004 / 4.0e-6
        pr = state.cp_j_kg_k * 4.0e-6 / 0.09
        nu = cv.taylor_hydrogen_nusselt(re, pr, 400.0 / 60.0, 50.0).value
        self.assertAlmostEqual(result.heat_transfer_coefficient_w_m2_k, nu * 0.09 / 0.004,
                               delta=1e-9 * nu)
        self.assertEqual(result.specific_heat_j_kg_k, state.cp_j_kg_k)
        with self.assertRaises(InputError):
            cv.taylor_hydrogen_coefficient(fluid="oxygen", **kw)

    def test_near_critical_region_screen(self) -> None:
        reg = DATA["taylor_tn_d4332"]["region_2"]
        t_star = cv.pseudocritical_temperature("parahydrogen", 2.0e6)
        f = reference_fluid("parahydrogen")
        cp_star = f.state_tp(t_star, 2.0e6).cp_j_kg_k
        for dt in (-0.5, 0.5):
            self.assertGreater(cp_star, f.state_tp(t_star + dt, 2.0e6).cp_j_kg_k)
        inside = cv.taylor_near_critical_region(reg["inlet_temperature_min_k"] + 5.0, 2.0e6)
        self.assertTrue(inside.in_region)
        self.assertAlmostEqual(inside.pseudocritical_temperature_k, t_star, delta=1e-4)
        self.assertEqual(inside.warnings[0].code, "taylor_near_critical_region")
        self.assertFalse(cv.taylor_near_critical_region(t_star + 5.0, 2.0e6).in_region)
        self.assertFalse(cv.taylor_near_critical_region(30.0, reg["inlet_pressure_max_pa"]
                                                        + 1.0e5).in_region)
        below = cv.taylor_near_critical_region(30.0, 1.0e6)
        self.assertFalse(below.in_region)
        self.assertIsNone(below.pseudocritical_temperature_k)


class HessKunzAndSp125(unittest.TestCase):
    def test_cl_table_points_and_range(self) -> None:
        t = DATA["cr678_cl_table"]
        for geometry in ("straight", "curved"):
            for rankine, cl in zip(t["rankine"], t[geometry], strict=True):
                if cl is None:
                    with self.assertRaises(DomainError):
                        cv.hess_kunz_cl_factor(rankine * RANKINE, geometry)
                    continue
                value = cv.hess_kunz_cl_factor(rankine * RANKINE, geometry).value
                self.assertAlmostEqual(value, cl, places=9)
        mid = cv.hess_kunz_cl_factor(72.5 * RANKINE, "straight").value
        self.assertAlmostEqual(mid, 0.5 * (1.07 + 0.93), places=9)

    def test_hess_kunz_form_and_wrapper(self) -> None:
        base = cv.hess_kunz_nusselt(1e5, 0.9, 2.0)
        self.assertAlmostEqual(base.value, 0.0208 * 1e5 ** 0.8 * 0.9 ** 0.4 * (1 + 0.02914),
                               delta=1e-9 * base.value)
        self.assertEqual(base.warnings[0].code, "range_not_stated")
        scaled = cv.hess_kunz_nusselt(1e5, 0.9, 2.0, cl_factor=0.85).value
        self.assertAlmostEqual(scaled / base.value, 0.85, places=12)
        res = cv.hess_kunz_coefficient(fluid="normal_hydrogen", bulk_temperature_k=50.0,
                                       wall_temperature_k=250.0, pressure_pa=7.0e6,
                                       mass_flux_kg_m2_s=3000.0, diameter_m=0.005,
                                       film_viscosity_pa_s=7.0e-6, film_conductivity_w_m_k=0.12,
                                       bulk_viscosity_pa_s=5.0e-6, wall_viscosity_pa_s=9.0e-6)
        f = reference_fluid("normal_hydrogen")
        b, m, w = (f.state_tp(t, 7.0e6) for t in (50.0, 150.0, 250.0))
        re_f = m.density_kg_m3 * (3000.0 / b.density_kg_m3) * 0.005 / 7.0e-6
        pr_f = m.cp_j_kg_k * 7.0e-6 / 0.12
        ratio = (9.0e-6 / w.density_kg_m3) / (5.0e-6 / b.density_kg_m3)
        nu = cv.hess_kunz_nusselt(re_f, pr_f, ratio).value
        self.assertAlmostEqual(res.heat_transfer_coefficient_w_m2_k, nu * 0.12 / 0.005,
                               delta=1e-9 * nu)
        self.assertEqual(res.reference_temperature_k, 150.0)

    def test_sp125_eq_4_25_nusselt_form(self) -> None:
        # Algebraic identity of Eq. (4-25): h d / k = 0.029 Re^0.8 Pr^(1/3) (Tco/Twc)^0.55.
        cp, mu, k, g, d = 1.4e4, 1.2e-5, 0.15, 2500.0, 0.003
        pr = cp * mu / k
        h = cv.sp125_supercritical_film_coefficient(cp, mu, pr, g, d, 100.0, 300.0).value
        re = g * d / mu
        self.assertAlmostEqual(h * d / k, 0.029 * re ** 0.8 * pr ** (1 / 3) * (1 / 3) ** 0.55,
                               delta=1e-9 * h)
        res = cv.sp125_supercritical_film_heat_transfer(
            fluid="parahydrogen", bulk_temperature_k=100.0, wall_temperature_k=300.0,
            pressure_pa=5.0e6, mass_flux_kg_m2_s=g, hydraulic_diameter_m=d, viscosity_pa_s=mu,
            conductivity_w_m_k=k)
        codes = [w.code for w in res.warnings]
        self.assertIn("bulk_reference_state_assumed", codes)
        self.assertNotIn("not_supercritical", codes)
        st = cv.sp125_sieder_tate_nusselt(1e5, 3.0, 1.5, 0.023).value
        self.assertAlmostEqual(st, 0.023 * 1e5 ** 0.8 * 3.0 ** 0.4 * 1.5 ** 0.14, delta=1e-9 * st)


class CurvatureFrictionAndScreens(unittest.TestCase):
    def test_curvature_factor_matches_ito_group(self) -> None:
        concave = cv.curvature_heat_transfer_factor(2e5, 0.001, 0.02, "concave", "gas")
        convex = cv.curvature_heat_transfer_factor(2e5, 0.001, 0.02, "convex", "gas").value
        ito = cv.ito_curvature_friction_factor(2e5, 0.001, 0.02).value
        self.assertAlmostEqual(concave.value, ito, places=14)
        self.assertAlmostEqual(concave.value * convex, 1.0, places=14)
        liquid = cv.curvature_heat_transfer_factor(2e5, 0.001, 0.02, "concave", "liquid")
        self.assertIn("curvature_factor_not_recommended", [w.code for w in liquid.warnings])

    def test_chen_and_haaland_agree(self) -> None:
        c = DATA["cross_checks"]
        worst = 0.0
        for re in c["chen_vs_haaland_grid"]["reynolds"]:
            for e in c["chen_vs_haaland_grid"]["relative_roughness"]:
                a = cv.haaland_friction_factor(re, e).value
                b = cv.chen_friction_factor(re, e).value
                worst = max(worst, rel(b, a))
        self.assertLessEqual(worst, c["measured_chen_vs_haaland_max_rel"] * 1.01)

    def test_gnielinski_forms_overlap(self) -> None:
        c = DATA["cross_checks"]
        worst = 0.0
        for re in c["gnielinski_grid"]["reynolds"]:
            for pr in c["gnielinski_grid"]["prandtl"]:
                worst = max(worst, rel(cv.gnielinski_power_law_nusselt(re, pr).value,
                                       cv.gnielinski_nusselt(re, pr).value))
        self.assertLessEqual(worst, c["measured_gnielinski_power_law_vs_7_41_max_rel"])
        # Eq. (7.41) differs from Eq. (7.35) only by (Re - 1000)/Re.
        for re in (1e4, 1e5, 1e6):
            ratio = (cv.gnielinski_nusselt(re, 5.0).value
                     / cv.smooth_pipe_analogy_nusselt(re, 5.0).value)
            self.assertAlmostEqual(ratio, (re - 1000.0) / re, places=12)

    def test_flat_plate_turbulent_forms_overlap_for_air(self) -> None:
        c = DATA["cross_checks"]
        for re in c["flat_plate_grid"]["reynolds"]:
            pr = c["flat_plate_grid"]["prandtl"]
            a = cv.flat_plate_turbulent_nusselt(re, pr).value
            b = cv.flat_plate_turbulent_power_law_nusselt(re, pr, "gas").value
            self.assertLessEqual(rel(a, b), c["measured_flat_plate_6_111_vs_6_112_air_max_rel"])

    def test_burnout_margin_screen(self) -> None:
        limit = DATA["sp8087"]["allowable_fraction"]
        ok = cv.burnout_margin_screen(0.7e7, 1.0e7)
        self.assertAlmostEqual(ok.value, 0.7)
        self.assertEqual(ok.warnings, ())
        bad = cv.burnout_margin_screen((limit + 0.05) * 1.0e7, 1.0e7)
        self.assertEqual(bad.warnings[0].severity, "critical")

    def test_fail_closed_and_extrapolation(self) -> None:
        with self.assertRaises(DomainError) as ctx:
            cv.gnielinski_nusselt(1000.0, 1.0)
        self.assertEqual(ctx.exception.code, "outside_correlation_range")
        warned = cv.gnielinski_nusselt(1000.0, 1.0, extrapolate=True)
        self.assertIn("correlation_extrapolated", [w.code for w in warned.warnings])
        with self.assertRaises(DomainError):
            cv.power_law_nusselt(5000.0, 1.0)
        with self.assertRaises(DomainError):
            cv.rough_pipe_nusselt(2e4, 1.0, 0.003)  # Re_eps < 70: not fully rough
        with self.assertRaises(DomainError):
            cv.property_ratio_nusselt_factor("gas", 3.0, heating=True)
        with self.assertRaises(DomainError):
            cv.reynolds_number(1.0, -1.0, 1.0, 1.0)
        with self.assertRaises(InputError):
            cv.power_law_nusselt(1e5, 1.0, "unknown")


class HatchPapellFilmCooling(unittest.TestCase):
    def setUp(self) -> None:
        e = DATA["hatch_papell_example"]
        self.e = e
        self.common = dict(
            heat_transfer_coefficient_w_m2_k=e["h_btu_ft2_s_r"] * BTU / (FT ** 2 * RANKINE),
            slot_length_m=e["slot_length_ft"] * FT, distance_m=e["distance_ft"] * FT,
            coolant_specific_heat_j_kg_k=e["coolant_cp_btu_lb_r"] * BTU / (LB * RANKINE),
            slot_height_m=e["slot_height_ft"] * FT, gas_velocity_m_s=e["gas_velocity_ft_s"] * FT,
            coolant_velocity_m_s=e["coolant_velocity_ft_s"] * FT,
            coolant_thermal_diffusivity_m2_s=e["coolant_diffusivity_ft2_s"] * FT ** 2,
            gas_temperature_k=e["gas_temperature_r"] * RANKINE,
            coolant_temperature_k=e["coolant_temperature_r"] * RANKINE)

    def test_velocity_function_and_limits(self) -> None:
        e = self.e
        f = cv.hatch_papell_velocity_ratio_function(e["gas_velocity_ft_s"] / e["coolant_velocity_ft_s"])
        self.assertAlmostEqual(f, e["velocity_function"], delta=0.005)
        self.assertEqual(cv.hatch_papell_velocity_ratio_function(1.0), 1.0)
        self.assertAlmostEqual(cv.hatch_papell_velocity_ratio_function(1.0 - 1e-9), 1.0, places=6)

    def test_worked_example_effectiveness(self) -> None:
        e = self.e
        eta_def = ((e["adiabatic_wall_temperature_r"] - e["wall_temperature_r"])
                   / (e["adiabatic_wall_temperature_r"] - e["coolant_temperature_r"]))
        self.assertAlmostEqual(eta_def, e["effectiveness"], delta=e["eta_definition_abs_tol"])
        with self.assertRaises(DomainError):  # 2000 ft/s is above the 1040 ft/s data range
            cv.hatch_papell_effectiveness(coolant_mass_flow_kg_s=e["first_iteration_flow_lb_s"]
                                          * LB, **self.common)
        r = cv.hatch_papell_effectiveness(coolant_mass_flow_kg_s=e["first_iteration_flow_lb_s"]
                                          * LB, extrapolate=True, **self.common)
        self.assertLess(rel(r.value, e["effectiveness"]), e["printed_rel_tol"])
        self.assertLess(rel(log(r.value), e["ln_effectiveness"]), e["printed_rel_tol"])
        self.assertEqual([w.field for w in r.warnings], ["gas_velocity_m_s"])
        wall = cv.film_cooled_wall_temperature(r.value, e["adiabatic_wall_temperature_r"] * RANKINE,
                                               e["coolant_temperature_r"] * RANKINE) / RANKINE
        self.assertAlmostEqual(wall, e["wall_temperature_r"], delta=2.0)

    def test_worked_example_required_flow_two_iterations(self) -> None:
        e = self.e
        eta = exp(e["ln_effectiveness"])
        w1 = cv.hatch_papell_required_coolant_flow(effectiveness=eta, extrapolate=True,
                                                   **self.common).value / LB
        self.assertLess(rel(w1, e["first_iteration_flow_lb_s"]), e["printed_rel_tol"])
        second = dict(self.common, slot_height_m=e["second_iteration_slot_ft"] * FT)
        w2 = cv.hatch_papell_required_coolant_flow(effectiveness=eta, extrapolate=True,
                                                   **second).value / LB
        self.assertLess(rel(w2, e["second_iteration_flow_lb_s"]), e["printed_rel_tol"])

    def test_no_cooling_length_gives_unit_effectiveness(self) -> None:
        r = cv.hatch_papell_effectiveness(coolant_mass_flow_kg_s=1.0e3, extrapolate=True,
                                          **self.common)
        self.assertEqual(r.value, 1.0)
        self.assertIn("Eq. (12a)", r.source)


if __name__ == "__main__":
    unittest.main()
