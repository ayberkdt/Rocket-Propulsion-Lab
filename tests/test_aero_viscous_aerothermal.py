"""Source-traceable checks of the viscous-flow, aerothermal and entry modules.

Reference values and their page locators are in
``tests/reference_data/aero_viscous.json``.
"""

from __future__ import annotations

import json
import math
import unittest
from pathlib import Path

from rocket_propulsion.aerodynamics.aerothermal import (
    SUTTON_GRAVES_TABLE_II,
    newtonian_stagnation_velocity_gradient,
    stagnation_point_heating_boundary_layer,
    sutton_graves_coefficient,
    sutton_graves_heating,
    tauber_laminar_cone_heating,
    tauber_laminar_flat_plate_heating,
    tauber_leading_edge_heating,
    tauber_stagnation_heating,
    tauber_swept_cylinder_heating,
    tauber_turbulent_cone_heating,
    tauber_turbulent_flat_plate_heating,
)
from rocket_propulsion.aerodynamics.boundary_layer import (
    SutherlandViscosity,
    adiabatic_wall_temperature,
    compressible_flat_plate_reference_temperature,
    flat_plate_regime,
    laminar_flat_plate_incompressible,
    reference_temperature_ratio,
    reynolds_analogy_stanton,
    solve_blasius,
    turbulent_flat_plate_incompressible,
    van_driest_ii,
)
from rocket_propulsion.aerodynamics.couette import compressible_couette, incompressible_couette
from rocket_propulsion.aerodynamics.entry import (
    USSA1976_MAX_GEOMETRIC_ALTITUDE_M,
    ExponentialAtmosphere,
    allen_eggers_entry,
    integrate_planar_entry,
    total_heat_input_reynolds_analogy,
    ussa1976_density,
    ussa1976_geopotential_altitude,
)
from rocket_propulsion.core.errors import DomainError, InputError

REFERENCE = json.loads((Path(__file__).parent / "reference_data" / "aero_viscous.json").read_text(encoding="utf-8"))


def _relative(actual: float, expected: float) -> float:
    return abs(actual / expected - 1.0)


class BlasiusTests(unittest.TestCase):
    def test_blasius_constants_match_anderson(self) -> None:
        case = REFERENCE["blasius"]
        solution = solve_blasius()
        self.assertLess(solution.edge_residual, 1.0e-10)
        for field, expected in case["expected"].items():
            with self.subTest(field=field):
                self.assertLessEqual(_relative(getattr(solution, field), expected), case["relative_tolerance"][field])

    def test_grid_refinement_changes_wall_curvature_below_1e9(self) -> None:
        coarse = solve_blasius(steps=2000).wall_curvature
        fine = solve_blasius(steps=4000).wall_curvature
        self.assertLess(abs(coarse - fine), 1.0e-9)

    def test_momentum_thickness_equals_twice_wall_curvature(self) -> None:
        # Momentum integral for a flat plate: theta sqrt(Re)/x = cf sqrt(Re)/2 = f''(0)*2/2*... = 2 f''(0).
        solution = solve_blasius()
        self.assertAlmostEqual(solution.momentum_coefficient, 2.0 * solution.wall_curvature, delta=2.0e-5)

    def test_rejects_bad_grid(self) -> None:
        with self.assertRaises(DomainError):
            solve_blasius(eta_max=3.0)
        with self.assertRaises(DomainError):
            solve_blasius(steps=10)


class AndersonFlatPlateExampleTests(unittest.TestCase):
    def _cases(self):
        return {case["id"]: case for case in REFERENCE["anderson_flat_plate_examples"]}

    def test_example_18_1a(self) -> None:
        case = self._cases()["example_18_1a_laminar_incompressible"]
        i = case["inputs"]
        reynolds = i["density"] * i["velocity"] * i["chord"] / i["viscosity"]
        result = laminar_flat_plate_incompressible(reynolds)
        drag = 2.0 * 0.5 * i["density"] * i["velocity"] ** 2 * i["area"] * result.average_skin_friction
        tolerance = case["relative_tolerance"]
        self.assertLessEqual(_relative(reynolds, case["expected"]["reynolds"]), tolerance)
        self.assertLessEqual(_relative(result.average_skin_friction, case["expected"]["average_skin_friction"]), tolerance)
        self.assertLessEqual(_relative(drag, case["expected"]["total_drag_n"]), tolerance)

    def test_example_19_1a(self) -> None:
        case = self._cases()["example_19_1a_turbulent_incompressible"]
        i = case["inputs"]
        result = turbulent_flat_plate_incompressible(i["reynolds"])
        drag = 2.0 * 0.5 * i["density"] * i["velocity"] ** 2 * i["area"] * result.average_skin_friction
        tolerance = case["relative_tolerance"]
        self.assertLessEqual(_relative(result.average_skin_friction, case["expected"]["average_skin_friction"]), tolerance)
        self.assertLessEqual(_relative(drag, case["expected"]["total_drag_n"]), tolerance)

    def test_reference_temperature_examples(self) -> None:
        for name in ("example_18_2_eckert_reference_temperature", "example_18_3_meador_smart_laminar",
                     "example_19_2_eckert_turbulent", "example_19_3_meador_smart_turbulent"):
            case = self._cases()[name]
            i = dict(case["inputs"])
            area = i.pop("area")
            viscosity = SutherlandViscosity(i.pop("mu0"), i.pop("t0"), 110.0)
            result = compressible_flat_plate_reference_temperature(viscosity=viscosity, **i)
            drag = 2.0 * area * result.friction_drag_per_area_pa
            with self.subTest(case=name):
                self.assertLessEqual(_relative(drag, case["expected"]["total_drag_n"]), case["relative_tolerance"])
                if "reference_temperature_k" in case["expected"]:
                    self.assertLessEqual(
                        _relative(result.reference_temperature_k, case["expected"]["reference_temperature_k"]),
                        case["relative_tolerance"],
                    )
                self.assertAlmostEqual(result.wall_heat_flux_w_m2, 0.0, delta=1.0e-9)

    def test_cold_wall_heats_and_hot_wall_cools(self) -> None:
        base = {"edge_temperature_k": 220.0, "edge_pressure_pa": 2000.0, "edge_velocity_m_s": 1500.0, "x_m": 1.0,
                    "gas_constant_j_kg_k": 287.0}
        cold = compressible_flat_plate_reference_temperature(wall_temperature_k=300.0, **base)
        hot = compressible_flat_plate_reference_temperature(wall_temperature_k=2.0 * cold.adiabatic_wall_temperature_k,
                                                            **base)
        self.assertGreater(cold.wall_heat_flux_w_m2, 0.0)
        self.assertLess(hot.wall_heat_flux_w_m2, 0.0)
        self.assertTrue(any(w.code == "calorically_perfect_limit" for w in hot.warnings))

    def test_reference_temperature_forms(self) -> None:
        self.assertAlmostEqual(reference_temperature_ratio(2.94, 2.467), 1.0 + 0.032 * 2.94**2 + 0.58 * 1.467)
        self.assertAlmostEqual(reference_temperature_ratio(2.94, 2.467, "meador_smart_laminar", 0.843),
                               0.45 + 0.55 * 2.467 + 0.16 * 0.843 * 0.2 * 2.94**2)
        with self.assertRaises(InputError):
            reference_temperature_ratio(2.0, 1.0, "meador_smart_laminar")
        with self.assertRaises(InputError):
            reference_temperature_ratio(2.0, 1.0, "unknown")

    def test_recovery_and_reynolds_analogy(self) -> None:
        self.assertAlmostEqual(adiabatic_wall_temperature(288.0, 2.94, 0.843) / 288.0, 1.0 + 0.843 * 0.2 * 2.94**2)
        self.assertAlmostEqual(reynolds_analogy_stanton(1.0e-3, 1.0), 5.0e-4)
        self.assertAlmostEqual(reynolds_analogy_stanton(1.0e-3, 0.71, "couette"), 5.0e-4 / 0.71)
        self.assertEqual(flat_plate_regime(1.0e5, 5.0e5), "laminar")
        self.assertEqual(flat_plate_regime(1.0e6, 5.0e5), "turbulent")


class VanDriestIITests(unittest.TestCase):
    def test_hopkins_table_2_ratios(self) -> None:
        table = REFERENCE["van_driest_ii_hopkins_table_2"]
        worst = 0.0
        for row in table["rows"]:
            taw_088 = row["te_k"] * (1.0 + 0.88 * 0.2 * row["mach"] ** 2)
            wall = row["tw_taw"] * taw_088
            computed = []
            for reynolds in (1.0e5, 1.0e9):
                common = {"mach": row["mach"], "edge_temperature_k": row["te_k"], "reynolds_x": reynolds,
                              "wall_temperature_k": wall}
                reference = van_driest_ii(recovery_factor=0.88, **common).local_skin_friction
                computed.append((reference, van_driest_ii(recovery_factor=1.0, **common).local_skin_friction,
                                 van_driest_ii(recovery_factor=1.0, viscosity_law="power", viscosity_exponent=0.76,
                                               **common).local_skin_friction))
            ratios = [computed[0][0] / computed[0][1], computed[1][0] / computed[1][1],
                      computed[0][0] / computed[0][2], computed[1][0] / computed[1][2]]
            for value, expected in zip(ratios, row["ratios"]):
                worst = max(worst, _relative(value, expected))
                with self.subTest(row=row, expected=expected):
                    self.assertLessEqual(_relative(value, expected), table["relative_tolerance"])
        self.assertLess(worst, 0.017)

    def test_incompressible_limit_reduces_to_karman_schoenherr(self) -> None:
        result = van_driest_ii(mach=0.0, edge_temperature_k=288.0, reynolds_x=1.0e7)
        self.assertAlmostEqual(result.f_c, 1.0)
        self.assertAlmostEqual(0.242 / math.sqrt(result.average_skin_friction),
                               math.log10(1.0e7 * result.average_skin_friction), places=9)

    def test_compressibility_lowers_skin_friction(self) -> None:
        low = van_driest_ii(mach=1.0, edge_temperature_k=220.0, reynolds_x=1.0e7).average_skin_friction
        high = van_driest_ii(mach=6.0, edge_temperature_k=220.0, reynolds_x=1.0e7).average_skin_friction
        self.assertLess(high, low)

    def test_rejects_out_of_range(self) -> None:
        with self.assertRaises(DomainError):
            van_driest_ii(mach=12.0, edge_temperature_k=220.0, reynolds_x=1.0e7)
        with self.assertRaises(InputError):
            van_driest_ii(mach=2.0, edge_temperature_k=220.0, reynolds_x=1.0e7, wall_temperature_k=300.0,
                          wall_to_adiabatic_temperature=0.5)


class CouetteTests(unittest.TestCase):
    def test_example_16_1(self) -> None:
        case = REFERENCE["anderson_couette_examples"][0]
        c = case["conversion_factors"]
        i = case["inputs_english"]
        ft, lbf, slug, rankine = c["ft_m"], c["lbf_n"], c["slug_kg"], c["R_K"]
        cp = i["cp_ft_lb_slug_R"] * ft * lbf / slug / rankine
        result = incompressible_couette(
            upper_plate_velocity_m_s=i["ue_ft_s"] * ft, gap_m=i["gap_in"] * c["in_m"],
            viscosity_pa_s=i["mu_slug_ft_s"] * slug / ft, density_kg_m3=1.2, specific_heat_j_kg_k=cp,
            prandtl=i["prandtl"], upper_wall_temperature_k=i["temperature_R"] * rankine,
            lower_wall_temperature_k=i["temperature_R"] * rankine,
        )
        computed = {
            "u_mid_ft_s": result.profile(0.5)[0] / ft,
            "shear_lb_ft2": result.shear_stress_pa / (lbf / ft**2),
            "max_temperature_R": result.maximum_temperature_k / rankine,
            "wall_heat_flux_ft_lb_ft2_s": result.lower_wall_heat_flux_w_m2 / (lbf * ft / ft**2),
            "adiabatic_wall_temperature_R": result.adiabatic_wall_temperature_k / rankine,
        }
        for field, expected in case["expected_english"].items():
            with self.subTest(field=field):
                self.assertLessEqual(_relative(computed[field], expected), case["relative_tolerance"][field])
        self.assertAlmostEqual(result.stanton / result.skin_friction, 0.5 / i["prandtl"])
        self.assertAlmostEqual(result.upper_wall_heat_flux_w_m2, -result.lower_wall_heat_flux_w_m2)

    def test_example_16_2(self) -> None:
        case = REFERENCE["anderson_couette_examples"][1]
        i = case["inputs"]
        ue = i["mach"] * math.sqrt(i["gamma"] * i["gas_constant_j_kg_k"] * i["temperature_k"])
        result = compressible_couette(
            upper_plate_velocity_m_s=ue, gap_m=i["gap_m"], pressure_pa=i["pressure_pa"],
            upper_wall_temperature_k=i["temperature_k"], lower_wall_temperature_k=i["temperature_k"],
            gas_constant_j_kg_k=i["gas_constant_j_kg_k"], gamma=i["gamma"], prandtl=i["prandtl"],
            viscosity=SutherlandViscosity(1.7894e-5, 288.0, 110.0),
        )
        tolerance = case["relative_tolerance"]
        self.assertLessEqual(_relative(ue, case["expected"]["upper_velocity_m_s"]), tolerance)
        flux_per_shear = result.lower_wall_heat_flux_w_m2 / result.shear_stress_pa
        self.assertAlmostEqual(flux_per_shear, ue / 2.0, delta=1.0e-9 * ue)
        self.assertLessEqual(_relative(i["given_shear_pa"] * flux_per_shear, case["expected"]["wall_heat_flux_w_m2"]),
                             tolerance)
        self.assertLessEqual(_relative(result.adiabatic_wall_temperature_k,
                                       case["expected"]["adiabatic_wall_temperature_k"]), tolerance)
        self.assertAlmostEqual(result.stanton / result.skin_friction, 0.5 / i["prandtl"], places=9)

    def test_compressible_reduces_to_incompressible_at_low_speed(self) -> None:
        law = SutherlandViscosity()
        result = compressible_couette(upper_plate_velocity_m_s=1.0, gap_m=1.0e-3, pressure_pa=1.0e5,
                                      upper_wall_temperature_k=300.0, lower_wall_temperature_k=300.0,
                                      gas_constant_j_kg_k=287.0, viscosity=law)
        self.assertLess(_relative(result.shear_stress_pa, law.viscosity(300.0) * 1.0 / 1.0e-3), 1.0e-6)

    def test_adiabatic_lower_wall(self) -> None:
        result = compressible_couette(upper_plate_velocity_m_s=2000.0, gap_m=1.0e-3, pressure_pa=1.0e5,
                                      upper_wall_temperature_k=250.0, lower_wall_temperature_k=None,
                                      gas_constant_j_kg_k=287.0)
        self.assertEqual(result.lower_wall_heat_flux_w_m2, 0.0)
        self.assertAlmostEqual(result.maximum_temperature_k, result.adiabatic_wall_temperature_k)
        self.assertAlmostEqual(result.profile[0][2], result.adiabatic_wall_temperature_k)
        self.assertAlmostEqual(result.profile[-1][0], 1.0)
        self.assertTrue(result.warnings)

    def test_quadrature_converged(self) -> None:
        kwargs = {"upper_plate_velocity_m_s": 3000.0, "gap_m": 1.0e-3, "pressure_pa": 1.0e5, "upper_wall_temperature_k": 250.0,
                      "lower_wall_temperature_k": 300.0, "gas_constant_j_kg_k": 287.0}
        default = compressible_couette(**kwargs).shear_stress_pa
        fine = compressible_couette(intervals=1600, **kwargs).shear_stress_pa
        self.assertLess(_relative(default, fine), 1.0e-9)


class StagnationHeatingTests(unittest.TestCase):
    def test_anderson_example_14_2(self) -> None:
        case = REFERENCE["stagnation_heating"][0]
        i = case["inputs"]
        v = i["velocity_m_s"]
        q = tauber_stagnation_heating(freestream_density_kg_m3=i["freestream_density_kg_m3"], velocity_m_s=v,
                                      nose_radius_m=i["nose_radius_m"],
                                      wall_enthalpy_j_kg=i["cp_j_kg_k"] * i["wall_temperature_k"],
                                      total_enthalpy_j_kg=0.5 * v * v)
        self.assertLessEqual(_relative(q / 1.0e4, case["expected"]["heat_flux_w_cm2"]), case["relative_tolerance"])

    def test_sutton_graves_mixture_rule_against_table_ii(self) -> None:
        data = REFERENCE["sutton_graves"]
        for case in data["cases"]:
            with self.subTest(gas=case["gas"]):
                self.assertEqual(SUTTON_GRAVES_TABLE_II[case["gas"]], case["table_ii_K"])
                k = sutton_graves_coefficient(case["mass_fractions"])
                self.assertLessEqual(_relative(k, case["table_ii_K"]), data["eq42_vs_table_ii_relative_tolerance"])

    def test_sutton_graves_heating_definition(self) -> None:
        q = sutton_graves_heating(stagnation_pressure_pa=101325.0 * 4.0, nose_radius_m=1.0,
                                  stagnation_enthalpy_j_kg=2.0e7, wall_enthalpy_j_kg=0.0, coefficient=0.1)
        self.assertAlmostEqual(q, 0.1 * 2.0 * 2.0e7)
        with self.assertRaises(DomainError):
            sutton_graves_coefficient({"O2": 0.5, "N2": 0.4})
        with self.assertRaises(InputError):
            sutton_graves_coefficient({"Xe": 1.0})

    def test_boundary_layer_stagnation_forms(self) -> None:
        gradient = newtonian_stagnation_velocity_gradient(0.5, 2.0e4, 1.0e3, 0.05)
        self.assertAlmostEqual(gradient, math.sqrt(2.0 * 1.9e4 / 0.05) / 0.5)
        common = {"prandtl": 0.71, "edge_density_kg_m3": 0.05, "edge_viscosity_pa_s": 1.0e-4,
                      "velocity_gradient_1_s": gradient, "adiabatic_wall_enthalpy_j_kg": 5.0e6, "wall_enthalpy_j_kg": 3.0e5}
        sphere = stagnation_point_heating_boundary_layer(geometry="sphere", **common)
        cylinder = stagnation_point_heating_boundary_layer(geometry="cylinder", **common)
        self.assertAlmostEqual(sphere / cylinder, 0.763 / 0.57)
        # q ~ R^-1/2 (Anderson Eq. (18.83)) through the Newtonian gradient.
        half = newtonian_stagnation_velocity_gradient(1.0, 2.0e4, 1.0e3, 0.05)
        common["velocity_gradient_1_s"] = half
        self.assertAlmostEqual(stagnation_point_heating_boundary_layer(geometry="sphere", **common) / sphere,
                               1.0 / math.sqrt(2.0))

    def test_tauber_closed_forms(self) -> None:
        common = {"freestream_density_kg_m3": 1.0e-3, "velocity_m_s": 3000.0}
        cone = tauber_laminar_cone_heating(distance_m=1.0, cone_half_angle_deg=10.0, wall_enthalpy_j_kg=3.0e5,
                                           freestream_enthalpy_j_kg=2.0e5, **common)
        plate = tauber_laminar_flat_plate_heating(distance_m=1.0, plate_angle_deg=10.0, wall_enthalpy_j_kg=3.0e5,
                                                  freestream_enthalpy_j_kg=2.0e5, **common)
        self.assertAlmostEqual(cone / plate, 4.03 / 2.42)
        haw = 2.0e5 + 0.40 * 3000.0**2
        self.assertAlmostEqual(
            cone, 4.03e-5 * math.sqrt(1.0e-3 * math.cos(math.radians(10.0))) * 3000.0**3.2
            * math.sin(math.radians(10.0)) * (1.0 - 3.0e5 / haw))
        cylinder = tauber_swept_cylinder_heating(cylinder_radius_m=0.1, sweep_deg=0.0, wall_enthalpy_j_kg=3.0e5,
                                                 freestream_enthalpy_j_kg=2.0e5, **common)
        self.assertAlmostEqual(tauber_leading_edge_heating(cylinder, plate, 0.0, 0.0), cylinder)
        turbulent = tauber_turbulent_flat_plate_heating(distance_from_transition_m=1.0, plate_angle_deg=10.0,
                                                        wall_temperature_k=555.0, wall_enthalpy_j_kg=0.0,
                                                        total_enthalpy_j_kg=4.7e6, **common)
        d = math.radians(10.0)
        expected = 3.72e-4 * (1.0e-3 * math.sin(d) ** 2 * math.cos(d) ** 2.22) ** 0.8 * 3000.0**3.37 * 0.9
        self.assertAlmostEqual(turbulent / expected, 1.0)
        self.assertAlmostEqual(tauber_turbulent_cone_heating(turbulent) / turbulent, 1.15)
        with self.assertRaises(DomainError):
            tauber_turbulent_flat_plate_heating(distance_from_transition_m=1.0, plate_angle_deg=10.0,
                                                wall_temperature_k=555.0, wall_enthalpy_j_kg=0.0,
                                                total_enthalpy_j_kg=4.7e6, freestream_density_kg_m3=1.0e-3,
                                                velocity_m_s=1000.0)


class AtmosphereTests(unittest.TestCase):
    def test_ussa1976_density_table_points(self) -> None:
        data = REFERENCE["ussa1976_density"]
        for point in data["points"]:
            with self.subTest(z=point["z_m"]):
                self.assertLessEqual(_relative(ussa1976_density(point["z_m"]), point["density_kg_m3"]),
                                     data["relative_tolerance"])
                self.assertAlmostEqual(ussa1976_geopotential_altitude(point["z_m"]), point["geopotential_m"],
                                       delta=data["geopotential_absolute_tolerance_m"])

    def test_ussa1976_ceiling_fails_closed(self) -> None:
        self.assertAlmostEqual(USSA1976_MAX_GEOMETRIC_ALTITUDE_M, 86000.0, delta=1.0)
        with self.assertRaises(DomainError):
            ussa1976_density(90000.0)
        with self.assertRaises(DomainError):
            integrate_planar_entry(entry_altitude_m=120000.0, entry_velocity_m_s=7000.0,
                                   entry_flight_path_angle_deg=-5.0, mass_kg=100.0, drag_coefficient=1.0,
                                   reference_area_m2=1.0)


def _allen_eggers_example():
    data = REFERENCE["allen_eggers"]["example_atmosphere_english"]
    ft, slug = 0.3048, 14.59390293720636
    atmosphere = ExponentialAtmosphere(data["rho0_slug_ft3"] * slug / ft**3, data["beta_1_per_ft"] / ft)
    diameter = data["sphere_diameter_ft"] * ft
    mass = 7870.0 * math.pi * diameter**3 / 6.0  # iron sphere; density chosen here, not from the report
    area = math.pi * diameter**2 / 4.0
    return data, atmosphere, mass, area, data["entry_speed_ft_s"] * ft, data["entry_altitude_miles"] * 1609.344


class AllenEggersTests(unittest.TestCase):
    def test_published_velocity_ratios(self) -> None:
        data, atmosphere, mass, area, speed, _ = _allen_eggers_example()
        entry = allen_eggers_entry(entry_velocity_m_s=speed, entry_angle_deg=90.0, mass_kg=mass,
                                   drag_coefficient=data["drag_coefficient"], reference_area_m2=area,
                                   atmosphere=atmosphere)
        ratios = REFERENCE["allen_eggers"]["expected_velocity_ratios"]
        tolerance = REFERENCE["allen_eggers"]["absolute_tolerance"]
        self.assertAlmostEqual(entry.peak_deceleration_velocity_m_s / speed, ratios["peak_deceleration"],
                               delta=tolerance)
        self.assertAlmostEqual(entry.peak_average_heating_velocity_m_s / speed, ratios["peak_average_heating"],
                               delta=tolerance)
        self.assertAlmostEqual(entry.peak_stagnation_heating_velocity_m_s / speed, ratios["peak_stagnation_heating"],
                               delta=tolerance)
        self.assertAlmostEqual(entry.velocity(entry.peak_deceleration_altitude_m), entry.peak_deceleration_velocity_m_s)

    def test_numerical_integration_reproduces_closed_form(self) -> None:
        data, atmosphere, mass, area, speed, start = _allen_eggers_example()
        for angle in (90.0, 30.0):
            entry = allen_eggers_entry(entry_velocity_m_s=speed, entry_angle_deg=angle, mass_kg=mass,
                                       drag_coefficient=data["drag_coefficient"], reference_area_m2=area,
                                       atmosphere=atmosphere, stagnation_heating_constant=1.83e-4,
                                       nose_radius_m=0.1524)
            trajectory = integrate_planar_entry(
                entry_altitude_m=start, entry_velocity_m_s=entry.velocity(start), entry_flight_path_angle_deg=-angle,
                mass_kg=mass, drag_coefficient=data["drag_coefficient"], reference_area_m2=area,
                atmosphere=atmosphere, gravity="none", curvature=False, time_step_s=0.005,
                final_altitude_m=0.0, final_velocity_m_s=0.3 * speed, stagnation_heating_constant=1.83e-4,
                nose_radius_m=0.1524,
            )
            with self.subTest(angle=angle):
                self.assertEqual(trajectory.termination, "velocity")
                worst = max(_relative(s.velocity_m_s, entry.velocity(s.altitude_m)) for s in trajectory.samples)
                self.assertLess(worst, 1.0e-6)
                self.assertLess(_relative(trajectory.peak_deceleration.deceleration_g, entry.peak_deceleration_g), 1.0e-3)
                self.assertLess(_relative(trajectory.peak_heating.stagnation_heating_w_m2,
                                          entry.peak_stagnation_heating_w_m2), 1.0e-3)
                self.assertAlmostEqual(trajectory.peak_deceleration.altitude_m, entry.peak_deceleration_altitude_m,
                                       delta=50.0)

    def test_heavy_body_peaks_at_sea_level(self) -> None:
        atmosphere = ExponentialAtmosphere(1.2, 1.0 / 7000.0)
        entry = allen_eggers_entry(entry_velocity_m_s=7000.0, entry_angle_deg=60.0, mass_kg=1.0e6,
                                   drag_coefficient=0.1, reference_area_m2=0.01, atmosphere=atmosphere)
        self.assertLess(entry.peak_deceleration_altitude_m, 0.0)
        self.assertAlmostEqual(entry.peak_deceleration_g, entry.deceleration_g(0.0))
        self.assertAlmostEqual(entry.peak_deceleration_velocity_m_s, entry.impact_velocity_m_s)

    def test_total_heat_input(self) -> None:
        self.assertAlmostEqual(total_heat_input_reynolds_analogy(mass_kg=10.0, entry_velocity_m_s=100.0,
                                                                 final_velocity_m_s=0.0, friction_to_drag_ratio=0.5),
                               0.5 * 0.5 * 0.5 * 10.0 * 100.0**2)

    def test_ussa1976_trajectory_with_gravity(self) -> None:
        trajectory = integrate_planar_entry(entry_altitude_m=85000.0, entry_velocity_m_s=7000.0,
                                            entry_flight_path_angle_deg=-20.0, mass_kg=500.0, drag_coefficient=1.0,
                                            reference_area_m2=1.0, time_step_s=0.05, final_velocity_m_s=500.0,
                                            nose_radius_m=0.5, stagnation_heating_constant=1.83e-4)
        self.assertEqual(trajectory.termination, "velocity")
        speeds = [s.velocity_m_s for s in trajectory.samples]
        self.assertGreater(trajectory.peak_deceleration.deceleration_g, 1.0)
        self.assertLess(trajectory.peak_heating.altitude_m, 85000.0)
        self.assertGreater(trajectory.integrated_heat_load_j_m2, 0.0)
        self.assertLess(speeds[-1], speeds[0])



class EntryTerminalStepTests(unittest.TestCase):
    def test_coarse_step_reaching_the_ground_terminates_at_the_crossing(self) -> None:
        from rocket_propulsion.aerodynamics.entry import integrate_planar_entry

        trajectory = integrate_planar_entry(
            entry_altitude_m=85_000.0, entry_velocity_m_s=7_000.0, entry_flight_path_angle_deg=-20.0,
            mass_kg=1_000.0, drag_coefficient=1.0, reference_area_m2=3.0, time_step_s=0.5,
            final_velocity_m_s=1.0,
        )
        self.assertEqual(trajectory.termination, "altitude")
        self.assertEqual(trajectory.samples[-1].altitude_m, 0.0)


if __name__ == "__main__":
    unittest.main()
