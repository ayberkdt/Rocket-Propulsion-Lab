"""Source-traceable checks of the launch-vehicle aerodynamics modules.

Reference values and page locators are in
``tests/reference_data/launch_vehicle_aero.json``.
"""

from __future__ import annotations

import json
import math
import unittest
from pathlib import Path

from rocket_propulsion.aerodynamics.aero_database import (
    RANGE_METHOD_D2,
    AeroDatabase,
    RangeUncertainty,
    UncertaintyTerm,
    combine_uncertainty_terms,
    coverage_bound,
    fit_linear,
    fit_moment_center_of_pressure,
    net_modeling_uncertainty,
    range_method_sigma,
    repeatability_uncertainty,
)
from rocket_propulsion.aerodynamics.launch_vehicle import (
    SP8035_A1_M,
    SP8035_A2_M,
    barrowman_base_drag,
    barrowman_skin_friction,
    body_friction_drag_coefficient,
    combine_wind_loads,
    dynamic_pressure_at_altitude,
    dynamic_pressure_from_mach,
    hama_turbulent_skin_friction,
    integrate_line_load,
    loads_from_coefficients,
    maximum_dynamic_pressure,
    schoenherr_turbulent_skin_friction,
    sp8035_discrete_gust_speed,
    subsonic_forebody_pressure_drag_coefficient,
    wind_induced_incidence,
)
from rocket_propulsion.aerodynamics.slender_body import (
    AxisymmetricBody,
    NormalForceComponent,
    TrapezoidalFinSet,
    barrowman_fin_center_of_pressure,
    barrowman_fin_normal_force_slope,
    body_geometry,
    combine_normal_force_components,
    cone_nose,
    conical_frustum,
    crossflow_line_load,
    cylinder_section,
    fin_in_presence_of_body,
    interference_factors,
    slender_body_normal_force,
    tangent_ogive_nose,
    viscous_crossflow_coefficients,
)
from rocket_propulsion.core.constants import (
    AIR_SPECIFIC_GAS_CONSTANT_J_KG_K,
    USSA1976_SEA_LEVEL_PRESSURE_PA,
    USSA1976_SEA_LEVEL_TEMPERATURE_K,
)
from rocket_propulsion.core.errors import DomainError, InputError

DATA = json.loads(
    (Path(__file__).parent / "reference_data" / "launch_vehicle_aero.json").read_text(encoding="utf-8")
)
INCH = 0.0254


def tr1048_body(cylinder_in: float) -> AxisymmetricBody:
    """Return an Allen-Perkins ogive-cylinder (NACA TR 1048 figure 6)."""

    inputs = DATA["tr1048_bodies"]["inputs_in"]
    radius = 0.5 * inputs["diameter"] * INCH
    return AxisymmetricBody(
        (
            tangent_ogive_nose(inputs["nose_length"] * INCH, radius),
            cylinder_section(cylinder_in * INCH, radius),
        )
    )


class Tr1048Tests(unittest.TestCase):
    def test_body_geometry_reproduces_figure_6(self) -> None:
        ref = DATA["tr1048_bodies"]
        inputs = ref["inputs_in"]
        tol = ref["tolerance"]
        for key, cylinder in (("I", inputs["cylinder_body_I"]), ("II", inputs["cylinder_body_II"])):
            geometry = body_geometry(tr1048_body(cylinder))
            self.assertAlmostEqual(
                geometry.fineness_ratio, ref["expected"][f"fineness_body_{key}"], delta=tol["fineness_abs"]
            )
        length = inputs["nose_length"]
        radius = 0.5 * inputs["diameter"]
        rho = (length * length + radius * radius) / (2.0 * radius)
        self.assertAlmostEqual(rho, ref["expected"]["generating_radius_in"], delta=tol["generating_radius_abs_in"])

    def test_figure_7_lift_at_8_degrees(self) -> None:
        ref = DATA["tr1048_fig7_body_I"]
        body = tr1048_body(DATA["tr1048_bodies"]["inputs_in"]["cylinder_body_I"])
        geometry = body_geometry(body)
        crossflow = ref["inputs"]["cdc"] * ref["inputs"]["eta"]
        kwargs = {
            "angle_of_attack_deg": ref["inputs"]["angle_of_attack_deg"],
            "crossflow_drag_coefficient": crossflow,
            "reference_area_m2": geometry.base_area_m2,
            "reference_length_m": DATA["tr1048_bodies"]["inputs_in"]["diameter"] * INCH,
        }
        full = viscous_crossflow_coefficients(body, **kwargs)
        small = viscous_crossflow_coefficients(body, small_angle=True, **kwargs)
        tolerance = ref["reading_tolerance_abs"]
        read = ref["figure_read"]
        self.assertAlmostEqual(full.lift_coefficient, read["dashed_lift_coefficient"], delta=tolerance)
        self.assertAlmostEqual(small.lift_coefficient, read["dashed_lift_coefficient"], delta=tolerance)
        self.assertAlmostEqual(
            small.potential_normal_force_coefficient, read["potential_lift_coefficient"], delta=tolerance
        )
        # Measured deviation from the experiment: documented, not asserted to agree.
        experiment = viscous_crossflow_coefficients(
            body, **{**kwargs, "angle_of_attack_deg": read["experiment_alpha_deg"]}
        )
        deviation = experiment.lift_coefficient / read["experiment_lift_coefficient"] - 1.0
        self.assertTrue(0.0 < deviation < 0.30, deviation)

    def test_equation_14_internal_consistency(self) -> None:
        body = tr1048_body(11.5)
        geometry = body_geometry(body)
        result = viscous_crossflow_coefficients(
            body,
            angle_of_attack_deg=15.0,
            crossflow_drag_coefficient=0.9,
            reference_area_m2=geometry.base_area_m2,
            reference_length_m=0.75 * INCH,
            moment_reference_m=0.3,
        )
        alpha = math.radians(15.0)
        # Lift is the normal force rotated with the potential part inclined at alpha/2.
        self.assertAlmostEqual(
            result.lift_coefficient,
            result.potential_normal_force_coefficient + result.viscous_normal_force_coefficient * math.cos(alpha),
            places=12,
        )
        # Integrating the printed cross-force distribution recovers C_N and the CP.
        stations = [geometry.length_m * i / 4000 for i in range(4001)]
        loads = [
            crossflow_line_load(
                body,
                x_m=x,
                angle_of_attack_deg=15.0,
                dynamic_pressure_pa=1.0,
                crossflow_drag_coefficient=0.9,
            )
            for x in stations
        ]
        distribution = integrate_line_load(stations_m=stations, line_load_n_m=loads)
        self.assertAlmostEqual(
            distribution.total_normal_force_n / geometry.base_area_m2, result.normal_force_coefficient, places=3
        )
        self.assertAlmostEqual(distribution.center_of_pressure_m, result.center_of_pressure_m, places=4)

    def test_angle_domain(self) -> None:
        body = tr1048_body(11.5)
        with self.assertRaises(DomainError):
            viscous_crossflow_coefficients(
                body, angle_of_attack_deg=95.0, crossflow_drag_coefficient=1.0,
                reference_area_m2=1.0, reference_length_m=1.0,
            )


class BarrowmanBodyTests(unittest.TestCase):
    def test_nose_slope_and_cone_cp(self) -> None:
        ref = DATA["barrowman"]
        cone = AxisymmetricBody((cone_nose(1.2, 0.1),))
        combined = combine_normal_force_components(
            slender_body_normal_force(cone, reference_area_m2=math.pi * 0.01)
        )
        self.assertAlmostEqual(combined.normal_force_slope_per_rad, ref["nose_normal_force_slope"]["expected_on_base_area"], places=12)
        self.assertAlmostEqual(
            combined.center_of_pressure_m / 1.2, ref["cone_center_of_pressure"]["expected_fraction_of_length"], places=12
        )

    def test_boattail_is_negative_and_body_is_shape_independent(self) -> None:
        body = AxisymmetricBody(
            (
                tangent_ogive_nose(0.5, 0.1),
                cylinder_section(2.0, 0.1),
                conical_frustum(0.3, 0.1, 0.07),
            )
        )
        parts = slender_body_normal_force(body, reference_area_m2=math.pi * 0.01)
        self.assertLess(parts[2].normal_force_slope_per_rad, 0.0)
        self.assertIsNone(parts[1].center_of_pressure_m)
        total = combine_normal_force_components(parts)
        # Equation 3-65: only the end areas matter.
        self.assertAlmostEqual(total.normal_force_slope_per_rad, 2.0 * 0.07**2 / 0.1**2, places=10)
        # Equation 3-89 for the whole body: X = l - V / A(l).
        geometry = body_geometry(body)
        self.assertAlmostEqual(
            total.center_of_pressure_m,
            geometry.length_m - geometry.volume_m3 / geometry.base_area_m2,
            places=9,
        )

    def test_equation_3_90_normalisation(self) -> None:
        expected_ratio = DATA["barrowman"]["ogive_volume_eq_3_90"]["observed_ratio_integrated_to_rhs_times_dN_A"]
        for length, radius in ((4.314, 0.375), (3.0, 0.5), (1.0, 0.5)):
            body = AxisymmetricBody((tangent_ogive_nose(length, radius),))
            volume = body_geometry(body).volume_m3
            diameter = 2.0 * radius
            f_n = length / diameter
            rhs = (
                f_n * (f_n**2 + 0.25) ** 2
                - f_n**3 / 3.0
                - (f_n**2 - 0.25) * (f_n**2 + 0.25) ** 2 * math.asin(f_n / (f_n**2 + 0.25))
            )
            self.assertAlmostEqual(volume / (rhs * diameter * math.pi * radius**2), expected_ratio, places=6)

    def test_discontinuous_body_is_refused(self) -> None:
        with self.assertRaises(InputError):
            AxisymmetricBody((cone_nose(1.0, 0.1), cylinder_section(1.0, 0.2)))
        with self.assertRaises(InputError):
            AxisymmetricBody((cylinder_section(1.0, 0.1),))


class FinTests(unittest.TestCase):
    def test_interference_factor_limits_and_figure_2(self) -> None:
        ref = DATA["tr1307_kwb_figure2"]
        for ratio, value in ref["figure_read"].items():
            r = float(ratio)
            factors = interference_factors(body_radius_m=r, exposed_semispan_m=1.0 - r)
            self.assertAlmostEqual(factors.tail_in_presence_of_body, value, delta=ref["reading_tolerance_abs"])
        small = interference_factors(body_radius_m=1.0e-7, exposed_semispan_m=1.0)
        large = interference_factors(body_radius_m=0.9999, exposed_semispan_m=1.0e-4)
        self.assertAlmostEqual(small.tail_in_presence_of_body, ref["limits"]["kwb_r0"], places=5)
        self.assertAlmostEqual(small.body_in_presence_of_tail, ref["limits"]["kbw_r0"], places=5)
        self.assertAlmostEqual(large.tail_in_presence_of_body, ref["limits"]["kwb_r1"], places=3)
        self.assertAlmostEqual(large.body_in_presence_of_tail, ref["limits"]["kbw_r1"], places=3)
        with self.assertRaises(DomainError):
            interference_factors(body_radius_m=1.0, exposed_semispan_m=1.0e-7)

    def test_fin_cp_and_slope(self) -> None:
        fins = TrapezoidalFinSet(4, 0.3, 0.1, 0.15, 0.2, 2.0, 0.05)
        x_cp, y_cp = barrowman_fin_center_of_pressure(fins)
        # Rectangular-fin check: x_t = 0, c_t = c_r gives the quarter chord.
        square = TrapezoidalFinSet(3, 0.2, 0.2, 0.1, 0.0, 1.0, 0.05)
        self.assertAlmostEqual(barrowman_fin_center_of_pressure(square)[0], 1.0 + 0.05, places=12)
        self.assertGreater(x_cp, 2.0)
        self.assertGreater(y_cp, 0.05)
        slope = barrowman_fin_normal_force_slope(fins, mach=0.5, aspect_ratio=1.5, reference_area_m2=0.01)
        beta = math.sqrt(0.75)
        ratio = beta * 1.5 / math.cos(fins.midchord_sweep_rad)
        expected = 4 * math.pi * 1.5 * (fins.panel_area_m2 / 0.01) / (2.0 + math.sqrt(4.0 + ratio**2))
        self.assertAlmostEqual(slope, expected, places=12)
        component = fin_in_presence_of_body(fins, mach=0.5, aspect_ratio=1.5, reference_area_m2=0.01)
        self.assertGreater(component.normal_force_slope_per_rad, slope)
        with self.assertRaises(DomainError):
            barrowman_fin_normal_force_slope(fins, mach=1.2, aspect_ratio=1.5, reference_area_m2=0.01)

    def test_combine_requires_cp(self) -> None:
        with self.assertRaises(InputError):
            combine_normal_force_components([NormalForceComponent("x", 1.0, None)])
        with self.assertRaises(DomainError):
            combine_normal_force_components([NormalForceComponent("x", 0.0, None)])


class DragTests(unittest.TestCase):
    def test_hama_within_two_percent_of_schoenherr(self) -> None:
        ref = DATA["barrowman"]["hama_vs_schoenherr"]
        for reynolds in ref["reynolds"]:
            hama = hama_turbulent_skin_friction(reynolds)
            schoenherr = schoenherr_turbulent_skin_friction(reynolds)
            self.assertLess(abs(hama / schoenherr - 1.0), ref["max_relative_difference"])
            # Schoenherr's equation 4-3 is satisfied by the solved value.
            self.assertAlmostEqual(math.log10(reynolds * schoenherr), 0.242 / math.sqrt(schoenherr), places=9)

    def test_skin_friction_branches(self) -> None:
        laminar = barrowman_skin_friction(reynolds_number=1.0e5, mach=0.3)
        self.assertEqual(laminar.regime, "laminar")
        self.assertAlmostEqual(laminar.compressible_coefficient, DATA["barrowman"]["laminar"]["expected"], places=6)
        transitional = barrowman_skin_friction(reynolds_number=1.0e7, mach=0.5)
        self.assertAlmostEqual(
            transitional.incompressible_coefficient, hama_turbulent_skin_friction(1.0e7) - 1700.0 / 1.0e7, places=15
        )
        self.assertAlmostEqual(
            transitional.compressible_coefficient, transitional.incompressible_coefficient * (1 - 0.09 * 0.25), places=15
        )
        supersonic = barrowman_skin_friction(reynolds_number=1.0e7, mach=2.0, fully_turbulent=True)
        self.assertAlmostEqual(
            supersonic.compressible_coefficient, hama_turbulent_skin_friction(1.0e7) / (1 + 0.15 * 4.0) ** 0.58, places=15
        )
        rough = barrowman_skin_friction(reynolds_number=1.0e8, mach=0.0, roughness_to_length=1.0e-4)
        self.assertEqual(rough.regime, "turbulent_rough")
        self.assertAlmostEqual(rough.incompressible_coefficient, 0.032 * 1.0e-4**0.2, places=15)

    def test_base_drag_sonic_values_and_branches(self) -> None:
        ref = DATA["barrowman"]["base_drag_sonic"]
        common = {"base_area_m2": 1.0, "reference_area_m2": 1.0, "wetted_area_m2": 40.0,
                      "compressible_skin_friction": 0.003, "sonic_skin_friction": 0.003}
        sonic = barrowman_base_drag(mach=1.0, **common)
        self.assertAlmostEqual(sonic.sonic_base_drag_coefficient, ref["expected_no_fins"], places=12)
        self.assertAlmostEqual(sonic.base_drag_coefficient, ref["expected_no_fins"], places=12)
        finned = barrowman_base_drag(mach=1.0, fin_trailing_edge_thickness_to_root_chord=0.1, **common)
        self.assertAlmostEqual(finned.sonic_base_drag_coefficient, ref["expected_h_over_c_0p1"], places=12)
        self.assertAlmostEqual(sonic.critical_mach, 0.892 / math.sqrt(0.185), places=12)
        # Subsonic constant K makes equation 4-50 equal C*_DB at M = 1 when C_fc = C*_fc.
        near_sonic = barrowman_base_drag(mach=0.999999, **common)
        self.assertAlmostEqual(near_sonic.base_drag_coefficient, 0.185, places=4)
        high = barrowman_base_drag(mach=3.0, **common)
        self.assertEqual(high.branch, "supersonic_inverse_square")
        self.assertAlmostEqual(high.base_drag_coefficient, 0.7 / 9.0, places=12)

    def test_forebody_pressure_drag_matches_sonic_fit(self) -> None:
        # K is chosen so that equation 4-41 at M = 1 equals the fit 4-43.
        common = {"nose_shape": "tangent_ogive", "nose_fineness_ratio": 3.0, "body_fineness_ratio": 15.0,
                      "wetted_area_m2": 30.0, "reference_area_m2": 0.5,
                      "compressible_skin_friction": 0.0025, "sonic_skin_friction": 0.0025}
        near = subsonic_forebody_pressure_drag_coefficient(mach=1.0 - 1.0e-12, **common)
        # 6.82 in equations 4-45/4-46 is 6/0.88 = 6.818 rounded: 0.03% offset.
        self.assertAlmostEqual(near / (0.88 / (3.0 + 1.0) ** 2.22), 1.0, delta=5.0e-4)
        with self.assertRaises(DomainError):
            subsonic_forebody_pressure_drag_coefficient(mach=1.2, **common)
        with self.assertRaises(DomainError):
            subsonic_forebody_pressure_drag_coefficient(mach=0.5, **{**common, "nose_fineness_ratio": 12.0})

    def test_body_friction(self) -> None:
        self.assertAlmostEqual(
            body_friction_drag_coefficient(compressible_skin_friction=0.003, body_fineness_ratio=10.0,
                                           wetted_area_m2=20.0, reference_area_m2=0.5),
            1.05 * 0.003 * 40.0,
            places=12,
        )


class LoadsTests(unittest.TestCase):
    def test_dynamic_pressure_identity(self) -> None:
        # q = rho u^2/2 with rho = p/(R T) and u = M sqrt(gamma R T).
        p = USSA1976_SEA_LEVEL_PRESSURE_PA
        t = USSA1976_SEA_LEVEL_TEMPERATURE_K
        r = AIR_SPECIFIC_GAS_CONSTANT_J_KG_K
        rho = p / (r * t)
        u = 0.8 * math.sqrt(1.4 * r * t)
        self.assertAlmostEqual(dynamic_pressure_from_mach(static_pressure_pa=p, mach=0.8), 0.5 * rho * u * u, places=6)
        self.assertAlmostEqual(dynamic_pressure_at_altitude(geopotential_altitude_m=0.0, mach=0.8), 0.5 * rho * u * u, places=6)
        peak = maximum_dynamic_pressure(times_s=[0, 1, 2], geopotential_altitudes_m=[0, 8000, 20000], machs=[0.2, 1.2, 2.5])
        self.assertEqual(peak.index, 1)

    def test_sp8035_gust(self) -> None:
        ref = DATA["sp8035_gust"]
        constants = ref["constants"]
        thickness = 130.0
        a1 = SP8035_A1_M
        self.assertEqual(sp8035_discrete_gust_speed(delta_altitude_m=60.0, thickness_m=thickness), constants["amplitude_m_s"])
        self.assertEqual(sp8035_discrete_gust_speed(delta_altitude_m=thickness + 30.0 - a1, thickness_m=thickness), 0.0)
        # Inflection points (half amplitude) are 'th' apart.
        self.assertAlmostEqual(sp8035_discrete_gust_speed(delta_altitude_m=15.0 - a1, thickness_m=thickness), 3.825, places=12)
        self.assertAlmostEqual(sp8035_discrete_gust_speed(delta_altitude_m=thickness + 15.0 - a1, thickness_m=thickness), 3.825, places=12)
        slope = sp8035_discrete_gust_speed(delta_altitude_m=SP8035_A2_M, thickness_m=thickness) / SP8035_A2_M
        self.assertAlmostEqual(slope, ref["tangent_slope_expected"], places=4)
        self.assertAlmostEqual(slope * 10.0, ref["printed_equation_1_slope"], places=3)
        right = sp8035_discrete_gust_speed(delta_altitude_m=SP8035_A2_M + 1e-9, thickness_m=thickness)
        self.assertLess(abs(right - slope * SP8035_A2_M), 1.0e-3)
        with self.assertRaises(DomainError):
            sp8035_discrete_gust_speed(delta_altitude_m=10.0, thickness_m=constants["thickness_max_m"] + 1.0)

    def test_load_combination_and_incidence(self) -> None:
        combined = combine_wind_loads(quasi_steady=-100.0, random_components=(30.0, 40.0))
        self.assertEqual((combined.random_rss, combined.design_load), (50.0, -150.0))
        tail = wind_induced_incidence(ground_speed_m_s=400.0, flight_path_angle_deg=60.0, horizontal_wind_m_s=50.0)
        head = wind_induced_incidence(ground_speed_m_s=400.0, flight_path_angle_deg=60.0, horizontal_wind_m_s=-50.0)
        self.assertGreater(tail.angle_of_attack_deg, 0.0)
        self.assertLess(head.angle_of_attack_deg, 0.0)
        loads = loads_from_coefficients(dynamic_pressure_pa=2.0e4, reference_area_m2=10.0, reference_length_m=8.0,
                                        normal=0.5, pitching=-0.1)
        self.assertEqual((loads.normal_force_n, loads.pitching_moment_n_m), (1.0e5, -1.6e5))

    def test_line_load_integration_exact_for_linear_load(self) -> None:
        stations = [0.0, 0.5, 1.0, 2.0]
        result = integrate_line_load(stations_m=stations, line_load_n_m=[1.0] * 4)
        self.assertAlmostEqual(result.total_normal_force_n, 2.0, places=12)
        self.assertAlmostEqual(result.bending_moment_n_m[-1], 2.0, places=12)
        self.assertAlmostEqual(result.center_of_pressure_m, 1.0, places=12)


class DatabaseTests(unittest.TestCase):
    def make_db(self) -> AeroDatabase:
        mach = (0.5, 1.0, 2.0)
        alpha = (-4.0, 0.0, 4.0)
        beta = (-4.0, 4.0)
        values = []
        for m in mach:
            for a in alpha:
                for b in beta:
                    values.append(0.1 * m + 0.05 * a + 0.01 * b)  # multilinear -> exact
        return AeroDatabase(
            axes=("mach", "alpha_deg", "beta_deg"),
            breakpoints=(mach, alpha, beta),
            coefficients={"CN": tuple(values), "CA": tuple(0.3 for _ in values)},
            range_uncertainty={"CN": RangeUncertainty("mach", ((0.5, 1.3, 0.02), (1.6, 2.0, 0.03)))},
        )

    def test_multilinear_exact_and_refusals(self) -> None:
        db = self.make_db()
        query = db.query({"mach": 0.8, "alpha_deg": 1.5, "beta_deg": -1.0})
        self.assertAlmostEqual(query.nominal["CN"], 0.1 * 0.8 + 0.05 * 1.5 - 0.01, places=12)
        self.assertEqual(query.uncertainty["CN"], 0.02)
        self.assertAlmostEqual(query.dispersed({"CN": -1.0})["CN"], query.lower["CN"], places=15)
        for bad in ({"mach": 2.5, "alpha_deg": 0.0, "beta_deg": 0.0},
                    {"mach": 1.4, "alpha_deg": 0.0, "beta_deg": 0.0}):
            with self.assertRaises(DomainError):
                db.query(bad)
        with self.assertRaises(InputError):
            db.query({"mach": 1.0, "alpha_deg": 0.0})
        with self.assertRaises(DomainError):
            query.dispersed({"CN": 1.5})
        with self.assertRaises(DomainError):
            query.dispersed({"CA": 0.5})

    def test_uncertainty_arithmetic(self) -> None:
        ref = DATA["sls_uncertainty"]
        self.assertEqual({str(k): v for k, v in RANGE_METHOD_D2.items()}, ref["d2"]["values"])
        sigma = range_method_sigma([0.01, 0.03], group_size=2)
        self.assertAlmostEqual(sigma, 0.02 / 1.128, places=15)
        self.assertAlmostEqual(repeatability_uncertainty(sigma, combined_terms=2), 3.0 * math.sqrt(2.0) * sigma, places=15)
        self.assertEqual(net_modeling_uncertainty(0.01, 0.02), 0.0)
        self.assertAlmostEqual(net_modeling_uncertainty(0.05, 0.03), 0.04, places=15)
        self.assertEqual(coverage_bound(list(range(1, 1001)), coverage=0.997), 997)
        u_bal = 0.01
        build = combine_uncertainty_terms(
            [UncertaintyTerm("dbm_net", 0.0), UncertaintyTerm("tst_net", 0.02), UncertaintyTerm("rpt", 0.01),
             UncertaintyTerm("bal", u_bal), UncertaintyTerm("g2f", 0.04), UncertaintyTerm("wt_corr", u_bal)],
            uncertainty_factor=ref["final_bound"]["uncertainty_factor"],
        )
        self.assertAlmostEqual(build.total, 1.05 * math.sqrt(0.0004 + 0.0001 + 0.0001 + 0.0016 + 0.0001), places=15)
        self.assertAlmostEqual(sum(build.contributions.values()), 1.0, places=12)
        doubled = combine_uncertainty_terms([UncertaintyTerm("cal", 0.01, weight=2.0), UncertaintyTerm("rep", 0.01)])
        self.assertAlmostEqual(doubled.total, math.sqrt(3.0) * 0.01, places=15)
        with self.assertRaises(InputError):
            combine_uncertainty_terms([UncertaintyTerm("a", 0.1)], uncertainty_factor=0.9)

    def test_hemsch_fits(self) -> None:
        alpha = [-4.0, -2.0, 0.0, 2.0, 4.0]
        cn = [0.01 + 0.04 * a for a in alpha]
        line = fit_linear(alpha, cn)
        self.assertAlmostEqual(line.slope, 0.04, places=12)
        self.assertAlmostEqual(line.intercept, 0.01, places=12)
        cm = [0.002 + (3.5 - 5.0) * c for c in cn]
        fit = fit_moment_center_of_pressure(cm, cn, moment_reference=5.0)
        self.assertAlmostEqual(fit.center_of_pressure, 3.5, places=10)


if __name__ == "__main__":
    unittest.main()
