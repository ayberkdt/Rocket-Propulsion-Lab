"""Inviscid incompressible and lifting aerodynamics, checked against primary sources.

Every reference value is read from ``tests/reference_data/aero_subsonic.json``,
which records the source and page of each number:

* J. D. Anderson Jr., Fundamentals of Aerodynamics, 6th ed., 2017 - worked
  Examples 3.16-3.19, 4.5-4.7, 5.1-5.2, Eqs. (3.101), (4.33)-(4.66),
  (5.43), (5.69)-(5.70), Fig. 5.20 (read-off with stated tolerance);
* NACA Reports 824 (ordinate tables, mean-line theory examples, measured
  trends), 460 (four-digit equations), 537 (Table I), 610 (230 camber), 411
  (Joukowski map, exact elliptic and Joukowski-section flows);
* DLMF Table 3.5.1 (Gauss-Legendre nodes).

Exactness tests: the source-panel cylinder, the elliptic lift distribution,
the conformal-map flows. Convergence tests: panel refinement on exact
elliptic and Joukowski-section flows. Limits: flat plate, symmetric section,
elliptic wing, zero circulation.
"""

from __future__ import annotations

import json
import unittest
from math import atan2, atanh, cos, cosh, degrees, hypot, pi, radians, sin, sinh, tanh
from pathlib import Path

from rocket_propulsion.aerodynamics._subsonic_numerics import gauss_legendre, integrate, solve_dense
from rocket_propulsion.aerodynamics.airfoil_geometry import naca_section
from rocket_propulsion.aerodynamics.lifting_line import (
    custom_wing,
    elliptic_wing,
    finite_wing_lift_slope,
    induced_drag_coefficient,
    lift_slope,
    solve_lifting_line,
    trapezoidal_wing,
)
from rocket_propulsion.aerodynamics.panel_methods import (
    influence_matrices,
    panel_geometry,
    solve_lifting_panels,
    solve_source_panels,
)
from rocket_propulsion.aerodynamics.potential_flow import (
    CircularCylinderFlow,
    DoubletFlow,
    JoukowskiAirfoil,
    SourceFlow,
    UniformFlow,
    VortexFlow,
    ellipse_surface_speed_ratio,
    kutta_joukowski_lift_per_span,
    pressure_coefficient,
    pressure_coefficient_from_pressure,
    superpose_stream_function,
    superpose_velocity,
)
from rocket_propulsion.aerodynamics.thin_airfoil import thin_airfoil, thin_airfoil_naca
from rocket_propulsion.core.errors import ConvergenceError, DomainError, InputError

REFERENCE = json.loads(
    (Path(__file__).parent / "reference_data" / "aero_subsonic.json").read_text(encoding="utf-8")
)


def _octagon() -> list[tuple[float, float]]:
    """Anderson Fig. 3.41: unit circle, panel 1 vertical facing the stream, clockwise."""

    nodes = [(cos(radians(202.5 - 45.0 * k)), sin(radians(202.5 - 45.0 * k))) for k in range(8)]
    nodes.append(nodes[0])
    return nodes


class QuadratureTests(unittest.TestCase):
    def test_gauss_legendre_matches_dlmf_table_3_5_1(self) -> None:
        ref = REFERENCE["quadrature"]["dlmf_table_3_5_1"]
        nodes, weights = gauss_legendre(5)
        pairs = sorted(zip(nodes, weights, strict=True))
        for (node, weight), (ref_node, ref_weight) in zip(
            pairs, [(-ref["nodes"][2], ref["weights"][2]), (-ref["nodes"][1], ref["weights"][1]),
                    (ref["nodes"][0], ref["weights"][0]), (ref["nodes"][1], ref["weights"][1]),
                    (ref["nodes"][2], ref["weights"][2])], strict=True):
            self.assertAlmostEqual(node, ref_node, delta=ref["tolerance"])
            self.assertAlmostEqual(weight, ref_weight, delta=ref["tolerance"])

    def test_composite_rule_and_linear_solver(self) -> None:
        self.assertAlmostEqual(integrate(sin, [0.0, pi]), 2.0, places=13)
        self.assertEqual(solve_dense([[2.0, 1.0], [1.0, 3.0]], [3.0, 5.0]), [0.8, 1.4])
        with self.assertRaises(ConvergenceError):
            solve_dense([[1.0, 2.0], [2.0, 4.0]], [1.0, 2.0])


class PotentialFlowTests(unittest.TestCase):
    def test_uniform_plus_doublet_is_the_nonlifting_cylinder(self) -> None:
        speed, radius = 3.0, 0.7
        flows = [UniformFlow(speed), DoubletFlow(2.0 * pi * speed * radius**2)]
        cylinder = CircularCylinderFlow(speed, radius)
        for r, theta in ((0.9, 0.3), (1.5, 2.0), (4.0, -1.2)):
            x, y = r * cos(theta), r * sin(theta)
            velocity = superpose_velocity(flows, x, y)
            v_r, v_t = cylinder.velocity_polar(r, theta)
            self.assertAlmostEqual(velocity.u, v_r * cos(theta) - v_t * sin(theta), places=12)
            self.assertAlmostEqual(velocity.v, v_r * sin(theta) + v_t * cos(theta), places=12)
            # psi = 0 on r = R (Eq. 3.97)
        self.assertAlmostEqual(superpose_stream_function(flows, radius * cos(1.1), radius * sin(1.1)), 0.0,
                               places=12)

    def test_vortex_source_signs_follow_anderson(self) -> None:
        vortex = VortexFlow(2.0 * pi)
        top = vortex.velocity(0.0, 1.0)  # clockwise positive: flow to +x above the vortex
        self.assertAlmostEqual(top.u, 1.0, places=14)
        source = SourceFlow(2.0 * pi)
        self.assertAlmostEqual(source.velocity(2.0, 0.0).u, 0.5, places=14)
        self.assertAlmostEqual(vortex.stream_function(1.0, 0.0), 0.0, places=14)
        with self.assertRaises(DomainError):
            source.velocity(0.0, 0.0)

    def test_nonlifting_cylinder_pressure_limits(self) -> None:
        ref = REFERENCE["potential_flow"]["nonlifting_cylinder"]
        cylinder = CircularCylinderFlow(10.0, 1.0)
        self.assertAlmostEqual(cylinder.surface_pressure_coefficient(0.0), ref["stagnation_cp"], places=14)
        self.assertAlmostEqual(cylinder.surface_pressure_coefficient(0.5 * pi), ref["minimum_cp"], places=14)
        self.assertEqual(cylinder.lift_coefficient, 0.0)

    def test_example_3_16_and_3_17(self) -> None:
        ref16 = REFERENCE["potential_flow"]["example_3_16"]
        ref17 = REFERENCE["potential_flow"]["example_3_17"]
        cylinder = CircularCylinderFlow(1.0, 1.0, ref16["lift_coefficient"])
        self.assertAlmostEqual(cylinder.lift_coefficient, ref16["lift_coefficient"], places=14)
        peak = cylinder.surface_pressure_coefficient(0.5 * pi)
        self.assertAlmostEqual(peak, ref16["peak_pressure_coefficient"], delta=ref16["printed_digits_tolerance"])
        stagnation = [degrees(theta) for _, theta in cylinder.stagnation_points()]
        for value, expected in zip(stagnation, ref17["stagnation_angles_deg"], strict=True):
            self.assertAlmostEqual(value, expected, delta=ref17["angle_tolerance_deg"])
        angles = [degrees(theta) for theta in cylinder.freestream_pressure_angles()]
        for value, expected in zip(angles, ref17["freestream_pressure_angles_deg"], strict=True):
            self.assertAlmostEqual(value, expected, delta=ref17["angle_tolerance_deg"])
        self.assertAlmostEqual(cylinder.surface_pressure_coefficient(1.5 * pi), ref17["bottom_pressure_coefficient"],
                               delta=ref17["bottom_tolerance"])

    def test_example_3_18_lift_per_span(self) -> None:
        ref = REFERENCE["potential_flow"]["example_3_18"]
        radius = 0.5 * ref["diameter_m"]
        speed = ref["freestream_speed_m_s"]
        # V_theta(90 deg) = -2V - Gamma/(2 pi R) = -V_max  (Eq. 3.125)
        circulation = -2.0 * pi * radius * (-ref["maximum_surface_speed_m_s"] + 2.0 * speed)
        cylinder = CircularCylinderFlow(speed, radius, circulation)
        self.assertAlmostEqual(-cylinder.surface_velocity(0.5 * pi), ref["maximum_surface_speed_m_s"], places=12)
        self.assertAlmostEqual(circulation, ref["circulation_m2_s"], delta=ref["relative_tolerance"] * 40.0)
        lift = kutta_joukowski_lift_per_span(ref["density_kg_m3"], speed, circulation)
        self.assertAlmostEqual(lift, ref["lift_per_span_n_m"], delta=ref["relative_tolerance"] * 900.0)

    def test_stagnation_regimes_and_pressure_definition(self) -> None:
        one = CircularCylinderFlow(1.0, 1.0, 4.0 * pi)
        self.assertEqual(len(one.stagnation_points()), 1)
        self.assertAlmostEqual(one.stagnation_points()[0][1], 1.5 * pi, places=14)
        off = CircularCylinderFlow(1.0, 1.0, 8.0 * pi).stagnation_points()
        self.assertGreater(off[0][0], 1.0)  # Eq. (3.124) outer root lies off the body
        self.assertAlmostEqual(pressure_coefficient(1.5, 1.0), -1.25, places=14)
        self.assertAlmostEqual(pressure_coefficient_from_pressure(101_000.0, 100_000.0, 1.2, 10.0),
                               1000.0 / 60.0, places=10)
        with self.assertRaises(DomainError):
            CircularCylinderFlow(1.0, -1.0)

    def test_joukowski_flat_plate_limit_is_two_pi_sin_alpha(self) -> None:
        # Circle of radius a centred at the origin maps to a flat plate of chord 4a.
        alpha = radians(6.0)
        plate = JoukowskiAirfoil(1.0, 0.0, 0.0, alpha)
        self.assertAlmostEqual(plate.chord, 4.0, places=9)
        self.assertAlmostEqual(plate.lift_coefficient, 2.0 * pi * sin(alpha), places=9)

    def test_joukowski_ellipse_matches_tr411_closed_form(self) -> None:
        a, psi, alpha = 1.0, 0.3, radians(7.0)
        section = JoukowskiAirfoil(a, 0.0, 0.0, alpha, radius=a * (cosh(psi) + sinh(psi)))
        self.assertAlmostEqual(section.chord, 4.0 * a * cosh(psi), places=9)
        self.assertAlmostEqual(section.lift_coefficient, 2.0 * pi * (1.0 + tanh(psi)) * sin(alpha), places=9)
        for theta in (0.4, 1.3, 2.5, 3.9, 5.5):
            ratio = ellipse_surface_speed_ratio(psi, alpha, theta)
            self.assertAlmostEqual(section.surface_speed(theta), ratio, places=10)
        self.assertAlmostEqual(ellipse_surface_speed_ratio(psi, alpha, 0.0), 0.0, places=14)

    def test_cusped_joukowski_kutta_circulation_is_tr411_eq_vii(self) -> None:
        a, cx, cy, alpha = 1.0, -0.1, 0.1, radians(5.0)
        section = JoukowskiAirfoil(a, cx, cy, alpha)
        radius = hypot(a - cx, cy)
        beta = atan2(cy, a - cx)
        self.assertTrue(section.is_cusped)
        self.assertAlmostEqual(section.circulation, 4.0 * pi * radius * sin(alpha + beta), places=12)
        self.assertAlmostEqual(section.surface_point(section.trailing_edge_angle)[0], 2.0 * a, places=12)
        # Kutta condition: speed vanishes approaching the trailing edge on both sides (finite, small)
        near = [section.surface_speed(section.trailing_edge_angle + s) for s in (1e-3, -1e-3)]
        self.assertLess(abs(near[0] - near[1]), 0.05)
        with self.assertRaises(DomainError):
            section.surface_speed(section.trailing_edge_angle)


class AirfoilGeometryTests(unittest.TestCase):
    def _check_table(self, designation: str, ref: dict) -> None:
        section = naca_section(designation)
        tol = ref["print_tolerance_percent"]
        for station, (xu, yu), (xl, yl) in zip(ref["nominal_stations_percent"], ref["upper"], ref["lower"],
                                               strict=True):
            ux, uy = section.upper_surface(station / 100.0)
            lx, ly = section.lower_surface(station / 100.0)
            if station < 100.0:  # TR-824 prints the trailing-edge station as 100.000 on both surfaces
                self.assertAlmostEqual(100.0 * ux, xu, delta=tol)
                self.assertAlmostEqual(100.0 * lx, xl, delta=tol)
            self.assertAlmostEqual(100.0 * uy, yu, delta=tol)
            self.assertAlmostEqual(100.0 * ly, yl, delta=tol)
        # Printed to two decimals (1.10 t^2 = 6.336 for t = 0.24 is printed 6.33).
        self.assertAlmostEqual(100.0 * section.leading_edge_radius, ref["leading_edge_radius_percent"], delta=0.01)
        self.assertAlmostEqual(section.mean_line(0.0)[1], ref["slope_of_radius_through_le"], delta=0.0006)

    def test_naca_1412_matches_tr824_table(self) -> None:
        self._check_table("1412", REFERENCE["airfoil_geometry"]["naca_1412_ordinates_percent"])

    def test_naca_23024_matches_tr824_table(self) -> None:
        self._check_table("23024", REFERENCE["airfoil_geometry"]["naca_23024_ordinates_percent"])

    def test_naca_0006_and_0009_thickness(self) -> None:
        ref = REFERENCE["airfoil_geometry"]["naca_0006_ordinates_percent"]
        section = naca_section("NACA 0006")
        for station, ordinate in ref["points"]:
            self.assertAlmostEqual(100.0 * section.thickness(station / 100.0), ordinate,
                                   delta=ref["print_tolerance_percent"])
        self.assertAlmostEqual(100.0 * section.leading_edge_radius, ref["leading_edge_radius_percent"],
                               delta=ref["print_tolerance_percent"])
        ref9 = REFERENCE["airfoil_geometry"]["naca_0009_leading_edge_radius_percent"]
        self.assertAlmostEqual(100.0 * naca_section("0009").leading_edge_radius, ref9["value"], delta=ref9["tolerance"])

    def test_five_digit_mean_lines(self) -> None:
        anderson = REFERENCE["thin_airfoil"]["example_4_6_naca_23012"]["mean_line_constants"]
        section = naca_section("23012")
        self.assertAlmostEqual(section.k1 / 6.0, anderson["k1_over_6"], delta=0.00005)
        self.assertAlmostEqual(3.0 * section.junction, anderson["three_m"], places=12)
        self.assertAlmostEqual(section.junction**2 * (3.0 - section.junction), anderson["m_squared_3_minus_m"],
                               delta=0.00005)
        self.assertAlmostEqual(-section.mean_line(0.5)[1], anderson["k1_m_cubed_over_6"], delta=0.000005)
        ref610 = REFERENCE["airfoil_geometry"]["tr_610_230_camber"]
        self.assertAlmostEqual(100.0 * section.max_camber, ref610["camber_percent"],
                               delta=ref610["print_tolerance_percent"])
        table = REFERENCE["airfoil_geometry"]["tr_537_table_i"]
        for line, column in table["reflexed_k2_over_k1_column"].items():
            ratio = naca_section(line + "12").k2_over_k1
            self.assertAlmostEqual(1.0e4 * ratio, column, delta=table["relative_tolerance"] * column + 0.01)
        for designation in ("22112", "23112", "24112", "25112", "21012", "25012", "43012"):
            reflexed = naca_section(designation)
            samples = [i / 4000.0 for i in range(4001)]
            peak = max(samples, key=lambda x: reflexed.mean_line(x)[0])
            self.assertAlmostEqual(peak, reflexed.camber_position, delta=0.0005)  # max camber at p
            self.assertAlmostEqual(reflexed.mean_line(1.0)[0], 0.0, places=12)

    def test_designation_validation(self) -> None:
        for bad in ("12", "0000", "2012", "21112", "33112", "ab12", "123456"):
            with self.assertRaises(InputError):
                naca_section(bad)
        with self.assertRaises(DomainError):
            naca_section("2412").thickness(1.5)

    def test_closed_panel_nodes(self) -> None:
        section = naca_section("2412")
        closed = section.panel_nodes(20)
        self.assertEqual(closed[0], closed[-1])
        opened = section.panel_nodes(20, close_trailing_edge=False)
        gap = hypot(opened[-1][0] - opened[0][0], opened[-1][1] - opened[0][1])
        self.assertAlmostEqual(gap, 2.0 * section.thickness(1.0), places=12)
        self.assertLessEqual(max(abs(p[1] - q[1]) for p, q in zip(closed, opened, strict=True)),
                             section.thickness(1.0) + 1e-15)


class ThinAirfoilTests(unittest.TestCase):
    def test_example_4_5_flat_plate(self) -> None:
        ref = REFERENCE["thin_airfoil"]["example_4_5_flat_plate"]
        result = thin_airfoil(lambda x: 0.0, radians(ref["alpha_deg"]))
        self.assertAlmostEqual(result.lift_coefficient, ref["lift_coefficient"], delta=ref["tolerance"])
        self.assertAlmostEqual(result.moment_coefficient_leading_edge, ref["moment_coefficient_leading_edge"],
                               delta=ref["tolerance"])
        self.assertAlmostEqual(result.moment_coefficient_quarter_chord, ref["moment_coefficient_quarter_chord"],
                               places=14)
        self.assertAlmostEqual(result.center_of_pressure, 0.25, places=14)
        self.assertEqual(thin_airfoil(lambda x: 0.0, 0.0).center_of_pressure, None)

    def test_example_4_6_naca_23012(self) -> None:
        ref = REFERENCE["thin_airfoil"]["example_4_6_naca_23012"]
        tol = ref["tolerances"]
        constants = ref["mean_line_constants"]

        def slope(x: float) -> float:
            if x < constants["junction"]:
                return constants["k1_over_6"] * (3.0 * x * x - 2.0 * constants["three_m"] * x
                                                 + constants["m_squared_3_minus_m"])
            return -constants["k1_m_cubed_over_6"]

        for result in (
            thin_airfoil(slope, radians(ref["alpha_deg"]), slope_breakpoints=(constants["junction"],)),
            thin_airfoil_naca(naca_section("23012"), radians(ref["alpha_deg"])),
        ):
            self.assertAlmostEqual(degrees(result.zero_lift_angle_rad), ref["zero_lift_angle_deg"],
                                   delta=tol["zero_lift_angle_deg"])
            self.assertAlmostEqual(result.lift_coefficient, ref["lift_coefficient"], delta=tol["lift_coefficient"])
            self.assertAlmostEqual(result.fourier_coefficients[0], ref["a1"], delta=tol["a1"])
            self.assertAlmostEqual(result.fourier_coefficients[1], ref["a2"], delta=tol["a2"])
            self.assertAlmostEqual(result.moment_coefficient_quarter_chord, ref["moment_coefficient_quarter_chord"],
                                   delta=tol["moment_coefficient_quarter_chord"])
            self.assertAlmostEqual(result.center_of_pressure, ref["center_of_pressure_x_over_c"],
                                   delta=tol["center_of_pressure_x_over_c"])
        # Theory vs the measured values Anderson tabulates (p. 364): within 0.03 in c_l, 0.05 deg.
        experiment = ref["experiment"]
        result = thin_airfoil_naca(naca_section("23012"), radians(ref["alpha_deg"]))
        self.assertAlmostEqual(degrees(result.zero_lift_angle_rad), experiment["zero_lift_angle_deg"], delta=0.05)
        self.assertAlmostEqual(result.lift_coefficient, experiment["lift_coefficient_at_4_deg"], delta=0.03)

    def test_relations_hold_for_any_camber(self) -> None:
        section = naca_section("4412")
        low = thin_airfoil_naca(section, radians(-2.0))
        high = thin_airfoil_naca(section, radians(6.0))
        self.assertAlmostEqual((high.lift_coefficient - low.lift_coefficient) / radians(8.0), 2.0 * pi, places=10)
        self.assertAlmostEqual(high.moment_coefficient_quarter_chord, low.moment_coefficient_quarter_chord,
                               places=13)  # aerodynamic centre at c/4 (Eq. 4.64)
        self.assertAlmostEqual(high.lift_coefficient, 2.0 * pi * (radians(6.0) - high.zero_lift_angle_rad),
                               places=10)  # Eq. (4.60)
        self.assertAlmostEqual(high.moment_coefficient_leading_edge,
                               -high.lift_coefficient / 4.0 + high.moment_coefficient_quarter_chord, places=12)
        self.assertEqual(high.warnings, ())
        self.assertIn("thin_airfoil_thick_section",
                      [w.code for w in thin_airfoil_naca(naca_section("4415"), 0.0).warnings])
        self.assertIn("thin_airfoil_large_incidence",
                      [w.code for w in thin_airfoil_naca(section, radians(12.0)).warnings])

    def test_naca_824_mean_line_examples(self) -> None:
        ref = REFERENCE["thin_airfoil"]["naca_824_mean_line_theory"]
        tol = ref["tolerances"]
        line64 = thin_airfoil_naca(naca_section("6412"), 0.0)
        self.assertAlmostEqual(degrees(line64.ideal_angle_of_attack_rad), ref["naca_64_mean_line"]["ideal_angle_deg"],
                               delta=tol["ideal_angle_deg"])
        self.assertAlmostEqual(line64.ideal_lift_coefficient, ref["naca_64_mean_line"]["ideal_lift_coefficient"],
                               delta=tol["ideal_lift_coefficient_relative"] * 0.76)
        self.assertAlmostEqual(line64.moment_coefficient_quarter_chord,
                               ref["naca_64_mean_line"]["moment_coefficient_quarter_chord"],
                               delta=tol["moment_coefficient_relative"] * 0.157)
        self.assertAlmostEqual(thin_airfoil_naca(naca_section("4415"), 0.0).moment_coefficient_quarter_chord,
                               ref["naca_4415_moment_coefficient_quarter_chord"],
                               delta=tol["moment_coefficient_relative"] * 0.105)
        line24 = thin_airfoil_naca(naca_section("2415"), 0.0)
        given = ref["naca_2415"]
        tr824_zero_lift = given["ideal_angle_deg"] - 57.3 * given["ideal_lift_coefficient"] / (2.0 * pi)
        self.assertAlmostEqual(degrees(line24.zero_lift_angle_rad), tr824_zero_lift, delta=tol["zero_lift_angle_deg"])
        # TR-824 Eq. (16) is the same thin-airfoil relation: alpha_L0 = alpha_i - c_li / (2 pi)
        self.assertAlmostEqual(line24.zero_lift_angle_rad,
                               line24.ideal_angle_of_attack_rad - line24.ideal_lift_coefficient / (2.0 * pi),
                               places=12)

    def test_naca_2412_against_measured_zero_lift_angle(self) -> None:
        ref = REFERENCE["thin_airfoil"]["naca_2412_experiment"]
        result = thin_airfoil_naca(naca_section("2412"), radians(4.0))
        self.assertAlmostEqual(degrees(result.zero_lift_angle_rad), ref["zero_lift_angle_deg"],
                               delta=ref["read_tolerance_deg"])
        # Thin-airfoil slope 2 pi /rad = 0.1097 /deg vs TR-824 "about 0.110 /deg" for 6-10 % sections
        slope = REFERENCE["thin_airfoil"]["naca_824_lift_slope"]["lift_slope_per_deg"]
        self.assertAlmostEqual(radians(result.lift_slope_per_rad), slope, delta=0.0005)


class PanelMethodTests(unittest.TestCase):
    def test_example_3_19_influence_coefficients_and_strengths(self) -> None:
        ref = REFERENCE["panel_methods"]["example_3_19"]
        normal, _ = influence_matrices(panel_geometry(_octagon()))
        for j, value in ref["I_4j"].items():
            self.assertAlmostEqual(normal[3][int(j) - 1], value, delta=ref["printed_digits_tolerance"])
        solution = solve_source_panels(_octagon())
        for value, expected in zip(solution.source_strengths, ref["lambda_over_2_pi_v"], strict=True):
            self.assertAlmostEqual(value / (2.0 * pi), expected, delta=ref["printed_digits_tolerance"])
        self.assertAlmostEqual(solution.source_closure, 0.0, places=12)  # Eq. (3.157)

    def test_closed_form_integrals_match_quadrature(self) -> None:
        geometry = panel_geometry(naca_section("4412").panel_nodes(12))
        normal, tangential = influence_matrices(geometry)
        nodes, weights = gauss_legendre(40)
        for i, j in ((3, 15), (10, 11), (0, 23), (12, 2)):
            xi, yi = geometry.control_points[i]
            xj, yj = geometry.nodes[j]
            phi_i, phi_j, length = geometry.angles[i], geometry.angles[j], geometry.lengths[j]
            beta = phi_i + 0.5 * pi
            n_sum = t_sum = 0.0
            for node, weight in zip(nodes, weights, strict=True):
                s = 0.5 * length * (node + 1.0)
                dx = xi - (xj + s * cos(phi_j))
                dy = yi - (yj + s * sin(phi_j))
                r2 = dx * dx + dy * dy
                n_sum += weight * (dx * cos(beta) + dy * sin(beta)) / r2   # Eq. (3.159)
                t_sum += weight * (dx * cos(phi_i) + dy * sin(phi_i)) / r2
            self.assertAlmostEqual(normal[i][j], 0.5 * length * n_sum, delta=1e-6)
            self.assertAlmostEqual(tangential[i][j], 0.5 * length * t_sum, delta=1e-6)

    def test_source_panels_reproduce_exact_cylinder_pressure(self) -> None:
        for count in (8, 32):
            step = 2.0 * pi / count
            nodes = [(cos(pi + 0.5 * step - k * step), sin(pi + 0.5 * step - k * step)) for k in range(count)]
            nodes.append(nodes[0])
            solution = solve_source_panels(nodes, freestream_speed=7.0)
            for cp, (x, y) in zip(solution.pressure_coefficients, solution.geometry.control_points, strict=True):
                self.assertAlmostEqual(cp, 1.0 - 4.0 * sin(atan2(y, x)) ** 2, places=11)  # Eq. (3.101)

    def test_lifting_panels_converge_to_exact_ellipse(self) -> None:
        thickness, alpha = 0.12, radians(5.0)
        psi = atanh(thickness)
        exact = 2.0 * pi * (1.0 + thickness) * sin(alpha)
        errors = []
        for count in (40, 80, 160):
            nodes = [(2.0 * cosh(psi) * cos(-2.0 * pi * k / count), 2.0 * sinh(psi) * sin(-2.0 * pi * k / count))
                     for k in range(count)]
            nodes.append(nodes[0])
            errors.append(abs(solve_lifting_panels(nodes, alpha).lift_coefficient - exact))
        self.assertGreater(errors[0] / errors[1], 1.8)
        self.assertGreater(errors[1] / errors[2], 1.8)
        self.assertLess(errors[2] / exact, 0.015)

    def test_lifting_panels_converge_to_exact_joukowski_section(self) -> None:
        alpha = radians(5.0)
        section = JoukowskiAirfoil(1.0, -0.08, 0.08, alpha, radius=1.3 * hypot(1.08, 0.08))
        exact = section.lift_coefficient
        results = [solve_lifting_panels(section.surface_nodes(count), alpha) for count in (40, 80, 160)]
        errors = [abs(r.lift_coefficient - exact) for r in results]
        self.assertGreater(errors[0] / errors[1], 1.8)
        self.assertGreater(errors[1] / errors[2], 1.8)
        self.assertLess(errors[2] / exact, 0.01)
        finest = results[-1]
        self.assertAlmostEqual(finest.pressure_lift_coefficient, exact, delta=0.01 * exact)
        self.assertLess(abs(finest.pressure_drag_coefficient), 0.005)
        # Surface pressure away from the stagnation regions within 0.02 of the exact value
        count = 160
        start = section.trailing_edge_angle
        for k in range(10, count - 10, 15):
            theta = start - 2.0 * pi * (k + 0.5) / count
            self.assertAlmostEqual(finest.pressure_coefficients[k], section.surface_pressure_coefficient(theta),
                                   delta=0.02)

    def test_cusped_trailing_edge_is_rejected(self) -> None:
        cusped = JoukowskiAirfoil(1.0, -0.1, 0.1, radians(5.0))
        with self.assertRaises(DomainError):
            solve_lifting_panels(cusped.surface_nodes(80), radians(5.0))

    def test_contour_validation(self) -> None:
        nodes = _octagon()
        with self.assertRaises(DomainError):
            solve_source_panels(list(reversed(nodes)))  # counterclockwise
        with self.assertRaises(DomainError):
            solve_source_panels(nodes[:-1])  # open
        with self.assertRaises(DomainError):
            solve_lifting_panels(naca_section("0012").panel_nodes(10, close_trailing_edge=False), 0.0)

    def test_symmetric_section_and_naca_2412_against_experiment(self) -> None:
        symmetric = solve_lifting_panels(naca_section("0012").panel_nodes(40), 0.0)
        self.assertAlmostEqual(symmetric.lift_coefficient, 0.0, places=8)
        ref = REFERENCE["thin_airfoil"]["naca_2412_experiment"]
        band = REFERENCE["panel_methods"]["naca_2412_experiment_comparison"]
        section = naca_section("2412")
        zero = solve_lifting_panels(section.panel_nodes(40), 0.0)
        four = solve_lifting_panels(section.panel_nodes(40), radians(4.0))
        slope = (four.lift_coefficient - zero.lift_coefficient) / 4.0
        zero_lift_deg = -zero.lift_coefficient / slope
        self.assertAlmostEqual(zero_lift_deg, ref["zero_lift_angle_deg"], delta=band["allowed_zero_lift_angle_difference_deg"])
        ratio = four.lift_coefficient / ref["lift_coefficient_at_4_deg"]
        self.assertGreater(ratio, band["documented_lift_overprediction_band"][0])
        self.assertLess(ratio, band["documented_lift_overprediction_band"][1])
        self.assertAlmostEqual(four.pressure_lift_coefficient, four.lift_coefficient, delta=0.03)


class LiftingLineTests(unittest.TestCase):
    def test_elliptic_wing_is_exact(self) -> None:
        wing = elliptic_wing(12.0, 1.5)
        alpha = radians(4.0)
        result = solve_lifting_line(wing, alpha, terms=15)
        aspect_ratio = wing.aspect_ratio
        expected = finite_wing_lift_slope(2.0 * pi, aspect_ratio) * alpha  # Eq. (5.69)
        self.assertAlmostEqual(result.lift_coefficient, expected, places=12)
        self.assertAlmostEqual(result.induced_drag_coefficient, induced_drag_coefficient(expected, aspect_ratio),
                               places=12)  # Eq. (5.43)
        self.assertAlmostEqual(result.induced_drag_factor, 0.0, places=14)
        self.assertAlmostEqual(result.span_efficiency, 1.0, places=14)
        induced = result.lift_coefficient / (pi * aspect_ratio)  # Eq. (5.42), constant along the span
        for value in result.induced_angles_rad:
            self.assertAlmostEqual(value, induced, places=12)
        for value in result.section_lift_coefficients:
            self.assertAlmostEqual(value, result.lift_coefficient, places=10)

    def test_examples_5_1_and_5_2(self) -> None:
        ref = REFERENCE["lifting_line"]["example_5_1"]
        rel = ref["relative_tolerance"]
        a = finite_wing_lift_slope(2.0 * pi, ref["aspect_ratio"], ref["tau"])
        self.assertAlmostEqual(a, ref["lift_slope_per_rad"], delta=rel * ref["lift_slope_per_rad"])
        lift = radians(a) * ref["alpha_deg"]
        self.assertAlmostEqual(lift, ref["lift_coefficient"], delta=rel * ref["lift_coefficient"])
        self.assertAlmostEqual(induced_drag_coefficient(ref["lift_coefficient"], ref["aspect_ratio"], ref["delta"]),
                               ref["induced_drag_coefficient"], delta=rel * ref["induced_drag_coefficient"])
        ref2 = REFERENCE["lifting_line"]["example_5_2"]
        rel2 = ref2["relative_tolerance"]
        lift6 = (6.0 * pi * ref2["induced_drag_ar6"] / 1.055) ** 0.5
        self.assertAlmostEqual(lift6, ref2["lift_coefficient_ar6"], delta=rel2 * 0.423)
        self.assertAlmostEqual(finite_wing_lift_slope(ref2["section_lift_slope_per_rad"], 6.0, 0.055),
                               ref2["lift_slope_ar6_per_rad"], delta=rel2 * 4.485)
        a10 = finite_wing_lift_slope(ref2["section_lift_slope_per_rad"], 10.0, 0.105)
        self.assertAlmostEqual(a10, ref2["lift_slope_ar10_per_rad"], delta=rel2 * 4.95)
        lift10 = radians(a10) * 5.4
        self.assertAlmostEqual(lift10, ref2["lift_coefficient_ar10"], delta=rel2 * 0.464)
        self.assertAlmostEqual(induced_drag_coefficient(lift10, 10.0, 0.105), ref2["induced_drag_ar10"], delta=0.0001)

    def test_series_convergence_and_twist(self) -> None:
        wing = trapezoidal_wing(10.0, 1.6, 0.4)
        coarse, middle, fine = (lift_slope(wing, terms=n) for n in (20, 40, 80))
        first = abs(middle.lift_slope_per_rad - coarse.lift_slope_per_rad)
        second = abs(fine.lift_slope_per_rad - middle.lift_slope_per_rad)
        self.assertGreater(first / second, 3.0)  # converging (the root chord kink limits the rate)
        self.assertLess(second / fine.lift_slope_per_rad, 2e-4)
        self.assertAlmostEqual(middle.induced_drag_factor, fine.induced_drag_factor, delta=1e-4)
        self.assertAlmostEqual(fine.span_efficiency, 1.0 / (1.0 + fine.induced_drag_factor), places=14)
        # Washout: zero net lift leaves a nonzero (basic) loading and positive induced drag.
        twisted = trapezoidal_wing(10.0, 1.6, 0.4, tip_twist_rad=radians(-4.0))
        flat = solve_lifting_line(twisted, 0.0)
        slope = solve_lifting_line(twisted, 1.0).lift_coefficient - flat.lift_coefficient
        zero_lift_alpha = -flat.lift_coefficient / slope
        balanced = solve_lifting_line(twisted, zero_lift_alpha)
        self.assertAlmostEqual(balanced.lift_coefficient, 0.0, places=10)
        self.assertGreater(balanced.induced_drag_coefficient, 0.0)
        self.assertGreater(max(balanced.circulation), 0.0)
        self.assertLess(min(balanced.circulation), 0.0)
        # A custom description of the same planform gives the same answer.
        same = custom_wing(10.0, lambda y: 1.6 * (1.0 - 0.6 * abs(y) / 5.0), breakpoints=(0.0,))
        self.assertAlmostEqual(same.area, wing.area, places=12)
        self.assertAlmostEqual(lift_slope(same).lift_slope_per_rad, lift_slope(wing).lift_slope_per_rad, places=12)

    def test_figure_5_20_trend_and_documented_discrepancy(self) -> None:
        ref = REFERENCE["lifting_line"]["figure_5_20"]
        computed = {}
        for aspect_ratio, figure_value in ref["delta_at_taper_1"].items():
            ar = float(aspect_ratio)
            delta = lift_slope(trapezoidal_wing(1.0, 1.0 / ar, 1.0)).induced_drag_factor
            computed[ar] = delta
            ratio = delta / figure_value
            self.assertGreater(ratio, ref["documented_ratio_band"][0])
            self.assertLess(ratio, ref["documented_ratio_band"][1])
        ordered = [computed[ar] for ar in sorted(computed)]
        self.assertEqual(ordered, sorted(ordered))  # delta grows with AR at taper 1, as in the figure
        tapers = [0.2 + 0.02 * k for k in range(16)]
        values = [lift_slope(trapezoidal_wing(1.0, 2.0 / (6.0 * (1.0 + t)), t), terms=30).induced_drag_factor
                  for t in tapers]
        best = tapers[values.index(min(values))]
        self.assertAlmostEqual(best, ref["minimum_taper_ratio_read"], delta=ref["minimum_taper_read_tolerance"])

    def test_input_validation(self) -> None:
        with self.assertRaises(DomainError):
            trapezoidal_wing(10.0, 1.0, 1.5)
        with self.assertRaises(DomainError):
            lift_slope(trapezoidal_wing(10.0, 1.0, 0.5, tip_twist_rad=0.01))
        with self.assertRaises(DomainError):
            solve_lifting_line(elliptic_wing(10.0, 1.0), float("nan"))
        with self.assertRaises(DomainError):
            induced_drag_coefficient(0.5, 8.0, -0.1)


if __name__ == "__main__":
    unittest.main()
