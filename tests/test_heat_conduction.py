"""Conduction layer checked against a textbook, DLMF identities and its own convergence.

Every printed reference value is read from ``tests/reference_data/heat_conduction.json``,
which records the source and locator of each number:

* J. H. Lienhard IV and J. H. Lienhard V, A Heat Transfer Textbook, 6th ed. (v6.00,
  2024) - Table 5.2 (one-term coefficients, all 41 rows, three shapes), Table 5.3
  (error function), eqn. (5.33) printed values, Examples 1.4-1.6, 2.1, 2.5-2.10, 4.8-4.10,
  5.2, 5.4, 5.5, 5.7 and the printed answers of Problems 2.15 and 2.20;
* Yovanovich, Culham and Lemczyk (1988) for the annular fin, checked here against an
  independent finite-difference solution and the straight-fin limit;
* NIST DLMF identities for the Bessel functions (Wronskian 10.28.2, derivatives 10.29.3,
  10.6.3, the power series 10.2.2 as an independent representation).

Convergence evidence: the finite-volume solver is compared with the series solutions
(observed spatial order), with itself under step halving from an eigenmode start
(observed temporal order for Crank-Nicolson and backward Euler), with the semi-infinite
flux solution, and with the independent Kirchhoff steady solver for a nonlinear
multilayer wall. Energy conservation of the solver is audited in every run.
"""

from __future__ import annotations

import json
import math
import unittest
from pathlib import Path

from rocket_propulsion.core.constants import CELSIUS_OFFSET_K
from rocket_propulsion.core.errors import ConvergenceError, DomainError, InputError
from rocket_propulsion.heat_transfer import conduction as c

DATA = json.loads((Path(__file__).resolve().parent / "reference_data" / "heat_conduction.json")
                  .read_text(encoding="utf-8"))


def kelvin(celsius: float) -> float:
    return celsius + CELSIUS_OFFSET_K


def bisect(function, low: float, high: float, iterations: int = 200) -> float:
    f_low = function(low)
    for _ in range(iterations):
        middle = 0.5 * (low + high)
        f_middle = function(middle)
        if (f_middle > 0.0) == (f_low > 0.0):
            low, f_low = middle, f_middle
        else:
            high = middle
    return 0.5 * (low + high)


def observed_orders(errors: list[float]) -> list[float]:
    return [math.log2(errors[i] / errors[i + 1]) for i in range(len(errors) - 1)]


class BesselFunctionTests(unittest.TestCase):
    ARGUMENTS = (1e-3, 0.1, 0.5, 1.0, 1.9, 2.0, 2.1, 5.0, 10.0, 30.0, 49.9, 50.0, 50.1, 100.0,
                 600.0)

    def test_wronskian_dlmf_10_28_2(self) -> None:
        for x in self.ARGUMENTS:
            value = x * (c.bessel_i0e(x) * c.bessel_k1e(x) + c.bessel_i1e(x) * c.bessel_k0e(x))
            self.assertAlmostEqual(value, 1.0, delta=5e-15, msg=x)

    def test_k_series_and_integral_representations_agree(self) -> None:
        # DLMF 10.31.1-2 (series) against DLMF 10.32.9 (integral, trapezoidal rule)
        for x in (0.05, 0.5, 1.0, 1.5, 2.0):
            for order in (0, 1):
                series = c._k_series(order, x) * math.exp(x)
                integral = c._scaled_k_integral(order, x)
                self.assertAlmostEqual(series / integral, 1.0, delta=5e-15, msg=(order, x))

    def test_derivative_identities_dlmf_10_29_3_and_10_6_3(self) -> None:
        h = 1e-5
        for x in (0.3, 1.0, 3.0, 12.0):
            d_i0 = (c.bessel_i0(x + h) - c.bessel_i0(x - h)) / (2 * h)
            d_k0 = (c.bessel_k0(x + h) - c.bessel_k0(x - h)) / (2 * h)
            d_j0 = (c.bessel_j0(x + h) - c.bessel_j0(x - h)) / (2 * h)
            self.assertAlmostEqual(d_i0 / c.bessel_i1(x), 1.0, delta=1e-8)
            self.assertAlmostEqual(d_k0 / -c.bessel_k1(x), 1.0, delta=1e-8)
            self.assertAlmostEqual(d_j0, -c.bessel_j1(x), delta=1e-9)

    def test_j_matches_the_independent_power_series_dlmf_10_2_2(self) -> None:
        for x in (0.0, 0.5, 2.0, 4.0, 6.0):
            for order in (0, 1):
                terms = [(-1) ** k * (0.25 * x * x) ** k
                         / (math.factorial(k) * math.factorial(k + order)) for k in range(40)]
                series = (0.5 * x) ** order * math.fsum(terms)
                function = c.bessel_j0 if order == 0 else c.bessel_j1
                self.assertAlmostEqual(function(x), series, delta=2e-15, msg=(order, x))

    def test_first_zero_of_j0_against_table_5_2(self) -> None:
        ref = DATA["bessel"]
        zero = bisect(c.bessel_j0, 2.0, 3.0)
        self.assertAlmostEqual(zero, ref["first_zero_j0"], delta=6e-6)
        self.assertAlmostEqual(2.0 / (zero * c.bessel_j1(zero)),
                               ref["cylinder_a1_infinite_biot"], delta=6e-5)

    def test_domain_and_overflow_fail_closed(self) -> None:
        with self.assertRaises(DomainError):
            c.bessel_k0(0.0)
        with self.assertRaises(DomainError):
            c.bessel_i0e(-1.0)
        with self.assertRaises(DomainError):
            c.bessel_i0(800.0)
        self.assertGreater(c.bessel_i0e(800.0), 0.0)
        with self.assertRaises(DomainError):
            c.bessel_j0(math.nan)


class SeriesSolutionTests(unittest.TestCase):
    def test_table_5_2_all_rows(self) -> None:
        ref = DATA["one_term_coefficients"]
        for row in ref["rows"]:
            biot = math.inf if row["bi"] is None else row["bi"]
            for shape in ("slab", "cylinder", "sphere"):
                lam, a1, d1 = row[shape]
                got = c.one_term_coefficients(shape, biot)
                self.assertAlmostEqual(got.eigenvalue, lam, delta=ref["lambda_tolerance"],
                                       msg=(shape, biot))
                self.assertAlmostEqual(got.amplitude, a1, delta=ref["coefficient_tolerance"],
                                       msg=(shape, biot))
                self.assertAlmostEqual(got.energy_coefficient, d1,
                                       delta=ref["coefficient_tolerance"], msg=(shape, biot))

    def test_eqn_5_33_printed_terms(self) -> None:
        ref = DATA["slab_first_kind_series"]
        for i, fo in enumerate(ref["fourier_numbers"]):
            printed = []
            for n, key in ((1, "t1"), (3, "t3"), (5, "t5")):
                exact = math.exp(-(n * math.pi / 2) ** 2 * fo) / n
                value = ref[key][i]
                if value is not None:
                    self.assertAlmostEqual(exact, value, delta=ref["tolerance"])
                printed.append(0.0 if value is None else value)
            if fo >= 0.1:  # higher terms negligible at the printed precision
                approx = 4 / math.pi * (printed[0] - printed[1] + printed[2])
                series = c.series_transient_solution("slab", math.inf, fo, 0.0)
                self.assertAlmostEqual(series.theta, approx, delta=2e-3)

    def test_energy_fraction_equals_volume_average_of_theta(self) -> None:
        # D_n (volume averages) checked against direct Simpson integration of the series.
        weights = {"slab": 0, "cylinder": 1, "sphere": 2}
        for shape, power in weights.items():
            n = 200
            values = []
            for j in range(n + 1):
                x = j / n
                values.append(c.series_transient_solution(shape, 2.0, 0.05, x).theta * x**power)
            simpson = (values[0] + values[-1] + 4 * sum(values[1:-1:2])
                       + 2 * sum(values[2:-1:2])) / (3 * n)
            mean = (power + 1) * simpson
            series = c.series_transient_solution(shape, 2.0, 0.05, 0.0)
            self.assertAlmostEqual(series.mean_theta, mean, delta=2e-7, msg=shape)
            self.assertAlmostEqual(series.energy_fraction, 1 - series.mean_theta, delta=1e-15)

    def test_example_5_2_apples(self) -> None:
        ref = DATA["series_examples"]["example_5_2"]
        radius = ref["diameter_m"] / 2
        alpha = ref["k_w_m_k"] / (ref["density_kg_m3"] * ref["specific_heat_j_kg_k"])
        fourier = alpha * ref["time_s"] / radius**2
        biot = ref["h_w_m2_k"] * radius / ref["k_w_m_k"]
        self.assertAlmostEqual(fourier, ref["fourier"], delta=5e-4)
        self.assertAlmostEqual(biot, ref["biot"], delta=5e-4)
        centre = c.series_transient_solution("sphere", biot, fourier, 0.0)
        self.assertAlmostEqual(centre.theta, ref["theta_centre_chart"],
                               delta=ref["chart_read_tolerance"])
        target = bisect(lambda fo: c.series_transient_solution("sphere", biot, fo, 0.0).theta
                        - ref["theta_target"], 0.3, 3.0, iterations=60)
        self.assertAlmostEqual(target, ref["fourier_to_target_chart"],
                               delta=ref["fourier_read_tolerance"])
        phi = c.series_transient_solution("sphere", 0.5, ref["fourier_to_target_chart"]).energy_fraction
        self.assertAlmostEqual(phi, ref["phi_chart"], delta=ref["phi_read_tolerance"])

    def test_small_fourier_needs_semi_infinite_solution(self) -> None:
        with self.assertRaises(ConvergenceError):
            c.series_transient_solution("slab", 1.0, 1e-7, 0.5, max_terms=50)
        with self.assertRaises(DomainError):
            c.series_transient_solution("slab", 0.0, 0.5)
        with self.assertRaises(InputError):
            c.transient_eigenvalues("cube", 1.0, 3)

    def test_eigenvalues_are_ordered_and_satisfy_table_5_1(self) -> None:
        for biot in (0.05, 1.0, 30.0):
            lams = c.transient_eigenvalues("cylinder", biot, 12)
            self.assertTrue(all(b - a > 2.0 for a, b in zip(lams, lams[1:], strict=False)))
            for lam in lams:
                residual = lam * c.bessel_j1(lam) - biot * c.bessel_j0(lam)
                self.assertLess(abs(residual), 1e-12 * max(1.0, biot))
            for lam in c.transient_eigenvalues("sphere", biot, 12):
                self.assertAlmostEqual(lam / math.tan(lam), 1 - biot, delta=1e-11)


class LumpedAndSemiInfiniteTests(unittest.TestCase):
    def test_example_1_4_thermocouple(self) -> None:
        ref = DATA["lumped"]["example_1_4"]
        d = ref["diameter_m"]
        result = c.lumped_capacitance_response(
            initial_temperature_k=kelvin(ref["initial_c"]), ambient_temperature_k=kelvin(ref["gas_c"]),
            heat_transfer_coefficient_w_m2_k=ref["h_w_m2_k"], surface_area_m2=math.pi * d * d,
            volume_m3=math.pi * d**3 / 6, density_kg_m3=ref["density_kg_m3"],
            specific_heat_j_kg_k=ref["specific_heat_j_kg_k"], times_s=[0.0, 1.116, 3.348],
            thermal_conductivity_w_m_k=ref["k_w_m_k"])
        self.assertAlmostEqual(result.time_constant_s, ref["time_constant_s"],
                               delta=ref["tolerance_s"])
        # the book's Biot number uses L = D/2; the library's uses V/A = D/6
        self.assertAlmostEqual(result.biot_number * 3, ref["biot_with_radius"], delta=5e-6)
        self.assertTrue(result.lumped_assumption_satisfied)
        self.assertAlmostEqual(result.temperatures_k[2] - kelvin(ref["gas_c"]),
                               -180 * math.exp(-3.348 / result.time_constant_s), delta=1e-9)

    def test_example_2_8_and_biot_warning(self) -> None:
        ref = DATA["lumped"]["example_2_8"]
        volume = math.pi * ref["length_m"] * ref["diameter_m"] ** 2 / 4
        kwargs = dict(initial_temperature_k=kelvin(72.3), ambient_temperature_k=kelvin(35.0),
                      heat_transfer_coefficient_w_m2_k=ref["h_eff_w_m2_k"],
                      surface_area_m2=ref["area_m2"], volume_m3=volume,
                      density_kg_m3=ref["density_kg_m3"],
                      specific_heat_j_kg_k=ref["specific_heat_j_kg_k"], times_s=[0.0])
        result = c.lumped_capacitance_response(**kwargs)
        self.assertAlmostEqual(result.time_constant_s, ref["time_constant_s"],
                               delta=ref["tolerance_s"])
        self.assertIn("lumped_biot_unchecked", [w.code for w in result.warnings])
        bad = c.lumped_capacitance_response(**kwargs, thermal_conductivity_w_m_k=1e-3)
        self.assertFalse(bad.lumped_assumption_satisfied)
        self.assertIn("lumped_biot_too_large", [w.code for w in bad.warnings])

    def test_table_5_3_through_the_temperature_step_solution(self) -> None:
        ref = DATA["error_function"]
        for z, erf_value, erfc_value in ref["rows"]:
            result = c.semi_infinite_surface_temperature_step(
                depth_m=2 * z, time_s=1.0, thermal_diffusivity_m2_s=1.0,
                thermal_conductivity_w_m_k=1.0, initial_temperature_k=400.0,
                surface_temperature_k=300.0)
            self.assertAlmostEqual((result.temperature_k - 300.0) / 100.0, erf_value,
                                   delta=ref["tolerance"])
            self.assertAlmostEqual(math.erfc(z), erfc_value, delta=ref["tolerance"])

    def test_penetration_depth_and_example_5_4(self) -> None:
        pen = DATA["semi_infinite"]["penetration"]
        self.assertAlmostEqual(math.erf(pen["erf_argument_for_0_99"]), 0.99, delta=1e-5)
        ref = DATA["semi_infinite"]["example_5_4"]
        time = bisect(lambda t: c.semi_infinite_surface_temperature_step(
            depth_m=0.0, time_s=t, thermal_diffusivity_m2_s=ref["diffusivity_m2_s"],
            thermal_conductivity_w_m_k=1.0, initial_temperature_k=400.0,
            surface_temperature_k=300.0).penetration_depth_99_m - ref["half_thickness_m"],
            1e-4, 1.0)
        self.assertAlmostEqual(time, ref["time_s"], delta=ref["tolerance_s"])

    def test_example_5_5_finger_in_flame(self) -> None:
        ref = DATA["semi_infinite"]["example_5_5"]
        self.assertAlmostEqual(ref["h_w_m2_k"] ** 2 * ref["diffusivity_m2_s"] / ref["k_w_m_k"] ** 2,
                               ref["beta_squared_per_s"], delta=5e-5)

        def theta(t: float) -> float:
            result = c.semi_infinite_convective_surface(
                depth_m=0.0, time_s=t, thermal_diffusivity_m2_s=ref["diffusivity_m2_s"],
                thermal_conductivity_w_m_k=ref["k_w_m_k"],
                heat_transfer_coefficient_w_m2_k=ref["h_w_m2_k"], initial_temperature_k=310.0,
                fluid_temperature_k=1073.0)
            return (result.surface_temperature_k - 1073.0) / (310.0 - 1073.0)

        time = bisect(lambda t: theta(t) - ref["theta_surface"], 0.01, 2.0)
        self.assertAlmostEqual(time, ref["time_s"], delta=ref["tolerance_s"])
        beta = (1 - ref["theta_surface"]) * math.sqrt(math.pi) / 2
        self.assertAlmostEqual(beta, ref["small_beta_estimate"], delta=5e-6)

    def test_convective_solution_limits(self) -> None:
        common = dict(depth_m=0.002, time_s=5.0, thermal_diffusivity_m2_s=1e-5,
                      thermal_conductivity_w_m_k=20.0, initial_temperature_k=300.0)
        step = c.semi_infinite_surface_temperature_step(surface_temperature_k=900.0, **common)
        strong = c.semi_infinite_convective_surface(heat_transfer_coefficient_w_m2_k=1e9,
                                                    fluid_temperature_k=900.0, **common)
        self.assertAlmostEqual(strong.temperature_k, step.temperature_k, delta=1e-3)
        far = c.semi_infinite_convective_surface(heat_transfer_coefficient_w_m2_k=1e3,
                                                 fluid_temperature_k=900.0,
                                                 **{**common, "depth_m": 1.0})
        self.assertAlmostEqual(far.temperature_k, 300.0, delta=1e-9)

    def test_flux_step_profile_obeys_eqn_5_55_and_5_56(self) -> None:
        kwargs = dict(time_s=2.0, thermal_diffusivity_m2_s=4e-6, thermal_conductivity_w_m_k=15.0,
                      initial_temperature_k=300.0, surface_heat_flux_w_m2=2e5)
        surface = c.semi_infinite_surface_flux_step(depth_m=0.0, **kwargs)
        self.assertAlmostEqual(surface.temperature_k, surface.surface_temperature_k, delta=1e-9)
        self.assertAlmostEqual(surface.surface_temperature_k - 300.0,
                               2 * 2e5 / 15.0 * math.sqrt(4e-6 * 2.0 / math.pi), delta=1e-9)
        for depth in (0.001, 0.003, 0.006):
            h = 1e-7
            gradient = (c.semi_infinite_surface_flux_step(depth_m=depth + h, **kwargs).temperature_k
                        - c.semi_infinite_surface_flux_step(depth_m=depth - h, **kwargs)
                        .temperature_k) / (2 * h)
            local = c.semi_infinite_surface_flux_step(depth_m=depth, **kwargs)
            self.assertAlmostEqual(-15.0 * gradient / local.heat_flux_w_m2, 1.0, delta=1e-6)

    def test_integrated_erfc_definition_dlmf_7_18_2(self) -> None:
        below = c.integrated_erfc(25.0)  # direct closed form
        above = c.integrated_erfc(math.nextafter(25.0, 30.0))  # asymptotic tail branch
        self.assertAlmostEqual(above / below, 1.0, delta=1e-10)
        for z in (0.0, 0.4, 1.5, 4.0, 25.5, 26.0):
            n = 4000
            upper = z + 12.0
            h = (upper - z) / n
            values = [math.erfc(z + j * h) for j in range(n + 1)]
            simpson = h / 3 * (values[0] + values[-1] + 4 * sum(values[1:-1:2])
                               + 2 * sum(values[2:-1:2]))
            got = c.integrated_erfc(z)
            if z > 25.0:
                # reference integral underflows; compare with the asymptotic leading term
                lead = math.exp(-z * z) / (2 * z * z * math.sqrt(math.pi))
                self.assertAlmostEqual(got / lead, 1.0, delta=2.0 / z**2)
            else:
                self.assertAlmostEqual(got, simpson, delta=5e-12)

    def test_example_5_7_annual_wave_depth(self) -> None:
        ref = DATA["semi_infinite"]["example_5_7"]
        omega = 2 * math.pi / (365 * 24 * 3600)

        def slope(x: float) -> float:
            h = 1e-4
            f = c.semi_infinite_periodic_surface
            return (f(depth_m=x + h, time_s=0.0, thermal_diffusivity_m2_s=ref["diffusivity_m2_s"],
                      mean_temperature_k=285.0, amplitude_k=10.0, angular_frequency_rad_s=omega)
                    - f(depth_m=x - h, time_s=0.0,
                        thermal_diffusivity_m2_s=ref["diffusivity_m2_s"], mean_temperature_k=285.0,
                        amplitude_k=10.0, angular_frequency_rad_s=omega)) / (2 * h)

        depth = bisect(slope, 1.0, 4.0, iterations=60)
        self.assertAlmostEqual(depth, ref["depth_m"], delta=ref["tolerance_m"])
        self.assertAlmostEqual(depth * math.sqrt(omega / (2 * ref["diffusivity_m2_s"])),
                               ref["xi_first_minimum"], delta=5e-4)

    def test_contact_temperature_eqn_5_60(self) -> None:
        self.assertAlmostEqual(c.contact_temperature(temperature_1_k=310.0, effusivity_1=1.0,
                                                     temperature_2_k=290.0, effusivity_2=1.0),
                               300.0, delta=1e-12)
        self.assertAlmostEqual(c.contact_temperature(temperature_1_k=310.0, effusivity_1=1.0,
                                                     temperature_2_k=290.0, effusivity_2=3.0),
                               295.0, delta=1e-12)


class SteadyWallTests(unittest.TestCase):
    def test_examples_1_5_and_1_6_radiating_thermocouple(self) -> None:
        for key in ("example_1_5", "example_1_6"):
            ref = DATA["steady"][key]
            offset = ref["kelvin_offset_used_in_book"]
            result = c.solve_steady_wall(
                c.WallLayer(0.001, 10.0), c.WallGeometry("plane"),
                c.SurfaceConvection(ref["h_w_m2_k"], ref["air_c"] + offset, ref["emissivity"],
                                    ref["wall_c"] + offset))
            self.assertAlmostEqual(result.outer_surface_temperature_k - offset,
                                   ref["thermocouple_c"], delta=ref["tolerance_c"], msg=key)
            self.assertLess(abs(result.boundary_residual_w), 1e-9)

    def test_example_2_1_parabolic_profile(self) -> None:
        q, k, length, wall = 5e4, 1.4, 0.2, 300.0
        result = c.solve_steady_wall(c.WallLayer(length, k, heat_generation_w_m3=q),
                                     c.WallGeometry("plane"), c.SurfaceTemperature(wall),
                                     c.SurfaceTemperature(wall), profile_points_per_layer=10)
        for x, t in zip(result.profile_positions_m, result.profile_temperatures_k, strict=True):
            xi = x / length
            self.assertAlmostEqual((t - wall) / (q * length**2 / k), 0.5 * (xi - xi * xi),
                                   delta=1e-12)
        self.assertAlmostEqual(result.inner_heat_flux_w_m2, -q * length / 2, delta=1e-6)
        self.assertAlmostEqual(result.maximum_temperature_position_m, length / 2, delta=1e-12)

    def test_example_2_5_cylinder_with_convection(self) -> None:
        ri, ro, k, h, length = 0.01, 0.02, 2.0, 50.0, 3.0
        result = c.solve_steady_wall(
            c.WallLayer(ro - ri, k),
            c.WallGeometry("cylinder", inner_radius_m=ri, length_m=length),
            c.SurfaceConvection(h, 300.0), c.SurfaceTemperature(400.0))
        expected = 100.0 / (c.convection_resistance(h, 2 * math.pi * ro * length)
                            + c.cylindrical_wall_resistance(ri, ro, k, length))
        self.assertAlmostEqual(result.inner_heat_rate_w / expected, 1.0, delta=1e-12)
        biot = h * ro / k
        for r, t in zip(result.profile_positions_m, result.profile_temperatures_k, strict=True):
            theta = math.log(r / ri) / (1 / biot + math.log(ro / ri))  # eqn (2.24)
            self.assertAlmostEqual((t - 400.0) / (300.0 - 400.0), theta, delta=1e-12)

    def test_example_2_6_critical_radius(self) -> None:
        ref = DATA["steady"]["example_2_6"]
        ri = ref["outer_diameter_m"] / 2
        k, h = ref["k_insulation_w_m_k"], ref["h_w_m2_k"]
        critical = c.critical_insulation_radius(k, h)
        self.assertAlmostEqual(critical, ref["critical_radius_m"], delta=5e-5)
        self.assertAlmostEqual(critical / ri, ref["critical_ratio"], delta=ref["ratio_tolerance"])
        bare = c.convection_resistance(h, 2 * math.pi * ri)

        def excess(ro: float) -> float:
            return (c.cylindrical_wall_resistance(ri, ro, k) + c.convection_resistance(
                h, 2 * math.pi * ro)) - bare

        break_even = bisect(excess, critical, 10 * ri)
        self.assertAlmostEqual(break_even / ri, ref["break_even_ratio"],
                               delta=ref["ratio_tolerance"])
        self.assertAlmostEqual(break_even, ref["break_even_radius_m"], delta=5e-5)

    def test_example_2_7_resistor_with_exact_radiation(self) -> None:
        ref = DATA["steady"]["example_2_7"]
        self.assertAlmostEqual(c.radiation_heat_transfer_coefficient(316.0, ref["emissivity"]),
                               ref["h_rad_at_316_k"], delta=0.005)
        self.assertAlmostEqual(c.radiation_heat_transfer_coefficient(327.0, ref["emissivity"]),
                               ref["h_rad_at_327_k"], delta=ref["h_rad_327_discrepancy_tolerance"])
        area = ref["area_m2"]
        result = c.solve_steady_wall(
            c.WallLayer(1e-4, 50.0), c.WallGeometry("plane", area_m2=area),
            c.SurfaceConvection(ref["h_w_m2_k"], kelvin(ref["air_c"]), ref["emissivity"]),
            c.SurfaceHeatFlux(ref["power_w"] / area))
        self.assertAlmostEqual(result.outer_surface_temperature_k - CELSIUS_OFFSET_K,
                               ref["resistor_c"], delta=ref["tolerance_c"])

    def test_example_2_9_contact_resistance_share(self) -> None:
        ref = DATA["steady"]["example_2_9"]
        layer = c.WallLayer(ref["thickness_m"], ref["k_w_m_k"])
        result = c.solve_steady_wall(
            [layer, layer], c.WallGeometry("plane"), c.SurfaceTemperature(300.0),
            c.SurfaceTemperature(400.0),
            contact_conductances_w_m2_k=[ref["contact_conductance_w_m2_k"]])
        jump = result.layer_outer_temperatures_k[0] - result.layer_inner_temperatures_k[1]
        self.assertAlmostEqual(jump, result.inner_heat_flux_w_m2
                               / ref["contact_conductance_w_m2_k"], delta=1e-9)
        self.assertAlmostEqual(jump / 100.0, ref["contact_share_printed"], delta=0.015)

    def test_example_2_10_tea_kettle(self) -> None:
        ref = DATA["steady"]["example_2_10"]
        result = c.solve_steady_wall(
            c.WallLayer(ref["thickness_m"], ref["k_aluminium_w_m_k"]), c.WallGeometry("plane"),
            c.SurfaceConvection(ref["h_boiling_w_m2_k"], 373.0),
            c.SurfaceConvection(ref["h_gas_w_m2_k"], 873.0))
        self.assertAlmostEqual(result.inner_heat_flux_w_m2 / 500.0, ref["overall_u_w_m2_k"],
                               delta=ref["tolerance"])

    def test_problems_2_15_and_2_20(self) -> None:
        ref = DATA["steady"]["problem_2_15"]
        rate = (ref["sphere_c"] - ref["far_c"]) / c.spherical_wall_resistance(
            ref["diameter_m"] / 2, math.inf, ref["k_w_m_k"])
        self.assertAlmostEqual(rate, ref["heat_rate_w"], delta=ref["tolerance_w"])
        ref = DATA["steady"]["problem_2_20"]
        result = c.solve_steady_wall(
            c.WallLayer(ref["diameter_m"] / 2, ref["k_w_m_k"],
                        heat_generation_w_m3=ref["generation_w_m3"]),
            c.WallGeometry("sphere"), c.SurfaceConvection(ref["h_w_m2_k"], kelvin(ref["fluid_c"])))
        self.assertAlmostEqual(result.inner_surface_temperature_k - CELSIUS_OFFSET_K,
                               ref["centre_c"], delta=ref["tolerance_c"])
        self.assertAlmostEqual(result.outer_heat_rate_w, result.generated_heat_rate_w, delta=1e-12)

    def test_linear_conductivity_gives_mean_k_flux_problem_2_1(self) -> None:
        # Problem 2.1 (printed p. 86): for k linear in T the plane-wall flux is k(T_mean) dT / L.
        table = c.TabulatedProperty((300.0, 900.0), (20.0, 35.0))
        result = c.solve_steady_wall(c.WallLayer(0.01, table), c.WallGeometry("plane"),
                                     c.SurfaceTemperature(350.0), c.SurfaceTemperature(800.0))
        expected = table.value_at(575.0) * 450.0 / 0.01
        self.assertAlmostEqual(result.inner_heat_flux_w_m2 / expected, 1.0, delta=1e-12)

    def test_invalid_problems_fail_closed(self) -> None:
        layer = c.WallLayer(0.01, 10.0)
        with self.assertRaises(InputError):
            c.solve_steady_wall(layer, c.WallGeometry("plane"), c.SurfaceHeatFlux(1.0))
        with self.assertRaises(InputError):
            c.solve_steady_wall(layer, c.WallGeometry("sphere"), c.SurfaceTemperature(300.0),
                                c.SurfaceTemperature(400.0))
        table = c.TabulatedProperty((300.0, 400.0), (10.0, 12.0))
        with self.assertRaises(DomainError):
            c.solve_steady_wall(c.WallLayer(0.01, table), c.WallGeometry("plane"),
                                c.SurfaceTemperature(300.0), c.SurfaceTemperature(500.0))
        with self.assertRaises(DomainError):
            c.WallGeometry("plane", inner_radius_m=0.1)
        with self.assertRaises(DomainError):
            c.SurfaceConvection(0.0, 300.0)
        with self.assertRaises(InputError):
            c.solve_steady_wall([layer, layer], c.WallGeometry("plane"),
                                c.SurfaceTemperature(300.0), c.SurfaceTemperature(400.0),
                                contact_conductances_w_m2_k=[1.0, 2.0])


def _slab_run(shape: str, cells: int, steps: int, theta: float = 0.5,
              fourier: float = 0.5, biot: float = 1.0, initial=400.0):
    geometry = c.WallGeometry({"slab": "plane"}.get(shape, shape))
    return c.solve_transient_wall(
        c.WallLayer(1.0, 1.0, 1.0), geometry, c.SurfaceConvection(biot, 300.0),
        initial_temperature_k=initial, time_step_s=fourier / steps, step_count=steps,
        cells_per_layer=cells, theta=theta)


class FiniteVolumeTests(unittest.TestCase):
    def test_spatial_convergence_against_series(self) -> None:
        minimum_order = {"slab": 1.95, "cylinder": 1.75, "sphere": 1.75}
        for shape, floor in minimum_order.items():
            field_errors, surface_errors = [], []
            for cells in (10, 20, 40):
                run = _slab_run(shape, cells, 2 * cells)
                field_errors.append(max(
                    abs((t - 300.0) / 100.0
                        - c.series_transient_solution(shape, 1.0, 0.5, x).theta)
                    for x, t in zip(run.cell_centers_m, run.final_temperatures_k, strict=True)))
                exact = c.series_transient_solution(shape, 1.0, 0.5, 1.0).theta
                surface_errors.append(abs((run.outer_surface_temperatures_k[-1] - 300.0) / 100.0
                                          - exact))
                self.assertLess(abs(run.energy_residual_j), 1e-9 * abs(run.stored_energy_change_j))
            for order in observed_orders(field_errors):
                self.assertGreater(order, floor, msg=(shape, field_errors))
            for order in observed_orders(surface_errors):
                self.assertGreater(order, 1.95, msg=(shape, surface_errors))
            self.assertLess(field_errors[-1], 5e-4)

    def test_temporal_order_from_an_eigenmode_start(self) -> None:
        lam = c.transient_eigenvalues("slab", 1.0, 1)[0]
        cells = 20
        start = [300.0 + 100.0 * math.cos(lam * (j + 0.5) / cells) for j in range(cells)]
        for theta, expected in ((0.5, 2.0), (1.0, 1.0)):
            fields = [_slab_run("slab", cells, steps, theta=theta, initial=start)
                      .final_temperatures_k for steps in (10, 20, 40, 80)]
            gaps = [max(abs(a - b) for a, b in zip(fields[i], fields[i + 1], strict=True))
                    for i in range(3)]
            self.assertAlmostEqual(observed_orders(gaps)[-1], expected, delta=0.06, msg=theta)

    def test_semi_infinite_flux_surface_temperature(self) -> None:
        exact = c.semi_infinite_surface_flux_step(
            depth_m=0.0, time_s=0.01, thermal_diffusivity_m2_s=1.0, thermal_conductivity_w_m_k=1.0,
            initial_temperature_k=300.0, surface_heat_flux_w_m2=1000.0).surface_temperature_k
        errors = []
        for cells in (50, 100, 200):
            run = c.solve_transient_wall(
                c.WallLayer(1.0, 1.0, 1.0), c.WallGeometry("plane"), c.SurfaceHeatFlux(0.0),
                c.SurfaceHeatFlux(1000.0), initial_temperature_k=300.0,
                time_step_s=0.01 / cells, step_count=cells, cells_per_layer=cells)
            errors.append(abs(run.inner_surface_temperatures_k[-1] - exact))
            self.assertAlmostEqual(run.energy_residual_j, 0.0, delta=1e-9)
            self.assertAlmostEqual(run.boundary_energy_in_j, 1000.0 * 0.01, delta=1e-9)
        for order in observed_orders(errors):
            self.assertGreater(order, 1.95)

    def test_nonlinear_multilayer_converges_to_kirchhoff_steady_solution(self) -> None:
        conductivity = c.TabulatedProperty((250.0, 500.0, 900.0, 1500.0),
                                           (390.0, 370.0, 340.0, 300.0))
        capacity = c.TabulatedProperty((250.0, 1500.0), (3.2e6, 4.2e6))
        layers = [c.WallLayer(0.002, conductivity, capacity, heat_generation_w_m3=2e7),
                  c.WallLayer(0.003, 20.0, 4.0e6)]
        geometry = c.WallGeometry("cylinder", inner_radius_m=0.02)
        inner = c.SurfaceConvection(20000.0, 1200.0)
        outer = c.SurfaceConvection(50.0, 300.0, 0.8, 280.0)
        steady = c.solve_steady_wall(layers, geometry, outer, inner,
                                     contact_conductances_w_m2_k=[2e4])
        errors = []
        for cells in (4, 8, 16):
            run = c.solve_transient_wall(
                layers, geometry, outer, inner, contact_conductances_w_m2_k=[2e4],
                initial_temperature_k=300.0, time_step_s=2.0, step_count=100,
                cells_per_layer=cells, theta=1.0)
            errors.append(abs(run.outer_surface_temperatures_k[-1]
                              - steady.outer_surface_temperature_k))
            self.assertLess(abs(run.energy_residual_j), 1e-9 * abs(run.stored_energy_change_j))
            self.assertAlmostEqual(run.inner_heat_rates_w[-1] / steady.inner_heat_rate_w, 1.0,
                                   delta=1e-5)
        for order in observed_orders(errors):
            self.assertAlmostEqual(order, 2.0, delta=0.1)

    def test_energy_audit_during_a_nonlinear_transient(self) -> None:
        conductivity = c.TabulatedProperty((200.0, 2000.0), (15.0, 30.0))
        capacity = c.TabulatedProperty((200.0, 2000.0), (3.5e6, 5.0e6))
        layers = [c.WallLayer(0.004, conductivity, capacity),
                  c.WallLayer(0.001, 1.0, 1.5e6, heat_generation_w_m3=-1e6)]
        run = c.solve_transient_wall(
            layers, c.WallGeometry("sphere", inner_radius_m=0.05),
            c.SurfaceConvection(10.0, 300.0, 0.9, 250.0), c.SurfaceHeatFlux(5e5),
            contact_conductances_w_m2_k=[5e3], initial_temperature_k=300.0,
            time_step_s=0.05, step_count=60, cells_per_layer=(10, 4), snapshot_every=20)
        self.assertGreater(run.stored_energy_change_j, 0.0)
        self.assertLess(abs(run.energy_residual_j), 1e-9 * run.stored_energy_change_j)
        self.assertEqual(len(run.snapshot_times_s), 4)
        self.assertGreater(run.max_picard_iterations, 1)

    def test_single_layer_kirchhoff_flux_is_exact_at_steady_state(self) -> None:
        table = c.TabulatedProperty((250.0, 600.0, 1000.0), (400.0, 330.0, 200.0))
        for kind, inner_radius in (("plane", 0.0), ("cylinder", 0.01), ("sphere", 0.01)):
            geometry = c.WallGeometry(kind, inner_radius_m=inner_radius)
            layer = c.WallLayer(0.01, table, 3e6)
            steady = c.solve_steady_wall(layer, geometry, c.SurfaceTemperature(300.0),
                                         c.SurfaceTemperature(900.0))
            run = c.solve_transient_wall(layer, geometry, c.SurfaceTemperature(300.0),
                                         c.SurfaceTemperature(900.0), initial_temperature_k=600.0,
                                         time_step_s=5.0, step_count=300, cells_per_layer=5,
                                         theta=1.0)
            self.assertAlmostEqual(run.inner_heat_rates_w[-1] / steady.inner_heat_rate_w, 1.0,
                                   delta=1e-12, msg=kind)

    def test_steady_state_reached_by_chaining_runs(self) -> None:
        layer = c.WallLayer(0.01, 20.0, 4e6)
        geometry = c.WallGeometry("plane")
        first = c.solve_transient_wall(layer, geometry, c.SurfaceTemperature(300.0),
                                       c.SurfaceTemperature(500.0), initial_temperature_k=300.0,
                                       time_step_s=1.0, step_count=200, cells_per_layer=10,
                                       theta=1.0)
        second = c.solve_transient_wall(layer, geometry, c.SurfaceTemperature(300.0),
                                        c.SurfaceTemperature(500.0),
                                        initial_temperature_k=first.final_temperatures_k,
                                        time_step_s=1.0, step_count=200, cells_per_layer=10,
                                        theta=1.0)
        for x, t in zip(second.cell_centers_m, second.final_temperatures_k, strict=True):
            self.assertAlmostEqual(t, 500.0 - 200.0 * x / 0.01, delta=1e-6)

    def test_invalid_runs_fail_closed(self) -> None:
        layer = c.WallLayer(0.01, 20.0, 4e6)
        kwargs = dict(initial_temperature_k=300.0, time_step_s=1.0, step_count=2)
        with self.assertRaises(DomainError):
            c.solve_transient_wall(layer, c.WallGeometry("plane"), c.SurfaceTemperature(300.0),
                                   theta=0.3, **kwargs)
        with self.assertRaises(InputError):
            c.solve_transient_wall(c.WallLayer(0.01, 20.0), c.WallGeometry("plane"),
                                   c.SurfaceTemperature(300.0), **kwargs)
        with self.assertRaises(InputError):
            c.solve_transient_wall(layer, c.WallGeometry("plane"), c.SurfaceTemperature(300.0),
                                   cells_per_layer=0, **kwargs)
        table = c.TabulatedProperty((250.0, 350.0), (10.0, 12.0))
        with self.assertRaises(DomainError):
            c.solve_transient_wall(c.WallLayer(0.01, table, 4e6), c.WallGeometry("plane"),
                                   c.SurfaceTemperature(600.0), **kwargs)


class FinTests(unittest.TestCase):
    def test_example_4_8_tip_temperatures(self) -> None:
        ref = DATA["fins"]["example_4_8"]
        common = dict(thermal_conductivity_w_m_k=ref["k_w_m_k"],
                      heat_transfer_coefficient_w_m2_k=ref["h_w_m2_k"],
                      diameter_m=ref["diameter_m"], length_m=ref["length_m"],
                      base_temperature_k=kelvin(ref["wall_c"]),
                      fluid_temperature_k=kelvin(ref["air_c"]))
        tipped = c.pin_fin(tip_heat_transfer_coefficient_w_m2_k=ref["h_w_m2_k"], **common)
        insulated = c.pin_fin(**common)
        excess = ref["wall_c"] - ref["air_c"]
        self.assertAlmostEqual(tipped.fin_parameter_per_m * ref["length_m"], ref["mL"], delta=5e-5)
        self.assertAlmostEqual(ref["h_w_m2_k"] * ref["length_m"] / ref["k_w_m_k"],
                               ref["bi_axial"], delta=5e-5)
        self.assertAlmostEqual((tipped.tip_temperature_k - kelvin(ref["air_c"])) / excess,
                               ref["theta_tip_convective"], delta=5e-4)
        self.assertAlmostEqual((insulated.tip_temperature_k - kelvin(ref["air_c"])) / excess,
                               ref["theta_tip_adiabatic"], delta=5e-5)
        self.assertAlmostEqual(tipped.tip_temperature_k - CELSIUS_OFFSET_K,
                               ref["tip_c_convective"], delta=0.05)
        self.assertAlmostEqual(insulated.tip_temperature_k - CELSIUS_OFFSET_K,
                               ref["tip_c_adiabatic"], delta=0.05)
        # eqn (4.55): effectiveness = efficiency * surface area / cross-section
        area = math.pi * ref["diameter_m"] ** 2 / 4
        self.assertAlmostEqual(insulated.effectiveness,
                               insulated.efficiency * insulated.convecting_area_m2 / area,
                               delta=1e-12)

    def test_example_4_9_lead_wires(self) -> None:
        ref = DATA["fins"]["example_4_9"]
        lead = c.pin_fin(thermal_conductivity_w_m_k=ref["k_w_m_k"],
                         heat_transfer_coefficient_w_m2_k=ref["h_eff_w_m2_k"],
                         diameter_m=ref["diameter_m"], length_m=1.0, base_temperature_k=330.0,
                         fluid_temperature_k=300.0)
        self.assertAlmostEqual(lead.thermal_resistance_k_w, ref["fin_resistance_k_w"], delta=5.0)
        conductance = (2 / lead.thermal_resistance_k_w
                       + ref["area_m2"] * (ref["h_rad_w_m2_k"] + ref["h_conv_w_m2_k"]))
        self.assertAlmostEqual(1 / conductance, ref["equivalent_resistance_k_w"], delta=0.5)
        self.assertAlmostEqual(ref["air_c"] + ref["power_w"] / conductance, ref["resistor_c"],
                               delta=0.05)

    def test_example_4_10_annular_fins_on_a_tube(self) -> None:
        ref = DATA["fins"]["example_4_10"]
        r1, r2 = ref["tube_diameter_m"] / 2, ref["fin_diameter_m"] / 2
        thickness = ref["thickness_m"]
        self.assertAlmostEqual(math.sqrt(ref["h_fin_w_m2_k"] * (r2 - r1) ** 3
                                         / (ref["k_w_m_k"] * (r2 - r1) * thickness)),
                               ref["chart_parameter"], delta=5e-4)
        fin = c.annular_fin(thermal_conductivity_w_m_k=ref["k_w_m_k"],
                            heat_transfer_coefficient_w_m2_k=ref["h_fin_w_m2_k"],
                            thickness_m=thickness, inner_radius_m=r1, outer_radius_m=r2,
                            base_temperature_k=kelvin(ref["water_c"]),
                            fluid_temperature_k=kelvin(ref["air_c"]))
        self.assertAlmostEqual(fin.efficiency, ref["efficiency_chart"],
                               delta=ref["chart_read_tolerance"])
        excess = ref["water_c"] - ref["air_c"]
        bare = math.pi * ref["tube_diameter_m"] * ref["h_tube_w_m2_k"] * excess
        self.assertAlmostEqual(bare, ref["heat_rate_bare_w_m"], delta=0.5)
        fins_per_m = 1 / ref["spacing_m"]
        total = (bare * (ref["spacing_m"] - thickness) / ref["spacing_m"]
                 + fins_per_m * fin.heat_rate_w)
        self.assertAlmostEqual(total, ref["heat_rate_finned_w_m"],
                               delta=ref["finned_tolerance_w_m"])

    def test_annular_fin_matches_independent_finite_differences(self) -> None:
        k, h, thickness, a, b = 180.0, 60.0, 0.002, 0.01, 0.05
        m2 = h / (k * thickness / 2)
        for h_edge in (0.0, 60.0):
            nodes = 4000
            dr = (b - a) / nodes
            radii = [a + i * dr for i in range(nodes + 1)]
            # (1/r)(r psi')' = m^2 psi, psi(a) = 1, psi'(b) = -(h_e/k) psi(b): second-order FD
            lower, diag, upper, rhs = ([0.0] * (nodes + 1) for _ in range(4))
            diag[0], rhs[0] = 1.0, 1.0
            for i in range(1, nodes):
                r = radii[i]
                lower[i] = (r - dr / 2) / (r * dr * dr)
                upper[i] = (r + dr / 2) / (r * dr * dr)
                diag[i] = -(lower[i] + upper[i]) - m2
            # ghost node: psi_{N+1} = psi_{N-1} - 2 dr (h_e/k) psi_N
            r = radii[-1]
            west = (r - dr / 2) / (r * dr * dr)
            east = (r + dr / 2) / (r * dr * dr)
            lower[-1] = west + east
            diag[-1] = -(west + east) - m2 - east * 2 * dr * h_edge / k
            psi = c._thomas(lower, diag, upper, rhs)
            gradient = (-3 * psi[0] + 4 * psi[1] - psi[2]) / (2 * dr)
            rate = -k * 2 * math.pi * a * thickness * 50.0 * gradient
            fin = c.annular_fin(thermal_conductivity_w_m_k=k, heat_transfer_coefficient_w_m2_k=h,
                                thickness_m=thickness, inner_radius_m=a, outer_radius_m=b,
                                base_temperature_k=350.0, fluid_temperature_k=300.0,
                                edge_heat_transfer_coefficient_w_m2_k=h_edge)
            self.assertAlmostEqual(fin.heat_rate_w / rate, 1.0, delta=2e-6, msg=h_edge)
            self.assertAlmostEqual((fin.tip_temperature_k - 300.0) / 50.0, psi[-1], delta=2e-6)

    def test_annular_fin_tends_to_straight_fin(self) -> None:
        k, h, thickness, length = 200.0, 100.0, 0.001, 0.02
        straight = c.straight_fin(thermal_conductivity_w_m_k=k, heat_transfer_coefficient_w_m2_k=h,
                                  perimeter_m=2.0, cross_section_area_m2=thickness,
                                  length_m=length, base_temperature_k=400.0,
                                  fluid_temperature_k=300.0)
        previous = math.inf
        for radius in (0.01, 1.0, 100.0, 1e4):
            fin = c.annular_fin(thermal_conductivity_w_m_k=k, heat_transfer_coefficient_w_m2_k=h,
                                thickness_m=thickness, inner_radius_m=radius,
                                outer_radius_m=radius + length, base_temperature_k=400.0,
                                fluid_temperature_k=300.0)
            gap = abs(fin.efficiency - straight.efficiency)
            self.assertLess(gap, previous)
            previous = gap
        self.assertLess(previous, 1e-6)
        # large m b must not overflow the scaled Bessel evaluation
        big = c.annular_fin(thermal_conductivity_w_m_k=k, heat_transfer_coefficient_w_m2_k=h,
                            thickness_m=thickness, inner_radius_m=20.0, outer_radius_m=21.0,
                            base_temperature_k=400.0, fluid_temperature_k=300.0)
        self.assertGreater(big.heat_rate_w, 0.0)

    def test_corrected_length_and_validation(self) -> None:
        self.assertAlmostEqual(c.corrected_fin_length(0.02, 0.002), 0.021, delta=1e-15)
        with self.assertRaises(DomainError):
            c.straight_fin(thermal_conductivity_w_m_k=1.0, heat_transfer_coefficient_w_m2_k=1.0,
                           perimeter_m=1.0, cross_section_area_m2=1.0, length_m=1.0,
                           base_temperature_k=300.0, fluid_temperature_k=300.0)
        thick = c.pin_fin(thermal_conductivity_w_m_k=1.0, heat_transfer_coefficient_w_m2_k=100.0,
                          diameter_m=0.02, length_m=0.05, base_temperature_k=350.0,
                          fluid_temperature_k=300.0)
        self.assertIn("fin_transverse_biot_large", [w.code for w in thick.warnings])


if __name__ == "__main__":
    unittest.main()
