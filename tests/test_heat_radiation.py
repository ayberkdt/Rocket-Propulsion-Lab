"""Source-traceable checks of ``rocket_propulsion.heat_transfer.radiation``.

Reference values and locators are in ``tests/reference_data/heat_radiation.json``.
Where the textbook used sigma = 5.67e-8 the comparison rescales by that value
or states the measured deviation in the assertion tolerance. Independent
oracles: composite Simpson quadrature of Planck's law (stdlib), the
algebraic identities of view factors, and energy conservation.
"""

from __future__ import annotations

import itertools
import json
import math
import unittest
from pathlib import Path

from rocket_propulsion.core.errors import DomainError, InputError
from rocket_propulsion.heat_transfer import radiation as rad
from rocket_propulsion.thermochemistry.rocket_performance import RocketStation

REFERENCE = json.loads((Path(__file__).parent / "reference_data" / "heat_radiation.json")
                       .read_text(encoding="utf-8"))
SIGMA = rad.STEFAN_BOLTZMANN_W_M2_K4
ATM = 101_325.0


def _relative(actual: float, expected: float) -> float:
    return abs(actual / expected - 1.0)


def _simpson_planck_fraction(temperature: float, low: float, high: float, n: int = 4000) -> float:
    """Independent oracle: Simpson's rule on Eq. (1.30) in log-wavelength."""

    a, b = math.log(low), math.log(high)
    h = (b - a) / n
    total = 0.0
    for i in range(n + 1):
        wavelength = math.exp(a + i * h)
        weight = 1 if i in (0, n) else (4 if i % 2 else 2)
        total += weight * rad.planck_spectral_emissive_power(wavelength, temperature) * wavelength
    return total * h / 3.0 / (SIGMA * temperature**4)


class ConstantTests(unittest.TestCase):
    def test_constants_match_nist_and_lienhard(self) -> None:
        ref = REFERENCE["constants"]
        self.assertLess(_relative(SIGMA, ref["nist_stefan_boltzmann_w_m2_k4"]), 1e-9)
        self.assertLess(_relative(rad.WIEN_DISPLACEMENT_CONSTANT_M_K, ref["nist_wien_wavelength_m_k"]), 1e-9)
        self.assertLess(_relative(rad.SECOND_RADIATION_CONSTANT_M_K, ref["nist_second_radiation_constant_m_k"]),
                        1e-9)
        self.assertLess(_relative(SIGMA, ref["lienhard_stefan_boltzmann_w_m2_k4"]), 1e-7)
        self.assertAlmostEqual(rad.WIEN_DISPLACEMENT_CONSTANT_M_K * 1e6, ref["lienhard_wien_um_k"], places=2)

    def test_planck_integral_reproduces_stefan_boltzmann(self) -> None:
        for temperature in (300.0, 1000.0, 3500.0):
            self.assertLess(abs(_simpson_planck_fraction(temperature, 1e-8, 1.0) - 1.0), 1e-9)

    def test_wien_peak_is_the_maximum(self) -> None:
        peak = rad.wien_peak_wavelength(2000.0)
        value = rad.planck_spectral_emissive_power(peak, 2000.0)
        for factor in (0.999, 1.001):
            self.assertLess(rad.planck_spectral_emissive_power(peak * factor, 2000.0), value)


class BandFractionTests(unittest.TestCase):
    def test_series_matches_quadrature_on_both_branches(self) -> None:
        for temperature, upper in ((1000.0, 1.5e-6), (1000.0, 3e-6), (1000.0, 10e-6), (1000.0, 13e-6),
                                   (1000.0, 15e-6), (300.0, 200e-6)):
            oracle = _simpson_planck_fraction(temperature, 1e-8, upper)
            self.assertAlmostEqual(rad.blackbody_fraction(upper * temperature), oracle, delta=1e-10)

    def test_branch_switch_is_continuous(self) -> None:
        switch = rad.SECOND_RADIATION_CONSTANT_M_K  # z = 1
        below = rad.blackbody_fraction(switch * (1 - 1e-13))
        above = rad.blackbody_fraction(switch * (1 + 1e-13))
        # dF/d(ln lambdaT) = 0.09 at z = 1, so a 2e-13 step moves F by ~1.8e-14.
        self.assertLess(abs(above - below), 5e-14)

    def test_limits_and_monotonicity(self) -> None:
        self.assertEqual(rad.blackbody_fraction(math.inf), 1.0)
        self.assertLess(rad.blackbody_fraction(1e-5), 1e-200)
        values = [rad.blackbody_fraction(x * 1e-6) for x in range(200, 100_000, 997)]
        self.assertTrue(all(b >= a for a, b in itertools.pairwise(values)))
        with self.assertRaises(DomainError):
            rad.blackbody_fraction(0.0)

    def test_lienhard_band_statements(self) -> None:
        ref = REFERENCE["blackbody_statements"]
        above = 1.0 - rad.blackbody_fraction(REFERENCE["constants"]["lienhard_wien_um_k"] * 1e-6)
        self.assertAlmostEqual(above, ref["fraction_above_peak_approx"], delta=ref["fraction_above_peak_tolerance"])
        for band in ref["bands_containing_ninety_percent"]:
            fraction = rad.blackbody_band_fraction(band["lower_um"] * 1e-6, band["upper_um"] * 1e-6,
                                                   band["temperature_k"])
            self.assertAlmostEqual(fraction, band["fraction"], delta=ref["band_tolerance"])


class ViewFactorTests(unittest.TestCase):
    def test_example_10_1_slit_algebra(self) -> None:
        ref = REFERENCE["example_10_1"]
        jet_area = math.pi * ref["jet_diameter_m"]
        slit = ref["shield_diameter_m"] * math.sin(math.radians(ref["slit_angle_deg"] / 2.0))
        self.assertAlmostEqual(slit, ref["slit_area_m2_per_m"], delta=5e-6)
        f_slit_jet = rad.reciprocal_view_factor(jet_area, ref["book_f_jet_room_used"], slit)
        self.assertAlmostEqual(f_slit_jet, ref["f_slit_jet"], delta=5e-6)
        self.assertAlmostEqual(1.0 - f_slit_jet, ref["f_slit_shield"], delta=5e-5)
        shield_area = math.pi * ref["shield_diameter_m"] * (360.0 - ref["slit_angle_deg"]) / 360.0
        f_shield_room = rad.reciprocal_view_factor(slit, 1.0 - f_slit_jet, shield_area)
        self.assertAlmostEqual(f_shield_room, ref["f_shield_room"], delta=2e-5)  # book rounds 0.01294, 0.9393
        sigma = ref["book_sigma"]
        q_room = shield_area * f_shield_room * sigma * (ref["shield_temperature_k"] ** 4 - ref["room_temperature_k"] ** 4)
        self.assertLess(_relative(q_room, ref["q_shield_room_w_m"]), 2e-3)  # measured 0.16 %
        q_jet = jet_area * (330 / 360) * sigma * (ref["jet_temperature_k"] ** 4 - ref["shield_temperature_k"] ** 4)
        self.assertLess(_relative(q_jet, ref["q_jet_shield_w_m"]), 1e-3)

    def test_example_10_2_against_chart_reading(self) -> None:
        ref = REFERENCE["example_10_2"]
        f14 = rad.view_factor_perpendicular_rectangles(1.0, 1.0, 1.0)
        combined = rad.view_factor_perpendicular_rectangles(0.5, 0.5, 1.0)
        self.assertAlmostEqual(f14, ref["f_1_4"], delta=5e-3)
        # Closed form gives 0.2406; the book read 0.245 off Fig. 10.9 (1.8 % chart-reading difference).
        self.assertAlmostEqual(combined, ref["f_combined_chart_read"], delta=5e-3)
        self.assertAlmostEqual(combined - f14, ref["f_1_2"], delta=5e-3)

    def test_examples_10_3_and_10_4_frustum(self) -> None:
        ref = REFERENCE["example_10_3_10_4"]
        enclosure = rad.frustum_enclosure_view_factors(ref["heater_radius_m"], ref["opening_radius_m"],
                                                       ref["height_m"])
        self.assertAlmostEqual(enclosure.view_factors[0][1], ref["f_heater_opening"], delta=5e-4)
        self.assertAlmostEqual(enclosure.view_factors[0][2], ref["f_heater_shield"], delta=5e-4)
        self.assertAlmostEqual(enclosure.areas_m2[2], ref["shield_area_m2"], delta=5e-6)
        self.assertAlmostEqual(enclosure.view_factors[2][0], ref["f_shield_heater"], delta=5e-5)
        q = (enclosure.areas_m2[0] * enclosure.view_factors[0][2] * ref["book_sigma"]
             * (ref["heater_temperature_k"] ** 4 - ref["shield_temperature_k"] ** 4))
        self.assertLess(_relative(q, ref["q_heater_shield_w"]), 1e-3)

    def test_frustum_algebra_matches_howell_c111_and_c78(self) -> None:
        for small, large, length in ((0.05, 0.1, 0.2), (0.1, 0.3, 0.5), (0.3, 1.2, 2.0), (0.2, 0.2, 1.0)):
            enclosure = rad.frustum_enclosure_view_factors(small, large, length)
            wall = enclosure.view_factors[2]
            self.assertAlmostEqual(wall[1], rad.view_factor_frustum_to_base(small, large, length), places=12)
            self.assertAlmostEqual(wall[0], rad.view_factor_frustum_to_base(large, small, length), places=12)
            summation, reciprocity = rad.view_factor_matrix_residuals(enclosure.areas_m2, enclosure.view_factors)
            self.assertLess(max(summation, reciprocity), 1e-14)
        cylinder = rad.frustum_enclosure_view_factors(0.2, 0.2, 1.0)
        self.assertAlmostEqual(cylinder.view_factors[2][2], rad.view_factor_cylinder_interior_to_itself(0.2, 1.0),
                               places=12)

    def test_disk_and_rectangle_identities(self) -> None:
        # Howell C-40 (equal disks) is C-41 with R1 = R2.
        for radius in (0.1, 1.0, 5.0):
            x = (2 * radius**2 + 1) / radius**2
            c40 = 0.5 * (x - math.sqrt(x * x - 4))
            self.assertAlmostEqual(rad.view_factor_coaxial_disks(radius, radius, 1.0), c40, places=13)
        for a, b, c in ((1.0, 2.0, 0.5), (3.0, 0.2, 1.0)):
            self.assertAlmostEqual(rad.view_factor_parallel_rectangles(a, b, c),
                                   rad.view_factor_parallel_rectangles(b, a, c), places=13)
        self.assertGreater(rad.view_factor_parallel_rectangles(1e3, 1e3, 1.0), 0.998)
        # Reciprocity of perpendicular rectangles: w F(w->h) = h F(h->w).
        self.assertAlmostEqual(2.0 * rad.view_factor_perpendicular_rectangles(2.0, 0.7, 1.3),
                               0.7 * rad.view_factor_perpendicular_rectangles(0.7, 2.0, 1.3), places=13)

    def test_two_dimensional_table(self) -> None:
        self.assertAlmostEqual(rad.view_factor_parallel_strips(1.0, 2.4), REFERENCE["example_10_9"]["f_12"],
                               places=12)
        self.assertAlmostEqual(rad.view_factor_hinged_strips(math.pi), 0.0, places=15)
        # Perpendicular strips of equal width: Table 10.2 #3 equals #2 at alpha = 90 deg.
        self.assertAlmostEqual(rad.view_factor_perpendicular_strips(1.0, 1.0),
                               rad.view_factor_hinged_strips(math.pi / 2), places=14)
        self.assertAlmostEqual(rad.view_factor_triangular_duct(0.5, 0.3, 0.4), REFERENCE["example_10_10"]["f_12"],
                               places=14)
        self.assertAlmostEqual(rad.view_factor_parallel_cylinders(1.0, 0.0), 0.5 - 1 / math.pi, places=14)
        self.assertLess(rad.view_factor_strip_to_cylinder(0.1, -1e6, 1e6, 1.0), 0.1 * math.pi / 2e6 + 1e-12)
        ((_, f12), (f21, f22)) = rad.view_factor_concentric_spheres(1.0, 2.0)
        self.assertEqual((f12, f21, f22), (1.0, 0.25, 0.75))
        with self.assertRaises(DomainError):
            rad.view_factor_triangular_duct(1.0, 1.0, 3.0)


class EnclosureTests(unittest.TestCase):
    def _strips(self, case: str) -> rad.EnclosureSolution:
        ref = REFERENCE["example_10_9"]
        f12 = rad.view_factor_parallel_strips(ref["strip_width_m"], ref["separation_m"])
        a3 = 2 * ref["separation_m"]
        f31 = (1 - f12) / a3
        matrix = [[0, f12, 1 - f12], [f12, 0, 1 - f12], [f31, f31, 1 - 2 * f31]]
        third = (rad.GraySurface("surroundings", a3, 1.0, temperature_k=ref["t_surroundings_k"]) if case == "a"
                 else rad.GraySurface("reflector", a3, 0.5, net_heat_flux_w_m2=0.0))
        return rad.solve_gray_enclosure([rad.GraySurface("1", 1.0, ref["e1"], temperature_k=ref["t1_k"]),
                                         rad.GraySurface("2", 1.0, ref["e2"], temperature_k=ref["t2_k"]),
                                         third], matrix)

    def test_example_10_9(self) -> None:
        ref = REFERENCE["example_10_9"]
        f12 = ref["f_12"]
        a = self._strips("a")
        # Book used sigma = 5.67e-8 and rounded e_b; measured deviations <= 0.2 %.
        self.assertLess(_relative(a.radiosities_w_m2[0], ref["case_a"]["b1"]), 2e-3)
        self.assertLess(_relative(a.radiosities_w_m2[1], ref["case_a"]["b2"]), 2e-3)
        q12 = (a.radiosities_w_m2[0] - a.radiosities_w_m2[1]) * f12
        self.assertLess(_relative(q12, ref["case_a"]["q_12_w_m"]), 3e-3)
        b = self._strips("b")
        for index, key in enumerate(("b1", "b2", "b3")):
            self.assertLess(_relative(b.radiosities_w_m2[index], ref["case_b"][key]), 2e-3)
        self.assertLess(_relative(b.net_heat_rates_w[0], ref["case_b"]["q_12_w_m"]), 3e-3)
        self.assertAlmostEqual(b.temperatures_k[2], ref["case_b"]["t3_k"], delta=0.5)
        for solution in (a, b):
            self.assertLess(abs(solution.energy_residual_w), 1e-9)

    def test_example_10_10(self) -> None:
        ref = REFERENCE["example_10_10"]
        w = ref["widths_m"]
        f = rad.view_factor_triangular_duct
        matrix = [[0, f(w[0], w[1], w[2]), f(w[0], w[2], w[1])],
                  [f(w[1], w[0], w[2]), 0, f(w[1], w[2], w[0])],
                  [f(w[2], w[0], w[1]), f(w[2], w[1], w[0]), 0]]
        surfaces = [rad.GraySurface(str(i), w[i], ref["emissivities"][i], temperature_k=ref["temperatures_k"][i])
                    for i in range(3)]
        solution = rad.solve_gray_enclosure(surfaces, matrix)
        q_book_sigma = solution.net_heat_rates_w[0] * ref["book_sigma"] / SIGMA
        self.assertLess(_relative(q_book_sigma, ref["q_1_w_m"]), 3e-3)  # measured 0.22 % (book rounding)
        self.assertLess(abs(solution.energy_residual_w), 1e-9)

    def test_enclosed_body_identity_and_two_surface_formula(self) -> None:
        a1, a2, e1, e2, t1, t2 = 0.4, 3.0, 0.3, 0.6, 900.0, 400.0
        f21 = a1 / a2
        solution = rad.solve_gray_enclosure(
            [rad.GraySurface("1", a1, e1, temperature_k=t1), rad.GraySurface("2", a2, e2, temperature_k=t2)],
            [[0.0, 1.0], [f21, 1.0 - f21]])
        transfer = rad.enclosed_body_transfer_factor(e1, e2, a1 / a2)
        expected = a1 * transfer * SIGMA * (t1**4 - t2**4)
        self.assertAlmostEqual(solution.net_heat_rates_w[0] / expected, 1.0, places=12)
        circuit = rad.two_surface_heat_rate(temperature_1_k=t1, temperature_2_k=t2, area_1_m2=a1, area_2_m2=a2,
                                            emissivity_1=e1, emissivity_2=e2, view_factor_12=1.0)
        self.assertAlmostEqual(circuit / expected, 1.0, places=12)
        self.assertAlmostEqual(rad.enclosed_body_transfer_factor(e1, e2, 1.0),
                               rad.parallel_plate_transfer_factor(e1, e2), places=15)

    def test_fail_closed_inputs(self) -> None:
        good = [rad.GraySurface("1", 1.0, 0.5, temperature_k=500.0), rad.GraySurface("2", 1.0, 0.5, temperature_k=300.0)]
        with self.assertRaises(DomainError):
            rad.solve_gray_enclosure(good, [[0.0, 0.9], [0.9, 0.0]])
        with self.assertRaises(InputError):
            rad.solve_gray_enclosure([rad.GraySurface("1", 1.0, 0.5, net_heat_flux_w_m2=1.0),
                                      rad.GraySurface("2", 1.0, 0.5, net_heat_flux_w_m2=-1.0)],
                                     [[0.0, 1.0], [1.0, 0.0]])
        with self.assertRaises(InputError):
            rad.solve_gray_enclosure([rad.GraySurface("1", 1.0, 0.5), good[1]], [[0.0, 1.0], [1.0, 0.0]])
        with self.assertRaises(DomainError):
            rad.solve_gray_enclosure([rad.GraySurface("1", 1.0, 0.0, temperature_k=1.0), good[1]],
                                     [[0.0, 1.0], [1.0, 0.0]])

    def test_example_10_8_shield(self) -> None:
        ref = REFERENCE["example_10_8"]
        scale = ref["book_sigma"] / SIGMA
        line_area = math.pi * ref["line_diameter_m"]
        bare = rad.shielded_heat_rate(inner_temperature_k=ref["line_temperature_k"],
                                      outer_temperature_k=ref["chamber_temperature_k"], inner_area_m2=line_area,
                                      inner_emissivity=ref["emissivity"], outer_area_m2=None, outer_emissivity=1.0)
        self.assertLess(_relative(-bare.heat_rate_w * scale, ref["q_gain_unshielded_w_m"]), 1e-3)
        shielded = rad.shielded_heat_rate(
            inner_temperature_k=ref["line_temperature_k"], outer_temperature_k=ref["chamber_temperature_k"],
            inner_area_m2=line_area, inner_emissivity=ref["emissivity"], outer_area_m2=None, outer_emissivity=1.0,
            shields=[rad.RadiationShield(math.pi * ref["shield_diameter_m"], ref["emissivity"], ref["emissivity"])])
        self.assertLess(_relative(-shielded.heat_rate_w * scale, ref["q_gain_shielded_w_m"]), 2e-3)
        self.assertAlmostEqual(shielded.shield_temperatures_k[0], ref["shield_temperature_k"], delta=0.5)

    def test_frustum_transfer_factor_matches_enclosure_solver(self) -> None:
        small, large, length, e = 0.12, 0.6, 1.1, 0.8
        t_small, t_large, _ = rad.frustum_interior_transfer_factors(small_radius_m=small, large_radius_m=large,
                                                                    length_m=length, wall_emissivity=e)
        enclosure = rad.frustum_enclosure_view_factors(small, large, length)
        a = enclosure.areas_m2
        solution = rad.solve_gray_enclosure([rad.GraySurface("throat", a[0], 1.0, temperature_k=900.0),
                                             rad.GraySurface("exit", a[1], 1.0, temperature_k=4.0),
                                             rad.GraySurface("wall", a[2], e, temperature_k=1400.0)],
                                            enclosure.view_factors)
        expected = SIGMA * (t_small * (1400.0**4 - 900.0**4) + t_large * (1400.0**4 - 4.0**4))
        self.assertAlmostEqual(solution.net_heat_fluxes_w_m2[2] / expected, 1.0, places=12)


class EquilibriumTests(unittest.TestCase):
    def test_lienhard_examples_1_5_and_1_6(self) -> None:
        ref = REFERENCE["examples_1_5_1_6"]
        offset = ref["kelvin_offset_used"]
        for emissivity, reading in ((1.0, ref["black_reading_c"]), (ref["gray_emissivity"], ref["gray_reading_c"])):
            result = rad.radiative_equilibrium_wall_temperature(
                sinks=[rad.RadiativeSink(emissivity, ref["wall_c"] + offset)],
                heat_transfer_coefficient_w_m2_k=ref["h_w_m2_k"], recovery_temperature_k=ref["air_c"] + offset)
            self.assertAlmostEqual(result.wall_temperature_k - offset, reading, delta=0.05)
            self.assertLess(abs(result.residual_w_m2), 1e-8)

    def test_pure_radiation_closed_form(self) -> None:
        q, e = 2.5e5, 0.85
        result = rad.radiative_equilibrium_wall_temperature(sinks=[rad.RadiativeSink(e, 0.0)],
                                                            imposed_heat_flux_w_m2=q)
        self.assertAlmostEqual(result.wall_temperature_k / (q / (e * SIGMA)) ** 0.25, 1.0, places=13)
        two_sided = rad.radiative_equilibrium_wall_temperature(
            sinks=[rad.RadiativeSink(e, 0.0), rad.RadiativeSink(e, 0.0)], imposed_heat_flux_w_m2=q)
        self.assertAlmostEqual(two_sided.wall_temperature_k / result.wall_temperature_k, 0.5**0.25, places=13)
        with self.assertRaises(InputError):
            rad.radiative_equilibrium_wall_temperature(sinks=[rad.RadiativeSink(e, 0.0)],
                                                       heat_transfer_coefficient_w_m2_k=10.0)


class GasRadiationTests(unittest.TestCase):
    def test_mean_beam_length_table_equals_four_v_over_a(self) -> None:
        d = 0.7
        cases = {
            "sphere": (math.pi * d**3 / 6, math.pi * d**2),
            "cube": (d**3, 6 * d**2),
            "cylinder_length_equals_diameter": (math.pi * d**3 / 4, 1.5 * math.pi * d**2),
            "infinite_cylinder": (math.pi * d**2 / 4, math.pi * d),
            "infinite_parallel_plates": (d, 2.0),
        }
        for geometry, (volume, area) in cases.items():
            self.assertAlmostEqual(rad.mean_beam_length_for_geometry(geometry, d),
                                   rad.mean_beam_length(volume, area), places=14)
        self.assertAlmostEqual(rad.mean_beam_length(1.0, 1.0, optically_thick=True), 3.6)

    def test_wsgg_weights_are_physical_over_the_fitted_range(self) -> None:
        for species in (rad.BORDBAR_2020_CO2, rad.BORDBAR_2020_H2O):
            for temperature in range(300, 2401, 25):
                weights = rad.wsgg_weights(species, float(temperature))
                self.assertTrue(all(w >= 0.0 for w in weights))
                self.assertAlmostEqual(sum(weights), 1.0, places=12)
        sums = [sum(rad.wsgg_weights(rad.SADEGHI_2021_SOOT, float(t))) for t in range(300, 3001, 25)]
        # Published Table B9 with the Eq. (10) clipping: measured 0.939 (2900 K) to 1.010 (825 K).
        self.assertTrue(0.93 < min(sums) and max(sums) < 1.011)
        with self.assertRaises(DomainError):
            rad.wsgg_weights(rad.BORDBAR_2020_H2O, 2500.0)

    def test_superposition_identity(self) -> None:
        t, length = 1500.0, 0.8
        mix = rad.wsgg_mixture_emissivity(temperature_k=t, path_length_m=length, total_pressure_pa=ATM,
                                          h2o_mole_fraction=0.2, co2_mole_fraction=0.1)
        e_w = rad.wsgg_species_emissivity(rad.BORDBAR_2020_H2O, t, 0.2 * length).value
        e_c = rad.wsgg_species_emissivity(rad.BORDBAR_2020_CO2, t, 0.1 * length).value
        self.assertAlmostEqual(mix.value, e_w + e_c - e_w * e_c, places=14)
        self.assertEqual(len(mix.gray_gas_weights), 25)
        self.assertAlmostEqual(sum(mix.gray_gas_weights), 1.0, places=12)

    def test_absorptivity_equals_emissivity_at_equal_temperatures(self) -> None:
        kwargs = {"path_length_m": 0.5, "total_pressure_pa": ATM, "h2o_mole_fraction": 0.3}
        emissivity = rad.wsgg_mixture_emissivity(temperature_k=1100.0, **kwargs).value
        absorptivity = rad.wsgg_mixture_absorptivity(gas_temperature_k=1100.0, wall_temperature_k=1100.0,
                                                     **kwargs).value
        self.assertAlmostEqual(emissivity, absorptivity, places=15)

    def test_example_10_11_cross_model(self) -> None:
        ref = REFERENCE["example_10_11"]
        length = rad.mean_beam_length_for_geometry("infinite_cylinder", ref["diameter_m"])
        self.assertEqual(length, ref["mean_beam_length_m"])
        x = ref["co2_pressure_atm"] / ref["total_pressure_atm"]
        eps = rad.wsgg_mixture_emissivity(temperature_k=ref["gas_temperature_k"], path_length_m=length,
                                          total_pressure_pa=ATM, co2_mole_fraction=x).value
        alpha = rad.wsgg_mixture_absorptivity(gas_temperature_k=ref["gas_temperature_k"],
                                              wall_temperature_k=ref["wall_temperature_k"], path_length_m=length,
                                              total_pressure_pa=ATM, co2_mole_fraction=x).value
        # Different models (pure-CO2 WSGG vs. Alberti chart, p_a -> 0): measured +1.0 % and +12.6 %.
        self.assertLess(_relative(eps, ref["eps_g_chart_read"]), 0.03)
        self.assertLess(_relative(alpha, ref["alpha_g"]), 0.15)
        q = math.pi * ref["diameter_m"] * rad.gas_to_black_wall_heat_flux(
            gas_emissivity=eps, gas_absorptivity=alpha, gas_temperature_k=ref["gas_temperature_k"],
            wall_temperature_k=ref["wall_temperature_k"])
        self.assertLess(_relative(q * ref["book_sigma"] / SIGMA, ref["q_w_m"]), 0.03)

    def test_fail_closed_pressure_path_and_soot(self) -> None:
        base = {"temperature_k": 1500.0, "path_length_m": 0.3, "h2o_mole_fraction": 0.5}
        with self.assertRaises(DomainError):
            rad.wsgg_mixture_emissivity(total_pressure_pa=20 * ATM, **base)
        flagged = rad.wsgg_mixture_emissivity(total_pressure_pa=20 * ATM, allow_pressure_extrapolation=True, **base)
        self.assertEqual(flagged.warnings[0].severity, "critical")
        with self.assertRaises(DomainError):
            rad.wsgg_mixture_emissivity(temperature_k=1500.0, path_length_m=0.01, total_pressure_pa=ATM,
                                        h2o_mole_fraction=0.1)
        sooty = rad.wsgg_mixture_emissivity(total_pressure_pa=ATM, soot_volume_fraction=1e-6, **base)
        clean = rad.wsgg_mixture_emissivity(total_pressure_pa=ATM, **base)
        self.assertGreater(sooty.value, clean.value)
        self.assertEqual(sooty.warnings[0].code, "wsgg_soot_path_unstated")
        self.assertGreater(rad.soot_planck_mean_absorption_coefficient(1500.0, 1e-6), 0.0)

    def test_rocket_station_bridge(self) -> None:
        station = RocketStation("exit", 40.0, 0.8 * ATM, 1800.0, 0.1, 0.02, 0.0, 0.0, 4000.0,
                                (("H2O", 0.6), ("H2", 0.35), ("CO2", 0.05)))
        result = rad.station_gas_emissivity(station, 0.5, allow_pressure_extrapolation=True)
        direct = rad.wsgg_mixture_emissivity(temperature_k=1800.0, path_length_m=0.5, total_pressure_pa=0.8 * ATM,
                                             h2o_mole_fraction=0.6, co2_mole_fraction=0.05,
                                             allow_pressure_extrapolation=True)
        self.assertEqual(result.value, direct.value)
        hot = RocketStation("chamber", None, 7e6, 3500.0, 1.0, 0.02, 0.0, 0.0, 0.0, (("H2O", 0.7),))
        with self.assertRaises(DomainError):
            rad.station_gas_emissivity(hot, 0.1, allow_pressure_extrapolation=True)


if __name__ == "__main__":
    unittest.main()
