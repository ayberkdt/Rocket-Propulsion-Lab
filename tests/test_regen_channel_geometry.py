"""Variable channel geometry, helical channels and the coolant momentum pressure drop.

* A profile with constant columns must reproduce ``CoolingChannels`` exactly.
* Every station-to-station pressure drop is recomputed here from the reported
  states: the Darcy drop of NASA/TM-1998-206313 (Wadel) Eq. (14) with the two
  stations' mean density, velocity ``G/rho`` and hydraulic diameter, plus the
  momentum change from a control volume on a stream tube of varying area,
  ``A dp = -mdot dv`` (Euler's ``dp = -rho v dv``), integrated with the mean
  ``rho v``, plus at a channel-count change the GFSSP sudden-area-change loss;
  densities come from the reference equation of state at the reported
  ``(T, p)``. For constant density the momentum term must reduce exactly to
  Bernoulli, and the GFSSP expansion loss plus Bernoulli to the momentum-theory
  (Borda-Carnot) static rise ``rho v2 (v1 - v2)``.
* A helix lengthens the path by ``1/cos(beta)`` and narrows the land.
* Narrowing the channels at the throat raises the coolant mass flux where the
  heat flux peaks: the peak hot-wall temperature falls and the pressure drop
  rises, the trade Wadel studies.
"""

from __future__ import annotations

import itertools
import math
import unittest
from functools import cache

from rocket_propulsion.core.errors import DomainError, InputError
from rocket_propulsion.geometry import generate_nozzle_contour
from rocket_propulsion.heat_transfer import convection
from rocket_propulsion.heat_transfer.regen_chamber import (
    COOLANT_MACH_WARNING,
    ChannelProfile,
    Coolant,
    CoolingChannels,
    hot_gas_from_propellant,
    solve_regen_chamber,
    sudden_area_change_coefficient,
)
from rocket_propulsion.thermodynamics.real_fluids import reference_fluid

PRESSURE = 7.0e6


@cache
def gas():
    return hot_gas_from_propellant(propellant_pair="LOX_CH4", chamber_pressure_pa=PRESSURE,
                                   throat_radius_of_curvature_m=0.05)


@cache
def contour():
    return generate_nozzle_contour(contour="conical", throat_area_m2=0.01, area_ratio=10.0, gamma=gas().gamma,
                                   chamber_radius_ratio=2.0, station_count=41)


def coolant(pressure: float = 1.5e7) -> Coolant:
    return Coolant(fluid="methane", inlet_temperature_k=120.0, inlet_pressure_pa=pressure,
                   mass_flow_kg_s=PRESSURE * 0.01 / gas().characteristic_velocity_m_s / 4.4)


def solve(channels):
    return solve_regen_chamber(contour=contour(), gas=gas(), channels=channels, coolant=coolant())


def span() -> tuple[float, float]:
    xs = [s.x_m for s in contour().stations]
    return min(xs), max(xs)


def throat_x() -> float:
    return min(contour().stations, key=lambda s: s.radius_m).x_m


CONSTANT = CoolingChannels(count=100, width_m=1.5e-3, height_m=5.0e-3, hot_wall_thickness_m=1.0e-3, roughness_m=1e-6)


def narrowed(helix: float = 0.0, throat_width: float = 1.0e-3) -> ChannelProfile:
    """1.5 mm x 5 mm channels narrowing to ``throat_width`` at the throat."""
    x0, x1 = span()
    xt = throat_x()
    xs = (x0, xt - 0.06, xt, xt + 0.06, x1)
    widths = (1.5e-3, 1.5e-3, throat_width, 1.5e-3, 1.5e-3)
    return ChannelProfile(x_m=xs, width_m=widths, height_m=(5e-3,) * 5, hot_wall_thickness_m=(1e-3,) * 5,
                          count=100, roughness_m=1e-6, helix_angle_deg=helix)


class ProfileEquivalence(unittest.TestCase):
    def test_constant_density_area_change_is_reversible(self) -> None:
        # trapezoidal integral of rho v dv with rho fixed is exactly Bernoulli's rho (v2^2 - v1^2)/2
        rho = 300.0
        for g_a, g_b in ((1e4, 2e4), (2e4, 1e4), (1.5e4, 1.5e4)):
            v_a, v_b = g_a / rho, g_b / rho
            self.assertAlmostEqual(0.5 * (g_a + g_b) * (v_b - v_a), rho * (v_b ** 2 - v_a ** 2) / 2,
                                   delta=1e-9 * rho * max(v_a, v_b) ** 2)

    def test_constant_profile_reproduces_cooling_channels(self) -> None:
        x0, x1 = span()
        profile = ChannelProfile(x_m=(x0, 0.5 * (x0 + x1), x1), width_m=(1.5e-3,) * 3, height_m=(5.0e-3,) * 3,
                                 hot_wall_thickness_m=(1.0e-3,) * 3, count=100, roughness_m=1e-6)
        a, b = solve(CONSTANT), solve(profile)
        self.assertEqual(a.stations, b.stations)
        self.assertEqual(a.total_heat_load_w, b.total_heat_load_w)

    def test_station_geometry_is_the_interpolated_profile(self) -> None:
        profile = narrowed()
        result = solve(profile)
        for station in result.stations:
            width, height, thickness = profile.geometry_at(station.x_m)
            self.assertEqual((station.channel_width_m, station.channel_height_m, station.hot_wall_thickness_m),
                             (width, height, thickness))
        narrowest = min(result.stations, key=lambda s: s.channel_width_m)
        self.assertAlmostEqual(narrowest.x_m, throat_x(), delta=0.02)


def count_at(channels, x: float) -> int:
    if isinstance(channels, ChannelProfile):
        return channels.count_at(x)
    return channels.count


class PressureDropRecomputed(unittest.TestCase):
    def check(self, channels, helix: float) -> None:
        result = solve(channels)
        fluid = reference_fluid("methane")
        flow = coolant().mass_flow_kg_s
        stations = list(reversed(result.stations))  # counterflow: exit to injector
        self.assertLess(abs(result.energy_balance_residual_w), 1e-9 * result.total_heat_load_w)
        momentum_total = 0.0
        for a, b in itertools.pairwise(stations):
            rho_a = fluid.state_tp(a.coolant_temperature_k, a.coolant_pressure_pa).density_kg_m3
            rho_b = fluid.state_tp(b.coolant_temperature_k, b.coolant_pressure_pa).density_kg_m3
            n_a, n_b = count_at(channels, a.x_m), count_at(channels, b.x_m)
            area_a, area_b = a.channel_width_m * a.channel_height_m, b.channel_width_m * b.channel_height_m
            d_a = 2 * area_a / (a.channel_width_m + a.channel_height_m)
            d_b = 2 * area_b / (b.channel_width_m + b.channel_height_m)
            g_a, g_b = flow / (n_a * area_a), flow / (n_b * area_b)
            v_a, v_b = g_a / rho_a, g_b / rho_b
            friction = convection.chen_friction_factor(a.coolant_reynolds, 1e-6 / d_a).value
            path = math.hypot(b.x_m - a.x_m, b.radius_m - a.radius_m) / math.cos(math.radians(helix))
            viscous = friction * path / (0.5 * (d_a + d_b)) * 0.5 * (rho_a + rho_b) * (0.5 * (v_a + v_b)) ** 2 / 2
            momentum = 0.5 * (g_a + g_b) * (v_b - v_a)
            junction = 0.0
            if n_a != n_b:
                k1 = sudden_area_change_coefficient((n_a * area_a) / (n_b * area_b), a.coolant_reynolds, friction)
                junction = k1 * rho_a * v_a ** 2 / 2
            momentum_total += momentum
            drop = a.coolant_pressure_pa - b.coolant_pressure_pa
            with self.subTest(x=b.x_m):
                self.assertLess(abs(drop - viscous - momentum - junction), 1e-6 * max(abs(drop), 1.0))
                self.assertLess(abs(a.coolant_velocity_m_s - v_a), 1e-9 * a.coolant_velocity_m_s)
        # heating lowers the density by about a factor of five: acceleration is a large share
        self.assertGreater(momentum_total, 0.2 * result.coolant_pressure_drop_pa)
        return result

    def test_bifurcated_channels(self) -> None:
        x0, x1 = span()
        xt = throat_x()
        # 100 channels through chamber and throat, split into 200 half-width ones in the
        # expansion where the circumference is large (narrow channels at the throat would choke)
        xb = xt + 0.15
        self.assertLess(xb, x1)
        profile = ChannelProfile(x_m=(x0, xb, x1), width_m=(1.5e-3, 0.75e-3, 0.75e-3),
                                 height_m=(5e-3, 5e-3, 5e-3), hot_wall_thickness_m=(1e-3,) * 3,
                                 count=(100, 200, 200), roughness_m=1e-6)
        self.assertEqual(profile.count_at(xb - 1e-9), 100)
        self.assertEqual(profile.count_at(xb), 200)
        result = solve_regen_chamber(contour=contour(), gas=gas(), channels=profile, coolant=coolant())
        self.assertIn("channel_count_change", {w.code for w in result.warnings})

    def test_constant_channels(self) -> None:
        self.check(CONSTANT, 0.0)

    def test_bifurcation_against_the_recomputation(self) -> None:
        x0, x1 = span()
        xb = throat_x() + 0.12  # merge 200 -> 100 downstream of the throat, seen in counterflow as a split
        profile = ChannelProfile(x_m=(x0, xb, x1), width_m=(0.75e-3, 1.5e-3, 1.5e-3),
                                 height_m=(5e-3,) * 3, hot_wall_thickness_m=(1e-3,) * 3,
                                 count=(200, 100, 100), roughness_m=1e-6)
        result = self.check(profile, 0.0)
        self.assertIn("channel_count_change", {w.code for w in result.warnings})

    def test_varying_helical_channels(self) -> None:
        self.check(narrowed(helix=30.0), 30.0)


class SuddenAreaChange(unittest.TestCase):
    def test_coefficients_against_gfssp_and_momentum_theory(self) -> None:
        self.assertEqual(sudden_area_change_coefficient(1.0, 1e5, 0.02), 0.0)
        # GFSSP Eq. (39), reduction A1/A2 = 2, turbulent
        self.assertAlmostEqual(sudden_area_change_coefficient(2.0, 1e5, 0.02), (0.6 + 0.48 * 0.02) * 2 * 1, places=14)
        # Eq. (38) laminar reduction and Eq. (40) laminar expansion
        self.assertAlmostEqual(sudden_area_change_coefficient(2.0, 1000.0, 0.064), (1.2 + 0.16) * 3, places=14)
        self.assertAlmostEqual(sudden_area_change_coefficient(0.5, 1000.0, 0.064), 2 * 0.75, places=14)
        for sigma in (0.1, 0.25, 0.5, 0.8, 0.99):
            with self.subTest(sigma=sigma):
                k1 = sudden_area_change_coefficient(sigma, 1e6, 0.0)
                self.assertAlmostEqual(k1, (1 - sigma) ** 2, places=14)  # Borda-Carnot
                # Bernoulli minus the loss equals the control-volume momentum result rho v2 (v1 - v2)
                rho, v1 = 400.0, 50.0
                v2 = sigma * v1
                static_rise = rho * (v1 ** 2 - v2 ** 2) / 2 - k1 * rho * v1 ** 2 / 2
                self.assertAlmostEqual(static_rise, rho * v2 * (v1 - v2), delta=1e-9 * rho * v1 ** 2)
        # losses grow with the mismatch and with friction
        expansions = [sudden_area_change_coefficient(s, 1e5, 0.02) for s in (0.9, 0.7, 0.5, 0.3)]
        reductions = [sudden_area_change_coefficient(s, 1e5, 0.02) for s in (1.1, 1.5, 2.0, 3.0)]
        self.assertEqual(expansions, sorted(expansions))
        self.assertEqual(reductions, sorted(reductions))
        self.assertGreater(sudden_area_change_coefficient(0.5, 1e5, 0.05), sudden_area_change_coefficient(0.5, 1e5, 0.01))

    def test_hostile_values_are_refused(self) -> None:
        for args in ((0.0, 1e5, 0.02), (-1.0, 1e5, 0.02), (math.nan, 1e5, 0.02), (0.5, 0.0, 0.02),
                     (0.5, math.inf, 0.02), (0.5, 1e5, -0.01), (0.5, 1e5, math.nan)):
            with self.subTest(args=args), self.assertRaises(DomainError):
                sudden_area_change_coefficient(*args)

    def test_profile_counts_are_validated(self) -> None:
        good = {"x_m": (0.0, 1.0), "width_m": (1e-3, 1e-3), "height_m": (3e-3, 3e-3),
                "hot_wall_thickness_m": (1e-3, 1e-3)}
        for counts in ((50,), (50, 0), (50, True), (50, 2.5), (50, 100, 200)):
            with self.subTest(counts=counts), self.assertRaises(DomainError):
                ChannelProfile(**good, count=counts)
        profile = ChannelProfile(**good, count=[40, 80])
        self.assertEqual(profile.count, (40, 80))
        self.assertEqual((profile.count_at(0.0), profile.count_at(0.999), profile.count_at(1.0)), (40, 40, 80))


class HelixAndHarccTrends(unittest.TestCase):
    def test_helix_lengthens_the_path_and_narrows_the_land(self) -> None:
        axial = solve(CONSTANT)
        helical = solve(CoolingChannels(count=100, width_m=1.5e-3, height_m=5.0e-3, hot_wall_thickness_m=1.0e-3,
                                        roughness_m=1e-6, helix_angle_deg=35.0))
        self.assertGreater(helical.coolant_pressure_drop_pa, axial.coolant_pressure_drop_pa)
        # same gas-side surface and nearly the same wall: the heat load changes by a few percent at most
        self.assertLess(abs(helical.total_heat_load_w / axial.total_heat_load_w - 1.0), 0.05)
        too_steep = CoolingChannels(count=100, width_m=1.5e-3, height_m=5.0e-3, hot_wall_thickness_m=1.0e-3,
                                    helix_angle_deg=72.0)
        with self.assertRaises(DomainError):
            solve(too_steep)  # the channel-normal pitch no longer fits 100 channels at the throat

    def test_narrower_throat_channels_cool_the_wall_and_cost_pressure(self) -> None:
        results = [solve(narrowed(throat_width=w)) for w in (1.5e-3, 1.25e-3, 1.0e-3)]
        peaks = [r.maximum_hot_wall_temperature_k for r in results]
        drops = [r.coolant_pressure_drop_pa for r in results]
        self.assertEqual(peaks, sorted(peaks, reverse=True))
        self.assertEqual(drops, sorted(drops))
        self.assertGreater(peaks[0] - peaks[-1], 30.0)
        # the 1.5 mm profile is the constant design
        self.assertEqual(results[0].maximum_hot_wall_temperature_k, solve(CONSTANT).maximum_hot_wall_temperature_k)

    def test_near_choking_design_is_flagged(self) -> None:
        tight = CoolingChannels(count=100, width_m=1.5e-3, height_m=3.0e-3, hot_wall_thickness_m=1.0e-3,
                                roughness_m=1e-6)
        result = solve_regen_chamber(contour=contour(), gas=gas(), channels=tight, coolant=coolant(2.0e7))
        self.assertGreater(max(s.coolant_mach for s in result.stations), COOLANT_MACH_WARNING)
        self.assertIn("coolant_high_mach", {w.code for w in result.warnings})
        with self.assertRaises(DomainError):
            solve_regen_chamber(contour=contour(), gas=gas(), channels=tight, coolant=coolant(1.5e7))


class TwoDimensionalWall(unittest.TestCase):
    def test_2d_wall_model_conserves_energy_and_bounds_the_fin_model(self) -> None:
        fin = solve(CONSTANT)
        wall = solve_regen_chamber(contour=contour(), gas=gas(), channels=CONSTANT, coolant=coolant(), wall_model="2d")
        self.assertEqual(wall.wall_model, "2d")
        self.assertLess(abs(wall.energy_balance_residual_w), 1e-9 * wall.total_heat_load_w)
        self.assertLess(abs(wall.total_heat_load_w / fin.total_heat_load_w - 1.0), 0.05)
        for a, b in zip(fin.stations, wall.stations, strict=True):
            with self.subTest(x=a.x_m):
                self.assertGreaterEqual(b.hot_wall_temperature_k, b.hot_wall_mean_temperature_k)
                self.assertGreater(b.hot_wall_mean_temperature_k, b.cold_wall_temperature_k)
                self.assertEqual(a.hot_wall_temperature_k, a.hot_wall_mean_temperature_k)  # fin: one value
                self.assertLess(abs(b.heat_flux_w_m2 / a.heat_flux_w_m2 - 1.0), 0.05)
        # the fin model's uniform hot wall hides the hot spot above the channel centre
        self.assertGreater(wall.maximum_hot_wall_temperature_k, fin.maximum_hot_wall_temperature_k)
        self.assertIn("two-dimensional wall conduction", next(w.message for w in wall.warnings
                                                               if w.code == "rte_class_model"))
        with self.assertRaises(InputError):
            solve_regen_chamber(contour=contour(), gas=gas(), channels=CONSTANT, coolant=coolant(), wall_model="3d")


class HostileProfiles(unittest.TestCase):
    def test_malformed_profiles_are_refused(self) -> None:
        good = {"x_m": (0.0, 1.0), "width_m": (1e-3, 1e-3), "height_m": (3e-3, 3e-3),
                "hot_wall_thickness_m": (1e-3, 1e-3), "count": 50}
        for change in ({"x_m": (1.0, 0.0)}, {"x_m": (0.0, 0.0)}, {"x_m": (0.0,)}, {"width_m": (1e-3,)},
                       {"width_m": (1e-3, -1e-3)}, {"height_m": (3e-3, math.nan)}, {"x_m": (0.0, math.inf)},
                       {"count": 0}, {"count": True}, {"count": 2.5}, {"roughness_m": -1e-6},
                       {"helix_angle_deg": 75.0}, {"helix_angle_deg": -1.0}, {"helix_angle_deg": math.nan}):
            with self.subTest(change=change), self.assertRaises(DomainError):
                ChannelProfile(**{**good, **change})

    def test_contour_outside_the_profile_is_refused(self) -> None:
        x0, x1 = span()
        short = ChannelProfile(x_m=(x0, 0.5 * (x0 + x1)), width_m=(1.5e-3,) * 2, height_m=(5e-3,) * 2,
                               hot_wall_thickness_m=(1e-3,) * 2, count=100)
        with self.assertRaises(DomainError):
            solve(short)
        with self.assertRaises(DomainError):
            solve(CoolingChannels(count=100, width_m=1.5e-3, height_m=5e-3, hot_wall_thickness_m=1e-3,
                                  helix_angle_deg=80.0))


if __name__ == "__main__":
    unittest.main()
