"""Variable channel geometry, helical channels and the coolant momentum pressure drop.

* A profile with constant columns must reproduce ``CoolingChannels`` exactly.
* Every station-to-station pressure drop is recomputed here from the reported
  states with NASA/TM-1998-206313 (Wadel) Eqs. (14)-(16): the Darcy drop with
  the two stations' mean density, velocity ``G/rho`` and hydraulic diameter,
  plus the momentum drop ``G_2^2/rho_2 - G_1^2/rho_1``, densities from the
  reference equation of state at the reported ``(T, p)``.
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

from rocket_propulsion.core.errors import DomainError
from rocket_propulsion.geometry import generate_nozzle_contour
from rocket_propulsion.heat_transfer import convection
from rocket_propulsion.heat_transfer.regen_chamber import (
    COOLANT_MACH_WARNING,
    ChannelProfile,
    Coolant,
    CoolingChannels,
    hot_gas_from_propellant,
    solve_regen_chamber,
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


class PressureDropRecomputed(unittest.TestCase):
    def check(self, channels, helix: float) -> None:
        result = solve(channels)
        fluid = reference_fluid("methane")
        flow = coolant().mass_flow_kg_s
        stations = list(reversed(result.stations))  # counterflow: exit to injector
        self.assertLess(abs(result.energy_balance_residual_w), 1e-9 * result.total_heat_load_w)
        for a, b in itertools.pairwise(stations):
            rho_a = fluid.state_tp(a.coolant_temperature_k, a.coolant_pressure_pa).density_kg_m3
            rho_b = fluid.state_tp(b.coolant_temperature_k, b.coolant_pressure_pa).density_kg_m3
            area_a, area_b = a.channel_width_m * a.channel_height_m, b.channel_width_m * b.channel_height_m
            d_a = 2 * area_a / (a.channel_width_m + a.channel_height_m)
            d_b = 2 * area_b / (b.channel_width_m + b.channel_height_m)
            g_a, g_b = flow / (100 * area_a), flow / (100 * area_b)
            friction = convection.chen_friction_factor(a.coolant_reynolds, 1e-6 / d_a).value
            path = math.hypot(b.x_m - a.x_m, b.radius_m - a.radius_m) / math.cos(math.radians(helix))
            mean_v = 0.5 * (g_a / rho_a + g_b / rho_b)
            viscous = friction * path / (0.5 * (d_a + d_b)) * 0.5 * (rho_a + rho_b) * mean_v ** 2 / 2
            momentum = g_b ** 2 / rho_b - g_a ** 2 / rho_a
            drop = a.coolant_pressure_pa - b.coolant_pressure_pa
            with self.subTest(x=b.x_m):
                self.assertLess(abs(drop - viscous - momentum), 1e-6 * max(abs(drop), 1.0))
                self.assertLess(abs(a.coolant_velocity_m_s - g_a / rho_a), 1e-9 * a.coolant_velocity_m_s)
        momentum_total = sum(
            (flow / (100 * b.channel_width_m * b.channel_height_m)) ** 2
            / fluid.state_tp(b.coolant_temperature_k, b.coolant_pressure_pa).density_kg_m3
            - (flow / (100 * a.channel_width_m * a.channel_height_m)) ** 2
            / fluid.state_tp(a.coolant_temperature_k, a.coolant_pressure_pa).density_kg_m3
            for a, b in itertools.pairwise(stations))
        # heating lowers the density by about a factor of five: acceleration is a large share
        self.assertGreater(momentum_total, 0.2 * result.coolant_pressure_drop_pa)

    def test_constant_channels(self) -> None:
        self.check(CONSTANT, 0.0)

    def test_varying_helical_channels(self) -> None:
        self.check(narrowed(helix=30.0), 30.0)


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
