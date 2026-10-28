"""Coupled regenerative-cooling solver: consistency, conservation, trends and refinement."""

from __future__ import annotations

import unittest
from functools import cache

from rocket_propulsion.core.errors import DomainError, InputError
from rocket_propulsion.geometry import generate_nozzle_contour
from rocket_propulsion.heat_transfer.regen_chamber import (
    Coolant,
    CoolingChannels,
    bartz_coefficient,
    hot_gas_from_propellant,
    solve_regen_chamber,
)
from rocket_propulsion.thermodynamics import calculate_bartz_heat_transfer
from rocket_propulsion.thermodynamics.materials import get_wall_material

PRESSURE = 7.0e6
CHANNELS = CoolingChannels(count=100, width_m=1.5e-3, height_m=5.0e-3, hot_wall_thickness_m=1.0e-3,
                           roughness_m=1.0e-6)


@cache
def _gas():
    return hot_gas_from_propellant(propellant_pair="LOX_CH4", chamber_pressure_pa=PRESSURE,
                                   throat_radius_of_curvature_m=0.05)


@cache
def _contour(stations: int = 41):
    return generate_nozzle_contour(contour="conical", throat_area_m2=0.01, area_ratio=10.0,
                                   gamma=_gas().gamma, chamber_radius_ratio=2.0, station_count=stations)


def _fuel_flow() -> float:
    return PRESSURE * 0.01 / _gas().characteristic_velocity_m_s / 4.4


@cache
def _solve(stations: int = 41, flow_scale: float = 1.0, material: str = "grcop-84"):
    return solve_regen_chamber(
        contour=_contour(stations), gas=_gas(), channels=CHANNELS, wall_material=material,
        coolant=Coolant(fluid="methane", inlet_temperature_k=120.0, inlet_pressure_pa=1.5e7,
                        mass_flow_kg_s=flow_scale * _fuel_flow()))


class GasSideConsistencyTests(unittest.TestCase):
    def test_station_bartz_equals_the_profile_implementation_at_uniform_wall_temperature(self) -> None:
        gas = _gas()
        profile = calculate_bartz_heat_transfer(
            contour=_contour(), chamber_pressure_pa=gas.chamber_pressure_pa,
            chamber_temperature_k=gas.chamber_temperature_k, wall_temperature_k=700.0,
            characteristic_velocity_m_s=gas.characteristic_velocity_m_s, dynamic_viscosity_pa_s=gas.viscosity_pa_s,
            specific_heat_j_kg_k=gas.specific_heat_j_kg_k, prandtl_number=gas.prandtl_number, gamma=gas.gamma,
            throat_radius_of_curvature_m=gas.throat_radius_of_curvature_m)
        for station, reference in zip(_contour().stations, profile.stations, strict=True):
            h, recovery = bartz_coefficient(gas, throat_radius_m=_contour().throat_radius_m,
                                            area_ratio=station.area_ratio, mach=station.mach,
                                            wall_temperature_k=700.0)
            with self.subTest(x=station.x_m):
                self.assertAlmostEqual(h / reference.gas_side_coefficient_w_m2_k, 1.0, delta=1.0e-12)
                self.assertAlmostEqual(recovery / reference.adiabatic_wall_temperature_k, 1.0, delta=1.0e-12)


class ConservationTests(unittest.TestCase):
    def test_coolant_enthalpy_rise_equals_absorbed_heat(self) -> None:
        result = _solve()
        self.assertLess(abs(result.energy_balance_residual_w), 1.0e-9 * result.total_heat_load_w)

    def test_every_station_satisfies_its_series_resistance_network(self) -> None:
        result = _solve()
        conductivity = get_wall_material("grcop-84").thermal_conductivity.value
        for station in result.stations:
            with self.subTest(x=station.x_m):
                conducted = conductivity * (station.hot_wall_temperature_k - station.cold_wall_temperature_k) \
                    / CHANNELS.hot_wall_thickness_m
                self.assertAlmostEqual(conducted / station.heat_flux_w_m2, 1.0, delta=1.0e-9)
                self.assertGreater(station.hot_wall_temperature_k, station.cold_wall_temperature_k)
                self.assertGreater(station.cold_wall_temperature_k, station.coolant_temperature_k)
                self.assertLess(station.hot_wall_temperature_k, station.adiabatic_wall_temperature_k)


class PhysicalTrendTests(unittest.TestCase):
    def test_counterflow_coolant_heats_and_loses_pressure_toward_the_injector(self) -> None:
        stations = _solve().stations
        temperatures = [item.coolant_temperature_k for item in stations]
        pressures = [item.coolant_pressure_pa for item in stations]
        self.assertEqual(temperatures, sorted(temperatures, reverse=True))
        self.assertEqual(pressures, sorted(pressures))

    def test_peak_heat_flux_sits_at_the_throat(self) -> None:
        peak = max(_solve().stations, key=lambda item: item.heat_flux_w_m2)
        self.assertAlmostEqual(peak.mach, 1.0, delta=0.15)

    def test_more_coolant_cools_the_wall_and_a_poor_conductor_heats_it(self) -> None:
        nominal = _solve().maximum_hot_wall_temperature_k
        self.assertLess(_solve(flow_scale=1.5).maximum_hot_wall_temperature_k, nominal)
        self.assertGreater(_solve(material="inconel-718").maximum_hot_wall_temperature_k, nominal)

    def test_wall_limit_margin_is_reported_against_the_material_limit(self) -> None:
        result = _solve()
        limit = get_wall_material("grcop-84").temperature_limit()[0].value
        self.assertEqual(result.wall_limit_temperature_k, limit)
        hottest = min(item.wall_limit_margin_k for item in result.stations)
        self.assertAlmostEqual(hottest, limit - result.maximum_hot_wall_temperature_k, delta=1.0e-9)


class RefinementTests(unittest.TestCase):
    def test_heat_load_converges_with_station_refinement(self) -> None:
        loads = [_solve(stations).total_heat_load_w for stations in (41, 81, 161)]
        coarse, fine = abs(loads[1] - loads[0]), abs(loads[2] - loads[1])
        self.assertLess(fine, coarse)
        self.assertLess(fine / loads[2], 0.02)



# Hydrogen design: 4 mm channels at a 20 MPa inlet keep the coolant Mach number near 0.4 (with
# the momentum pressure drop, 3 mm channels at 15 MPa choke), and its short channels exercise
# Taylor's out-of-range reporting for the entrance distance.
HYDROGEN_CHANNELS = CoolingChannels(count=100, width_m=1.5e-3, height_m=4.0e-3, hot_wall_thickness_m=1.0e-3,
                                    roughness_m=1.0e-6)


class HydrogenCoolantTests(unittest.TestCase):
    def test_taylor_is_the_default_for_hydrogen_and_runs_hotter_than_uncorrected_gnielinski(self) -> None:
        gas = hot_gas_from_propellant(propellant_pair="LOX_LH2", chamber_pressure_pa=PRESSURE,
                                      throat_radius_of_curvature_m=0.05)
        contour = generate_nozzle_contour(contour="conical", throat_area_m2=0.01, area_ratio=10.0,
                                          gamma=gas.gamma, chamber_radius_ratio=2.0, station_count=41)
        coolant = Coolant(fluid="parahydrogen", inlet_temperature_k=40.0, inlet_pressure_pa=2.0e7,
                          mass_flow_kg_s=PRESSURE * 0.01 / gas.characteristic_velocity_m_s / 6.55157)
        auto = solve_regen_chamber(contour=contour, gas=gas, channels=HYDROGEN_CHANNELS, coolant=coolant)
        taylor = solve_regen_chamber(contour=contour, gas=gas, channels=HYDROGEN_CHANNELS, coolant=coolant,
                                     coolant_correlation="taylor")
        gnielinski = solve_regen_chamber(contour=contour, gas=gas, channels=HYDROGEN_CHANNELS, coolant=coolant,
                                         coolant_correlation="gnielinski")
        self.assertEqual(auto.maximum_hot_wall_temperature_k, taylor.maximum_hot_wall_temperature_k)
        # Taylor's (Tw/Tb) exponent lowers h for strongly heated supercritical hydrogen.
        self.assertGreater(taylor.maximum_hot_wall_temperature_k, gnielinski.maximum_hot_wall_temperature_k)
        self.assertLess(abs(taylor.energy_balance_residual_w), 1.0e-9 * taylor.total_heat_load_w)
        # Taylor's T_w is the coolant-side wall: evaluated there, this design stays inside the
        # correlation's temperature-ratio range (1.1-27.6), so only out-of-range entrance
        # distances are reported, once per parameter.
        ratios = [item.cold_wall_temperature_k / item.coolant_temperature_k for item in taylor.stations]
        self.assertLess(max(ratios), 27.6)
        extrapolated = [item for item in taylor.warnings if item.code == "coolant_correlation_extrapolated"]
        self.assertTrue(extrapolated)
        self.assertFalse(any("wall_to_bulk_temperature_ratio" in item.message for item in extrapolated))
        self.assertEqual(len(extrapolated), len({item.message.split()[0] for item in extrapolated}))

class InputTests(unittest.TestCase):
    def test_invalid_designs_are_refused(self) -> None:
        coolant = Coolant(fluid="methane", inlet_temperature_k=120.0, inlet_pressure_pa=1.5e7,
                          mass_flow_kg_s=_fuel_flow())
        with self.assertRaises(InputError):
            solve_regen_chamber(contour=_contour(), gas=_gas(), channels=CHANNELS, coolant=coolant,
                                coolant_correlation="taylor")
        crowded = CoolingChannels(count=400, width_m=1.5e-3, height_m=3.0e-3, hot_wall_thickness_m=1.0e-3)
        with self.assertRaises(DomainError):
            solve_regen_chamber(contour=_contour(), gas=_gas(), channels=crowded, coolant=coolant)


if __name__ == "__main__":
    unittest.main()
