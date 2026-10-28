"""Conjugate axial wall conduction in the regenerative march.

The axial fluxes are recomputed here from the reported station states with
an independent implementation of the wall network: metal section per unit
circumference ``k (t + land share h) + k_c c``, conductance
``G = mean(k A') pi (r_a + r_b) / ds`` between the stations that carry march
area, mid-wall temperature ``(T_hot,mean + T_cold)/2`` and
``S_i A_i = sum_j G_ij (T_j - T_i)``. Conservation demands ``sum S_i A_i = 0``
exactly and the gas-side heat to equal the coolant's; physics demands that
axial conduction draw heat out of the throat peak.
"""

from __future__ import annotations

import itertools
import math
import unittest
from functools import cache

from rocket_propulsion.core.errors import DomainError, InputError
from rocket_propulsion.geometry import generate_nozzle_contour
from rocket_propulsion.heat_transfer.regen_chamber import (
    Coolant,
    CoolingChannels,
    hot_gas_from_propellant,
    solve_regen_chamber,
)
from rocket_propulsion.thermodynamics.materials import get_wall_material
from rocket_propulsion.thermodynamics.real_fluids import reference_fluid

PRESSURE = 7.0e6
CHANNELS = CoolingChannels(count=100, width_m=1.5e-3, height_m=5.0e-3, hot_wall_thickness_m=1.0e-3, roughness_m=1e-6)


@cache
def gas():
    return hot_gas_from_propellant(propellant_pair="LOX_CH4", chamber_pressure_pa=PRESSURE,
                                   throat_radius_of_curvature_m=0.05)


@cache
def contour(stations: int = 41):
    return generate_nozzle_contour(contour="conical", throat_area_m2=0.01, area_ratio=10.0, gamma=gas().gamma,
                                   chamber_radius_ratio=2.0, station_count=stations)


def coolant() -> Coolant:
    return Coolant(fluid="methane", inlet_temperature_k=120.0, inlet_pressure_pa=1.5e7,
                   mass_flow_kg_s=PRESSURE * 0.01 / gas().characteristic_velocity_m_s / 4.4)


@cache
def solved(axial: bool, wall_model: str = "fin", stations: int = 41, closeout: float = 0.0):
    return solve_regen_chamber(contour=contour(stations), gas=gas(), channels=CHANNELS, coolant=coolant(),
                               axial_conduction=axial, wall_model=wall_model, closeout_thickness_m=closeout)


def march_areas(result) -> dict[int, float]:
    st = result.stations
    order = list(reversed(range(len(st))))  # counterflow
    return {i: 2 * math.pi * st[i].radius_m * math.hypot(st[order[n + 1]].x_m - st[i].x_m,
                                                         st[order[n + 1]].radius_m - st[i].radius_m)
            for n, i in enumerate(order[:-1])}


def recomputed_sources(result, closeout: float = 0.0) -> dict[int, float]:
    k = get_wall_material("grcop-84").thermal_conductivity.value
    st = result.stations
    areas = march_areas(result)
    nodes = sorted(areas)

    def ka(i: int) -> float:
        s = st[i]
        pitch = 2 * math.pi * (s.radius_m + s.hot_wall_thickness_m) / 100
        return k * (s.hot_wall_thickness_m + (1 - s.channel_width_m / pitch) * s.channel_height_m) + k * closeout

    def mid(i: int) -> float:
        return 0.5 * (st[i].hot_wall_mean_temperature_k + st[i].cold_wall_temperature_k)
    net = dict.fromkeys(nodes, 0.0)
    for a, b in itertools.pairwise(nodes):
        g = 0.5 * (ka(a) + ka(b)) * math.pi * (st[a].radius_m + st[b].radius_m) / math.hypot(
            st[b].x_m - st[a].x_m, st[b].radius_m - st[a].radius_m)
        flow = g * (mid(b) - mid(a))
        net[a] += flow
        net[b] -= flow
    return {i: net[i] / areas[i] for i in nodes}


class ConjugateAxialConduction(unittest.TestCase):
    def check(self, result, closeout: float = 0.0) -> None:
        st = result.stations
        areas = march_areas(result)
        peak_gas = max(s.heat_flux_w_m2 for s in st)
        # exact zero sum: conduction only moves heat
        self.assertLess(abs(sum(st[i].axial_conduction_flux_w_m2 * a for i, a in areas.items())),
                        1e-9 * result.total_heat_load_w)
        # the gas-side heat equals the coolant's
        gas_heat = sum(st[i].heat_flux_w_m2 * a for i, a in areas.items())
        self.assertLess(abs(gas_heat - result.total_heat_load_w), 1e-8 * result.total_heat_load_w)
        self.assertLess(abs(result.energy_balance_residual_w), 1e-9 * result.total_heat_load_w)
        # converged: the fluxes the march used equal those of its own final temperatures
        expected = recomputed_sources(result, closeout)
        for i, value in expected.items():
            with self.subTest(x=st[i].x_m):
                self.assertLess(abs(st[i].axial_conduction_flux_w_m2 - value), 2e-5 * peak_gas)
        self.assertIn("axial_wall_conduction", {w.code for w in result.warnings})
        # segment by segment the coolant takes the gas flux plus the axial share
        fluid = reference_fluid("methane")
        flow = coolant().mass_flow_kg_s
        order = list(reversed(range(len(st))))
        for n, i in enumerate(order[:-1]):
            a_state, b_state = st[i], st[order[n + 1]]
            rise = (fluid.state_tp(b_state.coolant_temperature_k, b_state.coolant_pressure_pa).enthalpy_j_kg
                    - fluid.state_tp(a_state.coolant_temperature_k, a_state.coolant_pressure_pa).enthalpy_j_kg)
            expected = (a_state.heat_flux_w_m2 + a_state.axial_conduction_flux_w_m2) * areas[i] / flow
            with self.subTest(segment=n):
                self.assertLess(abs(rise - expected), 1e-7 * abs(expected) + 1e-3)

    def test_fin_wall(self) -> None:
        local, conjugate = solved(False), solved(True)
        self.check(conjugate)
        self.assertTrue(all(s.axial_conduction_flux_w_m2 == 0.0 for s in local.stations))
        # heat flows out of the throat peak along the wall: the peak falls, its neighbours gain
        self.assertLess(conjugate.maximum_hot_wall_temperature_k, local.maximum_hot_wall_temperature_k)
        hottest = max(range(len(local.stations)), key=lambda i: local.stations[i].hot_wall_temperature_k)
        self.assertLess(conjugate.stations[hottest].axial_conduction_flux_w_m2, 0.0)
        # redistribution only: the load moves by far less than the axial fluxes themselves
        self.assertLess(abs(conjugate.total_heat_load_w / local.total_heat_load_w - 1), 1e-3)

    def test_2d_wall_with_closeout(self) -> None:
        result = solved(True, "2d", 41, 1.0e-3)
        self.assertEqual(result.wall_model, "2d")
        self.check(result, closeout=1.0e-3)

    def test_wall_profile_is_smoother(self) -> None:
        def roughness(result) -> float:
            t = [s.hot_wall_mean_temperature_k for s in result.stations]
            return max(abs(t[i + 1] - 2 * t[i] + t[i - 1]) for i in range(1, len(t) - 1))
        self.assertLess(roughness(solved(True)), roughness(solved(False)))


class Refusals(unittest.TestCase):
    def test_options_are_validated(self) -> None:
        base = {"contour": contour(41), "gas": gas(), "channels": CHANNELS, "coolant": coolant()}
        with self.assertRaises(InputError):
            solve_regen_chamber(**base, closeout_thickness_m=1e-3)  # fin model has no closeout
        with self.assertRaises(InputError):
            solve_regen_chamber(**base, axial_conduction="yes")
        with self.assertRaises(DomainError):
            solve_regen_chamber(**base, wall_model="2d", closeout_thickness_m=-1e-3)
        with self.assertRaises(InputError):
            solve_regen_chamber(**base, wall_material=None, wall_thermal_conductivity_w_m_k=None,
                                axial_conduction=True)


if __name__ == "__main__":
    unittest.main()
