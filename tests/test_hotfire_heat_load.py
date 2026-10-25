"""Integrated gas-side heat load against the coolant enthalpy rise of NASA TN D-3532 firings.

The 500-lbf dump-cooled LOX/GH2 engine of TN D-3532 was cooled by liquid
hydrogen whose inlet and outlet temperature and pressure were measured
(Table II). The measured heat load is ``W_c [h(T_out, p_out) - h(T_in, p_in)]``
with the package's parahydrogen equation of state, the same bookkeeping as
the authors' Eq. (1) with enthalpy in place of ``c_p dT``.

The prediction integrates ``h_Bartz (T_aw - T_w)`` over the Fig. 3 wall
(1.881 in. chamber radius, 4.352 in. arcs, 15 deg exit cone). The gas state is
the CEA-equivalent equilibrium for GH2 (298.15 K) and LOX at the firing's
O/F and 114.7 psia, scaled to the authors' "calculated combustion chamber
total temperature" (Table I), i.e. ``T_c`` and ``c*`` with
``eta_c* = sqrt(T_table / T_ideal)``.

What it shows (a bracketing check, not an accuracy target): the wall
temperature, which the report gives only as Fig. 15 curves between about 900
and 1600 K, decides the answer. Uniform walls of 1100 K and 1500 K bracket the
measured heat load for firings 8-10, from about 1.36 to 0.92 times it; the
ideal (eta = 1) gas over-predicts by 1.5-1.9 times at 700-1000 K. The
convergence half-angle (not printed; 15-24.7 deg fit the drawing) changes the
result by under 1 %.
"""

from __future__ import annotations

import json
import unittest
from dataclasses import replace
from math import cos, pi, radians, sin, sqrt, tan
from pathlib import Path

from rocket_propulsion.compressible.isentropic import mach_from_isentropic
from rocket_propulsion.heat_transfer.regen_chamber import bartz_coefficient, hot_gas_from_reactants
from rocket_propulsion.thermochemistry.species import get_reactant
from rocket_propulsion.thermodynamics.real_fluids import reference_fluid

DATA = json.loads((Path(__file__).parent / "reference_data" / "hotfire_tn_d3532.json").read_text(encoding="utf-8"))
IN = 0.0254
PSI = 6894.757293168361
LB = 0.45359237
RANKINE = 5.0 / 9.0
G = DATA["geometry_in"]
R_ARC, RC, RT, RE = G["wall_arc_radius"] * IN, G["chamber_diameter"] / 2 * IN, G["throat_diameter"] / 2 * IN, G["exit_diameter"] / 2 * IN
XT, XE = G["injector_to_throat"] * IN, G["injector_to_exit"] * IN


def wall_radius(theta_deg: float):
    th = radians(theta_deg)
    dx_arc, dr_arc = R_ARC * sin(th), R_ARC * (1 - cos(th))
    dx_cone = ((RC - RT) - 2 * dr_arc) / tan(th)
    x0 = XT - 2 * dx_arc - dx_cone
    th_e = radians(G["exit_half_angle_deg"])
    dx_e, dr_e = R_ARC * sin(th_e), R_ARC * (1 - cos(th_e))

    def r(x: float) -> float:
        if x <= x0:
            return RC
        if x <= x0 + dx_arc:
            u = x - x0
            return RC - (R_ARC - sqrt(R_ARC ** 2 - u * u))
        if x <= x0 + dx_arc + dx_cone:
            return RC - dr_arc - (x - x0 - dx_arc) * tan(th)
        if x <= XT + dx_e:
            u = abs(x - XT)
            return RT + (R_ARC - sqrt(R_ARC ** 2 - u * u))
        return RT + dr_e + (x - XT - dx_e) * tan(th_e)
    return r


def predicted_load(gas, theta_deg: float, wall_k: float, n: int = 1500) -> float:
    r = wall_radius(theta_deg)
    xs = [XE * k / n for k in range(n + 1)]
    rs = [r(x) for x in xs]
    total = 0.0
    for i in range(n):
        x, radius = 0.5 * (xs[i] + xs[i + 1]), 0.5 * (rs[i] + rs[i + 1])
        ds = sqrt((xs[i + 1] - xs[i]) ** 2 + (rs[i + 1] - rs[i]) ** 2)
        ratio = (radius / RT) ** 2
        branch = "area_ratio_subsonic" if x < XT else "area_ratio_supersonic"
        mach = 1.0 if abs(ratio - 1.0) < 1e-9 else mach_from_isentropic(branch, ratio, gas.gamma)
        h, taw = bartz_coefficient(gas, throat_radius_m=RT, area_ratio=ratio, mach=mach, wall_temperature_k=wall_k)
        total += h * (taw - wall_k) * 2 * pi * radius * ds
    return total


def measured_load(firing: dict, fluid: str = "parahydrogen") -> float:
    f = reference_fluid(fluid)
    t_in = sum(firing["inlet_temperature_r"]) / 2 * RANKINE
    p_in = (sum(firing["inlet_pressure_psig"]) / 2 + 14.7) * PSI
    t_out = sum(firing["outlet_temperature_r"]) / 2 * RANKINE
    p_out = (sum(firing["outlet_pressure_psig"]) / 2 + 14.7) * PSI
    return firing["coolant_flow_lb_s"] * LB * (f.state_tp(t_out, p_out).enthalpy_j_kg - f.state_tp(t_in, p_in).enthalpy_j_kg)


def gas_state(firing: dict, *, efficiency_corrected: bool):
    h2, o2 = get_reactant("H2"), get_reactant("O2(L)")
    ratio = firing["oxidizer_fuel_ratio"]
    gas = hot_gas_from_reactants(
        reactants={h2.key: 1.0, o2.key: ratio * h2.molecular_mass_kg_mol / o2.molecular_mass_kg_mol},
        chamber_pressure_pa=(DATA["nominal_chamber_pressure_psig"] + 14.7) * PSI, throat_radius_of_curvature_m=R_ARC)
    if not efficiency_corrected:
        return gas, 1.0
    t_table = firing["calculated_chamber_temperature_r"] * RANKINE
    eta = sqrt(t_table / gas.chamber_temperature_k)
    return replace(gas, chamber_temperature_k=t_table,
                   characteristic_velocity_m_s=gas.characteristic_velocity_m_s * eta), eta


class TnD3532HeatLoad(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.cases = {}
        for key, firing in DATA["firings"].items():
            corrected, eta = gas_state(firing, efficiency_corrected=True)
            ideal, _ = gas_state(firing, efficiency_corrected=False)
            cls.cases[key] = {"measured": measured_load(firing), "eta": eta, "corrected": corrected, "ideal": ideal}

    def test_contour_matches_the_drawing(self) -> None:
        for theta in (15.0, 24.5):
            r = wall_radius(theta)
            self.assertAlmostEqual(r(0.0), RC, places=12)
            self.assertAlmostEqual(r(XT), RT, places=12)
            self.assertLess(abs(r(XE) - RE), 2e-4)  # 2.94 in. of divergence against the printed 10.94 - 8

    def test_measured_heat_load_and_spin_state(self) -> None:
        for key, case in self.cases.items():
            with self.subTest(firing=key):
                self.assertGreater(case["measured"], 2.0e5)
                self.assertLess(case["measured"], 2.7e5)
                normal = measured_load(DATA["firings"][key], "normal_hydrogen")
                self.assertLess(abs(normal / case["measured"] - 0.87), 0.02)

    def test_efficiency_implied_by_table_i(self) -> None:
        for case in self.cases.values():
            self.assertGreater(case["eta"], 0.95)
            self.assertLess(case["eta"], 0.98)

    def test_wall_temperature_band_brackets_the_measurement(self) -> None:
        for key, case in self.cases.items():
            cool_wall = predicted_load(case["corrected"], 20.0, 1100.0) / case["measured"]
            hot_wall = predicted_load(case["corrected"], 20.0, 1500.0) / case["measured"]
            with self.subTest(firing=key):
                self.assertGreater(cool_wall, 1.0)
                self.assertLess(cool_wall, 1.40)
                self.assertLess(hot_wall, 1.0)
                self.assertGreater(hot_wall, 0.90)

    def test_ideal_gas_state_overpredicts(self) -> None:
        for case in self.cases.values():
            ratio = predicted_load(case["ideal"], 20.0, 1000.0) / case["measured"]
            self.assertGreater(ratio, 1.45)

    def test_convergence_angle_is_immaterial(self) -> None:
        case = self.cases["9"]
        a = predicted_load(case["corrected"], 15.0, 1300.0)
        b = predicted_load(case["corrected"], 24.5, 1300.0)
        self.assertLess(abs(a / b - 1.0), 0.01)


if __name__ == "__main__":
    unittest.main()
