"""One-dimensional rocket performance against the conservation laws it must satisfy.

For every station the solver reports ``p, T, M, h, s, u`` and ``rho``; the
quantities below are recomputed here from those states alone, with the
equations of one-dimensional steady isentropic flow (Sutton and Biblarz,
*Rocket Propulsion Elements*, 9th ed., Ch. 3; NASA RP-1311 Part I, Sect. 6):

* ideal gas ``rho = p M / (R T)``;
* energy ``u = sqrt(2 (h_c - h))`` and isentropic expansion ``s = s_c``;
* continuity ``A/A* = (rho u)_t / (rho u)_e``; the throat carries the largest
  mass flux of every station;
* ``c* = p_c A_t / mdot = p_c / (rho u)_t``;
* ``C_F,vac = (u_e + p_e A_e / mdot) / c*`` and ``Isp_vac = c* C_F,vac / g0``;
* frozen flow keeps the chamber composition; shifting equilibrium releases
  recombination energy, so its vacuum impulse is not lower than frozen.
"""

from __future__ import annotations

import itertools
import math
import unittest
from typing import Any, ClassVar

from rocket_propulsion.core.constants import STANDARD_GRAVITY_M_S2
from rocket_propulsion.thermochemistry.rocket_equilibrium import (
    get_propellant_pair,
    reactants_for_mass_ratio,
)
from rocket_propulsion.thermochemistry.rocket_performance import solve_rocket_performance

R = 8.314510  # NASA/TP-2002-211556 Appendix A; independent model convention.
CASES = (("LOX_LH2", 6.0, 6.894757e6), ("LOX_CH4", 3.4, 1.0e7), ("LOX_LH2", 4.0, 2.0e6))
AREA_RATIOS = (3.0, 12.0, 40.0)


def mass_flux(station) -> float:
    return station.density_kg_m3 * station.velocity_m_s


class PerformanceIdentities(unittest.TestCase):
    results: ClassVar[dict[tuple[str, float, float, str], Any]] = {}

    @classmethod
    def setUpClass(cls) -> None:
        for key, ratio, pressure in CASES:
            reactants = reactants_for_mass_ratio(get_propellant_pair(key), ratio)
            for chemistry in ("equilibrium", "frozen"):
                cls.results[key, ratio, pressure, chemistry] = solve_rocket_performance(
                    reactants=reactants, chamber_pressure_pa=pressure, area_ratios=AREA_RATIOS, chemistry=chemistry)

    def each(self):
        for label, result in self.results.items():
            with self.subTest(case=label):
                yield label, result

    def test_ideal_gas_energy_and_entropy_at_every_station(self) -> None:
        for _, r in self.each():
            chamber = r.chamber
            self.assertEqual(chamber.velocity_m_s, 0.0)
            for station in (r.throat, *(e.station for e in r.exits)):
                rho = station.pressure_pa * station.molecular_mass_kg_mol / (R * station.temperature_k)
                self.assertLess(abs(station.density_kg_m3 - rho), 1e-12 * rho)
                u = math.sqrt(2.0 * (chamber.enthalpy_j_kg - station.enthalpy_j_kg))
                self.assertLess(abs(station.velocity_m_s - u), 1e-9 * u)
                self.assertLess(abs(station.entropy_j_kg_k - chamber.entropy_j_kg_k),
                                1e-7 * abs(chamber.entropy_j_kg_k))

    def test_continuity_fixes_every_area_ratio(self) -> None:
        for _, r in self.each():
            throat_flux = mass_flux(r.throat)
            for e in r.exits:
                self.assertLess(abs(throat_flux / mass_flux(e.station) - e.area_ratio), 1e-6 * e.area_ratio)
                self.assertGreater(throat_flux, mass_flux(e.station))

    def test_characteristic_velocity_thrust_coefficient_and_impulse(self) -> None:
        for _, r in self.each():
            throat_flux = mass_flux(r.throat)
            cstar = r.chamber_pressure_pa / throat_flux
            self.assertLess(abs(r.characteristic_velocity_m_s - cstar), 1e-9 * cstar)
            for e in r.exits:
                # p_e A_e / mdot = p_e (A_e/A_t) / (rho u)_t
                vacuum_velocity = e.station.velocity_m_s + e.station.pressure_pa * e.area_ratio / throat_flux
                self.assertLess(abs(e.effective_exhaust_velocity_m_s - e.station.velocity_m_s), 1e-12 * vacuum_velocity)
                self.assertLess(abs(e.vacuum_effective_exhaust_velocity_m_s - vacuum_velocity), 1e-6 * vacuum_velocity)
                self.assertLess(abs(e.vacuum_thrust_coefficient - vacuum_velocity / cstar), 1e-6)
                self.assertLess(abs(e.vacuum_specific_impulse_s * STANDARD_GRAVITY_M_S2 - vacuum_velocity),
                                1e-9 * vacuum_velocity)

    def test_expansion_is_monotone_in_area_ratio(self) -> None:
        for _, r in self.each():
            exits = sorted(r.exits, key=lambda e: e.area_ratio)
            for a, b in itertools.pairwise(exits):
                self.assertLess(b.station.pressure_pa, a.station.pressure_pa)
                self.assertLess(b.station.temperature_k, a.station.temperature_k)
                self.assertGreater(b.vacuum_specific_impulse_s, a.vacuum_specific_impulse_s)
            self.assertLess(r.throat.pressure_pa, r.chamber.pressure_pa)
            self.assertGreater(r.throat.pressure_pa, exits[0].station.pressure_pa)

    def test_frozen_keeps_composition_and_equilibrium_is_not_worse(self) -> None:
        for key, ratio, pressure in CASES:
            frozen = self.results[key, ratio, pressure, "frozen"]
            shifting = self.results[key, ratio, pressure, "equilibrium"]
            with self.subTest(case=(key, ratio, pressure)):
                chamber = dict(frozen.chamber.mole_fractions)
                for station in (frozen.throat, *(e.station for e in frozen.exits)):
                    # the frozen mixture omits species below 1e-10 mole fraction
                    # (rocket_performance._Expansion), a bounded ~n * 1e-10 drift
                    molar_mass = frozen.chamber.molecular_mass_kg_mol
                    self.assertLess(abs(station.molecular_mass_kg_mol - molar_mass), 1e-9 * molar_mass)
                    for species, fraction in station.mole_fractions:
                        self.assertLess(abs(fraction - chamber[species]), 1e-12)
                self.assertEqual(frozen.chamber, shifting.chamber)
                for f, s in zip(frozen.exits, shifting.exits, strict=True):
                    self.assertGreaterEqual(s.vacuum_specific_impulse_s, f.vacuum_specific_impulse_s)
                    self.assertGreaterEqual(s.station.temperature_k, f.station.temperature_k)


if __name__ == "__main__":
    unittest.main()
