"""Method-of-characteristics water hammer: reference, limit, conservation and errors.

Oracles and sources
-------------------
* Joukowsky ``dP = -rho a dV`` and wave speed
  ``a = sqrt(beta/rho) / sqrt(1 + beta D C / (E t))`` with the boundary factors
  of Figure 5: NASA/TM-20210022275, printed pp. 19-20,
  https://ntrs.nasa.gov/citations/20210022275
* Case study 8-1 of the same report (printed p. 22): 0.088 lb/s through a
  0.62 in^2 pipe, a = 39,370 in/s, rho = 0.036 lb/in^3, "instantaneous"
  closure; Joukowsky surge 14.5 psi, recovered by the method of
  characteristics.  Evaluated here in lbm-in-s units with
  g_c = 32.174 lbm ft/(lbf s^2) (Majumdar and Steadman, AIAA 99-0879 draft,
  nomenclature, https://ntrs.nasa.gov/citations/19980236009), so that no SI
  conversion constant is needed.
* Closures faster than the wave round trip give the Joukowsky surge; slower
  closures give smaller surges (same report, printed pp. 13 and 24).
"""

from __future__ import annotations

import itertools
import time
import unittest
from dataclasses import replace
from math import pi, sqrt

from rocket_propulsion.core.errors import DomainError
from rocket_propulsion.propulsion.burns.water_hammer import (
    LinearValveLaw,
    PipeAnchoring,
    WaterHammerConfiguration,
    WaterHammerPipe,
    joukowsky_surge_pa,
    pipe_boundary_factor,
    pipe_wave_speed_m_s,
    simulate_water_hammer,
)

GC_LBM_IN_PER_LBF_S2 = 32.174 * 12.0


def si_line(*, friction: float = 0.02, reaches: int = 20, closure_s: float = 0.0,
            duration_s: float = 0.5, acceleration: float = 0.0) -> WaterHammerConfiguration:
    """A 100 m, 20 mm line carrying a dense storable-like liquid (synthetic)."""

    pipe = WaterHammerPipe(length_m=100.0, inner_diameter_m=0.02, wave_speed_m_s=1200.0,
                           liquid_density_kg_m3=1140.0, darcy_friction_factor=friction,
                           axial_acceleration_m_s2=acceleration)
    return WaterHammerConfiguration(
        pipe=pipe, reservoir_pressure_pa=2.0e6, downstream_pressure_pa=1.8e6,
        valve_full_open_resistance_pa_s2_m6=1.0e12, valve_opening=LinearValveLaw(0.0, closure_s),
        reach_count=reaches, duration_s=duration_s)


class WaveSpeedAndJoukowskyTests(unittest.TestCase):
    def test_boundary_factors_match_figure_5(self) -> None:
        mu = 0.3
        thin = {"poisson_ratio": mu, "inner_diameter_m": 0.1, "wall_thickness_m": 0.001}
        self.assertAlmostEqual(pipe_boundary_factor(PipeAnchoring.UPSTREAM_ANCHORED, **thin),
                               1.0 - mu / 2.0)
        self.assertAlmostEqual(
            pipe_boundary_factor(PipeAnchoring.NO_LONGITUDINAL_MOVEMENT, **thin), 1.0 - mu * mu)
        self.assertEqual(pipe_boundary_factor(PipeAnchoring.EXPANSION_JOINTS, **thin), 1.0)
        self.assertEqual(pipe_boundary_factor(PipeAnchoring.RIGID, **thin), 0.0)
        d, t = 0.02, 0.002  # D/t = 10 <= 25: thick-walled forms
        thick = {"poisson_ratio": mu, "inner_diameter_m": d, "wall_thickness_m": t}
        base = 2.0 * t / d * (1.0 + mu)
        self.assertAlmostEqual(pipe_boundary_factor(PipeAnchoring.UPSTREAM_ANCHORED, **thick),
                               base + d / (d + t) * (1.0 - mu / 2.0))
        self.assertAlmostEqual(pipe_boundary_factor(PipeAnchoring.EXPANSION_JOINTS, **thick),
                               base + d / (d + t))
        self.assertEqual(pipe_boundary_factor(PipeAnchoring.RIGID, **thick), 0.0)

    def test_rigid_pipe_wave_speed_is_liquid_sound_speed_and_elasticity_slows_it(self) -> None:
        common = {"bulk_modulus_pa": 2.2e9, "density_kg_m3": 1000.0, "inner_diameter_m": 0.1,
                      "wall_thickness_m": 0.002, "youngs_modulus_pa": 7.0e10, "poisson_ratio": 0.33}
        rigid = pipe_wave_speed_m_s(anchoring=PipeAnchoring.RIGID, **common)
        self.assertAlmostEqual(rigid, sqrt(2.2e9 / 1000.0), places=9)
        elastic = pipe_wave_speed_m_s(anchoring=PipeAnchoring.EXPANSION_JOINTS, **common)
        expected = sqrt(2.2e9 / 1000.0) / sqrt(1.0 + 2.2e9 * 0.1 / (7.0e10 * 0.002))
        self.assertAlmostEqual(elastic, expected, places=9)
        self.assertLess(elastic, rigid)

    def test_case_study_joukowsky_value_is_reproduced(self) -> None:
        rho, area, mdot = 0.036, 0.62, 0.088  # lbm/in^3, in^2, lbm/s
        velocity = mdot / (rho * area)
        surge_psi = joukowsky_surge_pa(density_kg_m3=rho, wave_speed_m_s=39_370.0,
                                       velocity_change_m_s=-velocity) / GC_LBM_IN_PER_LBF_S2
        self.assertAlmostEqual(surge_psi, 14.5, delta=0.05)

    def test_moc_recovers_case_study_joukowsky_for_instantaneous_closure(self) -> None:
        rho, area, mdot = 0.036, 0.62, 0.088
        q = mdot / rho
        p_res, p_down = 508.0 * GC_LBM_IN_PER_LBF_S2, 500.0 * GC_LBM_IN_PER_LBF_S2
        pipe = WaterHammerPipe(length_m=39.4, inner_diameter_m=sqrt(4.0 * area / pi),
                               wave_speed_m_s=39_370.0, liquid_density_kg_m3=rho,
                               darcy_friction_factor=0.0)
        config = WaterHammerConfiguration(
            pipe=pipe, reservoir_pressure_pa=p_res, downstream_pressure_pa=p_down,
            valve_full_open_resistance_pa_s2_m6=(p_res - p_down) / q ** 2,
            valve_opening=LinearValveLaw(0.0, 0.0), reach_count=40, duration_s=0.004)
        result = simulate_water_hammer(config)
        self.assertAlmostEqual(result.initial_volume_flow_m3_s, q, places=9)
        self.assertAlmostEqual(result.valve_peak_surge_pa / GC_LBM_IN_PER_LBF_S2, 14.5, delta=0.05)
        self.assertAlmostEqual(result.surge_to_joukowsky_ratio, 1.0, places=12)


class MethodOfCharacteristicsTests(unittest.TestCase):
    def test_unchanged_valve_keeps_the_discrete_steady_state(self) -> None:
        config = replace(si_line(friction=0.03, acceleration=0.3),
                         valve_opening=LinearValveLaw(0.0, 0.0, 1.0, 1.0))
        result = simulate_water_hammer(config)
        for series in result.station_pressure_pa + result.station_volume_flow_m3_s:
            self.assertLess(max(series) - min(series), 1e-9 * max(abs(v) for v in series) + 1e-15)
        # Steady drop equals friction plus body force: rho f L V^2/(2D) + rho L alpha.
        pipe = config.pipe
        v = result.initial_velocity_m_s
        drop = result.station_pressure_pa[0][0] - result.valve_initial_pressure_pa
        self.assertGreater(v, 0.0)
        expected = (pipe.liquid_density_kg_m3 * pipe.darcy_friction_factor * pipe.length_m * v * v
                    / (2.0 * pipe.inner_diameter_m)
                    + pipe.liquid_density_kg_m3 * pipe.length_m * pipe.axial_acceleration_m_s2)
        self.assertAlmostEqual(drop, expected, delta=1e-9 * expected)

    def test_closure_faster_than_round_trip_converges_to_joukowsky(self) -> None:
        round_trip = 2.0 * 100.0 / 1200.0
        ratios = []
        for fraction in (4.0, 2.0, 1.0, 0.5, 0.0):
            result = simulate_water_hammer(si_line(friction=0.0, reaches=40,
                                                   closure_s=fraction * round_trip,
                                                   duration_s=5.0 * round_trip))
            ratios.append(result.surge_to_joukowsky_ratio)
        self.assertTrue(all(b >= a - 1e-12 for a, b in itertools.pairwise(ratios)),
                        ratios)
        self.assertLess(ratios[0], 0.5)
        for ratio in ratios[2:]:
            self.assertAlmostEqual(ratio, 1.0, places=10)

    def test_friction_line_packing_raises_peak_above_joukowsky_and_converges(self) -> None:
        values = [simulate_water_hammer(si_line(reaches=n, duration_s=0.3)).valve_peak_surge_pa
                  for n in (10, 20, 40, 80)]
        joukowsky = simulate_water_hammer(si_line(reaches=10, duration_s=0.01)).joukowsky_surge_pa
        self.assertTrue(all(v > joukowsky for v in values))
        differences = [abs(b - a) for a, b in itertools.pairwise(values)]
        self.assertLess(differences[-1], 0.6 * differences[0])

    def test_first_step_after_instantaneous_closure_is_exactly_joukowsky(self) -> None:
        result = simulate_water_hammer(si_line(reaches=20, duration_s=0.01))
        valve = result.station_pressure_pa[-1]
        self.assertAlmostEqual(valve[1] - valve[0], result.joukowsky_surge_pa,
                               delta=1e-9 * result.joukowsky_surge_pa)

    def test_wave_arrives_at_mid_pipe_after_half_the_travel_time(self) -> None:
        config = si_line(friction=0.0, reaches=20, duration_s=0.1)
        result = simulate_water_hammer(config)
        middle = result.station_pressure_pa[1]
        arrival = next(i for i, p in enumerate(middle) if p - middle[0] > 1.0)
        # The closure acts at the first step, t = dt, then travels 50 m.
        self.assertAlmostEqual(result.time_s[arrival], result.time_step_s + 50.0 / 1200.0,
                               places=12)

    def test_mass_and_energy_balances_close_and_converge(self) -> None:
        mass, energy = [], []
        for n in (10, 20, 40):
            result = simulate_water_hammer(si_line(reaches=n, closure_s=0.05, duration_s=0.5))
            stored_scale = (result.joukowsky_surge_pa * 0.25 * pi * 0.02 ** 2 * 100.0
                            / 1200.0 ** 2)
            mass.append(result.maximum_mass_balance_residual_kg / stored_scale)
            energy.append(abs(result.energy_balance_residual_j) / result.initial_acoustic_energy_j)
            self.assertGreater(result.friction_dissipation_j, 0.0)
        self.assertLess(mass[-1], 0.01)
        self.assertLess(energy[-1], 0.01)
        self.assertLess(mass[-1], 0.6 * mass[0])
        self.assertLess(energy[-1], 0.6 * energy[0])

    def test_lossless_closed_line_conserves_acoustic_energy_exactly(self) -> None:
        result = simulate_water_hammer(si_line(friction=0.0, reaches=25, duration_s=1.0))
        self.assertAlmostEqual(result.final_acoustic_energy_j, result.initial_acoustic_energy_j,
                               delta=1e-9 * result.initial_acoustic_energy_j)
        self.assertLess(abs(result.energy_balance_residual_j),
                        1e-9 * result.initial_acoustic_energy_j)

    def test_valve_opening_from_rest_accelerates_the_line(self) -> None:
        config = replace(si_line(friction=0.02, duration_s=1.0),
                         valve_opening=LinearValveLaw(0.0, 0.05, 0.0, 1.0))
        result = simulate_water_hammer(config)
        self.assertEqual(result.initial_volume_flow_m3_s, 0.0)
        valve_flow = result.station_volume_flow_m3_s[-1]
        self.assertGreater(valve_flow[-1], 0.0)
        self.assertLess(min(result.station_pressure_pa[-1]), 2.0e6)

    def test_runs_fast(self) -> None:
        start = time.perf_counter()
        simulate_water_hammer(si_line(reaches=50, duration_s=1.0))
        self.assertLess(time.perf_counter() - start, 2.0)


class WaterHammerErrorTests(unittest.TestCase):
    def test_invalid_inputs_fail_closed(self) -> None:
        base = si_line()
        for kwargs in (
            {"reach_count": 0}, {"duration_s": 0.0}, {"reservoir_pressure_pa": float("nan")},
            {"downstream_pressure_pa": -1.0}, {"valve_full_open_resistance_pa_s2_m6": 0.0},
            {"valve_opening": 1.0}, {"station_indices": (0, 0)},
            {"station_indices": (99,)}, {"schema": "other"}, {"vapor_pressure_pa": -1.0},
        ):
            with self.subTest(kwargs=kwargs), self.assertRaises(DomainError):
                replace(base, **kwargs)
        for kwargs in ({"length_m": 0.0}, {"darcy_friction_factor": -0.1},
                       {"wave_speed_m_s": float("inf")}, {"liquid_density_kg_m3": True}):
            with self.subTest(kwargs=kwargs), self.assertRaises(DomainError):
                replace(base.pipe, **kwargs)
        with self.assertRaises(DomainError):
            LinearValveLaw(0.0, 0.1, 1.5, 0.0)
        with self.assertRaises(DomainError):
            pipe_boundary_factor("rigid", poisson_ratio=0.3, inner_diameter_m=0.1,
                                 wall_thickness_m=0.01)
        with self.assertRaises(DomainError):
            pipe_boundary_factor(PipeAnchoring.RIGID, poisson_ratio=0.5, inner_diameter_m=0.1,
                                 wall_thickness_m=0.01)
        with self.assertRaises(DomainError):
            joukowsky_surge_pa(density_kg_m3=1.0, wave_speed_m_s=0.0, velocity_change_m_s=1.0)

    def test_opening_law_leaving_unit_interval_is_rejected(self) -> None:
        config = replace(si_line(), valve_opening=lambda t: 1.0 + t)
        with self.assertRaises(DomainError):
            simulate_water_hammer(config)

    def test_column_separation_is_refused(self) -> None:
        # Low reservoir pressure: the reflected expansion falls below vapour pressure.
        config = replace(si_line(friction=0.0, duration_s=0.5), reservoir_pressure_pa=1.3e6,
                         downstream_pressure_pa=1.1e6, vapor_pressure_pa=1.0e5)
        with self.assertRaises(DomainError):
            simulate_water_hammer(config)


if __name__ == "__main__":
    unittest.main()
