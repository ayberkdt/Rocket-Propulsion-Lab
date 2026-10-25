"""Radiation exchange against Monte Carlo ray tracing and independent linear algebra.

* View factors are checked by tracing cosine-weighted (diffuse) rays from a
  uniformly sampled emitter and counting hits, the definition of ``F1-2``.
  The frustum and the coaxial disks are convex enclosures, so a ray from one
  opening reaches the other opening exactly when its end point at that plane
  lies inside the disk. Tolerances are five binomial standard deviations of
  the estimate, ``5 sqrt(F (1 - F) / N)``, with seeded generators.
* The net-radiation method is re-solved here with numpy from the radiosity
  equations ``B_i - (1 - e_i) sum_j F_ij B_j = e_i sigma T_i^4`` (prescribed
  temperature) and ``B_k = sum_j F_kj B_j`` (adiabatic wall).
* The radiative-equilibrium wall temperature is re-solved with ``brentq``;
  the parallel-plate shield stack has the closed form
  ``q = sigma (T1^4 - T2^4) / ((N + 1)(2/e - 1))`` and, for one shield,
  ``Ts^4 = (T1^4 + T2^4)/2``.
"""

from __future__ import annotations

import itertools
import math
import unittest

import numpy as np
from scipy import optimize

from rocket_propulsion.core.constants import STEFAN_BOLTZMANN_W_M2_K4
from rocket_propulsion.core.errors import DomainError
from rocket_propulsion.heat_transfer import radiation as rad

SIGMA = STEFAN_BOLTZMANN_W_M2_K4
RAYS = 2_000_000


def diffuse_directions(rng: np.random.Generator, n: int) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Unit vectors with a cosine (Lambertian) distribution about +z."""
    u, v = rng.random(n), rng.random(n)
    sin_theta, phi = np.sqrt(u), 2.0 * math.pi * v
    return sin_theta * np.cos(phi), sin_theta * np.sin(phi), np.sqrt(1.0 - u)


def disk_points(rng: np.random.Generator, n: int, radius: float) -> tuple[np.ndarray, np.ndarray]:
    r, phi = radius * np.sqrt(rng.random(n)), 2.0 * math.pi * rng.random(n)
    return r * np.cos(phi), r * np.sin(phi)


def disk_to_coaxial_disk(radius_1: float, radius_2: float, separation: float, seed: int) -> float:
    rng = np.random.default_rng(seed)
    x, y = disk_points(rng, RAYS, radius_1)
    dx, dy, dz = diffuse_directions(rng, RAYS)
    t = separation / dz
    return float(np.mean((x + t * dx) ** 2 + (y + t * dy) ** 2 <= radius_2 * radius_2))


def tolerance(f: float) -> float:
    return 5.0 * math.sqrt(max(f * (1.0 - f), 1e-12) / RAYS)


class ViewFactorsByRayTracing(unittest.TestCase):
    def test_coaxial_disks(self) -> None:
        for seed, (r1, r2, h) in enumerate(((1.0, 1.0, 1.0), (0.3, 1.2, 0.5), (1.5, 0.4, 2.0), (1.0, 2.0, 0.05))):
            expected = rad.view_factor_coaxial_disks(r1, r2, h)
            traced = disk_to_coaxial_disk(r1, r2, h, seed)
            with self.subTest(r1=r1, r2=r2, h=h):
                self.assertLess(abs(traced - expected), tolerance(expected))

    def test_frustum_openings_and_wall(self) -> None:
        for seed, (r1, r2, length) in enumerate(((0.1, 0.3, 0.5), (0.2, 0.2, 0.4), (0.05, 0.6, 1.2)), start=10):
            enclosure = rad.frustum_enclosure_view_factors(r1, r2, length)
            small_to_large = disk_to_coaxial_disk(r1, r2, length, seed)
            f = enclosure.view_factors
            with self.subTest(r1=r1, r2=r2, length=length):
                self.assertLess(abs(f[0][1] - small_to_large), tolerance(f[0][1]))
                # the rest of the small opening's rays strike the wall (convex enclosure)
                self.assertLess(abs(f[0][2] - (1.0 - small_to_large)), tolerance(f[0][2]))
                self.assertEqual(f[0][0], 0.0)
                areas = enclosure.areas_m2
                self.assertAlmostEqual(areas[0], math.pi * r1 * r1, places=14)
                self.assertAlmostEqual(areas[2], math.pi * (r1 + r2) * math.hypot(r2 - r1, length), places=12)
                for i in range(3):
                    self.assertAlmostEqual(sum(f[i]), 1.0, places=12)
                    for j in range(3):
                        self.assertAlmostEqual(areas[i] * f[i][j], areas[j] * f[j][i], places=12)

    def test_parallel_rectangles(self) -> None:
        rng = np.random.default_rng(21)
        for a, b, c in ((1.0, 1.0, 1.0), (2.0, 0.5, 0.3), (0.4, 3.0, 2.5)):
            x, y = a * rng.random(RAYS), b * rng.random(RAYS)
            dx, dy, dz = diffuse_directions(rng, RAYS)
            t = c / dz
            xh, yh = x + t * dx, y + t * dy
            traced = float(np.mean((xh >= 0) & (xh <= a) & (yh >= 0) & (yh <= b)))
            expected = rad.view_factor_parallel_rectangles(a, b, c)
            with self.subTest(a=a, b=b, c=c):
                self.assertLess(abs(traced - expected), tolerance(expected))

    def test_perpendicular_rectangles_sharing_an_edge(self) -> None:
        # rectangle 1 in z = 0 with x in (0, w), y in (0, l); rectangle 2 in x = 0 with z in (0, h)
        rng = np.random.default_rng(31)
        for w, h, edge in ((1.0, 1.0, 1.0), (0.5, 2.0, 1.5), (2.0, 0.3, 0.8)):
            x, y = w * rng.random(RAYS), edge * rng.random(RAYS)
            dx, dy, dz = diffuse_directions(rng, RAYS)
            towards = dx < 0.0
            t = np.where(towards, -x / np.where(towards, dx, -1.0), np.inf)
            zh, yh = t * dz, y + t * dy
            traced = float(np.mean(towards & (zh <= h) & (yh >= 0) & (yh <= edge)))
            expected = rad.view_factor_perpendicular_rectangles(w, h, edge)
            with self.subTest(w=w, h=h, edge=edge):
                self.assertLess(abs(traced - expected), tolerance(expected))


class NetRadiationMethod(unittest.TestCase):
    def setUp(self) -> None:
        self.enclosure = rad.frustum_enclosure_view_factors(0.1, 0.35, 0.6)
        self.f = np.array(self.enclosure.view_factors)
        self.areas = np.array(self.enclosure.areas_m2)

    def test_all_temperatures_prescribed(self) -> None:
        emissivity = np.array([1.0, 0.6, 0.35])
        temperatures = np.array([2400.0, 300.0, 1300.0])
        surfaces = [rad.GraySurface(name, float(a), float(e), temperature_k=float(t))
                    for name, a, e, t in zip(("small", "large", "wall"), self.areas, emissivity, temperatures, strict=True)]
        solution = rad.solve_gray_enclosure(surfaces, self.enclosure.view_factors)
        matrix = np.eye(3) - (1.0 - emissivity)[:, None] * self.f
        radiosity = np.linalg.solve(matrix, emissivity * SIGMA * temperatures ** 4)
        net = radiosity - self.f @ radiosity
        np.testing.assert_allclose(solution.radiosities_w_m2, radiosity, rtol=1e-12)
        np.testing.assert_allclose(solution.net_heat_fluxes_w_m2, net, rtol=1e-10, atol=1e-9 * np.abs(net).max())
        self.assertLess(abs(sum(solution.net_heat_rates_w)), 1e-9 * float(np.abs(self.areas * net).max()))
        # the hottest surface loses heat, the coldest gains it
        self.assertGreater(solution.net_heat_rates_w[0], 0.0)
        self.assertLess(solution.net_heat_rates_w[1], 0.0)

    def test_adiabatic_wall(self) -> None:
        e_small, e_large, e_wall = 0.9, 0.5, 0.2
        t_small, t_large = 2000.0, 400.0
        surfaces = [rad.GraySurface("small", float(self.areas[0]), e_small, temperature_k=t_small),
                    rad.GraySurface("large", float(self.areas[1]), e_large, temperature_k=t_large),
                    rad.GraySurface("wall", float(self.areas[2]), e_wall, net_heat_flux_w_m2=0.0)]
        solution = rad.solve_gray_enclosure(surfaces, self.enclosure.view_factors)
        matrix = np.zeros((3, 3))
        rhs = np.zeros(3)
        for i, (e, t) in enumerate(((e_small, t_small), (e_large, t_large))):
            matrix[i] = -(1.0 - e) * self.f[i]
            matrix[i, i] += 1.0
            rhs[i] = e * SIGMA * t ** 4
        matrix[2] = -self.f[2]
        matrix[2, 2] += 1.0
        radiosity = np.linalg.solve(matrix, rhs)
        wall_temperature = (radiosity[2] / SIGMA) ** 0.25
        self.assertLess(abs(solution.temperatures_k[2] - wall_temperature), 1e-9 * wall_temperature)
        self.assertLess(abs(solution.net_heat_fluxes_w_m2[2]), 1e-9 * SIGMA * t_small ** 4)
        # a refractory wall does not depend on its own emittance
        other = rad.solve_gray_enclosure([*surfaces[:2], rad.GraySurface("wall", float(self.areas[2]), 0.8,
                                                                         net_heat_flux_w_m2=0.0)],
                                         self.enclosure.view_factors)
        self.assertLess(abs(other.temperatures_k[2] - solution.temperatures_k[2]), 1e-9 * wall_temperature)
        self.assertLess(abs(other.net_heat_rates_w[0] - solution.net_heat_rates_w[0]),
                        1e-9 * abs(solution.net_heat_rates_w[0]))

    def test_view_factor_errors_are_refused(self) -> None:
        broken = [list(row) for row in self.enclosure.view_factors]
        broken[0][1] += 1e-3
        surfaces = [rad.GraySurface(n, float(a), 0.5, temperature_k=1000.0) for n, a in zip("abc", self.areas, strict=True)]
        with self.assertRaises(DomainError):
            rad.solve_gray_enclosure(surfaces, broken)


class EquilibriumAndShields(unittest.TestCase):
    def test_radiative_equilibrium_matches_brentq(self) -> None:
        cases = (((0.85, 300.0),), ((0.6, 3.0), (0.25, 1500.0)), ((0.9, 250.0), (0.05, 2800.0)))
        for sinks in cases:
            for h, recovery, imposed in ((150.0, 2500.0, 0.0), (0.0, None, 4.0e5), (800.0, 3200.0, 2.0e4)):
                def balance(t: float, sinks=sinks, h=h, recovery=recovery, imposed=imposed) -> float:
                    convective = h * (recovery - t) if recovery is not None else 0.0
                    return imposed + convective - SIGMA * sum(f * (t ** 4 - ts ** 4) for f, ts in sinks)
                root = optimize.brentq(balance, 1.0, 1.0e4, xtol=1e-12, rtol=1e-15)
                result = rad.radiative_equilibrium_wall_temperature(
                    sinks=[rad.RadiativeSink(f, t) for f, t in sinks], heat_transfer_coefficient_w_m2_k=h,
                    recovery_temperature_k=recovery, imposed_heat_flux_w_m2=imposed)
                with self.subTest(sinks=sinks, h=h, recovery=recovery, imposed=imposed):
                    self.assertLess(abs(result.wall_temperature_k - root), 1e-9 * root)
                    self.assertLess(abs(result.residual_w_m2), 1e-6 * max(1.0, abs(result.radiated_heat_flux_w_m2)))

    def test_parallel_plate_shield_stack(self) -> None:
        t1, t2 = 1200.0, 300.0
        for emissivity in (0.05, 0.3, 0.9):
            for count in range(5):
                shields = [rad.RadiationShield(1.0, emissivity, emissivity)] * count
                result = rad.shielded_heat_rate(inner_temperature_k=t1, outer_temperature_k=t2, inner_area_m2=1.0,
                                                inner_emissivity=emissivity, outer_area_m2=1.0,
                                                outer_emissivity=emissivity, shields=shields)
                expected = SIGMA * (t1 ** 4 - t2 ** 4) / ((count + 1) * (2.0 / emissivity - 1.0))
                with self.subTest(emissivity=emissivity, shields=count):
                    self.assertLess(abs(result.heat_rate_w - expected), 1e-12 * expected)
                    # equal resistances: the shield T^4 values are equally spaced
                    fourth = [t1 ** 4, *(t ** 4 for t in result.shield_temperatures_k), t2 ** 4]
                    steps = [b - a for a, b in itertools.pairwise(fourth)]
                    for step in steps:
                        self.assertLess(abs(step - steps[0]), 1e-9 * abs(steps[0]))


if __name__ == "__main__":
    unittest.main()
