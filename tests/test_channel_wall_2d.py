"""Two-dimensional wall conduction across channel and land, against independent solutions.

* Lumped limit: as ``k -> inf`` the wall is isothermal and the heat balance
  ``h_g (T_aw - T) W/2 + q_rad W/2 = h_c (T - T_b)(w/2 + h)`` fixes ``T``.
* An independent discretization written here (uniform cells, vectorized
  assembly, ``scipy.sparse.linalg.spsolve``) must converge to the same heat
  per length as the library's banded-Cholesky solver.
* Every solve must conserve energy to round-off; the hottest point of the hot
  face must lie above the channel centre; stronger cooling, higher
  conductivity and wider lands must lower the peak.
* For a slender, conductive land the mean flux must approach the fin model
  ``q = (T_aw - T_b + q_rad/h_g) / (1/h_g + t/k + W/(h_c (w + 2 eta h)))``.
"""

from __future__ import annotations

import itertools
import math
import unittest
from unittest.mock import patch

import numpy as np
from scipy.sparse import lil_matrix
from scipy.sparse.linalg import spsolve

from rocket_propulsion.core.errors import DomainError, InputError
from rocket_propulsion.heat_transfer import channel_wall
from rocket_propulsion.heat_transfer.channel_wall import solve_channel_wall_section

BASE = {"pitch_m": 4e-3, "channel_width_m": 1.5e-3, "channel_height_m": 5e-3, "hot_wall_thickness_m": 1e-3,
        "thermal_conductivity_w_m_k": 300.0, "gas_coefficient_w_m2_k": 2.0e4, "recovery_temperature_k": 3500.0,
        "coolant_coefficient_w_m2_k": 6.0e4, "coolant_temperature_k": 200.0}


def independent_section(p: dict, cell: float) -> tuple[float, np.ndarray]:
    """Sparse solve with explicit hot-surface unknowns, not an eliminated Robin source.

    Uniform square cells; geometry must be a multiple of ``cell`` and the
    closeout has the same conductivity. At each hot face, the two equations
    separately impose Fourier conduction and the surface heat balance.
    """
    k, hg, hc = p["thermal_conductivity_w_m_k"], p["gas_coefficient_w_m2_k"], p["coolant_coefficient_w_m2_k"]
    taw, tb, qr = p["recovery_temperature_k"], p["coolant_temperature_k"], p.get("radiative_heat_flux_w_m2", 0.0)
    source = qr + p.get("imposed_heat_flux_w_m2", 0.0)
    ny = round(p["pitch_m"] / 2 / cell)
    nch = round(p["channel_width_m"] / 2 / cell)
    nw = round(p["hot_wall_thickness_m"] / cell)
    nh = round(p["channel_height_m"] / cell)
    nc = round(p.get("closeout_thickness_m", 0.0) / cell)
    solid = np.zeros((nw + nh + nc, ny), dtype=bool)
    solid[:nw, :] = True
    solid[nw:nw + nh, nch:] = True
    solid[nw + nh:, :] = True
    ids = -np.ones(solid.shape, dtype=int)
    ids[solid] = np.arange(solid.sum())
    size = int(solid.sum())
    a = lil_matrix((size + ny, size + ny))
    b = np.zeros(size + ny)
    g_half = 2 * k  # k * face area / half-cell distance
    g_cool = 1.0 / (cell / (2 * k) + 1 / hc) * cell
    for j, i in zip(*np.nonzero(solid), strict=True):
        n = ids[j, i]
        for dj, di in ((0, 1), (1, 0), (0, -1), (-1, 0)):
            jj, ii = j + dj, i + di
            if 0 <= jj < solid.shape[0] and 0 <= ii < ny and solid[jj, ii]:
                a[n, n] += k
                a[n, ids[jj, ii]] -= k
        if j == 0:
            face = size + i
            a[n, n] += g_half
            a[n, face] -= g_half
            a[face, n] -= g_half
            a[face, face] += g_half + hg * cell
            b[face] += (hg * taw + source) * cell
        if j == nw - 1 and i < nch:
            a[n, n] += g_cool
            b[n] += g_cool * tb
        if nw <= j < nw + nh and i == nch:
            a[n, n] += g_cool
            b[n] += g_cool * tb
        if nc and j == nw + nh and i < nch:  # channel ceiling
            a[n, n] += g_cool
            b[n] += g_cool * tb
    temps = spsolve(a.tocsr(), b)
    surface = temps[size:]
    return float(np.sum((hg * (taw - surface) + source) * cell)), surface


def independent_heat(p: dict, cell: float) -> float:
    return independent_section(p, cell)[0]


class ImposedSurfaceFlux(unittest.TestCase):
    def test_local_fourier_and_surface_balances_with_signed_sources(self) -> None:
        # Finite Biot number and one wall cell make a missing half-cell
        # source attenuation visible. A global energy ledger alone cannot.
        base = {**BASE, "channel_width_m": 2e-3, "channel_height_m": 3e-3,
                "thermal_conductivity_w_m_k": 30.0, "recovery_temperature_k": 2500.0,
                "coolant_temperature_k": 300.0, "coolant_coefficient_w_m2_k": 1e4,
                "cells_channel": 2, "cells_land": 1, "cells_wall": 1, "cells_height": 2}
        for radiation, imposed in ((0.0, 0.0), (2e6, 0.0), (0.0, 2e6), (0.0, -2e6), (3e6, -1e6)):
            captured = []
            original = channel_wall._banded_cholesky_solve

            def capture(diagonal, bands, rhs, original=original, captured=captured):
                temperatures = original(diagonal, bands, rhs)
                captured.extend(temperatures)
                return temperatures

            with self.subTest(radiation=radiation, imposed=imposed):
                with patch.object(channel_wall, "_banded_cholesky_solve", side_effect=capture):
                    result = solve_channel_wall_section(**base, radiative_heat_flux_w_m2=radiation,
                                                        imposed_heat_flux_w_m2=imposed)
                face_fluxes = []
                for i, (_, surface) in enumerate(result.hot_surface):
                    film = base["gas_coefficient_w_m2_k"] * (base["recovery_temperature_k"] - surface)
                    conduction = base["thermal_conductivity_w_m_k"] * (surface - captured[i]) / (
                        base["hot_wall_thickness_m"] / 2)
                    self.assertAlmostEqual(film + radiation + imposed, conduction, delta=1e-6)
                    face_fluxes.append(conduction)
                # The public heat result must represent those same physical faces.
                widths = [0.5e-3, 0.5e-3, 1e-3]
                self.assertAlmostEqual(result.heat_per_length_w_m,
                                       math.fsum(q * w for q, w in zip(face_fluxes, widths, strict=True)),
                                       delta=1e-8)

    def test_explicit_surface_unknown_oracle(self) -> None:
        for radiation, imposed, closeout in ((2e6, 0.0, 0.0), (0.0, -2e6, 0.0),
                                              (3e6, -1e6, 1e-3), (0.0, 0.0, 1e-3)):
            p = {**BASE, "channel_width_m": 2e-3, "channel_height_m": 3e-3,
                 "thermal_conductivity_w_m_k": 30.0, "radiative_heat_flux_w_m2": radiation,
                 "imposed_heat_flux_w_m2": imposed, "closeout_thickness_m": closeout}
            with self.subTest(radiation=radiation, imposed=imposed, closeout=closeout):
                heat, surfaces = independent_section(p, 0.5e-3)
                result = solve_channel_wall_section(**p, cells_channel=2, cells_land=2, cells_wall=2,
                                                    cells_height=6, cells_closeout=2)
                self.assertAlmostEqual(result.heat_per_length_w_m, heat, delta=1e-8)
                for (_, actual), expected in zip(result.hot_surface, surfaces, strict=True):
                    self.assertAlmostEqual(actual, expected, delta=1e-8)


class LimitsAndOracles(unittest.TestCase):
    def test_lumped_limit(self) -> None:
        for q_rad in (0.0, 3.0e6):
            p = {**BASE, "thermal_conductivity_w_m_k": 1e10, "radiative_heat_flux_w_m2": q_rad}
            r = solve_channel_wall_section(**p)
            w, half, h = 1.5e-3, 2e-3, 5e-3
            t = ((2e4 * 3500 + q_rad) * half + 6e4 * 200 * (w / 2 + h)) / (2e4 * half + 6e4 * (w / 2 + h))
            with self.subTest(q_rad=q_rad):
                self.assertLess(abs(r.mean_coolant_wall_temperature_k - t), 1e-5 * t)
                self.assertLess(r.peak_hot_wall_temperature_k - r.minimum_hot_wall_temperature_k, 1e-4)

    def test_independent_discretization_converges_to_the_same_heat(self) -> None:
        p = {**BASE, "pitch_m": 4.0e-3, "channel_width_m": 1.6e-3, "channel_height_m": 4.0e-3,
             "hot_wall_thickness_m": 0.8e-3, "radiative_heat_flux_w_m2": 1.0e6}
        reference = independent_heat(p, 0.05e-3)
        # identical 0.05 mm cells: two independent assemblies and solvers agree to round-off
        same = solve_channel_wall_section(**p, cells_channel=16, cells_land=24, cells_wall=16, cells_height=80)
        self.assertLess(abs(same.heat_per_length_w_m - reference), 1e-9 * reference)
        # and the library's finer grid stays within the discretization error
        fine = solve_channel_wall_section(**p, cells_channel=32, cells_land=48, cells_wall=32, cells_height=160)
        self.assertLess(abs(fine.heat_per_length_w_m - reference) / reference, 2e-3)

    def test_refinement_converges(self) -> None:
        heats = [solve_channel_wall_section(**BASE, cells_channel=6 * n, cells_land=4 * n, cells_wall=4 * n,
                                            cells_height=8 * n).heat_per_length_w_m for n in (1, 2, 4, 8)]
        steps = [abs(b - a) for a, b in itertools.pairwise(heats)]
        self.assertTrue(steps[0] > steps[1] > steps[2])
        self.assertLess(steps[-1] / heats[-1], 1e-3)


class Closeout(unittest.TestCase):
    def test_independent_assembly_with_a_closeout(self) -> None:
        p = {**BASE, "pitch_m": 4.0e-3, "channel_width_m": 1.6e-3, "channel_height_m": 4.0e-3,
             "hot_wall_thickness_m": 0.8e-3, "closeout_thickness_m": 1.0e-3}
        reference = independent_heat(p, 0.1e-3)
        same = solve_channel_wall_section(**p, cells_channel=8, cells_land=12, cells_wall=8, cells_height=40,
                                          cells_closeout=10)
        self.assertLess(abs(same.heat_per_length_w_m - reference), 1e-9 * reference)
        self.assertLess(abs(same.energy_residual_w_m), 1e-11 * reference)

    def test_thin_closeout_limit_and_trends(self) -> None:
        bare = solve_channel_wall_section(**BASE)
        thin = solve_channel_wall_section(**BASE, closeout_thickness_m=1e-8)
        self.assertLess(abs(thin.heat_per_length_w_m / bare.heat_per_length_w_m - 1), 1e-4)
        self.assertLess(abs(thin.peak_hot_wall_temperature_k - bare.peak_hot_wall_temperature_k), 0.05)
        peaks = [solve_channel_wall_section(**BASE, closeout_thickness_m=1e-3, closeout_conductivity_w_m_k=kc)
                 .peak_hot_wall_temperature_k for kc in (10.0, 100.0, 1000.0)]
        self.assertEqual(peaks, sorted(peaks, reverse=True))  # a better closeout path cools more
        self.assertLess(peaks[0], bare.peak_hot_wall_temperature_k + 1e-9)

    def test_outer_film_closes_the_energy_balance(self) -> None:
        r = solve_channel_wall_section(**BASE, closeout_thickness_m=1e-3, outer_coefficient_w_m2_k=50.0,
                                       outer_temperature_k=300.0)
        self.assertLess(abs(r.energy_residual_w_m), 1e-11 * r.heat_per_length_w_m)
        for change in ({"closeout_thickness_m": -1e-3}, {"closeout_thickness_m": math.nan},
                       {"closeout_thickness_m": 1e-3, "closeout_conductivity_w_m_k": 0.0},
                       {"outer_coefficient_w_m2_k": 10.0},
                       {"closeout_thickness_m": 1e-3, "outer_coefficient_w_m2_k": 10.0, "outer_temperature_k": -5.0},
                       {"imposed_heat_flux_w_m2": math.inf}):
            with self.subTest(change=change), self.assertRaises(DomainError):
                solve_channel_wall_section(**{**BASE, **change})
        # an imposed flux of either sign enters the balance
        for imposed in (-2.0e6, 3.0e6):
            q = solve_channel_wall_section(**BASE, imposed_heat_flux_w_m2=imposed)
            with self.subTest(imposed=imposed):
                self.assertLess(abs(q.energy_residual_w_m), 1e-11 * abs(q.heat_per_length_w_m))


class ConservationAndTrends(unittest.TestCase):
    def test_energy_is_conserved_over_random_designs(self) -> None:
        rng = np.random.default_rng(6)
        for _ in range(40):
            pitch = rng.uniform(2e-3, 8e-3)
            p = {**BASE, "pitch_m": pitch, "channel_width_m": rng.uniform(0.2, 0.8) * pitch,
                 "channel_height_m": rng.uniform(1e-3, 8e-3), "hot_wall_thickness_m": rng.uniform(0.3e-3, 2e-3),
                 "thermal_conductivity_w_m_k": rng.uniform(15.0, 400.0),
                 "gas_coefficient_w_m2_k": rng.uniform(1e3, 5e4), "coolant_coefficient_w_m2_k": rng.uniform(1e3, 2e5),
                 "radiative_heat_flux_w_m2": rng.uniform(0.0, 2e6)}
            r = solve_channel_wall_section(**p)
            with self.subTest(p=p):
                self.assertLess(abs(r.energy_residual_w_m), 1e-11 * r.heat_per_length_w_m)
                self.assertGreater(r.peak_hot_wall_temperature_k, r.mean_hot_wall_temperature_k - 1e-9)
                self.assertGreater(r.mean_hot_wall_temperature_k, r.mean_coolant_wall_temperature_k)
                self.assertGreater(r.mean_coolant_wall_temperature_k, p["coolant_temperature_k"])
                self.assertLess(r.peak_hot_wall_temperature_k, p["recovery_temperature_k"]
                                + p["radiative_heat_flux_w_m2"] / p["gas_coefficient_w_m2_k"])

    def test_peak_sits_over_the_channel_and_cooling_lowers_it(self) -> None:
        base = solve_channel_wall_section(**BASE)
        self.assertLess(base.peak_location_y_m, 0.5 * BASE["channel_width_m"])
        for key, factor in (("coolant_coefficient_w_m2_k", 2.0), ("thermal_conductivity_w_m_k", 2.0)):
            with self.subTest(change=key):
                better = solve_channel_wall_section(**{**BASE, key: BASE[key] * factor})
                self.assertLess(better.peak_hot_wall_temperature_k, base.peak_hot_wall_temperature_k)

    def test_slender_conductive_land_approaches_the_fin_model(self) -> None:
        p = {**BASE, "thermal_conductivity_w_m_k": 5000.0, "hot_wall_thickness_m": 0.3e-3, "channel_height_m": 3e-3,
             "pitch_m": 3.0e-3, "channel_width_m": 2.0e-3}
        r = solve_channel_wall_section(**p, cells_channel=24, cells_land=12, cells_wall=8, cells_height=48)
        k, hc, land, h = 5000.0, 6.0e4, 1.0e-3, 3e-3
        m = math.sqrt(hc * 2.0 / (k * land))
        eta = math.tanh(m * h) / (m * h)
        fin_q = (3500.0 - 200.0) / (1 / 2e4 + 0.3e-3 / k + 3.0e-3 / (hc * (2.0e-3 + 2 * eta * h)))
        self.assertLess(abs(r.mean_heat_flux_w_m2 - fin_q) / fin_q, 0.02)


class HostileInputs(unittest.TestCase):
    def test_wide_grid_attack_is_refused_before_matrix_solution(self) -> None:
        # Only 20,000 rectangular-grid cells, but ~100 million entries per
        # band matrix and ~10^12 factorization work units without a guard.
        with patch.object(channel_wall, "_banded_cholesky_solve") as solve:
            with self.assertRaises(InputError) as raised:
                solve_channel_wall_section(**BASE, cells_channel=9999, cells_land=1,
                                           cells_wall=1, cells_height=1)
            self.assertEqual(raised.exception.code, "wall_grid_budget")
            solve.assert_not_called()

    def test_malformed_scalars_are_not_coerced(self) -> None:
        for field in BASE:
            for value in (True, "1", None, [1.0]):
                with self.subTest(field=field, value=value), self.assertRaises(DomainError):
                    solve_channel_wall_section(**{**BASE, field: value})

    def test_unrepresentable_or_nonphysical_sourced_wall_is_refused(self) -> None:
        for change in ({"gas_coefficient_w_m2_k": 1e-308, "radiative_heat_flux_w_m2": 1e308},
                       {"thermal_conductivity_w_m_k": 1e-320},
                       {"hot_wall_thickness_m": 5e-324},
                       {"imposed_heat_flux_w_m2": -1e12}):
            with self.subTest(change=change), self.assertRaises(DomainError):
                solve_channel_wall_section(**{**BASE, **change})

    def test_refused(self) -> None:
        for change in ({"channel_width_m": 4e-3}, {"channel_width_m": 5e-3}, {"pitch_m": 0.0},
                       {"thermal_conductivity_w_m_k": -1.0}, {"gas_coefficient_w_m2_k": math.nan},
                       {"coolant_temperature_k": 4000.0}, {"radiative_heat_flux_w_m2": -1.0},
                       {"radiative_heat_flux_w_m2": math.inf}, {"hot_wall_thickness_m": math.inf}):
            with self.subTest(change=change), self.assertRaises(DomainError):
                solve_channel_wall_section(**{**BASE, **change})
        for change in ({"cells_channel": 0}, {"cells_land": True}, {"cells_wall": 2.0},
                       {"cells_height": 10_000}):
            with self.subTest(change=change), self.assertRaises(InputError):
                solve_channel_wall_section(**{**BASE, **change})


if __name__ == "__main__":
    unittest.main()
