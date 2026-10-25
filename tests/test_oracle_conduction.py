"""Conduction mathematics against independent implementations (scipy, mpmath).

The package implements its own Bessel functions, integrated error function,
transient eigenvalues and series solutions in pure Python. Here each is
compared with an independent library or with a different solution method:

* ``J_n``, ``I_0``, ``I_1``, ``K_0``, ``K_1`` and their scaled forms against
  ``scipy.special`` over ten decades of argument;
* the integrated complementary error function against mpmath at 50 digits;
* transient eigenvalues against brentq roots of the characteristic
  equations (slab ``l tan l = Bi``, cylinder ``l J1/J0 = Bi``, sphere
  ``1 - l cot l = Bi``);
* the series transient solution against a method-of-lines discretisation of
  the heat equation solved by scipy's stiff BDF integrator, for all three
  shapes;
* straight and annular fins against scipy's boundary-value solver on the fin
  differential equation;
* the semi-infinite convective-surface solution against mpmath's erfc with
  a cancellation-safe form, deep into the large-argument tail.
"""

from __future__ import annotations

import itertools
import math
import unittest

import mpmath
import numpy as np
from scipy import integrate, optimize, sparse, special

from rocket_propulsion.heat_transfer import conduction as c

mpmath.mp.dps = 50
ARGS = np.concatenate([np.logspace(-8, 2.6, 160), [0.0]])


def rel(a: float, b: float) -> float:
    return abs(a - b) / max(abs(b), 1e-300)


class BesselAgainstScipy(unittest.TestCase):
    def test_modified_bessel_functions(self) -> None:
        worst: dict[str, float] = {}
        for x in ARGS:
            x = float(x)
            pairs = [("i0e", c.bessel_i0e(x), special.i0e(x)), ("i1e", c.bessel_i1e(x), special.i1e(x))]
            if x > 0:
                pairs += [("k0e", c.bessel_k0e(x), special.k0e(x)), ("k1e", c.bessel_k1e(x), special.k1e(x))]
            if x < 700:
                pairs += [("i0", c.bessel_i0(x), special.i0(x)), ("i1", c.bessel_i1(x), special.i1(x))]
                if x > 0:
                    pairs += [("k0", c.bessel_k0(x), special.k0(x)), ("k1", c.bessel_k1(x), special.k1(x))]
            for name, mine, ref in pairs:
                if ref == 0.0:
                    self.assertEqual(mine, 0.0)
                    continue
                worst[name] = max(worst.get(name, 0.0), rel(mine, float(ref)))
        for name, err in worst.items():
            with self.subTest(function=name):
                self.assertLess(err, 5e-13)

    def test_bessel_j_of_integer_order_absolute(self) -> None:
        # J_n oscillates through zero, so the meaningful bound is absolute.
        for order in range(9):
            for x in np.linspace(-60.0, 60.0, 481):
                mine = c.bessel_jn(order, float(x))
                with self.subTest(order=order, x=float(x)):
                    self.assertLess(abs(mine - special.jv(order, x)), 2e-14)
        self.assertLess(abs(c.bessel_j0(7.3) - special.j0(7.3)), 1e-15)
        self.assertLess(abs(c.bessel_j1(7.3) - special.j1(7.3)), 1e-15)

    def test_recurrence_and_wronskian_identities(self) -> None:
        # J_{n-1} + J_{n+1} = (2n/x) J_n ; I_0 K_1 + I_1 K_0 = 1/x (DLMF 10.5, 10.28.2)
        for x in (0.3, 2.7, 9.1, 31.4):
            for n in range(1, 6):
                lhs = c.bessel_jn(n - 1, x) + c.bessel_jn(n + 1, x)
                self.assertLess(abs(lhs - 2 * n / x * c.bessel_jn(n, x)), 3e-14)
            wronskian = c.bessel_i0e(x) * c.bessel_k1e(x) + c.bessel_i1e(x) * c.bessel_k0e(x)
            self.assertLess(rel(wronskian, 1 / x), 1e-13)


class IntegratedErfcAgainstMpmath(unittest.TestCase):
    def test_against_fifty_digit_reference(self) -> None:
        for z in [0.0, 1e-6, 0.1, 0.5, 1.0, 2.0, 4.0, 6.0, 8.0, 10.0, 15.0, 25.0]:
            exact = mpmath.exp(-mpmath.mpf(z) ** 2) / mpmath.sqrt(mpmath.pi) - z * mpmath.erfc(z)
            with self.subTest(z=z):
                self.assertLess(rel(c.integrated_erfc(z), float(exact)), 1e-12)

    def test_semi_infinite_convective_surface_against_mpmath(self) -> None:
        alpha, k, h, ti, tf = 1.1e-5, 15.0, 800.0, 300.0, 1500.0
        for depth, time in ((0.0, 1.0), (0.002, 5.0), (0.01, 30.0), (0.05, 2.0), (0.0, 1e4)):
            r = c.semi_infinite_convective_surface(
                depth_m=depth, time_s=time, thermal_diffusivity_m2_s=alpha, thermal_conductivity_w_m_k=k,
                heat_transfer_coefficient_w_m2_k=h, initial_temperature_k=ti, fluid_temperature_k=tf)
            x, t = mpmath.mpf(depth), mpmath.mpf(time)
            eta = x / (2 * mpmath.sqrt(alpha * t))
            b = h * mpmath.sqrt(alpha * t) / k
            frac = mpmath.erfc(eta) - mpmath.exp(h * x / k + b * b) * mpmath.erfc(eta + b)
            exact = ti + (tf - ti) * frac
            with self.subTest(depth=depth, time=time):
                self.assertLess(abs(r.temperature_k - float(exact)), 1e-9 * (tf - ti))


def _characteristic(shape: str, bi: float):
    if shape == "slab":
        return lambda lam: lam * math.sin(lam) - bi * math.cos(lam)
    if shape == "cylinder":
        return lambda lam: lam * special.j1(lam) - bi * special.j0(lam)
    return lambda lam: 1 - lam / math.tan(lam) - bi if lam % math.pi > 1e-9 else float("nan")


class EigenvaluesAgainstBrentq(unittest.TestCase):
    def test_first_twenty_roots_for_every_shape(self) -> None:
        for shape in ("slab", "cylinder", "sphere"):
            for bi in (1e-4, 0.01, 0.1, 1.0, 5.0, 50.0, 1e3):
                roots = c.transient_eigenvalues(shape, bi, 20)
                f = _characteristic(shape, bi)
                for k, root in enumerate(roots):
                    lo, hi = root - 1e-6 * max(root, 1.0), root + 1e-6 * max(root, 1.0)
                    ref = optimize.brentq(f, lo, hi, xtol=1e-15, rtol=1e-15)
                    with self.subTest(shape=shape, bi=bi, k=k):
                        self.assertLess(abs(root - ref), 1e-11 * max(1.0, ref))
                # strictly increasing, one root per branch
                self.assertTrue(all(b > a for a, b in itertools.pairwise(roots)))


def _method_of_lines(shape: str, bi: float, fo: float, position: float, n: int = 400) -> float:
    """Dimensionless theta from the heat equation by finite volumes + BDF."""
    m = {"slab": 0, "cylinder": 1, "sphere": 2}[shape]
    edges = np.linspace(0.0, 1.0, n + 1)
    centres = 0.5 * (edges[1:] + edges[:-1])
    area = edges ** m
    volume = (edges[1:] ** (m + 1) - edges[:-1] ** (m + 1)) / (m + 1)

    def rhs(_t: float, th: np.ndarray) -> np.ndarray:
        flux = np.zeros(n + 1)
        flux[1:-1] = area[1:-1] * (th[1:] - th[:-1]) / (centres[1:] - centres[:-1])
        # Robin boundary through a half cell: -dtheta/dx = Bi theta at x = 1
        half = 1.0 - centres[-1]
        surface = th[-1] / (1.0 + bi * half)
        flux[-1] = -area[-1] * bi * surface
        return (flux[1:] - flux[:-1]) / volume

    pattern = sparse.diags([np.ones(n - 1), np.ones(n), np.ones(n - 1)], [-1, 0, 1])
    sol = integrate.solve_ivp(rhs, (0.0, fo), np.ones(n), method="BDF", rtol=1e-10, atol=1e-12,
                              jac_sparsity=pattern)
    return float(np.interp(position, centres, sol.y[:, -1]))


class SeriesAgainstMethodOfLines(unittest.TestCase):
    def test_three_shapes(self) -> None:
        for shape in ("slab", "cylinder", "sphere"):
            for bi, fo, pos in ((0.5, 0.1, 0.0), (5.0, 0.3, 0.6), (50.0, 0.05, 0.9)):
                series = c.series_transient_solution(shape, bi, fo, pos).theta
                numeric = _method_of_lines(shape, bi, fo, pos)
                with self.subTest(shape=shape, bi=bi, fo=fo, pos=pos):
                    # second-order finite volumes on 400 cells: agreement to ~1e-4
                    self.assertLess(abs(series - numeric), 3e-4)


class FinsAgainstBoundaryValueSolver(unittest.TestCase):
    def test_straight_fin_with_convective_tip(self) -> None:
        k, h, w, t, length = 200.0, 60.0, 0.05, 0.002, 0.04
        p, a = 2 * (w + t), w * t
        r = c.straight_fin(thermal_conductivity_w_m_k=k, heat_transfer_coefficient_w_m2_k=h, perimeter_m=p,
                           cross_section_area_m2=a, length_m=length, base_temperature_k=400.0,
                           fluid_temperature_k=300.0, tip_heat_transfer_coefficient_w_m2_k=h)
        m2 = h * p / (k * a)

        def ode(x: np.ndarray, y: np.ndarray) -> np.ndarray:
            return np.vstack([y[1], m2 * y[0]])

        def bc(ya: np.ndarray, yb: np.ndarray) -> np.ndarray:
            return np.array([ya[0] - 100.0, k * yb[1] + h * yb[0]])
        x = np.linspace(0, length, 200)
        sol = integrate.solve_bvp(ode, bc, x, np.vstack([np.full_like(x, 100.0), np.zeros_like(x)]), tol=1e-10)
        q_base = -k * a * sol.sol(0.0)[1]
        self.assertTrue(sol.success)
        self.assertLess(rel(r.heat_rate_w, q_base), 1e-6)

    def test_annular_fin_against_bvp(self) -> None:
        k, h, t, r1, r2 = 180.0, 120.0, 0.0015, 0.01, 0.035
        r = c.annular_fin(thermal_conductivity_w_m_k=k, heat_transfer_coefficient_w_m2_k=h, thickness_m=t,
                          inner_radius_m=r1, outer_radius_m=r2, base_temperature_k=450.0,
                          fluid_temperature_k=300.0)
        m2 = 2 * h / (k * t)

        def ode(rr: np.ndarray, y: np.ndarray) -> np.ndarray:
            return np.vstack([y[1], m2 * y[0] - y[1] / rr])

        def bc(ya: np.ndarray, yb: np.ndarray) -> np.ndarray:
            return np.array([ya[0] - 150.0, yb[1]])
        rr = np.linspace(r1, r2, 300)
        sol = integrate.solve_bvp(ode, bc, rr, np.vstack([np.full_like(rr, 150.0), np.zeros_like(rr)]), tol=1e-8,
                                  max_nodes=100000)
        q_base = -k * (2 * math.pi * r1 * t) * sol.sol(r1)[1]
        self.assertTrue(sol.success)
        self.assertLess(rel(r.heat_rate_w, q_base), 1e-6)
