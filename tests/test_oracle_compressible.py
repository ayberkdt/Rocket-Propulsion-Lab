"""Compressible-flow relations against independent physics, not against their own formulas.

Every check here rebuilds a quantity from a different principle than the code
uses, so a shared algebra error cannot pass both sides:

* normal and oblique shocks conserve mass, momentum and energy across the
  wave (Rankine-Hugoniot), computed from the returned ratios alone, and never
  decrease entropy;
* the oblique-shock wave angle satisfies the theta-beta-M relation and the
  normal Mach components reproduce the normal-shock solution;
* the Prandtl-Meyer function and the Fanno friction parameter are recomputed
  by adaptive quadrature (scipy) of their defining differential relations;
* the isentropic area ratio equals the mass-flux ratio built from the
  returned density and temperature ratios;
* every ``mach_from_isentropic`` input kind inverts ``calculate_isentropic``
  over random Mach numbers and ratios of specific heats;
* the Taylor-Maccoll cone solution is reproduced by integrating the ODE
  from the returned shock to the cone surface with scipy's RK45.

Random cases use fixed seeds so failures are reproducible.
"""

from __future__ import annotations

import math
import unittest

import numpy as np
from scipy import integrate

from rocket_propulsion.compressible.cone_flow import solve_cone_flow
from rocket_propulsion.compressible.duct_flow import calculate_fanno, calculate_rayleigh
from rocket_propulsion.compressible.isentropic import calculate_isentropic, mach_from_isentropic
from rocket_propulsion.compressible.normal_shock import calculate_normal_shock
from rocket_propulsion.compressible.oblique_shock import calculate_oblique_shock
from rocket_propulsion.core.errors import RocketPropulsionError

RNG = np.random.default_rng(20261025)
GAMMAS = (1.1, 1.2, 1.3, 1.4, 5.0 / 3.0)


def rel(a: float, b: float) -> float:
    return abs(a - b) / max(abs(b), 1e-300)


class NormalShockConservation(unittest.TestCase):
    def test_rankine_hugoniot_from_returned_ratios(self) -> None:
        for gamma in GAMMAS:
            for m1 in (1.0001, 1.05, 1.5, 2.0, 3.7, 8.0, 25.0):
                r = calculate_normal_shock(m1, gamma)
                m2 = r.downstream_mach
                # Mass: rho1 u1 = rho2 u2 with u = M a, a ~ sqrt(T)
                mass = r.density_ratio * (m2 / m1) * math.sqrt(r.temperature_ratio)
                # Momentum: p + rho u^2 = p (1 + gamma M^2)
                momentum = r.pressure_ratio * (1 + gamma * m2 * m2) / (1 + gamma * m1 * m1)
                # Energy: h0 constant, T0 = T (1 + (g-1)/2 M^2)
                energy = r.temperature_ratio * (1 + 0.5 * (gamma - 1) * m2 * m2) / (
                    1 + 0.5 * (gamma - 1) * m1 * m1)
                # Ideal gas: p = rho R T
                state = r.pressure_ratio / (r.density_ratio * r.temperature_ratio)
                with self.subTest(gamma=gamma, m1=m1):
                    self.assertLess(abs(mass - 1), 1e-12)
                    self.assertLess(abs(momentum - 1), 1e-12)
                    self.assertLess(abs(energy - 1), 1e-12)
                    self.assertLess(abs(state - 1), 1e-12)
                    self.assertLess(m2, 1.0)

    def test_entropy_never_decreases_and_matches_total_pressure_loss(self) -> None:
        for gamma in GAMMAS:
            for m1 in np.linspace(1.0001, 12.0, 45):
                r = calculate_normal_shock(float(m1), gamma)
                ds_over_r = (math.log(r.temperature_ratio) * gamma / (gamma - 1)
                             - math.log(r.pressure_ratio))
                with self.subTest(gamma=gamma, m1=float(m1)):
                    self.assertGreaterEqual(ds_over_r, -1e-13)
                    # s2 - s1 = -R ln(p02/p01)
                    self.assertLess(abs(ds_over_r + math.log(r.stagnation_pressure_ratio)), 1e-11)

    def test_subsonic_upstream_is_refused(self) -> None:
        for m1 in (0.0, 0.5, 0.999999):
            with self.assertRaises(RocketPropulsionError):
                calculate_normal_shock(m1, 1.4)


class ObliqueShockGeometry(unittest.TestCase):
    def test_theta_beta_mach_identity_and_normal_components(self) -> None:
        for gamma in (1.2, 1.4):
            for m1 in (1.6, 2.0, 3.0, 6.0):
                limit = calculate_oblique_shock(upstream_mach=m1, turning_angle_deg=0.5,
                                                gamma=gamma).maximum_turning_angle_deg
                for frac in (0.05, 0.3, 0.6, 0.95):
                    theta = frac * limit
                    for branch in ("weak", "strong"):
                        r = calculate_oblique_shock(upstream_mach=m1, turning_angle_deg=theta,
                                                    gamma=gamma, branch=branch)
                        beta = math.radians(r.wave_angle_deg)
                        rhs = 2 / math.tan(beta) * (m1 ** 2 * math.sin(beta) ** 2 - 1) / (
                            m1 ** 2 * (gamma + math.cos(2 * beta)) + 2)
                        normal = calculate_normal_shock(m1 * math.sin(beta), gamma)
                        with self.subTest(gamma=gamma, m1=m1, theta=theta, branch=branch):
                            self.assertLess(abs(math.tan(math.radians(theta)) - rhs), 1e-10)
                            self.assertLess(rel(r.upstream_normal_mach, m1 * math.sin(beta)), 1e-12)
                            self.assertLess(rel(r.pressure_ratio, normal.pressure_ratio), 1e-11)
                            self.assertLess(rel(r.stagnation_pressure_ratio,
                                                normal.stagnation_pressure_ratio), 1e-10)
                            self.assertLess(rel(r.downstream_mach * math.sin(beta - math.radians(theta)),
                                                normal.downstream_mach), 1e-9)

    def test_weak_branch_is_weaker_than_strong(self) -> None:
        for m1 in (1.8, 2.5, 4.0):
            theta = 0.5 * calculate_oblique_shock(upstream_mach=m1, turning_angle_deg=1.0).maximum_turning_angle_deg
            weak = calculate_oblique_shock(upstream_mach=m1, turning_angle_deg=theta, branch="weak")
            strong = calculate_oblique_shock(upstream_mach=m1, turning_angle_deg=theta, branch="strong")
            self.assertLess(weak.wave_angle_deg, strong.wave_angle_deg)
            self.assertLess(weak.pressure_ratio, strong.pressure_ratio)
            self.assertGreater(weak.stagnation_pressure_ratio, strong.stagnation_pressure_ratio)
            self.assertLess(strong.downstream_mach, 1.0)

    def test_detached_turning_is_refused(self) -> None:
        limit = calculate_oblique_shock(upstream_mach=2.0, turning_angle_deg=1.0).maximum_turning_angle_deg
        with self.assertRaises(RocketPropulsionError):
            calculate_oblique_shock(upstream_mach=2.0, turning_angle_deg=limit * 1.001)


class IsentropicByQuadrature(unittest.TestCase):
    def test_prandtl_meyer_equals_integral_of_its_ode(self) -> None:
        for gamma in GAMMAS:
            for m in (1.1, 1.5, 2.0, 3.0, 5.0, 10.0):
                def integrand(x: float, g: float = gamma) -> float:
                    return math.sqrt(x * x - 1) / (1 + 0.5 * (g - 1) * x * x) / x
                nu, _ = integrate.quad(integrand, 1.0, m, epsabs=1e-13, epsrel=1e-13, limit=200)
                with self.subTest(gamma=gamma, m=m):
                    self.assertLess(abs(math.degrees(nu) - calculate_isentropic(m, gamma).prandtl_meyer_deg), 1e-9)

    def test_area_ratio_is_the_mass_flux_ratio(self) -> None:
        for gamma in GAMMAS:
            star = calculate_isentropic(1.0, gamma)
            for m in (0.05, 0.3, 0.8, 1.0, 1.7, 4.0, 15.0):
                r = calculate_isentropic(m, gamma)
                flux = r.density_total_ratio * m * math.sqrt(r.temperature_total_ratio)
                flux_star = star.density_total_ratio * math.sqrt(star.temperature_total_ratio)
                with self.subTest(gamma=gamma, m=m):
                    self.assertLess(rel(r.area_critical_ratio, flux_star / flux), 1e-12)
                    self.assertLess(rel(r.pressure_total_ratio,
                                        r.density_total_ratio * r.temperature_total_ratio), 1e-13)

    def test_every_inverse_kind_round_trips(self) -> None:
        kinds_sub = ("pressure_total_ratio", "density_total_ratio", "temperature_total_ratio",
                     "pressure_critical_ratio", "density_critical_ratio", "temperature_critical_ratio",
                     "area_ratio_subsonic")
        kinds_sup = ("mach_angle_deg", "prandtl_meyer_deg", "area_ratio_supersonic") + kinds_sub[:-1]
        failures = []
        for _ in range(250):
            gamma = float(RNG.uniform(1.05, 1.67))
            sub = float(RNG.uniform(0.02, 0.98))
            sup = float(10 ** RNG.uniform(math.log10(1.02), math.log10(30.0)))
            for m, kinds in ((sub, kinds_sub), (sup, kinds_sup)):
                ref = calculate_isentropic(m, gamma)
                for kind in kinds:
                    attr = {"area_ratio_subsonic": "area_critical_ratio",
                            "area_ratio_supersonic": "area_critical_ratio"}.get(kind, kind)
                    value = getattr(ref, attr)
                    back = mach_from_isentropic(kind, value, gamma)
                    if rel(back, m) > 2e-12:
                        failures.append((kind, gamma, m, back))
        self.assertEqual(failures, [])


class DuctFlowByQuadrature(unittest.TestCase):
    def test_fanno_friction_parameter_is_the_integral_of_the_fanno_ode(self) -> None:
        for gamma in (1.2, 1.4):
            for m in (0.2, 0.5, 0.9, 1.3, 2.0, 3.5):
                def integrand(x: float, g: float = gamma) -> float:
                    # d(4fL/D) = (1 - M^2) / (g M^4 (1 + (g-1)/2 M^2)) d(M^2)
                    return (1 - x) / (g * x * x * (1 + 0.5 * (g - 1) * x)) * 1.0
                value, _ = integrate.quad(integrand, m * m, 1.0, epsabs=1e-13, epsrel=1e-12, limit=200)
                with self.subTest(gamma=gamma, m=m):
                    self.assertLess(abs(abs(value) - calculate_fanno(m, gamma).friction_parameter_to_sonic), 1e-9)

    def test_rayleigh_entropy_matches_temperature_and_pressure(self) -> None:
        for gamma in (1.2, 1.4):
            for m in (0.1, 0.5, 0.95, 1.5, 3.0):
                r = calculate_rayleigh(m, gamma)
                ds = gamma / (gamma - 1) * math.log(r.temperature_critical_ratio) - math.log(r.pressure_critical_ratio)
                with self.subTest(gamma=gamma, m=m):
                    self.assertLess(abs(abs(ds) - abs(r.entropy_difference_over_r)), 1e-11)
                    # Rayleigh line: p + rho u^2 constant
                    self.assertLess(rel(r.pressure_critical_ratio * (1 + gamma * m * m), 1 + gamma), 1e-12)

    def test_fanno_energy_is_constant(self) -> None:
        for m in (0.3, 0.8, 1.4, 2.6):
            r = calculate_fanno(m, 1.4)
            # T0 constant: T/T* (1 + 0.2 M^2) = 1.2
            self.assertLess(rel(r.temperature_critical_ratio * (1 + 0.2 * m * m), 1.2), 1e-13)
            # continuity: rho u = const
            self.assertLess(rel(r.density_critical_ratio * r.velocity_critical_ratio, 1.0), 1e-13)


class TaylorMaccollByIntegration(unittest.TestCase):
    def test_cone_surface_from_scipy_integration(self) -> None:
        gamma = 1.4
        for mach, cone in ((2.0, 10.0), (3.0, 20.0), (5.0, 15.0), (1.8, 5.0)):
            r = solve_cone_flow(mach=mach, cone_half_angle_deg=cone, gamma=gamma)
            beta = math.radians(r.shock_angle_deg)
            delta = math.radians(r.post_shock_deflection_deg)
            m2 = r.post_shock_mach
            # Nondimensional velocity V' = V / V_max
            v2 = (2 / ((gamma - 1) * m2 * m2) + 1) ** -0.5
            vr, vt = v2 * math.cos(beta - delta), -v2 * math.sin(beta - delta)

            def tm(theta: float, y: np.ndarray, g: float = gamma) -> list[float]:
                a, b = y
                c = 0.5 * (g - 1) * (1 - a * a - b * b)
                d2 = (b * b * a - c * (2 * a + b / math.tan(theta))) / (c - b * b)
                return [b, d2]

            def surface(theta: float, y: np.ndarray) -> float:
                return y[1]
            surface.terminal = True  # type: ignore[attr-defined]  # scipy event protocol attribute
            sol = integrate.solve_ivp(tm, (beta, 1e-4), [vr, vt], events=surface, rtol=1e-11, atol=1e-13)
            theta_c = math.degrees(sol.t_events[0][0])
            v_c = sol.y_events[0][0][0]
            m_c = math.sqrt(2 / ((gamma - 1) * (1 / (v_c * v_c) - 1)))
            with self.subTest(mach=mach, cone=cone):
                self.assertLess(abs(theta_c - cone), 1e-6)
                self.assertLess(rel(m_c, r.surface_mach), 1e-6)
