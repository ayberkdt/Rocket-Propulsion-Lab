"""Thermodynamic consistency of the reference equations of state.

A Helmholtz equation of state determines every property from one function,
so its outputs must satisfy exact thermodynamic identities. These tests check
the identities on random states of every reference fluid, with central
finite differences of the package's own ``state_rho_t``/``state_tp`` calls
as the independent route to each derivative:

* ``(dh/dT)_p = c_p``, ``(ds/dT)_p = c_p/T`` and ``(du/dT)_rho = c_v``;
* ``(dp/drho)_T`` and ``(dp/dT)_rho`` equal the analytic derivatives
  reported on the state;
* ``c_p - c_v = T (dp/dT)_rho^2 / (rho^2 (dp/drho)_T)`` and
  ``w^2 = (c_p/c_v)(dp/drho)_T``;
* Gibbs energy ``g = h - T s``; Joule-Thomson coefficient
  ``mu_JT = (T (dv/dT)_p - v)/c_p``;
* on the saturation line: equal Gibbs energies, Clausius-Clapeyron
  ``dp_sat/dT = h_fg/(T (v_v - v_l))`` from neighbouring saturation solves;
* ``state_tp`` inverts ``state_rho_t``.
"""

from __future__ import annotations

import unittest

import numpy as np

from rocket_propulsion.core.errors import RocketPropulsionError
from rocket_propulsion.thermodynamics.real_fluids import (
    REFERENCE_FLUID_KEYS,
    reference_fluid,
    reference_fluid_model,
)

RNG = np.random.default_rng(7123)


def rel(a: float, b: float) -> float:
    return abs(a - b) / max(abs(b), 1e-300)


def random_states(fluid: str, count: int) -> list[tuple[float, float]]:
    """Single-phase (T, p) pairs inside the model range, away from the critical point."""
    model = reference_fluid_model(fluid)
    states = []
    tc, pc = model.t_crit_k, model.p_crit_pa
    while len(states) < count:
        t = float(RNG.uniform(max(model.t_min_k * 1.05, 0.55 * tc), min(model.t_max_k, 6.0 * tc)))
        p = float(10 ** RNG.uniform(np.log10(2e4), np.log10(min(model.p_max_pa, 30e6))))
        if abs(t / tc - 1) < 0.04 and abs(p / pc - 1) < 0.2:
            continue
        states.append((t, p))
    return states


class FirstDerivatives(unittest.TestCase):
    def test_heat_capacities_and_pressure_derivatives(self) -> None:
        failures = []
        for fluid in REFERENCE_FLUID_KEYS:
            f = reference_fluid(fluid)
            for t, p in random_states(fluid, 12):
                try:
                    s0 = f.state_tp(t, p)
                except RocketPropulsionError:
                    continue
                rho = s0.density_kg_m3

                def richardson(g, step: float) -> float:
                    """Fourth-order central difference by Richardson extrapolation."""
                    d1 = (g(step) - g(-step)) / (2 * step)
                    d2 = (g(step / 2) - g(-step / 2)) / step
                    return (4 * d2 - d1) / 3

                dt = 2e-4 * t
                phases = {f.state_tp(t + sign * dt, p).phase for sign in (-1, 1)}
                if phases != {s0.phase}:
                    continue
                cp_fd = richardson(lambda e, f=f, t=t, p=p: f.state_tp(t + e, p).enthalpy_j_kg, dt)
                ds_fd = richardson(lambda e, f=f, t=t, p=p: f.state_tp(t + e, p).entropy_j_kg_k, dt)
                cv_fd = richardson(lambda e, f=f, t=t, rho=rho: f.state_rho_t(rho, t + e).internal_energy_j_kg, dt)
                dpdt_fd = richardson(lambda e, f=f, t=t, rho=rho: f.state_rho_t(rho, t + e).pressure_pa, dt)
                dpdrho_fd = richardson(lambda e, f=f, t=t, rho=rho: f.state_rho_t(rho + e, t).pressure_pa, 2e-5 * rho)
                checks = {
                    "cp": rel(cp_fd, s0.cp_j_kg_k),
                    "ds/dT": rel(ds_fd, s0.cp_j_kg_k / t),
                    "cv": rel(cv_fd, s0.cv_j_kg_k),
                    "dp/dT": rel(dpdt_fd, s0.dp_dt_rho_pa_k),
                    "dp/drho": rel(dpdrho_fd, s0.dp_drho_t_pa_m3_kg),
                }
                for name, err in checks.items():
                    if err > 2e-7:
                        failures.append((fluid, t, p, name, err))
        self.assertEqual(failures, [])


class ExactIdentities(unittest.TestCase):
    def test_mayer_sound_speed_gibbs_and_joule_thomson(self) -> None:
        failures = []
        for fluid in REFERENCE_FLUID_KEYS:
            f = reference_fluid(fluid)
            for t, p in random_states(fluid, 15):
                try:
                    s = f.state_tp(t, p)
                except RocketPropulsionError:
                    continue
                rho = s.density_kg_m3
                mayer = t * s.dp_dt_rho_pa_k ** 2 / (rho * rho * s.dp_drho_t_pa_m3_kg)
                sound = (s.cp_j_kg_k / s.cv_j_kg_k) * s.dp_drho_t_pa_m3_kg
                # (dv/dT)_p = (dp/dT)_rho / (rho^2 (dp/drho)_T)
                dvdt = s.dp_dt_rho_pa_k / (rho * rho * s.dp_drho_t_pa_m3_kg)
                jt = (t * dvdt - 1 / rho) / s.cp_j_kg_k
                checks = {
                    "mayer": rel(s.cp_j_kg_k - s.cv_j_kg_k, mayer),
                    "sound": rel(s.speed_of_sound_m_s ** 2, sound),
                    "gibbs": abs(s.gibbs_energy_j_kg - (s.enthalpy_j_kg - t * s.entropy_j_kg_k))
                    / max(abs(s.enthalpy_j_kg), t * abs(s.entropy_j_kg_k), 1.0),
                    "jt": abs(s.joule_thomson_k_pa - jt) / max(abs(jt), 1e-9),
                    "z": rel(s.compressibility_factor,
                             p / (rho * reference_fluid_model(fluid).gas_constant_j_mol_k * t
                                  / reference_fluid_model(fluid).molar_mass_kg_mol)),
                }
                for name, err in checks.items():
                    if err > 1e-9:
                        failures.append((fluid, t, p, name, err))
        self.assertEqual(failures, [])

    def test_state_tp_inverts_state_rho_t(self) -> None:
        for fluid in REFERENCE_FLUID_KEYS:
            f = reference_fluid(fluid)
            for t, p in random_states(fluid, 8):
                try:
                    s = f.state_tp(t, p)
                except RocketPropulsionError:
                    continue
                back = f.state_rho_t(s.density_kg_m3, t)
                with self.subTest(fluid=fluid, t=t, p=p):
                    self.assertLess(rel(back.pressure_pa, p), 1e-10)


class SaturationLine(unittest.TestCase):
    def test_phase_equilibrium_and_clausius_clapeyron(self) -> None:
        for fluid in REFERENCE_FLUID_KEYS:
            model = reference_fluid_model(fluid)
            f = reference_fluid(fluid)
            for frac in (0.6, 0.75, 0.9, 0.97):
                t = model.t_min_k + frac * (model.t_crit_k - model.t_min_k)
                try:
                    sat = f.saturation_t(t)
                except RocketPropulsionError:
                    continue
                dt = 1e-4 * t
                hi, lo = f.saturation_t(t + dt), f.saturation_t(t - dt)
                dpdt = (hi.pressure_pa - lo.pressure_pa) / (2 * dt)
                dv = 1 / sat.vapor.density_kg_m3 - 1 / sat.liquid.density_kg_m3
                clapeyron = sat.enthalpy_of_vaporization_j_kg / (t * dv)
                with self.subTest(fluid=fluid, t=t):
                    self.assertLess(abs(sat.liquid.gibbs_energy_j_kg - sat.vapor.gibbs_energy_j_kg)
                                    / abs(sat.enthalpy_of_vaporization_j_kg), 1e-10)
                    self.assertLess(rel(sat.liquid.pressure_pa, sat.vapor.pressure_pa), 1e-9)
                    self.assertLess(rel(dpdt, clapeyron), 2e-6)
                    self.assertGreater(sat.liquid.density_kg_m3, sat.vapor.density_kg_m3)
                    self.assertGreater(sat.enthalpy_of_vaporization_j_kg, 0.0)
