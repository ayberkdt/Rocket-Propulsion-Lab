"""Finite-rate kinetics against an independent stiff integrator and thermodynamics.

The H2-O2 GRI-Mech 3.0 subset is integrated at constant temperature and
volume and checked in four independent ways:

* the package's RODAS integrator and scipy's Radau IIA (a different implicit
  Runge-Kutta family) agree on the reactor history to the tolerance of both;
* hydrogen and oxygen atoms are conserved along the whole trajectory;
* the long-time state satisfies detailed balance, every reaction's forward
  rate equals its reverse rate, and matches the Gibbs-minimisation solver's
  TP equilibrium at the same temperature and final pressure (two codes that
  share only species data);
* the equilibrium constants equal ``exp(-sum nu g/RT) (p0/RT)^dnu`` rebuilt
  from the species data, and the analytic Jacobian equals a fourth-order
  finite difference of the production rates;
* a rate multiplier changes the path but never the end state.
"""

from __future__ import annotations

import math
import os
import unittest

import numpy as np
from scipy import integrate

from rocket_propulsion.core.constants import MOLAR_GAS_CONSTANT_J_MOL_K
from rocket_propulsion.thermochemistry import EquilibriumProblem, LocalEquilibriumProvider
from rocket_propulsion.thermochemistry import kinetics as k
from rocket_propulsion.thermochemistry.kinetics_mechanisms import gri30_h2o2_mechanism
from rocket_propulsion.thermochemistry.species import get_species

MECH = gri30_h2o2_mechanism()
SPECIES = MECH.species
R = MOLAR_GAS_CONSTANT_J_MOL_K
P0 = 1.0e5
T_REACTOR = 2800.0
P_START = 2.0e6


def omega(c: np.ndarray | list[float], multiplier: float = 1.0) -> np.ndarray:
    return np.array(k.evaluate_rates(MECH, T_REACTOR, list(np.maximum(c, 0.0)),
                                     rate_multiplier=multiplier).net_production_rates)


def start() -> list[float]:
    return list(k.concentrations_from_mole_fractions(MECH, T_REACTOR, P_START, {"H2": 0.6, "O2": 0.4}))


def atoms(c: np.ndarray | list[float]) -> dict[str, float]:
    out: dict[str, float] = {}
    for name, ci in zip(SPECIES, c, strict=True):
        for element, count in get_species(name).formula:
            out[element] = out.get(element, 0.0) + count * ci
    return out


class ReactorAgainstRadau(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.times = [1e-8, 1e-7, 1e-6, 1e-5, 1e-4, 1e-3]
        cls.mine = k.solve_stiff(lambda t, y: list(omega(y)), (0.0, 1e-3), start(), rtol=1e-8, atol=1e-12,
                                 jacobian=lambda t, y: k.production_rate_jacobian(MECH, T_REACTOR, list(np.maximum(y, 0.0))),
                                 autonomous=True, t_eval=cls.times)
        cls.ref = integrate.solve_ivp(lambda t, y: omega(y), (0.0, 1e-3), start(), method="Radau",
                                      rtol=1e-10, atol=1e-14, t_eval=cls.times)

    def test_histories_agree(self) -> None:
        mine = np.array(self.mine.y)
        ref = self.ref.y.T
        scale = np.max(np.abs(ref), axis=0) + 1e-12
        for i, t in enumerate(self.times):
            with self.subTest(t=t):
                self.assertLess(float(np.max(np.abs(mine[i] - ref[i]) / scale)), 2e-6)

    def test_atoms_conserved_along_trajectory(self) -> None:
        a0 = atoms(start())
        for t, y in zip(self.times, self.mine.y, strict=True):
            a = atoms(y)
            for element in a0:
                with self.subTest(t=t, element=element):
                    self.assertLess(abs(a[element] - a0[element]) / a0[element], 1e-9)


@unittest.skipUnless(os.environ.get("RPL_SLOW_TESTS") == "1", "slow (~4 min): set RPL_SLOW_TESTS=1")
class EquilibriumLimit(unittest.TestCase):
    def test_detailed_balance_and_agreement_with_gibbs_solver(self) -> None:
        sol = k.solve_stiff(lambda t, y: list(omega(y)), (0.0, 1.0), start(), rtol=1e-10, atol=1e-14,
                            jacobian=lambda t, y: k.production_rate_jacobian(MECH, T_REACTOR, list(np.maximum(y, 0.0))),
                            autonomous=True, t_eval=[1.0])
        c = np.maximum(np.array(sol.y[-1]), 0.0)
        rates = k.evaluate_rates(MECH, T_REACTOR, list(c))
        for j, (fwd, rev) in enumerate(zip(rates.forward_rates, rates.reverse_rates, strict=True)):
            with self.subTest(reaction=MECH.reactions[j].equation):
                self.assertLess(abs(fwd - rev), 1e-6 * max(fwd, rev, 1e-30) + 1e-18)
        total = float(c.sum())
        pressure = total * R * T_REACTOR
        x_kin = c / total
        moles = atoms(c)
        eq = LocalEquilibriumProvider().solve(EquilibriumProblem.tp(
            {"H2": moles["H"] / 2.0, "O2": moles["O"] / 2.0}, temperature_k=T_REACTOR, pressure_pa=pressure,
            candidate_species=SPECIES))
        x_eq = {s.species: s.mole_fraction for s in eq.species}
        for name, x in zip(SPECIES, x_kin, strict=True):
            with self.subTest(species=name):
                self.assertLess(abs(x - x_eq.get(name, 0.0)), 1e-6 * max(x, 1e-6))

    def test_rate_multiplier_changes_path_not_end_state(self) -> None:
        ends = []
        for mult in (0.1, 1.0, 10.0):
            sol = k.solve_stiff(lambda t, y, m=mult: list(omega(y, m)), (0.0, 5.0), start(), rtol=1e-10,
                                atol=1e-14, autonomous=True, t_eval=[5.0])
            ends.append(np.array(sol.y[-1]))
        for end in ends[1:]:
            self.assertLess(float(np.max(np.abs(end - ends[0]) / (np.abs(ends[0]) + 1e-9))), 1e-5)


class ThermodynamicConsistency(unittest.TestCase):
    def test_equilibrium_constants_from_species_data(self) -> None:
        for t in (1000.0, 2000.0, 3000.0, 4000.0):
            kc = k.equilibrium_constants(MECH, t)
            for j, reaction in enumerate(MECH.reactions):
                net: dict[str, int] = {}
                for s, n in reaction.products:
                    net[s] = net.get(s, 0) + n
                for s, n in reaction.reactants:
                    net[s] = net.get(s, 0) - n
                dg = 0.0
                for s, n in net.items():
                    prop = get_species(s).properties(t)
                    dg += n * (prop.enthalpy_molar_j_mol - t * prop.entropy_molar_j_mol_k) / (R * t)
                dnu = sum(net.values())
                expected = math.exp(-dg) * (P0 / (R * t)) ** dnu
                with self.subTest(t=t, reaction=reaction.equation):
                    self.assertLess(abs(kc[j] / expected - 1.0), 1e-10)

    def test_analytic_jacobian_against_fourth_order_differences(self) -> None:
        c = np.array(start()) * 0.7 + 1e-3
        jac = np.array(k.production_rate_jacobian(MECH, T_REACTOR, list(c)))
        for i in range(len(SPECIES)):
            h = 1e-3 * max(c[i], 1e-3)

            def f(e: float, i: int = i) -> np.ndarray:
                cc = c.copy()
                cc[i] += e
                return omega(cc)
            fd = (8 * (f(h / 2) - f(-h / 2)) - (f(h) - f(-h))) / (6 * h)
            scale = np.max(np.abs(jac)) + 1e-30
            with self.subTest(column=SPECIES[i]):
                self.assertLess(float(np.max(np.abs(jac[:, i] - fd))) / scale, 1e-8)
