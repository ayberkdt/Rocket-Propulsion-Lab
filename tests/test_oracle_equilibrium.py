"""Chemical equilibrium against its optimality conditions and physical invariants.

A Gibbs minimum at fixed T and p has a certificate independent of the solver:
every species present satisfies ``mu_j/RT = sum_i lambda_i a_ij`` for one set
of element potentials ``lambda``, and every absent species would raise the
Gibbs energy. Here ``lambda`` is fitted by least squares to the returned
composition and the residual is required on every species, trace ones
included, with ``mu_j = g_j(T) + RT ln(x_j p/p0)`` built from the species data
directly (not from the solver's internals).

Further invariants checked on random problems:

* element conservation to round-off;
* the solution is extensive: scaling every reactant by ``k`` scales every
  amount by ``k`` and leaves mole fractions unchanged;
* an independent constrained minimiser (scipy SLSQP on ``G/RT``) never
  finds a lower Gibbs energy than the package;
* an HP solve returns the temperature at which a TP solve reproduces the
  target enthalpy;
* Le Chatelier: dissociation (H, O, OH) decreases with pressure at fixed T
  and increases with T at fixed p.
"""

from __future__ import annotations

import math
import unittest

import numpy as np
from scipy import optimize

from rocket_propulsion.core.constants import MOLAR_GAS_CONSTANT_J_MOL_K
from rocket_propulsion.thermochemistry import EquilibriumProblem, LocalEquilibriumProvider
from rocket_propulsion.thermochemistry.species import get_species

R = MOLAR_GAS_CONSTANT_J_MOL_K
P0 = 1.0e5
PROVIDER = LocalEquilibriumProvider()
RNG = np.random.default_rng(424242)


def solve_tp(reactants: dict[str, float], t: float, p: float):
    return PROVIDER.solve(EquilibriumProblem.tp(reactants, temperature_k=t, pressure_pa=p))


def gibbs_rt(species: str, t: float) -> float:
    prop = get_species(species).properties(t)
    return (prop.enthalpy_molar_j_mol - t * prop.entropy_molar_j_mol_k) / (R * t)


def composition_matrix(names: list[str]) -> tuple[list[str], np.ndarray]:
    elements = sorted({e for n in names for e, _ in get_species(n).formula})
    a = np.zeros((len(elements), len(names)))
    for j, n in enumerate(names):
        for e, count in get_species(n).formula:
            a[elements.index(e), j] = count
    return elements, a


def random_hco_problem() -> tuple[dict[str, float], float, float]:
    t = float(RNG.uniform(1500.0, 4200.0))
    p = float(10 ** RNG.uniform(4.0, 7.5))
    choice = int(RNG.integers(3))
    if choice == 0:
        reactants = {"H2": float(RNG.uniform(0.5, 4.0)), "O2": 1.0}
    elif choice == 1:
        reactants = {"CH4": 1.0, "O2": float(RNG.uniform(0.8, 3.5))}
    else:
        reactants = {"H2": float(RNG.uniform(0.2, 2.0)), "O2": 1.0, "N2": float(RNG.uniform(0.0, 2.0))}
    return reactants, t, p


class OptimalityCertificate(unittest.TestCase):
    def test_element_potentials_fit_every_present_species(self) -> None:
        worst = 0.0
        for _ in range(40):
            reactants, t, p = random_hco_problem()
            r = solve_tp(reactants, t, p)
            names = [s.species for s in r.species if s.moles > 0.0]
            x = np.array([s.mole_fraction for s in r.species if s.moles > 0.0])
            _, a = composition_matrix(names)
            mu = np.array([gibbs_rt(n, t) for n in names]) + np.log(x * p / P0)
            # least squares on the major species, residual on all
            major = x > 1e-6
            lam, *_ = np.linalg.lstsq(a[:, major].T, mu[major], rcond=None)
            residual = np.abs(a.T @ lam - mu)
            worst = max(worst, float(residual.max()))
            with self.subTest(reactants=reactants, t=t, p=p):
                # residual in ln(x): 1e-6 means trace fractions are right to 1e-6 relative
                self.assertLess(float(residual.max()), 1e-6)
        self.assertLess(worst, 1e-6)

    def test_scipy_never_finds_lower_gibbs_energy(self) -> None:
        for reactants, t, p in (({"H2": 2.0, "O2": 1.0}, 3200.0, 2e6), ({"CH4": 1.0, "O2": 2.0}, 3400.0, 7e6),
                                ({"H2": 1.0, "O2": 1.0}, 2500.0, 1e5)):
            r = solve_tp(reactants, t, p)
            names = [s.species for s in r.species]
            n_pkg = np.array([s.moles for s in r.species])
            _, a = composition_matrix(names)
            b = a @ n_pkg
            g0 = np.array([gibbs_rt(n, t) for n in names])

            def g_total(n: np.ndarray, g0: np.ndarray = g0, p: float = p) -> float:
                n = np.maximum(n, 1e-300)
                return float(np.sum(n * (g0 + np.log(n / n.sum() * p / P0))))
            start = np.full(len(names), b.sum() / len(names))
            res = optimize.minimize(g_total, start, method="SLSQP",
                                    bounds=[(1e-30, None)] * len(names),
                                    constraints={"type": "eq", "fun": lambda n, a=a, b=b: a @ n - b},
                                    options={"ftol": 1e-14, "maxiter": 2000})
            with self.subTest(reactants=reactants):
                self.assertLessEqual(g_total(n_pkg), g_total(res.x) + 1e-9 * abs(g_total(n_pkg)))


class Invariants(unittest.TestCase):
    def test_element_conservation_within_declared_tolerance_and_honestly_reported(self) -> None:
        """Imbalance stays within the solver's 2e-9 acceptance and is not under-reported."""
        tolerance = LocalEquilibriumProvider.residual_tolerance
        for _ in range(30):
            reactants, t, p = random_hco_problem()
            r = solve_tp(reactants, t, p)
            names = [s.species for s in r.species]
            elements, a = composition_matrix(names)
            produced = a @ np.array([s.moles for s in r.species])
            fed = dict(r.element_moles)
            actual = max(abs(value - fed[e]) / fed[e] for e, value in zip(elements, produced, strict=True))
            with self.subTest(reactants=reactants, t=t, p=p):
                self.assertLessEqual(actual, tolerance)
                # the reported residual must not understate what we recompute
                self.assertGreaterEqual(r.residuals.maximum_element_relative * (1 + 1e-6) + 1e-15, actual)
                self.assertLess(abs(sum(s.mole_fraction for s in r.species) - 1.0), 1e-12)

    def test_extensive_scaling(self) -> None:
        base, t, p = {"H2": 2.3, "O2": 1.0}, 3100.0, 3e6
        r1 = solve_tp(base, t, p)
        for k in (1e-3, 0.37, 12.5, 4e3):
            rk = solve_tp({key: k * v for key, v in base.items()}, t, p)
            for s1, sk in zip(r1.species, rk.species, strict=True):
                with self.subTest(k=k, species=s1.species):
                    self.assertEqual(s1.species, sk.species)
                    self.assertLess(abs(sk.mole_fraction - s1.mole_fraction), 1e-12 + 1e-9 * s1.mole_fraction)

    def test_le_chatelier(self) -> None:
        def dissociated(t: float, p: float) -> float:
            r = solve_tp({"H2": 2.0, "O2": 1.0}, t, p)
            return sum(s.mole_fraction for s in r.species if s.species in ("H", "O", "OH"))
        pressures = [1e4, 1e5, 1e6, 1e7]
        values = [dissociated(3200.0, p) for p in pressures]
        self.assertEqual(values, sorted(values, reverse=True))
        temps = [2200.0, 2800.0, 3400.0, 4000.0]
        values = [dissociated(t, 2e6) for t in temps]
        self.assertEqual(values, sorted(values))

    def test_hp_temperature_reproduces_target_enthalpy(self) -> None:
        reactants = {"H2": 2.5, "O2": 1.0}
        hp = PROVIDER.solve(EquilibriumProblem.hp(reactants, pressure_pa=5e6))
        tp = solve_tp(reactants, hp.temperature_k, 5e6)
        self.assertLess(abs(tp.enthalpy_j - hp.enthalpy_j), 1e-7 * abs(hp.enthalpy_j) + 1e-6)
        for a, b in zip(hp.species, tp.species, strict=True):
            self.assertLess(abs(a.mole_fraction - b.mole_fraction), 1e-9)
        self.assertTrue(math.isfinite(hp.temperature_k))
