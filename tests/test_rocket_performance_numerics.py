"""Numerics behind the faster rocket-performance solver.

* ``_maximize_bounded`` is the bounded Brent method of SciPy ``fminbound``
  applied to ``-f``; on smooth, kinked, flat-topped and boundary maxima it must
  land where SciPy lands, within the shared tolerance.
* Isentrope warm starts reuse solved ``(ln p, ln T)`` points; a state reached
  after any history of earlier solves must equal a cold solve.
* The throat must be a mass-flux maximum: nearby pressures carry less flux.
"""

from __future__ import annotations

import math
import random
import unittest

from scipy import optimize

from rocket_propulsion.core.errors import ConvergenceError
from rocket_propulsion.thermochemistry.equilibrium import LocalEquilibriumProvider
from rocket_propulsion.thermochemistry.rocket_equilibrium import (
    get_propellant_pair,
    reactants_for_mass_ratio,
)
from rocket_propulsion.thermochemistry.rocket_performance import (
    _Expansion,
    _maximize_bounded,
    solve_rocket_performance,
)

TOLERANCE = 1e-5


class BoundedMaximizer(unittest.TestCase):
    def check(self, function, low: float, high: float) -> None:
        mine = _maximize_bounded(function, low, high, TOLERANCE)
        reference = optimize.fminbound(lambda x: -function(x), low, high, xtol=TOLERANCE)
        self.assertLessEqual(low, mine)
        self.assertLessEqual(mine, high)
        self.assertLess(abs(mine - reference), 4 * TOLERANCE)
        # the value is as good as SciPy's up to the curvature over the tolerance
        self.assertGreaterEqual(function(mine), function(reference) - 1e-8 * max(1.0, abs(function(reference))))

    def test_smooth_kinked_and_boundary_maxima(self) -> None:
        rng = random.Random(7)
        for _ in range(200):
            centre, scale = rng.uniform(-3.0, 3.0), rng.uniform(0.05, 20.0)
            low, high = centre - rng.uniform(0.1, 4.0), centre + rng.uniform(0.1, 4.0)
            with self.subTest(centre=centre, scale=scale, low=low, high=high):
                self.check(lambda x, c=centre, k=scale: -k * (x - c) ** 2 + math.sin(x) * 1e-3, low, high)
                self.check(lambda x, c=centre: -abs(x - c), low, high)
        # maximum on either bound
        self.check(lambda x: x, 0.0, 1.0)
        self.check(lambda x: -x, 0.0, 1.0)
        # quartic: flat top, parabolic steps are poor
        self.check(lambda x: -(x - 0.3) ** 4, -1.0, 2.0)

    def test_evaluation_count_beats_golden_section(self) -> None:
        calls = 0

        def counted(x: float) -> float:
            nonlocal calls
            calls += 1
            return -((x - 0.42) ** 2)
        _maximize_bounded(counted, 0.0, 1.0, TOLERANCE)
        # golden-section to the same tolerance needs ceil(ln(1e-5)/ln(0.618)) + 2 = 26 calls
        self.assertLess(calls, 13)

    def test_non_finite_values_cannot_hang(self) -> None:
        for bad in (math.nan, math.inf):
            with self.subTest(value=bad), self.assertRaises(ConvergenceError):
                _maximize_bounded(lambda x, bad=bad: bad, 0.0, 1.0, TOLERANCE)
        # a NaN that appears only partway through the search is caught too
        with self.assertRaises(ConvergenceError):
            _maximize_bounded(lambda x: math.nan if x > 0.6 else -((x - 0.9) ** 2), 0.0, 1.0, TOLERANCE)


class WarmStartedIsentrope(unittest.TestCase):
    def test_states_do_not_depend_on_solve_history(self) -> None:
        reactants = reactants_for_mass_ratio(get_propellant_pair("LOX_CH4"), 3.2)
        provider = LocalEquilibriumProvider()
        pc = 8.0e6
        rng = random.Random(3)
        pressures = [pc * 10 ** rng.uniform(-3.0, -0.1) for _ in range(12)]
        for chemistry in ("equilibrium", "frozen"):
            warm = _Expansion(reactants, provider, chemistry, pc)
            warm_states = [warm.isentrope(p) for p in pressures]
            for p, state in zip(pressures, warm_states, strict=True):
                cold = _Expansion(reactants, provider, chemistry, pc).isentrope(p)
                with self.subTest(chemistry=chemistry, pressure=p):
                    self.assertLess(abs(state.temperature_k - cold.temperature_k), 1e-8 * cold.temperature_k)
                    self.assertLess(abs(state.entropy_j_kg_k - warm.chamber.entropy_j_kg_k),
                                    1e-10 * abs(warm.chamber.entropy_j_kg_k))


class ThroatIsAMaximum(unittest.TestCase):
    def test_neighbouring_pressures_carry_less_flux(self) -> None:
        for key, ratio in (("LOX_LH2", 6.0), ("LOX_RP1", 2.6)):
            reactants = reactants_for_mass_ratio(get_propellant_pair(key), ratio)
            pc = 6.894757e6
            result = solve_rocket_performance(reactants=reactants, chamber_pressure_pa=pc, area_ratios=(5.0,))
            expansion = _Expansion(reactants, LocalEquilibriumProvider(), "equilibrium", pc)
            pt = result.throat.pressure_pa
            peak = result.throat.density_kg_m3 * result.throat.velocity_m_s
            with self.subTest(pair=key):
                for factor in (0.97, 0.99, 1.01, 1.03):
                    self.assertLess(expansion.mass_flux(pt * factor)[0], peak)
                # equilibrium p*/pc lies between the frozen-gamma limits 0.5 and 0.6
                self.assertTrue(0.5 < pt / pc < 0.6, pt / pc)


if __name__ == "__main__":
    unittest.main()
