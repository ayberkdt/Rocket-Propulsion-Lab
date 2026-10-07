"""Try to break the Lunaris adapter: odd timings, dry-mass edge cases, oracles.

Every expected value here is computed by code written in this file, not by
the adapter:

* Mass history: piecewise-linear flow integrates exactly to
  m(t) = m0 - Σ [q0 τ + (q1 - q0) τ² / (2 d)] (trapezoid areas; NASA Glenn,
  Ideal Rocket Equation, for the mass-flow/thrust relations,
  https://www1.grc.nasa.gov/beginners-guide-to-aeronautics/ideal-rocket-equation/).
  The dry-mass crossing is found by plain bisection on that function.
* Free-space Δv with a constant thrust/flow ratio c: Δv = c ln(m0/mf) (same
  source).
* Lunar point-mass gravity a = -μ r / |r|³ with the Lunaris constants
  ``MU_MOON`` and ``R_MOON``; the oracle ODE is integrated with SciPy's DOP853
  (Hairer, Nørsett and Wanner, Solving Ordinary Differential Equations I,
  2nd ed., Springer, 1993, bibliographic citation; SciPy ``solve_ivp``,
  https://docs.scipy.org/doc/scipy/reference/generated/scipy.integrate.solve_ivp.html),
  restarted at every discontinuity the test computes itself.
* Restart and one-sided-limit semantics: Orekit 13.1.5 EventDetector API,
  https://www.orekit.org/static/apidocs/org/orekit/propagation/events/EventDetector.html
"""

import importlib.util
import random
import unittest
from itertools import pairwise
from math import log, sqrt

from rocket_propulsion.core.errors import DomainError, InputError
from rocket_propulsion.integrations import (
    LunarisFiniteBurnDynamics,
    lunaris_burn_windows,
    lunaris_impulsive_equivalents,
)
from rocket_propulsion.propulsion.burns import (
    EstimableParameter,
    PropagatorPropulsionBridge,
    PropulsionTimeline,
    ScheduledBurn,
    TabulatedBurnArtifact,
    TabulatedBurnSegment,
)
from tests.test_lunaris_integration import _FreeSpaceDynamics

LUNARIS = importlib.util.find_spec("lunaris") is not None
SCIPY = importlib.util.find_spec("scipy") is not None
X = lambda burn_id, t, y: (1.0, 0.0, 0.0)


def bridge_from(rows, dry=None, tag="adv", bias_limit_s=None):
    """rows: (duration, f0, f1, q0, q1) per segment, contiguous from zero."""

    segments, time = [], 0.0
    for duration, f0, f1, q0, q1 in rows:
        segments.append(TabulatedBurnSegment(time, time + duration, f0, f1, f0, f1, q0, q1))
        time += duration
    artifact = TabulatedBurnArtifact(segments=tuple(segments), source_id=tag,
                                     source_sha256="d" * 64)
    extra = {}
    if bias_limit_s is not None:
        extra["parameters"] = (
            EstimableParameter("thrust_scale", 1.0, 1e-3, 0.01, 2.0, "1"),
            EstimableParameter("mass_flow_scale", 1.0, 1e-3, 0.01, 2.0, "1"),
            EstimableParameter("ignition_time_bias_s", 0.0, 1e-3, -bias_limit_s, bias_limit_s,
                               "s"),
        )
    return PropagatorPropulsionBridge(artifact, require_complete_stream_binding=False,
                                      protected_dry_mass_kg=dry, **extra)


def consumed_oracle(rows, scale, tau):
    used, time = 0.0, 0.0
    for duration, _, _, q0, q1 in rows:
        if tau <= time:
            break
        local = min(tau, time + duration) - time
        used += q0 * local + (q1 - q0) * local * local / (2.0 * duration)
        time += duration
    return used * scale


def first_crossing_oracle(rows, scale, budget):
    total = sum(d for d, *_ in rows)
    if consumed_oracle(rows, scale, total) < budget:
        return None
    low, high = 0.0, total
    for _ in range(200):
        middle = 0.5 * (low + high)
        if consumed_oracle(rows, scale, middle) >= budget:
            high = middle
        else:
            low = middle
    return high


def thrust_at(timeline, t, side, mass, overrides=None):
    value = timeline.evaluate(time_s=t, vehicle_mass_kg=mass,
                              directions={b.burn_id: (1.0, 0.0, 0.0) for b in timeline.burns},
                              side=side, parameter_overrides=overrides,
                              apply_inventory_limits=False)
    return value.acceleration_m_s2[0] * mass


class ReviewRegressionTests(unittest.TestCase):
    def test_ignition_inside_the_previous_burn_by_half_a_nanosecond_is_refused(self) -> None:
        timeline = PropulsionTimeline(burns=(
            ScheduledBurn("A", bridge_from(((60.0, 400, 400, 0.14, 0.14),), tag="a"), 100.0),
            ScheduledBurn("B", bridge_from(((30.0, 300, 300, 0.1, 0.1),), tag="b"), 160.0),
        ))
        overrides = {"B": {"ignition_time_bias_s": -5e-10}}
        for call in (
            lambda: lunaris_burn_windows(timeline, parameter_overrides=overrides),
            lambda: LunarisFiniteBurnDynamics(_FreeSpaceDynamics(), timeline, X,
                                              parameter_overrides=overrides),
            lambda: lunaris_impulsive_equivalents(
                timeline, initial_mass_kg=1000.0, parameter_overrides=overrides,
                directions={"A": ((1, 0, 0), "inertial"), "B": ((1, 0, 0), "inertial")}),
        ):
            with self.assertRaises(DomainError):
                call()

    def test_back_to_back_burns_keep_each_burns_own_limit(self) -> None:
        timeline = PropulsionTimeline(burns=(
            ScheduledBurn("A", bridge_from(((60.0, 400, 400, 0.14, 0.14),), tag="a"), 100.0),
            ScheduledBurn("B", bridge_from(((30.0, 250, 250, 0.1, 0.1),), tag="b"), 160.0),
        ))
        dynamics = LunarisFiniteBurnDynamics(_FreeSpaceDynamics(), timeline, X)
        self.assertEqual(dynamics.segmentation_times, (100.0, 160.0, 190.0))
        y = [0, 0, 0, 0, 0, 0, 1000.0]
        rhs = dynamics.build_rhs()
        rhs(100.0, y)  # segment start
        self.assertAlmostEqual(rhs(160.0, y)[3] * 1000.0, 400.0)  # A's left limit
        rhs = dynamics.build_rhs()
        self.assertAlmostEqual(rhs(160.0, y)[3] * 1000.0, 250.0)  # B's right limit

    def test_impulsive_path_validates_overrides_like_the_finite_path(self) -> None:
        timeline = PropulsionTimeline(burns=(
            ScheduledBurn("A", bridge_from(((30.0, 400, 400, 0.14, 0.14),), tag="a"), 100.0),
            ScheduledBurn("B", bridge_from(((30.0, 400, 400, 0.14, 0.14),), tag="b"), 140.0),
        ))
        directions = {"A": ((1, 0, 0), "inertial"), "B": ((1, 0, 0), "inertial")}
        bad = (
            {"A": {"ignition_time_bias_s": 15.0}},       # overlaps B after the bias
            {"ghost": {"thrust_scale": 1.0}},            # unknown burn
            5, "x", [("A", {})], {"A": 5}, {"A": {"thrust_scale": "1"}},
        )
        for overrides in bad:
            with self.assertRaises((DomainError, InputError), msg=repr(overrides)):
                lunaris_impulsive_equivalents(timeline, initial_mass_kg=1000.0,
                                              directions=directions,
                                              parameter_overrides=overrides)

    def test_a_bias_that_reorders_burns_chains_masses_in_the_new_order(self) -> None:
        rows_a = ((5.0, 400, 400, 0.8, 0.8),)
        rows_b = ((10.0, 300, 300, 0.1, 0.1),)
        timeline = PropulsionTimeline(burns=(
            ScheduledBurn("A", bridge_from(rows_a, tag="a", bias_limit_s=60.0), 100.0),
            ScheduledBurn("B", bridge_from(rows_b, tag="b"), 106.0),
        ))
        # A bridge with a +-60 s timing driver; a +18 s bias moves A after B.
        impulses = lunaris_impulsive_equivalents(
            timeline, initial_mass_kg=500.0,
            directions={"A": ((1, 0, 0), "inertial"), "B": ((1, 0, 0), "inertial")},
            parameter_overrides={"A": {"ignition_time_bias_s": 18.0}},
        )
        self.assertEqual([item.burn_id for item in impulses], ["B", "A"])
        self.assertAlmostEqual(impulses[0].final_mass_kg, 500.0 - 1.0, places=12)
        self.assertEqual(impulses[1].initial_mass_kg, impulses[0].final_mass_kg)
        self.assertAlmostEqual(impulses[1].delta_v_m_s, 500.0 * log(499.0 / 495.0), places=9)

    @unittest.skipUnless(LUNARIS, "lunaris is not installed")
    def test_partial_window_propagation_uses_only_markers_inside_it(self) -> None:
        from lunaris.common.type_defs import EventConfig, PropagatorConfig, TimeConfig
        from lunaris.core.propagation import propagate

        timeline = PropulsionTimeline(burns=(
            ScheduledBurn("A", bridge_from(((60.0, 400, 400, 0.14, 0.14),), tag="a"), 100.0),
            ScheduledBurn("B", bridge_from(((40.0, 400, 400, 0.14, 0.14),), tag="b"), 3600.0),
        ))
        dynamics = LunarisFiniteBurnDynamics(_FreeSpaceDynamics(), timeline, X)
        cfg = PropagatorConfig(method="DOP853", rtol=1e-12, atol=1e-13, verbose=False,
                               use_nyquist_max_step=False,
                               events=EventConfig(detect_impact=False,
                                                  enable_peri_apo_events=False))
        with self.assertRaises(ValueError):  # Lunaris refuses out-of-window markers
            propagate(dynamics, [2e6, 0, 0, 0, 1, 0, 1000.0], cfg,
                      time_cfg=TimeConfig(duration_s=1000.0, output_dt_s=100.0),
                      maneuver_plan=dynamics.segmentation_plan())
        result = propagate(dynamics, [2e6, 0, 0, 0, 1, 0, 1000.0], cfg,
                           time_cfg=TimeConfig(duration_s=1000.0, output_dt_s=100.0),
                           maneuver_plan=dynamics.segmentation_plan(0.0, 1000.0))
        self.assertAlmostEqual(float(result.y[-1][6]), 1000.0 - 8.4, places=9)
        with self.assertRaises(DomainError):
            dynamics.segmentation_plan(10.0, 5.0)


class DryMassBatteryTests(unittest.TestCase):
    def random_rows(self, generator):
        rows = []
        for _ in range(generator.randint(1, 5)):
            kind = generator.choice(("up", "down", "flat", "zero", "random"))
            q0, q1 = {
                "up": (0.0, generator.uniform(0.05, 0.5)),
                "down": (generator.uniform(0.05, 0.5), generator.uniform(0.0, 0.04)),
                "flat": (0.2, 0.2),
                "zero": (0.0, 0.0),
                "random": (generator.uniform(0, 0.5), generator.uniform(0, 0.5)),
            }[kind]
            f0, f1 = 2000.0 * q0 + 1.0, 2000.0 * q1 + 1.0
            rows.append((generator.uniform(0.5, 30.0), f0, f1, q0, q1))
        if all(q0 == q1 == 0.0 for *_, q0, q1 in rows):
            rows.append((5.0, 400.0, 400.0, 0.2, 0.2))
        return rows

    def test_cut_time_is_the_first_crossing_found_by_an_independent_oracle(self) -> None:
        generator = random.Random(424242)
        for case in range(400):
            rows = self.random_rows(generator)
            scale = generator.uniform(0.7, 1.3)
            total = consumed_oracle(rows, scale, sum(d for d, *_ in rows))
            mode = generator.choice(("inside", "boundary", "never", "at_stop"))
            if mode == "inside":
                budget = generator.uniform(0.01, 0.99) * total
            elif mode == "boundary":
                k = generator.randint(1, len(rows))
                budget = consumed_oracle(rows, scale, sum(d for d, *_ in rows[:k]))
            elif mode == "never":
                budget = total * 1.5
            else:
                budget = total
            if budget <= 0.0:
                continue
            m0 = 100.0
            dry = m0 - budget
            ignition = generator.uniform(0.0, 1e6)
            timeline = PropulsionTimeline(
                burns=(ScheduledBurn("b", bridge_from(rows, dry=dry, tag=f"c{case}"), ignition),)
            )
            (window,) = lunaris_burn_windows(
                timeline, initial_mass_kg=m0,
                parameter_overrides={"b": {"mass_flow_scale": scale}},
            )
            # Contract: the cut is the first time the mass is within 1e-12 kg of
            # the protected dry mass (the bridge clamp's margin); a cut within
            # the profile's boundary tolerance of the stop is the stop.
            expected = first_crossing_oracle(rows, scale, budget - 1e-12)
            duration = sum(d for d, *_ in rows)
            if expected is None or expected >= duration * (1 - 1e-12):
                # Reaching dry mass only at (or after) the stop is not a cut.
                self.assertIsNone(window.depletion_time_s, (case, mode))
                continue
            self.assertIsNotNone(window.depletion_time_s, (case, mode))
            relative = window.depletion_time_s - ignition
            # A plateau (zero flow) right after the crossing makes every point
            # of it a crossing; the cut must be the *first* one.
            self.assertLessEqual(relative, expected + 1e-9 * duration + 1e-9 * ignition,
                                 (case, mode))
            mass_at_cut = m0 - consumed_oracle(rows, scale, relative)
            self.assertLess(abs(mass_at_cut - dry), 1e-11 * m0, (case, mode))
            # Nothing between ignition and the cut is already at the limit.
            for k in range(1, 50):
                earlier = relative * k / 50 - 1e-9 * duration
                if earlier > 0:
                    self.assertGreater(m0 - consumed_oracle(rows, scale, earlier), dry + 1e-12,
                                       (case, mode))

    def test_vehicle_already_at_dry_mass_never_fires(self) -> None:
        timeline = PropulsionTimeline(burns=(
            ScheduledBurn("b", bridge_from(((10.0, 400, 400, 0.2, 0.2),), dry=100.0), 50.0),
        ))
        (window,) = lunaris_burn_windows(timeline, initial_mass_kg=100.0)
        self.assertEqual(window.depletion_time_s, 50.0)
        dynamics = LunarisFiniteBurnDynamics(_FreeSpaceDynamics(), timeline, X,
                                             initial_mass_kg=100.0)
        rhs = dynamics.build_rhs()
        self.assertEqual(rhs(50.0, [0, 0, 0, 0, 0, 0, 100.0])[3], 0.0)
        self.assertEqual(rhs(55.0, [0, 0, 0, 0, 0, 0, 100.0])[3], 0.0)

    def test_mass_check_tolerance_is_exactly_what_is_documented(self) -> None:
        timeline = PropulsionTimeline(burns=(
            ScheduledBurn("b", bridge_from(((10.0, 400, 400, 0.2, 0.2),), dry=90.0), 50.0),
        ))
        dynamics = LunarisFiniteBurnDynamics(_FreeSpaceDynamics(), timeline, X,
                                             initial_mass_kg=100.0, mass_tolerance=1e-9)
        dynamics.build_rhs()(0.0, [0, 0, 0, 0, 0, 0, 100.0 * (1 + 5e-10)])
        with self.assertRaises(DomainError):
            dynamics.build_rhs()(0.0, [0, 0, 0, 0, 0, 0, 100.0 * (1 + 2e-9)])
        # Mid-burn segment starts are checked against the predicted mass too.
        predicted = dynamics.predicted_mass_kg(53.0)
        self.assertAlmostEqual(predicted, 100.0 - 0.6, places=12)
        with self.assertRaises(DomainError):
            dynamics.build_rhs()(53.0, [0, 0, 0, 0, 0, 0, predicted + 1e-6])

    def test_predicted_mass_is_monotone_and_matches_the_oracle(self) -> None:
        generator = random.Random(99)
        for _ in range(60):
            rows_a = self.random_rows(generator)
            rows_b = self.random_rows(generator)
            ign_a = generator.uniform(0, 100)
            ign_b = ign_a + sum(d for d, *_ in rows_a) + generator.uniform(0, 50)
            timeline = PropulsionTimeline(burns=(
                ScheduledBurn("A", bridge_from(rows_a, tag="a"), ign_a),
                ScheduledBurn("B", bridge_from(rows_b, tag="b"), ign_b),
            ))
            m0 = 1000.0
            dynamics = LunarisFiniteBurnDynamics(_FreeSpaceDynamics(), timeline, X,
                                                 initial_mass_kg=m0)
            end = ign_b + sum(d for d, *_ in rows_b) + 10
            previous = m0
            for k in range(301):
                t = end * k / 300
                mass = dynamics.predicted_mass_kg(t)
                self.assertLessEqual(mass, previous + 1e-12)
                previous = mass
                expected = (m0 - consumed_oracle(rows_a, 1.0, max(0.0, t - ign_a))
                            - consumed_oracle(rows_b, 1.0, max(0.0, t - ign_b)))
                self.assertAlmostEqual(mass, expected, delta=1e-11 * m0)


@unittest.skipUnless(LUNARIS and SCIPY, "needs lunaris and scipy")
class IndependentOracleTests(unittest.TestCase):
    """Lunaris + adapter against an ODE written and integrated here."""

    def test_biased_two_burn_sequence_with_dry_mass_cut_matches_scipy(self) -> None:
        import numpy as np
        from lunaris.common.constants import MU_MOON, R_MOON
        from lunaris.common.type_defs import EventConfig, PropagatorConfig, TimeConfig
        from lunaris.core.propagation import propagate
        from scipy.integrate import solve_ivp

        rows_a = ((2.0, 0.0, 450.0, 0.0, 0.15), (50.0, 450.0, 450.0, 0.15, 0.15),
                  (1.0, 300.0, 300.0, 0.1, 0.1))
        rows_b = ((400.0, 380.0, 380.0, 0.13, 0.13),)
        ign_a, ign_b, bias_a, m0, dry = 600.0, 3000.0, 0.73, 1000.0, 975.0
        timeline = PropulsionTimeline(burns=(
            ScheduledBurn("A", bridge_from(rows_a, tag="oa"), ign_a),
            ScheduledBurn("B", bridge_from(rows_b, dry=dry, tag="ob"), ign_b),
        ))
        overrides = {"A": {"ignition_time_bias_s": bias_a}}
        mu = float(MU_MOON)

        class Moon:
            grav = ephem = flags = None

            def build_rhs(self):
                def rhs(_t, y):
                    state = np.asarray(y, dtype=float)
                    out = np.zeros_like(state)
                    out[:3] = state[3:6]
                    out[3:6] = -mu * state[:3] / np.linalg.norm(state[:3]) ** 3
                    return out

                return rhs

        def prograde(burn_id, t, y):
            return (float(y[3]), float(y[4]), float(y[5]))

        dynamics = LunarisFiniteBurnDynamics(Moon(), timeline, prograde,
                                             parameter_overrides=overrides, initial_mass_kg=m0)
        radius = float(R_MOON) + 100_000.0
        y0 = [radius, 0.0, 0.0, 0.0, sqrt(mu / radius), 0.0, m0]
        cfg = PropagatorConfig(method="DOP853", rtol=1e-12, atol=1e-12, verbose=False,
                               use_nyquist_max_step=False,
                               events=EventConfig(detect_impact=False,
                                                  enable_peri_apo_events=False))
        duration = 4000.0
        result = propagate(dynamics, y0, cfg,
                           time_cfg=TimeConfig(duration_s=duration, output_dt_s=100.0),
                           maneuver_plan=dynamics.segmentation_plan(0.0, duration))

        # --- independent oracle -------------------------------------------
        def branch(rows, tau, tau_mid):
            """Linear law of the row that contains tau_mid, evaluated at tau.

            Choosing the row by the interval midpoint gives the one-sided
            limit at both interval ends, as a restarted integrator needs.
            """

            time = 0.0
            for d, f0, f1, q0, q1 in rows:
                if tau_mid < time + d:
                    w = (tau - time) / d
                    return f0 + (f1 - f0) * w, q0 + (q1 - q0) * w
                time += d
            raise AssertionError("midpoint outside the profile")

        consumed_a = consumed_oracle(rows_a, 1.0, 1e9)
        cut_b = first_crossing_oracle(rows_b, 1.0, m0 - consumed_a - dry)
        self.assertIsNotNone(cut_b)
        burns = ((ign_a + bias_a, rows_a, sum(d for d, *_ in rows_a)), (ign_b, rows_b, cut_b))
        edges = {0.0, duration}
        for burn_start, rows, stop in burns:
            elapsed = 0.0
            edges.add(burn_start)
            for d, *_ in rows:
                elapsed += d
                if elapsed < stop:
                    edges.add(burn_start + elapsed)
            edges.add(burn_start + stop)
        cuts = sorted(edge for edge in edges if edge <= duration)

        def oracle_rhs(t, state, mid):
            out = np.zeros(7)
            out[:3] = state[3:6]
            out[3:6] = -mu * state[:3] / np.linalg.norm(state[:3]) ** 3
            for burn_start, rows, stop in burns:
                if 0.0 <= mid - burn_start < stop:
                    force, flow = branch(rows, t - burn_start, mid - burn_start)
                    out[3:6] += force / state[6] * state[3:6] / np.linalg.norm(state[3:6])
                    out[6] -= flow
            return out

        y = np.array(y0, dtype=float)
        for start, end in pairwise(cuts):
            mid = 0.5 * (start + end)
            solution = solve_ivp(lambda t, s, mid=mid: oracle_rhs(t, s, mid), (start, end), y,
                                 method="DOP853", rtol=1e-12, atol=1e-12)
            y = solution.y[:, -1]

        final = np.asarray(result.y[-1], dtype=float)
        self.assertAlmostEqual(y[6], dry, delta=1e-9)          # oracle reached dry mass
        self.assertAlmostEqual(final[6], y[6], delta=1e-9)      # adapter agrees on mass
        self.assertLess(np.linalg.norm(final[:3] - y[:3]), 1e-3)     # 1 mm
        self.assertLess(np.linalg.norm(final[3:6] - y[3:6]), 1e-6)   # 1 µm/s

if __name__ == "__main__":
    unittest.main()
