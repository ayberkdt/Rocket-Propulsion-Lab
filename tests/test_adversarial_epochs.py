"""Attack burn discontinuities with the rounding of large absolute epochs.

A propagator works in absolute time, often 1e5..1e8 s past its epoch.  A
consumer that asks for the left or right limit at a burn boundary forms that
time as ``(t_ignition + bias) + t_boundary``; at 1e7 s that sum is only known
to ~2e-9 s, while a tabulated profile's own boundary tolerance is 1e-12 s.
These tests generate profiles with a jump at every boundary and check that
the requested one-sided limit is returned at every boundary, for every epoch
size, and that times well inside a segment are never moved.

Oracles (independent of the library code):

* Left/right limits: the segment end/start values the test itself put into
  the profile; ``a = F s_F / m`` and ``dm/dt = -mdot s_m`` (Newton's second
  law with the propellant flow; NASA Glenn, Rocket Thrust Equation,
  https://www1.grc.nasa.gov/beginners-guide-to-aeronautics/rocket-thrust-equation/).
* Interior values: linear interpolation written out in the test.
* Rounding bound: a time built from a few additions of magnitude |t| is
  correct to a few units in the last place, below 16 eps |t| (D. Goldberg,
  What Every Computer Scientist Should Know About Floating-Point Arithmetic,
  ACM Computing Surveys 23(1), 1991, bibliographic citation).
* Integrated burns: constant thrust/flow ratio c gives Δv = c ln(m0/mf)
  exactly (NASA Glenn, Ideal Rocket Equation,
  https://www1.grc.nasa.gov/beginners-guide-to-aeronautics/ideal-rocket-equation/),
  and mf = m0 - Σ (q0 + q1) d / 2 per linear-flow segment.
* Left/right-limit contract: Orekit 13.1.5 EventDetector API,
  https://www.orekit.org/static/apidocs/org/orekit/propagation/events/EventDetector.html
"""

import random
import unittest
from math import isclose, log, ulp

from rocket_propulsion.core.errors import DomainError
from rocket_propulsion.integrations import LunarisFiniteBurnDynamics
from rocket_propulsion.propulsion.burns import (
    PropagatorPropulsionBridge,
    PropulsionTimeline,
    ScheduledBurn,
    TabulatedBurnArtifact,
    TabulatedBurnSegment,
    epoch_time_tolerance,
)
from tests.test_lunaris_integration import _FreeSpaceDynamics, integrate_like_lunaris

MASS = 1000.0


def jumpy_profile(generator: random.Random, segments: int) -> tuple[TabulatedBurnArtifact, list]:
    """Random profile with a thrust and flow jump at every interior boundary."""

    built, values, time = [], [], 0.0
    for _ in range(segments):
        duration = generator.uniform(0.2, 40.0)
        f0, f1 = generator.uniform(10.0, 900.0), generator.uniform(10.0, 900.0)
        q0, q1 = generator.uniform(0.01, 0.4), generator.uniform(0.01, 0.4)
        built.append(TabulatedBurnSegment(time, time + duration, f0, f1, f0, f1, q0, q1))
        values.append((time, time + duration, f0, f1, q0, q1))
        time += duration
    artifact = TabulatedBurnArtifact(segments=tuple(built), source_id="jumpy",
                                     source_sha256="e" * 64)
    return artifact, values


def timeline_for(artifact: TabulatedBurnArtifact, ignition: float) -> PropulsionTimeline:
    bridge = PropagatorPropulsionBridge(artifact, require_complete_stream_binding=False)
    return PropulsionTimeline(burns=(ScheduledBurn("b", bridge, ignition),))


def thrust_and_flow(timeline, time_s, side, bias, scale_f=1.0, scale_q=1.0):
    value = timeline.evaluate(
        time_s=time_s, vehicle_mass_kg=MASS, directions={"b": (1.0, 0.0, 0.0)}, side=side,
        parameter_overrides={"b": {"ignition_time_bias_s": bias, "thrust_scale": scale_f,
                                   "mass_flow_scale": scale_q}},
        apply_inventory_limits=False,
    )
    return value.acceleration_m_s2[0] * MASS, -value.mass_derivative_kg_s


class OneSidedLimitAttackTests(unittest.TestCase):
    def test_every_boundary_keeps_both_limits_at_every_epoch_size(self) -> None:
        generator = random.Random(20261007)
        failures = []
        checked = 0
        for case in range(600):
            artifact, values = jumpy_profile(generator, generator.randint(1, 5))
            ignition = 10.0 ** generator.uniform(0.0, 8.5) * generator.choice((1.0, 1.0, -1.0))
            ignition = abs(ignition)
            bias = generator.choice((0.0, 0.37, -0.0123, generator.uniform(-2.0, 2.0)))
            scale_f, scale_q = generator.uniform(0.9, 1.1), generator.uniform(0.9, 1.1)
            timeline = timeline_for(artifact, ignition)
            boundaries = [0.0] + [item[1] for item in values]
            for index, boundary in enumerate(boundaries):
                # Two ways a consumer may form the absolute time.
                for time_s in ((ignition + bias) + boundary, ignition + (bias + boundary)):
                    left_expected = (0.0, 0.0) if index == 0 else (
                        values[index - 1][3] * scale_f, values[index - 1][5] * scale_q)
                    right_expected = (0.0, 0.0) if index == len(values) else (
                        values[index][2] * scale_f, values[index][4] * scale_q)
                    for side, expected in (("left", left_expected), ("right", right_expected)):
                        checked += 1
                        actual = thrust_and_flow(timeline, time_s, side, bias, scale_f, scale_q)
                        if not all(isclose(a, e, rel_tol=1e-9, abs_tol=1e-9)
                                   for a, e in zip(actual, expected, strict=True)):
                            failures.append((case, ignition, bias, boundary, side, actual,
                                             expected))
        self.assertGreater(checked, 5000)
        self.assertEqual(failures[:5], [], f"{len(failures)} lost limits")

    def test_interior_times_are_never_snapped(self) -> None:
        generator = random.Random(7)
        for _ in range(400):
            artifact, values = jumpy_profile(generator, generator.randint(1, 4))
            ignition = 10.0 ** generator.uniform(0.0, 8.0)
            timeline = timeline_for(artifact, ignition)
            start, end, f0, f1, q0, q1 = generator.choice(values)
            fraction = generator.uniform(0.05, 0.95)
            time_s = ignition + start + fraction * (end - start)
            # Independent linear interpolation at the time actually represented.
            local = ((time_s - ignition) - start) / (end - start)
            expected = (f0 + (f1 - f0) * local, q0 + (q1 - q0) * local)
            for side in ("left", "right"):
                actual = thrust_and_flow(timeline, time_s, side, 0.0)
                for a, e in zip(actual, expected, strict=True):
                    self.assertTrue(isclose(a, e, rel_tol=1e-9), (a, e, ignition))

    def test_snapping_moves_a_time_only_onto_the_nearest_boundary_within_tolerance(self) -> None:
        generator = random.Random(11)
        for _ in range(300):
            artifact, values = jumpy_profile(generator, generator.randint(1, 6))
            boundaries = [0.0] + [item[1] for item in values]
            tolerance = generator.choice((None, 1e-10, 1e-7, 1e-4))
            used = max(artifact.boundary_tolerance_s, tolerance or 0.0)
            for _ in range(40):
                anchor = generator.choice(boundaries)
                offset = generator.choice((0.0, 0.3, 0.9, 1.1, 3.0, 1e3)) * used
                time_s = anchor + generator.choice((-1.0, 1.0)) * offset
                snapped = artifact.snap_to_boundary(time_s, tolerance)
                nearest = min(boundaries, key=lambda b, t=time_s: abs(b - t))
                if abs(nearest - time_s) <= used:
                    self.assertEqual(snapped, nearest)
                else:
                    self.assertEqual(snapped, time_s)

    def test_resolution_coarser_than_the_profile_fails_closed(self) -> None:
        # At 1e9 s a time is only known to ~3.6e-6 s; a 2 µs segment can then
        # not be told apart from its neighbours.  Refuse instead of guessing.
        segments = (
            TabulatedBurnSegment(0.0, 1.0, 100.0, 100.0, 100.0, 100.0, 0.1, 0.1),
            TabulatedBurnSegment(1.0, 1.000002, 500.0, 500.0, 500.0, 500.0, 0.5, 0.5),
            TabulatedBurnSegment(1.000002, 2.0, 100.0, 100.0, 100.0, 100.0, 0.1, 0.1),
        )
        artifact = TabulatedBurnArtifact(segments=segments, source_id="short",
                                         source_sha256="e" * 64)
        timeline = timeline_for(artifact, 1e9)
        self.assertGreater(epoch_time_tolerance(1e9), 0.5 * 2e-6)
        with self.assertRaises(DomainError):
            thrust_and_flow(timeline, 1e9 + 1.0, "right", 0.0)
        # The same profile at a small epoch is resolved normally.
        near = timeline_for(artifact, 10.0)
        self.assertAlmostEqual(thrust_and_flow(near, 11.0, "right", 0.0)[0], 500.0)
        self.assertAlmostEqual(thrust_and_flow(near, 11.0, "left", 0.0)[0], 100.0)

    def test_epoch_tolerance_bounds_the_actual_rounding(self) -> None:
        generator = random.Random(3)
        for _ in range(20000):
            ignition = 10.0 ** generator.uniform(0.0, 9.0)
            bias = generator.uniform(-5.0, 5.0)
            boundary = generator.uniform(0.0, 3000.0)
            formed = (ignition + bias) + boundary
            error = abs((formed - ignition) - bias - boundary)
            self.assertLessEqual(error, epoch_time_tolerance(formed, ignition, bias))


class LargeEpochIntegrationTests(unittest.TestCase):
    """Fixed-step RK4 takes the full error of any lost limit; it must see none."""

    def test_biased_burns_at_large_epochs_integrate_to_the_rocket_equation(self) -> None:
        c = 3000.0
        generator = random.Random(1)
        cases = [(ignition, bias) for ignition in (1e3, 123_456.789, 2_625_183.3548202747, 4.4e7)
                 for bias in (0.0, 0.37, -0.81)]
        # Random epochs: whether a given boundary sum rounds up or down is
        # effectively random, so a fixed handful of epochs can miss it.
        cases += [(generator.uniform(1e5, 1e7), generator.choice((0.0, 0.37)))
                  for _ in range(40)]
        for ignition, bias in cases:
            if True:
                segments = (
                    TabulatedBurnSegment(0.0, 1.3, 0.0, 420.0, 0.0, 420.0, 0.0, 420.0 / c),
                    TabulatedBurnSegment(1.3, 62.0, 420.0, 420.0, 420.0, 420.0, 0.14, 0.14),
                    # Jump down at 62 s and an abrupt cutoff at full thrust: both
                    # limits differ at every boundary, so a lost limit shows.
                    TabulatedBurnSegment(62.0, 63.1, 210.0, 210.0, 210.0, 210.0, 0.07, 0.07),
                )
                artifact = TabulatedBurnArtifact(
                    segments=tuple(
                        TabulatedBurnSegment(s.start_time_s, s.end_time_s,
                                             s.start_delivered_thrust_n, s.end_delivered_thrust_n,
                                             s.start_axial_thrust_n, s.end_axial_thrust_n,
                                             s.start_delivered_thrust_n / c,
                                             s.end_delivered_thrust_n / c)
                        for s in segments
                    ),
                    source_id="ramp", source_sha256="f" * 64,
                )
                timeline = timeline_for(artifact, ignition)
                overrides = {"b": {"ignition_time_bias_s": bias}}
                dynamics = LunarisFiniteBurnDynamics(
                    _FreeSpaceDynamics(), timeline, lambda burn_id, t, y: (1.0, 0.0, 0.0),
                    parameter_overrides=overrides,
                )
                final = integrate_like_lunaris(
                    dynamics, timeline, [0, 0, 0, 0, 0, 0, MASS], ignition - 5.0,
                    ignition + 70.0, 0.37, times=dynamics.segmentation_times,
                )
                consumed = sum((s.start_total_mass_flow_kg_s + s.end_total_mass_flow_kg_s)
                               * s.duration_s / 2 for s in artifact.segments)
                # The test integrator itself advances t by rounded additions: each
                # of its ~200 steps can misplace time by half an ulp of t, which
                # moves at most 0.14 kg/s of flow (or 0.43 m/s^2 of acceleration)
                # across a boundary.  That is the attainable bound at this epoch.
                time_slack = 200 * ulp(ignition + 70.0)
                self.assertLessEqual(abs(final[6] - (MASS - consumed)),
                                     1e-12 * MASS + 0.14 * time_slack,
                                     (ignition, bias, final[6], MASS - consumed))
                expected_dv = c * log(MASS / (MASS - consumed))
                self.assertLessEqual(abs(final[3] - expected_dv),
                                     1e-9 * expected_dv + 0.43 * time_slack, (ignition, bias))


if __name__ == "__main__":
    unittest.main()
