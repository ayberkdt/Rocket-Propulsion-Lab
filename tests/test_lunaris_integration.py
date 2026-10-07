"""Lunaris adapter: exact-mass impulsive equivalents and finite-burn dynamics.

The first part needs no Lunaris installation.  ``LunarisEndToEndTests`` runs
the adapter against the real Lunaris propagator when ``lunaris`` is
importable and is skipped otherwise.
"""

import importlib.util
import unittest
from itertools import pairwise
from math import exp, isclose, log, sqrt

from rocket_propulsion.core.errors import DomainError
from rocket_propulsion.integrations import (
    LunarisFiniteBurnDynamics,
    lunaris_burn_windows,
    lunaris_impulsive_equivalents,
    lunaris_segmentation_times,
)
from rocket_propulsion.propulsion.burns import (
    MassJettison,
    PropagatorPropulsionBridge,
    PropulsionTimeline,
    ScheduledBurn,
    StreamStateBinding,
    TabulatedBurnArtifact,
    TabulatedBurnSegment,
    free_space_burn_metrics,
)
from rocket_propulsion.propulsion.performance import STANDARD_GRAVITY_M_S2

LUNARIS_AVAILABLE = importlib.util.find_spec("lunaris") is not None


def mass_only_bridge(thrust: float, flow: float, duration: float, ramp: float = 0.0):
    segments = []
    if ramp > 0.0:
        segments.append(TabulatedBurnSegment(0.0, ramp, 0.0, thrust, 0.0, thrust, 0.0, flow,
                                             phase="ramp"))
    segments.append(TabulatedBurnSegment(ramp, ramp + duration, thrust, thrust, thrust, thrust,
                                         flow, flow, phase="steady"))
    artifact = TabulatedBurnArtifact(
        segments=tuple(segments), source_id=f"lunaris-{thrust}-{duration}",
        source_sha256="7" * 64, reference_ids=("NASA-20040000363",),
    )
    return PropagatorPropulsionBridge(artifact, require_complete_stream_binding=False)


def two_burns() -> PropulsionTimeline:
    return PropulsionTimeline(
        burns=(
            ScheduledBurn("lower", mass_only_bridge(400.0, 0.14, 60.0, ramp=2.0), 100.0),
            ScheduledBurn("raise", mass_only_bridge(400.0, 0.14, 40.0), 3600.0),
        )
    )


def lunaris_mass_law(mass: float, dv: tuple[float, float, float], isp: float) -> float:
    """Lunaris' apply_impulsive_maneuver mass update, reimplemented for the test."""

    return mass * exp(-sqrt(sum(v * v for v in dv)) / (STANDARD_GRAVITY_M_S2 * isp))


class ImpulsiveEquivalentTests(unittest.TestCase):
    def test_mass_and_delta_v_are_exact_and_sequenced(self) -> None:
        timeline = two_burns()
        impulses = lunaris_impulsive_equivalents(
            timeline,
            initial_mass_kg=1000.0,
            directions={"lower": ((0.0, -1.0, 0.0), "ric"), "raise": ((0.0, 2.0, 0.0), "ric")},
        )
        self.assertEqual([item.burn_id for item in impulses], ["lower", "raise"])
        mass = 1000.0
        for item, burn in zip(impulses, timeline.burns, strict=True):
            metrics = free_space_burn_metrics(burn.bridge.artifact, mass)
            self.assertAlmostEqual(item.delta_v_m_s, metrics.delta_v_m_s, places=12)
            self.assertAlmostEqual(sqrt(sum(v * v for v in item.dv_mps)), item.delta_v_m_s)
            self.assertTrue(isclose(lunaris_mass_law(mass, item.dv_mps, item.isp_s),
                                    item.final_mass_kg, rel_tol=1e-14))
            self.assertAlmostEqual(item.final_mass_kg, mass - burn.bridge.artifact.consumed_propellant_kg)
            self.assertGreater(item.t_burn_s, burn.ignition_time_s)
            self.assertLess(item.t_burn_s, burn.cutoff_time_s)
            mass = item.final_mass_kg
        self.assertEqual(impulses[0].dv_mps[1] < 0.0, True)

    def test_refusals(self) -> None:
        directions = {"lower": ((1.0, 0.0, 0.0), "inertial"), "raise": ((1.0, 0.0, 0.0), "inertial")}
        with self.assertRaises(DomainError):
            lunaris_impulsive_equivalents(two_burns(), initial_mass_kg=1000.0,
                                          directions={"lower": directions["lower"]})
        with self.assertRaises(DomainError):
            lunaris_impulsive_equivalents(
                two_burns(), initial_mass_kg=1000.0,
                directions={"lower": ((1.0, 0.0, 0.0), "vnb"), "raise": directions["raise"]},
            )
        jettison = PropulsionTimeline(
            burns=two_burns().burns, jettisons=(MassJettison("stage", 2000.0, 10.0),)
        )
        with self.assertRaises(DomainError):
            lunaris_impulsive_equivalents(jettison, initial_mass_kg=1000.0, directions=directions)
        with self.assertRaises(DomainError):
            lunaris_impulsive_equivalents(two_burns(), initial_mass_kg=-1.0, directions=directions)


class _FreeSpaceDynamics:
    """Lunaris-style dynamics without gravity, to isolate the propulsion term."""

    grav = None
    marker = "delegated"

    def build_rhs(self):
        def rhs(_t, y):
            return [y[3], y[4], y[5], 0.0, 0.0, 0.0, 0.0]

        return rhs


def integrate_like_lunaris(dynamics, timeline, y0, t0, tf, step, times=None):
    """Fresh RHS per segment, segments split at every discontinuity (as Lunaris does)."""

    times = lunaris_segmentation_times(timeline) if times is None else times
    cuts = [t for t in times if t0 < t < tf]
    bounds = [t0, *cuts, tf]
    y = list(y0)
    for start, end in pairwise(bounds):
        rhs = dynamics._build_solver_rhs()
        t = start
        while t < end - 1e-12:
            h = min(step, end - t)
            k1 = rhs(t, y)
            k2 = rhs(t + h / 2, [a + h / 2 * b for a, b in zip(y, k1, strict=True)])
            k3 = rhs(t + h / 2, [a + h / 2 * b for a, b in zip(y, k2, strict=True)])
            k4 = rhs(t + h, [a + h * b for a, b in zip(y, k3, strict=True)])
            y = [a + h / 6 * (b + 2 * c + 2 * d + e)
                 for a, b, c, d, e in zip(y, k1, k2, k3, k4, strict=True)]
            t += h
    return y


class FiniteBurnDynamicsTests(unittest.TestCase):
    def test_free_space_burns_reproduce_the_rocket_equation(self) -> None:
        timeline = two_burns()
        dynamics = LunarisFiniteBurnDynamics(
            _FreeSpaceDynamics(), timeline, lambda burn_id, t, y: (1.0, 0.0, 0.0)
        )
        self.assertEqual(dynamics.marker, "delegated")
        final = integrate_like_lunaris(
            dynamics, timeline, [0, 0, 0, 0, 0, 0, 1000.0], 0.0, 4000.0, 0.5
        )
        m1 = 1000.0 - timeline.burns[0].bridge.artifact.consumed_propellant_kg
        m2 = m1 - timeline.burns[1].bridge.artifact.consumed_propellant_kg
        self.assertAlmostEqual(final[6], m2, places=9)
        expected = (free_space_burn_metrics(timeline.burns[0].bridge.artifact, 1000.0).delta_v_m_s
                    + 400.0 / 0.14 * log(m1 / m2))
        self.assertTrue(isclose(final[3], expected, rel_tol=1e-10))

    def test_refuses_tank_states_conditions_and_jettisons(self) -> None:
        flows = (("ox", 1.0),)
        bound = PropagatorPropulsionBridge(
            TabulatedBurnArtifact(
                segments=(TabulatedBurnSegment(0.0, 1.0, 10, 10, 10, 10, 1, 1, flows, flows),),
                source_id="b", source_sha256="8" * 64,
            ),
            stream_bindings=(StreamStateBinding("ox", "m_ox"),),
        )
        with self.assertRaises(DomainError):
            LunarisFiniteBurnDynamics(
                _FreeSpaceDynamics(), PropulsionTimeline(burns=(ScheduledBurn("a", bound, 0.0),))
            )
        with self.assertRaises(DomainError):
            LunarisFiniteBurnDynamics(
                _FreeSpaceDynamics(),
                PropulsionTimeline(burns=two_burns().burns,
                                   jettisons=(MassJettison("s", 2000.0, 1.0),)),
            )
        with self.assertRaises(DomainError):
            LunarisFiniteBurnDynamics(object(), two_burns())
        dynamics = LunarisFiniteBurnDynamics(_FreeSpaceDynamics(), two_burns())
        with self.assertRaises(DomainError):
            dynamics.build_rhs()(0.0, [1, 0, 0, 0, 1, 0])  # 6-state has no mass


def protected_bridge(thrust, flow, duration, ramp, dry_mass):
    bridge = mass_only_bridge(thrust, flow, duration, ramp=ramp)
    return PropagatorPropulsionBridge(
        bridge.artifact, require_complete_stream_binding=False, protected_dry_mass_kg=dry_mass
    )


class BurnWindowTests(unittest.TestCase):
    """Bias-shifted events and exact dry-mass cuts become segment boundaries."""

    def test_ignition_bias_shifts_every_segmentation_time(self) -> None:
        timeline = two_burns()
        overrides = {"lower": {"ignition_time_bias_s": 3.25}, "raise": {"ignition_time_bias_s": -7.5}}
        nominal = lunaris_segmentation_times(timeline)
        shifted = lunaris_segmentation_times(timeline, parameter_overrides=overrides)
        self.assertEqual(shifted, (103.25, 105.25, 165.25, 3592.5, 3632.5))
        self.assertEqual(nominal, (100.0, 102.0, 162.0, 3600.0, 3640.0))
        dynamics = LunarisFiniteBurnDynamics(
            _FreeSpaceDynamics(), timeline, lambda burn_id, t, y: (1.0, 0.0, 0.0),
            parameter_overrides=overrides,
        )
        self.assertEqual(dynamics.segmentation_times, shifted)
        exact = integrate_like_lunaris(dynamics, timeline, [0, 0, 0, 0, 0, 0, 1000.0], 0.0,
                                       4000.0, 0.5, times=dynamics.segmentation_times)
        reference = integrate_like_lunaris(
            LunarisFiniteBurnDynamics(_FreeSpaceDynamics(), timeline,
                                      lambda burn_id, t, y: (1.0, 0.0, 0.0)),
            timeline, [0, 0, 0, 0, 0, 0, 1000.0], 0.0, 4000.0, 0.5,
        )
        # A shift in time does not change free-space delta-v or final mass.
        self.assertAlmostEqual(exact[6], reference[6], places=9)
        self.assertTrue(isclose(exact[3], reference[3], rel_tol=1e-10))
        # The pre-fix behaviour (nominal cuts with a biased burn) leaves the
        # ramp kink and both edges inside RK4 steps and loses accuracy.
        stale = integrate_like_lunaris(dynamics, timeline, [0, 0, 0, 0, 0, 0, 1000.0], 0.0,
                                       4000.0, 0.5, times=nominal)
        self.assertGreater(abs(stale[3] - reference[3]), 1e3 * abs(exact[3] - reference[3]))

    def test_bias_that_creates_an_overlap_is_refused(self) -> None:
        with self.assertRaises(DomainError):
            lunaris_segmentation_times(
                two_burns(), parameter_overrides={"raise": {"ignition_time_bias_s": -3490.0}}
            )
        with self.assertRaises(DomainError):
            lunaris_segmentation_times(two_burns(), parameter_overrides={"ghost": {}})

    def test_dry_mass_cut_time_is_exact_inside_the_ramp_and_steady_phase(self) -> None:
        # Ramp 0..4 s to q = 0.5 kg/s, then steady.  Ramp consumes q * 4 / 2 = 1 kg.
        for reserve, expected in ((0.5, 4.0 * sqrt(0.5 / 1.0)), (3.0, 4.0 + 2.0 / 0.5)):
            timeline = PropulsionTimeline(
                burns=(ScheduledBurn("burn", protected_bridge(1500.0, 0.5, 100.0, 4.0,
                                                               100.0 - reserve), 10.0),)
            )
            (window,) = lunaris_burn_windows(timeline, initial_mass_kg=100.0)
            self.assertAlmostEqual(window.depletion_time_s, 10.0 + expected, places=12)
            self.assertAlmostEqual(window.cutoff_time_s, 10.0 + expected, places=12)
            self.assertIn(window.depletion_time_s, lunaris_segmentation_times(
                timeline, initial_mass_kg=100.0))

    def test_finite_burn_stops_exactly_at_protected_dry_mass(self) -> None:
        dry = 940.0
        timeline = PropulsionTimeline(
            burns=(ScheduledBurn("burn", protected_bridge(1500.0, 0.5, 400.0, 4.0, dry), 10.0),)
        )
        for limits in (None, False):
            dynamics = LunarisFiniteBurnDynamics(
                _FreeSpaceDynamics(), timeline, lambda burn_id, t, y: (1.0, 0.0, 0.0),
                initial_mass_kg=1000.0, apply_inventory_limits=limits,
            )
            final = integrate_like_lunaris(dynamics, timeline, [0, 0, 0, 0, 0, 0, 1000.0], 0.0,
                                           600.0, 0.5, times=dynamics.segmentation_times)
            self.assertAlmostEqual(final[6], dry, places=9)
            # Thrust/flow ratio is constant (ramp included): c = 3000 m/s.
            self.assertTrue(isclose(final[3], 3000.0 * log(1000.0 / dry), rel_tol=1e-10))
            self.assertAlmostEqual(dynamics.predicted_mass_kg(600.0), dry, places=12)
        # Forcing the internal clamp on lets RK stages next to the cut see the
        # dry mass and stop the flow early: the documented bias, caught here.
        clamped = LunarisFiniteBurnDynamics(
            _FreeSpaceDynamics(), timeline, lambda burn_id, t, y: (1.0, 0.0, 0.0),
            initial_mass_kg=1000.0, apply_inventory_limits=True,
        )
        with self.assertRaises(DomainError):
            integrate_like_lunaris(clamped, timeline, [0, 0, 0, 0, 0, 0, 1000.0], 0.0, 600.0,
                                   0.5, times=clamped.segmentation_times)

    def test_protected_dry_mass_needs_consistent_initial_mass(self) -> None:
        timeline = PropulsionTimeline(
            burns=(ScheduledBurn("burn", protected_bridge(1500.0, 0.5, 400.0, 4.0, 940.0), 10.0),)
        )
        with self.assertRaises(DomainError):
            LunarisFiniteBurnDynamics(_FreeSpaceDynamics(), timeline)
        dynamics = LunarisFiniteBurnDynamics(_FreeSpaceDynamics(), timeline,
                                             initial_mass_kg=1000.0)
        with self.assertRaises(DomainError):
            dynamics.build_rhs()(0.0, [1, 0, 0, 0, 1, 0, 990.0])
        # Impulsive equivalents refuse a burn that its dry mass would truncate.
        with self.assertRaises(DomainError):
            lunaris_impulsive_equivalents(
                timeline, initial_mass_kg=1000.0, directions={"burn": ((1.0, 0, 0), "inertial")}
            )


@unittest.skipUnless(LUNARIS_AVAILABLE, "lunaris is not installed")
class LunarisEndToEndTests(unittest.TestCase):
    def setUp(self) -> None:
        from lunaris.common.constants import MU_MOON, R_MOON
        from lunaris.common.type_defs import EventConfig, PropagatorConfig, TimeConfig

        self.mu = MU_MOON
        r = R_MOON + 100_000.0
        self.y0 = [r, 0.0, 0.0, 0.0, sqrt(MU_MOON / r), 0.0, 1000.0]
        self.cfg = PropagatorConfig(
            method="DOP853", rtol=1e-12, atol=1e-13, verbose=False, use_nyquist_max_step=False,
            events=EventConfig(detect_impact=False, enable_peri_apo_events=False),
        )
        self.time_cfg = TimeConfig(duration_s=4000.0, output_dt_s=50.0)

        mu = MU_MOON

        class PointMassMoon:
            grav = None
            ephem = None
            flags = None

            def build_rhs(self):
                def rhs(_t, y):
                    import numpy as np

                    state = np.asarray(y, dtype=np.float64)
                    out = np.empty_like(state)
                    out[:3] = state[3:6]
                    rn = float(np.linalg.norm(state[:3]))
                    out[3:6] = -mu * state[:3] / rn**3
                    out[6] = 0.0
                    return out

                return rhs

        self.base = PointMassMoon()

    def semi_major_axis(self, y) -> float:
        r = sqrt(sum(value * value for value in y[:3]))
        v2 = sum(value * value for value in y[3:6])
        return 1.0 / (2.0 / r - v2 / self.mu)

    def test_finite_and_impulsive_paths_agree_on_mass_and_orbit(self) -> None:
        from lunaris.core.propagation import propagate

        from rocket_propulsion.integrations import (
            lunaris_segmentation_plan,
            to_lunaris_maneuver_plan,
            velocity_direction,
        )

        timeline = two_burns()
        finite = propagate(
            LunarisFiniteBurnDynamics(self.base, timeline, velocity_direction),
            self.y0, self.cfg, time_cfg=self.time_cfg,
            maneuver_plan=lunaris_segmentation_plan(timeline),
        )
        impulses = lunaris_impulsive_equivalents(
            timeline, initial_mass_kg=1000.0,
            directions={"lower": ((0.0, 1.0, 0.0), "ric"), "raise": ((0.0, 1.0, 0.0), "ric")},
        )
        impulsive = propagate(self.base, self.y0, self.cfg, time_cfg=self.time_cfg,
                              maneuver_plan=to_lunaris_maneuver_plan(impulses))
        expected_mass = impulses[-1].final_mass_kg
        self.assertAlmostEqual(float(finite.y[-1][6]), expected_mass, places=8)
        self.assertAlmostEqual(float(impulsive.y[-1][6]), expected_mass, places=10)
        a0 = self.semi_major_axis(self.y0)
        a_finite = self.semi_major_axis(finite.y[-1])
        a_impulsive = self.semi_major_axis(impulsive.y[-1])
        self.assertGreater(a_finite, a0)
        # A 62 s and a 40 s burn on a ~2 h orbit: finite-duration loss is tiny.
        self.assertLess(abs(a_finite - a_impulsive) / (a_finite - a0), 1e-3)
        self.assertEqual(len(impulsive.diagnostics["maneuvers_applied"]), 2)

    def test_dry_mass_cut_and_bias_match_an_explicitly_truncated_burn(self) -> None:
        from lunaris.core.propagation import propagate

        from rocket_propulsion.integrations import lunaris_segmentation_plan, velocity_direction

        # 600 s profile cut by its protected dry mass after 250 s (0.14 kg/s),
        # fired 12.5 s late through the timing driver ...
        protected = PropulsionTimeline(
            burns=(ScheduledBurn("burn", protected_bridge(400.0, 0.14, 600.0, 0.0, 965.0),
                                 300.0),)
        )
        overrides = {"burn": {"ignition_time_bias_s": 12.5}}
        dynamics = LunarisFiniteBurnDynamics(
            self.base, protected, velocity_direction, parameter_overrides=overrides,
            initial_mass_kg=1000.0,
        )
        self.assertAlmostEqual(dynamics.burn_windows[0].depletion_time_s, 562.5, places=9)
        cut = propagate(dynamics, self.y0, self.cfg, time_cfg=self.time_cfg,
                        maneuver_plan=dynamics.segmentation_plan())
        # ... must equal a 250 s burn at 312.5 s with no dry-mass logic at all.
        explicit = PropulsionTimeline(
            burns=(ScheduledBurn("burn", mass_only_bridge(400.0, 0.14, 250.0), 312.5),)
        )
        reference = propagate(
            LunarisFiniteBurnDynamics(self.base, explicit, velocity_direction),
            self.y0, self.cfg, time_cfg=self.time_cfg,
            maneuver_plan=lunaris_segmentation_plan(explicit),
        )
        self.assertAlmostEqual(float(cut.y[-1][6]), 965.0, places=8)
        self.assertAlmostEqual(float(reference.y[-1][6]), 965.0, places=8)
        for index in range(6):
            scale = 1.0 if index >= 3 else 1e3
            self.assertAlmostEqual(float(cut.y[-1][index]) / scale,
                                   float(reference.y[-1][index]) / scale, places=6)
        # A different start mass is refused instead of silently mis-cutting.
        with self.assertRaises(DomainError):
            propagate(dynamics, [*self.y0[:6], 999.0], self.cfg, time_cfg=self.time_cfg,
                      maneuver_plan=dynamics.segmentation_plan())


if __name__ == "__main__":
    unittest.main()
