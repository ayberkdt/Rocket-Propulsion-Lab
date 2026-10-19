"""Independent boundary ownership and binary epoch clock adversarial checks."""

import unittest
from dataclasses import replace
from decimal import Decimal, localcontext
from fractions import Fraction
from math import inf, ldexp, nextafter, ulp

from rocket_propulsion.core.errors import DomainError, InputError
from rocket_propulsion.integrations.lunaris import (
    lunaris_segmentation_times,
    LunarisFiniteBurnDynamics,
    lunaris_burn_windows,
)
from rocket_propulsion.propulsion.burns import (
    PropagatorPropulsionBridge,
    PropulsionTimeline,
    ScheduledBurn,
    TabulatedBurnArtifact,
    TabulatedBurnSegment,
)


def profile(h=1.0):
    return TabulatedBurnArtifact(
        (
            TabulatedBurnSegment(0.0, h, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, phase="first"),
            TabulatedBurnSegment(h, 2 * h, 2.0, 2.0, 2.0, 2.0, 0.0, 0.0, phase="second"),
        ),
        "clock oracle",
        "8" * 64,
    )


def bridge(artifact):
    return PropagatorPropulsionBridge(artifact, require_complete_stream_binding=False)


def evaluate(model, time, side="right", **kwargs):
    return model.evaluate(
        relative_time_s=time, vehicle_mass_kg=1.0, direction=(1.0, 0.0, 0.0), side=side, **kwargs
    )


def schedule(artifact, ignition=0.0, bias=0.0):
    model = bridge(artifact)
    if abs(bias) > 1.0:
        model = replace(
            model,
            parameters=tuple(
                replace(p, minimum=-2e16, maximum=2e16) if p.name == "ignition_time_bias_s" else p
                for p in model.parameters
            ),
        )
    return PropulsionTimeline((ScheduledBurn("burn", model, ignition),))


def timeline_value(artifact, time, ignition=0.0, bias=0.0, side="right"):
    timeline = schedule(artifact, ignition, bias)
    return timeline.evaluate(
        time_s=time,
        vehicle_mass_kg=1.0,
        directions={"burn": (1.0, 0.0, 0.0)},
        parameter_overrides={"burn": {"ignition_time_bias_s": bias}},
        side=side,
    )


class TabulatedTimeAttacks(unittest.TestCase):
    def test_real_near_stop_reserve_cut_is_retained(self):
        artifact = TabulatedBurnArtifact(
            (TabulatedBurnSegment(0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 2.0, 2.0),),
            "near stop margin",
            "8" * 64,
        )
        model = replace(bridge(artifact), protected_dry_mass_kg=2.0)
        timeline = PropulsionTimeline((ScheduledBurn("burn", model, 0.0),))
        window = lunaris_burn_windows(timeline, initial_mass_kg=4.0)[0]
        self.assertIsNotNone(window.depletion_time_s)
        self.assertLess(window.depletion_time_s, 1.0)
        self.assertLess(1.0 - window.depletion_time_s, 1e-12)
        self.assertEqual(window.depletion_time_s, float((Fraction(2) - Fraction(1e-12)) / 2))

    def test_epoch_saltation_has_independently_known_force_jump(self):
        model = bridge(profile())
        for bias in (-1e-5, 1e-5):
            jump = model.time_event_jump_at_epoch(
                time_s=1.0 + bias,
                ignition_time_s=0.0,
                vehicle_mass_kg=1.0,
                direction=(1.0, 0.0, 0.0),
                parameter_overrides={"ignition_time_bias_s": bias},
            )
            self.assertEqual(jump.acceleration_m_s2, (-1.0, 0.0, 0.0))
            self.assertEqual(jump.mass_derivative_kg_s, 0.0)

    def test_once_rounded_epoch_dates_match_paired_decimal_reference(self):
        cases = (
            (1e7, 0.37, 0.25),
            (ldexp(1.0, 30), 0.0, ulp(ldexp(1.0, 30))),
            (1e16, -1e16, 5e-21),
            (-1e16, 1e16, 5e-21),
        )
        for ignition, bias, h in cases:
            reference = []
            for precision in (1200, 1600):
                with localcontext() as context:
                    context.prec = precision
                    reference.append(
                        tuple(
                            float(
                                Decimal.from_float(ignition)
                                + Decimal.from_float(bias)
                                + Decimal.from_float(t)
                            )
                            for t in (0.0, h, 2 * h)
                        )
                    )
            self.assertEqual(reference[0], reference[1])
            timeline = schedule(profile(h), ignition, bias)
            actual = lunaris_segmentation_times(
                timeline, parameter_overrides={"burn": {"ignition_time_bias_s": bias}}
            )
            self.assertEqual(actual, reference[0])
            for date, side, force in (
                (actual[0], "right", 1.0),
                (actual[1], "left", 1.0),
                (actual[1], "right", 2.0),
                (actual[2], "left", 2.0),
            ):
                self.assertEqual(
                    timeline_value(profile(h), date, ignition, bias, side).acceleration_m_s2[0],
                    force,
                )

    def test_time_policy_is_required_and_bound_at_each_document_layer(self):
        artifact = profile()
        model = bridge(artifact)
        timeline = schedule(artifact)
        for value in (artifact, model, timeline):
            payload = value.to_dict()
            self.assertEqual(
                payload["time_numerics"], "exact_relative_clock_and_once_rounded_epoch_events_v1"
            )
            for operation in ("missing", "changed"):
                changed = dict(payload)
                if operation == "missing":
                    changed.pop("time_numerics")
                else:
                    changed["time_numerics"] = "unknown_policy"
                with (
                    self.subTest(kind=type(value).__name__, operation=operation),
                    self.assertRaises(InputError),
                ):
                    type(value).from_dict(changed)

    def test_small_real_overlaps_are_not_excused_by_absolute_tolerance(self):
        model = bridge(profile(5e-21))
        with self.assertRaises(DomainError):
            PropulsionTimeline((ScheduledBurn("a", model, 0.0), ScheduledBurn("b", model, 7.5e-21)))
        touching = PropulsionTimeline(
            (ScheduledBurn("a", model, 0.0), ScheduledBurn("b", model, 1e-20))
        )
        self.assertEqual(len(touching.burns), 2)

    def test_native_rhs_keeps_tiny_interiors_and_exact_one_sided_restarts(self):
        class ZeroDynamics:
            def build_rhs(self):
                return lambda t, y: [0.0] * 7

        for h in (5e-21, 1.0):
            dynamics = LunarisFiniteBurnDynamics(
                ZeroDynamics(), schedule(profile(h)), lambda *args: (1.0, 0.0, 0.0)
            )
            y = [0.0] * 6 + [1.0]
            rhs = dynamics.build_rhs()
            self.assertEqual(rhs(0.0, y)[3], 1.0)
            self.assertEqual(rhs(h, y)[3], 1.0)
            self.assertEqual(rhs(1.5 * h, y)[3], 2.0)
            self.assertEqual(rhs(2 * h, y)[3], 2.0)
            self.assertEqual(dynamics.build_rhs()(h, y)[3], 2.0)
            self.assertEqual(dynamics.build_rhs()(2 * h, y)[3], 0.0)

    def test_timeline_time_events_share_the_epoch_evaluation_clock(self):
        ignition = ldexp(1.0, 30)
        h = ulp(ignition)
        timeline = schedule(profile(h), ignition)
        dates = tuple(event.time_s for event in timeline.events() if event.time_s is not None)
        self.assertEqual(dates, (ignition, ignition + h, ignition + 2 * h))
        with self.assertRaises(DomainError):
            schedule(profile(h / 8), ignition).events()

    def test_shifted_fraction_is_not_rounded_onto_a_discontinuity(self):
        model = bridge(profile())
        for bias, force in ((-ldexp(1.0, -54), 2.0), (ldexp(1.0, -55), 1.0)):
            exact = Fraction(1.0) - Fraction(bias)
            self.assertNotEqual(exact, Fraction(1.0))
            self.assertEqual(float(exact), 1.0)
            for side in ("left", "right"):
                value = evaluate(
                    model, 1.0, side, parameter_overrides={"ignition_time_bias_s": bias}
                )
                self.assertEqual(value.acceleration_m_s2[0], force)
                self.assertEqual(value.partials.nonsmooth_parameters, ())

    def test_relative_bias_and_explicit_snap_do_not_change_unrelated_interiors(self):
        model = bridge(profile())
        for time in (0.1, 0.75, 1.25, 1.75):
            for side in ("left", "right"):
                value = evaluate(
                    model, time + 0.125, side, parameter_overrides={"ignition_time_bias_s": 0.125}
                )
                self.assertEqual(value.effective_profile_time_s, time)
                self.assertEqual(value.acceleration_m_s2[0], 1.0 if time < 1.0 else 2.0)
        self.assertEqual(
            evaluate(model, 1.0 + 1e-6, "left", boundary_tolerance_s=1e-5).acceleration_m_s2[0], 1.0
        )

    def test_epoch_aliases_do_not_override_canonical_dates_after_cancellation(self):
        # The right-associated additions all lose the tiny profile. Canonical
        # dates still distinguish its events; ambiguous aliases are discarded.
        h = 5e-21
        artifact = profile(h)
        self.assertEqual(1e16 + (-1e16 + h), 0.0)
        self.assertEqual(
            timeline_value(artifact, 0.0, 1e16, -1e16, "right").acceleration_m_s2[0], 1.0
        )
        self.assertEqual(
            timeline_value(artifact, h, 1e16, -1e16, "right").acceleration_m_s2[0], 2.0
        )

    def test_public_relative_entry_cannot_override_resolved_internal_time(self):
        with self.assertRaises(TypeError):
            evaluate(bridge(profile()), 0.5, _resolved_profile_time=Fraction(2))
        with self.assertRaises(DomainError):
            bridge(profile()).evaluate_at_epoch(
                time_s=0.5,
                ignition_time_s=0.0,
                vehicle_mass_kg=1.0,
                direction=(1.0, 0.0, 0.0),
                _resolved_profile_time=Fraction(2),
            )

    def test_tiny_interior_has_one_owner_and_remains_active(self):
        artifact = profile(5e-21)
        time = 7.5e-21
        self.assertEqual(artifact.snap_to_boundary(time), time)
        self.assertFalse(artifact.is_boundary_time(time))
        for side in ("left", "right"):
            self.assertEqual(artifact.segment_index_at(time, side=side), 1)
            self.assertEqual(artifact.state_at(time, side=side)["phase"], "second")
            value = evaluate(bridge(artifact), time, side)
            self.assertTrue(value.active)
            self.assertEqual(value.acceleration_m_s2, (2.0, 0.0, 0.0))

    def test_binary_neighbors_use_exact_half_open_ownership_across_scales(self):
        for exponent in (-1074, -1000, -500, -50, 0, 40, 500, 1000):
            h = ldexp(1.0, exponent)
            artifact = profile(h)
            times = {0.0, h, 2 * h, nextafter(h, 0.0), nextafter(h, inf), nextafter(2 * h, 0.0)}
            for time in sorted(times):
                for side in ("left", "right"):
                    with self.subTest(exponent=exponent, time=time, side=side):
                        expected = 0 if time < h or (time == h and side == "left") else 1
                        self.assertEqual(artifact.segment_index_at(time, side=side), expected)
                        self.assertEqual(artifact.snap_to_boundary(time), time)
                        self.assertEqual(artifact.is_boundary_time(time), time in (0.0, h, 2 * h))

    def test_only_exact_endpoints_control_bridge_activity_and_nonsmoothness(self):
        artifact = profile()
        model = bridge(artifact)
        for time in (
            nextafter(0.0, 1.0),
            nextafter(1.0, 0.0),
            nextafter(1.0, 2.0),
            nextafter(2.0, 0.0),
        ):
            for side in ("left", "right"):
                value = evaluate(model, time, side)
                self.assertTrue(value.active, (time, side))
                self.assertEqual(value.partials.nonsmooth_parameters, ())
        self.assertFalse(evaluate(model, 0.0, "left").active)
        self.assertFalse(evaluate(model, 2.0, "right").active)
        self.assertTrue(evaluate(model, 0.0, "right").active)
        self.assertTrue(evaluate(model, 2.0, "left").active)
        for side, force in (("left", 1.0), ("right", 2.0)):
            self.assertEqual(evaluate(model, 1.0, side).acceleration_m_s2[0], force)

    def test_explicit_snap_uses_exact_distance_and_no_hidden_floor(self):
        h = 5e-21
        artifact = profile(h)
        t = h + 0.1 * h
        tolerance = 0.2 * h
        self.assertEqual(artifact.snap_to_boundary(t, 0.0), t)
        self.assertEqual(artifact.snap_to_boundary(t, tolerance), h)
        self.assertEqual(artifact.snap_to_boundary(-0.1 * h, tolerance), 0.0)
        self.assertEqual(artifact.snap_to_boundary(2 * h + 0.1 * h, tolerance), 2 * h)
        self.assertEqual(artifact.snap_to_boundary(h + 0.3 * h, tolerance), h + 0.3 * h)
        for bad in (True, "0", -1.0, inf, float("nan"), 0.5 * h):
            with self.subTest(bad=bad), self.assertRaises(DomainError):
                artifact.snap_to_boundary(t, bad)

    def test_subnormal_zero_tolerance_does_not_compare_against_rounded_half(self):
        h = nextafter(0.0, 1.0)
        artifact = profile(h)
        for time in (0.0, h, 2 * h):
            self.assertEqual(artifact.snap_to_boundary(time, 0.0), time)
        with self.assertRaises(DomainError):
            artifact.snap_to_boundary(h, h)

    def test_declared_zero_origin_and_contiguity_are_exact(self):
        for origin in (5e-13, nextafter(0.0, 1.0)):
            with self.subTest(origin=origin), self.assertRaises(DomainError):
                TabulatedBurnArtifact(
                    (TabulatedBurnSegment(origin, 1.0, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0),),
                    "shifted",
                    "8" * 64,
                )
        for start in (nextafter(1.0, 0.0), nextafter(1.0, inf)):
            with self.subTest(start=start), self.assertRaises(DomainError):
                TabulatedBurnArtifact(
                    (
                        profile().segments[0],
                        TabulatedBurnSegment(start, 2.0, 2.0, 2.0, 2.0, 2.0, 0.0, 0.0),
                    ),
                    "gap or overlap",
                    "8" * 64,
                )
        self.assertEqual(profile(nextafter(0.0, 1.0)).duration_s, 2 * nextafter(0.0, 1.0))

    def test_lookup_is_logarithmic_without_overlap_walk(self):
        # Structural rejection keeps binary search meaningful; no backward walk
        # over thousands of overlapping intervals can be hidden in a lookup.
        rows = tuple(
            TabulatedBurnSegment(i * 1e-20, (i + 2) * 1e-20, 1.0, 1.0, 1.0, 1.0, 0.0, 0.0)
            for i in range(50)
        )
        with self.assertRaises(DomainError):
            TabulatedBurnArtifact(rows, "overlap", "8" * 64)

    def test_zero_epoch_preserves_tiny_and_subnormal_profile_clocks(self):
        for h in (nextafter(0.0, 1.0), 5e-21, 1.0):
            artifact = profile(h)
            for time, side, force in (
                (0.0, "right", 1.0),
                (h, "left", 1.0),
                (h, "right", 2.0),
                (2 * h, "left", 2.0),
            ):
                with self.subTest(h=h, time=time, side=side):
                    self.assertEqual(
                        timeline_value(artifact, time, side=side).acceleration_m_s2[0], force
                    )

    def test_absolute_clock_accepts_distinct_one_ulp_events(self):
        ignition = ldexp(1.0, 30)
        h = ulp(ignition)
        artifact = profile(h)
        for side, force in (("left", 1.0), ("right", 2.0)):
            self.assertEqual(
                timeline_value(artifact, ignition + h, ignition, side=side).acceleration_m_s2[0],
                force,
            )

    def test_actual_absolute_event_collisions_fail_closed(self):
        ignition = ldexp(1.0, 30)
        artifact = profile(ulp(ignition) / 8)
        with self.assertRaises(DomainError):
            timeline_value(artifact, ignition, ignition)
        timeline = PropulsionTimeline((ScheduledBurn("burn", bridge(artifact), ignition),))
        with self.assertRaises(DomainError):
            lunaris_segmentation_times(timeline)

    def test_epoch_neighbors_are_interiors_not_tolerance_wide_limits(self):
        ignition = 1e7
        time = nextafter(ignition + 1.0, inf)
        for side in ("left", "right"):
            self.assertEqual(
                timeline_value(profile(), time, ignition, side=side).acceleration_m_s2[0], 2.0
            )

    def test_exact_epoch_cancellation_preserves_recoverable_profile_time(self):
        ignition = 1e16
        bias = -1e16
        artifact = profile(5e-21)
        # Fraction arithmetic is an independent ledger of the supplied binary
        # clock, not subtraction of two prematurely rounded intermediate times.
        time = 7.5e-21
        self.assertEqual(Fraction(time) - Fraction(ignition) - Fraction(bias), Fraction(time))
        self.assertEqual(timeline_value(artifact, time, ignition, bias).acceleration_m_s2[0], 2.0)

    def test_native_segmentation_keeps_every_distinct_tiny_event(self):
        for ignition, h in (
            (0.0, 5e-21),
            (0.0, nextafter(0.0, 1.0)),
            (ldexp(1.0, 30), ulp(ldexp(1.0, 30))),
        ):
            timeline = PropulsionTimeline((ScheduledBurn("burn", bridge(profile(h)), ignition),))
            self.assertEqual(
                lunaris_segmentation_times(timeline), (ignition, ignition + h, ignition + 2 * h)
            )

    def test_ordinary_limits_and_serialized_roundtrip_control(self):
        artifact = profile(2.0)
        opened = TabulatedBurnArtifact.from_json(artifact.to_json())
        self.assertEqual(opened.artifact_hash, artifact.artifact_hash)
        self.assertEqual(opened.state_at(2.0, side="left")["phase"], "first")
        self.assertEqual(opened.state_at(2.0, side="right")["phase"], "second")
        self.assertEqual(evaluate(bridge(opened), 3.0).acceleration_m_s2[0], 2.0)
