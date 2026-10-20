"""Native cutoffs checked against independent exact mass and clock ledgers."""

import unittest
from dataclasses import replace
from decimal import Decimal, localcontext
from fractions import Fraction
from math import inf, ldexp, nextafter, ulp
from sys import float_info

from rocket_propulsion.core.errors import DomainError, InputError
from rocket_propulsion.integrations import lunaris as native
from rocket_propulsion.propulsion.burns import (
    PropagatorPropulsionBridge,
    PropulsionTimeline,
    ScheduledBurn,
    TabulatedBurnArtifact,
    TabulatedBurnSegment,
)


def model(rows, dry=None):
    segments = []
    for start, end, q0, q1 in rows:
        segments.append(TabulatedBurnSegment(start, end, 1.0, 1.0, 1.0, 1.0, q0, q1))
    artifact = TabulatedBurnArtifact(tuple(segments), "synthetic reserve oracle", "9" * 64)
    return PropagatorPropulsionBridge(
        artifact, protected_dry_mass_kg=dry, require_complete_stream_binding=False
    )


def mass_used(rows, time, scale=1.0):
    # Independent Bernstein endpoint integration, with exact input conversion.
    t = Fraction(time)
    total = Fraction(0)
    for start, end, q0, q1 in rows:
        a, b, q, r = map(Fraction, (start, end, q0, q1))
        if t <= a:
            break
        x = min((t - a) / (b - a), Fraction(1))
        total += (b - a) * (q * (x - x * x / 2) + r * x * x / 2)
    return total * Fraction(scale)


def timeline_for(value, ignition=0.0):
    return PropulsionTimeline((ScheduledBurn("burn", value, ignition),))


class ZeroDynamics:
    def build_rhs(self):
        return lambda t, y: [0.0] * 7


class NativeReserveAttacks(unittest.TestCase):
    def test_exact_inventory_clamp_policy_is_required_by_bridge_documents(self):
        value = model(((0.0, 1.0, 1.0, 1.0),), 1.5)
        payload = value.to_dict()
        policy = "exact_binary_state_and_limit_plus_declared_margin_v1"
        self.assertEqual(payload.get("inventory_clamp_numerics"), policy)
        self.assertEqual(PropagatorPropulsionBridge.from_dict(payload).model_hash, value.model_hash)
        for operation in ("missing", "changed"):
            changed = dict(payload)
            if operation == "missing":
                changed.pop("inventory_clamp_numerics")
            else:
                changed["inventory_clamp_numerics"] = "unknown_policy"
            with self.subTest(operation=operation), self.assertRaises(InputError):
                PropagatorPropulsionBridge.from_dict(changed)

    def test_public_initial_mass_normalization_matches_the_dynamics_ledger(self):
        initial = Fraction(2**53) + Fraction(3, 4)
        binary_initial = float(initial)
        dry = binary_initial - 1.0
        timeline = timeline_for(model(((0.0, 2.0, 1.0, 1.0),), dry))
        from_fraction = native.lunaris_burn_windows(timeline, initial_mass_kg=initial)[0]
        from_binary = native.lunaris_burn_windows(timeline, initial_mass_kg=binary_initial)[0]
        self.assertEqual(from_fraction, from_binary)
        dynamics = native.LunarisFiniteBurnDynamics(
            ZeroDynamics(), timeline, lambda *args: (1.0, 0.0, 0.0), initial_mass_kg=initial
        )
        self.assertGreaterEqual(dynamics.predicted_mass_kg(2.0), dry)
        self.assertEqual(dynamics.burn_windows[0], from_binary)

    def test_lattice_inverse_returns_the_last_safe_signed_float(self):
        from rocket_propulsion.integrations import native_mass

        for lower, upper, target in (
            (-float_info.max, float_info.max, Fraction(1, 3)),
            (-1e9, 1e9, Fraction(-1, 3)),
            (-1.0, 0.0, Fraction(-1, 3)),
            (0.0, 1.0, Fraction(1, 3)),
        ):
            calls = []

            def safe(date):
                calls.append(date)
                return Fraction(date) <= target

            cut = native_mass._greatest_safe(lower, upper, safe)
            self.assertLessEqual(Fraction(cut), target)
            self.assertGreater(Fraction(nextafter(cut, inf)), target)
            self.assertLessEqual(len(calls), 65)

    def test_base_failure_does_not_initialize_the_segment_clock(self):
        class FailsOnce:
            calls = 0

            def build_rhs(self):
                def rhs(t, y):
                    self.calls += 1
                    if self.calls == 1:
                        raise RuntimeError("synthetic base failure")
                    return [0.0] * 7

                return rhs

        artifact = TabulatedBurnArtifact(
            (
                TabulatedBurnSegment(0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.25, 0.25),
                TabulatedBurnSegment(1.0, 2.0, 2.0, 2.0, 2.0, 2.0, 0.25, 0.25),
            ),
            "base failure oracle",
            "6" * 64,
        )
        bridge = PropagatorPropulsionBridge(artifact, require_complete_stream_binding=False)
        dynamics = native.LunarisFiniteBurnDynamics(
            FailsOnce(), timeline_for(bridge), lambda *args: (1.0, 0.0, 0.0), initial_mass_kg=3.0
        )
        rhs = dynamics.build_rhs()
        with self.assertRaisesRegex(RuntimeError, "synthetic base failure"):
            rhs(0.25, [0.0] * 6 + [2.9375])
        self.assertEqual(rhs(1.0, [0.0] * 6 + [2.75])[3], 2.0 / 2.75)

    def test_initial_consistency_uses_unrounded_mass_difference(self):
        initial = ldexp(1.0, 53)
        timeline = timeline_for(model(((0.0, 1.0, 0.25, 0.25),)))
        dynamics = native.LunarisFiniteBurnDynamics(
            ZeroDynamics(),
            timeline,
            lambda *args: (1.0, 0.0, 0.0),
            initial_mass_kg=initial,
            mass_tolerance=ldexp(1.0, -56),
        )
        self.assertEqual(dynamics.predicted_mass_kg(1.0), initial)
        self.assertGreater(Fraction(1, 4), Fraction(ldexp(1.0, -56)) * Fraction(initial))
        with self.assertRaises(DomainError):
            dynamics.build_rhs()(1.0, [0.0] * 6 + [initial])
        tolerant = native.LunarisFiniteBurnDynamics(
            ZeroDynamics(),
            timeline,
            lambda *args: (1.0, 0.0, 0.0),
            initial_mass_kg=initial,
            mass_tolerance=ldexp(1.0, -54),
        )
        self.assertEqual(tolerant.build_rhs()(1.0, [0.0] * 6 + [initial])[6], 0.0)

    def test_signed_epoch_cutoffs_bracket_consumption_and_publish_exact_coordinates(self):
        rows = ((0.0, 1.0, 1.0, 1.0),)
        target = Fraction(0.5) - Fraction(1e-12)
        for ignition in (-1e9, -0.375, 0.0, 0.375, 1e9):
            with self.subTest(ignition=ignition):
                window = native.lunaris_burn_windows(
                    timeline_for(model(rows, 1.5), ignition), initial_mass_kg=2.0
                )[0]
                point = Fraction(window.depletion_time_s) - Fraction(ignition)
                used = mass_used(rows, point)
                self.assertEqual(Fraction(*window.cutoff_profile_time_ratio), point)
                self.assertLessEqual(used, target)
                self.assertGreater(
                    mass_used(
                        rows, Fraction(nextafter(window.depletion_time_s, inf)) - Fraction(ignition)
                    ),
                    target,
                )
                self.assertGreaterEqual(
                    Fraction(window.depletion_mass_shortfall_bound_kg), target - used
                )
                self.assertEqual(window.mass_at_cutoff_kg, float(Fraction(2) - used))
                self.assertEqual(
                    window.cutoff_numerics, "exact_mass_ledger_guarded_epoch_cutoff_v1"
                )

    def test_outward_mass_deficit_bounds_never_round_inwards(self):
        from rocket_propulsion.integrations import native_mass

        tiny = Fraction(nextafter(0.0, 1.0))
        for exact in (
            Fraction(0),
            tiny / 2,
            tiny * Fraction(5, 4),
            Fraction(1) + Fraction(1, 2**54),
            Fraction(7, 3),
        ):
            with self.subTest(exact=exact):
                bound = native_mass._upper_bound(exact)
                self.assertGreaterEqual(Fraction(bound), exact)
                if exact > 0:
                    self.assertLess(Fraction(nextafter(bound, -inf)), exact)
        for exact in (Fraction(-1), Fraction(float_info.max) * 2):
            with self.assertRaises(DomainError):
                native_mass._upper_bound(exact)

    def test_custom_zero_flow_scale_is_a_valid_closed_flow_control(self):
        value = model(((0.0, 1.0, 1.0, 1.0),), 1.5)
        value = replace(
            value,
            parameters=tuple(
                replace(p, minimum=0.0) if p.name == "mass_flow_scale" else p
                for p in value.parameters
            ),
        )
        dynamics = native.LunarisFiniteBurnDynamics(
            ZeroDynamics(),
            timeline_for(value),
            lambda *args: (1.0, 0.0, 0.0),
            initial_mass_kg=3.0,
            parameter_overrides={"burn": {"mass_flow_scale": 0.0}},
        )
        self.assertIsNone(dynamics.burn_windows[0].depletion_time_s)
        self.assertEqual(dynamics.predicted_mass_kg(1.0), 3.0)
        self.assertEqual(dynamics.build_rhs()(0.5, [0.0] * 6 + [3.0])[6], 0.0)

    def test_finite_cut_prefix_survives_a_nominal_end_beyond_clock_range(self):
        ignition = float_info.max / 2
        rows = ((0.0, float_info.max, 1e-308, 1e-308),)
        window = native.lunaris_burn_windows(
            timeline_for(model(rows, 1.5), ignition), initial_mass_kg=2.0
        )[0]
        cut = window.depletion_time_s
        self.assertIsNotNone(cut)
        self.assertLess(cut, float_info.max)
        used = mass_used(rows, Fraction(cut) - Fraction(ignition))
        self.assertLessEqual(used, Fraction(0.5) - Fraction(1e-12))
        self.assertGreater(
            mass_used(rows, Fraction(nextafter(cut, inf)) - Fraction(ignition)),
            Fraction(0.5) - Fraction(1e-12),
        )

    def test_actual_depletion_beyond_finite_clock_range_is_rejected(self):
        rows = ((0.0, float_info.max, 1e-308, 1e-308),)
        with self.assertRaises(DomainError):
            native.lunaris_burn_windows(
                timeline_for(model(rows, 1.5), float_info.max / 2), initial_mass_kg=3.0
            )

    def test_rejected_first_state_cannot_disable_consistency_validation(self):
        class CountingDynamics:
            calls = 0

            def build_rhs(self):
                def rhs(t, y):
                    self.calls += 1
                    return [0.0] * 7

                return rhs

        base = CountingDynamics()
        dynamics = native.LunarisFiniteBurnDynamics(
            base,
            timeline_for(model(((0.0, 1.0, 1.0, 1.0),))),
            lambda *args: (1.0, 0.0, 0.0),
            initial_mass_kg=3.0,
        )
        rhs = dynamics.build_rhs()
        for attempt in range(2):
            with self.subTest(attempt=attempt), self.assertRaises(DomainError):
                rhs(0.25, [0.0] * 6 + [2.0])
        self.assertEqual(base.calls, 0)
        self.assertEqual(rhs(0.25, [0.0] * 6 + [2.75])[6], -1.0)

    def test_failed_call_does_not_choose_the_first_successful_boundary_side(self):
        first = TabulatedBurnSegment(0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.25, 0.25)
        second = TabulatedBurnSegment(1.0, 2.0, 2.0, 2.0, 2.0, 2.0, 0.25, 0.25)
        artifact = TabulatedBurnArtifact((first, second), "restart oracle", "7" * 64)
        value = PropagatorPropulsionBridge(artifact, require_complete_stream_binding=False)
        dynamics = native.LunarisFiniteBurnDynamics(
            ZeroDynamics(), timeline_for(value), lambda *args: (1.0, 0.0, 0.0), initial_mass_kg=3.0
        )
        rhs = dynamics.build_rhs()
        with self.assertRaises(DomainError):
            rhs(0.25, [0.0] * 6 + [2.0])
        # The first accepted evaluation starts a new segment at the jump.
        self.assertEqual(rhs(1.0, [0.0] * 6 + [2.75])[3], 2.0 / 2.75)

    def test_window_copies_nested_inputs_and_rejects_inconsistent_records(self):
        dates = [0.0, 1.0]
        ratio = [2, 2]
        window = native.LunarisBurnWindow(
            "burn", 0.0, 1.0, 2.0, None, dates, cutoff_profile_time_ratio=ratio
        )
        dates[1] = 5.0
        ratio[0] = 9
        self.assertEqual(window.time_events_s, (0.0, 1.0))
        self.assertEqual(window.cutoff_profile_time_ratio, (1, 1))
        for change in (
            {"cutoff_time_s": -1.0},
            {"time_events_s": {0.0: "zero", 1.0: "one"}},
            {"time_events_s": (0.0, 0.0, 1.0)},
            {"depletion_time_s": 0.5},
            {"mass_at_cutoff_kg": 0.0},
            {"cutoff_profile_time_ratio": (True, 1)},
            {"depletion_mass_shortfall_bound_kg": 0.1},
            {"cutoff_numerics": "unsupported"},
        ):
            with self.subTest(change=change), self.assertRaises(DomainError):
                replace(window, **change)

    def test_parameter_overrides_are_frozen_with_the_precomputed_ledger(self):
        value = model(((0.0, 1.0, 1.0, 1.0),))
        timeline = timeline_for(value)
        overrides = {"burn": {"mass_flow_scale": 1.0}}
        dynamics = native.LunarisFiniteBurnDynamics(
            ZeroDynamics(),
            timeline,
            lambda *args: (1.0, 0.0, 0.0),
            initial_mass_kg=3.0,
            parameter_overrides=overrides,
        )
        overrides["burn"]["mass_flow_scale"] = 2.0
        self.assertEqual(dynamics.predicted_mass_kg(0.5), 2.5)
        self.assertEqual(dynamics.build_rhs()(0.5, [0.0] * 6 + [2.5])[6], -1.0)

    def test_rhs_rejects_invalid_mass_and_time_before_calling_base(self):
        class NeverCalled:
            def build_rhs(self):
                def rhs(t, y):
                    raise AssertionError("Invalid state reached base dynamics")

                return rhs

        dynamics = native.LunarisFiniteBurnDynamics(
            NeverCalled(),
            timeline_for(model(((0.0, 1.0, 1.0, 1.0),))),
            lambda *args: (1.0, 0.0, 0.0),
        )
        for time, mass in (
            (True, 2.0),
            ("0", 2.0),
            (inf, 2.0),
            (0.0, True),
            (0.0, "2"),
            (0.0, 0.0),
            (0.0, -1.0),
            (0.0, float("nan")),
        ):
            with self.subTest(time=time, mass=mass), self.assertRaises(DomainError):
                dynamics.build_rhs()(time, [0.0] * 6 + [mass])

    def test_constant_flow_public_cutoffs_survive_square_range_loss(self):
        for q, h in ((1.0, 1.0), (1e308, 1e-308), (1e-200, 1e200)):
            rows = ((0.0, h, q, q),)
            value = model(rows, 1.5)
            window = native.lunaris_burn_windows(timeline_for(value), initial_mass_kg=2.0)[0]
            self.assertIsNotNone(window.depletion_time_s)
            actual = window.depletion_time_s
            reference = []
            for precision in (1200, 1600):
                with localcontext() as context:
                    context.prec = precision
                    reference.append(
                        float((Decimal(".5") - Decimal.from_float(1e-12)) / Decimal.from_float(q))
                    )
            self.assertEqual(reference[0], reference[1])
            self.assertLessEqual(abs(actual - reference[0]), ulp(reference[0]))
            budget = Fraction(0.5) - Fraction(1e-12)
            self.assertLessEqual(mass_used(rows, actual), budget)
            self.assertGreater(mass_used(rows, nextafter(actual, inf)), budget)

    def test_increasing_and_decreasing_flows_bracket_the_exact_target(self):
        for exponent in (-1000, -600, -200, 0, 200, 600, 1000):
            q = ldexp(1.0, exponent)
            h = ldexp(1.0, -exponent)
            for a, b in ((0.0, 2.0), (2.0, 0.0), (1.0, 2.0), (2.0, 1.0)):
                rows = ((0.0, h, a * q, b * q),)
                value = model(rows, 1.5)
                with self.subTest(exponent=exponent, a=a, b=b):
                    cut = native.lunaris_burn_windows(timeline_for(value), initial_mass_kg=2.0)[
                        0
                    ].depletion_time_s
                    self.assertIsNotNone(cut)
                    target = Fraction(0.5) - Fraction(1e-12)
                    self.assertLessEqual(mass_used(rows, cut), target)
                    self.assertGreater(mass_used(rows, nextafter(cut, inf)), target)

    def test_partial_consumption_matches_paired_decimal_integrals(self):
        for q, h in ((1e308, 1e-308), (1e-200, 1e200), (3.0, 0.125)):
            rows = ((0.0, h, q, q / 2),)
            value = model(rows)
            for x in (0.25, 0.5, 1.0):
                t = x * h
                expected = []
                for precision in (1200, 1600):
                    with localcontext() as context:
                        context.prec = precision
                        d, f, g, tau = map(Decimal.from_float, (h, q, q / 2, t))
                        expected.append(float(f * tau + (g - f) * tau * tau / (2 * d)))
                self.assertEqual(expected[0], expected[1])
                self.assertEqual(native._consumed_by(value, 1.0, t), expected[0])

    def test_grouped_consumption_is_scaled_before_scalar_export(self):
        rows = ((0.0, 1.0, 1e-308, 1e-308),)
        value = model(rows)
        self.assertEqual(native._consumed_by(value, 1e308, 0.5), float(mass_used(rows, 0.5, 1e308)))
        tiny = nextafter(0.0, 1.0)
        # An unscaled partial is below binary64 range; the final scaled value is not.
        self.assertGreater(native._consumed_by(value, 1e308, tiny), 0.0)

    def test_unrepresentable_positive_cut_does_not_become_zero_success(self):
        h = 2 * nextafter(0.0, 1.0)
        value = model(((0.0, h, 1e308, 1e308),), 1e-12 - 2.5e-16)
        with self.assertRaises(DomainError):
            native.lunaris_burn_windows(timeline_for(value), initial_mass_kg=2e-12)

    def test_first_crossing_precedes_a_zero_flow_plateau(self):
        margin = 1e-12
        rows = ((0.0, 1.0, margin, margin), (1.0, 2.0, 0.0, 0.0), (2.0, 3.0, margin, margin))
        value = model(rows, 2 * margin)
        window = native.lunaris_burn_windows(timeline_for(value), initial_mass_kg=4 * margin)[0]
        self.assertEqual(window.depletion_time_s, 1.0)
        self.assertEqual(window.time_events_s, (0.0, 1.0))

    def test_absolute_cutoff_is_safe_on_the_actual_epoch_grid(self):
        ignition = 1e9
        rows = ((0.0, 1.0, 1.0, 1.0),)
        value = model(rows, 1.5)
        window = native.lunaris_burn_windows(timeline_for(value, ignition), initial_mass_kg=2.0)[0]
        cut = window.depletion_time_s
        relative = Fraction(cut) - Fraction(ignition)
        target = Fraction(0.5) - Fraction(1e-12)
        self.assertLessEqual(mass_used(rows, relative), target)
        self.assertGreater(
            mass_used(rows, Fraction(nextafter(cut, inf)) - Fraction(ignition)), target
        )
        self.assertGreaterEqual(
            window.depletion_mass_shortfall_bound_kg, float(target - mass_used(rows, relative))
        )

    def test_unused_unresolvable_tail_does_not_reject_a_valid_cut_prefix(self):
        rows = ((0.0, 1.0, 1.0, 1.0), (1.0, 1.00000001, 1.0, 1.0), (1.00000001, 2.0, 1.0, 1.0))
        value = model(rows, 1.5)
        timeline = timeline_for(value, 1e9)
        window = native.lunaris_burn_windows(timeline, initial_mass_kg=2.0)[0]
        self.assertLess(window.depletion_time_s, 1e9 + 1.0)
        dynamics = native.LunarisFiniteBurnDynamics(
            ZeroDynamics(), timeline, lambda *args: (1.0, 0.0, 0.0), initial_mass_kg=2.0
        )
        self.assertGreater(dynamics.build_rhs()(1e9 + 0.25, [0.0] * 6 + [1.75])[3], 0.0)

    def test_sequential_mass_ledger_retains_small_accumulated_consumption(self):
        value = model(((0.0, 1.0, 0.25, 0.25),))
        timeline = PropulsionTimeline(
            tuple(ScheduledBurn(str(i), value, float(i)) for i in range(4))
        )
        initial = ldexp(1.0, 53)
        windows = native.lunaris_burn_windows(timeline, initial_mass_kg=initial)
        self.assertEqual(windows[3].mass_at_ignition_kg, float(Fraction(initial) - Fraction(3, 4)))
        dynamics = native.LunarisFiniteBurnDynamics(
            ZeroDynamics(), timeline, lambda *args: (1.0, 0.0, 0.0), initial_mass_kg=initial
        )
        self.assertEqual(dynamics.predicted_mass_kg(4.0), initial - 1.0)

    def test_overlapping_unprotected_burns_use_mass_at_actual_ignition(self):
        value = model(((0.0, 2.0, 1.0, 1.0),))
        timeline = PropulsionTimeline(
            (ScheduledBurn("a", value, 0.0), ScheduledBurn("b", value, 1.0)),
            allow_overlapping_burns=True,
        )
        windows = native.lunaris_burn_windows(timeline, initial_mass_kg=10.0)
        self.assertEqual(windows[1].mass_at_ignition_kg, 9.0)
        dynamics = native.LunarisFiniteBurnDynamics(
            ZeroDynamics(), timeline, lambda *args: (1.0, 0.0, 0.0), initial_mass_kg=10.0
        )
        self.assertEqual(dynamics.predicted_mass_kg(1.5), 8.0)
        self.assertEqual(dynamics.predicted_mass_kg(3.0), 6.0)

    def test_nonpositive_mass_history_is_rejected_before_propagation(self):
        value = model(((0.0, 2.0, 1.0, 1.0),))
        with self.assertRaises(DomainError):
            native.lunaris_burn_windows(timeline_for(value), initial_mass_kg=2.0)

    def test_exact_clamp_comparison_does_not_round_a_limit_up_to_the_state(self):
        dry = 8192.0
        mass = nextafter(dry, inf)
        value = model(((0.0, 1.0, 1.0, 1.0),), dry)
        self.assertGreater(Fraction(mass), Fraction(dry) + Fraction(1e-12))
        result = value.evaluate(
            relative_time_s=0.5, vehicle_mass_kg=mass, direction=(1.0, 0.0, 0.0)
        )
        self.assertTrue(result.active)

    def test_invalid_private_scalar_inputs_are_typed(self):
        value = model(((0.0, 1.0, 1.0, 1.0),))
        for invalid in (True, "1", -1.0, inf, float("nan")):
            with self.subTest(invalid=invalid), self.assertRaises(DomainError):
                native._consumed_by(value, invalid, 0.5)

    def test_zero_flow_no_crossing_and_already_limited_controls(self):
        self.assertEqual(native._consumed_by(model(((0.0, 1.0, 1.0, 1.0),)), 0.0, 0.5), 0.0)
        value = model(((0.0, 1.0, 0.0, 0.0),), 1.5)
        window = native.lunaris_burn_windows(timeline_for(value), initial_mass_kg=2.0)[0]
        self.assertIsNone(window.depletion_time_s)
        self.assertEqual(window.cutoff_time_s, 1.0)
        limited = native.lunaris_burn_windows(timeline_for(value), initial_mass_kg=1.5)[0]
        self.assertEqual(limited.depletion_time_s, 0.0)

    def test_ordinary_public_dynamics_and_diagnostic_mass_control(self):
        value = model(((0.0, 1.0, 1.0, 1.0),), 1.5)
        timeline = timeline_for(value)
        dynamics = native.LunarisFiniteBurnDynamics(
            ZeroDynamics(), timeline, lambda *args: (1.0, 0.0, 0.0), initial_mass_kg=2.0
        )
        self.assertEqual(dynamics.predicted_mass_kg(0.25), 1.75)
        self.assertGreaterEqual(
            dynamics.predicted_mass_kg(1.0), float(Fraction(1.5) + Fraction(1e-12))
        )
        cut = dynamics.burn_windows[0].depletion_time_s
        self.assertEqual(
            dynamics.build_rhs()(cut, [0.0] * 6 + [dynamics.predicted_mass_kg(cut)])[3], 0.0
        )
