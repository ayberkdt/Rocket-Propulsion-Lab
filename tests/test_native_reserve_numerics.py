"""Native cutoffs checked against independent exact mass and clock ledgers."""

import unittest
from dataclasses import replace
from decimal import Decimal, localcontext
from fractions import Fraction
from math import inf, ldexp, nextafter, ulp

from rocket_propulsion.core.errors import DomainError
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
