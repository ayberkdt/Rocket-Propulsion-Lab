"""Independent L1 exposure/root, boundary and public burn-ledger attacks."""

import unittest
from dataclasses import replace
from decimal import Decimal, localcontext
from math import nextafter, ulp

from rocket_propulsion.core.errors import DomainError
from rocket_propulsion.propulsion.burns import (
    BurnTargetKind,
    CutoffReason,
    ProfiledPerformanceProvider,
    ThrottleSchedule,
    ThrottleSegment,
    burn_result_from_json,
    burn_result_to_json,
    integrated_throttle,
    simulate_burn,
    throttle_first_moment,
    time_for_integrated_throttle,
)
from rocket_propulsion.propulsion.burns.profiles import (
    ScheduleInterval,
    resolved_schedule_intervals,
)
from tests.test_burn_constant_numerics import inverse_series, mono


def schedule(*segments):
    return ThrottleSchedule(tuple(ThrottleSegment(*values) for values in segments))


def run(definition, point, command=None):
    command = command or schedule((10.0, 1.0, 1.0))
    return simulate_burn(definition, ProfiledPerformanceProvider(point, command))


def decimal_integrals(command, time, precision):
    with localcontext() as context:
        context.prec = precision
        stop = Decimal.from_float(time)
        cursor = exposure = moment = Decimal(0)
        for segment in command.segments:
            duration, a, b = map(
                Decimal.from_float,
                (segment.duration_s, segment.start_throttle, segment.end_throttle),
            )
            used = min(duration, max(Decimal(0), stop - cursor))
            local = a * used + (b - a) * used * used / (2 * duration)
            exposure += local
            moment += cursor * local + a * used * used / 2 + (b - a) * used**3 / (3 * duration)
            cursor += duration
        return float(exposure), float(moment)


def inverse_bisection(command, exposure, precision):
    with localcontext() as context:
        context.prec = precision
        target = Decimal.from_float(exposure)
        pieces = [
            tuple(map(Decimal.from_float, (s.duration_s, s.start_throttle, s.end_throttle)))
            for s in command.segments
        ]
        left = Decimal(0)
        right = sum((p[0] for p in pieces), Decimal(0))
        for _ in range(4096):
            middle = (left + right) / 2
            elapsed = value = Decimal(0)
            for duration, a, b in pieces:
                used = min(duration, max(Decimal(0), middle - elapsed))
                value += a * used + (b - a) * used * used / (2 * duration)
                elapsed += duration
            if value >= target:
                right = middle
            else:
                left = middle
            if float(left) == float(right):
                return float(right)
        raise AssertionError("independent inverse did not settle")


class ProfileBurnNumericalAttacks(unittest.TestCase):
    def close(self, actual, expected, count=4):
        self.assertLessEqual(abs(actual - expected), count * ulp(expected), (actual, expected))

    def test_positive_tiny_exposure_is_not_replaced_by_zero(self):
        command = schedule((1.0, 1.0, 1.0))
        for exposure in (1e-14, 1e-100, 1e-300):
            with self.subTest(exposure=exposure):
                self.assertEqual(time_for_integrated_throttle(command, exposure), exposure)
                self.assertEqual(integrated_throttle(command, exposure), exposure)

    def test_increasing_and_decreasing_roots_match_independent_bisection(self):
        for command in (
            schedule((1.0, 0.0, 1.0)),
            schedule((1.0, 1.0, 0.0)),
            schedule((2.0, 0.125, 0.875)),
            schedule((2.0, 0.875, 0.125)),
        ):
            for exposure in (1e-300, 1e-14, 0.125, 0.375):
                with self.subTest(command=command, exposure=exposure):
                    reference = inverse_bisection(command, exposure, 300)
                    self.assertEqual(reference, inverse_bisection(command, exposure, 200))
                    self.close(time_for_integrated_throttle(command, exposure), reference)

    def test_reachable_total_neighbors_and_idle_plateau_use_earliest_time(self):
        command = schedule((1.0, 1.0, 1.0), (2.0, 0.0, 0.0), (1.0, 1.0, 1.0))
        self.assertEqual(time_for_integrated_throttle(command, 1.0), 1.0)
        # The first positive increment after the plateau rounds to 3.0;
        # another exposure ulp gives a distinguishable time increment.
        self.assertEqual(time_for_integrated_throttle(command, nextafter(1.0, 2.0)), 3.0)
        self.assertGreater(
            time_for_integrated_throttle(command, nextafter(nextafter(1.0, 2.0), 2.0)), 3.0
        )
        self.assertEqual(time_for_integrated_throttle(command, 2.0), 4.0)
        self.assertIsNone(time_for_integrated_throttle(command, nextafter(2.0, 3.0)))
        self.assertEqual(time_for_integrated_throttle(command, 0.0), 0.0)

    def test_small_real_slope_over_a_long_segment_is_not_treated_as_constant(self):
        command = schedule((1e20, 0.25, 0.75))
        exposure = 2.5e19
        reference = inverse_bisection(command, exposure, 300)
        self.assertEqual(reference, inverse_bisection(command, exposure, 200))
        self.close(time_for_integrated_throttle(command, exposure), reference)

    def test_integrals_and_first_moments_match_independent_polynomials(self):
        command = schedule(
            (0.5, 0.0, 0.75), (1.25, 0.75, 0.125), (0.25, 0.0, 0.0), (0.75, 0.25, 1.0)
        )
        for time in (0.125, 0.5, nextafter(0.5, 1.0), 1.75, 2.0, 2.5, 3.0):
            with self.subTest(time=time):
                reference = decimal_integrals(command, time, 300)
                self.assertEqual(reference, decimal_integrals(command, time, 200))
                self.close(integrated_throttle(command, time), reference[0])
                self.close(throttle_first_moment(command, time), reference[1])

    def test_positive_segment_exposure_survives_early_half_underflow(self):
        tiny = float.fromhex("0x0.0000000000001p-1022")
        command = schedule((1e308, 0.0, tiny))
        expected = float(Decimal.from_float(tiny) * Decimal.from_float(1e308) / 2)
        self.assertEqual(command.segments[0].exposure_s, expected)
        self.assertEqual(command.total_exposure_s, expected)
        self.assertEqual(integrated_throttle(command, 1e308), expected)

    def test_unexportable_total_duration_does_not_prevent_early_finite_burn(self):
        command = schedule((1e308, 1.0, 1.0), (1e308, 1.0, 1.0))
        with self.assertRaises(DomainError):
            _ = command.total_duration_s
        definition, point = mono()
        result = run(definition, point, command)
        self.assertEqual(result.summary.duration_s, 1.0)
        self.assertTrue(result.summary.target_achieved)

    def test_true_moment_overflow_is_typed_but_centroid_can_remain_finite(self):
        command = schedule((1e200, 1.0, 1.0))
        with self.assertRaises(DomainError):
            throttle_first_moment(command, 1e200)
        definition, point = mono(flow=1e-200, force=1e-200, target=1e200, maximum=1e200)
        result = run(definition, point, command)
        self.close(result.summary.thrust_centroid_time_s, 5e199)
        self.close(result.summary.consumed_propellant_kg, 1.0)

    def test_malformed_and_boolean_schedule_queries_raise_typed_errors(self):
        command = schedule((1.0, 1.0, 1.0))
        for value in (None, True, "1", 1j, float("nan"), float("inf"), -1.0):
            for function in (
                integrated_throttle,
                throttle_first_moment,
                time_for_integrated_throttle,
            ):
                with (
                    self.subTest(function=function.__name__, value=value),
                    self.assertRaises(DomainError),
                ):
                    function(command, value)

    def test_small_public_delta_v_targets_match_independent_series(self):
        for target in (1e-10, 1e-14, 1e-100, 1e-300):
            with self.subTest(target=target):
                definition, point = mono(kind=BurnTargetKind.IDEAL_DELTA_V, target=target)
                reference = inverse_series(definition, point, 300)
                self.assertEqual(reference, inverse_series(definition, point, 200))
                result = run(definition, point)
                self.close(result.summary.duration_s, reference)
                self.close(result.summary.ideal_delta_v_m_s, target)
                self.assertTrue(result.summary.target_achieved)
                self.assertEqual(burn_result_from_json(burn_result_to_json(result)), result)

    def test_tiny_target_after_leading_idle_preserves_consumption(self):
        definition, point = mono(kind=BurnTargetKind.PROPELLANT_MASS, target=1e-14)
        result = run(definition, point, schedule((2.0, 0.0, 0.0), (1.0, 1.0, 1.0)))
        self.close(result.summary.duration_s, 2.0 + 1e-14)
        self.assertEqual(result.summary.consumed_propellant_kg, 1e-14)
        self.assertEqual(dict(result.summary.consumed_by_tank_kg)["tank"], 1e-14)
        self.assertTrue(result.summary.target_achieved)

    def test_profile_maximum_cutoff_is_not_overridden_by_target_tolerance(self):
        for kind in (BurnTargetKind.DURATION, BurnTargetKind.PROPELLANT_MASS):
            for maximum, reached in (
                (nextafter(1.0, 0.0), False),
                (1.0, True),
                (nextafter(1.0, 2.0), True),
            ):
                with self.subTest(maximum=maximum, kind=kind):
                    definition, point = mono(maximum=maximum, kind=kind)
                    summary = run(definition, point).summary
                    self.assertEqual(summary.duration_s, min(1.0, maximum))
                    self.assertEqual(summary.target_achieved, reached)

    def test_profile_schedule_end_cannot_be_reported_as_neighboring_target_success(self):
        definition, point = mono(target=nextafter(1.0, 2.0))
        summary = run(definition, point, schedule((1.0, 1.0, 1.0))).summary
        self.assertEqual(summary.duration_s, 1.0)
        self.assertFalse(summary.target_achieved)
        self.assertEqual(summary.cutoff_reason, CutoffReason.SCHEDULE_END)

    def test_reserve_inversion_on_ramp_preserves_exact_inventory(self):
        definition, point = mono(target=2.0, loaded=1.0, reserve=0.25)
        result = run(definition, point, schedule((2.0, 0.0, 1.0)))
        self.assertEqual(result.summary.cutoff_reason, CutoffReason.TANK_RESERVE)
        self.assertFalse(result.summary.target_achieved)
        self.assertEqual(dict(result.final_state.tank_masses_kg)["tank"], 0.25)
        self.assertEqual(result.summary.consumed_propellant_kg, 0.75)
        self.assertEqual(result.profile[-1].tank_masses_kg, result.final_state.tank_masses_kg)

    def test_small_profile_consumption_is_not_lost_in_rounded_endpoints(self):
        definition, point = mono(target=1e-20)
        result = run(definition, point)
        self.assertEqual(result.summary.final_mass_kg, definition.initial_mass_kg)
        self.assertEqual(dict(result.summary.consumed_by_tank_kg)["tank"], 1e-20)
        self.assertEqual(result.summary.mass_closure_error_kg, 0.0)

    def test_profile_recovers_finite_velocity_after_intermediate_ratio_loss(self):
        for arguments, expected in (
            ({"flow": 1e-200, "force": 1e100, "mass": 1e200, "loaded": 1.0, "dry": 0.0}, 1e-100),
            ({"flow": 1e-308, "force": 5.0}, 0.05),
        ):
            with self.subTest(arguments=arguments):
                definition, point = mono(**arguments)
                self.close(run(definition, point).summary.ideal_delta_v_m_s, expected)

    def test_tiny_nonconstant_profile_cannot_claim_exact_constant_export(self):
        definition, point = mono()
        result = run(definition, point, schedule((1.0, 0.5, nextafter(0.5, 1.0))))
        self.assertFalse(result.summary.exact_constant_sidera_export_available)
        control = run(definition, point, schedule((1.0, 0.5, 0.5)))
        self.assertTrue(control.summary.exact_constant_sidera_export_available)

    def test_declared_flow_discrepancy_remains_visible(self):
        definition, point = mono(target=0.125)
        point = replace(point, streams=(replace(point.streams[0], mass_flow_kg_s=1.0 + 2.0**-34),))
        summary = run(definition, point).summary
        self.assertEqual(summary.mass_closure_error_kg, 2.0**-37)

    def test_unreachable_large_targets_clip_without_finite_exposure_precondition(self):
        for kind in (
            BurnTargetKind.PROPELLANT_MASS,
            BurnTargetKind.TOTAL_IMPULSE,
            BurnTargetKind.IDEAL_DELTA_V,
        ):
            with self.subTest(kind=kind):
                definition, point = mono(
                    flow=1e-200, force=1e-200, kind=kind, target=1e200, maximum=1.0
                )
                summary = run(definition, point).summary
                self.assertFalse(summary.target_achieved)
                self.assertEqual(summary.duration_s, 1.0)
                self.assertEqual(summary.cutoff_reason, CutoffReason.MAXIMUM_DURATION)

    def test_zero_inventory_and_zero_throttle_prefix_keep_finite_cutoff_fields(self):
        definition, point = mono(loaded=1.0, reserve=1.0)
        result = run(definition, point, schedule((1.0, 0.0, 0.0), (1.0, 1.0, 1.0)))
        self.assertEqual(result.summary.duration_s, 0.0)
        self.assertFalse(result.summary.target_achieved)
        self.assertEqual(result.summary.cutoff_reason, CutoffReason.TANK_RESERVE)
        self.assertEqual(result.summary.ideal_delta_v_m_s, 0.0)

    def test_schedule_records_reject_malformed_fields_before_arithmetic(self):
        for field, value in (
            ("duration_s", True),
            ("start_throttle", None),
            ("end_throttle", "1"),
            ("phase", 3),
        ):
            with self.subTest(field=field), self.assertRaises(DomainError):
                ThrottleSegment(
                    **{

                            "duration_s": 1.0, "start_throttle": 0.5, "end_throttle": 0.5, "phase": "steady"
                        ,
                        field: value,
                    }
                )
        with self.assertRaises(DomainError):
            ThrottleSchedule((object(),))
        with self.assertRaises(DomainError):
            ScheduleInterval(0.0, 1.0, True, 0.5, "steady")
        with self.assertRaises(DomainError):
            ScheduleInterval(1.0, 0.0, 0.5, 0.5, "steady")

    def test_collapsed_public_interval_is_explicit_while_integral_survives(self):
        command = schedule((1e16, 0.0, 0.0), (1.0, 1.0, 1.0))
        self.assertEqual(integrated_throttle(command, 1e16 + 2.0), 1.0)
        with self.assertRaises(DomainError) as caught:
            resolved_schedule_intervals(command)
        self.assertEqual(caught.exception.field, "Interval time resolution")

    def test_impulse_closure_bound_covers_rounded_output_time(self):
        definition, point = mono(target=2.0, loaded=1.0, reserve=0.25)
        result = run(definition, point, schedule((2.0, 0.0, 1.0)))
        errors = []
        for precision in (200, 300):
            with localcontext() as context:
                context.prec = precision
                time = Decimal.from_float(result.summary.duration_s)
                errors.append(abs(Decimal("7.5") - Decimal(10) * time * time / 4))
        self.assertEqual(float(errors[0]), float(errors[1]))
        bound = result.summary.impulse_closure_error_n_s
        self.assertGreaterEqual(Decimal.from_float(bound), errors[1])
        self.assertLessEqual(bound, nextafter(float(errors[1]), float("inf")))

    def test_public_interval_properties_recover_finite_moment_and_reject_true_overflow(self):
        interval = ScheduleInterval(2.0, 5.0, 0.25, 0.75, "linear")
        self.assertEqual(interval.duration_s, 3.0)
        self.assertEqual(interval.exposure_s, 1.5)
        self.assertEqual(interval.throttle_first_moment_s2, 5.625)
        tiny = float.fromhex("0x0.0000000000001p-1022")
        interval = ScheduleInterval(0.0, 1e308, 0.0, tiny, "tiny ramp")
        references = []
        for precision in (200, 300):
            with localcontext() as context:
                context.prec = precision
                duration, throttle = map(Decimal.from_float, (1e308, tiny))
                references.append(float(throttle * duration * duration / 3))
        self.assertEqual(*references)
        self.assertEqual(interval.throttle_first_moment_s2, references[0])
        with self.assertRaises(DomainError):
            _ = ScheduleInterval(
                1e308, nextafter(1e308, float("inf")), 1.0, 1.0, "overflow"
            ).throttle_first_moment_s2

    def test_true_physical_output_range_errors_have_field_labels(self):
        for arguments, label in (
            ({"flow": 1e-200, "force": 1e-200, "target": 1e-200}, "Consumed propellant"),
            (
                {
                    "flow": 1e100,
                    "force": 1e308,
                    "mass": 1e308,
                    "loaded": 1e308,
                    "dry": 0.0,
                    "target": 10.0,
                    "maximum": 10.0,
                },
                "Delivered impulse",
            ),
        ):
            with self.subTest(arguments=arguments), self.assertRaises(DomainError) as caught:
                definition, point = mono(**arguments)
                run(definition, point)
            self.assertIn(label, caught.exception.field or "")
