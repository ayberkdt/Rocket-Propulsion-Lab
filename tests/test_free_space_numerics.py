"""Free-space integration attacks, independent exact/decimal analytic references."""

import unittest
from dataclasses import replace
from decimal import Decimal, localcontext
from fractions import Fraction
from math import inf, nextafter, ulp

from rocket_propulsion.core.errors import DomainError
from rocket_propulsion.propulsion.burns import burn_metrics as metrics_module
from rocket_propulsion.propulsion.burns import free_space_burn_metrics, gauss_legendre
from rocket_propulsion.propulsion.burns.propagator import PropagatorPropulsionBridge
from tests.test_tabulated_ledger_numerics import history, segment


def constant_reference(force, flow, duration, mass, precision):
    with localcontext() as context:
        context.prec = precision
        f, q, t, m = map(Decimal.from_float, (force, flow, duration, mass))
        consumed = q * t
        ratio = -(1 - consumed / m).ln() if q else Decimal(0)
        dv = f / q * ratio if q else f * t / m
        moment = f / q * (m / q * ratio - t) if q else f * t * t / (2 * m)
        return {
            "dv": dv,
            "centroid": moment / dv if dv else Decimal(0),
            "exhaust": f / q if q else None,
            "final": m - consumed,
            "consumed": consumed,
        }


class FreeSpaceNumericalAttacks(unittest.TestCase):
    def encloses(self, actual, reference):
        for field, bound, key in (
            ("delta_v_m_s", "delta_v_error_estimate_m_s", "dv"),
            ("impulsive_equivalent_time_s", "centroid_error_bound_s", "centroid"),
            (
                "effective_exhaust_velocity_m_s",
                "effective_exhaust_velocity_error_bound_m_s",
                "exhaust",
            ),
        ):
            if reference[key] is None:
                continue
            self.assertLessEqual(
                abs(Decimal.from_float(getattr(actual, field)) - reference[key]),
                Decimal.from_float(getattr(actual, bound)),
                (field, actual, reference),
            )
            if actual.quadrature_points == 32:
                self.close(getattr(actual, field), float(reference[key]))

    def close(self, actual, expected, ulps=8):
        self.assertLessEqual(abs(actual - expected), ulps * ulp(expected), (actual, expected))

    def test_small_acceleration_and_large_moment_are_rescued_before_export(self):
        for force, flow, duration, mass in (
            (1e-200, 1e-150, 1e150, 1e200),
            (1e-200, 1e-200, 1e200, 1e200),
            (1e200, 0.0, 1e200, 1e100),
        ):
            artifact = history(segment(end=duration, f0=force, q0=flow))
            a = constant_reference(force, flow, duration, mass, 800)
            b = constant_reference(force, flow, duration, mass, 1200)
            self.assertEqual(
                {k: float(v) if v is not None else None for k, v in a.items()},
                {k: float(v) if v is not None else None for k, v in b.items()},
            )
            leaves = tuple(
                metrics_module._mass_intervals(artifact, Fraction(mass), Fraction(1), Fraction(1))
            )
            raw_dv = sum((metrics_module._quadrature(leaf, 32)[0] for leaf in leaves), Fraction(0))
            self.close(float(raw_dv), float(b["dv"]))
            try:
                actual = free_space_burn_metrics(artifact, mass)
            except DomainError as error:
                self.fail(f"Known representable metric fixture was rejected: {error}")
            self.close(actual.delta_v_m_s, float(b["dv"]))
            self.close(actual.impulsive_equivalent_time_s, float(b["centroid"]))
            if b["exhaust"] is None:
                self.assertIsNone(actual.effective_exhaust_velocity_m_s)
            else:
                self.close(actual.effective_exhaust_velocity_m_s, float(b["exhaust"]))

    def test_quadrature_error_cannot_be_zero_for_a_wrong_positive_result(self):
        actual = free_space_burn_metrics(history(segment(end=1e150, f0=1e-200, q0=1e-150)), 1e200)
        reference = constant_reference(1e-200, 1e-150, 1e150, 1e200, 1200)["dv"]
        self.assertGreater(actual.delta_v_m_s, 0.0)
        error = abs(Decimal.from_float(actual.delta_v_m_s) - reference)
        self.assertLessEqual(error, Decimal.from_float(actual.delta_v_error_estimate_m_s))
        self.assertGreater(actual.delta_v_error_estimate_m_s, 0.0)

    def test_low_and_high_order_bounds_cover_closed_form_across_mass_ratios(self):
        for mass in (2.0, 1.01, nextafter(1.0, inf)):
            reference = constant_reference(2.0, 1.0, 1.0, mass, 300)["dv"]
            self.assertEqual(
                float(reference), float(constant_reference(2.0, 1.0, 1.0, mass, 500)["dv"])
            )
            for order in (1, 2, 16, 100):
                with self.subTest(mass=mass, order=order):
                    result = free_space_burn_metrics(
                        history(segment(f0=2.0)), mass, quadrature_points=order
                    )
                    self.assertLessEqual(
                        abs(Decimal.from_float(result.delta_v_m_s) - reference),
                        Decimal.from_float(result.delta_v_error_estimate_m_s),
                    )
                    self.assertGreater(result.impulsive_equivalent_time_s, 0.0)
                    self.assertLessEqual(result.impulsive_equivalent_time_s, 1.0)
                    self.assertEqual(result.quadrature_points, 2 * order)

    def test_zero_flow_and_zero_force_have_distinct_controls(self):
        interval = (
            Fraction(0),
            Fraction(2),
            Fraction(3),
            Fraction(3),
            Fraction(0),
            Fraction(0),
            Fraction(4),
        )
        quadrature_value, _ = metrics_module._quadrature(interval, 32)
        self.assertLessEqual(abs(quadrature_value - Fraction(3, 2)), Fraction(1, 10**90))
        result = free_space_burn_metrics(history(segment(end=2.0, f0=3.0, q0=0.0)), 4.0)
        self.assertEqual(result.delta_v_m_s, 1.5)
        self.assertEqual(result.impulsive_equivalent_time_s, 1.0)
        self.assertIsNone(result.effective_exhaust_velocity_m_s)
        zero = free_space_burn_metrics(history(segment(f0=0.0, q0=0.5)), 2.0)
        self.assertEqual(zero.delta_v_m_s, 0.0)
        self.assertEqual(zero.impulsive_equivalent_time_s, 0.0)
        self.assertEqual(zero.effective_exhaust_velocity_m_s, 0.0)

    def test_mass_decision_uses_unrounded_segment_consumption(self):
        # Trapezoid consumes 1 + 2**-53; a rounded consumption of 1 must
        # not permit the burn from exactly 1 kg.
        artifact = history(segment(q0=1.0, q1=nextafter(1.0, inf)))
        with self.assertRaises(DomainError):
            free_space_burn_metrics(artifact, 1.0)
        legal = free_space_burn_metrics(artifact, nextafter(1.0, inf))
        self.assertEqual(legal.final_mass_kg, ulp(1.0) / 2)

    def test_scaled_mass_flow_rescues_unrepresentable_unscaled_integral(self):
        artifact = history(segment(end=1e-100, f0=1.0, q0=1e-300))
        actual = free_space_burn_metrics(artifact, 10.0, mass_flow_scale=1e300)
        self.assertGreater(actual.consumed_propellant_kg, 0.0)
        self.close(actual.consumed_propellant_kg, 1e-100)
        self.close(actual.delta_v_m_s, 1e-101)
        self.close(actual.effective_exhaust_velocity_m_s, 1.0)

    def test_scaled_thrust_does_not_need_a_finite_intermediate_force(self):
        artifact = history(segment(end=1e-308, f0=1e308, q0=1e308))
        actual = free_space_burn_metrics(
            artifact, 1.7e308, thrust_scale=1e308, mass_flow_scale=1e308
        )
        with localcontext() as context:
            context.prec = 500
            q = Decimal.from_float(1e308) * Decimal.from_float(1e308)
            expected = float(
                -(1 - q * Decimal.from_float(1e-308) / Decimal.from_float(1.7e308)).ln()
            )
        self.close(actual.delta_v_m_s, expected)
        self.assertGreater(actual.impulsive_equivalent_time_s, 0.0)
        self.close(actual.effective_exhaust_velocity_m_s, 1.0)

    def test_true_unrepresentable_metric_outputs_fail_typed(self):
        cases = (
            (history(segment(end=0.25, f0=nextafter(0.0, 1.0), q0=0.0)), 1.0),
            (history(segment(end=2.0, f0=1e308, q0=0.0)), 1.0),
            (history(segment(end=1e100, f0=1e300, q0=1e-10)), 1e200),
        )
        for artifact, mass in cases:
            with self.subTest(mass=mass), self.assertRaises(DomainError) as caught:
                free_space_burn_metrics(artifact, mass)
            self.assertIsNotNone(caught.exception.field)

    def test_quadrature_count_validation_is_not_bypassed_by_warm_cache(self):
        gauss_legendre(2)
        for value in (True, False, 2.0, "2", None, 0, 201):
            with self.subTest(value=value), self.assertRaises(DomainError):
                gauss_legendre(value)
        for value in (True, 1.5, 0, 101):
            with self.subTest(metrics_order=value), self.assertRaises(DomainError):
                free_space_burn_metrics(history(segment()), 10.0, quadrature_points=value)

    def test_error_contract_and_partition_work_are_explicit(self):
        result = free_space_burn_metrics(history(segment()), nextafter(1.0, inf))
        self.assertEqual(
            result.error_estimate_kind, "directed_bernstein_enclosure_including_output_rounding"
        )
        self.assertGreater(result.integration_subintervals, 1)
        self.assertGreater(result.delta_v_error_estimate_m_s, 0.0)

    def test_ordinary_bridge_and_analytic_control_still_agree(self):
        artifact = history(segment(end=2.0, f0=12.0, q0=0.4))
        PropagatorPropulsionBridge(artifact, require_complete_stream_binding=False)
        actual = free_space_burn_metrics(artifact, 10.0)
        reference = constant_reference(12.0, 0.4, 2.0, 10.0, 200)
        self.close(actual.delta_v_m_s, float(reference["dv"]))
        self.close(actual.impulsive_equivalent_time_s, float(reference["centroid"]))
        self.close(actual.effective_exhaust_velocity_m_s, 30.0)
        self.assertEqual(actual.final_mass_kg, float(reference["final"]))

    def test_linear_force_constant_flow_bounds_cover_independent_log_moments(self):
        refs = []
        for precision in (200, 300):
            with localcontext() as context:
                context.prec = precision
                mass, q, t, f, s = map(Decimal, ("16", "1", "2", "3", "-1"))
                j0 = (mass / (mass - q * t)).ln() / q
                j1 = (mass * j0 - t) / q
                j2 = (mass * j1 - t * t / 2) / q
                dv = f * j0 + s * j1
                moment = f * j1 + s * j2
                refs.append(
                    {
                        "dv": dv,
                        "centroid": moment / dv,
                        "exhaust": dv / (mass / (mass - q * t)).ln(),
                    }
                )
        self.assertEqual(
            {k: float(v) for k, v in refs[0].items()}, {k: float(v) for k, v in refs[1].items()}
        )
        for order in (1, 16):
            self.encloses(
                free_space_burn_metrics(
                    history(segment(end=2.0, f0=3.0, f1=1.0)), 16.0, quadrature_points=order
                ),
                refs[1],
            )

    def test_increasing_flow_bounds_cover_independent_quadratic_mass_integrals(self):
        refs = []
        for precision in (200, 300):
            with localcontext() as context:
                context.prec = precision
                a, b, t, f, s = map(Decimal, ("4", "2", "2", "1", "2"))
                j0 = ((a + b) / (a - b)).ln() / (2 * a * b)
                j1 = (a * a / (a * a - b * b)).ln() / (2 * b * b)
                j2 = (a * a * j0 - 1) / (b * b)
                dv = t * (f * j0 + s * j1)
                moment = t * t * (f * j1 + s * j2) + dv  # main segment starts after an idle second
                refs.append(
                    {
                        "dv": dv,
                        "centroid": moment / dv,
                        "exhaust": dv / (a * a / (a * a - b * b)).ln(),
                    }
                )
        self.assertEqual(
            {k: float(v) for k, v in refs[0].items()}, {k: float(v) for k, v in refs[1].items()}
        )
        artifact = history(
            segment(f0=0.0, q0=0.0), segment(start=1.0, end=3.0, f0=1.0, f1=3.0, q0=0.0, q1=4.0)
        )
        for order in (1, 16):
            self.encloses(free_space_burn_metrics(artifact, 16.0, quadrature_points=order), refs[1])

    def test_decreasing_flow_bounds_cover_independent_square_mass_integrals(self):
        refs = []
        for precision in (200, 300):
            with localcontext() as context:
                context.prec = precision
                a, b, t, f, s = map(Decimal, ("4", "1", "2", "3", "-2"))
                k = f + s * a / b
                diff = 1 / (a - b) - 1 / a
                logarithm = (a / (a - b)).ln()
                dv = t / b * (k * diff - s / b * logarithm)
                moment = t * t / (b * b) * (a * k * diff - (k + a * s / b) * logarithm + s)
                refs.append({"dv": dv, "centroid": moment / dv, "exhaust": dv / (2 * logarithm)})
        self.assertEqual(
            {k: float(v) for k, v in refs[0].items()}, {k: float(v) for k, v in refs[1].items()}
        )
        for order in (1, 16):
            self.encloses(
                free_space_burn_metrics(
                    history(segment(end=2.0, f0=3.0, f1=1.0, q0=4.0, q1=3.0)),
                    16.0,
                    quadrature_points=order,
                ),
                refs[1],
            )

    def test_ambient_decimal_rounding_and_traps_do_not_change_metrics(self):
        from decimal import ROUND_CEILING, Context, Inexact

        artifact = history(segment(f0=3.0, q0=0.2))
        expected = free_space_burn_metrics(artifact, 10.0)
        metrics_module._bernstein_enclosure.cache_clear()
        hostile = Context(prec=2, rounding=ROUND_CEILING, Emax=5, Emin=-5)
        hostile.traps[Inexact] = True
        with localcontext(hostile):
            actual = free_space_burn_metrics(artifact, 10.0)
        self.assertEqual(actual, expected)

    def test_metric_record_rejects_malformed_scalars_counts_and_mass_order(self):
        base = free_space_burn_metrics(history(segment()), 10.0)
        for name, value in (
            ("delta_v_m_s", True),
            ("centroid_error_bound_s", float("nan")),
            ("quadrature_points", 1.5),
            ("integration_subintervals", True),
            ("final_mass_kg", 11.0),
            ("delta_v_error_estimate_m_s", -1.0),
        ):
            with self.subTest(field=name), self.assertRaises(DomainError):
                replace(base, **{name: value})

    def test_exact_mass_partitions_enforce_the_enclosure_precondition(self):
        for initial in (Fraction(2), Fraction.from_float(1.01)):
            intervals = list(
                metrics_module._mass_intervals(
                    history(segment(q0=1.0, q1=1.0)), initial, Fraction(1), Fraction(1)
                )
            )
            losses = []
            for _, h, _, _, q0, q1, mass in intervals:
                loss = h * (q0 + q1) / 2
                self.assertLessEqual(loss, mass / 2)
                losses.append(loss)
            self.assertEqual(sum(losses, Fraction(0)), Fraction(1))

    def test_direct_bernstein_bounds_cover_exact_rationals_and_log_tails(self):
        zero_flow = (
            Fraction(0),
            Fraction(1),
            Fraction(1),
            Fraction(1),
            Fraction(0),
            Fraction(0),
            Fraction(3),
        )
        lo, hi, ml, mh = metrics_module._bernstein_enclosure(zero_flow)
        self.assertLessEqual(lo, Fraction(1, 3))
        self.assertGreaterEqual(hi, Fraction(1, 3))
        self.assertLessEqual(ml, Fraction(1, 6))
        self.assertGreaterEqual(mh, Fraction(1, 6))
        draining = (
            Fraction(0),
            Fraction(1),
            Fraction(2),
            Fraction(2),
            Fraction(1),
            Fraction(1),
            Fraction(2),
        )
        lo, hi, ml, mh = metrics_module._bernstein_enclosure(draining)
        for precision in (300, 500):
            with localcontext() as context:
                context.prec = precision
                dv = 2 * Decimal(2).ln()
                moment = 2 * (2 * Decimal(2).ln() - 1)

                def convert(x):
                    return Decimal(x.numerator) / Decimal(x.denominator)

                self.assertLessEqual(convert(lo), dv)
                self.assertGreaterEqual(convert(hi), dv)
                self.assertLessEqual(convert(ml), moment)
                self.assertGreaterEqual(convert(mh), moment)

    def test_outward_bounds_preserve_subnormal_positive_residuals(self):
        for value in (Fraction(1, 3), Fraction(1, 10), Fraction(1, 2**2200)):
            upper = metrics_module._upper_output(value, "test bound")
            self.assertGreaterEqual(Fraction(upper), value)
            self.assertLess(Fraction(nextafter(upper, 0.0)), value)
        self.assertEqual(metrics_module._upper_output(Fraction(0), "zero bound"), 0.0)
        with self.assertRaises(DomainError):
            metrics_module._upper_output(Fraction(-1), "negative bound")

    def test_rules_cover_extreme_supported_orders_and_odd_moments(self):
        from math import fsum

        for order in (1, 3, 17, 100, 200):
            nodes, weights = gauss_legendre(order)
            self.assertTrue(all(-1.0 < x < 1.0 for x in nodes))
            self.assertTrue(all(w > 0.0 for w in weights))
            for degree in {0, 1, 2 * order - 2, 2 * order - 1}:
                expected = 0.0 if degree % 2 else 2.0 / (degree + 1)
                self.assertLessEqual(
                    abs(fsum(w * x**degree for x, w in zip(nodes, weights)) - expected), 1e-12
                )

    def test_representable_half_subnormal_increment_is_not_rejected_by_precision_loss(self):
        mass = float.fromhex("0x1p+1023")
        duration = float.fromhex("0x1p-52")
        flow = float.fromhex("0x1p-1022")
        result = free_space_burn_metrics(history(segment(end=duration, f0=1.0, q0=flow)), mass)
        self.assertEqual(result.delta_v_m_s, nextafter(0.0, 1.0))
        self.assertEqual(result.consumed_propellant_kg, nextafter(0.0, 1.0))
        self.assertEqual(result.effective_exhaust_velocity_m_s, float.fromhex("0x1p+1022"))
        self.assertGreater(result.impulsive_equivalent_time_s, 0.0)

    def test_exact_half_range_ties_are_classified_without_float_or_decimal_bias(self):
        tiny = nextafter(0.0, 1.0)
        with self.assertRaises(DomainError):
            free_space_burn_metrics(history(segment(end=0.5, f0=tiny, q0=0.0)), 1.0)
        # F/q is exactly half the smallest subnormal although delta-v rounds
        # to a representable positive number. The required exhaust field fails.
        with self.assertRaises(DomainError) as caught:
            free_space_burn_metrics(history(segment(f0=tiny, q0=2.0)), 3.0)
        self.assertEqual(caught.exception.field, "Effective exhaust velocity")
