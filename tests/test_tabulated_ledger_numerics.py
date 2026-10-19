"""Tabulated integrals, moments and producers against independent decimal ledgers."""

import unittest
from dataclasses import replace
from decimal import Decimal, localcontext
from fractions import Fraction
from math import inf, nextafter, ulp

from rocket_propulsion.core.errors import DomainError, InputError
from rocket_propulsion.propulsion.burns import TabulatedBurnArtifact, TabulatedBurnSegment
from rocket_propulsion.propulsion.burns.tabulated import (
    InterpolationPolicy,
    scale_tabulated_artifact,
    tabulated_artifact_from_csv,
    tabulated_artifact_from_profiled_burn,
)


def history(*segments):
    return TabulatedBurnArtifact(tuple(segments), "numerical history", "6" * 64)


def segment(start=0.0, end=1.0, f0=1.0, f1=None, q0=1.0, q1=None):
    return TabulatedBurnSegment(
        start,
        end,
        f0,
        f0 if f1 is None else f1,
        f0,
        f0 if f1 is None else f1,
        q0,
        q0 if q1 is None else q1,
    )


def decimal_ledger(segments, precision):
    with localcontext() as context:
        context.prec = precision
        impulse = mass = moment = Decimal(0)
        for s in segments:
            a, b, f, g, q, r = map(
                Decimal.from_float,
                (
                    s.start_time_s,
                    s.end_time_s,
                    s.start_delivered_thrust_n,
                    s.end_delivered_thrust_n,
                    s.start_total_mass_flow_kg_s,
                    s.end_total_mass_flow_kg_s,
                ),
            )
            h = b - a
            impulse += h * (f + g) / 2
            mass += h * (q + r) / 2
            moment += a * h * (f + g) / 2 + h * h * (f + 2 * g) / 6
        return impulse, mass, moment


class TabulatedLedgerAttacks(unittest.TestCase):
    def close(self, actual, expected, count=3):
        self.assertLessEqual(abs(actual - expected), count * ulp(expected), (actual, expected))

    def test_large_and_small_intermediates_preserve_representable_integrals(self):
        for force, duration in ((1e308, 1e-308), (1e-308, 1e308), (nextafter(0.0, 1.0), 2.0)):
            s = segment(end=duration, f0=force, q0=force)
            reference = decimal_ledger((s,), 1200)
            self.assertEqual(reference, decimal_ledger((s,), 1600))
            for name, expected in (
                ("delivered_impulse_n_s", reference[0]),
                ("axial_impulse_n_s", reference[0]),
                ("propellant_mass_kg", reference[1]),
            ):
                with self.subTest(force=force, field=name):
                    self.close(getattr(s, name), float(expected))

    def test_global_first_moments_use_positive_polynomial_coefficients(self):
        for s in (
            segment(end=1e200, f0=1e-200),
            segment(end=1e-200, f0=1e200),
            segment(start=1e200, end=nextafter(1e200, inf), f0=1e-200, f1=2e-200),
            segment(end=1.0, f0=1e-300, f1=0.0),
        ):
            expected = decimal_ledger((s,), 1200)[2]
            self.assertEqual(float(expected), float(decimal_ledger((s,), 1600)[2]))
            self.close(s.delivered_thrust_first_moment_n_s2, float(expected))

    def test_aggregate_terms_are_not_exported_before_the_sum(self):
        tiny = nextafter(0.0, 1.0)
        artifact = history(
            segment(end=0.25, f0=tiny, q0=tiny),
            segment(start=0.25, end=0.5, f0=tiny, q0=tiny),
            segment(start=0.5, end=1.0, f0=tiny, q0=tiny),
        )
        self.assertEqual(artifact.delivered_total_impulse_n_s, tiny)
        self.assertEqual(artifact.consumed_propellant_kg, tiny)
        self.assertEqual(artifact.thrust_centroid_time_s, 0.5)
        with self.assertRaises(DomainError):
            _ = artifact.segments[0].delivered_impulse_n_s

    def test_centroid_and_isp_do_not_require_representable_moment_or_g0_product(self):
        artifact = history(segment(end=1e200, f0=1e100, q0=1.0))
        self.close(artifact.thrust_centroid_time_s, 5e199)
        with self.assertRaises(DomainError):
            _ = artifact.segments[0].delivered_thrust_first_moment_n_s2
        artifact = history(segment(f0=1e308, q0=1e308))
        with localcontext() as context:
            context.prec = 100
            expected = float(Decimal(1) / Decimal.from_float(9.80665))
        self.close(artifact.system_equivalent_specific_impulse_s, expected)

    def test_grouped_large_plus_small_impulses_preserve_final_rounding(self):
        artifact = history(
            segment(end=1.0, f0=1e16),
            segment(start=1.0, end=2.0, f0=1.0),
            segment(start=2.0, end=3.0, f0=1.0),
        )
        self.assertEqual(artifact.delivered_total_impulse_n_s, 1e16 + 2.0)
        self.assertEqual(
            artifact.delivered_total_impulse_n_s, float(decimal_ledger(artifact.segments, 100)[0])
        )

    def test_query_types_and_small_positive_interpolation_fail_explicitly(self):
        s = segment(f0=0.0, f1=nextafter(0.0, 1.0), q0=0.0, q1=0.0)
        for bad in (True, None, "1", float("nan"), inf, -1.0):
            with self.subTest(time=bad), self.assertRaises(DomainError):
                s.state_at(bad)
        with self.assertRaises(DomainError):
            s.state_at(0.25)
        self.assertEqual(s._axial_at(0.25), Fraction(nextafter(0.0, 1.0)) / 4)
        self.assertEqual(s.state_at(0.0)["delivered_thrust_n"], 0.0)
        self.assertEqual(s.state_at(1.0)["delivered_thrust_n"], nextafter(0.0, 1.0))

    def test_large_endpoint_interpolation_stays_inside_its_convex_hull(self):
        s = segment(end=1e308, f0=1e308, f1=1e-308, q0=0.0, q1=0.0)
        with localcontext() as context:
            context.prec = 1200
            a, b, t = map(
                Decimal.from_float, (s.start_delivered_thrust_n, s.end_delivered_thrust_n, 5e307)
            )
            expected = float(
                a * (1 - t / Decimal.from_float(s.end_time_s))
                + b * t / Decimal.from_float(s.end_time_s)
            )
        self.assertEqual(s.state_at(5e307)["delivered_thrust_n"], expected)

    def test_true_integral_output_overflow_or_erasure_is_typed(self):
        for s in (segment(end=2.0, f0=1e308), segment(end=0.25, f0=nextafter(0.0, 1.0))):
            with self.assertRaises(DomainError) as caught:
                _ = s.delivered_impulse_n_s
            self.assertIsNotNone(caught.exception.field)

    def test_scale_ratio_and_chained_factors_are_not_rounded_early(self):
        base = history(segment(f0=1e-308, q0=1e-308))
        try:
            scaled = scale_tabulated_artifact(
                base,
                {"thrust_scale": 1e308, "specific_impulse_scale": 1e-308},
                source_id="scaled",
                provenance="synthetic",
            )
        except DomainError as error:
            self.fail(f"Known finite scaled endpoints were rejected: {error}")
        with localcontext() as context:
            context.prec = 100
            expected = float(
                Decimal.from_float(1e-308) * Decimal.from_float(1e308) / Decimal.from_float(1e-308)
            )
        self.close(scaled.segments[0].start_total_mass_flow_kg_s, expected)
        base = history(
            replace(
                segment(f0=1e-100, q0=0.0), start_axial_thrust_n=1e-300, end_axial_thrust_n=1e-300
            )
        )
        scaled = scale_tabulated_artifact(
            base,
            {"thrust_scale": 1e-100, "axial_efficiency_scale": 1e200},
            source_id="rescue",
            provenance="synthetic",
        )
        with localcontext() as context:
            context.prec = 100
            expected = float(
                Decimal.from_float(1e-300) * Decimal.from_float(1e-100) * Decimal.from_float(1e200)
            )
        self.assertEqual(scaled.segments[0].start_axial_thrust_n, expected)

    def test_scaling_binds_raw_source_when_derived_base_outputs_do_not_fit(self):
        base = history(segment(end=1e308, f0=1e308, q0=1e308))
        scaled = scale_tabulated_artifact(
            base,
            {"thrust_scale": 1e-308, "time_scale": 1e-308},
            source_id="finite scaled history",
            provenance={"case": 1},
        )
        self.assertGreater(scaled.delivered_total_impulse_n_s, 0.0)
        self.assertEqual(TabulatedBurnArtifact.from_json(scaled.to_json()), scaled)

    def test_boolean_scale_and_real_erased_scaled_outputs_are_rejected(self):
        base = history(segment())
        for scales in (
            {"thrust_scale": True},
            {"time_scale": "1"},
            {"specific_impulse_scale": float("nan")},
            {"thrust_scale": nextafter(0.0, 1.0), "specific_impulse_scale": 1e308},
        ):
            with self.subTest(scales=scales), self.assertRaises(DomainError):
                scale_tabulated_artifact(base, scales, source_id="invalid", provenance=None)

    def test_numerical_policy_is_required_and_bound_to_artifact_identity(self):
        artifact = history(segment())
        payload = artifact.to_dict()
        self.assertIn("integral_numerics", payload)
        for value in (None, "unreviewed"):
            changed = dict(payload)
            if value is None:
                del changed["integral_numerics"]
            else:
                changed["integral_numerics"] = value
            with self.assertRaises(InputError):
                TabulatedBurnArtifact.from_dict(changed)
        from rocket_propulsion.propulsion.burns import PropagatorPropulsionBridge

        bridge = PropagatorPropulsionBridge(artifact, require_complete_stream_binding=False)
        payload = bridge.to_dict()
        self.assertIn("tabulated_state_numerics", payload)
        payload["tabulated_state_numerics"] = "unreviewed"
        with self.assertRaises(InputError):
            PropagatorPropulsionBridge.from_dict(payload)

    def test_ordinary_zero_and_discontinuity_controls_roundtrip(self):
        artifact = history(
            segment(end=1.0, f0=0.0, f1=2.0, q0=0.0, q1=0.2),
            segment(start=1.0, end=2.0, f0=3.0, q0=0.3),
        )
        expected = decimal_ledger(artifact.segments, 100)
        self.close(artifact.delivered_total_impulse_n_s, float(expected[0]))
        self.close(artifact.thrust_centroid_time_s, float(expected[2] / expected[0]))
        self.assertEqual(artifact.state_at(1.0, side="left")["delivered_thrust_n"], 2.0)
        self.assertEqual(artifact.state_at(1.0, side="right")["delivered_thrust_n"], 3.0)
        self.assertEqual(TabulatedBurnArtifact.from_json(artifact.to_json()), artifact)
        zero = history(segment(f0=0.0, q0=0.0))
        self.assertEqual(zero.thrust_centroid_time_s, 0.0)
        self.assertIsNone(zero.system_equivalent_specific_impulse_s)

    def test_hold_policy_cannot_hide_a_real_endpoint_slope(self):
        changing = segment(f0=1.0, f1=nextafter(1.0, inf))
        with self.assertRaises(DomainError):
            TabulatedBurnArtifact((changing,), "hold", "6" * 64, InterpolationPolicy.HOLD)
        constant = history(segment())
        held = replace(constant, interpolation="hold")
        self.assertEqual(held.state_at(0.25)["delivered_thrust_n"], 1.0)

    def test_csv_recognized_hold_spelling_preserves_discontinuities(self):
        source = "time_s,delivered_thrust_n,total_mass_flow_kg_s\n0,1,0.1\n1,2,0.2\n2,0,0\n"
        literal = tabulated_artifact_from_csv(source, source_id="csv", interpolation="hold")
        enum = tabulated_artifact_from_csv(
            source, source_id="csv", interpolation=InterpolationPolicy.HOLD
        )
        self.assertEqual(literal, enum)
        self.assertEqual(literal.state_at(1.0, side="left")["delivered_thrust_n"], 1.0)
        self.assertEqual(literal.state_at(1.0, side="right")["delivered_thrust_n"], 2.0)

    def test_profile_conversion_groups_before_scaling_and_checks_result_authority(self):
        from rocket_propulsion.propulsion.burns import (
            ThrottleSchedule,
            ThrottleSegment,
            simulate_profiled_burn,
        )
        from tests.test_burn_constant_numerics import mono

        definition, point = mono(flow=1e16 + 2.0, force=1e16, mass=2e17, loaded=1e17, dry=0.0)
        streams = tuple(
            replace(point.streams[0], stream_id=f"path-{i}", mass_flow_kg_s=rate)
            for i, rate in enumerate((1e16, 1.0, 1.0))
        )
        point = replace(point, streams=streams)
        result = simulate_profiled_burn(
            definition, point, ThrottleSchedule((ThrottleSegment(1.0, 0.5, 0.5),))
        )
        artifact = tabulated_artifact_from_profiled_burn(result)
        self.assertEqual(
            dict(artifact.segments[0].start_stream_mass_flows_kg_s)["tank"], 5e15 + 1.0
        )
        self.assertEqual(artifact.consumed_propellant_kg, result.summary.consumed_propellant_kg)
        with self.assertRaises(DomainError):
            tabulated_artifact_from_profiled_burn(replace(result, result_hash=""))
        from rocket_propulsion.propulsion.burns.serialization import burn_result_hash

        forged = replace(
            result,
            summary=replace(
                result.summary, ideal_delta_v_m_s=result.summary.ideal_delta_v_m_s + 1.0
            ),
        )
        forged = replace(forged, result_hash=burn_result_hash(forged))
        with self.assertRaises(DomainError):
            tabulated_artifact_from_profiled_burn(forged)
