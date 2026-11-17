"""Independent scaled-map and quotient-tangent attacks, using synthetic data.

The pressure-map oracle eliminates the common force/flow scale analytically.
The cancellation oracle differentiates a high-precision linear-fractional path
by symmetric evaluations, not by copying the production quotient-rule formula.
Neither test is an experimental validation of a motor/feed configuration.
"""

import unittest
from dataclasses import replace
from decimal import Decimal, localcontext
from math import inf, ldexp, nextafter
from sys import float_info

from rocket_propulsion.core.errors import DomainError
from tests import test_burn_bipropellant as fixtures


def path_tangent(oxidizer, fuel, d_oxidizer, d_fuel, precision):
    """Differentiate the composite O/F path by Decimal symmetric evaluations."""
    with localcontext() as context:
        context.prec = precision
        ox, fu, dox, dfu = map(Decimal.from_float, (oxidizer, fuel, d_oxidizer, d_fuel))
        step = Decimal("1e-50")
        right = (ox + step * dox) / (fu + step * dfu)
        left = (ox - step * dox) / (fu - step * dfu)
        return float((right - left) / (2 * step))


class BipropellantNumericalAttacks(unittest.TestCase):
    def setUp(self):
        self.case = fixtures.BipropellantDynamicFeedTests()
        self.case.setUp()

    def scaled_coupling(self, scale):
        artifact = self.case.bridge.artifact
        segments = []
        for segment in artifact.segments:
            updates = {name: getattr(segment, name) * scale for name in (
                "start_delivered_thrust_n", "end_delivered_thrust_n",
                "start_axial_thrust_n", "end_axial_thrust_n",
                "start_total_mass_flow_kg_s", "end_total_mass_flow_kg_s")}
            for name in ("start_stream_mass_flows_kg_s", "end_stream_mass_flows_kg_s"):
                updates[name] = tuple((key, value * scale) for key, value in getattr(segment, name))
            segments.append(replace(segment, **updates))
        bridge = replace(self.case.bridge, artifact=replace(artifact, segments=tuple(segments)))
        return replace(self.case.coupling, bridge=bridge)

    def evaluate_scale(self, scale):
        return self.scaled_coupling(scale).evaluate(
            self.case.state, relative_time_s=0.5, vehicle_mass_kg=100.0, direction=(1, 0, 0))

    def evaluate_with_tolerance(self, state, tolerance, *, time=0.5):
        coupling = replace(self.case.coupling, configuration=replace(
            self.case.configuration, mass_closure_absolute_tolerance_kg_s=tolerance,
        ))
        return coupling.evaluate(state, relative_time_s=time, vehicle_mass_kg=100.0,
                                 direction=(1, 0, 0))

    def state_with_flows(self, oxidizer, fuel):
        state = self.case.state
        return replace(
            state, oxidizer_line=replace(state.oxidizer_line, mass_flow_kg_s=oxidizer),
            fuel_line=replace(state.fuel_line, mass_flow_kg_s=fuel),
        )

    def test_component_rounding_cannot_hide_error_above_tolerance(self):
        state = self.case.state
        state = replace(
            state, oxidizer_line=replace(state.oxidizer_line, manifold_pressure_pa=1.3e6),
            fuel_line=replace(state.fuel_line, manifold_pressure_pa=1.5e6),
        )
        result = self.case.evaluate(state)
        with localcontext() as context:
            context.prec = 1100
            expected = abs(sum(map(Decimal.from_float, (
                result.oxidizer_engine_flow_kg_s, result.fuel_engine_flow_kg_s,
                result.propulsion.mass_derivative_kg_s,
            ))))
        self.assertGreater(expected, Decimal("1e-16"))
        self.assertEqual(result.component_to_total_flow_closure_error_kg_s, float(expected))
        with self.assertRaises(DomainError) as error:
            self.evaluate_with_tolerance(state, 1e-16)
        self.assertEqual(error.exception.field, "Bipropellant component closure")

    def test_pressure_storage_rounding_does_not_false_reject_conservative_state(self):
        result = self.evaluate_with_tolerance(self.state_with_flows(0.11, 0.3), 2e-16)
        self.assertTrue(result.propulsion.active)
        self.assertLess(result.vehicle_mass_closure_error_kg_s, 1e-16)

    def test_vehicle_residual_audits_the_integrators_pressure_rates(self):
        with self.assertRaises(DomainError) as error:
            self.evaluate_with_tolerance(self.state_with_flows(0.0, 0.0), 1.5e-16)
        self.assertEqual(error.exception.field, "Bipropellant vehicle closure")

    def test_opposite_leg_errors_cannot_cancel_each_others_acceptance(self):
        state = self.state_with_flows(0.1, -0.1)
        result = self.evaluate_with_tolerance(state, 1e-17, time=3.0)
        self.assertEqual(result.vehicle_mass_closure_error_kg_s, 0.0)
        self.assertGreater(result.oxidizer_line.external_mass_closure_error_kg_s, 1e-18)
        self.assertGreater(result.fuel_line.external_mass_closure_error_kg_s, 1e-18)
        with self.assertRaises(DomainError) as error:
            self.evaluate_with_tolerance(state, 1e-18, time=3.0)
        self.assertEqual(error.exception.field, "Bipropellant oxidizer leg closure")

    def test_either_leg_can_fail_while_the_vehicle_total_is_within_tolerance(self):
        for flows, leg in (((0.25, 0.1), "oxidizer"), ((0.1, 0.25), "fuel")):
            with self.subTest(leg=leg):
                state = self.state_with_flows(*flows)
                valid = self.evaluate_with_tolerance(state, 1e-17, time=3.0)
                self.assertLess(valid.vehicle_mass_closure_error_kg_s, 4e-18)
                with self.assertRaises(DomainError) as error:
                    self.evaluate_with_tolerance(state, 4e-18, time=3.0)
                self.assertEqual(error.exception.field, f"Bipropellant {leg} leg closure")

    def test_exact_closure_boundary_is_inclusive(self):
        coupling = replace(self.case.coupling, oxidizer_line=replace(
            self.case.oxidizer_line, configuration=replace(
                self.case.oxidizer_line.configuration, compliance_kg_pa=1.5,
            ),
        ))
        self.case.coupling = coupling
        state = self.state_with_flows(1.0, 0.0)
        # C * rounded(1/C) − 1 is exactly −2**−54 for C = 1.5.
        boundary = ldexp(1.0, -54)
        result = self.evaluate_with_tolerance(state, boundary, time=3.0)
        self.assertEqual(result.vehicle_mass_closure_error_kg_s, boundary)
        with self.assertRaises(DomainError) as error:
            self.evaluate_with_tolerance(state, nextafter(boundary, 0.0), time=3.0)
        self.assertEqual(error.exception.field, "Bipropellant oxidizer leg closure")

    def test_closure_threshold_uses_unrounded_error(self):
        state = self.state_with_flows(1e-100, 0.3)
        result = self.case.evaluate(state)
        with localcontext() as context:
            context.prec = 1100
            compliance = Decimal(self.case.oxidizer_line.configuration.compliance_kg_pa)
            pressure_rate = Decimal(result.oxidizer_line.derivative.manifold_pressure_pa_s)
            exact = abs(
                compliance * pressure_rate
                + Decimal(result.oxidizer_feed.derivative.propellant_mass_kg_s)
                + Decimal(result.oxidizer_engine_flow_kg_s)
            )
        displayed = float(exact)
        self.assertGreater(exact, Decimal(displayed))
        with self.assertRaises(DomainError) as error:
            self.evaluate_with_tolerance(state, displayed)
        self.assertEqual(error.exception.field, "Bipropellant oxidizer leg closure")
        control = self.evaluate_with_tolerance(state, nextafter(displayed, inf))
        self.assertEqual(control.vehicle_mass_closure_error_kg_s,
                         result.vehicle_mass_closure_error_kg_s)

    def test_active_pressure_jacobian_cannot_erase_a_nonzero_gradient(self):
        scaled = self.scaled_coupling(1e-280)
        nominal = scaled.evaluate(
            self.case.state, relative_time_s=0.5, vehicle_mass_kg=100.0, direction=(1, 0, 0),
        )
        self.assertGreater(abs(nominal.partials.line_state_jacobian[1][1]), 0.0)
        coupling = replace(scaled, oxidizer_line=replace(
            scaled.oxidizer_line, configuration=replace(
                scaled.oxidizer_line.configuration, compliance_kg_pa=1e100,
            ),
        ))
        with self.assertRaises(DomainError) as error:
            coupling.evaluate(self.case.state, relative_time_s=0.5, vehicle_mass_kg=100.0,
                              direction=(1, 0, 0))
        self.assertEqual(error.exception.field, "Oxidizer manifold self-pressure sensitivity")

    def test_representable_subnormal_products_do_not_erase_pressure_gradient(self):
        value = self.evaluate_scale(1e-160)
        derivatives = dict(value.partials.engine_mixture_ratio_wrt_conditions)
        # qO/qF = 3 A(pO)/B(pF); scale has canceled before differentiation.
        expected = 3.0 * 2e-7 / 1.16
        self.assertAlmostEqual(derivatives["oxidizer_pressure_pa"], expected,
                               delta=abs(expected) * 5e-14)

    def test_common_scale_pressure_map_matches_analytic_response(self):
        expected = {"oxidizer_pressure_pa": 3.0 * 2e-7 / 1.16,
                    "fuel_pressure_pa": -3.0 * 1.06 * 1e-7 / 1.16**2}
        for scale in (1.0, 1e-20, 1e-100, 1e-150, 1e-160, 1e-170, 1e-200, 1e-250, 1e-280):
            with self.subTest(scale=scale):
                value = self.evaluate_scale(scale)
                self.assertAlmostEqual(value.engine_mixture_ratio, 3.0 * 1.06 / 1.16,
                                       delta=2e-14)
                actual = dict(value.partials.engine_mixture_ratio_wrt_conditions)
                for name, target in expected.items():
                    self.assertAlmostEqual(actual[name], target, delta=abs(target) * 5e-14)

    def test_factored_tangents_survive_large_small_and_equal_flows(self):
        from rocket_propulsion.propulsion.burns.bipropellant import _mixture_ratio_partial

        for exponent in (-1000, -600, -100, 0, 100, 600, 1000):
            for ratio_exponent in (-5, 0, 5):
                for a, b in ((1.0, nextafter(1.0, inf)), (-1.0, nextafter(-1.0, -inf)),
                             (0.5, 0.5), (-0.5, 0.75)):
                    with self.subTest(exponent=exponent, ratio=ratio_exponent, slopes=(a, b)):
                        # O(t)=2^(k+r)(1+a*t), F(t)=2^k(1+b*t).
                        # All input multiplications are exact binary scaling.
                        actual = _mixture_ratio_partial(
                            ldexp(1.0, exponent + ratio_exponent), ldexp(1.0, exponent),
                            ldexp(a, exponent + ratio_exponent), ldexp(b, exponent))
                        self.assertEqual(actual, ldexp(a - b, ratio_exponent))

    def test_nearly_correlated_tangents_match_independent_composite_path(self):
        from rocket_propulsion.propulsion.burns.bipropellant import _mixture_ratio_partial

        for oxidizer, fuel in ((3.0, 7.0), (1.1, 1.3), (17.0, 23.0), (0.3, 0.7)):
            for exponent in (-800, 0, 800):
                values = tuple(ldexp(value, exponent) for value in (
                    oxidizer, fuel, oxidizer, nextafter(fuel, inf)))
                with self.subTest(flows=(oxidizer, fuel), exponent=exponent):
                    first = path_tangent(*values, 160)
                    refined = path_tangent(*values, 240)
                    self.assertEqual(first, refined, "The independent oracle must settle first.")
                    self.assertNotEqual(refined, 0.0)
                    self.assertEqual(_mixture_ratio_partial(*values), refined)

    def test_unrepresentable_outputs_fail_typed_and_exact_zero_remains_valid(self):
        from rocket_propulsion.propulsion.burns.bipropellant import (
            _finite_quotient,
            _mixture_ratio_partial,
        )

        tiny = nextafter(0.0, 1.0)
        huge = float_info.max
        for numerator, denominator in ((huge, tiny), (tiny, huge), (1.0, 0.0)):
            with self.subTest(ratio=(numerator, denominator)), self.assertRaises(DomainError):
                _finite_quotient(numerator, denominator, "Test ratio")
        for values in ((1.0, tiny, 1.0, 0.0), (1.0, huge, tiny, 0.0)):
            with self.subTest(tangent=values), self.assertRaises(DomainError):
                _mixture_ratio_partial(*values)
        self.assertEqual(_mixture_ratio_partial(huge, huge, tiny, tiny), 0.0)
        self.assertEqual(_finite_quotient(0.0, tiny, "Exact zero"), 0.0)
        # Exercise the public line-ratio output too, not only private helpers.
        state = replace(self.case.state, fuel_line=replace(
            self.case.state.fuel_line, mass_flow_kg_s=tiny))
        with self.assertRaises(DomainError):
            self.case.evaluate(state)

    def test_nonfinite_composed_results_cannot_reach_a_propagator(self):
        value = self.case.evaluate()
        for field, replacement in (
            ("line_state_jacobian", ((inf, 0.0, 0.0, 0.0),) + value.partials.line_state_jacobian[1:]),
            ("line_derivative_wrt_feed_state", (("tank", (0.0, inf, 0.0, 0.0)),)),
            ("engine_component_flow_wrt_conditions", (("oxidizer", (("pressure", inf),)),)),
            ("engine_mixture_ratio_wrt_conditions", (("pressure", inf),)),
            ("acceleration_wrt_line_state", (("pressure", (0.0, inf, 0.0)),)),
        ):
            with self.subTest(partial=field), self.assertRaises(DomainError):
                replace(value.partials, **{field: replacement})
        for field in ("engine_mixture_ratio", "line_mixture_ratio",
                      "stored_liquid_mass_rate_kg_s", "vehicle_mass_closure_error_kg_s"):
            with self.subTest(output=field), self.assertRaises(DomainError):
                replace(value, **{field: inf})
        for field, malformed in (
            ("line_state_jacobian", value.partials.line_state_jacobian[:3]),
            ("line_state_order", ("q", "p", "q", "p")),
            ("line_derivative_wrt_feed_state", (("tank", (0.0, 0.0, 0.0)),)),
            ("engine_mixture_ratio_wrt_conditions", (("p", 1.0), ("p", 2.0))),
            ("acceleration_wrt_line_state", (("p", (0.0, 0.0)),)),
        ):
            with self.subTest(shape=field), self.assertRaises(DomainError):
                replace(value.partials, **{field: malformed})
        matrix = [list(row) for row in value.partials.line_state_jacobian]
        saved = replace(value.partials, line_state_jacobian=matrix)
        matrix[0][0] += 1.0
        self.assertEqual(saved.line_state_jacobian, value.partials.line_state_jacobian)
