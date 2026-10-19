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
