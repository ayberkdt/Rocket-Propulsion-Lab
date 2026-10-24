"""Synthetic PCHIP/tensor attacks; Decimal Bezier references, no hardware claim."""

import json
import random
import unittest
from bisect import bisect_right
from decimal import Decimal, localcontext
from itertools import pairwise, product
from math import isfinite, ldexp, nextafter, prod, ulp
from sys import float_info

from rocket_propulsion.core.errors import DomainError, InputError
from rocket_propulsion.propulsion.burns import (
    MonotoneCubicCurve,
    PerformanceAxis,
    PerformanceScalePoint,
    RectilinearPerformanceSurface,
    ThrottleableEngineModel,
)
from rocket_propulsion.propulsion.burns.interpolation import _hermite_segment, pchip_slopes
from tests.test_burn_pulse_mode import synthetic_thruster


def decimal_reference(x, y, query, precision):
    """PCHIP specification slopes and independent Bezier/de Casteljau evaluation."""
    with localcontext() as context:
        context.prec = precision
        xx = [Decimal.from_float(float(value)) for value in x]
        yy = [Decimal.from_float(float(value)) for value in y]
        q = Decimal.from_float(float(query))
        widths = [b - a for a, b in pairwise(xx)]
        secants = [(b - a) / h for a, b, h in zip(yy, yy[1:], widths)]

        def endpoint(h0, h1, d0, d1):
            slope = ((2*h0 + h1)*d0 - h0*d1) / (h0 + h1)
            if d0 == 0 or slope == 0 or (slope > 0) != (d0 > 0):
                return Decimal(0)
            if (d0 > 0) != (d1 > 0) and abs(slope) > 3*abs(d0):
                return 3*d0
            return slope

        if len(xx) == 2:
            slopes = [secants[0], secants[0]]
        else:
            slopes = [endpoint(widths[0], widths[1], secants[0], secants[1])]
            for i in range(1, len(xx)-1):
                a, b = secants[i-1:i+1]
                if a == 0 or b == 0 or (a > 0) != (b > 0):
                    slopes.append(Decimal(0))
                else:
                    w1, w2 = 2*widths[i]+widths[i-1], widths[i]+2*widths[i-1]
                    slopes.append((w1+w2)/(w1/a+w2/b))
            slopes.append(endpoint(widths[-1], widths[-2], secants[-1], secants[-2]))
        index = min(max(bisect_right(xx, q)-1, 0), len(xx)-2)
        h = widths[index]
        t = (q-xx[index])/h
        controls = [yy[index], yy[index]+h*slopes[index]/3,
                    yy[index+1]-h*slopes[index+1]/3, yy[index+1]]

        def de_casteljau(values):
            while len(values) > 1:
                values = [(1-t)*a+t*b for a, b in pairwise(values)]
            return values[0]

        value = de_casteljau(controls)
        derivative = de_casteljau([3*(b-a)/h for a, b in pairwise(controls)])
        return tuple(float(v) for v in slopes), float(value), float(derivative)


def surface(x, y):
    return RectilinearPerformanceSurface(
        (PerformanceAxis("x", "1", x),),
        tuple(PerformanceScalePoint((coordinate,), value, 2*value)
              for coordinate, value in zip(x, y, strict=True)),
        "synthetic-pchip", "a"*64, interpolation="cubic_hermite")


class PchipNumericalAttacks(unittest.TestCase):
    def assert_scalar(self, actual, expected):
        self.assertTrue(isfinite(actual))
        if expected == 0:
            self.assertEqual(actual, 0.)
        else:
            self.assertNotEqual(actual, 0.)
            self.assertAlmostEqual(actual, expected, delta=max(8*ulp(expected), abs(expected)*2e-14))

    def reference(self, x, y, query):
        # Extreme binary exponents can put exact endpoint slopes on a float
        # rounding midpoint. Preserve enough decimal digits to settle those
        # sampled dyadic cancellations before comparing production outputs.
        first = decimal_reference(x, y, query, 1200)
        second = decimal_reference(x, y, query, 1600)
        self.assertEqual(first, second, "Reference must settle independently of the implementation.")
        return second

    def assert_curve(self, x, y, query):
        reference = self.reference(x, y, query)
        value, derivative = MonotoneCubicCurve(x, y).evaluate(query)
        self.assert_scalar(value, reference[1])
        self.assert_scalar(derivative, reference[2])
        return value, derivative

    def assert_field(self, field, call):
        with self.assertRaises(DomainError) as caught:
            call()
        self.assertIn(field, caught.exception.field or "")

    def test_small_positive_secants_are_not_flattened(self):
        x, y = (.25, .5, .75), (1e-200, 2e-200, 3e-200)
        expected = self.reference(x, y, .3125)
        for actual, target in zip(pchip_slopes(x, y), expected[0], strict=True):
            self.assert_scalar(actual, target)
        self.assert_curve(x, y, .3125)

    def test_unreduced_compiled_ratios_preserve_signed_and_extreme_curve_outputs(self):
        for x, y, query in (
            ((.125,.875,1.125),(1e200,-1e200,3e199),.5),
            ((0.,1.,2.),(-1e-200,1e-200,-3e-200),.5),
            ((.01,.02,.05),(.001,.004,.04),.015),
            ((0.,1e308),(0.,1e308),1.),
        ):
            curve=MonotoneCubicCurve(x,y)
            index=min(max(bisect_right(x,query)-1,0),len(x)-2)
            n,d,denominator,*coefficients=curve._segments[index]
            expected=curve.evaluate(query)
            independent=self.reference(x,y,query)
            self.assert_scalar(expected[0],independent[1])
            self.assert_scalar(expected[1],independent[2])
            for exponent in (0,512,4096):
                with self.subTest(x=x,y=y,exponent=exponent):
                    factor=2**exponent
                    packed=(n,d,denominator*factor,*(value*factor for value in coefficients))
                    self.assertEqual(_hermite_segment(packed,query),expected)
                    self.assertEqual(_hermite_segment(packed,query,derivative=False),expected[0])

    def test_small_reversals_and_endpoint_limiting_follow_shape_rules(self):
        x = (0., .5, 1.)
        for y in ((0., 1e-200, -99e-200), (3e-200, 2e-200, 1.5e-200),
                  (1e-200, 2e-200, 1.5e-200), (1., 1., 2.)):
            with self.subTest(y=y):
                expected = self.reference(x, y, .125)
                for actual, target in zip(pchip_slopes(x, y), expected[0], strict=True):
                    self.assert_scalar(actual, target)
                self.assert_curve(x, y, .125)

    def test_overflowed_ordinate_difference_has_representable_slopes_and_samples(self):
        self.assert_curve((0., 2., 4.), (-1e308, 1e308, 1.5e308), 1.)

    def test_overflowed_coordinate_width_keeps_finite_value_and_gradient(self):
        self.assert_curve((-float_info.max, float_info.max), (1., 3.), 0.)

    def test_fraction_rounding_does_not_corrupt_a_small_value_near_a_large_cancellation(self):
        self.assert_curve((0., 3.), (-1e308, 1e308), nextafter(1.5, float("inf")))

    def test_constant_large_values_have_exact_zero_derivative_on_tiny_span(self):
        value, derivative = MonotoneCubicCurve((0., 1e-300), (1e308, 1e308)).evaluate(5e-301)
        self.assertEqual(value, 1e308)
        self.assertEqual(derivative, 0.)

    def test_subnormal_coordinates_and_ordinates_preserve_unit_derivative(self):
        tiny = float.fromhex("0x0.0000000000001p-1022")
        value, derivative = MonotoneCubicCurve((-tiny, tiny), (-tiny, tiny)).evaluate(0.)
        self.assertEqual(value, 0.)
        self.assertEqual(derivative, 1.)

    def test_fixed_seed_curves_are_covariant_under_extreme_binary_scaling(self):
        generator = random.Random(611)
        for sample in range(12):
            x = (0., .5, 1.5, 3.)
            y = tuple(generator.uniform(.25, 3.) for _ in x)
            for x_scale, y_scale in ((-500, -500), (500, 500), (-500, 500), (500, -500)):
                scaled_x = tuple(ldexp(v, x_scale) for v in x)
                scaled_y = tuple(ldexp(v, y_scale) for v in y)
                with self.subTest(sample=sample, x_scale=x_scale, y_scale=y_scale):
                    self.assert_curve(scaled_x, scaled_y, ldexp(.875, x_scale))

    def test_public_slope_range_failure_does_not_destroy_later_plateau_evaluation(self):
        x, y = (0., 1e-320, 1.), (0., 1., 1.)
        self.assert_field("PCHIP nodal slope", lambda: pchip_slopes(x, y))
        curve = MonotoneCubicCurve(x, y)
        self.assertEqual(curve.evaluate(1.), (1., 0.))
        self.assertEqual(curve._value_checked(0.), 0.)
        self.assert_field("Curve derivative", lambda: curve.evaluate(0.))

    def test_genuinely_unrepresentable_curve_outputs_fail_typed(self):
        tiny = float.fromhex("0x0.0000000000001p-1022")
        self.assert_field("Curve value", lambda: MonotoneCubicCurve((0., 2.), (0., tiny)).evaluate(1.))
        self.assert_field("Curve derivative", lambda: MonotoneCubicCurve(
            (0., 1e308), (0., tiny)).evaluate(1e308))
        self.assertEqual(MonotoneCubicCurve((0., 1e308), (0., 0.)).evaluate(1e308), (0., 0.))

    def test_direct_kernel_rejects_malformed_and_nonfinite_inputs(self):
        for x, y in (((), ()), ((0.,), (1.,)), ((0., 1.), (1.,)),
                     ((0., 0.), (1., 2.)), ((1., 0.), (1., 2.)),
                     ((0., float("inf")), (1., 2.)), ((0., 1.), (True, 2.)),
                     ((0., 1.), (1., float("nan")))):
            with self.subTest(x=x, y=y), self.assertRaises(DomainError):
                pchip_slopes(x, y)

    def test_throttle_consumer_preserves_the_scaled_curve_flow_and_tangent(self):
        curve = MonotoneCubicCurve((.25, .5, .75), (1e-200, 2e-200, 3e-200))
        model = ThrottleableEngineModel(150., 1e200, .25, .75, (("tank", 1.),), flow_curve=curve)
        value = model.evaluate(throttle=.3125, vehicle_mass_kg=100., direction=(1., 0., 0.),
                               additional_states={"tank": 10.})
        self.assert_scalar(value.mass_derivative_kg_s, -1.25)
        self.assert_scalar(value.partials.mass_rate_wrt_throttle, -4.)

    def test_one_dimensional_cubic_surface_matches_independent_bezier_reference(self):
        x, y, query = (.25, .5, .75), (1e-200, 2e-200, 3e-200), .3125
        reference = self.reference(x, y, query)
        value = surface(x, y).evaluate({"x": query})
        self.assert_scalar(value.force_scale, reference[1])
        self.assert_scalar(dict(value.force_scale_gradient)["x"], reference[2])
        self.assert_scalar(value.mass_flow_scale, 2*reference[1])
        self.assert_scalar(dict(value.mass_flow_scale_gradient)["x"], 2*reference[2])

    def test_cubic_surface_handles_width_beyond_the_float_range(self):
        x = (-float_info.max, float_info.max)
        value = surface(x, (1., 3.)).evaluate({"x": 0.})
        reference = self.reference(x, (1., 3.), 0.)
        self.assert_scalar(value.force_scale, 2.)
        self.assert_scalar(dict(value.force_scale_gradient)["x"], reference[2])

    def test_genuine_cubic_gradient_range_errors_are_typed_and_flat_maps_remain_valid(self):
        configurations = (
            (1e308, (1., nextafter(1., 2.)), (2., nextafter(2., 3.)), "Cubic force gradient"),
            (1e-320, (1., 2.), (2., 4.), "Cubic force gradient"),
            (1e308, (1., 1.), (1., nextafter(1., 2.)), "Cubic mass-flow gradient"),
        )
        for span, force, flow, field in configurations:
            model = RectilinearPerformanceSurface(
                (PerformanceAxis("x", "1", (0., span)),),
                tuple(PerformanceScalePoint((coordinate,), a, b)
                      for coordinate, a, b in zip((0., span), force, flow, strict=True)),
                "synthetic-range", "d"*64, interpolation="cubic_hermite")
            with self.subTest(field=field, span=span):
                self.assert_field(field, lambda*, model=model, span=span: model.evaluate({"x": span/2}))
        value = surface((0., 1e308), (1e-300, 1e-300)).evaluate({"x": 5e307})
        self.assertEqual(value.force_scale_gradient, (("x", 0.),))
        self.assertEqual(value.mass_flow_scale_gradient, (("x", 0.),))

    def test_unrepresentable_mixed_partials_do_not_erase_finite_tensor_outputs(self):
        for dimension in (2, 3):
            with self.subTest(dimension=dimension):
                h = 1e-180
                axes = tuple(PerformanceAxis(f"x{i}", "1", (0., h, 2*h)) for i in range(dimension))
                points = tuple(PerformanceScalePoint(tuple(i*h for i in index),
                               float(prod(i+1 for i in index)), float(2*prod(i+1 for i in index)))
                               for index in product(range(3), repeat=dimension))
                model = RectilinearPerformanceSurface(axes, points, "synthetic-tensor", "b"*64,
                                                      interpolation="cubic_hermite")
                value = model.evaluate({axis.name: h/4 for axis in axes})
                self.assert_scalar(value.force_scale, 1.25**dimension)
                self.assert_scalar(value.mass_flow_scale, 2*1.25**dimension)
                for _, gradient in value.force_scale_gradient:
                    self.assert_scalar(gradient, 1.25**(dimension-1)/h)

    def test_anisotropic_nonlinear_tensor_matches_separable_bezier_paths(self):
        x = tuple(ldexp(v, 500) for v in (0., 1., 2.))
        y = tuple(ldexp(v, -500) for v in (0., 1., 2.))
        a, b = (1., 2., 2.5), (1., 1.5, 3.)
        axes = (PerformanceAxis("x", "1", x), PerformanceAxis("y", "1", y))
        points = tuple(PerformanceScalePoint((xx, yy), a[i]*b[j], 2*a[i]*b[j])
                       for i, xx in enumerate(x) for j, yy in enumerate(y))
        model = RectilinearPerformanceSurface(axes, points, "synthetic-separable", "c"*64,
                                              interpolation="cubic_hermite")
        qx, qy = ldexp(.25, 500), ldexp(.75, -500)
        first, second = self.reference(x, a, qx), self.reference(y, b, qy)
        value = model.evaluate({"x": qx, "y": qy})
        self.assert_scalar(value.force_scale, first[1]*second[1])
        gradients = dict(value.force_scale_gradient)
        self.assert_scalar(gradients["x"], first[2]*second[1])
        self.assert_scalar(gradients["y"], first[1]*second[2])

    def test_cubic_surface_numerical_policy_is_bound_and_verified(self):
        model = surface((0., 1., 2.), (1., 2., 3.))
        document = json.loads(model.to_json())
        self.assertIn("interpolation_numerics", document)
        self.assertEqual(RectilinearPerformanceSurface.from_dict(document).surface_hash, model.surface_hash)
        for change in ("replace", "omit"):
            modified = dict(document)
            if change == "replace":
                modified["interpolation_numerics"] = "unsupported"
            else:
                del modified["interpolation_numerics"]
            with self.subTest(change=change), self.assertRaises(InputError):
                RectilinearPerformanceSurface.from_dict(modified)

    def test_pulse_consumer_curve_policy_is_bound_and_verified(self):
        model = synthetic_thruster()
        document = json.loads(model.to_json())
        self.assertIn("curve_numerics", document)
        self.assertEqual(type(model).from_dict(document).model_hash, model.model_hash)
        for change in ("replace", "omit"):
            modified = dict(document)
            if change == "replace":
                modified["curve_numerics"] = "unsupported"
            else:
                del modified["curve_numerics"]
            with self.subTest(change=change), self.assertRaises(InputError):
                type(model).from_dict(modified)
