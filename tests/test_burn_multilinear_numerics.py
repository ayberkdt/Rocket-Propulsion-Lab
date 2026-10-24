"""Synthetic map/bridge attacks with recursive Decimal interpolation oracles."""

import json
import random
import unittest
from decimal import Decimal, localcontext
from itertools import product
from math import fsum, isfinite, ldexp, nextafter, ulp
from sys import float_info

from rocket_propulsion.core.errors import DomainError, InputError
from rocket_propulsion.propulsion.burns import (
    PerformanceAxis,
    PerformanceScalePoint,
    PropagatorPropulsionBridge,
    RectilinearPerformanceSurface,
    RectilinearStreamPerformanceSurface,
    StreamPerformanceScalePoint,
    StreamStateBinding,
    TabulatedBurnArtifact,
    TabulatedBurnSegment,
)

TINY = float.fromhex("0x0.0000000000001p-1022")


def models(axes, data):
    coordinates = tuple(product(*(axis.values for axis in axes)))
    values = {point: data(point) for point in coordinates}
    return (
        RectilinearPerformanceSurface(axes, tuple(PerformanceScalePoint(point, *values[point][:2])
                                      for point in coordinates), "synthetic-total", "a"*64),
        RectilinearStreamPerformanceSurface(axes, tuple(StreamPerformanceScalePoint(
            point, values[point][0], (("fuel", values[point][2]), ("oxidizer", values[point][3])))
            for point in coordinates), "synthetic-stream", "b"*64),
    )


def channels(point):
    return ({"force": point.force_scale, "total": point.mass_flow_scale}
            if isinstance(point, PerformanceScalePoint)
            else {"force": point.force_scale, **dict(point.stream_flow_scales)})


def fields(value):
    if hasattr(value, "mass_flow_scale"):
        return ({"force": value.force_scale, "total": value.mass_flow_scale},
                {"force": dict(value.force_scale_gradient), "total": dict(value.mass_flow_scale_gradient)})
    return ({"force": value.force_scale, **dict(value.stream_flow_scales)},
            {"force": dict(value.force_scale_gradient),
             **{name: dict(gradient) for name, gradient in value.stream_flow_scale_gradients}})


def reference(model, query, precision):
    """Recursive linear blending and central differences of a fixed cell path."""
    with localcontext() as context:
        context.prec = precision
        nodes = {point.coordinates: {name: Decimal.from_float(value)
                 for name, value in channels(point).items()} for point in model.points}
        location = [Decimal.from_float(query[axis.name]) for axis in model.axes]
        cells = []
        for axis, value in zip(model.axes, location, strict=True):
            lower = 0
            for index, coordinate in enumerate(axis.values[:-1]):
                if value >= Decimal.from_float(coordinate):
                    lower = index
            cells.append((axis.values[lower], axis.values[lower+1]))

        def path(coordinates):
            def blend(depth, prefix):
                if depth == len(cells):
                    return nodes[tuple(prefix)]
                left, right = cells[depth]
                t = (coordinates[depth]-Decimal.from_float(left))/(Decimal.from_float(right)-Decimal.from_float(left))
                a, b = blend(depth+1, [*prefix, left]), blend(depth+1, [*prefix, right])
                return {name: a[name]+t*(b[name]-a[name]) for name in a}
            return blend(0, [])

        values = path(location)
        gradients = {name: {} for name in values}
        for column, (axis, (left, right)) in enumerate(zip(model.axes, cells, strict=True)):
            step = (Decimal.from_float(right)-Decimal.from_float(left))*Decimal("1e-40")
            plus, minus = list(location), list(location)
            plus[column] += step
            minus[column] -= step
            upper, lower = path(plus), path(minus)
            for name in values:
                gradients[name][axis.name] = float((upper[name]-lower[name])/(2*step))
        return {name: float(value) for name, value in values.items()}, gradients


def bridge(model):
    streams = (("fuel", 1.), ("oxidizer", 3.))
    artifact = TabulatedBurnArtifact((TabulatedBurnSegment(
        0., 1., 100., 100., 100., 100., 4., 4., streams, streams, "synthetic"),),
        "synthetic-map-bridge", "c"*64)
    return PropagatorPropulsionBridge(artifact,
        stream_bindings=(StreamStateBinding("fuel", "mf"), StreamStateBinding("oxidizer", "mo")),
        **({"performance_surface": model} if isinstance(model, RectilinearPerformanceSurface)
           else {"stream_performance_surface": model}))


class MultilinearNumericalAttacks(unittest.TestCase):
    def assert_scalar(self, actual, expected):
        self.assertTrue(isfinite(actual))
        if expected == 0 or abs(expected) < float_info.min:
            self.assertEqual(actual, expected)
        else:
            self.assertNotEqual(actual, 0.)
            self.assertAlmostEqual(actual, expected, delta=max(8*ulp(expected), abs(expected)*2e-14))

    def assert_map(self, model, query):
        expected = reference(model, query, 1600)
        self.assertEqual(reference(model, query, 1200), expected, "Reference must settle first.")
        value = model.evaluate(query)
        actual, gradients = fields(value)
        self.assertEqual(set(actual), set(expected[0]))
        for name, number in actual.items():
            self.assert_scalar(number, expected[0][name])
            for axis, partial in gradients[name].items():
                self.assert_scalar(partial, expected[1][name][axis])
        return value

    def assert_field(self, field, call):
        with self.assertRaises(DomainError) as caught:
            call()
        self.assertIn(field, caught.exception.field or "")

    def test_overflowed_width_preserves_the_finite_midcell_value_and_gradient(self):
        axes = (PerformanceAxis("x", "1", (-float_info.max, float_info.max)),)
        for model in models(axes, lambda point: (1., 2., 3., 4.) if point[0] < 0 else (3., 4., 5., 6.)):
            with self.subTest(model=type(model).__name__):
                self.assert_map(model, {"x": 0.})

    def test_underflowed_corner_weight_recovers_each_large_channel_contribution(self):
        axes = tuple(PerformanceAxis(name, "1", (0., 1.)) for name in ("x", "y"))
        for model in models(axes, lambda point: tuple(k*(1e300 if all(point) else 1e-300) for k in (1, 2, 3, 4))):
            with self.subTest(model=type(model).__name__):
                self.assert_map(model, {"x": 1e-200, "y": 1e-200})

    def test_constant_large_fields_have_zero_gradients_on_tiny_spans(self):
        axes = (PerformanceAxis("x", "1", (0., 1e-300)),)
        for model in models(axes, lambda point: (1e308, 1e308, 1e308, 1e308)):
            with self.subTest(model=type(model).__name__):
                self.assert_map(model, {"x": 5e-301})

    def test_unrepresentable_cell_fraction_can_still_change_a_subnormal_output(self):
        axes = (PerformanceAxis("x", "1", (0., 1e308)),)
        for model in models(axes, lambda point: (TINY, 2*TINY, 3*TINY, 4*TINY)
                            if point[0] == 0 else (1e308,)*4):
            with self.subTest(model=type(model).__name__):
                self.assert_map(model, {"x": TINY})

    def test_small_gradient_survives_large_term_cancellation_and_axis_order(self):
        for reversed_axes in (False, True):
            names = ("y", "x") if reversed_axes else ("x", "y")
            axes = tuple(PerformanceAxis(name, "1", (0., 1.)) for name in names)

            def data(point, *, reversed_axes=reversed_axes):
                x, y = reversed(point) if reversed_axes else point
                return tuple(k*(1e300 if y == 0 else 1+x) for k in (1, 2, 3, 4))

            for model in models(axes, data):
                with self.subTest(order=names, model=type(model).__name__):
                    value = self.assert_map(model, {"x": .25, "y": .5})
                    self.assertEqual(dict(value.force_scale_gradient)["x"], .5)

    def test_fixed_seed_nd_channel_maps_match_recursive_decimal_paths(self):
        generator = random.Random(9027)
        for dimension in (1, 2, 3):
            for exponent in (-500, 0, 500):
                axes = tuple(PerformanceAxis(f"x{i}", "1", (0., ldexp(1., exponent if i%2 == 0 else -exponent)))
                             for i in range(dimension))
                data = {point: tuple(ldexp(generator.uniform(.5, 4.), -exponent) for _ in range(4))
                        for point in product(*(axis.values for axis in axes))}
                query = {axis.name: axis.values[-1]*generator.uniform(.1, .9) for axis in axes}
                for model in models(axes, data.__getitem__):
                    with self.subTest(dimension=dimension, exponent=exponent, model=type(model).__name__):
                        self.assert_map(model, query)

    def test_grid_line_value_and_right_cell_derivative_remain_distinct(self):
        axes = (PerformanceAxis("x", "1", (0., .25, 1.)),)
        values = {0.: 1., .25: 2., 1.: 4.}
        for model in models(axes, lambda point: (values[point[0]],)*4):
            for query in (nextafter(.25, 0.), .25, nextafter(.25, 1.), 0., 1.):
                with self.subTest(model=type(model).__name__, query=query):
                    self.assert_map(model, {"x": query})

    def test_genuine_force_gradient_overflow_fails_with_an_output_field(self):
        axes = (PerformanceAxis("x", "1", (0., 1e-320)),)
        for model in models(axes, lambda point: (1.,)*4 if point[0] == 0 else (2.,)*4):
            with self.subTest(model=type(model).__name__):
                self.assert_field("Multilinear force gradient", lambda*, model=model: model.evaluate({"x": 5e-321}))

    def test_genuine_underflow_is_detected_per_total_or_named_stream_channel(self):
        axes = (PerformanceAxis("x", "1", (0., 1e308)),)
        for model in models(axes, lambda point: (1.,)*4 if point[0] == 0 else
                            (1., nextafter(1., 2.), nextafter(1., 2.), nextafter(1., 2.))):
            label = "mass-flow gradient" if isinstance(model, RectilinearPerformanceSurface) else "stream gradient: fuel"
            with self.subTest(model=type(model).__name__):
                self.assert_field(label, lambda*, model=model: model.evaluate({"x": 5e307}))

    def test_flat_tiny_maps_keep_genuine_zero_gradients_on_huge_spans(self):
        axes = (PerformanceAxis("x", "1", (-1e308, 1e308)),)
        for model in models(axes, lambda point: (TINY,)*4):
            with self.subTest(model=type(model).__name__):
                self.assert_map(model, {"x": 0.})

    def test_bridge_force_flow_and_condition_partials_follow_the_map_ledger(self):
        axes = tuple(PerformanceAxis(name, "1", (0., 1.)) for name in ("x", "y"))
        for model in models(axes, lambda point: tuple(k*(1e300 if all(point) else 1e-300) for k in (1, 2, 3, 4))):
            query = {"x": 1e-200, "y": 1e-200}
            expected, gradients = reference(model, query, 1600)
            value = bridge(model).evaluate(relative_time_s=.5, vehicle_mass_kg=100., direction=(1., 0., 0.),
                                            additional_states={"mf": 10., "mo": 30.}, operating_conditions=query)
            total = expected.get("total", None)
            fuel, oxidizer = (total, 3*total) if total is not None else (expected["fuel"], 3*expected["oxidizer"])
            with self.subTest(model=type(model).__name__):
                self.assert_scalar(value.axial_thrust_n, 100*expected["force"])
                self.assert_scalar(value.mass_derivative_kg_s, -(fuel+oxidizer))
                self.assert_scalar(dict(value.additional_state_derivatives)["mf"], -fuel)
                self.assert_scalar(dict(value.additional_state_derivatives)["mo"], -oxidizer)
                self.assert_scalar(value.mass_derivative_kg_s, fsum(rate for _, rate in value.additional_state_derivatives))
                for name, partial in value.partials.mass_rate_wrt_conditions:
                    expected_partial = -4*gradients["total"][name] if total is not None else -(gradients["fuel"][name]+3*gradients["oxidizer"][name])
                    self.assert_scalar(partial, expected_partial)

    def test_multilinear_numerical_policy_changes_or_omissions_are_rejected(self):
        axes = (PerformanceAxis("x", "1", (0., 1.)),)
        for model in models(axes, lambda point: (1+point[0],)*4):
            document = json.loads(model.to_json())
            with self.subTest(model=type(model).__name__):
                self.assertIn("interpolation_numerics", document)
                self.assertEqual(type(model).from_dict(document).surface_hash, model.surface_hash)
                for mode in ("replace", "omit"):
                    modified = dict(document)
                    if mode == "replace":
                        modified["interpolation_numerics"] = "unsupported"
                    else:
                        del modified["interpolation_numerics"]
                    with self.subTest(mode=mode), self.assertRaises(InputError):
                        type(model).from_dict(modified)

    def test_exact_binary_affine_coordinate_changes_preserve_values_and_transform_gradients(self):
        axes = (PerformanceAxis("x", "1", (0., 1.)), PerformanceAxis("y", "1", (0., 1.)))
        base = models(axes, lambda p: (1+p[0]+p[1]+p[0]*p[1], 2+p[0], 3+p[1], 4+p[0]*p[1]))
        scales, shifts = (ldexp(1., 500), ldexp(1., -500)), (ldexp(1., 500), ldexp(1., -500))
        changed_axes = tuple(PerformanceAxis(axis.name, "1", (shift, shift+scale))
                             for axis, scale, shift in zip(axes, scales, shifts, strict=True))
        changed = models(changed_axes, lambda p: (1+(p[0]-shifts[0])/scales[0]+(p[1]-shifts[1])/scales[1]+
                          ((p[0]-shifts[0])/scales[0])*((p[1]-shifts[1])/scales[1]),
                          2+(p[0]-shifts[0])/scales[0], 3+(p[1]-shifts[1])/scales[1],
                          4+((p[0]-shifts[0])/scales[0])*((p[1]-shifts[1])/scales[1])))
        query = {"x": .25, "y": .75}
        for ordinary, transformed in zip(base, changed, strict=True):
            a, da = fields(ordinary.evaluate(query))
            b, db = fields(transformed.evaluate({axis.name: shift+scale*query[axis.name]
                           for axis, scale, shift in zip(axes, scales, shifts, strict=True)}))
            for name in a:
                self.assert_scalar(b[name], a[name])
                for axis, scale in zip(axes, scales, strict=True):
                    self.assert_scalar(db[name][axis.name], da[name][axis.name]/scale)

    def test_varying_conditions_integrate_analytic_tank_and_vehicle_consumption(self):
        axes = tuple(PerformanceAxis(name, "1", (0., 1.)) for name in ("x", "y"))
        for model in models(axes, lambda p: (1+p[0]*p[1], 2+p[0]*p[1], 3+2*p[0]*p[1], 4+3*p[0]*p[1])):
            coupling = bridge(model)
            state = (100., 10., 30.)
            step = .02

            def rate(time, values, *, coupling=coupling):
                value = coupling.evaluate(relative_time_s=time, vehicle_mass_kg=values[0], direction=(1., 0., 0.),
                    additional_states={"mf": values[1], "mo": values[2]}, operating_conditions={"x": time, "y": time},
                    side="left" if time == 1 else "right")
                drains = dict(value.additional_state_derivatives)
                return value.mass_derivative_kg_s, drains["mf"], drains["mo"]

            def shifted(values, rates, factor):
                return tuple(v+factor*r for v, r in zip(values, rates, strict=True))

            for i in range(50):
                t = i*step
                a = rate(t, state)
                b = rate(t+step/2, shifted(state, a, step/2))
                c = rate(t+step/2, shifted(state, b, step/2))
                d = rate(t+step, shifted(state, c, step))
                state = tuple(v+step*(aa+2*bb+2*cc+dd)/6 for v, aa, bb, cc, dd in zip(state, a, b, c, d, strict=True))
            expected = (100-28/3, 10-7/3, 23.) if isinstance(model, RectilinearPerformanceSurface) else (100-56/3, 10-11/3, 15.)
            with self.subTest(model=type(model).__name__):
                for actual, target in zip(state, expected, strict=True):
                    self.assert_scalar(actual, target)
                self.assertAlmostEqual(state[0]-state[1]-state[2], 60., delta=2e-11)
