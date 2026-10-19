"""Hostile public-bridge mass ledgers; synthetic data, no hardware validation."""

import unittest
from dataclasses import replace
from decimal import Decimal, localcontext
from itertools import permutations
from math import fsum, nextafter

from rocket_propulsion.core.errors import DomainError
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
from tests import test_burn_bipropellant_numerics as coupled_fixtures


def fixture(flows, gradients, *, scale=1.0, end_flows=None, grouped=True, common=False):
    """Independent linear scale maps and constant/linear stream profiles."""
    end_flows = flows if end_flows is None else end_flows
    names = tuple(f"stream{index}" for index in range(len(flows)))
    artifact = TabulatedBurnArtifact(
        (TabulatedBurnSegment(0.0, 1.0, 0.0, 0.0, 0.0, 0.0, fsum(flows), fsum(end_flows),
                              tuple(zip(names, flows, strict=True)),
                              tuple(zip(names, end_flows, strict=True)), "synthetic"),),
        "synthetic-flow-ledger", "a" * 64)
    axes = (PerformanceAxis("x", "1", (0.0, 1.0)),)
    if common:
        surface = RectilinearPerformanceSurface(
            axes, tuple(PerformanceScalePoint((x,), 1.0, scale + gradients[0] * x)
                        for x in (0.0, 1.0)), "synthetic-common-map", "b" * 64)
        kwargs = {"performance_surface": surface}
    else:
        surface = RectilinearStreamPerformanceSurface(
            axes, tuple(StreamPerformanceScalePoint(
                (x,), 1.0, tuple((name, scale + gradient * x)
                                 for name, gradient in zip(names, gradients, strict=True)))
                        for x in (0.0, 1.0)), "synthetic-stream-map", "c" * 64)
        kwargs = {"stream_performance_surface": surface}
    return PropagatorPropulsionBridge(
        artifact, stream_bindings=tuple(StreamStateBinding(name, "tank" if grouped else name)
                                        for name in names), **kwargs)


def evaluate(model, *, time=0.5, side="right"):
    states = {binding.state_name: 1e300 for binding in model.stream_bindings}
    return model.evaluate(relative_time_s=time, vehicle_mass_kg=1e300, direction=(1, 0, 0),
                          additional_states=states, operating_conditions={"x": 0.0}, side=side)


class BridgeFlowNumericalAttacks(unittest.TestCase):
    def test_mass_profile_preserves_aggregate_interpolation_and_exact_endpoint(self):
        tiny = nextafter(0.0, 1.0)
        value = evaluate(fixture((0.0, 0.0), (0.0, 0.0), end_flows=(tiny, tiny)))
        self.assertEqual(value.mass_derivative_kg_s, -tiny)
        self.assertEqual(dict(value.additional_state_derivatives)["tank"], -tiny)
        end = evaluate(fixture((1e100,), (0.0,), end_flows=(1.0,)), time=1.0, side="left")
        self.assertTrue(end.active)
        self.assertEqual(end.mass_derivative_kg_s, -1.0)
        self.assertEqual(dict(end.additional_state_derivatives)["tank"], -1.0)

    def test_upstream_unrepresentable_nonzero_gradient_is_not_reported_as_insensitivity(self):
        case = coupled_fixtures.BipropellantNumericalAttacks()
        case.setUp()
        with self.assertRaises(DomainError):
            case.evaluate_scale(1e-320)

    def test_grouped_half_subnormal_condition_terms_are_rescued_before_rounding(self):
        tiny = nextafter(0.0, 1.0)
        value = evaluate(fixture((tiny, tiny), (0.5, 0.5)))
        self.assertEqual(dict(value.partials.mass_rate_wrt_conditions)["x"], -tiny)
        partials = dict(dict(value.partials.additional_state_rates_wrt_conditions)["x"])
        self.assertEqual(partials["tank"], -tiny)
        with self.assertRaises(DomainError):
            evaluate(fixture((tiny, tiny), (0.5, 0.5), grouped=False))
        # Exact cancellation is a legitimate zero, even if its terms are too
        # small to export individually; no individual channel is exposed here.
        canceled = evaluate(fixture((tiny, tiny), (0.5, -0.5)))
        self.assertEqual(dict(canceled.partials.mass_rate_wrt_conditions)["x"], 0.0)

    def test_grouped_half_subnormal_drains_and_parameter_rates_are_rescued(self):
        tiny = nextafter(0.0, 1.0)
        value = evaluate(fixture((tiny, tiny), (0.0, 0.0), scale=0.5, common=True))
        self.assertEqual(value.mass_derivative_kg_s, -tiny)
        self.assertEqual(dict(value.additional_state_derivatives)["tank"], -tiny)
        self.assertEqual(value.performance_mass_flow_scale, 0.5)
        self.assertEqual(dict(value.partials.mass_rate_wrt_parameters)["mass_flow_scale"], -tiny)
        grouped = dict(dict(value.partials.additional_state_rates_wrt_parameters)["mass_flow_scale"])
        self.assertEqual(grouped["tank"], -tiny)

    def test_condition_cancellation_is_order_independent_and_preserves_the_small_channel(self):
        # The large channels have opposite map slopes. Only the unit-flow
        # channel changes total flow: d[(A+C)*2 + B*(2+x)]/dx = 1.
        for channels in permutations(((1e16, 1.0), (1.0, 1.0), (1e16, -1.0))):
            with self.subTest(channels=channels):
                flows, gradients = zip(*channels, strict=True)
                value = evaluate(fixture(flows, gradients, scale=2.0))
                self.assertEqual(dict(value.partials.mass_rate_wrt_conditions)["x"], -1.0)
                states = dict(dict(value.partials.additional_state_rates_wrt_conditions)["x"])
                self.assertEqual(states["tank"], -1.0)

    def test_ignition_tangent_preserves_cancellation_in_the_stream_profile(self):
        # A(t)+C(t) is constant; B(t)=1+t, hence a delayed ignition changes
        # each negative drain at +1 kg/s per second of time bias.
        for endpoints in permutations(((1e16, 2e16), (1.0, 2.0), (2e16, 1e16))):
            with self.subTest(endpoints=endpoints):
                start, end = zip(*endpoints, strict=True)
                value = evaluate(fixture(start, (0.0, 0.0, 0.0), end_flows=end))
                self.assertEqual(dict(value.partials.mass_rate_wrt_parameters)["ignition_time_bias_s"], 1.0)
                states = dict(dict(value.partials.additional_state_rates_wrt_parameters)[
                    "ignition_time_bias_s"])
                self.assertEqual(states["tank"], 1.0)

    def test_driver_and_map_scaling_recover_representable_final_products(self):
        for flow, driver, map_scale in ((1e-250, 1e-100, 1e100), (1e250, 1e100, 1e-100)):
            with self.subTest(flow=flow, driver=driver):
                model = fixture((flow,), (map_scale,), scale=map_scale)
                parameters = tuple(replace(parameter, reference_value=driver,
                    minimum=0.0, maximum=max(2.0, 2 * driver))
                    if parameter.name == "mass_flow_scale" else parameter
                    for parameter in model.parameters)
                model = replace(model, parameters=parameters)
                references = []
                for precision in (160, 240):
                    with localcontext() as context:
                        context.prec = precision
                        product = Decimal.from_float(flow) * Decimal.from_float(driver)
                        references.append(float(-product * Decimal.from_float(map_scale)))
                self.assertEqual(references[0], references[1])
                expected = references[1]
                self.assertNotEqual(expected, 0.0)
                value = evaluate(model)
                self.assertEqual(value.mass_derivative_kg_s, expected)
                self.assertEqual(dict(value.additional_state_derivatives)["tank"], expected)
                self.assertEqual(dict(value.partials.mass_rate_wrt_conditions)["x"], expected)
                states = dict(dict(value.partials.additional_state_rates_wrt_conditions)["x"])
                self.assertEqual(states["tank"], expected)
