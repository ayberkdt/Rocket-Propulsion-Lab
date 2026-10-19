"""Direction-scale and tangent attacks with independent Decimal force paths.

All fixtures are synthetic. No pointing hardware or measured motor performance
is validated. Reference Jacobians differentiate the normalized force path by
two-sided Decimal evaluations rather than copying the implementation projector.
"""

import unittest
from dataclasses import replace
from decimal import Decimal, localcontext
from math import isfinite, ldexp, ulp
from sys import float_info

from rocket_propulsion.core.errors import DomainError
from rocket_propulsion.propulsion.burns import (
    PropagatorPropulsionBridge,
    StreamStateBinding,
    TabulatedBurnArtifact,
    TabulatedBurnSegment,
    ThrottleableEngineModel,
)
from tests import test_burn_electric as electric_fixtures


def bridge(thrust=150.0, area=None):
    segment = TabulatedBurnSegment(0, 1, thrust, thrust, thrust, thrust, 1, 1,
                                   (("propellant", 1),), (("propellant", 1),), "synthetic")
    return PropagatorPropulsionBridge(
        TabulatedBurnArtifact((segment,), "synthetic-direction", "a" * 64),
        stream_bindings=(StreamStateBinding("propellant", "tank"),), nozzle_exit_area_m2=area)


def engines():
    return {
        "bridge": bridge(),
        "throttle": ThrottleableEngineModel(150.0, 1.0, 0.1, 1.0, (("tank", 1.0),)),
        "electric": electric_fixtures.thrust_model(propellant_state_name="tank"),
    }


def evaluate(name, model, direction, *, mass=100.0, ambient=None, off=False):
    arguments = {"vehicle_mass_kg": mass, "direction": direction, "additional_states": {"tank": 1000.0}}
    if name == "bridge":
        return model.evaluate(relative_time_s=2.0 if off else 0.5,
                              ambient_pressure_pa=ambient, **arguments)
    if name == "throttle":
        return model.evaluate(throttle=0.0 if off else 1.0, ambient_pressure_pa=ambient, **arguments)
    return model.evaluate(level_index=-1 if off else 0, available_power_w=1000.0, **arguments)


def decimal_geometry(direction, thrust, mass, precision):
    with localcontext() as context:
        context.prec = precision
        d = tuple(Decimal.from_float(float(value)) for value in direction)
        force, vehicle = Decimal.from_float(float(thrust)), Decimal.from_float(float(mass))

        def acceleration(coordinates):
            radius = sum(value * value for value in coordinates).sqrt()
            return tuple(force * value / (vehicle * radius) for value in coordinates)

        value = acceleration(d)
        step = max(abs(value) for value in d) * Decimal("1e-40")
        columns = []
        for column in range(3):
            plus, minus = list(d), list(d)
            plus[column] += step
            minus[column] -= step
            right, left = acceleration(plus), acceleration(minus)
            columns.append(tuple((right_value - left_value) / (2 * step)
                                 for right_value, left_value in zip(right, left, strict=True)))
        matrix = tuple(tuple(float(columns[column][row]) for column in range(3)) for row in range(3))
        return tuple(float(component) for component in value), matrix


class DirectionNumericalAttacks(unittest.TestCase):
    def assert_scalar_close(self, actual, expected):
        self.assertTrue(isfinite(actual))
        if expected == 0.0:
            self.assertEqual(actual, 0.0)
        else:
            self.assertAlmostEqual(actual, expected, delta=max(8 * ulp(expected), abs(expected) * 2e-14))

    def assert_geometry(self, value, direction, thrust, mass=100.0, precisions=(160, 240)):
        first = decimal_geometry(direction, thrust, mass, precisions[0])
        reference = decimal_geometry(direction, thrust, mass, precisions[1])
        self.assertEqual(first, reference, "Independent force-path reference must settle first.")
        for actual, expected in zip(value.acceleration_m_s2, reference[0], strict=True):
            self.assert_scalar_close(actual, expected)
        for row, expected_row in zip(value.partials.acceleration_wrt_direction, reference[1], strict=True):
            for actual, expected in zip(row, expected_row, strict=True):
                self.assert_scalar_close(actual, expected)

    def test_huge_axis_scaling_preserves_force_in_every_shared_consumer(self):
        for name, model in engines().items():
            with self.subTest(model=name):
                baseline = evaluate(name, model, (1.0, 0.0, 0.0))
                scaled = evaluate(name, model, (1e308, 0.0, 0.0))
                self.assertEqual(scaled.force_vector_n, baseline.force_vector_n)
                self.assertEqual(scaled.acceleration_m_s2, baseline.acceleration_m_s2)

    def test_small_axis_scaling_with_representable_jacobian_is_accepted(self):
        for name, model in engines().items():
            with self.subTest(model=name):
                baseline = evaluate(name, model, (1.0, 0.0, 0.0))
                value = evaluate(name, model, (1e-300, 0.0, 0.0))
                self.assertEqual(value.force_vector_n, baseline.force_vector_n)
                self.assert_geometry(value, (1e-300, 0.0, 0.0), baseline.force_vector_n[0])

    def test_unrepresentable_length_does_not_discard_representable_vectors(self):
        direction = (float_info.max,) * 3
        for name, model in engines().items():
            with self.subTest(model=name):
                baseline = evaluate(name, model, (1.0, 0.0, 0.0))
                value = evaluate(name, model, direction)
                self.assert_geometry(value, direction, baseline.force_vector_n[0])

    def test_power_of_two_rescaling_preserves_force_and_inverse_tangent_scaling(self):
        for name, model in engines().items():
            nominal = evaluate(name, model, (1.0, 2.0, 2.0))
            force = evaluate(name, model, (1.0, 0.0, 0.0)).force_vector_n[0]
            for exponent in (-900, -200, 200, 900):
                direction = tuple(ldexp(component, exponent) for component in (1.0, 2.0, 2.0))
                with self.subTest(model=name, exponent=exponent):
                    value = evaluate(name, model, direction)
                    for actual, expected in zip(value.force_vector_n, nominal.force_vector_n, strict=True):
                        self.assert_scalar_close(actual, expected)
                    self.assert_geometry(value, direction, force)

    def test_near_axis_projector_keeps_representable_diagonal_sensitivity(self):
        direction = (1.0, 1e-300, 0.0)
        value = evaluate("bridge", bridge(1e300), direction)
        self.assertNotEqual(value.partials.acceleration_wrt_direction[0][0], 0.0)
        self.assert_geometry(value, direction, 1e300, precisions=(800, 1000))

    def test_signed_thrust_is_preserved_in_direction_tangents(self):
        models = engines()
        for name in ("bridge", "throttle"):
            model = replace(models[name], nozzle_exit_area_m2=1.0)
            value = evaluate(name, model, (1.0, 2.0, 2.0), ambient=200.0)
            self.assertLess(value.net_thrust_margin_n, 0.0)
            self.assert_geometry(value, (1.0, 2.0, 2.0), -50.0)

    def test_off_and_zero_force_keep_valid_tiny_directions(self):
        tiny = float.fromhex("0x0.0000000000001p-1022")
        for name, model in engines().items():
            with self.subTest(model=name):
                value = evaluate(name, model, (tiny, 0.0, 0.0), off=True)
                self.assertEqual(value.force_vector_n, (0.0, 0.0, 0.0))
        value = evaluate("bridge", bridge(0.0), (tiny, 0.0, 0.0))
        self.assertEqual(value.partials.acceleration_wrt_direction, ((0.0, 0.0, 0.0),) * 3)

    def test_projection_can_recover_an_unrepresentable_unit_component(self):
        from rocket_propulsion.propulsion.burns.direction import _direction_geometry

        geometry = _direction_geometry((1e308, 1e-308, 0.0))
        self.assertEqual(geometry.project(1e308, label="Force"), (1e308, 1e-308, 0.0))
        self.assertEqual(geometry.project(0.0, label="Exact zero"), (0.0, 0.0, 0.0))

    def test_true_unrepresentable_jacobian_has_a_typed_output_error(self):
        tiny = float.fromhex("0x0.0000000000001p-1022")
        with self.assertRaises(DomainError) as caught:
            evaluate("bridge", bridge(), (tiny, 0.0, 0.0))
        self.assertIn("Direction Jacobian", caught.exception.field or "")
        # A lower force makes every requested output representable again.
        value = evaluate("bridge", bridge(1e-300), (tiny, 0.0, 0.0))
        self.assertEqual(value.force_vector_n, (1e-300, 0.0, 0.0))
        self.assert_geometry(value, (tiny, 0.0, 0.0), 1e-300, precisions=(700, 900))
