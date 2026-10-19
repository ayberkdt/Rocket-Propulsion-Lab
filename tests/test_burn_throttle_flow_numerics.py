"""Synthetic throttle-flow attacks with independent Decimal path/ledger oracles.

No measured motor, feed hardware or physical uncertainty is validated here.
"""

import unittest
import json
from decimal import Decimal, localcontext
from itertools import permutations
from math import fsum, isfinite, log, ulp

from rocket_propulsion.core.errors import DomainError, InputError
from rocket_propulsion.propulsion.burns import MonotoneCubicCurve, ThrottleableEngineModel
from rocket_propulsion.propulsion.performance import STANDARD_GRAVITY_M_S2
from tests.test_burn_direction_numerics import engines as direction_consumers


TINY = float.fromhex("0x0.0000000000001p-1022")


def engine(flow=1.0, thrust=150.0, fractions=(("tank", 1.0),), curve=None):
    return ThrottleableEngineModel(thrust, flow, .25 if curve else .1,
                                  .75 if curve else 1.0, fractions, flow_curve=curve)


def evaluate(model, throttle=1.0, scale=1.0):
    return model.evaluate(throttle=throttle, vehicle_mass_kg=100.0, direction=(1., 0., 0.),
                          additional_states={name: 10.0 for name, _ in model.stream_fractions},
                          parameter_overrides={"mass_flow_scale": scale})


def decimal_ledger(model, throttle, scale, precision):
    """Differentiate a separate flow path; linear data need no cubic evaluator."""
    with localcontext() as context:
        context.prec = precision
        def number(value):
            return Decimal.from_float(float(value))
        rated, command, driver = map(number, (model.rated_mass_flow_kg_s, throttle, scale))
        weights = {name: number(value) for name, value in model.stream_fractions}
        total_weight = sum(weights.values())

        def path(u, s):
            if model.flow_curve is None:
                fraction = u
            else:
                x0, x1 = map(number, model.flow_curve.x)
                y0, y1 = map(number, model.flow_curve.y)
                fraction = y0 + (u - x0) * (y1 - y0) / (x1 - x0)
            return -rated * fraction * s

        step = Decimal("1e-30")
        rate = path(command, driver)
        command_partial = (path(command + step, driver) - path(command - step, driver)) / (2 * step)
        driver_partial = (path(command, driver + step) - path(command, driver - step)) / (2 * step)
        vacuum = number(model.rated_thrust_n) * command
        result = {"total": float(rate), "command": float(command_partial),
                  "driver": float(driver_partial),
                  "isp": float(vacuum / (-number(STANDARD_GRAVITY_M_S2) * rate))}
        for label, value in (("tank", rate), ("tank_command", command_partial),
                             ("tank_driver", driver_partial)):
            result[label] = {name: float(value * weight / total_weight)
                             for name, weight in weights.items()}
        return result


class ThrottleFlowNumericalAttacks(unittest.TestCase):
    def assert_scalar_close(self, actual, expected):
        self.assertTrue(isfinite(actual))
        if expected == 0:
            self.assertEqual(actual, 0.0)
        else:
            self.assertNotEqual(actual, 0.0)
            self.assertAlmostEqual(actual, expected, delta=max(8 * ulp(expected), abs(expected) * 2e-14))

    def assert_ledger(self, model, throttle=1.0, scale=1.0):
        reference = decimal_ledger(model, throttle, scale, 240)
        self.assertEqual(decimal_ledger(model, throttle, scale, 160), reference)
        value = evaluate(model, throttle, scale)
        partials = value.partials
        self.assert_scalar_close(value.mass_derivative_kg_s, reference["total"])
        self.assert_scalar_close(partials.mass_rate_wrt_throttle, reference["command"])
        self.assert_scalar_close(dict(partials.mass_rate_wrt_parameters)["mass_flow_scale"], reference["driver"])
        self.assert_scalar_close(value.vacuum_specific_impulse_s, reference["isp"])
        for actual, key in ((value.additional_state_derivatives, "tank"),
                            (partials.additional_state_rates_wrt_throttle, "tank_command"),
                            (dict(partials.additional_state_rates_wrt_parameters)["mass_flow_scale"],
                             "tank_driver")):
            self.assertEqual(set(dict(actual)), set(reference[key]))
            for name, rate in actual:
                self.assert_scalar_close(rate, reference[key][name])
        return value

    def assert_domain_field(self, field, call):
        with self.assertRaises(DomainError) as caught:
            call()
        self.assertIn(field, caught.exception.field or "")

    def test_unrepresentable_total_has_a_typed_mass_output_error(self):
        self.assert_domain_field("Total mass flow", lambda: evaluate(engine(TINY, 1e-300), .1))

    def test_positive_half_subnormal_streams_cannot_be_reported_as_zero(self):
        self.assert_domain_field("Tank mass flow", lambda: evaluate(
            engine(TINY, 1e-300, (("a", .5), ("b", .5)))))
        # A deliberately zero fraction remains a valid exact-zero stream.
        value = evaluate(engine(TINY, 1e-300, (("a", 1.), ("b", 0.))))
        self.assertEqual(dict(value.additional_state_derivatives), {"a": -TINY, "b": 0.})

    def test_near_unity_fractions_preserve_total_and_relative_stream_weights(self):
        fractions = (("a", .5), ("b", .5 - 2.0**-42))
        for order in (fractions, tuple(reversed(fractions))):
            with self.subTest(order=order):
                value = self.assert_ledger(engine(2.0**50, fractions=order))
                self.assertEqual(fsum(rate for _, rate in value.additional_state_derivatives),
                                 value.mass_derivative_kg_s)

    def test_command_partial_recovers_underflowed_intermediate_product(self):
        curve = MonotoneCubicCurve((.25, .75), (.875, 1.125))  # .75 + .5 u
        model = engine(TINY, 1e-300, curve=curve)
        value = self.assert_ledger(model, .5, 2.)
        self.assertEqual(value.partials.mass_rate_wrt_throttle, -TINY)
        self.assertEqual(dict(value.partials.additional_state_rates_wrt_throttle)["tank"], -TINY)

    def test_command_partial_recovers_overflowed_intermediate_product(self):
        curve = MonotoneCubicCurve((.25, .75), (.1, 1.1))
        value = self.assert_ledger(engine(1e308, curve=curve), .5, .01)
        self.assertGreater(abs(value.partials.mass_rate_wrt_throttle), 1e306)

    def test_isp_uses_exact_flow_before_subnormal_export_rounding(self):
        model = engine(1e-320, 1e-300)
        values = [self.assert_ledger(model, command).vacuum_specific_impulse_s
                  for command in (.1, .3, 1.)]
        self.assertEqual(values[0], values[1])
        self.assertEqual(values[1], values[2])

    def test_public_ledger_matches_decimal_command_paths_and_normalized_splits(self):
        for flow, thrust in ((9.5, 30000.), (1e308, 150.), (1e-300, 1e-280)):
            for command, scale in ((.1, .01), (.55, 1.), (.75, 2.)):
                if flow == 1e308 and scale == 2.:
                    continue  # covered as an explicit output range failure below
                with self.subTest(flow=flow, command=command, scale=scale):
                    self.assert_ledger(engine(flow, thrust, (("ox", .78), ("fuel", .22))),
                                       command, scale)

    def test_true_unrepresentable_command_and_driver_partials_fail_explicitly(self):
        curve = MonotoneCubicCurve((.25, .75), (.1, 1.1))
        self.assert_domain_field("Throttle mass-flow partial", lambda: evaluate(
            engine(1e308, curve=curve), .5))
        high = MonotoneCubicCurve((.25, .75), (2., 2.))
        self.assert_domain_field("Mass-flow scale partial", lambda: evaluate(
            engine(1e308, curve=high), .5, .01))

    def test_mass_closure_preserves_small_inventory_terms_in_every_order(self):
        items = (("large", .5), ("small_a", .25), ("small_b", .25))
        states = {"large": 1e16, "small_a": 1., "small_b": 1.}
        for order in permutations(items):
            with self.subTest(order=order):
                closure = engine(fractions=order).mass_closure(
                    vehicle_mass_kg=1e16 + 2., additional_states=states, non_propellant_mass_kg=0.)
                self.assertEqual(closure.tracked_propellant_kg, 1e16 + 2.)
                self.assertEqual(closure.residual_kg, 0.)

    def test_mass_closure_residual_is_computed_before_rounding_the_tracked_sum(self):
        closure = engine(fractions=(("large", .5), ("small", .5))).mass_closure(
            vehicle_mass_kg=1e16, additional_states={"large": 1e16, "small": 1.},
            non_propellant_mass_kg=0.)
        self.assertEqual(closure.tracked_propellant_kg, 1e16)
        self.assertEqual(closure.residual_kg, -1.)
        closure = engine().mass_closure(vehicle_mass_kg=1e16 + 2., additional_states={"tank": 1e16},
                                        non_propellant_mass_kg=1.)
        self.assertEqual(closure.residual_kg, 1.)

    def test_mass_closure_range_errors_are_typed_and_exact_zero_is_valid(self):
        self.assert_domain_field("Relative mass closure residual", lambda: engine().mass_closure(
            vehicle_mass_kg=1e308, additional_states={"tank": TINY}, non_propellant_mass_kg=1e308))
        self.assert_domain_field("Tracked propellant mass", lambda: engine(
            fractions=(("a", .5), ("b", .5))).mass_closure(
                vehicle_mass_kg=1e308, additional_states={"a": 1e308, "b": 1e308},
                non_propellant_mass_kg=0.))
        closure = engine().mass_closure(vehicle_mass_kg=1e308, additional_states={"tank": 0.},
                                        non_propellant_mass_kg=1e308)
        self.assertEqual(closure.residual_kg, 0.)
        self.assertEqual(closure.relative_residual, 0.)

    def test_off_state_does_not_evaluate_an_unrepresentable_firing_ledger(self):
        value = evaluate(engine(TINY, 1e-300, (("a", .5), ("b", .5))), 0.)
        self.assertFalse(value.active)
        self.assertEqual(value.mass_derivative_kg_s, 0.)
        self.assertTrue(all(rate == 0 for _, rate in value.additional_state_derivatives))
        self.assertIsNone(value.vacuum_specific_impulse_s)

    def test_serialized_numerical_policies_cannot_be_changed_or_omitted(self):
        for name, model in direction_consumers().items():
            document = json.loads(model.to_json())
            self.assertEqual(type(model).from_dict(document).model_hash, model.model_hash)
            for key in (key for key in document if key.endswith("_numerics")):
                for change in ("replace", "omit"):
                    modified = dict(document)
                    if change == "replace":
                        modified[key] = "unsupported_numerical_policy"
                    else:
                        del modified[key]
                    with self.subTest(model=name, field=key, change=change):
                        with self.assertRaises(InputError) as caught:
                            type(model).from_dict(modified)
                        self.assertEqual(caught.exception.field, key)

    def test_varying_command_integrates_analytic_flow_and_rocket_equation(self):
        model = engine(2., fractions=(("a", .5), ("b", .5 - 2.0**-42)))
        state = (100., 10., 10., 0.)
        step = .02

        def derivative(time, values):
            mass, a, b, velocity = values
            value = model.evaluate(throttle=.4 + .2 * time, vehicle_mass_kg=mass,
                                   direction=(1., 0., 0.), additional_states={"a": a, "b": b})
            drains = dict(value.additional_state_derivatives)
            return (value.mass_derivative_kg_s, drains["a"], drains["b"], value.acceleration_m_s2[0])

        def shifted(values, rates, factor):
            return tuple(value + factor * rate for value, rate in zip(values, rates, strict=True))

        for index in range(100):
            time = index * step
            first = derivative(time, state)
            second = derivative(time + step / 2, shifted(state, first, step / 2))
            third = derivative(time + step / 2, shifted(state, second, step / 2))
            fourth = derivative(time + step, shifted(state, third, step))
            state = tuple(value + step * (a + 2*b + 2*c + d) / 6
                          for value, a, b, c, d in zip(state, first, second, third, fourth, strict=True))
        self.assertAlmostEqual(state[0], 97.6, delta=2e-12)
        self.assertAlmostEqual(state[0] - state[1] - state[2], 80., delta=2e-12)
        with localcontext() as context:
            context.prec = 80
            weights = [Decimal.from_float(value) for _, value in model.stream_fractions]
            for actual, weight in zip(state[1:3], weights, strict=True):
                expected = float(Decimal(10) - Decimal("2.4") * weight / sum(weights))
                self.assertAlmostEqual(actual, expected, delta=2e-12)
        self.assertAlmostEqual(state[3], 75 * log(100 / 97.6), delta=2e-11)
