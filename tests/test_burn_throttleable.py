"""Throttle-commanded engine model and its monotone cubic throttle curves."""

import copy
import random
import unittest
from math import isclose, log

from rocket_propulsion.core.errors import DomainError, InputError
from rocket_propulsion.propulsion.burns import (
    AdditionalStateConstraint,
    MonotoneCubicCurve,
    ThrottleableEngineModel,
)
from rocket_propulsion.propulsion.performance import STANDARD_GRAVITY_M_S2


class MonotoneCubicCurveTests(unittest.TestCase):
    def test_interpolates_nodes_and_reproduces_linear_data(self) -> None:
        curve = MonotoneCubicCurve((0.0, 1.0, 3.0, 4.0), (1.0, 3.0, 7.0, 9.0))
        for x, y in zip(curve.x, curve.y, strict=True):
            self.assertAlmostEqual(curve.evaluate(x)[0], y)
        for x in (0.3, 1.7, 2.9, 3.5):
            value, slope = curve.evaluate(x)
            self.assertAlmostEqual(value, 1.0 + 2.0 * x)
            self.assertAlmostEqual(slope, 2.0)

    def test_monotone_data_never_overshoot_and_slope_is_continuous(self) -> None:
        generator = random.Random(4)
        for _ in range(30):
            x = [0.0]
            y = [generator.uniform(0, 1)]
            for _ in range(generator.randint(1, 8)):
                x.append(x[-1] + generator.uniform(0.05, 1.0))
                y.append(y[-1] + generator.choice((0.0, generator.uniform(0.0, 2.0))))
            curve = MonotoneCubicCurve(tuple(x), tuple(y))
            previous = None
            for k in range(2001):
                value, slope = curve.evaluate(x[0] + (x[-1] - x[0]) * k / 2000)
                self.assertGreaterEqual(slope, -1e-12)
                if previous is not None:
                    self.assertGreaterEqual(value, previous - 1e-12)
                previous = value
            self.assertLessEqual(previous, y[-1] + 1e-12)
            for knot in x[1:-1]:
                left = curve.evaluate(knot - 1e-9)[1]
                right = curve.evaluate(knot + 1e-9)[1]
                self.assertAlmostEqual(left, right, places=5)

    def test_derivative_matches_finite_difference_and_no_extrapolation(self) -> None:
        curve = MonotoneCubicCurve((0.4, 0.6, 0.8, 1.0, 1.1), (0.37, 0.58, 0.79, 1.0, 1.09))
        h = 1e-7
        for x in (0.45, 0.61, 0.79, 0.95, 1.05):
            numeric = (curve.evaluate(x + h)[0] - curve.evaluate(x - h)[0]) / (2 * h)
            self.assertAlmostEqual(curve.evaluate(x)[1], numeric, places=6)
        with self.assertRaises(DomainError):
            curve.evaluate(1.2)
        with self.assertRaises(DomainError):
            MonotoneCubicCurve((0.0, 0.0), (1.0, 2.0))
        with self.assertRaises(DomainError):
            MonotoneCubicCurve((0.0,), (1.0,))


def engine(**kwargs: object) -> ThrottleableEngineModel:
    arguments = {
        "rated_thrust_n": 30_000.0,
        "rated_mass_flow_kg_s": 9.5,
        "minimum_throttle": 0.4,
        "maximum_throttle": 1.05,
        "stream_fractions": (("m_ox", 0.78), ("m_fu", 0.22)),
        # Isp falls at deep throttle: thrust fraction below the flow fraction.
        "thrust_curve": MonotoneCubicCurve(
            (0.4, 0.6, 0.8, 1.0, 1.05), (0.37, 0.58, 0.79, 1.0, 1.05)
        ),
        "state_constraints": (AdditionalStateConstraint("m_ox", 5.0),),
        "protected_dry_mass_kg": 900.0,
        "nozzle_exit_area_m2": 0.3,
    }
    arguments.update(kwargs)
    return ThrottleableEngineModel(**arguments)  # type: ignore[arg-type]


STATES = {"m_ox": 400.0, "m_fu": 120.0}


class ThrottleableEngineTests(unittest.TestCase):
    def evaluate(self, model, **kwargs):
        arguments = {
            "throttle": 0.7,
            "vehicle_mass_kg": 1500.0,
            "direction": (0.2, -1.0, 2.0),
            "additional_states": STATES,
            "ambient_pressure_pa": 2000.0,
        }
        arguments.update(kwargs)
        return model.evaluate(**arguments)

    def test_force_flow_isp_and_drain_closure(self) -> None:
        model = engine()
        value = self.evaluate(model)
        fraction = model.thrust_curve.evaluate(0.7)[0]
        self.assertAlmostEqual(value.thrust_n, 30_000.0 * fraction - 2000.0 * 0.3)
        self.assertAlmostEqual(value.mass_derivative_kg_s, -9.5 * 0.7)
        drains = dict(value.additional_state_derivatives)
        self.assertAlmostEqual(sum(drains.values()), value.mass_derivative_kg_s)
        self.assertAlmostEqual(drains["m_ox"], -9.5 * 0.7 * 0.78)
        self.assertAlmostEqual(
            value.vacuum_specific_impulse_s,
            30_000.0 * fraction / (STANDARD_GRAVITY_M_S2 * 9.5 * 0.7),
        )
        rated = self.evaluate(model, throttle=1.0)
        self.assertGreater(rated.vacuum_specific_impulse_s, value.vacuum_specific_impulse_s)

    def test_every_partial_matches_central_finite_differences(self) -> None:
        model = engine()
        nominal = self.evaluate(model)
        partials = nominal.partials

        def check(vector, kwargs_plus, kwargs_minus, step, places=7):
            upper = self.evaluate(model, **kwargs_plus)
            lower = self.evaluate(model, **kwargs_minus)
            for index in range(3):
                self.assertAlmostEqual(
                    vector[index],
                    (upper.acceleration_m_s2[index] - lower.acceleration_m_s2[index]) / (2 * step),
                    places=places,
                )
            return upper, lower

        h = 1e-6
        upper, lower = check(partials.acceleration_wrt_throttle,
                             {"throttle": 0.7 + h}, {"throttle": 0.7 - h}, h)
        self.assertAlmostEqual(
            partials.mass_rate_wrt_throttle,
            (upper.mass_derivative_kg_s - lower.mass_derivative_kg_s) / (2 * h), places=6,
        )
        rates_u = dict(partials.additional_state_rates_wrt_throttle)
        for name in STATES:
            self.assertAlmostEqual(
                rates_u[name],
                (dict(upper.additional_state_derivatives)[name]
                 - dict(lower.additional_state_derivatives)[name]) / (2 * h), places=6,
            )
        check(partials.acceleration_wrt_mass, {"vehicle_mass_kg": 1500.001},
              {"vehicle_mass_kg": 1499.999}, 0.001, places=9)
        check(partials.acceleration_wrt_ambient_pressure, {"ambient_pressure_pa": 2010.0},
              {"ambient_pressure_pa": 1990.0}, 10.0, places=12)
        for column in range(3):
            plus = [0.2, -1.0, 2.0]
            minus = [0.2, -1.0, 2.0]
            plus[column] += h
            minus[column] -= h
            column_vector = tuple(partials.acceleration_wrt_direction[r][column] for r in range(3))
            check(column_vector, {"direction": tuple(plus)}, {"direction": tuple(minus)}, h,
                  places=6)
        acc_p = dict(partials.acceleration_wrt_parameters)
        mass_p = dict(partials.mass_rate_wrt_parameters)
        for name in ("thrust_scale", "mass_flow_scale"):
            upper, lower = check(acc_p[name], {"parameter_overrides": {name: 1 + h}},
                                 {"parameter_overrides": {name: 1 - h}}, h, places=5)
            self.assertAlmostEqual(
                mass_p[name], (upper.mass_derivative_kg_s - lower.mass_derivative_kg_s) / (2 * h),
                places=6,
            )

    def test_throttle_limits_off_state_and_inventory(self) -> None:
        model = engine()
        off = self.evaluate(model, throttle=0.0)
        self.assertFalse(off.active)
        self.assertEqual(off.acceleration_m_s2, (0.0, 0.0, 0.0))
        for bad in (0.2, -0.1, 1.2, float("nan"), True):
            with self.assertRaises(DomainError):
                self.evaluate(model, throttle=bad)
        self.assertEqual(model.validated_throttle(0.4 - 1e-13), 0.4)
        self.assertEqual(model.validated_throttle(1e-13), 0.0)
        with self.assertRaises(DomainError):
            self.evaluate(model, ambient_pressure_pa=None)
        no_area = engine(nozzle_exit_area_m2=None)
        with self.assertRaises(DomainError):
            self.evaluate(no_area)
        self.assertTrue(self.evaluate(no_area, ambient_pressure_pa=None).active)
        limited = self.evaluate(model, vehicle_mass_kg=900.0)
        self.assertEqual(limited.inhibited_by, ("protected_dry_mass_kg",))
        self.assertTrue(self.evaluate(model, vehicle_mass_kg=900.0,
                                      apply_inventory_limits=False).active)
        self.assertEqual([s.event_id for s in model.inventory_surfaces()],
                         ["engine:dry-mass", "engine:reserve:m_ox"])

    def test_construction_refusals(self) -> None:
        for kwargs in (
            {"minimum_throttle": 0.0},
            {"minimum_throttle": 1.1},
            {"stream_fractions": (("m_ox", 0.7), ("m_fu", 0.2))},
            {"stream_fractions": (("m_ox", 0.5), ("m_ox", 0.5))},
            {"thrust_curve": MonotoneCubicCurve((0.5, 1.05), (0.5, 1.05))},
            {"thrust_curve": MonotoneCubicCurve((0.4, 1.05), (0.0, 1.05))},
            {"state_constraints": (AdditionalStateConstraint("m_x", 1.0),)},
            {"rated_thrust_n": -1.0},
        ):
            with self.subTest(kwargs), self.assertRaises(DomainError):
                engine(**kwargs)

    def test_persistence_round_trip_and_tamper_detection(self) -> None:
        model = engine()
        reopened = ThrottleableEngineModel.from_json(model.to_json())
        self.assertEqual(reopened.model_hash, model.model_hash)
        self.assertEqual(self.evaluate(reopened), self.evaluate(model))
        document = model.manifest()
        for label, mutate in {
            "curve": lambda d: d["thrust_curve"]["y"].__setitem__(1, 0.59),
            "limit": lambda d: d.__setitem__("minimum_throttle", 0.35),
            "fractions": lambda d: d["stream_fractions"][0].__setitem__(1, 0.8),
        }.items():
            tampered = copy.deepcopy(document)
            mutate(tampered)
            with self.subTest(label), self.assertRaises(InputError):
                ThrottleableEngineModel.from_dict(tampered)

    def test_constant_throttle_integration_obeys_rocket_equation(self) -> None:
        model = engine(nozzle_exit_area_m2=None)
        mass, velocity, h = 1500.0, 0.0, 0.5
        u = 0.55

        def rate(m):
            value = model.evaluate(throttle=u, vehicle_mass_kg=m, direction=(1.0, 0.0, 0.0),
                                   additional_states=STATES)
            return value.acceleration_m_s2[0], value.mass_derivative_kg_s, value

        for _ in range(60):
            a1, m1, _ = rate(mass)
            a2, m2, _ = rate(mass + h / 2 * m1)
            a3, m3, _ = rate(mass + h / 2 * m2)
            a4, m4, _ = rate(mass + h * m3)
            velocity += h / 6 * (a1 + 2 * a2 + 2 * a3 + a4)
            mass += h / 6 * (m1 + 2 * m2 + 2 * m3 + m4)
        exhaust = rate(mass)[2].vacuum_specific_impulse_s * STANDARD_GRAVITY_M_S2
        self.assertAlmostEqual(mass, 1500.0 - 9.5 * u * 30.0, places=9)
        self.assertTrue(isclose(velocity, exhaust * log(1500.0 / mass), rel_tol=1e-10))


if __name__ == "__main__":
    unittest.main()
