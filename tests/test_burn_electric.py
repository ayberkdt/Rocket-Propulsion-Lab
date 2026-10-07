"""Power-limited electric propulsion: tables, throttle paths, events, partials.

The throttle table below is synthetic.  Its levels are chosen to be
physically consistent (total efficiency below one) and to exercise the path
construction rules; they are not data for any real thruster.
"""

import unittest
from math import isclose, log

from rocket_propulsion.core.errors import DomainError, InputError
from rocket_propulsion.propulsion.burns import (
    DiscreteThrottlePath,
    ElectricThrottleLevel,
    ElectricThrottleTable,
    FixedEfficiencyThrottle,
    PowerLimitedPropulsionModel,
)
from rocket_propulsion.propulsion.performance import STANDARD_GRAVITY_M_S2

G0 = STANDARD_GRAVITY_M_S2


def synthetic_table() -> ElectricThrottleTable:
    return ElectricThrottleTable(
        levels=(
            ElectricThrottleLevel("L1", 1000.0, 0.040, 1.5e-6),
            ElectricThrottleLevel("L2", 2000.0, 0.080, 2.6e-6),
            ElectricThrottleLevel("L2-high-isp", 2000.0, 0.070, 1.9e-6),
            ElectricThrottleLevel("L-dominated", 2500.0, 0.075, 3.5e-6),
            ElectricThrottleLevel("L3", 3000.0, 0.110, 3.2e-6),
        ),
        source_id="synthetic-electric-table",
        source_sha256="f" * 64,
    )


def thrust_model(**kwargs: object) -> PowerLimitedPropulsionModel:
    path = DiscreteThrottlePath.from_table(synthetic_table(), objective="thrust")
    return PowerLimitedPropulsionModel(path, **kwargs)  # type: ignore[arg-type]


class ThrottleTableTests(unittest.TestCase):
    def test_derived_isp_and_efficiency(self) -> None:
        level = synthetic_table().level("L2")
        self.assertAlmostEqual(level.specific_impulse_s, 0.080 / (G0 * 2.6e-6))
        self.assertAlmostEqual(level.total_efficiency, 0.080**2 / (2 * 2.6e-6 * 2000.0))
        self.assertLess(level.total_efficiency, 1.0)

    def test_inconsistent_levels_and_tables_are_refused(self) -> None:
        with self.assertRaises(DomainError):
            ElectricThrottleLevel("impossible", 1000.0, 1.0, 1e-6)
        with self.assertRaises(DomainError):
            ElectricThrottleLevel("negative", 1000.0, -0.01, 1e-6)
        with self.assertRaises(DomainError):
            ElectricThrottleTable(
                levels=(
                    ElectricThrottleLevel("A", 1000.0, 0.04, 1.5e-6),
                    ElectricThrottleLevel("A", 2000.0, 0.08, 2.6e-6),
                ),
                source_id="dup",
                source_sha256="f" * 64,
            )
        with self.assertRaises(DomainError):
            ElectricThrottleTable(
                levels=(ElectricThrottleLevel("A", 1000.0, 0.04, 1.5e-6),),
                source_id="bad-hash",
                source_sha256="F" * 64,
            )

    def test_hash_round_trip_and_tamper_detection(self) -> None:
        table = synthetic_table()
        self.assertEqual(ElectricThrottleTable.from_json(table.to_json()), table)
        payload = table.to_dict()
        payload["levels"][0]["thrust_n"] = 0.041  # type: ignore[index]
        with self.assertRaises(InputError):
            ElectricThrottleTable.from_dict(payload)
        with self.assertRaises(InputError):
            ElectricThrottleTable.from_json("[]")


class ThrottlePathTests(unittest.TestCase):
    def test_paths_never_select_a_worse_level_with_more_power(self) -> None:
        table = synthetic_table()
        self.assertEqual(
            DiscreteThrottlePath.from_table(table, objective="thrust").level_ids,
            ("L1", "L2", "L3"),
        )
        self.assertEqual(
            DiscreteThrottlePath.from_table(table, objective="specific_impulse").level_ids,
            ("L1", "L2-high-isp"),
        )
        with self.assertRaises(DomainError):
            DiscreteThrottlePath.from_table(table, objective="cost")
        with self.assertRaises(DomainError):
            DiscreteThrottlePath(table, ("L2", "L2-high-isp"))
        with self.assertRaises(DomainError):
            DiscreteThrottlePath(table, ("L3", "L1"))

    def test_initial_selection_and_hysteresis(self) -> None:
        plain = thrust_model()
        self.assertEqual(
            [plain.select_level(power) for power in (500.0, 1000.0, 1500.0, 2999.0, 5000.0)],
            [-1, 0, 0, 1, 2],
        )
        model = PowerLimitedPropulsionModel(
            DiscreteThrottlePath.from_table(
                synthetic_table(), objective="thrust", up_switch_hysteresis_w=100.0
            )
        )
        self.assertEqual(model.select_level(2050.0, 0), 0)
        self.assertEqual(model.select_level(2100.0, 0), 1)
        self.assertEqual(model.select_level(2050.0, 1), 1)
        self.assertEqual(model.select_level(1999.0, 1), 0)
        self.assertEqual(model.select_level(999.0, 1), -1)
        self.assertEqual(model.select_level(1050.0, -1), -1)
        self.assertEqual(model.select_level(1100.0, -1), 0)
        self.assertEqual(model.select_level(9000.0, -1), 2)
        with self.assertRaises(DomainError):
            model.select_level(1500.0, 3)

    def test_switching_surfaces_have_deterministic_targets(self) -> None:
        model = PowerLimitedPropulsionModel(
            DiscreteThrottlePath.from_table(
                synthetic_table(), objective="thrust", up_switch_hysteresis_w=100.0
            ),
            power_margin_w=50.0,
        )
        middle = model.switching_surfaces(1)
        self.assertEqual(
            [(item.crossing, item.threshold_w, item.from_level, item.to_level) for item in middle],
            [("decreasing", 2050.0, 1, 0), ("increasing", 3150.0, 1, 2)],
        )
        self.assertEqual(middle[0].value(2000.0), -50.0)
        self.assertEqual(
            [(item.crossing, item.to_level) for item in model.switching_surfaces(2)],
            [("decreasing", 1)],
        )
        self.assertEqual(
            [(item.crossing, item.threshold_w, item.to_level) for item in model.switching_surfaces(-1)],
            [("increasing", 1150.0, 0)],
        )
        self.assertEqual(model.select_level(2049.0), 0)


class EvaluationTests(unittest.TestCase):
    def setUp(self) -> None:
        self.states = {"xenon_mass_kg": 20.0}

    def evaluate(self, model: PowerLimitedPropulsionModel, **kwargs: object):
        arguments = {
            "level_index": 1,
            "available_power_w": 2500.0,
            "vehicle_mass_kg": 500.0,
            "direction": (0.0, 3.0, 4.0),
            "additional_states": self.states,
        }
        arguments.update(kwargs)
        return model.evaluate(**arguments)  # type: ignore[arg-type]

    def test_level_force_drain_power_and_isp(self) -> None:
        value = self.evaluate(thrust_model())
        self.assertTrue(value.active)
        self.assertEqual(value.level_id, "L2")
        self.assertAlmostEqual(value.thrust_n, 0.080)
        for actual, expected in zip(value.force_vector_n, (0.0, 0.048, 0.064), strict=True):
            self.assertAlmostEqual(actual, expected)
        self.assertAlmostEqual(value.mass_derivative_kg_s, -2.6e-6)
        self.assertEqual(dict(value.additional_state_derivatives)["xenon_mass_kg"], -2.6e-6)
        self.assertEqual(value.consumed_power_w, 2000.0)
        self.assertEqual(value.power_deficit_w, 0.0)
        self.assertAlmostEqual(value.specific_impulse_s, 0.080 / (G0 * 2.6e-6))
        self.assertEqual(value.partials.acceleration_wrt_available_power, (0.0, 0.0, 0.0))

    def test_rk_stage_past_a_root_keeps_level_and_reports_deficit(self) -> None:
        value = self.evaluate(thrust_model(), available_power_w=1990.0)
        self.assertTrue(value.active)
        self.assertEqual(value.level_id, "L2")
        self.assertAlmostEqual(value.power_deficit_w, 10.0)

    def test_duty_cycle_off_mode_and_inventory_limits(self) -> None:
        half = self.evaluate(thrust_model(duty_cycle=0.5))
        self.assertAlmostEqual(half.thrust_n, 0.040)
        self.assertAlmostEqual(half.mass_derivative_kg_s, -1.3e-6)
        self.assertAlmostEqual(half.consumed_power_w, 1000.0)
        self.assertAlmostEqual(half.specific_impulse_s, 0.080 / (G0 * 2.6e-6))
        off = self.evaluate(thrust_model(), level_index=-1)
        self.assertFalse(off.active)
        self.assertEqual(off.acceleration_m_s2, (0.0, 0.0, 0.0))
        limited = thrust_model(propellant_reserve_kg=20.0, protected_dry_mass_kg=500.0)
        inhibited = self.evaluate(limited)
        self.assertEqual(inhibited.inhibited_by, ("protected_dry_mass_kg", "xenon_mass_kg"))
        unclamped = self.evaluate(limited, apply_inventory_limits=False)
        self.assertTrue(unclamped.active)
        self.assertEqual(
            [surface.event_id for surface in limited.inventory_surfaces()],
            ["electric:dry-mass", "electric:reserve:xenon_mass_kg"],
        )
        with self.assertRaises(DomainError):
            self.evaluate(thrust_model(), level_index=3)
        with self.assertRaises(DomainError):
            self.evaluate(thrust_model(), additional_states={})
        with self.assertRaises(DomainError):
            thrust_model(duty_cycle=0.0)

    def test_partials_match_central_finite_differences(self) -> None:
        model = PowerLimitedPropulsionModel(
            FixedEfficiencyThrottle(0.6, 2500.0, 600.0, 4000.0)
        )
        nominal = self.evaluate(model, level_index=0, direction=(0.3, -1.2, 2.0))
        step_mass, step_power, step_dir, step_param = 1e-3, 1e-2, 1e-6, 1e-6

        def acceleration(**kwargs: object) -> tuple[float, float, float]:
            arguments: dict[str, object] = {"level_index": 0, "direction": (0.3, -1.2, 2.0)}
            arguments.update(kwargs)
            return self.evaluate(model, **arguments).acceleration_m_s2

        upper = acceleration(vehicle_mass_kg=500.0 + step_mass)
        lower = acceleration(vehicle_mass_kg=500.0 - step_mass)
        for index in range(3):
            self.assertAlmostEqual(
                nominal.partials.acceleration_wrt_mass[index],
                (upper[index] - lower[index]) / (2 * step_mass),
                places=12,
            )
        upper = acceleration(available_power_w=2500.0 + step_power)
        lower = acceleration(available_power_w=2500.0 - step_power)
        for index in range(3):
            self.assertAlmostEqual(
                nominal.partials.acceleration_wrt_available_power[index],
                (upper[index] - lower[index]) / (2 * step_power),
                places=12,
            )
        for column in range(3):
            plus = [0.3, -1.2, 2.0]
            minus = [0.3, -1.2, 2.0]
            plus[column] += step_dir
            minus[column] -= step_dir
            upper = acceleration(direction=tuple(plus))
            lower = acceleration(direction=tuple(minus))
            for row in range(3):
                self.assertAlmostEqual(
                    nominal.partials.acceleration_wrt_direction[row][column],
                    (upper[row] - lower[row]) / (2 * step_dir),
                    places=10,
                )
        flow_partials = dict(nominal.partials.mass_rate_wrt_parameters)
        acceleration_partials = dict(nominal.partials.acceleration_wrt_parameters)
        for name in ("thrust_scale", "mass_flow_scale"):
            plus = self.evaluate(
                model, level_index=0, direction=(0.3, -1.2, 2.0),
                parameter_overrides={name: 1.0 + step_param},
            )
            minus = self.evaluate(
                model, level_index=0, direction=(0.3, -1.2, 2.0),
                parameter_overrides={name: 1.0 - step_param},
            )
            self.assertAlmostEqual(
                flow_partials[name],
                (plus.mass_derivative_kg_s - minus.mass_derivative_kg_s) / (2 * step_param),
                places=12,
            )
            for index in range(3):
                self.assertAlmostEqual(
                    acceleration_partials[name][index],
                    (plus.acceleration_m_s2[index] - minus.acceleration_m_s2[index])
                    / (2 * step_param),
                    places=10,
                )
        mass_rate_power = (
            self.evaluate(model, level_index=0, available_power_w=2500.0 + step_power)
            .mass_derivative_kg_s
            - self.evaluate(model, level_index=0, available_power_w=2500.0 - step_power)
            .mass_derivative_kg_s
        ) / (2 * step_power)
        self.assertAlmostEqual(
            nominal.partials.mass_rate_wrt_available_power, mass_rate_power, places=14
        )

    def test_manifest_hash_binds_law_and_limits(self) -> None:
        base = thrust_model()
        self.assertEqual(base.model_hash, thrust_model().model_hash)
        self.assertEqual(base.manifest()["model_hash"], base.model_hash)
        self.assertNotEqual(base.model_hash, thrust_model(power_margin_w=10.0).model_hash)
        self.assertNotEqual(
            base.model_hash,
            PowerLimitedPropulsionModel(
                DiscreteThrottlePath.from_table(synthetic_table(), objective="specific_impulse")
            ).model_hash,
        )


class FixedEfficiencyTests(unittest.TestCase):
    law = FixedEfficiencyThrottle(0.6, 2500.0, 600.0, 4000.0, start_hysteresis_w=50.0)

    def test_jet_power_identity_saturation_and_shutdown(self) -> None:
        point = self.law.operating_point(0, 2000.0)
        exhaust = G0 * 2500.0
        self.assertAlmostEqual(0.5 * point.mass_flow_kg_s * exhaust**2, 0.6 * 2000.0)
        self.assertAlmostEqual(point.thrust_n / point.mass_flow_kg_s, exhaust)
        self.assertAlmostEqual(point.thrust_wrt_available_power, 2 * 0.6 / exhaust)
        saturated = self.law.operating_point(0, 5000.0)
        self.assertAlmostEqual(saturated.input_power_w, 4000.0)
        self.assertEqual(saturated.thrust_wrt_available_power, 0.0)
        self.assertEqual(self.law.select_level(599.0, 0), -1)
        self.assertEqual(self.law.select_level(620.0, -1), -1)
        self.assertEqual(self.law.select_level(650.0, -1), 0)
        probe = PowerLimitedPropulsionModel(self.law, power_margin_w=100.0).evaluate(
            level_index=0,
            available_power_w=50.0,
            vehicle_mass_kg=500.0,
            direction=(1.0, 0.0, 0.0),
            additional_states={"xenon_mass_kg": 20.0},
        )
        self.assertEqual(probe.thrust_n, 0.0)
        self.assertIsNone(probe.specific_impulse_s)
        self.assertAlmostEqual(probe.power_deficit_w, 650.0)
        with self.assertRaises(DomainError):
            FixedEfficiencyThrottle(1.2, 2500.0, 600.0, 4000.0)
        with self.assertRaises(DomainError):
            FixedEfficiencyThrottle(0.6, 2500.0, 600.0, 600.0)

    def test_on_mode_registers_shutdown_and_saturation_kink(self) -> None:
        model = PowerLimitedPropulsionModel(self.law, power_margin_w=100.0)
        self.assertEqual(
            [(s.event_id, s.threshold_w, s.to_level) for s in model.switching_surfaces(0)],
            [
                ("power:fixed-efficiency:down", 700.0, -1),
                ("power:fixed-efficiency:saturation", 4100.0, 0),
            ],
        )
        self.assertEqual(
            [(s.threshold_w, s.to_level) for s in model.switching_surfaces(-1)],
            [(750.0, 0)],
        )


class PowerDrivenIntegrationTests(unittest.TestCase):
    """Event-driven RK4 with a linearly decaying power supply."""

    initial_power_w = 3500.0
    power_slope_w_s = -1.0  # stands in for the consumer's power model
    initial_mass_kg = 500.0
    initial_xenon_kg = 20.0

    def power(self, time_s: float) -> float:
        return self.initial_power_w + self.power_slope_w_s * time_s

    def integrate(
        self, model: PowerLimitedPropulsionModel, end_s: float, step_s: float
    ) -> tuple[list[float], list[tuple[str, float]]]:
        direction = (1.0, 0.0, 0.0)
        time_s = 0.0
        y = [0.0, 0.0, 0.0, self.initial_mass_kg, self.initial_xenon_kg]
        level = model.select_level(self.power(0.0))
        log_events: list[tuple[str, float]] = []

        def rate(t: float, state: list[float]) -> list[float]:
            value = model.evaluate(
                level_index=level,
                available_power_w=self.power(t),
                vehicle_mass_kg=state[3],
                direction=direction,
                additional_states={"xenon_mass_kg": state[4]},
                apply_inventory_limits=False,
            )
            return [
                *value.acceleration_m_s2,
                value.mass_derivative_kg_s,
                dict(value.additional_state_derivatives)["xenon_mass_kg"],
            ]

        def rk4(t: float, state: list[float], h: float) -> list[float]:
            k1 = rate(t, state)
            k2 = rate(t + h / 2, [a + h / 2 * b for a, b in zip(state, k1, strict=True)])
            k3 = rate(t + h / 2, [a + h / 2 * b for a, b in zip(state, k2, strict=True)])
            k4 = rate(t + h, [a + h * b for a, b in zip(state, k3, strict=True)])
            return [
                a + h / 6 * (b + 2 * c + 2 * d + e)
                for a, b, c, d, e in zip(state, k1, k2, k3, k4, strict=True)
            ]

        while time_s < end_s - 1e-12:
            h = min(step_s, end_s - time_s)
            surfaces = model.switching_surfaces(level)
            start_power = self.power(time_s)
            end_power = self.power(time_s + h)

            def fires(surface, before: float, after: float) -> bool:
                if surface.crossing == "decreasing":
                    return before > 0.0 >= after
                if surface.crossing == "increasing":
                    return before < 0.0 <= after
                return (before > 0.0 >= after) or (before < 0.0 <= after)

            roots = []
            for surface in surfaces:
                before = surface.value(start_power)
                if not fires(surface, before, surface.value(end_power)):
                    continue
                low, high = 0.0, h
                for _ in range(80):
                    middle = 0.5 * (low + high)
                    if fires(surface, before, surface.value(self.power(time_s + middle))):
                        high = middle
                    else:
                        low = middle
                roots.append((high, surface))
            if roots:
                h, surface = min(roots, key=lambda item: item[0])
                y = rk4(time_s, y, h)
                time_s += h
                log_events.append((surface.event_id, time_s))
                level = surface.to_level
                continue
            y = rk4(time_s, y, h)
            time_s += h
        return y, log_events

    def test_discrete_levels_switch_at_analytic_times_and_close_mass(self) -> None:
        model = thrust_model()
        final, events = self.integrate(model, 3200.0, 37.0)
        self.assertEqual(
            [name for name, _ in events],
            ["power:L3:down", "power:L2:down", "power:L1:down"],
        )
        for (_, time_s), expected in zip(events, (500.0, 1500.0, 2500.0), strict=True):
            self.assertAlmostEqual(time_s, expected, places=8)
        levels = DiscreteThrottlePath.from_table(synthetic_table()).levels[::-1]
        durations = (500.0, 1000.0, 1000.0)
        mass = self.initial_mass_kg
        speed = 0.0
        for level, duration in zip(levels, durations, strict=True):
            end_mass = mass - level.mass_flow_kg_s * duration
            speed += level.thrust_n / level.mass_flow_kg_s * log(mass / end_mass)
            mass = end_mass
        self.assertAlmostEqual(final[3], mass, places=10)
        self.assertTrue(isclose(final[0], speed, rel_tol=1e-10))
        closure = model.mass_closure(
            vehicle_mass_kg=final[3],
            additional_states={"xenon_mass_kg": final[4]},
            non_propellant_mass_kg=self.initial_mass_kg - self.initial_xenon_kg,
        )
        self.assertLess(abs(closure.residual_kg), 1e-10)

    def test_continuous_law_follows_power_and_obeys_rocket_equation(self) -> None:
        # Saturated at 3000 W until t=500 s, follows power down to the 1000 W
        # shutdown at t=2500 s, then off.
        law = FixedEfficiencyThrottle(0.6, 2500.0, 1000.0, 3000.0)
        model = PowerLimitedPropulsionModel(law)
        final, events = self.integrate(model, 3000.0, 41.0)
        self.assertEqual(
            [name for name, _ in events],
            ["power:fixed-efficiency:saturation", "power:fixed-efficiency:down"],
        )
        self.assertAlmostEqual(events[0][1], 500.0, places=8)
        self.assertAlmostEqual(events[1][1], 2500.0, places=8)
        exhaust = G0 * 2500.0
        energy_j = 3000.0 * 500.0 + 0.5 * (3000.0 + 1000.0) * 2000.0
        consumed = 2.0 * 0.6 / exhaust**2 * energy_j
        self.assertAlmostEqual(final[3], self.initial_mass_kg - consumed, places=10)
        expected_speed = exhaust * log(self.initial_mass_kg / final[3])
        self.assertTrue(isclose(final[0], expected_speed, rel_tol=1e-10))


if __name__ == "__main__":
    unittest.main()
