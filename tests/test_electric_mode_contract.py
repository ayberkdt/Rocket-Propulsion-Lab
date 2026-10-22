"""Independent budget, relay, event-coordinate and frontier counterexamples."""

import sys
import unittest
from decimal import Decimal, localcontext
from fractions import Fraction
from math import ldexp, nextafter

from rocket_propulsion.core.errors import DomainError, InputError
from rocket_propulsion.propulsion.burns import (
    DiscreteThrottlePath,
    ElectricThrottleLevel,
    ElectricThrottleTable,
    FixedEfficiencyThrottle,
    PowerLimitedPropulsionModel,
    PowerSwitchSurface,
)

MAX = sys.float_info.max


def path(powers=(1.0, 2.0, 3.0), hysteresis=0.0):
    return DiscreteThrottlePath(
        ElectricThrottleTable(
            tuple(ElectricThrottleLevel(f"L{i}", p, p, p) for i, p in enumerate(powers)),
            "synthetic-mode",
            "a" * 64,
        ),
        tuple(f"L{i}" for i in range(len(powers))),
        hysteresis,
    )


def fixed(minimum=1.0, hysteresis=0.0):
    return FixedEfficiencyThrottle(0.6, 2500.0, minimum, MAX, hysteresis)


def requirements(law):
    return (
        ([level.input_power_w for level in law.levels], law.up_switch_hysteresis_w)
        if isinstance(law, DiscreteThrottlePath)
        else ([law.minimum_power_w], law.start_hysteresis_w)
    )


def reference(law, power, current=None, margin=0.0):
    # Enumerate possible destinations with an independent decimal ledger.
    with localcontext() as context:
        context.prec = 1700
        req, hysteresis = requirements(law)
        powers = [Decimal.from_float(x) for x in req]
        supplied = Decimal.from_float(power) - Decimal.from_float(margin)
        feasible = [i for i, demand in enumerate(powers) if demand <= supplied]
        if current is None or current >= 0 and powers[current] > supplied:
            return max(feasible, default=-1)
        upgrades = [
            i
            for i in range(current + 1, len(powers))
            if powers[i] + Decimal.from_float(hysteresis) <= supplied
        ]
        return max(upgrades, default=current)


def visible_power(threshold, crossing):
    value = float(threshold)
    if crossing == "increasing":
        return nextafter(value, float("inf")) if Fraction(value) < threshold else value
    return nextafter(value, -float("inf")) if Fraction(value) >= threshold else value


class ElectricModeAttacks(unittest.TestCase):
    def test_zero_power_cannot_fund_positive_demand_at_any_scale(self):
        for power in (ldexp(1.0, -1074), 1e-12, 1.0, 1e100):
            for law in (path((power,)), fixed(power)):
                for current in (None, -1, 0):
                    with self.subTest(power=power, law=type(law).__name__, current=current):
                        self.assertEqual(law.select_level(0.0, current), -1)

    def test_neighboring_power_inputs_use_exact_budget_without_absolute_floor(self):
        for demand in (ldexp(1.0, -1074), 1e-12, 1.0, 1e100):
            for law in (path((demand,)), fixed(demand)):
                for power in (nextafter(demand, 0.0), demand, nextafter(demand, float("inf"))):
                    for current in (None, -1, 0):
                        self.assertEqual(
                            law.select_level(power, current), reference(law, power, current)
                        )

    def test_fractional_hysteresis_sum_is_not_rounded_before_decision(self):
        for law in (path((ldexp(1.0, -54),), 1.0), fixed(ldexp(1.0, -54), 1.0)):
            self.assertEqual(law.select_level(1.0, -1), -1)
            self.assertEqual(law.select_level(nextafter(1.0, float("inf")), -1), 0)

    def test_model_margin_is_not_rounded_before_budget_comparison(self):
        for law in (path((1.0,)), fixed()):
            model = PowerLimitedPropulsionModel(law, power_margin_w=ldexp(1.0, -54))
            self.assertEqual(model.select_level(1.0), -1)
            self.assertEqual(model.select_level(nextafter(1.0, float("inf"))), 0)

    def test_decimal_enumeration_matches_every_state_and_margin(self):
        for hysteresis in (0.0, 0.25, 100.0):
            for law in (path(hysteresis=hysteresis), fixed(hysteresis=hysteresis)):
                for margin in (0.0, ldexp(1.0, -54), 1e16):
                    model = PowerLimitedPropulsionModel(law, power_margin_w=margin)
                    demands, h = requirements(law)
                    with localcontext() as context:
                        context.prec = 1700
                        thresholds = [
                            float(
                                Decimal.from_float(d)
                                + Decimal.from_float(margin)
                                + Decimal.from_float(extra)
                            )
                            for d in demands
                            for extra in (0.0, h)
                        ]
                    inputs = {0.0, -1.0, MAX}
                    for threshold in thresholds:
                        inputs.update(
                            (
                                nextafter(threshold, -float("inf")),
                                threshold,
                                nextafter(threshold, float("inf")),
                            )
                        )
                    for power in inputs:
                        for current in (None, *range(-1, len(demands))):
                            with self.subTest(
                                h=hysteresis, margin=margin, power=power, current=current
                            ):
                                self.assertEqual(
                                    model.select_level(power, current),
                                    reference(law, power, current, margin),
                                )

    def test_initialization_and_known_off_hysteresis_are_explicitly_different(self):
        for law in (path(hysteresis=100.0), fixed(hysteresis=100.0)):
            self.assertEqual(law.select_level(1.5, None), 0)
            self.assertEqual(law.select_level(1.5, -1), -1)
            self.assertEqual(law.select_level(1.5, 0), 0)
        self.assertEqual(path(hysteresis=100.0).select_level(1.5, 2), 0)

    def test_off_decisions_do_not_export_unrepresentable_effective_power(self):
        model = PowerLimitedPropulsionModel(fixed(), power_margin_w=MAX)
        self.assertEqual(model.select_level(-MAX), -1)
        with self.assertRaises(DomainError):
            model.effective_power_w(-MAX)
        value = model.evaluate(
            level_index=-1,
            available_power_w=-MAX,
            vehicle_mass_kg=1.0,
            direction=(1.0, 0.0, 0.0),
            additional_states={"xenon_mass_kg": 1.0},
        )
        self.assertFalse(value.active)
        self.assertEqual(value.consumed_power_w, 0.0)

    def test_exact_margin_prevents_false_saturation_in_held_mode(self):
        law = fixed()
        model = PowerLimitedPropulsionModel(law, power_margin_w=1.0)
        value = model.evaluate(
            level_index=0,
            available_power_w=MAX,
            vehicle_mass_kg=1.0,
            direction=(1.0, 0.0, 0.0),
            additional_states={"xenon_mass_kg": 1.0},
        )
        self.assertGreater(value.partials.acceleration_wrt_available_power[0], 0.0)
        self.assertEqual([event.crossing for event in model.switching_surfaces(0)], ["decreasing"])

    def test_event_residual_uses_true_sum_not_rounded_display(self):
        for law in (path((1.0,)), fixed()):
            model = PowerLimitedPropulsionModel(law, power_margin_w=1e16)
            event = model.switching_surfaces(-1)[0]
            self.assertEqual(event.threshold_w, 1e16)
            self.assertEqual(event.value(1e16), -1.0)
            self.assertEqual(event.value(nextafter(1e16, float("inf"))), 1.0)

    def test_root_metadata_and_decimal_signs_include_hysteresis_and_margin(self):
        for law in (path(hysteresis=0.25), fixed(hysteresis=0.25)):
            model = PowerLimitedPropulsionModel(law, power_margin_w=1e16)
            for current in range(-1, law.level_count):
                for event in model.switching_surfaces(current):
                    ratio = getattr(event, "threshold_ratio", None)
                    self.assertIsNotNone(ratio)
                    threshold = Fraction(*ratio)
                    for power in (
                        nextafter(event.threshold_w, 0.0),
                        event.threshold_w,
                        nextafter(event.threshold_w, float("inf")),
                    ):
                        with localcontext() as context:
                            context.prec = 1700
                            expected = Decimal.from_float(power) - Decimal(
                                threshold.numerator
                            ) / Decimal(threshold.denominator)
                        self.assertEqual(event.value(power), float(expected))

    def test_colocated_up_guards_coalesce_to_the_highest_reachable_mode(self):
        model = PowerLimitedPropulsionModel(path(hysteresis=0.25), power_margin_w=1e16)
        event = next(e for e in model.switching_surfaces(0) if e.crossing == "increasing")
        self.assertEqual(event.to_level, 2)
        threshold = Fraction(*event.threshold_ratio)
        power = visible_power(threshold, "increasing")
        self.assertEqual(power, 1e16 + 4.0)
        self.assertEqual(model.select_level(nextafter(power, 0.0), 0), 0)
        self.assertEqual(model.select_level(power, 0), event.to_level)

    def test_colocated_down_guards_do_not_leave_an_infeasible_intermediate_mode(self):
        model = PowerLimitedPropulsionModel(path(), power_margin_w=1e16)
        event = next(e for e in model.switching_surfaces(1) if e.crossing == "decreasing")
        self.assertEqual(event.to_level, -1)
        self.assertEqual(model.select_level(1e16, 1), -1)

    def test_every_directional_target_matches_first_visible_crossed_budget(self):
        for h in (0.0, 0.25, 100.0):
            for law in (path(hysteresis=h), fixed(hysteresis=h)):
                for margin in (0.0, ldexp(1.0, -54), 1e16):
                    model = PowerLimitedPropulsionModel(law, power_margin_w=margin)
                    for current in range(-1, law.level_count):
                        for event in model.switching_surfaces(current):
                            if event.crossing == "either":
                                self.assertEqual(event.to_level, current)
                                continue
                            power = visible_power(Fraction(*event.threshold_ratio), event.crossing)
                            self.assertEqual(event.to_level, reference(law, power, current, margin))
                            self.assertEqual(
                                model.select_level(power, event.to_level), event.to_level
                            )

    def test_sampled_bidirectional_trace_handles_zero_start_and_plateau_without_chatter(self):
        model = PowerLimitedPropulsionModel(path(hysteresis=0.25), power_margin_w=1e16)
        powers = [1e16 + delta for delta in (0.0, 2.0, 4.0, 6.0, 4.0, 2.0, 2.0, 0.0, 0.0)]
        current = model.select_level(powers[0])
        modes = [current]
        event_names = []
        for before, after in zip(powers, powers[1:]):
            fired = []
            for event in model.switching_surfaces(current):
                left, right = event.value(before), event.value(after)
                # Sampled budget convention: upward equality enters, downward
                # equality remains feasible until the next strictly lower sample.
                if (
                    event.crossing == "increasing"
                    and left < 0 <= right
                    or event.crossing == "decreasing"
                    and left >= 0 > right
                ):
                    fired.append(event)
            self.assertLessEqual(len(fired), 1)
            if fired:
                event_names.append(fired[0].event_id)
                current = fired[0].to_level
            modes.append(current)
        self.assertEqual(modes, [-1, 0, 2, 2, 2, 1, 1, -1, -1])
        self.assertEqual(
            event_names, ["power:L0:up", "power:L1:up", "power:L2:down", "power:L1:down"]
        )
        consumed = 0.0
        for power, mode in zip(powers, modes):
            value = model.evaluate(
                level_index=mode,
                available_power_w=power,
                vehicle_mass_kg=100.0,
                direction=(1.0, 0.0, 0.0),
                additional_states={"xenon_mass_kg": 20.0},
            )
            self.assertEqual(value.power_deficit_w, 0.0)
            self.assertEqual(
                value.mass_derivative_kg_s,
                dict(value.additional_state_derivatives)["xenon_mass_kg"],
            )
            self.assertLessEqual(value.consumed_power_w, power - 1e16)
            consumed -= value.mass_derivative_kg_s  # one-second samples
        self.assertEqual(consumed, 14.0)

    def test_unreachable_start_roots_are_omitted_without_false_activation(self):
        for law in (path((ldexp(1.0, -1074),), MAX), fixed(ldexp(1.0, -1074), MAX)):
            self.assertEqual(law.select_level(MAX, -1), -1)
            self.assertEqual(law.switching_surfaces(-1), ())
        model = PowerLimitedPropulsionModel(path((ldexp(1.0, -1074),)), power_margin_w=MAX)
        self.assertEqual(model.select_level(MAX), -1)
        self.assertEqual(model.switching_surfaces(-1), ())
        with self.assertRaises(DomainError):
            model.switching_surfaces(0)

    def test_root_residual_overflow_is_not_an_infinite_solver_value(self):
        root = PowerSwitchSurface("generic", -MAX, "increasing", -1, 0)
        with self.assertRaises(DomainError):
            root.value(MAX)

    def test_surface_literals_direction_and_exact_ratio_are_validated(self):
        class Label(str):
            pass

        event = PowerSwitchSurface(Label("event"), 1.0, Label("increasing"), -1, 0)
        self.assertIs(type(event.event_id), str)
        self.assertIs(type(event.crossing), str)
        for crossing, source, target in (
            ("increasing", 0, -1),
            ("decreasing", -1, 0),
            ("increasing", 0, 0),
        ):
            with self.subTest(crossing=crossing), self.assertRaises(DomainError):
                PowerSwitchSurface("bad", 1.0, crossing, source, target)
        for ratio in ((1, 0), (1, -1), (True, 1), (1.0, 1), (2, 1), [1, 1], (1, 1, 1)):
            with self.subTest(ratio=ratio), self.assertRaises(DomainError):
                PowerSwitchSurface("bad", 1.0, "increasing", -1, 0, threshold_ratio=ratio)

    def test_supported_exact_ratio_can_distinguish_a_decimal_display(self):
        root = PowerSwitchSurface("decimal", 0.1, "increasing", -1, 0, threshold_ratio=(1, 10))
        self.assertEqual(root.threshold_ratio, (1, 10))
        self.assertGreater(root.value(0.1), 0.0)
        self.assertLess(root.value(nextafter(0.1, 0.0)), 0.0)

    def test_hysteresis_and_bus_offset_do_not_change_semantic_event_identity(self):
        law = path(hysteresis=0.25)
        model = PowerLimitedPropulsionModel(law, power_margin_w=0.125)
        self.assertEqual(
            [e.event_id for e in law.switching_surfaces(0)],
            [e.event_id for e in model.switching_surfaces(0)],
        )
        for event in model.switching_surfaces(0):
            self.assertEqual(event.from_level, 0)

    def test_policy_is_hash_bound_and_legacy_omission_is_not_silent_upgrade(self):
        model = PowerLimitedPropulsionModel(fixed())
        raw = model.to_dict()
        self.assertIn("power_mode_numerics", raw)
        for bad in (None, "unsupported"):
            changed = dict(raw)
            if bad is None:
                changed.pop("power_mode_numerics")
            else:
                changed["power_mode_numerics"] = bad
            with self.assertRaises(InputError) as caught:
                type(model).from_dict(changed)
            self.assertEqual(caught.exception.field, "power_mode_numerics")
        self.assertEqual(type(model).from_json(model.to_json()).model_hash, model.model_hash)

    def test_modes_and_powers_refuse_boolean_nonfinite_and_nonintegral_inputs(self):
        for law in (path(), fixed()):
            for current in (True, False, 0.0, "0", -2, law.level_count):
                with self.subTest(current=current), self.assertRaises(DomainError):
                    law.select_level(1.0, current)
            for power in (None, True, "1", float("inf"), float("nan")):
                with self.subTest(power=repr(power)), self.assertRaises(DomainError):
                    law.select_level(power, None)


class ElectricFrontierAttacks(unittest.TestCase):
    def test_exact_isp_order_is_not_lost_when_displays_are_identical(self):
        low = ElectricThrottleLevel("low", 1000.0, 1.6613593041120234, 1.6177267169785168)
        high = ElectricThrottleLevel("high", 2000.0, 1.6613593041120231, 1.6177267169785166)
        self.assertEqual(low.specific_impulse_s, high.specific_impulse_s)
        with localcontext() as context:
            context.prec = 120
            self.assertGreater(
                Decimal.from_float(high.thrust_n) / Decimal.from_float(high.mass_flow_kg_s),
                Decimal.from_float(low.thrust_n) / Decimal.from_float(low.mass_flow_kg_s),
            )
        table = ElectricThrottleTable((low, high), "tie", "a" * 64)
        self.assertEqual(
            DiscreteThrottlePath.from_table(table, objective="specific_impulse").level_ids,
            ("low", "high"),
        )

    def test_unrepresentable_isp_need_not_be_exported_for_exact_ranking(self):
        levels = (
            ElectricThrottleLevel("a", 1e308, 1e-14, ldexp(1.0, -1074)),
            ElectricThrottleLevel(
                "b",
                nextafter(1e308, float("inf")),
                nextafter(1e-14, float("inf")),
                ldexp(1.0, -1074),
            ),
        )
        table = ElectricThrottleTable(levels, "ranking", "a" * 64)
        for level in levels:
            with self.assertRaises(DomainError):
                _ = level.specific_impulse_s
        try:
            result = DiscreteThrottlePath.from_table(table, objective="specific_impulse")
        except DomainError as error:
            self.fail(f"Ordering must not export unrepresentable Isp: {error}")
        self.assertEqual(result.level_ids, ("a", "b"))

    def test_invalid_table_or_forged_objective_is_typed(self):
        class Objective(str):
            def __hash__(self):
                return hash("thrust")

            def __eq__(self, other):
                return True

        with self.assertRaises(DomainError):
            DiscreteThrottlePath.from_table(path().table, objective=Objective("unrecognized"))
        for table in (None, {}, True):
            with self.subTest(table=table), self.assertRaises(DomainError):
                DiscreteThrottlePath.from_table(table)
