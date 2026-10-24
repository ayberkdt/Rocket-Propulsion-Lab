"""Inventory/decision attacks with independent high-precision Decimal ledgers.

All inputs are synthetic. Closure of a numerical ledger is not physical validity.
"""

import json
import unittest
from dataclasses import asdict
from decimal import Decimal, localcontext
from itertools import permutations
from math import nextafter
from sys import float_info

from rocket_propulsion.core.errors import DomainError, InputError
from rocket_propulsion.propulsion.burns import (
    PropagatorPropulsionBridge,
    PropellantMassClosure,
    StreamStateBinding,
    TabulatedBurnArtifact,
    TabulatedBurnSegment,
    ThrottleableEngineModel,
)
from tests.test_burn_direction_numerics import engines

TINY = float.fromhex("0x0.0000000000001p-1022")


def many_tanks(names=("a", "b"), *, alias=False):
    streams = tuple((name, 1.) for name in names)
    artifact = TabulatedBurnArtifact((TabulatedBurnSegment(
        0., 1., 150., 150., 150., 150., float(len(names)), float(len(names)),
        streams, streams, "synthetic"),), "synthetic-inventory", "a"*64)
    binding = tuple(StreamStateBinding(name, "shared" if alias else name) for name in names)
    model = PropagatorPropulsionBridge(artifact, stream_bindings=binding)
    return {"bridge": model, "throttle": ThrottleableEngineModel(
        150., 1., .1, 1., tuple((name, 1/len(names)) for name in names))}


def decimal_ledger(vehicle, other, inventories, precision):
    with localcontext() as context:
        context.prec = precision
        v, dry = map(Decimal.from_float, (vehicle, other))
        tracked = sum(map(Decimal.from_float, inventories), Decimal(0))
        residual = v-dry-tracked
        numerator, denominator = residual.as_integer_ratio()
        return (float(tracked), float(residual), float(residual/v),
                (str(numerator), str(denominator)))


class InventoryClosureAttacks(unittest.TestCase):
    def assert_ledger(self, model, vehicle, other, states, selected=None):
        inventories = tuple(states[name] for name in (selected or states))
        reference = decimal_ledger(vehicle, other, inventories, 2000)
        self.assertEqual(decimal_ledger(vehicle, other, inventories, 1600), reference)
        value = model.mass_closure(vehicle_mass_kg=vehicle, non_propellant_mass_kg=other,
                                   additional_states=states)
        self.assertEqual((value.tracked_propellant_kg, value.residual_kg, value.relative_residual), reference[:3])
        return value

    def assert_field(self, label, call):
        with self.assertRaises(DomainError) as caught:
            call()
        self.assertIn(label, caught.exception.field or "")

    def test_all_models_preserve_the_residual_and_reject_false_closed_decision(self):
        for name, model in engines().items():
            with self.subTest(model=name):
                value = self.assert_ledger(model, 1e16+2., 1., {"tank": 1e16})
                self.assertEqual(value.residual_kg, 1.)
                self.assertFalse(value.is_closed(absolute_tolerance_kg=.1))
                self.assertTrue(value.is_closed(absolute_tolerance_kg=1.))

    def test_inventory_sum_rounding_does_not_erase_a_small_residual(self):
        for name, model in many_tanks().items():
            with self.subTest(model=name):
                value = self.assert_ledger(model, 1e16, 0., {"a": 1e16, "b": 1.})
                self.assertEqual(value.tracked_propellant_kg, 1e16)
                self.assertEqual(value.residual_kg, -1.)

    def test_shared_bound_tank_is_counted_once_and_unbound_states_are_ignored(self):
        model = many_tanks(alias=True)["bridge"]
        value = self.assert_ledger(model, 10., 1., {"shared": 9., "ignored": 1e308}, ("shared",))
        self.assertTrue(value.is_closed(absolute_tolerance_kg=0.))
        for name, model in engines().items():
            with self.subTest(model=name):
                self.assert_ledger(model, 10., 1., {"tank": 9., "ignored": 1e308}, ("tank",))

    def test_signed_trial_inventory_sum_recovers_finite_result_after_intermediate_overflow(self):
        for order in permutations(("a", "b", "c")):
            for name, model in many_tanks(order).items():
                with self.subTest(model=name, order=order):
                    value = self.assert_ledger(model, 1e308, 0., {"a": 1e308, "b": 1e308, "c": -1e308})
                    self.assertEqual(value.residual_kg, 0.)

    def test_true_unrepresentable_inventory_and_residual_fields_fail_explicitly(self):
        for name, model in many_tanks().items():
            with self.subTest(model=name, field="tracked"):
                self.assert_field("Tracked propellant mass", lambda*, model=model: model.mass_closure(
                    vehicle_mass_kg=1e308, non_propellant_mass_kg=0., additional_states={"a": 1e308, "b": 1e308}))
        for name, model in engines().items():
            with self.subTest(model=name, field="residual"):
                self.assert_field("Mass closure residual", lambda*, model=model: model.mass_closure(
                    vehicle_mass_kg=1e308, non_propellant_mass_kg=0., additional_states={"tank": -1e308}))

    def test_relative_residual_range_errors_are_not_reported_as_zero_or_infinity(self):
        for name, model in engines().items():
            for vehicle, other, tracked in ((1e308, 1e308, TINY), (TINY, 0., 1e308)):
                with self.subTest(model=name, vehicle=vehicle):
                    self.assert_field("Relative mass closure residual", lambda*, model=model, vehicle=vehicle, other=other, tracked=tracked: model.mass_closure(
                        vehicle_mass_kg=vehicle, non_propellant_mass_kg=other, additional_states={"tank": tracked}))

    def test_decision_uses_unrounded_ratio_when_displayed_records_coincide(self):
        inventories = (
            ({"a": .25, "b": .25}, True),
            ({"a": nextafter(.5, 0.), "b": 3*2.**-56}, False),
        )
        for name, model in many_tanks().items():
            displays = []
            records = []
            for states, closed in inventories:
                with self.subTest(model=name, closed=closed):
                    value = self.assert_ledger(model, 2., .5, states)
                    displays.append((value.tracked_propellant_kg, value.residual_kg, value.relative_residual))
                    records.append(value)
                    self.assertEqual(value.is_closed(absolute_tolerance_kg=1.), closed)
            self.assertEqual(displays[0], displays[1])
            self.assertNotEqual(records[0], records[1])

    def test_exact_and_below_boundary_residuals_remain_closed(self):
        for name, model in engines().items():
            for other, tracked, closed in ((.5, .5, True), (1., 2.**-55, True),
                                          (nextafter(1., 0.), 2.**-55, False)):
                with self.subTest(model=name, other=other):
                    value = self.assert_ledger(model, 2., other, {"tank": tracked})
                    self.assertEqual(value.is_closed(absolute_tolerance_kg=1.), closed)

    def test_ratio_decision_does_not_require_a_finite_outward_float_bound(self):
        for name, model in engines().items():
            with self.subTest(model=name):
                value = self.assert_ledger(model, float_info.max, 0., {"tank": -TINY})
                self.assertEqual(value.residual_kg, float_info.max)
                self.assertFalse(value.is_closed(absolute_tolerance_kg=float_info.max))

    def test_tolerances_reject_boolean_malformed_and_nonfinite_values(self):
        value = engines()["bridge"].mass_closure(
            vehicle_mass_kg=10., non_propellant_mass_kg=1., additional_states={"tank": 9.})
        for tolerance in (True, False, None, "0.1", 1j, float("nan"), float("inf"), -1.):
            with self.subTest(tolerance=tolerance), self.assertRaises(DomainError):
                value.is_closed(absolute_tolerance_kg=tolerance)

    def test_public_records_validate_fields_and_copy_the_numerical_ratio(self):
        ordinary = PropellantMassClosure(2., .5, .5, 1., .5)
        self.assertTrue(ordinary.is_closed(absolute_tolerance_kg=1.))
        payload = asdict(ordinary)
        for field, bad in (("vehicle_mass_kg", True), ("vehicle_mass_kg", 0.),
                           ("non_propellant_mass_kg", -1.), ("tracked_propellant_kg", float("inf")),
                           ("residual_kg", float("nan")), ("relative_residual", .25)):
            with self.subTest(field=field), self.assertRaises(DomainError):
                PropellantMassClosure(**{**payload, field: bad})
        for ratio in (("1", "0"), ("nan", "1"), ("0", "1"), (True, 1), ["1"]):
            with self.subTest(ratio=ratio), self.assertRaises(DomainError):
                PropellantMassClosure(**{**payload, "residual_ratio": ratio})
        ratio = ["1", "1"]
        value = PropellantMassClosure(2., .5, .5, 1., .5, ratio)
        ratio[0] = "0"
        self.assertEqual(value.residual_ratio, ("1", "1"))
        self.assertFalse(value.is_closed(absolute_tolerance_kg=.5))

    def test_exact_decision_survives_strict_json_round_trip(self):
        model = many_tanks()["bridge"]
        value = model.mass_closure(vehicle_mass_kg=2., non_propellant_mass_kg=.5,
                                  additional_states={"a": nextafter(.5, 0.), "b": 3*2.**-56})
        expected = decimal_ledger(2., .5, (nextafter(.5, 0.), 3*2.**-56), 2000)
        self.assertEqual(value.residual_ratio, expected[3])
        raw = json.loads(json.dumps(asdict(value), allow_nan=False))
        reopened = PropellantMassClosure(**raw)
        self.assertEqual(value, reopened)
        self.assertFalse(reopened.is_closed(absolute_tolerance_kg=1.))

    def test_all_model_identities_bind_the_same_inventory_policy(self):
        policies = {model.manifest().get("mass_closure_numerics") for model in engines().values()}
        self.assertEqual(len(policies), 1)
        self.assertNotIn(None, policies)
        for name, model in engines().items():
            document = json.loads(model.to_json())
            self.assertEqual(type(model).from_dict(document).model_hash, model.model_hash)
            for mode in ("replace", "omit"):
                modified = dict(document)
                if mode == "replace":
                    modified["mass_closure_numerics"] = "unsupported"
                else:
                    del modified["mass_closure_numerics"]
                with self.subTest(model=name, mode=mode), self.assertRaises(InputError):
                    type(model).from_dict(modified)
