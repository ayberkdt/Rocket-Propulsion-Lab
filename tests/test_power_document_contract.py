"""Public raw schemas, literal identity and electric level arithmetic attacks."""

import json
import unittest
from dataclasses import replace
from decimal import Decimal, localcontext
from math import ldexp, nextafter
from random import Random

from rocket_propulsion.core.errors import DomainError, InputError
from rocket_propulsion.propulsion.burns import (
    DiscreteThrottlePath,
    ElectricThrottleLevel,
    ElectricThrottleTable,
    FixedEfficiencyThrottle,
    PowerLimitedPropulsionModel,
)
from rocket_propulsion.propulsion.performance import STANDARD_GRAVITY_M_S2
from tests.test_burn_model_persistence import configured_bridge, discrete_model
from tests.test_burn_throttleable import STATES, engine


def continuous():
    return PowerLimitedPropulsionModel(FixedEfficiencyThrottle(0.6, 2500.0, 600.0, 4000.0, 20.0))


def models():
    discrete = discrete_model()
    return (discrete.throttle_law.table, discrete, continuous(), engine())


def raw(model):
    return json.loads(model.to_json())


class PowerDocumentAttacks(unittest.TestCase):
    def refuse(self, model, mutate):
        document = raw(model)
        mutate(document)
        for encoding in ("mapping", "json"):
            with (
                self.subTest(model=type(model).__name__, encoding=encoding),
                self.assertRaises(InputError),
            ):
                if encoding == "mapping":
                    type(model).from_dict(document)
                else:
                    type(model).from_json(json.dumps(document))

    def test_unknown_root_fields_do_not_disappear(self):
        for model in models():
            self.refuse(model, lambda d: d.update(unmodeled={"units": "lbf"}))

    def test_unknown_law_fields_are_refused_for_both_laws(self):
        for model in (discrete_model(), continuous()):
            self.refuse(model, lambda d: d["throttle_law"].update(unmodeled=1))

    def test_fields_from_the_other_law_are_not_ignored(self):
        self.refuse(discrete_model(), lambda d: d["throttle_law"].update(total_efficiency=0.9))
        self.refuse(continuous(), lambda d: d["throttle_law"].update(table_hash="a" * 64))

    def test_unknown_table_level_fields_are_refused_standalone_and_nested(self):
        model = discrete_model()
        self.refuse(model.throttle_law.table, lambda d: d["levels"][0].update(thrust_unit="lbf"))
        self.refuse(model, lambda d: d["throttle_table"]["levels"][0].update(thrust_unit="lbf"))

    def test_unknown_parameter_and_constraint_fields_are_refused(self):
        for model in (discrete_model(), continuous(), engine()):
            self.refuse(model, lambda d: d["parameters"][0].update(prior="not-applied"))
        self.refuse(engine(), lambda d: d["state_constraints"][0].update(action="warn-only"))

    def test_unknown_curve_fields_are_refused(self):
        for key in ("thrust_curve", "flow_curve"):
            model = engine(flow_curve=engine().thrust_curve)
            self.refuse(model, lambda d, *, key=key: d[key].update(extrapolate=True))

    def test_declared_references_and_power_semantics_are_checked(self):
        for model in (discrete_model(), continuous(), engine()):
            self.refuse(model, lambda d: d.__setitem__("reference_ids", ["unreviewed-source"]))
        for key in ("force_semantics", "mass_semantics", "power_semantics"):
            for model in (discrete_model(), continuous()):
                with self.subTest(key=key):
                    self.refuse(model, lambda d, *, key=key: d.__setitem__(key, "different-semantics"))

    def test_required_derived_metadata_cannot_be_deleted(self):
        for model in (discrete_model(), continuous(), engine()):
            self.refuse(model, lambda d: d.pop("reference_ids"))
        for key in ("force_semantics", "mass_semantics", "power_semantics"):
            self.refuse(discrete_model(), lambda d, *, key=key: d.pop(key))

    def test_continuous_law_cannot_discard_a_supplied_table(self):
        table = discrete_model().throttle_law.table.to_dict()
        self.refuse(continuous(), lambda d: d.__setitem__("throttle_table", table))
        self.refuse(continuous(), lambda d: d.__setitem__("throttle_table", {"invalid": True}))

    def test_empty_parameter_list_cannot_invent_drivers_in_any_model_family(self):
        for model in (configured_bridge(), discrete_model(), continuous(), engine()):
            self.refuse(model, lambda d: d.__setitem__("parameters", []))

    def test_positional_stream_pairs_are_not_reconstructed_from_mapping_keys(self):
        model = engine()
        document = raw(model)
        document["stream_fractions"] = [
            {name: "discarded-value", value: "discarded-value"}
            for name, value in document["stream_fractions"]
        ]
        with self.assertRaises(InputError):
            type(model).from_dict(document)

    def test_json_array_fields_are_not_tuple_shortcuts(self):
        for model, path in (
            (discrete_model(), ("parameters",)),
            (discrete_model(), ("throttle_law", "level_ids")),
            (engine(), ("stream_fractions",)),
            (engine(), ("thrust_curve", "x")),
        ):
            document = raw(model)
            target = document
            for key in path[:-1]:
                target = target[key]
            target[path[-1]] = tuple(target[path[-1]])
            with self.subTest(path=path), self.assertRaises(InputError):
                type(model).from_dict(document)

    def test_plain_json_number_types_are_required(self):
        for value in (True, "1", float("inf"), float("nan")):
            document = raw(discrete_model())
            document["duty_cycle"] = value
            with self.subTest(value=repr(value)), self.assertRaises(InputError):
                PowerLimitedPropulsionModel.from_dict(document)

    def test_canonical_roundtrips_preserve_behavior_and_events(self):
        for model in models():
            opened = type(model).from_json(model.to_json())
            self.assertEqual(opened, model)
            self.assertEqual(type(model).from_dict(raw(model)), model)
        for model in (discrete_model(), continuous()):
            opened = type(model).from_dict(raw(model))
            kwargs = {
                "level_index": 0,
                "available_power_w": 1800.0,
                "vehicle_mass_kg": 400.0,
                "direction": (1.0, 0.0, 0.0),
                "additional_states": {"xenon_mass_kg": 20.0},
            }
            self.assertEqual(opened.evaluate(**kwargs), model.evaluate(**kwargs))
            self.assertEqual(opened.switching_surfaces(0), model.switching_surfaces(0))
        model = engine()
        opened = type(model).from_dict(raw(model))
        kwargs = {
            "throttle": 0.7,
            "vehicle_mass_kg": 1500.0,
            "direction": (1.0, 0.0, 0.0),
            "additional_states": STATES,
            "ambient_pressure_pa": 10000.0,
        }
        self.assertEqual(opened.evaluate(**kwargs), model.evaluate(**kwargs))

    def test_optional_null_fields_and_null_curves_remain_supported(self):
        for model in (
            continuous(),
            engine(
                thrust_curve=None,
                flow_curve=None,
                protected_dry_mass_kg=None,
                nozzle_exit_area_m2=None,
            ),
        ):
            self.assertEqual(type(model).from_dict(raw(model)), model)

    def test_parameter_and_numeric_representation_are_not_byte_authentication(self):
        for model in (discrete_model(), continuous(), engine()):
            document = raw(model)
            document = dict(reversed(list(document.items())))
            for parameter in document["parameters"]:
                parameter["reference_value"] = int(parameter["reference_value"])
            self.assertEqual(type(model).from_dict(document).to_json(), model.to_json())
        changed = replace(continuous(), model_id="another-valid-id")
        self.assertNotEqual(changed.model_hash, continuous().model_hash)
        self.assertEqual(type(changed).from_dict(raw(changed)), changed)

    def test_programmatic_metadata_is_stored_as_plain_literal_text(self):
        class Label(str):
            pass

        table = replace(
            discrete_model().throttle_law.table,
            source_id=Label(" Δv data "),
            source_sha256=Label("a" * 64),
            schema=Label("rocket_propulsion_electric_throttle_table_v1"),
            levels=(ElectricThrottleLevel(Label("L"), 1000.0, 0.04, 1.5e-6),),
        )
        path = DiscreteThrottlePath(table, (Label("L"),), path_id=Label("path"))
        law = FixedEfficiencyThrottle(0.6, 2500.0, 600.0, 4000.0, law_id=Label("law"))
        power = PowerLimitedPropulsionModel(
            law,
            model_id=Label("model"),
            propellant_state_name=Label("tank"),
            schema=Label("rocket_propulsion_power_limited_propulsion_v1"),
        )
        throttle = engine(
            model_id=Label("model"), schema=Label("rocket_propulsion_throttleable_engine_v1")
        )
        values = (
            table.source_id,
            table.source_sha256,
            table.schema,
            table.levels[0].level_id,
            path.path_id,
            path.level_ids[0],
            law.law_id,
            power.model_id,
            power.propellant_state_name,
            power.schema,
            throttle.model_id,
            throttle.schema,
        )
        for value in values:
            self.assertIs(type(value), str)
        self.assertEqual(table.source_id, " Δv data ")

    def test_programmatic_forged_schema_comparisons_cannot_change_literal_identity(self):
        class Forged(str):
            def __eq__(self, other):
                return True

            def __ne__(self, other):
                return False

        for model in (discrete_model().throttle_law.table, continuous(), engine()):
            with self.subTest(model=type(model).__name__), self.assertRaises(DomainError):
                replace(model, schema=Forged("wrong-schema"))

    def test_programmatic_digest_validation_reads_literal_characters(self):
        class Digest(str):
            def __iter__(self):
                return iter("a" * 64)

        with self.assertRaises(DomainError):
            replace(discrete_model().throttle_law.table, source_sha256=Digest("z" * 64))


class ElectricLevelNumericalAttacks(unittest.TestCase):
    def test_consistency_ceiling_neighbors_use_exact_binary_input_ratio(self):
        generator = Random(20261020)
        outcomes = []
        with localcontext() as context:
            context.prec = 120
            ceiling = Decimal.from_float(1.0 + 1e-12)
            for exponent in (-1000, -600, -1, 0, 400, 900):
                for _ in range(12):
                    power = ldexp(generator.uniform(1.0, 1.8), exponent)
                    flow = ldexp(generator.uniform(1.0, 1.8), exponent)
                    candidate = float(
                        (2 * ceiling * Decimal.from_float(power) * Decimal.from_float(flow)).sqrt()
                    )
                    for thrust in (
                        nextafter(candidate, 0.0),
                        candidate,
                        nextafter(candidate, float("inf")),
                    ):
                        ratio = Decimal.from_float(thrust) ** 2 / (
                            2 * Decimal.from_float(power) * Decimal.from_float(flow)
                        )
                        expected = ratio <= ceiling
                        try:
                            ElectricThrottleLevel("neighbor", power, thrust, flow)
                        except DomainError:
                            accepted = False
                        else:
                            accepted = True
                        self.assertEqual(accepted, expected, (exponent, power, thrust, flow))
                        outcomes.append(accepted)
        self.assertIn(True, outcomes)
        self.assertIn(False, outcomes)

    def test_efficiency_and_isp_are_invariant_under_exact_binary_unit_scaling(self):
        expected = 1.0 / STANDARD_GRAVITY_M_S2
        for exponent in (-1074, -1022, -900, -600, -10, 0, 200, 512, 900, 1023):
            value = ldexp(1.0, exponent)
            with self.subTest(exponent=exponent):
                level = ElectricThrottleLevel("scale", value, value, value)
                self.assertEqual(level.total_efficiency, 0.5)
                self.assertEqual(level.specific_impulse_s, expected)

    def test_large_denominator_does_not_turn_nonzero_isp_into_zero(self):
        level = ElectricThrottleLevel("subnormal-isp", 1e308, 1.0, 1e308)
        with localcontext() as context:
            context.prec = 100
            expected = float(
                Decimal(1) / (Decimal.from_float(STANDARD_GRAVITY_M_S2) * Decimal.from_float(1e308))
            )
        self.assertGreater(expected, 0.0)
        self.assertEqual(level.specific_impulse_s, expected)

    def test_consistency_decision_is_made_before_rounded_efficiency_export(self):
        with self.assertRaises(DomainError):
            ElectricThrottleLevel("impossible", 1e-100, 1e-200, 1e-308)
        # A legitimate tiny ratio need not be representable as a float until requested.
        level = ElectricThrottleLevel("tiny-efficiency", 1.0, 1e-200, 1.0)
        with self.assertRaises(DomainError):
            _ = level.total_efficiency

    def test_unrepresentable_isp_outputs_fail_typed_without_physical_floors(self):
        for power, thrust, flow in (
            (1.0, ldexp(1.0, -1074), 1e308),
            (1e308, 1e-14, ldexp(1.0, -1074)),
        ):
            level = ElectricThrottleLevel("export", power, thrust, flow)
            with self.subTest(thrust=thrust), self.assertRaises(DomainError):
                _ = level.specific_impulse_s


class ContinuousPowerNumericalAttacks(unittest.TestCase):
    def test_model_recovers_representable_outputs_across_private_point_boundary(self):
        law = FixedEfficiencyThrottle(ldexp(1.0, -1074), 2.5 / STANDARD_GRAVITY_M_S2, 1.0, 1e308)
        # A complete public point cannot export its tiny unscaled flow slope.
        with self.assertRaises(DomainError):
            law.operating_point(0, 1e307)
        model = PowerLimitedPropulsionModel(law)
        try:
            value = model.evaluate(
                level_index=0,
                available_power_w=1e307,
                vehicle_mass_kg=1.0,
                direction=(1.0, 0.0, 0.0),
                additional_states={"xenon_mass_kg": 1.0},
                parameter_overrides={"mass_flow_scale": 2.0},
            )
        except DomainError as error:
            self.fail(
                f"All final outputs have representable Decimal references; premature point export: {error}"
            )
        with localcontext() as context:
            context.prec = 120
            efficiency = Decimal.from_float(law.total_efficiency)
            velocity = Decimal.from_float(STANDARD_GRAVITY_M_S2) * Decimal.from_float(
                law.specific_impulse_s
            )
            gain = 2 * efficiency / velocity
            flow_gain = 2 * gain / velocity
            self.assertEqual(value.thrust_n, float(gain * Decimal.from_float(1e307)))
            self.assertEqual(
                value.mass_derivative_kg_s, float(-flow_gain * Decimal.from_float(1e307))
            )
            self.assertEqual(value.partials.mass_rate_wrt_available_power, float(-flow_gain))
            self.assertEqual(
                value.partials.acceleration_wrt_available_power, (float(gain), 0.0, 0.0)
            )
        self.assertNotEqual(value.partials.mass_rate_wrt_available_power, 0.0)

    def test_negative_python_indices_do_not_select_last_thruster(self):
        law = discrete_model().throttle_law
        for index in (-1, -2):
            with self.subTest(index=index), self.assertRaises(DomainError):
                law.operating_point(index, 1800.0)

    def expected(self, law, power):
        with localcontext() as context:
            context.prec = 120
            velocity = Decimal.from_float(STANDARD_GRAVITY_M_S2) * Decimal.from_float(
                law.specific_impulse_s
            )
            gain = 2 * Decimal.from_float(law.total_efficiency) / velocity
            thrust = gain * Decimal.from_float(power)
            return float(thrust), float(thrust / velocity)

    def test_continuous_large_isp_avoids_false_zero_thrust_and_mass_flow(self):
        law = FixedEfficiencyThrottle(1.0, 1e308, 1.0, 1e308)
        point = law.operating_point(0, 1e308)
        thrust, flow = self.expected(law, 1e308)
        self.assertGreater(thrust, 0.0)
        self.assertGreater(flow, 0.0)
        self.assertEqual(point.thrust_n, thrust)
        self.assertEqual(point.mass_flow_kg_s, flow)
        self.assertEqual(point.thrust_wrt_available_power, 0.0)

    def test_unrepresentable_gain_can_recover_at_saturation_before_export(self):
        law = FixedEfficiencyThrottle(ldexp(1.0, -1074), 1.0, 1.0, 1e308)
        with self.assertRaises(DomainError):
            _ = law.thrust_per_watt_n_w
        point = law.operating_point(0, 1e308)
        thrust, flow = self.expected(law, 1e308)
        self.assertEqual(point.thrust_n, thrust)
        self.assertEqual(point.mass_flow_kg_s, flow)
        self.assertGreater(point.thrust_n, 0.0)
        self.assertEqual(point.mass_flow_wrt_available_power, 0.0)

    def test_unrepresentable_required_unsaturated_slope_fails_typed(self):
        law = FixedEfficiencyThrottle(ldexp(1.0, -1074), 1.0, 1.0, 1e308)
        with self.assertRaises(DomainError):
            law.operating_point(0, 1e307)

    def test_clipped_negative_power_has_zero_local_derivatives(self):
        law = continuous().throttle_law
        point = law.operating_point(0, -10.0)
        self.assertEqual(point.thrust_n, 0.0)
        self.assertEqual(point.mass_flow_kg_s, 0.0)
        self.assertEqual(point.thrust_wrt_available_power, 0.0)
        self.assertEqual(point.mass_flow_wrt_available_power, 0.0)
        self.assertEqual(
            law.operating_point(0, 0.0).thrust_wrt_available_power, law.thrust_per_watt_n_w
        )

    def test_public_operating_points_validate_mode_and_power_inputs(self):
        for law in (discrete_model().throttle_law, continuous().throttle_law):
            for mode in (False, True, -1, law.level_count, 0.0, None, "0"):
                with (
                    self.subTest(law=type(law).__name__, mode=mode),
                    self.assertRaises(DomainError),
                ):
                    law.operating_point(mode, 1800.0)
            for power in (True, None, "1800", float("nan"), float("inf")):
                with (
                    self.subTest(law=type(law).__name__, power=repr(power)),
                    self.assertRaises(DomainError),
                ):
                    law.operating_point(0, power)

    def test_model_does_not_turn_positive_tank_drain_into_zero(self):
        level = ElectricThrottleLevel("tiny-flow", 1e308, 1e-160, ldexp(1.0, -1074))
        table = ElectricThrottleTable((level,), "synthetic", "a" * 64)
        model = PowerLimitedPropulsionModel(
            DiscreteThrottlePath(table, ("tiny-flow",)), duty_cycle=0.5
        )
        with self.assertRaises(DomainError):
            model.evaluate(
                level_index=0,
                available_power_w=1e308,
                vehicle_mass_kg=1.0,
                direction=(1.0, 0.0, 0.0),
                additional_states={"xenon_mass_kg": 1.0},
            )

    def test_model_isp_uses_unrounded_flow_and_exhaust_denominator(self):
        level = ElectricThrottleLevel("scale", 1e308, 1e308, 1e308)
        table = ElectricThrottleTable((level,), "synthetic", "a" * 64)
        model = PowerLimitedPropulsionModel(DiscreteThrottlePath(table, ("scale",)))
        value = model.evaluate(
            level_index=0,
            available_power_w=1e308,
            vehicle_mass_kg=1.0,
            direction=(1.0, 0.0, 0.0),
            additional_states={"xenon_mass_kg": 1.0},
        )
        self.assertEqual(value.specific_impulse_s, 1.0 / STANDARD_GRAVITY_M_S2)
        self.assertEqual(value.mass_derivative_kg_s, -1e308)
