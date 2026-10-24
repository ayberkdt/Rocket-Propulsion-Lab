"""Adversarial programmatic normalization: no erased values or reordered data."""

import unittest
from array import array
from collections.abc import Mapping
from dataclasses import dataclass, field, make_dataclass
from fractions import Fraction
from math import copysign, inf, nextafter
from typing import Optional, Union

from rocket_propulsion.core import validation as v
from rocket_propulsion.core.errors import DomainError, InputError
from rocket_propulsion.propulsion.burns import TabulatedBurnArtifact, TabulatedBurnSegment


@dataclass(frozen=True)
class Record:
    value: float


class LegacySequence:
    def __init__(self, values):
        self.values = values

    def __getitem__(self, index):
        return self.values[index]


class BrokenIterator:
    def __len__(self):
        return 3

    def __iter__(self):
        raise TypeError("broken data iterator")


ANNOTATION_CALLS = []


def annotation_probe():
    ANNOTATION_CALLS.append(True)
    return float


@dataclass(frozen=True)
class OpaqueAnnotation:
    value: "annotation_probe()"


class CoreValidationAttacks(unittest.TestCase):
    def test_stream_metadata_uses_literal_normalization_for_non_depleting_paths(self):
        from dataclasses import replace

        from tests.test_burn_constant_numerics import mono

        definition, point = mono()
        base = point.streams[0]

        class ForgedEncoding(str):
            def encode(self, *args, **kwargs):
                return b"fabricated"

        class ForgedRole(str):
            def __hash__(self):
                return hash(base.role.value)

            def __eq__(self, other):
                return True

            def __ne__(self, other):
                return False

        for change in (
            {"tank_depleting": False, "tank_id": ForgedEncoding("\ud800")},
            {"role": ForgedRole("unsupported")},
        ):
            with self.subTest(change=repr(change)), self.assertRaises(DomainError):
                replace(base, **change)

        class Label(str):
            def __str__(self):
                return "fabricated"

        value = replace(base, tank_depleting=False, tank_id=Label(""), stream_id=Label("actual"))
        self.assertIs(type(value.tank_id), str)
        self.assertIs(type(value.stream_id), str)
        self.assertEqual(value.tank_id, "")
        self.assertEqual(value.stream_id, "actual")
        from rocket_propulsion.propulsion.burns import simulate_constant_burn

        result = simulate_constant_burn(definition, point)
        summary = replace(result.summary, binding_constraint=Label("duration_target"))
        self.assertIs(type(summary.binding_constraint), str)
        self.assertEqual(summary.binding_constraint, "duration_target")

    def test_text_normalizers_use_literal_data_not_subclass_overrides(self):
        class ForgedText(str):
            def __str__(self):
                return "fabricated"

            def __float__(self):
                return 99.0

            def strip(self, *args):
                return "fabricated"

            def encode(self, *args, **kwargs):
                return b"fabricated"

        value = v.text(ForgedText("actual"), "name")
        self.assertIs(type(value), str)
        self.assertEqual(value, "actual")
        for bad in (ForgedText(" "), ForgedText("\ud800")):
            with self.subTest(bad=repr(bad)), self.assertRaises(DomainError):
                v.text(bad, "name")
        names = v.strings([ForgedText("actual")], "names")
        self.assertEqual(names, ("actual",))
        self.assertIs(type(names[0]), str)
        self.assertEqual(v.finite_float_text(ForgedText("1.25"), "value"), 1.25)
        with self.assertRaises(InputError):
            v.finite_float_text(ForgedText("not-a-number"), "value")
        source = v.json_text(ForgedText('{"actual":1}'), "document")
        self.assertIs(type(source), str)
        self.assertEqual(source, '{"actual":1}')

    def test_tabulated_metadata_normalization_preserves_literal_model_identity(self):
        segment = TabulatedBurnSegment(0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0)

        class ForgedDigest(str):
            def lower(self):
                return "a" * 64

        class ForgedSchema(str):
            def __eq__(self, other):
                return True

            def __ne__(self, other):
                return False

        for options in (
            {"source_sha256": ForgedDigest("z" * 64)},
            {"schema": ForgedSchema("unsupported")},
            {"interpolation": ForgedSchema("unsupported")},
        ):
            fields = {
                "segments": (segment,),
                "source_id": "source",
                "source_sha256": "a" * 64,
                **options,
            }
            with self.subTest(options=options), self.assertRaises(DomainError):
                TabulatedBurnArtifact(**fields)

        class Label(str):
            def __str__(self):
                return "fabricated"

        segment = TabulatedBurnSegment(0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, phase=Label("phase"))
        value = TabulatedBurnArtifact(
            (segment,), Label("source"), Label("a" * 64), reference_ids=(Label("reference"),)
        )
        self.assertIs(type(value.source_id), str)
        self.assertIs(type(value.source_sha256), str)
        self.assertIs(type(value.segments[0].phase), str)
        self.assertEqual(value.source_id, "source")
        self.assertEqual(TabulatedBurnArtifact.from_json(value.to_json()), value)

    def test_existing_domain_errors_keep_type_and_structured_context(self):
        declared = InputError(
            "source sample is invalid", field="sensor.flow", details={"unit": "kg/s"}
        )

        class DomainFloat(float):
            def __float__(self):
                raise declared

        class DomainInteger(int):
            def __int__(self):
                raise declared

        class DomainIterator:
            def __iter__(self):
                raise declared

        class DomainMapping(dict):
            def items(self):
                raise declared

        for reader, value in (
            (v.real, DomainFloat(1.0)),
            (v.integer, DomainInteger(1)),
            (v.real_tuple, DomainIterator()),
            (v.real_mapping, DomainMapping()),
        ):
            with self.subTest(reader=reader):
                try:
                    reader(value, "outer")
                except DomainError as actual:
                    self.assertIs(type(actual), type(declared))
                    self.assertEqual(actual.as_dict(), declared.as_dict())
                else:
                    self.fail("Structured source error was suppressed")

    def domain_failure(self, call, label=None):
        """Assert the declared error type for recognized malformed data protocols."""
        try:
            call()
        except DomainError as error:
            if label is not None:
                self.assertEqual(error.field, label)
        except (TypeError, ValueError, OverflowError, KeyError, IndexError) as error:
            self.fail(f"Recognized bad-data protocol escaped as {type(error).__name__}: {error}")
        else:
            self.fail("Malformed input was accepted")

    def test_mapping_protocol_cannot_drop_duplicates_or_leak_data_errors(self):
        class InvalidMapping(Mapping):
            def __iter__(self):
                return iter(("flow",))

            def __len__(self):
                return 1

            def __getitem__(self, key):
                return 1.0

            def items(self):
                return [("flow", 1.0), ("flow", 2.0)]

        self.domain_failure(lambda: v.real_mapping(InvalidMapping(), "values"), "values")

        class BrokenMapping(InvalidMapping):
            def items(self):
                raise TypeError("invalid data mapping")

        self.domain_failure(lambda: v.real_mapping(BrokenMapping(), "values"), "values")

    def test_iterator_value_failures_are_typed_without_partial_output(self):
        for error in (TypeError, ValueError, KeyError, IndexError, OverflowError):

            def broken(*, error=error):
                yield 1.0
                raise error("bad data iteration")

            with self.subTest(error=error):
                self.domain_failure(lambda: v.real_tuple(broken(), "values"), "values")

    def test_dataclass_scalar_types_optional_none_and_literal_identity(self):
        cls = make_dataclass(
            "Mixed",
            [("mass", float), ("count", int), ("name", str), ("limit", Optional[float])],  # noqa: UP045 - Exercise legacy typing.Optional compatibility.
            frozen=True,
        )
        record = cls(1, 2, " ölçüm ", None)
        v.coerce_fields(record, "record")
        self.assertIs(type(record.mass), float)
        self.assertIs(type(record.count), int)
        self.assertEqual(record.name, " ölçüm ")
        self.assertIsNone(record.limit)
        for bad in (cls(1, True, "name", None), cls(1, 2, "\ud800", None)):
            with self.assertRaises(DomainError):
                v.coerce_fields(bad, "record")

    def test_nonzero_real_underflow_is_not_a_valid_zero(self):
        tiny = Fraction(nextafter(0.0, 1.0))
        for value in (tiny / 4, -tiny / 4, tiny / 2, -tiny / 2):
            with self.subTest(value=value):
                with self.assertRaises(DomainError) as raised:
                    v.real(value, "mass")
                self.assertEqual(raised.exception.field, "mass")

    def test_model_cannot_receive_a_false_zero_flow_from_fraction_input(self):
        tiny = Fraction(nextafter(0.0, 1.0)) / 4
        with self.assertRaises(DomainError):
            TabulatedBurnSegment(0.0, 1.0, 1.0, 1.0, 1.0, 1.0, tiny, tiny)

    def test_real_accepts_representable_subnormals_and_preserves_signed_zero(self):
        tiny = Fraction(nextafter(0.0, 1.0))
        for value, expected in (
            (tiny, float(tiny)),
            (tiny * Fraction(3, 4), float(tiny)),
            (-tiny, -float(tiny)),
            (Fraction(1, 8), 0.125),
            (7, 7.0),
        ):
            self.assertEqual(v.real(value, "value"), expected)
        self.assertEqual(copysign(1.0, v.real(-0.0, "zero")), -1.0)

    def test_malformed_scalar_values_fail_with_field_context(self):
        for value in (None, True, "1", [], {}, inf, -inf, float("nan"), 2**2000):
            with self.subTest(value=repr(value)[:70]):
                with self.assertRaises(DomainError) as raised:
                    v.real(value, "pressure")
                self.assertEqual(raised.exception.field, "pressure")

    def test_recognized_real_conversion_errors_are_typed(self):
        for error in (TypeError, ValueError, OverflowError):

            class BrokenFloat(float):
                def __float__(self, *, error=error):
                    raise error("conversion failed")

            with self.subTest(error=error):
                self.domain_failure(lambda: v.real(BrokenFloat(1.0), "flow"), "flow")

    def test_integer_keeps_large_exact_values_and_refuses_boolean_or_fraction(self):
        self.assertEqual(v.integer(2**2000, "seed"), 2**2000)
        for bad in (True, 1.0, Fraction(1), "1", None):
            with self.subTest(bad=bad), self.assertRaises(DomainError):
                v.integer(bad, "seed")

    def test_recognized_integer_conversion_errors_are_typed(self):
        for error in (TypeError, ValueError, OverflowError):

            class BrokenInteger(int):
                def __int__(self, *, error=error):
                    raise error("conversion failed")

            with self.subTest(error=error):
                self.domain_failure(lambda: v.integer(BrokenInteger(1), "count"), "count")

    def test_text_preserves_literal_unicode_but_rejects_invalid_encoding(self):
        text = "  ölçüm 🚀  "
        self.assertEqual(v.text(text, "name"), text)
        for bad in ("\ud800", "\udfff", None, "  ", True):
            with self.subTest(bad=repr(bad)), self.assertRaises(DomainError):
                v.text(bad, "name")

    def test_string_lists_validate_unicode_even_when_blanks_are_allowed(self):
        self.assertEqual(v.strings(["", "ölçüm"], "warnings", allow_blank=True), ("", "ölçüm"))
        for flag in (0, 1, None, "false"):
            with self.subTest(flag=flag), self.assertRaises(DomainError):
                v.strings([], "warnings", allow_blank=flag)
        with self.assertRaises(DomainError):
            v.strings(["\ud800"], "warnings", allow_blank=True)

    def test_unordered_containers_cannot_define_component_or_record_order(self):
        cases = (
            (v.vector3, {1.0, 2.0, 3.0}),
            (v.real_tuple, {1.0, 2.0}),
            (v.real_matrix, {(1.0, 2.0), (3.0, 4.0)}),
            (lambda x, label: v.records(x, int, label), {1, 2}),
            (v.strings, {"a", "b"}),
            (v.name_value_pairs, {("a", 1.0), ("b", 2.0)}),
        )
        for reader, value in cases:
            for bad in (value, frozenset(value)):
                with self.subTest(reader=reader, bad=bad), self.assertRaises(DomainError):
                    reader(bad, "data")

    def test_binary_buffers_are_not_implicit_numeric_sequences(self):
        readers = (
            v.vector3,
            v.real_tuple,
            v.real_matrix,
            v.strings,
            v.name_value_pairs,
            lambda x, label: v.records(x, int, label),
        )
        for bad in (b"abc", bytearray((1, 2, 3)), memoryview(b"abc")):
            for reader in readers:
                with (
                    self.subTest(reader=reader, bad=type(bad).__name__),
                    self.assertRaises(DomainError),
                ):
                    reader(bad, "data")

    def test_bad_iterator_protocol_does_not_escape_as_type_error(self):
        readers = (
            v.vector3,
            v.real_tuple,
            v.real_matrix,
            v.strings,
            v.name_value_pairs,
            lambda x, label: v.records(x, int, label),
        )
        for reader in readers:
            with self.subTest(reader=reader):
                self.domain_failure(lambda*, reader=reader: reader(BrokenIterator(), "data"), "data")

    def test_vector_cardinality_check_reads_at_most_four_items(self):
        seen = []

        def endless():
            while True:
                seen.append(1)
                yield 1.0

        with self.assertRaises(DomainError):
            v.vector3(endless(), "direction")
        self.assertEqual(len(seen), 4)

    def test_pair_cardinality_check_reads_at_most_three_items(self):
        seen = []

        def endless():
            for value in ("fuel", 1.0):
                seen.append(value)
                yield value
            while True:
                seen.append(2.0)
                yield 2.0

        with self.assertRaises(DomainError):
            v.name_value_pairs([endless()], "flow")
        self.assertEqual(len(seen), 3)

    def test_valid_ordered_iterators_and_legacy_sequences_keep_order(self):
        expected = (3.0, 1.0, 2.0)
        for data in (
            [3, 1, 2],
            iter((3, 1, 2)),
            LegacySequence((3, 1, 2)),
            array("d", (3.0, 1.0, 2.0)),
        ):
            self.assertEqual(v.vector3(data, "vector"), expected)
        self.assertEqual(v.real_tuple(LegacySequence((2, 1)), "data"), (2.0, 1.0))

    def test_nested_mapping_and_set_pairs_are_not_positional_data(self):
        for bad in (
            {0: "fuel", 1: 1.0},
            {"name": "fuel", "value": 1.0},
            {"fuel": 99.0, 1.0: 88.0},
            {"fuel", 1.0},
        ):
            with self.subTest(bad=bad), self.assertRaises(DomainError):
                v.name_value_pairs([bad], "flow")

    def test_name_value_pairs_keep_duplicates_and_snapshot_mutable_inputs(self):
        raw = [["fuel", 1], ["fuel", 2]]
        result = v.name_value_pairs(raw, "flow")
        raw[0][1] = 99
        self.assertEqual(result, (("fuel", 1.0), ("fuel", 2.0)))
        self.assertEqual(v.name_value_pairs([iter(("fuel", 1))], "flow"), (("fuel", 1.0),))

    def test_pair_fields_and_cardinality_are_checked_without_truncation(self):
        for bad in (["fuel"], ["fuel", 1, 2], ["fuel", True], [" ", 1], ["\ud800", 1]):
            with self.subTest(bad=repr(bad)), self.assertRaises(DomainError):
                v.name_value_pairs([bad], "flow")

    def test_matrix_is_rectangular_and_detached_from_mutable_rows(self):
        raw = [[1, 2], [3, 4]]
        matrix = v.real_matrix(raw, "matrix")
        raw[0][0] = 99
        self.assertEqual(matrix, ((1.0, 2.0), (3.0, 4.0)))
        self.assertEqual(v.real_matrix([], "matrix"), ())
        self.assertEqual(v.real_matrix([[], []], "matrix"), ((), ()))
        with self.assertRaises(DomainError):
            v.real_matrix([[1, 2], [3]], "matrix")

    def test_mapping_values_and_keys_use_the_same_scalar_and_text_contract(self):
        raw = {"flow": 1}
        result = v.real_mapping(raw, "values")
        raw["flow"] = 99
        self.assertEqual(result, {"flow": 1.0})
        self.assertEqual(v.real_mapping(None, "values"), {})
        for bad in ({"flow": True}, {"\ud800": 1}, {"flow": Fraction(nextafter(0.0, 1.0)) / 4}, []):
            with self.subTest(bad=repr(bad)), self.assertRaises(DomainError):
                v.real_mapping(bad, "values")

    def test_records_snapshot_and_preserve_explicit_allowed_type_contract(self):
        raw = [Record(1.0), Record(2.0)]
        result = v.records(raw, Record, "rows")
        raw.clear()
        self.assertEqual(result, (Record(1.0), Record(2.0)))
        self.assertEqual(v.records([1, "a"], (int, str), "rows"), (1, "a"))
        with self.assertRaises(DomainError):
            v.records([Record(1.0), "wrong"], Record, "rows")
        for bad_kind in ([], "Record", ()):
            with self.subTest(bad_kind=bad_kind), self.assertRaises(TypeError):
                v.records([], bad_kind, "rows")

    def test_optional_scalar_annotations_are_validated_without_evaluation(self):
        for annotation in (
            Optional[float],  # noqa: UP045 - Legacy annotations are deliberate validation fixtures.
            Union[None, float],  # noqa: UP007 - Preserve typing.Union coverage beside modern syntax.
            "typing.Optional[float]",
            "float | None",
            "'float'",
            "typing.Union[float, None]",
        ):
            cls = make_dataclass("Scalar", [("value", annotation)], frozen=True)
            bad = cls(True)
            with self.subTest(annotation=annotation), self.assertRaises(DomainError):
                v.coerce_fields(bad, "record")
            good = cls(2)
            v.coerce_fields(good, "record")
            self.assertIs(type(good.value), float)

    def test_failed_dataclass_normalization_does_not_mutate_earlier_fields(self):
        cls = make_dataclass("Pair", [("first", float), ("second", float)], frozen=True)
        record = cls(1, "bad")
        with self.assertRaises(DomainError):
            v.coerce_fields(record, "record")
        self.assertIs(type(record.first), int)
        self.assertEqual(record.second, "bad")

    def test_dataclass_class_is_not_an_instance_and_opaque_fields_are_untouched(self):
        for bad in (None, object(), Record):
            with self.subTest(bad=bad), self.assertRaises(TypeError):
                v.coerce_fields(bad, "record")

        @dataclass(frozen=True)
        class Deferred:
            value: float
            opaque: object = field(init=False)

        record = Deferred(1)
        v.coerce_fields(record, "record")
        self.assertIs(type(record.value), float)
        self.assertFalse(hasattr(record, "opaque"))
        ANNOTATION_CALLS.clear()
        opaque = OpaqueAnnotation("not a scalar")
        v.coerce_fields(opaque, "record")
        self.assertEqual(ANNOTATION_CALLS, [])
        self.assertEqual(opaque.value, "not a scalar")

    def test_unrelated_callback_failures_are_not_misreported_as_bad_data(self):
        class FailsUnexpectedly(float):
            def __float__(self):
                raise RuntimeError("application failure")

        with self.assertRaisesRegex(RuntimeError, "application failure"):
            v.real(FailsUnexpectedly(1.0), "value")
