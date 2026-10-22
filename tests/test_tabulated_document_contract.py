"""Reject ambiguous or silently discarded trace data before model consumption."""

import hashlib
import json
import unittest
from copy import deepcopy
from dataclasses import replace

from rocket_propulsion.core.errors import DomainError, InputError
from rocket_propulsion.core.validation import read_json_document
from rocket_propulsion.propulsion.burns import (
    TabulatedBurnArtifact,
    TabulatedBurnSegment,
    tabulated_artifact_from_csv,
)

REFUSED = (DomainError, InputError)
HEADER = "time_s,delivered_thrust_n,total_mass_flow_kg_s"


def artifact():
    segment = TabulatedBurnSegment(
        0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, (("fuel", 1.0),), (("fuel", 1.0),), "trace"
    )
    return TabulatedBurnArtifact((segment,), "raw", "a" * 64, reference_ids=("ref",))


def csv_model(source, policy="linear"):
    return tabulated_artifact_from_csv(source, source_id="trace", interpolation=policy)


class TabulatedDocumentAttacks(unittest.TestCase):
    def test_nonzero_decimal_underflow_cannot_become_an_accepted_zero(self):
        for token in ("1e-9999", "-1e-9999"):
            with self.subTest(token=token), self.assertRaises(InputError):
                read_json_document('{"value":' + token + "}", "test", lambda data: data)
        self.assertEqual(
            read_json_document('{"value":0e-999999}', "test", lambda data: data), {"value": 0.0}
        )

    def test_unknown_root_fields_cannot_disappear_behind_model_hash(self):
        payload = artifact().to_dict()
        for name, value in (("units", "lbf"), ("approved", True), ("zz_unknown", 1)):
            bad = {**payload, name: value}
            with self.subTest(name=name), self.assertRaises(REFUSED):
                TabulatedBurnArtifact.from_dict(bad)

    def test_unknown_segment_fields_cannot_disappear(self):
        bad = artifact().to_dict()
        bad["segments"][0]["flow_units"] = "lbm/s"
        with self.assertRaises(REFUSED):
            TabulatedBurnArtifact.from_dict(bad)

    def test_summary_booleans_are_not_numerical_agreement(self):
        for field in (
            "duration_s",
            "delivered_total_impulse_n_s",
            "axial_total_impulse_n_s",
            "consumed_propellant_kg",
        ):
            bad = artifact().to_dict()
            self.assertEqual(bad["summary"][field], 1.0)
            bad["summary"][field] = True
            with self.subTest(field=field), self.assertRaises(REFUSED):
                TabulatedBurnArtifact.from_dict(bad)

    def test_raw_arrays_cannot_be_python_only_containers(self):
        payload = artifact().to_dict()
        for key in ("segments", "warnings", "reference_ids"):
            bad = deepcopy(payload)
            bad[key] = tuple(bad[key])
            with self.subTest(key=key), self.assertRaises(REFUSED):
                TabulatedBurnArtifact.from_dict(bad)

    def test_flow_pairs_cannot_be_indexed_mappings_or_tuples(self):
        for pair in ({0: "fuel", 1: 1.0}, ("fuel", 1.0)):
            bad = artifact().to_dict()
            bad["segments"][0]["start_stream_mass_flows_kg_s"] = [pair]
            with self.subTest(pair=pair), self.assertRaises(REFUSED):
                TabulatedBurnArtifact.from_dict(bad)

    def test_required_metadata_cannot_be_recreated_by_defaults(self):
        bad = artifact().to_dict()
        bad.pop("warnings")
        with self.assertRaises(REFUSED):
            TabulatedBurnArtifact.from_dict(bad)

    def test_duplicate_json_members_are_rejected_at_every_depth(self):
        source = json.dumps(artifact().to_dict(), separators=(",", ":"))
        for old, new in (
            ('"source_id":"raw"', '"source_id":"fake","source_id":"raw"'),
            ('"phase":"trace"', '"phase":"fake","phase":"trace"'),
        ):
            bad = source.replace(old, new)
            self.assertNotEqual(source, bad)
            with self.subTest(member=old), self.assertRaises(REFUSED):
                TabulatedBurnArtifact.from_json(bad)

    def test_shared_json_parser_rejects_duplicates_before_builder(self):
        calls = []
        for source in ('{"a":1,"a":1}', '{"nested":{"a":0,"a":1}}'):
            with self.subTest(source=source), self.assertRaises(InputError):
                read_json_document(source, "test", lambda data: calls.append(data))
        self.assertEqual(calls, [])

    def test_nonfinite_json_tokens_and_overflow_literals_never_reach_builder(self):
        calls = []
        for token in ("NaN", "Infinity", "-Infinity", "1e309", "-1e309"):
            for source in (f'{{"ignored":{token}}}', f'{{"nested":[{token}]}}'):
                with self.subTest(source=source), self.assertRaises(InputError):
                    read_json_document(source, "test", lambda data: calls.append(data))
        self.assertEqual(calls, [])

    def test_unencodable_json_keys_and_values_never_reach_builder(self):
        for source in ('{"bad":"\\ud800"}', '{"\\udfff":1}', '{"nested":["\\ud800"]}'):
            with self.subTest(source=source), self.assertRaises(InputError):
                read_json_document(source, "test", lambda data: data)

    def test_model_metadata_must_be_utf8_encodable_before_hashing(self):
        value = artifact()
        for change in (
            {"source_id": "\ud800"},
            {"warnings": ("\ud800",)},
            {"reference_ids": ("\ud800",)},
        ):
            with self.subTest(change=change), self.assertRaises(DomainError):
                replace(value, **change)
        with self.assertRaises(DomainError):
            replace(value.segments[0], phase="\ud800")
        with self.assertRaises(DomainError):
            replace(value.segments[0], start_stream_mass_flows_kg_s=(("\ud800", 1.0),))

    def test_valid_unicode_and_literal_source_digest_roundtrip(self):
        value = replace(
            artifact(), source_id="ölçüm 🚀", source_sha256="A" * 64, warnings=("ölçüm koşulları",)
        )
        opened = TabulatedBurnArtifact.from_json(value.to_json())
        self.assertEqual(opened, value)
        self.assertEqual(opened.to_json(), value.to_json())
        self.assertNotEqual(
            value.artifact_hash, replace(value, source_sha256="a" * 64).artifact_hash
        )

    def test_missing_or_corrupt_summary_is_rejected(self):
        for replacement in (None, {}, {"duration_s": 1.0}, "summary"):
            bad = artifact().to_dict()
            bad["summary"] = replacement
            with self.subTest(replacement=replacement), self.assertRaises(REFUSED):
                TabulatedBurnArtifact.from_dict(bad)

    def test_optional_empty_reference_list_has_defined_semantics(self):
        value = replace(artifact(), reference_ids=())
        payload = value.to_dict()
        payload["reference_ids"] = []
        self.assertEqual(TabulatedBurnArtifact.from_dict(payload), value)

    def test_csv_duplicate_header_names_are_never_last_wins(self):
        for extra in ("time_s", "delivered_thrust_n", "stream:fuel"):
            header = HEADER + (",stream:fuel" if extra == "stream:fuel" else "") + "," + extra
            width = len(header.split(","))
            first_extra = ["0"] if extra == "time_s" else ["1"] * (width - 3)
            source = (
                header
                + "\n"
                + ",".join(["0", "1", "1"] + first_extra)
                + "\n"
                + ",".join(["1", "1", "1"] + ["1"] * (width - 3))
                + "\n"
            )
            with self.subTest(extra=extra), self.assertRaises(REFUSED):
                csv_model(source)

    def test_csv_unknown_columns_cannot_silently_drop_physics(self):
        for extra in ("thrust_units", "total_mass_flow_lb_s", "pressure_pa"):
            source = HEADER + "," + extra + "\n0,1,1,1\n1,1,1,1\n"
            with self.subTest(extra=extra), self.assertRaises(REFUSED):
                csv_model(source)

    def test_csv_each_record_must_match_the_declared_width(self):
        for source in (
            HEADER + "\n0,1,1,extra\n1,1,1\n",
            HEADER + ",phase\n0,1,1,a\n1,1,1\n",
            HEADER + "\n0,1,1\n1,1,1,\n",
        ):
            with self.subTest(source=source), self.assertRaises(REFUSED):
                csv_model(source)

    def test_csv_declared_optional_numeric_cells_cannot_default(self):
        for field in ("axial_thrust_n", "stream:fuel"):
            source = HEADER + "," + field + "\n0,1,1,\n1,1,1,1\n"
            with self.subTest(field=field), self.assertRaises(REFUSED):
                csv_model(source)

    def test_csv_hold_terminal_sample_is_still_validated(self):
        cases = (
            ("", "1,-1"),
            ("", "1,nan"),
            (",axial_thrust_n", "1,1,2"),
            (",stream:fuel", "1,1,0.5"),
            (",phase", "1,1,"),
        )
        for suffix, tail in cases:
            first = "0,1,1" + (
                ",1" if suffix and suffix != ",phase" else ",start" if suffix else ""
            )
            source = HEADER + suffix + "\n" + first + "\n1," + tail + "\n"
            with self.subTest(source=source), self.assertRaises(REFUSED):
                csv_model(source, "hold")

    def test_csv_nonfinite_and_negative_physics_are_refused(self):
        for row in ("0,nan,1", "0,1,inf", "0,-1,1", "0,1,-1", "0,1e-9999,1", "0,1,1e-9999"):
            with self.subTest(row=row), self.assertRaises(REFUSED):
                csv_model(HEADER + "\n" + row + "\n1,1,1\n")

    def test_csv_defaults_only_apply_when_columns_are_absent(self):
        value = csv_model(HEADER + "\n0,1,.25\n1,3,.75\n2,2,.5\n")
        self.assertEqual(value.axial_total_impulse_n_s, 4.5)
        self.assertEqual(value.delivered_total_impulse_n_s, 4.5)
        self.assertEqual(value.consumed_propellant_kg, 1.125)
        self.assertEqual(value.segments[0].phase, "imported-trace")
        self.assertEqual(value.segments[0].start_stream_mass_flows_kg_s, ())
        held = csv_model(HEADER + "\n0,1,.25\n1,3,.75\n2,2,.5\n", "hold")
        self.assertEqual(held.delivered_total_impulse_n_s, 4.0)
        self.assertEqual(held.consumed_propellant_kg, 1.0)

    def test_csv_source_hash_binds_literal_bytes_and_quoted_phase(self):
        source = HEADER + ',phase\n0,1,1,"başlangıç, faz"\n1,1,1,bitiş\n'
        a, b = csv_model(source), csv_model(source.replace("\n", "\r\n"))
        self.assertEqual(a.source_sha256, hashlib.sha256(source.encode("utf-8")).hexdigest())
        self.assertNotEqual(a.source_sha256, b.source_sha256)
        self.assertEqual(a.segments, b.segments)
        self.assertEqual(a.segments[0].phase, "başlangıç, faz")

    def test_csv_parser_and_encoding_errors_are_typed(self):
        for source in (
            None,
            b"text",
            7,
            HEADER + ',phase\n0,1,1,"start"junk\n1,1,1,end\n',
            HEADER + ",phase\n0,1,1,\ud800\n1,1,1,end\n",
        ):
            with self.subTest(source=repr(source)[:100]), self.assertRaises(REFUSED):
                csv_model(source)

    def test_reader_snapshot_is_detached_from_mutable_raw_document(self):
        payload = artifact().to_dict()
        opened = TabulatedBurnArtifact.from_dict(payload)
        before = opened.to_json()
        payload["segments"][0]["start_stream_mass_flows_kg_s"][0][1] = 99.0
        payload["reference_ids"].append("fake")
        self.assertEqual(opened.to_json(), before)
