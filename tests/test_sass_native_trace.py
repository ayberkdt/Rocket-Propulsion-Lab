"""Exercise native CSV channel and row binding through the normal validation path."""

from dataclasses import replace
from fractions import Fraction
from hashlib import sha256
from math import inf, nextafter
from unittest.mock import patch

import pytest

from rocket_propulsion.core.errors import DomainError
from rocket_propulsion.studies.sass import validation
from rocket_propulsion.studies.sass.native_trace import NativeCsvTrace, native_csv_trace_status
from rocket_propulsion.studies.sass.validation import (
    NativeAcquisitionFile,
    acquisition_partition_sha256,
)
from tests.test_sass_acquisition_partition import setup_case


def approve_plan(record, plan):
    return replace(plan, reviewed_acquisition_partition_sha256=acquisition_partition_sha256(record, plan))


def csv_case(tmp_path, *, header="tick,pressure,spare", value="10", errors=None):
    record, plan, model, run = setup_case(tmp_path)
    condition = plan.targets[0].condition_record
    origin = 2**55
    rows = [(origin - 4, "12", "10"), (origin, value, "11"), (origin + 4, "10", "11"), (origin + 8, "13", "10")]
    suffixes = ("", "", "", "") if errors is None else tuple(f",{error}" for error in errors)
    if errors is not None:
        header += ",error"
    csv = (header + "\n" + "".join(f"{tick},{center},{spare}{suffix}\n"
                                   for (tick, center, spare), suffix in zip(rows, suffixes, strict=True))).encode()
    path = tmp_path / "native.csv"
    path.write_bytes(csv)
    digest = sha256(csv).hexdigest()
    partition = plan.acquisition_partition
    original = next(row for row in partition.slices if row.condition_record == condition)
    scopes = tuple(replace(row, native_sha256=digest, past_end_index=4) if row is original else row
                   for row in partition.slices)
    sources = tuple(source for source in plan.native_acquisitions if source.native_sha256 != original.native_sha256)
    mapping = NativeCsvTrace(condition, run.channels[0].trace.channel_id, digest, "pressure", "Pa",
                             1, 3, -origin, "reviewed original pressure column and time origin")
    plan = replace(plan, acquisition_partition=replace(partition, slices=scopes, native_csv_traces=(mapping,)),
                   native_acquisitions=(*sources, NativeAcquisitionFile(digest, path)))
    return record, approve_plan(record, plan), model, run, path, mapping


def with_mapping(record, plan, mapping):
    updated = replace(plan, acquisition_partition=replace(plan.acquisition_partition, native_csv_traces=(mapping,)))
    return approve_plan(record, updated)


def test_wrong_column_or_rows_fail_despite_matching_file_and_reviewed_digest(tmp_path):
    record, plan, model, run, _, mapping = csv_case(tmp_path)
    valid = validation.assess_validation(record, plan=plan, model_run=model, transient_run=run)
    assert valid["comparison_status"] == "pass"
    assert valid["native_trace_review"]["status"] == "verified"
    for edits, reason in (
        ({"value_column": "spare"}, "native_csv_sample_mismatch"),
        ({"first_index": 0, "past_end_index": 2}, "native_csv_sample_mismatch"),
        ({"tick_offset": mapping.tick_offset + 1}, "native_csv_sample_mismatch"),
        ({"unit": "N"}, "native_csv_unit_mismatch"),
        ({"first_index": 3, "past_end_index": 5}, "native_csv_outside_reviewed_slice"),
        ({"first_index": 0}, "native_csv_sample_count_mismatch"),
    ):
        attacked = with_mapping(record, plan, replace(mapping, **edits))
        report = validation.assess_validation(record, plan=attacked, model_run=model, transient_run=run)
        assert report["acquisition_partition_review"]["status"] == "verified"
        assert report["transient_replay_review"]["status"] == "verified"
        assert report["comparison_status"] == "unassessed"
        assert report["samples"][0]["reason"] == reason
        assert report["physical_acceptance_review"]["status"] == "unassessed"


def test_wrong_header_and_nonfinite_sample_cannot_enter_trace(tmp_path):
    for kwargs, reason in (
        ({"header": "tick,pressure,pressure"}, "native_csv_invalid_header"),
        ({"header": "time,pressure,spare"}, "native_csv_invalid_header"),
        ({"value": "1e999"}, "native_csv_invalid_sample"),
        ({"value": "NaN"}, "native_csv_invalid_sample"),
        ({"value": "10.000000000000002"}, "native_csv_sample_mismatch"),
    ):
        record, plan, model, run, _, _ = csv_case(tmp_path, **kwargs)
        report = validation.assess_validation(record, plan=plan, model_run=model, transient_run=run)
        assert report["comparison_status"] == "unassessed"
        assert report["samples"][0]["reason"] == reason


def test_csv_bytes_must_match_while_parsing_even_after_initial_hash_check(tmp_path):
    record, plan, model, run, path, _ = csv_case(tmp_path)
    original_hash = validation._file_sha256

    def swap_after_hash(file, label):
        digest = original_hash(file, label)
        if file == path:
            # An unselected record changes. The selected pressure centers stay identical.
            path.write_bytes(path.read_bytes().replace(b",13,10", b",14,10"))
        return digest

    with patch.object(validation, "_file_sha256", side_effect=swap_after_hash):
        result = validation.assess_validation(record, plan=plan, model_run=model, transient_run=run)
    assert result["acquisition_partition_review"]["native_acquisition_review"]["status"] == "verified"
    assert result["comparison_status"] == "unassessed"
    assert result["samples"][0]["reason"] == "native_csv_digest_mismatch"


def test_csv_mapping_cannot_reuse_old_approval_or_move_to_another_channel(tmp_path):
    record, plan, model, run, _, mapping = csv_case(tmp_path)
    unsigned = replace(plan, acquisition_partition=replace(
        plan.acquisition_partition, native_csv_traces=(replace(mapping, value_column="spare"),)
    ))
    result = validation.assess_validation(record, plan=unsigned, model_run=model, transient_run=run)
    assert result["samples"][0].get("reason") == "native_csv_mapping_unapproved"
    attacked = with_mapping(record, plan, replace(mapping, channel_id="some other channel"))
    result = validation.assess_validation(record, plan=attacked, model_run=model, transient_run=run)
    assert result["samples"][0]["reason"] == "native_csv_mapping_missing"
    assert result["native_trace_review"]["unplanned_channels"] == [(mapping.condition_record, "some other channel")]


def test_narrowed_native_error_withholds_measurement_acceptance(tmp_path):
    record, plan, model, run, _, mapping = csv_case(tmp_path, errors=("0", "0", "0.001", "0"))
    declared = replace(mapping, absolute_error_column="error")
    unsigned = replace(plan, acquisition_partition=replace(plan.acquisition_partition, native_csv_traces=(declared,)))
    rejected = validation.assess_validation(record, plan=unsigned, model_run=model, transient_run=run)
    assert rejected["samples"][0]["reason"] == "native_csv_mapping_unapproved"

    narrowed = with_mapping(record, plan, declared)
    result = validation.assess_validation(record, plan=narrowed, model_run=model, transient_run=run)
    assert result["acquisition_partition_review"]["status"] == "verified"
    assert result["transient_replay_review"]["status"] == "verified"
    assert result["comparison_status"] == "unassessed"
    assert result["physical_acceptance_review"]["status"] == "unassessed"
    assert result["samples"][0]["reason"] == "native_csv_uncertainty_narrowed"

    record, plan, model, run, _, mapping = csv_case(tmp_path, errors=("invalid", "0", "0.0", "invalid"))
    approved = with_mapping(record, plan, replace(mapping, absolute_error_column="error"))
    result = validation.assess_validation(record, plan=approved, model_run=model, transient_run=run)
    assert result["comparison_status"] == "pass"
    assert result["physical_acceptance_review"]["status"] == "pass"
    assert result["native_trace_review"]["status"] == "verified"


@pytest.mark.parametrize(("center", "radius", "bound", "accepted"), [
    ("1.00000000000000001", "0.125", 0.125, False),
    ("0.99999999999999999", "0.125", 0.125, False),
    ("1.00000000000000001", "0.125", nextafter(0.125, inf), True),
    ("1", "0.12500000000000001", 0.125, False),
    ("1", "0.12500000000000001", nextafter(0.125, inf), True),
    ("1", "0.125", 0.125, True),
    ("1", "1e-325", 0.0, False),
    ("1", "1e-325", nextafter(0.0, inf), True),
])
def test_native_decimal_interval_must_fit_both_binary64_endpoints(tmp_path, center, radius, bound, accepted):
    _, _, _, run, path, mapping = csv_case(tmp_path, value=center, errors=("0", radius, "0", "0"))
    trace = run.channels[0].trace
    sample = replace(trace.samples[0], value=float(center), absolute_error=bound)
    trace = replace(trace, samples=(sample, trace.samples[1]))
    x, u, c, r = Fraction(center), Fraction(radius), Fraction(sample.value), Fraction(bound)
    assert (c - r <= x - u and x + u <= c + r) is accepted
    expected = "verified" if accepted else "native_csv_uncertainty_narrowed"
    assert native_csv_trace_status(path, replace(mapping, absolute_error_column="error"), trace) == expected


@pytest.mark.parametrize("radius", ["-1e-325", "1e-999999999", "1e+999999999", "inf", "", "0_1", " 0"])
def test_malformed_or_unbounded_native_radius_fails_closed(tmp_path, radius):
    _, _, _, run, path, mapping = csv_case(tmp_path, errors=("0", radius, "0", "0"))
    scoped = replace(mapping, absolute_error_column="error")
    assert native_csv_trace_status(path, scoped, run.channels[0].trace) == "native_csv_invalid_absolute_error"


def test_uncertainty_column_cannot_alias_value_or_clock_or_go_missing(tmp_path):
    _, _, _, run, path, mapping = csv_case(tmp_path)
    trace = run.channels[0].trace
    with pytest.raises(DomainError, match="own column"):
        replace(mapping, absolute_error_column="tick")
    with pytest.raises(DomainError, match="own column"):
        replace(mapping, absolute_error_column=mapping.value_column)
    assert native_csv_trace_status(path, replace(mapping, absolute_error_column="error"), trace) == (
        "native_csv_invalid_header"
    )


def test_native_clock_scale_requires_approval_and_matches_pipeline_clock(tmp_path):
    record, plan, model, run, _, mapping = csv_case(tmp_path)
    changed_period = nextafter(1.0, inf)
    wrong_scope = replace(mapping, tick_period_s=changed_period)
    unsigned = replace(plan, acquisition_partition=replace(plan.acquisition_partition, native_csv_traces=(wrong_scope,)))
    result = validation.assess_validation(record, plan=unsigned, model_run=model, transient_run=run)
    assert result["samples"][0]["reason"] == "native_csv_mapping_unapproved"

    approved = with_mapping(record, plan, wrong_scope)
    result = validation.assess_validation(record, plan=approved, model_run=model, transient_run=run)
    assert result["transient_replay_review"]["status"] == "verified"
    assert result["acquisition_partition_review"]["status"] == "verified"
    assert result["samples"][0].get("reason") == "native_csv_tick_period_mismatch"
    assert result["comparison_status"] == "unassessed"
    assert result["physical_acceptance_review"]["status"] == "unassessed"

    correct = with_mapping(record, plan, replace(mapping, tick_period_s=1.0))
    result = validation.assess_validation(record, plan=correct, model_run=model, transient_run=run)
    assert result["native_trace_review"]["status"] == "verified"
    assert result["comparison_status"] == "pass"
    assert result["physical_acceptance_review"]["status"] == "pass"


@pytest.mark.parametrize("clock", [0.0, -1.0, inf, float("nan"), True, "1"])
def test_native_clock_period_must_be_positive_finite_numeric(tmp_path, clock):
    *_, mapping = csv_case(tmp_path)
    with pytest.raises(DomainError):
        replace(mapping, tick_period_s=clock)
