"""Challenge calibration/holdout separation independently of channel and sample names."""

from dataclasses import replace
from hashlib import sha256

import pytest

from rocket_propulsion.core.errors import DomainError
from rocket_propulsion.studies.sass.validation import (
    AcquisitionPartition,
    NativeAcquisitionFile,
    NativeAcquisitionSlice,
    acquisition_partition_sha256,
    assess_validation,
)
from tests.test_sass_physical_acceptance import approve, limits
from tests.test_sass_transient_validation import (
    RECIPE,
    archived,
    pressure_row,
    synthetic_native_files,
    with_fixed_prediction,
)


def native(condition, firing, file, start, end):
    return NativeAcquisitionSlice(condition, firing, sha256(file.encode()).hexdigest(), start, end)


def reviewed(record, plan, *slices, calibration_samples=None):
    calibration, = plan.calibration_condition_records
    if calibration_samples is None:
        calibration_samples = tuple((identifier, calibration) for identifier in record.calibration_sample_ids)
    names = ("DAQ A", "DAQ B", "DAQ C")
    expected = {item.native_sha256 for item in slices}
    sources = synthetic_native_files(record, (name for name in names if sha256(name.encode()).hexdigest() in expected))
    plan = replace(plan, acquisition_partition=AcquisitionPartition("independent DAQ mapping", slices,
                                                                   calibration_samples), native_acquisitions=sources)
    return replace(plan, reviewed_acquisition_partition_sha256=acquisition_partition_sha256(record, plan))


def setup_case(tmp_path):
    run = archived()
    row = pressure_row()
    record, plan, model, _ = with_fixed_prediction(tmp_path, run, (row, RECIPE))
    target = approve(record, plan.targets[0], row, limits(high=14))
    return record, replace(plan, targets=(target,)), model, run


def test_new_names_and_reapproval_cannot_hide_shared_firing_or_native_indices(tmp_path):
    record, plan, model, run = setup_case(tmp_path)
    calibration, = plan.calibration_condition_records
    holdout = plan.targets[0].condition_record
    start = 2**58
    used = native(calibration, "firing A", "DAQ A", start, start + 8)
    for replayed, reason in (
        (native(holdout, "firing A", "DAQ B", 0, 4), "shared_firing"),
        (native(holdout, "renamed firing", "DAQ A", start + 7, start + 12), "overlapping_native_slice"),
        (native(holdout, "renamed firing", "DAQ A", start - 1, start + 1), "overlapping_native_slice"),
    ):
        approved = reviewed(record, plan, used, replayed)
        report = assess_validation(record, plan=approved, model_run=model, transient_run=run)
        assert report["comparison_status"] == "unassessed"
        assert report["samples"][0].get("reason") == "calibration_acquisition_overlap"
        assert report["acquisition_partition_review"]["conflicts"] == [{
            "holdout_condition": holdout, "calibration_condition": calibration, "reasons": [reason],
        }]
        assert report["physical_acceptance_review"]["status"] == "pass"


def test_adjacent_native_slices_can_be_independent_if_the_firings_differ(tmp_path):
    record, plan, model, run = setup_case(tmp_path)
    calibration, = plan.calibration_condition_records
    holdout = plan.targets[0].condition_record
    start = 2**58
    used = native(calibration, "firing A", "DAQ A", start, start + 8)
    following = native(holdout, "firing B", "DAQ A", start + 8, start + 16)
    secondary = native(holdout, "firing B", "DAQ B", 0, 4)
    approved = reviewed(record, plan, used, following, secondary)
    reordered = reviewed(record, plan, secondary, following, used)
    assert reordered == approved
    report = assess_validation(record, plan=approved, model_run=model, transient_run=run)
    assert report["comparison_status"] == "pass"
    assert report["acquisition_partition_review"]["status"] == "verified"


def test_model_bytes_and_sample_review_cannot_replace_partition_approval(tmp_path):
    record, plan, model, run = setup_case(tmp_path)
    calibration, = plan.calibration_condition_records
    holdout = plan.targets[0].condition_record
    used = native(calibration, "firing A", "DAQ A", 0, 8)
    following = native(holdout, "firing B", "DAQ B", 0, 8)
    approved = reviewed(record, plan, used, following)
    for attacked, reason in (
        (replace(approved, reviewed_acquisition_partition_sha256=None), "missing_reviewed_acquisition_partition"),
        (replace(approved, acquisition_partition=None), "missing_reviewed_acquisition_partition"),
        (replace(approved, acquisition_partition=replace(approved.acquisition_partition,
                                                         slices=(used, replace(following, past_end_index=7)))),
         "acquisition_partition_sha256_mismatch"),
        (replace(approved, calibration_condition_records=("new calibration label",)),
         "acquisition_partition_sha256_mismatch"),
        (reviewed(record, plan, following), "incomplete_acquisition_partition"),
        (reviewed(record, plan, used, following, calibration_samples=()), "incomplete_acquisition_partition"),
        (reviewed(record, plan, used, following, calibration_samples=(("renamed calibration ID", calibration),)),
         "incomplete_acquisition_partition"),
        (reviewed(record, plan, used, following,
                  calibration_samples=((record.calibration_sample_ids[0], holdout),)), "incomplete_acquisition_partition"),
        (reviewed(record, plan, used, following, native("unplanned", "firing C", "DAQ C", 0, 8)),
         "incomplete_acquisition_partition"),
    ):
        result = assess_validation(record, plan=attacked, model_run=model, transient_run=run)
        assert result["model_run_review"]["status"] == "verified"
        assert result["transient_replay_review"]["status"] == "verified"
        assert result["samples"][0].get("reason") == reason
        assert result["comparison_status"] == "unassessed"
        assert result["physical_acceptance_review"]["status"] == "pass"


def test_native_slice_uses_exact_half_open_indices_and_sha256():
    valid = native("condition", "firing", "native", 0, 8)
    for changes in (
        {"first_index": True}, {"past_end_index": 8.0}, {"first_index": -1}, {"past_end_index": 0},
        {"native_sha256": "A" * 64}, {"firing_id": " "},
    ):
        with pytest.raises(DomainError):
            replace(valid, **changes)
    with pytest.raises(DomainError, match="unique native slices"):
        AcquisitionPartition("review", (valid, valid))
    with pytest.raises(DomainError, match="unique acquisition mapping"):
        AcquisitionPartition("review", (valid,), (("same ID", "first firing"), ("same ID", "second firing")))


def test_partition_approval_binds_raw_model_and_calibration_identifiers(tmp_path):
    record, plan, _, _ = setup_case(tmp_path)
    digest = acquisition_partition_sha256(record, plan)
    for changes in (
        {"raw_data_sha256": "0" * 64}, {"model_sha256": "0" * 64},
        {"model_configuration_id": "edited model configuration"},
        {"calibration_sample_ids": ("different calibration sample",)},
        {"independence_record": "different review"},
    ):
        assert acquisition_partition_sha256(replace(record, **changes), plan) != digest


def test_complete_native_file_content_is_checked_after_relocation_and_after_edit(tmp_path):
    record, plan, model, run = setup_case(tmp_path)
    scope = plan.acquisition_partition.slices[0]
    original = b"A" * 1_048_576 + b"last DAQ frame"
    path = tmp_path / "original.daq"
    path.write_bytes(original)
    native_digest = sha256(original).hexdigest()
    scopes = tuple(replace(row, native_sha256=native_digest) if row is scope else row
                   for row in plan.acquisition_partition.slices)
    sources = tuple(source for source in plan.native_acquisitions if source.native_sha256 != scope.native_sha256)
    updated = replace(plan, acquisition_partition=replace(plan.acquisition_partition, slices=scopes),
                      native_acquisitions=(*sources, NativeAcquisitionFile(native_digest, path)))
    approved = replace(updated, reviewed_acquisition_partition_sha256=acquisition_partition_sha256(record, updated))
    relocated = tmp_path / "relocated.daq"
    path.rename(relocated)
    moved = replace(approved, native_acquisitions=(*sources, NativeAcquisitionFile(native_digest, relocated)))
    assert acquisition_partition_sha256(record, moved) == approved.reviewed_acquisition_partition_sha256
    valid = assess_validation(record, plan=moved, model_run=model, transient_run=run)
    assert valid["comparison_status"] == "pass"
    assert native_digest in valid["acquisition_partition_review"]["native_acquisition_review"]["verified_sha256"]
    changed = original[:-1] + b"!"
    relocated.write_bytes(changed)
    invalid = assess_validation(record, plan=moved, model_run=model, transient_run=run)
    assert invalid["comparison_status"] == "unassessed"
    assert invalid["samples"][0].get("reason") == "native_acquisition_archive_unverified"
    check = invalid["acquisition_partition_review"]["native_acquisition_review"]
    assert check["sha256_mismatches"] == {native_digest: sha256(changed).hexdigest()}
    assert invalid["physical_acceptance_review"]["status"] == "pass"


def test_missing_unreadable_and_unplanned_sources_withhold_model_agreement(tmp_path):
    record, plan, model, run = setup_case(tmp_path)
    selected, other = plan.native_acquisitions
    extra, = synthetic_native_files(record, ("unexpected DAQ file",))
    for sources, key, digest in (
        ((other,), "missing_sha256", selected.native_sha256),
        ((other, replace(selected, path=tmp_path / "missing.daq")), "unreadable_sha256", selected.native_sha256),
        ((other, replace(selected, path=tmp_path)), "unreadable_sha256", selected.native_sha256),
        ((selected, other, extra), "unplanned_sha256", extra.native_sha256),
    ):
        report = assess_validation(record, plan=replace(plan, native_acquisitions=sources),
                                   model_run=model, transient_run=run)
        assert report["comparison_status"] == "unassessed"
        assert report["samples"][0].get("reason") == "native_acquisition_archive_unverified"
        assert report["acquisition_partition_review"]["native_acquisition_review"][key] == [digest]
    with pytest.raises(DomainError, match="unique file binding"):
        replace(plan, native_acquisitions=(selected, other, selected))


def test_model_exports_and_empty_file_cannot_stand_in_for_original_daq(tmp_path):
    record, plan, model, run = setup_case(tmp_path)
    scope = plan.acquisition_partition.slices[0]
    empty = tmp_path / "empty.daq"
    empty.write_bytes(b"")
    sources = tuple(source for source in plan.native_acquisitions if source.native_sha256 != scope.native_sha256)
    for path, digest, reason in (
        (record.raw_data_path, record.raw_data_sha256, "repurposed_sha256"),
        (model.model_path, record.model_sha256, "repurposed_sha256"),
        (model.predictions_path, plan.reviewed_predictions_sha256, "repurposed_sha256"),
        (empty, sha256(b"").hexdigest(), "empty_sha256"),
    ):
        scopes = tuple(replace(row, native_sha256=digest) if row is scope else row
                       for row in plan.acquisition_partition.slices)
        updated = replace(plan, acquisition_partition=replace(plan.acquisition_partition, slices=scopes),
                          native_acquisitions=(*sources, NativeAcquisitionFile(digest, path)))
        approved = replace(updated, reviewed_acquisition_partition_sha256=acquisition_partition_sha256(record, updated))
        report = assess_validation(record, plan=approved, model_run=model, transient_run=run)
        assert report["comparison_status"] == "unassessed"
        assert report["acquisition_partition_review"]["native_acquisition_review"][reason] == [digest]
