"""Pressure data and model referents must match the independently reviewed target."""

from dataclasses import asdict, replace
from hashlib import sha256
from json import dumps, loads

import pytest

from rocket_propulsion.core.errors import DomainError, InputError
from rocket_propulsion.studies.sass.transient_validation import (
    PressureReferent,
    transient_archive_bytes,
    transient_archive_from_bytes,
)
from rocket_propulsion.studies.sass.validation import assess_validation, validation_sample_sha256
from tests.test_sass_physical_acceptance import limits
from tests.test_sass_shutdown_replay import fixture
from tests.test_sass_transient_validation import (
    RECIPE,
    approve_synthetic_partition,
    archived,
    channel,
    fixed_prediction_row,
    pressure_row,
)
from tests.test_sass_validation import approve_prediction_bytes, archive_fixture

ABSOLUTE = PressureReferent("feed manifold", "absolute")
GAUGE = PressureReferent("feed manifold", "gauge", "test chamber")
DIFFERENTIAL = PressureReferent("feed manifold", "differential", "injector inlet")


def configure(tmp_path, *, measured=ABSOLUTE, target=ABSOLUTE, predicted=ABSOLUTE, legacy=False, v5=False):
    run = archived()
    run = replace(run, channels=(replace(run.channels[0], pressure_referent=measured),))
    raw = tmp_path / "acquisition.json"
    raw.write_bytes(transient_archive_bytes(run))
    row = pressure_row()
    record, plan, model, export = archive_fixture(tmp_path, row)
    record = replace(record, raw_data_path=raw, raw_data_sha256=sha256(raw.read_bytes()).hexdigest())
    recipe = replace(RECIPE, pressure_referent=target)
    criteria = limits(high=14)
    target_row = replace(plan.targets[0], transient_target=recipe, physical_criteria=criteria,
                         reviewed_sample_sha256=validation_sample_sha256(
                             record, row, transient_target=recipe, physical_criteria=criteria))
    plan = replace(plan, targets=(target_row,))
    if v5:
        export["schema"] = "rocket-propulsion-prediction-run-v5"
        export["predictions"][0].update(fixed_prediction_row(
            replace(recipe, pressure_referent=predicted), run.channels[0].trace.tick_period.value
        ))
    elif not legacy:
        export["schema"] = "rocket-propulsion-prediction-run-v3"
        export["predictions"][0]["pressure_window"] = None if predicted is None else {
            "start_tick": recipe.start_tick, "end_tick": recipe.end_tick, "metric": recipe.metric,
            "tick_period_s": run.channels[0].trace.tick_period.value,
            "window_review_record": recipe.window_review_record,
        }
        export["predictions"][0]["pressure_referent"] = None if predicted is None else asdict(predicted)
    plan = approve_prediction_bytes(plan, model, dumps(export).encode("utf-8"))
    return record, approve_synthetic_partition(record, plan), model, run


def test_pressure_referent_roundtrip_matches_archived_model_and_physical_measurement(tmp_path):
    record, plan, model, run = configure(tmp_path, v5=True)
    data = record.raw_data_path.read_bytes()
    assert loads(data)["schema"] == "rocket-propulsion-transient-run-v3"
    reopened = transient_archive_from_bytes(data)
    assert reopened == run
    result = assess_validation(record, plan=plan, model_run=model, transient_run=reopened)
    assert result["comparison_status"] == "pass"
    assert result["physical_acceptance_review"]["status"] == "pass"
    assert result["transient_replay_review"]["samples"][0]["pressure_referent_status"] == "matched"


def test_same_values_with_wrong_measured_location_or_datum_never_supply_a_result(tmp_path):
    for referent in (
        GAUGE, PressureReferent("injector inlet", "absolute"), None,
    ):
        record, plan, model, run = configure(tmp_path, measured=referent)
        result = assess_validation(record, plan=plan, model_run=model, transient_run=run)
        assert result["samples"][0].get("reason") == "pressure_measurement_referent_mismatch"
        assert result["comparison_status"] == "unassessed"
        assert result["physical_acceptance_review"]["status"] == "unassessed"


def test_model_location_and_datum_are_checked_even_with_identical_numbers(tmp_path):
    for wrong in (
        GAUGE, PressureReferent("injector inlet", "absolute"),
        PressureReferent("test chamber", "differential", "feed manifold"),
    ):
        record, plan, model, run = configure(tmp_path, predicted=wrong)
        result = assess_validation(record, plan=plan, model_run=model, transient_run=run)
        assert result["samples"][0].get("reason") == "pressure_prediction_referent_mismatch"
        assert result["comparison_status"] == "unassessed"
        assert result["physical_acceptance_review"]["status"] == "pass"


def test_reference_node_and_differential_polarity_cannot_change(tmp_path):
    for matching, wrong in (
        (GAUGE, replace(GAUGE, reference_location="room ambient")),
        (DIFFERENTIAL, PressureReferent("injector inlet", "differential", "feed manifold")),
    ):
        record, plan, model, run = configure(tmp_path, measured=matching, target=matching, predicted=wrong)
        report = assess_validation(record, plan=plan, model_run=model, transient_run=run)
        assert report["samples"][0].get("reason") == "pressure_prediction_referent_mismatch"
        assert report["physical_acceptance_review"]["status"] == "pass"


def test_legacy_prediction_or_missing_archive_cannot_clear_typed_pressure(tmp_path):
    for legacy in (True, False):
        record, plan, model, run = configure(tmp_path, legacy=legacy, predicted=None)
        for kwargs in ({"model_run": model}, {}):
            result = assess_validation(record, plan=plan, transient_run=run, **kwargs)
            assert result["samples"][0].get("reason") == "missing_pressure_prediction_referent"
            assert result["physical_acceptance_review"]["status"] == "pass"


def test_downgrading_pressure_plan_does_not_match_typed_acquisition(tmp_path):
    record, plan, model, run = configure(tmp_path, target=None, predicted=None)
    result = assess_validation(record, plan=plan, model_run=model, transient_run=run)
    assert result["comparison_status"] == "unassessed"
    assert result["samples"][0].get("reason") == "pressure_measurement_referent_mismatch"
    # Legacy measurements remain replayable and retain their pressure scope label.
    record, plan, model, run = configure(tmp_path, target=None, measured=None, predicted=None, legacy=True)
    result = assess_validation(record, plan=plan, model_run=model, transient_run=run)
    assert result["comparison_status"] == "unassessed"
    assert result["samples"][0].get("reason") == "missing_fixed_prediction_scope"
    assert result["transient_replay_review"]["samples"][0]["pressure_referent_status"] == "legacy_not_declared"


def test_independent_review_digest_binds_location_and_acquisition_bytes(tmp_path):
    record, plan, model, run = configure(tmp_path)
    item = plan.targets[0]
    changed = replace(item, transient_target=replace(item.transient_target, pressure_referent=GAUGE))
    result = assess_validation(record, plan=replace(plan, targets=(changed,)), model_run=model,
                               transient_run=run)
    assert result["samples"][0].get("reason") == "reviewed_sample_sha256_mismatch"
    forged = replace(run, channels=(replace(run.channels[0], pressure_referent=GAUGE),))
    try:
        assess_validation(record, plan=plan, model_run=model, transient_run=forged)
    except DomainError as error:
        assert "reconstruction SHA-256 mismatch" in str(error)
    else:
        raise AssertionError("Altered pressure reference replayed against original acquisition bytes")


def test_absolute_and_relative_pressure_require_explicit_reference_and_channel_unit(tmp_path):
    for args in (("tap", "gauge"), ("tap", "differential", "tap"), ("tap", "absolute", "room"),
                 ("tap", "unspecified"), (" ", "absolute")):
        with pytest.raises(DomainError):
            PressureReferent(*args)
    with pytest.raises(DomainError):
        replace(channel((0, 1, 0), (4, 1, 0), unit="N"), pressure_referent=ABSOLUTE)
    record, plan, model, run = configure(tmp_path)
    export = loads(model.predictions_path.read_bytes())
    del export["predictions"][0]["pressure_referent"]["reference_location"]
    invalid = approve_prediction_bytes(plan, model, dumps(export).encode("utf-8"))
    with pytest.raises((DomainError, InputError), match="Pressure referent"):
        assess_validation(record, plan=invalid, model_run=model, transient_run=run)



def test_pressure_and_valve_evidence_share_the_v3_archive(tmp_path):
    run, _, _ = fixture()
    _, _, _, pressure = configure(tmp_path)
    joint = replace(run, channels=(*run.channels, *pressure.channels))
    data = transient_archive_bytes(joint)
    assert loads(data)["schema"] == "rocket-propulsion-transient-run-v3"
    assert "valve_acquisitions" in loads(data)
    try:
        reopened = transient_archive_from_bytes(data)
    except (DomainError, InputError):
        reopened = None
    assert reopened == joint
    assert transient_archive_bytes(reopened) == data
    payload = loads(data)
    payload["schema"] = "rocket-propulsion-transient-run-v2"
    edited = dumps(payload, sort_keys=True, separators=(",", ":")).encode("ascii")
    with pytest.raises(DomainError, match="canonical"):
        transient_archive_from_bytes(edited)
