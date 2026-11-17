"""Attack archive identity, recipe binding and uncertainty replay end to end."""

from dataclasses import asdict, replace
from fractions import Fraction
from hashlib import sha256
from json import dumps, loads
from math import nextafter, ulp

import pytest

from rocket_propulsion.core.errors import DomainError, InputError
from rocket_propulsion.studies.sass import (
    Citation,
    CitedScalar,
    TransientSample,
    TransientTrace,
    assumption,
)
from rocket_propulsion.studies.sass.transient_validation import (
    ArchivedTransientChannel,
    ArchivedTransientRun,
    TransientMeasurement,
    TransientValidationTarget,
    transient_archive_bytes,
    transient_archive_from_bytes,
    transient_measurement,
)
from rocket_propulsion.studies.sass.validation import (
    AcquisitionPartition,
    NativeAcquisitionFile,
    NativeAcquisitionSlice,
    ValidationPlan,
    ValidationTarget,
    acquisition_partition_sha256,
    assess_validation,
    validation_sample_sha256,
)
from tests.test_sass_adversarial import q
from tests.test_sass_validation import approve_prediction_bytes, archive_fixture, evidence, sample

RUN = "synthetic shared conditions"
ARTICLE = "synthetic-article"
RECIPE = TransientValidationTarget("feed pressure", 0, 4, "peak", "synthetic window review")


def acquisition_value(value, unit):
    # Synthetic substitute for a reviewed acquisition record, not empirical data.
    return CitedScalar(value, unit, (Citation("NASA-7009B", "synthetic test fixture"),), "derived")


def channel(*rows, condition=RUN, unit="Pa", name="feed pressure", gap=4, slew=2):
    trace = TransientTrace(name, unit, acquisition_value(1, "s"), tuple(TransientSample(*r) for r in rows),
                           "synthetic DAQ", "synthetic acquisition review", gap,
                           None if slew is None else acquisition_value(slew, f"{unit}/s"))
    return ArchivedTransientChannel(condition, trace)


def archived(*channels):
    return ArchivedTransientRun(ARTICLE, channels or (channel((0, 10, 0), (4, 10, 0)),))


def planned(tmp_path, run, *pairs):
    raw_path = tmp_path / "acquisition.json"
    raw_path.write_bytes(transient_archive_bytes(run))
    items = tuple(row for row, _ in pairs)
    records = evidence(raw_path, items)
    targets = tuple(ValidationTarget(
        row.identifier, row.condition_record, row.allowable_error,
        validation_sample_sha256(records, row, transient_target=recipe), recipe,
    ) for row, recipe in pairs)
    return records, ValidationPlan(ARTICLE, "synthetic plan review", targets, ("synthetic calibration",))


def pressure_row(**kwargs):
    row = sample(measured=q(12, "Pa"), predicted=q(12, "Pa"), measurement_bound=q(2, "Pa"),
                 prediction_bound=q(0, "Pa"), allowable_error=q(3, "Pa"))
    return replace(row, **kwargs)


def fixed_prediction_row(recipe, clock):
    return {
        "pressure_referent": None if recipe.pressure_referent is None else asdict(recipe.pressure_referent),
        "shutdown_scope": None,
        "fixed_window": {
            "channel_id": recipe.channel_id, "start_tick": recipe.start_tick, "end_tick": recipe.end_tick,
            "tick_period_s": clock, "metric": recipe.metric, "window_review_record": recipe.window_review_record,
        },
    }


def approve_synthetic_partition(record, plan):
    conditions = set(plan.calibration_condition_records) | {target.condition_record for target in plan.targets}
    calibration, = plan.calibration_condition_records
    partition = AcquisitionPartition("synthetic partition review", tuple(
        NativeAcquisitionSlice(condition, condition, sha256(condition.encode()).hexdigest(), 0, 8)
        for condition in sorted(conditions)
    ), tuple((identifier, calibration) for identifier in record.calibration_sample_ids))
    plan = replace(plan, acquisition_partition=partition,
                   native_acquisitions=synthetic_native_files(record, sorted(conditions)))
    return replace(plan, reviewed_acquisition_partition_sha256=acquisition_partition_sha256(record, plan))


def synthetic_native_files(record, labels):
    sources = []
    for label in labels:
        payload = label.encode()
        digest = sha256(payload).hexdigest()
        path = record.raw_data_path.with_name(f"synthetic-native-{digest}.daq")
        path.write_bytes(payload)
        sources.append(NativeAcquisitionFile(digest, path))
    return tuple(sources)


def with_fixed_prediction(tmp_path, run, *pairs):
    rows = tuple(row for row, _ in pairs)
    base, model_plan, model_run, export = archive_fixture(tmp_path, *rows)
    observed, observed_plan = planned(tmp_path, run, *pairs)
    record = replace(base, raw_data_path=observed.raw_data_path, raw_data_sha256=observed.raw_data_sha256)
    targets = tuple(replace(target, reviewed_sample_sha256=validation_sample_sha256(
        record, row, transient_target=recipe
    )) for target, (row, recipe) in zip(observed_plan.targets, pairs, strict=True))
    export["schema"] = "rocket-propulsion-prediction-run-v5"
    clocks = {(item.condition_record, item.trace.channel_id): item.trace.tick_period.value for item in run.channels}
    for model_row, (row, recipe) in zip(export["predictions"], pairs, strict=True):
        model_row.update(fixed_prediction_row(recipe, clocks[(row.condition_record, recipe.channel_id)]))
    plan = approve_prediction_bytes(replace(model_plan, targets=targets), model_run, dumps(export).encode())
    return record, approve_synthetic_partition(record, plan), model_run, export


def test_valid_pressure_and_fixed_impulse_replay_on_same_run(tmp_path):
    # The archived pressure tent permits peaks from 10 to 14 Pa.
    # Signed force reverses from -2 to +2 N, leaving 2 N*s forward impulse.
    impulse = replace(RECIPE, channel_id="load cell", metric="positive_impulse")
    run = archived(channel((0, 10, 0), (4, 10, 0)),
                   channel((0, -2, 0), (4, 2, 0), unit="N", name="load cell", slew=1))
    thrust_row = sample(identifier="fixed-window-forward", measured=q(2, "N*s"),
                        predicted=q(2, "N*s"), measurement_bound=q(ulp(0.0), "N*s"),
                        prediction_bound=q(0, "N*s"), allowable_error=q(0.1, "N*s"))
    submission, plan, model, _ = with_fixed_prediction(tmp_path, run, (pressure_row(), RECIPE), (thrust_row, impulse))
    reopened = transient_archive_from_bytes(submission.raw_data_path.read_bytes())
    assert transient_archive_bytes(reopened) == submission.raw_data_path.read_bytes()
    report = assess_validation(submission, plan=plan, transient_run=reopened, model_run=model)
    assert report["comparison_status"] == "pass"
    replay = report["transient_replay_review"]
    assert replay["status"] == "verified"
    assert [row["enclosure"] for row in replay["samples"]] == [(10, 14), (2, 2)]
    assert report["hardware_readiness"] == "not_assessed"
    assert report["qualification_evidence"] is False


def test_no_archive_cannot_be_treated_as_a_verified_transcription(tmp_path):
    submission, plan = planned(tmp_path, archived(), (pressure_row(), RECIPE))
    result = assess_validation(submission, plan=plan)
    assert result["comparison_status"] == "unassessed"
    assert result["samples"][0]["reason"] == "missing_transient_archive"


def test_reapproved_understatement_still_fails_raw_replay(tmp_path):
    run = archived()
    for bad in (pressure_row(measurement_bound=q(.1, "Pa"), allowable_error=q(1, "Pa")),
                pressure_row(measured=q(11, "Pa"), allowable_error=q(4, "Pa")),
                pressure_row(measurement_bound=None)):
        submission, plan = planned(tmp_path, run, (bad, RECIPE))
        result = assess_validation(submission, plan=plan, transient_run=run)
        assert result["sample_review"]["status"] == "verified"
        assert result["comparison_status"] == "unassessed"
        assert result["samples"][0].get("reason") == "replayed_transient_measurement_mismatch"


def test_synchronized_channel_name_cannot_cross_physical_runs(tmp_path):
    run = archived(channel((0, 10, 0), (4, 10, 0), condition="different test"))
    submission, plan = planned(tmp_path, run, (pressure_row(), RECIPE))
    result = assess_validation(submission, plan=plan, transient_run=run)
    assert result["comparison_status"] == "unassessed"
    assert result["samples"][0]["reason"] == "missing_archived_transient_channel"


def test_original_bytes_and_reconstructed_trace_are_both_bound(tmp_path):
    run = archived()
    submission, plan = planned(tmp_path, run, (pressure_row(), RECIPE))
    record = run.channels[0]
    for replacement in (
        replace(record.trace, samples=(TransientSample(0, 10, 0), TransientSample(4, 10, .1))),
        replace(record.trace, acquisition_review_record="different calibration"),
    ):
        forged = archived(replace(record, trace=replacement))
        try:
            assess_validation(submission, plan=plan, transient_run=forged)
        except DomainError as error:
            assert "reconstruction SHA-256 mismatch" in str(error)
        else:
            raise AssertionError("altered replay trace accepted against original raw bytes")

    submission.raw_data_path.write_bytes(submission.raw_data_path.read_bytes() + b"\n")
    try:
        assess_validation(submission, plan=plan, transient_run=run)
    except DomainError as error:
        assert "Raw experimental evidence SHA-256 mismatch" in str(error)
    else:
        raise AssertionError("altered raw bytes accepted against original review digest")


def test_same_numeric_value_cannot_reuse_a_digest_after_recipe_change(tmp_path):
    run = archived(channel((0, 12, 0), (2, 12, 0), (4, 12, 0), slew=0))
    row = pressure_row(measurement_bound=q(ulp(0.0), "Pa"))
    submission, plan = planned(tmp_path, run, (row, RECIPE))
    target = plan.targets[0]
    for changed in (replace(RECIPE, start_tick=2),
                    replace(RECIPE, window_review_record="substituted window review")):
        attacked = replace(plan, targets=(replace(target, transient_target=changed),))
        result = assess_validation(submission, plan=attacked, transient_run=run)
        assert result["comparison_status"] == "unassessed"
        assert result["samples"][0]["reason"] == "reviewed_sample_sha256_mismatch"


def test_incomplete_acquisition_cannot_supply_a_measurement(tmp_path):
    good = channel((0, 10, 0), (4, 10, 0))
    for capture, reason in (
        (channel((0, 10, 0), (3, 10, 0)), "missing_window_boundary"),
        (channel((0, 10, 0), (4, 10, 0), gap=3), "sample_gap"),
        (channel((0, 10, 0), (4, 10, 0), slew=None), "missing_slew"),
        (channel((0, 0, 0), (1, 50, 50), (4, 100, 0)), "inconsistent_samples"),
        (replace(good, trace=replace(good.trace, max_absolute_slew=assumption(2, "Pa/s", "guess"))),
         "assumed_time_or_slew"),
        (replace(good, trace=replace(good.trace, tick_period=assumption(1, "s", "guess"))),
         "assumed_time_or_slew"),
    ):
        run = archived(capture)
        submission, plan = planned(tmp_path, run, (pressure_row(), RECIPE))
        report = assess_validation(submission, plan=plan, transient_run=run)
        assert report["comparison_status"] == "unassessed"
        assert report["samples"][0]["reason"] == "transient_" + reason


def test_metric_dimensions_and_sign_do_not_cross(tmp_path):
    pressure = archived()
    bad = replace(RECIPE, metric="positive_impulse")
    submission, plan = planned(tmp_path, pressure, (pressure_row(), bad))
    result = assess_validation(submission, plan=plan, transient_run=pressure)
    assert result["comparison_status"] == "unassessed"
    assert result["samples"][0]["reason"] == "transient_metric_unit_mismatch"
    for metric, unit, enclosure in (
        ("peak", "Pa", (10, 14)),
        ("signed_integral", "Pa*s", (32, 48)),
    ):
        measurement = transient_measurement(pressure.channels[0].trace, replace(RECIPE, metric=metric))
        assert isinstance(measurement, TransientMeasurement)
        assert (measurement.unit, measurement.enclosure) == (unit, enclosure)


def test_midpoint_rounding_cannot_shrink_the_transient_interval():
    observed = channel((0, 0, 1), (1, .5, nextafter(.5, float("inf"))), slew=0)
    result = transient_measurement(observed.trace, replace(RECIPE, end_tick=1))
    assert isinstance(result, TransientMeasurement)
    lo, hi = map(Fraction, result.enclosure)
    assert (lo, hi) == (-Fraction(1, 2**53), 1)
    center, bound = Fraction(result.value), Fraction(result.absolute_error)
    assert center - bound <= lo and center + bound >= hi
    narrowed = Fraction(nextafter(result.absolute_error, 0))
    assert center - narrowed > lo or center + narrowed < hi


def test_archive_identity_keeps_signed_zero_and_deterministic_channel_order():
    first = channel((0, 0.0, 0), (4, 0, 0))
    second = channel((0, 0, 0), (4, 0, 0), condition="different run")
    assert transient_archive_bytes(archived(first, second)) == transient_archive_bytes(archived(second, first))
    negative = replace(first, trace=replace(
        first.trace, samples=(TransientSample(0, -0.0, 0), first.trace.samples[1])))
    assert transient_archive_bytes(archived(negative, second)) != transient_archive_bytes(archived(first, second))
    with pytest.raises(DomainError, match="unique run and channel"):
        archived(first, negative)



def test_replay_and_model_archive_are_both_required_when_planned(tmp_path):
    pressure = pressure_row()
    run = archived()
    record, plan, model_run, _ = with_fixed_prediction(tmp_path, run, (pressure, RECIPE))
    ok = assess_validation(record, plan=plan, model_run=model_run, transient_run=run)
    assert ok["comparison_status"] == "pass"
    assert ok["model_run_review"]["status"] == ok["transient_replay_review"]["status"] == "verified"
    assert assess_validation(record, plan=plan, model_run=model_run)["comparison_status"] == "unassessed"
    assert assess_validation(record, plan=plan, transient_run=run)["comparison_status"] == "unassessed"


def test_archive_reader_rejects_noncanonical_and_ambiguous_records():
    data = transient_archive_bytes(archived())
    assert transient_archive_bytes(transient_archive_from_bytes(data)) == data
    extra = loads(data)
    extra["channels"][0]["trace"]["ignored"] = True
    invalid_tick = loads(data)
    invalid_tick["channels"][0]["trace"]["samples"][0][0] = True
    for bad in (
        data + b" ",
        data.replace(b'"channel_id":', b'"channel_id":"override","channel_id":', 1),
        dumps(extra, sort_keys=True, separators=(",", ":")).encode("ascii"),
        dumps(invalid_tick, sort_keys=True, separators=(",", ":")).encode("ascii"),
        data.replace(b"0x1.4000000000000p+3", b"inf", 1),
    ):
        try:
            transient_archive_from_bytes(bad)
        except (DomainError, InputError):
            continue
        raise AssertionError("an ambiguous or noncanonical acquisition export was accepted")
