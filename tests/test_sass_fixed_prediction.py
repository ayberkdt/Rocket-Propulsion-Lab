"""Reject numerically identical predictions that use the wrong acquisition scope."""

from dataclasses import replace
from json import dumps
from math import nextafter, ulp

from rocket_propulsion.core.errors import DomainError, InputError
from rocket_propulsion.studies.sass.validation import assess_validation
from tests.test_sass_adversarial import q
from tests.test_sass_pressure_referent import ABSOLUTE
from tests.test_sass_shutdown_replay import fixture, shutdown_prediction_scope
from tests.test_sass_transient_validation import (
    RECIPE,
    acquisition_value,
    archived,
    channel,
    pressure_row,
    with_fixed_prediction,
)
from tests.test_sass_validation import approve_prediction_bytes, sample


def force_run(tmp_path):
    force = channel((0, 2, 0), (4, 2, 0), (8, 2, 0), name="force", unit="N", slew=0)
    other = replace(force, trace=replace(force.trace, channel_id="spare force"))
    # Same channel label in another test at a different clock scale.
    earlier = replace(force, condition_record="other firing",
                      trace=replace(force.trace, tick_period=acquisition_value(.5, "s")))
    run = archived(force, other, earlier)
    recipe = replace(RECIPE, channel_id="force", metric="positive_impulse")
    row = sample(measured=q(8, "N*s"), predicted=q(8, "N*s"), measurement_bound=q(ulp(0.0), "N*s"),
                 prediction_bound=q(0, "N*s"), allowable_error=q(.1, "N*s"))
    evidence, plan, model, export = with_fixed_prediction(tmp_path, run, (row, recipe))
    return evidence, plan, model, export, run


def test_same_impulse_on_other_channel_window_clock_or_metric_cannot_match(tmp_path):
    evidence, plan, model, export, run = force_run(tmp_path)
    report = assess_validation(evidence, plan=plan, model_run=model, transient_run=run)
    assert report["comparison_status"] == "pass"
    assert report["transient_replay_review"]["samples"][0]["enclosure"] == (8, 8)
    row = export["predictions"][0]
    original = row["fixed_window"]
    for changes, reason in (
        ({"channel_id": "spare force"}, "fixed_prediction_scope_mismatch"),
        ({"start_tick": 1}, "fixed_prediction_scope_mismatch"),
        ({"end_tick": 5}, "fixed_prediction_scope_mismatch"),
        ({"start_tick": 4, "end_tick": 8}, "fixed_prediction_scope_mismatch"),
        ({"metric": "signed_integral"}, "fixed_prediction_scope_mismatch"),
        ({"window_review_record": "independent time origin"}, "fixed_prediction_scope_mismatch"),
        ({"tick_period_s": .5}, "fixed_prediction_clock_mismatch"),
        ({"tick_period_s": nextafter(1, 2)}, "fixed_prediction_clock_mismatch"),
    ):
        row["fixed_window"] = {**original, **changes}
        approved = approve_prediction_bytes(plan, model, dumps(export).encode())
        report = assess_validation(evidence, plan=approved, model_run=model, transient_run=run)
        assert report["comparison_status"] == "unassessed"
        assert report["samples"][0].get("reason") == reason
        assert report["model_run_review"]["fixed_prediction_scope_exclusions"] == {row["sample_id"]: reason}
        assert report["transient_replay_review"]["status"] == "verified"


def test_valid_typed_pressure_uses_exact_ticks_beyond_binary64_integer_precision(tmp_path):
    start = -(2**55)
    capture = channel((start, 10, 0), (start + 4, 10, 0), slew=0)
    capture = replace(capture, pressure_referent=ABSOLUTE,
                      trace=replace(capture.trace, tick_period=acquisition_value(.25, "s")))
    run = archived(capture)
    recipe = replace(RECIPE, start_tick=start, end_tick=start + 4, metric="signed_integral",
                     pressure_referent=ABSOLUTE)
    row = pressure_row(measured=q(10, "Pa*s"), predicted=q(10, "Pa*s"),
                       measurement_bound=q(ulp(0.0), "Pa*s"), prediction_bound=q(0, "Pa*s"),
                       allowable_error=q(.1, "Pa*s"))
    evidence, plan, model, export = with_fixed_prediction(tmp_path, run, (row, recipe))
    assert assess_validation(evidence, plan=plan, model_run=model, transient_run=run)["comparison_status"] == "pass"
    scope = export["predictions"][0]["fixed_window"]
    assert float(start) == float(start + 1)
    scope["start_tick"] += 1
    approved = approve_prediction_bytes(plan, model, dumps(export).encode())
    report = assess_validation(evidence, plan=approved, model_run=model, transient_run=run)
    assert report["samples"][0].get("reason") == "fixed_prediction_scope_mismatch"


def test_legacy_numeric_prediction_and_null_scope_leave_valid_measurement_unassessed(tmp_path):
    evidence, plan, model, export, run = force_run(tmp_path)
    without_model = assess_validation(evidence, plan=plan, transient_run=run)
    assert without_model["samples"][0].get("reason") == "missing_fixed_prediction_scope"
    row = export["predictions"][0]
    row["fixed_window"] = None
    for version in (5, 4, 1):
        export["schema"] = f"rocket-propulsion-prediction-run-v{version}"
        if version == 4:
            del row["fixed_window"]
            row["pressure_window"] = None
        if version == 1:
            for key in ("pressure_window", "pressure_referent", "shutdown_scope"):
                del row[key]
        approved = approve_prediction_bytes(plan, model, dumps(export).encode())
        report = assess_validation(evidence, plan=approved, model_run=model, transient_run=run)
        assert report["comparison_status"] == "unassessed"
        assert report["samples"][0].get("reason") == "missing_fixed_prediction_scope"
        assert report["transient_replay_review"]["status"] == "verified"


def test_invalid_fixed_metric_clock_or_dual_window_cannot_enter_archive(tmp_path):
    evidence, plan, model, export, run = force_run(tmp_path)
    row = export["predictions"][0]
    original = row["fixed_window"]
    for changes in (
        {"metric": "peak"}, {"metric": "post_closure_positive_impulse"}, {"channel_id": " "},
        {"start_tick": False}, {"start_tick": 0.0}, {"end_tick": 0},
        {"tick_period_s": 0}, {"tick_period_s": True}, {"window_review_record": ""},
    ):
        row["fixed_window"] = {**original, **changes}
        approved = approve_prediction_bytes(plan, model, dumps(export).encode())
        try:
            assess_validation(evidence, plan=approved, model_run=model, transient_run=run)
        except (DomainError, InputError):
            continue
        raise AssertionError("Invalid fixed-window scope was admitted")
    _, shutdown, _ = fixture()
    row["fixed_window"] = original
    row["shutdown_scope"] = shutdown_prediction_scope(shutdown)
    approved = approve_prediction_bytes(plan, model, dumps(export).encode())
    reason = None
    try:
        assess_validation(evidence, plan=approved, model_run=model, transient_run=run)
    except (DomainError, InputError) as error:
        reason = str(error)
    assert reason is not None and "both fixed and closure-relative windows" in reason
