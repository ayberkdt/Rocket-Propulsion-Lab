"""Exercise acquisition-aligned model prediction windows end to end."""

from dataclasses import asdict, replace
from hashlib import sha256
from json import dumps, loads

from rocket_propulsion.core.errors import DomainError, InputError
from rocket_propulsion.studies.sass.transient_validation import transient_archive_bytes
from rocket_propulsion.studies.sass.validation import assess_validation, validation_sample_sha256
from tests.test_sass_adversarial import q
from tests.test_sass_pressure_referent import ABSOLUTE, configure
from tests.test_sass_transient_validation import RECIPE, pressure_row
from tests.test_sass_validation import approve_prediction_bytes, archive_fixture


def test_reapproved_different_prediction_window_cannot_match_same_numeric_peak(tmp_path):
    record, plan, model, run = configure(tmp_path)
    original = loads(model.predictions_path.read_bytes())
    for edits in (
        {"start_tick": 1}, {"end_tick": 5},
        {"start_tick": 4, "end_tick": 8},
        {"window_review_record": "different time origin alignment"},
    ):
        export = loads(dumps(original))
        export["predictions"][0]["pressure_window"].update(edits)
        approved = approve_prediction_bytes(plan, model, dumps(export).encode("utf-8"))
        result = assess_validation(record, plan=approved, model_run=model, transient_run=run)
        assert result["comparison_status"] == "unassessed"
        assert result["samples"][0].get("reason") == "pressure_prediction_window_mismatch"
        assert result["physical_acceptance_review"]["status"] == "pass"


def test_same_tick_numbers_with_different_clock_scale_cannot_match(tmp_path):
    record, plan, model, run = configure(tmp_path)
    export = loads(model.predictions_path.read_bytes())
    export["predictions"][0]["pressure_window"]["tick_period_s"] = 0.5
    approved = approve_prediction_bytes(plan, model, dumps(export).encode("utf-8"))
    report = assess_validation(record, plan=approved, model_run=model, transient_run=run)
    assert report["samples"][0].get("reason") == "pressure_prediction_clock_mismatch"
    assert report["physical_acceptance_review"]["status"] == "pass"


def test_v2_archive_can_identify_pressure_but_cannot_establish_its_window(tmp_path):
    record, plan, model, run = configure(tmp_path)
    export = loads(model.predictions_path.read_bytes())
    export["schema"] = "rocket-propulsion-prediction-run-v2"
    del export["predictions"][0]["pressure_window"]
    approved = approve_prediction_bytes(plan, model, dumps(export).encode("utf-8"))
    report = assess_validation(record, plan=approved, model_run=model, transient_run=run)
    assert report["samples"][0].get("reason") == "missing_pressure_prediction_window"
    assert report["comparison_status"] == "unassessed"


def test_metric_with_wrong_unit_and_invalid_clock_cannot_enter_archive(tmp_path):
    record, plan, model, run = configure(tmp_path)
    original = loads(model.predictions_path.read_bytes())
    for edits in (
        {"metric": "signed_integral"}, {"metric": "maximum_sample"},
        {"start_tick": True}, {"start_tick": 0.0}, {"end_tick": 0},
        {"tick_period_s": 0.0}, {"tick_period_s": True},
        {"window_review_record": " "},
    ):
        export = loads(dumps(original))
        export["predictions"][0]["pressure_window"].update(edits)
        approved = approve_prediction_bytes(plan, model, dumps(export).encode("utf-8"))
        try:
            assess_validation(record, plan=approved, model_run=model, transient_run=run)
        except (DomainError, InputError):
            continue
        raise AssertionError("Invalid pressure metric or time base entered the prediction archive")


def test_partial_v3_pressure_metadata_is_rejected(tmp_path):
    record, plan, model, run = configure(tmp_path)
    original = loads(model.predictions_path.read_bytes())
    for field in ("pressure_window", "pressure_referent"):
        export = loads(dumps(original))
        export["predictions"][0][field] = None
        approved = approve_prediction_bytes(plan, model, dumps(export).encode("utf-8"))
        try:
            assess_validation(record, plan=approved, model_run=model, transient_run=run)
        except (DomainError, InputError) as error:
            assert "declared together" in str(error)
        else:
            raise AssertionError("Partial v3 pressure metadata was accepted")


def test_valid_integral_with_large_exact_ticks_and_fractional_clock(tmp_path):
    _, _, _, run = configure(tmp_path)
    channel = run.channels[0]
    origin = 10**20
    trace = replace(channel.trace, samples=tuple(replace(x, tick=x.tick + origin)
                                                 for x in channel.trace.samples),
                    tick_period=replace(channel.trace.tick_period, value=.25))
    run = replace(run, channels=(replace(channel, trace=trace),))
    raw = tmp_path / "large-clock.json"
    raw.write_bytes(transient_archive_bytes(run))
    # Four ticks at 0.25 s. The pressure is bounded by two equal-endpoint
    # cones with a slope limit of 2 Pa/s, giving an integral of 10 +/- 0.5 Pa*s.
    row = pressure_row(measured=q(10, "Pa*s"), predicted=q(10, "Pa*s"),
                       measurement_bound=q(.5, "Pa*s"), prediction_bound=q(0, "Pa*s"),
                       allowable_error=q(1, "Pa*s"))
    record, plan, model, export = archive_fixture(tmp_path, row)
    record = replace(record, raw_data_path=raw, raw_data_sha256=sha256(raw.read_bytes()).hexdigest())
    recipe = replace(RECIPE, start_tick=origin, end_tick=origin + 4, metric="signed_integral",
                     pressure_referent=ABSOLUTE)
    plan = replace(plan, targets=(replace(plan.targets[0], transient_target=recipe,
                                        reviewed_sample_sha256=validation_sample_sha256(
                                            record, row, transient_target=recipe)),))
    export["schema"] = "rocket-propulsion-prediction-run-v3"
    model_row = export["predictions"][0]
    model_row["pressure_referent"] = asdict(ABSOLUTE)
    model_row["pressure_window"] = {
        "start_tick": origin, "end_tick": origin + 4, "tick_period_s": .25,
        "metric": "signed_integral", "window_review_record": recipe.window_review_record,
    }
    approved = approve_prediction_bytes(plan, model, dumps(export).encode("utf-8"))
    result = assess_validation(record, plan=approved, model_run=model, transient_run=run)
    assert result["comparison_status"] == "unassessed"
    assert result["samples"][0].get("reason") == "missing_fixed_prediction_scope"
    model_row["pressure_window"]["start_tick"] += 1
    approved = approve_prediction_bytes(plan, model, dumps(export).encode("utf-8"))
    result = assess_validation(record, plan=approved, model_run=model, transient_run=run)
    assert result["samples"][0].get("reason") == "pressure_prediction_window_mismatch"

