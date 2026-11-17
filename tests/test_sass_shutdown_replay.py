"""Adversarial replay of closure-relative impulse from archived channels."""

from dataclasses import asdict, replace
from hashlib import sha256
from json import dumps, loads

import pytest

from rocket_propulsion.core.errors import DomainError, InputError
from rocket_propulsion.studies.sass import ShutdownWindow, ValveCommand, assumption
from rocket_propulsion.studies.sass.transient_validation import (
    ArchivedTransientRun,
    ArchivedValveAcquisition,
    ShutdownValidationTarget,
    shutdown_measurement,
    transient_archive_bytes,
    transient_archive_from_bytes,
)
from rocket_propulsion.studies.sass.validation import (
    PhysicalAcceptanceCriteria,
    assess_validation,
    validation_sample_sha256,
)
from tests.test_sass_adversarial import q
from tests.test_sass_shutdown_observation import histories
from tests.test_sass_transient_validation import (
    ARTICLE,
    RECIPE,
    RUN,
    acquisition_value,
    approve_synthetic_partition,
    channel,
    fixed_prediction_row,
    planned,
    pressure_row,
)
from tests.test_sass_validation import approve_prediction_bytes, archive_fixture, sample
from tests.test_sass_valve_response import limits, observed


def fixture(*, late=False):
    commanded, measured, selection = histories(duration=20)
    period = acquisition_value(1, "s")
    commanded = replace(commanded, tick_period=period)
    measured = replace(measured, tick_period=period)
    valve = ArchivedValveAcquisition(RUN, commanded, observed(measured, error=1))
    force = channel((3, 2, 0), (7, 2, 0), (11, 2, 0), (15, 2, 0), unit="N",
                    name="force", slew=0)
    run = ArchivedTransientRun(ARTICLE, (force,), (valve,))
    recipe = ShutdownValidationTarget(
        "force", selection, limits(5, 1 if late else 5),
        ShutdownWindow(15, 1, "force sync review", "tail attribution review"),
    )
    row = sample(measured=q(14, "N*s"), predicted=q(14, "N*s"), measurement_bound=q(4, "N*s"),
                 prediction_bound=q(0, "N*s"), allowable_error=q(5, "N*s"))
    return run, recipe, row


def shutdown_prediction_scope(recipe):
    return {
        "metric": "post_closure_positive_impulse",
        "channel_id": recipe.channel_id,
        "closing_commands": [{"head_id": c.head_id, "tick": c.tick} for c in recipe.close_commands],
        "window": asdict(recipe.window),
        "final_closure_ticks": [6, 10],
        "tick_period_s": 1,
    }


def with_shutdown_prediction(tmp_path, run, rows_and_recipes):
    rows, recipes = zip(*rows_and_recipes, strict=True)
    record, plan, model, export = archive_fixture(tmp_path, *rows)
    raw = tmp_path / "acquisition.json"
    raw.write_bytes(transient_archive_bytes(run))
    record = replace(record, raw_data_path=raw, raw_data_sha256=sha256(raw.read_bytes()).hexdigest())
    targets = tuple(replace(target, transient_target=recipe,
                            reviewed_sample_sha256=validation_sample_sha256(record, row, transient_target=recipe))
                    for target, row, recipe in zip(plan.targets, rows, recipes, strict=True))
    mixed = any(not isinstance(recipe, ShutdownValidationTarget) for recipe in recipes)
    export["schema"] = "rocket-propulsion-prediction-run-v5" if mixed else "rocket-propulsion-prediction-run-v4"
    for model_row, recipe in zip(export["predictions"], recipes, strict=True):
        if isinstance(recipe, ShutdownValidationTarget):
            model_row.update(pressure_referent=None, shutdown_scope=shutdown_prediction_scope(recipe))
            model_row["fixed_window" if mixed else "pressure_window"] = None
        else:
            clock = next(item.trace.tick_period.value for item in run.channels
                         if (item.condition_record, item.trace.channel_id)
                         == (model_row["condition_record"], recipe.channel_id))
            model_row.update(fixed_prediction_row(recipe, clock))
    plan = approve_prediction_bytes(replace(plan, targets=targets), model, dumps(export).encode())
    return record, approve_synthetic_partition(record, plan), model, export


def test_shutdown_archive_roundtrip_preserves_actual_closed_window(tmp_path):
    run, recipe, row = fixture()
    record, plan = planned(tmp_path, run, (row, recipe))
    archive = record.raw_data_path.read_bytes()
    assert loads(archive)["schema"] == "rocket-propulsion-transient-run-v2"
    reopened = transient_archive_from_bytes(archive)
    assert reopened == run
    assert transient_archive_bytes(reopened) == archive
    report = assess_validation(record, plan=plan, transient_run=reopened)
    assert report["comparison_status"] == "unassessed"
    assert report["samples"][0]["reason"] == "missing_shutdown_prediction_scope"
    assert report["transient_replay_review"]["status"] == "verified"
    item = report["transient_replay_review"]["samples"][0]
    assert item["enclosure"] == (10, 18)
    assert item["final_closure_ticks"] == (6, 10)
    assert item["valve_timing_status"] == "pass"


def test_no_valve_record_or_other_condition_cannot_supply_closure(tmp_path):
    run, recipe, row = fixture()
    valve = run.valve_acquisitions[0]
    for valves in ((), (replace(valve, condition_record="different firing"),)):
        changed = replace(run, valve_acquisitions=valves)
        record, plan = planned(tmp_path, changed, (row, recipe))
        report = assess_validation(record, plan=plan, transient_run=changed)
        assert report["comparison_status"] == "unassessed"
        assert report["samples"][0].get("reason") == "missing_archived_valve_acquisition"


def test_review_digest_binds_recipe_even_for_unchanged_numerical_result(tmp_path):
    run, recipe, row = fixture()
    record, plan = planned(tmp_path, run, (row, recipe))
    target = plan.targets[0]
    # Dropping A does not change final closure or the numeric result.
    for new in (
        replace(recipe, close_commands=recipe.close_commands[1:]),
        replace(recipe, window=replace(recipe.window, clock_alignment_record="new clock review")),
        replace(recipe, limits=replace(recipe.limits, approval_record="new limit review")),
    ):
        original = shutdown_measurement(run.channels[0].trace, recipe, run.valve_acquisitions[0])
        assert shutdown_measurement(run.channels[0].trace, new, run.valve_acquisitions[0]) == original
        changed = replace(plan, targets=(replace(target, transient_target=new),))
        result = assess_validation(record, plan=changed, transient_run=run)
        assert result["samples"][0].get("reason") == "reviewed_sample_sha256_mismatch"
        assert result["comparison_status"] == "unassessed"


def test_archived_alignment_and_observed_position_are_hash_bound(tmp_path):
    run, recipe, row = fixture()
    record, plan = planned(tmp_path, run, (row, recipe))
    valve = run.valve_acquisitions[0]
    event = ValveCommand(9, "B", False)
    history = valve.observed.history
    for edited in (
        replace(valve.observed, clock_alignment_record="edited clock calibration"),
        replace(valve.observed, history=replace(history, commands=(history.commands[0], event))),
        replace(valve.observed, max_relative_timing_error_ticks=0),
    ):
        mutated = replace(run, valve_acquisitions=(replace(valve, observed=edited),))
        try:
            assess_validation(record, plan=plan, transient_run=mutated)
        except DomainError as error:
            assert "reconstruction SHA-256 mismatch" in str(error)
        else:
            raise AssertionError("Changed valve evidence replayed against the original archive")


def test_reapproved_input_does_not_replace_missing_closure_or_justify_assumed_clock(tmp_path):
    run, recipe, row = fixture()
    valve = run.valve_acquisitions[0]
    history = valve.observed.history
    edits = (
        (replace(valve, observed=replace(valve.observed, history=replace(
            history, commands=history.commands[:1]))), "shutdown_closure_unresolved"),
        (replace(valve, observed=replace(valve.observed, verified_head_ids=("A",))),
         "shutdown_closure_unresolved"),
        (replace(valve, observed=replace(valve.observed, history=replace(
            history, tick_period=assumption(1, "s", "guessed clock")))), "transient_assumed_time_or_slew"),
        (replace(valve, observed=replace(valve.observed, history=replace(
            history, tick_period=acquisition_value(.5, "s")))), "shutdown_clock_scale_mismatch"),
        (replace(valve, observed=replace(valve.observed, history=replace(
            history, commands=(*history.commands, ValveCommand(17, "B", True))))),
         "shutdown_possible_valve_reopening"),
    )
    for edited, reason in edits:
        changed = replace(run, valve_acquisitions=(edited,))
        record, plan = planned(tmp_path, changed, (row, recipe))
        report = assess_validation(record, plan=plan, transient_run=changed)
        assert report["sample_review"]["status"] == "verified"
        assert report["comparison_status"] == "unassessed"
        assert report["samples"][0].get("reason") == reason


def test_reviewed_recipe_cannot_select_an_unissued_close(tmp_path):
    run, recipe, row = fixture()
    modified = replace(recipe, close_commands=(ValveCommand(6, "B", False),))
    record, plan = planned(tmp_path, run, (row, modified))
    report = assess_validation(record, plan=plan, transient_run=run)
    assert report["comparison_status"] == "unassessed"
    assert report["samples"][0].get("reason") == "shutdown_close_not_archived"


def test_late_closure_can_agree_with_prediction_but_keeps_timing_failure(tmp_path):
    run, recipe, row = fixture(late=True)
    record, plan, model, _ = with_shutdown_prediction(tmp_path, run, ((row, recipe),))
    report = assess_validation(record, plan=plan, transient_run=run, model_run=model)
    assert report["comparison_status"] == "pass"
    assert report["transient_replay_review"]["samples"][0]["valve_timing_status"] == "fail"
    assert report["qualification_evidence"] is False


def test_both_archive_versions_reject_downgrade_and_duplicate_valve_runs():
    run, _, _ = fixture()
    data = transient_archive_bytes(run)
    stripped = replace(run, valve_acquisitions=())
    assert loads(transient_archive_bytes(stripped))["schema"] == "rocket-propulsion-transient-run-v1"
    assert transient_archive_from_bytes(transient_archive_bytes(stripped)) == stripped
    with pytest.raises(DomainError, match="unique valve condition"):
        replace(run, valve_acquisitions=run.valve_acquisitions * 2)
    for key, value in (("schema", "rocket-propulsion-transient-run-v1"), ("valve_acquisitions", [])):
        attack = loads(data)
        attack[key] = value
        invalid = dumps(attack, sort_keys=True, separators=(",", ":")).encode("ascii")
        try:
            transient_archive_from_bytes(invalid)
        except (DomainError, InputError):
            continue
        raise AssertionError("An archive version change hid valve records or omitted required evidence")


def test_recipe_digest_changes_with_force_alignment_even_if_measurement_reapproved(tmp_path):
    run, recipe, row = fixture()
    record, plan = planned(tmp_path, run, (row, recipe))
    recipe2 = replace(recipe, window=replace(recipe.window, force_to_command_error_ticks=0))
    assert validation_sample_sha256(record, row, transient_target=recipe2) != plan.targets[0].reviewed_sample_sha256
    new_record, new_plan = planned(tmp_path, run, (row, recipe2))
    report = assess_validation(new_record, plan=new_plan, transient_run=run)
    assert report["comparison_status"] == "unassessed"
    assert report["samples"][0].get("reason") == "replayed_transient_measurement_mismatch"




def test_fixed_pressure_and_closed_impulse_compare_with_archived_model(tmp_path):
    run, shutdown, force = fixture()
    force = replace(force, identifier="shutdown-1")
    pressure = pressure_row(identifier="manifold-1")
    run = replace(run, channels=(*run.channels, channel((0, 10, 0), (4, 10, 0))))
    evidence, plan, model, _ = with_shutdown_prediction(tmp_path, run, ((pressure, RECIPE), (force, shutdown)))
    reopened = transient_archive_from_bytes(evidence.raw_data_path.read_bytes())
    report = assess_validation(evidence, plan=plan, model_run=model,
                               transient_run=reopened)
    assert report["comparison_status"] == "pass"
    assert report["model_run_review"]["status"] == report["transient_replay_review"]["status"] == "verified"
    pressure_result, shutdown_result = report["transient_replay_review"]["samples"]
    assert pressure_result["enclosure"] == (10, 14)
    assert pressure_result.get("final_closure_ticks") is None
    assert shutdown_result["final_closure_ticks"] == (6, 10)


def test_reapproved_shutdown_prediction_cannot_change_channel_commands_or_window(tmp_path):
    run, recipe, row = fixture()
    evidence, plan, model, export = with_shutdown_prediction(tmp_path, run, ((row, recipe),))
    criteria = PhysicalAcceptanceCriteria("approved force tail limit", maximum=acquisition_value(18, "N*s"))
    target = replace(plan.targets[0], physical_criteria=criteria,
                     reviewed_sample_sha256=validation_sample_sha256(
                         evidence, row, transient_target=recipe, physical_criteria=criteria
                     ))
    plan = replace(plan, targets=(target,))
    valid = assess_validation(evidence, plan=plan, model_run=model, transient_run=run)
    assert valid["comparison_status"] == "pass"
    original = export["predictions"][0]["shutdown_scope"]
    for changes, reason in (
        ({"channel_id": "different force channel"}, "shutdown_prediction_scope_mismatch"),
        ({"closing_commands": original["closing_commands"][1:]}, "shutdown_prediction_scope_mismatch"),
        ({"closing_commands": [{**original["closing_commands"][0], "tick": 0},
                               original["closing_commands"][1]]}, "shutdown_prediction_scope_mismatch"),
        ({"window": {**original["window"], "tail_end_tick": 14}}, "shutdown_prediction_scope_mismatch"),
        ({"window": {**original["window"], "force_to_command_error_ticks": 2}},
         "shutdown_prediction_scope_mismatch"),
        ({"window": {**original["window"], "clock_alignment_record": "different reviewed clock origin"}},
         "shutdown_prediction_scope_mismatch"),
        ({"window": {**original["window"], "tail_window_review_record": "different attribution"}},
         "shutdown_prediction_scope_mismatch"),
        ({"final_closure_ticks": [6, 9]}, "shutdown_prediction_closure_mismatch"),
        ({"final_closure_ticks": [7, 10]}, "shutdown_prediction_closure_mismatch"),
        ({"tick_period_s": .5}, "shutdown_prediction_clock_mismatch"),
    ):
        export["predictions"][0]["shutdown_scope"] = {**original, **changes}
        approved = approve_prediction_bytes(plan, model, dumps(export).encode())
        report = assess_validation(evidence, plan=approved, model_run=model, transient_run=run)
        assert report["samples"][0].get("reason") == reason
        assert report["model_run_review"]["shutdown_prediction_scope_exclusions"] == {row.identifier: reason}
        assert report["comparison_status"] == "unassessed"
        assert report["transient_replay_review"]["status"] == "verified"
        assert report["physical_acceptance_review"]["status"] == "pass"

    export["predictions"][0]["shutdown_scope"] = {
        **original, "closing_commands": list(reversed(original["closing_commands"])),
    }
    approved = approve_prediction_bytes(plan, model, dumps(export).encode())
    assert assess_validation(evidence, plan=approved, model_run=model, transient_run=run)["comparison_status"] == "pass"


def test_legacy_or_explicit_null_prediction_cannot_establish_shutdown_scope(tmp_path):
    run, recipe, row = fixture()
    evidence, plan, model, export = with_shutdown_prediction(tmp_path, run, ((row, recipe),))
    for schema, include_key in (("rocket-propulsion-prediction-run-v4", True),
                                ("rocket-propulsion-prediction-run-v3", False)):
        export["schema"] = schema
        model_row = export["predictions"][0]
        if include_key:
            model_row["shutdown_scope"] = None
        else:
            del model_row["shutdown_scope"]
        approved = approve_prediction_bytes(plan, model, dumps(export).encode())
        report = assess_validation(evidence, plan=approved, model_run=model, transient_run=run)
        assert report["comparison_status"] == "unassessed"
        assert report["samples"][0].get("reason") == "missing_shutdown_prediction_scope"


def test_prediction_cannot_use_signed_impulse_or_malformed_closure_scope(tmp_path):
    run, recipe, row = fixture()
    evidence, plan, model, export = with_shutdown_prediction(tmp_path, run, ((row, recipe),))
    original = export["predictions"][0]["shutdown_scope"]
    for changes in (
        {"metric": "signed_integral"},
        {"metric": "positive_impulse"},
        {"final_closure_ticks": [True, 10]},
        {"final_closure_ticks": [6, 10.0]},
        {"final_closure_ticks": [10, 6]},
        {"final_closure_ticks": [6, recipe.window.tail_end_tick]},
        {"closing_commands": original["closing_commands"][:1] * 2},
        {"tick_period_s": 0},
        {"window": {**original["window"], "clock_alignment_record": " "}},
    ):
        export["predictions"][0]["shutdown_scope"] = {**original, **changes}
        approved = approve_prediction_bytes(plan, model, dumps(export).encode())
        try:
            assess_validation(evidence, plan=approved, model_run=model, transient_run=run)
        except (DomainError, InputError):
            continue
        raise AssertionError("Malformed closure scope or incompatible impulse metric entered the archive")


def test_shutdown_model_bytes_without_observed_closure_cannot_verify_run(tmp_path):
    run, recipe, row = fixture()
    evidence, plan, model, _ = with_shutdown_prediction(tmp_path, run, ((row, recipe),))
    report = assess_validation(evidence, plan=plan, model_run=model)
    assert report["comparison_status"] == "unassessed"
    assert report["samples"][0].get("reason") == "missing_transient_archive"
    assert report["model_run_review"]["shutdown_prediction_scope_exclusions"] == {
        row.identifier: "shutdown_prediction_replay_unavailable",
    }


def test_shutdown_scoped_prediction_cannot_be_reused_for_an_unscoped_impulse(tmp_path):
    _, recipe, row = fixture()
    evidence, plan, model, export = archive_fixture(tmp_path, row)
    export["schema"] = "rocket-propulsion-prediction-run-v4"
    export["predictions"][0].update(pressure_referent=None, pressure_window=None,
                                    shutdown_scope=shutdown_prediction_scope(recipe))
    approved = approve_prediction_bytes(plan, model, dumps(export).encode())
    report = assess_validation(evidence, plan=approved, model_run=model)
    assert report["comparison_status"] == "unassessed"
    assert report["samples"][0].get("reason") == "unexpected_shutdown_prediction_scope"
