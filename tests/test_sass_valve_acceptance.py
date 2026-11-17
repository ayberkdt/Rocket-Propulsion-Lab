"""Keep independently verified valve timing visible when impulse evidence fails."""

from dataclasses import replace

from rocket_propulsion.studies.sass import ValveCommand, assumption
from rocket_propulsion.studies.sass.transient_validation import (
    transient_archive_bytes,
    transient_archive_from_bytes,
)
from rocket_propulsion.studies.sass.validation import assess_validation
from tests.test_sass_adversarial import q
from tests.test_sass_physical_acceptance import approve, evaluate, limits
from tests.test_sass_shutdown_replay import fixture
from tests.test_sass_transient_validation import acquisition_value, planned


def test_confirmed_late_close_survives_force_channel_absence(tmp_path):
    run, recipe, row = fixture(late=True)
    archive = replace(run, channels=())
    reopened = transient_archive_from_bytes(transient_archive_bytes(archive))
    assert reopened == archive
    criteria = limits(high=18, unit="N*s", timing=True)
    result = evaluate(tmp_path, reopened, recipe, row, criteria)
    assert result["comparison_status"] == "unassessed"
    assert result["physical_acceptance_review"]["status"] == "fail"
    physical = result["physical_acceptance_review"]["samples"][0]
    assert physical["reason"] == "missing_archived_transient_channel"
    assert physical["measurement_band_status"] == "unassessed"
    assert physical["valve_timing_status"] == "fail"
    failed = [entry for entry in physical["valve_timing_review"]["responses"]
              if entry["status"] == "late_response"]
    assert [(entry["command"]["head_id"], entry["latency_bounds_ticks"]) for entry in failed] == [
        ("B", (2, 4)),
    ]


def test_passing_timing_only_criterion_does_not_supply_an_impulse(tmp_path):
    run, recipe, row = fixture()
    archive = replace(run, channels=())
    for criteria, expected, band in (
        (limits(timing=True), "pass", "not_declared"),
        (limits(high=18, unit="N*s", timing=True), "unassessed", "unassessed"),
        (limits(high=18, unit="N*s"), "unassessed", "unassessed"),
    ):
        result = evaluate(tmp_path, archive, recipe, row, criteria)
        assert result["comparison_status"] == "unassessed"
        assert result["physical_acceptance_review"]["status"] == expected
        physical = result["physical_acceptance_review"]["samples"][0]
        assert physical["measurement_band_status"] == band
        assert physical["valve_timing_status"] == ("pass" if criteria.require_valve_timing else "not_declared")
        assert physical["enclosure"] is None


def test_unexpected_reopening_remains_failure_even_when_tail_is_unassessable(tmp_path):
    run, recipe, row = fixture()
    valve = run.valve_acquisitions[0]
    history = valve.observed.history
    reopened = replace(history, commands=(*history.commands, ValveCommand(17, "B", True)))
    valve = replace(valve, observed=replace(valve.observed, history=reopened))
    run = replace(run, valve_acquisitions=(valve,))
    report = evaluate(tmp_path, run, recipe, row, limits(timing=True))
    assert report["samples"][0].get("reason") == "shutdown_possible_valve_reopening"
    physical = report["physical_acceptance_review"]["samples"][0]
    assert physical["status"] == "fail"
    assert physical["valve_timing_review"]["observation_issues"] == (
        {"observation": {"tick": 17, "head_id": "B", "opened": True},
         "status": "unexpected_observation"},
    )


def test_unreviewed_clock_or_substituted_condition_cannot_supply_timing(tmp_path):
    run, recipe, row = fixture(late=True)
    valve = run.valve_acquisitions[0]
    observed = valve.observed
    for edited, reason in (
        (replace(valve, condition_record="another firing"), "missing_archived_valve_acquisition"),
        (replace(valve, commanded=replace(
            valve.commanded, tick_period=assumption(1, "s", "guess"))), "shutdown_assumed_valve_clock"),
        (replace(valve, observed=replace(observed, history=replace(
            observed.history, tick_period=acquisition_value(.5, "s")))), "shutdown_valve_clock_scale_mismatch"),
    ):
        archive = replace(run, channels=(), valve_acquisitions=(edited,))
        report = evaluate(tmp_path, archive, recipe, row, limits(timing=True))
        physical = report["physical_acceptance_review"]["samples"][0]
        assert physical["status"] == "unassessed"
        assert physical.get("valve_timing_reason") == reason
        assert physical.get("valve_timing_review") is None


def test_force_transcription_mismatch_cannot_erase_verified_timing(tmp_path):
    run, recipe, row = fixture(late=True)
    erroneous = replace(row, measurement_bound=q(.01, "N*s"))
    report = evaluate(tmp_path, run, recipe, erroneous, limits(high=18, unit="N*s", timing=True))
    physical = report["physical_acceptance_review"]["samples"][0]
    assert report["comparison_status"] == "unassessed"
    assert physical["status"] == "fail"
    assert physical["measurement_band_status"] == "unassessed"
    assert physical.get("reason") == "replayed_transient_measurement_mismatch"


def test_required_timing_still_requires_reviewed_recipe_and_selected_commands(tmp_path):
    run, recipe, row = fixture(late=True)
    archive = replace(run, channels=())
    record, plan = planned(tmp_path, archive, (row, recipe))
    approved = approve(record, plan.targets[0], row, limits(timing=True))
    invented = replace(recipe, close_commands=(ValveCommand(6, "B", False),))
    forged = replace(approved, transient_target=invented)
    report = assess_validation(record, plan=replace(plan, targets=(forged,)), transient_run=archive)
    physical = report["physical_acceptance_review"]["samples"][0]
    assert physical["status"] == "unassessed"
    assert physical.get("valve_timing_reason") == "reviewed_sample_sha256_mismatch"

    report = evaluate(tmp_path, archive, invented, row, limits(timing=True))
    physical = report["physical_acceptance_review"]["samples"][0]
    assert physical["status"] == "unassessed"
    assert physical.get("valve_timing_reason") == "shutdown_close_not_archived"

