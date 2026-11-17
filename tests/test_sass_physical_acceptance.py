"""Adversarial checks for approved physical limits on replayed measurements."""

from dataclasses import replace
from fractions import Fraction
from math import nextafter

import pytest

from rocket_propulsion.core.errors import DomainError
from rocket_propulsion.studies.sass import assumption
from rocket_propulsion.studies.sass.transient_validation import (
    TransientMeasurement,
    transient_measurement,
)
from rocket_propulsion.studies.sass.validation import (
    PhysicalAcceptanceCriteria,
    assess_validation,
    validation_sample_sha256,
)
from tests.test_sass_adversarial import q
from tests.test_sass_shutdown_replay import fixture, with_shutdown_prediction
from tests.test_sass_transient_validation import (
    RECIPE,
    acquisition_value,
    archived,
    channel,
    planned,
    pressure_row,
    with_fixed_prediction,
)


def limits(*, low=None, high=None, unit="Pa", timing=False):
    return PhysicalAcceptanceCriteria(
        "synthetic physical requirements approval",
        None if low is None else acquisition_value(low, unit),
        None if high is None else acquisition_value(high, unit), timing,
    )


def approve(record, target, row, criteria):
    return replace(target, physical_criteria=criteria, reviewed_sample_sha256=validation_sample_sha256(
        record, row, transient_target=target.transient_target, physical_criteria=criteria
    ))


def evaluate(tmp_path, run, recipe, row, criteria):
    record, plan = planned(tmp_path, run, (row, recipe))
    target = approve(record, plan.targets[0], row, criteria)
    return assess_validation(record, plan=replace(plan, targets=(target,)), transient_run=run)


def test_inclusive_limits_and_one_ulp_straddles_in_both_directions(tmp_path):
    cases = (
        (limits(high=nextafter(10, 0)), "fail"),
        (limits(high=10), "indeterminate"),
        (limits(high=nextafter(14, 0)), "indeterminate"),
        (limits(high=14), "pass"),
        (limits(low=nextafter(14, 20)), "fail"),
        (limits(low=14), "indeterminate"),
        (limits(low=nextafter(10, 20)), "indeterminate"),
        (limits(low=10, high=14), "pass"),
        (limits(low=11, high=13), "indeterminate"),
    )
    for criteria, expected in cases:
        report = evaluate(tmp_path, archived(), RECIPE, pressure_row(), criteria)
        physical = report["physical_acceptance_review"]
        assert report["comparison_status"] == "unassessed"
        assert physical["samples"][0]["measurement_band_status"] == expected
        assert physical["status"] == ("unassessed" if expected == "indeterminate" else expected)


def test_physical_band_uses_original_outward_interval_not_symmetrized_row(tmp_path):
    next_one = nextafter(1.0, 2.0)
    width = next_one - 1
    run = archived(channel((0, 1, width), (4, next_one, width), slew=0))
    measurement = transient_measurement(run.channels[0].trace, RECIPE)
    assert isinstance(measurement, TransientMeasurement)
    assert measurement.enclosure == (1, next_one)
    assert Fraction(measurement.value) - Fraction(measurement.absolute_error) < 1
    row = pressure_row(measured=q(measurement.value, "Pa"),
                       measurement_bound=q(measurement.absolute_error, "Pa"), predicted=q(1, "Pa"))
    report = evaluate(tmp_path, run, RECIPE, row, limits(low=1, high=next_one))
    assert report["physical_acceptance_review"]["status"] == "pass"


def test_matching_prediction_and_missing_model_archive_cannot_hide_physical_failure(tmp_path):
    row = pressure_row()
    run = archived()
    record, plan, model_run, _ = with_fixed_prediction(tmp_path, run, (row, RECIPE))
    criteria = limits(high=9)
    target = replace(plan.targets[0], transient_target=RECIPE)
    plan = replace(plan, targets=(approve(record, target, row, criteria),))
    for model, expected in ((None, "unassessed"), (model_run, "pass")):
        report = assess_validation(record, plan=plan, transient_run=run, model_run=model)
        assert report["comparison_status"] == expected
        assert report["physical_acceptance_review"]["status"] == "fail"


def test_model_outside_domain_or_failed_prediction_does_not_change_physical_pass(tmp_path):
    for row, expected in ((pressure_row(in_validation_domain=False), "unassessed"),
                          (pressure_row(predicted=q(40, "Pa")), "fail")):
        run = archived()
        record, plan, model_run, _ = with_fixed_prediction(tmp_path, run, (row, RECIPE))
        target = approve(record, plan.targets[0], row, limits(high=14))
        result = assess_validation(record, plan=replace(plan, targets=(target,)), transient_run=run, model_run=model_run)
        assert result["comparison_status"] == expected
        assert result["physical_acceptance_review"]["status"] == "pass"
        assert result["hardware_readiness"] == "not_assessed"


def test_edited_physical_limit_or_approval_cannot_reuse_prior_review(tmp_path):
    run = archived()
    row = pressure_row()
    record, plan = planned(tmp_path, run, (row, RECIPE))
    criteria = limits(high=9)
    target = approve(record, plan.targets[0], row, criteria)
    for changed in (limits(high=14), replace(criteria, approval_record="substitute"), None):
        attacked = replace(plan, targets=(replace(target, physical_criteria=changed),))
        report = assess_validation(record, plan=attacked, transient_run=run)
        assert report["comparison_status"] == "unassessed"
        assert report["samples"][0].get("reason") == "reviewed_sample_sha256_mismatch"
        assert report["physical_acceptance_review"]["status"] in ("unassessed", "not_declared")


def test_timing_failure_dominates_passing_band_and_matching_prediction(tmp_path):
    run, recipe, row = fixture(late=True)
    criteria = limits(high=18, unit="N*s", timing=True)
    record, plan, model, _ = with_shutdown_prediction(tmp_path, run, ((row, recipe),))
    target = approve(record, plan.targets[0], row, criteria)
    result = assess_validation(record, plan=replace(plan, targets=(target,)), transient_run=run, model_run=model)
    physical = result["physical_acceptance_review"]
    assert result["comparison_status"] == "pass"
    assert physical["status"] == "fail"
    assert physical["samples"][0]["measurement_band_status"] == "pass"
    assert physical["samples"][0]["valve_timing_status"] == "fail"


def test_absent_required_measurement_cannot_turn_a_partial_pass_into_complete_pass(tmp_path):
    run = archived()
    measured = pressure_row(identifier="available")
    missing = pressure_row(identifier="required-but-missing")
    record, plan = planned(tmp_path, run, (measured, RECIPE), (missing, RECIPE))
    for limit, expected in ((14, "unassessed"), (9, "fail")):
        criteria = limits(high=limit)
        approved = tuple(approve(record, target, row, criteria)
                         for target, row in zip(plan.targets, (measured, missing)))
        attacked = replace(record, samples=(measured,))
        result = assess_validation(attacked, plan=replace(plan, targets=approved), transient_run=run)
        physical = result["physical_acceptance_review"]
        assert physical["status"] == expected
        assert physical["samples"][1]["status"] == "unassessed"


def test_missing_or_reapproved_understated_transcription_has_no_physical_result(tmp_path):
    run = archived()
    criteria = limits(high=9)
    row = pressure_row()
    record, plan = planned(tmp_path, run, (row, RECIPE))
    approved = approve(record, plan.targets[0], row, criteria)
    no_archive = assess_validation(record, plan=replace(plan, targets=(approved,)))
    assert no_archive["physical_acceptance_review"]["samples"][0].get("reason") == "missing_transient_archive"
    bad = pressure_row(measurement_bound=q(.1, "Pa"))
    report = evaluate(tmp_path, run, RECIPE, bad, criteria)
    assert report["physical_acceptance_review"]["status"] == "unassessed"
    assert (report["physical_acceptance_review"]["samples"][0].get("reason")
            == "replayed_transient_measurement_mismatch")


def test_criteria_require_approved_units_and_applicable_replay(tmp_path):
    _, plan = planned(tmp_path, archived(), (pressure_row(), RECIPE))
    target = plan.targets[0]
    for kwargs in (
        {},
        {"require_valve_timing": 1},
        {"maximum": assumption(20, "Pa", "unapproved limit")},
        {"minimum": acquisition_value(10, "N"), "maximum": acquisition_value(20, "Pa")},
        {"minimum": acquisition_value(20, "Pa"), "maximum": acquisition_value(10, "Pa")},
    ):
        with pytest.raises(DomainError):
            PhysicalAcceptanceCriteria("synthetic", **kwargs)
    for kwargs in (
        {"physical_criteria": limits(high=20, unit="N")},
        {"physical_criteria": limits(timing=True)},
        {"physical_criteria": limits(high=20), "transient_target": None},
    ):
        with pytest.raises(DomainError):
            replace(target, **kwargs)



def test_timing_only_criterion_retains_indeterminate_observed_response(tmp_path):
    run, recipe, row = fixture()
    criteria = limits(timing=True)
    assert evaluate(tmp_path, run, recipe, row, criteria)["physical_acceptance_review"]["status"] == "pass"
    ambiguous = replace(recipe, limits=replace(recipe.limits, closing_ticks=3))
    report = evaluate(tmp_path, run, ambiguous, row, criteria)
    physical = report["physical_acceptance_review"]
    assert report["comparison_status"] == "unassessed"
    assert physical["status"] == "unassessed"
    assert physical["samples"][0]["valve_timing_status"] == "unassessed"
