"""Attack independence, uncertainty, provenance and result-write failure modes.

All inputs here are synthetic. Passing demonstrates software verification,
not experimental validation of any motor, feed system or scientific model.
"""

import hashlib
import json
import math
import random
import sys
import tempfile
import unittest
from dataclasses import replace
from fractions import Fraction
from pathlib import Path
from unittest.mock import patch

from rocket_propulsion.core.errors import DomainError
from rocket_propulsion.studies.sass.artifacts import distinct_artifact_paths, write_artifact_pair
from rocket_propulsion.studies.sass.validation import (
    ArchivedModelRun,
    ValidationEvidence,
    ValidationPlan,
    ValidationSample,
    ValidationTarget,
    assess_validation,
    compare_validation_sample,
    validation_sample_sha256,
)
from tests.test_sass_adversarial import q

ROOT = Path(__file__).resolve().parents[1]


def sample(**changes):
    """NASA-7009B: synthetic, explicitly bounded comparison for test isolation."""
    return replace(ValidationSample("holdout-1", "synthetic shared conditions", q(10, "N"),
                                   q(10.1, "N"), q(.1, "N"), q(.1, "N"),
                                   q(.4, "N"), True), **changes)


def evidence(path, samples):
    """NASA-7009B: synthetic evidence records; no authentic experiment asserted."""
    return ValidationEvidence("synthetic-article", "synthetic-model", path,
                              hashlib.sha256(path.read_bytes()).hexdigest(), "a" * 64,
                              ("calibration-1",), samples, "synthetic independence review",
                              "synthetic transcription review", "synthetic absolute bounds",
                              "synthetic acceptance approval")


class ValidationAttacks(unittest.TestCase):
    def test_reported_residual_interval_contains_exact_corner_errors(self):
        record = sample(measured=q(0, "N"), predicted=q(3.0941676495717276, "N"),
                        measurement_bound=q(5.050647017735236, "N"),
                        prediction_bound=q(7.246650282258111, "N"),
                        allowable_error=q(15.391464949565075, "N"))
        actual_high = sum(map(Fraction, (record.predicted.value, record.measurement_bound.value,
                                        record.prediction_bound.value)))
        result = compare_validation_sample(record)
        self.assertGreaterEqual(Fraction(result["maximum_error"]), actual_high)
        self.assertEqual(result["status"], "indeterminate")

    def test_sampled_bounds_and_decisions_are_sound_against_exact_interval_corners(self):
        rng = random.Random(20261018)
        statuses = set()
        for index in range(256):
            scale = 10. ** rng.randrange(-250, 251)
            measured, predicted = (rng.uniform(-10, 10) * scale for _ in range(2))
            mb, pb = ((rng.random() + .01) * scale for _ in range(2))
            m, p, um, up = map(Fraction, (measured, predicted, mb, pb))
            corners = tuple(y - x for x in (m - um, m + um) for y in (p - up, p + up))
            high = max(map(abs, corners))
            low = Fraction(0) if min(corners) <= 0 <= max(corners) else min(map(abs, corners))
            allowance = float(high * 2 if index % 3 == 0 else low / 2 if index % 3 == 1 else high)
            result = compare_validation_sample(sample(
                measured=q(measured, "N"), predicted=q(predicted, "N"),
                measurement_bound=q(mb, "N"), prediction_bound=q(pb, "N"),
                allowable_error=q(allowance, "N")))
            statuses.add(result["status"])
            with self.subTest(index=index):
                self.assertGreaterEqual(Fraction(result["maximum_error"]), high)
                self.assertLessEqual(Fraction(result["minimum_error"]), low)
                if result["status"] == "pass":
                    self.assertLess(high, Fraction(allowance))
                if result["status"] == "fail":
                    self.assertGreater(low, Fraction(allowance))
        self.assertEqual(statuses, {"pass", "fail", "indeterminate"})

    def test_declared_numerical_review_guard_is_not_silently_removed(self):
        record = sample(measured=q(0, "N"), predicted=q(1, "N"),
                        measurement_bound=q(1e-10, "N"), prediction_bound=q(0, "N"),
                        allowable_error=q(1.000000000100001, "N"))
        self.assertEqual(compare_validation_sample(record)["status"], "indeterminate")

    def test_exact_upper_limit_beyond_float_max_cannot_be_published_as_finite(self):
        record = sample(measured=q(0, "N"), predicted=q(sys.float_info.max, "N"),
                        measurement_bound=q(math.nextafter(0., 1.), "N"), prediction_bound=q(0, "N"),
                        allowable_error=q(sys.float_info.max, "N"))
        with self.assertRaises(DomainError):
            compare_validation_sample(record)

    def test_subnormal_guard_is_an_upper_enclosure_instead_of_exact_zero(self):
        smallest = math.nextafter(0., 1.)
        result = compare_validation_sample(sample(measured=q(0, "N"), predicted=q(0, "N"),
                                           measurement_bound=q(smallest, "N"), prediction_bound=q(0, "N"),
                                           allowable_error=q(smallest, "N")))
        self.assertEqual(result["maximum_error"], smallest)
        exact_guard = 16 * Fraction(sys.float_info.epsilon) * Fraction(smallest)
        self.assertGreaterEqual(Fraction(result["numerical_guard"]), exact_guard)
        self.assertEqual(result["status"], "indeterminate")

    def test_uncertainty_cannot_make_a_bad_prediction_pass(self):
        self.assertEqual(compare_validation_sample(sample())["status"], "pass")
        self.assertEqual(compare_validation_sample(sample(predicted=q(11, "N")))["status"], "fail")
        self.assertEqual(compare_validation_sample(sample(measurement_bound=q(.4, "N")))["status"],
                         "indeterminate")
        self.assertEqual(compare_validation_sample(sample(allowable_error=q(.3, "N")))["status"],
                         "indeterminate")

    def test_missing_bounds_and_outside_domain_cannot_pass(self):
        for change in ({"measurement_bound": None}, {"prediction_bound": None},
                       {"in_validation_domain": False}):
            self.assertEqual(compare_validation_sample(sample(**change))["status"], "unassessed")

    def test_bounds_have_units_and_measured_bound_cannot_be_zero(self):
        for change in ({"measurement_bound": q(0, "N")}, {"measurement_bound": q(-1, "N")},
                       {"prediction_bound": q(1, "Pa")}, {"predicted": q(10, "s")},
                       {"in_validation_domain": 1}):
            with self.subTest(change=change), self.assertRaises(DomainError):
                sample(**change)

    def test_no_sample_reuse_or_duplicate_holdout_identity(self):
        with tempfile.TemporaryDirectory() as folder:
            path = Path(folder) / "raw.csv"
            path.write_text("synthetic raw evidence", encoding="utf-8")
            for samples in ((sample(), sample()), (sample(identifier="calibration-1"),)):
                with self.assertRaises(DomainError):
                    evidence(path, samples)

    def test_changed_file_cannot_reuse_a_pass(self):
        with tempfile.TemporaryDirectory() as folder:
            path = Path(folder) / "raw.csv"
            path.write_text("synthetic raw evidence", encoding="utf-8")
            record = evidence(path, (sample(),))
            result = assess_validation(record)
            self.assertEqual(result["comparison_status"], "pass")
            self.assertFalse(result["qualification_evidence"])
            self.assertEqual(result["hardware_readiness"], "not_assessed")
            path.write_text("altered raw evidence", encoding="utf-8")
            with self.assertRaises(DomainError):
                assess_validation(record)

    def test_one_failed_sample_dominates_unknown_and_passed_samples(self):
        with tempfile.TemporaryDirectory() as folder:
            path = Path(folder) / "raw.csv"
            path.write_text("synthetic raw evidence", encoding="utf-8")
            samples = (sample(), sample(identifier="holdout-2", predicted=q(99, "N")),
                       sample(identifier="holdout-3", prediction_bound=None))
            self.assertEqual(assess_validation(evidence(path, samples))["comparison_status"], "fail")


def reviewed_plan(record):
    return ValidationPlan(
        "synthetic-article", "test-only coverage review",
        tuple(ValidationTarget(item.identifier, item.condition_record, item.allowable_error,
                               validation_sample_sha256(record, item))
              for item in record.samples),
        ("calibration-run",),
    )


class ValidationCoverageAttacks(unittest.TestCase):
    def test_missing_local_sensor_cannot_inherit_aggregate_peak_agreement(self):
        # Williams, Ares I cold-flow thesis, Table 3.2, WH5b. Only the
        # measured and modeled pressures are from the publication. Bounds
        # and a 60 psia allowance are test inputs, with no hardware claim.
        pressures = (
            ("101", 640, 638), ("103", 746, 766), ("304", 515, 529),
            ("306", 521, 525), ("310", 420, 526), ("312", 437, 522),
            ("401", 662, 636), ("402", 637, 627), ("403", 776, 801), ("407", 614, 644),
        )
        samples = tuple(ValidationSample(
            f"WH5b:HfP-{sensor}", "WH5b", q(measured, "psia"), q(predicted, "psia"),
            q(1, "psia"), q(1, "psia"), q(60, "psia"), True,
        ) for sensor, measured, predicted in pressures)
        self.assertEqual(max(row[2] for row in pressures) - max(row[1] for row in pressures), 25)
        with tempfile.TemporaryDirectory() as folder:
            path = Path(folder) / "raw.csv"
            path.write_text("\n".join(f"{site},{obs},{pred}" for site, obs, pred in pressures),
                            encoding="utf-8")
            record = evidence(path, samples)
            plan = reviewed_plan(record)
            complete = assess_validation(record, plan=plan)
            self.assertEqual(complete["coverage"]["status"], "complete")
            self.assertEqual(complete["sample_review"]["status"], "verified")
            self.assertEqual(complete["comparison_status"], "fail")
            self.assertEqual({row["sample_id"] for row in complete["samples"] if row["status"] == "fail"},
                             {"WH5b:HfP-310", "WH5b:HfP-312"})
            cherry_picked = samples[:4] + samples[6:]
            legacy = assess_validation(evidence(path, cherry_picked))
            self.assertEqual(legacy["comparison_status"], "pass")
            self.assertEqual(legacy["coverage"]["status"], "not_declared")
            incomplete = assess_validation(evidence(path, cherry_picked), plan=plan)
            self.assertEqual(incomplete["comparison_status"], "unassessed")
            self.assertEqual(incomplete["coverage"]["status"], "incomplete")
            self.assertEqual(incomplete["coverage"]["missing_sample_ids"],
                             ["WH5b:HfP-310", "WH5b:HfP-312"])
            still_failing = assess_validation(evidence(path, samples[:-1]), plan=plan)
            self.assertEqual(still_failing["comparison_status"], "fail")

    def test_scope_and_allowance_cannot_be_changed_by_result_rows(self):
        original = sample()
        with tempfile.TemporaryDirectory() as folder:
            path = Path(folder) / "raw.csv"
            path.write_text("test-only measurements", encoding="utf-8")
            record = evidence(path, (original,))
            plan = reviewed_plan(record)
            accepted = assess_validation(record, plan=plan)
            self.assertEqual(accepted["comparison_status"], "pass")
            for item, reason in (
                (replace(original, condition_record="different firing"), "condition_record"),
                (replace(original, allowable_error=q(99, "N")), "allowable_error"),
                (replace(original, measured=q(10, "Pa"), predicted=q(10.1, "Pa"),
                         measurement_bound=q(.1, "Pa"), prediction_bound=q(.1, "Pa"),
                         allowable_error=q(.4, "Pa")), "unit"),
            ):
                with self.subTest(reason=reason):
                    report = assess_validation(evidence(path, (item,)), plan=plan)
                    self.assertEqual(report["comparison_status"], "unassessed")
                    self.assertEqual(report["samples"][0]["status"], "unassessed")
                    self.assertEqual(report["coverage"]["target_mismatches"],
                                     {original.identifier: [reason]})
            extra = sample(identifier="unplanned", predicted=q(99, "N"))
            with_extra = assess_validation(evidence(path, (original, extra)), plan=plan)
            self.assertEqual(with_extra["comparison_status"], "unassessed")
            self.assertEqual(with_extra["coverage"]["status"], "incomplete")
            self.assertEqual(with_extra["coverage"]["unexpected_sample_ids"], ["unplanned"])
            other_article = replace(plan, test_article_id="different-hardware")
            with self.assertRaisesRegex(DomainError, "test article"):
                assess_validation(evidence(path, (original,)), plan=other_article)

    def test_same_run_cannot_cross_calibration_boundary_under_a_new_sensor_id(self):
        original = sample(condition_record="firing-A")
        with tempfile.TemporaryDirectory() as folder:
            path = Path(folder) / "raw.csv"
            path.write_text("test-only measurements", encoding="utf-8")
            record = evidence(path, (original,))
            plan = reviewed_plan(record)
            with self.assertRaisesRegex(DomainError, "used for calibration"):
                replace(plan, calibration_condition_records=("firing-A",))
            with self.assertRaisesRegex(DomainError, "uniquely identified"):
                replace(plan, targets=(plan.targets[0], plan.targets[0]))
            incomplete_plan = replace(plan, calibration_condition_records=())
            report = assess_validation(record, plan=incomplete_plan)
            self.assertEqual(report["comparison_status"], "unassessed")
            self.assertTrue(report["coverage"]["calibration_conditions_unreported"])


class ValidationTranscriptAttacks(unittest.TestCase):
    def test_archived_file_hash_alone_does_not_validate_an_edited_row(self):
        with tempfile.TemporaryDirectory() as folder:
            path = Path(folder) / "raw.csv"
            path.write_text("holdout-1,10,11\n", encoding="utf-8")
            original = sample(predicted=q(11, "N"))
            record = evidence(path, (original,))
            plan = reviewed_plan(record)
            self.assertEqual(assess_validation(record, plan=plan)["comparison_status"], "fail")
            changed = replace(record, samples=(replace(original, predicted=q(10, "N")),))
            self.assertEqual(assess_validation(changed)["comparison_status"], "pass")
            for row in (changed.samples[0],
                        replace(original, measured=q(11, "N")),
                        replace(original, measured=q(math.nextafter(10.0, 11.0), "N")),
                        replace(original, measurement_bound=q(0.11, "N")),
                        replace(original, prediction_bound=None),
                        replace(original, in_validation_domain=False)):
                with self.subTest(row=row):
                    report = assess_validation(replace(record, samples=(row,)), plan=plan)
                    self.assertEqual(report["raw_data_sha256"], record.raw_data_sha256)
                    self.assertEqual(report["coverage"]["status"], "complete")
                    self.assertEqual(report["comparison_status"], "unassessed")
                    self.assertEqual(report["samples"][0]["reason"], "reviewed_sample_sha256_mismatch")
                    check = report["sample_review"]["samples"][0]
                    self.assertEqual(check["expected_sha256"], plan.targets[0].reviewed_sample_sha256)
                    self.assertNotEqual(check["submitted_sha256"], check["expected_sha256"])

    def test_independent_row_review_cannot_follow_a_replaced_archive_or_model(self):
        with tempfile.TemporaryDirectory() as folder:
            path = Path(folder) / "raw.csv"
            path.write_text("holdout-1,10,10.1\n", encoding="utf-8")
            original = sample()
            record = evidence(path, (original,))
            plan = reviewed_plan(record)
            accepted = assess_validation(record, plan=plan)
            self.assertEqual(accepted["comparison_status"], "pass")
            unreviewed = replace(plan, targets=(replace(plan.targets[0], reviewed_sample_sha256=None),))
            report = assess_validation(record, plan=unreviewed)
            self.assertEqual(report["comparison_status"], "unassessed")
            self.assertEqual(report["samples"][0]["reason"], "missing_reviewed_sample_sha256")
            other_archive = path.with_name("other.csv")
            other_archive.write_text("holdout-1,10,10.1\nextra,0,0\n", encoding="utf-8")
            for edited in (
                replace(record, raw_data_path=other_archive,
                        raw_data_sha256=hashlib.sha256(other_archive.read_bytes()).hexdigest()),
                replace(record, model_sha256="b" * 64),
                replace(record, model_configuration_id="different model revision"),
                replace(record, calibration_sample_ids=("calibration-other",)),
                replace(record, transcription_review_record="different transcription review"),
            ):
                with self.subTest(edited=edited):
                    report = assess_validation(edited, plan=plan)
                    self.assertEqual(report["comparison_status"], "unassessed")
                    self.assertEqual(report["samples"][0]["reason"], "reviewed_sample_sha256_mismatch")


def archive_fixture(folder, *samples):
    raw = folder / "measurements.csv"
    raw.write_bytes(b"test-only measurements\n")
    model = folder / "model.json"
    model.write_bytes(b'{"configuration":"independent archive fixture"}\n')
    record = replace(evidence(raw, samples), model_sha256=hashlib.sha256(model.read_bytes()).hexdigest())
    predictions = {
        "schema": "rocket-propulsion-prediction-run-v1",
        "test_article_id": record.test_article_id,
        "model_configuration_id": record.model_configuration_id,
        "model_sha256": record.model_sha256,
        "predictions": [
            {"sample_id": item.identifier, "condition_record": item.condition_record,
             "unit": item.predicted.unit, "predicted": item.predicted.value,
             "prediction_bound": None if item.prediction_bound is None else item.prediction_bound.value}
            for item in samples
        ],
    }
    archives = ArchivedModelRun(model, folder / "predictions.json")
    plan = approve_prediction_bytes(reviewed_plan(record), archives, json.dumps(predictions).encode("utf-8"))
    return record, plan, archives, predictions


def approve_prediction_bytes(plan, archives, data):
    archives.predictions_path.write_bytes(data)
    return replace(plan, reviewed_predictions_sha256=hashlib.sha256(data).hexdigest())


class ArchivedPredictionAttacks(unittest.TestCase):
    def test_pass_requires_declared_archives_and_actual_matching_bytes(self):
        with tempfile.TemporaryDirectory() as folder:
            record, plan, archives, _ = archive_fixture(Path(folder), sample())
            verified = assess_validation(record, plan=plan, model_run=archives)
            self.assertEqual(verified["comparison_status"], "pass")
            self.assertEqual(verified["model_run_review"]["status"], "verified")
            self.assertEqual(verified["model_run_review"]["prediction_archive_sha256"],
                             plan.reviewed_predictions_sha256)
            self.assertEqual(verified["model_configuration_verification"],
                             "archived_model_bytes_match_declared_sha256")
            self.assertFalse(verified["qualification_evidence"])
            missing = assess_validation(record, plan=plan)
            self.assertEqual(missing["comparison_status"], "unassessed")
            self.assertEqual(missing["samples"][0]["reason"], "missing_model_run_archives")
            with self.assertRaisesRegex(DomainError, "independently reviewed"):
                assess_validation(record, plan=reviewed_plan(record), model_run=archives)
            original = archives.model_path.read_bytes()
            archives.model_path.write_bytes(original + b" ")
            with self.assertRaisesRegex(DomainError, "Archived model SHA-256 mismatch"):
                assess_validation(record, plan=plan, model_run=archives)
            archives.model_path.write_bytes(original)
            # A whitespace edit changes no JSON values but invalidates the archived bytes.
            archives.predictions_path.write_bytes(archives.predictions_path.read_bytes() + b" ")
            with self.assertRaisesRegex(DomainError, "Prediction archive SHA-256 mismatch"):
                assess_validation(record, plan=plan, model_run=archives)

    def test_approved_row_cannot_replace_prediction_or_bound_from_same_run(self):
        original = sample(predicted=q(11, "N"))
        with tempfile.TemporaryDirectory() as folder:
            record, plan, archives, _ = archive_fixture(Path(folder), original)
            self.assertEqual(assess_validation(record, plan=plan, model_run=archives)["comparison_status"],
                             "fail")
            for edited in (
                replace(original, predicted=q(10, "N")),
                replace(original, predicted=q(math.nextafter(11.0, 12.0), "N")),
                replace(original, prediction_bound=q(0.0, "N")),
                replace(original, prediction_bound=None),
                replace(original, condition_record="another firing"),
            ):
                with self.subTest(edited=edited):
                    changed = replace(record, samples=(edited,))
                    independently_approved = replace(
                        reviewed_plan(changed), reviewed_predictions_sha256=plan.reviewed_predictions_sha256)
                    report = assess_validation(changed, plan=independently_approved, model_run=archives)
                    self.assertEqual(report["sample_review"]["status"], "verified")
                    self.assertEqual(report["comparison_status"], "unassessed")
                    self.assertEqual(report["samples"][0]["reason"], "archived_prediction_mismatch")
                    self.assertEqual(report["model_run_review"]["mismatched_submitted_sample_ids"],
                                     [edited.identifier])

    def test_missing_archive_row_preserves_a_matched_failure(self):
        failing = sample(predicted=q(11, "N"))
        passing = sample(identifier="holdout-2")
        with tempfile.TemporaryDirectory() as folder:
            record, plan, archives, export = archive_fixture(Path(folder), failing, passing)
            export["predictions"].pop()
            plan = approve_prediction_bytes(plan, archives, json.dumps(export).encode("utf-8"))
            report = assess_validation(record, plan=plan, model_run=archives)
            self.assertEqual(report["comparison_status"], "fail")
            self.assertEqual(report["model_run_review"]["missing_planned_sample_ids"], [passing.identifier])
            self.assertEqual([row["status"] for row in report["samples"]], ["fail", "unassessed"])
            self.assertEqual(report["samples"][1]["reason"], "missing_archived_prediction")

    def test_hash_match_does_not_override_schema_or_configuration_identity(self):
        with tempfile.TemporaryDirectory() as folder:
            record, plan, archives, export = archive_fixture(Path(folder), sample())
            wrong_header = {**export, "model_configuration_id": "another revision"}
            duplicates = {**export, "predictions": export["predictions"] * 2}
            valid = json.dumps(export)
            for data, reason in (
                (json.dumps(wrong_header).encode(), "header does not match"),
                (json.dumps(duplicates).encode(), "duplicate sample ID"),
                (valid.replace('"predicted": 10.1', '"predicted": 99, "predicted": 10.1').encode(),
                 "JSON contract"),
                (valid.replace('"predicted": 10.1', '"predicted": 1e-999').encode(), "nonzero finite binary range"),
                (valid.replace('"predicted": 10.1', '"predicted": true').encode(), "real number"),
                (b"\xff", "UTF-8"),
            ):
                with self.subTest(reason=reason):
                    approved = approve_prediction_bytes(plan, archives, data)
                    with self.assertRaisesRegex(DomainError, reason):
                        assess_validation(record, plan=approved, model_run=archives)


class ArtifactFailureTests(unittest.TestCase):
    def test_existing_samefile_alias_is_rejected(self):
        with tempfile.TemporaryDirectory() as folder:
            first, second = Path(folder) / "input.json", Path(folder) / "alias.json"
            first.write_text("source", encoding="utf-8")
            second.write_text("source", encoding="utf-8")
            with patch.object(Path, "samefile", return_value=True), self.assertRaises(DomainError):
                distinct_artifact_paths((first, second))

    def test_hard_link_to_request_cannot_be_an_output(self):
        with tempfile.TemporaryDirectory() as folder:
            first, second = Path(folder) / "input.json", Path(folder) / "alias.json"
            first.write_text("source", encoding="utf-8")
            try:
                second.hardlink_to(first)
            except PermissionError:
                self.skipTest("Host sandbox prohibits creating hard links; samefile branch tested separately.")
            with self.assertRaises(DomainError):
                distinct_artifact_paths((first, second))
            self.assertEqual(first.read_text(encoding="utf-8"), "source")

    def test_second_replace_failure_restores_first_artifact(self):
        parent = ROOT / 'output' / 'assurance' / 'test_artifacts' / 'sass_transaction'
        parent.mkdir(parents=True, exist_ok=True)
        with tempfile.TemporaryDirectory(dir=parent) as folder:
            first, second = Path(folder) / "study.json", Path(folder) / "study.html"
            first.write_text("previous JSON", encoding="utf-8")
            second.write_text("previous HTML", encoding="utf-8")
            original_replace = Path.replace
            def fail_second(path, target):
                if target == second:
                    raise OSError("Injected second-artifact failure")
                return original_replace(path, target)
            with patch.object(Path, "replace", fail_second), self.assertRaises(OSError):
                write_artifact_pair(((first, "new JSON"), (second, "new HTML")))
            self.assertEqual(first.read_text(encoding="utf-8"), "previous JSON")
            self.assertEqual(second.read_text(encoding="utf-8"), "previous HTML")
            self.assertEqual(sorted(x.name for x in Path(folder).iterdir()),
                             ["study.html", "study.json"])


if __name__ == "__main__":
    unittest.main()
