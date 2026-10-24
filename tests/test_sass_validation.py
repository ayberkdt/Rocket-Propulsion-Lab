"""Attack independence, uncertainty, provenance and result-write failure modes.

All inputs here are synthetic. Passing demonstrates software verification,
not experimental validation of any motor, feed system or scientific model.
"""

import hashlib
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
    ValidationEvidence,
    ValidationSample,
    assess_validation,
    compare_validation_sample,
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
