"""Challenge verification-report truthfulness with altered and omitted evidence."""

import json
import subprocess
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from scripts import summarize_verification as reporting
from scripts import run_bounded_regression as runner
from scripts.summarize_verification import summarize, summarize_mutations
from scripts.run_bounded_regression import verification_inputs


class VerificationEvidenceAttacks(unittest.TestCase):
    def test_collection_failure_replaces_prior_success_and_retains_partial_output(self):
        cases = (
            subprocess.TimeoutExpired(["pytest"], 30, output=b"partial ids\n", stderr=b"stack\n"),
            subprocess.CompletedProcess(["pytest"], 2, b"partial ids\n", b"import failed\n"),
        )
        for process_result in cases:
            with self.subTest(failure=type(process_result).__name__), tempfile.TemporaryDirectory() as folder:
                root = Path(folder)
                (root / "src").mkdir()
                (root / "tests").mkdir()
                (root / "README.md").write_text("fixture", encoding="utf-8")
                (root / "pyproject.toml").write_text("", encoding="utf-8")
                (root / "tests/test_model.py").write_text("", encoding="utf-8")
                output = root / "regression"
                output.mkdir()
                (output / "summary.json").write_text(
                    json.dumps({"complete_suite_passed": True}), encoding="utf-8")
                run_options = ({"side_effect": process_result} if isinstance(process_result, Exception)
                               else {"return_value": process_result})
                exit_status = None
                with patch.object(runner, "ROOT", root), patch.object(runner, "OUT", output), \
                     patch.dict(runner.os.environ), \
                     patch.object(runner.sys, "argv", ["runner"]), \
                     patch.object(runner, "verification_inputs", side_effect=lambda: verification_inputs(root)), \
                     patch.object(runner.importlib.util, "find_spec", return_value=None), \
                     patch.object(runner.subprocess, "run", **run_options):
                    try:
                        exit_status = runner.main()
                    except (subprocess.SubprocessError, RuntimeError):
                        pass
                current = json.loads((output / "summary.json").read_text())
                self.assertFalse(current["complete_suite_passed"])
                self.assertEqual(exit_status, 1)
                self.assertIn(current["collection"]["status"], ("timeout", "failed"))
                self.assertIn("partial ids", (output / "collection.txt").read_text())
                self.assertEqual(current["results"], [])

    def test_surviving_invalid_or_omitted_mutants_cannot_produce_overall_success(self):
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            (root / "src").mkdir()
            (root / "scripts").mkdir()
            (root / "README.md").write_text("fixture", encoding="utf-8")
            (root / "pyproject.toml").write_text("", encoding="utf-8")
            for name in ("run_sass_mutations.py", "run_feed_mutations.py", "run_lunaris_impulse_mutations.py"):
                (root / "scripts" / name).write_text(
                    "MUTATIONS = (('fixture_fault', None),)\n", encoding="utf-8")
            native_root = root / "provider"
            native_root.mkdir()
            (native_root / "model.py").write_text("fixture", encoding="utf-8")
            snapshot = verification_inputs(root)
            regression = root / "regression"
            regression.mkdir()
            report = {"source_sha256": snapshot, "source_unchanged_during_run": True,
                      "collected_test_ids": ["tests/test_model.py::Checks::test_physics"],
                      "explicit_sidera_test_path": None,
                      "external_dependency_roots": [str(native_root)],
                      "external_dependency_sha256": runner.dependency_inputs([native_root]),
                      "results": [{"module": "tests/test_model.py", "status": "passed"}]}
            (regression / "summary.json").write_text(json.dumps(report), encoding="utf-8")
            (regression / "test_model.xml").write_text(
                '<testsuites><testsuite><testcase classname="tests.test_model.Checks" '
                'name="test_physics"/></testsuite></testsuites>', encoding="utf-8")
            killed = {"mutation": "fixture_fault", "status": "killed",
                      "assertion_failures": 1, "errors": 0,
                      "baseline_passed":True,"post_control_passed":True,"execution_valid":True}
            for status, row, expected_exit in (
                ("survived", {**killed, "status": "survived", "assertion_failures": 0}, 1),
                ("invalid", {**killed, "status": "invalid", "errors": 1}, 1),
                ("failed_baseline", {**killed,"status":"invalid","baseline_passed":False},1),
                ("failed_post_control", {**killed,"status":"invalid","post_control_passed":False},1),
                ("invalid_execution", {**killed,"status":"invalid","execution_valid":False},1),
                ("forged_kill_after_failed_baseline", {**killed,"baseline_passed":False},1),
                ("killed", killed, 0),
            ):
                with self.subTest(status=status):
                    for name in ("mutations.json", "feed_mutations.json", "native_lunaris_mutations.json"):
                        data = {"scope": "synthetic report fixture", "source_sha256": snapshot, "mutations": [row],
                                "runtime_policy":reporting.MUTATION_RUNTIME_POLICY}
                        if name == "native_lunaris_mutations.json":
                            data.update(external_dependency_roots=[str(native_root)],
                                        external_dependency_sha256=runner.dependency_inputs([native_root]))
                        (root / name).write_text(json.dumps(data), encoding="utf-8")
                    with patch.object(reporting, "ROOT", root), patch.object(reporting, "OUT", root):
                        exit_status = reporting.main()
                    self.assertEqual(exit_status, expected_exit)
                    current = json.loads((root / "verification_summary.json").read_text())
                    self.assertEqual(current["complete_suite_passed"],status!='forged_kill_after_failed_baseline')
                    self.assertEqual(current.get("software_checks_passed", False), expected_exit == 0)
                    self.assertEqual(current["verification_status"],
                                     'incomplete_or_invalid' if status=='forged_kill_after_failed_baseline'
                                     else "complete" if expected_exit == 0 else "failed")
            # Another valid native-copy report must not silently qualify a
            # different dependency context from the one used by regression.
            other = root / "other_provider"
            other.mkdir()
            (other / "model.py").write_text("fixture")
            native_path = root / "native_lunaris_mutations.json"
            native_data = json.loads(native_path.read_text())
            native_data.update(external_dependency_roots=[str(other)],
                               external_dependency_sha256=runner.dependency_inputs([other]))
            native_path.write_text(json.dumps(native_data))
            self.assertTrue(summarize_mutations(root, native_path)['passed'])
            with patch.object(reporting, "ROOT", root), patch.object(reporting, "OUT", root):
                self.assertEqual(reporting.main(), 1)
            self.assertIn("regression dependency context", json.loads(
                (root / "verification_summary.json").read_text())["reason"])
            for rows in ([], [killed, killed], [{**killed, "mutation": "unconfigured"}],
                         [{**killed, "assertion_failures": 0}], [{**killed, "status": "unknown"}],
                         [{**killed,"baseline_passed":1}], [{**killed,"execution_valid":"true"}]):
                with self.subTest(rows=rows):
                    path = root / "mutations.json"
                    path.write_text(json.dumps({"scope": "synthetic report fixture",
                        "source_sha256": snapshot, "mutations": rows,
                        "runtime_policy":reporting.MUTATION_RUNTIME_POLICY}), encoding="utf-8")
                    with self.assertRaises(ValueError):
                        summarize_mutations(root, path)

    def test_native_mutation_dependency_changes_cannot_reuse_kills(self):
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            (root / "scripts").mkdir()
            (root / "src").mkdir()
            (root / "README.md").write_text("fixture")
            (root / "pyproject.toml").write_text("")
            (root / "scripts/run_lunaris_impulse_mutations.py").write_text("MUTATIONS = (('fixture_fault', None),)\n")
            native = root / "provider"
            native.mkdir()
            module = native / "model.py"
            module.write_text("original")
            snapshot = verification_inputs(root)
            report = {"scope":"synthetic provider fixture", "source_sha256":snapshot,
                      "runtime_policy":reporting.MUTATION_RUNTIME_POLICY,
                      "external_dependency_roots":[str(native)],
                      "external_dependency_sha256":runner.dependency_inputs([native]),
                      "mutations":[{"mutation":"fixture_fault", "status":"killed", "assertion_failures":1,"errors":0,
                                    "baseline_passed":True,"post_control_passed":True,"execution_valid":True}]}
            path = root / "native_lunaris_mutations.json"
            path.write_text(json.dumps(report))
            self.assertTrue(summarize_mutations(root, path)['passed'])
            # UI build dependencies are outside the native source snapshot;
            # changing their Python helper must not invalidate this report.
            build = native / "ui/node_modules/tooling"
            build.mkdir(parents=True)
            (build / "helper.py").write_text("unrelated build helper")
            self.assertTrue(summarize_mutations(root, path)['passed'])
            module.write_text("changed")
            self.assertEqual(snapshot, verification_inputs(root))
            with self.assertRaisesRegex(ValueError, "dependency.*stale"):
                summarize_mutations(root, path)
            module.write_text("original")
            for field in ('external_dependency_roots', 'external_dependency_sha256'):
                modified = dict(report)
                del modified[field]
                path.write_text(json.dumps(modified))
                with self.subTest(field=field), self.assertRaisesRegex(ValueError, "dependency"):
                    summarize_mutations(root, path)

    def test_missing_crashed_module_evidence_replaces_previous_green_summary(self):
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            (root / "src").mkdir()
            (root / "README.md").write_text("fixture", encoding="utf-8")
            (root / "pyproject.toml").write_text("", encoding="utf-8")
            report = {"source_sha256": verification_inputs(root),
                      "source_unchanged_during_run": True,
                      "collected_test_ids": ["tests/test_model.py::Checks::test_physics"],
                      "explicit_sidera_test_path": None,
                      "results": [{"module": "tests/test_model.py", "status": "failed"}]}
            regression = root / "regression"
            regression.mkdir()
            (regression / "summary.json").write_text(json.dumps(report), encoding="utf-8")
            output = root / "verification_summary.json"
            output.write_text(json.dumps({"complete_suite_passed": True}), encoding="utf-8")
            exit_status = None
            with patch.object(reporting, "ROOT", root), patch.object(reporting, "OUT", root):
                try:
                    exit_status = reporting.main()
                except (OSError, ValueError):
                    pass
            current = json.loads(output.read_text(encoding="utf-8"))
            self.assertFalse(current["complete_suite_passed"])
            self.assertEqual(exit_status, 1)
            self.assertEqual(current["verification_status"], "incomplete_or_invalid")
            self.assertIn("FileNotFoundError", current["reason"])
            self.assertIn("test_model.xml", current["reason"])

    def test_changed_inputs_cannot_reuse_mutation_kills(self):
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            (root / "src").mkdir()
            (root / "README.md").write_text("fixture", encoding="utf-8")
            (root / "pyproject.toml").write_text("", encoding="utf-8")
            source = root / "src/model.py"
            source.write_text("original", encoding="utf-8")
            (root / "scripts").mkdir()
            (root / "scripts/run_sass_mutations.py").write_text(
                "MUTATIONS = (('fixture_fault', None),)\n", encoding="utf-8")
            path = root / "mutations.json"
            path.write_text(json.dumps({"scope": "synthetic fixture", "mutations": [
                {"mutation": "fixture_fault", "status": "killed", "assertion_failures": 1,
                 "errors": 0,"baseline_passed":True,"post_control_passed":True,"execution_valid":True}],
                 "runtime_policy":reporting.MUTATION_RUNTIME_POLICY,
                 "source_sha256": verification_inputs(root)}), encoding="utf-8")
            self.assertEqual(summarize_mutations(root, path)["killed"], 1)
            valid=path.read_text();missing=json.loads(valid);missing.pop('runtime_policy')
            path.write_text(json.dumps(missing))
            with self.assertRaisesRegex(ValueError,'runtime policy'):summarize_mutations(root,path)
            path.write_text(valid)
            source.write_text("changed", encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "stale"):
                summarize_mutations(root, path)

    def test_changed_source_or_failed_subtest_cannot_reuse_a_pass(self):
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            (root / "src").mkdir()
            (root / "README.md").write_text("fixture", encoding="utf-8")
            (root / "pyproject.toml").write_text("", encoding="utf-8")
            source = root / "src/model.py"
            source.write_text("original", encoding="utf-8")
            report = {"source_sha256": verification_inputs(root),
                      "source_unchanged_during_run": True,
                      "collected_test_ids": ["tests/test_model.py::Checks::test_physics"],
                      "explicit_sidera_test_path": None,
                      "results": [{"module": "tests/test_model.py", "status": "passed"}]}
            (root / "summary.json").write_text(json.dumps(report), encoding="utf-8")
            xml = root / "test_model.xml"
            case = '<testcase classname="tests.test_model.Checks" name="test_physics"/>'
            xml.write_text('<testsuites><testsuite>' + case + '</testsuite></testsuites>', encoding="utf-8")
            self.assertTrue(summarize(root, root)["complete_suite_passed"])
            source.write_text("changed", encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "stale"):
                summarize(root, root)
            source.write_text("original", encoding="utf-8")
            added = root / "src/new.py"
            added.write_text("new code", encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "stale"):
                summarize(root, root)
            added.unlink()
            failed = case.replace('/>', '><failure message="physical mismatch"/></testcase>')
            xml.write_text('<testsuites><testsuite>' + case + failed + '</testsuite></testsuites>',
                           encoding="utf-8")
            result = summarize(root, root)
            self.assertFalse(result["complete_suite_passed"])
            self.assertEqual(result["results"]["failed"], 1)
            xml.write_text('<testsuites><testsuite>' + case + '</testsuite></testsuites>', encoding="utf-8")
            report["results"][0]["status"] = "failed"
            (root / "summary.json").write_text(json.dumps(report), encoding="utf-8")
            self.assertFalse(summarize(root, root)["complete_suite_passed"])

    def test_missing_collected_test_cannot_disappear_from_report(self):
        with tempfile.TemporaryDirectory() as folder:
            root = Path(folder)
            (root / "src").mkdir()
            (root / "README.md").write_text("fixture", encoding="utf-8")
            (root / "pyproject.toml").write_text("", encoding="utf-8")
            source = root / "src/model.py"
            source.write_text("original", encoding="utf-8")
            report = {"source_sha256": verification_inputs(root),
                      "source_unchanged_during_run": True,
                      "collected_test_ids": ["tests/test_model.py::Checks::test_physics"],
                      "explicit_sidera_test_path": None,
                      "results": [{"module": "tests/test_model.py", "status": "passed"}]}
            (root / "summary.json").write_text(json.dumps(report), encoding="utf-8")
            (root / "test_model.xml").write_text('<testsuites><testsuite/></testsuites>', encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "identities"):
                summarize(root, root)


if __name__ == "__main__":
    unittest.main()
