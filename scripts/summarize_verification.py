"""Derive current verification evidence from collected IDs, JUnit and source digests.

Missing/changed inputs, omitted tests and failed subtests cannot become passes.
This report makes no claim about test adequacy, experiment authenticity or TRL.
References: NASA-STD-7009B configuration, verification evidence and caveat reporting.
"""

from __future__ import annotations

import ast
import hashlib
import json
import xml.etree.ElementTree as ET
from pathlib import Path

from scripts.mutation_runtime import MUTATION_RUNTIME_POLICY
from scripts.run_bounded_regression import dependency_inputs, verification_inputs

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "output/assurance"


def summarize(root: Path, regression: Path) -> dict:
    """Require current sources and terminal JUnit evidence for every collected identity."""
    manifest_bytes = (regression / "summary.json").read_bytes()
    report = json.loads(manifest_bytes)
    if report.get("collection", {}).get("status", "passed") != "passed":
        raise ValueError("Test collection failed or did not finish; current verification is incomplete.")
    fingerprints = report["source_sha256"]
    if not fingerprints or not report["source_unchanged_during_run"]:
        raise ValueError("Verification source snapshot is missing or changed during the run.")
    if verification_inputs(root) != fingerprints:
        raise ValueError("Verification evidence is stale: input file set/content changed.")
    for name, digest in fingerprints.items():
        path = (root / name).resolve()
        if not path.is_relative_to(root.resolve()) or not path.is_file():
            raise ValueError("Verification input is outside the project or missing.")
        if hashlib.sha256(path.read_bytes()).hexdigest() != digest:
            raise ValueError(f"Verification evidence is stale: {name}.")
    external = report.get("external_dependency_sha256", {})
    if dependency_inputs([Path(p) for p in report.get("external_dependency_roots", [])]) != external:
        raise ValueError("External test dependency file set/content changed.")
    for name, digest in external.items():
        path = Path(name)
        if not path.is_file() or hashlib.sha256(path.read_bytes()).hexdigest() != digest:
            raise ValueError(f"External test dependency changed: {name}.")
    collected = report["collected_test_ids"]
    expected = set()
    for node_id in collected:
        module, *parts = node_id.split("::")
        expected.add(".".join((module.removesuffix(".py").replace("/", "."), *parts)))
    if not expected or len(expected) != len(collected):
        raise ValueError("Collected method identities are missing or ambiguous.")
    states = {}
    junit_rows = 0
    evidence_files = {}
    for row in report["results"]:
        if row["status"] not in {"passed", "failed"}:
            raise ValueError("A timed out or unfinished module cannot close verification.")
        xml = regression / (Path(row["module"]).stem + ".xml")
        xml_bytes = xml.read_bytes()
        evidence_files[str(xml.relative_to(root))] = hashlib.sha256(xml_bytes).hexdigest()
        for case in ET.fromstring(xml_bytes).iter("testcase"):
            identity = case.attrib["classname"] + "." + case.attrib["name"]
            state = ("failed" if case.find("failure") is not None or case.find("error") is not None
                     else "skipped" if case.find("skipped") is not None else "passed")
            states.setdefault(identity, set()).add(state)
            junit_rows += 1
    if set(states) != expected:
        raise ValueError("JUnit identities do not match collected tests; omitted/extra tests exist.")
    totals = {"passed": 0, "failed": 0, "skipped": 0}
    for state_set in states.values():
        status = "failed" if "failed" in state_set else "skipped" if "skipped" in state_set else "passed"
        totals[status] += 1
    evidence_files[str((regression / "summary.json").relative_to(root))] = hashlib.sha256(
        manifest_bytes).hexdigest()
    return {"scope": "software verification only; adequacy and physical validation not inferred",
            "collected_test_methods": len(expected), "results": totals,
            "junit_rows_including_subtests": junit_rows,
            "source_inputs_verified": len(fingerprints),
            "external_dependency_sources_verified": len(report.get("external_dependency_sha256", {})),
            "failed_modules": [row["module"] for row in report["results"] if row["status"] != "passed"],
            "complete_suite_passed": totals["failed"] == totals["skipped"] == 0 and all(
                row["status"] == "passed" for row in report["results"]),
            "explicit_sidera_test_path": report["explicit_sidera_test_path"],
            "resource_isolated_modules": report.get("resource_isolated_modules", []),
            "runtime_performance_under_concurrent_load": report.get(
                "runtime_performance_under_concurrent_load", "not_assessed"),
            "native_sidera_scope": "constant-burn public adapter; advanced provider consumption not demonstrated",
            "evidence_sha256": evidence_files,
            "hardware_readiness": "not_assessed", "qualification_evidence": False}


def _declared_mutations(root: Path, path: Path) -> set[str]:
    """Read fault identities from the source-bound literal registry without executing it."""
    script = {"mutations.json": "run_sass_mutations.py",
              "feed_mutations.json": "run_feed_mutations.py",
              "native_lunaris_mutations.json": "run_lunaris_impulse_mutations.py"}.get(path.name)
    if script is None:
        raise ValueError("Unknown curated mutation report.")
    tree = ast.parse((root / "scripts" / script).read_text(encoding="utf-8"))
    tables = [node.value for node in tree.body if isinstance(node, ast.Assign)
              and any(isinstance(target, ast.Name) and target.id == "MUTATIONS"
                      for target in node.targets)]
    if len(tables) != 1 or not isinstance(tables[0], ast.Tuple):
        raise ValueError("Curated mutation registry must be one explicit tuple.")
    identifiers = []
    for row in tables[0].elts:
        if not isinstance(row, ast.Tuple) or not row.elts or not isinstance(row.elts[0], ast.Constant):
            raise ValueError("Curated fault identifiers must be literal strings.")
        name = row.elts[0].value
        if not isinstance(name, str) or not name.strip():
            raise ValueError("Curated fault identifiers must be nonempty strings.")
        identifiers.append(name)
    if not identifiers or len(identifiers) != len(set(identifiers)):
        raise ValueError("Curated mutation registry must be nonempty and unique.")
    return set(identifiers)


def summarize_mutations(root: Path, path: Path) -> dict:
    """Require current sources and one coherent outcome per configured curated fault."""
    raw = path.read_bytes()
    data = json.loads(raw)
    if data.get('runtime_policy')!=MUTATION_RUNTIME_POLICY:
        raise ValueError('Mutation runtime policy is missing or unsupported.')
    if not data.get("source_sha256") or data["source_sha256"] != verification_inputs(root):
        raise ValueError("Mutation evidence is stale: input file set/content changed.")
    if path.name == "native_lunaris_mutations.json":
        external = data.get("external_dependency_sha256")
        roots = data.get("external_dependency_roots")
        if not external or not roots or dependency_inputs([Path(p) for p in roots]) != external:
            raise ValueError("Native mutation dependency evidence is missing or stale.")
    outcomes = data["mutations"]
    if not isinstance(outcomes, list) or not outcomes:
        raise ValueError("Curated mutation evidence must be a nonempty list.")
    names = []
    for row in outcomes:
        if not isinstance(row, dict) or not isinstance(row.get("mutation"), str):
            raise ValueError("Curated mutation outcome lacks a fault identifier.")  # noqa: TRY004 - Malformed JSON evidence retains its ValueError schema contract.
        failures, errors = row.get("assertion_failures"), row.get("errors")
        if type(failures) is not int or type(errors) is not int or min(failures, errors) < 0:
            raise ValueError("Curated outcome counts must be nonnegative integers.")
        controls=[row.get(key) for key in ('baseline_passed','post_control_passed','execution_valid')]
        if any(type(value) is not bool for value in controls):
            raise ValueError('Mutation control outcomes must be explicit booleans.')
        expected = "invalid" if not all(controls) or errors else "killed" if failures else "survived"
        if row.get("status") != expected:
            raise ValueError("Curated status contradicts its assertion/error evidence.")
        names.append(row["mutation"])
    declared = _declared_mutations(root, path)
    if len(names) != len(set(names)) or set(names) != declared:
        raise ValueError("Curated outcomes omit, duplicate or add configured fault identities.")
    return {"scope": data["scope"], "total": len(outcomes),
            "killed": sum(row["status"] == "killed" for row in outcomes),
            "invalid": sum(row["status"] == "invalid" for row in outcomes),
            "survived": sum(row["status"] == "survived" for row in outcomes),
            "declared_mutation_ids": sorted(declared),
            "passed": all(row["status"] == "killed" for row in outcomes),
            "external_dependency_sources_verified": len(data.get("external_dependency_sha256", {})),
            "sha256": hashlib.sha256(raw).hexdigest()}


def main():
    """Replace old success before evaluation; missing/crashed evidence stays failing."""
    OUT.mkdir(parents=True, exist_ok=True)
    output = OUT / "verification_summary.json"
    pending = {"verification_status": "evaluation_in_progress", "complete_suite_passed": False,
               "software_checks_passed": False,
               "hardware_readiness": "not_assessed", "qualification_evidence": False}
    output.write_text(json.dumps(pending, indent=2) + "\n", encoding="utf-8")
    try:
        result = summarize(ROOT, OUT / "regression")
        result["curated_mutation_runs"] = {
            name: summarize_mutations(ROOT, OUT / name)
            for name in ("mutations.json", "feed_mutations.json", "native_lunaris_mutations.json")}
        native_data = json.loads((OUT / "native_lunaris_mutations.json").read_text())
        regression_data = json.loads((OUT / "regression/summary.json").read_text())
        if not native_data["external_dependency_sha256"].items() <= regression_data.get(
                "external_dependency_sha256", {}).items():
            raise ValueError("Native mutation sources do not match the regression dependency context.")
        result["native_mutation_dependency_alignment_verified"] = True
        result["curated_mutation_runs_passed"] = all(
            row["passed"] for row in result["curated_mutation_runs"].values())
        result["software_checks_passed"] = (
            result["complete_suite_passed"] and result["curated_mutation_runs_passed"])
        result["verification_status"] = "complete" if result["software_checks_passed"] else "failed"
    except Exception as error:  # noqa: BLE001 - evidence evaluation must not leave stale green output
        result = {**pending, "verification_status": "incomplete_or_invalid",
                  "reason": f"{type(error).__name__}: {error}"}
    output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({key: result.get(key) for key in ("verification_status",
                     "collected_test_methods", "results", "source_inputs_verified",
                     "complete_suite_passed", "software_checks_passed", "reason")}, indent=2))
    return int(not result["software_checks_passed"])


if __name__ == "__main__":
    raise SystemExit(main())
