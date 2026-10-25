"""Run independent test modules in child processes with recorded time budgets.

Timeouts and deferred tests are open results, never passes. Each child has its
own Python process so a hung native integration cannot prevent other checks.
References: NASA-7009B, preserve verification results and unresolved concerns.
"""

import argparse
import concurrent.futures
import hashlib
import importlib.metadata
import importlib.util
import json
import os
import subprocess
import sys
import time
import xml.etree.ElementTree as ET
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "output/assurance/regression"
DEFERRED = (
    "tests/test_lunaris_integration.py::LunarisEndToEndTests",
    ("tests/test_adversarial_pulse_mode.py::DeltaVPlanAttackTests::"
    "test_random_thrusters_and_hostile_requests_against_enumeration"),
    ("tests/test_adversarial_pulse_mode.py::ImpulsePlanAttackTests::"
    "test_random_impulse_requests_against_enumeration"),
)
MODULE_TIMEOUT_SECONDS = 120
COLLECTION_TIMEOUT_SECONDS = 30


def verification_inputs(root=ROOT):
    """Hash executable sources, tests and their repository data/document inputs."""
    suffixes = {".py", ".json", ".csv", ".txt", ".md", ".html", ".css", ".js",
                ".mjs", ".cjs", ".yml", ".yaml", ".toml", ".patch"}
    files = {root / "pyproject.toml", root / "README.md"}
    files.update(p for name in ('app.cmd', 'LICENSE') if (p := root / name).is_file())
    for folder in ("src", "tests", "scripts", "docs", "examples", ".github/workflows"):
        files.update(p for p in (root / folder).rglob("*") if p.is_file() and p.suffix in suffixes)
    return {str(p.relative_to(root)): hashlib.sha256(p.read_bytes()).hexdigest()
            for p in sorted(files)}


def dependency_inputs(roots):
    """Bind native Python sources, excluding the snapshot's build/cache trees."""
    result = {}
    ignored = {".git", "__pycache__", "node_modules"}
    for directory in roots:
        for current, directories, filenames in os.walk(directory):
            directories[:] = [name for name in directories if name.casefold() not in ignored]
            for name in filenames:
                if Path(name).suffix == ".py":
                    path = (Path(current) / name).resolve()
                    result[str(path)] = hashlib.sha256(path.read_bytes()).hexdigest()
    return result


def run_module(path):
    """Retain stdout, JUnit and explicit timeout/error status for one module."""
    stem = path.stem
    xml = OUT / (stem + ".xml")
    xml.unlink(missing_ok=True)
    cmd = [sys.executable, "-m", "pytest", path.relative_to(ROOT).as_posix(), "-q",
           "--tb=short", "-o", "faulthandler_timeout=0", "-p", "scripts.pytest_diagnostics",
           "--diagnostic-timeout=20", f"--junitxml={xml}"]
    cmd.extend("--deselect=" + name for name in DEFERRED)
    started = time.monotonic()
    try:
        process = subprocess.run(cmd, check=False, cwd=ROOT, capture_output=True, timeout=MODULE_TIMEOUT_SECONDS)
        log = process.stdout + process.stderr
        status = "passed" if process.returncode == 0 else "failed"
        return_code = process.returncode
    except subprocess.TimeoutExpired as error:
        log = (error.stdout or b"") + (error.stderr or b"")
        status = "timeout"
        return_code = None
    (OUT / (stem + ".txt")).write_bytes(log)
    result = {"module": path.relative_to(ROOT).as_posix(), "status": status,
              "return_code": return_code,
              "seconds": round(time.monotonic() - started, 2), "command": cmd}
    if xml.exists():
        suites = ET.parse(xml).getroot().findall("testsuite")
        result["junit"] = {key: sum(int(s.attrib.get(key, 0)) for s in suites)
                           for key in ("tests", "failures", "errors", "skipped")}
    print(path.name, status, flush=True)
    return result


def main():
    """Summarize every module without upgrading missing work to success."""
    global DEFERRED
    parser = argparse.ArgumentParser()
    parser.add_argument("--include-deferred", action="store_true")
    parser.add_argument("--sidera-path", type=Path)
    options = parser.parse_args()
    if options.include_deferred:
        DEFERRED = ()
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "summary.json").write_text(json.dumps({
        "complete_suite_passed": False, "collection": {"status": "not_started"},
        "results": [], "scope": "current run in progress; previous success invalidated",
    }, indent=2) + "\n", encoding="utf-8")
    paths = sorted((ROOT / "tests").glob("test_*.py"))
    os.environ["PYTHONPATH"] = str(ROOT / "src") + os.pathsep + str(ROOT)
    if options.sidera_path is not None:
        explicit_path = options.sidera_path.resolve(strict=True)
        if not explicit_path.is_dir():
            raise ValueError("The explicit Sidera test path must be a directory.")
        os.environ["PYTHONPATH"] += os.pathsep + str(explicit_path)
    os.environ["PYTHONDONTWRITEBYTECODE"] = "1"
    cache = OUT / "numba-cache"
    cache.mkdir(exist_ok=True)
    os.environ["NUMBA_CACHE_DIR"] = str(cache)
    fingerprints = verification_inputs()
    dependency_roots = [] if options.sidera_path is None else [explicit_path]
    lunaris_spec = importlib.util.find_spec("lunaris")
    if lunaris_spec is not None and lunaris_spec.origin is not None:
        dependency_roots.append(Path(lunaris_spec.origin).parent)
    dependency_hashes = dependency_inputs(dependency_roots)
    collect_command = [sys.executable, "-m", "pytest", "--collect-only", "-q"]
    started = time.monotonic()
    try:
        process = subprocess.run(collect_command, check=False, cwd=ROOT, capture_output=True,
                                 timeout=COLLECTION_TIMEOUT_SECONDS)
        collect_stdout, collect_stderr = process.stdout, process.stderr
        collection = {"status": "passed" if process.returncode == 0 else "failed",
                      "return_code": process.returncode}
    except subprocess.TimeoutExpired as error:
        collect_stdout, collect_stderr = error.stdout or b"", error.stderr or b""
        collection = {"status": "timeout", "return_code": None}
    except OSError as error:
        collect_stdout, collect_stderr = b"", str(error).encode("utf-8")
        collection = {"status": "failed", "return_code": None}
    collection.update({"seconds": round(time.monotonic() - started, 2),
                       "timeout_seconds": COLLECTION_TIMEOUT_SECONDS, "command": collect_command})
    (OUT / "collection.txt").write_bytes(collect_stdout + collect_stderr)
    collected_ids = [line for line in collect_stdout.decode("utf-8", errors="replace").splitlines()
                     if line.startswith("tests/") and "::" in line]
    if collection["status"] == "passed" and (
            not collected_ids or len(collected_ids) != len(set(collected_ids))):
        collection["status"] = "invalid"
        collection["reason"] = "Could not determine unique collected test identities."
    results = []
    if collection["status"] == "passed":
        with concurrent.futures.ThreadPoolExecutor(max_workers=4) as executor:
            results = list(executor.map(run_module, paths))
    unchanged = (verification_inputs() == fingerprints
                 and dependency_inputs(dependency_roots) == dependency_hashes)
    versions = {}
    for package in ("pytest", "numpy", "scipy", "numba", "coverage", "ruff"):
        try:
            versions[package] = importlib.metadata.version(package)
        except importlib.metadata.PackageNotFoundError:
            versions[package] = "not_installed"
    report = {"module_timeout_seconds": MODULE_TIMEOUT_SECONDS, "results": results,
              "collection": collection,
              "expected_modules": [path.relative_to(ROOT).as_posix() for path in paths],
              "diagnostic_snapshot": "20s best-effort Python timer callback with GIL; fatal handler retained",
              "resource_isolated_modules": [],
              "other_module_workers": 4,
              "runtime_performance_under_concurrent_load": "not_assessed",
              "deferred": list(DEFERRED),
              "deferred_reason": "previous runs did not complete; not counted as passed",
              "source_sha256": fingerprints,
              "fingerprint_scope": "sources, tests, scripts, documents, examples and package data",
              "external_dependency_sha256": dependency_hashes,
              "external_dependency_roots": [str(path) for path in dependency_roots],
              "runtime": {"python": sys.version, "packages": versions},
              "collected_test_ids": collected_ids,
              "explicit_sidera_test_path": None if options.sidera_path is None
              else str(explicit_path),
              "source_unchanged_during_run": unchanged,
              "complete_suite_passed": collection["status"] == "passed" and unchanged
              and not DEFERRED and all(
                  r["status"] == "passed" and r.get("junit", {}).get("skipped", 1) == 0
                  for r in results),
              "environment": "Windows; explicit workspace NUMBA_CACHE_DIR; recorded host permissions"}
    (OUT / "summary.json").write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"collection_status": collection["status"], "modules": len(results),
        "expected_modules": len(paths), "passed_modules": sum(
        r["status"] == "passed" for r in results), "open_modules": [
            r["module"] for r in results if r["status"] != "passed"]}))
    return int(not report["complete_suite_passed"])


if __name__ == "__main__":
    raise SystemExit(main())
