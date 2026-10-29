"""Check selected numeric faults in isolated copies, with before/after controls.

The case file describes implementation or oracle faults, not physical validation.
Timeouts, crashes, skips, collection errors and surviving faults are open results.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import re
import shutil
import subprocess
import sys
import time
import xml.etree.ElementTree as ET
from pathlib import Path
from zipfile import ZIP_DEFLATED, ZipFile

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from scripts.run_bounded_regression import verification_inputs

CASE_FILE = ROOT / 'tests/fixtures/numeric_admission_mutations.json'
BOOTSTRAP = '''
import json,sys
from pathlib import Path
sandbox=Path(sys.argv[1]).resolve()
sys.path[:0]=[str(sandbox/'src'),str(sandbox)]
import rocket_propulsion, tests
assert Path(rocket_propulsion.__file__).resolve().is_relative_to(sandbox/'src')
test_paths=[Path(p).resolve() for p in tests.__path__]
assert test_paths and all(p == sandbox/'tests' for p in test_paths)
print(json.dumps({'loaded_source':rocket_propulsion.__file__,
                  'loaded_tests':[str(p) for p in test_paths]}),flush=True)
import pytest
import _pytest.unittest
_pytest.unittest.TestCaseFunction.failfast=True
raise SystemExit(pytest.main(sys.argv[2:]))
'''


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def copied_inventory(sandbox):
    return {p.relative_to(sandbox).as_posix(): hashlib.sha256(p.read_bytes()).hexdigest()
            for name in ('src', 'tests') for p in (sandbox / name).rglob('*')
            if p.is_file() and '__pycache__' not in p.parts}


def archive_baseline(out, sandbox, inventory):
    archive = out / 'source_test_baseline.zip'
    with ZipFile(archive, 'w', compression=ZIP_DEFLATED, compresslevel=6) as z:
        for relative in inventory:
            z.write(sandbox / relative, relative)
        z.writestr('inventory.json', json.dumps(inventory, indent=2))
    with ZipFile(archive) as z:
        require(z.testzip() is None, 'Baseline archive CRC failed')
        for relative, digest in inventory.items():
            require(hashlib.sha256(z.read(relative)).hexdigest() == digest, relative)
    return hashlib.sha256(archive.read_bytes()).hexdigest()


def run_test(out, sandbox, selector, tag):
    xml = out / (tag + '.xml')
    filename, _, suffix = selector.partition('::')
    target = str(sandbox / filename) + ('::' + suffix if suffix else '')
    args = [sys.executable, '-c', BOOTSTRAP, str(sandbox), '-c', str(ROOT / 'pyproject.toml'),
            '-o', 'pythonpath=' + str(sandbox / 'src'), '--rootdir=' + str(sandbox), target,
            '-q', '--tb=short', '--maxfail=1', '--junitxml=' + str(xml)]
    env = {**os.environ, 'PYTHONPATH': str(sandbox / 'src') + os.pathsep + str(sandbox),
           'PYTHONDONTWRITEBYTECODE': '1', 'OPENBLAS_NUM_THREADS': '1', 'OMP_NUM_THREADS': '1'}
    started = time.monotonic()
    try:
        process = subprocess.run(args, cwd=sandbox, env=env, capture_output=True,
                                 timeout=90, check=False)
    except subprocess.TimeoutExpired as exc:
        (out / (tag + '.txt')).write_bytes((exc.stdout or b'') + (exc.stderr or b''))
        return {'status': 'timeout', 'seconds': time.monotonic() - started}
    (out / (tag + '.txt')).write_bytes(process.stdout + process.stderr)
    tree = ET.parse(xml).getroot() if xml.exists() else ET.Element('missing_junit')
    suites = tree.findall('testsuite')
    counts = {key: sum(int(s.attrib.get(key, 0)) for s in suites)
              for key in ('tests', 'failures', 'errors', 'skipped')}
    failures = list(tree.iter('failure'))
    assertion_only = bool(failures) and all(
        item.attrib.get('message', '').startswith(('AssertionError', 'assert '))
        or re.search(r'^E\s+AssertionError\b', item.text or '', re.MULTILINE)
        for item in failures)
    return {'return_code': process.returncode, 'junit': counts,
            'assertion_failures_only': bool(assertion_only),
            'seconds': round(time.monotonic() - started, 2),
            'source_origin_assertion': process.stdout.decode(errors='replace').splitlines()[0]
            if process.stdout else ''}


def passing_control(result):
    c = result.get('junit', {})
    return (result.get('return_code') == 0 and c.get('tests', 0) > 0
            and all(c.get(k) == 0 for k in ('failures', 'errors', 'skipped')))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output-dir', type=Path,
                        default=ROOT / 'output/assurance/numeric_mutations')
    parser.add_argument('--case', action='append', dest='names', help='Select a named fault; repeatable')
    parser.add_argument('--list', action='store_true', help='Preflight and list cases without running')
    args = parser.parse_args()
    cases = json.loads(CASE_FILE.read_text(encoding='utf-8'))
    names = [row[0] for row in cases]
    require(len(names) == len(set(names)), 'Duplicate fault name')
    for name, relative, original, changed, selector in cases:
        require(bool(re.fullmatch('[a-z0-9_]+', name)), 'Unsafe fault name')
        source = (ROOT / relative).resolve()
        require(source.is_relative_to(ROOT / 'src') or source.is_relative_to(ROOT / 'tests'),
                'Mutation target outside source/test tree')
        target = (ROOT / selector.split('::')[0]).resolve()
        require(target.is_relative_to(ROOT / 'tests') and target.is_file(), 'Invalid test target')
        text = source.read_text(encoding='utf-8')
        require(text.count(original) == 1 and original != changed, 'Stale fault anchor: ' + name)
        compile(text.replace(original, changed), relative, 'exec')
    require(not args.names or set(args.names) <= set(names), 'Unknown selected fault')
    cases = [c for c in cases if not args.names or c[0] in args.names]
    if args.list:
        print(json.dumps([{'name': c[0], 'file': c[1], 'selector': c[4]} for c in cases], indent=2))
        return 0
    out = args.output_dir.resolve()
    require(out.is_relative_to(ROOT / 'output'), 'Output directory must be inside checkout/output')
    out.mkdir(parents=True, exist_ok=False)
    before = verification_inputs(ROOT)
    result = {'complete': False, 'selected_faults': len(cases), 'results': [],
              'source_sha256': before,
              'case_file_sha256': hashlib.sha256(CASE_FILE.read_bytes()).hexdigest(),
              'scope': 'Selected equation, refusal, ownership and test-oracle faults; '
                       'not exhaustive mutation coverage or hardware qualification'}
    snapshot_inventory = None
    try:
        for name, relative, original, changed, selector in cases:
            sandbox = out / name
            for directory in ('src', 'tests'):
                shutil.copytree(ROOT / directory, sandbox / directory,
                                ignore=shutil.ignore_patterns('__pycache__'))
            marker = sandbox / 'tests/__init__.py'
            require(not marker.exists(), 'Test namespace changed; review isolation bootstrap')
            marker.write_text('"""Isolated test namespace marker; no test logic."""\n')
            inventory = copied_inventory(sandbox)
            if snapshot_inventory is None:
                snapshot_inventory = inventory
                result['baseline_archive_sha256'] = archive_baseline(out, sandbox, inventory)
            require(inventory == snapshot_inventory, 'Source/test inputs changed between faults')
            path = sandbox / relative
            original_bytes = path.read_bytes()
            text = path.read_text(encoding='utf-8')
            baseline = run_test(out, sandbox, selector, name + '_before')
            require(passing_control(baseline), 'Original control failed: ' + name)
            mutated = text.replace(original, changed)
            try:
                path.write_text(mutated, encoding='utf-8', newline='')
                attack = run_test(out, sandbox, selector, name + '_mutant')
            finally:
                path.write_bytes(original_bytes)
            control = run_test(out, sandbox, selector, name + '_after')
            counts = attack.get('junit', {})
            killed = (attack.get('return_code') == 1 and counts.get('failures', 0) > 0
                      and counts.get('errors') == counts.get('skipped') == 0
                      and attack.get('assertion_failures_only') and passing_control(control)
                      and baseline['junit']['tests'] == control['junit']['tests'])
            row = {'name': name, 'file': relative, 'selector': selector, 'before': baseline,
                   'mutant': attack, 'after': control, 'killed': bool(killed),
                   'original_sha256': hashlib.sha256(original_bytes).hexdigest(),
                   'restored_sha256': hashlib.sha256(path.read_bytes()).hexdigest(),
                   'mutated_sha256': hashlib.sha256(mutated.encode('utf-8')).hexdigest()}
            result['results'].append(row)
            print(name, 'killed' if killed else 'OPEN', flush=True)
            require(copied_inventory(sandbox) == inventory, 'Copied files not restored byte-exactly')
            require(killed, 'Fault survived or control failed: ' + name)
            for directory in ('src', 'tests'):
                target = (sandbox / directory).resolve()
                require(target.parent == sandbox.resolve() and target.is_relative_to(out)
                        and not target.is_symlink(), 'Unsafe temporary tree cleanup')
                shutil.rmtree(target)
            (out / 'results.json').write_text(json.dumps(result, indent=2) + '\n')
        result.update(complete=True, killed=len(cases), passing_original_controls=2 * len(cases))
        return 0
    finally:
        result['production_inputs_unchanged'] = before == verification_inputs(ROOT)
        (out / 'results.json').write_text(json.dumps(result, indent=2) + '\n')
        require(result['production_inputs_unchanged'], 'Production inputs changed during run')
        print(json.dumps({k: v for k, v in result.items()
                          if k not in ('results', 'source_sha256')}), flush=True)


if __name__ == '__main__':
    raise SystemExit(main())
