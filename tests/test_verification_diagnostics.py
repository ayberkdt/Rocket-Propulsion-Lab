"""Exercise the real pytest diagnostic hook without a C watchdog."""

import os
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


class VerificationDiagnosticTests(unittest.TestCase):
    def test_python_timer_records_blocked_test_and_keeps_its_result(self):
        with tempfile.TemporaryDirectory(dir=ROOT) as folder:
            fixture = Path(folder) / "test_snapshot_fixture.py"
            fixture.write_text("import time\ndef test_wait_for_snapshot():\n    time.sleep(.25)\n"
                               "    assert 6 * 7 == 42\n", encoding="utf-8")
            env = dict(os.environ, PYTHONPATH=str(ROOT))
            command = [sys.executable, "-m", "pytest", str(fixture), "-q", "-o",
                       "faulthandler_timeout=0", "-p", "scripts.pytest_diagnostics",
                       "--diagnostic-timeout=.02"]
            child = subprocess.run(command, cwd=ROOT, env=env, capture_output=True, timeout=30)
            log = (child.stdout + child.stderr).decode(errors="replace")
            self.assertEqual(child.returncode, 0, log)
            self.assertIn("Diagnostic snapshot after", log)
            self.assertIn("Python thread holds GIL", log)
            self.assertIn("test_wait_for_snapshot", log)
