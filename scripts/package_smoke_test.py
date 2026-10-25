"""Release smoke test: build the wheel, install it into an empty venv and use it.

The Phase 7 release gate asks for installation, a first calculation, saving and
reopening a workspace and exporting on a clean machine. This script runs the
same steps in a freshly created virtual environment (a clean *environment*,
not a clean machine) from outside the source tree, so nothing is imported
from ``src/``:

1. build the wheel with build isolation;
2. create an empty venv and install only that wheel;
3. solve a CEA-gated rocket point and an ideal nozzle;
4. reopen ``examples/two_case.rplab.json`` and check the save/open round trip;
5. produce CSV/SVG/DXF/STL geometry exports and a transient result;
6. start the ``rocket-propulsion`` console script and fetch the web shell.

Usage::

    python scripts/package_smoke_test.py [--keep]
"""

from __future__ import annotations

import argparse
import json
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]

CHECKS = r'''
import json, pathlib, subprocess, sys, time, urllib.request
import rocket_propulsion
from rocket_propulsion.api.routes import dispatch_calculation
from rocket_propulsion.studies import workspace_from_json, workspace_to_json
assert "site-packages" in rocket_propulsion.__file__, rocket_propulsion.__file__
record = {"version": rocket_propulsion.__version__}
result = dispatch_calculation("/api/v1/rocket-performance", {"propellant_pair": "LOX_RP1", "area_ratios": [16]})
record["lox_rp1_vacuum_isp_s"] = result["result"]["exits"][0]["vacuum_specific_impulse_s"]
record["nozzle_thrust_n"] = dispatch_calculation("/api/v1/nozzle", {
    "chamber_pressure_pa": 7e6, "chamber_temperature_k": 3500, "throat_area_m2": 0.01, "area_ratio": 25,
    "ambient_pressure_pa": 0, "gamma": 1.22, "gas_constant_j_kg_k": 355})["thrust_n"]
saved = workspace_to_json(workspace_from_json(pathlib.Path(sys.argv[1]).read_text(encoding="utf-8")))
assert workspace_to_json(workspace_from_json(saved)) == saved
record["workspace_round_trip"] = True
exports = dispatch_calculation("/api/v1/nozzle-contour", {
    "contour": "rao", "throat_area_m2": 0.01, "area_ratio": 25, "gamma": 1.22, "half_angle_deg": 15,
    "chamber_radius_ratio": 2.5, "bell_length_fraction": 0.8, "characteristic_count": 20,
    "throat_angle_deg": 30, "exit_angle_deg": 8, "station_count": 61, "axisymmetric_correction": True,
    "include_exports": True, "include_comparison": False, "angular_segments": 16})["exports"]
assert all(exports[key] for key in ("csv", "svg", "dxf", "stl"))
record["exports"] = ["csv", "svg", "dxf", "stl"]
record["blowdown_duration_s"] = dispatch_calculation("/api/v1/transient/blowdown", {})["summary"]["duration_s"]
script = pathlib.Path(sys.executable).with_name("rocket-propulsion.exe" if sys.platform == "win32" else "rocket-propulsion")
server = subprocess.Popen([str(script), "--port", sys.argv[2]], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
try:
    base = f"http://127.0.0.1:{sys.argv[2]}"
    for _ in range(100):
        try:
            record["health"] = json.load(urllib.request.urlopen(base + "/api/v1/health", timeout=1))["status"]
            break
        except OSError:
            time.sleep(0.2)
    page = urllib.request.urlopen(base + "/", timeout=5).read().decode("utf-8")
    record["web_shell"] = 'id="performance"' in page and 'id="feed-transient"' in page
finally:
    server.terminate()
    server.wait(timeout=10)
assert record.get("health") == "ok" and record["web_shell"], record
print(json.dumps(record))
'''


def run(command: list[str], **options) -> subprocess.CompletedProcess[str]:
    return subprocess.run(command, check=True, text=True, capture_output=True, **options)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--keep", action="store_true", help="keep the temporary build directory")
    parser.add_argument("--port", default="8799", help="port for the console-script check")
    arguments = parser.parse_args()
    work = Path(tempfile.mkdtemp(prefix="rplab-smoke-"))
    try:
        run([sys.executable, "-m", "pip", "wheel", str(ROOT), "--no-deps", "-w", str(work / "wheel"), "-q"])
        wheel = next((work / "wheel").glob("*.whl"))
        run([sys.executable, "-m", "venv", str(work / "venv")])
        python = work / "venv" / ("Scripts" if sys.platform == "win32" else "bin") / "python"
        run([str(python), "-m", "pip", "install", "-q", str(wheel)])
        (work / "checks.py").write_text(CHECKS, encoding="utf-8")
        result = run([str(python), "-I", str(work / "checks.py"),
                      str(ROOT / "examples" / "two_case.rplab.json"), arguments.port], cwd=work)
        record = json.loads(result.stdout.strip().splitlines()[-1])
        record["wheel"] = wheel.name
        print(json.dumps(record, indent=1))
        return 0
    except subprocess.CalledProcessError as error:
        print(error.stdout, error.stderr, sep="\n")
        return 1
    finally:
        if arguments.keep:
            print(f"kept {work}")
        else:
            shutil.rmtree(work, ignore_errors=True)


if __name__ == "__main__":
    sys.exit(main())
