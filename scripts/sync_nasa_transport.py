"""Extract the bundled NASA CEA transport-property coefficients from a pinned ``trans.inp``.

NASA's CEA distributes least-squares coefficients of viscosity and thermal
conductivity for individual gases, and of the viscosity interaction parameter
for some binary pairs, in ``data/trans.inp`` (https://github.com/nasa/cea,
Apache-2.0). Each fit has the form of RP-1311 Part I, equation (5.1) [1]_::

    ln(eta) or ln(lambda) = A ln T + B / T + C / T**2 + D

with ``eta`` in micropoise and ``lambda`` in microwatt per centimetre-kelvin.
The fixed-column layout read here is the one of CEA's own reader
(``source/database_compile.f90``, subroutine ``utran``, formats ``99001`` and
``99002``) at the same commit [2]_: a header ``(2a16, 2x, a1, i1, a1, i1)``
with one or two species names and the number of viscosity (``V``) and
conductivity (``C``) intervals, then one line ``(1x, a1, 2f9.2, 4e15.8)`` per
interval.

Usage::

    python scripts/sync_nasa_transport.py --source path/to/trans.inp
    python scripts/sync_nasa_transport.py --source path/to/trans.inp --check

``--check`` re-extracts and fails when the bundled JSON differs.

References
----------
.. [1] NASA RP-1311, "Computer Program for Calculation of Complex Chemical
       Equilibrium Compositions and Applications, Part I: Analysis", section
       5.1, equation (5.1). https://ntrs.nasa.gov/citations/19950013764
.. [2] NASA-CEA-TRANS. NASA Glenn Research Center, CEA transport database
       ``data/trans.inp`` and reader ``source/database_compile.f90``, commit
       2f79a647e85737742fb08a127da00d1ea1d069da. https://github.com/nasa/cea
"""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
TARGET = ROOT / "src" / "rocket_propulsion" / "thermochemistry" / "data" / "nasa_cea_transport.json"
REPOSITORY = "https://github.com/nasa/cea"
COMMIT = "2f79a647e85737742fb08a127da00d1ea1d069da"
EXPECTED_SHA256 = "379c514a7f7638371d8d9254a1d653cd078084a9371061aa0b45c5a5c3334d41"


def _fortran_real(field: str) -> float:
    """Parse an ``e15.8``/``f9.2`` field; Fortran ignores blanks (``0.61E 00``)."""
    text = field.replace(" ", "").replace("D", "E").replace("d", "e")
    if not text:
        raise SystemExit(f"blank numeric field {field!r}")
    return float(text)


def _interval(line: str) -> tuple[str, dict]:
    """Read one ``(1x, a1, 2f9.2, 4e15.8)`` coefficient line."""
    kind = line[1:2]
    if kind not in {"V", "C"}:
        raise SystemExit(f"unexpected interval line {line!r}")
    record = {
        "minimum_temperature_k": _fortran_real(line[2:11]),
        "maximum_temperature_k": _fortran_real(line[11:20]),
        "coefficients": [_fortran_real(line[20 + 15 * i:35 + 15 * i]) for i in range(4)],
    }
    return kind, record


def parse_trans_inp(text: str) -> list[dict]:
    """Return every record of ``trans.inp`` in file order."""
    lines = text.splitlines()
    start = next(i for i, line in enumerate(lines) if line.lower().startswith("tran")) + 1
    records = []
    index = start
    while index < len(lines):
        header = lines[index]
        first, second = header[0:16].strip(), header[16:32].strip()
        if first in {"end", "LAST"}:
            break
        if header[34:35] != "V" or header[36:37] != "C":
            raise SystemExit(f"malformed header {header!r}")
        viscosity_count, conductivity_count = int(header[35:36]), int(header[37:38])
        viscosity, conductivity = [], []
        for offset in range(viscosity_count + conductivity_count):
            kind, interval = _interval(lines[index + 1 + offset])
            (viscosity if kind == "V" else conductivity).append(interval)
        if len(viscosity) != viscosity_count or len(conductivity) != conductivity_count:
            raise SystemExit(f"interval count mismatch for {first} {second}")
        records.append({
            "species": [first] if not second else [first, second],
            "comment": header[38:].strip(),
            "viscosity": viscosity,
            "thermal_conductivity": conductivity,
        })
        index += 1 + viscosity_count + conductivity_count
    return records


def build_document(text: str, sha256: str) -> dict:
    """Assemble the bundled JSON document from parsed trans.inp records."""
    records = parse_trans_inp(text)
    species = []
    binaries = []
    for record in records:
        if len(record["species"]) == 1:
            species.append({
                "name": record["species"][0],
                "trans_inp_comment": record["comment"],
                "viscosity": record["viscosity"],
                "thermal_conductivity": record["thermal_conductivity"],
            })
        else:
            binaries.append({
                "pair": record["species"],
                "trans_inp_comment": record["comment"],
                "viscosity": record["viscosity"],
                "thermal_conductivity": record["thermal_conductivity"],
            })
    return {
        "schema_version": 1,
        "source": {
            "title": "NASA CEA transport property database (trans.inp)",
            "authors": "NASA Glenn Research Center; data from Svehla (1995) and the per-record sources",
            "url": f"{REPOSITORY}/blob/{COMMIT}/data/trans.inp",
            "repository": REPOSITORY,
            "commit": COMMIT,
            "sha256": sha256,
            "license": "Apache-2.0",
            "format_reference": (
                "NASA RP-1311 Part I eq. (5.1); CEA source/database_compile.f90 "
                "subroutine utran, formats 99001 and 99002"
            ),
        },
        "model": {
            "name": "ln(property) = A ln T + B/T + C/T^2 + D",
            "coefficient_order": ["A", "B", "C", "D"],
            "viscosity_unit": "micropoise",
            "thermal_conductivity_unit": "microwatt/(cm K)",
        },
        "species": species,
        "binary_interactions": binaries,
    }


def render(document: dict) -> str:
    """Serialize deterministically with a trailing newline."""
    return json.dumps(document, indent=1, ensure_ascii=False) + "\n"


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--source", required=True, type=Path, help="path to trans.inp")
    parser.add_argument("--check", action="store_true", help="fail when the bundle differs")
    arguments = parser.parse_args()
    raw = arguments.source.read_bytes()
    digest = hashlib.sha256(raw).hexdigest()
    if digest != EXPECTED_SHA256:
        print(f"trans.inp sha256 {digest} differs from pinned {EXPECTED_SHA256}")
        return 1
    content = render(build_document(raw.decode("utf-8"), digest))
    if arguments.check:
        current = TARGET.read_text(encoding="utf-8") if TARGET.is_file() else ""
        if current != content:
            print(f"{TARGET} is out of date")
            return 1
        print("bundled NASA CEA transport data matches trans.inp")
        return 0
    TARGET.write_text(content, encoding="utf-8", newline="\n")
    print(f"wrote {TARGET} ({len(content)} bytes)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
