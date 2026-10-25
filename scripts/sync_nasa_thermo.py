"""Extract the bundled NASA CEA thermodynamic subset from a pinned ``thermo.inp``.

The local equilibrium provider needs more product species than the ten
NASA/TP-2002-211556 records in ``nasa_glenn_2002.json`` and needs the
assigned-enthalpy condensed reactants (liquid oxygen, liquid hydrogen, RP-1,
NTO, MMH, ...) that NASA CEA uses. Both come from NASA's own ``thermo.inp``
(https://github.com/nasa/cea, Apache-2.0), read in the fixed-column format
documented in NASA RP-1311 Part II, Appendix A [1]_ and NASA/TP-2002-211556 [2]_.

Usage::

    python scripts/sync_nasa_thermo.py --source path/to/thermo.inp
    python scripts/sync_nasa_thermo.py --source path/to/thermo.inp --check

``--check`` re-extracts and fails when the bundled JSON differs.

References
----------
.. [1] NASA RP-1311, "Computer Program for Calculation of Complex Chemical
       Equilibrium Compositions and Applications, Part II: Users Manual and
       Program Description". https://ntrs.nasa.gov/citations/19960044559
.. [2] NASA/TP-2002-211556. NASA Glenn, "NASA Glenn Coefficients for
       Calculating Thermodynamic Properties of Individual Species".
       https://ntrs.nasa.gov/citations/20020085330
"""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
TARGET = ROOT / "src" / "rocket_propulsion" / "thermochemistry" / "data" / "nasa_cea_thermo.json"
REPOSITORY = "https://github.com/nasa/cea"
COMMIT = "2f79a647e85737742fb08a127da00d1ea1d069da"
EXPECTED_SHA256 = "fa7746572952d74e249e818a82a35c113829742fb421a308e167185528884363"

# Gas-phase products added to the ten NASA/TP-2002-211556 species. Names are
# the thermo.inp names; the bundled key is the same name.
PRODUCTS = (
    "HO2", "H2O2", "O3",
    "N", "NO", "NO2", "N2O", "NH", "NH2", "NH3", "HNO", "N2H2",
    "C", "CH", "CH2", "CH3", "HCO", "HCHO,formaldehy", "COOH", "HCOOH",
    "C2H2,acetylene", "C2H4", "HCN", "CN", "HNCO", "NCO",
)

# Condensed reactants with an assigned enthalpy at one temperature.
REACTANTS = (
    "H2(L)", "O2(L)", "CH4(L)", "RP-1", "N2O4(L)", "CH6N2(L)",
)

DESCRIPTIONS = {
    "H2(L)": "Liquid hydrogen", "O2(L)": "Liquid oxygen", "CH4(L)": "Liquid methane",
    "RP-1": "RP-1 kerosene", "N2O4(L)": "Nitrogen tetroxide",
    "CH6N2(L)": "Monomethylhydrazine",
}


def _number(text: str) -> float:
    """Parse a Fortran real, accepting the ``D`` exponent marker; blank reads as zero."""
    text = text.strip()
    return float(text.replace("D", "E").replace("d", "e")) if text else 0.0


def _formula(line: str) -> dict[str, float]:
    """Read the five (symbol, count) fields in columns 11-50 of record line 2."""
    formula: dict[str, float] = {}
    for start in range(10, 50, 8):
        symbol = line[start:start + 2].strip()
        count = line[start + 2:start + 8].strip()
        if symbol and count and float(count) != 0.0:
            formula[symbol.capitalize() if len(symbol) == 2 else symbol] = float(count)
    return formula


def parse_thermo_inp(text: str) -> tuple[dict[str, dict], dict[str, dict]]:
    """Return ``(products, reactants)`` records keyed by thermo.inp name."""
    lines = text.splitlines()
    start = next(index for index, line in enumerate(lines) if line.startswith("thermo")) + 2
    products: dict[str, dict] = {}
    reactants: dict[str, dict] = {}
    section = products
    index = start
    while index < len(lines):
        line = lines[index]
        if line.startswith("END PRODUCTS"):
            section = reactants
            index += 1
            continue
        if line.startswith("END REACTANTS"):
            break
        if not line.strip() or line.startswith("!"):
            index += 1
            continue
        name = line[:18].strip()
        comment = line[18:].strip()
        header = lines[index + 1]
        intervals_count = int(header[0:2])
        record = {
            "name": name,
            "comment": comment,
            "reference_date_code": header[3:9].strip(),
            "formula": _formula(header),
            "phase": int(header[50:52]),
            "molecular_mass_g_mol": _number(header[52:65]),
            "heat_of_formation_j_mol": _number(header[65:80]),
        }
        if intervals_count == 0:
            assigned = lines[index + 2]
            record["assigned_temperature_k"] = _number(assigned[0:11])
            record["intervals"] = []
            index += 3
        else:
            intervals = []
            cursor = index + 2
            for _ in range(intervals_count):
                bounds, first, second = lines[cursor], lines[cursor + 1], lines[cursor + 2]
                values_first = [_number(first[i:i + 16]) for i in range(0, 80, 16)]
                values_second = [_number(second[i:i + 16]) for i in range(0, 80, 16)]
                intervals.append({
                    "minimum_temperature_k": _number(bounds[0:11]),
                    "maximum_temperature_k": _number(bounds[11:21]),
                    "coefficients": values_first + values_second[0:2] + values_second[3:5],
                })
                cursor += 3
            record["intervals"] = intervals
            index = cursor
        section[name] = record
    return products, reactants


def build_document(text: str, sha256: str) -> dict:
    """Assemble the bundled JSON document from parsed thermo.inp records."""
    products, reactants = parse_thermo_inp(text)
    missing = [name for name in PRODUCTS if name not in products]
    missing += [name for name in REACTANTS if name not in reactants]
    if missing:
        raise SystemExit(f"thermo.inp lacks required records: {missing}")
    species = []
    for name in PRODUCTS:
        record = products[name]
        if record["phase"] != 0:
            raise SystemExit(f"{name} is not a gas-phase record")
        species.append({
            "key": name,
            "name": name,
            "thermo_inp_comment": record["comment"],
            "reference_date_code": record["reference_date_code"],
            "formula": record["formula"],
            "molecular_mass_g_mol": record["molecular_mass_g_mol"],
            "heat_of_formation_j_mol_298": record["heat_of_formation_j_mol"],
            "intervals": record["intervals"],
        })
    condensed = []
    for name in REACTANTS:
        record = reactants[name]
        if record["intervals"]:
            raise SystemExit(f"{name} is not an assigned-enthalpy reactant")
        condensed.append({
            "key": name,
            "name": DESCRIPTIONS[name],
            "thermo_inp_comment": record["comment"],
            "reference_date_code": record["reference_date_code"],
            "formula": record["formula"],
            "molecular_mass_g_mol": record["molecular_mass_g_mol"],
            "assigned_enthalpy_j_mol": record["heat_of_formation_j_mol"],
            "assigned_temperature_k": record["assigned_temperature_k"],
            "phase": "condensed",
        })
    return {
        "schema_version": 1,
        "source": {
            "title": "NASA CEA thermodynamic database (thermo.inp)",
            "authors": "NASA Glenn Research Center",
            "url": f"{REPOSITORY}/blob/{COMMIT}/data/thermo.inp",
            "repository": REPOSITORY,
            "commit": COMMIT,
            "sha256": sha256,
            "license": "Apache-2.0",
            "format_reference": "NASA RP-1311 Part II Appendix A; NASA/TP-2002-211556",
        },
        "model": {
            "name": "NASA Glenn seven-term cp polynomial with two integration constants",
            "coefficient_order": ["a1", "a2", "a3", "a4", "a5", "a6", "a7", "b1", "b2"],
            "reference_pressure_pa": 100000,
        },
        "species": species,
        "reactants": condensed,
    }


def render(document: dict) -> str:
    """Serialize deterministically with a trailing newline."""
    return json.dumps(document, indent=1, ensure_ascii=False) + "\n"


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--source", required=True, type=Path, help="path to thermo.inp")
    parser.add_argument("--check", action="store_true", help="fail when the bundle differs")
    arguments = parser.parse_args()
    raw = arguments.source.read_bytes()
    digest = hashlib.sha256(raw).hexdigest()
    if digest != EXPECTED_SHA256:
        print(f"thermo.inp sha256 {digest} differs from pinned {EXPECTED_SHA256}")
        return 1
    content = render(build_document(raw.decode("utf-8"), digest))
    if arguments.check:
        current = TARGET.read_text(encoding="utf-8") if TARGET.is_file() else ""
        if current != content:
            print(f"{TARGET} is out of date")
            return 1
        print("bundled NASA CEA subset matches thermo.inp")
        return 0
    TARGET.write_text(content, encoding="utf-8", newline="\n")
    print(f"wrote {TARGET} ({len(content)} bytes)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
