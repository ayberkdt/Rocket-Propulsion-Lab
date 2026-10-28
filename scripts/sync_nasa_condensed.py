"""Extract the NASA CEA condensed-phase product subset from a pinned ``thermo.inp``.

The condensed-phase equilibrium provider
(:mod:`rocket_propulsion.thermochemistry.condensed_equilibrium`) needs, in
addition to the gas species already bundled in ``nasa_glenn_2002.json`` and
``nasa_cea_thermo.json``:

* every condensed product of NASA's ``thermo.inp`` made only of the elements
  Al, Cl, C, H, O and N (graphite, ice, liquid water, aluminium, alumina,
  aluminium chloride, aluminium nitride, ...);
* every non-ionic gas product of the same element set that is not bundled
  yet (Al, Cl, Al-Cl-H-O-N compounds and the C/H/O/N gases that NASA CEA also
  considers for those elements);
* two reactants given by their enthalpy at 298.15 K: ``AL(cr)`` and
  ``NH4CLO4(I)``.

The records are read from NASA's ``thermo.inp`` (https://github.com/nasa/cea,
Apache-2.0) with the fixed-column reader of ``scripts/sync_nasa_thermo.py``
(format of NASA RP-1311 Part II, Appendix A [1]_). Before anything is written
the script checks, for every record, three identities and stores the measured
residuals in the JSON file:

1. heat of formation: ``R T (H/RT)(298.15 K)`` from the polynomial equals the
   header heat of formation for records whose first interval covers 298.15 K
   or starts at 300 K, with ``H(298.15) = Delta_f H(298.15)`` by definition
   (RP-1311 Part I, eqs. (4.1)-(4.2) [2]_) and R = 8314.51 J/(kmol K) as in
   RP-1311 Part I section 2.1 [2]_;
2. phase transitions: the Gibbs energies of two consecutive condensed phases
   of one substance are equal at their shared transition temperature, and the
   enthalpy of transition is positive;
3. continuity of ``H/RT`` and ``S/R`` across interval boundaries (``cp/R``
   may step between separately fitted intervals; its largest step is
   recorded). One record, ``ALN(L)``, carries a step in ``H`` and ``S`` at
   2700 K in the source data; it is recorded, not altered.

Usage::

    python scripts/sync_nasa_condensed.py --source path/to/thermo.inp
    python scripts/sync_nasa_condensed.py --source path/to/thermo.inp --check

``--check`` re-extracts and fails when the bundled JSON differs.

References
----------
.. [1] NASA RP-1311, "Computer Program for Calculation of Complex Chemical
       Equilibrium Compositions and Applications, Part II: Users Manual and
       Program Description". https://ntrs.nasa.gov/citations/19960044559
.. [2] NASA-RP1311-P1-CONDENSED. Gordon and McBride, NASA RP-1311 Part I
       (1994). https://ntrs.nasa.gov/citations/19950013764
"""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
from math import log
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import itertools

from sync_nasa_thermo import (
    COMMIT,
    EXPECTED_SHA256,
    REPOSITORY,
    parse_thermo_inp,
)

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "src" / "rocket_propulsion" / "thermochemistry" / "data"
TARGET = DATA / "nasa_cea_condensed.json"
EXISTING = (DATA / "nasa_glenn_2002.json", DATA / "nasa_cea_thermo.json")

ELEMENTS = ("Al", "C", "Cl", "H", "N", "O")
# Counts fixed by the pinned thermo.inp; a different count means the selection
# rule no longer selects what was reviewed.
EXPECTED_GAS_COUNT = 186
EXPECTED_CONDENSED_COUNT = 16
REACTANTS = {
    "AL(cr)": "Aluminium (cubic) reactant at 298.15 K",
    "NH4CLO4(I)": "Ammonium perchlorate (I) reactant at 298.15 K",
}
REACTANT_TEMPERATURE_K = 298.15
# RP-1311 Part I section 2.1: R = 8314.51 J/(kg-mole K); used only for the
# heat-of-formation identity, because the header enthalpies are in J/mol.
CEA_GAS_CONSTANT_J_MOL_K = 8.31451
HEAT_OF_FORMATION_TOLERANCE_J_MOL = 0.5
TRANSITION_GIBBS_TOLERANCE = 1.0e-3
# Absolute step in H/RT or S/R; the largest measured in the pinned file is
# 1.24e-3 (ALOCL at 1000 K), the next 5.2e-5 (H2O(L) at 373.15 K).
CONTINUITY_TOLERANCE = 2.0e-3
# Steps in H/RT and S/R present in the source data itself (record, boundary K).
KNOWN_DISCONTINUITIES = {("ALN(L)", 2700.0)}


def _h_over_rt(coefficients: list[float], temperature: float) -> float:
    a1, a2, a3, a4, a5, a6, a7, b1, _ = coefficients
    t = temperature
    return (-a1 / t**2 + a2 * log(t) / t + a3 + a4 * t / 2.0 + a5 * t**2 / 3.0
            + a6 * t**3 / 4.0 + a7 * t**4 / 5.0 + b1 / t)


def _s_over_r(coefficients: list[float], temperature: float) -> float:
    a1, a2, a3, a4, a5, a6, a7, _, b2 = coefficients
    t = temperature
    return (-a1 / t**2 / 2.0 - a2 / t + a3 * log(t) + a4 * t + a5 * t**2 / 2.0
            + a6 * t**3 / 3.0 + a7 * t**4 / 4.0 + b2)


def _cp_over_r(coefficients: list[float], temperature: float) -> float:
    a1, a2, a3, a4, a5, a6, a7, _, _ = coefficients
    t = temperature
    return a1 / t**2 + a2 / t + a3 + a4 * t + a5 * t**2 + a6 * t**3 + a7 * t**4


def _phase_group(name: str) -> str:
    """Return the substance of a condensed record: ``AL2O3(a)`` -> ``AL2O3``."""
    return name[: name.rfind("(")] if name.endswith(")") else name


def _existing_keys() -> set[str]:
    keys: set[str] = set()
    for path in EXISTING:
        document = json.loads(path.read_text(encoding="utf-8"))
        keys.update(record["key"].upper() for record in document["species"])
    return keys


def _selected(record: dict) -> bool:
    formula = record["formula"]
    return "E" not in formula and set(formula) <= set(ELEMENTS)


def _heat_of_formation_residual(record: dict) -> float | None:
    """Return ``H_poly(298.15) - Delta_f H`` in J/mol, or None if not applicable."""
    first = record["intervals"][0]
    if not (first["minimum_temperature_k"] <= 300.0 and first["maximum_temperature_k"] >= 298.15):
        return None
    enthalpy = _h_over_rt(first["coefficients"], 298.15) * CEA_GAS_CONSTANT_J_MOL_K * 298.15
    return enthalpy - record["heat_of_formation_j_mol"]


def _continuity(record: dict) -> tuple[float, float, list[dict]]:
    """Return the largest relative H/S step, the largest cp step and the known steps."""
    worst = 0.0
    worst_cp = 0.0
    known: list[dict] = []
    for lower, upper in zip(record["intervals"], record["intervals"][1:], strict=False):
        t = lower["maximum_temperature_k"]
        steps = {}
        for label, function in (("h_over_rt", _h_over_rt), ("s_over_r", _s_over_r)):
            a = function(lower["coefficients"], t)
            b = function(upper["coefficients"], t)
            steps[label] = b - a
        cp_low, cp_high = _cp_over_r(lower["coefficients"], t), _cp_over_r(upper["coefficients"], t)
        worst_cp = max(worst_cp, abs(cp_high - cp_low) / abs(cp_low))
        if (record["name"], t) in KNOWN_DISCONTINUITIES:
            known.append({"key": record["name"], "temperature_k": t,
                          "h_over_rt_step": round(steps["h_over_rt"], 6),
                          "s_over_r_step": round(steps["s_over_r"], 6)})
            continue
        worst = max(worst, max(abs(value) for value in steps.values()))
    return worst, worst_cp, known


def _species_record(record: dict) -> dict:
    return {
        "key": record["name"],
        "name": record["name"],
        "thermo_inp_comment": record["comment"],
        "reference_date_code": record["reference_date_code"],
        "formula": record["formula"],
        "molecular_mass_g_mol": record["molecular_mass_g_mol"],
        "heat_of_formation_j_mol_298": record["heat_of_formation_j_mol"],
        "intervals": record["intervals"],
    }


def build_document(text: str, sha256: str) -> dict:
    """Select, verify and assemble the condensed-phase data document."""
    products, reactants = parse_thermo_inp(text)
    existing = _existing_keys()
    gases = [r for name, r in products.items()
             if r["phase"] == 0 and _selected(r) and name.upper() not in existing]
    condensed = [r for r in products.values() if r["phase"] != 0 and _selected(r)]
    if len(gases) != EXPECTED_GAS_COUNT or len(condensed) != EXPECTED_CONDENSED_COUNT:
        raise SystemExit(f"selection changed: {len(gases)} gases, {len(condensed)} condensed")
    for record in gases + condensed:
        if not record["intervals"]:
            raise SystemExit(f"{record['name']} has no polynomial intervals")

    identity: list[dict] = []
    continuity_worst = 0.0
    cp_step_worst = 0.0
    known_steps: list[dict] = []
    for record in gases + condensed:
        step, cp_step, known = _continuity(record)
        continuity_worst = max(continuity_worst, step)
        cp_step_worst = max(cp_step_worst, cp_step)
        known_steps.extend(known)
        residual = _heat_of_formation_residual(record)
        if residual is not None:
            identity.append({"key": record["name"], "residual_j_mol": round(residual, 6)})
            if abs(residual) > HEAT_OF_FORMATION_TOLERANCE_J_MOL:
                raise SystemExit(f"{record['name']}: heat-of-formation identity off by {residual} J/mol")
    if continuity_worst > CONTINUITY_TOLERANCE or len(known_steps) != len(KNOWN_DISCONTINUITIES):
        raise SystemExit(f"interval continuity residual {continuity_worst}")

    groups: dict[str, list[dict]] = {}
    for record in condensed:
        groups.setdefault(_phase_group(record["name"]), []).append(record)
    transitions = []
    for group, records in sorted(groups.items()):
        records.sort(key=lambda item: item["intervals"][0]["minimum_temperature_k"])
        for lower, upper in itertools.pairwise(records):
            t = lower["intervals"][-1]["maximum_temperature_k"]
            if t != upper["intervals"][0]["minimum_temperature_k"]:
                raise SystemExit(f"{group}: phases {lower['name']} and {upper['name']} are not adjacent")
            low_c, high_c = lower["intervals"][-1]["coefficients"], upper["intervals"][0]["coefficients"]
            gibbs = (_h_over_rt(high_c, t) - _s_over_r(high_c, t)) - (_h_over_rt(low_c, t) - _s_over_r(low_c, t))
            enthalpy = (_h_over_rt(high_c, t) - _h_over_rt(low_c, t)) * CEA_GAS_CONSTANT_J_MOL_K * t
            if abs(gibbs) > TRANSITION_GIBBS_TOLERANCE or enthalpy < 0.0:
                raise SystemExit(f"{group}: inconsistent transition at {t} K")
            transitions.append({
                "lower": lower["name"], "upper": upper["name"], "temperature_k": t,
                "gibbs_difference_over_rt": float(f"{gibbs:.6e}"),
                "transition_enthalpy_j_mol": round(enthalpy, 3),
            })

    reactant_records = []
    for name, description in REACTANTS.items():
        record = products.get(name) or reactants.get(name)
        if record is None or not record["intervals"]:
            raise SystemExit(f"thermo.inp lacks a polynomial record for {name}")
        first = record["intervals"][0]
        polynomial = _h_over_rt(first["coefficients"], REACTANT_TEMPERATURE_K) * (
            CEA_GAS_CONSTANT_J_MOL_K * REACTANT_TEMPERATURE_K)
        residual = polynomial - record["heat_of_formation_j_mol"]
        if abs(residual) > HEAT_OF_FORMATION_TOLERANCE_J_MOL:
            raise SystemExit(f"{name}: assigned enthalpy identity off by {residual} J/mol")
        reactant_records.append({
            "key": name,
            "name": description,
            "thermo_inp_comment": record["comment"],
            "reference_date_code": record["reference_date_code"],
            "formula": record["formula"],
            "molecular_mass_g_mol": record["molecular_mass_g_mol"],
            "assigned_enthalpy_j_mol": record["heat_of_formation_j_mol"],
            "assigned_temperature_k": REACTANT_TEMPERATURE_K,
            "polynomial_enthalpy_residual_j_mol": round(residual, 6),
            "phase": "condensed",
        })

    condensed_records = []
    for record in condensed:
        entry = _species_record(record)
        entry["phase_index"] = record["phase"]
        entry["phase_group"] = _phase_group(record["name"])
        condensed_records.append(entry)

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
            "generator": "scripts/sync_nasa_condensed.py",
        },
        "selection": {
            "elements": list(ELEMENTS),
            "gas_rule": ("phase-0 product records of thermo.inp made only of these elements, "
                         "no electron (ion) entry, not already bundled in nasa_glenn_2002.json "
                         "or nasa_cea_thermo.json"),
            "condensed_rule": "phase>0 product records of thermo.inp made only of these elements",
            "reactants": list(REACTANTS),
        },
        "model": {
            "name": "NASA Glenn seven-term cp polynomial with two integration constants",
            "coefficient_order": ["a1", "a2", "a3", "a4", "a5", "a6", "a7", "b1", "b2"],
            "reference_pressure_pa": 100000,
        },
        "verification": {
            "gas_constant_for_identities_j_mol_k": CEA_GAS_CONSTANT_J_MOL_K,
            "heat_of_formation_tolerance_j_mol": HEAT_OF_FORMATION_TOLERANCE_J_MOL,
            "heat_of_formation_identity": identity,
            "maximum_interval_step_h_s": float(f"{continuity_worst:.3e}"),
            "maximum_interval_step_cp_relative": float(f"{cp_step_worst:.3e}"),
            "known_source_discontinuities": known_steps,
            "phase_transitions": transitions,
        },
        "gas_species": [_species_record(record) for record in gases],
        "condensed_species": condensed_records,
        "reactants": reactant_records,
    }


def render(document: dict) -> str:
    """Serialize deterministically: one species record per line, trailing newline."""
    head = {key: value for key, value in document.items()
            if key not in ("gas_species", "condensed_species", "reactants")}
    text = json.dumps(head, indent=1, ensure_ascii=False)[:-2]
    parts = [text]
    for section in ("gas_species", "condensed_species", "reactants"):
        rows = ",\n  ".join(json.dumps(item, ensure_ascii=False) for item in document[section])
        parts.append(f',\n "{section}": [\n  {rows}\n ]')
    return "".join(parts) + "\n}\n"


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
    document = build_document(raw.decode("utf-8"), digest)
    content = render(document)
    if json.loads(content) != document:
        raise SystemExit("render round-trip failed")
    if arguments.check:
        current = TARGET.read_text(encoding="utf-8") if TARGET.is_file() else ""
        if current != content:
            print(f"{TARGET} is out of date")
            return 1
        print("bundled NASA CEA condensed subset matches thermo.inp")
        return 0
    TARGET.write_text(content, encoding="utf-8", newline="\n")
    verification = document["verification"]
    worst = max(abs(item["residual_j_mol"]) for item in verification["heat_of_formation_identity"])
    print(f"wrote {TARGET} ({len(content)} bytes): {len(document['gas_species'])} gases, "
          f"{len(document['condensed_species'])} condensed, {len(document['reactants'])} reactants; "
          f"{len(verification['heat_of_formation_identity'])} heat-of-formation identities, "
          f"worst {worst:.4f} J/mol; {len(verification['phase_transitions'])} transitions")
    return 0


if __name__ == "__main__":
    sys.exit(main())
