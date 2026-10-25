"""Record NASA CEA reference points for equilibrium with condensed products.

The points exercise the condensed-phase provider of
``rocket_propulsion.thermochemistry.condensed_equilibrium``:

* fuel-rich LOX/CH4 and LOX/RP-1 where graphite, C(gr), forms (rocket, IAC);
* aluminium with ammonium perchlorate (20/80 by mass, CEA reactants
  ``AL(cr)`` and ``NH4CLO4(I)``, no binder) producing liquid alumina, with
  area ratios that cross the 2327 K melting point of AL2O3 (rocket, IAC), and
  aluminium with oxygen (TP below and above the melting point, HP);
* hydrogen with oxygen (3:1 by moles) at 250, 300 and 350 K producing ice
  or liquid water;
* fuel-rich CH4/O2 at 1000 K producing C(gr) (TP).

NASA CEA (the ``cea`` Python package, an implementation independent of this
repository) is run with its full product database for the reactant elements
(``products_from_reactants=True``); for the elements Al, Cl, C, H, O and N that
set equals the provider's bundled set. The script writes the inputs, CEA's
version, the SHA-256 of CEA's ``thermo.lib``, of the pinned ``thermo.inp`` the
bundled data came from and of ``nasa_cea_condensed.json``, and the results to
``tests/reference_data/cea_condensed_points.json``. The comparison lives in
``tests/test_condensed_equilibrium.py`` and runs without CEA installed.

Usage (in an environment with ``pip install cea``)::

    python scripts/generate_cea_condensed_points.py

References
----------
.. [1] NASA-RP1311-P1-CONDENSED. Gordon and McBride, NASA RP-1311 Part I (1994).
       https://ntrs.nasa.gov/citations/19950013764
.. [2] NASA-CEA-PY. NASA Glenn Research Center, "CEA: Chemical Equilibrium with
       Applications", Python package 3.3.4. https://pypi.org/project/cea/
"""

from __future__ import annotations

import datetime
import hashlib
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
TARGET = ROOT / "tests" / "reference_data" / "cea_condensed_points.json"
DATA = ROOT / "src" / "rocket_propulsion" / "thermochemistry" / "data" / "nasa_cea_condensed.json"
THERMO_INP_SHA256 = "fa7746572952d74e249e818a82a35c113829742fb421a308e167185528884363"

CHAMBER_PRESSURE_BAR = 68.9476  # 1000 psia, as the release gate
REACTANT_TEMPERATURE_K = {"H2(L)": 20.27, "O2(L)": 90.17, "CH4(L)": 111.643, "RP-1": 298.15,
                          "AL(cr)": 298.15, "NH4CLO4(I)": 298.15, "O2": 298.15, "H2": 298.15}
G0 = 9.80665
GAS_MOLE_FRACTION_FLOOR = 1.0e-4

ROCKETS = (
    ("lox-ch4-of-0.8", ["CH4(L)", "O2(L)"], 0.8, None, [10.0, 40.0]),
    ("lox-ch4-of-1.2", ["CH4(L)", "O2(L)"], 1.2, None, [40.0]),
    ("lox-rp1-of-1.0", ["RP-1", "O2(L)"], 1.0, None, [40.0]),
    ("al-ap-20-80", ["AL(cr)", "NH4CLO4(I)"], None, [0.2, 0.8], [10.0, 40.0, 60.0, 80.0, 100.0]),
)
# (id, reactants, mass fractions or O/F, kind, temperature K or None, pressure bar)
EQUILIBRIA = (
    ("h2-o2-tp-250", ["H2", "O2"], "moles:3:1", "TP", 250.0, 1.0),
    ("h2-o2-tp-300", ["H2", "O2"], "moles:3:1", "TP", 300.0, 1.0),
    ("h2-o2-tp-350", ["H2", "O2"], "moles:3:1", "TP", 350.0, 1.0),
    ("al-o2-tp-2000", ["AL(cr)", "O2"], [0.5, 0.5], "TP", 2000.0, 1.0),
    ("al-o2-tp-3000", ["AL(cr)", "O2"], [0.5, 0.5], "TP", 3000.0, 1.0),
    ("al-o2-hp-20bar", ["AL(cr)", "O2"], [0.5, 0.5], "HP", None, 20.0),
    ("ch4-o2-of-1.0-tp-1000", ["CH4(L)", "O2(L)"], "of:1.0", "TP", 1000.0, 10.0),
)


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _condensed_keys() -> set[str]:
    document = json.loads(DATA.read_text(encoding="utf-8"))
    return {record["key"] for record in document["condensed_species"]}


def _composition(mass_fractions: dict, mole_fractions: dict, index: int | None, condensed: set[str]):
    def value(table, key):
        item = table[key]
        return float(item if index is None else item[index])

    solid = {key: value(mass_fractions, key) for key in mass_fractions
             if key in condensed and value(mass_fractions, key) > 0.0}
    gases = {key: value(mole_fractions, key) for key in mole_fractions
             if key not in condensed and value(mole_fractions, key) >= GAS_MOLE_FRACTION_FLOOR}
    return solid, gases


def main() -> int:
    import cea
    import numpy as np

    condensed = _condensed_keys()
    thermo_lib = Path(cea.__file__).parent / "data" / "thermo.lib"
    rockets = []
    for identifier, names, ratio, weights, area_ratios in ROCKETS:
        reactants = cea.Mixture(names)
        products = cea.Mixture(names, products_from_reactants=True)
        solver = cea.RocketSolver(products, reactants=reactants)
        solution = cea.RocketSolution(solver)
        if weights is None:
            mass = reactants.of_ratio_to_weights(np.array([0.0, 1.0]), np.array([1.0, 0.0]), ratio)
        else:
            mass = np.array(weights)
        temperatures = np.array([REACTANT_TEMPERATURE_K[name] for name in names])
        enthalpy = reactants.calc_property(cea.ENTHALPY, mass, temperatures) / cea.R
        solver.solve(solution, mass, CHAMBER_PRESSURE_BAR, supar=area_ratios, hc=enthalpy, iac=True)
        if not solution.converged:
            raise SystemExit(f"CEA did not converge for {identifier}")
        stations = []
        labels = ["chamber", "throat"] + ["exit"] * len(area_ratios)
        for index, label in enumerate(labels):
            solid, gases = _composition(solution.mass_fractions, solution.mole_fractions, index, condensed)
            station = {
                "name": label,
                "area_ratio": None if label == "chamber" else float(solution.ae_at[index]),
                "temperature_k": float(solution.T[index]),
                "pressure_pa": float(solution.P[index]) * 1.0e5,
                "molecular_mass_g_mol": float(solution.M[index]),
                "mixture_molecular_mass_g_mol": float(solution.MW[index]),
                "condensed_mass_fractions": solid,
                "condensed_mass_fraction": sum(solid.values()),
                "gas_mole_fractions": gases,
            }
            if label == "exit":
                station["vacuum_specific_impulse_s"] = float(solution.Isp_vacuum[index]) / G0
            stations.append(station)
        rockets.append({
            "id": identifier,
            "reactants": names,
            "mass_fractions": [float(value) for value in mass],
            "oxidizer_fuel_ratio": ratio,
            "reactant_temperature_k": [REACTANT_TEMPERATURE_K[name] for name in names],
            "chamber_pressure_pa": CHAMBER_PRESSURE_BAR * 1.0e5,
            "characteristic_velocity_m_s": float(solution.c_star[-1]),
            "stations": stations,
        })
    equilibria = []
    for identifier, names, amount, kind, temperature, pressure in EQUILIBRIA:
        reactants = cea.Mixture(names)
        products = cea.Mixture(names, products_from_reactants=True)
        solver = cea.EqSolver(products, reactants=reactants)
        solution = cea.EqSolution(solver)
        if isinstance(amount, str) and amount.startswith("moles:"):
            mass = reactants.moles_to_weights(np.array([float(x) for x in amount.split(":")[1:]]))
        elif isinstance(amount, str) and amount.startswith("of:"):
            mass = reactants.of_ratio_to_weights(np.array([0.0, 1.0]), np.array([1.0, 0.0]),
                                                 float(amount.split(":")[1]))
        else:
            mass = np.array(amount)
        if kind == "TP":
            solver.solve(solution, cea.TP, temperature, pressure, mass)
        else:
            temperatures = np.array([REACTANT_TEMPERATURE_K[name] for name in names])
            enthalpy = reactants.calc_property(cea.ENTHALPY, mass, temperatures) / cea.R
            solver.solve(solution, cea.HP, enthalpy, pressure, mass)
        if not solution.converged:
            raise SystemExit(f"CEA did not converge for {identifier}")
        solid, gases = _composition(solution.mass_fractions, solution.mole_fractions, None, condensed)
        equilibria.append({
            "id": identifier,
            "problem_type": kind,
            "reactants": names,
            "mass_fractions": [float(value) for value in mass],
            "reactant_temperature_k": [REACTANT_TEMPERATURE_K[name] for name in names],
            "pressure_pa": pressure * 1.0e5,
            "temperature_k": float(solution.T),
            "molecular_mass_g_mol": float(solution.M),
            "mixture_molecular_mass_g_mol": float(solution.MW),
            "enthalpy_j_kg": float(solution.enthalpy) * 1.0e3,
            "entropy_j_kg_k": float(solution.entropy) * 1.0e3,
            "condensed_mass_fractions": solid,
            "condensed_mass_fraction": sum(solid.values()),
            "gas_mole_fractions": gases,
        })
    document = {
        "schema_version": 1,
        "generated": datetime.date.today().isoformat(),
        "generator": "scripts/generate_cea_condensed_points.py",
        "reference": {
            "program": "NASA CEA (cea Python package)",
            "version": cea.__version__,
            "url": "https://pypi.org/project/cea/",
            "repository": "https://github.com/nasa/cea",
            "thermo_lib_sha256": _sha256(thermo_lib),
            "product_set": "products_from_reactants=True (full CEA database for the reactant elements)",
            "units_note": "CEA enthalpy and entropy are converted from kJ/kg and kJ/(kg K)",
        },
        "local_data": {
            "thermo_inp_sha256": THERMO_INP_SHA256,
            "nasa_cea_condensed_json_sha256": _sha256(DATA),
        },
        "conditions": {
            "combustor": "infinite area (IAC)",
            "expansion": "shifting equilibrium, zero velocity and thermal lag",
            "standard_gravity_m_s2": G0,
            "gas_mole_fraction_floor": GAS_MOLE_FRACTION_FLOOR,
        },
        "rocket_points": rockets,
        "equilibrium_points": equilibria,
    }
    TARGET.write_text(json.dumps(document, indent=1) + "\n", encoding="utf-8", newline="\n")
    print(f"wrote {len(rockets)} rocket and {len(equilibria)} equilibrium CEA points to {TARGET}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
