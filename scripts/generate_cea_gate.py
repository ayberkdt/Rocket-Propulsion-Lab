"""Record the nine-point NASA CEA reference set for the 1.0 release gate.

The plan's Phase 2 acceptance criterion asks for chamber temperature, ``c*``
and vacuum Isp to agree with NASA CEA at nine or more points across LOX/LH2,
LOX/RP-1 and LOX/CH4 within tolerances declared in advance. This script runs
NASA's own CEA (the ``cea`` Python package, an independent implementation
from NASA Glenn) and writes the reference values together with the exact
inputs, the CEA version and the hash of its thermodynamic database to
``tests/reference_data/cea_gate_points.json``. The comparison itself lives in
``tests/test_cea_gate.py`` and runs without CEA installed.

Usage (in an environment with ``pip install cea``)::

    python scripts/generate_cea_gate.py

References
----------
.. [1] NASA RP-1311, "Computer Program for Calculation of Complex Chemical
       Equilibrium Compositions and Applications, Part I: Analysis".
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
TARGET = ROOT / "tests" / "reference_data" / "cea_gate_points.json"

# Declared before any local result was compared (2026-10-25).
CHAMBER_PRESSURE_BAR = 68.9476  # 1000 psia
AREA_RATIO = 40.0
TOLERANCE = {"chamber_temperature_k": 0.01, "characteristic_velocity_m_s": 0.01,
             "vacuum_specific_impulse_s": 0.01}
REACTANT_TEMPERATURE_K = {"H2(L)": 20.27, "O2(L)": 90.17, "CH4(L)": 111.643, "RP-1": 298.15,
                          "CH6N2(L)": 298.15, "N2O4(L)": 298.15}
POINTS = (
    ("LOX_LH2", "H2(L)", "O2(L)", (4.5, 5.5, 6.5)),
    ("LOX_CH4", "CH4(L)", "O2(L)", (2.8, 3.4, 4.0)),
    ("LOX_RP1", "RP-1", "O2(L)", (2.2, 2.6, 3.0)),
)
G0 = 9.80665
# Domain evidence beyond the nine gate points: every bundled pair at the ends
# of the chamber-pressure range, reference O/F, same area ratio.
DOMAIN_PRESSURES_BAR = (5.0, 250.0)
DOMAIN_PAIRS = (
    ("LOX_LH2", "H2(L)", "O2(L)", 5.55157),
    ("LOX_CH4", "CH4(L)", "O2(L)", 3.4),
    ("LOX_RP1", "RP-1", "O2(L)", 2.7),
    ("NTO_MMH", "CH6N2(L)", "N2O4(L)", 1.65),
)
# Values printed by NASA's own regression output for RP-1311 example 8.
NASA_EXAMPLE8 = {
    "source": "https://github.com/nasa/cea/blob/2f79a647e85737742fb08a127da00d1ea1d069da/"
              "test/main_interface/reference_output/example8.out",
    "conditions": {"fuel": "H2(L)", "oxidizer": "O2(L)", "oxidizer_fuel_ratio": 5.55157,
                   "chamber_pressure_pa": 5331720.0},
    "expected": {
        "chamber_temperature_k": 3383.84462594510,
        "throat_temperature_k": 3185.67317308071,
        "characteristic_velocity_m_s": 2332.33560769645,
        "area_ratio_25_vacuum_effective_exhaust_velocity_m_s": 4348.51382282217,
        "area_ratio_25_effective_exhaust_velocity_m_s": 4124.41484261268,
        "area_ratio_75_vacuum_effective_exhaust_velocity_m_s": 4554.91082019847,
        "area_ratio_75_effective_exhaust_velocity_m_s": 4399.11825930155,
    },
}


def main() -> int:
    import cea
    import numpy as np

    data = Path(cea.__file__).parent / "data" / "thermo.lib"
    thermo_sha256 = hashlib.sha256(data.read_bytes()).hexdigest()
    records = []
    for pair, fuel, oxidizer, ratios in POINTS:
        reactants = cea.Mixture([fuel, oxidizer])
        products = cea.Mixture([fuel, oxidizer], products_from_reactants=True)
        solver = cea.RocketSolver(products, reactants=reactants)
        temperatures = np.array([REACTANT_TEMPERATURE_K[fuel], REACTANT_TEMPERATURE_K[oxidizer]])
        for ratio in ratios:
            solution = cea.RocketSolution(solver)
            weights = reactants.of_ratio_to_weights(np.array([0.0, 1.0]), np.array([1.0, 0.0]), ratio)
            enthalpy = reactants.calc_property(cea.ENTHALPY, weights, temperatures) / cea.R
            solver.solve(solution, weights, CHAMBER_PRESSURE_BAR, supar=[AREA_RATIO], hc=enthalpy, iac=True)
            if not solution.converged:
                raise SystemExit(f"CEA did not converge for {pair} O/F={ratio}")
            records.append({
                "id": f"cea-{pair.lower()}-of-{ratio:g}",
                "pair": pair,
                "fuel": fuel,
                "oxidizer": oxidizer,
                "oxidizer_fuel_ratio": ratio,
                "expected": {
                    "chamber_temperature_k": float(solution.T[0]),
                    "characteristic_velocity_m_s": float(solution.c_star[-1]),
                    "vacuum_specific_impulse_s": float(solution.Isp_vacuum[-1]) / G0,
                    "throat_pressure_pa": float(solution.P[1]) * 1.0e5,
                    "exit_pressure_pa": float(solution.P[-1]) * 1.0e5,
                    "chamber_molecular_mass_g_mol": float(solution.MW[0]),
                },
            })
    domain = []
    for pair, fuel, oxidizer, ratio in DOMAIN_PAIRS:
        reactants = cea.Mixture([fuel, oxidizer])
        products = cea.Mixture([fuel, oxidizer], products_from_reactants=True)
        solver = cea.RocketSolver(products, reactants=reactants)
        temperatures = np.array([REACTANT_TEMPERATURE_K[fuel], REACTANT_TEMPERATURE_K[oxidizer]])
        for pressure in DOMAIN_PRESSURES_BAR:
            solution = cea.RocketSolution(solver)
            weights = reactants.of_ratio_to_weights(np.array([0.0, 1.0]), np.array([1.0, 0.0]), ratio)
            enthalpy = reactants.calc_property(cea.ENTHALPY, weights, temperatures) / cea.R
            solver.solve(solution, weights, pressure, supar=[AREA_RATIO], hc=enthalpy, iac=True)
            if not solution.converged:
                raise SystemExit(f"CEA did not converge for {pair} p={pressure} bar")
            domain.append({
                "pair": pair, "fuel": fuel, "oxidizer": oxidizer, "oxidizer_fuel_ratio": ratio,
                "chamber_pressure_pa": pressure * 1.0e5,
                "expected": {
                    "chamber_temperature_k": float(solution.T[0]),
                    "characteristic_velocity_m_s": float(solution.c_star[-1]),
                    "vacuum_specific_impulse_s": float(solution.Isp_vacuum[-1]) / G0,
                },
            })
    document = {
        "schema_version": 1,
        "generated": datetime.date.today().isoformat(),
        "generator": "scripts/generate_cea_gate.py",
        "reference": {
            "program": "NASA CEA (cea Python package)",
            "version": cea.__version__,
            "url": "https://pypi.org/project/cea/",
            "repository": "https://github.com/nasa/cea",
            "thermo_lib_sha256": thermo_sha256,
        },
        "conditions": {
            "chamber_pressure_pa": CHAMBER_PRESSURE_BAR * 1.0e5,
            "area_ratio": AREA_RATIO,
            "combustor": "infinite area (IAC)",
            "expansion": "shifting equilibrium",
            "reactant_temperature_k": REACTANT_TEMPERATURE_K,
            "standard_gravity_m_s2": G0,
        },
        "declared_relative_tolerance": TOLERANCE,
        "declared_on": "2026-10-25",
        "points": records,
        "nasa_reference_output_example8": NASA_EXAMPLE8,
        "pressure_domain_points": domain,
    }
    TARGET.write_text(json.dumps(document, indent=1) + "\n", encoding="utf-8", newline="\n")
    print(f"wrote {len(records)} CEA points to {TARGET}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
