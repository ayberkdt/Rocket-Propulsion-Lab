"""Record NASA CEA finite-area-combustor (FAC) reference points.

Runs NASA's own CEA (the ``cea`` Python package 3.3.4) in its FAC rocket mode
(``RocketSolver.solve(..., iac=False, ac_at=...)`` or ``mdot=...``) for
LOX/LH2, LOX/CH4 and LOX/RP-1 and writes the inputs, the CEA version, the hash
of its ``thermo.lib`` and the station values to
``tests/reference_data/cea_fac_points.json``. CEA's FAC stations are, in order:
injector, fictitious infinite-area point, combustor end, throat, exits
(RP-1311 Part I, section 6 [1]_). The comparison lives in
``tests/test_combustion_efficiency.py`` and runs without CEA installed.

The file also carries the values printed by NASA's CEA2 regression outputs of
RP-1311 examples 9 (``acat=1.58``) and 10 (``mdot/Ac = 1333.9``) [2]_ as
published, independent anchors.

Usage (in an environment with ``pip install cea``)::

    python scripts/generate_cea_fac_points.py

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
TARGET = ROOT / "tests" / "reference_data" / "cea_fac_points.json"

# Declared before any local result was compared (2026-10-25).
TOLERANCE = {
    "stagnation_pressure_ratio": 2.0e-4,
    "combustor_end_pressure_pa": 5.0e-4,
    "combustor_end_temperature_k": 1.0e-3,
    "combustor_end_mach_number": 2.0e-3,
    "throat_pressure_pa": 5.0e-4,
    "characteristic_velocity_m_s": 1.0e-3,
    "vacuum_effective_exhaust_velocity_m_s": 1.0e-3,
    "contraction_ratio": 2.0e-4,
}
REACTANT_TEMPERATURE_K = {"H2(L)": 20.27, "O2(L)": 90.17, "CH4(L)": 111.643, "RP-1": 298.15}
# (id, pair, fuel, oxidizer, O/F, injector pressure bar, option, value, supersonic area ratios)
POINTS = (
    ("rp1311-ex9", "LOX_LH2", "H2(L)", "O2(L)", 5.55157, 53.3172, "ac_at", 1.58, (25.0, 50.0, 75.0)),
    ("rp1311-ex10", "LOX_LH2", "H2(L)", "O2(L)", 5.55157, 53.3172, "mdot", 1333.9, (25.0, 50.0, 75.0)),
    ("lox_lh2-cr1.3", "LOX_LH2", "H2(L)", "O2(L)", 6.0, 68.9476, "ac_at", 1.3, (40.0,)),
    ("lox_lh2-cr3", "LOX_LH2", "H2(L)", "O2(L)", 6.0, 68.9476, "ac_at", 3.0, (40.0,)),
    ("lox_ch4-cr2", "LOX_CH4", "CH4(L)", "O2(L)", 3.4, 68.9476, "ac_at", 2.0, (40.0,)),
    ("lox_rp1-cr2", "LOX_RP1", "RP-1", "O2(L)", 2.7, 68.9476, "ac_at", 2.0, (40.0,)),
)
CEA_REPOSITORY_COMMIT = "2f79a647e85737742fb08a127da00d1ea1d069da"
NASA_OUTPUTS = {
    "example9": {
        "source": f"https://github.com/nasa/cea/blob/{CEA_REPOSITORY_COMMIT}/"
                  "test/main_interface/reference_output/example9.out",
        "sha256": "0fabc39640cadb5edead23275dcd269ea4f8fa9cc9971521f50a9ef784157c27",
        "locator": "lines 126-160 (first table)",
        "conditions": {"oxidizer_fuel_ratio": 5.55157, "injector_pressure_pa": 5331720.0,
                       "contraction_ratio": 1.58},
        "expected": {
            "injector_to_infinite_pressure_ratio": 1.084767,
            "combustor_end_pressure_pa": 0.446032972915175e7,
            "combustor_end_temperature_k": 0.334106999141470e4,
            "combustor_end_mach_number": 0.413196100264602,
            "throat_pressure_pa": 0.282669122299647e7,
            "characteristic_velocity_m_s": 0.233099548251667e4,
            "area_ratio_25_vacuum_effective_exhaust_velocity_m_s": 0.434767250611254e4,
            "area_ratio_50_vacuum_effective_exhaust_velocity_m_s": 0.448663302519288e4,
        },
    },
    "example10": {
        "source": f"https://github.com/nasa/cea/blob/{CEA_REPOSITORY_COMMIT}/"
                  "test/main_interface/reference_output/example10.out",
        "sha256": "9d3a2274529f01cee33cf2f88419173efb66496a3f265296b653a85d04fe2663",
        "locator": "lines 160-194 (first table)",
        "conditions": {"oxidizer_fuel_ratio": 5.55157, "injector_pressure_pa": 5331720.0,
                       "mass_flux_per_chamber_area_kg_m2_s": 1333.9},
        "expected": {
            "injector_to_infinite_pressure_ratio": 1.084668,
            "contraction_ratio": 0.158090453944294e1,
            "combustor_end_pressure_pa": 0.446133802307749e7,
            "combustor_end_mach_number": 0.412907602531960,
            "characteristic_velocity_m_s": 0.233099699394121e4,
        },
    },
}


def main() -> int:
    import cea
    import numpy as np

    data = Path(cea.__file__).parent / "data" / "thermo.lib"
    thermo_sha256 = hashlib.sha256(data.read_bytes()).hexdigest()
    records = []
    for point_id, pair, fuel, oxidizer, ratio, pressure_bar, option, value, supar in POINTS:
        reactants = cea.Mixture([fuel, oxidizer])
        products = cea.Mixture([fuel, oxidizer], products_from_reactants=True)
        solver = cea.RocketSolver(products, reactants=reactants)
        solution = cea.RocketSolution(solver)
        temperatures = np.array([REACTANT_TEMPERATURE_K[fuel], REACTANT_TEMPERATURE_K[oxidizer]])
        weights = reactants.of_ratio_to_weights(np.array([0.0, 1.0]), np.array([1.0, 0.0]), ratio)
        enthalpy = reactants.calc_property(cea.ENTHALPY, weights, temperatures) / cea.R
        solver.solve(solution, weights, pressure_bar, supar=list(supar), hc=enthalpy, iac=False,
                     **{option: value})
        if not solution.converged:
            raise SystemExit(f"CEA FAC did not converge for {point_id}")
        exits = []
        for offset, area_ratio in enumerate(supar):
            index = 4 + offset
            exits.append({
                "area_ratio": area_ratio,
                "pressure_pa": float(solution.P[index]) * 1.0e5,
                "temperature_k": float(solution.T[index]),
                "effective_exhaust_velocity_m_s": float(solution.Isp[index]),
                "vacuum_effective_exhaust_velocity_m_s": float(solution.Isp_vacuum[index]),
            })
        records.append({
            "id": point_id,
            "pair": pair,
            "fuel": fuel,
            "oxidizer": oxidizer,
            "oxidizer_fuel_ratio": ratio,
            "injector_pressure_pa": pressure_bar * 1.0e5,
            "option": option,
            "option_value": value,
            "expected": {
                "injector_temperature_k": float(solution.T[0]),
                "infinite_area_pressure_pa": float(solution.P[1]) * 1.0e5,
                "stagnation_pressure_ratio": float(solution.P[1]) / pressure_bar,
                "infinite_area_entropy_j_kg_k": float(solution.entropy[1]) * 1000.0,
                "injector_entropy_j_kg_k": float(solution.entropy[0]) * 1000.0,
                "contraction_ratio": float(solution.ae_at[2]),
                "combustor_end_pressure_pa": float(solution.P[2]) * 1.0e5,
                "combustor_end_temperature_k": float(solution.T[2]),
                "combustor_end_mach_number": float(solution.Mach[2]),
                "throat_pressure_pa": float(solution.P[3]) * 1.0e5,
                "throat_temperature_k": float(solution.T[3]),
                "characteristic_velocity_m_s": float(solution.c_star[3]),
                "exits": exits,
            },
        })
    document = {
        "schema_version": 1,
        "generated": datetime.date.today().isoformat(),
        "generator": "scripts/generate_cea_fac_points.py",
        "reference": {
            "program": "NASA CEA (cea Python package)",
            "version": cea.__version__,
            "url": "https://pypi.org/project/cea/",
            "repository": "https://github.com/nasa/cea",
            "thermo_lib_sha256": thermo_sha256,
            "station_order": ["injector", "infinite_area", "combustor_end", "throat", "exits..."],
        },
        "conditions": {
            "combustor": "finite area (FAC), RP-1311 Part I section 6.4",
            "expansion": "shifting equilibrium",
            "reactant_temperature_k": REACTANT_TEMPERATURE_K,
        },
        "declared_relative_tolerance": TOLERANCE,
        "declared_on": "2026-10-25",
        "points": records,
        "nasa_reference_outputs": NASA_OUTPUTS,
    }
    TARGET.write_text(json.dumps(document, indent=1) + "\n", encoding="utf-8", newline="\n")
    print(f"wrote {len(records)} CEA FAC points to {TARGET}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
