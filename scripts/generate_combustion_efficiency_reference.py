"""Record NASA CEA composition, stoichiometry and O/F-optimum references.

Runs NASA's own CEA (``cea`` Python package 3.3.4) and writes, together with
the CEA version and the hash of its ``thermo.lib``, to
``tests/reference_data/combustion_efficiency.json``:

* equilibrium mole and mass fractions, including trace radicals, at the
  chamber, throat and an area-ratio-40 exit of infinite-area-combustor rocket
  problems for LOX/LH2 (three O/F, three pressures), LOX/CH4 and LOX/RP-1;
* stoichiometric O/F (chemical equivalence ratio ``r = 1``, RP-1311 Part I
  eq. (9.18)) and O/F at assigned ``r`` and ``phi`` (eq. (9.19)) from CEA's
  ``Mixture.chem_eq_ratio_to_of_ratio`` and ``weight_eq_ratio_to_of_ratio``;
* the O/F that maximises CEA vacuum specific impulse, and the O/F that
  maximises ``rho_0 I_vac`` with declared test densities, found by golden
  section on CEA itself.

It also carries published values read from NASA documents (with page
locators) that the tests use as independent anchors. The comparison lives in
``tests/test_combustion_efficiency.py`` and runs without CEA installed.

Usage (in an environment with ``pip install cea``)::

    python scripts/generate_combustion_efficiency_reference.py

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
import math
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
TARGET = ROOT / "tests" / "reference_data" / "combustion_efficiency.json"

AREA_RATIO = 40.0
REACTANT_TEMPERATURE_K = {"H2(L)": 20.27, "O2(L)": 90.17, "CH4(L)": 111.643, "RP-1": 298.15}
PAIRS = {"LOX_LH2": ("H2(L)", "O2(L)"), "LOX_CH4": ("CH4(L)", "O2(L)"), "LOX_RP1": ("RP-1", "O2(L)")}
# (id, pair, O/F, chamber pressure in bar)
COMPOSITION_CASES = (
    ("lox_lh2-of4-p68.9", "LOX_LH2", 4.0, 68.9476),
    ("lox_lh2-of6-p68.9", "LOX_LH2", 6.0, 68.9476),
    ("lox_lh2-of8-p68.9", "LOX_LH2", 8.0, 68.9476),
    ("lox_lh2-of6-p10", "LOX_LH2", 6.0, 10.0),
    ("lox_lh2-of6-p200", "LOX_LH2", 6.0, 200.0),
    ("lox_ch4-of3.4-p68.9", "LOX_CH4", 3.4, 68.9476),
    ("lox_rp1-of2.7-p68.9", "LOX_RP1", 2.7, 68.9476),
)
# Declared before any local result was compared (2026-10-25). Relative
# tolerance by mole-fraction band; species below the last band are not compared.
COMPOSITION_TOLERANCE = {
    "bands": [
        {"minimum_mole_fraction": 1.0e-2, "relative": 2.0e-3},
        {"minimum_mole_fraction": 1.0e-6, "relative": 1.0e-2},
        {"minimum_mole_fraction": 1.0e-10, "relative": 5.0e-2},
    ],
    "temperature_relative": 1.0e-3,
    "molecular_mass_relative": 1.0e-3,
}
OPTIMUM_BRACKETS = {"LOX_LH2": (3.5, 7.0), "LOX_CH4": (2.5, 4.0), "LOX_RP1": (2.0, 3.3)}
DENSITY_BRACKETS = {"LOX_LH2": (5.0, 11.0), "LOX_CH4": (3.0, 5.0), "LOX_RP1": (2.0, 3.5)}
OPTIMUM_PRESSURE_BAR = 68.9476
# Declared test inputs for the density-impulse objective, kg/m^3. They exercise
# the optimiser and are NOT property data for design use.
TEST_DENSITIES_KG_M3 = {"H2(L)": 70.0, "O2(L)": 1140.0, "CH4(L)": 420.0, "RP-1": 800.0}
OPTIMUM_TOLERANCE = {"oxidizer_fuel_ratio_absolute": 0.03, "objective_relative": 2.0e-4}
G0 = 9.80665

PUBLISHED = {
    "sp8089_table_iv_mixing_factor": {
        "source": "NASA SP-8089, Liquid Rocket Engine Injectors (1976), Table IV, printed p. 20 "
                  "(PDF p. 29); equation (1), printed p. 19 (PDF p. 28)",
        "url": "https://ntrs.nasa.gov/citations/19760023196",
        "values": {"1-on-1": 1.0, "2-on-2": 1.0, "2-on-1": 1.6, "3-on-1": 3.5, "4-on-1": 9.4,
                   "5-on-1": 27.5},
    },
    "sp125_table_4_1_characteristic_length_in": {
        "source": "NASA SP-125 (Huzel and Huang, 1967), Table 4-1, printed p. 87 (PDF p. 96)",
        "url": "https://ntrs.nasa.gov/citations/19710019929",
        "values": {
            "LOX_LH2_GH2_injection": [22.0, 28.0],
            "LOX_LH2_LH2_injection": [30.0, 40.0],
            "LOX_RP1": [40.0, 50.0],
            "NTO_HYDRAZINE_BASE": [30.0, 35.0],
        },
        "stay_time_range_s": [0.002, 0.040],
        "stay_time_locator": "printed p. 87, first column: L* 15-120 in for stay times 0.002-0.040 s",
    },
    "sp125_table_1_7_density_impulse_rows": {
        "source": "NASA SP-125, Table 1-7, printed p. 25 (PDF p. 34); columns Is (lb-sec/lb), "
                  "d (g/cc), Isd (sec g/cc); definition Isd = Is d, equation (1-45), printed p. 20",
        "url": "https://ntrs.nasa.gov/citations/19710019929",
        "rows": [
            {"oxidizer": "IRFNA", "fuel": "UDMH", "specific_impulse_s": 276.0, "bulk_density_g_cc": 1.26,
             "density_impulse_s_g_cc": 348.0},
            {"oxidizer": "IRFNA", "fuel": "RP-1", "specific_impulse_s": 268.0, "bulk_density_g_cc": 1.35,
             "density_impulse_s_g_cc": 362.0},
            {"oxidizer": "95% H2O2", "fuel": "RP-1", "specific_impulse_s": 273.0, "bulk_density_g_cc": 1.30,
             "density_impulse_s_g_cc": 355.0},
        ],
    },
}


def _solver(cea, fuel, oxidizer):
    reactants = cea.Mixture([fuel, oxidizer])
    products = cea.Mixture([fuel, oxidizer], products_from_reactants=True)
    return reactants, cea.RocketSolver(products, reactants=reactants)


def _solve(cea, np, reactants, solver, fuel, oxidizer, ratio, pressure_bar):
    solution = cea.RocketSolution(solver)
    temperatures = np.array([REACTANT_TEMPERATURE_K[fuel], REACTANT_TEMPERATURE_K[oxidizer]])
    weights = reactants.of_ratio_to_weights(np.array([0.0, 1.0]), np.array([1.0, 0.0]), ratio)
    enthalpy = reactants.calc_property(cea.ENTHALPY, weights, temperatures) / cea.R
    solver.solve(solution, weights, pressure_bar, supar=[AREA_RATIO], hc=enthalpy, iac=True)
    if not solution.converged:
        raise SystemExit(f"CEA did not converge: {fuel}/{oxidizer} O/F={ratio} p={pressure_bar}")
    return solution


def _golden(function, low, high, tolerance):
    golden = (math.sqrt(5.0) - 1.0) / 2.0
    first, second = high - golden * (high - low), low + golden * (high - low)
    f_first, f_second = function(first), function(second)
    while high - low > tolerance:
        if f_first > f_second:
            high, second, f_second = second, first, f_first
            first = high - golden * (high - low)
            f_first = function(first)
        else:
            low, first, f_first = first, second, f_second
            second = low + golden * (high - low)
            f_second = function(second)
    best = 0.5 * (low + high)
    return best, function(best)


def main() -> int:
    import cea
    import numpy as np

    thermo_sha256 = hashlib.sha256((Path(cea.__file__).parent / "data" / "thermo.lib").read_bytes()).hexdigest()
    compositions = []
    for case_id, pair, ratio, pressure in COMPOSITION_CASES:
        fuel, oxidizer = PAIRS[pair]
        reactants, solver = _solver(cea, fuel, oxidizer)
        solution = _solve(cea, np, reactants, solver, fuel, oxidizer, ratio, pressure)
        stations = {}
        for index, name in enumerate(("chamber", "throat", "exit")):
            stations[name] = {
                "pressure_pa": float(solution.P[index]) * 1.0e5,
                "temperature_k": float(solution.T[index]),
                "molecular_mass_g_mol": float(solution.MW[index]),
                "mole_fractions": {key: float(values[index]) for key, values in solution.mole_fractions.items()
                                   if float(values[index]) >= 1.0e-12},
                "mass_fractions": {key: float(values[index]) for key, values in solution.mass_fractions.items()
                                   if float(values[index]) >= 1.0e-12},
            }
        compositions.append({"id": case_id, "pair": pair, "fuel": fuel, "oxidizer": oxidizer,
                             "oxidizer_fuel_ratio": ratio, "chamber_pressure_pa": pressure * 1.0e5,
                             "area_ratio": AREA_RATIO, "stations": stations})

    stoichiometry = []
    for pair, (fuel, oxidizer) in PAIRS.items():
        reactants = cea.Mixture([fuel, oxidizer])
        oxidant_weights, fuel_weights = np.array([0.0, 1.0]), np.array([1.0, 0.0])
        stoichiometry.append({
            "pair": pair, "fuel": fuel, "oxidizer": oxidizer,
            "stoichiometric_oxidizer_fuel_ratio": float(
                reactants.chem_eq_ratio_to_of_ratio(oxidant_weights, fuel_weights, 1.0)),
            "oxidizer_fuel_ratio_at_r_0.9": float(
                reactants.chem_eq_ratio_to_of_ratio(oxidant_weights, fuel_weights, 0.9)),
            "oxidizer_fuel_ratio_at_phi_0.8": float(
                reactants.weight_eq_ratio_to_of_ratio(oxidant_weights, fuel_weights, 0.8)),
            "oxidizer_fuel_ratio_at_phi_1.25": float(
                reactants.weight_eq_ratio_to_of_ratio(oxidant_weights, fuel_weights, 1.25)),
        })

    optima = []
    for pair, (fuel, oxidizer) in PAIRS.items():
        reactants, solver = _solver(cea, fuel, oxidizer)
        low, high = OPTIMUM_BRACKETS[pair]

        def isp(ratio, reactants=reactants, solver=solver, fuel=fuel, oxidizer=oxidizer):
            solution = _solve(cea, np, reactants, solver, fuel, oxidizer, ratio, OPTIMUM_PRESSURE_BAR)
            return float(solution.Isp_vacuum[-1]) / G0

        def density_isp(ratio, isp=isp, fuel=fuel, oxidizer=oxidizer):
            rho_o, rho_f = TEST_DENSITIES_KG_M3[oxidizer], TEST_DENSITIES_KG_M3[fuel]
            bulk = (ratio + 1.0) * rho_o * rho_f / (rho_o + ratio * rho_f)  # RP-1311 eq. (9.13)
            return bulk * isp(ratio)

        best_isp = _golden(isp, low, high, 1.0e-4)
        best_density = _golden(density_isp, *DENSITY_BRACKETS[pair], 1.0e-4)
        optima.append({
            "pair": pair, "fuel": fuel, "oxidizer": oxidizer, "bracket": [low, high],
            "chamber_pressure_pa": OPTIMUM_PRESSURE_BAR * 1.0e5, "area_ratio": AREA_RATIO,
            "maximum_vacuum_specific_impulse": {"oxidizer_fuel_ratio": best_isp[0], "objective": best_isp[1]},
            "density_bracket": list(DENSITY_BRACKETS[pair]),
            "maximum_density_vacuum_impulse": {"oxidizer_fuel_ratio": best_density[0],
                                               "objective": best_density[1],
                                               "objective_unit": "kg s / m^3"},
        })

    document = {
        "schema_version": 1,
        "generated": datetime.date.today().isoformat(),
        "generator": "scripts/generate_combustion_efficiency_reference.py",
        "reference": {
            "program": "NASA CEA (cea Python package)",
            "version": cea.__version__,
            "url": "https://pypi.org/project/cea/",
            "repository": "https://github.com/nasa/cea",
            "thermo_lib_sha256": thermo_sha256,
        },
        "conditions": {"combustor": "infinite area (IAC)", "expansion": "shifting equilibrium",
                       "reactant_temperature_k": REACTANT_TEMPERATURE_K,
                       "standard_gravity_m_s2": G0, "test_densities_kg_m3": TEST_DENSITIES_KG_M3},
        "declared_on": "2026-10-25",
        "composition_tolerance": COMPOSITION_TOLERANCE,
        "optimum_tolerance": OPTIMUM_TOLERANCE,
        "compositions": compositions,
        "stoichiometry": stoichiometry,
        "optima": optima,
        "published": PUBLISHED,
    }
    TARGET.write_text(json.dumps(document, indent=1) + "\n", encoding="utf-8", newline="\n")
    print(f"wrote {len(compositions)} compositions, {len(stoichiometry)} stoichiometries and "
          f"{len(optima)} optima to {TARGET}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
