"""Record NASA CEA transport properties and equilibrium derivatives as test references.

This script runs NASA's own CEA (the ``cea`` Python package, an independent
implementation from NASA Glenn) with ``transport=True`` and writes, for every
state, the exact inputs, the full product mole fractions and CEA's viscosity,
frozen and equilibrium thermal conductivity, Prandtl numbers, frozen and
equilibrium specific heats, ``gamma_s`` and sound speed to
``tests/reference_data/cea_transport_points.json``, together with the CEA
version and the SHA-256 of the ``thermo.lib`` and ``trans.lib`` files it used.
The comparison lives in ``tests/test_transport_and_derivatives.py`` and runs
without CEA installed.

Three groups of states are recorded:

* ``pure_species``: one gas at a time (CEA ``only`` product list) to check the
  species fits of RP-1311 Part I equation (5.1) [1]_;
* ``mixtures``: small fixed product sets solved at assigned T and p, to check
  the mixture rules (5.3)--(5.7) and the reaction terms (5.8)--(5.18);
* ``rocket``: chamber and throat of the nine LOX/LH2, LOX/CH4 and LOX/RP-1
  points of ``tests/reference_data/cea_gate_points.json`` (68.9476 bar,
  infinite-area combustor, area ratio 40).

Usage (in an environment with ``pip install cea``)::

    python scripts/generate_cea_transport_points.py

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
TARGET = ROOT / "tests" / "reference_data" / "cea_transport_points.json"

CHAMBER_PRESSURE_BAR = 68.9476  # 1000 psia, as in scripts/generate_cea_gate.py
AREA_RATIO = 40.0
REACTANT_TEMPERATURE_K = {"H2(L)": 20.27, "O2(L)": 90.17, "CH4(L)": 111.643, "RP-1": 298.15}
ROCKET_POINTS = (
    ("LOX_LH2", "H2(L)", "O2(L)", (4.5, 5.5, 6.5)),
    ("LOX_CH4", "CH4(L)", "O2(L)", (2.8, 3.4, 4.0)),
    ("LOX_RP1", "RP-1", "O2(L)", (2.2, 2.6, 3.0)),
)
PURE_SPECIES = (
    ("N2", (300.0, 1000.0, 3000.0)),
    ("O2", (1000.0, 3000.0)),
    ("H2", (300.0, 2000.0)),
    ("H2O", (500.0, 1500.0, 3000.0)),
    ("CO", (2000.0,)),
    ("CO2", (1000.0, 3000.0)),
    ("OH", (3000.0,)),
    ("H", (3000.0,)),
    ("O", (3000.0,)),
    ("CH4", (800.0,)),
)
# (id, product species, reactant moles, temperature K, pressure bar)
MIXTURES = (
    ("air-like-n2-o2-1000k", ("N2", "O2"), (0.79, 0.21), 1000.0, 1.0),
    ("h2-h2o-2000k", ("H2", "H2O"), (0.3, 0.7), 2000.0, 10.0),
    ("co-co2-o2-3000k", ("CO", "CO2", "O2"), (0.3, 0.6, 0.1), 3000.0, 10.0),
    ("h-o-system-3500k", ("H2", "O2", "H2O", "OH", "H", "O"), (0.4, 0.1, 0.5, 0.0, 0.0, 0.0),
     3500.0, 20.0),
)
LIB_UNITS = {
    "temperature": "K", "pressure": "bar", "viscosity": "millipoise",
    "thermal_conductivity": "milliwatt/(cm K)", "specific_heat": "kJ/(kg K)",
    "sonic_velocity": "m/s", "molecular_weight": "g/mol",
}


def _properties(solution, index=None) -> dict:
    """Collect the transport and derivative outputs of one CEA state."""
    def pick(name):
        value = getattr(solution, name)
        return float(value if index is None else value[index])

    names = ("T", "P", "viscosity", "conductivity_fr", "conductivity_eq", "Pr_fr", "Pr_eq",
             "cp_fr", "cp_eq", "cv_fr", "cv_eq", "gamma_s", "M", "MW")
    values = {name: pick(name) for name in names}
    if index is not None:
        values["sonic_velocity"] = pick("sonic_velocity")
    fractions = solution.mole_fractions
    values["mole_fractions"] = {
        key: float(value if index is None else value[index])
        for key, value in sorted(fractions.items())
        if float(value if index is None else value[index]) > 0.0
    }
    return values


def main() -> int:
    import cea
    import numpy as np

    data = Path(cea.__file__).parent / "data"
    digests = {
        name: hashlib.sha256((data / name).read_bytes()).hexdigest()
        for name in ("thermo.lib", "trans.lib")
    }

    def tp_state(species, moles, temperature, pressure):
        products = cea.Mixture(list(species))
        reactants = cea.Mixture(list(species))
        solver = cea.EqSolver(products, reactants=reactants, transport=True)
        solution = cea.EqSolution(solver)
        weights = reactants.moles_to_weights(np.array(moles, dtype=float))
        solver.solve(solution, cea.TP, temperature, pressure, weights)
        if not solution.converged:
            raise SystemExit(f"CEA did not converge for {species} at {temperature} K")
        return _properties(solution)

    pure = []
    for species, temperatures in PURE_SPECIES:
        for temperature in temperatures:
            record = tp_state((species,), (1.0,), temperature, 1.0)
            pure.append({"id": f"pure-{species}-{temperature:g}k", "species": species,
                         "inputs": {"temperature_k": temperature, "pressure_bar": 1.0},
                         "cea": record})
    mixtures = []
    for identifier, species, moles, temperature, pressure in MIXTURES:
        mixtures.append({
            "id": identifier,
            "inputs": {"products": list(species), "reactant_moles": list(moles),
                       "temperature_k": temperature, "pressure_bar": pressure},
            "cea": tp_state(species, moles, temperature, pressure),
        })
    rocket = []
    for pair, fuel, oxidizer, ratios in ROCKET_POINTS:
        reactants = cea.Mixture([fuel, oxidizer])
        products = cea.Mixture([fuel, oxidizer], products_from_reactants=True)
        solver = cea.RocketSolver(products, reactants=reactants, transport=True)
        temperatures = np.array([REACTANT_TEMPERATURE_K[fuel], REACTANT_TEMPERATURE_K[oxidizer]])
        for ratio in ratios:
            solution = cea.RocketSolution(solver)
            weights = reactants.of_ratio_to_weights(np.array([0.0, 1.0]), np.array([1.0, 0.0]), ratio)
            enthalpy = reactants.calc_property(cea.ENTHALPY, weights, temperatures) / cea.R
            solver.solve(solution, weights, CHAMBER_PRESSURE_BAR, supar=[AREA_RATIO], hc=enthalpy,
                         iac=True)
            if not solution.converged:
                raise SystemExit(f"CEA did not converge for {pair} O/F={ratio}")
            for index, station in ((0, "chamber"), (1, "throat")):
                rocket.append({
                    "id": f"cea-{pair.lower()}-of-{ratio:g}-{station}",
                    "pair": pair, "fuel": fuel, "oxidizer": oxidizer,
                    "oxidizer_fuel_ratio": ratio, "station": station,
                    "cea": _properties(solution, index),
                })
    document = {
        "schema_version": 1,
        "generated": datetime.date.today().isoformat(),
        "generator": "scripts/generate_cea_transport_points.py",
        "reference": {
            "program": "NASA CEA (cea Python package)",
            "version": cea.__version__,
            "url": "https://pypi.org/project/cea/",
            "repository": "https://github.com/nasa/cea",
            "thermo_lib_sha256": digests["thermo.lib"],
            "trans_lib_sha256": digests["trans.lib"],
            "universal_gas_constant_j_kmol_k": float(cea.R),
        },
        "cea_output_units": LIB_UNITS,
        "rocket_conditions": {
            "chamber_pressure_bar": CHAMBER_PRESSURE_BAR,
            "area_ratio": AREA_RATIO,
            "combustor": "infinite area (IAC)",
            "expansion": "shifting equilibrium",
            "reactant_temperature_k": REACTANT_TEMPERATURE_K,
        },
        "pure_species": pure,
        "mixtures": mixtures,
        "rocket": rocket,
    }
    TARGET.write_text(json.dumps(document, indent=1) + "\n", encoding="utf-8", newline="\n")
    print(f"wrote {len(pure)} pure, {len(mixtures)} mixture and {len(rocket)} rocket states to {TARGET}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
