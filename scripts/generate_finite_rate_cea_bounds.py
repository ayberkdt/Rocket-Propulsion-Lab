"""Record NASA CEA frozen and shifting-equilibrium bounds for the finite-rate nozzle case.

The finite-rate H2-O2 nozzle test case of ``tests/test_finite_rate_kinetics.py``
(LOX/LH2, O/F 6.0, 1000 psia, A_e/A_t = 40) must lie between frozen and
shifting-equilibrium performance. This script runs NASA's own CEA (the ``cea``
Python package) for that case with the product set restricted to the eight
species of the GRI-Mech 3.0 H2-O2 subset, and prints the JSON block stored
under ``cea_bounds`` in ``tests/reference_data/finite_rate_kinetics.json``.
Frozen composition at the chamber is ``n_frz=1`` and at the throat
``n_frz=2`` (CEA's RP-1311 example 12 uses ``n_frz=2`` for "frozen composition
from the throat").

Usage (in an environment with ``pip install cea``)::

    python scripts/generate_finite_rate_cea_bounds.py

References
----------
.. [1] NASA RP-1311, "Computer Program for Calculation of Complex Chemical
       Equilibrium Compositions and Applications, Part I: Analysis".
       https://ntrs.nasa.gov/citations/19950013764
.. [2] NASA-CEA-PY. NASA Glenn Research Center, "CEA: Chemical Equilibrium with
       Applications", Python package 3.3.4. https://pypi.org/project/cea/
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

CHAMBER_PRESSURE_BAR = 68.9476
AREA_RATIO = 40.0
OXIDIZER_FUEL_RATIO = 6.0
REACTANT_TEMPERATURE_K = {"H2(L)": 20.27, "O2(L)": 90.17}
PRODUCTS = ["H2", "H", "O", "O2", "OH", "H2O", "HO2", "H2O2"]
G0 = 9.80665


def main() -> int:
    import cea
    import numpy as np

    data = Path(cea.__file__).parent / "data" / "thermo.lib"
    reactants = cea.Mixture(["H2(L)", "O2(L)"])
    products = cea.Mixture(PRODUCTS)
    solver = cea.RocketSolver(products, reactants=reactants)
    temperatures = np.array([REACTANT_TEMPERATURE_K["H2(L)"], REACTANT_TEMPERATURE_K["O2(L)"]])
    weights = reactants.of_ratio_to_weights(np.array([0.0, 1.0]), np.array([1.0, 0.0]),
                                            OXIDIZER_FUEL_RATIO)
    enthalpy = reactants.calc_property(cea.ENTHALPY, weights, temperatures) / cea.R
    results = {}
    for label, frozen_at in (("shifting_equilibrium", None), ("frozen_at_chamber", 1),
                             ("frozen_at_throat", 2)):
        solution = cea.RocketSolution(solver)
        solver.solve(solution, weights, CHAMBER_PRESSURE_BAR, supar=[AREA_RATIO], iac=True,
                     hc=enthalpy, n_frz=frozen_at)
        if not solution.converged:
            raise SystemExit(f"CEA did not converge for {label}")
        results[label] = {
            "vacuum_specific_impulse_s": float(solution.Isp_vacuum[-1]) / G0,
            "characteristic_velocity_m_s": float(solution.c_star[-1]),
            "exit_temperature_k": float(solution.T[-1]),
            "chamber_temperature_k": float(solution.T[0]),
        }
    block = {
        "program": "NASA CEA (cea Python package)",
        "version": cea.__version__,
        "thermo_lib_sha256": hashlib.sha256(data.read_bytes()).hexdigest(),
        "generator": "scripts/generate_finite_rate_cea_bounds.py",
        "products": PRODUCTS,
        "reactant_temperature_k": REACTANT_TEMPERATURE_K,
        "chamber_pressure_pa": CHAMBER_PRESSURE_BAR * 1.0e5,
        "oxidizer_fuel_ratio": OXIDIZER_FUEL_RATIO,
        "area_ratio": AREA_RATIO,
        "results": results,
    }
    print(json.dumps(block, indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
