"""Optional NASA CEA rocket solver through NASA's ``cea`` Python package.

NASA Glenn publishes CEA as an Apache-2.0 package with a Python API [1]_. When
it is installed (``pip install "rocket-propulsion-lab[cea]"``) this module runs
CEA's infinite-area-combustor rocket problem [2]_ and returns the same
:class:`~rocket_propulsion.thermochemistry.rocket_performance.RocketPerformance`
as the local solver, so the two can be compared station by station. Nothing
here is imported by the core package unless asked for; without the package
:func:`cea_library_info` reports it as unavailable.

References
----------
.. [1] NASA-CEA-PY. NASA Glenn Research Center, "CEA: Chemical Equilibrium with
       Applications", Python package. https://pypi.org/project/cea/
.. [2] NASA RP-1311, "Computer Program for Calculation of Complex Chemical
       Equilibrium Compositions and Applications, Part I: Analysis".
       https://ntrs.nasa.gov/citations/19950013764
"""

from __future__ import annotations

import hashlib
from collections.abc import Sequence
from dataclasses import dataclass
from functools import lru_cache
from math import isfinite
from pathlib import Path
from typing import Any

from rocket_propulsion.core.constants import STANDARD_GRAVITY_M_S2
from rocket_propulsion.core.errors import ConvergenceError, DomainError, FeatureUnavailableError
from rocket_propulsion.core.metadata import WarningMessage

from .rocket_performance import RocketExit, RocketPerformance, RocketStation
from .species import CondensedReactant, get_reactant


@dataclass(frozen=True, slots=True)
class CeaLibraryInfo:
    """Availability and provenance of the installed NASA CEA package [1]_.

    Attributes
    ----------
    available : bool
        Whether ``import cea`` succeeds.
    version : str or None
        Package version.
    thermo_sha256 : str or None
        SHA-256 of the packaged ``thermo.lib`` database.
    detail : str
        Installation hint or import error.

    References
    ----------
    .. [1] NASA-CEA-PY. NASA Glenn Research Center, "CEA: Chemical Equilibrium
           with Applications", Python package. https://pypi.org/project/cea/
    """

    available: bool
    version: str | None
    thermo_sha256: str | None
    detail: str


def _import_cea() -> Any:
    import cea  # noqa: PLC0415 - optional dependency

    return cea


@lru_cache(maxsize=1)
def cea_library_info() -> CeaLibraryInfo:
    """Report whether NASA's ``cea`` package is importable, with provenance [1]_.

    Returns
    -------
    CeaLibraryInfo
        Version and database hash when available; the import error otherwise.

    References
    ----------
    .. [1] NASA-CEA-PY. NASA Glenn Research Center, "CEA: Chemical Equilibrium
           with Applications", Python package. https://pypi.org/project/cea/
    """

    try:
        cea = _import_cea()
    except Exception as error:  # noqa: BLE001 - any import failure means unavailable
        return CeaLibraryInfo(False, None, None,
                              f"Install with: pip install \"rocket-propulsion-lab[cea]\" ({error})")
    database = Path(cea.__file__).parent / "data" / "thermo.lib"
    digest = hashlib.sha256(database.read_bytes()).hexdigest() if database.is_file() else None
    return CeaLibraryInfo(True, str(cea.__version__), digest, str(Path(cea.__file__).parent))


def _station(name: str, area_ratio: float | None, solution: Any, index: int) -> RocketStation:
    fractions = {
        str(key): float(values[index])
        for key, values in solution.mole_fractions.items()
        if float(values[index]) >= 1.0e-5
    }
    return RocketStation(
        name=name,
        area_ratio=area_ratio,
        pressure_pa=float(solution.P[index]) * 1.0e5,
        temperature_k=float(solution.T[index]),
        density_kg_m3=float(solution.density[index]),
        molecular_mass_kg_mol=float(solution.MW[index]) / 1000.0,
        enthalpy_j_kg=float(solution.enthalpy[index]) * 1000.0,
        entropy_j_kg_k=float(solution.entropy[index]) * 1000.0,
        velocity_m_s=float(solution.Mach[index]) * float(solution.sonic_velocity[index]),
        mole_fractions=tuple(sorted(fractions.items(), key=lambda item: -item[1])),
    )


def solve_cea_rocket(
    *,
    fuel: str,
    oxidizer: str,
    oxidizer_fuel_ratio: float,
    chamber_pressure_pa: float,
    area_ratios: Sequence[float] = (10.0, 25.0, 50.0),
    chemistry: str = "equilibrium",
) -> RocketPerformance:
    """Run NASA CEA's infinite-area-combustor rocket problem [1]_, [2]_.

    Parameters
    ----------
    fuel, oxidizer : str
        Bundled reactant keys, which are also CEA names (``"RP-1"``,
        ``"O2(L)"``). Condensed reactants enter at their assigned
        temperature; gas reactants at 298.15 K.
    oxidizer_fuel_ratio : float
        Oxidizer-to-fuel mass ratio.
    chamber_pressure_pa : float
        Chamber pressure in Pa.
    area_ratios : sequence of float, default (10, 25, 50)
        Supersonic exit-to-throat area ratios.
    chemistry : {"equilibrium", "frozen"}, default "equilibrium"
        Shifting equilibrium, or frozen at the chamber (``n_frz=1``).

    Returns
    -------
    RocketPerformance
        Result with ``provider="nasa-cea"`` and the CEA version.

    Raises
    ------
    FeatureUnavailableError
        If the ``cea`` package is not installed.
    DomainError
        If the ratio, pressure or an area ratio is out of range.
    ConvergenceError
        If CEA reports a non-converged solution.

    References
    ----------
    .. [1] NASA-CEA-PY. NASA Glenn Research Center, "CEA: Chemical Equilibrium
           with Applications", Python package. https://pypi.org/project/cea/
    .. [2] NASA RP-1311, "Computer Program for Calculation of Complex Chemical
           Equilibrium Compositions and Applications, Part I: Analysis".
           https://ntrs.nasa.gov/citations/19950013764
    """

    info = cea_library_info()
    if not info.available:
        raise FeatureUnavailableError("NASA CEA Python package is not installed.",
                                      details={"install": info.detail})
    if not isfinite(oxidizer_fuel_ratio) or oxidizer_fuel_ratio <= 0.0:
        raise DomainError("O/F ratio must be positive.", field="oxidizer_fuel_ratio")
    if not isfinite(chamber_pressure_pa) or chamber_pressure_pa <= 0.0:
        raise DomainError("Chamber pressure must be positive.", field="chamber_pressure_pa")
    ratios = [float(value) for value in area_ratios]
    if not ratios or any(not isfinite(value) or value <= 1.0 for value in ratios):
        raise DomainError("Rocket area ratios must all exceed one.", field="area_ratios")
    if chemistry not in ("equilibrium", "frozen"):
        raise DomainError("Rocket chemistry must be equilibrium or frozen.", field="chemistry")
    cea = _import_cea()
    import numpy  # noqa: PLC0415 - installed with cea

    records = [get_reactant(fuel), get_reactant(oxidizer)]
    names = [record.key for record in records]
    temperatures = numpy.array([
        record.temperature_k if isinstance(record, CondensedReactant) else 298.15
        for record in records
    ])
    reactants = cea.Mixture(names)
    products = cea.Mixture(names, products_from_reactants=True)
    solver = cea.RocketSolver(products, reactants=reactants)
    solution = cea.RocketSolution(solver)
    weights = reactants.of_ratio_to_weights(numpy.array([0.0, 1.0]), numpy.array([1.0, 0.0]),
                                            oxidizer_fuel_ratio)
    enthalpy = reactants.calc_property(cea.ENTHALPY, weights, temperatures) / cea.R
    options: dict[str, Any] = {"supar": ratios, "hc": enthalpy, "iac": True}
    if chemistry == "frozen":
        options["n_frz"] = 1
    solver.solve(solution, weights, chamber_pressure_pa / 1.0e5, **options)
    if not solution.converged:
        raise ConvergenceError("NASA CEA did not converge.", code="cea_not_converged")
    exits = []
    for offset, ratio in enumerate(ratios):
        index = 2 + offset
        station = _station("exit", ratio, solution, index)
        vacuum = float(solution.Isp_vacuum[index])
        exits.append(RocketExit(
            area_ratio=ratio,
            station=station,
            effective_exhaust_velocity_m_s=float(solution.Isp[index]),
            vacuum_effective_exhaust_velocity_m_s=vacuum,
            vacuum_specific_impulse_s=vacuum / STANDARD_GRAVITY_M_S2,
            vacuum_thrust_coefficient=float(solution.coefficient_of_thrust[index])
            * vacuum / float(solution.Isp[index]),
        ))
    fuel_record, oxidizer_record = records
    oxidizer_moles = oxidizer_fuel_ratio * fuel_record.molecular_mass_kg_mol / oxidizer_record.molecular_mass_kg_mol
    return RocketPerformance(
        provider="nasa-cea",
        provider_version=str(info.version),
        chemistry=chemistry,
        reactants=tuple(sorted(((fuel_record.key.upper(), 1.0), (oxidizer_record.key.upper(), oxidizer_moles)))),
        chamber_pressure_pa=chamber_pressure_pa,
        chamber=_station("chamber", None, solution, 0),
        throat=_station("throat", 1.0, solution, 1),
        exits=tuple(exits),
        characteristic_velocity_m_s=float(solution.c_star[-1]),
        warnings=(WarningMessage(
            code="nasa_cea_reference",
            message=f"NASA CEA {info.version}; thermo.lib sha256 {info.thermo_sha256}.",
        ),),
    )
