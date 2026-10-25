"""HTTP routes for CEA-equivalent rocket performance and chemistry-driven nozzles.

Handlers only translate JSON to library calls; the physics lives in
:mod:`rocket_propulsion.thermochemistry.rocket_performance`,
:mod:`rocket_propulsion.thermochemistry.cea_library` and
:mod:`rocket_propulsion.thermochemistry.chemistry_nozzle`.

References
----------
.. [1] NASA RP-1311, "Computer Program for Calculation of Complex Chemical
       Equilibrium Compositions and Applications, Part I: Analysis".
       https://ntrs.nasa.gov/citations/19950013764
"""

from __future__ import annotations

from dataclasses import asdict
from math import isfinite
from typing import Any

from rocket_propulsion.core.constants import USSA1976_SEA_LEVEL_PRESSURE
from rocket_propulsion.core.errors import InputError
from rocket_propulsion.core.metadata import ModelMetadata
from rocket_propulsion.thermochemistry import (
    cea_library_info,
    get_propellant_pair,
    reactants_for_mass_ratio,
    solve_cea_rocket,
    solve_rocket_performance,
)
from rocket_propulsion.thermochemistry.chemistry_nozzle import solve_chemistry_coupled_nozzle

MAXIMUM_AREA_RATIOS = 8


def _number(payload: dict[str, Any], key: str, default: float) -> float:
    raw = payload.get(key, default)
    if raw in (None, ""):
        raw = default
    if isinstance(raw, bool):
        raise InputError(f"Field {key} must be a number.", field=key)
    try:
        value = float(raw)
    except (TypeError, ValueError, OverflowError) as error:
        raise InputError(f"Field {key} must be a number.", field=key) from error
    if not isfinite(value):
        raise InputError(f"Field {key} must be finite.", field=key)
    return value


def _area_ratios(payload: dict[str, Any]) -> tuple[float, ...]:
    raw = payload.get("area_ratios", [5.0, 10.0, 25.0, 40.0, 80.0])
    if isinstance(raw, str):
        raw = [part for part in raw.replace(";", ",").split(",") if part.strip()]
    if not isinstance(raw, list) or not raw:
        raise InputError("Field area_ratios must be a non-empty list.", field="area_ratios")
    if len(raw) > MAXIMUM_AREA_RATIOS:
        raise InputError(f"At most {MAXIMUM_AREA_RATIOS} area ratios per request.", field="area_ratios")
    values = []
    for item in raw:
        values.append(_number({"value": item}, "value", 0.0))
    return tuple(sorted(set(values)))


def _relative(left: float, right: float) -> float:
    return left / right - 1.0


def rocket_performance_route(payload: dict[str, Any]) -> dict[str, Any]:
    """Solve the CEA-equivalent infinite-area-combustor rocket problem [1]_, [2]_.

    Parameters
    ----------
    payload : dict
        propellant_pair : str, default "LOX_LH2"
            Preset key.
        oxidizer_fuel_ratio : float, optional
            Mass ratio; the preset reference ratio by default.
        chamber_pressure_pa : float, default 7.0e6
            Chamber pressure, Pa.
        area_ratios : list of float, default [5, 10, 25, 40, 80]
            Up to eight supersonic area ratios.
        chemistry : {"equilibrium", "frozen"}, default "equilibrium"
            Expansion chemistry.
        provider : {"local", "cea", "compare"}, default "local"
            Local solver, NASA CEA (optional package) or both with relative
            differences.

    Returns
    -------
    dict
        ``result`` (local or CEA), and for ``compare`` also ``reference``
        and ``differences``; ``cea`` availability and the pair preset.

    Raises
    ------
    InputError
        For an unknown provider, pair or malformed field.
    DomainError
        For values outside the model domain.
    FeatureUnavailableError
        If ``cea`` or ``compare`` is requested without the package.

    References
    ----------
    .. [1] NASA RP-1311, "Computer Program for Calculation of Complex Chemical
           Equilibrium Compositions and Applications, Part I: Analysis".
           https://ntrs.nasa.gov/citations/19950013764
    .. [2] NASA-CEA-PY. NASA Glenn Research Center, "CEA: Chemical Equilibrium
           with Applications", Python package. https://pypi.org/project/cea/
    """

    pair = get_propellant_pair(str(payload.get("propellant_pair", "LOX_LH2")))
    ratio = _number(payload, "oxidizer_fuel_ratio", pair.reference_oxidizer_fuel_ratio)
    pressure = _number(payload, "chamber_pressure_pa", 7.0e6)
    areas = _area_ratios(payload)
    chemistry = str(payload.get("chemistry", "equilibrium")).strip().lower()
    provider = str(payload.get("provider", "local")).strip().lower()
    if provider not in {"local", "cea", "compare"}:
        raise InputError("Provider must be local, cea or compare.", field="provider",
                         details={"supported": ["local", "cea", "compare"]})
    info = cea_library_info()
    response: dict[str, Any] = {
        "pair": asdict(pair),
        "oxidizer_fuel_ratio": ratio,
        "cea": {"available": info.available, "version": info.version,
                "thermo_sha256": info.thermo_sha256},
    }

    def cea() -> Any:
        return solve_cea_rocket(fuel=pair.fuel_reactant, oxidizer=pair.oxidizer_reactant,
                                oxidizer_fuel_ratio=ratio, chamber_pressure_pa=pressure,
                                area_ratios=areas, chemistry=chemistry)

    if provider == "cea":
        response["result"] = asdict(cea())
        return response
    local = solve_rocket_performance(reactants=reactants_for_mass_ratio(pair, ratio),
                                     chamber_pressure_pa=pressure, area_ratios=areas,
                                     chemistry=chemistry)  # type: ignore[arg-type]
    response["result"] = asdict(local)
    if provider == "compare":
        reference = cea()
        response["reference"] = asdict(reference)
        response["differences"] = {
            "chamber_temperature": _relative(local.chamber.temperature_k, reference.chamber.temperature_k),
            "characteristic_velocity": _relative(local.characteristic_velocity_m_s,
                                                 reference.characteristic_velocity_m_s),
            "vacuum_specific_impulse": [
                _relative(mine.vacuum_specific_impulse_s, theirs.vacuum_specific_impulse_s)
                for mine, theirs in zip(local.exits, reference.exits, strict=True)
            ],
        }
    return response


def chemistry_nozzle_route(payload: dict[str, Any]) -> dict[str, Any]:
    """Solve nozzle regimes with properties calibrated to propellant chemistry [1]_.

    Parameters
    ----------
    payload : dict
        propellant_pair : str, default "LOX_RP1"
        oxidizer_fuel_ratio : float, optional
        chamber_pressure_pa : float, default 7.0e6
        area_ratio : float, default 16.0
        ambient_pressure_pa : float, default 101325.0
        throat_area_m2 : float, default 0.01
        maximum_altitude_m : float, default 50000.0
        point_count : int, default 26

    Returns
    -------
    dict
        Equilibrium design point, calibrated gas, regime solution and a
        standard-atmosphere sweep.

    Raises
    ------
    InputError
        For an unknown pair or malformed field.
    DomainError
        For values outside the model domain.

    References
    ----------
    .. [1] NASA RP-1311, "Computer Program for Calculation of Complex Chemical
           Equilibrium Compositions and Applications, Part I: Analysis".
           https://ntrs.nasa.gov/citations/19950013764
    """

    pair = get_propellant_pair(str(payload.get("propellant_pair", "LOX_RP1")))
    point_count = _number(payload, "point_count", 26.0)
    if not point_count.is_integer():
        raise InputError("Field point_count must be an integer.", field="point_count")
    result = solve_chemistry_coupled_nozzle(
        propellant_pair=pair.key,
        oxidizer_fuel_ratio=_number(payload, "oxidizer_fuel_ratio", pair.reference_oxidizer_fuel_ratio),
        chamber_pressure_pa=_number(payload, "chamber_pressure_pa", 7.0e6),
        area_ratio=_number(payload, "area_ratio", 16.0),
        ambient_pressure_pa=_number(payload, "ambient_pressure_pa", USSA1976_SEA_LEVEL_PRESSURE.value),
        throat_area_m2=_number(payload, "throat_area_m2", 0.01),
        include_altitude_sweep=True,
        maximum_altitude_m=_number(payload, "maximum_altitude_m", 50_000.0),
        point_count=int(point_count),
    )
    return asdict(result)


CHEMISTRY_ROUTES = {
    "/api/v1/rocket-performance": rocket_performance_route,
    "/api/v1/nozzle-chemistry": chemistry_nozzle_route,
}

CHEMISTRY_ROUTE_METADATA = {
    "/api/v1/rocket-performance": ModelMetadata(
        model="infinite-area-combustor rocket performance",
        version="1.0",
        assumptions=(
            "Assigned-enthalpy chamber equilibrium",
            "One-dimensional isentropic expansion, shifting or frozen chemistry",
            "Throat at maximum mass flux; ideal gas products",
            "Gated against NASA CEA 3.3.4 at nine points and 0.5-25 MPa",
        ),
        units={"pressure": "Pa", "temperature": "K", "velocity": "m/s", "specific_impulse": "s"},
        references=("NASA RP-1311", "NASA CEA 3.3.4 (cea Python package)"),
    ),
    "/api/v1/nozzle-chemistry": ModelMetadata(
        model="chemistry-calibrated off-design nozzle",
        version="1.0",
        assumptions=(
            "Perfect gas calibrated to the shifting-equilibrium design point",
            "Mass flow and design exit pressure match equilibrium; other stations approximate",
            "Summerfield separation screen; external waves not resolved",
        ),
        units={"pressure": "Pa", "temperature": "K", "thrust": "N", "specific_impulse": "s"},
        references=("NASA RP-1311", "NACA Report 1135", "NASA SP-8120 p.17"),
    ),
}
