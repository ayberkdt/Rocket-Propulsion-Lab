"""JSON request adapters that connect the HTTP API to the calculation domain.

Each ``*_route`` function takes a decoded JSON request body, converts its
fields with strict number, integer and boolean readers, calls one domain
calculation and returns a JSON-ready dictionary. ``ROUTES`` maps versioned
API paths to these adapters and ``ROUTE_METADATA`` records the model,
assumptions, units and source labels returned with every response. The
physical models and their primary sources are documented on the domain
functions each route calls.
"""

from __future__ import annotations

import json
import os
from dataclasses import asdict
from decimal import Decimal, InvalidOperation
from fractions import Fraction
from math import isfinite, log
from pathlib import Path
from typing import Any

from rocket_propulsion.compressible import (
    calculate_fanno,
    calculate_ideal_nozzle,
    calculate_isentropic,
    calculate_normal_shock,
    calculate_oblique_shock,
    calculate_off_design_nozzle,
    calculate_rayleigh,
    generate_isentropic_sweep,
    generate_nozzle_altitude_sweep,
    mach_from_isentropic,
)
from rocket_propulsion.core.constants import (
    AIR_SPECIFIC_GAS_CONSTANT_J_KG_K,
    STANDARD_GRAVITY_M_S2,
)
from rocket_propulsion.core.errors import DomainError, InputError
from rocket_propulsion.core.metadata import ModelMetadata, ValidityRange
from rocket_propulsion.geometry import (
    NozzleContour,
    compare_nozzle_contours,
    export_nozzle_contour,
    generate_moc_nozzle,
    generate_nozzle_contour,
    generate_rao_contour,
)
from rocket_propulsion.integrations import export_sidera_artifact, inspect_sidera_capabilities
from rocket_propulsion.propulsion import (
    BurnDefinition,
    BurnTarget,
    BurnTargetKind,
    ConstantPerformanceProvider,
    EngineOperatingPoint,
    FlowDestination,
    ProfiledPerformanceProvider,
    PropellantFlow,
    PropellantRole,
    StageDefinition,
    TankDefinition,
    ThrottleSchedule,
    ThrottleSegment,
    build_burn_html_report,
    burn_result_to_dict,
    calculate_engine_cycle,
    calculate_loss_budget,
    calculate_rocket_equation,
    calculate_thrust,
    calculate_two_phase_band,
    list_propellants,
    optimize_staging,
    simulate_constant_burn,
    simulate_profiled_burn,
)
from rocket_propulsion.propulsion.burns import burn_result_from_dict
from rocket_propulsion.propulsion.burns.numeric_outputs import _finite_rational
from rocket_propulsion.studies import (
    ParameterRange,
    build_report_bundle,
    compare_cases,
    propagate_uncertainty,
    workspace_from_json,
    workspace_to_json,
)
from rocket_propulsion.thermochemistry import (
    PROPELLANT_PAIRS,
    CeaEquilibriumProvider,
    CeaSubprocessAdapter,
    EquilibriumProblem,
    LocalEquilibriumProvider,
    Mixture,
    cea_library_info,
    generate_oxidizer_fuel_sweep,
    get_species,
    species_catalog,
)
from rocket_propulsion.thermodynamics import (
    calculate_bartz_heat_transfer,
    calculate_isentropic_process,
    calculate_polytropic_process,
    calculate_regenerative_cooling,
    calculate_stagnation_state,
    complete_ideal_gas_state,
)

from .aero_routes import AERO_ROUTE_METADATA, AERO_ROUTES
from .capabilities import SideraCapabilityProbe
from .chemistry_routes import CHEMISTRY_ROUTE_METADATA, CHEMISTRY_ROUTES
from .thermal_routes import THERMAL_ROUTE_METADATA, THERMAL_ROUTES


def _number(payload: dict[str, Any], key: str, *, default: float | None = None) -> float:
    """Read a required finite number from a request payload.

    ``bool`` is refused even though it converts to a number. Text that
    ``float`` accepts (for example form values) is converted. ``default`` is
    used when the key is absent; a missing value without a default raises
    :class:`InputError` with code ``"missing_field"``.
    """

    raw_value = payload.get(key, default)
    if raw_value is None:
        raise InputError(f"Missing required field: {key}.", code="missing_field", field=key)
    if isinstance(raw_value, bool):
        raise InputError(f"Field {key} must be a number.", field=key)
    try:
        value = float(raw_value)
    except (TypeError, ValueError, OverflowError) as error:
        raise InputError(f"Field {key} must be a number.", field=key) from error
    if not isfinite(value):
        raise InputError(f"Field {key} must be finite.", field=key)
    return value


def _optional_number(payload: dict[str, Any], key: str) -> float | None:
    """Read an optional finite number; ``None`` or an empty string means absent."""

    raw_value = payload.get(key)
    if raw_value in (None, ""):
        return None
    return _number(payload, key)


def _integer(payload: dict[str, Any], key: str, *, default: int) -> int:
    """Read an integral number as ``int``; fractional values are refused, not truncated."""

    value = _number(payload, key, default=float(default))
    if not value.is_integer():
        raise InputError(f"Field {key} must be an integer.", field=key)
    return int(value)


def _boolean(payload: dict[str, Any], key: str, *, default: bool) -> bool:
    """Read an explicit boolean from JSON or form text.

    JSON booleans are returned as is. Text is matched case-insensitively
    against ``true/1/yes/on`` and ``false/0/no/off``; anything else raises
    :class:`InputError`, so Python truthiness (e.g. of ``"false"``) never
    applies.
    """

    raw_value = payload.get(key, default)
    if isinstance(raw_value, bool):
        return raw_value
    normalized = str(raw_value).strip().lower()
    if normalized in {"true", "1", "yes", "on"}:
        return True
    if normalized in {"false", "0", "no", "off"}:
        return False
    raise InputError(f"Field {key} must be true or false.", field=key)


def calculate_isentropic_route(payload: dict[str, Any]) -> dict[str, Any]:
    """Resolve one isentropic property to Mach number and return all isentropic ratios.

    Parameters
    ----------
    payload : dict
        Decoded JSON request body with the fields:

        input_kind : str, default "mach"
            Which property ``value`` gives (as accepted by
            :func:`~rocket_propulsion.compressible.mach_from_isentropic`).
        value : float
            Value of that property.
        gamma : float, default 1.4
            Ratio of specific heats.

    Returns
    -------
    dict
        Fields of the isentropic result for the resolved Mach number.

    Raises
    ------
    InputError
        If a required field is missing or a numeric field is not a finite
        number.
    DomainError
        If the values are outside the model domain.

    Notes
    -----
    Perfect-gas, adiabatic and reversible one-dimensional flow [1]_.

    References
    ----------
    .. [1] NACA-TR-1135. Ames Research Staff, "Equations, Tables, and Charts for
           Compressible Flow", NACA Report 1135, 1953,
           https://ntrs.nasa.gov/citations/19930091059
    """

    input_kind = str(payload.get("input_kind", "mach"))
    gamma = _number(payload, "gamma", default=1.4)
    mach = mach_from_isentropic(input_kind, _number(payload, "value"), gamma)
    return asdict(calculate_isentropic(mach, gamma))


def calculate_normal_shock_route(payload: dict[str, Any]) -> dict[str, Any]:
    """Return normal-shock properties for a supersonic upstream Mach number.

    Parameters
    ----------
    payload : dict
        Decoded JSON request body with the fields:

        upstream_mach : float
            Upstream Mach number; must exceed one.
        gamma : float, default 1.4
            Ratio of specific heats.

    Returns
    -------
    dict
        Fields of the normal-shock result.

    Raises
    ------
    InputError
        If a required field is missing or a numeric field is not a finite
        number.
    DomainError
        If the values are outside the model domain.

    Notes
    -----
    Perfect-gas Rankine-Hugoniot relations for a planar normal shock [1]_.

    References
    ----------
    .. [1] NACA-TR-1135. Ames Research Staff, "Equations, Tables, and Charts for
           Compressible Flow", NACA Report 1135, 1953,
           https://ntrs.nasa.gov/citations/19930091059
    """

    gamma = _number(payload, "gamma", default=1.4)
    return asdict(calculate_normal_shock(_number(payload, "upstream_mach"), gamma))


def calculate_oblique_shock_route(payload: dict[str, Any]) -> dict[str, Any]:
    """Solve the weak or strong attached oblique-shock relations.

    Parameters
    ----------
    payload : dict
        Decoded JSON request body with the fields:

        upstream_mach : float
            Upstream Mach number; must exceed one.
        turning_angle_deg : float
            Flow deflection angle, deg.
        gamma : float, default 1.4
            Ratio of specific heats.
        branch : {"weak", "strong"}, default "weak"
            Shock-angle branch of the theta-beta-M relation.

    Returns
    -------
    dict
        Fields of the oblique-shock result.

    Raises
    ------
    InputError
        If a required field is missing or a numeric field is not a finite
        number.
    DomainError
        If the values are outside the model domain.

    Notes
    -----
    Planar attached shock from the theta-beta-M relation [1]_.

    References
    ----------
    .. [1] NACA-TR-1135. Ames Research Staff, "Equations, Tables, and Charts for
           Compressible Flow", NACA Report 1135, 1953,
           https://ntrs.nasa.gov/citations/19930091059
    """

    return asdict(
        calculate_oblique_shock(
            upstream_mach=_number(payload, "upstream_mach"),
            turning_angle_deg=_number(payload, "turning_angle_deg"),
            gamma=_number(payload, "gamma", default=1.4),
            branch=str(payload.get("branch", "weak")),
        )
    )


def calculate_isentropic_sweep_route(payload: dict[str, Any]) -> dict[str, Any]:
    """Return sampled isentropic curves for plotting.

    Parameters
    ----------
    payload : dict
        Decoded JSON request body with the fields:

        gamma : float, default 1.4
            Ratio of specific heats.
        minimum_mach, maximum_mach : float, default 0.05 and 5.0
            Sampled Mach-number interval.
        point_count : int, default 121
            Number of samples.

    Returns
    -------
    dict
        ``{"points": [...]}`` with one isentropic record per sample.

    Raises
    ------
    InputError
        If a required field is missing or a numeric field is not a finite
        number.
    DomainError
        If the values are outside the model domain.

    References
    ----------
    .. [1] NACA-TR-1135. Ames Research Staff, "Equations, Tables, and Charts for
           Compressible Flow", NACA Report 1135, 1953,
           https://ntrs.nasa.gov/citations/19930091059
    """

    points = generate_isentropic_sweep(
        gamma=_number(payload, "gamma", default=1.4),
        minimum_mach=_number(payload, "minimum_mach", default=0.05),
        maximum_mach=_number(payload, "maximum_mach", default=5.0),
        point_count=_integer(payload, "point_count", default=121),
    )
    return {"points": [asdict(point) for point in points]}


def calculate_ideal_gas_route(payload: dict[str, Any]) -> dict[str, Any]:
    """Complete an ideal-gas state from exactly two state properties.

    Parameters
    ----------
    payload : dict
        Decoded JSON request body with the fields:

        pressure_pa, temperature_k, density_kg_m3 : float, optional
            Exactly two must be given; ``null`` or ``""`` means absent.
        gas_constant_j_kg_k : float, optional
            Specific gas constant, J/(kg K); defaults to
            ``AIR_SPECIFIC_GAS_CONSTANT_J_KG_K``.
        gamma : float, default 1.4
            Ratio of specific heats.

    Returns
    -------
    dict
        Fields of :class:`~rocket_propulsion.thermodynamics.IdealGasState`.

    Raises
    ------
    InputError
        If a required field is missing or a numeric field is not a finite
        number.
    DomainError
        If the values are outside the model domain.

    Notes
    -----
    ``p = rho R T`` [1]_ with constant ``cp`` and ``cv`` [2]_.

    References
    ----------
    .. [1] NASA-GAS. NASA Glenn Research Center, "Equation of State (Ideal Gas)",
           https://www1.grc.nasa.gov/beginners-guide-to-aeronautics/equation-of-state-ideal-gas-2/
    .. [2] NASA Glenn Research Center, "Specific Heats",
           https://www.grc.nasa.gov/www/k-12/airplane/specheat.html
    """

    state = complete_ideal_gas_state(
        pressure_pa=_optional_number(payload, "pressure_pa"),
        temperature_k=_optional_number(payload, "temperature_k"),
        density_kg_m3=_optional_number(payload, "density_kg_m3"),
        gas_constant_j_kg_k=_number(payload, "gas_constant_j_kg_k", default=AIR_SPECIFIC_GAS_CONSTANT_J_KG_K),
        gamma=_number(payload, "gamma", default=1.4),
    )
    return asdict(state)


def calculate_stagnation_route(payload: dict[str, Any]) -> dict[str, Any]:
    """Connect a static flow state to its stagnation (total) state.

    Parameters
    ----------
    payload : dict
        Decoded JSON request body with the fields:

        mach : float
            Local Mach number.
        static_pressure_pa : float
            Static pressure, Pa.
        static_temperature_k : float
            Static temperature, K.
        gamma : float, default 1.4
            Ratio of specific heats.
        gas_constant_j_kg_k : float, optional
            Specific gas constant, J/(kg K); defaults to the registered dry-air
            value.

    Returns
    -------
    dict
        Fields of :class:`~rocket_propulsion.thermodynamics.StagnationState`.

    Raises
    ------
    InputError
        If a required field is missing or a numeric field is not a finite
        number.
    DomainError
        If the values are outside the model domain.

    Notes
    -----
    Adiabatic, reversible perfect-gas deceleration [1]_.

    References
    ----------
    .. [1] NACA-TR-1135. Ames Research Staff, "Equations, Tables, and Charts for
           Compressible Flow", NACA Report 1135, 1953,
           https://ntrs.nasa.gov/citations/19930091059
    """

    return asdict(
        calculate_stagnation_state(
            mach=_number(payload, "mach"),
            static_pressure_pa=_number(payload, "static_pressure_pa"),
            static_temperature_k=_number(payload, "static_temperature_k"),
            gamma=_number(payload, "gamma", default=1.4),
            gas_constant_j_kg_k=_number(payload, "gas_constant_j_kg_k", default=AIR_SPECIFIC_GAS_CONSTANT_J_KG_K),
        )
    )


def calculate_process_route(payload: dict[str, Any]) -> dict[str, Any]:
    """Calculate ideal and efficiency-corrected compression or expansion.

    Parameters
    ----------
    payload : dict
        Decoded JSON request body with the fields:

        inlet_pressure_pa, outlet_pressure_pa : float
            End pressures, Pa.
        inlet_temperature_k : float
            Inlet temperature, K.
        efficiency : float, default 1.0
            Isentropic efficiency in (0, 1].
        gamma : float, default 1.4
            Ratio of specific heats.
        gas_constant_j_kg_k : float, optional
            Specific gas constant, J/(kg K); defaults to the registered dry-air
            value.

    Returns
    -------
    dict
        Fields of :class:`~rocket_propulsion.thermodynamics.IsentropicProcessResult`.

    Raises
    ------
    InputError
        If a required field is missing or a numeric field is not a finite
        number.
    DomainError
        If the values are outside the model domain.

    Notes
    -----
    Isentropic temperature ratio [1]_ with compressor or turbine efficiency
    definitions [2]_.

    References
    ----------
    .. [1] NACA-TR-1135. Ames Research Staff, "Equations, Tables, and Charts for
           Compressible Flow", NACA Report 1135, 1953,
           https://ntrs.nasa.gov/citations/19930091059
    .. [2] NASA-ENGINESIM-THEORY. NASA Glenn Research Center, ideal gas turbine
           temperature/pressure and shaft-work relations,
           https://www.grc.nasa.gov/www/k-12/airplane/EngineTheory.pdf
    """

    return asdict(
        calculate_isentropic_process(
            inlet_pressure_pa=_number(payload, "inlet_pressure_pa"),
            outlet_pressure_pa=_number(payload, "outlet_pressure_pa"),
            inlet_temperature_k=_number(payload, "inlet_temperature_k"),
            efficiency=_number(payload, "efficiency", default=1.0),
            gamma=_number(payload, "gamma", default=1.4),
            gas_constant_j_kg_k=_number(payload, "gas_constant_j_kg_k", default=AIR_SPECIFIC_GAS_CONSTANT_J_KG_K),
        )
    )


def calculate_polytropic_route(payload: dict[str, Any]) -> dict[str, Any]:
    """Return a polytropic process end state, first-law terms and diagram path.

    Parameters
    ----------
    payload : dict
        Decoded JSON request body with the fields:

        inlet_pressure_pa, outlet_pressure_pa : float
            End pressures, Pa.
        inlet_temperature_k : float
            Inlet temperature, K.
        exponent : float
            Polytropic exponent ``n``.
        gamma : float, default 1.4
            Ratio of specific heats.
        gas_constant_j_kg_k : float, optional
            Specific gas constant, J/(kg K); defaults to the registered dry-air
            value.
        point_count : int, default 61
            Number of path samples.

    Returns
    -------
    dict
        Fields of :class:`~rocket_propulsion.thermodynamics.PolytropicResult`.

    Raises
    ------
    InputError
        If a required field is missing or a numeric field is not a finite
        number.
    DomainError
        If the values are outside the model domain.

    Notes
    -----
    Quasistatic closed-system path ``p v^n = const`` for an ideal gas [1]_
    with constant heat capacities [2]_.

    References
    ----------
    .. [1] NASA-GAS. NASA Glenn Research Center, "Equation of State (Ideal Gas)",
           https://www1.grc.nasa.gov/beginners-guide-to-aeronautics/equation-of-state-ideal-gas-2/
    .. [2] NASA Glenn Research Center, "Specific Heats",
           https://www.grc.nasa.gov/www/k-12/airplane/specheat.html
    """

    return asdict(
        calculate_polytropic_process(
            inlet_pressure_pa=_number(payload, "inlet_pressure_pa"),
            outlet_pressure_pa=_number(payload, "outlet_pressure_pa"),
            inlet_temperature_k=_number(payload, "inlet_temperature_k"),
            exponent=_number(payload, "exponent"),
            gamma=_number(payload, "gamma", default=1.4),
            gas_constant_j_kg_k=_number(payload, "gas_constant_j_kg_k", default=AIR_SPECIFIC_GAS_CONSTANT_J_KG_K),
            point_count=_integer(payload, "point_count", default=61),
        )
    )


def calculate_duct_flow_route(payload: dict[str, Any]) -> dict[str, Any]:
    """Return Fanno or Rayleigh line properties for a local Mach number.

    Parameters
    ----------
    payload : dict
        Decoded JSON request body with the fields:

        model : {"fanno", "rayleigh"}, default "fanno"
            Constant-area duct model.
        mach : float
            Local Mach number.
        gamma : float, default 1.4
            Ratio of specific heats.

    Returns
    -------
    dict
        ``{"model": ..., **fields}`` with the sonic-reference ratios.

    Raises
    ------
    InputError
        If a required field is missing or a numeric field is not a finite
        number.
    DomainError
        If the values are outside the model domain.
        Also raised when ``model`` is not ``"fanno"`` or ``"rayleigh"``.

    Notes
    -----
    Constant-area adiabatic flow with friction (Fanno) or frictionless flow
    with heat addition (Rayleigh) [1]_.

    References
    ----------
    .. [1] NACA-TN-1336. B. L. Hicks, D. J. Montgomery and R. H. Wasserman, "On the
           One-Dimensional Theory of Steady Compressible Fluid Flow in Ducts with
           Friction and Heat Addition", NACA TN 1336, 1947,
           https://ntrs.nasa.gov/citations/19930081963
    """

    model = str(payload.get("model", "fanno")).strip().lower()
    mach = _number(payload, "mach")
    gamma = _number(payload, "gamma", default=1.4)
    if model == "fanno":
        return {"model": model, **asdict(calculate_fanno(mach, gamma))}
    if model == "rayleigh":
        return {"model": model, **asdict(calculate_rayleigh(mach, gamma))}
    raise DomainError("Duct-flow model must be fanno or rayleigh.")


def calculate_nozzle_route(payload: dict[str, Any]) -> dict[str, Any]:
    """Calculate ideal choked converging-diverging nozzle performance.

    Parameters
    ----------
    payload : dict
        Decoded JSON request body with the fields:

        chamber_pressure_pa, chamber_temperature_k : float
            Chamber stagnation state, Pa and K.
        throat_area_m2 : float
            Throat area, m^2.
        area_ratio : float
            Exit-to-throat area ratio, at least one.
        ambient_pressure_pa : float, default 0.0
            Ambient pressure, Pa.
        gamma : float, default 1.22
            Ratio of specific heats.
        gas_constant_j_kg_k : float, default 355.0
            Specific gas constant of the products, J/(kg K).

    Returns
    -------
    dict
        Fields of the ideal nozzle result (mass flow, exit state, thrust, Isp).

    Raises
    ------
    InputError
        If a required field is missing or a numeric field is not a finite
        number.
    DomainError
        If the values are outside the model domain.

    Notes
    -----
    Isentropic, inviscid, frozen-composition quasi-one-dimensional flow
    [1]_, [2]_; thrust is the momentum plus pressure thrust [3]_.

    References
    ----------
    .. [1] NACA-TR-1135. Ames Research Staff, "Equations, Tables, and Charts for
           Compressible Flow", NACA Report 1135, 1953,
           https://ntrs.nasa.gov/citations/19930091059
    .. [2] SUTTON-BIBLARZ-2017. G. P. Sutton and O. Biblarz, "Rocket Propulsion
           Elements", 9th ed., Wiley, 2017.
    .. [3] NASA-STEADY-ROCKET-THRUST. NASA Glenn Research Center, momentum and
           pressure thrust equation,
           https://www1.grc.nasa.gov/beginners-guide-to-aeronautics/rocket-thrust/
    """

    return asdict(
        calculate_ideal_nozzle(
            chamber_pressure_pa=_number(payload, "chamber_pressure_pa"),
            chamber_temperature_k=_number(payload, "chamber_temperature_k"),
            throat_area_m2=_number(payload, "throat_area_m2"),
            area_ratio=_number(payload, "area_ratio"),
            ambient_pressure_pa=_number(payload, "ambient_pressure_pa", default=0.0),
            gamma=_number(payload, "gamma", default=1.22),
            gas_constant_j_kg_k=_number(payload, "gas_constant_j_kg_k", default=355.0),
        )
    )


def calculate_off_design_nozzle_route(payload: dict[str, Any]) -> dict[str, Any]:
    """Solve the back-pressure-driven nozzle regime and optionally an altitude sweep.

    Parameters
    ----------
    payload : dict
        Decoded JSON request body with the fields:

        chamber_pressure_pa, chamber_temperature_k : float
            Chamber stagnation state, Pa and K.
        throat_area_m2 : float
            Throat area, m^2.
        area_ratio : float
            Exit-to-throat area ratio.
        ambient_pressure_pa : float
            Back pressure, Pa.
        gamma : float, default 1.22
            Ratio of specific heats.
        gas_constant_j_kg_k : float, default 355.0
            Specific gas constant, J/(kg K).
        station_count : int, default 81
            Number of axial stations.
        include_altitude_sweep : {"true", "false"}, default "true"
            Attach a U.S. Standard Atmosphere 1976 sweep.
        maximum_altitude_m : float, default 50000.0
            Upper altitude of the sweep, m.
        altitude_point_count : int, default 31
            Number of sweep points.

    Returns
    -------
    dict
        Fields of the off-design result, plus ``altitude_sweep`` when requested.

    Raises
    ------
    InputError
        If a required field is missing or a numeric field is not a finite
        number.
    DomainError
        If the values are outside the model domain.
        Also raised when ``include_altitude_sweep`` is not true or false.

    Notes
    -----
    Quasi-one-dimensional flow with a normal shock for internal compression
    [1]_, [2]_; separation is screened empirically [3]_ and ambient pressure
    follows the standard atmosphere [4]_.

    References
    ----------
    .. [1] SUTTON-BIBLARZ-2017. G. P. Sutton and O. Biblarz, "Rocket Propulsion
           Elements", 9th ed., Wiley, 2017.
    .. [2] NACA-TR-1135. Ames Research Staff, "Equations, Tables, and Charts for
           Compressible Flow", NACA Report 1135, 1953,
           https://ntrs.nasa.gov/citations/19930091059
    .. [3] NASA-SP8120-SEPARATION. NASA SP-8120, "Liquid Rocket Engine Nozzles",
           1976, p. 17, https://ntrs.nasa.gov/citations/19770009165
    .. [4] NOAA/NASA/USAF, "U.S. Standard Atmosphere, 1976", NOAA-S/T 76-1562 /
           NASA-TM-X-74335, https://ntrs.nasa.gov/citations/19770009539
    """

    inputs = {
        "chamber_pressure_pa": _number(payload, "chamber_pressure_pa"),
        "chamber_temperature_k": _number(payload, "chamber_temperature_k"),
        "throat_area_m2": _number(payload, "throat_area_m2"),
        "area_ratio": _number(payload, "area_ratio"),
        "ambient_pressure_pa": _number(payload, "ambient_pressure_pa"),
        "gamma": _number(payload, "gamma", default=1.22),
        "gas_constant_j_kg_k": _number(
            payload, "gas_constant_j_kg_k", default=355.0
        ),
        "station_count": _integer(payload, "station_count", default=81),
    }
    result = asdict(calculate_off_design_nozzle(**inputs))
    include_sweep = str(payload.get("include_altitude_sweep", "true")).lower()
    if include_sweep not in {"true", "false"}:
        raise InputError(
            "include_altitude_sweep must be true or false.",
            field="include_altitude_sweep",
        )
    if include_sweep == "true":
        result["altitude_sweep"] = asdict(
            generate_nozzle_altitude_sweep(
                chamber_pressure_pa=inputs["chamber_pressure_pa"],
                chamber_temperature_k=inputs["chamber_temperature_k"],
                throat_area_m2=inputs["throat_area_m2"],
                area_ratio=inputs["area_ratio"],
                gamma=inputs["gamma"],
                gas_constant_j_kg_k=inputs["gas_constant_j_kg_k"],
                maximum_altitude_m=_number(
                    payload, "maximum_altitude_m", default=50_000.0
                ),
                point_count=_integer(payload, "altitude_point_count", default=31),
            )
        )
    return result


def calculate_nozzle_contour_route(payload: dict[str, Any]) -> dict[str, Any]:
    """Generate a conical, preliminary, Rao or method-of-characteristics nozzle contour.

    Parameters
    ----------
    payload : dict
        Decoded JSON request body with the fields:

        contour : {"preliminary", "bell", "conical", "rao", "moc"}
            Contour family; ``"bell"`` is an alias of ``"preliminary"``.
        throat_area_m2 : float
            Throat area, m^2.
        area_ratio : float
            Exit-to-throat area ratio.
        gamma : float, default 1.22
            Ratio of specific heats.
        half_angle_deg, chamber_radius_ratio, bell_length_fraction : float
            Conical/preliminary shape parameters (defaults 15, 2.5 and 0.8).
        throat_angle_deg, exit_angle_deg : float
            Rao wall angles (defaults 30 and 8), deg.
        station_count : int, default 121
            Number of contour stations.
        characteristic_count : int, default 20
            Number of characteristics for the MOC design.
        axisymmetric_correction : {"true", "false"}, default "true"
            Apply the axisymmetric correction to the MOC design.
        include_exports, include_comparison : {"true", "false"}, default "false"
            Attach CSV/SVG/DXF/STL exports, or a comparison of all four families.
        angular_segments : int, default 32
            Revolution segments of the STL export.

    Returns
    -------
    dict
        Fields of the nozzle contour plus ``characteristic_lines``; MOC evidence,
        ``exports`` and ``comparison`` when applicable.

    Raises
    ------
    InputError
        If a required field is missing or a numeric field is not a finite
        number.
    DomainError
        If the values are outside the model domain.
        Also raised for an unknown contour family or a malformed boolean flag.

    Notes
    -----
    Rao thrust-optimized parabola [1]_; preliminary and conical shapes [2]_;
    MOC characteristic and compatibility equations [3]_ with averaged
    axisymmetric coefficients [4]_. Exports are a sampled fluid domain in
    metres (DXF [5]_, ASCII STL [6]_), not a manufactured wall.

    References
    ----------
    .. [1] RAO-1958. G. V. R. Rao, "Exhaust Nozzle Contour for Optimum Thrust",
           Jet Propulsion 28(6), 377-382, 1958, doi:10.2514/8.7324
    .. [2] SUTTON-BIBLARZ-2017. G. P. Sutton and O. Biblarz, "Rocket Propulsion
           Elements", 9th ed., Wiley, 2017.
    .. [3] NASA-MOC-CARTESIAN-1993. Aldo and Argrow, Cartesian characteristic and
           compatibility equations (1a,b), https://ntrs.nasa.gov/citations/19940019183
    .. [4] NASA-TRR33-AXISYMMETRIC. Guentert and Neumann, 1959, radius source and
           averaged characteristic coefficients,
           https://ntrs.nasa.gov/citations/19980232191
    .. [5] AUTODESK-LWPOLYLINE. Autodesk, lightweight-polyline subclass group codes,
           https://help.autodesk.com/cloudhelp/2018/ENU/AutoCAD-DXF/files/GUID-748FC305-F3F2-4F74-825A-61F04D757A50.htm
    .. [6] LOC-STL-ASCII. Library of Congress, ASCII STL records and metadata
           limitations, https://www.loc.gov/preservation/digital/formats/fdd/fdd000506.shtml
    """

    model = str(payload.get("contour", "preliminary")).strip().lower()
    if model == "bell":
        model = "preliminary"
    throat_area = _number(payload, "throat_area_m2")
    area_ratio = _number(payload, "area_ratio")
    gamma = _number(payload, "gamma", default=1.22)
    characteristic_lines: list[list[dict[str, Any]]] = []
    evidence: dict[str, Any] = {}
    contour: NozzleContour
    if model in {"conical", "preliminary"}:
        contour = generate_nozzle_contour(
            throat_area_m2=throat_area,
            area_ratio=area_ratio,
            gamma=gamma,
            contour=model,
            half_angle_deg=_number(payload, "half_angle_deg", default=15.0),
            chamber_radius_ratio=_number(payload, "chamber_radius_ratio", default=2.5),
            bell_length_fraction=_number(payload, "bell_length_fraction", default=0.8),
            station_count=_integer(payload, "station_count", default=121),
        )
    elif model == "rao":
        contour = generate_rao_contour(
            throat_area_m2=throat_area,
            area_ratio=area_ratio,
            gamma=gamma,
            length_fraction=_number(payload, "bell_length_fraction", default=0.8),
            throat_angle_deg=_number(payload, "throat_angle_deg", default=30.0),
            exit_angle_deg=_number(payload, "exit_angle_deg", default=8.0),
            chamber_radius_ratio=_number(payload, "chamber_radius_ratio", default=2.5),
            station_count=_integer(payload, "station_count", default=121),
        )
    elif model == "moc":
        correction = str(payload.get("axisymmetric_correction", "true")).lower()
        if correction not in {"true", "false"}:
            raise InputError(
                "axisymmetric_correction must be true or false.",
                field="axisymmetric_correction",
            )
        design = generate_moc_nozzle(
            throat_area_m2=throat_area,
            area_ratio=area_ratio,
            gamma=gamma,
            characteristic_count=_integer(payload, "characteristic_count", default=20),
            axisymmetric_correction=correction == "true",
        )
        contour = design.contour
        characteristic_lines = [
            [asdict(point) for point in line] for line in design.characteristic_lines
        ]
        evidence = {
            "geometry_mode": design.geometry_mode,
            "exit_mach_target": design.exit_mach_target,
            "exit_mach_computed": design.exit_mach_computed,
            "characteristic_count": design.characteristic_count,
            "residuals": asdict(design.residuals),
            "refined_length_change_relative": design.refined_length_change_relative,
            "source_url": design.source_url,
        }
    else:
        raise DomainError(
            "Nozzle contour must be conical, preliminary, rao, or moc.",
            field="contour",
        )

    result: dict[str, Any] = asdict(contour)
    result["characteristic_lines"] = characteristic_lines
    result.update(evidence)
    include_exports = str(payload.get("include_exports", "false")).lower()
    if include_exports not in {"true", "false"}:
        raise InputError("include_exports must be true or false.", field="include_exports")
    if include_exports == "true":
        result["exports"] = asdict(
            export_nozzle_contour(
                contour,
                angular_segments=_integer(payload, "angular_segments", default=32),
            )
        )
    include_comparison = str(payload.get("include_comparison", "false")).lower()
    if include_comparison not in {"true", "false"}:
        raise InputError(
            "include_comparison must be true or false.", field="include_comparison"
        )
    if include_comparison == "true":
        station_count = _integer(payload, "station_count", default=121)
        chamber_ratio = _number(payload, "chamber_radius_ratio", default=2.5)
        length_fraction = _number(payload, "bell_length_fraction", default=0.8)
        contours = (
            generate_nozzle_contour(
                throat_area_m2=throat_area,
                area_ratio=area_ratio,
                gamma=gamma,
                contour="conical",
                half_angle_deg=_number(payload, "half_angle_deg", default=15.0),
                chamber_radius_ratio=chamber_ratio,
                station_count=station_count,
            ),
            generate_nozzle_contour(
                throat_area_m2=throat_area,
                area_ratio=area_ratio,
                gamma=gamma,
                contour="preliminary",
                chamber_radius_ratio=chamber_ratio,
                bell_length_fraction=length_fraction,
                station_count=station_count,
            ),
            generate_rao_contour(
                throat_area_m2=throat_area,
                area_ratio=area_ratio,
                gamma=gamma,
                length_fraction=length_fraction,
                throat_angle_deg=_number(payload, "throat_angle_deg", default=30.0),
                exit_angle_deg=_number(payload, "exit_angle_deg", default=8.0),
                chamber_radius_ratio=chamber_ratio,
                station_count=station_count,
            ),
            generate_moc_nozzle(
                throat_area_m2=throat_area,
                area_ratio=area_ratio,
                gamma=gamma,
                characteristic_count=_integer(
                    payload, "characteristic_count", default=20
                ),
                axisymmetric_correction=True,
            ).contour,
        )
        result["comparison"] = [
            asdict(record) for record in compare_nozzle_contours(contours)
        ]
    return result


def calculate_rocket_equation_route(payload: dict[str, Any]) -> dict[str, Any]:
    """Solve the Tsiolkovsky rocket equation for its missing quantity.

    Parameters
    ----------
    payload : dict
        Decoded JSON request body with the fields:

        initial_mass_kg : float
            Initial vehicle mass, kg.
        final_mass_kg : float, optional
            Final vehicle mass, kg; give it or ``delta_v_m_s``.
        delta_v_m_s : float, optional
            Velocity increment, m/s.
        specific_impulse_s, effective_exhaust_velocity_m_s : float, optional
            Performance input, s or m/s; exactly one is required.

    Returns
    -------
    dict
        Fields of the rocket-equation result.

    Raises
    ------
    InputError
        If a required field is missing or a numeric field is not a finite
        number.
    DomainError
        If the values are outside the model domain.

    Notes
    -----
    ``dv = c ln(m0/mf)`` with ``c = Isp g_n`` [1]_, [2]_; the domain function
    decides which combination of inputs is accepted.

    References
    ----------
    .. [1] NASA-ROCKET. NASA Glenn Research Center, "Ideal Rocket Equation",
           https://www1.grc.nasa.gov/beginners-guide-to-aeronautics/ideal-rocket-equation/
    .. [2] NASA-ISP. NASA Glenn Research Center, "Specific Impulse",
           https://www1.grc.nasa.gov/beginners-guide-to-aeronautics/specific-impulse/
    """

    return asdict(
        calculate_rocket_equation(
            initial_mass_kg=_optional_number(payload, "initial_mass_kg"),
            final_mass_kg=_optional_number(payload, "final_mass_kg"),
            delta_v_m_s=_optional_number(payload, "delta_v_m_s"),
            specific_impulse_s=_optional_number(payload, "specific_impulse_s"),
            effective_exhaust_velocity_m_s=_optional_number(
                payload, "effective_exhaust_velocity_m_s"
            ),
        )
    )


def calculate_thrust_route(payload: dict[str, Any]) -> dict[str, Any]:
    """Calculate thrust, equivalent exhaust velocity, Isp and optional total impulse.

    Parameters
    ----------
    payload : dict
        Decoded JSON request body with the fields:

        mass_flow_kg_s : float
            Propellant mass flow, kg/s.
        exit_velocity_m_s : float
            Exit-plane velocity, m/s.
        exit_area_m2 : float
            Exit area, m^2.
        exit_pressure_pa : float
            Exit-plane static pressure, Pa.
        ambient_pressure_pa : float, default 0.0
            Ambient pressure, Pa.
        burn_time_s : float, optional
            Burn duration for total impulse, s.

    Returns
    -------
    dict
        Fields of the thrust result.

    Raises
    ------
    InputError
        If a required field is missing or a numeric field is not a finite
        number.
    DomainError
        If the values are outside the model domain.

    Notes
    -----
    ``F = mdot Ve + (pe - pa) Ae`` [1]_ and ``Isp = F / (mdot g_n)`` [2]_.

    References
    ----------
    .. [1] NASA-STEADY-ROCKET-THRUST. NASA Glenn Research Center, momentum and
           pressure thrust equation,
           https://www1.grc.nasa.gov/beginners-guide-to-aeronautics/rocket-thrust/
    .. [2] NASA-ISP. NASA Glenn Research Center, "Specific Impulse",
           https://www1.grc.nasa.gov/beginners-guide-to-aeronautics/specific-impulse/
    """

    return asdict(
        calculate_thrust(
            mass_flow_kg_s=_number(payload, "mass_flow_kg_s"),
            exit_velocity_m_s=_number(payload, "exit_velocity_m_s"),
            exit_area_m2=_number(payload, "exit_area_m2"),
            exit_pressure_pa=_number(payload, "exit_pressure_pa"),
            ambient_pressure_pa=_number(payload, "ambient_pressure_pa", default=0.0),
            burn_time_s=_optional_number(payload, "burn_time_s"),
        )
    )


def _object(value: Any, field: str) -> dict[str, Any]:
    """Return ``value`` if it is a JSON object, else raise :class:`InputError`."""
    if not isinstance(value, dict):
        raise InputError(f"Field {field} must be an object.", field=field)
    return value


def _objects(value: Any, field: str) -> list[dict[str, Any]]:
    """Return ``value`` if it is a non-empty JSON array of objects."""
    if not isinstance(value, list) or not value:
        raise InputError(f"Field {field} must be a non-empty list.", field=field)
    if not all(isinstance(item, dict) for item in value):
        raise InputError(f"Every item in {field} must be an object.", field=field)
    return value


def _enum_value(enum_type: Any, value: Any, field: str) -> Any:
    """Convert ``value`` to ``enum_type``; the error lists the allowed values."""
    try:
        return enum_type(str(value))
    except ValueError as error:
        choices = ", ".join(item.value for item in enum_type)
        raise InputError(
            f"Field {field} must be one of: {choices}.", field=field
        ) from error


def simulate_burn_route(payload: dict[str, Any]) -> dict[str, Any]:
    """Evaluate a standalone constant (L0) or throttle-profiled (L1) propulsion burn.

    Parameters
    ----------
    payload : dict
        Decoded JSON request body with the fields:

        definition : object
            Burn definition: ``name``, ``initial_mass_kg``,
            ``protected_dry_mass_kg``, ``tanks`` (non-empty list), ``target``
            (``kind`` and ``value``), ``maximum_duration_s`` and
            ``cant_efficiency`` (default 1).
        operating_point : object
            Rated operating point: throttles, engine count, ideal and delivered
            thrust, chamber and system Isp, mass flows, mixture ratio, pressures
            and a non-empty ``streams`` list.
        schedule : object, optional
            Piecewise-linear throttle schedule with non-empty ``segments``; when
            absent a constant burn is simulated.

    Returns
    -------
    dict
        Serialized, hash-bearing burn result.

    Raises
    ------
    InputError
        If a required field is missing or a numeric field is not a finite
        number.
    DomainError
        If the values are outside the model domain.
        Also raised for a malformed object, list or enumeration value.

    Notes
    -----
    Rated delivered thrust and explicit tank flows; a throttle schedule scales
    thrust and all streams. The velocity increment is the ideal one-dimensional
    value [1]_ with constant system specific impulse [2]_; gravity, drag,
    steering and orbital propagation are not modelled.

    References
    ----------
    .. [1] NASA-ROCKET. NASA Glenn Research Center, "Ideal Rocket Equation",
           https://www1.grc.nasa.gov/beginners-guide-to-aeronautics/ideal-rocket-equation/
    .. [2] NASA-ISP. NASA Glenn Research Center, "Specific Impulse",
           https://www1.grc.nasa.gov/beginners-guide-to-aeronautics/specific-impulse/
    """

    definition_data = _object(payload.get("definition"), "definition")
    target_data = _object(definition_data.get("target"), "definition.target")
    tank_items = _objects(definition_data.get("tanks"), "definition.tanks")
    definition = BurnDefinition(
        name=str(definition_data.get("name", "constant burn")),
        initial_mass_kg=_number(definition_data, "initial_mass_kg"),
        protected_dry_mass_kg=_number(definition_data, "protected_dry_mass_kg"),
        tanks=tuple(
            TankDefinition(
                tank_id=str(item.get("tank_id", "")),
                propellant_key=str(item.get("propellant_key", "")),
                role=_enum_value(
                    PropellantRole, item.get("role"), f"definition.tanks[{index}].role"
                ),
                loaded_mass_kg=_number(item, "loaded_mass_kg"),
                reserve_mass_kg=_number(item, "reserve_mass_kg", default=0.0),
            )
            for index, item in enumerate(tank_items)
        ),
        target=BurnTarget(
            kind=_enum_value(
                BurnTargetKind, target_data.get("kind"), "definition.target.kind"
            ),
            value=_number(target_data, "value"),
        ),
        maximum_duration_s=_number(definition_data, "maximum_duration_s"),
        cant_efficiency=_number(definition_data, "cant_efficiency", default=1.0),
    )

    point_data = _object(payload.get("operating_point"), "operating_point")
    stream_items = _objects(point_data.get("streams"), "operating_point.streams")
    point = EngineOperatingPoint(
        name=str(point_data.get("name", "constant operating point")),
        commanded_throttle=_number(point_data, "commanded_throttle", default=1.0),
        realized_throttle=_number(point_data, "realized_throttle", default=1.0),
        active_engine_count=_integer(point_data, "active_engine_count", default=1),
        ideal_thrust_n=_number(point_data, "ideal_thrust_n"),
        delivered_thrust_n=_number(point_data, "delivered_thrust_n"),
        chamber_specific_impulse_s=_number(point_data, "chamber_specific_impulse_s"),
        system_specific_impulse_s=_number(point_data, "system_specific_impulse_s"),
        chamber_mass_flow_kg_s=_number(point_data, "chamber_mass_flow_kg_s"),
        total_tank_flow_kg_s=_number(point_data, "total_tank_flow_kg_s"),
        oxidizer_mass_flow_kg_s=_number(
            point_data, "oxidizer_mass_flow_kg_s", default=0.0
        ),
        fuel_mass_flow_kg_s=_number(point_data, "fuel_mass_flow_kg_s", default=0.0),
        mixture_ratio=_optional_number(point_data, "mixture_ratio"),
        chamber_pressure_pa=_number(point_data, "chamber_pressure_pa"),
        ambient_pressure_pa=_number(point_data, "ambient_pressure_pa", default=0.0),
        streams=tuple(
            PropellantFlow(
                stream_id=str(item.get("stream_id", "")),
                tank_id=str(item.get("tank_id", "")),
                propellant_key=str(item.get("propellant_key", "")),
                role=_enum_value(
                    PropellantRole,
                    item.get("role"),
                    f"operating_point.streams[{index}].role",
                ),
                destination=_enum_value(
                    FlowDestination,
                    item.get("destination"),
                    f"operating_point.streams[{index}].destination",
                ),
                mass_flow_kg_s=_number(item, "mass_flow_kg_s"),
                tank_depleting=_boolean(item, "tank_depleting", default=True),
                contributes_to_delivered_thrust=_boolean(
                    item, "contributes_to_delivered_thrust", default=True
                ),
            )
            for index, item in enumerate(stream_items)
        ),
        model_id=str(point_data.get("model_id", "constant_performance_v1")),
        source=str(point_data.get("source", "api")),
        warnings=tuple(str(item) for item in point_data.get("warnings", ())),
    )
    schedule_payload = payload.get("schedule")
    if schedule_payload is None:
        return burn_result_to_dict(simulate_constant_burn(definition, point))
    schedule_data = _object(schedule_payload, "schedule")
    segment_items = _objects(schedule_data.get("segments"), "schedule.segments")
    schedule = ThrottleSchedule(
        name=str(schedule_data.get("name", "piecewise-linear throttle schedule")),
        segments=tuple(
            ThrottleSegment(
                duration_s=_number(item, "duration_s"),
                start_throttle=_number(item, "start_throttle"),
                end_throttle=_number(item, "end_throttle"),
                phase=str(item.get("phase", f"segment-{index + 1}")),
            )
            for index, item in enumerate(segment_items)
        ),
    )
    provider = ProfiledPerformanceProvider(point, schedule)
    return burn_result_to_dict(
        simulate_profiled_burn(definition, provider.operating_point(), provider.schedule)
    )


def calculate_engine_operating_point_route(payload: dict[str, Any]) -> dict[str, Any]:
    """Resolve a constant-performance shortcut into an explicit operating point.

    Parameters
    ----------
    payload : dict
        Decoded JSON request body with the fields:

        family : {"bipropellant", "monopropellant"}, default "bipropellant"
            Provider family.
        delivered_thrust_n, system_specific_impulse_s, chamber_pressure_pa : float
            Delivered performance and chamber pressure.
        ideal_thrust_n, chamber_specific_impulse_s : float, optional
            Optional ideal thrust and chamber Isp.
        commanded_throttle, realized_throttle : float, default 1.0
            Throttle settings.
        active_engine_count : int, default 1
            Number of firing engines.
        ambient_pressure_pa : float, default 0.0
            Ambient pressure, Pa.
        oxidizer_fuel_ratio : float
            Bipropellant mixture ratio (bipropellant only).
        dump_mass_flow_kg_s, dump_role, tank and propellant keys : optional
            Further stream bookkeeping fields.

    Returns
    -------
    dict
        Fields of the resolved operating point, including explicit streams.

    Raises
    ------
    InputError
        If a required field is missing or a numeric field is not a finite
        number.
    DomainError
        If the values are outside the model domain.
        Also raised for an unknown family.

    Notes
    -----
    Total tank flow follows ``mdot = F / (Isp g_n)`` [1]_; a bipropellant
    split uses a constant oxidizer-to-fuel ratio [2]_.

    References
    ----------
    .. [1] NASA-ISP. NASA Glenn Research Center, "Specific Impulse",
           https://www1.grc.nasa.gov/beginners-guide-to-aeronautics/specific-impulse/
    .. [2] SUTTON-BIBLARZ-2017. G. P. Sutton and O. Biblarz, "Rocket Propulsion
           Elements", 9th ed., Wiley, 2017.
    """

    family = str(payload.get("family", "bipropellant")).strip().lower()
    common = {
        "delivered_thrust_n": _number(payload, "delivered_thrust_n"),
        "system_specific_impulse_s": _number(payload, "system_specific_impulse_s"),
        "ideal_thrust_n": _optional_number(payload, "ideal_thrust_n"),
        "commanded_throttle": _number(payload, "commanded_throttle", default=1.0),
        "realized_throttle": _number(payload, "realized_throttle", default=1.0),
        "active_engine_count": _integer(payload, "active_engine_count", default=1),
        "chamber_pressure_pa": _number(payload, "chamber_pressure_pa"),
        "ambient_pressure_pa": _number(payload, "ambient_pressure_pa", default=0.0),
        "name": str(payload.get("name", "constant operating point")),
        "source": str(payload.get("source", "constant provider API")),
    }
    if family == "bipropellant":
        provider = ConstantPerformanceProvider.bipropellant(
            **common,
            chamber_specific_impulse_s=_optional_number(
                payload, "chamber_specific_impulse_s"
            ),
            oxidizer_fuel_ratio=_number(payload, "oxidizer_fuel_ratio"),
            dump_mass_flow_kg_s=_number(payload, "dump_mass_flow_kg_s", default=0.0),
            dump_role=_enum_value(
                PropellantRole, payload.get("dump_role", "fuel"), "dump_role"
            ),
            oxidizer_tank_id=str(payload.get("oxidizer_tank_id", "oxidizer")),
            fuel_tank_id=str(payload.get("fuel_tank_id", "fuel")),
            oxidizer_key=str(payload.get("oxidizer_key", "oxidizer")),
            fuel_key=str(payload.get("fuel_key", "fuel")),
        )
    elif family == "monopropellant":
        provider = ConstantPerformanceProvider.monopropellant(
            **common,
            tank_id=str(payload.get("tank_id", "propellant")),
            propellant_key=str(payload.get("propellant_key", "monopropellant")),
        )
    else:
        raise InputError(
            "Field family must be bipropellant or monopropellant.", field="family"
        )
    return asdict(provider.operating_point())


def build_burn_report_route(payload: dict[str, Any]) -> dict[str, Any]:
    """Render a verified burn artifact as a standalone HTML report.

    Parameters
    ----------
    payload : dict
        Request body with ``burn_result``, a serialized burn result object.

    Returns
    -------
    dict
        ``{"html": ..., "result_hash": ...}``.

    Raises
    ------
    InputError
        If ``burn_result`` is not an object or fails its integrity checks.
    """

    result = burn_result_from_dict(_object(payload.get("burn_result"), "burn_result"))
    return {"html": build_burn_html_report(result), "result_hash": result.result_hash}


def export_sidera_burn_route(payload: dict[str, Any]) -> dict[str, Any]:
    """Transform a verified burn result into a Sidera finite-burn hand-off document.

    Parameters
    ----------
    payload : dict
        Decoded JSON request body with the fields:

        burn_result : object
            Serialized burn result; re-verified before use.
        direction : list of 3 numbers
            Thrust direction in ``frame``.
        t_start_s : float
            Burn start time, s.
        frame : {"inertial", "ric", "vnb"}, default "inertial"
            Direction frame.
        mode : {"exact", "equivalent_constant"}, default "exact"
            Export mode (see
            :func:`~rocket_propulsion.integrations.export_sidera_artifact`).

    Returns
    -------
    dict
        The ``sidera_maneuver_plan_v1`` document.

    Raises
    ------
    InputError
        If a required field is missing or a numeric field is not a finite
        number.
    DomainError
        If the values are outside the model domain.
        Also raised for a malformed direction, frame or mode.

    Notes
    -----
    Thrust is the net axial thrust and ``isp_s`` represents total tank drain
    [1]_; the document is versioned and unit-explicit [2]_.

    References
    ----------
    .. [1] NASA-ISP. NASA Glenn Research Center, "Specific Impulse",
           https://www1.grc.nasa.gov/beginners-guide-to-aeronautics/specific-impulse/
    .. [2] CCSDS 502.0-B-3. CCSDS, "Orbit Data Messages", Recommended Standard
           502.0-B-3, https://ccsds.org/Pubs/502x0b3e1.pdf
    """

    result_payload = _object(payload.get("burn_result"), "burn_result")
    raw_direction = payload.get("direction")
    if not isinstance(raw_direction, (list, tuple)) or len(raw_direction) != 3:
        raise InputError("Field direction must contain three numbers.", field="direction")
    direction_data = {f"component_{index}": value for index, value in enumerate(raw_direction)}
    direction = tuple(
        _number(direction_data, f"component_{index}") for index in range(3)
    )
    return export_sidera_artifact(
        burn_result_from_dict(result_payload),
        t_start_s=_number(payload, "t_start_s"),
        direction=direction,  # type: ignore[arg-type]
        frame=str(payload.get("frame", "inertial")),  # type: ignore[arg-type]
        mode=str(payload.get("mode", "exact")),  # type: ignore[arg-type]
    )


_SIDERA_CAPABILITY_PROBE = SideraCapabilityProbe(lambda: inspect_sidera_capabilities())


def sidera_capabilities_route(payload: dict[str, Any]) -> dict[str, Any]:
    """Return the Sidera capability snapshot without waiting for optional imports.

    Parameters
    ----------
    payload : dict
        Ignored; accepted for a uniform route signature.

    Returns
    -------
    dict
        The completed capability record with ``probe_complete=True``, or a
        ``"checking"`` record with ``probe_complete=False`` while the
        background inspection runs.

    References
    ----------
    .. [1] Python 3.12, "threading - Thread-based parallelism",
           https://docs.python.org/3.12/library/threading.html
    """

    del payload
    return _SIDERA_CAPABILITY_PROBE.snapshot()


def calculate_loss_budget_route(payload: dict[str, Any]) -> dict[str, Any]:
    """Evaluate a declared-basis performance-loss budget with first-order covariance.

    Parameters
    ----------
    payload : dict
        Decoded JSON request body with the fields:

        ideal_thrust_n, ideal_specific_impulse_s : float
            Ideal reference performance.
        discharge_coefficient, divergence_efficiency : float
            Loss factors (defaults 0.98 and 0.98).
        boundary_layer_efficiency, combustion_efficiency : float
            Loss factors (defaults 0.985 and 0.97).
        enable_discharge_loss, enable_divergence_loss : bool, default true
            Switches for the individual factors.
        enable_boundary_layer_loss, enable_combustion_loss : bool, default true
            Switches for the individual factors.
        relative_uncertainty : float, default 0.02
            Relative standard uncertainty of each factor.
        operating_basis : str, default "fixed-mass-flow-lumped"
            Declared operating basis of the factor ratios.
        correlation_matrix : 4x4 nested list, default identity
            Factor correlation in the order combustion, discharge, boundary layer,
            divergence.

    Returns
    -------
    dict
        Fields of the loss-budget result, including the retained basis and
        covariance terms.

    Raises
    ------
    InputError
        If a required field is missing or a numeric field is not a finite
        number.
    DomainError
        If the values are outside the model domain.

    Notes
    -----
    Discharge-coefficient and nozzle-loss definitions [1]_, [2]_; the
    intrinsic c*/discharge convention is that of [3]_. Uncertainty is
    propagated to first order with signed covariances [4]_. Sequential
    attribution does not prove that the physical losses do not overlap.

    References
    ----------
    .. [1] SUTTON-BIBLARZ-2017. G. P. Sutton and O. Biblarz, "Rocket Propulsion
           Elements", 9th ed., Wiley, 2017.
    .. [2] NASA-SP8120-DISCHARGE. NASA SP-8120, "Liquid Rocket Engine Nozzles", 1976,
           discharge coefficient and flow effects, printed pp. 10, 19-20 and 93,
           https://ntrs.nasa.gov/api/citations/19770009165/downloads/19770009165.pdf
    .. [3] ROCKETISP-CSTAR-CD-CONVENTION. RocketIsp 0.1.13 documentation,
           intrinsic c*/discharge convention (not independent NASA qualification),
           https://rocketisp.readthedocs.io/en/latest/efficiencies.html
    .. [4] NIST-TN1297-COVARIANCE. NIST TN 1297, Appendix A, law of propagation of
           uncertainty with sensitivities and covariances,
           https://www.nist.gov/pml/nist-technical-note-1297/nist-tn-1297-appendix-law-propagation-uncertainty
    """

    return asdict(
        calculate_loss_budget(
            ideal_thrust_n=_number(payload, "ideal_thrust_n"),
            ideal_specific_impulse_s=_number(payload, "ideal_specific_impulse_s"),
            discharge_coefficient=_number(
                payload, "discharge_coefficient", default=0.98
            ),
            divergence_efficiency=_number(
                payload, "divergence_efficiency", default=0.98
            ),
            boundary_layer_efficiency=_number(
                payload, "boundary_layer_efficiency", default=0.985
            ),
            combustion_efficiency=_number(
                payload, "combustion_efficiency", default=0.97
            ),
            enable_discharge_loss=_boolean(
                payload, "enable_discharge_loss", default=True
            ),
            enable_divergence_loss=_boolean(
                payload, "enable_divergence_loss", default=True
            ),
            enable_boundary_layer_loss=_boolean(
                payload, "enable_boundary_layer_loss", default=True
            ),
            enable_combustion_loss=_boolean(
                payload, "enable_combustion_loss", default=True
            ),
            relative_uncertainty=_number(
                payload, "relative_uncertainty", default=0.02
            ),
            operating_basis=payload.get("operating_basis", "fixed-mass-flow-lumped"),
            correlation_matrix=payload.get("correlation_matrix", (
                (1., 0., 0., 0.), (0., 1., 0., 0.), (0., 0., 1., 0.), (0., 0., 0., 1.))),
        )
    )


def _thermal_contour(payload: dict[str, Any]) -> NozzleContour:
    """Build the nozzle contour for the thermal route (Rao by default)."""
    model = str(payload.get("contour", "rao")).strip().lower()
    common = {
        "throat_area_m2": _number(payload, "throat_area_m2"),
        "area_ratio": _number(payload, "area_ratio"),
        "gamma": _number(payload, "gamma", default=1.22),
    }
    if model == "rao":
        return generate_rao_contour(
            **common,
            length_fraction=_number(payload, "bell_length_fraction", default=0.8),
            station_count=_integer(payload, "station_count", default=121),
        )
    if model == "moc":
        return generate_moc_nozzle(
            **common,
            characteristic_count=_integer(payload, "characteristic_count", default=20),
        ).contour
    return generate_nozzle_contour(
        **common,
        contour=model,
        half_angle_deg=_number(payload, "half_angle_deg", default=15.0),
        bell_length_fraction=_number(payload, "bell_length_fraction", default=0.8),
        station_count=_integer(payload, "station_count", default=121),
    )


def calculate_nozzle_thermal_route(payload: dict[str, Any]) -> dict[str, Any]:
    """Calculate a prescribed Bartz wall heat-flux profile and an optional coolant balance.

    Parameters
    ----------
    payload : dict
        Decoded JSON request body with the fields:

        contour and its shape fields : optional
            As for the nozzle-contour route; Rao by default.
        chamber_pressure_pa, chamber_temperature_k, wall_temperature_k : float
            Chamber stagnation state and uniform wall temperature.
        characteristic_velocity_m_s, dynamic_viscosity_pa_s : float
            ``c*`` and gas viscosity.
        specific_heat_j_kg_k, prandtl_number : float
            Constant gas properties.
        gamma : float, default 1.22
            Ratio of specific heats.
        throat_radius_of_curvature_m : float
            Throat wall radius of curvature, m.
        include_cooling : bool, default false
            Add the bulk coolant balance driven by the total heat load.
        coolant_* and channel fields : float
            Coolant properties and channel geometry (``channel_count`` default 100,
            ``laminar_darcy_poiseuille_number`` default 64).

    Returns
    -------
    dict
        Fields of the Bartz result, plus ``cooling`` when requested.

    Raises
    ------
    InputError
        If a required field is missing or a numeric field is not a finite
        number.
    DomainError
        If the values are outside the model domain.
        Also raised with code ``"unsupported_nonpositive_bulk_heat"`` when
        cooling is requested for a nonpositive net heat load.

    Notes
    -----
    Bartz coefficient [1]_ with the recovery convention of [2]_; the coolant
    balance uses a supplied laminar ``f Re`` [3]_ and the Swamee-Jain formula
    [4]_. Gas properties and wall temperature are prescribed; no coupled
    wall/coolant solve and no physical qualification.

    References
    ----------
    .. [1] NASA-SP125-BARTZ-PAGES. NASA SP-125, "Design of Liquid Propellant Rocket
           Engines", printed pp. 100-101, equations 4-10, 4-13 and 4-14,
           https://ntrs.nasa.gov/api/citations/19710019929/downloads/19710019929.pdf
    .. [2] NASA-HOPPER-2012-THERMAL. Linne et al., AIAA 2012-0806, printed p. 8,
           equations 1, 2 and 4,
           https://ntrs.nasa.gov/api/citations/20120004103/downloads/20120004103.pdf
    .. [3] NASA-GFSSP6-NONCIRCULAR. NASA GFSSP v6 user manual, printed p. 40
           (PDF p. 80), equation 28,
           https://www.nasa.gov/wp-content/uploads/2024/04/gfssp-v6-usermanual.pdf
    .. [4] SWAMEE-JAIN-1976-IDENTITY. Swamee and Jain, "Explicit Equations for
           Pipe-Flow Problems", J. Hydraulics Division 102(5), 657-664, 1976,
           https://doi.org/10.1061/JYCEAJ.0004542
    """

    heat = calculate_bartz_heat_transfer(
        contour=_thermal_contour(payload),
        chamber_pressure_pa=_number(payload, "chamber_pressure_pa"),
        chamber_temperature_k=_number(payload, "chamber_temperature_k"),
        wall_temperature_k=_number(payload, "wall_temperature_k"),
        characteristic_velocity_m_s=_number(
            payload, "characteristic_velocity_m_s"
        ),
        dynamic_viscosity_pa_s=_number(payload, "dynamic_viscosity_pa_s"),
        specific_heat_j_kg_k=_number(payload, "specific_heat_j_kg_k"),
        prandtl_number=_number(payload, "prandtl_number"),
        gamma=_number(payload, "gamma", default=1.22),
        throat_radius_of_curvature_m=_number(
            payload, "throat_radius_of_curvature_m"
        ),
    )
    result = asdict(heat)
    if _boolean(payload, "include_cooling", default=False):
        if heat.total_heat_load_w <= 0:
            raise DomainError("The current bulk coolant model requires positive net absorbed heat; signed gas-side heat is not rectified.",
                              code="unsupported_nonpositive_bulk_heat", field="include_cooling")
        wall_material = _catalog_key(payload, "wall_material")
        result["cooling"] = asdict(
            calculate_regenerative_cooling(
                absorbed_heat_w=heat.total_heat_load_w,
                coolant_mass_flow_kg_s=_number(payload, "coolant_mass_flow_kg_s"),
                coolant_specific_heat_j_kg_k=_optional_number(
                    payload, "coolant_specific_heat_j_kg_k"
                ),
                coolant_inlet_temperature_k=_number(
                    payload, "coolant_inlet_temperature_k"
                ),
                coolant_density_kg_m3=_optional_number(payload, "coolant_density_kg_m3"),
                coolant_dynamic_viscosity_pa_s=_optional_number(
                    payload, "coolant_dynamic_viscosity_pa_s"
                ),
                channel_count=_integer(payload, "channel_count", default=100),
                channel_flow_area_m2=_number(payload, "channel_flow_area_m2"),
                hydraulic_diameter_m=_number(payload, "hydraulic_diameter_m"),
                channel_length_m=_number(payload, "channel_length_m"),
                roughness_m=_number(payload, "roughness_m"),
                maximum_coolant_temperature_k=_number(
                    payload, "maximum_coolant_temperature_k"
                ),
                laminar_darcy_poiseuille_number=_number(
                    payload, "laminar_darcy_poiseuille_number", default=64.0
                ),
                coolant=_catalog_key(payload, "coolant"),
                wall_material=wall_material,
                wall_thickness_m=_optional_number(payload, "wall_thickness_m") if wall_material else None,
                wall_heat_flux_w_m2=heat.peak_heat_flux_w_m2 if wall_material else None,
                coolant_side_heat_transfer_coefficient_w_m2_k=(
                    _optional_number(payload, "coolant_side_heat_transfer_coefficient_w_m2_k")
                    if wall_material else None
                ),
            )
        )
    return result


def _catalog_key(payload: dict[str, Any], key: str) -> str | None:
    """Read an optional catalog key; absent, ``None`` or blank text means none."""

    raw = payload.get(key)
    if raw in (None, ""):
        return None
    if not isinstance(raw, str):
        raise InputError(f"Field {key} must be a catalog key.", field=key)
    return raw.strip()


def calculate_two_phase_band_route(payload: dict[str, Any]) -> dict[str, Any]:
    """Return conditional gas/particle exit-state scenarios with an explicit coupling definition.

    Parameters
    ----------
    payload : dict
        Decoded JSON request body with the fields:

        ideal_thrust_n, ideal_specific_impulse_s : float
            Ideal reference performance.
        condensed_mass_fraction : float
            Condensed-phase mass fraction.
        minimum_particle_coupling_efficiency : float
            Lower scenario endpoint of the coupling efficiency.
        maximum_particle_coupling_efficiency : float
            Upper scenario endpoint of the coupling efficiency.
        reference_pressure_thrust_n : float, default 0.0
            Separate pressure-thrust term, N.
        coupling_definition : str, default "exit-velocity-ratio"
            Meaning of the coupling efficiency.

    Returns
    -------
    dict
        Fields of the two-phase band result.

    Raises
    ------
    InputError
        If a required field is missing or a numeric field is not a finite
        number.
    DomainError
        If the values are outside the model domain.

    Notes
    -----
    Separate gas and particle momentum with pressure thrust [1]_ and the
    momentum debit of coupled gas/particle flow [2]_. Scenario endpoints are
    not validated physical bounds or statistical coverage.

    References
    ----------
    .. [1] DREWNIAK-TWO-PHASE-2022. Drewniak, Boyer and Yetter, 2022, equation (3),
           separate gas/particle momentum and pressure,
           https://ntrs.nasa.gov/citations/20220006969
    .. [2] KIBBEY-TWO-PHASE-2012. Kibbey, 2012, coupled gas/particle flow and
           momentum debit, https://ntrs.nasa.gov/citations/20120015678
    """

    return asdict(
        calculate_two_phase_band(
            ideal_thrust_n=_number(payload, "ideal_thrust_n"),
            ideal_specific_impulse_s=_number(payload, "ideal_specific_impulse_s"),
            condensed_mass_fraction=_number(payload, "condensed_mass_fraction"),
            minimum_particle_coupling_efficiency=_number(
                payload, "minimum_particle_coupling_efficiency"
            ),
            maximum_particle_coupling_efficiency=_number(
                payload, "maximum_particle_coupling_efficiency"
            ),
            reference_pressure_thrust_n=_number(payload, "reference_pressure_thrust_n", default=0.),
            coupling_definition=payload.get("coupling_definition", "exit-velocity-ratio"),
        )
    )


def calculate_staging_route(payload: dict[str, Any]) -> dict[str, Any]:
    """Optimize multi-stage delta-v allocation for minimum initial mass.

    Parameters
    ----------
    payload : dict
        Decoded JSON request body with the fields:

        stages : list or JSON text
            Stage objects with ``name``, ``specific_impulse_s`` and
            ``structural_fraction``. Alternatively the two-stage form fields
            ``stage1_*`` and ``stage2_*``.
        target_delta_v_m_s : float
            Required total velocity increment, m/s.
        payload_mass_kg : float
            Payload mass, kg.
        resolution : int, default 240
            Grid resolution of the allocation search.

    Returns
    -------
    dict
        Fields of the staging result.

    Raises
    ------
    InputError
        If a required field is missing or a numeric field is not a finite
        number.
    DomainError
        If the values are outside the model domain.
        Also raised when ``stages`` is not valid JSON or not a list of objects.

    Notes
    -----
    Tsiolkovsky equation per stage [1]_ with instantaneous staging and no
    gravity or drag loss [2]_; the optimum is global on the declared grid only.

    References
    ----------
    .. [1] NASA-ROCKET. NASA Glenn Research Center, "Ideal Rocket Equation",
           https://www1.grc.nasa.gov/beginners-guide-to-aeronautics/ideal-rocket-equation/
    .. [2] SUTTON-BIBLARZ-2017. G. P. Sutton and O. Biblarz, "Rocket Propulsion
           Elements", 9th ed., Wiley, 2017.
    """

    raw_stages = payload.get("stages")
    if isinstance(raw_stages, str):
        try:
            raw_stages = json.loads(raw_stages)
        except (ValueError, RecursionError) as error:  # malformed, oversized or too deep
            raise InputError("stages must be valid JSON.", field="stages") from error
    if raw_stages is None and "stage1_specific_impulse_s" in payload:
        raw_stages = [
            {
                "name": payload.get("stage1_name", "Booster"),
                "specific_impulse_s": payload.get("stage1_specific_impulse_s"),
                "structural_fraction": payload.get("stage1_structural_fraction"),
            },
            {
                "name": payload.get("stage2_name", "Upper"),
                "specific_impulse_s": payload.get("stage2_specific_impulse_s"),
                "structural_fraction": payload.get("stage2_structural_fraction"),
            },
        ]
    if not isinstance(raw_stages, list):
        raise InputError("stages must be a list.", field="stages")
    stages: list[StageDefinition] = []
    for index, raw_stage in enumerate(raw_stages):
        if not isinstance(raw_stage, dict):
            raise InputError("Each stage must be an object.", field=f"stages[{index}]")
        stages.append(
            StageDefinition(
                name=raw_stage.get("name", f"Stage {index + 1}"),
                specific_impulse_s=_number(raw_stage, "specific_impulse_s"),
                structural_fraction=_number(raw_stage, "structural_fraction"),
            )
        )
    return asdict(
        optimize_staging(
            stages=tuple(stages),
            target_delta_v_m_s=_number(payload, "target_delta_v_m_s"),
            payload_mass_kg=_number(payload, "payload_mass_kg"),
            resolution=_integer(payload, "resolution", default=240),
        )
    )


def calculate_engine_cycle_route(payload: dict[str, Any]) -> dict[str, Any]:
    """Close a preliminary feed-system pump and turbine energy balance.

    Parameters
    ----------
    payload : dict
        Decoded JSON request body with the fields:

        cycle : str, default "pressure-fed"
            Engine cycle.
        chamber_pressure_pa, chamber_mass_flow_kg_s, inlet_pressure_pa : float
            Chamber and feed conditions.
        injector_pressure_drop_fraction : float
            Injector pressure drop as a fraction of chamber pressure.
        propellant_density_kg_m3, pump_efficiency, turbine_efficiency : float, optional
            Pump-side inputs.
        turbine_inlet_temperature_k, turbine_specific_heat_j_kg_k : float, optional
            Turbine-side inputs.
        turbine_gamma, turbine_pressure_ratio : float, optional
            Turbine-side inputs.
        pump_feed_policy : str, default "shared-pump"
            Whether generator flow passes through the main pumps.

    Returns
    -------
    dict
        Fields of the engine-cycle result.

    Raises
    ------
    InputError
        If a required field is missing or a numeric field is not a finite
        number.
    DomainError
        If the values are outside the model domain.

    Notes
    -----
    Cycle power balances with whole pumped mass flows [1]_ and ideal-gas
    turbine temperature, pressure and shaft-work relations [2]_. No
    turbomachinery maps or cavitation model.

    References
    ----------
    .. [1] NASA-TMX64749-CYCLE. R. H. Schmucker, 1973, equivalent cycle power
           balances and whole pumped mass flows,
           https://ntrs.nasa.gov/citations/19730016059
    .. [2] NASA-ENGINESIM-THEORY. NASA Glenn Research Center, ideal gas turbine
           temperature/pressure and shaft-work relations,
           https://www.grc.nasa.gov/www/k-12/airplane/EngineTheory.pdf
    """

    return asdict(
        calculate_engine_cycle(
            cycle=payload.get("cycle", "pressure-fed"),
            chamber_pressure_pa=_number(payload, "chamber_pressure_pa"),
            chamber_mass_flow_kg_s=_number(payload, "chamber_mass_flow_kg_s"),
            propellant_density_kg_m3=_optional_number(payload, "propellant_density_kg_m3"),
            inlet_pressure_pa=_number(payload, "inlet_pressure_pa"),
            injector_pressure_drop_fraction=_number(
                payload, "injector_pressure_drop_fraction"
            ),
            pump_efficiency=_optional_number(payload, "pump_efficiency"),
            turbine_efficiency=_optional_number(payload, "turbine_efficiency"),
            turbine_inlet_temperature_k=_optional_number(
                payload, "turbine_inlet_temperature_k"
            ),
            turbine_specific_heat_j_kg_k=_optional_number(
                payload, "turbine_specific_heat_j_kg_k"
            ),
            turbine_gamma=_optional_number(payload, "turbine_gamma"),
            turbine_pressure_ratio=_optional_number(payload, "turbine_pressure_ratio"),
            pump_feed_policy=payload.get("pump_feed_policy", "shared-pump"),
        )
    )


def _study_request_integer(payload, key, default):
    """Read an integer field exactly, including large seeds given as text.

    JSON integers are returned unchanged; floats and numeric text must be
    finite and integral. ``bool`` and other types raise :class:`InputError`.
    """
    value = payload.get(key, default)
    if isinstance(value, bool) or not isinstance(value, (int, float, str)):
        raise InputError(f'{key} must be an integer.', field=key)
    if isinstance(value, int):
        return value
    try:
        number = Decimal(value) if isinstance(value, str) else Decimal.from_float(value)
        if not number.is_finite() or number != number.to_integral_value():
            raise InputError(f'{key} must be a finite integer.', field=key)
        return int(number)
    except (InvalidOperation, ValueError, OverflowError) as error:
        raise InputError(f'{key} must be an integer.', field=key) from error


def _uncertainty_delta_v(sample):
    """Return ``g_n Isp ln(MR)`` for one uncertainty sample, rounded once.

    The product of ``STANDARD_GRAVITY_M_S2``, the specific impulse and the
    platform ``math.log`` of the mass ratio is formed exactly and rounded
    once [1]_, [2]_.

    References
    ----------
    .. [1] NASA-ROCKET. NASA Glenn Research Center, "Ideal Rocket Equation",
           https://www1.grc.nasa.gov/beginners-guide-to-aeronautics/ideal-rocket-equation/
    .. [2] BIPM-SI. BIPM, "The International System of Units (SI Brochure)",
           9th ed., 2019, https://www.bipm.org/en/publications/si-brochure
    """
    return _finite_rational(Fraction(STANDARD_GRAVITY_M_S2)
                            * Fraction(sample['specific_impulse_s'])
                            * Fraction(log(sample['mass_ratio'])), 'Uncertainty delta-v')


def _validate_uncertainty_physical_ranges(parameters):
    """Require positive Isp and a mass ratio of at least one over each whole range."""
    for parameter in parameters:
        if ((parameter.name == 'specific_impulse_s' and parameter.minimum <= 0)
                or (parameter.name == 'mass_ratio' and parameter.minimum < 1)):
            raise InputError('Isp must be positive and mass ratio at least one over the entire range.',
                             field=parameter.name)


def calculate_uncertainty_route(payload: dict[str, Any]) -> dict[str, Any]:
    """Propagate bounded Isp and mass-ratio uncertainty through the Tsiolkovsky delta-v.

    Parameters
    ----------
    payload : dict
        Decoded JSON request body with the fields:

        parameters : list or JSON text
            Exactly the ranges ``specific_impulse_s`` and ``mass_ratio`` as objects
            with ``name``, ``minimum`` and ``maximum``; alternatively the form
            fields ``minimum_/maximum_specific_impulse_s`` and
            ``minimum_/maximum_mass_ratio``.
        sample_count : int, default 200
            Number of samples.
        seed : int, default 42
            Random seed (kept exact for large values).
        method : str, default "latin-hypercube"
            Sampling method.

    Returns
    -------
    dict
        Fields of the uncertainty result (samples, summaries, correlations).

    Raises
    ------
    InputError
        If a required field is missing or a numeric field is not a finite
        number.
    DomainError
        If the values are outside the model domain.
        Also raised for malformed or missing parameter ranges, or Isp not
        positive and mass ratio below one anywhere in the ranges.

    Notes
    -----
    Independent uniform marginals sampled by Latin hypercube [1]_ through
    ``dv = g_n Isp ln(MR)`` [2]_; Pearson coefficients [3]_ are linear
    screening only. Summaries are empirical, without qualified confidence
    coverage.

    References
    ----------
    .. [1] MCKAY-1979. M. D. McKay, R. J. Beckman and W. J. Conover, "A Comparison
           of Three Methods for Selecting Values of Input Variables in the Analysis
           of Output from a Computer Code", Technometrics 21(2), 239-245, 1979,
           doi:10.2307/1268522
    .. [2] NASA-ROCKET. NASA Glenn Research Center, "Ideal Rocket Equation",
           https://www1.grc.nasa.gov/beginners-guide-to-aeronautics/ideal-rocket-equation/
    .. [3] NIST-PEARSON. NIST Dataplot, centered-product correlation definition,
           https://www.itl.nist.gov/div898/software/dataplot/refman2/auxillar/correlat.htm
    """

    raw_parameters = payload.get("parameters")
    if isinstance(raw_parameters, str):
        try:
            raw_parameters = json.loads(raw_parameters)
        except (ValueError, RecursionError) as error:  # malformed, oversized or too deep
            raise InputError("parameters must be valid JSON.", field="parameters") from error
    if raw_parameters is None and "minimum_specific_impulse_s" in payload:
        raw_parameters = [
            {
                "name": "specific_impulse_s",
                "minimum": payload.get("minimum_specific_impulse_s"),
                "maximum": payload.get("maximum_specific_impulse_s"),
            },
            {
                "name": "mass_ratio",
                "minimum": payload.get("minimum_mass_ratio"),
                "maximum": payload.get("maximum_mass_ratio"),
            },
        ]
    if not isinstance(raw_parameters, list):
        raise InputError("parameters must be a list.", field="parameters")
    parameters = tuple(
        ParameterRange(
            name=item.get("name", ""),
            minimum=_number(item, "minimum"),
            maximum=_number(item, "maximum"),
        )
        for item in raw_parameters
        if isinstance(item, dict)
    )
    if len(parameters) != len(raw_parameters):
        raise InputError("Each uncertainty parameter must be an object.")
    required = {"specific_impulse_s", "mass_ratio"}
    if {parameter.name for parameter in parameters} != required:
        raise InputError(
            "rocket-delta-v uncertainty requires specific_impulse_s and mass_ratio ranges.",
            field="parameters",
        )
    _validate_uncertainty_physical_ranges(parameters)
    result = propagate_uncertainty(
        parameters,
        evaluator=_uncertainty_delta_v,
        sample_count=_study_request_integer(payload, "sample_count", 200),
        seed=_study_request_integer(payload, "seed", 42),
        method=payload.get("method", "latin-hypercube"),
    )
    return asdict(result)


def validate_workspace_route(payload: dict[str, Any]) -> dict[str, Any]:
    """Migrate and normalize a ``.rplab.json`` workspace; optionally build reports.

    Parameters
    ----------
    payload : dict
        Request body with ``workspace_json`` (JSON text) and optional
        ``include_report`` (bool, default false).

    Returns
    -------
    dict
        ``workspace`` (normalized fields), ``normalized_json`` and, when
        requested, ``report``.

    Raises
    ------
    InputError
        If ``workspace_json`` is not a string or not a valid workspace
        document.

    References
    ----------
    .. [1] IETF RFC 8259, "The JavaScript Object Notation (JSON) Data
           Interchange Format", https://www.rfc-editor.org/rfc/rfc8259
    """

    raw_text = payload.get("workspace_json")
    if not isinstance(raw_text, str):
        raise InputError("workspace_json must be a JSON string.", field="workspace_json")
    workspace = workspace_from_json(raw_text)
    result: dict[str, Any] = {
        "workspace": asdict(workspace),
        "normalized_json": workspace_to_json(workspace),
    }
    if _boolean(payload, "include_report", default=False):
        result["report"] = asdict(build_report_bundle(workspace))
    return result


def compare_workspace_cases_route(payload: dict[str, Any]) -> dict[str, Any]:
    """Compare two cases selected from a validated workspace.

    Parameters
    ----------
    payload : dict
        Request body with ``workspace_json`` (JSON text),
        ``baseline_case_id``, ``candidate_case_id`` and optional ``metrics``
        (list of metric names).

    Returns
    -------
    dict
        Fields of the case comparison (baseline-relative differences of
        shared finite scalar metrics).

    Raises
    ------
    InputError
        If the workspace is invalid, a case id does not exist, or
        ``metrics`` is not a list of strings.
    """

    raw_text = payload.get("workspace_json")
    if not isinstance(raw_text, str):
        raise InputError("workspace_json must be a JSON string.", field="workspace_json")
    workspace = workspace_from_json(raw_text)
    baseline_id = str(payload.get("baseline_case_id", ""))
    candidate_id = str(payload.get("candidate_case_id", ""))
    cases = {case.case_id: case for case in workspace.cases}
    if baseline_id not in cases or candidate_id not in cases:
        raise InputError("Baseline and candidate case ids must exist in the workspace.")
    raw_metrics = payload.get("metrics")
    metrics = None
    if raw_metrics is not None:
        if not isinstance(raw_metrics, list) or not all(
            isinstance(item, str) for item in raw_metrics
        ):
            raise InputError("metrics must be a list of strings.", field="metrics")
        metrics = tuple(raw_metrics)
    return asdict(compare_cases(cases[baseline_id], cases[candidate_id], metrics))


def list_propellants_route(payload: dict[str, Any]) -> dict[str, Any]:
    """Return preliminary propellant screening records for a category.

    Parameters
    ----------
    payload : dict
        Request body with optional ``category`` (default ``"all"``).

    Returns
    -------
    dict
        ``{"category": ..., "count": ..., "records": [...]}``.

    Raises
    ------
    DomainError
        If the category is not recognized.

    Notes
    -----
    Records are broad preliminary-design ranges [1]_, [2]_, not equilibrium
    results or test data.

    References
    ----------
    .. [1] NASA KSC Engineering, "Examples of Specific Impulse",
           https://ntrs.nasa.gov/api/citations/20100036816/downloads/20100036816.pdf
    .. [2] NASA Technology Taxonomy, "TX01 Propulsion Systems",
           https://www.nasa.gov/wp-content/uploads/2015/03/2020_nasa_technology_taxonomy_lowres.pdf
    """

    category = str(payload.get("category", "all"))
    records = [asdict(record) for record in list_propellants(category)]
    return {"category": category.lower(), "count": len(records), "records": records}


def list_species_route(payload: dict[str, Any]) -> dict[str, Any]:
    """Return the source-labelled variable-property species catalog.

    Parameters
    ----------
    payload : dict
        Ignored; accepted for a uniform route signature.

    Returns
    -------
    dict
        ``{"count": ..., "records": [...]}`` with key, name, formula, molar
        mass, polynomial temperature range and source of each species.

    References
    ----------
    .. [1] NASA/TP-2002-211556. NASA Glenn, "NASA Glenn Coefficients for Calculating
           Thermodynamic Properties of Individual Species",
           https://ntrs.nasa.gov/citations/20020085330
    """

    del payload
    records = [
        {
            "key": species.key,
            "name": species.name,
            "formula": dict(species.formula),
            "molecular_mass_kg_mol": species.molecular_mass_kg_mol,
            "minimum_temperature_k": species.minimum_temperature_k,
            "maximum_temperature_k": species.maximum_temperature_k,
            "source": species.source,
            "source_url": species.source_url,
        }
        for species in species_catalog().values()
    ]
    return {"count": len(records), "records": records}


def _composition(payload: dict[str, Any]) -> dict[str, float] | None:
    """Read the optional ``composition`` object as ``{species: finite fraction}``."""
    raw_composition = payload.get("composition")
    if raw_composition is None:
        return None
    if not isinstance(raw_composition, dict):
        raise InputError("Field composition must be an object.", field="composition")
    result: dict[str, float] = {}
    for key, raw_value in raw_composition.items():
        if isinstance(raw_value, bool):
            raise InputError(
                f"Composition fraction for {key} must be a number.", field="composition"
            )
        try:
            value = float(raw_value)
        except (TypeError, ValueError) as error:
            raise InputError(
                f"Composition fraction for {key} must be a number.", field="composition"
            ) from error
        if not isfinite(value):
            raise InputError(
                f"Composition fraction for {key} must be finite.", field="composition"
            )
        result[str(key)] = value
    return result


def calculate_thermal_properties_route(payload: dict[str, Any]) -> dict[str, Any]:
    """Evaluate a species or frozen ideal mixture and return a temperature sweep.

    Parameters
    ----------
    payload : dict
        Decoded JSON request body with the fields:

        temperature_k : float
            Evaluation temperature, K.
        pressure_pa : float, default 100000.0
            Mixture pressure, Pa.
        species : str, default "H2O"
            Species key, used when ``composition`` is absent.
        composition : object, optional
            Species fractions of a frozen ideal mixture.
        basis : {"mole", "mass"}, default "mole"
            Basis of ``composition``.
        minimum_temperature_k, maximum_temperature_k : float, optional
            Sweep interval; defaults to the common polynomial range.
        point_count : int, default 81
            Number of sweep points, 21 to 401.

    Returns
    -------
    dict
        ``model``, ``valid_temperature_range_k``, ``properties``, ``curve``
        and a calorically perfect reference evaluated at 298.15 K.

    Raises
    ------
    InputError
        If a field is missing or malformed.
    DomainError
        If the sweep leaves the common polynomial range (code
        ``"property_temperature_out_of_range"``), is empty, or has an
        invalid point count.

    Notes
    -----
    Properties use the NASA Glenn seven-term polynomials [1]_ with ideal
    mixing; they are not extrapolated.

    References
    ----------
    .. [1] NASA/TP-2002-211556. NASA Glenn, "NASA Glenn Coefficients for Calculating
           Thermodynamic Properties of Individual Species",
           https://ntrs.nasa.gov/citations/20020085330
    """

    temperature_k = _number(payload, "temperature_k")
    pressure_pa = _number(payload, "pressure_pa", default=100_000.0)
    composition = _composition(payload)
    if composition is None:
        species = get_species(str(payload.get("species", "H2O")))
        properties = asdict(species.properties(temperature_k))
        minimum_temperature = species.minimum_temperature_k
        maximum_temperature = species.maximum_temperature_k
        model = {
            "kind": "species",
            "label": species.key,
            "composition": {species.key: 1.0},
            "basis": "mole",
            "source": species.source,
            "source_url": species.source_url,
        }

        def evaluate(sample_temperature: float) -> dict[str, Any]:
            return asdict(species.properties(sample_temperature))

    else:
        mixture = Mixture.from_fractions(
            composition, basis=str(payload.get("basis", "mole"))
        )
        properties = asdict(mixture.properties(temperature_k, pressure_pa=pressure_pa))
        minimum_temperature = mixture.minimum_temperature_k
        maximum_temperature = mixture.maximum_temperature_k
        model = {
            "kind": "mixture",
            "label": " + ".join(component.species.key for component in mixture.components),
            "composition": {
                component.species.key: component.mole_fraction
                for component in mixture.components
            },
            "basis": "mole",
            "source": "NASA/TP-2002-211556 species data; ideal mixing rules",
        }

        def evaluate(sample_temperature: float) -> dict[str, Any]:
            return asdict(
                mixture.properties(sample_temperature, pressure_pa=pressure_pa)
            )

    sweep_minimum = _number(
        payload, "minimum_temperature_k", default=minimum_temperature
    )
    sweep_maximum = _number(
        payload, "maximum_temperature_k", default=maximum_temperature
    )
    point_count = _integer(payload, "point_count", default=81)
    if sweep_minimum < minimum_temperature or sweep_maximum > maximum_temperature:
        raise DomainError(
            "Requested sweep exceeds the common polynomial temperature range.",
            code="property_temperature_out_of_range",
            details={
                "minimum_temperature_k": minimum_temperature,
                "maximum_temperature_k": maximum_temperature,
            },
        )
    if sweep_maximum <= sweep_minimum:
        raise DomainError("Maximum sweep temperature must exceed the minimum.")
    if not 21 <= point_count <= 401:
        raise DomainError("Point count must be between 21 and 401.", field="point_count")
    step = (sweep_maximum - sweep_minimum) / (point_count - 1)
    curve = [evaluate(sweep_minimum + index * step) for index in range(point_count)]
    reference_temperature = 298.15
    reference = evaluate(reference_temperature)
    calorically_perfect_curve = [
        {
            "temperature_k": point["temperature_k"],
            "cp_j_kg_k": reference["cp_j_kg_k"],
            "cv_j_kg_k": reference["cv_j_kg_k"],
            "gamma": reference["gamma"],
            "enthalpy_j_kg": reference["enthalpy_j_kg"]
            + reference["cp_j_kg_k"]
            * (point["temperature_k"] - reference_temperature),
        }
        for point in curve
    ]
    return {
        "model": model,
        "valid_temperature_range_k": [minimum_temperature, maximum_temperature],
        "properties": properties,
        "curve": curve,
        "calorically_perfect_reference": {
            "temperature_k": reference_temperature,
            "properties": reference,
            "curve": calorically_perfect_curve,
        },
    }


def _reactant_amounts(payload: dict[str, Any]) -> dict[str, float]:
    """Read the non-empty ``reactants`` object as ``{species: finite amount}``."""
    raw_reactants = payload.get("reactants")
    if not isinstance(raw_reactants, dict) or not raw_reactants:
        raise InputError("Field reactants must be a non-empty object.", field="reactants")
    result: dict[str, float] = {}
    for key, raw_value in raw_reactants.items():
        try:
            value = float(raw_value)
        except (TypeError, ValueError) as error:
            raise InputError(
                f"Reactant amount for {key} must be numeric.", field="reactants"
            ) from error
        if not isfinite(value):
            raise InputError(
                f"Reactant amount for {key} must be finite.", field="reactants"
            )
        result[str(key)] = value
    return result


def _candidate_species(payload: dict[str, Any]) -> tuple[str, ...] | None:
    """Read the optional non-empty ``candidate_species`` array."""
    raw_candidates = payload.get("candidate_species")
    if raw_candidates is None:
        return None
    if not isinstance(raw_candidates, list) or not raw_candidates:
        raise InputError(
            "Field candidate_species must be a non-empty array.",
            field="candidate_species",
        )
    return tuple(str(value) for value in raw_candidates)


def _equilibrium_provider(payload: dict[str, Any]):
    """Return the local provider or a CEA provider configured from the environment.

    ``RPLAB_CEA_EXECUTABLE`` names the CEA executable and
    ``RPLAB_CEA_CACHE`` an optional cache directory.
    """
    name = str(payload.get("provider", "local")).strip().lower()
    if name == "local":
        return LocalEquilibriumProvider()
    if name == "cea":
        executable = os.environ.get("RPLAB_CEA_EXECUTABLE")
        cache = os.environ.get("RPLAB_CEA_CACHE")
        adapter = CeaSubprocessAdapter(
            executable,
            cache_directory=Path(cache) if cache else None,
        )
        return CeaEquilibriumProvider(adapter)
    raise InputError(
        "Equilibrium provider must be local or cea.",
        field="provider",
        details={"supported": ["local", "cea"]},
    )


def calculate_equilibrium_route(payload: dict[str, Any]) -> dict[str, Any]:
    """Solve a provider-neutral TP, HP or UV chemical-equilibrium problem.

    Parameters
    ----------
    payload : dict
        Decoded JSON request body with the fields:

        problem_type : {"TP", "HP", "UV"}, default "HP"
            Fixed-state pair of the problem.
        reactants : object
            Non-empty ``{species: amount}`` mapping.
        candidate_species : list of str, optional
            Product species to consider.
        temperature_k, pressure_pa : float
            State for TP; pressure also for HP.
        reactant_temperature_k : float, default 298.15
            Reactant temperature for HP and UV, K.
        target_enthalpy_j, target_internal_energy_j : float, optional
            Explicit energy targets for HP and UV.
        volume_m3 : float
            Volume for UV, m^3.
        provider : {"local", "cea"}, default "local"
            Equilibrium provider.

    Returns
    -------
    dict
        Fields of the equilibrium result.

    Raises
    ------
    InputError
        If a required field is missing or a numeric field is not a finite
        number.
    DomainError
        If the values are outside the model domain.
        Also raised for an unknown problem type or provider.

    Notes
    -----
    Element-conserving ideal-gas equilibrium [1]_ with NASA Glenn species
    data [2]_; the local provider is experimental.

    References
    ----------
    .. [1] NASA RP-1311, "Computer Program for Calculation of Complex Chemical
           Equilibrium Compositions and Applications, Part I: Analysis",
           https://ntrs.nasa.gov/citations/19950013764
    .. [2] NASA/TP-2002-211556. NASA Glenn, "NASA Glenn Coefficients for Calculating
           Thermodynamic Properties of Individual Species",
           https://ntrs.nasa.gov/citations/20020085330
    """

    problem_type = str(payload.get("problem_type", "HP")).strip().upper()
    reactants = _reactant_amounts(payload)
    candidates = _candidate_species(payload)
    if problem_type == "TP":
        problem = EquilibriumProblem.tp(
            reactants,
            temperature_k=_number(payload, "temperature_k"),
            pressure_pa=_number(payload, "pressure_pa"),
            candidate_species=candidates,
        )
    elif problem_type == "HP":
        problem = EquilibriumProblem.hp(
            reactants,
            pressure_pa=_number(payload, "pressure_pa"),
            reactant_temperature_k=_number(
                payload, "reactant_temperature_k", default=298.15
            ),
            target_enthalpy_j=_optional_number(payload, "target_enthalpy_j"),
            candidate_species=candidates,
        )
    elif problem_type == "UV":
        problem = EquilibriumProblem.uv(
            reactants,
            volume_m3=_number(payload, "volume_m3"),
            reactant_temperature_k=_number(
                payload, "reactant_temperature_k", default=298.15
            ),
            target_internal_energy_j=_optional_number(
                payload, "target_internal_energy_j"
            ),
            candidate_species=candidates,
        )
    else:
        raise InputError(
            "Equilibrium problem type must be TP, HP, or UV.",
            field="problem_type",
        )
    return asdict(_equilibrium_provider(payload).solve(problem))


def list_chemistry_capabilities_route(payload: dict[str, Any]) -> dict[str, Any]:
    """List chemistry providers and source-labelled propellant presets.

    Parameters
    ----------
    payload : dict
        Ignored; accepted for a uniform route signature.

    Returns
    -------
    dict
        ``providers`` (key, availability, status and supported problems)
        and ``propellant_pairs``. CEA availability reflects whether
        ``RPLAB_CEA_EXECUTABLE`` is set in the current process.

    References
    ----------
    .. [1] NASA RP-1311, "Computer Program for Calculation of Complex Chemical
           Equilibrium Compositions and Applications, Part II",
           https://ntrs.nasa.gov/citations/19960044559
    """

    del payload
    cea_path = os.environ.get("RPLAB_CEA_EXECUTABLE")
    library = cea_library_info()
    return {
        "providers": [
            {
                "key": "local",
                "available": True,
                "status": "experimental",
                "supports": ["TP", "HP", "UV", "O/F sweep"],
            },
            {
                "key": "cea",
                "available": bool(cea_path),
                "status": "validated external adapter",
                "supports": ["TP", "HP", "UV", "rocket"],
                "configuration": "RPLAB_CEA_EXECUTABLE",
            },
            {
                "key": "cea-library",
                "available": library.available,
                "status": "NASA CEA Python package (optional)",
                "supports": ["rocket"],
                "version": library.version,
                "thermo_sha256": library.thermo_sha256,
                "configuration": 'pip install "rocket-propulsion-lab[cea]"',
            },
        ],
        "propellant_pairs": [asdict(pair) for pair in PROPELLANT_PAIRS.values()],
    }


def calculate_oxidizer_fuel_sweep_route(payload: dict[str, Any]) -> dict[str, Any]:
    """Generate experimental equilibrium and frozen-expansion curves over mixture ratio.

    Parameters
    ----------
    payload : dict
        Decoded JSON request body with the fields:

        propellant_pair : str, default "LOX_LH2"
            Propellant preset key.
        chamber_pressure_pa : float
            Chamber pressure, Pa.
        area_ratio : float, default 20.0
            Nozzle area ratio.
        minimum_ratio, maximum_ratio : float, optional
            Oxidizer-to-fuel ratio interval; defaults to the preset range.
        point_count : int, default 5
            Number of mixture ratios.
        provider : {"local"}, default "local"
            Only the local provider supports this route.

    Returns
    -------
    dict
        Fields of the O/F sweep result.

    Raises
    ------
    InputError
        If a required field is missing or a numeric field is not a finite
        number.
    DomainError
        If the values are outside the model domain.
        Also raised with code ``"provider_not_supported_for_route"`` for a
        non-local provider.

    Notes
    -----
    Ideal-gas equilibrium chamber [1]_, [2]_ with NASA Glenn species data
    [3]_; frozen chamber-property throat mass flux and one-dimensional
    inviscid expansion.

    References
    ----------
    .. [1] NASA RP-1311, "Computer Program for Calculation of Complex Chemical
           Equilibrium Compositions and Applications, Part I: Analysis",
           https://ntrs.nasa.gov/citations/19950013764
    .. [2] NASA RP-1311, "Computer Program for Calculation of Complex Chemical
           Equilibrium Compositions and Applications, Part II",
           https://ntrs.nasa.gov/citations/19960044559
    .. [3] NASA/TP-2002-211556. NASA Glenn, "NASA Glenn Coefficients for Calculating
           Thermodynamic Properties of Individual Species",
           https://ntrs.nasa.gov/citations/20020085330
    """

    provider = _equilibrium_provider(payload)
    if not isinstance(provider, LocalEquilibriumProvider):
        raise InputError(
            "The structured CEA rocket sweep requires a configured CEA rocket adapter; "
            "use the local provider in this release.",
            code="provider_not_supported_for_route",
            field="provider",
        )
    return asdict(
        generate_oxidizer_fuel_sweep(
            propellant_pair=str(payload.get("propellant_pair", "LOX_LH2")),
            chamber_pressure_pa=_number(payload, "chamber_pressure_pa"),
            area_ratio=_number(payload, "area_ratio", default=20.0),
            minimum_ratio=_optional_number(payload, "minimum_ratio"),
            maximum_ratio=_optional_number(payload, "maximum_ratio"),
            point_count=_integer(payload, "point_count", default=5),
            provider=provider,
        )
    )


ROUTES = {
    "/api/v1/isentropic": calculate_isentropic_route,
    "/api/v1/normal-shock": calculate_normal_shock_route,
    "/api/v1/oblique-shock": calculate_oblique_shock_route,
    "/api/v1/isentropic-sweep": calculate_isentropic_sweep_route,
    "/api/v1/ideal-gas": calculate_ideal_gas_route,
    "/api/v1/stagnation": calculate_stagnation_route,
    "/api/v1/thermo-process": calculate_process_route,
    "/api/v1/polytropic": calculate_polytropic_route,
    "/api/v1/duct-flow": calculate_duct_flow_route,
    "/api/v1/nozzle": calculate_nozzle_route,
    "/api/v1/nozzle-off-design": calculate_off_design_nozzle_route,
    "/api/v1/nozzle-contour": calculate_nozzle_contour_route,
    "/api/v1/rocket-equation": calculate_rocket_equation_route,
    "/api/v1/thrust": calculate_thrust_route,
    "/api/v1/burn/simulate": simulate_burn_route,
    "/api/v1/burn/report": build_burn_report_route,
    "/api/v1/burn/export/sidera": export_sidera_burn_route,
    "/api/v1/engine-operating-point": calculate_engine_operating_point_route,
    "/api/v1/integrations/sidera/capabilities": sidera_capabilities_route,
    "/api/v1/loss-budget": calculate_loss_budget_route,
    "/api/v1/nozzle-thermal": calculate_nozzle_thermal_route,
    "/api/v1/two-phase-band": calculate_two_phase_band_route,
    "/api/v1/staging": calculate_staging_route,
    "/api/v1/engine-cycle": calculate_engine_cycle_route,
    "/api/v1/uncertainty": calculate_uncertainty_route,
    "/api/v1/workspace": validate_workspace_route,
    "/api/v1/case-comparison": compare_workspace_cases_route,
    "/api/v1/propellants": list_propellants_route,
    "/api/v1/species": list_species_route,
    "/api/v1/thermal-properties": calculate_thermal_properties_route,
    "/api/v1/equilibrium": calculate_equilibrium_route,
    "/api/v1/chemistry-capabilities": list_chemistry_capabilities_route,
    "/api/v1/of-sweep": calculate_oxidizer_fuel_sweep_route,
}


_PERFECT_GAS = (
    "calorically perfect gas",
    "constant specific-heat ratio",
    "one-dimensional steady flow",
)

ROUTE_METADATA = {
    "/api/v1/isentropic": ModelMetadata(
        model="perfect-gas isentropic relations",
        version="1.0",
        assumptions=_PERFECT_GAS + ("adiabatic reversible flow",),
        units={"mach": "1", "*_ratio": "1", "*_deg": "deg"},
        validity=(ValidityRange("mach", minimum=0.0, inclusive_minimum=False),),
        references=("NASA Glenn isentropic flow equations",),
    ),
    "/api/v1/normal-shock": ModelMetadata(
        model="perfect-gas normal shock",
        version="1.0",
        assumptions=_PERFECT_GAS + ("planar normal discontinuity",),
        units={"mach": "1", "*_ratio": "1"},
        validity=(ValidityRange("upstream_mach", minimum=1.0, inclusive_minimum=False),),
        references=("NASA Glenn normal shock equations",),
    ),
    "/api/v1/oblique-shock": ModelMetadata(
        model="theta-beta-M attached oblique shock",
        version="1.0",
        assumptions=_PERFECT_GAS + ("planar attached shock",),
        units={"mach": "1", "*_ratio": "1", "*_deg": "deg"},
        validity=(ValidityRange("upstream_mach", minimum=1.0, inclusive_minimum=False),),
        references=("theta-beta-M relation",),
    ),
    "/api/v1/isentropic-sweep": ModelMetadata(
        model="perfect-gas isentropic sweep",
        version="1.0",
        assumptions=_PERFECT_GAS + ("adiabatic reversible flow",),
        units={"mach": "1", "*_ratio": "1", "*_deg": "deg"},
    ),
    "/api/v1/ideal-gas": ModelMetadata(
        model="calorically perfect ideal gas",
        version="1.0",
        assumptions=("ideal-gas equation of state", "constant cp and cv"),
        units={
            "pressure_pa": "Pa",
            "temperature_k": "K",
            "density_kg_m3": "kg/m3",
            "*_j_kg_k": "J/kg/K",
            "speed_of_sound_m_s": "m/s",
        },
    ),
    "/api/v1/stagnation": ModelMetadata(
        model="perfect-gas static-to-stagnation state",
        version="1.0",
        assumptions=_PERFECT_GAS + ("adiabatic reversible deceleration",),
        units={"*_pressure_pa": "Pa", "*_temperature_k": "K", "*_m_s": "m/s"},
    ),
    "/api/v1/thermo-process": ModelMetadata(
        model="efficiency-corrected compression or expansion",
        version="1.0",
        assumptions=("calorically perfect gas", "adiabatic steady component"),
        units={"*_pressure_pa": "Pa", "*_temperature_k": "K", "*_j_kg_k": "J/kg/K"},
    ),
    "/api/v1/polytropic": ModelMetadata(
        model="closed-system polytropic process",
        version="1.0",
        assumptions=("ideal gas", "quasistatic path", "constant exponent"),
        units={"*_pa": "Pa", "*_k": "K", "*_j_kg": "J/kg", "*_m3_kg": "m3/kg"},
    ),
    "/api/v1/duct-flow": ModelMetadata(
        model="perfect-gas Fanno or Rayleigh line",
        version="1.0",
        assumptions=_PERFECT_GAS + ("constant-area duct",),
        units={"mach": "1", "*_ratio": "1", "friction_parameter_to_sonic": "1"},
    ),
    "/api/v1/nozzle": ModelMetadata(
        model="ideal choked converging-diverging nozzle",
        version="1.0",
        assumptions=_PERFECT_GAS + ("inviscid", "isentropic", "frozen properties"),
        units={
            "*_pressure_pa": "Pa",
            "*_temperature_k": "K",
            "*_area_m2": "m2",
            "*_kg_s": "kg/s",
            "*_n": "N",
            "specific_impulse_s": "s",
        },
        validity=(ValidityRange("area_ratio", minimum=1.0),),
    ),
    "/api/v1/nozzle-off-design": ModelMetadata(
        model="perfect-gas back-pressure-driven C-D nozzle",
        version="1.0",
        assumptions=(
            "calorically perfect gas",
            "quasi-one-dimensional steady flow",
            "normal shock for internal compression",
            "inviscid core flow",
            "empirical separation screening only",
        ),
        units={
            "*_pressure_pa": "Pa",
            "*_temperature_k": "K",
            "*_area_m2": "m2",
            "*_m_s": "m/s",
            "*_kg_s": "kg/s",
            "*_n": "N",
            "specific_impulse_s": "s",
            "altitude_m": "m",
        },
        validity=(ValidityRange("area_ratio", minimum=1.0),),
        references=(
            "quasi-one-dimensional normal-shock nozzle relations",
            "NASA SP-8120 p.17 exit-pressure-ratio separation screen",
            "U.S. Standard Atmosphere 1976",
        ),
    ),
    "/api/v1/nozzle-contour": ModelMetadata(
        model="multi-fidelity axisymmetric nozzle contour",
        version="2.1",
        assumptions=(
            "axisymmetric equivalent-area geometry",
            "prescribed isentropic station states; construction and flow applicability remain under review",
            "mapped planar MOC field; full axisymmetric source-term solve not implemented",
            "exported closed sampled fluid domain, not a manufactured nozzle wall; STL coordinates in metres",
            "owned export replay and binary64 coordinate round trips do not establish physical qualification",
        ),
        units={"*_m": "m", "*_m2": "m2", "mach": "1", "*_ratio": "1"},
        references=(
            "Rao thrust-optimized parabola construction",
            "NASA method-of-characteristics nozzle design literature",
        ),
    ),
    "/api/v1/rocket-equation": ModelMetadata(
        model="Tsiolkovsky ideal rocket equation",
        version="2.0",
        assumptions=("constant effective exhaust velocity", "no gravity or drag losses",
                     "retained solve mode and performance input; independent checked exports",
                     "finite precision directed log/exp enclosures; no hardware validation"),
        units={"*_kg": "kg", "*_m_s": "m/s", "specific_impulse_s": "s", "*_ratio": "1",
               "propellant_mass_fraction": "1"},
    ),
    "/api/v1/thrust": ModelMetadata(
        model="steady rocket thrust equation",
        version="2.0",
        assumptions=("steady one-dimensional exit plane", "uniform exit properties",
                     "exact binary-input pressure/momentum ledger; retained impulse duration"),
        units={"*_kg_s": "kg/s", "*_m_s": "m/s", "*_m2": "m2", "*_pa": "Pa", "*_n": "N",
               "specific_impulse_s": "s", "burn_time_s": "s", "total_impulse_n_s": "N s"},
    ),
    "/api/v1/burn/simulate": ModelMetadata(
        model="standalone analytic constant/profiled propulsion burn",
        version="0.11.0-l1",
        assumptions=(
            "rated delivered thrust and explicit tank flows",
            "optional piecewise-linear throttle scales thrust and all streams",
            "constant system specific impulse over L1 schedules",
            "one-dimensional ideal velocity increment",
            "no gravity, drag, steering, attitude, or orbital propagation",
        ),
        units={
            "*_kg": "kg",
            "*_kg_s": "kg/s",
            "*_n": "N",
            "*_n_s": "N s",
            "*_m_s": "m/s",
            "*_pa": "Pa",
            "*_s": "s",
        },
        references=("NASA SP-125", "Tsiolkovsky ideal rocket equation"),
    ),
    "/api/v1/engine-operating-point": ModelMetadata(
        model="constant propulsion operating-point provider",
        version="1.0",
        assumptions=(
            "constant delivered thrust and system specific impulse",
            "constant total oxidizer/fuel ratio for bipropellant mode",
        ),
        units={"*_n": "N", "*_s": "s", "*_kg_s": "kg/s", "*_pa": "Pa"},
    ),
    "/api/v1/burn/report": ModelMetadata(
        model="standalone propulsion-burn HTML report",
        version="1.0",
        assumptions=("source artifact passes hash and physical reproduction checks",),
        units={},
    ),
    "/api/v1/burn/export/sidera": ModelMetadata(
        model="exact or explicitly reduced constant-burn Sidera adapter",
        version="1.1",
        assumptions=(
            "source is a verified rocket_propulsion_burn_v1 or v2 result",
            "delivered thrust is already throttle-resolved",
            "system Isp represents total tank drain",
            "equivalent_constant mode preserves duration, impulse, and propellant only",
        ),
        units={"t_start_s": "s", "duration_s": "s", "thrust_n": "N", "isp_s": "s"},
    ),
    "/api/v1/integrations/sidera/capabilities": ModelMetadata(
        model="Sidera public finite-burn capability probe",
        version="1.1",
        assumptions=("only installed public modules are inspected",
                     "checking means unknown; flags are omitted until inspection completes",
                     "process-local snapshot; restart server after installation changes"),
        units={},
    ),
    "/api/v1/loss-budget": ModelMetadata(
        model="declared operating-basis factor ratios and first-order covariance proxy",
        version="2.0",
        assumptions=(
            "default fixed-mass-flow-lumped mode retains legacy multipliers; discharge is not a mass-flow coefficient",
            "fixed-pressure-area mode separates intrinsic cstar, mass-flow Cd and net thrust-coefficient factors",
            "factor correlation matrix order: combustion, discharge, boundary_layer, divergence",
            "first-order signed covariance propagation; ideal reference fixed; no qualified coverage",
            "sequential attribution does not prove nonoverlapping physical losses; no coupled gas/nozzle closure",
        ),
        units={"*_n": "N", "*_s": "s", "*_efficiency": "1", "*_residual_*": "varies"},
    ),
    "/api/v1/nozzle-thermal": ModelMetadata(
        model="prescribed signed Bartz profile and optional experimental bulk coolant balance",
        version="2.0",
        assumptions=(
            "steady turbulent hot-gas boundary layer",
            "user-supplied frozen gas properties",
            "single-phase constant-property coolant when cooling is enabled",
            "nested coolant model 2.0: exact temperature-limit/Reynolds decisions; default laminar Darcy f*Re=64 assumes a circular equivalent",
            "coolant cross-section/transition/friction domain and boiling applicability are unqualified",
            "prescribed uniform wall temperature; wall/coolant coupling, radiation and physical applicability are unqualified",
            "piecewise linear heat-flux/radius reconstruction; no continuous-profile error bound",
            "zero gas-side residual is an integral ledger identity, not measured energy closure",
        ),
        units={"*_w_m2": "W/m2", "*_w": "W", "*_k": "K", "*_pa": "Pa", "*_m": "m"},
        references=("Bartz correlation, NASA SP-125", "Darcy-Weisbach pressure loss"),
    ),
    "/api/v1/two-phase-band": ModelMetadata(
        model="conditional gas/particle exit-state scenarios with separated pressure thrust",
        version="2.0",
        assumptions=("fixed total mass flow and prescribed exit gas state; gas response is not solved",
            "particle/gas exit velocity ratio, or explicitly selected per-unit-mass kinetic-energy ratio",
            "reference pressure thrust is separate; default zero is an assumption",
            "scenario endpoints are not validated physical bounds or statistical coverage",
            "nonpositive conditional thrust is retained; actual nozzle compatibility is not solved"),
        units={"*_n": "N", "*_s": "s", "*_fraction": "1", "*_efficiency*": "1", "*_ratio*": "1"},
    ),
    "/api/v1/staging": ModelMetadata(
        model="conditional finite-grid multi-stage initial-mass optimizer",
        version="2.0",
        assumptions=(
            "Tsiolkovsky equation per stage",
            "constant stage Isp and structural fraction",
            "instantaneous staging without gravity or drag loss",
            "global optimum on the declared grid; exact grid increments and certified cost comparisons",
            "two-grid mass difference is diagnostic; no continuous-optimum or physical qualification",
        ),
        units={"*_m_s": "m/s", "*_kg": "kg", "*_s": "s", "*_fraction": "1"},
    ),
    "/api/v1/engine-cycle": ModelMetadata(
        model="preliminary engine-cycle component energy balance",
        version="2.0",
        assumptions=(
            "single equivalent propellant density",
            "constant component efficiencies",
            "no turbomachinery maps or cavitation model",
            "shared-pump GG includes generator flow; external GG feed is explicitly excluded",
            "pressure-fed requirements sizing; hot-gas pressure and hardware feasibility not closed",
        ),
        units={"*_pa": "Pa", "*_w": "W", "*_kg_s": "kg/s", "*_j_kg": "J/kg", "*_fraction": "1"},
    ),
    "/api/v1/uncertainty": ModelMetadata(
        model="seeded bounded uncertainty propagation",
        version="2.0",
        assumptions=("declared independent uniform input marginals; LHS rows dependent",
                     "empirical sample summaries; no qualified confidence coverage",
                     "Pearson linear screening; undefined for constant vectors",
                     "ideal Tsiolkovsky relation; platform math.log and checked product"),
        units={"output_*": "m/s", "specific_impulse_s": "s", "mass_ratio": "1"},
    ),
    "/api/v1/workspace": ModelMetadata(
        model="versioned Rocket Propulsion Lab workspace schema",
        version="2",
        assumptions=("JSON finite-number contract", "deterministic sorted serialization"),
        units={},
    ),
    "/api/v1/case-comparison": ModelMetadata(
        model="workspace scalar result comparison",
        version="1.0",
        assumptions=("shared finite scalar metrics", "baseline-relative difference"),
        units={"relative_difference": "1"},
    ),
    "/api/v1/propellants": ModelMetadata(
        model="preliminary propellant screening catalog",
        version="1.0",
        assumptions=("representative performance ranges", "not an equilibrium calculation"),
        units={"vacuum_isp_*_s": "s"},
        references=("NASA propulsion screening references listed in docs/REFERENCES.md",),
    ),
    "/api/v1/species": ModelMetadata(
        model="NASA Glenn thermodynamic species catalog",
        version="1.0",
        assumptions=("ideal-gas standard state", "no polynomial extrapolation"),
        units={"molecular_mass_kg_mol": "kg/mol", "*_temperature_k": "K"},
        references=("NASA/TP-2002-211556",),
    ),
    "/api/v1/thermal-properties": ModelMetadata(
        model="NASA Glenn variable-property ideal gas",
        version="1.0",
        assumptions=(
            "ideal gas",
            "frozen composition",
            "NASA Glenn piecewise polynomial properties",
        ),
        units={
            "temperature_k": "K",
            "pressure_pa": "Pa",
            "*_j_kg_k": "J/kg/K",
            "enthalpy_j_kg": "J/kg",
            "molecular_mass_kg_mol": "kg/mol",
        },
        validity=(ValidityRange("temperature", minimum=200.0, maximum=6000.0, unit="K"),),
        references=("NASA/TP-2002-211556",),
    ),
    "/api/v1/equilibrium": ModelMetadata(
        model="provider-neutral ideal-gas chemical equilibrium",
        version="0.2-experimental",
        assumptions=(
            "element conservation",
            "ideal-gas products",
            "gas-phase bundled species for local provider",
        ),
        units={
            "temperature_k": "K",
            "pressure_pa": "Pa",
            "*_j": "J",
            "*_j_k": "J/K",
            "molecular_mass_kg_mol": "kg/mol",
        },
        validity=(ValidityRange("temperature", minimum=200.0, maximum=6000.0, unit="K"),),
        references=("NASA RP-1311 Part I", "NASA/TP-2002-211556"),
    ),
    "/api/v1/chemistry-capabilities": ModelMetadata(
        model="chemistry provider capability registry",
        version="1.0",
        assumptions=("availability reflects current process configuration",),
        units={},
        references=("NASA RP-1311 Part II",),
    ),
    "/api/v1/of-sweep": ModelMetadata(
        model="experimental equilibrium chamber with frozen/equilibrium nozzle comparison",
        version="0.2-experimental",
        assumptions=(
            "gas reactants at 298.15 K for local provider",
            "ideal-gas equilibrium chamber",
            "frozen chamber-property throat mass flux",
            "one-dimensional inviscid nozzle",
        ),
        units={
            "chamber_pressure_pa": "Pa",
            "adiabatic_flame_temperature_k": "K",
            "characteristic_velocity_m_s": "m/s",
            "*_specific_impulse_s": "s",
        },
        references=("NASA RP-1311", "NASA/TP-2002-211556"),
    ),
}


def calculation_metadata(path: str) -> dict[str, Any]:
    """Return the model metadata registered for an API path.

    Parameters
    ----------
    path : str
        Versioned API path, e.g. ``"/api/v1/nozzle"``.

    Returns
    -------
    dict
        Model name, version, assumptions, units, validity ranges and
        source labels.

    Raises
    ------
    KeyError
        If the path is not registered.
    """

    try:
        return asdict(ROUTE_METADATA[path])
    except KeyError as error:
        raise KeyError(path) from error

ROUTES.update(CHEMISTRY_ROUTES)
ROUTE_METADATA.update(CHEMISTRY_ROUTE_METADATA)
ROUTES.update(AERO_ROUTES)
ROUTE_METADATA.update(AERO_ROUTE_METADATA)
ROUTES.update(THERMAL_ROUTES)
ROUTE_METADATA.update(THERMAL_ROUTE_METADATA)


def dispatch_calculation(path: str, payload: dict[str, Any]) -> dict[str, Any]:
    """Dispatch a versioned API path to its calculation adapter.

    Parameters
    ----------
    path : str
        Versioned API path, a key of ``ROUTES``.
    payload : dict
        Decoded JSON request body.

    Returns
    -------
    dict
        JSON-ready result of the route.

    Raises
    ------
    KeyError
        If the path is not registered.
    RocketPropulsionError
        Any input, domain, convergence or availability error raised by the
        route.
    """

    route = ROUTES.get(path)
    if route is None:
        raise KeyError(path)
    return route(payload)


# --- Transient-core and SASS routes (defined in transient_routes.py) ---------
from .transient_routes import TRANSIENT_ROUTE_METADATA, TRANSIENT_ROUTES

ROUTES.update(TRANSIENT_ROUTES)
ROUTE_METADATA.update(TRANSIENT_ROUTE_METADATA)
# --- End of transient-core and SASS routes -----------------------------------

