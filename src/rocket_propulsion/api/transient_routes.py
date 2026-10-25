"""JSON request adapters for the isolated transient propulsion core and SASS sizing.

Each ``*_route`` function takes a decoded JSON request body, reads its fields
with the strict readers of :mod:`.routes`, builds the library records, calls
one library model and returns a JSON-ready dictionary with

``model``
    Model label of the library function that produced the result.
``inputs``
    Every resolved scalar input, so a result can be reproduced.
``summary``
    Scalar results.
``series``
    Plot-ready groups ``{"x", "source_point_count", "point_count",
    "columns"}``; every column of a group has the same length, at most
    ``max_points`` (default and ceiling 400). Longer histories are decimated
    deterministically to evenly spaced indices that keep the first and last
    point.
``reference_ids``, ``warnings``
    Source identifiers and model warnings reported by the library.

Every field has a default taken from the matching script in ``examples/``, so
an empty request returns the example case. Unknown fields are refused, numeric
controls that set the amount of work are bounded, and all physics lives in the
library functions each route calls. ``TRANSIENT_ROUTES`` and
``TRANSIENT_ROUTE_METADATA`` are merged into the registries of :mod:`.routes`.
"""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import asdict
from math import ceil
from typing import Any

from rocket_propulsion.core.errors import DomainError, InputError
from rocket_propulsion.core.metadata import ModelMetadata
from rocket_propulsion.propulsion.burns import (
    BlowdownTankModel,
    CalibrationScope,
    ConstantPerformanceProvider,
    DiscreteThrottlePath,
    DistributionKind,
    DynamicFeedLine,
    DynamicFeedLineConfiguration,
    DynamicFeedLineState,
    ElectricThrottleLevel,
    ElectricThrottleTable,
    EngineClusterMember,
    EngineTransientCommand,
    HotFireScaleMeasurement,
    IdealPressurantGas,
    MassJettison,
    PerformanceAxis,
    PerformanceScalePoint,
    PowerLimitedPropulsionModel,
    PressurePerformanceLaw,
    PropagatorPropulsionBridge,
    PropulsionScenario,
    PropulsionScenarioSet,
    PropulsionTimeline,
    PulseModeThruster,
    RectilinearPerformanceSurface,
    RegulatedFeedConfiguration,
    RegulatedFeedSystem,
    ScenarioEvent,
    ScenarioEventType,
    ScheduledBurn,
    StreamStateBinding,
    TabulatedBurnArtifact,
    TabulatedBurnSegment,
    ThrottleSchedule,
    ThrottleSegment,
    UncertainParameter,
    UncertaintyClass,
    build_propulsion_scenario_ensemble,
    calibrate_hot_fire_scale,
    evaluate_coupled_feed_propulsion,
    evaluate_dynamic_feed_propulsion,
    free_space_burn_metrics,
    simulate_blowdown_uncertainty,
    simulate_pressure_fed_blowdown,
)
from rocket_propulsion.propulsion.performance import STANDARD_GRAVITY_M_S2
from rocket_propulsion.studies.global_sensitivity import rank_sensitivity
from rocket_propulsion.studies.sass import project_from_dict

from .routes import _boolean, _enum_value, _integer, _number, _object, _objects

MAX_SERIES_POINTS = 400
"""Ceiling (and default) of the number of points returned per series column."""

MAX_UNCERTAINTY_SAMPLES = 256
"""Ceiling of the uncertainty sample count accepted over HTTP."""

MAX_BLOWDOWN_MESH_INTERVALS = 2000
"""Ceiling of ``initial_propellant_mass_kg / maximum_mass_step_kg`` accepted over HTTP."""

_MAX_LIST_ITEMS = 400


# ---------------------------------------------------------------------------
# Strict readers and series helpers (no physics).
# ---------------------------------------------------------------------------


def _reject_unknown(payload: dict[str, Any], allowed: set[str], label: str = "") -> None:
    """Refuse fields that the route does not read, so typos never fall back to defaults."""

    unknown = sorted(str(key) for key in payload if key not in allowed)
    if unknown:
        prefix = f"{label}." if label else ""
        raise InputError(
            f"Unknown field(s): {', '.join(prefix + key for key in unknown)}.",
            code="unknown_field",
            field=prefix + unknown[0],
        )


def _text(payload: dict[str, Any], key: str, *, default: str) -> str:
    """Read a non-blank string; anything else raises :class:`InputError`."""

    value = payload.get(key, default)
    if not isinstance(value, str) or not value.strip():
        raise InputError(f"Field {key} must be a non-empty string.", field=key)
    return value


def _bounded_integer(
    payload: dict[str, Any], key: str, *, default: int, minimum: int, maximum: int
) -> int:
    """Read an integer and refuse values outside ``[minimum, maximum]``."""

    value = _integer(payload, key, default=default)
    if not minimum <= value <= maximum:
        raise InputError(
            f"Field {key} must be an integer from {minimum} to {maximum}.", field=key
        )
    return value


def _optional_number_or(payload: dict[str, Any], key: str, default: float | None) -> float | None:
    """Read a finite number; an absent key gives ``default`` and ``null`` gives ``None``."""

    if key not in payload:
        return default
    if payload[key] in (None, ""):
        return None
    return _number(payload, key)


def _numbers(
    payload: dict[str, Any],
    key: str,
    *,
    default: tuple[float, ...],
    minimum_length: int = 1,
    maximum_length: int = _MAX_LIST_ITEMS,
) -> tuple[float, ...]:
    """Read a JSON array of finite numbers with a bounded length."""

    raw = payload.get(key, default)
    if not isinstance(raw, (list, tuple)) or not minimum_length <= len(raw) <= maximum_length:
        raise InputError(
            f"Field {key} must be a list of {minimum_length} to {maximum_length} numbers.",
            field=key,
        )
    return tuple(
        _number({f"{key}[{index}]": value}, f"{key}[{index}]") for index, value in enumerate(raw)
    )


def _items(
    payload: dict[str, Any], key: str, *, default: tuple[dict[str, Any], ...]
) -> list[dict[str, Any]]:
    """Read a non-empty, bounded JSON array of objects (``default`` when absent)."""

    raw = payload[key] if key in payload else list(default)
    items = _objects(raw, key)
    if len(items) > _MAX_LIST_ITEMS:
        raise InputError(f"Field {key} accepts at most {_MAX_LIST_ITEMS} items.", field=key)
    return items


def _max_points(payload: dict[str, Any]) -> int:
    """Read ``max_points``, the per-series point ceiling (2 to 400, default 400)."""

    return _bounded_integer(
        payload, "max_points", default=MAX_SERIES_POINTS, minimum=2, maximum=MAX_SERIES_POINTS
    )


def _sweep_points(payload: dict[str, Any], default: int = 101) -> int:
    """Read ``sweep_points``, the number of evaluations of a sweep (2 to 400)."""

    return _bounded_integer(
        payload, "sweep_points", default=default, minimum=2, maximum=MAX_SERIES_POINTS
    )


def _decimated_indices(count: int, limit: int) -> list[int]:
    """Return evenly spaced indices of ``range(count)`` that keep both end points.

    The spacing ``(count - 1) / (limit - 1)`` is at least one, so the integer
    floors are strictly increasing and the selection is deterministic.
    """

    if count <= limit:
        return list(range(count))
    return [(index * (count - 1)) // (limit - 1) for index in range(limit)]


def _series(x: str, columns: dict[str, list[Any]], limit: int) -> dict[str, Any]:
    """Pack equal-length columns into one plot-ready, decimated series group."""

    count = len(columns[x])
    indices = _decimated_indices(count, limit)
    return {
        "x": x,
        "source_point_count": count,
        "point_count": len(indices),
        "columns": {name: [values[index] for index in indices] for name, values in columns.items()},
    }


def _linspace(start: float, stop: float, count: int) -> list[float]:
    """Return ``count`` evenly spaced sweep abscissae from ``start`` to ``stop``."""

    return [start + (stop - start) * index / (count - 1) for index in range(count)]


def _artifact_series(artifact: TabulatedBurnArtifact, limit: int) -> dict[str, Any]:
    """Tabulate the segment end points of a burn artifact (two points per segment)."""

    columns: dict[str, list[Any]] = {
        "time_s": [],
        "delivered_thrust_n": [],
        "axial_thrust_n": [],
        "total_mass_flow_kg_s": [],
    }
    for segment in artifact.segments:
        for end in ("start", "end"):
            columns["time_s"].append(getattr(segment, f"{end}_time_s"))
            columns["delivered_thrust_n"].append(getattr(segment, f"{end}_delivered_thrust_n"))
            columns["axial_thrust_n"].append(getattr(segment, f"{end}_axial_thrust_n"))
            columns["total_mass_flow_kg_s"].append(getattr(segment, f"{end}_total_mass_flow_kg_s"))
    return _series("time_s", columns, limit)


def _artifact_summary(artifact: TabulatedBurnArtifact) -> dict[str, Any]:
    """Scalar integrals that a tabulated burn artifact reports about itself."""

    return {
        "duration_s": artifact.duration_s,
        "delivered_impulse_n_s": artifact.delivered_total_impulse_n_s,
        "axial_impulse_n_s": artifact.axial_total_impulse_n_s,
        "consumed_propellant_kg": artifact.consumed_propellant_kg,
        "thrust_centroid_time_s": artifact.thrust_centroid_time_s,
        "system_equivalent_specific_impulse_s": artifact.system_equivalent_specific_impulse_s,
        "segment_count": len(artifact.segments),
        "artifact_hash": artifact.artifact_hash,
    }


def _histogram(values: list[float], bins: int) -> dict[str, list[float] | list[int]]:
    """Equal-width histogram of sample values; a constant sample gives one bin."""

    low, high = min(values), max(values)
    if high <= low:
        return {"edges": [low, high], "counts": [len(values)]}
    edges = _linspace(low, high, bins + 1)
    counts = [0] * bins
    for value in values:
        counts[min(bins - 1, int((value - low) / (high - low) * bins))] += 1
    return {"edges": edges, "counts": counts}


# ---------------------------------------------------------------------------
# Pressure-fed blowdown (examples/pressure_fed_blowdown.py)
# ---------------------------------------------------------------------------

_BLOWDOWN_DEFAULTS: dict[str, float | None] = {
    "delivered_thrust_n": 100.0,
    "system_specific_impulse_s": 220.0,
    "chamber_pressure_pa": 1.5e6,
    "initial_pressure_pa": 2.0e6,
    "initial_ullage_volume_m3": 0.01,
    "initial_propellant_mass_kg": 10.0,
    "propellant_density_kg_m3": 1000.0,
    "polytropic_exponent": 1.2,
    "reference_supply_pressure_pa": 2.0e6,
    "minimum_supply_pressure_pa": 1.0e6,
    "maximum_supply_pressure_pa": None,
    "thrust_pressure_exponent": 1.0,
    "specific_impulse_pressure_exponent": 0.02,
    "chamber_pressure_exponent": 1.0,
    "reserve_mass_kg": 0.0,
    "cant_efficiency": 1.0,
}
_BLOWDOWN_MODEL_FIELDS = tuple(_BLOWDOWN_DEFAULTS)


def _blowdown_records(payload: dict[str, Any], defaults: dict[str, float | None]):
    """Build the rated point, tank and pressure law shared by blowdown routes."""

    values: dict[str, Any] = {
        key: (
            _optional_number_or(payload, key, default)
            if key == "maximum_supply_pressure_pa"
            else _number(payload, key, default=default)
        )
        for key, default in defaults.items()
    }
    rated = ConstantPerformanceProvider.monopropellant(
        delivered_thrust_n=values["delivered_thrust_n"],
        system_specific_impulse_s=values["system_specific_impulse_s"],
        chamber_pressure_pa=values["chamber_pressure_pa"],
        source="illustrative rated point",
    ).operating_point()
    tank = BlowdownTankModel(
        initial_pressure_pa=values["initial_pressure_pa"],
        initial_ullage_volume_m3=values["initial_ullage_volume_m3"],
        initial_propellant_mass_kg=values["initial_propellant_mass_kg"],
        propellant_density_kg_m3=values["propellant_density_kg_m3"],
        polytropic_exponent=values["polytropic_exponent"],
    )
    law = PressurePerformanceLaw(
        reference_supply_pressure_pa=values["reference_supply_pressure_pa"],
        minimum_supply_pressure_pa=values["minimum_supply_pressure_pa"],
        thrust_pressure_exponent=values["thrust_pressure_exponent"],
        specific_impulse_pressure_exponent=values["specific_impulse_pressure_exponent"],
        chamber_pressure_exponent=values["chamber_pressure_exponent"],
        maximum_supply_pressure_pa=values["maximum_supply_pressure_pa"],
    )
    return values, rated, tank, law


def _mass_step(payload: dict[str, Any], values: dict[str, Any], default: float) -> float:
    """Read the blowdown mesh step and bound the initial mesh size."""

    step = _number(payload, "maximum_mass_step_kg", default=default)
    mass = values["initial_propellant_mass_kg"]
    if step > 0.0 and mass > 0.0 and ceil(mass / step) > MAX_BLOWDOWN_MESH_INTERVALS:
        raise InputError(
            "Field maximum_mass_step_kg is too small for an HTTP request: "
            f"initial_propellant_mass_kg / maximum_mass_step_kg must not exceed "
            f"{MAX_BLOWDOWN_MESH_INTERVALS}.",
            code="work_limit",
            field="maximum_mass_step_kg",
        )
    return step


def _tolerance(payload: dict[str, Any], key: str, *, default: float, minimum: float) -> float:
    """Read a relative tolerance and refuse values tighter than ``minimum``."""

    value = _number(payload, key, default=default)
    if value < minimum:
        raise InputError(
            f"Field {key} must be at least {minimum:g} for an HTTP request.",
            code="work_limit",
            field=key,
        )
    return value


def simulate_blowdown_route(payload: dict[str, Any]) -> dict[str, Any]:
    """Simulate an unregulated pressure-fed burn to reserve or pressure cutoff.

    Parameters
    ----------
    payload : dict
        Decoded JSON request body; every field is optional:

        delivered_thrust_n, system_specific_impulse_s, chamber_pressure_pa : float
            Rated monopropellant point at the reference pressure; defaults
            100 N, 220 s and 1.5 MPa.
        initial_pressure_pa, initial_ullage_volume_m3, initial_propellant_mass_kg,
        propellant_density_kg_m3, polytropic_exponent : float
            Rigid tank; defaults 2.0 MPa, 0.01 m3, 10 kg, 1000 kg/m3 and 1.2.
        reference_supply_pressure_pa, minimum_supply_pressure_pa : float
            Pressure-law reference and cutoff; defaults 2.0 MPa and 1.0 MPa.
        maximum_supply_pressure_pa : float or null, default null
            Upper validity of the pressure law (null: the reference pressure).
        thrust_pressure_exponent, specific_impulse_pressure_exponent,
        chamber_pressure_exponent : float
            Pressure-law exponents; defaults 1.0, 0.02 and 1.0.
        reserve_mass_kg, cant_efficiency : float
            Unburned reserve and axial fraction; defaults 0 kg and 1.
        maximum_mass_step_kg, relative_tolerance : float
            Mesh controls; defaults 0.5 kg and 1e-6 (the integration controls
            of ``examples/blowdown_uncertainty.py``); the tolerance must be at
            least 1e-9 and the initial mesh at most 2000 intervals.
        maximum_refinements : int, default 12
            Step halvings, 0 to 12.
        max_points : int, default 400
            Series point ceiling.

    Returns
    -------
    dict
        ``summary`` (duration, impulse, propellant, centroid, equivalent Isp,
        initial/final pressure and thrust, cutoff reason, artifact hash),
        ``evidence`` (mesh convergence record), ``series.history`` (time
        versus remaining mass, ullage volume, supply pressure, thrust, Isp
        and mass flow), ``reference_ids`` and ``warnings``.

    Raises
    ------
    InputError
        If a field is unknown, not a finite number, or exceeds a work limit.
    DomainError
        If the library refuses the inputs, the initial pressure lies outside
        the pressure law, or the mesh does not converge.

    Notes
    -----
    Calls
    :func:`~rocket_propulsion.propulsion.burns.simulate_pressure_fed_blowdown`:
    polytropic ullage expansion [1]_, pressure-correlated thrust and Isp [2]_,
    and the unregulated pressure decay documented in [3]_.

    References
    ----------
    .. [1] NASA-SP-8112. NASA, "Pressurization Systems for Liquid Rockets", 1975.
           https://ntrs.nasa.gov/citations/19760015212
    .. [2] NASA-CR-122347. NASA, "Monopropellant Hydrazine Resisto Jet".
           https://ntrs.nasa.gov/citations/19720011120
    .. [3] NASA-CR-131400. NASA, "Reliability Model of a Monopropellant
           Auxiliary Propulsion System". https://ntrs.nasa.gov/citations/19730012094
    """

    _reject_unknown(
        payload,
        set(_BLOWDOWN_MODEL_FIELDS)
        | {"maximum_mass_step_kg", "relative_tolerance", "maximum_refinements", "max_points"},
    )
    limit = _max_points(payload)
    values, rated, tank, law = _blowdown_records(payload, _BLOWDOWN_DEFAULTS)
    values["maximum_mass_step_kg"] = _mass_step(payload, values, 0.5)
    values["relative_tolerance"] = _tolerance(
        payload, "relative_tolerance", default=1e-6, minimum=1e-9
    )
    values["maximum_refinements"] = _bounded_integer(
        payload, "maximum_refinements", default=12, minimum=0, maximum=12
    )
    result = simulate_pressure_fed_blowdown(
        tank,
        rated,
        law,
        reserve_mass_kg=values["reserve_mass_kg"],
        cant_efficiency=values["cant_efficiency"],
        maximum_mass_step_kg=values["maximum_mass_step_kg"],
        relative_tolerance=values["relative_tolerance"],
        maximum_refinements=values["maximum_refinements"],
    )
    samples = result.samples
    columns: dict[str, list[Any]] = {
        name: [getattr(sample, name) for sample in samples]
        for name in (
            "time_s",
            "remaining_propellant_mass_kg",
            "ullage_volume_m3",
            "supply_pressure_pa",
            "delivered_thrust_n",
            "axial_thrust_n",
            "system_specific_impulse_s",
            "total_mass_flow_kg_s",
        )
    }
    return {
        "model": "pressure-fed blowdown, preliminary L2",
        "inputs": values,
        "summary": {
            **_artifact_summary(result.artifact),
            "cutoff_reason": result.cutoff_reason,
            "initial_supply_pressure_pa": samples[0].supply_pressure_pa,
            "final_supply_pressure_pa": samples[-1].supply_pressure_pa,
            "initial_delivered_thrust_n": samples[0].delivered_thrust_n,
            "final_delivered_thrust_n": samples[-1].delivered_thrust_n,
            "final_remaining_propellant_kg": samples[-1].remaining_propellant_mass_kg,
        },
        "evidence": asdict(result.evidence),
        "series": {"history": _series("time_s", columns, limit)},
        "reference_ids": list(result.reference_ids),
        "warnings": list(result.warnings),
    }


# ---------------------------------------------------------------------------
# Regulated feed and dynamic feed line
# (examples/coupled_feed_propagator.py, examples/dynamic_feed_line_propagator.py)
# ---------------------------------------------------------------------------

_FEED_CONFIGURATION_DEFAULTS: dict[str, float] = {
    "tank_internal_volume_m3": 0.02,
    "propellant_density_kg_m3": 1_000.0,
    "bottle_volume_m3": 0.01,
    "regulator_set_pressure_pa": 2.1e6,
    "regulator_full_open_error_pa": 0.2e6,
    "regulator_maximum_area_m2": 1.0e-7,
    "regulator_discharge_coefficient": 0.8,
    "valve_time_constant_s": 0.05,
    "feed_line_resistance_pa_s2_kg2": 1_000.0,
    "tank_heat_transfer_w_k": 0.0,
    "tank_wall_temperature_k": 300.0,
    "bottle_heat_transfer_w_k": 0.0,
    "bottle_wall_temperature_k": 300.0,
    "minimum_propellant_mass_kg": 1.0,
    "minimum_injector_pressure_pa": 1.0e6,
    "minimum_bottle_pressure_margin_pa": 0.1e6,
    "maximum_tank_pressure_pa": 2.3e6,
}
_FEED_STATE_DEFAULTS: dict[str, float] = {
    "propellant_mass_kg": 10.0,
    "tank_pressure_pa": 2.0e6,
    "tank_gas_temperature_k": 300.0,
    "bottle_pressure_pa": 20.0e6,
    "bottle_gas_temperature_k": 300.0,
    "regulator_opening": 0.5,
}
_PRESSURANT_DEFAULTS: dict[str, float] = {
    "pressurant_gas_constant_j_kg_k": 2_077.1,
    "pressurant_gamma": 1.667,
}
_ENGINE_CONDITION_DEFAULTS: dict[str, float] = {
    "relative_time_s": 2.0,
    "vehicle_mass_kg": 100.0,
    "ambient_pressure_pa": 50_000.0,
}
_FEED_FIELDS = (
    set(_FEED_CONFIGURATION_DEFAULTS)
    | set(_FEED_STATE_DEFAULTS)
    | set(_PRESSURANT_DEFAULTS)
    | set(_ENGINE_CONDITION_DEFAULTS)
    | {"pressurant_name"}
)


def _feed_example_bridge(source_id: str, surface_sha: str, profile_id: str, profile_sha: str):
    """Pressure-dependent 500 N profile and grid of the feed examples (fixed data)."""

    surface = RectilinearPerformanceSurface(
        axes=(
            PerformanceAxis("supply_pressure_pa", "Pa", (1.0e6, 2.0e6)),
            PerformanceAxis("ambient_pressure_pa", "Pa", (0.0, 100_000.0)),
        ),
        points=(
            PerformanceScalePoint((1.0e6, 0.0), 0.8, 0.9),
            PerformanceScalePoint((1.0e6, 100_000.0), 0.7, 0.95),
            PerformanceScalePoint((2.0e6, 0.0), 1.0, 1.1),
            PerformanceScalePoint((2.0e6, 100_000.0), 0.9, 1.15),
        ),
        source_id=source_id,
        source_sha256=surface_sha,
    )
    artifact = TabulatedBurnArtifact(
        segments=(
            TabulatedBurnSegment(0.0, 20.0, 500.0, 500.0, 490.0, 490.0, 0.2, 0.2, phase="steady"),
        ),
        source_id=profile_id,
        source_sha256=profile_sha,
    )
    return PropagatorPropulsionBridge(artifact, performance_surface=surface)


def _feed_records(payload: dict[str, Any], line_resistance_default: float):
    """Build the regulated feed system and its initial state from the request."""

    defaults = {
        **_PRESSURANT_DEFAULTS,
        **_FEED_CONFIGURATION_DEFAULTS,
        "feed_line_resistance_pa_s2_kg2": line_resistance_default,
        **_FEED_STATE_DEFAULTS,
        **_ENGINE_CONDITION_DEFAULTS,
    }
    values: dict[str, Any] = {
        key: _number(payload, key, default=default) for key, default in defaults.items()
    }
    values["pressurant_name"] = _text(payload, "pressurant_name", default="helium")
    feed = RegulatedFeedSystem(
        RegulatedFeedConfiguration(
            gas=IdealPressurantGas(
                values["pressurant_gas_constant_j_kg_k"],
                values["pressurant_gamma"],
                values["pressurant_name"],
            ),
            **{key: values[key] for key in _FEED_CONFIGURATION_DEFAULTS},
        )
    )
    state = feed.initial_state(**{key: values[key] for key in _FEED_STATE_DEFAULTS})
    return values, feed, state


def _feed_evaluation_summary(evaluation) -> dict[str, Any]:
    """Scalar fields of a :class:`RegulatedFeedEvaluation`."""

    return {
        "tank_pressure_pa": evaluation.tank_pressure_pa,
        "bottle_pressure_pa": evaluation.bottle_pressure_pa,
        "injector_pressure_pa": evaluation.injector_pressure_pa,
        "feed_line_pressure_drop_pa": evaluation.feed_line_pressure_drop_pa,
        "regulator_command": evaluation.regulator_command,
        "regulator_mass_flow_kg_s": evaluation.regulator_mass_flow_kg_s,
        "regulator_choked": evaluation.regulator_choked,
        "tank_heat_flow_w": evaluation.tank_heat_flow_w,
        "bottle_heat_flow_w": evaluation.bottle_heat_flow_w,
        "derivative": asdict(evaluation.derivative),
    }


def _propulsion_summary(evaluation) -> dict[str, Any]:
    """Scalar fields of a :class:`PropulsionDynamicsEvaluation`."""

    return {
        "active": evaluation.active,
        "phase": evaluation.phase,
        "axial_thrust_n": evaluation.axial_thrust_n,
        "force_vector_n": list(evaluation.force_vector_n),
        "engine_mass_flow_kg_s": 0.0 - evaluation.mass_derivative_kg_s,
        "performance_force_scale": evaluation.performance_force_scale,
        "performance_mass_flow_scale": evaluation.performance_mass_flow_scale,
        "operating_conditions": dict(evaluation.operating_conditions),
        "inhibited_by": list(evaluation.inhibited_by),
    }


def _feed_sweeps(feed, values: dict[str, Any], flow: float, count: int, maximum_flow: float):
    """Evaluate the feed at swept demand and swept tank pressure (state evaluations only)."""

    state = feed.initial_state(**{key: values[key] for key in _FEED_STATE_DEFAULTS})
    demand: dict[str, list[Any]] = {
        "propellant_mass_flow_kg_s": [],
        "injector_pressure_pa": [],
        "feed_line_pressure_drop_pa": [],
        "tank_pressurant_mass_rate_kg_s": [],
    }
    for demand_flow in _linspace(0.0, maximum_flow, count):
        evaluation = feed.evaluate(state, demand_flow)
        demand["propellant_mass_flow_kg_s"].append(demand_flow)
        demand["injector_pressure_pa"].append(evaluation.injector_pressure_pa)
        demand["feed_line_pressure_drop_pa"].append(evaluation.feed_line_pressure_drop_pa)
        demand["tank_pressurant_mass_rate_kg_s"].append(
            evaluation.derivative.tank_pressurant_mass_kg_s
        )
    regulator: dict[str, list[Any]] = {
        "tank_pressure_pa": [],
        "regulator_command": [],
        "regulator_mass_flow_kg_s": [],
        "regulator_choked": [],
    }
    high = values["regulator_set_pressure_pa"]
    low = high - values["regulator_full_open_error_pa"]
    for pressure in _linspace(low, high, count):
        swept_state = feed.initial_state(
            **{**{key: values[key] for key in _FEED_STATE_DEFAULTS}, "tank_pressure_pa": pressure}
        )
        evaluation = feed.evaluate(swept_state, flow)
        regulator["tank_pressure_pa"].append(pressure)
        regulator["regulator_command"].append(evaluation.regulator_command)
        regulator["regulator_mass_flow_kg_s"].append(evaluation.regulator_mass_flow_kg_s)
        regulator["regulator_choked"].append(evaluation.regulator_choked)
    return demand, regulator


def evaluate_regulated_feed_route(payload: dict[str, Any]) -> dict[str, Any]:
    """Close a regulated pressure-fed engine operating point and sweep the feed.

    Parameters
    ----------
    payload : dict
        Decoded JSON request body; every field is optional (defaults from
        ``examples/coupled_feed_propagator.py``):

        pressurant_gas_constant_j_kg_k, pressurant_gamma, pressurant_name
            Pressurant gas; defaults 2077.1 J/kg/K, 1.667, ``"helium"``.
        tank_internal_volume_m3, propellant_density_kg_m3, bottle_volume_m3 : float
            Defaults 0.02 m3, 1000 kg/m3, 0.01 m3.
        regulator_set_pressure_pa, regulator_full_open_error_pa,
        regulator_maximum_area_m2, regulator_discharge_coefficient,
        valve_time_constant_s : float
            Defaults 2.1 MPa, 0.2 MPa, 1e-7 m2, 0.8, 0.05 s.
        feed_line_resistance_pa_s2_kg2 : float, default 1000
            Quadratic line resistance, Pa s2/kg2.
        tank_heat_transfer_w_k, tank_wall_temperature_k, bottle_heat_transfer_w_k,
        bottle_wall_temperature_k : float
            Defaults 0 W/K, 300 K, 0 W/K, 300 K.
        minimum_propellant_mass_kg, minimum_injector_pressure_pa,
        minimum_bottle_pressure_margin_pa, maximum_tank_pressure_pa : float
            Event limits; defaults 1 kg, 1.0 MPa, 0.1 MPa, 2.3 MPa.
        propellant_mass_kg, tank_pressure_pa, tank_gas_temperature_k,
        bottle_pressure_pa, bottle_gas_temperature_k, regulator_opening : float
            Feed state; defaults 10 kg, 2.0 MPa, 300 K, 20 MPa, 300 K, 0.5.
        relative_time_s, vehicle_mass_kg, ambient_pressure_pa : float
            Engine evaluation point; defaults 2 s, 100 kg, 50 kPa.
        initial_mass_flow_kg_s, relaxation : float
            Coupling iteration start and relaxation; defaults 0.2 kg/s, 0.8.
        demand_sweep_maximum_flow_kg_s : float, default 0.5
            Upper end of the demand sweep (the line flow bound of
            ``examples/dynamic_feed_line_propagator.py``).
        sweep_points : int, default 101
            Evaluations per sweep, 2 to 400.

    Returns
    -------
    dict
        ``summary`` with the coupled engine point (thrust, flow, injector
        pressure, regulator state, ODE right-hand side) and coupling
        ``evidence``; ``series.demand_sweep`` (injector pressure and line drop
        versus propellant flow at the declared state) and
        ``series.regulator_sweep`` (regulator command and gas flow versus tank
        pressure across the regulator error band at the coupled flow).

    Raises
    ------
    InputError
        If a field is unknown or not a finite number.
    DomainError
        If a record is outside its domain, the injector pressure leaves the
        performance grid, or the coupling iteration does not converge.

    Notes
    -----
    Calls :func:`~rocket_propulsion.propulsion.burns.evaluate_coupled_feed_propulsion`
    and :meth:`RegulatedFeedSystem.evaluate`. The library supplies the
    right-hand side for a propagator; it does not integrate a feed time
    history, so none is returned. The pressure-dependent 500 N profile and its
    supply/ambient-pressure grid are the fixed synthetic data of the example.
    Single-node pressurization [1]_, dynamic regulator [2]_, nodal feed network
    [3]_ and choked restriction flow [4]_.

    References
    ----------
    .. [1] NASA-SP-8112. NASA, "Pressurization Systems for Liquid Rockets", 1975.
           https://ntrs.nasa.gov/citations/19760015212
    .. [2] NASA-SP-8080. NASA, "Liquid Rocket Pressure Regulators and Valves".
           https://ntrs.nasa.gov/citations/19740002611
    .. [3] NASA-20240003493. NASA, "Nodal Modeling of Liquid Propellant Feed and
           Pressurization System", 2024. https://ntrs.nasa.gov/citations/20240003493
    .. [4] NASA-GRC-MASS-FLOW-CHOKING. NASA Glenn Research Center, "Mass Flow
           Choking". https://www.grc.nasa.gov/www/k-12/BGP/mflchk.html
    """

    _reject_unknown(
        payload,
        _FEED_FIELDS
        | {"initial_mass_flow_kg_s", "relaxation", "demand_sweep_maximum_flow_kg_s", "sweep_points"},
    )
    count = _sweep_points(payload)
    values, feed, state = _feed_records(payload, 1_000.0)
    values["initial_mass_flow_kg_s"] = _number(payload, "initial_mass_flow_kg_s", default=0.2)
    values["relaxation"] = _number(payload, "relaxation", default=0.8)
    values["demand_sweep_maximum_flow_kg_s"] = _number(
        payload, "demand_sweep_maximum_flow_kg_s", default=0.5
    )
    bridge = _feed_example_bridge(
        "example-feed-performance-grid", "4" * 64, "example-feed-coupled-profile", "5" * 64
    )
    result = evaluate_coupled_feed_propulsion(
        bridge,
        feed,
        state,
        relative_time_s=values["relative_time_s"],
        vehicle_mass_kg=values["vehicle_mass_kg"],
        direction=(1.0, 0.0, 0.0),
        pressure_condition_name="supply_pressure_pa",
        other_operating_conditions={"ambient_pressure_pa": values["ambient_pressure_pa"]},
        initial_mass_flow_kg_s=values["initial_mass_flow_kg_s"],
        relaxation=values["relaxation"],
    )
    flow = 0.0 - result.propulsion.mass_derivative_kg_s
    demand, regulator = _feed_sweeps(
        feed, values, flow, count, values["demand_sweep_maximum_flow_kg_s"]
    )
    return {
        "model": "regulated single-node feed coupled to a pressure-dependent engine",
        "inputs": values,
        "summary": {
            "propulsion": _propulsion_summary(result.propulsion),
            "feed": _feed_evaluation_summary(result.feed),
            "feed_model_hash": feed.model_hash,
        },
        "evidence": asdict(result.evidence),
        "series": {
            "demand_sweep": _series("propellant_mass_flow_kg_s", demand, MAX_SERIES_POINTS),
            "regulator_sweep": _series("tank_pressure_pa", regulator, MAX_SERIES_POINTS),
        },
        "reference_ids": [
            "NASA-SP-8112",
            "NASA-SP-8080",
            "NASA-20240003493",
            "NASA-GRC-MASS-FLOW-CHOKING",
        ],
        "warnings": [
            ("State evaluation only: the library returns the feed ODE right-hand side; "
            "no time history is integrated."),
            "Performance grid and burn profile are synthetic example data.",
        ],
    }


_LINE_DEFAULTS: dict[str, float] = {
    "line_inertance_pa_s2_kg": 1.0e6,
    "line_compliance_kg_pa": 1.0e-8,
    # examples/dynamic_feed_line_propagator.py: (2.0e6 - 1.8e6) / 0.217**2.
    "line_resistance_pa_s2_kg2": (2.0e6 - 1.8e6) / 0.217**2,
    "minimum_manifold_pressure_pa": 1.0e6,
    "maximum_manifold_pressure_pa": 2.3e6,
    "maximum_absolute_flow_kg_s": 0.5,
    "line_mass_flow_kg_s": 0.217,
    "manifold_pressure_pa": 1.8e6,
}


def evaluate_feed_line_route(payload: dict[str, Any]) -> dict[str, Any]:
    """Evaluate a regulated feed, an inertial feed line and the engine together.

    Parameters
    ----------
    payload : dict
        Decoded JSON request body; every field is optional (defaults from
        ``examples/dynamic_feed_line_propagator.py``). It accepts every
        regulated-feed field of ``/api/v1/transient/regulated-feed`` (same
        defaults, except ``feed_line_resistance_pa_s2_kg2`` = 0 because the
        dynamic line owns the resistance) and:

        line_inertance_pa_s2_kg, line_compliance_kg_pa : float
            Defaults 1e6 Pa s2/kg and 1e-8 kg/Pa.
        line_resistance_pa_s2_kg2 : float
            Default (2.0e6 - 1.8e6) / 0.217**2 Pa s2/kg2.
        minimum_manifold_pressure_pa, maximum_manifold_pressure_pa : float
            Defaults 1.0 MPa and 2.3 MPa.
        maximum_absolute_flow_kg_s : float, default 0.5
            Line flow bound; also the end of the line-flow sweep.
        line_mass_flow_kg_s, manifold_pressure_pa : float
            Line state; defaults 0.217 kg/s and 1.8 MPa.
        sweep_points : int, default 101
            Evaluations of the line-flow sweep, 2 to 400.

    Returns
    -------
    dict
        ``summary`` with the engine point, the regulated feed state and the
        line evaluation (derivatives, natural frequency, damping ratio, mass
        closure error), and ``series.line_flow_sweep``: line flow
        acceleration, manifold pressure rate, friction drop and damping ratio
        versus line flow at the declared upstream pressure and engine demand.

    Raises
    ------
    InputError
        If a field is unknown or not a finite number.
    DomainError
        If a record is outside its domain or the manifold pressure leaves the
        performance grid.

    Notes
    -----
    Calls :func:`~rocket_propulsion.propulsion.burns.evaluate_dynamic_feed_propulsion`
    and :meth:`DynamicFeedLine.evaluate`: the first inertance/compliance mode
    of the line [1]_, not a water-hammer solver [2]_, with coefficients meant
    to be calibrated against a nodal network [3]_. No time history is
    integrated.

    References
    ----------
    .. [1] NASA-19740028545. NASA, "Liquid Rocket Propellant Feedline Dynamics",
           1973. https://ntrs.nasa.gov/citations/19740028545
    .. [2] NASA-20220006583. NASA NESC Technical Bulletin 22-03, "Treatment of
           Transient Pressure Events", 2022. https://ntrs.nasa.gov/citations/20220006583
    .. [3] NASA-20240003493. NASA, "Nodal Modeling of Liquid Propellant Feed and
           Pressurization System", 2024. https://ntrs.nasa.gov/citations/20240003493
    """

    _reject_unknown(payload, _FEED_FIELDS | set(_LINE_DEFAULTS) | {"sweep_points"})
    count = _sweep_points(payload)
    values, feed, feed_state = _feed_records(payload, 0.0)
    for key, default in _LINE_DEFAULTS.items():
        values[key] = _number(payload, key, default=default)
    line = DynamicFeedLine(
        DynamicFeedLineConfiguration(
            inertance_pa_s2_kg=values["line_inertance_pa_s2_kg"],
            compliance_kg_pa=values["line_compliance_kg_pa"],
            resistance_pa_s2_kg2=values["line_resistance_pa_s2_kg2"],
            minimum_manifold_pressure_pa=values["minimum_manifold_pressure_pa"],
            maximum_manifold_pressure_pa=values["maximum_manifold_pressure_pa"],
            maximum_absolute_flow_kg_s=values["maximum_absolute_flow_kg_s"],
        )
    )
    line_state = DynamicFeedLineState(
        mass_flow_kg_s=values["line_mass_flow_kg_s"],
        manifold_pressure_pa=values["manifold_pressure_pa"],
    )
    bridge = _feed_example_bridge(
        "example-dynamic-feed-grid", "6" * 64, "example-dynamic-feed-profile", "7" * 64
    )
    result = evaluate_dynamic_feed_propulsion(
        bridge,
        feed,
        feed_state,
        line,
        line_state,
        relative_time_s=values["relative_time_s"],
        vehicle_mass_kg=values["vehicle_mass_kg"],
        direction=(1.0, 0.0, 0.0),
        pressure_condition_name="supply_pressure_pa",
        other_operating_conditions={"ambient_pressure_pa": values["ambient_pressure_pa"]},
    )
    line_result = result.feed_line
    sweep: dict[str, list[Any]] = {
        "line_mass_flow_kg_s": [],
        "line_flow_rate_kg_s2": [],
        "manifold_pressure_rate_pa_s": [],
        "friction_pressure_drop_pa": [],
        "damping_ratio": [],
    }
    for line_flow in _linspace(0.0, values["maximum_absolute_flow_kg_s"], count):
        evaluation = line.evaluate(
            DynamicFeedLineState(
                mass_flow_kg_s=line_flow, manifold_pressure_pa=values["manifold_pressure_pa"]
            ),
            upstream_pressure_pa=line_result.upstream_pressure_pa,
            engine_mass_flow_kg_s=line_result.engine_mass_flow_kg_s,
        )
        sweep["line_mass_flow_kg_s"].append(line_flow)
        sweep["line_flow_rate_kg_s2"].append(evaluation.derivative.mass_flow_kg_s2)
        sweep["manifold_pressure_rate_pa_s"].append(evaluation.derivative.manifold_pressure_pa_s)
        sweep["friction_pressure_drop_pa"].append(evaluation.friction_pressure_drop_pa)
        sweep["damping_ratio"].append(evaluation.damping_ratio)
    return {
        "model": "lumped inertance/compliance feed line with regulated feed and engine",
        "inputs": values,
        "summary": {
            "propulsion": _propulsion_summary(result.propulsion),
            "feed": _feed_evaluation_summary(result.regulated_feed),
            "line": {
                "upstream_pressure_pa": line_result.upstream_pressure_pa,
                "engine_mass_flow_kg_s": line_result.engine_mass_flow_kg_s,
                "friction_pressure_drop_pa": line_result.friction_pressure_drop_pa,
                "stored_liquid_mass_rate_kg_s": line_result.stored_liquid_mass_rate_kg_s,
                "external_mass_closure_error_kg_s": line_result.external_mass_closure_error_kg_s,
                "natural_angular_frequency_rad_s": line_result.natural_angular_frequency_rad_s,
                "damping_ratio": line_result.damping_ratio,
                "derivative": asdict(line_result.derivative),
            },
            "reverse_flow_event_extension": result.reverse_flow_event_extension,
            "line_model_hash": line.model_hash,
            "feed_model_hash": feed.model_hash,
        },
        "series": {"line_flow_sweep": _series("line_mass_flow_kg_s", sweep, MAX_SERIES_POINTS)},
        "reference_ids": ["NASA-19740028545", "NASA-20220006583", "NASA-20240003493"],
        "warnings": [
            ("State evaluation only: the library returns the line ODE right-hand side; "
            "no time history is integrated."),
            "First inertance/compliance mode only; not a water-hammer solver.",
            "Performance grid and burn profile are synthetic example data.",
        ],
    }


# ---------------------------------------------------------------------------
# Power-limited electric propulsion (examples/electric_propulsion.py)
# ---------------------------------------------------------------------------

_ELECTRIC_LEVELS = (
    {"level_id": "low", "input_power_w": 1000.0, "thrust_n": 0.040, "mass_flow_kg_s": 1.5e-6},
    {"level_id": "mid", "input_power_w": 2000.0, "thrust_n": 0.080, "mass_flow_kg_s": 2.6e-6},
    {"level_id": "high", "input_power_w": 3000.0, "thrust_n": 0.110, "mass_flow_kg_s": 3.2e-6},
)
_ELECTRIC_NUMBERS: dict[str, float] = {
    "up_switch_hysteresis_w": 50.0,
    "propellant_reserve_kg": 2.0,
    "power_margin_w": 100.0,
    "duty_cycle": 1.0,
    "available_power_w": 2400.0,
    "vehicle_mass_kg": 500.0,
    "propellant_mass_kg": 20.0,
}


def evaluate_electric_route(payload: dict[str, Any]) -> dict[str, Any]:
    """Select an electric-thruster level from available power and sweep power.

    Parameters
    ----------
    payload : dict
        Decoded JSON request body; every field is optional (defaults from
        ``examples/electric_propulsion.py``):

        levels : list of object
            Throttle table rows ``{level_id, input_power_w, thrust_n,
            mass_flow_kg_s}``; default the synthetic low/mid/high table
            (1000/2000/3000 W, 0.040/0.080/0.110 N, 1.5/2.6/3.2 mg/s).
        source_id, source_sha256 : str
            Table provenance; defaults ``"example-synthetic-throttle-table"``
            and 64 zeros.
        objective : {"thrust", "specific_impulse"}, default "thrust"
            Quantity the throttle path maximizes.
        up_switch_hysteresis_w, power_margin_w : float
            Defaults 50 W and 100 W.
        propellant_reserve_kg, duty_cycle : float
            Defaults 2 kg and 1.
        available_power_w, vehicle_mass_kg, propellant_mass_kg : float
            Evaluation point; defaults 2400 W, 500 kg, 20 kg xenon.
        sweep_minimum_power_w, sweep_maximum_power_w : float
            Power sweep; defaults 0 W and twice ``available_power_w`` (a
            presentation range, not example data).
        sweep_points : int, default 101
            Evaluations of the sweep, 2 to 400.

    Returns
    -------
    dict
        ``summary`` with the selected level, thrust, Isp, consumed power,
        propellant rate and power deficit; ``switching_surfaces`` of the
        selected level; ``levels`` with derived Isp and total efficiency;
        ``series.power_sweep`` with thrust, Isp, consumed power, mass flow and
        the selected level versus available power (``level_index`` from a
        fresh selection, ``level_index_increasing`` and
        ``level_index_decreasing`` carrying hysteresis memory along the sweep).

    Raises
    ------
    InputError
        If a field is unknown, malformed or not finite.
    DomainError
        If the table or model is invalid (for example a level whose jet power
        exceeds its input power).

    Notes
    -----
    Calls :class:`~rocket_propulsion.propulsion.burns.PowerLimitedPropulsionModel`
    ``select_level`` and ``evaluate``. Throttle tables take the form of [1]_;
    the efficiency check is the jet-to-input power ratio of [2]_; the force,
    mass-derivative and event contract follows [3]_. The default table is
    synthetic.

    References
    ----------
    .. [1] NASA-20090004685. NASA, "NASA's Evolutionary Xenon Thruster (NEXT) Ion
           Propulsion System Information Summary",
           https://ntrs.nasa.gov/citations/20090004685.
    .. [2] NASA-CR-2004-213099. NASA, "Development and Characterization of
           High-Efficiency, High-Specific Impulse Xenon Hall Thrusters".
    .. [3] OREKIT-13.1.5-PROPULSION-MODEL. Orekit 13.1.5, "PropulsionModel /
           EventDetector" API.
    """

    _reject_unknown(
        payload,
        set(_ELECTRIC_NUMBERS)
        | {
            "levels",
            "source_id",
            "source_sha256",
            "objective",
            "sweep_minimum_power_w",
            "sweep_maximum_power_w",
            "sweep_points",
        },
    )
    count = _sweep_points(payload)
    values: dict[str, Any] = {
        key: _number(payload, key, default=default) for key, default in _ELECTRIC_NUMBERS.items()
    }
    values["source_id"] = _text(payload, "source_id", default="example-synthetic-throttle-table")
    values["source_sha256"] = _text(payload, "source_sha256", default="0" * 64)
    values["objective"] = _text(payload, "objective", default="thrust")
    values["sweep_minimum_power_w"] = _number(payload, "sweep_minimum_power_w", default=0.0)
    values["sweep_maximum_power_w"] = _number(
        payload, "sweep_maximum_power_w", default=2.0 * values["available_power_w"]
    )
    if values["sweep_maximum_power_w"] <= values["sweep_minimum_power_w"]:
        raise InputError(
            "Field sweep_maximum_power_w must exceed sweep_minimum_power_w.",
            field="sweep_maximum_power_w",
        )
    level_items = _items(payload, "levels", default=_ELECTRIC_LEVELS)
    levels = []
    for index, item in enumerate(level_items):
        _reject_unknown(
            item, {"level_id", "input_power_w", "thrust_n", "mass_flow_kg_s"}, f"levels[{index}]"
        )
        levels.append(
            ElectricThrottleLevel(
                _text(item, "level_id", default=""),
                _number(item, "input_power_w"),
                _number(item, "thrust_n"),
                _number(item, "mass_flow_kg_s"),
            )
        )
    table = ElectricThrottleTable(
        levels=tuple(levels),
        source_id=values["source_id"],
        source_sha256=values["source_sha256"],
    )
    model = PowerLimitedPropulsionModel(
        DiscreteThrottlePath.from_table(
            table,
            objective=values["objective"],
            up_switch_hysteresis_w=values["up_switch_hysteresis_w"],
        ),
        propellant_reserve_kg=values["propellant_reserve_kg"],
        power_margin_w=values["power_margin_w"],
        duty_cycle=values["duty_cycle"],
    )

    def evaluate(level: int, power: float):
        """Evaluate the model at one commanded level and available power."""
        return model.evaluate(
            level_index=level,
            available_power_w=power,
            vehicle_mass_kg=values["vehicle_mass_kg"],
            direction=(0.0, 1.0, 0.0),
            additional_states={model.propellant_state_name: values["propellant_mass_kg"]},
        )

    level = model.select_level(values["available_power_w"])
    evaluation = evaluate(level, values["available_power_w"])
    powers = _linspace(values["sweep_minimum_power_w"], values["sweep_maximum_power_w"], count)
    sweep: dict[str, list[Any]] = {
        "available_power_w": powers,
        "level_index": [],
        "thrust_n": [],
        "specific_impulse_s": [],
        "consumed_power_w": [],
        "mass_flow_kg_s": [],
        "level_index_increasing": [],
        "level_index_decreasing": [],
    }
    for power in powers:
        selected = model.select_level(power)
        point = evaluate(selected, power)
        sweep["level_index"].append(selected)
        sweep["thrust_n"].append(point.thrust_n)
        sweep["specific_impulse_s"].append(point.specific_impulse_s)
        sweep["consumed_power_w"].append(point.consumed_power_w)
        sweep["mass_flow_kg_s"].append(0.0 - point.mass_derivative_kg_s)
    current = None
    for power in powers:
        current = model.select_level(power, current)
        sweep["level_index_increasing"].append(current)
    current = None
    decreasing = []
    for power in reversed(powers):
        current = model.select_level(power, current)
        decreasing.append(current)
    sweep["level_index_decreasing"] = decreasing[::-1]
    throttle_law = model.throttle_law
    if not isinstance(throttle_law, DiscreteThrottlePath):  # constructed from one above
        raise DomainError("Electric throttle sweep requires a discrete throttle path.")
    path_levels = throttle_law.levels
    return {
        "model": "power-limited electric propulsion, discrete throttle path",
        "inputs": {**values, "levels": [asdict(item) for item in levels]},
        "summary": {
            "level_index": level,
            "level_id": evaluation.level_id,
            "active": evaluation.active,
            "thrust_n": evaluation.thrust_n,
            "specific_impulse_s": evaluation.specific_impulse_s,
            "available_power_w": evaluation.available_power_w,
            "effective_power_w": model.effective_power_w(values["available_power_w"]),
            "consumed_power_w": evaluation.consumed_power_w,
            "power_deficit_w": evaluation.power_deficit_w,
            "propellant_mass_rate_kg_s": 0.0 - evaluation.mass_derivative_kg_s,
            "acceleration_m_s2": list(evaluation.acceleration_m_s2),
            "inhibited_by": list(evaluation.inhibited_by),
            "model_hash": model.model_hash,
            "table_hash": table.table_hash,
        },
        "levels": [
            {
                "level_index": index,
                "level_id": item.level_id,
                "input_power_w": item.input_power_w,
                "thrust_n": item.thrust_n,
                "mass_flow_kg_s": item.mass_flow_kg_s,
                "specific_impulse_s": item.specific_impulse_s,
                "total_efficiency": item.total_efficiency,
            }
            for index, item in enumerate(path_levels)
        ],
        "switching_surfaces": [
            {
                "event_id": surface.event_id,
                "crossing": surface.crossing,
                "threshold_w": surface.threshold_w,
                "from_level": surface.from_level,
                "to_level": surface.to_level,
            }
            for surface in model.switching_surfaces(level)
        ],
        "series": {"power_sweep": _series("available_power_w", sweep, MAX_SERIES_POINTS)},
        "reference_ids": ["NASA-20090004685", "OREKIT-13.1.5-PROPULSION-MODEL"],
        "warnings": [
            ("Throttle-table values are caller-supplied; the model is not validated against "
            "thruster hardware."),
        ],
    }


# ---------------------------------------------------------------------------
# Pulse-mode thruster (examples/pulse_mode_trim.py)
# ---------------------------------------------------------------------------

_PULSE_ON_TIMES = (0.010, 0.020, 0.050, 0.100, 0.200)
# Synthetic 1 N tables exactly as built in examples/pulse_mode_trim.py.
_PULSE_STEADY_FLOW = 1.0 / (STANDARD_GRAVITY_M_S2 * 220.0)
_PULSE_IMPULSE_BITS = tuple(1.0 * (t - 0.006) + 0.002 for t in _PULSE_ON_TIMES)
_PULSE_PROPELLANT = tuple(_PULSE_STEADY_FLOW * (t + 0.003) for t in _PULSE_ON_TIMES)
_PULSE_NUMBERS: dict[str, float] = {
    "steady_thrust_n": 1.0,
    "steady_mass_flow_kg_s": _PULSE_STEADY_FLOW,
    "minimum_off_time_s": 0.02,
    "maximum_on_time_s": 2.0,
    "initial_mass_kg": 12.0,
    "propellant_budget_kg": 2.0,
    "delta_v_m_s": 0.05,
    "pulse_off_time_s": 0.05,
}


def plan_pulse_mode_route(payload: dict[str, Any]) -> dict[str, Any]:
    """Plan a quantized pulse train for a Δv request and tabulate the impulse-bit curve.

    Parameters
    ----------
    payload : dict
        Decoded JSON request body; every field is optional (defaults from
        ``examples/pulse_mode_trim.py``):

        on_times_s, impulse_bits_n_s, propellant_per_pulse_kg : list of float
            Characterization table; default the synthetic 1 N table at 10, 20,
            50, 100 and 200 ms.
        steady_thrust_n, steady_mass_flow_kg_s : float
            Defaults 1 N and 1 / (g0 * 220 s) kg/s.
        minimum_off_time_s, maximum_on_time_s : float
            Defaults 0.02 s and 2 s.
        source_id, source_sha256 : str
            Defaults ``"synthetic-1N"`` and 64 ``"c"``.
        initial_mass_kg, propellant_budget_kg, delta_v_m_s : float
            Δv request; defaults 12 kg, 2 kg and 0.05 m/s.
        quantization : {"nearest", "no_overshoot"}, default "nearest"
            Train choice when no integer train meets the request.
        pulse_off_time_s : float, default 0.05
            Off-time added to the on-time to form the train period.
        quantization_requests_m_s : list of float
            Further Δv requests tabulated as plans; default 0.0002, 0.0004,
            0.05 and 1.0 m/s.
        sweep_points : int, default 101
            Geometrically spaced on-times of the impulse-bit curve, 2 to 400.
        max_points : int, default 400
            Point ceiling of the pulse-train history.

    Returns
    -------
    dict
        ``summary`` (minimum impulse bit, Isp at the minimum bit, plan and
        free-space burn metrics of the train), ``quantization`` (one plan per
        request), ``series.impulse_bit_curve`` (impulse, propellant and Isp
        versus on-time) and ``series.pulse_train`` (thrust and mass flow
        versus time).

    Raises
    ------
    InputError
        If a field is unknown, malformed or not finite.
    DomainError
        If the table is invalid, the budget cannot deliver a request, or the
        period is too short.
    ConvergenceError
        If the free-space burn metrics of the train do not converge.

    Notes
    -----
    Calls :class:`~rocket_propulsion.propulsion.burns.PulseModeThruster`
    (``pulse``, ``plan_delta_v``, ``pulse_train_artifact``) and
    :func:`~rocket_propulsion.propulsion.burns.free_space_burn_metrics`.
    Quantization follows from the minimum impulse bit [1]_; short-pulse
    impulse and Isp vary with on-time [2]_; the train Δv uses the ideal rocket
    equation [3]_.

    References
    ----------
    .. [1] NASA-20070032005. J. M. Parker and M. J. Wilson, "The Minimum Impulse
           Thruster", 2005. https://ntrs.nasa.gov/citations/20070032005
    .. [2] NASA-CR-122347. NASA, "Monopropellant Hydrazine Resisto Jet".
           https://ntrs.nasa.gov/citations/19720011120
    .. [3] NASA-ROCKET. NASA Glenn Research Center, "Ideal Rocket Equation",
           https://www1.grc.nasa.gov/beginners-guide-to-aeronautics/ideal-rocket-equation/
    """

    _reject_unknown(
        payload,
        set(_PULSE_NUMBERS)
        | {
            "on_times_s",
            "impulse_bits_n_s",
            "propellant_per_pulse_kg",
            "source_id",
            "source_sha256",
            "quantization",
            "quantization_requests_m_s",
            "sweep_points",
            "max_points",
        },
    )
    count = _sweep_points(payload)
    limit = _max_points(payload)
    values: dict[str, Any] = {
        key: _number(payload, key, default=default) for key, default in _PULSE_NUMBERS.items()
    }
    values["on_times_s"] = list(_numbers(payload, "on_times_s", default=_PULSE_ON_TIMES))
    values["impulse_bits_n_s"] = list(
        _numbers(payload, "impulse_bits_n_s", default=_PULSE_IMPULSE_BITS)
    )
    values["propellant_per_pulse_kg"] = list(
        _numbers(payload, "propellant_per_pulse_kg", default=_PULSE_PROPELLANT)
    )
    values["source_id"] = _text(payload, "source_id", default="synthetic-1N")
    values["source_sha256"] = _text(payload, "source_sha256", default="c" * 64)
    values["quantization"] = _text(payload, "quantization", default="nearest")
    values["quantization_requests_m_s"] = list(
        _numbers(
            payload,
            "quantization_requests_m_s",
            default=(0.0002, 0.0004, 0.05, 1.0),
            maximum_length=16,
        )
    )
    thruster = PulseModeThruster(
        on_times_s=tuple(values["on_times_s"]),
        impulse_bits_n_s=tuple(values["impulse_bits_n_s"]),
        propellant_per_pulse_kg=tuple(values["propellant_per_pulse_kg"]),
        steady_thrust_n=values["steady_thrust_n"],
        steady_mass_flow_kg_s=values["steady_mass_flow_kg_s"],
        minimum_off_time_s=values["minimum_off_time_s"],
        source_id=values["source_id"],
        source_sha256=values["source_sha256"],
        maximum_on_time_s=values["maximum_on_time_s"],
    )

    def plan_for(request: float):
        """Plan the fewest equal pulses for one delta-v request."""
        return thruster.plan_delta_v(
            request,
            values["initial_mass_kg"],
            propellant_budget_kg=values["propellant_budget_kg"],
            quantization=values["quantization"],
        )

    def plan_row(plan) -> dict[str, Any]:
        """Flatten a pulse plan into JSON scalars."""
        return {
            "requested_m_s": plan.requested,
            "delivered_m_s": plan.delivered,
            "residual_m_s": plan.residual,
            "quantized": plan.quantized,
            "pulse_count": plan.pulse_count,
            "on_time_s": None if plan.pulse is None else plan.pulse.on_time_s,
            "delivered_impulse_n_s": plan.delivered_impulse_n_s,
            "propellant_kg": plan.propellant_kg,
            "final_mass_kg": plan.final_mass_kg,
            "effective_specific_impulse_s": plan.effective_specific_impulse_s,
        }

    plan = plan_for(values["delta_v_m_s"])
    train: dict[str, Any] | None = None
    train_series: dict[str, Any] | None = None
    if plan.pulse is not None and plan.pulse_count > 0:
        artifact = thruster.pulse_train_artifact(
            plan, period_s=plan.pulse.on_time_s + values["pulse_off_time_s"]
        )
        metrics = free_space_burn_metrics(artifact, values["initial_mass_kg"])
        train = {
            **_artifact_summary(artifact),
            "free_space_delta_v_m_s": metrics.delta_v_m_s,
            "free_space_final_mass_kg": metrics.final_mass_kg,
            "delta_v_error_estimate_m_s": metrics.delta_v_error_estimate_m_s,
        }
        train_series = _artifact_series(artifact, limit)
    shortest = thruster.minimum_on_time_s
    longest = thruster.longest_on_time_s
    curve: dict[str, list[Any]] = {
        "on_time_s": [],
        "impulse_n_s": [],
        "propellant_kg": [],
        "specific_impulse_s": [],
    }
    for index in range(count):
        on_time = (
            longest if index == count - 1 else shortest * (longest / shortest) ** (index / (count - 1))
        )
        bit = thruster.pulse(on_time)
        curve["on_time_s"].append(bit.on_time_s)
        curve["impulse_n_s"].append(bit.impulse_n_s)
        curve["propellant_kg"].append(bit.propellant_kg)
        curve["specific_impulse_s"].append(bit.specific_impulse_s)
    series = {"impulse_bit_curve": _series("on_time_s", curve, MAX_SERIES_POINTS)}
    if train_series is not None:
        series["pulse_train"] = train_series
    return {
        "model": "pulse-mode thruster with minimum-impulse-bit quantization",
        "inputs": values,
        "summary": {
            "minimum_impulse_bit_n_s": thruster.minimum_impulse_bit_n_s,
            "minimum_on_time_s": shortest,
            "longest_on_time_s": longest,
            "specific_impulse_at_minimum_bit_s": thruster.pulse(shortest).specific_impulse_s,
            "plan": plan_row(plan),
            "pulse_train": train,
            "model_hash": thruster.model_hash,
        },
        "quantization": [plan_row(plan_for(request)) for request in values["quantization_requests_m_s"]],
        "series": series,
        "reference_ids": list(thruster.reference_ids),
        "warnings": [
            ("The impulse-bit and propellant tables are caller-supplied; the default table is "
            "synthetic."),
            ("Pulse-train segments are rectangles over the electrical on-time; rise and tail-off "
            "are not represented."),
        ],
    }


# ---------------------------------------------------------------------------
# Seeded blowdown uncertainty ensemble (examples/blowdown_uncertainty.py)
# ---------------------------------------------------------------------------

_UNCERTAINTY_MODEL_DEFAULTS: dict[str, float | None] = {
    **_BLOWDOWN_DEFAULTS,
    "maximum_supply_pressure_pa": 2.2e6,
}
_UNCERTAINTY_PARAMETERS = (
    {
        "name": "initial_pressure_pa",
        "lower": 1.9e6,
        "upper": 2.1e6,
        "distribution": "truncated-normal",
        "uncertainty_class": "aleatory",
        "standard_deviation": 40_000.0,
        "source": "illustrative tank-pressure acceptance band",
        "rationale": "demonstrate bounded realization-to-realization variation",
    },
    {
        "name": "rated_delivered_thrust_n",
        "lower": 97.0,
        "upper": 103.0,
        "distribution": "triangular",
        "uncertainty_class": "epistemic",
        "source": "illustrative hot-fire calibration band",
        "rationale": "demonstrate reducible rated-performance uncertainty",
    },
)
_UNCERTAINTY_METRICS = (
    "duration_s",
    "delivered_impulse_n_s",
    "consumed_propellant_kg",
    "thrust_centroid_time_s",
    "final_supply_pressure_pa",
)


def _rank_sensitivities(parameters: Sequence[Any], realizations: dict[str, list[Any]]) -> list[dict[str, Any]]:
    """Spearman, PRCC and SRRC of every metric to every sampled input [1]_.

    Computed from the same seeded sample as the Pearson screening records;
    a metric whose values are not all finite, or a sample too small for
    PRCC, is reported with the reason instead.

    References
    ----------
    .. [1] SAND2001-0417. Helton and Davis, Sect. 6, https://doi.org/10.2172/806696
    """

    inputs = {parameter.name: [float(value) for value in realizations[f"input_{parameter.name}"]]
              for parameter in parameters}
    records: list[dict[str, Any]] = []
    for metric in _UNCERTAINTY_METRICS:
        values = realizations[metric]
        try:
            if any(value is None for value in values):
                raise DomainError(f"{metric} is undefined in some realizations.")
            records.append(asdict(rank_sensitivity(inputs, [float(value) for value in values], output_name=metric)))
        except DomainError as error:
            records.append({"output": metric, "unavailable": str(error)})
    return records


def simulate_uncertainty_route(payload: dict[str, Any]) -> dict[str, Any]:
    """Propagate sourced input uncertainty through seeded blowdown realizations.

    Parameters
    ----------
    payload : dict
        Decoded JSON request body; every field is optional (defaults from
        ``examples/blowdown_uncertainty.py``). It accepts every model field of
        ``/api/v1/transient/blowdown`` (same defaults, except
        ``maximum_supply_pressure_pa`` = 2.2 MPa) and:

        parameters : list of object
            Uncertain inputs ``{name, lower, upper, distribution,
            uncertainty_class, standard_deviation?, source, rationale,
            nominal?}``; ``nominal`` defaults to the resolved model input.
            Default: initial pressure 1.9-2.1 MPa truncated normal
            (sigma 40 kPa, aleatory) and rated thrust 97-103 N triangular
            (epistemic).
        sample_count : int, default 16
            Realizations, 4 to 256 (the example uses 32; 16 keeps the default
            request near one second).
        seed : int, default 42
            Sampling seed.
        method : {"latin-hypercube", "monte-carlo"}, default "latin-hypercube"
        convergence_tolerance : float, default 0.05
            Half/full sample stability tolerance.
        maximum_mass_step_kg : float, default 0.5
            Mesh step of each realization.
        integration_relative_tolerance : float, default 1e-3
            Mesh tolerance of each realization (the example uses 1e-6, about
            1.4 s per realization); at least 1e-6.
        histogram_bins : int, default 10
            Equal-width bins per metric histogram, 1 to 50.

    Returns
    -------
    dict
        ``statistics`` (mean, standard deviation, P05/P50/P95, extremes and
        mean interval per metric), ``sensitivity`` (Pearson screening),
        ``convergence``, ``histograms`` per metric, ``series.realizations``
        (sampled inputs and metric values per realization), ``reference_ids``
        and ``warnings``.

    Raises
    ------
    InputError
        If a field is unknown, malformed, not finite or exceeds a work limit.
    DomainError
        If a parameter, bound or control is invalid, or a realization fails.

    Notes
    -----
    Calls :func:`~rocket_propulsion.propulsion.burns.simulate_blowdown_uncertainty`.
    Uncertain inputs are classified as aleatory or epistemic [1]_; Latin
    hypercube sampling follows [2]_; the outcome summaries are engineering
    estimates within a probabilistic risk framework [3]_, not flight
    certification.

    References
    ----------
    .. [1] NASA-STD-7009B, "Standard for Models and Simulations".
           https://standards.nasa.gov/standard/NASA/NASA-STD-7009
    .. [2] MCKAY-1979. M. D. McKay, R. J. Beckman and W. J. Conover, "A
           Comparison of Three Methods for Selecting Values of Input Variables in
           the Analysis of Output from a Computer Code", Technometrics 21(2),
           239-245, 1979.
    .. [3] NASA-SP-2011-3421. NASA, "Probabilistic Risk Assessment Procedures
           Guide for NASA Managers and Practitioners".
           https://ntrs.nasa.gov/citations/20120001369
    """

    _reject_unknown(
        payload,
        set(_BLOWDOWN_MODEL_FIELDS)
        | {
            "parameters",
            "sample_count",
            "seed",
            "method",
            "convergence_tolerance",
            "maximum_mass_step_kg",
            "integration_relative_tolerance",
            "histogram_bins",
        },
    )
    values, rated, tank, law = _blowdown_records(payload, _UNCERTAINTY_MODEL_DEFAULTS)
    values["sample_count"] = _bounded_integer(
        payload, "sample_count", default=16, minimum=4, maximum=MAX_UNCERTAINTY_SAMPLES
    )
    values["seed"] = _integer(payload, "seed", default=42)
    values["method"] = _text(payload, "method", default="latin-hypercube")
    values["convergence_tolerance"] = _number(payload, "convergence_tolerance", default=0.05)
    values["maximum_mass_step_kg"] = _mass_step(payload, values, 0.5)
    values["integration_relative_tolerance"] = _tolerance(
        payload, "integration_relative_tolerance", default=1e-3, minimum=1e-6
    )
    values["histogram_bins"] = _bounded_integer(
        payload, "histogram_bins", default=10, minimum=1, maximum=50
    )
    nominal_values = {
        "rated_delivered_thrust_n": values["delivered_thrust_n"],
        "rated_system_specific_impulse_s": values["system_specific_impulse_s"],
        **values,
    }
    parameters = []
    for index, item in enumerate(_items(payload, "parameters", default=_UNCERTAINTY_PARAMETERS)):
        label = f"parameters[{index}]"
        _reject_unknown(
            item,
            {
                "name",
                "nominal",
                "lower",
                "upper",
                "distribution",
                "uncertainty_class",
                "standard_deviation",
                "source",
                "rationale",
            },
            label,
        )
        name = _text(item, "name", default="")
        if "nominal" in item:
            nominal = _number(item, "nominal")
        elif isinstance(nominal_values.get(name), float):
            nominal = nominal_values[name]
        else:
            raise InputError(
                f"Field {label}.nominal is required for parameter {name!r}.",
                code="missing_field",
                field=f"{label}.nominal",
            )
        parameters.append(
            UncertainParameter(
                name=name,
                nominal=nominal,
                lower=_number(item, "lower"),
                upper=_number(item, "upper"),
                distribution=_enum_value(
                    DistributionKind, item.get("distribution"), f"{label}.distribution"
                ),
                uncertainty_class=_enum_value(
                    UncertaintyClass, item.get("uncertainty_class"), f"{label}.uncertainty_class"
                ),
                source=_text(item, "source", default=""),
                rationale=_text(item, "rationale", default=""),
                standard_deviation=_optional_number_or(item, "standard_deviation", None),
            )
        )
    result = simulate_blowdown_uncertainty(
        tank,
        rated,
        law,
        tuple(parameters),
        reserve_mass_kg=values["reserve_mass_kg"],
        cant_efficiency=values["cant_efficiency"],
        sample_count=values["sample_count"],
        seed=values["seed"],
        method=values["method"],
        convergence_tolerance=values["convergence_tolerance"],
        maximum_mass_step_kg=values["maximum_mass_step_kg"],
        integration_relative_tolerance=values["integration_relative_tolerance"],
    )
    realizations: dict[str, list[Any]] = {"index": [sample.index for sample in result.samples]}
    for parameter in result.parameters:
        realizations[f"input_{parameter.name}"] = [
            dict(sample.inputs)[parameter.name] for sample in result.samples
        ]
    for metric in _UNCERTAINTY_METRICS:
        realizations[metric] = [getattr(sample, metric) for sample in result.samples]
    realizations["cutoff_reason"] = [sample.cutoff_reason for sample in result.samples]
    return {
        "model": "seeded pressure-fed blowdown uncertainty ensemble",
        "inputs": {
            **values,
            "parameters": [
                {
                    **asdict(parameter),
                    "distribution": parameter.distribution.value,
                    "uncertainty_class": parameter.uncertainty_class.value,
                }
                for parameter in result.parameters
            ],
        },
        "summary": {
            "method": result.method,
            "seed": result.seed,
            "sample_count": result.sample_count,
            "aleatory_parameter_count": result.aleatory_parameter_count,
            "epistemic_parameter_count": result.epistemic_parameter_count,
            "sampling_stable": result.convergence.converged,
            "sampling_numerics": result.sampling_numerics,
            "summary_numerics": result.summary_numerics,
        },
        "statistics": [asdict(item) for item in result.statistics],
        "sensitivity": [asdict(item) for item in result.sensitivity],
        "global_sensitivity": _rank_sensitivities(result.parameters, realizations),
        "convergence": asdict(result.convergence),
        "histograms": {
            metric: _histogram(realizations[metric], values["histogram_bins"])
            for metric in _UNCERTAINTY_METRICS
        },
        "series": {"realizations": _series("index", realizations, MAX_SERIES_POINTS)},
        "reference_ids": list(result.reference_ids),
        "warnings": list(result.warnings),
    }


# ---------------------------------------------------------------------------
# Discrete propulsion scenarios (examples/propulsion_scenarios.py)
# ---------------------------------------------------------------------------

_SCENARIO_SOURCE = "illustrative event-tree fixture"
_SCENARIO_EVENT_SOURCE = "illustrative reliability fixture"
_SCENARIOS = (
    {
        "scenario_id": "nominal",
        "probability": 0.90,
        "events": [],
        "source": _SCENARIO_SOURCE,
        "rationale": "both engines complete the command",
    },
    {
        "scenario_id": "engine-a-out",
        "probability": 0.08,
        "events": [
            {
                "event_type": "engine-disabled",
                "engine_ids": ["engine-a"],
                "value": None,
                "source": _SCENARIO_EVENT_SOURCE,
                "rationale": "engine-a fails to ignite",
            }
        ],
        "source": _SCENARIO_SOURCE,
        "rationale": "single-engine contingency",
    },
    {
        "scenario_id": "common-cause-no-ignition",
        "probability": 0.02,
        "events": [
            {
                "event_type": "engine-disabled",
                "engine_ids": ["engine-a", "engine-b"],
                "value": None,
                "source": _SCENARIO_EVENT_SOURCE,
                "rationale": "shared initiation system unavailable",
            }
        ],
        "source": _SCENARIO_SOURCE,
        "rationale": "dependent loss of both engines",
    },
)


def _strings(item: dict[str, Any], key: str, label: str) -> tuple[str, ...]:
    """Read a non-empty list of non-blank strings."""

    raw = item.get(key)
    if (
        not isinstance(raw, list)
        or not raw
        or len(raw) > _MAX_LIST_ITEMS
        or not all(isinstance(value, str) and value.strip() for value in raw)
    ):
        raise InputError(f"Field {label} must be a non-empty list of strings.", field=label)
    return tuple(raw)


def build_scenarios_route(payload: dict[str, Any]) -> dict[str, Any]:
    """Build probability-weighted engine-cluster histories for discrete outcomes.

    Parameters
    ----------
    payload : dict
        Decoded JSON request body; every field is optional (defaults from
        ``examples/propulsion_scenarios.py``):

        delivered_thrust_n, system_specific_impulse_s, chamber_pressure_pa : float
            Per-engine monopropellant point; defaults 100 N, 220 s, 1.0 MPa.
        engine_ids : list of str, default ["engine-a", "engine-b"]
            Cluster members, each commanded by the same schedule.
        burn_duration_s, throttle : float
            One steady command segment; defaults 10 s and 1.
        scenarios : list of object
            Exhaustive, mutually exclusive scenarios ``{scenario_id,
            probability, events, source, rationale}``; each event is
            ``{event_type, engine_ids, value, source, rationale}`` with
            ``event_type`` one of engine-disabled, ignition-delay,
            early-cutoff, thrust-scale, specific-impulse-scale,
            response-time-scale. Default: nominal 0.90, engine-a out 0.08,
            common-cause loss of both engines 0.02 (illustrative, not failure
            rates).
        objective : str
            Default ``"illustrate propulsion consequence weighting"``.
        maximum_step_s : float, default 0.25
            Sampling step of each cluster profile.
        max_points : int, default 400
            Point ceiling of each scenario history.

    Returns
    -------
    dict
        ``summary`` (zero-thrust probability), ``statistics`` (probability
        weighted expected value, P05/P50/P95 and extremes per metric),
        ``scenarios`` (probability and consequences per scenario),
        ``series.consequences`` (one bar per scenario) and
        ``scenario_histories`` (thrust and mass flow versus time per
        scenario).

    Raises
    ------
    InputError
        If a field is unknown, malformed or not finite.
    DomainError
        If probabilities are not exhaustive, an event names an unknown
        engine, or a value violates its event rule.

    Notes
    -----
    Calls :func:`~rocket_propulsion.propulsion.burns.build_propulsion_scenario_ensemble`.
    Discrete outcomes are weighted as in probabilistic risk assessment [1]_;
    engine-out, start and cutoff outcomes are reliability drivers [2]_;
    several engines in one event express a declared common-cause failure
    [3]_. No failure rate or independence is assumed.

    References
    ----------
    .. [1] NASA-SP-2011-3421. NASA, "Probabilistic Risk Assessment Procedures
           Guide for NASA Managers and Practitioners".
           https://ntrs.nasa.gov/citations/20120001369
    .. [2] NASA-20050207429. NASA, "Key Reliability Drivers of Liquid Propulsion
           Engines and a Reliability Model for Sensitivity Analysis".
           https://ntrs.nasa.gov/citations/20050207429
    .. [3] NASA-20160007073. NASA, "Common Cause Failure Modeling in Space Launch
           Vehicles". https://ntrs.nasa.gov/citations/20160007073
    """

    _reject_unknown(
        payload,
        {
            "delivered_thrust_n",
            "system_specific_impulse_s",
            "chamber_pressure_pa",
            "engine_ids",
            "burn_duration_s",
            "throttle",
            "scenarios",
            "objective",
            "maximum_step_s",
            "max_points",
        },
    )
    limit = _max_points(payload)
    values: dict[str, Any] = {
        "delivered_thrust_n": _number(payload, "delivered_thrust_n", default=100.0),
        "system_specific_impulse_s": _number(payload, "system_specific_impulse_s", default=220.0),
        "chamber_pressure_pa": _number(payload, "chamber_pressure_pa", default=1.0e6),
        "engine_ids": list(
            _strings({"engine_ids": ["engine-a", "engine-b"], **payload}, "engine_ids", "engine_ids")
        ),
        "burn_duration_s": _number(payload, "burn_duration_s", default=10.0),
        "throttle": _number(payload, "throttle", default=1.0),
        "objective": _text(
            payload, "objective", default="illustrate propulsion consequence weighting"
        ),
        "maximum_step_s": _number(payload, "maximum_step_s", default=0.25),
    }
    point = ConstantPerformanceProvider.monopropellant(
        delivered_thrust_n=values["delivered_thrust_n"],
        system_specific_impulse_s=values["system_specific_impulse_s"],
        chamber_pressure_pa=values["chamber_pressure_pa"],
        source="illustrative engine point",
    ).operating_point()
    members = tuple(EngineClusterMember(engine_id, point) for engine_id in values["engine_ids"])
    schedule = ThrottleSchedule(
        (ThrottleSegment(values["burn_duration_s"], values["throttle"], values["throttle"], "steady"),)
    )
    commands = tuple(EngineTransientCommand(engine_id, schedule) for engine_id in values["engine_ids"])
    scenario_items = _items(payload, "scenarios", default=_SCENARIOS)
    scenarios = []
    for index, item in enumerate(scenario_items):
        label = f"scenarios[{index}]"
        _reject_unknown(
            item, {"scenario_id", "probability", "events", "source", "rationale"}, label
        )
        raw_events = item.get("events", [])
        if not isinstance(raw_events, list) or len(raw_events) > _MAX_LIST_ITEMS:
            raise InputError(f"Field {label}.events must be a list.", field=f"{label}.events")
        events = []
        for event_index, raw_event in enumerate(raw_events):
            event_label = f"{label}.events[{event_index}]"
            event = _object(raw_event, event_label)
            _reject_unknown(
                event, {"event_type", "engine_ids", "value", "source", "rationale"}, event_label
            )
            events.append(
                ScenarioEvent(
                    _enum_value(
                        ScenarioEventType, event.get("event_type"), f"{event_label}.event_type"
                    ),
                    _strings(event, "engine_ids", f"{event_label}.engine_ids"),
                    _optional_number_or(event, "value", None),
                    _text(event, "source", default=""),
                    _text(event, "rationale", default=""),
                )
            )
        scenarios.append(
            PropulsionScenario(
                _text(item, "scenario_id", default=""),
                _number(item, "probability"),
                tuple(events),
                _text(item, "source", default=""),
                _text(item, "rationale", default=""),
            )
        )
    result = build_propulsion_scenario_ensemble(
        members,
        commands,
        PropulsionScenarioSet(scenarios=tuple(scenarios), objective=values["objective"]),
        maximum_step_s=values["maximum_step_s"],
    )
    rows = []
    for realization in result.realizations:
        artifact = realization.artifact
        rows.append(
            {
                "scenario_id": realization.scenario_id,
                "probability": realization.probability,
                "duration_s": artifact.duration_s,
                "delivered_impulse_n_s": artifact.delivered_total_impulse_n_s,
                "axial_impulse_n_s": artifact.axial_total_impulse_n_s,
                "consumed_propellant_kg": artifact.consumed_propellant_kg,
                "thrust_centroid_time_s": artifact.thrust_centroid_time_s,
                "artifact_hash": artifact.artifact_hash,
            }
        )
    consequence_columns = (
        "scenario_id",
        "probability",
        "delivered_impulse_n_s",
        "consumed_propellant_kg",
        "duration_s",
    )
    return {
        "model": "probability-weighted discrete engine-cluster scenarios",
        "inputs": {
            **values,
            "scenarios": [
                {
                    "scenario_id": scenario.scenario_id,
                    "probability": scenario.probability,
                    "events": [
                        {
                            "event_type": event.event_type.value,
                            "engine_ids": list(event.engine_ids),
                            "value": event.value,
                            "source": event.source,
                            "rationale": event.rationale,
                        }
                        for event in scenario.events
                    ],
                    "source": scenario.source,
                    "rationale": scenario.rationale,
                }
                for scenario in result.scenarios
            ],
        },
        "summary": {
            "objective": result.objective,
            "scenario_count": len(result.realizations),
            "zero_thrust_probability": result.zero_thrust_probability,
        },
        "statistics": [asdict(item) for item in result.statistics],
        "scenarios": rows,
        "series": {
            "consequences": _series(
                "scenario_id",
                {name: [row[name] for row in rows] for name in consequence_columns},
                MAX_SERIES_POINTS,
            )
        },
        "scenario_histories": {
            realization.scenario_id: _artifact_series(realization.artifact, limit)
            for realization in result.realizations
        },
        "reference_ids": list(result.reference_ids),
        "warnings": list(result.warnings),
    }


# ---------------------------------------------------------------------------
# Multi-burn timeline with jettison (examples/propulsion_timeline.py)
# ---------------------------------------------------------------------------

_TIMELINE_BURNS = (
    {
        "burn_id": "perigee",
        "ignition_time_s": 10.0,
        "duration_s": 120.0,
        "thrust_n": 30_000.0,
        "oxidizer_mass_flow_kg_s": 7.8,
        "fuel_mass_flow_kg_s": 2.2,
    },
    {
        "burn_id": "apogee",
        "ignition_time_s": 400.0,
        "duration_s": 40.0,
        "thrust_n": 30_000.0,
        "oxidizer_mass_flow_kg_s": 7.8,
        "fuel_mass_flow_kg_s": 2.2,
    },
)
_TIMELINE_JETTISONS = (
    {"event_id": "adapter", "time_s": 200.0, "mass_kg": 150.0, "label": "payload adapter"},
)


def evaluate_timeline_route(payload: dict[str, Any]) -> dict[str, Any]:
    """Place constant-thrust burns and jettisons on one clock and tabulate the dynamics.

    Parameters
    ----------
    payload : dict
        Decoded JSON request body; every field is optional (defaults from
        ``examples/propulsion_timeline.py``):

        burns : list of object
            ``{burn_id, ignition_time_s, duration_s, thrust_n,
            oxidizer_mass_flow_kg_s, fuel_mass_flow_kg_s}``; defaults perigee
            (10 s, 120 s) and apogee (400 s, 40 s), each 30 kN with 7.8 kg/s
            oxidizer and 2.2 kg/s fuel (c = F / mdot = 3000 m/s).
        jettisons : list of object
            ``{event_id, time_s, mass_kg, label}``; default a 150 kg payload
            adapter at 200 s. An empty list is allowed.
        initial_vehicle_mass_kg, oxidizer_mass_kg, fuel_mass_kg : float
            Defaults 6000 kg, 2400 kg and 700 kg.
        end_time_s : float, default 500
            End of the tabulated window.
        allow_overlapping_burns : bool, default false
        max_points : int, default 400
            Point ceiling of the timeline history.

    Returns
    -------
    dict
        ``summary`` (total impulse, propellant, jettisoned mass and the
        resulting mass ledger, timeline hash), ``burns`` (per-burn artifact
        integrals), ``events`` (burn start, boundary, stop and jettison
        roots) and ``series.history`` (summed axial thrust, mass flow and
        active-burn count versus time, with left and right limits at every
        timed event).

    Raises
    ------
    InputError
        If a field is unknown, malformed or not finite.
    DomainError
        If burns overlap without permission, identifiers collide, or a record
        is outside its domain.

    Notes
    -----
    Calls :class:`~rocket_propulsion.propulsion.burns.PropulsionTimeline`
    ``evaluate``, ``events`` and ``jettisoned_mass_kg``. The timeline computes
    no orbit and integrates no state: the consumer propagates mass and
    velocity, so acceleration and Δv histories are not returned. The final
    mass is the ledger initial mass minus every burn's consumed propellant and
    every jettison. Event-driven propulsion contract [1]_; jettisons are
    instantaneous staging mass drops [2]_; tabulated burn profiles follow the
    transient-modelling context of [3]_.

    References
    ----------
    .. [1] OREKIT-13.1.5-PROPULSION-MODEL. Orekit 13.1.5, "PropulsionModel /
           EventDetector" API.
    .. [2] SUTTON-BIBLARZ-2017. G. P. Sutton and O. Biblarz, "Rocket Propulsion
           Elements", 9th ed., Wiley, 2017.
    .. [3] NASA-20040000363. NASA, "Transient Mathematical Modeling for Liquid
           Rocket Engine Systems". https://ntrs.nasa.gov/citations/20040000363
    """

    _reject_unknown(
        payload,
        {
            "burns",
            "jettisons",
            "initial_vehicle_mass_kg",
            "oxidizer_mass_kg",
            "fuel_mass_kg",
            "end_time_s",
            "allow_overlapping_burns",
            "max_points",
        },
    )
    limit = _max_points(payload)
    values: dict[str, Any] = {
        "initial_vehicle_mass_kg": _number(payload, "initial_vehicle_mass_kg", default=6000.0),
        "oxidizer_mass_kg": _number(payload, "oxidizer_mass_kg", default=2400.0),
        "fuel_mass_kg": _number(payload, "fuel_mass_kg", default=700.0),
        "end_time_s": _number(payload, "end_time_s", default=500.0),
        "allow_overlapping_burns": _boolean(payload, "allow_overlapping_burns", default=False),
    }
    if values["end_time_s"] <= 0.0:
        raise InputError("Field end_time_s must be positive.", field="end_time_s")
    burn_fields = (
        "ignition_time_s",
        "duration_s",
        "thrust_n",
        "oxidizer_mass_flow_kg_s",
        "fuel_mass_flow_kg_s",
    )
    burns = []
    burn_inputs = []
    for index, item in enumerate(_items(payload, "burns", default=_TIMELINE_BURNS)):
        label = f"burns[{index}]"
        _reject_unknown(item, {"burn_id", *burn_fields}, label)
        burn_id = _text(item, "burn_id", default="")
        numbers = {key: _number(item, key) for key in burn_fields}
        burn: dict[str, str | float] = {"burn_id": burn_id, **numbers}
        burn_inputs.append(burn)
        streams = (("ox", numbers["oxidizer_mass_flow_kg_s"]), ("fu", numbers["fuel_mass_flow_kg_s"]))
        flow = numbers["oxidizer_mass_flow_kg_s"] + numbers["fuel_mass_flow_kg_s"]
        thrust = numbers["thrust_n"]
        artifact = TabulatedBurnArtifact(
            segments=(
                TabulatedBurnSegment(
                    0.0, numbers["duration_s"], thrust, thrust, thrust, thrust,
                    flow, flow, streams, streams, "steady",
                ),
            ),
            source_id=f"{burn_id}-burn",
            source_sha256="0" * 64,
            reference_ids=("NASA-20040000363",),
        )
        bridge = PropagatorPropulsionBridge(
            artifact,
            stream_bindings=(StreamStateBinding("ox", "m_ox"), StreamStateBinding("fu", "m_fu")),
        )
        burns.append(ScheduledBurn(burn_id, bridge, numbers["ignition_time_s"]))
    jettison_items = payload.get("jettisons", list(_TIMELINE_JETTISONS))
    if not isinstance(jettison_items, list) or len(jettison_items) > _MAX_LIST_ITEMS:
        raise InputError("Field jettisons must be a list.", field="jettisons")
    jettisons = []
    for index, raw in enumerate(jettison_items):
        label = f"jettisons[{index}]"
        item = _object(raw, label)
        _reject_unknown(item, {"event_id", "time_s", "mass_kg", "label"}, label)
        jettisons.append(
            MassJettison(
                _text(item, "event_id", default=""),
                _number(item, "time_s"),
                _number(item, "mass_kg"),
                _text(item, "label", default="jettison"),
            )
        )
    timeline = PropulsionTimeline(
        burns=tuple(burns),
        jettisons=tuple(jettisons),
        allow_overlapping_burns=values["allow_overlapping_burns"],
    )
    events = timeline.events()
    timed = sorted(
        {
            event.time_s
            for event in events
            if event.time_s is not None and 0.0 <= event.time_s <= values["end_time_s"]
        }
    )
    uniform_count = max(2, MAX_SERIES_POINTS - 2 * len(timed))
    samples = sorted(
        {(time, "right") for time in _linspace(0.0, values["end_time_s"], uniform_count)}
        | {(time, side) for time in timed for side in ("left", "right")}
    )
    directions = {burn.burn_id: (1.0, 0.0, 0.0) for burn in burns}
    states = {"m_ox": values["oxidizer_mass_kg"], "m_fu": values["fuel_mass_kg"]}
    history: dict[str, list[Any]] = {
        "time_s": [],
        "axial_thrust_n": [],
        "mass_flow_kg_s": [],
        "active_burn_count": [],
        "jettisoned_mass_kg": [],
    }
    for time, side in samples:
        jettisoned = timeline.jettisoned_mass_kg(time, include_at_time=side == "right")
        evaluation = timeline.evaluate(
            time_s=time,
            vehicle_mass_kg=values["initial_vehicle_mass_kg"] - jettisoned,
            directions=directions,
            additional_states=states,
            side=side,
            apply_inventory_limits=False,
        )
        history["time_s"].append(time)
        history["axial_thrust_n"].append(
            sum(
                (
                    burn_evaluation.axial_thrust_n
                    for _, burn_evaluation in evaluation.burn_evaluations
                    if burn_evaluation.active
                ),
                0.0,
            )
        )
        history["mass_flow_kg_s"].append(0.0 - evaluation.mass_derivative_kg_s)
        history["active_burn_count"].append(len(evaluation.active_burns))
        history["jettisoned_mass_kg"].append(jettisoned)
    burn_rows = [
        {
            "burn_id": burn.burn_id,
            "ignition_time_s": burn.ignition_time_s,
            "cutoff_time_s": burn.cutoff_time_s,
            **_artifact_summary(burn.bridge.artifact),
        }
        for burn in timeline.burns
    ]
    total_propellant = sum(row["consumed_propellant_kg"] for row in burn_rows)
    total_jettison = sum(jettison.mass_kg for jettison in timeline.jettisons)
    return {
        "model": "multi-burn propulsion timeline with discrete mass jettison",
        "inputs": {
            **values,
            "burns": burn_inputs,
            "jettisons": [asdict(jettison) for jettison in timeline.jettisons],
        },
        "summary": {
            "burn_count": len(timeline.burns),
            "jettison_count": len(timeline.jettisons),
            "total_delivered_impulse_n_s": sum(row["delivered_impulse_n_s"] for row in burn_rows),
            "total_consumed_propellant_kg": total_propellant,
            "total_jettisoned_mass_kg": total_jettison,
            "final_vehicle_mass_kg": values["initial_vehicle_mass_kg"]
            - total_propellant
            - total_jettison,
            "timeline_hash": timeline.timeline_hash,
        },
        "burns": burn_rows,
        "events": [
            {
                "event_id": event.event_id,
                "kind": event.kind.value,
                "owner_id": event.owner_id,
                "time_s": event.time_s,
                "stops_burn": event.stops_burn,
                "jettison_mass_kg": event.jettison_mass_kg,
            }
            for event in events
        ],
        "series": {"history": _series("time_s", history, limit)},
        "reference_ids": ["OREKIT-13.1.5-PROPULSION-MODEL", "SUTTON-BIBLARZ-2017",
                          "NASA-20040000363"],
        "warnings": [
            "The timeline computes no orbit; mass and velocity are integrated by the consumer.",
            "Inventory limits are not applied to the tabulated history (consumer-located roots).",
        ],
    }


# ---------------------------------------------------------------------------
# Hot-fire random-effects calibration (examples/hot_fire_calibration.py)
# ---------------------------------------------------------------------------

_CALIBRATION_SOURCE = "synthetic qualification campaign"
_CALIBRATION_MEASUREMENTS = tuple(
    {
        "engine_id": f"engine-{engine_index + 1}",
        "run_id": f"run-{run_index + 1}",
        "condition_id": "rated-steady-state",
        "value": 1.008 + offset + deviation,
        "standard_uncertainty": 0.001,
        "source": _CALIBRATION_SOURCE,
    }
    for engine_index, offset in enumerate((-0.025, -0.015, -0.005, 0.005, 0.015, 0.025))
    for run_index, deviation in enumerate((-0.004, 0.0, 0.004))
)


def calibrate_hot_fire_route(payload: dict[str, Any]) -> dict[str, Any]:
    """Separate engine, run and measurement variance in repeated hot-fire scales.

    Parameters
    ----------
    payload : dict
        Decoded JSON request body; every field is optional (defaults from
        ``examples/hot_fire_calibration.py``):

        parameter_name : str, default "thrust_scale"
            One of thrust_scale, specific_impulse_scale, time_scale,
            axial_efficiency_scale.
        measurements : list of object
            ``{engine_id, run_id, condition_id, value, standard_uncertainty,
            source}``; default six synthetic engines with three runs each
            (scale 1.008 + engine offset + run deviation, u = 0.001).
        source, rationale : str
            Calibration provenance; defaults ``"synthetic qualification
            campaign"`` and ``"demonstrate engine/run/measurement variance
            separation"``.
        relative_tolerance : float, default 1e-10
        maximum_iterations : int, default 200
            Mandel-Paule solver controls (iterations 1 to 1000).
        coverage_factor : float, default 3.0
            Truncation half-width of the residual parameters.

    Returns
    -------
    dict
        ``summary`` (consensus scale and uncertainty, variance components
        and fractions, calibration hash), ``evidence`` (solver record),
        ``residual_parameters`` (population and same-engine truncated-normal
        scales), ``engines`` (per-engine summaries), ``series.engines`` (mean,
        minimum and maximum scale per engine) and ``series.measurements``.

    Raises
    ------
    InputError
        If a field is unknown, malformed or not finite.
    DomainError
        If the data cannot support the random-effects model (for example one
        engine, one run per engine or mixed conditions).

    Notes
    -----
    Calls :func:`~rocket_propulsion.propulsion.burns.calibrate_hot_fire_scale`:
    a one-way random-effects model [1]_ with declared standard uncertainties
    [2]_ and Mandel-Paule consensus weighting [3]_.

    References
    ----------
    .. [1] NIST-SEMATECH-7.4.4. NIST/SEMATECH e-Handbook of Statistical Methods,
           section 7.4.4, "Variance Components".
           https://www.itl.nist.gov/div898/handbook/prc/section4/prc44.htm
    .. [2] NIST-TN-1297. NIST TN 1297, "Guidelines for Evaluating and Expressing
           the Uncertainty of NIST Measurement Results", 1994.
           https://www.nist.gov/pml/nist-technical-note-1297
    .. [3] PAULE-MANDEL-1982. R. C. Paule and J. Mandel, "Consensus Values and
           Weighting Factors", J. Res. Natl. Bur. Stand. 87(5), 377-385, 1982.
    """

    _reject_unknown(
        payload,
        {
            "parameter_name",
            "measurements",
            "source",
            "rationale",
            "relative_tolerance",
            "maximum_iterations",
            "coverage_factor",
        },
    )
    values: dict[str, Any] = {
        "parameter_name": _text(payload, "parameter_name", default="thrust_scale"),
        "source": _text(payload, "source", default=_CALIBRATION_SOURCE),
        "rationale": _text(
            payload, "rationale", default="demonstrate engine/run/measurement variance separation"
        ),
        "relative_tolerance": _number(payload, "relative_tolerance", default=1e-10),
        "maximum_iterations": _bounded_integer(
            payload, "maximum_iterations", default=200, minimum=1, maximum=1000
        ),
        "coverage_factor": _number(payload, "coverage_factor", default=3.0),
    }
    measurements = []
    for index, item in enumerate(_items(payload, "measurements", default=_CALIBRATION_MEASUREMENTS)):
        label = f"measurements[{index}]"
        _reject_unknown(
            item,
            {"engine_id", "run_id", "condition_id", "value", "standard_uncertainty", "source"},
            label,
        )
        measurements.append(
            HotFireScaleMeasurement(
                engine_id=_text(item, "engine_id", default=""),
                run_id=_text(item, "run_id", default=""),
                condition_id=_text(item, "condition_id", default=""),
                value=_number(item, "value"),
                standard_uncertainty=_number(item, "standard_uncertainty"),
                source=_text(item, "source", default=""),
            )
        )
    result = calibrate_hot_fire_scale(
        values["parameter_name"],
        tuple(measurements),
        source=values["source"],
        rationale=values["rationale"],
        relative_tolerance=values["relative_tolerance"],
        maximum_iterations=values["maximum_iterations"],
    )
    residuals = {}
    for scope in (CalibrationScope.POPULATION, CalibrationScope.SAME_ENGINE):
        parameter = result.residual_uncertain_parameter(
            scope=scope, coverage_factor=values["coverage_factor"]
        )
        residuals[scope.value] = {
            "name": parameter.name,
            "nominal": parameter.nominal,
            "lower": parameter.lower,
            "upper": parameter.upper,
            "standard_deviation": parameter.standard_deviation,
            "distribution": parameter.distribution.value,
            "uncertainty_class": parameter.uncertainty_class.value,
        }
    engines = [asdict(summary) for summary in result.engine_summaries]
    engine_index = {row["engine_id"]: index for index, row in enumerate(engines)}
    return {
        "model": "one-way random-effects hot-fire scale calibration (Mandel-Paule)",
        "inputs": {
            **values,
            "measurements": [asdict(measurement) for measurement in result.measurements],
        },
        "summary": {
            "parameter_name": result.parameter_name,
            "condition_id": result.condition_id,
            "method": result.method,
            "observation_count": result.observation_count,
            "engine_count": result.engine_count,
            "within_degrees_of_freedom": result.within_degrees_of_freedom,
            "consensus_scale": result.consensus_scale,
            "consensus_standard_uncertainty": result.consensus_standard_uncertainty,
            "observed_within_engine_variance": result.observed_within_engine_variance,
            "measurement_within_variance_contribution": (
                result.measurement_within_variance_contribution
            ),
            "within_engine_process_variance": result.within_engine_process_variance,
            "between_engine_process_variance": result.between_engine_process_variance,
            "population_process_variance": result.population_process_variance,
            "within_engine_fraction": result.within_engine_fraction,
            "between_engine_fraction": result.between_engine_fraction,
            "calibration_hash": result.calibration_hash,
        },
        "evidence": asdict(result.evidence),
        "residual_parameters": residuals,
        "engines": engines,
        "series": {
            "engines": _series(
                "engine_id",
                {
                    "engine_id": [row["engine_id"] for row in engines],
                    "mean_scale": [row["mean_scale"] for row in engines],
                    "minimum_scale": [row["minimum_scale"] for row in engines],
                    "maximum_scale": [row["maximum_scale"] for row in engines],
                },
                MAX_SERIES_POINTS,
            ),
            "measurements": _series(
                "engine_index",
                {
                    "engine_index": [
                        engine_index[item.engine_id] for item in result.measurements
                    ],
                    "engine_id": [item.engine_id for item in result.measurements],
                    "run_id": [item.run_id for item in result.measurements],
                    "value": [item.value for item in result.measurements],
                    "standard_uncertainty": [
                        item.standard_uncertainty for item in result.measurements
                    ],
                },
                MAX_SERIES_POINTS,
            ),
        },
        "reference_ids": list(result.reference_ids),
        "warnings": list(result.warnings),
    }


# ---------------------------------------------------------------------------
# SASS spacecraft propulsion sizing and trade (examples/sass_mission.json)
# ---------------------------------------------------------------------------

_SASS_NOTE = "Synthetic teaching scenario, not an ITU assignment or hardware measurement: "
_SASS_THRESHOLD = "analyst-selected preliminary feasibility threshold"
_SASS_VERIFICATION = "preliminary numerical analysis only; assumptions require later validation"


def _assumed(value: float, unit: str, note: str) -> dict[str, Any]:
    """One analyst-assumption scalar in the cited-scalar form of the SASS request."""

    return {
        "value": value,
        "unit": unit,
        "citations": [{"source_id": "NASA-SE-6", "locator": "6.8: analyst assumption"}],
        "status": "assumption",
        "note": _SASS_NOTE + note,
    }


def _sass_item(name: str, delta_v: float, margin: float, count: int, note: str) -> dict[str, Any]:
    """One delta-v budget item of the SASS example request."""

    return {
        "name": name,
        "nominal_delta_v": _assumed(delta_v, "m/s", note),
        "margin_fraction": _assumed(margin, "1", "chosen delta-v contingency"),
        "count": count,
    }


def _sass_candidate(
    identifier: str,
    catalog_key: str,
    other_dry_mass: float,
    installed: int,
    firing: int,
    power: float | None,
) -> dict[str, Any]:
    """One propulsion candidate of the SASS example request."""

    return {
        "identifier": identifier,
        "catalog_key": catalog_key,
        "other_propulsion_dry_mass": _assumed(
            other_dry_mass,
            "kg",
            "chosen tank/feed/pressurant/electronics/mounting allocation; excludes catalogued "
            "thruster mass",
        ),
        "installed_count": installed,
        "firing_count": firing,
        "axial_efficiency": _assumed(
            1.0, "1", "parallel thrust along the selected translation direction"
        ),
        "duty_fraction": _assumed(
            0.5, "1", "50 percent operational availability, not a trajectory solution"
        ),
        "firing_bus_power": None
        if power is None
        else _assumed(
            power,
            "W",
            "complete propulsion bus firing allocation; not a published product power",
        ),
    }


def _sass_requirement(
    identifier: str, description: str, metric: str, value: float, unit: str
) -> dict[str, Any]:
    """One at-most requirement of the SASS example request."""

    return {
        "identifier": identifier,
        "description": description,
        "metric": metric,
        "relation": "at_most",
        "limit": _assumed(value, unit, _SASS_THRESHOLD),
        "verification": _SASS_VERIFICATION,
    }


def _sass_criterion(metric: str, best: float, worst: float, weight: float, unit: str):
    """One weighted trade criterion of the SASS example request."""

    return {
        "metric": metric,
        "best": _assumed(best, unit, "fixed preferred utility anchor"),
        "worst": _assumed(worst, unit, "fixed least-preferred utility anchor"),
        "weight": _assumed(weight, "1", "illustrative preference weight"),
    }


def default_sass_study() -> dict[str, Any]:
    """Return the ``examples/sass_mission.json`` request as a fresh dictionary.

    Returns
    -------
    dict
        The synthetic five-candidate UZM 508E teaching study; every scalar is
        an analyst assumption cited to NASA-SE-6 section 6.8.
    """

    return {
        "schema": "rocket_propulsion_sass_study_v1",
        "title": "UZM 508E provisional propulsion trade - synthetic small spacecraft",
        "bus_dry_mass_excluding_propulsion": _assumed(
            50.0, "kg", "chosen bus and payload dry allocation, excluding propulsion"
        ),
        "residual_fraction": _assumed(
            0.02,
            "1",
            "unusable fraction of loaded propellant; not a prescribed industry margin",
        ),
        "fixed_reserve": _assumed(
            0.2, "kg", "additional retained contingency mass, distinct from unusable fraction"
        ),
        "items": [
            _sass_item(
                "injection_correction", 20.0, 0.1, 1,
                "allocation supplied by the assumed mission design team",
            ),
            _sass_item(
                "collision_avoidance", 0.1, 0.3, 20,
                "20 hypothetical avoidance maneuvers; not an encounter forecast",
            ),
            _sass_item(
                "maintenance", 1.0, 0.2, 3,
                "three hypothetical annual allocations; no atmospheric density model",
            ),
            _sass_item(
                "end_of_life", 50.0, 0.1, 1,
                "placeholder disposal allocation; no regulatory compliance or lifetime claim",
            ),
        ],
        "candidates": [
            _sass_candidate("cold_gas_nominal", "moog-058e151", 4.0, 4, 2, 25.0),
            _sass_candidate("hydrazine_nominal", "ariane-1n", 4.0, 2, 1, 25.0),
            _sass_candidate("bipropellant_nominal", "ariane-s10-13", 6.0, 2, 1, 30.0),
            _sass_candidate("green_map_needed", "ecaps-1n-hpgp", 4.0, 2, 1, None),
            _sass_candidate("ion_map_needed", "busek-bit3-page", 3.0, 1, 1, None),
        ],
        "requirements": [
            _sass_requirement(
                "PREP-MASS", "Wet spacecraft mass allocation", "wet_mass", 70.0, "kg"
            ),
            _sass_requirement(
                "PREP-TIME",
                "Total allocated maneuver elapsed-time ceiling",
                "elapsed_time",
                86400.0,
                "s",
            ),
            _sass_requirement(
                "PREP-POWER", "Firing bus power allocation", "firing_bus_power", 50.0, "W"
            ),
        ],
        "criteria": [
            _sass_criterion("wet_mass", 50.0, 80.0, 0.6, "kg"),
            _sass_criterion("elapsed_time", 1000.0, 60000.0, 0.4, "s"),
        ],
    }


_SASS_OVERRIDES = {
    "bus_dry_mass_kg": "bus_dry_mass_excluding_propulsion",
    "residual_fraction": "residual_fraction",
    "fixed_reserve_kg": "fixed_reserve",
}
_SASS_METRICS = (
    "wet_mass",
    "dry_mass",
    "loaded_propellant",
    "consumed_propellant",
    "elapsed_time",
    "on_time",
    "total_delta_v",
    "firing_bus_power",
)


def _cited_source_ids(value: Any) -> set[str]:
    """Collect every ``source_id`` cited anywhere in a study document."""

    if isinstance(value, dict):
        found = {value["source_id"]} if isinstance(value.get("source_id"), str) else set()
        for item in value.values():
            found |= _cited_source_ids(item)
        return found
    if isinstance(value, list):
        found = set()
        for item in value:
            found |= _cited_source_ids(item)
        return found
    return set()


def size_sass_mission_route(payload: dict[str, Any]) -> dict[str, Any]:
    """Size and trade spacecraft propulsion candidates for a cited mission budget.

    Parameters
    ----------
    payload : dict
        Decoded JSON request body; every field is optional:

        study : object
            A complete ``rocket_propulsion_sass_study_v1`` request (schema of
            ``examples/sass_mission.json``); default that example.
        bus_dry_mass_kg, residual_fraction, fixed_reserve_kg : float
            Optional replacements of the study's bus dry mass (default 50 kg),
            unusable propellant fraction (0.02) and fixed reserve (0.2 kg);
            the replaced scalar keeps its analyst-assumption citation and its
            note records the API override.
        include_document : bool, default false
            Also return the complete hash-bearing study document (several
            hundred kB with every citation).

    Returns
    -------
    dict
        ``summary`` (title, nominal and allocated delta-v, study hash),
        ``candidates`` (sized metrics, eligibility, score, requirement checks
        and unresolved items per candidate), ``series.trade`` (one bar per
        candidate: wet mass, consumed propellant, elapsed time, score) and,
        on request, ``document``.

    Raises
    ------
    InputError
        If a field is unknown, malformed or not finite.
    DomainError
        If the study request is invalid or the analysis fails.

    Notes
    -----
    Calls :func:`~rocket_propulsion.studies.sass.project_from_dict`: the
    inverse ideal rocket equation [1]_ with Isp-based flow identities [2]_,
    analyst-defined margins and interfaces [3]_, and the trade-study structure
    of [4]_. Results are preliminary design screens, not qualification
    evidence.

    References
    ----------
    .. [1] NASA-ROCKET. NASA Glenn Research Center, "Ideal Rocket Equation".
    .. [2] NASA-ISP. NASA Glenn Research Center, "Specific Impulse".
    .. [3] NASA-SE-6. NASA, "Systems Engineering Handbook: Crosscutting Technical
           Management".
    .. [4] NASA-SE-APP. NASA, "Systems Engineering Handbook: Appendix"; Trade
           Study Report.
    """

    _reject_unknown(payload, {"study", "include_document", *_SASS_OVERRIDES})
    include_document = _boolean(payload, "include_document", default=False)
    study_request = (
        _object(payload["study"], "study") if "study" in payload else default_sass_study()
    )
    overrides = {}
    for key, field in _SASS_OVERRIDES.items():
        if key in payload:
            value = _number(payload, key)
            overrides[key] = value
            original = _object(study_request.get(field), f"study.{field}")
            study_request = {
                **study_request,
                field: {
                    **original,
                    "value": value,
                    "note": f"{original.get('note', '')} [value overridden through the HTTP API]",
                },
            }
    study = project_from_dict(study_request)
    document = study.to_dict()
    scores = {score.name: score for score in study.scores}
    rows = []
    candidate_values: list[dict[str, float]] = []
    for item in study.candidates:
        metrics = dict(item.metrics)
        candidate_values.append(
            {name: metrics[name].value for name in _SASS_METRICS if name in metrics}
        )
        score = scores.get(item.candidate.identifier)
        rows.append(
            {
                "identifier": item.candidate.identifier,
                "catalog_key": item.candidate.catalog_key,
                "metrics": {
                    name: {"value": metrics[name].value, "unit": metrics[name].unit}
                    for name in _SASS_METRICS
                    if name in metrics
                },
                "eligibility": None if score is None else score.eligibility,
                "score": None if score is None or score.score is None else score.score.value,
                "missing_metrics": [] if score is None else list(score.missing_metrics),
                "checks": []
                if score is None
                else [
                    {
                        "requirement": check.requirement.identifier,
                        "metric": check.requirement.metric,
                        "status": check.status,
                        "measured": None if check.measured is None else check.measured.value,
                        "margin": None if check.margin is None else check.margin.value,
                    }
                    for check in score.checks
                ],
                "unresolved": list(item.unresolved),
            }
        )

    def metric_column(name: str) -> list[float | None]:
        """Candidate values of one sized metric; null where a candidate is unsized."""
        return [values.get(name) for values in candidate_values]

    result: dict[str, Any] = {
        "model": "cited preliminary spacecraft propulsion sizing and weighted trade",
        "inputs": {"study": study.input_data, "overrides": overrides},
        "summary": {
            "title": study.title,
            "nominal_delta_v_m_s": study.budget.get("nominal_delta_v").value,
            "allocated_delta_v_m_s": study.budget.get("allocated_delta_v").value,
            "candidate_count": len(rows),
            "eligible_count": sum(row["eligibility"] == "eligible" for row in rows),
            "study_hash": document["study_hash"],
            "assurance": document["assurance"],
        },
        "candidates": rows,
        "series": {
            "trade": _series(
                "identifier",
                {
                    "identifier": [row["identifier"] for row in rows],
                    "wet_mass_kg": metric_column("wet_mass"),
                    "consumed_propellant_kg": metric_column("consumed_propellant"),
                    "elapsed_time_s": metric_column("elapsed_time"),
                    "score": [row["score"] for row in rows],
                },
                MAX_SERIES_POINTS,
            )
        },
        "reference_ids": sorted(_cited_source_ids(document)),
        "warnings": [
            ("Preliminary design screening only: no hardware readiness, experimental "
            "validation or qualification evidence."),
            f"Catalogue extrapolation: {document['extrapolation_policy']}.",
        ],
    }
    if include_document:
        result["document"] = document
    return result


# ---------------------------------------------------------------------------
# Registries merged into routes.ROUTES / routes.ROUTE_METADATA.
# ---------------------------------------------------------------------------

TRANSIENT_ROUTES = {
    "/api/v1/transient/blowdown": simulate_blowdown_route,
    "/api/v1/transient/regulated-feed": evaluate_regulated_feed_route,
    "/api/v1/transient/feed-line": evaluate_feed_line_route,
    "/api/v1/transient/electric": evaluate_electric_route,
    "/api/v1/transient/pulse-mode": plan_pulse_mode_route,
    "/api/v1/transient/uncertainty": simulate_uncertainty_route,
    "/api/v1/transient/scenarios": build_scenarios_route,
    "/api/v1/transient/timeline": evaluate_timeline_route,
    "/api/v1/transient/calibration": calibrate_hot_fire_route,
    "/api/v1/sass/mission": size_sass_mission_route,
}

_HISTORY_UNITS = {
    "*_s": "s",
    "*_n": "N",
    "*_n_s": "N s",
    "*_kg": "kg",
    "*_kg_s": "kg/s",
    "*_pa": "Pa",
    "*_m3": "m3",
}

TRANSIENT_ROUTE_METADATA = {
    "/api/v1/transient/blowdown": ModelMetadata(
        model="pressure-fed blowdown with mesh-converged tabulation",
        version="1.0",
        assumptions=(
            "rigid tank with user-declared polytropic ullage",
            "pressure-power-law thrust, Isp and chamber pressure; no extrapolation",
            "no regulator, vapor pressure, dissolved gas or line thermal state",
            "synthetic example inputs; not a qualification dataset",
        ),
        units={**_HISTORY_UNITS, "ullage_volume_m3": "m3"},
        references=("NASA-SP-8112", "NASA-CR-131400", "NASA-CR-122347"),
    ),
    "/api/v1/transient/regulated-feed": ModelMetadata(
        model="regulated single-node feed coupled to a pressure-dependent engine",
        version="1.0",
        assumptions=(
            "ideal pressurant, rigid tank and bottle, first-order regulator",
            "quadratic feed-line loss; compressible regulator flow with choking",
            "state evaluation and sweeps only; no time integration",
            "fixed synthetic performance grid and 500 N profile from the example",
        ),
        units={**_HISTORY_UNITS, "*_k": "K", "*_w": "W", "regulator_command": "1"},
        references=("NASA-SP-8112", "NASA-SP-8080", "NASA-20240003493"),
    ),
    "/api/v1/transient/feed-line": ModelMetadata(
        model="lumped inertance/compliance feed line with regulated feed and engine",
        version="1.0",
        assumptions=(
            "first inertance/compliance mode only; not a water-hammer solver",
            "coefficients require calibration against plumbing or a nodal model",
            "state evaluation and sweeps only; no time integration",
            "fixed synthetic performance grid and 500 N profile from the example",
        ),
        units={
            **_HISTORY_UNITS,
            "*_kg_s2": "kg/s2",
            "*_pa_s": "Pa/s",
            "*_rad_s": "rad/s",
            "damping_ratio": "1",
        },
        references=("NASA-19740028545", "NASA-20220006583", "NASA-20240003493"),
    ),
    "/api/v1/transient/electric": ModelMetadata(
        model="power-limited electric propulsion with a discrete throttle path",
        version="1.0",
        assumptions=(
            "caller-supplied throttle table; not validated against hardware",
            "highest level whose input power fits, with up-switch hysteresis",
            "duty cycle averages on/off firing",
        ),
        units={"*_w": "W", "*_n": "N", "*_kg_s": "kg/s", "*_s": "s", "total_efficiency": "1"},
        references=("NASA-20090004685", "OREKIT-13.1.5-PROPULSION-MODEL"),
    ),
    "/api/v1/transient/pulse-mode": ModelMetadata(
        model="pulse-mode thruster impulse-bit table with minimum-impulse-bit quantization",
        version="1.0",
        assumptions=(
            "caller-supplied impulse-bit and propellant tables with PCHIP interpolation",
            "identical pulses; free-space ideal rocket equation per train",
            "rectangular on-time pulses in the train artifact",
        ),
        units={**_HISTORY_UNITS, "*_m_s": "m/s"},
        references=("NASA-20070032005", "NASA-CR-122347"),
    ),
    "/api/v1/transient/uncertainty": ModelMetadata(
        model="seeded pressure-fed blowdown uncertainty ensemble",
        version="1.0",
        assumptions=(
            "declared bounded marginals classified aleatory or epistemic",
            "every realization is a converged blowdown; failures abort the study",
            "empirical statistics; mean interval without established coverage",
            "Pearson screening sensitivity; default mesh tolerance relaxed to 1e-3",
        ),
        units={**_HISTORY_UNITS, "*_m3": "m3"},
        references=("NASA-STD-7009B", "NASA-SP-2011-3421", "MCKAY-1979"),
    ),
    "/api/v1/transient/scenarios": ModelMetadata(
        model="probability-weighted discrete engine-cluster scenarios",
        version="1.0",
        assumptions=(
            "exhaustive mutually exclusive analyst-declared scenarios",
            "no failure rates or independence are inferred",
            "illustrative example probabilities, not reliability data",
        ),
        units={**_HISTORY_UNITS, "probability": "1"},
        references=("NASA-SP-2011-3421", "NASA-20050207429", "NASA-20160007073"),
    ),
    "/api/v1/transient/timeline": ModelMetadata(
        model="multi-burn propulsion timeline with discrete mass jettison",
        version="1.0",
        assumptions=(
            "constant-thrust tabulated burns on one clock",
            "instantaneous jettison mass drops",
            "no orbit; mass and velocity integration belong to the consumer",
        ),
        units=_HISTORY_UNITS,
        references=("OREKIT-13.1.5-PROPULSION-MODEL", "SUTTON-BIBLARZ-2017", "NASA-20040000363"),
    ),
    "/api/v1/transient/calibration": ModelMetadata(
        model="one-way random-effects hot-fire scale calibration",
        version="1.0",
        assumptions=(
            "engine effect, run-to-run process and declared measurement error kept separate",
            "Mandel-Paule heteroscedastic consensus mean",
            "excludes variance-component sampling uncertainty and model-form discrepancy",
        ),
        units={"*_scale": "1", "*_variance": "1", "*_fraction": "1", "standard_uncertainty": "1"},
        references=("NIST-SEMATECH-7.4.4", "NIST-TN-1297", "PAULE-MANDEL-1982"),
    ),
    "/api/v1/sass/mission": ModelMetadata(
        model="cited preliminary spacecraft propulsion sizing and weighted trade",
        version="1.0",
        assumptions=(
            "inverse ideal rocket equation per ordered budget item",
            "published nominal catalogue points; no extrapolation",
            "hard requirement screens then affine weighted utilities",
            "preliminary design screening only; no qualification evidence",
        ),
        units={"*_kg": "kg", "*_s": "s", "*_m_s": "m/s", "score": "1"},
        references=("NASA-ROCKET", "NASA-ISP", "NASA-SE-6", "NASA-SE-APP"),
    ),
}
