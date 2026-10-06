"""JSON-shaped adapters for the calculation domain."""

from __future__ import annotations

import json
import os
from dataclasses import asdict
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


def _number(payload: dict[str, Any], key: str, *, default: float | None = None) -> float:
    """Read a required finite-convertible number from a request payload."""

    raw_value = payload.get(key, default)
    if raw_value is None:
        raise InputError(f"Missing required field: {key}.", code="missing_field", field=key)
    if isinstance(raw_value, bool):
        raise InputError(f"Field {key} must be a number.", field=key)
    try:
        value = float(raw_value)
    except (TypeError, ValueError) as error:
        raise InputError(f"Field {key} must be a number.", field=key) from error
    if not isfinite(value):
        raise InputError(f"Field {key} must be finite.", field=key)
    return value


def _optional_number(payload: dict[str, Any], key: str) -> float | None:
    """Read an optional numeric field, treating an empty string as absent."""

    raw_value = payload.get(key)
    if raw_value in (None, ""):
        return None
    return _number(payload, key)


def _integer(payload: dict[str, Any], key: str, *, default: int) -> int:
    """Read an integer without silently truncating a fractional request value."""

    value = _number(payload, key, default=float(default))
    if not value.is_integer():
        raise InputError(f"Field {key} must be an integer.", field=key)
    return int(value)


def _boolean(payload: dict[str, Any], key: str, *, default: bool) -> bool:
    """Read an explicit JSON/form boolean without Python truthiness surprises."""

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
    """Resolve an input property to Mach number and return all isentropic ratios."""

    input_kind = str(payload.get("input_kind", "mach"))
    gamma = _number(payload, "gamma", default=1.4)
    mach = mach_from_isentropic(input_kind, _number(payload, "value"), gamma)
    return asdict(calculate_isentropic(mach, gamma))


def calculate_normal_shock_route(payload: dict[str, Any]) -> dict[str, Any]:
    """Return normal-shock properties for an upstream supersonic Mach number."""

    gamma = _number(payload, "gamma", default=1.4)
    return asdict(calculate_normal_shock(_number(payload, "upstream_mach"), gamma))


def calculate_oblique_shock_route(payload: dict[str, Any]) -> dict[str, Any]:
    """Solve weak or strong attached oblique-shock relations."""

    return asdict(
        calculate_oblique_shock(
            upstream_mach=_number(payload, "upstream_mach"),
            turning_angle_deg=_number(payload, "turning_angle_deg"),
            gamma=_number(payload, "gamma", default=1.4),
            branch=str(payload.get("branch", "weak")),
        )
    )


def calculate_isentropic_sweep_route(payload: dict[str, Any]) -> dict[str, Any]:
    """Return sampled isentropic curves for browser plots."""

    points = generate_isentropic_sweep(
        gamma=_number(payload, "gamma", default=1.4),
        minimum_mach=_number(payload, "minimum_mach", default=0.05),
        maximum_mach=_number(payload, "maximum_mach", default=5.0),
        point_count=_integer(payload, "point_count", default=121),
    )
    return {"points": [asdict(point) for point in points]}


def calculate_ideal_gas_route(payload: dict[str, Any]) -> dict[str, Any]:
    """Complete an ideal-gas state from exactly two state properties."""

    state = complete_ideal_gas_state(
        pressure_pa=_optional_number(payload, "pressure_pa"),
        temperature_k=_optional_number(payload, "temperature_k"),
        density_kg_m3=_optional_number(payload, "density_kg_m3"),
        gas_constant_j_kg_k=_number(payload, "gas_constant_j_kg_k", default=287.05),
        gamma=_number(payload, "gamma", default=1.4),
    )
    return asdict(state)


def calculate_stagnation_route(payload: dict[str, Any]) -> dict[str, Any]:
    """Connect static pressure and temperature to total-flow properties."""

    return asdict(
        calculate_stagnation_state(
            mach=_number(payload, "mach"),
            static_pressure_pa=_number(payload, "static_pressure_pa"),
            static_temperature_k=_number(payload, "static_temperature_k"),
            gamma=_number(payload, "gamma", default=1.4),
            gas_constant_j_kg_k=_number(payload, "gas_constant_j_kg_k", default=287.05),
        )
    )


def calculate_process_route(payload: dict[str, Any]) -> dict[str, Any]:
    """Calculate efficiency-corrected compression or expansion."""

    return asdict(
        calculate_isentropic_process(
            inlet_pressure_pa=_number(payload, "inlet_pressure_pa"),
            outlet_pressure_pa=_number(payload, "outlet_pressure_pa"),
            inlet_temperature_k=_number(payload, "inlet_temperature_k"),
            efficiency=_number(payload, "efficiency", default=1.0),
            gamma=_number(payload, "gamma", default=1.4),
            gas_constant_j_kg_k=_number(payload, "gas_constant_j_kg_k", default=287.05),
        )
    )


def calculate_polytropic_route(payload: dict[str, Any]) -> dict[str, Any]:
    """Return polytropic end-state, energy terms, and diagram path."""

    return asdict(
        calculate_polytropic_process(
            inlet_pressure_pa=_number(payload, "inlet_pressure_pa"),
            outlet_pressure_pa=_number(payload, "outlet_pressure_pa"),
            inlet_temperature_k=_number(payload, "inlet_temperature_k"),
            exponent=_number(payload, "exponent"),
            gamma=_number(payload, "gamma", default=1.4),
            gas_constant_j_kg_k=_number(payload, "gas_constant_j_kg_k", default=287.05),
            point_count=_integer(payload, "point_count", default=61),
        )
    )


def calculate_duct_flow_route(payload: dict[str, Any]) -> dict[str, Any]:
    """Return Fanno or Rayleigh line properties for a local Mach number."""

    model = str(payload.get("model", "fanno")).strip().lower()
    mach = _number(payload, "mach")
    gamma = _number(payload, "gamma", default=1.4)
    if model == "fanno":
        return {"model": model, **asdict(calculate_fanno(mach, gamma))}
    if model == "rayleigh":
        return {"model": model, **asdict(calculate_rayleigh(mach, gamma))}
    raise DomainError("Duct-flow model must be fanno or rayleigh.")


def calculate_nozzle_route(payload: dict[str, Any]) -> dict[str, Any]:
    """Calculate ideal choked C-D nozzle performance."""

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
    """Solve a nozzle regime and optionally attach a standard-atmosphere sweep."""

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
    """Generate a conical, preliminary, Rao, or MOC nozzle design."""

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
    """Solve forward or inverse Tsiolkovsky rocket equation inputs."""

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
    """Calculate thrust, equivalent velocity, Isp, and optional impulse."""

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
    if not isinstance(value, dict):
        raise InputError(f"Field {field} must be an object.", field=field)
    return value


def _objects(value: Any, field: str) -> list[dict[str, Any]]:
    if not isinstance(value, list) or not value:
        raise InputError(f"Field {field} must be a non-empty list.", field=field)
    if not all(isinstance(item, dict) for item in value):
        raise InputError(f"Every item in {field} must be an object.", field=field)
    return value


def _enum_value(enum_type: Any, value: Any, field: str) -> Any:
    try:
        return enum_type(str(value))
    except ValueError as error:
        choices = ", ".join(item.value for item in enum_type)
        raise InputError(
            f"Field {field} must be one of: {choices}.", field=field
        ) from error


def simulate_burn_route(payload: dict[str, Any]) -> dict[str, Any]:
    """Evaluate a standalone L0 constant or L1 profiled propulsion burn."""

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
    """Resolve a constant provider shortcut into explicit propulsion streams."""

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
    """Render a verified burn artifact as standalone HTML."""

    result = burn_result_from_dict(_object(payload.get("burn_result"), "burn_result"))
    return {"html": build_burn_html_report(result), "result_hash": result.result_hash}


def export_sidera_burn_route(payload: dict[str, Any]) -> dict[str, Any]:
    """Transform a verified constant burn into a Sidera hand-off artifact."""

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


def sidera_capabilities_route(payload: dict[str, Any]) -> dict[str, Any]:
    """Return the installed Sidera capability record without changing state."""

    del payload
    return asdict(inspect_sidera_capabilities())


def calculate_loss_budget_route(payload: dict[str, Any]) -> dict[str, Any]:
    """Bridge ideal and delivered performance with individually switchable terms."""

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
        )
    )


def _thermal_contour(payload: dict[str, Any]) -> NozzleContour:
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
    """Calculate Bartz wall loading and an optional regenerative bulk balance."""

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
        result["cooling"] = asdict(
            calculate_regenerative_cooling(
                absorbed_heat_w=heat.total_heat_load_w,
                coolant_mass_flow_kg_s=_number(payload, "coolant_mass_flow_kg_s"),
                coolant_specific_heat_j_kg_k=_number(
                    payload, "coolant_specific_heat_j_kg_k"
                ),
                coolant_inlet_temperature_k=_number(
                    payload, "coolant_inlet_temperature_k"
                ),
                coolant_density_kg_m3=_number(payload, "coolant_density_kg_m3"),
                coolant_dynamic_viscosity_pa_s=_number(
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
            )
        )
    return result


def calculate_two_phase_band_route(payload: dict[str, Any]) -> dict[str, Any]:
    """Return an explicit condensed-phase performance uncertainty band."""

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
        )
    )


def calculate_staging_route(payload: dict[str, Any]) -> dict[str, Any]:
    """Optimize multi-stage delta-v allocation for minimum initial mass."""

    raw_stages = payload.get("stages")
    if isinstance(raw_stages, str):
        try:
            raw_stages = json.loads(raw_stages)
        except json.JSONDecodeError as error:
            raise InputError("stages must be valid JSON.", field="stages") from error
    if raw_stages is None and "stage1_specific_impulse_s" in payload:
        raw_stages = [
            {
                "name": str(payload.get("stage1_name", "Booster")),
                "specific_impulse_s": payload.get("stage1_specific_impulse_s"),
                "structural_fraction": payload.get("stage1_structural_fraction"),
            },
            {
                "name": str(payload.get("stage2_name", "Upper")),
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
                name=str(raw_stage.get("name", f"Stage {index + 1}")),
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
    """Close a preliminary feed-system pump/turbine energy balance."""

    return asdict(
        calculate_engine_cycle(
            cycle=str(payload.get("cycle", "pressure-fed")),
            chamber_pressure_pa=_number(payload, "chamber_pressure_pa"),
            chamber_mass_flow_kg_s=_number(payload, "chamber_mass_flow_kg_s"),
            propellant_density_kg_m3=_number(payload, "propellant_density_kg_m3"),
            inlet_pressure_pa=_number(payload, "inlet_pressure_pa"),
            injector_pressure_drop_fraction=_number(
                payload, "injector_pressure_drop_fraction"
            ),
            pump_efficiency=_number(payload, "pump_efficiency"),
            turbine_efficiency=_number(payload, "turbine_efficiency"),
            turbine_inlet_temperature_k=_number(
                payload, "turbine_inlet_temperature_k"
            ),
            turbine_specific_heat_j_kg_k=_number(
                payload, "turbine_specific_heat_j_kg_k"
            ),
            turbine_gamma=_number(payload, "turbine_gamma"),
            turbine_pressure_ratio=_number(payload, "turbine_pressure_ratio"),
        )
    )


def calculate_uncertainty_route(payload: dict[str, Any]) -> dict[str, Any]:
    """Propagate bounded Isp/mass-ratio uncertainty through Tsiolkovsky delta-v."""

    raw_parameters = payload.get("parameters")
    if isinstance(raw_parameters, str):
        try:
            raw_parameters = json.loads(raw_parameters)
        except json.JSONDecodeError as error:
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
            name=str(item.get("name", "")),
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
    result = propagate_uncertainty(
        parameters,
        evaluator=lambda sample: 9.80665
        * sample["specific_impulse_s"]
        * log(sample["mass_ratio"]),
        sample_count=_integer(payload, "sample_count", default=200),
        seed=_integer(payload, "seed", default=42),
        method=str(payload.get("method", "latin-hypercube")),
    )
    return asdict(result)


def validate_workspace_route(payload: dict[str, Any]) -> dict[str, Any]:
    """Migrate and normalize a `.rplab.json` document; optionally build reports."""

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
    """Compare two cases selected from a validated workspace."""

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
    """Return study records for a requested propulsion-system category."""

    category = str(payload.get("category", "all"))
    records = [asdict(record) for record in list_propellants(category)]
    return {"category": category.lower(), "count": len(records), "records": records}


def list_species_route(payload: dict[str, Any]) -> dict[str, Any]:
    """Return the source-labelled variable-property species catalog."""

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
    """Evaluate one species or frozen ideal mixture and return a plot sweep."""

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
    """Solve a provider-neutral TP, HP, or UV equilibrium problem."""

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
    """List chemistry providers and source-labelled propellant presets."""

    del payload
    cea_path = os.environ.get("RPLAB_CEA_EXECUTABLE")
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
        ],
        "propellant_pairs": [asdict(pair) for pair in PROPELLANT_PAIRS.values()],
    }


def calculate_oxidizer_fuel_sweep_route(payload: dict[str, Any]) -> dict[str, Any]:
    """Generate local experimental equilibrium and frozen-expansion O/F curves."""

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
            "NASA SP-8041 Summerfield separation screen",
            "U.S. Standard Atmosphere 1976",
        ),
    ),
    "/api/v1/nozzle-contour": ModelMetadata(
        model="multi-fidelity axisymmetric nozzle contour",
        version="2.0",
        assumptions=(
            "axisymmetric equivalent-area geometry",
            "quasi-one-dimensional isentropic station properties",
            "planar compatibility equations for MOC",
        ),
        units={"*_m": "m", "*_m2": "m2", "mach": "1", "*_ratio": "1"},
        references=(
            "Rao thrust-optimized parabola construction",
            "NASA method-of-characteristics nozzle design literature",
        ),
    ),
    "/api/v1/rocket-equation": ModelMetadata(
        model="Tsiolkovsky ideal rocket equation",
        version="1.0",
        assumptions=("constant effective exhaust velocity", "no gravity or drag losses"),
        units={"*_kg": "kg", "*_m_s": "m/s", "specific_impulse_s": "s", "*_ratio": "1"},
    ),
    "/api/v1/thrust": ModelMetadata(
        model="steady rocket thrust equation",
        version="1.0",
        assumptions=("steady one-dimensional exit plane", "uniform exit properties"),
        units={"*_kg_s": "kg/s", "*_m_s": "m/s", "*_m2": "m2", "*_pa": "Pa", "*_n": "N"},
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
        version="1.0",
        assumptions=("only installed public modules are inspected",),
        units={},
    ),
    "/api/v1/loss-budget": ModelMetadata(
        model="sequential preliminary rocket performance loss budget",
        version="1.0",
        assumptions=(
            "independent lumped efficiencies",
            "sequential attribution without double counting",
            "root-sum-square relative uncertainty",
        ),
        units={"*_n": "N", "*_s": "s", "*_efficiency": "1", "*_residual_*": "varies"},
    ),
    "/api/v1/nozzle-thermal": ModelMetadata(
        model="Bartz gas-side heat transfer and experimental regenerative balance",
        version="1.0-experimental",
        assumptions=(
            "steady turbulent hot-gas boundary layer",
            "user-supplied frozen gas properties",
            "single-phase constant-property coolant when cooling is enabled",
        ),
        units={"*_w_m2": "W/m2", "*_w": "W", "*_k": "K", "*_pa": "Pa", "*_m": "m"},
        references=("Bartz correlation, NASA SP-125", "Darcy-Weisbach pressure loss"),
    ),
    "/api/v1/two-phase-band": ModelMetadata(
        model="condensed-phase particle-coupling performance bracket",
        version="0.1-experimental",
        assumptions=("user-supplied condensed fraction", "bounded particle coupling efficiency"),
        units={"*_n": "N", "*_s": "s", "*_fraction": "1", "*_efficiency*": "1"},
    ),
    "/api/v1/staging": ModelMetadata(
        model="discrete constrained multi-stage mass optimizer",
        version="1.0",
        assumptions=(
            "Tsiolkovsky equation per stage",
            "constant stage Isp and structural fraction",
            "instantaneous staging without gravity or drag loss",
        ),
        units={"*_m_s": "m/s", "*_kg": "kg", "*_s": "s", "*_fraction": "1"},
    ),
    "/api/v1/engine-cycle": ModelMetadata(
        model="preliminary engine-cycle component energy balance",
        version="1.0",
        assumptions=(
            "single equivalent propellant density",
            "constant component efficiencies",
            "no turbomachinery maps or cavitation model",
        ),
        units={"*_pa": "Pa", "*_w": "W", "*_kg_s": "kg/s", "*_j_kg": "J/kg"},
    ),
    "/api/v1/uncertainty": ModelMetadata(
        model="seeded bounded uncertainty propagation",
        version="1.0",
        assumptions=("uniform independent inputs", "Pearson linear sensitivity ranking"),
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
        version="0.1-experimental",
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
        version="0.1-experimental",
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
    """Return JSON-shaped provenance for a registered calculation route."""

    try:
        return asdict(ROUTE_METADATA[path])
    except KeyError as error:
        raise KeyError(path) from error


def dispatch_calculation(path: str, payload: dict[str, Any]) -> dict[str, Any]:
    """Dispatch a versioned API path to its calculation adapter."""

    route = ROUTES.get(path)
    if route is None:
        raise KeyError(path)
    return route(payload)

