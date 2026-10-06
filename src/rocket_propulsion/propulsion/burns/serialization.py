"""Canonical JSON serialization for burn inputs and results."""

from __future__ import annotations

import hashlib
import json
from dataclasses import fields, is_dataclass
from enum import Enum
from typing import Any

from rocket_propulsion.core.errors import InputError


def _primitive(value: Any) -> Any:
    if is_dataclass(value):
        return {field.name: _primitive(getattr(value, field.name)) for field in fields(value)}
    if isinstance(value, Enum):
        return value.value
    if isinstance(value, tuple):
        return [_primitive(item) for item in value]
    if isinstance(value, list):
        return [_primitive(item) for item in value]
    if isinstance(value, dict):
        return {str(key): _primitive(item) for key, item in value.items()}
    return value


def canonical_json_bytes(value: Any) -> bytes:
    """Encode a public value with deterministic ordering and float spelling."""

    return json.dumps(
        _primitive(value),
        ensure_ascii=False,
        allow_nan=False,
        sort_keys=True,
        separators=(",", ":"),
    ).encode("utf-8")


def _sha256(value: Any) -> str:
    return hashlib.sha256(canonical_json_bytes(value)).hexdigest()


def burn_input_hash(
    definition: Any, operating_point: Any, schedule: Any | None = None
) -> str:
    """Hash normalized physical inputs and model identity."""

    payload = {"definition": definition, "operating_point": operating_point}
    if schedule is not None:
        payload["schedule"] = schedule
    return _sha256(payload)


def burn_result_hash(result: Any) -> str:
    """Hash an authoritative result while excluding its own hash field."""

    payload = _primitive(result)
    payload["result_hash"] = ""
    return _sha256(payload)


def burn_result_to_dict(result: Any) -> dict[str, Any]:
    """Return the JSON-compatible canonical artifact object."""

    payload = _primitive(result)
    if payload.get("schema") not in {
        "rocket_propulsion_burn_v1",
        "rocket_propulsion_burn_v2",
    }:
        raise ValueError("Unsupported burn result schema.")
    return payload


def burn_result_to_json(result: Any, *, indent: int | None = 2) -> str:
    """Serialize a burn result for persistence or transport."""

    return json.dumps(
        burn_result_to_dict(result),
        ensure_ascii=False,
        allow_nan=False,
        sort_keys=True,
        indent=indent,
    )


def _pairs(value: Any, field: str) -> tuple[tuple[str, float], ...]:
    try:
        return tuple((str(item[0]), float(item[1])) for item in value)
    except (TypeError, ValueError, IndexError) as error:
        raise InputError(f"Field {field} must contain name/value pairs.", field=field) from error


def burn_result_from_dict(payload: dict[str, Any]) -> Any:
    """Reconstruct and integrity-check a ``rocket_propulsion_burn_v1`` result."""

    from .models import (
        BurnDefinition,
        BurnProfilePoint,
        BurnResult,
        BurnSummary,
        BurnTarget,
        BurnTargetKind,
        CutoffReason,
        EngineOperatingPoint,
        FlowDestination,
        IntegrationEvidence,
        ProfiledBurnPoint,
        ProfiledBurnResult,
        PropellantFlow,
        PropellantRole,
        PropulsionState,
        TankDefinition,
        ThrottleSchedule,
        ThrottleSegment,
    )

    if not isinstance(payload, dict):
        raise InputError("Burn artifact root must be an object.")
    schema = payload.get("schema")
    if schema not in {"rocket_propulsion_burn_v1", "rocket_propulsion_burn_v2"}:
        raise InputError("Unsupported burn result schema.", field="schema")
    try:
        definition_data = payload["definition"]
        target_data = definition_data["target"]
        definition = BurnDefinition(
            name=str(definition_data["name"]),
            initial_mass_kg=float(definition_data["initial_mass_kg"]),
            protected_dry_mass_kg=float(definition_data["protected_dry_mass_kg"]),
            tanks=tuple(
                TankDefinition(
                    tank_id=str(item["tank_id"]),
                    propellant_key=str(item["propellant_key"]),
                    role=PropellantRole(item["role"]),
                    loaded_mass_kg=float(item["loaded_mass_kg"]),
                    reserve_mass_kg=float(item["reserve_mass_kg"]),
                )
                for item in definition_data["tanks"]
            ),
            target=BurnTarget(BurnTargetKind(target_data["kind"]), float(target_data["value"])),
            maximum_duration_s=float(definition_data["maximum_duration_s"]),
            cant_efficiency=float(definition_data["cant_efficiency"]),
        )
        point_data = payload["operating_point"]
        operating_point = EngineOperatingPoint(
            name=str(point_data["name"]),
            commanded_throttle=float(point_data["commanded_throttle"]),
            realized_throttle=float(point_data["realized_throttle"]),
            active_engine_count=int(point_data["active_engine_count"]),
            ideal_thrust_n=float(point_data["ideal_thrust_n"]),
            delivered_thrust_n=float(point_data["delivered_thrust_n"]),
            chamber_specific_impulse_s=float(point_data["chamber_specific_impulse_s"]),
            system_specific_impulse_s=float(point_data["system_specific_impulse_s"]),
            chamber_mass_flow_kg_s=float(point_data["chamber_mass_flow_kg_s"]),
            total_tank_flow_kg_s=float(point_data["total_tank_flow_kg_s"]),
            oxidizer_mass_flow_kg_s=float(point_data["oxidizer_mass_flow_kg_s"]),
            fuel_mass_flow_kg_s=float(point_data["fuel_mass_flow_kg_s"]),
            mixture_ratio=(
                None if point_data["mixture_ratio"] is None else float(point_data["mixture_ratio"])
            ),
            chamber_pressure_pa=float(point_data["chamber_pressure_pa"]),
            ambient_pressure_pa=float(point_data["ambient_pressure_pa"]),
            streams=tuple(
                PropellantFlow(
                    stream_id=str(item["stream_id"]),
                    tank_id=str(item["tank_id"]),
                    propellant_key=str(item["propellant_key"]),
                    role=PropellantRole(item["role"]),
                    destination=FlowDestination(item["destination"]),
                    mass_flow_kg_s=float(item["mass_flow_kg_s"]),
                    tank_depleting=bool(item["tank_depleting"]),
                    contributes_to_delivered_thrust=bool(
                        item["contributes_to_delivered_thrust"]
                    ),
                )
                for item in point_data["streams"]
            ),
            model_id=str(point_data["model_id"]),
            source=str(point_data["source"]),
            warnings=tuple(str(item) for item in point_data["warnings"]),
        )

        def state(item: dict[str, Any], field: str) -> PropulsionState:
            return PropulsionState(
                time_s=float(item["time_s"]),
                vehicle_mass_kg=float(item["vehicle_mass_kg"]),
                tank_masses_kg=_pairs(item["tank_masses_kg"], field),
            )

        if schema == "rocket_propulsion_burn_v1":
            schedule = None
            profile = tuple(
                BurnProfilePoint(
                    time_s=float(item["time_s"]),
                    delivered_thrust_n=float(item["delivered_thrust_n"]),
                    axial_thrust_n=float(item["axial_thrust_n"]),
                    system_specific_impulse_s=float(item["system_specific_impulse_s"]),
                    total_mass_flow_kg_s=float(item["total_mass_flow_kg_s"]),
                    vehicle_mass_kg=float(item["vehicle_mass_kg"]),
                    tank_masses_kg=_pairs(
                        item["tank_masses_kg"], "profile.tank_masses_kg"
                    ),
                    phase=str(item["phase"]),
                )
                for item in payload["profile"]
            )
        else:
            schedule_data = payload["schedule"]
            schedule = ThrottleSchedule(
                name=str(schedule_data["name"]),
                segments=tuple(
                    ThrottleSegment(
                        duration_s=float(item["duration_s"]),
                        start_throttle=float(item["start_throttle"]),
                        end_throttle=float(item["end_throttle"]),
                        phase=str(item["phase"]),
                    )
                    for item in schedule_data["segments"]
                ),
            )
            profile = tuple(
                ProfiledBurnPoint(
                    time_s=float(item["time_s"]),
                    commanded_throttle=float(item["commanded_throttle"]),
                    realized_throttle=float(item["realized_throttle"]),
                    delivered_thrust_n=float(item["delivered_thrust_n"]),
                    axial_thrust_n=float(item["axial_thrust_n"]),
                    system_specific_impulse_s=float(item["system_specific_impulse_s"]),
                    total_mass_flow_kg_s=float(item["total_mass_flow_kg_s"]),
                    oxidizer_mass_flow_kg_s=float(item["oxidizer_mass_flow_kg_s"]),
                    fuel_mass_flow_kg_s=float(item["fuel_mass_flow_kg_s"]),
                    vehicle_mass_kg=float(item["vehicle_mass_kg"]),
                    tank_masses_kg=_pairs(
                        item["tank_masses_kg"], "profile.tank_masses_kg"
                    ),
                    chamber_pressure_pa=float(item["chamber_pressure_pa"]),
                    phase=str(item["phase"]),
                )
                for item in payload["profile"]
            )
        summary_data = payload["summary"]
        summary = BurnSummary(
            requested_target_kind=BurnTargetKind(summary_data["requested_target_kind"]),
            requested_target_value=float(summary_data["requested_target_value"]),
            target_achieved=bool(summary_data["target_achieved"]),
            cutoff_reason=CutoffReason(summary_data["cutoff_reason"]),
            binding_constraint=(
                None
                if summary_data["binding_constraint"] is None
                else str(summary_data["binding_constraint"])
            ),
            duration_s=float(summary_data["duration_s"]),
            initial_mass_kg=float(summary_data["initial_mass_kg"]),
            final_mass_kg=float(summary_data["final_mass_kg"]),
            consumed_propellant_kg=float(summary_data["consumed_propellant_kg"]),
            consumed_by_tank_kg=_pairs(
                summary_data["consumed_by_tank_kg"], "summary.consumed_by_tank_kg"
            ),
            delivered_total_impulse_n_s=float(summary_data["delivered_total_impulse_n_s"]),
            axial_total_impulse_n_s=float(summary_data["axial_total_impulse_n_s"]),
            average_delivered_thrust_n=float(summary_data["average_delivered_thrust_n"]),
            system_equivalent_specific_impulse_s=float(
                summary_data["system_equivalent_specific_impulse_s"]
            ),
            ideal_delta_v_m_s=float(summary_data["ideal_delta_v_m_s"]),
            rocket_equivalent_specific_impulse_s=float(
                summary_data["rocket_equivalent_specific_impulse_s"]
            ),
            thrust_centroid_time_s=float(summary_data["thrust_centroid_time_s"]),
            mass_closure_error_kg=float(summary_data["mass_closure_error_kg"]),
            impulse_closure_error_n_s=float(summary_data["impulse_closure_error_n_s"]),
            delta_v_closure_error_m_s=float(summary_data["delta_v_closure_error_m_s"]),
            exact_constant_sidera_export_available=bool(
                summary_data["exact_constant_sidera_export_available"]
            ),
        )
        evidence_data = payload["evidence"]
        common = {
            "schema": str(payload["schema"]),
            "definition": definition,
            "operating_point": operating_point,
            "initial_state": state(
                payload["initial_state"], "initial_state.tank_masses_kg"
            ),
            "final_state": state(payload["final_state"], "final_state.tank_masses_kg"),
            "profile": profile,
            "summary": summary,
            "evidence": IntegrationEvidence(
                method=str(evidence_data["method"]),
                absolute_tolerance=float(evidence_data["absolute_tolerance"]),
                relative_tolerance=float(evidence_data["relative_tolerance"]),
                segment_count=int(evidence_data["segment_count"]),
                accepted_steps=int(evidence_data["accepted_steps"]),
                rejected_steps=int(evidence_data["rejected_steps"]),
            ),
            "warnings": tuple(str(item) for item in payload["warnings"]),
            "input_hash": str(payload["input_hash"]),
            "result_hash": str(payload["result_hash"]),
        }
        if schema == "rocket_propulsion_burn_v1":
            result = BurnResult(**common)
        else:
            result = ProfiledBurnResult(schedule=schedule, **common)
    except (KeyError, TypeError, ValueError) as error:
        if isinstance(error, InputError):
            raise
        raise InputError("Malformed burn result artifact.") from error

    declared_schedule = result.schedule if isinstance(result, ProfiledBurnResult) else None
    if result.input_hash != burn_input_hash(
        result.definition, result.operating_point, declared_schedule
    ):
        raise InputError("Burn artifact input hash does not match its inputs.", field="input_hash")
    if result.result_hash != burn_result_hash(result):
        raise InputError("Burn artifact result hash does not match its contents.", field="result_hash")
    # Re-evaluation makes persistence validation stronger than a
    # self-consistent, recomputed hash.
    if isinstance(result, ProfiledBurnResult):
        from .profile_integration import simulate_profiled_burn

        expected = simulate_profiled_burn(
            result.definition, result.operating_point, result.schedule
        )
    else:
        from .integration import simulate_constant_burn

        expected = simulate_constant_burn(result.definition, result.operating_point)
    if canonical_json_bytes(result) != canonical_json_bytes(expected):
        raise InputError("Burn artifact does not reproduce from its declared inputs.")
    return result


def burn_result_from_json(text: str) -> Any:
    """Parse and integrity-check a serialized burn result."""

    try:
        payload = json.loads(text)
    except json.JSONDecodeError as error:
        raise InputError("Burn artifact is not valid JSON.") from error
    return burn_result_from_dict(payload)

