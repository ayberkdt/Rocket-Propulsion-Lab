"""Canonical JSON serialization for burn inputs and results."""

from __future__ import annotations

import hashlib
import json
from dataclasses import fields, is_dataclass
from enum import Enum
from typing import Any

from rocket_propulsion.core.errors import DomainError, InputError
from rocket_propulsion.core.validation import DOCUMENT_ERRORS, json_text


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


def canonical_sha256(value: Any) -> str:
    """Return the SHA-256 hex digest of :func:`canonical_json_bytes`."""

    return hashlib.sha256(canonical_json_bytes(value)).hexdigest()


def _require_numerical_policies(payload, manifest) -> None:
    """Reject missing or changed numerical policies recognized by this model."""
    for key, value in manifest.items():
        if key.endswith("_numerics") and payload.get(key) != value:
            raise InputError(f"Model numerical policy {key} is missing or unsupported.", field=key)


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
    """Export only a typed, hash-consistent result reproducible from its inputs."""

    from .models import BurnResult, ProfiledBurnResult

    if not isinstance(result, (BurnResult, ProfiledBurnResult)):
        raise InputError('Burn result must be a supported result record.', field='result')
    payload = _primitive(result)
    burn_result_from_dict(payload)
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


def _record_data(value, kind, label):
    """Require precisely the declared record fields; never discard unhashed data."""
    if not isinstance(value, dict):
        raise InputError(f"{label} must be an object.", field=label)
    expected = {item.name for item in fields(kind)}
    if set(value) != expected:
        raise InputError(f"{label} has missing or unsupported fields.", field=label)
    return dict(value)


def _array(value, label):
    """Accept JSON-compatible arrays without iterating strings or mappings."""
    if not isinstance(value, (list, tuple)):
        raise InputError(f"{label} must be an array.", field=label)
    return tuple(value)


def burn_result_from_dict(payload: dict[str, Any]) -> Any:
    """Validate exact v1/v2 structure, raw types, hashes and solver reproduction."""
    from .models import (
        BurnDefinition,
        BurnProfilePoint,
        BurnResult,
        BurnSummary,
        BurnTarget,
        EngineOperatingPoint,
        IntegrationEvidence,
        ProfiledBurnPoint,
        ProfiledBurnResult,
        PropellantFlow,
        PropulsionState,
        TankDefinition,
        ThrottleSchedule,
        ThrottleSegment,
    )

    if not isinstance(payload, dict):
        raise InputError("Burn artifact root must be an object.", field="artifact")
    schema = payload.get("schema")
    if not isinstance(schema, str) or schema not in {"rocket_propulsion_burn_v1", "rocket_propulsion_burn_v2"}:
        raise InputError("Unsupported burn result schema.", field="schema")
    result_kind = BurnResult if schema == "rocket_propulsion_burn_v1" else ProfiledBurnResult
    point_kind = BurnProfilePoint if result_kind is BurnResult else ProfiledBurnPoint

    def record(kind, value, label):
        return kind(**_record_data(value, kind, label))

    try:
        data = _record_data(payload, result_kind, "artifact")
        definition_data = _record_data(data["definition"], BurnDefinition, "definition")
        definition = BurnDefinition(**{
            **definition_data,
            "target": record(BurnTarget, definition_data["target"], "definition.target"),
            "tanks": tuple(record(TankDefinition, item, "definition.tanks")
                           for item in _array(definition_data["tanks"], "definition.tanks")),
        })
        point_data = _record_data(data["operating_point"], EngineOperatingPoint, "operating_point")
        operating_point = EngineOperatingPoint(**{
            **point_data,
            "streams": tuple(record(PropellantFlow, item, "operating_point.streams")
                             for item in _array(point_data["streams"], "operating_point.streams")),
        })
        converted = {
            **data, "definition": definition, "operating_point": operating_point,
            "initial_state": record(PropulsionState, data["initial_state"], "initial_state"),
            "final_state": record(PropulsionState, data["final_state"], "final_state"),
            "summary": record(BurnSummary, data["summary"], "summary"),
            "evidence": record(IntegrationEvidence, data["evidence"], "evidence"),
            "profile": tuple(record(point_kind, item, "profile")
                             for item in _array(data["profile"], "profile")),
        }
        if result_kind is ProfiledBurnResult:
            schedule_data = _record_data(data["schedule"], ThrottleSchedule, "schedule")
            converted["schedule"] = ThrottleSchedule(**{
                **schedule_data,
                "segments": tuple(record(ThrottleSegment, item, "schedule.segments")
                                  for item in _array(schedule_data["segments"], "schedule.segments")),
            })
        result = result_kind(**converted)
        declared_schedule = result.schedule if isinstance(result, ProfiledBurnResult) else None
        if result.input_hash != burn_input_hash(result.definition, result.operating_point, declared_schedule):
            raise InputError("Burn artifact input hash does not match its inputs.", field="input_hash")
        if result.result_hash != burn_result_hash(result):
            raise InputError("Burn artifact result hash does not match its contents.", field="result_hash")
        # A recomputed hash is not evidence that the output follows its inputs.
        if isinstance(result, ProfiledBurnResult):
            from .profile_integration import simulate_profiled_burn
            expected = simulate_profiled_burn(result.definition, result.operating_point, result.schedule)
        else:
            from .integration import simulate_constant_burn
            expected = simulate_constant_burn(result.definition, result.operating_point)
        if canonical_json_bytes(result) != canonical_json_bytes(expected):
            raise InputError("Burn artifact does not reproduce from its declared inputs.", field="result")
        return result
    except InputError:
        raise
    except DomainError as error:
        raise InputError(f"Malformed burn result artifact: {error}", field=error.field) from error
    except DOCUMENT_ERRORS + (UnicodeError,) as error:
        raise InputError("Malformed burn result artifact.", field="artifact") from error


def _unique_object(pairs):
    """Reject ambiguous object members at every JSON depth before model coercion."""
    result = {}
    for name, value in pairs:
        if name in result:
            raise InputError(f"Duplicate burn artifact field {name!r}.", field=name)
        result[name] = value
    return result


def _reject_constant(value):
    """Keep nonstandard NaN and infinity tokens out of the whole artifact."""
    raise InputError(f"Nonstandard JSON constant {value!r}.", field="JSON constant")


def burn_result_from_json(text: str) -> Any:
    """Parse strict JSON and validate the complete burn artifact."""
    source = json_text(text, "Burn artifact")
    try:
        payload = json.loads(source, object_pairs_hook=_unique_object, parse_constant=_reject_constant)
    except InputError:
        raise
    except DOCUMENT_ERRORS as error:
        raise InputError("Burn artifact is not valid JSON.", field="artifact") from error
    return burn_result_from_dict(payload)
