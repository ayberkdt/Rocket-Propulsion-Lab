"""Canonical JSON serialization and integrity checks for burn inputs and results.

Canonical bytes are UTF-8 JSON with sorted keys, compact separators and no
``NaN`` or infinity tokens [1]_ [2]_; SHA-256 digests of these bytes are
reproducible content identities, not authentication. Reading a burn result
checks the closed record structure, both hashes, and that the result reproduces
from its declared inputs.

References
----------
.. [1] IETF RFC 8259, "The JavaScript Object Notation (JSON) Data Interchange
       Format". https://www.rfc-editor.org/rfc/rfc8259
.. [2] PYTHON-3.12-JSON. Python 3.12, ``json``: encoder/decoder hooks,
       duplicate-member and nonstandard numeric token behavior.
"""

from __future__ import annotations

import hashlib
import json
from dataclasses import fields, is_dataclass
from enum import Enum
from typing import Any

from rocket_propulsion.core.errors import DomainError, InputError
from rocket_propulsion.core.validation import DOCUMENT_ERRORS, json_text


def _primitive(value: Any) -> Any:
    """Convert dataclasses, enums, tuples, lists and dicts to JSON-compatible values.

    Dataclasses become objects of their fields, enums their values, tuples
    arrays and dictionary keys strings; other values are returned unchanged.
    """
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
    """Encode a value as canonical UTF-8 JSON bytes.

    Keys are sorted, separators are compact, non-ASCII text is kept as Unicode
    and floats use Python's shortest round-trip spelling [1]_ [2]_.

    Parameters
    ----------
    value : object
        Dataclass record or JSON-compatible value.

    Returns
    -------
    bytes
        Canonical encoding.

    Raises
    ------
    ValueError
        If the value contains ``NaN`` or an infinity.
    TypeError
        If the value contains an object that is not JSON serializable.

    References
    ----------
    .. [1] IETF RFC 8259, "The JavaScript Object Notation (JSON) Data Interchange
           Format". https://www.rfc-editor.org/rfc/rfc8259
    .. [2] PYTHON-3.12-JSON. Python 3.12, ``json``: encoder/decoder hooks,
           duplicate-member and nonstandard numeric token behavior.
    """

    return json.dumps(
        _primitive(value),
        ensure_ascii=False,
        allow_nan=False,
        sort_keys=True,
        separators=(",", ":"),
    ).encode("utf-8")


def canonical_sha256(value: Any) -> str:
    """Return the SHA-256 hex digest of :func:`canonical_json_bytes`.

    Parameters
    ----------
    value : object
        Value to identify.

    Returns
    -------
    str
        Lower-case hexadecimal digest; a reproducible content identity, not an
        authentication.
    """

    return hashlib.sha256(canonical_json_bytes(value)).hexdigest()


def _require_numerical_policies(payload, manifest) -> None:
    """Require every ``*_numerics`` policy of ``manifest`` to match ``payload``.

    Raises
    ------
    InputError
        If a policy identifier is missing from the payload or differs.
    """
    for key, value in manifest.items():
        if key.endswith("_numerics") and payload.get(key) != value:
            raise InputError(f"Model numerical policy {key} is missing or unsupported.", field=key)


def _sha256(value: Any) -> str:
    """Return the SHA-256 hex digest of the canonical JSON of ``value``."""
    return hashlib.sha256(canonical_json_bytes(value)).hexdigest()


def burn_input_hash(
    definition: Any, operating_point: Any, schedule: Any | None = None
) -> str:
    """Return the SHA-256 identity of a burn's normalized physical inputs.

    Parameters
    ----------
    definition : BurnDefinition
        Burn definition.
    operating_point : EngineOperatingPoint
        Engine operating point.
    schedule : ThrottleSchedule or None, optional
        Throttle schedule of a profiled burn; omitted from the hash when
        ``None``.

    Returns
    -------
    str
        Lower-case hexadecimal digest.
    """

    payload = {"definition": definition, "operating_point": operating_point}
    if schedule is not None:
        payload["schedule"] = schedule
    return _sha256(payload)


def burn_result_hash(result: Any) -> str:
    """Return the SHA-256 identity of a burn result.

    The ``result_hash`` field is blanked before hashing so the digest does not
    depend on itself.

    Parameters
    ----------
    result : BurnResult or ProfiledBurnResult
        Result record.

    Returns
    -------
    str
        Lower-case hexadecimal digest.
    """

    payload = _primitive(result)
    payload["result_hash"] = ""
    return _sha256(payload)


def burn_result_to_dict(result: Any) -> dict[str, Any]:
    """Export a burn result as a JSON-compatible document after verifying it.

    The exported document is read back with :func:`burn_result_from_dict`, so
    only a typed, hash-consistent result that reproduces from its inputs is
    exported.

    Parameters
    ----------
    result : BurnResult or ProfiledBurnResult
        Result record.

    Returns
    -------
    dict
        JSON-compatible document.

    Raises
    ------
    InputError
        If ``result`` is not a supported record or does not verify.
    """

    from .models import BurnResult, ProfiledBurnResult

    if not isinstance(result, (BurnResult, ProfiledBurnResult)):
        raise InputError('Burn result must be a supported result record.', field='result')
    payload = _primitive(result)
    burn_result_from_dict(payload)
    return payload


def burn_result_to_json(result: Any, *, indent: int | None = 2) -> str:
    """Serialize a verified burn result as JSON text with sorted keys.

    Non-finite numbers are rejected, so the output is standard JSON [1]_ [2]_.

    Parameters
    ----------
    result : BurnResult or ProfiledBurnResult
        Result record.
    indent : int or None, optional
        Indentation passed to ``json.dumps``.

    Returns
    -------
    str
        JSON text of :func:`burn_result_to_dict`.

    Raises
    ------
    InputError
        If the result does not verify.

    References
    ----------
    .. [1] IETF RFC 8259, "The JavaScript Object Notation (JSON) Data Interchange
           Format". https://www.rfc-editor.org/rfc/rfc8259
    .. [2] PYTHON-3.12-JSON. Python 3.12, ``json``: encoder/decoder hooks,
           duplicate-member and nonstandard numeric token behavior.
    """

    return json.dumps(
        burn_result_to_dict(result),
        ensure_ascii=False,
        allow_nan=False,
        sort_keys=True,
        indent=indent,
    )


def _record_data(value, kind, label):
    """Return a copy of ``value`` if it has exactly the fields of ``kind``.

    Extra members are rejected rather than discarded, so no data outside the
    hashed record can be silently dropped.

    Raises
    ------
    InputError
        If ``value`` is not an object or its member set differs.
    """
    if not isinstance(value, dict):
        raise InputError(f"{label} must be an object.", field=label)
    expected = {item.name for item in fields(kind)}
    if set(value) != expected:
        raise InputError(f"{label} has missing or unsupported fields.", field=label)
    return dict(value)


def _array(value, label):
    """Return a list or tuple as a tuple; strings and mappings are rejected.

    Raises
    ------
    InputError
        If ``value`` is not a list or tuple.
    """
    if not isinstance(value, (list, tuple)):
        raise InputError(f"{label} must be an array.", field=label)
    return tuple(value)


def burn_result_from_dict(payload: dict[str, Any]) -> Any:
    """Reopen a burn-result document and verify that it reproduces.

    Parameters
    ----------
    payload : dict
        Document with schema ``rocket_propulsion_burn_v1`` (constant burn) or
        ``rocket_propulsion_burn_v2`` (profiled burn).

    Returns
    -------
    BurnResult or ProfiledBurnResult
        The validated result.

    Raises
    ------
    InputError
        If the schema is unsupported, the record structure is not exactly the
        declared one, the input or result hash does not match, or rerunning the
        solver on the declared inputs does not reproduce the stored result.

    Notes
    -----
    A matching hash only shows that the content is unaltered. Re-simulation
    is required to show that the outputs follow from the declared inputs.
    """
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
        expected: BurnResult | ProfiledBurnResult
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
    """Build an object from decoded members, rejecting duplicate names.

    RFC 8259 leaves the meaning of duplicate member names undefined [1]_ and
    the Python decoder keeps the last one by default [2]_; duplicates are
    therefore rejected at every depth before model coercion.

    Raises
    ------
    InputError
        If a member name repeats.

    References
    ----------
    .. [1] IETF RFC 8259, "The JavaScript Object Notation (JSON) Data Interchange
           Format". https://www.rfc-editor.org/rfc/rfc8259
    .. [2] PYTHON-3.12-JSON. Python 3.12, ``json``: encoder/decoder hooks,
           duplicate-member and nonstandard numeric token behavior.
    """
    result = {}
    for name, value in pairs:
        if name in result:
            raise InputError(f"Duplicate burn artifact field {name!r}.", field=name)
        result[name] = value
    return result


def _reject_constant(value):
    """Reject the nonstandard ``NaN``, ``Infinity`` and ``-Infinity`` tokens.

    These tokens are not JSON [1]_ but are accepted by the Python decoder by
    default [2]_.

    Raises
    ------
    InputError
        Always.

    References
    ----------
    .. [1] IETF RFC 8259, "The JavaScript Object Notation (JSON) Data Interchange
           Format". https://www.rfc-editor.org/rfc/rfc8259
    .. [2] PYTHON-3.12-JSON. Python 3.12, ``json``: encoder/decoder hooks,
           duplicate-member and nonstandard numeric token behavior.
    """
    raise InputError(f"Nonstandard JSON constant {value!r}.", field="JSON constant")


def burn_result_from_json(text: str) -> Any:
    """Parse strict JSON text and verify the complete burn result.

    Duplicate member names and nonstandard numeric tokens are rejected [1]_
    [2]_ before :func:`burn_result_from_dict` validates the document.

    Parameters
    ----------
    text : str
        JSON document text.

    Returns
    -------
    BurnResult or ProfiledBurnResult
        The validated result.

    Raises
    ------
    InputError
        If the text is not strict JSON or the document does not verify.

    References
    ----------
    .. [1] IETF RFC 8259, "The JavaScript Object Notation (JSON) Data Interchange
           Format". https://www.rfc-editor.org/rfc/rfc8259
    .. [2] PYTHON-3.12-JSON. Python 3.12, ``json``: encoder/decoder hooks,
           duplicate-member and nonstandard numeric token behavior.
    """
    source = json_text(text, "Burn artifact")
    try:
        payload = json.loads(source, object_pairs_hook=_unique_object, parse_constant=_reject_constant)
    except InputError:
        raise
    except DOCUMENT_ERRORS as error:
        raise InputError("Burn artifact is not valid JSON.", field="artifact") from error
    return burn_result_from_dict(payload)
