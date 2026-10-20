"""Closed raw tabulated schema; normalized model hashes are not authentication.

JSON arrays are lists at this boundary. Constructors accept their documented
iterables separately. Empty references mean absence; JSON integers normalize
through the model. No unknown extension field can disappear during a read.
"""

from rocket_propulsion.core.errors import DomainError, InputError
from rocket_propulsion.core.validation import real, text

_NUMBERS = frozenset(
    (
        "start_time_s",
        "end_time_s",
        "start_delivered_thrust_n",
        "end_delivered_thrust_n",
        "start_axial_thrust_n",
        "end_axial_thrust_n",
        "start_total_mass_flow_kg_s",
        "end_total_mass_flow_kg_s",
    )
)
_FLOWS = ("start_stream_mass_flows_kg_s", "end_stream_mass_flows_kg_s")
_SUMMARY = frozenset(
    (
        "duration_s",
        "delivered_total_impulse_n_s",
        "axial_total_impulse_n_s",
        "consumed_propellant_kg",
        "thrust_centroid_time_s",
        "system_equivalent_specific_impulse_s",
    )
)
_ROOT = frozenset(
    (
        "schema",
        "integral_numerics",
        "time_numerics",
        "source_id",
        "source_sha256",
        "interpolation",
        "segments",
        "warnings",
        "summary",
        "artifact_hash",
    )
)


def _utf8(value, label, *, allow_blank=False):
    if allow_blank:
        if not isinstance(value, str):
            raise DomainError(f"{label} must be text.", field=label)
        value = str.__str__(value)
    else:
        value = text(value, label)
    try:
        value.encode("utf-8")
    except UnicodeError as error:
        raise DomainError(f"{label} must be UTF-8 encodable.", field=label) from error
    return value


def _object(value, required, optional, label):
    if (
        not isinstance(value, dict)
        or not required <= value.keys()
        or value.keys() - required - optional
    ):
        raise InputError(f"{label} has missing, unsupported or invalid fields.", field=label)


def _array(value, label):
    if not isinstance(value, list):
        raise InputError(f"{label} must be a JSON array.", field=label)
    return value


def _number(value, label, *, nullable=False):
    if value is None and nullable:
        return
    if type(value) not in (int, float):
        raise InputError(f"{label} must be a JSON number.", field=label)
    real(value, label)


def _validate_document(payload):
    """Check all declared structure before reconstruction or hash comparison."""
    _object(payload, _ROOT, {"reference_ids"}, "Tabulated artifact")
    for key in (
        "schema",
        "integral_numerics",
        "time_numerics",
        "source_id",
        "source_sha256",
        "interpolation",
        "artifact_hash",
    ):
        _utf8(payload[key], key)
    for key in ("warnings", "reference_ids"):
        for item in _array(payload.get(key, []), key):
            _utf8(item, key, allow_blank=key == "warnings")
    for segment in _array(payload["segments"], "segments"):
        _object(segment, _NUMBERS | set(_FLOWS) | {"phase"}, set(), "Tabulated segment")
        for key in _NUMBERS:
            _number(segment[key], key)
        _utf8(segment["phase"], "phase")
        for key in _FLOWS:
            for pair in _array(segment[key], key):
                if len(_array(pair, key)) != 2:
                    raise InputError("Stream flow pairs require two entries.", field=key)
                _utf8(pair[0], "Stream identifier")
                _number(pair[1], "Stream flow")
    _object(payload["summary"], _SUMMARY, set(), "Tabulated summary")
    for key in _SUMMARY:
        _number(
            payload["summary"][key], key, nullable=key == "system_equivalent_specific_impulse_s"
        )
