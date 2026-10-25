"""Closed JSON schemas and redundant-metadata consistency for exchange models.

The checks run on raw decoded JSON: objects must have exactly their declared
members with text keys, and each member must have the declared JSON type [1]_
[2]_. Redundant metadata declared in a document must equal the value of the
reconstructed model. Content hashes are reproducible identities, not source
authentication or physical validation. Programmatic constructors keep their
separate iterable and numeric contracts.

References
----------
.. [1] IETF RFC 8259, "The JavaScript Object Notation (JSON) Data Interchange
       Format". https://www.rfc-editor.org/rfc/rfc8259
.. [2] PYTHON-3.12-JSON. Python 3.12, ``json``: encoder/decoder hooks,
       duplicate-member and nonstandard numeric token behavior.
"""

from collections.abc import Callable
from typing import TypeVar, cast

from rocket_propulsion.core.errors import DomainError, InputError
from rocket_propulsion.core.validation import read_json_document, real

from .tabulated_documents import _validate_document as _tabulated_document

_Model = TypeVar("_Model")


def _read_document(
    value: object, label: str, build: Callable[[dict[str, object]], _Model]
) -> _Model:
    """Typed front end of :func:`~rocket_propulsion.core.validation.read_json_document`.

    Parameters
    ----------
    value : object
        Document text, passed unchanged to ``read_json_document``.
    label : str
        Document name used in error messages.
    build : callable
        Reader that turns the decoded root object into a model.

    Returns
    -------
    object
        The model returned by ``build``.

    Raises
    ------
    InputError, DomainError
        Exactly as raised by ``read_json_document``; the root-object check
        below repeats its own and cannot fire on its output.
    """

    def _build(payload: object) -> object:
        if not isinstance(payload, dict):
            raise InputError(f"{label} root must be a JSON object.", field=label)
        return build(payload)

    # read_json_document returns whatever ``_build`` returned, i.e. ``build``'s model.
    return cast(_Model, read_json_document(value, label, _build))


def _object(value, required, label, optional=()):
    """Require a JSON object with text keys and exactly the allowed members.

    Raises
    ------
    InputError
        If ``value`` is not a ``dict`` with ``str`` keys, a required member is
        missing, or an unsupported member is present.
    """
    if type(value) is not dict or any(type(key) is not str for key in value):
        raise InputError(f"{label} must be a JSON object with text keys.", field=label)
    if set(value) - set(required) - set(optional) or set(required) - set(value):
        missing = set(required) - set(value)
        field = next(iter(missing)) if len(missing) == 1 else label
        raise InputError(f"{label} has missing or unsupported fields.", field=field)


def _array(value, label):
    """Return ``value`` if it is a JSON array (a Python ``list``).

    Raises
    ------
    InputError
        If ``value`` is not a ``list``.
    """
    if type(value) is not list:
        raise InputError(f"{label} must be a JSON array.", field=label)
    return value


def _string(value, label, *, nullable=False):
    """Require UTF-8 encodable JSON text, or ``null`` when ``nullable``.

    Raises
    ------
    InputError
        If ``value`` is not a ``str`` or contains unpaired surrogates.
    """
    if value is None and nullable:
        return
    if type(value) is not str:
        raise InputError(f"{label} must be JSON text.", field=label)
    try:
        value.encode("utf-8")
    except UnicodeError as error:
        raise InputError(f"{label} must be UTF-8 encodable.", field=label) from error


def _number(value, label, *, nullable=False):
    """Require a finite JSON number (``int`` or ``float``, not ``bool``).

    ``null`` is accepted when ``nullable``. The value must also pass the
    library's real-number validation.

    Raises
    ------
    InputError
        If ``value`` is not a number, is not finite, or is a nonzero value that
        converts to ``0.0``.
    """
    if value is None and nullable:
        return
    if type(value) not in (int, float):
        raise InputError(f"{label} must be a JSON number.", field=label)
    try:
        real(value, label)
    except DomainError as error:
        raise InputError(error.message, field=label) from error


def _strings(value, label):
    """Require a JSON array whose items are all JSON text."""
    for item in _array(value, label):
        _string(item, label)


def _boolean(value, label):
    """Require a JSON boolean.

    Raises
    ------
    InputError
        If ``value`` is not a ``bool``.
    """
    if type(value) is not bool:
        raise InputError(f"{label} must be a JSON boolean.", field=label)


def _same(left, right):
    """Return JSON equality of two decoded values.

    Integers and floats compare by numeric value (``1 == 1.0``), but booleans
    are never equal to numbers; objects compare by key set and member values,
    arrays element by element.
    """
    if type(left) in (int, float) and type(right) in (int, float):
        return left == right
    if type(left) is not type(right):
        return False
    if type(left) is dict:
        return left.keys() == right.keys() and all(_same(left[k], right[k]) for k in left)
    if type(left) is list:
        return len(left) == len(right) and all(_same(a, b) for a, b in zip(left, right))
    return left == right


def _metadata(payload, manifest, keys):
    """Require each declared metadata key to equal the reconstructed manifest.

    Raises
    ------
    InputError
        If a declared value differs from the manifest value; an absent key
        reads as ``null``.
    """
    for key in keys:
        if not _same(payload.get(key), manifest[key]):
            raise InputError(f"Declared {key} does not match the reconstructed model.", field=key)


def _surface_document(payload, *, streams=False):
    """Check the closed raw schema of a performance-surface document.

    With ``streams`` the points carry per-stream flow scales instead of one
    mass-flow scale.

    Raises
    ------
    InputError
        If the document or a nested axis or point is malformed.
    """
    root = {
        "schema",
        "source_id",
        "source_sha256",
        "axes",
        "points",
        "reference_ids",
        "interpolation_numerics",
        "surface_hash",
    }
    _object(payload, root, "Performance surface", () if streams else ("interpolation",))
    for key in ("schema", "source_id", "source_sha256", "interpolation_numerics", "surface_hash"):
        _string(payload[key], key)
    if "interpolation" in payload:
        _string(payload["interpolation"], "interpolation")
    _strings(payload["reference_ids"], "reference_ids")
    for axis in _array(payload["axes"], "axes"):
        _object(axis, {"name", "unit", "values"}, "Performance axis")
        for key in ("name", "unit"):
            _string(axis[key], key)
        for value in _array(axis["values"], "Axis values"):
            _number(value, "Axis coordinate")
    for point in _array(payload["points"], "points"):
        flow = "stream_flow_scales" if streams else "mass_flow_scale"
        _object(point, {"coordinates", "force_scale", flow}, "Performance point")
        for value in _array(point["coordinates"], "Point coordinates"):
            _number(value, "Point coordinate")
        _number(point["force_scale"], "force_scale")
        if streams:
            for pair in _array(point[flow], flow):
                if len(_array(pair, flow)) != 2:
                    raise InputError("Stream scales require two-entry arrays.", field=flow)
                _string(pair[0], "Stream identifier")
                _number(pair[1], "Stream scale")
        else:
            _number(point[flow], flow)


_BRIDGE_TEXT = {
    "schema",
    "model_id",
    "artifact_hash",
    "force_semantics",
    "mass_semantics",
    "time_semantics",
    "mass_flow_numerics",
    "tabulated_state_numerics",
    "mass_closure_numerics",
    "inventory_clamp_numerics",
    "direction_numerics",
    "time_numerics",
    "model_hash",
}
_BRIDGE_FIELDS = _BRIDGE_TEXT | {
    "performance_surface_hash",
    "operating_condition_axes",
    "stream_bindings",
    "state_constraints",
    "protected_dry_mass_kg",
    "require_complete_stream_binding",
    "parameters",
    "reference_ids",
    "artifact",
    "performance_surface",
}
_BRIDGE_METADATA = (
    "parameters",
    "artifact_hash",
    "performance_surface_hash",
    "stream_performance_surface_hash",
    "operating_condition_axes",
    "force_semantics",
    "mass_semantics",
    "time_semantics",
    "reference_ids",
)


def _bridge_document(payload):
    """Check the closed raw schema of a propagator-bridge document.

    The nested tabulated artifact and performance surfaces are checked with
    their own schemas.

    Raises
    ------
    InputError
        If the document or a nested record is malformed.
    """
    _object(
        payload,
        _BRIDGE_FIELDS,
        "Propagator bridge",
        ("stream_performance_surface", "stream_performance_surface_hash", "nozzle_exit_area_m2"),
    )
    for key in _BRIDGE_TEXT:
        _string(payload[key], key)
    for key in ("performance_surface_hash", "stream_performance_surface_hash"):
        _string(payload.get(key), key, nullable=True)
    _number(payload["protected_dry_mass_kg"], "protected_dry_mass_kg", nullable=True)
    if "nozzle_exit_area_m2" in payload:
        _number(payload["nozzle_exit_area_m2"], "nozzle_exit_area_m2")
    _boolean(payload["require_complete_stream_binding"], "require_complete_stream_binding")
    _strings(payload["reference_ids"], "reference_ids")
    try:
        _tabulated_document(payload["artifact"])
    except InputError:
        raise
    except DomainError as error:
        raise InputError(error.message, field=error.field, details=error.details) from error
    if payload["performance_surface"] is not None:
        _surface_document(payload["performance_surface"])
    if payload.get("stream_performance_surface") is not None:
        _surface_document(payload["stream_performance_surface"], streams=True)
    for item in _array(payload["operating_condition_axes"], "operating_condition_axes"):
        _object(item, {"name", "unit", "minimum", "maximum"}, "Operating condition axis")
        _string(item["name"], "Axis name")
        _string(item["unit"], "Axis unit")
        _number(item["minimum"], "Axis minimum")
        _number(item["maximum"], "Axis maximum")
    for item in _array(payload["stream_bindings"], "stream_bindings"):
        _object(item, {"stream_id", "state_name"}, "Stream binding")
        _string(item["stream_id"], "stream_id")
        _string(item["state_name"], "state_name")
    for item in _array(payload["state_constraints"], "state_constraints"):
        _object(item, {"state_name", "minimum_value", "label"}, "State constraint")
        _string(item["state_name"], "state_name")
        _string(item["label"], "Constraint label")
        _number(item["minimum_value"], "minimum_value")
    for item in _array(payload["parameters"], "parameters"):
        _object(
            item, {"name", "reference_value", "scale", "minimum", "maximum", "unit"}, "Parameter"
        )
        _string(item["name"], "Parameter name")
        _string(item["unit"], "Parameter unit")
        for key in ("reference_value", "scale", "minimum", "maximum"):
            _number(item[key], key, nullable=key in ("minimum", "maximum"))


def _timeline_document(payload):
    """Check the closed raw schema of a propulsion-timeline document.

    The ``bridges`` object must have exactly one member per declared burn.

    Raises
    ------
    InputError
        If the document or a nested record is malformed.
    """
    _object(
        payload,
        {
            "schema",
            "time_numerics",
            "timeline_id",
            "allow_overlapping_burns",
            "burns",
            "jettisons",
            "timeline_hash",
            "bridges",
        },
        "Propulsion timeline",
    )
    for key in ("schema", "time_numerics", "timeline_id", "timeline_hash"):
        _string(payload[key], key)
    _boolean(payload["allow_overlapping_burns"], "allow_overlapping_burns")
    identifiers = []
    for item in _array(payload["burns"], "burns"):
        _object(item, {"burn_id", "ignition_time_s", "bridge_model_hash"}, "Scheduled burn")
        _string(item["burn_id"], "burn_id")
        _string(item["bridge_model_hash"], "bridge_model_hash")
        _number(item["ignition_time_s"], "ignition_time_s")
        identifiers.append(item["burn_id"])
    # Every serialized bridge must be used exactly by a declared burn. Duplicate
    # burn IDs are diagnosed by the timeline constructor, before acceptance.
    _object(payload["bridges"], set(identifiers), "Timeline bridges")
    for bridge in payload["bridges"].values():
        _bridge_document(bridge)
    for item in _array(payload["jettisons"], "jettisons"):
        _object(item, {"event_id", "time_s", "mass_kg", "label"}, "Mass jettison")
        _string(item["event_id"], "event_id")
        _string(item["label"], "Jettison label")
        _number(item["time_s"], "time_s")
        _number(item["mass_kg"], "mass_kg")
