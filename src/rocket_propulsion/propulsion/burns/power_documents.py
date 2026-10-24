"""Closed electric/throttle document structures and declared model metadata.

Source identifiers and hashes remain caller-supplied metadata. These structures
neither authenticate a test campaign nor establish physical qualification.
"""

from rocket_propulsion.core.errors import InputError

from .exchange_documents import _array, _number, _object, _string, _strings


def _parameters(values):
    for item in _array(values, "parameters"):
        _object(
            item, {"name", "reference_value", "scale", "minimum", "maximum", "unit"}, "Parameter"
        )
        _string(item["name"], "Parameter name")
        _string(item["unit"], "Parameter unit")
        for key in ("reference_value", "scale", "minimum", "maximum"):
            _number(item[key], key, nullable=key in ("minimum", "maximum"))


def _table_document(payload):
    _object(
        payload,
        {"schema", "source_id", "source_sha256", "reference_ids", "levels", "table_hash"},
        "Electric throttle table",
    )
    for key in ("schema", "source_id", "source_sha256", "table_hash"):
        _string(payload[key], key)
    _strings(payload["reference_ids"], "reference_ids")
    for item in _array(payload["levels"], "levels"):
        _object(item, {"level_id", "input_power_w", "thrust_n", "mass_flow_kg_s"}, "Throttle level")
        _string(item["level_id"], "level_id")
        for key in ("input_power_w", "thrust_n", "mass_flow_kg_s"):
            _number(item[key], key)


_POWER_TEXT = {
    "schema",
    "model_id",
    "model_hash",
    "propellant_state_name",
    "mass_closure_numerics",
    "power_mode_numerics",
    "direction_numerics",
    "force_semantics",
    "mass_semantics",
    "power_semantics",
}
_POWER_FIELDS = _POWER_TEXT | {
    "throttle_law",
    "propellant_reserve_kg",
    "protected_dry_mass_kg",
    "power_margin_w",
    "duty_cycle",
    "parameters",
    "reference_ids",
}
_DISCRETE_FIELDS = {"law", "path_id", "table_hash", "level_ids", "up_switch_hysteresis_w"}
_FIXED_FIELDS = {
    "law",
    "law_id",
    "total_efficiency",
    "specific_impulse_s",
    "minimum_power_w",
    "maximum_power_w",
    "start_hysteresis_w",
}
_POWER_METADATA = (
    "force_semantics",
    "mass_semantics",
    "power_semantics",
    "reference_ids",
    "parameters",
)


def _power_document(payload):
    _object(payload, _POWER_FIELDS, "Power-limited model", ("throttle_table",))
    for key in _POWER_TEXT:
        _string(payload[key], key)
    _strings(payload["reference_ids"], "reference_ids")
    for key in ("propellant_reserve_kg", "protected_dry_mass_kg", "power_margin_w", "duty_cycle"):
        _number(
            payload[key], key, nullable=key in ("propellant_reserve_kg", "protected_dry_mass_kg")
        )
    _parameters(payload["parameters"])
    law = payload["throttle_law"]
    _object(law, {"law"}, "Throttle law", _DISCRETE_FIELDS | _FIXED_FIELDS)
    _string(law["law"], "throttle_law.law")
    if law["law"] == "discrete_throttle_path":
        _object(law, _DISCRETE_FIELDS, "Discrete throttle law")
        for key in ("path_id", "table_hash"):
            _string(law[key], key)
        _strings(law["level_ids"], "level_ids")
        _number(law["up_switch_hysteresis_w"], "up_switch_hysteresis_w")
        if "throttle_table" not in payload:
            raise InputError("A discrete law requires its throttle table.", field="throttle_table")
        _table_document(payload["throttle_table"])
    elif law["law"] == "fixed_efficiency":
        _object(law, _FIXED_FIELDS, "Fixed-efficiency throttle law")
        _string(law["law_id"], "law_id")
        for key in _FIXED_FIELDS - {"law", "law_id"}:
            _number(law[key], key)
        if "throttle_table" in payload:
            raise InputError(
                "A fixed-efficiency law does not consume a throttle table.", field="throttle_table"
            )
    else:
        raise InputError("Unknown throttle law in model document.", field="throttle_law.law")


def _curve_document(value, label):
    if value is None:
        return
    _object(value, {"x", "y"}, label)
    for key in ("x", "y"):
        for number in _array(value[key], label + "." + key):
            _number(number, label + "." + key)


_THROTTLE_TEXT = {
    "schema",
    "model_id",
    "model_hash",
    "direction_numerics",
    "mass_flow_numerics",
    "mass_closure_numerics",
    "curve_numerics",
}
_THROTTLE_FIELDS = _THROTTLE_TEXT | {
    "rated_thrust_n",
    "rated_mass_flow_kg_s",
    "minimum_throttle",
    "maximum_throttle",
    "stream_fractions",
    "thrust_curve",
    "flow_curve",
    "state_constraints",
    "protected_dry_mass_kg",
    "nozzle_exit_area_m2",
    "parameters",
    "reference_ids",
}


def _throttle_document(payload):
    _object(payload, _THROTTLE_FIELDS, "Throttleable engine")
    for key in _THROTTLE_TEXT:
        _string(payload[key], key)
    _strings(payload["reference_ids"], "reference_ids")
    for key in (
        "rated_thrust_n",
        "rated_mass_flow_kg_s",
        "minimum_throttle",
        "maximum_throttle",
        "protected_dry_mass_kg",
        "nozzle_exit_area_m2",
    ):
        _number(payload[key], key, nullable=key in ("protected_dry_mass_kg", "nozzle_exit_area_m2"))
    _parameters(payload["parameters"])
    for pair in _array(payload["stream_fractions"], "stream_fractions"):
        if len(_array(pair, "Stream fraction pair")) != 2:
            raise InputError("Stream fractions require two-entry arrays.", field="stream_fractions")
        _string(pair[0], "Stream fraction name")
        _number(pair[1], "Stream fraction value")
    for key in ("thrust_curve", "flow_curve"):
        _curve_document(payload[key], key)
    for item in _array(payload["state_constraints"], "state_constraints"):
        _object(item, {"state_name", "minimum_value", "label"}, "State constraint")
        _string(item["state_name"], "state_name")
        _string(item["label"], "Constraint label")
        _number(item["minimum_value"], "minimum_value")
