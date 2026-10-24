"""Closed raw JSON records for hierarchical ensembles and their Sidera handoff.

Nested artifacts keep their own reader/identity contract. Source identifiers
and reconstructed content hashes do not authenticate measurements or campaigns.
"""

from rocket_propulsion.core.errors import InputError

from .exchange_documents import _array, _boolean, _number, _object, _string, _strings


def _integer(value, label):
    if type(value) is not int:
        raise InputError(f"{label} must be a JSON integer.", field=label)


def _record(
    value,
    label,
    *,
    texts=(),
    numbers=(),
    integers=(),
    booleans=(),
    arrays=(),
    opaque=(),
    optional=(),
):
    required = set(texts) | set(numbers) | set(integers) | set(booleans) | set(arrays) | set(opaque)
    _object(value, required, label, optional)
    for key in texts:
        _string(value[key], key)
    for key in numbers:
        _number(value[key], key)
    for key in integers:
        _integer(value[key], key)
    for key in booleans:
        _boolean(value[key], key)
    for key in arrays:
        _array(value[key], key)


def _definition(value):
    _record(
        value,
        "Conditional definition",
        texts=("scenario_id", "source", "rationale"),
        integers=("sample_count",),
        arrays=("parameters",),
        optional=("correlation",),
    )
    for parameter in value["parameters"]:
        _record(
            parameter,
            "Uncertain parameter",
            texts=("name", "distribution", "uncertainty_class", "source", "rationale"),
            numbers=("nominal", "lower", "upper"),
            optional=("standard_deviation",),
        )
        if "standard_deviation" in parameter:
            _number(parameter["standard_deviation"], "standard_deviation", nullable=True)
    correlation = value.get("correlation")
    if correlation is not None:
        _record(
            correlation,
            "Conditional correlation",
            texts=("source", "rationale"),
            arrays=("parameter_names", "matrix"),
        )
        _strings(correlation["parameter_names"], "Correlation names")
        for row in correlation["matrix"]:
            for number in _array(row, "Correlation row"):
                _number(number, "Correlation coefficient")


def _ensemble_document(payload):
    _record(
        payload,
        "Hierarchical ensemble",
        texts=("schema", "method", "ensemble_hash", "statistics_numerics", "sampling_numerics"),
        integers=("base_seed", "scenario_count", "realization_count"),
        numbers=(
            "zero_thrust_probability",
            "joint_probability_sum",
            "joint_probability_closure_error",
        ),
        arrays=(
            "definitions",
            "realizations",
            "statistics",
            "variance_decomposition",
            "contributions",
            "sampling_evidence",
            "reference_ids",
            "warnings",
        ),
    )
    _strings(payload["reference_ids"], "reference_ids")
    _strings(payload["warnings"], "warnings")
    for value in payload["definitions"]:
        _definition(value)
    for value in payload["realizations"]:
        _record(
            value,
            "Conditional realization",
            texts=("scenario_id",),
            integers=("conditional_index",),
            numbers=("scenario_probability", "conditional_probability", "joint_probability"),
            arrays=("inputs",),
            opaque=("artifact",),
        )
        for pair in value["inputs"]:
            if len(_array(pair, "Realization input")) != 2:
                raise InputError("Realization inputs require two-entry arrays.", field="inputs")
            _string(pair[0], "Input name")
            _number(pair[1], "Input value")
    for value in payload["statistics"]:
        _record(
            value,
            "Hierarchical statistic",
            texts=("metric",),
            numbers=(
                "expected_value",
                "standard_deviation",
                "minimum",
                "p05",
                "p50",
                "p95",
                "maximum",
            ),
        )
    for value in payload["variance_decomposition"]:
        _record(
            value,
            "Variance decomposition",
            texts=("metric",),
            numbers=(
                "total_variance",
                "within_scenario_variance",
                "between_scenario_variance",
                "closure_error",
                "within_fraction",
                "between_fraction",
            ),
        )
    for value in payload["contributions"]:
        _record(
            value,
            "Scenario contribution",
            texts=("scenario_id", "metric"),
            numbers=(
                "scenario_probability",
                "conditional_expected_value",
                "joint_expected_value_contribution",
            ),
        )
    for value in payload["sampling_evidence"]:
        _record(
            value,
            "Sampling evidence",
            texts=("scenario_id",),
            integers=("sample_count", "primary_seed", "audit_seed"),
            booleans=("converged",),
            numbers=(
                "impulse_mean_relative_difference",
                "impulse_standard_deviation_relative_difference",
                "requested_tolerance",
            ),
        )


def _handoff_document(payload):
    _record(
        payload,
        "Sidera ensemble handoff",
        texts=(
            "schema",
            "source_ensemble_schema",
            "source_ensemble_hash",
            "force_channel",
            "mass_flow_channel",
            "direction_semantics",
            "direction_numerics",
            "handoff_hash",
        ),
        numbers=("joint_probability_sum",),
        arrays=("required_capabilities", "jobs", "warnings"),
    )
    for key in ("required_capabilities", "warnings"):
        _strings(payload[key], key)
    for value in payload["jobs"]:
        _record(
            value,
            "Sidera ensemble job",
            texts=("scenario_id", "frame", "artifact_hash"),
            integers=("conditional_index",),
            arrays=("direction",),
            opaque=("artifact",),
            numbers=(
                "scenario_probability",
                "conditional_probability",
                "joint_probability",
                "t_start_s",
                "t_end_s",
            ),
        )
        if len(value["direction"]) != 3:
            raise InputError("Sidera direction requires three coordinates.", field="direction")
        for number in value["direction"]:
            _number(number, "Direction coordinate")
