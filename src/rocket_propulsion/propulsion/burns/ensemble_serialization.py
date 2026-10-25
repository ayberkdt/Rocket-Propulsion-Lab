"""Persistence and integrity checks for hierarchical propulsion ensembles.

The exchange document is target-neutral: it carries propulsion histories, joint
probabilities, evidence and provenance, but no epoch, frame, attitude or orbital
state, so a trajectory tool such as Sidera can consume it without becoming a
dependency of the propulsion core. Documents are standard JSON [1]_ [2]_.

Reading a document recomputes every derived quantity (probability closure,
weighted statistics, variance decomposition, scenario contributions and seed
links) and rejects any mismatch [3]_ [4]_. The content hash is a reproducible
identity, not an authentication of the data.

References
----------
.. [1] IETF RFC 8259, "The JavaScript Object Notation (JSON) Data Interchange
       Format". https://www.rfc-editor.org/rfc/rfc8259
.. [2] PYTHON-3.12-JSON. Python 3.12, ``json``: encoder/decoder hooks,
       duplicate-member and nonstandard numeric token behavior.
.. [3] NASA-STD-7009B, "Standard for Models and Simulations".
       https://standards.nasa.gov/standard/NASA/NASA-STD-7009
.. [4] NASA/SP-2011-3421, "Probabilistic Risk Assessment Procedures Guide for
       NASA Managers and Practitioners". https://ntrs.nasa.gov/citations/20120001369
"""

from __future__ import annotations

import hashlib
import json
from dataclasses import asdict
from fractions import Fraction
from math import isclose, isfinite
from typing import Any

from rocket_propulsion.core.errors import DomainError, InputError
from rocket_propulsion.core.validation import read_json_document

from .ensemble import (
    ConditionalPropulsionRealization,
    ConditionalSamplingEvidence,
    HierarchicalMetricStatistics,
    HierarchicalPropulsionEnsemble,
    HierarchicalVarianceDecomposition,
    ScenarioConditionalUncertainty,
    ScenarioMetricContribution,
    _derived_seed,
    _variance_decomposition,
    _weighted_statistics,
)
from .ensemble_documents import _ensemble_document
from .marginal_numerics import _SAMPLING_POLICY
from .numeric_outputs import _finite_rational, _fraction
from .probability_numerics import _PROBABILITY_POLICY, _weight_ledger
from .serialization import _require_numerical_policies, canonical_json_bytes
from .tabulated import TabulatedBurnArtifact
from .uncertainty import (
    CorrelationModel,
    DistributionKind,
    UncertainParameter,
    UncertaintyClass,
)

HIERARCHICAL_ENSEMBLE_SCHEMA = "rocket_propulsion_hierarchical_ensemble_v1"

_METRIC_GETTERS = {
    "duration_s": lambda artifact: artifact.duration_s,
    "delivered_impulse_n_s": lambda artifact: artifact.delivered_total_impulse_n_s,
    "axial_impulse_n_s": lambda artifact: artifact.axial_total_impulse_n_s,
    "consumed_propellant_kg": lambda artifact: artifact.consumed_propellant_kg,
    "thrust_centroid_time_s": lambda artifact: artifact.thrust_centroid_time_s,
}


def _parameter_to_dict(parameter: UncertainParameter) -> dict[str, Any]:
    """Return the JSON-compatible record of an uncertain parameter."""
    return {
        "name": parameter.name,
        "nominal": parameter.nominal,
        "lower": parameter.lower,
        "upper": parameter.upper,
        "distribution": parameter.distribution.value,
        "uncertainty_class": parameter.uncertainty_class.value,
        "source": parameter.source,
        "rationale": parameter.rationale,
        "standard_deviation": parameter.standard_deviation,
    }


def _parameter_from_dict(payload: dict[str, Any]) -> UncertainParameter:
    """Rebuild and validate an uncertain parameter from its record."""
    return UncertainParameter(
        name=str(payload["name"]),
        nominal=float(payload["nominal"]),
        lower=float(payload["lower"]),
        upper=float(payload["upper"]),
        distribution=DistributionKind(payload["distribution"]),
        uncertainty_class=UncertaintyClass(payload["uncertainty_class"]),
        source=str(payload["source"]),
        rationale=str(payload["rationale"]),
        standard_deviation=(
            None
            if payload.get("standard_deviation") is None
            else float(payload["standard_deviation"])
        ),
    )


def _correlation_to_dict(correlation: CorrelationModel | None) -> dict[str, Any] | None:
    """Return the JSON-compatible record of a correlation model, or ``None``."""
    if correlation is None:
        return None
    return {
        "parameter_names": list(correlation.parameter_names),
        "matrix": [list(row) for row in correlation.matrix],
        "source": correlation.source,
        "rationale": correlation.rationale,
    }


def _correlation_from_dict(payload: Any) -> CorrelationModel | None:
    """Rebuild and validate a correlation model from its record or ``None``.

    Raises
    ------
    InputError
        If the record is neither an object nor ``None``.
    """
    if payload is None:
        return None
    if not isinstance(payload, dict):
        raise InputError("Conditional correlation must be an object or null.")
    return CorrelationModel(
        parameter_names=tuple(str(item) for item in payload["parameter_names"]),
        matrix=tuple(
            tuple(float(value) for value in row) for row in payload["matrix"]
        ),
        source=str(payload["source"]),
        rationale=str(payload["rationale"]),
    )


def _definition_to_dict(
    definition: ScenarioConditionalUncertainty,
) -> dict[str, Any]:
    """Return the JSON-compatible record of a conditional definition."""
    return {
        "scenario_id": definition.scenario_id,
        "parameters": [
            _parameter_to_dict(parameter) for parameter in definition.parameters
        ],
        "sample_count": definition.sample_count,
        "source": definition.source,
        "rationale": definition.rationale,
        "correlation": _correlation_to_dict(definition.correlation),
    }


def _definition_from_dict(payload: dict[str, Any]) -> ScenarioConditionalUncertainty:
    """Rebuild and validate a conditional definition from its record."""
    return ScenarioConditionalUncertainty(
        scenario_id=str(payload["scenario_id"]),
        parameters=tuple(
            _parameter_from_dict(parameter) for parameter in payload["parameters"]
        ),
        sample_count=int(payload["sample_count"]),
        source=str(payload["source"]),
        rationale=str(payload["rationale"]),
        correlation=_correlation_from_dict(payload.get("correlation")),
    )


def hierarchical_ensemble_to_dict(
    ensemble: HierarchicalPropulsionEnsemble,
    *,
    include_hash: bool = True,
) -> dict[str, Any]:
    """Return a deterministic JSON-compatible ensemble document.

    Parameters
    ----------
    ensemble : HierarchicalPropulsionEnsemble
        Ensemble to export.
    include_hash : bool, optional
        Add ``ensemble_hash``, the SHA-256 of the canonical JSON of the
        document without that field.

    Returns
    -------
    dict
        Document with schema, numerical-policy identifiers, definitions,
        realizations (each with its nested artifact document), statistics,
        evidence, references and warnings [1]_ [2]_.

    References
    ----------
    .. [1] NASA-STD-7009B, "Standard for Models and Simulations".
           https://standards.nasa.gov/standard/NASA/NASA-STD-7009
    .. [2] NASA/SP-2011-3421, "Probabilistic Risk Assessment Procedures Guide for
           NASA Managers and Practitioners". https://ntrs.nasa.gov/citations/20120001369
    """

    payload: dict[str, Any] = {
        "schema": HIERARCHICAL_ENSEMBLE_SCHEMA,
        "statistics_numerics": _PROBABILITY_POLICY,
        "sampling_numerics": _SAMPLING_POLICY,
        "method": ensemble.method,
        "base_seed": ensemble.base_seed,
        "scenario_count": ensemble.scenario_count,
        "realization_count": ensemble.realization_count,
        "definitions": [
            _definition_to_dict(definition) for definition in ensemble.definitions
        ],
        "realizations": [
            {
                "scenario_id": realization.scenario_id,
                "scenario_probability": realization.scenario_probability,
                "conditional_index": realization.conditional_index,
                "conditional_probability": realization.conditional_probability,
                "joint_probability": realization.joint_probability,
                "inputs": [list(item) for item in realization.inputs],
                "artifact": realization.artifact.to_dict(),
            }
            for realization in ensemble.realizations
        ],
        "statistics": [asdict(item) for item in ensemble.statistics],
        "variance_decomposition": [
            asdict(item) for item in ensemble.variance_decomposition
        ],
        "contributions": [asdict(item) for item in ensemble.contributions],
        "sampling_evidence": [asdict(item) for item in ensemble.sampling_evidence],
        "zero_thrust_probability": ensemble.zero_thrust_probability,
        "joint_probability_sum": ensemble.joint_probability_sum,
        "joint_probability_closure_error": ensemble.joint_probability_closure_error,
        "reference_ids": list(ensemble.reference_ids),
        "warnings": list(ensemble.warnings),
    }
    if include_hash:
        payload["ensemble_hash"] = hashlib.sha256(
            canonical_json_bytes(payload)
        ).hexdigest()
    return payload


def hierarchical_ensemble_hash(ensemble: HierarchicalPropulsionEnsemble) -> str:
    """Return the SHA-256 content identity of the complete propulsion ensemble.

    The digest is taken over the canonical JSON [1]_ of the document without
    its hash field. It is a reproducible identity, not an authentication.

    Returns
    -------
    str
        Lower-case hexadecimal digest.

    References
    ----------
    .. [1] IETF RFC 8259, "The JavaScript Object Notation (JSON) Data Interchange
           Format". https://www.rfc-editor.org/rfc/rfc8259
    """

    return hashlib.sha256(
        canonical_json_bytes(hierarchical_ensemble_to_dict(ensemble, include_hash=False))
    ).hexdigest()


def hierarchical_ensemble_to_json(
    ensemble: HierarchicalPropulsionEnsemble,
    *,
    indent: int | None = 2,
) -> str:
    """Serialize a hierarchical ensemble as JSON text with sorted keys.

    Non-ASCII text is kept as Unicode and non-finite numbers are rejected, so
    the output is standard JSON [1]_ [2]_.

    Parameters
    ----------
    ensemble : HierarchicalPropulsionEnsemble
        Ensemble to export.
    indent : int or None, optional
        Indentation passed to ``json.dumps``.

    Returns
    -------
    str
        JSON text of :func:`hierarchical_ensemble_to_dict` including the hash.

    References
    ----------
    .. [1] IETF RFC 8259, "The JavaScript Object Notation (JSON) Data Interchange
           Format". https://www.rfc-editor.org/rfc/rfc8259
    .. [2] PYTHON-3.12-JSON. Python 3.12, ``json``: encoder/decoder hooks,
           duplicate-member and nonstandard numeric token behavior.
    """

    return json.dumps(
        hierarchical_ensemble_to_dict(ensemble),
        ensure_ascii=False,
        allow_nan=False,
        sort_keys=True,
        indent=indent,
    )


def _close(left: float, right: float) -> bool:
    """Return whether two floats agree within ``1e-12`` relative or absolute."""
    return isclose(left, right, rel_tol=1e-12, abs_tol=1e-12)


def _required_bool(value: Any, field: str) -> bool:
    """Return ``value`` if it is a JSON boolean.

    Raises
    ------
    InputError
        If ``value`` is not a ``bool``.
    """
    if not isinstance(value, bool):
        raise InputError(f"Hierarchical ensemble field {field!r} must be boolean.")
    return value


def _validate_semantics(ensemble: HierarchicalPropulsionEnsemble) -> None:
    """Check that every derived quantity of a reopened ensemble reproduces.

    Verified are: method and counts; one definition per scenario and complete
    conditional indices; probabilities in ``[0, 1]`` with equal conditional
    weights ``1 / N_s`` and joint weights ``P_s / N_s``; exact closure of the
    joint weights to one within ``1e-12``; zero-thrust probability; weighted
    statistics [1]_ and the law-of-total-variance decomposition [2]_; scenario
    contributions; derived seeds and the replicate convergence flags; and
    unique reference identifiers.

    Raises
    ------
    InputError
        If any check fails.

    References
    ----------
    .. [1] NIST-WEIGHTED-MEAN. NIST Dataplot, "Weighted Mean", equation 2-20,
           normalized weighted mean.
    .. [2] MIT-TOTAL-VARIANCE. MIT RES.6-012, lecture 13.7, conditional and
           between-group variance identity (law of total variance).
    """
    if ensemble.method not in {"latin-hypercube", "monte-carlo"}:
        raise InputError("Unsupported hierarchical ensemble sampling method.")
    if ensemble.scenario_count != len(ensemble.definitions):
        raise InputError("Hierarchical ensemble scenario count does not reproduce.")
    if ensemble.realization_count != len(ensemble.realizations):
        raise InputError("Hierarchical ensemble realization count does not reproduce.")
    if not ensemble.definitions or not ensemble.realizations:
        raise InputError("Hierarchical ensemble must contain conditional realizations.")
    if not 0.0 <= ensemble.zero_thrust_probability <= 1.0:
        raise InputError("Hierarchical probability evidence has an invalid domain.")
    definition_map = {
        definition.scenario_id: definition for definition in ensemble.definitions
    }
    if len(definition_map) != len(ensemble.definitions):
        raise InputError("Hierarchical ensemble definitions are not unique.")
    scenario_ids = {item.scenario_id for item in ensemble.realizations}
    if scenario_ids != definition_map.keys():
        raise InputError("Hierarchical ensemble definition coverage does not reproduce.")

    for scenario_id, definition in definition_map.items():
        selected = tuple(
            item for item in ensemble.realizations if item.scenario_id == scenario_id
        )
        if len(selected) != definition.sample_count:
            raise InputError("Conditional realization count does not reproduce.")
        if tuple(sorted(item.conditional_index for item in selected)) != tuple(
            range(definition.sample_count)
        ):
            raise InputError("Conditional realization indices are incomplete.")
        probability = selected[0].scenario_probability
        if not isfinite(probability) or not 0.0 <= probability <= 1.0:
            raise InputError("Scenario probability must be finite and in [0, 1].")
        if any(not _close(item.scenario_probability, probability) for item in selected):
            raise InputError("Scenario probability changes inside a conditional ensemble.")
        if any(
            not isfinite(item.conditional_probability)
            or not 0.0 <= item.conditional_probability <= 1.0
            or not isfinite(item.joint_probability)
            or not 0.0 <= item.joint_probability <= 1.0
            for item in selected
        ):
            raise InputError("Conditional and joint probabilities must lie in [0, 1].")
        if any(
            not isfinite(value)
            for item in selected
            for _name, value in item.inputs
        ):
            raise InputError("Conditional realization inputs must be finite.")
        parameters = {parameter.name: parameter for parameter in definition.parameters}
        for item in selected:
            inputs = dict(item.inputs)
            if len(inputs) != len(item.inputs) or inputs.keys() != parameters.keys():
                raise InputError("Conditional input identities do not match the definition.", field="inputs")
            if any(not parameters[name].lower <= value <= parameters[name].upper for name, value in inputs.items()):
                raise InputError("Conditional input lies outside its declared bounds.", field="inputs")
            if item.conditional_probability != 1.0 / definition.sample_count:
                raise InputError("Conditional draws must have the declared equal weight.", field="conditional_probability")
        if not _close(sum(item.conditional_probability for item in selected), 1.0):
            raise InputError("Conditional probabilities do not close to one.")
        for item in selected:
            if not _close(
                item.joint_probability,
                item.scenario_probability * item.conditional_probability,
            ):
                raise InputError("A joint realization probability does not reproduce.")

    total_weight = sum((_fraction(item.joint_probability) for item in ensemble.realizations), Fraction(0))
    if abs(total_weight-1) > _fraction(1e-12) or ensemble.joint_probability_sum != _finite_rational(total_weight, "Joint probability sum"):
        raise InputError("Hierarchical joint probability does not reproduce.")
    expected_closure = _finite_rational(abs(total_weight-1), "Probability closure error")
    if ensemble.joint_probability_closure_error != expected_closure:
        raise InputError("Hierarchical probability closure evidence does not reproduce.")
    zero_probability = _finite_rational(sum((
        _fraction(item.joint_probability)
        for item in ensemble.realizations
        if item.artifact.delivered_total_impulse_n_s == 0.0
    ), Fraction(0))/total_weight, "Zero-thrust probability")
    if zero_probability != ensemble.zero_thrust_probability:
        raise InputError("Zero-thrust probability does not reproduce.")

    statistic_map = {item.metric: item for item in ensemble.statistics}
    decomposition_map = {
        item.metric: item for item in ensemble.variance_decomposition
    }
    expected_metrics = set(_METRIC_GETTERS)
    if len(statistic_map) != len(ensemble.statistics) or len(decomposition_map) != len(ensemble.variance_decomposition):
        raise InputError("Hierarchical metric identities must be unique.", field="metric")
    if set(statistic_map) != expected_metrics or set(decomposition_map) != expected_metrics:
        raise InputError("Hierarchical metric coverage does not reproduce.")
    for metric, getter in _METRIC_GETTERS.items():
        weighted_values = tuple(
            (getter(item.artifact), item.joint_probability)
            for item in ensemble.realizations
        )
        expected_statistic = _weighted_statistics(metric, weighted_values)
        actual_statistic = statistic_map[metric]
        if any(
            getattr(actual_statistic, field) != getattr(expected_statistic, field)
            for field in (
                "expected_value",
                "standard_deviation",
                "minimum",
                "p05",
                "p50",
                "p95",
                "maximum",
            )
        ):
            raise InputError(f"Hierarchical statistic {metric!r} does not reproduce.")
        expected_decomposition = _variance_decomposition(
            metric, ensemble.realizations, getter
        )
        actual_decomposition = decomposition_map[metric]
        if any(getattr(actual_decomposition, name) < 0.0 for name in
               ("total_variance", "within_scenario_variance", "between_scenario_variance", "closure_error")):
            raise InputError("Hierarchical variance evidence must be nonnegative.", field="variance_decomposition")
        if any(getattr(actual_decomposition, name) < 0.0 for name in ("within_fraction", "between_fraction")):
            raise InputError("Hierarchical variance fractions must be nonnegative.", field="variance_decomposition")
        if any(
            getattr(actual_decomposition, field) != getattr(expected_decomposition, field)
            for field in (
                "total_variance",
                "within_scenario_variance",
                "between_scenario_variance",
                "closure_error",
                "within_fraction",
                "between_fraction",
            )
        ):
            raise InputError(
                f"Hierarchical variance decomposition {metric!r} does not reproduce."
            )

    contribution_map = {
        (item.scenario_id, item.metric): item for item in ensemble.contributions
    }
    if len(contribution_map) != len(ensemble.contributions):
        raise InputError("Hierarchical contribution identities must be unique.", field="contributions")
    if len(contribution_map) != len(definition_map) * len(expected_metrics):
        raise InputError("Hierarchical scenario contribution coverage is incomplete.")
    for scenario_id in definition_map:
        selected = tuple(
            item for item in ensemble.realizations if item.scenario_id == scenario_id
        )
        scenario_probability = selected[0].scenario_probability
        for metric, getter in _METRIC_GETTERS.items():
            actual = contribution_map.get((scenario_id, metric))
            if actual is None:
                raise InputError("Hierarchical scenario contribution is missing.")
            group = _weight_ledger((getter(item.artifact), item.joint_probability) for item in selected)
            total = _weight_ledger((0.0, item.joint_probability) for item in ensemble.realizations).total
            conditional_mean = _finite_rational(group.mean, "Conditional mean")
            if not all(
                (
                    actual.scenario_probability == scenario_probability,
                    actual.conditional_expected_value == conditional_mean,
                    actual.joint_expected_value_contribution == _finite_rational(group.mean*group.total/total, "Scenario mean contribution"),
                )
            ):
                raise InputError("Hierarchical scenario contribution does not reproduce.")

    evidence_map = {item.scenario_id: item for item in ensemble.sampling_evidence}
    if len(evidence_map) != len(ensemble.sampling_evidence):
        raise InputError("Conditional evidence identities must be unique.", field="sampling_evidence")
    if set(evidence_map) != set(definition_map):
        raise InputError("Conditional sampling evidence coverage is incomplete.")
    for scenario_id, definition in definition_map.items():
        evidence = evidence_map[scenario_id]
        differences = (
            evidence.impulse_mean_relative_difference,
            evidence.impulse_standard_deviation_relative_difference,
        )
        if evidence.sample_count != definition.sample_count:
            raise InputError("Conditional sampling evidence count does not reproduce.")
        if evidence.primary_seed == evidence.audit_seed:
            raise InputError("Conditional primary and audit seeds must differ.")
        if (evidence.primary_seed != _derived_seed(ensemble.base_seed, scenario_id, "primary")
                or evidence.audit_seed != _derived_seed(ensemble.base_seed, scenario_id, "audit")):
            raise InputError("Conditional seeds do not match the base seed and scenario.", field="sampling_evidence")
        if (
            not isfinite(evidence.requested_tolerance)
            or evidence.requested_tolerance <= 0.0
        ):
            raise InputError("Conditional sampling tolerance must be positive and finite.")
        if any(not isfinite(value) or value < 0.0 for value in differences):
            raise InputError("Conditional sampling differences must be finite and nonnegative.")
        if evidence.converged != (
            max(differences) <= evidence.requested_tolerance
        ):
            raise InputError("Conditional sampling convergence flag does not reproduce.")
    if len(ensemble.reference_ids) != len(set(ensemble.reference_ids)):
        raise InputError("Hierarchical ensemble reference identifiers must be unique.")


def hierarchical_ensemble_from_dict(
    payload: dict[str, Any],
) -> HierarchicalPropulsionEnsemble:
    """Reopen an ensemble document and verify hashes and derived quantities.

    The closed document schema and numerical-policy identifiers are checked
    first, nested artifacts are reopened with their own integrity checks, every
    derived quantity is recomputed, and the ensemble hash must match [1]_ [2]_.

    Parameters
    ----------
    payload : dict
        Document produced by :func:`hierarchical_ensemble_to_dict`.

    Returns
    -------
    HierarchicalPropulsionEnsemble
        The validated ensemble.

    Raises
    ------
    InputError
        If the document is malformed, uses another schema or numerical policy,
        any derived quantity does not reproduce, or the hash does not match.

    References
    ----------
    .. [1] NASA-STD-7009B, "Standard for Models and Simulations".
           https://standards.nasa.gov/standard/NASA/NASA-STD-7009
    .. [2] NASA/SP-2011-3421, "Probabilistic Risk Assessment Procedures Guide for
           NASA Managers and Practitioners". https://ntrs.nasa.gov/citations/20120001369
    """

    _ensemble_document(payload)
    _require_numerical_policies(payload, {"statistics_numerics": _PROBABILITY_POLICY})
    _require_numerical_policies(payload, {"sampling_numerics": _SAMPLING_POLICY})
    if payload.get("schema") != HIERARCHICAL_ENSEMBLE_SCHEMA:
        raise InputError("Unsupported hierarchical ensemble schema.")
    try:
        ensemble = HierarchicalPropulsionEnsemble(
            method=str(payload["method"]),
            base_seed=int(payload["base_seed"]),
            scenario_count=int(payload["scenario_count"]),
            realization_count=int(payload["realization_count"]),
            definitions=tuple(
                _definition_from_dict(item) for item in payload["definitions"]
            ),
            realizations=tuple(
                ConditionalPropulsionRealization(
                    scenario_id=str(item["scenario_id"]),
                    scenario_probability=float(item["scenario_probability"]),
                    conditional_index=int(item["conditional_index"]),
                    conditional_probability=float(item["conditional_probability"]),
                    joint_probability=float(item["joint_probability"]),
                    inputs=tuple(
                        (str(pair[0]), float(pair[1])) for pair in item["inputs"]
                    ),
                    artifact=TabulatedBurnArtifact.from_dict(item["artifact"]),
                )
                for item in payload["realizations"]
            ),
            statistics=tuple(
                HierarchicalMetricStatistics(
                    metric=str(item["metric"]),
                    expected_value=float(item["expected_value"]),
                    standard_deviation=float(item["standard_deviation"]),
                    minimum=float(item["minimum"]),
                    p05=float(item["p05"]),
                    p50=float(item["p50"]),
                    p95=float(item["p95"]),
                    maximum=float(item["maximum"]),
                )
                for item in payload["statistics"]
            ),
            variance_decomposition=tuple(
                HierarchicalVarianceDecomposition(
                    metric=str(item["metric"]),
                    total_variance=float(item["total_variance"]),
                    within_scenario_variance=float(
                        item["within_scenario_variance"]
                    ),
                    between_scenario_variance=float(
                        item["between_scenario_variance"]
                    ),
                    closure_error=float(item["closure_error"]),
                    within_fraction=float(item["within_fraction"]),
                    between_fraction=float(item["between_fraction"]),
                )
                for item in payload["variance_decomposition"]
            ),
            contributions=tuple(
                ScenarioMetricContribution(
                    scenario_id=str(item["scenario_id"]),
                    metric=str(item["metric"]),
                    scenario_probability=float(item["scenario_probability"]),
                    conditional_expected_value=float(
                        item["conditional_expected_value"]
                    ),
                    joint_expected_value_contribution=float(
                        item["joint_expected_value_contribution"]
                    ),
                )
                for item in payload["contributions"]
            ),
            sampling_evidence=tuple(
                ConditionalSamplingEvidence(
                    scenario_id=str(item["scenario_id"]),
                    sample_count=int(item["sample_count"]),
                    primary_seed=int(item["primary_seed"]),
                    audit_seed=int(item["audit_seed"]),
                    impulse_mean_relative_difference=float(
                        item["impulse_mean_relative_difference"]
                    ),
                    impulse_standard_deviation_relative_difference=float(
                        item["impulse_standard_deviation_relative_difference"]
                    ),
                    requested_tolerance=float(item["requested_tolerance"]),
                    converged=_required_bool(item["converged"], "converged"),
                )
                for item in payload["sampling_evidence"]
            ),
            zero_thrust_probability=float(payload["zero_thrust_probability"]),
            joint_probability_sum=float(payload["joint_probability_sum"]),
            joint_probability_closure_error=float(
                payload["joint_probability_closure_error"]
            ),
            reference_ids=tuple(str(item) for item in payload["reference_ids"]),
            warnings=tuple(str(item) for item in payload["warnings"]),
        )
    except InputError:
        raise
    except (KeyError, TypeError, ValueError, IndexError) as error:
        raise InputError("Malformed hierarchical propulsion ensemble.") from error
    try:
        _validate_semantics(ensemble)
    except InputError:
        raise
    except DomainError as error:
        raise InputError(error.message, field=error.field) from error
    except (TypeError, ValueError, OverflowError, ZeroDivisionError) as error:
        raise InputError("Hierarchical semantic calculation could not be represented.") from error
    if payload.get("ensemble_hash") != hierarchical_ensemble_hash(ensemble):
        raise InputError("Hierarchical ensemble hash does not match its contents.")
    return ensemble


def hierarchical_ensemble_from_json(text: str) -> HierarchicalPropulsionEnsemble:
    """Parse strict JSON text [1]_ [2]_ and verify the hierarchical ensemble.

    Parameters
    ----------
    text : str
        JSON document text.

    Returns
    -------
    HierarchicalPropulsionEnsemble
        The validated ensemble.

    Raises
    ------
    InputError
        If the text is not valid JSON or fails
        :func:`hierarchical_ensemble_from_dict`.

    References
    ----------
    .. [1] IETF RFC 8259, "The JavaScript Object Notation (JSON) Data Interchange
           Format". https://www.rfc-editor.org/rfc/rfc8259
    .. [2] PYTHON-3.12-JSON. Python 3.12, ``json``: encoder/decoder hooks,
           duplicate-member and nonstandard numeric token behavior.
    """

    return read_json_document(text, "Hierarchical ensemble", hierarchical_ensemble_from_dict)

