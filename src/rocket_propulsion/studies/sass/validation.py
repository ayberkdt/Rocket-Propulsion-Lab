"""Conservative comparison of independent measurements and model predictions.

This module verifies supplied evidence and compares bounded intervals. It does
not authenticate the origin of files or certify a physical system. Acceptance
limits and uncertainty bounds must be established outside this implementation.

References
----------
.. [1] NASA-7009B. NASA, "Standard for Models and Simulations, NASA-STD-7009B";
       4.2.5-4.2.7 empirical referents, validation domain and uncertainty.
.. [2] NASA-TRL. NASA ESTO, "Technology Readiness Levels"; hardware demonstration is
       separate from model agreement.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass
from fractions import Fraction
from hashlib import sha256
from json import dumps
from pathlib import Path
from sys import float_info
from typing import cast

from rocket_propulsion.core.errors import DomainError
from rocket_propulsion.core.validation import (
    integer,
    read_json_document,
    real,
    records,
    strings,
    text,
)

from .evidence import CitedScalar, quantity
from .native_trace import NativeCsvTrace, native_csv_trace_status
from .rational_numerics import outward_endpoint
from .shutdown_observation import ShutdownWindow
from .transient_validation import (
    ArchivedTransientRun,
    PressureReferent,
    ShutdownValidationTarget,
    TransientMeasurement,
    TransientValidationTarget,
    _pressure_referent_from_payload,
    shutdown_measurement,
    shutdown_valve_timing,
    transient_archive_bytes,
    transient_measurement,
)
from .valve_response import ValveResponseReview


@dataclass(frozen=True, slots=True)
class ValidationSample:
    """One holdout measurement and prediction under a documented common condition.

    Bounds are absolute error bounds in the measurand's units, not standard deviations
    or confidence intervals; None means unknown, never exact zero. Domain membership
    requires external review of the condition record; ``in_validation_domain`` records
    that review and is not inferred from a numeric range [1]_.

    Attributes
    ----------
    identifier : str
        Unique sample identifier.
    condition_record : str
        Joint operating-condition record.
    measured, predicted : CitedScalar
        Measured referent and model prediction in the same unit.
    measurement_bound, prediction_bound : CitedScalar or None
        Absolute error bounds of measurement (positive) and prediction (non-negative).
    allowable_error : CitedScalar
        Approved absolute error allowance; non-negative.
    in_validation_domain : bool
        Outcome of the external domain review.

    Raises
    ------
    DomainError
        If units differ, a bound is negative, or the domain review is not a bool.

    References
    ----------
    .. [1] NASA-7009B. NASA, "Standard for Models and Simulations, NASA-STD-7009B";
           4.2.5 validation referents/domain; 4.2.7 uncertainty rationale.
    """

    identifier: str
    condition_record: str
    measured: CitedScalar
    predicted: CitedScalar
    measurement_bound: CitedScalar | None
    prediction_bound: CitedScalar | None
    allowable_error: CitedScalar
    in_validation_domain: bool

    def __post_init__(self) -> None:
        """Reject mismatched units, negative bounds and an implicit domain review.

        References
        ----------
        .. [1] NASA-7009B. NASA, "Standard for Models and Simulations, NASA-STD-7009B";
               4.2.1 units; 4.2.5/4.2.7 validation and uncertainty records.
        """
        text(self.identifier, "validation sample identifier")
        text(self.condition_record, "joint operating-condition record")
        if not isinstance(self.measured, CitedScalar):
            raise DomainError("A measurement must retain its source evidence.")
        unit = self.measured.unit
        quantity(self.predicted, unit)
        quantity(self.allowable_error, unit, minimum=0)
        if self.measurement_bound is not None:
            quantity(self.measurement_bound, unit, strictly_positive=True)
        if self.prediction_bound is not None:
            quantity(self.prediction_bound, unit, minimum=0)
        if not isinstance(self.in_validation_domain, bool):
            raise DomainError("Domain review must be an explicit boolean.")


@dataclass(frozen=True, slots=True)
class ValidationEvidence:
    """A reproducible, independent holdout comparison with explicit pedigree.

    The raw-file digest binds the archived bytes. A reviewed plan can also
    bind the transcribed values to that archive and the declared model.
    A transcription-review record is mandatory. Calibration
    sample identifiers must not overlap holdout identifiers. Dataset authorship,
    applicability and sensor calibration are human-reviewed records [1]_.

    Attributes
    ----------
    test_article_id, model_configuration_id : str
        Identifiers of the tested hardware and the model configuration.
    raw_data_path : pathlib.Path
        Archived raw data file.
    raw_data_sha256, model_sha256 : str
        Lowercase hexadecimal SHA-256 digests.
    calibration_sample_ids : tuple of str
        Samples used for calibration.
    samples : tuple of ValidationSample
        Non-empty, uniquely identified holdout samples.
    independence_record, transcription_review_record : str
        Review records.
    uncertainty_basis, acceptance_approval_record : str
        Basis of the bounds and approval of the allowances.

    Raises
    ------
    DomainError
        If a record is empty, a digest is malformed, the path is not a ``Path``, or
        samples are missing, duplicated or overlap calibration.

    References
    ----------
    .. [1] NASA-7009B. NASA, "Standard for Models and Simulations, NASA-STD-7009B";
           4.2.1 data pedigree, 4.2.5 empirical data and 4.2.7 uncertainty.
    """

    test_article_id: str
    model_configuration_id: str
    raw_data_path: Path
    raw_data_sha256: str
    model_sha256: str
    calibration_sample_ids: tuple[str, ...]
    samples: tuple[ValidationSample, ...]
    independence_record: str
    transcription_review_record: str
    uncertainty_basis: str
    acceptance_approval_record: str

    def __post_init__(self) -> None:
        """Freeze sample identities and require traceable review records.

        References
        ----------
        .. [1] NASA-7009B. NASA, "Standard for Models and Simulations, NASA-STD-7009B";
               4.1.1.5 project acceptance criteria; 4.2.5 data referents.
        """
        for field in ("test_article_id", "model_configuration_id", "independence_record",
                      "transcription_review_record", "uncertainty_basis", "acceptance_approval_record"):
            text(getattr(self, field), field)
        for digest in (self.raw_data_sha256, self.model_sha256):
            if not isinstance(digest, str) or len(digest) != 64 or any(
                char not in "0123456789abcdef" for char in digest
            ):
                raise DomainError("Evidence digests must be lowercase SHA-256 hex strings.")
        if not isinstance(self.raw_data_path, Path):
            raise DomainError("Raw data path must be an explicit pathlib.Path.")
        samples = records(self.samples, ValidationSample, "validation samples")
        identifiers = {sample.identifier for sample in samples}
        calibration = strings(self.calibration_sample_ids, "calibration sample IDs")
        if not samples or len(identifiers) != len(samples):
            raise DomainError("Holdout samples must be non-empty and uniquely identified.")
        if identifiers.intersection(calibration):
            raise DomainError("Calibration and validation holdout samples overlap.")
        object.__setattr__(self, "samples", samples)
        object.__setattr__(self, "calibration_sample_ids", calibration)


@dataclass(frozen=True, slots=True)
class PhysicalAcceptanceCriteria:
    """Approved inclusive bounds on a replayed measurement and optional valve timing.

    References
    ----------
    .. [1] NASA-7009B, 4.1.1.5 and 4.2.7, approved acceptance criteria and uncertainty.
    """

    approval_record: str
    minimum: CitedScalar | None = None
    maximum: CitedScalar | None = None
    require_valve_timing: bool = False

    def __post_init__(self) -> None:
        """Require an approved criterion, coherent units and non-assumed limits.

        References
        ----------
        .. [1] NASA-7009B, 4.1.1.5 and 4.2.1, criteria and input traceability.
        """
        text(self.approval_record, "physical criteria approval")
        if not isinstance(self.require_valve_timing, bool):
            raise DomainError("The valve timing requirement must be an explicit boolean.")
        if self.minimum is None and self.maximum is None and not self.require_valve_timing:
            raise DomainError("Physical acceptance requires a measurement limit or valve timing.")
        for limit in (self.minimum, self.maximum):
            if limit is not None and (not isinstance(limit, CitedScalar) or limit.assumption_notes):
                raise DomainError("Physical acceptance limits require cited values without assumptions.")
        if (self.minimum is not None and self.maximum is not None
                and (self.minimum.unit != self.maximum.unit or self.minimum.value > self.maximum.value)):
            raise DomainError("Physical acceptance bounds must share units and be ordered.")


@dataclass(frozen=True, slots=True)
class ValidationTarget:
    """One required holdout identity, operating condition and approved allowance.

    The unit on `allowable_error` also fixes the required measurand's unit.
    Sample identifiers should identify the observation, including its sensor
    and event where applicable. Condition records should identify the whole
    test run, shared by observations made during that run. Set the reviewed
    digest after checking the transcription and its evidence context.

    References
    ----------
    .. [1] NASA-7009B, 4.1.1.5 and 4.2.5, acceptance criteria and empirical referents.
    """

    sample_id: str
    condition_record: str
    allowable_error: CitedScalar
    reviewed_sample_sha256: str | None = None
    transient_target: TransientValidationTarget | ShutdownValidationTarget | None = None
    physical_criteria: PhysicalAcceptanceCriteria | None = None

    def __post_init__(self) -> None:
        """Require a cited, nonnegative allowance and stable target identity.

        References
        ----------
        .. [1] NASA-7009B, 4.1.1.5 and 4.2.1, acceptance criteria and units.
        """
        text(self.sample_id, "planned holdout sample ID")
        text(self.condition_record, "planned operating-condition record")
        if not isinstance(self.allowable_error, CitedScalar):
            raise DomainError("A planned allowance requires cited units and approval.")
        quantity(self.allowable_error, self.allowable_error.unit, minimum=0)
        digest = self.reviewed_sample_sha256
        if digest is not None and (not isinstance(digest, str) or len(digest) != 64
                                   or any(char not in "0123456789abcdef" for char in digest)):
            raise DomainError("Reviewed sample digest must be lowercase SHA-256 hex.")
        if self.transient_target is not None and not isinstance(
            self.transient_target, (TransientValidationTarget, ShutdownValidationTarget)
        ):
            raise DomainError("A planned replay recipe must be a transient or shutdown target.")
        _check_physical_criteria(self.physical_criteria, self.transient_target, self.allowable_error.unit)


def _check_physical_criteria(
    criteria: PhysicalAcceptanceCriteria | None,
    recipe: TransientValidationTarget | ShutdownValidationTarget | None, unit: str,
) -> None:
    """Check the criterion's scope against the measurement and replay method.

    References
    ----------
    .. [1] NASA-7009B, 4.1.1.5 and 4.2.5, applicable acceptance criteria and referents.
    """
    if criteria is None:
        return
    if not isinstance(criteria, PhysicalAcceptanceCriteria) or recipe is None:
        raise DomainError("Physical acceptance requires approved criteria and acquisition replay.")
    if criteria.require_valve_timing and not isinstance(recipe, ShutdownValidationTarget):
        raise DomainError("Required valve timing needs a shutdown replay target.")
    for limit in (criteria.minimum, criteria.maximum):
        if limit is not None:
            quantity(limit, unit)


@dataclass(frozen=True, slots=True, order=True)
class NativeAcquisitionSlice:
    """Reviewed portion of one native acquisition, indexed before export or resampling.

    References
    ----------
    .. [1] NASA-7009B, 4.2.1 and 4.2.5, pedigree and identity of empirical data.
    """

    condition_record: str
    firing_id: str
    native_sha256: str
    first_index: int
    past_end_index: int

    def __post_init__(self) -> None:
        """Require a traceable firing identity and nonempty half-open native index range.

        References
        ----------
        .. [1] NASA-7009B, 4.2.1, reviewed input representation.
        """
        text(self.condition_record, "native acquisition condition")
        text(self.firing_id, "native acquisition firing ID")
        if not isinstance(self.native_sha256, str) or len(self.native_sha256) != 64 or any(
            char not in "0123456789abcdef" for char in self.native_sha256
        ):
            raise DomainError("A native acquisition digest must be lowercase SHA-256 hex.")
        first = integer(self.first_index, "first native acquisition index")
        end = integer(self.past_end_index, "past-end native acquisition index")
        if not 0 <= first < end:
            raise DomainError("A native acquisition slice requires 0 <= first index < past-end index.")
        object.__setattr__(self, "first_index", first)
        object.__setattr__(self, "past_end_index", end)


@dataclass(frozen=True, slots=True)
class NativeAcquisitionFile:
    """Local location of a content-addressed original acquisition file.

    References
    ----------
    .. [1] NASA-7009B, 4.2.1 and 4.2.5, empirical input pedigree.
    """

    native_sha256: str
    path: Path

    def __post_init__(self) -> None:
        """Validate the digest and path without opening the file during planning.

        References
        ----------
        .. [1] NASA-7009B, 4.2.1, reviewed input identities.
        """
        if not isinstance(self.native_sha256, str) or len(self.native_sha256) != 64 or any(
            char not in "0123456789abcdef" for char in self.native_sha256
        ):
            raise DomainError("A native acquisition digest must be lowercase SHA-256 hex.")
        if not isinstance(self.path, Path):
            raise DomainError("A native acquisition file requires an explicit pathlib.Path.")


@dataclass(frozen=True, slots=True)
class AcquisitionPartition:
    """Reviewed map from calibration and holdout conditions to native acquisition slices.

    References
    ----------
    .. [1] NASA-7009B, 4.2.1 and 4.2.5, reviewed acquisition and use.
    """

    review_record: str
    slices: tuple[NativeAcquisitionSlice, ...]
    calibration_samples: tuple[tuple[str, str], ...] = ()
    native_csv_traces: tuple[NativeCsvTrace, ...] = ()

    def __post_init__(self) -> None:
        """Freeze the reviewed mapping in canonical order.

        References
        ----------
        .. [1] NASA-7009B, 4.2.1, deterministic input records.
        """
        text(self.review_record, "native acquisition partition review")
        slices = records(self.slices, NativeAcquisitionSlice, "native acquisition slices")
        if not slices or len(set(slices)) != len(slices):
            raise DomainError("An acquisition partition requires nonempty, unique native slices.")
        samples = records(self.calibration_samples, tuple, "calibration sample-to-condition mappings")
        sample_ids: set[str] = set()
        normalized = []
        for mapping in samples:
            if len(mapping) != 2:
                raise DomainError("A calibration mapping needs a sample ID and condition record.")
            identifier = text(mapping[0], "calibration sample mapping ID")
            condition = text(mapping[1], "calibration sample mapping condition")
            if identifier in sample_ids:
                raise DomainError("A calibration sample requires a unique acquisition mapping.")
            sample_ids.add(identifier)
            normalized.append((identifier, condition))
        object.__setattr__(self, "slices", tuple(sorted(slices)))
        object.__setattr__(self, "calibration_samples", tuple(sorted(normalized)))
        conversions: tuple[NativeCsvTrace, ...] = records(self.native_csv_traces, NativeCsvTrace, "native CSV mappings")
        if len({(row.condition_record, row.channel_id) for row in conversions}) != len(conversions):
            raise DomainError("A native CSV mapping requires a unique condition and channel.")
        object.__setattr__(self, "native_csv_traces", tuple(sorted(conversions)))


@dataclass(frozen=True, slots=True)
class ValidationPlan:
    """Reviewed holdout coverage for a named test article.

    Prepare the target list and allowances independently of the submitted
    results. List the whole conditions used for calibration, even when a
    holdout observation has a different sensor ID. The review record identifies
    the external approval. When supplied, `reviewed_predictions_sha256`
    requires the archived model and prediction export for the comparison.

    References
    ----------
    .. [1] NASA-7009B, 4.1.1.5 and 4.2.5, acceptance criteria and validation.
    """

    test_article_id: str
    review_record: str
    targets: tuple[ValidationTarget, ...]
    calibration_condition_records: tuple[str, ...]
    reviewed_predictions_sha256: str | None = None
    acquisition_partition: AcquisitionPartition | None = None
    reviewed_acquisition_partition_sha256: str | None = None
    native_acquisitions: tuple[NativeAcquisitionFile, ...] = ()

    def __post_init__(self) -> None:
        """Freeze the required targets and reject calibration-run overlap.

        References
        ----------
        .. [1] NASA-7009B, 4.2.5, empirical validation referents and domain.
        """
        text(self.test_article_id, "planned test article ID")
        text(self.review_record, "holdout plan review record")
        targets = records(self.targets, ValidationTarget, "planned holdout targets")
        if not targets or len({target.sample_id for target in targets}) != len(targets):
            raise DomainError("Planned holdout targets must be non-empty and uniquely identified.")
        calibration = strings(self.calibration_condition_records, "calibration condition records")
        if {target.condition_record for target in targets}.intersection(calibration):
            raise DomainError("A planned holdout condition was used for calibration.")
        digest = self.reviewed_predictions_sha256
        if digest is not None and (not isinstance(digest, str) or len(digest) != 64
                                   or any(char not in "0123456789abcdef" for char in digest)):
            raise DomainError("Reviewed prediction archive digest must be lowercase SHA-256 hex.")
        if self.acquisition_partition is not None and not isinstance(self.acquisition_partition, AcquisitionPartition):
            raise DomainError("Expected a typed acquisition partition.")
        partition_digest = self.reviewed_acquisition_partition_sha256
        if partition_digest is not None and (not isinstance(partition_digest, str) or len(partition_digest) != 64
                                             or any(char not in "0123456789abcdef" for char in partition_digest)):
            raise DomainError("Reviewed acquisition partition digest must be lowercase SHA-256 hex.")
        sources: tuple[NativeAcquisitionFile, ...] = records(
            self.native_acquisitions, NativeAcquisitionFile, "native acquisition files"
        )
        if len({source.native_sha256 for source in sources}) != len(sources):
            raise DomainError("Each native acquisition digest requires a unique file binding.")
        object.__setattr__(self, "targets", targets)
        object.__setattr__(self, "calibration_condition_records", calibration)
        object.__setattr__(self, "native_acquisitions", sources)


def acquisition_partition_sha256(evidence: ValidationEvidence, plan: ValidationPlan) -> str:
    """Identify the mapping and its evidence context for independent approval.

    References
    ----------
    .. [1] NASA-7009B, 4.2.1 and 4.2.5, traceability of validation datasets.
    """
    if (not isinstance(evidence, ValidationEvidence) or not isinstance(plan, ValidationPlan)
            or evidence.test_article_id != plan.test_article_id or plan.acquisition_partition is None):
        raise DomainError("An acquisition partition digest requires a plan and evidence for the same test article.")
    partition = asdict(plan.acquisition_partition)
    if not plan.acquisition_partition.native_csv_traces:
        partition.pop("native_csv_traces")
    for mapping in partition.get("native_csv_traces", ()):
        if mapping["absolute_error_column"] is None:
            mapping.pop("absolute_error_column")
        if mapping["tick_period_s"] is None:
            mapping.pop("tick_period_s")
    payload = {
        "schema": "rocket-propulsion-acquisition-partition-v1",
        "test_article_id": evidence.test_article_id,
        "raw_data_sha256": evidence.raw_data_sha256,
        "model_configuration_id": evidence.model_configuration_id,
        "model_sha256": evidence.model_sha256,
        "independence_record": evidence.independence_record,
        "plan_review_record": plan.review_record,
        "calibration_sample_ids": sorted(set(evidence.calibration_sample_ids)),
        "calibration_conditions": sorted(set(plan.calibration_condition_records)),
        "holdout_targets": sorted((target.sample_id, target.condition_record) for target in plan.targets),
        "acquisition_partition": partition,
    }
    return sha256(dumps(payload, sort_keys=True, ensure_ascii=True, separators=(",", ":")).encode("ascii")).hexdigest()


@dataclass(frozen=True, slots=True)
class ArchivedModelRun:
    """Paths to the archived model configuration and its prediction export.

    The model file may be a packaged configuration with its dependencies.
    Assessment checks its bytes and reads the prediction JSON as data.

    References
    ----------
    .. [1] NASA-7009B, 4.2.1, data pedigree and model configuration.
    """

    model_path: Path
    predictions_path: Path

    def __post_init__(self) -> None:
        """Require explicit paths for the two evidence files.

        References
        ----------
        .. [1] NASA-7009B, 4.2.1, traceable inputs and configuration.
        """
        if not isinstance(self.model_path, Path) or not isinstance(self.predictions_path, Path):
            raise DomainError("Archived model and prediction paths must be explicit pathlib.Path values.")


def validation_sample_sha256(
    evidence: ValidationEvidence, sample: ValidationSample, *,
    transient_target: TransientValidationTarget | ShutdownValidationTarget | None = None,
    physical_criteria: PhysicalAcceptanceCriteria | None = None,
) -> str:
    """Identify the reviewed numeric row and the evidence it refers to.

    After inspecting the source-to-row transcription, record this digest in
    the separately controlled plan target. It is an identity for the reviewed
    content. The function performs no file reads and grants no approval.
    The canonical payload uses sorted JSON keys, ASCII escapes, compact
    separators, exact binary64 hex values, and a versioned schema. A replayed
    transient row also binds the reviewed channel, window and statistic.
    Physical criteria, when present, are bound with their approval record.

    References
    ----------
    .. [1] NASA-7009B, 4.2.1 and 4.2.5, data pedigree and empirical referents.
    .. [2] NIST FIPS 180-4, SHA-256, https://csrc.nist.gov/pubs/fips/180-4/upd1/final
    """
    if not isinstance(evidence, ValidationEvidence) or not isinstance(sample, ValidationSample):
        raise DomainError("A sample digest requires validation evidence and a validation sample.")
    if transient_target is not None and not isinstance(
        transient_target, (TransientValidationTarget, ShutdownValidationTarget)
    ):
        raise DomainError("A reviewed replay digest requires a transient or shutdown target.")
    _check_physical_criteria(physical_criteria, transient_target, sample.measured.unit)
    sample_data = asdict(sample)
    for field in ("measured", "predicted", "measurement_bound", "prediction_bound", "allowable_error"):
        scalar = sample_data[field]
        if scalar is not None:
            scalar["value"] = scalar["value"].hex()
    payload: dict[str, object] = {
        "schema": "rocket-propulsion-validation-sample-v1",
        "test_article_id": evidence.test_article_id,
        "raw_data_sha256": evidence.raw_data_sha256,
        "model_configuration_id": evidence.model_configuration_id,
        "model_sha256": evidence.model_sha256,
        "calibration_sample_ids": sorted(set(evidence.calibration_sample_ids)),
        "independence_record": evidence.independence_record,
        "transcription_review_record": evidence.transcription_review_record,
        "uncertainty_basis": evidence.uncertainty_basis,
        "acceptance_approval_record": evidence.acceptance_approval_record,
        "sample": sample_data,
    }
    if transient_target is not None:
        payload["schema"] = ("rocket-propulsion-validation-sample-v3"
                             if isinstance(transient_target, ShutdownValidationTarget)
                             else "rocket-propulsion-validation-sample-v2")
        payload["transient_target"] = asdict(transient_target)
        if isinstance(transient_target, TransientValidationTarget) and transient_target.pressure_referent is None:
            cast(dict, payload["transient_target"]).pop("pressure_referent", None)
    if physical_criteria is not None:
        payload["schema"] = "rocket-propulsion-validation-sample-v4"
        limits = asdict(physical_criteria)
        for endpoint in ("minimum", "maximum"):
            if limits[endpoint] is not None:
                limits[endpoint]["value"] = limits[endpoint]["value"].hex()
        payload["physical_criteria"] = limits
    return sha256(dumps(payload, sort_keys=True, ensure_ascii=True, separators=(",", ":"),
                        allow_nan=False).encode("ascii")).hexdigest()


def compare_validation_sample(sample: ValidationSample) -> dict[str, object]:
    """Bound the prediction error with the supplied error intervals and classify it.

    Parameters
    ----------
    sample : ValidationSample
        Reviewed holdout sample.

    Returns
    -------
    dict
        ``status`` is ``'unassessed'`` (with a reason) outside the reviewed domain or
        without both bounds; otherwise ``'pass'``, ``'fail'`` or ``'indeterminate'``
        with the discrepancy, minimum and maximum error, allowance and numerical
        guard.

    Raises
    ------
    DomainError
        If ``sample`` is not a ``ValidationSample`` or an interval endpoint is not
        finite.

    Notes
    -----
    max error = |prediction - measurement| + measurement_bound + prediction_bound and
    min error = max(0, |prediction - measurement| - both bounds), by the triangle
    inequality on like-unit intervals [2]_. Pass requires the whole interval inside
    the approved allowance; a partial overlap is indeterminate, with no Gaussian
    assumption or fitted tolerance [1]_. This conservative interval rule is the
    method defined by this module, not a NASA formula. Exact rational intervals are
    rounded outward for publication and decisions use those bounds plus an
    outward-rounded guard of 16 machine epsilons of the reported scale.

    References
    ----------
    .. [1] NASA-7009B. NASA, "Standard for Models and Simulations, NASA-STD-7009B";
           4.2.5 and 4.2.7 require justified metrics and uncertainty.
    .. [2] BIPM-SI. BIPM, "SI Brochure, 9th edition"; like-unit interval arithmetic;
           triangle inequality applied in this module.
    """
    if not isinstance(sample, ValidationSample):
        raise DomainError("Expected a reviewed validation sample.")
    result = {"sample_id": sample.identifier, "unit": sample.measured.unit,
              "condition_record": sample.condition_record, "status": "unassessed"}
    if not sample.in_validation_domain:
        return {**result, "reason": "outside_reviewed_validation_domain"}
    if sample.measurement_bound is None or sample.prediction_bound is None:
        return {**result, "reason": "missing_uncertainty_bound"}
    discrepancy = abs(Fraction(sample.predicted.value) - Fraction(sample.measured.value))
    bound = Fraction(sample.measurement_bound.value) + Fraction(sample.prediction_bound.value)
    high = outward_endpoint(discrepancy + bound, upward=True)
    low = outward_endpoint(max(Fraction(0), discrepancy - bound), upward=False)
    allowance = Fraction(sample.allowable_error.value)
    # Expose a near-boundary result for review rather than hiding roundoff in
    # a pass. Publish an upper bound for 16 epsilons of the reported scale.
    guard_value = 16 * Fraction(float_info.epsilon) * max(Fraction(high), allowance)
    guard = outward_endpoint(guard_value, upward=True)
    status = ("pass" if Fraction(high) + Fraction(guard) < allowance else
              "fail" if Fraction(low) > allowance + Fraction(guard) else "indeterminate")
    return {**result, "status": status, "discrepancy": float(discrepancy),
            "minimum_error": low, "maximum_error": high,
            "allowable_error": sample.allowable_error.value, "numerical_guard": guard,
            "numerical_method": "exact_intervals_with_outward_binary_endpoints"}


def _compare_coverage(
    plan: ValidationPlan,
    evidence: ValidationEvidence,
) -> tuple[set[str], dict[str, object]]:
    """Match planned observations by ID, condition, unit and allowance.

    References
    ----------
    .. [1] NASA-7009B, 4.1.1.5 and 4.2.5, acceptance criteria and validation domain.
    """
    targets = {target.sample_id: target for target in plan.targets}
    samples = {sample.identifier: sample for sample in evidence.samples}
    missing = sorted(targets.keys() - samples.keys())
    unexpected = sorted(samples.keys() - targets.keys())
    mismatches = {}
    matched = set()
    calibration_conditions_unreported = bool(
        evidence.calibration_sample_ids and not plan.calibration_condition_records
    )
    for sample_id in sorted(targets.keys() & samples.keys()):
        target, sample = targets[sample_id], samples[sample_id]
        reasons = []
        if sample.condition_record != target.condition_record:
            reasons.append("condition_record")
        if sample.measured.unit != target.allowable_error.unit:
            reasons.append("unit")
        if sample.allowable_error.value != target.allowable_error.value:
            reasons.append("allowable_error")
        if reasons:
            mismatches[sample_id] = reasons
        else:
            matched.add(sample_id)
    return matched, {
        "status": ("incomplete" if missing or unexpected or mismatches
                   or calibration_conditions_unreported else "complete"),
        "plan_review_record": plan.review_record,
        "required_sample_ids": sorted(targets),
        "required_targets": [
            {"sample_id": target.sample_id, "condition_record": target.condition_record,
             "unit": target.allowable_error.unit, "allowable_error": target.allowable_error.value,
             "reviewed_sample_sha256": target.reviewed_sample_sha256,
             **({"transient_target": asdict(target.transient_target)}
                if target.transient_target is not None else {}),
             **({"physical_criteria": asdict(target.physical_criteria)}
                if target.physical_criteria is not None else {})}
            for target in plan.targets
        ],
        "missing_sample_ids": missing,
        "unexpected_sample_ids": unexpected,
        "target_mismatches": mismatches,
        "calibration_conditions_unreported": calibration_conditions_unreported,
        "calibration_condition_records": list(plan.calibration_condition_records),
    }


def _check_sample_reviews(
    plan: ValidationPlan, evidence: ValidationEvidence, matched: set[str],
) -> tuple[dict[str, str], dict[str, object]]:
    """Compare the submitted rows with the digests from independent review.

    References
    ----------
    .. [1] NASA-7009B, 4.2.5, empirical referents and their validation records.
    """
    targets = {target.sample_id: target for target in plan.targets}
    excluded: dict[str, str] = {}
    reviewed = []
    for sample in evidence.samples:
        if sample.identifier not in matched:
            continue
        target = targets[sample.identifier]
        expected = target.reviewed_sample_sha256
        actual = validation_sample_sha256(
            evidence, sample, transient_target=target.transient_target,
            physical_criteria=target.physical_criteria,
        )
        if expected is None:
            excluded[sample.identifier] = "missing_reviewed_sample_sha256"
        elif actual != expected:
            excluded[sample.identifier] = "reviewed_sample_sha256_mismatch"
        reviewed.append({"sample_id": sample.identifier, "expected_sha256": expected,
                         "submitted_sha256": actual,
                         "status": excluded.get(sample.identifier, "matched")})
    status = "verified" if len(reviewed) == len(targets) and not excluded else "incomplete"
    return excluded, {"status": status, "samples": reviewed}


def _file_sha256(path: Path, label: str) -> str:
    """Hash a regular archive file incrementally.

    References
    ----------
    .. [1] NASA-7009B, 4.2.1, data pedigree and configuration records.
    .. [2] NIST FIPS 180-4, SHA-256, https://csrc.nist.gov/pubs/fips/180-4/upd1/final
    """
    digest = sha256()
    try:
        if not path.is_file():
            raise DomainError(f"{label} must be a regular archive file.")
        with path.open("rb") as stream:
            for block in iter(lambda: stream.read(1_048_576), b""):
                digest.update(block)
    except OSError as error:
        raise DomainError(f"{label} cannot be read.") from error
    return digest.hexdigest()


@dataclass(frozen=True, slots=True)
class _PressurePredictionWindow:
    """A pressure prediction's statistic and acquisition-aligned time support.

    References
    ----------
    .. [1] NASA-7009B, 4.2.5, validation referent and condition definition.
    """

    start_tick: int
    end_tick: int
    tick_period_s: float
    window_review_record: str
    metric: str


def _pressure_prediction_window(raw: object, unit: str) -> _PressurePredictionWindow | None:
    """Decode a pressure window with exact ticks and a matching statistic unit.

    References
    ----------
    .. [1] NASA-7009B, 4.2.1 and 4.2.5, documented representation and referent.
    """
    if raw is None:
        return None
    if not isinstance(raw, dict) or raw.keys() != {
        "start_tick", "end_tick", "tick_period_s", "window_review_record", "metric",
    }:
        raise DomainError("Pressure prediction window has an invalid schema.")
    start = integer(raw["start_tick"], "pressure prediction start tick")
    end = integer(raw["end_tick"], "pressure prediction end tick")
    period = real(raw["tick_period_s"], "pressure prediction tick period in seconds")
    review = text(raw["window_review_record"], "pressure prediction alignment review")
    metric = text(raw["metric"], "pressure prediction metric")
    if period <= 0 or end <= start:
        raise DomainError("Pressure prediction requires a positive tick period and nonempty window.")
    if (metric, unit) not in (("peak", "Pa"), ("signed_integral", "Pa*s")):
        raise DomainError("Pressure prediction metric and unit do not match.")
    return _PressurePredictionWindow(start, end, period, review, metric)


@dataclass(frozen=True, slots=True)
class _ShutdownPredictionScope:
    """Force, commands and time support of a predicted post-closure impulse.

    References
    ----------
    .. [1] NASA-7009B, 4.2.5 and 4.2.7, validation referents and timing uncertainty.
    """

    channel_id: str
    closing_commands: tuple[tuple[str, int], ...]
    window: ShutdownWindow
    final_closure_ticks: tuple[int, int]
    tick_period_s: float


def _shutdown_prediction_scope(raw: object, unit: str) -> _ShutdownPredictionScope | None:
    """Require an explicit closure-relative impulse and its complete time support.

    References
    ----------
    .. [1] NASA-7009B, 4.2.1 and 4.2.5, representation and referent definition.
    """
    if raw is None:
        return None
    if not isinstance(raw, dict) or raw.keys() != {
        "metric", "channel_id", "closing_commands", "window", "final_closure_ticks", "tick_period_s",
    }:
        raise DomainError("Shutdown prediction scope has an invalid schema.")
    if raw["metric"] != "post_closure_positive_impulse" or unit != "N*s":
        raise DomainError("Shutdown prediction requires post-closure positive impulse in N*s.")
    channel = text(raw["channel_id"], "shutdown prediction force channel")
    period = real(raw["tick_period_s"], "shutdown prediction tick period in seconds")
    window = raw["window"]
    if not isinstance(window, dict) or window.keys() != {
        "tail_end_tick", "force_to_command_error_ticks", "clock_alignment_record", "tail_window_review_record",
    }:
        raise DomainError("Shutdown prediction window has an invalid schema.")
    window = ShutdownWindow(**window)
    commands = raw["closing_commands"]
    if not isinstance(commands, list) or not commands:
        raise DomainError("Shutdown prediction requires a nonempty array of closing commands.")
    closing = {}
    for command in commands:
        if not isinstance(command, dict) or command.keys() != {"head_id", "tick"}:
            raise DomainError("Shutdown prediction close requires head_id and tick.")
        head = text(command["head_id"], "shutdown prediction closing valve")
        tick = integer(command["tick"], "shutdown prediction closing command tick")
        if head in closing or not 0 <= tick < window.tail_end_tick:
            raise DomainError("Shutdown prediction needs one closing command per valve before the tail ends.")
        closing[head] = tick
    bounds = raw["final_closure_ticks"]
    if not isinstance(bounds, list) or len(bounds) != 2:
        raise DomainError("Shutdown prediction requires two final closure endpoints.")
    low, high = (integer(tick, "shutdown prediction final closure endpoint") for tick in bounds)
    if period <= 0 or not 0 <= low <= high < window.tail_end_tick:
        raise DomainError("Shutdown prediction requires ordered closure endpoints before the tail end and positive tick period.")
    return _ShutdownPredictionScope(channel, tuple(sorted(closing.items())), window, (low, high), period)


@dataclass(frozen=True, slots=True)
class _FixedPredictionScope:
    """A predicted statistic tied to an acquisition channel and time base.

    References
    ----------
    .. [1] NASA-7009B, 4.2.5 and 4.2.6, validation metrics and referents.
    """

    channel_id: str
    start_tick: int
    end_tick: int
    tick_period_s: float
    window_review_record: str
    metric: str


def _fixed_prediction_scope(raw: object, unit: str) -> _FixedPredictionScope | None:
    """Read a channel and fixed window with exact ticks and explicit metric.

    References
    ----------
    .. [1] NASA-7009B, 4.2.1 and 4.2.5, documented input and metric scope.
    """
    if raw is None:
        return None
    if not isinstance(raw, dict) or raw.keys() != {
        "channel_id", "start_tick", "end_tick", "tick_period_s", "window_review_record", "metric",
    }:
        raise DomainError("Fixed-window prediction scope has an invalid schema.")
    channel = text(raw["channel_id"], "fixed prediction acquisition channel")
    start = integer(raw["start_tick"], "fixed prediction start tick")
    end = integer(raw["end_tick"], "fixed prediction end tick")
    period = real(raw["tick_period_s"], "fixed prediction tick period in seconds")
    review = text(raw["window_review_record"], "fixed prediction alignment review")
    metric = text(raw["metric"], "fixed prediction metric")
    if period <= 0 or end <= start:
        raise DomainError("Fixed prediction requires a positive tick period and nonempty window.")
    if (metric, unit) not in (
        ("peak", "N"), ("peak", "Pa"), ("signed_integral", "N*s"),
        ("signed_integral", "Pa*s"), ("positive_impulse", "N*s"),
    ):
        raise DomainError("Fixed prediction metric and unit do not match.")
    return _FixedPredictionScope(channel, start, end, period, review, metric)


_PredictionRows = tuple[
    dict[str, tuple[str, str, float, float | None]],
    dict[str, PressureReferent],
    dict[str, _PressurePredictionWindow],
    dict[str, _ShutdownPredictionScope],
    dict[str, _FixedPredictionScope],
]


def _parse_prediction_rows(
    payload: object, evidence: ValidationEvidence,
) -> _PredictionRows:
    """Read the versioned export and require unique, finite prediction rows.

    References
    ----------
    .. [1] NASA-7009B, 4.2.1, documented data units, sources and pedigree.
    .. [2] RFC 8259, JSON, https://www.rfc-editor.org/rfc/rfc8259
    """
    header = {
        "test_article_id": evidence.test_article_id,
        "model_configuration_id": evidence.model_configuration_id,
        "model_sha256": evidence.model_sha256,
    }
    if not isinstance(payload, dict) or payload.keys() != (header.keys() | {"schema", "predictions"}):
        raise DomainError("Prediction archive has an invalid document schema.")
    if any(payload[key] != value for key, value in header.items()):
        raise DomainError("Prediction archive header does not match the declared model and test article.")
    schema = payload["schema"]
    if schema not in ("rocket-propulsion-prediction-run-v1", "rocket-propulsion-prediction-run-v2",
                      "rocket-propulsion-prediction-run-v3", "rocket-propulsion-prediction-run-v4",
                      "rocket-propulsion-prediction-run-v5"):
        raise DomainError("Unknown prediction archive schema.")
    rows = payload["predictions"]
    if not isinstance(rows, list) or not rows:
        raise DomainError("Prediction archive requires a nonempty predictions array.")
    result: dict[str, tuple[str, str, float, float | None]] = {}
    pressure_referents: dict[str, PressureReferent] = {}
    keys = {"sample_id", "condition_record", "unit", "predicted", "prediction_bound"}
    pressure_windows: dict[str, _PressurePredictionWindow] = {}
    shutdown_scopes: dict[str, _ShutdownPredictionScope] = {}
    fixed_scopes: dict[str, _FixedPredictionScope] = {}
    if schema != "rocket-propulsion-prediction-run-v1":
        keys.add("pressure_referent")
    has_pressure_window = schema in ("rocket-propulsion-prediction-run-v3", "rocket-propulsion-prediction-run-v4")
    if has_pressure_window:
        keys.add("pressure_window")
    if schema in ("rocket-propulsion-prediction-run-v4", "rocket-propulsion-prediction-run-v5"):
        keys.add("shutdown_scope")
    if schema == "rocket-propulsion-prediction-run-v5":
        keys.add("fixed_window")
    for row in rows:
        if not isinstance(row, dict) or row.keys() != keys:
            raise DomainError("Prediction archive row has an invalid schema.")
        identifier = text(row["sample_id"], "archived prediction sample ID")
        if identifier in result:
            raise DomainError("Prediction archive contains a duplicate sample ID.")
        condition = text(row["condition_record"], "archived prediction condition")
        unit = text(row["unit"], "archived prediction unit")
        predicted = real(row["predicted"], "archived prediction")
        bound = None if row["prediction_bound"] is None else real(
            row["prediction_bound"], "archived prediction bound")
        if bound is not None and bound < 0.0:
            raise DomainError("Archived prediction bound must be nonnegative.")
        result[identifier] = (condition, unit, predicted, bound)
        pressure = _pressure_referent_from_payload(row.get("pressure_referent"))
        if pressure is not None:
            if unit not in ("Pa", "Pa*s"):
                raise DomainError("A pressure prediction referent requires Pa or Pa*s.")
            pressure_referents[identifier] = pressure
        if has_pressure_window:
            window = _pressure_prediction_window(row["pressure_window"], unit)
            if (window is None) != (pressure is None):
                raise DomainError("Pressure prediction referent and window must be declared together.")
            if window is not None:
                pressure_windows[identifier] = window
        shutdown = _shutdown_prediction_scope(row.get("shutdown_scope"), unit)
        fixed = _fixed_prediction_scope(row.get("fixed_window"), unit)
        if schema == "rocket-propulsion-prediction-run-v5" and pressure is not None and fixed is None:
            raise DomainError("A v5 pressure prediction referent requires a fixed window.")
        if shutdown is not None and fixed is not None:
            raise DomainError("A prediction cannot declare both fixed and closure-relative windows.")
        if shutdown is not None:
            shutdown_scopes[identifier] = shutdown
        if fixed is not None:
            fixed_scopes[identifier] = fixed
    return result, pressure_referents, pressure_windows, shutdown_scopes, fixed_scopes


def _check_model_run(
    evidence: ValidationEvidence,
    plan: ValidationPlan | None,
    model_run: ArchivedModelRun | None,
    matched: set[str],
    transient_run: ArchivedTransientRun | None,
    replayed: dict[str, TransientMeasurement],
) -> tuple[dict[str, str], dict[str, object]]:
    """Verify archived bytes and compare the reviewed predictions row by row.

    References
    ----------
    .. [1] NASA-7009B, 4.2.1 and 4.2.5, model pedigree and data comparisons.
    """
    approved_digest = None if plan is None else plan.reviewed_predictions_sha256
    expected_fixed = {
        target.sample_id: target.transient_target for target in plan.targets
        if isinstance(target.transient_target, TransientValidationTarget)
    } if plan is not None else {}
    expected_shutdown = {
        target.sample_id: target.transient_target for target in plan.targets
        if isinstance(target.transient_target, ShutdownValidationTarget)
    } if plan is not None else {}
    if model_run is None:
        status = "not_declared" if approved_digest is None else "missing_archives"
        excluded = {} if approved_digest is None else dict.fromkeys(matched, "missing_model_run_archives")
        if expected_fixed:
            excluded.update({
                key: ("missing_fixed_prediction_scope" if recipe.pressure_referent is None
                      else "missing_pressure_prediction_referent")
                for key, recipe in expected_fixed.items() if key in matched
            })
            status = "missing_archives"
        if expected_shutdown:
            excluded.update(dict.fromkeys(matched & expected_shutdown.keys(), "missing_shutdown_prediction_scope"))
            status = "missing_archives"
        return excluded, {"status": status}
    if not isinstance(model_run, ArchivedModelRun) or plan is None or approved_digest is None:
        raise DomainError("An archived model run requires an independently reviewed prediction archive digest.")
    model_digest = _file_sha256(model_run.model_path, "Archived model")
    if model_digest != evidence.model_sha256:
        raise DomainError("Archived model SHA-256 mismatch.")
    label = "Prediction archive"
    try:
        if not model_run.predictions_path.is_file():
            raise DomainError("Prediction archive must be a regular archive file.")
        with model_run.predictions_path.open("rb") as stream:
            data = stream.read(2_097_153)
    except OSError as error:
        raise DomainError("Prediction archive cannot be read.") from error
    if len(data) > 2_097_152:
        raise DomainError("Prediction archive exceeds the 2 MiB JSON limit.")
    output_digest = sha256(data).hexdigest()
    if output_digest != approved_digest:
        raise DomainError("Prediction archive SHA-256 mismatch.")
    try:
        source = data.decode("utf-8")
    except UnicodeError as error:
        raise DomainError("Prediction archive must use UTF-8.") from error
    archived, pressure_referents, pressure_windows, shutdown_scopes, fixed_scopes = cast(
        _PredictionRows,
        read_json_document(source, label, lambda payload: _parse_prediction_rows(payload, evidence)),
    )
    clocks = {(row.condition_record, row.trace.channel_id): row.trace.tick_period.value
              for row in transient_run.channels} if isinstance(transient_run, ArchivedTransientRun) else {}
    excluded = {}
    for sample in evidence.samples:
        if sample.identifier not in matched:
            continue
        row = archived.get(sample.identifier)
        bound = None if sample.prediction_bound is None else sample.prediction_bound.value
        if row is None:
            excluded[sample.identifier] = "missing_archived_prediction"
        elif row != (sample.condition_record, sample.predicted.unit, sample.predicted.value, bound):
            excluded[sample.identifier] = "archived_prediction_mismatch"
        else:
            recipe = expected_fixed.get(sample.identifier)
            reference = None if recipe is None else recipe.pressure_referent
            if pressure_referents.get(sample.identifier) != reference:
                excluded[sample.identifier] = ("pressure_prediction_referent_mismatch"
                                               if sample.identifier in pressure_referents
                                               else "missing_pressure_prediction_referent")
                continue
            fixed = fixed_scopes.get(sample.identifier)
            if reference is not None and recipe is not None and fixed is None:
                window = pressure_windows.get(sample.identifier)
                clock = clocks.get((sample.condition_record, recipe.channel_id))
                if window is None:
                    excluded[sample.identifier] = "missing_pressure_prediction_window"
                elif ((window.start_tick, window.end_tick, window.metric, window.window_review_record)
                      != (recipe.start_tick, recipe.end_tick, recipe.metric, recipe.window_review_record)):
                    excluded[sample.identifier] = "pressure_prediction_window_mismatch"
                elif clock is not None and window.tick_period_s != clock:
                    excluded[sample.identifier] = "pressure_prediction_clock_mismatch"
                if sample.identifier in excluded:
                    continue
            if recipe is not None or fixed is not None:
                clock = None if recipe is None else clocks.get((sample.condition_record, recipe.channel_id))
                reason = _fixed_prediction_mismatch(
                    recipe, fixed, clock, sample.identifier in replayed, sample.identifier in shutdown_scopes,
                )
                if reason is not None:
                    excluded[sample.identifier] = reason
            elif sample.identifier in expected_shutdown or sample.identifier in shutdown_scopes:
                shutdown_recipe = expected_shutdown.get(sample.identifier)
                clock = None if shutdown_recipe is None else clocks.get(
                    (sample.condition_record, shutdown_recipe.channel_id)
                )
                reason = _shutdown_prediction_mismatch(
                    shutdown_recipe, shutdown_scopes.get(sample.identifier), replayed.get(sample.identifier), clock,
                )
                if reason is not None:
                    excluded[sample.identifier] = reason
    required = {target.sample_id for target in plan.targets}
    missing = sorted(required - archived.keys())
    return excluded, {
        "status": "verified" if not missing and not excluded and matched == required else "incomplete",
        "model_sha256": model_digest,
        "prediction_archive_sha256": output_digest,
        "missing_planned_sample_ids": missing,
        "unplanned_archived_sample_ids": sorted(archived.keys() - required),
        "pressure_prediction_referent_exclusions": {
            key: reason for key, reason in excluded.items()
            if reason in ("pressure_prediction_referent_mismatch", "missing_pressure_prediction_referent")
        },
        "pressure_prediction_window_exclusions": {
            key: reason for key, reason in excluded.items()
            if reason in ("missing_pressure_prediction_window", "pressure_prediction_window_mismatch",
                          "pressure_prediction_clock_mismatch")
        },
        "shutdown_prediction_scope_exclusions": {
            key: reason for key, reason in excluded.items() if "shutdown_prediction_" in reason
        },
        "fixed_prediction_scope_exclusions": {
            key: reason for key, reason in excluded.items() if "fixed_prediction_" in reason
        },
        "mismatched_submitted_sample_ids": sorted(
            identifier for identifier, reason in excluded.items() if reason == "archived_prediction_mismatch"
        ),
    }


def _fixed_prediction_mismatch(
    recipe: TransientValidationTarget | None, scope: _FixedPredictionScope | None,
    clock: float | None, replay_verified: bool, shutdown_scope_present: bool,
) -> str | None:
    """Require a channel-aligned model result and successfully replayed measurement.

    References
    ----------
    .. [1] NASA-7009B, 4.2.5 and 4.2.6, matched metric, acquisition and time support.
    """
    if shutdown_scope_present or recipe is None:
        return "fixed_prediction_scope_mismatch"
    if scope is None:
        return "missing_fixed_prediction_scope"
    if (scope.channel_id, scope.start_tick, scope.end_tick, scope.metric, scope.window_review_record) != (
        recipe.channel_id, recipe.start_tick, recipe.end_tick, recipe.metric, recipe.window_review_record,
    ):
        return "fixed_prediction_scope_mismatch"
    if clock is None or not replay_verified:
        return "fixed_prediction_replay_unavailable"
    if scope.tick_period_s != clock:
        return "fixed_prediction_clock_mismatch"
    return None


def _shutdown_prediction_mismatch(
    recipe: ShutdownValidationTarget | None, scope: _ShutdownPredictionScope | None,
    measurement: TransientMeasurement | None, clock: float | None,
) -> str | None:
    """Match the prediction to the reviewed shutdown and replayed closure interval.

    References
    ----------
    .. [1] NASA-7009B, 4.2.5 and 4.2.7, matched experimental scope.
    """
    if scope is None:
        return "missing_shutdown_prediction_scope"
    if recipe is None:
        return "unexpected_shutdown_prediction_scope"
    commands = tuple(sorted((command.head_id, command.tick) for command in recipe.close_commands))
    if (scope.channel_id, scope.closing_commands, scope.window) != (recipe.channel_id, commands, recipe.window):
        return "shutdown_prediction_scope_mismatch"
    if clock is None or measurement is None:
        return "shutdown_prediction_replay_unavailable"
    if scope.tick_period_s != clock:
        return "shutdown_prediction_clock_mismatch"
    if scope.final_closure_ticks != measurement.final_closure_ticks:
        return "shutdown_prediction_closure_mismatch"
    return None


@dataclass(frozen=True, slots=True)
class _TransientReplay:
    """Verified measurements, independent valve reviews and replay diagnostics.

    References
    ----------
    .. [1] NASA-7009B, 4.2.5 and 4.3.8, evidence and reporting.
    """

    exclusions: dict[str, str]
    report: dict[str, object]
    measurements: dict[str, TransientMeasurement]
    valve_timing: dict[str, ValveResponseReview | str]


def _check_transient_replay(
    evidence: ValidationEvidence, plan: ValidationPlan | None,
    run: ArchivedTransientRun | None, matched: set[str],
) -> _TransientReplay:
    """Require planned measurements to match independent replay of the archived input.

    References
    ----------
    .. [1] NASA-7009B, 4.2.1, 4.2.5 and 4.2.7, referent pedigree and uncertainty.
    """
    targets = {target.sample_id: target for target in plan.targets
               if target.transient_target is not None} if plan is not None else {}
    if run is None:
        if not targets:
            return _TransientReplay({}, {"status": "not_declared"}, {}, {})
        return _TransientReplay(
            dict.fromkeys(matched & targets.keys(), "missing_transient_archive"),
            {"status": "missing_archive", "required_sample_ids": sorted(targets)}, {}, {},
        )
    if not targets or not isinstance(run, ArchivedTransientRun):
        raise DomainError("A transient archive requires a validation plan with replay targets.")
    if run.test_article_id != evidence.test_article_id:
        raise DomainError("Transient archive is for a different test article.")
    if sha256(transient_archive_bytes(run)).hexdigest() != evidence.raw_data_sha256:
        raise DomainError("Transient archive reconstruction SHA-256 mismatch.")
    channels = {(row.condition_record, row.trace.channel_id): row.trace for row in run.channels}
    pressure_channels = {(row.condition_record, row.trace.channel_id): row.pressure_referent
                         for row in run.channels if row.pressure_referent is not None}
    valve_acquisitions = {row.condition_record: row for row in run.valve_acquisitions}
    exclusions = {}
    rows = []
    measurements = {}
    valve_timing = {}
    for sample in evidence.samples:
        planned = targets.get(sample.identifier)
        if planned is None or sample.identifier not in matched:
            continue
        target = planned.transient_target
        assert target is not None
        trace = channels.get((sample.condition_record, target.channel_id))
        row: dict[str, object] = {"sample_id": sample.identifier}
        result: TransientMeasurement | str
        valve = valve_acquisitions.get(sample.condition_record)
        criteria = planned.physical_criteria
        if isinstance(target, ShutdownValidationTarget) and criteria is not None and criteria.require_valve_timing:
            valve_timing[sample.identifier] = (
                "missing_archived_valve_acquisition" if valve is None
                else shutdown_valve_timing(target, valve)
            )
        if trace is None:
            result = "missing_archived_transient_channel"
        elif isinstance(target, ShutdownValidationTarget):
            result = ("missing_archived_valve_acquisition" if valve is None
                      else shutdown_measurement(trace, target, valve))
        else:
            acquired_pressure = pressure_channels.get((sample.condition_record, target.channel_id))
            if acquired_pressure != target.pressure_referent:
                result = "pressure_measurement_referent_mismatch"
            else:
                result = transient_measurement(trace, target)
            if trace.unit == "Pa":
                row["pressure_referent_status"] = (
                    "mismatch" if acquired_pressure != target.pressure_referent else
                    "legacy_not_declared" if acquired_pressure is None else "matched"
                )
                if acquired_pressure is not None:
                    row["acquired_pressure_referent"] = asdict(acquired_pressure)
        if isinstance(result, str):
            exclusions[sample.identifier] = result
        else:
            row.update(unit=result.unit, value=result.value, absolute_error=result.absolute_error,
                       enclosure=result.enclosure)
            if result.final_closure_ticks is not None:
                row.update(final_closure_ticks=result.final_closure_ticks,
                           valve_timing_status=result.valve_timing_status)
            if result.unit != sample.measured.unit:
                exclusions[sample.identifier] = "transient_metric_unit_mismatch"
            elif (sample.measurement_bound is None or sample.measured.value != result.value
                  or sample.measurement_bound.value != result.absolute_error):
                exclusions[sample.identifier] = "replayed_transient_measurement_mismatch"
            else:
                measurements[sample.identifier] = result
        row["status"] = exclusions.get(sample.identifier, "matched")
        rows.append(row)
    return _TransientReplay(exclusions, {
        "status": "verified" if len(rows) == len(targets) and not exclusions else "incomplete",
        "required_sample_ids": sorted(targets),
        "samples": rows,
    }, measurements, valve_timing)


def _physical_row(
    criteria: PhysicalAcceptanceCriteria, measurement: TransientMeasurement | None,
    timing_review: ValveResponseReview | str | None, unavailable_reason: str,
) -> dict[str, object]:
    """Assess the outward enclosure against inclusive physical bounds.

    References
    ----------
    .. [1] NASA-7009B, 4.1.1.5 and 4.2.7, acceptance criteria and uncertainty.
    """
    band = "not_declared"
    if criteria.minimum is not None or criteria.maximum is not None:
        band = "unassessed"
        if measurement is not None:
            lo, hi = measurement.enclosure
            minimum = -float("inf") if criteria.minimum is None else criteria.minimum.value
            maximum = float("inf") if criteria.maximum is None else criteria.maximum.value
            band = ("fail" if hi < minimum or lo > maximum else
                    "pass" if minimum <= lo and hi <= maximum else "indeterminate")
    timing = "not_declared"
    if criteria.require_valve_timing:
        timing = timing_review.status if isinstance(timing_review, ValveResponseReview) else "unassessed"
    states = {band, timing} - {"not_declared"}
    status = ("fail" if "fail" in states else "pass" if states == {"pass"} else
              "indeterminate" if "indeterminate" in states else "unassessed")
    result: dict[str, object] = {
        "status": status, "measurement_band_status": band, "valve_timing_status": timing,
        "enclosure": None if measurement is None else measurement.enclosure,
    }
    if measurement is None:
        result["reason"] = unavailable_reason
    if isinstance(timing_review, ValveResponseReview):
        result["valve_timing_review"] = asdict(timing_review)
    elif criteria.require_valve_timing:
        result["valve_timing_reason"] = timing_review or unavailable_reason
    return result


def _check_physical_acceptance(
    plan: ValidationPlan | None, replay: _TransientReplay, measurement_exclusions: dict[str, str],
) -> dict[str, object]:
    """Report each declared physical criterion even if a model comparison is unavailable.

    References
    ----------
    .. [1] NASA-7009B, 4.1.1.5 and 4.3.8, approved requirements and reporting.
    """
    rows = []
    for target in (() if plan is None else plan.targets):
        criteria = target.physical_criteria
        if criteria is None:
            continue
        result = _physical_row(
            criteria, replay.measurements.get(target.sample_id), replay.valve_timing.get(target.sample_id),
            measurement_exclusions.get(target.sample_id, "missing_reviewed_measurement"),
        )
        rows.append({"sample_id": target.sample_id, "condition_record": target.condition_record,
                     "unit": target.allowable_error.unit, "criteria": asdict(criteria), **result})
    states = {row["status"] for row in rows}
    status = ("not_declared" if not rows else "fail" if "fail" in states else
              "pass" if states == {"pass"} else "unassessed")
    return {"status": status, "samples": rows}


def _check_acquisition_partition(
    evidence: ValidationEvidence, plan: ValidationPlan | None, matched: set[str],
) -> tuple[dict[str, str], dict[str, object]]:
    """Detect shared firings or native data across calibration and holdout.

    References
    ----------
    .. [1] NASA-7009B, 4.2.1 and 4.2.5, input pedigree and empirical datasets.
    """
    if plan is None or (plan.acquisition_partition is None and plan.reviewed_acquisition_partition_sha256 is None
                        and not plan.native_acquisitions
                        and all(target.transient_target is None for target in plan.targets)):
        return {}, {"status": "not_declared"}
    if plan.acquisition_partition is None or plan.reviewed_acquisition_partition_sha256 is None:
        return dict.fromkeys(matched, "missing_reviewed_acquisition_partition"), {"status": "missing_review"}
    digest = acquisition_partition_sha256(evidence, plan)
    if digest != plan.reviewed_acquisition_partition_sha256:
        return dict.fromkeys(matched, "acquisition_partition_sha256_mismatch"), {
            "status": "review_mismatch", "submitted_sha256": digest,
        }
    calibration = set(plan.calibration_condition_records)
    holdout = {target.condition_record for target in plan.targets}
    scopes = plan.acquisition_partition.slices
    declared = {scope.condition_record for scope in scopes}
    missing = sorted((calibration | holdout) - declared)
    unexpected = sorted(declared - (calibration | holdout))
    mapped_ids = {identifier for identifier, _ in plan.acquisition_partition.calibration_samples}
    mapped_conditions = {condition for _, condition in plan.acquisition_partition.calibration_samples}
    calibration_ids = set(evidence.calibration_sample_ids)
    missing_ids = sorted(calibration_ids - mapped_ids)
    unexpected_ids = sorted(mapped_ids - calibration_ids)
    mapping_matches = mapped_conditions == calibration
    report: dict[str, object] = {
        "reviewed_sha256": digest,
        "review_record": plan.acquisition_partition.review_record,
        "missing_conditions": missing,
        "unplanned_conditions": unexpected,
        "missing_calibration_sample_ids": missing_ids,
        "unplanned_calibration_sample_ids": unexpected_ids,
        "calibration_mapping_matches_declared_conditions": mapping_matches,
    }
    if missing or unexpected or missing_ids or unexpected_ids or not mapping_matches:
        return dict.fromkeys(matched, "incomplete_acquisition_partition"), {**report, "status": "incomplete"}
    calibration_by_file: dict[str, list[NativeAcquisitionSlice]] = {}
    calibration_by_firing: dict[str, set[str]] = {}
    for scope in scopes:
        if scope.condition_record in calibration:
            calibration_by_file.setdefault(scope.native_sha256, []).append(scope)
            calibration_by_firing.setdefault(scope.firing_id, set()).add(scope.condition_record)
    conflicts: dict[tuple[str, str], set[str]] = {}
    for scope in scopes:
        if scope.condition_record not in holdout:
            continue
        for condition in calibration_by_firing.get(scope.firing_id, ()):
            conflicts.setdefault((scope.condition_record, condition), set()).add("shared_firing")
        for other in calibration_by_file.get(scope.native_sha256, ()):
            if max(scope.first_index, other.first_index) < min(scope.past_end_index, other.past_end_index):
                conflicts.setdefault((scope.condition_record, other.condition_record), set()).add("overlapping_native_slice")
    conflicting_conditions = {condition for condition, _ in conflicts}
    native_review = _check_native_acquisitions(evidence, plan, {scope.native_sha256 for scope in scopes})
    report["native_acquisition_review"] = native_review
    report["status"] = ("overlap" if conflicts else "verified" if native_review["status"] == "verified"
                        else "native_archive_unverified")
    report["conflicts"] = [
        {"holdout_condition": condition, "calibration_condition": source, "reasons": sorted(reasons)}
        for (condition, source), reasons in sorted(conflicts.items())
    ]
    exclusions = {
        target.sample_id: "calibration_acquisition_overlap" for target in plan.targets
        if target.sample_id in matched and target.condition_record in conflicting_conditions
    }
    if native_review["status"] != "verified":
        for identifier in matched:
            exclusions.setdefault(identifier, "native_acquisition_archive_unverified")
    return exclusions, report


def _check_native_acquisitions(
    evidence: ValidationEvidence, plan: ValidationPlan, expected: set[str],
) -> dict[str, object]:
    """Verify native files against the hashes in the approved partition.

    References
    ----------
    .. [1] NASA-7009B, 4.2.1 and 4.2.5, input identity and pedigree.
    """
    sources = {source.native_sha256: source.path for source in plan.native_acquisitions}
    missing = sorted(expected - sources.keys())
    unplanned = sorted(sources.keys() - expected)
    repurposed = expected & {
        evidence.raw_data_sha256, evidence.model_sha256, plan.reviewed_predictions_sha256,
    }
    empty = expected & {sha256(b"").hexdigest()}
    verified = []
    unreadable = []
    mismatches = {}
    for digest in sorted((expected & sources.keys()) - repurposed - empty):
        try:
            actual = _file_sha256(sources[digest], "Native acquisition file")
        except DomainError:
            unreadable.append(digest)
            continue
        if actual != digest:
            mismatches[digest] = actual
        else:
            verified.append(digest)
    return {
        "status": ("verified" if not (missing or unplanned or repurposed or empty or unreadable or mismatches)
                   else "incomplete"),
        "verified_sha256": verified,
        "missing_sha256": missing,
        "unplanned_sha256": unplanned,
        "repurposed_sha256": sorted(repurposed),
        "empty_sha256": sorted(empty),
        "unreadable_sha256": unreadable,
        "sha256_mismatches": mismatches,
    }


def _check_native_trace_replay(
    plan: ValidationPlan | None, run: ArchivedTransientRun | None,
    partition_review: dict[str, object],
) -> tuple[dict[str, str], dict[str, object]]:
    """Verify approved CSV conversions without treating an alias as a source.

    References
    ----------
    .. [1] NASA-7009B, 4.2.5, empirical acquisition and processing.
    """
    partition = None if plan is None else plan.acquisition_partition
    if plan is None or partition is None or not partition.native_csv_traces:
        return {}, {"status": "not_declared"}
    targets = {
        target.sample_id: (target.condition_record, target.transient_target.channel_id)
        for target in plan.targets if target.transient_target is not None
    }
    mappings = {(row.condition_record, row.channel_id): row for row in partition.native_csv_traces}
    channels = {(row.condition_record, row.trace.channel_id): row.trace for row in (() if run is None else run.channels)}
    sources = {source.native_sha256: source.path for source in plan.native_acquisitions}
    source_review = partition_review.get("native_acquisition_review")
    verified_sources = set(source_review.get("verified_sha256", ())) if isinstance(source_review, dict) else set()
    slices: dict[tuple[str, str], list[tuple[int, int]]] = {}
    for item in partition.slices:
        slices.setdefault((item.condition_record, item.native_sha256), []).append((item.first_index, item.past_end_index))
    outcomes = {}
    for condition, channel in sorted(set(targets.values())):
        key = (condition, channel)
        mapping = mappings.get(key)
        trace = channels.get(key)
        if "reviewed_sha256" not in partition_review:
            outcome = "native_csv_mapping_unapproved"
        elif mapping is None:
            outcome = "native_csv_mapping_missing"
        elif trace is None:
            outcome = "native_csv_interpreted_channel_missing"
        elif not any(lo <= mapping.first_index and mapping.past_end_index <= hi
                     for lo, hi in slices.get((condition, mapping.native_sha256), ())):
            outcome = "native_csv_outside_reviewed_slice"
        elif mapping.native_sha256 not in sources:
            outcome = "native_csv_source_missing"
        elif mapping.native_sha256 not in verified_sources:
            outcome = "native_csv_source_unverified"
        else:
            outcome = native_csv_trace_status(sources[mapping.native_sha256], mapping, trace)
        outcomes[key] = outcome
    unexpected = sorted(mappings.keys() - set(targets.values()))
    exclusions = {identifier: outcomes[key] for identifier, key in targets.items() if outcomes[key] != "verified"}
    return exclusions, {
        "status": "incomplete" if exclusions or unexpected else "verified",
        "channels": [
            {"condition_record": condition, "channel_id": channel, "status": outcome}
            for (condition, channel), outcome in outcomes.items()
        ],
        "unplanned_channels": unexpected,
    }


def assess_validation(
    evidence: ValidationEvidence, *, plan: ValidationPlan | None = None,
    model_run: ArchivedModelRun | None = None,
    transient_run: ArchivedTransientRun | None = None,
) -> dict[str, object]:
    """Verify the archived raw bytes and summarize all holdout comparisons.

    Parameters
    ----------
    evidence : ValidationEvidence
        Supplied evidence record.
    plan : ValidationPlan or None
        Independently reviewed coverage, allowances and numeric row digests.
        When omitted, the comparison is scoped to submitted intervals.
    model_run : ArchivedModelRun or None
        Reviewed model and prediction archives. A comparison with a replay
        target requires an archive with matching metric and acquisition scope.
    transient_run : ArchivedTransientRun or None
        Acquisition export required when the plan declares transient targets.

    Returns
    -------
    dict
        ``comparison_status`` (``'fail'`` if any applicable sample fails,
        ``'pass'`` if all pass and any supplied plan has complete coverage and
        matched review digests, otherwise ``'unassessed'``), the per-sample
        rows, coverage, sample review and evidence records. A separate
        physical acceptance review compares replayed measurements with the
        plan's declared criteria and does not depend on model agreement.
        Model comparisons with replay targets also require a reviewed map
        of calibration and holdout conditions to native acquisition slices
        and the corresponding digest-verified original files.
        Hardware readiness and qualification are always reported as not established.

    Raises
    ------
    DomainError
        If the evidence type is wrong, the raw file cannot be read, or its SHA-256
        does not match, or a supplied plan is for a different test article.

    Notes
    -----
    The row digest binds submitted values to the reviewed transcription.
    With archived model evidence, the prediction export must identify the
    same model and match the submitted prediction rows. This verifies the
    supplied archive bytes and numeric records. Execution and provenance
    of the export require separate records [1]_.

    References
    ----------
    .. [1] NASA-7009B. NASA, "Standard for Models and Simulations, NASA-STD-7009B";
           4.3.8 report unresolved evidence and uncertainty.
    .. [2] NASA-TRL. NASA ESTO, "Technology Readiness Levels"; system
           hardware/environment demonstration cannot be inferred from this check.
    """
    if not isinstance(evidence, ValidationEvidence):
        raise DomainError("Expected complete validation evidence.")
    if plan is not None and (not isinstance(plan, ValidationPlan)
                             or plan.test_article_id != evidence.test_article_id):
        raise DomainError("A validation plan must match the evidence's test article.")
    raw_digest = _file_sha256(evidence.raw_data_path, "Raw experimental evidence")
    if raw_digest != evidence.raw_data_sha256:
        raise DomainError("Raw experimental evidence SHA-256 mismatch.")
    coverage: dict[str, object] = {"status": "not_declared"}
    sample_review: dict[str, object] = {"status": "not_declared"}
    exclusions: dict[str, str] = {}
    matched = {sample.identifier for sample in evidence.samples}
    if plan is not None:
        matched, coverage = _compare_coverage(plan, evidence)
        if coverage["calibration_conditions_unreported"]:
            matched.clear()
        exclusions, sample_review = _check_sample_reviews(plan, evidence, matched)
        matched.difference_update(exclusions)
    replay = _check_transient_replay(evidence, plan, transient_run, matched)
    partition_exclusions, partition_review = _check_acquisition_partition(evidence, plan, matched)
    native_exclusions, native_trace_review = _check_native_trace_replay(plan, transient_run, partition_review)
    archive_exclusions, archive_review = _check_model_run(
        evidence, plan, model_run, matched, transient_run, replay.measurements
    )
    exclusions.update(replay.exclusions)
    if native_exclusions:
        for identifier, reason in native_exclusions.items():
            exclusions.setdefault(identifier, reason)
        replay = _TransientReplay(replay.exclusions, replay.report, {
            key: value for key, value in replay.measurements.items() if key not in native_exclusions
        }, replay.valve_timing)
    physical_review = _check_physical_acceptance(plan, replay, exclusions)
    for identifier, reason in archive_exclusions.items():
        exclusions.setdefault(identifier, reason)
    for identifier, reason in partition_exclusions.items():
        exclusions.setdefault(identifier, reason)
    matched.difference_update(replay.exclusions)
    matched.difference_update(archive_exclusions)
    matched.difference_update(partition_exclusions)
    matched.difference_update(native_exclusions)
    rows = [
        compare_validation_sample(sample) if sample.identifier in matched else {
            "sample_id": sample.identifier, "condition_record": sample.condition_record,
            "unit": sample.measured.unit, "status": "unassessed",
            "reason": exclusions.get(sample.identifier, "outside_reviewed_plan"),
        }
        for sample in evidence.samples
    ]
    states = {row["status"] for row in rows}
    covered = plan is None or (coverage["status"] == "complete" and sample_review["status"] == "verified")
    covered = (covered and archive_review["status"] in ("not_declared", "verified")
               and replay.report["status"] in ("not_declared", "verified")
               and partition_review["status"] in ("not_declared", "verified")
               and native_trace_review["status"] in ("not_declared", "verified"))
    status = "fail" if "fail" in states else "pass" if covered and states == {"pass"} else "unassessed"
    return {"comparison_status": status, "samples": rows,
            "comparison_scope": ("supplied_numeric_intervals_only" if plan is None
                                 else "reviewed_holdout_targets"),
            "coverage": coverage,
            "sample_review": sample_review,
            "model_run_review": archive_review,
            "transient_replay_review": replay.report,
            "physical_acceptance_review": physical_review,
            "acquisition_partition_review": partition_review,
            "native_trace_review": native_trace_review,
            "empirical_validation_status": "not_established_by_this_numeric_check",
            "test_article_id": evidence.test_article_id,
            "model_configuration_id": evidence.model_configuration_id,
            "raw_data_sha256": raw_digest, "model_sha256": evidence.model_sha256,
            "model_configuration_verification": (
                "archived_model_bytes_match_declared_sha256" if "model_sha256" in archive_review
                else "supplied_identifier_and_digest_only_not_independently_verified"
            ),
            "independence_record": evidence.independence_record,
            "transcription_review_record": evidence.transcription_review_record,
            "uncertainty_basis": evidence.uncertainty_basis,
            "acceptance_approval_record": evidence.acceptance_approval_record,
            "calibration_sample_ids": list(evidence.calibration_sample_ids),
            "hardware_readiness": "not_assessed",
            "qualification_evidence": False,
            "evidence_authentication": "analyst_records_not_authenticated"}
