"""Conservative comparison of independent measurements and model predictions.

This module verifies supplied evidence and compares bounded intervals. It does
not authenticate the origin of files or certify a physical system. Acceptance
limits and uncertainty bounds must be established outside this implementation.

References
----------
NASA-7009B: 4.2.5-4.2.7 empirical referents, validation domain and uncertainty.
NASA-TRL: hardware demonstration is separate from model agreement.
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from hashlib import sha256
from math import isfinite, nextafter
from pathlib import Path
from sys import float_info

from rocket_propulsion.core.errors import DomainError
from rocket_propulsion.core.validation import records, strings, text

from .evidence import CitedScalar, quantity


@dataclass(frozen=True, slots=True)
class ValidationSample:
    """One holdout datum and prediction under a documented common condition.

    Bounds are absolute error bounds in the measurand's units, not standard
    deviations or confidence intervals. None means unknown, never exact zero.
    Domain membership requires external review of the condition record; the
    boolean records that review and is not inferred from one numeric range.

    References
    ----------
    NASA-7009B: 4.2.5 validation referents/domain; 4.2.7 uncertainty rationale.
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
        NASA-7009B: 4.2.1 units; 4.2.5/4.2.7 validation and uncertainty records.
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

    The raw file digest binds the reviewed snapshot, not the analyst's numeric
    transcription to it. A transcription-review record remains mandatory.
    Calibration sample IDs cannot overlap holdout IDs. Dataset authorship,
    applicability and sensor calibration are human-reviewed records.

    References
    ----------
    NASA-7009B: 4.2.1 data pedigree, 4.2.5 empirical data and 4.2.7 uncertainty.
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
        NASA-7009B: 4.1.1.5 project acceptance criteria; 4.2.5 data referents.
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


def _outward_endpoint(value: Fraction, *, upward: bool) -> float:
    """Publish a finite enclosure endpoint without rounding into the interval.

    References
    ----------
    NASA-7009B: declared numerical representation limits and error reporting.
    BIPM-SI: preserve the units of bounded errors; rounding method derived here.
    """
    try:
        binary = float(value)
    except OverflowError as error:
        raise DomainError("Validation interval exceeds the finite binary domain.") from error
    if not isfinite(binary):
        raise DomainError("Validation interval exceeds the finite binary domain.")
    if upward and Fraction(binary) < value:
        binary = nextafter(binary, float("inf"))
    elif not upward and Fraction(binary) > value:
        binary = nextafter(binary, -float("inf"))
    if not isfinite(binary):
        raise DomainError("Validation interval has no finite outward binary bound.")
    return binary


def compare_validation_sample(sample: ValidationSample) -> dict[str, object]:
    """Bound |prediction-truth| by the supplied bounded error intervals.

    max error = |prediction-measurement| + measurement_bound + prediction_bound.
    min error = max(0, |prediction-measurement| - both bounds). Pass requires
    the complete residual interval inside the approved allowance. A partial
    overlap is indeterminate; no Gaussian assumptions or fitted tolerances.
    This conservative interval rule is our stated method, not a NASA formula.
    Exact rational intervals are rounded outward for publication. Decisions
    use those reported bounds and an outward-rounded numerical review guard.

    References
    ----------
    NASA-7009B: 4.2.5 and 4.2.7 require justified metrics and uncertainty.
    BIPM-SI: like-unit interval arithmetic; triangle inequality derived here.
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
    high = _outward_endpoint(discrepancy + bound, upward=True)
    low = _outward_endpoint(max(Fraction(0), discrepancy - bound), upward=False)
    allowance = Fraction(sample.allowable_error.value)
    # Expose a near-boundary result for review rather than hiding roundoff in
    # a pass. Publish an upper bound for 16 epsilons of the reported scale.
    guard_value = 16 * Fraction(float_info.epsilon) * max(Fraction(high), allowance)
    guard = _outward_endpoint(guard_value, upward=True)
    status = ("pass" if Fraction(high) + Fraction(guard) < allowance else
              "fail" if Fraction(low) > allowance + Fraction(guard) else "indeterminate")
    return {**result, "status": status, "discrepancy": float(discrepancy),
            "minimum_error": low, "maximum_error": high,
            "allowable_error": sample.allowable_error.value, "numerical_guard": guard,
            "numerical_method": "exact_intervals_with_outward_binary_endpoints"}


def assess_validation(evidence: ValidationEvidence) -> dict[str, object]:
    """Verify the archived bytes and summarize all holdout samples without TRL promotion.

    A supplied digest of model code is recorded, not independently verified by
    this function. Hashes identify content, not authorship or test authenticity.

    References
    ----------
    NASA-7009B: 4.3.8 report unresolved evidence and uncertainty.
    NASA-TRL: system hardware/environment demonstration cannot be inferred here.
    """
    if not isinstance(evidence, ValidationEvidence):
        raise DomainError("Expected complete validation evidence.")
    digest = sha256()
    try:
        with evidence.raw_data_path.open("rb") as stream:
            for block in iter(lambda: stream.read(1_048_576), b""):
                digest.update(block)
    except OSError as error:
        raise DomainError("Raw experimental evidence cannot be read.") from error
    if digest.hexdigest() != evidence.raw_data_sha256:
        raise DomainError("Raw experimental evidence SHA-256 mismatch.")
    rows = [compare_validation_sample(sample) for sample in evidence.samples]
    states = {row["status"] for row in rows}
    status = "fail" if "fail" in states else "pass" if states == {"pass"} else "unassessed"
    return {"comparison_status": status, "samples": rows,
            "comparison_scope": "supplied_numeric_intervals_only",
            "empirical_validation_status": "not_established_by_this_numeric_check",
            "test_article_id": evidence.test_article_id,
            "model_configuration_id": evidence.model_configuration_id,
            "raw_data_sha256": digest.hexdigest(), "model_sha256": evidence.model_sha256,
            "model_configuration_verification": "supplied_identifier_and_digest_only_not_independently_verified",
            "independence_record": evidence.independence_record,
            "transcription_review_record": evidence.transcription_review_record,
            "uncertainty_basis": evidence.uncertainty_basis,
            "acceptance_approval_record": evidence.acceptance_approval_record,
            "calibration_sample_ids": list(evidence.calibration_sample_ids),
            "hardware_readiness": "not_assessed",
            "qualification_evidence": False,
            "evidence_authentication": "analyst_records_not_authenticated"}
