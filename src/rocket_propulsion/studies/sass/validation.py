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

    The raw-file digest binds the reviewed data snapshot, not the analyst's numeric
    transcription of it, so a transcription-review record is mandatory. Calibration
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


def _outward_endpoint(value: Fraction, *, upward: bool) -> float:
    """Publish a finite enclosure endpoint without rounding into the interval.

    References
    ----------
    .. [1] NASA-7009B. NASA, "Standard for Models and Simulations, NASA-STD-7009B";
           declared numerical representation limits and error reporting.
    .. [2] BIPM-SI. BIPM, "SI Brochure, 9th edition"; preserve the units of bounded
           errors; rounding method defined in this module.
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
    """Verify the archived raw bytes and summarize all holdout comparisons.

    Parameters
    ----------
    evidence : ValidationEvidence
        Complete evidence record.

    Returns
    -------
    dict
        ``comparison_status`` (``'fail'`` if any sample fails, ``'pass'`` only if all
        pass, otherwise ``'unassessed'``), the per-sample rows and the evidence
        records. Hardware readiness and qualification are always reported as not
        established.

    Raises
    ------
    DomainError
        If the evidence type is wrong, the raw file cannot be read, or its SHA-256
        does not match.

    Notes
    -----
    The supplied model digest is recorded, not independently verified. Hashes
    identify content, not authorship or test authenticity [1]_; agreement does not
    promote hardware readiness [2]_.

    References
    ----------
    .. [1] NASA-7009B. NASA, "Standard for Models and Simulations, NASA-STD-7009B";
           4.3.8 report unresolved evidence and uncertainty.
    .. [2] NASA-TRL. NASA ESTO, "Technology Readiness Levels"; system
           hardware/environment demonstration cannot be inferred from this check.
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
