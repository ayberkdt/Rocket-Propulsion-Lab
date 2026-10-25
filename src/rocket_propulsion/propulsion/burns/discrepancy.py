"""Model-form and qualification-to-flight propulsion discrepancy evidence.

Hot-fire repeatability does not cover bias between a model and reality, nor the
change between qualification and flight application domains [1]_ [2]_. This
module keeps those epistemic effects separate from aleatory engine/run
variation, calibrates them from independent validation-ratio evidence [3]_ with
the Mandel-Paule random-effects estimator [4]_ [5]_, and composes them with the
hot-fire calibration so that no physical scale is corrected twice.

A discrepancy calibration is valid only inside its declared applicability
domain; it is a constant factor, not an extrapolating response surface [6]_.

References
----------
.. [1] NASA-STD-7009B, "Standard for Models and Simulations".
       https://standards.nasa.gov/standard/NASA/NASA-STD-7009
.. [2] NASA-HDBK-7009, "NASA Handbook for Models and Simulations".
       https://ntrs.nasa.gov/citations/20140002378
.. [3] NASA, "Calibration of Predictor Models Using Multiple Validation
       Experiments". https://ntrs.nasa.gov/citations/20150006032
.. [4] PAULE-MANDEL-1982. R. C. Paule and J. Mandel, "Consensus Values and
       Weighting Factors", J. Res. Natl. Bur. Stand. 87(5), 377-385, 1982.
.. [5] NIST Dataplot, "Consensus Mean" (heteroscedastic random effects and
       Mandel-Paule estimation).
       https://www.itl.nist.gov/div898/software/dataplot/refman1/auxillar/consmean.htm
.. [6] NASA/SP-2009-569, "Bayesian Inference for NASA Probabilistic Risk and
       Reliability Analysis". https://ntrs.nasa.gov/citations/20090023159
"""

from __future__ import annotations

import hashlib
import json
from collections import Counter
from dataclasses import dataclass
from enum import Enum
from math import isfinite, sqrt
from typing import Any

from rocket_propulsion.core.errors import DomainError, InputError
from rocket_propulsion.core.validation import read_json_document

from .calibration import CalibrationScope, MandelPauleEvidence, _mandel_paule
from .ensemble import ScenarioConditionalUncertainty
from .multivariate_calibration import HotFireMultivariateCalibration
from .serialization import canonical_json_bytes
from .tabulated import TabulatedBurnArtifact, scale_tabulated_artifact
from .uncertainty import (
    CorrelationModel,
    DistributionKind,
    UncertainParameter,
    UncertaintyClass,
)

DISCREPANCY_CALIBRATION_SCHEMA = "rocket_propulsion_discrepancy_calibration_v1"
DISCREPANCY_REFERENCE_IDS = (
    "NASA-STD-7009B",
    "NASA-HDBK-7009",
    "NASA-20150006032",
    "NASA-SP-2009-569",
    "NIST-DATAPLOT-CONSENSUS-MEAN",
)
_SCALE_COMPONENT_SEPARATOR = "::"
_ALLOWED_TARGET_SCALES = {
    "thrust_scale",
    "specific_impulse_scale",
    "time_scale",
    "axial_efficiency_scale",
}


class DiscrepancyKind(str, Enum):
    """Epistemic physical meaning assigned to validation-ratio evidence [1]_ [2]_.

    The kinds are kept distinct so that one data set cannot silently serve both
    roles.

    Attributes
    ----------
    MODEL_FORM
        Measured/model residual structure remaining after parameter
        calibration.
    QUALIFICATION_TO_FLIGHT
        Ratio of flight-domain to qualification-domain performance.

    References
    ----------
    .. [1] NASA-STD-7009B, "Standard for Models and Simulations".
           https://standards.nasa.gov/standard/NASA/NASA-STD-7009
    .. [2] NASA-HDBK-7009, "NASA Handbook for Models and Simulations".
           https://ntrs.nasa.gov/citations/20140002378
    """

    MODEL_FORM = "model-form"
    QUALIFICATION_TO_FLIGHT = "qualification-to-flight"


@dataclass(frozen=True, slots=True)
class DiscrepancyRatioObservation:
    """One positive target/reference ratio with its standard uncertainty.

    Parameters
    ----------
    case_id : str
        Nonempty validation-case identifier.
    applicability_id : str
        Nonempty applicability domain; all observations pooled in one
        calibration must share it.
    ratio : float
        Positive finite ratio. For model-form evidence it is observed divided
        by predicted performance [1]_; for qualification transfer it is
        flight-representative divided by qualification-domain performance.
    standard_uncertainty : float
        Finite, nonnegative standard uncertainty of the referent ratio [2]_.
    source : str
        Nonempty provenance.

    Raises
    ------
    DomainError
        If any constraint above is violated.

    References
    ----------
    .. [1] NASA, "Calibration of Predictor Models Using Multiple Validation
           Experiments". https://ntrs.nasa.gov/citations/20150006032
    .. [2] NIST TN 1297, "Guidelines for Evaluating and Expressing the Uncertainty
           of NIST Measurement Results", 1994.
           https://www.nist.gov/pml/nist-technical-note-1297
    """

    case_id: str
    applicability_id: str
    ratio: float
    standard_uncertainty: float
    source: str

    def __post_init__(self) -> None:
        if not self.case_id.strip() or not self.applicability_id.strip():
            raise DomainError("Discrepancy case and applicability identifiers are required.")
        if not isfinite(self.ratio) or self.ratio <= 0.0:
            raise DomainError("Discrepancy ratio must be positive and finite.")
        if not isfinite(self.standard_uncertainty) or self.standard_uncertainty < 0.0:
            raise DomainError(
                "Discrepancy measurement standard uncertainty must be finite and nonnegative."
            )
        if not self.source.strip():
            raise DomainError("Discrepancy observation source is required.")


def qualified_scale_component_name(
    target_scale: str,
    *,
    uncertainty_class: UncertaintyClass,
    component_id: str,
) -> str:
    """Return the ensemble name of one factor acting on a physical scale.

    The name ``"<target>::<class>:<component>"`` lets several factors multiply
    the same scale while each keeps its aleatory or epistemic classification
    [1]_ [2]_.

    Parameters
    ----------
    target_scale : str
        ``thrust_scale``, ``specific_impulse_scale``, ``time_scale`` or
        ``axial_efficiency_scale``.
    uncertainty_class : UncertaintyClass
        Classification encoded in the name.
    component_id : str
        Nonempty component identifier without ``"::"``.

    Returns
    -------
    str
        Qualified component name.

    Raises
    ------
    DomainError
        If the target, component identifier or class is invalid.

    References
    ----------
    .. [1] NASA-STD-7009B, "Standard for Models and Simulations".
           https://standards.nasa.gov/standard/NASA/NASA-STD-7009
    .. [2] NASA/SP-2009-569, "Bayesian Inference for NASA Probabilistic Risk and
           Reliability Analysis". https://ntrs.nasa.gov/citations/20090023159
    """

    normalized_target = target_scale.strip()
    normalized_component = component_id.strip()
    if normalized_target not in _ALLOWED_TARGET_SCALES:
        raise DomainError("Unsupported discrepancy target scale.")
    if not normalized_component or _SCALE_COMPONENT_SEPARATOR in normalized_component:
        raise DomainError("Discrepancy component ID is empty or contains a reserved separator.")
    if not isinstance(uncertainty_class, UncertaintyClass):
        raise DomainError("A valid uncertainty class is required for a scale component.")
    return (
        f"{normalized_target}{_SCALE_COMPONENT_SEPARATOR}"
        f"{uncertainty_class.value}:{normalized_component}"
    )


@dataclass(frozen=True, slots=True)
class ScaleDiscrepancyCalibration:
    """Random-effects estimate of one epistemic propulsion scale discrepancy.

    The consensus ratio is a deterministic correction. The predictive epistemic
    variance combines the between-case discrepancy variance and the squared
    standard uncertainty of the consensus correction. Referent measurement
    uncertainty informs the fit but is not relabelled as physical engine
    variability [1]_ [2]_. Instances are produced by
    :func:`calibrate_scale_discrepancy`.

    Parameters
    ----------
    target_scale, component_id : str
        Corrected scale and component identifier.
    kind : DiscrepancyKind
        Physical meaning of the evidence.
    applicability_id : str
        Common applicability domain of the observations.
    evidence_id : str
        Identifier of the evidence set, used to detect double counting.
    method, source, rationale : str
        Method identifier, provenance and justification.
    observations : tuple of DiscrepancyRatioObservation
        Input ratios.
    consensus_scale : float
        Mandel-Paule weighted mean ratio ``mu`` [3]_ [4]_.
    consensus_standard_uncertainty : float
        ``(sum_i w_i)**-0.5`` with ``w_i = 1 / (tau**2 + u_i**2)``.
    between_case_variance : float
        Mandel-Paule between-case variance ``tau**2``.
    predictive_epistemic_variance : float
        ``tau**2 + (sum_i w_i)**-1``.
    evidence : MandelPauleEvidence
        Solver evidence.
    relative_tolerance : float
        Mandel-Paule tolerance.
    maximum_iterations : int
        Mandel-Paule iteration limit.
    reference_ids, warnings : tuple of str
        Source identifiers and interpretation warnings.

    References
    ----------
    .. [1] NASA, "Calibration of Predictor Models Using Multiple Validation
           Experiments". https://ntrs.nasa.gov/citations/20150006032
    .. [2] NASA-HDBK-7009, "NASA Handbook for Models and Simulations".
           https://ntrs.nasa.gov/citations/20140002378
    .. [3] PAULE-MANDEL-1982. R. C. Paule and J. Mandel, "Consensus Values and
           Weighting Factors", J. Res. Natl. Bur. Stand. 87(5), 377-385, 1982.
    .. [4] NIST Dataplot, "Consensus Mean" (heteroscedastic random effects and
           Mandel-Paule estimation).
           https://www.itl.nist.gov/div898/software/dataplot/refman1/auxillar/consmean.htm
    """

    target_scale: str
    component_id: str
    kind: DiscrepancyKind
    applicability_id: str
    evidence_id: str
    method: str
    source: str
    rationale: str
    observations: tuple[DiscrepancyRatioObservation, ...]
    consensus_scale: float
    consensus_standard_uncertainty: float
    between_case_variance: float
    predictive_epistemic_variance: float
    evidence: MandelPauleEvidence
    relative_tolerance: float
    maximum_iterations: int
    reference_ids: tuple[str, ...]
    warnings: tuple[str, ...]

    @property
    def calibration_hash(self) -> str:
        """SHA-256 hex digest of the canonical JSON of the complete record.

        The digest covers evidence, assumptions, results and sources; it is a
        reproducible content identity, not an authentication of the data.
        """

        return hashlib.sha256(canonical_json_bytes(self)).hexdigest()

    def residual_uncertain_parameter(
        self,
        *,
        coverage_factor: float = 3.0,
    ) -> UncertainParameter:
        """Return the bounded epistemic residual factor after consensus correction.

        Apply ``consensus_scale`` to the deterministic artifact first. The
        returned truncated-normal factor is centered at one with relative
        standard deviation ``r = sqrt(predictive_epistemic_variance) / mu`` and
        bounds ``[max(1e-12, 1 - k r), 1 + k r]`` for coverage factor ``k``
        [1]_. It is classified as epistemic [2]_ [3]_.

        Parameters
        ----------
        coverage_factor : float, optional
            Positive finite truncation half-width in standard deviations.

        Returns
        -------
        UncertainParameter
            Epistemic factor with a qualified component name.

        Raises
        ------
        DomainError
            If ``coverage_factor`` is invalid or the predictive variance is zero.

        References
        ----------
        .. [1] NIST TN 1297, "Guidelines for Evaluating and Expressing the Uncertainty
               of NIST Measurement Results", 1994.
               https://www.nist.gov/pml/nist-technical-note-1297
        .. [2] NASA-HDBK-7009, "NASA Handbook for Models and Simulations".
               https://ntrs.nasa.gov/citations/20140002378
        .. [3] NASA/SP-2009-569, "Bayesian Inference for NASA Probabilistic Risk and
               Reliability Analysis". https://ntrs.nasa.gov/citations/20090023159
        """

        if not isfinite(coverage_factor) or coverage_factor <= 0.0:
            raise DomainError("Discrepancy coverage factor must be positive and finite.")
        relative_sigma = sqrt(self.predictive_epistemic_variance) / self.consensus_scale
        if relative_sigma <= 0.0:
            raise DomainError("Discrepancy predictive variance is zero.")
        return UncertainParameter(
            name=qualified_scale_component_name(
                self.target_scale,
                uncertainty_class=UncertaintyClass.EPISTEMIC,
                component_id=f"{self.kind.value}:{self.component_id}",
            ),
            nominal=1.0,
            lower=max(1e-12, 1.0 - coverage_factor * relative_sigma),
            upper=1.0 + coverage_factor * relative_sigma,
            distribution=DistributionKind.TRUNCATED_NORMAL,
            uncertainty_class=UncertaintyClass.EPISTEMIC,
            source=f"{self.source}; discrepancy sha256={self.calibration_hash}",
            rationale=(
                f"Predictive {self.kind.value} residual for {self.applicability_id} after "
                f"applying consensus correction {self.consensus_scale:.17g}; {self.rationale}"
            ),
            standard_deviation=relative_sigma,
        )

    def to_dict(self, *, include_hash: bool = True) -> dict[str, Any]:
        """Return the complete versioned discrepancy-calibration document.

        Parameters
        ----------
        include_hash : bool, optional
            Add ``calibration_hash`` to the document.

        Returns
        -------
        dict
            JSON-compatible document with ``schema`` and every record field.
        """

        payload = json.loads(canonical_json_bytes(self).decode("utf-8"))
        payload = {"schema": DISCREPANCY_CALIBRATION_SCHEMA, **payload}
        if include_hash:
            payload["calibration_hash"] = self.calibration_hash
        return payload

    def to_json(self, *, indent: int | None = 2) -> str:
        """Serialize the discrepancy calibration as JSON with sorted keys.

        Non-finite numbers are rejected, so the output is standard JSON [1]_
        [2]_.

        Parameters
        ----------
        indent : int or None, optional
            Indentation passed to ``json.dumps``.

        Returns
        -------
        str
            JSON text of :meth:`to_dict`.

        References
        ----------
        .. [1] IETF RFC 8259, "The JavaScript Object Notation (JSON) Data Interchange
               Format". https://www.rfc-editor.org/rfc/rfc8259
        .. [2] PYTHON-3.12-JSON. Python 3.12, ``json``: encoder/decoder hooks,
               duplicate-member and nonstandard numeric token behavior.
        """

        return json.dumps(
            self.to_dict(),
            ensure_ascii=False,
            allow_nan=False,
            sort_keys=True,
            indent=indent,
        )

    @classmethod
    def from_dict(cls, payload: dict[str, Any]) -> ScaleDiscrepancyCalibration:
        """Rebuild a discrepancy calibration and verify that it reproduces.

        The calibration is recomputed from the stored observations and
        controls; the stored hash and every stored field must match.

        Parameters
        ----------
        payload : dict
            Document produced by :meth:`to_dict`.

        Returns
        -------
        ScaleDiscrepancyCalibration
            The recomputed calibration.

        Raises
        ------
        InputError
            If the document is malformed, has another schema, or its hash or
            content does not reproduce.
        """

        if not isinstance(payload, dict):
            raise InputError("Discrepancy calibration root must be an object.")
        if payload.get("schema") != DISCREPANCY_CALIBRATION_SCHEMA:
            raise InputError("Unsupported discrepancy calibration schema.")
        try:
            observations = tuple(
                DiscrepancyRatioObservation(
                    case_id=str(item["case_id"]),
                    applicability_id=str(item["applicability_id"]),
                    ratio=float(item["ratio"]),
                    standard_uncertainty=float(item["standard_uncertainty"]),
                    source=str(item["source"]),
                )
                for item in payload["observations"]
            )
            rebuilt = calibrate_scale_discrepancy(
                str(payload["target_scale"]),
                str(payload["component_id"]),
                DiscrepancyKind(payload["kind"]),
                observations,
                evidence_id=str(payload["evidence_id"]),
                source=str(payload["source"]),
                rationale=str(payload["rationale"]),
                relative_tolerance=float(payload["relative_tolerance"]),
                maximum_iterations=int(payload["maximum_iterations"]),
            )
        except (KeyError, TypeError, ValueError, DomainError) as error:
            raise InputError("Malformed discrepancy calibration.") from error
        if payload.get("calibration_hash") != rebuilt.calibration_hash:
            raise InputError("Discrepancy calibration hash does not reproduce.")
        if canonical_json_bytes(payload) != canonical_json_bytes(rebuilt.to_dict()):
            raise InputError("Discrepancy calibration semantics do not reproduce.")
        return rebuilt

    @classmethod
    def from_json(cls, text: str) -> ScaleDiscrepancyCalibration:
        """Parse strict JSON text [1]_ [2]_ and verify it with :meth:`from_dict`.

        Raises
        ------
        InputError
            If the text is not valid JSON or the document does not reproduce.

        References
        ----------
        .. [1] IETF RFC 8259, "The JavaScript Object Notation (JSON) Data Interchange
               Format". https://www.rfc-editor.org/rfc/rfc8259
        .. [2] PYTHON-3.12-JSON. Python 3.12, ``json``: encoder/decoder hooks,
               duplicate-member and nonstandard numeric token behavior.
        """

        return read_json_document(text, "Discrepancy calibration", cls.from_dict)


def calibrate_scale_discrepancy(
    target_scale: str,
    component_id: str,
    kind: DiscrepancyKind,
    observations: tuple[DiscrepancyRatioObservation, ...],
    *,
    evidence_id: str,
    source: str,
    rationale: str,
    relative_tolerance: float = 1e-10,
    maximum_iterations: int = 200,
) -> ScaleDiscrepancyCalibration:
    """Calibrate a consensus correction and predictive epistemic discrepancy.

    Parameters
    ----------
    target_scale : str
        ``thrust_scale``, ``specific_impulse_scale``, ``time_scale`` or
        ``axial_efficiency_scale``.
    component_id : str
        Nonempty component identifier without ``"::"``.
    kind : DiscrepancyKind
        Physical meaning of the evidence.
    observations : tuple of DiscrepancyRatioObservation
        At least three independent cases with unique ``case_id`` and one
        common ``applicability_id``.
    evidence_id, source, rationale : str
        Nonempty evidence identifier, provenance and justification.
    relative_tolerance : float, optional
        Positive tolerance of the Mandel-Paule solve.
    maximum_iterations : int, optional
        Iteration limit of the Mandel-Paule solve, at least ten.

    Returns
    -------
    ScaleDiscrepancyCalibration
        Consensus correction, between-case and predictive variances, solver
        evidence and warnings.

    Raises
    ------
    DomainError
        If the inputs violate the constraints above or the variance solve does
        not converge.

    Notes
    -----
    Independent validation ratios ``x_i`` with referent uncertainties ``u_i``
    [1]_ are combined with the heteroscedastic Mandel-Paule model [2]_ [3]_:
    ``tau**2`` solves ``sum_i w_i (x_i - mu)**2 = k - 1`` with
    ``w_i = 1 / (tau**2 + u_i**2)`` and weighted mean
    ``mu = sum_i w_i x_i / sum_i w_i`` [4]_. The predictive variance of the
    discrepancy in a future application within the same domain is
    ``tau**2 + (sum_i w_i)**-1``, the between-case variance plus the squared
    inverse-variance standard uncertainty of the consensus. This is a constant
    discrepancy within one declared applicability domain, not an extrapolating
    response surface [5]_. With fewer than six cases the random-effects
    uncertainty can be too narrow [3]_.

    References
    ----------
    .. [1] NASA, "Calibration of Predictor Models Using Multiple Validation
           Experiments". https://ntrs.nasa.gov/citations/20150006032
    .. [2] PAULE-MANDEL-1982. R. C. Paule and J. Mandel, "Consensus Values and
           Weighting Factors", J. Res. Natl. Bur. Stand. 87(5), 377-385, 1982.
    .. [3] NIST Dataplot, "Consensus Mean" (heteroscedastic random effects and
           Mandel-Paule estimation).
           https://www.itl.nist.gov/div898/software/dataplot/refman1/auxillar/consmean.htm
    .. [4] NIST-WEIGHTED-MEAN. NIST Dataplot, "Weighted Mean", equation 2-20,
           normalized weighted mean.
    .. [5] NASA-HDBK-7009, "NASA Handbook for Models and Simulations".
           https://ntrs.nasa.gov/citations/20140002378
    """

    observations = tuple(observations)
    normalized_target = target_scale.strip()
    normalized_component = component_id.strip()
    if normalized_target not in _ALLOWED_TARGET_SCALES:
        raise DomainError("Unsupported discrepancy target scale.")
    if not normalized_component or _SCALE_COMPONENT_SEPARATOR in normalized_component:
        raise DomainError("Discrepancy component ID is empty or contains a reserved separator.")
    if not isinstance(kind, DiscrepancyKind):
        raise DomainError("Discrepancy kind must be a DiscrepancyKind value.")
    if not evidence_id.strip() or not source.strip() or not rationale.strip():
        raise DomainError("Discrepancy evidence ID, source, and rationale are required.")
    if not isfinite(relative_tolerance) or relative_tolerance <= 0.0:
        raise DomainError("Discrepancy tolerance must be positive and finite.")
    if maximum_iterations < 10:
        raise DomainError("Discrepancy iteration limit must be at least ten.")
    if len(observations) < 3:
        raise DomainError("Discrepancy calibration requires at least three independent cases.")
    case_ids = tuple(item.case_id for item in observations)
    if len(case_ids) != len(set(case_ids)):
        raise DomainError("Discrepancy case identifiers must be unique.")
    applicability_ids = {item.applicability_id for item in observations}
    if len(applicability_ids) != 1:
        raise DomainError("Discrepancy observations cannot mix applicability domains.")
    values = tuple(item.ratio for item in observations)
    variances = tuple(item.standard_uncertainty**2 for item in observations)
    consensus, between_variance, evidence = _mandel_paule(
        values,
        variances,
        relative_tolerance=relative_tolerance,
        maximum_iterations=maximum_iterations,
    )
    effective_variances = tuple(max(value, evidence.variance_floor) for value in variances)
    weights = tuple(1.0 / (between_variance + value) for value in effective_variances)
    consensus_standard_uncertainty = sqrt(1.0 / sum(weights))
    predictive_variance = between_variance + consensus_standard_uncertainty**2
    warnings = [
        "Apply consensus_scale to the deterministic artifact before residual sampling.",
        "The discrepancy is valid only for the declared applicability domain; no extrapolation is implied.",
        "Predictive epistemic variance combines between-case discrepancy and consensus uncertainty.",
    ]
    if len(observations) < 6:
        warnings.append(
            "Fewer than six independent cases: random-effects uncertainty can be too narrow."
        )
    if between_variance == 0.0:
        warnings.append("Between-case discrepancy variance is on its nonnegative boundary.")
    if all(item.standard_uncertainty == 0.0 for item in observations):
        warnings.append("No referent measurement standard uncertainty was supplied.")
    return ScaleDiscrepancyCalibration(
        target_scale=normalized_target,
        component_id=normalized_component,
        kind=kind,
        applicability_id=next(iter(applicability_ids)),
        evidence_id=evidence_id.strip(),
        method="heteroscedastic-mandel-paule-predictive-discrepancy",
        source=source,
        rationale=rationale,
        observations=observations,
        consensus_scale=consensus,
        consensus_standard_uncertainty=consensus_standard_uncertainty,
        between_case_variance=between_variance,
        predictive_epistemic_variance=predictive_variance,
        evidence=evidence,
        relative_tolerance=relative_tolerance,
        maximum_iterations=maximum_iterations,
        reference_ids=DISCREPANCY_REFERENCE_IDS,
        warnings=tuple(warnings),
    )


@dataclass(frozen=True, slots=True)
class LayeredCalibrationPlan:
    """Deterministic corrections plus separated aleatory and epistemic factors.

    Deterministic consensus factors are applied once through
    :meth:`corrected_artifact`; the conditional definition then samples only the
    residual factors [1]_ [2]_ [3]_. Instances are produced by
    :func:`build_layered_calibration_plan`.

    Parameters
    ----------
    scenario_id, applicability_id : str
        Scenario and applicability domain of the plan.
    hot_fire_evidence_id : str
        Identifier of the hot-fire evidence set.
    hot_fire_calibration_hash : str
        Hash of the paired hot-fire calibration.
    discrepancy_hashes : tuple of str
        Hashes of the discrepancy calibrations.
    deterministic_scales : tuple of (str, float)
        Unique, supported scale names with positive finite products of all
        consensus factors.
    conditional_definition : ScenarioConditionalUncertainty
        Residual definition for the same scenario.
    reference_ids, warnings : tuple of str
        Source identifiers and interpretation warnings.

    Raises
    ------
    DomainError
        If the scenario does not match the definition, or the deterministic
        scales are duplicated, unsupported, nonpositive or non-finite.

    References
    ----------
    .. [1] NASA-STD-7009B, "Standard for Models and Simulations".
           https://standards.nasa.gov/standard/NASA/NASA-STD-7009
    .. [2] NASA-HDBK-7009, "NASA Handbook for Models and Simulations".
           https://ntrs.nasa.gov/citations/20140002378
    .. [3] NASA/SP-2009-569, "Bayesian Inference for NASA Probabilistic Risk and
           Reliability Analysis". https://ntrs.nasa.gov/citations/20090023159
    """

    scenario_id: str
    applicability_id: str
    hot_fire_evidence_id: str
    hot_fire_calibration_hash: str
    discrepancy_hashes: tuple[str, ...]
    deterministic_scales: tuple[tuple[str, float], ...]
    conditional_definition: ScenarioConditionalUncertainty
    reference_ids: tuple[str, ...]
    warnings: tuple[str, ...]

    def __post_init__(self) -> None:
        if self.conditional_definition.scenario_id != self.scenario_id:
            raise DomainError("Layered plan scenario does not match its conditional definition.")
        names = tuple(name for name, _value in self.deterministic_scales)
        if len(names) != len(set(names)) or set(names) - _ALLOWED_TARGET_SCALES:
            raise DomainError("Layered deterministic scale names are invalid or duplicated.")
        if any(not isfinite(value) or value <= 0.0 for _name, value in self.deterministic_scales):
            raise DomainError("Layered deterministic scales must be positive and finite.")

    @property
    def plan_hash(self) -> str:
        """SHA-256 hex digest of the canonical JSON of the complete plan.

        The digest is a reproducible provenance identity, not an
        authentication of the evidence.
        """

        return hashlib.sha256(canonical_json_bytes(self)).hexdigest()

    def corrected_artifact(self, artifact: TabulatedBurnArtifact) -> TabulatedBurnArtifact:
        """Apply all deterministic consensus corrections exactly once.

        The thrust, specific-impulse, time and axial-efficiency scales are
        applied through ``scale_tabulated_artifact``, which keeps delivered
        thrust [1]_, specific impulse and mass flow ``F / (g0 Isp)`` [2]_
        consistent. The plan hash and source calibration hashes are recorded as
        provenance [3]_.

        Parameters
        ----------
        artifact : TabulatedBurnArtifact
            Uncorrected deterministic artifact.

        Returns
        -------
        TabulatedBurnArtifact
            Corrected artifact.

        References
        ----------
        .. [1] NASA Glenn, "Rocket Thrust Equation".
               https://www1.grc.nasa.gov/beginners-guide-to-aeronautics/rocket-thrust-equation/
        .. [2] NASA Glenn, "Specific Impulse".
               https://www1.grc.nasa.gov/beginners-guide-to-aeronautics/specific-impulse/
        .. [3] NASA-STD-7009B, "Standard for Models and Simulations".
               https://standards.nasa.gov/standard/NASA/NASA-STD-7009
        """

        return scale_tabulated_artifact(
            artifact,
            dict(self.deterministic_scales),
            source_id=f"layered-calibration:{self.scenario_id}",
            provenance={
                "plan_hash": self.plan_hash,
                "hot_fire_calibration_hash": self.hot_fire_calibration_hash,
                "discrepancy_hashes": self.discrepancy_hashes,
            },
            reference_ids=self.reference_ids,
            warnings=(
                "Deterministic hot-fire and discrepancy consensus corrections applied once.",
            ),
        )


def _renamed_parameter(parameter: UncertainParameter, name: str) -> UncertainParameter:
    """Return a copy of ``parameter`` with a new name and all other fields kept."""
    return UncertainParameter(
        name=name,
        nominal=parameter.nominal,
        lower=parameter.lower,
        upper=parameter.upper,
        distribution=parameter.distribution,
        uncertainty_class=parameter.uncertainty_class,
        source=parameter.source,
        rationale=parameter.rationale,
        standard_deviation=parameter.standard_deviation,
    )


def build_layered_calibration_plan(
    calibration: HotFireMultivariateCalibration,
    discrepancies: tuple[ScaleDiscrepancyCalibration, ...],
    *,
    scenario_id: str,
    applicability_id: str,
    hot_fire_evidence_id: str,
    sample_count: int,
    source: str,
    rationale: str,
    discrepancy_correlation: CorrelationModel | None = None,
) -> LayeredCalibrationPlan:
    """Compose hot-fire variability with independent epistemic transfer evidence.

    Hot-fire population residuals remain aleatory; model-form and
    qualification-to-flight residuals remain epistemic [1]_ [2]_ [3]_. Qualified
    component names let several factors multiply the same physical scale
    without collapsing their classifications.

    Parameters
    ----------
    calibration : HotFireMultivariateCalibration
        Paired hot-fire calibration; its population residual model is used.
    discrepancies : tuple of ScaleDiscrepancyCalibration
        Discrepancy calibrations for ``applicability_id`` with unique component
        names.
    scenario_id, applicability_id : str
        Nonempty scenario and applicability identifiers.
    hot_fire_evidence_id : str
        Nonempty identifier of the hot-fire evidence; it must not be reused as
        discrepancy evidence (double counting).
    sample_count : int
        Conditional sample count.
    source, rationale : str
        Nonempty provenance and justification.
    discrepancy_correlation : CorrelationModel or None, optional
        Correlation of the discrepancy factors in their order; required when
        two discrepancies share an ``evidence_id``.

    Returns
    -------
    LayeredCalibrationPlan
        Deterministic scales and the residual conditional definition.

    Raises
    ------
    DomainError
        If identifiers are missing, applicability or evidence rules are
        violated, component names repeat, or the correlation is missing or out
        of order.

    Notes
    -----
    Each deterministic scale is the hot-fire consensus scale multiplied by the
    consensus of every discrepancy acting on it. The combined copula matrix is
    block diagonal: the hot-fire block, the discrepancy block (identity unless
    ``discrepancy_correlation`` is given) and zero cross-block correlation,
    which is an explicit independence assumption.

    References
    ----------
    .. [1] NASA-STD-7009B, "Standard for Models and Simulations".
           https://standards.nasa.gov/standard/NASA/NASA-STD-7009
    .. [2] NASA-HDBK-7009, "NASA Handbook for Models and Simulations".
           https://ntrs.nasa.gov/citations/20140002378
    .. [3] NASA/SP-2009-569, "Bayesian Inference for NASA Probabilistic Risk and
           Reliability Analysis". https://ntrs.nasa.gov/citations/20090023159
    """

    if not scenario_id.strip() or not applicability_id.strip():
        raise DomainError("Layered plan scenario and applicability IDs are required.")
    if not hot_fire_evidence_id.strip() or not source.strip() or not rationale.strip():
        raise DomainError("Layered plan evidence ID, source, and rationale are required.")
    if any(item.applicability_id != applicability_id for item in discrepancies):
        raise DomainError("Every discrepancy must match the plan applicability domain.")
    if any(item.evidence_id == hot_fire_evidence_id for item in discrepancies):
        raise DomainError("Hot-fire evidence cannot be reused as discrepancy evidence.")
    discrepancy_names = tuple(item.residual_uncertain_parameter().name for item in discrepancies)
    if len(discrepancy_names) != len(set(discrepancy_names)):
        raise DomainError("Layered discrepancy component names must be unique.")
    evidence_counts = Counter(item.evidence_id for item in discrepancies)
    if any(count > 1 for count in evidence_counts.values()) and discrepancy_correlation is None:
        raise DomainError("Shared discrepancy evidence requires an explicit correlation model.")
    if discrepancy_correlation is not None and (
        discrepancy_correlation.parameter_names != discrepancy_names
    ):
        raise DomainError("Discrepancy correlation order must match discrepancy components.")

    hot_fire_parameters, hot_fire_correlation = calibration.residual_model(
        scope=CalibrationScope.POPULATION
    )
    renamed_hot_fire = tuple(
        _renamed_parameter(
            parameter,
            qualified_scale_component_name(
                parameter.name,
                uncertainty_class=UncertaintyClass.ALEATORY,
                component_id="hot-fire",
            ),
        )
        for parameter in hot_fire_parameters
    )
    discrepancy_parameters = tuple(item.residual_uncertain_parameter() for item in discrepancies)
    parameters = renamed_hot_fire + discrepancy_parameters
    size = len(parameters)
    hot_size = len(renamed_hot_fire)
    matrix = [[0.0] * size for _ in range(size)]
    for index in range(size):
        matrix[index][index] = 1.0
    for row in range(hot_size):
        for column in range(hot_size):
            matrix[row][column] = hot_fire_correlation.matrix[row][column]
    if discrepancy_correlation is not None:
        for row in range(len(discrepancies)):
            for column in range(len(discrepancies)):
                matrix[hot_size + row][hot_size + column] = discrepancy_correlation.matrix[row][
                    column
                ]
    correlation = CorrelationModel(
        parameter_names=tuple(item.name for item in parameters),
        matrix=tuple(tuple(row) for row in matrix),
        source=(
            f"{source}; hot-fire sha256={calibration.calibration_hash}; "
            f"discrepancy sha256={','.join(item.calibration_hash for item in discrepancies)}"
        ),
        rationale=(
            "Hot-fire paired correlation and discrepancy block correlation; zero cross-block "
            f"correlation is an explicit independence assumption; {rationale}"
        ),
    )
    deterministic = {
        marginal.parameter_name: marginal.consensus_scale for marginal in calibration.marginals
    }
    for item in discrepancies:
        deterministic[item.target_scale] = (
            deterministic.get(item.target_scale, 1.0) * item.consensus_scale
        )
    definition = ScenarioConditionalUncertainty(
        scenario_id=scenario_id,
        parameters=parameters,
        sample_count=sample_count,
        source=source,
        rationale=rationale,
        correlation=correlation,
    )
    warnings = [
        "Aleatory hot-fire and epistemic discrepancy factors remain separately named.",
        "Hot-fire/discrepancy cross-block correlation is assumed zero.",
        "Apply corrected_artifact before the conditional definition; do not apply consensus factors again.",
    ]
    return LayeredCalibrationPlan(
        scenario_id=scenario_id,
        applicability_id=applicability_id,
        hot_fire_evidence_id=hot_fire_evidence_id,
        hot_fire_calibration_hash=calibration.calibration_hash,
        discrepancy_hashes=tuple(item.calibration_hash for item in discrepancies),
        deterministic_scales=tuple(sorted(deterministic.items())),
        conditional_definition=definition,
        reference_ids=tuple(
            dict.fromkeys(
                calibration.reference_ids
                + tuple(
                    reference_id for item in discrepancies for reference_id in item.reference_ids
                )
            )
        ),
        warnings=tuple(warnings),
    )
