"""Hot-fire random-effects calibration for propulsion scale uncertainty.

The model separates three quantities that must not be collapsed into one
arbitrary percentage: between-engine variability, within-engine burn-to-burn
variability, and declared measurement standard uncertainty. Engine means are
combined with the heteroscedastic Mandel-Paule random-effects procedure.

References
----------
NIST TN 1297: https://www.nist.gov/pml/nist-technical-note-1297
NIST/SEMATECH variance components:
https://www.itl.nist.gov/div898/handbook/prc/section4/prc44.htm
NIST Dataplot consensus mean:
https://www.itl.nist.gov/div898/software/dataplot/refman1/auxillar/consmean.htm
NASA-SP-2009-569: https://ntrs.nasa.gov/citations/20090023159
NASA-STD-7009B: https://standards.nasa.gov/standard/NASA/NASA-STD-7009
"""

from __future__ import annotations

import hashlib
from dataclasses import dataclass
from enum import Enum
from math import isfinite, sqrt
from statistics import fmean

from rocket_propulsion.core.errors import DomainError

from .serialization import canonical_json_bytes
from .uncertainty import (
    DistributionKind,
    UncertainParameter,
    UncertaintyClass,
)

HOT_FIRE_CALIBRATION_REFERENCE_IDS = (
    "NIST-TN-1297",
    "NIST-SEMATECH-7.4.4",
    "NIST-DATAPLOT-CONSENSUS-MEAN",
    "NASA-SP-2009-569",
    "NASA-STD-7009B",
)

_SUPPORTED_SCALE_NAMES = {
    "thrust_scale",
    "specific_impulse_scale",
    "time_scale",
    "axial_efficiency_scale",
}


class CalibrationScope(str, Enum):
    """Physical population represented by a calibrated residual distribution.

    ``POPULATION`` predicts a burn from a randomly selected engine and includes
    both between-engine and within-engine process variance. ``SAME_ENGINE``
    predicts another burn from an engine whose persistent engine effect is
    already known/corrected and includes only within-engine process variance.

    References
    ----------
    NIST/SEMATECH variance components:
    https://www.itl.nist.gov/div898/handbook/prc/section4/prc44.htm
    """

    POPULATION = "population"
    SAME_ENGINE = "same-engine"


@dataclass(frozen=True, slots=True)
class HotFireScaleMeasurement:
    """One positive measured/model scale with standard measurement uncertainty.

    Measurements pooled in one calibration must use the same ``condition_id``;
    pressure, mixture ratio, throttle, thermal state, or instrumentation changes
    belong in separate strata unless an explicit regression model removes them.

    References
    ----------
    NIST TN 1297: https://www.nist.gov/pml/nist-technical-note-1297
    NASA-STD-7009B: https://standards.nasa.gov/standard/NASA/NASA-STD-7009
    """

    engine_id: str
    run_id: str
    condition_id: str
    value: float
    standard_uncertainty: float
    source: str

    def __post_init__(self) -> None:
        if not self.engine_id.strip() or not self.run_id.strip():
            raise DomainError("Hot-fire engine and run identifiers must not be empty.")
        if not self.condition_id.strip():
            raise DomainError("Hot-fire operating-condition identifier is required.")
        if not isfinite(self.value) or self.value <= 0.0:
            raise DomainError("Hot-fire scale measurement must be positive and finite.")
        if not isfinite(self.standard_uncertainty) or self.standard_uncertainty < 0.0:
            raise DomainError(
                "Hot-fire measurement standard uncertainty must be finite and nonnegative."
            )
        if not self.source.strip():
            raise DomainError("Hot-fire measurement source is required.")


@dataclass(frozen=True, slots=True)
class EngineCalibrationSummary:
    """Per-engine sufficient statistics retained by the calibration.

    References
    ----------
    NIST Dataplot consensus mean:
    https://www.itl.nist.gov/div898/software/dataplot/refman1/auxillar/consmean.htm
    """

    engine_id: str
    observation_count: int
    mean_scale: float
    minimum_scale: float
    maximum_scale: float
    sample_variance: float
    measurement_variance_of_mean: float
    modeled_variance_of_mean: float


@dataclass(frozen=True, slots=True)
class MandelPauleEvidence:
    """Deterministic solver evidence for the between-engine variance estimate.

    References
    ----------
    NIST Dataplot consensus mean:
    https://www.itl.nist.gov/div898/software/dataplot/refman1/auxillar/consmean.htm
    """

    target_degrees_of_freedom: int
    iterations: int
    equation_value: float
    equation_residual: float
    relative_tolerance: float
    variance_floor: float
    boundary_solution: bool
    converged: bool


@dataclass(frozen=True, slots=True)
class HotFireScaleCalibration:
    """Traceable random-effects calibration of one propulsion scale.

    The consensus scale is a correction to the deterministic base model. The
    variance components describe residual physical variability about that
    corrected mean. Measurement uncertainty is reported and removed from the
    pooled within-engine process variance; it is not silently re-labelled as
    engine variability.

    References
    ----------
    NIST TN 1297: https://www.nist.gov/pml/nist-technical-note-1297
    NIST/SEMATECH variance components:
    https://www.itl.nist.gov/div898/handbook/prc/section4/prc44.htm
    NIST Dataplot consensus mean:
    https://www.itl.nist.gov/div898/software/dataplot/refman1/auxillar/consmean.htm
    NASA-SP-2009-569: https://ntrs.nasa.gov/citations/20090023159
    NASA-STD-7009B: https://standards.nasa.gov/standard/NASA/NASA-STD-7009
    """

    parameter_name: str
    condition_id: str
    method: str
    source: str
    rationale: str
    measurements: tuple[HotFireScaleMeasurement, ...]
    engine_summaries: tuple[EngineCalibrationSummary, ...]
    observation_count: int
    engine_count: int
    within_degrees_of_freedom: int
    consensus_scale: float
    consensus_standard_uncertainty: float
    observed_within_engine_variance: float
    measurement_within_variance_contribution: float
    within_engine_process_variance: float
    between_engine_process_variance: float
    population_process_variance: float
    within_engine_fraction: float
    between_engine_fraction: float
    evidence: MandelPauleEvidence
    reference_ids: tuple[str, ...]
    warnings: tuple[str, ...]

    @property
    def calibration_hash(self) -> str:
        """SHA-256 identity covering data, method, evidence, and sources.

        References
        ----------
        NASA-STD-7009B: https://standards.nasa.gov/standard/NASA/NASA-STD-7009
        """

        return hashlib.sha256(canonical_json_bytes(self)).hexdigest()

    def residual_uncertain_parameter(
        self,
        *,
        scope: CalibrationScope = CalibrationScope.POPULATION,
        coverage_factor: float = 3.0,
    ) -> UncertainParameter:
        """Return a bounded residual scale for the conditional ensemble layer.

        Apply ``consensus_scale`` to the deterministic base artifact first.
        The returned parameter is centered at one and represents remaining
        physical variability. It deliberately excludes sampling uncertainty in
        the estimated mean/variance components and model-form discrepancy.

        References
        ----------
        NIST TN 1297: https://www.nist.gov/pml/nist-technical-note-1297
        NIST/SEMATECH variance components:
        https://www.itl.nist.gov/div898/handbook/prc/section4/prc44.htm
        """

        if not isfinite(coverage_factor) or coverage_factor <= 0.0:
            raise DomainError("Calibration coverage factor must be positive and finite.")
        if not isinstance(scope, CalibrationScope):
            raise DomainError("Calibration scope must be a CalibrationScope value.")
        variance = (
            self.population_process_variance
            if scope is CalibrationScope.POPULATION
            else self.within_engine_process_variance
        )
        relative_standard_deviation = sqrt(variance) / self.consensus_scale
        if relative_standard_deviation <= 0.0:
            raise DomainError(
                "Calibrated residual variance is zero; no stochastic parameter is needed."
            )
        lower = max(
            1e-12,
            1.0 - coverage_factor * relative_standard_deviation,
        )
        upper = 1.0 + coverage_factor * relative_standard_deviation
        return UncertainParameter(
            name=self.parameter_name,
            nominal=1.0,
            lower=lower,
            upper=upper,
            distribution=DistributionKind.TRUNCATED_NORMAL,
            uncertainty_class=UncertaintyClass.ALEATORY,
            source=(
                f"{self.source}; hot-fire calibration sha256={self.calibration_hash}"
            ),
            rationale=(
                f"{scope.value} residual variability after applying consensus "
                f"correction scale {self.consensus_scale:.17g}; "
                f"{self.rationale}"
            ),
            standard_deviation=relative_standard_deviation,
        )


def _group_measurements(
    measurements: tuple[HotFireScaleMeasurement, ...],
) -> tuple[tuple[str, tuple[HotFireScaleMeasurement, ...]], ...]:
    engine_ids = tuple(dict.fromkeys(item.engine_id for item in measurements))
    return tuple(
        (
            engine_id,
            tuple(item for item in measurements if item.engine_id == engine_id),
        )
        for engine_id in engine_ids
    )


def _mandel_paule(
    means: tuple[float, ...],
    variances_of_mean: tuple[float, ...],
    *,
    relative_tolerance: float,
    maximum_iterations: int,
) -> tuple[float, float, MandelPauleEvidence]:
    target = len(means) - 1
    reference_scale = max(max(abs(value) for value in means), 1.0)
    variance_floor = (reference_scale * 1e-12) ** 2
    effective_variances = tuple(
        max(value, variance_floor) for value in variances_of_mean
    )

    def evaluate(between_variance: float) -> tuple[float, float]:
        weights = tuple(
            1.0 / (between_variance + value) for value in effective_variances
        )
        consensus = sum(
            weight * value for weight, value in zip(weights, means, strict=True)
        ) / sum(weights)
        equation = sum(
            weight * (value - consensus) ** 2
            for weight, value in zip(weights, means, strict=True)
        )
        return consensus, equation

    consensus_at_zero, equation_at_zero = evaluate(0.0)
    if equation_at_zero <= target:
        evidence = MandelPauleEvidence(
            target_degrees_of_freedom=target,
            iterations=0,
            equation_value=equation_at_zero,
            equation_residual=0.0,
            relative_tolerance=relative_tolerance,
            variance_floor=variance_floor,
            boundary_solution=True,
            converged=True,
        )
        return consensus_at_zero, 0.0, evidence

    mean_of_means = fmean(means)
    high = max(
        sum((value - mean_of_means) ** 2 for value in means) / target,
        max(effective_variances),
    )
    _consensus, high_equation = evaluate(high)
    expansion_count = 0
    while high_equation > target and expansion_count < maximum_iterations:
        high *= 4.0
        _consensus, high_equation = evaluate(high)
        expansion_count += 1
    if high_equation > target:
        raise DomainError("Mandel-Paule variance bracket did not converge.")

    low = 0.0
    iterations = expansion_count
    consensus = consensus_at_zero
    equation = equation_at_zero
    for _ in range(maximum_iterations - expansion_count):
        between_variance = 0.5 * (low + high)
        consensus, equation = evaluate(between_variance)
        iterations += 1
        residual = abs(equation - target)
        if residual <= relative_tolerance * max(float(target), 1.0):
            evidence = MandelPauleEvidence(
                target_degrees_of_freedom=target,
                iterations=iterations,
                equation_value=equation,
                equation_residual=residual,
                relative_tolerance=relative_tolerance,
                variance_floor=variance_floor,
                boundary_solution=False,
                converged=True,
            )
            return consensus, between_variance, evidence
        if equation > target:
            low = between_variance
        else:
            high = between_variance
    raise DomainError("Mandel-Paule variance solve exceeded its iteration limit.")


def calibrate_hot_fire_scale(
    parameter_name: str,
    measurements: tuple[HotFireScaleMeasurement, ...],
    *,
    source: str,
    rationale: str,
    relative_tolerance: float = 1e-10,
    maximum_iterations: int = 200,
) -> HotFireScaleCalibration:
    """Calibrate correction and residual variance from repeated hot-fire data.

    A common within-engine process variance is estimated after subtracting the
    contribution of each observation's declared measurement standard
    uncertainty. Engine means may have unequal replicate counts and unequal
    measurement uncertainties. Mandel-Paule then estimates the nonnegative
    between-engine variance and consensus correction.

    References
    ----------
    NIST TN 1297: https://www.nist.gov/pml/nist-technical-note-1297
    NIST/SEMATECH variance components:
    https://www.itl.nist.gov/div898/handbook/prc/section4/prc44.htm
    NIST Dataplot consensus mean:
    https://www.itl.nist.gov/div898/software/dataplot/refman1/auxillar/consmean.htm
    NASA-SP-2009-569: https://ntrs.nasa.gov/citations/20090023159
    NASA-STD-7009B: https://standards.nasa.gov/standard/NASA/NASA-STD-7009
    """

    normalized_name = parameter_name.strip()
    if normalized_name not in _SUPPORTED_SCALE_NAMES:
        raise DomainError("Unsupported hot-fire propulsion scale parameter.")
    if not source.strip() or not rationale.strip():
        raise DomainError("Hot-fire calibration source and rationale are required.")
    if not isfinite(relative_tolerance) or relative_tolerance <= 0.0:
        raise DomainError("Calibration tolerance must be positive and finite.")
    if maximum_iterations < 10:
        raise DomainError("Calibration iteration limit must be at least ten.")
    if len(measurements) < 3:
        raise DomainError("Hot-fire calibration requires at least three observations.")
    identities = tuple((item.engine_id, item.run_id) for item in measurements)
    if len(identities) != len(set(identities)):
        raise DomainError("Hot-fire engine/run identities must be unique.")
    condition_ids = {item.condition_id for item in measurements}
    if len(condition_ids) != 1:
        raise DomainError(
            "Hot-fire measurements from different operating conditions cannot be pooled."
        )
    groups = _group_measurements(measurements)
    if len(groups) < 2:
        raise DomainError("Hot-fire calibration requires at least two engines.")
    within_degrees_of_freedom = len(measurements) - len(groups)
    if within_degrees_of_freedom <= 0:
        raise DomainError(
            "At least one engine requires a repeated hot-fire observation."
        )

    raw_summaries: list[tuple[str, int, float, float, float, float, float]] = []
    within_sum_squares = 0.0
    measurement_within_sum_squares = 0.0
    for engine_id, items in groups:
        count = len(items)
        mean = fmean(item.value for item in items)
        sum_squares = sum((item.value - mean) ** 2 for item in items)
        sample_variance = sum_squares / (count - 1) if count > 1 else 0.0
        measurement_variance_of_mean = (
            sum(item.standard_uncertainty**2 for item in items) / count**2
        )
        within_sum_squares += sum_squares
        measurement_within_sum_squares += (1.0 - 1.0 / count) * sum(
            item.standard_uncertainty**2 for item in items
        )
        raw_summaries.append(
            (
                engine_id,
                count,
                mean,
                min(item.value for item in items),
                max(item.value for item in items),
                sample_variance,
                measurement_variance_of_mean,
            )
        )

    observed_within_variance = within_sum_squares / within_degrees_of_freedom
    measurement_contribution = (
        measurement_within_sum_squares / within_degrees_of_freedom
    )
    within_residual = observed_within_variance - measurement_contribution
    variance_scale = max(
        observed_within_variance,
        measurement_contribution,
        1e-30,
    )
    within_process_variance = (
        0.0
        if within_residual <= 1e-12 * variance_scale
        else within_residual
    )
    means = tuple(item[2] for item in raw_summaries)
    modeled_variances_of_mean = tuple(
        within_process_variance / item[1] + item[6] for item in raw_summaries
    )
    consensus, between_variance, evidence = _mandel_paule(
        means,
        modeled_variances_of_mean,
        relative_tolerance=relative_tolerance,
        maximum_iterations=maximum_iterations,
    )
    effective_variances = tuple(
        max(value, evidence.variance_floor) for value in modeled_variances_of_mean
    )
    weights = tuple(
        1.0 / (between_variance + value) for value in effective_variances
    )
    consensus_standard_uncertainty = sqrt(1.0 / sum(weights))
    summaries = tuple(
        EngineCalibrationSummary(
            engine_id=item[0],
            observation_count=item[1],
            mean_scale=item[2],
            minimum_scale=item[3],
            maximum_scale=item[4],
            sample_variance=item[5],
            measurement_variance_of_mean=item[6],
            modeled_variance_of_mean=modeled,
        )
        for item, modeled in zip(raw_summaries, modeled_variances_of_mean, strict=True)
    )
    population_variance = within_process_variance + between_variance
    if population_variance > 0.0:
        within_fraction = within_process_variance / population_variance
        between_fraction = between_variance / population_variance
    else:
        within_fraction = 0.0
        between_fraction = 0.0
    warnings = [
        "Variance-component sampling uncertainty and model-form discrepancy are not included.",
        "Apply consensus_scale to the deterministic base model before residual sampling.",
    ]
    if len(groups) < 6:
        warnings.append(
            "Fewer than six engines: NIST notes Mandel-Paule uncertainty intervals can be too narrow."
        )
    if within_process_variance == 0.0 and observed_within_variance > 0.0:
        warnings.append(
            "Declared measurement uncertainty consumed the observed within-engine variance."
        )
    if between_variance == 0.0:
        warnings.append("Between-engine variance estimate is on its nonnegative boundary.")
    if all(item.standard_uncertainty == 0.0 for item in measurements):
        warnings.append("No Type B measurement standard uncertainty was supplied.")
    return HotFireScaleCalibration(
        parameter_name=normalized_name,
        condition_id=next(iter(condition_ids)),
        method="heteroscedastic-mandel-paule-random-effects",
        source=source,
        rationale=rationale,
        measurements=measurements,
        engine_summaries=summaries,
        observation_count=len(measurements),
        engine_count=len(groups),
        within_degrees_of_freedom=within_degrees_of_freedom,
        consensus_scale=consensus,
        consensus_standard_uncertainty=consensus_standard_uncertainty,
        observed_within_engine_variance=observed_within_variance,
        measurement_within_variance_contribution=measurement_contribution,
        within_engine_process_variance=within_process_variance,
        between_engine_process_variance=between_variance,
        population_process_variance=population_variance,
        within_engine_fraction=within_fraction,
        between_engine_fraction=between_fraction,
        evidence=evidence,
        reference_ids=HOT_FIRE_CALIBRATION_REFERENCE_IDS,
        warnings=tuple(warnings),
    )

