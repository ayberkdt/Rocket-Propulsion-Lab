"""Hot-fire random-effects calibration of propulsion scale uncertainty.

Repeated hot-fire measurements of one dimensionless scale (for example measured
over predicted thrust) are described by a one-way random-effects model [1]_

    x_ij = mu + b_i + w_ij + e_ij,

with consensus correction ``mu``, persistent engine effect ``b_i`` of variance
``sigma_b**2``, burn-to-burn process deviation ``w_ij`` of variance
``sigma_w**2`` and measurement error ``e_ij`` with the declared standard
uncertainty ``u_ij`` [2]_. These three sources are kept separate rather than
collapsed into one arbitrary percentage. Engine means are combined with the
heteroscedastic Mandel-Paule procedure [3]_ [4]_.

The result is a statistical calibration of supplied data. It does not validate
the propulsion model against hardware outside the calibrated operating
condition, and it excludes the sampling uncertainty of the variance components
and model-form discrepancy [5]_ [6]_.

References
----------
.. [1] NIST/SEMATECH e-Handbook of Statistical Methods, section 7.4.4,
       "Variance Components" (one-way random-effects model).
       https://www.itl.nist.gov/div898/handbook/prc/section4/prc44.htm
.. [2] NIST TN 1297, "Guidelines for Evaluating and Expressing the Uncertainty
       of NIST Measurement Results", 1994.
       https://www.nist.gov/pml/nist-technical-note-1297
.. [3] PAULE-MANDEL-1982. R. C. Paule and J. Mandel, "Consensus Values and
       Weighting Factors", J. Res. Natl. Bur. Stand. 87(5), 377-385, 1982.
.. [4] NIST Dataplot, "Consensus Mean" (heteroscedastic random effects and
       Mandel-Paule estimation).
       https://www.itl.nist.gov/div898/software/dataplot/refman1/auxillar/consmean.htm
.. [5] NASA/SP-2009-569, "Bayesian Inference for NASA Probabilistic Risk and
       Reliability Analysis". https://ntrs.nasa.gov/citations/20090023159
.. [6] NASA-STD-7009B, "Standard for Models and Simulations".
       https://standards.nasa.gov/standard/NASA/NASA-STD-7009
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

    Attributes
    ----------
    POPULATION
        Predicts a burn of a randomly selected engine; the residual variance is
        ``sigma_w**2 + sigma_b**2`` (within-engine plus between-engine process
        variance) [1]_.
    SAME_ENGINE
        Predicts another burn of an engine whose persistent effect is already
        known and corrected; the residual variance is ``sigma_w**2`` only.

    References
    ----------
    .. [1] NIST/SEMATECH e-Handbook of Statistical Methods, section 7.4.4,
           "Variance Components" (one-way random-effects model).
           https://www.itl.nist.gov/div898/handbook/prc/section4/prc44.htm
    """

    POPULATION = "population"
    SAME_ENGINE = "same-engine"


@dataclass(frozen=True, slots=True)
class HotFireScaleMeasurement:
    """One positive measured/model scale with its standard measurement uncertainty.

    Parameters
    ----------
    engine_id, run_id : str
        Nonempty engine and run identifiers.
    condition_id : str
        Nonempty operating-condition stratum. Measurements pooled in one
        calibration must share it: pressure, mixture ratio, throttle, thermal
        state or instrumentation changes belong in separate strata unless an
        explicit regression model removes them [2]_.
    value : float
        Positive finite measured/model scale.
    standard_uncertainty : float
        Finite, nonnegative standard measurement uncertainty of ``value`` in
        the same units [1]_.
    source : str
        Nonempty provenance.

    Raises
    ------
    DomainError
        If any constraint above is violated.

    References
    ----------
    .. [1] NIST TN 1297, "Guidelines for Evaluating and Expressing the Uncertainty
           of NIST Measurement Results", 1994.
           https://www.nist.gov/pml/nist-technical-note-1297
    .. [2] NASA-STD-7009B, "Standard for Models and Simulations".
           https://standards.nasa.gov/standard/NASA/NASA-STD-7009
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

    Parameters
    ----------
    engine_id : str
        Engine identifier.
    observation_count : int
        Number of burns ``n_i``.
    mean_scale, minimum_scale, maximum_scale : float
        Engine mean ``x_i`` and extrema.
    sample_variance : float
        ``n_i - 1`` sample variance; zero for a single burn.
    measurement_variance_of_mean : float
        ``sum_j u_ij**2 / n_i**2``.
    modeled_variance_of_mean : float
        ``v_i = sigma_w**2 / n_i + sum_j u_ij**2 / n_i**2``, the within-engine
        variance of the engine mean used as the Mandel-Paule input [1]_.

    References
    ----------
    .. [1] NIST Dataplot, "Consensus Mean" (heteroscedastic random effects and
           Mandel-Paule estimation).
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
    """Deterministic solver evidence for the Mandel-Paule between-group variance.

    Parameters
    ----------
    target_degrees_of_freedom : int
        ``k - 1`` for ``k`` groups, the target of the Mandel-Paule equation
        [1]_ [2]_.
    iterations : int
        Bracket expansions plus bisection steps.
    equation_value : float
        Weighted sum of squares at the returned variance.
    equation_residual : float
        ``|equation_value - (k - 1)|``; zero for a boundary solution.
    relative_tolerance : float
        Acceptance tolerance relative to ``max(k - 1, 1)``.
    variance_floor : float
        Lower floor ``(1e-12 max(max |x_i|, 1))**2`` applied to each
        within-group variance so that zero-variance groups keep finite weight.
    boundary_solution : bool
        ``True`` when the weighted sum of squares at zero between-group
        variance already falls at or below ``k - 1``; the estimate is then zero.
    converged : bool
        ``True`` for every returned record.

    References
    ----------
    .. [1] PAULE-MANDEL-1982. R. C. Paule and J. Mandel, "Consensus Values and
           Weighting Factors", J. Res. Natl. Bur. Stand. 87(5), 377-385, 1982.
    .. [2] NIST Dataplot, "Consensus Mean" (heteroscedastic random effects and
           Mandel-Paule estimation).
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
    corrected mean [1]_. Declared measurement uncertainty [2]_ is reported and
    removed from the pooled within-engine process variance; it is not relabelled
    as engine variability. Variance-component sampling uncertainty and
    model-form discrepancy are not included [7]_ [8]_. Instances are produced
    by :func:`calibrate_hot_fire_scale`.

    Parameters
    ----------
    parameter_name : str
        Calibrated scale name.
    condition_id : str
        Common operating-condition stratum of all measurements.
    method : str
        Method identifier.
    source, rationale : str
        Provenance and justification.
    measurements : tuple of HotFireScaleMeasurement
        Input data.
    engine_summaries : tuple of EngineCalibrationSummary
        Per-engine statistics.
    observation_count, engine_count : int
        ``N`` burns on ``k`` engines.
    within_degrees_of_freedom : int
        ``N - k``.
    consensus_scale : float
        Mandel-Paule weighted mean ``mu = sum_i w_i x_i / sum_i w_i`` with
        ``w_i = 1 / (sigma_b**2 + v_i)`` [3]_ [4]_ [5]_.
    consensus_standard_uncertainty : float
        Inverse-variance weighted-mean standard uncertainty
        ``(sum_i w_i)**-0.5``. This is not the Rukhin-Vangel interval reported
        by [4]_.
    observed_within_engine_variance : float
        Pooled within-engine sample variance.
    measurement_within_variance_contribution : float
        Expected measurement contribution to that pooled variance.
    within_engine_process_variance : float
        ``sigma_w**2``, nonnegative.
    between_engine_process_variance : float
        ``sigma_b**2``, nonnegative.
    population_process_variance : float
        ``sigma_w**2 + sigma_b**2`` [6]_.
    within_engine_fraction, between_engine_fraction : float
        Shares of the population variance; both zero when it is zero.
    evidence : MandelPauleEvidence
        Solver evidence.
    reference_ids, warnings : tuple of str
        Source identifiers and interpretation warnings.

    References
    ----------
    .. [1] NIST/SEMATECH e-Handbook of Statistical Methods, section 7.4.4,
           "Variance Components" (one-way random-effects model).
           https://www.itl.nist.gov/div898/handbook/prc/section4/prc44.htm
    .. [2] NIST TN 1297, "Guidelines for Evaluating and Expressing the Uncertainty
           of NIST Measurement Results", 1994.
           https://www.nist.gov/pml/nist-technical-note-1297
    .. [3] PAULE-MANDEL-1982. R. C. Paule and J. Mandel, "Consensus Values and
           Weighting Factors", J. Res. Natl. Bur. Stand. 87(5), 377-385, 1982.
    .. [4] NIST Dataplot, "Consensus Mean" (heteroscedastic random effects and
           Mandel-Paule estimation).
           https://www.itl.nist.gov/div898/software/dataplot/refman1/auxillar/consmean.htm
    .. [5] NIST-WEIGHTED-MEAN. NIST Dataplot, "Weighted Mean", equation 2-20,
           normalized weighted mean.
    .. [6] MIT-TOTAL-VARIANCE. MIT RES.6-012, lecture 13.7, conditional and
           between-group variance identity (law of total variance).
    .. [7] NASA/SP-2009-569, "Bayesian Inference for NASA Probabilistic Risk and
           Reliability Analysis". https://ntrs.nasa.gov/citations/20090023159
    .. [8] NASA-STD-7009B, "Standard for Models and Simulations".
           https://standards.nasa.gov/standard/NASA/NASA-STD-7009
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
        """SHA-256 hex digest of the canonical JSON of the complete record.

        The digest covers data, method, evidence and sources; it is a
        reproducible content identity, not an authentication of the data.

        Returns
        -------
        str
            Lower-case hexadecimal digest.
        """

        return hashlib.sha256(canonical_json_bytes(self)).hexdigest()

    def residual_uncertain_parameter(
        self,
        *,
        scope: CalibrationScope = CalibrationScope.POPULATION,
        coverage_factor: float = 3.0,
    ) -> UncertainParameter:
        """Return a bounded residual scale for the conditional ensemble layer.

        Apply ``consensus_scale`` to the deterministic base artifact first. The
        returned truncated-normal parameter is centered at one with relative
        standard deviation ``r = sigma / mu``, where ``sigma**2`` is the residual
        variance of ``scope`` [1]_, and bounds ``[max(1e-12, 1 - k r), 1 + k r]``
        for coverage factor ``k`` [2]_. It represents remaining physical
        variability and excludes the sampling uncertainty of the estimated mean
        and variance components and model-form discrepancy.

        Parameters
        ----------
        scope : CalibrationScope, optional
            Population or same-engine residual variance.
        coverage_factor : float, optional
            Positive finite truncation half-width in standard deviations.

        Returns
        -------
        UncertainParameter
            Aleatory truncated-normal scale named ``parameter_name``.

        Raises
        ------
        DomainError
            If ``coverage_factor`` or ``scope`` is invalid, or the selected
            residual variance is zero.

        References
        ----------
        .. [1] NIST/SEMATECH e-Handbook of Statistical Methods, section 7.4.4,
               "Variance Components" (one-way random-effects model).
               https://www.itl.nist.gov/div898/handbook/prc/section4/prc44.htm
        .. [2] NIST TN 1297, "Guidelines for Evaluating and Expressing the Uncertainty
               of NIST Measurement Results", 1994.
               https://www.nist.gov/pml/nist-technical-note-1297
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
    """Group measurements by engine, in order of first appearance."""
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
    """Solve the Mandel-Paule equation for the between-group variance.

    For group means ``x_i`` with within-group variances of the mean ``v_i``, the
    between-group variance ``tau**2 >= 0`` satisfies [1]_ [2]_

        sum_i w_i (x_i - mu)**2 = k - 1,
        w_i = 1 / (tau**2 + v_i),  mu = sum_i w_i x_i / sum_i w_i.

    If the left side at ``tau**2 = 0`` is already at most ``k - 1`` the boundary
    solution ``tau**2 = 0`` is returned. Otherwise the root is bracketed by
    repeated fourfold expansion and found by bisection on the monotonically
    decreasing left side.

    Parameters
    ----------
    means : tuple of float
        Group means ``x_i``.
    variances_of_mean : tuple of float
        Nonnegative ``v_i``; each is floored at a small scale-relative value.
    relative_tolerance : float
        Acceptance tolerance relative to ``max(k - 1, 1)``.
    maximum_iterations : int
        Combined limit on bracket expansions and bisection steps.

    Returns
    -------
    tuple
        ``(mu, tau**2, evidence)``.

    Raises
    ------
    DomainError
        If the bracket or the bisection does not converge within the limit.

    References
    ----------
    .. [1] PAULE-MANDEL-1982. R. C. Paule and J. Mandel, "Consensus Values and
           Weighting Factors", J. Res. Natl. Bur. Stand. 87(5), 377-385, 1982.
    .. [2] NIST Dataplot, "Consensus Mean" (heteroscedastic random effects and
           Mandel-Paule estimation).
           https://www.itl.nist.gov/div898/software/dataplot/refman1/auxillar/consmean.htm
    """
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
    """Calibrate a consensus correction and residual variances from hot-fire data.

    Parameters
    ----------
    parameter_name : str
        One of ``thrust_scale``, ``specific_impulse_scale``, ``time_scale`` or
        ``axial_efficiency_scale``.
    measurements : tuple of HotFireScaleMeasurement
        At least three observations on at least two engines, with unique
        engine/run identities, one common ``condition_id`` and at least one
        repeated engine.
    source, rationale : str
        Nonempty provenance and justification.
    relative_tolerance : float, optional
        Positive tolerance of the Mandel-Paule solve.
    maximum_iterations : int, optional
        Iteration limit of the Mandel-Paule solve, at least ten.

    Returns
    -------
    HotFireScaleCalibration
        Consensus correction, variance components, per-engine summaries,
        solver evidence and warnings.

    Raises
    ------
    DomainError
        If the inputs violate the constraints above or the variance solve does
        not converge.

    Notes
    -----
    For engine ``i`` with ``n_i`` burns, mean ``x_i`` and declared measurement
    uncertainties ``u_ij``, and ``N`` burns on ``k`` engines:

    1. Pooled within-engine variance [1]_:
       ``s_w**2 = sum_i sum_j (x_ij - x_i)**2 / (N - k)``.
    2. Expected contribution of independent measurement errors with declared
       standard uncertainties ``u_ij`` [2]_ to ``s_w**2``:
       ``s_e**2 = sum_i (1 - 1/n_i) sum_j u_ij**2 / (N - k)``.
    3. Within-engine process variance ``sigma_w**2 = s_w**2 - s_e**2``, set to
       zero when it is not above ``1e-12`` times the larger of the two.
    4. Variance of each engine mean
       ``v_i = sigma_w**2 / n_i + sum_j u_ij**2 / n_i**2``.
    5. Mandel-Paule between-engine variance ``sigma_b**2`` and weighted
       consensus ``mu = sum_i w_i x_i / sum_i w_i`` with
       ``w_i = 1 / (sigma_b**2 + v_i)`` [3]_ [4]_ [5]_, reported with the
       inverse-variance standard uncertainty ``(sum_i w_i)**-0.5``.

    Engines may have unequal replicate counts and unequal measurement
    uncertainties. With fewer than six engines the Mandel-Paule uncertainty can
    be too narrow [4]_; the result carries a warning. The population residual
    variance is ``sigma_w**2 + sigma_b**2`` [6]_. Variance-component sampling
    uncertainty and model-form discrepancy are not included [7]_ [8]_.

    References
    ----------
    .. [1] NIST/SEMATECH e-Handbook of Statistical Methods, section 7.4.4,
           "Variance Components" (one-way random-effects model).
           https://www.itl.nist.gov/div898/handbook/prc/section4/prc44.htm
    .. [2] NIST TN 1297, "Guidelines for Evaluating and Expressing the Uncertainty
           of NIST Measurement Results", 1994.
           https://www.nist.gov/pml/nist-technical-note-1297
    .. [3] PAULE-MANDEL-1982. R. C. Paule and J. Mandel, "Consensus Values and
           Weighting Factors", J. Res. Natl. Bur. Stand. 87(5), 377-385, 1982.
    .. [4] NIST Dataplot, "Consensus Mean" (heteroscedastic random effects and
           Mandel-Paule estimation).
           https://www.itl.nist.gov/div898/software/dataplot/refman1/auxillar/consmean.htm
    .. [5] NIST-WEIGHTED-MEAN. NIST Dataplot, "Weighted Mean", equation 2-20,
           normalized weighted mean.
    .. [6] MIT-TOTAL-VARIANCE. MIT RES.6-012, lecture 13.7, conditional and
           between-group variance identity (law of total variance).
    .. [7] NASA/SP-2009-569, "Bayesian Inference for NASA Probabilistic Risk and
           Reliability Analysis". https://ntrs.nasa.gov/citations/20090023159
    .. [8] NASA-STD-7009B, "Standard for Models and Simulations".
           https://standards.nasa.gov/standard/NASA/NASA-STD-7009
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

