"""Paired hot-fire calibration of cross-parameter propulsion uncertainty.

The scalar calibration module estimates each marginal correction and variance
component.  This module adds the information that separate fits cannot retain:
within-engine and between-engine cross-covariance from measurements made on the
same engine/run.  Declared measurement covariance is removed before physical
correlation is estimated.

References
----------
JCGM 100:2008, sections 5.2 and C.3.6:
https://www.bipm.org/documents/20126/2071204/JCGM_100_2008_E.pdf
NIST TN 1297, section 5:
https://www.nist.gov/pml/nist-technical-note-1297/nist-tn-1297-5-combined-standard-uncertainty
NIST/SEMATECH variance components:
https://www.itl.nist.gov/div898/handbook/prc/section4/prc44.htm
NASA-STD-7009B: https://standards.nasa.gov/standard/NASA/NASA-STD-7009
"""

from __future__ import annotations

import hashlib
import json
from dataclasses import dataclass
from math import isfinite, sqrt
from statistics import fmean
from typing import Any

from rocket_propulsion.core.errors import DomainError, InputError

from .calibration import (
    CalibrationScope,
    HotFireScaleCalibration,
    HotFireScaleMeasurement,
    calibrate_hot_fire_scale,
)
from .serialization import canonical_json_bytes
from .uncertainty import CorrelationModel, UncertainParameter

HOT_FIRE_MULTIVARIATE_SCHEMA = "rocket_propulsion_hot_fire_multivariate_calibration_v1"
HOT_FIRE_MULTIVARIATE_REFERENCE_IDS = (
    "JCGM-100-2008",
    "NIST-TN-1297",
    "NIST-SEMATECH-7.4.4",
    "NASA-STD-7009B",
)

_MATRIX_TOLERANCE = 1e-12


def _matrix(values: Any, *, size: int, field: str) -> tuple[tuple[float, ...], ...]:
    try:
        matrix = tuple(tuple(float(item) for item in row) for row in values)
    except (TypeError, ValueError) as error:
        raise DomainError(f"{field} must be a finite square matrix.") from error
    if len(matrix) != size or any(len(row) != size for row in matrix):
        raise DomainError(f"{field} dimensions must match parameter_names.")
    if any(not isfinite(item) for row in matrix for item in row):
        raise DomainError(f"{field} must contain only finite values.")
    return matrix


def _is_positive_semidefinite(matrix: tuple[tuple[float, ...], ...]) -> bool:
    """Check symmetric positive semidefiniteness with a pivoted-tolerance LDL form."""

    size = len(matrix)
    scale = max(max(abs(value) for value in row) for row in matrix) if size else 1.0
    tolerance = _MATRIX_TOLERANCE * max(scale, 1.0)
    lower = [[0.0] * size for _ in range(size)]
    diagonal = [0.0] * size
    for column in range(size):
        pivot = matrix[column][column] - sum(
            lower[column][prior] ** 2 * diagonal[prior] for prior in range(column)
        )
        if pivot < -tolerance:
            return False
        diagonal[column] = 0.0 if abs(pivot) <= tolerance else pivot
        lower[column][column] = 1.0
        for row in range(column + 1, size):
            residual = matrix[row][column] - sum(
                lower[row][prior] * lower[column][prior] * diagonal[prior]
                for prior in range(column)
            )
            if diagonal[column] == 0.0:
                if abs(residual) > tolerance:
                    return False
                lower[row][column] = 0.0
            else:
                lower[row][column] = residual / diagonal[column]
    return True


def _is_positive_definite(matrix: tuple[tuple[float, ...], ...]) -> bool:
    size = len(matrix)
    lower = [[0.0] * size for _ in range(size)]
    for row in range(size):
        for column in range(row + 1):
            residual = matrix[row][column] - sum(
                lower[row][prior] * lower[column][prior] for prior in range(column)
            )
            if row == column:
                if residual <= _MATRIX_TOLERANCE:
                    return False
                lower[row][column] = sqrt(residual)
            else:
                lower[row][column] = residual / lower[column][column]
    return True


@dataclass(frozen=True, slots=True)
class HotFireVectorMeasurement:
    """One paired vector measurement and its standard-uncertainty covariance.

    Parameter order is part of the data contract.  The covariance matrix is in
    squared scale units and may include off-diagonal instrument or reduction
    covariance.  Errors from different runs are assumed independent.

    References
    ----------
    JCGM 100:2008, section 5.2:
    https://www.bipm.org/documents/20126/2071204/JCGM_100_2008_E.pdf
    NIST TN 1297, section 5:
    https://www.nist.gov/pml/nist-technical-note-1297/nist-tn-1297-5-combined-standard-uncertainty
    """

    engine_id: str
    run_id: str
    condition_id: str
    parameter_names: tuple[str, ...]
    values: tuple[float, ...]
    measurement_covariance: tuple[tuple[float, ...], ...]
    source: str

    def __post_init__(self) -> None:
        if not self.engine_id.strip() or not self.run_id.strip():
            raise DomainError("Hot-fire engine and run identifiers must not be empty.")
        if not self.condition_id.strip() or not self.source.strip():
            raise DomainError("Hot-fire condition and measurement source are required.")
        names = tuple(str(name).strip() for name in self.parameter_names)
        try:
            values = tuple(float(value) for value in self.values)
        except (TypeError, ValueError) as error:
            raise DomainError("Hot-fire vector values must be numeric.") from error
        size = len(names)
        if size < 2 or len(set(names)) != size:
            raise DomainError("Vector calibration needs at least two unique parameters.")
        if any(not name for name in names):
            raise DomainError("Hot-fire vector parameter names must not be empty.")
        if len(values) != size or any(not isfinite(value) or value <= 0.0 for value in values):
            raise DomainError("Hot-fire vector values must be positive finite scales.")
        covariance = _matrix(
            self.measurement_covariance,
            size=size,
            field="measurement_covariance",
        )
        for row in range(size):
            if covariance[row][row] < 0.0:
                raise DomainError("Measurement covariance diagonal must be nonnegative.")
            for column in range(row):
                if abs(covariance[row][column] - covariance[column][row]) > 1e-12:
                    raise DomainError("Measurement covariance must be symmetric.")
        if not _is_positive_semidefinite(covariance):
            raise DomainError("Measurement covariance must be positive semidefinite.")
        object.__setattr__(self, "parameter_names", names)
        object.__setattr__(self, "values", values)
        object.__setattr__(self, "measurement_covariance", covariance)


@dataclass(frozen=True, slots=True)
class CorrelationRegularizationEvidence:
    """Recorded transformation from a raw estimate to a valid copula matrix.

    Pairwise coefficients are first clipped to the open physical interval and,
    only when necessary, off-diagonal terms are uniformly shrunk toward the
    identity until strict positive definiteness is obtained.  The transformation
    is evidence, not hidden numerical cleanup.

    References
    ----------
    JCGM 100:2008, section C.3.6:
    https://www.bipm.org/documents/20126/2071204/JCGM_100_2008_E.pdf
    NASA-STD-7009B: https://standards.nasa.gov/standard/NASA/NASA-STD-7009
    """

    scope: str
    raw_matrix: tuple[tuple[float, ...], ...]
    effective_matrix: tuple[tuple[float, ...], ...]
    coefficient_clip_count: int
    identity_shrinkage: float
    iterations: int
    regularized: bool


def _regularized_correlation(
    covariance: tuple[tuple[float, ...], ...],
    *,
    scope: str,
) -> CorrelationRegularizationEvidence | None:
    size = len(covariance)
    variances = tuple(covariance[index][index] for index in range(size))
    if any(value <= 0.0 for value in variances):
        return None
    raw_rows = [[0.0] * size for _ in range(size)]
    clipped_rows = [[0.0] * size for _ in range(size)]
    clip_count = 0
    limit = 1.0 - 1e-12
    for row in range(size):
        raw_rows[row][row] = 1.0
        clipped_rows[row][row] = 1.0
        for column in range(row):
            coefficient = covariance[row][column] / sqrt(variances[row] * variances[column])
            raw_rows[row][column] = coefficient
            raw_rows[column][row] = coefficient
            clipped = max(-limit, min(limit, coefficient))
            if clipped != coefficient:
                clip_count += 1
            clipped_rows[row][column] = clipped
            clipped_rows[column][row] = clipped
    raw = tuple(tuple(row) for row in raw_rows)
    clipped = tuple(tuple(row) for row in clipped_rows)
    if _is_positive_definite(clipped):
        return CorrelationRegularizationEvidence(
            scope=scope,
            raw_matrix=raw,
            effective_matrix=clipped,
            coefficient_clip_count=clip_count,
            identity_shrinkage=0.0,
            iterations=0,
            regularized=clip_count > 0,
        )

    low = 0.0
    high = 1.0
    iterations = 0
    for _ in range(80):
        shrinkage = 0.5 * (low + high)
        candidate = tuple(
            tuple(
                1.0 if row == column else (1.0 - shrinkage) * clipped[row][column]
                for column in range(size)
            )
            for row in range(size)
        )
        iterations += 1
        if _is_positive_definite(candidate):
            high = shrinkage
        else:
            low = shrinkage
    effective = tuple(
        tuple(
            1.0 if row == column else (1.0 - high) * clipped[row][column] for column in range(size)
        )
        for row in range(size)
    )
    if not _is_positive_definite(effective):
        raise DomainError("Correlation regularization did not reach positive definiteness.")
    return CorrelationRegularizationEvidence(
        scope=scope,
        raw_matrix=raw,
        effective_matrix=effective,
        coefficient_clip_count=clip_count,
        identity_shrinkage=high,
        iterations=iterations,
        regularized=True,
    )


@dataclass(frozen=True, slots=True)
class HotFireMultivariateCalibration:
    """Traceable paired calibration and ensemble-ready correlation model.

    Marginal variances use the scalar Mandel-Paule estimator.  Cross-covariance
    uses a paired method-of-moments estimate after subtracting the declared
    measurement covariance.  This hybrid keeps the reviewed nonnegative
    marginal estimator while preserving same-run cross-parameter evidence.

    References
    ----------
    JCGM 100:2008, sections 5.2 and C.3.6:
    https://www.bipm.org/documents/20126/2071204/JCGM_100_2008_E.pdf
    NIST/SEMATECH variance components:
    https://www.itl.nist.gov/div898/handbook/prc/section4/prc44.htm
    NASA-STD-7009B: https://standards.nasa.gov/standard/NASA/NASA-STD-7009
    """

    parameter_names: tuple[str, ...]
    condition_id: str
    method: str
    source: str
    rationale: str
    measurements: tuple[HotFireVectorMeasurement, ...]
    marginals: tuple[HotFireScaleCalibration, ...]
    within_engine_process_covariance: tuple[tuple[float, ...], ...]
    between_engine_process_covariance: tuple[tuple[float, ...], ...]
    population_process_covariance: tuple[tuple[float, ...], ...]
    population_correlation: CorrelationRegularizationEvidence
    same_engine_correlation: CorrelationRegularizationEvidence | None
    relative_tolerance: float
    maximum_iterations: int
    reference_ids: tuple[str, ...]
    warnings: tuple[str, ...]

    @property
    def calibration_hash(self) -> str:
        """Return a SHA-256 over inputs, estimates, solver evidence, and sources.

        References
        ----------
        NASA-STD-7009B: https://standards.nasa.gov/standard/NASA/NASA-STD-7009
        """

        return hashlib.sha256(canonical_json_bytes(self)).hexdigest()

    def residual_model(
        self,
        *,
        scope: CalibrationScope = CalibrationScope.POPULATION,
        coverage_factor: float = 3.0,
    ) -> tuple[tuple[UncertainParameter, ...], CorrelationModel]:
        """Return calibrated marginal residuals and their Gaussian copula.

        Apply every marginal ``consensus_scale`` to its deterministic base
        quantity before sampling these residuals.  Correlation coefficients are
        invariant under the positive consensus rescaling.

        References
        ----------
        JCGM 100:2008, section C.3.6:
        https://www.bipm.org/documents/20126/2071204/JCGM_100_2008_E.pdf
        NASA-SP-2011-3421: https://ntrs.nasa.gov/citations/20120001369
        """

        if not isinstance(scope, CalibrationScope):
            raise DomainError("Calibration scope must be a CalibrationScope value.")
        evidence = (
            self.population_correlation
            if scope is CalibrationScope.POPULATION
            else self.same_engine_correlation
        )
        if evidence is None:
            raise DomainError(
                "Same-engine correlation is undefined because at least one within-engine "
                "residual variance is zero."
            )
        parameters = tuple(
            marginal.residual_uncertain_parameter(
                scope=scope,
                coverage_factor=coverage_factor,
            )
            for marginal in self.marginals
        )
        correlation = CorrelationModel(
            parameter_names=self.parameter_names,
            matrix=evidence.effective_matrix,
            source=(f"{self.source}; paired hot-fire calibration sha256={self.calibration_hash}"),
            rationale=(
                f"{scope.value} paired residual correlation after measurement-covariance "
                f"removal; {self.rationale}"
            ),
        )
        return parameters, correlation

    def to_dict(self, *, include_hash: bool = True) -> dict[str, Any]:
        """Return the complete versioned calibration document.

        References
        ----------
        NASA-STD-7009B: https://standards.nasa.gov/standard/NASA/NASA-STD-7009
        """

        payload = json.loads(canonical_json_bytes(self).decode("utf-8"))
        payload = {"schema": HOT_FIRE_MULTIVARIATE_SCHEMA, **payload}
        if include_hash:
            payload["calibration_hash"] = self.calibration_hash
        return payload

    def to_json(self, *, indent: int | None = 2) -> str:
        """Serialize the versioned calibration with deterministic field ordering.

        References
        ----------
        NASA-STD-7009B: https://standards.nasa.gov/standard/NASA/NASA-STD-7009
        """

        return json.dumps(
            self.to_dict(),
            ensure_ascii=False,
            allow_nan=False,
            sort_keys=True,
            indent=indent,
        )

    @classmethod
    def from_dict(cls, payload: dict[str, Any]) -> HotFireMultivariateCalibration:
        """Recompute and integrity-check a persisted calibration document.

        References
        ----------
        NASA-STD-7009B: https://standards.nasa.gov/standard/NASA/NASA-STD-7009
        """

        if not isinstance(payload, dict):
            raise InputError("Hot-fire multivariate calibration root must be an object.")
        if payload.get("schema") != HOT_FIRE_MULTIVARIATE_SCHEMA:
            raise InputError("Unsupported hot-fire multivariate calibration schema.")
        declared_hash = payload.get("calibration_hash")
        try:
            measurements = tuple(
                HotFireVectorMeasurement(
                    engine_id=str(item["engine_id"]),
                    run_id=str(item["run_id"]),
                    condition_id=str(item["condition_id"]),
                    parameter_names=tuple(str(name) for name in item["parameter_names"]),
                    values=tuple(float(value) for value in item["values"]),
                    measurement_covariance=tuple(
                        tuple(float(value) for value in row)
                        for row in item["measurement_covariance"]
                    ),
                    source=str(item["source"]),
                )
                for item in payload["measurements"]
            )
            rebuilt = calibrate_hot_fire_multivariate(
                measurements,
                source=str(payload["source"]),
                rationale=str(payload["rationale"]),
                relative_tolerance=float(payload["relative_tolerance"]),
                maximum_iterations=int(payload["maximum_iterations"]),
            )
        except (KeyError, TypeError, ValueError, DomainError) as error:
            raise InputError("Malformed hot-fire multivariate calibration.") from error
        if declared_hash != rebuilt.calibration_hash:
            raise InputError("Hot-fire multivariate calibration hash does not reproduce.")
        if canonical_json_bytes(payload) != canonical_json_bytes(rebuilt.to_dict()):
            raise InputError("Hot-fire multivariate calibration semantics do not reproduce.")
        return rebuilt

    @classmethod
    def from_json(cls, text: str) -> HotFireMultivariateCalibration:
        """Parse and integrity-check a calibration JSON document.

        References
        ----------
        NASA-STD-7009B: https://standards.nasa.gov/standard/NASA/NASA-STD-7009
        """

        try:
            payload = json.loads(text)
        except (TypeError, json.JSONDecodeError) as error:
            raise InputError("Hot-fire multivariate calibration is not valid JSON.") from error
        return cls.from_dict(payload)


def _group_vectors(
    measurements: tuple[HotFireVectorMeasurement, ...],
) -> tuple[tuple[str, tuple[HotFireVectorMeasurement, ...]], ...]:
    engine_ids = tuple(dict.fromkeys(item.engine_id for item in measurements))
    return tuple(
        (engine_id, tuple(item for item in measurements if item.engine_id == engine_id))
        for engine_id in engine_ids
    )


def _covariance_sum(
    items: tuple[HotFireVectorMeasurement, ...],
    row: int,
    column: int,
) -> float:
    return sum(item.measurement_covariance[row][column] for item in items)


def calibrate_hot_fire_multivariate(
    measurements: tuple[HotFireVectorMeasurement, ...],
    *,
    source: str,
    rationale: str,
    relative_tolerance: float = 1e-10,
    maximum_iterations: int = 200,
) -> HotFireMultivariateCalibration:
    """Fit marginal random effects and paired residual correlations.

    The within-engine cross-covariance is the pooled paired covariance minus
    the expected declared measurement covariance.  The between-engine
    cross-covariance is the sample covariance of engine means minus the average
    covariance of an engine mean.  Marginal diagonals come from the reviewed
    scalar Mandel-Paule fits; off-diagonal regularization, if required for a
    Gaussian copula, is fully retained as evidence.

    References
    ----------
    JCGM 100:2008, sections 5.2 and C.3.6:
    https://www.bipm.org/documents/20126/2071204/JCGM_100_2008_E.pdf
    NIST/SEMATECH variance components:
    https://www.itl.nist.gov/div898/handbook/prc/section4/prc44.htm
    NIST Dataplot consensus mean:
    https://www.itl.nist.gov/div898/software/dataplot/refman1/auxillar/consmean.htm
    NASA-STD-7009B: https://standards.nasa.gov/standard/NASA/NASA-STD-7009
    """

    if not source.strip() or not rationale.strip():
        raise DomainError("Multivariate calibration source and rationale are required.")
    if len(measurements) < 3:
        raise DomainError("Multivariate calibration requires at least three observations.")
    parameter_names = measurements[0].parameter_names
    if any(item.parameter_names != parameter_names for item in measurements):
        raise DomainError("All paired measurements must use identical parameter order.")
    identities = tuple((item.engine_id, item.run_id) for item in measurements)
    if len(identities) != len(set(identities)):
        raise DomainError("Paired hot-fire engine/run identities must be unique.")
    conditions = {item.condition_id for item in measurements}
    if len(conditions) != 1:
        raise DomainError("Paired hot-fire measurements cannot mix operating conditions.")
    groups = _group_vectors(measurements)
    if len(groups) < 2 or len(measurements) - len(groups) <= 0:
        raise DomainError(
            "Paired calibration needs at least two engines and one repeated observation."
        )

    marginals: list[HotFireScaleCalibration] = []
    for index, parameter_name in enumerate(parameter_names):
        marginal_measurements = tuple(
            HotFireScaleMeasurement(
                engine_id=item.engine_id,
                run_id=item.run_id,
                condition_id=item.condition_id,
                value=item.values[index],
                standard_uncertainty=sqrt(item.measurement_covariance[index][index]),
                source=item.source,
            )
            for item in measurements
        )
        marginals.append(
            calibrate_hot_fire_scale(
                parameter_name,
                marginal_measurements,
                source=source,
                rationale=f"{rationale}; paired marginal {parameter_name}",
                relative_tolerance=relative_tolerance,
                maximum_iterations=maximum_iterations,
            )
        )

    size = len(parameter_names)
    within_df = len(measurements) - len(groups)
    within = [[0.0] * size for _ in range(size)]
    between = [[0.0] * size for _ in range(size)]
    for index, marginal in enumerate(marginals):
        within[index][index] = marginal.within_engine_process_variance
        between[index][index] = marginal.between_engine_process_variance

    engine_means: list[tuple[float, ...]] = []
    for _engine_id, items in groups:
        engine_means.append(
            tuple(fmean(item.values[index] for item in items) for index in range(size))
        )

    for row in range(size):
        for column in range(row):
            observed_within_sum = 0.0
            measurement_within_sum = 0.0
            for _engine_id, items in groups:
                count = len(items)
                row_mean = fmean(item.values[row] for item in items)
                column_mean = fmean(item.values[column] for item in items)
                observed_within_sum += sum(
                    (item.values[row] - row_mean) * (item.values[column] - column_mean)
                    for item in items
                )
                measurement_within_sum += (1.0 - 1.0 / count) * _covariance_sum(items, row, column)
            within_covariance = (observed_within_sum - measurement_within_sum) / within_df
            within[row][column] = within_covariance
            within[column][row] = within_covariance

            row_mean_of_means = fmean(item[row] for item in engine_means)
            column_mean_of_means = fmean(item[column] for item in engine_means)
            sample_engine_covariance = sum(
                (item[row] - row_mean_of_means) * (item[column] - column_mean_of_means)
                for item in engine_means
            ) / (len(groups) - 1)
            mean_noise_covariances = []
            for _engine_id, items in groups:
                count = len(items)
                mean_noise_covariances.append(
                    within_covariance / count + _covariance_sum(items, row, column) / count**2
                )
            between_covariance = sample_engine_covariance - fmean(mean_noise_covariances)
            between[row][column] = between_covariance
            between[column][row] = between_covariance

    within_matrix = tuple(tuple(row) for row in within)
    between_matrix = tuple(tuple(row) for row in between)
    population_matrix = tuple(
        tuple(within[row][column] + between[row][column] for column in range(size))
        for row in range(size)
    )
    population_correlation = _regularized_correlation(
        population_matrix,
        scope=CalibrationScope.POPULATION.value,
    )
    if population_correlation is None:
        raise DomainError("Every paired parameter requires positive population residual variance.")
    same_engine_correlation = _regularized_correlation(
        within_matrix,
        scope=CalibrationScope.SAME_ENGINE.value,
    )
    warnings = [
        "Cross-covariance is a paired method-of-moments estimate; marginal variances use Mandel-Paule.",
        "Run-to-run measurement errors are assumed independent after applying each declared covariance.",
        "Sampling uncertainty, model-form discrepancy, drift, and qualification-to-flight transfer are excluded.",
        "Apply every marginal consensus correction before sampling residual scales.",
    ]
    if population_correlation.regularized:
        warnings.append(
            "Population correlation required recorded clipping or shrinkage to become positive definite."
        )
    if same_engine_correlation is None:
        warnings.append(
            "Same-engine correlation is unavailable because at least one within-engine variance is zero."
        )
    elif same_engine_correlation.regularized:
        warnings.append(
            "Same-engine correlation required recorded clipping or shrinkage to become positive definite."
        )
    return HotFireMultivariateCalibration(
        parameter_names=parameter_names,
        condition_id=next(iter(conditions)),
        method="paired-random-effects-covariance-with-mandel-paule-marginals",
        source=source,
        rationale=rationale,
        measurements=measurements,
        marginals=tuple(marginals),
        within_engine_process_covariance=within_matrix,
        between_engine_process_covariance=between_matrix,
        population_process_covariance=population_matrix,
        population_correlation=population_correlation,
        same_engine_correlation=same_engine_correlation,
        relative_tolerance=relative_tolerance,
        maximum_iterations=maximum_iterations,
        reference_ids=HOT_FIRE_MULTIVARIATE_REFERENCE_IDS,
        warnings=tuple(warnings),
    )
