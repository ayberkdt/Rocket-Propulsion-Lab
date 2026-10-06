"""Traceable uncertainty propagation for pressure-fed propulsion histories.

The module keeps uncertainty definition separate from astrodynamics. Each
realization is a physically closed propulsion artifact that a future Sidera
ensemble adapter can propagate without resampling or altering engine physics.

References
----------
NASA-STD-7009B: https://standards.nasa.gov/standard/NASA/NASA-STD-7009
NASA-SP-2011-3421: https://ntrs.nasa.gov/citations/20120001369
NASA-SP-2009-569: https://ntrs.nasa.gov/citations/20090023159
"""

from __future__ import annotations

from dataclasses import dataclass, replace
from enum import Enum
from math import exp, isfinite, log, sqrt
from random import Random
from statistics import NormalDist, fmean

from rocket_propulsion.core.errors import DomainError

from .blowdown import (
    BlowdownTankModel,
    PressurePerformanceLaw,
    _scaled_operating_point,
    simulate_pressure_fed_blowdown,
)
from .models import EngineOperatingPoint
from .tabulated import TabulatedBurnArtifact

UNCERTAINTY_REFERENCE_IDS = (
    "NASA-STD-7009B",
    "NASA-SP-2011-3421",
    "NASA-SP-2009-569",
    "NASA-SP-8112",
    "NASA-CR-131400",
)


class UncertaintyClass(str, Enum):
    """Source class assigned to every uncertain input.

    ``ALEATORY`` denotes modeled variability between realizations;
    ``EPISTEMIC`` denotes incomplete knowledge that may be reduced with data.

    References
    ----------
    NASA-STD-7009B: https://standards.nasa.gov/standard/NASA/NASA-STD-7009
    NASA-SP-2009-569: https://ntrs.nasa.gov/citations/20090023159
    """

    ALEATORY = "aleatory"
    EPISTEMIC = "epistemic"


class DistributionKind(str, Enum):
    """Supported bounded marginal distributions.

    References
    ----------
    NASA-SP-2011-3421: https://ntrs.nasa.gov/citations/20120001369
    """

    UNIFORM = "uniform"
    TRIANGULAR = "triangular"
    TRUNCATED_NORMAL = "truncated-normal"
    LOG_UNIFORM = "log-uniform"


@dataclass(frozen=True, slots=True)
class UncertainParameter:
    """One bounded, sourced propulsion uncertainty input.

    ``nominal`` is also the triangular mode. ``standard_deviation`` is required
    only for a truncated-normal marginal. The source and rationale are
    mandatory so an arbitrary engineering percentage cannot enter unnoticed.

    References
    ----------
    NASA-STD-7009B: https://standards.nasa.gov/standard/NASA/NASA-STD-7009
    NASA-SP-2011-3421: https://ntrs.nasa.gov/citations/20120001369
    """

    name: str
    nominal: float
    lower: float
    upper: float
    distribution: DistributionKind
    uncertainty_class: UncertaintyClass
    source: str
    rationale: str
    standard_deviation: float | None = None

    def __post_init__(self) -> None:
        if not self.name.strip():
            raise DomainError("Uncertain parameter name must not be empty.")
        if not all(isfinite(value) for value in (self.nominal, self.lower, self.upper)):
            raise DomainError("Uncertain parameter values must be finite.", field=self.name)
        if not self.lower < self.upper:
            raise DomainError("Uncertain parameter lower bound must be below upper bound.")
        if not self.lower <= self.nominal <= self.upper:
            raise DomainError("Uncertain parameter nominal must lie within its bounds.")
        if not self.source.strip() or not self.rationale.strip():
            raise DomainError("Uncertain parameter source and rationale are required.")
        if self.distribution is DistributionKind.TRUNCATED_NORMAL:
            if (
                self.standard_deviation is None
                or not isfinite(self.standard_deviation)
                or self.standard_deviation <= 0.0
            ):
                raise DomainError("Truncated-normal uncertainty requires positive sigma.")
        elif self.standard_deviation is not None:
            raise DomainError("Sigma is only valid for truncated-normal uncertainty.")
        if self.distribution is DistributionKind.LOG_UNIFORM and self.lower <= 0.0:
            raise DomainError("Log-uniform uncertainty requires positive bounds.")

    def value_at_probability(self, probability: float) -> float:
        """Transform a unit-interval probability into this marginal.

        References
        ----------
        NASA-SP-2011-3421: https://ntrs.nasa.gov/citations/20120001369
        """

        if not isfinite(probability) or not 0.0 < probability < 1.0:
            raise DomainError("Sampling probability must be finite and in (0, 1).")
        span = self.upper - self.lower
        if self.distribution is DistributionKind.UNIFORM:
            return self.lower + probability * span
        if self.distribution is DistributionKind.LOG_UNIFORM:
            return exp(log(self.lower) + probability * (log(self.upper) - log(self.lower)))
        if self.distribution is DistributionKind.TRIANGULAR:
            mode_fraction = (self.nominal - self.lower) / span
            if probability <= mode_fraction:
                return self.lower + sqrt(probability * span * (self.nominal - self.lower))
            return self.upper - sqrt(
                (1.0 - probability) * span * (self.upper - self.nominal)
            )
        assert self.standard_deviation is not None
        normal = NormalDist(mu=self.nominal, sigma=self.standard_deviation)
        lower_probability = normal.cdf(self.lower)
        upper_probability = normal.cdf(self.upper)
        transformed = lower_probability + probability * (
            upper_probability - lower_probability
        )
        return normal.inv_cdf(transformed)


@dataclass(frozen=True, slots=True)
class CorrelationModel:
    """Gaussian-copula correlation model for ordered uncertain inputs.

    The matrix must be symmetric, positive definite, and have a unit diagonal.
    Its source is mandatory because guessed correlations can dominate ensemble
    tails even when every marginal is individually reasonable.

    References
    ----------
    NASA-STD-7009B: https://standards.nasa.gov/standard/NASA/NASA-STD-7009
    NASA-SP-2011-3421: https://ntrs.nasa.gov/citations/20120001369
    """

    parameter_names: tuple[str, ...]
    matrix: tuple[tuple[float, ...], ...]
    source: str
    rationale: str

    def __post_init__(self) -> None:
        size = len(self.parameter_names)
        if size < 2 or len(set(self.parameter_names)) != size:
            raise DomainError("Correlation model needs at least two unique parameters.")
        if len(self.matrix) != size or any(len(row) != size for row in self.matrix):
            raise DomainError("Correlation matrix dimensions must match parameter names.")
        if not self.source.strip() or not self.rationale.strip():
            raise DomainError("Correlation source and rationale are required.")
        for row in self.matrix:
            if any(not isfinite(value) or abs(value) > 1.0 for value in row):
                raise DomainError("Correlation coefficients must be finite and in [-1, 1].")
        for index in range(size):
            if abs(self.matrix[index][index] - 1.0) > 1e-12:
                raise DomainError("Correlation matrix diagonal must equal one.")
            for other in range(index):
                if abs(self.matrix[index][other] - self.matrix[other][index]) > 1e-12:
                    raise DomainError("Correlation matrix must be symmetric.")
        _cholesky(self.matrix)


@dataclass(frozen=True, slots=True)
class UncertaintyDraw:
    """One deterministic, indexed uncertainty realization.

    References
    ----------
    NASA-SP-2011-3421: https://ntrs.nasa.gov/citations/20120001369
    """

    index: int
    values: tuple[tuple[str, float], ...]

    def as_dict(self) -> dict[str, float]:
        """Return this immutable draw as a fresh name/value mapping.

        References
        ----------
        NASA-SP-2011-3421: https://ntrs.nasa.gov/citations/20120001369
        """

        return dict(self.values)


@dataclass(frozen=True, slots=True)
class PropulsionMetricStatistics:
    """Empirical statistics and a mean standard-error interval for one metric.

    The 95% interval is a normal-approximation interval for the sample mean; it
    is not a tolerance interval for individual engine realizations.

    References
    ----------
    NASA-SP-2011-3421: https://ntrs.nasa.gov/citations/20120001369
    """

    metric: str
    mean: float
    standard_deviation: float
    mean_standard_error: float
    mean_ci95_lower: float
    mean_ci95_upper: float
    minimum: float
    p05: float
    p50: float
    p95: float
    maximum: float


@dataclass(frozen=True, slots=True)
class PropulsionSensitivity:
    """Pearson screening sensitivity for one parameter/metric pair.

    Correlation is a screening statistic, not proof of causality or a complete
    global-sensitivity decomposition.

    References
    ----------
    NASA-SP-2011-3421: https://ntrs.nasa.gov/citations/20120001369
    """

    metric: str
    parameter: str
    correlation: float
    absolute_correlation: float
    rank: int


@dataclass(frozen=True, slots=True)
class SamplingConvergenceEvidence:
    """First-half versus full-ensemble stability diagnostic.

    This diagnostic does not certify probabilistic convergence. It makes sample
    adequacy visible and allows a workflow to demand more realizations.

    References
    ----------
    NASA-STD-7009B: https://standards.nasa.gov/standard/NASA/NASA-STD-7009
    NASA-SP-2011-3421: https://ntrs.nasa.gov/citations/20120001369
    """

    metric: str
    half_sample_count: int
    full_sample_count: int
    mean_relative_change: float
    standard_deviation_relative_change: float
    requested_tolerance: float
    converged: bool


@dataclass(frozen=True, slots=True)
class BlowdownUncertaintySample:
    """One physically closed blowdown realization and its exchange artifact.

    References
    ----------
    NASA-SP-8112: https://ntrs.nasa.gov/citations/19760015212
    NASA-CR-131400: https://ntrs.nasa.gov/citations/19730012094
    """

    index: int
    inputs: tuple[tuple[str, float], ...]
    duration_s: float
    delivered_impulse_n_s: float
    consumed_propellant_kg: float
    thrust_centroid_time_s: float
    final_supply_pressure_pa: float
    cutoff_reason: str
    artifact: TabulatedBurnArtifact


@dataclass(frozen=True, slots=True)
class BlowdownUncertaintyResult:
    """Traceable propulsion ensemble ready for later trajectory propagation.

    References
    ----------
    NASA-STD-7009B: https://standards.nasa.gov/standard/NASA/NASA-STD-7009
    NASA-SP-2011-3421: https://ntrs.nasa.gov/citations/20120001369
    NASA-SP-8112: https://ntrs.nasa.gov/citations/19760015212
    """

    method: str
    seed: int
    sample_count: int
    aleatory_parameter_count: int
    epistemic_parameter_count: int
    parameters: tuple[UncertainParameter, ...]
    correlation: CorrelationModel | None
    samples: tuple[BlowdownUncertaintySample, ...]
    statistics: tuple[PropulsionMetricStatistics, ...]
    sensitivity: tuple[PropulsionSensitivity, ...]
    convergence: SamplingConvergenceEvidence
    reference_ids: tuple[str, ...]
    warnings: tuple[str, ...]


def _cholesky(matrix: tuple[tuple[float, ...], ...]) -> tuple[tuple[float, ...], ...]:
    size = len(matrix)
    lower = [[0.0 for _ in range(size)] for _ in range(size)]
    for row in range(size):
        for column in range(row + 1):
            residual = matrix[row][column] - sum(
                lower[row][index] * lower[column][index] for index in range(column)
            )
            if row == column:
                if residual <= 1e-12:
                    raise DomainError("Correlation matrix must be positive definite.")
                lower[row][column] = sqrt(residual)
            else:
                lower[row][column] = residual / lower[column][column]
    return tuple(tuple(row) for row in lower)


def _unit_samples(
    dimension: int, *, sample_count: int, seed: int, method: str
) -> list[list[float]]:
    random = Random(seed)
    if method == "latin-hypercube":
        columns: list[list[float]] = []
        for _ in range(dimension):
            column = [(index + random.random()) / sample_count for index in range(sample_count)]
            random.shuffle(column)
            columns.append(column)
        return [[columns[column][row] for column in range(dimension)] for row in range(sample_count)]
    if method == "monte-carlo":
        return [[random.random() for _ in range(dimension)] for _ in range(sample_count)]
    raise DomainError("Uncertainty method must be latin-hypercube or monte-carlo.")


def sample_uncertain_parameters(
    parameters: tuple[UncertainParameter, ...],
    *,
    sample_count: int,
    seed: int,
    method: str = "latin-hypercube",
    correlation: CorrelationModel | None = None,
) -> tuple[UncertaintyDraw, ...]:
    """Generate reproducible bounded draws, optionally through a Gaussian copula.

    Latin-hypercube stratification is exact for independent marginals. When a
    correlation model is supplied, the Gaussian-copula transform prioritizes
    the declared dependency structure; the final marginal ranks are therefore
    not claimed to remain an exact Latin hypercube.

    References
    ----------
    NASA-SP-2011-3421: https://ntrs.nasa.gov/citations/20120001369
    NASA-STD-7009B: https://standards.nasa.gov/standard/NASA/NASA-STD-7009
    """

    if not parameters:
        raise DomainError("At least one propulsion uncertainty parameter is required.")
    if len({parameter.name for parameter in parameters}) != len(parameters):
        raise DomainError("Uncertainty parameter names must be unique.")
    if not 4 <= sample_count <= 10_000:
        raise DomainError("Propulsion uncertainty sample count must be between 4 and 10000.")
    normalized_method = method.strip().lower()
    unit_rows = _unit_samples(
        len(parameters), sample_count=sample_count, seed=seed, method=normalized_method
    )
    if correlation is not None:
        expected_names = tuple(parameter.name for parameter in parameters)
        if correlation.parameter_names != expected_names:
            raise DomainError("Correlation parameter order must match uncertainty parameters.")
        lower = _cholesky(correlation.matrix)
        standard_normal = NormalDist()
        epsilon = 1e-14
        for row_index, row in enumerate(unit_rows):
            independent = [
                standard_normal.inv_cdf(min(1.0 - epsilon, max(epsilon, value)))
                for value in row
            ]
            correlated = [
                sum(lower[index][column] * independent[column] for column in range(index + 1))
                for index in range(len(parameters))
            ]
            unit_rows[row_index] = [standard_normal.cdf(value) for value in correlated]
    return tuple(
        UncertaintyDraw(
            index=index,
            values=tuple(
                (parameter.name, parameter.value_at_probability(row[column]))
                for column, parameter in enumerate(parameters)
            ),
        )
        for index, row in enumerate(unit_rows)
    )


def _percentile(values: list[float], probability: float) -> float:
    ordered = sorted(values)
    position = probability * (len(ordered) - 1)
    lower = int(position)
    upper = min(lower + 1, len(ordered) - 1)
    fraction = position - lower
    return ordered[lower] + fraction * (ordered[upper] - ordered[lower])


def _statistics(metric: str, values: list[float]) -> PropulsionMetricStatistics:
    mean = fmean(values)
    deviation = sqrt(sum((value - mean) ** 2 for value in values) / (len(values) - 1))
    standard_error = deviation / sqrt(len(values))
    return PropulsionMetricStatistics(
        metric=metric,
        mean=mean,
        standard_deviation=deviation,
        mean_standard_error=standard_error,
        mean_ci95_lower=mean - 1.96 * standard_error,
        mean_ci95_upper=mean + 1.96 * standard_error,
        minimum=min(values),
        p05=_percentile(values, 0.05),
        p50=_percentile(values, 0.50),
        p95=_percentile(values, 0.95),
        maximum=max(values),
    )


def _correlation(left: list[float], right: list[float]) -> float:
    left_mean = fmean(left)
    right_mean = fmean(right)
    numerator = sum(
        (x - left_mean) * (y - right_mean) for x, y in zip(left, right, strict=True)
    )
    denominator = sqrt(
        sum((value - left_mean) ** 2 for value in left)
        * sum((value - right_mean) ** 2 for value in right)
    )
    return 0.0 if denominator == 0.0 else numerator / denominator


def _relative_change(previous: float, current: float) -> float:
    return abs(current - previous) / max(abs(current), 1e-15)


def simulate_blowdown_uncertainty(
    tank: BlowdownTankModel,
    rated_point: EngineOperatingPoint,
    pressure_law: PressurePerformanceLaw,
    parameters: tuple[UncertainParameter, ...],
    *,
    reserve_mass_kg: float = 0.0,
    cant_efficiency: float = 1.0,
    sample_count: int = 128,
    seed: int = 42,
    method: str = "latin-hypercube",
    correlation: CorrelationModel | None = None,
    convergence_tolerance: float = 0.05,
    maximum_mass_step_kg: float = 0.1,
    integration_relative_tolerance: float = 1e-7,
) -> BlowdownUncertaintyResult:
    """Propagate sourced input uncertainty into closed blowdown histories.

    Supported input names are ``initial_pressure_pa``,
    ``initial_ullage_volume_m3``, ``initial_propellant_mass_kg``,
    ``propellant_density_kg_m3``, ``polytropic_exponent``,
    ``rated_delivered_thrust_n``, ``rated_system_specific_impulse_s``,
    ``minimum_supply_pressure_pa``, ``thrust_pressure_exponent``,
    ``specific_impulse_pressure_exponent``, ``chamber_pressure_exponent``,
    ``reserve_mass_kg`` and ``cant_efficiency``.

    Every realization recomputes stream mass flows from delivered thrust and
    system Isp, executes the deterministic blowdown convergence gate, and
    retains its complete tabulated artifact. Invalid realizations abort the
    study with their sample index; they are never silently discarded.

    References
    ----------
    NASA-STD-7009B: https://standards.nasa.gov/standard/NASA/NASA-STD-7009
    NASA-SP-2011-3421: https://ntrs.nasa.gov/citations/20120001369
    NASA-SP-8112: https://ntrs.nasa.gov/citations/19760015212
    NASA-CR-131400: https://ntrs.nasa.gov/citations/19730012094
    """

    if not isfinite(convergence_tolerance) or convergence_tolerance <= 0.0:
        raise DomainError("Uncertainty convergence tolerance must be positive.")
    nominal = {
        "initial_pressure_pa": tank.initial_pressure_pa,
        "initial_ullage_volume_m3": tank.initial_ullage_volume_m3,
        "initial_propellant_mass_kg": tank.initial_propellant_mass_kg,
        "propellant_density_kg_m3": tank.propellant_density_kg_m3,
        "polytropic_exponent": tank.polytropic_exponent,
        "rated_delivered_thrust_n": rated_point.delivered_thrust_n,
        "rated_system_specific_impulse_s": rated_point.system_specific_impulse_s,
        "minimum_supply_pressure_pa": pressure_law.minimum_supply_pressure_pa,
        "thrust_pressure_exponent": pressure_law.thrust_pressure_exponent,
        "specific_impulse_pressure_exponent": pressure_law.specific_impulse_pressure_exponent,
        "chamber_pressure_exponent": pressure_law.chamber_pressure_exponent,
        "reserve_mass_kg": reserve_mass_kg,
        "cant_efficiency": cant_efficiency,
    }
    positive_names = {
        "initial_pressure_pa",
        "initial_ullage_volume_m3",
        "initial_propellant_mass_kg",
        "propellant_density_kg_m3",
        "polytropic_exponent",
        "rated_delivered_thrust_n",
        "rated_system_specific_impulse_s",
        "minimum_supply_pressure_pa",
    }
    nonnegative_names = {
        "thrust_pressure_exponent",
        "chamber_pressure_exponent",
    }
    for parameter in parameters:
        if parameter.name not in nominal:
            raise DomainError(f"Unsupported blowdown uncertainty input: {parameter.name}.")
        scale = max(1.0, abs(nominal[parameter.name]))
        if abs(parameter.nominal - nominal[parameter.name]) > 1e-10 * scale:
            raise DomainError(
                f"Uncertainty nominal for {parameter.name} does not match the model input."
            )
        if parameter.name in positive_names and parameter.lower <= 0.0:
            raise DomainError(f"Uncertainty bounds for {parameter.name} must remain positive.")
        if parameter.name in nonnegative_names and parameter.lower < 0.0:
            raise DomainError(
                f"Uncertainty bounds for {parameter.name} must remain non-negative."
            )
        if parameter.name == "reserve_mass_kg" and parameter.lower < 0.0:
            raise DomainError("Reserve-mass uncertainty cannot include negative values.")
        if parameter.name == "cant_efficiency" and (
            parameter.lower <= 0.0 or parameter.upper > 1.0
        ):
            raise DomainError("Cant-efficiency uncertainty must remain in (0, 1].")
    parameter_map = {parameter.name: parameter for parameter in parameters}
    pressure_bounds = parameter_map.get("initial_pressure_pa")
    minimum_bounds = parameter_map.get("minimum_supply_pressure_pa")
    lowest_initial_pressure = (
        pressure_bounds.lower if pressure_bounds is not None else tank.initial_pressure_pa
    )
    highest_initial_pressure = (
        pressure_bounds.upper if pressure_bounds is not None else tank.initial_pressure_pa
    )
    highest_minimum_pressure = (
        minimum_bounds.upper
        if minimum_bounds is not None
        else pressure_law.minimum_supply_pressure_pa
    )
    if highest_initial_pressure > pressure_law.maximum_pressure_pa:
        raise DomainError("Initial-pressure uncertainty exceeds pressure-law validity.")
    if highest_minimum_pressure >= lowest_initial_pressure:
        raise DomainError("Minimum-pressure uncertainty can eliminate the usable pressure range.")
    lowest_inventory = (
        parameter_map["initial_propellant_mass_kg"].lower
        if "initial_propellant_mass_kg" in parameter_map
        else tank.initial_propellant_mass_kg
    )
    highest_reserve = (
        parameter_map["reserve_mass_kg"].upper
        if "reserve_mass_kg" in parameter_map
        else reserve_mass_kg
    )
    if highest_reserve >= lowest_inventory:
        raise DomainError("Reserve uncertainty can consume the complete propellant inventory.")

    draws = sample_uncertain_parameters(
        parameters,
        sample_count=sample_count,
        seed=seed,
        method=method,
        correlation=correlation,
    )
    realizations: list[BlowdownUncertaintySample] = []
    for draw in draws:
        values = nominal | draw.as_dict()
        sampled_tank = BlowdownTankModel(
            initial_pressure_pa=values["initial_pressure_pa"],
            initial_ullage_volume_m3=values["initial_ullage_volume_m3"],
            initial_propellant_mass_kg=values["initial_propellant_mass_kg"],
            propellant_density_kg_m3=values["propellant_density_kg_m3"],
            polytropic_exponent=values["polytropic_exponent"],
        )
        sampled_point = _scaled_operating_point(
            rated_point,
            supply_pressure_pa=pressure_law.reference_supply_pressure_pa,
            delivered_thrust_n=values["rated_delivered_thrust_n"],
            system_specific_impulse_s=values["rated_system_specific_impulse_s"],
            chamber_pressure_pa=rated_point.chamber_pressure_pa,
            model_id="uncertain_rated_point_v1",
            source=f"{rated_point.source}; uncertainty realization {draw.index}",
            warning="Rated performance contains a sourced uncertainty realization.",
        )
        sampled_law = PressurePerformanceLaw(
            reference_supply_pressure_pa=pressure_law.reference_supply_pressure_pa,
            minimum_supply_pressure_pa=values["minimum_supply_pressure_pa"],
            maximum_supply_pressure_pa=pressure_law.maximum_pressure_pa,
            thrust_pressure_exponent=values["thrust_pressure_exponent"],
            specific_impulse_pressure_exponent=values[
                "specific_impulse_pressure_exponent"
            ],
            chamber_pressure_exponent=values["chamber_pressure_exponent"],
        )
        try:
            result = simulate_pressure_fed_blowdown(
                sampled_tank,
                sampled_point,
                sampled_law,
                reserve_mass_kg=values["reserve_mass_kg"],
                cant_efficiency=values["cant_efficiency"],
                maximum_mass_step_kg=maximum_mass_step_kg,
                relative_tolerance=integration_relative_tolerance,
            )
        except DomainError as error:
            raise DomainError(
                f"Blowdown uncertainty realization {draw.index} is invalid: {error}"
            ) from error
        artifact = replace(
            result.artifact,
            reference_ids=tuple(
                dict.fromkeys(result.artifact.reference_ids + UNCERTAINTY_REFERENCE_IDS)
            ),
            warnings=result.artifact.warnings
            + ("This artifact is one sourced propulsion uncertainty realization.",),
        )
        realizations.append(
            BlowdownUncertaintySample(
                index=draw.index,
                inputs=draw.values,
                duration_s=artifact.duration_s,
                delivered_impulse_n_s=artifact.delivered_total_impulse_n_s,
                consumed_propellant_kg=artifact.consumed_propellant_kg,
                thrust_centroid_time_s=artifact.thrust_centroid_time_s,
                final_supply_pressure_pa=result.samples[-1].supply_pressure_pa,
                cutoff_reason=result.cutoff_reason,
                artifact=artifact,
            )
        )

    metric_values = {
        "duration_s": [sample.duration_s for sample in realizations],
        "delivered_impulse_n_s": [sample.delivered_impulse_n_s for sample in realizations],
        "consumed_propellant_kg": [sample.consumed_propellant_kg for sample in realizations],
        "thrust_centroid_time_s": [sample.thrust_centroid_time_s for sample in realizations],
        "final_supply_pressure_pa": [sample.final_supply_pressure_pa for sample in realizations],
    }
    statistics = tuple(
        _statistics(metric, values) for metric, values in metric_values.items()
    )
    sensitivity: list[PropulsionSensitivity] = []
    for metric, outputs in metric_values.items():
        ranked = sorted(
            (
                (
                    parameter.name,
                    _correlation(
                        [dict(sample.inputs)[parameter.name] for sample in realizations],
                        outputs,
                    ),
                )
                for parameter in parameters
            ),
            key=lambda item: abs(item[1]),
            reverse=True,
        )
        sensitivity.extend(
            PropulsionSensitivity(metric, name, value, abs(value), rank)
            for rank, (name, value) in enumerate(ranked, start=1)
        )

    primary = metric_values["delivered_impulse_n_s"]
    half_count = len(primary) // 2
    half_mean = fmean(primary[:half_count])
    full_mean = fmean(primary)
    half_deviation = sqrt(
        sum((value - half_mean) ** 2 for value in primary[:half_count])
        / (half_count - 1)
    )
    full_deviation = sqrt(
        sum((value - full_mean) ** 2 for value in primary) / (len(primary) - 1)
    )
    mean_change = _relative_change(half_mean, full_mean)
    deviation_change = _relative_change(half_deviation, full_deviation)
    convergence = SamplingConvergenceEvidence(
        metric="delivered_impulse_n_s",
        half_sample_count=half_count,
        full_sample_count=len(primary),
        mean_relative_change=mean_change,
        standard_deviation_relative_change=deviation_change,
        requested_tolerance=convergence_tolerance,
        converged=max(mean_change, deviation_change) <= convergence_tolerance,
    )
    warnings = [
        "Ensemble statistics are engineering uncertainty estimates, not flight certification.",
        "Pearson sensitivity is a screening measure and can miss nonlinear interactions.",
    ]
    if not convergence.converged:
        warnings.append(
            "Half/full ensemble stability did not meet the requested tolerance; increase sample count."
        )
    return BlowdownUncertaintyResult(
        method=method.strip().lower(),
        seed=seed,
        sample_count=sample_count,
        aleatory_parameter_count=sum(
            parameter.uncertainty_class is UncertaintyClass.ALEATORY
            for parameter in parameters
        ),
        epistemic_parameter_count=sum(
            parameter.uncertainty_class is UncertaintyClass.EPISTEMIC
            for parameter in parameters
        ),
        parameters=parameters,
        correlation=correlation,
        samples=tuple(realizations),
        statistics=statistics,
        sensitivity=tuple(sensitivity),
        convergence=convergence,
        reference_ids=UNCERTAINTY_REFERENCE_IDS,
        warnings=tuple(warnings),
    )
