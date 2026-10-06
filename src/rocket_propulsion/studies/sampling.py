"""Deterministic parameter sweeps and uncertainty propagation."""

from __future__ import annotations

from collections.abc import Callable, Mapping
from dataclasses import dataclass
from math import isfinite, sqrt
from random import Random
from statistics import fmean

from rocket_propulsion.core.errors import DomainError


@dataclass(frozen=True, slots=True)
class ParameterRange:
    name: str
    minimum: float
    maximum: float


@dataclass(frozen=True, slots=True)
class SweepPoint:
    inputs: Mapping[str, float]
    output: float


@dataclass(frozen=True, slots=True)
class SensitivityRecord:
    parameter: str
    correlation: float
    absolute_correlation: float
    rank: int


@dataclass(frozen=True, slots=True)
class UncertaintyResult:
    method: str
    seed: int
    sample_count: int
    samples: tuple[SweepPoint, ...]
    output_mean: float
    output_standard_deviation: float
    output_percentiles: Mapping[str, float]
    sensitivity: tuple[SensitivityRecord, ...]


def _validate_ranges(parameters: tuple[ParameterRange, ...]) -> None:
    if not parameters:
        raise DomainError("At least one parameter range is required.")
    names: set[str] = set()
    for parameter in parameters:
        if not parameter.name or parameter.name in names:
            raise DomainError("Parameter names must be non-empty and unique.")
        if not all(isfinite(value) for value in (parameter.minimum, parameter.maximum)):
            raise DomainError("Parameter bounds must be finite.", field=parameter.name)
        if parameter.minimum >= parameter.maximum:
            raise DomainError("Parameter minimum must be below maximum.", field=parameter.name)
        names.add(parameter.name)


def latin_hypercube_samples(
    parameters: tuple[ParameterRange, ...], *, sample_count: int, seed: int
) -> tuple[Mapping[str, float], ...]:
    """Generate a reproducible centered-jitter Latin hypercube."""

    _validate_ranges(parameters)
    if not 2 <= sample_count <= 100_000:
        raise DomainError("Sample count must be between 2 and 100000.")
    random = Random(seed)
    columns: dict[str, list[float]] = {}
    for parameter in parameters:
        values = [
            parameter.minimum
            + (index + random.random())
            / sample_count
            * (parameter.maximum - parameter.minimum)
            for index in range(sample_count)
        ]
        random.shuffle(values)
        columns[parameter.name] = values
    return tuple(
        {parameter.name: columns[parameter.name][index] for parameter in parameters}
        for index in range(sample_count)
    )


def monte_carlo_samples(
    parameters: tuple[ParameterRange, ...], *, sample_count: int, seed: int
) -> tuple[Mapping[str, float], ...]:
    """Generate reproducible independent uniform samples over declared bounds."""

    _validate_ranges(parameters)
    if not 2 <= sample_count <= 100_000:
        raise DomainError("Sample count must be between 2 and 100000.")
    random = Random(seed)
    return tuple(
        {
            parameter.name: random.uniform(parameter.minimum, parameter.maximum)
            for parameter in parameters
        }
        for _ in range(sample_count)
    )


def parameter_sweep(
    parameter: ParameterRange,
    *,
    point_count: int,
    evaluator: Callable[[Mapping[str, float]], float],
    fixed_inputs: Mapping[str, float] | None = None,
) -> tuple[SweepPoint, ...]:
    """Evaluate a deterministic inclusive one-parameter sweep."""

    _validate_ranges((parameter,))
    if not 2 <= point_count <= 10_000:
        raise DomainError("Point count must be between 2 and 10000.")
    base = dict(fixed_inputs or {})
    points: list[SweepPoint] = []
    for index in range(point_count):
        value = parameter.minimum + index / (point_count - 1) * (
            parameter.maximum - parameter.minimum
        )
        inputs = {**base, parameter.name: value}
        output = float(evaluator(inputs))
        if not isfinite(output):
            raise DomainError("Sweep evaluator returned a non-finite result.")
        points.append(SweepPoint(inputs, output))
    return tuple(points)


def _percentile(sorted_values: list[float], probability: float) -> float:
    position = probability * (len(sorted_values) - 1)
    lower = int(position)
    upper = min(lower + 1, len(sorted_values) - 1)
    fraction = position - lower
    return sorted_values[lower] * (1.0 - fraction) + sorted_values[upper] * fraction


def _correlation(left: list[float], right: list[float]) -> float:
    left_mean = fmean(left)
    right_mean = fmean(right)
    covariance = sum(
        (x - left_mean) * (y - right_mean) for x, y in zip(left, right)
    )
    left_sum = sum((value - left_mean) ** 2 for value in left)
    right_sum = sum((value - right_mean) ** 2 for value in right)
    denominator = sqrt(left_sum * right_sum)
    return 0.0 if denominator == 0.0 else covariance / denominator


def propagate_uncertainty(
    parameters: tuple[ParameterRange, ...],
    *,
    evaluator: Callable[[Mapping[str, float]], float],
    sample_count: int = 200,
    seed: int = 42,
    method: str = "latin-hypercube",
) -> UncertaintyResult:
    """Evaluate samples and rank linear input/output sensitivity."""

    normalized_method = method.strip().lower()
    if normalized_method == "latin-hypercube":
        inputs = latin_hypercube_samples(parameters, sample_count=sample_count, seed=seed)
    elif normalized_method == "monte-carlo":
        inputs = monte_carlo_samples(parameters, sample_count=sample_count, seed=seed)
    else:
        raise DomainError("Uncertainty method must be latin-hypercube or monte-carlo.")
    samples: list[SweepPoint] = []
    for sample in inputs:
        output = float(evaluator(sample))
        if not isfinite(output):
            raise DomainError("Uncertainty evaluator returned a non-finite result.")
        samples.append(SweepPoint(sample, output))
    outputs = [sample.output for sample in samples]
    mean = fmean(outputs)
    standard_deviation = sqrt(
        sum((value - mean) ** 2 for value in outputs) / (len(outputs) - 1)
    )
    correlations = [
        (
            parameter.name,
            _correlation(
                [float(sample.inputs[parameter.name]) for sample in samples], outputs
            ),
        )
        for parameter in parameters
    ]
    correlations.sort(key=lambda item: abs(item[1]), reverse=True)
    sensitivity = tuple(
        SensitivityRecord(name, value, abs(value), rank)
        for rank, (name, value) in enumerate(correlations, start=1)
    )
    ordered = sorted(outputs)
    return UncertaintyResult(
        method=normalized_method,
        seed=seed,
        sample_count=sample_count,
        samples=tuple(samples),
        output_mean=mean,
        output_standard_deviation=standard_deviation,
        output_percentiles={
            "p05": _percentile(ordered, 0.05),
            "p50": _percentile(ordered, 0.50),
            "p95": _percentile(ordered, 0.95),
        },
        sensitivity=sensitivity,
    )
