"""Reproducible bounded studies with checked empirical arithmetic.

Uniform draws are exact affine transforms of the binary RNG fraction, rounded
once to binary64. LHS stratification belongs to latent draws; finite grids may
collapse distinct draws. Summary exactness is relative to retained binary data,
not measured physics, confidence coverage or a global sensitivity analysis.
"""

from __future__ import annotations

from collections.abc import Callable, Mapping
from dataclasses import dataclass
from fractions import Fraction
from random import Random

from rocket_propulsion.core.errors import DomainError
from rocket_propulsion.core.validation import integer, real, real_mapping, records, text
from rocket_propulsion.propulsion.burns.numeric_outputs import _finite_rational
from rocket_propulsion.propulsion.burns.probability_numerics import _sqrt_output
from rocket_propulsion.propulsion.burns.sample_statistics import (
    _SUMMARY_POLICY, _linear_quantile, _pearson, _sample_ledger,
)

_SAMPLING_POLICY = 'exact_binary_affine_uniform_latent_strata_v1'


class _FrozenValues(dict):
    """Owned scalar mapping compatible with dataclasses.asdict and JSON.

    Normal mapping mutations are disabled, including reinitialization. Like a
    frozen dataclass, this is an ownership contract, not a hostile-code sandbox.
    """
    __slots__ = ()

    def __new__(cls, values=()):
        instance = dict.__new__(cls)
        dict.update(instance, values)
        return instance

    def __init__(self, values=()):
        pass

    def _immutable(self, *args, **kwargs):
        raise TypeError('Study value mappings are immutable.')

    __setitem__ = __delitem__ = clear = pop = popitem = setdefault = update = __ior__ = _immutable

    def copy(self):
        return self

    def __deepcopy__(self, memo):
        return self

    def __reduce_ex__(self, protocol):
        return type(self), (dict(self),)


def _owned_values(value, label):
    if not isinstance(value, Mapping):
        raise DomainError(f'{label} must be a mapping.', field=label)
    return _FrozenValues(real_mapping(value, label))


@dataclass(frozen=True, slots=True)
class ParameterRange:
    name: str
    minimum: float
    maximum: float

    def __post_init__(self):
        object.__setattr__(self, 'name', text(self.name, 'Parameter name'))
        object.__setattr__(self, 'minimum', real(self.minimum, 'Parameter minimum'))
        object.__setattr__(self, 'maximum', real(self.maximum, 'Parameter maximum'))
        if self.minimum >= self.maximum:
            raise DomainError('Parameter minimum must be below maximum.', field=self.name)


@dataclass(frozen=True, slots=True)
class SweepPoint:
    inputs: Mapping[str, float]
    output: float

    def __post_init__(self):
        object.__setattr__(self, 'inputs', _owned_values(self.inputs, 'Sample inputs'))
        object.__setattr__(self, 'output', real(self.output, 'Evaluator output'))


@dataclass(frozen=True, slots=True)
class SensitivityRecord:
    parameter: str
    correlation: float | None
    absolute_correlation: float | None
    rank: int | None
    undefined_reason: str | None = None

    def __post_init__(self):
        object.__setattr__(self, 'parameter', text(self.parameter, 'Sensitivity parameter'))
        if self.correlation is None:
            if (self.absolute_correlation is not None or self.rank is not None
                    or self.undefined_reason != 'constant_input_or_output'):
                raise DomainError('Undefined Pearson requires null magnitude/rank and a reason.', field='Sensitivity')
            return
        coefficient = real(self.correlation, 'Pearson correlation')
        magnitude = real(self.absolute_correlation, 'Pearson magnitude')
        rank = integer(self.rank, 'Sensitivity rank')
        if (not -1 <= coefficient <= 1 or magnitude != abs(coefficient)
                or rank < 1 or self.undefined_reason is not None):
            raise DomainError('Inconsistent Pearson sensitivity record.', field='Sensitivity')
        object.__setattr__(self, 'correlation', coefficient)
        object.__setattr__(self, 'absolute_correlation', magnitude)
        object.__setattr__(self, 'rank', rank)


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
    parameters: tuple[ParameterRange, ...] = ()
    sampling_numerics: str = _SAMPLING_POLICY
    summary_numerics: str = _SUMMARY_POLICY
    warnings: tuple[str, ...] = ()

    def __post_init__(self):
        method = _method(self.method)
        count, seed = _controls(self.sample_count, self.seed)
        parameters = _validate_ranges(self.parameters)
        samples = records(self.samples, SweepPoint, 'Samples')
        sensitivity = records(self.sensitivity, SensitivityRecord, 'Sensitivity')
        mean = real(self.output_mean, 'Output mean')
        deviation = real(self.output_standard_deviation, 'Output standard deviation')
        percentiles = _owned_values(self.output_percentiles, 'Output percentiles')
        if self.sampling_numerics != _SAMPLING_POLICY or self.summary_numerics != _SUMMARY_POLICY:
            raise DomainError('Unknown study numerical policy.', field='Numerical policy')
        if len(samples) != count:
            raise DomainError('Retained sample count is inconsistent.', field='Samples')
        names = {parameter.name for parameter in parameters}
        for sample in samples:
            if set(sample.inputs) != names or any(
                not parameter.minimum <= sample.inputs[parameter.name] <= parameter.maximum
                for parameter in parameters
            ):
                raise DomainError('Retained inputs do not match declared parameter ranges.', field='Samples')
        expected_mean, expected_deviation, expected_percentiles = _summaries(samples)
        expected_sensitivity = _sensitivities(parameters, samples)
        if ((mean, deviation, percentiles) != (expected_mean, expected_deviation, expected_percentiles)
                or sensitivity != expected_sensitivity):
            raise DomainError('Retained summaries do not reproduce samples.', field='Study summaries')
        warnings = _warnings(method)
        if self.warnings and self.warnings != warnings:
            raise DomainError('Study warnings do not match the declared method.', field='Warnings')
        for name, value in dict(method=method, seed=seed, sample_count=count, parameters=parameters,
                                samples=samples, sensitivity=sensitivity, output_mean=mean,
                                output_standard_deviation=deviation, output_percentiles=percentiles,
                                warnings=warnings).items():
            object.__setattr__(self, name, value)


def _validate_ranges(parameters):
    parameters = records(parameters, ParameterRange, 'Parameter ranges')
    if not parameters or len({p.name for p in parameters}) != len(parameters):
        raise DomainError('Parameter ranges must be nonempty with unique names.', field='Parameter ranges')
    return parameters


def _count(value, maximum, label):
    value = integer(value, label)
    if not 2 <= value <= maximum:
        raise DomainError(f'{label} must be between 2 and {maximum}.', field=label)
    return value


def _controls(sample_count, seed):
    return _count(sample_count, 100000, 'Sample count'), integer(seed, 'Seed')


def _method(value):
    value = text(value, 'Uncertainty method').strip().lower()
    if value not in ('latin-hypercube', 'monte-carlo'):
        raise DomainError('Uncertainty method must be latin-hypercube or monte-carlo.', field='Uncertainty method')
    return value


def _affine(parameter, probability):
    # A convex combination cannot overflow its finite bounds. Round once; an
    # adjacent-to-zero grid cell may legitimately round a nonzero draw to zero.
    return float((1-probability)*Fraction(parameter.minimum)+probability*Fraction(parameter.maximum))


def latin_hypercube_samples(
    parameters: tuple[ParameterRange, ...], *, sample_count: int, seed: int
) -> tuple[Mapping[str, float], ...]:
    """Generate one random-jitter LHS, with exact latent half-open strata."""
    parameters = _validate_ranges(parameters)
    sample_count, seed = _controls(sample_count, seed)
    random = Random(seed)
    columns = {}
    for parameter in parameters:
        values = [_affine(parameter, (Fraction(index)+Fraction(random.random()))/sample_count)
                  for index in range(sample_count)]
        random.shuffle(values)
        columns[parameter.name] = values
    return tuple(_FrozenValues((p.name, columns[p.name][index]) for p in parameters)
                 for index in range(sample_count))


def monte_carlo_samples(
    parameters: tuple[ParameterRange, ...], *, sample_count: int, seed: int
) -> tuple[Mapping[str, float], ...]:
    """Generate seeded uniform draws over the declared finite bounds."""
    parameters = _validate_ranges(parameters)
    sample_count, seed = _controls(sample_count, seed)
    random = Random(seed)
    return tuple(_FrozenValues((p.name, _affine(p, Fraction(random.random()))) for p in parameters)
                 for _ in range(sample_count))


def _evaluator(value):
    if not callable(value):
        raise DomainError('Evaluator must be callable.', field='Evaluator')
    return value


def parameter_sweep(
    parameter: ParameterRange, *, point_count: int,
    evaluator: Callable[[Mapping[str, float]], float],
    fixed_inputs: Mapping[str, float] | None = None,
) -> tuple[SweepPoint, ...]:
    """Evaluate an inclusive sweep; the swept parameter overrides its fixed value."""
    parameter, = _validate_ranges((parameter,))
    point_count = _count(point_count, 10000, 'Point count')
    evaluator = _evaluator(evaluator)
    base = {} if fixed_inputs is None else _owned_values(fixed_inputs, 'Fixed inputs')
    points = []
    for index in range(point_count):
        inputs = _FrozenValues({**base, parameter.name: _affine(parameter, Fraction(index, point_count-1))})
        points.append(SweepPoint(inputs, real(evaluator(inputs), 'Evaluator output')))
    return tuple(points)


def _percentile(sorted_values, probability):
    return _linear_quantile(sorted_values, probability)


def _correlation(left, right):
    return _pearson(left, right)


def _summaries(samples):
    outputs = tuple(sample.output for sample in samples)
    ledger = _sample_ledger(outputs)
    return (_finite_rational(ledger.mean, 'Output mean'),
            _sqrt_output(ledger.variance, 'Output standard deviation'),
            _FrozenValues((key, _percentile(outputs, probability)) for key, probability in
                          (('p05',Fraction(1,20)), ('p50',Fraction(1,2)), ('p95',Fraction(19,20)))))


def _sensitivities(parameters, samples):
    outputs = tuple(sample.output for sample in samples)
    correlations = [(parameter.name, _correlation(tuple(point.inputs[parameter.name] for point in samples), outputs))
                    for parameter in parameters]
    # Stable ties follow declared parameter order. Undefined values follow all
    # defined values, including a defined zero; an undefined value has no rank.
    correlations.sort(key=lambda item: (item[1] is None, -abs(item[1]) if item[1] is not None else 0))
    return tuple(SensitivityRecord(name, value, abs(value) if value is not None else None,
                                   index if value is not None else None,
                                   'constant_input_or_output' if value is None else None)
                 for index, (name, value) in enumerate(correlations, start=1))


def _warnings(method):
    common = ('Empirical sample moments and percentiles do not establish confidence coverage or physical validation.',
              'Pearson is a linear screening coefficient, not a global sensitivity measure; constant vectors are undefined.')
    return common + (('LHS rows are dependent; IID standard-error formulas are not qualified for this design.',)
                     if method == 'latin-hypercube' else ())


def propagate_uncertainty(
    parameters: tuple[ParameterRange, ...], *, evaluator: Callable[[Mapping[str, float]], float],
    sample_count: int = 200, seed: int = 42, method: str = 'latin-hypercube',
) -> UncertaintyResult:
    """Evaluate owned samples and rank defined linear Pearson sensitivities."""
    parameters = _validate_ranges(parameters)
    sample_count, seed = _controls(sample_count, seed)
    method = _method(method)
    evaluator = _evaluator(evaluator)
    sampler = latin_hypercube_samples if method == 'latin-hypercube' else monte_carlo_samples
    inputs = sampler(parameters, sample_count=sample_count, seed=seed)
    samples = []
    for row in inputs:
        owned = _owned_values(row, 'Sample inputs')
        samples.append(SweepPoint(owned, real(evaluator(owned), 'Evaluator output')))
    samples = tuple(samples)
    mean, deviation, percentiles = _summaries(samples)
    return UncertaintyResult(method, seed, sample_count, samples, mean, deviation, percentiles,
                             _sensitivities(parameters, samples), parameters=parameters)
