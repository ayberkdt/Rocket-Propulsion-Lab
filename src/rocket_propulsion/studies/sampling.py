"""Reproducible bounded parameter studies with checked empirical summaries.

Uniform draws are exact affine maps of the binary random fraction onto the declared
bounds, rounded once to binary64. Latin hypercube stratification [1]_ applies to
the latent draws; a finite binary64 grid may map distinct draws to the same value.
Summary exactness is relative to the retained binary samples. The summaries are
empirical screening statistics [2]_ [3]_ [4]_; they do not establish measured
physics, confidence coverage or a global sensitivity analysis.

References
----------
.. [1] MCKAY-1979. M. D. McKay, R. J. Beckman and W. J. Conover, "A Comparison of
       Three Methods for Selecting Values of Input Variables in the Analysis of
       Output from a Computer Code", Technometrics 21(2), 239-245, 1979.
.. [2] NIST TN 1297, "Guidelines for Evaluating and Expressing the Uncertainty of
       NIST Measurement Results", 1994 (Type A evaluation: sample mean and standard
       deviation).
.. [3] NUMPY-2.2-LINEAR-QUANTILE. NumPy 2.2, numpy.quantile: linear (Hyndman-Fan
       type 7) sample quantile definition.
.. [4] NIST-PEARSON. NIST Dataplot reference manual, "Correlation": centered-product
       (Pearson) correlation definition.
.. [5] NASA/SP-2011-3421, "Probabilistic Risk Assessment Procedures Guide for NASA
       Managers and Practitioners", NASA, 2011 (Monte Carlo and LHS uncertainty
       propagation).
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
    _SUMMARY_POLICY,
    _linear_quantile,
    _pearson,
    _sample_ledger,
)

_SAMPLING_POLICY = 'exact_binary_affine_uniform_latent_strata_v1'


class _FrozenValues(dict):
    """Owned scalar mapping compatible with ``dataclasses.asdict`` and JSON.

    Mutating mapping methods raise ``TypeError`` and reinitialization through
    ``__init__`` leaves the contents unchanged. Like a frozen dataclass, this is an
    ownership contract, not a sandbox against deliberately hostile code.
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
        """Return this mapping itself; it is immutable, so no copy is needed."""
        return self

    def __deepcopy__(self, memo):
        return self

    def __reduce_ex__(self, protocol):
        return type(self), (dict(self),)


def _owned_values(value, label):
    """Return an immutable copy of a mapping of names to finite reals.

    Parameters
    ----------
    value : Mapping
        Candidate name-to-value mapping.
    label : str
        Label used in errors.

    Returns
    -------
    _FrozenValues

    Raises
    ------
    DomainError
        If ``value`` is not a mapping or a value is not a finite real.
    """
    if not isinstance(value, Mapping):
        raise DomainError(f'{label} must be a mapping.', field=label)
    return _FrozenValues(real_mapping(value, label))


@dataclass(frozen=True, slots=True)
class ParameterRange:
    """Closed bounds of one uniformly sampled study parameter.

    Attributes
    ----------
    name : str
        Parameter name, used as the sample mapping key.
    minimum, maximum : float
        Finite bounds with ``minimum < maximum``.

    Raises
    ------
    DomainError
        If the name is empty or the bounds are not finite and ordered.

    References
    ----------
    .. [1] NASA/SP-2011-3421, "Probabilistic Risk Assessment Procedures Guide for NASA
           Managers and Practitioners", NASA, 2011 (Monte Carlo and LHS uncertainty
           propagation).
    """
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
    """One evaluated study point.

    Attributes
    ----------
    inputs : Mapping of str to float
        Immutable input values passed to the evaluator.
    output : float
        Finite evaluator output.

    Raises
    ------
    DomainError
        If an input or the output is not a finite real.
    """
    inputs: Mapping[str, float]
    output: float

    def __post_init__(self):
        object.__setattr__(self, 'inputs', _owned_values(self.inputs, 'Sample inputs'))
        object.__setattr__(self, 'output', real(self.output, 'Evaluator output'))


@dataclass(frozen=True, slots=True)
class SensitivityRecord:
    """Pearson correlation of one parameter with the study output and its rank.

    Attributes
    ----------
    parameter : str
        Parameter name.
    correlation : float or None
        Pearson coefficient in [-1, 1] [1]_; None when undefined because the input
        or the output is constant [2]_.
    absolute_correlation : float or None
        abs(correlation), or None when undefined.
    rank : int or None
        1 for the largest magnitude; None when undefined.
    undefined_reason : str or None
        ``'constant_input_or_output'`` exactly when the coefficient is undefined.

    Raises
    ------
    DomainError
        If the fields are inconsistent.

    Notes
    -----
    Pearson is a linear screening coefficient, not a global sensitivity measure.

    References
    ----------
    .. [1] NIST-PEARSON. NIST Dataplot reference manual, "Correlation": centered-product
           (Pearson) correlation definition.
    .. [2] SCIPY-1.14.1-PEARSON. SciPy 1.14.1, scipy.stats.pearsonr: Pearson coefficient
           and its undefined value for constant input.
    """
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
    """Sampled uncertainty study, replayable from its retained samples.

    On construction the record validates the method, controls and parameter ranges,
    checks that every sample lies inside the ranges, recomputes the summaries and
    sensitivities from the samples and refuses any value they do not reproduce.

    Attributes
    ----------
    method : {'latin-hypercube', 'monte-carlo'}
        Sampling design [1]_ [2]_.
    seed : int
        Seed of the pseudo-random generator.
    sample_count : int
        Number of samples, in [2, 100000].
    samples : tuple of SweepPoint
        Evaluated samples.
    output_mean : float
        Sample mean of the outputs [3]_.
    output_standard_deviation : float
        Sample standard deviation with the n - 1 divisor [3]_.
    output_percentiles : Mapping of str to float
        ``'p05'``, ``'p50'`` and ``'p95'`` linear sample quantiles [4]_.
    sensitivity : tuple of SensitivityRecord
        Pearson records ranked by magnitude; undefined records last [5]_.
    parameters : tuple of ParameterRange
        Sampled parameter ranges.
    sampling_numerics, summary_numerics : str
        Identifiers of the sampling and summary arithmetic policies.
    warnings : tuple of str
        Method-scope warnings; derived from ``method`` when empty.

    Raises
    ------
    DomainError
        If a control or range is invalid, a sample lies outside its range, or a
        summary, sensitivity or warning does not reproduce.

    References
    ----------
    .. [1] MCKAY-1979. M. D. McKay, R. J. Beckman and W. J. Conover, "A Comparison of
           Three Methods for Selecting Values of Input Variables in the Analysis of
           Output from a Computer Code", Technometrics 21(2), 239-245, 1979.
    .. [2] NASA/SP-2011-3421, "Probabilistic Risk Assessment Procedures Guide for NASA
           Managers and Practitioners", NASA, 2011 (Monte Carlo and LHS uncertainty
           propagation).
    .. [3] NIST TN 1297, "Guidelines for Evaluating and Expressing the Uncertainty of
           NIST Measurement Results", 1994 (Type A evaluation: sample mean and standard
           deviation).
    .. [4] NUMPY-2.2-LINEAR-QUANTILE. NumPy 2.2, numpy.quantile: linear (Hyndman-Fan
           type 7) sample quantile definition.
    .. [5] NIST-PEARSON. NIST Dataplot reference manual, "Correlation": centered-product
           (Pearson) correlation definition.
    """
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
        for name, value in {'method': method, 'seed': seed, 'sample_count': count, 'parameters': parameters,
                                'samples': samples, 'sensitivity': sensitivity, 'output_mean': mean,
                                'output_standard_deviation': deviation, 'output_percentiles': percentiles,
                                'warnings': warnings}.items():
            object.__setattr__(self, name, value)


def _validate_ranges(parameters):
    """Return parameter ranges as a non-empty tuple with unique names.

    Parameters
    ----------
    parameters : iterable of ParameterRange

    Returns
    -------
    tuple of ParameterRange

    Raises
    ------
    DomainError
        If the collection is empty, has a wrong type or duplicate names.
    """
    parameters = records(parameters, ParameterRange, 'Parameter ranges')
    if not parameters or len({p.name for p in parameters}) != len(parameters):
        raise DomainError('Parameter ranges must be nonempty with unique names.', field='Parameter ranges')
    return parameters


def _count(value, maximum, label):
    """Return an integer count in ``[2, maximum]``.

    Parameters
    ----------
    value : int
        Candidate count.
    maximum : int
        Upper limit.
    label : str
        Label used in errors.

    Returns
    -------
    int

    Raises
    ------
    DomainError
        If ``value`` is not an integer in range.
    """
    value = integer(value, label)
    if not 2 <= value <= maximum:
        raise DomainError(f'{label} must be between 2 and {maximum}.', field=label)
    return value


def _controls(sample_count, seed):
    """Return the validated ``(sample_count, seed)`` pair.

    Parameters
    ----------
    sample_count : int
        Number of samples, in [2, 100000].
    seed : int
        Integer seed.

    Returns
    -------
    tuple of int

    Raises
    ------
    DomainError
        If either value is invalid.
    """
    return _count(sample_count, 100000, 'Sample count'), integer(seed, 'Seed')


def _method(value):
    """Return the normalized sampling method name.

    Parameters
    ----------
    value : str
        ``'latin-hypercube'`` or ``'monte-carlo'``; case and surrounding whitespace
        are ignored.

    Returns
    -------
    str

    Raises
    ------
    DomainError
        If the method is unknown.
    """
    value = text(value, 'Uncertainty method').strip().lower()
    if value not in ('latin-hypercube', 'monte-carlo'):
        raise DomainError('Uncertainty method must be latin-hypercube or monte-carlo.', field='Uncertainty method')
    return value


def _affine(parameter, probability):
    # A convex combination cannot overflow its finite bounds. Round once; an
    # adjacent-to-zero grid cell may legitimately round a nonzero draw to zero.
    """Map a probability in [0, 1] onto the parameter bounds, rounding once.

    The convex combination (1 - p) * minimum + p * maximum is evaluated exactly and
    cannot overflow its finite bounds. Near zero, a nonzero draw may legitimately
    round to zero.

    Parameters
    ----------
    parameter : ParameterRange
    probability : Fraction

    Returns
    -------
    float
    """
    return float((1-probability)*Fraction(parameter.minimum)+probability*Fraction(parameter.maximum))


def latin_hypercube_samples(
    parameters: tuple[ParameterRange, ...], *, sample_count: int, seed: int
) -> tuple[Mapping[str, float], ...]:
    """Generate one Latin hypercube design with random jitter inside each stratum.

    Each parameter's range is divided into ``sample_count`` equal-probability
    half-open latent strata; one uniform point is drawn in every stratum and the
    column is randomly permuted [1]_.

    Parameters
    ----------
    parameters : tuple of ParameterRange
        Non-empty ranges with unique names.
    sample_count : int
        Number of samples, in [2, 100000].
    seed : int
        Seed of the pseudo-random generator.

    Returns
    -------
    tuple of Mapping
        One immutable name-to-value mapping per sample.

    Raises
    ------
    DomainError
        If a range, the count or the seed is invalid.

    Notes
    -----
    LHS rows are dependent, so independent-sample standard-error formulas are not
    qualified for this design.

    References
    ----------
    .. [1] MCKAY-1979. M. D. McKay, R. J. Beckman and W. J. Conover, "A Comparison of
           Three Methods for Selecting Values of Input Variables in the Analysis of
           Output from a Computer Code", Technometrics 21(2), 239-245, 1979.
    """
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
    """Generate seeded independent uniform draws over the declared bounds.

    Parameters
    ----------
    parameters : tuple of ParameterRange
        Non-empty ranges with unique names.
    sample_count : int
        Number of samples, in [2, 100000].
    seed : int
        Seed of the pseudo-random generator.

    Returns
    -------
    tuple of Mapping
        One immutable name-to-value mapping per sample.

    Raises
    ------
    DomainError
        If a range, the count or the seed is invalid.

    Notes
    -----
    Simple random (Monte Carlo) sampling of the input ranges [1]_.

    References
    ----------
    .. [1] NASA/SP-2011-3421, "Probabilistic Risk Assessment Procedures Guide for NASA
           Managers and Practitioners", NASA, 2011 (Monte Carlo and LHS uncertainty
           propagation).
    """
    parameters = _validate_ranges(parameters)
    sample_count, seed = _controls(sample_count, seed)
    random = Random(seed)
    return tuple(_FrozenValues((p.name, _affine(p, Fraction(random.random()))) for p in parameters)
                 for _ in range(sample_count))


def _evaluator(value):
    """Return ``value`` if it is callable.

    Raises
    ------
    DomainError
        If ``value`` is not callable.
    """
    if not callable(value):
        raise DomainError('Evaluator must be callable.', field='Evaluator')
    return value


def parameter_sweep(
    parameter: ParameterRange, *, point_count: int,
    evaluator: Callable[[Mapping[str, float]], float],
    fixed_inputs: Mapping[str, float] | None = None,
) -> tuple[SweepPoint, ...]:
    """Evaluate a model at evenly spaced points of one parameter, endpoints included.

    Parameters
    ----------
    parameter : ParameterRange
        Swept parameter.
    point_count : int
        Number of points, in [2, 10000].
    evaluator : callable
        ``evaluator(inputs) -> float`` called with an immutable input mapping.
    fixed_inputs : Mapping of str to float, optional
        Other inputs held constant; the swept parameter overrides a fixed value of
        the same name.

    Returns
    -------
    tuple of SweepPoint

    Raises
    ------
    DomainError
        If an argument is invalid or the evaluator returns a non-finite value.
    """
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
    """Return the linear sample quantile of ``sorted_values`` at ``probability``.

    References
    ----------
    .. [1] NUMPY-2.2-LINEAR-QUANTILE. NumPy 2.2, numpy.quantile: linear (Hyndman-Fan
           type 7) sample quantile definition.
    """
    return _linear_quantile(sorted_values, probability)


def _correlation(left, right):
    """Return the Pearson correlation of two samples, or None when undefined.

    References
    ----------
    .. [1] NIST-PEARSON. NIST Dataplot reference manual, "Correlation": centered-product
           (Pearson) correlation definition.
    .. [2] SCIPY-1.14.1-PEARSON. SciPy 1.14.1, scipy.stats.pearsonr: Pearson coefficient
           and its undefined value for constant input.
    """
    return _pearson(left, right)


def _summaries(samples):
    """Return the output mean, sample standard deviation and 5/50/95 % quantiles.

    Parameters
    ----------
    samples : tuple of SweepPoint

    Returns
    -------
    tuple
        ``(mean, standard_deviation, percentiles)``; the deviation uses the n - 1
        divisor [1]_ and the quantiles the linear definition [2]_.

    References
    ----------
    .. [1] NIST TN 1297, "Guidelines for Evaluating and Expressing the Uncertainty of
           NIST Measurement Results", 1994 (Type A evaluation: sample mean and standard
           deviation).
    .. [2] NUMPY-2.2-LINEAR-QUANTILE. NumPy 2.2, numpy.quantile: linear (Hyndman-Fan
           type 7) sample quantile definition.
    """
    outputs = tuple(sample.output for sample in samples)
    ledger = _sample_ledger(outputs)
    return (_finite_rational(ledger.mean, 'Output mean'),
            _sqrt_output(ledger.variance, 'Output standard deviation'),
            _FrozenValues((key, _percentile(outputs, probability)) for key, probability in
                          (('p05',Fraction(1,20)), ('p50',Fraction(1,2)), ('p95',Fraction(19,20)))))


def _sensitivities(parameters, samples):
    """Return Pearson sensitivity records ranked by decreasing magnitude.

    Ties keep the declared parameter order. Undefined coefficients follow all defined
    ones, including a defined zero, and carry no rank.

    Parameters
    ----------
    parameters : tuple of ParameterRange
    samples : tuple of SweepPoint

    Returns
    -------
    tuple of SensitivityRecord

    References
    ----------
    .. [1] NIST-PEARSON. NIST Dataplot reference manual, "Correlation": centered-product
           (Pearson) correlation definition.
    .. [2] SCIPY-1.14.1-PEARSON. SciPy 1.14.1, scipy.stats.pearsonr: Pearson coefficient
           and its undefined value for constant input.
    """
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
    """Return the scope warnings for the given sampling method.

    Parameters
    ----------
    method : str
        Normalized method name.

    Returns
    -------
    tuple of str
    """
    common = ('Empirical sample moments and percentiles do not establish confidence coverage or physical validation.',
              'Pearson is a linear screening coefficient, not a global sensitivity measure; constant vectors are undefined.')
    return common + (('LHS rows are dependent; IID standard-error formulas are not qualified for this design.',)
                     if method == 'latin-hypercube' else ())


def propagate_uncertainty(
    parameters: tuple[ParameterRange, ...], *, evaluator: Callable[[Mapping[str, float]], float],
    sample_count: int = 200, seed: int = 42, method: str = 'latin-hypercube',
) -> UncertaintyResult:
    """Propagate uniform input ranges through a model by sampling.

    Parameters
    ----------
    parameters : tuple of ParameterRange
        Non-empty ranges with unique names.
    evaluator : callable
        ``evaluator(inputs) -> float`` called with an immutable input mapping.
    sample_count : int, default 200
        Number of samples, in [2, 100000].
    seed : int, default 42
        Seed of the pseudo-random generator.
    method : {'latin-hypercube', 'monte-carlo'}, default 'latin-hypercube'
        Sampling design [1]_ [2]_.

    Returns
    -------
    UncertaintyResult
        Samples, empirical mean, standard deviation and percentiles [3]_ [4]_, and
        ranked Pearson sensitivities [5]_.

    Raises
    ------
    DomainError
        If an argument is invalid or the evaluator returns a non-finite value.

    Notes
    -----
    Empirical moments and percentiles do not establish confidence coverage or
    physical validation. Pearson is a linear screening coefficient, not a global
    sensitivity measure, and is undefined for constant vectors.

    References
    ----------
    .. [1] MCKAY-1979. M. D. McKay, R. J. Beckman and W. J. Conover, "A Comparison of
           Three Methods for Selecting Values of Input Variables in the Analysis of
           Output from a Computer Code", Technometrics 21(2), 239-245, 1979.
    .. [2] NASA/SP-2011-3421, "Probabilistic Risk Assessment Procedures Guide for NASA
           Managers and Practitioners", NASA, 2011 (Monte Carlo and LHS uncertainty
           propagation).
    .. [3] NIST TN 1297, "Guidelines for Evaluating and Expressing the Uncertainty of
           NIST Measurement Results", 1994 (Type A evaluation: sample mean and standard
           deviation).
    .. [4] NUMPY-2.2-LINEAR-QUANTILE. NumPy 2.2, numpy.quantile: linear (Hyndman-Fan
           type 7) sample quantile definition.
    .. [5] NIST-PEARSON. NIST Dataplot reference manual, "Correlation": centered-product
           (Pearson) correlation definition.
    """
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
