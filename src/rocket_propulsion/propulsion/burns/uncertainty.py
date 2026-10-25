"""Traceable uncertainty propagation for pressure-fed propulsion histories.

Bounded, sourced input marginals are sampled by Latin hypercube [1]_ or plain
Monte Carlo [3]_, optionally coupled through a Gaussian copula, and propagated through
the deterministic blowdown model. Uncertainty definition stays separate from
astrodynamics: every realization is a physically closed propulsion artifact that
a downstream trajectory propagator (for example Sidera) can consume without
resampling or altering engine physics. Every input is classified as aleatory
or epistemic [4]_.

Ensemble summaries are empirical statistics of the retained realizations. They
are engineering uncertainty estimates, not flight certification, and the models
are not validated against hardware by this module [2]_.

References
----------
.. [1] MCKAY-1979. M. D. McKay, R. J. Beckman and W. J. Conover, "A Comparison
       of Three Methods for Selecting Values of Input Variables in the Analysis
       of Output from a Computer Code", Technometrics 21(2), 239-245, 1979.
.. [2] NASA-STD-7009B, "Standard for Models and Simulations".
       https://standards.nasa.gov/standard/NASA/NASA-STD-7009
.. [3] NASA/SP-2011-3421, "Probabilistic Risk Assessment Procedures Guide for
       NASA Managers and Practitioners". https://ntrs.nasa.gov/citations/20120001369
.. [4] NASA/SP-2009-569, "Bayesian Inference for NASA Probabilistic Risk and
       Reliability Analysis". https://ntrs.nasa.gov/citations/20090023159
"""

from __future__ import annotations

from dataclasses import dataclass, replace
from enum import Enum
from fractions import Fraction
from math import inf, isfinite, nextafter
from random import Random
from statistics import NormalDist

from rocket_propulsion.core.errors import ConvergenceError, DomainError
from rocket_propulsion.core.validation import (
    coerce_fields,
    integer,
    name_value_pairs,
    real,
    real_matrix,
    records,
    strings,
    text,
)

from .blowdown import (
    BlowdownTankModel,
    PressurePerformanceLaw,
    _scaled_operating_point,
    simulate_pressure_fed_blowdown,
)
from .marginal_numerics import (
    _SAMPLING_POLICY,
    _GaussianProbability,
    _normal_quantile,
    marginal_value,
)
from .models import EngineOperatingPoint
from .numeric_outputs import _finite_rational, _fraction
from .probability_numerics import _sqrt_output
from .sample_statistics import (
    _INTERVAL_KIND,
    _NORMAL_PROXY_MULTIPLIER,
    _SUMMARY_POLICY,
    _affine_sqrt_output,
    _deviation_relative_change,
    _linear_quantile,
    _pearson,
    _relative_exact,
    _sample_ledger,
)
from .tabulated import TabulatedBurnArtifact

UNCERTAINTY_REFERENCE_IDS = (
    "NASA-STD-7009B",
    "NASA-SP-2011-3421",
    "NASA-SP-2009-569",
    "NASA-SP-8112",
    "NASA-CR-131400",
)


class UncertaintyClass(str, Enum):
    """Source class assigned to every uncertain input [1]_ [2]_.

    Attributes
    ----------
    ALEATORY
        Modeled variability between realizations (irreducible at the stated
        level of modeling).
    EPISTEMIC
        Incomplete knowledge that may be reduced with additional data.

    References
    ----------
    .. [1] NASA-STD-7009B, "Standard for Models and Simulations".
           https://standards.nasa.gov/standard/NASA/NASA-STD-7009
    .. [2] NASA/SP-2009-569, "Bayesian Inference for NASA Probabilistic Risk and
           Reliability Analysis". https://ntrs.nasa.gov/citations/20090023159
    """

    ALEATORY = "aleatory"
    EPISTEMIC = "epistemic"


class DistributionKind(str, Enum):
    """Supported bounded marginal distributions [3]_.

    Attributes
    ----------
    UNIFORM
        Constant density on ``[lower, upper]``.
    TRIANGULAR
        Triangular density with mode ``nominal`` [1]_.
    TRUNCATED_NORMAL
        Normal density with mean ``nominal`` and ``standard_deviation``,
        renormalized on ``[lower, upper]``.
    LOG_UNIFORM
        Uniform density of ``ln x`` on positive bounds [2]_.

    References
    ----------
    .. [1] NIST-TRIPPF. NIST Dataplot, triangular percent point function.
    .. [2] SCIPY-1.14.1-LOGUNIFORM. SciPy, log-uniform CDF and quantile formulas.
    .. [3] NASA/SP-2011-3421, "Probabilistic Risk Assessment Procedures Guide for
           NASA Managers and Practitioners". https://ntrs.nasa.gov/citations/20120001369
    """

    UNIFORM = "uniform"
    TRIANGULAR = "triangular"
    TRUNCATED_NORMAL = "truncated-normal"
    LOG_UNIFORM = "log-uniform"


@dataclass(frozen=True, slots=True)
class UncertainParameter:
    """One bounded, sourced propulsion uncertainty input.

    Parameters
    ----------
    name : str
        Model input name; nonempty.
    nominal : float
        Nominal value inside ``[lower, upper]``. It is also the triangular mode
        and the truncated-normal mean.
    lower, upper : float
        Finite bounds with ``lower < upper``. Log-uniform marginals require
        ``lower > 0``.
    distribution : DistributionKind
        Marginal distribution [2]_.
    uncertainty_class : UncertaintyClass
        Aleatory or epistemic classification [1]_.
    source : str
        Nonempty provenance of the bounds.
    rationale : str
        Nonempty justification. Source and rationale are mandatory so that an
        arbitrary engineering percentage cannot enter unnoticed.
    standard_deviation : float or None, optional
        Positive finite sigma, required for and only allowed with a
        truncated-normal marginal.

    Raises
    ------
    DomainError
        If any of the above constraints is violated.

    References
    ----------
    .. [1] NASA-STD-7009B, "Standard for Models and Simulations".
           https://standards.nasa.gov/standard/NASA/NASA-STD-7009
    .. [2] NASA/SP-2011-3421, "Probabilistic Risk Assessment Procedures Guide for
           NASA Managers and Practitioners". https://ntrs.nasa.gov/citations/20120001369
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
        coerce_fields(self, "Uncertain parameter")
        if not isinstance(self.distribution, DistributionKind) or not isinstance(self.uncertainty_class, UncertaintyClass):
            raise DomainError("Uncertain parameter requires declared distribution and uncertainty enums.")
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
        """Return the marginal quantile at ``probability`` (inverse-CDF transform).

        Parameters
        ----------
        probability : float
            Probability in the open interval ``(0, 1)``.

        Returns
        -------
        float
            Correctly rounded quantile in ``[lower, upper]``; see
            :func:`~rocket_propulsion.propulsion.burns.marginal_numerics.marginal_value`
            for the distribution formulas [1]_ [2]_ [3]_. Inverse-CDF
            sampling of the marginals follows [4]_.

        Raises
        ------
        DomainError
            If ``probability`` is not finite or not in ``(0, 1)``.

        References
        ----------
        .. [1] NIST-TRIPPF. NIST Dataplot, triangular percent point function.
        .. [2] SCIPY-1.14.1-LOGUNIFORM. SciPy, log-uniform CDF and quantile formulas.
        .. [3] NIST-DLMF-7.8. NIST DLMF, error-function and Mills ratio bounds.
        .. [4] NASA/SP-2011-3421, "Probabilistic Risk Assessment Procedures Guide for
               NASA Managers and Practitioners". https://ntrs.nasa.gov/citations/20120001369
        """

        probability = real(probability, "Sampling probability")
        if not 0.0 < probability < 1.0:
            raise DomainError("Sampling probability must be finite and in (0, 1).")
        return marginal_value(self, probability)


@dataclass(frozen=True, slots=True)
class CorrelationModel:
    """Gaussian-copula correlation model for ordered uncertain inputs.

    Parameters
    ----------
    parameter_names : tuple of str
        At least two unique input names, in matrix order.
    matrix : tuple of tuple of float
        Correlation matrix of the latent standard normal variates. It must be
        square, symmetric, strictly positive definite, have a unit diagonal and
        off-diagonal entries in ``[-1, 1]``.
    source : str
        Nonempty provenance of the coefficients.
    rationale : str
        Nonempty justification. Source and rationale are mandatory because
        guessed correlations can dominate ensemble tails even when every
        marginal is individually reasonable [1]_ [2]_.

    Raises
    ------
    DomainError
        If the names or matrix violate the constraints above.

    Notes
    -----
    Positive definiteness is decided with an exact rational ``L D L^T``
    factorization, so the decision is free of rounding error.

    References
    ----------
    .. [1] NASA-STD-7009B, "Standard for Models and Simulations".
           https://standards.nasa.gov/standard/NASA/NASA-STD-7009
    .. [2] NASA/SP-2011-3421, "Probabilistic Risk Assessment Procedures Guide for
           NASA Managers and Practitioners". https://ntrs.nasa.gov/citations/20120001369
    """

    parameter_names: tuple[str, ...]
    matrix: tuple[tuple[float, ...], ...]
    source: str
    rationale: str

    def __post_init__(self) -> None:
        coerce_fields(self, "Correlation model")
        object.__setattr__(self,"parameter_names",strings(self.parameter_names,"Correlation parameter names"))
        object.__setattr__(self,"matrix",real_matrix(self.matrix,"Correlation matrix"))
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
            if self.matrix[index][index] != 1.0:
                raise DomainError("Correlation matrix diagonal must equal one.")
            for other in range(index):
                if self.matrix[index][other] != self.matrix[other][index]:
                    raise DomainError("Correlation matrix must be symmetric.")
        _cholesky(self.matrix)


@dataclass(frozen=True, slots=True)
class UncertaintyDraw:
    """One deterministic, indexed uncertainty realization [1]_.

    Parameters
    ----------
    index : int
        Nonnegative realization index.
    values : tuple of (str, float)
        Sampled value of each named input; names must be unique.

    Raises
    ------
    DomainError
        If the index is negative or a name repeats.

    References
    ----------
    .. [1] NASA/SP-2011-3421, "Probabilistic Risk Assessment Procedures Guide for
           NASA Managers and Practitioners". https://ntrs.nasa.gov/citations/20120001369
    """

    index: int
    values: tuple[tuple[str, float], ...]

    def __post_init__(self) -> None:
        coerce_fields(self,"Uncertainty draw")
        object.__setattr__(self,"values",name_value_pairs(self.values,"Uncertainty draw values"))
        if self.index < 0 or len({name for name,value in self.values}) != len(self.values):
            raise DomainError("Uncertainty draw requires a nonnegative index and unique names.")

    def as_dict(self) -> dict[str, float]:
        """Return this immutable draw as a fresh name-to-value mapping.

        Returns
        -------
        dict of str to float
            A new dictionary; mutating it does not affect the draw.
        """

        return dict(self.values)


@dataclass(frozen=True, slots=True)
class PropulsionMetricStatistics:
    """Empirical statistics of one metric and an unqualified IID-normal proxy.

    The statistics summarize sampling-based propagation output [5]_.

    Parameters
    ----------
    metric : str
        Metric name.
    mean, standard_deviation : float
        Sample mean and ``n - 1`` sample standard deviation ``s`` [1]_.
    mean_standard_error : float
        ``s / sqrt(n)`` [2]_.
    mean_ci95_lower, mean_ci95_upper : float
        ``mean -/+ 1.96 s / sqrt(n)``. The normal multiplier assumes a known
        standard deviation and independent draws [3]_; with an estimated ``s``
        and dependent Latin-hypercube rows it establishes neither 95 % coverage
        nor a design-specific sampling error.
    minimum, maximum : float
        Sample extrema.
    p05, p50, p95 : float
        Linear-interpolation sample percentiles [4]_.
    mean_interval_kind : str, optional
        Fixed label of the proxy interval kind.
    coverage_established : bool, optional
        Always ``False``.

    Raises
    ------
    DomainError
        If dispersion is negative, percentiles are not ordered within the
        extrema, the proxy does not bracket the mean, or the record claims
        established coverage.

    References
    ----------
    .. [1] NIST TN 1297, "Guidelines for Evaluating and Expressing the Uncertainty
           of NIST Measurement Results", 1994 (Type A evaluation).
           https://www.nist.gov/pml/nist-technical-note-1297
    .. [2] NIST-MEAN-INTERVAL. NIST, confidence limits for the mean.
    .. [3] NIST-KNOWN-SIGMA-INTERVAL. NIST, known and estimated population
           deviation in confidence intervals.
    .. [4] NUMPY-2.2-LINEAR-QUANTILE. NumPy, ``numpy.quantile``, linear
           (Hyndman-Fan type 7) sample quantile definition.
    .. [5] NASA/SP-2011-3421, "Probabilistic Risk Assessment Procedures Guide for
           NASA Managers and Practitioners". https://ntrs.nasa.gov/citations/20120001369
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
    mean_interval_kind: str = _INTERVAL_KIND
    coverage_established: bool = False

    def __post_init__(self) -> None:
        coerce_fields(self,"Sample metric statistics")
        if self.standard_deviation<0 or self.mean_standard_error<0:
            raise DomainError("Sample dispersion must be nonnegative.")
        if not self.minimum<=self.p05<=self.p50<=self.p95<=self.maximum:
            raise DomainError("Sample percentiles must be ordered within the extrema.")
        if not self.minimum<=self.mean<=self.maximum:
            raise DomainError("Sample mean must lie within the extrema.")
        if not self.mean_ci95_lower<=self.mean<=self.mean_ci95_upper:
            raise DomainError("Mean proxy endpoints must bracket the mean.")
        if self.mean_interval_kind!=_INTERVAL_KIND or self.coverage_established is not False:
            raise DomainError("Sample summary does not establish confidence coverage.")


@dataclass(frozen=True, slots=True)
class PropulsionSensitivity:
    """Pearson screening sensitivity for one parameter/metric pair.

    Correlation is a screening statistic [1]_ [2]_. It is not proof of causality nor
    a complete global-sensitivity decomposition, and it can miss nonlinear
    interactions.

    Parameters
    ----------
    metric, parameter : str
        Output metric and input parameter names.
    correlation : float or None
        Pearson coefficient in ``[-1, 1]``; ``None`` when either vector is
        constant.
    absolute_correlation : float or None
        ``abs(correlation)``; ``None`` when undefined.
    rank : int or None
        Positive rank by absolute correlation within the metric; ``None`` when
        undefined.
    undefined_reason : str or None, optional
        ``"constant_input_or_output"`` exactly when ``correlation`` is ``None``.

    Raises
    ------
    DomainError
        If the defined and undefined fields are inconsistent.

    References
    ----------
    .. [1] NIST-PEARSON. NIST Dataplot, centered-product correlation definition.
    .. [2] NASA/SP-2011-3421, "Probabilistic Risk Assessment Procedures Guide for
           NASA Managers and Practitioners". https://ntrs.nasa.gov/citations/20120001369
    """

    metric: str
    parameter: str
    correlation: float | None
    absolute_correlation: float | None
    rank: int | None
    undefined_reason: str | None = None

    def __post_init__(self) -> None:
        coerce_fields(self,"Propulsion sensitivity")
        if self.undefined_reason is not None:
            object.__setattr__(self,"undefined_reason",text(self.undefined_reason,"Undefined sensitivity reason"))
        if self.correlation is None:
            if self.absolute_correlation is not None or self.rank is not None or self.undefined_reason!='constant_input_or_output':
                raise DomainError("Undefined Pearson sensitivity requires null magnitude/rank and a constant-vector reason.")
        else:
            object.__setattr__(self,"rank",integer(self.rank,"Sensitivity rank"))
            if not -1<=self.correlation<=1 or self.absolute_correlation!=abs(self.correlation):
                raise DomainError("Pearson coefficient and absolute magnitude disagree.")
            if self.rank<1 or self.undefined_reason is not None:
                raise DomainError("Defined sensitivity needs a positive rank and no undefined reason.")


@dataclass(frozen=True, slots=True)
class SamplingConvergenceEvidence:
    """First-half versus full-ensemble stability diagnostic.

    The diagnostic does not certify probabilistic convergence. It makes sample
    adequacy visible and allows a workflow to demand more realizations [1]_ [2]_.

    Parameters
    ----------
    metric : str
        Metric under test.
    half_sample_count, full_sample_count : int
        ``full_sample_count >= 4`` and ``half_sample_count`` equal to its floor
        half.
    mean_relative_change : float
        ``|m_half - m_full| / max(|m_half|, |m_full|)``.
    standard_deviation_relative_change : float
        ``|s_half - s_full| / max(s_half, s_full)``.
    requested_tolerance : float
        Positive acceptance tolerance.
    converged : bool
        ``True`` when both changes are within the tolerance.

    Raises
    ------
    DomainError
        If counts, changes or tolerance are out of range, or the flag
        contradicts the displayed changes.

    References
    ----------
    .. [1] NASA-STD-7009B, "Standard for Models and Simulations".
           https://standards.nasa.gov/standard/NASA/NASA-STD-7009
    .. [2] NASA/SP-2011-3421, "Probabilistic Risk Assessment Procedures Guide for
           NASA Managers and Practitioners". https://ntrs.nasa.gov/citations/20120001369
    """

    metric: str
    half_sample_count: int
    full_sample_count: int
    mean_relative_change: float
    standard_deviation_relative_change: float
    requested_tolerance: float
    converged: bool

    def __post_init__(self) -> None:
        coerce_fields(self,"Sampling stability")
        if self.full_sample_count<4 or self.half_sample_count!=self.full_sample_count//2:
            raise DomainError("Sampling stability requires half/full sample counts of at least two/four.")
        changes=(self.mean_relative_change,self.standard_deviation_relative_change)
        if any(not 0<=value<=2 for value in changes) or self.requested_tolerance<=0:
            raise DomainError("Sampling stability requires bounded nonnegative changes and positive tolerance.")
        if type(self.converged) is not bool:
            raise DomainError("Sampling stability flag must be boolean.")
        # Display equality can hide either exact threshold side; the owning
        # result replays the decision from its retained sample values.
        shown=max(changes)
        if self.converged and shown>self.requested_tolerance or not self.converged and shown<self.requested_tolerance:
            raise DomainError("Sampling stability flag contradicts its displayed evidence.")


@dataclass(frozen=True, slots=True)
class BlowdownUncertaintySample:
    """One physically closed blowdown realization and its exchange artifact.

    The realization is one run of the pressure-fed blowdown model [1]_ [2]_.

    Parameters
    ----------
    index : int
        Nonnegative realization index.
    inputs : tuple of (str, float)
        Sampled inputs with unique names.
    duration_s, delivered_impulse_n_s, consumed_propellant_kg : float
        Burn metrics; must equal the artifact's values.
    thrust_centroid_time_s : float
        Thrust-weighted mean time; must equal the artifact's value.
    final_supply_pressure_pa : float
        Positive supply pressure at cutoff [Pa].
    cutoff_reason : str
        Reason the blowdown simulation terminated.
    artifact : TabulatedBurnArtifact
        Complete thrust and mass-flow history of the realization.

    Raises
    ------
    DomainError
        If the record is inconsistent with its artifact or a constraint above
        is violated.

    References
    ----------
    .. [1] NASA SP-8112, "Pressurization Systems for Liquid Rockets".
           https://ntrs.nasa.gov/citations/19760015212
    .. [2] NASA CR-131400, "Reliability Model of a Monopropellant Auxiliary
           Propulsion System". https://ntrs.nasa.gov/citations/19730012094
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

    def __post_init__(self) -> None:
        coerce_fields(self,"Blowdown uncertainty sample")
        object.__setattr__(self,"inputs",name_value_pairs(self.inputs,"Blowdown sample inputs"))
        if self.index<0 or len(dict(self.inputs))!=len(self.inputs):
            raise DomainError("Blowdown sample needs a nonnegative index and unique inputs.")
        if not isinstance(self.artifact,TabulatedBurnArtifact) or self.final_supply_pressure_pa<=0:
            raise DomainError("Blowdown sample requires an artifact and positive final pressure.")
        expected=(self.artifact.duration_s,self.artifact.delivered_total_impulse_n_s,
                  self.artifact.consumed_propellant_kg,self.artifact.thrust_centroid_time_s)
        if (self.duration_s,self.delivered_impulse_n_s,self.consumed_propellant_kg,self.thrust_centroid_time_s)!=expected:
            raise DomainError("Blowdown sample metrics do not reproduce its artifact.")


@dataclass(frozen=True, slots=True)
class BlowdownUncertaintyResult:
    """Traceable blowdown uncertainty ensemble ready for trajectory propagation.

    Construction recomputes every statistic, sensitivity and stability record
    from the retained samples and rejects a record that does not reproduce, so a
    persisted result cannot carry summaries inconsistent with its realizations.
    The statistics are engineering estimates for credibility assessment [2]_ of
    sampling-based propagation [3]_ through pressure-fed blowdown histories
    [4]_.

    Parameters
    ----------
    method : {"latin-hypercube", "monte-carlo"}
        Sampling method [1]_.
    seed : int
        Random seed of the sampling stream.
    sample_count : int
        Number of realizations, 4 to 10000.
    aleatory_parameter_count, epistemic_parameter_count : int
        Counts by uncertainty class; must reproduce from ``parameters``.
    parameters : tuple of UncertainParameter
        Supported blowdown inputs with unique names.
    correlation : CorrelationModel or None
        Optional Gaussian-copula model in parameter order.
    samples : tuple of BlowdownUncertaintySample
        One sample per index ``0 .. sample_count - 1``, inside declared bounds.
    statistics : tuple of PropulsionMetricStatistics
        Per-metric empirical statistics.
    sensitivity : tuple of PropulsionSensitivity
        Pearson screening sensitivities.
    convergence : SamplingConvergenceEvidence
        Half/full stability of delivered impulse.
    reference_ids, warnings : tuple of str
        Source identifiers and interpretation warnings.
    sampling_numerics, summary_numerics : str, optional
        Fixed numerical-policy identifiers.

    Raises
    ------
    DomainError
        If any record does not reproduce from the retained samples or a
        constraint above is violated.

    References
    ----------
    .. [1] MCKAY-1979. M. D. McKay, R. J. Beckman and W. J. Conover, "A Comparison
           of Three Methods for Selecting Values of Input Variables in the Analysis
           of Output from a Computer Code", Technometrics 21(2), 239-245, 1979.
    .. [2] NASA-STD-7009B, "Standard for Models and Simulations".
           https://standards.nasa.gov/standard/NASA/NASA-STD-7009
    .. [3] NASA/SP-2011-3421, "Probabilistic Risk Assessment Procedures Guide for
           NASA Managers and Practitioners". https://ntrs.nasa.gov/citations/20120001369
    .. [4] NASA SP-8112, "Pressurization Systems for Liquid Rockets".
           https://ntrs.nasa.gov/citations/19760015212
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
    sampling_numerics: str = _SAMPLING_POLICY
    summary_numerics: str = _SUMMARY_POLICY

    def __post_init__(self) -> None:
        coerce_fields(self,"Blowdown uncertainty result")
        for name,kind in (("parameters",UncertainParameter),("samples",BlowdownUncertaintySample),
                          ("statistics",PropulsionMetricStatistics),("sensitivity",PropulsionSensitivity)):
            object.__setattr__(self,name,records(getattr(self,name),kind,f"Blowdown {name}"))
        for name in ("reference_ids","warnings"):
            object.__setattr__(self,name,strings(getattr(self,name),f"Blowdown {name}"))
        if self.sampling_numerics!=_SAMPLING_POLICY or self.summary_numerics!=_SUMMARY_POLICY:
            raise DomainError("Unsupported blowdown uncertainty numerical policy.")
        _validate_uncertainty_result(self)


def _cholesky(matrix: tuple[tuple[float, ...], ...]) -> tuple[tuple[float, ...], ...]:
    """Return the lower Cholesky factor of a correlation matrix.

    The ``L D L^T`` factorization is computed in exact rational arithmetic, so
    the positive-definiteness decision is exact. Only the final coefficients
    ``L_ij sqrt(D_j)`` are rounded, each once.

    Raises
    ------
    DomainError
        If the matrix is not strictly positive definite.
    """
    size = len(matrix)
    # Exact LDL^T decisions; only the final Cholesky coefficients are rounded.
    lower = [[Fraction(int(row == column)) for column in range(size)] for row in range(size)]
    diagonal=[]
    for row in range(size):
        for column in range(row):
            residual=_fraction(matrix[row][column])-sum((lower[row][k]*lower[column][k]*diagonal[k] for k in range(column)),Fraction(0))
            lower[row][column]=residual/diagonal[column]
        residual=_fraction(matrix[row][row])-sum((lower[row][k]**2*diagonal[k] for k in range(row)),Fraction(0))
        if residual <= 0:raise DomainError("Correlation matrix must be positive definite.")
        diagonal.append(residual)
    return tuple(tuple(
        (1 if lower[row][column] >= 0 else -1)*_sqrt_output(lower[row][column]**2*diagonal[column],"Cholesky coefficient")
        if column <= row else 0. for column in range(size)) for row in range(size))


def _open_random(random):
    """Return a variate from ``random.random()`` in the open interval ``(0, 1)``.

    The zero endpoint is redrawn up to 1024 times.

    Raises
    ------
    DomainError
        If the generator returns a value outside ``[0, 1)``.
    ConvergenceError
        If every draw returns zero.
    """
    for _ in range(1024):
        value=real(random.random(),"Unit random variate")
        if 0 < value < 1:return value
        if value != 0:raise DomainError("Random variate is outside [0, 1).")
    raise ConvergenceError("Random source repeatedly returned the zero endpoint.",field="Unit random variate")


def _stratum_probability(index, value, count):
    """Return the binary64 probability ``(index + value) / count`` inside its stratum.

    The result lies in the open unit interval and in the half-open stratum
    ``[index / count, (index + 1) / count)``; it is moved by one ulp when
    rounding would leave the stratum.

    Raises
    ------
    DomainError
        If the stratum contains no representable open probability.
    """
    lower,upper=Fraction(index,count),Fraction(index+1,count)
    exact=(index+_fraction(value))/count
    # This is an intermediate; zero may need one directed step into the cell.
    point=float(exact)
    if _fraction(point) >= upper or point >= 1:point=nextafter(point,-inf)
    if _fraction(point) < lower or point <= 0:point=nextafter(point,inf)
    if not 0 < point < 1 or not lower <= _fraction(point) < upper:
        raise DomainError("Stratum has no representable open probability.",field="Unit stratum probability")
    return point


def _unit_samples(
    dimension: int, *, sample_count: int, seed: int, method: str
) -> list[list[float]]:
    """Return ``sample_count`` rows of ``dimension`` open unit-interval variates.

    ``"latin-hypercube"`` places one variate in each of ``sample_count`` equal
    strata per column and permutes each column independently [1]_;
    ``"monte-carlo"`` draws independent variates.

    Raises
    ------
    DomainError
        If the method is not supported.

    References
    ----------
    .. [1] MCKAY-1979. M. D. McKay, R. J. Beckman and W. J. Conover, "A Comparison
           of Three Methods for Selecting Values of Input Variables in the Analysis
           of Output from a Computer Code", Technometrics 21(2), 239-245, 1979.
    """
    random = Random(seed)
    if method == "latin-hypercube":
        columns: list[list[float]] = []
        for _ in range(dimension):
            column = [_stratum_probability(index,_open_random(random),sample_count) for index in range(sample_count)]
            random.shuffle(column)
            columns.append(column)
        return [[columns[column][row] for column in range(dimension)] for row in range(sample_count)]
    if method == "monte-carlo":
        return [[_open_random(random) for _ in range(dimension)] for _ in range(sample_count)]
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

    Unit-interval variates are produced by Latin hypercube [1]_ or Monte Carlo
    sampling and mapped through each marginal's inverse CDF. With a correlation
    model, each row is mapped to independent standard normal variates
    ``u = Phi^-1(p)``, correlated as ``z = L u`` with the lower Cholesky factor
    ``L`` of the matrix, and passed to the marginals as the latent ``Phi(z)``.

    Parameters
    ----------
    parameters : tuple of UncertainParameter
        Nonempty inputs with unique names.
    sample_count : int
        Number of draws, 4 to 10000.
    seed : int
        Seed of the ``random.Random`` stream; equal seeds give equal draws.
    method : {"latin-hypercube", "monte-carlo"}, optional
        Sampling method (case-insensitive).
    correlation : CorrelationModel or None, optional
        Gaussian-copula model whose ``parameter_names`` match ``parameters``
        in order.

    Returns
    -------
    tuple of UncertaintyDraw
        Draws indexed ``0 .. sample_count - 1``.

    Raises
    ------
    DomainError
        If the inputs, count, method or correlation order are invalid.
    ConvergenceError
        If a certified quantile comparison cannot be resolved.

    Notes
    -----
    Latin-hypercube stratification is exact for independent marginals. With a
    correlation model the copula transform prioritizes the declared dependency
    structure, so the final marginal ranks are not claimed to remain an exact
    Latin hypercube. Sampling-based propagation follows [2]_; the draws carry
    no credibility claim beyond their declared sources [3]_.

    References
    ----------
    .. [1] MCKAY-1979. M. D. McKay, R. J. Beckman and W. J. Conover, "A Comparison
           of Three Methods for Selecting Values of Input Variables in the Analysis
           of Output from a Computer Code", Technometrics 21(2), 239-245, 1979.
    .. [2] NASA/SP-2011-3421, "Probabilistic Risk Assessment Procedures Guide for
           NASA Managers and Practitioners". https://ntrs.nasa.gov/citations/20120001369
    .. [3] NASA-STD-7009B, "Standard for Models and Simulations".
           https://standards.nasa.gov/standard/NASA/NASA-STD-7009
    """

    parameters=records(parameters,UncertainParameter,"Uncertainty parameters")
    sample_count=integer(sample_count,"Uncertainty sample count")
    seed=integer(seed,"Uncertainty seed")
    normalized_method=text(method,"Uncertainty method").strip().lower()
    if correlation is not None and not isinstance(correlation,CorrelationModel):
        raise DomainError("Uncertainty correlation requires a CorrelationModel.")
    if not parameters:
        raise DomainError("At least one propulsion uncertainty parameter is required.")
    if len({parameter.name for parameter in parameters}) != len(parameters):
        raise DomainError("Uncertainty parameter names must be unique.")
    if not 4 <= sample_count <= 10_000:
        raise DomainError("Propulsion uncertainty sample count must be between 4 and 10000.")
    unit_rows = _unit_samples(
        len(parameters), sample_count=sample_count, seed=seed, method=normalized_method
    )
    if correlation is not None:
        expected_names = tuple(parameter.name for parameter in parameters)
        if correlation.parameter_names != expected_names:
            raise DomainError("Correlation parameter order must match uncertainty parameters.")
        lower = _cholesky(correlation.matrix)
        for row_index, row in enumerate(unit_rows):
            independent = [
                _normal_quantile(value)
                for value in row
            ]
            correlated = [
                _finite_rational(sum((_fraction(lower[index][column])*_fraction(independent[column]) for column in range(index+1)),Fraction(0)),"Correlated Gaussian variate")
                for index in range(len(parameters))
            ]
            unit_rows[row_index] = [_GaussianProbability(value) for value in correlated]
    return tuple(
        UncertaintyDraw(
            index=index,
            values=tuple(
                (parameter.name, marginal_value(parameter,row[column]) if type(parameter) is UncertainParameter else _custom_parameter_value(parameter,row[column]))
                for column, parameter in enumerate(parameters)
            ),
        )
        for index, row in enumerate(unit_rows)
    )


def _custom_parameter_value(parameter, probability):
    """Evaluate a subclass ``value_at_probability`` and check its bounds.

    A Gaussian-copula latent value is converted to a float probability with
    ``statistics.NormalDist`` because custom callbacks take a float; a
    probability that rounds to 0 or 1 is rejected instead of being clipped.

    Raises
    ------
    DomainError
        If the probability is not representable in ``(0, 1)`` or the returned
        value is outside the declared bounds.
    """
    if isinstance(probability,_GaussianProbability):
        # Custom callbacks retain their public float interface; no silent tail clip.
        value=NormalDist().cdf(float(probability.z))
        if not 0 < value < 1:
            raise DomainError("Custom marginal callback requires a representable open probability.",field="Sampling probability")
        probability=value
    value=real(parameter.value_at_probability(probability),"Custom marginal value")
    if not parameter.lower <= value <= parameter.upper:
        raise DomainError("Custom marginal value is outside its declared bounds.")
    return value


def _percentile(values: list[float], probability: float) -> float:
    """Return the linear-interpolation sample percentile [1]_.

    References
    ----------
    .. [1] NUMPY-2.2-LINEAR-QUANTILE. NumPy, ``numpy.quantile``, linear
           (Hyndman-Fan type 7) sample quantile definition.
    """
    return _linear_quantile(values,probability)


def _statistics(metric: str, values: list[float]) -> PropulsionMetricStatistics:
    """Return exact-input empirical statistics of one metric.

    The mean, ``n - 1`` standard deviation [1]_, standard error ``s / sqrt(n)``
    and the proxy endpoints ``mean -/+ 1.96 s / sqrt(n)`` [2]_ are each rounded
    once from exact rationals; percentiles use linear interpolation [3]_.

    References
    ----------
    .. [1] NIST TN 1297, "Guidelines for Evaluating and Expressing the Uncertainty
           of NIST Measurement Results", 1994 (Type A evaluation).
           https://www.nist.gov/pml/nist-technical-note-1297
    .. [2] NIST-MEAN-INTERVAL. NIST, confidence limits for the mean.
    .. [3] NUMPY-2.2-LINEAR-QUANTILE. NumPy, ``numpy.quantile``, linear
           (Hyndman-Fan type 7) sample quantile definition.
    """
    ledger=_sample_ledger(values)
    mean = _finite_rational(ledger.mean,"Sample mean")
    deviation = _sqrt_output(ledger.variance,"Sample standard deviation")
    standard_error = _sqrt_output(ledger.variance/len(ledger.values),"Mean standard-error proxy")
    width_squared = ledger.variance/len(ledger.values)*_NORMAL_PROXY_MULTIPLIER**2
    return PropulsionMetricStatistics(
        metric=metric,
        mean=mean,
        standard_deviation=deviation,
        mean_standard_error=standard_error,
        mean_ci95_lower=_affine_sqrt_output(ledger.mean,width_squared,-1,"Mean proxy lower endpoint"),
        mean_ci95_upper=_affine_sqrt_output(ledger.mean,width_squared,1,"Mean proxy upper endpoint"),
        minimum=float(min(ledger.values)),
        p05=_percentile(ledger.values,Fraction(1,20)),
        p50=_percentile(ledger.values,Fraction(1,2)),
        p95=_percentile(ledger.values,Fraction(19,20)),
        maximum=float(max(ledger.values)),
    )


def _correlation(left: list[float], right: list[float]) -> float | None:
    """Return the Pearson coefficient [1]_, or ``None`` for a constant vector.

    References
    ----------
    .. [1] NIST-PEARSON. NIST Dataplot, centered-product correlation definition.
    """
    return _pearson(left,right)


def _relative_change(previous: float, current: float) -> float:
    """Return ``|previous - current| / max(|previous|, |current|)``."""
    return _relative_exact(previous,current)


def _sampling_convergence(values,tolerance):
    """Compare the first-half and full-sample mean and standard deviation.

    Both decisions are made on exact rationals: the mean passes when its
    relative change is within ``tolerance``; the standard deviation passes when
    ``min(v_half, v_full) / max(v_half, v_full) >= (1 - tolerance)**2``, which is
    equivalent to a relative standard-deviation change within ``tolerance``.

    Raises
    ------
    DomainError
        If fewer than four values are supplied or the tolerance is not positive.
    """
    tolerance=real(tolerance,"Sampling convergence tolerance")
    full=_sample_ledger(values)
    if len(full.values)<4 or tolerance<=0:
        raise DomainError("Sampling stability needs at least four values and positive tolerance.")
    half=_sample_ledger(full.values[:len(full.values)//2])
    scale=max(abs(half.mean),abs(full.mean))
    mean_ratio=abs(half.mean-full.mean)/scale if scale else Fraction(0)
    variance_scale=max(half.variance,full.variance)
    variance_ratio=min(half.variance,full.variance)/variance_scale if variance_scale else Fraction(1)
    exact_tolerance=_fraction(tolerance)
    deviation_pass=exact_tolerance>=1 or variance_ratio>=(1-exact_tolerance)**2
    return SamplingConvergenceEvidence('delivered_impulse_n_s',len(half.values),len(full.values),
        _finite_rational(mean_ratio,'Relative sample difference'),
        _deviation_relative_change(half.variance,full.variance),tolerance,
        mean_ratio<=exact_tolerance and deviation_pass)


_BLOWDOWN_METRICS=('duration_s','delivered_impulse_n_s','consumed_propellant_kg',
                   'thrust_centroid_time_s','final_supply_pressure_pa')
_BLOWDOWN_INPUT_NAMES={
    'initial_pressure_pa','initial_ullage_volume_m3','initial_propellant_mass_kg',
    'propellant_density_kg_m3','polytropic_exponent','rated_delivered_thrust_n',
    'rated_system_specific_impulse_s','minimum_supply_pressure_pa',
    'thrust_pressure_exponent','specific_impulse_pressure_exponent','chamber_pressure_exponent',
    'reserve_mass_kg','cant_efficiency',
}


def _sensitivity_records(samples,parameters):
    """Rank Pearson screening sensitivities [1]_ of every metric to every input.

    Within each metric, inputs are ranked by absolute correlation; undefined
    coefficients (constant input or output) are listed last with no rank.

    References
    ----------
    .. [1] NIST-PEARSON. NIST Dataplot, centered-product correlation definition.
    """
    result=[]
    for metric in _BLOWDOWN_METRICS:
        output=[getattr(row,metric) for row in samples]
        ranked=[(p.name,_correlation([dict(row.inputs)[p.name] for row in samples],output)) for p in parameters]
        ranked.sort(key=lambda item:abs(item[1]) if item[1] is not None else -1,reverse=True)
        for rank,(name,value) in enumerate(ranked,start=1):
            result.append(PropulsionSensitivity(metric,name,value,None if value is None else abs(value),
                None if value is None else rank,'constant_input_or_output' if value is None else None))
    return tuple(result)


def _validate_uncertainty_result(result):
    """Check that a blowdown result reproduces from its retained samples.

    Samples are reordered by index, and statistics and sensitivities are
    replaced by their recomputed canonical values after the comparison.

    Raises
    ------
    DomainError
        If counts, coverage, bounds, statistics, sensitivities or stability
        evidence do not reproduce.
    """
    if result.method not in {'latin-hypercube','monte-carlo'} or not 4<=result.sample_count<=10000:
        raise DomainError("Unsupported blowdown sampling method or sample count.")
    names=tuple(p.name for p in result.parameters)
    if not names or len(set(names))!=len(names) or not set(names)<=_BLOWDOWN_INPUT_NAMES or len(result.samples)!=result.sample_count:
        raise DomainError("Blowdown result needs unique parameters and complete sample coverage.")
    if result.correlation is not None and (
        not isinstance(result.correlation,CorrelationModel) or result.correlation.parameter_names!=names
    ):
        raise DomainError("Blowdown correlation order does not match the parameters.")
    if result.aleatory_parameter_count!=sum(p.uncertainty_class is UncertaintyClass.ALEATORY for p in result.parameters) or result.epistemic_parameter_count!=sum(p.uncertainty_class is UncertaintyClass.EPISTEMIC for p in result.parameters):
        raise DomainError("Blowdown uncertainty class counts do not reproduce.")
    samples=tuple(sorted(result.samples,key=lambda row:row.index))
    if tuple(row.index for row in samples)!=tuple(range(result.sample_count)):
        raise DomainError("Blowdown sample indices must cover each realization exactly once.")
    for row in samples:
        values=dict(row.inputs)
        if set(values)!=set(names) or any(not p.lower<=values[p.name]<=p.upper for p in result.parameters):
            raise DomainError("Blowdown sample inputs do not match their declared bounded parameters.")
    expected_statistics=tuple(_statistics(metric,[getattr(row,metric) for row in samples]) for metric in _BLOWDOWN_METRICS)
    actual_statistics={row.metric:row for row in result.statistics}
    if len(actual_statistics)!=len(result.statistics) or actual_statistics!={row.metric:row for row in expected_statistics}:
        raise DomainError("Blowdown statistics do not reproduce the retained samples.")
    expected_sensitivity=_sensitivity_records(samples,result.parameters)
    actual_sensitivity={(row.metric,row.parameter):row for row in result.sensitivity}
    if len(actual_sensitivity)!=len(result.sensitivity) or actual_sensitivity!={(row.metric,row.parameter):row for row in expected_sensitivity}:
        raise DomainError("Blowdown sensitivities do not reproduce the retained samples.")
    if not isinstance(result.convergence,SamplingConvergenceEvidence) or result.convergence!=_sampling_convergence(
            [row.delivered_impulse_n_s for row in samples],result.convergence.requested_tolerance):
        raise DomainError("Blowdown stability evidence does not reproduce the retained samples.")
    object.__setattr__(result,'samples',samples)
    object.__setattr__(result,'statistics',expected_statistics)
    object.__setattr__(result,'sensitivity',expected_sensitivity)


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

    Every realization recomputes stream mass flows from delivered thrust and
    system specific impulse, runs the deterministic pressure-fed blowdown model
    with its convergence gate [1]_ [2]_, and retains the complete tabulated
    artifact. An invalid realization aborts the study with its sample index; it
    is never silently discarded.

    Parameters
    ----------
    tank : BlowdownTankModel
        Nominal tank state.
    rated_point : EngineOperatingPoint
        Nominal rated operating point.
    pressure_law : PressurePerformanceLaw
        Nominal supply-pressure performance law.
    parameters : tuple of UncertainParameter
        Uncertain inputs. Supported names are ``initial_pressure_pa``,
        ``initial_ullage_volume_m3``, ``initial_propellant_mass_kg``,
        ``propellant_density_kg_m3``, ``polytropic_exponent``,
        ``rated_delivered_thrust_n``, ``rated_system_specific_impulse_s``,
        ``minimum_supply_pressure_pa``, ``thrust_pressure_exponent``,
        ``specific_impulse_pressure_exponent``, ``chamber_pressure_exponent``,
        ``reserve_mass_kg`` and ``cant_efficiency``. Each ``nominal`` must equal
        the corresponding model input.
    reserve_mass_kg : float, optional
        Nominal unusable reserve [kg], nonnegative.
    cant_efficiency : float, optional
        Nominal axial cant efficiency in ``(0, 1]``.
    sample_count : int, optional
        Number of realizations, 4 to 10000.
    seed : int, optional
        Sampling seed.
    method : {"latin-hypercube", "monte-carlo"}, optional
        Sampling method [3]_.
    correlation : CorrelationModel or None, optional
        Gaussian-copula model in parameter order.
    convergence_tolerance : float, optional
        Positive tolerance of the half/full stability diagnostic.
    maximum_mass_step_kg : float, optional
        Positive maximum propellant-mass integration step [kg].
    integration_relative_tolerance : float, optional
        Positive relative tolerance of the blowdown integration.

    Returns
    -------
    BlowdownUncertaintyResult
        Retained realizations with statistics, Pearson screening
        sensitivities and stability evidence.

    Raises
    ------
    DomainError
        If a record type, control, parameter name, nominal or bound is invalid;
        if the uncertainty bounds can exceed the pressure-law validity, remove
        the usable pressure range or consume the whole inventory; or if any
        realization fails.

    Notes
    -----
    The reported mean interval is an IID normal-formula proxy without
    established 95 % coverage; Latin-hypercube rows are dependent, so it is not
    a design-specific confidence interval. Pearson sensitivity is a screening
    measure. Statistics are engineering uncertainty estimates, not flight
    certification [4]_ [5]_.

    References
    ----------
    .. [1] NASA SP-8112, "Pressurization Systems for Liquid Rockets".
           https://ntrs.nasa.gov/citations/19760015212
    .. [2] NASA CR-131400, "Reliability Model of a Monopropellant Auxiliary
           Propulsion System". https://ntrs.nasa.gov/citations/19730012094
    .. [3] MCKAY-1979. M. D. McKay, R. J. Beckman and W. J. Conover, "A Comparison
           of Three Methods for Selecting Values of Input Variables in the Analysis
           of Output from a Computer Code", Technometrics 21(2), 239-245, 1979.
    .. [4] NASA-STD-7009B, "Standard for Models and Simulations".
           https://standards.nasa.gov/standard/NASA/NASA-STD-7009
    .. [5] NASA/SP-2011-3421, "Probabilistic Risk Assessment Procedures Guide for
           NASA Managers and Practitioners". https://ntrs.nasa.gov/citations/20120001369
    """

    if not isinstance(tank,BlowdownTankModel) or not isinstance(rated_point,EngineOperatingPoint) or not isinstance(pressure_law,PressurePerformanceLaw):
        raise DomainError("Blowdown uncertainty requires declared tank, operating point and pressure-law records.")
    parameters=records(parameters,UncertainParameter,"Blowdown uncertainty parameters")
    sample_count=integer(sample_count,"Blowdown sample count")
    seed=integer(seed,"Blowdown seed")
    method=text(method,"Blowdown sampling method").strip().lower()
    if method not in {'latin-hypercube','monte-carlo'} or not 4<=sample_count<=10000:
        raise DomainError("Unsupported blowdown sampling method or sample count.")
    convergence_tolerance=real(convergence_tolerance,"Uncertainty convergence tolerance")
    reserve_mass_kg=real(reserve_mass_kg,"Uncertainty reserve mass")
    cant_efficiency=real(cant_efficiency,"Uncertainty cant efficiency")
    maximum_mass_step_kg=real(maximum_mass_step_kg,"Uncertainty maximum mass step")
    integration_relative_tolerance=real(integration_relative_tolerance,"Uncertainty integration tolerance")
    if reserve_mass_kg<0 or not 0<cant_efficiency<=1 or maximum_mass_step_kg<=0 or integration_relative_tolerance<=0:
        raise DomainError("Blowdown uncertainty controls are outside their physical domains.")
    if convergence_tolerance <= 0.0:
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
        if parameter.nominal != nominal[parameter.name]:
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
    sensitivity = _sensitivity_records(tuple(realizations),parameters)
    convergence = _sampling_convergence(metric_values["delivered_impulse_n_s"],convergence_tolerance)
    warnings = [
        "Ensemble statistics are engineering uncertainty estimates, not flight certification.",
        "Pearson sensitivity is a screening measure and can miss nonlinear interactions.",
        "The reported mean interval is an IID normal-formula proxy; 95% coverage is not established.",
    ]
    if method=='latin-hypercube':
        warnings.append("Latin-hypercube rows are dependent; the IID-formula interval is not a design-specific confidence interval.")
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
