"""Hierarchical discrete-scenario and conditional-uncertainty ensembles.

Probability-labelled propulsion scenarios (for example nominal operation,
engine-out or degraded modes [1]_ [2]_) are composed with sourced continuous
uncertainty conditional on each scenario [3]_. Scenario probability is applied
exactly once, through the joint weight of each realization. Every retained
realization is a complete thrust and mass-flow artifact that a downstream
trajectory propagator (for example Sidera) can consume directly.

Ensemble statistics use the joint weights and report the law-of-total-variance
split into within- and between-scenario parts [4]_. They are engineering
uncertainty estimates, not flight certification [5]_.

References
----------
.. [1] NASA, "Key Reliability Drivers of Liquid Propulsion Engines and a
       Reliability Model for Sensitivity Analysis".
       https://ntrs.nasa.gov/citations/20050207429
.. [2] NASA, "Common Cause Failure Modeling in Space Launch Vehicles".
       https://ntrs.nasa.gov/citations/20160007073
.. [3] NASA/SP-2011-3421, "Probabilistic Risk Assessment Procedures Guide for
       NASA Managers and Practitioners". https://ntrs.nasa.gov/citations/20120001369
.. [4] MIT-TOTAL-VARIANCE. MIT RES.6-012, lecture 13.7, conditional and
       between-group variance identity (law of total variance).
.. [5] NASA-STD-7009B, "Standard for Models and Simulations".
       https://standards.nasa.gov/standard/NASA/NASA-STD-7009
"""

from __future__ import annotations

import hashlib
from collections.abc import Callable
from dataclasses import dataclass
from fractions import Fraction
from math import isfinite

from rocket_propulsion.core.errors import DomainError
from rocket_propulsion.core.validation import (
    coerce_fields,
    integer,
    name_value_pairs,
    real,
    records,
    strings,
    text,
)

from .marginal_numerics import _SAMPLING_POLICY
from .numeric_outputs import _finite_rational, _fraction
from .probability_numerics import (
    _joint_weight,
    _quantiles,
    _relative_change,
    _sample_moments,
    _sqrt_output,
    _weight_ledger,
)
from .scenarios import PropulsionScenarioEnsemble
from .tabulated import TabulatedBurnArtifact, scale_tabulated_artifact
from .uncertainty import (
    UNCERTAINTY_REFERENCE_IDS,
    CorrelationModel,
    UncertainParameter,
    UncertaintyDraw,
    sample_uncertain_parameters,
)

HIERARCHICAL_ENSEMBLE_REFERENCE_IDS = tuple(
    dict.fromkeys(
        UNCERTAINTY_REFERENCE_IDS
        + (
            "NASA-20050207429",
            "NASA-20160007073",
            "NASA-TM-107318",
        )
    )
)

_ALLOWED_SCALE_NAMES = {
    "thrust_scale",
    "specific_impulse_scale",
    "time_scale",
    "axial_efficiency_scale",
}
_SCALE_COMPONENT_SEPARATOR = "::"


def _target_scale_name(parameter_name: str) -> str:
    """Return the physical scale of a plain or qualified component name."""
    return parameter_name.split(_SCALE_COMPONENT_SEPARATOR, 1)[0]


def _validate_scale_component_name(parameter: UncertainParameter) -> None:
    """Check that a qualified name encodes the parameter's class and a component.

    Plain names without ``"::"`` are accepted unchanged.

    Raises
    ------
    DomainError
        If a qualified name does not have the form
        ``"<scale>::<class>:<component>"`` with the parameter's own class.
    """
    if _SCALE_COMPONENT_SEPARATOR not in parameter.name:
        return
    _target, qualifier = parameter.name.split(_SCALE_COMPONENT_SEPARATOR, 1)
    declared_class, separator, component_id = qualifier.partition(":")
    if (
        not separator
        or not component_id.strip()
        or declared_class != parameter.uncertainty_class.value
    ):
        raise DomainError(
            "Qualified scale names must encode their uncertainty class and component ID."
        )


@dataclass(frozen=True, slots=True)
class ScenarioConditionalUncertainty:
    """Continuous uncertainty definition conditional on one scenario.

    The conditional distribution is sampled by Latin hypercube or Monte Carlo
    [2]_.

    Parameters
    ----------
    scenario_id : str
        Nonempty scenario identifier.
    parameters : tuple of UncertainParameter
        Multiplicative artifact scales (``thrust_scale``,
        ``specific_impulse_scale``, ``time_scale``, ``axial_efficiency_scale``
        or qualified components of them) with unique names, nominal one and
        positive lower bounds. An empty tuple defines a deterministic scenario.
    sample_count : int
        One for a deterministic scenario; 4 to 10000 otherwise.
    source, rationale : str
        Nonempty provenance and the reason the conditional distribution applies
        to this scenario [1]_.
    correlation : CorrelationModel or None, optional
        Gaussian-copula model in parameter order; not allowed for a
        deterministic scenario.

    Raises
    ------
    DomainError
        If any constraint above is violated.

    References
    ----------
    .. [1] NASA-STD-7009B, "Standard for Models and Simulations".
           https://standards.nasa.gov/standard/NASA/NASA-STD-7009
    .. [2] NASA/SP-2011-3421, "Probabilistic Risk Assessment Procedures Guide for
           NASA Managers and Practitioners". https://ntrs.nasa.gov/citations/20120001369
    """

    scenario_id: str
    parameters: tuple[UncertainParameter, ...]
    sample_count: int
    source: str
    rationale: str
    correlation: CorrelationModel | None = None

    def __post_init__(self) -> None:
        coerce_fields(self, "Conditional definition")
        object.__setattr__(self, "parameters", records(self.parameters, UncertainParameter, "Conditional parameters"))
        if self.correlation is not None and not isinstance(self.correlation, CorrelationModel):
            raise DomainError("Conditional correlation requires a CorrelationModel.")
        if not self.scenario_id.strip():
            raise DomainError("Conditional uncertainty scenario ID must not be empty.")
        if not self.source.strip() or not self.rationale.strip():
            raise DomainError("Conditional uncertainty source and rationale are required.")
        if not self.parameters:
            if self.sample_count != 1 or self.correlation is not None:
                raise DomainError(
                    "A deterministic conditional scenario requires sample_count=1 and no correlation."
                )
            return
        if not 4 <= self.sample_count <= 10_000:
            raise DomainError("Conditional uncertainty sample count must be between 4 and 10000.")
        names = tuple(parameter.name for parameter in self.parameters)
        if len(names) != len(set(names)):
            raise DomainError("Conditional uncertainty parameter names must be unique.")
        unknown = {name for name in names if _target_scale_name(name) not in _ALLOWED_SCALE_NAMES}
        if unknown:
            raise DomainError(
                f"Unsupported conditional artifact scales: {', '.join(sorted(unknown))}."
            )
        for parameter in self.parameters:
            _validate_scale_component_name(parameter)
            if abs(parameter.nominal - 1.0) > 1e-12:
                raise DomainError("Conditional artifact scale nominals must equal one.")
            if parameter.lower <= 0.0:
                raise DomainError("Conditional artifact scale bounds must remain positive.")
        if self.correlation is not None and self.correlation.parameter_names != names:
            raise DomainError(
                "Conditional correlation order must match the declared parameter order."
            )


@dataclass(frozen=True, slots=True)
class ConditionalPropulsionRealization:
    """One joint scenario and continuous realization with explicit weights.

    Parameters
    ----------
    scenario_id : str
        Scenario identifier.
    scenario_probability : float
        Scenario probability ``P_s`` in ``(0, 1]``.
    conditional_index : int
        Draw index within the scenario.
    conditional_probability : float
        Conditional weight ``1 / N_s`` in ``(0, 1]``.
    joint_probability : float
        ``P_s / N_s``; must equal the correctly rounded product of the two
        probabilities [1]_.
    inputs : tuple of (str, float)
        Sampled scale inputs.
    artifact : TabulatedBurnArtifact
        Perturbed propulsion history.

    Raises
    ------
    DomainError
        If the artifact type, probabilities or joint product are invalid.

    References
    ----------
    .. [1] NASA/SP-2011-3421, "Probabilistic Risk Assessment Procedures Guide for
           NASA Managers and Practitioners". https://ntrs.nasa.gov/citations/20120001369
    """

    scenario_id: str
    scenario_probability: float
    conditional_index: int
    conditional_probability: float
    joint_probability: float
    inputs: tuple[tuple[str, float], ...]
    artifact: TabulatedBurnArtifact

    def __post_init__(self) -> None:
        coerce_fields(self, "Conditional realization")
        object.__setattr__(self, "inputs", name_value_pairs(self.inputs, "Conditional inputs"))
        if not isinstance(self.artifact, TabulatedBurnArtifact):
            raise DomainError("Conditional realization requires a TabulatedBurnArtifact.")
        if not 0 < self.scenario_probability <= 1 or not 0 < self.conditional_probability <= 1:
            raise DomainError("Retained realization probabilities must lie in (0, 1].")
        if self.joint_probability != _joint_weight(self.scenario_probability, self.conditional_probability):
            raise DomainError("Joint probability does not match its checked binary product.", field="joint_probability")


@dataclass(frozen=True, slots=True)
class HierarchicalMetricStatistics:
    """Joint-probability-weighted statistics of one propulsion metric.

    The statistics summarize probabilistic propagation output [2]_.

    Parameters
    ----------
    metric : str
        Metric name.
    expected_value : float
        Normalized weighted mean [1]_.
    standard_deviation : float
        Square root of the normalized weighted variance.
    minimum, maximum : float
        Smallest and largest retained values.
    p05, p50, p95 : float
        Weighted percentiles: the smallest value whose cumulative joint weight
        reaches the probability.

    References
    ----------
    .. [1] NIST-WEIGHTED-MEAN. NIST Dataplot, "Weighted Mean", equation 2-20,
           normalized weighted mean.
    .. [2] NASA/SP-2011-3421, "Probabilistic Risk Assessment Procedures Guide for
           NASA Managers and Practitioners". https://ntrs.nasa.gov/citations/20120001369
    """

    metric: str
    expected_value: float
    standard_deviation: float
    minimum: float
    p05: float
    p50: float
    p95: float
    maximum: float

    def __post_init__(self) -> None:
        coerce_fields(self, "Hierarchical statistic")


@dataclass(frozen=True, slots=True)
class ScenarioMetricContribution:
    """Conditional mean and joint expected-value contribution of one scenario.

    Parameters
    ----------
    scenario_id, metric : str
        Scenario and metric.
    scenario_probability : float
        Scenario probability ``P_s`` [2]_.
    conditional_expected_value : float
        Weighted mean of the metric within the scenario, ``E[X | s]``.
    joint_expected_value_contribution : float
        ``P_s E[X | s]`` normalized by the total joint weight; the
        contributions of all scenarios sum to the ensemble mean [1]_.

    References
    ----------
    .. [1] MIT-TOTAL-VARIANCE. MIT RES.6-012, lecture 13.7, conditional and
           between-group variance identity (law of total variance).
    .. [2] NASA/SP-2011-3421, "Probabilistic Risk Assessment Procedures Guide for
           NASA Managers and Practitioners". https://ntrs.nasa.gov/citations/20120001369
    """

    scenario_id: str
    metric: str
    scenario_probability: float
    conditional_expected_value: float
    joint_expected_value_contribution: float

    def __post_init__(self) -> None:
        coerce_fields(self, "Scenario contribution")


@dataclass(frozen=True, slots=True)
class HierarchicalVarianceDecomposition:
    """Law-of-total-variance decomposition of one propulsion metric.

    ``Var(X) = E[Var(X | s)] + Var(E[X | s])`` [1]_. The split attributes the
    output variance to the discrete scenarios and to the continuous
    uncertainty within them, a sensitivity measure for the ensemble [2]_ [3]_.

    Parameters
    ----------
    metric : str
        Metric name.
    total_variance : float
        Joint-weighted variance ``Var(X)``.
    within_scenario_variance : float
        Probability-weighted conditional variance ``E[Var(X | s)]``.
    between_scenario_variance : float
        Variance of the conditional scenario means ``Var(E[X | s])``.
    closure_error : float
        ``|total - within - between|`` of the three exported floats; the exact
        rational parts close identically.
    within_fraction, between_fraction : float
        Shares of the total variance; both zero when it is zero.

    References
    ----------
    .. [1] MIT-TOTAL-VARIANCE. MIT RES.6-012, lecture 13.7, conditional and
           between-group variance identity (law of total variance).
    .. [2] NASA/SP-2011-3421, "Probabilistic Risk Assessment Procedures Guide for
           NASA Managers and Practitioners". https://ntrs.nasa.gov/citations/20120001369
    .. [3] NASA-STD-7009B, "Standard for Models and Simulations".
           https://standards.nasa.gov/standard/NASA/NASA-STD-7009
    """

    metric: str
    total_variance: float
    within_scenario_variance: float
    between_scenario_variance: float
    closure_error: float
    within_fraction: float
    between_fraction: float

    def __post_init__(self) -> None:
        coerce_fields(self, "Variance decomposition")


@dataclass(frozen=True, slots=True)
class ConditionalSamplingEvidence:
    """Independent-replicate stability evidence for one conditional ensemble.

    Primary and audit samples use separate, deterministically derived seeds but
    the same distributions and sample count. The relative differences of their
    delivered-impulse mean and standard deviation are a diagnostic of sample
    adequacy [1]_ [2]_, not a proof that tail probabilities have converged.

    Persisted diagnostics are supplied records: the reader checks identities,
    seed links and flag consistency, but does not authenticate a campaign or
    rerun the primary and audit physical histories.

    Parameters
    ----------
    scenario_id : str
        Scenario identifier.
    sample_count : int
        Conditional sample count of both replicates.
    primary_seed, audit_seed : int
        Derived seeds of the two replicates.
    impulse_mean_relative_difference : float
        ``|m_p - m_a| / max(|m_p|, |m_a|)`` of delivered impulse.
    impulse_standard_deviation_relative_difference : float
        The same relative difference of the ``n - 1`` sample standard
        deviations [3]_.
    requested_tolerance : float
        Acceptance tolerance.
    converged : bool
        ``True`` when both differences are within the tolerance.

    Raises
    ------
    DomainError
        If ``converged`` is not a boolean.

    References
    ----------
    .. [1] NASA-STD-7009B, "Standard for Models and Simulations".
           https://standards.nasa.gov/standard/NASA/NASA-STD-7009
    .. [2] NASA/SP-2011-3421, "Probabilistic Risk Assessment Procedures Guide for
           NASA Managers and Practitioners". https://ntrs.nasa.gov/citations/20120001369
    .. [3] NIST TN 1297, "Guidelines for Evaluating and Expressing the Uncertainty
           of NIST Measurement Results", 1994.
           https://www.nist.gov/pml/nist-technical-note-1297
    """

    scenario_id: str
    sample_count: int
    primary_seed: int
    audit_seed: int
    impulse_mean_relative_difference: float
    impulse_standard_deviation_relative_difference: float
    requested_tolerance: float
    converged: bool

    def __post_init__(self) -> None:
        coerce_fields(self, "Conditional sampling evidence")
        if type(self.converged) is not bool:
            raise DomainError("Conditional convergence must be boolean.", field="converged")


@dataclass(frozen=True, slots=True)
class HierarchicalPropulsionEnsemble:
    """Scenario-weighted conditional propulsion ensemble for trajectory propagation.

    Produced by :func:`build_hierarchical_propulsion_ensemble` and consumed by
    downstream propagators such as Sidera. The statistics are engineering
    estimates, not flight certification [4]_ [5]_.

    Parameters
    ----------
    method : str
        Conditional sampling method [1]_.
    base_seed : int
        Seed from which every scenario seed is derived.
    scenario_count, realization_count : int
        Number of scenarios and retained realizations.
    definitions : tuple of ScenarioConditionalUncertainty
        One conditional definition per scenario.
    realizations : tuple of ConditionalPropulsionRealization
        Retained joint realizations.
    statistics : tuple of HierarchicalMetricStatistics
        Joint-weighted metric statistics.
    variance_decomposition : tuple of HierarchicalVarianceDecomposition
        Law-of-total-variance split per metric [2]_.
    contributions : tuple of ScenarioMetricContribution
        Per-scenario mean contributions.
    sampling_evidence : tuple of ConditionalSamplingEvidence
        Replicate stability per scenario.
    zero_thrust_probability : float
        Normalized joint weight of realizations with zero delivered impulse,
        for example engine-out scenarios [3]_.
    joint_probability_sum, joint_probability_closure_error : float
        Sum of the joint weights and its exact distance from one.
    reference_ids, warnings : tuple of str
        Source identifiers and interpretation warnings.

    References
    ----------
    .. [1] MCKAY-1979. M. D. McKay, R. J. Beckman and W. J. Conover, "A Comparison
           of Three Methods for Selecting Values of Input Variables in the Analysis
           of Output from a Computer Code", Technometrics 21(2), 239-245, 1979.
    .. [2] MIT-TOTAL-VARIANCE. MIT RES.6-012, lecture 13.7, conditional and
           between-group variance identity (law of total variance).
    .. [3] NASA, "Key Reliability Drivers of Liquid Propulsion Engines and a
           Reliability Model for Sensitivity Analysis".
           https://ntrs.nasa.gov/citations/20050207429
    .. [4] NASA-STD-7009B, "Standard for Models and Simulations".
           https://standards.nasa.gov/standard/NASA/NASA-STD-7009
    .. [5] NASA/SP-2011-3421, "Probabilistic Risk Assessment Procedures Guide for
           NASA Managers and Practitioners". https://ntrs.nasa.gov/citations/20120001369
    """

    method: str
    base_seed: int
    scenario_count: int
    realization_count: int
    definitions: tuple[ScenarioConditionalUncertainty, ...]
    realizations: tuple[ConditionalPropulsionRealization, ...]
    statistics: tuple[HierarchicalMetricStatistics, ...]
    variance_decomposition: tuple[HierarchicalVarianceDecomposition, ...]
    contributions: tuple[ScenarioMetricContribution, ...]
    sampling_evidence: tuple[ConditionalSamplingEvidence, ...]
    zero_thrust_probability: float
    joint_probability_sum: float
    joint_probability_closure_error: float
    reference_ids: tuple[str, ...]
    warnings: tuple[str, ...]

    def __post_init__(self) -> None:
        coerce_fields(self, "Hierarchical ensemble")
        for name, kind in (("definitions", ScenarioConditionalUncertainty),
                           ("realizations", ConditionalPropulsionRealization),
                           ("statistics", HierarchicalMetricStatistics),
                           ("variance_decomposition", HierarchicalVarianceDecomposition),
                           ("contributions", ScenarioMetricContribution),
                           ("sampling_evidence", ConditionalSamplingEvidence)):
            object.__setattr__(self, name, records(getattr(self, name), kind, name))
        object.__setattr__(self, "reference_ids", strings(self.reference_ids, "Ensemble references"))
        object.__setattr__(self, "warnings", strings(self.warnings, "Ensemble warnings"))


def _derived_seed(base_seed: int, scenario_id: str, channel: str) -> int:
    """Derive a 64-bit seed from ``"<base_seed>:<scenario_id>:<channel>"``.

    The first eight bytes of the SHA-256 digest, read big-endian, give a
    reproducible seed per scenario and channel (``"primary"`` or ``"audit"``).
    """
    digest = hashlib.sha256(f"{base_seed}:{scenario_id}:{channel}".encode()).digest()
    return int.from_bytes(digest[:8], "big", signed=False)


def _draws(
    definition: ScenarioConditionalUncertainty,
    *,
    seed: int,
    method: str,
) -> tuple[UncertaintyDraw, ...]:
    """Sample a conditional definition; a deterministic scenario gives one empty draw."""
    if not definition.parameters:
        return (UncertaintyDraw(index=0, values=()),)
    return sample_uncertain_parameters(
        definition.parameters,
        sample_count=definition.sample_count,
        seed=seed,
        method=method,
        correlation=definition.correlation,
    )


def _scale_ledger(draw: UncertaintyDraw):
    """Return the exact product of all sampled factors acting on each scale.

    Unsampled scales are one; qualified components multiply their target scale.
    """
    values = {
        "thrust_scale": Fraction(1),
        "specific_impulse_scale": Fraction(1),
        "time_scale": Fraction(1),
        "axial_efficiency_scale": Fraction(1),
    }
    for component_name, value in draw.values:
        values[_target_scale_name(component_name)] *= _fraction(value)
    return values


def _scale_values(draw: UncertaintyDraw) -> dict[str, float]:
    """Return the scale products of ``draw`` rounded once to binary64."""
    return {name: _finite_rational(value, f"Artifact scale {name}")
            for name, value in _scale_ledger(draw).items()}


def _validate_axial_scale(
    artifact: TabulatedBurnArtifact, definition: ScenarioConditionalUncertainty
) -> None:
    """Check that the largest axial-efficiency factor keeps axial <= delivered thrust.

    The worst case multiplies the largest axial/delivered thrust ratio of the
    artifact by the product of every axial-efficiency upper bound.

    Raises
    ------
    DomainError
        If that product can exceed one (beyond a ``1e-12`` tolerance).
    """
    parameters = tuple(
        item
        for item in definition.parameters
        if _target_scale_name(item.name) == "axial_efficiency_scale"
    )
    if not parameters:
        return
    maximum_ratio = max(
        (
            axial / delivered
            for segment in artifact.segments
            for axial, delivered in (
                (segment.start_axial_thrust_n, segment.start_delivered_thrust_n),
                (segment.end_axial_thrust_n, segment.end_delivered_thrust_n),
            )
            if delivered > 0.0
        ),
        default=0.0,
    )
    maximum_scale = 1.0
    for parameter in parameters:
        maximum_scale *= parameter.upper
    if maximum_ratio * maximum_scale > 1.0 + 1e-12:
        raise DomainError(
            "Axial-efficiency uncertainty can make axial thrust exceed delivered thrust."
        )


def _perturb_artifact(
    artifact: TabulatedBurnArtifact,
    *,
    scenario_id: str,
    definition: ScenarioConditionalUncertainty,
    draw: UncertaintyDraw,
    primary_seed: int,
    method: str,
) -> TabulatedBurnArtifact:
    """Scale a scenario artifact by one draw and bind the draw as provenance.

    The scenario probability is not stored in the artifact; it is applied once
    through the realization's joint weight.
    """
    values = _scale_values(draw)
    return scale_tabulated_artifact(
        artifact,
        values,
        source_id=f"conditional:{scenario_id}:{draw.index}",
        provenance={
            "scenario_id": scenario_id,
            "conditional_definition": definition,
            "conditional_index": draw.index,
            "inputs": draw.values,
            "primary_seed": primary_seed,
            "method": method,
            "sampling_numerics": _SAMPLING_POLICY,
        },
        reference_ids=tuple(dict.fromkeys(HIERARCHICAL_ENSEMBLE_REFERENCE_IDS)),
        warnings=(
            "Artifact contains a conditional continuous uncertainty realization.",
            "Scenario probability is stored outside the artifact and must be applied once.",
        ),
    )


def _weighted_percentile(values: tuple[tuple[float, float], ...], probability: float) -> float:
    """Return the smallest value whose cumulative weight reaches ``probability``."""
    return _weight_ledger(values).quantile(probability)


def _weighted_statistics(
    metric: str, values: tuple[tuple[float, float], ...]
) -> HierarchicalMetricStatistics:
    """Return joint-weighted statistics from ``(value, weight)`` pairs.

    The mean is the normalized weighted mean [1]_ and the variance the
    normalized weighted second central moment; each export is rounded once.

    References
    ----------
    .. [1] NIST-WEIGHTED-MEAN. NIST Dataplot, "Weighted Mean", equation 2-20,
           normalized weighted mean.
    """
    ledger = _weight_ledger(values)
    p05, p50, p95 = _quantiles(ledger)
    return HierarchicalMetricStatistics(
        metric=metric,
        expected_value=_finite_rational(ledger.mean, "Weighted mean"),
        standard_deviation=_sqrt_output(ledger.variance, "Weighted standard deviation"),
        minimum=_finite_rational(ledger.rows[0][0], "Weighted minimum"),
        p05=p05, p50=p50, p95=p95,
        maximum=_finite_rational(ledger.rows[-1][0], "Weighted maximum"),
    )


def _variance_decomposition(
    metric: str,
    realizations: tuple[ConditionalPropulsionRealization, ...],
    getter: Callable[[TabulatedBurnArtifact], float],
) -> HierarchicalVarianceDecomposition:
    """Split the joint-weighted variance of a metric by scenario.

    With total weight ``W``, scenario weight ``W_s``, scenario mean ``m_s``,
    scenario variance ``v_s`` and overall mean ``m``, the within part is
    ``sum_s W_s v_s / W`` and the between part ``sum_s W_s (m_s - m)**2 / W``;
    their exact sum equals the total variance [1]_.

    References
    ----------
    .. [1] MIT-TOTAL-VARIANCE. MIT RES.6-012, lecture 13.7, conditional and
           between-group variance identity (law of total variance).
    """
    ledger = _weight_ledger((getter(item.artifact), item.joint_probability) for item in realizations)
    scenario_ids = tuple(dict.fromkeys(item.scenario_id for item in realizations))
    within_variance = Fraction(0)
    between_variance = Fraction(0)
    for scenario_id in scenario_ids:
        selected = tuple(item for item in realizations if item.scenario_id == scenario_id)
        group = _weight_ledger((getter(item.artifact), item.joint_probability) for item in selected)
        within_variance += group.total * group.variance / ledger.total
        between_variance += group.total * (group.mean - ledger.mean)**2 / ledger.total
    total_output = _finite_rational(ledger.variance, "Total variance")
    within_output = _finite_rational(within_variance, "Within-scenario variance")
    between_output = _finite_rational(between_variance, "Between-scenario variance")
    closure_error = abs(_fraction(total_output)-_fraction(within_output)-_fraction(between_output))
    if ledger.variance > 0:
        within_fraction = within_variance / ledger.variance
        between_fraction = between_variance / ledger.variance
    else:
        within_fraction = 0.0
        between_fraction = 0.0
    return HierarchicalVarianceDecomposition(
        metric=metric,
        total_variance=total_output,
        within_scenario_variance=within_output,
        between_scenario_variance=between_output,
        closure_error=_finite_rational(closure_error, "Variance display closure"),
        within_fraction=_finite_rational(_fraction(within_fraction), "Within variance fraction"),
        between_fraction=_finite_rational(_fraction(between_fraction), "Between variance fraction"),
    )


def _mean_and_deviation(values: tuple[float, ...]) -> tuple[float, float]:
    """Return the sample mean and ``n - 1`` sample standard deviation."""
    return _sample_moments(values)


def _relative_difference(left: float, right: float) -> float:
    """Return ``|left - right| / max(|left|, |right|)``, or zero if both are zero."""
    return _relative_change(left, right)


def _impulses_for_draws(
    base_impulse_n_s: float, draws: tuple[UncertaintyDraw, ...]
) -> tuple[float, ...]:
    """Return the delivered impulse of each draw without building artifacts.

    Delivered impulse scales with ``thrust_scale * time_scale``; specific
    impulse and axial efficiency do not change it.
    """
    impulses = []
    for draw in draws:
        scales = _scale_ledger(draw)
        impulses.append(_finite_rational(_fraction(base_impulse_n_s)*scales["thrust_scale"]*scales["time_scale"], "Audit impulse"))
    return tuple(impulses)


def build_hierarchical_propulsion_ensemble(
    scenarios: PropulsionScenarioEnsemble,
    definitions: tuple[ScenarioConditionalUncertainty, ...],
    *,
    seed: int = 42,
    method: str = "latin-hypercube",
    replicate_tolerance: float = 0.05,
) -> HierarchicalPropulsionEnsemble:
    """Compose scenario probabilities with conditional continuous uncertainty.

    Parameters
    ----------
    scenarios : PropulsionScenarioEnsemble
        Discrete, probability-labelled scenarios with their artifacts.
    definitions : tuple of ScenarioConditionalUncertainty
        Exactly one conditional definition per scenario.
    seed : int, optional
        Base seed; per-scenario seeds are derived from it.
    method : {"latin-hypercube", "monte-carlo"}, optional
        Conditional sampling method [1]_.
    replicate_tolerance : float, optional
        Positive tolerance of the replicate stability check.

    Returns
    -------
    HierarchicalPropulsionEnsemble
        Joint realizations, weighted statistics, variance decomposition,
        scenario contributions and replicate evidence.

    Raises
    ------
    DomainError
        If inputs are invalid, definitions do not cover the scenarios one to
        one, an axial-efficiency factor can exceed delivered thrust, or the
        joint probabilities do not close to one within ``1e-12``.

    Notes
    -----
    For scenario probability ``P_s`` and ``N_s`` equally weighted conditional
    draws, each realization receives the joint probability ``P_s / N_s``, so the
    scenario probability is applied exactly once [2]_. A second, independently
    seeded ensemble of the same size is used only to check the stability of
    the delivered-impulse mean and standard deviation; its samples are not
    mixed into the primary distribution. Scenario sets typically represent
    engine-out, degraded or common-cause cases [3]_ [4]_. Statistics use the
    law of total variance [5]_. Trajectory, pointing, epoch and attitude
    dispersion are outside this model [6]_.

    References
    ----------
    .. [1] MCKAY-1979. M. D. McKay, R. J. Beckman and W. J. Conover, "A Comparison
           of Three Methods for Selecting Values of Input Variables in the Analysis
           of Output from a Computer Code", Technometrics 21(2), 239-245, 1979.
    .. [2] NASA/SP-2011-3421, "Probabilistic Risk Assessment Procedures Guide for
           NASA Managers and Practitioners". https://ntrs.nasa.gov/citations/20120001369
    .. [3] NASA, "Key Reliability Drivers of Liquid Propulsion Engines and a
           Reliability Model for Sensitivity Analysis".
           https://ntrs.nasa.gov/citations/20050207429
    .. [4] NASA, "Common Cause Failure Modeling in Space Launch Vehicles".
           https://ntrs.nasa.gov/citations/20160007073
    .. [5] MIT-TOTAL-VARIANCE. MIT RES.6-012, lecture 13.7, conditional and
           between-group variance identity (law of total variance).
    .. [6] NASA-STD-7009B, "Standard for Models and Simulations".
           https://standards.nasa.gov/standard/NASA/NASA-STD-7009
    """

    if not isinstance(scenarios, PropulsionScenarioEnsemble):
        raise DomainError("Hierarchical builder requires a PropulsionScenarioEnsemble.")
    definitions = records(definitions, ScenarioConditionalUncertainty, "Conditional definitions")
    seed = integer(seed, "Hierarchical base seed")
    replicate_tolerance = real(replicate_tolerance, "Replicate tolerance")
    normalized_method = text(method, "Hierarchical method").strip().lower()
    if normalized_method not in {"latin-hypercube", "monte-carlo"}:
        raise DomainError("Hierarchical method must be latin-hypercube or monte-carlo.")
    if not isfinite(replicate_tolerance) or replicate_tolerance <= 0.0:
        raise DomainError("Replicate convergence tolerance must be positive and finite.")
    definition_map = {definition.scenario_id: definition for definition in definitions}
    if len(definition_map) != len(definitions):
        raise DomainError("Conditional scenario definitions must be unique.")
    scenario_ids = {item.scenario_id for item in scenarios.realizations}
    if definition_map.keys() != scenario_ids:
        raise DomainError("Every discrete scenario requires exactly one conditional definition.")

    realizations: list[ConditionalPropulsionRealization] = []
    evidence: list[ConditionalSamplingEvidence] = []
    for scenario in scenarios.realizations:
        definition = definition_map[scenario.scenario_id]
        _validate_axial_scale(scenario.artifact, definition)
        primary_seed = _derived_seed(seed, scenario.scenario_id, "primary")
        audit_seed = _derived_seed(seed, scenario.scenario_id, "audit")
        primary_draws = _draws(definition, seed=primary_seed, method=normalized_method)
        audit_draws = _draws(definition, seed=audit_seed, method=normalized_method)
        conditional_probability = 1.0 / len(primary_draws)
        for draw in primary_draws:
            artifact = _perturb_artifact(
                scenario.artifact,
                scenario_id=scenario.scenario_id,
                definition=definition,
                draw=draw,
                primary_seed=primary_seed,
                method=normalized_method,
            )
            realizations.append(
                ConditionalPropulsionRealization(
                    scenario_id=scenario.scenario_id,
                    scenario_probability=scenario.probability,
                    conditional_index=draw.index,
                    conditional_probability=conditional_probability,
                    joint_probability=_joint_weight(scenario.probability, conditional_probability),
                    inputs=draw.values,
                    artifact=artifact,
                )
            )
        primary_impulses = _impulses_for_draws(
            scenario.artifact.delivered_total_impulse_n_s, primary_draws
        )
        audit_impulses = _impulses_for_draws(
            scenario.artifact.delivered_total_impulse_n_s, audit_draws
        )
        primary_mean, primary_deviation = _mean_and_deviation(primary_impulses)
        audit_mean, audit_deviation = _mean_and_deviation(audit_impulses)
        mean_difference = _relative_difference(primary_mean, audit_mean)
        deviation_difference = _relative_difference(primary_deviation, audit_deviation)
        evidence.append(
            ConditionalSamplingEvidence(
                scenario_id=scenario.scenario_id,
                sample_count=len(primary_draws),
                primary_seed=primary_seed,
                audit_seed=audit_seed,
                impulse_mean_relative_difference=mean_difference,
                impulse_standard_deviation_relative_difference=deviation_difference,
                requested_tolerance=replicate_tolerance,
                converged=max(mean_difference, deviation_difference) <= replicate_tolerance,
            )
        )

    total_weight = sum((_fraction(item.joint_probability) for item in realizations), Fraction(0))
    exact_closure = abs(total_weight-1)
    joint_probability_sum = _finite_rational(total_weight, "Joint probability sum")
    closure_error = _finite_rational(exact_closure, "Probability closure error")
    if exact_closure > _fraction(1e-12):
        raise DomainError("Hierarchical joint probabilities do not close to one.")
    metric_getters = {
        "duration_s": lambda artifact: artifact.duration_s,
        "delivered_impulse_n_s": lambda artifact: artifact.delivered_total_impulse_n_s,
        "axial_impulse_n_s": lambda artifact: artifact.axial_total_impulse_n_s,
        "consumed_propellant_kg": lambda artifact: artifact.consumed_propellant_kg,
        "thrust_centroid_time_s": lambda artifact: artifact.thrust_centroid_time_s,
    }
    weighted_values = {
        metric: tuple((getter(item.artifact), item.joint_probability) for item in realizations)
        for metric, getter in metric_getters.items()
    }
    contributions: list[ScenarioMetricContribution] = []
    for scenario in scenarios.realizations:
        selected = tuple(item for item in realizations if item.scenario_id == scenario.scenario_id)
        for metric, getter in metric_getters.items():
            group = _weight_ledger((getter(item.artifact), item.joint_probability) for item in selected)
            conditional_mean = _finite_rational(group.mean, "Conditional mean")
            contributions.append(
                ScenarioMetricContribution(
                    scenario_id=scenario.scenario_id,
                    metric=metric,
                    scenario_probability=scenario.probability,
                    conditional_expected_value=conditional_mean,
                    joint_expected_value_contribution=_finite_rational(group.mean*group.total/total_weight, "Scenario mean contribution"),
                )
            )
    warnings = [
        "Scenario probabilities are applied once; conditional draws are equally weighted within scenario.",
        "Trajectory, pointing, epoch, and attitude dispersion remain Sidera responsibilities.",
    ]
    if any(not item.converged for item in evidence):
        warnings.append(
            "At least one independent replicate check failed; increase its conditional sample count."
        )
    return HierarchicalPropulsionEnsemble(
        method=normalized_method,
        base_seed=seed,
        scenario_count=len(scenarios.realizations),
        realization_count=len(realizations),
        definitions=definitions,
        realizations=tuple(realizations),
        statistics=tuple(
            _weighted_statistics(metric, values) for metric, values in weighted_values.items()
        ),
        variance_decomposition=tuple(
            _variance_decomposition(metric, tuple(realizations), getter)
            for metric, getter in metric_getters.items()
        ),
        contributions=tuple(contributions),
        sampling_evidence=tuple(evidence),
        zero_thrust_probability=_finite_rational(sum((
            _fraction(item.joint_probability)
            for item in realizations
            if item.artifact.delivered_total_impulse_n_s == 0.0
        ), Fraction(0))/total_weight, "Zero-thrust probability"),
        joint_probability_sum=joint_probability_sum,
        joint_probability_closure_error=closure_error,
        reference_ids=HIERARCHICAL_ENSEMBLE_REFERENCE_IDS,
        warnings=tuple(warnings),
    )
