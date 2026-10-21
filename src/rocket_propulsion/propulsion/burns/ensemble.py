"""Hierarchical discrete-scenario and conditional-uncertainty ensembles.

The contract composes probability-labelled propulsion scenarios with sourced
continuous uncertainty without double-counting scenario probability. Every
retained realization remains a complete thrust/mass-flow artifact suitable for
later Sidera propagation.

References
----------
NASA-STD-7009B: https://standards.nasa.gov/standard/NASA/NASA-STD-7009
NASA-SP-2011-3421: https://ntrs.nasa.gov/citations/20120001369
NASA-20050207429: https://ntrs.nasa.gov/citations/20050207429
NASA-20160007073: https://ntrs.nasa.gov/citations/20160007073
"""

from __future__ import annotations

import hashlib
from collections.abc import Callable
from dataclasses import dataclass
from math import isfinite

from rocket_propulsion.core.errors import DomainError
from rocket_propulsion.core.validation import coerce_fields, name_value_pairs, records, strings, integer, real, text

from .scenarios import PropulsionScenarioEnsemble
from fractions import Fraction
from .numeric_outputs import _finite_rational, _fraction
from .probability_numerics import _weight_ledger, _quantiles, _sqrt_output, _sample_moments, _relative_change, _joint_weight
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
    return parameter_name.split(_SCALE_COMPONENT_SEPARATOR, 1)[0]


def _validate_scale_component_name(parameter: UncertainParameter) -> None:
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

    A deterministic scenario uses no parameters and ``sample_count=1``. A
    stochastic definition needs at least four samples. Source and rationale
    describe why the conditional distribution is applicable to this scenario.

    References
    ----------
    NASA-STD-7009B: https://standards.nasa.gov/standard/NASA/NASA-STD-7009
    NASA-SP-2011-3421: https://ntrs.nasa.gov/citations/20120001369
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
    """One joint scenario/continuous realization with explicit weights.

    References
    ----------
    NASA-SP-2011-3421: https://ntrs.nasa.gov/citations/20120001369
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
    """Joint-probability statistics for one propulsion metric.

    References
    ----------
    NASA-SP-2011-3421: https://ntrs.nasa.gov/citations/20120001369
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
    """Conditional mean and joint expected-value contribution by scenario.

    References
    ----------
    NASA-SP-2011-3421: https://ntrs.nasa.gov/citations/20120001369
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
    """Law-of-total-variance decomposition for one propulsion metric.

    ``within_scenario_variance`` is the probability-weighted conditional
    variance. ``between_scenario_variance`` is the variance of conditional
    scenario means. Their sum must reproduce ``total_variance``.

    References
    ----------
    NASA-SP-2011-3421: https://ntrs.nasa.gov/citations/20120001369
    NASA-STD-7009B: https://standards.nasa.gov/standard/NASA/NASA-STD-7009
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

    Primary and audit samples use separate deterministically derived seeds but
    the same distributions and sample count. This is a diagnostic, not a proof
    that tail probabilities have converged.

    Persisted diagnostics are supplied records: the reader checks identities,
    seed links and flag consistency, but does not authenticate a campaign or
    independently rerun the primary/audit physical histories.

    References
    ----------
    NASA-STD-7009B: https://standards.nasa.gov/standard/NASA/NASA-STD-7009
    NASA-SP-2011-3421: https://ntrs.nasa.gov/citations/20120001369
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
    """Scenario-weighted conditional propulsion ensemble for Sidera.

    References
    ----------
    NASA-STD-7009B: https://standards.nasa.gov/standard/NASA/NASA-STD-7009
    NASA-SP-2011-3421: https://ntrs.nasa.gov/citations/20120001369
    NASA-20050207429: https://ntrs.nasa.gov/citations/20050207429
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
    digest = hashlib.sha256(f"{base_seed}:{scenario_id}:{channel}".encode()).digest()
    return int.from_bytes(digest[:8], "big", signed=False)


def _draws(
    definition: ScenarioConditionalUncertainty,
    *,
    seed: int,
    method: str,
) -> tuple[UncertaintyDraw, ...]:
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
    return {name: _finite_rational(value, f"Artifact scale {name}")
            for name, value in _scale_ledger(draw).items()}


def _validate_axial_scale(
    artifact: TabulatedBurnArtifact, definition: ScenarioConditionalUncertainty
) -> None:
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
        },
        reference_ids=tuple(dict.fromkeys(HIERARCHICAL_ENSEMBLE_REFERENCE_IDS)),
        warnings=(
            "Artifact contains a conditional continuous uncertainty realization.",
            "Scenario probability is stored outside the artifact and must be applied once.",
        ),
    )


def _weighted_percentile(values: tuple[tuple[float, float], ...], probability: float) -> float:
    return _weight_ledger(values).quantile(probability)


def _weighted_statistics(
    metric: str, values: tuple[tuple[float, float], ...]
) -> HierarchicalMetricStatistics:
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
    return _sample_moments(values)


def _relative_difference(left: float, right: float) -> float:
    return _relative_change(left, right)


def _impulses_for_draws(
    base_impulse_n_s: float, draws: tuple[UncertaintyDraw, ...]
) -> tuple[float, ...]:
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

    For scenario probability ``P_s`` and ``N_s`` equally weighted conditional
    samples, each retained realization receives joint probability
    ``P_s / N_s``. A second independently seeded ensemble of the same size is
    used only to audit impulse mean/standard-deviation stability; its samples
    are not mixed into the primary distribution.

    References
    ----------
    NASA-STD-7009B: https://standards.nasa.gov/standard/NASA/NASA-STD-7009
    NASA-SP-2011-3421: https://ntrs.nasa.gov/citations/20120001369
    NASA-20050207429: https://ntrs.nasa.gov/citations/20050207429
    NASA-20160007073: https://ntrs.nasa.gov/citations/20160007073
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
