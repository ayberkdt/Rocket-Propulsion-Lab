"""Explicit constraints and fixed-scale weighted trade studies.

References
----------
NASA-SE-6: 6.8.1.2.3-6, weighted matrices, operational scales and decision reports.
NASA-SE-APP: Appendix D, requirement/source/verification traceability.
"""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass
from fractions import Fraction
from math import fsum, isclose, isfinite
from typing import Literal

from rocket_propulsion.core.errors import DomainError
from rocket_propulsion.core.validation import records, text

from .evidence import Citation, CitedScalar, derived, quantity


@dataclass(frozen=True, slots=True)
class Requirement:
    """An analyst-selected numerical screen, never a fabricated course requirement.

    References
    ----------
    NASA-SE-APP: Appendix D, identifier, source and verification approach.
    """

    identifier: str
    description: str
    metric: str
    relation: Literal["at_most", "at_least"]
    limit: CitedScalar
    verification: str = "preliminary numerical analysis"

    def __post_init__(self) -> None:
        """Require an identifiable source-bearing threshold and evaluation method.

        References
        ----------
        NASA-SE-APP: Appendix D, traceable requirement and verification method.
        """
        for name in ("identifier", "description", "metric", "verification"):
            text(getattr(self, name), name)
        if not isinstance(self.relation, str) or self.relation not in {"at_most", "at_least"}:
            raise DomainError("Requirement relation must be at_most or at_least.")
        if not isinstance(self.limit, CitedScalar):
            raise DomainError("Requirement thresholds must retain their source or assumption.")


@dataclass(frozen=True, slots=True)
class RequirementCheck:
    """A screen result with an evidence-bearing signed margin, or explicit unknown.

    References
    ----------
    NASA-SE-APP: Appendix D, verification matrix and analysis evidence.
    """

    requirement: Requirement
    status: Literal["pass", "fail", "unknown"]
    measured: CitedScalar | None
    margin: CitedScalar | None


def check_requirements(metrics: Mapping[str, CitedScalar], requirements: tuple[Requirement, ...]
                       ) -> tuple[RequirementCheck, ...]:
    """Evaluate all constraints, leaving absent data unknown and visible.

    Positive margin is feasible. Passing means only that the provided numbers
    pass this preliminary analysis; tests/hardware compliance are not inferred.

    References
    ----------
    NASA-SE-APP: Appendix D, source-linked requirements verification matrix.
    """
    if not isinstance(metrics, Mapping):
        raise DomainError("Metrics must be a name-to-cited-quantity mapping.")
    entries = records(requirements, Requirement, "requirements")
    if len({x.identifier for x in entries}) != len(entries):
        raise DomainError("Requirement identifiers must be unique.")
    checks = []
    for req in entries:
        value = metrics.get(req.metric)
        if value is None:
            checks.append(RequirementCheck(req, "unknown", None, None))
            continue
        actual = quantity(value, req.limit.unit)
        margin = req.limit.value - actual if req.relation == "at_most" else actual - req.limit.value
        evidence = derived(margin, req.limit.unit, (value, req.limit), (
            Citation("NASA-SE-APP", "Appendix D: explicit numerical requirement check"),
        ))
        checks.append(RequirementCheck(req, "pass" if margin >= 0 else "fail", value, evidence))
    return tuple(checks)


@dataclass(frozen=True, slots=True)
class Criterion:
    """Fixed best/worst utility anchors and a user-sourced weight.

    Anchors define preferences, not an empirical data-extrapolation model.
    They are fixed across alternatives so adding a dominated candidate cannot
    rescale the existing scores. Values outside the anchors are clipped.

    References
    ----------
    NASA-SE-6: 6.8.1.2.3-4, weighted matrices and operational normalization scales.
    """

    metric: str
    best: CitedScalar
    worst: CitedScalar
    weight: CitedScalar

    def __post_init__(self) -> None:
        """Reject equal anchors or unsupported weight values.

        References
        ----------
        NASA-SE-6: 6.8.1.2.4, repeatable criteria scales and explicit weights.
        """
        text(self.metric, "trade metric")
        if not isinstance(self.best, CitedScalar):
            raise DomainError("Trade utility anchors must retain evidence.")
        quantity(self.worst, self.best.unit)
        if self.best.value == self.worst.value:
            raise DomainError("Best and worst anchors must differ.")
        if not isfinite(self.best.value - self.worst.value):
            raise DomainError("Utility anchor separation exceeds the finite numerical domain.")
        if quantity(self.weight, "1", minimum=0) > 1:
            raise DomainError("Criterion weight must lie in [0,1].")


@dataclass(frozen=True, slots=True)
class AlternativeScore:
    """Complete alternatives have a score; failed or incomplete ones remain visible.

    References
    ----------
    NASA-SE-6: 6.8.1.2.4-6, constraints before recommendation and retain alternatives.
    """

    name: str
    eligibility: Literal["eligible", "infeasible", "incomplete"]
    score: CitedScalar | None
    utilities: tuple[tuple[str, CitedScalar], ...]
    checks: tuple[RequirementCheck, ...]
    missing_metrics: tuple[str, ...]


def weighted_trade(alternatives: Mapping[str, Mapping[str, CitedScalar]],
                    criteria: tuple[Criterion, ...], requirements: tuple[Requirement, ...]
                    ) -> tuple[AlternativeScore, ...]:
    """Compute sum(w*clip((x-worst)/(best-worst),0,1)) after hard screens.

    This explicit affine preference rule is this implementation's choice of
    the cited normalization method; NASA does not prescribe these anchors or
    weights. Missing power, mass, cost or heritage data never becomes zero.
    Returned order is input order; the caller can inspect ties without a
    hidden tie-break or a claim of a universally best propulsion technology.
    Affine utilities and weighted terms use exact rational input arithmetic;
    displayed utilities and the final sum are rounded separately. The final
    sum never consumes already rounded utility/contribution display values.

    References
    ----------
    NASA-SE-6: 6.8.1.2.3-4, weighted matrices and repeatable normalization.
    NASA-SE-APP: Trade Study Report, state selection rule and all alternatives.
    """
    if not isinstance(alternatives, Mapping) or not alternatives:
        raise DomainError("Trade study needs a non-empty mapping of alternatives.")
    entries = records(criteria, Criterion, "trade criteria")
    if not entries or len({x.metric for x in entries}) != len(entries):
        raise DomainError("Trade criteria must be non-empty and unique.")
    if not isclose(fsum(x.weight.value for x in entries), 1.0, rel_tol=0, abs_tol=1e-12):
        raise DomainError("Trade weights must sum to one; no silent normalization.")
    results = []
    for name, metrics in alternatives.items():
        text(name, "alternative name")
        checks = check_requirements(metrics, requirements)
        missing = tuple(x.metric for x in entries if metrics.get(x.metric) is None)
        failed = any(x.status == "fail" for x in checks)
        incomplete = bool(missing) or any(x.status == "unknown" for x in checks)
        eligibility = "infeasible" if failed else "incomplete" if incomplete else "eligible"
        utilities = []
        terms = []
        score_inputs = []
        for criterion in entries:
            if criterion.metric in missing:
                continue
            value = metrics[criterion.metric]
            actual = quantity(value, criterion.best.unit)
            best, worst, actual = map(Fraction, (criterion.best.value, criterion.worst.value, actual))
            utility = max(Fraction(0), min(Fraction(1), (actual - worst) / (best - worst)))
            q = derived(float(utility), "1", (value, criterion.best, criterion.worst, criterion.weight), (
                Citation("NASA-SE-6", "6.8.1.2.4: explicitly chosen clipped affine utility"),
            ))
            utilities.append((criterion.metric, q))
            terms.append(Fraction(criterion.weight.value) * utility)
            score_inputs.extend((value, criterion.best, criterion.worst, criterion.weight))
        score = None
        if eligibility == "eligible":
            score = derived(float(sum(terms, Fraction(0))), "1", tuple(score_inputs), (
                Citation("NASA-SE-6", "6.8: exact affine weighted sum, final binary conversion"),
            ), "Preference arithmetic only; neither input uncertainty nor hardware fitness is quantified.")
        results.append(AlternativeScore(name, eligibility, score, tuple(utilities), checks, missing))
    return tuple(results)
