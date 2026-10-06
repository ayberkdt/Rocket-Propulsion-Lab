"""Cross-case numerical result comparison."""

from __future__ import annotations

from dataclasses import dataclass
from math import isfinite

from rocket_propulsion.core.errors import DomainError

from .workspace import StudyCase


@dataclass(frozen=True, slots=True)
class ResultDifference:
    metric: str
    baseline_value: float
    candidate_value: float
    absolute_difference: float
    relative_difference: float | None


@dataclass(frozen=True, slots=True)
class CaseComparison:
    baseline_case_id: str
    candidate_case_id: str
    differences: tuple[ResultDifference, ...]


def compare_cases(
    baseline: StudyCase,
    candidate: StudyCase,
    metrics: tuple[str, ...] | None = None,
) -> CaseComparison:
    """Compare shared finite scalar outputs and retain baseline-relative deltas."""

    selected = metrics or tuple(
        sorted(set(baseline.results).intersection(candidate.results))
    )
    differences: list[ResultDifference] = []
    for metric in selected:
        try:
            left = float(baseline.results[metric])
            right = float(candidate.results[metric])
        except (KeyError, TypeError, ValueError) as error:
            raise DomainError(
                f"Metric {metric} is not a shared numeric result.", field=metric
            ) from error
        if not isfinite(left) or not isfinite(right):
            raise DomainError(f"Metric {metric} must be finite.", field=metric)
        difference = right - left
        differences.append(
            ResultDifference(
                metric=metric,
                baseline_value=left,
                candidate_value=right,
                absolute_difference=difference,
                relative_difference=None if left == 0.0 else difference / abs(left),
            )
        )
    return CaseComparison(baseline.case_id, candidate.case_id, tuple(differences))
