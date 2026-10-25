"""Numerical comparison of scalar results shared by two study cases."""

from __future__ import annotations

from dataclasses import dataclass
from math import isfinite

from rocket_propulsion.core.errors import DomainError

from .workspace import StudyCase


@dataclass(frozen=True, slots=True)
class ResultDifference:
    """Difference of one scalar metric between a baseline and a candidate case.

    Attributes
    ----------
    metric : str
        Result name.
    baseline_value, candidate_value : float
        Finite metric values of the two cases.
    absolute_difference : float
        candidate_value - baseline_value.
    relative_difference : float or None
        absolute_difference / abs(baseline_value); None when the baseline is zero.
    """
    metric: str
    baseline_value: float
    candidate_value: float
    absolute_difference: float
    relative_difference: float | None


@dataclass(frozen=True, slots=True)
class CaseComparison:
    """Ordered per-metric differences between two study cases.

    Attributes
    ----------
    baseline_case_id, candidate_case_id : str
        Identifiers of the compared cases.
    differences : tuple of ResultDifference
        One record per compared metric.
    """
    baseline_case_id: str
    candidate_case_id: str
    differences: tuple[ResultDifference, ...]


def compare_cases(
    baseline: StudyCase,
    candidate: StudyCase,
    metrics: tuple[str, ...] | None = None,
) -> CaseComparison:
    """Compare finite scalar outputs of two cases relative to the baseline.

    Parameters
    ----------
    baseline, candidate : StudyCase
        Cases whose ``results`` are compared.
    metrics : tuple of str, optional
        Metrics to compare, in order. By default every result name shared by both
        cases is compared in sorted order.

    Returns
    -------
    CaseComparison
        Absolute and baseline-relative differences.

    Raises
    ------
    DomainError
        If a selected metric is missing, not numeric, or not finite in either case.
    """

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
