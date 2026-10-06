"""Common comparison metrics for nozzle contour alternatives."""

from __future__ import annotations

from dataclasses import dataclass
from math import cos, radians

from .nozzle_contour import NozzleContour


@dataclass(frozen=True, slots=True)
class ContourComparison:
    """Decision metrics derived from the common :class:`NozzleContour` contract."""

    contour: str
    model_fidelity: str
    total_length_m: float
    diverging_length_m: float
    surface_area_m2: float | None
    estimated_divergence_efficiency: float | None
    estimated_divergence_loss_fraction: float | None
    exit_flow_angle_deg: float
    exit_uniformity_estimate: float
    convergence_delta_exit_mach: float | None


def _exit_angle(contour: NozzleContour) -> float:
    if contour.contour == "conical":
        return contour.half_angle_deg
    if contour.contour == "preliminary":
        return 8.0
    if contour.contour == "rao":
        return 8.0
    if contour.contour.startswith("moc"):
        return 0.0
    return contour.half_angle_deg


def compare_nozzle_contours(
    contours: tuple[NozzleContour, ...] | list[NozzleContour],
) -> tuple[ContourComparison, ...]:
    """Normalize geometry, divergence, and one-dimensional exit metrics."""

    records: list[ContourComparison] = []
    for contour in contours:
        exit_angle = _exit_angle(contour)
        efficiency = contour.estimated_divergence_efficiency
        records.append(
            ContourComparison(
                contour=contour.contour,
                model_fidelity=contour.model_fidelity,
                total_length_m=contour.total_length_m,
                diverging_length_m=contour.diverging_length_m,
                surface_area_m2=contour.surface_area_m2,
                estimated_divergence_efficiency=efficiency,
                estimated_divergence_loss_fraction=(
                    None if efficiency is None else 1.0 - efficiency
                ),
                exit_flow_angle_deg=exit_angle,
                exit_uniformity_estimate=0.5 * (1.0 + cos(radians(exit_angle))),
                convergence_delta_exit_mach=contour.convergence_delta_exit_mach,
            )
        )
    return tuple(records)
