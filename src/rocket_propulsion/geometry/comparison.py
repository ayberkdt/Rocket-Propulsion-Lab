"""Common comparison metrics for nozzle contour alternatives.

Divergence metrics use the conical divergence factor ``(1 + cos(alpha))/2``
[1]_.

References
----------
.. [1] SUTTON-BIBLARZ-2017. G. P. Sutton and O. Biblarz, "Rocket Propulsion
       Elements", 9th ed., Wiley, 2017.
"""

from __future__ import annotations

from dataclasses import dataclass
from math import cos, radians

from .nozzle_contour import NozzleContour


@dataclass(frozen=True, slots=True)
class ContourComparison:
    """Decision metrics derived from the common :class:`NozzleContour` contract.

    Attributes
    ----------
    contour, model_fidelity : str
        Contour family and fidelity label copied from the source contour.
    total_length_m, diverging_length_m : float
        Overall and diverging-section lengths in m.
    surface_area_m2 : float or None
        Wetted wall area from the sampled stations, if available.
    estimated_divergence_efficiency : float or None
        Divergence factor reported by the contour.
    estimated_divergence_loss_fraction : float or None
        ``1 - estimated_divergence_efficiency``.
    exit_flow_angle_deg : float
        Nominal exit flow angle assumed for the contour family.
    exit_uniformity_estimate : float
        ``(1 + cos(exit_flow_angle))/2`` evaluated at the nominal exit angle
        [1]_.
    convergence_delta_exit_mach : float or None
        Exit-Mach discrepancy reported by the contour generator, if any.

    References
    ----------
    .. [1] SUTTON-BIBLARZ-2017. G. P. Sutton and O. Biblarz, "Rocket Propulsion
           Elements", 9th ed., Wiley, 2017.
    """

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
    """Return the nominal exit flow angle in degrees for a contour family.

    Conical contours use their half-angle. The preliminary Hermite and Rao
    parabola contours are assigned a nominal 8 degree exit angle (the default
    exit slope of both constructions), independent of a non-default Rao exit
    angle. Method-of-characteristics contours are designed for axial exit flow
    (0 degrees). Unknown families fall back to the reported half-angle.
    """
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
    """Normalize geometry, divergence, and one-dimensional exit metrics.

    Parameters
    ----------
    contours : sequence of NozzleContour
        Contours to compare, in the order the records are returned.

    Returns
    -------
    tuple of ContourComparison
        One record per input contour.

    Notes
    -----
    The exit-uniformity estimate applies the conical divergence factor
    ``lambda = (1 + cos(alpha))/2`` [1]_ to the nominal exit flow angle of each
    family. It is a one-dimensional screening metric, not a computed
    exit-plane flow profile.

    References
    ----------
    .. [1] SUTTON-BIBLARZ-2017. G. P. Sutton and O. Biblarz, "Rocket Propulsion
           Elements", 9th ed., Wiley, 2017.
    """

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
