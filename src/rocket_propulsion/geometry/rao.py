"""Rao-style thrust-optimized parabolic bell nozzle construction."""

from __future__ import annotations

from itertools import pairwise
from math import cos, isfinite, pi, radians, sqrt, tan

from rocket_propulsion.compressible import calculate_isentropic, mach_from_isentropic
from rocket_propulsion.core.errors import DomainError

from .nozzle_contour import NozzleContour, NozzleStation


def generate_rao_contour(
    *,
    throat_area_m2: float,
    area_ratio: float,
    gamma: float = 1.22,
    length_fraction: float = 0.8,
    throat_angle_deg: float = 30.0,
    exit_angle_deg: float = 8.0,
    chamber_radius_ratio: float = 2.5,
    station_count: int = 121,
) -> NozzleContour:
    """Generate a quadratic Rao/TOP approximation with explicit end angles.

    The control point is the intersection of the throat and exit tangent lines.
    This is the standard thrust-optimized-parabola construction, not a digitized
    reproduction of a particular Rao optimization chart.
    """

    values = (
        (throat_area_m2, "throat_area_m2"),
        (area_ratio, "area_ratio"),
        (length_fraction, "length_fraction"),
        (chamber_radius_ratio, "chamber_radius_ratio"),
    )
    for value, field in values:
        if not isfinite(value) or value <= 0.0:
            raise DomainError(f"{field} must be finite and positive.", field=field)
    if area_ratio <= 1.0:
        raise DomainError("Rao area ratio must exceed one.", field="area_ratio")
    if not 0.6 <= length_fraction <= 1.0:
        raise DomainError("Rao length fraction must be in [0.6, 1.0].")
    if not 15.0 <= throat_angle_deg <= 40.0:
        raise DomainError("Rao throat angle must be in [15, 40] degrees.")
    if not 0.0 <= exit_angle_deg < throat_angle_deg:
        raise DomainError("Rao exit angle must be non-negative and below throat angle.")
    if chamber_radius_ratio <= 1.0:
        raise DomainError("Chamber radius ratio must exceed one.")
    if not 41 <= station_count <= 401:
        raise DomainError("Station count must be between 41 and 401.")

    throat_radius = sqrt(throat_area_m2 / pi)
    exit_radius = throat_radius * sqrt(area_ratio)
    chamber_radius = chamber_radius_ratio * throat_radius
    converging_length = 1.5 * throat_radius
    cone_length = (exit_radius - throat_radius) / tan(radians(15.0))
    diverging_length = length_fraction * cone_length
    throat_slope = tan(radians(throat_angle_deg))
    exit_slope = tan(radians(exit_angle_deg))
    control_x = (
        exit_radius - throat_radius - exit_slope * diverging_length
    ) / (throat_slope - exit_slope)
    control_y = throat_radius + throat_slope * control_x
    if not 0.0 < control_x < diverging_length:
        raise DomainError(
            "Rao tangent lines do not intersect inside the nozzle length.",
            code="invalid_rao_control_point",
            details={"control_x_m": control_x, "diverging_length_m": diverging_length},
        )

    converging_count = max(15, round(station_count * 0.3))
    diverging_count = station_count - converging_count + 1
    stations: list[NozzleStation] = []

    def append(x_m: float, radius_m: float, branch: str) -> None:
        area = pi * radius_m**2
        local_ratio = area / throat_area_m2
        mach = (
            1.0
            if abs(local_ratio - 1.0) <= 1.0e-12
            else mach_from_isentropic(f"area_ratio_{branch}", local_ratio, gamma)
        )
        flow = calculate_isentropic(mach, gamma)
        stations.append(
            NozzleStation(
                x_m=x_m,
                radius_m=radius_m,
                area_m2=area,
                area_ratio=local_ratio,
                mach=mach,
                pressure_total_ratio=flow.pressure_total_ratio,
                temperature_total_ratio=flow.temperature_total_ratio,
            )
        )

    for index in range(converging_count):
        fraction = index / (converging_count - 1)
        radius = chamber_radius + (throat_radius - chamber_radius) * (
            0.5 - 0.5 * cos(pi * fraction)
        )
        append(-converging_length * (1.0 - fraction), radius, "subsonic")
    for index in range(1, diverging_count):
        fraction = index / (diverging_count - 1)
        one_minus = 1.0 - fraction
        x_m = 2.0 * one_minus * fraction * control_x + fraction**2 * diverging_length
        radius = (
            one_minus**2 * throat_radius
            + 2.0 * one_minus * fraction * control_y
            + fraction**2 * exit_radius
        )
        append(x_m, radius, "supersonic")

    surface_area = sum(
        pi
        * (left.radius_m + right.radius_m)
        * sqrt((right.x_m - left.x_m) ** 2 + (right.radius_m - left.radius_m) ** 2)
        for left, right in pairwise(stations)
    )
    effective_angle = 0.5 * (throat_angle_deg + exit_angle_deg)
    return NozzleContour(
        contour="rao",
        throat_radius_m=throat_radius,
        chamber_radius_m=chamber_radius,
        exit_radius_m=exit_radius,
        converging_length_m=converging_length,
        diverging_length_m=diverging_length,
        total_length_m=converging_length + diverging_length,
        throat_area_m2=throat_area_m2,
        exit_area_m2=throat_area_m2 * area_ratio,
        area_ratio=area_ratio,
        half_angle_deg=effective_angle,
        stations=tuple(stations),
        model_fidelity="Rao thrust-optimized parabola approximation",
        wall_tangency_residual_deg=0.0,
        estimated_divergence_efficiency=0.5 * (1.0 + cos(radians(effective_angle))),
        surface_area_m2=surface_area,
        source="Rao-style quadratic tangent construction; NASA SP-8120 context",
    )

