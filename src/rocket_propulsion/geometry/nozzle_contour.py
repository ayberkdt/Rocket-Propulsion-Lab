"""Parametric axisymmetric converging-diverging nozzle contours."""

from dataclasses import dataclass
from itertools import pairwise
from math import cos, isfinite, pi, radians, sqrt, tan

from rocket_propulsion.compressible.isentropic import (
    calculate_isentropic,
    mach_from_isentropic,
)
from rocket_propulsion.core.errors import DomainError


@dataclass(frozen=True, slots=True)
class NozzleStation:
    """Geometry and ideal-flow properties at one axial location."""

    x_m: float
    radius_m: float
    area_m2: float
    area_ratio: float
    mach: float
    pressure_total_ratio: float
    temperature_total_ratio: float


@dataclass(frozen=True, slots=True)
class NozzleContour:
    """Dimensions and sampled stations for an axisymmetric C-D nozzle."""

    contour: str
    throat_radius_m: float
    chamber_radius_m: float
    exit_radius_m: float
    converging_length_m: float
    diverging_length_m: float
    total_length_m: float
    throat_area_m2: float
    exit_area_m2: float
    area_ratio: float
    half_angle_deg: float
    stations: tuple[NozzleStation, ...]
    model_fidelity: str = "preliminary"
    wall_tangency_residual_deg: float | None = None
    centerline_symmetry_residual_deg: float | None = None
    convergence_delta_exit_mach: float | None = None
    estimated_divergence_efficiency: float | None = None
    surface_area_m2: float | None = None
    source: str | None = None


def _positive(value: float, label: str) -> None:
    if not isfinite(value) or value <= 0.0:
        raise DomainError(f"{label} must be finite and greater than zero.")


def _hermite_radius(
    fraction: float,
    length: float,
    throat_radius: float,
    exit_radius: float,
) -> float:
    """Return a smooth bell-like radius using endpoint slope constraints."""

    h00 = 2.0 * fraction**3 - 3.0 * fraction**2 + 1.0
    h10 = fraction**3 - 2.0 * fraction**2 + fraction
    h01 = -2.0 * fraction**3 + 3.0 * fraction**2
    h11 = fraction**3 - fraction**2
    inlet_slope = tan(radians(30.0))
    exit_slope = tan(radians(8.0))
    return (
        h00 * throat_radius
        + h10 * length * inlet_slope
        + h01 * exit_radius
        + h11 * length * exit_slope
    )


def generate_nozzle_contour(
    *,
    throat_area_m2: float,
    area_ratio: float,
    gamma: float = 1.22,
    contour: str = "preliminary",
    half_angle_deg: float = 15.0,
    chamber_radius_ratio: float = 2.5,
    bell_length_fraction: float = 0.8,
    station_count: int = 121,
) -> NozzleContour:
    """Generate an axisymmetric nozzle contour and ideal station solution.

    The converging side uses a cosine contraction. A conical diverging side uses
    the selected half-angle; the preliminary option uses a cubic-Hermite curve
    with fixed 30-degree throat and 8-degree exit slopes. Use
    :func:`generate_rao_contour` when an explicit Rao/TOP construction is needed.
    """

    for value, label in (
        (throat_area_m2, "Throat area"),
        (area_ratio, "Area ratio"),
        (half_angle_deg, "Half-angle"),
        (chamber_radius_ratio, "Chamber radius ratio"),
        (bell_length_fraction, "Bell length fraction"),
    ):
        _positive(value, label)
    if area_ratio < 1.0:
        raise DomainError("Exit-to-throat area ratio must be at least one.")
    if not 1.0 < gamma:
        raise DomainError("Gamma must be greater than one.")
    normalized_contour = contour.strip().lower()
    if normalized_contour == "bell":
        normalized_contour = "preliminary"
    if normalized_contour not in {"conical", "preliminary"}:
        raise DomainError("Nozzle contour must be conical or preliminary.")
    if not 2.0 <= half_angle_deg < 45.0:
        raise DomainError("Conical half-angle must be in [2, 45) degrees.")
    if chamber_radius_ratio <= 1.0:
        raise DomainError("Chamber radius ratio must be greater than one.")
    if not 0.5 <= bell_length_fraction <= 1.2:
        raise DomainError("Bell length fraction must be in [0.5, 1.2].")
    if not 41 <= station_count <= 401:
        raise DomainError("Station count must be between 41 and 401.")

    throat_radius = sqrt(throat_area_m2 / pi)
    exit_radius = throat_radius * sqrt(area_ratio)
    chamber_radius = chamber_radius_ratio * throat_radius
    converging_length = 1.5 * throat_radius
    reference_cone_length = (exit_radius - throat_radius) / tan(radians(15.0))
    if normalized_contour == "conical":
        diverging_length = (exit_radius - throat_radius) / tan(radians(half_angle_deg))
    else:
        diverging_length = bell_length_fraction * reference_cone_length

    converging_count = max(15, round(station_count * 0.3))
    diverging_count = station_count - converging_count + 1
    stations: list[NozzleStation] = []

    def append_station(x_position: float, radius: float, branch: str) -> None:
        station_area = pi * radius**2
        local_area_ratio = station_area / throat_area_m2
        if abs(local_area_ratio - 1.0) < 1.0e-12:
            mach = 1.0
        else:
            mach = mach_from_isentropic(f"area_ratio_{branch}", local_area_ratio, gamma)
        ratios = calculate_isentropic(mach, gamma)
        stations.append(
            NozzleStation(
                x_m=x_position,
                radius_m=radius,
                area_m2=station_area,
                area_ratio=local_area_ratio,
                mach=mach,
                pressure_total_ratio=ratios.pressure_total_ratio,
                temperature_total_ratio=ratios.temperature_total_ratio,
            )
        )

    for index in range(converging_count):
        fraction = index / (converging_count - 1)
        x_position = -converging_length * (1.0 - fraction)
        radius = throat_radius + 0.5 * (chamber_radius - throat_radius) * (
            1.0 + cos(pi * fraction)
        )
        append_station(x_position, radius, "subsonic")

    for index in range(1, diverging_count):
        fraction = index / (diverging_count - 1)
        x_position = fraction * diverging_length
        if normalized_contour == "conical":
            radius = throat_radius + fraction * (exit_radius - throat_radius)
        else:
            radius = _hermite_radius(fraction, diverging_length, throat_radius, exit_radius)
        append_station(x_position, radius, "supersonic")

    surface_area = 0.0
    for left, right in pairwise(stations):
        segment = sqrt((right.x_m - left.x_m) ** 2 + (right.radius_m - left.radius_m) ** 2)
        surface_area += pi * (left.radius_m + right.radius_m) * segment
    effective_angle = half_angle_deg if normalized_contour == "conical" else 19.0
    return NozzleContour(
        contour=normalized_contour,
        throat_radius_m=throat_radius,
        chamber_radius_m=chamber_radius,
        exit_radius_m=exit_radius,
        converging_length_m=converging_length,
        diverging_length_m=diverging_length,
        total_length_m=converging_length + diverging_length,
        throat_area_m2=throat_area_m2,
        exit_area_m2=area_ratio * throat_area_m2,
        area_ratio=area_ratio,
        half_angle_deg=half_angle_deg,
        stations=tuple(stations),
        model_fidelity=(
            "analytical conical" if normalized_contour == "conical" else "preliminary Hermite"
        ),
        estimated_divergence_efficiency=0.5 * (1.0 + cos(radians(effective_angle))),
        surface_area_m2=surface_area,
        source=(
            "analytical cone"
            if normalized_contour == "conical"
            else "cubic Hermite preliminary contour; not Rao optimized"
        ),
    )

