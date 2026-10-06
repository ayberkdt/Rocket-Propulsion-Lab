"""Planar minimum-length nozzle method of characteristics.

The characteristic compatibility equations are solved in the planar field.
An optional first-order axisymmetric area mapping preserves the requested
throat and exit radii; it is deliberately labelled as a correction rather than
a full axisymmetric source-term MOC solution.
"""

from __future__ import annotations

from bisect import bisect_right
from dataclasses import dataclass
from itertools import pairwise
from math import asin, atan2, cos, degrees, isfinite, pi, radians, sqrt, tan

from rocket_propulsion.compressible import calculate_isentropic, mach_from_isentropic
from rocket_propulsion.core.errors import ConvergenceError, DomainError

from .nozzle_contour import NozzleContour, NozzleStation


@dataclass(frozen=True, slots=True)
class CharacteristicPoint:
    """One solved node in the planar characteristic mesh."""

    row: int
    column: int
    x_m: float
    radius_m: float
    flow_angle_deg: float
    prandtl_meyer_deg: float
    mach: float
    mach_angle_deg: float
    k_plus_deg: float
    k_minus_deg: float
    boundary: str


@dataclass(frozen=True, slots=True)
class MocResiduals:
    """Boundary and target residuals for an MOC design."""

    wall_tangency_max_deg: float
    centerline_symmetry_max_deg: float
    exit_mach_absolute: float


@dataclass(frozen=True, slots=True)
class MocDesign:
    """Common nozzle contour plus the characteristic net and evidence."""

    contour: NozzleContour
    characteristic_lines: tuple[tuple[CharacteristicPoint, ...], ...]
    geometry_mode: str
    exit_mach_target: float
    exit_mach_computed: float
    characteristic_count: int
    residuals: MocResiduals
    refined_length_change_relative: float
    source: str
    source_url: str


@dataclass(slots=True)
class _Point:
    theta: float
    nu: float
    mach: float
    mu: float
    k_plus: float
    k_minus: float
    x: float = 0.0
    y: float = 0.0
    boundary: str = "interior"


def _intersection(
    x1: float,
    y1: float,
    slope1: float,
    x2: float,
    y2: float,
    slope2: float,
) -> tuple[float, float]:
    denominator = slope1 - slope2
    if abs(denominator) < 1.0e-12:
        raise ConvergenceError(
            "MOC characteristic lines became parallel.",
            code="moc_parallel_characteristics",
        )
    x = (y2 - y1 + slope1 * x1 - slope2 * x2) / denominator
    return x, y1 + slope1 * (x - x1)


def _state_gamma(theta: float, nu: float, gamma: float) -> _Point:
    if nu < 0.0:
        raise ConvergenceError(
            "MOC produced a negative Prandtl-Meyer angle.",
            code="moc_invalid_state",
            details={"theta_deg": theta, "nu_deg": nu},
        )
    mach = mach_from_isentropic("prandtl_meyer_deg", max(0.0, nu), gamma)
    mu = degrees(asin(1.0 / mach))
    return _Point(theta, nu, mach, mu, theta - nu, theta + nu)


def _normalized_mesh(
    *, exit_mach: float, gamma: float, characteristic_count: int
) -> tuple[list[list[_Point]], list[_Point]]:
    exit_nu = calculate_isentropic(exit_mach, gamma).prandtl_meyer_deg
    assert exit_nu is not None
    maximum_wall_angle = 0.5 * exit_nu
    angle_step = maximum_wall_angle / characteristic_count
    rows: list[list[_Point]] = []
    wall: list[_Point] = []
    throat_wall = _state_gamma(maximum_wall_angle, maximum_wall_angle, gamma)
    throat_wall.x = 0.0
    throat_wall.y = 1.0
    throat_wall.boundary = "wall"
    wall.append(throat_wall)

    for row_index in range(characteristic_count):
        fan_index = row_index + 1
        fan_angle = fan_index * angle_step
        length = characteristic_count + 2 - fan_index
        row: list[_Point] = []
        for column in range(length):
            if column == 0:
                point = _state_gamma(0.0, 2.0 * fan_angle, gamma)
                point.boundary = "centerline"
            elif column == length - 1:
                incoming = row[column - 1]
                theta = maximum_wall_angle - fan_angle
                nu = theta - incoming.k_plus
                point = _state_gamma(theta, nu, gamma)
                point.boundary = "wall"
            else:
                lower = row[column - 1]
                initial_angle = (fan_index + column) * angle_step
                k_plus = lower.k_plus
                k_minus = 2.0 * initial_angle
                point = _state_gamma(
                    0.5 * (k_plus + k_minus),
                    0.5 * (k_minus - k_plus),
                    gamma,
                )

            if column == 0:
                if row_index == 0:
                    source = _state_gamma(fan_angle, fan_angle, gamma)
                    source_x, source_y = 0.0, 1.0
                else:
                    source = rows[row_index - 1][1]
                    source_x, source_y = source.x, source.y
                slope = tan(
                    radians(
                        0.5
                        * (
                            source.theta
                            - source.mu
                            + point.theta
                            - point.mu
                        )
                    )
                )
                if abs(slope) < 1.0e-12:
                    raise ConvergenceError("MOC centerline characteristic has zero slope.")
                point.x = source_x - source_y / slope
                point.y = 0.0
                point.boundary = "centerline"
            elif column == length - 1:
                lower = row[column - 1]
                previous_wall = wall[-1]
                characteristic_slope = tan(
                    radians(
                        0.5
                        * (
                            lower.theta
                            + lower.mu
                            + point.theta
                            + point.mu
                        )
                    )
                )
                wall_slope = tan(
                    radians(0.5 * (previous_wall.theta + point.theta))
                )
                point.x, point.y = _intersection(
                    lower.x,
                    lower.y,
                    characteristic_slope,
                    previous_wall.x,
                    previous_wall.y,
                    wall_slope,
                )
                wall.append(point)
            else:
                lower = row[column - 1]
                if row_index == 0:
                    initial_angle = (fan_index + column) * angle_step
                    upper_state = _state_gamma(initial_angle, initial_angle, gamma)
                    upper_x, upper_y = 0.0, 1.0
                else:
                    upper_state = rows[row_index - 1][column + 1]
                    upper_x, upper_y = upper_state.x, upper_state.y
                positive_slope = tan(
                    radians(
                        0.5
                        * (
                            lower.theta
                            + lower.mu
                            + point.theta
                            + point.mu
                        )
                    )
                )
                negative_slope = tan(
                    radians(
                        0.5
                        * (
                            upper_state.theta
                            - upper_state.mu
                            + point.theta
                            - point.mu
                        )
                    )
                )
                point.x, point.y = _intersection(
                    lower.x,
                    lower.y,
                    positive_slope,
                    upper_x,
                    upper_y,
                    negative_slope,
                )
            row.append(point)
        rows.append(row)
    return rows, wall


def _build_design(
    *,
    throat_area_m2: float,
    area_ratio: float,
    gamma: float,
    characteristic_count: int,
    axisymmetric_correction: bool,
) -> tuple[NozzleContour, tuple[tuple[CharacteristicPoint, ...], ...], MocResiduals]:
    exit_mach = mach_from_isentropic("area_ratio_supersonic", area_ratio, gamma)
    rows, wall = _normalized_mesh(
        exit_mach=exit_mach,
        gamma=gamma,
        characteristic_count=characteristic_count,
    )
    throat_radius = sqrt(throat_area_m2 / pi)
    desired_exit_ratio = sqrt(area_ratio)
    raw_exit_ratio = wall[-1].y
    if raw_exit_ratio <= 1.0:
        raise ConvergenceError(
            "MOC wall did not expand beyond the throat.",
            code="moc_invalid_contour",
        )

    def mapped_wall_radius(normalized_radius: float) -> float:
        fraction = (normalized_radius - 1.0) / (raw_exit_ratio - 1.0)
        return throat_radius * (1.0 + fraction * (desired_exit_ratio - 1.0))

    x_scale = throat_radius * (
        (desired_exit_ratio - 1.0) / (raw_exit_ratio - 1.0)
    )
    wall_x = [point.x for point in wall]

    def wall_height_at(x_position: float) -> float:
        index = min(max(1, bisect_right(wall_x, x_position)), len(wall_x) - 1)
        left = wall[index - 1]
        right = wall[index]
        if abs(right.x - left.x) < 1.0e-12:
            return right.y
        fraction = (x_position - left.x) / (right.x - left.x)
        return left.y + fraction * (right.y - left.y)

    def mapped_mesh_radius(point: _Point) -> float:
        if point.boundary == "centerline":
            return 0.0
        local_wall_height = max(wall_height_at(point.x), point.y, 1.0e-12)
        local_wall_radius = mapped_wall_radius(local_wall_height)
        return local_wall_radius * point.y / local_wall_height
    characteristic_lines: list[tuple[CharacteristicPoint, ...]] = []
    for row_index, row in enumerate(rows):
        points = tuple(
            CharacteristicPoint(
                row=row_index,
                column=column,
                x_m=point.x * x_scale,
                radius_m=mapped_mesh_radius(point),
                flow_angle_deg=point.theta,
                prandtl_meyer_deg=point.nu,
                mach=point.mach,
                mach_angle_deg=point.mu,
                k_plus_deg=point.k_plus,
                k_minus_deg=point.k_minus,
                boundary=point.boundary,
            )
            for column, point in enumerate(row)
        )
        characteristic_lines.append(points)

    chamber_radius = 2.5 * throat_radius
    converging_length = 1.5 * throat_radius
    stations: list[NozzleStation] = []
    converging_count = max(15, characteristic_count)
    for index in range(converging_count):
        fraction = index / (converging_count - 1)
        radius = chamber_radius + (throat_radius - chamber_radius) * (
            0.5 - 0.5 * cos(pi * fraction)
        )
        local_area = pi * radius**2
        local_ratio = local_area / throat_area_m2
        mach = (
            1.0
            if index == converging_count - 1
            else mach_from_isentropic("area_ratio_subsonic", local_ratio, gamma)
        )
        flow = calculate_isentropic(mach, gamma)
        stations.append(
            NozzleStation(
                x_m=-converging_length * (1.0 - fraction),
                radius_m=radius,
                area_m2=local_area,
                area_ratio=local_ratio,
                mach=mach,
                pressure_total_ratio=flow.pressure_total_ratio,
                temperature_total_ratio=flow.temperature_total_ratio,
            )
        )
    for index, point in enumerate(wall[1:], start=1):
        radius = mapped_wall_radius(point.y)
        local_ratio = (radius / throat_radius) ** 2
        mach = point.mach
        flow = calculate_isentropic(mach, gamma)
        stations.append(
            NozzleStation(
                x_m=point.x * x_scale,
                radius_m=radius,
                area_m2=pi * radius**2,
                area_ratio=local_ratio,
                mach=mach,
                pressure_total_ratio=flow.pressure_total_ratio,
                temperature_total_ratio=flow.temperature_total_ratio,
            )
        )
    wall_residuals = []
    mapped_wall = [
        (point.x * x_scale, mapped_wall_radius(point.y), point.theta) for point in wall
    ]
    for left, right in pairwise(mapped_wall):
        segment_angle = degrees(
            atan2(right[1] - left[1], right[0] - left[0])
        )
        wall_residuals.append(abs(segment_angle - 0.5 * (left[2] + right[2])))
    residuals = MocResiduals(
        wall_tangency_max_deg=max(wall_residuals, default=0.0),
        centerline_symmetry_max_deg=max(
            abs(point.theta)
            for row in rows
            for point in row
            if point.boundary == "centerline"
        ),
        exit_mach_absolute=abs(wall[-1].mach - exit_mach),
    )
    surface_area = sum(
        pi
        * (left.radius_m + right.radius_m)
        * sqrt((right.x_m - left.x_m) ** 2 + (right.radius_m - left.radius_m) ** 2)
        for left, right in pairwise(stations)
    )
    diverging_length = stations[-1].x_m
    contour = NozzleContour(
        contour="moc-axisymmetric-corrected" if axisymmetric_correction else "moc-planar",
        throat_radius_m=throat_radius,
        chamber_radius_m=chamber_radius,
        exit_radius_m=throat_radius * desired_exit_ratio,
        converging_length_m=converging_length,
        diverging_length_m=diverging_length,
        total_length_m=converging_length + diverging_length,
        throat_area_m2=throat_area_m2,
        exit_area_m2=throat_area_m2 * area_ratio,
        area_ratio=area_ratio,
        half_angle_deg=0.5
        * (calculate_isentropic(exit_mach, gamma).prandtl_meyer_deg or 0.0),
        stations=tuple(stations),
        model_fidelity=(
            "planar MOC compatibility with first-order axisymmetric area mapping"
            if axisymmetric_correction
            else "planar minimum-length MOC with equivalent circular-area contour"
        ),
        wall_tangency_residual_deg=residuals.wall_tangency_max_deg,
        centerline_symmetry_residual_deg=residuals.centerline_symmetry_max_deg,
        convergence_delta_exit_mach=residuals.exit_mach_absolute,
        estimated_divergence_efficiency=1.0,
        surface_area_m2=surface_area,
        source="planar characteristic compatibility; NASA MOC design context",
    )
    return contour, tuple(characteristic_lines), residuals


def generate_moc_nozzle(
    *,
    throat_area_m2: float,
    area_ratio: float,
    gamma: float = 1.22,
    characteristic_count: int = 20,
    axisymmetric_correction: bool = True,
) -> MocDesign:
    """Generate a minimum-length characteristic contour and refined evidence."""

    if not isfinite(throat_area_m2) or throat_area_m2 <= 0.0:
        raise DomainError("Throat area must be finite and positive.")
    if not isfinite(area_ratio) or area_ratio <= 1.0:
        raise DomainError("MOC area ratio must exceed one.")
    if not isfinite(gamma) or gamma <= 1.0:
        raise DomainError("Gamma must exceed one.")
    if not 5 <= characteristic_count <= 80:
        raise DomainError("Characteristic count must be between 5 and 80.")
    contour, lines, residuals = _build_design(
        throat_area_m2=throat_area_m2,
        area_ratio=area_ratio,
        gamma=gamma,
        characteristic_count=characteristic_count,
        axisymmetric_correction=axisymmetric_correction,
    )
    refined, _, _ = _build_design(
        throat_area_m2=throat_area_m2,
        area_ratio=area_ratio,
        gamma=gamma,
        characteristic_count=min(80, characteristic_count * 2),
        axisymmetric_correction=axisymmetric_correction,
    )
    relative_change = abs(refined.diverging_length_m - contour.diverging_length_m) / max(
        refined.diverging_length_m, 1.0e-12
    )
    return MocDesign(
        contour=contour,
        characteristic_lines=lines,
        geometry_mode=(
            "axisymmetric-area-corrected" if axisymmetric_correction else "planar"
        ),
        exit_mach_target=mach_from_isentropic(
            "area_ratio_supersonic", area_ratio, gamma
        ),
        exit_mach_computed=contour.stations[-1].mach,
        characteristic_count=characteristic_count,
        residuals=residuals,
        refined_length_change_relative=relative_change,
        source="Method of Characteristics; NASA nozzle-design literature",
        source_url="https://www.nasa.gov/glenn/research/inlets-and-nozzles/design-analysis-software/",
    )

