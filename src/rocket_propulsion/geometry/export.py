"""Portable text exports for axisymmetric nozzle contours."""

from __future__ import annotations

from dataclasses import dataclass
from itertools import pairwise
from math import cos, isfinite, pi, sin

from rocket_propulsion.core.errors import DomainError

from .nozzle_contour import NozzleContour, NozzleStation


@dataclass(frozen=True, slots=True)
class ContourExportBundle:
    """Self-contained engineering exchange files encoded as UTF-8 text."""

    csv: str
    svg: str
    dxf: str
    stl: str


def _csv(contour: NozzleContour) -> str:
    header = (
        "x_m,radius_m,area_m2,area_ratio,mach,"
        "pressure_total_ratio,temperature_total_ratio"
    )
    rows = [header]
    rows.extend(
        ",".join(
            f"{value:.12g}"
            for value in (
                station.x_m,
                station.radius_m,
                station.area_m2,
                station.area_ratio,
                station.mach,
                station.pressure_total_ratio,
                station.temperature_total_ratio,
            )
        )
        for station in contour.stations
    )
    return "\n".join(rows) + "\n"


def _svg(contour: NozzleContour) -> str:
    stations = contour.stations
    minimum_x = stations[0].x_m
    maximum_x = stations[-1].x_m
    maximum_radius = max(station.radius_m for station in stations)
    width, height, margin = 1200.0, 500.0, 30.0
    x_scale = (width - 2.0 * margin) / max(maximum_x - minimum_x, 1.0e-12)
    y_scale = (height - 2.0 * margin) / max(2.0 * maximum_radius, 1.0e-12)

    def point(station: NozzleStation, sign: float) -> str:
        x = margin + (station.x_m - minimum_x) * x_scale
        y = height / 2.0 - sign * station.radius_m * y_scale
        return f"{x:.3f},{y:.3f}"

    upper = " ".join(point(station, 1.0) for station in stations)
    lower = " ".join(point(station, -1.0) for station in reversed(stations))
    return (
        '<?xml version="1.0" encoding="UTF-8"?>\n'
        f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {width:.0f} {height:.0f}" '
        'role="img" aria-label="Axisymmetric nozzle contour">\n'
        '  <rect width="100%" height="100%" fill="white"/>\n'
        f'  <polygon points="{upper} {lower}" fill="#e7f1f4" stroke="#0b4f6c" '
        'stroke-width="3"/>\n'
        f'  <line x1="{margin:.3f}" y1="{height / 2.0:.3f}" '
        f'x2="{width - margin:.3f}" y2="{height / 2.0:.3f}" '
        'stroke="#6b7280" stroke-dasharray="8 6"/>\n'
        f'  <metadata>model={contour.contour}; units=m; area_ratio={contour.area_ratio:.12g}</metadata>\n'
        '</svg>\n'
    )


def _dxf(contour: NozzleContour) -> str:
    points = [(station.x_m, station.radius_m) for station in contour.stations]
    points.extend((station.x_m, -station.radius_m) for station in reversed(contour.stations))
    output = ["0", "SECTION", "2", "HEADER", "9", "$INSUNITS", "70", "6", "0", "ENDSEC"]
    output.extend(["0", "SECTION", "2", "ENTITIES", "0", "LWPOLYLINE", "8", "NOZZLE"])
    output.extend(["90", str(len(points)), "70", "1"])
    for x_m, radius_m in points:
        output.extend(["10", f"{x_m:.12g}", "20", f"{radius_m:.12g}"])
    output.extend(["0", "ENDSEC", "0", "EOF"])
    return "\n".join(output) + "\n"


def _normal_and_triangle(
    a: tuple[float, float, float],
    b: tuple[float, float, float],
    c: tuple[float, float, float],
) -> str:
    ab = tuple(right - left for left, right in zip(a, b))
    ac = tuple(right - left for left, right in zip(a, c))
    normal = (
        ab[1] * ac[2] - ab[2] * ac[1],
        ab[2] * ac[0] - ab[0] * ac[2],
        ab[0] * ac[1] - ab[1] * ac[0],
    )
    magnitude = sum(value * value for value in normal) ** 0.5
    unit = tuple(value / magnitude for value in normal) if magnitude else (0.0, 0.0, 0.0)
    lines = [f"  facet normal {unit[0]:.9g} {unit[1]:.9g} {unit[2]:.9g}", "    outer loop"]
    lines.extend(f"      vertex {x:.12g} {y:.12g} {z:.12g}" for x, y, z in (a, b, c))
    lines.extend(["    endloop", "  endfacet"])
    return "\n".join(lines)


def _stl(contour: NozzleContour, angular_segments: int) -> str:
    rings = [
        [
            (
                station.x_m,
                station.radius_m * cos(2.0 * pi * index / angular_segments),
                station.radius_m * sin(2.0 * pi * index / angular_segments),
            )
            for index in range(angular_segments)
        ]
        for station in contour.stations
    ]
    triangles: list[str] = []
    for left, right in pairwise(rings):
        for index in range(angular_segments):
            following = (index + 1) % angular_segments
            triangles.append(_normal_and_triangle(left[index], right[index], right[following]))
            triangles.append(_normal_and_triangle(left[index], right[following], left[following]))
    inlet_center = (contour.stations[0].x_m, 0.0, 0.0)
    exit_center = (contour.stations[-1].x_m, 0.0, 0.0)
    for index in range(angular_segments):
        following = (index + 1) % angular_segments
        triangles.append(_normal_and_triangle(inlet_center, rings[0][following], rings[0][index]))
        triangles.append(_normal_and_triangle(exit_center, rings[-1][index], rings[-1][following]))
    return "solid rocket_propulsion_nozzle\n" + "\n".join(triangles) + "\nendsolid rocket_propulsion_nozzle\n"


def export_nozzle_contour(
    contour: NozzleContour, *, angular_segments: int = 32
) -> ContourExportBundle:
    """Return CSV, SVG, ASCII DXF, and watertight ASCII STL representations."""

    if not 8 <= angular_segments <= 128:
        raise DomainError("STL angular segment count must be between 8 and 128.")
    if not contour.stations:
        raise DomainError("A contour must contain at least one station.")
    if any(
        not isfinite(value)
        for station in contour.stations
        for value in (station.x_m, station.radius_m)
    ):
        raise DomainError("Contour stations must be finite before export.")
    return ContourExportBundle(
        csv=_csv(contour),
        svg=_svg(contour),
        dxf=_dxf(contour),
        stl=_stl(contour, angular_segments),
    )
