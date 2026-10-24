"""Replayed preliminary text exchange of a closed sampled fluid domain."""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass, fields
from fractions import Fraction
from functools import lru_cache
from itertools import pairwise
from xml.sax.saxutils import escape

from rocket_propulsion.core.errors import DomainError
from rocket_propulsion.core.validation import _ordered_items, integer, real, text

from .nozzle_contour import NozzleContour, NozzleStation
from .export_numerics import _check_geometry, _ring, _unit_normal

_SERIALIZATION_POLICY = 'binary64_round_trip_17_digits_exact_vertex_cross_normals_v1'
_MESH_ROLE = 'closed_fluid_domain_visualization'


def _owned(value, kind, label):
    """Copy retained records, typing malformed protocols without hiding callbacks."""
    try:
        if isinstance(value, Mapping):
            return kind(**dict(value))
        if isinstance(value, kind):
            return kind(**{item.name: getattr(value, item.name) for item in fields(kind)})
    except (TypeError, ValueError, KeyError, IndexError, AttributeError, OverflowError) as error:
        raise DomainError('Malformed retained export record.', field=label) from error
    raise DomainError('Typed retained export record required.', field=label)


@dataclass(frozen=True, slots=True)
class ContourExportInputs:
    """Owned quantities consumed by export; no authentication of a flow solver."""
    contour: str
    throat_area_m2: float
    area_ratio: float
    stations: tuple[NozzleStation, ...]
    angular_segments: int = 32

    def __post_init__(self):
        name = text(self.contour, 'Export contour name')
        if any((ord(c) < 32 and c not in '\t\n\r') or ord(c) in (65534, 65535) for c in name):
            raise DomainError('Contour name contains characters forbidden in XML 1.0.', field='Export contour name')
        object.__setattr__(self, 'contour', name)
        for key in ('throat_area_m2', 'area_ratio'):
            value = real(getattr(self, key), key)
            if value <= 0 or (key == 'area_ratio' and value < 1):
                raise DomainError('Throat area must be positive and exit area ratio at least one.', field=key)
            object.__setattr__(self, key, value)
        count = integer(self.angular_segments, 'angular_segments')
        if not 8 <= count <= 128:
            raise DomainError('STL angular segment count must be between 8 and 128.', field='angular_segments')
        object.__setattr__(self, 'angular_segments', count)
        points = []
        for station in _ordered_items(self.stations, 'Export stations'):
            station = _owned(station, NozzleStation, 'Export station')
            values = {item.name: real(getattr(station, item.name), item.name) for item in fields(NozzleStation)}
            station = NozzleStation(**values)
            if any(values[key] <= 0 for key in ('radius_m', 'area_m2', 'area_ratio', 'pressure_total_ratio', 'temperature_total_ratio')) or station.mach < 0:
                raise DomainError('Export station contains invalid positive quantities or negative Mach.', field='Export station')
            if points and station.x_m <= points[-1].x_m:
                raise DomainError('Export axial stations must increase strictly.', field='Export stations')
            points.append(station)
        if len(points) < 2:
            raise DomainError('An export needs at least two axial stations.', field='Export stations')
        points = tuple(points)
        _check_geometry(points, self.throat_area_m2, self.area_ratio)
        object.__setattr__(self, 'stations', points)


@dataclass(frozen=True, slots=True)
class ContourExportBundle:
    """Replayed files for visualization/exchange, not a nozzle wall design."""

    csv: str
    svg: str
    dxf: str
    stl: str
    inputs: ContourExportInputs | None = None
    coordinate_unit: str = 'm'
    mesh_role: str = _MESH_ROLE
    manufacturing_qualified: bool = False
    physical_model_qualified: bool = False
    serialization_policy: str = _SERIALIZATION_POLICY

    def __post_init__(self):
        source = _owned(self.inputs, ContourExportInputs, 'Export inputs')
        object.__setattr__(self, 'inputs', source)
        for name in ('csv', 'svg', 'dxf', 'stl', 'coordinate_unit', 'mesh_role', 'serialization_policy'):
            object.__setattr__(self, name, text(getattr(self, name), name))
        for name in ('manufacturing_qualified', 'physical_model_qualified'):
            if type(getattr(self, name)) is not bool:
                raise DomainError('Export qualification flags must be booleans.', field=name)
        if (self.coordinate_unit != 'm' or self.mesh_role != _MESH_ROLE
                or self.serialization_policy != _SERIALIZATION_POLICY
                or self.manufacturing_qualified or self.physical_model_qualified):
            raise DomainError('Invalid export policy or qualification claim.', field='Export bundle')
        expected = _render(source)
        if any(getattr(self, name) != value for name, value in zip(('csv', 'svg', 'dxf', 'stl'), expected)):
            raise DomainError('Export text does not replay from owned inputs.', field='Export bundle')


def _csv(contour: ContourExportInputs) -> str:
    header = (
        "x_m,radius_m,area_m2,area_ratio,mach,"
        "pressure_total_ratio,temperature_total_ratio"
    )
    rows = [header]
    rows.extend(
        ",".join(
            f"{value:.17g}"
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


def _svg(contour: ContourExportInputs) -> str:
    stations = contour.stations
    minimum_x = stations[0].x_m
    maximum_x = stations[-1].x_m
    maximum_radius = max(station.radius_m for station in stations)
    width, height, margin = 1200.0, 500.0, 30.0
    x_span = Fraction(maximum_x) - Fraction(minimum_x)
    radius = Fraction(maximum_radius)

    def point(station: NozzleStation, sign: float) -> str:
        x = float(30 + 1140 * (Fraction(station.x_m) - Fraction(minimum_x)) / x_span)
        y = float(250 - Fraction(sign) * 220 * Fraction(station.radius_m) / radius)
        return f"{x:.3f},{y:.3f}"

    upper = " ".join(point(station, 1.0) for station in stations)
    lower = " ".join(point(station, -1.0) for station in reversed(stations))
    name = escape(contour.contour, {'\r': '&#13;'})
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
        f'  <metadata>model={name}; units=m; area_ratio={contour.area_ratio:.17g}; role={_MESH_ROLE}; manufacturing_qualified=false</metadata>\n'
        '</svg>\n'
    )


def _dxf(contour: ContourExportInputs) -> str:
    points = [(station.x_m, station.radius_m) for station in contour.stations]
    points.extend((station.x_m, -station.radius_m) for station in reversed(contour.stations))
    output = ["0", "SECTION", "2", "HEADER", "9", "$ACADVER", "1", "AC1015", "9", "$INSUNITS", "70", "6", "0", "ENDSEC"]
    output.extend(["0", "SECTION", "2", "ENTITIES", "0", "LWPOLYLINE", "100", "AcDbEntity", "8", "NOZZLE", "100", "AcDbPolyline"])
    output.extend(["90", str(len(points)), "70", "1"])
    for x_m, radius_m in points:
        output.extend(["10", f"{x_m:.17g}", "20", f"{radius_m:.17g}"])
    output.extend(["0", "ENDSEC", "0", "EOF"])
    return "\n".join(output) + "\n"


def _normal_and_triangle(
    a: tuple[float, float, float],
    b: tuple[float, float, float],
    c: tuple[float, float, float],
) -> str:
    unit = _unit_normal(a, b, c)
    lines = [f"  facet normal {unit[0]:.17g} {unit[1]:.17g} {unit[2]:.17g}", "    outer loop"]
    lines.extend(f"      vertex {x:.17g} {y:.17g} {z:.17g}" for x, y, z in (a, b, c))
    lines.extend(["    endloop", "  endfacet"])
    return "\n".join(lines)


def _stl(contour: ContourExportInputs, angular_segments: int) -> str:
    rings = [_ring(station, angular_segments) for station in contour.stations]
    triangles: list[str] = []
    for left, right in pairwise(rings):
        for index in range(angular_segments):
            following = (index + 1) % angular_segments
            triangles.append(_normal_and_triangle(left[index], right[following], right[index]))
            triangles.append(_normal_and_triangle(left[index], left[following], right[following]))
    inlet_center = (contour.stations[0].x_m, 0.0, 0.0)
    exit_center = (contour.stations[-1].x_m, 0.0, 0.0)
    for index in range(angular_segments):
        following = (index + 1) % angular_segments
        triangles.append(_normal_and_triangle(inlet_center, rings[0][following], rings[0][index]))
        triangles.append(_normal_and_triangle(exit_center, rings[-1][index], rings[-1][following]))
    return "solid rocket_propulsion_nozzle\n" + "\n".join(triangles) + "\nendsolid rocket_propulsion_nozzle\n"


@lru_cache(maxsize=1)
def _render(inputs):
    """Bound retained mesh memory while allowing factory/replay to share rendering."""
    return _csv(inputs), _svg(inputs), _dxf(inputs), _stl(inputs, inputs.angular_segments)


def export_nozzle_contour(
    contour: NozzleContour, *, angular_segments: int = 32
) -> ContourExportBundle:
    """Export a sampled closed fluid domain; no wall/flow/manufacturing qualification.

    CSV/DXF/STL use 17 significant digits for binary64 coordinate round trips.
    STL itself carries no standardized unit: consume ``coordinate_unit='m'``
    from the enclosing bundle. SVG is an independently scaled viewport view.
    """
    if not isinstance(contour, NozzleContour):
        raise DomainError('A NozzleContour record is required.', field='Export contour')
    inputs = ContourExportInputs(contour.contour, contour.throat_area_m2,
                                 contour.area_ratio, contour.stations, angular_segments)
    return ContourExportBundle(*_render(inputs), inputs=inputs)
