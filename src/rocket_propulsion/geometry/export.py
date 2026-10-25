"""Text exports of a sampled axisymmetric nozzle fluid domain.

The exported CSV station table, SVG side view, DXF closed polyline and ASCII
STL surface describe the closed fluid domain bounded by the sampled wall and
the inlet and exit planes. They are intended for visualization and data
exchange; they are neither a qualified wall design nor a manufacturing model.
Every bundle is regenerated from its retained inputs on construction, so its
text always matches those inputs. DXF output follows [1]_ and [2]_, ASCII STL
output [3]_ and [4]_, and contour names obey the XML 1.0 character rules [5]_.

References
----------
.. [1] AUTODESK-LWPOLYLINE. Autodesk, DXF lightweight-polyline subclass,
       vertex/count and closure group codes.
.. [2] AUTODESK-INSUNITS. Autodesk, drawing insertion unit convention.
.. [3] LOC-STL-FAMILY. Library of Congress, STL format family: oriented
       triangular surface and format limitations.
.. [4] LOC-STL-ASCII. Library of Congress, ASCII STL records and metadata
       limitations.
.. [5] W3C-XML10. W3C, "Extensible Markup Language (XML) 1.0", fifth edition,
       permitted characters and literal character data.
"""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass, fields
from fractions import Fraction
from functools import lru_cache
from itertools import pairwise
from xml.sax.saxutils import escape

from rocket_propulsion.core.errors import DomainError
from rocket_propulsion.core.validation import _ordered_items, integer, real, text

from .export_numerics import _check_geometry, _ring, _unit_normal
from .nozzle_contour import NozzleContour, NozzleStation

_SERIALIZATION_POLICY = 'binary64_round_trip_17_digits_exact_vertex_cross_normals_v1'
_MESH_ROLE = 'closed_fluid_domain_visualization'


def _owned(value, kind, label):
    """Return a fresh ``kind`` instance copied from a mapping or a ``kind`` record.

    Raises
    ------
    DomainError
        If ``value`` is neither a mapping nor a ``kind`` instance, or if
        construction from it fails with a type, value, key, index, attribute
        or overflow error. Other exceptions propagate unchanged.
    """
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
    """Validated, owned copy of the contour quantities consumed by the export.

    The record does not verify that the stations come from a particular flow
    solver; it checks only their internal geometric consistency.

    Attributes
    ----------
    contour : str
        Contour name; it must not contain characters forbidden in XML 1.0
        [1]_.
    throat_area_m2 : float
        Positive throat area in m^2.
    area_ratio : float
        Exit-to-throat area ratio, at least one.
    stations : tuple of NozzleStation
        At least two stations with strictly increasing ``x_m``, positive
        radius, area, area ratio and stagnation ratios, and non-negative Mach.
    angular_segments : int, default 32
        Number of STL facets around the circumference, between 8 and 128.

    Raises
    ------
    DomainError
        If any field violates the constraints above or the redundant station
        geometry disagrees beyond binary64 representation cells.

    References
    ----------
    .. [1] W3C-XML10. W3C, "Extensible Markup Language (XML) 1.0", fifth edition,
           permitted characters and literal character data.
    """
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
    """Export texts for visualization and exchange, not a nozzle wall design.

    Construction re-renders all four texts from ``inputs`` and rejects any
    mismatch, so a bundle always reproduces from its retained inputs.

    Attributes
    ----------
    csv : str
        Station table with 17 significant digits per value.
    svg : str
        Independently scaled side view with model metadata.
    dxf : str
        Closed ``LWPOLYLINE`` of the upper and mirrored lower wall in metres
        (``$INSUNITS`` = 6) [1]_, [2]_.
    stl : str
        ASCII STL surface of the closed fluid domain [3]_, [4]_.
    inputs : ContourExportInputs
        Retained source quantities.
    coordinate_unit : str, default "m"
        Coordinate unit of all texts; STL itself carries no unit.
    mesh_role : str
        Fixed statement that the mesh is a closed fluid-domain visualization.
    manufacturing_qualified, physical_model_qualified : bool
        Always ``False``.
    serialization_policy : str
        Identifier of the numeric serialization policy.

    Raises
    ------
    DomainError
        If a policy field differs from its fixed value, a qualification flag
        is set, or a text does not reproduce from ``inputs``.

    References
    ----------
    .. [1] AUTODESK-LWPOLYLINE. Autodesk, DXF lightweight-polyline subclass,
           vertex/count and closure group codes.
    .. [2] AUTODESK-INSUNITS. Autodesk, drawing insertion unit convention.
    .. [3] LOC-STL-FAMILY. Library of Congress, STL format family: oriented
           triangular surface and format limitations.
    .. [4] LOC-STL-ASCII. Library of Congress, ASCII STL records and metadata
           limitations.
    """

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
    """Return the station table as comma-separated text with a header row.

    Values use 17 significant digits so every binary64 value round-trips
    exactly; every line, including the last, ends with a newline.
    """
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
    """Return a 1200 x 500 SVG side view of the contour and its centreline.

    The view is scaled independently in x and r with exact rational
    arithmetic and is not a dimensioned drawing. The contour name is
    XML-escaped, including carriage returns, before it is written into the
    metadata element [1]_.

    References
    ----------
    .. [1] W3C-XML10. W3C, "Extensible Markup Language (XML) 1.0", fifth edition,
           permitted characters and literal character data.
    """
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
    """Return a minimal DXF drawing with one closed lightweight polyline.

    The header sets ``$INSUNITS`` to 6 (metres) [2]_. The ``LWPOLYLINE`` entity
    lists the upper wall followed by the mirrored lower wall, with group code
    90 giving the vertex count and group code 70 set to 1 (closed) [1]_.
    Coordinates use 17 significant digits.

    References
    ----------
    .. [1] AUTODESK-LWPOLYLINE. Autodesk, DXF lightweight-polyline subclass,
           vertex/count and closure group codes.
    .. [2] AUTODESK-INSUNITS. Autodesk, drawing insertion unit convention.
    """
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
    """Return one ASCII STL ``facet`` record for the triangle ``(a, b, c)``.

    The facet normal is computed from the serialized vertices and the vertex
    order follows the right-hand rule [1]_.

    References
    ----------
    .. [1] LOC-STL-ASCII. Library of Congress, ASCII STL records and metadata
           limitations.
    """

    unit = _unit_normal(a, b, c)
    lines = [f"  facet normal {unit[0]:.17g} {unit[1]:.17g} {unit[2]:.17g}", "    outer loop"]
    lines.extend(f"      vertex {x:.17g} {y:.17g} {z:.17g}" for x, y, z in (a, b, c))
    lines.extend(["    endloop", "  endfacet"])
    return "\n".join(lines)


def _stl(contour: ContourExportInputs, angular_segments: int) -> str:
    """Return a closed ASCII STL surface of the sampled fluid domain.

    Adjacent station rings are joined by two triangles per angular segment,
    and the inlet and exit planes are closed by triangle fans about the axis,
    with outward-facing normals. STL carries no unit; coordinates are in
    metres as declared by the enclosing bundle [1]_, [2]_.

    References
    ----------
    .. [1] LOC-STL-FAMILY. Library of Congress, STL format family: oriented
           triangular surface and format limitations.
    .. [2] LOC-STL-ASCII. Library of Congress, ASCII STL records and metadata
           limitations.
    """
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
    """Return ``(csv, svg, dxf, stl)`` for ``inputs``.

    A one-entry cache lets :func:`export_nozzle_contour` and the bundle's
    reproduction check share one rendering while bounding retained memory.
    """
    return _csv(inputs), _svg(inputs), _dxf(inputs), _stl(inputs, inputs.angular_segments)


def export_nozzle_contour(
    contour: NozzleContour, *, angular_segments: int = 32
) -> ContourExportBundle:
    """Export a sampled closed fluid domain as CSV, SVG, DXF and ASCII STL text.

    Parameters
    ----------
    contour : NozzleContour
        Contour whose stations define the wall.
    angular_segments : int, default 32
        Number of STL facets around the circumference, between 8 and 128.

    Returns
    -------
    ContourExportBundle
        The four texts with their retained inputs and fixed policy fields.

    Raises
    ------
    DomainError
        If ``contour`` is not a :class:`NozzleContour` or its stations fail
        the checks of :class:`ContourExportInputs`.

    Notes
    -----
    The export carries no wall, flow or manufacturing qualification. CSV, DXF
    and STL use 17 significant digits for binary64 coordinate round trips.
    STL itself carries no standardized unit [3]_: consume
    ``coordinate_unit='m'`` from the enclosing bundle. DXF declares metres
    through ``$INSUNITS`` [2]_ and stores the outline as a closed lightweight
    polyline [1]_. SVG is an independently scaled viewport view.

    References
    ----------
    .. [1] AUTODESK-LWPOLYLINE. Autodesk, DXF lightweight-polyline subclass,
           vertex/count and closure group codes.
    .. [2] AUTODESK-INSUNITS. Autodesk, drawing insertion unit convention.
    .. [3] LOC-STL-FAMILY. Library of Congress, STL format family: oriented
           triangular surface and format limitations.
    """
    if not isinstance(contour, NozzleContour):
        raise DomainError('A NozzleContour record is required.', field='Export contour')
    inputs = ContourExportInputs(contour.contour, contour.throat_area_m2,
                                 contour.area_ratio, contour.stations, angular_segments)
    return ContourExportBundle(*_render(inputs), inputs=inputs)
