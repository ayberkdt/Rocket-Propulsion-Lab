"""Bounded engine maps with independent force and propellant-stream scales.

Unlike a single total-flow multiplier, this surface preserves oxidizer/fuel or
other named stream response independently.  That makes mixture-ratio drift and
component-specific feed sensitivities visible to a numerical propagator.
"""

from __future__ import annotations

import json
from bisect import bisect_right
from collections.abc import Mapping
from dataclasses import dataclass
from fractions import Fraction
from itertools import product
from math import isfinite

from rocket_propulsion.core.errors import DomainError, InputError
from .exchange_documents import _surface_document
from rocket_propulsion.core.validation import (
    name_value_pairs,
    read_json_document,
    real,
    real_mapping,
    real_tuple,
    records,
    strings,
    text,
)

from .performance_surface import (
    PerformanceAxis,
    PerformanceGridLineSurface,
    PerformanceValiditySurface,
)
from .interpolation import _MULTILINEAR_POLICY, _multilinear_terms
from .numeric_outputs import _finite_rational, _fraction
from .serialization import _require_numerical_policies, canonical_sha256

STREAM_PERFORMANCE_SURFACE_SCHEMA = "rocket_propulsion_stream_performance_surface_v1"
STREAM_PERFORMANCE_SURFACE_REFERENCE_IDS = (
    "NASA-20000031654",
    "NASA-20050206420",
    "NASA-20040000363",
    "NASA-STD-7009B",
)


@dataclass(frozen=True, slots=True)
class StreamPerformanceScalePoint:
    """One grid node with independent force and named-stream multipliers.

    References
    ----------
    NASA computer-controlled mixture-ratio study:
    https://ntrs.nasa.gov/citations/20050206420
    NASA aerospike parametric engine model:
    https://ntrs.nasa.gov/citations/20000031654
    """

    coordinates: tuple[float, ...]
    force_scale: float
    stream_flow_scales: tuple[tuple[str, float], ...]

    def __post_init__(self) -> None:
        object.__setattr__(
            self,
            "coordinates",
            real_tuple(self.coordinates, "Stream-performance coordinates"),
        )
        object.__setattr__(
            self,
            "force_scale",
            real(self.force_scale, "Stream-performance force scale"),
        )
        object.__setattr__(
            self,
            "stream_flow_scales",
            name_value_pairs(self.stream_flow_scales, "Stream-performance flow scales"),
        )
        if not self.coordinates:
            raise DomainError("Stream-performance coordinates must be finite and non-empty.")
        if self.force_scale <= 0.0:
            raise DomainError("Stream-performance force scale must be positive and finite.")
        if not self.stream_flow_scales:
            raise DomainError("A stream-performance point requires named stream scales.")
        names = [name for name, _ in self.stream_flow_scales]
        if any(not name.strip() for name in names) or len(names) != len(set(names)):
            raise DomainError("Stream-performance names must be unique and non-empty.")
        if any(not isfinite(value) or value <= 0.0 for _, value in self.stream_flow_scales):
            raise DomainError("Stream-performance flow scales must be positive and finite.")
        object.__setattr__(self, "stream_flow_scales", tuple(sorted(self.stream_flow_scales)))


@dataclass(frozen=True, slots=True)
class StreamPerformanceSurfaceEvaluation:
    """Interpolated force/stream scales and their analytic cell gradients.

    References
    ----------
    NASA computer-controlled mixture-ratio study:
    https://ntrs.nasa.gov/citations/20050206420
    NASA-STD-7009B: https://standards.nasa.gov/standard/NASA/NASA-STD-7009
    """

    conditions: tuple[tuple[str, float], ...]
    force_scale: float
    stream_flow_scales: tuple[tuple[str, float], ...]
    force_scale_gradient: tuple[tuple[str, float], ...]
    stream_flow_scale_gradients: tuple[tuple[str, tuple[tuple[str, float], ...]], ...]
    cell_lower_coordinates: tuple[float, ...]
    cell_upper_coordinates: tuple[float, ...]
    boundary_axes: tuple[str, ...]


@dataclass(frozen=True, slots=True)
class RectilinearStreamPerformanceSurface:
    """Complete N-D map of force and independent named-stream scales.

    Every Cartesian node is required exactly once and every node must expose
    the same stream identifiers.  Multilinear interpolation is performed on
    force and each stream independently. Values beyond the endpoint tolerance
    are refused; accepted round-off offsets are snapped to the endpoint.

    References
    ----------
    NASA computer-controlled mixture-ratio study:
    https://ntrs.nasa.gov/citations/20050206420
    NASA aerospike parametric engine model:
    https://ntrs.nasa.gov/citations/20000031654
    NASA transient-modeling overview:
    https://ntrs.nasa.gov/citations/20040000363
    NASA-STD-7009B: https://standards.nasa.gov/standard/NASA/NASA-STD-7009
    """

    axes: tuple[PerformanceAxis, ...]
    points: tuple[StreamPerformanceScalePoint, ...]
    source_id: str
    source_sha256: str
    reference_ids: tuple[str, ...] = STREAM_PERFORMANCE_SURFACE_REFERENCE_IDS
    schema: str = STREAM_PERFORMANCE_SURFACE_SCHEMA

    def __post_init__(self) -> None:
        object.__setattr__(
            self, "axes", records(self.axes, PerformanceAxis, "Stream-performance axes")
        )
        object.__setattr__(
            self,
            "points",
            records(self.points, StreamPerformanceScalePoint, "Stream-performance points"),
        )
        object.__setattr__(
            self,
            "reference_ids",
            strings(self.reference_ids, "Stream-performance references"),
        )
        text(self.source_id, "Stream-performance source identifier")
        text(self.source_sha256, "Stream-performance source hash")
        text(self.schema, "Stream-performance schema")
        if self.schema != STREAM_PERFORMANCE_SURFACE_SCHEMA:
            raise DomainError("Unsupported stream-performance surface schema.")
        if not self.axes:
            raise DomainError("A stream-performance surface requires at least one axis.")
        axis_names = [axis.name for axis in self.axes]
        if len(axis_names) != len(set(axis_names)):
            raise DomainError("Stream-performance axis names must be unique.")
        if not self.source_id.strip():
            raise DomainError("Stream-performance source identifier must not be empty.")
        if len(self.source_sha256) != 64 or any(
            character not in "0123456789abcdef" for character in self.source_sha256.lower()
        ):
            raise DomainError("Stream-performance source hash must be a SHA-256 digest.")
        if not self.reference_ids or any(not item.strip() for item in self.reference_ids):
            raise DomainError("Stream-performance references must not be empty.")
        if len(self.reference_ids) != len(set(self.reference_ids)):
            raise DomainError("Stream-performance references must be unique.")

        expected = set(product(*(axis.values for axis in self.axes)))
        coordinates = [point.coordinates for point in self.points]
        if any(len(item) != len(self.axes) for item in coordinates):
            raise DomainError("Stream-performance point dimensionality must match axes.")
        if len(coordinates) != len(set(coordinates)):
            raise DomainError("Stream-performance grid coordinates must be unique.")
        actual = set(coordinates)
        if actual != expected:
            raise DomainError(
                "Stream-performance surface must contain the complete Cartesian grid "
                f"(missing={len(expected - actual)}, extra={len(actual - expected)})."
            )
        stream_sets = {tuple(name for name, _ in point.stream_flow_scales) for point in self.points}
        if len(stream_sets) != 1:
            raise DomainError("Every stream-performance node must expose identical streams.")

    @property
    def stream_ids(self) -> tuple[str, ...]:
        """Return the stable sorted output-stream order.

        References
        ----------
        NASA computer-controlled mixture-ratio study:
        https://ntrs.nasa.gov/citations/20050206420
        """

        return tuple(name for name, _ in self.points[0].stream_flow_scales)

    @property
    def surface_hash(self) -> str:
        """Return content identity for grid, channels, and provenance.

        References
        ----------
        NASA-STD-7009B: https://standards.nasa.gov/standard/NASA/NASA-STD-7009
        """

        return canonical_sha256(self.to_dict(include_hash=False))

    def to_dict(self, *, include_hash: bool = True) -> dict[str, object]:
        """Return deterministic JSON-compatible surface content.

        References
        ----------
        NASA-STD-7009B: https://standards.nasa.gov/standard/NASA/NASA-STD-7009
        CCSDS-502.0-B-3: https://ccsds.org/Pubs/502x0b3e1.pdf
        """

        payload: dict[str, object] = {
            "schema": self.schema,
            "source_id": self.source_id,
            "source_sha256": self.source_sha256,
            "axes": [
                {"name": axis.name, "unit": axis.unit, "values": list(axis.values)}
                for axis in self.axes
            ],
            "points": [
                {
                    "coordinates": list(point.coordinates),
                    "force_scale": point.force_scale,
                    "stream_flow_scales": [list(item) for item in point.stream_flow_scales],
                }
                for point in sorted(self.points, key=lambda item: item.coordinates)
            ],
            "reference_ids": list(self.reference_ids),
            "interpolation_numerics": _MULTILINEAR_POLICY,
        }
        if include_hash:
            payload["surface_hash"] = canonical_sha256(payload)
        return payload

    def to_json(self, *, indent: int | None = 2) -> str:
        """Serialize the complete hash-bound stream surface.

        References
        ----------
        NASA-STD-7009B: https://standards.nasa.gov/standard/NASA/NASA-STD-7009
        """

        return json.dumps(
            self.to_dict(),
            allow_nan=False,
            ensure_ascii=False,
            indent=indent,
            sort_keys=True,
        )

    @classmethod
    def from_dict(cls, payload: Mapping[str, object]) -> RectilinearStreamPerformanceSurface:
        """Reopen a stream surface and reject semantic or hash changes.

        References
        ----------
        NASA-STD-7009B: https://standards.nasa.gov/standard/NASA/NASA-STD-7009
        """

        _surface_document(payload, streams=True)
        try:
            axes_payload = payload["axes"]
            points_payload = payload["points"]
            if not isinstance(axes_payload, list) or not isinstance(points_payload, list):
                raise TypeError
            surface = cls(
                schema=payload["schema"],
                source_id=payload["source_id"],
                source_sha256=payload["source_sha256"],
                reference_ids=payload["reference_ids"],
                axes=tuple(
                    PerformanceAxis(
                        name=item["name"],
                        unit=item["unit"],
                        values=item["values"],
                    )
                    for item in axes_payload
                ),
                points=tuple(
                    StreamPerformanceScalePoint(
                        coordinates=item["coordinates"],
                        force_scale=item["force_scale"],
                        stream_flow_scales=item["stream_flow_scales"],
                    )
                    for item in points_payload
                ),
            )
        except (KeyError, TypeError, ValueError, IndexError, DomainError) as error:
            raise InputError("Malformed stream-performance artifact.") from error
        if payload.get("surface_hash") != surface.surface_hash:
            raise InputError("Stream-performance hash does not match its contents.")
        _require_numerical_policies(payload, surface.to_dict(include_hash=False))
        return surface

    @classmethod
    def from_json(cls, text: str) -> RectilinearStreamPerformanceSurface:
        """Deserialize and integrity-check a stream surface.

        References
        ----------
        NASA-STD-7009B: https://standards.nasa.gov/standard/NASA/NASA-STD-7009
        """

        return read_json_document(text, "Stream-performance artifact", cls.from_dict)  # type: ignore[return-value]

    def validity_surfaces(self) -> tuple[PerformanceValiditySurface, ...]:
        """Return lower and upper event roots for every input coordinate.

        References
        ----------
        Orekit EventDetectorsProvider 13.1.5 API:
        https://www.orekit.org/static/apidocs/org/orekit/propagation/events/EventDetectorsProvider.html
        """

        return tuple(
            surface
            for axis in self.axes
            for surface in (
                PerformanceValiditySurface(
                    f"stream-performance:{axis.name}:minimum",
                    axis.name,
                    "lower",
                    axis.minimum,
                    axis.unit,
                ),
                PerformanceValiditySurface(
                    f"stream-performance:{axis.name}:maximum",
                    axis.name,
                    "upper",
                    axis.maximum,
                    axis.unit,
                ),
            )
        )

    def grid_line_surfaces(self) -> tuple[PerformanceGridLineSurface, ...]:
        """Return non-stopping roots at interior interpolation-grid lines.

        A state-transition or variational-equation consumer should restart at
        these derivative kinks even though the interpolated force and stream
        flows remain continuous.

        References
        ----------
        Orekit EventDetectorsProvider 13.1.5 API:
        https://www.orekit.org/static/apidocs/org/orekit/propagation/events/EventDetectorsProvider.html
        NASA-STD-7009B: https://standards.nasa.gov/standard/NASA/NASA-STD-7009
        """

        return tuple(
            PerformanceGridLineSurface(
                event_id=f"stream-performance:{axis.name}:grid:{index}",
                condition_name=axis.name,
                coordinate=value,
                unit=axis.unit,
            )
            for axis in self.axes
            for index, value in enumerate(axis.values[1:-1], start=1)
        )

    def evaluate(self, conditions: Mapping[str, float]) -> StreamPerformanceSurfaceEvaluation:
        """Interpolate force/stream scales and analytic gradients.

        References
        ----------
        NASA computer-controlled mixture-ratio study:
        https://ntrs.nasa.gov/citations/20050206420
        NASA-STD-7009B: https://standards.nasa.gov/standard/NASA/NASA-STD-7009
        """

        supplied = real_mapping(conditions, "Operating conditions")
        expected_names = {axis.name for axis in self.axes}
        supplied_names = set(supplied)
        if supplied_names != expected_names:
            raise DomainError(
                "Operating conditions must match stream-performance axes exactly "
                f"(missing={sorted(expected_names - supplied_names)}, "
                f"extra={sorted(supplied_names - expected_names)})."
            )

        normalized: dict[str, float] = {}
        brackets: list[tuple[int, int, Fraction, Fraction]] = []
        boundary_axes: list[str] = []
        for axis in self.axes:
            value = supplied[axis.name]
            tolerance = max(1e-12, max(abs(axis.minimum), abs(axis.maximum)) * 1e-12)
            if value < axis.minimum - tolerance or value > axis.maximum + tolerance:
                raise DomainError(
                    f"Operating condition {axis.name!r} lies outside "
                    f"[{axis.minimum}, {axis.maximum}] {axis.unit}; extrapolation refused."
                )
            value = min(axis.maximum, max(axis.minimum, value))
            normalized[axis.name] = value
            if any(abs(value - coordinate) <= tolerance for coordinate in axis.values):
                boundary_axes.append(axis.name)
            upper = min(bisect_right(axis.values, value), len(axis.values) - 1)
            lower = max(0, upper - 1)
            if upper == lower:
                lower = max(0, upper - 1)
            width = _fraction(axis.values[upper]) - _fraction(axis.values[lower])
            brackets.append((lower, upper, (_fraction(value) - _fraction(axis.values[lower])) / width, width))

        point_by_coordinate = {point.coordinates: point for point in self.points}
        force_scale = Fraction(0)
        stream_scales = {name: Fraction(0) for name in self.stream_ids}
        force_gradient = [Fraction(0)] * len(self.axes)
        stream_gradients = {name: [Fraction(0)] * len(self.axes) for name in self.stream_ids}
        for corner, weight, derivative_weights in _multilinear_terms(brackets):
            coordinate = tuple(
                axis.values[brackets[index][corner[index]]] for index, axis in enumerate(self.axes)
            )
            point = point_by_coordinate[coordinate]
            point_streams = {name: _fraction(value) for name, value in point.stream_flow_scales}
            point_force = _fraction(point.force_scale)
            force_scale += weight * point_force
            for stream_id in self.stream_ids:
                stream_scales[stream_id] += weight * point_streams[stream_id]
            for derivative_axis, derivative_weight in enumerate(derivative_weights):
                force_gradient[derivative_axis] += derivative_weight * point_force
                for stream_id in self.stream_ids:
                    stream_gradients[stream_id][derivative_axis] += (
                        derivative_weight * point_streams[stream_id]
                    )

        return StreamPerformanceSurfaceEvaluation(
            conditions=tuple((axis.name, normalized[axis.name]) for axis in self.axes),
            force_scale=_finite_rational(force_scale, "Multilinear force scale"),
            stream_flow_scales=tuple((name, _finite_rational(value, f"Multilinear stream scale: {name}"))
                                    for name, value in stream_scales.items()),
            force_scale_gradient=tuple(
                (axis.name, _finite_rational(force_gradient[index], f"Multilinear force gradient: {axis.name}"))
                for index, axis in enumerate(self.axes)
            ),
            stream_flow_scale_gradients=tuple(
                (
                    stream_id,
                    tuple(
                        (axis.name, _finite_rational(stream_gradients[stream_id][index],
                         f"Multilinear stream gradient: {stream_id} / {axis.name}"))
                        for index, axis in enumerate(self.axes)
                    ),
                )
                for stream_id in self.stream_ids
            ),
            cell_lower_coordinates=tuple(
                axis.values[brackets[index][0]] for index, axis in enumerate(self.axes)
            ),
            cell_upper_coordinates=tuple(
                axis.values[brackets[index][1]] for index, axis in enumerate(self.axes)
            ),
            boundary_axes=tuple(boundary_axes),
        )
