"""Bounded multidimensional propulsion performance surfaces.

The surface maps explicit operating conditions (for example supply pressure,
ambient pressure, throttle, mixture ratio, or available electrical power) to
force and tank-flow multipliers.  It uses a complete rectilinear grid,
multilinear interpolation, and analytic cell gradients.  Extrapolation is
never performed.
"""

from __future__ import annotations

import hashlib
import json
from bisect import bisect_right
from collections.abc import Mapping
from dataclasses import dataclass
from itertools import product
from math import isfinite, prod

from rocket_propulsion.core.errors import DomainError, InputError

PERFORMANCE_SURFACE_SCHEMA = "rocket_propulsion_performance_surface_v1"
PERFORMANCE_SURFACE_REFERENCE_IDS = (
    "NASA-20000031654",
    "NASA-STD-7009B",
    "NASA-20040000363",
)


def _canonical_hash(payload: object) -> str:
    encoded = json.dumps(
        payload,
        allow_nan=False,
        ensure_ascii=False,
        separators=(",", ":"),
        sort_keys=True,
    ).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()


@dataclass(frozen=True, slots=True)
class PerformanceAxis:
    """One strictly increasing operating-condition coordinate axis.

    References
    ----------
    NASA aerospike parametric engine model:
    https://ntrs.nasa.gov/citations/20000031654
    NASA-STD-7009B: https://standards.nasa.gov/standard/NASA/NASA-STD-7009
    """

    name: str
    unit: str
    values: tuple[float, ...]

    def __post_init__(self) -> None:
        if not self.name.strip() or not self.unit.strip():
            raise DomainError("Performance-axis name and unit must not be empty.")
        if len(self.values) < 2:
            raise DomainError("A performance axis requires at least two coordinates.")
        if any(not isfinite(value) for value in self.values):
            raise DomainError("Performance-axis coordinates must be finite.")
        if any(end <= start for start, end in zip(self.values, self.values[1:], strict=False)):
            raise DomainError("Performance-axis coordinates must be strictly increasing.")

    @property
    def minimum(self) -> float:
        """Return the inclusive lower validation-domain coordinate.

        References
        ----------
        NASA-STD-7009B: https://standards.nasa.gov/standard/NASA/NASA-STD-7009
        """

        return self.values[0]

    @property
    def maximum(self) -> float:
        """Return the inclusive upper validation-domain coordinate.

        References
        ----------
        NASA-STD-7009B: https://standards.nasa.gov/standard/NASA/NASA-STD-7009
        """

        return self.values[-1]


@dataclass(frozen=True, slots=True)
class PerformanceScalePoint:
    """One measured or validated grid node of force and flow multipliers.

    Coordinates follow the surface axis order.  Force and total tank-flow
    multipliers are independent because a measured operating map can contain
    changing effective specific impulse.

    References
    ----------
    NASA aerospike parametric engine model:
    https://ntrs.nasa.gov/citations/20000031654
    NASA transient-modeling overview:
    https://ntrs.nasa.gov/citations/20040000363
    """

    coordinates: tuple[float, ...]
    force_scale: float
    mass_flow_scale: float

    def __post_init__(self) -> None:
        if not self.coordinates or any(not isfinite(value) for value in self.coordinates):
            raise DomainError("Performance-point coordinates must be finite and non-empty.")
        for value, label in (
            (self.force_scale, "Performance force scale"),
            (self.mass_flow_scale, "Performance mass-flow scale"),
        ):
            if not isfinite(value) or value <= 0.0:
                raise DomainError(f"{label} must be finite and greater than zero.")


@dataclass(frozen=True, slots=True)
class PerformanceSurfaceEvaluation:
    """Interpolated scales, analytic gradients, and active grid cell.

    References
    ----------
    NASA aerospike parametric engine model:
    https://ntrs.nasa.gov/citations/20000031654
    NASA-STD-7009B: https://standards.nasa.gov/standard/NASA/NASA-STD-7009
    """

    conditions: tuple[tuple[str, float], ...]
    force_scale: float
    mass_flow_scale: float
    force_scale_gradient: tuple[tuple[str, float], ...]
    mass_flow_scale_gradient: tuple[tuple[str, float], ...]
    cell_lower_coordinates: tuple[float, ...]
    cell_upper_coordinates: tuple[float, ...]
    boundary_axes: tuple[str, ...]


@dataclass(frozen=True, slots=True)
class PerformanceValiditySurface:
    """Signed lower or upper validity-domain event for one condition.

    A positive value is inside the declared domain, zero is the boundary, and
    a negative value is outside.  A propagator can register both surfaces for
    each state-dependent condition and stop the burn before extrapolation.

    References
    ----------
    NASA-STD-7009B: https://standards.nasa.gov/standard/NASA/NASA-STD-7009
    Orekit EventDetector 13.1.5 API:
    https://www.orekit.org/static/apidocs/org/orekit/propagation/events/EventDetector.html
    """

    event_id: str
    condition_name: str
    bound: str
    boundary_value: float
    unit: str
    stops_burn: bool = True

    def __post_init__(self) -> None:
        if not self.event_id.strip() or not self.condition_name.strip() or not self.unit.strip():
            raise DomainError("Performance validity event names and unit must not be empty.")
        if self.bound not in {"lower", "upper"}:
            raise DomainError("Performance validity bound must be lower or upper.")
        if not isfinite(self.boundary_value):
            raise DomainError("Performance validity boundary must be finite.")

    def value(self, conditions: Mapping[str, float]) -> float:
        """Evaluate positive-inside signed distance in axis units."""

        if self.condition_name not in conditions:
            raise DomainError(f"Missing operating condition {self.condition_name!r}.")
        try:
            value = float(conditions[self.condition_name])
        except (TypeError, ValueError) as error:
            raise DomainError(
                f"Operating condition {self.condition_name!r} must be numeric."
            ) from error
        if not isfinite(value):
            raise DomainError(f"Operating condition {self.condition_name!r} must be finite.")
        if self.bound == "lower":
            return value - self.boundary_value
        return self.boundary_value - value


@dataclass(frozen=True, slots=True)
class RectilinearPerformanceSurface:
    """Complete non-extrapolating N-D force and tank-flow scale map.

    Every Cartesian grid node must be supplied exactly once.  Evaluation uses
    tensor-product multilinear interpolation; analytic gradients are constant
    with respect to each coordinate only in the one-dimensional case and vary
    multilinearly with the other coordinates in higher dimensions.  The
    surface records a source hash and stable reference identifiers.

    This is a validated-domain interpolation mechanism, not a fitted physical
    law.  The input nodes must come from applicable test data or a separately
    verified model.

    References
    ----------
    NASA aerospike parametric engine model:
    https://ntrs.nasa.gov/citations/20000031654
    NASA-STD-7009B: https://standards.nasa.gov/standard/NASA/NASA-STD-7009
    NASA transient-modeling overview:
    https://ntrs.nasa.gov/citations/20040000363
    """

    axes: tuple[PerformanceAxis, ...]
    points: tuple[PerformanceScalePoint, ...]
    source_id: str
    source_sha256: str
    reference_ids: tuple[str, ...] = PERFORMANCE_SURFACE_REFERENCE_IDS
    schema: str = PERFORMANCE_SURFACE_SCHEMA

    def __post_init__(self) -> None:
        if self.schema != PERFORMANCE_SURFACE_SCHEMA:
            raise DomainError("Unsupported performance-surface schema.")
        if not self.axes:
            raise DomainError("A performance surface requires at least one axis.")
        names = [axis.name for axis in self.axes]
        if len(names) != len(set(names)):
            raise DomainError("Performance-axis names must be unique.")
        if not self.source_id.strip():
            raise DomainError("Performance-surface source identifier must not be empty.")
        if len(self.source_sha256) != 64 or any(
            character not in "0123456789abcdef"
            for character in self.source_sha256.lower()
        ):
            raise DomainError("Performance-surface source hash must be a SHA-256 digest.")
        if not self.reference_ids or any(not item.strip() for item in self.reference_ids):
            raise DomainError("Performance-surface references must not be empty.")
        if len(self.reference_ids) != len(set(self.reference_ids)):
            raise DomainError("Performance-surface references must be unique.")

        expected_coordinates = set(product(*(axis.values for axis in self.axes)))
        actual_coordinates = [point.coordinates for point in self.points]
        if any(len(coordinates) != len(self.axes) for coordinates in actual_coordinates):
            raise DomainError("Performance-point dimensionality must match the axes.")
        if len(actual_coordinates) != len(set(actual_coordinates)):
            raise DomainError("Performance grid coordinates must be unique.")
        actual_set = set(actual_coordinates)
        if actual_set != expected_coordinates:
            missing = len(expected_coordinates - actual_set)
            extra = len(actual_set - expected_coordinates)
            raise DomainError(
                "Performance surface must contain the complete Cartesian grid "
                f"(missing={missing}, extra={extra})."
            )

    @property
    def surface_hash(self) -> str:
        """Return content identity for nodes, units, domain, and provenance.

        References
        ----------
        NASA-STD-7009B: https://standards.nasa.gov/standard/NASA/NASA-STD-7009
        """

        return _canonical_hash(self.to_dict(include_hash=False))

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
                    "mass_flow_scale": point.mass_flow_scale,
                }
                for point in sorted(self.points, key=lambda item: item.coordinates)
            ],
            "reference_ids": list(self.reference_ids),
        }
        if include_hash:
            payload["surface_hash"] = _canonical_hash(payload)
        return payload

    def to_json(self, *, indent: int | None = 2) -> str:
        """Serialize the complete hash-bound surface.

        References
        ----------
        NASA-STD-7009B: https://standards.nasa.gov/standard/NASA/NASA-STD-7009
        CCSDS-502.0-B-3: https://ccsds.org/Pubs/502x0b3e1.pdf
        """

        return json.dumps(
            self.to_dict(),
            allow_nan=False,
            ensure_ascii=False,
            indent=indent,
            sort_keys=True,
        )

    @classmethod
    def from_dict(cls, payload: Mapping[str, object]) -> RectilinearPerformanceSurface:
        """Reopen a surface and reject malformed or hash-altered content.

        References
        ----------
        NASA-STD-7009B: https://standards.nasa.gov/standard/NASA/NASA-STD-7009
        CCSDS-502.0-B-3: https://ccsds.org/Pubs/502x0b3e1.pdf
        """

        if not isinstance(payload, Mapping):
            raise InputError("Performance-surface artifact root must be an object.")
        try:
            axes_payload = payload["axes"]
            points_payload = payload["points"]
            if not isinstance(axes_payload, list) or not isinstance(points_payload, list):
                raise TypeError
            surface = cls(
                schema=str(payload["schema"]),
                source_id=str(payload["source_id"]),
                source_sha256=str(payload["source_sha256"]),
                reference_ids=tuple(str(item) for item in payload["reference_ids"]),
                axes=tuple(
                    PerformanceAxis(
                        name=str(item["name"]),
                        unit=str(item["unit"]),
                        values=tuple(float(value) for value in item["values"]),
                    )
                    for item in axes_payload
                ),
                points=tuple(
                    PerformanceScalePoint(
                        coordinates=tuple(float(value) for value in item["coordinates"]),
                        force_scale=float(item["force_scale"]),
                        mass_flow_scale=float(item["mass_flow_scale"]),
                    )
                    for item in points_payload
                ),
            )
        except (KeyError, TypeError, ValueError, IndexError, DomainError) as error:
            raise InputError("Malformed performance-surface artifact.") from error
        if payload.get("surface_hash") != surface.surface_hash:
            raise InputError("Performance-surface hash does not match its contents.")
        return surface

    @classmethod
    def from_json(cls, text: str) -> RectilinearPerformanceSurface:
        """Deserialize and integrity-check a JSON surface artifact.

        References
        ----------
        NASA-STD-7009B: https://standards.nasa.gov/standard/NASA/NASA-STD-7009
        CCSDS-502.0-B-3: https://ccsds.org/Pubs/502x0b3e1.pdf
        """

        try:
            payload = json.loads(text)
        except json.JSONDecodeError as error:
            raise InputError("Performance-surface artifact is not valid JSON.") from error
        if not isinstance(payload, dict):
            raise InputError("Performance-surface artifact root must be an object.")
        return cls.from_dict(payload)

    def validity_surfaces(self) -> tuple[PerformanceValiditySurface, ...]:
        """Return lower and upper event surfaces for every input coordinate.

        References
        ----------
        NASA-STD-7009B: https://standards.nasa.gov/standard/NASA/NASA-STD-7009
        Orekit EventDetectorsProvider 13.1.5 API:
        https://www.orekit.org/static/apidocs/org/orekit/propagation/events/EventDetectorsProvider.html
        """

        return tuple(
            surface
            for axis in self.axes
            for surface in (
                PerformanceValiditySurface(
                    f"performance:{axis.name}:minimum",
                    axis.name,
                    "lower",
                    axis.minimum,
                    axis.unit,
                ),
                PerformanceValiditySurface(
                    f"performance:{axis.name}:maximum",
                    axis.name,
                    "upper",
                    axis.maximum,
                    axis.unit,
                ),
            )
        )

    def evaluate(self, conditions: Mapping[str, float]) -> PerformanceSurfaceEvaluation:
        """Interpolate scales and analytic gradients without extrapolation.

        Conditions must match the axis names exactly.  A floating-point value
        within ``1e-12`` of a domain endpoint (relative to the axis magnitude)
        is snapped to that endpoint; values farther outside fail closed.

        References
        ----------
        NASA aerospike parametric engine model:
        https://ntrs.nasa.gov/citations/20000031654
        NASA-STD-7009B: https://standards.nasa.gov/standard/NASA/NASA-STD-7009
        """

        expected_names = {axis.name for axis in self.axes}
        supplied_names = set(conditions)
        if supplied_names != expected_names:
            missing = expected_names - supplied_names
            extra = supplied_names - expected_names
            raise DomainError(
                "Operating conditions must match surface axes exactly "
                f"(missing={sorted(missing)}, extra={sorted(extra)})."
            )

        normalized: dict[str, float] = {}
        brackets: list[tuple[int, int, float, float]] = []
        boundary_axes: list[str] = []
        for axis in self.axes:
            try:
                value = float(conditions[axis.name])
            except (TypeError, ValueError) as error:
                raise DomainError(
                    f"Operating condition {axis.name!r} must be numeric."
                ) from error
            if not isfinite(value):
                raise DomainError(f"Operating condition {axis.name!r} must be finite.")
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
            width = axis.values[upper] - axis.values[lower]
            fraction = (value - axis.values[lower]) / width
            brackets.append((lower, upper, fraction, width))

        point_by_coordinate = {point.coordinates: point for point in self.points}
        force_scale = 0.0
        flow_scale = 0.0
        force_gradient = [0.0] * len(self.axes)
        flow_gradient = [0.0] * len(self.axes)

        for corner in product((0, 1), repeat=len(self.axes)):
            coordinate = tuple(
                axis.values[brackets[index][corner[index]]]
                for index, axis in enumerate(self.axes)
            )
            point = point_by_coordinate[coordinate]
            factors = tuple(
                fraction if corner[index] else 1.0 - fraction
                for index, (_, _, fraction, _) in enumerate(brackets)
            )
            weight = prod(factors)
            force_scale += weight * point.force_scale
            flow_scale += weight * point.mass_flow_scale
            for derivative_axis, (_, _, _, width) in enumerate(brackets):
                sign = 1.0 if corner[derivative_axis] else -1.0
                other_weight = prod(
                    factor
                    for index, factor in enumerate(factors)
                    if index != derivative_axis
                )
                derivative_weight = sign * other_weight / width
                force_gradient[derivative_axis] += derivative_weight * point.force_scale
                flow_gradient[derivative_axis] += derivative_weight * point.mass_flow_scale

        return PerformanceSurfaceEvaluation(
            conditions=tuple((axis.name, normalized[axis.name]) for axis in self.axes),
            force_scale=force_scale,
            mass_flow_scale=flow_scale,
            force_scale_gradient=tuple(
                (axis.name, force_gradient[index])
                for index, axis in enumerate(self.axes)
            ),
            mass_flow_scale_gradient=tuple(
                (axis.name, flow_gradient[index])
                for index, axis in enumerate(self.axes)
            ),
            cell_lower_coordinates=tuple(
                axis.values[brackets[index][0]] for index, axis in enumerate(self.axes)
            ),
            cell_upper_coordinates=tuple(
                axis.values[brackets[index][1]] for index, axis in enumerate(self.axes)
            ),
            boundary_axes=tuple(boundary_axes),
        )

