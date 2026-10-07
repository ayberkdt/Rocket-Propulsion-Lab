"""Bounded multidimensional propulsion performance surfaces.

The surface maps explicit operating conditions (for example supply pressure,
ambient pressure, throttle, mixture ratio, or available electrical power) to
force and tank-flow multipliers.  It uses a complete rectilinear grid,
multilinear interpolation, and analytic cell gradients.  Extrapolation is
never performed.
"""

from __future__ import annotations

import json
from bisect import bisect_left, bisect_right
from collections.abc import Mapping
from dataclasses import dataclass, field
from itertools import combinations, product
from math import isfinite, prod
from types import MappingProxyType

from rocket_propulsion.core.errors import DomainError, InputError
from rocket_propulsion.core.validation import (
    read_json_document,
    real,
    real_tuple,
    records,
    strings,
    text,
)

from .interpolation import pchip_slopes
from .serialization import canonical_sha256

PERFORMANCE_SURFACE_SCHEMA = "rocket_propulsion_performance_surface_v1"
INTERPOLATION_MODES = ("multilinear", "cubic_hermite")
PERFORMANCE_SURFACE_REFERENCE_IDS = (
    "NASA-20000031654",
    "NASA-STD-7009B",
    "NASA-20040000363",
)


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
        text(self.name, "Performance-axis name")
        text(self.unit, "Performance-axis unit")
        object.__setattr__(self, "values", real_tuple(self.values, "Performance-axis values"))
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
        object.__setattr__(
            self, "coordinates", real_tuple(self.coordinates, "Performance-point coordinates")
        )
        object.__setattr__(self, "force_scale", real(self.force_scale, "Performance force scale"))
        object.__setattr__(
            self, "mass_flow_scale", real(self.mass_flow_scale, "Performance mass-flow scale")
        )
        if not self.coordinates:
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

        if not isinstance(conditions, Mapping):
            raise DomainError("Operating conditions must be a mapping.")
        if self.condition_name not in conditions:
            raise DomainError(f"Missing operating condition {self.condition_name!r}.")
        value = real(conditions[self.condition_name], f"Operating condition {self.condition_name!r}")
        if self.bound == "lower":
            return value - self.boundary_value
        return self.boundary_value - value


@dataclass(frozen=True, slots=True)
class PerformanceGridLineSurface:
    """Signed non-stopping event on one interior grid coordinate.

    Multilinear interpolation is continuous, but the derivative along an axis
    jumps when that condition crosses an interior grid coordinate.  A
    variable-step integrator that steps across such a kink without restarting
    loses its error estimate, and gradient-based optimizers see a
    discontinuous Jacobian.  Registering these surfaces lets the consumer
    land on each kink and restart, exactly as it does for profile boundaries.
    The value is ``condition - coordinate``; it never inhibits the burn.

    Exactly on the coordinate, :meth:`RectilinearPerformanceSurface.evaluate`
    returns the gradient of the upper cell (the right-hand limit along that
    axis), except at the axis maximum, where only the lower cell exists.

    References
    ----------
    Hairer, Norsett and Wanner, Solving Ordinary Differential Equations I,
    section II.6 (discontinuities and dense-output event location).
    Orekit EventDetector 13.1.5 API:
    https://www.orekit.org/static/apidocs/org/orekit/propagation/events/EventDetector.html
    """

    event_id: str
    condition_name: str
    coordinate: float
    unit: str
    stops_burn: bool = False

    def __post_init__(self) -> None:
        if not self.event_id.strip() or not self.condition_name.strip() or not self.unit.strip():
            raise DomainError("Performance grid-line event names and unit must not be empty.")
        if not isfinite(self.coordinate):
            raise DomainError("Performance grid-line coordinate must be finite.")
        if self.stops_burn:
            raise DomainError("Performance grid-line events never stop the burn.")

    def value(self, conditions: Mapping[str, float]) -> float:
        """Evaluate the signed distance ``condition - coordinate`` in axis units."""

        if not isinstance(conditions, Mapping):
            raise DomainError("Operating conditions must be a mapping.")
        if self.condition_name not in conditions:
            raise DomainError(f"Missing operating condition {self.condition_name!r}.")
        value = real(conditions[self.condition_name], f"Operating condition {self.condition_name!r}")
        return value - self.coordinate


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
    interpolation: str = "multilinear"
    _point_index: Mapping[tuple[float, ...], PerformanceScalePoint] = field(
        init=False, repr=False, compare=False
    )
    _hermite: Mapping[tuple[int, ...], Mapping[tuple[int, ...], tuple[float, float]]] = field(
        init=False, repr=False, compare=False
    )

    def __post_init__(self) -> None:
        if self.schema != PERFORMANCE_SURFACE_SCHEMA:
            raise DomainError("Unsupported performance-surface schema.")
        for name, kind in (("axes", PerformanceAxis), ("points", PerformanceScalePoint)):
            object.__setattr__(
                self, name, records(getattr(self, name), kind, f"Performance-surface {name}")
            )
        text(self.source_id, "Performance-surface source identifier")
        text(self.source_sha256, "Performance-surface source hash")
        object.__setattr__(
            self,
            "reference_ids",
            strings(self.reference_ids, "Performance-surface references", allow_blank=True),
        )
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
        if self.interpolation not in INTERPOLATION_MODES:
            raise DomainError(f"Surface interpolation must be one of {INTERPOLATION_MODES}.")
        # Built once: evaluation is called at every propagator stage.
        object.__setattr__(
            self,
            "_point_index",
            MappingProxyType({point.coordinates: point for point in self.points}),
        )
        object.__setattr__(
            self,
            "_hermite",
            self._hermite_derivatives() if self.interpolation == "cubic_hermite" else {},
        )

    def _hermite_derivatives(
        self,
    ) -> Mapping[tuple[int, ...], Mapping[tuple[int, ...], tuple[float, float]]]:
        """Nodal values and mixed partials for tensor-product cubic Hermite.

        ``D[S]`` is the mixed partial over the axes in ``S`` (a sorted tuple).
        It is built by applying the PCHIP slope operator along the largest
        axis of ``S`` to ``D[S minus that axis]``, one grid line at a time.
        """

        dimension = len(self.axes)
        shape = [len(axis.values) for axis in self.axes]
        base: dict[tuple[int, ...], tuple[float, float]] = {}
        for index in product(*(range(size) for size in shape)):
            coordinate = tuple(self.axes[a].values[i] for a, i in enumerate(index))
            point = self._point_index[coordinate]
            base[index] = (point.force_scale, point.mass_flow_scale)
        derivatives: dict[tuple[int, ...], Mapping[tuple[int, ...], tuple[float, float]]] = {
            (): MappingProxyType(base)
        }
        subsets = sorted(
            (subset for size in range(1, dimension + 1) for subset in combinations(range(dimension), size)),
            key=len,
        )
        for subset in subsets:
            axis_index = subset[-1]
            source = derivatives[subset[:-1]]
            values = self.axes[axis_index].values
            result: dict[tuple[int, ...], tuple[float, float]] = {}
            others = [range(size) for a, size in enumerate(shape) if a != axis_index]
            for fixed in product(*others):
                line = []
                for i in range(shape[axis_index]):
                    index = list(fixed)
                    index.insert(axis_index, i)
                    line.append(tuple(index))
                force = pchip_slopes(values, tuple(source[item][0] for item in line))
                flow = pchip_slopes(values, tuple(source[item][1] for item in line))
                for position, item in enumerate(line):
                    result[item] = (force[position], flow[position])
            derivatives[subset] = MappingProxyType(result)
        return MappingProxyType(derivatives)

    @property
    def surface_hash(self) -> str:
        """Return content identity for nodes, units, domain, and provenance.

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
                    "mass_flow_scale": point.mass_flow_scale,
                }
                for point in sorted(self.points, key=lambda item: item.coordinates)
            ],
            "reference_ids": list(self.reference_ids),
        }
        if self.interpolation != "multilinear":
            # Present only when used, so multilinear surface identities are unchanged.
            payload["interpolation"] = self.interpolation
        if include_hash:
            payload["surface_hash"] = canonical_sha256(payload)
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
                schema=payload["schema"],
                source_id=payload["source_id"],
                source_sha256=payload["source_sha256"],
                reference_ids=payload["reference_ids"],
                interpolation=payload.get("interpolation", "multilinear"),
                axes=tuple(
                    PerformanceAxis(
                        name=item["name"],
                        unit=item["unit"],
                        values=item["values"],
                    )
                    for item in axes_payload
                ),
                points=tuple(
                    PerformanceScalePoint(
                        coordinates=item["coordinates"],
                        force_scale=item["force_scale"],
                        mass_flow_scale=item["mass_flow_scale"],
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

        return read_json_document(text, "Performance-surface artifact", cls.from_dict)

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

    def grid_line_surfaces(self) -> tuple[PerformanceGridLineSurface, ...]:
        """Return one non-stopping kink event per interior axis coordinate.

        Domain endpoints are excluded; they are covered by the stopping
        :meth:`validity_surfaces`.  A ``cubic_hermite`` surface is C1, so it
        has no derivative kinks and returns no surfaces.

        References
        ----------
        Orekit EventDetectorsProvider 13.1.5 API:
        https://www.orekit.org/static/apidocs/org/orekit/propagation/events/EventDetectorsProvider.html
        """

        if self.interpolation == "cubic_hermite":
            return ()
        return tuple(
            PerformanceGridLineSurface(
                f"performance:{axis.name}:grid:{index}",
                axis.name,
                coordinate,
                axis.unit,
            )
            for axis in self.axes
            for index, coordinate in enumerate(axis.values[1:-1], start=1)
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

        if not isinstance(conditions, Mapping):
            raise DomainError("Operating conditions must be a mapping.")
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
            value = real(conditions[axis.name], f"Operating condition {axis.name!r}")
            tolerance = max(1e-12, max(abs(axis.minimum), abs(axis.maximum)) * 1e-12)
            if value < axis.minimum - tolerance or value > axis.maximum + tolerance:
                raise DomainError(
                    f"Operating condition {axis.name!r} lies outside "
                    f"[{axis.minimum}, {axis.maximum}] {axis.unit}; extrapolation refused."
                )
            value = min(axis.maximum, max(axis.minimum, value))
            normalized[axis.name] = value
            nearest = bisect_left(axis.values, value - tolerance)
            if nearest < len(axis.values) and axis.values[nearest] <= value + tolerance:
                boundary_axes.append(axis.name)
            upper = min(bisect_right(axis.values, value), len(axis.values) - 1)
            lower = max(0, upper - 1)
            if upper == lower:
                lower = max(0, upper - 1)
            width = axis.values[upper] - axis.values[lower]
            fraction = (value - axis.values[lower]) / width
            brackets.append((lower, upper, fraction, width))

        if self.interpolation == "cubic_hermite":
            force_scale, flow_scale, force_gradient, flow_gradient = self._evaluate_hermite(
                brackets
            )
            if force_scale <= 0.0 or flow_scale <= 0.0:
                raise DomainError(
                    "Cubic surface produced a non-positive scale inside the grid; "
                    "use multilinear interpolation for these data."
                )
            return self._evaluation(
                normalized, brackets, boundary_axes,
                force_scale, flow_scale, force_gradient, flow_gradient,
            )
        point_by_coordinate = self._point_index
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

        return self._evaluation(
            normalized, brackets, boundary_axes,
            force_scale, flow_scale, force_gradient, flow_gradient,
        )

    def _evaluate_hermite(
        self, brackets: list[tuple[int, int, float, float]]
    ) -> tuple[float, float, list[float], list[float]]:
        dimension = len(self.axes)
        force = 0.0
        flow = 0.0
        force_gradient = [0.0] * dimension
        flow_gradient = [0.0] * dimension
        basis_cache = []
        for _, _, t, width in brackets:
            t2, t3 = t * t, t * t * t
            # (value basis, its d/dx) for [corner 0 value, corner 1 value,
            # corner 0 slope, corner 1 slope]
            basis_cache.append(
                (
                    (2 * t3 - 3 * t2 + 1, (6 * t2 - 6 * t) / width),
                    (-2 * t3 + 3 * t2, (-6 * t2 + 6 * t) / width),
                    ((t3 - 2 * t2 + t) * width, 3 * t2 - 4 * t + 1),
                    ((t3 - t2) * width, 3 * t2 - 2 * t),
                )
            )
        for subset, table in self._hermite.items():
            members = set(subset)
            for corner in product((0, 1), repeat=dimension):
                index = tuple(brackets[a][corner[a]] for a in range(dimension))
                d_force, d_flow = table[index]
                values = []
                slopes = []
                for a in range(dimension):
                    value, slope = basis_cache[a][corner[a] + (2 if a in members else 0)]
                    values.append(value)
                    slopes.append(slope)
                weight = prod(values)
                force += weight * d_force
                flow += weight * d_flow
                for j in range(dimension):
                    partial = slopes[j] * prod(values[a] for a in range(dimension) if a != j)
                    force_gradient[j] += partial * d_force
                    flow_gradient[j] += partial * d_flow
        return force, flow, force_gradient, flow_gradient

    def _evaluation(
        self,
        normalized: dict[str, float],
        brackets: list[tuple[int, int, float, float]],
        boundary_axes: list[str],
        force_scale: float,
        flow_scale: float,
        force_gradient: list[float],
        flow_gradient: list[float],
    ) -> PerformanceSurfaceEvaluation:
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

