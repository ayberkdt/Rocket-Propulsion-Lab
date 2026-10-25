"""Bounded engine maps with independent force and propellant-stream scales.

Unlike a single total-flow multiplier, this surface carries the oxidizer, fuel
or other named stream responses independently, so mixture-ratio drift [1]_ and
component-specific feed sensitivities remain visible to a numerical
propagator.  Engine performance is tabulated over operating parameters in the
manner of the trajectory-facing tables of [2]_ and interpolated multilinearly
inside the grid [3]_; extrapolation is refused.

References
----------
.. [1] NASA, "Simulation Analysis of Computer-Controlled Pressurization for
       Mixture Ratio Control". https://ntrs.nasa.gov/citations/20050206420
.. [2] NASA, "Parametric Model of an Aerospike Rocket Engine".
       https://ntrs.nasa.gov/citations/20000031654
.. [3] SCIPY-1.14.1-REGULAR-GRID. SciPy 1.14.1, "RegularGridInterpolator",
       rectilinear grids and linear interpolation.
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

from .exchange_documents import _surface_document
from .interpolation import _MULTILINEAR_POLICY, _multilinear_terms
from .numeric_outputs import _finite_rational, _fraction
from .performance_surface import (
    PerformanceAxis,
    PerformanceGridLineSurface,
    PerformanceValiditySurface,
)
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

    Parameters
    ----------
    coordinates : tuple of float
        Non-empty node coordinates, one per surface axis.
    force_scale : float
        Positive multiplier on rated delivered thrust.
    stream_flow_scales : tuple of (str, float)
        Positive multiplier per uniquely named propellant stream; stored
        sorted by name.

    Raises
    ------
    DomainError
        If a field has the wrong type, the coordinates are empty, a scale is
        non-finite or not strictly positive, or a stream name is blank or
        repeated.

    Notes
    -----
    Separate stream multipliers let a map express mixture-ratio shifts [1]_
    as well as the power-level dependence of a parametric engine table [2]_.

    References
    ----------
    .. [1] NASA, "Simulation Analysis of Computer-Controlled Pressurization for
           Mixture Ratio Control". https://ntrs.nasa.gov/citations/20050206420
    .. [2] NASA, "Parametric Model of an Aerospike Rocket Engine".
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
    """Interpolated force and stream scales with their analytic cell gradients.

    Independent per-stream scales represent a mixture-ratio response of the
    kind studied in [1]_; the reported cell and boundary axes record where in
    the declared domain the evaluation lies, for the model-use reporting of
    [2]_.

    Attributes
    ----------
    conditions : tuple of (str, float)
        Operating conditions in axis order, after endpoint snapping.
    force_scale : float
        Interpolated thrust multiplier.
    stream_flow_scales : tuple of (str, float)
        Interpolated multiplier per stream, sorted by stream name.
    force_scale_gradient : tuple of (str, float)
        Partial derivative of the thrust multiplier per axis, inside the cell.
    stream_flow_scale_gradients : tuple of (str, tuple of (str, float))
        Partial derivatives of each stream multiplier per axis.
    cell_lower_coordinates, cell_upper_coordinates : tuple of float
        Bounds of the interpolation cell that produced the values.
    boundary_axes : tuple of str
        Axes on which the condition lies on a grid line (within tolerance),
        where the gradient is one-sided.

    References
    ----------
    .. [1] NASA, "Simulation Analysis of Computer-Controlled Pressurization for
           Mixture Ratio Control". https://ntrs.nasa.gov/citations/20050206420
    .. [2] NASA-STD-7009B. NASA, "Standard for Models and Simulations", 2024.
           https://standards.nasa.gov/standard/NASA/NASA-STD-7009
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
    """Complete N-dimensional map of force and independent named-stream scales.

    Parameters
    ----------
    axes : tuple of PerformanceAxis
        At least one uniquely named, unit-labelled grid axis.
    points : tuple of StreamPerformanceScalePoint
        Exactly one node per point of the Cartesian product of the axes.
    source_id : str
        Non-blank identifier of the data source.
    source_sha256 : str
        SHA-256 hex digest of the source data.
    reference_ids : tuple of str, default STREAM_PERFORMANCE_SURFACE_REFERENCE_IDS
        Non-empty, unique source identifiers.
    schema : str, default STREAM_PERFORMANCE_SURFACE_SCHEMA
        Document schema; only the current schema is accepted.

    Raises
    ------
    DomainError
        If a field has the wrong type, or the schema, axes, provenance fields
        or grid violate the contract: missing or extra nodes, mismatched
        dimensionality, or nodes exposing different stream identifiers.

    Notes
    -----
    Every Cartesian node is required exactly once and every node exposes the
    same stream identifiers.  Force and each stream are interpolated
    multilinearly and independently [3]_, so mixture-ratio variation across
    the map [1]_ is preserved; the grid layout follows the parameterized
    engine tables of [2]_.  Conditions beyond a relative endpoint tolerance of
    ``1e-12`` are refused; accepted round-off offsets are snapped to the
    endpoint.  [4]_ and [5]_ give the transient-modelling and model-credibility
    context.

    References
    ----------
    .. [1] NASA, "Simulation Analysis of Computer-Controlled Pressurization for
           Mixture Ratio Control". https://ntrs.nasa.gov/citations/20050206420
    .. [2] NASA, "Parametric Model of an Aerospike Rocket Engine".
           https://ntrs.nasa.gov/citations/20000031654
    .. [3] SCIPY-1.14.1-REGULAR-GRID. SciPy 1.14.1, "RegularGridInterpolator",
           rectilinear grids and linear interpolation.
    .. [4] NASA, "Transient Mathematical Modeling for Liquid Rocket Engine
           Systems". https://ntrs.nasa.gov/citations/20040000363
    .. [5] NASA-STD-7009B. NASA, "Standard for Models and Simulations", 2024.
           https://standards.nasa.gov/standard/NASA/NASA-STD-7009
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
        """Stream identifiers in the stable, sorted output order [1]_.

        References
        ----------
        .. [1] NASA, "Simulation Analysis of Computer-Controlled Pressurization for
               Mixture Ratio Control". https://ntrs.nasa.gov/citations/20050206420
        """

        return tuple(name for name, _ in self.points[0].stream_flow_scales)

    @property
    def surface_hash(self) -> str:
        """SHA-256 content identity of the grid, channels and provenance.

        The hash covers the canonical document without its own hash field, so
        any change to values, units or provenance changes the identity used
        for traceable model reporting [1]_.

        References
        ----------
        .. [1] NASA-STD-7009B. NASA, "Standard for Models and Simulations", 2024.
               https://standards.nasa.gov/standard/NASA/NASA-STD-7009
        """

        return canonical_sha256(self.to_dict(include_hash=False))

    def to_dict(self, *, include_hash: bool = True) -> dict[str, object]:
        """Return the deterministic, JSON-compatible surface document.

        Parameters
        ----------
        include_hash : bool, default True
            Whether to add the ``surface_hash`` field.

        Returns
        -------
        dict
            Schema, provenance, unit-labelled axes, nodes sorted by
            coordinates, reference identifiers and interpolation policy.

        Notes
        -----
        The document is versioned and unit-explicit in the manner of [2]_ and
        supports the reporting requirements of [1]_.

        References
        ----------
        .. [1] NASA-STD-7009B. NASA, "Standard for Models and Simulations", 2024.
               https://standards.nasa.gov/standard/NASA/NASA-STD-7009
        .. [2] CCSDS 502.0-B-3. CCSDS, "Orbit Data Messages".
               https://ccsds.org/Pubs/502x0b3e1.pdf
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
        """Serialize the complete hash-bound stream surface as JSON text.

        Parameters
        ----------
        indent : int or None, default 2
            JSON indentation; ``None`` gives the compact form.

        Returns
        -------
        str
            Key-sorted JSON without NaN or infinity tokens [1]_.

        References
        ----------
        .. [1] IETF RFC 8259, "The JavaScript Object Notation (JSON) Data Interchange
               Format". https://www.rfc-editor.org/rfc/rfc8259
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

        Parameters
        ----------
        payload : Mapping
            Document produced by `to_dict`.

        Returns
        -------
        RectilinearStreamPerformanceSurface
            The reconstructed, validated surface.

        Raises
        ------
        InputError
            If the document is malformed, its ``surface_hash`` does not match
            its contents, or its numerical policies differ from the current
            ones.

        Notes
        -----
        Hash verification keeps a reopened map traceable to the reported
        model [1]_.

        References
        ----------
        .. [1] NASA-STD-7009B. NASA, "Standard for Models and Simulations", 2024.
               https://standards.nasa.gov/standard/NASA/NASA-STD-7009
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
        """Deserialize and integrity-check a stream surface from JSON text.

        Parameters
        ----------
        text : str
            JSON document produced by `to_json`.

        Returns
        -------
        RectilinearStreamPerformanceSurface
            The reconstructed, validated surface.

        Raises
        ------
        InputError
            If the text is not a conforming JSON document [1]_ or fails the
            checks of `from_dict`.

        References
        ----------
        .. [1] IETF RFC 8259, "The JavaScript Object Notation (JSON) Data Interchange
               Format". https://www.rfc-editor.org/rfc/rfc8259
        """

        return read_json_document(text, "Stream-performance artifact", cls.from_dict)  # type: ignore[return-value]

    def validity_surfaces(self) -> tuple[PerformanceValiditySurface, ...]:
        """Return lower and upper validity event surfaces for every axis.

        Returns
        -------
        tuple of PerformanceValiditySurface
            Two surfaces per axis, at its minimum and maximum, for a
            propagator's event detection [1]_.

        References
        ----------
        .. [1] Orekit 13.1.5, "PropulsionModel / EventDetector API", the
               propagation-time force, parameter-driver and event contract.
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
        """Return non-stopping event surfaces at interior grid lines.

        Returns
        -------
        tuple of PerformanceGridLineSurface
            One surface per interior grid value of every axis.

        Notes
        -----
        Multilinear interpolation [2]_ is continuous, but its gradient jumps
        across grid lines.  A state-transition or variational-equation consumer
        restarts integration at these derivative kinks, detected as propagator
        events [1]_, even though the interpolated force and stream flows remain
        continuous.

        References
        ----------
        .. [1] Orekit 13.1.5, "PropulsionModel / EventDetector API", the
               propagation-time force, parameter-driver and event contract.
        .. [2] SCIPY-1.14.1-REGULAR-GRID. SciPy 1.14.1, "RegularGridInterpolator",
               rectilinear grids and linear interpolation.
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
        """Interpolate force and stream scales with their analytic gradients.

        Parameters
        ----------
        conditions : Mapping of str to float
            One value per surface axis, keyed by axis name.

        Returns
        -------
        StreamPerformanceSurfaceEvaluation
            Interpolated multipliers, per-axis gradients and cell bounds.

        Raises
        ------
        DomainError
            If the condition names differ from the axis names, or a value lies
            outside its axis by more than the endpoint tolerance
            ``max(1e-12, 1e-12 * max(|min|, |max|))``.

        Notes
        -----
        On the cell bracketing the condition, each output is the multilinear
        blend ``sum_c w_c * v_c`` of its ``2**N`` corner values, where ``w_c``
        is the product of the per-axis linear weights [2]_; the gradients are
        the exact derivatives of the same polynomial.  Weights and sums are
        accumulated in exact rational arithmetic and rounded once on export.
        Independent stream multipliers preserve mixture-ratio variation [1]_.
        Values within the endpoint tolerance are snapped to the endpoint.

        References
        ----------
        .. [1] NASA, "Simulation Analysis of Computer-Controlled Pressurization for
               Mixture Ratio Control". https://ntrs.nasa.gov/citations/20050206420
        .. [2] SCIPY-1.14.1-REGULAR-GRID. SciPy 1.14.1, "RegularGridInterpolator",
               rectilinear grids and linear interpolation.
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
