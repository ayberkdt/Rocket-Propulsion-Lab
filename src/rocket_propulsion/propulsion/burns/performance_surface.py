"""Bounded multidimensional propulsion performance surfaces.

A surface maps explicit operating conditions (for example supply pressure,
ambient pressure, throttle, mixture ratio or available electrical power) to
force and tank-flow multipliers, in the manner of trajectory-facing engine
performance tables [1]_.  It uses a complete rectilinear grid with
multilinear interpolation [2]_ or a tensor-product monotone cubic Hermite
(PCHIP) interpolant [3]_, [4]_, together with analytic gradients.
Extrapolation is never performed.

References
----------
.. [1] NASA, "Parametric Model of an Aerospike Rocket Engine".
       https://ntrs.nasa.gov/citations/20000031654
.. [2] SCIPY-1.14.1-REGULAR-GRID. SciPy 1.14.1 documentation,
       ``scipy.interpolate.RegularGridInterpolator``: rectilinear grids and
       linear interpolation.
.. [3] F. N. Fritsch and R. E. Carlson, "Monotone Piecewise Cubic
       Interpolation", SIAM J. Numer. Anal. 17(2), 1980.
.. [4] SCIPY-1.14.1-PCHIP. SciPy 1.14.1 documentation,
       ``scipy.interpolate.PchipInterpolator``: algorithm and implementation.
"""

from __future__ import annotations

import json
from bisect import bisect_left, bisect_right
from collections.abc import Mapping
from dataclasses import dataclass, field
from fractions import Fraction
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

from .exchange_documents import _surface_document
from .interpolation import _MULTILINEAR_POLICY, _multilinear_terms, _pchip_rational_slopes
from .numeric_outputs import _finite_rational, _fraction
from .serialization import _require_numerical_policies, canonical_sha256

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

    Parameters
    ----------
    name : str
        Non-blank condition name; unique within a surface.
    unit : str
        Non-blank unit label of the coordinate values.
    values : tuple of float
        At least two finite, strictly increasing grid coordinates.

    Raises
    ------
    DomainError
        If a label is blank or the coordinates are too few, non-finite or not
        strictly increasing.

    Notes
    -----
    The axis range is the declared validation domain of the surface along
    this condition [2]_; performance tables of [1]_ are parameterized by such
    operating conditions.

    References
    ----------
    .. [1] NASA, "Parametric Model of an Aerospike Rocket Engine".
           https://ntrs.nasa.gov/citations/20000031654
    .. [2] NASA-STD-7009B. NASA, "Standard for Models and Simulations", 2024.
           https://standards.nasa.gov/standard/NASA/NASA-STD-7009
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
        """Inclusive lower validation-domain coordinate, in axis units [1]_.

        References
        ----------
        .. [1] NASA-STD-7009B. NASA, "Standard for Models and Simulations",
               2024. https://standards.nasa.gov/standard/NASA/NASA-STD-7009
        """

        return self.values[0]

    @property
    def maximum(self) -> float:
        """Inclusive upper validation-domain coordinate, in axis units [1]_.

        References
        ----------
        .. [1] NASA-STD-7009B. NASA, "Standard for Models and Simulations",
               2024. https://standards.nasa.gov/standard/NASA/NASA-STD-7009
        """

        return self.values[-1]


@dataclass(frozen=True, slots=True)
class PerformanceScalePoint:
    """One measured or validated grid node of force and flow multipliers.

    Parameters
    ----------
    coordinates : tuple of float
        Non-empty node coordinates in the surface axis order.
    force_scale : float
        Positive finite thrust multiplier.
    mass_flow_scale : float
        Positive finite total tank-flow multiplier.

    Raises
    ------
    DomainError
        If the coordinates are empty or a multiplier is not positive and
        finite.

    Notes
    -----
    Force and flow multipliers are independent, so a measured operating map
    can carry a varying effective specific impulse ``F / (g0 * mdot)``, with
    ``g0`` the registered ``STANDARD_GRAVITY_M_S2``.  The
    nodes come from test data or a separately verified model such as the
    parametric tables of [1]_; [2]_ gives the transient-modelling context.

    References
    ----------
    .. [1] NASA, "Parametric Model of an Aerospike Rocket Engine".
           https://ntrs.nasa.gov/citations/20000031654
    .. [2] NASA, "Transient Mathematical Modeling for Liquid Rocket Engine
           Systems". https://ntrs.nasa.gov/citations/20040000363
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
    """Interpolated scales, analytic gradients and active grid cell.

    The multipliers come from a trajectory-facing performance table of the
    kind described in [1]_; the active cell and boundary axes record where in
    the declared domain the evaluation lies, for the model-use reporting of
    [2]_.

    Attributes
    ----------
    conditions : tuple of (str, float)
        Evaluated conditions in axis order, after endpoint snapping.
    force_scale, mass_flow_scale : float
        Interpolated thrust and total tank-flow multipliers.
    force_scale_gradient, mass_flow_scale_gradient : tuple of (str, float)
        Partial derivatives of each multiplier per axis, per axis unit.
    cell_lower_coordinates, cell_upper_coordinates : tuple of float
        Corners of the active grid cell.
    boundary_axes : tuple of str
        Axes on which the condition lies on a grid coordinate within the
        snapping tolerance.

    References
    ----------
    .. [1] NASA, "Parametric Model of an Aerospike Rocket Engine".
           https://ntrs.nasa.gov/citations/20000031654
    .. [2] NASA-STD-7009B. NASA, "Standard for Models and Simulations", 2024.
           https://standards.nasa.gov/standard/NASA/NASA-STD-7009
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

    Parameters
    ----------
    event_id : str
        Non-blank event identifier.
    condition_name : str
        Non-blank name of the operating condition.
    bound : {"lower", "upper"}
        Which end of the domain the event guards.
    boundary_value : float
        Finite domain limit, in axis units.
    unit : str
        Non-blank unit label.
    stops_burn : bool, default True
        Whether crossing the event stops the burn.

    Raises
    ------
    DomainError
        If a label is blank, ``bound`` is invalid or the limit is non-finite.

    Notes
    -----
    The event value is positive inside the declared validation domain [1]_,
    zero on its boundary and negative outside.  A propagator registers both
    events for each state-dependent condition as event detectors [2]_ and
    stops the burn before the surface would have to extrapolate.

    References
    ----------
    .. [1] NASA-STD-7009B. NASA, "Standard for Models and Simulations", 2024.
           https://standards.nasa.gov/standard/NASA/NASA-STD-7009
    .. [2] Orekit 13.1.5, "PropulsionModel / EventDetector" API
           documentation, ``org.orekit.forces.maneuvers.propulsion``
           (URL listed in docs/REFERENCES.md).
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
        """Return the positive-inside signed distance to the bound.

        Parameters
        ----------
        conditions : Mapping of str to float
            Operating conditions; must contain ``condition_name``.

        Returns
        -------
        float
            ``condition - boundary_value`` for a lower bound and
            ``boundary_value - condition`` for an upper bound, in axis units.

        Raises
        ------
        DomainError
            If ``conditions`` is not a mapping, lacks the condition, or holds
            a non-real value for it.
        """

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

    Parameters
    ----------
    event_id : str
        Non-blank event identifier.
    condition_name : str
        Non-blank name of the operating condition.
    coordinate : float
        Finite interior grid coordinate, in axis units.
    unit : str
        Non-blank unit label.
    stops_burn : bool, default False
        Must be ``False``; grid-line events never stop the burn.

    Raises
    ------
    DomainError
        If a label is blank, the coordinate is non-finite, or ``stops_burn``
        is ``True``.

    Notes
    -----
    Multilinear interpolation is continuous, but the derivative along an axis
    jumps where that condition crosses an interior grid coordinate.  A
    variable-step integrator that steps across such a kink without
    restarting degrades its error estimate, and gradient-based optimizers see
    a discontinuous Jacobian.  Registering these events as event detectors
    [1]_ lets the consumer land on each kink and restart, as it does at
    profile boundaries.  The event value is ``condition - coordinate``.

    Exactly on the coordinate, :meth:`RectilinearPerformanceSurface.evaluate`
    returns the gradient of the upper cell (the right-hand limit along that
    axis), except at the axis maximum, where only the lower cell exists.

    References
    ----------
    .. [1] Orekit 13.1.5, "PropulsionModel / EventDetector" API
           documentation, ``org.orekit.forces.maneuvers.propulsion``
           (URL listed in docs/REFERENCES.md).
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
        """Return the signed distance ``condition - coordinate``, in axis units.

        Parameters
        ----------
        conditions : Mapping of str to float
            Operating conditions; must contain ``condition_name``.

        Returns
        -------
        float
            Signed distance from the grid line.

        Raises
        ------
        DomainError
            If ``conditions`` is not a mapping, lacks the condition, or holds
            a non-real value for it.
        """

        if not isinstance(conditions, Mapping):
            raise DomainError("Operating conditions must be a mapping.")
        if self.condition_name not in conditions:
            raise DomainError(f"Missing operating condition {self.condition_name!r}.")
        value = real(conditions[self.condition_name], f"Operating condition {self.condition_name!r}")
        return value - self.coordinate


@dataclass(frozen=True, slots=True)
class RectilinearPerformanceSurface:
    """Complete, non-extrapolating N-dimensional force and tank-flow scale map.

    Parameters
    ----------
    axes : tuple of PerformanceAxis
        At least one axis, with unique names.
    points : tuple of PerformanceScalePoint
        Every node of the Cartesian product of the axis coordinates, each
        exactly once.
    source_id : str
        Non-blank identifier of the data or model source.
    source_sha256 : str
        SHA-256 hex digest of the source content.
    reference_ids : tuple of str, default PERFORMANCE_SURFACE_REFERENCE_IDS
        Non-empty, unique, non-blank source identifiers.
    schema : str, default PERFORMANCE_SURFACE_SCHEMA
        Document schema; no other value is accepted.
    interpolation : {"multilinear", "cubic_hermite"}, default "multilinear"
        Interpolation mode.

    Raises
    ------
    DomainError
        If the schema or mode is unsupported, a label or hash is invalid, the
        reference list is empty or repeats an entry, axis names repeat, or the
        nodes do not form the complete Cartesian grid.

    Notes
    -----
    The default mode is tensor-product multilinear interpolation on the
    rectilinear grid [2]_.  Its analytic gradient is constant within a cell
    only in one dimension; in higher dimensions each partial derivative
    varies multilinearly with the other coordinates.

    The ``cubic_hermite`` mode builds nodal values and mixed partial
    derivatives from the Fritsch-Carlson PCHIP slope operator [3]_, [4]_,
    applied axis by axis, and keeps them and the tensor sums in exact
    rational arithmetic until the checked value and gradient are exported.
    The interpolant is C1.  One-dimensional monotonicity does not establish
    positivity in the N-dimensional interior, so evaluation refuses a
    non-positive interpolated scale.

    The surface is a declared-domain interpolation mechanism, not a fitted
    physical law.  Its nodes come from applicable test data or a separately
    verified model, such as the parametric engine tables of [1]_, and its axis
    ranges are the declared validation domain [5]_.  The surface records a
    source hash and stable reference identifiers; [6]_ gives the
    transient-modelling context.

    References
    ----------
    .. [1] NASA, "Parametric Model of an Aerospike Rocket Engine".
           https://ntrs.nasa.gov/citations/20000031654
    .. [2] SCIPY-1.14.1-REGULAR-GRID. SciPy 1.14.1 documentation,
           ``scipy.interpolate.RegularGridInterpolator``: rectilinear grids
           and linear interpolation.
    .. [3] F. N. Fritsch and R. E. Carlson, "Monotone Piecewise Cubic
           Interpolation", SIAM J. Numer. Anal. 17(2), 1980.
    .. [4] SCIPY-1.14.1-PCHIP. SciPy 1.14.1 documentation,
           ``scipy.interpolate.PchipInterpolator``: algorithm and
           implementation.
    .. [5] NASA-STD-7009B. NASA, "Standard for Models and Simulations", 2024.
           https://standards.nasa.gov/standard/NASA/NASA-STD-7009
    .. [6] NASA, "Transient Mathematical Modeling for Liquid Rocket Engine
           Systems". https://ntrs.nasa.gov/citations/20040000363
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
    _hermite: Mapping[tuple[int, ...], Mapping[tuple[int, ...], tuple[Fraction, Fraction]]] = field(
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
    ) -> Mapping[tuple[int, ...], Mapping[tuple[int, ...], tuple[Fraction, Fraction]]]:
        """Return nodal values and mixed partials for tensor-product cubic Hermite.

        ``D[S]`` is the mixed partial over the axes in ``S`` (a sorted tuple);
        ``D[()]`` holds the nodal values.  It is built by applying the PCHIP
        slope operator [1]_, [2]_ along the largest axis of ``S`` to
        ``D[S minus that axis]``, one grid line at a time.  All entries are
        exact rationals ``(force, flow)`` keyed by node index.

        References
        ----------
        .. [1] F. N. Fritsch and R. E. Carlson, "Monotone Piecewise Cubic
               Interpolation", SIAM J. Numer. Anal. 17(2), 1980.
        .. [2] F. N. Fritsch and J. Butland, "A Method for Constructing Local
               Monotone Piecewise Cubic Interpolants", SIAM J. Sci. Stat.
               Comput. 5(2), 1984.
        """

        dimension = len(self.axes)
        shape = [len(axis.values) for axis in self.axes]
        base: dict[tuple[int, ...], tuple[Fraction, Fraction]] = {}
        for index in product(*(range(size) for size in shape)):
            coordinate = tuple(self.axes[a].values[i] for a, i in enumerate(index))
            point = self._point_index[coordinate]
            base[index] = (_fraction(point.force_scale), _fraction(point.mass_flow_scale))
        derivatives: dict[tuple[int, ...], Mapping[tuple[int, ...], tuple[Fraction, Fraction]]] = {
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
            result: dict[tuple[int, ...], tuple[Fraction, Fraction]] = {}
            others = [range(size) for a, size in enumerate(shape) if a != axis_index]
            for fixed in product(*others):
                line = []
                for i in range(shape[axis_index]):
                    index = list(fixed)
                    index.insert(axis_index, i)
                    line.append(tuple(index))
                force = _pchip_rational_slopes(values, tuple(source[item][0] for item in line))
                flow = _pchip_rational_slopes(values, tuple(source[item][1] for item in line))
                for position, item in enumerate(line):
                    result[item] = (force[position], flow[position])
            derivatives[subset] = MappingProxyType(result)
        return MappingProxyType(derivatives)

    @property
    def surface_hash(self) -> str:
        """SHA-256 content identity of nodes, units, domain and provenance.

        The digest of the canonical document without its hash field
        identifies content; it is not authentication.  [1]_ gives the
        model-credibility context.

        References
        ----------
        .. [1] NASA-STD-7009B. NASA, "Standard for Models and Simulations",
               2024. https://standards.nasa.gov/standard/NASA/NASA-STD-7009
        """

        return canonical_sha256(self.to_dict(include_hash=False))

    def to_dict(self, *, include_hash: bool = True) -> dict[str, object]:
        """Return the deterministic JSON-compatible surface document.

        Parameters
        ----------
        include_hash : bool, default True
            Whether to add ``surface_hash``.

        Returns
        -------
        dict
            Schema, source, axes, nodes sorted by coordinates, reference
            identifiers and the numerical-policy identifier; ``interpolation``
            is written only for a non-default mode.

        Notes
        -----
        The document is versioned, unit-explicit and records its source
        provenance.  [1]_ gives the model-credibility context and [2]_ the
        context for versioned, unit-explicit exchange documents.

        References
        ----------
        .. [1] NASA-STD-7009B. NASA, "Standard for Models and Simulations",
               2024. https://standards.nasa.gov/standard/NASA/NASA-STD-7009
        .. [2] CCSDS 502.0-B-3. CCSDS, "Orbit Data Messages", 2023.
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
                    "mass_flow_scale": point.mass_flow_scale,
                }
                for point in sorted(self.points, key=lambda item: item.coordinates)
            ],
            "reference_ids": list(self.reference_ids),
        }
        if self.interpolation != "multilinear":
            payload["interpolation"] = self.interpolation
            payload["interpolation_numerics"] = "rational_pchip_tensor_hermite_export_v1"
        else:
            payload["interpolation_numerics"] = _MULTILINEAR_POLICY
        if include_hash:
            payload["surface_hash"] = canonical_sha256(payload)
        return payload

    def to_json(self, *, indent: int | None = 2) -> str:
        """Serialize the complete hash-bound surface as JSON.

        Parameters
        ----------
        indent : int or None, default 2
            Indentation passed to `json.dumps`.

        Returns
        -------
        str
            Sorted JSON without NaN or Infinity tokens [3]_.

        Notes
        -----
        See `to_dict` for the document content and its model-credibility [1]_
        and exchange-document [2]_ context.

        References
        ----------
        .. [1] NASA-STD-7009B. NASA, "Standard for Models and Simulations",
               2024. https://standards.nasa.gov/standard/NASA/NASA-STD-7009
        .. [2] CCSDS 502.0-B-3. CCSDS, "Orbit Data Messages", 2023.
               https://ccsds.org/Pubs/502x0b3e1.pdf
        .. [3] IETF RFC 8259, "The JavaScript Object Notation (JSON) Data
               Interchange Format". https://www.rfc-editor.org/rfc/rfc8259
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
        """Reconstruct a surface and reject malformed or hash-altered content.

        Parameters
        ----------
        payload : Mapping
            Document produced by `to_dict`.

        Returns
        -------
        RectilinearPerformanceSurface
            The reconstructed surface.

        Raises
        ------
        InputError
            If the document is malformed, a field violates the surface
            contract, the hash does not match, or a numerical-policy
            identifier is missing or differs.

        Notes
        -----
        Every read verifies the content hash, so a reopened surface is the
        recorded one.  [1]_ gives the model-credibility context and [2]_ the
        context for versioned, unit-explicit exchange documents.

        References
        ----------
        .. [1] NASA-STD-7009B. NASA, "Standard for Models and Simulations",
               2024. https://standards.nasa.gov/standard/NASA/NASA-STD-7009
        .. [2] CCSDS 502.0-B-3. CCSDS, "Orbit Data Messages", 2023.
               https://ccsds.org/Pubs/502x0b3e1.pdf
        """

        _surface_document(payload)
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
        _require_numerical_policies(payload, surface.to_dict(include_hash=False))
        return surface

    @classmethod
    def from_json(cls, text: str) -> RectilinearPerformanceSurface:
        """Parse a strict JSON surface document [3]_ and check it with `from_dict`.

        Parameters
        ----------
        text : str
            JSON text produced by `to_json`.

        Returns
        -------
        RectilinearPerformanceSurface
            The reconstructed surface.

        Raises
        ------
        InputError
            If the text is not a valid JSON object or `from_dict` rejects it.

        Notes
        -----
        See `from_dict` for the integrity checks and their model-credibility
        [1]_ and exchange-document [2]_ context.

        References
        ----------
        .. [1] NASA-STD-7009B. NASA, "Standard for Models and Simulations",
               2024. https://standards.nasa.gov/standard/NASA/NASA-STD-7009
        .. [2] CCSDS 502.0-B-3. CCSDS, "Orbit Data Messages", 2023.
               https://ccsds.org/Pubs/502x0b3e1.pdf
        .. [3] IETF RFC 8259, "The JavaScript Object Notation (JSON) Data
               Interchange Format". https://www.rfc-editor.org/rfc/rfc8259
        """

        return read_json_document(text, "Performance-surface artifact", cls.from_dict)

    def validity_surfaces(self) -> tuple[PerformanceValiditySurface, ...]:
        """Return stopping lower and upper validity events for every axis.

        Returns
        -------
        tuple of PerformanceValiditySurface
            For each axis in order, the ``"lower"`` event at its minimum and
            the ``"upper"`` event at its maximum, with identifiers
            ``performance:<axis>:minimum`` and ``performance:<axis>:maximum``.

        Notes
        -----
        The events bound the declared validation domain [1]_ and are meant to
        be registered as propagator event detectors [2]_.

        References
        ----------
        .. [1] NASA-STD-7009B. NASA, "Standard for Models and Simulations",
               2024. https://standards.nasa.gov/standard/NASA/NASA-STD-7009
        .. [2] Orekit 13.1.5, "PropulsionModel / EventDetector" API
               documentation, ``org.orekit.forces.maneuvers.propulsion``
               (URL listed in docs/REFERENCES.md).
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

        Returns
        -------
        tuple of PerformanceGridLineSurface
            Events with identifiers ``performance:<axis>:grid:<index>``; empty
            for a ``cubic_hermite`` surface.

        Notes
        -----
        Domain endpoints are excluded; the stopping :meth:`validity_surfaces`
        cover them.  A ``cubic_hermite`` surface is C1, so it has no
        derivative kinks.  The events are meant to be registered as
        propagator event detectors [1]_.

        References
        ----------
        .. [1] Orekit 13.1.5, "PropulsionModel / EventDetector" API
               documentation, ``org.orekit.forces.maneuvers.propulsion``
               (URL listed in docs/REFERENCES.md).
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

        Parameters
        ----------
        conditions : Mapping of str to float
            One value per axis, keyed exactly by the axis names.

        Returns
        -------
        PerformanceSurfaceEvaluation
            Interpolated multipliers, their gradients and the active cell.

        Raises
        ------
        DomainError
            If ``conditions`` is not a mapping, its keys differ from the axis
            names, a value is not real, a value lies outside an axis range
            beyond the snapping tolerance, or the ``cubic_hermite`` interpolant
            is non-positive at the query.

        Notes
        -----
        A value within ``max(1e-12, 1e-12 * max(|min|, |max|))`` of a domain
        endpoint is snapped to that endpoint; values farther outside are
        refused, keeping evaluation inside the declared validation domain
        [2]_.  Multilinear weights [3]_ and cubic Hermite bases built on PCHIP
        derivatives [4]_ are evaluated in exact rational arithmetic and
        rounded once on export.  [1]_ gives an example of trajectory-facing
        engine performance tables parameterized by operating conditions.

        References
        ----------
        .. [1] NASA, "Parametric Model of an Aerospike Rocket Engine".
               https://ntrs.nasa.gov/citations/20000031654
        .. [2] NASA-STD-7009B. NASA, "Standard for Models and Simulations",
               2024. https://standards.nasa.gov/standard/NASA/NASA-STD-7009
        .. [3] SCIPY-1.14.1-REGULAR-GRID. SciPy 1.14.1 documentation,
               ``scipy.interpolate.RegularGridInterpolator``: rectilinear grids
               and linear interpolation.
        .. [4] F. N. Fritsch and R. E. Carlson, "Monotone Piecewise Cubic
               Interpolation", SIAM J. Numer. Anal. 17(2), 1980.
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
        brackets: list[tuple[int, int, float | Fraction, float | Fraction]] = []
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
            width = _fraction(axis.values[upper]) - _fraction(axis.values[lower])
            fraction = (_fraction(value) - _fraction(axis.values[lower])) / width
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
        force_scale = Fraction(0)
        flow_scale = Fraction(0)
        force_gradient = [Fraction(0)] * len(self.axes)
        flow_gradient = [Fraction(0)] * len(self.axes)

        for corner, weight, derivative_weights in _multilinear_terms(brackets):
            coordinate = tuple(
                axis.values[brackets[index][corner[index]]]
                for index, axis in enumerate(self.axes)
            )
            point = point_by_coordinate[coordinate]
            force, flow = _fraction(point.force_scale), _fraction(point.mass_flow_scale)
            force_scale += weight * force
            flow_scale += weight * flow
            for derivative_axis, derivative_weight in enumerate(derivative_weights):
                force_gradient[derivative_axis] += derivative_weight * force
                flow_gradient[derivative_axis] += derivative_weight * flow

        return self._evaluation(
            normalized, brackets, boundary_axes,
            _finite_rational(force_scale, "Multilinear force scale"),
            _finite_rational(flow_scale, "Multilinear mass-flow scale"),
            [_finite_rational(value, f"Multilinear force gradient: {axis.name}")
             for axis, value in zip(self.axes, force_gradient, strict=True)],
            [_finite_rational(value, f"Multilinear mass-flow gradient: {axis.name}")
             for axis, value in zip(self.axes, flow_gradient, strict=True)],
        )

    def _evaluate_hermite(
        self, brackets: list[tuple[int, int, float | Fraction, float | Fraction]]
    ) -> tuple[float, float, list[float], list[float]]:
        """Evaluate the tensor-product cubic Hermite interpolant and its gradient.

        On a cell with local coordinate ``t`` in ``[0, 1]`` and width ``h``,
        the one-dimensional cubic Hermite bases are ``2t^3 - 3t^2 + 1`` and
        ``-2t^3 + 3t^2`` for the corner values and ``h (t^3 - 2t^2 + t)`` and
        ``h (t^3 - t^2)`` for the corner slopes.  Each stored mixed partial
        ``D[S]`` is weighted by the slope bases on the axes in ``S`` and the
        value bases elsewhere, using PCHIP derivatives [1]_.

        Returns
        -------
        tuple
            ``(force_scale, flow_scale, force_gradient, flow_gradient)``,
            each rounded once from its exact value.

        References
        ----------
        .. [1] F. N. Fritsch and R. E. Carlson, "Monotone Piecewise Cubic
               Interpolation", SIAM J. Numer. Anal. 17(2), 1980.
        """
        dimension = len(self.axes)
        force = Fraction(0)
        flow = Fraction(0)
        force_gradient = [Fraction(0)] * dimension
        flow_gradient = [Fraction(0)] * dimension
        basis_cache = []
        for _, _, t, width in brackets:
            t, width = _fraction(t), _fraction(width)
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
        return (
            _finite_rational(force, "Cubic force scale"),
            _finite_rational(flow, "Cubic mass-flow scale"),
            [_finite_rational(value, f"Cubic force gradient: {axis.name}")
             for axis, value in zip(self.axes, force_gradient, strict=True)],
            [_finite_rational(value, f"Cubic mass-flow gradient: {axis.name}")
             for axis, value in zip(self.axes, flow_gradient, strict=True)],
        )

    def _evaluation(
        self,
        normalized: dict[str, float],
        brackets: list[tuple[int, int, float | Fraction, float | Fraction]],
        boundary_axes: list[str],
        force_scale: float,
        flow_scale: float,
        force_gradient: list[float],
        flow_gradient: list[float],
    ) -> PerformanceSurfaceEvaluation:
        """Assemble the evaluation record in axis order."""
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

