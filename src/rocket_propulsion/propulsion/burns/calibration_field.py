"""Hot-fire calibration fields over a rectilinear operating-condition surface.

A :class:`HotFireCalibrationField` places one random-effects hot-fire
calibration [1]_ [2]_ at every node of a complete rectilinear grid of
operating conditions (for example chamber or supply pressure by mixture
ratio) and interpolates, inside that grid only, the consensus correction
scale, its standard uncertainty and the within- and between-engine process
variances. Interpolation weights are the tensor-product multilinear weights
of the performance surfaces [3]_; the domain check, snapping tolerance and
validity events follow
:class:`~rocket_propulsion.propulsion.burns.performance_surface.RectilinearPerformanceSurface`
so that a propagator stops before either layer would extrapolate [4]_.

Notes
-----
Uncertainty interpolation rule. For non-negative multilinear weights ``w_c``
(summing to one) on the corner nodes ``c`` of the active cell:

* consensus scale ``mu(x) = sum_c w_c mu_c``;
* consensus standard uncertainty ``u(x) = sum_c w_c u_c``. ``mu(x)`` is a
  linear function of the node estimates with sensitivity coefficients
  ``w_c``, so its combined variance follows the law of propagation of
  uncertainty with correlation coefficients ``r`` [5]_. For ``r = +1`` that
  law reduces to the linear sum ``sum_c w_c u_c`` (JCGM 100:2008, 5.2.2
  NOTE 1); since ``|r| <= 1`` and ``w_c >= 0`` it is also an upper bound for
  every correlation between node estimates. The independent-node value
  ``sqrt(sum_c w_c**2 u_c**2)`` (JCGM 100:2008, equation 10) is reported
  alongside as the lower, uncorrelated case;
* process variances ``sigma_w**2(x) = sum_c w_c sigma_w,c**2`` and
  ``sigma_b**2(x) = sum_c w_c sigma_b,c**2``, i.e. multilinear interpolation
  of the variances themselves. No propagation law applies to a variance
  *parameter*; this is this library's stated engineering convention. It
  reproduces every node, stays non-negative, keeps the population identity
  ``sigma_w**2 + sigma_b**2`` [6]_ exact at every point, and by the
  Cauchy-Schwarz inequality ``sum_c w_c s_c**2 >= (sum_c w_c s_c)**2`` is
  never smaller than linear interpolation of the standard deviations.

The field interpolates calibrated evidence; it does not estimate the
interpolation (model-form) error between nodes, the sampling uncertainty of
the variance components, or covariance across nodes beyond the envelope
above [7]_.

References
----------
.. [1] NIST/SEMATECH e-Handbook of Statistical Methods, section 7.4.4,
       "Variance Components" (one-way random-effects model).
       https://www.itl.nist.gov/div898/handbook/prc/section4/prc44.htm
.. [2] NIST Dataplot, "Consensus Mean" (heteroscedastic random effects and
       Mandel-Paule estimation).
       https://www.itl.nist.gov/div898/software/dataplot/refman1/auxillar/consmean.htm
.. [3] SCIPY-1.14.1-REGULAR-GRID. SciPy 1.14.1 documentation,
       ``scipy.interpolate.RegularGridInterpolator``: rectilinear grids and
       linear interpolation.
.. [4] NASA-STD-7009B, "Standard for Models and Simulations".
       https://standards.nasa.gov/standard/NASA/NASA-STD-7009
.. [5] JCGM-100-2008-SEC5. JCGM 100:2008, sections 5.1.2 (equation 10) and
       5.2.2 (equations 13-16, NOTE 1).
       https://www.bipm.org/documents/20126/2071204/JCGM_100_2008_E.pdf
.. [6] MIT-TOTAL-VARIANCE. MIT RES.6-012, lecture 13.7, conditional and
       between-group variance identity (law of total variance).
.. [7] NASA-HDBK-7009, "NASA Handbook for Models and Simulations".
       https://ntrs.nasa.gov/citations/20140002378
"""

from __future__ import annotations

import json
from bisect import bisect_left, bisect_right
from collections.abc import Mapping
from dataclasses import dataclass, field
from fractions import Fraction
from math import isfinite, sqrt
from types import MappingProxyType
from typing import Any

from rocket_propulsion.core.errors import DomainError, InputError
from rocket_propulsion.core.validation import (
    DOCUMENT_ERRORS,
    real,
    real_tuple,
    records,
    strings,
    text,
)

from .calibration import (
    CalibrationScope,
    HotFireScaleCalibration,
    HotFireScaleMeasurement,
    calibrate_hot_fire_scale,
)
from .discrepancy import qualified_scale_component_name
from .exchange_documents import _read_document
from .interpolation import _MULTILINEAR_POLICY, _multilinear_terms
from .numeric_outputs import _finite_rational, _fraction
from .performance_surface import (
    PerformanceAxis,
    PerformanceGridLineSurface,
    PerformanceScalePoint,
    PerformanceSurfaceEvaluation,
    PerformanceValiditySurface,
    RectilinearPerformanceSurface,
)
from .serialization import _require_numerical_policies, canonical_json_bytes, canonical_sha256
from .uncertainty import DistributionKind, UncertainParameter, UncertaintyClass

CALIBRATION_FIELD_SCHEMA = "rocket_propulsion_hot_fire_calibration_field_v1"
CALIBRATION_FIELD_UNCERTAINTY_POLICY = (
    "consensus_u_linear_envelope_variance_components_multilinear_v1"
)
CALIBRATION_FIELD_REFERENCE_IDS = (
    "NIST-SEMATECH-7.4.4",
    "NIST-DATAPLOT-CONSENSUS-MEAN",
    "JCGM-100-2008",
    "NASA-STD-7009B",
    "NASA-HDBK-7009",
)
_SUPPORTED_SCALE_NAMES = (
    "thrust_scale",
    "specific_impulse_scale",
    "time_scale",
    "axial_efficiency_scale",
)
_FIELD_WARNINGS = (
    ("Between grid nodes the calibration is interpolated, not measured; interpolation "
    "(model-form) error is not included."),
    ("Consensus standard uncertainty is the linear envelope sum_c w_c u_c, an upper bound "
    "for any correlation between node estimates."),
    "Process variances are interpolated linearly in variance (library convention).",
    "Variance-component sampling uncertainty and model-form discrepancy are not included.",
)


@dataclass(frozen=True, slots=True)
class CalibrationFieldNode:
    """One hot-fire calibration located at one grid node.

    Parameters
    ----------
    coordinates : tuple of float
        Node coordinates in the field axis order.
    calibration : HotFireScaleCalibration
        Random-effects calibration of the measurements taken at this
        operating condition [1]_.

    Raises
    ------
    DomainError
        If the coordinates are empty or non-finite, or ``calibration`` is not
        a :class:`HotFireScaleCalibration`.

    References
    ----------
    .. [1] NIST/SEMATECH e-Handbook of Statistical Methods, section 7.4.4,
           "Variance Components" (one-way random-effects model).
           https://www.itl.nist.gov/div898/handbook/prc/section4/prc44.htm
    """

    coordinates: tuple[float, ...]
    calibration: HotFireScaleCalibration

    def __post_init__(self) -> None:
        object.__setattr__(
            self, "coordinates", real_tuple(self.coordinates, "Calibration-field node coordinates")
        )
        if not self.coordinates:
            raise DomainError("Calibration-field node coordinates must not be empty.")
        if not isinstance(self.calibration, HotFireScaleCalibration):
            raise DomainError("Calibration-field node requires a HotFireScaleCalibration.")


@dataclass(frozen=True, slots=True)
class CalibrationFieldEvaluation:
    """Interpolated calibration at one in-domain operating point.

    Attributes
    ----------
    parameter_name : str
        Calibrated scale name.
    conditions : tuple of (str, float)
        Evaluated conditions in axis order, after endpoint snapping.
    consensus_scale : float
        ``sum_c w_c mu_c`` with multilinear weights ``w_c`` [1]_.
    consensus_standard_uncertainty : float
        ``sum_c w_c u_c``: the combined standard uncertainty for fully
        correlated node estimates and an upper bound for any correlation
        [2]_.
    independent_nodes_standard_uncertainty : float
        ``sqrt(sum_c w_c**2 u_c**2)``, the uncorrelated-node case [2]_.
    within_engine_process_variance, between_engine_process_variance : float
        Multilinearly interpolated variance components (library convention).
    population_process_variance : float
        Their sum [3]_.
    consensus_scale_gradient : tuple of (str, float)
        Partial derivative of ``consensus_scale`` per axis, per axis unit.
    cell_lower_coordinates, cell_upper_coordinates : tuple of float
        Corners of the active grid cell.
    boundary_axes : tuple of str
        Axes on which the condition lies on a grid coordinate.
    node_weights : tuple of (tuple of float, float)
        Corner nodes with a nonzero interpolation weight and that weight.

    References
    ----------
    .. [1] SCIPY-1.14.1-REGULAR-GRID. SciPy 1.14.1 documentation,
           ``scipy.interpolate.RegularGridInterpolator``: rectilinear grids
           and linear interpolation.
    .. [2] JCGM-100-2008-SEC5. JCGM 100:2008, sections 5.1.2 (equation 10)
           and 5.2.2 (equations 13-16, NOTE 1).
           https://www.bipm.org/documents/20126/2071204/JCGM_100_2008_E.pdf
    .. [3] MIT-TOTAL-VARIANCE. MIT RES.6-012, lecture 13.7, conditional and
           between-group variance identity (law of total variance).
    """

    parameter_name: str
    conditions: tuple[tuple[str, float], ...]
    consensus_scale: float
    consensus_standard_uncertainty: float
    independent_nodes_standard_uncertainty: float
    within_engine_process_variance: float
    between_engine_process_variance: float
    population_process_variance: float
    consensus_scale_gradient: tuple[tuple[str, float], ...]
    cell_lower_coordinates: tuple[float, ...]
    cell_upper_coordinates: tuple[float, ...]
    boundary_axes: tuple[str, ...]
    node_weights: tuple[tuple[tuple[float, ...], float], ...]

    @property
    def relative_standard_uncertainty(self) -> float:
        """``consensus_standard_uncertainty / consensus_scale`` [1]_.

        References
        ----------
        .. [1] JCGM-100-2008-SEC5. JCGM 100:2008, section 5.1.6 (equation 12),
               relative standard uncertainty.
               https://www.bipm.org/documents/20126/2071204/JCGM_100_2008_E.pdf
        """

        return self.consensus_standard_uncertainty / self.consensus_scale


def _cell_brackets(
    axes: tuple[PerformanceAxis, ...], conditions: Mapping[str, float]
) -> tuple[dict[str, float], list[tuple[int, int, Fraction, Fraction]], list[str]]:
    """Locate the active grid cell without extrapolation.

    Applies the same rules as ``RectilinearPerformanceSurface.evaluate``:
    exact key match, endpoint snapping within
    ``max(1e-12, 1e-12 * max(|min|, |max|))`` and refusal beyond it, keeping
    evaluation inside the declared validation domain [1]_.

    References
    ----------
    .. [1] NASA-STD-7009B, "Standard for Models and Simulations".
           https://standards.nasa.gov/standard/NASA/NASA-STD-7009
    """
    if not isinstance(conditions, Mapping):
        raise DomainError("Operating conditions must be a mapping.")
    expected = {axis.name for axis in axes}
    supplied = set(conditions)
    if supplied != expected:
        raise DomainError(
            "Operating conditions must match calibration-field axes exactly "
            f"(missing={sorted(expected - supplied)}, extra={sorted(supplied - expected)})."
        )
    normalized: dict[str, float] = {}
    brackets: list[tuple[int, int, Fraction, Fraction]] = []
    boundary_axes: list[str] = []
    for axis in axes:
        value = real(conditions[axis.name], f"Operating condition {axis.name!r}")
        tolerance = max(1e-12, max(abs(axis.minimum), abs(axis.maximum)) * 1e-12)
        if value < axis.minimum - tolerance or value > axis.maximum + tolerance:
            raise DomainError(
                f"Operating condition {axis.name!r} lies outside the calibrated domain "
                f"[{axis.minimum}, {axis.maximum}] {axis.unit}; extrapolation refused."
            )
        value = min(axis.maximum, max(axis.minimum, value))
        normalized[axis.name] = value
        nearest = bisect_left(axis.values, value - tolerance)
        if nearest < len(axis.values) and axis.values[nearest] <= value + tolerance:
            boundary_axes.append(axis.name)
        upper = min(bisect_right(axis.values, value), len(axis.values) - 1)
        lower = max(0, upper - 1)
        width = _fraction(axis.values[upper]) - _fraction(axis.values[lower])
        brackets.append(
            (lower, upper, (_fraction(value) - _fraction(axis.values[lower])) / width, width)
        )
    return normalized, brackets, boundary_axes


@dataclass(frozen=True, slots=True)
class HotFireCalibrationField:
    """Complete, non-extrapolating field of one hot-fire calibrated scale.

    Parameters
    ----------
    parameter_name : str
        ``thrust_scale``, ``specific_impulse_scale``, ``time_scale`` or
        ``axial_efficiency_scale``; every node calibration must calibrate it.
    axes : tuple of PerformanceAxis
        At least one strictly increasing operating-condition axis with a
        unique name.
    nodes : tuple of CalibrationFieldNode
        Exactly one calibration at every node of the Cartesian product of the
        axis coordinates. Missing nodes, duplicate nodes, coordinates outside
        the grid or between grid lines, and a ``condition_id`` used at two
        nodes are refused.
    source, rationale : str
        Non-blank provenance and justification.
    reference_ids : tuple of str, default CALIBRATION_FIELD_REFERENCE_IDS
        Non-empty, unique source identifiers.
    schema : str, default CALIBRATION_FIELD_SCHEMA
        Document schema; no other value is accepted.

    Raises
    ------
    DomainError
        If any constraint above is violated.

    Notes
    -----
    Each node is an independent stratum calibrated by the heteroscedastic
    Mandel-Paule random-effects model [1]_ [2]_. Between nodes the consensus
    scale and process variances are interpolated with multilinear weights
    [3]_ and the consensus standard uncertainty with the linear envelope of
    JCGM 100:2008, 5.2.2 NOTE 1 [4]_; see the module docstring for the full
    rule and its justification. The axis ranges are the declared calibration
    domain [5]_ and :meth:`evaluate` refuses extrapolation.

    References
    ----------
    .. [1] NIST/SEMATECH e-Handbook of Statistical Methods, section 7.4.4,
           "Variance Components" (one-way random-effects model).
           https://www.itl.nist.gov/div898/handbook/prc/section4/prc44.htm
    .. [2] NIST Dataplot, "Consensus Mean" (heteroscedastic random effects and
           Mandel-Paule estimation).
           https://www.itl.nist.gov/div898/software/dataplot/refman1/auxillar/consmean.htm
    .. [3] SCIPY-1.14.1-REGULAR-GRID. SciPy 1.14.1 documentation,
           ``scipy.interpolate.RegularGridInterpolator``: rectilinear grids
           and linear interpolation.
    .. [4] JCGM-100-2008-SEC5. JCGM 100:2008, sections 5.1.2 (equation 10)
           and 5.2.2 (equations 13-16, NOTE 1).
           https://www.bipm.org/documents/20126/2071204/JCGM_100_2008_E.pdf
    .. [5] NASA-STD-7009B, "Standard for Models and Simulations".
           https://standards.nasa.gov/standard/NASA/NASA-STD-7009
    """

    parameter_name: str
    axes: tuple[PerformanceAxis, ...]
    nodes: tuple[CalibrationFieldNode, ...]
    source: str
    rationale: str
    reference_ids: tuple[str, ...] = CALIBRATION_FIELD_REFERENCE_IDS
    schema: str = CALIBRATION_FIELD_SCHEMA
    _node_index: Mapping[tuple[float, ...], HotFireScaleCalibration] = field(
        init=False, repr=False, compare=False
    )

    def __post_init__(self) -> None:
        if self.schema != CALIBRATION_FIELD_SCHEMA:
            raise DomainError("Unsupported calibration-field schema.")
        text(self.parameter_name, "Calibration-field parameter name")
        if self.parameter_name not in _SUPPORTED_SCALE_NAMES:
            raise DomainError("Unsupported calibration-field scale parameter.")
        text(self.source, "Calibration-field source")
        text(self.rationale, "Calibration-field rationale")
        object.__setattr__(self, "axes", records(self.axes, PerformanceAxis, "Calibration-field axes"))
        object.__setattr__(
            self, "nodes", records(self.nodes, CalibrationFieldNode, "Calibration-field nodes")
        )
        object.__setattr__(
            self,
            "reference_ids",
            strings(self.reference_ids, "Calibration-field references", allow_blank=True),
        )
        if not self.reference_ids or any(not item.strip() for item in self.reference_ids):
            raise DomainError("Calibration-field references must not be empty.")
        if len(self.reference_ids) != len(set(self.reference_ids)):
            raise DomainError("Calibration-field references must be unique.")
        if not self.axes:
            raise DomainError("A calibration field requires at least one axis.")
        names = [axis.name for axis in self.axes]
        if len(names) != len(set(names)):
            raise DomainError("Calibration-field axis names must be unique.")

        index: dict[tuple[float, ...], HotFireScaleCalibration] = {}
        condition_ids: set[str] = set()
        for node in self.nodes:
            if len(node.coordinates) != len(self.axes):
                raise DomainError("Calibration-field node dimensionality must match the axes.")
            for axis, value in zip(self.axes, node.coordinates, strict=True):
                if value < axis.minimum or value > axis.maximum:
                    raise DomainError(
                        f"Calibration node lies outside the grid on axis {axis.name!r} "
                        f"({value} not in [{axis.minimum}, {axis.maximum}] {axis.unit})."
                    )
                if value not in axis.values:
                    raise DomainError(
                        f"Calibration node is not on a grid coordinate of axis {axis.name!r} "
                        f"({value} {axis.unit})."
                    )
            if node.coordinates in index:
                raise DomainError(f"Duplicate calibration-field node at {node.coordinates}.")
            calibration = node.calibration
            if calibration.parameter_name != self.parameter_name:
                raise DomainError(
                    "Every node calibration must calibrate the field parameter "
                    f"{self.parameter_name!r}."
                )
            if calibration.condition_id in condition_ids:
                raise DomainError(
                    f"Operating condition {calibration.condition_id!r} is assigned to two nodes."
                )
            condition_ids.add(calibration.condition_id)
            index[node.coordinates] = calibration
        expected = len(self._grid_coordinates())
        if len(index) != expected:
            raise DomainError(
                "Calibration field must contain the complete Cartesian grid "
                f"(missing={expected - len(index)} of {expected} nodes)."
            )
        object.__setattr__(self, "_node_index", MappingProxyType(index))

    def _grid_coordinates(self) -> list[tuple[float, ...]]:
        """Return every grid node in lexicographic axis order."""
        coordinates: list[tuple[float, ...]] = [()]
        for axis in self.axes:
            coordinates = [prefix + (value,) for prefix in coordinates for value in axis.values]
        return coordinates

    def node_calibration(self, coordinates: tuple[float, ...]) -> HotFireScaleCalibration:
        """Return the calibration stored at exactly ``coordinates``.

        Parameters
        ----------
        coordinates : tuple of float
            Grid node in axis order.

        Returns
        -------
        HotFireScaleCalibration
            Node calibration [1]_.

        Raises
        ------
        DomainError
            If ``coordinates`` is not a grid node.

        References
        ----------
        .. [1] NIST/SEMATECH e-Handbook of Statistical Methods, section 7.4.4,
               "Variance Components" (one-way random-effects model).
               https://www.itl.nist.gov/div898/handbook/prc/section4/prc44.htm
        """

        key = real_tuple(coordinates, "Calibration-field node coordinates")
        if key not in self._node_index:
            raise DomainError(f"{key} is not a calibration-field grid node.")
        return self._node_index[key]

    @property
    def warnings(self) -> tuple[str, ...]:
        """Field interpretation warnings followed by distinct node warnings.

        The field warnings record the interpolation assumptions; node
        warnings (for example fewer than six engines [1]_) are carried over
        once each.

        References
        ----------
        .. [1] NIST Dataplot, "Consensus Mean" (heteroscedastic random effects
               and Mandel-Paule estimation).
               https://www.itl.nist.gov/div898/software/dataplot/refman1/auxillar/consmean.htm
        """

        node_warnings = (
            warning for node in self.nodes for warning in node.calibration.warnings
        )
        return tuple(dict.fromkeys((*_FIELD_WARNINGS, *node_warnings)))

    @property
    def field_hash(self) -> str:
        """SHA-256 content identity of axes, node calibrations and provenance.

        The digest of the canonical document without its hash field is a
        reproducible content identity, not an authentication [1]_.

        References
        ----------
        .. [1] NASA-STD-7009B, "Standard for Models and Simulations".
               https://standards.nasa.gov/standard/NASA/NASA-STD-7009
        """

        return canonical_sha256(self.to_dict(include_hash=False))

    def validity_surfaces(self) -> tuple[PerformanceValiditySurface, ...]:
        """Return stopping lower and upper calibrated-domain events per axis.

        Returns
        -------
        tuple of PerformanceValiditySurface
            For each axis, events ``calibration-field:<parameter>:<axis>:minimum``
            and ``...:maximum``; positive inside, zero on and negative outside
            the calibrated domain [1]_, for registration as propagator event
            detectors [2]_.

        References
        ----------
        .. [1] NASA-STD-7009B, "Standard for Models and Simulations".
               https://standards.nasa.gov/standard/NASA/NASA-STD-7009
        .. [2] Orekit 13.1.5, "PropulsionModel / EventDetector" API
               documentation, ``org.orekit.forces.maneuvers.propulsion``
               (URL listed in docs/REFERENCES.md).
        """

        prefix = f"calibration-field:{self.parameter_name}"
        return tuple(
            surface
            for axis in self.axes
            for surface in (
                PerformanceValiditySurface(
                    f"{prefix}:{axis.name}:minimum", axis.name, "lower", axis.minimum, axis.unit
                ),
                PerformanceValiditySurface(
                    f"{prefix}:{axis.name}:maximum", axis.name, "upper", axis.maximum, axis.unit
                ),
            )
        )

    def grid_line_surfaces(self) -> tuple[PerformanceGridLineSurface, ...]:
        """Return one non-stopping kink event per interior axis coordinate.

        Returns
        -------
        tuple of PerformanceGridLineSurface
            Events ``calibration-field:<parameter>:<axis>:grid:<index>``; the
            multilinear field is continuous but its gradient jumps there [1]_.

        References
        ----------
        .. [1] Orekit 13.1.5, "PropulsionModel / EventDetector" API
               documentation, ``org.orekit.forces.maneuvers.propulsion``
               (URL listed in docs/REFERENCES.md).
        """

        return tuple(
            PerformanceGridLineSurface(
                f"calibration-field:{self.parameter_name}:{axis.name}:grid:{index}",
                axis.name,
                coordinate,
                axis.unit,
            )
            for axis in self.axes
            for index, coordinate in enumerate(axis.values[1:-1], start=1)
        )

    def evaluate(self, conditions: Mapping[str, float]) -> CalibrationFieldEvaluation:
        """Interpolate the calibration at an in-domain operating point.

        Parameters
        ----------
        conditions : Mapping of str to float
            One value per axis, keyed exactly by the axis names.

        Returns
        -------
        CalibrationFieldEvaluation
            Interpolated consensus, uncertainties, variance components,
            gradient and active cell.

        Raises
        ------
        DomainError
            If the keys differ from the axis names, a value is not real, or a
            value lies outside an axis range beyond the snapping tolerance
            (extrapolation refused).

        Notes
        -----
        Weights ``w_c`` are the exact rational multilinear weights [1]_; all
        sums are formed exactly and rounded once on export. At a node the
        weight of that node is one, so every output equals the node
        calibration exactly. The uncertainty rule (linear envelope
        ``sum_c w_c u_c`` from JCGM 100:2008 5.2.2 NOTE 1 [2]_ and
        variance-linear process variances) is stated in the module docstring.

        References
        ----------
        .. [1] SCIPY-1.14.1-REGULAR-GRID. SciPy 1.14.1 documentation,
               ``scipy.interpolate.RegularGridInterpolator``: rectilinear
               grids and linear interpolation.
        .. [2] JCGM-100-2008-SEC5. JCGM 100:2008, sections 5.1.2 (equation
               10) and 5.2.2 (equations 13-16, NOTE 1).
               https://www.bipm.org/documents/20126/2071204/JCGM_100_2008_E.pdf
        """

        normalized, brackets, boundary_axes = _cell_brackets(self.axes, conditions)
        consensus = Fraction(0)
        envelope = Fraction(0)
        independent_variance = Fraction(0)
        within = Fraction(0)
        between = Fraction(0)
        gradient = [Fraction(0)] * len(self.axes)
        weights: list[tuple[tuple[float, ...], float]] = []
        for corner, weight, derivative_weights in _multilinear_terms(brackets):
            coordinate = tuple(
                # ``[:2]`` keeps the integer (lower, upper) grid indices of each bracket.
                axis.values[brackets[index][:2][corner[index]]]
                for index, axis in enumerate(self.axes)
            )
            calibration = self._node_index[coordinate]
            mu = _fraction(calibration.consensus_scale)
            u = _fraction(calibration.consensus_standard_uncertainty)
            consensus += weight * mu
            envelope += weight * u
            independent_variance += weight * weight * u * u
            within += weight * _fraction(calibration.within_engine_process_variance)
            between += weight * _fraction(calibration.between_engine_process_variance)
            for axis_index, derivative in enumerate(derivative_weights):
                gradient[axis_index] += derivative * mu
            if weight != 0:
                weights.append((coordinate, _finite_rational(weight, "Calibration-field weight")))
        return CalibrationFieldEvaluation(
            parameter_name=self.parameter_name,
            conditions=tuple((axis.name, normalized[axis.name]) for axis in self.axes),
            consensus_scale=_finite_rational(consensus, "Calibration-field consensus scale"),
            consensus_standard_uncertainty=_finite_rational(
                envelope, "Calibration-field consensus standard uncertainty"
            ),
            independent_nodes_standard_uncertainty=sqrt(
                _finite_rational(independent_variance, "Calibration-field independent variance")
            ),
            within_engine_process_variance=_finite_rational(
                within, "Calibration-field within-engine variance"
            ),
            between_engine_process_variance=_finite_rational(
                between, "Calibration-field between-engine variance"
            ),
            population_process_variance=_finite_rational(
                within + between, "Calibration-field population variance"
            ),
            consensus_scale_gradient=tuple(
                (axis.name, _finite_rational(value, f"Calibration-field gradient: {axis.name}"))
                for axis, value in zip(self.axes, gradient, strict=True)
            ),
            cell_lower_coordinates=tuple(
                axis.values[brackets[index][0]] for index, axis in enumerate(self.axes)
            ),
            cell_upper_coordinates=tuple(
                axis.values[brackets[index][1]] for index, axis in enumerate(self.axes)
            ),
            boundary_axes=tuple(boundary_axes),
            node_weights=tuple(weights),
        )

    def residual_uncertain_parameter(
        self,
        conditions: Mapping[str, float],
        *,
        scope: CalibrationScope = CalibrationScope.POPULATION,
        coverage_factor: float = 3.0,
        component_id: str | None = None,
    ) -> UncertainParameter:
        """Return the aleatory residual scale at an operating point.

        Apply the interpolated ``consensus_scale`` deterministically first.
        The returned truncated-normal factor is centered at one with relative
        standard deviation ``r = sigma(x) / mu(x)``, where ``sigma(x)**2`` is
        the interpolated population or within-engine process variance [1]_,
        and bounds ``[max(1e-12, 1 - k r), 1 + k r]`` — the same construction
        and defaults as ``HotFireScaleCalibration.residual_uncertain_parameter``
        [2]_.

        Parameters
        ----------
        conditions : Mapping of str to float
            In-domain operating point.
        scope : CalibrationScope, optional
            Population or same-engine residual variance.
        coverage_factor : float, optional
            Positive finite truncation half-width ``k`` in standard deviations.
        component_id : str or None, optional
            When given, the name is qualified as
            ``"<scale>::aleatory:<component_id>"`` so that it can coexist with
            other components of the same scale; otherwise the plain scale
            name is used.

        Returns
        -------
        UncertainParameter
            Aleatory truncated-normal scale factor.

        Raises
        ------
        DomainError
            If the point is outside the domain, ``scope`` or
            ``coverage_factor`` is invalid, or the residual variance is zero.

        References
        ----------
        .. [1] NIST/SEMATECH e-Handbook of Statistical Methods, section 7.4.4,
               "Variance Components" (one-way random-effects model).
               https://www.itl.nist.gov/div898/handbook/prc/section4/prc44.htm
        .. [2] NIST TN 1297, "Guidelines for Evaluating and Expressing the
               Uncertainty of NIST Measurement Results", 1994.
               https://www.nist.gov/pml/nist-technical-note-1297
        """

        if not isinstance(scope, CalibrationScope):
            raise DomainError("Calibration scope must be a CalibrationScope value.")
        coverage = _coverage(coverage_factor)
        evaluation = self.evaluate(conditions)
        variance = (
            evaluation.population_process_variance
            if scope is CalibrationScope.POPULATION
            else evaluation.within_engine_process_variance
        )
        relative = sqrt(variance) / evaluation.consensus_scale
        if relative <= 0.0:
            raise DomainError(
                "Interpolated residual variance is zero; no stochastic parameter is needed."
            )
        name = (
            self.parameter_name
            if component_id is None
            else qualified_scale_component_name(
                self.parameter_name,
                uncertainty_class=UncertaintyClass.ALEATORY,
                component_id=component_id,
            )
        )
        return UncertainParameter(
            name=name,
            nominal=1.0,
            lower=max(1e-12, 1.0 - coverage * relative),
            upper=1.0 + coverage * relative,
            distribution=DistributionKind.TRUNCATED_NORMAL,
            uncertainty_class=UncertaintyClass.ALEATORY,
            source=f"{self.source}; hot-fire calibration field sha256={self.field_hash}",
            rationale=(
                f"{scope.value} residual variability at {_conditions_text(evaluation)} after "
                f"applying interpolated consensus correction "
                f"{evaluation.consensus_scale:.17g}; {self.rationale}"
            ),
            standard_deviation=relative,
        )

    def consensus_uncertain_parameter(
        self,
        conditions: Mapping[str, float],
        *,
        coverage_factor: float = 3.0,
        component_id: str = "hot-fire-field-consensus",
    ) -> UncertainParameter:
        """Return the epistemic factor for the interpolated consensus uncertainty.

        The factor is centered at one with relative standard deviation
        ``u(x) / mu(x)`` (the linear-envelope uncertainty [1]_) and the same
        truncation as :meth:`residual_uncertain_parameter`. The uncertainty
        of an estimated correction is reducible with more data, so it is
        classified as epistemic [2]_ and named
        ``"<scale>::epistemic:<component_id>"``.

        Parameters
        ----------
        conditions : Mapping of str to float
            In-domain operating point.
        coverage_factor : float, optional
            Positive finite truncation half-width in standard deviations.
        component_id : str, optional
            Component identifier without ``"::"``.

        Returns
        -------
        UncertainParameter
            Epistemic truncated-normal scale factor.

        Raises
        ------
        DomainError
            If the point is outside the domain or an argument is invalid.

        References
        ----------
        .. [1] JCGM-100-2008-SEC5. JCGM 100:2008, sections 5.1.6 (equation 12)
               and 5.2.2 (NOTE 1).
               https://www.bipm.org/documents/20126/2071204/JCGM_100_2008_E.pdf
        .. [2] NASA/SP-2009-569, "Bayesian Inference for NASA Probabilistic
               Risk and Reliability Analysis".
               https://ntrs.nasa.gov/citations/20090023159
        """

        coverage = _coverage(coverage_factor)
        evaluation = self.evaluate(conditions)
        relative = evaluation.relative_standard_uncertainty
        return UncertainParameter(
            name=qualified_scale_component_name(
                self.parameter_name,
                uncertainty_class=UncertaintyClass.EPISTEMIC,
                component_id=component_id,
            ),
            nominal=1.0,
            lower=max(1e-12, 1.0 - coverage * relative),
            upper=1.0 + coverage * relative,
            distribution=DistributionKind.TRUNCATED_NORMAL,
            uncertainty_class=UncertaintyClass.EPISTEMIC,
            source=f"{self.source}; hot-fire calibration field sha256={self.field_hash}",
            rationale=(
                f"Consensus standard uncertainty (linear envelope over grid nodes) at "
                f"{_conditions_text(evaluation)}; {self.rationale}"
            ),
            standard_deviation=relative,
        )

    def ensemble_parameters(
        self,
        conditions: Mapping[str, float],
        *,
        scope: CalibrationScope = CalibrationScope.POPULATION,
        coverage_factor: float = 3.0,
    ) -> tuple[UncertainParameter, UncertainParameter]:
        """Return the aleatory residual and epistemic consensus factors.

        Parameters
        ----------
        conditions : Mapping of str to float
            In-domain operating point.
        scope : CalibrationScope, optional
            Residual population.
        coverage_factor : float, optional
            Positive finite truncation half-width in standard deviations.

        Returns
        -------
        tuple of UncertainParameter
            ``"<scale>::aleatory:hot-fire-field"`` and
            ``"<scale>::epistemic:hot-fire-field-consensus"``, ready for a
            ``ScenarioConditionalUncertainty`` whose artifact already carries
            the interpolated consensus correction. The two factors are
            classified separately [1]_ [2]_.

        Raises
        ------
        DomainError
            As :meth:`residual_uncertain_parameter`.

        References
        ----------
        .. [1] NASA-STD-7009B, "Standard for Models and Simulations".
               https://standards.nasa.gov/standard/NASA/NASA-STD-7009
        .. [2] NASA/SP-2009-569, "Bayesian Inference for NASA Probabilistic
               Risk and Reliability Analysis".
               https://ntrs.nasa.gov/citations/20090023159
        """

        return (
            self.residual_uncertain_parameter(
                conditions,
                scope=scope,
                coverage_factor=coverage_factor,
                component_id="hot-fire-field",
            ),
            self.consensus_uncertain_parameter(conditions, coverage_factor=coverage_factor),
        )

    def to_dict(self, *, include_hash: bool = True) -> dict[str, Any]:
        """Return the deterministic JSON-compatible field document.

        Parameters
        ----------
        include_hash : bool, default True
            Whether to add ``field_hash``.

        Returns
        -------
        dict
            Schema, parameter, provenance, axes, nodes sorted by coordinates
            (each with its complete calibration record and calibration hash),
            reference identifiers and the numerical-policy identifiers.

        Notes
        -----
        The document is versioned and unit-explicit and records its evidence
        [1]_; its encoding is standard JSON [2]_.

        References
        ----------
        .. [1] NASA-STD-7009B, "Standard for Models and Simulations".
               https://standards.nasa.gov/standard/NASA/NASA-STD-7009
        .. [2] IETF RFC 8259, "The JavaScript Object Notation (JSON) Data
               Interchange Format". https://www.rfc-editor.org/rfc/rfc8259
        """

        payload: dict[str, Any] = {
            "schema": self.schema,
            "parameter_name": self.parameter_name,
            "source": self.source,
            "rationale": self.rationale,
            "axes": [
                {"name": axis.name, "unit": axis.unit, "values": list(axis.values)}
                for axis in self.axes
            ],
            "nodes": [
                {
                    "coordinates": list(node.coordinates),
                    "calibration": json.loads(canonical_json_bytes(node.calibration)),
                    "calibration_hash": node.calibration.calibration_hash,
                }
                for node in sorted(self.nodes, key=lambda item: item.coordinates)
            ],
            "reference_ids": list(self.reference_ids),
            "interpolation_numerics": _MULTILINEAR_POLICY,
            "uncertainty_numerics": CALIBRATION_FIELD_UNCERTAINTY_POLICY,
        }
        if include_hash:
            payload["field_hash"] = canonical_sha256(payload)
        return payload

    def to_json(self, *, indent: int | None = 2) -> str:
        """Serialize the complete hash-bound field as JSON.

        Parameters
        ----------
        indent : int or None, default 2
            Indentation passed to ``json.dumps``.

        Returns
        -------
        str
            Sorted JSON without NaN or Infinity tokens [1]_ [2]_.

        References
        ----------
        .. [1] IETF RFC 8259, "The JavaScript Object Notation (JSON) Data
               Interchange Format". https://www.rfc-editor.org/rfc/rfc8259
        .. [2] PYTHON-3.12-JSON. Python 3.12, ``json``: encoder/decoder hooks,
               duplicate-member and nonstandard numeric token behavior.
        """

        return json.dumps(
            self.to_dict(), allow_nan=False, ensure_ascii=False, indent=indent, sort_keys=True
        )

    @classmethod
    def from_dict(cls, payload: Mapping[str, Any]) -> HotFireCalibrationField:
        """Rebuild a field, rerun every node calibration and verify the hash.

        Parameters
        ----------
        payload : Mapping
            Document produced by :meth:`to_dict`.

        Returns
        -------
        HotFireCalibrationField
            The reconstructed field.

        Raises
        ------
        InputError
            If the document is malformed or has another schema, a node
            calibration does not reproduce from its stored measurements, the
            field hash does not match, the document differs from the rebuilt
            one, or a numerical-policy identifier is missing or differs.

        Notes
        -----
        Every node is recalibrated with ``calibrate_hot_fire_scale`` [1]_
        from its stored measurements, source, rationale and tolerance, so a
        falsified derived result is refused even when its hash is updated.

        References
        ----------
        .. [1] NIST Dataplot, "Consensus Mean" (heteroscedastic random effects
               and Mandel-Paule estimation).
               https://www.itl.nist.gov/div898/software/dataplot/refman1/auxillar/consmean.htm
        """

        if not isinstance(payload, Mapping):
            raise InputError("Calibration-field root must be an object.")
        if payload.get("schema") != CALIBRATION_FIELD_SCHEMA:
            raise InputError("Unsupported calibration-field schema.")
        try:
            axes = tuple(
                PerformanceAxis(name=item["name"], unit=item["unit"], values=item["values"])
                for item in payload["axes"]
            )
            nodes = tuple(
                CalibrationFieldNode(
                    coordinates=item["coordinates"],
                    calibration=_reproduce_calibration(item["calibration"]),
                )
                for item in payload["nodes"]
            )
            rebuilt = cls(
                parameter_name=payload["parameter_name"],
                axes=axes,
                nodes=nodes,
                source=payload["source"],
                rationale=payload["rationale"],
                reference_ids=payload["reference_ids"],
                schema=payload["schema"],
            )
        except (*DOCUMENT_ERRORS, DomainError) as error:
            raise InputError("Malformed calibration-field artifact.") from error
        if payload.get("field_hash") != rebuilt.field_hash:
            raise InputError("Calibration-field hash does not match its contents.")
        _require_numerical_policies(payload, rebuilt.to_dict(include_hash=False))
        if canonical_json_bytes(dict(payload)) != canonical_json_bytes(rebuilt.to_dict()):
            raise InputError("Calibration-field node calibrations do not reproduce.")
        return rebuilt

    @classmethod
    def from_json(cls, text: str) -> HotFireCalibrationField:
        """Parse a strict JSON field document [1]_ and check it with :meth:`from_dict`.

        Parameters
        ----------
        text : str
            JSON text produced by :meth:`to_json`.

        Returns
        -------
        HotFireCalibrationField
            The reconstructed field.

        Raises
        ------
        InputError
            If the text is not a valid JSON object or :meth:`from_dict`
            rejects it.

        References
        ----------
        .. [1] IETF RFC 8259, "The JavaScript Object Notation (JSON) Data
               Interchange Format". https://www.rfc-editor.org/rfc/rfc8259
        """

        return _read_document(text, "Calibration-field artifact", cls.from_dict)


def _coverage(coverage_factor: float) -> float:
    """Validate a positive finite coverage factor."""
    value = real(coverage_factor, "Calibration-field coverage factor")
    if not isfinite(value) or value <= 0.0:
        raise DomainError("Calibration-field coverage factor must be positive and finite.")
    return value


def _conditions_text(evaluation: CalibrationFieldEvaluation) -> str:
    """Return ``name=value`` pairs of the evaluated conditions."""
    return ", ".join(f"{name}={value:.17g}" for name, value in evaluation.conditions)


def _reproduce_calibration(record: Mapping[str, Any]) -> HotFireScaleCalibration:
    """Rerun one stored node calibration and require identical content.

    The Mandel-Paule iteration limit is not stored; any limit at least the
    recorded iteration count reproduces the same bracket and bisection path.
    """
    measurements = tuple(
        HotFireScaleMeasurement(
            engine_id=item["engine_id"],
            run_id=item["run_id"],
            condition_id=item["condition_id"],
            value=item["value"],
            standard_uncertainty=item["standard_uncertainty"],
            source=item["source"],
        )
        for item in record["measurements"]
    )
    evidence = record["evidence"]
    rebuilt = calibrate_hot_fire_scale(
        record["parameter_name"],
        measurements,
        source=record["source"],
        rationale=record["rationale"],
        relative_tolerance=float(evidence["relative_tolerance"]),
        maximum_iterations=max(200, int(evidence["iterations"]) + 1),
    )
    if canonical_json_bytes(rebuilt) != canonical_json_bytes(dict(record)):
        raise DomainError("Calibration-field node calibration does not reproduce.")
    return rebuilt


def calibrate_hot_fire_field(
    parameter_name: str,
    axes: tuple[PerformanceAxis, ...],
    condition_coordinates: Mapping[str, tuple[float, ...]],
    measurements: tuple[HotFireScaleMeasurement, ...],
    *,
    source: str,
    rationale: str,
    relative_tolerance: float = 1e-10,
    maximum_iterations: int = 200,
) -> HotFireCalibrationField:
    """Calibrate every operating-condition stratum and assemble the field.

    Parameters
    ----------
    parameter_name : str
        Calibrated scale name.
    axes : tuple of PerformanceAxis
        Grid axes.
    condition_coordinates : Mapping of str to tuple of float
        Grid-node coordinates of every ``condition_id``.
    measurements : tuple of HotFireScaleMeasurement
        Hot-fire measurements; each ``condition_id`` is one stratum [1]_.
    source, rationale : str
        Non-blank provenance and justification, recorded on every node.
    relative_tolerance : float, optional
        Mandel-Paule tolerance passed to ``calibrate_hot_fire_scale`` [2]_.
    maximum_iterations : int, optional
        Mandel-Paule iteration limit passed to ``calibrate_hot_fire_scale``.

    Returns
    -------
    HotFireCalibrationField
        The complete field.

    Raises
    ------
    DomainError
        If a measurement's condition has no coordinates, a mapped condition
        has no measurements, two conditions share coordinates, a stratum
        cannot be calibrated, or the nodes do not form the complete grid
        (missing, off-grid or out-of-grid nodes).

    References
    ----------
    .. [1] NIST/SEMATECH e-Handbook of Statistical Methods, section 7.4.4,
           "Variance Components" (one-way random-effects model).
           https://www.itl.nist.gov/div898/handbook/prc/section4/prc44.htm
    .. [2] NIST Dataplot, "Consensus Mean" (heteroscedastic random effects and
           Mandel-Paule estimation).
           https://www.itl.nist.gov/div898/software/dataplot/refman1/auxillar/consmean.htm
    """

    if not isinstance(condition_coordinates, Mapping):
        raise DomainError("Condition coordinates must be a mapping of condition_id to node.")
    measurements = records(measurements, HotFireScaleMeasurement, "Calibration-field measurements")
    coordinates = {
        text(condition_id, "Calibration-field condition_id"): real_tuple(
            value, "Calibration-field condition coordinates"
        )
        for condition_id, value in condition_coordinates.items()
    }
    unmapped = sorted({item.condition_id for item in measurements} - set(coordinates))
    if unmapped:
        raise DomainError(f"Measurements at unmapped operating conditions: {unmapped}.")
    nodes_by_coordinate: dict[tuple[float, ...], str] = {}
    for condition_id, node in coordinates.items():
        if node in nodes_by_coordinate:
            raise DomainError(
                f"Duplicate calibration-field node {node} for conditions "
                f"{nodes_by_coordinate[node]!r} and {condition_id!r}."
            )
        nodes_by_coordinate[node] = condition_id
    nodes = []
    for condition_id, node in coordinates.items():
        stratum = tuple(item for item in measurements if item.condition_id == condition_id)
        if not stratum:
            raise DomainError(f"Operating condition {condition_id!r} has no measurements.")
        nodes.append(
            CalibrationFieldNode(
                coordinates=node,
                calibration=calibrate_hot_fire_scale(
                    parameter_name,
                    stratum,
                    source=source,
                    rationale=rationale,
                    relative_tolerance=relative_tolerance,
                    maximum_iterations=maximum_iterations,
                ),
            )
        )
    return HotFireCalibrationField(
        parameter_name=parameter_name.strip(),
        axes=tuple(axes),
        nodes=tuple(nodes),
        source=source,
        rationale=rationale,
    )


@dataclass(frozen=True, slots=True)
class CalibratedSurfaceEvaluation:
    """Corrected performance-surface evaluation with field uncertainties.

    Attributes
    ----------
    surface : PerformanceSurfaceEvaluation
        Evaluation of the corrected surface.
    thrust, specific_impulse : CalibrationFieldEvaluation or None
        Field evaluations at the same point, when a field was applied.
    force_scale_relative_standard_uncertainty : float
        ``u_F / mu_F`` of the thrust field; zero without one [1]_.
    specific_impulse_relative_standard_uncertainty : float
        ``u_I / mu_I`` of the specific-impulse field; zero without one.
    mass_flow_scale_relative_standard_uncertainty : float
        First-order relative uncertainty of ``mu_F / mu_I``. With unknown
        correlation between the two fields the bound
        ``u_F / mu_F + u_I / mu_I`` is reported (JCGM 100:2008 5.1.6 and
        5.2.2 NOTE 4 with ``|r| <= 1``) [1]_.

    References
    ----------
    .. [1] JCGM-100-2008-SEC5. JCGM 100:2008, sections 5.1.6 (equation 12,
           NOTE 2) and 5.2.2 (NOTE 4).
           https://www.bipm.org/documents/20126/2071204/JCGM_100_2008_E.pdf
    """

    surface: PerformanceSurfaceEvaluation
    thrust: CalibrationFieldEvaluation | None
    specific_impulse: CalibrationFieldEvaluation | None
    force_scale_relative_standard_uncertainty: float
    specific_impulse_relative_standard_uncertainty: float
    mass_flow_scale_relative_standard_uncertainty: float


@dataclass(frozen=True, slots=True)
class CalibratedPerformanceSurface:
    """Performance surface corrected node-wise by calibration fields.

    Instances are produced by :func:`apply_calibration_fields`.

    Parameters
    ----------
    base_surface_hash : str
        ``surface_hash`` of the uncorrected surface.
    surface : RectilinearPerformanceSurface
        Corrected surface, interpolated with the base surface's own rule.
    thrust_field, specific_impulse_field : HotFireCalibrationField or None
        Applied fields; at least one is present and both share the surface
        axes.

    Raises
    ------
    DomainError
        If no field is present, a field calibrates another parameter, or a
        field's axes differ from the surface axes.

    Notes
    -----
    Thrust scales with the thrust consensus and mass flow with
    ``thrust / specific impulse`` because ``mdot = F / (g0 Isp)`` [1]_. The
    correction is applied at the nodes and the corrected nodes are
    interpolated; between nodes this differs at second order from the product
    of the two interpolants. The validated domain is the shared grid [2]_.

    References
    ----------
    .. [1] NASA Glenn, "Specific Impulse".
           https://www1.grc.nasa.gov/beginners-guide-to-aeronautics/specific-impulse/
    .. [2] NASA-STD-7009B, "Standard for Models and Simulations".
           https://standards.nasa.gov/standard/NASA/NASA-STD-7009
    """

    base_surface_hash: str
    surface: RectilinearPerformanceSurface
    thrust_field: HotFireCalibrationField | None
    specific_impulse_field: HotFireCalibrationField | None

    def __post_init__(self) -> None:
        if not isinstance(self.surface, RectilinearPerformanceSurface):
            raise DomainError("Calibrated surface requires a RectilinearPerformanceSurface.")
        _check_fields(self.surface.axes, self.thrust_field, self.specific_impulse_field)

    def validity_surfaces(self) -> tuple[PerformanceValiditySurface, ...]:
        """Return the stopping domain events of the corrected surface [1]_.

        References
        ----------
        .. [1] NASA-STD-7009B, "Standard for Models and Simulations".
               https://standards.nasa.gov/standard/NASA/NASA-STD-7009
        """

        return self.surface.validity_surfaces()

    def evaluate(self, conditions: Mapping[str, float]) -> CalibratedSurfaceEvaluation:
        """Evaluate the corrected surface and the field uncertainties.

        Parameters
        ----------
        conditions : Mapping of str to float
            In-domain operating point.

        Returns
        -------
        CalibratedSurfaceEvaluation
            Corrected scales, field evaluations and relative uncertainties
            [1]_.

        Raises
        ------
        DomainError
            If the point is outside the shared domain.

        References
        ----------
        .. [1] JCGM-100-2008-SEC5. JCGM 100:2008, sections 5.1.6 (equation
               12) and 5.2.2 (NOTE 4).
               https://www.bipm.org/documents/20126/2071204/JCGM_100_2008_E.pdf
        """

        surface = self.surface.evaluate(conditions)
        thrust = None if self.thrust_field is None else self.thrust_field.evaluate(conditions)
        isp = (
            None
            if self.specific_impulse_field is None
            else self.specific_impulse_field.evaluate(conditions)
        )
        force_relative = 0.0 if thrust is None else thrust.relative_standard_uncertainty
        isp_relative = 0.0 if isp is None else isp.relative_standard_uncertainty
        return CalibratedSurfaceEvaluation(
            surface=surface,
            thrust=thrust,
            specific_impulse=isp,
            force_scale_relative_standard_uncertainty=force_relative,
            specific_impulse_relative_standard_uncertainty=isp_relative,
            mass_flow_scale_relative_standard_uncertainty=force_relative + isp_relative,
        )

    def ensemble_parameters(
        self,
        conditions: Mapping[str, float],
        *,
        scope: CalibrationScope = CalibrationScope.POPULATION,
        coverage_factor: float = 3.0,
    ) -> tuple[UncertainParameter, ...]:
        """Return the residual and consensus factors of every applied field.

        Parameters
        ----------
        conditions : Mapping of str to float
            In-domain operating point.
        scope : CalibrationScope, optional
            Residual population.
        coverage_factor : float, optional
            Positive finite truncation half-width in standard deviations.

        Returns
        -------
        tuple of UncertainParameter
            Thrust factors then specific-impulse factors, as returned by
            ``HotFireCalibrationField.ensemble_parameters``; their mutual
            correlation is not modelled here (independence by default in
            the ensemble) [1]_.

        Raises
        ------
        DomainError
            As ``HotFireCalibrationField.residual_uncertain_parameter``.

        References
        ----------
        .. [1] NASA-STD-7009B, "Standard for Models and Simulations".
               https://standards.nasa.gov/standard/NASA/NASA-STD-7009
        """

        return tuple(
            parameter
            for item in (self.thrust_field, self.specific_impulse_field)
            if item is not None
            for parameter in item.ensemble_parameters(
                conditions, scope=scope, coverage_factor=coverage_factor
            )
        )


def _check_fields(
    axes: tuple[PerformanceAxis, ...],
    thrust_field: HotFireCalibrationField | None,
    specific_impulse_field: HotFireCalibrationField | None,
) -> None:
    """Require at least one field, the right parameters and identical axes."""
    if thrust_field is None and specific_impulse_field is None:
        raise DomainError("At least one calibration field is required.")
    for item, expected in (
        (thrust_field, "thrust_scale"),
        (specific_impulse_field, "specific_impulse_scale"),
    ):
        if item is None:
            continue
        if not isinstance(item, HotFireCalibrationField) or item.parameter_name != expected:
            raise DomainError(f"The {expected} field must be a {expected} calibration field.")
        if item.axes != axes:
            raise DomainError(
                f"The {expected} field axes (names, units and coordinates) must equal the "
                "performance-surface axes."
            )


def apply_calibration_fields(
    surface: RectilinearPerformanceSurface,
    *,
    thrust_field: HotFireCalibrationField | None = None,
    specific_impulse_field: HotFireCalibrationField | None = None,
) -> CalibratedPerformanceSurface:
    """Correct a performance surface with thrust and/or Isp calibration fields.

    Parameters
    ----------
    surface : RectilinearPerformanceSurface
        Uncorrected force and tank-flow multiplier map.
    thrust_field, specific_impulse_field : HotFireCalibrationField or None
        ``thrust_scale`` and ``specific_impulse_scale`` fields on exactly the
        surface axes; at least one is required.

    Returns
    -------
    CalibratedPerformanceSurface
        Corrected surface (node force ``s_F mu_F``, node flow
        ``s_mdot mu_F / mu_I``, with ``mdot = F / (g0 Isp)`` [1]_), the
        fields, and the base surface hash. The corrected surface keeps the
        base interpolation mode; its source hash binds the base surface hash
        and the field hashes [2]_.

    Raises
    ------
    DomainError
        If the surface or fields are invalid or the axes differ.

    References
    ----------
    .. [1] NASA Glenn, "Specific Impulse".
           https://www1.grc.nasa.gov/beginners-guide-to-aeronautics/specific-impulse/
    .. [2] NASA-STD-7009B, "Standard for Models and Simulations".
           https://standards.nasa.gov/standard/NASA/NASA-STD-7009
    """

    if not isinstance(surface, RectilinearPerformanceSurface):
        raise DomainError("apply_calibration_fields requires a RectilinearPerformanceSurface.")
    _check_fields(surface.axes, thrust_field, specific_impulse_field)
    points = []
    for point in surface.points:
        thrust = (
            Fraction(1)
            if thrust_field is None
            else _fraction(thrust_field.node_calibration(point.coordinates).consensus_scale)
        )
        isp = (
            Fraction(1)
            if specific_impulse_field is None
            else _fraction(
                specific_impulse_field.node_calibration(point.coordinates).consensus_scale
            )
        )
        points.append(
            PerformanceScalePoint(
                coordinates=point.coordinates,
                force_scale=_finite_rational(
                    _fraction(point.force_scale) * thrust, "Calibrated force scale"
                ),
                mass_flow_scale=_finite_rational(
                    _fraction(point.mass_flow_scale) * thrust / isp, "Calibrated mass-flow scale"
                ),
            )
        )
    base_hash = surface.surface_hash
    provenance = {
        "base_surface_hash": base_hash,
        "thrust_field_hash": None if thrust_field is None else thrust_field.field_hash,
        "specific_impulse_field_hash": (
            None if specific_impulse_field is None else specific_impulse_field.field_hash
        ),
    }
    references = surface.reference_ids + tuple(
        reference
        for item in (thrust_field, specific_impulse_field)
        if item is not None
        for reference in item.reference_ids
    )
    corrected = RectilinearPerformanceSurface(
        axes=surface.axes,
        points=tuple(points),
        source_id=f"calibrated:{surface.source_id}",
        source_sha256=canonical_sha256(provenance),
        reference_ids=tuple(dict.fromkeys(references)),
        interpolation=surface.interpolation,
    )
    return CalibratedPerformanceSurface(
        base_surface_hash=base_hash,
        surface=corrected,
        thrust_field=thrust_field,
        specific_impulse_field=specific_impulse_field,
    )
