"""Gridded aerodynamic force-and-moment database with an SLS-style uncertainty model.

The NASA Space Launch System (SLS) ascent force-and-moment database is
delivered as look-up tables of six body-axis coefficients at breakpoints in
Mach number, angle of attack and sideslip, built from response surfaces fitted
to wind-tunnel and CFD data [1]_. Each query returns a nominal value ``C_Q``
and a centred uncertainty bound ``U`` such that
``C_Q - U <= C_true <= C_Q + U`` [2]_ (equation 1), a +-3-sigma-equivalent
(99.7%) bound [1]_ obtained by root-sum-squaring independent terms and
multiplying by an uncertainty factor [1]_ [2]_ [3]_.

This module provides

* :class:`AeroDatabase` - an N-dimensional rectilinear table with
  multilinear interpolation that refuses queries outside the grid;
* uncertainty attached either as a gridded table or as piecewise-constant
  ranges of one axis (the SLS ascent uncertainty is a function of Mach range
  only [1]_ [2]_);
* the uncertainty build-up arithmetic of the cited SLS papers: range-method
  standard deviation, 3-sigma repeatability, 99.7% coverage bounds, net
  modelling uncertainty and the factored root-sum-square;
* the linear slender-body-theory fits Hemsch used to compare CFD and
  wind-tunnel databases [4]_.

Interpolation choice
--------------------
The cited papers state that the trajectory codes interpolate the delivered
tables, and that this interpolation can differ from the response surface
(an uncertainty term where significant) [3]_, but do not state the
interpolant. Multilinear interpolation is this module's choice; the
response-surface construction itself (natural-neighbour interpolation of
scattered data, kriging) is not implemented.

References
----------
.. [1] SHEA-2020-SLS-ASCENT-DB. P. R. Shea, J. T. Pinier, H. P. Houlden,
       A. L. Favaregh, M. J. Hemsch, D. J. Dalle, S. E. Rogers, J. G. Meeroff
       and H. C. Lee, "Ascent Aerodynamic Force and Moment Database
       Development for the Space Launch System", NTRS 20200002382.
.. [2] FAVAREGH-2016-SLS-UQ. A. L. Favaregh, H. P. Houlden and J. T. Pinier,
       "Quantification of the Uncertainties for the Space Launch System
       Liftoff/Transition and Ascent Databases", AIAA 2016-0795; NTRS
       20160007692.
.. [3] LEE-2024-ARTEMIS-UQ. M. W. Lee, P. R. Shea, J. T. Pinier and
       D. T. Chan, "Uncertainty Quantification of Artemis I Space Launch
       System Integrated Aerodynamics Databases", Journal of Spacecraft and
       Rockets; NTRS 20240000054.
.. [4] HEMSCH-2016-SLS-CFD-WT. M. J. Hemsch, "Global Comparison of CFD and
       Wind-Tunnel-Derived Force and Moment Databases for the Space Launch
       System", AIAA 2016-0794; NTRS 20160007695.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass, field
from math import ceil, isfinite, sqrt

from rocket_propulsion.core.errors import DomainError, InputError
from rocket_propulsion.core.metadata import WarningMessage

#: Range-method bias-correction factors ``d2`` for samples of n = 2..5
#: (Lee et al., NTRS 20240000054, below equation 2, citing Montgomery).
RANGE_METHOD_D2 = {2: 1.128, 3: 1.693, 4: 2.059, 5: 2.326}

#: Coverage used by the SLS papers for 3-sigma-equivalent bounds.
SLS_COVERAGE = 0.997


def _finite(value: float, field_name: str) -> float:
    """Return ``value`` as a finite float.

    References
    ----------
    .. [1] SHEA-2020-SLS-ASCENT-DB. P. R. Shea et al., NTRS 20200002382.
    """

    if not isinstance(value, (int, float)) or isinstance(value, bool):
        raise InputError(f"{field_name} must be a number.", field=field_name)
    number = float(value)
    if not isfinite(number):
        raise InputError(f"{field_name} must be finite.", field=field_name)
    return number


@dataclass(frozen=True, slots=True)
class RangeUncertainty:
    """Piecewise-constant uncertainty over ranges of one database axis.

    Attributes
    ----------
    axis : str
        Name of the axis the ranges refer to (for example ``"mach"``).
    ranges : tuple of (float, float, float)
        ``(lower, upper, bound)`` triples, inclusive at both ends, ordered
        and non-overlapping. A query outside every range is refused: the SLS
        ascent uncertainties were developed for Mach 0.5-1.3 and 1.6-2.5
        separately [1]_.

    References
    ----------
    .. [1] SHEA-2020-SLS-ASCENT-DB. P. R. Shea et al., NTRS 20200002382,
           section VII.
    """

    axis: str
    ranges: tuple[tuple[float, float, float], ...]

    def __post_init__(self) -> None:
        if not isinstance(self.axis, str) or not self.axis:
            raise InputError("axis must be a non-empty string.", field="axis")
        if not isinstance(self.ranges, tuple) or not self.ranges:
            raise InputError("ranges must be a non-empty tuple.", field="ranges")
        previous_upper = None
        for item in self.ranges:
            if not isinstance(item, tuple) or len(item) != 3:
                raise InputError("Each range is (lower, upper, bound).", field="ranges")
            lower, upper, bound = (_finite(v, "ranges") for v in item)
            if upper < lower or bound < 0.0:
                raise InputError("Ranges need lower <= upper and bound >= 0.", field="ranges")
            if previous_upper is not None and lower <= previous_upper:
                raise InputError("Ranges must be ordered and non-overlapping.", field="ranges")
            previous_upper = upper

    def bound_at(self, value: float) -> float:
        """Return the bound of the range containing ``value``.

        Raises
        ------
        DomainError
            If ``value`` is outside every range.

        References
        ----------
        .. [1] SHEA-2020-SLS-ASCENT-DB. P. R. Shea et al., NTRS 20200002382.
        """

        number = _finite(value, self.axis)
        for lower, upper, bound in self.ranges:
            if lower <= number <= upper:
                return float(bound)
        raise DomainError(
            f"{self.axis} = {number:g} lies outside every uncertainty range.", field=self.axis
        )


@dataclass(frozen=True, slots=True)
class AeroQuery:
    """Result of one database query.

    Attributes
    ----------
    point : dict of str to float
        The query coordinates.
    nominal : dict of str to float
        Interpolated nominal coefficients ``C_Q``.
    uncertainty : dict of str to float
        Bounds ``U`` (only for coefficients that have one).
    lower, upper : dict of str to float
        ``C_Q - U`` and ``C_Q + U`` (equation 1 of [1]_).
    warnings : tuple of WarningMessage
        Notes, e.g. on interpolated uncertainty.

    References
    ----------
    .. [1] FAVAREGH-2016-SLS-UQ. A. L. Favaregh et al., NTRS 20160007692,
           equation 1.
    """

    point: dict[str, float]
    nominal: dict[str, float]
    uncertainty: dict[str, float]
    lower: dict[str, float]
    upper: dict[str, float]
    warnings: tuple[WarningMessage, ...] = ()

    def dispersed(self, multipliers: Mapping[str, float]) -> dict[str, float]:
        """Return ``C_Q + k U`` for user-chosen multipliers ``-1 <= k <= 1``.

        Parameters
        ----------
        multipliers : mapping of str to float
            One ``k`` per coefficient to disperse; others stay nominal.

        Returns
        -------
        dict of str to float
            Dispersed coefficients.

        Raises
        ------
        DomainError
            If ``|k| > 1`` (outside the stated bound) or a coefficient without
            uncertainty is dispersed.

        Notes
        -----
        The papers define the bound [1]_ and state that the databases feed
        Monte Carlo trajectory simulations [2]_; they do not state how ``k``
        is sampled or whether coefficients are dispersed together, so the
        caller chooses ``k``.

        References
        ----------
        .. [1] FAVAREGH-2016-SLS-UQ. A. L. Favaregh et al., NTRS 20160007692,
               equation 1.
        .. [2] LEE-2024-ARTEMIS-UQ. M. W. Lee et al., NTRS 20240000054,
               section I.
        """

        result = dict(self.nominal)
        for name, value in multipliers.items():
            k = _finite(value, "multipliers")
            if abs(k) > 1.0:
                raise DomainError("Dispersion multipliers must satisfy |k| <= 1.", field=name)
            if name not in self.uncertainty:
                raise DomainError(f"Coefficient {name!r} has no uncertainty.", field=name)
            result[name] = self.nominal[name] + k * self.uncertainty[name]
        return result


@dataclass(frozen=True, slots=True)
class AeroDatabase:
    """Rectilinear N-dimensional coefficient tables with uncertainty [1]_.

    Attributes
    ----------
    axes : tuple of str
        Axis names in table order, e.g. ``("mach", "alpha_deg", "beta_deg")``.
    breakpoints : tuple of tuple of float
        Strictly increasing breakpoints per axis (at least two each).
    coefficients : dict of str to tuple of float
        Row-major (last axis fastest) nominal values per coefficient.
    gridded_uncertainty : dict of str to tuple of float
        Optional row-major uncertainty bounds on the same grid.
    range_uncertainty : dict of str to RangeUncertainty
        Optional piecewise-constant bounds over one axis.
    reference_area_m2, reference_length_m : float or None
        Reference quantities the coefficients are normalised with.
    moment_reference : str
        Description of the moment reference point.

    References
    ----------
    .. [1] SHEA-2020-SLS-ASCENT-DB. P. R. Shea et al., NTRS 20200002382,
           section V.B (look-up tables at Mach, alpha, beta breakpoints).
    """

    axes: tuple[str, ...]
    breakpoints: tuple[tuple[float, ...], ...]
    coefficients: dict[str, tuple[float, ...]]
    gridded_uncertainty: dict[str, tuple[float, ...]] = field(default_factory=dict)
    range_uncertainty: dict[str, RangeUncertainty] = field(default_factory=dict)
    reference_area_m2: float | None = None
    reference_length_m: float | None = None
    moment_reference: str = ""

    def __post_init__(self) -> None:
        if not isinstance(self.axes, tuple) or not self.axes:
            raise InputError("axes must be a non-empty tuple.", field="axes")
        if len(set(self.axes)) != len(self.axes):
            raise InputError("Axis names must be unique.", field="axes")
        if not isinstance(self.breakpoints, tuple) or len(self.breakpoints) != len(self.axes):
            raise InputError("One breakpoint tuple per axis is required.", field="breakpoints")
        size = 1
        for name, points in zip(self.axes, self.breakpoints, strict=True):
            if not isinstance(points, tuple) or len(points) < 2:
                raise InputError(f"Axis {name!r} needs at least two breakpoints.", field="breakpoints")
            values = [_finite(v, "breakpoints") for v in points]
            if any(b <= a for a, b in zip(values, values[1:], strict=False)):
                raise InputError(f"Axis {name!r} breakpoints must increase.", field="breakpoints")
            size *= len(points)
        if not self.coefficients:
            raise InputError("At least one coefficient table is required.", field="coefficients")
        for table_name, table in (("coefficients", self.coefficients), ("gridded_uncertainty", self.gridded_uncertainty)):
            for name, values in table.items():
                if not isinstance(values, tuple) or len(values) != size:
                    raise InputError(f"{table_name}[{name!r}] needs {size} values.", field=table_name)
                for value in values:
                    number = _finite(value, table_name)
                    if table_name == "gridded_uncertainty" and number < 0.0:
                        raise InputError("Uncertainty bounds must be non-negative.", field=table_name)
        for name in list(self.gridded_uncertainty) + list(self.range_uncertainty):
            if name not in self.coefficients:
                raise InputError(f"Uncertainty given for unknown coefficient {name!r}.", field="coefficients")
        for name, ranged in self.range_uncertainty.items():
            if name in self.gridded_uncertainty:
                raise InputError(f"Coefficient {name!r} has two uncertainty models.", field=name)
            if not isinstance(ranged, RangeUncertainty) or ranged.axis not in self.axes:
                raise InputError("Range uncertainty must name a database axis.", field=name)
        for value in (self.reference_area_m2, self.reference_length_m):
            if value is not None and (not isfinite(value) or value <= 0.0):
                raise InputError("Reference quantities must be positive.", field="reference_area_m2")

    def _weights(self, point: Mapping[str, float]) -> list[tuple[int, float]]:
        """Return flat indices and multilinear weights of the enclosing cell.

        References
        ----------
        .. [1] SHEA-2020-SLS-ASCENT-DB. P. R. Shea et al., NTRS 20200002382.
        """

        missing = [name for name in self.axes if name not in point]
        extra = [name for name in point if name not in self.axes]
        if missing or extra:
            raise InputError(f"Query needs exactly the axes {self.axes}.", field="point")
        per_axis = []
        for name, points in zip(self.axes, self.breakpoints, strict=True):
            value = _finite(point[name], name)
            if value < points[0] or value > points[-1]:
                raise DomainError(
                    f"{name} = {value:g} lies outside the database [{points[0]:g}, {points[-1]:g}].",
                    field=name,
                )
            upper = 1
            while upper < len(points) - 1 and points[upper] < value:
                upper += 1
            lower = upper - 1
            fraction = (value - points[lower]) / (points[upper] - points[lower])
            per_axis.append(((lower, 1.0 - fraction), (upper, fraction)))
        strides = []
        stride = 1
        for points in reversed(self.breakpoints):
            strides.append(stride)
            stride *= len(points)
        strides.reverse()
        corners: list[tuple[int, float]] = [(0, 1.0)]
        for axis_index, pair in enumerate(per_axis):
            corners = [
                (flat + index * strides[axis_index], weight * w)
                for flat, weight in corners
                for index, w in pair
                if w != 0.0
            ]
        return corners

    def query(self, point: Mapping[str, float]) -> AeroQuery:
        """Interpolate every coefficient and its uncertainty at ``point``.

        Parameters
        ----------
        point : mapping of str to float
            One value per axis.

        Returns
        -------
        AeroQuery
            Nominal values, bounds and ``C_Q -+ U``.

        Raises
        ------
        DomainError
            If any coordinate lies outside the grid (no extrapolation), or
            outside every range of a range uncertainty.

        Notes
        -----
        Multilinear interpolation between the ``2^N`` corners of the
        enclosing cell; exact at breakpoints. Gridded uncertainty is
        interpolated the same way (the papers deliver uncertainty as bounds;
        interpolating them is this implementation's choice, flagged with an
        ``info`` warning).

        References
        ----------
        .. [1] SHEA-2020-SLS-ASCENT-DB. P. R. Shea et al., NTRS 20200002382.
        .. [2] FAVAREGH-2016-SLS-UQ. A. L. Favaregh et al., NTRS 20160007692,
               equation 1.
        """

        corners = self._weights(point)
        nominal = {
            name: sum(weight * values[index] for index, weight in corners)
            for name, values in self.coefficients.items()
        }
        uncertainty = {
            name: sum(weight * values[index] for index, weight in corners)
            for name, values in self.gridded_uncertainty.items()
        }
        for name, ranged in self.range_uncertainty.items():
            uncertainty[name] = ranged.bound_at(point[ranged.axis])
        warnings = ()
        if self.gridded_uncertainty and len(corners) > 1:
            warnings = (
                WarningMessage(
                    "interpolated_uncertainty",
                    "Uncertainty bounds were interpolated between breakpoints.",
                    severity="info",
                ),
            )
        return AeroQuery(
            point={name: float(point[name]) for name in self.axes},
            nominal=nominal,
            uncertainty=uncertainty,
            lower={name: nominal[name] - u for name, u in uncertainty.items()},
            upper={name: nominal[name] + u for name, u in uncertainty.items()},
            warnings=warnings,
        )


# ---------------------------------------------------------------------------
# Uncertainty build-up
# ---------------------------------------------------------------------------


def range_method_sigma(ranges: Sequence[float], *, group_size: int) -> float:
    """Estimate a standard deviation from ranges, ``sigma = R_bar / d2``.

    Parameters
    ----------
    ranges : sequence of float
        Ranges (max - min) of repeat groups, all of the same size.
    group_size : int
        Number of measurements per group, 2 to 5.

    Returns
    -------
    float
        ``mean(R) / d2(n)`` with ``d2 = 1.128, 1.693, 2.059, 2.326`` for
        ``n = 2..5`` [1]_ (equation 2); Favaregh et al. apply the same ratio
        to each repeat set (``sigma_rep = R / d2``, equation 3 of [2]_).

    References
    ----------
    .. [1] LEE-2024-ARTEMIS-UQ. M. W. Lee et al., NTRS 20240000054, equation 2.
    .. [2] FAVAREGH-2016-SLS-UQ. A. L. Favaregh et al., NTRS 20160007692,
           equation 3.
    """

    if group_size not in RANGE_METHOD_D2 or isinstance(group_size, bool):
        raise InputError("group_size must be 2, 3, 4 or 5.", field="group_size")
    values = tuple(_finite(v, "ranges") for v in ranges)
    if not values or any(v < 0.0 for v in values):
        raise InputError("ranges must be non-empty and non-negative.", field="ranges")
    return sum(values) / len(values) / RANGE_METHOD_D2[group_size]


def repeatability_uncertainty(sigma: float, *, combined_terms: int = 1) -> float:
    """Return the 3-sigma experimental bound ``3 sqrt(n) sigma_rep``.

    Parameters
    ----------
    sigma : float
        Repeatability standard deviation ``sigma_rep``.
    combined_terms : int, default 1
        Number ``n`` of test-data terms combined (by addition or subtraction)
        in the database quantity; ``n = 2`` for increments formed as a
        difference of two measurements.

    Returns
    -------
    float
        ``sqrt(n (3 sigma)^2) = 3 sqrt(n) sigma`` [1]_ (equations 4 and 5);
        ``n = 1`` gives the ``3 sigma_RM`` of [2]_ (equation 3).

    References
    ----------
    .. [1] FAVAREGH-2016-SLS-UQ. A. L. Favaregh et al., NTRS 20160007692,
           equations 4 and 5.
    .. [2] LEE-2024-ARTEMIS-UQ. M. W. Lee et al., NTRS 20240000054, equation 3.
    """

    value = _finite(sigma, "sigma")
    if value < 0.0:
        raise InputError("sigma must be non-negative.", field="sigma")
    if not isinstance(combined_terms, int) or isinstance(combined_terms, bool) or combined_terms < 1:
        raise InputError("combined_terms must be a positive integer.", field="combined_terms")
    return 3.0 * sqrt(combined_terms) * value


def coverage_bound(residuals: Sequence[float], *, coverage: float = SLS_COVERAGE) -> float:
    """Return the smallest bound that covers ``coverage`` of ``|residuals|``.

    Parameters
    ----------
    residuals : sequence of float
        Residuals or ranges (absolute values are taken).
    coverage : float, default 0.997
        Fraction to cover, ``0 < coverage <= 1``.

    Returns
    -------
    float
        The ``ceil(coverage n)``-th smallest absolute residual.

    Notes
    -----
    The SLS papers define bounds "to cover 99.7%" of residuals or ranges
    [1]_ [2]_ without naming a quantile estimator; the order statistic used
    here guarantees that at least the requested fraction lies at or below the
    bound.

    References
    ----------
    .. [1] SHEA-2020-SLS-ASCENT-DB. P. R. Shea et al., NTRS 20200002382,
           section VII.A-B.
    .. [2] FAVAREGH-2016-SLS-UQ. A. L. Favaregh et al., NTRS 20160007692,
           sections IV-V.
    """

    values = sorted(abs(_finite(v, "residuals")) for v in residuals)
    if not values:
        raise InputError("residuals must not be empty.", field="residuals")
    fraction = _finite(coverage, "coverage")
    if not 0.0 < fraction <= 1.0:
        raise InputError("coverage must lie in (0, 1].", field="coverage")
    rank = max(1, ceil(fraction * len(values) - 1.0e-12))
    return values[rank - 1]


def net_modeling_uncertainty(modeling: float, repeatability: float) -> float:
    """Return ``sqrt(U_DBM^2 - U_rpt^2)``, or zero if repeatability exceeds it.

    Removes the repeatability already counted in the experimental term from a
    modelling bound computed from residuals against test data [1]_ (equations
    2 and 3), [2]_ (equations 7, 8 and 11).

    References
    ----------
    .. [1] SHEA-2020-SLS-ASCENT-DB. P. R. Shea et al., NTRS 20200002382,
           equations 2 and 3.
    .. [2] FAVAREGH-2016-SLS-UQ. A. L. Favaregh et al., NTRS 20160007692,
           equation 7.
    """

    model = _finite(modeling, "modeling")
    rep = _finite(repeatability, "repeatability")
    if model < 0.0 or rep < 0.0:
        raise InputError("Uncertainty bounds must be non-negative.", field="modeling")
    return sqrt(model * model - rep * rep) if model > rep else 0.0


@dataclass(frozen=True, slots=True)
class UncertaintyTerm:
    """One independent uncertainty source in a root-sum-square build-up.

    Attributes
    ----------
    name : str
        Label (e.g. ``"balance"``, ``"repeatability"``, ``"ground_to_flight"``).
    bound : float
        3-sigma-equivalent bound ``u`` (non-negative).
    weight : float, default 1.0
        Multiplier ``w`` of ``u^2`` inside the root-sum-square: 2 for the
        doubled calibration term of the Artemis I full-stack model, or an
        inner uncertainty factor such as ``F_exp`` and ``F_db`` [1]_
        (equations 5-7).

    References
    ----------
    .. [1] LEE-2024-ARTEMIS-UQ. M. W. Lee et al., NTRS 20240000054,
           equations 5-7.
    """

    name: str
    bound: float
    weight: float = 1.0

    def __post_init__(self) -> None:
        if _finite(self.bound, "bound") < 0.0:
            raise InputError("bound must be non-negative.", field="bound")
        if _finite(self.weight, "weight") < 0.0:
            raise InputError("weight must be non-negative.", field="weight")


@dataclass(frozen=True, slots=True)
class UncertaintyBuildUp:
    """Total bound and the contribution of each term.

    Attributes
    ----------
    total : float
        ``UF sqrt(sum(w_i u_i^2))``.
    uncertainty_factor : float
        The outer factor ``UF``.
    contributions : dict of str to float
        ``w_i u_i^2 / sum(w u^2)`` per term (zero when all terms vanish).

    References
    ----------
    .. [1] SHEA-2020-SLS-ASCENT-DB. P. R. Shea et al., NTRS 20200002382,
           equation 4.
    """

    total: float
    uncertainty_factor: float
    contributions: dict[str, float]


def combine_uncertainty_terms(
    terms: Sequence[UncertaintyTerm], *, uncertainty_factor: float = 1.0
) -> UncertaintyBuildUp:
    """Root-sum-square independent terms and apply an uncertainty factor.

    Parameters
    ----------
    terms : sequence of UncertaintyTerm
        Independent, uncorrelated sources.
    uncertainty_factor : float, default 1.0
        Outer factor ``UF >= 1`` for unmodelled effects.

    Returns
    -------
    UncertaintyBuildUp
        ``U = UF sqrt(sum(w_i u_i^2))``.

    Notes
    -----
    Special cases in the cited papers:

    * SLS Block 1B ascent [1]_ (equation 4), ``UF = 1.05``:
      ``U = UF sqrt(U_DBM,net^2 + U_tst,net^2 + U_rpt^2 + U_bal^2 + U_G2F^2
      + U_WTcorr^2)`` with ``U_WTcorr = U_bal``.
    * SLS Block 1 ascent and transition [2]_ (equations 2, 9, 12):
      ``U = UF sqrt(U_DBMnet^2 + U_EXP^2 + U_bal^2)``.
    * Artemis I full stack [3]_ (equation 7):
      ``F sqrt(2 u_cal^2 + u_rep^2 + u_mod^2 + u_G2F^2 + u_T2T^2)``; the
      transition model (equation 5) puts ``F_exp`` and ``F_db`` inside the
      root as weights with ``UF = 1``.

    The papers assume the sources are independent and uncorrelated [2]_; a
    source perfectly correlated with another is handled by weighting (the
    Artemis I calibration term is doubled for that reason [3]_).

    References
    ----------
    .. [1] SHEA-2020-SLS-ASCENT-DB. P. R. Shea et al., NTRS 20200002382,
           equation 4.
    .. [2] FAVAREGH-2016-SLS-UQ. A. L. Favaregh et al., NTRS 20160007692,
           equations 2, 9-12.
    .. [3] LEE-2024-ARTEMIS-UQ. M. W. Lee et al., NTRS 20240000054,
           equations 5-11.
    """

    items = tuple(terms)
    if not items or not all(isinstance(t, UncertaintyTerm) for t in items):
        raise InputError("terms must be a non-empty sequence of UncertaintyTerm.", field="terms")
    names = [t.name for t in items]
    if len(set(names)) != len(names):
        raise InputError("Term names must be unique.", field="terms")
    factor = _finite(uncertainty_factor, "uncertainty_factor")
    if factor < 1.0:
        raise InputError("uncertainty_factor must be at least 1.", field="uncertainty_factor")
    squares = {t.name: t.weight * t.bound * t.bound for t in items}
    total_square = sum(squares.values())
    contributions = {
        name: (value / total_square if total_square > 0.0 else 0.0) for name, value in squares.items()
    }
    return UncertaintyBuildUp(factor * sqrt(total_square), factor, contributions)


# ---------------------------------------------------------------------------
# Slender-body-theory linear fits (Hemsch 2016)
# ---------------------------------------------------------------------------


@dataclass(frozen=True, slots=True)
class LinearFit:
    """Least-squares straight line ``y = intercept + slope x``.

    Attributes
    ----------
    intercept, slope : float
        Fitted coefficients.
    residual_rms : float
        Root-mean-square of the residuals.
    count : int
        Number of points.

    References
    ----------
    .. [1] HEMSCH-2016-SLS-CFD-WT. M. J. Hemsch, NTRS 20160007695,
           equations 1-4.
    """

    intercept: float
    slope: float
    residual_rms: float
    count: int


def fit_linear(x: Sequence[float], y: Sequence[float]) -> LinearFit:
    """Fit ``y = y0 + b x`` by ordinary least squares.

    Hemsch [1]_ models each Mach number separately with ``C_N = C_N0 +
    C_N_alpha alpha`` (equation 1) and ``C_Y = C_Y0 + C_Y_beta beta``
    (equation 2), "least squares is used to choose the coefficients".

    References
    ----------
    .. [1] HEMSCH-2016-SLS-CFD-WT. M. J. Hemsch, NTRS 20160007695,
           equations 1 and 2.
    """

    xs = tuple(_finite(v, "x") for v in x)
    ys = tuple(_finite(v, "y") for v in y)
    if len(xs) < 2 or len(xs) != len(ys):
        raise InputError("Need at least two (x, y) pairs of equal length.", field="x")
    mean_x = sum(xs) / len(xs)
    mean_y = sum(ys) / len(ys)
    sxx = sum((v - mean_x) ** 2 for v in xs)
    if sxx == 0.0:
        raise DomainError("x values are all equal; slope undefined.", field="x")
    slope = sum((a - mean_x) * (b - mean_y) for a, b in zip(xs, ys, strict=True)) / sxx
    intercept = mean_y - slope * mean_x
    rms = sqrt(sum((b - intercept - slope * a) ** 2 for a, b in zip(xs, ys, strict=True)) / len(xs))
    return LinearFit(intercept, slope, rms, len(xs))


@dataclass(frozen=True, slots=True)
class MomentFit:
    """Hemsch's moment-versus-force fit and the implied centre of pressure.

    Attributes
    ----------
    fit : LinearFit
        ``C_m = C_m0 + (x_CN - x_BMC) C_N`` fitted as a line in ``C_N``.
    center_of_pressure : float
        ``x_CN = x_BMC + slope`` in the units of ``x_BMC`` (Hemsch normalises
        by the core diameter).

    References
    ----------
    .. [1] HEMSCH-2016-SLS-CFD-WT. M. J. Hemsch, NTRS 20160007695,
           equations 3 and 4.
    """

    fit: LinearFit
    center_of_pressure: float


def fit_moment_center_of_pressure(
    moment_coefficients: Sequence[float],
    force_coefficients: Sequence[float],
    *,
    moment_reference: float,
) -> MomentFit:
    """Fit ``C_m = C_m0 + (x_CN - x_BMC) C_N`` and return ``x_CN``.

    Parameters
    ----------
    moment_coefficients, force_coefficients : sequence of float
        Paired ``C_m``/``C_N`` (or ``C_n``/``C_Y``) values at one Mach number.
    moment_reference : float
        ``x_BMC``, the normalised moment-reference location.

    Returns
    -------
    MomentFit
        Fit and centre of pressure.

    Notes
    -----
    Slender-body theory predicts a centre of pressure independent of
    ``alpha``, ``beta`` and ``M``, hence a moment linear in the force [1]_
    (equations 3 and 4, as printed). The excerpt does not state the direction
    of the ``x`` axis; the sign of ``x_CN - x_BMC`` follows the user's
    moment sign convention.

    References
    ----------
    .. [1] HEMSCH-2016-SLS-CFD-WT. M. J. Hemsch, NTRS 20160007695,
           equations 3 and 4.
    """

    reference = _finite(moment_reference, "moment_reference")
    line = fit_linear(force_coefficients, moment_coefficients)
    return MomentFit(line, reference + line.slope)
