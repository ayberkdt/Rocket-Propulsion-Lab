"""Throttle-commanded chemical engine model for trajectory optimization.

The propagator bridge replays a *prescribed* thrust history.  A trajectory
optimizer or a closed-loop guidance law instead chooses the throttle ``u`` as
a control and needs ``a(u)``, ``mdot(u)`` and their derivatives with respect to
``u`` [4]_.  :class:`ThrottleableEngineModel` provides that::

    F(u)    = F_rated * f_T(u) * s_F - p_a * A_e
    mdot(u) = mdot_rated * f_m(u) * s_m

with ``u = 0`` (off) or ``u`` in ``[u_min, u_max]``.  The ``- p_a * A_e`` term
is the ambient-pressure correction of the rocket thrust equation [1]_.  The
thrust and flow fractions ``f_T`` and ``f_m`` default to ``u``; measured
throttle characteristics (for example specific impulse falling at deep
throttle) are given as :class:`MonotoneCubicCurve` instances.  These monotone
piecewise cubic Hermite interpolants [2]_, [3]_ are C1-continuous, so a
gradient-based optimizer sees no kinks, and they never overshoot monotone
data.

Throttle values strictly between zero and ``u_min`` are refused: a
deep-throttle gap is a hard engine limit, not something to interpolate.

References
----------
.. [1] NASA Glenn Research Center, "Rocket Thrust Equation".
       https://www1.grc.nasa.gov/beginners-guide-to-aeronautics/rocket-thrust-equation/
.. [2] F. N. Fritsch and R. E. Carlson, "Monotone Piecewise Cubic
       Interpolation", SIAM J. Numer. Anal. 17(2), 1980.
.. [3] F. N. Fritsch and J. Butland, "A Method for Constructing Local
       Monotone Piecewise Cubic Interpolants", SIAM J. Sci. Stat. Comput.
       5(2), 1984.
.. [4] Orekit 13.1.5, "PropulsionModel / EventDetector API", API documentation of
       org.orekit.forces.maneuvers.propulsion.PropulsionModel.
"""

from __future__ import annotations

import json
from bisect import bisect_right
from collections.abc import Mapping
from dataclasses import dataclass, field
from fractions import Fraction
from itertools import pairwise
from typing import cast

from rocket_propulsion.core.errors import DomainError, InputError
from rocket_propulsion.core.validation import (
    name_value_pairs,
    real,
    real_mapping,
    real_tuple,
    records,
    text,
)
from rocket_propulsion.propulsion.performance import STANDARD_GRAVITY_M_S2

from .direction import _direction_geometry
from .exchange_documents import _metadata, _read_document
from .interpolation import _hermite_coefficients, _hermite_segment, _pchip_rational_slopes
from .mass_closure import _MASS_CLOSURE_POLICY, _inventory_ledger
from .numeric_outputs import _finite_rational, _fraction
from .power_documents import _throttle_document
from .propagator import (
    AdditionalStateConstraint,
    BurnEventKind,
    BurnEventSurface,
    EstimableParameter,
    PropellantMassClosure,
    _document_parameters,
)
from .serialization import _require_numerical_policies, canonical_sha256

THROTTLEABLE_ENGINE_SCHEMA = "rocket_propulsion_throttleable_engine_v1"
THROTTLEABLE_REFERENCE_IDS = ("NASA-20040000363", "OREKIT-13.1.5-PROPULSION-MODEL")
_THROTTLE_TOLERANCE = 1e-12






@dataclass(frozen=True, slots=True)
class MonotoneCubicCurve:
    """Shape-preserving piecewise cubic Hermite curve (PCHIP) without extrapolation.

    Parameters
    ----------
    x : tuple of float
        At least two strictly increasing finite abscissae.
    y : tuple of float
        Finite ordinates, one per abscissa.

    Raises
    ------
    DomainError
        If fewer than two nodes are given, the lengths differ, a value is
        non-finite, or the abscissae are not strictly increasing.

    Notes
    -----
    Nodal slopes follow Fritsch-Carlson [1]_ with the Fritsch-Butland
    weighted-harmonic-mean interior slopes [2]_, the scheme also used by
    [3]_.  The curve passes through every node, is continuously
    differentiable, and is monotone on every interval where the data are
    monotone (no overshoot).  Slopes and segment coefficients are held as
    exact rationals; evaluation rounds the value and the analytic derivative
    once at output.

    References
    ----------
    .. [1] F. N. Fritsch and R. E. Carlson, "Monotone Piecewise Cubic
           Interpolation", SIAM J. Numer. Anal. 17(2), 1980.
    .. [2] F. N. Fritsch and J. Butland, "A Method for Constructing Local
           Monotone Piecewise Cubic Interpolants", SIAM J. Sci. Stat. Comput.
           5(2), 1984.
    .. [3] SCIPY-1.14.1-PCHIP. SciPy 1.14.1, documentation and implementation of
           "scipy.interpolate.PchipInterpolator".
    """

    x: tuple[float, ...]
    y: tuple[float, ...]
    _slopes: tuple[Fraction, ...] = field(init=False, repr=False, compare=False)
    _segments: tuple[tuple[int, ...], ...] = field(init=False, repr=False, compare=False)

    def __post_init__(self) -> None:
        object.__setattr__(self, "x", real_tuple(self.x, "Curve abscissae"))
        object.__setattr__(self, "y", real_tuple(self.y, "Curve ordinates"))
        if len(self.x) < 2 or len(self.x) != len(self.y):
            raise DomainError("A curve needs at least two nodes and matching x/y lengths.")
        if any(b <= a for a, b in pairwise(self.x)):
            raise DomainError("Curve abscissae must be strictly increasing.")
        object.__setattr__(self, "_slopes", _pchip_rational_slopes(self.x, self.y))
        object.__setattr__(self, "_segments", tuple(
            _hermite_coefficients(self.x[i], self.x[i+1], self.y[i], self.y[i+1],
                                  self._slopes[i], self._slopes[i+1])
            for i in range(len(self.x)-1)))

    def evaluate(self, value: float) -> tuple[float, float]:
        """Return the curve value and its derivative at ``value``.

        Parameters
        ----------
        value : float
            Abscissa inside ``[x[0], x[-1]]``.  Values within a relative
            ``1e-12`` (absolute ``1e-12`` near zero) of an end node are
            clamped to it.

        Returns
        -------
        tuple of float
            ``(y(value), dy/dvalue)`` of the monotone cubic Hermite
            interpolant [1]_.

        Raises
        ------
        DomainError
            If ``value`` is non-finite or outside the node range beyond the
            clamping tolerance; extrapolation is refused.

        References
        ----------
        .. [1] F. N. Fritsch and R. E. Carlson, "Monotone Piecewise Cubic
               Interpolation", SIAM J. Numer. Anal. 17(2), 1980.
        """

        value = real(value, "Curve argument")
        tolerance = 1e-12 * max(1.0, abs(self.x[0]), abs(self.x[-1]))
        if value < self.x[0] - tolerance or value > self.x[-1] + tolerance:
            raise DomainError(
                f"Curve argument {value} lies outside [{self.x[0]}, {self.x[-1]}]; "
                "extrapolation refused."
            )
        return self._evaluate_checked(min(max(value, self.x[0]), self.x[-1]))

    def _evaluate_checked(self, value: float) -> tuple[float, float]:
        """Return ``(value, derivative)`` for an argument already validated and clamped."""

        index = min(max(bisect_right(self.x, value) - 1, 0), len(self.x) - 2)
        return _hermite_segment(self._segments[index], value)

    def _value_checked(self, value: float) -> float:
        """Return the value only for an argument already validated and clamped."""
        index = min(max(bisect_right(self.x, value)-1, 0), len(self.x)-2)
        return _hermite_segment(self._segments[index], value, derivative=False)

    def to_dict(self) -> dict[str, object]:
        """Return the JSON-ready node document ``{"x": [...], "y": [...]}``."""
        return {"x": list(self.x), "y": list(self.y)}


@dataclass(frozen=True, slots=True)
class ThrottleablePartials:
    """Analytic local derivatives of one throttle-commanded evaluation.

    Attributes
    ----------
    acceleration_wrt_mass : tuple of float
        ``da/dm = -F * d_hat / m**2``, m/s^2 per kg.
    acceleration_wrt_direction : tuple of tuple of float
        3 x 3 Jacobian of ``a = F * d / (m * |d|)`` with respect to the
        unnormalized direction ``d``.
    acceleration_wrt_throttle : tuple of float
        ``(dF/du) * d_hat / m``.
    mass_rate_wrt_throttle : float
        ``-dmdot/du``, kg/s.
    additional_state_rates_wrt_throttle : tuple of (str, float)
        Per-tank share of ``-dmdot/du``.
    acceleration_wrt_ambient_pressure : tuple of float
        ``-A_e * d_hat / m``; zero without a nozzle exit area.
    acceleration_wrt_parameters : tuple of (str, tuple of float)
        Acceleration sensitivities to ``thrust_scale`` and
        ``mass_flow_scale``.
    mass_rate_wrt_parameters : tuple of (str, float)
        Vehicle mass-rate sensitivities to the scale parameters.
    additional_state_rates_wrt_parameters : tuple of (str, tuple of (str, float))
        Per-tank mass-rate sensitivities to the scale parameters.

    Notes
    -----
    The blocks are the acceleration, mass-derivative and parameter-driver
    partials of a propagation-time propulsion model [1]_; all are zero while
    the engine is off or inhibited.

    References
    ----------
    .. [1] Orekit 13.1.5, "PropulsionModel / EventDetector API", API documentation of
           org.orekit.forces.maneuvers.propulsion.PropulsionModel.
    """

    acceleration_wrt_mass: tuple[float, float, float]
    acceleration_wrt_direction: tuple[
        tuple[float, float, float],
        tuple[float, float, float],
        tuple[float, float, float],
    ]
    acceleration_wrt_throttle: tuple[float, float, float]
    mass_rate_wrt_throttle: float
    additional_state_rates_wrt_throttle: tuple[tuple[str, float], ...]
    acceleration_wrt_ambient_pressure: tuple[float, float, float]
    acceleration_wrt_parameters: tuple[tuple[str, tuple[float, float, float]], ...]
    mass_rate_wrt_parameters: tuple[tuple[str, float], ...]
    additional_state_rates_wrt_parameters: tuple[tuple[str, tuple[tuple[str, float], ...]], ...]


@dataclass(frozen=True, slots=True)
class ThrottleableEvaluation:
    """Force, tank drains and partials for one throttle command.

    Attributes
    ----------
    throttle : float
        Validated throttle; ``0.0`` when off.
    active : bool
        Whether the engine fires.
    thrust_n : float
        Axial thrust after the ambient-pressure correction, N.
    force_vector_n, acceleration_m_s2 : tuple of float
        Thrust force, N, and thrust acceleration, m/s^2.
    mass_derivative_kg_s : float
        Vehicle mass rate, ``-mdot``, kg/s.
    additional_state_derivatives : tuple of (str, float)
        Tank mass rates, kg/s.
    vacuum_specific_impulse_s : float or None
        ``F_vac / (g0 * mdot)`` with ``g0`` the registered
        ``STANDARD_GRAVITY_M_S2``; ``None`` when off.
    net_thrust_margin_n : float or None
        Net thrust with ambient back-pressure, N.
    parameters : tuple of (str, float)
        Resolved parameter values.
    inhibited_by : tuple of str
        Inventory limits that inhibited firing.
    partials : ThrottleablePartials
        Analytic local derivatives.

    Notes
    -----
    ``net_thrust_margin_n`` is defined only while the engine fires and has a
    nozzle exit area; otherwise it is ``None``, so a consumer monitors it only
    while firing.  Thrust and its ambient-pressure correction follow [1]_ and
    specific impulse follows [2]_; [3]_ gives the transient-modelling context.

    References
    ----------
    .. [1] NASA Glenn Research Center, "Rocket Thrust Equation".
           https://www1.grc.nasa.gov/beginners-guide-to-aeronautics/rocket-thrust-equation/
    .. [2] NASA Glenn Research Center, "Specific Impulse".
           https://www1.grc.nasa.gov/beginners-guide-to-aeronautics/specific-impulse/
    .. [3] NASA, "Transient Mathematical Modeling for Liquid Rocket Engine
           Systems". https://ntrs.nasa.gov/citations/20040000363
    """

    throttle: float
    active: bool
    thrust_n: float
    force_vector_n: tuple[float, float, float]
    acceleration_m_s2: tuple[float, float, float]
    mass_derivative_kg_s: float
    additional_state_derivatives: tuple[tuple[str, float], ...]
    vacuum_specific_impulse_s: float | None
    net_thrust_margin_n: float | None
    parameters: tuple[tuple[str, float], ...]
    inhibited_by: tuple[str, ...]
    partials: ThrottleablePartials


_ZERO3 = (0.0, 0.0, 0.0)


@dataclass(frozen=True, slots=True)
class ThrottleableEngineModel:
    """Throttle-commanded engine with smooth throttle characteristics.

    Parameters
    ----------
    rated_thrust_n : float
        Positive rated vacuum thrust ``F_rated``, N.
    rated_mass_flow_kg_s : float
        Positive rated total mass flow ``mdot_rated``, kg/s.
    minimum_throttle, maximum_throttle : float
        Throttle limits with ``0 < minimum <= maximum``.
    stream_fractions : tuple of (str, float)
        Non-negative split of total flow into uniquely named tank states, for
        example ``(("m_ox", 0.78), ("m_fu", 0.22))``.
    thrust_curve, flow_curve : MonotoneCubicCurve or None
        Positive fractions of rated thrust and flow as functions of throttle;
        ``None`` uses the fraction ``u``.
    state_constraints : tuple of AdditionalStateConstraint
        Tank reserves; each must refer to a drained state.
    protected_dry_mass_kg : float or None
        Positive vehicle mass below which firing is inhibited, kg.
    nozzle_exit_area_m2 : float or None
        Positive nozzle exit area ``A_e`` for the ambient correction, m^2.
    parameters : tuple of EstimableParameter
        Exactly ``thrust_scale`` and ``mass_flow_scale``; empty selects both
        with reference value 1.
    model_id : str, default "throttleable-engine"
        Model identifier.
    schema : str
        Must equal ``THROTTLEABLE_ENGINE_SCHEMA``.

    Raises
    ------
    DomainError
        If the schema is unsupported or any field violates the constraints
        above.

    Notes
    -----
    Thrust is ``F = F_rated * f_T(u) * s_F - p_a * A_e`` [1]_ and flow
    ``mdot = mdot_rated * f_m(u) * s_m``.  Stream fractions must sum to one
    within ``1e-12``; accepted fractions are normalized by their exact binary
    sum, which preserves relative weights.  The internal total/tank ledger
    closes exactly; exported fields round independently.  ``thrust_curve``
    and ``flow_curve`` must cover ``[minimum_throttle, maximum_throttle]``.
    The model exposes the acceleration, mass-derivative, parameter-driver and
    event contract of [2]_; [3]_ gives the transient-modelling context.

    References
    ----------
    .. [1] NASA Glenn Research Center, "Rocket Thrust Equation".
           https://www1.grc.nasa.gov/beginners-guide-to-aeronautics/rocket-thrust-equation/
    .. [2] Orekit 13.1.5, "PropulsionModel / EventDetector API", API documentation of
           org.orekit.forces.maneuvers.propulsion.PropulsionModel.
    .. [3] NASA, "Transient Mathematical Modeling for Liquid Rocket Engine
           Systems". https://ntrs.nasa.gov/citations/20040000363
    """

    rated_thrust_n: float
    rated_mass_flow_kg_s: float
    minimum_throttle: float
    maximum_throttle: float
    stream_fractions: tuple[tuple[str, float], ...]
    thrust_curve: MonotoneCubicCurve | None = None
    flow_curve: MonotoneCubicCurve | None = None
    state_constraints: tuple[AdditionalStateConstraint, ...] = ()
    protected_dry_mass_kg: float | None = None
    nozzle_exit_area_m2: float | None = None
    parameters: tuple[EstimableParameter, ...] = ()
    model_id: str = "throttleable-engine"
    schema: str = THROTTLEABLE_ENGINE_SCHEMA

    def __post_init__(self) -> None:
        object.__setattr__(self, "schema", text(self.schema, "Throttleable model schema"))
        if self.schema != THROTTLEABLE_ENGINE_SCHEMA:
            raise DomainError("Unsupported throttleable engine schema.")
        object.__setattr__(self, "model_id", text(self.model_id, "Model identifier"))
        for name in ("rated_thrust_n", "rated_mass_flow_kg_s", "minimum_throttle",
                     "maximum_throttle"):
            object.__setattr__(self, name, real(getattr(self, name), name))
        if self.rated_thrust_n <= 0.0 or self.rated_mass_flow_kg_s <= 0.0:
            raise DomainError("Rated thrust and mass flow must be greater than zero.")
        if not 0.0 < self.minimum_throttle <= self.maximum_throttle:
            raise DomainError("Throttle limits must satisfy 0 < minimum <= maximum.")
        fractions = list(name_value_pairs(self.stream_fractions, "Stream fractions"))
        if any(fraction < 0.0 for _, fraction in fractions):
            raise DomainError("Stream fractions must be non-negative.")
        if not fractions or len({name for name, _ in fractions}) != len(fractions):
            raise DomainError("Stream fractions need unique, non-empty state names.")
        fraction_sum = sum((_fraction(value) for _, value in fractions), Fraction(0))
        if abs(fraction_sum - 1) > _fraction(1e-12):
            raise DomainError("Stream fractions must sum to one.")
        object.__setattr__(self, "stream_fractions", tuple(fractions))
        for curve_name in ("thrust_curve", "flow_curve"):
            curve = getattr(self, curve_name)
            if curve is None:
                continue
            if not isinstance(curve, MonotoneCubicCurve):
                raise DomainError(f"{curve_name} must be a MonotoneCubicCurve.")
            if curve.x[0] > self.minimum_throttle + 1e-12 or (
                curve.x[-1] < self.maximum_throttle - 1e-12
            ):
                raise DomainError(f"{curve_name} must cover the full throttle range.")
            if any(value <= 0.0 for value in curve.y):
                raise DomainError(f"{curve_name} fractions must be greater than zero.")
        for name in ("protected_dry_mass_kg", "nozzle_exit_area_m2"):
            value = getattr(self, name)
            if value is not None:
                value = real(value, name)
                if value <= 0.0:
                    raise DomainError(f"{name} must be greater than zero.")
                object.__setattr__(self, name, value)
        for name, kind in (("state_constraints", AdditionalStateConstraint),
                           ("parameters", EstimableParameter)):
            object.__setattr__(self, name, records(getattr(self, name), kind, name))
        states = {name for name, _ in self.stream_fractions}
        if any(constraint.state_name not in states for constraint in self.state_constraints):
            raise DomainError("State constraints must refer to a drained tank state.")
        if not self.parameters:
            object.__setattr__(
                self,
                "parameters",
                (
                    EstimableParameter("thrust_scale", 1.0, 1e-3, 0.01, 2.0, "1"),
                    EstimableParameter("mass_flow_scale", 1.0, 1e-3, 0.01, 2.0, "1"),
                ),
            )
        if sorted(parameter.name for parameter in self.parameters) != [
            "mass_flow_scale",
            "thrust_scale",
        ]:
            raise DomainError("Throttleable model requires exactly thrust_scale and mass_flow_scale.")

    def _fraction(self, curve: MonotoneCubicCurve | None, throttle: float) -> tuple[float, float]:
        """Return ``(fraction, d fraction / du)``; ``(u, 1)`` when no curve is set."""
        if curve is None:
            return throttle, 1.0
        return curve.evaluate(throttle)

    def validated_throttle(self, throttle: float) -> float:
        """Return a throttle that is zero or inside the limits.

        Parameters
        ----------
        throttle : float
            Requested throttle.

        Returns
        -------
        float
            ``0.0`` when ``|throttle| <= 1e-12``; otherwise the throttle
            clamped to ``[minimum_throttle, maximum_throttle]`` when it lies
            within ``1e-12`` of that interval.

        Raises
        ------
        DomainError
            If the throttle is non-finite or lies in neither region, including
            the deep-throttle gap between zero and the minimum.
        """

        value = real(throttle, "Throttle")
        if abs(value) <= _THROTTLE_TOLERANCE:
            return 0.0
        if self.minimum_throttle - _THROTTLE_TOLERANCE <= value <= (
            self.maximum_throttle + _THROTTLE_TOLERANCE
        ):
            return min(max(value, self.minimum_throttle), self.maximum_throttle)
        raise DomainError(
            f"Throttle {value} must be 0 or within [{self.minimum_throttle}, "
            f"{self.maximum_throttle}]."
        )

    def inventory_surfaces(self) -> tuple[BurnEventSurface, ...]:
        """Return the burn-stopping dry-mass and tank-reserve event surfaces.

        Returns
        -------
        tuple of BurnEventSurface
            A dry-mass surface when ``protected_dry_mass_kg`` is set, followed
            by one reserve surface per state constraint, for an event
            detector [1]_.

        References
        ----------
        .. [1] Orekit 13.1.5, "PropulsionModel / EventDetector API", API documentation of
               org.orekit.forces.maneuvers.propulsion.PropulsionModel.
        """

        surfaces = []
        if self.protected_dry_mass_kg is not None:
            surfaces.append(
                BurnEventSurface(
                    "engine:dry-mass", BurnEventKind.DRY_MASS, True,
                    threshold=self.protected_dry_mass_kg,
                )
            )
        surfaces.extend(
            BurnEventSurface(
                f"engine:reserve:{item.state_name}", BurnEventKind.TANK_RESERVE, True,
                state_name=item.state_name, threshold=item.minimum_value,
            )
            for item in self.state_constraints
        )
        return tuple(surfaces)

    def resolved_parameters(self, overrides: Mapping[str, float] | None = None) -> dict[str, float]:
        """Resolve named overrides against the bounded parameter drivers.

        Parameters
        ----------
        overrides : Mapping[str, float] or None
            Values for ``thrust_scale`` or ``mass_flow_scale``.

        Returns
        -------
        dict of str to float
            Every parameter, using its reference value when not overridden.

        Raises
        ------
        DomainError
            If a name is unknown or a value violates its driver bounds.

        Notes
        -----
        The drivers are estimable propulsion parameters in the sense of [1]_.

        References
        ----------
        .. [1] Orekit 13.1.5, "PropulsionModel / EventDetector API", API documentation of
               org.orekit.forces.maneuvers.propulsion.PropulsionModel.
        """

        supplied = real_mapping(overrides, "Parameter overrides")
        unknown = set(supplied) - {parameter.name for parameter in self.parameters}
        if unknown:
            raise DomainError("Unknown propulsion parameters: " + ", ".join(sorted(unknown)))
        return {
            parameter.name: parameter.validate(supplied.get(parameter.name, parameter.reference_value))
            for parameter in self.parameters
        }

    def evaluate(
        self,
        *,
        throttle: float,
        vehicle_mass_kg: float,
        direction: tuple[float, float, float],
        additional_states: Mapping[str, float],
        ambient_pressure_pa: float | None = None,
        parameter_overrides: Mapping[str, float] | None = None,
        apply_inventory_limits: bool = True,
    ) -> ThrottleableEvaluation:
        """Evaluate force, tank drains and analytic partials for one throttle command.

        Parameters
        ----------
        throttle : float
            Throttle command, zero or inside the limits.
        vehicle_mass_kg : float
            Positive vehicle mass ``m``, kg.
        direction : tuple of float
            Non-zero thrust direction ``d``; it need not be normalized.
        additional_states : Mapping[str, float]
            Tank states; every drained state must be present.
        ambient_pressure_pa : float or None
            Non-negative ambient pressure ``p_a``, Pa; required while firing
            when a nozzle exit area is set and refused when none is set.
        parameter_overrides : Mapping[str, float] or None
            Overrides of ``thrust_scale`` and ``mass_flow_scale``.
        apply_inventory_limits : bool, default True
            Whether dry-mass and reserve limits (with ``1e-12`` margin)
            inhibit firing.

        Returns
        -------
        ThrottleableEvaluation
            Force, acceleration, drains, vacuum specific impulse and partials;
            all zero when the throttle is zero or a limit inhibits firing.

        Raises
        ------
        DomainError
            If an input is outside its domain, a tank state is missing, the
            ambient pressure is inconsistent with the nozzle exit area, or the
            active flow fraction is not positive.

        Notes
        -----
        Thrust is ``F = F_rated * f_T(u) * s_F - p_a * A_e`` [1]_, flow is
        ``mdot = mdot_rated * f_m(u) * s_m`` split by the normalized stream
        fractions, acceleration is ``a = F * d / (m * |d|)``, and vacuum
        specific impulse is ``F_vac / (g0 * mdot)`` [2]_ with ``g0`` the
        registered ``STANDARD_GRAVITY_M_S2``.  Intermediate quantities are
        exact rationals rounded once per exported field.  The outputs follow
        the propagation-time propulsion contract of [3]_; [4]_ gives the
        transient-modelling context.

        References
        ----------
        .. [1] NASA Glenn Research Center, "Rocket Thrust Equation".
               https://www1.grc.nasa.gov/beginners-guide-to-aeronautics/rocket-thrust-equation/
        .. [2] NASA Glenn Research Center, "Specific Impulse".
               https://www1.grc.nasa.gov/beginners-guide-to-aeronautics/specific-impulse/
        .. [3] Orekit 13.1.5, "PropulsionModel / EventDetector API", API documentation of
               org.orekit.forces.maneuvers.propulsion.PropulsionModel.
        .. [4] NASA, "Transient Mathematical Modeling for Liquid Rocket Engine
               Systems". https://ntrs.nasa.gov/citations/20040000363
        """

        throttle = self.validated_throttle(throttle)
        vehicle_mass_kg = real(vehicle_mass_kg, "Vehicle mass")
        if vehicle_mass_kg <= 0.0:
            raise DomainError("Vehicle mass must be greater than zero.")
        geometry = _direction_geometry(direction)
        states = real_mapping(additional_states, "Additional states")
        names = [name for name, _ in self.stream_fractions]
        missing = set(names) - set(states)
        if missing:
            raise DomainError("Missing tank states: " + ", ".join(sorted(missing)))
        if not isinstance(apply_inventory_limits, bool):
            raise DomainError("apply_inventory_limits must be a boolean.")
        if ambient_pressure_pa is not None:
            if self.nozzle_exit_area_m2 is None:
                raise DomainError("Ambient pressure was supplied but no nozzle exit area is set.")
            ambient_pressure_pa = real(ambient_pressure_pa, "Ambient pressure")
            if ambient_pressure_pa < 0.0:
                raise DomainError("Ambient pressure must be non-negative.")
        parameters = self.resolved_parameters(parameter_overrides)
        thrust_scale = parameters["thrust_scale"]
        flow_scale = parameters["mass_flow_scale"]

        inhibited: list[str] = []
        if apply_inventory_limits:
            if self.protected_dry_mass_kg is not None and (
                vehicle_mass_kg <= self.protected_dry_mass_kg + 1e-12
            ):
                inhibited.append("protected_dry_mass_kg")
            for constraint in self.state_constraints:
                if states[constraint.state_name] <= constraint.minimum_value + 1e-12:
                    inhibited.append(constraint.state_name)

        if throttle == 0.0 or inhibited:
            zero_states = tuple((name, 0.0) for name in names)
            return ThrottleableEvaluation(
                throttle=throttle, active=False, thrust_n=0.0, force_vector_n=_ZERO3,
                acceleration_m_s2=_ZERO3, mass_derivative_kg_s=0.0,
                additional_state_derivatives=zero_states, vacuum_specific_impulse_s=None,
                net_thrust_margin_n=None, parameters=tuple(parameters.items()),
                inhibited_by=tuple(inhibited),
                partials=ThrottleablePartials(
                    acceleration_wrt_mass=_ZERO3,
                    acceleration_wrt_direction=(_ZERO3, _ZERO3, _ZERO3),
                    acceleration_wrt_throttle=_ZERO3,
                    mass_rate_wrt_throttle=0.0,
                    additional_state_rates_wrt_throttle=zero_states,
                    acceleration_wrt_ambient_pressure=_ZERO3,
                    acceleration_wrt_parameters=tuple((n, _ZERO3) for n in parameters),
                    mass_rate_wrt_parameters=tuple((n, 0.0) for n in parameters),
                    additional_state_rates_wrt_parameters=tuple(
                        (n, zero_states) for n in parameters
                    ),
                ),
            )

        if self.nozzle_exit_area_m2 is not None and ambient_pressure_pa is None:
            raise DomainError("Ambient pressure is required while the engine is firing.")
        thrust_fraction, thrust_slope = self._fraction(self.thrust_curve, throttle)
        flow_fraction, flow_slope = self._fraction(self.flow_curve, throttle)
        vacuum = _fraction(self.rated_thrust_n) * _fraction(thrust_fraction) * _fraction(thrust_scale)
        pressure_loss = (
            0.0 if self.nozzle_exit_area_m2 is None
            else _fraction(ambient_pressure_pa) * _fraction(self.nozzle_exit_area_m2)
        )
        exact_thrust = vacuum - _fraction(pressure_loss)
        thrust = _finite_rational(exact_thrust, "Axial thrust")
        base_flow = _fraction(self.rated_mass_flow_kg_s) * _fraction(flow_fraction)
        if base_flow <= 0:
            raise DomainError("Active throttle flow fraction must be positive.",
                              field="Throttle flow fraction")
        flow = base_flow * _fraction(flow_scale)
        fraction_sum = sum((_fraction(value) for _, value in self.stream_fractions), Fraction(0))
        tank_weights = tuple((name, _fraction(value) / fraction_sum)
                             for name, value in self.stream_fractions)
        force = geometry.project(exact_thrust, label="Force")
        acceleration = geometry.project(exact_thrust, vehicle_mass_kg, label="Acceleration")
        d_thrust_du = _fraction(self.rated_thrust_n) * _fraction(thrust_slope) * _fraction(thrust_scale)
        d_flow_du = _fraction(self.rated_mass_flow_kg_s) * _fraction(flow_slope) * _fraction(flow_scale)
        return ThrottleableEvaluation(
            throttle=throttle,
            active=True,
            thrust_n=thrust,
            force_vector_n=force,
            acceleration_m_s2=acceleration,
            mass_derivative_kg_s=_finite_rational(-flow, "Total mass flow"),
            additional_state_derivatives=tuple(
                (name, _finite_rational(-flow * weight, f"Tank mass flow: {name}"))
                for name, weight in tank_weights
            ),
            vacuum_specific_impulse_s=_finite_rational(
                vacuum / (_fraction(STANDARD_GRAVITY_M_S2) * _fraction(flow)), "Vacuum specific impulse"),
            net_thrust_margin_n=None if self.nozzle_exit_area_m2 is None else thrust,
            parameters=tuple(parameters.items()),
            inhibited_by=(),
            partials=ThrottleablePartials(
                acceleration_wrt_mass=geometry.project(
                    -exact_thrust, _fraction(vehicle_mass_kg)**2, label="Mass acceleration partial"),
                acceleration_wrt_direction=geometry.jacobian(exact_thrust, vehicle_mass_kg),
                acceleration_wrt_throttle=geometry.project(
                    d_thrust_du, vehicle_mass_kg, label="Throttle acceleration partial"),
                mass_rate_wrt_throttle=_finite_rational(-d_flow_du, "Throttle mass-flow partial"),
                additional_state_rates_wrt_throttle=tuple(
                    (name, _finite_rational(-d_flow_du * weight,
                                            f"Tank throttle mass-flow partial: {name}"))
                    for name, weight in tank_weights
                ),
                acceleration_wrt_ambient_pressure=(
                    _ZERO3
                    if self.nozzle_exit_area_m2 is None
                    else geometry.project(-self.nozzle_exit_area_m2, vehicle_mass_kg,
                                          label="Ambient pressure acceleration partial")
                ),
                acceleration_wrt_parameters=(
                    (
                        "thrust_scale",
                        geometry.project(_fraction(self.rated_thrust_n) * _fraction(thrust_fraction),
                                         vehicle_mass_kg, label="Thrust scale acceleration partial"),
                    ),
                    ("mass_flow_scale", _ZERO3),
                ),
                mass_rate_wrt_parameters=(
                    ("thrust_scale", 0.0),
                    ("mass_flow_scale", _finite_rational(-base_flow, "Mass-flow scale partial"))),
                additional_state_rates_wrt_parameters=(
                    ("thrust_scale", tuple((name, 0.0) for name in names)),
                    (
                        "mass_flow_scale",
                        tuple(
                            (name, _finite_rational(-base_flow * weight,
                                                    f"Tank mass-flow scale partial: {name}"))
                            for name, weight in tank_weights
                        ),
                    ),
                ),
            ),
        )

    def mass_closure(
        self,
        *,
        vehicle_mass_kg: float,
        additional_states: Mapping[str, float],
        non_propellant_mass_kg: float,
    ) -> PropellantMassClosure:
        """Compare vehicle mass with the drained tank states.

        Parameters
        ----------
        vehicle_mass_kg : float
            Positive vehicle mass, kg.
        additional_states : Mapping[str, float]
            Tank states; every drained state must be present.
        non_propellant_mass_kg : float
            Non-negative mass not held in the drained tanks, kg.

        Returns
        -------
        PropellantMassClosure
            Vehicle, non-propellant and tank-inventory ledger with its
            closure residual.

        Raises
        ------
        DomainError
            If a mass is outside its domain or a tank state is missing.

        Notes
        -----
        Vehicle mass should equal non-propellant mass plus the propellant held
        in the tanks of the pressurization and feed system [1]_.

        References
        ----------
        .. [1] NASA SP-8112. NASA, "Pressurization Systems for Liquid Rockets", 1975.
               https://ntrs.nasa.gov/citations/19760015212
        """

        vehicle_mass_kg = real(vehicle_mass_kg, "Vehicle mass")
        non_propellant_mass_kg = real(non_propellant_mass_kg, "Non-propellant mass")
        if vehicle_mass_kg <= 0.0 or non_propellant_mass_kg < 0.0:
            raise DomainError("Masses must be positive (vehicle) and non-negative (other).")
        states = real_mapping(additional_states, "Additional states")
        names = [name for name, _ in self.stream_fractions]
        missing = set(names) - set(states)
        if missing:
            raise DomainError("Missing tank states: " + ", ".join(sorted(missing)))
        return PropellantMassClosure(
            vehicle_mass_kg, non_propellant_mass_kg,
            *_inventory_ledger(vehicle_mass_kg, non_propellant_mass_kg, (states[name] for name in names)),
        )

    def manifest(self, *, include_hash: bool = True) -> dict[str, object]:
        """Return the complete deterministic model document.

        Parameters
        ----------
        include_hash : bool, default True
            Whether to add ``model_hash``, the canonical SHA-256 of the
            document without it.

        Returns
        -------
        dict
            Schema, numerics policies, rated values, limits, stream fractions,
            curves, constraints, parameters and reference identifiers.
        """

        payload: dict[str, object] = {
            "schema": self.schema,
            "direction_numerics": "scaled_hypot_rational_projection_and_signed_tangents_v1",
            "mass_flow_numerics": "exact_binary_flow_ledger_normalized_stream_weights_v1",
            "mass_closure_numerics": _MASS_CLOSURE_POLICY,
            "curve_numerics": "rational_pchip_and_hermite_export_v1",
            "model_id": self.model_id,
            "rated_thrust_n": self.rated_thrust_n,
            "rated_mass_flow_kg_s": self.rated_mass_flow_kg_s,
            "minimum_throttle": self.minimum_throttle,
            "maximum_throttle": self.maximum_throttle,
            "stream_fractions": [list(item) for item in self.stream_fractions],
            "thrust_curve": None if self.thrust_curve is None else self.thrust_curve.to_dict(),
            "flow_curve": None if self.flow_curve is None else self.flow_curve.to_dict(),
            "state_constraints": [
                {"state_name": c.state_name, "minimum_value": c.minimum_value, "label": c.label}
                for c in self.state_constraints
            ],
            "protected_dry_mass_kg": self.protected_dry_mass_kg,
            "nozzle_exit_area_m2": self.nozzle_exit_area_m2,
            "parameters": [
                {
                    "name": p.name, "reference_value": p.reference_value, "scale": p.scale,
                    "minimum": p.minimum, "maximum": p.maximum, "unit": p.unit,
                }
                for p in self.parameters
            ],
            "reference_ids": list(THROTTLEABLE_REFERENCE_IDS),
        }
        if include_hash:
            payload["model_hash"] = canonical_sha256(payload)
        return payload

    @property
    def model_hash(self) -> str:
        """Canonical SHA-256 content identity of the engine model."""

        return canonical_sha256(self.manifest(include_hash=False))

    def to_json(self, *, indent: int | None = 2) -> str:
        """Serialize the complete model document as sorted, NaN-free JSON.

        Parameters
        ----------
        indent : int or None, default 2
            JSON indentation; ``None`` gives a compact document.

        Returns
        -------
        str
            JSON text of :meth:`manifest`.
        """

        return json.dumps(
            self.manifest(), allow_nan=False, ensure_ascii=False, indent=indent, sort_keys=True
        )

    @classmethod
    def from_dict(cls, payload: Mapping[str, object]) -> ThrottleableEngineModel:
        """Rebuild a model from a document and verify its integrity.

        Parameters
        ----------
        payload : Mapping[str, object]
            Document produced by :meth:`manifest`.

        Returns
        -------
        ThrottleableEngineModel
            Reconstructed model.

        Raises
        ------
        InputError
            If the document is malformed or violates the model constraints,
            its hash does not match its contents, or its numerical policies
            or metadata differ from the rebuilt model.
        """

        _throttle_document(payload)
        try:
            def curve(item):
                return None if item is None else MonotoneCubicCurve(item["x"], item["y"])

            # _throttle_document has checked the JSON types of every member read
            # below; each stream fraction is a two-entry [name, value] array.
            model = cls(
                rated_thrust_n=cast(float, payload["rated_thrust_n"]),
                rated_mass_flow_kg_s=cast(float, payload["rated_mass_flow_kg_s"]),
                minimum_throttle=cast(float, payload["minimum_throttle"]),
                maximum_throttle=cast(float, payload["maximum_throttle"]),
                stream_fractions=tuple(
                    (cast(str, pair[0]), cast(float, pair[1]))
                    for pair in cast(list[list[object]], payload["stream_fractions"])
                ),
                thrust_curve=curve(payload["thrust_curve"]),
                flow_curve=curve(payload["flow_curve"]),
                state_constraints=tuple(
                    AdditionalStateConstraint(
                        cast(str, i["state_name"]),
                        cast(float, i["minimum_value"]),
                        cast(str, i["label"]),
                    )
                    for i in cast(list[dict[str, object]], payload["state_constraints"])
                ),
                protected_dry_mass_kg=cast(float | None, payload["protected_dry_mass_kg"]),
                nozzle_exit_area_m2=cast(float | None, payload["nozzle_exit_area_m2"]),
                parameters=_document_parameters(payload["parameters"]),
                model_id=cast(str, payload["model_id"]),
                schema=cast(str, payload["schema"]),
            )
        except (KeyError, TypeError, ValueError, AttributeError) as error:
            raise InputError("Malformed throttleable engine document.") from error
        if payload.get("model_hash") != model.model_hash:
            raise InputError("Throttleable engine hash does not match its contents.")
        _require_numerical_policies(payload, model.manifest(include_hash=False))
        _metadata(payload, model.manifest(include_hash=False), ('reference_ids','parameters'))
        return model

    @classmethod
    def from_json(cls, text_value: str) -> ThrottleableEngineModel:
        """Deserialize and integrity-check a JSON model document.

        Parameters
        ----------
        text_value : str
            JSON text produced by :meth:`to_json`.

        Returns
        -------
        ThrottleableEngineModel
            Reconstructed model.

        Raises
        ------
        InputError
            If the text is not a valid model document.
        """

        return _read_document(text_value, "Throttleable engine document", cls.from_dict)
