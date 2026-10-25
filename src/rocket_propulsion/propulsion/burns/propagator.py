"""Propagation-ready force, mass-flow, event, and sensitivity contract.

This module stops at the propulsion/astrodynamics boundary.  It does not own
epochs, frames, attitude laws, or orbit integration.  A numerical propagator
supplies relative time, current mass, a resolved thrust direction, and any
tracked tank states; the bridge returns simultaneous translational and mass
derivatives, analytic partials, and the event surfaces that make
discontinuities explicit, following the propulsion-model boundary of [1]_.
The tabulated profile is a transient engine history of the kind described in
[2]_; it is caller-supplied and not validated against hardware here.

References
----------
.. [1] OREKIT-13.1.5-PROPULSION-MODEL. Orekit 13.1.5, "PropulsionModel /
       EventDetector" API: force, mass-derivative, parameter-driver and
       event contract.
.. [2] NASA-20040000363. NASA, "Transient Mathematical Modeling for Liquid
       Rocket Engine Systems", https://ntrs.nasa.gov/citations/20040000363.
"""

from __future__ import annotations

import json
from collections.abc import Mapping
from dataclasses import dataclass, field
from enum import Enum
from fractions import Fraction
from math import isfinite
from types import MappingProxyType

from rocket_propulsion.core.errors import DomainError, InputError
from rocket_propulsion.core.validation import (
    read_json_document,
    real,
    real_mapping,
    records,
    text,
)

from .direction import _direction_geometry
from .exchange_documents import _BRIDGE_METADATA, _bridge_document, _metadata
from .mass_closure import _MASS_CLOSURE_POLICY, _closure_ratio, _inventory_ledger
from .numeric_outputs import _finite_rational, _fraction
from .performance_surface import (
    PerformanceGridLineSurface,
    PerformanceValiditySurface,
    RectilinearPerformanceSurface,
)
from .serialization import _require_numerical_policies, canonical_sha256
from .stream_performance import RectilinearStreamPerformanceSurface
from .tabulated import InterpolationPolicy, TabulatedBurnArtifact, TabulatedBurnSegment
from .time_coordinates import _TIME_POLICY, _epoch_profile_time

# Absolute margin (kg, or the bound state's unit) within which a mass or tank
# state counts as having reached its protected limit.  Shared by the bridge's
# clamp and by consumers that precompute the cut (the Lunaris adapter).
INVENTORY_LIMIT_TOLERANCE = 1e-12

PROPAGATOR_BRIDGE_SCHEMA = "rocket_propulsion_propagator_bridge_v1"
AMBIENT_PRESSURE_CONDITION = "ambient_pressure_pa"
PROPAGATOR_BRIDGE_REFERENCE_IDS = (
    "OREKIT-13.1.5-PROPULSION-MODEL",
    "NASA-20040000363",
    "CCSDS-502.0-B-3",
)


class BurnEventKind(str, Enum):
    """Kinds of root surfaces exposed to a numerical propagator.

    ``START``, ``PROFILE_BOUNDARY`` and ``STOP`` are time events of the
    tabulated profile; ``DRY_MASS`` and ``TANK_RESERVE`` are inventory events
    [1]_.

    References
    ----------
    .. [1] OREKIT-13.1.5-PROPULSION-MODEL. Orekit 13.1.5, "PropulsionModel /
           EventDetector" API: force, mass-derivative, parameter-driver and
           event contract.
    """

    START = "start"
    PROFILE_BOUNDARY = "profile_boundary"
    STOP = "stop"
    DRY_MASS = "dry_mass"
    TANK_RESERVE = "tank_reserve"


@dataclass(frozen=True, slots=True)
class EstimableParameter:
    """One bounded parameter-driver definition for orbit estimation.

    Parameters
    ----------
    name : str
        Parameter name; models address parameters by name, so an estimator
        cannot silently reorder a positional parameter vector [1]_.
    reference_value : float
        Nominal value, inside ``[minimum, maximum]``.
    scale : float
        Positive finite-difference/normalization scale; it does not alter the
        physical value.
    minimum, maximum : float
        Finite bounds with ``minimum < maximum``.
    unit : str
        Non-empty unit label (``"1"`` for dimensionless scales).

    Raises
    ------
    DomainError
        If a value is non-finite, the scale is not positive, the bounds are
        empty, the reference value is out of bounds, or the unit is blank.

    References
    ----------
    .. [1] OREKIT-13.1.5-PROPULSION-MODEL. Orekit 13.1.5, "PropulsionModel /
           EventDetector" API: force, mass-derivative, parameter-driver and
           event contract.
    """

    name: str
    reference_value: float
    scale: float
    minimum: float
    maximum: float
    unit: str

    def __post_init__(self) -> None:
        text(self.name, "Estimable parameter name")
        text(self.unit, "Estimable parameter unit")
        for field_name in ("reference_value", "scale", "minimum", "maximum"):
            object.__setattr__(
                self,
                field_name,
                real(getattr(self, field_name), f"Estimable parameter {field_name}"),
            )
        values = (self.reference_value, self.scale, self.minimum, self.maximum)
        if any(not isfinite(value) for value in values):
            raise DomainError("Estimable parameter values and bounds must be finite.")
        if self.scale <= 0.0 or self.minimum >= self.maximum:
            raise DomainError("Estimable parameter scale and bounds are invalid.")
        if not self.minimum <= self.reference_value <= self.maximum:
            raise DomainError("Estimable parameter reference value lies outside its bounds.")
        if not self.unit.strip():
            raise DomainError("Estimable parameter unit must not be empty.")

    def validate(self, value: float) -> float:
        """Return ``value`` if it is finite and inside the parameter bounds.

        Parameters
        ----------
        value : float
            Candidate parameter value.

        Returns
        -------
        float
            The validated value.

        Raises
        ------
        DomainError
            If the value is non-finite or outside ``[minimum, maximum]``.
        """

        numeric = real(value, f"Parameter {self.name!r}")
        if not self.minimum <= numeric <= self.maximum:
            raise DomainError(
                f"Parameter {self.name!r} must be finite and in [{self.minimum}, {self.maximum}]."
            )
        return numeric


@dataclass(frozen=True, slots=True)
class StreamStateBinding:
    """Map one artifact stream to one propagator additional-mass state.

    Parameters
    ----------
    stream_id : str
        Identifier of a stream in the tabulated artifact.
    state_name : str
        Name of the propagator's additional state (a tank mass) that the
        stream drains.

    Notes
    -----
    Several streams may drain the same state; their derivatives are summed,
    so each tank's mass derivative is minus the total flow drawn from it [1]_.

    References
    ----------
    .. [1] NASA-SP-8112. NASA, "Pressurization Systems for Liquid Rockets",
           https://ntrs.nasa.gov/citations/19760015212.
    """

    stream_id: str
    state_name: str

    def __post_init__(self) -> None:
        text(self.stream_id, "Stream identifier")
        text(self.state_name, "Additional-state name")


@dataclass(frozen=True, slots=True)
class AdditionalStateConstraint:
    """Minimum permitted value (reserve) of a tracked tank or propellant state.

    Parameters
    ----------
    state_name : str
        Name of a bound additional state.
    minimum_value : float
        Finite, non-negative reserve in the state's unit.
    label : str, optional
        Human-readable label of the limit.

    Raises
    ------
    DomainError
        If the minimum is negative or non-finite.

    Notes
    -----
    Reserve propellant is held back from the burn, as in pressure-fed
    propellant budgeting [1]_.

    References
    ----------
    .. [1] NASA-SP-8112. NASA, "Pressurization Systems for Liquid Rockets",
           https://ntrs.nasa.gov/citations/19760015212.
    """

    state_name: str
    minimum_value: float
    label: str = "tank reserve"

    def __post_init__(self) -> None:
        text(self.state_name, "State-constraint state name")
        text(self.label, "State-constraint label")
        object.__setattr__(
            self, "minimum_value", real(self.minimum_value, "State-constraint minimum")
        )
        if self.minimum_value < 0.0:
            raise DomainError("State-constraint minimum must be finite and non-negative.")


@dataclass(frozen=True, slots=True)
class BurnEventSurface:
    """A signed root function with explicit propulsion-stop semantics.

    Parameters
    ----------
    event_id : str
        Stable identifier of the event.
    kind : BurnEventKind
        Event category.
    stops_burn : bool
        Whether crossing the root inhibits this propulsion model.  It does not
        request termination of the orbit propagation.
    relative_time_s : float or None, optional
        Profile time of a time event, in s.
    state_name : str or None, optional
        Additional state of a tank-reserve event.
    threshold : float or None, optional
        Limit value of a dry-mass or reserve event.

    Notes
    -----
    Time surfaces use ``g = t - (t_event + ignition_time_bias)``; mass
    surfaces use ``g = current - minimum``.  A zero means the propagator lands
    exactly on the discontinuity and recomputes the right-hand derivative
    [1]_.  Roots are located from sign changes of ``g``, so a step-wise
    detector can miss multiple roots inside one step [2]_.

    References
    ----------
    .. [1] OREKIT-13.1.5-PROPULSION-MODEL. Orekit 13.1.5, "PropulsionModel /
           EventDetector" API: force, mass-derivative, parameter-driver and
           event contract.
    .. [2] SCIPY-1.14.1-IVP. SciPy 1.14.1, "solve_ivp": event sign-change
           detection; multiple roots within one step may be missed.
    """

    event_id: str
    kind: BurnEventKind
    stops_burn: bool
    relative_time_s: float | None = None
    state_name: str | None = None
    threshold: float | None = None

    def value(
        self,
        *,
        relative_time_s: float,
        vehicle_mass_kg: float,
        additional_states: Mapping[str, float],
        ignition_time_bias_s: float = 0.0,
    ) -> float:
        """Evaluate the signed event function ``g``; the surface holds no state.

        Parameters
        ----------
        relative_time_s : float
            Propagation time relative to nominal ignition, in s.
        vehicle_mass_kg : float
            Current vehicle mass, in kg.
        additional_states : Mapping of str to float
            Current additional-state values; must contain ``state_name`` for
            a reserve event.
        ignition_time_bias_s : float, optional
            Ignition-time bias added to a time event, in s.

        Returns
        -------
        float
            ``t - (t_event + bias)`` for a time event, ``m - threshold`` for a
            dry-mass event, or ``state - threshold`` for a reserve event.

        Raises
        ------
        DomainError
            If an input is non-finite, the surface is malformed, or the
            required state is missing.
        """

        relative_time_s = real(relative_time_s, "Event evaluation time")
        vehicle_mass_kg = real(vehicle_mass_kg, "Event evaluation vehicle mass")
        ignition_time_bias_s = real(ignition_time_bias_s, "Event ignition-time bias")
        if self.relative_time_s is not None:
            return relative_time_s - (self.relative_time_s + ignition_time_bias_s)
        if self.kind is BurnEventKind.DRY_MASS:
            assert self.threshold is not None
            return vehicle_mass_kg - self.threshold
        if self.state_name is None or self.threshold is None:
            raise DomainError("Malformed state event surface.")
        states = real_mapping(additional_states, "Additional states")
        if self.state_name not in states:
            raise DomainError(f"Missing additional state {self.state_name!r}.")
        return states[self.state_name] - self.threshold


@dataclass(frozen=True, slots=True)
class PropulsionPartials:
    """Analytic local derivatives for variational equations and estimation.

    Attributes
    ----------
    acceleration_wrt_mass : tuple of float
        ``da/dm``, in m/(s**2 kg).
    acceleration_wrt_parameters, mass_rate_wrt_parameters : tuple
        Named partials with respect to the estimable parameters.
    additional_state_rates_wrt_parameters : tuple
        Per-parameter partials of each bound state's drain.
    acceleration_wrt_conditions, mass_rate_wrt_conditions : tuple
        Named partials with respect to performance-surface conditions.
    additional_state_rates_wrt_conditions : tuple
        Per-condition partials of each bound state's drain.
    nonsmooth_parameters : tuple of str
        Parameters whose derivative does not exist at this point.
    acceleration_wrt_direction : tuple of tuple of float
        3x3 Jacobian ``J[i][j] = da_i/dd_j``.
    acceleration_wrt_ambient_pressure : tuple of float
        ``da/dp_a``, in m/(s**2 Pa).

    Notes
    -----
    Acceleration partials are expressed in the same resolved frame as the
    supplied direction [1]_.  At a command discontinuity the ignition-time
    partial is reported as non-smooth instead of an invented derivative.

    ``acceleration_wrt_direction`` is taken with respect to the *supplied*
    (not necessarily unit) direction vector ``d``.  Because the bridge
    normalizes ``d``, ``J = (T/m) (I - u u^T) / |d|``, where ``T`` is the
    signed axial thrust and ``u = d / |d|``: only the component of a
    perturbation perpendicular to the thrust axis changes the acceleration.
    A consumer that parameterizes pointing by angles chains this with
    ``dd/dangle``.

    References
    ----------
    .. [1] OREKIT-13.1.5-PROPULSION-MODEL. Orekit 13.1.5, "PropulsionModel /
           EventDetector" API: force, mass-derivative, parameter-driver and
           event contract.
    """

    acceleration_wrt_mass: tuple[float, float, float]
    acceleration_wrt_parameters: tuple[tuple[str, tuple[float, float, float]], ...]
    mass_rate_wrt_parameters: tuple[tuple[str, float], ...]
    additional_state_rates_wrt_parameters: tuple[tuple[str, tuple[tuple[str, float], ...]], ...]
    acceleration_wrt_conditions: tuple[tuple[str, tuple[float, float, float]], ...] = ()
    mass_rate_wrt_conditions: tuple[tuple[str, float], ...] = ()
    additional_state_rates_wrt_conditions: tuple[
        tuple[str, tuple[tuple[str, float], ...]], ...
    ] = ()
    nonsmooth_parameters: tuple[str, ...] = ()
    acceleration_wrt_direction: tuple[
        tuple[float, float, float],
        tuple[float, float, float],
        tuple[float, float, float],
    ] = ((0.0, 0.0, 0.0), (0.0, 0.0, 0.0), (0.0, 0.0, 0.0))
    acceleration_wrt_ambient_pressure: tuple[float, float, float] = (0.0, 0.0, 0.0)


@dataclass(frozen=True, slots=True)
class PropulsionDynamicsEvaluation:
    """Simultaneous force and inventory derivative at one propagator call.

    Attributes
    ----------
    relative_time_s : float
        Queried time relative to nominal ignition, in s.
    effective_profile_time_s : float
        Profile time after the ignition-time bias, in s.
    active : bool
        Whether the burn produces thrust at this point.
    phase : str
        Profile phase label, or ``"off"`` / ``"inhibited"``.
    axial_thrust_n : float
        Signed net axial thrust, in N.
    force_vector_n, acceleration_m_s2 : tuple of float
        Force (N) and acceleration (m/s**2) along the resolved direction.
    mass_derivative_kg_s : float
        Vehicle mass derivative (negative drain), in kg/s.
    additional_state_derivatives : tuple of (str, float)
        Derivative of each bound state, in its unit per second.
    parameters : tuple of (str, float)
        Resolved estimable-parameter values.
    inhibited_by : tuple of str
        Inventory limits that inhibited the burn.
    partials : PropulsionPartials
        Analytic local derivatives.
    operating_conditions : tuple of (str, float)
        Performance-surface conditions used.
    performance_force_scale, performance_mass_flow_scale : float
        Multipliers applied by the performance surface.
    ambient_pressure_pa : float or None
        Ambient pressure used for the pressure-thrust term, in Pa.
    pressure_thrust_loss_n : float
        ``p_a A_e`` subtracted from the vacuum thrust, in N.
    net_thrust_margin_n : float or None
        Net-thrust root monitored by the consumer while firing, in N.

    Notes
    -----
    The record follows the propulsion-model boundary of [1]_ for a transient
    engine history [2]_.

    References
    ----------
    .. [1] OREKIT-13.1.5-PROPULSION-MODEL. Orekit 13.1.5, "PropulsionModel /
           EventDetector" API: force, mass-derivative, parameter-driver and
           event contract.
    .. [2] NASA-20040000363. NASA, "Transient Mathematical Modeling for Liquid
           Rocket Engine Systems", https://ntrs.nasa.gov/citations/20040000363.
    """

    relative_time_s: float
    effective_profile_time_s: float
    active: bool
    phase: str
    axial_thrust_n: float
    force_vector_n: tuple[float, float, float]
    acceleration_m_s2: tuple[float, float, float]
    mass_derivative_kg_s: float
    additional_state_derivatives: tuple[tuple[str, float], ...]
    parameters: tuple[tuple[str, float], ...]
    inhibited_by: tuple[str, ...]
    partials: PropulsionPartials
    operating_conditions: tuple[tuple[str, float], ...] = ()
    performance_force_scale: float = 1.0
    performance_mass_flow_scale: float = 1.0
    ambient_pressure_pa: float | None = None
    pressure_thrust_loss_n: float = 0.0
    net_thrust_margin_n: float | None = None


@dataclass(frozen=True, slots=True)
class PropellantMassClosure:
    """Vehicle-mass versus tracked-propellant bookkeeping residual.

    Parameters
    ----------
    vehicle_mass_kg : float
        Positive propagated vehicle mass, in kg.
    non_propellant_mass_kg : float
        Non-negative mass not represented by a bound state, in kg.
    tracked_propellant_kg : float
        Sum of the bound propellant states, in kg.
    residual_kg : float
        ``vehicle_mass - non_propellant_mass - tracked_propellant``, in kg.
    relative_residual : float
        ``residual_kg / vehicle_mass_kg``.
    residual_ratio : tuple of (str, str) or None, optional
        Exact residual as integer numerator and denominator text [3]_; when
        omitted, the float ``residual_kg`` is taken as exact.

    Raises
    ------
    DomainError
        If a mass is invalid, the ratio is malformed or contradicts
        ``residual_kg``, or ``relative_residual`` contradicts the ratio.

    Notes
    -----
    With complete stream binding and consistent initial conditions the
    residual stays at integration round-off for the whole burn, because the
    total drain equals the sum of the bound stream drains [1]_, [2]_.  A
    growing residual reveals inconsistent initial states, an unbound stream,
    or a consumer that integrates the vehicle mass and tank states with
    different derivatives.

    References
    ----------
    .. [1] NASA-SP-8112. NASA, "Pressurization Systems for Liquid Rockets",
           https://ntrs.nasa.gov/citations/19760015212.
    .. [2] OREKIT-13.1.5-PROPULSION-MODEL. Orekit 13.1.5, "PropulsionModel /
           EventDetector" API: force, mass-derivative, parameter-driver and
           event contract.
    .. [3] PYTHON-3.12-FRACTION. Python 3.12, "fractions": float inputs retain
           their exact binary value as a rational,
           https://docs.python.org/3.12/library/fractions.html.
    """

    vehicle_mass_kg: float
    non_propellant_mass_kg: float
    tracked_propellant_kg: float
    residual_kg: float
    relative_residual: float
    residual_ratio: tuple[str, str] | None = field(default=None, repr=False)

    def __post_init__(self) -> None:
        for name, label in (("vehicle_mass_kg", "Vehicle mass"),
                            ("non_propellant_mass_kg", "Non-propellant mass"),
                            ("tracked_propellant_kg", "Tracked propellant mass"),
                            ("residual_kg", "Mass closure residual"),
                            ("relative_residual", "Relative mass closure residual")):
            object.__setattr__(self, name, real(getattr(self, name), label))
        if self.vehicle_mass_kg <= 0 or self.non_propellant_mass_kg < 0:
            raise DomainError("Vehicle mass must be positive and other mass non-negative.")
        ratio = _closure_ratio(self.residual_ratio, self.residual_kg)
        expected = _finite_rational(ratio / _fraction(self.vehicle_mass_kg), "Relative mass closure residual")
        if self.relative_residual != expected:
            raise DomainError("Relative residual contradicts the numerical ratio.",
                              field="Relative mass closure residual")
        object.__setattr__(self, "residual_ratio", (str(ratio.numerator), str(ratio.denominator)))

    def is_closed(self, *, absolute_tolerance_kg: float) -> bool:
        """Return whether the exact residual satisfies ``|residual| <= tolerance``.

        Parameters
        ----------
        absolute_tolerance_kg : float
            Finite, non-negative tolerance, in kg.

        Returns
        -------
        bool
            Comparison of the exact, unrounded residual with the tolerance.

        Raises
        ------
        DomainError
            If the tolerance is negative or non-finite.
        """

        absolute_tolerance_kg = real(absolute_tolerance_kg, "Mass-closure tolerance")
        if absolute_tolerance_kg < 0.0:
            raise DomainError("Mass-closure tolerance must be non-negative.", field="Mass-closure tolerance")
        exact_residual = Fraction(int(self.residual_ratio[0]), int(self.residual_ratio[1]))
        return abs(exact_residual) <= _fraction(absolute_tolerance_kg)


def _validated_states(additional_states: Mapping[str, float] | None) -> dict[str, float]:
    """Return the additional states as a validated ``{name: finite float}`` dict."""
    return real_mapping(additional_states, "Additional states")


def _mass_profile(segment: TabulatedBurnSegment, time_s: float, hold: bool):
    """Return exact total and per-stream flows and their time slopes on a segment.

    Each channel is linearly interpolated between the segment end points,
    ``mdot(t) = mdot_0 + (t - t_0) / (t_1 - t_0) * (mdot_1 - mdot_0)``, or held
    at its start value with zero slope when ``hold`` is set.  Values stay
    exact rationals [1]_ so grouping and scaling round only once.

    References
    ----------
    .. [1] PYTHON-3.12-FRACTION. Python 3.12, "fractions": float inputs retain
           their exact binary value as a rational,
           https://docs.python.org/3.12/library/fractions.html.
    """
    start, end = _fraction(segment.start_time_s), _fraction(segment.end_time_s)
    duration = end - start
    fraction = Fraction(0) if hold else (_fraction(time_s) - start) / duration

    def channel(first, last):
        initial, final = _fraction(first), _fraction(last)
        return (initial + fraction * (final - initial),
                Fraction(0) if hold else (final - initial) / duration)

    total, slope = channel(segment.start_total_mass_flow_kg_s, segment.end_total_mass_flow_kg_s)
    initial = dict(segment.start_stream_mass_flows_kg_s)
    final = dict(segment.end_stream_mass_flows_kg_s)
    samples = {name: channel(initial.get(name, 0.0), final.get(name, 0.0))
               for name in sorted(initial.keys() | final.keys())}
    return total, {name: pair[0] for name, pair in samples.items()}, slope, {
        name: pair[1] for name, pair in samples.items()}


@dataclass(frozen=True, slots=True)
class PropulsionDerivativeJump:
    """Left-minus-right propulsion derivatives at one time event.

    Attributes
    ----------
    relative_time_s : float
        Event time relative to nominal ignition, in s.
    acceleration_m_s2 : tuple of float
        ``a(left) - a(right)``, in m/s**2.
    mass_derivative_kg_s : float
        ``mdot(left) - mdot(right)``, in kg/s.
    additional_state_derivatives : tuple of (str, float)
        Left-minus-right derivative of each bound state.

    Notes
    -----
    For a consumer integrating sensitivities ``S = dy/dp``, an event whose
    time ``tau`` depends on a parameter requires the saltation update
    ``S+ = S- + (f- - f+) dtau/dp``; maneuver execution errors in navigation
    filtering are treated in [1]_.  Profile time events move one-for-one with
    ``ignition_time_bias_s`` (``dtau/db = 1``) and not with the scale drivers,
    so the propulsion part of ``f- - f+`` is exactly this jump.

    References
    ----------
    .. [1] NASA-20180003657. NASA, "Navigation Filter Best Practices"
           (NASA/TP-2018-219822), https://ntrs.nasa.gov/citations/20180003657.
    """

    relative_time_s: float
    acceleration_m_s2: tuple[float, float, float]
    mass_derivative_kg_s: float
    additional_state_derivatives: tuple[tuple[str, float], ...]


@dataclass(frozen=True, slots=True)
class PropagatorPropulsionBridge:
    """Immutable tabulated propulsion model for numerical orbit propagation.

    Parameters
    ----------
    artifact : TabulatedBurnArtifact
        Hash-verified thrust and mass-flow history; re-opened through
        :meth:`TabulatedBurnArtifact.from_dict` before it is accepted.
    performance_surface : RectilinearPerformanceSurface or None, optional
        Total-flow multipliers over operating conditions.
    stream_performance_surface : RectilinearStreamPerformanceSurface or None, optional
        Stream-resolved multipliers; exclusive with ``performance_surface``.
        Its stream identifiers must equal the artifact's.
    stream_bindings : tuple of StreamStateBinding, optional
        At most one binding per artifact stream.
    state_constraints : tuple of AdditionalStateConstraint, optional
        At most one reserve per bound state.
    protected_dry_mass_kg : float or None, optional
        Finite, positive dry-mass limit, in kg.
    require_complete_stream_binding : bool, optional
        Require every artifact stream to be bound.
    parameters : tuple of EstimableParameter, optional
        Exactly ``thrust_scale``, ``mass_flow_scale`` and
        ``ignition_time_bias_s``; defaults are created when empty.
    nozzle_exit_area_m2 : float or None, optional
        Finite, positive axial-projected total nozzle exit area, in m**2.
    model_id : str, optional
        Model identifier.
    schema : str, optional
        Must equal ``PROPAGATOR_BRIDGE_SCHEMA``.

    Raises
    ------
    DomainError
        If a field is invalid, both performance surfaces are given, bindings
        or constraints are duplicated or reference unknown streams or states,
        complete binding is required but missing, ambient pressure is also a
        performance-surface axis, or the parameter set is not exactly the
        three expected drivers.
    InputError
        If the artifact fails its integrity check.

    Notes
    -----
    The returned force is the artifact's *axial* force: cant and divergence
    already represented by the producer are not applied twice.  Vehicle and
    per-stream mass derivatives are negative drains.  ``thrust_scale`` and
    ``mass_flow_scale`` are independent so orbit determination can estimate
    force and propellant bookkeeping without imposing an Isp constraint;
    calibrated upstream realizations normally leave both at one.

    With ``nozzle_exit_area_m2`` set, the artifact force is the vacuum thrust
    and the ambient-pressure term is subtracted at evaluation time,
    ``F = F_vac - p_a A_e`` [3]_.

    ``ignition_time_bias_s`` shifts the entire profile and all date events by
    the same amount.  Its derivative is analytic inside linear segments and is
    marked non-smooth at every declared boundary.  The force, mass-derivative,
    parameter-driver and event contract follows [1]_; the profile is a
    caller-supplied transient history [2]_.

    References
    ----------
    .. [1] OREKIT-13.1.5-PROPULSION-MODEL. Orekit 13.1.5, "PropulsionModel /
           EventDetector" API: force, mass-derivative, parameter-driver and
           event contract.
    .. [2] NASA-20040000363. NASA, "Transient Mathematical Modeling for Liquid
           Rocket Engine Systems", https://ntrs.nasa.gov/citations/20040000363.
    .. [3] NASA-STEADY-ROCKET-THRUST. NASA Glenn Research Center, momentum and
           pressure thrust equation,
           https://www1.grc.nasa.gov/beginners-guide-to-aeronautics/rocket-thrust/.
    """

    artifact: TabulatedBurnArtifact
    performance_surface: RectilinearPerformanceSurface | None = None
    stream_performance_surface: RectilinearStreamPerformanceSurface | None = None
    stream_bindings: tuple[StreamStateBinding, ...] = ()
    state_constraints: tuple[AdditionalStateConstraint, ...] = ()
    protected_dry_mass_kg: float | None = None
    require_complete_stream_binding: bool = True
    parameters: tuple[EstimableParameter, ...] = ()
    nozzle_exit_area_m2: float | None = None
    model_id: str = "tabulated-propulsion-dynamics"
    schema: str = PROPAGATOR_BRIDGE_SCHEMA

    def __post_init__(self) -> None:
        if self.schema != PROPAGATOR_BRIDGE_SCHEMA:
            raise DomainError("Unsupported propagator bridge schema.")
        if not isinstance(self.artifact, TabulatedBurnArtifact):
            raise DomainError("Propagator bridge requires a TabulatedBurnArtifact.")
        if self.performance_surface is not None and not isinstance(
            self.performance_surface, RectilinearPerformanceSurface
        ):
            raise DomainError("Performance surface must be a RectilinearPerformanceSurface.")
        if self.stream_performance_surface is not None and not isinstance(
            self.stream_performance_surface, RectilinearStreamPerformanceSurface
        ):
            raise DomainError(
                "Stream performance surface must be a RectilinearStreamPerformanceSurface."
            )
        if self.performance_surface is not None and self.stream_performance_surface is not None:
            raise DomainError(
                "A bridge may use either a total-flow or stream-resolved performance "
                "surface, not both."
            )
        for field_name, kind in (
            ("stream_bindings", StreamStateBinding),
            ("state_constraints", AdditionalStateConstraint),
            ("parameters", EstimableParameter),
        ):
            object.__setattr__(
                self, field_name, records(getattr(self, field_name), kind, f"Bridge {field_name}")
            )
        for field_name in ("protected_dry_mass_kg", "nozzle_exit_area_m2"):
            if getattr(self, field_name) is not None:
                object.__setattr__(
                    self, field_name, real(getattr(self, field_name), f"Bridge {field_name}")
                )
        if not isinstance(self.require_complete_stream_binding, bool):
            raise DomainError("require_complete_stream_binding must be a boolean.")
        text(self.model_id, "Propagator bridge model identifier")
        # Re-open at the public integrity boundary before accepting the source.
        TabulatedBurnArtifact.from_dict(self.artifact.to_dict())
        if self.protected_dry_mass_kg is not None and (
            not isfinite(self.protected_dry_mass_kg) or self.protected_dry_mass_kg <= 0.0
        ):
            raise DomainError("Protected dry mass must be finite and greater than zero.")
        if self.nozzle_exit_area_m2 is not None:
            if not isfinite(self.nozzle_exit_area_m2) or self.nozzle_exit_area_m2 <= 0.0:
                raise DomainError("Nozzle exit area must be finite and greater than zero.")
            active_surface = self.performance_surface or self.stream_performance_surface
            if active_surface is not None and any(
                axis.name == AMBIENT_PRESSURE_CONDITION for axis in active_surface.axes
            ):
                raise DomainError(
                    "Ambient pressure is applied through the nozzle exit area; it must not "
                    "also be a performance-surface axis."
                )

        binding_ids = [binding.stream_id for binding in self.stream_bindings]
        if len(binding_ids) != len(set(binding_ids)):
            raise DomainError("Each artifact stream may have only one state binding.")
        constraint_names = [constraint.state_name for constraint in self.state_constraints]
        if len(constraint_names) != len(set(constraint_names)):
            raise DomainError("Each additional state may have only one constraint.")
        bound_states = {binding.state_name for binding in self.stream_bindings}
        missing_constraint_states = set(constraint_names) - bound_states
        if missing_constraint_states:
            raise DomainError(
                "State constraints require a stream binding: "
                + ", ".join(sorted(missing_constraint_states))
            )

        artifact_streams = {
            stream_id
            for segment in self.artifact.segments
            for stream_id, _ in (
                segment.start_stream_mass_flows_kg_s + segment.end_stream_mass_flows_kg_s
            )
        }
        unknown_bindings = set(binding_ids) - artifact_streams
        if unknown_bindings:
            raise DomainError(
                "Bindings reference unknown artifact streams: "
                + ", ".join(sorted(unknown_bindings))
            )
        if self.require_complete_stream_binding and artifact_streams != set(binding_ids):
            missing = artifact_streams - set(binding_ids)
            raise DomainError(
                "Every artifact stream must be bound when complete binding is required: "
                + ", ".join(sorted(missing))
            )
        if self.stream_performance_surface is not None and (
            set(self.stream_performance_surface.stream_ids) != artifact_streams
        ):
            raise DomainError("Stream-performance outputs must match the artifact streams exactly.")

        if not self.parameters:
            limit = max(1.0, self.artifact.duration_s)
            object.__setattr__(
                self,
                "parameters",
                (
                    EstimableParameter("thrust_scale", 1.0, 1e-3, 0.01, 2.0, "1"),
                    EstimableParameter("mass_flow_scale", 1.0, 1e-3, 0.01, 2.0, "1"),
                    EstimableParameter("ignition_time_bias_s", 0.0, 1e-3, -limit, limit, "s"),
                ),
            )
        names = [parameter.name for parameter in self.parameters]
        expected = {"thrust_scale", "mass_flow_scale", "ignition_time_bias_s"}
        if set(names) != expected or len(names) != len(expected):
            raise DomainError(
                "Propagator bridge requires exactly thrust_scale, mass_flow_scale, "
                "and ignition_time_bias_s parameter definitions."
            )

    @property
    def model_hash(self) -> str:
        """SHA-256 content identity of the source artifact and propagation semantics."""

        return canonical_sha256(self.manifest(include_hash=False))

    def manifest(self, *, include_hash: bool = True) -> dict[str, object]:
        """Return the deterministic model manifest recorded with integrations.

        Parameters
        ----------
        include_hash : bool, optional
            Add the ``model_hash`` member computed over the other members.

        Returns
        -------
        dict
            Source hashes, semantics and numerical-policy labels, bindings,
            constraints, parameters and reference identifiers.
        """

        payload: dict[str, object] = {
            "schema": self.schema,
            "model_id": self.model_id,
            "artifact_hash": self.artifact.artifact_hash,
            "performance_surface_hash": (
                None if self.performance_surface is None else self.performance_surface.surface_hash
            ),
            "stream_performance_surface_hash": (
                None
                if self.stream_performance_surface is None
                else self.stream_performance_surface.surface_hash
            ),
            "operating_condition_axes": (
                []
                if self.performance_surface is None and self.stream_performance_surface is None
                else [
                    {
                        "name": axis.name,
                        "unit": axis.unit,
                        "minimum": axis.minimum,
                        "maximum": axis.maximum,
                    }
                    for axis in (
                        self.performance_surface.axes
                        if self.performance_surface is not None
                        else self.stream_performance_surface.axes
                    )
                ]
            ),
            "force_semantics": "net_axial_thrust_along_resolved_direction",
            "mass_semantics": "negative_total_and_bound_stream_tank_drain",
            "mass_flow_numerics": "exact_linear_mass_sampling_and_grouped_outputs_v1",
            "tabulated_state_numerics": "exact_binary_axial_and_grouped_flow_channels_v1",
            "mass_closure_numerics": _MASS_CLOSURE_POLICY,
            "inventory_clamp_numerics": "exact_binary_state_and_limit_plus_declared_margin_v1",
            "direction_numerics": "scaled_hypot_rational_projection_and_signed_tangents_v1",
            "time_semantics": "relative_seconds_with_global_ignition_bias",
            "time_numerics": _TIME_POLICY,
            "stream_bindings": [
                {"stream_id": item.stream_id, "state_name": item.state_name}
                for item in self.stream_bindings
            ],
            "state_constraints": [
                {
                    "state_name": item.state_name,
                    "minimum_value": item.minimum_value,
                    "label": item.label,
                }
                for item in self.state_constraints
            ],
            "protected_dry_mass_kg": self.protected_dry_mass_kg,
            "require_complete_stream_binding": self.require_complete_stream_binding,
            "parameters": [
                {
                    "name": item.name,
                    "reference_value": item.reference_value,
                    "scale": item.scale,
                    "minimum": item.minimum,
                    "maximum": item.maximum,
                    "unit": item.unit,
                }
                for item in self.parameters
            ],
            "reference_ids": list(PROPAGATOR_BRIDGE_REFERENCE_IDS),
        }
        if self.nozzle_exit_area_m2 is not None:
            # Present only when used, so a vacuum-only bridge identity does not
            # depend on this field.
            payload["nozzle_exit_area_m2"] = self.nozzle_exit_area_m2
            payload["force_semantics"] = (
                "vacuum_axial_thrust_minus_ambient_pressure_times_exit_area"
            )
        if include_hash:
            payload["model_hash"] = canonical_sha256(payload)
        return payload

    def to_dict(self) -> dict[str, object]:
        """Return the complete, reopenable bridge document.

        Returns
        -------
        dict
            The manifest (with ``model_hash``) extended with the full source
            artifact and optional performance surfaces, each carrying its own
            content hash.

        Notes
        -----
        :meth:`from_dict` rebuilds the bridge and refuses any document whose
        recomputed ``model_hash`` differs, so a recorded model is reproducible
        as required for simulation records [1]_.

        References
        ----------
        .. [1] NASA-STD-7009B. NASA, "Standard for Models and Simulations",
               https://standards.nasa.gov/standard/NASA/NASA-STD-7009.
        """

        document = self.manifest()
        document["artifact"] = self.artifact.to_dict()
        document["performance_surface"] = (
            None if self.performance_surface is None else self.performance_surface.to_dict()
        )
        document["stream_performance_surface"] = (
            None
            if self.stream_performance_surface is None
            else self.stream_performance_surface.to_dict()
        )
        return document

    def to_json(self, *, indent: int | None = 2) -> str:
        """Serialize the complete bridge document as JSON text.

        Parameters
        ----------
        indent : int or None, optional
            Indentation passed to :func:`json.dumps`; ``None`` is compact.

        Returns
        -------
        str
            Key-sorted JSON without ``NaN`` or ``Infinity`` tokens [1]_.

        References
        ----------
        .. [1] IETF RFC 8259, "The JavaScript Object Notation (JSON) Data
               Interchange Format", https://www.rfc-editor.org/rfc/rfc8259.
        """

        return json.dumps(
            self.to_dict(), allow_nan=False, ensure_ascii=False, indent=indent, sort_keys=True
        )

    @classmethod
    def from_dict(cls, payload: Mapping[str, object]) -> PropagatorPropulsionBridge:
        """Rebuild a bridge and verify every nested hash and the model hash.

        Parameters
        ----------
        payload : Mapping
            Document produced by :meth:`to_dict`.

        Returns
        -------
        PropagatorPropulsionBridge
            The reconstructed bridge.

        Raises
        ------
        InputError
            If the document is malformed, a nested or model hash does not
            match, or its numerical-policy or metadata members differ from the
            rebuilt manifest.
        DomainError
            If a reconstructed field violates the bridge constraints.

        Notes
        -----
        Hash verification makes a recorded model reproducible [1]_.

        References
        ----------
        .. [1] NASA-STD-7009B. NASA, "Standard for Models and Simulations",
               https://standards.nasa.gov/standard/NASA/NASA-STD-7009.
        """

        _bridge_document(payload)
        try:
            surface_payload = payload["performance_surface"]
            stream_surface_payload = payload.get("stream_performance_surface")
            bridge = cls(
                artifact=TabulatedBurnArtifact.from_dict(payload["artifact"]),
                performance_surface=(
                    None
                    if surface_payload is None
                    else RectilinearPerformanceSurface.from_dict(surface_payload)
                ),
                stream_performance_surface=(
                    None
                    if stream_surface_payload is None
                    else RectilinearStreamPerformanceSurface.from_dict(stream_surface_payload)
                ),
                stream_bindings=tuple(
                    StreamStateBinding(item["stream_id"], item["state_name"])
                    for item in payload["stream_bindings"]
                ),
                state_constraints=tuple(
                    AdditionalStateConstraint(
                        item["state_name"], item["minimum_value"], item["label"]
                    )
                    for item in payload["state_constraints"]
                ),
                protected_dry_mass_kg=payload["protected_dry_mass_kg"],
                require_complete_stream_binding=payload["require_complete_stream_binding"],
                parameters=tuple(
                    EstimableParameter(
                        item["name"],
                        item["reference_value"],
                        item["scale"],
                        item["minimum"],
                        item["maximum"],
                        item["unit"],
                    )
                    for item in payload["parameters"]
                ),
                nozzle_exit_area_m2=payload.get("nozzle_exit_area_m2"),
                model_id=payload["model_id"],
                schema=payload["schema"],
            )
        except InputError:
            raise
        except (KeyError, TypeError, ValueError, AttributeError) as error:
            raise InputError("Malformed propagator bridge document.") from error
        if payload.get("model_hash") != bridge.model_hash:
            raise InputError("Propagator bridge model hash does not match its contents.")
        _require_numerical_policies(payload, bridge.manifest(include_hash=False))
        _metadata(payload, bridge.manifest(include_hash=False), _BRIDGE_METADATA)
        return bridge

    @classmethod
    def from_json(cls, text: str) -> PropagatorPropulsionBridge:
        """Parse JSON text and rebuild the bridge with :meth:`from_dict`.

        Parameters
        ----------
        text : str
            RFC 8259 JSON document [1]_.

        Returns
        -------
        PropagatorPropulsionBridge
            The reconstructed, hash-verified bridge.

        Raises
        ------
        InputError
            If the text is not a valid document or fails a hash check.

        References
        ----------
        .. [1] IETF RFC 8259, "The JavaScript Object Notation (JSON) Data
               Interchange Format", https://www.rfc-editor.org/rfc/rfc8259.
        """

        return read_json_document(text, "Propagator bridge document", cls.from_dict)

    def resolved_parameters(
        self, overrides: Mapping[str, float] | None = None
    ) -> Mapping[str, float]:
        """Resolve named overrides against the bounded parameter definitions.

        Parameters
        ----------
        overrides : Mapping of str to float or None, optional
            Values replacing the reference values of named parameters.

        Returns
        -------
        Mapping of str to float
            Read-only mapping of every parameter to its validated value.

        Raises
        ------
        DomainError
            If a name is unknown or a value is non-finite or out of bounds.

        References
        ----------
        .. [1] OREKIT-13.1.5-PROPULSION-MODEL. Orekit 13.1.5, "PropulsionModel /
               EventDetector" API: force, mass-derivative, parameter-driver and
               event contract.
        """

        supplied = real_mapping(overrides, "Parameter overrides")
        known = {parameter.name for parameter in self.parameters}
        unknown = set(supplied) - known
        if unknown:
            raise DomainError("Unknown propulsion parameters: " + ", ".join(sorted(unknown)))
        resolved = {
            parameter.name: parameter.validate(
                supplied.get(parameter.name, parameter.reference_value)
            )
            for parameter in self.parameters
        }
        return MappingProxyType(resolved)

    def event_surfaces(self) -> tuple[BurnEventSurface, ...]:
        """Return every time and inventory discontinuity exactly once.

        Returns
        -------
        tuple of BurnEventSurface
            Start, each interior profile boundary, stop, and, when declared,
            the dry-mass and per-state reserve surfaces [1]_.

        References
        ----------
        .. [1] OREKIT-13.1.5-PROPULSION-MODEL. Orekit 13.1.5, "PropulsionModel /
               EventDetector" API: force, mass-derivative, parameter-driver and
               event contract.
        """

        surfaces = [BurnEventSurface("burn:start", BurnEventKind.START, False, relative_time_s=0.0)]
        for index, segment in enumerate(self.artifact.segments[:-1], start=1):
            surfaces.append(
                BurnEventSurface(
                    f"burn:boundary:{index}",
                    BurnEventKind.PROFILE_BOUNDARY,
                    False,
                    relative_time_s=segment.end_time_s,
                )
            )
        surfaces.append(
            BurnEventSurface(
                "burn:stop",
                BurnEventKind.STOP,
                True,
                relative_time_s=self.artifact.duration_s,
            )
        )
        if self.protected_dry_mass_kg is not None:
            surfaces.append(
                BurnEventSurface(
                    "burn:dry-mass",
                    BurnEventKind.DRY_MASS,
                    True,
                    threshold=self.protected_dry_mass_kg,
                )
            )
        surfaces.extend(
            BurnEventSurface(
                f"burn:reserve:{constraint.state_name}",
                BurnEventKind.TANK_RESERVE,
                True,
                state_name=constraint.state_name,
                threshold=constraint.minimum_value,
            )
            for constraint in self.state_constraints
        )
        return tuple(surfaces)

    def operating_condition_surfaces(self) -> tuple[PerformanceValiditySurface, ...]:
        """Return the validity-domain roots of the optional performance surface.

        Returns
        -------
        tuple of PerformanceValiditySurface
            Roots at the edges of the surface's declared condition domain;
            empty without a surface.

        Notes
        -----
        A model is used only inside its stated domain of validity [1]_; the
        roots are registered with the propagator's event detection [2]_.

        References
        ----------
        .. [1] NASA-STD-7009B. NASA, "Standard for Models and Simulations",
               https://standards.nasa.gov/standard/NASA/NASA-STD-7009.
        .. [2] OREKIT-13.1.5-PROPULSION-MODEL. Orekit 13.1.5, "PropulsionModel /
               EventDetector" API: force, mass-derivative, parameter-driver and
               event contract.
        """

        if self.performance_surface is None and self.stream_performance_surface is None:
            return ()
        if self.performance_surface is not None:
            return self.performance_surface.validity_surfaces()
        assert self.stream_performance_surface is not None
        return self.stream_performance_surface.validity_surfaces()

    def mass_closure(
        self,
        *,
        vehicle_mass_kg: float,
        additional_states: Mapping[str, float],
        non_propellant_mass_kg: float,
    ) -> PropellantMassClosure:
        """Compare propagated vehicle mass with the bound tank states.

        Parameters
        ----------
        vehicle_mass_kg : float
            Finite, positive propagated vehicle mass, in kg.
        additional_states : Mapping of str to float
            Current value of every bound state, in kg.
        non_propellant_mass_kg : float
            Finite, non-negative mass not represented by a bound state (dry
            structure, pressurant, unbound or untracked propellant), in kg.

        Returns
        -------
        PropellantMassClosure
            Residual ``m_vehicle - m_other - sum(bound states)``; each bound
            state is counted once even if several streams drain it.

        Raises
        ------
        DomainError
            If a mass is invalid or a bound state is missing.

        References
        ----------
        .. [1] NASA-SP-8112. NASA, "Pressurization Systems for Liquid Rockets",
               https://ntrs.nasa.gov/citations/19760015212.
        """

        vehicle_mass_kg = real(vehicle_mass_kg, "Vehicle mass")
        non_propellant_mass_kg = real(non_propellant_mass_kg, "Non-propellant mass")
        for value, label in (
            (vehicle_mass_kg, "Vehicle mass"),
            (non_propellant_mass_kg, "Non-propellant mass"),
        ):
            if value < 0.0:
                raise DomainError(f"{label} must be finite and non-negative.")
        if vehicle_mass_kg <= 0.0:
            raise DomainError("Vehicle mass must be greater than zero.")
        states = _validated_states(additional_states)
        bound_states = sorted({binding.state_name for binding in self.stream_bindings})
        missing = set(bound_states) - set(states)
        if missing:
            raise DomainError("Missing bound propulsion states: " + ", ".join(sorted(missing)))
        return PropellantMassClosure(
            vehicle_mass_kg, non_propellant_mass_kg,
            *_inventory_ledger(vehicle_mass_kg, non_propellant_mass_kg,
                               (states[name] for name in bound_states)),
        )

    def time_event_jump(
        self,
        *,
        relative_time_s: float,
        vehicle_mass_kg: float,
        direction: tuple[float, float, float],
        additional_states: Mapping[str, float] | None = None,
        operating_conditions: Mapping[str, float] | None = None,
        parameter_overrides: Mapping[str, float] | None = None,
        ambient_pressure_pa: float | None = None,
    ) -> PropulsionDerivativeJump:
        """Return ``f(left) - f(right)`` for the saltation update at a time event.

        Parameters
        ----------
        relative_time_s : float
            Event time relative to nominal ignition, in s.
        vehicle_mass_kg : float
            Finite, positive current vehicle mass, in kg.
        direction : tuple of float
            Non-zero thrust direction in the propagator's resolved frame; it is
            normalized internally.
        additional_states : Mapping of str to float, optional
            Current value of every bound additional state.
        operating_conditions : Mapping of str to float, optional
            Performance-surface conditions; only allowed with a surface.
        parameter_overrides : Mapping of str to float, optional
            Values of ``thrust_scale``, ``mass_flow_scale`` and
            ``ignition_time_bias_s`` that replace the reference values.
        ambient_pressure_pa : float or None, optional
            Non-negative ambient pressure, in Pa, as for :meth:`evaluate`.

        Returns
        -------
        PropulsionDerivativeJump
            Left-minus-right acceleration and drains.

        Raises
        ------
        DomainError
            If an input is invalid (see :meth:`evaluate`).

        Notes
        -----
        Inventory limits are not applied, so the jump is the profile's own
        discontinuity.  Its use in sensitivity propagation is described in
        :class:`PropulsionDerivativeJump` and [1]_.

        References
        ----------
        .. [1] NASA-20180003657. NASA, "Navigation Filter Best Practices"
               (NASA/TP-2018-219822), https://ntrs.nasa.gov/citations/20180003657.
        """

        common = {
            "relative_time_s": relative_time_s,
            "vehicle_mass_kg": vehicle_mass_kg,
            "direction": direction,
            "additional_states": additional_states,
            "operating_conditions": operating_conditions,
            "parameter_overrides": parameter_overrides,
            "ambient_pressure_pa": ambient_pressure_pa,
            "apply_inventory_limits": False,
        }
        left = self.evaluate(side="left", **common)  # type: ignore[arg-type]
        right = self.evaluate(side="right", **common)  # type: ignore[arg-type]
        return self._jump_from_evaluations(left, right)

    def time_event_jump_at_epoch(self, *, time_s: float, ignition_time_s: float,
                                  **kwargs) -> PropulsionDerivativeJump:
        """Return the saltation jump at a represented canonical epoch event.

        Parameters
        ----------
        time_s : float
            Absolute propagation time of the event, in s.
        ignition_time_s : float
            Absolute nominal ignition time, in s.
        **kwargs
            Remaining keyword arguments of :meth:`evaluate_at_epoch`, except
            ``side`` and ``apply_inventory_limits``.

        Returns
        -------
        PropulsionDerivativeJump
            Left-minus-right acceleration and drains.

        Raises
        ------
        DomainError
            If an input is invalid (see :meth:`evaluate_at_epoch`).

        Notes
        -----
        Inventory limits are excluded, as for :meth:`time_event_jump`.  The
        sensitivity update is described in :class:`PropulsionDerivativeJump`
        and [1]_; the absolute clock is resolved exactly [2]_.

        References
        ----------
        .. [1] NASA-20180003657. NASA, "Navigation Filter Best Practices"
               (NASA/TP-2018-219822), https://ntrs.nasa.gov/citations/20180003657.
        .. [2] PYTHON-3.12-FRACTION. Python 3.12, "fractions": float inputs retain
               their exact binary value as a rational,
               https://docs.python.org/3.12/library/fractions.html.
        """
        kwargs = dict(kwargs)
        kwargs['apply_inventory_limits'] = False
        left = self.evaluate_at_epoch(time_s=time_s, ignition_time_s=ignition_time_s,
                                      side='left', **kwargs)
        right = self.evaluate_at_epoch(time_s=time_s, ignition_time_s=ignition_time_s,
                                       side='right', **kwargs)
        return self._jump_from_evaluations(left, right)

    @staticmethod
    def _jump_from_evaluations(left, right) -> PropulsionDerivativeJump:
        """Subtract two evaluations; shared by the relative and epoch jump methods."""
        right_states = dict(right.additional_state_derivatives)
        return PropulsionDerivativeJump(
            relative_time_s=left.relative_time_s,
            acceleration_m_s2=tuple(
                a - b for a, b in zip(left.acceleration_m_s2, right.acceleration_m_s2, strict=True)
            ),  # type: ignore[arg-type]
            mass_derivative_kg_s=left.mass_derivative_kg_s - right.mass_derivative_kg_s,
            additional_state_derivatives=tuple(
                (name, rate - right_states.get(name, 0.0))
                for name, rate in left.additional_state_derivatives
            ),
        )

    def operating_condition_grid_surfaces(self) -> tuple[PerformanceGridLineSurface, ...]:
        """Return non-stopping derivative-kink roots of the performance surface.

        Returns
        -------
        tuple of PerformanceGridLineSurface
            One root per interior grid line; empty without a surface.

        Notes
        -----
        A consumer registers these next to :meth:`operating_condition_surfaces`
        and restarts integration on each crossing, as for profile boundaries
        [1]_.

        References
        ----------
        .. [1] OREKIT-13.1.5-PROPULSION-MODEL. Orekit 13.1.5, "PropulsionModel /
               EventDetector" API: force, mass-derivative, parameter-driver and
               event contract.
        """

        if self.performance_surface is None and self.stream_performance_surface is None:
            return ()
        if self.performance_surface is not None:
            return self.performance_surface.grid_line_surfaces()
        assert self.stream_performance_surface is not None
        return self.stream_performance_surface.grid_line_surfaces()

    def _segment_at(
        self, effective_time_s: float, side: str
    ) -> tuple[TabulatedBurnSegment | None, bool]:
        """Return the segment owning an exact profile time on ``side`` and a boundary flag."""
        artifact = self.artifact
        effective_time_s = _fraction(effective_time_s)
        at_boundary = artifact._is_boundary_time(effective_time_s)
        if effective_time_s < 0.0 or effective_time_s > artifact.duration_s:
            return None, at_boundary
        if effective_time_s == 0.0 and side == "left":
            return None, True
        if effective_time_s == artifact.duration_s and side == "right":
            return None, True
        index = artifact._segment_index_at(effective_time_s, side=side)
        segments = artifact.segments
        if (
            side == "left"
            and index + 1 < len(segments)
            and effective_time_s == segments[index + 1].start_time_s
        ):
            at_boundary = True
        return segments[index], at_boundary

    def evaluate(self, *, relative_time_s: float, vehicle_mass_kg: float,
                 direction: tuple[float, float, float],
                 additional_states: Mapping[str, float] | None = None,
                 operating_conditions: Mapping[str, float] | None = None,
                 parameter_overrides: Mapping[str, float] | None = None,
                 side: str = 'right', apply_inventory_limits: bool = True,
                 ambient_pressure_pa: float | None = None,
                 boundary_tolerance_s: float | None = None) -> PropulsionDynamicsEvaluation:
        """Evaluate force, drains, event-safe activity and analytic partials.

        Parameters
        ----------
        relative_time_s : float
            Finite propagation time relative to nominal ignition, in s.
        vehicle_mass_kg : float
            Finite, positive current vehicle mass, in kg.
        direction : tuple of float
            Non-zero thrust direction in the propagator's resolved frame; it is
            normalized internally.
        additional_states : Mapping of str to float, optional
            Current value of every bound additional state.
        operating_conditions : Mapping of str to float, optional
            Performance-surface conditions; only allowed with a surface.
        parameter_overrides : Mapping of str to float, optional
            Values of ``thrust_scale``, ``mass_flow_scale`` and
            ``ignition_time_bias_s`` that replace the reference values.
        side : {'right', 'left'}, optional
            Limit taken exactly at a discontinuity.
        apply_inventory_limits : bool, optional
            Zero the output at or below a dry-mass or reserve limit.
        ambient_pressure_pa : float or None, optional
            Non-negative ambient pressure, in Pa; required while firing when
            ``nozzle_exit_area_m2`` is set and refused otherwise.
        boundary_tolerance_s : float or None, optional
            Caller-owned window within which a time snaps to a profile
            boundary; ``None`` uses exact relative coordinates.

        Returns
        -------
        PropulsionDynamicsEvaluation
            Force, acceleration, drains and partials; exact zeros outside the
            shifted firing interval or when an inventory limit inhibits the
            burn.

        Raises
        ------
        DomainError
            If an input is invalid, a bound state is missing, conditions are
            supplied without a surface, or ambient pressure is missing while
            firing with ``nozzle_exit_area_m2`` (or supplied without it).

        Notes
        -----
        ``side`` selects the limiting value exactly at discontinuities: normal
        integration uses ``'right'`` after an event reset; event localization
        may query both sides.  Only exact boundary coordinates select a
        one-sided limit.  Outside the shifted firing interval the method
        returns exact zeros rather than extrapolating the profile.

        ``apply_inventory_limits=True`` zeroes the output when the vehicle
        mass or a constrained state is at or below its limit plus
        ``INVENTORY_LIMIT_TOLERANCE``.  This is a safety net for consumers
        without event handling.  It makes the right-hand side discontinuous
        *inside* an integration step: a Runge-Kutta stage that probes just
        past a reserve root sees zero drain, so the step under-drains and a
        root located on that step is biased by a fraction of the step.  A
        consumer that registers the dry-mass and reserve
        :meth:`event_surfaces`, locates their roots and switches the burn off
        itself passes ``False``.  The profile is never extrapolated in time
        either way.

        With ``nozzle_exit_area_m2`` set, the artifact's axial thrust is the
        *vacuum* thrust and ``F = F_vac - p_a A_e``, ``dF/dp_a = -A_e`` [2]_.
        The net thrust is not clamped; ``net_thrust_margin_n`` is the root a
        consumer monitors to stop the burn before the net force reverses.  It
        is defined only while the burn is firing (``None`` otherwise), so the
        consumer enables that root between the burn's start and stop events.

        Absolute-clock consumers use :meth:`evaluate_at_epoch`;
        :class:`PropulsionTimeline` does so.  Shifted coordinates keep their
        exact rational position [3]_ through segment ownership and
        interpolation, even when the displayed time rounds onto a boundary.
        Numerical clock rules do not establish sensor accuracy.  The
        evaluation contract follows [1]_.

        References
        ----------
        .. [1] OREKIT-13.1.5-PROPULSION-MODEL. Orekit 13.1.5, "PropulsionModel /
               EventDetector" API: force, mass-derivative, parameter-driver and
               event contract.
        .. [2] NASA-STEADY-ROCKET-THRUST. NASA Glenn Research Center, momentum and
               pressure thrust equation,
               https://www1.grc.nasa.gov/beginners-guide-to-aeronautics/rocket-thrust/.
        .. [3] PYTHON-3.12-FRACTION. Python 3.12, "fractions": float inputs retain
               their exact binary value as a rational,
               https://docs.python.org/3.12/library/fractions.html.
        """
        arguments = locals().copy()
        arguments.pop('self')
        return self._evaluate(**arguments)

    def evaluate_at_epoch(self, *, time_s: float, ignition_time_s: float,
                          vehicle_mass_kg: float, direction: tuple[float,float,float],
                          **kwargs) -> PropulsionDynamicsEvaluation:
        """Evaluate at an absolute binary64 clock reading without a snap window.

        Parameters
        ----------
        time_s : float
            Finite absolute propagation time, in s.
        ignition_time_s : float
            Finite absolute nominal ignition time, in s.
        vehicle_mass_kg : float
            Finite, positive current vehicle mass, in kg.
        direction : tuple of float
            Non-zero thrust direction in the resolved frame.
        **kwargs
            Remaining keyword arguments of :meth:`evaluate` other than
            ``relative_time_s``.

        Returns
        -------
        PropulsionDynamicsEvaluation
            As for :meth:`evaluate`, with ``relative_time_s`` the exact
            difference ``time_s - ignition_time_s`` rounded once.

        Raises
        ------
        DomainError
            If an input is invalid, ``relative_time_s`` is supplied, or the
            absolute clock cannot distinguish successive profile boundaries.

        Notes
        -----
        Each profile boundary maps to the once-rounded absolute date
        ``ignition + bias + t_k``; distinct boundaries require distinct
        dates.  A reading equal to such a date, or to an unambiguous
        addition-order alias of it, selects that boundary's one-sided limit;
        any other reading maps to its exact profile time [1]_.

        References
        ----------
        .. [1] PYTHON-3.12-FRACTION. Python 3.12, "fractions": float inputs retain
               their exact binary value as a rational,
               https://docs.python.org/3.12/library/fractions.html.
        """
        time_s = real(time_s, 'Absolute propagation time')
        if '_resolved_profile_time' in kwargs or 'relative_time_s' in kwargs:
            raise DomainError('Epoch evaluation resolves its own profile coordinate.')
        ignition_time_s = real(ignition_time_s, 'Ignition time')
        parameters = self.resolved_parameters(kwargs.get('parameter_overrides'))
        effective = _epoch_profile_time(self.artifact.segments, time_s, ignition_time_s,
                                        parameters['ignition_time_bias_s'])
        relative = _finite_rational(_fraction(time_s)-_fraction(ignition_time_s),
                                    'Relative propagation time')
        return self._evaluate(relative_time_s=relative, vehicle_mass_kg=vehicle_mass_kg,
                              direction=direction, _resolved_profile_time=effective, **kwargs)

    def _evaluate(
        self,
        *,
        relative_time_s: float,
        vehicle_mass_kg: float,
        direction: tuple[float, float, float],
        additional_states: Mapping[str, float] | None = None,
        operating_conditions: Mapping[str, float] | None = None,
        parameter_overrides: Mapping[str, float] | None = None,
        side: str = "right",
        apply_inventory_limits: bool = True,
        ambient_pressure_pa: float | None = None,
        boundary_tolerance_s: float | None = None,
        _resolved_profile_time: Fraction | None = None,
    ) -> PropulsionDynamicsEvaluation:
        """Shared implementation of :meth:`evaluate` and :meth:`evaluate_at_epoch`.

        ``_resolved_profile_time`` is the exact profile coordinate already
        resolved from an absolute clock; when ``None`` it is
        ``relative_time_s - ignition_time_bias_s``.  The axial thrust is
        ``F = F_vac - p_a A_e`` when a nozzle exit area is declared [2]_, and
        the tabulated transient history [3]_ is never extrapolated.  The full
        contract is documented on :meth:`evaluate` and follows [1]_.

        References
        ----------
        .. [1] OREKIT-13.1.5-PROPULSION-MODEL. Orekit 13.1.5, "PropulsionModel /
               EventDetector" API: force, mass-derivative, parameter-driver and
               event contract.
        .. [2] NASA-STEADY-ROCKET-THRUST. NASA Glenn Research Center, momentum and
               pressure thrust equation,
               https://www1.grc.nasa.gov/beginners-guide-to-aeronautics/rocket-thrust/.
        .. [3] NASA-20040000363. NASA, "Transient Mathematical Modeling for Liquid
               Rocket Engine Systems", https://ntrs.nasa.gov/citations/20040000363.
        """

        if side not in ("left", "right"):
            raise DomainError("Propagation boundary side must be 'left' or 'right'.")
        relative_time_s = real(relative_time_s, "Relative propagation time")
        vehicle_mass_kg = real(vehicle_mass_kg, "Vehicle mass")
        if not isinstance(apply_inventory_limits, bool):
            raise DomainError("apply_inventory_limits must be a boolean.")
        if not isfinite(relative_time_s):
            raise DomainError("Relative propagation time must be finite.")
        if not isfinite(vehicle_mass_kg) or vehicle_mass_kg <= 0.0:
            raise DomainError("Vehicle mass must be finite and greater than zero.")
        geometry = _direction_geometry(direction)
        states = _validated_states(additional_states)
        required_states = {binding.state_name for binding in self.stream_bindings}
        missing_states = required_states - set(states)
        if missing_states:
            raise DomainError(
                "Missing bound propulsion states: " + ", ".join(sorted(missing_states))
            )
        if operating_conditions is not None and not isinstance(operating_conditions, Mapping):
            raise DomainError("Operating conditions must be a mapping.")
        supplied_conditions = dict(operating_conditions or {})
        if (
            self.performance_surface is None
            and self.stream_performance_surface is None
            and supplied_conditions
        ):
            raise DomainError(
                "Operating conditions were supplied but this bridge has no performance surface."
            )
        if ambient_pressure_pa is not None:
            if self.nozzle_exit_area_m2 is None:
                raise DomainError(
                    "Ambient pressure was supplied but this bridge has no nozzle exit area."
                )
            ambient_pressure_pa = real(ambient_pressure_pa, "Ambient pressure")
            if ambient_pressure_pa < 0.0:
                raise DomainError("Ambient pressure must be finite and non-negative.")

        parameter_values = self.resolved_parameters(parameter_overrides)
        thrust_scale = parameter_values["thrust_scale"]
        flow_scale = parameter_values["mass_flow_scale"]
        time_bias = parameter_values["ignition_time_bias_s"]
        exact_time = (_fraction(relative_time_s) - _fraction(time_bias)
                      if _resolved_profile_time is None else _resolved_profile_time)
        exact_time = self.artifact._snap_to_boundary(exact_time, boundary_tolerance_s)
        segment, at_boundary = self._segment_at(exact_time, side)
        effective_time = _finite_rational(exact_time, 'Effective profile time')

        inhibited: list[str] = []
        if apply_inventory_limits:
            if self.protected_dry_mass_kg is not None and _fraction(vehicle_mass_kg) <= (
                _fraction(self.protected_dry_mass_kg) + _fraction(INVENTORY_LIMIT_TOLERANCE)
            ):
                inhibited.append("protected_dry_mass_kg")
            for constraint in self.state_constraints:
                if (
                    _fraction(states[constraint.state_name])
                    <= _fraction(constraint.minimum_value) + _fraction(INVENTORY_LIMIT_TOLERANCE)
                ):
                    inhibited.append(constraint.state_name)
        active = segment is not None and not inhibited

        if not active:
            zero = (0.0, 0.0, 0.0)
            return PropulsionDynamicsEvaluation(
                relative_time_s=relative_time_s,
                effective_profile_time_s=effective_time,
                active=False,
                phase="inhibited" if inhibited else "off",
                axial_thrust_n=0.0,
                force_vector_n=zero,
                acceleration_m_s2=zero,
                mass_derivative_kg_s=0.0,
                additional_state_derivatives=tuple((name, 0.0) for name in sorted(required_states)),
                parameters=tuple(parameter_values.items()),
                inhibited_by=tuple(inhibited),
                partials=PropulsionPartials(
                    acceleration_wrt_mass=zero,
                    acceleration_wrt_parameters=tuple((name, zero) for name in parameter_values),
                    mass_rate_wrt_parameters=tuple((name, 0.0) for name in parameter_values),
                    additional_state_rates_wrt_parameters=tuple(
                        (name, tuple((state, 0.0) for state in sorted(required_states)))
                        for name in parameter_values
                    ),
                    nonsmooth_parameters=("ignition_time_bias_s",) if at_boundary else (),
                ),
            )

        assert segment is not None
        surface_evaluation = (
            None
            if self.performance_surface is None
            else self.performance_surface.evaluate(supplied_conditions)
        )
        stream_surface_evaluation = (
            None
            if self.stream_performance_surface is None
            else self.stream_performance_surface.evaluate(supplied_conditions)
        )
        surface_force_scale = (
            surface_evaluation.force_scale
            if surface_evaluation is not None
            else (
                stream_surface_evaluation.force_scale
                if stream_surface_evaluation is not None
                else 1.0
            )
        )
        base_thrust = segment._axial_at(
            min(max(exact_time, segment.start_time_s), segment.end_time_s)
        )
        base_flow, stream_flows, flow_slope, stream_slopes = _mass_profile(
            segment, min(max(exact_time, segment.start_time_s), segment.end_time_s),
            self.artifact.interpolation is InterpolationPolicy.HOLD)
        flow_driver = _fraction(flow_scale)
        common_surface_flow_scale = (
            _fraction(1.0 if surface_evaluation is None else surface_evaluation.mass_flow_scale)
        )
        stream_surface_scales = (
            {}
            if stream_surface_evaluation is None
            else dict(stream_surface_evaluation.stream_flow_scales)
        )
        stream_scale_by_id = {
            stream_id: _fraction(stream_surface_scales.get(stream_id, common_surface_flow_scale))
            for stream_id in stream_flows
        }
        surface_scaled_total_flow = (
            sum((flow * stream_scale_by_id[stream_id] for stream_id, flow in stream_flows.items()),
                Fraction(0))
            if stream_flows
            else base_flow * common_surface_flow_scale
        )
        surface_flow_scale = (_finite_rational(surface_scaled_total_flow / base_flow,
                                               "Performance mass-flow scale")
                              if base_flow > 0 else 1.0)
        vacuum_thrust = _fraction(base_thrust) * _fraction(thrust_scale) * _fraction(surface_force_scale)
        pressure_loss = 0.0
        if self.nozzle_exit_area_m2 is not None:
            if ambient_pressure_pa is None:
                raise DomainError(
                    "This bridge has a nozzle exit area; ambient pressure is required "
                    "while the burn is active."
                )
            pressure_loss = _fraction(ambient_pressure_pa) * _fraction(self.nozzle_exit_area_m2)
        # Kept on the linear branch even when negative: a consumer stops the burn
        # at the net_thrust_margin_n root instead of seeing a clamp inside a step.
        exact_thrust = vacuum_thrust - pressure_loss
        thrust = _finite_rational(exact_thrust, "Axial thrust")
        pressure_loss = _finite_rational(_fraction(pressure_loss), "Pressure thrust loss")
        total_flow = _finite_rational(surface_scaled_total_flow * flow_driver, "Total mass flow")
        force = geometry.project(exact_thrust, label="Force")
        acceleration = geometry.project(exact_thrust, vehicle_mass_kg, label="Acceleration")

        binding_by_stream = {
            binding.stream_id: binding.state_name for binding in self.stream_bindings
        }
        state_drains = {state_name: Fraction(0) for state_name in required_states}
        for stream_id, flow in stream_flows.items():
            if stream_id in binding_by_stream:
                state_drains[binding_by_stream[stream_id]] -= (
                    flow * flow_driver * stream_scale_by_id[stream_id]
                )

        duration = segment.duration_s
        if self.artifact.interpolation is InterpolationPolicy.HOLD:
            thrust_slope = 0.0
        else:
            thrust_slope = (segment.end_axial_thrust_n - segment.start_axial_thrust_n) / duration
        surface_scaled_flow_slope = (
            sum((slope * stream_scale_by_id[stream_id] for stream_id, slope in stream_slopes.items()),
                Fraction(0))
            if stream_slopes
            else flow_slope * common_surface_flow_scale
        )

        zero = (0.0, 0.0, 0.0)
        thrust_scale_partial = geometry.project(
            _fraction(base_thrust) * _fraction(surface_force_scale), vehicle_mass_kg,
            label="Thrust scale acceleration partial")
        time_partial = geometry.project(
            -_fraction(thrust_scale) * _fraction(surface_force_scale) * _fraction(thrust_slope),
            vehicle_mass_kg, label="Ignition time acceleration partial")
        mass_partials = (
            ("thrust_scale", 0.0),
            ("mass_flow_scale", _finite_rational(-surface_scaled_total_flow,
                                                "Total mass-flow scale partial")),
            (
                "ignition_time_bias_s",
                _finite_rational(flow_driver * surface_scaled_flow_slope,
                                 "Total ignition-time flow partial"),
            ),
        )

        def state_parameter_partials(parameter_name: str) -> tuple[tuple[str, float], ...]:
            values = {state_name: Fraction(0) for state_name in required_states}
            if parameter_name == "mass_flow_scale":
                for stream_id, flow in stream_flows.items():
                    if stream_id in binding_by_stream:
                        values[binding_by_stream[stream_id]] -= flow * stream_scale_by_id[stream_id]
            elif parameter_name == "ignition_time_bias_s":
                for stream_id, slope in stream_slopes.items():
                    if stream_id in binding_by_stream:
                        values[binding_by_stream[stream_id]] += (
                            flow_driver * stream_scale_by_id[stream_id] * slope
                        )
            return tuple((name, _finite_rational(value, f"Tank {parameter_name} partial: {name}"))
                         for name, value in sorted(values.items()))

        force_condition_gradients = (
            dict(surface_evaluation.force_scale_gradient)
            if surface_evaluation is not None
            else (
                dict(stream_surface_evaluation.force_scale_gradient)
                if stream_surface_evaluation is not None
                else {}
            )
        )
        common_flow_condition_gradients = (
            dict(surface_evaluation.mass_flow_scale_gradient)
            if surface_evaluation is not None
            else {}
        )
        stream_flow_condition_gradients = (
            {}
            if stream_surface_evaluation is None
            else {
                stream_id: dict(gradients)
                for stream_id, gradients in (stream_surface_evaluation.stream_flow_scale_gradients)
            }
        )
        condition_names = tuple(force_condition_gradients)
        total_flow_condition_derivatives = {
            name: (
                sum((flow * _fraction(common_flow_condition_gradients[name]
                     if surface_evaluation is not None else stream_flow_condition_gradients[stream_id][name])
                     for stream_id, flow in stream_flows.items()), Fraction(0))
                if stream_flows else base_flow * _fraction(common_flow_condition_gradients[name])
            )
            for name in condition_names
        }
        stream_flow_condition_derivative = {
            (stream_id, name): (
                _fraction(common_flow_condition_gradients[name]
                          if surface_evaluation is not None
                          else stream_flow_condition_gradients[stream_id][name])
            )
            for stream_id in stream_flows
            for name in condition_names
        }
        acceleration_condition_partials = tuple(
            (
                name,
                geometry.project(_fraction(base_thrust) * _fraction(thrust_scale) * _fraction(gradient),
                                 vehicle_mass_kg, label=f"Condition acceleration partial: {name}"),
            )
            for name, gradient in force_condition_gradients.items()
        )
        mass_condition_partials = tuple(
            (name, _finite_rational(-flow_driver * total_flow_condition_derivatives[name],
                                   f"Total mass-flow condition partial: {name}"))
            for name in condition_names
        )

        def state_condition_partials(condition_name: str) -> tuple[tuple[str, float], ...]:
            values = {state_name: Fraction(0) for state_name in required_states}
            for stream_id, flow in stream_flows.items():
                if stream_id in binding_by_stream:
                    values[binding_by_stream[stream_id]] -= (
                        flow
                        * flow_driver
                        * stream_flow_condition_derivative[(stream_id, condition_name)]
                    )
            return tuple((name, _finite_rational(value, f"Tank condition partial: {name}"))
                         for name, value in sorted(values.items()))

        return PropulsionDynamicsEvaluation(
            relative_time_s=relative_time_s,
            effective_profile_time_s=effective_time,
            active=True,
            phase=segment.phase,
            axial_thrust_n=thrust,
            force_vector_n=force,  # type: ignore[arg-type]
            acceleration_m_s2=acceleration,  # type: ignore[arg-type]
            mass_derivative_kg_s=-total_flow,
            additional_state_derivatives=tuple(
                (name, _finite_rational(value, f"Tank mass drain: {name}"))
                for name, value in sorted(state_drains.items())),
            parameters=tuple(parameter_values.items()),
            inhibited_by=(),
            partials=PropulsionPartials(
                acceleration_wrt_mass=geometry.project(
                    -exact_thrust, _fraction(vehicle_mass_kg)**2, label="Mass acceleration partial"),
                acceleration_wrt_parameters=(
                    ("thrust_scale", thrust_scale_partial),  # type: ignore[arg-type]
                    ("mass_flow_scale", zero),
                    ("ignition_time_bias_s", time_partial),  # type: ignore[arg-type]
                ),
                mass_rate_wrt_parameters=mass_partials,
                additional_state_rates_wrt_parameters=tuple(
                    (name, state_parameter_partials(name)) for name in parameter_values
                ),
                acceleration_wrt_conditions=acceleration_condition_partials,
                mass_rate_wrt_conditions=mass_condition_partials,
                additional_state_rates_wrt_conditions=tuple(
                    (name, state_condition_partials(name)) for name in condition_names
                ),
                nonsmooth_parameters=("ignition_time_bias_s",) if at_boundary else (),
                acceleration_wrt_direction=geometry.jacobian(exact_thrust, vehicle_mass_kg),
                acceleration_wrt_ambient_pressure=(
                    (0.0, 0.0, 0.0)
                    if self.nozzle_exit_area_m2 is None
                    else geometry.project(-self.nozzle_exit_area_m2, vehicle_mass_kg,
                                          label="Ambient pressure acceleration partial")
                ),  # type: ignore[arg-type]
            ),
            operating_conditions=(
                surface_evaluation.conditions
                if surface_evaluation is not None
                else (
                    stream_surface_evaluation.conditions
                    if stream_surface_evaluation is not None
                    else ()
                )
            ),
            performance_force_scale=surface_force_scale,
            performance_mass_flow_scale=surface_flow_scale,
            ambient_pressure_pa=ambient_pressure_pa,
            pressure_thrust_loss_n=pressure_loss,
            net_thrust_margin_n=(None if self.nozzle_exit_area_m2 is None else thrust),
        )
