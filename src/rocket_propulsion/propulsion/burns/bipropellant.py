"""Dual-feed bipropellant coupling for numerical orbit propagation.

The model joins two regulated tanks, two inertance/compliance liquid lines
[2]_, [3]_, a stream-resolved engine map and the propagation bridge [4]_.
Oxidizer and fuel remain independent mass channels, so mixture-ratio drift
under separate feed pressures [1]_, first-component depletion and
cross-feed-pressure sensitivities are never hidden inside one equivalent
specific impulse.

References
----------
.. [1] NASA, "Simulation Analysis of Computer-Controlled Pressurization for
       Mixture Ratio Control". https://ntrs.nasa.gov/citations/20050206420
.. [2] NASA, "Nodal Modeling of Liquid Propellant Feed and Pressurization
       System". https://ntrs.nasa.gov/citations/20240003493
.. [3] NASA, "Liquid Rocket Propellant Feedline Dynamics".
       https://ntrs.nasa.gov/citations/19740028545
.. [4] Orekit 13.1.5, "PropulsionModel / EventDetector API", API documentation of
       org.orekit.forces.maneuvers.propulsion.PropulsionModel.
"""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass
from enum import Enum
from fractions import Fraction
from math import fsum, inf, isfinite

from rocket_propulsion.core.errors import DomainError
from rocket_propulsion.core.validation import (
    coerce_fields,
    name_value_pairs,
    real,
    real_mapping,
    real_matrix,
    real_tuple,
    records,
    strings,
    text,
)

from .feed_line import (
    DynamicFeedLine,
    DynamicFeedLineEvaluation,
    DynamicFeedLineState,
)
from .feed_system import (
    RegulatedFeedEvaluation,
    RegulatedFeedState,
    RegulatedFeedSystem,
)
from .numeric_outputs import _finite_rational
from .propagator import PropagatorPropulsionBridge, PropulsionDynamicsEvaluation
from .serialization import canonical_sha256

BIPROPELLANT_DYNAMIC_FEED_SCHEMA = "rocket_propulsion_bipropellant_dynamic_feed_v1"
BIPROPELLANT_DYNAMIC_FEED_REFERENCE_IDS = (
    "NASA-20050206420",
    "NASA-19730002016",
    "NASA-20240003493",
    "NASA-19740028545",
    "OREKIT-13.1.5-PROPULSION-MODEL",
)


def _finite_quotient(numerator: float, denominator: float, label: str) -> float:
    """Return ``numerator / denominator`` rounded once from the exact rational quotient.

    Raises
    ------
    DomainError
        If either input is non-finite, the denominator is not positive, or the
        quotient is not representable as a finite float.
    """
    top = real(numerator, label + " numerator")
    bottom = real(denominator, label + " denominator")
    if bottom <= 0.0:
        raise DomainError(f"{label} denominator must be positive.", field=label)
    return _finite_rational(Fraction(top) / Fraction(bottom), label)


def _mixture_ratio_partial(oxidizer: float, fuel: float,
                           d_oxidizer: float, d_fuel: float) -> float:
    """Return the sensitivity of the mixture ratio ``O/F`` to one condition.

    With oxidizer flow ``O``, fuel flow ``F`` and their partials ``dO`` and
    ``dF``, the quotient rule gives ``d(O/F) = (dO * F - O * dF) / F**2``,
    where ``O/F`` is the oxidizer-to-fuel mass-flow ratio [1]_.  Products and
    differences are formed exactly on the supplied binary values and rounded
    once at output; this bounds numerical loss and says nothing about physical
    uncertainty.

    Raises
    ------
    DomainError
        If an input is non-finite, the fuel flow is not positive, or the
        result is not a finite float.

    References
    ----------
    .. [1] SUTTON-BIBLARZ-2017. G. P. Sutton and O. Biblarz, "Rocket
           Propulsion Elements", 9th ed., Wiley, 2017.
    """
    ox = Fraction(real(oxidizer, "Oxidizer engine flow"))
    fu = Fraction(real(fuel, "Fuel engine flow"))
    dox = Fraction(real(d_oxidizer, "Oxidizer flow partial"))
    dfu = Fraction(real(d_fuel, "Fuel flow partial"))
    if fu <= 0:
        raise DomainError("Mixture-ratio partial requires positive fuel flow.")
    exact = (dox * fu - ox * dfu) / (fu * fu)
    return _finite_rational(exact, "Engine mixture-ratio partial")


def _finite_columns(value, size: int, label: str):
    """Return ``(name, vector)`` pairs of finite ``size``-vectors with unique names."""
    columns = []
    for item in records(value, (list, tuple), label):
        if len(item) != 2:
            raise DomainError(f"{label} requires named vectors.")
        name = text(item[0], label + " name")
        vector = real_tuple(item[1], label)
        if len(vector) != size:
            raise DomainError(f"{label} vectors require {size} components.")
        columns.append((name, vector))
    if len({name for name, _ in columns}) != len(columns):
        raise DomainError(f"{label} names must be unique.")
    return tuple(columns)


def _finite_pairs(value, label: str):
    """Return ``(name, value)`` pairs of finite scalars with unique names."""
    pairs = name_value_pairs(value, label)
    if len({name for name, _ in pairs}) != len(pairs):
        raise DomainError(f"{label} names must be unique.")
    return pairs


@dataclass(frozen=True, slots=True)
class BipropellantDynamicFeedConfiguration:
    """State names, mixture-ratio limits and closure tolerance of two feed legs.

    Parameters
    ----------
    oxidizer_state_name, fuel_state_name : str
        Distinct bridge additional-mass states drained by every stream
        assigned to the oxidizer and fuel legs.
    oxidizer_pressure_condition_name, fuel_pressure_condition_name : str
        Distinct operating-condition names of the two manifold pressures in
        the stream-resolved engine map.
    minimum_engine_mixture_ratio, maximum_engine_mixture_ratio : float
        Finite, positive, strictly ordered engine mixture-ratio limits.
    mass_closure_absolute_tolerance_kg_s : float, default 1.0e-10
        Finite, positive absolute tolerance of the mass-closure checks, kg/s.
    model_id : str, default "dual-regulated-dynamic-feed"
        Non-blank model identifier.
    schema : str
        Must equal ``BIPROPELLANT_DYNAMIC_FEED_SCHEMA``.

    Raises
    ------
    DomainError
        If the schema is unsupported, a name is blank, the two legs share a
        state or pressure-condition name, the limits are not finite, positive
        and ordered, or the tolerance is not finite and positive.

    Notes
    -----
    The engine mixture ratio is oxidizer engine flow divided by fuel engine
    flow.  It drifts when the two feeds are pressurized independently [1]_, and
    short pulses show transient flow-ratio excursions [2]_; the limits bound
    the operating region in which the engine map is used.

    References
    ----------
    .. [1] NASA, "Simulation Analysis of Computer-Controlled Pressurization for
           Mixture Ratio Control". https://ntrs.nasa.gov/citations/20050206420
    .. [2] NASA, "Space Simulation Experiments on Reaction Control System Thruster
           Plumes". https://ntrs.nasa.gov/citations/19730002016
    """

    oxidizer_state_name: str
    fuel_state_name: str
    oxidizer_pressure_condition_name: str
    fuel_pressure_condition_name: str
    minimum_engine_mixture_ratio: float
    maximum_engine_mixture_ratio: float
    mass_closure_absolute_tolerance_kg_s: float = 1.0e-10
    model_id: str = "dual-regulated-dynamic-feed"
    schema: str = BIPROPELLANT_DYNAMIC_FEED_SCHEMA

    def __post_init__(self) -> None:
        coerce_fields(self, "Bipropellant configuration")
        if self.schema != BIPROPELLANT_DYNAMIC_FEED_SCHEMA:
            raise DomainError("Unsupported bipropellant dynamic-feed schema.")
        names = (
            self.oxidizer_state_name,
            self.fuel_state_name,
            self.oxidizer_pressure_condition_name,
            self.fuel_pressure_condition_name,
            self.model_id,
        )
        if any(not name.strip() for name in names):
            raise DomainError("Bipropellant state, condition, and model names must not be empty.")
        if self.oxidizer_state_name == self.fuel_state_name:
            raise DomainError("Oxidizer and fuel must use distinct mass states.")
        if self.oxidizer_pressure_condition_name == self.fuel_pressure_condition_name:
            raise DomainError("Oxidizer and fuel must use distinct pressure conditions.")
        if (
            not isfinite(self.minimum_engine_mixture_ratio)
            or not isfinite(self.maximum_engine_mixture_ratio)
            or self.minimum_engine_mixture_ratio <= 0.0
            or self.maximum_engine_mixture_ratio <= self.minimum_engine_mixture_ratio
        ):
            raise DomainError(
                "Bipropellant mixture-ratio limits must be finite, positive, and ordered."
            )
        if (
            not isfinite(self.mass_closure_absolute_tolerance_kg_s)
            or self.mass_closure_absolute_tolerance_kg_s <= 0.0
        ):
            raise DomainError("Bipropellant closure tolerance must be finite and positive.")


@dataclass(frozen=True, slots=True)
class BipropellantDynamicFeedState:
    """Two regulated-feed states and two liquid-line states.

    Parameters
    ----------
    oxidizer_feed, fuel_feed : RegulatedFeedState
        Tank and pressurant states of each leg.
    oxidizer_line, fuel_line : DynamicFeedLineState
        Line mass-flow and manifold-pressure states of each leg.

    Raises
    ------
    DomainError
        If a field has the wrong type.

    Notes
    -----
    The states are the nodal mass and pressure variables of a coupled feed
    network [1]_.

    References
    ----------
    .. [1] NASA, "Nodal Modeling of Liquid Propellant Feed and Pressurization
           System". https://ntrs.nasa.gov/citations/20240003493
    """

    oxidizer_feed: RegulatedFeedState
    oxidizer_line: DynamicFeedLineState
    fuel_feed: RegulatedFeedState
    fuel_line: DynamicFeedLineState

    def __post_init__(self) -> None:
        """Require typed tank and line records before any equation is evaluated."""
        for name, kind in (("oxidizer_feed", RegulatedFeedState),
                           ("fuel_feed", RegulatedFeedState),
                           ("oxidizer_line", DynamicFeedLineState),
                           ("fuel_line", DynamicFeedLineState)):
            if not isinstance(getattr(self, name), kind):
                raise DomainError(f"Bipropellant {name} requires a {kind.__name__}.")


@dataclass(frozen=True, slots=True)
class ComponentDepletionPrediction:
    """Instantaneous reserve-crossing prediction for component-limited cutoff.

    Attributes
    ----------
    oxidizer_time_to_reserve_s, fuel_time_to_reserve_s : float
        Usable mass above the tank reserve divided by the current
        tank-to-line flow, s; ``inf`` for zero or reverse flow.
    limiting_component : str or None
        ``"oxidizer"`` or ``"fuel"``, whichever reaches its reserve first
        (oxidizer on a tie); ``None`` when neither time is finite.

    Notes
    -----
    The values are local linear predictions from the current flow.  The
    authoritative cutoff is the event root on each propagated tank state.
    Separately pressurized legs deplete at rates set by their own feed
    pressures [1]_.

    References
    ----------
    .. [1] NASA, "Simulation Analysis of Computer-Controlled Pressurization for
           Mixture Ratio Control". https://ntrs.nasa.gov/citations/20050206420
    """

    oxidizer_time_to_reserve_s: float
    fuel_time_to_reserve_s: float
    limiting_component: str | None


class BipropellantMixtureRatioEventKind(str, Enum):
    """Engine mixture-ratio validity boundary kinds.

    Members
    -------
    MINIMUM
        Lower engine mixture-ratio limit.
    MAXIMUM
        Upper engine mixture-ratio limit.

    Notes
    -----
    Measured short-pulse bipropellant flow-ratio transients [1]_ show why the
    mixture ratio is monitored as an event rather than assumed constant.

    References
    ----------
    .. [1] NASA, "Space Simulation Experiments on Reaction Control System Thruster
           Plumes". https://ntrs.nasa.gov/citations/19730002016
    """

    MINIMUM = "minimum_engine_mixture_ratio"
    MAXIMUM = "maximum_engine_mixture_ratio"


@dataclass(frozen=True, slots=True)
class BipropellantMixtureRatioEventSurface:
    """Active-engine mixture-ratio event function, positive inside its limit.

    Parameters
    ----------
    event_id : str
        Event identifier.
    kind : BipropellantMixtureRatioEventKind or str
        Boundary kind; a string is normalized to the enum.
    boundary_value : float
        Positive mixture-ratio boundary.
    stops_burn : bool, default True
        Whether a crossing terminates the burn.

    Raises
    ------
    DomainError
        If the kind is unknown, the boundary is not positive, or
        ``stops_burn`` is not a boolean.

    Notes
    -----
    The sign convention follows a propagator event detector [1]_: the
    function is positive inside the valid region and crosses zero at the
    boundary.  The monitored ratio responds to the separate feed pressures
    [2]_.

    References
    ----------
    .. [1] Orekit 13.1.5, "PropulsionModel / EventDetector API", API documentation of
           org.orekit.forces.maneuvers.propulsion.PropulsionModel.
    .. [2] NASA, "Simulation Analysis of Computer-Controlled Pressurization for
           Mixture Ratio Control". https://ntrs.nasa.gov/citations/20050206420
    """

    event_id: str
    kind: BipropellantMixtureRatioEventKind
    boundary_value: float
    stops_burn: bool = True

    def __post_init__(self) -> None:
        """Normalize the event kind and reject malformed cutoff metadata."""
        coerce_fields(self, "Bipropellant mixture-ratio event")
        try:
            kind = BipropellantMixtureRatioEventKind(self.kind)
        except (TypeError, ValueError) as error:
            raise DomainError("Unknown bipropellant mixture-ratio event kind.") from error
        object.__setattr__(self, "kind", kind)
        if self.boundary_value <= 0.0:
            raise DomainError("Mixture-ratio event boundary must be positive.")
        if not isinstance(self.stops_burn, bool):
            raise DomainError("Mixture-ratio event stops_burn must be a boolean.")

    def value(self, engine_mixture_ratio: float | None) -> float:
        """Evaluate the event function.

        Parameters
        ----------
        engine_mixture_ratio : float or None
            Active-engine mixture ratio, or ``None`` for an inactive engine.

        Returns
        -------
        float
            ``MR - boundary`` for a minimum, ``boundary - MR`` for a maximum,
            and ``1.0`` for an inactive engine, which therefore never
            triggers the event [1]_.

        Raises
        ------
        DomainError
            If an active mixture ratio is non-finite or not positive.

        References
        ----------
        .. [1] Orekit 13.1.5, "PropulsionModel / EventDetector API", API documentation of
               org.orekit.forces.maneuvers.propulsion.PropulsionModel.
        """

        if engine_mixture_ratio is None:
            return 1.0
        engine_mixture_ratio = real(engine_mixture_ratio, "Active engine mixture ratio")
        if engine_mixture_ratio <= 0.0:
            raise DomainError("Active engine mixture ratio must be finite and positive.")
        if self.kind is BipropellantMixtureRatioEventKind.MINIMUM:
            return engine_mixture_ratio - self.boundary_value
        return self.boundary_value - engine_mixture_ratio


@dataclass(frozen=True, slots=True)
class BipropellantCoupledPartials:
    """Four-state line Jacobian and feed and trajectory sensitivity blocks.

    Attributes
    ----------
    line_state_order : tuple of str
        Four distinct state names: oxidizer line flow, oxidizer manifold
        pressure, fuel line flow and fuel manifold pressure.
    line_state_jacobian : tuple of tuple of float
        4 x 4 Jacobian of the line-state derivatives.
    line_derivative_wrt_feed_state : tuple of (str, 4-tuple of float)
        Line-derivative sensitivities to named feed states.
    engine_component_flow_wrt_conditions : tuple of (str, tuple of (str, float))
        Component engine-flow sensitivities to each operating condition.
    engine_mixture_ratio_wrt_conditions : tuple of (str, float)
        Mixture-ratio sensitivities to each operating condition.
    acceleration_wrt_line_state : tuple of (str, 3-tuple of float)
        Acceleration sensitivities to each line state.

    Raises
    ------
    DomainError
        If a block has the wrong shape, a value is non-finite, or names are
        duplicated.

    Notes
    -----
    The line states are those of an inertance/compliance feedline model [1]_;
    the Jacobian supports implicit integration of the coupled transient [2]_.

    References
    ----------
    .. [1] NASA, "Liquid Rocket Propellant Feedline Dynamics".
           https://ntrs.nasa.gov/citations/19740028545
    .. [2] NASA, "Transient Mathematical Modeling for Liquid Rocket Engine
           Systems". https://ntrs.nasa.gov/citations/20040000363
    """

    line_state_order: tuple[str, str, str, str]
    line_state_jacobian: tuple[
        tuple[float, float, float, float],
        tuple[float, float, float, float],
        tuple[float, float, float, float],
        tuple[float, float, float, float],
    ]
    line_derivative_wrt_feed_state: tuple[tuple[str, tuple[float, float, float, float]], ...]
    engine_component_flow_wrt_conditions: tuple[tuple[str, tuple[tuple[str, float], ...]], ...]
    engine_mixture_ratio_wrt_conditions: tuple[tuple[str, float], ...]
    acceleration_wrt_line_state: tuple[tuple[str, tuple[float, float, float]], ...]

    def __post_init__(self) -> None:
        """Require finite, uniquely named and correctly shaped sensitivity blocks."""
        order = strings(self.line_state_order, "Coupled line state order")
        if len(order) != 4 or len(set(order)) != 4:
            raise DomainError("Coupled line state order requires four distinct names.")
        object.__setattr__(self, "line_state_order", order)
        matrix = real_matrix(self.line_state_jacobian, "Coupled line Jacobian")
        if len(matrix) != 4 or any(len(row) != 4 for row in matrix):
            raise DomainError("Coupled line Jacobian must have shape 4 by 4.")
        object.__setattr__(self, "line_state_jacobian", matrix)
        object.__setattr__(self, "line_derivative_wrt_feed_state", _finite_columns(
            self.line_derivative_wrt_feed_state, 4, "Coupled feed-state partials"))
        object.__setattr__(self, "acceleration_wrt_line_state", _finite_columns(
            self.acceleration_wrt_line_state, 3, "Coupled acceleration partials"))
        object.__setattr__(self, "engine_mixture_ratio_wrt_conditions", _finite_pairs(
            self.engine_mixture_ratio_wrt_conditions, "Mixture-ratio condition partials"))
        components = []
        for item in records(self.engine_component_flow_wrt_conditions, (list, tuple),
                            "Engine component partials"):
            if len(item) != 2:
                raise DomainError("Engine component partials require named condition tables.")
            components.append((text(item[0], "Engine partial component"),
                               _finite_pairs(item[1], "Engine component condition partials")))
        if len({name for name, _ in components}) != len(components):
            raise DomainError("Engine partial component names must be unique.")
        object.__setattr__(self, "engine_component_flow_wrt_conditions", tuple(components))


@dataclass(frozen=True, slots=True)
class BipropellantDynamicFeedEvaluation:
    """Simultaneous two-leg feed, engine, conservation and Jacobian result.

    Attributes
    ----------
    propulsion : PropulsionDynamicsEvaluation
        Bridge evaluation of force, mass rate and state drains.
    oxidizer_feed, fuel_feed : RegulatedFeedEvaluation
        Feed evaluations at the current line flows.
    oxidizer_line, fuel_line : DynamicFeedLineEvaluation
        Line evaluations.
    oxidizer_engine_flow_kg_s, fuel_engine_flow_kg_s : float
        Component engine flows, kg/s.
    engine_mixture_ratio : float or None
        Engine O/F; ``None`` when the engine is inactive.
    line_mixture_ratio : float or None
        Line-flow O/F; ``None`` unless the fuel line flow is positive.
    stored_liquid_mass_rate_kg_s : float
        Rate of change of liquid stored in both compliant lines, kg/s.
    component_to_total_flow_closure_error_kg_s : float
        ``|O + F - total engine flow|``, kg/s.
    vehicle_mass_closure_error_kg_s : float
        Residual of tank, line-storage and engine mass rates, kg/s.
    depletion : ComponentDepletionPrediction
        Local reserve-crossing prediction.
    partials : BipropellantCoupledPartials
        Composed Jacobian and sensitivity blocks.

    Raises
    ------
    DomainError
        If a scalar diagnostic is non-finite.

    Notes
    -----
    The mass balances follow a nodal feed-network formulation [2]_ with
    independently pressurized oxidizer and fuel legs [1]_.

    References
    ----------
    .. [1] NASA, "Simulation Analysis of Computer-Controlled Pressurization for
           Mixture Ratio Control". https://ntrs.nasa.gov/citations/20050206420
    .. [2] NASA, "Nodal Modeling of Liquid Propellant Feed and Pressurization
           System". https://ntrs.nasa.gov/citations/20240003493
    """

    propulsion: PropulsionDynamicsEvaluation
    oxidizer_feed: RegulatedFeedEvaluation
    oxidizer_line: DynamicFeedLineEvaluation
    fuel_feed: RegulatedFeedEvaluation
    fuel_line: DynamicFeedLineEvaluation
    oxidizer_engine_flow_kg_s: float
    fuel_engine_flow_kg_s: float
    engine_mixture_ratio: float | None
    line_mixture_ratio: float | None
    stored_liquid_mass_rate_kg_s: float
    component_to_total_flow_closure_error_kg_s: float
    vehicle_mass_closure_error_kg_s: float
    depletion: ComponentDepletionPrediction
    partials: BipropellantCoupledPartials

    def __post_init__(self) -> None:
        """Reject non-finite scalar flow, ratio and conservation diagnostics."""
        coerce_fields(self, "Bipropellant evaluation")

    @property
    def reverse_flow_event_extension(self) -> bool:
        """Whether either feed leg is evaluated beyond its zero-flow cutoff."""
        return (self.oxidizer_feed.reverse_flow_event_extension
                or self.fuel_feed.reverse_flow_event_extension)


@dataclass(frozen=True, slots=True)
class BipropellantDynamicFeedCoupling:
    """Two-tank, two-line, stream-resolved propulsion model for a propagator.

    Parameters
    ----------
    configuration : BipropellantDynamicFeedConfiguration
        State names, pressure-condition names and mixture-ratio limits.
    bridge : PropagatorPropulsionBridge
        Propagation bridge with a stream-resolved performance surface and
        complete stream binding.
    oxidizer_feed, fuel_feed : RegulatedFeedSystem
        Regulated feeds; their static feed-line resistance must be zero.
    oxidizer_line, fuel_line : DynamicFeedLine
        Dynamic inertance/compliance liquid lines.

    Raises
    ------
    DomainError
        If a field has the wrong type, the bridge lacks a stream-resolved
        surface or complete binding, a pressure condition is not a surface
        axis, the bound states differ from the configured oxidizer and fuel
        states, or a feed declares static line resistance.

    Notes
    -----
    The bridge uses a `RectilinearStreamPerformanceSurface` whose stream
    bindings collapse to exactly the configured oxidizer and fuel states.
    Several artifact streams may drain one leg, so coolant or turbine flows
    stay assigned to their physical source.  Oxidizer and fuel are pressurized
    independently [1]_, short-pulse flow-ratio transients are measured in
    [2]_, the network follows a nodal feed formulation [3]_, and the
    propagation boundary follows [4]_.

    References
    ----------
    .. [1] NASA, "Simulation Analysis of Computer-Controlled Pressurization for
           Mixture Ratio Control". https://ntrs.nasa.gov/citations/20050206420
    .. [2] NASA, "Space Simulation Experiments on Reaction Control System Thruster
           Plumes". https://ntrs.nasa.gov/citations/19730002016
    .. [3] NASA, "Nodal Modeling of Liquid Propellant Feed and Pressurization
           System". https://ntrs.nasa.gov/citations/20240003493
    .. [4] Orekit 13.1.5, "PropulsionModel / EventDetector API", API documentation of
           org.orekit.forces.maneuvers.propulsion.PropulsionModel.
    """

    configuration: BipropellantDynamicFeedConfiguration
    bridge: PropagatorPropulsionBridge
    oxidizer_feed: RegulatedFeedSystem
    oxidizer_line: DynamicFeedLine
    fuel_feed: RegulatedFeedSystem
    fuel_line: DynamicFeedLine

    def __post_init__(self) -> None:
        for name, kind in (("configuration", BipropellantDynamicFeedConfiguration),
                           ("bridge", PropagatorPropulsionBridge),
                           ("oxidizer_feed", RegulatedFeedSystem),
                           ("fuel_feed", RegulatedFeedSystem),
                           ("oxidizer_line", DynamicFeedLine),
                           ("fuel_line", DynamicFeedLine)):
            if not isinstance(getattr(self, name), kind):
                raise DomainError(f"Bipropellant {name} requires a {kind.__name__}.")
        config = self.configuration
        if self.bridge.stream_performance_surface is None:
            raise DomainError(
                "Bipropellant coupling requires a stream-resolved performance surface."
            )
        axis_names = {axis.name for axis in self.bridge.stream_performance_surface.axes}
        required_axes = {
            config.oxidizer_pressure_condition_name,
            config.fuel_pressure_condition_name,
        }
        if not required_axes.issubset(axis_names):
            raise DomainError("Bipropellant pressure conditions must be performance axes.")
        if not self.bridge.require_complete_stream_binding:
            raise DomainError("Bipropellant coupling requires complete stream binding.")
        bound_states = {binding.state_name for binding in self.bridge.stream_bindings}
        if bound_states != {config.oxidizer_state_name, config.fuel_state_name}:
            raise DomainError(
                "Bridge streams must bind exclusively to the configured oxidizer and fuel states."
            )
        for feed in (self.oxidizer_feed, self.fuel_feed):
            if feed.configuration.feed_line_resistance_pa_s2_kg2 != 0.0:
                raise DomainError(
                    "Static feed resistance must be zero when dynamic lines are used."
                )

    @property
    def model_hash(self) -> str:
        """SHA-256 identity of the complete dual-feed coupling manifest.

        Model identity supports the configuration-traceability requirements of
        [1]_.

        References
        ----------
        .. [1] NASA-STD-7009B. NASA, "Standard for Models and Simulations", 2024.
               https://standards.nasa.gov/standard/NASA/NASA-STD-7009
        """

        return canonical_sha256(self.manifest(include_hash=False))

    def manifest(self, *, include_hash: bool = True) -> dict[str, object]:
        """Return a deterministic, hash-bound integration manifest.

        Parameters
        ----------
        include_hash : bool, default True
            Whether to add ``model_hash``, the canonical SHA-256 of the
            manifest without it.

        Returns
        -------
        dict
            Schema, model identifier, component model hashes, configuration,
            line-state order, numerics policy and reference identifiers.

        Notes
        -----
        The manifest records model identity and configuration for model
        traceability [1]_; the referenced mixture-ratio study is [2]_.

        References
        ----------
        .. [1] NASA-STD-7009B. NASA, "Standard for Models and Simulations", 2024.
               https://standards.nasa.gov/standard/NASA/NASA-STD-7009
        .. [2] NASA, "Simulation Analysis of Computer-Controlled Pressurization for
               Mixture Ratio Control". https://ntrs.nasa.gov/citations/20050206420
        """

        config = self.configuration
        payload: dict[str, object] = {
            "schema": config.schema,
            "model_id": config.model_id,
            "bridge_model_hash": self.bridge.model_hash,
            "oxidizer_feed_model_hash": self.oxidizer_feed.model_hash,
            "oxidizer_line_model_hash": self.oxidizer_line.model_hash,
            "fuel_feed_model_hash": self.fuel_feed.model_hash,
            "fuel_line_model_hash": self.fuel_line.model_hash,
            "configuration": {
                "oxidizer_state_name": config.oxidizer_state_name,
                "fuel_state_name": config.fuel_state_name,
                "oxidizer_pressure_condition_name": (config.oxidizer_pressure_condition_name),
                "fuel_pressure_condition_name": config.fuel_pressure_condition_name,
                "minimum_engine_mixture_ratio": config.minimum_engine_mixture_ratio,
                "maximum_engine_mixture_ratio": config.maximum_engine_mixture_ratio,
                "mass_closure_absolute_tolerance_kg_s": (
                    config.mass_closure_absolute_tolerance_kg_s
                ),
            },
            "line_state_order": [
                "oxidizer_mass_flow_kg_s",
                "oxidizer_manifold_pressure_pa",
                "fuel_mass_flow_kg_s",
                "fuel_manifold_pressure_pa",
            ],
            "mixture_ratio_numerics": "exact_binary_rational_then_finite_output_v1",
            "reference_ids": list(BIPROPELLANT_DYNAMIC_FEED_REFERENCE_IDS),
        }
        if include_hash:
            payload["model_hash"] = canonical_sha256(payload)
        return payload

    def mixture_ratio_surfaces(
        self,
    ) -> tuple[BipropellantMixtureRatioEventSurface, ...]:
        """Return the minimum and maximum active-engine O/F event functions.

        Returns
        -------
        tuple of BipropellantMixtureRatioEventSurface
            Burn-stopping surfaces at the configured lower and upper limits,
            in the event-detector convention of [1]_.  Short-pulse flow-ratio
            transients are measured in [2]_.

        References
        ----------
        .. [1] Orekit 13.1.5, "PropulsionModel / EventDetector API", API documentation of
               org.orekit.forces.maneuvers.propulsion.PropulsionModel.
        .. [2] NASA, "Space Simulation Experiments on Reaction Control System Thruster
               Plumes". https://ntrs.nasa.gov/citations/19730002016
        """

        config = self.configuration
        return (
            BipropellantMixtureRatioEventSurface(
                "bipropellant:mixture-ratio:minimum",
                BipropellantMixtureRatioEventKind.MINIMUM,
                config.minimum_engine_mixture_ratio,
            ),
            BipropellantMixtureRatioEventSurface(
                "bipropellant:mixture-ratio:maximum",
                BipropellantMixtureRatioEventKind.MAXIMUM,
                config.maximum_engine_mixture_ratio,
            ),
        )

    def _depletion_prediction(
        self, state: BipropellantDynamicFeedState
    ) -> ComponentDepletionPrediction:
        """Predict each leg's time to its reserve from the current line flow.

        Each time is the usable mass above the tank reserve divided by the
        tank-to-line flow, a first-order mass balance of the nodal network
        [1]_; zero or reverse flow gives ``inf``.

        References
        ----------
        .. [1] NASA, "Nodal Modeling of Liquid Propellant Feed and Pressurization
               System". https://ntrs.nasa.gov/citations/20240003493
        """
        oxidizer_usable = max(
            0.0,
            state.oxidizer_feed.propellant_mass_kg
            - self.oxidizer_feed.configuration.minimum_propellant_mass_kg,
        )
        fuel_usable = max(
            0.0,
            state.fuel_feed.propellant_mass_kg
            - self.fuel_feed.configuration.minimum_propellant_mass_kg,
        )
        oxidizer_time = (
            oxidizer_usable / state.oxidizer_line.mass_flow_kg_s
            if state.oxidizer_line.mass_flow_kg_s > 0.0
            else inf
        )
        fuel_time = (
            fuel_usable / state.fuel_line.mass_flow_kg_s
            if state.fuel_line.mass_flow_kg_s > 0.0
            else inf
        )
        limiting = None
        if isfinite(oxidizer_time) or isfinite(fuel_time):
            limiting = "oxidizer" if oxidizer_time <= fuel_time else "fuel"
        return ComponentDepletionPrediction(oxidizer_time, fuel_time, limiting)

    def evaluate(
        self,
        state: BipropellantDynamicFeedState,
        *,
        relative_time_s: float,
        vehicle_mass_kg: float,
        direction: tuple[float, float, float],
        other_operating_conditions: Mapping[str, float] | None = None,
        other_additional_states: Mapping[str, float] | None = None,
        parameter_overrides: Mapping[str, float] | None = None,
        ambient_pressure_pa: float | None = None,
        side: str = "right",
        apply_inventory_limits: bool = True,
        boundary_tolerance_s: float | None = None,
    ) -> BipropellantDynamicFeedEvaluation:
        """Evaluate both feeds, the stream-resolved engine and the composed Jacobian.

        Parameters
        ----------
        state : BipropellantDynamicFeedState
            Feed and line states of both legs.
        relative_time_s : float
            Time relative to the burn command, s.
        vehicle_mass_kg : float
            Current vehicle mass, kg.
        direction : tuple of float
            Thrust direction in the propagation frame.
        other_operating_conditions : Mapping[str, float] or None
            Operating conditions other than the two manifold pressures, which
            the lines own.
        other_additional_states : Mapping[str, float] or None
            Additional states other than the two tank masses, which the feeds
            own.
        parameter_overrides : Mapping[str, float] or None
            Bridge parameter overrides.
        ambient_pressure_pa : float or None
            Ambient pressure, Pa.
        side : str, default "right"
            Side of a schedule discontinuity at which to evaluate.
        apply_inventory_limits : bool, default True
            Whether bridge inventory limits gate the engine.
        boundary_tolerance_s : float or None
            Schedule-boundary tolerance, s.

        Returns
        -------
        BipropellantDynamicFeedEvaluation
            Feed, line and engine evaluations, mixture ratios, closure errors,
            depletion prediction and composed partials.

        Raises
        ------
        DomainError
            If the state has the wrong type, a caller supplies an owned
            condition or state, component flows do not close to total engine
            flow, an active engine has non-positive fuel flow, or tank,
            line-storage and vehicle mass rates do not close within tolerance.

        Notes
        -----
        The line manifold pressures become the engine-map pressure conditions,
        the tank masses become bridge additional states, and the engine demand
        drains each line.  Each manifold obeys
        ``C * dp/dt = mdot_line - mdot_engine(p_ox, p_fuel)`` [3]_, so the
        pressure rows of the Jacobian hold ``-d(mdot_engine)/dp / C``, which
        couples the oxidizer and fuel legs through the engine map [1]_.  Mass
        conservation is checked over tanks, line storage and the vehicle
        [2]_.

        Negative line states are mathematical continuations used to locate the
        zero-flow cutoff and are marked ``reverse_flow_event_extension``; they
        are not predictions of sustained reverse flow in hardware.

        References
        ----------
        .. [1] NASA, "Simulation Analysis of Computer-Controlled Pressurization for
               Mixture Ratio Control". https://ntrs.nasa.gov/citations/20050206420
        .. [2] NASA, "Nodal Modeling of Liquid Propellant Feed and Pressurization
               System". https://ntrs.nasa.gov/citations/20240003493
        .. [3] NASA, "Liquid Rocket Propellant Feedline Dynamics".
               https://ntrs.nasa.gov/citations/19740028545
        """

        if not isinstance(state, BipropellantDynamicFeedState):
            raise DomainError("Bipropellant evaluation requires a dynamic-feed state.")
        config = self.configuration
        conditions = real_mapping(other_operating_conditions, "Other operating conditions")
        owned_conditions = {
            config.oxidizer_pressure_condition_name,
            config.fuel_pressure_condition_name,
        }
        conflict = owned_conditions & set(conditions)
        if conflict:
            raise DomainError(
                "Bipropellant lines own pressure conditions: " + ", ".join(sorted(conflict))
            )
        conditions[config.oxidizer_pressure_condition_name] = (
            state.oxidizer_line.manifold_pressure_pa
        )
        conditions[config.fuel_pressure_condition_name] = state.fuel_line.manifold_pressure_pa

        additional_states = real_mapping(other_additional_states, "Other additional states")
        owned_states = {config.oxidizer_state_name, config.fuel_state_name}
        state_conflict = owned_states & set(additional_states)
        if state_conflict:
            raise DomainError(
                "Bipropellant feeds own additional states: " + ", ".join(sorted(state_conflict))
            )
        additional_states[config.oxidizer_state_name] = state.oxidizer_feed.propellant_mass_kg
        additional_states[config.fuel_state_name] = state.fuel_feed.propellant_mass_kg

        propulsion = self.bridge.evaluate(
            relative_time_s=relative_time_s,
            vehicle_mass_kg=vehicle_mass_kg,
            direction=direction,
            additional_states=additional_states,
            operating_conditions=conditions,
            parameter_overrides=parameter_overrides,
            ambient_pressure_pa=ambient_pressure_pa,
            side=side,
            apply_inventory_limits=apply_inventory_limits,
            boundary_tolerance_s=boundary_tolerance_s,
        )
        state_drains = dict(propulsion.additional_state_derivatives)
        oxidizer_engine_flow = -state_drains[config.oxidizer_state_name]
        fuel_engine_flow = -state_drains[config.fuel_state_name]
        total_engine_flow = -propulsion.mass_derivative_kg_s
        component_closure = abs(oxidizer_engine_flow + fuel_engine_flow - total_engine_flow)
        if component_closure > config.mass_closure_absolute_tolerance_kg_s:
            raise DomainError("Bipropellant component flows do not close to total engine flow.")
        if propulsion.active and fuel_engine_flow <= 0.0:
            raise DomainError("An active bipropellant engine requires positive fuel flow.")
        engine_mixture_ratio = (
            _finite_quotient(oxidizer_engine_flow, fuel_engine_flow, "Engine mixture ratio")
            if propulsion.active else None
        )

        oxidizer_feed = self.oxidizer_feed.evaluate(
            state.oxidizer_feed, state.oxidizer_line.mass_flow_kg_s,
            for_event_location=state.oxidizer_line.mass_flow_kg_s < 0.0,
        )
        fuel_feed = self.fuel_feed.evaluate(
            state.fuel_feed, state.fuel_line.mass_flow_kg_s,
            for_event_location=state.fuel_line.mass_flow_kg_s < 0.0,
        )
        oxidizer_line = self.oxidizer_line.evaluate(
            state.oxidizer_line,
            upstream_pressure_pa=oxidizer_feed.tank_pressure_pa,
            engine_mass_flow_kg_s=oxidizer_engine_flow,
        )
        fuel_line = self.fuel_line.evaluate(
            state.fuel_line,
            upstream_pressure_pa=fuel_feed.tank_pressure_pa,
            engine_mass_flow_kg_s=fuel_engine_flow,
        )
        line_mixture_ratio = (
            _finite_quotient(state.oxidizer_line.mass_flow_kg_s,
                             state.fuel_line.mass_flow_kg_s, "Line mixture ratio")
            if state.fuel_line.mass_flow_kg_s > 0.0
            else None
        )
        stored_rate = (
            oxidizer_line.stored_liquid_mass_rate_kg_s + fuel_line.stored_liquid_mass_rate_kg_s
        )
        vehicle_closure = abs(fsum((
            oxidizer_feed.derivative.propellant_mass_kg_s,
            fuel_feed.derivative.propellant_mass_kg_s,
            stored_rate,
            total_engine_flow,
        )))
        if vehicle_closure > config.mass_closure_absolute_tolerance_kg_s:
            raise DomainError("Bipropellant tank, line-storage, and vehicle mass do not close.")

        condition_names = tuple(conditions)
        zero_by_condition = {name: 0.0 for name in condition_names}
        component_flow_partials = {
            config.oxidizer_state_name: dict(zero_by_condition),
            config.fuel_state_name: dict(zero_by_condition),
        }
        acceleration_by_condition = {name: (0.0, 0.0, 0.0) for name in condition_names}
        if propulsion.active:
            acceleration_by_condition.update(dict(propulsion.partials.acceleration_wrt_conditions))
            for (
                condition_name,
                state_values,
            ) in propulsion.partials.additional_state_rates_wrt_conditions:
                for state_name, derivative in state_values:
                    if state_name in component_flow_partials:
                        component_flow_partials[state_name][condition_name] = -derivative

        oxidizer_pressure = config.oxidizer_pressure_condition_name
        fuel_pressure = config.fuel_pressure_condition_name
        oxidizer_demand = component_flow_partials[config.oxidizer_state_name]
        fuel_demand = component_flow_partials[config.fuel_state_name]
        oxidizer_local = dict(oxidizer_line.partials.derivative_wrt_state)
        fuel_local = dict(fuel_line.partials.derivative_wrt_state)
        oxidizer_compliance = self.oxidizer_line.configuration.compliance_kg_pa
        fuel_compliance = self.fuel_line.configuration.compliance_kg_pa
        jacobian = (
            (
                oxidizer_local["mass_flow_kg_s"][0],
                oxidizer_local["manifold_pressure_pa"][0],
                0.0,
                0.0,
            ),
            (
                oxidizer_local["mass_flow_kg_s"][1],
                -oxidizer_demand[oxidizer_pressure] / oxidizer_compliance,
                0.0,
                -oxidizer_demand[fuel_pressure] / oxidizer_compliance,
            ),
            (
                0.0,
                0.0,
                fuel_local["mass_flow_kg_s"][0],
                fuel_local["manifold_pressure_pa"][0],
            ),
            (
                0.0,
                -fuel_demand[oxidizer_pressure] / fuel_compliance,
                fuel_local["mass_flow_kg_s"][1],
                -fuel_demand[fuel_pressure] / fuel_compliance,
            ),
        )

        feed_partials: list[tuple[str, tuple[float, float, float, float]]] = []
        oxidizer_upstream_gain = oxidizer_line.partials.derivative_wrt_upstream_pressure[0]
        for name, derivative in oxidizer_feed.pressure_partials.injector_pressure_wrt_state:
            feed_partials.append(
                (f"oxidizer.{name}", (oxidizer_upstream_gain * derivative, 0.0, 0.0, 0.0))
            )
        fuel_upstream_gain = fuel_line.partials.derivative_wrt_upstream_pressure[0]
        for name, derivative in fuel_feed.pressure_partials.injector_pressure_wrt_state:
            feed_partials.append((f"fuel.{name}", (0.0, 0.0, fuel_upstream_gain * derivative, 0.0)))

        zero_acceleration = (0.0, 0.0, 0.0)
        state_order = (
            "oxidizer_mass_flow_kg_s",
            "oxidizer_manifold_pressure_pa",
            "fuel_mass_flow_kg_s",
            "fuel_manifold_pressure_pa",
        )
        acceleration_wrt_line_state = (
            (state_order[0], zero_acceleration),
            (state_order[1], acceleration_by_condition[oxidizer_pressure]),
            (state_order[2], zero_acceleration),
            (state_order[3], acceleration_by_condition[fuel_pressure]),
        )
        mixture_ratio_partials: tuple[tuple[str, float], ...] = ()
        if engine_mixture_ratio is not None:
            mixture_ratio_partials = tuple(
                (
                    name,
                    _mixture_ratio_partial(oxidizer_engine_flow, fuel_engine_flow,
                                           oxidizer_demand[name], fuel_demand[name]),
                )
                for name in condition_names
            )

        return BipropellantDynamicFeedEvaluation(
            propulsion=propulsion,
            oxidizer_feed=oxidizer_feed,
            oxidizer_line=oxidizer_line,
            fuel_feed=fuel_feed,
            fuel_line=fuel_line,
            oxidizer_engine_flow_kg_s=oxidizer_engine_flow,
            fuel_engine_flow_kg_s=fuel_engine_flow,
            engine_mixture_ratio=engine_mixture_ratio,
            line_mixture_ratio=line_mixture_ratio,
            stored_liquid_mass_rate_kg_s=stored_rate,
            component_to_total_flow_closure_error_kg_s=component_closure,
            vehicle_mass_closure_error_kg_s=vehicle_closure,
            depletion=self._depletion_prediction(state),
            partials=BipropellantCoupledPartials(
                line_state_order=state_order,
                line_state_jacobian=jacobian,
                line_derivative_wrt_feed_state=tuple(feed_partials),
                engine_component_flow_wrt_conditions=(
                    (
                        "oxidizer",
                        tuple((name, oxidizer_demand[name]) for name in condition_names),
                    ),
                    (
                        "fuel",
                        tuple((name, fuel_demand[name]) for name in condition_names),
                    ),
                ),
                engine_mixture_ratio_wrt_conditions=mixture_ratio_partials,
                acceleration_wrt_line_state=acceleration_wrt_line_state,
            ),
        )
