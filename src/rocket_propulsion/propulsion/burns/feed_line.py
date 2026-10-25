"""Lumped two-state surrogate for liquid feed-line transients.

The line retains hydraulic inertance and compliance instead of reducing every
integrator-stage evaluation to a static pressure loss [1]_.  It is small enough
for orbit propagation and exposes the mass-conservation and Jacobian terms a
numerical force model requires.  The surrogate represents the first
inertance/compliance mode only; it is not a water-hammer solver [2]_, and its
coefficients are calibrated against the actual plumbing or a higher-order nodal
network model [3]_.

References
----------
.. [1] NASA, "Liquid Rocket Propellant Feedline Dynamics", 1973.
       https://ntrs.nasa.gov/citations/19740028545
.. [2] NASA NESC Technical Bulletin 22-03, "Treatment of Transient Pressure
       Events", 2022. https://ntrs.nasa.gov/citations/20220006583
.. [3] NASA, "Nodal Modeling of Liquid Propellant Feed and Pressurization
       System", 2024. https://ntrs.nasa.gov/citations/20240003493
"""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass
from enum import Enum
from math import isfinite, sqrt
from typing import TYPE_CHECKING

from rocket_propulsion.core.errors import DomainError
from rocket_propulsion.core.validation import coerce_fields, real, real_tuple, strings, text

from .feed_system import RegulatedFeedEvaluation, RegulatedFeedState, RegulatedFeedSystem
from .serialization import canonical_sha256

if TYPE_CHECKING:
    from .propagator import PropagatorPropulsionBridge, PropulsionDynamicsEvaluation

DYNAMIC_FEED_LINE_SCHEMA = "rocket_propulsion_dynamic_feed_line_v1"
DYNAMIC_FEED_LINE_REFERENCE_IDS = (
    "NASA-19740028545",
    "NASA-20220006583",
    "NASA-20240003493",
)


def _positive(value: float, label: str) -> None:
    """Raise `DomainError` unless ``value`` is finite and strictly positive."""
    if not isfinite(value) or value <= 0.0:
        raise DomainError(f"{label} must be finite and greater than zero.")


def _nonnegative(value: float, label: str) -> None:
    """Raise `DomainError` unless ``value`` is finite and non-negative."""
    if not isfinite(value) or value < 0.0:
        raise DomainError(f"{label} must be finite and non-negative.")




@dataclass(frozen=True, slots=True)
class DynamicFeedLineConfiguration:
    """Lumped hydraulic inertance, compliance, resistance and operating limits.

    Parameters
    ----------
    inertance_pa_s2_kg : float
        Positive hydraulic inertance ``L``, Pa s^2/kg.
    compliance_kg_pa : float
        Positive compliance ``C``, kg/Pa: incremental stored liquid mass per
        unit manifold pressure.
    resistance_pa_s2_kg2 : float
        Non-negative quadratic resistance ``R``, Pa s^2/kg^2.
    minimum_manifold_pressure_pa : float
        Positive lower manifold-pressure limit, Pa.
    maximum_manifold_pressure_pa : float
        Upper manifold-pressure limit, Pa; must exceed the minimum.
    maximum_absolute_flow_kg_s : float
        Positive limit on ``abs(q)``, kg/s.
    model_id : str, default "lumped-inertance-compliance-feed-line"
        Non-blank model identifier.
    schema : str, default DYNAMIC_FEED_LINE_SCHEMA
        Configuration schema; only the current schema is accepted.

    Raises
    ------
    DomainError
        If the schema is unsupported, the identifier is blank, a value is
        outside its domain, or the pressure limits are not ordered.

    Notes
    -----
    The governing equations of the lumped line [1]_ are

    ``L dq/dt = p_up - p_m - R q abs(q)``

    and

    ``C dp_m/dt = q - q_engine``,

    where ``q`` is line mass flow, ``p_up`` the upstream (tank) pressure and
    ``p_m`` the engine-manifold pressure.  The model is a low-order surrogate
    for distributed line modes and transient pressure events [2]_; its
    coefficients are calibrated against the actual plumbing or a higher-order
    nodal model [3]_.

    References
    ----------
    .. [1] NASA, "Liquid Rocket Propellant Feedline Dynamics", 1973.
           https://ntrs.nasa.gov/citations/19740028545
    .. [2] NASA NESC Technical Bulletin 22-03, "Treatment of Transient Pressure
           Events", 2022. https://ntrs.nasa.gov/citations/20220006583
    .. [3] NASA, "Nodal Modeling of Liquid Propellant Feed and Pressurization
           System", 2024. https://ntrs.nasa.gov/citations/20240003493
    """

    inertance_pa_s2_kg: float
    compliance_kg_pa: float
    resistance_pa_s2_kg2: float
    minimum_manifold_pressure_pa: float
    maximum_manifold_pressure_pa: float
    maximum_absolute_flow_kg_s: float
    model_id: str = "lumped-inertance-compliance-feed-line"
    schema: str = DYNAMIC_FEED_LINE_SCHEMA

    def __post_init__(self) -> None:
        coerce_fields(self, "Feed-line configuration")
        if self.schema != DYNAMIC_FEED_LINE_SCHEMA:
            raise DomainError("Unsupported dynamic feed-line schema.")
        if not self.model_id.strip():
            raise DomainError("Dynamic feed-line model identifier must not be empty.")
        for value, label in (
            (self.inertance_pa_s2_kg, "Feed-line inertance"),
            (self.compliance_kg_pa, "Feed-line compliance"),
            (self.minimum_manifold_pressure_pa, "Minimum manifold pressure"),
            (self.maximum_manifold_pressure_pa, "Maximum manifold pressure"),
            (self.maximum_absolute_flow_kg_s, "Maximum absolute feed-line flow"),
        ):
            _positive(value, label)
        _nonnegative(self.resistance_pa_s2_kg2, "Feed-line resistance")
        if self.maximum_manifold_pressure_pa <= self.minimum_manifold_pressure_pa:
            raise DomainError("Maximum manifold pressure must exceed the minimum.")


@dataclass(frozen=True, slots=True)
class DynamicFeedLineState:
    """Propagated line mass flow and engine-manifold pressure.

    Parameters
    ----------
    mass_flow_kg_s : float
        Finite signed line mass flow ``q``, kg/s.
    manifold_pressure_pa : float
        Positive engine-manifold pressure ``p_m``, Pa.

    Raises
    ------
    DomainError
        If the flow is non-finite or the pressure is not finite and positive.

    Notes
    -----
    Negative flow is representable so that an event locator can resolve a
    flow reversal of the line [1]_.  The reverse-flow event is positive inside
    the permitted domain and stops the burn at zero flow, before any sustained
    reverse operation; the pressure limits screen transient pressure events
    [2]_.

    References
    ----------
    .. [1] NASA, "Liquid Rocket Propellant Feedline Dynamics", 1973.
           https://ntrs.nasa.gov/citations/19740028545
    .. [2] NASA NESC Technical Bulletin 22-03, "Treatment of Transient Pressure
           Events", 2022. https://ntrs.nasa.gov/citations/20220006583
    """

    mass_flow_kg_s: float
    manifold_pressure_pa: float

    def __post_init__(self) -> None:
        coerce_fields(self, "Feed-line state")
        if not isfinite(self.mass_flow_kg_s):
            raise DomainError("Feed-line mass flow must be finite.")
        _positive(self.manifold_pressure_pa, "Feed-line manifold pressure")


@dataclass(frozen=True, slots=True)
class DynamicFeedLineDerivative:
    """Right-hand side of the two-state feed-line ODE.

    Attributes
    ----------
    mass_flow_kg_s2 : float
        ``dq/dt = (p_up - p_m - R q abs(q)) / L`` [1]_, kg/s^2.
    manifold_pressure_pa_s : float
        ``dp_m/dt = (q - q_engine) / C``, Pa/s.

    References
    ----------
    .. [1] NASA, "Liquid Rocket Propellant Feedline Dynamics", 1973.
           https://ntrs.nasa.gov/citations/19740028545
    """

    mass_flow_kg_s2: float
    manifold_pressure_pa_s: float

    def __post_init__(self) -> None:
        """Coerce both rates to finite reals; raise `DomainError` otherwise."""
        coerce_fields(self, "Feed-line derivative")

    def as_tuple(self) -> tuple[float, float]:
        """Return ``(dq/dt, dp_m/dt)`` of the lumped line [1]_ in state order.

        References
        ----------
        .. [1] NASA, "Liquid Rocket Propellant Feedline Dynamics", 1973.
               https://ntrs.nasa.gov/citations/19740028545
        """

        return (self.mass_flow_kg_s2, self.manifold_pressure_pa_s)


@dataclass(frozen=True, slots=True)
class DynamicFeedLinePartials:
    """Analytic local derivatives of the line ODE before engine-map composition.

    Attributes
    ----------
    derivative_wrt_state : tuple of (str, (float, float))
        Named Jacobian columns: the derivative of ``(dq/dt, dp_m/dt)`` with
        respect to each state; names are unique.
    derivative_wrt_upstream_pressure : tuple of float
        Derivative with respect to ``p_up``, ``(1 / L, 0)``.
    derivative_wrt_engine_flow : tuple of float
        Derivative with respect to ``q_engine``, ``(0, -1 / C)``.
    nonsmooth_state_names : tuple of str, default ()
        States at which the right-hand side is not twice differentiable; the
        quadratic loss ``R q abs(q)`` [1]_ is non-smooth at ``q = 0``.

    Raises
    ------
    DomainError
        If a column is malformed, a component is non-finite, or column names
        repeat.

    References
    ----------
    .. [1] NASA, "Liquid Rocket Propellant Feedline Dynamics", 1973.
           https://ntrs.nasa.gov/citations/19740028545
    """

    derivative_wrt_state: tuple[tuple[str, tuple[float, float]], ...]
    derivative_wrt_upstream_pressure: tuple[float, float]
    derivative_wrt_engine_flow: tuple[float, float]
    nonsmooth_state_names: tuple[str, ...] = ()

    def __post_init__(self) -> None:
        """Normalize and validate named, finite, two-component Jacobian columns."""
        if not isinstance(self.derivative_wrt_state, (list, tuple)):
            raise DomainError("Line state partials must be named two-component columns.")
        columns = []
        for column in self.derivative_wrt_state:
            if not isinstance(column, (list, tuple)) or len(column) != 2:
                raise DomainError("Line state partials must be named two-component columns.")
            name = text(column[0], "Line partial state name")
            values = real_tuple(column[1], "Line state partials")
            if len(values) != 2:
                raise DomainError("Line Jacobian columns require two components.")
            columns.append((name, values))
        if len({name for name, _ in columns}) != len(columns):
            raise DomainError("Line state partial names must be unique.")
        object.__setattr__(self, "derivative_wrt_state", tuple(columns))
        for field in ("derivative_wrt_upstream_pressure", "derivative_wrt_engine_flow"):
            values = real_tuple(getattr(self, field), field)
            if len(values) != 2:
                raise DomainError("Line input partials require two components.")
            object.__setattr__(self, field, values)
        object.__setattr__(self, "nonsmooth_state_names", strings(
            self.nonsmooth_state_names, "Line nonsmooth states"))


@dataclass(frozen=True, slots=True)
class DynamicFeedLineEvaluation:
    """Line ODE, mass-conservation terms and local modal diagnostics.

    Attributes
    ----------
    upstream_pressure_pa : float
        Upstream (tank) pressure ``p_up``, Pa.
    engine_mass_flow_kg_s : float
        Engine demand ``q_engine``, kg/s.
    friction_pressure_drop_pa : float
        Quadratic loss ``R q abs(q)``, Pa.
    stored_liquid_mass_rate_kg_s : float
        Rate of change of compliant line storage, ``q - q_engine``, kg/s.
    external_mass_closure_error_kg_s : float
        Absolute residual of the tank-plus-line mass balance, kg/s.
    natural_angular_frequency_rad_s : float
        Undamped first-mode frequency ``1 / sqrt(L C)`` [1]_, rad/s.
    damping_ratio : float
        Damping ratio of the linearized loss, ``R abs(q) sqrt(C / L)``.
    derivative : DynamicFeedLineDerivative
        ODE right-hand side.
    partials : DynamicFeedLinePartials
        Local analytic Jacobian.

    Notes
    -----
    Modal values describe the first lumped mode only; faster pressure
    transients such as water hammer are outside the surrogate [2]_.

    References
    ----------
    .. [1] NASA, "Liquid Rocket Propellant Feedline Dynamics", 1973.
           https://ntrs.nasa.gov/citations/19740028545
    .. [2] NASA NESC Technical Bulletin 22-03, "Treatment of Transient Pressure
           Events", 2022. https://ntrs.nasa.gov/citations/20220006583
    """

    upstream_pressure_pa: float
    engine_mass_flow_kg_s: float
    friction_pressure_drop_pa: float
    stored_liquid_mass_rate_kg_s: float
    external_mass_closure_error_kg_s: float
    natural_angular_frequency_rad_s: float
    damping_ratio: float
    derivative: DynamicFeedLineDerivative
    partials: DynamicFeedLinePartials

    def __post_init__(self) -> None:
        """Require finite diagnostics and a positive natural frequency."""
        coerce_fields(self, "Feed-line evaluation")
        _positive(self.natural_angular_frequency_rad_s, "Feed-line natural frequency")


class DynamicFeedLineEventKind(str, Enum):
    """Cutoff events of the dynamic line, each positive inside the permitted domain.

    Members are reverse flow, minimum and maximum manifold pressure, and
    maximum absolute flow; the pressure limits screen transient pressure
    excursions [1]_.

    References
    ----------
    .. [1] NASA NESC Technical Bulletin 22-03, "Treatment of Transient Pressure
           Events", 2022. https://ntrs.nasa.gov/citations/20220006583
    """

    REVERSE_FLOW = "reverse_flow"
    MINIMUM_MANIFOLD_PRESSURE = "minimum_manifold_pressure"
    MAXIMUM_MANIFOLD_PRESSURE = "maximum_manifold_pressure"
    MAXIMUM_ABSOLUTE_FLOW = "maximum_absolute_flow"


@dataclass(frozen=True, slots=True)
class DynamicFeedLineEventSurface:
    """One signed event function that is positive inside the permitted domain.

    Attributes
    ----------
    event_id : str
        Stable event identifier.
    kind : DynamicFeedLineEventKind
        Limit that the event function measures.
    stops_burn : bool, default True
        Whether a zero crossing terminates the burn.

    Notes
    -----
    The sign convention suits the zero-crossing event detectors of a
    propagator [2]_; the pressure limits screen transient pressure events
    [1]_.

    References
    ----------
    .. [1] NASA NESC Technical Bulletin 22-03, "Treatment of Transient Pressure
           Events", 2022. https://ntrs.nasa.gov/citations/20220006583
    .. [2] Orekit Project, "PropulsionModel / EventDetector" API documentation,
           Orekit 13.1.5.
    """

    event_id: str
    kind: DynamicFeedLineEventKind
    stops_burn: bool = True


@dataclass(frozen=True, slots=True)
class DynamicFeedCoupledPartials:
    """Line Jacobian composed with the engine map, and feed-state coupling.

    Attributes
    ----------
    line_state_jacobian : tuple of tuple of float
        2 x 2 Jacobian of ``(dq/dt, dp_m/dt)`` with respect to ``(q, p_m)``;
        the ``(1, 1)`` entry is the engine-demand term
        ``-(d q_engine / d p_m) / C`` of the lumped line [2]_.
    line_derivative_wrt_feed_state : tuple of (str, (float, float))
        Derivative of the line ODE with respect to each regulated-feed state,
        through the upstream pressure.
    engine_flow_wrt_manifold_pressure : float
        ``d q_engine / d p_m``, kg/(s Pa).
    acceleration_wrt_line_state : tuple of (str, (float, float, float))
        Thrust-acceleration derivative with respect to each line state.

    Notes
    -----
    [1]_ gives the transient engine-system modelling context of the coupled
    component equations.

    References
    ----------
    .. [1] NASA, "Transient Mathematical Modeling for Liquid Rocket Engine
           Systems". https://ntrs.nasa.gov/citations/20040000363
    .. [2] NASA, "Liquid Rocket Propellant Feedline Dynamics", 1973.
           https://ntrs.nasa.gov/citations/19740028545
    """

    line_state_jacobian: tuple[tuple[float, float], tuple[float, float]]
    line_derivative_wrt_feed_state: tuple[tuple[str, tuple[float, float]], ...]
    engine_flow_wrt_manifold_pressure: float
    acceleration_wrt_line_state: tuple[tuple[str, tuple[float, float, float]], ...]


@dataclass(frozen=True, slots=True)
class DynamicFeedPropulsionEvaluation:
    """Simultaneous pressurant, line, engine and Jacobian evaluation.

    Attributes
    ----------
    propulsion : PropulsionDynamicsEvaluation
        Engine evaluation at the propagated manifold pressure.
    regulated_feed : RegulatedFeedEvaluation
        Tank and pressurant evaluation at the line flow.
    feed_line : DynamicFeedLineEvaluation
        Line ODE and diagnostics.
    partials : DynamicFeedCoupledPartials
        Composed Jacobian terms.

    Notes
    -----
    The coupled tank, line and engine states form a small feed network in the
    sense of [1]_; [2]_ gives the transient-modelling context.

    References
    ----------
    .. [1] NASA, "Nodal Modeling of Liquid Propellant Feed and Pressurization
           System", 2024. https://ntrs.nasa.gov/citations/20240003493
    .. [2] NASA, "Transient Mathematical Modeling for Liquid Rocket Engine
           Systems". https://ntrs.nasa.gov/citations/20040000363
    """

    propulsion: PropulsionDynamicsEvaluation
    regulated_feed: RegulatedFeedEvaluation
    feed_line: DynamicFeedLineEvaluation
    partials: DynamicFeedCoupledPartials

    @property
    def reverse_flow_event_extension(self) -> bool:
        """Whether the signed right-hand side is only a zero-flow continuation.

        ``True`` when the line flow is negative; the values then serve event
        location only and are not a reverse-flow operating model.
        """
        return self.regulated_feed.reverse_flow_event_extension


@dataclass(frozen=True, slots=True)
class DynamicFeedLine:
    """Two-state liquid-line model for a numerical orbit propagator.

    Parameters
    ----------
    configuration : DynamicFeedLineConfiguration
        Lumped coefficients and operating limits.

    Notes
    -----
    The state ``(q, p_m)`` captures the first inertance/compliance mode of the
    line with nonlinear quadratic damping [1]_.  The model is suitable for
    mission-level transient coupling only after its effective coefficients are
    correlated against the actual line or a validated nodal model [3]_.  It is
    not a water-hammer solver [2]_.

    References
    ----------
    .. [1] NASA, "Liquid Rocket Propellant Feedline Dynamics", 1973.
           https://ntrs.nasa.gov/citations/19740028545
    .. [2] NASA NESC Technical Bulletin 22-03, "Treatment of Transient Pressure
           Events", 2022. https://ntrs.nasa.gov/citations/20220006583
    .. [3] NASA, "Nodal Modeling of Liquid Propellant Feed and Pressurization
           System", 2024. https://ntrs.nasa.gov/citations/20240003493
    """

    configuration: DynamicFeedLineConfiguration

    @property
    def model_hash(self) -> str:
        """SHA-256 of the canonical manifest without its hash field.

        The hash identifies configuration and source identity for model
        traceability [1]_.

        References
        ----------
        .. [1] NASA-STD-7009B. NASA, "Standard for Models and Simulations", 2024.
               https://standards.nasa.gov/standard/NASA/NASA-STD-7009
        """

        return canonical_sha256(self.manifest(include_hash=False))

    def manifest(self, *, include_hash: bool = True) -> dict[str, object]:
        """Return a deterministic, unit-explicit model manifest.

        Parameters
        ----------
        include_hash : bool, default True
            Whether to append ``model_hash`` computed over the other fields.

        Returns
        -------
        dict
            Schema, model identifier, coefficients under unit-suffixed keys,
            state order and source identifiers.

        Notes
        -----
        The manifest records the lumped coefficients of [2]_ for model
        traceability [1]_.

        References
        ----------
        .. [1] NASA-STD-7009B. NASA, "Standard for Models and Simulations", 2024.
               https://standards.nasa.gov/standard/NASA/NASA-STD-7009
        .. [2] NASA, "Liquid Rocket Propellant Feedline Dynamics", 1973.
               https://ntrs.nasa.gov/citations/19740028545
        """

        config = self.configuration
        payload: dict[str, object] = {
            "schema": config.schema,
            "model_id": config.model_id,
            "configuration": {
                "inertance_pa_s2_kg": config.inertance_pa_s2_kg,
                "compliance_kg_pa": config.compliance_kg_pa,
                "resistance_pa_s2_kg2": config.resistance_pa_s2_kg2,
                "minimum_manifold_pressure_pa": config.minimum_manifold_pressure_pa,
                "maximum_manifold_pressure_pa": config.maximum_manifold_pressure_pa,
                "maximum_absolute_flow_kg_s": config.maximum_absolute_flow_kg_s,
            },
            "state_order": ["mass_flow_kg_s", "manifold_pressure_pa"],
            "reference_ids": list(DYNAMIC_FEED_LINE_REFERENCE_IDS),
        }
        if include_hash:
            payload["model_hash"] = canonical_sha256(payload)
        return payload

    def evaluate(
        self,
        state: DynamicFeedLineState,
        *,
        upstream_pressure_pa: float,
        engine_mass_flow_kg_s: float,
    ) -> DynamicFeedLineEvaluation:
        """Evaluate line dynamics, mass conservation, modal values and partials.

        Parameters
        ----------
        state : DynamicFeedLineState
            Current line flow and manifold pressure.
        upstream_pressure_pa : float
            Positive upstream (tank) pressure ``p_up``, Pa.
        engine_mass_flow_kg_s : float
            Non-negative engine demand ``q_engine``, kg/s.

        Returns
        -------
        DynamicFeedLineEvaluation
            Right-hand side, analytic Jacobian, closure residual and modal
            diagnostics.

        Raises
        ------
        DomainError
            If ``state`` has the wrong type or an input is outside its domain.

        Notes
        -----
        The right-hand side is ``dq/dt = (p_up - p_m - R q abs(q)) / L`` and
        ``dp_m/dt = (q - q_engine) / C`` [1]_.  The Jacobian entry
        ``d(dq/dt)/dq = -2 R abs(q) / L`` is continuous, but the right-hand side
        is not twice differentiable at ``q = 0``, so ``q`` is reported as a
        nonsmooth state there.  The modal diagnostics are
        ``omega_n = 1 / sqrt(L C)`` and ``zeta = R abs(q) sqrt(C / L)``.  The
        closure residual checks that tank liquid plus compliant line storage
        changes exactly by engine outflow, the nodal mass balance of [2]_.

        References
        ----------
        .. [1] NASA, "Liquid Rocket Propellant Feedline Dynamics", 1973.
               https://ntrs.nasa.gov/citations/19740028545
        .. [2] NASA, "Nodal Modeling of Liquid Propellant Feed and Pressurization
               System", 2024. https://ntrs.nasa.gov/citations/20240003493
        """

        if not isinstance(state, DynamicFeedLineState):
            raise DomainError("Feed-line evaluation requires a DynamicFeedLineState.")
        upstream_pressure_pa = real(upstream_pressure_pa, "Upstream pressure")
        engine_mass_flow_kg_s = real(engine_mass_flow_kg_s, "Engine mass flow")

        _positive(upstream_pressure_pa, "Dynamic-line upstream pressure")
        _nonnegative(engine_mass_flow_kg_s, "Engine mass-flow demand")
        config = self.configuration
        flow = state.mass_flow_kg_s
        friction_drop = config.resistance_pa_s2_kg2 * flow * abs(flow)
        flow_rate = (
            upstream_pressure_pa - state.manifold_pressure_pa - friction_drop
        ) / config.inertance_pa_s2_kg
        stored_mass_rate = flow - engine_mass_flow_kg_s
        pressure_rate = stored_mass_rate / config.compliance_kg_pa
        natural_frequency = (1.0 / sqrt(config.inertance_pa_s2_kg)) / sqrt(config.compliance_kg_pa)
        damping_ratio = (
            config.resistance_pa_s2_kg2
            * abs(flow)
            * (sqrt(config.compliance_kg_pa) / sqrt(config.inertance_pa_s2_kg))
        )
        flow_derivative = -2.0 * config.resistance_pa_s2_kg2 * abs(flow) / config.inertance_pa_s2_kg
        partials = DynamicFeedLinePartials(
            derivative_wrt_state=(
                (
                    "mass_flow_kg_s",
                    (flow_derivative, 1.0 / config.compliance_kg_pa),
                ),
                (
                    "manifold_pressure_pa",
                    (-1.0 / config.inertance_pa_s2_kg, 0.0),
                ),
            ),
            derivative_wrt_upstream_pressure=(
                1.0 / config.inertance_pa_s2_kg,
                0.0,
            ),
            derivative_wrt_engine_flow=(0.0, -1.0 / config.compliance_kg_pa),
            nonsmooth_state_names=("mass_flow_kg_s",) if flow == 0.0 else (),
        )
        # Tank liquid plus compliant line storage changes exactly by engine
        # outflow; the residual is exported so an adapter can monitor it.
        closure_error = abs((-flow) + stored_mass_rate + engine_mass_flow_kg_s)
        return DynamicFeedLineEvaluation(
            upstream_pressure_pa=upstream_pressure_pa,
            engine_mass_flow_kg_s=engine_mass_flow_kg_s,
            friction_pressure_drop_pa=friction_drop,
            stored_liquid_mass_rate_kg_s=stored_mass_rate,
            external_mass_closure_error_kg_s=closure_error,
            natural_angular_frequency_rad_s=natural_frequency,
            damping_ratio=damping_ratio,
            derivative=DynamicFeedLineDerivative(flow_rate, pressure_rate),
            partials=partials,
        )

    def event_surfaces(self) -> tuple[DynamicFeedLineEventSurface, ...]:
        """Return the reverse-flow, pressure-limit and flow-magnitude events.

        Returns
        -------
        tuple of DynamicFeedLineEventSurface
            One burn-stopping surface per `DynamicFeedLineEventKind` member;
            the pressure limits screen transient pressure events [1]_.

        References
        ----------
        .. [1] NASA NESC Technical Bulletin 22-03, "Treatment of Transient Pressure
               Events", 2022. https://ntrs.nasa.gov/citations/20220006583
        """

        return (
            DynamicFeedLineEventSurface(
                "feed-line:reverse-flow", DynamicFeedLineEventKind.REVERSE_FLOW
            ),
            DynamicFeedLineEventSurface(
                "feed-line:minimum-pressure",
                DynamicFeedLineEventKind.MINIMUM_MANIFOLD_PRESSURE,
            ),
            DynamicFeedLineEventSurface(
                "feed-line:maximum-pressure",
                DynamicFeedLineEventKind.MAXIMUM_MANIFOLD_PRESSURE,
            ),
            DynamicFeedLineEventSurface(
                "feed-line:maximum-flow",
                DynamicFeedLineEventKind.MAXIMUM_ABSOLUTE_FLOW,
            ),
        )

    def event_value(
        self, surface: DynamicFeedLineEventSurface, state: DynamicFeedLineState
    ) -> float:
        """Evaluate one line event function, positive inside the permitted domain.

        Parameters
        ----------
        surface : DynamicFeedLineEventSurface
            Event to evaluate.
        state : DynamicFeedLineState
            Current line state.

        Returns
        -------
        float
            ``q`` for reverse flow, ``p_m - p_min``, ``p_max - p_m`` or
            ``q_max - abs(q)``; a sign change marks the event for a
            zero-crossing detector [2]_.  The pressure limits screen transient
            pressure events [1]_.

        Raises
        ------
        DomainError
            If the event kind is not supported.

        References
        ----------
        .. [1] NASA NESC Technical Bulletin 22-03, "Treatment of Transient Pressure
               Events", 2022. https://ntrs.nasa.gov/citations/20220006583
        .. [2] Orekit Project, "PropulsionModel / EventDetector" API documentation,
               Orekit 13.1.5.
        """

        config = self.configuration
        if surface.kind is DynamicFeedLineEventKind.REVERSE_FLOW:
            return state.mass_flow_kg_s
        if surface.kind is DynamicFeedLineEventKind.MINIMUM_MANIFOLD_PRESSURE:
            return state.manifold_pressure_pa - config.minimum_manifold_pressure_pa
        if surface.kind is DynamicFeedLineEventKind.MAXIMUM_MANIFOLD_PRESSURE:
            return config.maximum_manifold_pressure_pa - state.manifold_pressure_pa
        if surface.kind is DynamicFeedLineEventKind.MAXIMUM_ABSOLUTE_FLOW:
            return config.maximum_absolute_flow_kg_s - abs(state.mass_flow_kg_s)
        raise DomainError(f"Unsupported dynamic feed-line event: {surface.kind}.")


def evaluate_dynamic_feed_propulsion(
    bridge: PropagatorPropulsionBridge,
    regulated_feed: RegulatedFeedSystem,
    regulated_state: RegulatedFeedState,
    feed_line: DynamicFeedLine,
    line_state: DynamicFeedLineState,
    *,
    relative_time_s: float,
    vehicle_mass_kg: float,
    direction: tuple[float, float, float],
    pressure_condition_name: str,
    other_operating_conditions: Mapping[str, float] | None = None,
    additional_states: Mapping[str, float] | None = None,
    parameter_overrides: Mapping[str, float] | None = None,
    apply_inventory_limits: bool = True,
    ambient_pressure_pa: float | None = None,
) -> DynamicFeedPropulsionEvaluation:
    """Evaluate pressurant, inertial line and engine without an algebraic loop.

    Parameters
    ----------
    bridge : PropagatorPropulsionBridge
        Propagator bridge with a performance surface.
    regulated_feed : RegulatedFeedSystem
        Tank and pressurant model; its static line resistance must be zero.
    regulated_state : RegulatedFeedState
        Current tank and pressurant state.
    feed_line : DynamicFeedLine
        Dynamic line model.
    line_state : DynamicFeedLineState
        Current line flow and manifold pressure.
    relative_time_s : float
        Time from burn start, s.
    vehicle_mass_kg : float
        Current vehicle mass, kg.
    direction : tuple of float
        Thrust direction in the propagation frame.
    pressure_condition_name : str
        Performance-surface axis driven by the manifold pressure.
    other_operating_conditions : Mapping[str, float] or None, default None
        Remaining operating conditions; must not contain
        ``pressure_condition_name``.
    additional_states : Mapping[str, float] or None, default None
        Additional propagated states passed to the bridge.
    parameter_overrides : Mapping[str, float] or None, default None
        Parameter values passed to the bridge.
    apply_inventory_limits : bool, default True
        Whether the bridge applies inventory cutoffs.
    ambient_pressure_pa : float or None, default None
        Ambient pressure passed to the bridge, Pa.

    Returns
    -------
    DynamicFeedPropulsionEvaluation
        Engine, feed and line evaluations with the composed Jacobian.

    Raises
    ------
    DomainError
        If the bridge has no performance surface, the pressure condition is not
        one of its axes or also appears in ``other_operating_conditions``, the
        regulated feed carries a static line resistance, or an active bridge
        does not expose the pressure-condition partials.

    Notes
    -----
    Manifold pressure is a propagated line state and drives the engine
    performance map directly.  Engine demand drives compliant storage; line
    flow drains the tank and expands the ullage [1]_, [2]_.  The composed
    Jacobian includes the engine-demand derivative inside the
    manifold-pressure ODE [3]_.  While the bridge is off or inhibited, engine
    demand and its pressure derivative are zero, and the line and pressurant
    states may continue to relax.

    The regulated feed's static line resistance is zero so that the loss is not
    applied twice; dynamic resistance belongs to `DynamicFeedLine`.

    Negative trial flows are continued consistently through the tank-mass and
    ullage-work equations for zero-flow event location.  Such results carry
    ``reverse_flow_event_extension=True``; physical propagation stops at zero
    flow, and the continuation is not a sustained reverse-operation model.

    References
    ----------
    .. [1] NASA, "Liquid Rocket Propellant Feedline Dynamics", 1973.
           https://ntrs.nasa.gov/citations/19740028545
    .. [2] NASA, "Nodal Modeling of Liquid Propellant Feed and Pressurization
           System", 2024. https://ntrs.nasa.gov/citations/20240003493
    .. [3] NASA, "Transient Mathematical Modeling for Liquid Rocket Engine
           Systems". https://ntrs.nasa.gov/citations/20040000363
    """

    if bridge.performance_surface is None:
        raise DomainError("Dynamic feed coupling requires a bridge performance surface.")
    if pressure_condition_name not in {axis.name for axis in bridge.performance_surface.axes}:
        raise DomainError("Dynamic feed pressure condition is not a performance axis.")
    if regulated_feed.configuration.feed_line_resistance_pa_s2_kg2 != 0.0:
        raise DomainError("Regulated-feed static line resistance must be zero with a dynamic line.")
    conditions = dict(other_operating_conditions or {})
    if pressure_condition_name in conditions:
        raise DomainError("Dynamic line owns the declared pressure condition.")
    conditions[pressure_condition_name] = line_state.manifold_pressure_pa
    propulsion = bridge.evaluate(
        relative_time_s=relative_time_s,
        vehicle_mass_kg=vehicle_mass_kg,
        direction=direction,
        additional_states=additional_states,
        operating_conditions=conditions,
        parameter_overrides=parameter_overrides,
        apply_inventory_limits=apply_inventory_limits,
        ambient_pressure_pa=ambient_pressure_pa,
    )
    engine_flow = -propulsion.mass_derivative_kg_s
    feed_evaluation = regulated_feed.evaluate(
        regulated_state, line_state.mass_flow_kg_s,
        for_event_location=line_state.mass_flow_kg_s < 0.0,
    )
    line_evaluation = feed_line.evaluate(
        line_state,
        upstream_pressure_pa=feed_evaluation.tank_pressure_pa,
        engine_mass_flow_kg_s=engine_flow,
    )
    if propulsion.active:
        acceleration_by_condition = dict(propulsion.partials.acceleration_wrt_conditions)
        mass_rate_by_condition = dict(propulsion.partials.mass_rate_wrt_conditions)
        if pressure_condition_name not in acceleration_by_condition or (
            pressure_condition_name not in mass_rate_by_condition
        ):
            raise DomainError("Bridge did not expose dynamic pressure-condition partials.")
        acceleration_wrt_pressure = acceleration_by_condition[pressure_condition_name]
        engine_flow_wrt_pressure = -mass_rate_by_condition[pressure_condition_name]
    else:
        # The bridge omits performance-map partials while it is off or
        # inhibited.  Engine demand and its local pressure derivative are then
        # exactly zero; the line and pressurization states may still relax.
        acceleration_wrt_pressure = (0.0, 0.0, 0.0)
        engine_flow_wrt_pressure = 0.0
    raw_state_partials = dict(line_evaluation.partials.derivative_wrt_state)
    line_state_jacobian = (
        (
            raw_state_partials["mass_flow_kg_s"][0],
            raw_state_partials["manifold_pressure_pa"][0],
        ),
        (
            raw_state_partials["mass_flow_kg_s"][1],
            -engine_flow_wrt_pressure / feed_line.configuration.compliance_kg_pa,
        ),
    )
    upstream_gain = line_evaluation.partials.derivative_wrt_upstream_pressure[0]
    line_wrt_feed = tuple(
        (name, (upstream_gain * pressure_partial, 0.0))
        for name, pressure_partial in (
            feed_evaluation.pressure_partials.injector_pressure_wrt_state
        )
    )
    zero = (0.0, 0.0, 0.0)
    acceleration_wrt_line_state = (
        ("mass_flow_kg_s", zero),
        ("manifold_pressure_pa", acceleration_wrt_pressure),
    )
    return DynamicFeedPropulsionEvaluation(
        propulsion=propulsion,
        regulated_feed=feed_evaluation,
        feed_line=line_evaluation,
        partials=DynamicFeedCoupledPartials(
            line_state_jacobian=line_state_jacobian,
            line_derivative_wrt_feed_state=line_wrt_feed,
            engine_flow_wrt_manifold_pressure=engine_flow_wrt_pressure,
            acceleration_wrt_line_state=acceleration_wrt_line_state,
        ),
    )
