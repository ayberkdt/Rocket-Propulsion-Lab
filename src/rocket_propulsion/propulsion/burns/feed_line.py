"""Distributed-effect surrogate for liquid feed-line transients.

The two-state lumped line retains hydraulic inertance and compliance rather
than collapsing every stage evaluation to a static pressure loss.  It is small
enough for orbit propagation while exposing the conservation and Jacobian
terms required by a numerical force-model implementation.
"""

from __future__ import annotations

import hashlib
import json
from collections.abc import Mapping
from dataclasses import dataclass
from enum import Enum
from math import isfinite, sqrt
from typing import TYPE_CHECKING

from rocket_propulsion.core.errors import DomainError

from .feed_system import RegulatedFeedEvaluation, RegulatedFeedState, RegulatedFeedSystem

if TYPE_CHECKING:
    from .propagator import PropagatorPropulsionBridge, PropulsionDynamicsEvaluation

DYNAMIC_FEED_LINE_SCHEMA = "rocket_propulsion_dynamic_feed_line_v1"
DYNAMIC_FEED_LINE_REFERENCE_IDS = (
    "NASA-19740028545",
    "NASA-20220006583",
    "NASA-20240003493",
)


def _positive(value: float, label: str) -> None:
    if not isfinite(value) or value <= 0.0:
        raise DomainError(f"{label} must be finite and greater than zero.")


def _nonnegative(value: float, label: str) -> None:
    if not isfinite(value) or value < 0.0:
        raise DomainError(f"{label} must be finite and non-negative.")


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
class DynamicFeedLineConfiguration:
    """Lumped hydraulic inertance, compliance, resistance, and limits.

    The governing equations are

    ``L dq/dt = p_up - p_m - R q abs(q)``

    and

    ``C dp_m/dt = q - q_engine``.

    ``C`` therefore has units kg/Pa and represents incremental stored liquid
    mass per pressure.  The model is a low-order surrogate for distributed line
    modes and must be calibrated against the actual plumbing or a higher-order
    nodal model.

    References
    ----------
    NASA liquid-rocket feedline dynamics:
    https://ntrs.nasa.gov/citations/19740028545
    NASA NESC transient-pressure bulletin:
    https://ntrs.nasa.gov/citations/20220006583
    NASA nodal feed/pressurization modeling:
    https://ntrs.nasa.gov/citations/20240003493
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
    """Propagated line flow and engine-manifold pressure.

    Negative flow is retained so an event locator can resolve a flow reversal;
    the configured reverse-flow event is positive-safe and stops the burn at
    zero before sustained reverse operation.

    References
    ----------
    NASA liquid-rocket feedline dynamics:
    https://ntrs.nasa.gov/citations/19740028545
    NASA NESC transient-pressure bulletin:
    https://ntrs.nasa.gov/citations/20220006583
    """

    mass_flow_kg_s: float
    manifold_pressure_pa: float

    def __post_init__(self) -> None:
        if not isfinite(self.mass_flow_kg_s):
            raise DomainError("Feed-line mass flow must be finite.")
        _positive(self.manifold_pressure_pa, "Feed-line manifold pressure")


@dataclass(frozen=True, slots=True)
class DynamicFeedLineDerivative:
    """Two-state feed-line ODE right-hand side.

    References
    ----------
    NASA liquid-rocket feedline dynamics:
    https://ntrs.nasa.gov/citations/19740028545
    """

    mass_flow_kg_s2: float
    manifold_pressure_pa_s: float

    def as_tuple(self) -> tuple[float, float]:
        """Return derivatives in state-field order.

        References
        ----------
        NASA liquid-rocket feedline dynamics:
        https://ntrs.nasa.gov/citations/19740028545
        """

        return (self.mass_flow_kg_s2, self.manifold_pressure_pa_s)


@dataclass(frozen=True, slots=True)
class DynamicFeedLinePartials:
    """Analytic local line derivatives before engine-map composition.

    References
    ----------
    NASA liquid-rocket feedline dynamics:
    https://ntrs.nasa.gov/citations/19740028545
    """

    derivative_wrt_state: tuple[tuple[str, tuple[float, float]], ...]
    derivative_wrt_upstream_pressure: tuple[float, float]
    derivative_wrt_engine_flow: tuple[float, float]
    nonsmooth_state_names: tuple[str, ...] = ()


@dataclass(frozen=True, slots=True)
class DynamicFeedLineEvaluation:
    """Line ODE, conservation terms, and local modal diagnostics.

    References
    ----------
    NASA liquid-rocket feedline dynamics:
    https://ntrs.nasa.gov/citations/19740028545
    NASA NESC transient-pressure bulletin:
    https://ntrs.nasa.gov/citations/20220006583
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


class DynamicFeedLineEventKind(str, Enum):
    """Positive-safe cutoff events for the dynamic line.

    References
    ----------
    NASA NESC transient-pressure bulletin:
    https://ntrs.nasa.gov/citations/20220006583
    """

    REVERSE_FLOW = "reverse_flow"
    MINIMUM_MANIFOLD_PRESSURE = "minimum_manifold_pressure"
    MAXIMUM_MANIFOLD_PRESSURE = "maximum_manifold_pressure"
    MAXIMUM_ABSOLUTE_FLOW = "maximum_absolute_flow"


@dataclass(frozen=True, slots=True)
class DynamicFeedLineEventSurface:
    """One signed root that is positive inside the permitted line domain.

    References
    ----------
    NASA NESC transient-pressure bulletin:
    https://ntrs.nasa.gov/citations/20220006583
    Orekit EventDetector 13.1.5 API:
    https://www.orekit.org/static/apidocs/org/orekit/propagation/events/EventDetector.html
    """

    event_id: str
    kind: DynamicFeedLineEventKind
    stops_burn: bool = True


@dataclass(frozen=True, slots=True)
class DynamicFeedCoupledPartials:
    """Engine-composed line Jacobian and feed-state coupling.

    References
    ----------
    NASA transient-modeling overview:
    https://ntrs.nasa.gov/citations/20040000363
    NASA liquid-rocket feedline dynamics:
    https://ntrs.nasa.gov/citations/19740028545
    """

    line_state_jacobian: tuple[tuple[float, float], tuple[float, float]]
    line_derivative_wrt_feed_state: tuple[tuple[str, tuple[float, float]], ...]
    engine_flow_wrt_manifold_pressure: float
    acceleration_wrt_line_state: tuple[tuple[str, tuple[float, float, float]], ...]


@dataclass(frozen=True, slots=True)
class DynamicFeedPropulsionEvaluation:
    """Simultaneous pressurant, line, engine, and Jacobian evaluation.

    References
    ----------
    NASA nodal feed/pressurization modeling:
    https://ntrs.nasa.gov/citations/20240003493
    NASA transient-modeling overview:
    https://ntrs.nasa.gov/citations/20040000363
    """

    propulsion: PropulsionDynamicsEvaluation
    regulated_feed: RegulatedFeedEvaluation
    feed_line: DynamicFeedLineEvaluation
    partials: DynamicFeedCoupledPartials


@dataclass(frozen=True, slots=True)
class DynamicFeedLine:
    """Two-state liquid-line model for a numerical orbit propagator.

    The state captures the first inertance/compliance mode and nonlinear
    quadratic damping.  It is appropriate for mission-level transient coupling
    only after its effective coefficients have been correlated against the
    actual line or a validated nodal model.  It is not a water-hammer solver.

    References
    ----------
    NASA liquid-rocket feedline dynamics:
    https://ntrs.nasa.gov/citations/19740028545
    NASA NESC transient-pressure bulletin:
    https://ntrs.nasa.gov/citations/20220006583
    NASA nodal feed/pressurization modeling:
    https://ntrs.nasa.gov/citations/20240003493
    """

    configuration: DynamicFeedLineConfiguration

    @property
    def model_hash(self) -> str:
        """Return configuration and reference identity.

        References
        ----------
        NASA-STD-7009B: https://standards.nasa.gov/standard/NASA/NASA-STD-7009
        """

        return _canonical_hash(self.manifest(include_hash=False))

    def manifest(self, *, include_hash: bool = True) -> dict[str, object]:
        """Return a deterministic, unit-explicit model manifest.

        References
        ----------
        NASA-STD-7009B: https://standards.nasa.gov/standard/NASA/NASA-STD-7009
        NASA liquid-rocket feedline dynamics:
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
            payload["model_hash"] = _canonical_hash(payload)
        return payload

    def evaluate(
        self,
        state: DynamicFeedLineState,
        *,
        upstream_pressure_pa: float,
        engine_mass_flow_kg_s: float,
    ) -> DynamicFeedLineEvaluation:
        """Evaluate line dynamics, conservation, modes, and partials.

        References
        ----------
        NASA liquid-rocket feedline dynamics:
        https://ntrs.nasa.gov/citations/19740028545
        NASA nodal feed/pressurization modeling:
        https://ntrs.nasa.gov/citations/20240003493
        """

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
        natural_frequency = sqrt(1.0 / (config.inertance_pa_s2_kg * config.compliance_kg_pa))
        damping_ratio = (
            config.resistance_pa_s2_kg2
            * abs(flow)
            * sqrt(config.compliance_kg_pa / config.inertance_pa_s2_kg)
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
        # outflow.  This residual is explicit so an adapter can monitor it.
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
        """Return reverse-flow, pressure, and magnitude limit roots.

        References
        ----------
        NASA NESC transient-pressure bulletin:
        https://ntrs.nasa.gov/citations/20220006583
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
        """Evaluate one positive-safe line event function.

        References
        ----------
        NASA NESC transient-pressure bulletin:
        https://ntrs.nasa.gov/citations/20220006583
        Orekit EventDetector 13.1.5 API:
        https://www.orekit.org/static/apidocs/org/orekit/propagation/events/EventDetector.html
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
) -> DynamicFeedPropulsionEvaluation:
    """Evaluate pressurant, inertial line, and engine without an algebraic loop.

    Manifold pressure is a propagated line state and directly drives the engine
    performance map.  Engine demand drives compliant storage; line flow drains
    the tank and expands the ullage.  The returned composed Jacobian includes
    the engine demand derivative inside the manifold-pressure ODE.

    The regulated feed's static line resistance must be zero to prevent double
    application; dynamic resistance belongs to `DynamicFeedLine`.

    References
    ----------
    NASA liquid-rocket feedline dynamics:
    https://ntrs.nasa.gov/citations/19740028545
    NASA nodal feed/pressurization modeling:
    https://ntrs.nasa.gov/citations/20240003493
    NASA transient-modeling overview:
    https://ntrs.nasa.gov/citations/20040000363
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
    )
    engine_flow = -propulsion.mass_derivative_kg_s
    feed_evaluation = regulated_feed.evaluate(regulated_state, max(0.0, line_state.mass_flow_kg_s))
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
        # The bridge deliberately omits performance-map partials while it is
        # off or inhibited.  The physical engine demand and its local pressure
        # derivative are both exactly zero in that state; the line and
        # pressurization states may nevertheless continue to relax.
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
