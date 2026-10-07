"""Propagation-ready force, mass-flow, event, and sensitivity contract.

This module deliberately stops at the propulsion/astrodynamics boundary.  It
does not own epochs, frames, attitude laws, or orbit integration.  A numerical
propagator supplies relative time, current mass, a resolved thrust direction,
and any tracked tank states; the bridge returns simultaneous translational and
mass derivatives plus the event surfaces that make discontinuities explicit.
"""

from __future__ import annotations

import hashlib
import json
from collections.abc import Mapping
from dataclasses import dataclass
from enum import Enum
from math import isfinite, sqrt
from types import MappingProxyType

from rocket_propulsion.core.errors import DomainError

from .performance_surface import (
    PerformanceValiditySurface,
    RectilinearPerformanceSurface,
)
from .tabulated import InterpolationPolicy, TabulatedBurnArtifact, TabulatedBurnSegment

PROPAGATOR_BRIDGE_SCHEMA = "rocket_propulsion_propagator_bridge_v1"
PROPAGATOR_BRIDGE_REFERENCE_IDS = (
    "OREKIT-13.1.5-PROPULSION-MODEL",
    "NASA-20040000363",
    "CCSDS-502.0-B-3",
)


class BurnEventKind(str, Enum):
    """Kinds of root surfaces exposed to a numerical propagator.

    References
    ----------
    Orekit EventDetectorsProvider 13.1.5 API:
    https://www.orekit.org/static/apidocs/org/orekit/propagation/events/EventDetectorsProvider.html
    """

    START = "start"
    PROFILE_BOUNDARY = "profile_boundary"
    STOP = "stop"
    DRY_MASS = "dry_mass"
    TANK_RESERVE = "tank_reserve"


@dataclass(frozen=True, slots=True)
class EstimableParameter:
    """One bounded parameter-driver definition for orbit estimation.

    ``scale`` is a recommended finite-difference/normalization scale; it does
    not alter the physical value.  The bridge uses immutable parameter names so
    an estimator cannot silently reorder a positional parameter vector.

    References
    ----------
    Orekit ParameterDriversProvider 13.1.5 API:
    https://www.orekit.org/static/apidocs/org/orekit/utils/ParameterDriversProvider.html
    """

    name: str
    reference_value: float
    scale: float
    minimum: float
    maximum: float
    unit: str

    def __post_init__(self) -> None:
        if not self.name.strip():
            raise DomainError("Estimable parameter name must not be empty.")
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
        """Return a finite in-bounds value or fail closed."""

        try:
            numeric = float(value)
        except (TypeError, ValueError) as error:
            raise DomainError(f"Parameter {self.name!r} must be numeric.") from error
        if not isfinite(numeric) or not self.minimum <= numeric <= self.maximum:
            raise DomainError(
                f"Parameter {self.name!r} must be finite and in "
                f"[{self.minimum}, {self.maximum}]."
            )
        return numeric


@dataclass(frozen=True, slots=True)
class StreamStateBinding:
    """Map one artifact stream to one propagator additional-mass state.

    Several streams may drain the same state; their derivatives are summed.

    References
    ----------
    NASA-SP-8112: https://ntrs.nasa.gov/citations/19760015212
    """

    stream_id: str
    state_name: str

    def __post_init__(self) -> None:
        if not self.stream_id.strip() or not self.state_name.strip():
            raise DomainError("Stream and additional-state names must not be empty.")


@dataclass(frozen=True, slots=True)
class AdditionalStateConstraint:
    """Minimum permitted value for a tracked tank or propellant state.

    References
    ----------
    NASA-SP-8112: https://ntrs.nasa.gov/citations/19760015212
    """

    state_name: str
    minimum_value: float
    label: str = "tank reserve"

    def __post_init__(self) -> None:
        if not self.state_name.strip() or not self.label.strip():
            raise DomainError("State-constraint names must not be empty.")
        if not isfinite(self.minimum_value) or self.minimum_value < 0.0:
            raise DomainError("State-constraint minimum must be finite and non-negative.")


@dataclass(frozen=True, slots=True)
class BurnEventSurface:
    """A signed root function with explicit propulsion-stop semantics.

    Time surfaces use ``g = t - (t_event + ignition_time_bias)``.  Mass
    surfaces use ``g = current - minimum``.  A zero therefore means the
    propagator must land exactly on the discontinuity and recompute the right
    hand derivative.  ``stops_burn`` does not request termination of the orbit
    propagation; it only inhibits this propulsion model.

    References
    ----------
    Orekit EventDetector 13.1.5 API:
    https://www.orekit.org/static/apidocs/org/orekit/propagation/events/EventDetector.html
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
        """Evaluate the signed event function without hidden state."""

        if not isfinite(relative_time_s) or not isfinite(vehicle_mass_kg):
            raise DomainError("Event evaluation time and mass must be finite.")
        if not isfinite(ignition_time_bias_s):
            raise DomainError("Event ignition-time bias must be finite.")
        if self.relative_time_s is not None:
            return relative_time_s - (self.relative_time_s + ignition_time_bias_s)
        if self.kind is BurnEventKind.DRY_MASS:
            assert self.threshold is not None
            return vehicle_mass_kg - self.threshold
        if self.state_name is None or self.threshold is None:
            raise DomainError("Malformed state event surface.")
        if self.state_name not in additional_states:
            raise DomainError(f"Missing additional state {self.state_name!r}.")
        value = float(additional_states[self.state_name])
        if not isfinite(value):
            raise DomainError(f"Additional state {self.state_name!r} must be finite.")
        return value - self.threshold


@dataclass(frozen=True, slots=True)
class PropulsionPartials:
    """Analytic local derivatives for variational equations and estimation.

    Acceleration partials are expressed in the same resolved frame as the
    supplied direction.  At a command discontinuity the ignition-time partial
    is intentionally reported as non-smooth instead of inventing a derivative.

    References
    ----------
    Orekit PropulsionModel 13.1.5 API:
    https://www.orekit.org/static/apidocs/org/orekit/forces/maneuvers/propulsion/PropulsionModel.html
    """

    acceleration_wrt_mass: tuple[float, float, float]
    acceleration_wrt_parameters: tuple[
        tuple[str, tuple[float, float, float]], ...
    ]
    mass_rate_wrt_parameters: tuple[tuple[str, float], ...]
    additional_state_rates_wrt_parameters: tuple[
        tuple[str, tuple[tuple[str, float], ...]], ...
    ]
    acceleration_wrt_conditions: tuple[
        tuple[str, tuple[float, float, float]], ...
    ] = ()
    mass_rate_wrt_conditions: tuple[tuple[str, float], ...] = ()
    additional_state_rates_wrt_conditions: tuple[
        tuple[str, tuple[tuple[str, float], ...]], ...
    ] = ()
    nonsmooth_parameters: tuple[str, ...] = ()


@dataclass(frozen=True, slots=True)
class PropulsionDynamicsEvaluation:
    """Simultaneous force and inventory derivative at one propagator call.

    References
    ----------
    NASA-20040000363: https://ntrs.nasa.gov/citations/20040000363
    Orekit PropulsionModel 13.1.5 API:
    https://www.orekit.org/static/apidocs/org/orekit/forces/maneuvers/propulsion/PropulsionModel.html
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


def _normalized(vector: tuple[float, float, float]) -> tuple[float, float, float]:
    if len(vector) != 3 or any(not isfinite(value) for value in vector):
        raise DomainError("Thrust direction must contain three finite components.")
    norm = sqrt(sum(value * value for value in vector))
    if norm <= 0.0:
        raise DomainError("Thrust direction must be non-zero.")
    return tuple(value / norm for value in vector)  # type: ignore[return-value]


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
class PropagatorPropulsionBridge:
    """Immutable tabulated propulsion model for numerical orbit propagation.

    The returned force is the artifact's *axial* force: cant/divergence already
    represented by the producer is not applied twice.  Vehicle and per-stream
    mass derivatives are negative drains.  ``thrust_scale`` and
    ``mass_flow_scale`` are independent so orbit determination can estimate
    force and propellant bookkeeping without fabricating an Isp constraint;
    upstream calibrated realizations should normally leave both at one.

    ``ignition_time_bias_s`` shifts the entire profile and all date events by
    the same amount.  Its derivative is analytic inside linear segments and is
    marked non-smooth at every declared boundary.

    References
    ----------
    Orekit PropulsionModel 13.1.5 API:
    https://www.orekit.org/static/apidocs/org/orekit/forces/maneuvers/propulsion/PropulsionModel.html
    Orekit Maneuver 13.1.5 API:
    https://www.orekit.org/static/apidocs/org/orekit/forces/maneuvers/Maneuver.html
    NASA-20040000363: https://ntrs.nasa.gov/citations/20040000363
    """

    artifact: TabulatedBurnArtifact
    performance_surface: RectilinearPerformanceSurface | None = None
    stream_bindings: tuple[StreamStateBinding, ...] = ()
    state_constraints: tuple[AdditionalStateConstraint, ...] = ()
    protected_dry_mass_kg: float | None = None
    require_complete_stream_binding: bool = True
    parameters: tuple[EstimableParameter, ...] = ()
    model_id: str = "tabulated-propulsion-dynamics"
    schema: str = PROPAGATOR_BRIDGE_SCHEMA

    def __post_init__(self) -> None:
        if self.schema != PROPAGATOR_BRIDGE_SCHEMA:
            raise DomainError("Unsupported propagator bridge schema.")
        if not self.model_id.strip():
            raise DomainError("Propagator bridge model identifier must not be empty.")
        # Re-open at the public integrity boundary before accepting the source.
        TabulatedBurnArtifact.from_dict(self.artifact.to_dict())
        if self.protected_dry_mass_kg is not None and (
            not isfinite(self.protected_dry_mass_kg) or self.protected_dry_mass_kg <= 0.0
        ):
            raise DomainError("Protected dry mass must be finite and greater than zero.")

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
                segment.start_stream_mass_flows_kg_s
                + segment.end_stream_mass_flows_kg_s
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

        if not self.parameters:
            limit = max(1.0, self.artifact.duration_s)
            object.__setattr__(
                self,
                "parameters",
                (
                    EstimableParameter("thrust_scale", 1.0, 1e-3, 0.01, 2.0, "1"),
                    EstimableParameter("mass_flow_scale", 1.0, 1e-3, 0.01, 2.0, "1"),
                    EstimableParameter(
                        "ignition_time_bias_s", 0.0, 1e-3, -limit, limit, "s"
                    ),
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
        """Content identity binding source artifact and propagation semantics."""

        return _canonical_hash(self.manifest(include_hash=False))

    def manifest(self, *, include_hash: bool = True) -> dict[str, object]:
        """Return a deterministic, integration-reviewable model manifest."""

        payload: dict[str, object] = {
            "schema": self.schema,
            "model_id": self.model_id,
            "artifact_hash": self.artifact.artifact_hash,
            "performance_surface_hash": (
                None
                if self.performance_surface is None
                else self.performance_surface.surface_hash
            ),
            "operating_condition_axes": (
                []
                if self.performance_surface is None
                else [
                    {
                        "name": axis.name,
                        "unit": axis.unit,
                        "minimum": axis.minimum,
                        "maximum": axis.maximum,
                    }
                    for axis in self.performance_surface.axes
                ]
            ),
            "force_semantics": "net_axial_thrust_along_resolved_direction",
            "mass_semantics": "negative_total_and_bound_stream_tank_drain",
            "time_semantics": "relative_seconds_with_global_ignition_bias",
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
        if include_hash:
            payload["model_hash"] = _canonical_hash(payload)
        return payload

    def resolved_parameters(
        self, overrides: Mapping[str, float] | None = None
    ) -> Mapping[str, float]:
        """Resolve named overrides against immutable driver definitions."""

        supplied = dict(overrides or {})
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
        """Return every time and inventory discontinuity exactly once."""

        surfaces = [
            BurnEventSurface("burn:start", BurnEventKind.START, False, relative_time_s=0.0)
        ]
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
        """Return state-domain roots for the optional performance surface.

        References
        ----------
        NASA-STD-7009B: https://standards.nasa.gov/standard/NASA/NASA-STD-7009
        Orekit EventDetectorsProvider 13.1.5 API:
        https://www.orekit.org/static/apidocs/org/orekit/propagation/events/EventDetectorsProvider.html
        """

        if self.performance_surface is None:
            return ()
        return self.performance_surface.validity_surfaces()

    def _segment_at(
        self, effective_time_s: float, side: str
    ) -> tuple[TabulatedBurnSegment | None, bool]:
        tolerance = max(1e-12, self.artifact.duration_s * 1e-12)
        at_boundary = any(
            abs(effective_time_s - value) <= tolerance
            for value in (
                0.0,
                *(segment.end_time_s for segment in self.artifact.segments),
            )
        )
        if effective_time_s < 0.0 or effective_time_s > self.artifact.duration_s:
            return None, at_boundary
        if abs(effective_time_s) <= tolerance and side == "left":
            return None, True
        if abs(effective_time_s - self.artifact.duration_s) <= tolerance and side == "right":
            return None, True
        for index, segment in enumerate(self.artifact.segments):
            if segment.start_time_s <= effective_time_s < segment.end_time_s:
                if (
                    side == "left"
                    and index > 0
                    and abs(effective_time_s - segment.start_time_s) <= tolerance
                ):
                    return self.artifact.segments[index - 1], True
                return segment, at_boundary
        return self.artifact.segments[-1], at_boundary

    def evaluate(
        self,
        *,
        relative_time_s: float,
        vehicle_mass_kg: float,
        direction: tuple[float, float, float],
        additional_states: Mapping[str, float] | None = None,
        operating_conditions: Mapping[str, float] | None = None,
        parameter_overrides: Mapping[str, float] | None = None,
        side: str = "right",
    ) -> PropulsionDynamicsEvaluation:
        """Evaluate force, drains, event-safe activity, and analytic partials.

        ``side`` selects the limiting value exactly at discontinuities.  Normal
        integration should use ``right`` after an event reset; event localization
        may query both sides.  Outside the shifted firing interval the method
        returns exact zeros rather than extrapolating the profile.

        References
        ----------
        Orekit PropulsionModel 13.1.5 API:
        https://www.orekit.org/static/apidocs/org/orekit/forces/maneuvers/propulsion/PropulsionModel.html
        NASA-20040000363: https://ntrs.nasa.gov/citations/20040000363
        """

        if side not in {"left", "right"}:
            raise DomainError("Propagation boundary side must be 'left' or 'right'.")
        if not isfinite(relative_time_s):
            raise DomainError("Relative propagation time must be finite.")
        if not isfinite(vehicle_mass_kg) or vehicle_mass_kg <= 0.0:
            raise DomainError("Vehicle mass must be finite and greater than zero.")
        unit_direction = _normalized(direction)
        states = dict(additional_states or {})
        for name, value in states.items():
            try:
                numeric_value = float(value)
            except (TypeError, ValueError) as error:
                raise DomainError(
                    "Additional-state names and values must be valid and finite."
                ) from error
            if not isinstance(name, str) or not name.strip() or not isfinite(numeric_value):
                raise DomainError("Additional-state names and values must be valid and finite.")
            states[name] = numeric_value
        required_states = {binding.state_name for binding in self.stream_bindings}
        missing_states = required_states - set(states)
        if missing_states:
            raise DomainError(
                "Missing bound propulsion states: " + ", ".join(sorted(missing_states))
            )
        supplied_conditions = dict(operating_conditions or {})
        if self.performance_surface is None and supplied_conditions:
            raise DomainError(
                "Operating conditions were supplied but this bridge has no performance surface."
            )

        parameter_values = self.resolved_parameters(parameter_overrides)
        thrust_scale = parameter_values["thrust_scale"]
        flow_scale = parameter_values["mass_flow_scale"]
        time_bias = parameter_values["ignition_time_bias_s"]
        effective_time = relative_time_s - time_bias
        segment, at_boundary = self._segment_at(effective_time, side)

        inhibited: list[str] = []
        if self.protected_dry_mass_kg is not None and vehicle_mass_kg <= (
            self.protected_dry_mass_kg + 1e-12
        ):
            inhibited.append("protected_dry_mass_kg")
        for constraint in self.state_constraints:
            if states[constraint.state_name] <= constraint.minimum_value + 1e-12:
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
                additional_state_derivatives=tuple(
                    (name, 0.0) for name in sorted(required_states)
                ),
                parameters=tuple(parameter_values.items()),
                inhibited_by=tuple(inhibited),
                partials=PropulsionPartials(
                    acceleration_wrt_mass=zero,
                    acceleration_wrt_parameters=tuple(
                        (name, zero) for name in parameter_values
                    ),
                    mass_rate_wrt_parameters=tuple(
                        (name, 0.0) for name in parameter_values
                    ),
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
        surface_force_scale = (
            1.0 if surface_evaluation is None else surface_evaluation.force_scale
        )
        surface_flow_scale = (
            1.0 if surface_evaluation is None else surface_evaluation.mass_flow_scale
        )
        sample = segment.state_at(
            min(max(effective_time, segment.start_time_s), segment.end_time_s)
        )
        base_thrust = float(sample["axial_thrust_n"])
        base_flow = float(sample["total_mass_flow_kg_s"])
        stream_flows = dict(sample["stream_mass_flows_kg_s"])
        thrust = base_thrust * thrust_scale * surface_force_scale
        total_flow = base_flow * flow_scale * surface_flow_scale
        force = tuple(thrust * value for value in unit_direction)
        acceleration = tuple(value / vehicle_mass_kg for value in force)

        binding_by_stream = {
            binding.stream_id: binding.state_name for binding in self.stream_bindings
        }
        state_drains = {state_name: 0.0 for state_name in required_states}
        for stream_id, flow in stream_flows.items():
            if stream_id in binding_by_stream:
                state_drains[binding_by_stream[stream_id]] -= (
                    flow * flow_scale * surface_flow_scale
                )

        duration = segment.duration_s
        if self.artifact.interpolation is InterpolationPolicy.HOLD:
            thrust_slope = 0.0
            flow_slope = 0.0
            stream_slopes = {stream_id: 0.0 for stream_id in stream_flows}
        else:
            thrust_slope = (
                segment.end_axial_thrust_n - segment.start_axial_thrust_n
            ) / duration
            flow_slope = (
                segment.end_total_mass_flow_kg_s - segment.start_total_mass_flow_kg_s
            ) / duration
            start_streams = dict(segment.start_stream_mass_flows_kg_s)
            end_streams = dict(segment.end_stream_mass_flows_kg_s)
            stream_slopes = {
                stream_id: (
                    end_streams.get(stream_id, 0.0) - start_streams.get(stream_id, 0.0)
                )
                / duration
                for stream_id in start_streams.keys() | end_streams.keys()
            }

        zero = (0.0, 0.0, 0.0)
        thrust_scale_partial = tuple(
            base_thrust * surface_force_scale * value / vehicle_mass_kg
            for value in unit_direction
        )
        time_partial = tuple(
            -thrust_scale
            * surface_force_scale
            * thrust_slope
            * value
            / vehicle_mass_kg
            for value in unit_direction
        )
        mass_partials = (
            ("thrust_scale", 0.0),
            ("mass_flow_scale", -base_flow * surface_flow_scale),
            (
                "ignition_time_bias_s",
                flow_scale * surface_flow_scale * flow_slope,
            ),
        )

        def state_parameter_partials(parameter_name: str) -> tuple[tuple[str, float], ...]:
            values = {state_name: 0.0 for state_name in required_states}
            if parameter_name == "mass_flow_scale":
                for stream_id, flow in stream_flows.items():
                    if stream_id in binding_by_stream:
                        values[binding_by_stream[stream_id]] -= flow * surface_flow_scale
            elif parameter_name == "ignition_time_bias_s":
                for stream_id, slope in stream_slopes.items():
                    if stream_id in binding_by_stream:
                        values[binding_by_stream[stream_id]] += (
                            flow_scale * surface_flow_scale * slope
                        )
            return tuple(sorted(values.items()))

        force_condition_gradients = (
            {}
            if surface_evaluation is None
            else dict(surface_evaluation.force_scale_gradient)
        )
        flow_condition_gradients = (
            {}
            if surface_evaluation is None
            else dict(surface_evaluation.mass_flow_scale_gradient)
        )
        acceleration_condition_partials = tuple(
            (
                name,
                tuple(
                    base_thrust
                    * thrust_scale
                    * gradient
                    * value
                    / vehicle_mass_kg
                    for value in unit_direction
                ),
            )
            for name, gradient in force_condition_gradients.items()
        )
        mass_condition_partials = tuple(
            (name, -base_flow * flow_scale * gradient)
            for name, gradient in flow_condition_gradients.items()
        )

        def state_condition_partials(condition_name: str) -> tuple[tuple[str, float], ...]:
            values = {state_name: 0.0 for state_name in required_states}
            gradient = flow_condition_gradients[condition_name]
            for stream_id, flow in stream_flows.items():
                if stream_id in binding_by_stream:
                    values[binding_by_stream[stream_id]] -= flow * flow_scale * gradient
            return tuple(sorted(values.items()))

        return PropulsionDynamicsEvaluation(
            relative_time_s=relative_time_s,
            effective_profile_time_s=effective_time,
            active=True,
            phase=str(sample["phase"]),
            axial_thrust_n=thrust,
            force_vector_n=force,  # type: ignore[arg-type]
            acceleration_m_s2=acceleration,  # type: ignore[arg-type]
            mass_derivative_kg_s=-total_flow,
            additional_state_derivatives=tuple(sorted(state_drains.items())),
            parameters=tuple(parameter_values.items()),
            inhibited_by=(),
            partials=PropulsionPartials(
                acceleration_wrt_mass=tuple(
                    -value / vehicle_mass_kg for value in acceleration
                ),  # type: ignore[arg-type]
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
                    (name, state_condition_partials(name))
                    for name in flow_condition_gradients
                ),
                nonsmooth_parameters=("ignition_time_bias_s",) if at_boundary else (),
            ),
            operating_conditions=(
                () if surface_evaluation is None else surface_evaluation.conditions
            ),
            performance_force_scale=surface_force_scale,
            performance_mass_flow_scale=surface_flow_scale,
        )

