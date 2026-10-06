"""Immutable public contracts for standalone propulsion-burn analysis.

The records in this module intentionally contain no position, velocity, orbit,
reference-frame, or guidance quantities.  They describe propulsion hardware,
propellant inventory, a burn target, and the resulting one-dimensional
propulsion bookkeeping only.
"""

from __future__ import annotations

from dataclasses import dataclass
from enum import Enum
from math import isfinite

from rocket_propulsion.core.errors import DomainError


class PropellantRole(str, Enum):
    """Functional role of a tank-depleting propellant stream."""

    FUEL = "fuel"
    OXIDIZER = "oxidizer"
    MONOPROPELLANT = "monopropellant"
    OTHER = "other"


class FlowDestination(str, Enum):
    """Destination of a propellant stream after it leaves its tank."""

    CHAMBER = "chamber"
    GAS_GENERATOR_DUMP = "gas_generator_dump"
    COOLANT_DUMP = "coolant_dump"
    VENT = "vent"
    OTHER = "other"


class BurnTargetKind(str, Enum):
    """Supported primary termination targets for the constant-burn solver."""

    DURATION = "duration"
    PROPELLANT_MASS = "propellant_mass"
    TOTAL_IMPULSE = "total_impulse"
    IDEAL_DELTA_V = "ideal_delta_v"


class CutoffReason(str, Enum):
    """Why a burn simulation stopped."""

    TARGET_REACHED = "target_reached"
    TANK_RESERVE = "tank_reserve"
    DRY_MASS_FLOOR = "dry_mass_floor"
    MAXIMUM_DURATION = "maximum_duration"
    SCHEDULE_END = "schedule_end"


def _finite_nonnegative(value: float, label: str) -> None:
    if not isfinite(value) or value < 0.0:
        raise DomainError(f"{label} must be finite and non-negative.")


def _finite_positive(value: float, label: str) -> None:
    if not isfinite(value) or value <= 0.0:
        raise DomainError(f"{label} must be finite and greater than zero.")


@dataclass(frozen=True, slots=True)
class PropellantFlow:
    """One explicit mass-flow path leaving a named tank."""

    stream_id: str
    tank_id: str
    propellant_key: str
    role: PropellantRole
    destination: FlowDestination
    mass_flow_kg_s: float
    tank_depleting: bool = True
    contributes_to_delivered_thrust: bool = True

    def __post_init__(self) -> None:
        if not self.stream_id.strip():
            raise DomainError("Stream identifier must not be empty.")
        if self.tank_depleting and not self.tank_id.strip():
            raise DomainError("A tank-depleting stream requires a tank identifier.")
        if not self.propellant_key.strip():
            raise DomainError("Propellant key must not be empty.")
        _finite_nonnegative(self.mass_flow_kg_s, "Stream mass flow")


@dataclass(frozen=True, slots=True)
class EngineOperatingPoint:
    """Delivered steady performance with explicit tank-depleting streams.

    Thrust and flows are totals for ``active_engine_count``.  Thus downstream
    solvers never multiply them by the engine count a second time.
    """

    name: str
    commanded_throttle: float
    realized_throttle: float
    active_engine_count: int
    ideal_thrust_n: float
    delivered_thrust_n: float
    chamber_specific_impulse_s: float
    system_specific_impulse_s: float
    chamber_mass_flow_kg_s: float
    total_tank_flow_kg_s: float
    oxidizer_mass_flow_kg_s: float
    fuel_mass_flow_kg_s: float
    mixture_ratio: float | None
    chamber_pressure_pa: float
    ambient_pressure_pa: float
    streams: tuple[PropellantFlow, ...]
    model_id: str = "constant_performance_v1"
    source: str = "user"
    warnings: tuple[str, ...] = ()

    def __post_init__(self) -> None:
        from rocket_propulsion.propulsion.performance import STANDARD_GRAVITY_M_S2

        if not self.name.strip():
            raise DomainError("Operating-point name must not be empty.")
        for value, label in (
            (self.commanded_throttle, "Commanded throttle"),
            (self.realized_throttle, "Realized throttle"),
        ):
            if not isfinite(value) or not 0.0 < value <= 1.0:
                raise DomainError(f"{label} must be finite and in (0, 1].")
        if self.active_engine_count < 1:
            raise DomainError("Active engine count must be at least one.")
        for value, label in (
            (self.ideal_thrust_n, "Ideal thrust"),
            (self.delivered_thrust_n, "Delivered thrust"),
            (self.chamber_specific_impulse_s, "Chamber specific impulse"),
            (self.system_specific_impulse_s, "System specific impulse"),
            (self.chamber_mass_flow_kg_s, "Chamber mass flow"),
            (self.total_tank_flow_kg_s, "Total tank flow"),
            (self.chamber_pressure_pa, "Chamber pressure"),
        ):
            _finite_positive(value, label)
        for value, label in (
            (self.oxidizer_mass_flow_kg_s, "Oxidizer mass flow"),
            (self.fuel_mass_flow_kg_s, "Fuel mass flow"),
            (self.ambient_pressure_pa, "Ambient pressure"),
        ):
            _finite_nonnegative(value, label)
        if self.delivered_thrust_n > self.ideal_thrust_n * (1.0 + 1e-12):
            raise DomainError("Delivered thrust cannot exceed ideal thrust.")
        if not self.streams:
            raise DomainError("At least one propellant stream is required.")

        tank_flow = sum(flow.mass_flow_kg_s for flow in self.streams if flow.tank_depleting)
        chamber_flow = sum(
            flow.mass_flow_kg_s
            for flow in self.streams
            if flow.destination is FlowDestination.CHAMBER
        )
        oxidizer_flow = sum(
            flow.mass_flow_kg_s
            for flow in self.streams
            if flow.tank_depleting and flow.role is PropellantRole.OXIDIZER
        )
        fuel_flow = sum(
            flow.mass_flow_kg_s
            for flow in self.streams
            if flow.tank_depleting and flow.role is PropellantRole.FUEL
        )
        for actual, declared, label in (
            (tank_flow, self.total_tank_flow_kg_s, "Total tank flow"),
            (chamber_flow, self.chamber_mass_flow_kg_s, "Chamber mass flow"),
            (oxidizer_flow, self.oxidizer_mass_flow_kg_s, "Oxidizer mass flow"),
            (fuel_flow, self.fuel_mass_flow_kg_s, "Fuel mass flow"),
        ):
            tolerance = max(1e-12, abs(declared) * 1e-10)
            if abs(actual - declared) > tolerance:
                raise DomainError(f"{label} does not close against explicit streams.")

        expected_isp = self.delivered_thrust_n / (
            STANDARD_GRAVITY_M_S2 * self.total_tank_flow_kg_s
        )
        if abs(expected_isp - self.system_specific_impulse_s) > max(
            1e-10, expected_isp * 1e-10
        ):
            raise DomainError(
                "System specific impulse must use delivered thrust and total tank flow."
            )
        if self.mixture_ratio is not None:
            _finite_positive(self.mixture_ratio, "Mixture ratio")
            if self.fuel_mass_flow_kg_s <= 0.0:
                raise DomainError("Mixture ratio requires a positive fuel flow.")
            expected_ratio = self.oxidizer_mass_flow_kg_s / self.fuel_mass_flow_kg_s
            if abs(expected_ratio - self.mixture_ratio) > max(1e-10, expected_ratio * 1e-10):
                raise DomainError("Mixture ratio does not close against oxidizer and fuel flow.")


@dataclass(frozen=True, slots=True)
class TankDefinition:
    """Initial and protected inventory for one propellant tank."""

    tank_id: str
    propellant_key: str
    role: PropellantRole
    loaded_mass_kg: float
    reserve_mass_kg: float = 0.0

    def __post_init__(self) -> None:
        if not self.tank_id.strip():
            raise DomainError("Tank identifier must not be empty.")
        if not self.propellant_key.strip():
            raise DomainError("Tank propellant key must not be empty.")
        _finite_nonnegative(self.loaded_mass_kg, "Loaded tank mass")
        _finite_nonnegative(self.reserve_mass_kg, "Tank reserve mass")
        if self.reserve_mass_kg > self.loaded_mass_kg:
            raise DomainError("Tank reserve cannot exceed loaded mass.")

    @property
    def expendable_mass_kg(self) -> float:
        """Mass available to this burn after protecting the reserve."""

        return self.loaded_mass_kg - self.reserve_mass_kg


@dataclass(frozen=True, slots=True)
class BurnTarget:
    """Exactly one primary target for a standalone burn."""

    kind: BurnTargetKind
    value: float

    def __post_init__(self) -> None:
        _finite_positive(self.value, "Burn target")


@dataclass(frozen=True, slots=True)
class BurnDefinition:
    """Normalized propulsion-only request for a constant burn."""

    name: str
    initial_mass_kg: float
    protected_dry_mass_kg: float
    tanks: tuple[TankDefinition, ...]
    target: BurnTarget
    maximum_duration_s: float
    cant_efficiency: float = 1.0

    def __post_init__(self) -> None:
        if not self.name.strip():
            raise DomainError("Burn name must not be empty.")
        _finite_positive(self.initial_mass_kg, "Initial vehicle mass")
        _finite_nonnegative(self.protected_dry_mass_kg, "Protected dry mass")
        _finite_positive(self.maximum_duration_s, "Maximum burn duration")
        if not isfinite(self.cant_efficiency) or not 0.0 < self.cant_efficiency <= 1.0:
            raise DomainError("Cant efficiency must be finite and in (0, 1].")
        if not self.tanks:
            raise DomainError("At least one tank is required.")
        tank_ids = [tank.tank_id for tank in self.tanks]
        if len(tank_ids) != len(set(tank_ids)):
            raise DomainError("Tank identifiers must be unique.")
        if self.protected_dry_mass_kg + sum(t.loaded_mass_kg for t in self.tanks) > (
            self.initial_mass_kg * (1.0 + 1e-12)
        ):
            raise DomainError("Protected dry mass plus tank loads exceeds initial mass.")


@dataclass(frozen=True, slots=True)
class PropulsionState:
    """Propulsion inventory and vehicle mass at one instant."""

    time_s: float
    vehicle_mass_kg: float
    tank_masses_kg: tuple[tuple[str, float], ...]


@dataclass(frozen=True, slots=True)
class BurnProfilePoint:
    """One exact endpoint from a constant burn."""

    time_s: float
    delivered_thrust_n: float
    axial_thrust_n: float
    system_specific_impulse_s: float
    total_mass_flow_kg_s: float
    vehicle_mass_kg: float
    tank_masses_kg: tuple[tuple[str, float], ...]
    phase: str


@dataclass(frozen=True, slots=True)
class ThrottleSegment:
    """One exact piecewise-linear throttle segment.

    Segments are ordered and contiguous by construction in
    :class:`ThrottleSchedule`. A discontinuous command is represented by an
    end throttle that differs from the next segment's start throttle.
    """

    duration_s: float
    start_throttle: float
    end_throttle: float
    phase: str = "commanded"

    def __post_init__(self) -> None:
        _finite_positive(self.duration_s, "Throttle-segment duration")
        for value, label in (
            (self.start_throttle, "Segment start throttle"),
            (self.end_throttle, "Segment end throttle"),
        ):
            if not isfinite(value) or not 0.0 <= value <= 1.0:
                raise DomainError(f"{label} must be finite and in [0, 1].")
        if not self.phase.strip():
            raise DomainError("Throttle-segment phase must not be empty.")

    @property
    def exposure_s(self) -> float:
        """Return the exact full-throttle-equivalent duration."""

        return 0.5 * (self.start_throttle + self.end_throttle) * self.duration_s


@dataclass(frozen=True, slots=True)
class ThrottleSchedule:
    """Ordered dependency-free L1 throttle/start/stop command history."""

    segments: tuple[ThrottleSegment, ...]
    name: str = "piecewise-linear throttle schedule"

    def __post_init__(self) -> None:
        if not self.name.strip():
            raise DomainError("Throttle-schedule name must not be empty.")
        if not isinstance(self.segments, tuple) or not self.segments:
            raise DomainError("Throttle schedule requires at least one segment.")
        if self.total_exposure_s <= 0.0:
            raise DomainError("Throttle schedule must contain positive commanded exposure.")

    @property
    def total_duration_s(self) -> float:
        """Return command time from the first segment to schedule end."""

        return sum(segment.duration_s for segment in self.segments)

    @property
    def total_exposure_s(self) -> float:
        """Return ``integral(throttle dt)`` over the complete schedule."""

        return sum(segment.exposure_s for segment in self.segments)


@dataclass(frozen=True, slots=True)
class ProfiledBurnPoint:
    """One exact L1 endpoint with throttle, force, flow, and inventory state."""

    time_s: float
    commanded_throttle: float
    realized_throttle: float
    delivered_thrust_n: float
    axial_thrust_n: float
    system_specific_impulse_s: float
    total_mass_flow_kg_s: float
    oxidizer_mass_flow_kg_s: float
    fuel_mass_flow_kg_s: float
    vehicle_mass_kg: float
    tank_masses_kg: tuple[tuple[str, float], ...]
    chamber_pressure_pa: float
    phase: str


@dataclass(frozen=True, slots=True)
class BurnSummary:
    """Conservation-aware summary of a standalone burn."""

    requested_target_kind: BurnTargetKind
    requested_target_value: float
    target_achieved: bool
    cutoff_reason: CutoffReason
    binding_constraint: str | None
    duration_s: float
    initial_mass_kg: float
    final_mass_kg: float
    consumed_propellant_kg: float
    consumed_by_tank_kg: tuple[tuple[str, float], ...]
    delivered_total_impulse_n_s: float
    axial_total_impulse_n_s: float
    average_delivered_thrust_n: float
    system_equivalent_specific_impulse_s: float
    ideal_delta_v_m_s: float
    rocket_equivalent_specific_impulse_s: float
    thrust_centroid_time_s: float
    mass_closure_error_kg: float
    impulse_closure_error_n_s: float
    delta_v_closure_error_m_s: float
    exact_constant_sidera_export_available: bool


@dataclass(frozen=True, slots=True)
class IntegrationEvidence:
    """Numerical-method evidence carried with every result."""

    method: str
    absolute_tolerance: float
    relative_tolerance: float
    segment_count: int
    accepted_steps: int
    rejected_steps: int


@dataclass(frozen=True, slots=True)
class BurnResult:
    """Top-level, Sidera-independent burn result contract."""

    schema: str
    definition: BurnDefinition
    operating_point: EngineOperatingPoint
    initial_state: PropulsionState
    final_state: PropulsionState
    profile: tuple[BurnProfilePoint, ...]
    summary: BurnSummary
    evidence: IntegrationEvidence
    warnings: tuple[str, ...]
    input_hash: str
    result_hash: str


@dataclass(frozen=True, slots=True)
class ProfiledBurnResult:
    """Top-level result for an exact piecewise-linear L1 throttle schedule."""

    schema: str
    definition: BurnDefinition
    operating_point: EngineOperatingPoint
    schedule: ThrottleSchedule
    initial_state: PropulsionState
    final_state: PropulsionState
    profile: tuple[ProfiledBurnPoint, ...]
    summary: BurnSummary
    evidence: IntegrationEvidence
    warnings: tuple[str, ...]
    input_hash: str
    result_hash: str

