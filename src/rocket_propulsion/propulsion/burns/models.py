"""Immutable public contracts for standalone propulsion-burn analysis.

The records in this module contain no position, velocity, orbit,
reference-frame or guidance quantities.  They describe propulsion hardware,
propellant inventory, a burn target and the resulting one-dimensional
propulsion bookkeeping only: total impulse and specific impulse [1]_ and the
ideal, external-force-free rocket equation [2]_.

References
----------
.. [1] NASA-ISP. NASA Glenn Research Center, "Specific Impulse",
       https://www1.grc.nasa.gov/beginners-guide-to-aeronautics/specific-impulse/.
.. [2] NASA-ROCKET. NASA Glenn Research Center, "Ideal Rocket Equation",
       https://www1.grc.nasa.gov/beginners-guide-to-aeronautics/ideal-rocket-equation/.
"""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass, fields
from enum import Enum
from math import isfinite

from rocket_propulsion.core.errors import DomainError
from rocket_propulsion.core.validation import (
    coerce_fields,
    name_value_pairs,
    real,
    records,
    strings,
    text,
)

from .numeric_outputs import _finite_rational, _fraction


class PropellantRole(str, Enum):
    """Functional role of a tank-depleting propellant stream.

    Oxidizer and fuel roles define the mixture ratio ``mdot_ox / mdot_fuel``
    of a bipropellant operating point [1]_.

    References
    ----------
    .. [1] SUTTON-BIBLARZ-2017. G. P. Sutton and O. Biblarz, "Rocket Propulsion
           Elements", 9th ed., Wiley, 2017.
    """

    FUEL = "fuel"
    OXIDIZER = "oxidizer"
    MONOPROPELLANT = "monopropellant"
    OTHER = "other"


class FlowDestination(str, Enum):
    """Destination of a propellant stream after it leaves its tank.

    Flows that bypass the chamber (gas-generator or coolant dumps, vents) still
    deplete tanks and therefore lower the system specific impulse relative to
    the chamber value [1]_.

    References
    ----------
    .. [1] SUTTON-BIBLARZ-2017. G. P. Sutton and O. Biblarz, "Rocket Propulsion
           Elements", 9th ed., Wiley, 2017.
    """

    CHAMBER = "chamber"
    GAS_GENERATOR_DUMP = "gas_generator_dump"
    COOLANT_DUMP = "coolant_dump"
    VENT = "vent"
    OTHER = "other"


class BurnTargetKind(str, Enum):
    """Supported primary termination targets for the constant-burn solver.

    ``IDEAL_DELTA_V`` is the external-force-free velocity change of the ideal
    rocket equation [1]_; ``TOTAL_IMPULSE`` is the time integral of
    delivered thrust [2]_.

    References
    ----------
    .. [1] NASA-ROCKET. NASA Glenn Research Center, "Ideal Rocket Equation",
           https://www1.grc.nasa.gov/beginners-guide-to-aeronautics/ideal-rocket-equation/.
    .. [2] NASA-ISP. NASA Glenn Research Center, "Specific Impulse",
           https://www1.grc.nasa.gov/beginners-guide-to-aeronautics/specific-impulse/.
    """

    DURATION = "duration"
    PROPELLANT_MASS = "propellant_mass"
    TOTAL_IMPULSE = "total_impulse"
    IDEAL_DELTA_V = "ideal_delta_v"


class CutoffReason(str, Enum):
    """Reason a burn simulation stopped."""

    TARGET_REACHED = "target_reached"
    TANK_RESERVE = "tank_reserve"
    DRY_MASS_FLOOR = "dry_mass_floor"
    MAXIMUM_DURATION = "maximum_duration"
    SCHEDULE_END = "schedule_end"


def _finite_nonnegative(value: float, label: str) -> None:
    """Raise ``DomainError`` unless ``value`` is a finite, non-negative real number."""
    if real(value, label) < 0.0:
        raise DomainError(f"{label} must be finite and non-negative.", field=label)


def _finite_positive(value: float, label: str) -> None:
    """Raise ``DomainError`` unless ``value`` is a finite real number greater than zero."""
    if real(value, label) <= 0.0:
        raise DomainError(f"{label} must be finite and greater than zero.", field=label)


def _boolean(value, label):
    """Return ``value`` if it is a ``bool``; text and numbers are rejected."""
    if not isinstance(value, bool):
        raise DomainError(f"{label} must be boolean.", field=label)
    return value


def _enum(value, kind, label):
    """Return the ``kind`` member for a recognized string value, or raise ``DomainError``."""
    if not isinstance(value, str):
        raise DomainError(f"{label} must be a recognized {kind.__name__}.", field=label)
    try:
        return kind(text(value, label))
    except ValueError as error:
        raise DomainError(f"{label} must be a recognized {kind.__name__}.", field=label) from error


def _utf8(value, label):
    """Require text encodable as UTF-8, as needed for canonical content hashes."""
    try:
        value.encode('utf-8')
    except UnicodeError as error:
        raise DomainError(f'{label} contains an unpaired Unicode surrogate.', field=label) from error


def _model_fields(record, label):
    """Coerce the scalar fields of ``record`` and require UTF-8-encodable text fields."""
    coerce_fields(record, label)
    for item in fields(record):
        value = getattr(record, item.name)
        if isinstance(value, str):
            _utf8(value, item.name)


def _model_strings(value, label):
    """Return a tuple snapshot of strings, each encodable as UTF-8."""
    result = strings(value, label)
    for item in result:
        _utf8(item, label)
    return result


def _inventory_pairs(value, label):
    """Return ``(name, mass)`` pairs with unique names and non-negative finite masses."""
    if isinstance(value, (str, bytes, Mapping)) or not hasattr(value, '__iter__'):
        raise DomainError(f"{label} must be a sequence of pairs.", field=label)
    items = tuple(value)
    if any(isinstance(item, (str, bytes, Mapping)) or not hasattr(item, '__len__')
           or len(item) != 2 for item in items):
        raise DomainError(f"{label} must contain exactly two-item pairs.", field=label)
    pairs = name_value_pairs(items, label)
    if len({name for name, _ in pairs}) != len(pairs):
        raise DomainError(f"{label} identifiers must be unique.", field=label)
    for name, amount in pairs:
        _utf8(name, label)
        _finite_nonnegative(amount, label)
    return pairs


def _point_contract(point):
    """Validate a profile sample, including zero-thrust (idle) samples.

    Float fields are finite and non-negative, the vehicle mass is positive,
    axial thrust does not exceed delivered thrust, and tank inventories are
    snapshotted.
    """
    _model_fields(point, type(point).__name__)
    for item in fields(point):
        if item.type == 'float':
            _finite_nonnegative(getattr(point, item.name), item.name)
    _finite_positive(point.vehicle_mass_kg, 'vehicle_mass_kg')
    if point.axial_thrust_n > point.delivered_thrust_n:
        raise DomainError('Axial thrust cannot exceed delivered thrust.', field='axial_thrust_n')
    object.__setattr__(point, 'tank_masses_kg', _inventory_pairs(point.tank_masses_kg, 'tank_masses_kg'))


@dataclass(frozen=True, slots=True)
class PropellantFlow:
    """One explicit mass-flow path leaving a named tank.

    Parameters
    ----------
    stream_id : str
        Non-empty stream identifier.
    tank_id : str
        Source tank; may be blank only for a stream that does not deplete a
        tank.
    propellant_key : str
        Non-empty propellant identifier.
    role : PropellantRole
        Functional role of the stream.
    destination : FlowDestination
        Where the flow goes after leaving the tank.
    mass_flow_kg_s : float
        Finite, non-negative stream mass flow, in kg/s.
    tank_depleting : bool, optional
        Whether the stream draws from ``tank_id``.
    contributes_to_delivered_thrust : bool, optional
        Whether the stream is part of the delivered-thrust accounting.

    Raises
    ------
    DomainError
        If a field has the wrong type, an identifier is blank or not
        UTF-8-encodable, or the flow is negative or non-finite.
    """

    stream_id: str
    tank_id: str
    propellant_key: str
    role: PropellantRole
    destination: FlowDestination
    mass_flow_kg_s: float
    tank_depleting: bool = True
    contributes_to_delivered_thrust: bool = True

    def __post_init__(self) -> None:
        object.__setattr__(self, 'stream_id', text(self.stream_id, 'stream_id'))
        object.__setattr__(self, 'propellant_key', text(self.propellant_key, 'propellant_key'))
        _boolean(self.tank_depleting, 'tank_depleting')
        _boolean(self.contributes_to_delivered_thrust, 'contributes_to_delivered_thrust')
        object.__setattr__(self, 'tank_id', strings((self.tank_id,), 'tank_id',
                           allow_blank=not self.tank_depleting)[0])
        for name in ('stream_id', 'propellant_key', 'tank_id'):
            _utf8(getattr(self, name), name)
        object.__setattr__(self, 'role', _enum(self.role, PropellantRole, 'role'))
        object.__setattr__(self, 'destination', _enum(self.destination, FlowDestination, 'destination'))
        object.__setattr__(self, 'mass_flow_kg_s', real(self.mass_flow_kg_s, 'mass_flow_kg_s'))
        _finite_nonnegative(self.mass_flow_kg_s, "Stream mass flow")


@dataclass(frozen=True, slots=True)
class EngineOperatingPoint:
    """Delivered steady performance with explicit tank-depleting streams.

    Parameters
    ----------
    name : str
        Operating-point name.
    commanded_throttle, realized_throttle : float
        Throttle fractions in ``(0, 1]``.
    active_engine_count : int
        Number of firing engines, at least one.
    ideal_thrust_n, delivered_thrust_n : float
        Positive ideal and delivered thrust, in N; delivered thrust does not
        exceed ideal thrust.
    chamber_specific_impulse_s, system_specific_impulse_s : float
        Positive chamber and system specific impulse, in s.
    chamber_mass_flow_kg_s, total_tank_flow_kg_s : float
        Positive flow into the chamber and total tank-depleting flow, in kg/s.
    oxidizer_mass_flow_kg_s, fuel_mass_flow_kg_s : float
        Non-negative tank-depleting oxidizer and fuel flows, in kg/s.
    mixture_ratio : float or None
        Oxidizer-to-fuel mass-flow ratio, or ``None``.
    chamber_pressure_pa : float
        Positive chamber pressure, in Pa.
    ambient_pressure_pa : float
        Non-negative ambient pressure, in Pa.
    streams : tuple of PropellantFlow
        At least one explicit stream.
    model_id, source : str, optional
        Model and provenance identifiers.
    warnings : tuple of str, optional
        Interpretation warnings.

    Raises
    ------
    DomainError
        If a scalar is out of range, there is no stream, a declared flow total
        does not close against the streams, the system specific impulse is
        inconsistent, or the mixture ratio does not close.

    Notes
    -----
    Thrust and flows are totals for ``active_engine_count``; downstream
    solvers never multiply them by the engine count a second time.  Declared
    flow totals must match the stream sums within a relative tolerance of
    ``1e-10`` (absolute floor ``1e-12``).  The system specific impulse is
    ``Isp_sys = F_delivered / (g0 mdot_tank)`` [1]_ with ``g0`` the registered
    ``STANDARD_GRAVITY_M_S2``, so flows that bypass the chamber lower it
    below the chamber value [2]_.  The mixture ratio, when given, is
    ``mdot_ox / mdot_fuel`` [2]_.

    References
    ----------
    .. [1] NASA-ISP. NASA Glenn Research Center, "Specific Impulse",
           https://www1.grc.nasa.gov/beginners-guide-to-aeronautics/specific-impulse/.
    .. [2] SUTTON-BIBLARZ-2017. G. P. Sutton and O. Biblarz, "Rocket Propulsion
           Elements", 9th ed., Wiley, 2017.
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

        _model_fields(self, 'Operating point')
        object.__setattr__(self, 'streams', records(self.streams, PropellantFlow, 'streams'))
        object.__setattr__(self, 'warnings', _model_strings(self.warnings, 'warnings'))
        for value, label in (
            (self.commanded_throttle, "Commanded throttle"),
            (self.realized_throttle, "Realized throttle"),
        ):
            if not isfinite(value) or not 0.0 < value <= 1.0:
                raise DomainError(f"{label} must be finite and in (0, 1].", field=label)
        if self.active_engine_count < 1:
            raise DomainError("Active engine count must be at least one.", field='active_engine_count')
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
        if _fraction(self.delivered_thrust_n) > _fraction(self.ideal_thrust_n) * _fraction(1.0 + 1e-12):
            raise DomainError("Delivered thrust cannot exceed ideal thrust.", field='delivered_thrust_n')
        if not self.streams:
            raise DomainError("At least one propellant stream is required.", field='streams')

        tank_flow = sum((_fraction(flow.mass_flow_kg_s) for flow in self.streams if flow.tank_depleting), _fraction(0.))
        chamber_flow = sum(
            (_fraction(flow.mass_flow_kg_s)
            for flow in self.streams
            if flow.destination is FlowDestination.CHAMBER), _fraction(0.)
        )
        oxidizer_flow = sum(
            (_fraction(flow.mass_flow_kg_s)
            for flow in self.streams
            if flow.tank_depleting and flow.role is PropellantRole.OXIDIZER), _fraction(0.)
        )
        fuel_flow = sum(
            (_fraction(flow.mass_flow_kg_s)
            for flow in self.streams
            if flow.tank_depleting and flow.role is PropellantRole.FUEL), _fraction(0.)
        )
        for actual, declared, label in (
            (tank_flow, self.total_tank_flow_kg_s, "Total tank flow"),
            (chamber_flow, self.chamber_mass_flow_kg_s, "Chamber mass flow"),
            (oxidizer_flow, self.oxidizer_mass_flow_kg_s, "Oxidizer mass flow"),
            (fuel_flow, self.fuel_mass_flow_kg_s, "Fuel mass flow"),
        ):
            tolerance = max(_fraction(1e-12), abs(_fraction(declared)) * _fraction(1e-10))
            if abs(actual - _fraction(declared)) > tolerance:
                raise DomainError(f"{label} does not close against explicit streams.", field=label)

        expected_isp = _fraction(self.delivered_thrust_n) / (_fraction(STANDARD_GRAVITY_M_S2) * _fraction(self.total_tank_flow_kg_s))
        if abs(expected_isp - _fraction(self.system_specific_impulse_s)) > max(_fraction(1e-10), expected_isp * _fraction(1e-10)):
            raise DomainError(
                "System specific impulse must use delivered thrust and total tank flow.", field='system_specific_impulse_s'
            )
        if self.mixture_ratio is not None:
            _finite_positive(self.mixture_ratio, "Mixture ratio")
            if self.fuel_mass_flow_kg_s <= 0.0:
                raise DomainError("Mixture ratio requires a positive fuel flow.", field='mixture_ratio')
            expected_ratio = _fraction(self.oxidizer_mass_flow_kg_s) / _fraction(self.fuel_mass_flow_kg_s)
            if abs(expected_ratio - _fraction(self.mixture_ratio)) > max(_fraction(1e-10), expected_ratio * _fraction(1e-10)):
                raise DomainError("Mixture ratio does not close against oxidizer and fuel flow.", field='mixture_ratio')


@dataclass(frozen=True, slots=True)
class TankDefinition:
    """Initial and protected inventory for one propellant tank.

    Parameters
    ----------
    tank_id : str
        Tank identifier.
    propellant_key : str
        Propellant identifier.
    role : PropellantRole
        Functional role of the stored propellant.
    loaded_mass_kg : float
        Finite, non-negative loaded propellant mass, in kg.
    reserve_mass_kg : float, optional
        Protected reserve, in kg; between zero and ``loaded_mass_kg``.

    Raises
    ------
    DomainError
        If a mass is negative or non-finite, or the reserve exceeds the load.
    """

    tank_id: str
    propellant_key: str
    role: PropellantRole
    loaded_mass_kg: float
    reserve_mass_kg: float = 0.0

    def __post_init__(self) -> None:
        _model_fields(self, 'Tank')
        object.__setattr__(self, 'role', _enum(self.role, PropellantRole, 'role'))
        _finite_nonnegative(self.loaded_mass_kg, "Loaded tank mass")
        _finite_nonnegative(self.reserve_mass_kg, "Tank reserve mass")
        if self.reserve_mass_kg > self.loaded_mass_kg:
            raise DomainError("Tank reserve cannot exceed loaded mass.", field='reserve_mass_kg')

    @property
    def expendable_mass_kg(self) -> float:
        """Mass available to this burn after protecting the reserve, in kg."""

        return self.loaded_mass_kg - self.reserve_mass_kg


@dataclass(frozen=True, slots=True)
class BurnTarget:
    """Exactly one primary target for a standalone burn.

    Parameters
    ----------
    kind : BurnTargetKind
        Target quantity.
    value : float
        Finite, positive target value in the unit of ``kind`` (s, kg, N s or
        m/s).

    Raises
    ------
    DomainError
        If ``kind`` is not recognized or ``value`` is not finite and positive.
    """

    kind: BurnTargetKind
    value: float

    def __post_init__(self) -> None:
        _model_fields(self, 'Burn target')
        object.__setattr__(self, 'kind', _enum(self.kind, BurnTargetKind, 'kind'))
        _finite_positive(self.value, "Burn target")


@dataclass(frozen=True, slots=True)
class BurnDefinition:
    """Normalized propulsion-only request for a constant burn.

    Parameters
    ----------
    name : str
        Burn name.
    initial_mass_kg : float
        Positive initial vehicle mass, in kg.
    protected_dry_mass_kg : float
        Non-negative mass floor the burn must not cross, in kg.
    tanks : tuple of TankDefinition
        At least one tank, with unique identifiers.
    target : BurnTarget
        Primary termination target.
    maximum_duration_s : float
        Positive upper bound on burn duration, in s.
    cant_efficiency : float, optional
        Ratio of axial to delivered thrust, in ``(0, 1]``.

    Raises
    ------
    DomainError
        If a field is out of range, tanks are missing or duplicated, or the
        protected dry mass plus tank loads exceeds the initial mass by more
        than a relative ``1e-12``.

    Notes
    -----
    ``cant_efficiency`` is the cosine loss of canted thrusters: only the
    thrust component along the vehicle axis contributes to axial impulse
    [1]_.

    References
    ----------
    .. [1] NASA SP-8059, "Spacecraft Attitude Control During Thrusting
           Maneuvers", https://ntrs.nasa.gov/citations/19710016722.
    """

    name: str
    initial_mass_kg: float
    protected_dry_mass_kg: float
    tanks: tuple[TankDefinition, ...]
    target: BurnTarget
    maximum_duration_s: float
    cant_efficiency: float = 1.0

    def __post_init__(self) -> None:
        _model_fields(self, 'Burn definition')
        object.__setattr__(self, 'tanks', records(self.tanks, TankDefinition, 'tanks'))
        if not isinstance(self.target, BurnTarget):
            raise DomainError('target must be a BurnTarget.', field='target')
        _finite_positive(self.initial_mass_kg, "Initial vehicle mass")
        _finite_nonnegative(self.protected_dry_mass_kg, "Protected dry mass")
        _finite_positive(self.maximum_duration_s, "Maximum burn duration")
        if not isfinite(self.cant_efficiency) or not 0.0 < self.cant_efficiency <= 1.0:
            raise DomainError("Cant efficiency must be finite and in (0, 1].", field='cant_efficiency')
        if not self.tanks:
            raise DomainError("At least one tank is required.", field='tanks')
        tank_ids = [tank.tank_id for tank in self.tanks]
        if len(tank_ids) != len(set(tank_ids)):
            raise DomainError("Tank identifiers must be unique.", field='tanks')
        if _fraction(self.protected_dry_mass_kg) + sum((_fraction(t.loaded_mass_kg) for t in self.tanks), _fraction(0.)) > (
            _fraction(self.initial_mass_kg) * _fraction(1.0 + 1e-12)
        ):
            raise DomainError("Protected dry mass plus tank loads exceeds initial mass.", field='initial_mass_kg')


@dataclass(frozen=True, slots=True)
class PropulsionState:
    """Propulsion inventory and vehicle mass at one instant.

    Parameters
    ----------
    time_s : float
        Non-negative burn time, in s.
    vehicle_mass_kg : float
        Positive vehicle mass, in kg.
    tank_masses_kg : tuple of (str, float)
        Unique tank identifiers with non-negative masses, in kg.

    Raises
    ------
    DomainError
        If a field is out of range or tank identifiers repeat.
    """

    time_s: float
    vehicle_mass_kg: float
    tank_masses_kg: tuple[tuple[str, float], ...]

    def __post_init__(self) -> None:
        _model_fields(self, 'Propulsion state')
        _finite_nonnegative(self.time_s, 'time_s')
        _finite_positive(self.vehicle_mass_kg, 'vehicle_mass_kg')
        object.__setattr__(self, 'tank_masses_kg', _inventory_pairs(self.tank_masses_kg, 'tank_masses_kg'))


@dataclass(frozen=True, slots=True)
class BurnProfilePoint:
    """One binary64 sample of a constant burn.

    Parameters
    ----------
    time_s : float
        Burn time, in s.
    delivered_thrust_n, axial_thrust_n : float
        Delivered and axial thrust, in N; axial does not exceed delivered.
    system_specific_impulse_s : float
        System specific impulse, in s.
    total_mass_flow_kg_s : float
        Total tank-depleting flow, in kg/s.
    vehicle_mass_kg : float
        Positive vehicle mass, in kg.
    tank_masses_kg : tuple of (str, float)
        Tank inventories, in kg.
    phase : str
        Phase label.

    Raises
    ------
    DomainError
        If a float field is negative or non-finite, the vehicle mass is not
        positive, or axial thrust exceeds delivered thrust.
    """

    time_s: float
    delivered_thrust_n: float
    axial_thrust_n: float
    system_specific_impulse_s: float
    total_mass_flow_kg_s: float
    vehicle_mass_kg: float
    tank_masses_kg: tuple[tuple[str, float], ...]
    phase: str

    def __post_init__(self) -> None:
        _point_contract(self)


@dataclass(frozen=True, slots=True)
class ThrottleSegment:
    """One piecewise-linear throttle segment defined by supplied binary values.

    Parameters
    ----------
    duration_s : float
        Finite, positive segment duration, in s.
    start_throttle, end_throttle : float
        Throttle fractions in ``[0, 1]`` at the segment ends.
    phase : str, optional
        Non-empty phase label.

    Raises
    ------
    DomainError
        If the duration is not positive or a throttle is outside ``[0, 1]``.

    Notes
    -----
    Segments are ordered and contiguous by construction in
    :class:`ThrottleSchedule`.  A discontinuous command is an end throttle
    that differs from the next segment's start throttle.
    """

    duration_s: float
    start_throttle: float
    end_throttle: float
    phase: str = "commanded"

    def __post_init__(self) -> None:
        _model_fields(self, 'Throttle segment')
        _finite_positive(self.duration_s, "Throttle-segment duration")
        for value, label in (
            (self.start_throttle, "Segment start throttle"),
            (self.end_throttle, "Segment end throttle"),
        ):
            if not isfinite(value) or not 0.0 <= value <= 1.0:
                raise DomainError(f"{label} must be finite and in [0, 1].", field=label)
        text(self.phase, "Throttle-segment phase")

    @property
    def exposure_s(self) -> float:
        """Full-throttle exposure ``(start + end) duration / 2``, in s.

        The integral of the linear throttle is evaluated exactly from the
        binary inputs and rounded once on export.
        """

        return _finite_rational(
            (_fraction(self.start_throttle) + _fraction(self.end_throttle))
            * _fraction(self.duration_s)
            / 2,
            "Throttle segment exposure",
        )


@dataclass(frozen=True, slots=True)
class ThrottleSchedule:
    """Ordered, self-contained piecewise-linear throttle command history.

    Parameters
    ----------
    segments : tuple of ThrottleSegment
        At least one segment; together they must command positive exposure.
    name : str, optional
        Schedule name.

    Raises
    ------
    DomainError
        If there is no segment or every throttle value is zero.
    """

    segments: tuple[ThrottleSegment, ...]
    name: str = "piecewise-linear throttle schedule"

    def __post_init__(self) -> None:
        _model_fields(self, 'Throttle schedule')
        object.__setattr__(self, 'segments', records(self.segments, ThrottleSegment, 'segments'))
        if not self.segments:
            raise DomainError("Throttle schedule requires at least one segment.", field='segments')
        if not any(
            segment.start_throttle > 0 or segment.end_throttle > 0 for segment in self.segments
        ):
            raise DomainError("Throttle schedule must contain positive commanded exposure.", field='segments')

    @property
    def total_duration_s(self) -> float:
        """Command time from the start of the first segment to the schedule end, in s."""

        return _finite_rational(
            sum((_fraction(segment.duration_s) for segment in self.segments), _fraction(0.0)),
            "Schedule duration",
        )

    @property
    def total_exposure_s(self) -> float:
        """Full-throttle exposure ``integral(throttle dt)`` over the schedule, in s."""

        return _finite_rational(
            sum(
                (
                    (_fraction(segment.start_throttle) + _fraction(segment.end_throttle))
                    * _fraction(segment.duration_s)
                    / 2
                    for segment in self.segments
                ),
                _fraction(0.0),
            ),
            "Schedule exposure",
        )


@dataclass(frozen=True, slots=True)
class ProfiledBurnPoint:
    """One sample of a profiled burn with rounded throttle, force, flow and inventory.

    Fields follow :class:`BurnProfilePoint`, plus commanded and realized
    throttle in ``[0, 1]``, oxidizer and fuel flows in kg/s and chamber
    pressure in Pa.

    Raises
    ------
    DomainError
        If a float field is negative or non-finite, the vehicle mass is not
        positive, axial thrust exceeds delivered thrust, or a throttle exceeds
        one.
    """

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

    def __post_init__(self) -> None:
        _point_contract(self)
        for name in ('commanded_throttle', 'realized_throttle'):
            if getattr(self, name) > 1:
                raise DomainError(f'{name} must be in [0, 1].', field=name)


@dataclass(frozen=True, slots=True)
class BurnSummary:
    """Conservation-checked summary of a standalone burn.

    Attributes include the requested target and whether it was achieved, the
    cutoff reason and binding constraint, duration, initial and final mass,
    consumed propellant (total and per tank), delivered and axial total
    impulse, average delivered thrust, equivalent specific impulses, ideal
    delta-v, thrust-centroid time, the mass, impulse and delta-v closure
    errors, and whether an exact constant-thrust export to Sidera is
    available.

    Raises
    ------
    DomainError
        If an enum or flag has the wrong type, a float field is negative or
        non-finite, the target value or a mass is not positive, the final mass
        exceeds the initial mass, or axial impulse exceeds delivered impulse.

    Notes
    -----
    Total impulse is ``integral F dt`` and the system-equivalent specific
    impulse is ``I / (g0 m_consumed)`` [1]_.  The ideal delta-v and the
    rocket-equivalent specific impulse follow the external-force-free rocket
    equation ``dv = g0 Isp ln(m0 / mf)`` [2]_, with ``g0`` the registered
    ``STANDARD_GRAVITY_M_S2``.  The closure errors are numerical
    consistency checks, not physical uncertainty.

    References
    ----------
    .. [1] NASA-ISP. NASA Glenn Research Center, "Specific Impulse",
           https://www1.grc.nasa.gov/beginners-guide-to-aeronautics/specific-impulse/.
    .. [2] NASA-ROCKET. NASA Glenn Research Center, "Ideal Rocket Equation",
           https://www1.grc.nasa.gov/beginners-guide-to-aeronautics/ideal-rocket-equation/.
    """

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

    def __post_init__(self) -> None:
        _model_fields(self, 'Burn summary')
        object.__setattr__(self, 'requested_target_kind', _enum(self.requested_target_kind, BurnTargetKind, 'requested_target_kind'))
        object.__setattr__(self, 'cutoff_reason', _enum(self.cutoff_reason, CutoffReason, 'cutoff_reason'))
        _boolean(self.target_achieved, 'target_achieved')
        _boolean(self.exact_constant_sidera_export_available, 'exact_constant_sidera_export_available')
        if self.binding_constraint is not None:
            object.__setattr__(self, 'binding_constraint',
                               text(self.binding_constraint, 'binding_constraint'))
        for item in fields(self):
            if item.type == 'float':
                _finite_nonnegative(getattr(self, item.name), item.name)
        for name in ('requested_target_value', 'initial_mass_kg', 'final_mass_kg'):
            _finite_positive(getattr(self, name), name)
        if self.final_mass_kg > self.initial_mass_kg:
            raise DomainError('Burn final mass exceeds initial mass.', field='final_mass_kg')
        if self.axial_total_impulse_n_s > self.delivered_total_impulse_n_s:
            raise DomainError('Axial impulse exceeds delivered impulse.', field='axial_total_impulse_n_s')
        object.__setattr__(self, 'consumed_by_tank_kg', _inventory_pairs(self.consumed_by_tank_kg, 'consumed_by_tank_kg'))


@dataclass(frozen=True, slots=True)
class IntegrationEvidence:
    """Numerical-method evidence carried with every result.

    Parameters
    ----------
    method : str
        Integration method name.
    absolute_tolerance, relative_tolerance : float
        Non-negative solver tolerances.
    segment_count, accepted_steps, rejected_steps : int
        Non-negative step accounting.

    Raises
    ------
    DomainError
        If a tolerance or count is negative.

    Notes
    -----
    The record documents numerical verification evidence as required for
    model reporting by [1]_; it is not validation evidence.

    References
    ----------
    .. [1] NASA-STD-7009B, "Standard for Models and Simulations",
           https://standards.nasa.gov/standard/NASA/NASA-STD-7009.
    """

    method: str
    absolute_tolerance: float
    relative_tolerance: float
    segment_count: int
    accepted_steps: int
    rejected_steps: int

    def __post_init__(self) -> None:
        _model_fields(self, 'Integration evidence')
        for name in ('absolute_tolerance', 'relative_tolerance'):
            _finite_nonnegative(getattr(self, name), name)
        for name in ('segment_count', 'accepted_steps', 'rejected_steps'):
            if getattr(self, name) < 0:
                raise DomainError(f'{name} must be nonnegative.', field=name)


@dataclass(frozen=True, slots=True)
class BurnResult:
    """Top-level result of a constant burn, independent of any trajectory tool.

    Parameters
    ----------
    schema : str
        Must equal ``"rocket_propulsion_burn_v1"``.
    definition : BurnDefinition
        The normalized request.
    operating_point : EngineOperatingPoint
        The engine operating point used.
    initial_state, final_state : PropulsionState
        Inventory at burn start and end.
    profile : tuple of BurnProfilePoint
        Sampled burn history.
    summary : BurnSummary
        Conservation-checked summary.
    evidence : IntegrationEvidence
        Numerical-method evidence.
    warnings : tuple of str
        Interpretation warnings.
    input_hash, result_hash : str
        Lower-case SHA-256 digests; ``result_hash`` may be empty while the
        record is provisional.

    Raises
    ------
    DomainError
        If the schema or a component type is wrong or a digest is malformed.
    """

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

    def __post_init__(self) -> None:
        _result_contract(self, 'rocket_propulsion_burn_v1', BurnProfilePoint)


@dataclass(frozen=True, slots=True)
class ProfiledBurnResult:
    """Top-level result of a burn under a piecewise-linear throttle schedule.

    Fields follow :class:`BurnResult`, with ``schema`` equal to
    ``"rocket_propulsion_burn_v2"``, an added ``schedule`` and
    :class:`ProfiledBurnPoint` samples.

    Raises
    ------
    DomainError
        If the schema or a component type is wrong or a digest is malformed.
    """

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

    def __post_init__(self) -> None:
        _result_contract(self, 'rocket_propulsion_burn_v2', ProfiledBurnPoint)
        if not isinstance(self.schedule, ThrottleSchedule):
            raise DomainError('schedule must be a ThrottleSchedule.', field='schedule')


def _result_contract(result, schema, point_kind):
    """Check a result's schema, component types and digests.

    Reproduction of the hashes from content is verified when a document is
    read, not here.
    """
    if result.schema != schema:
        raise DomainError('Result schema does not match its record type.', field='schema')
    for name, kind in (('definition', BurnDefinition), ('operating_point', EngineOperatingPoint),
                       ('initial_state', PropulsionState), ('final_state', PropulsionState),
                       ('summary', BurnSummary), ('evidence', IntegrationEvidence)):
        if not isinstance(getattr(result, name), kind):
            raise DomainError(f'{name} must be a {kind.__name__}.', field=name)
    object.__setattr__(result, 'profile', records(result.profile, point_kind, 'profile'))
    object.__setattr__(result, 'warnings', _model_strings(result.warnings, 'warnings'))
    for name in ('input_hash', 'result_hash'):
        digest = getattr(result, name)
        # Solvers construct a provisional record before computing its own hash.
        if name == 'result_hash' and digest == '':
            continue
        if not isinstance(digest, str) or len(digest) != 64 or any(c not in '0123456789abcdef' for c in digest):
            raise DomainError(f'{name} must be a lowercase SHA-256 digest.', field=name)
