"""Regulated pressurant and propellant-feed dynamics for orbit propagation.

The model is a compact zero-dimensional network [1]_, [3]_: a finite pressurant
bottle feeds a growing tank ullage through a dynamic regulator [2]_; liquid
withdrawal creates ullage volume, and a quadratic feed-line loss sets the engine
inlet pressure.  Regulator gas flow follows the compressible-restriction
relation including choking [4]_.  The model supplies additional states for a
numerical propagator; it does not replace validated nodal or water-hammer
tools.

References
----------
.. [1] NASA SP-8112. NASA, "Pressurization Systems for Liquid Rockets", 1975.
       https://ntrs.nasa.gov/citations/19760015212
.. [2] NASA SP-8080. NASA, "Liquid Rocket Pressure Regulators and Valves".
       https://ntrs.nasa.gov/citations/19740002611
.. [3] NASA, "Nodal Modeling of Liquid Propellant Feed and Pressurization
       System", 2024. https://ntrs.nasa.gov/citations/20240003493
.. [4] NASA Glenn Research Center, "Mass Flow Choking".
       https://www.grc.nasa.gov/www/k-12/BGP/mflchk.html
"""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass
from enum import Enum
from math import exp, expm1, isfinite, log, log1p, sqrt
from typing import TYPE_CHECKING

from rocket_propulsion.core.errors import DomainError
from rocket_propulsion.core.validation import coerce_fields, name_value_pairs, real

from .serialization import canonical_sha256

if TYPE_CHECKING:
    from .propagator import PropagatorPropulsionBridge, PropulsionDynamicsEvaluation

FEED_SYSTEM_SCHEMA = "rocket_propulsion_regulated_feed_system_v1"
FEED_SYSTEM_REFERENCE_IDS = (
    "NASA-SP-8112",
    "NASA-SP-8080",
    "NASA-20240003493",
    "NASA-GRC-MASS-FLOW-CHOKING",
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
class IdealPressurantGas:
    """Calorically perfect pressurant properties shared by both gas nodes.

    Parameters
    ----------
    specific_gas_constant_j_kg_k : float
        Positive specific gas constant ``R``, J/(kg K).
    gamma : float
        Ratio of specific heats, finite and greater than one.
    name : str, default "pressurant"
        Non-blank gas label.

    Raises
    ------
    DomainError
        If ``R`` or the derived heat capacities are not finite and positive,
        ``gamma <= 1``, or the name is blank.

    Notes
    -----
    The heat capacities follow ``cp = cv + R`` and ``gamma = cp / cv`` [2]_.
    The pressurant is a single ideal gas in both the bottle and the tank
    ullage, as in pressurization-system sizing [1]_ and nodal feed modelling
    [3]_.

    References
    ----------
    .. [1] NASA SP-8112. NASA, "Pressurization Systems for Liquid Rockets", 1975.
           https://ntrs.nasa.gov/citations/19760015212
    .. [2] NASA Glenn Research Center, "Specific Heats".
           https://www.grc.nasa.gov/www/k-12/airplane/specheat.html
    .. [3] NASA, "Nodal Modeling of Liquid Propellant Feed and Pressurization
           System", 2024. https://ntrs.nasa.gov/citations/20240003493
    """

    specific_gas_constant_j_kg_k: float
    gamma: float
    name: str = "pressurant"

    def __post_init__(self) -> None:
        coerce_fields(self, "Pressurant gas")
        _positive(self.specific_gas_constant_j_kg_k, "Pressurant gas constant")
        if not isfinite(self.gamma) or self.gamma <= 1.0:
            raise DomainError("Pressurant gamma must be finite and greater than one.")
        _positive(self.specific_heat_cv_j_kg_k, "Pressurant constant-volume heat capacity")
        _positive(self.specific_heat_cp_j_kg_k, "Pressurant constant-pressure heat capacity")
        if not self.name.strip():
            raise DomainError("Pressurant name must not be empty.")

    @property
    def specific_heat_cv_j_kg_k(self) -> float:
        """Constant-volume heat capacity ``cv = R / (gamma - 1)`` [1]_, J/(kg K).

        References
        ----------
        .. [1] NASA Glenn Research Center, "Specific Heats".
               https://www.grc.nasa.gov/www/k-12/airplane/specheat.html
        """

        return self.specific_gas_constant_j_kg_k / (self.gamma - 1.0)

    @property
    def specific_heat_cp_j_kg_k(self) -> float:
        """Constant-pressure heat capacity ``cp = gamma R / (gamma - 1)`` [1]_, J/(kg K).

        References
        ----------
        .. [1] NASA Glenn Research Center, "Specific Heats".
               https://www.grc.nasa.gov/www/k-12/airplane/specheat.html
        """

        return self.gamma * self.specific_heat_cv_j_kg_k


@dataclass(frozen=True, slots=True)
class RegulatedFeedConfiguration:
    """Immutable geometry, control, thermal and cutoff configuration.

    Parameters
    ----------
    gas : IdealPressurantGas
        Pressurant properties.
    tank_internal_volume_m3 : float
        Positive rigid internal tank volume, m^3.
    propellant_density_kg_m3 : float
        Positive constant liquid density, kg/m^3.
    bottle_volume_m3 : float
        Positive rigid pressurant-bottle volume, m^3.
    regulator_set_pressure_pa : float
        Positive regulated tank pressure, Pa.
    regulator_full_open_error_pa : float
        Positive tank-pressure error at which the command saturates, Pa.
    regulator_maximum_area_m2 : float
        Positive fully open regulator flow area, m^2.
    regulator_discharge_coefficient : float
        Discharge coefficient in ``(0, 1]``.
    valve_time_constant_s : float
        Positive first-order regulator time constant, s.
    feed_line_resistance_pa_s2_kg2 : float
        Non-negative quadratic liquid-line resistance, Pa s^2/kg^2.
    tank_heat_transfer_w_k : float
        Non-negative tank wall-to-gas conductance, W/K.
    tank_wall_temperature_k : float
        Positive tank wall temperature, K.
    bottle_heat_transfer_w_k : float
        Non-negative bottle wall-to-gas conductance, W/K.
    bottle_wall_temperature_k : float
        Positive bottle wall temperature, K.
    minimum_propellant_mass_kg : float, default 0.0
        Non-negative propellant-reserve cutoff, kg.
    minimum_injector_pressure_pa : float, default 0.0
        Non-negative injector-pressure cutoff, Pa.
    minimum_bottle_pressure_margin_pa : float, default 0.0
        Non-negative bottle-minus-tank pressure margin cutoff, Pa.
    maximum_tank_pressure_pa : float or None, default None
        Optional tank overpressure cutoff, Pa; must exceed the set pressure.
    model_id : str, default "regulated-single-node-feed"
        Non-blank model identifier.
    schema : str, default FEED_SYSTEM_SCHEMA
        Configuration schema; only the current schema is accepted.

    Raises
    ------
    DomainError
        If the schema is unsupported, the identifier is blank, a value is
        outside its domain, or the overpressure limit does not exceed the
        set pressure.

    Notes
    -----
    The regulator command is proportional to the tank-pressure error and
    clipped to ``[0, 1]``; the opening follows it with a first-order lag.
    Flow through the effective regulator area ``Cd * A_max * opening``
    follows the isentropic compressible-restriction relation, including
    choking [3]_.  The liquid line loss is ``delta_p = R_line * mdot**2``.
    [1]_ and [2]_ give the pressurization-system and regulator design context.

    References
    ----------
    .. [1] NASA SP-8112. NASA, "Pressurization Systems for Liquid Rockets", 1975.
           https://ntrs.nasa.gov/citations/19760015212
    .. [2] NASA SP-8080. NASA, "Liquid Rocket Pressure Regulators and Valves".
           https://ntrs.nasa.gov/citations/19740002611
    .. [3] NASA Glenn Research Center, "Mass Flow Choking".
           https://www.grc.nasa.gov/www/k-12/BGP/mflchk.html
    """

    gas: IdealPressurantGas
    tank_internal_volume_m3: float
    propellant_density_kg_m3: float
    bottle_volume_m3: float
    regulator_set_pressure_pa: float
    regulator_full_open_error_pa: float
    regulator_maximum_area_m2: float
    regulator_discharge_coefficient: float
    valve_time_constant_s: float
    feed_line_resistance_pa_s2_kg2: float
    tank_heat_transfer_w_k: float
    tank_wall_temperature_k: float
    bottle_heat_transfer_w_k: float
    bottle_wall_temperature_k: float
    minimum_propellant_mass_kg: float = 0.0
    minimum_injector_pressure_pa: float = 0.0
    minimum_bottle_pressure_margin_pa: float = 0.0
    maximum_tank_pressure_pa: float | None = None
    model_id: str = "regulated-single-node-feed"
    schema: str = FEED_SYSTEM_SCHEMA

    def __post_init__(self) -> None:
        coerce_fields(self, "Feed configuration")
        if self.schema != FEED_SYSTEM_SCHEMA:
            raise DomainError("Unsupported regulated-feed schema.")
        if not self.model_id.strip():
            raise DomainError("Feed-system model identifier must not be empty.")
        for value, label in (
            (self.tank_internal_volume_m3, "Tank internal volume"),
            (self.propellant_density_kg_m3, "Propellant density"),
            (self.bottle_volume_m3, "Pressurant bottle volume"),
            (self.regulator_set_pressure_pa, "Regulator set pressure"),
            (self.regulator_full_open_error_pa, "Regulator full-open error"),
            (self.regulator_maximum_area_m2, "Regulator maximum area"),
            (self.valve_time_constant_s, "Valve time constant"),
            (self.tank_wall_temperature_k, "Tank wall temperature"),
            (self.bottle_wall_temperature_k, "Bottle wall temperature"),
        ):
            _positive(value, label)
        if (
            not isfinite(self.regulator_discharge_coefficient)
            or not 0.0 < self.regulator_discharge_coefficient <= 1.0
        ):
            raise DomainError("Regulator discharge coefficient must be in (0, 1].")
        for value, label in (
            (self.feed_line_resistance_pa_s2_kg2, "Feed-line resistance"),
            (self.tank_heat_transfer_w_k, "Tank heat-transfer conductance"),
            (self.bottle_heat_transfer_w_k, "Bottle heat-transfer conductance"),
            (self.minimum_propellant_mass_kg, "Minimum propellant mass"),
            (self.minimum_injector_pressure_pa, "Minimum injector pressure"),
            (
                self.minimum_bottle_pressure_margin_pa,
                "Minimum bottle pressure margin",
            ),
        ):
            _nonnegative(value, label)
        if self.maximum_tank_pressure_pa is not None:
            _positive(self.maximum_tank_pressure_pa, "Maximum tank pressure")
            if self.maximum_tank_pressure_pa <= self.regulator_set_pressure_pa:
                raise DomainError("Maximum tank pressure must exceed regulator set pressure.")


@dataclass(frozen=True, slots=True)
class RegulatedFeedState:
    """Six additional states propagated with the spacecraft.

    Parameters
    ----------
    propellant_mass_kg : float
        Positive liquid mass in the tank, kg.
    tank_pressurant_mass_kg : float
        Positive pressurant mass in the tank ullage, kg.
    tank_gas_temperature_k : float
        Positive ullage gas temperature, K.
    bottle_pressurant_mass_kg : float
        Positive pressurant mass in the bottle, kg.
    bottle_gas_temperature_k : float
        Positive bottle gas temperature, K.
    regulator_opening : float
        Regulator opening fraction in ``[0, 1]``.

    Raises
    ------
    DomainError
        If a mass or temperature is not finite and positive, or the opening
        lies outside ``[0, 1]``.

    Notes
    -----
    Masses and temperatures are the conserved node variables of a lumped
    feed network [1]_; pressures are derived algebraically from them [2]_.

    References
    ----------
    .. [1] NASA, "Nodal Modeling of Liquid Propellant Feed and Pressurization
           System", 2024. https://ntrs.nasa.gov/citations/20240003493
    .. [2] NASA SP-8112. NASA, "Pressurization Systems for Liquid Rockets", 1975.
           https://ntrs.nasa.gov/citations/19760015212
    """

    propellant_mass_kg: float
    tank_pressurant_mass_kg: float
    tank_gas_temperature_k: float
    bottle_pressurant_mass_kg: float
    bottle_gas_temperature_k: float
    regulator_opening: float

    def __post_init__(self) -> None:
        coerce_fields(self, "Feed state")
        for value, label in (
            (self.propellant_mass_kg, "Propellant mass"),
            (self.tank_pressurant_mass_kg, "Tank pressurant mass"),
            (self.tank_gas_temperature_k, "Tank gas temperature"),
            (self.bottle_pressurant_mass_kg, "Bottle pressurant mass"),
            (self.bottle_gas_temperature_k, "Bottle gas temperature"),
        ):
            _positive(value, label)
        if not isfinite(self.regulator_opening) or not 0.0 <= self.regulator_opening <= 1.0:
            raise DomainError("Regulator opening must be finite and in [0, 1].")


@dataclass(frozen=True, slots=True)
class RegulatedFeedDerivative:
    """Right-hand side of the feed ODE in `RegulatedFeedState` field order.

    Attributes
    ----------
    propellant_mass_kg_s : float
        Liquid mass rate, the negative of the demanded flow, kg/s.
    tank_pressurant_mass_kg_s : float
        Regulator inflow to the tank ullage, kg/s.
    tank_gas_temperature_k_s : float
        Ullage gas temperature rate from the tank energy balance [1]_, K/s.
    bottle_pressurant_mass_kg_s : float
        Bottle mass rate, the negative of the regulator flow, kg/s.
    bottle_gas_temperature_k_s : float
        Bottle gas temperature rate from the bottle energy balance, K/s.
    regulator_opening_s : float
        Opening rate ``(command - opening) / tau``, 1/s.

    References
    ----------
    .. [1] NASA, "Nodal Modeling of Liquid Propellant Feed and Pressurization
           System", 2024. https://ntrs.nasa.gov/citations/20240003493
    """

    propellant_mass_kg_s: float
    tank_pressurant_mass_kg_s: float
    tank_gas_temperature_k_s: float
    bottle_pressurant_mass_kg_s: float
    bottle_gas_temperature_k_s: float
    regulator_opening_s: float

    def __post_init__(self) -> None:
        """Coerce every rate to a finite real; raise `DomainError` otherwise."""
        coerce_fields(self, "Feed derivative")

    def as_tuple(self) -> tuple[float, ...]:
        """Return the six node rates [1]_ in `RegulatedFeedState` field order.

        References
        ----------
        .. [1] NASA, "Nodal Modeling of Liquid Propellant Feed and Pressurization
               System", 2024. https://ntrs.nasa.gov/citations/20240003493
        """

        return (
            self.propellant_mass_kg_s,
            self.tank_pressurant_mass_kg_s,
            self.tank_gas_temperature_k_s,
            self.bottle_pressurant_mass_kg_s,
            self.bottle_gas_temperature_k_s,
            self.regulator_opening_s,
        )


@dataclass(frozen=True, slots=True)
class FeedPressurePartials:
    """Analytic injector-pressure derivatives for Jacobian composition.

    Attributes
    ----------
    injector_pressure_wrt_state : tuple of (str, float)
        Derivative of injector pressure with respect to each named state.
        Through the ideal-gas ullage pressure ``p_t = m_t R T_t / V_u`` [2]_
        only the propellant mass, tank pressurant mass and tank temperature
        contribute.
    injector_pressure_wrt_propellant_flow : float
        ``d p_inj / d mdot = -2 R_line mdot``, Pa s/kg.

    Notes
    -----
    These partials feed the coupled Jacobian of the lumped network [1]_.

    References
    ----------
    .. [1] NASA, "Nodal Modeling of Liquid Propellant Feed and Pressurization
           System", 2024. https://ntrs.nasa.gov/citations/20240003493
    .. [2] NASA-GAS. NASA Glenn Research Center, "Equation of State (Ideal Gas)".
           https://www1.grc.nasa.gov/beginners-guide-to-aeronautics/equation-of-state-ideal-gas-2/
    """

    injector_pressure_wrt_state: tuple[tuple[str, float], ...]
    injector_pressure_wrt_propellant_flow: float

    def __post_init__(self) -> None:
        """Coerce entries to finite reals and detach the named-pair sequence."""
        coerce_fields(self, "Feed pressure partials")
        object.__setattr__(self, "injector_pressure_wrt_state", name_value_pairs(
            self.injector_pressure_wrt_state, "Injector-pressure state partials"))


@dataclass(frozen=True, slots=True)
class RegulatedFeedEvaluation:
    """Algebraic feed variables and the simultaneous ODE right-hand side.

    Attributes
    ----------
    tank_ullage_volume_m3 : float
        Ullage volume ``V_tank - m_p / rho`` [1]_, m^3.
    tank_pressure_pa : float
        Ideal-gas ullage pressure, Pa.
    bottle_pressure_pa : float
        Ideal-gas bottle pressure, Pa.
    injector_pressure_pa : float
        Tank pressure minus the quadratic line loss, Pa.
    feed_line_pressure_drop_pa : float
        ``R_line * mdot**2``, Pa.
    regulator_command : float
        Proportional command clipped to ``[0, 1]``.
    regulator_mass_flow_kg_s : float
        Bottle-to-tank pressurant flow, kg/s.
    regulator_choked : bool
        Whether the regulator flow is choked.
    tank_heat_flow_w : float
        Wall-to-ullage heat flow, W.
    bottle_heat_flow_w : float
        Wall-to-bottle-gas heat flow, W.
    derivative : RegulatedFeedDerivative
        Six-state ODE right-hand side [2]_.
    pressure_partials : FeedPressurePartials
        Analytic injector-pressure derivatives.
    reverse_flow_event_extension : bool, default False
        ``True`` when a negative liquid flow was evaluated as a signed
        continuation for zero-flow event location only.

    References
    ----------
    .. [1] NASA SP-8112. NASA, "Pressurization Systems for Liquid Rockets", 1975.
           https://ntrs.nasa.gov/citations/19760015212
    .. [2] NASA, "Nodal Modeling of Liquid Propellant Feed and Pressurization
           System", 2024. https://ntrs.nasa.gov/citations/20240003493
    """

    tank_ullage_volume_m3: float
    tank_pressure_pa: float
    bottle_pressure_pa: float
    injector_pressure_pa: float
    feed_line_pressure_drop_pa: float
    regulator_command: float
    regulator_mass_flow_kg_s: float
    regulator_choked: bool
    tank_heat_flow_w: float
    bottle_heat_flow_w: float
    derivative: RegulatedFeedDerivative
    pressure_partials: FeedPressurePartials
    reverse_flow_event_extension: bool = False

    def __post_init__(self) -> None:
        """Coerce every algebraic output to a finite value; raise otherwise."""
        coerce_fields(self, "Feed evaluation")


class FeedEventKind(str, Enum):
    """Protective cutoff limits owned by the feed system.

    Members are the propellant reserve, minimum injector pressure, minimum
    bottle-to-tank pressure margin and, when configured, tank overpressure;
    [1]_ and [2]_ give the pressurization and regulator design context for
    these limits.

    References
    ----------
    .. [1] NASA SP-8112. NASA, "Pressurization Systems for Liquid Rockets", 1975.
           https://ntrs.nasa.gov/citations/19760015212
    .. [2] NASA SP-8080. NASA, "Liquid Rocket Pressure Regulators and Valves".
           https://ntrs.nasa.gov/citations/19740002611
    """

    PROPELLANT_RESERVE = "propellant_reserve"
    INJECTOR_PRESSURE = "injector_pressure"
    BOTTLE_PRESSURE_MARGIN = "bottle_pressure_margin"
    TANK_OVERPRESSURE = "tank_overpressure"


@dataclass(frozen=True, slots=True)
class FeedSystemEventSurface:
    """Signed event function for one feed-system limit, positive inside it.

    Attributes
    ----------
    event_id : str
        Stable event identifier.
    kind : FeedEventKind
        Limit that the event function measures.
    stops_burn : bool, default True
        Whether a zero crossing terminates the burn.

    Notes
    -----
    The sign convention suits the zero-crossing event detectors of a
    propagator [2]_; the limits come from the pressurization-system design
    [1]_.

    References
    ----------
    .. [1] NASA SP-8112. NASA, "Pressurization Systems for Liquid Rockets", 1975.
           https://ntrs.nasa.gov/citations/19760015212
    .. [2] Orekit Project, "PropulsionModel / EventDetector" API documentation,
           Orekit 13.1.5.
    """

    event_id: str
    kind: FeedEventKind
    stops_burn: bool = True


@dataclass(frozen=True, slots=True)
class FeedCouplingEvidence:
    """Fixed-point convergence evidence for the feed/performance closure.

    Attributes
    ----------
    method : str
        Identifier of the iteration scheme.
    iterations : int
        Iterations performed, including the accepted one.
    converged : bool
        Whether the flow residual met the tolerance.
    relative_tolerance : float
        Relative flow tolerance (scale floor 1 kg/s).
    relaxation : float
        Under-relaxation factor in ``(0, 1]``.
    mass_flow_residual_kg_s : float
        Absolute demanded-minus-assumed flow residual, kg/s.
    injector_pressure_change_pa : float
        Injector-pressure change over the last iteration, Pa.

    Notes
    -----
    [1]_ and [2]_ give the transient and nodal feed-network modelling context
    of the closure.

    References
    ----------
    .. [1] NASA, "Transient Mathematical Modeling for Liquid Rocket Engine
           Systems". https://ntrs.nasa.gov/citations/20040000363
    .. [2] NASA, "Nodal Modeling of Liquid Propellant Feed and Pressurization
           System", 2024. https://ntrs.nasa.gov/citations/20240003493
    """

    method: str
    iterations: int
    converged: bool
    relative_tolerance: float
    relaxation: float
    mass_flow_residual_kg_s: float
    injector_pressure_change_pa: float


@dataclass(frozen=True, slots=True)
class CoupledFeedPartials:
    """Implicit derivatives through the feed-pressure / engine-flow loop.

    Attributes
    ----------
    algebraic_denominator : float
        ``1 - P_q Q_p``; its magnitude measures the distance from an algebraic
        singularity of the loop.
    injector_pressure_wrt_feed_state : tuple of (str, float)
        Closed derivative ``dp/dx`` for each feed state.
    propellant_flow_wrt_feed_state : tuple of (str, float)
        ``Q_p dp/dx`` for each feed state.
    acceleration_wrt_feed_state : tuple of (str, (float, float, float))
        Thrust-acceleration derivative for each feed state.

    Notes
    -----
    If engine demand is ``q = Q(p)`` and injector pressure is ``p = P(x, q)``,
    the closed derivative is

    ``dp/dx = P_x / (1 - P_q Q_p)``.

    The denominator is reported explicitly so that proximity to an algebraic
    singularity is visible to the caller rather than implicit in the
    fixed-point iteration.  [1]_ and [2]_ give the coupled-component
    modelling context.

    References
    ----------
    .. [1] NASA, "Transient Mathematical Modeling for Liquid Rocket Engine
           Systems". https://ntrs.nasa.gov/citations/20040000363
    .. [2] NASA, "Nodal Modeling of Liquid Propellant Feed and Pressurization
           System", 2024. https://ntrs.nasa.gov/citations/20240003493
    """

    algebraic_denominator: float
    injector_pressure_wrt_feed_state: tuple[tuple[str, float], ...]
    propellant_flow_wrt_feed_state: tuple[tuple[str, float], ...]
    acceleration_wrt_feed_state: tuple[
        tuple[str, tuple[float, float, float]], ...
    ]


@dataclass(frozen=True, slots=True)
class CoupledFeedPropulsionEvaluation:
    """Closed feed-system and propagator-force evaluation at one stage.

    Attributes
    ----------
    propulsion : PropulsionDynamicsEvaluation
        Engine evaluation at the converged injector pressure.
    feed : RegulatedFeedEvaluation
        Feed evaluation at the converged flow.
    evidence : FeedCouplingEvidence
        Convergence evidence.
    partials : CoupledFeedPartials
        Implicit closed-loop derivatives.

    Notes
    -----
    [1]_ and [2]_ give the coupled transient and feed-network modelling
    context.

    References
    ----------
    .. [1] NASA, "Transient Mathematical Modeling for Liquid Rocket Engine
           Systems". https://ntrs.nasa.gov/citations/20040000363
    .. [2] NASA, "Nodal Modeling of Liquid Propellant Feed and Pressurization
           System", 2024. https://ntrs.nasa.gov/citations/20240003493
    """

    propulsion: PropulsionDynamicsEvaluation
    feed: RegulatedFeedEvaluation
    evidence: FeedCouplingEvidence
    partials: CoupledFeedPartials


@dataclass(frozen=True, slots=True)
class RegulatedFeedSystem:
    """Single-ullage-node regulated feed model for propagation.

    Parameters
    ----------
    configuration : RegulatedFeedConfiguration
        Geometry, control, thermal and cutoff configuration.

    Notes
    -----
    The tank energy equation includes inlet enthalpy, ullage boundary work
    and wall heat transfer; the fixed-volume bottle includes outlet enthalpy
    and wall heat transfer [3]_.  Pressurant mass is conserved exactly between
    the two nodes.  [1]_ and [2]_ give the pressurization-system and regulator
    design context.

    The model omits line inertance, acoustic and water-hammer modes,
    stratified ullage, phase change, dissolved gas, regulator hysteresis and
    real-fluid properties.  Fast pressure events may require a higher-order
    validated nodal model [3]_, [4]_.

    References
    ----------
    .. [1] NASA SP-8112. NASA, "Pressurization Systems for Liquid Rockets", 1975.
           https://ntrs.nasa.gov/citations/19760015212
    .. [2] NASA SP-8080. NASA, "Liquid Rocket Pressure Regulators and Valves".
           https://ntrs.nasa.gov/citations/19740002611
    .. [3] NASA, "Nodal Modeling of Liquid Propellant Feed and Pressurization
           System", 2024. https://ntrs.nasa.gov/citations/20240003493
    .. [4] NASA NESC Technical Bulletin 22-03, "Treatment of Transient Pressure
           Events", 2022. https://ntrs.nasa.gov/citations/20220006583
    """

    configuration: RegulatedFeedConfiguration

    @property
    def model_hash(self) -> str:
        """SHA-256 binding configuration, event policy and source identifiers.

        The digest covers the manifest without its hash field, for model
        traceability [1]_; it is not a hash of the executable implementation.

        References
        ----------
        .. [1] NASA-STD-7009B. NASA, "Standard for Models and Simulations", 2024.
               https://standards.nasa.gov/standard/NASA/NASA-STD-7009
        """

        return canonical_sha256(self.manifest(include_hash=False))

    def manifest(self, *, include_hash: bool = True) -> dict[str, object]:
        """Return the deterministic integration manifest.

        Parameters
        ----------
        include_hash : bool, default True
            Whether to append ``model_hash`` computed over the other fields.

        Returns
        -------
        dict
            Schema, model identifier, gas properties, configuration values,
            state order, source identifiers and the reverse-flow event policy.

        Notes
        -----
        The manifest records the node configuration of [2]_ for model
        traceability [1]_.

        References
        ----------
        .. [1] NASA-STD-7009B. NASA, "Standard for Models and Simulations", 2024.
               https://standards.nasa.gov/standard/NASA/NASA-STD-7009
        .. [2] NASA, "Nodal Modeling of Liquid Propellant Feed and Pressurization
               System", 2024. https://ntrs.nasa.gov/citations/20240003493
        """

        config = self.configuration
        payload: dict[str, object] = {
            "schema": config.schema,
            "model_id": config.model_id,
            "gas": {
                "name": config.gas.name,
                "specific_gas_constant_j_kg_k": config.gas.specific_gas_constant_j_kg_k,
                "gamma": config.gas.gamma,
            },
            "configuration": {
                name: getattr(config, name)
                for name in (
                    "tank_internal_volume_m3",
                    "propellant_density_kg_m3",
                    "bottle_volume_m3",
                    "regulator_set_pressure_pa",
                    "regulator_full_open_error_pa",
                    "regulator_maximum_area_m2",
                    "regulator_discharge_coefficient",
                    "valve_time_constant_s",
                    "feed_line_resistance_pa_s2_kg2",
                    "tank_heat_transfer_w_k",
                    "tank_wall_temperature_k",
                    "bottle_heat_transfer_w_k",
                    "bottle_wall_temperature_k",
                    "minimum_propellant_mass_kg",
                    "minimum_injector_pressure_pa",
                    "minimum_bottle_pressure_margin_pa",
                    "maximum_tank_pressure_pa",
                )
            },
            "state_order": [
                "propellant_mass_kg",
                "tank_pressurant_mass_kg",
                "tank_gas_temperature_k",
                "bottle_pressurant_mass_kg",
                "bottle_gas_temperature_k",
                "regulator_opening",
            ],
            "reference_ids": list(FEED_SYSTEM_REFERENCE_IDS),
            "event_location_policy": {
                "negative_liquid_flow": "signed_continuation_with_zero_static_resistance",
                "physical_cutoff": "zero_liquid_flow",
                "sustained_reverse_operation": "unsupported",
            },
        }
        if include_hash:
            payload["model_hash"] = canonical_sha256(payload)
        return payload

    def initial_state(
        self,
        *,
        propellant_mass_kg: float,
        tank_pressure_pa: float,
        tank_gas_temperature_k: float,
        bottle_pressure_pa: float,
        bottle_gas_temperature_k: float,
        regulator_opening: float = 0.0,
    ) -> RegulatedFeedState:
        """Build a mass-consistent state from declared node pressures.

        Parameters
        ----------
        propellant_mass_kg : float
            Positive liquid mass leaving positive ullage, kg.
        tank_pressure_pa : float
            Positive initial ullage pressure, Pa.
        tank_gas_temperature_k : float
            Positive initial ullage temperature, K.
        bottle_pressure_pa : float
            Positive initial bottle pressure, Pa.
        bottle_gas_temperature_k : float
            Positive initial bottle temperature, K.
        regulator_opening : float, default 0.0
            Initial opening fraction in ``[0, 1]``.

        Returns
        -------
        RegulatedFeedState
            State whose pressurant masses are ``m = p V / (R T)`` [3]_ for the
            ullage volume ``V_tank - m_p / rho`` [1]_ and the bottle volume.

        Raises
        ------
        DomainError
            If an input is outside its domain or the propellant leaves no
            positive ullage.

        Notes
        -----
        Node masses rather than pressures are the conserved states of the
        lumped network [2]_.

        References
        ----------
        .. [1] NASA SP-8112. NASA, "Pressurization Systems for Liquid Rockets", 1975.
               https://ntrs.nasa.gov/citations/19760015212
        .. [2] NASA, "Nodal Modeling of Liquid Propellant Feed and Pressurization
               System", 2024. https://ntrs.nasa.gov/citations/20240003493
        .. [3] NASA-GAS. NASA Glenn Research Center, "Equation of State (Ideal Gas)",
               Beginner's Guide to Aeronautics, accessed 2026-10-18.
        """

        propellant_mass_kg = real(propellant_mass_kg, "Initial propellant mass")
        tank_pressure_pa = real(tank_pressure_pa, "Initial tank pressure")
        tank_gas_temperature_k = real(tank_gas_temperature_k, "Initial tank gas temperature")
        bottle_pressure_pa = real(bottle_pressure_pa, "Initial bottle pressure")
        bottle_gas_temperature_k = real(bottle_gas_temperature_k, "Initial bottle gas temperature")
        regulator_opening = real(regulator_opening, "Initial regulator opening")

        _positive(propellant_mass_kg, "Initial propellant mass")
        _positive(tank_pressure_pa, "Initial tank pressure")
        _positive(tank_gas_temperature_k, "Initial tank gas temperature")
        _positive(bottle_pressure_pa, "Initial bottle pressure")
        _positive(bottle_gas_temperature_k, "Initial bottle gas temperature")
        ullage = self._ullage_volume(propellant_mass_kg)
        gas_constant = self.configuration.gas.specific_gas_constant_j_kg_k
        return RegulatedFeedState(
            propellant_mass_kg=propellant_mass_kg,
            tank_pressurant_mass_kg=(
                tank_pressure_pa * ullage / (gas_constant * tank_gas_temperature_k)
            ),
            tank_gas_temperature_k=tank_gas_temperature_k,
            bottle_pressurant_mass_kg=(
                bottle_pressure_pa
                * self.configuration.bottle_volume_m3
                / (gas_constant * bottle_gas_temperature_k)
            ),
            bottle_gas_temperature_k=bottle_gas_temperature_k,
            regulator_opening=regulator_opening,
        )

    def _ullage_volume(self, propellant_mass_kg: float) -> float:
        """Return the rigid-tank ullage volume ``V_tank - m_p / rho`` [1]_, m^3.

        Raises `DomainError` unless the propellant leaves a positive,
        representable ullage below tank capacity.

        References
        ----------
        .. [1] NASA SP-8112. NASA, "Pressurization Systems for Liquid Rockets", 1975.
               https://ntrs.nasa.gov/citations/19760015212
        """
        config = self.configuration
        maximum_propellant = (
            config.tank_internal_volume_m3 * config.propellant_density_kg_m3
        )
        if not isfinite(propellant_mass_kg) or not 0.0 < propellant_mass_kg < maximum_propellant:
            raise DomainError(
                "Propellant mass must leave positive ullage and remain below tank capacity."
            )
        ullage = config.tank_internal_volume_m3 - (
            propellant_mass_kg / config.propellant_density_kg_m3
        )
        _positive(ullage, "Representable positive tank ullage")
        return ullage

    def _pressures(self, state: RegulatedFeedState) -> tuple[float, float, float]:
        """Return ``(ullage_m3, tank_pa, bottle_pa)`` from the ideal-gas law.

        Each node pressure is ``p = m R T / V`` [1]_; non-positive or
        non-representable pressures raise `DomainError`.

        References
        ----------
        .. [1] NASA-GAS. NASA Glenn Research Center, "Equation of State (Ideal Gas)",
               Beginner's Guide to Aeronautics, accessed 2026-10-18.
        """
        config = self.configuration
        ullage = self._ullage_volume(state.propellant_mass_kg)
        gas_constant = config.gas.specific_gas_constant_j_kg_k
        tank_pressure = (
            state.tank_pressurant_mass_kg
            * gas_constant
            * state.tank_gas_temperature_k
            / ullage
        )
        bottle_pressure = (
            state.bottle_pressurant_mass_kg
            * gas_constant
            * state.bottle_gas_temperature_k
            / config.bottle_volume_m3
        )
        _positive(tank_pressure, "Representable positive tank pressure")
        _positive(bottle_pressure, "Representable positive bottle pressure")
        return ullage, tank_pressure, bottle_pressure

    def _regulator_flow(
        self,
        *,
        upstream_pressure_pa: float,
        upstream_temperature_k: float,
        downstream_pressure_pa: float,
        opening: float,
    ) -> tuple[float, bool]:
        """Return the compressible regulator mass flow and whether it is choked.

        With pressure ratio ``r = p_down / p_up``, area
        ``A = A_max * opening`` and ``G = Cd A p_up / sqrt(R T_up)``, the flow
        is [1]_

        ``mdot = G sqrt(2 gamma / (gamma - 1) * (r**(2/gamma) - r**((gamma+1)/gamma)))``

        for ``r`` above the critical ratio ``(2 / (gamma + 1))**(gamma/(gamma-1))``
        and the choked value
        ``mdot = G sqrt(gamma) (1 + (gamma - 1) / 2)**(-(gamma+1) / (2 (gamma-1)))``
        otherwise.  ``log r`` is formed with ``log1p`` near equilibrium and the
        unchoked bracket is evaluated as ``r**(2/gamma) * -expm1(...)`` to
        avoid cancellation.  Zero flow is returned when the regulator is
        closed or the pressure difference is not forward.

        References
        ----------
        .. [1] NASA Glenn Research Center, "Mass Flow Choking".
               https://www.grc.nasa.gov/www/k-12/BGP/mflchk.html
        """
        config = self.configuration
        gas = config.gas
        if upstream_pressure_pa <= downstream_pressure_pa or opening <= 0.0:
            return 0.0, False
        area = config.regulator_maximum_area_m2 * opening
        delta_gamma = gas.gamma - 1.0
        # log1p keeps the exact pressure gap near equilibrium; rounding
        # p_down/p_up first and subtracting two powers would lose the flow.
        log_ratio = (log1p((downstream_pressure_pa - upstream_pressure_pa) / upstream_pressure_pa)
                     if downstream_pressure_pa > .5 * upstream_pressure_pa else
                     log(downstream_pressure_pa) - log(upstream_pressure_pa))
        critical_log = -(gas.gamma / delta_gamma) * log1p(.5 * delta_gamma)
        common = (
            config.regulator_discharge_coefficient
            * area
            * upstream_pressure_pa
            / (sqrt(gas.specific_gas_constant_j_kg_k) * sqrt(upstream_temperature_k))
        )
        choked = log_ratio <= critical_log
        if choked:
            factor = sqrt(gas.gamma) * exp(
                -(.5 + 1.0 / delta_gamma) * log1p(.5 * delta_gamma))
        else:
            # r^(2/gamma) * (1-r^((gamma-1)/gamma)); no near-equal power subtraction.
            radicand = (2.0 * (gas.gamma / delta_gamma) * exp(2.0 * log_ratio / gas.gamma)
                        * -expm1((delta_gamma / gas.gamma) * log_ratio))
            factor = sqrt(radicand)
        flow = real(common * factor, "Regulator mass flow")
        _positive(flow, "Forward regulator flow in the representable numerical domain")
        return flow, choked

    def evaluate(
        self,
        state: RegulatedFeedState,
        propellant_mass_flow_kg_s: float,
        *,
        for_event_location: bool = False,
    ) -> RegulatedFeedEvaluation:
        """Return the algebraic pressures and the six-state ODE right-hand side.

        Parameters
        ----------
        state : RegulatedFeedState
            Current feed state.
        propellant_mass_flow_kg_s : float
            Liquid demand; positive values drain the tank, kg/s.
        for_event_location : bool, default False
            Permit a negative demand as a signed continuation for locating a
            zero-flow cutoff.

        Returns
        -------
        RegulatedFeedEvaluation
            Pressures, regulator flow, heat flows, ODE right-hand side and
            injector-pressure partials.

        Raises
        ------
        DomainError
            If ``state`` has the wrong type, ``for_event_location`` is not a
            bool, a negative demand is supplied without
            ``for_event_location`` or with a non-zero static line resistance,
            or a derived quantity is not finite and positive.

        Notes
        -----
        Node pressures follow the ideal-gas law on the ullage
        ``V_u = V_tank - m_p / rho`` [1]_ and the bottle.  The gas-node energy
        balances [2]_ are

        ``d(m_t cv T_t)/dt = mdot_r cp T_b + Q_t - p_t dV_u/dt``

        and

        ``d(m_b cv T_b)/dt = -mdot_r cp T_b + Q_b``,

        with wall heat flows ``Q = hA (T_wall - T_gas)``.  The regulator flow
        ``mdot_r`` follows the compressible-restriction relation with choking
        [3]_, and the injector pressure is ``p_t - R_line mdot**2``.

        Negative liquid flow is refused by default.  ``for_event_location``
        permits a signed continuation, only with zero static line resistance,
        for locating the dynamic line's zero-flow cutoff.  Tank mass and ullage
        work retain the signed flow so that the trial right-hand side conserves
        mass and gas energy.  Such a result is marked
        ``reverse_flow_event_extension``; it does not describe sustained
        reverse operation, and physical propagation terminates at the zero-flow
        event.

        References
        ----------
        .. [1] NASA SP-8112. NASA, "Pressurization Systems for Liquid Rockets", 1975.
               https://ntrs.nasa.gov/citations/19760015212
        .. [2] NASA, "Nodal Modeling of Liquid Propellant Feed and Pressurization
               System", 2024. https://ntrs.nasa.gov/citations/20240003493
        .. [3] NASA Glenn Research Center, "Mass Flow Choking".
               https://www.grc.nasa.gov/www/k-12/BGP/mflchk.html
        """

        if not isinstance(state, RegulatedFeedState):
            raise DomainError("Feed evaluation requires a RegulatedFeedState.")
        propellant_mass_flow_kg_s = real(propellant_mass_flow_kg_s, "Propellant mass flow")

        config = self.configuration
        if not isinstance(for_event_location, bool):
            raise DomainError("for_event_location must be a boolean.")
        if propellant_mass_flow_kg_s < 0.0:
            if not for_event_location:
                raise DomainError("Propellant mass-flow demand must be non-negative.")
            if config.feed_line_resistance_pa_s2_kg2 != 0.0:
                raise DomainError("Reverse event continuation requires zero static line resistance.")
        gas = config.gas
        ullage, tank_pressure, bottle_pressure = self._pressures(state)
        pressure_error = config.regulator_set_pressure_pa - tank_pressure
        command = min(1.0, max(0.0, pressure_error / config.regulator_full_open_error_pa))
        regulator_flow, choked = self._regulator_flow(
            upstream_pressure_pa=bottle_pressure,
            upstream_temperature_k=state.bottle_gas_temperature_k,
            downstream_pressure_pa=tank_pressure,
            opening=state.regulator_opening,
        )
        ullage_rate = propellant_mass_flow_kg_s / config.propellant_density_kg_m3
        tank_heat = config.tank_heat_transfer_w_k * (
            config.tank_wall_temperature_k - state.tank_gas_temperature_k
        )
        bottle_heat = config.bottle_heat_transfer_w_k * (
            config.bottle_wall_temperature_k - state.bottle_gas_temperature_k
        )
        cv = gas.specific_heat_cv_j_kg_k
        tank_capacity = state.tank_pressurant_mass_kg * cv
        bottle_capacity = state.bottle_pressurant_mass_kg * cv
        _positive(tank_capacity, "Tank gas thermal capacity")
        _positive(bottle_capacity, "Bottle gas thermal capacity")
        inlet_thermal_flow = (regulator_flow * gas.specific_gas_constant_j_kg_k
                              * state.bottle_gas_temperature_k)
        tank_temperature_rate = (
            regulator_flow * (state.bottle_gas_temperature_k - state.tank_gas_temperature_k)
            / state.tank_pressurant_mass_kg
        ) + (
            inlet_thermal_flow
            + tank_heat
            - tank_pressure * ullage_rate
        ) / tank_capacity
        bottle_temperature_rate = (
            bottle_heat - inlet_thermal_flow
        ) / bottle_capacity
        pressure_drop = (
            config.feed_line_resistance_pa_s2_kg2 * propellant_mass_flow_kg_s
            * propellant_mass_flow_kg_s
        )
        injector_pressure = tank_pressure - pressure_drop
        derivative = RegulatedFeedDerivative(
            propellant_mass_kg_s=-propellant_mass_flow_kg_s,
            tank_pressurant_mass_kg_s=regulator_flow,
            tank_gas_temperature_k_s=tank_temperature_rate,
            bottle_pressurant_mass_kg_s=-regulator_flow,
            bottle_gas_temperature_k_s=bottle_temperature_rate,
            regulator_opening_s=(command - state.regulator_opening)
            / config.valve_time_constant_s,
        )
        gas_constant = gas.specific_gas_constant_j_kg_k
        pressure_partials = FeedPressurePartials(
            injector_pressure_wrt_state=(
                (
                    "propellant_mass_kg",
                    state.tank_pressurant_mass_kg
                    * gas_constant
                    * state.tank_gas_temperature_k
                    / (config.propellant_density_kg_m3 * ullage**2),
                ),
                (
                    "tank_pressurant_mass_kg",
                    gas_constant * state.tank_gas_temperature_k / ullage,
                ),
                (
                    "tank_gas_temperature_k",
                    state.tank_pressurant_mass_kg * gas_constant / ullage,
                ),
                ("bottle_pressurant_mass_kg", 0.0),
                ("bottle_gas_temperature_k", 0.0),
                ("regulator_opening", 0.0),
            ),
            injector_pressure_wrt_propellant_flow=(
                -2.0 * (config.feed_line_resistance_pa_s2_kg2 * propellant_mass_flow_kg_s)
            ),
        )
        return RegulatedFeedEvaluation(
            tank_ullage_volume_m3=ullage,
            tank_pressure_pa=tank_pressure,
            bottle_pressure_pa=bottle_pressure,
            injector_pressure_pa=injector_pressure,
            feed_line_pressure_drop_pa=pressure_drop,
            regulator_command=command,
            regulator_mass_flow_kg_s=regulator_flow,
            regulator_choked=choked,
            tank_heat_flow_w=tank_heat,
            bottle_heat_flow_w=bottle_heat,
            derivative=derivative,
            pressure_partials=pressure_partials,
            reverse_flow_event_extension=propellant_mass_flow_kg_s < 0.0,
        )

    def event_surfaces(self) -> tuple[FeedSystemEventSurface, ...]:
        """Return every configured feed cutoff event.

        Returns
        -------
        tuple of FeedSystemEventSurface
            Propellant-reserve, injector-pressure and bottle-margin events,
            plus tank overpressure when ``maximum_tank_pressure_pa`` is set.
            [1]_ and [2]_ give the design context of these limits.

        References
        ----------
        .. [1] NASA SP-8112. NASA, "Pressurization Systems for Liquid Rockets", 1975.
               https://ntrs.nasa.gov/citations/19760015212
        .. [2] NASA SP-8080. NASA, "Liquid Rocket Pressure Regulators and Valves".
               https://ntrs.nasa.gov/citations/19740002611
        """

        surfaces = [
            FeedSystemEventSurface(
                "feed:propellant-reserve", FeedEventKind.PROPELLANT_RESERVE
            ),
            FeedSystemEventSurface(
                "feed:injector-pressure", FeedEventKind.INJECTOR_PRESSURE
            ),
            FeedSystemEventSurface(
                "feed:bottle-pressure-margin", FeedEventKind.BOTTLE_PRESSURE_MARGIN
            ),
        ]
        if self.configuration.maximum_tank_pressure_pa is not None:
            surfaces.append(
                FeedSystemEventSurface(
                    "feed:tank-overpressure", FeedEventKind.TANK_OVERPRESSURE
                )
            )
        return tuple(surfaces)

    def event_value(
        self,
        surface: FeedSystemEventSurface,
        state: RegulatedFeedState,
        propellant_mass_flow_kg_s: float,
    ) -> float:
        """Evaluate one feed cutoff function, positive inside the permitted domain.

        Parameters
        ----------
        surface : FeedSystemEventSurface
            Event to evaluate.
        state : RegulatedFeedState
            Current feed state.
        propellant_mass_flow_kg_s : float
            Non-negative liquid demand, kg/s.

        Returns
        -------
        float
            Reserve margin ``m_p - m_min``, injector margin
            ``p_inj - p_inj_min``, bottle margin ``p_b - p_t - margin_min``, or
            overpressure margin ``p_max - p_t``; a sign change marks the event
            for a zero-crossing detector [2]_.  The limits follow the
            pressurization-system design [1]_.

        Raises
        ------
        DomainError
            If the event kind is unsupported or the feed evaluation fails.

        References
        ----------
        .. [1] NASA SP-8112. NASA, "Pressurization Systems for Liquid Rockets", 1975.
               https://ntrs.nasa.gov/citations/19760015212
        .. [2] Orekit Project, "PropulsionModel / EventDetector" API documentation,
               Orekit 13.1.5.
        """

        config = self.configuration
        evaluation = self.evaluate(state, propellant_mass_flow_kg_s)
        if surface.kind is FeedEventKind.PROPELLANT_RESERVE:
            return state.propellant_mass_kg - config.minimum_propellant_mass_kg
        if surface.kind is FeedEventKind.INJECTOR_PRESSURE:
            return evaluation.injector_pressure_pa - config.minimum_injector_pressure_pa
        if surface.kind is FeedEventKind.BOTTLE_PRESSURE_MARGIN:
            return (
                evaluation.bottle_pressure_pa
                - evaluation.tank_pressure_pa
                - config.minimum_bottle_pressure_margin_pa
            )
        if surface.kind is FeedEventKind.TANK_OVERPRESSURE:
            assert config.maximum_tank_pressure_pa is not None
            return config.maximum_tank_pressure_pa - evaluation.tank_pressure_pa
        raise DomainError(f"Unsupported feed event kind: {surface.kind}.")


def evaluate_coupled_feed_propulsion(
    bridge: PropagatorPropulsionBridge,
    feed_system: RegulatedFeedSystem,
    feed_state: RegulatedFeedState,
    *,
    relative_time_s: float,
    vehicle_mass_kg: float,
    direction: tuple[float, float, float],
    pressure_condition_name: str,
    other_operating_conditions: Mapping[str, float] | None = None,
    additional_states: Mapping[str, float] | None = None,
    parameter_overrides: Mapping[str, float] | None = None,
    initial_mass_flow_kg_s: float = 0.0,
    apply_inventory_limits: bool = True,
    ambient_pressure_pa: float | None = None,
    relative_tolerance: float = 1e-10,
    relaxation: float = 0.5,
    maximum_iterations: int = 100,
) -> CoupledFeedPropulsionEvaluation:
    """Close engine flow demand against the feed-line-dependent inlet pressure.

    Parameters
    ----------
    bridge : PropagatorPropulsionBridge
        Propagator bridge with a performance surface.
    feed_system : RegulatedFeedSystem
        Regulated feed model.
    feed_state : RegulatedFeedState
        Current feed state.
    relative_time_s : float
        Time from burn start, s.
    vehicle_mass_kg : float
        Current vehicle mass, kg.
    direction : tuple of float
        Thrust direction in the propagation frame.
    pressure_condition_name : str
        Non-blank performance-surface axis driven by injector pressure.
    other_operating_conditions : Mapping[str, float] or None, default None
        Remaining operating conditions; must not contain
        ``pressure_condition_name``.
    additional_states : Mapping[str, float] or None, default None
        Additional propagated states passed to the bridge.
    parameter_overrides : Mapping[str, float] or None, default None
        Parameter values passed to the bridge.
    initial_mass_flow_kg_s : float, default 0.0
        Non-negative initial flow guess, kg/s.
    apply_inventory_limits : bool, default True
        Whether the bridge applies inventory cutoffs.
    ambient_pressure_pa : float or None, default None
        Ambient pressure passed to the bridge, Pa.
    relative_tolerance : float, default 1e-10
        Positive relative flow tolerance (scale floor 1 kg/s).
    relaxation : float, default 0.5
        Under-relaxation factor in ``(0, 1]``.
    maximum_iterations : int, default 100
        Positive iteration limit.

    Returns
    -------
    CoupledFeedPropulsionEvaluation
        Converged engine and feed evaluations, convergence evidence and the
        implicit closed-loop partials.

    Raises
    ------
    DomainError
        If an argument is outside its domain, the pressure condition is not a
        performance axis or is also supplied in ``other_operating_conditions``,
        the bridge does not expose the pressure-condition partials, the
        algebraic Jacobian denominator is within ``1e-10`` of zero, an iterate
        is invalid, or the iteration does not converge.

    Notes
    -----
    The performance surface produces tank-flow demand ``Q(p)`` from injector
    pressure, while the feed line produces injector pressure ``P(x, q)`` from
    that flow.  The relaxed fixed-point iteration
    ``q <- q + omega (Q(P(x, q)) - q)`` resolves this algebraic loop and
    returns explicit residual evidence; non-convergence raises `DomainError`.
    At convergence the state derivatives are closed through
    ``dp/dx = P_x / (1 - P_q Q_p)``.  [1]_ and [2]_ give the coupled
    feed-network modelling context and [3]_ a pressurization-control study of
    the feed-pressure / flow coupling.

    References
    ----------
    .. [1] NASA, "Transient Mathematical Modeling for Liquid Rocket Engine
           Systems". https://ntrs.nasa.gov/citations/20040000363
    .. [2] NASA, "Nodal Modeling of Liquid Propellant Feed and Pressurization
           System", 2024. https://ntrs.nasa.gov/citations/20240003493
    .. [3] NASA, "Simulation Analysis of Computer-Controlled Pressurization for
           Mixture Ratio Control", 2005. https://ntrs.nasa.gov/citations/20050206420
    """

    if bridge.performance_surface is None:
        raise DomainError("Coupled feed evaluation requires a bridge performance surface.")
    axis_names = {axis.name for axis in bridge.performance_surface.axes}
    if pressure_condition_name not in axis_names:
        raise DomainError("Pressure condition is not an axis of the bridge performance surface.")
    if not pressure_condition_name.strip():
        raise DomainError("Pressure condition name must not be empty.")
    _nonnegative(initial_mass_flow_kg_s, "Initial coupled mass-flow guess")
    _positive(relative_tolerance, "Coupled relative tolerance")
    if not isfinite(relaxation) or not 0.0 < relaxation <= 1.0:
        raise DomainError("Coupling relaxation must be finite and in (0, 1].")
    if maximum_iterations < 1:
        raise DomainError("Coupling maximum iterations must be positive.")
    conditions = dict(other_operating_conditions or {})
    if pressure_condition_name in conditions:
        raise DomainError("Pressure condition is owned by the feed-system coupling.")

    flow = initial_mass_flow_kg_s
    previous_pressure: float | None = None
    pressure_change = 0.0
    propulsion: PropulsionDynamicsEvaluation | None = None
    feed: RegulatedFeedEvaluation | None = None
    residual = float("inf")
    for iteration in range(1, maximum_iterations + 1):
        feed = feed_system.evaluate(feed_state, flow)
        active_conditions = {
            **conditions,
            pressure_condition_name: feed.injector_pressure_pa,
        }
        propulsion = bridge.evaluate(
            relative_time_s=relative_time_s,
            vehicle_mass_kg=vehicle_mass_kg,
            direction=direction,
            additional_states=additional_states,
            operating_conditions=active_conditions,
            parameter_overrides=parameter_overrides,
            apply_inventory_limits=apply_inventory_limits,
            ambient_pressure_pa=ambient_pressure_pa,
        )
        demanded_flow = -propulsion.mass_derivative_kg_s
        residual = demanded_flow - flow
        pressure_change = (
            0.0
            if previous_pressure is None
            else abs(feed.injector_pressure_pa - previous_pressure)
        )
        scale = max(1.0, abs(demanded_flow), abs(flow))
        if abs(residual) <= relative_tolerance * scale:
            acceleration_condition_partials = dict(
                propulsion.partials.acceleration_wrt_conditions
            )
            mass_rate_condition_partials = dict(
                propulsion.partials.mass_rate_wrt_conditions
            )
            if pressure_condition_name not in acceleration_condition_partials or (
                pressure_condition_name not in mass_rate_condition_partials
            ):
                raise DomainError(
                    "Bridge did not expose the pressure-condition partials required "
                    "for coupled feed linearization."
                )
            acceleration_wrt_pressure = acceleration_condition_partials[
                pressure_condition_name
            ]
            demanded_flow_wrt_pressure = -mass_rate_condition_partials[
                pressure_condition_name
            ]
            injector_pressure_wrt_flow = (
                feed.pressure_partials.injector_pressure_wrt_propellant_flow
            )
            algebraic_denominator = (
                1.0
                - injector_pressure_wrt_flow * demanded_flow_wrt_pressure
            )
            if abs(algebraic_denominator) <= 1e-10:
                raise DomainError(
                    "Feed/performance algebraic Jacobian is singular or ill-conditioned."
                )
            closed_pressure_partials = tuple(
                (name, partial / algebraic_denominator)
                for name, partial in feed.pressure_partials.injector_pressure_wrt_state
            )
            flow_state_partials = tuple(
                (name, demanded_flow_wrt_pressure * partial)
                for name, partial in closed_pressure_partials
            )
            acceleration_state_partials = tuple(
                (
                    name,
                    tuple(component * partial for component in acceleration_wrt_pressure),
                )
                for name, partial in closed_pressure_partials
            )
            evidence = FeedCouplingEvidence(
                method="relaxed_fixed_point_feed_performance_v1",
                iterations=iteration,
                converged=True,
                relative_tolerance=relative_tolerance,
                relaxation=relaxation,
                mass_flow_residual_kg_s=abs(residual),
                injector_pressure_change_pa=pressure_change,
            )
            return CoupledFeedPropulsionEvaluation(
                propulsion,
                feed,
                evidence,
                CoupledFeedPartials(
                    algebraic_denominator=algebraic_denominator,
                    injector_pressure_wrt_feed_state=closed_pressure_partials,
                    propellant_flow_wrt_feed_state=flow_state_partials,
                    acceleration_wrt_feed_state=acceleration_state_partials,
                ),
            )
        previous_pressure = feed.injector_pressure_pa
        flow += relaxation * residual
        if not isfinite(flow) or flow < 0.0:
            raise DomainError("Feed/performance coupling produced an invalid flow iterate.")
    raise DomainError(
        "Feed/performance coupling did not converge: "
        f"iterations={maximum_iterations}, residual_kg_s={abs(residual):.9g}, "
        f"pressure_change_pa={pressure_change:.9g}."
    )

