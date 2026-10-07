"""Propagation-ready regulated pressurant and propellant-feed dynamics.

The model is a compact, auditable zero-dimensional network: a finite
pressurant bottle feeds a growing tank ullage through a dynamic regulator;
liquid withdrawal creates ullage volume and a quadratic feed-line loss sets the
engine inlet pressure.  It is intended as a numerical propagator additional-
state model, not as a replacement for validated nodal or water-hammer tools.
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
    if not isfinite(value) or value <= 0.0:
        raise DomainError(f"{label} must be finite and greater than zero.")


def _nonnegative(value: float, label: str) -> None:
    if not isfinite(value) or value < 0.0:
        raise DomainError(f"{label} must be finite and non-negative.")


def _canonical_hash(payload: object) -> str:
    return hashlib.sha256(
        json.dumps(
            payload,
            allow_nan=False,
            ensure_ascii=False,
            separators=(",", ":"),
            sort_keys=True,
        ).encode("utf-8")
    ).hexdigest()


@dataclass(frozen=True, slots=True)
class IdealPressurantGas:
    """Calorically perfect pressurant properties used by both gas nodes.

    References
    ----------
    NASA-SP-8112: https://ntrs.nasa.gov/citations/19760015212
    NASA nodal feed/pressurization modeling:
    https://ntrs.nasa.gov/citations/20240003493
    """

    specific_gas_constant_j_kg_k: float
    gamma: float
    name: str = "pressurant"

    def __post_init__(self) -> None:
        _positive(self.specific_gas_constant_j_kg_k, "Pressurant gas constant")
        if not isfinite(self.gamma) or self.gamma <= 1.0:
            raise DomainError("Pressurant gamma must be finite and greater than one.")
        if not self.name.strip():
            raise DomainError("Pressurant name must not be empty.")

    @property
    def specific_heat_cv_j_kg_k(self) -> float:
        """Return ``cv = R/(gamma-1)``.

        References
        ----------
        NASA nodal feed/pressurization modeling:
        https://ntrs.nasa.gov/citations/20240003493
        """

        return self.specific_gas_constant_j_kg_k / (self.gamma - 1.0)

    @property
    def specific_heat_cp_j_kg_k(self) -> float:
        """Return ``cp = gamma*R/(gamma-1)``.

        References
        ----------
        NASA nodal feed/pressurization modeling:
        https://ntrs.nasa.gov/citations/20240003493
        """

        return self.gamma * self.specific_heat_cv_j_kg_k


@dataclass(frozen=True, slots=True)
class RegulatedFeedConfiguration:
    """Immutable geometry, control, thermal, and cutoff configuration.

    The regulator command is proportional to tank-pressure error and clipped
    to ``[0, 1]``.  Flow through the effective regulator area follows the
    isentropic compressible-orifice relation, including choking.  Liquid line
    loss is ``delta_p = resistance * mdot**2``.

    References
    ----------
    NASA-SP-8112: https://ntrs.nasa.gov/citations/19760015212
    NASA-SP-8080: https://ntrs.nasa.gov/citations/19740002611
    NASA mass-flow choking relation:
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

    References
    ----------
    NASA nodal feed/pressurization modeling:
    https://ntrs.nasa.gov/citations/20240003493
    NASA-SP-8112: https://ntrs.nasa.gov/citations/19760015212
    """

    propellant_mass_kg: float
    tank_pressurant_mass_kg: float
    tank_gas_temperature_k: float
    bottle_pressurant_mass_kg: float
    bottle_gas_temperature_k: float
    regulator_opening: float

    def __post_init__(self) -> None:
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
    """ODE right-hand side in the exact state order.

    References
    ----------
    NASA nodal feed/pressurization modeling:
    https://ntrs.nasa.gov/citations/20240003493
    """

    propellant_mass_kg_s: float
    tank_pressurant_mass_kg_s: float
    tank_gas_temperature_k_s: float
    bottle_pressurant_mass_kg_s: float
    bottle_gas_temperature_k_s: float
    regulator_opening_s: float

    def as_tuple(self) -> tuple[float, ...]:
        """Return derivatives in `RegulatedFeedState` field order.

        References
        ----------
        NASA nodal feed/pressurization modeling:
        https://ntrs.nasa.gov/citations/20240003493
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

    References
    ----------
    NASA nodal feed/pressurization modeling:
    https://ntrs.nasa.gov/citations/20240003493
    """

    injector_pressure_wrt_state: tuple[tuple[str, float], ...]
    injector_pressure_wrt_propellant_flow: float


@dataclass(frozen=True, slots=True)
class RegulatedFeedEvaluation:
    """Algebraic feed state and simultaneous ODE derivative.

    References
    ----------
    NASA-SP-8112: https://ntrs.nasa.gov/citations/19760015212
    NASA nodal feed/pressurization modeling:
    https://ntrs.nasa.gov/citations/20240003493
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


class FeedEventKind(str, Enum):
    """Protective cutoff surfaces owned by the feed system.

    References
    ----------
    NASA-SP-8112: https://ntrs.nasa.gov/citations/19760015212
    NASA-SP-8080: https://ntrs.nasa.gov/citations/19740002611
    """

    PROPELLANT_RESERVE = "propellant_reserve"
    INJECTOR_PRESSURE = "injector_pressure"
    BOTTLE_PRESSURE_MARGIN = "bottle_pressure_margin"
    TANK_OVERPRESSURE = "tank_overpressure"


@dataclass(frozen=True, slots=True)
class FeedSystemEventSurface:
    """Positive-safe signed root function for one feed-system limit.

    References
    ----------
    NASA-SP-8112: https://ntrs.nasa.gov/citations/19760015212
    Orekit EventDetector 13.1.5 API:
    https://www.orekit.org/static/apidocs/org/orekit/propagation/events/EventDetector.html
    """

    event_id: str
    kind: FeedEventKind
    stops_burn: bool = True


@dataclass(frozen=True, slots=True)
class FeedCouplingEvidence:
    """Fixed-point convergence evidence for feed/performance closure.

    References
    ----------
    NASA transient-modeling overview:
    https://ntrs.nasa.gov/citations/20040000363
    NASA nodal feed/pressurization modeling:
    https://ntrs.nasa.gov/citations/20240003493
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
    """Implicit derivatives through the feed-pressure/engine-flow loop.

    If engine demand is ``q = Q(p)`` and injector pressure is ``p = P(x,q)``,
    the closed derivative is

    ``dp/dx = P_x / (1 - P_q Q_p)``.

    The reported denominator makes proximity to an algebraic singularity
    auditable instead of hiding it inside a fixed-point iteration.

    References
    ----------
    NASA transient-modeling overview:
    https://ntrs.nasa.gov/citations/20040000363
    NASA nodal feed/pressurization modeling:
    https://ntrs.nasa.gov/citations/20240003493
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

    References
    ----------
    NASA transient-modeling overview:
    https://ntrs.nasa.gov/citations/20040000363
    NASA nodal feed/pressurization modeling:
    https://ntrs.nasa.gov/citations/20240003493
    """

    propulsion: PropulsionDynamicsEvaluation
    feed: RegulatedFeedEvaluation
    evidence: FeedCouplingEvidence
    partials: CoupledFeedPartials


@dataclass(frozen=True, slots=True)
class RegulatedFeedSystem:
    """Single-ullage-node regulated feed model for propagation.

    The tank energy equation includes inlet enthalpy, ullage boundary work and
    wall heat transfer.  The fixed-volume bottle includes outlet enthalpy and
    wall heat transfer.  Pressurant mass is conserved exactly between nodes.

    The model omits line inertance, acoustic/water-hammer modes, stratified
    ullage, phase change, dissolved gas, regulator hysteresis and real-fluid
    properties.  Those omissions are explicit because fast pressure events may
    require a higher-order validated nodal model.

    References
    ----------
    NASA-SP-8112: https://ntrs.nasa.gov/citations/19760015212
    NASA-SP-8080: https://ntrs.nasa.gov/citations/19740002611
    NASA nodal feed/pressurization modeling:
    https://ntrs.nasa.gov/citations/20240003493
    NASA NESC transient-pressure bulletin:
    https://ntrs.nasa.gov/citations/20220006583
    """

    configuration: RegulatedFeedConfiguration

    @property
    def model_hash(self) -> str:
        """Bind configuration, units, equations, and source identifiers.

        References
        ----------
        NASA-STD-7009B: https://standards.nasa.gov/standard/NASA/NASA-STD-7009
        """

        return _canonical_hash(self.manifest(include_hash=False))

    def manifest(self, *, include_hash: bool = True) -> dict[str, object]:
        """Return the deterministic integration manifest.

        References
        ----------
        NASA-STD-7009B: https://standards.nasa.gov/standard/NASA/NASA-STD-7009
        NASA nodal feed/pressurization modeling:
        https://ntrs.nasa.gov/citations/20240003493
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
        }
        if include_hash:
            payload["model_hash"] = _canonical_hash(payload)
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

        References
        ----------
        NASA-SP-8112: https://ntrs.nasa.gov/citations/19760015212
        NASA nodal feed/pressurization modeling:
        https://ntrs.nasa.gov/citations/20240003493
        """

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
        config = self.configuration
        maximum_propellant = (
            config.tank_internal_volume_m3 * config.propellant_density_kg_m3
        )
        if not isfinite(propellant_mass_kg) or not 0.0 < propellant_mass_kg < maximum_propellant:
            raise DomainError(
                "Propellant mass must leave positive ullage and remain below tank capacity."
            )
        return config.tank_internal_volume_m3 - (
            propellant_mass_kg / config.propellant_density_kg_m3
        )

    def _pressures(self, state: RegulatedFeedState) -> tuple[float, float, float]:
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
        return ullage, tank_pressure, bottle_pressure

    def _regulator_flow(
        self,
        *,
        upstream_pressure_pa: float,
        upstream_temperature_k: float,
        downstream_pressure_pa: float,
        opening: float,
    ) -> tuple[float, bool]:
        config = self.configuration
        gas = config.gas
        if upstream_pressure_pa <= downstream_pressure_pa or opening <= 0.0:
            return 0.0, False
        area = config.regulator_maximum_area_m2 * opening
        ratio = downstream_pressure_pa / upstream_pressure_pa
        critical_ratio = (2.0 / (gas.gamma + 1.0)) ** (
            gas.gamma / (gas.gamma - 1.0)
        )
        common = (
            config.regulator_discharge_coefficient
            * area
            * upstream_pressure_pa
            / sqrt(gas.specific_gas_constant_j_kg_k * upstream_temperature_k)
        )
        if ratio <= critical_ratio:
            factor = sqrt(gas.gamma) * (2.0 / (gas.gamma + 1.0)) ** (
                (gas.gamma + 1.0) / (2.0 * (gas.gamma - 1.0))
            )
            return common * factor, True
        radicand = (
            2.0
            * gas.gamma
            / (gas.gamma - 1.0)
            * (ratio ** (2.0 / gas.gamma) - ratio ** ((gas.gamma + 1.0) / gas.gamma))
        )
        return common * sqrt(max(0.0, radicand)), False

    def evaluate(
        self,
        state: RegulatedFeedState,
        propellant_mass_flow_kg_s: float,
    ) -> RegulatedFeedEvaluation:
        """Return algebraic pressures and the six-state ODE right-hand side.

        Positive ``propellant_mass_flow_kg_s`` drains the liquid tank.  The gas
        node energy balances are

        ``d(m_t cv T_t)/dt = mdot_r cp T_b + Q_t - p_t dV_u/dt``

        and

        ``d(m_b cv T_b)/dt = -mdot_r cp T_b + Q_b``.

        References
        ----------
        NASA-SP-8112: https://ntrs.nasa.gov/citations/19760015212
        NASA nodal feed/pressurization modeling:
        https://ntrs.nasa.gov/citations/20240003493
        NASA mass-flow choking relation:
        https://www.grc.nasa.gov/www/k-12/BGP/mflchk.html
        """

        _nonnegative(propellant_mass_flow_kg_s, "Propellant mass-flow demand")
        config = self.configuration
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
        cp = gas.specific_heat_cp_j_kg_k
        tank_temperature_rate = (
            regulator_flow
            * (cp * state.bottle_gas_temperature_k - cv * state.tank_gas_temperature_k)
            + tank_heat
            - tank_pressure * ullage_rate
        ) / (state.tank_pressurant_mass_kg * cv)
        bottle_temperature_rate = (
            bottle_heat
            - regulator_flow
            * (cp - cv)
            * state.bottle_gas_temperature_k
        ) / (state.bottle_pressurant_mass_kg * cv)
        pressure_drop = (
            config.feed_line_resistance_pa_s2_kg2 * propellant_mass_flow_kg_s**2
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
                -2.0
                * config.feed_line_resistance_pa_s2_kg2
                * propellant_mass_flow_kg_s
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
        )

    def event_surfaces(self) -> tuple[FeedSystemEventSurface, ...]:
        """Return all configured positive-safe feed cutoff roots.

        References
        ----------
        NASA-SP-8112: https://ntrs.nasa.gov/citations/19760015212
        NASA-SP-8080: https://ntrs.nasa.gov/citations/19740002611
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
        """Evaluate one positive-safe cutoff root.

        References
        ----------
        NASA-SP-8112: https://ntrs.nasa.gov/citations/19760015212
        Orekit EventDetector 13.1.5 API:
        https://www.orekit.org/static/apidocs/org/orekit/propagation/events/EventDetector.html
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
    relative_tolerance: float = 1e-10,
    relaxation: float = 0.5,
    maximum_iterations: int = 100,
) -> CoupledFeedPropulsionEvaluation:
    """Close engine flow demand against feed-line-dependent inlet pressure.

    The bridge's performance surface produces tank-flow demand from injector
    pressure, while the feed line produces injector pressure from that same
    flow.  A relaxed fixed-point iteration resolves this algebraic loop and
    returns explicit residual evidence.  Non-convergence fails closed.

    References
    ----------
    NASA transient-modeling overview:
    https://ntrs.nasa.gov/citations/20040000363
    NASA nodal feed/pressurization modeling:
    https://ntrs.nasa.gov/citations/20240003493
    NASA computer-controlled pressurization study:
    https://ntrs.nasa.gov/citations/20050206420
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

