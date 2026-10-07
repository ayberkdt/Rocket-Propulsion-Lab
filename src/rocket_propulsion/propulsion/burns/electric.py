"""Power-limited electric propulsion for numerical orbit propagation.

An electric thruster's operating point is set by the electrical power the
spacecraft can deliver, not by a time schedule.  The astrodynamics framework
owns that power (solar distance, array degradation, eclipse, bus load); this
module receives the resulting *available* power and returns thrust, propellant
drain, consumed power, analytic partials, and the power-switching roots a
propagator must register.

Two throttle laws are provided:

``DiscreteThrottlePath``
    An ordered subset of a hash-verified tabulated throttle table (e.g. an
    ion/Hall thruster qualification table).  The highest level whose required
    input power fits the available power is selected, with an optional
    up-switch hysteresis band so that power noise near a threshold does not
    chatter between levels.  Within a level thrust and flow are constant, so
    the acceleration does not depend on available power; level changes are
    explicit events.
``FixedEfficiencyThrottle``
    A continuous law with constant total efficiency ``eta`` and specific
    impulse.  From jet power ``P_jet = T g0 Isp / 2 = eta P``::

        T = 2 eta P / (g0 Isp),     mdot = T / (g0 Isp)

    between a minimum (start) and a maximum (saturation) input power.  This
    gives a non-zero, analytic ``dT/dP`` for low-thrust trajectory
    optimization.

The discrete level index is a propagator-owned mode state, as in an
event-driven maneuver model: the consumer initializes it with
:meth:`PowerLimitedPropulsionModel.select_level`, registers
:meth:`~PowerLimitedPropulsionModel.switching_surfaces` for the current level,
and after locating a root switches to that surface's ``to_level``.  The target
is fixed by the crossing direction, so a root located a round-off on either
side of the threshold cannot leave the mode unchanged and re-trigger the same
event.  ``evaluate`` never changes the level by itself, so Runge-Kutta
stages that probe slightly past a switching root stay on a smooth branch; the
evaluation reports any ``power_deficit_w`` instead.
"""

from __future__ import annotations

import hashlib
import json
from collections.abc import Mapping
from dataclasses import dataclass, field
from itertools import pairwise
from math import isfinite

from rocket_propulsion.core.errors import DomainError, InputError
from rocket_propulsion.propulsion.performance import STANDARD_GRAVITY_M_S2

from .propagator import (
    BurnEventKind,
    BurnEventSurface,
    EstimableParameter,
    PropellantMassClosure,
    _normalized_with_norm,
    _validated_states,
)

ELECTRIC_THROTTLE_TABLE_SCHEMA = "rocket_propulsion_electric_throttle_table_v1"
POWER_LIMITED_PROPULSION_SCHEMA = "rocket_propulsion_power_limited_propulsion_v1"
ELECTRIC_PROPULSION_REFERENCE_IDS = (
    "NASA-20090004685",
    "OREKIT-13.1.5-PROPULSION-MODEL",
)
OFF_LEVEL = -1
_POWER_TOLERANCE_W = 1e-9


def _canonical_hash(payload: object) -> str:
    encoded = json.dumps(
        payload,
        allow_nan=False,
        ensure_ascii=False,
        separators=(",", ":"),
        sort_keys=True,
    ).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()


def _positive(value: float, label: str) -> None:
    if not isfinite(value) or value <= 0.0:
        raise DomainError(f"{label} must be finite and greater than zero.")


def _nonnegative(value: float, label: str) -> None:
    if not isfinite(value) or value < 0.0:
        raise DomainError(f"{label} must be finite and non-negative.")


def _sha256_hex(value: str, label: str) -> None:
    if len(value) != 64 or any(character not in "0123456789abcdef" for character in value):
        raise DomainError(f"{label} must be a lower-case SHA-256 hex digest.")


@dataclass(frozen=True, slots=True)
class ElectricThrottleLevel:
    """One qualified electric-thruster operating point.

    ``input_power_w`` is the power drawn from the spacecraft bus by the
    thruster string (power-processing unit input), which is the quantity the
    power model must supply.  The derived total efficiency
    ``T^2 / (2 mdot P)`` is the ratio of jet power to input power and must not
    exceed one; a table that violates this is rejected as inconsistent.

    References
    ----------
    NASA-20090004685: https://ntrs.nasa.gov/citations/20090004685
    """

    level_id: str
    input_power_w: float
    thrust_n: float
    mass_flow_kg_s: float

    def __post_init__(self) -> None:
        if not self.level_id.strip():
            raise DomainError("Throttle level identifier must not be empty.")
        _positive(self.input_power_w, "Throttle level input power")
        _positive(self.thrust_n, "Throttle level thrust")
        _positive(self.mass_flow_kg_s, "Throttle level mass flow")
        if self.total_efficiency > 1.0 + 1e-12:
            raise DomainError(
                f"Throttle level {self.level_id!r} implies jet power above input power."
            )

    @property
    def specific_impulse_s(self) -> float:
        """Return ``T / (g0 mdot)`` in seconds."""

        return self.thrust_n / (STANDARD_GRAVITY_M_S2 * self.mass_flow_kg_s)

    @property
    def total_efficiency(self) -> float:
        """Return jet power over input power, ``T^2 / (2 mdot P)``."""

        return self.thrust_n**2 / (2.0 * self.mass_flow_kg_s * self.input_power_w)


@dataclass(frozen=True, slots=True)
class ElectricThrottleTable:
    """Hash-verified tabulated throttle table with source provenance.

    The table is evidence: values come from a qualification or acceptance
    test campaign identified by ``source_id`` and ``source_sha256``.  No level
    is interpolated or extrapolated.

    References
    ----------
    NASA-20090004685: https://ntrs.nasa.gov/citations/20090004685
    CCSDS-502.0-B-3: https://ccsds.org/Pubs/502x0b3e1.pdf
    """

    levels: tuple[ElectricThrottleLevel, ...]
    source_id: str
    source_sha256: str
    reference_ids: tuple[str, ...] = ("NASA-20090004685",)
    schema: str = ELECTRIC_THROTTLE_TABLE_SCHEMA

    def __post_init__(self) -> None:
        if self.schema != ELECTRIC_THROTTLE_TABLE_SCHEMA:
            raise DomainError("Unsupported electric throttle-table schema.")
        if not self.levels:
            raise DomainError("An electric throttle table requires at least one level.")
        identifiers = [level.level_id for level in self.levels]
        if len(identifiers) != len(set(identifiers)):
            raise DomainError("Throttle level identifiers must be unique.")
        if not self.source_id.strip():
            raise DomainError("Throttle-table source identifier must not be empty.")
        _sha256_hex(self.source_sha256, "Throttle-table source hash")
        if not self.reference_ids or any(not item.strip() for item in self.reference_ids):
            raise DomainError("Throttle-table references must not be empty.")
        if len(self.reference_ids) != len(set(self.reference_ids)):
            raise DomainError("Throttle-table references must be unique.")

    def level(self, level_id: str) -> ElectricThrottleLevel:
        """Return one level by identifier or fail closed."""

        for level in self.levels:
            if level.level_id == level_id:
                return level
        raise DomainError(f"Unknown throttle level {level_id!r}.")

    @property
    def table_hash(self) -> str:
        """Return content identity of levels and provenance."""

        return _canonical_hash(self.to_dict(include_hash=False))

    def to_dict(self, *, include_hash: bool = True) -> dict[str, object]:
        """Return deterministic JSON-compatible table content."""

        payload: dict[str, object] = {
            "schema": self.schema,
            "source_id": self.source_id,
            "source_sha256": self.source_sha256,
            "reference_ids": list(self.reference_ids),
            "levels": [
                {
                    "level_id": level.level_id,
                    "input_power_w": level.input_power_w,
                    "thrust_n": level.thrust_n,
                    "mass_flow_kg_s": level.mass_flow_kg_s,
                }
                for level in self.levels
            ],
        }
        if include_hash:
            payload["table_hash"] = _canonical_hash(payload)
        return payload

    def to_json(self, *, indent: int | None = 2) -> str:
        """Serialize the complete hash-bound table."""

        return json.dumps(
            self.to_dict(), allow_nan=False, ensure_ascii=False, indent=indent, sort_keys=True
        )

    @classmethod
    def from_dict(cls, payload: Mapping[str, object]) -> ElectricThrottleTable:
        """Reopen a table and reject malformed or hash-altered content."""

        if not isinstance(payload, Mapping):
            raise InputError("Throttle-table artifact root must be an object.")
        try:
            levels_payload = payload["levels"]
            if not isinstance(levels_payload, list):
                raise TypeError
            table = cls(
                schema=str(payload["schema"]),
                source_id=str(payload["source_id"]),
                source_sha256=str(payload["source_sha256"]),
                reference_ids=tuple(str(item) for item in payload["reference_ids"]),
                levels=tuple(
                    ElectricThrottleLevel(
                        level_id=str(item["level_id"]),
                        input_power_w=float(item["input_power_w"]),
                        thrust_n=float(item["thrust_n"]),
                        mass_flow_kg_s=float(item["mass_flow_kg_s"]),
                    )
                    for item in levels_payload
                ),
            )
        except (KeyError, TypeError, ValueError, DomainError) as error:
            raise InputError("Malformed electric throttle-table artifact.") from error
        if payload.get("table_hash") != table.table_hash:
            raise InputError("Electric throttle-table hash does not match its contents.")
        return table

    @classmethod
    def from_json(cls, text: str) -> ElectricThrottleTable:
        """Deserialize and integrity-check a JSON throttle table."""

        try:
            payload = json.loads(text)
        except json.JSONDecodeError as error:
            raise InputError("Electric throttle-table artifact is not valid JSON.") from error
        if not isinstance(payload, dict):
            raise InputError("Electric throttle-table artifact root must be an object.")
        return cls.from_dict(payload)


@dataclass(frozen=True, slots=True)
class LevelOperatingPoint:
    """Thrust, flow and power of the active level at one available power.

    References
    ----------
    NASA-20090004685: https://ntrs.nasa.gov/citations/20090004685
    """

    level_index: int
    level_id: str
    thrust_n: float
    mass_flow_kg_s: float
    input_power_w: float
    required_power_w: float
    thrust_wrt_available_power: float
    mass_flow_wrt_available_power: float


@dataclass(frozen=True, slots=True)
class PowerSwitchSurface:
    """Signed root on available power that changes the throttle level.

    ``value = available_power_w - threshold_w``.  A ``decreasing`` surface
    fires when the value crosses zero from above (power fell below what the
    current level needs); an ``increasing`` surface fires when it crosses from
    below (power now supports the next level plus hysteresis).  After the
    root the consumer switches to ``to_level``.  An ``either`` surface has
    ``to_level == from_level``: it only marks a derivative kink, fires in both
    directions, and asks for an integrator restart.

    References
    ----------
    Orekit EventDetector 13.1.5 API:
    https://www.orekit.org/static/apidocs/org/orekit/propagation/events/EventDetector.html
    """

    event_id: str
    threshold_w: float
    crossing: str
    from_level: int
    to_level: int

    def __post_init__(self) -> None:
        if not self.event_id.strip():
            raise DomainError("Power-switch event identifier must not be empty.")
        if self.from_level < OFF_LEVEL or self.to_level < OFF_LEVEL:
            raise DomainError("Power-switch levels must be -1 (off) or a level index.")
        if not isfinite(self.threshold_w):
            raise DomainError("Power-switch threshold must be finite.")
        if self.crossing not in {"increasing", "decreasing", "either"}:
            raise DomainError("Power-switch crossing must be increasing, decreasing or either.")
        if self.crossing == "either" and self.to_level != self.from_level:
            raise DomainError("Only a mode-preserving kink surface may fire in either direction.")

    def value(self, available_power_w: float) -> float:
        """Evaluate ``available_power_w - threshold_w``."""

        if not isfinite(available_power_w):
            raise DomainError("Available power must be finite.")
        return available_power_w - self.threshold_w


@dataclass(frozen=True, slots=True)
class DiscreteThrottlePath:
    """Ordered throttle levels with strictly increasing required power.

    ``up_switch_hysteresis_w`` raises only the up-switch threshold: from level
    ``k`` the model moves to ``k+1`` when power reaches
    ``P_{k+1} + hysteresis``, but drops back only when power falls below
    ``P_k``.  Starting from off, the lowest level likewise needs
    ``P_0 + hysteresis``.

    References
    ----------
    NASA-20090004685: https://ntrs.nasa.gov/citations/20090004685
    """

    table: ElectricThrottleTable
    level_ids: tuple[str, ...]
    up_switch_hysteresis_w: float = 0.0
    path_id: str = "discrete-throttle-path"
    _levels: tuple[ElectricThrottleLevel, ...] = field(
        init=False, repr=False, compare=False
    )

    def __post_init__(self) -> None:
        if not self.level_ids:
            raise DomainError("A throttle path requires at least one level.")
        if len(self.level_ids) != len(set(self.level_ids)):
            raise DomainError("A throttle path may use each level only once.")
        if not self.path_id.strip():
            raise DomainError("Throttle path identifier must not be empty.")
        _nonnegative(self.up_switch_hysteresis_w, "Up-switch hysteresis")
        levels = tuple(self.table.level(level_id) for level_id in self.level_ids)
        if any(
            upper.input_power_w <= lower.input_power_w for lower, upper in pairwise(levels)
        ):
            raise DomainError("Throttle-path required power must be strictly increasing.")
        # Resolved once: operating_point runs at every propagator stage.
        object.__setattr__(self, "_levels", levels)

    @classmethod
    def from_table(
        cls,
        table: ElectricThrottleTable,
        *,
        objective: str = "thrust",
        up_switch_hysteresis_w: float = 0.0,
    ) -> DiscreteThrottlePath:
        """Build the path that maximizes ``thrust`` or ``specific_impulse``.

        Levels are visited in increasing required power; a level enters the
        path only if it strictly improves the objective over every cheaper
        level already on the path, so more power never selects a worse level.
        At equal power the better level wins.

        References
        ----------
        NASA-20090004685: https://ntrs.nasa.gov/citations/20090004685
        """

        if objective not in {"thrust", "specific_impulse"}:
            raise DomainError("Throttle-path objective must be thrust or specific_impulse.")

        def score(level: ElectricThrottleLevel) -> float:
            return level.thrust_n if objective == "thrust" else level.specific_impulse_s

        ordered = sorted(table.levels, key=lambda level: (level.input_power_w, -score(level)))
        selected: list[ElectricThrottleLevel] = []
        for level in ordered:
            if selected and level.input_power_w == selected[-1].input_power_w:
                continue
            if not selected or score(level) > score(selected[-1]):
                selected.append(level)
        return cls(
            table,
            tuple(level.level_id for level in selected),
            up_switch_hysteresis_w,
            f"max-{objective}",
        )

    @property
    def levels(self) -> tuple[ElectricThrottleLevel, ...]:
        """Return the path levels in increasing required power."""

        return self._levels

    def select_level(self, available_power_w: float, current_level: int | None) -> int:
        """Return the level index for ``available_power_w`` (``-1`` = off)."""

        levels = self.levels
        feasible = OFF_LEVEL
        for index, level in enumerate(levels):
            if level.input_power_w <= available_power_w + _POWER_TOLERANCE_W:
                feasible = index
        if current_level is None:
            return feasible
        if current_level < OFF_LEVEL or current_level >= len(levels):
            raise DomainError("Current throttle level is outside the path.")
        if current_level != OFF_LEVEL and (
            available_power_w < levels[current_level].input_power_w - _POWER_TOLERANCE_W
        ):
            return feasible
        level = current_level
        while (
            level + 1 < len(levels)
            and available_power_w + _POWER_TOLERANCE_W
            >= levels[level + 1].input_power_w + self.up_switch_hysteresis_w
        ):
            level += 1
        return level

    def switching_surfaces(self, current_level: int) -> tuple[PowerSwitchSurface, ...]:
        """Return the down- and up-switch roots for the current level."""

        levels = self.levels
        if current_level < OFF_LEVEL or current_level >= len(levels):
            raise DomainError("Current throttle level is outside the path.")
        surfaces: list[PowerSwitchSurface] = []
        if current_level != OFF_LEVEL:
            surfaces.append(
                PowerSwitchSurface(
                    f"power:{levels[current_level].level_id}:down",
                    levels[current_level].input_power_w,
                    "decreasing",
                    current_level,
                    current_level - 1,
                )
            )
        if current_level + 1 < len(levels):
            following = levels[current_level + 1]
            surfaces.append(
                PowerSwitchSurface(
                    f"power:{following.level_id}:up",
                    following.input_power_w + self.up_switch_hysteresis_w,
                    "increasing",
                    current_level,
                    current_level + 1,
                )
            )
        return tuple(surfaces)

    def operating_point(self, level_index: int, available_power_w: float) -> LevelOperatingPoint:
        """Return the constant operating point of ``level_index``."""

        level = self.levels[level_index]
        return LevelOperatingPoint(
            level_index=level_index,
            level_id=level.level_id,
            thrust_n=level.thrust_n,
            mass_flow_kg_s=level.mass_flow_kg_s,
            input_power_w=level.input_power_w,
            required_power_w=level.input_power_w,
            thrust_wrt_available_power=0.0,
            mass_flow_wrt_available_power=0.0,
        )

    @property
    def level_count(self) -> int:
        return len(self.level_ids)

    def manifest(self) -> dict[str, object]:
        return {
            "law": "discrete_throttle_path",
            "path_id": self.path_id,
            "table_hash": self.table.table_hash,
            "level_ids": list(self.level_ids),
            "up_switch_hysteresis_w": self.up_switch_hysteresis_w,
        }


@dataclass(frozen=True, slots=True)
class FixedEfficiencyThrottle:
    """Continuous constant-efficiency, constant-Isp power-to-thrust law.

    Between ``minimum_power_w`` and ``maximum_power_w`` the thruster uses all
    available power: ``T = 2 eta P / (g0 Isp)`` and ``mdot = T / (g0 Isp)``.
    Above the maximum it saturates (``dT/dP = 0``); below the minimum it is
    off.  Level ``0`` is on and ``-1`` is off.  The start threshold is
    ``minimum_power_w + start_hysteresis_w``.

    References
    ----------
    NASA-20090004685: https://ntrs.nasa.gov/citations/20090004685
    """

    total_efficiency: float
    specific_impulse_s: float
    minimum_power_w: float
    maximum_power_w: float
    start_hysteresis_w: float = 0.0
    law_id: str = "fixed-efficiency"

    def __post_init__(self) -> None:
        if not isfinite(self.total_efficiency) or not 0.0 < self.total_efficiency <= 1.0:
            raise DomainError("Total efficiency must be in (0, 1].")
        _positive(self.specific_impulse_s, "Specific impulse")
        _positive(self.minimum_power_w, "Minimum input power")
        _positive(self.maximum_power_w, "Maximum input power")
        if self.maximum_power_w <= self.minimum_power_w:
            raise DomainError("Maximum input power must exceed the minimum.")
        _nonnegative(self.start_hysteresis_w, "Start hysteresis")
        if not self.law_id.strip():
            raise DomainError("Throttle law identifier must not be empty.")

    @property
    def level_count(self) -> int:
        return 1

    @property
    def thrust_per_watt_n_w(self) -> float:
        """Return ``2 eta / (g0 Isp)``."""

        return 2.0 * self.total_efficiency / (STANDARD_GRAVITY_M_S2 * self.specific_impulse_s)

    def select_level(self, available_power_w: float, current_level: int | None) -> int:
        """Return ``0`` when the thruster runs and ``-1`` when it is off."""

        if current_level is not None and current_level not in {OFF_LEVEL, 0}:
            raise DomainError("Current throttle level is outside the continuous law.")
        threshold = self.minimum_power_w
        if current_level == OFF_LEVEL:
            threshold += self.start_hysteresis_w
        return 0 if available_power_w + _POWER_TOLERANCE_W >= threshold else OFF_LEVEL

    def switching_surfaces(self, current_level: int) -> tuple[PowerSwitchSurface, ...]:
        """Return the shutdown root when on, or the start root when off."""

        if current_level == 0:
            return (
                PowerSwitchSurface(
                    f"power:{self.law_id}:down",
                    self.minimum_power_w,
                    "decreasing",
                    0,
                    OFF_LEVEL,
                ),
            )
        if current_level == OFF_LEVEL:
            return (
                PowerSwitchSurface(
                    f"power:{self.law_id}:up",
                    self.minimum_power_w + self.start_hysteresis_w,
                    "increasing",
                    OFF_LEVEL,
                    0,
                ),
            )
        raise DomainError("Current throttle level is outside the continuous law.")

    def kink_surfaces(self) -> tuple[PowerSwitchSurface, ...]:
        """Return the non-mode-changing saturation root where ``dT/dP`` jumps."""

        return (
            PowerSwitchSurface(
                f"power:{self.law_id}:saturation", self.maximum_power_w, "either", 0, 0
            ),
        )

    def operating_point(self, level_index: int, available_power_w: float) -> LevelOperatingPoint:
        """Return thrust and flow using ``min(P, P_max)`` and their power slopes."""

        if level_index != 0:
            raise DomainError("Continuous law operating point requires level 0.")
        saturated = available_power_w >= self.maximum_power_w
        power = self.maximum_power_w if saturated else max(available_power_w, 0.0)
        gain = self.thrust_per_watt_n_w
        exhaust_velocity = STANDARD_GRAVITY_M_S2 * self.specific_impulse_s
        thrust = gain * power
        thrust_slope = 0.0 if saturated else gain
        return LevelOperatingPoint(
            level_index=0,
            level_id=self.law_id,
            thrust_n=thrust,
            mass_flow_kg_s=thrust / exhaust_velocity,
            input_power_w=power,
            required_power_w=self.minimum_power_w,
            thrust_wrt_available_power=thrust_slope,
            mass_flow_wrt_available_power=thrust_slope / exhaust_velocity,
        )

    def manifest(self) -> dict[str, object]:
        return {
            "law": "fixed_efficiency",
            "law_id": self.law_id,
            "total_efficiency": self.total_efficiency,
            "specific_impulse_s": self.specific_impulse_s,
            "minimum_power_w": self.minimum_power_w,
            "maximum_power_w": self.maximum_power_w,
            "start_hysteresis_w": self.start_hysteresis_w,
        }


@dataclass(frozen=True, slots=True)
class ElectricPropulsionPartials:
    """Analytic local derivatives of one power-limited evaluation.

    ``acceleration_wrt_direction`` is ``da_i/dd_j`` for the supplied, not
    necessarily unit, direction ``d``.  Available-power partials let a
    consumer chain the thrust law with its own power model
    (``dP/dr`` from solar distance, eclipse factors); they are zero inside a
    discrete level and above continuous saturation.

    References
    ----------
    Orekit PropulsionModel 13.1.5 API:
    https://www.orekit.org/static/apidocs/org/orekit/forces/maneuvers/propulsion/PropulsionModel.html
    """

    acceleration_wrt_mass: tuple[float, float, float]
    acceleration_wrt_direction: tuple[
        tuple[float, float, float],
        tuple[float, float, float],
        tuple[float, float, float],
    ]
    acceleration_wrt_available_power: tuple[float, float, float]
    mass_rate_wrt_available_power: float
    acceleration_wrt_parameters: tuple[tuple[str, tuple[float, float, float]], ...]
    mass_rate_wrt_parameters: tuple[tuple[str, float], ...]


@dataclass(frozen=True, slots=True)
class ElectricPropulsionEvaluation:
    """Force, drain, power use and partials at one propagator call.

    ``power_deficit_w`` is positive when the commanded level needs more power
    than is available, which happens legitimately in Runge-Kutta stages just
    past a down-switch root.  A deficit after an accepted step means the
    consumer missed a switching event.

    References
    ----------
    NASA-20090004685: https://ntrs.nasa.gov/citations/20090004685
    """

    active: bool
    level_index: int
    level_id: str | None
    thrust_n: float
    force_vector_n: tuple[float, float, float]
    acceleration_m_s2: tuple[float, float, float]
    mass_derivative_kg_s: float
    additional_state_derivatives: tuple[tuple[str, float], ...]
    available_power_w: float
    consumed_power_w: float
    power_deficit_w: float
    specific_impulse_s: float | None
    parameters: tuple[tuple[str, float], ...]
    inhibited_by: tuple[str, ...]
    partials: ElectricPropulsionPartials


_ZERO3 = (0.0, 0.0, 0.0)


@dataclass(frozen=True, slots=True)
class PowerLimitedPropulsionModel:
    """Immutable power-limited electric propulsion model for propagation.

    ``power_margin_w`` is subtracted from the available power before level
    selection (for example a housekeeping reserve).  ``duty_cycle`` scales
    thrust, flow and consumed power to represent an averaged on/off firing
    pattern; it is not a substitute for modelling individual pulses.

    References
    ----------
    NASA-20090004685: https://ntrs.nasa.gov/citations/20090004685
    Orekit PropulsionModel 13.1.5 API:
    https://www.orekit.org/static/apidocs/org/orekit/forces/maneuvers/propulsion/PropulsionModel.html
    """

    throttle_law: DiscreteThrottlePath | FixedEfficiencyThrottle
    propellant_state_name: str = "xenon_mass_kg"
    propellant_reserve_kg: float | None = None
    protected_dry_mass_kg: float | None = None
    power_margin_w: float = 0.0
    duty_cycle: float = 1.0
    parameters: tuple[EstimableParameter, ...] = ()
    model_id: str = "power-limited-electric-propulsion"
    schema: str = POWER_LIMITED_PROPULSION_SCHEMA

    def __post_init__(self) -> None:
        if self.schema != POWER_LIMITED_PROPULSION_SCHEMA:
            raise DomainError("Unsupported power-limited propulsion schema.")
        if not self.model_id.strip() or not self.propellant_state_name.strip():
            raise DomainError("Model and propellant-state names must not be empty.")
        if self.propellant_reserve_kg is not None:
            _nonnegative(self.propellant_reserve_kg, "Propellant reserve")
        if self.protected_dry_mass_kg is not None:
            _positive(self.protected_dry_mass_kg, "Protected dry mass")
        _nonnegative(self.power_margin_w, "Power margin")
        if not isfinite(self.duty_cycle) or not 0.0 < self.duty_cycle <= 1.0:
            raise DomainError("Duty cycle must be in (0, 1].")
        if not self.parameters:
            object.__setattr__(
                self,
                "parameters",
                (
                    EstimableParameter("thrust_scale", 1.0, 1e-3, 0.01, 2.0, "1"),
                    EstimableParameter("mass_flow_scale", 1.0, 1e-3, 0.01, 2.0, "1"),
                ),
            )
        names = [parameter.name for parameter in self.parameters]
        if sorted(names) != ["mass_flow_scale", "thrust_scale"]:
            raise DomainError(
                "Power-limited model requires exactly thrust_scale and mass_flow_scale."
            )

    def effective_power_w(self, available_power_w: float) -> float:
        """Return available power minus the declared margin."""

        if not isfinite(available_power_w):
            raise DomainError("Available power must be finite.")
        return available_power_w - self.power_margin_w

    def select_level(self, available_power_w: float, current_level: int | None = None) -> int:
        """Initialize (``current_level=None``) or update the throttle mode.

        References
        ----------
        NASA-20090004685: https://ntrs.nasa.gov/citations/20090004685
        """

        return self.throttle_law.select_level(
            self.effective_power_w(available_power_w), current_level
        )

    def switching_surfaces(self, current_level: int) -> tuple[PowerSwitchSurface, ...]:
        """Return mode-changing roots expressed in *available* power.

        References
        ----------
        Orekit EventDetector 13.1.5 API:
        https://www.orekit.org/static/apidocs/org/orekit/propagation/events/EventDetector.html
        """

        surfaces = self.throttle_law.switching_surfaces(current_level)
        if isinstance(self.throttle_law, FixedEfficiencyThrottle) and current_level == 0:
            surfaces = surfaces + self.throttle_law.kink_surfaces()
        return tuple(
            PowerSwitchSurface(
                surface.event_id,
                surface.threshold_w + self.power_margin_w,
                surface.crossing,
                surface.from_level,
                surface.to_level,
            )
            for surface in surfaces
        )

    def inventory_surfaces(self) -> tuple[BurnEventSurface, ...]:
        """Return stopping dry-mass and propellant-reserve roots.

        References
        ----------
        Orekit EventDetector 13.1.5 API:
        https://www.orekit.org/static/apidocs/org/orekit/propagation/events/EventDetector.html
        """

        surfaces: list[BurnEventSurface] = []
        if self.protected_dry_mass_kg is not None:
            surfaces.append(
                BurnEventSurface(
                    "electric:dry-mass",
                    BurnEventKind.DRY_MASS,
                    True,
                    threshold=self.protected_dry_mass_kg,
                )
            )
        if self.propellant_reserve_kg is not None:
            surfaces.append(
                BurnEventSurface(
                    f"electric:reserve:{self.propellant_state_name}",
                    BurnEventKind.TANK_RESERVE,
                    True,
                    state_name=self.propellant_state_name,
                    threshold=self.propellant_reserve_kg,
                )
            )
        return tuple(surfaces)

    def resolved_parameters(
        self, overrides: Mapping[str, float] | None = None
    ) -> dict[str, float]:
        """Resolve named overrides against the bounded driver definitions."""

        supplied = dict(overrides or {})
        known = {parameter.name for parameter in self.parameters}
        unknown = set(supplied) - known
        if unknown:
            raise DomainError("Unknown propulsion parameters: " + ", ".join(sorted(unknown)))
        return {
            parameter.name: parameter.validate(
                supplied.get(parameter.name, parameter.reference_value)
            )
            for parameter in self.parameters
        }

    def manifest(self, *, include_hash: bool = True) -> dict[str, object]:
        """Return a deterministic, integration-reviewable model manifest."""

        payload: dict[str, object] = {
            "schema": self.schema,
            "model_id": self.model_id,
            "throttle_law": self.throttle_law.manifest(),
            "propellant_state_name": self.propellant_state_name,
            "propellant_reserve_kg": self.propellant_reserve_kg,
            "protected_dry_mass_kg": self.protected_dry_mass_kg,
            "power_margin_w": self.power_margin_w,
            "duty_cycle": self.duty_cycle,
            "force_semantics": "thrust_along_resolved_direction",
            "mass_semantics": "negative_propellant_drain",
            "power_semantics": "available_bus_power_minus_margin_selects_level",
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
            "reference_ids": list(ELECTRIC_PROPULSION_REFERENCE_IDS),
        }
        if include_hash:
            payload["model_hash"] = _canonical_hash(payload)
        return payload

    @property
    def model_hash(self) -> str:
        """Content identity of the throttle law, limits and semantics."""

        return _canonical_hash(self.manifest(include_hash=False))

    def mass_closure(
        self,
        *,
        vehicle_mass_kg: float,
        additional_states: Mapping[str, float],
        non_propellant_mass_kg: float,
    ) -> PropellantMassClosure:
        """Compare vehicle mass with the tracked propellant state.

        References
        ----------
        NASA-SP-8112: https://ntrs.nasa.gov/citations/19760015212
        """

        _positive(vehicle_mass_kg, "Vehicle mass")
        _nonnegative(non_propellant_mass_kg, "Non-propellant mass")
        states = _validated_states(additional_states)
        if self.propellant_state_name not in states:
            raise DomainError(f"Missing propellant state {self.propellant_state_name!r}.")
        tracked = states[self.propellant_state_name]
        residual = vehicle_mass_kg - non_propellant_mass_kg - tracked
        return PropellantMassClosure(
            vehicle_mass_kg=vehicle_mass_kg,
            non_propellant_mass_kg=non_propellant_mass_kg,
            tracked_propellant_kg=tracked,
            residual_kg=residual,
            relative_residual=residual / vehicle_mass_kg,
        )

    def evaluate(
        self,
        *,
        level_index: int,
        available_power_w: float,
        vehicle_mass_kg: float,
        direction: tuple[float, float, float],
        additional_states: Mapping[str, float],
        parameter_overrides: Mapping[str, float] | None = None,
        apply_inventory_limits: bool = True,
    ) -> ElectricPropulsionEvaluation:
        """Evaluate thrust, drain and partials for the commanded mode.

        ``level_index`` is the propagator-owned mode (``-1`` = off).  The
        model does not switch it; see :meth:`select_level`.  As in
        :meth:`PropagatorPropulsionBridge.evaluate`, pass
        ``apply_inventory_limits=False`` when dry-mass and reserve roots are
        located by the consumer.

        References
        ----------
        NASA-20090004685: https://ntrs.nasa.gov/citations/20090004685
        Orekit PropulsionModel 13.1.5 API:
        https://www.orekit.org/static/apidocs/org/orekit/forces/maneuvers/propulsion/PropulsionModel.html
        """

        if level_index < OFF_LEVEL or level_index >= self.throttle_law.level_count:
            raise DomainError("Throttle level index is outside the throttle law.")
        _positive(vehicle_mass_kg, "Vehicle mass")
        effective_power = self.effective_power_w(available_power_w)
        unit, norm = _normalized_with_norm(direction)
        states = _validated_states(additional_states)
        if self.propellant_state_name not in states:
            raise DomainError(f"Missing propellant state {self.propellant_state_name!r}.")
        parameters = self.resolved_parameters(parameter_overrides)
        thrust_scale = parameters["thrust_scale"]
        flow_scale = parameters["mass_flow_scale"]

        inhibited: list[str] = []
        if apply_inventory_limits:
            if self.protected_dry_mass_kg is not None and vehicle_mass_kg <= (
                self.protected_dry_mass_kg + 1e-12
            ):
                inhibited.append("protected_dry_mass_kg")
            if self.propellant_reserve_kg is not None and states[
                self.propellant_state_name
            ] <= (self.propellant_reserve_kg + 1e-12):
                inhibited.append(self.propellant_state_name)

        if level_index == OFF_LEVEL or inhibited:
            return ElectricPropulsionEvaluation(
                active=False,
                level_index=level_index,
                level_id=None,
                thrust_n=0.0,
                force_vector_n=_ZERO3,
                acceleration_m_s2=_ZERO3,
                mass_derivative_kg_s=0.0,
                additional_state_derivatives=((self.propellant_state_name, 0.0),),
                available_power_w=available_power_w,
                consumed_power_w=0.0,
                power_deficit_w=0.0,
                specific_impulse_s=None,
                parameters=tuple(parameters.items()),
                inhibited_by=tuple(inhibited),
                partials=ElectricPropulsionPartials(
                    acceleration_wrt_mass=_ZERO3,
                    acceleration_wrt_direction=(_ZERO3, _ZERO3, _ZERO3),
                    acceleration_wrt_available_power=_ZERO3,
                    mass_rate_wrt_available_power=0.0,
                    acceleration_wrt_parameters=tuple((name, _ZERO3) for name in parameters),
                    mass_rate_wrt_parameters=tuple((name, 0.0) for name in parameters),
                ),
            )

        point = self.throttle_law.operating_point(level_index, effective_power)
        duty = self.duty_cycle
        base_thrust = point.thrust_n * duty
        base_flow = point.mass_flow_kg_s * duty
        thrust = base_thrust * thrust_scale
        flow = base_flow * flow_scale
        force = tuple(thrust * value for value in unit)
        acceleration = tuple(value / vehicle_mass_kg for value in force)
        direction_gain = thrust / (vehicle_mass_kg * norm)
        thrust_power_slope = point.thrust_wrt_available_power * duty * thrust_scale
        flow_power_slope = point.mass_flow_wrt_available_power * duty * flow_scale
        return ElectricPropulsionEvaluation(
            active=True,
            level_index=level_index,
            level_id=point.level_id,
            thrust_n=thrust,
            force_vector_n=force,  # type: ignore[arg-type]
            acceleration_m_s2=acceleration,  # type: ignore[arg-type]
            mass_derivative_kg_s=-flow,
            additional_state_derivatives=((self.propellant_state_name, -flow),),
            available_power_w=available_power_w,
            consumed_power_w=point.input_power_w * duty,
            power_deficit_w=max(0.0, point.required_power_w - effective_power),
            specific_impulse_s=(
                thrust / (STANDARD_GRAVITY_M_S2 * flow) if flow > 0.0 else None
            ),
            parameters=tuple(parameters.items()),
            inhibited_by=(),
            partials=ElectricPropulsionPartials(
                acceleration_wrt_mass=tuple(
                    -value / vehicle_mass_kg for value in acceleration
                ),  # type: ignore[arg-type]
                acceleration_wrt_direction=tuple(
                    tuple(
                        direction_gain
                        * ((1.0 if row == column else 0.0) - unit[row] * unit[column])
                        for column in range(3)
                    )
                    for row in range(3)
                ),  # type: ignore[arg-type]
                acceleration_wrt_available_power=tuple(
                    thrust_power_slope * value / vehicle_mass_kg for value in unit
                ),  # type: ignore[arg-type]
                mass_rate_wrt_available_power=-flow_power_slope,
                acceleration_wrt_parameters=(
                    (
                        "thrust_scale",
                        tuple(base_thrust * value / vehicle_mass_kg for value in unit),
                    ),  # type: ignore[arg-type]
                    ("mass_flow_scale", _ZERO3),
                ),
                mass_rate_wrt_parameters=(
                    ("thrust_scale", 0.0),
                    ("mass_flow_scale", -base_flow),
                ),
            ),
        )
