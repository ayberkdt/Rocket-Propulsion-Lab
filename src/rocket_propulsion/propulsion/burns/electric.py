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

import json
from collections.abc import Mapping
from dataclasses import dataclass, field
from fractions import Fraction
from itertools import pairwise
from math import isfinite

from rocket_propulsion.core.errors import DomainError, InputError
from rocket_propulsion.core.validation import (
    integer,
    read_json_document,
    real,
    real_mapping,
    records,
    strings,
    text,
)
from rocket_propulsion.propulsion.performance import STANDARD_GRAVITY_M_S2

from .direction import _direction_geometry
from .mass_closure import _MASS_CLOSURE_POLICY, _inventory_ledger
from .numeric_outputs import _finite_rational, _fraction
from .power_documents import _table_document, _power_document, _POWER_METADATA
from .exchange_documents import _metadata
from .power_modes import _POWER_MODE_POLICY, _mode_roots, _select_mode
from .propagator import (
    BurnEventKind,
    BurnEventSurface,
    EstimableParameter,
    PropellantMassClosure,
    _validated_states,
)
from .serialization import _require_numerical_policies, canonical_sha256

ELECTRIC_THROTTLE_TABLE_SCHEMA = "rocket_propulsion_electric_throttle_table_v1"
POWER_LIMITED_PROPULSION_SCHEMA = "rocket_propulsion_power_limited_propulsion_v1"
ELECTRIC_PROPULSION_REFERENCE_IDS = (
    "NASA-20090004685",
    "OREKIT-13.1.5-PROPULSION-MODEL",
)

OFF_LEVEL = -1




def _positive(value: float, label: str) -> float:
    number = real(value, label)
    if number <= 0.0:
        raise DomainError(f"{label} must be finite and greater than zero.")
    return number


def _nonnegative(value: float, label: str) -> float:
    number = real(value, label)
    if number < 0.0:
        raise DomainError(f"{label} must be finite and non-negative.")
    return number


def _strings_blank(value: object, label: str) -> tuple[str, ...]:
    return strings(value, label, allow_blank=True)


def _set(instance: object, name: str, value: object) -> None:
    object.__setattr__(instance, name, value)




def _sha256_hex(value: str, label: str) -> None:
    text(value, label)
    if len(value) != 64 or any(character not in "0123456789abcdef" for character in value):
        raise DomainError(f"{label} must be a lower-case SHA-256 hex digest.")


def _level_efficiency(level):
    """Exact jet/input power ratio of already validated binary level scalars."""
    return _fraction(level.thrust_n)**2 / (
        2 * _fraction(level.mass_flow_kg_s) * _fraction(level.input_power_w))


@dataclass(frozen=True, slots=True)
class ElectricThrottleLevel:
    """One supplied electric-thruster operating point.

    ``input_power_w`` is the power drawn from the spacecraft bus by the
    thruster string (power-processing unit input), which is the quantity the
    power model must supply.  The derived total efficiency
    ``T^2 / (2 mdot P)`` is the ratio of jet power to input power and must not
    exceed the retained numerical consistency ceiling ``1.0 + 1e-12``.
    This check is not a calibrated efficiency uncertainty or qualification
    record. Ratios use exact binary-input arithmetic before checked output.

    References
    ----------
    NASA-20090004685: https://ntrs.nasa.gov/citations/20090004685
    """

    level_id: str
    input_power_w: float
    thrust_n: float
    mass_flow_kg_s: float

    def __post_init__(self) -> None:
        _set(self, "level_id", text(self.level_id, "Throttle level identifier"))
        _set(self, "input_power_w", _positive(self.input_power_w, "Throttle level input power"))
        _set(self, "thrust_n", _positive(self.thrust_n, "Throttle level thrust"))
        _set(self, "mass_flow_kg_s", _positive(self.mass_flow_kg_s, "Throttle level mass flow"))
        if _level_efficiency(self) > _fraction(1.0 + 1e-12):
            raise DomainError(
                f"Throttle level {self.level_id!r} implies jet power above input power."
            )

    @property
    def specific_impulse_s(self) -> float:
        """Return ``T / (g0 mdot)`` in seconds."""

        return _finite_rational(
            _fraction(self.thrust_n) / (_fraction(STANDARD_GRAVITY_M_S2) * _fraction(self.mass_flow_kg_s)),
            "Throttle level specific impulse")

    @property
    def total_efficiency(self) -> float:
        """Return jet power over input power, ``T^2 / (2 mdot P)``."""

        return _finite_rational(_level_efficiency(self), "Throttle level total efficiency")


@dataclass(frozen=True, slots=True)
class ElectricThrottleTable:
    """Tabulated points with caller-supplied source identifiers and digest.

    Content identities support reproducible reopening. Source metadata alone
    does not authenticate data or prove a qualification/acceptance campaign.
    No level is interpolated or extrapolated.

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
        _set(self, "schema", text(self.schema, "Throttle-table schema"))
        if self.schema != ELECTRIC_THROTTLE_TABLE_SCHEMA:
            raise DomainError("Unsupported electric throttle-table schema.")
        _set(self, "levels", records(self.levels, ElectricThrottleLevel, "Throttle levels"))
        _set(self, "reference_ids", _strings_blank(self.reference_ids, "Throttle-table references"))
        _set(self, "source_id", text(self.source_id, "Throttle-table source identifier"))
        _set(self, "source_sha256", text(self.source_sha256, "Throttle-table source hash"))
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

        return canonical_sha256(self.to_dict(include_hash=False))

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
            payload["table_hash"] = canonical_sha256(payload)
        return payload

    def to_json(self, *, indent: int | None = 2) -> str:
        """Serialize the complete hash-bound table."""

        return json.dumps(
            self.to_dict(), allow_nan=False, ensure_ascii=False, indent=indent, sort_keys=True
        )

    @classmethod
    def from_dict(cls, payload: Mapping[str, object]) -> ElectricThrottleTable:
        """Reopen a table and reject malformed or hash-altered content."""

        _table_document(payload)
        try:
            levels_payload = payload["levels"]
            if not isinstance(levels_payload, list):
                raise TypeError
            table = cls(
                schema=payload["schema"],
                source_id=payload["source_id"],
                source_sha256=payload["source_sha256"],
                reference_ids=payload["reference_ids"],
                levels=tuple(
                    ElectricThrottleLevel(
                        level_id=item["level_id"],
                        input_power_w=item["input_power_w"],
                        thrust_n=item["thrust_n"],
                        mass_flow_kg_s=item["mass_flow_kg_s"],
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

        return read_json_document(text, "Electric throttle-table artifact", cls.from_dict)


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
class _ExactLevelPoint:
    """Private binary-input ledger; export only at a public record boundary."""

    level_index: int
    level_id: str
    thrust_n: Fraction
    mass_flow_kg_s: Fraction
    input_power_w: Fraction
    required_power_w: Fraction
    thrust_wrt_available_power: Fraction
    mass_flow_wrt_available_power: Fraction

    def export(self) -> LevelOperatingPoint:
        return LevelOperatingPoint(
            self.level_index, self.level_id,
            _finite_rational(self.thrust_n, "Continuous thrust"),
            _finite_rational(self.mass_flow_kg_s, "Continuous mass flow"),
            _finite_rational(self.input_power_w, "Throttle input power"),
            _finite_rational(self.required_power_w, "Throttle required power"),
            _finite_rational(self.thrust_wrt_available_power, "Continuous thrust/power slope"),
            _finite_rational(self.mass_flow_wrt_available_power, "Continuous flow/power slope"),
        )


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
    threshold_ratio: tuple[int, int] | None = field(default=None, repr=False)

    def __post_init__(self) -> None:
        _set(self, "event_id", text(self.event_id, "Power-switch event identifier"))
        _set(self, "crossing", text(self.crossing, "Power-switch crossing"))
        _set(self, "threshold_w", real(self.threshold_w, "Power-switch threshold"))
        _set(self, "from_level", integer(self.from_level, "Power-switch from_level"))
        _set(self, "to_level", integer(self.to_level, "Power-switch to_level"))
        if self.from_level < OFF_LEVEL or self.to_level < OFF_LEVEL:
            raise DomainError("Power-switch levels must be -1 (off) or a level index.")
        if not isfinite(self.threshold_w):
            raise DomainError("Power-switch threshold must be finite.")
        if self.crossing not in {"increasing", "decreasing", "either"}:
            raise DomainError("Power-switch crossing must be increasing, decreasing or either.")
        if self.crossing == "either" and self.to_level != self.from_level:
            raise DomainError("Only a mode-preserving kink surface may fire in either direction.")
        if self.crossing == "increasing" and self.to_level <= self.from_level:
            raise DomainError("An increasing power event must raise the mode index.", field="Power-switch target")
        if self.crossing == "decreasing" and self.to_level >= self.from_level:
            raise DomainError("A decreasing power event must lower the mode index.", field="Power-switch target")
        if self.threshold_ratio is None:
            exact = _fraction(self.threshold_w)
        else:
            if type(self.threshold_ratio) is not tuple or len(self.threshold_ratio) != 2:
                raise DomainError("Power-switch threshold ratio requires two integers.", field="Power-switch threshold ratio")
            numerator, denominator = (integer(value, "Power-switch threshold ratio") for value in self.threshold_ratio)
            if denominator <= 0:
                raise DomainError("Power-switch threshold denominator must be positive.", field="Power-switch threshold ratio")
            exact = Fraction(numerator, denominator)
            if _finite_rational(exact, "Power-switch threshold") != self.threshold_w:
                raise DomainError("Power-switch display does not match its exact threshold.", field="Power-switch threshold")
        _set(self, "threshold_ratio", (exact.numerator, exact.denominator))

    @property
    def threshold_exact_w(self) -> Fraction:
        """Exact threshold; threshold_w is its nearest-float display."""
        return Fraction(*self.threshold_ratio)

    def value(self, available_power_w: float) -> float:
        """Evaluate ``available_power_w - threshold_w``."""

        return _finite_rational(_fraction(real(available_power_w, "Available power")) - self.threshold_exact_w,
                                "Power-switch residual")


def _mode_data(law):
    """Declared inputs of the two built-in relay laws, without float export."""
    if isinstance(law, DiscreteThrottlePath):
        return (tuple(_fraction(level.input_power_w) for level in law.levels),
                _fraction(law.up_switch_hysteresis_w), tuple(level.level_id for level in law.levels))
    return ((_fraction(law.minimum_power_w),), _fraction(law.start_hysteresis_w), (law.law_id,))


def _law_surfaces(law, current_level, margin=Fraction(0), *, include_kink=False):
    required, hysteresis, identifiers = _mode_data(law)
    saturation = _fraction(law.maximum_power_w) if include_kink and isinstance(law, FixedEfficiencyThrottle) else None
    specs = _mode_roots(required, hysteresis, current_level, identifiers, margin, saturation=saturation)
    return tuple(PowerSwitchSurface(name, _finite_rational(threshold, "Power-switch threshold"),
                    crossing, source, target, threshold_ratio=(threshold.numerator, threshold.denominator))
                 for name, threshold, crossing, source, target in specs)


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
        if not isinstance(self.table, ElectricThrottleTable):
            raise DomainError("Throttle path requires an ElectricThrottleTable.")
        _set(self, "level_ids", _strings_blank(self.level_ids, "Throttle-path level identifiers"))
        _set(self, "path_id", text(self.path_id, "Throttle path identifier"))
        if not self.level_ids:
            raise DomainError("A throttle path requires at least one level.")
        if len(self.level_ids) != len(set(self.level_ids)):
            raise DomainError("A throttle path may use each level only once.")
        if not self.path_id.strip():
            raise DomainError("Throttle path identifier must not be empty.")
        _set(
            self,
            "up_switch_hysteresis_w",
            _nonnegative(self.up_switch_hysteresis_w, "Up-switch hysteresis"),
        )
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

        if not isinstance(table, ElectricThrottleTable):
            raise DomainError("Throttle path requires an ElectricThrottleTable.", field="Throttle table")
        objective = text(objective, "Throttle-path objective")
        if objective not in {"thrust", "specific_impulse"}:
            raise DomainError("Throttle-path objective must be thrust or specific_impulse.")

        def score(level: ElectricThrottleLevel) -> Fraction:
            return _fraction(level.thrust_n) if objective == "thrust" else _fraction(level.thrust_n) / _fraction(level.mass_flow_kg_s)

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

        required, hysteresis, _ = _mode_data(self)
        return _select_mode(_fraction(real(available_power_w, "Available power")), required, hysteresis, current_level)

    def switching_surfaces(self, current_level: int) -> tuple[PowerSwitchSurface, ...]:
        """Return the down- and up-switch roots for the current level."""

        return _law_surfaces(self, current_level)

    def operating_point(self, level_index: int, available_power_w: float) -> LevelOperatingPoint:
        """Return the constant operating point of ``level_index``."""

        return _operating_point_ledger(self, level_index, real(available_power_w, "Available power")).export()

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


def _fixed_gain(law):
    """Retain the exact binary-input thrust/power ratio until each output."""
    return 2 * _fraction(law.total_efficiency) / (
        _fraction(STANDARD_GRAVITY_M_S2) * _fraction(law.specific_impulse_s))


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
        _set(self, "total_efficiency", real(self.total_efficiency, "Total efficiency"))
        if not 0.0 < self.total_efficiency <= 1.0:
            raise DomainError("Total efficiency must be in (0, 1].")
        _set(self, "specific_impulse_s", _positive(self.specific_impulse_s, "Specific impulse"))
        _set(self, "minimum_power_w", _positive(self.minimum_power_w, "Minimum input power"))
        _set(self, "maximum_power_w", _positive(self.maximum_power_w, "Maximum input power"))
        if self.maximum_power_w <= self.minimum_power_w:
            raise DomainError("Maximum input power must exceed the minimum.")
        _set(
            self, "start_hysteresis_w", _nonnegative(self.start_hysteresis_w, "Start hysteresis")
        )
        _set(self, "law_id", text(self.law_id, "Throttle law identifier"))

    @property
    def level_count(self) -> int:
        return 1

    @property
    def thrust_per_watt_n_w(self) -> float:
        """Return ``2 eta / (g0 Isp)``."""

        return _finite_rational(_fixed_gain(self), "Continuous thrust/power gain")

    def select_level(self, available_power_w: float, current_level: int | None) -> int:
        """Return ``0`` when the thruster runs and ``-1`` when it is off."""

        required, hysteresis, _ = _mode_data(self)
        return _select_mode(_fraction(real(available_power_w, "Available power")), required, hysteresis, current_level)

    def switching_surfaces(self, current_level: int) -> tuple[PowerSwitchSurface, ...]:
        """Return the shutdown root when on, or the start root when off."""

        return _law_surfaces(self, current_level)

    def kink_surfaces(self) -> tuple[PowerSwitchSurface, ...]:
        """Return the non-mode-changing saturation root where ``dT/dP`` jumps."""

        return tuple(surface for surface in _law_surfaces(self, 0, include_kink=True) if surface.crossing == "either")

    def operating_point(self, level_index: int, available_power_w: float) -> LevelOperatingPoint:
        """Return thrust and flow using ``min(P, P_max)`` and their power slopes."""

        return _operating_point_ledger(self, level_index, real(available_power_w, "Available power")).export()

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


def _operating_point_ledger(law, level_index, available_power_w):
    """Keep built-in law values exact across duty/driver scaling in a model."""
    level_index = integer(level_index, "Throttle level index")
    available_power_w = _fraction(available_power_w, "Available power")
    if not 0 <= level_index < law.level_count:
        raise DomainError("Throttle level index is outside the law.", field="Throttle level index")
    if isinstance(law, DiscreteThrottlePath):
        level = law.levels[level_index]
        return _ExactLevelPoint(level_index, level.level_id, _fraction(level.thrust_n),
            _fraction(level.mass_flow_kg_s), _fraction(level.input_power_w),
            _fraction(level.input_power_w), Fraction(0), Fraction(0))
    saturated = available_power_w >= law.maximum_power_w
    power = _fraction(law.maximum_power_w if saturated else max(available_power_w, 0.0))
    gain = _fixed_gain(law)
    velocity = _fraction(STANDARD_GRAVITY_M_S2) * _fraction(law.specific_impulse_s)
    # At zero use the right derivative; below zero clamped power is constant.
    slope = Fraction(0) if saturated or available_power_w < 0 else gain
    return _ExactLevelPoint(0, law.law_id, gain * power, gain * power / velocity,
        power, _fraction(law.minimum_power_w), slope, slope / velocity)


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
        _set(self, "schema", text(self.schema, "Power-limited model schema"))
        if self.schema != POWER_LIMITED_PROPULSION_SCHEMA:
            raise DomainError("Unsupported power-limited propulsion schema.")
        if not isinstance(self.throttle_law, (DiscreteThrottlePath, FixedEfficiencyThrottle)):
            raise DomainError("Throttle law must be a DiscreteThrottlePath or FixedEfficiencyThrottle.")
        _set(self, "model_id", text(self.model_id, "Model identifier"))
        _set(self, "propellant_state_name", text(self.propellant_state_name, "Propellant state name"))
        if self.propellant_reserve_kg is not None:
            _set(
                self,
                "propellant_reserve_kg",
                _nonnegative(self.propellant_reserve_kg, "Propellant reserve"),
            )
        if self.protected_dry_mass_kg is not None:
            _set(
                self,
                "protected_dry_mass_kg",
                _positive(self.protected_dry_mass_kg, "Protected dry mass"),
            )
        _set(self, "power_margin_w", _nonnegative(self.power_margin_w, "Power margin"))
        _set(self, "duty_cycle", real(self.duty_cycle, "Duty cycle"))
        _set(
            self, "parameters", records(self.parameters, EstimableParameter, "Parameters")
        )
        if not 0.0 < self.duty_cycle <= 1.0:
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

        return _finite_rational(
            _fraction(available_power_w, "Available power") - _fraction(self.power_margin_w),
            "Effective power")

    def select_level(self, available_power_w: float, current_level: int | None = None) -> int:
        """Initialize (``current_level=None``) or update the throttle mode.

        References
        ----------
        NASA-20090004685: https://ntrs.nasa.gov/citations/20090004685
        """

        if type(self.throttle_law) in (DiscreteThrottlePath, FixedEfficiencyThrottle):
            required, hysteresis, _ = _mode_data(self.throttle_law)
            power = _fraction(real(available_power_w, "Available power")) - _fraction(self.power_margin_w)
            return _select_mode(power, required, hysteresis, current_level)
        return self.throttle_law.select_level(self.effective_power_w(available_power_w), current_level)

    def switching_surfaces(self, current_level: int) -> tuple[PowerSwitchSurface, ...]:
        """Return mode-changing roots expressed in *available* power.

        References
        ----------
        Orekit EventDetector 13.1.5 API:
        https://www.orekit.org/static/apidocs/org/orekit/propagation/events/EventDetector.html
        """

        if type(self.throttle_law) in (DiscreteThrottlePath, FixedEfficiencyThrottle):
            return _law_surfaces(self.throttle_law, current_level, _fraction(self.power_margin_w), include_kink=True)
        surfaces = self.throttle_law.switching_surfaces(current_level)
        if isinstance(self.throttle_law, FixedEfficiencyThrottle) and current_level == 0:
            surfaces = surfaces + self.throttle_law.kink_surfaces()
        return tuple(
            PowerSwitchSurface(
                surface.event_id,
                _finite_rational(surface.threshold_exact_w + _fraction(self.power_margin_w), "Power-switch threshold"),
                surface.crossing,
                surface.from_level,
                surface.to_level,
                threshold_ratio=((surface.threshold_exact_w + _fraction(self.power_margin_w)).numerator,
                                 (surface.threshold_exact_w + _fraction(self.power_margin_w)).denominator),
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

        supplied = real_mapping(overrides, "Parameter overrides")
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
            "mass_closure_numerics": _MASS_CLOSURE_POLICY,
            "power_mode_numerics": _POWER_MODE_POLICY,
            "direction_numerics": "scaled_hypot_rational_projection_and_signed_tangents_v1",
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
            payload["model_hash"] = canonical_sha256(payload)
        return payload

    def to_dict(self) -> dict[str, object]:
        """Return the complete, reopenable model document.

        The manifest (with ``model_hash``) is extended with the full throttle
        table of a discrete path.  :meth:`from_dict` refuses a document whose
        recomputed hashes differ.

        References
        ----------
        NASA-20090004685: https://ntrs.nasa.gov/citations/20090004685
        """

        document = self.manifest()
        if isinstance(self.throttle_law, DiscreteThrottlePath):
            document["throttle_table"] = self.throttle_law.table.to_dict()
        return document

    def to_json(self, *, indent: int | None = 2) -> str:
        """Serialize the complete model document."""

        return json.dumps(
            self.to_dict(), allow_nan=False, ensure_ascii=False, indent=indent, sort_keys=True
        )

    @classmethod
    def from_dict(cls, payload: Mapping[str, object]) -> PowerLimitedPropulsionModel:
        """Rebuild a model and verify the table and model hashes.

        References
        ----------
        NASA-20090004685: https://ntrs.nasa.gov/citations/20090004685
        """

        _power_document(payload)
        try:
            law_payload = payload["throttle_law"]
            if law_payload["law"] == "discrete_throttle_path":
                table = ElectricThrottleTable.from_dict(payload["throttle_table"])
                if table.table_hash != law_payload["table_hash"]:
                    raise InputError("Throttle-path table hash does not match the table.")
                law: DiscreteThrottlePath | FixedEfficiencyThrottle = DiscreteThrottlePath(
                    table,
                    tuple(law_payload["level_ids"]),
                    law_payload["up_switch_hysteresis_w"],
                    law_payload["path_id"],
                )
            elif law_payload["law"] == "fixed_efficiency":
                law = FixedEfficiencyThrottle(
                    law_payload["total_efficiency"],
                    law_payload["specific_impulse_s"],
                    law_payload["minimum_power_w"],
                    law_payload["maximum_power_w"],
                    law_payload["start_hysteresis_w"],
                    law_payload["law_id"],
                )
            else:
                raise InputError("Unknown throttle law in model document.")
            model = cls(
                throttle_law=law,
                propellant_state_name=payload["propellant_state_name"],
                propellant_reserve_kg=payload["propellant_reserve_kg"],
                protected_dry_mass_kg=payload["protected_dry_mass_kg"],
                power_margin_w=payload["power_margin_w"],
                duty_cycle=payload["duty_cycle"],
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
                model_id=payload["model_id"],
                schema=payload["schema"],
            )
        except InputError:
            raise
        except (KeyError, TypeError, ValueError, AttributeError) as error:
            raise InputError("Malformed power-limited model document.") from error
        if payload.get("model_hash") != model.model_hash:
            raise InputError("Power-limited model hash does not match its contents.")
        _require_numerical_policies(payload, model.manifest(include_hash=False))
        _metadata(payload, model.manifest(include_hash=False), _POWER_METADATA)
        return model

    @classmethod
    def from_json(cls, text: str) -> PowerLimitedPropulsionModel:
        """Deserialize and integrity-check a model document."""

        return read_json_document(text, "Power-limited model document", cls.from_dict)

    @property
    def model_hash(self) -> str:
        """Content identity of the throttle law, limits and semantics."""

        return canonical_sha256(self.manifest(include_hash=False))

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

        vehicle_mass_kg = _positive(vehicle_mass_kg, "Vehicle mass")
        non_propellant_mass_kg = _nonnegative(non_propellant_mass_kg, "Non-propellant mass")
        states = _validated_states(additional_states)
        if self.propellant_state_name not in states:
            raise DomainError(f"Missing propellant state {self.propellant_state_name!r}.")
        return PropellantMassClosure(
            vehicle_mass_kg, non_propellant_mass_kg,
            *_inventory_ledger(vehicle_mass_kg, non_propellant_mass_kg,
                               (states[self.propellant_state_name],)),
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

        level_index = integer(level_index, "Throttle level index")
        if level_index < OFF_LEVEL or level_index >= self.throttle_law.level_count:
            raise DomainError("Throttle level index is outside the throttle law.")
        vehicle_mass_kg = _positive(vehicle_mass_kg, "Vehicle mass")
        available_power_w = real(available_power_w, "Available power")
        if not isinstance(apply_inventory_limits, bool):
            raise DomainError("apply_inventory_limits must be a boolean.")
        effective_power = (_fraction(available_power_w) - _fraction(self.power_margin_w)
                           if type(self.throttle_law) in (DiscreteThrottlePath, FixedEfficiencyThrottle)
                           else self.effective_power_w(available_power_w))
        geometry = _direction_geometry(direction)
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

        if type(self.throttle_law) in (DiscreteThrottlePath, FixedEfficiencyThrottle):
            point = _operating_point_ledger(self.throttle_law, level_index, effective_power)
        else:
            # Retain subclass callback behavior at its existing public boundary.
            point = self.throttle_law.operating_point(level_index, effective_power)
        duty = self.duty_cycle
        base_thrust = _fraction(point.thrust_n) * _fraction(duty)
        base_flow = _fraction(point.mass_flow_kg_s) * _fraction(duty)
        exact_thrust = base_thrust * _fraction(thrust_scale)
        thrust = _finite_rational(exact_thrust, "Axial thrust")
        flow = base_flow * _fraction(flow_scale)
        drain = _finite_rational(-flow, "Electric propellant drain")
        force = geometry.project(exact_thrust, label="Force")
        acceleration = geometry.project(exact_thrust, vehicle_mass_kg, label="Acceleration")
        thrust_power_slope = (_fraction(point.thrust_wrt_available_power)
                              * _fraction(duty) * _fraction(thrust_scale))
        flow_power_slope = (_fraction(point.mass_flow_wrt_available_power)
                            * _fraction(duty) * _fraction(flow_scale))
        return ElectricPropulsionEvaluation(
            active=True,
            level_index=level_index,
            level_id=point.level_id,
            thrust_n=thrust,
            force_vector_n=force,  # type: ignore[arg-type]
            acceleration_m_s2=acceleration,  # type: ignore[arg-type]
            mass_derivative_kg_s=drain,
            additional_state_derivatives=((self.propellant_state_name, drain),),
            available_power_w=available_power_w,
            consumed_power_w=_finite_rational(_fraction(point.input_power_w) * _fraction(duty), "Consumed electric power"),
            power_deficit_w=_finite_rational(max(0, _fraction(point.required_power_w) - _fraction(effective_power)), "Electric power deficit"),
            specific_impulse_s=(
                _finite_rational(exact_thrust / (_fraction(STANDARD_GRAVITY_M_S2) * flow),
                                 "Electric specific impulse") if flow > 0 else None
            ),
            parameters=tuple(parameters.items()),
            inhibited_by=(),
            partials=ElectricPropulsionPartials(
                acceleration_wrt_mass=geometry.project(
                    -exact_thrust, _fraction(vehicle_mass_kg)**2, label="Mass acceleration partial"),
                acceleration_wrt_direction=geometry.jacobian(exact_thrust, vehicle_mass_kg),
                acceleration_wrt_available_power=geometry.project(
                    thrust_power_slope, vehicle_mass_kg, label="Power acceleration partial"),
                mass_rate_wrt_available_power=_finite_rational(-flow_power_slope, "Electric flow/power partial"),
                acceleration_wrt_parameters=(
                    (
                        "thrust_scale",
                        geometry.project(base_thrust, vehicle_mass_kg,
                                         label="Thrust scale acceleration partial"),
                    ),  # type: ignore[arg-type]
                    ("mass_flow_scale", _ZERO3),
                ),
                mass_rate_wrt_parameters=(
                    ("thrust_scale", 0.0),
                    ("mass_flow_scale", _finite_rational(-base_flow, "Electric flow scale partial")),
                ),
            ),
        )
