"""Validated constant-flow ledgers, exact for supplied binary inputs only.

These helpers use explicit stream rates, not rounded declared engine totals.
Arithmetic exactness does not establish physical model or measurement validity.
"""

from __future__ import annotations

from fractions import Fraction
from math import inf, nextafter
from sys import float_info

from rocket_propulsion.core.errors import DomainError
from rocket_propulsion.core.validation import real, records, text

from .models import EngineOperatingPoint, PropellantFlow, PropellantRole, TankDefinition
from .numeric_outputs import _finite_rational


def _role(value, label):
    """Compare recognized role spellings without altering caller records."""
    if not isinstance(value, str):
        raise DomainError(f"{label} must be a recognized propellant role.", field=label)
    try:
        return PropellantRole(value)
    except ValueError as error:
        raise DomainError(f"{label} must be a recognized propellant role.", field=label) from error


def _checked_streams(operating_point):
    """Validate ledger inputs and retain exact rates for tracked paths only."""
    if not isinstance(operating_point, EngineOperatingPoint):
        raise DomainError("Operating point must be an EngineOperatingPoint.", field="Operating point")
    result = []
    for stream in records(operating_point.streams, PropellantFlow, "Operating point streams"):
        name = text(stream.stream_id, "Stream identifier")
        rate = real(stream.mass_flow_kg_s, f"Stream {name!r} mass flow")
        if rate < 0:
            raise DomainError(f"Stream {name!r} mass flow must be non-negative.", field="Stream mass flow")
        if not isinstance(stream.tank_depleting, bool):
            raise DomainError(f"Stream {name!r} tank_depleting must be boolean.", field="tank_depleting")
        if stream.tank_depleting:
            result.append((name, text(stream.tank_id, f"Stream {name!r} tank identifier"),
                           text(stream.propellant_key, f"Stream {name!r} propellant key"),
                           _role(stream.role, f"Stream {name!r} role"), Fraction(rate)))
    return tuple(result)


def _checked_tanks(tanks):
    """Copy and validate ordered inventory, rejecting ambiguous identifiers."""
    result = []
    seen = set()
    for tank in records(tanks, TankDefinition, "Tanks"):
        name = text(tank.tank_id, "Tank identifier")
        if name in seen:
            raise DomainError("Tank identifiers must be unique.", field="Tank identifier")
        seen.add(name)
        propellant = text(tank.propellant_key, f"Tank {name!r} propellant key")
        role = _role(tank.role, f"Tank {name!r} role")
        loaded = real(tank.loaded_mass_kg, f"Tank {name!r} loaded mass")
        reserve = real(tank.reserve_mass_kg, f"Tank {name!r} reserve mass")
        if not 0 <= reserve <= loaded:
            raise DomainError(f"Tank {name!r} requires 0 <= reserve <= loaded mass.", field="Tank inventory")
        result.append((name, propellant, role, Fraction(loaded), Fraction(reserve)))
    return tuple(result)


def _group_rates(streams):
    """Sum all paths before any binary64 rate export."""
    result = {}
    for _, tank, _, _, rate in streams:
        result[tank] = result.get(tank, Fraction(0)) + rate
    return result


def _checked_inventory(tanks, operating_point):
    """Resolve every tracked path, including zero flow, to one compatible tank."""
    inventory = _checked_tanks(tanks)
    streams = _checked_streams(operating_point)
    index = {row[0]: row for row in inventory}
    for name, tank_id, propellant, role, _ in streams:
        tank = index.get(tank_id)
        if tank is None:
            raise DomainError(f"Stream {name!r} references unknown tank {tank_id!r}.", field="Tank coverage")
        if propellant != tank[1]:
            raise DomainError(f"Stream {name!r} propellant does not match tank {tank_id!r}.", field="Tank coverage")
        if role != tank[2]:
            raise DomainError(f"Stream {name!r} role does not match its tank.", field="Tank coverage")
    return inventory, _group_rates(streams)


def flow_by_tank(operating_point: EngineOperatingPoint) -> dict[str, float]:
    """Export grouped tracked rates once; reject a non-finite grouped output."""
    return {name: _finite_rational(rate, f"Tank {name!r} flow")
            for name, rate in _group_rates(_checked_streams(operating_point)).items()}


def validate_tank_coverage(
    tanks: tuple[TankDefinition, ...], operating_point: EngineOperatingPoint
) -> None:
    """Require unambiguous, valid tanks for all tank-depleting streams."""
    _checked_inventory(tanks, operating_point)


def tank_duration_limit(
    tanks: tuple[TankDefinition, ...], operating_point: EngineOperatingPoint
) -> tuple[float, str | None]:
    """Return the greatest safe binary64 forward time and first binding tank.

    A finite physical cutoff rounds downward when needed, so a query at the
    returned finite time cannot cross a reserve. A cutoff below the smallest
    positive float returns zero; one above the largest finite float returns
    infinity with its binding tank. No positive drain returns (infinity, None).
    Exact ties retain the supplied tank order.
    """
    tanks, rates = _checked_inventory(tanks, operating_point)
    shortest = None
    binding = None
    for name, _, _, loaded, reserve in tanks:
        rate = rates.get(name, Fraction(0))
        if rate == 0:
            continue
        duration = (loaded - reserve) / rate
        if shortest is None or duration < shortest:
            shortest, binding = duration, name
    if shortest is None or shortest > Fraction(float_info.max):
        return inf, binding
    horizon = float(shortest)
    if Fraction(horizon) > shortest:
        horizon = nextafter(horizon, 0.)
    return horizon, binding


def tank_masses_at(
    tanks: tuple[TankDefinition, ...], operating_point: EngineOperatingPoint, time_s: float
) -> tuple[tuple[str, float], ...]:
    """Export inventories after a finite nonnegative constant-flow time.

    Reserve decisions precede rounding, with no tolerance admitting a real
    crossing. Nonzero masses outside binary64 output range fail explicitly.
    """
    time = real(time_s, "Inventory time")
    if time < 0:
        raise DomainError("Inventory time must be non-negative.", field="Inventory time")
    tanks, rates = _checked_inventory(tanks, operating_point)
    result = []
    for name, _, _, loaded, reserve in tanks:
        remaining = loaded - rates.get(name, Fraction(0)) * Fraction(time)
        if remaining < reserve:
            raise DomainError(f"Tank {name!r} crossed its protected reserve.", field="Tank inventory")
        result.append((name, _finite_rational(remaining, f"Tank {name!r} remaining mass")))
    return tuple(result)

