"""Tank-ledger helpers for propulsion-only burn analysis."""

from __future__ import annotations

from math import inf

from rocket_propulsion.core.errors import DomainError

from .models import EngineOperatingPoint, TankDefinition


def flow_by_tank(operating_point: EngineOperatingPoint) -> dict[str, float]:
    """Aggregate tank-depleting streams by tank identifier."""

    result: dict[str, float] = {}
    for stream in operating_point.streams:
        if stream.tank_depleting:
            result[stream.tank_id] = result.get(stream.tank_id, 0.0) + stream.mass_flow_kg_s
    return result


def validate_tank_coverage(
    tanks: tuple[TankDefinition, ...], operating_point: EngineOperatingPoint
) -> None:
    """Require every depleting stream to resolve to a compatible tank."""

    tank_index = {tank.tank_id: tank for tank in tanks}
    for stream in operating_point.streams:
        if not stream.tank_depleting:
            continue
        tank = tank_index.get(stream.tank_id)
        if tank is None:
            raise DomainError(f"Stream {stream.stream_id!r} references unknown tank {stream.tank_id!r}.")
        if tank.propellant_key != stream.propellant_key:
            raise DomainError(
                f"Stream {stream.stream_id!r} propellant does not match tank {tank.tank_id!r}."
            )
        if tank.role is not stream.role:
            raise DomainError(f"Stream {stream.stream_id!r} role does not match its tank.")


def tank_duration_limit(
    tanks: tuple[TankDefinition, ...], operating_point: EngineOperatingPoint
) -> tuple[float, str | None]:
    """Return the first reserve-limited burn time and binding tank."""

    rates = flow_by_tank(operating_point)
    shortest = inf
    binding: str | None = None
    for tank in tanks:
        rate = rates.get(tank.tank_id, 0.0)
        if rate <= 0.0:
            continue
        duration = tank.expendable_mass_kg / rate
        if duration < shortest:
            shortest = duration
            binding = tank.tank_id
    if any(tank_id not in rates for tank_id in (flow.tank_id for flow in operating_point.streams if flow.tank_depleting)):
        raise DomainError("One or more tank-depleting streams are not covered by inventory.")
    return shortest, binding


def tank_masses_at(
    tanks: tuple[TankDefinition, ...], operating_point: EngineOperatingPoint, time_s: float
) -> tuple[tuple[str, float], ...]:
    """Return ordered tank inventory after an exact constant-flow interval."""

    rates = flow_by_tank(operating_point)
    masses: list[tuple[str, float]] = []
    for tank in tanks:
        remaining = tank.loaded_mass_kg - rates.get(tank.tank_id, 0.0) * time_s
        tolerance = max(1e-12, tank.loaded_mass_kg * 1e-12)
        if remaining < tank.reserve_mass_kg - tolerance:
            raise DomainError(f"Tank {tank.tank_id!r} crossed its protected reserve.")
        masses.append((tank.tank_id, max(remaining, tank.reserve_mass_kg)))
    return tuple(masses)

