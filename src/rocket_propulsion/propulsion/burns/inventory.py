"""Constant-flow tank inventory ledgers with exact binary-input arithmetic.

Each tank-depleting stream of an :class:`EngineOperatingPoint` drains one
declared tank at its own mass-flow rate, so a tank's remaining mass after a
time ``t`` at constant flow is ``m(t) = m_loaded - t * sum(mdot_i)`` [1]_.
The helpers use the explicit stream rates, not rounded declared engine totals,
and evaluate every sum and reserve decision on the exact rational values of
the supplied binary64 inputs [2]_ before a single checked float export.

Arithmetic exactness does not establish physical-model or measurement
validity; the ledgers assume constant flow and no unaccounted losses.

References
----------
.. [1] SUTTON-BIBLARZ-2017. G. P. Sutton and O. Biblarz, "Rocket Propulsion
       Elements", 9th ed., Wiley, 2017.
.. [2] PYTHON-3.12-FRACTION. Python 3.12, "fractions": float inputs retain their
       exact binary value as a rational,
       https://docs.python.org/3.12/library/fractions.html.
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
    """Return ``value`` as a :class:`PropellantRole` without altering the caller's record."""
    if not isinstance(value, str):
        raise DomainError(f"{label} must be a recognized propellant role.", field=label)
    try:
        return PropellantRole(value)
    except ValueError as error:
        raise DomainError(f"{label} must be a recognized propellant role.", field=label) from error


def _checked_streams(operating_point):
    """Validate the operating-point streams and return exact rates of tank-depleting paths."""
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
    """Validate tanks in order, rejecting duplicate identifiers and invalid reserves."""
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
    """Sum the exact stream rates per tank before any binary64 export."""
    result = {}
    for _, tank, _, _, rate in streams:
        result[tank] = result.get(tank, Fraction(0)) + rate
    return result


def _checked_inventory(tanks, operating_point):
    """Resolve every tracked path, including zero flow, to one compatible tank.

    Raises
    ------
    DomainError
        If a stream names an unknown tank or its propellant or role differs
        from the tank's.
    """
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
    """Return the total tank-depleting mass flow drawn from each tank.

    Parameters
    ----------
    operating_point : EngineOperatingPoint
        Operating point whose ``tank_depleting`` streams are grouped by
        ``tank_id``.

    Returns
    -------
    dict of str to float
        Mass flow per tank, in kg/s.  The exact per-tank sum of stream rates
        [1]_ is rounded once on export.

    Raises
    ------
    DomainError
        If a stream is invalid or a grouped rate is not representable as a
        finite float.

    References
    ----------
    .. [1] SUTTON-BIBLARZ-2017. G. P. Sutton and O. Biblarz, "Rocket Propulsion
           Elements", 9th ed., Wiley, 2017.
    """
    return {name: _finite_rational(rate, f"Tank {name!r} flow")
            for name, rate in _group_rates(_checked_streams(operating_point)).items()}


def validate_tank_coverage(
    tanks: tuple[TankDefinition, ...], operating_point: EngineOperatingPoint
) -> None:
    """Require a valid, unambiguous tank for every tank-depleting stream.

    Parameters
    ----------
    tanks : tuple of TankDefinition
        Tanks with unique identifiers and ``0 <= reserve <= loaded`` mass.
    operating_point : EngineOperatingPoint
        Operating point whose ``tank_depleting`` streams are checked.

    Raises
    ------
    DomainError
        If a tank or stream is invalid, a tank identifier is duplicated, or a
        stream names an unknown tank or one with a different propellant or
        role.
    """
    _checked_inventory(tanks, operating_point)


def tank_duration_limit(
    tanks: tuple[TankDefinition, ...], operating_point: EngineOperatingPoint
) -> tuple[float, str | None]:
    """Return the greatest safe binary64 burn time and the first binding tank.

    Parameters
    ----------
    tanks : tuple of TankDefinition
        Tanks with loaded and reserve masses, in kg.
    operating_point : EngineOperatingPoint
        Constant-flow operating point drawing on the tanks.

    Returns
    -------
    duration_s : float
        Time, in s, at which the first tank reaches its reserve at constant
        flow, ``min_i (m_loaded,i - m_reserve,i) / mdot_i`` [1]_.  ``inf`` when
        the cutoff exceeds the largest finite float.
    binding_tank : str or None
        Identifier of the tank that sets the limit; ``None`` when no tank has
        a positive drain (the duration is then ``inf``).

    Raises
    ------
    DomainError
        If the tanks or streams are invalid or inconsistent.

    Notes
    -----
    The exact rational cutoff is rounded downward to the adjacent
    representable float when needed [2]_, so a query at the returned finite
    time cannot cross a reserve.  A cutoff below the smallest positive float
    returns zero.  Exact ties retain the supplied tank order.

    References
    ----------
    .. [1] SUTTON-BIBLARZ-2017. G. P. Sutton and O. Biblarz, "Rocket Propulsion
           Elements", 9th ed., Wiley, 2017.
    .. [2] PYTHON-3.12-NEXTAFTER. Python 3.12, ``math.nextafter``: adjacent
           representable floating-point values,
           https://docs.python.org/3.12/library/math.html#math.nextafter.
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
    """Return each tank's remaining mass after a constant-flow time.

    Parameters
    ----------
    tanks : tuple of TankDefinition
        Tanks with loaded and reserve masses, in kg.
    operating_point : EngineOperatingPoint
        Constant-flow operating point drawing on the tanks.
    time_s : float
        Finite, non-negative elapsed burn time, in s.

    Returns
    -------
    tuple of (str, float)
        ``(tank_id, m_loaded - mdot * t)`` [1]_ in kg, in the supplied tank
        order.

    Raises
    ------
    DomainError
        If ``time_s`` is negative or non-finite, the inputs are invalid, a
        tank falls below its reserve, or a non-zero mass is outside the
        binary64 output range.

    Notes
    -----
    Reserve decisions use the exact remaining mass before rounding; no
    tolerance admits a real reserve crossing.

    References
    ----------
    .. [1] SUTTON-BIBLARZ-2017. G. P. Sutton and O. Biblarz, "Rocket Propulsion
           Elements", 9th ed., Wiley, 2017.
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

