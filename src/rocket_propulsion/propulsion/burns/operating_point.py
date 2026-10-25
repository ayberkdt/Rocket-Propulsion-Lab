"""Performance-provider protocol and constant-performance constructors.

A provider supplies the delivered operating point consumed by the burn
solvers.  The constant (L0) provider builds stream-closed points from delivered
thrust and system specific impulse through ``mdot = F / (g0 * Isp)`` [1]_,
with ``g0`` the registered ``STANDARD_GRAVITY_M_S2``.

References
----------
.. [1] NASA Glenn Research Center, "Specific Impulse".
       https://www1.grc.nasa.gov/beginners-guide-to-aeronautics/specific-impulse/
"""

from __future__ import annotations

from dataclasses import dataclass
from math import isfinite
from typing import Protocol

from rocket_propulsion.core.errors import DomainError
from rocket_propulsion.propulsion.performance import STANDARD_GRAVITY_M_S2

from .models import (
    BurnDefinition,
    BurnResult,
    EngineOperatingPoint,
    FlowDestination,
    ProfiledBurnResult,
    PropellantFlow,
    PropellantRole,
)


class PropulsionProvider(Protocol):
    """Structural contract of a performance provider consumed by burn solvers.

    Attributes
    ----------
    provider_id : str
        Stable provider identifier.
    provider_version : str
        Provider version string.

    Notes
    -----
    A provider resolves one delivered operating point whose total tank flow
    satisfies ``mdot = F / (g0 * Isp)`` [1]_ and reports the schedule
    discontinuities a solver must step across.

    References
    ----------
    .. [1] NASA Glenn Research Center, "Specific Impulse".
           https://www1.grc.nasa.gov/beginners-guide-to-aeronautics/specific-impulse/
    """

    provider_id: str
    provider_version: str

    def operating_point(self) -> EngineOperatingPoint:
        """Return the resolved delivered operating point."""

    def discontinuities_s(self, definition: BurnDefinition) -> tuple[float, ...]:
        """Return schedule discontinuity times relative to the command time, s."""


def _positive(value: float, label: str) -> None:
    """Raise `DomainError` unless ``value`` is finite and strictly positive."""
    if not isfinite(value) or value <= 0.0:
        raise DomainError(f"{label} must be finite and greater than zero.")


@dataclass(frozen=True, slots=True)
class ConstantPerformanceProvider:
    """L0 provider that returns one immutable delivered operating point.

    Parameters
    ----------
    point : EngineOperatingPoint
        Delivered operating point returned for every query.
    provider_id : str, default "constant-performance"
        Provider identifier.
    provider_version : str, default "1.0"
        Provider version string.

    Notes
    -----
    The classmethod constructors derive total tank flow from delivered thrust
    and system specific impulse, ``mdot = F / (g0 * Isp)`` [1]_, with ``g0`` the
    registered ``STANDARD_GRAVITY_M_S2``.

    References
    ----------
    .. [1] NASA Glenn Research Center, "Specific Impulse".
           https://www1.grc.nasa.gov/beginners-guide-to-aeronautics/specific-impulse/
    """

    point: EngineOperatingPoint
    provider_id: str = "constant-performance"
    provider_version: str = "1.0"

    def operating_point(self) -> EngineOperatingPoint:
        """Return the immutable operating point held by this provider."""

        return self.point

    def discontinuities_s(self, definition: BurnDefinition) -> tuple[float, ...]:
        """Return an empty tuple: a constant provider has no discontinuities.

        Parameters
        ----------
        definition : BurnDefinition
            Burn definition; unused.

        Returns
        -------
        tuple of float
            Always ``()``.
        """

        del definition
        return ()

    @classmethod
    def bipropellant(
        cls,
        *,
        delivered_thrust_n: float,
        system_specific_impulse_s: float,
        oxidizer_fuel_ratio: float,
        ideal_thrust_n: float | None = None,
        chamber_specific_impulse_s: float | None = None,
        dump_mass_flow_kg_s: float = 0.0,
        dump_role: PropellantRole = PropellantRole.FUEL,
        commanded_throttle: float = 1.0,
        realized_throttle: float = 1.0,
        active_engine_count: int = 1,
        chamber_pressure_pa: float = 1.0,
        ambient_pressure_pa: float = 0.0,
        oxidizer_tank_id: str = "oxidizer",
        fuel_tank_id: str = "fuel",
        oxidizer_key: str = "oxidizer",
        fuel_key: str = "fuel",
        name: str = "constant bipropellant point",
        source: str = "constant provider",
    ) -> ConstantPerformanceProvider:
        """Build a stream-closed bipropellant point from system performance.

        Parameters
        ----------
        delivered_thrust_n : float
            Delivered thrust ``F``, N; finite and positive.
        system_specific_impulse_s : float
            System specific impulse referenced to total tank drain, s; finite
            and positive.
        oxidizer_fuel_ratio : float
            Oxidizer-to-fuel mass ratio of total tank drain; finite and
            positive.
        ideal_thrust_n : float or None, default None
            Ideal thrust, N; ``None`` uses ``delivered_thrust_n``.
        chamber_specific_impulse_s : float or None, default None
            Chamber specific impulse, s; ``None`` derives
            ``F / (g0 * mdot_chamber)``.
        dump_mass_flow_kg_s : float, default 0.0
            Non-negative open-cycle dump flow that bypasses the chamber, kg/s.
        dump_role : PropellantRole, default PropellantRole.FUEL
            Component that feeds the dump; fuel or oxidizer.
        commanded_throttle, realized_throttle : float, default 1.0
            Commanded and realized throttle fractions.
        active_engine_count : int, default 1
            Number of active engines.
        chamber_pressure_pa : float, default 1.0
            Chamber pressure, Pa.
        ambient_pressure_pa : float, default 0.0
            Ambient pressure, Pa.
        oxidizer_tank_id, fuel_tank_id : str
            Source tank identifiers.
        oxidizer_key, fuel_key : str
            Propellant keys.
        name : str
            Operating-point name.
        source : str
            Provenance label.

        Returns
        -------
        ConstantPerformanceProvider
            Provider holding the stream-closed operating point.

        Raises
        ------
        DomainError
            If thrust, specific impulse or mixture ratio is not finite and
            positive, the dump flow is negative or non-finite, the dump role is
            neither fuel nor oxidizer, or the dump flow is not smaller than its
            component's tank flow.

        Notes
        -----
        Total tank flow is ``mdot = F / (g0 * Isp_sys)`` [1]_ with ``g0`` the
        registered ``STANDARD_GRAVITY_M_S2``; fuel flow is ``mdot / (1 + O/F)``
        [2]_ and the oxidizer receives the remainder.  ``oxidizer_fuel_ratio``
        applies to total tank drain.  An optional dump flow is subtracted from
        the selected component's main-chamber stream and carried as a
        non-thrust-contributing stream, so the total mixture ratio and the
        system-Isp mass relation are both preserved.

        References
        ----------
        .. [1] NASA Glenn Research Center, "Specific Impulse".
               https://www1.grc.nasa.gov/beginners-guide-to-aeronautics/specific-impulse/
        .. [2] SUTTON-BIBLARZ-2017. G. P. Sutton and O. Biblarz, "Rocket
               Propulsion Elements", 9th ed., Wiley, 2017.
        """

        _positive(delivered_thrust_n, "Delivered thrust")
        _positive(system_specific_impulse_s, "System specific impulse")
        _positive(oxidizer_fuel_ratio, "Oxidizer/fuel ratio")
        if not isfinite(dump_mass_flow_kg_s) or dump_mass_flow_kg_s < 0.0:
            raise DomainError("Dump mass flow must be finite and non-negative.")
        if dump_role not in {PropellantRole.FUEL, PropellantRole.OXIDIZER}:
            raise DomainError("Bipropellant dump role must be fuel or oxidizer.")

        total_flow = delivered_thrust_n / (
            STANDARD_GRAVITY_M_S2 * system_specific_impulse_s
        )
        fuel_flow = total_flow / (1.0 + oxidizer_fuel_ratio)
        oxidizer_flow = total_flow - fuel_flow
        role_flow = fuel_flow if dump_role is PropellantRole.FUEL else oxidizer_flow
        if dump_mass_flow_kg_s >= role_flow:
            raise DomainError("Dump flow must be smaller than its component tank flow.")
        chamber_fuel = fuel_flow - (
            dump_mass_flow_kg_s if dump_role is PropellantRole.FUEL else 0.0
        )
        chamber_oxidizer = oxidizer_flow - (
            dump_mass_flow_kg_s if dump_role is PropellantRole.OXIDIZER else 0.0
        )
        chamber_flow = chamber_fuel + chamber_oxidizer
        resolved_chamber_isp = chamber_specific_impulse_s
        if resolved_chamber_isp is None:
            resolved_chamber_isp = delivered_thrust_n / (
                STANDARD_GRAVITY_M_S2 * chamber_flow
            )
        streams = [
            PropellantFlow(
                "main-oxidizer",
                oxidizer_tank_id,
                oxidizer_key,
                PropellantRole.OXIDIZER,
                FlowDestination.CHAMBER,
                chamber_oxidizer,
            ),
            PropellantFlow(
                "main-fuel",
                fuel_tank_id,
                fuel_key,
                PropellantRole.FUEL,
                FlowDestination.CHAMBER,
                chamber_fuel,
            ),
        ]
        if dump_mass_flow_kg_s > 0.0:
            is_fuel = dump_role is PropellantRole.FUEL
            streams.append(
                PropellantFlow(
                    "open-cycle-dump",
                    fuel_tank_id if is_fuel else oxidizer_tank_id,
                    fuel_key if is_fuel else oxidizer_key,
                    dump_role,
                    FlowDestination.GAS_GENERATOR_DUMP,
                    dump_mass_flow_kg_s,
                    contributes_to_delivered_thrust=False,
                )
            )
        return cls(
            EngineOperatingPoint(
                name=name,
                commanded_throttle=commanded_throttle,
                realized_throttle=realized_throttle,
                active_engine_count=active_engine_count,
                ideal_thrust_n=(delivered_thrust_n if ideal_thrust_n is None else ideal_thrust_n),
                delivered_thrust_n=delivered_thrust_n,
                chamber_specific_impulse_s=resolved_chamber_isp,
                system_specific_impulse_s=system_specific_impulse_s,
                chamber_mass_flow_kg_s=chamber_flow,
                total_tank_flow_kg_s=total_flow,
                oxidizer_mass_flow_kg_s=oxidizer_flow,
                fuel_mass_flow_kg_s=fuel_flow,
                mixture_ratio=oxidizer_fuel_ratio,
                chamber_pressure_pa=chamber_pressure_pa,
                ambient_pressure_pa=ambient_pressure_pa,
                streams=tuple(streams),
                source=source,
            )
        )

    @classmethod
    def monopropellant(
        cls,
        *,
        delivered_thrust_n: float,
        system_specific_impulse_s: float,
        ideal_thrust_n: float | None = None,
        commanded_throttle: float = 1.0,
        realized_throttle: float = 1.0,
        active_engine_count: int = 1,
        chamber_pressure_pa: float = 1.0,
        ambient_pressure_pa: float = 0.0,
        tank_id: str = "propellant",
        propellant_key: str = "monopropellant",
        name: str = "constant monopropellant point",
        source: str = "constant provider",
    ) -> ConstantPerformanceProvider:
        """Build a single-stream monopropellant operating point.

        Parameters
        ----------
        delivered_thrust_n : float
            Delivered thrust ``F``, N; finite and positive.
        system_specific_impulse_s : float
            System specific impulse, s; finite and positive.  Chamber and
            system specific impulse are equal for a single chamber stream.
        ideal_thrust_n : float or None, default None
            Ideal thrust, N; ``None`` uses ``delivered_thrust_n``.
        commanded_throttle, realized_throttle : float, default 1.0
            Commanded and realized throttle fractions.
        active_engine_count : int, default 1
            Number of active engines.
        chamber_pressure_pa : float, default 1.0
            Chamber pressure, Pa.
        ambient_pressure_pa : float, default 0.0
            Ambient pressure, Pa.
        tank_id : str
            Source tank identifier.
        propellant_key : str
            Propellant key.
        name : str
            Operating-point name.
        source : str
            Provenance label.

        Returns
        -------
        ConstantPerformanceProvider
            Provider holding the single-stream operating point.

        Raises
        ------
        DomainError
            If thrust or specific impulse is not finite and positive.

        Notes
        -----
        Mass flow is ``mdot = F / (g0 * Isp)`` [1]_ with ``g0`` the registered
        ``STANDARD_GRAVITY_M_S2``.

        References
        ----------
        .. [1] NASA Glenn Research Center, "Specific Impulse".
               https://www1.grc.nasa.gov/beginners-guide-to-aeronautics/specific-impulse/
        """

        _positive(delivered_thrust_n, "Delivered thrust")
        _positive(system_specific_impulse_s, "System specific impulse")
        mass_flow = delivered_thrust_n / (
            STANDARD_GRAVITY_M_S2 * system_specific_impulse_s
        )
        return cls(
            EngineOperatingPoint(
                name=name,
                commanded_throttle=commanded_throttle,
                realized_throttle=realized_throttle,
                active_engine_count=active_engine_count,
                ideal_thrust_n=(delivered_thrust_n if ideal_thrust_n is None else ideal_thrust_n),
                delivered_thrust_n=delivered_thrust_n,
                chamber_specific_impulse_s=system_specific_impulse_s,
                system_specific_impulse_s=system_specific_impulse_s,
                chamber_mass_flow_kg_s=mass_flow,
                total_tank_flow_kg_s=mass_flow,
                oxidizer_mass_flow_kg_s=0.0,
                fuel_mass_flow_kg_s=0.0,
                mixture_ratio=None,
                chamber_pressure_pa=chamber_pressure_pa,
                ambient_pressure_pa=ambient_pressure_pa,
                streams=(
                    PropellantFlow(
                        "main",
                        tank_id,
                        propellant_key,
                        PropellantRole.MONOPROPELLANT,
                        FlowDestination.CHAMBER,
                        mass_flow,
                    ),
                ),
                source=source,
            )
        )


def simulate_burn(
    definition: BurnDefinition, provider: PropulsionProvider
) -> BurnResult | ProfiledBurnResult:
    """Evaluate a burn with a supported provider through the common entry point.

    Parameters
    ----------
    definition : BurnDefinition
        Burn command, duration and vehicle state.
    provider : PropulsionProvider
        Performance provider.  A ``ProfiledPerformanceProvider`` is routed to
        the profiled solver; any other provider must report no
        discontinuities and is routed to the constant (L0) solver.

    Returns
    -------
    BurnResult or ProfiledBurnResult
        Constant-burn result, or profiled-burn result for a profiled provider.

    Raises
    ------
    DomainError
        If a non-profiled provider reports schedule discontinuities, or the
        selected solver rejects the inputs.

    Notes
    -----
    Both solvers integrate thrust and propellant consumption linked by
    ``mdot = F / (g0 * Isp)`` [1]_.

    References
    ----------
    .. [1] NASA Glenn Research Center, "Specific Impulse".
           https://www1.grc.nasa.gov/beginners-guide-to-aeronautics/specific-impulse/
    """

    from .profile_integration import simulate_profiled_burn
    from .profiles import ProfiledPerformanceProvider

    if isinstance(provider, ProfiledPerformanceProvider):
        return simulate_profiled_burn(
            definition, provider.operating_point(), provider.schedule
        )

    from .integration import simulate_constant_burn

    if provider.discontinuities_s(definition):
        raise DomainError("L0 simulate_burn requires a provider without discontinuities.")
    return simulate_constant_burn(definition, provider.operating_point())

