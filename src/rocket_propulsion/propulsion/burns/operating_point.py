"""Provider protocol and convenience constructors for constant performance."""

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
    """Small performance-provider boundary consumed by future burn solvers."""

    provider_id: str
    provider_version: str

    def operating_point(self) -> EngineOperatingPoint:
        """Return the resolved delivered performance law for L0."""

    def discontinuities_s(self, definition: BurnDefinition) -> tuple[float, ...]:
        """Return exact schedule discontinuities relative to command time."""


def _positive(value: float, label: str) -> None:
    if not isfinite(value) or value <= 0.0:
        raise DomainError(f"{label} must be finite and greater than zero.")


@dataclass(frozen=True, slots=True)
class ConstantPerformanceProvider:
    """L0 provider that returns one immutable delivered operating point."""

    point: EngineOperatingPoint
    provider_id: str = "constant-performance"
    provider_version: str = "1.0"

    def operating_point(self) -> EngineOperatingPoint:
        """Return the immutable point represented by this provider."""

        return self.point

    def discontinuities_s(self, definition: BurnDefinition) -> tuple[float, ...]:
        """A constant provider has no internal discontinuities."""

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

        ``oxidizer_fuel_ratio`` applies to total tank drain. Optional dump flow
        is subtracted from the selected role's main-chamber stream, preserving
        both the total mixture ratio and the system-Isp mass law.
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
        """Build a single-stream monopropellant operating point."""

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
    """Evaluate a supported provider through the common provider-facing API."""

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

