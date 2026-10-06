"""Normal-shock relations for a calorically perfect gas."""

from dataclasses import dataclass
from math import isfinite

from rocket_propulsion.core.errors import DomainError


@dataclass(frozen=True, slots=True)
class NormalShockResult:
    """Flow-property ratios across a stationary normal shock."""

    upstream_mach: float
    downstream_mach: float
    pressure_ratio: float
    density_ratio: float
    temperature_ratio: float
    stagnation_pressure_ratio: float


def calculate_normal_shock(upstream_mach: float, gamma: float = 1.4) -> NormalShockResult:
    """Calculate Rankine-Hugoniot relations across a normal shock.

    Args:
        upstream_mach: Supersonic Mach number immediately before the shock.
        gamma: Ratio of specific heats, strictly greater than one.

    Returns:
        Downstream Mach number plus static and stagnation property ratios.

    Raises:
        DomainError: If the state cannot contain a normal shock.
    """

    if not isfinite(gamma) or gamma <= 1.0:
        raise DomainError("Gamma must be finite and greater than 1.")
    if not isfinite(upstream_mach) or upstream_mach <= 1.0:
        raise DomainError("Upstream Mach number must be greater than 1.")

    mach_squared = upstream_mach**2
    downstream_squared = (
        1.0 + 0.5 * (gamma - 1.0) * mach_squared
    ) / (gamma * mach_squared - 0.5 * (gamma - 1.0))
    pressure_ratio = (2.0 * gamma * mach_squared - (gamma - 1.0)) / (gamma + 1.0)
    density_ratio = ((gamma + 1.0) * mach_squared) / (
        (gamma - 1.0) * mach_squared + 2.0
    )
    stagnation_pressure_ratio = (
        ((gamma + 1.0) * mach_squared)
        / ((gamma - 1.0) * mach_squared + 2.0)
    ) ** (gamma / (gamma - 1.0)) * (
        (gamma + 1.0) / (2.0 * gamma * mach_squared - (gamma - 1.0))
    ) ** (1.0 / (gamma - 1.0))

    return NormalShockResult(
        upstream_mach=upstream_mach,
        downstream_mach=downstream_squared**0.5,
        pressure_ratio=pressure_ratio,
        density_ratio=density_ratio,
        temperature_ratio=pressure_ratio / density_ratio,
        stagnation_pressure_ratio=stagnation_pressure_ratio,
    )

