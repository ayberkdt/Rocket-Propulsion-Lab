"""Fanno and Rayleigh one-dimensional perfect-gas relations."""

from dataclasses import dataclass
from math import isfinite, log, sqrt

from rocket_propulsion.core.errors import DomainError


@dataclass(frozen=True, slots=True)
class FannoResult:
    """Local Fanno-flow properties referenced to the sonic state."""

    mach: float
    temperature_critical_ratio: float
    pressure_critical_ratio: float
    density_critical_ratio: float
    velocity_critical_ratio: float
    stagnation_pressure_critical_ratio: float
    friction_parameter_to_sonic: float
    entropy_difference_over_r: float


@dataclass(frozen=True, slots=True)
class RayleighResult:
    """Local Rayleigh-flow properties referenced to the sonic state."""

    mach: float
    stagnation_temperature_critical_ratio: float
    temperature_critical_ratio: float
    pressure_critical_ratio: float
    density_critical_ratio: float
    velocity_critical_ratio: float
    stagnation_pressure_critical_ratio: float
    entropy_difference_over_r: float


def _validate(mach: float, gamma: float) -> None:
    if not isfinite(mach) or mach <= 0.0:
        raise DomainError("Mach number must be finite and greater than zero.")
    if not isfinite(gamma) or gamma <= 1.0:
        raise DomainError("Gamma must be finite and greater than one.")


def calculate_fanno(mach: float, gamma: float = 1.4) -> FannoResult:
    """Calculate adiabatic constant-area flow with wall friction.

    ``4 f L*/D`` is the remaining Darcy-style friction parameter required to
    reach the sonic state. Starred properties share the same mass flux and
    stagnation temperature.
    """

    _validate(mach, gamma)
    mach_squared = mach**2
    denominator = 2.0 + (gamma - 1.0) * mach_squared
    temperature = (gamma + 1.0) / denominator
    pressure = sqrt((gamma + 1.0) / denominator) / mach
    density = sqrt(denominator / (gamma + 1.0)) / mach
    velocity = mach * sqrt(temperature)
    stagnation_pressure = (
        ((denominator / (gamma + 1.0)) ** ((gamma + 1.0) / (2.0 * (gamma - 1.0))))
        / mach
    )
    friction = (1.0 - mach_squared) / (gamma * mach_squared) + (
        (gamma + 1.0) / (2.0 * gamma)
    ) * log(((gamma + 1.0) * mach_squared) / denominator)

    return FannoResult(
        mach=mach,
        temperature_critical_ratio=temperature,
        pressure_critical_ratio=pressure,
        density_critical_ratio=density,
        velocity_critical_ratio=velocity,
        stagnation_pressure_critical_ratio=stagnation_pressure,
        friction_parameter_to_sonic=max(0.0, friction),
        entropy_difference_over_r=log(stagnation_pressure),
    )


def calculate_rayleigh(mach: float, gamma: float = 1.4) -> RayleighResult:
    """Calculate frictionless constant-area flow with heat transfer.

    Results are referenced to the sonic Rayleigh-line state for the same mass
    and momentum flux. Positive heat addition drives either branch toward
    choking, subject to the maximum stagnation-temperature point.
    """

    _validate(mach, gamma)
    mach_squared = mach**2
    momentum_factor = 1.0 + gamma * mach_squared
    temperature = ((gamma + 1.0) ** 2 * mach_squared) / momentum_factor**2
    pressure = (gamma + 1.0) / momentum_factor
    density = momentum_factor / ((gamma + 1.0) * mach_squared)
    velocity = ((gamma + 1.0) * mach_squared) / momentum_factor
    stagnation_temperature = (
        (gamma + 1.0) * mach_squared * (2.0 + (gamma - 1.0) * mach_squared)
    ) / momentum_factor**2
    stagnation_pressure = pressure * (
        (2.0 + (gamma - 1.0) * mach_squared) / (gamma + 1.0)
    ) ** (gamma / (gamma - 1.0))
    entropy_difference = log(pressure) - (gamma / (gamma - 1.0)) * log(temperature)

    return RayleighResult(
        mach=mach,
        stagnation_temperature_critical_ratio=stagnation_temperature,
        temperature_critical_ratio=temperature,
        pressure_critical_ratio=pressure,
        density_critical_ratio=density,
        velocity_critical_ratio=velocity,
        stagnation_pressure_critical_ratio=stagnation_pressure,
        entropy_difference_over_r=entropy_difference,
    )

