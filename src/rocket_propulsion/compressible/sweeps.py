"""Sampled compressible-flow curves for plots and parameter studies."""

from dataclasses import dataclass
from math import isfinite

from rocket_propulsion.core.errors import DomainError

from .isentropic import calculate_isentropic


@dataclass(frozen=True, slots=True)
class IsentropicSweepPoint:
    """Selected nondimensional isentropic properties at one Mach number."""

    mach: float
    pressure_total_ratio: float
    density_total_ratio: float
    temperature_total_ratio: float
    area_critical_ratio: float


def generate_isentropic_sweep(
    *,
    gamma: float = 1.4,
    minimum_mach: float = 0.05,
    maximum_mach: float = 5.0,
    point_count: int = 121,
) -> tuple[IsentropicSweepPoint, ...]:
    """Sample isentropic property curves on a uniform Mach grid."""

    if not isfinite(gamma) or gamma <= 1.0:
        raise DomainError("Gamma must be finite and greater than one.")
    if not isfinite(minimum_mach) or minimum_mach <= 0.0:
        raise DomainError("Minimum Mach must be finite and positive.")
    if not isfinite(maximum_mach) or maximum_mach <= minimum_mach:
        raise DomainError("Maximum Mach must be greater than minimum Mach.")
    if not 21 <= point_count <= 401:
        raise DomainError("Point count must be between 21 and 401.")

    points: list[IsentropicSweepPoint] = []
    for index in range(point_count):
        fraction = index / (point_count - 1)
        mach = minimum_mach + fraction * (maximum_mach - minimum_mach)
        result = calculate_isentropic(mach, gamma)
        points.append(
            IsentropicSweepPoint(
                mach=mach,
                pressure_total_ratio=result.pressure_total_ratio,
                density_total_ratio=result.density_total_ratio,
                temperature_total_ratio=result.temperature_total_ratio,
                area_critical_ratio=result.area_critical_ratio,
            )
        )
    return tuple(points)

