"""Attached oblique-shock relations for a calorically perfect gas."""

from dataclasses import dataclass
from math import asin, atan, cos, degrees, isfinite, pi, radians, sin, tan

from rocket_propulsion.core.errors import DomainError
from rocket_propulsion.core.numerics import bisect_monotonic

from .normal_shock import calculate_normal_shock


@dataclass(frozen=True, slots=True)
class ObliqueShockResult:
    """Flow state across an attached weak or strong oblique shock."""

    branch: str
    upstream_mach: float
    downstream_mach: float
    turning_angle_deg: float
    wave_angle_deg: float
    maximum_turning_angle_deg: float
    upstream_normal_mach: float
    downstream_normal_mach: float
    pressure_ratio: float
    density_ratio: float
    temperature_ratio: float
    stagnation_pressure_ratio: float


def _turning_angle_rad(upstream_mach: float, wave_angle_rad: float, gamma: float) -> float:
    """Return theta from the theta-beta-M relation for a candidate beta."""

    mach_normal_squared = upstream_mach**2 * sin(wave_angle_rad) ** 2
    numerator = 2.0 * (mach_normal_squared - 1.0)
    denominator = tan(wave_angle_rad) * (
        upstream_mach**2 * (gamma + cos(2.0 * wave_angle_rad)) + 2.0
    )
    return atan(numerator / denominator)


def _maximum_turn_state(upstream_mach: float, gamma: float) -> tuple[float, float]:
    """Return wave angle and maximum attached turning angle in radians."""

    lower = asin(1.0 / upstream_mach) + 1.0e-10
    upper = 0.5 * pi - 1.0e-10
    for _ in range(120):
        left = lower + (upper - lower) / 3.0
        right = upper - (upper - lower) / 3.0
        if _turning_angle_rad(upstream_mach, left, gamma) < _turning_angle_rad(
            upstream_mach, right, gamma
        ):
            lower = left
        else:
            upper = right
    beta = 0.5 * (lower + upper)
    return beta, _turning_angle_rad(upstream_mach, beta, gamma)


def calculate_oblique_shock(
    *,
    upstream_mach: float,
    turning_angle_deg: float,
    gamma: float = 1.4,
    branch: str = "weak",
) -> ObliqueShockResult:
    """Solve an attached oblique shock using the theta-beta-M relation.

    Args:
        upstream_mach: Supersonic Mach number before the wave.
        turning_angle_deg: Positive flow-deflection angle in degrees.
        gamma: Ratio of specific heats.
        branch: ``weak`` or ``strong`` attached solution.

    Raises:
        DomainError: If no attached solution exists at the requested deflection.
    """

    if not isfinite(upstream_mach) or upstream_mach <= 1.0:
        raise DomainError("Upstream Mach number must be greater than one.")
    if not isfinite(gamma) or gamma <= 1.0:
        raise DomainError("Gamma must be finite and greater than one.")
    if not isfinite(turning_angle_deg) or turning_angle_deg <= 0.0:
        raise DomainError("Turning angle must be finite and greater than zero.")
    normalized_branch = branch.strip().lower()
    if normalized_branch not in {"weak", "strong"}:
        raise DomainError("Oblique-shock branch must be weak or strong.")

    mach_angle = asin(1.0 / upstream_mach)
    peak_beta, maximum_theta = _maximum_turn_state(upstream_mach, gamma)
    target_theta = radians(turning_angle_deg)
    if target_theta > maximum_theta + 1.0e-12:
        raise DomainError(
            f"Detached shock: maximum attached turn is {degrees(maximum_theta):.6g} degrees."
        )

    if normalized_branch == "weak":
        lower, upper = mach_angle + 1.0e-10, peak_beta
    else:
        lower, upper = peak_beta, 0.5 * pi - 1.0e-10
    wave_angle = bisect_monotonic(
        lambda candidate: _turning_angle_rad(upstream_mach, candidate, gamma),
        target_theta,
        lower,
        upper,
    )
    upstream_normal = upstream_mach * sin(wave_angle)
    normal = calculate_normal_shock(upstream_normal, gamma)
    downstream_mach = normal.downstream_mach / sin(wave_angle - target_theta)

    return ObliqueShockResult(
        branch=normalized_branch,
        upstream_mach=upstream_mach,
        downstream_mach=downstream_mach,
        turning_angle_deg=turning_angle_deg,
        wave_angle_deg=degrees(wave_angle),
        maximum_turning_angle_deg=degrees(maximum_theta),
        upstream_normal_mach=upstream_normal,
        downstream_normal_mach=normal.downstream_mach,
        pressure_ratio=normal.pressure_ratio,
        density_ratio=normal.density_ratio,
        temperature_ratio=normal.temperature_ratio,
        stagnation_pressure_ratio=normal.stagnation_pressure_ratio,
    )

