"""Isentropic relations for a calorically perfect gas."""

from __future__ import annotations

from dataclasses import dataclass
from math import asin, atan, degrees, isfinite, radians, sin, sqrt

from rocket_propulsion.core.errors import DomainError
from rocket_propulsion.core.numerics import bisect_monotonic


@dataclass(frozen=True, slots=True)
class IsentropicResult:
    """Dimensionless isentropic properties associated with a Mach number."""

    mach: float
    mach_angle_deg: float | None
    prandtl_meyer_deg: float | None
    pressure_total_ratio: float
    density_total_ratio: float
    temperature_total_ratio: float
    pressure_critical_ratio: float
    density_critical_ratio: float
    temperature_critical_ratio: float
    area_critical_ratio: float


def _validate_gamma(gamma: float) -> None:
    if not isfinite(gamma) or gamma <= 1.0:
        raise DomainError("Gamma must be finite and greater than 1.")


def _validate_mach(mach: float) -> None:
    if not isfinite(mach) or mach <= 0.0:
        raise DomainError("Mach number must be finite and greater than 0.")


def calculate_isentropic(mach: float, gamma: float = 1.4) -> IsentropicResult:
    """Calculate perfect-gas isentropic flow relations.

    Uses ``T/T0 = 1 / (1 + (gamma-1) M²/2)`` and its pressure and density
    power laws. Starred properties are referenced to the sonic state at the
    same stagnation condition. Mach and Prandtl-Meyer angles are defined only
    for ``M >= 1`` and are returned as ``None`` below sonic speed.

    Args:
        mach: Local Mach number, strictly positive.
        gamma: Ratio of specific heats, strictly greater than one.

    Returns:
        An immutable set of flow ratios and angles in degrees.
    """

    _validate_gamma(gamma)
    _validate_mach(mach)

    energy_factor = 1.0 + 0.5 * (gamma - 1.0) * mach**2
    temperature_total = energy_factor**-1.0
    pressure_total = energy_factor ** (-gamma / (gamma - 1.0))
    density_total = energy_factor ** (-1.0 / (gamma - 1.0))
    critical_factor = (gamma + 1.0) / (2.0 * energy_factor)
    area_ratio = (
        (1.0 / mach)
        * (2.0 * energy_factor / (gamma + 1.0))
        ** ((gamma + 1.0) / (2.0 * (gamma - 1.0)))
    )

    if mach >= 1.0:
        root = sqrt(max(0.0, mach**2 - 1.0))
        mach_angle = degrees(asin(1.0 / mach))
        prandtl_meyer = degrees(
            sqrt((gamma + 1.0) / (gamma - 1.0))
            * atan(sqrt((gamma - 1.0) / (gamma + 1.0)) * root)
            - atan(root)
        )
    else:
        mach_angle = None
        prandtl_meyer = None

    return IsentropicResult(
        mach=mach,
        mach_angle_deg=mach_angle,
        prandtl_meyer_deg=prandtl_meyer,
        pressure_total_ratio=pressure_total,
        density_total_ratio=density_total,
        temperature_total_ratio=temperature_total,
        pressure_critical_ratio=critical_factor ** (gamma / (gamma - 1.0)),
        density_critical_ratio=critical_factor ** (1.0 / (gamma - 1.0)),
        temperature_critical_ratio=critical_factor,
        area_critical_ratio=area_ratio,
    )


def mach_from_isentropic(
    input_kind: str,
    value: float,
    gamma: float = 1.4,
) -> float:
    """Recover Mach number from an isentropic property.

    Supported inputs are ``mach``, ``mach_angle_deg``, ``prandtl_meyer_deg``,
    ``pressure_total_ratio``, ``density_total_ratio``,
    ``temperature_total_ratio``, ``pressure_critical_ratio``,
    ``density_critical_ratio``, ``temperature_critical_ratio``,
    ``area_ratio_subsonic``, and ``area_ratio_supersonic``.

    Area ratio is double-valued and therefore requires an explicit branch.
    Numerical inversions use robust monotonic bisection.
    """

    _validate_gamma(gamma)
    if not isfinite(value):
        raise DomainError("Input value must be finite.")
    if input_kind == "mach":
        _validate_mach(value)
        return value
    if input_kind == "mach_angle_deg":
        if not 0.0 < value <= 90.0:
            raise DomainError("Mach angle must be in the interval (0, 90].")
        return 1.0 / sin(radians(value))
    if input_kind == "prandtl_meyer_deg":
        maximum = degrees(0.5 * 3.141592653589793 * (sqrt((gamma + 1) / (gamma - 1)) - 1))
        if not 0.0 <= value < maximum:
            raise DomainError(f"Prandtl-Meyer angle must be in [0, {maximum:g}).")
        return bisect_monotonic(
            lambda candidate: calculate_isentropic(candidate, gamma).prandtl_meyer_deg or 0.0,
            value,
            1.0,
            1000.0,
        )

    field_by_input = {
        "pressure_total_ratio": "pressure_total_ratio",
        "density_total_ratio": "density_total_ratio",
        "temperature_total_ratio": "temperature_total_ratio",
        "pressure_critical_ratio": "pressure_critical_ratio",
        "density_critical_ratio": "density_critical_ratio",
        "temperature_critical_ratio": "temperature_critical_ratio",
    }
    if input_kind in field_by_input:
        if value <= 0.0:
            raise DomainError("Flow ratios must be greater than zero.")
        field = field_by_input[input_kind]
        return bisect_monotonic(
            lambda candidate: getattr(calculate_isentropic(candidate, gamma), field),
            value,
            1.0e-8,
            1000.0,
        )

    if input_kind in {"area_ratio_subsonic", "area_ratio_supersonic"}:
        if value < 1.0:
            raise DomainError("A/A* must be greater than or equal to 1.")
        if value == 1.0:
            return 1.0
        bounds = (1.0e-8, 1.0) if input_kind.endswith("subsonic") else (1.0, 1000.0)
        return bisect_monotonic(
            lambda candidate: calculate_isentropic(candidate, gamma).area_critical_ratio,
            value,
            *bounds,
        )

    raise DomainError(f"Unsupported isentropic input kind: {input_kind!r}.")

