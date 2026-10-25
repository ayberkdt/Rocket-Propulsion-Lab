"""Isentropic relations for a calorically perfect gas.

The module evaluates the one-dimensional isentropic ratios of NACA Report 1135
[1]_ and inverts them for Mach number.

References
----------
.. [1] NACA-TR-1135. Ames Research Staff, "Equations, Tables, and Charts for
       Compressible Flow", NACA Report 1135, 1953.
"""

from __future__ import annotations

from dataclasses import dataclass
from math import asin, atan, degrees, isfinite, radians, sin, sqrt

from rocket_propulsion.core.errors import DomainError
from rocket_propulsion.core.numerics import bisect_monotonic


@dataclass(frozen=True, slots=True)
class IsentropicResult:
    """Dimensionless isentropic properties associated with a Mach number [1]_.

    Attributes
    ----------
    mach : float
        Local Mach number.
    mach_angle_deg : float or None
        Mach angle ``asin(1/M)`` in degrees; ``None`` for subsonic flow.
    prandtl_meyer_deg : float or None
        Prandtl-Meyer function in degrees; ``None`` for subsonic flow.
    pressure_total_ratio, density_total_ratio, temperature_total_ratio : float
        Static-to-stagnation ratios ``p/p0``, ``rho/rho0`` and ``T/T0``.
    pressure_critical_ratio, density_critical_ratio, temperature_critical_ratio : float
        Ratios ``p/p*``, ``rho/rho*`` and ``T/T*`` to the sonic state at the
        same stagnation condition.
    area_critical_ratio : float
        Area ratio ``A/A*`` for the same mass flow.

    References
    ----------
    .. [1] NACA-TR-1135. Ames Research Staff, "Equations, Tables, and Charts for
           Compressible Flow", NACA Report 1135, 1953.
    """

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
    """Raise :class:`DomainError` unless ``gamma`` is finite and greater than one."""
    if not isfinite(gamma) or gamma <= 1.0:
        raise DomainError("Gamma must be finite and greater than 1.")


def _validate_mach(mach: float) -> None:
    """Raise :class:`DomainError` unless ``mach`` is finite and positive."""
    if not isfinite(mach) or mach <= 0.0:
        raise DomainError("Mach number must be finite and greater than 0.")


def calculate_isentropic(mach: float, gamma: float = 1.4) -> IsentropicResult:
    """Calculate perfect-gas isentropic flow relations at one Mach number.

    Parameters
    ----------
    mach : float
        Local Mach number, finite and strictly positive.
    gamma : float, default 1.4
        Ratio of specific heats, finite and strictly greater than one.

    Returns
    -------
    IsentropicResult
        Immutable set of flow ratios and, for ``M >= 1``, angles in degrees.

    Raises
    ------
    DomainError
        If ``mach`` or ``gamma`` is outside its domain.

    Notes
    -----
    With ``X = 1 + (gamma - 1) M^2 / 2`` the stagnation ratios are
    ``T/T0 = 1/X``, ``p/p0 = X^(-gamma/(gamma-1))`` and
    ``rho/rho0 = X^(-1/(gamma-1))`` [1]_. Starred properties are referenced to
    the sonic state at the same stagnation condition, and
    ``A/A* = (1/M) (2 X/(gamma+1))^((gamma+1)/(2(gamma-1)))``. The Mach angle
    ``mu = asin(1/M)`` and the Prandtl-Meyer function
    ``nu = sqrt((gamma+1)/(gamma-1)) atan(sqrt((gamma-1)/(gamma+1) (M^2-1)))
    - atan(sqrt(M^2-1))`` [1]_ are defined only for ``M >= 1`` and are
    returned as ``None`` below sonic speed.

    References
    ----------
    .. [1] NACA-TR-1135. Ames Research Staff, "Equations, Tables, and Charts for
           Compressible Flow", NACA Report 1135, 1953.
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
    """Recover Mach number from one isentropic property.

    Parameters
    ----------
    input_kind : str
        Name of the supplied property: ``mach``, ``mach_angle_deg``,
        ``prandtl_meyer_deg``, ``pressure_total_ratio``,
        ``density_total_ratio``, ``temperature_total_ratio``,
        ``pressure_critical_ratio``, ``density_critical_ratio``,
        ``temperature_critical_ratio``, ``area_ratio_subsonic`` or
        ``area_ratio_supersonic``.
    value : float
        Finite value of that property (angles in degrees).
    gamma : float, default 1.4
        Ratio of specific heats, finite and strictly greater than one.

    Returns
    -------
    float
        Mach number that reproduces ``value`` under
        :func:`calculate_isentropic`.

    Raises
    ------
    DomainError
        If ``input_kind`` is unsupported, ``value`` is not finite, or
        ``value`` lies outside the range of the selected property: Mach angle
        in ``(0, 90]`` degrees, Prandtl-Meyer angle below its ``M -> infinity``
        limit, positive flow ratios, or ``A/A* >= 1``.

    Notes
    -----
    The Mach angle is inverted in closed form, ``M = 1/sin(mu)``. All other
    properties are inverted by monotonic bisection of the forward relations of
    [1]_ on a bracket inside ``1e-8 <= M <= 1000``. The area ratio is
    double-valued, so the subsonic or supersonic branch must be named
    explicitly; ``A/A* = 1`` returns ``M = 1`` on either branch.

    References
    ----------
    .. [1] NACA-TR-1135. Ames Research Staff, "Equations, Tables, and Charts for
           Compressible Flow", NACA Report 1135, 1953.
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

