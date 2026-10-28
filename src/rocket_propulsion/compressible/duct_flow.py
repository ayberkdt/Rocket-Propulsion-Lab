"""Fanno and Rayleigh one-dimensional perfect-gas relations.

Constant-area, steady, one-dimensional flow of a calorically perfect gas with
wall friction (Fanno) or heat transfer (Rayleigh) [1]_.

References
----------
.. [1] NACA-TN-1336. B. L. Hicks, D. J. Montgomery and R. H. Wasserman, "On the
       One-Dimensional Theory of Steady Compressible Fluid Flow in Ducts with
       Friction and Heat Addition", NACA TN 1336, 1947.
"""

from dataclasses import dataclass
from math import isfinite, log, sqrt

from rocket_propulsion.core.errors import DomainError
from rocket_propulsion.core.validation import numeric_range_guard


@dataclass(frozen=True, slots=True)
class FannoResult:
    """Local Fanno-flow properties referenced to the sonic state [1]_.

    Attributes
    ----------
    mach : float
        Local Mach number.
    temperature_critical_ratio, pressure_critical_ratio : float
        Static ratios ``T/T*`` and ``p/p*``.
    density_critical_ratio, velocity_critical_ratio : float
        Ratios ``rho/rho*`` and ``V/V*``.
    stagnation_pressure_critical_ratio : float
        Total-pressure ratio ``p0/p0*``.
    friction_parameter_to_sonic : float
        Friction parameter ``4 f L*/D`` still required to reach the sonic
        state, clipped at zero.
    entropy_difference_over_r : float
        ``(s* - s)/R = ln(p0/p0*)``.

    References
    ----------
    .. [1] NACA-TN-1336. B. L. Hicks, D. J. Montgomery and R. H. Wasserman, "On the
           One-Dimensional Theory of Steady Compressible Fluid Flow in Ducts with
           Friction and Heat Addition", NACA TN 1336, 1947.
    """

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
    """Local Rayleigh-flow properties referenced to the sonic state [1]_.

    Attributes
    ----------
    mach : float
        Local Mach number.
    stagnation_temperature_critical_ratio : float
        Total-temperature ratio ``T0/T0*``.
    temperature_critical_ratio, pressure_critical_ratio : float
        Static ratios ``T/T*`` and ``p/p*``.
    density_critical_ratio, velocity_critical_ratio : float
        Ratios ``rho/rho*`` and ``V/V*``.
    stagnation_pressure_critical_ratio : float
        Total-pressure ratio ``p0/p0*``.
    entropy_difference_over_r : float
        ``(s* - s)/R = ln(p/p*) - gamma/(gamma-1) ln(T/T*)``, zero at
        the sonic state.

    References
    ----------
    .. [1] NACA-TN-1336. B. L. Hicks, D. J. Montgomery and R. H. Wasserman, "On the
           One-Dimensional Theory of Steady Compressible Fluid Flow in Ducts with
           Friction and Heat Addition", NACA TN 1336, 1947.
    """

    mach: float
    stagnation_temperature_critical_ratio: float
    temperature_critical_ratio: float
    pressure_critical_ratio: float
    density_critical_ratio: float
    velocity_critical_ratio: float
    stagnation_pressure_critical_ratio: float
    entropy_difference_over_r: float


def _validate(mach: float, gamma: float) -> None:
    """Raise :class:`DomainError` unless Mach is positive and gamma exceeds one."""
    if not isfinite(mach) or mach <= 0.0:
        raise DomainError("Mach number must be finite and greater than zero.")
    if not isfinite(gamma) or gamma <= 1.0:
        raise DomainError("Gamma must be finite and greater than one.")


@numeric_range_guard("mach")
def calculate_fanno(mach: float, gamma: float = 1.4) -> FannoResult:
    """Calculate adiabatic constant-area flow with wall friction (Fanno flow).

    Parameters
    ----------
    mach : float
        Local Mach number, finite and strictly positive.
    gamma : float, default 1.4
        Ratio of specific heats, finite and strictly greater than one.

    Returns
    -------
    FannoResult
        Property ratios to the sonic state and the friction parameter to choke.

    Raises
    ------
    DomainError
        If ``mach`` or ``gamma`` is outside its domain.

    Notes
    -----
    Starred properties share the same mass flux and stagnation temperature.
    With ``D = 2 + (gamma-1) M^2`` the relations of [1]_ are
    ``T/T* = (gamma+1)/D``, ``p/p* = sqrt((gamma+1)/D)/M`` and
    ``p0/p0* = (D/(gamma+1))^((gamma+1)/(2(gamma-1)))/M``. The remaining
    Darcy-style friction parameter to the sonic state is
    ``4 f L*/D = (1 - M^2)/(gamma M^2)
    + (gamma+1)/(2 gamma) ln((gamma+1) M^2 / D)``.

    References
    ----------
    .. [1] NACA-TN-1336. B. L. Hicks, D. J. Montgomery and R. H. Wasserman, "On the
           One-Dimensional Theory of Steady Compressible Fluid Flow in Ducts with
           Friction and Heat Addition", NACA TN 1336, 1947.
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


@numeric_range_guard("mach")
def calculate_rayleigh(mach: float, gamma: float = 1.4) -> RayleighResult:
    """Calculate frictionless constant-area flow with heat transfer (Rayleigh flow).

    Parameters
    ----------
    mach : float
        Local Mach number, finite and strictly positive.
    gamma : float, default 1.4
        Ratio of specific heats, finite and strictly greater than one.

    Returns
    -------
    RayleighResult
        Property ratios to the sonic Rayleigh-line state.

    Raises
    ------
    DomainError
        If ``mach`` or ``gamma`` is outside its domain.

    Notes
    -----
    Results are referenced to the sonic Rayleigh-line state for the same mass
    and momentum flux. With ``B = 1 + gamma M^2`` the relations of [1]_ are
    ``p/p* = (gamma+1)/B``, ``T/T* = (gamma+1)^2 M^2 / B^2`` and
    ``T0/T0* = (gamma+1) M^2 (2 + (gamma-1) M^2) / B^2``. Positive heat
    addition drives either branch toward choking, subject to the maximum
    stagnation-temperature point.

    References
    ----------
    .. [1] NACA-TN-1336. B. L. Hicks, D. J. Montgomery and R. H. Wasserman, "On the
           One-Dimensional Theory of Steady Compressible Fluid Flow in Ducts with
           Friction and Heat Addition", NACA TN 1336, 1947.
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

