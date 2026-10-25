"""Normal-shock relations for a calorically perfect gas.

Steady, adiabatic, inviscid Rankine-Hugoniot jump relations [1]_.

References
----------
.. [1] NACA-TR-1135. Ames Research Staff, "Equations, Tables, and Charts for
       Compressible Flow", NACA Report 1135, 1953.
"""

from dataclasses import dataclass
from math import isfinite

from rocket_propulsion.core.errors import DomainError


@dataclass(frozen=True, slots=True)
class NormalShockResult:
    """Flow-property ratios across a stationary normal shock [1]_.

    Attributes
    ----------
    upstream_mach, downstream_mach : float
        Mach numbers immediately before and after the shock.
    pressure_ratio, density_ratio, temperature_ratio : float
        Downstream-to-upstream static ratios ``p2/p1``, ``rho2/rho1`` and
        ``T2/T1``.
    stagnation_pressure_ratio : float
        Total-pressure recovery ``p02/p01``.

    References
    ----------
    .. [1] NACA-TR-1135. Ames Research Staff, "Equations, Tables, and Charts for
           Compressible Flow", NACA Report 1135, 1953.
    """

    upstream_mach: float
    downstream_mach: float
    pressure_ratio: float
    density_ratio: float
    temperature_ratio: float
    stagnation_pressure_ratio: float


def calculate_normal_shock(upstream_mach: float, gamma: float = 1.4) -> NormalShockResult:
    """Calculate Rankine-Hugoniot relations across a normal shock.

    Parameters
    ----------
    upstream_mach : float
        Mach number immediately before the shock, finite and greater than one.
    gamma : float, default 1.4
        Ratio of specific heats, finite and strictly greater than one.

    Returns
    -------
    NormalShockResult
        Downstream Mach number plus static and stagnation property ratios.

    Raises
    ------
    DomainError
        If ``gamma`` is invalid or the upstream flow is not supersonic, so the
        state cannot contain a normal shock.

    Notes
    -----
    For a calorically perfect gas in steady, adiabatic, inviscid flow [1]_:

    * ``M2^2 = (1 + (gamma-1) M1^2/2) / (gamma M1^2 - (gamma-1)/2)``
    * ``p2/p1 = (2 gamma M1^2 - (gamma-1)) / (gamma+1)``
    * ``rho2/rho1 = (gamma+1) M1^2 / ((gamma-1) M1^2 + 2)``
    * ``T2/T1 = (p2/p1) / (rho2/rho1)``
    * ``p02/p01 = (rho2/rho1)^(gamma/(gamma-1)) (p1/p2)^(1/(gamma-1))``

    References
    ----------
    .. [1] NACA-TR-1135. Ames Research Staff, "Equations, Tables, and Charts for
           Compressible Flow", NACA Report 1135, 1953.
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

