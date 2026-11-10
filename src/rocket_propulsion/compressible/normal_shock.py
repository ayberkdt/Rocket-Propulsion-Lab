"""Normal-shock relations for a calorically perfect gas.

Steady, adiabatic, inviscid Rankine-Hugoniot jump relations [1]_.

References
----------
.. [1] NACA-TR-1135. Ames Research Staff, "Equations, Tables, and Charts for
       Compressible Flow", NACA Report 1135, 1953.
"""

from dataclasses import dataclass
from decimal import (
    ROUND_HALF_EVEN,
    Context,
    Decimal,
    DivisionByZero,
    InvalidOperation,
    Overflow,
    localcontext,
)
from math import isfinite

from rocket_propulsion.core.errors import DomainError
from rocket_propulsion.core.validation import real


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


def _positive_export(value: Decimal, field: str) -> float:
    result = float(value)
    if not isfinite(result) or result <= 0.0:
        raise DomainError("Normal shock result is outside the positive float range.",
                          code="numeric_range", field=field)
    return result


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
        If the inputs are invalid, the upstream flow is not supersonic, or
        an output has no positive finite binary64 representation.

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

    mach = real(upstream_mach, "upstream_mach")
    specific_heat_ratio = real(gamma, "gamma")
    if specific_heat_ratio <= 1.0:
        raise DomainError("Gamma must be greater than 1.", field="gamma")
    if mach <= 1.0:
        raise DomainError("Upstream Mach number must be greater than 1.", field="upstream_mach")

    # Evaluate complete ratios before export. In particular, the two individual
    # recovery powers may overflow and underflow while their product is finite.
    context = Context(prec=80, rounding=ROUND_HALF_EVEN,
                      Emin=-999999, Emax=999999, capitals=1, clamp=0, flags=[],
                      traps=[InvalidOperation, DivisionByZero, Overflow])
    with localcontext(context):
        g = Decimal.from_float(specific_heat_ratio)
        gm = g - 1
        m = Decimal.from_float(mach)
        mach_squared = m * m
        pressure_numerator = 2 * g * mach_squared - gm
        density_denominator = gm * mach_squared + 2
        pressure = pressure_numerator / (g + 1)
        density = (g + 1) * mach_squared / density_denominator
        downstream = (density_denominator / pressure_numerator).sqrt()
        pressure_ratio = _positive_export(pressure, "pressure_ratio")
        density_ratio = _positive_export(density, "density_ratio")
        temperature_ratio = _positive_export(pressure / density, "temperature_ratio")
        log_recovery = (g * density.ln() - pressure.ln()) / gm
        # Below this log the exact result rounds to zero even as a subnormal.
        recovery = log_recovery.exp() if log_recovery > -800 else Decimal(0)
        return NormalShockResult(
            upstream_mach=mach,
            downstream_mach=_positive_export(downstream, "downstream_mach"),
            pressure_ratio=pressure_ratio,
            density_ratio=density_ratio,
            temperature_ratio=temperature_ratio,
            stagnation_pressure_ratio=_positive_export(recovery, "stagnation_pressure_ratio"),
        )

