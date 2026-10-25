"""Ideal choked converging-diverging rocket-nozzle performance.

Quasi-one-dimensional, isentropic, frozen-property flow of a calorically
perfect gas [1]_, [2]_ with the choked-throat mass flow [3]_.

References
----------
.. [1] SUTTON-BIBLARZ-2017. G. P. Sutton and O. Biblarz, "Rocket Propulsion
       Elements", 9th ed., Wiley, 2017.
.. [2] NACA-TR-1135. Ames Research Staff, "Equations, Tables, and Charts for
       Compressible Flow", NACA Report 1135, 1953.
.. [3] NASA-GLENN-MASS-FLOW. NASA Glenn Research Center, isentropic mass-flow
       conservation and density/velocity context.
       https://www1.grc.nasa.gov/beginners-guide-to-aeronautics/mass-flow-rate-equations/
"""

from dataclasses import dataclass
from math import isfinite, sqrt

from rocket_propulsion.core.constants import STANDARD_GRAVITY_M_S2
from rocket_propulsion.core.errors import DomainError

from .isentropic import calculate_isentropic, mach_from_isentropic


@dataclass(frozen=True, slots=True)
class NozzleResult:
    """One-dimensional design-point state for an ideal choked nozzle [1]_.

    Attributes
    ----------
    throat_area_m2, exit_area_m2 : float
        Throat and exit areas in m^2.
    area_ratio : float
        Exit-to-throat area ratio.
    exit_mach, exit_temperature_k, exit_pressure_pa, exit_velocity_m_s : float
        Exit-plane state on the supersonic isentropic branch.
    mass_flow_kg_s : float
        Choked mass flow.
    momentum_thrust_n, pressure_thrust_n, thrust_n : float
        ``mdot V_e``, ``(p_e - p_a) A_e`` and their sum, in N.
    characteristic_velocity_m_s : float
        ``c* = p_c A_t / mdot``.
    thrust_coefficient : float
        ``C_F = F / (p_c A_t)``.
    specific_impulse_s : float
        ``F / (mdot g0)`` with ``g0`` from ``STANDARD_GRAVITY_M_S2``.
    pressure_regime : str
        ``"adapted"`` when ``|p_e - p_a| <= 0.01 max(p_e, 1 Pa)``, otherwise
        ``"underexpanded"`` or ``"overexpanded"``.

    References
    ----------
    .. [1] SUTTON-BIBLARZ-2017. G. P. Sutton and O. Biblarz, "Rocket Propulsion
           Elements", 9th ed., Wiley, 2017.
    """

    throat_area_m2: float
    exit_area_m2: float
    area_ratio: float
    exit_mach: float
    exit_temperature_k: float
    exit_pressure_pa: float
    exit_velocity_m_s: float
    mass_flow_kg_s: float
    momentum_thrust_n: float
    pressure_thrust_n: float
    thrust_n: float
    characteristic_velocity_m_s: float
    thrust_coefficient: float
    specific_impulse_s: float
    pressure_regime: str


def _positive(value: float, label: str) -> None:
    """Raise :class:`DomainError` unless ``value`` is finite and positive."""
    if not isfinite(value) or value <= 0.0:
        raise DomainError(f"{label} must be finite and greater than zero.")


def calculate_ideal_nozzle(
    *,
    chamber_pressure_pa: float,
    chamber_temperature_k: float,
    throat_area_m2: float,
    area_ratio: float,
    ambient_pressure_pa: float,
    gamma: float = 1.22,
    gas_constant_j_kg_k: float = 355.0,
) -> NozzleResult:
    """Calculate ideal, frozen-property, supersonic C-D nozzle performance.

    Parameters
    ----------
    chamber_pressure_pa : float
        Chamber stagnation pressure ``p_c`` in Pa, finite and positive.
    chamber_temperature_k : float
        Chamber stagnation temperature ``T_c`` in K, finite and positive.
    throat_area_m2 : float
        Throat area ``A_t`` in m^2, finite and positive.
    area_ratio : float
        Exit-to-throat area ratio ``A_e/A_t``, at least one.
    ambient_pressure_pa : float
        Ambient pressure ``p_a`` in Pa, finite and non-negative.
    gamma : float, default 1.22
        Frozen ratio of specific heats, finite and greater than one.
    gas_constant_j_kg_k : float, default 355.0
        Specific gas constant ``R`` in J/(kg K), finite and positive.

    Returns
    -------
    NozzleResult
        Exit state, mass flow, thrust components and performance figures.

    Raises
    ------
    DomainError
        If an input is outside its domain or the resulting thrust is not
        positive.

    Notes
    -----
    Chamber properties are treated as stagnation conditions. The throat is
    choked, so ``mdot = p_c A_t sqrt(gamma/(R T_c))
    (2/(gamma+1))^((gamma+1)/(2(gamma-1)))`` [3]_, and the exit Mach number is
    taken from the supersonic branch of the isentropic area ratio [2]_. Thrust
    is ``F = mdot V_e + (p_e - p_a) A_e`` [1]_. Boundary-layer, divergence,
    chemistry, separation and discharge losses are not included.

    References
    ----------
    .. [1] SUTTON-BIBLARZ-2017. G. P. Sutton and O. Biblarz, "Rocket Propulsion
           Elements", 9th ed., Wiley, 2017.
    .. [2] NACA-TR-1135. Ames Research Staff, "Equations, Tables, and Charts for
           Compressible Flow", NACA Report 1135, 1953.
    .. [3] NASA-GLENN-MASS-FLOW. NASA Glenn Research Center, isentropic mass-flow
           conservation and density/velocity context.
           https://www1.grc.nasa.gov/beginners-guide-to-aeronautics/mass-flow-rate-equations/
    """

    for value, label in (
        (chamber_pressure_pa, "Chamber pressure"),
        (chamber_temperature_k, "Chamber temperature"),
        (throat_area_m2, "Throat area"),
        (area_ratio, "Area ratio"),
        (gas_constant_j_kg_k, "Specific gas constant"),
    ):
        _positive(value, label)
    if area_ratio < 1.0:
        raise DomainError("Exit-to-throat area ratio must be at least one.")
    if not isfinite(ambient_pressure_pa) or ambient_pressure_pa < 0.0:
        raise DomainError("Ambient pressure must be finite and non-negative.")
    if not isfinite(gamma) or gamma <= 1.0:
        raise DomainError("Gamma must be finite and greater than one.")

    exit_mach = mach_from_isentropic("area_ratio_supersonic", area_ratio, gamma)
    exit_ratios = calculate_isentropic(exit_mach, gamma)
    exit_temperature = chamber_temperature_k * exit_ratios.temperature_total_ratio
    exit_pressure = chamber_pressure_pa * exit_ratios.pressure_total_ratio
    exit_velocity = exit_mach * sqrt(gamma * gas_constant_j_kg_k * exit_temperature)
    mass_flow = (
        chamber_pressure_pa
        * throat_area_m2
        / sqrt(chamber_temperature_k)
        * sqrt(gamma / gas_constant_j_kg_k)
        * (2.0 / (gamma + 1.0)) ** ((gamma + 1.0) / (2.0 * (gamma - 1.0)))
    )
    exit_area = area_ratio * throat_area_m2
    momentum_thrust = mass_flow * exit_velocity
    pressure_thrust = (exit_pressure - ambient_pressure_pa) * exit_area
    thrust = momentum_thrust + pressure_thrust
    if thrust <= 0.0:
        raise DomainError("The ideal nozzle state produces non-positive thrust.")
    characteristic_velocity = chamber_pressure_pa * throat_area_m2 / mass_flow
    if abs(exit_pressure - ambient_pressure_pa) <= 0.01 * max(exit_pressure, 1.0):
        regime = "adapted"
    elif exit_pressure > ambient_pressure_pa:
        regime = "underexpanded"
    else:
        regime = "overexpanded"

    return NozzleResult(
        throat_area_m2=throat_area_m2,
        exit_area_m2=exit_area,
        area_ratio=area_ratio,
        exit_mach=exit_mach,
        exit_temperature_k=exit_temperature,
        exit_pressure_pa=exit_pressure,
        exit_velocity_m_s=exit_velocity,
        mass_flow_kg_s=mass_flow,
        momentum_thrust_n=momentum_thrust,
        pressure_thrust_n=pressure_thrust,
        thrust_n=thrust,
        characteristic_velocity_m_s=characteristic_velocity,
        thrust_coefficient=thrust / (chamber_pressure_pa * throat_area_m2),
        specific_impulse_s=thrust / (mass_flow * STANDARD_GRAVITY_M_S2),
        pressure_regime=regime,
    )

