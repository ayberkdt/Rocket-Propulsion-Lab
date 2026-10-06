"""Calorically perfect ideal-gas state relations."""

from dataclasses import dataclass
from math import isfinite, sqrt

from rocket_propulsion.core.errors import DomainError


@dataclass(frozen=True, slots=True)
class IdealGasState:
    """A complete SI-unit state for a calorically perfect ideal gas."""

    pressure_pa: float
    temperature_k: float
    density_kg_m3: float
    gas_constant_j_kg_k: float
    gamma: float
    cp_j_kg_k: float
    cv_j_kg_k: float
    enthalpy_j_kg: float
    internal_energy_j_kg: float
    speed_of_sound_m_s: float


def complete_ideal_gas_state(
    *,
    pressure_pa: float | None = None,
    temperature_k: float | None = None,
    density_kg_m3: float | None = None,
    gas_constant_j_kg_k: float = 287.05,
    gamma: float = 1.4,
) -> IdealGasState:
    """Complete an ideal-gas state from exactly two of ``p``, ``T``, and ``rho``.

    The missing state property follows from ``p = rho R T``. Calorically perfect
    heat capacities use ``cv = R/(gamma-1)`` and ``cp = gamma*cv``.

    Args:
        pressure_pa: Absolute static pressure in pascals.
        temperature_k: Absolute static temperature in kelvin.
        density_kg_m3: Static density in kilograms per cubic metre.
        gas_constant_j_kg_k: Specific gas constant in J/(kg K).
        gamma: Ratio of specific heats.

    Returns:
        A complete immutable thermodynamic state in SI units.
    """

    if not isfinite(gas_constant_j_kg_k) or gas_constant_j_kg_k <= 0.0:
        raise DomainError("Specific gas constant must be finite and positive.")
    if not isfinite(gamma) or gamma <= 1.0:
        raise DomainError("Gamma must be finite and greater than 1.")

    supplied = [pressure_pa is not None, temperature_k is not None, density_kg_m3 is not None]
    if sum(supplied) != 2:
        raise DomainError("Provide exactly two of pressure, temperature, and density.")

    for label, candidate in (
        ("Pressure", pressure_pa),
        ("Temperature", temperature_k),
        ("Density", density_kg_m3),
    ):
        if candidate is not None and (not isfinite(candidate) or candidate <= 0.0):
            raise DomainError(f"{label} must be finite and positive.")

    if pressure_pa is None:
        assert temperature_k is not None and density_kg_m3 is not None
        pressure_pa = density_kg_m3 * gas_constant_j_kg_k * temperature_k
    elif temperature_k is None:
        assert density_kg_m3 is not None
        temperature_k = pressure_pa / (density_kg_m3 * gas_constant_j_kg_k)
    else:
        density_kg_m3 = pressure_pa / (gas_constant_j_kg_k * temperature_k)

    assert pressure_pa is not None and temperature_k is not None and density_kg_m3 is not None
    cv = gas_constant_j_kg_k / (gamma - 1.0)
    cp = gamma * cv
    return IdealGasState(
        pressure_pa=pressure_pa,
        temperature_k=temperature_k,
        density_kg_m3=density_kg_m3,
        gas_constant_j_kg_k=gas_constant_j_kg_k,
        gamma=gamma,
        cp_j_kg_k=cp,
        cv_j_kg_k=cv,
        enthalpy_j_kg=cp * temperature_k,
        internal_energy_j_kg=cv * temperature_k,
        speed_of_sound_m_s=sqrt(gamma * gas_constant_j_kg_k * temperature_k),
    )

