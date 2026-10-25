"""Calorically perfect ideal-gas state relations.

References
----------
.. [1] NASA-GAS. NASA Glenn Research Center, "Equation of State (Ideal Gas)",
       https://www1.grc.nasa.gov/beginners-guide-to-aeronautics/equation-of-state-ideal-gas-2/
.. [2] NASA Glenn Research Center, "Specific Heats",
       https://www.grc.nasa.gov/www/k-12/airplane/specheat.html
"""

from dataclasses import dataclass
from math import isfinite, sqrt

from rocket_propulsion.core.constants import AIR_SPECIFIC_GAS_CONSTANT_J_KG_K
from rocket_propulsion.core.errors import DomainError


@dataclass(frozen=True, slots=True)
class IdealGasState:
    """Complete SI state of a calorically perfect ideal gas.

    Attributes
    ----------
    pressure_pa : float
        Absolute static pressure, Pa.
    temperature_k : float
        Absolute static temperature, K.
    density_kg_m3 : float
        Static density, kg/m^3.
    gas_constant_j_kg_k : float
        Specific gas constant R, J/(kg K).
    gamma : float
        Ratio of specific heats cp/cv.
    cp_j_kg_k, cv_j_kg_k : float
        Constant-pressure and constant-volume specific heats, J/(kg K).
    enthalpy_j_kg : float
        Sensible enthalpy ``cp T`` referenced to 0 K, J/kg.
    internal_energy_j_kg : float
        Sensible internal energy ``cv T`` referenced to 0 K, J/kg.
    speed_of_sound_m_s : float
        ``sqrt(gamma R T)``, m/s.

    References
    ----------
    .. [1] NASA-GAS. NASA Glenn Research Center, "Equation of State (Ideal Gas)",
           https://www1.grc.nasa.gov/beginners-guide-to-aeronautics/equation-of-state-ideal-gas-2/
    .. [2] NASA Glenn Research Center, "Specific Heats",
           https://www.grc.nasa.gov/www/k-12/airplane/specheat.html
    .. [3] NACA-TR-1135. Ames Research Staff, "Equations, Tables, and Charts for
           Compressible Flow", NACA Report 1135, 1953,
           https://ntrs.nasa.gov/citations/19930091059
    """

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
    gas_constant_j_kg_k: float = AIR_SPECIFIC_GAS_CONSTANT_J_KG_K,
    gamma: float = 1.4,
) -> IdealGasState:
    """Complete an ideal-gas state from exactly two of ``p``, ``T`` and ``rho``.

    Parameters
    ----------
    pressure_pa : float, optional
        Absolute static pressure, Pa.
    temperature_k : float, optional
        Absolute static temperature, K.
    density_kg_m3 : float, optional
        Static density, kg/m^3.
    gas_constant_j_kg_k : float, optional
        Specific gas constant, J/(kg K). Defaults to the registered dry-air
        value ``AIR_SPECIFIC_GAS_CONSTANT_J_KG_K``.
    gamma : float, optional
        Ratio of specific heats; must exceed one.

    Returns
    -------
    IdealGasState
        Complete immutable state in SI units.

    Raises
    ------
    DomainError
        If the gas constant or ``gamma`` is out of range, if not exactly two
        state properties are supplied, or if a supplied property is not
        finite and positive.

    Notes
    -----
    The missing property follows from the thermal equation of state
    ``p = rho R T`` [1]_. Calorically perfect heat capacities use
    ``cp - cv = R`` and ``gamma = cp/cv`` [2]_, giving ``cv = R/(gamma - 1)``
    and ``cp = gamma cv``. Enthalpy and internal energy are sensible values
    referenced to 0 K, and the speed of sound is ``sqrt(gamma R T)`` [3]_. Real-gas
    effects, dissociation and temperature-dependent heat capacities are not
    modelled.

    References
    ----------
    .. [1] NASA-GAS. NASA Glenn Research Center, "Equation of State (Ideal Gas)",
           https://www1.grc.nasa.gov/beginners-guide-to-aeronautics/equation-of-state-ideal-gas-2/
    .. [2] NASA Glenn Research Center, "Specific Heats",
           https://www.grc.nasa.gov/www/k-12/airplane/specheat.html
    .. [3] NACA-TR-1135. Ames Research Staff, "Equations, Tables, and Charts for
           Compressible Flow", NACA Report 1135, 1953,
           https://ntrs.nasa.gov/citations/19930091059
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

