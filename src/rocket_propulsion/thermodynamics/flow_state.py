"""Connections between static, stagnation and process thermodynamic states.

All relations assume a calorically perfect gas (constant ``cp``, ``cv`` and
``gamma``).

References
----------
.. [1] NACA-TR-1135. Ames Research Staff, "Equations, Tables, and Charts for
       Compressible Flow", NACA Report 1135, 1953,
       https://ntrs.nasa.gov/citations/19930091059
.. [2] NASA-ENGINESIM-THEORY. NASA Glenn Research Center, ideal gas turbine
       temperature/pressure and shaft-work relations,
       https://www.grc.nasa.gov/www/k-12/airplane/EngineTheory.pdf
"""

from dataclasses import dataclass
from math import isfinite, log, sqrt

from rocket_propulsion.core.constants import AIR_SPECIFIC_GAS_CONSTANT_J_KG_K
from rocket_propulsion.core.errors import DomainError


@dataclass(frozen=True, slots=True)
class StagnationState:
    """Static and stagnation properties of adiabatic perfect-gas flow.

    Attributes
    ----------
    mach : float
        Local Mach number.
    static_pressure_pa, stagnation_pressure_pa : float
        Static and isentropic stagnation pressure, Pa.
    static_temperature_k, stagnation_temperature_k : float
        Static and stagnation temperature, K.
    static_density_kg_m3, stagnation_density_kg_m3 : float
        Static and isentropic stagnation density, kg/m^3.
    speed_of_sound_m_s : float
        Static speed of sound, m/s.
    velocity_m_s : float
        Flow speed ``M a``, m/s.
    dynamic_pressure_pa : float
        Conventional dynamic pressure ``rho V^2 / 2``, Pa.
    static_enthalpy_j_kg, stagnation_enthalpy_j_kg : float
        ``cp T`` and ``cp T0`` referenced to 0 K, J/kg.

    References
    ----------
    .. [1] NACA-TR-1135. Ames Research Staff, "Equations, Tables, and Charts for
           Compressible Flow", NACA Report 1135, 1953,
           https://ntrs.nasa.gov/citations/19930091059
    """

    mach: float
    static_pressure_pa: float
    static_temperature_k: float
    static_density_kg_m3: float
    stagnation_pressure_pa: float
    stagnation_temperature_k: float
    stagnation_density_kg_m3: float
    speed_of_sound_m_s: float
    velocity_m_s: float
    dynamic_pressure_pa: float
    static_enthalpy_j_kg: float
    stagnation_enthalpy_j_kg: float


@dataclass(frozen=True, slots=True)
class IsentropicProcessResult:
    """Ideal and efficiency-corrected compression or expansion state.

    Attributes
    ----------
    process : {"compression", "expansion"}
        Direction implied by the pressure ratio.
    inlet_pressure_pa, outlet_pressure_pa : float
        End pressures, Pa.
    inlet_temperature_k : float
        Inlet temperature, K.
    ideal_outlet_temperature_k : float
        Isentropic outlet temperature, K.
    actual_outlet_temperature_k : float
        Outlet temperature after applying the efficiency, K.
    efficiency : float
        Isentropic efficiency in (0, 1].
    specific_work_j_kg : float
        ``cp (T2 - T1)``; positive for compression (work input), negative for
        expansion, J/kg.
    entropy_change_j_kg_k : float
        Specific entropy change between inlet and actual outlet, J/(kg K).

    References
    ----------
    .. [1] NACA-TR-1135. Ames Research Staff, "Equations, Tables, and Charts for
           Compressible Flow", NACA Report 1135, 1953,
           https://ntrs.nasa.gov/citations/19930091059
    .. [2] NASA-ENGINESIM-THEORY. NASA Glenn Research Center, ideal gas turbine
           temperature/pressure and shaft-work relations,
           https://www.grc.nasa.gov/www/k-12/airplane/EngineTheory.pdf
    """

    process: str
    inlet_pressure_pa: float
    outlet_pressure_pa: float
    inlet_temperature_k: float
    ideal_outlet_temperature_k: float
    actual_outlet_temperature_k: float
    efficiency: float
    specific_work_j_kg: float
    entropy_change_j_kg_k: float


def _positive(value: float, label: str) -> None:
    """Raise :class:`DomainError` unless ``value`` is finite and positive."""
    if not isfinite(value) or value <= 0.0:
        raise DomainError(f"{label} must be finite and greater than zero.")


def calculate_stagnation_state(
    *,
    mach: float,
    static_pressure_pa: float,
    static_temperature_k: float,
    gamma: float = 1.4,
    gas_constant_j_kg_k: float = AIR_SPECIFIC_GAS_CONSTANT_J_KG_K,
) -> StagnationState:
    """Connect a static perfect-gas flow state to its stagnation state.

    Parameters
    ----------
    mach : float
        Local Mach number; must be finite and positive.
    static_pressure_pa : float
        Absolute static pressure, Pa.
    static_temperature_k : float
        Absolute static temperature, K.
    gamma : float, optional
        Ratio of specific heats; must exceed one.
    gas_constant_j_kg_k : float, optional
        Specific gas constant, J/(kg K). Defaults to the registered dry-air
        value ``AIR_SPECIFIC_GAS_CONSTANT_J_KG_K``.

    Returns
    -------
    StagnationState
        Static and stagnation properties.

    Raises
    ------
    DomainError
        If any input is non-finite or non-positive, or ``gamma <= 1``.

    Notes
    -----
    With ``f = 1 + (gamma - 1) M^2 / 2`` the adiabatic, reversible
    relations [1]_ are ``T0 = T f``, ``p0 = p f^(gamma/(gamma-1))`` and
    ``rho0 = rho f^(1/(gamma-1))``, consistent with ``h0 = h + V^2/2``. The
    dynamic pressure is the conventional ``rho V^2 / 2``; at finite Mach
    number it is not equal to ``p0 - p``.

    References
    ----------
    .. [1] NACA-TR-1135. Ames Research Staff, "Equations, Tables, and Charts for
           Compressible Flow", NACA Report 1135, 1953,
           https://ntrs.nasa.gov/citations/19930091059
    """

    for value, label in (
        (mach, "Mach number"),
        (static_pressure_pa, "Static pressure"),
        (static_temperature_k, "Static temperature"),
        (gas_constant_j_kg_k, "Specific gas constant"),
    ):
        _positive(value, label)
    if not isfinite(gamma) or gamma <= 1.0:
        raise DomainError("Gamma must be finite and greater than one.")

    factor = 1.0 + 0.5 * (gamma - 1.0) * mach**2
    static_density = static_pressure_pa / (gas_constant_j_kg_k * static_temperature_k)
    stagnation_temperature = static_temperature_k * factor
    stagnation_pressure = static_pressure_pa * factor ** (gamma / (gamma - 1.0))
    stagnation_density = static_density * factor ** (1.0 / (gamma - 1.0))
    speed_of_sound = sqrt(gamma * gas_constant_j_kg_k * static_temperature_k)
    velocity = mach * speed_of_sound
    cp = gamma * gas_constant_j_kg_k / (gamma - 1.0)
    return StagnationState(
        mach=mach,
        static_pressure_pa=static_pressure_pa,
        static_temperature_k=static_temperature_k,
        static_density_kg_m3=static_density,
        stagnation_pressure_pa=stagnation_pressure,
        stagnation_temperature_k=stagnation_temperature,
        stagnation_density_kg_m3=stagnation_density,
        speed_of_sound_m_s=speed_of_sound,
        velocity_m_s=velocity,
        dynamic_pressure_pa=0.5 * static_density * velocity**2,
        static_enthalpy_j_kg=cp * static_temperature_k,
        stagnation_enthalpy_j_kg=cp * stagnation_temperature,
    )


def calculate_isentropic_process(
    *,
    inlet_pressure_pa: float,
    outlet_pressure_pa: float,
    inlet_temperature_k: float,
    efficiency: float = 1.0,
    gamma: float = 1.4,
    gas_constant_j_kg_k: float = AIR_SPECIFIC_GAS_CONSTANT_J_KG_K,
) -> IsentropicProcessResult:
    """Calculate ideal and efficiency-corrected compression or expansion.

    Parameters
    ----------
    inlet_pressure_pa, outlet_pressure_pa : float
        Absolute end pressures, Pa; they must differ. A higher outlet
        pressure means compression.
    inlet_temperature_k : float
        Absolute inlet temperature, K.
    efficiency : float, optional
        Isentropic efficiency in (0, 1].
    gamma : float, optional
        Ratio of specific heats; must exceed one.
    gas_constant_j_kg_k : float, optional
        Specific gas constant, J/(kg K). Defaults to the registered dry-air
        value ``AIR_SPECIFIC_GAS_CONSTANT_J_KG_K``.

    Returns
    -------
    IsentropicProcessResult
        Ideal and actual outlet temperatures, specific work and entropy
        change.

    Raises
    ------
    DomainError
        If an input is non-finite or non-positive, ``gamma <= 1``, the
        efficiency is outside (0, 1], or the pressures are equal.

    Notes
    -----
    The isentropic outlet temperature is ``T2s = T1 (p2/p1)^((gamma-1)/gamma)``
    [1]_. Compressor efficiency is ``(T2s - T1)/(T2 - T1)`` and
    turbine/nozzle efficiency is ``(T1 - T2)/(T1 - T2s)`` [2]_. Specific work
    ``cp (T2 - T1)`` is positive for compression and negative for expansion
    (work-input sign convention). The entropy change is
    ``cp ln(T2/T1) - R ln(p2/p1)``. The component is adiabatic and steady;
    heat loss and kinetic-energy changes are not modelled.

    References
    ----------
    .. [1] NACA-TR-1135. Ames Research Staff, "Equations, Tables, and Charts for
           Compressible Flow", NACA Report 1135, 1953,
           https://ntrs.nasa.gov/citations/19930091059
    .. [2] NASA-ENGINESIM-THEORY. NASA Glenn Research Center, ideal gas turbine
           temperature/pressure and shaft-work relations,
           https://www.grc.nasa.gov/www/k-12/airplane/EngineTheory.pdf
    """

    for value, label in (
        (inlet_pressure_pa, "Inlet pressure"),
        (outlet_pressure_pa, "Outlet pressure"),
        (inlet_temperature_k, "Inlet temperature"),
        (gas_constant_j_kg_k, "Specific gas constant"),
    ):
        _positive(value, label)
    if not isfinite(gamma) or gamma <= 1.0:
        raise DomainError("Gamma must be finite and greater than one.")
    if not isfinite(efficiency) or not 0.0 < efficiency <= 1.0:
        raise DomainError("Efficiency must be in the interval (0, 1].")
    if outlet_pressure_pa == inlet_pressure_pa:
        raise DomainError("Inlet and outlet pressures must differ.")

    pressure_ratio = outlet_pressure_pa / inlet_pressure_pa
    ideal_temperature = inlet_temperature_k * pressure_ratio ** ((gamma - 1.0) / gamma)
    if pressure_ratio > 1.0:
        process = "compression"
        actual_temperature = inlet_temperature_k + (
            ideal_temperature - inlet_temperature_k
        ) / efficiency
    else:
        process = "expansion"
        actual_temperature = inlet_temperature_k - efficiency * (
            inlet_temperature_k - ideal_temperature
        )

    cp = gamma * gas_constant_j_kg_k / (gamma - 1.0)
    return IsentropicProcessResult(
        process=process,
        inlet_pressure_pa=inlet_pressure_pa,
        outlet_pressure_pa=outlet_pressure_pa,
        inlet_temperature_k=inlet_temperature_k,
        ideal_outlet_temperature_k=ideal_temperature,
        actual_outlet_temperature_k=actual_temperature,
        efficiency=efficiency,
        specific_work_j_kg=cp * (actual_temperature - inlet_temperature_k),
        entropy_change_j_kg_k=cp * log(actual_temperature / inlet_temperature_k)
        - gas_constant_j_kg_k * log(pressure_ratio),
    )

