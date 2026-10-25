"""Closed-system polytropic ideal-gas process paths and energy accounting.

References
----------
.. [1] NASA-GAS. NASA Glenn Research Center, "Equation of State (Ideal Gas)",
       https://www1.grc.nasa.gov/beginners-guide-to-aeronautics/equation-of-state-ideal-gas-2/
.. [2] NASA Glenn Research Center, "Specific Heats",
       https://www.grc.nasa.gov/www/k-12/airplane/specheat.html
"""

from dataclasses import dataclass
from math import isfinite, log

from rocket_propulsion.core.constants import AIR_SPECIFIC_GAS_CONSTANT_J_KG_K
from rocket_propulsion.core.errors import DomainError


@dataclass(frozen=True, slots=True)
class PolytropicPoint:
    """One sampled point along a quasistatic polytropic process path.

    Attributes
    ----------
    fraction : float
        Position along the path in [0, 1]; pressure is sampled
        geometrically, ``p = p1 (p2/p1)^fraction``.
    pressure_pa : float
        Absolute pressure, Pa.
    temperature_k : float
        Absolute temperature, K.
    specific_volume_m3_kg : float
        ``R T / p``, m^3/kg.
    entropy_change_j_kg_k : float
        Specific entropy relative to the inlet state, J/(kg K).

    References
    ----------
    .. [1] NASA-GAS. NASA Glenn Research Center, "Equation of State (Ideal Gas)",
           https://www1.grc.nasa.gov/beginners-guide-to-aeronautics/equation-of-state-ideal-gas-2/
    .. [2] NASA Glenn Research Center, "Specific Heats",
           https://www.grc.nasa.gov/www/k-12/airplane/specheat.html
    """

    fraction: float
    pressure_pa: float
    temperature_k: float
    specific_volume_m3_kg: float
    entropy_change_j_kg_k: float


@dataclass(frozen=True, slots=True)
class PolytropicResult:
    """End state, first-law terms and sampled path for ``p v^n = const``.

    Attributes
    ----------
    process : {"compression", "expansion"}
        Direction implied by the pressure ratio.
    exponent : float
        Polytropic exponent ``n``.
    inlet_pressure_pa, outlet_pressure_pa : float
        End pressures, Pa.
    inlet_temperature_k, outlet_temperature_k : float
        End temperatures, K.
    inlet_specific_volume_m3_kg, outlet_specific_volume_m3_kg : float
        End specific volumes, m^3/kg.
    boundary_work_by_gas_j_kg : float
        Specific boundary work done by the gas, J/kg.
    internal_energy_change_j_kg, enthalpy_change_j_kg : float
        ``cv (T2 - T1)`` and ``cp (T2 - T1)``, J/kg.
    heat_into_gas_j_kg : float
        Specific heat added to the gas from the first law, J/kg.
    entropy_change_j_kg_k : float
        Specific entropy change, J/(kg K).
    path : tuple of PolytropicPoint
        Sampled path from inlet to outlet.

    References
    ----------
    .. [1] NASA-GAS. NASA Glenn Research Center, "Equation of State (Ideal Gas)",
           https://www1.grc.nasa.gov/beginners-guide-to-aeronautics/equation-of-state-ideal-gas-2/
    .. [2] NASA Glenn Research Center, "Specific Heats",
           https://www.grc.nasa.gov/www/k-12/airplane/specheat.html
    """

    process: str
    exponent: float
    inlet_pressure_pa: float
    outlet_pressure_pa: float
    inlet_temperature_k: float
    outlet_temperature_k: float
    inlet_specific_volume_m3_kg: float
    outlet_specific_volume_m3_kg: float
    boundary_work_by_gas_j_kg: float
    internal_energy_change_j_kg: float
    enthalpy_change_j_kg: float
    heat_into_gas_j_kg: float
    entropy_change_j_kg_k: float
    path: tuple[PolytropicPoint, ...]


def calculate_polytropic_process(
    *,
    inlet_pressure_pa: float,
    outlet_pressure_pa: float,
    inlet_temperature_k: float,
    exponent: float,
    gamma: float = 1.4,
    gas_constant_j_kg_k: float = AIR_SPECIFIC_GAS_CONSTANT_J_KG_K,
    point_count: int = 61,
) -> PolytropicResult:
    """Calculate a quasistatic closed-system ideal-gas polytropic process.

    Parameters
    ----------
    inlet_pressure_pa, outlet_pressure_pa : float
        Absolute end pressures, Pa; they must differ.
    inlet_temperature_k : float
        Absolute inlet temperature, K.
    exponent : float
        Polytropic exponent ``n``; must be finite and positive.
    gamma : float, optional
        Ratio of specific heats; must exceed one.
    gas_constant_j_kg_k : float, optional
        Specific gas constant, J/(kg K). Defaults to the registered dry-air
        value ``AIR_SPECIFIC_GAS_CONSTANT_J_KG_K``.
    point_count : int, optional
        Number of sampled path points, between 21 and 401 inclusive.

    Returns
    -------
    PolytropicResult
        End states, first-law terms and the sampled path.

    Raises
    ------
    DomainError
        If an input is non-finite or non-positive, the pressures are equal,
        ``gamma <= 1`` or ``point_count`` is out of range.

    Notes
    -----
    The path follows ``p v^n = constant`` with the ideal-gas equation of
    state ``p v = R T`` [1]_, so ``T2 = T1 (p2/p1)^((n-1)/n)``. Boundary work
    by the gas is ``(p2 v2 - p1 v1)/(1 - n)``; for ``abs(n - 1) < 1e-10`` the
    isothermal limit ``R T1 ln(v2/v1)`` is used. Heat follows the first law
    with the sign convention ``Q - W_by = delta U``, using ``cv = R/(gamma-1)``
    and ``cp = gamma cv`` [2]_. The entropy change is
    ``cp ln(T2/T1) - R ln(p2/p1)``. The model connects thermodynamic diagrams
    with energy terms; it does not replace component-efficiency models for
    compressors or turbines.

    References
    ----------
    .. [1] NASA-GAS. NASA Glenn Research Center, "Equation of State (Ideal Gas)",
           https://www1.grc.nasa.gov/beginners-guide-to-aeronautics/equation-of-state-ideal-gas-2/
    .. [2] NASA Glenn Research Center, "Specific Heats",
           https://www.grc.nasa.gov/www/k-12/airplane/specheat.html
    """

    for value, label in (
        (inlet_pressure_pa, "Inlet pressure"),
        (outlet_pressure_pa, "Outlet pressure"),
        (inlet_temperature_k, "Inlet temperature"),
        (exponent, "Polytropic exponent"),
        (gas_constant_j_kg_k, "Specific gas constant"),
    ):
        if not isfinite(value) or value <= 0.0:
            raise DomainError(f"{label} must be finite and greater than zero.")
    if inlet_pressure_pa == outlet_pressure_pa:
        raise DomainError("Inlet and outlet pressures must differ.")
    if gamma <= 1.0:
        raise DomainError("Gamma must be greater than one.")
    if not 21 <= point_count <= 401:
        raise DomainError("Point count must be between 21 and 401.")

    pressure_ratio = outlet_pressure_pa / inlet_pressure_pa
    outlet_temperature = inlet_temperature_k * pressure_ratio ** ((exponent - 1.0) / exponent)
    inlet_volume = gas_constant_j_kg_k * inlet_temperature_k / inlet_pressure_pa
    outlet_volume = gas_constant_j_kg_k * outlet_temperature / outlet_pressure_pa
    if abs(exponent - 1.0) < 1.0e-10:
        boundary_work = gas_constant_j_kg_k * inlet_temperature_k * log(
            outlet_volume / inlet_volume
        )
    else:
        boundary_work = (
            outlet_pressure_pa * outlet_volume - inlet_pressure_pa * inlet_volume
        ) / (1.0 - exponent)

    cv = gas_constant_j_kg_k / (gamma - 1.0)
    cp = gamma * cv
    delta_u = cv * (outlet_temperature - inlet_temperature_k)
    delta_h = cp * (outlet_temperature - inlet_temperature_k)
    delta_s = cp * log(outlet_temperature / inlet_temperature_k) - gas_constant_j_kg_k * log(
        pressure_ratio
    )
    heat = delta_u + boundary_work

    path: list[PolytropicPoint] = []
    for index in range(point_count):
        fraction = index / (point_count - 1)
        pressure = inlet_pressure_pa * pressure_ratio**fraction
        temperature = inlet_temperature_k * (pressure / inlet_pressure_pa) ** (
            (exponent - 1.0) / exponent
        )
        volume = gas_constant_j_kg_k * temperature / pressure
        entropy = cp * log(temperature / inlet_temperature_k) - gas_constant_j_kg_k * log(
            pressure / inlet_pressure_pa
        )
        path.append(
            PolytropicPoint(
                fraction=fraction,
                pressure_pa=pressure,
                temperature_k=temperature,
                specific_volume_m3_kg=volume,
                entropy_change_j_kg_k=entropy,
            )
        )

    return PolytropicResult(
        process="compression" if pressure_ratio > 1.0 else "expansion",
        exponent=exponent,
        inlet_pressure_pa=inlet_pressure_pa,
        outlet_pressure_pa=outlet_pressure_pa,
        inlet_temperature_k=inlet_temperature_k,
        outlet_temperature_k=outlet_temperature,
        inlet_specific_volume_m3_kg=inlet_volume,
        outlet_specific_volume_m3_kg=outlet_volume,
        boundary_work_by_gas_j_kg=boundary_work,
        internal_energy_change_j_kg=delta_u,
        enthalpy_change_j_kg=delta_h,
        heat_into_gas_j_kg=heat,
        entropy_change_j_kg_k=delta_s,
        path=tuple(path),
    )

