"""Preliminary propellant storage and isothermal pressurant inventories.

References
----------
.. [1] NASA-GAS. NASA Glenn Research Center, "Equation of State (Ideal Gas)";
       density=m/V and ideal-gas equation.
.. [2] NASA-SP8112. NASA, "Pressurization Systems for Liquid Rockets"; preliminary
       blowdown/regulated pressurization context.
"""

from __future__ import annotations

from fractions import Fraction

from rocket_propulsion.core.errors import DomainError

from .evidence import AnalysisResult, Citation, CitedScalar, derived, quantity
from .rational_numerics import binary_ledger


def liquid_tank_volume(propellant_mass: CitedScalar, density: CitedScalar,
                       ullage_fraction: CitedScalar) -> AnalysisResult:
    """Size a liquid propellant tank from mass, density and ullage fraction.

    Parameters
    ----------
    propellant_mass : CitedScalar
        Stored liquid mass, kg; non-negative.
    density : CitedScalar
        Liquid density at the design storage condition, kg/m^3; positive. The value
        must carry its own temperature, pressure and source.
    ullage_fraction : CitedScalar
        Ullage as a fraction of TOTAL tank volume, in [0, 1); not a percentage added
        to the liquid volume.

    Returns
    -------
    AnalysisResult
        ``liquid_volume``, ``tank_volume`` and ``ullage_volume``, m^3.

    Raises
    ------
    DomainError
        If a quantity is out of domain or an output is not representable.

    Notes
    -----
    V_liquid = m / rho [1]_ and V_tank = V_liquid / (1 - ullage_fraction); the
    geometric ullage rearrangement is stated in this module [2]_. Wall, liner, insulation,
    propellant management device and structural mass are not estimated.

    References
    ----------
    .. [1] NASA-GAS. NASA Glenn Research Center, "Equation of State (Ideal Gas)";
           density is mass divided by volume.
    .. [2] NASA-SP8112. NASA, "Pressurization Systems for Liquid Rockets"; ullage and
           stored-liquid tank pressurization.
    """
    m_value = quantity(propellant_mass, "kg", minimum=0)
    rho_value = quantity(density, "kg/m^3", strictly_positive=True)
    ullage_value = quantity(ullage_fraction, "1", minimum=0)
    if ullage_value >= 1:
        raise DomainError("Ullage fraction must be less than one.")
    m, rho, ullage = map(Fraction, (m_value, rho_value, ullage_value))
    liquid = m / rho
    tank = liquid / (1 - ullage)
    inputs = (propellant_mass, density, ullage_fraction)
    methods = (Citation("NASA-GAS", "Density=m/V; liquid density is supplied separately"),
               Citation("NASA-SP8112", "Ullage bookkeeping; geometric rearrangement stated here"))
    values = binary_ledger((("liquid_volume", liquid, "m^3"), ("tank_volume", tank, "m^3"),
                            ("ullage_volume", tank * ullage, "m^3")), "Liquid storage")
    return AnalysisResult("constant_density_liquid_storage", tuple(
        (name, derived(value, unit, inputs, methods)) for name, value, unit in values
    ), ("No vessel wall, liner, insulation, PMD or structural qualification mass estimate.",
        "No density extrapolation: supplied density must apply to the design storage condition.",
        "Exact rational input arithmetic and final binary rounding; physical uncertainty is not quantified."))


def ideal_cold_gas_storage(usable_mass: CitedScalar, initial_pressure: CitedScalar,
                           cutoff_pressure: CitedScalar, temperature: CitedScalar,
                           specific_gas_constant: CitedScalar) -> AnalysisResult:
    """Size a rigid isothermal ideal-gas tank from its usable pressure interval.

    Parameters
    ----------
    usable_mass : CitedScalar
        Gas mass expelled between the initial and cutoff pressures, kg.
    initial_pressure : CitedScalar
        Initial absolute tank pressure p0, Pa; positive.
    cutoff_pressure : CitedScalar
        Absolute pressure at which blowdown ends, pf, Pa; positive and below p0.
    temperature : CitedScalar
        Absolute gas temperature, K; positive.
    specific_gas_constant : CitedScalar
        Specific gas constant R of the gas, J/(kg*K); positive.

    Returns
    -------
    AnalysisResult
        ``tank_volume`` (m^3), ``loaded_gas`` and ``remaining_gas`` (kg).

    Raises
    ------
    DomainError
        If a quantity is out of domain, p0 <= pf, or an output is not representable.

    Notes
    -----
    V = m_usable * R * T / (p0 - pf), m_loaded = p0 * V / (R * T) and
    m_remaining = pf * V / (R * T), from p * V = m * R * T at fixed V and T [1]_ and
    the blowdown cutoff inventory [2]_. This is an analytical ideal-gas model, not an
    extrapolation of a measured density or thruster map. Real-gas and thermal effects
    at high pressure are not validated, and unregulated thrust and Isp generally vary
    with pressure.

    References
    ----------
    .. [1] NASA-GAS. NASA Glenn Research Center, "Equation of State (Ideal Gas)"; p=rho
           R T, absolute pressure/temperature, ideal-gas domain.
    .. [2] NASA-SP8112. NASA, "Pressurization Systems for Liquid Rockets"; blowdown
           inventory and cutoff boundary.
    """
    m_value = quantity(usable_mass, "kg", minimum=0)
    p0_value = quantity(initial_pressure, "Pa", strictly_positive=True)
    pf_value = quantity(cutoff_pressure, "Pa", strictly_positive=True)
    t_value = quantity(temperature, "K", strictly_positive=True)
    gas_r_value = quantity(specific_gas_constant, "J/(kg*K)", strictly_positive=True)
    if p0_value <= pf_value:
        raise DomainError("Initial absolute pressure must exceed cutoff pressure.")
    m, p0, pf, t, gas_r = map(Fraction, (m_value, p0_value, pf_value, t_value, gas_r_value))
    thermal = gas_r * t
    volume = m * thermal / (p0 - pf)
    loaded = m * p0 / (p0 - pf)
    remaining = m * pf / (p0 - pf)
    inputs = (usable_mass, initial_pressure, cutoff_pressure, temperature, specific_gas_constant)
    methods = (Citation("NASA-GAS", "pV=mRT; subtract endpoint inventories at fixed V,T"),
               Citation("NASA-SP8112", "Preliminary blowdown cutoff and remaining gas"))
    values = binary_ledger((("tank_volume", volume, "m^3"), ("loaded_gas", loaded, "kg"),
                            ("remaining_gas", remaining, "kg")), "Cold gas storage")
    return AnalysisResult("isothermal_ideal_gas_blowdown_storage", tuple(
        (name, derived(value, unit, inputs, methods)) for name, value, unit in values
    ), ("Ideal gas and isothermal only; high-pressure real-gas and thermal effects unvalidated.",
        "Unregulated thrust/Isp generally varies with pressure; this is only an inventory model.",
        "No pressure-vessel hardware mass or pressure qualification is inferred.",
        "Exact rational input arithmetic and final binary rounding; physical uncertainty is not quantified."))


def ideal_regulated_pressurant(displaced_volume: CitedScalar, tank_pressure: CitedScalar,
                               bottle_initial_pressure: CitedScalar,
                               bottle_final_pressure: CitedScalar,
                               temperature: CitedScalar,
                               specific_gas_constant: CitedScalar) -> AnalysisResult:
    """Size an ideal isothermal pressurant bottle for a displaced liquid volume.

    Parameters
    ----------
    displaced_volume : CitedScalar
        Liquid volume to be displaced at regulated pressure, m^3; non-negative.
    tank_pressure : CitedScalar
        Regulated absolute propellant-tank pressure p_t, Pa; positive.
    bottle_initial_pressure : CitedScalar
        Initial absolute bottle pressure p_b0, Pa.
    bottle_final_pressure : CitedScalar
        Final absolute bottle pressure p_bf, Pa; must satisfy p_b0 > p_bf > p_t.
    temperature : CitedScalar
        Absolute gas temperature, K; positive.
    specific_gas_constant : CitedScalar
        Specific gas constant R of the gas, J/(kg*K); positive.

    Returns
    -------
    AnalysisResult
        ``bottle_volume`` (m^3), ``loaded_pressurant``, ``delivered_pressurant`` and
        ``remaining_pressurant`` (kg).

    Raises
    ------
    DomainError
        If a quantity is out of domain, the pressures are not ordered, or an output
        is not representable.

    Notes
    -----
    V_b = p_t * V_displaced / (p_b0 - p_bf) from the ideal-gas law and conservation of
    the transferred gas [1]_ in a regulated bottle-to-ullage architecture [2]_. Bottle
    and receiving ullage share the supplied temperature; initial tank gas, lines and
    regulator losses are excluded, and transient dynamics need a separate model.

    References
    ----------
    .. [1] NASA-GAS. NASA Glenn Research Center, "Equation of State (Ideal Gas)"; pV=mRT
           and conservation of transferred gas.
    .. [2] NASA-SP8112. NASA, "Pressurization Systems for Liquid Rockets"; regulated
           bottle-to-ullage pressurization architecture.
    """
    displaced_value = quantity(displaced_volume, "m^3", minimum=0)
    pt_value = quantity(tank_pressure, "Pa", strictly_positive=True)
    p0_value = quantity(bottle_initial_pressure, "Pa", strictly_positive=True)
    pf_value = quantity(bottle_final_pressure, "Pa", strictly_positive=True)
    t_value = quantity(temperature, "K", strictly_positive=True)
    gas_r_value = quantity(specific_gas_constant, "J/(kg*K)", strictly_positive=True)
    if not p0_value > pf_value > pt_value:
        raise DomainError("Require bottle initial > bottle final > regulated tank pressure.")
    displaced, pt, p0, pf, t, gas_r = map(Fraction, (
        displaced_value, pt_value, p0_value, pf_value, t_value, gas_r_value))
    thermal = gas_r * t
    volume = pt * displaced / (p0 - pf)
    used = pt * displaced / thermal
    loaded = used * p0 / (p0 - pf)
    remaining = used * pf / (p0 - pf)
    inputs = (displaced_volume, tank_pressure, bottle_initial_pressure,
              bottle_final_pressure, temperature, specific_gas_constant)
    methods = (Citation("NASA-GAS", "pV=mRT and bottle-to-ullage mass balance"),
               Citation("NASA-SP8112", "Preliminary regulated pressurization"))
    values = binary_ledger((("bottle_volume", volume, "m^3"),
                            ("loaded_pressurant", loaded, "kg"),
                            ("delivered_pressurant", used, "kg"),
                            ("remaining_pressurant", remaining, "kg")), "Regulated pressurant")
    return AnalysisResult("isothermal_ideal_regulated_pressurant_inventory", tuple(
        (name, derived(value, unit, inputs, methods)) for name, value, unit in values
    ), ("Same constant gas temperature in bottle and ullage; no real-gas/heat-transfer solution.",
        "Additional initial tank gas, regulator margins and non-expelled pressurant need budgeting.",
        "Transient dynamics require a separate model with independently validated parameters.",
        "Exact rational input arithmetic and final binary rounding; physical uncertainty is not quantified."))
