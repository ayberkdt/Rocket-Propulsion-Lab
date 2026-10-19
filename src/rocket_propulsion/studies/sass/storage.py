"""Preliminary propellant storage and isothermal pressurant inventories.

References
----------
NASA-GAS: density=m/V and ideal-gas equation.
NASA-SP8112: preliminary blowdown/regulated pressurization context.
"""

from __future__ import annotations

from fractions import Fraction

from rocket_propulsion.core.errors import DomainError

from .evidence import AnalysisResult, Citation, CitedScalar, derived, quantity
from .rational_numerics import binary_ledger


def liquid_tank_volume(propellant_mass: CitedScalar, density: CitedScalar,
                       ullage_fraction: CitedScalar) -> AnalysisResult:
    """Use V_liquid=m/rho and V_tank=V_liquid/(1-ullage_fraction).

    Density must carry its own temperature/pressure/source. Ullage is a
    fraction of TOTAL tank volume, not a percentage added to liquid volume.

    References
    ----------
    NASA-GAS: density is mass divided by volume.
    NASA-SP8112: ullage and stored-liquid tank pressurization.
    """
    m = quantity(propellant_mass, "kg", minimum=0)
    rho = quantity(density, "kg/m^3", strictly_positive=True)
    ullage = quantity(ullage_fraction, "1", minimum=0)
    if ullage >= 1:
        raise DomainError("Ullage fraction must be less than one.")
    m, rho, ullage = map(Fraction, (m, rho, ullage))
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
    """Size a rigid isothermal IDEAL-gas tank from its usable pressure interval.

    V=m_usable*R*T/(p0-pf); m_loaded=p0*V/(R*T). This is an analytical
    ideal-gas model, not extrapolation of a measured gas density or thruster map.

    References
    ----------
    NASA-GAS: p=rho R T, absolute pressure/temperature, ideal-gas domain.
    NASA-SP8112: blowdown inventory and cutoff boundary.
    """
    m = quantity(usable_mass, "kg", minimum=0)
    p0 = quantity(initial_pressure, "Pa", strictly_positive=True)
    pf = quantity(cutoff_pressure, "Pa", strictly_positive=True)
    t = quantity(temperature, "K", strictly_positive=True)
    gas_r = quantity(specific_gas_constant, "J/(kg*K)", strictly_positive=True)
    if p0 <= pf:
        raise DomainError("Initial absolute pressure must exceed cutoff pressure.")
    m, p0, pf, t, gas_r = map(Fraction, (m, p0, pf, t, gas_r))
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
    """Size an ideal isothermal bottle for displaced liquid volume.

    V_b=p_t*V_displaced/(p_b0-p_bf). Bottle and receiving ullage have the same
    supplied temperature; initial tank pressurant, lines and regulator losses
    are excluded. Bottle final pressure must stay above regulated tank pressure.

    References
    ----------
    NASA-GAS: pV=mRT and conservation of transferred gas.
    NASA-SP8112: regulated bottle-to-ullage pressurization architecture.
    """
    displaced = quantity(displaced_volume, "m^3", minimum=0)
    pt = quantity(tank_pressure, "Pa", strictly_positive=True)
    p0 = quantity(bottle_initial_pressure, "Pa", strictly_positive=True)
    pf = quantity(bottle_final_pressure, "Pa", strictly_positive=True)
    t = quantity(temperature, "K", strictly_positive=True)
    gas_r = quantity(specific_gas_constant, "J/(kg*K)", strictly_positive=True)
    if not p0 > pf > pt:
        raise DomainError("Require bottle initial > bottle final > regulated tank pressure.")
    displaced, pt, p0, pf, t, gas_r = map(Fraction, (displaced, pt, p0, pf, t, gas_r))
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
