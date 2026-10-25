"""Off-design nozzle regimes driven by equilibrium propellant chemistry.

The regime, internal-shock and separation solver of
:mod:`rocket_propulsion.compressible.off_design_nozzle` is a calorically
perfect gas model [1]_. This module feeds it from a propellant pair instead of
hand-entered ``gamma`` and ``R``: the shifting-equilibrium rocket solution at
the design area ratio [2]_ calibrates a two-parameter equivalent perfect gas
that reproduces exactly

* the equilibrium exit-to-chamber pressure ratio at the design area ratio
  (fixing ``gamma_eq``), and
* the equilibrium characteristic velocity ``c*`` with the equilibrium chamber
  temperature (fixing ``R_eq`` through
  ``c* = sqrt(R T_c) / (sqrt(gamma) (2/(gamma+1))^((gamma+1)/(2(gamma-1))))``
  [3]_).

Mass flow and design-point exit pressure therefore match the equilibrium
solution, so the regime thresholds, shock location and separation screen are
consistent with the chemistry at the design point. The remaining vacuum-Isp
difference is reported as the calibration residual. Away from the design
area ratio the equivalent gas is an approximation, which the result states.

References
----------
.. [1] NACA-TR-1135. Ames Research Staff, "Equations, Tables, and Charts for
       Compressible Flow", NACA Report 1135, 1953.
       https://ntrs.nasa.gov/citations/19930091059
.. [2] NASA RP-1311, "Computer Program for Calculation of Complex Chemical
       Equilibrium Compositions and Applications, Part I: Analysis".
       https://ntrs.nasa.gov/citations/19950013764
.. [3] SUTTON-BIBLARZ-2017. G. P. Sutton and O. Biblarz, "Rocket Propulsion
       Elements", 9th ed., Wiley, 2017.
"""

from __future__ import annotations

from dataclasses import dataclass
from math import isfinite, sqrt

from rocket_propulsion.compressible.isentropic import calculate_isentropic, mach_from_isentropic
from rocket_propulsion.compressible.nozzle import calculate_ideal_nozzle
from rocket_propulsion.compressible.off_design_nozzle import (
    NozzleAltitudeSweep,
    OffDesignNozzleResult,
    calculate_off_design_nozzle,
    generate_nozzle_altitude_sweep,
)
from rocket_propulsion.core.errors import ConvergenceError, DomainError
from rocket_propulsion.core.metadata import WarningMessage

from .equilibrium import EquilibriumProvider
from .rocket_equilibrium import get_propellant_pair, reactants_for_mass_ratio
from .rocket_performance import RocketPerformance, solve_rocket_performance


@dataclass(frozen=True, slots=True)
class EquivalentPerfectGas:
    """Perfect gas calibrated to one shifting-equilibrium design point [1]_, [2]_.

    Attributes
    ----------
    gamma : float
        Ratio of specific heats reproducing the equilibrium ``p_e/p_c`` at the
        design area ratio.
    gas_constant_j_kg_k : float
        Gas constant reproducing the equilibrium ``c*`` with ``T_c``.
    chamber_temperature_k : float
        Equilibrium chamber temperature.
    design_area_ratio : float
        Area ratio of the calibration.
    equilibrium_vacuum_specific_impulse_s : float
        Shifting-equilibrium vacuum Isp at the design area ratio.
    equivalent_vacuum_specific_impulse_s : float
        Perfect-gas vacuum Isp with the calibrated properties.
    vacuum_isp_relative_residual : float
        ``I_eq_gas / I_equilibrium - 1``.

    References
    ----------
    .. [1] NASA RP-1311, "Computer Program for Calculation of Complex Chemical
           Equilibrium Compositions and Applications, Part I: Analysis".
           https://ntrs.nasa.gov/citations/19950013764
    .. [2] SUTTON-BIBLARZ-2017. G. P. Sutton and O. Biblarz, "Rocket Propulsion
           Elements", 9th ed., Wiley, 2017.
    """

    gamma: float
    gas_constant_j_kg_k: float
    chamber_temperature_k: float
    design_area_ratio: float
    equilibrium_vacuum_specific_impulse_s: float
    equivalent_vacuum_specific_impulse_s: float
    vacuum_isp_relative_residual: float


@dataclass(frozen=True, slots=True)
class ChemistryCoupledNozzle:
    """Off-design nozzle solution with chemistry-calibrated gas properties.

    Attributes
    ----------
    propellant_pair : str
        Preset key.
    oxidizer_fuel_ratio, chamber_pressure_pa, area_ratio, throat_area_m2 : float
        Operating point and geometry.
    equilibrium : RocketPerformance
        Shifting-equilibrium design-point solution.
    equivalent_gas : EquivalentPerfectGas
        Calibrated perfect gas used by the regime solver.
    off_design : OffDesignNozzleResult
        Regime solution at ``ambient_pressure_pa``.
    altitude_sweep : NozzleAltitudeSweep or None
        Standard-atmosphere sweep when requested.
    warnings : tuple of WarningMessage
        Model statements.
    """

    propellant_pair: str
    oxidizer_fuel_ratio: float
    chamber_pressure_pa: float
    area_ratio: float
    throat_area_m2: float
    equilibrium: RocketPerformance
    equivalent_gas: EquivalentPerfectGas
    off_design: OffDesignNozzleResult
    altitude_sweep: NozzleAltitudeSweep | None
    warnings: tuple[WarningMessage, ...]


def _pressure_ratio(gamma: float, area_ratio: float) -> float:
    mach = mach_from_isentropic("area_ratio_supersonic", area_ratio, gamma)
    return calculate_isentropic(mach, gamma).pressure_total_ratio


def calibrate_equivalent_gas(performance: RocketPerformance, exit_index: int = 0) -> EquivalentPerfectGas:
    """Fit ``gamma`` and ``R`` of a perfect gas to an equilibrium design point [1]_, [2]_.

    Parameters
    ----------
    performance : RocketPerformance
        Shifting-equilibrium (or frozen) rocket solution.
    exit_index : int, default 0
        Which exit station of ``performance`` is the design point.

    Returns
    -------
    EquivalentPerfectGas
        Calibrated properties and the vacuum-Isp residual.

    Raises
    ------
    ConvergenceError
        If no ``gamma`` in ``(1.05, 1.67)`` reproduces the pressure ratio.

    Notes
    -----
    ``gamma`` is bisected on the supersonic isentropic relation
    ``p_e/p_c(gamma, A_e/A_t)`` [1]_; ``R`` then follows from the perfect-gas
    ``c*`` relation [2]_.

    References
    ----------
    .. [1] NACA-TR-1135. Ames Research Staff, "Equations, Tables, and Charts for
           Compressible Flow", NACA Report 1135, 1953.
           https://ntrs.nasa.gov/citations/19930091059
    .. [2] SUTTON-BIBLARZ-2017. G. P. Sutton and O. Biblarz, "Rocket Propulsion
           Elements", 9th ed., Wiley, 2017.
    """

    design = performance.exits[exit_index]
    area_ratio = design.area_ratio
    target = design.station.pressure_pa / performance.chamber_pressure_pa
    low, high = 1.05, 1.67
    f_low = _pressure_ratio(low, area_ratio) - target
    f_high = _pressure_ratio(high, area_ratio) - target
    if f_low * f_high > 0.0:
        raise ConvergenceError("No perfect-gas gamma reproduces the equilibrium pressure ratio.",
                               code="equivalent_gamma_not_bracketed")
    for _ in range(200):
        middle = 0.5 * (low + high)
        value = _pressure_ratio(middle, area_ratio) - target
        if f_low * value <= 0.0:
            high = middle
        else:
            low, f_low = middle, value
        if high - low < 1.0e-13:
            break
    gamma = 0.5 * (low + high)
    temperature = performance.chamber.temperature_k
    flow_function = sqrt(gamma) * (2.0 / (gamma + 1.0)) ** ((gamma + 1.0) / (2.0 * (gamma - 1.0)))
    gas_constant = (performance.characteristic_velocity_m_s * flow_function) ** 2 / temperature
    ideal = calculate_ideal_nozzle(
        chamber_pressure_pa=performance.chamber_pressure_pa,
        chamber_temperature_k=temperature,
        throat_area_m2=1.0,
        area_ratio=area_ratio,
        ambient_pressure_pa=0.0,
        gamma=gamma,
        gas_constant_j_kg_k=gas_constant,
    )
    equilibrium_isp = design.vacuum_specific_impulse_s
    return EquivalentPerfectGas(
        gamma=gamma,
        gas_constant_j_kg_k=gas_constant,
        chamber_temperature_k=temperature,
        design_area_ratio=area_ratio,
        equilibrium_vacuum_specific_impulse_s=equilibrium_isp,
        equivalent_vacuum_specific_impulse_s=ideal.specific_impulse_s,
        vacuum_isp_relative_residual=ideal.specific_impulse_s / equilibrium_isp - 1.0,
    )


def solve_chemistry_coupled_nozzle(
    *,
    propellant_pair: str,
    chamber_pressure_pa: float,
    area_ratio: float,
    ambient_pressure_pa: float,
    throat_area_m2: float = 0.01,
    oxidizer_fuel_ratio: float | None = None,
    include_altitude_sweep: bool = False,
    maximum_altitude_m: float = 50_000.0,
    point_count: int = 51,
    provider: EquilibriumProvider | None = None,
) -> ChemistryCoupledNozzle:
    """Solve nozzle regimes at back pressure with properties from propellant chemistry.

    Parameters
    ----------
    propellant_pair : str
        Preset key such as ``"LOX_RP1"``.
    chamber_pressure_pa : float
        Chamber pressure in Pa.
    area_ratio : float
        Exit-to-throat area ratio, greater than one.
    ambient_pressure_pa : float
        Back pressure in Pa.
    throat_area_m2 : float, default 0.01
        Throat area in m^2.
    oxidizer_fuel_ratio : float, optional
        Mass ratio; the preset reference ratio by default.
    include_altitude_sweep : bool, default False
        Also return a standard-atmosphere sweep with the same gas.
    maximum_altitude_m, point_count
        Sweep settings, as in
        :func:`~rocket_propulsion.compressible.off_design_nozzle.generate_nozzle_altitude_sweep`.
    provider : EquilibriumProvider, optional
        Equilibrium provider; the local provider by default.

    Returns
    -------
    ChemistryCoupledNozzle
        Equilibrium design point, calibrated gas, regime solution and sweep.

    Raises
    ------
    InputError
        If the preset is unknown.
    DomainError
        If an operating input is out of range.
    ConvergenceError
        If the equilibrium or calibration solve fails.

    References
    ----------
    .. [1] NASA RP-1311, "Computer Program for Calculation of Complex Chemical
           Equilibrium Compositions and Applications, Part I: Analysis".
           https://ntrs.nasa.gov/citations/19950013764
    .. [2] NACA-TR-1135. Ames Research Staff, "Equations, Tables, and Charts for
           Compressible Flow", NACA Report 1135, 1953.
           https://ntrs.nasa.gov/citations/19930091059
    .. [3] SUTTON-BIBLARZ-2017. G. P. Sutton and O. Biblarz, "Rocket Propulsion
           Elements", 9th ed., Wiley, 2017.
    """

    pair = get_propellant_pair(propellant_pair)
    ratio = pair.reference_oxidizer_fuel_ratio if oxidizer_fuel_ratio is None else oxidizer_fuel_ratio
    if not isfinite(area_ratio) or area_ratio <= 1.0:
        raise DomainError("Area ratio must exceed one.", field="area_ratio")
    performance = solve_rocket_performance(
        reactants=reactants_for_mass_ratio(pair, ratio),
        chamber_pressure_pa=chamber_pressure_pa,
        area_ratios=(area_ratio,),
        provider=provider,
    )
    gas = calibrate_equivalent_gas(performance)
    off_design = calculate_off_design_nozzle(
        chamber_pressure_pa=chamber_pressure_pa,
        chamber_temperature_k=gas.chamber_temperature_k,
        throat_area_m2=throat_area_m2,
        area_ratio=area_ratio,
        ambient_pressure_pa=ambient_pressure_pa,
        gamma=gas.gamma,
        gas_constant_j_kg_k=gas.gas_constant_j_kg_k,
    )
    sweep = None
    if include_altitude_sweep:
        sweep = generate_nozzle_altitude_sweep(
            chamber_pressure_pa=chamber_pressure_pa,
            chamber_temperature_k=gas.chamber_temperature_k,
            throat_area_m2=throat_area_m2,
            area_ratio=area_ratio,
            gamma=gas.gamma,
            gas_constant_j_kg_k=gas.gas_constant_j_kg_k,
            maximum_altitude_m=maximum_altitude_m,
            point_count=point_count,
        )
    return ChemistryCoupledNozzle(
        propellant_pair=pair.key,
        oxidizer_fuel_ratio=ratio,
        chamber_pressure_pa=chamber_pressure_pa,
        area_ratio=area_ratio,
        throat_area_m2=throat_area_m2,
        equilibrium=performance,
        equivalent_gas=gas,
        off_design=off_design,
        altitude_sweep=sweep,
        warnings=(WarningMessage(
            code="equivalent_perfect_gas",
            message=(
                "Regimes use a perfect gas calibrated to the shifting-equilibrium design point "
                f"(gamma {gas.gamma:.4f}, R {gas.gas_constant_j_kg_k:.1f} J/kg/K); mass flow and "
                "design exit pressure match the equilibrium solution, other stations are approximate."
            ),
        ),),
    )
