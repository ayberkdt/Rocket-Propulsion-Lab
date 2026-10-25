"""Infinite-area-combustor rocket performance with shifting or frozen chemistry.

The chamber is an assigned-enthalpy, assigned-pressure equilibrium state [1]_.
Expansion is isentropic. With shifting equilibrium the composition is
re-equilibrated at every station; with frozen chemistry it keeps the chamber
composition [1]_. The throat is the pressure of maximum mass flux ``rho u``
along the isentrope, which is the sonic condition of an isentropic expansion
[2]_. Each exit station is the supersonic pressure whose mass flux gives the
requested area ratio. Performance follows the usual definitions [2]_::

    c* = p_c / (rho_t u_t)
    I_vac = u_e + p_e (A_e / A_t) / (rho_t u_t)
    C_F,vac = (rho_t u_t u_e + p_e A_e / A_t) / p_c

This is the same problem NASA CEA solves in its ``rocket`` mode with an
infinite-area combustor [1]_, so the result can be compared point for point
with CEA through :mod:`rocket_propulsion.thermochemistry.cea_library`.

References
----------
.. [1] NASA RP-1311, "Computer Program for Calculation of Complex Chemical
       Equilibrium Compositions and Applications, Part I: Analysis".
       https://ntrs.nasa.gov/citations/19950013764
.. [2] SUTTON-BIBLARZ-2017. G. P. Sutton and O. Biblarz, "Rocket Propulsion
       Elements", 9th ed., Wiley, 2017.
"""

from __future__ import annotations

from collections.abc import Callable, Mapping, Sequence
from dataclasses import dataclass
from math import exp, isfinite, log, sqrt
from typing import Literal

from rocket_propulsion.core.constants import STANDARD_GRAVITY_M_S2
from rocket_propulsion.core.errors import ConvergenceError, DomainError, InputError
from rocket_propulsion.core.metadata import WarningMessage

from .equilibrium import EquilibriumProblem, EquilibriumProvider, LocalEquilibriumProvider
from .mixtures import Mixture

Chemistry = Literal["equilibrium", "frozen"]
_GOLDEN = (sqrt(5.0) - 1.0) / 2.0


@dataclass(frozen=True, slots=True)
class RocketStation:
    """Thermodynamic state and flow at one nozzle station.

    Attributes
    ----------
    name : {"chamber", "throat", "exit"}
        Station kind.
    area_ratio : float or None
        Local area over throat area; ``None`` for the infinite-area chamber.
    pressure_pa, temperature_k : float
        Static pressure and temperature.
    density_kg_m3 : float
        Ideal-gas density ``p / (R T)``.
    molecular_mass_kg_mol : float
        Mixture molar mass.
    enthalpy_j_kg, entropy_j_kg_k : float
        Specific enthalpy and entropy of the mixture.
    velocity_m_s : float
        Flow velocity ``sqrt(2 (h_c - h))``.
    mole_fractions : tuple of (str, float)
        Species with mole fraction of at least 1e-5, largest first.
    """

    name: str
    area_ratio: float | None
    pressure_pa: float
    temperature_k: float
    density_kg_m3: float
    molecular_mass_kg_mol: float
    enthalpy_j_kg: float
    entropy_j_kg_k: float
    velocity_m_s: float
    mole_fractions: tuple[tuple[str, float], ...]


@dataclass(frozen=True, slots=True)
class RocketExit:
    """Vacuum performance at one supersonic area ratio.

    Attributes
    ----------
    area_ratio : float
        Exit-to-throat area ratio.
    station : RocketStation
        Exit state.
    effective_exhaust_velocity_m_s : float
        Exit velocity ``u_e`` (CEA ``Isp`` for matched expansion).
    vacuum_effective_exhaust_velocity_m_s : float
        ``u_e + p_e A_e / mdot``.
    vacuum_specific_impulse_s : float
        Vacuum effective exhaust velocity over standard gravity.
    vacuum_thrust_coefficient : float
        ``F_vac / (p_c A_t)``.
    """

    area_ratio: float
    station: RocketStation
    effective_exhaust_velocity_m_s: float
    vacuum_effective_exhaust_velocity_m_s: float
    vacuum_specific_impulse_s: float
    vacuum_thrust_coefficient: float


@dataclass(frozen=True, slots=True)
class RocketPerformance:
    """Chamber, throat and exit results of one rocket problem.

    Attributes
    ----------
    provider, provider_version : str
        Solver that produced the states.
    chemistry : {"equilibrium", "frozen"}
        Expansion chemistry.
    reactants : tuple of (str, float)
        Reactant keys and amounts in mol.
    chamber_pressure_pa : float
        Chamber pressure.
    chamber, throat : RocketStation
        Infinite-area chamber and sonic throat.
    exits : tuple of RocketExit
        Requested supersonic exits in input order.
    characteristic_velocity_m_s : float
        ``c* = p_c / (rho_t u_t)``.
    warnings : tuple of WarningMessage
        Model statements.
    """

    provider: str
    provider_version: str
    chemistry: str
    reactants: tuple[tuple[str, float], ...]
    chamber_pressure_pa: float
    chamber: RocketStation
    throat: RocketStation
    exits: tuple[RocketExit, ...]
    characteristic_velocity_m_s: float
    warnings: tuple[WarningMessage, ...]


@dataclass(frozen=True, slots=True)
class _State:
    pressure_pa: float
    temperature_k: float
    density_kg_m3: float
    molecular_mass_kg_mol: float
    enthalpy_j_kg: float
    entropy_j_kg_k: float
    mole_fractions: tuple[tuple[str, float], ...]


def _illinois(function: Callable[[float], float], low: float, high: float,
              f_low: float, f_high: float, tolerance: float) -> float:
    """Return a root of ``function`` bracketed by ``[low, high]`` (Illinois method)."""
    side = 0
    for _ in range(100):
        root = (low * f_high - high * f_low) / (f_high - f_low)
        value = function(root)
        if abs(value) <= tolerance or abs(high - low) <= tolerance * max(1.0, abs(root)):
            return root
        if value * f_high > 0.0:
            high, f_high = root, value
            if side == -1:
                f_low *= 0.5
            side = -1
        else:
            low, f_low = root, value
            if side == 1:
                f_high *= 0.5
            side = 1
    raise ConvergenceError("Rocket station root did not converge.", code="rocket_root_not_converged")


def _bracket(function: Callable[[float], float], guess: float, lower_limit: float,
             upper_limit: float) -> tuple[float, float, float, float]:
    """Grow a multiplicative bracket around ``guess`` until ``function`` changes sign."""
    low, high = max(lower_limit, guess / 1.05), min(upper_limit, guess * 1.05)
    f_low, f_high = function(low), function(high)
    for _ in range(60):
        if f_low * f_high <= 0.0:
            return low, high, f_low, f_high
        if abs(f_low) < abs(f_high) and low > lower_limit:
            low = max(lower_limit, low / 1.2)
            f_low = function(low)
        elif high < upper_limit:
            high = min(upper_limit, high * 1.2)
            f_high = function(high)
        else:
            break
    raise DomainError("Rocket station is not bracketed by the model range.",
                      code="rocket_station_not_bracketed")


class _Expansion:
    """Isentropic states from one chamber, shifting or frozen."""

    def __init__(self, reactants: Mapping[str, float], provider: EquilibriumProvider,
                 chemistry: Chemistry, chamber_pressure_pa: float) -> None:
        self.reactants = dict(reactants)
        self.provider = provider
        self.chemistry = chemistry
        chamber = provider.solve(EquilibriumProblem.hp(self.reactants, pressure_pa=chamber_pressure_pa))
        mass = chamber.total_mass_kg
        fractions = {item.species: item.mole_fraction for item in chamber.species
                     if item.mole_fraction > 1.0e-10}
        total = sum(fractions.values())
        self.frozen_mixture = Mixture.from_fractions({key: value / total for key, value in fractions.items()})
        self.chamber = _State(
            pressure_pa=chamber.pressure_pa,
            temperature_k=chamber.temperature_k,
            density_kg_m3=chamber.pressure_pa / (chamber.gas_constant_j_kg_k * chamber.temperature_k),
            molecular_mass_kg_mol=chamber.molecular_mass_kg_mol,
            enthalpy_j_kg=chamber.enthalpy_j / mass,
            entropy_j_kg_k=chamber.entropy_j_k / mass,
            mole_fractions=_leading(fractions),
        )
        self.gamma = chamber.gamma

    def state(self, temperature_k: float, pressure_pa: float) -> _State:
        """Return the state at ``(T, p)`` with this expansion's chemistry."""
        if self.chemistry == "frozen":
            properties = self.frozen_mixture.properties(temperature_k, pressure_pa=pressure_pa)
            return _State(pressure_pa, temperature_k,
                          pressure_pa / (properties.gas_constant_j_kg_k * temperature_k),
                          properties.molecular_mass_kg_mol, properties.enthalpy_j_kg,
                          properties.entropy_j_kg_k, self.chamber.mole_fractions)
        result = self.provider.solve(EquilibriumProblem.tp(self.reactants, temperature_k=temperature_k,
                                                           pressure_pa=pressure_pa))
        mass = result.total_mass_kg
        return _State(pressure_pa, temperature_k,
                      pressure_pa / (result.gas_constant_j_kg_k * temperature_k),
                      result.molecular_mass_kg_mol, result.enthalpy_j / mass, result.entropy_j_k / mass,
                      _leading({item.species: item.mole_fraction for item in result.species}))

    def isentrope(self, pressure_pa: float) -> _State:
        """Return the state on the chamber isentrope at ``pressure_pa``."""
        target = self.chamber.entropy_j_kg_k
        chamber = self.chamber
        guess = chamber.temperature_k * (pressure_pa / chamber.pressure_pa) ** ((self.gamma - 1.0) / self.gamma)
        cache: dict[float, _State] = {}

        def residual(temperature: float) -> float:
            cache[temperature] = self.state(temperature, pressure_pa)
            return (cache[temperature].entropy_j_kg_k - target) / abs(target)

        low, high, f_low, f_high = _bracket(residual, min(guess, chamber.temperature_k), 200.0,
                                            chamber.temperature_k)
        temperature = _illinois(residual, low, high, f_low, f_high, 1.0e-11)
        return cache.get(temperature) or self.state(temperature, pressure_pa)

    def velocity(self, state: _State) -> float:
        drop = self.chamber.enthalpy_j_kg - state.enthalpy_j_kg
        if drop <= 0.0:
            raise DomainError("Isentropic expansion produced no enthalpy drop.",
                              code="invalid_rocket_expansion")
        return sqrt(2.0 * drop)

    def mass_flux(self, pressure_pa: float) -> tuple[float, _State]:
        state = self.isentrope(pressure_pa)
        return state.density_kg_m3 * self.velocity(state), state


def _leading(fractions: Mapping[str, float]) -> tuple[tuple[str, float], ...]:
    """Return species with mole fraction >= 1e-5, largest first."""
    return tuple(sorted(((key, value) for key, value in fractions.items() if value >= 1.0e-5),
                        key=lambda item: -item[1]))


def _station(name: str, area_ratio: float | None, state: _State, velocity: float) -> RocketStation:
    return RocketStation(name, area_ratio, state.pressure_pa, state.temperature_k, state.density_kg_m3,
                         state.molecular_mass_kg_mol, state.enthalpy_j_kg, state.entropy_j_kg_k,
                         velocity, state.mole_fractions)


def solve_rocket_performance(
    *,
    reactants: Mapping[str, float],
    chamber_pressure_pa: float,
    area_ratios: Sequence[float] = (10.0, 25.0, 50.0),
    chemistry: Chemistry = "equilibrium",
    provider: EquilibriumProvider | None = None,
) -> RocketPerformance:
    """Solve an infinite-area-combustor rocket problem [1]_, [2]_.

    Parameters
    ----------
    reactants : mapping of str to float
        Reactant amounts in mol; condensed reactants such as ``"O2(L)"`` or
        ``"RP-1"`` enter with their assigned enthalpy.
    chamber_pressure_pa : float
        Chamber pressure in Pa.
    area_ratios : sequence of float, default (10, 25, 50)
        Supersonic exit-to-throat area ratios, each greater than one.
    chemistry : {"equilibrium", "frozen"}, default "equilibrium"
        Shifting equilibrium, or composition frozen at the chamber.
    provider : EquilibriumProvider, optional
        Equilibrium provider; the local provider by default.

    Returns
    -------
    RocketPerformance
        Chamber, throat and exit stations with ``c*``, vacuum Isp and
        vacuum thrust coefficient.

    Raises
    ------
    InputError
        If ``chemistry`` is unknown.
    DomainError
        If the pressure or an area ratio is out of range, or a station is not
        bracketed inside 200 K to the chamber temperature.
    ConvergenceError
        If a station root does not converge.

    Notes
    -----
    The throat pressure maximizes ``rho u`` by golden-section search on
    ``ln p`` between ``p_c/2.5`` and ``p_c/1.2`` to 1e-5; the exits solve
    ``(rho u)_t / (rho u)_e = A_e/A_t`` on the supersonic branch by the
    Illinois method on ``ln p``. Every state lies on the chamber isentrope
    [1]_. Definitions of ``c*``, vacuum Isp and ``C_F`` follow [2]_.

    References
    ----------
    .. [1] NASA RP-1311, "Computer Program for Calculation of Complex Chemical
           Equilibrium Compositions and Applications, Part I: Analysis".
           https://ntrs.nasa.gov/citations/19950013764
    .. [2] SUTTON-BIBLARZ-2017. G. P. Sutton and O. Biblarz, "Rocket Propulsion
           Elements", 9th ed., Wiley, 2017.
    """

    if chemistry not in ("equilibrium", "frozen"):
        raise InputError("Rocket chemistry must be equilibrium or frozen.", field="chemistry")
    if not isfinite(chamber_pressure_pa) or chamber_pressure_pa <= 0.0:
        raise DomainError("Chamber pressure must be positive.", field="chamber_pressure_pa")
    ratios = tuple(float(value) for value in area_ratios)
    if not ratios or any(not isfinite(value) or value <= 1.0 for value in ratios):
        raise DomainError("Rocket area ratios must all exceed one.", field="area_ratios")
    active = provider or LocalEquilibriumProvider()
    expansion = _Expansion(reactants, active, chemistry, chamber_pressure_pa)
    pc = chamber_pressure_pa

    # Golden-section maximization of mass flux on ln p.
    low, high = log(pc / 2.5), log(pc / 1.2)
    first = high - _GOLDEN * (high - low)
    second = low + _GOLDEN * (high - low)
    flux_first = expansion.mass_flux(exp(first))[0]
    flux_second = expansion.mass_flux(exp(second))[0]
    while high - low > 1.0e-5:
        if flux_first > flux_second:
            high, second, flux_second = second, first, flux_first
            first = high - _GOLDEN * (high - low)
            flux_first = expansion.mass_flux(exp(first))[0]
        else:
            low, first, flux_first = first, second, flux_second
            second = low + _GOLDEN * (high - low)
            flux_second = expansion.mass_flux(exp(second))[0]
    throat_flux, throat_state = expansion.mass_flux(exp(0.5 * (low + high)))
    throat_velocity = expansion.velocity(throat_state)
    log_throat = log(throat_state.pressure_pa)

    exits: list[RocketExit] = []
    for ratio in ratios:
        states: dict[float, _State] = {}

        def residual(log_pressure: float, ratio: float = ratio, *, states=states) -> float:
            flux, state = expansion.mass_flux(exp(log_pressure))
            states[log_pressure] = state
            return log(throat_flux / flux) - log(ratio)

        guess = log_throat - log(ratio) * expansion.gamma
        lo, hi, f_lo, f_hi = _bracket(residual, guess, log_throat - 40.0,
                                      log_throat - 1.0e-3)
        root = _illinois(residual, lo, hi, f_lo, f_hi, 1.0e-12)
        state = states.get(root) or expansion.isentrope(exp(root))
        velocity = expansion.velocity(state)
        vacuum = velocity + state.pressure_pa * ratio / throat_flux
        exits.append(RocketExit(
            area_ratio=ratio,
            station=_station("exit", ratio, state, velocity),
            effective_exhaust_velocity_m_s=velocity,
            vacuum_effective_exhaust_velocity_m_s=vacuum,
            vacuum_specific_impulse_s=vacuum / STANDARD_GRAVITY_M_S2,
            vacuum_thrust_coefficient=(throat_flux * velocity + state.pressure_pa * ratio) / pc,
        ))
    return RocketPerformance(
        provider=active.name,
        provider_version=active.version,
        chemistry=chemistry,
        reactants=tuple(sorted((key.upper(), float(value)) for key, value in reactants.items())),
        chamber_pressure_pa=pc,
        chamber=_station("chamber", None, expansion.chamber, 0.0),
        throat=_station("throat", 1.0, throat_state, throat_velocity),
        exits=tuple(exits),
        characteristic_velocity_m_s=pc / throat_flux,
        warnings=(WarningMessage(
            code="ideal_rocket_assumptions",
            message=("Infinite-area combustor, adiabatic one-dimensional isentropic flow; "
                     "no kinetics, friction, divergence or condensed-phase lag."),
        ),),
    )
