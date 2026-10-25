"""Steady quasi-one-dimensional finite-rate expansion through a prescribed nozzle.

The flow is steady, inviscid, adiabatic and one-dimensional, with the
assigned-area equations of NASA RP-1328 section 2.4 [1]_: continuity
``rho V A = mdot`` (eq. (2.58)), momentum ``dp + rho V dV = 0`` (eq. (2.60)),
energy ``dh + V dV = 0`` (eq. (2.61) without heat loss) and species
``d sigma_i/dx = omega_i/(rho V)`` (eqs. (2.33), (2.57)), where
``sigma_i`` is mol of species ``i`` per kg of mixture. Eliminating ``V dV``
gives RP-1328 eq. (2.65); with the ideal-mixture entropy of NASA RP-1311
eqs. (2.16)--(2.17) [2]_ this is ``T ds = -sum_i mu_i d sigma_i``. The marched
state is therefore ``(sigma_i, s)``; at every ``x`` the temperature and
pressure follow algebraically from the entropy, the total enthalpy and the
mass flux ``mdot/A(x)`` on the subsonic or supersonic branch of the local
fixed-composition isentrope. That closure is singular at the frozen sonic
point (``M_f = 1``, the ``gamma M^2 - 1`` denominator of RP-1328 eqs.
(2.80)--(2.85)), a saddle of the equations.

Two start modes are provided. ``"equilibrium_throat"`` (default) assumes
shifting equilibrium from the infinite-area chamber through the throat to a
supersonic start station and integrates finite rates downstream of it; very
fast rates then reproduce shifting equilibrium. ``"chamber"`` integrates
finite rates from the convergent section with the chamber composition, finds
the choking mass flow by shooting and switches branch near the frozen sonic
point; zero rates then reproduce expansion frozen at the chamber, but the
supersonic continuation is reliable only when the throat chemistry is slow.
Vacuum specific impulse is ``(V_e + p_e A_e/mdot)/g0`` (RP-1311 eqs.
(6.8)--(6.9) [2]_).

The Bray sudden-freezing estimate follows K. N. C. Bray's criterion
``-(d alpha/dt)_e = K r_De`` (thesis eq. (36), with ``K`` of order one and
experimentally 0.42 to 1.25) [3]_: equilibrium expansion up to the freezing
point, frozen expansion after it.

References
----------
.. [1] NASA-RP-1328. K. Radhakrishnan, "LSENS -- A General Chemical Kinetics
       and Sensitivity Analysis Code for Homogeneous Gas-Phase Reactions.
       I. Theory and Numerical Solution Procedures", NASA RP-1328, 1994,
       section 2.4, pp. 16-24. https://ntrs.nasa.gov/citations/19950008642
.. [2] NASA RP-1311, "Computer Program for Calculation of Complex Chemical
       Equilibrium Compositions and Applications, Part I: Analysis".
       https://ntrs.nasa.gov/citations/19950013764
.. [3] BRAY-1960-THESIS. K. N. C. Bray, "Atomic recombination in a hypersonic
       wind-tunnel nozzle", PhD thesis, University of Southampton, 1960,
       section 5, pp. 35-41. https://eprints.soton.ac.uk/458375/1/38084.pdf
"""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass
from math import exp, isfinite, log, pi, radians, sqrt, tan
from typing import Literal

from rocket_propulsion.core.constants import MOLAR_GAS_CONSTANT_J_MOL_K, STANDARD_GRAVITY_M_S2
from rocket_propulsion.core.errors import ConvergenceError, DomainError, InputError
from rocket_propulsion.core.metadata import WarningMessage

from .equilibrium import EquilibriumProblem, EquilibriumProvider, LocalEquilibriumProvider
from .kinetics import Mechanism, StepRejected, _compile, solve_stiff
from .kinetics_mechanisms import gri30_h2o2_mechanism
from .mixtures import REFERENCE_PRESSURE_PA
from .species import get_reactant

_R = MOLAR_GAS_CONSTANT_J_MOL_K
_SQRT_EPS = 1.4901161193847656e-08
_MINIMUM_TEMPERATURE_K = 300.0
_SONIC_CLAMP_RELATIVE = 1.0e-4


@dataclass(frozen=True, slots=True)
class HyperbolicNozzle:
    """Axisymmetric nozzle with hyperbolic convergent and divergent walls.

    ``r(x)^2 = r_t^2 + (x tan(theta))^2`` with ``theta`` the convergent half
    angle for ``x < 0`` and the divergent half angle for ``x > 0``, so
    ``A/A_t = 1 + (x tan(theta)/r_t)^2``; each side is asymptotic to a cone
    and the wall curvature radius at the throat is ``r_t/tan(theta)^2``.

    Attributes
    ----------
    throat_radius_m : float
        Throat radius ``r_t``; sets the residence time and hence the kinetic
        loss.
    contraction_area_ratio : float
        Area ratio at the start of the convergent section (> 1).
    exit_area_ratio : float
        Exit-to-throat area ratio (> 1).
    convergent_half_angle_deg, divergent_half_angle_deg : float
        Asymptotic half angles in degrees, in (0, 90).

    Raises
    ------
    DomainError
        If a dimension is out of range.
    """

    throat_radius_m: float
    contraction_area_ratio: float
    exit_area_ratio: float
    convergent_half_angle_deg: float
    divergent_half_angle_deg: float

    def __post_init__(self) -> None:
        if not isfinite(self.throat_radius_m) or self.throat_radius_m <= 0.0:
            raise DomainError("Throat radius must be positive.", field="throat_radius_m")
        for name in ("contraction_area_ratio", "exit_area_ratio"):
            value = getattr(self, name)
            if not isfinite(value) or value <= 1.0:
                raise DomainError(f"{name} must exceed one.", field=name)
        for name in ("convergent_half_angle_deg", "divergent_half_angle_deg"):
            value = getattr(self, name)
            if not isfinite(value) or not 0.0 < value < 90.0:
                raise DomainError(f"{name} must be in (0, 90) degrees.", field=name)

    @property
    def throat_area_m2(self) -> float:
        """Throat area ``pi r_t^2``."""
        return pi * self.throat_radius_m ** 2

    def _slope(self, x_m: float) -> float:
        angle = self.convergent_half_angle_deg if x_m < 0.0 else self.divergent_half_angle_deg
        return tan(radians(angle))

    def area_ratio(self, x_m: float) -> float:
        """Return ``A(x)/A_t``; ``x`` is measured from the throat."""
        return 1.0 + (x_m * self._slope(x_m) / self.throat_radius_m) ** 2

    def position(self, area_ratio: float, *, divergent: bool) -> float:
        """Return ``x`` of an area ratio on the chosen side of the throat."""
        if area_ratio < 1.0:
            raise DomainError("Area ratio must be at least one.", field="area_ratio")
        sign = 1.0 if divergent else -1.0
        slope = self._slope(sign)
        return sign * self.throat_radius_m * sqrt(area_ratio - 1.0) / slope

    @property
    def start_x_m(self) -> float:
        """``x`` of the contraction area ratio (negative)."""
        return self.position(self.contraction_area_ratio, divergent=False)

    @property
    def exit_x_m(self) -> float:
        """``x`` of the exit area ratio."""
        return self.position(self.exit_area_ratio, divergent=True)


@dataclass(frozen=True, slots=True)
class NozzleStation:
    """Flow state at one axial station.

    Attributes
    ----------
    x_m, area_ratio : float
        Position from the throat and ``A/A_t``.
    temperature_k, pressure_pa, density_kg_m3, velocity_m_s : float
        Static state and velocity.
    frozen_mach : float
        ``V`` over the frozen sound speed ``sqrt(gamma_f R T / M)``.
    mole_fractions : tuple of (str, float)
        All mechanism species in mechanism order.
    """

    x_m: float
    area_ratio: float
    temperature_k: float
    pressure_pa: float
    density_kg_m3: float
    velocity_m_s: float
    frozen_mach: float
    mole_fractions: tuple[tuple[str, float], ...]


@dataclass(frozen=True, slots=True)
class NozzleKineticsResult:
    """Finite-rate nozzle expansion result.

    Attributes
    ----------
    mechanism : str
        Mechanism name.
    rate_multiplier : float
        Common forward-rate factor used.
    start : str
        ``"equilibrium_throat"`` or ``"chamber"`` (see
        :func:`solve_nozzle_kinetics`).
    chamber_temperature_k, chamber_pressure_pa : float
        Equilibrium chamber state.
    mass_flux_throat_kg_m2_s : float
        ``mdot / A_t``.
    characteristic_velocity_m_s : float
        ``p_c A_t / mdot``.
    start_x_m : float
        Start of the supersonic finite-rate march: the equilibrium start
        station, or the branch switch near the frozen sonic point.
    exit : NozzleStation
        Exit station.
    vacuum_specific_impulse_s : float
        ``(V_e + p_e A_e/mdot)/g0``.
    stations : tuple of NozzleStation
        Accepted integration points from the start to the exit.
    mass_flow_iterations : int
        Subsonic marches used to find the choking mass flow (``0`` in the
        ``"equilibrium_throat"`` mode).
    accepted_steps : int
        Accepted integrator steps of the final subsonic and supersonic marches.
    maximum_element_residual : float
        Largest relative change of an element's mol/kg between chamber and
        exit.
    warnings : tuple of WarningMessage
        Model statements.
    """

    mechanism: str
    rate_multiplier: float
    start: str
    chamber_temperature_k: float
    chamber_pressure_pa: float
    mass_flux_throat_kg_m2_s: float
    characteristic_velocity_m_s: float
    start_x_m: float
    exit: NozzleStation
    vacuum_specific_impulse_s: float
    stations: tuple[NozzleStation, ...]
    mass_flow_iterations: int
    accepted_steps: int
    maximum_element_residual: float
    warnings: tuple[WarningMessage, ...]


@dataclass(frozen=True, slots=True)
class BrayFreezingResult:
    """Bray sudden-freezing estimate along the same nozzle.

    Attributes
    ----------
    freezing_constant : float
        ``K`` of the criterion.
    freezing_x_m, freezing_area_ratio : float
        Freezing point (the exit when the criterion is never met).
    freezing_temperature_k, freezing_pressure_pa : float
        Equilibrium state at the freezing point.
    equilibrium_mass_flux_throat_kg_m2_s : float
        Shifting-equilibrium choking flux ``mdot/A_t``.
    vacuum_specific_impulse_s : float
        Equilibrium to the freezing point, frozen after it.
    equilibrium_vacuum_specific_impulse_s : float
        Shifting equilibrium to the exit over the mechanism species.
    criterion : tuple of (float, float)
        ``(A/A_t, -(d sigma_m/dt)_e / r_De)`` along the divergent section.
    warnings : tuple of WarningMessage
        Model statements.
    """

    freezing_constant: float
    freezing_x_m: float
    freezing_area_ratio: float
    freezing_temperature_k: float
    freezing_pressure_pa: float
    equilibrium_mass_flux_throat_kg_m2_s: float
    vacuum_specific_impulse_s: float
    equilibrium_vacuum_specific_impulse_s: float
    criterion: tuple[tuple[float, float], ...]
    warnings: tuple[WarningMessage, ...]


class _Choked(StepRejected):  # noqa: N818 - control-flow signal
    """The target mass flux exceeds the local frozen choking flux."""

    def __init__(self, excess_flux: float) -> None:
        super().__init__()
        self.local_maximum = excess_flux


class _MarchComplete(Exception):  # noqa: N818 - control-flow signal
    """The subsonic march has passed its minimum choking margin."""


@dataclass(frozen=True, slots=True)
class _Closure:
    temperature_k: float
    pressure_pa: float
    density_kg_m3: float
    velocity_m_s: float
    frozen_sound_speed_m_s: float
    h_rt: tuple[float, ...]
    s_r: tuple[float, ...]
    log_mole_fractions: tuple[float, ...]


def _illinois(function, low, high, f_low, f_high, tolerance):
    """Bracketed Illinois root of ``function`` on ``[low, high]``."""
    side = 0
    for _ in range(200):
        root = (low * f_high - high * f_low) / (f_high - f_low)
        value = function(root)
        if value == 0.0 or abs(high - low) <= tolerance * max(1.0, abs(root)):
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
    raise ConvergenceError("Nozzle closure root did not converge.", code="nozzle_root_not_converged")


class _Expansion:
    """Chamber state, mixture closure and right-hand side (internal)."""

    def __init__(self, reactants: Mapping[str, float], chamber_pressure_pa: float,
                 mechanism: Mechanism, provider: EquilibriumProvider) -> None:
        if not isfinite(chamber_pressure_pa) or chamber_pressure_pa <= 0.0:
            raise DomainError("Chamber pressure must be positive.", field="chamber_pressure_pa")
        elements = set(mechanism.elements)
        for key in reactants:
            reactant_elements = {element for element, _ in get_reactant(key).formula}
            if not reactant_elements <= elements:
                raise InputError(f"Reactant {key} has elements outside the mechanism "
                                 f"({sorted(elements)}).", field="reactants")
        self.mechanism = mechanism
        self.compiled = _compile(mechanism)
        self.size = len(mechanism.species)
        chamber = provider.solve(EquilibriumProblem.hp(
            dict(reactants), pressure_pa=chamber_pressure_pa,
            candidate_species=mechanism.species))
        mass = chamber.total_mass_kg
        moles = {item.species: item.moles for item in chamber.species}
        self.sigma0 = tuple(moles.get(key, 0.0) / mass for key in mechanism.species)
        self.h0 = chamber.enthalpy_j / mass
        self.s0 = chamber.entropy_j_k / mass
        self.chamber_temperature_k = chamber.temperature_k
        self.chamber_pressure_pa = chamber_pressure_pa
        self.reactants = dict(reactants)
        self.provider = provider
        self.element_matrix = mechanism.element_matrix()
        self.rate_multiplier = 1.0
        self.last_temperature: dict[str, float] = {}

    # -- fixed-composition thermodynamics -----------------------------------
    def _mixture(self, sigma):
        total = sum(sigma)
        if not total > 0.0:
            raise DomainError("Mixture has no moles.", code="nozzle_invalid_state")
        log_x = tuple(log(value / total) if value > 0.0 else -690.0 for value in sigma)
        mixing = sum(value * lx for value, lx in zip(sigma, log_x, strict=True) if value > 0.0)
        return total, log_x, mixing

    def _evaluate(self, temperature, sigma, total, mixing, entropy):
        cp_r, h_rt, s_r = self.compiled.thermo(temperature)
        enthalpy = _R * temperature * sum(s * h for s, h in zip(sigma, h_rt, strict=True))
        cp = _R * sum(s * c for s, c in zip(sigma, cp_r, strict=True))
        standard = _R * (sum(s * e for s, e in zip(sigma, s_r, strict=True)) - mixing)
        log_p = (standard - entropy) / (_R * total)
        v2 = 2.0 * (self.h0 - enthalpy)
        gas_constant = _R * total
        gamma = cp / (cp - gas_constant)
        return cp, log_p, v2, gamma, gas_constant, h_rt, s_r

    def _upper_temperature(self, sigma, total, mixing, entropy):
        temperature = max(self.chamber_temperature_k, _MINIMUM_TEMPERATURE_K + 1.0)
        for _ in range(60):
            v2 = self._evaluate(temperature, sigma, total, mixing, entropy)[2]
            if v2 < 0.0:
                return temperature
            temperature = min(temperature * 1.05, 6000.0)
            if temperature >= 6000.0:
                if self._evaluate(6000.0, sigma, total, mixing, entropy)[2] < 0.0:
                    return 6000.0
                break
        raise DomainError("Stagnation temperature exceeds the thermodynamic data range.",
                          code="nozzle_invalid_state")

    def sonic(self, sigma, entropy, guess=None):
        """Return ``(T*, G_max)`` on the fixed-composition isentrope."""
        total, _, mixing = self._mixture(sigma)

        def phi(temperature):
            _, _, v2, gamma, gas_constant, _, _ = self._evaluate(temperature, sigma, total,
                                                                  mixing, entropy)
            return v2 - gamma * gas_constant * temperature

        if guess is not None:
            low, high = guess * 0.995, guess * 1.005
            f_low, f_high = phi(low), phi(high)
        else:
            low, high = _MINIMUM_TEMPERATURE_K, self._upper_temperature(sigma, total, mixing,
                                                                        entropy)
            f_low, f_high = phi(low), phi(high)
        for _ in range(60):
            if f_low > 0.0 > f_high:
                break
            if f_low <= 0.0:
                if low <= _MINIMUM_TEMPERATURE_K:
                    raise DomainError("Sonic point lies below the thermodynamic data range.",
                                      code="nozzle_invalid_state")
                low = max(_MINIMUM_TEMPERATURE_K, low * 0.97)
                f_low = phi(low)
            if f_high >= 0.0:
                high = self._upper_temperature(sigma, total, mixing, entropy)
                f_high = phi(high)
        else:
            raise ConvergenceError("Sonic point not bracketed.", code="nozzle_root_not_converged")
        temperature = _illinois(phi, low, high, f_low, f_high, 1.0e-13)
        _, log_p, v2, _, gas_constant, _, _ = self._evaluate(temperature, sigma, total, mixing,
                                                             entropy)
        density = REFERENCE_PRESSURE_PA * exp(log_p) / (gas_constant * temperature)
        return temperature, density * sqrt(max(v2, 0.0))

    def closure(self, sigma, entropy, flux, branch, *, clamp=False):
        """Return the state with mass flux ``flux`` on ``branch`` ('sub' or 'sup')."""
        total, log_x, mixing = self._mixture(sigma)
        log_target = log(flux)

        def log_flux(temperature):
            cp, log_p, v2, gamma, gas_constant, h_rt, s_r = self._evaluate(
                temperature, sigma, total, mixing, entropy)
            if v2 <= 0.0:
                return None
            value = log_p + log(REFERENCE_PRESSURE_PA / (gas_constant * temperature)) + 0.5 * log(v2)
            derivative = (cp / gas_constant - 1.0) / temperature - cp / v2
            return value, derivative, v2 - gamma * gas_constant * temperature

        temperature = None
        guess = self.last_temperature.get(branch)
        if guess is not None:
            trial = guess
            for _ in range(25):
                if trial < _MINIMUM_TEMPERATURE_K or trial > 6000.0:
                    trial = None
                    break
                evaluated = log_flux(trial)
                if evaluated is None or evaluated[1] == 0.0:
                    trial = None
                    break
                step = -(evaluated[0] - log_target) / evaluated[1]
                step = max(-0.05 * trial, min(0.05 * trial, step))
                trial += step
                if abs(step) <= 1.0e-12 * trial:
                    final = log_flux(trial)
                    if final is not None and (final[2] > 0.0) == (branch == "sup") \
                            and abs(final[0] - log_target) < 1.0e-10:
                        temperature = trial
                    break
        if temperature is None:
            sonic_temperature, maximum = self.sonic(sigma, entropy)
            if flux > maximum:
                if clamp and flux <= maximum * (1.0 + _SONIC_CLAMP_RELATIVE):
                    temperature = sonic_temperature
                else:
                    raise _Choked(maximum)
            else:
                if branch == "sup":
                    low, high = _MINIMUM_TEMPERATURE_K, sonic_temperature
                else:
                    low = sonic_temperature
                    high = self._upper_temperature(sigma, total, mixing, entropy)
                    while log_flux(high) is not None:
                        high *= 1.01

                def residual(value):
                    evaluated = log_flux(value)
                    if evaluated is None:
                        return -50.0 if branch == "sub" else 50.0
                    return evaluated[0] - log_target

                f_low, f_high = residual(low), residual(high)
                if f_low * f_high > 0.0:
                    if branch == "sup" and f_low > 0.0:
                        raise DomainError("Supersonic state lies below the thermodynamic data "
                                          "range.", code="nozzle_invalid_state")
                    temperature = sonic_temperature
                else:
                    temperature = _illinois(residual, low, high, f_low, f_high, 1.0e-13)
        cp, log_p, v2, gamma, gas_constant, h_rt, s_r = self._evaluate(
            temperature, sigma, total, mixing, entropy)
        self.last_temperature[branch] = temperature
        pressure = REFERENCE_PRESSURE_PA * exp(log_p)
        density = pressure / (gas_constant * temperature)
        return _Closure(temperature, pressure, density, sqrt(max(v2, 0.0)),
                        sqrt(gamma * gas_constant * temperature), tuple(h_rt), tuple(s_r), log_x)

    # -- right-hand side -----------------------------------------------------
    def derivatives(self, sigma, entropy, flux, branch, clamp=False):
        state = self.closure(sigma, entropy, flux, branch, clamp=clamp)
        if self.rate_multiplier == 0.0:
            return [0.0] * (self.size + 1), state
        concentrations = [state.density_kg_m3 * max(value, 0.0) for value in sigma]
        constants = self.compiled.equilibrium_constants(state.temperature_k, state.h_rt,
                                                        state.s_r)
        production = self.compiled.rates(state.temperature_k, concentrations,
                                         self.rate_multiplier, constants)[3]
        dsigma = [value / flux for value in production]
        log_pressure = log(state.pressure_pa / REFERENCE_PRESSURE_PA)
        potential = [h - s + lx + log_pressure for h, s, lx in
                     zip(state.h_rt, state.s_r, state.log_mole_fractions, strict=True)]
        dentropy = -_R * sum(m * d for m, d in zip(potential, dsigma, strict=True))
        return dsigma + [dentropy], state

    def jacobian(self, sigma, entropy, flux, branch, clamp=False):
        """Analytic ``d(rhs)/d(sigma, s)`` at fixed ``x`` by implicit differentiation.

        The closure ``ln G(T, sigma, s) = ln(mdot/A)`` gives ``dT/dq``; then
        ``rho = G/V`` gives ``d ln rho/dq = (dh/dq)/V^2``, the concentrations
        ``C = rho sigma`` and the rates follow by the chain rule with the
        analytic ``d omega/d C`` and a one-sided difference in ``T``.
        """
        size = self.size
        n = size + 1
        if self.rate_multiplier == 0.0:
            return [[0.0] * n for _ in range(n)]
        state = self.closure(sigma, entropy, flux, branch, clamp=clamp)
        temperature, density = state.temperature_k, state.density_kg_m3
        total = sum(sigma)
        gas_constant = _R * total
        cp_r, h_rt, s_r = self.compiled.thermo(temperature)
        cp = _R * sum(sg * c for sg, c in zip(sigma, cp_r, strict=True))
        v2 = state.velocity_m_s ** 2
        log_p = log(state.pressure_pa / REFERENCE_PRESSURE_PA)
        log_x = state.log_mole_fractions
        dlnp_dt = cp / (gas_constant * temperature)
        dlng_dt = (cp / gas_constant - 1.0) / temperature - cp / v2
        # At the sonic point d ln G/dT vanishes (the closure is singular there);
        # keep the branch sign and a floor of 1e-3 cp/V^2 so the matrix stays finite.
        floor = 1.0e-3 * cp / v2
        if abs(dlng_dt) < floor:
            dlng_dt = floor if branch == "sup" else -floor
        # partial derivatives at fixed T: index k < size is sigma_k, size is s
        dlnp_dq = [(s_r[k] - log_x[k] - log_p) / total for k in range(size)] + [-1.0 / gas_constant]
        dh_dq = [_R * temperature * h_rt[k] for k in range(size)] + [0.0]
        dlng_dq = [dlnp_dq[k] - 1.0 / total - dh_dq[k] / v2 for k in range(size)] \
            + [dlnp_dq[size]]
        dt_dq = [-value / dlng_dt for value in dlng_dq]
        dlnrho_dq = [(dh_dq[q] + cp * dt_dq[q]) / v2 for q in range(n)]
        concentrations = [density * max(value, 0.0) for value in sigma]
        constants = self.compiled.equilibrium_constants(temperature, h_rt, s_r)
        production = self.compiled.rates(temperature, concentrations, self.rate_multiplier,
                                         constants)[3]
        jac_c = self.compiled.rates_jacobian(temperature, concentrations, self.rate_multiplier,
                                             constants)
        step = 1.0e-6 * temperature
        _, h_rt_t, s_r_t = self.compiled.thermo(temperature + step)
        constants_t = self.compiled.equilibrium_constants(temperature + step, h_rt_t, s_r_t)
        production_t = self.compiled.rates(temperature + step, concentrations,
                                           self.rate_multiplier, constants_t)[3]
        domega_dt = [(a - b) / step for a, b in zip(production_t, production, strict=True)]
        # d C_k / d q
        dc_dq = [[concentrations[k] * dlnrho_dq[q] for q in range(n)] for k in range(size)]
        for k in range(size):
            dc_dq[k][k] += density
        domega_dq = [[sum(jac_c[i][k] * dc_dq[k][q] for k in range(size))
                      + domega_dt[i] * dt_dq[q] for q in range(n)] for i in range(size)]
        matrix = [[value / flux for value in row] for row in domega_dq]
        potential = [h_rt[k] - s_r[k] + log_x[k] + log_p for k in range(size)]
        entropy_row = []
        for q in range(n):
            dlnp_total = dlnp_dq[q] + dlnp_dt * dt_dq[q]
            value = 0.0
            for k in range(size):
                dpotential = -h_rt[k] / temperature * dt_dq[q] + dlnp_total
                if q < size and sigma[k] > 0.0:
                    dpotential += (1.0 if q == k else 0.0) / sigma[k] - 1.0 / total
                elif q < size:
                    dpotential -= 1.0 / total
                value += potential[k] * domega_dq[k][q] + production[k] * dpotential
            entropy_row.append(-_R * value / flux)
        matrix.append(entropy_row)
        return matrix

    def station(self, x_m, area_ratio, y, flux, branch, clamp=False):
        sigma = y[:self.size]
        state = self.closure(sigma, y[self.size], flux, branch, clamp=clamp)
        total = sum(sigma)
        return NozzleStation(
            x_m=x_m, area_ratio=area_ratio, temperature_k=state.temperature_k,
            pressure_pa=state.pressure_pa, density_kg_m3=state.density_kg_m3,
            velocity_m_s=state.velocity_m_s,
            frozen_mach=state.velocity_m_s / state.frozen_sound_speed_m_s,
            mole_fractions=tuple((key, max(value, 0.0) / total)
                                 for key, value in zip(self.mechanism.species, sigma,
                                                       strict=True)))

    def element_residual(self, sigma):
        worst = 0.0
        for row in self.element_matrix:
            reference = sum(a * s for a, s in zip(row, self.sigma0, strict=True))
            current = sum(a * s for a, s in zip(row, sigma, strict=True))
            if reference > 0.0:
                worst = max(worst, abs(current - reference) / reference)
        return worst


def _finite_difference_jacobian(function, x, y):
    """Forward-difference Jacobian with increments ``sqrt(eps) max(|y|, floor)``."""
    size = len(y)
    base = function(x, y)
    matrix = [[0.0] * size for _ in range(size)]
    state = list(y)
    for k in range(size):
        saved = state[k]
        floor = 1.0e-8 if k < size - 1 else 1.0
        delta = _SQRT_EPS * max(abs(saved), floor)
        state[k] = saved + delta
        try:
            shifted = function(x, state)
        finally:
            state[k] = saved
        for i in range(size):
            matrix[i][k] = (shifted[i] - base[i]) / delta
    return matrix


def solve_nozzle_kinetics(
    *,
    reactants: Mapping[str, float],
    chamber_pressure_pa: float,
    nozzle: HyperbolicNozzle,
    mechanism: Mechanism | None = None,
    rate_multiplier: float = 1.0,
    start: Literal["equilibrium_throat", "chamber"] = "equilibrium_throat",
    start_area_ratio: float = 1.1,
    rtol: float = 1.0e-6,
    provider: EquilibriumProvider | None = None,
    mass_flow_tolerance: float = 1.0e-6,
) -> NozzleKineticsResult:
    """Integrate a finite-rate nozzle expansion to the exit [1]_.

    Parameters
    ----------
    reactants : mapping of str to float
        Reactant amounts in mol (gas species or condensed reactants such as
        ``"H2(L)"``); their elements must belong to the mechanism.
    chamber_pressure_pa : float
        Chamber pressure in Pa.
    nozzle : HyperbolicNozzle
        Geometry, including the throat radius that sets the residence time.
    mechanism : Mechanism, optional
        Reaction mechanism; the GRI-Mech 3.0 H2-O2 subset by default [2]_.
    rate_multiplier : float, default 1.0
        Common factor on every forward rate; ``0`` gives frozen flow.
    start : {"equilibrium_throat", "chamber"}, default "equilibrium_throat"
        ``"equilibrium_throat"``: shifting equilibrium from the chamber to
        ``start_area_ratio`` on the divergent side, finite rates after it.
        ``"chamber"``: finite rates from the contraction-ratio station with
        the chamber composition, through the throat (mass flow by shooting).
    start_area_ratio : float, default 1.1
        Supersonic start station of the ``"equilibrium_throat"`` mode; the
        frozen Mach number there must exceed one.
    rtol : float, default 1e-6
        Relative tolerance of the RODAS integrator; absolute tolerances are
        ``1e-12`` mol/kg for species and ``rtol`` J/(kg K) for entropy.
    provider : EquilibriumProvider, optional
        Equilibrium provider; the local provider by default.
    mass_flow_tolerance : float, default 1e-6
        Relative bracket of the choking mass flow (``"chamber"`` mode).

    Returns
    -------
    NozzleKineticsResult
        Exit state, vacuum specific impulse, histories and diagnostics.

    Raises
    ------
    InputError
        If a reactant has elements outside the mechanism or ``start`` is
        unknown.
    DomainError
        If inputs are out of range or the flow leaves the thermodynamic data
        range.
    ConvergenceError
        If the mass-flow search or the integration fails; codes
        ``nozzle_mass_flow_not_converged``, ``nozzle_transonic_unresolved``
        (``"chamber"`` mode could not continue past the frozen sonic point)
        and ``nozzle_rechoked``.

    Notes
    -----
    Marched state ``(sigma_1..sigma_N, s)``; ``d sigma_i/dx = omega_i/G`` with
    ``G = mdot/A(x)`` and ``ds/dx = -(1/T) sum mu_i d sigma_i/dx`` (RP-1328
    eqs. (2.33), (2.57), (2.58), (2.60), (2.61), (2.65); RP-1311 eqs.
    (2.16)--(2.17)). The temperature at each ``x`` solves
    ``rho(T) V(T) = G`` on the frozen isentrope of the current composition with
    ``V = sqrt(2 (h0 - h))``. The Jacobian is analytic (implicit
    differentiation through this closure).

    ``"chamber"`` mode: the mass flow is the largest value for which the
    subsonic march survives, found from the shifting-equilibrium choking flux
    by a fixed-point step on ``min A(x) G_max(x)`` followed by bisection; the
    march switches to the supersonic branch where the smallest choking flux
    failed. Because the frozen sonic point is a saddle of these equations,
    the supersonic continuation is reliable only when the chemistry near the
    throat is slow (or frozen); within ``0.25 r_t`` of the switch a closure
    whose flux exceeds the local frozen maximum by less than ``1e-4``
    (relative) uses the sonic state, and a larger excess raises
    ``nozzle_transonic_unresolved``. Zero rates reproduce expansion frozen at
    the chamber composition.

    ``"equilibrium_throat"`` mode: mass flow and start state come from the
    shifting-equilibrium isentrope over the mechanism species (golden-section
    maximum of ``rho V`` for the throat). Very fast rates reproduce shifting
    equilibrium; zero rates give expansion frozen at the start station.

    References
    ----------
    .. [1] NASA-RP-1328. K. Radhakrishnan, NASA RP-1328, 1994, section 2.4.
           https://ntrs.nasa.gov/citations/19950008642
    .. [2] GRI-MECH-3.0. G. P. Smith et al., "GRI-Mech 3.0".
           http://combustion.berkeley.edu/gri-mech/version30/text30.html
    """

    if not isfinite(rate_multiplier) or rate_multiplier < 0.0:
        raise DomainError("Rate multiplier must be finite and non-negative.",
                          field="rate_multiplier")
    if not isfinite(rtol) or not 0.0 < rtol < 1.0e-2:
        raise DomainError("rtol must be in (0, 1e-2).", field="rtol")
    if start not in ("equilibrium_throat", "chamber"):
        raise InputError("start must be 'equilibrium_throat' or 'chamber'.", field="start")
    if not isfinite(start_area_ratio) or not 1.0 < start_area_ratio < nozzle.exit_area_ratio:
        raise DomainError("start_area_ratio must lie between 1 and the exit area ratio.",
                          field="start_area_ratio")
    if not isfinite(mass_flow_tolerance) or not 0.0 < mass_flow_tolerance < 1.0e-2:
        raise DomainError("mass_flow_tolerance must be in (0, 1e-2).",
                          field="mass_flow_tolerance")
    active_mechanism = mechanism or gri30_h2o2_mechanism()
    expansion = _Expansion(reactants, chamber_pressure_pa, active_mechanism,
                           provider or LocalEquilibriumProvider())
    expansion.rate_multiplier = rate_multiplier
    size = expansion.size
    radius = nozzle.throat_radius_m
    atol = [1.0e-12] * size + [rtol]
    if start == "chamber":
        flux, x_switch, y_switch, stations, iterations, clamp_limit = _transonic_from_chamber(
            expansion, nozzle, chamber_pressure_pa, rtol, atol, mass_flow_tolerance)
        accepted = len(stations) - 1
    else:
        path = _EquilibriumPath(expansion)
        throat_pressure, flux = path.throat(chamber_pressure_pa)
        log_start = path.log_pressure_at_area(log(throat_pressure), flux, start_area_ratio)
        temperature, density, velocity, sigma = path.isentrope(exp(log_start))
        x_switch = nozzle.position(start_area_ratio, divergent=True)
        y_switch = tuple(sigma) + (expansion.s0,)
        expansion.last_temperature.clear()
        first = expansion.station(x_switch, start_area_ratio, list(y_switch),
                                  flux / start_area_ratio, "sup")
        if first.frozen_mach <= 1.0:
            raise DomainError("Frozen Mach number at start_area_ratio is not supersonic.",
                              field="start_area_ratio",
                              details={"frozen_mach": first.frozen_mach})
        stations = [first]
        iterations = 0
        accepted = 0
        clamp_limit = x_switch - 1.0

    x_exit = nozzle.exit_x_m
    supersonic_points: list[tuple[float, tuple[float, ...]]] = []
    expansion.last_temperature.clear()

    def rhs_supersonic(x, y):
        area = nozzle.area_ratio(x)
        try:
            return expansion.derivatives(y[:size], y[size], flux / area, "sup",
                                         clamp=x <= clamp_limit)[0]
        except _Choked:
            if x <= clamp_limit:
                raise
            raise ConvergenceError(
                "Supersonic flow re-choked." if start == "equilibrium_throat" else
                "Supersonic continuation from the frozen sonic point failed; use "
                "start='equilibrium_throat' when the throat chemistry is fast.",
                code="nozzle_rechoked" if start == "equilibrium_throat"
                else "nozzle_transonic_unresolved", details={"x_m": x}) from None

    def jacobian_supersonic(x, y):
        return expansion.jacobian(y[:size], y[size], flux / nozzle.area_ratio(x), "sup",
                                  clamp=True)

    try:
        solution = solve_stiff(rhs_supersonic, (x_switch, x_exit), list(y_switch), rtol=rtol,
                               atol=atol, jacobian=jacobian_supersonic,
                               first_step=1.0e-3 * radius, max_step=0.5 * radius,
                               max_steps=20_000,
                               step_callback=lambda x, y: supersonic_points.append((x, y)))
    except ConvergenceError as error:
        if start == "chamber" and error.code.startswith("stiff_"):
            raise ConvergenceError(
                "Supersonic continuation from the frozen sonic point failed; use "
                "start='equilibrium_throat' when the throat chemistry is fast.",
                code="nozzle_transonic_unresolved", details=error.details) from error
        raise
    accepted += solution.accepted_steps
    for x, y in supersonic_points:
        stations.append(expansion.station(x, nozzle.area_ratio(x), list(y),
                                          flux / nozzle.area_ratio(x), "sup",
                                          clamp=x <= clamp_limit))
    exit_station = stations[-1]
    y_exit = supersonic_points[-1][1]
    vacuum = (exit_station.velocity_m_s
              + exit_station.pressure_pa * nozzle.exit_area_ratio / flux) / STANDARD_GRAVITY_M_S2
    warnings = [
        WarningMessage(
            code="finite_rate_nozzle_assumptions",
            message=("Steady inviscid adiabatic quasi-one-dimensional flow from an "
                     "infinite-area equilibrium chamber; no friction, heat loss, "
                     "divergence or condensed phases."),
            severity="info"),
        WarningMessage(
            code="mechanism_extrapolated",
            message=("GRI-Mech 3.0 was optimized for natural-gas combustion; rocket "
                     "chamber temperatures and pressures lie outside its targets."),
            severity="warning"),
    ]
    if start == "equilibrium_throat":
        warnings.append(WarningMessage(
            code="equilibrium_to_start_station",
            message=(f"Shifting equilibrium assumed up to A/A_t = {start_area_ratio:g}; "
                     "finite rates only downstream."), severity="info"))
    return NozzleKineticsResult(
        mechanism=active_mechanism.name,
        rate_multiplier=rate_multiplier,
        start=start,
        chamber_temperature_k=expansion.chamber_temperature_k,
        chamber_pressure_pa=chamber_pressure_pa,
        mass_flux_throat_kg_m2_s=flux,
        characteristic_velocity_m_s=chamber_pressure_pa / flux,
        start_x_m=x_switch,
        exit=exit_station,
        vacuum_specific_impulse_s=vacuum,
        stations=tuple(stations),
        mass_flow_iterations=iterations,
        accepted_steps=accepted,
        maximum_element_residual=expansion.element_residual(y_exit[:size]),
        warnings=tuple(warnings),
    )


def _transonic_from_chamber(expansion, nozzle, chamber_pressure_pa, rtol, atol,
                            mass_flow_tolerance):
    """Shoot the choking mass flow from the chamber composition (``"chamber"`` mode)."""
    size = expansion.size
    radius = nozzle.throat_radius_m
    x_start = nozzle.start_x_m
    x_window = min(nozzle.position(1.5, divergent=True), nozzle.exit_x_m)
    y0 = list(expansion.sigma0) + [expansion.s0]

    def subsonic_march(flux_throat, stop_x=None):
        samples: list[tuple[float, float, tuple[float, ...]]] = []
        failures: list[float] = []
        expansion.last_temperature.clear()

        def rhs(x, y):
            area = nozzle.area_ratio(x)
            try:
                return expansion.derivatives(y[:size], y[size], flux_throat / area, "sub")[0]
            except _Choked as choked:
                failures.append(area * choked.local_maximum)
                raise

        def jacobian(x, y):
            area = nozzle.area_ratio(x)
            return expansion.jacobian(y[:size], y[size], flux_throat / area, "sub", clamp=True)

        def record(x, y):
            value = nozzle.area_ratio(x) * expansion.sonic(y[:size], y[size])[1]
            samples.append((x, value, tuple(y)))
            best_x, best = min(((sx, sf) for sx, sf, _ in samples), key=lambda item: item[1])
            if stop_x is None and x > 0.0 and x > best_x + 0.2 * radius \
                    and value > best * (1.0 + 2.0e-3):
                raise _MarchComplete

        record(x_start, tuple(y0))
        end = x_window if stop_x is None else stop_x
        targets = (0.0, end) if end > 0.0 else (end,)
        try:
            solve_stiff(rhs, (x_start, end), y0, rtol=rtol, atol=atol, jacobian=jacobian,
                        t_eval=targets, first_step=1.0e-3 * radius,
                        max_step=0.05 * radius, step_callback=record)
        except _MarchComplete:
            pass
        except (ConvergenceError, _Choked):
            if not failures:
                raise
            return False, samples
        return True, samples

    expansion.last_temperature.clear()
    frozen_flux = expansion.sonic(expansion.sigma0, expansion.s0)[1]
    if expansion.rate_multiplier == 0.0:
        flux = frozen_flux * (1.0 - 1.0e-3)
    else:
        flux = _EquilibriumPath(expansion).throat(chamber_pressure_pa, 1.0e-3)[1] \
            * (1.0 - 1.0e-3)
    low, high = 0.0, float("inf")
    best_samples: list = []
    failure_x = None
    iterations = 0
    converged = False
    for iterations in range(1, 41):
        survived, samples = subsonic_march(flux)
        if survived:
            low, best_samples = flux, samples
            index = min(range(len(samples)), key=lambda i: samples[i][1])
            estimate = _refined_minimum(samples, index)
            if estimate - flux <= mass_flow_tolerance * flux \
                    or high - low <= mass_flow_tolerance * low:
                converged = True
                break
            proposal = estimate * (1.0 - 0.5 * mass_flow_tolerance)
        else:
            high = flux
            failure_x = samples[-1][0]
            if low > 0.0 and high - low <= mass_flow_tolerance * low:
                converged = True
                break
            proposal = 0.5 * (low + high) if low > 0.0 else 0.98 * flux
        if not low < proposal < high:
            proposal = 0.5 * (low + high)
        flux = proposal
    if not converged:
        raise ConvergenceError("Choking mass flow did not converge.",
                               code="nozzle_mass_flow_not_converged",
                               details={"low": low, "high": high})
    flux = low
    minimum_index = min(range(len(best_samples)), key=lambda i: best_samples[i][1])
    if failure_x is not None:
        _, best_samples = subsonic_march(flux, stop_x=failure_x)
        minimum_index = len(best_samples) - 1
    x_switch, _, y_switch = best_samples[minimum_index]
    expansion.last_temperature.clear()
    stations = [expansion.station(x, nozzle.area_ratio(x), list(y), flux / nozzle.area_ratio(x),
                                  "sub")
                for x, _, y in best_samples[:minimum_index + 1]]
    return flux, x_switch, y_switch, stations, iterations, x_switch + 0.25 * radius


def _parabola_vertex(points):
    """Return the vertex value of the parabola through three ``(x, f)`` points, or ``None``."""
    (x0, f0), (x1, f1), (x2, f2) = points
    denominator = (x0 - x1) * (x0 - x2) * (x1 - x2)
    if denominator == 0.0:
        return None
    a = (x2 * (f1 - f0) + x1 * (f0 - f2) + x0 * (f2 - f1)) / denominator
    b = (x2 * x2 * (f0 - f1) + x1 * x1 * (f2 - f0) + x0 * x0 * (f1 - f2)) / denominator
    if a <= 0.0:
        return None
    vertex = -b / (2.0 * a)
    if not min(x0, x2) <= vertex <= max(x0, x2):
        return None
    c = f0 - a * x0 * x0 - b * x0
    return a * vertex * vertex + b * vertex + c


def _refined_minimum(samples, index):
    """Refine the sampled minimum of ``f(x) = A G_max`` by parabolas.

    The area law changes curvature at the throat (``x = 0``, always a sample),
    so a parabola never spans it: away from the throat the three points around
    the minimum are used; at the throat each side is fitted separately.
    """
    points = [(x, f) for x, f, _ in samples]
    best = points[index][1]
    candidates = []
    if points[index][0] == 0.0:
        if index >= 2:
            candidates.append(points[index - 2:index + 1])
        if index + 2 < len(points):
            candidates.append(points[index:index + 3])
    elif 0 < index < len(points) - 1:
        window = points[index - 1:index + 2]
        if not window[0][0] < 0.0 < window[2][0]:
            candidates.append(window)
    for window in candidates:
        value = _parabola_vertex(window)
        if value is not None:
            best = min(best, value)
    return best


# ---------------------------------------------------------------------------
# Shifting equilibrium along the nozzle and the Bray sudden-freezing estimate
# ---------------------------------------------------------------------------


class _EquilibriumPath:
    """Shifting-equilibrium isentrope over the mechanism species (internal)."""

    def __init__(self, expansion: _Expansion) -> None:
        self.expansion = expansion
        self.cache: dict[float, tuple] = {}
        self.last_temperature = expansion.chamber_temperature_k

    def state_tp(self, temperature, pressure):
        expansion = self.expansion
        result = expansion.provider.solve(EquilibriumProblem.tp(
            expansion.reactants, temperature_k=temperature, pressure_pa=pressure,
            candidate_species=expansion.mechanism.species))
        mass = result.total_mass_kg
        moles = {item.species: item.moles for item in result.species}
        sigma = tuple(moles.get(key, 0.0) / mass for key in expansion.mechanism.species)
        return result.entropy_j_k / mass, result.enthalpy_j / mass, sigma

    def isentrope(self, pressure):
        if pressure in self.cache:
            return self.cache[pressure]
        target = self.expansion.s0
        states: dict[float, tuple] = {}

        def residual(temperature):
            states[temperature] = self.state_tp(temperature, pressure)
            return (states[temperature][0] - target) / abs(target)

        guess = self.last_temperature
        low, high = guess * 0.97, min(guess * 1.01, 6000.0)
        f_low, f_high = residual(low), residual(high)
        for _ in range(60):
            if f_low <= 0.0 <= f_high:
                break
            if f_low > 0.0:
                low = max(_MINIMUM_TEMPERATURE_K, low * 0.95)
                f_low = residual(low)
            if f_high < 0.0:
                high = min(6000.0, high * 1.05)
                f_high = residual(high)
        else:
            raise ConvergenceError("Equilibrium isentrope not bracketed.",
                                   code="nozzle_root_not_converged")
        temperature = _illinois(residual, low, high, f_low, f_high, 1.0e-12)
        state = states.get(temperature) or self.state_tp(temperature, pressure)
        _, enthalpy, sigma = state
        v2 = 2.0 * (self.expansion.h0 - enthalpy)
        if v2 <= 0.0:
            raise DomainError("Equilibrium expansion produced no enthalpy drop.",
                              code="nozzle_invalid_state")
        density = pressure / (_R * sum(sigma) * temperature)
        value = (temperature, density, sqrt(v2), sigma)
        self.cache[pressure] = value
        self.last_temperature = temperature
        return value

    def flux(self, pressure):
        temperature, density, velocity, _ = self.isentrope(pressure)
        return density * velocity

    def throat(self, chamber_pressure_pa, tolerance=1.0e-6):
        """Return ``(p_t, G_max)``: golden-section maximum of ``rho V`` on ``ln p``."""
        golden = (sqrt(5.0) - 1.0) / 2.0
        low, high = log(chamber_pressure_pa / 2.5), log(chamber_pressure_pa / 1.2)
        first, second = high - golden * (high - low), low + golden * (high - low)
        f_first, f_second = self.flux(exp(first)), self.flux(exp(second))
        while high - low > tolerance:
            if f_first > f_second:
                high, second, f_second = second, first, f_first
                first = high - golden * (high - low)
                f_first = self.flux(exp(first))
            else:
                low, first, f_first = first, second, f_second
                second = low + golden * (high - low)
                f_second = self.flux(exp(second))
        pressure = exp(0.5 * (low + high))
        return pressure, self.flux(pressure)

    def log_pressure_at_area(self, log_throat, flux, area_ratio):
        """Return ``ln p`` of the supersonic equilibrium station with ``A/A_t = area_ratio``."""

        def residual(value):
            return log(flux / self.flux(exp(value)) / area_ratio)

        high = log_throat - 1.0e-3
        low = log_throat - 1.6 * log(area_ratio) - 0.2
        while residual(low) < 0.0:
            low -= 0.5
        return _illinois(residual, low, high, residual(low), residual(high), 1.0e-12)


def bray_sudden_freezing(
    *,
    reactants: Mapping[str, float],
    chamber_pressure_pa: float,
    nozzle: HyperbolicNozzle,
    mechanism: Mechanism | None = None,
    freezing_constant: float = 1.0,
    rate_multiplier: float = 1.0,
    provider: EquilibriumProvider | None = None,
    stations: int = 60,
) -> BrayFreezingResult:
    """Estimate the exit state by Bray's sudden-freezing approximation [1]_.

    Parameters
    ----------
    reactants, chamber_pressure_pa, nozzle, mechanism, provider
        As for :func:`solve_nozzle_kinetics`.
    freezing_constant : float, default 1.0
        ``K`` of the criterion; Bray reports experimental values from 0.42 to
        1.25 for his single-reaction gas [1]_.
    rate_multiplier : float, default 1.0
        Forward-rate factor applied to the dissociation rate ``r_D``.
    stations : int, default 60
        Pressure levels between the throat and the exit (log-spaced).

    Returns
    -------
    BrayFreezingResult
        Freezing point, sudden-freezing and shifting-equilibrium vacuum
        specific impulse.

    Raises
    ------
    DomainError
        If inputs are out of range.
    ConvergenceError
        If an equilibrium isentrope point does not converge.

    Notes
    -----
    Bray's criterion compares the net rate ``-(d alpha/dt)_e`` required to
    stay in equilibrium with the one-way dissociation rate ``r_De``, both on
    the equilibrium solution (thesis eq. (36)) [1]_. For a multi-reaction gas
    this implementation takes ``alpha`` as the total mol per kg ``sigma_m``,
    which only reactions with ``dn != 0`` change, and
    ``r_De = (1/rho) sum_j |dn_j| R_j,diss`` with ``R_j,diss`` the one-way rate
    in the mole-increasing direction. This generalization is the present
    model's choice, not Bray's. The equilibrium isentrope uses the
    repository's equilibrium provider over the mechanism species; the throat
    maximizes ``rho V`` (golden section on ``ln p``), the divergent stations
    follow from ``A = mdot/(rho V)``, and ``d sigma_m/dt = V d sigma_m/dx`` by
    central differences. After the freezing point the composition is frozen
    and the expansion isentropic (Bray's ``Phi = 0`` solution, eq. (30)).

    References
    ----------
    .. [1] BRAY-1960-THESIS. K. N. C. Bray, "Atomic recombination in a
           hypersonic wind-tunnel nozzle", PhD thesis, University of
           Southampton, 1960, section 5, eqs. (34)-(36), pp. 35-41.
           https://eprints.soton.ac.uk/458375/1/38084.pdf
    """

    if not isfinite(freezing_constant) or freezing_constant <= 0.0:
        raise DomainError("Freezing constant must be positive.", field="freezing_constant")
    if stations < 10:
        raise DomainError("At least 10 stations are required.", field="stations")
    active_mechanism = mechanism or gri30_h2o2_mechanism()
    expansion = _Expansion(reactants, chamber_pressure_pa, active_mechanism,
                           provider or LocalEquilibriumProvider())
    path = _EquilibriumPath(expansion)
    throat_pressure, flux = path.throat(chamber_pressure_pa)

    log_throat = log(throat_pressure)
    log_exit = path.log_pressure_at_area(log_throat, flux, nozzle.exit_area_ratio)
    levels = [log_throat - 1.0e-3 + (log_exit - log_throat + 1.0e-3) * i / (stations - 1)
              for i in range(stations)]
    compiled = expansion.compiled
    rows = []
    for level in levels:
        temperature, density, velocity, sigma = path.isentrope(exp(level))
        area = flux / (density * velocity)
        x = nozzle.position(max(area, 1.0), divergent=True)
        concentrations = [density * value for value in sigma]
        _, h_rt, s_r = compiled.thermo(temperature)
        constants = compiled.equilibrium_constants(temperature, h_rt, s_r)
        _, forward, reverse, _ = compiled.rates(temperature, concentrations, rate_multiplier,
                                                constants)
        dissociation = 0.0
        for j, change in enumerate(compiled.delta_moles):
            if change > 0:
                dissociation += change * forward[j]
            elif change < 0:
                dissociation += -change * reverse[j]
        rows.append((x, area, velocity, sum(sigma), dissociation / density, temperature,
                     exp(level)))
    criterion = []
    freezing_index = None
    for i in range(1, len(rows) - 1):
        (x0, _, _, m0, _, _, _), (_, area1, v1, _, r1, _, _), (x2, _, _, m2, _, _, _) = \
            rows[i - 1], rows[i], rows[i + 1]
        rate = -v1 * (m2 - m0) / (x2 - x0)
        ratio = rate / r1 if r1 > 0.0 else float("inf")
        criterion.append((area1, ratio))
        if freezing_index is None and ratio >= freezing_constant:
            freezing_index = i
    warnings = [WarningMessage(
        code="bray_generalized_criterion",
        message=("Bray's single-reaction criterion applied to total moles per kg; "
                 "the generalization is a model choice."), severity="info")]
    exit_temperature, exit_density, exit_velocity, _ = path.isentrope(exp(log_exit))
    equilibrium_isp = (exit_velocity + exp(log_exit) * nozzle.exit_area_ratio / flux) \
        / STANDARD_GRAVITY_M_S2
    if freezing_index is None:
        warnings.append(WarningMessage(code="bray_no_freezing",
                                       message="Criterion not met before the exit.",
                                       severity="info"))
        x_f, area_f, temperature_f, pressure_f = (rows[-1][0], rows[-1][1], rows[-1][5],
                                                  rows[-1][6])
        isp = equilibrium_isp
    else:
        # criterion[i - 1] belongs to rows[i]; interpolate ln(ratio) between the
        # last station below K and the first at or above it.
        level = levels[freezing_index]
        if freezing_index >= 2:
            ratio_a = criterion[freezing_index - 2][1]
            ratio_b = criterion[freezing_index - 1][1]
            if 0.0 < ratio_a < freezing_constant <= ratio_b:
                weight = (log(freezing_constant) - log(ratio_a)) / (log(ratio_b) - log(ratio_a))
                level = levels[freezing_index - 1] + weight * (
                    levels[freezing_index] - levels[freezing_index - 1])
        pressure_f = exp(level)
        temperature_f, density_f, velocity_f, sigma_f = path.isentrope(pressure_f)
        area_f = flux / (density_f * velocity_f)
        x_f = nozzle.position(max(area_f, 1.0), divergent=True)
        expansion.last_temperature.clear()
        state = expansion.closure(sigma_f, expansion.s0, flux / nozzle.exit_area_ratio, "sup")
        isp = (state.velocity_m_s + state.pressure_pa * nozzle.exit_area_ratio / flux) \
            / STANDARD_GRAVITY_M_S2
    del exit_temperature, exit_density
    return BrayFreezingResult(
        freezing_constant=freezing_constant,
        freezing_x_m=x_f,
        freezing_area_ratio=area_f,
        freezing_temperature_k=temperature_f,
        freezing_pressure_pa=pressure_f,
        equilibrium_mass_flux_throat_kg_m2_s=flux,
        vacuum_specific_impulse_s=isp,
        equilibrium_vacuum_specific_impulse_s=equilibrium_isp,
        criterion=tuple(criterion),
        warnings=tuple(warnings),
    )
