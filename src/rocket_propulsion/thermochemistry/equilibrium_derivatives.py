"""Frozen and equilibrium (shifting) thermodynamic derivatives of an ideal-gas mixture.

For a gas mixture in chemical equilibrium the derivatives
``(dln V/dln T)_p``, ``(dln V/dln p)_T`` and the equilibrium specific heat
``c_p,e`` follow from the linearized equilibrium conditions, RP-1311 Part I
section 2.5, equations (2.50)--(2.66) [1]_. The remaining derivatives
``c_v``, ``gamma = c_p/c_v``, ``gamma_s = (dln p/dln rho)_s`` and the sound
speed follow from section 2.6, equations (2.67)--(2.74). Species properties
are the NASA Glenn polynomials of :mod:`.species` [2]_.

Only gaseous species are treated (the local equilibrium provider has no
condensed phases), so the sums over ``j = NG+1, ..., NS`` of RP-1311 vanish.

References
----------
.. [1] NASA RP-1311, "Computer Program for Calculation of Complex Chemical
       Equilibrium Compositions and Applications, Part I: Analysis", sections
       2.5 and 2.6. https://ntrs.nasa.gov/citations/19950013764
.. [2] NASA/TP-2002-211556. NASA Glenn, "NASA Glenn Coefficients for
       Calculating Thermodynamic Properties of Individual Species".
       https://ntrs.nasa.gov/citations/20020085330
"""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass
from math import exp, isfinite, log, sqrt

from rocket_propulsion.core.errors import ConvergenceError, DomainError, InputError

from .equilibrium import EquilibriumProblem, EquilibriumProvider, EquilibriumResult
from .nasa_polynomials import UNIVERSAL_GAS_CONSTANT_J_MOL_K as _R
from .species import get_species


@dataclass(frozen=True, slots=True)
class EquilibriumDerivatives:
    """Frozen and equilibrium derivatives at one state [1]_.

    Attributes
    ----------
    temperature_k, pressure_pa : float
        State.
    molecular_mass_kg_mol : float
        Mixture molar mass ``m/n``.
    frozen_cp_j_kg_k, frozen_cv_j_kg_k, frozen_gamma : float
        Fixed-composition specific heats (equation (2.49b)) and their ratio.
    frozen_sound_speed_m_s : float
        ``sqrt(gamma_fr n R T / m)``.
    reaction_cp_j_kg_k : float
        Reaction contribution of equation (2.49c).
    equilibrium_cp_j_kg_k : float
        ``c_p,e`` of equation (2.59).
    dlnv_dlnt_p, dlnv_dlnp_t : float
        Equations (2.50) and (2.51).
    equilibrium_cv_j_kg_k : float
        Equation (2.70).
    equilibrium_gamma : float
        ``c_p/c_v``, equation (2.72).
    gamma_s : float
        Isentropic exponent ``(dln p/dln rho)_s``, equation (2.73).
    equilibrium_sound_speed_m_s : float
        Equation (2.74).
    dlnp_dlnt_s : float
        Equation (2.75).
    species : tuple of str
        Species included.

    References
    ----------
    .. [1] NASA RP-1311, "Computer Program for Calculation of Complex Chemical
           Equilibrium Compositions and Applications, Part I: Analysis",
           sections 2.5 and 2.6. https://ntrs.nasa.gov/citations/19950013764
    """

    temperature_k: float
    pressure_pa: float
    molecular_mass_kg_mol: float
    frozen_cp_j_kg_k: float
    frozen_cv_j_kg_k: float
    frozen_gamma: float
    frozen_sound_speed_m_s: float
    reaction_cp_j_kg_k: float
    equilibrium_cp_j_kg_k: float
    dlnv_dlnt_p: float
    dlnv_dlnp_t: float
    equilibrium_cv_j_kg_k: float
    equilibrium_gamma: float
    gamma_s: float
    equilibrium_sound_speed_m_s: float
    dlnp_dlnt_s: float
    species: tuple[str, ...]


def _solve(matrix: list[list[float]], vector: list[float]) -> list[float]:
    """Gaussian elimination with partial pivoting; raises on a singular matrix."""
    size = len(vector)
    rows = [row[:] + [value] for row, value in zip(matrix, vector, strict=True)]
    scale = max(abs(value) for row in matrix for value in row)
    for column in range(size):
        pivot = max(range(column, size), key=lambda index: abs(rows[index][column]))
        if abs(rows[pivot][column]) <= 1.0e-14 * scale:
            raise DomainError(
                "Singular derivative matrix (linearly dependent element basis).",
                code="derivative_singular_matrix",
            )
        rows[column], rows[pivot] = rows[pivot], rows[column]
        for index in range(column + 1, size):
            factor = rows[index][column] / rows[column][column]
            if factor:
                for k in range(column, size + 1):
                    rows[index][k] -= factor * rows[column][k]
    solution = [0.0] * size
    for index in range(size - 1, -1, -1):
        total = rows[index][size] - sum(rows[index][k] * solution[k] for k in range(index + 1, size))
        solution[index] = total / rows[index][index]
    return solution


def derivatives_from_composition(
    moles: Mapping[str, float], temperature_k: float, pressure_pa: float
) -> EquilibriumDerivatives:
    """Frozen and equilibrium derivatives of a gas mixture assumed in equilibrium [1]_.

    Parameters
    ----------
    moles : mapping of str to float
        Amount of each gaseous species (species-catalog key) in any consistent
        unit (mol, mol/kg); only ratios and the mass basis matter.
    temperature_k : float
        Temperature in K.
    pressure_pa : float
        Pressure in Pa.

    Returns
    -------
    EquilibriumDerivatives
        Derivatives per unit mass.

    Raises
    ------
    InputError
        If the state or an amount is invalid, or no species is positive.
    DomainError
        If a species polynomial does not cover the temperature or the element
        matrix is singular.

    Notes
    -----
    The composition must be the equilibrium composition at the state; the
    derivatives describe how it shifts. With ``h_j = H_j/(R T)`` and gas
    amounts ``n_j``, the temperature derivatives ``(d pi_i/d ln T)_p`` and
    ``(d ln n/d ln T)_p`` solve equations (2.56) and (2.58), and the pressure
    derivatives solve equations (2.64) and (2.66) (no condensed species).
    Then ``c_p,e`` is equation (2.59), ``(dln V/dln T)_p = 1 + (dln n/dln T)_p``
    (2.50), ``(dln V/dln p)_T = -1 + (dln n/dln p)_T`` (2.51), ``c_v`` is
    (2.70) with ``pV/T = n R``, ``gamma_s = -gamma/(dln V/dln p)_T`` (2.73),
    the sound speed ``a = sqrt(n R T gamma_s)`` per unit mass (2.74) and
    ``(dln p/dln T)_s = c_p/(n R (dln V/dln T)_p)`` (2.75). The frozen values
    keep composition fixed: ``(dln V/dln T)_p = 1``, ``(dln V/dln p)_T = -1``.

    References
    ----------
    .. [1] NASA RP-1311, "Computer Program for Calculation of Complex Chemical
           Equilibrium Compositions and Applications, Part I: Analysis",
           equations (2.49)--(2.75). https://ntrs.nasa.gov/citations/19950013764
    """
    if not isfinite(temperature_k) or temperature_k <= 0.0:
        raise InputError("Temperature must be finite and positive.", field="temperature_k")
    if not isfinite(pressure_pa) or pressure_pa <= 0.0:
        raise InputError("Pressure must be finite and positive.", field="pressure_pa")
    amounts: dict[str, float] = {}
    for key, value in moles.items():
        value = float(value)
        if not isfinite(value) or value < 0.0:
            raise InputError(f"Amount of {key} must be finite and non-negative.", field="moles")
        if value > 0.0:
            amounts[str(key)] = amounts.get(str(key), 0.0) + value
    if not amounts:
        raise InputError("At least one species amount must be positive.", field="moles")
    names = tuple(amounts)
    species = [get_species(name) for name in names]
    n = [amounts[name] for name in names]
    total_moles = sum(n)
    mass = sum(amount * item.molecular_mass_kg_mol for amount, item in zip(n, species, strict=True))
    rt = _R * temperature_k
    properties = [item.properties(temperature_k) for item in species]
    h = [prop.enthalpy_molar_j_mol / rt for prop in properties]
    cp_over_r = [prop.cp_molar_j_mol_k / _R for prop in properties]
    formulas = [dict(item.formula) for item in species]
    elements = sorted({element for formula in formulas for element in formula})
    a = [[formula.get(element, 0.0) for formula in formulas] for element in elements]
    count = len(elements)

    # Left-hand side shared by equations (2.56)/(2.58) and (2.64)/(2.66).
    matrix = [[0.0] * (count + 1) for _ in range(count + 1)]
    for k in range(count):
        for i in range(count):
            matrix[k][i] = sum(a[k][j] * a[i][j] * n[j] for j in range(len(n)))
        matrix[k][count] = sum(a[k][j] * n[j] for j in range(len(n)))
    for i in range(count):
        matrix[count][i] = sum(a[i][j] * n[j] for j in range(len(n)))
    rhs_temperature = [-sum(a[k][j] * n[j] * h[j] for j in range(len(n))) for k in range(count)]
    rhs_temperature.append(-sum(n[j] * h[j] for j in range(len(n))))
    rhs_pressure = [sum(a[k][j] * n[j] for j in range(len(n))) for k in range(count)]
    rhs_pressure.append(total_moles)
    temperature_solution = _solve(matrix, rhs_temperature)
    pressure_solution = _solve(matrix, rhs_pressure)
    dpi_dlnt = temperature_solution[:count]
    dlnn_dlnt = temperature_solution[count]
    dlnn_dlnp = pressure_solution[count]

    # Equation (2.59): c_p,e / R on the n_j basis.
    frozen_cp_over_r = sum(nj * cp for nj, cp in zip(n, cp_over_r, strict=True))
    reaction_cp_over_r = (
        sum(sum(a[i][j] * n[j] * h[j] for j in range(len(n))) * dpi_dlnt[i] for i in range(count))
        + sum(nj * hj for nj, hj in zip(n, h, strict=True)) * dlnn_dlnt
        + sum(nj * hj * hj for nj, hj in zip(n, h, strict=True))
    )
    gas_constant = _R * total_moles / mass
    frozen_cp = _R * frozen_cp_over_r / mass
    reaction_cp = _R * reaction_cp_over_r / mass
    equilibrium_cp = frozen_cp + reaction_cp
    dlnv_dlnt = 1.0 + dlnn_dlnt
    dlnv_dlnp = -1.0 + dlnn_dlnp
    if dlnv_dlnp >= 0.0 or equilibrium_cp <= 0.0:
        raise DomainError("Derivatives violate thermodynamic stability.",
                          code="derivative_unstable_state")
    equilibrium_cv = equilibrium_cp + gas_constant * dlnv_dlnt ** 2 / dlnv_dlnp
    frozen_cv = frozen_cp - gas_constant
    if equilibrium_cv <= 0.0 or frozen_cv <= 0.0:
        raise DomainError("Non-positive specific heat at constant volume.",
                          code="derivative_unstable_state")
    equilibrium_gamma = equilibrium_cp / equilibrium_cv
    gamma_s = -equilibrium_gamma / dlnv_dlnp
    frozen_gamma = frozen_cp / frozen_cv
    return EquilibriumDerivatives(
        temperature_k=temperature_k,
        pressure_pa=pressure_pa,
        molecular_mass_kg_mol=mass / total_moles,
        frozen_cp_j_kg_k=frozen_cp,
        frozen_cv_j_kg_k=frozen_cv,
        frozen_gamma=frozen_gamma,
        frozen_sound_speed_m_s=sqrt(frozen_gamma * gas_constant * temperature_k),
        reaction_cp_j_kg_k=reaction_cp,
        equilibrium_cp_j_kg_k=equilibrium_cp,
        dlnv_dlnt_p=dlnv_dlnt,
        dlnv_dlnp_t=dlnv_dlnp,
        equilibrium_cv_j_kg_k=equilibrium_cv,
        equilibrium_gamma=equilibrium_gamma,
        gamma_s=gamma_s,
        equilibrium_sound_speed_m_s=sqrt(gamma_s * gas_constant * temperature_k),
        dlnp_dlnt_s=equilibrium_cp / (gas_constant * dlnv_dlnt),
        species=names,
    )


def equilibrium_derivatives(result: EquilibriumResult) -> EquilibriumDerivatives:
    """Derivatives of an :class:`.EquilibriumResult` from a gas-phase provider [1]_.

    Parameters
    ----------
    result : EquilibriumResult
        Equilibrium state whose species keys are in the bundled species
        catalog (the local Gibbs provider's results satisfy this).

    Returns
    -------
    EquilibriumDerivatives
        Frozen and shifting-equilibrium derivatives at the result's state.

    Raises
    ------
    DomainError, InputError
        As :func:`derivatives_from_composition`.

    References
    ----------
    .. [1] NASA RP-1311, "Computer Program for Calculation of Complex Chemical
           Equilibrium Compositions and Applications, Part I: Analysis",
           sections 2.5 and 2.6. https://ntrs.nasa.gov/citations/19950013764
    """
    return derivatives_from_composition(
        {amount.species: amount.moles for amount in result.species},
        result.temperature_k,
        result.pressure_pa,
    )


@dataclass(frozen=True, slots=True)
class FiniteDifferenceDerivatives:
    """Derivatives estimated by re-solving the equilibrium at perturbed states.

    Attributes
    ----------
    equilibrium_cp_j_kg_k : float
        ``(h(T+dT) - h(T-dT)) / (2 dT)`` at constant pressure.
    dlnv_dlnt_p, dlnv_dlnp_t : float
        Central differences of ``ln(n T / p)``.
    equilibrium_sound_speed_m_s : float
        ``sqrt(dp/drho)`` with both states on the isentrope through the result.
    relative_step : float
        Relative perturbation used.

    References
    ----------
    .. [1] NASA RP-1311, "Computer Program for Calculation of Complex Chemical
           Equilibrium Compositions and Applications, Part I: Analysis",
           equations (2.50), (2.51) and (2.67). https://ntrs.nasa.gov/citations/19950013764
    """

    equilibrium_cp_j_kg_k: float
    dlnv_dlnt_p: float
    dlnv_dlnp_t: float
    equilibrium_sound_speed_m_s: float
    relative_step: float


def finite_difference_derivatives(
    result: EquilibriumResult,
    provider: EquilibriumProvider,
    *,
    relative_step: float = 1.0e-4,
    entropy_tolerance: float = 1.0e-12,
) -> FiniteDifferenceDerivatives:
    """Estimate equilibrium derivatives numerically with ``provider`` (a check) [1]_.

    Parameters
    ----------
    result : EquilibriumResult
        Reference state.
    provider : EquilibriumProvider
        Provider re-solved at perturbed temperatures and pressures with the
        same element inventory (assigned-T-and-p problems whose reactants and
        candidate products are the result's retained species).
    relative_step : float, default 1e-4
        Relative step in ``T`` and ``p``.
    entropy_tolerance : float, default 1e-12
        Relative entropy residual accepted on the isentrope.

    Returns
    -------
    FiniteDifferenceDerivatives
        Central-difference estimates, independent of the linearized equations
        of :func:`derivatives_from_composition`.

    Raises
    ------
    InputError
        If the step is not in ``(0, 0.1]``.
    ConvergenceError
        If the isentropic temperature is not found within 60 secant steps.

    Notes
    -----
    The sound speed uses the definition ``a**2 = (dp/drho)_s``, equation
    (2.67) of [1]_: at ``p (1 +/- step)`` the temperature with the reference
    mixture entropy is found by a secant iteration on assigned-T-and-p
    solutions, and ``rho = p M / (R T)``.

    References
    ----------
    .. [1] NASA RP-1311, "Computer Program for Calculation of Complex Chemical
           Equilibrium Compositions and Applications, Part I: Analysis",
           equations (2.50), (2.51) and (2.67). https://ntrs.nasa.gov/citations/19950013764
    """
    if not isfinite(relative_step) or not 0.0 < relative_step <= 0.1:
        raise InputError("relative_step must be in (0, 0.1].", field="relative_step")
    reactants = {amount.species: amount.moles for amount in result.species}
    candidates = tuple(amount.species for amount in result.species)
    mass = result.total_mass_kg

    def solve(temperature: float, pressure: float) -> EquilibriumResult:
        return provider.solve(EquilibriumProblem.tp(
            reactants, temperature_k=temperature, pressure_pa=pressure,
            candidate_species=candidates,
        ))

    t0, p0 = result.temperature_k, result.pressure_pa
    dt, dp = t0 * relative_step, p0 * relative_step
    hot, cold = solve(t0 + dt, p0), solve(t0 - dt, p0)
    high, low = solve(t0, p0 + dp), solve(t0, p0 - dp)
    cp = (hot.enthalpy_j - cold.enthalpy_j) / (2.0 * dt * mass)

    def log_volume(state: EquilibriumResult) -> float:
        return log(state.total_moles * state.temperature_k / state.pressure_pa)

    dlnv_dlnt = (log_volume(hot) - log_volume(cold)) / (log(t0 + dt) - log(t0 - dt))
    dlnv_dlnp = (log_volume(high) - log_volume(low)) / (log(p0 + dp) - log(p0 - dp))
    entropy = result.entropy_j_k

    def isentropic_density(pressure: float) -> float:
        # Secant iteration in ln T on the entropy residual.
        guess = t0 * exp(log(pressure / p0) * (1.0 / dlnp_dlnt_guess))
        x0, x1 = log(guess), log(guess) + 1.0e-6
        f0 = solve(exp(x0), pressure).entropy_j_k - entropy
        for _ in range(60):
            state = solve(exp(x1), pressure)
            f1 = state.entropy_j_k - entropy
            if abs(f1) <= entropy_tolerance * abs(entropy):
                return pressure * (state.total_mass_kg / state.total_moles) / (
                    _R * state.temperature_k
                )
            if f1 == f0:
                break
            x0, x1, f0 = x1, x1 - f1 * (x1 - x0) / (f1 - f0), f1
        raise ConvergenceError("Isentropic temperature not found.", code="isentrope_not_converged")

    # (dln p/dln T)_s from the difference estimates, only to start the secant.
    dlnp_dlnt_guess = cp * mass / (_R * result.total_moles * dlnv_dlnt)
    rho_high = isentropic_density(p0 + dp)
    rho_low = isentropic_density(p0 - dp)
    sound_speed = sqrt(2.0 * dp / (rho_high - rho_low))
    return FiniteDifferenceDerivatives(
        equilibrium_cp_j_kg_k=cp,
        dlnv_dlnt_p=dlnv_dlnt,
        dlnv_dlnp_t=dlnv_dlnp,
        equilibrium_sound_speed_m_s=sound_speed,
        relative_step=relative_step,
    )
