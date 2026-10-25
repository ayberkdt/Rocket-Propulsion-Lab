"""Experimental equilibrium-chamber and nozzle-expansion comparisons.

A chamber state from an HP equilibrium solve [1]_ is expanded to the same
exit pressure twice: with frozen composition and constant properties
(ideal-nozzle relations [2]_), and isentropically with shifting equilibrium
composition. Both use the frozen-chamber throat mass flow. The results are
screening values and are not validated against NASA CEA or engine data.

References
----------
.. [1] NASA RP-1311, "Computer Program for Calculation of Complex Chemical
       Equilibrium Compositions and Applications, Part I: Analysis".
       https://ntrs.nasa.gov/citations/19950013764
.. [2] SUTTON-BIBLARZ-2017. G. P. Sutton and O. Biblarz, "Rocket Propulsion
       Elements", 9th ed., Wiley, 2017.
"""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass
from math import isfinite, sqrt

from rocket_propulsion.compressible.nozzle import (
    STANDARD_GRAVITY_M_S2,
    NozzleResult,
    calculate_ideal_nozzle,
)
from rocket_propulsion.core.errors import DomainError, InputError
from rocket_propulsion.core.metadata import WarningMessage

from .equilibrium import (
    EquilibriumProblem,
    EquilibriumProvider,
    EquilibriumResult,
    LocalEquilibriumProvider,
)
from .species import get_reactant


@dataclass(frozen=True, slots=True)
class PropellantPair:
    """Source-labelled preset for local screening and optional CEA runs.

    Attributes
    ----------
    key, name : str
        Preset key and display name.
    fuel_reactant, oxidizer_reactant : str
        Bundled reactant keys; condensed reactants (``"H2(L)"``, ``"RP-1"``)
        carry the assigned storage-state enthalpy of NASA CEA's database.
    minimum_oxidizer_fuel_ratio, maximum_oxidizer_fuel_ratio : float
        Default sweep range of the oxidizer-to-fuel mass ratio.
    reference_oxidizer_fuel_ratio : float
        Reference mass ratio of the preset.
    cea_fuel_name, cea_oxidizer_name : str
        Reactant names for CEA input decks.
    source : str
        Basis of the preset values.
    verified_chamber_pressure_range_pa : tuple of float
        Chamber pressures at whose ends the local solver was compared with
        NASA CEA 3.3.4 at the reference O/F and area ratio 40
        (``tests/reference_data/cea_gate_points.json``, pressure-domain points).
    """

    key: str
    name: str
    fuel_reactant: str
    oxidizer_reactant: str
    minimum_oxidizer_fuel_ratio: float
    maximum_oxidizer_fuel_ratio: float
    reference_oxidizer_fuel_ratio: float
    cea_fuel_name: str
    cea_oxidizer_name: str
    source: str
    verified_chamber_pressure_range_pa: tuple[float, float] = (5.0e5, 2.5e7)


PROPELLANT_PAIRS: dict[str, PropellantPair] = {
    "LOX_LH2": PropellantPair(
        key="LOX_LH2",
        name="LOX / LH2",
        fuel_reactant="H2(L)",
        oxidizer_reactant="O2(L)",
        minimum_oxidizer_fuel_ratio=4.0,
        maximum_oxidizer_fuel_ratio=8.0,
        reference_oxidizer_fuel_ratio=5.55157,
        cea_fuel_name="H2(L)",
        cea_oxidizer_name="O2(L)",
        source="NASA RP-1311 Part II, Example 8",
    ),
    "LOX_CH4": PropellantPair(
        key="LOX_CH4",
        name="LOX / CH4",
        fuel_reactant="CH4(L)",
        oxidizer_reactant="O2(L)",
        minimum_oxidizer_fuel_ratio=2.0,
        maximum_oxidizer_fuel_ratio=5.0,
        reference_oxidizer_fuel_ratio=3.4,
        cea_fuel_name="CH4(L)",
        cea_oxidizer_name="O2(L)",
        source="NASA CEA thermodynamic species convention",
    ),
    "LOX_RP1": PropellantPair(
        key="LOX_RP1",
        name="LOX / RP-1",
        fuel_reactant="RP-1",
        oxidizer_reactant="O2(L)",
        minimum_oxidizer_fuel_ratio=2.0,
        maximum_oxidizer_fuel_ratio=3.5,
        reference_oxidizer_fuel_ratio=2.7,
        cea_fuel_name="RP-1",
        cea_oxidizer_name="O2(L)",
        source="NASA CEA thermo.inp RP-1 record (CH1.95; Mehta et al., AIAA 95-2962)",
    ),
    "NTO_MMH": PropellantPair(
        key="NTO_MMH",
        name="NTO / MMH",
        fuel_reactant="CH6N2(L)",
        oxidizer_reactant="N2O4(L)",
        minimum_oxidizer_fuel_ratio=1.2,
        maximum_oxidizer_fuel_ratio=2.6,
        reference_oxidizer_fuel_ratio=1.65,
        cea_fuel_name="CH6N2(L)",
        cea_oxidizer_name="N2O4(L)",
        source="Screening preset; reactants from NASA CEA thermo.inp N2O4(L) and CH6N2(L)",
    ),
}


@dataclass(frozen=True, slots=True)
class EquilibriumExpansionResult:
    """Isentropic-equilibrium exit using a frozen chamber/throat mass flow.

    The exit state is a shifting-equilibrium expansion [1]_; thrust uses the
    momentum and pressure terms [2]_.

    Attributes
    ----------
    exit_state : EquilibriumResult
        Equilibrium state at the exit pressure with the chamber entropy.
    exit_velocity_m_s : float
        ``sqrt(2 (h_c - h_e))``.
    thrust_n : float
        ``mdot V_e + (p_e - p_a) A_e``.
    thrust_coefficient : float
        ``F / (p_c A_t)``.
    specific_impulse_s : float
        ``F / (mdot g0)``.
    warnings : tuple of WarningMessage
        Statement of the hybrid frozen/equilibrium approximation.

    References
    ----------
    .. [1] NASA RP-1311, "Computer Program for Calculation of Complex Chemical
           Equilibrium Compositions and Applications, Part I: Analysis".
           https://ntrs.nasa.gov/citations/19950013764
    .. [2] SUTTON-BIBLARZ-2017. G. P. Sutton and O. Biblarz, "Rocket Propulsion
           Elements", 9th ed., Wiley, 2017.
    """

    exit_state: EquilibriumResult
    exit_velocity_m_s: float
    thrust_n: float
    thrust_coefficient: float
    specific_impulse_s: float
    warnings: tuple[WarningMessage, ...]


@dataclass(frozen=True, slots=True)
class ChemistryNozzleComparison:
    """Frozen and equilibrium expansion from one equilibrium chamber state.

    Frozen expansion uses ideal-nozzle relations [2]_; equilibrium expansion
    is isentropic with shifting composition [1]_.

    Attributes
    ----------
    chamber : EquilibriumResult
        HP equilibrium chamber state.
    frozen : NozzleResult
        Ideal nozzle with frozen chamber ``gamma`` and gas constant.
    equilibrium : EquilibriumExpansionResult
        Isentropic shifting-equilibrium expansion to the frozen exit pressure.

    References
    ----------
    .. [1] NASA RP-1311, "Computer Program for Calculation of Complex Chemical
           Equilibrium Compositions and Applications, Part I: Analysis".
           https://ntrs.nasa.gov/citations/19950013764
    .. [2] SUTTON-BIBLARZ-2017. G. P. Sutton and O. Biblarz, "Rocket Propulsion
           Elements", 9th ed., Wiley, 2017.
    """

    chamber: EquilibriumResult
    frozen: NozzleResult
    equilibrium: EquilibriumExpansionResult


@dataclass(frozen=True, slots=True)
class OxidizerFuelSweepPoint:
    """One O/F point with chamber and both expansion-mode performance.

    Values come from :func:`compare_nozzle_chemistry` [1]_, [2]_.

    Attributes
    ----------
    oxidizer_fuel_ratio : float
        Oxidizer-to-fuel mass ratio.
    adiabatic_flame_temperature_k, molecular_mass_kg_mol, gamma : float
        Chamber equilibrium temperature, molar mass and frozen ``gamma``.
    characteristic_velocity_m_s : float
        Frozen ``c*``.
    frozen_thrust_coefficient, frozen_vacuum_specific_impulse_s : float
        Frozen-expansion vacuum performance.
    equilibrium_thrust_coefficient, equilibrium_vacuum_specific_impulse_s : float
        Equilibrium-expansion vacuum performance.
    maximum_element_relative_residual : float
        Chamber element-conservation residual.

    References
    ----------
    .. [1] NASA RP-1311, "Computer Program for Calculation of Complex Chemical
           Equilibrium Compositions and Applications, Part I: Analysis".
           https://ntrs.nasa.gov/citations/19950013764
    .. [2] SUTTON-BIBLARZ-2017. G. P. Sutton and O. Biblarz, "Rocket Propulsion
           Elements", 9th ed., Wiley, 2017.
    """

    oxidizer_fuel_ratio: float
    adiabatic_flame_temperature_k: float
    molecular_mass_kg_mol: float
    gamma: float
    characteristic_velocity_m_s: float
    frozen_thrust_coefficient: float
    frozen_vacuum_specific_impulse_s: float
    equilibrium_thrust_coefficient: float
    equilibrium_vacuum_specific_impulse_s: float
    maximum_element_relative_residual: float


@dataclass(frozen=True, slots=True)
class OxidizerFuelSweep:
    """Source-labelled O/F screening result.

    Points are computed as in :func:`generate_oxidizer_fuel_sweep` [1]_, [2]_.

    Attributes
    ----------
    pair : PropellantPair
        Propellant preset.
    chamber_pressure_pa, area_ratio : float
        Fixed chamber pressure and nozzle area ratio.
    points : tuple of OxidizerFuelSweepPoint
        Results in increasing O/F.
    provider, provider_version : str
        Equilibrium provider used.
    warnings : tuple of WarningMessage
        Statement of the reactant reference state and the hybrid expansion.

    References
    ----------
    .. [1] NASA RP-1311, "Computer Program for Calculation of Complex Chemical
           Equilibrium Compositions and Applications, Part I: Analysis".
           https://ntrs.nasa.gov/citations/19950013764
    .. [2] SUTTON-BIBLARZ-2017. G. P. Sutton and O. Biblarz, "Rocket Propulsion
           Elements", 9th ed., Wiley, 2017.
    """

    pair: PropellantPair
    chamber_pressure_pa: float
    area_ratio: float
    points: tuple[OxidizerFuelSweepPoint, ...]
    provider: str
    provider_version: str
    warnings: tuple[WarningMessage, ...]


def get_propellant_pair(key: str) -> PropellantPair:
    """Return a propellant preset by a case-insensitive key.

    Parameters
    ----------
    key : str
        Preset key; case, ``/`` and ``-`` separators are normalized, so
        ``"lox/lh2"`` selects ``"LOX_LH2"``.

    Returns
    -------
    PropellantPair
        Preset record.

    Raises
    ------
    InputError
        If the key is unknown (code ``unknown_propellant_pair``).
    """

    normalized = key.strip().upper().replace("/", "_").replace("-", "_")
    try:
        return PROPELLANT_PAIRS[normalized]
    except KeyError as error:
        raise InputError(
            f"Unknown propellant pair: {key!r}.",
            code="unknown_propellant_pair",
            field="propellant_pair",
            details={"supported": sorted(PROPELLANT_PAIRS)},
        ) from error


def reactants_for_mass_ratio(pair: PropellantPair, ratio: float) -> dict[str, float]:
    """Build one-mole-fuel reactants from an oxidizer/fuel mass ratio.

    Parameters
    ----------
    pair : PropellantPair
        Propellant preset.
    ratio : float
        Oxidizer-to-fuel mass ratio, finite and positive.

    Returns
    -------
    dict of str to float
        ``{fuel: 1.0, oxidizer: ratio M_fuel / M_oxidizer}`` in mol.

    Raises
    ------
    DomainError
        If ``ratio`` is not finite and positive.
    """

    if not isfinite(ratio) or ratio <= 0.0:
        raise DomainError("O/F ratio must be positive.", field="oxidizer_fuel_ratio")
    fuel = get_reactant(pair.fuel_reactant)
    oxidizer = get_reactant(pair.oxidizer_reactant)
    oxidizer_moles = ratio * fuel.molecular_mass_kg_mol / oxidizer.molecular_mass_kg_mol
    return {fuel.key: 1.0, oxidizer.key: oxidizer_moles}


def _isentropic_equilibrium_exit(
    *,
    chamber: EquilibriumResult,
    problem: EquilibriumProblem,
    provider: EquilibriumProvider,
    exit_pressure_pa: float,
    exit_area_m2: float,
    ambient_pressure_pa: float,
    mass_flow_kg_s: float,
    chamber_pressure_pa: float,
    throat_area_m2: float,
) -> EquilibriumExpansionResult:
    """Expand isentropically in shifting equilibrium to ``exit_pressure_pa``.

    Raises
    ------
    DomainError
        If the chamber entropy is not bracketed or the enthalpy drop is not
        positive.

    Notes
    -----
    The exit temperature is bisected in ``[200 K, T_c]`` until the
    equilibrium entropy at the exit pressure equals the chamber entropy [1]_.
    The exit velocity follows from the energy equation,
    ``V_e = sqrt(2 (h_c - h_e))``, and thrust from
    ``F = mdot V_e + (p_e - p_a) A_e`` [2]_, with the frozen-chamber throat
    mass flow.

    References
    ----------
    .. [1] NASA RP-1311, "Computer Program for Calculation of Complex Chemical
           Equilibrium Compositions and Applications, Part I: Analysis".
           https://ntrs.nasa.gov/citations/19950013764
    .. [2] SUTTON-BIBLARZ-2017. G. P. Sutton and O. Biblarz, "Rocket Propulsion
           Elements", 9th ed., Wiley, 2017.
    """
    target_entropy = chamber.entropy_j_k
    low = 200.0
    high = chamber.temperature_k
    latest: EquilibriumResult | None = None

    def residual(temperature: float) -> float:
        nonlocal latest
        latest = provider.solve(
            EquilibriumProblem.tp(
                {item.species: item.moles for item in problem.reactants},
                temperature_k=temperature,
                pressure_pa=exit_pressure_pa,
                candidate_species=problem.candidate_species,
            )
        )
        return (latest.entropy_j_k - target_entropy) / max(abs(target_entropy), 1.0)

    low_residual = residual(low)
    high_residual = residual(high)
    if low_residual * high_residual > 0.0:
        raise DomainError(
            "Equilibrium exit entropy is not bracketed by the local model range.",
            code="equilibrium_expansion_not_bracketed",
            details={"low_residual": low_residual, "high_residual": high_residual},
        )
    for _ in range(80):
        midpoint = 0.5 * (low + high)
        value = residual(midpoint)
        if abs(value) <= 1.0e-9 or (high - low) <= 1.0e-7:
            break
        if low_residual * value <= 0.0:
            high = midpoint
        else:
            low = midpoint
            low_residual = value
    assert latest is not None
    enthalpy_drop_j_kg = (chamber.enthalpy_j - latest.enthalpy_j) / chamber.total_mass_kg
    if enthalpy_drop_j_kg <= 0.0:
        raise DomainError(
            "Equilibrium expansion produced a non-positive enthalpy drop.",
            code="invalid_equilibrium_expansion",
        )
    exit_velocity = sqrt(2.0 * enthalpy_drop_j_kg)
    thrust = mass_flow_kg_s * exit_velocity + (
        (exit_pressure_pa - ambient_pressure_pa) * exit_area_m2
    )
    return EquilibriumExpansionResult(
        exit_state=latest,
        exit_velocity_m_s=exit_velocity,
        thrust_n=thrust,
        thrust_coefficient=thrust / (chamber_pressure_pa * throat_area_m2),
        specific_impulse_s=thrust / (mass_flow_kg_s * STANDARD_GRAVITY_M_S2),
        warnings=(
            WarningMessage(
                code="hybrid_equilibrium_nozzle",
                message=(
                    "Equilibrium exit chemistry is solved isentropically, while throat mass "
                    "flow uses frozen chamber properties. Use NASA CEA for final validation."
                ),
            ),
        ),
    )


def compare_nozzle_chemistry(
    *,
    reactants: Mapping[str, float],
    chamber_pressure_pa: float,
    area_ratio: float,
    ambient_pressure_pa: float = 0.0,
    throat_area_m2: float = 1.0,
    reactant_temperature_k: float = 298.15,
    provider: EquilibriumProvider | None = None,
) -> ChemistryNozzleComparison:
    """Compare frozen and experimental equilibrium nozzle expansion.

    Parameters
    ----------
    reactants : mapping of str to float
        Reactant amounts in mol per species key.
    chamber_pressure_pa : float
        Chamber pressure in Pa.
    area_ratio : float
        Exit-to-throat area ratio.
    ambient_pressure_pa : float, default 0.0
        Ambient pressure in Pa.
    throat_area_m2 : float, default 1.0
        Throat area in m^2.
    reactant_temperature_k : float, default 298.15
        Reactant temperature in K.
    provider : EquilibriumProvider, optional
        Equilibrium provider; :class:`LocalEquilibriumProvider` by default.

    Returns
    -------
    ChemistryNozzleComparison
        Chamber state with frozen and equilibrium expansions.

    Raises
    ------
    InputError, DomainError, ConvergenceError
        Propagated from the equilibrium provider and the nozzle relations.

    Notes
    -----
    The chamber is an HP equilibrium state [1]_. The frozen expansion uses the
    chamber ``gamma`` and gas constant in the ideal-nozzle relations [2]_; the
    equilibrium expansion is isentropic to the same exit pressure with the
    frozen throat mass flow. Kinetic, two-phase, boundary-layer and divergence
    losses are not included.

    References
    ----------
    .. [1] NASA RP-1311, "Computer Program for Calculation of Complex Chemical
           Equilibrium Compositions and Applications, Part I: Analysis".
           https://ntrs.nasa.gov/citations/19950013764
    .. [2] SUTTON-BIBLARZ-2017. G. P. Sutton and O. Biblarz, "Rocket Propulsion
           Elements", 9th ed., Wiley, 2017.
    """

    active_provider = provider or LocalEquilibriumProvider()
    chamber_problem = EquilibriumProblem.hp(
        reactants,
        pressure_pa=chamber_pressure_pa,
        reactant_temperature_k=reactant_temperature_k,
    )
    chamber = active_provider.solve(chamber_problem)
    frozen = calculate_ideal_nozzle(
        chamber_pressure_pa=chamber_pressure_pa,
        chamber_temperature_k=chamber.temperature_k,
        throat_area_m2=throat_area_m2,
        area_ratio=area_ratio,
        ambient_pressure_pa=ambient_pressure_pa,
        gamma=chamber.gamma,
        gas_constant_j_kg_k=chamber.gas_constant_j_kg_k,
    )
    equilibrium = _isentropic_equilibrium_exit(
        chamber=chamber,
        problem=chamber_problem,
        provider=active_provider,
        exit_pressure_pa=frozen.exit_pressure_pa,
        exit_area_m2=frozen.exit_area_m2,
        ambient_pressure_pa=ambient_pressure_pa,
        mass_flow_kg_s=frozen.mass_flow_kg_s,
        chamber_pressure_pa=chamber_pressure_pa,
        throat_area_m2=throat_area_m2,
    )
    return ChemistryNozzleComparison(
        chamber=chamber,
        frozen=frozen,
        equilibrium=equilibrium,
    )


def generate_oxidizer_fuel_sweep(
    *,
    propellant_pair: str,
    chamber_pressure_pa: float,
    area_ratio: float,
    minimum_ratio: float | None = None,
    maximum_ratio: float | None = None,
    point_count: int = 9,
    provider: EquilibriumProvider | None = None,
) -> OxidizerFuelSweep:
    """Generate chamber, c-star, coefficient, and Isp curves over O/F.

    Parameters
    ----------
    propellant_pair : str
        Preset key accepted by :func:`get_propellant_pair`.
    chamber_pressure_pa : float
        Chamber pressure in Pa.
    area_ratio : float
        Exit-to-throat area ratio.
    minimum_ratio, maximum_ratio : float, optional
        Sweep limits; the preset range is used when omitted.
    point_count : int, default 9
        Number of evenly spaced O/F points, between 3 and 41.
    provider : EquilibriumProvider, optional
        Equilibrium provider; :class:`LocalEquilibriumProvider` by default.

    Returns
    -------
    OxidizerFuelSweep
        Per-point chamber and vacuum performance.

    Raises
    ------
    InputError
        If the preset is unknown or needs the CEA provider.
    DomainError
        If the limits are not positive and increasing or the point count is
        out of range.

    Notes
    -----
    Each point calls :func:`compare_nozzle_chemistry` with vacuum ambient
    pressure [1]_, [2]_. Reactants enter at the assigned storage-state
    enthalpy of NASA CEA's database (for example O2(L) at 90.17 K). The sweep
    keeps the fast hybrid expansion of :func:`compare_nozzle_chemistry`; use
    :func:`~rocket_propulsion.thermochemistry.rocket_performance.solve_rocket_performance`
    for the CEA-equivalent shifting-equilibrium throat and exit.

    References
    ----------
    .. [1] NASA RP-1311, "Computer Program for Calculation of Complex Chemical
           Equilibrium Compositions and Applications, Part I: Analysis".
           https://ntrs.nasa.gov/citations/19950013764
    .. [2] SUTTON-BIBLARZ-2017. G. P. Sutton and O. Biblarz, "Rocket Propulsion
           Elements", 9th ed., Wiley, 2017.
    """

    pair = get_propellant_pair(propellant_pair)
    low = pair.minimum_oxidizer_fuel_ratio if minimum_ratio is None else minimum_ratio
    high = pair.maximum_oxidizer_fuel_ratio if maximum_ratio is None else maximum_ratio
    if not isfinite(low) or not isfinite(high) or low <= 0.0 or high <= low:
        raise DomainError("O/F sweep limits must be positive and increasing.")
    if not 3 <= point_count <= 41:
        raise DomainError("O/F sweep point count must be between 3 and 41.")
    active_provider = provider or LocalEquilibriumProvider()
    step = (high - low) / (point_count - 1)
    points: list[OxidizerFuelSweepPoint] = []
    for index in range(point_count):
        ratio = low + index * step
        comparison = compare_nozzle_chemistry(
            reactants=reactants_for_mass_ratio(pair, ratio),
            chamber_pressure_pa=chamber_pressure_pa,
            area_ratio=area_ratio,
            provider=active_provider,
        )
        points.append(
            OxidizerFuelSweepPoint(
                oxidizer_fuel_ratio=ratio,
                adiabatic_flame_temperature_k=comparison.chamber.temperature_k,
                molecular_mass_kg_mol=comparison.chamber.molecular_mass_kg_mol,
                gamma=comparison.chamber.gamma,
                characteristic_velocity_m_s=comparison.frozen.characteristic_velocity_m_s,
                frozen_thrust_coefficient=comparison.frozen.thrust_coefficient,
                frozen_vacuum_specific_impulse_s=comparison.frozen.specific_impulse_s,
                equilibrium_thrust_coefficient=comparison.equilibrium.thrust_coefficient,
                equilibrium_vacuum_specific_impulse_s=(
                    comparison.equilibrium.specific_impulse_s
                ),
                maximum_element_relative_residual=(
                    comparison.chamber.residuals.maximum_element_relative
                ),
            )
        )
    return OxidizerFuelSweep(
        pair=pair,
        chamber_pressure_pa=chamber_pressure_pa,
        area_ratio=area_ratio,
        points=tuple(points),
        provider=active_provider.name,
        provider_version=active_provider.version,
        warnings=(
            WarningMessage(
                code="hybrid_sweep_expansion",
                message=(
                    "Reactants use NASA CEA storage-state enthalpies; the sweep's throat "
                    "flow is frozen at the chamber. Use the rocket-performance solve for "
                    "CEA-equivalent equilibrium values."
                ),
            ),
        ),
    )

