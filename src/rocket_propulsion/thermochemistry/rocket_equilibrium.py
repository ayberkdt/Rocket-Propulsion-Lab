"""Experimental equilibrium-chamber and nozzle-expansion comparisons."""

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
from .species import get_species


@dataclass(frozen=True, slots=True)
class PropellantPair:
    """Source-labelled preset for local screening and optional CEA runs."""

    key: str
    name: str
    fuel_species: str | None
    oxidizer_species: str
    minimum_oxidizer_fuel_ratio: float
    maximum_oxidizer_fuel_ratio: float
    reference_oxidizer_fuel_ratio: float
    cea_fuel_name: str
    cea_oxidizer_name: str
    source: str


PROPELLANT_PAIRS: dict[str, PropellantPair] = {
    "LOX_LH2": PropellantPair(
        key="LOX_LH2",
        name="LOX / LH2",
        fuel_species="H2",
        oxidizer_species="O2",
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
        fuel_species="CH4",
        oxidizer_species="O2",
        minimum_oxidizer_fuel_ratio=2.0,
        maximum_oxidizer_fuel_ratio=5.0,
        reference_oxidizer_fuel_ratio=3.4,
        cea_fuel_name="CH4(L)",
        cea_oxidizer_name="O2(L)",
        source="NASA CEA thermodynamic species convention",
    ),
    "LOX_RP1": PropellantPair(
        key="LOX_RP1",
        name="LOX / RP-1 surrogate",
        fuel_species=None,
        oxidizer_species="O2",
        minimum_oxidizer_fuel_ratio=2.0,
        maximum_oxidizer_fuel_ratio=3.5,
        reference_oxidizer_fuel_ratio=2.7,
        cea_fuel_name="C12H26(L)",
        cea_oxidizer_name="O2(L)",
        source="C12H26 liquid surrogate; CEA validation required",
    ),
}


@dataclass(frozen=True, slots=True)
class EquilibriumExpansionResult:
    """Isentropic-equilibrium exit using a frozen chamber/throat mass flow."""

    exit_state: EquilibriumResult
    exit_velocity_m_s: float
    thrust_n: float
    thrust_coefficient: float
    specific_impulse_s: float
    warnings: tuple[WarningMessage, ...]


@dataclass(frozen=True, slots=True)
class ChemistryNozzleComparison:
    """Frozen and equilibrium expansion from one equilibrium chamber state."""

    chamber: EquilibriumResult
    frozen: NozzleResult
    equilibrium: EquilibriumExpansionResult


@dataclass(frozen=True, slots=True)
class OxidizerFuelSweepPoint:
    """One O/F point with chamber and both expansion-mode performance."""

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
    """Source-labelled O/F screening result."""

    pair: PropellantPair
    chamber_pressure_pa: float
    area_ratio: float
    points: tuple[OxidizerFuelSweepPoint, ...]
    provider: str
    provider_version: str
    warnings: tuple[WarningMessage, ...]


def get_propellant_pair(key: str) -> PropellantPair:
    """Return a propellant preset by a case-insensitive key."""

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
    """Build one-mole-fuel reactants from an oxidizer/fuel mass ratio."""

    if pair.fuel_species is None:
        raise InputError(
            f"{pair.name} requires the NASA CEA provider; no local RP-1 species is bundled.",
            code="provider_required",
            field="propellant_pair",
        )
    if not isfinite(ratio) or ratio <= 0.0:
        raise DomainError("O/F ratio must be positive.", field="oxidizer_fuel_ratio")
    fuel = get_species(pair.fuel_species)
    oxidizer = get_species(pair.oxidizer_species)
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
    """Compare frozen and experimental equilibrium nozzle expansion."""

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
    """Generate chamber, c-star, coefficient, and Isp curves over O/F."""

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
                code="gaseous_reactant_reference",
                message=(
                    "Local sweeps use bundled gas species at 298.15 K. Select the CEA "
                    "adapter with real liquid reactants for release-grade comparisons."
                ),
            ),
        ),
    )

