"""Multi-stage delta-v accounting and constrained mass optimization."""

from __future__ import annotations

from dataclasses import dataclass
from math import exp, isfinite, log

from rocket_propulsion.core.errors import ConvergenceError, DomainError

from .performance import STANDARD_GRAVITY_M_S2


@dataclass(frozen=True, slots=True)
class StageDefinition:
    """A stage propulsion/structure model ordered from lower to upper stage."""

    name: str
    specific_impulse_s: float
    structural_fraction: float


@dataclass(frozen=True, slots=True)
class StageSolution:
    name: str
    delta_v_m_s: float
    mass_ratio: float
    initial_mass_kg: float
    propellant_mass_kg: float
    dry_mass_kg: float
    upper_stack_mass_kg: float
    active_constraint: str | None


@dataclass(frozen=True, slots=True)
class StagingOptimizationResult:
    target_delta_v_m_s: float
    achieved_delta_v_m_s: float
    payload_mass_kg: float
    initial_mass_kg: float
    payload_fraction: float
    stages: tuple[StageSolution, ...]
    discretization_step_m_s: float
    convergence_delta_initial_mass_relative: float
    active_constraints: tuple[str, ...]
    feasible: bool


def _validate(stages: tuple[StageDefinition, ...], target: float, payload: float) -> None:
    if not 1 <= len(stages) <= 5:
        raise DomainError("Stage count must be between 1 and 5.")
    if not isfinite(target) or target <= 0.0:
        raise DomainError("Target delta-v must be finite and positive.")
    if not isfinite(payload) or payload <= 0.0:
        raise DomainError("Payload mass must be finite and positive.")
    for stage in stages:
        if not stage.name.strip():
            raise DomainError("Stage name cannot be empty.")
        if not isfinite(stage.specific_impulse_s) or stage.specific_impulse_s <= 0.0:
            raise DomainError("Stage specific impulse must be positive.")
        if not isfinite(stage.structural_fraction) or not 0.0 < stage.structural_fraction < 1.0:
            raise DomainError("Stage structural fraction must be in (0, 1).")


def _mass_factor(stage: StageDefinition, delta_v_m_s: float) -> float:
    mass_ratio = exp(delta_v_m_s / (STANDARD_GRAVITY_M_S2 * stage.specific_impulse_s))
    denominator = 1.0 - mass_ratio * stage.structural_fraction
    if denominator <= 0.0:
        return float("inf")
    return mass_ratio * (1.0 - stage.structural_fraction) / denominator


def _optimize_units(
    stages: tuple[StageDefinition, ...], target: float, resolution: int
) -> tuple[list[int], float]:
    step = target / resolution
    infinity = float("inf")
    costs = [infinity] * (resolution + 1)
    costs[0] = 0.0
    predecessors: list[list[int]] = []
    for stage in stages:
        next_costs = [infinity] * (resolution + 1)
        selected = [-1] * (resolution + 1)
        factors = [_mass_factor(stage, units * step) for units in range(resolution + 1)]
        for total_units in range(resolution + 1):
            best = infinity
            best_units = -1
            for stage_units in range(total_units + 1):
                prior = costs[total_units - stage_units]
                factor = factors[stage_units]
                if prior == infinity or factor == infinity:
                    continue
                candidate = prior + log(factor)
                if candidate < best:
                    best = candidate
                    best_units = stage_units
            next_costs[total_units] = best
            selected[total_units] = best_units
        costs = next_costs
        predecessors.append(selected)
    if costs[resolution] == infinity:
        raise DomainError(
            "Target delta-v is infeasible for the supplied structural fractions and Isp values.",
            code="infeasible_staging_target",
        )
    allocations = [0] * len(stages)
    remaining = resolution
    for index in range(len(stages) - 1, -1, -1):
        stage_units = predecessors[index][remaining]
        if stage_units < 0:
            raise ConvergenceError("Staging optimizer predecessor chain is incomplete.")
        allocations[index] = stage_units
        remaining -= stage_units
    return allocations, step


def optimize_staging(
    *,
    stages: tuple[StageDefinition, ...],
    target_delta_v_m_s: float,
    payload_mass_kg: float,
    resolution: int = 240,
) -> StagingOptimizationResult:
    """Minimize initial mass using deterministic discrete delta-v allocation."""

    _validate(stages, target_delta_v_m_s, payload_mass_kg)
    if not 40 <= resolution <= 800:
        raise DomainError("Staging optimization resolution must be in [40, 800].")
    allocations, step = _optimize_units(stages, target_delta_v_m_s, resolution)
    refined_resolution = min(800, resolution * 2)
    refined_allocations, refined_step = _optimize_units(
        stages, target_delta_v_m_s, refined_resolution
    )

    def initial_mass(units: list[int], unit_step: float) -> float:
        mass = payload_mass_kg
        for stage, count in zip(reversed(stages), reversed(units)):
            mass *= _mass_factor(stage, count * unit_step)
        return mass

    coarse_mass = initial_mass(allocations, step)
    refined_mass = initial_mass(refined_allocations, refined_step)
    upper_mass = payload_mass_kg
    solutions_reversed: list[StageSolution] = []
    constraints: list[str] = []
    for stage, count in zip(reversed(stages), reversed(allocations)):
        delta_v = count * step
        mass_ratio = exp(delta_v / (STANDARD_GRAVITY_M_S2 * stage.specific_impulse_s))
        denominator = 1.0 - mass_ratio * stage.structural_fraction
        stage_mass = upper_mass * (mass_ratio - 1.0) / denominator
        dry_mass = stage.structural_fraction * stage_mass
        propellant_mass = stage_mass - dry_mass
        limit = -STANDARD_GRAVITY_M_S2 * stage.specific_impulse_s * log(
            stage.structural_fraction
        )
        active = None
        if delta_v <= step:
            active = "minimum delta-v allocation"
        elif limit - delta_v <= 2.0 * step:
            active = "structural-fraction mass-ratio limit"
        if active:
            constraints.append(f"{stage.name}: {active}")
        initial = upper_mass + stage_mass
        solutions_reversed.append(
            StageSolution(
                name=stage.name,
                delta_v_m_s=delta_v,
                mass_ratio=mass_ratio,
                initial_mass_kg=initial,
                propellant_mass_kg=propellant_mass,
                dry_mass_kg=dry_mass,
                upper_stack_mass_kg=upper_mass,
                active_constraint=active,
            )
        )
        upper_mass = initial
    solutions = tuple(reversed(solutions_reversed))
    return StagingOptimizationResult(
        target_delta_v_m_s=target_delta_v_m_s,
        achieved_delta_v_m_s=sum(stage.delta_v_m_s for stage in solutions),
        payload_mass_kg=payload_mass_kg,
        initial_mass_kg=coarse_mass,
        payload_fraction=payload_mass_kg / coarse_mass,
        stages=solutions,
        discretization_step_m_s=step,
        convergence_delta_initial_mass_relative=abs(refined_mass - coarse_mass)
        / refined_mass,
        active_constraints=tuple(constraints),
        feasible=True,
    )
