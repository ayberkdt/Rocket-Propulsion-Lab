"""Traceable preliminary rocket performance loss budgets."""

from __future__ import annotations

from dataclasses import dataclass
from math import isfinite, sqrt

from rocket_propulsion.core.errors import DomainError


@dataclass(frozen=True, slots=True)
class LossTerm:
    """One sequential, switchable efficiency contribution."""

    key: str
    label: str
    enabled: bool
    efficiency: float
    thrust_loss_n: float
    specific_impulse_loss_s: float
    source: str


@dataclass(frozen=True, slots=True)
class PerformanceLossBudget:
    """Ideal-to-delivered performance bridge without double-counted losses."""

    ideal_thrust_n: float
    delivered_thrust_n: float
    ideal_specific_impulse_s: float
    delivered_specific_impulse_s: float
    total_efficiency: float
    total_thrust_loss_n: float
    total_specific_impulse_loss_s: float
    terms: tuple[LossTerm, ...]
    thrust_balance_residual_n: float
    specific_impulse_balance_residual_s: float
    one_sigma_relative: float
    delivered_thrust_band_n: tuple[float, float]
    delivered_specific_impulse_band_s: tuple[float, float]
    assumptions: tuple[str, ...]


def calculate_loss_budget(
    *,
    ideal_thrust_n: float,
    ideal_specific_impulse_s: float,
    discharge_coefficient: float = 0.98,
    divergence_efficiency: float = 0.98,
    boundary_layer_efficiency: float = 0.985,
    combustion_efficiency: float = 0.97,
    enable_discharge_loss: bool = True,
    enable_divergence_loss: bool = True,
    enable_boundary_layer_loss: bool = True,
    enable_combustion_loss: bool = True,
    relative_uncertainty: float = 0.02,
) -> PerformanceLossBudget:
    """Apply named efficiencies sequentially and expose each performance delta."""

    if not isfinite(ideal_thrust_n) or ideal_thrust_n <= 0.0:
        raise DomainError("Ideal thrust must be finite and positive.", field="ideal_thrust_n")
    if not isfinite(ideal_specific_impulse_s) or ideal_specific_impulse_s <= 0.0:
        raise DomainError(
            "Ideal specific impulse must be finite and positive.",
            field="ideal_specific_impulse_s",
        )
    settings = (
        (
            "combustion",
            "Yanma verimi",
            enable_combustion_loss,
            combustion_efficiency,
            "characteristic-velocity efficiency",
        ),
        (
            "discharge",
            "Deşarj katsayısı",
            enable_discharge_loss,
            discharge_coefficient,
            "effective throat mass-flow coefficient",
        ),
        (
            "boundary_layer",
            "Sınır tabakası",
            enable_boundary_layer_loss,
            boundary_layer_efficiency,
            "lumped viscous momentum efficiency",
        ),
        (
            "divergence",
            "Iraksama",
            enable_divergence_loss,
            divergence_efficiency,
            "axisymmetric momentum projection estimate",
        ),
    )
    for key, _, _, efficiency, _ in settings:
        if not isfinite(efficiency) or not 0.0 < efficiency <= 1.0:
            raise DomainError(
                f"{key} efficiency must be finite and in (0, 1].", field=key
            )
    if not isfinite(relative_uncertainty) or not 0.0 <= relative_uncertainty <= 0.25:
        raise DomainError("Relative uncertainty must be in [0, 0.25].")

    current_thrust = ideal_thrust_n
    current_isp = ideal_specific_impulse_s
    terms: list[LossTerm] = []
    enabled_count = 0
    for key, label, enabled, efficiency, source in settings:
        applied = efficiency if enabled else 1.0
        next_thrust = current_thrust * applied
        next_isp = current_isp * applied
        terms.append(
            LossTerm(
                key=key,
                label=label,
                enabled=enabled,
                efficiency=efficiency,
                thrust_loss_n=current_thrust - next_thrust,
                specific_impulse_loss_s=current_isp - next_isp,
                source=source,
            )
        )
        current_thrust, current_isp = next_thrust, next_isp
        enabled_count += int(enabled)
    total_efficiency = current_thrust / ideal_thrust_n
    thrust_loss = sum(term.thrust_loss_n for term in terms)
    isp_loss = sum(term.specific_impulse_loss_s for term in terms)
    combined_uncertainty = relative_uncertainty * sqrt(max(enabled_count, 1))
    return PerformanceLossBudget(
        ideal_thrust_n=ideal_thrust_n,
        delivered_thrust_n=current_thrust,
        ideal_specific_impulse_s=ideal_specific_impulse_s,
        delivered_specific_impulse_s=current_isp,
        total_efficiency=total_efficiency,
        total_thrust_loss_n=thrust_loss,
        total_specific_impulse_loss_s=isp_loss,
        terms=tuple(terms),
        thrust_balance_residual_n=abs(ideal_thrust_n - current_thrust - thrust_loss),
        specific_impulse_balance_residual_s=abs(
            ideal_specific_impulse_s - current_isp - isp_loss
        ),
        one_sigma_relative=combined_uncertainty,
        delivered_thrust_band_n=(
            current_thrust * (1.0 - combined_uncertainty),
            current_thrust * (1.0 + combined_uncertainty),
        ),
        delivered_specific_impulse_band_s=(
            current_isp * (1.0 - combined_uncertainty),
            current_isp * (1.0 + combined_uncertainty),
        ),
        assumptions=(
            "sequential independent lumped efficiencies",
            "uncertainty terms combined by root-sum-square",
            "preliminary performance accounting, not viscous CFD",
        ),
    )
