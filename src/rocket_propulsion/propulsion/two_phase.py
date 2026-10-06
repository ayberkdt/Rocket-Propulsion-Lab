"""Preliminary performance bands for condensed-phase exhaust."""

from __future__ import annotations

from dataclasses import dataclass
from math import isfinite

from rocket_propulsion.core.errors import DomainError
from rocket_propulsion.core.metadata import WarningMessage


@dataclass(frozen=True, slots=True)
class TwoPhasePerformanceBand:
    """Bound ideal performance using a user-supplied particle coupling band."""

    condensed_mass_fraction: float
    particle_coupling_efficiency_range: tuple[float, float]
    performance_efficiency_range: tuple[float, float]
    thrust_range_n: tuple[float, float]
    specific_impulse_range_s: tuple[float, float]
    warnings: tuple[WarningMessage, ...]
    model: str


def calculate_two_phase_band(
    *,
    ideal_thrust_n: float,
    ideal_specific_impulse_s: float,
    condensed_mass_fraction: float,
    minimum_particle_coupling_efficiency: float,
    maximum_particle_coupling_efficiency: float,
) -> TwoPhasePerformanceBand:
    """Return a transparent bracket, not a hidden monodisperse particle model."""

    for field, value in {
        "ideal_thrust_n": ideal_thrust_n,
        "ideal_specific_impulse_s": ideal_specific_impulse_s,
    }.items():
        if not isfinite(value) or value <= 0.0:
            raise DomainError(f"{field} must be finite and positive.", field=field)
    if not isfinite(condensed_mass_fraction) or not 0.0 <= condensed_mass_fraction <= 1.0:
        raise DomainError("Condensed mass fraction must be in [0, 1].")
    minimum = minimum_particle_coupling_efficiency
    maximum = maximum_particle_coupling_efficiency
    if not all(isfinite(value) and 0.0 <= value <= 1.0 for value in (minimum, maximum)):
        raise DomainError("Particle coupling efficiencies must be in [0, 1].")
    if minimum > maximum:
        raise DomainError("Minimum particle coupling efficiency cannot exceed maximum.")
    lower_efficiency = 1.0 - condensed_mass_fraction * (1.0 - minimum)
    upper_efficiency = 1.0 - condensed_mass_fraction * (1.0 - maximum)
    return TwoPhasePerformanceBand(
        condensed_mass_fraction=condensed_mass_fraction,
        particle_coupling_efficiency_range=(minimum, maximum),
        performance_efficiency_range=(lower_efficiency, upper_efficiency),
        thrust_range_n=(ideal_thrust_n * lower_efficiency, ideal_thrust_n * upper_efficiency),
        specific_impulse_range_s=(
            ideal_specific_impulse_s * lower_efficiency,
            ideal_specific_impulse_s * upper_efficiency,
        ),
        warnings=(
            WarningMessage(
                code="two_phase_band_only",
                message=(
                    "This is a coupling-efficiency band; particle size, residence time, "
                    "phase change, and slip are not resolved."
                ),
                severity="warning",
            ),
        ),
        model="condensed-mass-fraction particle-coupling bracket",
    )
