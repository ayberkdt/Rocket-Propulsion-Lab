"""Public capability records shared by optional external adapters."""

from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True, slots=True)
class SideraCapabilities:
    """Observed public Sidera finite-burn capabilities."""

    status: str
    version: str | None
    finite_burn_available: bool
    maneuver_plan_available: bool
    constructor_compatible: bool
    supported_frames: tuple[str, ...]
    constant_burn_exact: bool
    tabulated_burn_available: bool
    independent_mass_flow_available: bool
    reason: str | None = None
    weighted_maneuver_ensemble_available: bool = False
    exact_weighted_tabulated_ensemble: bool = False
    propulsion_provider_available: bool = False
    additional_states_available: bool = False
    event_detectors_available: bool = False
    analytic_partials_available: bool = False
    body_frame_force_moment_available: bool = False
