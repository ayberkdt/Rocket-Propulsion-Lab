"""Public capability records shared by optional external adapters."""

from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True, slots=True)
class SideraCapabilities:
    """Capabilities observed on the installed public Sidera interface.

    Every flag records what inspection of the installed package found; a
    False flag means the capability was not observed, not that it is known
    to be absent upstream.

    Attributes
    ----------
    status : str
        ``"supported"``, ``"incompatible"`` or ``"unavailable"``.
    version : str or None
        Reported Sidera version, ``"unknown"`` if it exposes none, ``None``
        when Sidera is not importable.
    finite_burn_available, maneuver_plan_available : bool
        Public ``FiniteBurn`` and ``ManeuverPlan`` types were found.
    constructor_compatible : bool
        ``FiniteBurn`` accepts every field the adapter writes.
    supported_frames : tuple of str
        Frames the adapter may emit; empty when incompatible.
    constant_burn_exact : bool
        A constant burn can be handed over without approximation.
    tabulated_burn_available : bool
        A public ``TabulatedFiniteBurn`` type was found.
    independent_mass_flow_available : bool
        That type accepts a mass flow independent of thrust and Isp.
    reason : str or None
        Explanation when a capability is missing or inspection failed.
    weighted_maneuver_ensemble_available : bool
        A public ``WeightedManeuverEnsemble`` type was found.
    exact_weighted_tabulated_ensemble : bool
        Tabulated burns with independent mass flow and weighted ensembles are
        all available.
    propulsion_provider_available : bool
        A public ``PropulsionProvider`` type was found.
    additional_states_available, event_detectors_available,
    analytic_partials_available, body_frame_force_moment_available : bool
        The provider constructor exposes the corresponding parameter.
    """

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
