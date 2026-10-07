"""Optional adapters that keep external frameworks out of the physics core."""

from .contracts import SideraCapabilities
from .lunaris import (
    LunarisFiniteBurnDynamics,
    LunarisImpulsiveEquivalent,
    lunaris_impulsive_equivalents,
    lunaris_segmentation_plan,
    lunaris_segmentation_times,
    to_lunaris_maneuver_plan,
    velocity_direction,
)
from .sidera import (
    export_sidera_artifact,
    inspect_sidera_capabilities,
    to_native_finite_burn,
    to_native_maneuver_plan,
)
from .sidera_ensemble import (
    SIDERA_ENSEMBLE_HANDOFF_SCHEMA,
    SIDERA_ENSEMBLE_REQUIRED_CAPABILITIES,
    SideraEnsembleBurnJob,
    SideraEnsembleHandoff,
    export_sidera_ensemble_handoff,
    missing_sidera_ensemble_capabilities,
    require_sidera_ensemble_capabilities,
)

__all__ = [
    "SIDERA_ENSEMBLE_HANDOFF_SCHEMA",
    "SIDERA_ENSEMBLE_REQUIRED_CAPABILITIES",
    "LunarisFiniteBurnDynamics",
    "LunarisImpulsiveEquivalent",
    "SideraCapabilities",
    "SideraEnsembleBurnJob",
    "SideraEnsembleHandoff",
    "export_sidera_artifact",
    "export_sidera_ensemble_handoff",
    "inspect_sidera_capabilities",
    "lunaris_impulsive_equivalents",
    "lunaris_segmentation_plan",
    "lunaris_segmentation_times",
    "missing_sidera_ensemble_capabilities",
    "require_sidera_ensemble_capabilities",
    "to_lunaris_maneuver_plan",
    "to_native_finite_burn",
    "to_native_maneuver_plan",
    "velocity_direction",
]

