"""Source-cited preliminary propulsion design tools for provisional UZM 508E projects.

References
----------
ITU-UZM508: course identity only. NASA-SE-APP: engineering study traceability.
See docs/SASS_PROJECT_PREPARATION.md for confirmed facts versus inferred scope.
"""

from .catalog import (
    THRUSTERS,
    CatalogField,
    ThrusterRecord,
    list_thrusters,
    nominal_operating_point,
    thruster,
)
from .evidence import AnalysisResult, Citation, CitedScalar, assumption
from .maneuvers import (
    burn_duration_screen,
    circular_spiral_delta_v,
    hohmann_transfer,
    perigee_lowering_delta_v,
)
from .references import REFERENCES, Reference, reference
from .project import Candidate, CandidateStudy, ProjectStudy, analyze_project, project_from_dict, project_from_json
from .report import build_project_html
from .sensitivity import DesignCase, design_sensitivity
from .sizing import (
    DeltaVItem,
    Maneuver,
    MissionSizing,
    bipropellant_split,
    delta_v_budget,
    power_resources,
    size_mission,
)
from .storage import (
    ideal_cold_gas_storage,
    ideal_regulated_pressurant,
    liquid_tank_volume,
)
from .trade import (
    AlternativeScore,
    Criterion,
    Requirement,
    RequirementCheck,
    check_requirements,
    weighted_trade,
)

__all__ = [
    "REFERENCES", "THRUSTERS", "AlternativeScore", "AnalysisResult", "CatalogField", "Citation",
    "CitedScalar", "Criterion", "DeltaVItem", "Maneuver", "MissionSizing", "Reference", "Requirement",
    "RequirementCheck", "ThrusterRecord", "assumption", "bipropellant_split", "burn_duration_screen",
    "check_requirements", "circular_spiral_delta_v", "delta_v_budget", "hohmann_transfer",
    "ideal_cold_gas_storage", "ideal_regulated_pressurant", "liquid_tank_volume", "list_thrusters",
    "nominal_operating_point", "perigee_lowering_delta_v", "power_resources", "reference",
    "size_mission", "thruster", "weighted_trade",
    "Candidate", "CandidateStudy", "ProjectStudy", "analyze_project", "project_from_dict",
    "project_from_json", "build_project_html",
    "DesignCase", "design_sensitivity",
]
