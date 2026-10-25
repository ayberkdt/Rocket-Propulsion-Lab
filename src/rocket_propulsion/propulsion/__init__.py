"""Rocket-system performance equations, engine-cycle, staging and loss models.

The package provides the ideal rocket and steady thrust equations
(:mod:`.performance`), preliminary engine-cycle balances (:mod:`.cycles`),
finite-grid staging optimization (:mod:`.staging`), declared performance-factor
budgets (:mod:`.losses`), conditional two-phase exit-state ranges
(:mod:`.two_phase`), propellant screening records (:mod:`.propellants`) and the
finite-burn interfaces re-exported from :mod:`.burns`. Each module documents its
model assumptions and references.
"""

from .burns import (
    BurnDefinition,
    BurnResult,
    BurnTarget,
    BurnTargetKind,
    ConstantBurnReduction,
    ConstantPerformanceProvider,
    EngineOperatingPoint,
    FlowDestination,
    ProfiledBurnPoint,
    ProfiledBurnResult,
    ProfiledPerformanceProvider,
    PropellantFlow,
    PropellantRole,
    PropulsionProvider,
    TankDefinition,
    ThrottleSchedule,
    ThrottleSegment,
    build_burn_html_report,
    burn_result_from_dict,
    burn_result_from_json,
    burn_result_to_dict,
    burn_result_to_json,
    reduce_to_constant,
    simulate_burn,
    simulate_constant_burn,
    simulate_profiled_burn,
)
from .cycles import EngineCycleInputs, EngineCycleResult, calculate_engine_cycle
from .losses import LossBudgetInputs, LossTerm, PerformanceLossBudget, calculate_loss_budget
from .performance import (
    STANDARD_GRAVITY_M_S2,
    RocketEquationResult,
    ThrustResult,
    calculate_rocket_equation,
    calculate_thrust,
)
from .propellants import PROPELLANTS, PropellantRecord, list_propellants
from .staging import (
    StageDefinition,
    StageSolution,
    StagingInputs,
    StagingOptimizationResult,
    optimize_staging,
)
from .two_phase import TwoPhaseInputs, TwoPhasePerformanceBand, calculate_two_phase_band

__all__ = [
    "PROPELLANTS",
    "STANDARD_GRAVITY_M_S2",
    "BurnDefinition",
    "BurnResult",
    "BurnTarget",
    "BurnTargetKind",
    "ConstantBurnReduction",
    "ConstantPerformanceProvider",
    "EngineCycleInputs",
    "EngineCycleResult",
    "EngineOperatingPoint",
    "FlowDestination",
    "LossBudgetInputs",
    "LossTerm",
    "PerformanceLossBudget",
    "ProfiledBurnPoint",
    "ProfiledBurnResult",
    "ProfiledPerformanceProvider",
    "PropellantFlow",
    "PropellantRecord",
    "PropellantRole",
    "PropulsionProvider",
    "RocketEquationResult",
    "StageDefinition",
    "StageSolution",
    "StagingInputs",
    "StagingOptimizationResult",
    "TankDefinition",
    "ThrottleSchedule",
    "ThrottleSegment",
    "ThrustResult",
    "TwoPhaseInputs",
    "TwoPhasePerformanceBand",
    "build_burn_html_report",
    "burn_result_from_dict",
    "burn_result_from_json",
    "burn_result_to_dict",
    "burn_result_to_json",
    "calculate_engine_cycle",
    "calculate_loss_budget",
    "calculate_rocket_equation",
    "calculate_thrust",
    "calculate_two_phase_band",
    "list_propellants",
    "optimize_staging",
    "reduce_to_constant",
    "simulate_burn",
    "simulate_constant_burn",
    "simulate_profiled_burn",
]
