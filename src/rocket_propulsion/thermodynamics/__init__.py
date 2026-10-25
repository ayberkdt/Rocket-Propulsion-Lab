"""Thermodynamic state relations used in propulsion calculations."""

from . import real_fluids
from .flow_state import (
    IsentropicProcessResult,
    StagnationState,
    calculate_isentropic_process,
    calculate_stagnation_state,
)
from .heat_transfer import (
    BartzContour,
    BartzContourStation,
    BartzHeatTransferResult,
    BartzInputs,
    HeatFluxStation,
    calculate_bartz_heat_transfer,
)
from .ideal_gas import IdealGasState, complete_ideal_gas_state
from .materials import (
    CitedProperty,
    CoolantRecord,
    PropertyCitation,
    WallMaterialRecord,
    coolant_catalog,
    get_coolant,
    get_wall_material,
    wall_material_catalog,
)
from .polytropic import PolytropicPoint, PolytropicResult, calculate_polytropic_process
from .regenerative import (
    CoolantSelection,
    RegenerativeCoolingResult,
    RegenerativeInputs,
    calculate_regenerative_cooling,
)
from .wall_conduction import (
    WallConductionInputs,
    WallConductionResult,
    calculate_wall_conduction,
)

__all__ = [
    "real_fluids",
    "BartzContour",
    "BartzContourStation",
    "BartzHeatTransferResult",
    "BartzInputs",
    "CitedProperty",
    "CoolantRecord",
    "CoolantSelection",
    "HeatFluxStation",
    "IdealGasState",
    "IsentropicProcessResult",
    "PolytropicPoint",
    "PolytropicResult",
    "PropertyCitation",
    "RegenerativeCoolingResult",
    "RegenerativeInputs",
    "StagnationState",
    "WallConductionInputs",
    "WallConductionResult",
    "WallMaterialRecord",
    "calculate_bartz_heat_transfer",
    "calculate_isentropic_process",
    "calculate_polytropic_process",
    "calculate_regenerative_cooling",
    "calculate_stagnation_state",
    "calculate_wall_conduction",
    "complete_ideal_gas_state",
    "coolant_catalog",
    "get_coolant",
    "get_wall_material",
    "wall_material_catalog",
]
