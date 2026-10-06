"""Thermodynamic state relations used in propulsion calculations."""

from .flow_state import (
    IsentropicProcessResult,
    StagnationState,
    calculate_isentropic_process,
    calculate_stagnation_state,
)
from .heat_transfer import (
    BartzHeatTransferResult,
    HeatFluxStation,
    calculate_bartz_heat_transfer,
)
from .ideal_gas import IdealGasState, complete_ideal_gas_state
from .polytropic import PolytropicPoint, PolytropicResult, calculate_polytropic_process
from .regenerative import RegenerativeCoolingResult, calculate_regenerative_cooling

__all__ = [
    "BartzHeatTransferResult",
    "HeatFluxStation",
    "IdealGasState",
    "IsentropicProcessResult",
    "PolytropicPoint",
    "PolytropicResult",
    "RegenerativeCoolingResult",
    "StagnationState",
    "calculate_bartz_heat_transfer",
    "calculate_isentropic_process",
    "calculate_polytropic_process",
    "calculate_regenerative_cooling",
    "calculate_stagnation_state",
    "complete_ideal_gas_state",
]
