"""Compressible-flow relations for calorically perfect gases.

The package provides one-dimensional isentropic, normal-shock, oblique-shock,
Prandtl-Meyer [1]_, Fanno and Rayleigh relations, ideal and off-design
converging-diverging nozzle performance [2]_, and the U.S. Standard Atmosphere
1976 pressure used for altitude sweeps. Each public callable documents its own
assumptions and sources.

References
----------
.. [1] NACA-TR-1135. Ames Research Staff, "Equations, Tables, and Charts for
       Compressible Flow", NACA Report 1135, 1953.
.. [2] SUTTON-BIBLARZ-2017. G. P. Sutton and O. Biblarz, "Rocket Propulsion
       Elements", 9th ed., Wiley, 2017.
"""

from .duct_flow import FannoResult, RayleighResult, calculate_fanno, calculate_rayleigh
from .isentropic import IsentropicResult, calculate_isentropic, mach_from_isentropic
from .normal_shock import NormalShockResult, calculate_normal_shock
from .nozzle import NozzleResult, calculate_ideal_nozzle
from .oblique_shock import ObliqueShockResult, calculate_oblique_shock
from .off_design_nozzle import (
    NozzleAltitudePoint,
    NozzleAltitudeSweep,
    NozzleFlowStation,
    NozzleRegimeThresholds,
    OffDesignNozzleResult,
    SeparationAssessment,
    ShockLocation,
    calculate_off_design_nozzle,
    generate_nozzle_altitude_sweep,
    nozzle_regime_thresholds,
    standard_atmosphere_pressure_pa,
)
from .sweeps import IsentropicSweepPoint, generate_isentropic_sweep

__all__ = [
    "FannoResult",
    "IsentropicResult",
    "IsentropicSweepPoint",
    "NormalShockResult",
    "NozzleAltitudePoint",
    "NozzleAltitudeSweep",
    "NozzleFlowStation",
    "NozzleRegimeThresholds",
    "NozzleResult",
    "ObliqueShockResult",
    "OffDesignNozzleResult",
    "RayleighResult",
    "SeparationAssessment",
    "ShockLocation",
    "calculate_fanno",
    "calculate_ideal_nozzle",
    "calculate_isentropic",
    "calculate_normal_shock",
    "calculate_oblique_shock",
    "calculate_off_design_nozzle",
    "calculate_rayleigh",
    "generate_isentropic_sweep",
    "generate_nozzle_altitude_sweep",
    "mach_from_isentropic",
    "nozzle_regime_thresholds",
    "standard_atmosphere_pressure_pa",
]

