"""Real-fluid thermodynamics for propellants and pressurants.

Reference equations of state explicit in the reduced Helmholtz energy, with
literature coefficients for parahydrogen and normal hydrogen [1]_, oxygen
[2]_, nitrogen [3]_, methane [4]_ and helium [5]_, plus Peng-Robinson and
Soave-Redlich-Kwong cubic equations [6]_ as a lower-accuracy fallback for
fluids without a reference equation.

Example
-------
>>> from rocket_propulsion.thermodynamics.real_fluids import reference_fluid
>>> lox = reference_fluid("oxygen")
>>> state = lox.state_tp(90.0, 1.0e6)
>>> round(state.density_kg_m3, 1)  # NIST WebBook (accessed 2026-10-25): 1144.0956
1144.1

The submodules are :mod:`.helmholtz` (framework and solvers),
:mod:`.reference_fluids` (coefficients and validity ranges) and
:mod:`.cubic` (cubic fallbacks). See ``docs/REAL_FLUIDS.md`` for sources,
validity ranges and the measured agreement with NIST WebBook tables.

References
----------
.. [1] LEACHMAN-2009-H2. Leachman et al., J. Phys. Chem. Ref. Data 38,
       721-748 (2009), https://srd.nist.gov/jpcrdreprint/1.3160306.pdf
.. [2] STEWART-1991-O2. Stewart, Jacobsen and Wagner, J. Phys. Chem. Ref.
       Data 20, 917-1021 (1991), https://srd.nist.gov/jpcrdreprint/1.555897.pdf
.. [3] SPAN-2000-N2. Span et al., J. Phys. Chem. Ref. Data 29, 1361-1433
       (2000), https://srd.nist.gov/jpcrdreprint/1.1349047.pdf
.. [4] SETZMANN-WAGNER-1991-CH4. Setzmann and Wagner, J. Phys. Chem. Ref.
       Data 20, 1061-1155 (1991),
       https://srd.nist.gov/jpcrdreprint/1.555898.pdf
.. [5] ORTIZVEGA-2023-HE. Ortiz Vega et al., NIST IR 8474 (2023),
       https://doi.org/10.6028/NIST.IR.8474
.. [6] BELL-JAGER-2016-CUBIC. Bell and Jaeger, J. Res. NIST 121, 238-263
       (2016), https://doi.org/10.6028/jres.121.011
"""

from .cubic import CubicResidual, constant_cp_ideal, cubic_counterpart, cubic_model
from .helmholtz import (
    DepartureFunctions,
    FluidModel,
    FluidState,
    IdealHelmholtz,
    LiquidDensityAncillary,
    RealFluid,
    ResidualHelmholtz,
    ResidualTerm,
    SaturationState,
    VaporPressureAncillary,
    with_reference_state,
)
from .reference_fluids import (
    REFERENCE_FLUID_KEYS,
    ReferenceFluidInfo,
    reference_fluid,
    reference_fluid_info,
    reference_fluid_model,
)

from .transport import (
    TRANSPORT_FLUID_KEYS,
    TransportProperties,
    TransportValidity,
    OlchowySengersParameters,
    transport_tp,
    transport_rho_t,
    transport_from_state,
    thermal_conductivity_rho_t,
    viscosity_rho_t,
    dilute_gas_transport,
    transport_validity,
)

__all__ = [
    "transport_validity",
    "dilute_gas_transport",
    "viscosity_rho_t",
    "thermal_conductivity_rho_t",
    "transport_from_state",
    "transport_rho_t",
    "transport_tp",
    "OlchowySengersParameters",
    "TransportValidity",
    "TransportProperties",
    "TRANSPORT_FLUID_KEYS",
    "REFERENCE_FLUID_KEYS",
    "CubicResidual",
    "DepartureFunctions",
    "FluidModel",
    "FluidState",
    "IdealHelmholtz",
    "LiquidDensityAncillary",
    "RealFluid",
    "ReferenceFluidInfo",
    "ResidualHelmholtz",
    "ResidualTerm",
    "SaturationState",
    "VaporPressureAncillary",
    "constant_cp_ideal",
    "cubic_counterpart",
    "cubic_model",
    "reference_fluid",
    "reference_fluid_info",
    "reference_fluid_model",
    "with_reference_state",
]
