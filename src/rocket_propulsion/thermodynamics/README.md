# Thermodynamics

The first module models a calorically perfect ideal gas.

`complete_ideal_gas_state` solves `p = rho R T` from any two state variables and
derives speed of sound, specific enthalpy, and specific internal energy. Future
modules can add variable heat capacity, mixtures, combustion, and equilibrium
adapters without changing the API boundary.

`flow_state.py` links static and stagnation properties through Mach number and
models efficiency-corrected isentropic compression or expansion. This creates a
shared thermodynamic layer for inlets, compressors, turbines, chambers, and
nozzles.

`polytropic.py` adds a closed-system `p v^n = constant` path with boundary work,
heat transfer, internal-energy, enthalpy, entropy, and sampled P-v/T-s data for
thermodynamic diagrams.

`heat_transfer.py` evaluates a source-labelled Bartz gas-side heat-transfer
correlation over a shared nozzle contour. Gas properties are explicit inputs;
the model does not silently invent viscosity, heat capacity, or Prandtl number.
`regenerative.py` is an experimental constant-property, single-phase bulk
coolant energy balance with Darcy-Weisbach pressure loss. Required coolant and
channel data must all be supplied.

