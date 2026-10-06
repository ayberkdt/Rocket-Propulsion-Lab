# Physics and product gap analysis

This file is the explicit boundary between what the current application solves
and what would require a higher-fidelity model.

## Compressible aerodynamics

| Capability | Status | Current boundary |
|---|---|---|
| Isentropic perfect-gas flow | Implemented | Constant gamma, one-dimensional |
| Area-Mach inversion | Implemented | Explicit subsonic/supersonic branch |
| Normal and oblique shocks | Implemented | Attached planar waves, perfect gas |
| Fanno and Rayleigh lines | Implemented | Constant area and constant gamma |
| Choked C-D nozzle | Implemented | Inviscid, frozen-property design point |
| Nozzle station solution | Implemented | Quasi-one-dimensional, no internal shock |
| Shock trains / nozzle separation | Partial | Internal normal shock and Summerfield screen; no shock train or viscous separation field |
| Boundary layer / heat transfer | Partial | Bartz and lumped boundary-layer loss; no conjugate heat transfer or CFD |
| Method of characteristics | Implemented | Planar compatibility with first-order axisymmetric area mapping |

## Thermodynamics

| Capability | Status | Current boundary |
|---|---|---|
| Ideal-gas state and caloric properties | Implemented | Constant cp, cv, gamma |
| Static/stagnation connection | Implemented | Adiabatic and reversible |
| Component efficiency process | Implemented | Lumped compressor/turbine/nozzle efficiency |
| Polytropic P-v path and first law | Implemented | Quasistatic closed system |
| Variable heat capacity | Implemented | NASA Glenn subset, 200–6000 K, no extrapolation |
| Frozen ideal-gas mixtures | Implemented | Ten H/C/N/O combustion species, ideal mixing |
| Combustion equilibrium | Partial | Local experimental Gibbs solver and external CEA adapter; nine-point CEA gate remains |
| Dissociation / finite-rate chemistry | Open | Required at rocket chamber temperatures |
| Two-phase flow | Partial | Condensed-fraction coupling band only; no particle-resolved flow |

## Geometry and visualization

| Capability | Status | Current boundary |
|---|---|---|
| Conical nozzle | Implemented | Analytical axisymmetric contour |
| Smooth bell-like nozzle | Implemented | Hermite preliminary approximation |
| Station-by-station flow | Implemented | Isentropic branch selected from geometry |
| Engineering plots | Implemented | Responsive SVG, local and dependency-free |
| Rao optimum / MOC contour | Implemented | Rao/TOP approximation plus MOC mesh and refinement evidence |
| CAD/mesh export | Implemented | CSV/SVG/DXF/STL; preliminary engineering exchange only |

## Highest-value next steps

1. Add temperature-dependent species and mixture properties.
2. Integrate a NASA CEA-compatible equilibrium adapter.
3. Add back-pressure nozzle regimes, internal shocks, and separation warnings.
4. Replace the bell approximation with an MOC/Rao contour option.
5. Add study comparison, export, and uncertainty propagation.

