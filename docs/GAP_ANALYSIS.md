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
| Off-design with chemistry | Implemented | Perfect gas calibrated to the equilibrium design point; exact mass flow and design exit pressure |
| Shock trains / nozzle separation | Partial | Internal normal shock; Summerfield screen and Schmucker's criterion (NASA TM-77396), which locates incipient separation on the 1-D wall pressure; no shock train, side loads or viscous separation field |
| Boundary layer / heat transfer | Implemented (1-D) | Coupled gas/wall/coolant regenerative solve (Bartz at the local wall temperature, fin lands, Taylor or Gnielinski coolant side, real-fluid coolant, optional WSGG gas radiation); transient multilayer conduction; convection and radiation libraries (see the heat-transfer section); no axial or circumferential wall conduction, no two-phase coolant, no CFD |
| Method of characteristics | Partial | Local interior, boundary, axis-departure and straight-wall field processes have targeted checks; the legacy full-nozzle generator uses axisymmetric area mapping and has unresolved characteristic-geometry defects (`GEOMETRY_REVIEW.md`) |

## Thermodynamics

| Capability | Status | Current boundary |
|---|---|---|
| Ideal-gas state and caloric properties | Implemented | Constant cp, cv, gamma |
| Static/stagnation connection | Implemented | Adiabatic and reversible |
| Component efficiency process | Implemented | Lumped compressor/turbine/nozzle efficiency |
| Polytropic P-v path and first law | Implemented | Quasistatic closed system |
| Variable heat capacity | Implemented | 36 NASA species (TP-2002-211556 plus pinned CEA `thermo.inp`), no extrapolation |
| Frozen ideal-gas mixtures | Implemented | 36 H/C/N/O species, ideal mixing |
| Combustion equilibrium | Implemented | Local Gibbs solver with condensed liquid reactants, gated against NASA CEA 3.3.4 (nine points plus 0.5–25 MPa); optional NASA CEA provider; gas-phase products only |
| Rocket performance (IAC) | Implemented | Shifting or frozen expansion, maximum-mass-flux throat; matches CEA within 2e-5 |
| Condensed-phase products (Al₂O₃, carbon) | Implemented | Graphite, ice/water and liquid/α-alumina in TP/HP equilibrium and IAC rocket expansion, within 9e-6 of CEA 3.3.4; no particle velocity or thermal lag, no FAC with condensed phases |
| Dissociation / finite-rate chemistry | Implemented | Dissociation energy, radicals and recombination heat per station; finite-area combustor; quasi-1-D finite-rate nozzle with the GRI-Mech 3.0 H₂–O₂ subset and Bray freezing (LOX/LH₂ only; rate set not validated at rocket conditions) |
| Two-phase flow | Partial | Condensed-fraction coupling band only; no particle-resolved flow |
| Real-fluid propellant properties | Implemented | Reference Helmholtz EOS for H2 (para/normal), O2, N2, CH4, He within 11 ppm of NIST; PR/SRK cubic fallback; no melting line or mixtures |
| Transport properties | Implemented | NASA CEA trans.inp fits and mixture rules (CEA 3.3.4 reproduced); Sutherland and Chapman–Enskog fallbacks |
| Equilibrium derivatives | Implemented | Frozen/equilibrium cp, cv, γs, sound speed per RP-1311 Part I; matches CEA to ~1e-8 |

## Geometry and visualization

| Capability | Status | Current boundary |
|---|---|---|
| Conical nozzle | Implemented | Analytical axisymmetric contour |
| Smooth bell-like nozzle | Implemented | Hermite preliminary approximation |
| Station-by-station flow | Implemented | Isentropic branch selected from geometry |
| Engineering plots | Implemented | Responsive SVG, local and dependency-free |
| Rao optimum / MOC contour | Partial | Rao/TOP approximation and legacy mapped mesh; no verified optimized full-nozzle characteristic field; mesh self-refinement is not independent validation |
| CAD/mesh export | Implemented | CSV/SVG/DXF/STL; preliminary engineering exchange only |

## External aerodynamics

| Capability | Status | Current boundary |
|---|---|---|
| Thermally perfect gas flow | Implemented | Isentropic, normal and oblique shocks with vibrational excitation; no dissociation |
| Cone flow | Implemented | Taylor–Maccoll, zero incidence, matches NASA SP-3004 |
| Supersonic airfoils | Implemented | Shock-expansion, Ackeret and Busemann; attached shocks |
| Compressibility corrections | Implemented | Prandtl–Glauert, Kármán–Tsien, Laitone, critical Mach |
| Hypersonic inviscid | Implemented | Newtonian, modified Newtonian, strong-shock limits |
| Subsonic airfoils and wings | Implemented | Thin airfoil, lifting panel method, lifting line; inviscid (≈14 % above measured section lift) |
| Boundary layers and heating | Implemented | Blasius, reference temperature, Van Driest II, Couette, Tauber and Sutton–Graves; no transition prediction, no Fay–Riddell, no radiative heating |
| Atmospheric entry | Implemented | Allen–Eggers and planar point-mass integration below 86 km |
| Launch-vehicle aerodynamics | Implemented | Barrowman, Allen–Perkins crossflow, fin–body interference, SP-8035 loads, SLS-style database and uncertainty structure; vehicle coefficients from CFD/wind tunnel are user data |
| CFD-level flow fields | Open | Out of scope; supply CFD or wind-tunnel databases to `aero_database` |

## Propellant data and propulsion system

| Capability | Status | Current boundary |
|---|---|---|
| Propellant catalogue | Implemented | Seven families; four linked to computable pairs, three state why they are not computed (condensed products, HTPB and paraffin absent from NASA data) |
| Turbopump | Implemented | Affinity laws, specific/suction speed, NPSH margin and a homologous head–flow map; coefficients are caller data |
| Distributed water hammer | Implemented | Single-pipe MOC with friction, Joukowsky check; no column separation or fluid–structure interaction |
| Line priming | Open | Sources read give only qualitative priming guidance |
| Regulator hysteresis | Implemented | Quasi-static force balance with Coulomb-friction band and lockup |
| Ullage stratification | Implemented | Layered ullage with natural-convection wall/interface heat transfer; no condensation |
| Calibration fields | Implemented | Hot-fire calibrations interpolated over the operating-condition surface with a conservative uncertainty envelope |
| Transient HTTP/UI surfaces | Implemented | Ten routes and five browser modules |
| Global sensitivity | Implemented | Spearman, PRCC, SRRC with rank-regression R² and Sobol first-order/total indices (Saltelli design, bootstrap intervals) in `studies.global_sensitivity`, exposed on the Δv uncertainty route; the persisted burn-uncertainty documents keep their Pearson screening records; correlated-input Sobol indices and quasi-random designs open; the blowdown uncertainty route also returns Spearman, PRCC and SRRC per metric |

## Heat transfer and combustion

| Capability | Status | Current boundary |
|---|---|---|
| Conduction | Implemented | 1-D steady and transient (finite-volume θ-method) multilayer walls with k(T), contact resistance, generation, convective and radiative boundaries; fins with exact Bessel functions; no 2-D ribs or axial conduction |
| Convection | Implemented | Internal, entry, rough, liquid-metal, external and natural convection; Taylor and Hess–Kunz hydrogen, RTE coolant form with curvature; boiling, peak flux and film condensation; Hatch–Papell film cooling; SP-8087 burnout-margin screen with subcooled flow-boiling burnout for ethanol and MMH; no free RP-1, N₂O₄ or hydrazine burnout correlation found |
| Transpiration cooling | Implemented | Film theory (TN 3208) with the wall energy balance, wall efficiency and radiation; Eckert–Livingood and Rannie relations; two-temperature porous wall; inside the spread of four published methods on Kelley's case; no rocket-firing comparison, no foreign-gas diffusion term |
| Radiation | Implemented | Planck, band fractions, closed-form view factors, gray enclosures, radiation-cooled extensions; H₂O–CO₂ WSGG (Bordbar 2020) and soot WSGG (Sadeghi 2021) by superposition, fitted at 1 atm; no alumina particle radiation, no spectral line-by-line model |
| Regenerative cooling | Implemented (1-D gas and coolant, fin or 2-D wall) | Coupled station march with axially varying, helical or bifurcated channels, RTE friction plus Euler momentum pressure drop and GFSSP sudden-area-change losses, and an optional two-dimensional wall cross-section at every station, verified against its components and conservation; gas side compared with NASA TN D-6396 firings (Bartz 7–30 % high in the chamber, 13–61 % at the throat) and total heat load with TN D-3532 (bracketed for 1100–1500 K walls); station-by-station coolant side not compared with a firing |
| Combustion products and dissociation | Implemented | Composition to 1e-14, radicals, dissociation energy and recombination heat, within 2.5e-4 of CEA down to 1e-10 mole fraction |
| Finite-area combustor | Implemented | RP-1311 FAC within 1.3e-4 of CEA; gas phase only |
| Combustion efficiency | Implemented | Measured efficiencies, vaporisation-limited c*, Rupe mixing and JANNAF ERE bookkeeping, SP-8089 mixing, L*; Priem–Heidmann generalized length and curve as tabulated in NASA CR-168066 (TR R-67 itself not on NTRS) |
| Finite-rate kinetics | Implemented (H₂–O₂) | GRI-Mech 3.0 subset with RODAS; LOX/LH₂ only; chamber-start mode does not resolve fast throat chemistry; no published ODK/TDK case compared |
| Combustion instability | Implemented (linear, closed form) | SP-194 chamber modes, chugging and Crocco–Cheng limits, n–τ neutral curve, resonator tuning, response factor, Q and decay rate; feed-coupled chugging (lines, manifolds, orifices, tanks; Nyquist); no pumps, baffles, liners in the eigenvalue or nonlinear instability |

## Highest-value next steps

1. Hot-fire coolant-side and wall-temperature data (for example NASA TN D-7207 or a calorimeter chamber) to validate the coupled regenerative solver beyond the gas side.
2. Sourced HTPB and paraffin thermochemistry to compute the hybrid families.
3. Line priming once a primary model source is available.
4. Clean-machine (not only clean-environment) package smoke test for the 1.0 tag.
