# Roadmap

This page is the compact milestone view. Sequencing, effort, dependencies,
acceptance gates, and the first three sprints are specified in
[`DEVELOPMENT_PLAN.md`](DEVELOPMENT_PLAN.md).

## Milestone 1 — foundation (implemented)

- Isentropic perfect-gas relations and inverse conversions
- Normal-shock relations and ideal-gas state calculator
- Dependency-free local web server and responsive UI

## Milestone 2 — connected gas dynamics (implemented)

- Static/stagnation stations and efficiency-corrected compression/expansion
- Fanno and Rayleigh lines
- Choked C-D nozzle mass flow, exit state, c*, thrust coefficient, and Isp
- Weak/strong oblique shocks and attached-shock detachment limit

## Milestone 3 — thermodynamics and rocket-system layer (implemented)

- Forward/inverse Tsiolkovsky equation and complete thrust decomposition
- Closed-system polytropic P-v path and first-law energy accounting
- Common liquid-bipropellant, solid, and hybrid screening catalog

## Milestone 4 — geometry and plots (implemented)

- Parametric conical and smooth bell-like axisymmetric nozzle contours
- Subsonic throat approach and supersonic diverging station solution
- Isentropic property, area-Mach, nozzle contour, Mach, pressure, temperature,
  P-v, and T-s plots
- Plot-ready sampled API responses with no duplicated physics in JavaScript

## Milestone 5 — higher-fidelity physics (implemented with declared gates)

- Temperature-dependent heat capacities and frozen species/mixture properties
  (implemented in 0.4.0)
- External NASA CEA adapter and experimental local equilibrium provider
- Frozen/equilibrium nozzle comparison; finite-rate chemistry remains out of scope
- Switchable boundary-layer, divergence, discharge, combustion, and two-phase bands
- Internal shocks and a source-labelled nozzle-separation screening envelope

## Milestone 6 — engineering studies (implemented in 0.9.0)

- Multi-stage delta-v and payload-fraction optimization
- Ambient-pressure/altitude sweeps
- CSV/SVG/DXF/STL export, saved cases, comparisons, and uncertainty bands
- Rao/MOC geometry, Bartz thermal loads, regenerative bulk cooling balance
- Pressure-fed, gas-generator, and staged-combustion component energy balances

## Milestone 7 — 1.0 release gate (solver gate closed 2026-10-25)

- Closed: nine NASA CEA 3.3.4 reference points across LOX/LH2, LOX/RP-1 and
  LOX/CH4, recorded with the CEA version and database hash; the local solver
  agrees within 2e-5 against a declared 1 % tolerance
  (`tests/reference_data/cea_gate_points.json`)
- Closed: tolerance/result matrix from the same CEA version
  ([`TOLERANCES.md`](TOLERANCES.md), [`VALIDATION_REPORT.md`](VALIDATION_REPORT.md))
- Remaining before the 1.0 tag: clean-machine package smoke test

## Milestone 8 — standalone burn analysis and optional Sidera hand-off (in progress)

- Implemented: standalone Python/API/CLI constant-burn analysis, explicit
  propellant-stream inventory closure, analytic targeting and total impulse
- Implemented: system Isp from total tank drain, reserve/dry-mass/duration
  cutoffs, versioned hashed artifacts and standalone HTML reports
- Implemented: exact constant-burn artifact/native Sidera adapter with
  capability negotiation and no base dependency on Sidera
- Implemented: constant bipropellant/monopropellant provider constructors,
  browser Burn workbench, versioned workspace persistence and HTML/JSON export
- Implemented: dependency-free L1 piecewise-linear start/throttle/stop
  schedules, analytic integrals, v2 persistence, thrust/mass plots and explicit
  equivalent-constant Sidera reduction evidence
- Implemented in isolated Python core: analytic first-order response lag,
  ignition delay, restart/cooldown/minimum-pulse validation, independently
  commanded engine-cluster composition and convergence-qualified tabulation
- Implemented in isolated Python core: discontinuity-safe
  `rocket_propulsion_tabulated_burn_v1` exchange artifact plus unsmoothed CSV
  thrust/flow import with explicit linear/hold interpolation and source hashes
- Implemented in isolated Python core: reference-backed pressure-fed blowdown
  with polytropic ullage decay, pressure/reserve cutoff, bounded correlation or
  measured pressure map, variable thrust/Isp/flow and mesh-convergence evidence
- Implemented: immutable primary-source registry, citation-bearing artifact
  hashes and automated `References`-section coverage for the public transient API
- Implemented in isolated Python core: seeded propulsion uncertainty ensembles,
  bounded marginal distributions, aleatory/epistemic separation, sourced
  correlation, percentile/statistical summaries, stability diagnostics and
  per-output screening sensitivity; every realization retains its burn artifact
- Implemented in isolated Python core: exhaustive probability-labelled
  propulsion scenarios for engine-out/common-cause loss, ignition delay, early
  cutoff, response change and thrust/Isp degradation, including exact
  zero-thrust outcomes and weighted propulsion consequence percentiles
- Implemented in isolated Python core: hierarchical discrete/continuous
  propulsion ensembles with per-scenario conditional definitions, exact joint
  weights, physical thrust/Isp/time/axial scaling, scenario contributions,
  retained artifacts and independently seeded replicate diagnostics
- Implemented: hashed `rocket_propulsion_hierarchical_ensemble_v1` persistence
  with nested artifact verification, semantic reproduction checks, and
  within-scenario/between-scenario law-of-total-variance attribution
- Implemented in the isolated adapter: hashed weighted Sidera hand-off jobs,
  explicit axial-force/total-tank-flow channel semantics, zero-thrust retention,
  and fail-closed native capability negotiation
- Implemented in the isolated propulsion core: hot-fire scale calibration with
  unbalanced/heteroscedastic Mandel-Paule random effects, explicit measurement
  uncertainty removal, engine/run variance components, solver evidence and
  ensemble-ready residual parameters
- Implemented in the isolated propulsion core: paired multivariate hot-fire
  calibration with complete measurement covariance, within/between-engine
  cross-covariance, recorded positive-definite regularization, ensemble-ready
  copulas and semantically revalidated v1 persistence
- Implemented in the isolated propulsion core: model-form and
  qualification-to-flight random-effects discrepancy calibration, predictive
  epistemic variance, application-domain/evidence-reuse gates, qualified
  multi-component scales, and deterministic-correction/ensemble composition
- Implemented in the isolated propulsion core: a propagation-ready transient
  force model with simultaneous axial acceleration, total/tank mass drains,
  explicit left/right discontinuities, start/boundary/stop/reserve event
  surfaces, bounded estimation drivers and analytic local Jacobians
- Implemented in the isolated propulsion core: complete rectilinear
  operating-condition performance maps with N-D multilinear interpolation,
  analytic condition Jacobians, hash-verified persistence, validation-domain
  event surfaces and exact propagator-bridge composition
- Implemented in the isolated propulsion core: six-state regulated pressurant
  and feed dynamics with bottle/ullage energy balances, choked regulator flow,
  valve lag, liquid line loss, protective events, converged feed/performance
  closure and implicit coupled Jacobians
- Implemented in the isolated propulsion core: a two-state dynamic liquid-line
  surrogate with hydraulic inertance/compliance, quadratic damping, explicit
  stored-mass closure, pressure/flow events and engine/feed-composed Jacobians
- Implemented in the isolated propulsion core: power-limited electric
  propulsion with hash-verified throttle tables, maximum-thrust/maximum-Isp
  throttle paths, up-switch hysteresis, a fixed-efficiency continuous law,
  deterministic power-switching roots and analytic power/direction partials
- Implemented in the isolated propulsion core: ambient back-pressure thrust in
  the propagation bridge (`F = F_vac - p_a A_e`, analytic partial, net-thrust
  margin root), navigation-facing driver priors linearized from hot-fire
  calibration, the Gates execution-error model, first-order Gauss-Markov
  thrust noise, and body-frame thruster force/torque with propellant-dependent
  centre of mass and inertia
- Implemented: fail-closed input validation and canonical numbers across the
  propagation interface (fuzz- and property-tested), reopenable documents for
  bridges, electric models and timelines, multi-burn timelines with jettison,
  and bounded-variable thrust allocation
- Implemented: throttle-commanded engine model with smooth (PCHIP) throttle
  characteristics, free-space burn metrics and impulsive-equivalent time,
  saltation support (`time_event_jump`) with an end-to-end variational test,
  and fixes plus regression tests for every finding of a code review
- Implemented: Lunaris adapter (exact-mass impulsive equivalents, finite burns
  in Lunaris' propagator, segmentation markers) verified end-to-end, C1 cubic
  Hermite performance surfaces, minimum-effort balanced thrust allocation, and
  fail-closed validation of the feed models and Sidera adapters
- Implemented (2026-10-25): native Sidera wrapper flying the propagation
  bridge through Sidera's public `propagate` (`integrations/sidera_native.py`),
  HTTP routes and browser modules for the transient core and SASS sizing,
  distributed single-pipe water hammer, regulator hysteresis with lockup,
  layered ullage stratification, and hot-fire calibration fields over the
  operating-condition surface
- Remaining: line priming (no primary model source read yet), named-tank
  states and operating-condition inputs in the Sidera wrapper (Sidera exposes
  only a 7-state mass model), Sidera-side trajectory/pointing/epoch dispersion
  (owned by Sidera) and broader validation grids
- Keep trajectory integration, maneuver frames, guidance and orbital results in
  Sidera
- Add transient-profile coupling only through a versioned thrust/mass-flow
  provider contract, with no silent constant-burn approximation

## Milestone 9 — thermodynamics and aerodynamics depth (implemented 2026-10-25)

- Real-fluid reference equations of state for H2, O2, N2, CH4 and He
  (NIST-verified) with cubic fallbacks
- NASA CEA transport properties and equilibrium thermodynamic derivatives
- Thermally perfect gas flow, Taylor–Maccoll cones, shock-expansion and
  linearized supersonic airfoils, compressibility corrections, hypersonic
  Newtonian methods
- Potential flow, NACA sections, thin airfoil, lifting panels and lifting line
- Boundary layers, Couette flow, stagnation and flat-plate heating, entry
- Launch-vehicle normal force, crossflow, drag build-up, loads and an
  SLS-style aerodynamic database with uncertainty structure
- Remaining: transition prediction, Fay–Riddell and radiative heating (sources
  not yet read), CFD-level fields (external), and the hardware evidence that
  TRL 5 requires ([`TRL_ASSESSMENT.md`](TRL_ASSESSMENT.md))

## Milestone 10 — heat transfer and combustion depth (implemented 2026-10-25)

- Conduction: steady and transient multilayer walls, fins, analytic transient
  solutions, an energy-audited finite-volume solver
- Convection: internal, external, natural, boiling and condensation
  correlations; Taylor and Hess–Kunz hydrogen; Hatch–Papell film cooling
- Radiation: Planck and band fractions, view factors, gray enclosures,
  radiation-cooled extensions, H₂O–CO₂–soot WSGG
- Coupled regenerative-cooling solver with real-fluid coolant transport
- Combustion: finite-area combustor, products and radicals, dissociation
  energy, efficiency models, Isp loss chain, O/F optimisation, condensed-phase
  equilibrium, finite-rate H₂–O₂ nozzle kinetics
- Six HTTP routes and three browser views; cross-model consistency tests
- Added: transpiration cooling (film theory, Eckert–Livingood, Rannie, porous
  wall), a gas-side comparison with NASA TN D-6396 firings, the Priem–Heidmann
  vaporization length and curve, and SP-194 combustion-instability relations
- Added: ethanol and MMH flow-boiling burnout, feed-coupled chugging, a
  TN D-3532 total heat-load comparison, and HTTP routes and workbench panes
  for transpiration, stability and vaporization
- Remaining: station-by-station coolant-side hot-fire validation, RP-1
  burnout, hydrocarbon kinetics, pump-coupled and nonlinear combustion
  instability

The reviewed interface, ownership boundary, phases and acceptance gates are in
[`SIDERA_BURN_INTEGRATION_PLAN.md`](SIDERA_BURN_INTEGRATION_PLAN.md).

