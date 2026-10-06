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

## Milestone 7 — 1.0 release gate

- Close nine external CEA reference points across LOX/LH2, LOX/RP-1, and LOX/CH4
- Publish the final tolerance/result matrix from the same CEA version
- Tag 1.0 only after the external gate and clean-machine package smoke test pass

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
- Remaining: public HTTP API/UI surfaces for the isolated transient core,
  feed-line/regulator/thermal pressurant fidelity, model-form and
  operating-condition calibration, Sidera-side trajectory/pointing/epoch dispersion and broader
  validation grids
- Keep trajectory integration, maneuver frames, guidance and orbital results in
  Sidera
- Add transient-profile coupling only through a versioned thrust/mass-flow
  provider contract, with no silent constant-burn approximation

The reviewed interface, ownership boundary, phases and acceptance gates are in
[`SIDERA_BURN_INTEGRATION_PLAN.md`](SIDERA_BURN_INTEGRATION_PLAN.md).

