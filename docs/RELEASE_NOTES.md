# Release notes

## Unreleased — standalone burn vertical slice

- Added immutable burn, tank, explicit propellant-stream, operating-point,
  profile, summary, and integration-evidence contracts.
- Added exact constant-burn solutions for duration, propellant, impulse, and
  ideal one-dimensional delta-v targets with protective inventory cutoffs.
- Added deterministic `rocket_propulsion_burn_v1` persistence with hash and
  physical reproduction validation.
- Added Python, HTTP, CLI, and standalone HTML-report surfaces.
- Added closed bipropellant/monopropellant operating-point providers, including
  explicit open-cycle dump streams and chamber-versus-system Isp accounting.
- Added a browser Burn workbench with inventory/mass history, canonical JSON,
  printable report, `.rplab.json` workspace and Sidera export flows.
- Added optional exact constant-burn Sidera artifact/native export and a public
  capability probe; Sidera remains outside base dependencies.
- Added L1 piecewise-linear ignition, throttle and shutdown schedules with
  analytic exposure, impulse, tank drain, delta-v and force-centroid integrals.
- Added reproducible `rocket_propulsion_burn_v2` artifacts, CLI/API/browser
  profile flows and mass/thrust charts.
- Added fail-closed transient Sidera export plus an explicit
  `equivalent_constant` mode that preserves duration, impulse and propellant
  while carrying approximation metrics and warnings.
- Added an isolated first-order engine-response solver with analytic exposure
  and force-centroid integrals, including ignition delay and shutdown settling.
- Added minimum pulse, restart-count, overlap and cooldown validation for
  independently commanded engines.
- Added engine-cluster transient composition with cant efficiency, staggered
  starts, explicit tank drains and convergence evidence against exact response
  integrals.
- Added the target-neutral `rocket_propulsion_tabulated_burn_v1` contract. It
  preserves discontinuities through independent segment endpoints and carries
  deterministic source/artifact hashes without importing or modifying Sidera.
- Added unsmoothed CSV thrust/flow replay with explicit linear or hold
  interpolation for solid-motor and qualification-test histories.
- Added a preliminary pressure-fed blowdown solver with polytropic ullage
  expansion, pressure/reserve cutoff, variable thrust/Isp/tank flow and
  convergence-qualified tabulation.
- Added interchangeable pressure-performance power laws and measured maps;
  both refuse extrapolation and preserve explicit stream/system-Isp closure.
- Added a reviewed primary-source registry, citation identifiers inside
  tabulated artifact hashes, and an automated public-docstring reference gate.
- Added propulsion-specific uncertainty propagation with bounded uniform,
  triangular, truncated-normal and log-uniform marginals; explicit
  aleatory/epistemic labels; and optional sourced Gaussian-copula correlation.
- Added deterministic LHS/Monte Carlo ensembles in which every draw reruns the
  blowdown convergence gate and retains a complete thrust/mass-flow artifact.
- Added mean/standard-error intervals, empirical P05/P50/P95 ranges,
  half/full-sample stability evidence and per-output Pearson sensitivity ranks.
- Added a discrete propulsion scenario layer for ignition/engine-out,
  common-cause loss, delay, early cutoff, response change and thrust/Isp
  degradation without implicit independence assumptions.
- Added exhaustive probability validation, exact all-engine-out artifacts and
  probability-weighted propulsion consequence percentiles for later Sidera
  trajectory ensembles.
- Added hierarchical composition of discrete outcomes and scenario-conditional
  continuous uncertainty. Each realization carries an exact joint weight and
  complete tabulated artifact; zero-thrust mass is preserved and probabilities
  are checked against double application.
- Added physical thrust/Isp/time/axial scaling, per-scenario expected-value
  contributions, deterministic provenance including uncertainty evidence, and
  independently seeded impulse-stability replicates.
- Added the hashed `rocket_propulsion_hierarchical_ensemble_v1` exchange
  document with full JSON round trip, nested artifact integrity checks and
  semantic reproduction of probabilities, metrics and evidence.
- Added law-of-total-variance attribution so continuous within-scenario scatter
  and discrete between-scenario risk are reported separately with a closure
  residual.
- Added `rocket_propulsion_sidera_ensemble_handoff_v1`, binding every weighted
  artifact to Sidera start time, normalized direction and frame without
  modifying or importing Sidera.
- Fixed constant Sidera export for canted propulsion: the adapter now sends net
  axial thrust with rocket-equivalent Isp, preserving both trajectory force and
  total tank drain. Previous delivered-thrust export overstated net force when
  cant efficiency was below one.
- Added capability gates for tabulated burns, independent mass flow and weighted
  maneuver ensembles; current Sidera correctly fails closed for this future
  native path.
- Added hot-fire scale calibration using a heteroscedastic Mandel-Paule
  random-effects model. Repeated runs now separate engine-to-engine,
  within-engine and declared measurement uncertainty instead of collapsing them
  into an unsupported percentage.
- Added deterministic calibration hashes, solver convergence evidence,
  operating-condition pooling refusal, NIST reference reproduction, and
  population/same-engine residual distributions for conditional ensembles.
- Added paired multivariate hot-fire calibration. Full per-run measurement
  covariance is removed from process cross-covariance; marginal Mandel-Paule
  fits and ordered Gaussian-copula correlation now come from the same campaign.
- Added explicit correlation clipping/shrinkage evidence and the semantically
  recomputed `rocket_propulsion_hot_fire_multivariate_calibration_v1` schema.
- Added model-form and qualification-to-flight discrepancy calibration from
  independent validation ratios, including consensus correction, predictive
  epistemic variance and a semantically recomputed v1 artifact.
- Added qualified scale components so aleatory hot-fire and epistemic
  discrepancy factors can act on the same thrust/Isp/time channel without
  losing classification. Evidence reuse and uncorrelated shared evidence now
  fail closed.
- Centralized deterministic and sampled artifact scaling so thrust, Isp, time,
  axial force, total tank drain and named stream-flow identities use one path.
- Added `rocket_propulsion_propagator_bridge_v1`: a framework-neutral transient
  force model returning acceleration, total and named tank mass derivatives,
  exact event surfaces, bounded estimation drivers, and analytic Jacobians.
- Added dry-mass/tank-reserve inhibition and explicit left/right evaluation at
  every profile discontinuity; no force or flow is extrapolated outside a burn.
- Added hash-verified N-dimensional operating-condition performance surfaces
  for supply pressure, ambient pressure, throttle, mixture ratio, power or
  other named coordinates, with a complete-grid gate and no extrapolation.
- Coupled surface force/flow scales into the propagator bridge and added exact
  condition Jacobians plus lower/upper validation-domain event surfaces.
- Added `rocket_propulsion_regulated_feed_system_v1`: finite pressurant bottle,
  growing thermal ullage, choked/unchoked regulator, valve lag, wall heat
  transfer, propellant inventory and quadratic feed-line loss ODEs.
- Added positive-safe feed cutoff roots, conservation/energy evidence, analytic
  inlet-pressure derivatives and a fail-closed iterative engine/feed closure.
- Added implicit coupled Jacobians through the pressure↔mass-flow algebraic
  loop, including an exposed singularity denominator for estimator safety.
- Added `rocket_propulsion_dynamic_feed_line_v1`, a propagated two-state
  inertance/compliance liquid-line model with quadratic damping, compliant
  storage-mass closure, reverse-flow/pressure/flow roots, modal diagnostics,
  deterministic model identity and analytic local Jacobians.
- Composed the dynamic line with regulated pressurization, bounded
  pressure-dependent engine performance and the propagator bridge. Engine
  demand enters the manifold-pressure ODE, line flow drains the tank, and
  analytic feed-state/line-state/acceleration derivatives are returned without
  a hidden algebraic loop.
- Added central-finite-difference checks of the fully composed dynamic-line
  Jacobian, inactive-engine relaxation behavior, double-resistance refusal,
  reference/docstring audit coverage, documentation and a runnable example.

## 0.9.0 — engineering preview

- Added variable-property thermochemistry, external CEA adaptation, and an
  experimental local equilibrium provider.
- Added off-design nozzle regimes, internal shocks, altitude performance, and
  separation screening.
- Added Rao and MOC contours, comparison evidence, and CSV/SVG/DXF/STL export.
- Added traceable performance losses, Bartz heat flux, experimental regenerative
  cooling, and a two-phase performance band.
- Added versioned workspaces, seeded studies, staging optimization, engine-cycle
  balances, and printable reports.
- Expanded the browser workbench, API metadata, CI, validation docs, responsive
  keyboard navigation, and print output.

Known release blocker: the nine-point external NASA CEA comparison matrix is
not closed on this machine. See `VALIDATION_REPORT.md`.

