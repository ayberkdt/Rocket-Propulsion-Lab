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

### Propagator bridge hardening

- Segment lookup in `TabulatedBurnArtifact` and `PropagatorPropulsionBridge`
  is now a bisection (`segment_index_at`, `is_boundary_time`) instead of a
  linear scan, and the performance surface builds its node index once instead
  of on every call. Measured on one machine: 2000-segment profile 427 -> 44
  us per evaluation; with a 20^3 surface 1104 -> 91 us. Left-side queries
  within boundary tolerance no longer raise when the previous segment ends a
  round-off earlier.
- Added `PerformanceGridLineSurface`, `grid_line_surfaces()` and
  `operating_condition_grid_surfaces()`: non-stopping roots at every interior
  grid coordinate, where the multilinear gradient jumps.
- Added `PropulsionPartials.acceleration_wrt_direction`, the 3x3 Jacobian with
  respect to the supplied direction vector.
- Added `PropagatorPropulsionBridge.mass_closure()` and `PropellantMassClosure`.
- Added `evaluate(..., apply_inventory_limits=False)` (also forwarded by both
  feed couplings) for consumers that locate dry-mass and reserve roots
  themselves; the default internal clamp biases RK root location.
- Added `tests/test_burn_propagator_integration.py`: an event-aware RK4 loop
  that checks the rocket equation, mass and tank closure, a closed-form steady
  position, the analytic reserve cut time, ignition-bias invariance, the
  direction Jacobian, the grid-line events and lookup equivalence with the
  previous linear scan.

### Power-limited electric propulsion

- Added `ElectricThrottleTable` (hash-verified, provenance-bound throttle
  levels with Isp and total-efficiency consistency checks),
  `DiscreteThrottlePath` (maximum-thrust or maximum-Isp path, up-switch
  hysteresis) and `FixedEfficiencyThrottle` (`T = 2 eta P / (g0 Isp)` between
  start and saturation power).
- Added `PowerLimitedPropulsionModel` with propagator-owned mode state,
  `PowerSwitchSurface` roots carrying a deterministic target level, saturation
  kink roots, power margin, duty cycle, inventory roots, mass closure and
  analytic partials with respect to mass, direction, available power and the
  thrust/flow drivers.
- Added `tests/test_burn_electric.py`, including event-driven RK4 runs under
  decaying power that reproduce analytic switch times, the piecewise rocket
  equation and the energy-to-propellant identity; and
  `examples/electric_propulsion.py`.

### Navigation uncertainty, ambient pressure and vehicle geometry

- `PropagatorPropulsionBridge` accepts `nozzle_exit_area_m2`; evaluations then
  take `ambient_pressure_pa` and return `F = F_vac - p_a A_e`,
  `acceleration_wrt_ambient_pressure`, `pressure_thrust_loss_n` and the
  `net_thrust_margin_n` root. Double counting with a surface axis named
  `ambient_pressure_pa` is refused; vacuum-only manifests are unchanged. Both
  feed couplings forward the pressure.
- Added `execution_errors.py`: `PropulsionParameterPrior` (hash-verified,
  PSD-checked, from sigmas or linearized from a paired hot-fire calibration),
  `GatesExecutionErrorModel` (covariance, sampling, proportional magnitude
  from a prior via `delta_v_log_sensitivities`) and `FirstOrderGaussMarkov`.
- Added `vehicle_geometry.py`: `ThrusterLayout` (resultant force/torque,
  6xN influence matrix, `d tau/d c`, cosine efficiency) and
  `PropellantMassModel` (centre of mass and inertia versus tank masses with
  analytic derivatives).
- Registered NTRS 20180003657 (Navigation Filter Best Practices).
- New documentation: astrodynamics integration guide, execution errors,
  vehicle geometry; rewritten top-level README.
- New tests verify the calibration mapping against a nonlinear Monte Carlo,
  Gates sampling against its covariance, Gauss-Markov autocorrelation, all
  new partials against finite differences, and inertia against an independent
  parallel-axis sum.

### Robustness, persistence, timelines and thrust allocation

- Added `core/validation.py` and applied it at every public boundary of the
  propagation interface: malformed values now raise `DomainError` /
  `InputError` naming the field. A seeded fuzzer found 257 bare
  `TypeError`/`ValueError` crashes before the change and none after.
- Fixed: artifacts built with integer values (tabulated burns, performance
  surfaces, throttle tables, priors) failed their own hash check on reopening.
  Numbers are now normalized to float at construction.
- Fixed: a non-string parameter-override key crashed while formatting its
  error message.
- `PropagatorPropulsionBridge`, `PowerLimitedPropulsionModel`: complete
  `to_dict`/`to_json`/`from_dict`/`from_json` documents that re-verify every
  nested hash and the model hash.
- Added `timeline.py`: `PropulsionTimeline`, `ScheduledBurn`, `MassJettison`
  with overlap checks, namespaced events, jettison mass handling, per-burn
  directions/conditions/drivers, concurrent-burn summation and a verified
  timeline document.
- Added `ThrusterLayout.allocate`: bounded-variable least squares
  (0 ≤ T ≤ T_max) by an active-set Lawson–Hanson / Stark–Parker method with
  Householder QR, reporting residuals and binding limits.
- New tests: fuzzing, integer round trips, bool rejection, property tests on
  random profiles/tables/geometries, persistence tamper tests, a two-burn
  RK4 timeline against the rocket equation, and bounded allocation against
  brute-force optimality.

### Optimization inputs, verified sensitivities and review fixes

- Added `ThrottleableEngineModel` (throttle as a control, deep-throttle gap
  refused, stream split, ambient pressure, partials with respect to throttle,
  mass, direction, pressure and drivers, verified document) and
  `MonotoneCubicCurve` (PCHIP: C1, no overshoot, no extrapolation).
- Added `free_space_burn_metrics` / `gauss_legendre`: free-space delta-v,
  effective exhaust velocity and acceleration-weighted impulsive-equivalent
  time with a quadrature error estimate.
- Added `PropagatorPropulsionBridge.time_event_jump` for the saltation update.
  `tests/test_burn_variational.py` matches integrated sensitivities to finite
  differences of whole propagations for all bridge drivers and for the
  electric power chain rule.
- Code-review fixes, each with a regression test: non-string identifiers now
  fail closed; nested hash failures keep their diagnosis when a bridge is
  reopened; artifact, surface, throttle-table and prior readers no longer
  coerce text to numbers; timeline event ids cannot collide; the timeline
  validates the timing driver once even for inactive burns; the quadrature
  error estimate can no longer be a fake zero; `net_thrust_margin_n`
  semantics are documented.
- Consolidated duplicated validation into `core.validation` (`records`,
  `strings`, `name_value_pairs`, `real_matrix`) and the content hash into
  `serialization.canonical_sha256`; hash values are unchanged.

### Lunaris adapter, C1 surfaces, balanced allocation, wider validation

- Added `integrations/lunaris.py`: `lunaris_impulsive_equivalents`,
  `to_lunaris_maneuver_plan`, `lunaris_segmentation_plan`/`_times`,
  `LunarisFiniteBurnDynamics`, `velocity_direction`. Verified end-to-end with
  the real Lunaris propagator (test skipped when Lunaris is absent).
- `RectilinearPerformanceSurface(interpolation="cubic_hermite")`: C1
  tensor-product cubic Hermite with PCHIP nodal derivatives and analytic
  gradients; multilinear hashes unchanged.
- `ThrusterLayout.allocate(balance_weight=...)`: exact minimum-effort
  allocation; KKT optimality certified in tests.
- `core.validation.coerce_fields`; regulated feed system, dynamic feed line
  and Sidera adapters now fail closed (a fuzzer found 273 bare errors in the
  feed models; none remain).
- Shared `interpolation.pchip_slopes`.

### Lunaris finite burns: ignition bias and exact dry-mass cuts

- Fixed: `LunarisFiniteBurnDynamics` and `lunaris_segmentation_plan` used
  nominal burn times even when `ignition_time_bias_s` was overridden, so the
  biased discontinuities fell inside integration steps. Segmentation now uses
  the resolved bias; the wrapper exposes `segmentation_times` and
  `segmentation_plan()` so the plan always matches the dynamics.
- Added `lunaris_burn_windows` / `LunarisBurnWindow`: exact dry-mass cut times
  for mass-only burns become segment boundaries (Lunaris events cannot see
  mass). The wrapper checks the state mass against the timeline prediction at
  every segment start and fails closed on a mismatch.
- Impulsive equivalents now refuse burns their protected dry mass would
  truncate.

### Pulse-mode thrusters

- Added `PulseModeThruster` (`pulse_mode.py`): measured impulse-bit and
  propellant-per-pulse tables (PCHIP, exact derivatives), minimum impulse bit,
  `plan_impulse` / `plan_delta_v` with the fewest equal pulses and an explicit
  quantization residual (`nearest` or `no_overshoot`), and
  `pulse_train_artifact` for the propagator bridge, timelines and Lunaris.
  Plans are checked against brute-force enumeration and the exact
  rocket equation.
- Registered NTRS 20070032005 (*The Minimum Impulse Thruster*).

### Adversarial test suites, traceable constants and the defects they exposed

- Added `test_adversarial_epochs.py`, `test_adversarial_lunaris.py`,
  `test_adversarial_pulse_mode.py`, `test_adversarial_json_readers.py` and
  `test_constants_and_references.py`, with independent oracles and mutation
  checks; see `docs/ADVERSARIAL_TESTING.md`.
- Fixed: one-sided limits lost at large absolute epochs
  (`epoch_time_tolerance`, `TabulatedBurnArtifact.snap_to_boundary`, bridge
  `boundary_tolerance_s`); Lunaris boundary tolerance now matches the bridge.
- Fixed: Lunaris dry-mass cut uses the bridge clamp margin
  (`INVENTORY_LIMIT_TOLERANCE`) and treats a cut at the stop as the stop;
  bias-overlapping and reordered burns are validated identically on the
  impulsive and finite paths; `segmentation_plan(start_s, end_s)` for partial
  propagation windows.
- Fixed: cancellation in ln(m0/mf) for small burns (`log1p`) in pulse plans,
  burn metrics, both constant/profiled integrators and execution-error
  sensitivities.
- Pulse mode: `plan_delta_v` requires `propellant_budget_kg` and searches
  budget-limited shorter-pulse trains; the pulse limit applies to the returned
  plan; forged plans are refused; a 1e-9 relative tolerance contract.
- All 14 document readers go through `core.validation.read_json_document`:
  non-text input, oversized integers, NaN/Infinity, deep nesting and broken
  documents fail closed with `InputError`.
- `core.constants` registers every physical constant with its source; the
  U.S. Standard Atmosphere 1976 now uses its own R* = 8.31432 J/(mol K)
  (removes a ~3e-5 pressure jump at layer bases) and documents geopotential
  altitude.

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

