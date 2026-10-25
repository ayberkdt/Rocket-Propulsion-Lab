# Release notes

Raised three simplified parts to engineering models. The regenerative-cooling
solver accepts channels that vary along the chamber (`ChannelProfile`) and
helical channels, and its coolant pressure drop now carries the RTE momentum
term next to friction (Wadel, NASA/TM-1998-206313, eqs. (14)-(16)); for the
methane example the momentum term is half of the 2.86 MPa drop, and the
earlier 3 mm default channel was found to choke, so the default depth is now
5 mm. Stations report coolant velocity and Mach number, with a warning above
Mach 0.5. A new `studies.global_sensitivity` module adds rank-based measures
(Spearman, PRCC, SRRC, rank R²; Helton and Davis, SAND2001-0417) and Sobol
first-order and total indices with bootstrap intervals (Saltelli et al. 2010),
returned by the Δv uncertainty route; the persisted burn-uncertainty documents
keep their Pearson records. Rocket performance finds the throat with a bounded
Brent search and warm-starts each isentrope, cutting equilibrium solves per
call from 481 to 268 with c* and Isp unchanged to 3e-9; a non-finite mass
flux now fails closed instead of returning a point.

Hardened the existing systems against independent oracles and hostile input.
Eight new oracle suites check compressible flow, conduction, real-fluid
equations of state, Gibbs equilibrium, finite-rate kinetics, radiation,
rocket performance and incompressible aerodynamics against SciPy/mpmath
solutions, Monte Carlo ray tracing, KKT certificates and conservation laws
written in the test files; an AST-driven fuzzer sends type garbage and
extreme values to every field of every HTTP route. Defects they found and
fixed: the vortex-panel influence integrals used an unsigned distance, so
cambered sections were wrong (NACA 4412 zero-lift cl 1.44 instead of 0.517);
isentropic inverses are now closed-form or bracketed Newton (0.017 ms instead
of 2.6 ms, to a few ulp); ideal-gas equilibrium refuses pressures above 1 GPa
instead of hanging; 41 overflow, division and time-out failures under extreme
inputs now fail closed with named codes (`numeric_range`,
`sample_budget_exceeded`); the propagator kept thrust exact when no nozzle
exit area is set; a bare 101 325 Pa literal now uses the registered
constant; the SP-8120 docstring link pointed at SP-8041. Nozzle separation
gained Schmucker's Mach-number criterion (NASA TM-77396), which locates the
incipient separation point next to the Summerfield screen. `mypy` now runs in
CI with `check_untyped_defs` and `warn_unused_ignores`.

Added heat-transfer and combustion depth. A new `heat_transfer` package covers
conduction (multilayer k(T) walls, fins, analytic and energy-audited
finite-volume transients), convection (pipe, entry, external, natural,
boiling and condensation correlations, Taylor and Hess–Kunz hydrogen,
Hatch–Papell film cooling) and radiation (Planck, view factors, gray
enclosures, radiation-cooled extensions, H₂O–CO₂–soot WSGG), each checked
against the worked examples of Lienhard's *A Heat Transfer Textbook* v6.00. A
coupled regenerative-cooling solver marches the Bartz gas side, wall
conduction with fin lands and the real-fluid coolant (new NIST-verified
viscosity and conductivity for H₂, O₂, N₂, CH₄ and He) with an energy residual
below 1e-9 of the load. Thermochemistry gained a finite-area combustor,
product composition with radicals and dissociation energy, combustion
efficiency models and the TP-1999-208522 Isp loss chain (all against NASA CEA
3.3.4), condensed-phase equilibrium (graphite, water, alumina) and finite-rate
H₂–O₂ nozzle kinetics with a RODAS integrator. The Stefan–Boltzmann, Planck,
Boltzmann and light-speed constants are registered. Six HTTP routes and three
browser views expose the new models; `tests/test_cross_validation.py` checks
that independent models agree where their physics overlap. The README was
rewritten around a capability map and the verification evidence.

Added a literature-verified thermodynamics and aerodynamics layer. Real-fluid
reference equations of state for hydrogen, oxygen, nitrogen, methane and
helium match NIST WebBook within 11 ppm; NASA CEA transport properties and
equilibrium derivatives reproduce CEA 3.3.4. New compressible modules cover
thermally perfect gas flow, Taylor–Maccoll cones (NASA SP-3004), shock-expansion
and Busemann airfoils, compressibility corrections and Newtonian hypersonics. A
new `aerodynamics` package adds potential flow, NACA sections, thin airfoil,
lifting panels and lifting line, boundary layers (Blasius, reference
temperature, Van Driest II), Couette flow, Tauber/Sutton–Graves heating,
Allen–Eggers and point-mass entry, and launch-vehicle methods (Barrowman,
Allen–Perkins crossflow, NACA TR 1307 interference, SP-8035 loads, an
SLS-style aerodynamic database). Seven HTTP routes and six browser modules
expose them. Source misprints found while verifying (Anderson, Barrowman,
SP-8035) are recorded beside the tests. `TRL_ASSESSMENT.md` states what TRL 5
still needs from hardware.

Closed the Phase 2 NASA CEA release gate. Nine points (LOX/LH2, LOX/CH4,
LOX/RP-1 × three O/F, 6.895 MPa, area ratio 40) were produced by NASA CEA
3.3.4 and recorded with its version and database hash; the new
infinite-area-combustor solver (`thermochemistry.rocket_performance`) matches
them within 2e-5 against a declared 1 % tolerance, reproduces RP-1311 example 8
at its published station, and agrees at 0.5 and 25 MPa for every pair. The
local equilibrium solver gained 26 gas species and six condensed reactants
from NASA CEA `thermo.inp` (pinned commit, `scripts/sync_nasa_thermo.py`), a
warm-started Newton that no longer falls back to a zero-trace composition for
carbon systems, and NTO/MMH and local LOX/RP-1 pairs; liquid reactants now
enter at their storage-state enthalpy, which changes local sweep values
relative to 0.9.0's 298 K gas reactants. NASA's `cea` package is an optional
provider (`pip install ".[cea]"`). Off-design nozzle regimes can now be driven
by a perfect gas calibrated to the equilibrium design point. The propellant
catalogue links each family to its computable pair or states why it cannot be
computed, with page-cited NASA values. Added cited coolant and wall-material
catalogues, wall conduction with a material-limit check, centrifugal pump
affinity, suction and map relations, HTTP routes and browser modules for the
transient core (feed, blowdown, electric, pulse mode, uncertainty, scenarios,
calibration, timeline, SASS), distributed water hammer, regulator hysteresis,
layered ullage stratification, hot-fire calibration fields over the
operating-condition surface, a native Sidera wrapper for the propagation
bridge, and a clean-environment package smoke test
(`scripts/package_smoke_test.py`). See `VALIDATION_REPORT.md`.

Added transpiration cooling (`heat_transfer.transpiration`: Mickley film
theory, the wall energy balance with efficiency and radiation,
Eckert–Livingood and Rannie relations, the two-temperature porous wall), the
first hot-fire comparison (Bartz against NASA TN D-6396: 7–30 % high in the
chamber, 13–61 % at the throat), the Priem–Heidmann vaporization length
(`thermochemistry.vaporization`) and NASA SP-194 combustion-instability
relations (`thermochemistry.combustion_instability`). `hot_gas_from_reactants`
builds the chamber gas for any reactant mixture. Also added ethanol and MMH
flow-boiling burnout correlations, feed-coupled chugging with SP-194 feed
admittances and a Nyquist count, a total heat-load comparison with the NASA
TN D-3532 firings, and the `/api/v1/transpiration-cooling`,
`/api/v1/combustion-stability` and `/api/v1/propellant-vaporization` routes
with workbench panes.

Repaired relative-time boundary ownership, early bridge cutoff, shifted-time
rounding and native event merging. Exact origins/endpoints and zero implicit
snap windows preserve arbitrarily short representable relative profiles.
Epoch events are rounded once from complete sums; actual collisions fail
explicitly while distinct one-ULP events remain valid. Added 26 independent
adversarial methods and 16 selected faults, including time-translation invariance
and a shared epoch saltation path. Artifact/bridge/timeline documents
require a new hashed time policy; regenerate older documents. Full raw-document
and native reserve/consumption arithmetic reviews remain separate.
See `TABULATED_TIME_CONTRACT_REVIEW.md`.

Repaired free-space Δv/centroid erasure and tabulated integral/moment range loss.
Exact ledgers and mass partitions precede normalized quadrature; independent
directed Bernstein enclosures include final output rounding, with exact identity
and boundary-refinement paths. Added centroid/exhaust numerical bound metadata,
38 adversarial methods and 21 selected faults. Tabulated artifacts and bridge
documents require new hashed numerical-policy fields; regenerate older documents
from their source samples. See `FREE_SPACE_TABULATED_NUMERICS_REVIEW.md`.

Burn v1/v2 records now normalize valid physical scalars and recognized enums,
snapshot nested collections and reject malformed flags, counts, inventories
and UTF-8 identities. Artifact readers reject unknown/duplicate fields and
type-changing coercions; writers also verify solver reproduction. Range-safe
closure/budget comparisons preserve the existing declared thresholds. Added
25 record/artifact methods, two explicit AST inventory methods and 15 selected
faults. See `BURN_RECORD_CONTRACT_REVIEW.md` for compatibility and scope.

Direct constant-flow inventory helpers now reject missing/ambiguous tank
coverage, malformed time inputs and actual reserve crossings. Exact grouped
stream ledgers preserve small paths before decisions; the duration helper
returns the greatest safe binary64 horizon with an unrounded binding decision.
Added 16 adversarial methods and eight selected faults. See
`DIRECT_INVENTORY_HELPERS_REVIEW.md` for sentinel and output-range semantics.

## Unreleased — motor/feed assurance review

- Applied a byte-bound native Lunaris impulse correction: norm and mass-law
  factors stay in Decimal96 through final export, recoverable exponential range
  loss is preserved, and true norm/mass/velocity output range errors fail clearly.
  Added 15 native adversarial methods, six native fault injections and a stale
  dependency-evidence attack. Candidate/live application status and scope are in
  `NATIVE_LUNARIS_IMPULSE_REVIEW.md` and its patch manifest.

- Repaired raw-norm overflow/underflow in Sidera constant/ensemble directions and
  Lunaris impulse projections. Shared checked final-vector arithmetic retains
  signed components and recoverable small delta-v components. Ensemble reopening
  validates and preserves serialized coordinates instead of normalizing them
  again; its hashed numerical-policy field changes the handoff identity. Added
  14 adversarial methods and five selected faults, including pinned native
  constant-burn/impulse boundary checks. See `ADAPTER_DIRECTION_NUMERICS_REVIEW.md`.

- Repaired L1 positive-target zero-output success, dimensional slope clipping,
  near-boundary cutoff overrides and erased consumption. Exact command integrals
  and exposure decisions precede Decimal1100 nonlinear roots and checked exports.
  Impulse closure includes rounded output-time error; exact constant export now
  requires genuinely equal command values. Added 25 adversarial methods and nine
  selected faults. See `PROFILE_BURN_NUMERICS_REVIEW.md` for scope and limits.

- Repaired L0 small-delta-v inversion/forward errors, target-tolerance cutoff
  overrides and erased tank consumption. Rational ledgers and Decimal1100
  transcendentals precede checked exports; zero inventory is a finite cutoff.
  Added 16 adversarial methods and six selected faults with independent series
  references. Method/result identities change; the L1 profile path remains a
  separate review. See `CONSTANT_BURN_NUMERICS_REVIEW.md`.

- Repaired bridge/electric inventory residual cancellation and rounded-residual
  false closure decisions across all three models. Shared exact ledgers retain
  a JSON-portable residual ratio; public records validate finite fields and
  relative consistency. Added 13 adversarial methods and five selected faults;
  two existing inventory faults now target the shared helper. Model policies
  and identities change. See `INVENTORY_CLOSURE_NUMERICS_REVIEW.md`.

- Repaired width/fraction range loss, erased corner contributions and
  non-finite/cancelled gradients in both multilinear surface classes. Shared
  rational weights and per-channel ledgers precede checked output conversion.
  Added 14 adversarial methods and six selected faults, with downstream bridge
  conservation checks; changed surface identities bind numerical policies.
  See `MULTILINEAR_NUMERICS_REVIEW.md`.

- Repaired small-secant flattening, intermediate range loss and cancellation
  in the shared PCHIP kernel, throttle curves and cubic tensor performance
  surfaces. Rational nodal/mixed derivatives and polynomial sums precede
  checked field exports. Added 19 adversarial methods and seven selected faults;
  numerical policy fields bind changed identities. See `PCHIP_NUMERICS_REVIEW.md`.
  Exact power coefficients are compiled once for repeated curve queries;
  pulse searches skip unused derivatives and retain final PulseBit slopes.

- Repaired throttle flow/tank ledgers, command/driver derivatives, subnormal
  Isp errors and inventory residual cancellation. Accepted stream fractions
  normalize by their exact binary sum; final outputs fail explicitly when
  genuinely unrepresentable. Added 14 methods and seven curated faults.
  Bridge/throttle/electric readers now verify known numerical policy fields;
  see `THROTTLE_FLOW_NUMERICS_REVIEW.md` for verification limits.

- Repaired scale-dependent direction normalization in the bridge, throttleable
  and electric models. Scaled norm estimation and checked rational projections
  preserve representable force and signed direction tangents; near-axis
  diagonal cancellation no longer erases sensitivity. Added nine adversarial
  methods and five curated faults; the integrated transverse-speed oracle now
  uses a cross product with its original tolerance. See `DIRECTION_NUMERICS_REVIEW.md`.

- Repaired upstream bridge flow/tangent erasure, order-dependent grouped
  cancellation and recoverable intermediate overflow. Mass profiles, stream
  scaling and tank/total ledgers now remain exact through grouping and export
  once with explicit range failures. Added seven public-bridge oracle methods
  and five curated faults; force/direction arithmetic has a separate scoped review.

- Repaired lost/undefined mixture-ratio pressure derivatives at small common
  flow scales and cancellation in correlated tangents. Exact arithmetic on
  supplied binary values precedes one finite output conversion. Composed
  sensitivity records now reject non-finite/malformed values and copy mutable
  vectors. Added analytical scaled-map and independent Decimal path oracles
  with four curated fault injections; physical uncertainty remains unassessed.
- Overall software-check status now requires the complete declared mutation
  registries as well as pytest success. Survivors, invalid experiments and
  inconsistent/omitted outcomes cannot produce a successful overall report;
  the pytest result remains separately visible. Archived the prior green run
  with repository input bytes and JUnit/log evidence before this review.
- Test-collection timeouts/errors now retain partial output and replace old
  successful manifests with current failing state. Collection/module deadlines
  remain 30/120 seconds with four workers; startup latency variability and
  startup latency variability remains open. The later bridge-flow review
  addresses the previously observed component-gradient loss.

- Repaired signed tank/line conservation through zero-flow event-location
  trials in single and bipropellant couplings. Negative trial results are
  explicitly marked; physical propagation must stop at zero. The vehicle
  closure now uses the actual returned tank derivatives. Added independent
  mass/energy ledgers, an integrated cutoff-root check and three targeted faults.
- Hardened bipropellant configuration, state and event boundaries. Boolean
  physical limits and malformed records now fail with typed domain errors;
  known string event kinds retain the correct minimum/maximum semantics and
  unknown kinds are rejected. Added four hostile methods and two targeted faults.
- Reduced redundant pulse width-search interpolation with mass-only sampling,
  a bounded per-model mass inverse cache and per-call pulse reuse. The random
  enumeration workload, tolerances and regression deadline are retained.

- Reproduced a native diagnostic crash in a propulsion-free control process;
  replaced the C watchdog snapshot with a GIL-held Python callback while
  retaining fatal handlers, test cases and the 120-second process deadline.
- Incomplete/crashed evidence can no longer leave an earlier green summary;
  current incomplete/failed reports exit nonzero and module exit codes are recorded.

- Repaired inward-rounded residual bounds in measurement/model comparisons.
  Exact intervals now publish outward endpoints and a conservative numerical
  review guard. Independent corner checks exercise 256 fixed-seed cases plus
  extreme output limits; numerical agreement still makes no physical/TRL claim.

- Repaired spurious ordering of equal affine decision scores and lost summed
  contributions using exact weighted arithmetic. Positive scenario scaling
  can no longer silently erase a mission allocation. Ranking robustness under
  physical input/model uncertainty remains explicitly unassessed.

- Repaired the HTTP capability timeout at its cause: optional Sidera imports
  now run once in a background inspection. API metadata v1.1 exposes checking
  with unknown flags omitted, then a cached terminal record. Added actual HTTP
  gated-import/concurrency/failure tests and restored the fully parallel run.

- Repaired cold-gas and regulated-pressurant intermediate overflow/underflow
  and silent positive ullage loss. Storage and power ledgers share exact
  rational arithmetic/final conversion; independent endpoint-inventory roots
  and two additional curated semantic mutations exercise the repaired limits.

- Repaired doubled/lost standby energy near continuous duty and lost summed
  subnormal phase energy using exact rational evaluation before final output
  rounding. Both reproduced defects are included in curated mutations.

- Repaired residual-boundary propellant errors, overflowing intermediate mass
  ratios and silent positive-burn underflow. Mass sizing now records isolated
  Decimal precision agreement, conditioning and separate binary/working
  closure residuals. Added independent forward-root boundary checks.
- Curated mutation reports now reject stale source/test/tool snapshots.

- Hardened regulator flow near equilibrium and gamma=1, gas heat capacities,
  thermal energy evaluation, dynamic-line frequencies, finite ODE/Jacobian
  outputs and rounded-zero ullage handling. Added independent Decimal and
  conservation oracles with eight additional curated physics mutations.
- Verification summaries now derive results from collected test identities,
  JUnit and current source/data/dependency digests; stale or omitted evidence
  cannot silently preserve an earlier passing summary.

- Reproduced and repaired nine numerical/input integrity defect classes in
  the SASS addition. Added a 75-digit independent mass/orbit oracle and eight
  curated semantic mutation experiments with recorded assertion failures.
- Added reviewed holdout/uncertainty comparisons, source byte digests,
  staged artifact writes and a per-symbol audit inventory. Experimental
  measurements remain distinct from simulations and corrected estimates.
- Reassessed five catalogue entries; none provides validated joint operating
  data for qualification. JSON/HTML now carry explicit hardware/uncertainty
  limitations. Hardware design is unselected and TRL5 evidence gates are open.
- See `PROPULSION_ASSURANCE.md` and generated `output/assurance/` records for
  scope, unresolved checks and independent-data research.

## Unreleased — cited UZM 508E propulsion project preparation

- Researched the official UZM 508E course record and separated verified
  course identity from inferred project scope. The public course description
  does not provide an assignment brief.
- Added evidence-bearing quantities, field-level thruster specifications,
  ordered mission/residual propellant closure, impulsive/low-thrust analytic
  budgets, storage/pressurant inventories and phase-based bus energy budgets.
- Added explicit requirement screens, fixed-anchor weighted trades and
  deterministic mission-input sensitivity studies; missing data stays unknown.
- Added a strict citation-bearing JSON request, printable offline HTML/JSON
  report, synthetic example and source/function audit tests. Catalogue ranges
  and maxima cannot create nominal points, and extrapolation is disabled.

## Unreleased — standalone burn vertical slice

- Consolidated the advanced propulsion work on one Sidera-first branch while
  retaining Lunaris only as a tested legacy compatibility adapter.
- Added the hash-bound `rocket_propulsion_sidera_provider_handoff_v1` contract:
  absolute/relative time mapping, explicit named-state order, event inventory,
  operating inputs, estimable drivers, nested bridge integrity, and optional
  body-frame force/moment closure.
- Added non-extrapolating N-D engine maps with independent named propellant-
  stream scales and analytic force, total-flow, and component-flow gradients.
- Added dual regulated bipropellant feed/line coupling with engine and line O/F,
  first-component depletion, positive-safe O/F events, mass/storage closure,
  and a composed four-state Jacobian.
- Added ordered two-axis gimballed engine clusters with body-frame force and
  moment about a moving centre of mass, vector-cancellation loss, analytic
  thrust/gimbal partials, and integrity-checked persistence.
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
- Added implicit coupled Jacobians through the pressureâ†”mass-flow algebraic
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
  (0 â‰¤ T â‰¤ T_max) by an active-set Lawson–Hanson / Stark–Parker method with
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

