# Documentation

- [Tabulated relative and epoch clocks](TABULATED_TIME_CONTRACT_REVIEW.md):
  exact ownership, explicit snap windows, canonical event dates and actual clock collisions.

- [Free-space and tabulated numerics](FREE_SPACE_TABULATED_NUMERICS_REVIEW.md):
  exact integral ledgers, mass partitions, directed numerical bounds and output-range classification.

- [Burn record and artifact contracts](BURN_RECORD_CONTRACT_REVIEW.md):
  strict raw types, immutable snapshots, reproducible output and complete explicit AST inventory.

- [Native Lunaris impulse numerics](NATIVE_LUNARIS_IMPULSE_REVIEW.md):
  controlled external patch, independent mass-law series and native source binding.

- [Adapter direction and serialized coordinates](ADAPTER_DIRECTION_NUMERICS_REVIEW.md):
  Sidera unit vectors, Lunaris final-vector projection and native boundary checks.

- [Profile exposure and cutoff numerics](PROFILE_BURN_NUMERICS_REVIEW.md):
  earliest command-time inversion, independent moments and protected ledgers.

- [Constant burn target and cutoff numerics](CONSTANT_BURN_NUMERICS_REVIEW.md):
  stable inversion, protected limits and direct tank-consumption ledgers.

- [Tank inventory closure decisions](INVENTORY_CLOSURE_NUMERICS_REVIEW.md):
  exact residual arithmetic, portable decision ratios and tolerance boundaries.

- [Multilinear maps and channel numerics](MULTILINEAR_NUMERICS_REVIEW.md):
  exact cell geometry, independent gradient ledgers and downstream bridge checks.

- [PCHIP and cubic tensor numerics](PCHIP_NUMERICS_REVIEW.md): scale/sign,
  cancellation and mixed-derivative attacks with independent Bezier references.

- [Throttle flow and inventory numerics](THROTTLE_FLOW_NUMERICS_REVIEW.md):
  independent flow-path/closure checks and recognized numerical-policy integrity.

- [Diagnostic evidence review](DIAGNOSTIC_EVIDENCE_REVIEW.md): reproduced
  interpreter diagnostic failure, unchanged hard deadlines and failing current reports.

- [Interval validation review](INTERVAL_VALIDATION_REVIEW.md): outward error
  enclosures, exact-corner soundness checks and explicit numerical/physical scope.

- [Decision numerics review](DECISION_NUMERICS_REVIEW.md): exact affine ties,
  summed contributions, scenario input preservation and unassessed ranking uncertainty.

- [HTTP capability probe](HTTP_CAPABILITY_PROBE.md): nonblocking optional
  inspection, explicit unknown states and concurrency/failure checks.

- [Storage numerics review](STORAGE_NUMERICS_REVIEW.md): intermediate scaling
  defects, positive-output representation limits and independent inventory roots.

- [Phase energy review](POWER_LEDGER_REVIEW.md): duty-boundary cancellation,
  summed-energy underflow and exact rational phase accounting.

- [Mission mass closure review](MASS_CLOSURE_REVIEW.md): feasibility-boundary
  defects, precision-controlled closure and independent forward root checks.

- [Feed physics review](FEED_PHYSICS_REVIEW.md): regulator/thermal/line defects,
  independent conservation oracles and source-bound verification records.

- [Motor and feed-system assurance review](PROPULSION_ASSURANCE.md): actual defects,
  adversarial verification, independent-data research and unclosed hardware TRL5 gates.

- `SASS_PROJECT_PREPARATION.md` records verified UZM 508E course identity,
  inferred propulsion project branches, field-level references, preliminary
  sizing/trade tools, strict JSON/HTML workflow and the no-extrapolation policy.

- `ARCHITECTURE.md` explains dependency direction and extension points.
- `ROADMAP.md` records the staged engineering scope.
- `DEVELOPMENT_PLAN.md` turns the open gaps into a prioritized, test-gated
  implementation plan from version 0.3 to 1.0.
- `REFERENCES.md` records primary equation and screening-data sources.
- `GAP_ANALYSIS.md` distinguishes covered physics from the next fidelity limits.
- `USER_GUIDE.md` describes the end-to-end browser workflow and fidelity labels.
- `API.md` documents the JSON envelope and Python/PowerShell examples.
- `TOLERANCES.md` centralizes numerical acceptance thresholds.
- `VALIDATION_REPORT.md` records bounded software evidence and open external validation gates.
- `RELEASE_NOTES.md` summarizes packaged engineering-preview changes.
- `HOT_FIRE_CALIBRATION.md` documents scalar and paired test-data calibration.
- `MODEL_DISCREPANCY_TRANSFER.md` separates model-form and
  qualification-to-flight epistemic effects from hot-fire variability.
- `PROPAGATOR_PROPULSION_BRIDGE.md` defines the event-safe force, total/tank
  mass derivative, estimator-driver, and analytic-Jacobian boundary for Sidera.
- `OPERATING_CONDITION_SURFACES.md` documents bounded multidimensional
  supply-pressure/ambient/throttle/power performance maps and their gradients.
- `REGULATED_FEED_DYNAMICS.md` defines the pressurant bottle, dynamic
  regulator, thermal ullage, line-loss ODE and implicit feed/engine closure.
- `ASTRODYNAMICS_INTEGRATION.md` is the entry point for propagator, orbit-
  determination and trajectory-optimization developers: module map, state
  layout, event protocol and partials.
- `EXECUTION_ERRORS.md` defines driver priors (including the hot-fire
  calibration mapping), the Gates execution-error model and Gauss-Markov
  thrust process noise.
- `VEHICLE_GEOMETRY.md` defines thruster force/torque geometry and
  propellant-dependent centre of mass and inertia with analytic derivatives.
- `TRAJECTORY_OPTIMIZATION.md` covers the throttle-commanded engine, smooth
  monotone throttle curves, impulsive equivalents of finite burns and the
  verified variational-equation recipe with saltation at burn events.
- `LUNARIS_INTEGRATION.md` connects timelines to the Lunaris lunar-orbit
  propagator: exact-mass impulsive equivalents, finite burns inside Lunaris,
  bias-shifted segmentation, exact dry-mass cuts and the end-to-end
  verification.
- `ADVERSARIAL_TESTING.md` explains the attack suites, their independent
  oracles, the mutation checks and the defects they found.
- `PULSE_MODE.md` defines pulse-mode thrusters: measured impulse bits, the
  minimum impulse bit, fewest-pulse plans with reported quantization and
  pulse-train artifacts.
- `PROPULSION_TIMELINE.md` defines multi-burn schedules, jettison events,
  concurrent-burn summation and verified timeline documents.
- `ELECTRIC_PROPULSION.md` defines the power-limited electric thruster model:
  hash-verified throttle tables, throttle paths with hysteresis, the
  fixed-efficiency law, power-switching roots and power/direction partials.
- `DYNAMIC_FEED_LINE.md` defines the propagated liquid-flow/manifold-pressure
  states, conservation law, cutoff roots and composed Jacobians used to add a
  calibrated first feed-line mode to an orbit propagator.

Implementation docstrings also describe the equations. That proximity does not
automatically guarantee agreement between the implementation and documentation.

`DIRECT_INVENTORY_HELPERS_REVIEW.md` covers direct constant-flow tank helpers,
their safe binary64 duration horizon, exact reserve decisions and attack scope.

