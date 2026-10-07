# Standalone burn analysis

This package owns propulsion-side finite-burn bookkeeping. It intentionally
does not propagate position, velocity, attitude, or orbital frames.

The first implementation level (`L0`) evaluates a constant delivered operating
point analytically. Its explicit stream ledger distinguishes chamber flow from
total tank drain, so gas-generator dump, coolant dump, and vent flows cannot be
hidden inside an ambiguous mass-flow value. System specific impulse is always
closed as `F / (g0 * total tank flow)`.

`ConstantPerformanceProvider.bipropellant()` and `.monopropellant()` construct
closed, immutable operating points from delivered thrust, system Isp and the
declared stream split. `simulate_burn(definition, provider)` is the provider-
oriented entry point. It dispatches the L0 constant provider or an L1
`ProfiledPerformanceProvider` without coupling either solver to an engine family.

L1 uses an ordered `ThrottleSchedule` of exact piecewise-linear
`ThrottleSegment` records. Start ramps, steady operation, throttle steps and
shutdown ramps are split at every boundary. Thrust and every tank stream scale
proportionally with throttle while system Isp remains constant. Exposure
`integral(throttle dt)`, impulse, tank drain, force-time centroid and ideal
one-dimensional delta-v are integrated analytically; presentation samples do
not drive the authoritative result.

`simulate_constant_burn()` supports duration, total propellant, delivered
impulse, and ideal one-dimensional delta-v targets. Tank reserves, protected
dry mass, and maximum duration act as independent protective cutoffs. An
unreachable target returns a valid result with `target_achieved = false` and
the binding constraint; it is not presented as a successful shortened burn.

L0 results use `rocket_propulsion_burn_v1`; L1 results use
`rocket_propulsion_burn_v2` and persist the complete schedule. Both have
deterministic input/result hashes and are physically reproduced when reopened.
They remain independent of Sidera.

## Isolated transient and cluster layer

The next Sidera-facing capabilities are developed here first, without importing
or modifying Sidera:

- `response.py` solves a first-order throttle response to every linear command
  segment analytically. Ignition delay, exponential rise/decay, exact exposure,
  and force-time centroid are retained.
- `sequence.py` validates per-engine minimum pulse width, cooldown, overlap,
  and maximum start count before a profile can be built.
- `cluster.py` composes independently commanded engine members, cant losses,
  shared tank drains, staggered ignition, and engine-out histories. Exponential
  response is converted to a piecewise-linear exchange history only after an
  impulse, mass, axial-impulse, and centroid convergence check.
- `tabulated.py` defines `rocket_propulsion_tabulated_burn_v1`. Its independent
  segment endpoints preserve instantaneous jumps without duplicate timestamps.
  It carries delivered thrust, axial thrust, total tank flow, optional per-tank
  flows, source identity, SHA-256 provenance, interpolation policy, integrals,
  warnings, and a deterministic artifact hash.
- imported CSV thrust/flow curves are replayed with explicit `linear` or `hold`
  interpolation. Samples are never silently smoothed.
- `blowdown.py` adds a preliminary pressure-fed L2 model: rigid-tank
  polytropic ullage expansion, exact pressure/reserve cutoff, variable thrust,
  Isp and mass flow, and convergence-qualified exchange tabulation. Performance
  may come from a bounded power correlation or a non-extrapolating measured
  pressure map. See `docs/BLOWDOWN_MODEL.md` for equations and limits.
- `references.py` provides stable primary-source identifiers. Produced
  artifacts include those identifiers in their content hash, and automated
  documentation tests reject public transient APIs without a `References`
  section.
- `uncertainty.py` generates seeded propulsion ensembles with sourced bounded
  marginals, explicit aleatory/epistemic classification, optional Gaussian-
  copula correlation, empirical percentiles, sampling-stability evidence and
  ranked screening sensitivities. Every draw reruns the deterministic physics
  and retains its complete tabulated artifact for later Sidera propagation.
  See `docs/PROPULSION_UNCERTAINTY.md`.
- `scenarios.py` represents sourced mutually exclusive ignition failure,
  engine-out/common-cause, delay, early-cutoff and performance-degradation
  outcomes. It preserves each outcome as a probability-labelled artifact and
  reports weighted propulsion consequence percentiles without inventing
  component independence. See `docs/PROPULSION_SCENARIOS.md`.
- `ensemble.py` composes every discrete outcome with its own sourced continuous
  uncertainty definition. Joint realization weights are exactly
  `scenario_probability / conditional_sample_count`; complete artifacts,
  scenario contributions and an independently seeded replicate diagnostic are
  retained. See `docs/HIERARCHICAL_PROPULSION_ENSEMBLE.md`.
- `ensemble_serialization.py` defines the integrity-checked
  `rocket_propulsion_hierarchical_ensemble_v1` exchange document. Reopening
  verifies nested burn hashes and independently reproduces weights, statistics,
  total-variance decomposition, scenario contributions and convergence flags.
- `calibration.py` fits a sourced heteroscedastic Mandel-Paule random-effects
  model to repeated hot-fire scale measurements. It separates between-engine,
  within-engine and measurement variance, retains solver evidence and produces
  population or same-engine residual parameters for the conditional ensemble.
  See `docs/HOT_FIRE_CALIBRATION.md`.
- `multivariate_calibration.py` preserves paired thrust/Isp/time evidence,
  subtracts declared measurement covariance, records any positive-definite
  correlation regularization, emits population/same-engine Gaussian-copula
  models, and persists a semantically revalidated
  `rocket_propulsion_hot_fire_multivariate_calibration_v1` document.
- `discrepancy.py` calibrates independent model-form and
  qualification-to-flight ratios as epistemic predictive factors. Its layered
  plan preserves aleatory/epistemic names, refuses evidence reuse, applies all
  consensus corrections once, and feeds the existing hierarchical ensemble.
  See `docs/MODEL_DISCREPANCY_TRANSFER.md`.
- `propagator.py` turns a verified tabulated artifact into simultaneous axial
  acceleration, vehicle/tank mass derivatives, event roots and analytic
  estimator Jacobians. It keeps attitude, frames, epochs and orbit integration
  outside the propulsion core. See `docs/PROPAGATOR_PROPULSION_BRIDGE.md`.
- `performance_surface.py` supplies complete, non-extrapolating N-dimensional
  operating maps with analytic force/flow gradients and validity event roots.
  The propagator bridge composes these scales without losing tank-flow closure.
  See `docs/OPERATING_CONDITION_SURFACES.md`.
- `feed_system.py` propagates liquid inventory, bottle/ullage pressurant mass
  and temperature, plus regulator opening. It closes inlet pressure against
  pressure-dependent engine demand and returns implicit coupled Jacobians.
  See `docs/REGULATED_FEED_DYNAMICS.md`.
- `feed_line.py` replaces the static line-loss closure with propagated liquid
  flow and engine-manifold pressure states. It supplies inertance/compliance
  dynamics, nonlinear damping, positive-safe pressure/flow roots, explicit
  storage-mass closure, and engine-composed analytic Jacobians. See
  `docs/DYNAMIC_FEED_LINE.md`.

The exchange artifact intentionally has no epoch, direction, frame, position,
velocity, central body, or orbital result. A future Sidera adapter will combine
it with maneuver timing and guidance while Sidera performs force propagation.
See `docs/TABULATED_BURN_CONTRACT.md` for the integration contract and refusal
rules.

Current Sidera accepts constant finite burns. A genuinely constant L0/L1 result
maps exactly. A varying L1 profile fails closed unless the caller explicitly
selects `equivalent_constant`; that reduction preserves duration, total impulse
and total propellant, reports its variation/centroid metrics, and does not claim
to preserve trajectory response.

The HTTP surfaces are `POST /api/v1/engine-operating-point`,
`POST /api/v1/burn/simulate`, and `POST /api/v1/burn/report`. The browser Burn
workbench uses those routes, stores the complete canonical artifact inside a
versioned `.rplab.json` case, and offers exact constant-burn Sidera hand-off.
The same request object can be saved to disk and run without a server:

```powershell
rocket-propulsion burn simulate examples/constant_burn_input.json `
  --output burn-result.json --report burn-report.html
rocket-propulsion burn validate burn-result.json
rocket-propulsion burn export burn-result.json --target sidera `
  --output sidera-plan.json --start 100 --direction 1 0 0 --frame inertial
rocket-propulsion burn simulate examples/profiled_burn_input.json `
  --output profile-result.json --report profile-report.html
rocket-propulsion burn export profile-result.json --target sidera `
  --mode equivalent_constant --output sidera-equivalent.json `
  --start 100 --direction 1 0 0
rocket-propulsion burn tabulate-trace examples/solid_thrust_trace.csv `
  --source-id static-fire-001 --output solid-profile.json
```

The model is preliminary engineering software, not flight or safety-critical
command generation. Gravity, drag, steering, and finite-orbit effects require
an external trajectory tool such as Sidera.
