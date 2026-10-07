# Astrodynamics integration guide

This guide is for the developer of an orbit propagator, orbit-determination
filter or trajectory optimizer (Sidera, Orekit-style tools, or your own) who
wants to take propulsion models from Rocket Propulsion Lab. It lists what each
propagation-facing module provides, how to wire it into an integrator, and
what remains the consumer's responsibility.

Everything here is pure Python with no dependencies and no import of any
astrodynamics package. All quantities are SI.

## 1. Which module for which job

| Need | Module | Main type | Detail |
|---|---|---|---|
| Chemical burn with a known thrust/flow history | `propagator.py` | `PropagatorPropulsionBridge` | [Bridge](PROPAGATOR_PROPULSION_BRIDGE.md) |
| Thrust that depends on tank/feed pressure, throttle, power … | `performance_surface.py` | `RectilinearPerformanceSurface` | [Surfaces](OPERATING_CONDITION_SURFACES.md) |
| Ascent or landing: ambient back-pressure | `propagator.py` | `nozzle_exit_area_m2`, `ambient_pressure_pa` | [Bridge](PROPAGATOR_PROPULSION_BRIDGE.md#ambient-pressure) |
| Pressure-fed liquid stage with regulator dynamics | `feed_system.py`, `feed_line.py` | `RegulatedFeedSystem`, `DynamicFeedLine` | [Feed](REGULATED_FEED_DYNAMICS.md), [Line](DYNAMIC_FEED_LINE.md) |
| Ion / Hall thruster limited by available power | `electric.py` | `PowerLimitedPropulsionModel` | [Electric](ELECTRIC_PROPULSION.md) |
| Priors / consider covariance for OD, LinCov | `execution_errors.py` | `PropulsionParameterPrior` | [Errors](EXECUTION_ERRORS.md) |
| Impulsive maneuver execution errors | `execution_errors.py` | `GatesExecutionErrorModel` | [Errors](EXECUTION_ERRORS.md) |
| Slowly wandering thrust on long arcs | `execution_errors.py` | `FirstOrderGaussMarkov` | [Errors](EXECUTION_ERRORS.md) |
| Torque, cosine loss, CoM/inertia drift for 6-DOF | `vehicle_geometry.py` | `ThrusterLayout`, `PropellantMassModel` | [Geometry](VEHICLE_GEOMETRY.md) |
| Bounded (0 ≤ T ≤ T_max) thruster allocation for a force/torque demand | `vehicle_geometry.py` | `ThrusterLayout.allocate` | [Geometry](VEHICLE_GEOMETRY.md#thrust-allocation) |
| Several burns, coasts and stage/fairing jettison on one clock | `timeline.py` | `PropulsionTimeline` | [Timeline](PROPULSION_TIMELINE.md) |
| Throttle as an optimization control, smooth Isp(u) | `throttleable.py` | `ThrottleableEngineModel`, `MonotoneCubicCurve` | [Optimization](TRAJECTORY_OPTIMIZATION.md) |
| Impulsive equivalent (Δv, t*) of a finite burn | `burn_metrics.py` | `free_space_burn_metrics` | [Optimization](TRAJECTORY_OPTIMIZATION.md#3-impulsive-equivalent-of-a-finite-burn) |
| Saltation for sensitivities at burn events | `propagator.py` | `time_event_jump` | [Optimization](TRAJECTORY_OPTIMIZATION.md#4-variational-equations-across-burn-discontinuities) |
| Monte Carlo realizations with probabilities | `ensemble.py`, `scenarios.py` | hierarchical ensemble | [Ensemble](HIERARCHICAL_PROPULSION_ENSEMBLE.md) |

All public types are importable from `rocket_propulsion.propulsion.burns`.

## 2. Ownership boundary

| Rocket Propulsion Lab returns | The astrodynamics framework owns |
|---|---|
| net force magnitude along a supplied direction | epoch and time scales |
| vehicle and per-tank mass derivatives | attitude, guidance and frame transforms |
| event (root) functions and what each one means | event location and integrator restarts |
| analytic partials with respect to its inputs | variational equations and state ordering |
| available-power → thrust law | the power model (solar distance, eclipse) |
| propulsion uncertainty (priors, Gates magnitude terms) | pointing statistics and trajectory dispersion |

No module integrates position or velocity, and none clamps its output inside
an integration step when a correctly registered event would do the job.

## 3. State vector the consumer integrates

A typical augmented state for a chemical stage:

```text
r (3), v (3), m                       translational + vehicle mass
m_ox, m_fuel, ...                     one entry per StreamStateBinding target
[feed states]                         only with RegulatedFeedSystem / DynamicFeedLine
[x_GM]                                optional Gauss-Markov thrust-scale deviation
```

For electric propulsion: `r, v, m, m_xenon`, plus the discrete throttle
`level` kept outside the continuous state.

## 4. One integration step

```python
value = bridge.evaluate(
    relative_time_s=t - t_ignition,
    vehicle_mass_kg=m,
    direction=thrust_direction_in_propagation_frame,
    additional_states={"oxidizer_mass_kg": m_ox, "fuel_mass_kg": m_fuel},
    ambient_pressure_pa=p_atm,          # only if nozzle_exit_area_m2 is set
    apply_inventory_limits=False,       # because inventory roots are registered
    side="right",                       # "left" for a stage that ends on a boundary
)
dv_dt += value.acceleration_m_s2
dm_dt = value.mass_derivative_kg_s
dm_tank_dt = dict(value.additional_state_derivatives)
```

Rules that the tests in `tests/test_burn_propagator_integration.py` and
`tests/test_burn_electric.py` demonstrate with an RK4 loop:

1. **Register every root once at setup**: `event_surfaces()` (start, segment
   boundaries, stop, dry mass, reserves), `operating_condition_surfaces()`
   (stopping) and `operating_condition_grid_surfaces()` (restart only) when a
   performance surface is attached; `net_thrust_margin_n` when ambient pressure
   is used; `switching_surfaces(level)` and `inventory_surfaces()` for electric
   propulsion.
2. **Land on time events.** Use the right-hand limit at the start of a step and
   the left-hand limit at its end.
3. **Locate state roots, then switch the burn off yourself**, and pass
   `apply_inventory_limits=False`. The built-in clamp is only a safety net for
   consumers without event handling; inside a Runge-Kutta step it biases root
   location (measured: 6.7e-4 s with a 0.01 s step).
4. **Electric mode changes:** after a `PowerSwitchSurface` root, set
   `level = surface.to_level`. Never recompute the level from the power at the
   located root.
5. **Check mass consistency** with `mass_closure()`; the residual should stay
   at round-off.

## 5. Partials for variational equations and estimation

| Partial | Bridge | Electric |
|---|---|---|
| ∂a/∂m | `acceleration_wrt_mass` | `acceleration_wrt_mass` |
| ∂a/∂d (3×3, supplied direction) | `acceleration_wrt_direction` | `acceleration_wrt_direction` |
| ∂a/∂p_ambient | `acceleration_wrt_ambient_pressure` | — |
| ∂a/∂P_available | — | `acceleration_wrt_available_power` |
| ∂a, ∂ṁ, ∂ṁ_tank / ∂conditions | `*_wrt_conditions` | — |
| ∂a, ∂ṁ, ∂ṁ_tank / ∂drivers | `*_wrt_parameters` | `*_wrt_parameters` |

Drivers (`thrust_scale`, `mass_flow_scale`, `ignition_time_bias_s`) are named
and bounded, so an estimator cannot reorder them silently. Use
`PropulsionParameterPrior` for their a-priori covariance and
`prior.check_against(bridge.parameters)` before a filter run.

Sensitivities across events: profile time events move with
`ignition_time_bias_s`, so the sensitivity matrix needs the saltation update
`S+ = S- + (f- - f+) dtau/dp` at each one. `bridge.time_event_jump(...)`
returns the propulsion part of `f- - f+`. `tests/test_burn_variational.py`
integrates the full recipe and matches finite differences of whole
propagations to six digits.

Partials are local. Where a quantity is non-smooth (a time boundary for the
ignition bias, a grid line of a performance map, a throttle level change), the
bridge reports it (`nonsmooth_parameters`) or exposes a root instead of
returning a misleading derivative.

## 6. Data provenance

Every exchanged artifact (tabulated burn, performance surface, throttle table,
parameter prior, bridge manifest) is canonical JSON with a SHA-256 content hash
and reference identifiers from `references.py`. Re-opening a modified artifact
fails. Store the `model_hash` with every trajectory result so a propulsion
change can never pass silently into an orbit product.

## 7. Robustness guarantees

- **Fail-closed inputs:** every public boundary of these modules validates its
  inputs through `core/validation.py`. `None`, text, `bool`, containers, NaN
  and infinity raise `DomainError` or `InputError` naming the field, never a
  bare `TypeError`. A seeded fuzz test (`tests/test_burn_input_robustness.py`)
  enforces this.
- **Canonical numbers:** integers are normalized to floats at construction, so
  an artifact built with `1000` hashes exactly like one reopened from JSON.
  Before this, such artifacts could be written but not reopened.
- **Reopenable models:** `PropagatorPropulsionBridge`,
  `PowerLimitedPropulsionModel` and `PropulsionTimeline` serialize to complete
  JSON documents. Reopening re-verifies every nested hash. Ship the document
  with a trajectory product so the exact propulsion model can be reproduced.
- **Code review:** a systematic review of the whole interface found ten
  issues. Each is fixed and has a regression test (`ReviewRegressionTests`),
  for example readers that coerced `"1.0"` to a number, event-id namespace
  collisions in timelines, and a quadrature error estimate that could be a
  fake zero.
- **Property tests:** random multi-segment profiles, throttle tables and
  geometries are checked for continuity, drain closure, linearity, RK4
  propellant and mass closure, monotone throttle paths and bounded-allocation
  optimality against brute-force enumeration.

## 8. Known limits

- Steering, attitude dynamics and the decision when to start or stop a burn
  belong to the consumer. `PropulsionTimeline` sequences burns and jettisons
  but does not target them.
- Thrust direction, steering laws, gimbal dynamics and attitude control belong
  to the consumer; `vehicle_geometry.py` provides the forces, torques and mass
  properties those models need.
- Performance maps interpolate supplied evidence; no thruster or engine data
  are bundled.
- Ambient pressure enters as `p_a A_e` with a fixed exit area; nozzle flow
  separation at high back-pressure is not modelled in the bridge (see the
  off-design nozzle screen in `compressible/off_design_nozzle.py`).
