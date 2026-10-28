# JSON API guide

All calculations use `POST /api/v1/...` with a JSON object. Success returns:

```json
{
  "data": {"mach": 2.0},
  "meta": {
    "model": "perfect-gas isentropic relations",
    "version": "1.0",
    "assumptions": ["calorically perfect gas"],
    "units": {"mach": "1"}
  }
}
```

Failures use a stable `error` object. Input/domain errors are HTTP 400,
optional feature errors HTTP 501, and convergence failures HTTP 422.

## Examples

```powershell
$body = @{ input_kind = "mach"; value = 2; gamma = 1.4 } | ConvertTo-Json
Invoke-RestMethod http://127.0.0.1:8765/api/v1/isentropic `
  -Method Post -ContentType application/json -Body $body
```

```python
from rocket_propulsion.api.routes import dispatch_calculation

result = dispatch_calculation(
    "/api/v1/nozzle-off-design",
    {
        "chamber_pressure_pa": 5_000_000,
        "chamber_temperature_k": 3500,
        "throat_area_m2": 0.01,
        "area_ratio": 25,
        "ambient_pressure_pa": 101_325,
    },
)
print(result["regime"], result["thrust_n"])
```

The authoritative endpoint inventory and per-route metadata live in
`rocket_propulsion.api.routes.ROUTES` and `ROUTE_METADATA`. Nested arrays are
accepted directly in JSON; form-friendly staging and uncertainty aliases are
also supported by the browser interface.

## Ideal rocket equation and steady thrust

`/api/v1/rocket-equation` and `/api/v1/thrust` use model metadata version 2.0.
Rocket results add `solve_mode`, `performance_input` and `numerical_policy`;
thrust results add `burn_time_s` and `numerical_policy`. These retain the source
quantities needed to reproduce each output. In inverse solves rounded final and
initial masses can coincide while positive propellant consumption remains.
Preserve the original mode/input tags when reconstructing a record. An impulse
record also requires its duration. Named DomainError denotes an invalid input or
unrepresentable mandatory output; ConvergenceError denotes unresolved rounding,
not proof of unrepresentability. See `PERFORMANCE_ARITHMETIC_REVIEW.md`.

## Engine-cycle component sizing

`POST /api/v1/engine-cycle` has model metadata version 2.0. Its default
`pump_feed_policy=shared-pump` includes generator flow in GG pump duty and solves
the feedback balance jointly. `external-gg-feed` applies only to GG and explicitly
excludes its separate supply energy. Staged turbine flow returns to the chamber.
Turbine input temperature and pressure ratio use total quantities; the ratio is
inlet/outlet. Pumped modes require their component inputs; pressure-fed sizing
can omit them and reports whether inlet pressure meets the declared requirement.

Results add retained `inputs`, `pump_mass_flow_kg_s`, `energy_balance_scope`,
`numerical_policy`, `inlet_pressure_margin_pa`, its representability flag and
nullable `pressure_fed_supply_sufficient`. An unrepresentable optional margin is
null, while the exact sign still decides supply sufficiency. Warnings and scope
must be retained: algebraic zero residual does not validate integrated hardware
or hot-gas pressure closure. See `CYCLE_COMPONENT_REVIEW.md`.

## Generic delta-v uncertainty

`POST /api/v1/uncertainty` has model metadata version 2.0. It accepts unique
`specific_impulse_s` and `mass_ratio` uniform ranges, sample count 2..100000,
an integer seed, and `latin-hypercube` or `monte-carlo`. Form aliases remain
supported. Isp must be positive and mass ratio at least one over each full range.
JSON integer seeds and integral form text retain their exact integer value.

Results add retained `parameters`, `sampling_numerics`, `summary_numerics` and
inference `warnings`. A constant input/output makes Pearson `correlation`,
`absolute_correlation` and `rank` null, with `undefined_reason` set to
`constant_input_or_output`. Clients must handle that null value. P05/P50/P95
are empirical linear sample percentiles, not confidence intervals. See
`STUDY_SAMPLING_REVIEW.md` for arithmetic and validation limits.

## Loss-factor budget

`POST /api/v1/loss-budget` has model metadata version 2.0. Existing calls default
to `operating_basis="fixed-mass-flow-lumped"`; its discharge input is a legacy
performance multiplier, explicitly warned as not physical throat Cd.
`operating_basis="fixed-pressure-area"` selects separate intrinsic cstar,
mass-flow Cd and net Cf factors. It exports distinct thrust, Isp and flow ratios;
signed Isp changes can be negative under independently varied coefficients.
The model does not close actual gas/nozzle compatibility or validate coefficients.

Optional `correlation_matrix` is a numeric 4x4 JSON array, ordered combustion,
discharge, boundary_layer, divergence. Default identity declares zero correlation.
The matrix must be exactly symmetric, unit diagonal and positive semidefinite;
singular cases are allowed. `relative_uncertainty` is the same relative standard
input uncertainty per enabled factor; ideal reference uncertainty is excluded.
Disabled factors supply zero sensitivity. Distinct first-order relative proxies
are returned for thrust, Isp and flow. Bands carry `coverage_qualified=false`
and `uncertainty_kind="first_order_relative_factor_covariance_proxy"`.

Inputs, ordered terms, basis and numerical policy are retained and replayed.
Directly constructing old result/term records without retained inputs is a typed
migration error. Ordinary calculation calls retain their signature and defaults.
See [the loss-factor review](LOSS_FACTOR_REVIEW.md) for equations and limits.

## Two-phase exit-state scenarios

`POST /api/v1/two-phase-band` has metadata version 2.0. Optional
`reference_pressure_thrust_n` defaults to zero, a warned assumption, and is kept
separate from momentum thrust. `coupling_definition` defaults to
`exit-velocity-ratio`; `specific-kinetic-energy-ratio` means particle/gas kinetic
energy **per unit mass**, whose square root supplies the velocity ratio. Existing
coupling-efficiency field names remain aliases for the selected ratio in [0,1].

Reference total thrust/Isp and reference momentum
`ideal_thrust_n-reference_pressure_thrust_n` must be positive. Signed pressure
and nonpositive conditional endpoint thrust are retained with warnings; actual
gas/nozzle compatibility is unsolved. Separate particle velocity, momentum and
total-thrust-ratio ranges are returned. The legacy `performance_efficiency_range`
is a signed total-thrust-ratio alias. Results retain typed inputs and replayed
fields, with `band_kind="conditional_exit_state_scenario_range"`,
`physical_bound_qualified=false` and `statistical_coverage_qualified=false`.
See [the two-phase review](TWO_PHASE_REVIEW.md) for scope and validation gaps.

## Prescribed signed nozzle heat transfer

`POST /api/v1/nozzle-thermal` metadata version 2.0 retains consumed `inputs`
and replayed station records. Heat flux is signed positive into the wall;
`contains_negative_heat_flux` reports reversed local heat flow. The signed
surface integral uses a piecewise linear flux/radius reconstruction.
`physical_model_qualified=false` and `quadrature_error_bound_w=null` describe
its qualification. Zero `energy_balance_residual_w` is labelled by `residual_kind`
as an integral ledger identity, not measured energy closure.

If optional cooling is requested with a nonpositive net gas-side load, the
current positive-absorption bulk model fails with `unsupported_nonpositive_bulk_heat`.
Gas/wall properties remain prescribed; neither the wall temperature nor a
conjugate wall/coolant solution is predicted here. See [the Bartz review](BARTZ_REVIEW.md);
the coupled solve is `/api/v1/regen-chamber` ([REGENERATIVE_COOLING.md](REGENERATIVE_COOLING.md)).

The nested `cooling` result has `model_version="2.0"`, owned inputs and replayed
outputs. `laminar_darcy_poiseuille_number` defaults to 64 (a circular-equivalent
assumption) and sets Darcy f*Re in the laminar branch. Exact Re determines the
reported `flow_regime`; transition retains a labelled unqualified extension.
`physical_model_qualified` and `friction_correlation_domain_qualified` are false.

`temperature_limit_exceeded` and its critical warning use the exact energy
ledger, even if the exported outlet equals the limit. The optional
`temperature_limit_margin_k` and `displayed_temperature_energy_residual_w` may
be null with their corresponding `*_unrepresentable` flags. Zero
`energy_balance_residual_w` is a labelled bulk-ledger identity. These are not
boiling, coupled-cooling or hardware qualification claims. See
[the regenerative review](REGENERATIVE_REVIEW.md).

## Conditional finite-grid staging

`POST /api/v1/staging` metadata version 2.0 retains validated `inputs`, exact
integer `allocation_units`, `refined_allocation_units`, and `refined_resolution`.
The optimizer minimizes the constant-Isp/structural-fraction initial mass on the
requested exact delta-v grid. `infeasible_staging_grid` concerns this grid and
does not alone prove continuous-model infeasibility. Stage names require text.

`refinement_performed` and `refinement_grid_nested` qualify the legacy
`convergence_delta_initial_mass_relative` diagnostic. Its value is `null` at
resolution 800, where no different comparison grid exists. The difference
between grids is not a continuous-optimum error bound, and
`continuous_optimum_qualified=false`. Stage values and complete results replay
from owned inputs; rounded displayed stage delta-v values can sum differently
from the exact target ledger. See [the staging review](STAGING_REVIEW.md).

## Standalone constant and profiled burns

`POST /api/v1/engine-operating-point` builds a validated L0 operating point.
Set `kind` to `bipropellant` with delivered thrust, system Isp and total O/F,
or to `monopropellant` with delivered thrust and system Isp. Bipropellant
requests may declare an explicit dumped fuel or oxidizer stream; the response
keeps chamber flow and total tank drain distinct.

`POST /api/v1/burn/simulate` accepts a `definition` and an explicit
`operating_point`. The operating point contains named tank streams; its system
specific impulse must close against delivered thrust and total tank drain.
Without `schedule`, the response is a canonical `rocket_propulsion_burn_v1`
L0 artifact. An optional `schedule.segments` array selects L1; every segment
contains `duration_s`, `start_throttle`, `end_throttle`, and `phase`, and the
response is `rocket_propulsion_burn_v2`. Both include conservation residuals,
deterministic hashes and enough inputs for physical reproduction.

`POST /api/v1/burn/report` accepts `{ "burn_result": <canonical artifact> }` and
returns a printable, self-contained HTML engineering report. It validates the
artifact and its hashes before rendering it. `POST /api/v1/burn/export/sidera`
maps eligible constant results to the Sidera contract. Its default
`mode="exact"` rejects varying profiles. Explicit
`mode="equivalent_constant"` preserves duration, total impulse and total
propellant while recording variation, centroid offset and approximation
warnings in the manifest. Frame, direction and start time are export-only
inputs and never enter propulsion physics.

`GET /api/v1/integrations/sidera/capabilities` reports whether native Sidera
is installed and whether its public constant-burn contract is compatible.
Standalone analysis, JSON/workspace persistence and reports do not require it.
Metadata version 1.1 returns `status="checking"`, `probe_complete=false`
while inspection runs in the background; capability flags are omitted until
known. Poll again for a terminal `supported`, `incompatible` or `unavailable`
record with `probe_complete=true`. The completed snapshot lasts for this
server process; restart after installation/path changes. See
[capability inspection behavior](HTTP_CAPABILITY_PROBE.md).

The complete request shape is available in
[`examples/constant_burn_input.json`](../examples/constant_burn_input.json) and
[`examples/profiled_burn_input.json`](../examples/profiled_burn_input.json).
The same file works through the CLI:

```powershell
rocket-propulsion burn simulate examples/constant_burn_input.json `
  --output burn-result.json --report burn-report.html
rocket-propulsion burn validate burn-result.json
rocket-propulsion burn export burn-result.json --target sidera `
  --output sidera-plan.json --start 100 --direction 1 0 0 --frame inertial
```

## Transient propulsion core and SASS sizing

Ten `POST` routes in `rocket_propulsion.api.transient_routes` expose the
isolated transient core and the SASS spacecraft sizing study. Every field is
optional: `{}` runs the case of the named example script with the same
numbers. Unknown fields are refused with `InputError` code `unknown_field`
(HTTP 400); so are values above a work limit (code `work_limit`). Library
domain errors, including non-converged blowdown meshes and feed coupling
iterations (the library raises `DomainError` for these), are HTTP 400; a
`ConvergenceError` (for example from the pulse-train burn metrics) is HTTP 422.
All units are SI.

Every response `data` has the same frame:

| Key | Content |
| --- | --- |
| `model` | Label of the library model that produced the result. |
| `inputs` | Every resolved input, defaults included, for reproduction. |
| `summary` | Scalar results. |
| `series` | Plot groups `{x, source_point_count, point_count, columns}`. Every column has `point_count` entries (at most `max_points`, default and ceiling 400). Longer histories keep evenly spaced indices including the first and last point. |
| `reference_ids`, `warnings` | Source identifiers and model warnings from the library. |

The routes return library results only. The regulated-feed, feed-line and
timeline routes evaluate state right-hand sides; the library integrates no
feed or vehicle time history, so these routes return state evaluations and
sweeps, not integrated feed or Δv histories.

### `/api/v1/transient/blowdown`

Example: `examples/pressure_fed_blowdown.py`. Calls
`simulate_pressure_fed_blowdown`.

| Field | Default | Unit |
| --- | --- | --- |
| `delivered_thrust_n`, `system_specific_impulse_s`, `chamber_pressure_pa` | 100, 220, 1.5e6 | N, s, Pa |
| `initial_pressure_pa`, `initial_ullage_volume_m3`, `initial_propellant_mass_kg` | 2.0e6, 0.01, 10 | Pa, m3, kg |
| `propellant_density_kg_m3`, `polytropic_exponent` | 1000, 1.2 | kg/m3, 1 |
| `reference_supply_pressure_pa`, `minimum_supply_pressure_pa` | 2.0e6, 1.0e6 | Pa |
| `maximum_supply_pressure_pa` | null (reference pressure) | Pa |
| `thrust_pressure_exponent`, `specific_impulse_pressure_exponent`, `chamber_pressure_exponent` | 1.0, 0.02, 1.0 | 1 |
| `reserve_mass_kg`, `cant_efficiency` | 0, 1 | kg, 1 |
| `maximum_mass_step_kg`, `relative_tolerance` | 0.5, 1e-6 | kg, 1 |
| `maximum_refinements`, `max_points` | 12 (0-12), 400 (2-400) | - |

The mesh defaults are the integration controls of
`examples/blowdown_uncertainty.py`. The library defaults (0.1 kg, 1e-7) take
about 3 s. `relative_tolerance` must be at least 1e-9, and
`initial_propellant_mass_kg / maximum_mass_step_kg` at most 2000.
Response: `summary` (duration, delivered and axial impulse, consumed
propellant, thrust centroid, equivalent Isp, cutoff reason, initial and final
supply pressure and thrust, artifact hash), `evidence` (mesh convergence
record), `series.history` (x `time_s`; `remaining_propellant_mass_kg`,
`ullage_volume_m3`, `supply_pressure_pa`, `delivered_thrust_n`,
`axial_thrust_n`, `system_specific_impulse_s`, `total_mass_flow_kg_s`).

### `/api/v1/transient/regulated-feed`

Example: `examples/coupled_feed_propagator.py`. Calls
`evaluate_coupled_feed_propulsion` and `RegulatedFeedSystem.evaluate`.

| Field | Default | Unit |
| --- | --- | --- |
| `pressurant_gas_constant_j_kg_k`, `pressurant_gamma`, `pressurant_name` | 2077.1, 1.667, `"helium"` | J/kg/K, 1 |
| `tank_internal_volume_m3`, `propellant_density_kg_m3`, `bottle_volume_m3` | 0.02, 1000, 0.01 | m3, kg/m3, m3 |
| `regulator_set_pressure_pa`, `regulator_full_open_error_pa` | 2.1e6, 0.2e6 | Pa |
| `regulator_maximum_area_m2`, `regulator_discharge_coefficient`, `valve_time_constant_s` | 1e-7, 0.8, 0.05 | m2, 1, s |
| `feed_line_resistance_pa_s2_kg2` | 1000 | Pa s2/kg2 |
| `tank_heat_transfer_w_k`, `tank_wall_temperature_k`, `bottle_heat_transfer_w_k`, `bottle_wall_temperature_k` | 0, 300, 0, 300 | W/K, K |
| `minimum_propellant_mass_kg`, `minimum_injector_pressure_pa`, `minimum_bottle_pressure_margin_pa`, `maximum_tank_pressure_pa` | 1, 1.0e6, 0.1e6, 2.3e6 | kg, Pa |
| `propellant_mass_kg`, `tank_pressure_pa`, `tank_gas_temperature_k` | 10, 2.0e6, 300 | kg, Pa, K |
| `bottle_pressure_pa`, `bottle_gas_temperature_k`, `regulator_opening` | 20e6, 300, 0.5 | Pa, K, 1 |
| `relative_time_s`, `vehicle_mass_kg`, `ambient_pressure_pa` | 2, 100, 50000 | s, kg, Pa |
| `initial_mass_flow_kg_s`, `relaxation` | 0.2, 0.8 | kg/s, 1 |
| `demand_sweep_maximum_flow_kg_s`, `sweep_points` | 0.5, 101 (2-400) | kg/s, - |

The 500 N pressure-dependent profile and its supply/ambient-pressure grid are
the fixed synthetic data of the example. The 0.5 kg/s sweep bound is the line
flow bound of `examples/dynamic_feed_line_propagator.py`.
Response: `summary.propulsion` (axial thrust, force vector, engine mass flow,
performance scales, operating conditions), `summary.feed` (tank, bottle and
injector pressure, line drop, regulator command, flow and choking, heat
flows, six-state `derivative`), `evidence` (coupling iterations and
residuals), `series.demand_sweep` (x `propellant_mass_flow_kg_s`;
`injector_pressure_pa`, `feed_line_pressure_drop_pa`,
`tank_pressurant_mass_rate_kg_s`) and `series.regulator_sweep` (x
`tank_pressure_pa` across `[set - full_open_error, set]`; `regulator_command`,
`regulator_mass_flow_kg_s`, `regulator_choked`).

### `/api/v1/transient/feed-line`

Example: `examples/dynamic_feed_line_propagator.py`. Calls
`evaluate_dynamic_feed_propulsion` and `DynamicFeedLine.evaluate`. It
accepts every regulated-feed model field above, with the same defaults except
`feed_line_resistance_pa_s2_kg2 = 0` (the dynamic line owns the resistance).
The coupling and demand-sweep fields are not accepted.

| Field | Default | Unit |
| --- | --- | --- |
| `line_inertance_pa_s2_kg`, `line_compliance_kg_pa` | 1e6, 1e-8 | Pa s2/kg, kg/Pa |
| `line_resistance_pa_s2_kg2` | (2.0e6 - 1.8e6) / 0.217^2 | Pa s2/kg2 |
| `minimum_manifold_pressure_pa`, `maximum_manifold_pressure_pa` | 1.0e6, 2.3e6 | Pa |
| `maximum_absolute_flow_kg_s` | 0.5 | kg/s |
| `line_mass_flow_kg_s`, `manifold_pressure_pa` | 0.217, 1.8e6 | kg/s, Pa |
| `sweep_points` | 101 (2-400) | - |

Response: `summary.propulsion`, `summary.feed`, `summary.line` (upstream
pressure, engine demand, friction drop, line storage rate, mass closure
error, natural angular frequency, damping ratio, `derivative`), model hashes,
and `series.line_flow_sweep` (x `line_mass_flow_kg_s` from 0 to the flow
bound; `line_flow_rate_kg_s2`, `manifold_pressure_rate_pa_s`,
`friction_pressure_drop_pa`, `damping_ratio` at the evaluated upstream
pressure and engine demand).

### `/api/v1/transient/electric`

Example: `examples/electric_propulsion.py`. Calls
`PowerLimitedPropulsionModel.select_level` and `evaluate`.

| Field | Default | Unit |
| --- | --- | --- |
| `levels` | low/mid/high: 1000/2000/3000 W, 0.040/0.080/0.110 N, 1.5e-6/2.6e-6/3.2e-6 kg/s | objects `{level_id, input_power_w, thrust_n, mass_flow_kg_s}` |
| `source_id`, `source_sha256` | `"example-synthetic-throttle-table"`, 64 zeros | - |
| `objective` | `"thrust"` (or `"specific_impulse"`) | - |
| `up_switch_hysteresis_w`, `power_margin_w` | 50, 100 | W |
| `propellant_reserve_kg`, `duty_cycle` | 2, 1 | kg, 1 |
| `available_power_w`, `vehicle_mass_kg`, `propellant_mass_kg` | 2400, 500, 20 | W, kg, kg |
| `sweep_minimum_power_w`, `sweep_maximum_power_w` | 0, 2 x `available_power_w` | W |
| `sweep_points` | 101 (2-400) | - |

The sweep range is a display choice; the example has no sweep. Response:
`summary` (selected level, thrust, Isp, effective and consumed power, power
deficit, propellant rate, acceleration, model and table hashes), `levels`
(with Isp and total efficiency), `switching_surfaces` of the selected level,
and `series.power_sweep` (x `available_power_w`; `level_index`, `thrust_n`,
`specific_impulse_s` (null when off), `consumed_power_w`, `mass_flow_kg_s`,
and the hysteresis-aware `level_index_increasing`/`level_index_decreasing`).

### `/api/v1/transient/pulse-mode`

Example: `examples/pulse_mode_trim.py`. Calls `PulseModeThruster` and
`free_space_burn_metrics`.

| Field | Default | Unit |
| --- | --- | --- |
| `on_times_s` | 0.010, 0.020, 0.050, 0.100, 0.200 | s |
| `impulse_bits_n_s`, `propellant_per_pulse_kg` | example tables `t - 0.004` and `(t + 0.003) / (g0 220)` | N s, kg |
| `steady_thrust_n`, `steady_mass_flow_kg_s` | 1, 1 / (g0 220) | N, kg/s |
| `minimum_off_time_s`, `maximum_on_time_s` | 0.02, 2 | s |
| `source_id`, `source_sha256` | `"synthetic-1N"`, 64 `c` | - |
| `initial_mass_kg`, `propellant_budget_kg`, `delta_v_m_s` | 12, 2, 0.05 | kg, kg, m/s |
| `quantization` | `"nearest"` (or `"no_overshoot"`) | - |
| `pulse_off_time_s` | 0.05 (period = on-time + off-time) | s |
| `quantization_requests_m_s` | 0.0002, 0.0004, 0.05, 1.0 (at most 16) | m/s |
| `sweep_points`, `max_points` | 101, 400 | - |

Response: `summary` (minimum impulse bit, minimum and longest on-time, Isp at
the minimum bit, `plan` and `pulse_train` with free-space Δv), `quantization`
(one plan per request: delivered, residual, quantized flag, pulse count,
on-time, propellant), `series.impulse_bit_curve` (x `on_time_s` geometrically
spaced; `impulse_n_s`, `propellant_kg`, `specific_impulse_s`) and
`series.pulse_train` (x `time_s`; thrust and mass flow at segment ends).

### `/api/v1/transient/uncertainty`

Example: `examples/blowdown_uncertainty.py`. Calls
`simulate_blowdown_uncertainty`. Accepts every blowdown model field (same
defaults, except `maximum_supply_pressure_pa = 2.2e6`) and:

| Field | Default | Unit |
| --- | --- | --- |
| `parameters` | initial pressure 1.9e6-2.1e6 truncated normal, sigma 40000, aleatory; rated thrust 97-103 triangular, epistemic | objects `{name, lower, upper, distribution, uncertainty_class, standard_deviation?, source, rationale, nominal?}` |
| `sample_count` | 16 (4-256) | - |
| `seed` | 42 | - |
| `method` | `"latin-hypercube"` (or `"monte-carlo"`) | - |
| `convergence_tolerance` | 0.05 | 1 |
| `maximum_mass_step_kg` | 0.5 | kg |
| `integration_relative_tolerance` | 1e-3 (at least 1e-6) | 1 |
| `histogram_bins` | 10 (1-50) | - |

A parameter's `nominal` defaults to the resolved model input of that name.
The example runs 32 samples at a 1e-6 mesh tolerance, about 45 s. The route
defaults (16 samples, 1e-3) answer in about 1 s, and the realizations are
still reproducible for a given seed. Response: `summary` (method, seed,
sample count, aleatory/epistemic counts, `sampling_stable`), `statistics`
(per metric: mean, standard deviation, mean interval, minimum, P05, P50, P95,
maximum), `sensitivity` (Pearson screening), `convergence`, `histograms`
(`edges`, `counts` per metric) and `series.realizations` (x `index`; sampled
`input_<name>` columns, `duration_s`, `delivered_impulse_n_s`,
`consumed_propellant_kg`, `thrust_centroid_time_s`,
`final_supply_pressure_pa`, `cutoff_reason`).

### `/api/v1/transient/scenarios`

Example: `examples/propulsion_scenarios.py`. Calls
`build_propulsion_scenario_ensemble`.

| Field | Default | Unit |
| --- | --- | --- |
| `delivered_thrust_n`, `system_specific_impulse_s`, `chamber_pressure_pa` | 100, 220, 1.0e6 | N, s, Pa |
| `engine_ids` | `["engine-a", "engine-b"]` | - |
| `burn_duration_s`, `throttle` | 10, 1 | s, 1 |
| `scenarios` | nominal 0.90; engine-a out 0.08; common-cause loss of both 0.02 | objects `{scenario_id, probability, events, source, rationale}` |
| `objective` | `"illustrate propulsion consequence weighting"` | - |
| `maximum_step_s`, `max_points` | 0.25, 400 | s, - |

Each event is `{event_type, engine_ids, value, source, rationale}` with
`event_type` one of `engine-disabled` (value null), `ignition-delay`,
`early-cutoff` (s), `thrust-scale`, `specific-impulse-scale`,
`response-time-scale` (factor). Scenarios and events must state `source` and
`rationale`. The example probabilities are illustrations, not failure rates.
Response: `summary` (zero-thrust probability), `statistics` (probability
weighted expected value, P05/P50/P95, extremes per metric), `scenarios`
(probability, duration, impulse, propellant, centroid, artifact hash),
`series.consequences` (x `scenario_id`; `probability`,
`delivered_impulse_n_s`, `consumed_propellant_kg`, `duration_s`) and
`scenario_histories` (one thrust/mass-flow series per scenario).

### `/api/v1/transient/timeline`

Example: `examples/propulsion_timeline.py`. Calls `PropulsionTimeline`.

| Field | Default | Unit |
| --- | --- | --- |
| `burns` | perigee at 10 s for 120 s, apogee at 400 s for 40 s; each 30000 N, 7.8 kg/s oxidizer, 2.2 kg/s fuel | objects `{burn_id, ignition_time_s, duration_s, thrust_n, oxidizer_mass_flow_kg_s, fuel_mass_flow_kg_s}` |
| `jettisons` | 150 kg payload adapter at 200 s (`[]` allowed) | objects `{event_id, time_s, mass_kg, label}` |
| `initial_vehicle_mass_kg`, `oxidizer_mass_kg`, `fuel_mass_kg` | 6000, 2400, 700 | kg |
| `end_time_s` | 500 | s |
| `allow_overlapping_burns`, `max_points` | false, 400 | - |

Response: `summary` (total impulse, propellant, jettisoned mass, ledger final
mass, timeline hash), `burns` (per-burn ignition, cutoff and artifact
integrals), `events` (burn start/stop and jettison roots) and
`series.history` (x `time_s`, with left and right limits at every timed
event; `axial_thrust_n`, `mass_flow_kg_s`, `active_burn_count`,
`jettisoned_mass_kg`). The timeline computes no orbit, so no acceleration or
Δv history is returned.

### `/api/v1/transient/calibration`

Example: `examples/hot_fire_calibration.py`. Calls
`calibrate_hot_fire_scale`.

| Field | Default | Unit |
| --- | --- | --- |
| `parameter_name` | `"thrust_scale"` (or `specific_impulse_scale`, `time_scale`, `axial_efficiency_scale`) | - |
| `measurements` | six synthetic engines x three runs: 1.008 + engine offset + run deviation, u = 0.001 | objects `{engine_id, run_id, condition_id, value, standard_uncertainty, source}` |
| `source`, `rationale` | example strings | - |
| `relative_tolerance`, `maximum_iterations` | 1e-10, 200 (1-1000) | 1, - |
| `coverage_factor` | 3 | 1 |

Response: `summary` (consensus scale and standard uncertainty, observed and
process variance components, within/between fractions, calibration hash),
`evidence` (Mandel-Paule solver record), `residual_parameters` (`population`
and `same_engine` truncated-normal scales with bounds), `engines`,
`series.engines` (x `engine_id`; mean, minimum, maximum scale) and
`series.measurements` (x `engine_index`; `value`, `standard_uncertainty`).

### `/api/v1/sass/mission`

Example: `examples/sass_mission.json`. Calls
`rocket_propulsion.studies.sass.project_from_dict`.

| Field | Default | Unit |
| --- | --- | --- |
| `study` | the complete `rocket_propulsion_sass_study_v1` request of the example | object |
| `bus_dry_mass_kg`, `residual_fraction`, `fixed_reserve_kg` | none (50, 0.02, 0.2 in the example) | kg, 1, kg |
| `include_document` | false | - |

The three overrides replace the value of the corresponding cited scalar. Its
analyst-assumption citation is kept and its note records the override.
Response: `summary` (title, nominal and allocated Δv, eligible count, study
hash, assurance block), `candidates` (sized `metrics` with units,
eligibility, score, requirement `checks`, missing metrics, unresolved items),
`series.trade` (x `identifier`; `wet_mass_kg`, `consumed_propellant_kg`,
`elapsed_time_s`, `score`; null for unsized candidates) and, with
`include_document`, the complete hash-bearing study `document`
(about 650 kB).

## Heat transfer, regenerative cooling and combustion products

Six routes expose the heat-transfer and combustion modules. Every field has a
default, so `{}` returns a complete worked case; unknown propellant pairs and
nonpositive durations are refused with `InputError`/`DomainError`.

| Route | Main inputs (defaults) | Returns |
|---|---|---|
| `/api/v1/heat-conduction` | `hot_wall_thickness_m` (1e-3), `hot_wall_conductivity_w_m_k` (305), `hot_wall_heat_capacity_j_m3_k` (3.5e6), gas and coolant `*_heat_transfer_coefficient_w_m2_k` and `*_temperature_k`, `initial_temperature_k`, `duration_s` (0.2) | Crank–Nicolson wall history, four through-wall profiles, stored energy and the energy residual |
| `/api/v1/heat-convection` | `fluid` (methane), `temperature_k`, `pressure_pa`, `hydraulic_diameter_m` | Real-fluid transport state and h against mass flux for Gnielinski, Dittus–Boelter and Colburn |
| `/api/v1/heat-radiation` | `inlet_radius_m`, `exit_radius_m`, `length_m`, `recovery_temperature_k`, `emissivity` | Frustum view-factor matrix, radiation-cooled wall temperature and its sweep over emissivity |
| `/api/v1/regen-chamber` | `propellant_pair` (LOX_CH4), `chamber_pressure_pa`, `area_ratio`, `coolant`, `coolant_inlet_temperature_k`, `coolant_inlet_pressure_pa`, `coolant_mass_flow_kg_s`, `channel_count`, `channel_width_m`, `channel_height_m`, `hot_wall_thickness_m`, `roughness_m`, `wall_material` | Station table (gas, wall and coolant temperatures, heat flux, coolant pressure) and totals with the energy-balance residual |
| `/api/v1/combustion-products` | `propellant_pair` (LOX_LH2), `oxidizer_fuel_ratio` (pair reference), `chamber_pressure_pa`, `area_ratio` | Composition, radicals and residuals at chamber, throat and exit; dissociation energy split; stoichiometry |
| `/api/v1/combustion-losses` | the same plus `contraction_ratio` (2.5), `compute_kinetics` (LOX_LH2 only), `throat_radius_m`, `divergence_efficiency`, `boundary_layer_decrement_s`, `characteristic_velocity_efficiency` | Frozen/equilibrium limits, finite-area combustor, finite-rate nozzle (H, O, OH, H₂O along A/A*) and the ODE → delivered Isp chain |

Model details and verification: [HEAT_CONDUCTION.md](HEAT_CONDUCTION.md),
[HEAT_CONVECTION.md](HEAT_CONVECTION.md), [HEAT_RADIATION.md](HEAT_RADIATION.md),
[REGENERATIVE_COOLING.md](REGENERATIVE_COOLING.md),
[COMBUSTION_EFFICIENCY.md](COMBUSTION_EFFICIENCY.md) and
[FINITE_RATE_KINETICS.md](FINITE_RATE_KINETICS.md).

Three more routes cover transpiration cooling, combustion stability and
propellant vaporization (`api/combustion_dynamics_routes.py`); `{}` again
returns a worked case.

| Route | Main inputs (defaults) | Returns |
|---|---|---|
| `/api/v1/transpiration-cooling` | `unblown_heat_transfer_coefficient_w_m2_k` (2e4), `recovery_temperature_k` (3300), `target_wall_temperature_k` (900), `coolant_temperature_k` (100), `coolant_specific_heat_j_kg_k` (14300), `wall_efficiency` (1), `radiative_heat_flux_w_m2` (0) | Required coolant flux and wall state; wall temperature and h/h* over 0-2.5 × that flux |
| `/api/v1/combustion-stability` | `propellant_pair` (LOX_CH4), `oxidizer_fuel_ratio`, `chamber_pressure_pa` (7e6), `chamber_radius_m` (0.10), `chamber_length_m` (0.30), `contraction_ratio` (3), `oxidizer_drop_fraction`, `fuel_drop_fraction` (0.2), `combustion_time_lag_s` (8e-4) | Chamber state, acoustic modes, chugging gain and limit, Nyquist verdict and the f(iω) curve |
| `/api/v1/propellant-vaporization` | `propellant_pair` (LOX_CH4 or LOX_LH2), `chamber_pressure_pa`, `chamber_length_m` (0.25), `contraction_ratio` (3), `drop_radius_m` (7.62e-5), `injection_velocity_m_s` (30.48) | Liquid properties, generalized lengths, percent vaporized, c* efficiency and a chamber-length sweep |

Details: [TRANSPIRATION_COOLING.md](TRANSPIRATION_COOLING.md),
[COMBUSTION_INSTABILITY.md](COMBUSTION_INSTABILITY.md) and
[PROPELLANT_VAPORIZATION.md](PROPELLANT_VAPORIZATION.md).

The aerodynamics and real-fluid routes are `/api/v1/real-fluid`,
`/api/v1/airfoil`, `/api/v1/wing`, `/api/v1/supersonic-body`,
`/api/v1/flat-plate`, `/api/v1/entry` and `/api/v1/launch-vehicle-aero`; their
models are described in [REAL_FLUIDS.md](REAL_FLUIDS.md),
[AERODYNAMICS_SUBSONIC.md](AERODYNAMICS_SUBSONIC.md),
[COMPRESSIBLE_ADVANCED.md](COMPRESSIBLE_ADVANCED.md),
[AEROTHERMODYNAMICS.md](AEROTHERMODYNAMICS.md) and
[LAUNCH_VEHICLE_AERODYNAMICS.md](LAUNCH_VEHICLE_AERODYNAMICS.md).

## Near-throat Hall approximation (Python only)

`HallThroatModel(throat_radius_m=..., throat_curvature_radius_m=..., gamma=...,
axisymmetric=...)` supplies geometry/gas assumptions for
`evaluate_hall_throat(model=..., axial_parameter=..., radius_fraction=...,
radial_basis='throat', order=2)`.
`locate_hall_sonic_point(model=..., wall_fraction=..., order=2)` returns the
retained polynomial's sonic root with outward parameter, function, x and
radius bounds. Results retain nonlinear PDE/slip defects and explicit display
loss. This is an approximation with no physical error bound or complete
sonic-to-MOC startup; see [Hall review](HALL_THROAT_REVIEW.md).

## Interior characteristic intersection (Python only)

`solve_moc_field(axis_result=..., inlet_points=..., wall_segment=...,
marching_steps=...)` now assembles the local kernels across a vertical
supersonic inlet from axis to finite straight wall. `MocFieldInputs` and
`MocFieldResult` retain both staggered fronts, complete cell audits, relative
mass flux, cell imbalance and separate solved/unenforced residual diagnostics.
It needs existing supersonic axis history; it supplies no sonic initialization
or nozzle optimization and is not wired to the HTTP designer. See the
[field review](MOC_FIELD_REVIEW.md).

`solve_moc_axis_departure(axis_result=..., minus_parent=...)` carries an owned
axis-boundary result into the first off-axis point. Gamma and geometry mode
are inherited from the axis result. Cylindrical flow uses a second-order
radial Taylor relation and C- compatibility; the separately exposed
`averaged_plus_compatibility_defect` is not an enforced equation. Planar
flow retains the two-family interior process. Inputs/results are
`MocAxisDepartureInputs` and `MocAxisDepartureResult`. This local process is
not connected to the HTTP designer; see the [departure review](MOC_AXIS_DEPARTURE_REVIEW.md).

`solve_moc_axis_point` additionally takes `incoming` and exactly two
`axis_history` states. `solve_moc_wall_point` takes `incoming` and a
`MocWallSegment`; its accepted x is in `(x0,x1]`. Their owned
`MocBoundaryResult` includes signed normalized source diagnostics and always
retains false full-nozzle/physical qualification flags. These local boundary
processes are not yet connected to the HTTP nozzle designer; see the
[boundary review](MOC_BOUNDARY_REVIEW.md).

Python geometry also exposes `MocFlowPoint`, `MocUnitInputs`, `MocUnitResult`
and `solve_moc_interior_point` for a separately verified interior C+/C-
intersection. It accepts planar or axisymmetric equations with explicit
positive-radius restrictions for the latter. Owned results expose raw versus
coordinate-rounding-adjusted equation residuals and never claim full nozzle
or physical qualification. This kernel is not connected to the current HTTP
MOC nozzle generator; its construction findings remain open. See the
[unit-process review](MOC_UNIT_PROCESS_REVIEW.md).

## Sampled geometry exports (contour model metadata 2.1)

With `include_exports=true`, `/api/v1/nozzle-contour` returns CSV/SVG/DXF/STL
text plus owned consumed `inputs`, `coordinate_unit="m"`,
`mesh_role="closed_fluid_domain_visualization"`, `serialization_policy`, and
actual false `manufacturing_qualified`/`physical_model_qualified` flags. CSV,
DXF and STL preserve binary64 numbers with 17 significant digits. Reopening a
Python `ContourExportBundle` requires retained inputs and replays every file.
Legacy four-text-only construction cannot establish this replay contract.
STL carries no standardized units: select metres in the receiving application.
It represents the closed sampled fluid domain, with no nozzle wall thickness.
SVG independently scales axes; it is not a manufacturing drawing. Contour
construction and MOC flow findings remain open; see [geometry review](GEOMETRY_REVIEW.md).

