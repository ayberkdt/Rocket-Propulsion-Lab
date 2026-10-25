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
conjugate wall/coolant solution is predicted. See [the Bartz review](BARTZ_REVIEW.md).

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

