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

