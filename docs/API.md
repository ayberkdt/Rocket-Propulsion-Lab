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

