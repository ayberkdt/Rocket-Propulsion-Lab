# API

The local server exposes versioned JSON endpoints grouped by domain:

- `POST /api/v1/isentropic`
- `POST /api/v1/normal-shock`
- `POST /api/v1/oblique-shock`
- `POST /api/v1/isentropic-sweep`
- `POST /api/v1/ideal-gas`
- `POST /api/v1/stagnation`
- `POST /api/v1/thermo-process`
- `POST /api/v1/polytropic`
- `POST /api/v1/duct-flow`
- `POST /api/v1/nozzle`
- `POST /api/v1/nozzle-off-design`
- `POST /api/v1/nozzle-contour`
- `POST /api/v1/loss-budget`
- `POST /api/v1/nozzle-thermal`
- `POST /api/v1/two-phase-band`
- `POST /api/v1/rocket-equation`
- `POST /api/v1/thrust`
- `POST /api/v1/staging`
- `POST /api/v1/engine-cycle`
- `POST /api/v1/uncertainty`
- `POST /api/v1/workspace`
- `POST /api/v1/case-comparison`
- `POST /api/v1/propellants`
- `POST /api/v1/species`
- `POST /api/v1/thermal-properties`
- `POST /api/v1/equilibrium`
- `POST /api/v1/chemistry-capabilities`
- `POST /api/v1/of-sweep`

`routes.py` translates JSON-shaped dictionaries to domain function calls.
`server.py` owns HTTP concerns and static assets. Keeping those responsibilities
separate makes routing testable without opening a network port.

Successful HTTP calculations return `{ "data": ..., "meta": ... }`. The
metadata names the model, assumptions, units, validity ranges, and references.
Failures use a stable `{ "error": { "code", "message", "field?", "details?" } }`
envelope. Input/domain errors return HTTP 400; numerical non-convergence returns
HTTP 422.

