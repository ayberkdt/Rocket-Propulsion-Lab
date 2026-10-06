# Rocket Propulsion Lab

[![engineering-quality](https://github.com/ayberkdt/Rocket-Propulsion-Lab/actions/workflows/ci.yml/badge.svg)](https://github.com/ayberkdt/Rocket-Propulsion-Lab/actions/workflows/ci.yml)
![Python 3.11+](https://img.shields.io/badge/python-3.11%2B-blue)
![Dependencies: none](https://img.shields.io/badge/dependencies-none-success)
![License: MIT](https://img.shields.io/badge/license-MIT-lightgrey)

A dependency-free Python toolkit for preliminary rocket-propulsion engineering:
compressible flow, thermodynamics, thermochemistry, nozzle geometry, finite-burn
analysis and uncertainty studies. Everything runs locally from the standard
library and is exposed three ways:

- a **Python API** of pure, deterministic functions returning immutable dataclasses;
- a **versioned JSON HTTP API** (`/api/v1/...`) served by a built-in local server;
- a **browser workbench** (single-page UI, Turkish interface) on top of that API.

> This is an engineering study tool, not flight-certified analysis software.
> Every result carries explicit model metadata and assumptions; the fidelity
> boundary of each model is documented in [docs/GAP_ANALYSIS.md](docs/GAP_ANALYSIS.md).

## Contents

- [Quick start](#quick-start)
- [Usage](#usage)
- [Capabilities](#capabilities)
- [Project layout](#project-layout)
- [Documentation](#documentation)
- [Development](#development)
- [Conventions and scope](#conventions-and-scope)
- [License](#license)

## Quick start

Requirements: Python 3.11 or newer. There are no runtime dependencies.

```bash
git clone https://github.com/ayberkdt/Rocket-Propulsion-Lab.git
cd Rocket-Propulsion-Lab
python -m pip install -e ".[dev]"
rocket-propulsion
```

Then open <http://127.0.0.1:8765>. Use `--host` and `--port` to change the bind
address.

Without installing anything, the PowerShell launchers set `PYTHONPATH=src`
and run the same entry points:

```powershell
.\scripts\run.ps1     # start the local server and browser UI
.\scripts\test.ps1    # run the test suite
.\scripts\build.ps1   # build a wheel into dist/
```

## Usage

### Python

```python
from rocket_propulsion.compressible import calculate_off_design_nozzle
from rocket_propulsion.geometry import generate_rao_contour
from rocket_propulsion.propulsion import calculate_loss_budget

contour = generate_rao_contour(throat_area_m2=0.01, area_ratio=25.0)
nozzle = calculate_off_design_nozzle(
    chamber_pressure_pa=5_000_000.0,
    chamber_temperature_k=3500.0,
    throat_area_m2=contour.throat_area_m2,
    area_ratio=contour.area_ratio,
    ambient_pressure_pa=101_325.0,
    gamma=1.22,
    gas_constant_j_kg_k=355.0,
)
losses = calculate_loss_budget(
    ideal_thrust_n=nozzle.thrust_n,
    ideal_specific_impulse_s=nozzle.specific_impulse_s,
    divergence_efficiency=contour.estimated_divergence_efficiency or 1.0,
)
print(nozzle.regime, nozzle.thrust_n, losses.delivered_thrust_n)
```

More complete, reference-labelled workflows live in [`examples/`](examples/):

| Example | What it shows |
|---|---|
| `quick_start.py` | Rao contour, off-design nozzle and loss budget in a dozen lines |
| `pressure_fed_blowdown.py` | Preliminary pressure-fed blowdown with polytropic tank decay |
| `blowdown_uncertainty.py` | Sourced LHS/Monte Carlo ensemble over a blowdown case |
| `propulsion_scenarios.py` | Nominal, engine-out and common-cause loss histories |
| `hierarchical_propulsion_ensemble.py` | Discrete outcomes combined with conditional continuous variation |
| `hot_fire_calibration.py` | Random-effects calibration from repeated hot-fire data |
| `hot_fire_multivariate_calibration.py` | Paired thrust/Isp/duration residual covariance model |
| `constant_burn_input.json`, `profiled_burn_input.json` | Request bodies for the burn API and CLI |
| `two_case.rplab.json` | A versioned two-case workspace |

### HTTP API

Every calculation is a `POST` with a JSON object. Success returns `data` plus a
`meta` block naming the model, its version, assumptions and units. Failures
return a stable `error` envelope: HTTP 400 for input/domain errors, 422 for
convergence failures and 501 when an optional feature is unavailable.

```bash
curl -s http://127.0.0.1:8765/api/v1/isentropic -H "Content-Type: application/json" -d "{\"input_kind\": \"mach\", \"value\": 2, \"gamma\": 1.4}"
```

`GET /api/v1/health` reports the service version. The full endpoint list and
envelope contract are in [docs/API.md](docs/API.md).

### Command line

The burn-analysis slice is also available as a CLI:

```bash
rocket-propulsion burn simulate examples/constant_burn_input.json --report burn.html
rocket-propulsion burn validate result.json
rocket-propulsion burn tabulate-trace examples/solid_thrust_trace.csv --output trace.json
rocket-propulsion burn export result.json --target sidera --start 0 --direction 1 0 0 --output handoff.json
```

## Capabilities

**Compressible flow**
- Isentropic perfect-gas relations with forward and inverse conversions, and parameter sweeps.
- Normal shocks; weak/strong oblique-shock branches with detachment limits.
- Fanno and Rayleigh constant-area duct flow.
- Ideal choked converging-diverging nozzle performance.
- Off-design nozzles: back-pressure regimes, internal normal shocks, altitude sweeps and a source-labelled Summerfield separation screen.

**Thermodynamics and thermochemistry**
- Ideal-gas states, static/stagnation connections, compression/expansion processes and polytropic relations.
- NASA Glenn temperature-dependent properties for ten combustion species, with frozen ideal-gas mixtures and inverse enthalpy/entropy temperature solutions.
- An adapter for an externally installed NASA CEA executable (`RPLAB_CEA_EXECUTABLE`; CEA is not bundled) and an explicitly experimental local Gibbs equilibrium provider with O/F sweeps.

**Nozzle geometry**
- Conical, preliminary Hermite, Rao/TOP and method-of-characteristics contours with station-by-station flow solutions and MOC residual/mesh-refinement evidence.
- CSV, SVG, DXF and STL export.

**Propulsion systems**
- Switchable performance loss budgets, Bartz wall heat flux, experimental regenerative channel balance and condensed-phase (two-phase) performance bands.
- Tsiolkovsky forward/inverse calculations, thrust and multi-stage mass optimisation.
- Pressure-fed, gas-generator and staged-combustion component energy balances.
- A screening catalog of common liquid, solid and hybrid propellant systems.

**Finite-burn analysis**
- L0 constant and L1 piecewise-linear burns with explicit tank streams, system-Isp closure and exact analytic segment integrals.
- Transient development: first-order engine response, restart/cooldown constraints, engine clusters, imported solid/test thrust curves and a hashed, discontinuity-safe tabulated thrust/mass-flow contract.
- Preliminary pressure-fed blowdown with polytropic ullage decay, pressure/reserve cutoff, power-law or measured pressure-performance maps and convergence evidence.

**Uncertainty, scenarios and calibration**
- Seeded Monte Carlo and Latin-hypercube ensembles with bounded marginals, aleatory/epistemic separation, optional Gaussian-copula correlation, convergence diagnostics and sensitivity ranking.
- Discrete propulsion scenarios (engine-out, common-cause loss, ignition delay, early cutoff, degradation) with exhaustive probabilities.
- Hierarchical ensembles that combine discrete outcomes with scenario-conditional continuous variation, exact joint weights and within/between-scenario variance attribution.
- Hot-fire random-effects and paired covariance calibration that turn repeated test data into hash-linked residual ensemble parameters.

**Studies and integration**
- Versioned `.rplab.json` workspaces with migrations, case comparison and JSON/CSV/SVG/HTML reports.
- Optional, fail-closed hand-off of burn artifacts to the Sidera astrodynamics framework. Sidera is never imported by the core packages and is not a dependency.

## Project layout

```text
src/rocket_propulsion/
├── api/              HTTP transport, route registry and static-file server
├── compressible/     Isentropic, shock, duct and nozzle flow physics
├── core/             Errors, metadata, units and numerical solvers
├── geometry/         Parametric nozzle contours, MOC and export
├── integrations/     Optional, fail-closed Sidera hand-off adapters
├── propulsion/       Rocket equations, cycles, losses, propellants, burns/
├── studies/          Workspaces, comparison, sampling and reports
├── thermochemistry/  NASA species data, mixtures, equilibrium, CEA adapter
├── thermodynamics/   Perfect-gas states, processes and heat transfer
└── web/              Browser UI (index.html, app.js, styles.css)
tests/                Unit, reference-data and HTTP-level tests
docs/                 Architecture, models, validation and contracts
examples/             Runnable workflows and sample inputs
scripts/              PowerShell launchers for run, test and build
```

The architecture is a simple ports-and-adapters boundary: browser UI to JSON
API to domain packages to core utilities. Domain packages never import the API
or web layers. Every package has its own README, and public functions document
their equations, parameters and domain restrictions in docstrings so the code
and documentation cannot drift apart.

## Documentation

| Document | Purpose |
|---|---|
| [USER_GUIDE.md](docs/USER_GUIDE.md) | End-to-end browser workflow and fidelity labels |
| [API.md](docs/API.md) | JSON envelope, error codes, Python and PowerShell examples |
| [ARCHITECTURE.md](docs/ARCHITECTURE.md) | Dependency direction, design rules and extension points |
| [GAP_ANALYSIS.md](docs/GAP_ANALYSIS.md) | What is solved versus what needs higher-fidelity models |
| [REFERENCES.md](docs/REFERENCES.md) | Primary equation and data sources |
| [TOLERANCES.md](docs/TOLERANCES.md) | Numerical acceptance thresholds used by the tests |
| [VALIDATION_REPORT.md](docs/VALIDATION_REPORT.md) | Closed evidence and the remaining external CEA gate |
| [BLOWDOWN_MODEL.md](docs/BLOWDOWN_MODEL.md) | Pressure-fed equations, limits and bibliography |
| [TABULATED_BURN_CONTRACT.md](docs/TABULATED_BURN_CONTRACT.md) | Discontinuity-safe thrust/mass-flow exchange format |
| [PROPULSION_UNCERTAINTY.md](docs/PROPULSION_UNCERTAINTY.md) | Ensemble contract and ownership boundary |
| [PROPULSION_SCENARIOS.md](docs/PROPULSION_SCENARIOS.md) | Probability-labelled failure and execution scenarios |
| [HIERARCHICAL_PROPULSION_ENSEMBLE.md](docs/HIERARCHICAL_PROPULSION_ENSEMBLE.md) | Conditional composition of scenarios and continuous uncertainty |
| [HOT_FIRE_CALIBRATION.md](docs/HOT_FIRE_CALIBRATION.md) | Repeated-test calibration equations and references |
| [SIDERA_BURN_INTEGRATION_PLAN.md](docs/SIDERA_BURN_INTEGRATION_PLAN.md), [SIDERA_ENSEMBLE_HANDOFF.md](docs/SIDERA_ENSEMBLE_HANDOFF.md) | Hand-off contracts toward the Sidera framework |
| [ROADMAP.md](docs/ROADMAP.md), [DEVELOPMENT_PLAN.md](docs/DEVELOPMENT_PLAN.md), [RELEASE_NOTES.md](docs/RELEASE_NOTES.md) | Scope, staged plan and change history |

## Development

```bash
python -m pip install -e ".[dev]"
python -m ruff check .
python -m unittest discover -s tests -v
```

`python -m pytest` works as well; the pytest configuration in `pyproject.toml`
already points at `src` and `tests`.

The CI workflow ([`.github/workflows/ci.yml`](.github/workflows/ci.yml)) runs
lint, tests and a distribution build on Python 3.11 and 3.12. One test is
skipped unless a compatible Sidera installation is present on the explicit
test path.

Design rules that every contribution should respect are listed in
[docs/ARCHITECTURE.md](docs/ARCHITECTURE.md): keep equations pure and
deterministic, validate physical domains at module boundaries, add
reference-value tests with every new relation, keep the HTTP API versioned,
and version saved workspaces with explicit migrations.

## Conventions and scope

- SI units for all dimensional quantities; angles in degrees.
- `gamma` must be greater than one; Mach numbers must be positive.
- Results are labelled with model fidelity. Screening results (for example the
  Summerfield separation check or propellant catalog ranges) are warnings and
  starting points, not design answers.
- Nothing is sent to an external service; the server binds to localhost by default.

## License

MIT. See [LICENSE](LICENSE).
