# Rocket Propulsion Lab

A dependency-free Python toolkit for preliminary rocket-propulsion analysis,
and a propulsion provider for astrodynamics software.

It covers two jobs:

1. **Engineering studies:** compressible flow, thermochemistry, nozzle
   geometry, engine cycles, finite-burn sizing, hot-fire calibration and
   uncertainty, available from Python, a JSON API, a CLI and a browser
   workbench.
2. **Propulsion for orbit propagation:** event-safe force and mass-flow models
   with analytic partials, which an orbit propagator, orbit-determination
   filter or trajectory optimizer can call at every integration step. These
   cover chemical, pressure-fed and power-limited electric propulsion.

The library computes no orbits. Epochs, frames, attitude, gravity and
trajectory integration belong to the consuming astrodynamics framework (the
companion project is Sidera). This line is enforced in code: the propulsion
core never imports an astrodynamics package.

> Engineering study software, not flight-certified. Python 3.11+, MIT licence,
> no runtime dependencies.

---

## Contents

- [Highlights](#highlights)
- [Installation](#installation)
- [Quick start](#quick-start)
- [Using it from an orbit propagator](#using-it-from-an-orbit-propagator)
- [Project layout](#project-layout)
- [Documentation](#documentation)
- [Testing and quality](#testing-and-quality)
- [Validation status and limits](#validation-status-and-limits)
- [Conventions](#conventions)

## Highlights

### Propagation interface for astrodynamics

| Capability | What a propagator gets |
|---|---|
| **Transient burn bridge** | Thrust and vehicle/tank mass derivatives from a hash-verified thrust/flow history. Explicit left/right limits at discontinuities. Start, boundary, stop, dry-mass and reserve roots. O(log n) lookup. |
| **Operating-condition maps** | N-D performance surfaces over pressure, throttle, mixture ratio or power, multilinear or C¹ cubic Hermite. Analytic gradients, no extrapolation, stopping domain roots and restart roots at gradient kinks. |
| **Lunaris adapter** | Lunar-orbit propagation in [Lunaris](https://github.com/ayberkdt/lunaris): exact-mass impulsive equivalents or finite burns inside Lunaris' propagator, with bias-shifted segmentation and exact dry-mass cuts, verified end-to-end against the real Lunaris code. |
| **Ambient back-pressure** | `F = F_vac − p_a·A_e` with `∂a/∂p_a` for ascent and landing, plus a net-thrust margin root. |
| **Pressure-fed liquid stages** | Pressurant bottle, regulator, thermal ullage and feed-line ODEs, plus an inertance/compliance line model with composed Jacobians. |
| **Power-limited electric propulsion** | Hash-verified throttle tables, maximum-thrust or maximum-Isp throttle paths with hysteresis, a fixed-efficiency continuous law, and power-switching roots with deterministic mode targets. |
| **Navigation uncertainty** | Driver priors linearized from hot-fire calibration, the Gates execution-error covariance, and Gauss-Markov thrust process noise. |
| **Vehicle geometry** | Thruster force/torque influence matrix, cosine loss, and centre of mass and inertia versus tank mass with analytic derivatives. |
| **Thrust allocation** | Bounded least squares (0 ≤ T ≤ T_max, Lawson–Hanson / Stark–Parker) with an optional minimum-effort balance term, dependency-free, checked by brute force and KKT certificates. |
| **Optimization-ready engine** | Throttle as a control with ∂a/∂u, smooth monotone (PCHIP) throttle characteristics, deep-throttle gap enforced. |
| **Pulse-mode thrusters** | Measured impulse-bit tables, minimum impulse bit, fewest-pulse plans for an impulse or Δv with reported quantization residual (nearest or no-overshoot), and propagator-ready pulse-train artifacts. |
| **Impulsive equivalents** | Free-space Δv, effective exhaust velocity and impulse placement time t* of any tabulated burn (Gauss–Legendre with error estimate). |
| **Verified sensitivities** | Variational equations with saltation at burn events reproduce finite differences of whole propagations; `time_event_jump` supplies the jump. |
| **Multi-burn timelines** | Several bridges on one clock with overlap checks, stage/fairing jettison events, summed concurrent burns and one verified document. |
| **Robust by construction** | Fail-closed input validation at every boundary (fuzz-tested), canonical numbers for stable hashes, and reopenable model documents for bridges, electric models and timelines. |
| **Partials** | Acceleration and mass-rate partials with respect to mass, direction, ambient pressure, available power, operating conditions and bounded estimator drivers, all checked against finite differences. |

### Propulsion physics and studies

- Isentropic, normal/oblique shock, Fanno and Rayleigh relations.
- Ideal and off-design C-D nozzles, internal shocks, altitude sweeps and a
  separation screen.
- NASA Glenn species data, frozen mixtures, an external NASA CEA adapter and an
  experimental local Gibbs equilibrium solver.
- Conical, Hermite, Rao/TOP and MOC contours, exported as CSV/SVG/DXF/STL.
- Loss budgets, Bartz heat flux, regenerative cooling balance and two-phase
  bands.
- Tsiolkovsky sizing, multi-stage optimization and pressure-fed,
  gas-generator and staged-combustion cycle balances.
- Constant and profiled finite burns with explicit propellant-stream accounting
  and system Isp from total tank drain.
- Pressure-fed blowdown, engine clusters with first-order response, and
  imported test traces.
- Seeded uncertainty ensembles, discrete failure scenarios and hierarchical
  ensembles with exact joint weights.
- Hot-fire random-effects calibration, paired covariance calibration and
  model-form / qualification-to-flight discrepancy.
- Versioned workspaces, comparisons and JSON/CSV/SVG/HTML reports.

## Installation

```bash
git clone https://github.com/ayberkdt/Rocket-Propulsion-Lab.git
cd Rocket-Propulsion-Lab
python -m pip install -e .            # no runtime dependencies
python -m pip install -e ".[dev]"     # optional: pytest and ruff
```

Without installing, set `PYTHONPATH=src` (the PowerShell scripts in `scripts/`
do this for you).

## Quick start

### Browser workbench and JSON API

```bash
rocket-propulsion            # or: .\scripts\run.ps1 on Windows
```

Then open <http://127.0.0.1:8765>. The same calculations are available from
`POST /api/v1/...` (see [docs/API.md](docs/API.md)) and from the
`rocket-propulsion burn simulate` CLI.

### Engine performance in Python

```python
from rocket_propulsion.compressible import calculate_off_design_nozzle

nozzle = calculate_off_design_nozzle(
    chamber_pressure_pa=5.0e6, chamber_temperature_k=3500.0,
    throat_area_m2=0.01, area_ratio=25.0, ambient_pressure_pa=101_325.0,
    gamma=1.22, gas_constant_j_kg_k=355.0,
)
print(nozzle.regime, round(nozzle.thrust_n), "N", round(nozzle.specific_impulse_s), "s")
```

More in [`examples/`](examples/): blowdown, calibration, ensembles, feed
systems, electric propulsion, navigation uncertainty, vehicle geometry and
thrust allocation, pulse-mode trim maneuvers, and a two-burn timeline with
jettison checked against the rocket equation.

## Using it from an orbit propagator

A tabulated burn becomes a force model with one call per integration stage:

```python
from rocket_propulsion.propulsion.burns import (
    AdditionalStateConstraint, PropagatorPropulsionBridge, StreamStateBinding,
    TabulatedBurnArtifact, TabulatedBurnSegment,
)

flows = (("oxidizer", 0.24), ("fuel", 0.06))
artifact = TabulatedBurnArtifact(
    segments=(TabulatedBurnSegment(0.0, 120.0, 1000.0, 1000.0, 1000.0, 1000.0,
                                   0.3, 0.3, flows, flows, "steady"),),
    source_id="engine-acceptance-test-42", source_sha256="0" * 64,
    reference_ids=("NASA-20040000363",),
)
bridge = PropagatorPropulsionBridge(
    artifact,
    stream_bindings=(StreamStateBinding("oxidizer", "m_ox"),
                     StreamStateBinding("fuel", "m_fuel")),
    state_constraints=(AdditionalStateConstraint("m_ox", 2.0),),
    protected_dry_mass_kg=400.0,
    nozzle_exit_area_m2=0.004,
)

# Inside the propagator's right-hand side:
value = bridge.evaluate(
    relative_time_s=30.0, vehicle_mass_kg=520.0, direction=(0.0, 0.0, 1.0),
    additional_states={"m_ox": 60.0, "m_fuel": 18.0},
    ambient_pressure_pa=2_000.0, apply_inventory_limits=False,
)
print(value.acceleration_m_s2, value.mass_derivative_kg_s)
print([event.event_id for event in bridge.event_surfaces()])   # register once
```

The consumer registers the roots, lands on time events, switches the burn off
at inventory roots, and maps the analytic partials into its variational
equations. The full protocol is in the
**[astrodynamics integration guide](docs/ASTRODYNAMICS_INTEGRATION.md)**.
The tests run the same protocol with an event-aware RK4 loop that reproduces
the rocket equation, analytic cut-off times and mass closure.

## Project layout

```text
src/rocket_propulsion/
├── api/              HTTP transport, JSON API and static-file server
├── compressible/     Isentropic, shock, duct and nozzle flow
├── core/             Validation (strict input coercion), errors, units, solvers
├── geometry/         Nozzle contours, MOC and CAD export
├── integrations/     Optional, fail-closed Sidera and Lunaris adapters
├── propulsion/       Performance, losses, cycles, staging, propellants
│   └── burns/        Finite burns and the propagation interface:
│       ├── propagator.py          transient burn bridge (events, partials, ambient pressure)
│       ├── performance_surface.py N-D operating maps and kink/domain roots
│       ├── feed_system.py         regulated pressurant and feed ODEs
│       ├── feed_line.py           dynamic liquid line (inertance/compliance)
│       ├── electric.py            power-limited electric propulsion
│       ├── execution_errors.py    priors, Gates errors, Gauss-Markov noise
│       ├── vehicle_geometry.py    thruster force/torque, allocation, mass properties
│       ├── timeline.py            multi-burn schedules and mass jettison
│       ├── throttleable.py        throttle-commanded engine, monotone cubic curves
│       ├── burn_metrics.py        free-space delta-v and impulsive-equivalent time
│       ├── tabulated.py           hash-verified thrust/flow exchange artifact
│       └── ...                    sizing, profiles, blowdown, calibration, ensembles
├── studies/          Workspaces, comparison, sampling and reports
├── thermochemistry/  NASA species data, mixtures and equilibrium
├── thermodynamics/   Ideal-gas, polytropic and heat-transfer relations
└── web/              Browser workbench
tests/                Unit, API and end-to-end integration tests
examples/             Runnable scripts
docs/                 Model, contract and validation documentation
```

Every package has its own README, and public functions document equations,
units, validity limits and references in their docstrings.

## Documentation

| Topic | Documents |
|---|---|
| Start here | [User guide](docs/USER_GUIDE.md) · [API](docs/API.md) · [Architecture](docs/ARCHITECTURE.md) |
| Astrodynamics interface | [Integration guide](docs/ASTRODYNAMICS_INTEGRATION.md) · [Propagator bridge](docs/PROPAGATOR_PROPULSION_BRIDGE.md) · [Operating surfaces](docs/OPERATING_CONDITION_SURFACES.md) · [Electric propulsion](docs/ELECTRIC_PROPULSION.md) · [Execution errors](docs/EXECUTION_ERRORS.md) · [Vehicle geometry](docs/VEHICLE_GEOMETRY.md) · [Timelines](docs/PROPULSION_TIMELINE.md) · [Optimization & sensitivities](docs/TRAJECTORY_OPTIMIZATION.md) · [Lunaris](docs/LUNARIS_INTEGRATION.md) · [Pulse mode](docs/PULSE_MODE.md) |
| Feed systems | [Regulated feed](docs/REGULATED_FEED_DYNAMICS.md) · [Dynamic feed line](docs/DYNAMIC_FEED_LINE.md) · [Blowdown](docs/BLOWDOWN_MODEL.md) |
| Uncertainty | [Propulsion uncertainty](docs/PROPULSION_UNCERTAINTY.md) · [Scenarios](docs/PROPULSION_SCENARIOS.md) · [Hierarchical ensemble](docs/HIERARCHICAL_PROPULSION_ENSEMBLE.md) · [Hot-fire calibration](docs/HOT_FIRE_CALIBRATION.md) · [Model discrepancy](docs/MODEL_DISCREPANCY_TRANSFER.md) |
| Exchange contracts | [Tabulated burn](docs/TABULATED_BURN_CONTRACT.md) · [Sidera plan](docs/SIDERA_BURN_INTEGRATION_PLAN.md) · [Sidera ensemble hand-off](docs/SIDERA_ENSEMBLE_HANDOFF.md) |
| Evidence and status | [Validation report](docs/VALIDATION_REPORT.md) · [Tolerances](docs/TOLERANCES.md) · [Gap analysis](docs/GAP_ANALYSIS.md) · [References](docs/REFERENCES.md) · [Roadmap](docs/ROADMAP.md) · [Release notes](docs/RELEASE_NOTES.md) |

## Testing and quality

```bash
python -m unittest discover -s tests      # or: pytest, or .\scripts\test.ps1
python -m ruff check .
```

CI runs ruff, the full test suite and a package build on Python 3.11 and 3.12.

The test suite includes:
- analytic and reference-data checks for the flow and thermochemistry modules;
- central finite-difference checks of every analytic partial;
- event-aware RK4 integrations that reproduce the rocket equation, closed-form
  positions, analytic reserve and power-switch times, and mass closure;
- Monte Carlo checks of the covariance models;
- a seeded fuzz suite that sends malformed values to every public
  propagation-facing boundary and requires a typed, fail-closed error;
- seeded property tests over random profiles, throttle tables, geometries and
  allocation problems (bounded optimality against brute-force enumeration);
- tamper tests on every persisted document;
- variational-equation integrations (with saltation) checked against finite
  differences of whole propagations;
- a regression test for every finding of a systematic code review;
- an end-to-end Lunaris propagation (when Lunaris is installed) comparing
  finite and impulsive paths;
- an automated gate that every public propulsion-transient type and function
  cites its references.

## Validation status and limits

- **Release gate still open:** the nine-point external NASA CEA comparison
  (LOX/LH2, LOX/RP-1, LOX/CH4) is still open, and the local Gibbs equilibrium
  solver is labelled experimental. See
  [docs/VALIDATION_REPORT.md](docs/VALIDATION_REPORT.md) and
  [docs/DEVELOPMENT_PLAN.md](docs/DEVELOPMENT_PLAN.md).
- **No hardware data:** no engine, thruster or throttle data are bundled.
  Performance maps, throttle tables and test traces are user-supplied evidence,
  bound to their source by hash.
- **Not modelled:** finite-rate chemistry, particle-resolved two-phase flow,
  water hammer / acoustic feed modes, slosh, and nozzle separation inside the
  propagation bridge.
- **Responsibility of the consuming framework:** burn targeting, steering and
  attitude dynamics. `PropulsionTimeline` sequences burns and jettisons but does
  not decide when they happen.

The detailed fidelity boundary is in [docs/GAP_ANALYSIS.md](docs/GAP_ANALYSIS.md).

## Conventions

- SI units throughout; angles in the flow/geometry APIs are in degrees.
- Names carry their unit (`thrust_n`, `ambient_pressure_pa`, `time_s`).
- Invalid inputs fail closed with `DomainError` / `InputError`. Nothing is
  silently extrapolated, clamped or approximated; an approximation must be
  requested explicitly and carries evidence.
- Every exchange artifact is canonical JSON with a SHA-256 content hash.

## Licence

MIT. See [LICENSE](LICENSE).
