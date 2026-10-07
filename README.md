# Rocket Propulsion Lab

Rocket propulsion studies for compressible aerodynamics and thermodynamics in a
small, dependency-free Python application. The current release contains:

- isentropic perfect-gas relations with forward and inverse conversions;
- normal-shock relations;
- ideal-gas thermodynamic state calculations;
- NASA Glenn temperature-dependent properties for ten combustion species;
- frozen ideal-gas mixtures with mole/mass fraction conversion and inverse
  enthalpy/entropy temperature solutions;
- static/stagnation state and compression/expansion connections;
- Fanno and Rayleigh constant-area duct-flow relations;
- weak/strong oblique-shock branches with detachment limits;
- ideal choked converging-diverging nozzle performance;
- back-pressure nozzle regimes, internal normal shocks, altitude sweeps, and a
  source-labelled empirical separation screen;
- conical, preliminary Hermite, Rao/TOP, and MOC nozzle contours with
  station-by-station flow solutions and CSV/SVG/DXF/STL export;
- external NASA CEA adaptation plus an explicitly experimental local Gibbs
  equilibrium provider and O/F performance sweeps;
- switchable performance loss budgets, Bartz wall heat flux, experimental
  regenerative channel balance, and condensed-phase performance bands;
- isentropic, nozzle, P-v, and T-s engineering plots;
- Tsiolkovsky forward/inverse mission calculations and rocket thrust;
- multi-stage mass optimization and pressure-fed, gas-generator, and
  staged-combustion component energy balances;
- versioned `.rplab.json` workspaces, case comparison, seeded Monte Carlo/Latin
  hypercube uncertainty, sensitivity ranking, and JSON/CSV/SVG/HTML reports;
- standalone L0 constant and L1 piecewise-linear burn analysis with explicit
  tank streams, system-Isp closure, browser/workspace/report flows, and
  fail-closed exact or explicitly reduced Sidera hand-off;
- isolated transient-burn development for eventual Sidera integration:
  first-order engine response, restart/cooldown constraints, engine clusters,
  imported solid/test curves, and a hashed discontinuity-safe tabulated
  thrust/mass-flow contract;
- reference-backed preliminary pressure-fed blowdown with polytropic tank
  decay, variable performance, measured pressure maps, convergence evidence,
  and no Sidera source modification;
- sourced propulsion uncertainty ensembles with bounded distributions,
  aleatory/epistemic separation, optional correlated sampling, convergence
  diagnostics and one reusable thrust/mass-flow artifact per realization;
- explicit discrete propulsion scenarios for engine-out/common-cause loss,
  ignition delay, early cutoff and performance degradation, with exhaustive
  probabilities and one reusable artifact per outcome;
- hierarchical propulsion ensembles that combine those discrete outcomes with
  scenario-conditional continuous variation, exact joint weights, independent
  replicate diagnostics, within/between-scenario variance attribution, a
  complete artifact for every realization, and integrity-checked JSON
  persistence;
- an isolated weighted Sidera hand-off contract that binds every realization
  to execution timing/direction while preserving axial force, independent total
  tank drain, zero-thrust probability and nested provenance;
- hot-fire random-effects calibration that separates motor-to-motor,
  burn-to-burn and measurement uncertainty, retains solver evidence and turns
  test results into hash-linked residual ensemble parameters;
- paired hot-fire covariance calibration for thrust/Isp/time scales, including
  full measurement-covariance removal, recorded correlation regularization,
  ensemble-ready copulas, and recomputed versioned persistence;
- model-form and qualification-to-flight discrepancy calibration with explicit
  application domains, predictive epistemic variance, evidence-reuse refusal,
  and layered aleatory/epistemic ensemble components;
- a propagation-ready transient bridge with simultaneous force, total/tank
  mass derivatives, exact discontinuity and reserve event surfaces, bounded
  estimator drivers, analytic Jacobians, and a hash-bound integration manifest;
- multidimensional state-dependent performance surfaces with complete-grid
  validation, analytic condition gradients, no extrapolation and propagator
  events at every declared operating-domain boundary;
- regulated feed/pressurant ODEs with finite bottle inventory, thermal ullage,
  choked valve flow, regulator dynamics, feed-line loss, protective events and
  implicitly closed engine-flow/pressure Jacobians;
- a two-state dynamic liquid feed line with hydraulic inertance, compliance,
  nonlinear damping, manifold-pressure coupling, protective roots, mass
  closure and analytic composed Jacobians for orbit propagation;
- a screening catalog for common liquid, solid, and hybrid propellant systems;
- a responsive browser interface backed by a versioned JSON API.

## Quick start

```powershell
.\scripts\run.ps1
```

Then open <http://127.0.0.1:8765>. The launcher searches for an installed Python
3.11+ and falls back to the Python runtime bundled with Codex Desktop.

## Test

```powershell
.\scripts\test.ps1
```

No package installation is required. Both scripts set `PYTHONPATH` to `src`.

## Project map

```text
src/rocket_propulsion/
├── api/              HTTP transport and static-file server
├── compressible/     Isentropic and shock-flow physics
├── core/             Shared validation and numerical solvers
├── geometry/         Parametric nozzle contours and flow stations
├── integrations/     Optional, fail-closed Sidera hand-off adapters
├── propulsion/       Rocket equations, burn analysis, thrust, and propellants
├── studies/          Workspaces, comparison, sampling, and reports
├── thermochemistry/  NASA species data and frozen ideal-mixture properties
├── thermodynamics/   Calorically perfect-gas relations
└── web/              Browser UI
tests/                Unit and API-level tests
docs/                 Architecture and roadmap
```

Every package owns a focused README. Public functions include equations,
parameters, return values, and domain restrictions in their docstrings.

The first standalone finite-burn slice is available through the Python API,
browser Burn workbench, `POST /api/v1/engine-operating-point`,
`POST /api/v1/burn/simulate`, and the `rocket-propulsion burn simulate`
command. It performs analytic constant-thrust/constant-flow inventory,
impulse, and ideal one-dimensional delta-v accounting without depending on
Sidera. L1 ramp/steady/shutdown schedules use exact analytic segment integrals;
eligible constant results can be handed off exactly, while varying profiles
require an explicit, evidence-carrying equivalent reduction. See
[`examples/constant_burn_input.json`](examples/constant_burn_input.json) and
[`examples/profiled_burn_input.json`](examples/profiled_burn_input.json).

## Scope and conventions

- SI units are used for dimensional quantities.
- Angles are entered and returned in degrees.
- `gamma` must be greater than one; Mach numbers must be positive.
- This is an engineering study tool, not flight-certified analysis software.

The current fidelity boundary is explicit in
[docs/GAP_ANALYSIS.md](docs/GAP_ANALYSIS.md); data provenance and model
references are recorded in [docs/REFERENCES.md](docs/REFERENCES.md). Use the
[user guide](docs/USER_GUIDE.md), [API guide](docs/API.md), and
[validation report](docs/VALIDATION_REPORT.md) for reproducible engineering
work. The remaining external CEA release gate is tracked in
[docs/DEVELOPMENT_PLAN.md](docs/DEVELOPMENT_PLAN.md). The planned standalone
burn-analysis library and its optional hand-off to the existing Sidera
astrodynamics framework are specified in
[docs/SIDERA_BURN_INTEGRATION_PLAN.md](docs/SIDERA_BURN_INTEGRATION_PLAN.md).
The isolated transient exchange format and future consumer rules are in
[docs/TABULATED_BURN_CONTRACT.md](docs/TABULATED_BURN_CONTRACT.md).
The pressure-fed equations, applicability limits and primary bibliography are
in [docs/BLOWDOWN_MODEL.md](docs/BLOWDOWN_MODEL.md).
The uncertainty/ensemble contract and its Sidera ownership boundary are in
[docs/PROPULSION_UNCERTAINTY.md](docs/PROPULSION_UNCERTAINTY.md).
The probability-labelled failure/execution scenario contract is in
[docs/PROPULSION_SCENARIOS.md](docs/PROPULSION_SCENARIOS.md).
Their implemented conditional composition and future Sidera ensemble contract
are in
[docs/HIERARCHICAL_PROPULSION_ENSEMBLE.md](docs/HIERARCHICAL_PROPULSION_ENSEMBLE.md).
The producer-side weighted transient hand-off and exact capability requirements
for future Sidera consumption are in
[docs/SIDERA_ENSEMBLE_HANDOFF.md](docs/SIDERA_ENSEMBLE_HANDOFF.md).
The repeated-test calibration equations, assumptions and NIST/NASA references
are in [docs/HOT_FIRE_CALIBRATION.md](docs/HOT_FIRE_CALIBRATION.md).
Model validation discrepancy, qualification-to-flight transfer and their
layered ensemble workflow are in
[docs/MODEL_DISCREPANCY_TRANSFER.md](docs/MODEL_DISCREPANCY_TRANSFER.md).

