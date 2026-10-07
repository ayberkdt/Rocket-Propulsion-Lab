# Standalone burn analysis and Sidera integration plan

## 1. Decision

Rocket Propulsion Lab will own propulsion performance during a burn. Sidera
will continue to own astrodynamics, reference frames, guidance direction,
events, and trajectory integration.

No orbit propagator, gravity model, RIC/VNB transform, attitude law, or
trajectory optimizer will be duplicated in this project.

The burn subsystem is a first-class Rocket Propulsion Lab capability. It must
remain fully usable when Sidera is not installed. Sidera export is one optional
adapter over the public burn result, not the subsystem's internal data model or
primary execution path.

The integration boundary will be:

```text
Rocket Propulsion Lab                         Sidera
---------------------                         ------
engine/nozzle operating point                 epoch and maneuver ordering
thrust and propellant-flow history    --->    force integration
start/stop and throttle transients            inertial/RIC/VNB direction
mixture-ratio inventory accounting            gravity and perturbations
burn sizing and constraint checks             event handling
equivalent constant-burn reduction             trajectory and orbital results
```

This plan was prepared against the local `D:\sidera` working tree on
2026-10-05, at base commit `6dcb8f63eadbadd25183c49956c2c77ae1d447a9`
(`0.1.0rc3.dev0`). The working tree contains user changes, so the integration
tests must record the actual Sidera version and burn-schema capability rather
than relying only on the commit identifier.

### Implementation status — 2026-10-06

The first executable vertical slice is now present:

- immutable L0 burn, tank, stream, operating-point, profile, summary and
  evidence records;
- analytic duration, propellant, impulse and ideal-delta-v target solutions;
- reserve, dry-mass and maximum-duration protective cutoffs;
- explicit total-tank-flow/system-Isp closure, including dump streams;
- deterministic `rocket_propulsion_burn_v1` serialization, re-open,
  hash verification and physical reproduction check;
- Python, HTTP and CLI execution plus standalone HTML report;
- artifact-only and optional native exact-constant Sidera export;
- capability inspection through Sidera's public propagation facade;
- a cross-project native contract test against the local Sidera environment;
- closed bipropellant and monopropellant provider convenience constructors;
- browser Burn workbench with tank/mass results, JSON/report and Sidera export;
- versioned workspace save/re-open of the complete canonical burn artifact;
- L1 piecewise-linear throttle/start/stop schedules with analytic exposure,
  impulse, tank drain, delta-v and force-centroid integrals;
- reproducible `rocket_propulsion_burn_v2` persistence and Python/API/HTTP/CLI
  parity, including browser thrust/mass histories;
- fail-closed transient Sidera export and explicit equivalent-constant
  reduction with impulse/mass identities and approximation evidence.
- isolated analytic first-order response, ignition delay, restart/cooldown and
  minimum-pulse validation, engine-cluster composition, imported trace replay,
  and a discontinuity-safe `rocket_propulsion_tabulated_burn_v1` prototype;
- convergence evidence between exact cluster response integrals and the
  piecewise-linear artifact intended for a future Sidera propagator adapter.
- sourced continuous propulsion UQ, explicit discrete/common-cause outcomes,
  and their hierarchical composition into exact joint probabilities;
- one complete tabulated artifact per joint realization, including deterministic
  provenance and an independently seeded conditional sampling audit.
- a target-neutral, hashed `rocket_propulsion_hierarchical_ensemble_v1`
  document whose reader reproduces joint weights, propulsion statistics,
  within/between-scenario variance and nested artifact integrity before any
  future Sidera propagation.
- an integrity-checked `rocket_propulsion_sidera_ensemble_handoff_v1` adapter
  artifact that binds each realization to start time, frame and direction,
  explicitly selects axial thrust and total tank drain, and fails capability
  negotiation until Sidera exposes exact tabulated/weighted support.
- paired hot-fire covariance plus independent model-form and
  qualification-to-flight discrepancy layers; qualified scale components keep
  aleatory and epistemic inputs separate while producing the same complete
  weighted artifacts for future Sidera propagation.

Still open before the broader B1 evidence target is complete: expanded
reference fixtures/property grids. The B2 analytic and exchange foundation is
implemented in Rocket Propulsion Lab; its API/CLI/UI product surfaces and the
future Sidera-side consumer remain intentionally open. The preliminary
pressure-fed blowdown provider has started B3; coupled regulator/feed-line and
thermal-pressurant dynamics remain open.

### 1.1 Product modes

The same burn core supports four independent workflows:

1. **Standalone performance analysis** — size or simulate a burn and report
   thrust, total impulse, tank drain, mixture consumption, effective Isp,
   ideal one-dimensional velocity increment and constraints.
2. **Engine and propulsion-system trade study** — compare propellant pairs,
   cycles, engine counts, throttle schedules, tank reserves and uncertainty
   without requiring any trajectory software.
3. **Profile/artifact generation** — export a versioned thrust and mass-flow
   history to CSV/JSON for test planning, controls work or another simulator.
4. **Sidera hand-off** — turn a compatible result into Sidera's maneuver
   contract and let Sidera compute the orbit.
5. **Propulsion risk ensemble** — combine sourced discrete outcomes with
   scenario-conditional continuous variation, preserving one joint weight and
   one full artifact per realization before any trajectory propagation.

The standalone modes are release-blocking. An integration failure or absent
Sidera installation must never disable burn sizing, simulation or reporting.

### 1.2 Required fidelity levels

Every run selects and reports one fidelity level; the library never upgrades
or downgrades it silently.

| Level | Model | Intended use |
| --- | --- | --- |
| L0 | Constant delivered thrust and system Isp | quick sizing and exact current-Sidera export |
| L1 | Piecewise throttle/start/stop profile | total impulse, transient losses and short-burn studies |
| L2 | State-dependent provider such as blow-down or power-limited thrust | spacecraft propulsion-system trades |
| L3 | Imported test/qualification thrust and flow histories | measured-profile replay and correlation |

L0 and L1 are dependency-free. L2 providers may use existing thermochemistry,
nozzle, cycle and loss models. L3 treats imported data as evidence with units,
source identity and interpolation policy; it is not silently smoothed.

### 1.3 Success criteria

The subsystem is successful only when all of the following are true:

- `rocket_propulsion.propulsion.burns` can be imported and used without Sidera;
- a complete burn can be evaluated from Python without the HTTP server or UI;
- fuel, oxidizer, pressurant/dump streams and vehicle mass close numerically;
- constant and variable performance use unambiguous delta-v definitions;
- model provenance and approximation warnings survive JSON/CSV/report export;
- a supported constant result maps exactly to Sidera;
- an unsupported profile fails closed or requires explicit approximation;
- existing steady thrust/nozzle APIs remain backward compatible.

## 2. What Sidera already provides

Sidera's current public maneuver surface is sufficient for a first constant
finite-burn integration:

- `FiniteBurn(t_start_s, duration_s, thrust_n, isp_s, direction, frame,
  throttle)`;
- `ImpulsiveManeuver` and an immutable, ordered `ManeuverPlan`;
- inertial, RIC and VNB finite-burn directions, rebuilt from the current state
  during integration for local frames;
- seven-state propagation with mass at `y[6]`;
- ideal-rocket mass flow, `mdot = throttle * thrust / (g0 * Isp)`;
- segmented propagation at ignition and cutoff boundaries;
- terminal-event precedence, non-overlap validation and deterministic plan
  hashing;
- JSON round-trip through `sidera_maneuver_plan_v1`;
- CCSDS OPM maneuver conversion;
- body-centred propagation that obtains initial mass from
  `SpacecraftProps.mass_kg` when a maneuver changes mass.

The relevant Sidera implementation is in:

- `src/sidera/core/propagation/plans.py`;
- `src/sidera/core/propagation/propagator.py`;
- `src/sidera/core/body_engine/orbits.py`;
- `src/sidera/core/propagation/ccsds_maneuvers.py`;
- `tests/core/test_maneuver_propagation.py`.

### Current integration ceiling

Sidera's current `FiniteBurn` is intentionally simple. It does not yet carry:

- time-varying thrust, specific impulse or propellant flow;
- independent thrust and total-system propellant flow;
- ignition/shutdown ramps, pulsing or duty-cycle histories;
- separate fuel and oxidizer inventories;
- body-frame finite-burn steering;
- maneuver checkpoint/resume or maneuver telemetry;
- a maneuver-aware two-body comparison baseline.

Therefore a constant operating point can be exported exactly today. A
transient or blow-down burn must not be silently flattened into a constant
burn and presented as exact.

## 3. Existing Rocket Propulsion Lab gap

The current project can calculate a steady thrust and optionally multiply it
by a burn duration. It does not yet model a burn as a first-class engineering
object. In particular, it lacks:

- ignition-to-cutoff state and constraints;
- initial/final vehicle mass and propellant inventory closure;
- fuel/oxidizer flow split from mixture ratio;
- open-cycle dumped flow in the vehicle mass depletion rate;
- start/stop transients and throttle schedules;
- equivalent delivered performance over a varying profile;
- a versioned burn artifact suitable for another program;
- compatibility checks against Sidera's installed burn contract.

The existing `calculate_thrust`, loss budget, engine-cycle, equilibrium and
nozzle models are inputs to the new burn layer; they are not replaced.

## 4. Ownership rules

### Rocket Propulsion Lab owns

- delivered thrust magnitude and its uncertainty;
- chamber/nozzle operating point and vacuum performance;
- total vehicle propellant flow, including dumped open-cycle flow;
- oxidizer/fuel consumption and residual margins;
- engine count, cant efficiency, throttle limits and duty cycle;
- ignition delay, ramp-up, steady operation and ramp-down;
- burn duration or target-delta-v sizing;
- total impulse, equivalent Isp and ideal-rocket delta-v bookkeeping;
- conversion to a stable, versioned propulsion burn artifact;
- the optional Sidera adapter.

### Sidera owns

- absolute epoch and time-system interpretation;
- burn direction, reference frame and attitude/guidance law;
- position, velocity and orbital-element evolution;
- gravity, atmosphere and other perturbing accelerations;
- gravity, steering and finite-duration trajectory losses;
- terminal events, impact/SOI logic and numerical integration;
- CCSDS orbit-message semantics.

The propulsion core will never import Sidera. Only the adapter package may
import both domains.

## 5. Target package structure

```text
src/rocket_propulsion/
├── propulsion/
│   └── burns/
│       ├── __init__.py
│       ├── models.py          immutable burn inputs, points and results
│       ├── operating_point.py engine/nozzle/cycle performance assembly
│       ├── sizing.py          duration, propellant and target-delta-v solves
│       ├── profiles.py        throttle and start/stop histories
│       ├── inventory.py       fuel/oxidizer/residual mass accounting
│       ├── integration.py     propulsion-state integration and event roots
│       ├── reduction.py       equivalent constant-burn calculation
│       ├── uncertainty.py     burn-specific uncertain-input mapping
│       └── README.md          equations, validity and worked examples
├── integrations/
│   ├── __init__.py
│   ├── contracts.py           target-neutral export capability protocol
│   └── sidera.py              optional adapter; no orbital equations
├── studies/
│   └── burn_studies.py        sweeps, comparisons and seeded uncertainty
└── api/
    └── routes.py              versioned burn and export endpoints
```

Every new module receives its own docstrings and the `burns/` package receives
a README defining equations, units, fidelity and examples.

### 5.1 Dependency direction

```text
web -> api -> integrations -> propulsion.burns -> existing propulsion models
                    |                 |
                    |                 +-> thermochemistry / compressible /
                    |                     thermodynamics through public APIs
                    +-> optional Sidera public API

studies -> propulsion.burns
reports -> burn result/artifact
```

Rules:

- `propulsion.burns` imports neither `integrations` nor Sidera;
- providers depend on existing public domain APIs, not HTTP routes;
- the browser never calculates thrust, flow, impulse or inventory;
- adapters consume immutable public results and cannot mutate a burn;
- Sidera is not added to base dependencies; if a packaged version becomes
  available it is exposed through an optional `sidera` extra;
- production code never discovers Sidera through the local `D:\sidera` path;
  the local path is allowed only in an explicitly configured integration test.

### 5.2 Provider-oriented design

The burn solver receives performance through a small protocol rather than
containing a large engine-type conditional:

```python
class PropulsionProvider(Protocol):
    provider_id: str
    provider_version: str

    def evaluate(
        self,
        *,
        elapsed_s: float,
        command: BurnCommand,
        state: PropulsionState,
        environment: BurnEnvironment,
    ) -> EngineOperatingPoint: ...

    def discontinuities_s(self, definition: BurnDefinition) -> tuple[float, ...]: ...
```

Initial providers:

- `ConstantPerformanceProvider` for L0 sizing and Sidera parity;
- `ProfiledPerformanceProvider` for prescribed L1 schedules;
- `RocketModelProvider` that assembles equilibrium/nozzle/loss/cycle results;
- `ImportedTraceProvider` for test data with declared interpolation;
- `PressureFedBlowdownProvider` as the first L2 state-dependent model.

The protocol returns actual delivered thrust and explicit stream flows. It
does not return orbital acceleration, because division by vehicle mass and
force-frame transformation belong to the consuming dynamics program.

## 6. Domain contracts

### 6.1 `EngineOperatingPoint`

One delivered, steady operating point:

- commanded and realized throttle;
- active engine count;
- ideal and delivered thrust [N];
- chamber-only and system-effective Isp [s];
- chamber flow and total tank-drain flow [kg/s];
- oxidizer and fuel flow [kg/s];
- mixture ratio;
- chamber pressure and ambient pressure [Pa];
- loss-budget/cycle/equilibrium model identifiers;
- warnings, validity ranges and source metadata.

It also exposes `streams`, a tuple of `PropellantFlow` records. Each stream
contains a stable species/tank key, mass flow, destination (`chamber`,
`gas_generator_dump`, `coolant_dump`, `vent` or `other`) and whether its mass
contributes to delivered nozzle momentum. This prevents a single ambiguous
`mass_flow_kg_s` from representing both chamber flow and tank depletion.

For gas-generator engines, Sidera-facing mass depletion must use total tank
drain. Consequently the exported effective Isp is

```text
Isp_system = delivered_thrust / (g0 * total_tank_drain)
```

not the chamber-only nozzle Isp. This prevents the trajectory from silently
under-consuming propellant.

The following identities are checked at every operating point:

```text
mdot_total     = sum(all tank-depleting stream flows)
mdot_chamber   = sum(streams delivered to the main chamber/nozzle)
mdot_oxidizer  = sum(oxidizer-family streams)
mdot_fuel      = sum(fuel-family streams)
Isp_system     = F_delivered / (g0 mdot_total)
```

The system Isp identity may be undefined only at a true zero-thrust/zero-flow
point. A state with thrust but no positive mass flow is rejected unless the
provider explicitly declares a non-propellant momentum source.

### 6.2 `BurnDefinition`

The propulsion-owned burn request contains:

- initial wet mass and protected dry-mass floor [kg];
- available fuel and oxidizer [kg];
- engine/operating-point selection;
- throttle schedule;
- ignition delay and ramp durations;
- exactly one primary termination target: duration, propellant expenditure,
  total impulse, or ideal delta-v;
- reserve policy and maximum burn duration;
- optional engine count/cant efficiency.

It deliberately contains no position, velocity, orbital frame or direction.

`BurnDefinition` references immutable equipment rather than embedding copied
performance numbers:

- `EngineDefinition`: name, family, rated limits, restart/minimum-pulse data
  and provider selection;
- `TankDefinition`: propellant identity, usable load, trapped residual,
  pressurant and pressure model identifier;
- `EngineClusterDefinition`: engines, enabled set, body-frame mount direction,
  cant efficiency and optional moment-arm metadata;
- `ReservePolicy`: absolute/relative reserve by tank plus protected dry mass.

Mount vectors are propulsion/vehicle hardware data. The standalone library may
sum body-frame force and torque, but it does not propagate attitude. Sidera's
current three-frame `FiniteBurn` adapter consumes only the scalar effective
thrust and export-time direction selected by the caller.

### 6.3 `BurnCommand`, `BurnEnvironment` and `PropulsionState`

`BurnCommand` is what the propulsion system is asked to do:

- commanded throttle or per-engine throttles;
- enabled engine identifiers;
- optional gimbal commands expressed only in engine/body hardware axes;
- start, stop and pulse commands;
- command timestamp relative to the burn artifact.

`BurnEnvironment` contains only quantities that affect propulsion performance:

- ambient pressure [Pa];
- optional inlet/back pressure [Pa];
- optional available electrical power [W];
- optional declared thermal boundary values;
- environment source and interpolation policy.

It contains no central body, position, velocity or gravity. A constant vacuum
environment is the default for orbital chemical burns; it is always serialized
rather than implied.

`PropulsionState` carries the state integrated by this library:

- remaining mass in each tank [kg];
- total vehicle mass used for ideal velocity-increment bookkeeping [kg];
- tank/pressurant state required by the selected provider;
- engine phase, elapsed on-time, restart count and cooldown timers;
- provider-specific state in a versioned, namespaced record.

Provider state is never an unrestricted dictionary in the public result; it
must be represented by a declared schema/version so saved studies remain
reproducible.

### 6.4 `BurnProfilePoint` and `BurnProfile`

Each sampled point records:

- time since command and time since ignition [s];
- commanded/realized throttle;
- thrust [N], system-effective Isp [s] and total flow [kg/s];
- oxidizer/fuel flow and remaining inventory [kg];
- vehicle mass [kg], chamber pressure [Pa] and phase label;
- any active constraint or warning.

The profile carries a schema version, deterministic content hash, model
versions and units. Sampling is presentation/export resolution; integration
uses analytic segments where available and a declared quadrature tolerance
otherwise.

### 6.5 `BurnSummary`

Required outputs:

- achieved duration and cutoff reason;
- initial/final mass and consumed propellant by component;
- delivered total impulse;
- time-averaged and impulse-weighted thrust;
- system-effective equivalent Isp;
- ideal-rocket delta-v from actual mass ratio;
- thrust centroid time;
- residual mass, impulse and delta-v closure errors;
- uncertainty bands and fidelity warnings;
- whether exact constant-Sidera export is available.

The summary distinguishes quantities that are often incorrectly conflated:

```text
total impulse:                   J = integral(F dt)
consumed propellant:            Delta_m = integral(mdot_total dt)
impulse-weighted/system Isp:    Isp_J = J / (g0 Delta_m)
ideal 1-D velocity increment:   Delta_v_1D = integral(F / m_vehicle dt)
rocket-equivalent Isp:          Isp_dv = Delta_v_1D /
                                         (g0 ln(m0 / mf))
```

For constant system Isp, `Isp_J`, `Isp_dv` and the operating-point system Isp
agree. With varying Isp they need not agree, so all are named explicitly. The
library must not compute a variable-performance delta-v by inserting an
averaged Isp into Tsiolkovsky's constant-Isp form.

### 6.6 Standalone result contract

`BurnResult` is the top-level return value and contains:

- normalized `BurnDefinition` and initial/final `PropulsionState`;
- `BurnProfile` and `BurnSummary`;
- all warnings and binding constraints;
- convergence and integration evidence;
- model/dataset versions and external-source provenance;
- deterministic input and result hashes;
- an `export_capabilities` record that reports exact, approximate or
  unsupported status for each installed adapter.

The result has no Sidera types. It serializes as
`rocket_propulsion_burn_v1`; adapters are pure transformations from that
artifact.

### 6.7 Supported propulsion families

The architecture supports the following without forcing one model onto all
systems:

| Family | Initial implementation | Later fidelity |
| --- | --- | --- |
| Liquid bipropellant | constant/profiled O/F, separate tanks, cycle total flow | regulated/blow-down tanks, valve/feed transients |
| Liquid monopropellant | one propellant stream and imported/constant performance | pressure and catalyst-bed warm-up dependence |
| Solid rocket motor | imported thrust and mass-flow curve | grain regression provider and nozzle erosion |
| Hybrid | prescribed oxidizer flow and O/F history | regression-coupled fuel flow |
| Cold gas | pressure-dependent thrust/flow provider | thermal real-gas blow-down |
| Electric propulsion | explicit thrust/Isp or imported map | available-power and efficiency coupling |

Only liquid constant/profiled chemical propulsion is required for B1. Solid
curve replay and monopropellant blow-down enter B2/B3. Electric propulsion is
supported by the contracts but is not claimed validated until its own provider
and reference matrix exist.

### 6.8 Propulsion-only state equations

The standalone solver integrates propulsion state, not trajectory state. For
each tank or tracked inventory `i`:

```text
d m_i / dt      = -sum(mdot_stream,i)
d m_vehicle/dt  = -sum(all external tank-depleting streams)
d J / dt        = F_delivered
d Delta_v_1D/dt = F_axial / m_vehicle
```

Optional provider states add their own equations, for example an isothermal
ideal-gas blow-down pressure state. They must declare equations, assumptions,
validity ranges and invariants through metadata.

The default velocity-increment channel uses delivered axial thrust after
declared cant/divergence losses. A separate gross-thrust integral is retained
when gross and axial thrust differ. Neither channel includes gravity, drag or
steering losses; those are trajectory effects.

### 6.9 Numerical execution model

1. Validate and normalize all definitions, units and termination targets.
2. Resolve the schedule's exact discontinuities: command, ignition, throttle
   knots, engine events, shutdown and provider events.
3. Split the run at every discontinuity; never integrate across a step.
4. Use closed forms for constant and linear prescribed segments.
5. Use a deterministic adaptive integrator only for state-dependent providers.
6. Locate inventory, mass-floor, impulse and delta-v cutoffs with a bounded
   event root solve inside the active segment.
7. Evaluate both sides of a discontinuity and record pre/post semantics.
8. Recompute independent conservation residuals from the completed result.

The L0/L1 implementation remains dependency-free. A built-in embedded RK pair
may be added for L2 only after comparison against a reference integrator; until
then an L2 provider that cannot supply a closed segment solution is marked
experimental and unavailable rather than evaluated with an ad hoc time step.

Every numerical result reports:

- absolute and relative tolerances;
- accepted/rejected step counts where applicable;
- event time and active cutoff constraint;
- maximum local mass/impulse residual;
- profile-resolution convergence for exported samples;
- whether an analytic or numerical path was used.

### 6.10 Hard invariants and failure policy

The solver fails closed when any of these invariants is violated:

- vehicle and tank masses are finite and never negative;
- vehicle mass equals protected dry mass plus tracked inventories plus declared
  untracked consumables within tolerance;
- no tank crosses its usable reserve;
- total tank drain equals the sum of external stream drains;
- thrust, flow and pressure remain in the provider's validity envelope;
- throttle and engine count remain within hardware limits;
- time is strictly increasing and schedule discontinuities are ordered;
- total impulse and ideal velocity increment are monotone for non-negative
  axial thrust;
- an engine cannot ignite during mandatory cooldown or beyond restart count;
- an imported trace has monotone time, explicit units and no NaN/Inf values.

Warnings are reserved for valid-but-low-fidelity operation. Violations of
conservation, inventories, schema compatibility or required model data are
errors, not warnings.

### 6.11 Termination and feasibility

Exactly one primary target is selected:

- duration;
- total impulse;
- ideal one-dimensional velocity increment;
- total propellant expenditure;
- named tank reserve/depletion.

All other limits are protective constraints. The first reached limit ends the
burn and the result identifies `target_reached`, `constraint_limited`,
`commanded_shutdown`, `provider_invalid` or `numerical_failure`. A target that
cannot be reached with available propellant, thrust limits or maximum duration
returns an infeasible result with the attainable bound; it is never reported
as a successful shortened burn.

Sizing uses bracketed monotonic root solves where the selected provider makes
the target monotonic. Non-monotonic schedules require forward simulation and
explicit optimization in `studies`; the low-level solver does not guess a
branch.

## 7. Sidera adapter contract

### 7.1 Exact constant-burn path

`integrations.sidera.to_finite_burn(...)` will accept:

- a constant `BurnSummary`/operating point;
- `t_start_s`;
- direction and one of `inertial`, `ric`, or `vnb`.

It will return either a native `sidera.core.propagation.FiniteBurn` when Sidera
is importable, or a validated `sidera_maneuver_plan_v1` JSON-compatible object
when artifact-only export is requested.

Mapping:

| Rocket Propulsion Lab | Sidera `FiniteBurn` |
| --- | --- |
| start time | `t_start_s` |
| achieved duration | `duration_s` |
| actual constant delivered thrust | `thrust_n` |
| system-effective Isp | `isp_s` |
| export-only direction | `direction` |
| export-only frame | `frame` |
| unity (throttle is already resolved by the propulsion model) | `throttle` |

The original commanded and realized throttle remain in the propulsion artifact
provenance. Passing actual delivered thrust and the same throttle again would
double-apply throttle in Sidera. A future adapter may preserve Sidera's native
throttle field only when a validated linear full-thrust reference is available.

The adapter checks Sidera's public constructor and schema capability at
runtime and fails with a typed `FeatureUnavailableError` on incompatibility.
It will not guess field meanings or use Sidera private functions.

### 7.2 Varying-profile path

The first implementation will expose two explicit choices:

1. `exact`: refuse export unless thrust and system-effective Isp are constant
   within declared tolerances;
2. `equivalent_constant`: export the impulse/mass-matched reduction with a
   prominent approximation warning and the original profile hash.

Equivalent reduction is useful for short chemical burns but is not declared
trajectory-equivalent. Its report must include variation ranges, centroid
shift and a duration-to-local-orbital-period ratio supplied by Sidera or the
caller.

For exact transient propagation, propose a small Sidera contract extension:

```text
TabulatedFiniteBurn(
    t_start_s,
    time_offsets_s,
    thrust_n,
    mass_flow_kg_s,
    direction,
    frame,
    interpolation="linear" | "previous",
)
```

Sidera would only evaluate this supplied propulsion law inside its existing
segmented thrust wrapper. That is an extension of the force input, not a new
astrodynamics implementation. Independent `mass_flow_kg_s` is required so an
open-cycle engine and future non-ideal propulsion systems do not have to encode
tank drain indirectly through Isp.

Adjacent piecewise `FiniteBurn` objects are not a safe workaround: the current
`ManeuverPlan` rejects a following maneuver whose start is equal to the prior
burn's end, and artificial gaps would change the trajectory.

### 7.3 Compatibility matrix

| Rocket Propulsion Lab result | Current Sidera mapping | Status |
| --- | --- | --- |
| Constant thrust + constant system Isp | one `FiniteBurn` | exact |
| Constant thrust with open-cycle total drain represented by system Isp | one `FiniteBurn` | exact for scalar force/mass |
| Impulsive approximation requested by user | one `ImpulsiveManeuver` | approximate, explicit |
| Start/stop or throttle profile | equivalent constant only | approximate until tabulated contract |
| Solid/imported thrust curve | none | requires tabulated contract |
| Blow-down thrust/Isp/flow history | none | requires tabulated contract |
| Body-frame engine cluster force | no finite body-frame support | unsupported |
| Body torque/gimbal history | no current maneuver torque contract | unsupported |
| Pulsed duty cycle | averaged constant only when explicitly accepted | approximate |
| Electric power-limited profile | none | requires tabulated contract |

`exact` refers only to transferring the propulsion force/mass law represented
by Sidera. It does not claim that either program is flight-qualified or that a
scalar thrust model captures plume, slosh, attitude or feed-system dynamics.

### 7.4 Adapter capability negotiation

The adapter checks capabilities rather than branching only on a version
string:

- importability and distribution version;
- public availability and constructor signature of `FiniteBurn` and
  `ManeuverPlan`;
- accepted maneuver schema identifier;
- supported frames;
- whether independent mass-flow/tabulated burn types exist;
- serialization round-trip and a small numerical parity probe.

The result is a typed `SideraCapabilities` record persisted with exports. A
known version with an unexpected signature fails closed. An unknown version is
not accepted merely because an import succeeds.

### 7.5 Artifact boundary

Two artifacts remain distinct:

1. `rocket_propulsion_burn_v1` — authoritative propulsion result containing
   the full profile, streams, inventories, equations and provenance;
2. `sidera_maneuver_plan_v1` or its future successor — derived trajectory input.

The Sidera artifact includes a sidecar manifest with:

- source burn hash;
- reduction mode and tolerances;
- fields omitted by the target schema;
- Rocket Propulsion Lab and Sidera versions;
- exported thrust, Isp/mass-flow semantics and frame;
- exact/approximate status.

Round-tripping through Sidera must not replace or mutate the authoritative burn
artifact. Trajectory outputs are linked by the source burn hash.

### 7.6 Cross-project boundary test

The integration smoke test performs the following without editing Sidera:

1. build an L0 standalone burn in Rocket Propulsion Lab;
2. export it through the public Sidera adapter;
3. construct/round-trip the public `ManeuverPlan`;
4. propagate Sidera free-space and point-mass reference cases;
5. compare final mass and ideal velocity increment with the standalone result;
6. record Sidera diagnostics, plan hash and package version.

The normal Rocket Propulsion Lab test suite uses a contract fixture and does
not require Sidera. A separate optional integration job runs against supported
Sidera versions or the explicitly configured local checkout.

## 8. Delivery phases

Each phase must leave a usable standalone library. Optional Sidera work cannot
block the standalone release train.

### Phase B0 — contracts, equations and traceability

**Deliver**

- frozen public models for definitions, commands, environments, streams,
  states, profile points, summaries and results;
- `PropulsionProvider` and target-neutral export capability protocols;
- `rocket_propulsion_burn_v1` JSON schema and canonical encoder;
- stable warning/error/cutoff enumerations;
- explicit definitions for total impulse, system Isp, ideal 1-D delta-v and
  rocket-equivalent Isp;
- algebraic conservation/residual utilities;
- Sidera constant-burn contract fixture and reviewed-capability record;
- source/reference table for every equation and correlation.

**Acceptance**

- no Sidera import below `integrations/`;
- no mutable collections in public result objects;
- all dimensional values are SI, finite and documented;
- invalid inventory, dry-mass floor and conflicting targets fail before work;
- serialization is byte-deterministic and rejects NaN/Infinity;
- JSON schema round-trip preserves every public field and content hash;
- import-time dependency test proves the burn package works with Sidera and
  optional scientific packages absent.

**Exit artifact**

- one hand-authored constant bipropellant case and its canonical JSON result,
  used as the compatibility golden for B1.

### Phase B1 — standalone constant-burn release

**Deliver**

- L0 constant provider and closed-form sizing by duration, propellant, impulse
  or ideal 1-D delta-v;
- separate fuel/oxidizer inventories and O/F closure;
- explicit chamber, bypass/dump and total tank-drain streams;
- multi-engine scalar cluster, availability and cant-efficiency handling;
- dry-mass, reserve, maximum-duration and engine-limit cutoffs;
- Python API, `/api/v1/burn/simulate`, `/api/v1/burn/size` and CLI command;
- JSON/CSV/HTML standalone report and minimal Burn UI;
- `.rplab.json` workspace schema migration for burn cases;
- exact current-Sidera constant export as an optional feature.

**Acceptance**

- all four termination modes recover the same constant reference case;
- forward duration-to-delta-v and inverse delta-v-to-duration round-trip within
  `1e-12` relative for well-conditioned cases;
- total tank drain equals the sum of stream and tank drains;
- gas-generator mass depletion uses total, not chamber, flow;
- enabled engine thrust sums correctly and cant loss is independently closed;
- no burn crosses dry mass, per-tank reserve or provider validity limits;
- standalone Python/API/CLI/UI paths return the same canonical result hash;
- Rocket Propulsion Lab and Sidera produce matching final mass and ideal
  delta-v to `1e-10` relative for the parity case;
- Sidera export round-trips and preserves its maneuver-plan hash.

**Recommended version**

- `0.10.0` if the external CEA gate still keeps the project below 1.0;
  otherwise the next stable feature minor, provisionally `1.1.0`.

### Phase B2 — standalone transient/profile release (in progress)

**Deliver**

- ignition delay and named `off`, `starting`, `steady`, `stopping`, `cooldown`
  phases;
- linear/tabulated start and shutdown ramps;
- piecewise-constant and piecewise-linear throttle schedules;
- optional first-order throttle response lag;
- minimum throttle, minimum on/off time, restart count and cooldown rules;
- profile import with units, duplicate-time policy and interpolation choice;
- solid-motor/imported thrust and mass-flow curve replay;
- exact analytic segment integration and convergence-controlled sampling;
- thrust, flow, inventory, mass, impulse and ideal-delta-v plots;
- exact and equivalent-constant reduction reports.

**Acceptance**

- a zero-ramp unit-throttle profile degenerates bit-for-bit to B1;
- analytic constant/linear profiles match independent closed forms;
- sample refinement changes impulse, drain and delta-v below declared
  tolerances while the authoritative integral remains unchanged;
- every discontinuity appears once with documented pre/post semantics;
- cutoff occurs at the earliest binding constraint and is named;
- imported curve units/errors fail closed;
- total impulse equals integrated delivered thrust and mass loss equals
  integrated tank-drain flow;
- an approximate Sidera reduction cannot be produced without an explicit mode.

**Recommended version**

- `0.11.0` before 1.0, otherwise provisionally `1.2.0`.

### Phase B3 — state-dependent propulsion-system providers

**Current isolated implementation status (not integrated into Sidera)**

- implemented: rigid-tank polytropic blowdown pressure state and exact inverse
  minimum-pressure cutoff;
- implemented: pressure-correlated thrust/Isp/chamber-pressure law with an
  explicit validity range and no extrapolation;
- implemented: ordered measured pressure-performance table with bounded linear
  interpolation and source identity;
- implemented: mass-domain Simpson integration, stream/system-Isp closure,
  pressure/reserve cutoff and convergence-qualified tabulated output;
- implemented: stable primary-source registry, citation-bearing artifact hash,
  public docstring reference audit and analytic/monotonic/refusal tests;
- implemented: sourced continuous-parameter uncertainty with bounded
  distributions, aleatory/epistemic separation, optional correlated sampling,
  empirical statistics, sensitivity ranking and one retained propulsion
  artifact per realization;
- implemented: mutually exclusive discrete propulsion scenarios for
  engine-out/common-cause loss, ignition delay, early cutoff, response change
  and thrust/Isp degradation, with exhaustive probability checks and a complete
  artifact for every outcome;
- implemented: scenario-conditional continuous ensembles, persistence,
  hot-fire calibration, model/flight discrepancy separation and weighted
  Sidera hand-off artifacts;
- implemented: a framework-neutral propagation bridge returning axial force,
  total and named tank-state drains, exact start/boundary/stop/inventory event
  surfaces, bounded estimation drivers and analytic local Jacobians;
- implemented: hash-verified multidimensional operating-condition maps with
  independent force/flow response, analytic condition partials and explicit
  lower/upper validation-domain event roots;
- implemented: a six-state regulated feed provider and converged algebraic
  coupling to the pressure-dependent bridge, including positive-safe cutoff
  roots and implicit state/flow/acceleration derivatives;
- implemented: a two-state dynamic liquid feed line with inertance,
  compliance, nonlinear damping, stored-mass closure, positive-safe roots and
  engine/feed-composed analytic Jacobians;
- still open: higher-order distributed line/priming physics, regulator
  hysteresis, ullage stratification, bipropellant tank-pair depletion, the
  native Sidera wrapper around the bridge, Sidera-side pointing/epoch
  dispersion, external qualification-data closure and HTTP/UI.

**Deliver**

- pressure-fed regulated and isothermal blow-down provider interfaces;
- a validated preliminary monopropellant blow-down implementation;
- bipropellant tank-pair depletion with component-limited cutoff;
- rocket-model provider joining equilibrium, nozzle, loss and cycle outputs;
- ambient-pressure schedule support for propulsion-only altitude/test profiles;
- provider state schema/versioning and reproducible initial-state builder;
- test/qualification trace correlation metrics;
- burn sweeps, case comparison and seeded uncertainty using `studies`.

**Acceptance**

- constant-pressure limits converge to the B1 result;
- blow-down pressure, thrust and Isp evolve monotonically where the selected
  model requires it;
- analytic tank-law references and independent numerical integration agree;
- provider invalidity halts the run with the last valid state preserved;
- mass, energy/pressure model residuals and uncertainty seed are reported;
- standalone operation remains independent of Sidera.

**Recommended version**

- `0.12.0` before 1.0, otherwise provisionally `1.3.0`, with preliminary
  providers marked experimental until their external validation rows close.

### Phase B4 — optional exact transient Sidera coupling

**Dependency**

- the isolated `rocket_propulsion_propagator_bridge_v1` producer is now
  implemented and verified;
- Sidera still needs a released/public force-model wrapper able to consume its
  tabulated thrust, independent total/tank mass-flow, events and Jacobians.

**Deliver**

- capability-negotiated adapter for the new Sidera contract;
- native mapping of bridge event surfaces and named parameter partials;
- named mapping of propagated pressure/thermal/power states into performance
  conditions and composition of their analytic Jacobian columns;
- registration of feed-system additional states and their simultaneous ODE
  derivatives in the native Sidera numerical propagator;
- registration of dynamic line-flow/manifold-pressure states, replacement of
  the static injector-pressure root, and mapping of the composed two-state
  Jacobian;
- end-to-end trajectory example using a Rocket Propulsion Lab profile;
- constant/native versus constant/tabulated and transient comparisons;
- CCSDS OPM comparison for constant burns;
- trajectory-result linkage by source burn hash;
- a documented unsupported result if the installed Sidera lacks the contract.

**Acceptance**

- free-space transient cases reproduce analytic mass, velocity and position;
- a constant tabulated burn matches native Sidera `FiniteBurn` within solver
  tolerance;
- open-cycle cases close vehicle mass using independent total tank drain;
- adapter tests cover the supported Sidera version/capability matrix;
- unknown or changed Sidera schemas fail closed;
- every linked result records both versions, schemas, tolerances and hashes.

B4 is not a prerequisite for B1-B3 standalone releases.

### Phase B5 — operational, cluster and advanced providers

Deferred until the underlying providers are validated:

- multi-engine-out timelines, redundancy and restart reliability;
- body-frame force/moment history for mounted/gimballed clusters;
- thermal soak, cumulative on-time and cooldown availability;
- regulated pressurization and feed-line pressure-drop dynamics;
- hybrid regression coupling and solid grain/nozzle erosion providers;
- electric-propulsion power-limited performance maps;
- qualification-data residual dashboards and Bayesian calibration;
- burn targeting/optimization, which remains in Sidera or another mission
  layer and consumes this project's propulsion profiles.

## 9. API and UI plan

### 9.1 Python API is authoritative

The Python API is the primary reusable library surface; HTTP, CLI and UI call
the same functions.

```python
from rocket_propulsion.propulsion.burns import (
    BurnDefinition,
    BurnTarget,
    BurnTargetKind,
    ConstantPerformanceProvider,
    PropellantRole,
    TankDefinition,
    simulate_burn,
)

provider = ConstantPerformanceProvider.bipropellant(
    delivered_thrust_n=24_000.0,
    system_specific_impulse_s=328.0,
    oxidizer_fuel_ratio=3.4,
)
definition = BurnDefinition(
    name="420 m/s sizing burn",
    initial_mass_kg=2_400.0,
    protected_dry_mass_kg=1_400.0,
    tanks=(
        TankDefinition(
            "oxidizer", "oxidizer", PropellantRole.OXIDIZER, 760.0, 20.0
        ),
        TankDefinition("fuel", "fuel", PropellantRole.FUEL, 240.0, 8.0),
    ),
    target=BurnTarget(BurnTargetKind.IDEAL_DELTA_V, 420.0),
    maximum_duration_s=300.0,
)
result = simulate_burn(definition, provider)

print(result.summary.duration_s)
print(result.summary.total_impulse_n_s)
print(result.summary.ideal_velocity_increment_m_s)
print(result.summary.cutoff_reason)
```

Public functions:

- `simulate_burn(definition, provider, *, environment=...) -> BurnResult`;
- `size_burn(definition, provider, *, target=...) -> BurnSizingResult`;
- `reduce_to_constant(result, *, tolerance=...) -> ConstantBurnReduction`;
- `compare_burns(results, *, baseline=...) -> BurnComparison`;
- `burn_result_to_json` / `burn_result_from_json`;
- `export_burn(result, target, **options) -> ExportBundle`.

Provider construction shortcuts are conveniences over explicit immutable
models; saved artifacts always contain the resolved inputs.

### 9.2 HTTP API

`POST /api/v1/engine-operating-point`

- evaluates one engine/nozzle/cycle state without performing a burn;
- returns every propellant stream and chamber/system Isp distinction.

`POST /api/v1/burn/simulate`

- runs a complete standalone burn from a definition and provider;
- returns profile, summary, warnings, evidence and closure residuals.

`POST /api/v1/burn/size`

- solves duration/propellant for a selected target;
- returns feasibility bounds and the completed burn when feasible.

`POST /api/v1/burn/reduce`

- constructs and qualifies an equivalent constant burn;
- reports exactness, variation metrics and reduction error indicators.

`POST /api/v1/burn/export/sidera`

- accepts a burn result plus export-only start/frame/direction;
- returns exact native-compatible schema or a typed refusal;
- requires an explicit `mode="equivalent_constant"` to permit approximation.

`GET /api/v1/integrations/sidera/capabilities`

- reports unavailable/supported/incompatible without changing state;
- includes the detected version, schema and supported burn types.

No endpoint will propagate an orbit. The exported artifact is the hand-off to
Sidera.

All responses use the existing API envelope and carry model metadata. Large
profiles may be downsampled only for the response when the authoritative
integrals and full artifact remain unchanged; downsampling is reported.

### 9.3 CLI

Planned dependency-free commands:

```text
rocket-propulsion burn simulate burn-input.json --output burn-result.json
rocket-propulsion burn size burn-input.json --target delta-v:420
rocket-propulsion burn compare case-a.json case-b.json --output report.html
rocket-propulsion burn export case.json --target sidera --output plan.json
rocket-propulsion burn validate burn-result.json
```

CLI output defaults to a concise human summary; `--json` is machine-readable.
Failures have stable non-zero exit codes for input, infeasible target,
convergence, provider invalidity and integration incompatibility.

### 9.4 UI

Add one **Burn** workbench with:

- engine/propellant/operating-point selection;
- initial mass and separate propellant inventories;
- termination target and throttle/transient controls;
- summary cards for duration, impulse, delta-v, final mass and margins;
- thrust, mass-flow, vehicle-mass and inventory histories;
- an export drawer for Sidera start time, frame and direction;
- exact/approximate compatibility badge and all warnings.

The workbench has four tabs:

1. **System** — engine provider, cycle, tanks, inventories, reserve and cluster;
2. **Command** — target, throttle, start/stop profile and limits;
3. **Results** — summary, closures, binding constraint and histories;
4. **Export** — standalone artifact/report plus optional Sidera hand-off.

The result view must keep chamber Isp and system Isp visually distinct. It must
show gross thrust, delivered thrust and axial effective thrust separately when
losses or cant are enabled.

The UI must label ideal delta-v as bookkeeping. Achieved orbital changes are
displayed only when imported back from Sidera.

### 9.5 Standalone reporting

Every burn can be exported without Sidera as:

- canonical JSON with complete inputs/results/provenance;
- tidy CSV profile plus a separate summary CSV;
- SVG charts;
- printable HTML engineering report;
- optional compact `thrust_n, mass_flow_kg_s, time_s` exchange CSV.

The HTML report includes equations used, fidelity level, validity warnings,
provider/reference versions, inventory ledger, closure residuals, uncertainty
settings and the exact canonical result hash.

## 10. Validation matrix

### 10.1 Evidence layers

Every production model needs evidence at all applicable layers:

1. **Equation tests** — one formula against hand/independent values;
2. **Invariant tests** — mass, flow, impulse and monotonicity properties;
3. **Limit/metamorphic tests** — zero ramp, constant-pressure, zero cant,
   engine-count scaling and profile refinement;
4. **Reference cases** — published data or an independently identified tool;
5. **Cross-interface parity** — Python, API, CLI, workspace and report;
6. **Adapter tests** — only when an external integration is installed;
7. **Schema migration tests** — every persisted historical version;
8. **UI smoke tests** — input, calculation, warning, save/open and export.

An internal golden is a regression aid, not independent validation. The
validation report identifies which rows are analytic, external, cross-tool or
regression-only.

### 10.2 Minimum case matrix

| Case | Standalone evidence | Sidera evidence |
| --- | --- | --- |
| Constant thrust, constant Isp | closed-form mass, impulse and `integral(F/m dt)` | free-space/point-mass parity |
| Constant bipropellant | O/F split, reserve and inverse-sizing closure | one exact `FiniteBurn` |
| Gas-generator engine | chamber versus dump versus total tank drain | system-Isp mass parity |
| Pressure-fed constant limit | regulated provider collapses to L0 | exact constant export |
| Blow-down monopropellant | tank-law and thrust/Isp trend reference | tabulated path when available |
| O/F-limited bipropellant | first component cutoff and remaining inventory | final mass only |
| Multi-engine cluster | force sum, cant and engine-out metamorphic tests | scalar effective-thrust export |
| Linear start/stop ramps | analytic impulse, drain and centroid | transient path when available |
| Throttle step profile | exact discontinuities and refinement convergence | transient path when available |
| Solid/imported curve | source reconstruction and integrated impulse | transient path when available |
| Short-burn reduction | standalone reduction metrics | native versus reduced trajectory |
| Long-burn reduction | warning/refusal threshold | trajectory-error demonstration |
| Infeasible delta-v | attainable bound and binding inventory | no export |
| Schema/version mismatch | canonical typed error | fail-closed adapter behavior |

Each case tests equations, boundary validation, canonical serialization and at
least one user-facing path. Sidera columns marked unavailable remain open gates;
they do not block standalone validation.

### 10.3 Provisional tolerance policy

Numerical tolerances live in `docs/TOLERANCES.md` and reference-data files,
never only in test code. Initial targets, subject to conditioning evidence:

| Quantity | Initial acceptance |
| --- | --- |
| Closed-form L0 mass/impulse/delta-v | `1e-12` relative or scale-aware absolute floor |
| Inventory and stream mass closure | `1e-10 kg` absolute plus `1e-12` of initial mass |
| Analytic L1 linear-segment integrals | `1e-11` relative |
| State-dependent provider integration | tolerance study; initially `1e-8` relative |
| Inverse-sizing forward/back closure | `1e-10` relative in declared conditioned range |
| Canonical JSON round-trip/hash | exact bytes/hash |
| Constant Sidera mass/delta-v parity | `1e-10` relative |
| Tabulated/native Sidera parity | tied to Sidera solver tolerance and convergence study |

Near zero values use absolute tolerances derived from physical scale. A test
must not divide by a near-zero reference to manufacture a relative error.

### 10.4 Independent reference sources

- NASA SP-125 provides the liquid-engine design/performance foundation:
  <https://ntrs.nasa.gov/citations/19710019929>.
- NASA's rocket thrust equation defines momentum and pressure thrust:
  <https://www1.grc.nasa.gov/beginners-guide-to-aeronautics/rocket-thrust-equation/>.
- NASA CR-184099 documents a transient engine simulation that covers steady
  balance, throttle transients, starts and shutdowns:
  <https://ntrs.nasa.gov/citations/19910011919>.
- NASA's transient-modeling overview motivates explicit startup/shutdown
  models and redline/validity handling:
  <https://ntrs.nasa.gov/citations/20040000363>.
- NASA CR-131400 supplies a spacecraft monopropellant blow-down reference
  context: <https://ntrs.nasa.gov/citations/19730012094>.
- NASA SP-8112 supplies the pressure-fed/blowdown pressurization architecture
  and ullage-decay basis: <https://ntrs.nasa.gov/citations/19760015212>.
- NASA CR-122347 supplies measured supply-pressure versus thrust/Isp and pulse
  transient context: <https://ntrs.nasa.gov/citations/19720011120>.
- NASA-STD-7009B supplies model-credibility, validation and uncertainty-
  qualification requirements:
  <https://standards.nasa.gov/standard/NASA/NASA-STD-7009>.
- NASA/SP-2011-3421 supplies Monte Carlo/LHS uncertainty propagation and
  sensitivity-analysis guidance:
  <https://ntrs.nasa.gov/citations/20120001369>.
- NASA's liquid-engine reliability-driver study identifies engine-out,
  start/cutoff, restart and clustered-engine reliability concerns:
  <https://ntrs.nasa.gov/citations/20050207429>.
- NASA's launch-vehicle common-cause study motivates explicit dependent
  scenarios instead of multiplying independent engine probabilities:
  <https://ntrs.nasa.gov/citations/20160007073>.
- CCSDS 502.0-B-3 is the interchange authority for OPM maneuver blocks:
  <https://ccsds.org/Pubs/502x0b3e1.pdf>.
- The inspected Sidera tests remain the cross-project implementation oracle,
  not an independent propulsion-physics source.

Reference selection for a provider is recorded in its metadata. A provider
without an external reference is labelled preliminary or experimental in
Python, API, UI and reports using the same wording.

### 10.5 Test organization

```text
tests/
├── test_burn_models.py
├── test_burn_constant.py
├── test_burn_profiles.py
├── test_burn_inventory.py
├── test_burn_sizing.py
├── test_burn_providers.py
├── test_burn_serialization.py
├── test_burn_api.py
├── test_burn_workspace.py
├── test_burn_reports.py
├── integration/
│   └── test_sidera_burn_adapter.py
└── reference_data/
    └── burns/
        ├── constant_cases.json
        ├── transient_cases.json
        ├── blowdown_cases.json
        └── source_manifest.json
```

The normal suite skips no standalone burn tests. Only the external Sidera
integration module is optional and reports a clear skip reason/version.

## 11. Explicit non-goals

- implementing orbital mechanics in Rocket Propulsion Lab;
- computing RIC/VNB/body transforms here;
- hiding gravity or steering losses inside propulsion Isp;
- claiming an equivalent constant burn is exact for a long transient burn;
- silently replacing total system flow with chamber flow;
- importing Sidera into the propulsion domain package;
- silently inferring ambient pressure, available power, tank pressure or O/F;
- performing guidance, burn targeting or orbit optimization in the burn core;
- treating body-frame force/moment output as an attitude simulation;
- claiming preliminary blow-down or transient correlations are qualification
  models without external data;
- accepting negative propellant mass to help an optimizer converge;
- adding CFD, finite-rate combustion or detailed turbopump maps as a
  prerequisite for the first burn release;
- flight software, command generation, certification or safety-critical use.

## 12. File-by-file implementation backlog

### Slice 1 — model foundation

- add `propulsion/burns/models.py` with frozen public records and enums;
- add `propulsion/burns/inventory.py` with stream/tank ledgers;
- add canonical `rocket_propulsion_burn_v1` encoder/decoder;
- add schema and fixtures under `tests/reference_data/burns/`;
- extend `core.errors` only with stable burn-specific codes, not new exception
  families unless behavior differs from existing `InputError`, `DomainError`
  and `ConvergenceError`;
- document quantities and tolerances before implementing equations.

### Slice 2 — standalone L0 solver

- implement `ConstantPerformanceProvider`;
- implement closed-form target solvers and cutoff feasibility;
- implement independent conservation residual calculation;
- publish the Python API from `propulsion.burns.__init__`;
- add equation, inverse, limit and property-grid tests.

### Slice 3 — product surfaces

- add routes and route metadata;
- add CLI subcommands without changing the server command;
- add workspace schema migration and deterministic burn cases;
- add JSON/CSV/HTML report sections;
- build the Burn UI from returned Python profile data;
- test Python/API/CLI/workspace parity by canonical result hash.

### Slice 4 — current Sidera adapter

- add `integrations/contracts.py` and `integrations/sidera.py`;
- add capability inspection and exact L0 mapping;
- add artifact-only export when native Sidera import is unavailable;
- add optional local/candidate Sidera integration test;
- document the supported capability tuple, not only a version range.

### Slice 5 — L1 profiles

- implement schedule normalization and discontinuity splitting;
- implement analytic constant/linear segment integrals;
- implement trace import and reconstruction report;
- implement reduction metrics and explicit approximation policy;
- add profile charts, comparison and uncertainty mapping.

No slice starts by modifying Sidera. A Sidera tabulated-burn proposal is opened
only after the standalone L1 artifact and constant parity case are stable.

## 13. Schema, compatibility and migration policy

- Public burn Python models and `rocket_propulsion_burn_v1` use semantic
  versioning independent of Sidera schemas.
- Additive optional JSON fields do not change the schema identifier; changes
  to meaning, required fields or units do.
- Readers reject unknown major schemas and preserve unknown additive metadata
  only when lossless pass-through is explicitly supported.
- `.rplab.json` advances from schema v2 through a migration that adds burn
  cases without reinterpreting existing cases.
- Result hashes cover normalized physical inputs, provider identity/version,
  solver policy and authoritative results; presentation sampling and chart
  preferences use a separate hash.
- External adapter artifacts carry their own target schema and source hash.
- Deprecated Python fields remain readable for one minor release with a clear
  migration warning; physical semantics are never silently changed.

## 14. Risk register and design decisions

| Risk | Consequence | Control / decision gate |
| --- | --- | --- |
| Chamber Isp used as system Isp | under-predicted propellant drain | explicit streams and mandatory system-Isp closure |
| Double-applied throttle | wrong thrust and mass flow in Sidera | export actual thrust with Sidera throttle fixed to unity |
| Averaged Isp inserted into Tsiolkovsky | wrong variable-profile delta-v | integrate `F/m`; report distinct Isp definitions |
| Profile sampled too coarsely | impulse/cutoff bias | analytic authoritative integrals and convergence evidence |
| Discontinuity crossed by a solver | transient smearing | mandatory segment splitting |
| One propellant depletes first | negative component inventory | per-tank terminal roots and component-limited cutoff |
| External provider leaves validity range | plausible but invalid answer | terminate at last valid state, fail closed |
| Sidera API changes | corrupt export | capability probe, schema fixture and supported matrix |
| Local Sidera checkout leaks into package | non-portable release | adapter optional extra; path only in explicit test config |
| Burn core grows into GNC/trajectory | duplicated physics and ownership | non-goals and import-layer tests |
| Preliminary model mistaken for flight model | unsafe interpretation | fidelity badge, source/validity report and consistent warnings |
| External CEA gate still open | mixed validation status | burn releases remain pre-1.0 and identify independent open gates |

Decision gates:

1. B0 model review before persistence/API code;
2. L0 conservation and inverse sizing before any Sidera adapter;
3. Python/API/CLI parity before UI publication;
4. standalone L1 reference cases before proposing tabulated Sidera input;
5. external blow-down reference closure before removing `experimental`;
6. no stable 1.0 burn claim while project-wide required validation gates are
   open.

## 15. Definition of done for a burn capability

A burn feature is complete only when:

- equation, inputs, outputs, units and validity are documented;
- public data models are immutable and canonical-serializable;
- domain, limit, inverse and conservation tests pass;
- at least one independent or analytic reference row is recorded;
- Python, API and saved-workspace results agree;
- reports expose warnings and residuals without requiring debug output;
- unsupported adapters refuse rather than approximate silently;
- all new files have focused module/package documentation;
- Ruff, full tests, wheel installation and clean import pass;
- the D-drive canonical project and desktop junction copy are synchronized.

For a Sidera-connected feature, add:

- public API capability probe;
- version/schema/hashes in the artifact;
- cross-project parity test;
- exact/approximate/unsupported classification in the UI and report.

## 16. Recommended first implementation slice

Start with B0 and the standalone part of B1. The first milestone must be useful
even if Sidera is unavailable.

The first standalone vertical demonstration should be:

1. select a vacuum engine operating point already computed by Rocket
   Propulsion Lab;
2. size a constant burn from initial mass and requested delta-v;
3. close fuel, oxidizer, dump-flow, total mass and loss metadata;
4. display duration, impulse, system Isp, ideal 1-D delta-v and margins;
5. save/reopen the canonical burn artifact and generate an HTML report;
6. reproduce the same result through Python, API and CLI.

The second demonstration adds the adapter:

1. export the verified result as a native-compatible
   `sidera_maneuver_plan_v1` artifact;
2. load it in Sidera and verify identical mass depletion and ideal delta-v;
3. let Sidera alone report trajectory/orbit changes;
4. link the trajectory result to the authoritative burn hash.

Only after this parity test is green should transient profiles or a Sidera
tabulated-burn extension begin.
