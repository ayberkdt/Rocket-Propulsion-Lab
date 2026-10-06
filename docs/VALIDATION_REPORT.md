# Validation report — 0.9.0

| Capability | Evidence | Status |
|---|---|---|
| Ideal gas, isentropic, shocks | source-labelled reference fixtures and inverse/limit tests | closed |
| NASA Glenn species and mixtures | published coefficients, 298 K checks, range and inverse tests | closed for bundled species |
| Equilibrium provider contract | TP/HP/UV element and energy residual tests | closed as experimental |
| NASA CEA comparison | adapter/parser/deck tests; official example recorded | **external nine-point gate open** |
| Off-design nozzle | analytical regime order, internal-shock mass conservation, atmosphere sweep | closed for perfect gas |
| Separation | NASA SP-8041/Summerfield screen and validity warning | screening only |
| Rao/MOC geometry | common contract, wall/centerline residuals, 2× mesh refinement | closed for documented approximation |
| CAD exchange | CSV/SVG and parsed DXF/STL dimension tests | closed for preliminary exchange |
| Loss and thermal models | switch-off tests, budget/energy residuals, domain warnings | closed for preliminary scope |
| Workspaces and uncertainty | exact round trip, v1 migration, seeded reproducibility | closed |
| Staging and cycles | target/constraint tests and component power balance | closed for constant-property scope |
| Constant propulsion burn | analytic inverse cases, stream/tank conservation, cutoff, provider construction, persistence and Python/API/CLI/workspace parity tests | closed for L0 preliminary scope |
| Profiled propulsion burn | analytic ramp/step integrals, exact cutoff inversion, open-cycle scaling, v2 reproduction and Python/API/HTTP/CLI parity | closed for proportional-throttle L1 scope |
| Transient engine response | analytic first-order step exposure, first moment, zero-lag limit and convergence-qualified sampling | closed for declared first-order model |
| Engine sequence and cluster | minimum-pulse/cooldown/restart refusal, staggered ignition, cant, disabled-engine and force/mass composition tests | closed for isolated propulsion scope |
| Tabulated burn contract | discontinuity-side semantics, stream closure, source/artifact hashes, tamper refusal, CSV linear/hold replay and CLI round trip | closed for `rocket_propulsion_tabulated_burn_v1` prototype |
| Pressure-fed blowdown | analytic polytropic endpoints/inverse cutoff, monotonic pressure/thrust/mass history, system-Isp flow identity, pressure-map interpolation/refusal, Simpson-to-artifact residuals | closed for declared preliminary L2 assumptions |
| Burn source/docstring traceability | unique primary-source registry, official-domain URL checks, citation-bearing artifact hashes and public transient API `References` audit | closed for isolated transient modules |
| Propulsion uncertainty | seeded reproducibility, bounded marginal transforms, singular-correlation refusal, empirical copula correlation, nominal/model consistency, pressure-envelope refusal and per-realization artifact closure | closed for preliminary continuous-parameter UQ scope |
| Discrete propulsion scenarios | probability completeness, single-engine-out, common-cause all-engine loss, exact zero-thrust history, early cutoff, delay, performance degradation and conflicting-event refusal | closed for propulsion-consequence scope; trajectory risk remains in Sidera |
| Hierarchical propulsion ensemble | exact scenario/conditional/joint weight closure, deterministic artifact reproduction, evidence-sensitive hashes, thrust/Isp/time physical identities, axial-force refusal, zero-thrust retention, scenario-contribution closure and independent audit seeds | closed for conditional scale-based propulsion scope; trajectory ensemble remains in Sidera |
| Hierarchical ensemble persistence | complete JSON round trip, outer SHA-256, nested burn hash refusal, recomputed probability/statistic/contribution/evidence checks, and within/between/total variance identity | closed for `rocket_propulsion_hierarchical_ensemble_v1`; Sidera consumer remains open |
| Hot-fire scale calibration | NIST §7.4.4 variance-component reproduction, unbalanced engine groups, heteroscedastic measurement uncertainty removal, nonnegative boundary behavior, condition/identity refusal, deterministic hashes and population/same-engine residual mapping | closed for the declared one-factor random-effects model; model discrepancy remains open |
| Paired hot-fire covariance calibration | full positive-semidefinite measurement covariance validation, analytic cross-covariance removal, scalar/component diagonal identity, bounded correlated ensemble draws, recorded positive-definite regularization, complete JSON round trip and derived-result tamper refusal | closed for paired method-of-moments covariance with Mandel-Paule marginals; multivariate REML/Bayesian inference and model discrepancy remain open |
| Sidera weighted ensemble hand-off | parent/nested/outer hash chain, normalized frame direction, exact joint closure, zero-thrust retention, axial-force/total-flow semantics, semantic-tamper refusal and explicit missing-capability report | producer contract closed; native Sidera consumer remains intentionally absent |
| Sidera hand-off | artifact force/mass identity and native public-contract test against local Sidera environment | closed for exact constant burns |
| Canted constant Sidera hand-off | axial impulse maps to native force while rocket-equivalent Isp reproduces total tank drain | closed for constant cant efficiency |
| Equivalent Sidera reduction | duration/impulse/propellant identities, variation and centroid evidence, explicit opt-in requirement | closed as declared approximation |
| Browser UI | desktop visual QA, Burn target/inventory rendering, clean console, keyboard focus/skip link, responsive and print CSS | closed locally |
| Python packaging | 3.11/3.12 CI matrix, wheel/sdist job | configured; local wheel verified at release |

## Release decision

Version `0.9.0` is a feature-complete local preview. It is not labelled
`1.0.0` because the plan requires nine independent CEA reference points across
LOX/LH2, LOX/RP-1, and LOX/CH4. The adapter deliberately requires a
user-supplied CEA installation, so that gate cannot be replaced by self-generated
local-solver values.

