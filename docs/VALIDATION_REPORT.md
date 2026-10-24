# Software verification evidence — 0.9.0 preview

This table describes software checks and their intended model scope. It is
not a physical validation or qualification report. No selected motor/feed
hardware has accepted relevant-environment measurement, calibration or
holdout evidence in this review. Hardware readiness is not assessed.
Estimator tests and synthetic hot-fire fixtures do not establish that a
physical engine has been calibrated. Source references/hashes do not establish
measurement authenticity or model validity.

Current test identities, source/dependency digests and results are generated
in `output/assurance/verification_summary.json`. The per-symbol review remains
partial (`objective_audit.json`); a passing suite does not establish the
adequacy of every check or a universal correctness guarantee. The table's
scope labels must be read within these limits. See
[motor/feed assurance gaps](PROPULSION_ASSURANCE.md).

| Capability | Listed software evidence | Verification scope / open work |
|---|---|---|
| Bounded measurement/model comparisons | exact interval corners, outward endpoint containment, guarded decision checks and archived-file integrity | numerical protocol checks only; physical holdout/calibration and model-source authentication unestablished |
| Ideal gas, isentropic, shocks | source-labelled reference fixtures and inverse/limit tests | listed software checks |
| NASA Glenn species and mixtures | published coefficients, 298 K checks, range and inverse tests | software checks cover bundled species |
| Equilibrium provider contract | TP/HP/UV element and energy residual tests | software checks cover experimental |
| NASA CEA comparison | adapter/parser/deck tests; official example recorded | **external nine-point gate open** |
| Off-design nozzle | analytical regime order, internal-shock mass conservation, atmosphere sweep | software checks cover perfect gas |
| Separation | NASA SP-8041/Summerfield screen and validity warning | screening only |
| Rao/MOC geometry | legacy boundary checks plus independent raw/mapped characteristic-chord diagnostic; see [geometry review](GEOMETRY_REVIEW.md) | mapped interior flow, ineffective axisymmetric toggle, maximum-count self-refinement, generic boundaries and construction claims remain open |
| Interior MOC unit process | independent constant/radial/spherical flows, local refinement, Cartesian balances, range/coordinate-cell/hostile replay checks and controlled faults; see [unit-process review](MOC_UNIT_PROCESS_REVIEW.md) | constant-gamma interior intersection only; axis/sonic/wall processes, complete field/exit validation and physical qualification remain open; legacy mapped solver unchanged |
| MOC axis and prescribed-wall processes | independent signed radial/spherical boundary flows, refinement, chained axis history, containment/topology and controlled faults; see [boundary review](MOC_BOUNDARY_REVIEW.md) | local constant-gamma boundary processes only; sonic initialization, full characteristic-field assembly, optimized wall and physical/hardware validation remain open |
| CAD exchange | independent directed-edge/normal/analytic-volume oracles, exact serialized coordinate round trips, extreme-scale/hostile records and controlled faults; see [geometry review](GEOMETRY_REVIEW.md) | sampled closed fluid-domain software exchange only; continuous geometry/flow, third-party CAD conformance, nozzle-wall manufacture and hardware remain unqualified |
| Loss-factor model | declared operating bases, exact ledgers, replay and first-order covariance tests | conditional coefficient model; coupled nozzle physics and statistical coverage remain unqualified |
| Prescribed Bartz profile | independent signed profile equations, Gaussian surface-integration oracle, guarded geometry, retained replay and controlled faults; see [Bartz review](BARTZ_REVIEW.md) | integral ledger identity is not measured energy closure; prescribed wall/gas states, empirical domain, unsampled-profile error and hardware applicability remain unqualified |
| Bulk regenerative cooling | exact energy/laminar ledgers, enclosed turbulent formula, exact limit/regime decisions, explicit cross-section factor, owned replay and controlled faults; see [regenerative review](REGENERATIVE_REVIEW.md) | prescribed single-phase model; transition/cross-section/correlation applicability, wall coupling, changing properties, boiling and hardware validation remain open |
| Workspaces and uncertainty | exact round trip, v1 migration, seeded reproducibility | listed software checks |
| Staging and cycles | target/constraint tests and component power balance | software checks cover constant-property scope |
| Constant propulsion burn | analytic inverse cases, stream/tank conservation, cutoff, provider construction, persistence and Python/API/CLI/workspace parity tests | software checks cover L0 preliminary scope |
| Profiled propulsion burn | analytic ramp/step integrals, exact cutoff inversion, open-cycle scaling, v2 reproduction and Python/API/HTTP/CLI parity | software checks cover proportional-throttle L1 scope |
| Transient engine response | analytic first-order step exposure, first moment, zero-lag limit and convergence-qualified sampling | software checks cover declared first-order model |
| Engine sequence and cluster | minimum-pulse/cooldown/restart refusal, staggered ignition, cant, disabled-engine and force/mass composition tests | software checks cover isolated propulsion scope |
| Tabulated burn contract | discontinuity-side semantics, stream closure, source/artifact hashes, tamper refusal, CSV linear/hold replay and CLI round trip | software checks cover `rocket_propulsion_tabulated_burn_v1` prototype |
| Pressure-fed blowdown | analytic polytropic endpoints/inverse cutoff, monotonic pressure/thrust/mass history, system-Isp flow identity, pressure-map interpolation/refusal, Simpson-to-artifact residuals | software checks cover declared preliminary L2 assumptions |
| Burn source/docstring traceability | unique primary-source registry, official-domain URL checks, citation-bearing artifact hashes and public transient API `References` audit | software checks cover isolated transient modules |
| Propulsion uncertainty | seeded reproducibility, bounded marginal transforms, singular-correlation refusal, empirical copula correlation, nominal/model consistency, pressure-envelope refusal and per-realization artifact closure | software checks cover preliminary continuous-parameter UQ scope |
| Discrete propulsion scenarios | probability completeness, single-engine-out, common-cause all-engine loss, exact zero-thrust history, early cutoff, delay, performance degradation and conflicting-event refusal | software checks cover propulsion-consequence scope; trajectory risk remains in Sidera |
| Hierarchical propulsion ensemble | exact scenario/conditional/joint weight closure, deterministic artifact reproduction, evidence-sensitive hashes, thrust/Isp/time physical identities, axial-force refusal, zero-thrust retention, scenario-contribution closure and independent audit seeds | software checks cover conditional scale-based propulsion scope; trajectory ensemble remains in Sidera |
| Hierarchical ensemble persistence | complete JSON round trip, outer SHA-256, nested burn hash refusal, recomputed probability/statistic/contribution/evidence checks, and within/between/total variance identity | software checks cover `rocket_propulsion_hierarchical_ensemble_v1`; Sidera consumer remains open |
| Hot-fire scale calibration | NIST §7.4.4 variance-component reproduction, unbalanced engine groups, heteroscedastic measurement uncertainty removal, nonnegative boundary behavior, condition/identity refusal, deterministic hashes and population/same-engine residual mapping | software checks cover the declared one-factor random-effects model; variance-component sampling distributions remain open |
| Paired hot-fire covariance calibration | full positive-semidefinite measurement covariance validation, analytic cross-covariance removal, scalar/component diagonal identity, bounded correlated ensemble draws, recorded positive-definite regularization, complete JSON round trip and derived-result tamper refusal | software checks cover paired method-of-moments covariance with Mandel-Paule marginals; multivariate REML/Bayesian inference remains open |
| Model/flight discrepancy calibration | independent ratio evidence, application-domain pooling refusal, consensus-plus-predictive epistemic variance, hash and semantic round trip, evidence-reuse refusal, aleatory/epistemic component retention, and thrust/Isp/flow identities through a complete hierarchical ensemble | software checks cover constant in-domain random-effects discrepancy; response surfaces and extrapolation remain open |
| Propagator propulsion bridge | verified source-artifact reopening, force/total-mass/tank-state closure, explicit left/right jump limits, time/dry-mass/reserve roots, zero outside/inhibited behavior, bounded named drivers, analytic mass/parameter Jacobians checked by central finite differences, and configuration-sensitive manifest hash | software checks cover the framework-neutral bridge; native Sidera force-model wrapper remains intentionally isolated |
| Operating-condition performance surface | complete Cartesian-grid validation, arbitrary-dimensional multilinear interpolation, analytic force/flow condition gradients checked against central differences, strict no-extrapolation refusal, lower/upper domain event roots, hash-verified JSON round trip, and exact bridge composition | software checks cover deterministic rectilinear maps; node covariance and condition-dependent discrepancy remain open |
| Regulated feed dynamics | ideal-gas pressure reproduction, bottle-to-ullage pressurant mass closure, independent open-system energy residuals, choked-flow branch, valve command/lag, feed-line loss, positive-safe inventory/pressure events, analytic pressure partials, converged feed/engine algebraic closure, and implicit coupled derivatives checked against fully re-solved finite differences | software checks cover preliminary single-node ideal-gas scope; multi-node priming dynamics, real fluids and stratification remain open |
| Dynamic liquid feed line | steady operating-point reproduction, compliant storage/tank/engine mass closure, natural-frequency and damping diagnostics, positive-safe reverse-flow/pressure/flow roots, configuration hash, inactive-engine relaxation, double-loss refusal, and local/composed/feed-state/acceleration Jacobians checked against central finite differences | software checks cover a first-mode lumped surrogate; calibration for selected hardware is unaccepted; distributed water hammer, priming, cavitation, two-phase flow and structural coupling remain open |
| Sidera weighted ensemble hand-off | parent/nested/outer hash chain, normalized frame direction, exact joint closure, zero-thrust retention, axial-force/total-flow semantics, semantic-tamper refusal and explicit missing-capability report | producer contract checks listed; native Sidera consumer remains intentionally absent |
| Sidera hand-off | artifact force/mass identity and native public-contract test against local Sidera environment | software checks cover exact constant burns |
| Canted constant Sidera hand-off | axial impulse maps to native force while rocket-equivalent Isp reproduces total tank drain | software checks cover constant cant efficiency |
| Equivalent Sidera reduction | duration/impulse/propellant identities, variation and centroid evidence, explicit opt-in requirement | software checks cover declared approximation |
| Browser UI | HTTP/HTML contract checks and current JavaScript syntax check; earlier visual QA has not been revalidated for the current changes | current browser visual check open |
| Python packaging | 3.11/3.12 CI matrix and wheel/sdist job configuration | configured; current package build not part of this review |

## Release decision

Version `0.9.0` labels a local engineering preview. The plan still requires
nine independent CEA reference points across
LOX/LH2, LOX/RP-1, and LOX/CH4. The adapter deliberately requires a
user-supplied CEA installation, so that gate cannot be replaced by self-generated
local-solver values. Closing this solver-comparison gate alone would not
establish motor/feed hardware maturity, quantified physical uncertainty or
qualification. The separate hardware evidence gates remain open.

## References

The evidence rows above cite their sources by name; full entries and URLs are
in [`REFERENCES.md`](REFERENCES.md). The principal ones: NASA Glenn
coefficients, NASA/TP-2002-211556, <https://ntrs.nasa.gov/citations/20020085330>;
NASA RP-1311 (CEA), <https://ntrs.nasa.gov/citations/19950013764>; NASA SP-8041
(nozzle separation), <https://ntrs.nasa.gov/citations/19710021390>; NIST/SEMATECH
e-Handbook §7.4.4, <https://www.itl.nist.gov/div898/handbook/prc/section4/prc44.htm>. Physical constants and their defining
sources are registered in `rocket_propulsion.core.constants` and listed under
*Physical constants* there.
