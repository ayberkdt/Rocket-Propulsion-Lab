# Review resumed on the current upstream tree

The upstream baseline is commit
`fea6a324e9b9e8dc7f88872e6abfddbad8b4842a` (merged PR #9).
Relative to the previous published baseline `59c51c8`, it contains 16 commits
and changes 320 files (115232 inserted and 4807 deleted lines). Line counts
include documentation and generated/catalogue data; they are not reviewed
function counts.

## Added scope

- A new browser workbench, vendored display presets and a launcher.
- A larger NASA thermochemistry catalogue, optional NASA CEA 3.3.4 provider,
  solver comparison fixtures and rocket performance calculations.
- Reference and cubic real-fluid equations of state, transport properties
  and equilibrium derivatives.
- Turbopump maps, water hammer, regulator hysteresis, ullage stratification
  and calibration fields.
- Advanced compressible flow, external and launch-vehicle aerodynamics.
- Conduction, convection, radiation and a coupled regenerative-cooling march.
- Combustion-efficiency bookkeeping, finite-area combustion, condensed
  equilibrium and finite-rate nozzle kinetics.

The new capabilities are review inputs. Their presence and existing tests do
not establish that every function, range guard or source claim has been audited.

## Reconciliation policy

The previous dirty checkout is retained at the sibling directory
`advanced-propulsion-sidera`. Its NASA9 and mixture changes were not copied over
the new 36-species catalogue and equilibrium implementation. Those changes
need a semantic port and renewed comparison, not a wholesale file replacement.
The shared-root repair is transferred separately; its scope is documented in
[SHARED_ROOT_REVIEW.md](SHARED_ROOT_REVIEW.md).

The historical `post_mixture_forward_checkpoint` proof applies to its archived
source/dependency snapshot. It does not qualify this tree. In particular,
historical test totals and the 716 curated faults are not current results.
The current review keeps a source-bound regression report and a separate root
mutation report. Missing dependencies, skips, timeouts or errors remain open
results; the per-module 120 s / collection 30 s / four-worker policy is retained.
The source fingerprint now includes browser `.mjs`/`.cjs` modules, CI workflows,
the launcher and the license. Mutating those executable/configuration inputs
must invalidate prior results as well. Missing-provider tests simulate import
absence explicitly so they run with the real optional provider installed.
CI installs the CEA extra, enables the opt-in slow comparisons and checks the
NASA catalogue source against its pinned SHA-256 before running tests.
The first run with all slow cases enabled exposed a retained 120 s timeout.
Repeated finite-difference column work was then removed from the equilibrium
Newton loop; see [EQUILIBRIUM_JACOBIAN_REVIEW.md](EQUILIBRIUM_JACOBIAN_REVIEW.md).
This is a scoped numerical change, not a renewed review of the entire solver.

## Current priorities and limits

1. Repair the reproduced Ruff 0.17.2 diagnostics and the shared-root defects;
   test them against the current tree and actual optional NASA CEA provider.
2. Reconcile prior thermochemistry corrections with the expanded catalogue,
   transport and new equilibrium paths without overwriting upstream work.
3. Review new numerical models and their claimed domains with independent
   reference cases and adversarial boundary tests, including coupled cooling
   and finite-rate chemistry. Preserve counterexamples and compare source
   versions explicitly.
4. Renew the prior mutation and integration evidence on the reconciled tree.

Hardware selection, mission requirements, relevant-environment measurements,
uncertainty/holdout evidence and scaling remain open. Software changes and
reference-code agreement do not confer motor/feed TRL 5 or certification.

## References

- [Upstream baseline fea6a32](https://github.com/ayberkdt/Rocket-Propulsion-Lab/commit/fea6a324e9b9e8dc7f88872e6abfddbad8b4842a),
  retained locally as a full Git archive before the review changes.
- NASA NPR 7123.1D, [Appendix E](https://nodis3.gsfc.nasa.gov/displayDir.cfm?Internal_ID=N_PR_7123_001D_&page_name=AppendixE),
  hardware and software readiness criteria; no TRL is assigned by this review.
