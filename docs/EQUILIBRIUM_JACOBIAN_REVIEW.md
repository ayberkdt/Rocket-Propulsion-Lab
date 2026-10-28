# Element-potential residual and finite-difference reuse

The opt-in slow combustion-efficiency tests exposed a 120.03 s module timeout
in the four-worker bounded regression. The source-frozen failed run is retained
at `output/assurance/history/current_main_slow_regression_timeout/manifest.json`.
The timeout is incomplete verification, not a pass.

`LocalEquilibriumProvider._composition_at_tp` previously recomputed the same
shifted residual vector separately for every Jacobian row. Each vector is a
pure calculation of the same captured coefficients and shifted potentials.
For E elements, the matrix construction therefore requested E squared shifted
vectors per Newton iteration, although E vectors contain every required entry.
The repair evaluates each column once and uses its corresponding row entries.
The finite-difference step, equations, floating arithmetic per matrix entry,
linear solve, damping and convergence thresholds are unchanged.

Column reuse alone did not close the four-worker deadline; that failed run,
including two unchanged MOC module timeouts, is retained at
`history/current_main_post_column_reuse_timeout/manifest.json`. A separate
profile of the unchanged RP-1 density-optimum reference case found 8.8 million
atom-dictionary lookups and millions of generator resumptions. The solver now
prepares dense atom rows/columns once per TP call and uses ordered multiply/sum
iterations. Operand order, sum order and coefficients are retained; no global
cache or cross-provider state is introduced.
For that one instrumented case the complete returned result agreed exactly
before and after the change; profiler elapsed times were 41.03 s and 20.73 s.
Instrumentation and host load affect those times; no general speed ratio or
concurrent-runtime guarantee is inferred.

The new independent check evaluates the normalization and element-ratio
residuals with 100-digit mpmath arithmetic. It checks the initial finite-
difference matrix for three reference-element orderings. Its 2e-7 absolute
matrix-entry tolerance was fixed and passed against the original calculation
before the optimization. It is a test of this binary-input finite-difference
protocol, not a proof of the underlying unrounded thermodynamic model.
The independent residual vector is also checked at an absolute 1e-10 tolerance,
fixed and passed before the dense-iteration change.

Five curated faults exercise repeated columns, transposition, difference sign,
scale and Newton right-hand-side sign. The existing mutation runtime requires passing original controls
before and after each trial; only assertion failures count as kills.
Current source-bound results are in `equilibrium_jacobian_mutations.json` and
`regression/summary.json`. The separate full regression must still complete
with the original collection 30 s / module 120 s / four-worker policy, with
slow comparisons enabled. No case or tolerance is removed to fit the budget.

This scope covers column reuse and its matrix check. It does not qualify all
Gibbs-solver inputs, equilibrium branches, catalogue coefficients, inverse
enclosures, convergence over a domain, or motor/feed hardware.

## References

- [Upstream equilibrium implementation](https://github.com/ayberkdt/Rocket-Propulsion-Lab/blob/fea6a324e9b9e8dc7f88872e6abfddbad8b4842a/src/rocket_propulsion/thermochemistry/equilibrium.py),
  retained with the timeout source snapshot before the optimization.
- NASA NPR 7123.1D, [Appendix E](https://nodis3.gsfc.nasa.gov/displayDir.cfm?Internal_ID=N_PR_7123_001D_&page_name=AppendixE),
  relevant-environment and scaling evidence is separate from numerical checks.
