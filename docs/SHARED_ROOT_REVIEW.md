# Shared bisection: sampled evidence and representability limits

This review applies to `core/numerics.py`, not every inverse solver in the
library. The current source and test digests are recorded by
`scripts/run_bounded_regression.py` and `scripts/run_root_mutations.py`.
Their result files, including failures, are under `output/assurance/`.

## Reproduced problems on upstream fea6a32

The previous adversarial root tests were run against the new upstream code
before transferring the repair. The retained baseline is
`current_main/root_baseline.{txt,xml}`. Its failed assertions include missing
ownership/evidence fields as well as numerical counterexamples; its failure
count must not be described as a count of independent defects.

Concrete numerical problems include acceptance of infinite tolerances and
boolean controls, an expanded admission bracket that can exclude the actual
target, midpoint overflow on a finite same-sign interval, and an acceptance
decision changed by rounding the residual subtraction. The original upstream
tree is retained in `current_main/upstream_fea6a32.zip`.

## Repaired contract

- Normalize finite real controls and observations; reject booleans, text,
  nonfinite values and erased nonzero inputs before doing numerical work.
  The midpoint-evaluation budget is an integer in 1..4096.
- Admit a target only inside the actual observed endpoint range. Exact
  endpoint matches and matching singleton brackets need no interior call.
- Compute the midpoint and residual decision using exact rational arithmetic
  on the normalized binary inputs. The unchanged acceptance rule is
  `abs(observation-target) <= tolerance*max(1,abs(target))`.
- Reject interior observations that contradict observed monotonic ordering.
  An adjacent-float bracket that cannot meet the requested residual raises
  `ConvergenceError(code="root_float_resolution")`. Exhaustion of the
  declared evaluation budget remains `non_convergence`. Neither relaxes the
  tolerance; both retain finite, strict-JSON-safe diagnostics.
  Failed solves report their failure termination and `acceptance_satisfied=false`,
  rather than inheriting the success packet's residual-tolerance label.
- Own the ordered `RootEvaluation` records and replay the requested sample
  coordinates, derived fields and acceptance decision without calling the
  external evaluator again. A coherent invented trace cannot authenticate
  the observations; this is consistency evidence, not a truth attestation.

`RootResult` now requires the input controls, observations and evidence fields
as well as the former scalar fields. Code manually constructing the old
seven-field record must migrate; obtaining results through either solver
function preserves the float-returning wrapper and its residual convention.
`RootEvaluation` is exported from `rocket_propulsion.core`.

## Independent and adversarial checks

`tests/test_root_contract.py` contains exact-rational decision counterexamples,
subnormal and near-maximum finite values, signed affine relations, independently
evaluated nonlinear roots (180-digit mpmath and SciPy Brent), callback failures,
discontinuities, ordered observation ownership and one-ULP record tampering.
The tests compare the owned observations with the actual callback calls.

The root mutation campaign injects 21 curated faults. It uses the existing
live-module-global runtime and requires a passing original control both before
and after each mutated test. Only assertion failures count as kills; execution
errors, skips and failed controls cannot count. This is a targeted assessment,
not an exhaustive fault-space or mutation-score claim.

## What remains open

Finite samples do not establish continuity or global monotonicity. Their
bracket does not certify a root enclosure for an unrounded mathematical model.
A small residual on a flat function does not establish coordinate accuracy.
`root_enclosure_qualified` and `physical_model_qualified` are therefore always
false, and `coordinate_error_bound` must remain `None`.

Thermochemical inverses still need model-specific interval signs, interval
boundary/gap handling, derivative/conditioning evidence and representability
analysis. Changes to this shared utility do not supply those, validate the
new 36-species catalogue, or renew the historical NASA9/mixture review on a
different source tree. Motor/feed hardware and its relevant-environment data
remain unselected and unqualified.

## References

- Python 3.12, [fractions — rational numbers](https://docs.python.org/3.12/library/fractions.html):
  construction from a float preserves that float's exact value.
- NASA NPR 7123.1D, [Appendix E](https://nodis3.gsfc.nasa.gov/displayDir.cfm?Internal_ID=N_PR_7123_001D_&page_name=AppendixE):
  hardware TRL 5 includes relevant-environment demonstration and documented
  scaling, which the sampled numerical protocol does not supply.
