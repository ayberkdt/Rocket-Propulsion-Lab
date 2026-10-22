# Blowdown sample summaries and mutation runtime review

This is a bounded arithmetic, retained-record and verification-instrument review.
It does not qualify uncertainty distributions, confidence coverage, a physical
blowdown model or selected motor/feed hardware. The original global review and
all nine hardware evidence gates remain incomplete. Tests are not a TRL certificate.

## Reproduced summary failures

The previous sample summary reported zero standard deviation for [0, 2e-200],
raised an untyped fsum overflow for the representable constant mean of
[1e308, 1e308], erased Pearson correlation for nonconstant tiny vectors, and
used an absolute 1e-15 floor in the stability diagnostic. Ordinary examples
also exposed one-ULP standard-error errors, rounded-center variance errors and
a P05 of 0.30000000000000004 instead of 0.3 for seven consecutive integers.

Sample values now enter an exact binary-input rational ledger. The empirical
mean and n-1 dispersion are accumulated before conversion. Each required mean,
standard deviation and s/sqrt(n) proxy is rounded once at its own boundary.
The n-1 expression is the declared empirical sample formula; dependence between
LHS rows prevents interpreting it automatically as an unbiased IID variance
estimate. No minimum variance, perturbation or manufactured positive number is
introduced. A mandatory nonzero result that cannot be represented raises a
typed error rather than silently becoming zero.

Sample percentiles retain linear interpolation (Hyndman-Fan type 7): with sorted
values, h=(n-1)*p determines adjacent positions and their affine weights. The
named P05/P50/P95 use rational 1/20, 1/2 and 19/20. Affine interpolation avoids an
overflowing span between opposite extremes. These sample percentiles remain a
different convention from the weighted-distribution inverse CDF used by scenario
ensembles; changing this helper does not silently replace either convention.

## Mean interval and inference boundary

The legacy mean_standard_error and mean_ci95 fields remain available, but their
meaning is explicit: the IID formula s/sqrt(n) and mean +/- 1.96*s/sqrt(n).
The conventional decimal multiplier is represented as 49/25. The metadata
`mean_interval_kind=iid_normal_formula_proxy_unqualified_coverage` and
`coverage_established=False` prevent the record from claiming qualified coverage.

NIST's standard mean-interval treatment distinguishes a known population sigma
from an estimated sample deviation and states the assumptions behind normal/t
inference. The engine outputs are nonlinear transforms of declared uncertain
inputs; this review establishes neither normal output distributions nor a
finite-sample coverage guarantee. In particular, one LHS design has dependent
rows: in a single dimension two rows cannot both occupy the same stratum, while
independent rows could. The IID formula therefore is not a design-specific LHS
standard error. The result emits that limitation as a warning. Proper replicated
designs and validated inference remain additional work, not an achieved capability.

Interval endpoints are computed as algebraic mean +/- sqrt(variance/n*(49/25)^2),
without first exporting the mean or the root. Squared rational comparisons locate
the final binary rounding cell and destination-even ties. Explicit overflow
midpoints distinguish a last finite output from infinity. Exact cancellation to
zero is permitted; an erased nonzero endpoint is a typed boundary failure. No
interval endpoint is clipped to a physically convenient support limit.

## Pearson screening and stability

Centered sums and cross-products remain rational. The signed square root of
Sxy^2/(Sxx*Syy) produces the coefficient without overflowing or erasing the
dispersion factors. An exactly constant input or output makes Pearson undefined:
the record retains correlation, magnitude and rank as null, with reason
`constant_input_or_output`. A defined zero coefficient remains a ranked zero.
Defined sensitivities are ordered by magnitude with declared parameter order
breaking ties; undefined entries remain present after the defined entries.
These are linear screening statistics, not causality or global sensitivity.

Half/full mean changes use abs(a-b)/max(abs(a),abs(b)), with zero/zero defined as
zero. Deviation changes use 1-sqrt(min(v1,v2)/max(v1,v2)), with both zero
variances handled explicitly. Raw sample moments determine these changes before
display rounding. Threshold decisions also use exact comparisons: a change can
round to the displayed tolerance while still being strictly above it. The owning
result recomputes the decision from retained samples. This remains an order-based
stability diagnostic, not proof that a probability tail or physical model converged.

## Retained records and configuration

The metric, sensitivity, sample, stability and result records normalize their
owned scalars and collections. Sample metrics must reproduce their artifact;
inputs must be unique, declared and within their bounds. Results require complete
sample indices, metric and metric/parameter coverage, uncertainty-class counts,
supported numerical policies, and reproduced summary/ranking/stability fields.
Storage reordering is canonicalized by semantic identities, including sample
indices. The summary policy is
`exact_binary_sample_moments_linear_rational_quantiles_and_pearson_v1`.

The uncertain nominal must equal the declared deterministic input's binary value;
an absolute tolerance no longer permits an uncertainty distribution centered on
a different configuration. Malformed controls fail before executing a solver.
This internal consistency does not authenticate a campaign, regenerate sample
draws from a seed, rerun the physical solver when opening a record, or establish
the truth of supplied final pressures/source labels. Mocked solver outputs in
assembly tests are explicitly synthetic interface fixtures.

## Fault-injection runtime correction

The previous three mutation runners compiled replacements with copied globals.
A later controlled patch to a dependency could then be ignored. This can produce
a false assertion kill even when the intended faulty branch cannot execute, or
send an injected function into an uncontrolled real solver instead of its fixture.
Historical reports retain that limitation. Prior archive bytes are preserved;
their reported kills are not substituted for a current replay under the repaired
runtime.

Replacement definitions now use a separate locals dictionary and the original
live module globals. Recursive calls, controlled dependencies, cached functions,
properties and class methods retain the intended lookup context. Imported aliases
are resolved after an original-code control and restored after injection, including
aliases created while the fault was active. Numerical caches are cleared between
executions.

Every injected variant requires a passing original-code control before and after
the fault. A failing/skipped/expected-failure control, an invalid execution, or an
unexpected exception invalidates the experiment. Only an actual assertion failure
between passing controls counts as a kill. The required report policy is
`live_module_globals_and_passing_before_after_controls_v1`. The verifier rejects
legacy/missing policies, non-boolean control outcomes, and a claimed kill whose
controls failed. Source-bound original fault IDs are retained; no variant is
removed to improve the reported result.

The infrastructure attacks exercise all three real runners with a controlled
branch that must survive when inactive and fail when active. Additional attacks
cover recursion, cached/property/class-method lookups, failing controls, leaked
state and aliases imported by controls or during injection. Synthetic runtime
fixture reports are isolated in temporary output roots and are not native physics
evidence.

## Evidence limitations and retained failed attempts

Independent Decimal ledgers check moments, interval endpoints and correlations;
exact identities check zero/constant vectors and scale changes. Tests cover both
representable subnormal outputs and genuinely unrepresentable mandatory outputs.
The first post-repair test incorrectly called [0, 2*TINY] wholly representable
even though its P05 is not. A valid subnormal control and a separate typed-failure
attack replace that mistaken assumption. An owner-replay attack was also changed
from a constant-duration mean (rejected by an earlier local guard) to a varying
impulse mean. Production bytes were retained before those oracle corrections.

Local failed-attempt archives include `history/summary_first_oracle_attempt`
(360 files) and `history/summary_first_mutation_attempt` (359 files), under
`output/assurance`. The latter records 325 reported kills and two invalid results
from the old copied-global runtime, not a passing final experiment. Setup-error
logs are retained separately and do not count as product defects or fault kills.
Current full regression, historical replays and controlled fault runs must be
evaluated against their own source hashes before closing this checkpoint.

Remaining work includes calibrated/design-specific inference, extreme latent
Gaussian behavior, other sampling/study APIs, the rest of the source/catalogue
review, and motor/feed experimental evidence. None is inferred from these checks.

## References

- NIST, confidence limits for the mean and their interpretation:
  <https://www.itl.nist.gov/div898/handbook/eda/section3/eda352.htm>.
- NIST, known-sigma and estimated-sigma confidence-interval distinction:
  <https://www.itl.nist.gov/div898/handbook/prc/section1/prc14.htm>.
- NumPy 2.2, declared linear sample-quantile convention (not a runtime dependency):
  <https://numpy.org/doc/2.2/reference/generated/numpy.quantile.html>.
- NIST Dataplot, Pearson centered-product coefficient:
  <https://www.itl.nist.gov/div898/software/dataplot/refman2/auxillar/correlat.htm>.
- SciPy 1.14.1, constant-vector Pearson is undefined (not a runtime dependency):
  <https://docs.scipy.org/doc/scipy-1.14.1/reference/generated/scipy.stats.pearsonr.html>.
