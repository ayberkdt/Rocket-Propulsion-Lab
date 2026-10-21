# Bounded marginal and sampler review in progress

Follow-up: BLOWDOWN_SUMMARY_AND_MUTATION_RUNTIME_REVIEW.md records the subsequent
blowdown summary checks and the copied-global mutation-runtime limitation.
The earlier evidence retains its own source bytes and historical scope.

This review addresses finite-range inverse transforms, typed uncertainty and
correlation inputs, binary stratification, latent Gaussian transport, and sampler
identity in hierarchical ensemble documents. It remains a partial review. The
precision budget is finite and arbitrary extreme latent tails have not been
qualified. Blowdown summary statistics/sensitivity and other
uncertainty/scaling paths still require separate numerical review. All nine
motor/feed hardware evidence gates remain unmet. No TRL qualification or NASA
endorsement follows from these software checks.

## Counterexamples and output contract

Finite uniform bounds [-MAX, MAX] previously overflowed their span, returning
infinity at the median instead of zero. Triangular inverse products overflowed
or underflowed before the square root. Log-uniform log/exp round trips changed
representable geometric medians. Truncated normal CDF subtraction lost narrow
bounded intervals. Correlation entries within a tolerance were silently ignored,
and an absolute Cholesky pivot floor rejected a genuinely positive definite
binary matrix with rho = nextafter(1, 0).

Mandatory marginal outputs retain the existing checked binary64 convention:
an exact zero is permitted, while a nonzero mathematical quantile that rounds
to zero raises DomainError with field `Marginal value`. An unresolved enclosure
raises the distinct ConvergenceError with field `Marginal enclosure`; that is
not accepted as proof that an output is unrepresentable. Tests distinguish both.
The finite precision schedule is a numerical limit, not a claim that every
accepted distribution can be inverted within that budget.

## Inversion and enclosure argument

Uniform affine combinations retain exact rational intermediates. Triangular
CDF comparisons use squared rational distances, so neither the span product nor
its square root is first exported as a float. The ordered binary64 search checks
the final rounding cell and uses destination-even midpoint ties, including
subnormals. Log-uniform inversion compares logarithmic CDF enclosures, with an
exact rational power comparison for the small-denominator branch. Search hints
cannot establish a result; they must pass the same cell comparison.

For Gaussian primitives G(x) = integral from 0 to x of exp(-t*t/2) dt, the
implementation uses the positive series obtained by rescaling DLMF 7.6.2:
G(x) = exp(-x*x/2) sum over n >= 0 of x**(2*n+1)/(2*n+1)!!.
The successive term ratio is x*x/(2*n+3). After that ratio falls below one,
subsequent ratios decrease, so the first omitted term divided by one minus the
current ratio bounds the entire positive remainder. Directed Decimal arithmetic
encloses the sum, products and remainder. Correctly rounded Decimal exp/ln
endpoint evaluations are expanded by an adjacent Decimal value before use.
The calculation uses explicit contexts, independent of a caller's Decimal context.

The Gaussian normalizing integral uses R = isqrt(6*(precision+14))+1.
Integration by parts gives tail(R) <= exp(-R*R/2)/R, which is less than
10**(-precision-12) at these positive radii. This bound encloses discarded
integration work; it does not clip an input probability. Truncated normal
mass is assembled from centered primitives rather than subtracting two rounded
float CDFs near one half. Cancellation can still require more precision than
the configured schedule, particularly for extreme latent Gaussian targets.

The first implementation used an alternating series with large cancelling
intermediates. Two unchanged 120-second targeted runs timed out. Their complete
353-input source snapshot and logs are preserved locally in
`output/assurance/history/marginal_first_timeout/manifest.json` (357 files).
Replacing that series with the positive enclosure produced a 48-test passing
targeted run at the same budget. Later identity/input changes have their own
source-bound runs; this earlier pass does not qualify those later bytes.

## Sampler and identity boundaries

Random zero endpoints are resampled; a repeatedly broken source terminates with
a distinct convergence error after at most 1024 calls. No arbitrary positive
probability floor is introduced. Latin-hypercube arithmetic forms (i+u)/N exactly
and directs its final binary quantization into [i/N, (i+1)/N) and the global open
unit interval. Its intermediate may round to zero before that directed step.
The smallest positive jitter therefore has a valid first cell even at N=10000.

Correlation matrices are immutable snapshots with exact unit diagonal and
symmetry. Exact rational LDL-transpose decisions determine positive definiteness;
only final Cholesky coefficients are rounded. Copula products/sums retain rational
intermediates and export one latent variate. Built-in marginals receive that
latent Gaussian value, avoiding a premature float CDF of zero or one. Custom
callbacks retain their public float-probability interface and fail typed if
an open probability cannot be represented. This does not certify a user callback.

Independent copula normal inverses also certify their binary output cells.
NormalDist supplies a hint only: at the smallest positive probability and at
nextafter(1,0), its output differed by one ULP from the independent references.
The unit-normal search interval [-64,64] encloses all open binary64 probability
quantiles by the Mills tail bound; it is a search bracket, not a clipping rule.
When a hint misses its first rounding cell, up to three adjacent cells are also
certified before the full ordered-grid fallback. Each candidate needs the same
directed comparisons; there is no assumed accuracy bound on the hint. An attack
supplies both nearby incorrect hints and deliberately distant endpoint hints.

Hierarchical documents require the hash-bound `sampling_numerics` marker
`open_binary_strata_latent_gaussian_and_certified_bounded_quantiles_v1`.
Missing or unsupported markers require an explicit migration and are rejected.
The same marker enters conditional artifact provenance, so a sampler policy
change changes the conditional artifact identity even when numeric draws coincide.
Blowdown uncertainty results expose the sampler marker. A declared policy,
content hash, source string or seed does not authenticate measurement data or
prove that arbitrary caller-supplied records were generated by that sampler.

## Independent evidence and remaining work

`scripts/generate_marginal_references.py` reproduces the 22 exported reference
cases using mpmath 1.3.0 at two decimal precisions per case. It is an offline
generator, not a runtime or test dependency. The fixture covers full-range
log-uniform inputs, bounded/shifted normal inputs, smallest positive probability,
and a normal span [0, smallest-subnormal] with sigma=MAX. Exact binary input
values are retained as hexadecimal strings. Agreement at two high precisions
is finite reference evidence, not a universal correct-rounding proof.

Further attacks cover independent Decimal triangular ledgers, subnormal midpoint
ties, real near-singular SPD matrices, ignored matrix fields, caller mutation,
malformed controls, zero RNG endpoints, both ends of Latin-hypercube cells,
copula clipping, endpoint CDF loss, and changed sampler identity. Curated injected
faults count as killed only by assertions; exceptions and survivors remain open.
Full current-source regression, historical replay, curated faults and remaining
extreme-tail concerns must be assessed before sealing this review stage.

The first curated run produced 294 assertion kills and two survivors among 296
variants. Its 355 input files and seven reports/logs are retained in
`output/assurance/history/marginal_first_mutation_attempt/manifest.json`.
The diagonal attack used an above-one coefficient that another guard already
rejected; a below-one near-unit diagonal was added to exercise the unit-diagonal
rule itself. The sampling-policy attack now checks the typed policy field so a
later generic hash rejection cannot mask removal of the version guard. No old
variant was deleted or counted as killed in that failed run.

The first full regression at the unchanged 30-second collection/120-second
module budgets and four workers passed 92 of 93 modules, with
`test_burn_uncertainty.py` timing out. Its complete 552-file source/evidence
snapshot is retained in `output/assurance/history/marginal_first_full_timeout`.
Diagnostics identified full-grid Gaussian searches after a one-ULP hint miss,
alongside the unchanged blowdown workload. The local-cell search change retains
the same comparisons and full fallback; subsequent runs must establish its
effect under the original budgets. A timed-out module is never a pass.
Nonnegative interval products also use monotonic endpoint multiplication,
retaining outward rounding while avoiding the four-corner calculation needed
for general signed intervals. The independent reference values remain unchanged.

## References

- NIST DLMF 7.6.2, positive error-function series:
  <https://dlmf.nist.gov/7.6#E2>.
- NIST DLMF 7.8, Gaussian tail/Mills bounds:
  <https://dlmf.nist.gov/7.8>.
- NIST Dataplot TRIPPF, triangular inverse distribution formulas:
  <https://www.itl.nist.gov/div898/software/dataplot/refman2/ch8/trippf.pdf>.
- SciPy 1.14.1, log-uniform CDF and quantile:
  <https://docs.scipy.org/doc/scipy-1.14.1/tutorial/stats/continuous_loguniform.html>.
- Python 3.12 Decimal, explicit directed contexts and correctly rounded exp/ln:
  <https://docs.python.org/3.12/library/decimal.html>.
