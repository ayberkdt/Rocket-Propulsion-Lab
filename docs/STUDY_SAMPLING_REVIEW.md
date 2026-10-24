# Generic study sampling and result review

This review covers `studies.sampling`, the delta-v uncertainty route and shared
line-chart coordinates. It is software implementation evidence, not a physical
model, inference-coverage or motor/feed qualification. The global review remains
partial and all nine hardware gates remain unmet.

## Reproduced failures and arithmetic

The previous sampler overflowed finite `[-MAX, MAX]` spans. Summaries overflowed
the representable constant mean of `1e308`, erased tiny dispersion and Pearson,
and rounded five-percent positions. The route rounded large integer seeds through
float, and `g0 * Isp * log(ratio)` overflowed before a finite final delta-v.

Uniform transport now computes `(1-q)*lower + q*upper` from exact binary input
fractions before conversion. LHS retains exact `(i+u)/n`, the existing RNG draw
and shuffle order. Monte Carlo retains row/parameter RNG order. The policy is
`exact_binary_affine_uniform_latent_strata_v1`. Strata are half-open in the latent
variable; the finite output grid may merge draws or round to endpoints. Adjacent
zero/subnormal bounds remain supported: nearest rounding to zero can be a
legitimate sample draw. Inclusive sweeps use exact `i/(n-1)`.

Empirical mean, n-1 standard deviation, linear sample percentiles and Pearson
reuse the exact kernels in `BLOWDOWN_SUMMARY_AND_MUTATION_RUNTIME_REVIEW.md`.
Named P05/P50/P95 use rational 1/20, 1/2 and 19/20 positions. A mandatory nonzero
summary that cannot be represented fails with a named DomainError; exact zero
remains valid. Constant input/output yields null Pearson/magnitude/rank and a
reason. Defined zero remains ranked; stable ties retain declared parameter order.
The n-1 expression is an empirical convention, not an unbiasedness claim for
dependent LHS rows. Pearson screens linear association, not global sensitivity.
Sample percentiles are not confidence intervals.

## Ownership and API

Ranges, outputs, counts, seeds, methods and record collections have strict typed
boundaries. Bool/numeric-text callback outputs are rejected; unrelated callback
exceptions propagate. Parameters and fixed inputs are snapshotted. Callbacks
receive immutable owned scalar mappings. `asdict`, JSON, deepcopy and pickle
remain usable. Deliberate low-level bypasses of frozen Python objects are outside
this ownership contract.

The owning result retains parameter ranges and numerical policies, checks sample
cardinality/keys/bounds, and recomputes every summary and sensitivity. It does
not regenerate draws from the claimed seed, authenticate a campaign, rerun an
evaluator or qualify physics. The declared 2..100000 samples and 2..10000 sweep
points remain available. Maximum-count generation is exercised, not qualified
for concurrent-service performance or every evaluator/parameter count.

Uncertainty API metadata advances to 2.0. JSON integer seeds and integral form
strings retain their exact value. The entire Isp range must be positive and the
mass-ratio range at least one before sampling. Delta-v uses a checked exact
product of binary g0, Isp and the platform `math.log` result; logarithm enclosure
is not established. Clients must handle null Pearson and added range/policy/
warning fields. The UI labels sample percentiles and undefined sensitivities;
workspace export retains numerical policies and inference warnings.

## Chart and evidence scope

Four instrumented executions of the previous chart failed: opposite finite
extremes, constant MAX, visible tiny dispersion, and a NaN gap. Axes now normalize
before differences/padding, avoid argument-spread extrema and break invalid
paths. Approximate display ticks remain finite; model data are not clipped.
The Node standard-library harness executes actual chart/result functions with an
instrumented DOM, including 100000 points and literal parameter names. It does
not prove browser layout or every frontend flow.

Independent Decimal moments, exact affine draw oracles, known rational sample
percentiles, public API attacks, record forgery and controlled mutations support
the bounded review. Source-bound baseline/later reports remain locally under
ignored `output/assurance/`. Subtest rows, expected future-field errors and bare
historical exceptions are not unique methods or assertion kills. Design-specific
inference, extreme latent tails, global review and hardware evidence remain open.

Development fault attempts retain an equivalent bool-count coercion (0/1 still
fails the minimum count) and an uninjected handler replacement: the ROUTES table
already held the original callable. Neither is an assertion kill. The unchanged
dispatch test now attacks the identical physical guard extracted into a live
module helper; the retained handler resolves this helper normally. The runtime
does not claim to rewrite arbitrary callback registries.

## References

- Python 3.12, exact binary float-to-rational conversion:
  <https://docs.python.org/3.12/library/fractions.html>.
- NumPy 2.2, linear/type-7 sample quantiles:
  <https://numpy.org/doc/2.2/reference/generated/numpy.quantile.html>.
- SciPy 1.14.1, Pearson and constant inputs:
  <https://docs.scipy.org/doc/scipy-1.14.1/reference/generated/scipy.stats.pearsonr.html>.
