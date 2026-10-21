# Probability, moment and variance numerical review

This is a bounded review of the discrete/hierarchical weighted statistic paths,
conditional-weight export, grouped variance/contributions, sample diagnostics and
their impulse/scale intermediates. It does not validate the supplied scenario
probabilities, uncertainty distributions/copulas, physical engine maps, tail
convergence, native weighted-orbit execution or all other model functions.
All nine motor/feed hardware evidence gates remain open; no TRL is established.

## Reproduced failures and output boundaries

A smallest-subnormal positive scenario probability with four conditional draws
previously became zero in every retained joint weight; variance then raised an
untyped ZeroDivisionError. Other accepted ranges lost weighted products before
summing, overflowed or erased squared dispersion, admitted zero-mass outliers
into extrema, biased moments when total weight was not exactly one, or used an
absolute 1e-15 floor in quantile/relative diagnostics. Qualified scale products
and audit impulses could overflow an intermediate even when the final value
was representable.

Retained joint probabilities are mandatory binary64 outputs. Their exact
binary-input product is rounded once with a checked boundary; a positive product
that cannot be represented fails with a typed DomainError naming the joint
probability. A representable 1e-20 event remains present. No probability floor,
fabricated positive weight or silent scenario deletion is used. This output
contract does not promise that every mathematical positive weight can be stored
in the existing float field.

A smallest-subnormal event with one deterministic draw is also retained when
every required metric variance/contribution remains representable. Tests choose
a sufficient impulse/mass separation to meet that independent output condition;
a representable weight alone does not guarantee every mandatory statistic fits.

Scale components and audit impulse products retain exact rational intermediates.
Scale factors required by the existing artifact-scaling API are checked when
exported, and audit impulses are checked after the complete product. This does
not re-audit every artifact scaling, axial-bound or uncertainty-draw function.

## Weighted distribution and grouped identities

For retained binary weights w_i, W is their exact sum. Moments use w_i/W:
mean = sum(w_i*x_i)/W, variance = sum(w_i*(x_i-mean)^2)/W. Zero-weight points
do not enter the support. Empty, all-zero and negative-weight distributions fail
typed. This is a normalized empirical distribution of the supplied weights;
it is not an estimate of the reliability of those weights.

The existing aggregate admission tolerance, the represented value of 1e-12,
remains unchanged. The actual sum and its exact residual from one are separately
exported with checked rounding. Normalizing moments does not rewrite the stored
weights or invent perfect closure. Discrete and hierarchical zero-thrust
probabilities use their corresponding normalized retained mass.

Within each group, conditional moments normalize by that group's retained mass.
Within-group variance is weighted by group mass/W; between-group variance uses
the squared difference of the group's and overall exact means. Contributions
use the same normalized group weights. Total, within and between variance are
checked separately at their mandatory float boundaries. Fractions use the exact
variance ratios; closure_error reports the exact residual of the three displayed
float variances, rounded once, so a display discrepancy is not hidden as zero.

Weighted variance here is population dispersion of the empirical probability
distribution. Independent replicate sample deviation keeps its separate n-1
denominator. Neither is a measurement uncertainty budget or a certified coverage
interval. The diagnostic relative change is abs(a-b)/max(abs(a),abs(b)), with zero
only for two zero inputs; it has no unit-dependent absolute floor.

## Quantiles and square-root export

The 5/50/95 percentiles use exact rational probabilities 1/20, 1/2 and 19/20.
The first sorted positive-mass value whose cumulative mass reaches p*W is chosen;
there is no absolute probability slack. This is the local empirical inverse-CDF
definition, not interpolation or a claim of tail convergence.

Variance need not be exported before standard deviation. The binary64 square-root
export finds the integer floor of the scaled rational root on its destination
grid, then compares the squared exact midpoint. Midpoint ties select the even
significand, including subnormals and the normal/subnormal boundary. A true
nonzero output that rounds to zero, or an overflowing mandatory output, fails
typed. The independent tests include huge/small dispersion, cancellation,
destination-grid ties and the finite-range boundaries. Exactness is relative to
the supplied binary operands, not to measured physics.

## Persisted policy and evidence

The hierarchical manifest requires statistics_numerics with value
exact_binary_normalized_joint_weights_and_rational_percentiles_v1, included in
its content hash. Legacy payloads without it are rejected explicitly; migrating
an old mathematical policy requires deliberate reconstruction and a new hash.
Derived statistics, decomposition, contributions, aggregate sum/closure and
zero-thrust displays must reproduce deterministically. One output ULP is not
silently discarded under an arbitrary metadata tolerance. Supplied probability
factorization is additionally checked by the retained record's exact product.

Two existing curated negative-diagnostic variants are retargeted to the new
deterministic equality guards: those guards now enforce nonnegative expected
outputs as well as their values. They are not counted as new variants. Historical
evidence and unchanged final-test replay distinguish actual semantic failures
from controls. Decimal tests use an independent 2200-digit normalized ledger;
square-root tie cases use constructed exact rational midpoints. Selected mutant
kills are not exhaustive correctness, hardware validation or certification.

The first probability-stage mutation attempt preserved two surviving older
variants. Its equal-weight oracle now uses positive unequal power-of-two weights
and proves the baseline reopens, so a newer positive-product guard cannot mask
the target policy. The root scalar-normalization variant now uses the direct
programmatic boolean-input oracle; writer-side normalization had independently
fixed the old single-scenario symptom. These are oracle updates, not new kills
added to the variant count or a relaxation of any acceptance limit.

## References

- NIST Dataplot, Weighted Mean, equation 2-20: division by the total weight.
  <https://www.itl.nist.gov/div898/software/dataplot/refman2/ch2/weigmean.pdf>
- MIT RES.6-012, lecture 13.7 derivation transcript: decomposition of total
  variance into conditional dispersion and dispersion of conditional means.
  <https://ocw.mit.edu/courses/res-6-012-introduction-to-probability-spring-2018/c025b4baa48ff7801b96c3915ee62e5c_mHonq7Gjjqg.pdf>
- Python 3.12 fractions: exact rational representation of supplied float values.
  <https://docs.python.org/3.12/library/fractions.html>
- Python 3.12 math: integer square root and binary exponent scaling primitives.
  The local midpoint algorithm is specified above, not endorsed by Python/NIST.
  <https://docs.python.org/3.12/library/math.html>
