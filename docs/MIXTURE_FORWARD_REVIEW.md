# Frozen ideal-gas mixture forward contracts

This review covers resolved species identity, mole/mass composition conversion,
forward ideal-mixture properties and owned result records. Species retain the
NASA9 fits reviewed separately. Ideal gases and frozen composition are model
assumptions. This is not a real-gas/chemical-equilibrium model, measurement
validation, inverse-root certificate, motor/feed qualification or TRL5 evidence.
The original complete function/catalogue/claim review and nine hardware gates
remain open.

## Equations and independent comparison

For normalized mole fractions X, species molar masses M_i and the NASA fit R:

```text
M = sum X_i M_i,  Y_i = X_i M_i/M,  R_mix = R/M
cp = sum X_i cp_i / M,  cv = cp - R_mix,  gamma = cp/cv
h = sum X_i h_i / M
s = [sum X_i s_i - R sum X_i ln(X_i) - R ln(P/P_ref)]/M
```

For mass-basis raw fractions, molar amounts are proportional to Y_i/M_i.
The entropy relation is equivalently the mole-weighted sum of standard species
entropies evaluated at each species' partial pressure. P_ref is one bar.
No partial-pressure product or pressure ratio is rounded before its logarithm.
The defining NASA R is retained consistently; mixing does not select a newer R.

Cantera 3.2.0 was installed in a separate ignored toolchain. Its native binary
SHA256 is `00f00749cc90b76a7266656fd02e6a7d40be82e7d7e1389472ad285123064e9d`;
its wheel SHA256 is `258b95a4e29f9afd50788eab284e88440174302123c8b3a4a3a493904b40e238`.
The optional generator `scripts/generate_mixture_cantera_reference.py` uses the
pinned primary NASA extraction, not the production evaluator. It freezes 150
states across mole/mass compositions, pure gases and all ten gases, temperatures
200/298.15/1000/3500/6000 K, and pressures 100/100000/20000000 Pa.

Comparisons use cp/R, h/R and s/R. This removes the different defining gas
constants from the comparison; mass-basis inputs use the historical NASA molar
masses. Cantera's modern atomic weights are not silently substituted into that
conversion. Tolerances were set before running the comparison: relative 2e-12,
absolute cp/R 5e-12, h/R 3e-8 K, s/R 5e-11. They are cross-implementation
screening tolerances, not uncertainty estimates. Source/generator/binary hashes
and vectors are in `assurance/mixture_cantera_reference.json`. The unit suite
uses this frozen file and does not need Cantera installed.

These evaluators share the same coefficients. Their agreement is independent
implementation evidence, not independent physical data or validation of an
engine, real-gas domain, equilibrium chemistry or heat-transfer model.
An additional 180-digit equation oracle checks binary64 forward exports, while
700-digit attacks cover subnormal composition and extreme-pressure cases.

## Composition admission, ownership and migration

Input maps resolve case/whitespace aliases to species identities. Merged onto
current main (2026-10-29), aliases of one resolved species keep main's
contract: their raw amounts are combined exactly into one component before any
mixing term, so H2/h2 cannot invent entropy of mixing; every raw alias entry
stays in `input_fractions`. Direct constructors refuse two components with one
canonical species key.

Positive finite real raw fractions must sum to one within the existing 1e-9
admission tolerance, now evaluated with exact binary-input sums. This tolerance
is not composition uncertainty. Accepted raw fractions and their basis remain
in `input_fractions` and `composition_basis`; their fractions are normalized
exactly once under `exact_normalization_of_retained_fraction_inputs_1e-9_v1`.
That policy is explicit and replay-checked. Inputs outside the tolerance are
refused; arbitrary relative weights are not accepted. The former code admitted
near-one totals while using unnormalized mole fractions in its equations.

Component fraction fields are binary64 displays of exact retained composition
ledgers. Calculations do not reconstruct exact fractions from those displays.
Mandatory nonzero mole/mass fraction displays that genuinely underflow are
refused, never clipped to zero or silently removed. A minimum-subnormal mass
fraction may remain representable while the corresponding mole-basis request
requires an unrepresentable mass fraction; these are different inputs.

Components own reconstructed Species records. Mixtures own reconstructed
components and their raw inputs, with exact deterministic display checks.
Direct construction without explicit raw inputs first requires both display
sets to sum to one within 1e-12 and agree with each other (main's direct
contract), then derives the raw inputs from the declared basis' displays and
re-exports both displays from the exact ledger. Reopening a produced mixture should retain all its
raw inputs and policy; dropping them may define a different composition.
The maximum is 64 species/components. Variable input iteration reads at most
65 items, and fixed name/value pairs at most three. Empty, unordered, malformed
or inconsistent records fail with typed errors.

MixtureProperties now requires its owning Mixture and replays all seven
derived property fields. The owning mixture does not take part in record
equality, and the API reports state fields without it. Scalar fields must be finite, and
physical_model_qualified must be the actual False boolean. Its owned inputs
are JSON-safe; no private Fraction is serialized. Old standalone property
construction without mixture inputs is no longer the result contract.
Defining T/P may describe a different valid state when all derived displays
still replay; a one-ULP change to a derived property cannot bypass replay.

## Numerical scope and retained failures

Species polynomial terms and mixture sums retain exact normalized binary inputs
before separate property exports. Each constituent must have positive retained
cp-R; a positive aggregate cannot hide a constituent with invalid cv. Mixture
properties can remain finite when an individual unweighted property would
overflow; the forward path does not prematurely export that component property.
Positive cv is computed before display subtraction, including when displayed
cp and R_mix round equal.

Log arguments remain rational. Numerators/denominators are split by GCD into
pairwise-coprime integer bases and equal-base coefficients are combined.
This cancels identities such as 3 ln(2)-ln(8) exactly without general integer
factorization. Splitting has a hard 4096-step resource budget. Logarithms use
isolated contexts with explicit rounding, exponent limits, flags and traps.
Each correctly rounded Decimal.ln receives a conservative decimal-ULP bound.
Exact weighted bounds must round to the same binary64 export; precision is
bounded to 48/96/192/384/768/1536/3072 digits. An unresolved rounding fails.
The integer-log cache is bounded to 512 immutable entries.

Mandatory overflow and erased nonzero outputs fail typed; exact zero is allowed.
Exact log-basis cancellation distinguishes zero entropy from numerical
underflow. These are rounding bounds on the retained model, not fit, gas,
measurement or physical-uncertainty bounds. Resource limits are not validated
operating limits. Standard species NASA fits still limit temperature to
200--6000 K; positive finite pressure is an algebraic input domain, not proof
that the ideal-gas assumption is physically accurate there.

The initial unrepaired state and attack suite are retained in
`output/assurance/history/mixture_original_failures/manifest.json` (409 files).
Some initial record attacks required the new ownership API and are not distinct
physical defects. A derived-field forgery test was corrected to exclude T/P,
which are defining inputs rather than derived properties.

A later 86-method connected green run, including 150 Cantera states, still
missed refusal of a valid common endpoint. Fits [200,1000] and [1000,6000]
both permit forward evaluation at exactly 1000 K. That green state and probe
are retained in `history/mixture_green_before_common_endpoint_attack/manifest.json`
(413 files). The new test fails the old constructor; its source and failure are
retained in `history/mixture_common_endpoint_failure/manifest.json` (412 files).
Only the disjoint-range comparison changed from >= to >. A singleton common
temperature supports forward evaluation; it does not supply an inverse bracket.

Twenty-seven meaningful attack methods cover these contracts, and 27 curated
faults target weighting, thermodynamic signs, dependent logs, premature rounding,
typed boundaries, ownership, qualification claims and stored result consistency.
Each fault needs passing original controls before and after its assertion kill.
This is not an exhaustive mutation score. Inverse methods add only strict
target/pressure admission and retain shared bisection, without root enclosures,
global continuity/monotonicity certificates or new accuracy claims.

## Merge onto current main: solver path

Merged on 2026-10-29. Exported mixtures (`Mixture.from_fractions`,
`Mixture.properties`) follow these contracts. Equilibrium results and the
frozen rocket expansion use binary64 helpers (`_binary64_frozen`,
`Mixture._binary64_state`) that reproduce the previous composition and mixing
arithmetic without building owned records; full rocket-performance results
were bit-identical to pre-merge main for three propellant pairs. Equal exact
composition requests share one cached ledger.

## References

- Cantera 3.2, IdealGasPhase, standard and partial molar properties:
  <https://cantera.org/3.2/cxx/d7/dfa/classCantera_1_1IdealGasPhase.html>
- McBride, Zehe and Gordon, NASA/TP-2002-211556, species NASA9 equations,
  standard states and defining gas constant:
  <https://ntrs.nasa.gov/citations/20020085330>
