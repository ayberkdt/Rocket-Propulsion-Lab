# Inverting the supplied isentropic flow ratio

Parent: `504dbbcc3a35051d267134aec4fd41a33efce7a3`.
This review concerns six static/total and static/sonic temperature, pressure and
density inverses. It does not claim that every angle/area inverse, real-gas
state or manufactured motor/feed system is qualified.

## Requirements and reproduced defects

Interpret the supplied binary64 value itself. Return a finite Mach within the
existing derived-inverse bracket [1e-8,1000], or refuse with the correct typed
physical-domain / inverse-range error and `value` field. A rounded generating
Mach is not an oracle: a forward plateau can discard that information.

Four independently confirmed counterexamples were recorded by the forward
review. At gamma=1e308, supplied T/T0=5e-309 and T/T*=0.25 imply M≈2; at
gamma=nextafter(1,+inf), supplied T/T*=1 and p/p*=1 imply exactly M=1.
The old implementation rejected all four as outside [1e-8,1000]. The new test's
initial red JUnit records those assertion failures before the repair.

The failure mechanisms were a private 700 logarithm cutoff before a balanced
division, early rounding of `(gamma+1)/2`, and an implicit promise that critical
conversion was exact. The public docstring also conflated direct Mach identity
with derived-inverse bounds. These claims are corrected, without restricting the
forward input domain or silently clamping a root into the inverse bracket.

## Complete calculation

With supplied ratio r, let a be 1 for temperature, (gamma-1)/gamma for pressure,
or gamma-1 for density. Set L=-a ln(r), adding ln((gamma+1)/2) for a critical
ratio. Complete Decimal expressions retain both gamma±1 terms before these
operations. L<0 is a physical-domain refusal; otherwise compare L with
ln(1+(gamma-1)Mmin²/2) and ln(1+(gamma-1)Mmax²/2). A root outside that numerical
bracket is refused before computing exp(L). An admitted root is
sqrt(2(exp(L)-1)/(gamma-1)), exported once after the complete evaluation.

The context owns precision, rounding, exponent limits, traps and flags.
Caller contexts and mutable Decimal defaults cannot select solver behavior.
No runtime mpmath dependency is introduced. Ratio inversion does not require
unrelated forward outputs to be representable: near-unit gamma can have an
admitted pressure-ratio inverse even when the corresponding complete forward
state cannot export its area ratio.

## Independent controls and falsification

The 1200-digit mpmath reference first obtains temperature from entropy and the
sonic temperature reference, then uses total enthalpy to recover velocity/sound
speed. It does not invoke production forward code, invert a production result
or use the production logarithmic bracket algorithm. Selected valid Mach outputs
have a two-ULP screening budget; this finite sample is not a rigorous global
rounding enclosure or measured physical uncertainty.

Controls include exact-sonic critical ratios across gamma scales, raw supplied
ratios down to the smallest positive float, six quantity/reference combinations,
neighboring floats at both inverse boundaries, distinct pressure/density powers,
physical-domain classification, finite roots outside the bracket, rejection
before enormous private powers, caller/default context independence and the real
HTTP route. Nonrepresentable boundary ratios are not supplied as fake float
inputs; the loops explicitly omit them rather than recording skipped tests.

Seventeen selected faults target gamma forwarding, sonic normalization, entropy
powers, log sign, domain/range classification, both bracket ends, premature
energy export, private context ownership, precision, eager exponentiation,
computed error fields, exact-float ownership and rounding before bracket tests.
Only assertion failures with passing original/restored controls and byte-exact
restoration count as kills. A survivor, timeout, crash or collection error stays
open. The selected set is not exhaustive mutation coverage.

## Remaining scope

Angle inversion and numerical Newton/area brackets remain separately reviewable.
This repair does not make all forward-rounded properties identify the original
Mach, provide an all-input mathematical enclosure, establish perfect-gas
applicability or validate a relevant hardware environment. The motor/feed design,
requirements, test article and empirical evidence are still required for TRL 5.

Read-only follow-up checks found three valid `area_ratio_supersonic=2` cases
(gamma=nextafter(1,+inf), 1.001 and 1.01) refused because the unrelated M=1000
endpoint cannot export its area. Independent mass-flux roots are approximately
1.9216228894, 1.9222697684 and 1.9281007423. Also, the two smallest positive
Mach-angle inputs produce an untyped ZeroDivisionError instead of the expected
inverse-range refusal. One ordinary area case and one larger tiny-angle case
remain passing controls. Those five defects were open at this scalar-ratio
checkpoint; their subsequent repair is tracked in
[AREA_ANGLE_INVERSE_REVIEW.md](AREA_ANGLE_INVERSE_REVIEW.md). The original scalar
receipts remain historical evidence for their archived bytes.

## Measured checkpoint

The new independent module passed **14 methods and 1312 subtests**. The final
regression passed **21 modules, 271 methods and 12382 subtests**, with no failures,
errors or skips. The 17 selected faults all failed by assertions; 34 original
and byte-restored control runs passed. The float-before-bracket fault was caught
specifically at T/T0=1.999960000799982e-05, gamma=1.1: the complete root is
1000.00000000000004676…, while its rounded float is 1000. This is a boundary
classification defect, not an eager-power crash counted under another name.

The wheel and source distribution built from a copied source snapshot. An
installed-wheel smoke, outside the checkout and without mpmath, checked eight
independently generated expected cases, context ownership and the real HTTP
route. Ruff and mypy passed (207 source files). These are software checks on
the same host; no clean-machine or hardware qualification is claimed.

Compact source-bound evidence is in
`assurance/isentropic_ratio_inverse_verification.json`; raw receipts and preserved
initial failures are in `output/assurance/isentropic_ratio_inverse/`. Only this
report and the compact receipt are enriched after the runs; executable and test
bytes remain matched to their snapshots. The hardware ledger still has nine
unmet gates, so the user's complete TRL 5 objective remains unachieved.

## References

- NACA Report 1135: https://ntrs.nasa.gov/citations/19930091059
- NASA Glenn isentropic equations:
  https://www.grc.nasa.gov/www/k-12/airplane/isentrop.html
- Python Decimal exact float conversion and context ownership:
  https://docs.python.org/3.12/library/decimal.html
