# Complete area and Mach-angle inverses

Parent: `cc089cf` (the scalar flow-ratio repair, now committed and pushed).
This review addresses area-ratio and Mach-angle inversion. Perfect-gas model
applicability, Prandtl–Meyer accuracy, empirical flow validation and motor/feed
hardware readiness remain distinct requirements.

## Reproduced failures

Five failing assertions were retained before the repair. For A/A*=2 on the
supersonic branch, gamma=nextafter(1,+inf), 1.001 and 1.01 have admitted roots
near 1.9216228894, 1.9222697684 and 1.9281007423. The old solver refused them
because evaluating the unrelated M=1000 bracket endpoint required an
unrepresentable area. This did not make the requested root unrepresentable.

The two smallest positive angle floats were instead passed through binary64
degree-to-radian conversion, rounded to zero and divided into one. They raised
ZeroDivisionError rather than the required typed inverse-range refusal.

An additional boundary trap is independently documented: the supplied angle
0.0572957890623832 degrees implies M=1000.0000000000000184966… . Its float
export is 1000, so checking the bracket after float export falsely admits it.
The numerical bracket [1e-8,1000] is preserved; it is not a gas-model validity
range or a new physical restriction.

## Area calculation and termination

All admitted area ratios use the complete logarithmic relation on the selected
subsonic or supersonic branch. Endpoint comparison uses logarithmic values and
does not require the endpoint's binary64 area or other forward outputs.
The solver uses a safeguarded Decimal Newton step with

```text
L(M) = (gamma+1)/(2(gamma-1)) ln[(2+(gamma-1)M²)/(gamma+1)] − ln M
dL/dM = (M²−1) / [M(1+(gamma-1)M²/2)]
```

The bracket is updated by the branch's monotonic sign. A Newton step outside
that bracket or stalled at its current point is replaced by a geometric or
arithmetic midpoint. Computed zero or bounds that export to the same float end
the iteration. Exhausting the iteration budget raises a typed convergence
failure; no last iterate is exported as a solution. These computed signs and
finite precision are numerical controls, not a rigorous mathematical enclosure.

The unused native area derivative wrapper and the obsolete post-export Mach
bound wrapper are removed. The PM Newton/bisection path is retained separately;
this change does not silently qualify it.

## Mach-angle calculation

The input angle is interpreted as the supplied float in degrees. A private
Decimal context evaluates pi=4 atan(1), a converged sine series and 1/sin(mu).
The complete Mach is compared with the numerical bracket before float export.
Gamma does not enter this angle relation; the enclosing public function still
validates the supplied gamma under its existing contract. Exactly 90 degrees
returns the sonic value 1. No near-90 tolerance is substituted for the equation.

The arctangent and sine iteration failures propagate as typed Mach-angle
convergence errors with `value` ownership. Caller/default Decimal settings and
flags remain separate. No runtime mpmath dependency is introduced.

## Independent controls and falsification

Area references use sonic temperature, entropy and mass-flux conservation,
followed by 140 independent bisection steps at 600 decimal digits. Angle
references use mpmath sine and pi at 1200 digits. Neither inverts production
forward output. The supplied rounded input itself determines the expected
root and bracket classification. Valid outputs have a two-ULP screening budget;
the finite matrix is not a universal rounding guarantee or physical error bar.

Controls cover both branches, wide area ratios, nearly unit/large gamma,
neighboring boundary inputs, degree conversion, gamma independence, near-sonic
angles, typed domain/range failures, absence of private endpoint output
requirements, caller/default context ownership on success and refusal, forced
budget exhaustion, pi-failure propagation and the real calculation route.
Nonrepresentable boundary ratios are explicitly omitted from float-input
generation rather than recorded as skipped tests.

The selected fault set contains two prior invariants and 21 new faults. An
existing area-dispatch fault is explicitly moved from the removed near-sonic
threshold to the complete area dispatch; it restores the old rounded
callback/endpoint algorithm and keeps its independent near-sonic selector.
An existing ratio error-field anchor is expanded because the new angle helper
adds another range error in the same file. Both changes are recorded in
`assurance/area_angle_inverse_mutation_transition.json`. Historical results
remain bound to their original source bytes. A new final run is required for
any current kill claim.

Faults target branch choice, private endpoint requirements, sonic normalization,
native exponent overflow, continuity, bracket direction, safeguards, absolute
tolerance, budget exhaustion, error classification, native radians, degree
scale, root rounding before bounds, missing angle bounds, context ownership,
partial sine, swallowed pi failure, premature sonic clamping and a global
small-angle approximation. Only assertion failures with passing original and
byte-restored controls count as kills. Survivors, crashes, collection errors
and timeouts stay open.

Arithmetic/domain exceptions on independently admitted inputs are explicit
output-contract assertion failures in these tests. A crashed test process,
import or collection failure is not counted as a kill.

## Measured checkpoint and remaining work

The independent module passed 15 methods and 168 subtests. The source-bound
regression passed 22 modules, 286 methods and 12550 subtests, with no failures,
errors or skips. All 23 selected faults were killed by assertion failures, with
46 passing original/restored control runs and byte-identical source/test
restoration. In particular, the float-before-bracket angle fault was caught at
0.0572957890623832 degrees by a missing expected refusal, not a process crash.

The wheel and source distribution built. An isolated installed-wheel smoke,
outside the checkout and without mpmath, checked 11 independently generated
area/angle cases, context ownership and the real calculation route. These are
same-host software checks, not clean-machine or physical hardware validation.
The compact receipt is `assurance/area_angle_inverse_verification.json`; raw
results and the five original failing assertions remain under
`output/assurance/area_angle_inverse/`. Only this report and the compact receipt
are enriched after the measured runs; runtime and test bytes are matched to
their run snapshots.

The retained PM path is still open. An independent 24-case probe first rounds
the supplied angle, then inverts that exact float with a 400-digit reference.
It finds two valid inputs falsely refused at gamma=1000 and 100000 near the
upper Mach bracket. Seven valid outputs exceed a 64-ULP exploratory screen;
that screen is an arithmetic diagnostic, not a physical uncertainty criterion.
For example, gamma=1.4 and nu=130.16759847813037 degrees imply
M≈999.9999999999694, while the retained native callback returns 1000.
These observations include input quantisation and do not blame lost generating
Mach information for callback error. They are retained for the next review.

The motor/feed hardware ledger still has nine unmet evidence gates. The full
TRL 5 objective remains unachieved.

## References

- NACA Report 1135: https://ntrs.nasa.gov/citations/19930091059
- NASA Glenn isentropic relations and area branches:
  https://www.grc.nasa.gov/www/k-12/airplane/isentrop.html
- NIST DLMF sine series: https://dlmf.nist.gov/4.19.E1
- NIST DLMF arctangent series: https://dlmf.nist.gov/4.24.E3
- Decimal contexts and exact float conversion:
  https://docs.python.org/3.12/library/decimal.html
