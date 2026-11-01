# Complete forward isentropic ratios

Arithmetic baseline: `a799830519ad0dc0a0d6be2c78e87d704e3e1426`.
Current parent: `f62ad44` (intervening Sentinel research documents).
The admitted mathematical inputs are finite M>0 and gamma>1. This audit does
not treat extreme synthetic gamma as a physically realizable perfect gas.

## Reproduced failures

- M=2, gamma=1e308: an intermediate energy factor overflowed, although all
  requested ratios fit. T/T0 and p/p0 round to 5e-309, density ratio to 1 and
  T/T*, p/p* to 0.25. The old implementation refused the whole finite case.
- M=1e-160, gamma=1e308: `2*(gamma-1)` overflowed inside the area exponent.
  A/A* was exported as 1e160 instead of approximately 1414213.56237345.
  This was a successful-looking finite answer, not a detected overflow.
- M=10, gamma=nextafter(1,+inf): rounding `gamma+1` and sonic normalization
  before applying the large exponent yielded A/A*≈5.18470553e20 instead of
  ≈3.14468286e20. The static/sonic pressure and density were also wrong.
- A normalized power can overflow before division by M even when the complete
  area fits, e.g. M=3e15, gamma=1.1. Area-only evaluation must not require a
  representable unrelated stagnation-pressure output.

The initial probe is retained at
`output/assurance/isentropic_forward_range/initial_probe.json`. These are
counterexamples to implementation arithmetic, not new physical range limits.

## Repair and output contract

The ordinary well-conditioned path retains its native arithmetic. A conservative
dispatch uses complete Decimal expressions for gamma<1.1, gamma>1e6,
M<1e-80 or M>1000. Within that bounded native candidate range it also switches
when `max(abs(log(X)), abs(log((gamma+1)/(2*X))))*gamma/(gamma-1) > 16`.
This screens amplification of rounded power exponents. Invalid scalar exports
also trigger checked complete evaluation. These dispatch thresholds select an
algorithm; they do not reject input or define gas-model applicability.

The first 176-case scale sweep found 96–112 ULP pressure-ratio errors at
M=1e40, gamma=4 despite passing the first fixed grid. The original scale receipt
and a new failing large-power test are retained. The screening bound was not
relaxed: large-power cases have a two-ULP control. The previous native overflow
catch is no longer a reachable admitted path after complete dispatch; its fault
was replaced by a direct amplified-power-dispatch fault rather than counted
as an untested successful kill.

The complete expressions retain X=1+(gamma-1)M²/2, the sonic factor and both
gamma±1 terms before logarithms/powers. Static/total and static/sonic ratios
are exported independently with their computed field names. Positive subnormals
remain valid. A true result that rounds to zero or infinity raises a typed
`numeric_range` error; successful output never silently contains those ratios.
Area-only evaluation exports only area, rather than evaluating unrelated fields.

Each context explicitly owns precision, exponent limits, rounding, traps and
flags. Precision is conservatively increased for input scales. Active caller
contexts and mutable `decimal.DefaultContext` cannot choose solver behavior.
No process-global Decimal setting changes. No mpmath runtime dependency is added.

## Independent tests and falsification

`tests/test_isentropic_forward_range.py` builds a 1200-digit reference using
static/total temperature, entropy and continuity at the sonic station. Its area
comes from mass-flux conservation rather than the production normalized-power
expression. Complete-expression selected cases use a two-ULP screening bound;
the ordinary grid uses a declared 64-ULP arithmetic screening budget. These
finite grids are not rigorous all-input rounding enclosures or empirical error
bars for a real flow.

Controls include the smallest positive Mach paired with high gamma, retained
sonic normalization as gamma approaches one, subnormal rounding, genuinely
unrepresentable outputs, separate area evaluation, poisoned caller/default
contexts, strict JSON and preservation of the ordinary native path.

Fourteen explicit mutations remove complete-expression or amplified-power dispatch, reintroduce
the private gamma overflow, round the sonic reference early, corrupt pressure
or density exponents, export private infinity/false zero, reject subnormals,
inherit caller/default contexts, reduce precision, lose computed error fields,
or require unrelated outputs for area. Originals, assertion-only mutant
failures and byte-restored passing controls must be inspected before any kill
claim. Import/collection failures, skips, timeouts and surviving mutants remain
OPEN. The fixture is a selected falsification set, not exhaustive coverage.

The first isolated-wheel smoke also failed because its hand-entered near-unit
gamma reference was mistyped. The failed trace is retained. The corrected
installed check uses the independent oracle's rounded value and a two-ULP bound;
it is not a tolerance relaxation or evidence that the original smoke passed.

## Remaining scope

The broad inverse documentation no longer promises universal machine-precision
recovery: conditioning, forward rounding plateaus and different inverse
algorithms need separate review. This forward repair does not establish every
inverse, shock/real-gas relation, nozzle discretization, continuous flow or
manufactured hardware correct. Relevant-environment motor/feed evidence is
still required for the user's TRL 5 target.

Four inverse false refusals are confirmed independently by inverting the exact
supplied floats, rather than a forward/inverse round trip. At gamma=1e308,
T/T0=5e-309 and T/T*=0.25 each imply M≈2. At gamma=nextafter(1,+inf),
T/T*=1 and p/p*=1 each imply exactly M=1. The current inverse rejects all four
as outside [1e-8,1000]. Its 700 log cutoff and rounded sonic normalization
remain open defects. The docstring's claim of "exact" critical-ratio conversion
is also unproven and contradicted by these cases; it is not endorsed by this
forward review. These are synthetic arithmetic cases, not a real-gas envelope.

## Measured checkpoint

The final source-bound regression covered **20 modules, 257 test methods and
11070 subtests**, with no failures, errors or skips. The independent forward
module contributed 11 methods / 885 subtests. The 176-case scale matrix had
157 finite successful cases and 19 typed unrepresentable cases, no unexpected
result, and a largest observed ratio discrepancy of 11 ULP across its mixed
two-/64-ULP controls. That observation is not an all-input error guarantee.

All 14 final selected faults were killed by assertion failures, with 28 passing
original/restored control runs and byte-identical restored source/test snapshots.
The two within-review fault changes are recorded in
`assurance/isentropic_forward_range_mutation_review.json`; the retired unreachable
catch is not counted as a final kill. The wheel and source distribution built,
and an isolated no-dependency wheel smoke verified complete ratios, subnormal
exports, caller-context independence and the real HTTP route.

The compact, source-hashed receipt is
`assurance/isentropic_forward_range_verification.json`; raw results and preserved
failures are under `output/assurance/isentropic_forward_range/`. Only this report
and the compact receipt are enriched after the measured runs; executable and test
bytes are checked against those run snapshots. The hardware evidence audit still
has nine unmet gates. The full objective remains unachieved.

## References

- NACA Report 1135: https://ntrs.nasa.gov/citations/19930091059
- NASA Glenn isentropic flow equations:
  https://www.grc.nasa.gov/www/k-12/airplane/isentrop.html
- Python Decimal contexts and exact float conversion:
  https://docs.python.org/3.12/library/decimal.html
