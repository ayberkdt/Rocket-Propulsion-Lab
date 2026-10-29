# Frozen gas inputs and Bartz arithmetic: scoped review

Baseline: `191a9cd497ed71e3a627d72994fe742f60f022a6`. The reviewed implementation is
`heat_transfer/regen_chamber.py`: `HotGas.__post_init__`, `bartz_coefficient`, its two
arithmetic helpers, and pressure/radius admission in `hot_gas_from_reactants`.
This is input-contract and numerical evidence, not a review of the complete cooling
solver, a new physical validity envelope, or hardware TRL evidence.

## Reproduced defects

- `HotGas` had no constructor validation. A boolean specific heat produced a positive
  film coefficient, and H2O/CO2 fractions of 0.8 each were accepted.
- The public station function accepted a negative wall temperature. With the ordinary
  fixture it returned a positive coefficient around 32786 W/m²/K for −10 K.
- A finite pressure/c* ratio beyond binary64 range returned `inf`; very large finite
  Mach returned a zero coefficient and `NaN` recovery temperature, although both
  complete expressions can have finite positive results.
- Ratios and partial products in the subnormal range could lose significant digits
  before later factors restored normal magnitude. Checking only final finiteness
  would preserve these wrong finite answers.
- Recovering the kinetic term as `C−1` erased a small Mach contribution before a large
  recovery factor amplified it. The corrected expression retains `(γ−1)M²/2` directly.
- Invalid throat curvature could reach expensive chemistry before failing, or be
  stored as an invalid gas state.

The first test run was deliberately made against the unchanged implementation;
the before XML records refusal, arithmetic and ownership failures. Later tests
also cover the small kinetic contribution and caller Decimal context.

## Corrected contract

Frozen scalar fields are owned native finite floats, with the shared admission
policy refusing booleans, strings, nonfinite values and erased nonzero conversions.
Pressure, temperature, c*, cp, viscosity, Prandtl number and curvature are positive;
gamma exceeds one. Each recorded radiating-species fraction lies in [0,1], and
their floating-point sum cannot exceed one. These are necessary domain checks,
not thermochemical consistency or a validated combustion state.

Station radius and wall temperature are positive, Mach is nonnegative and geometric
A/A* is at least one. Zero Mach remains valid and recovers the stagnation temperature.
The reactant factory admits pressure and curvature before invoking chemistry.

Ordinary values use the existing product order. Nonfinite, zero or subnormal ratio/
product intermediates trigger evaluation of the complete positive expressions in
a private 90-digit Decimal context. Only the complete outputs are converted back
to float. Unrepresentable positive outputs raise `DomainError(code='numeric_range')`;
they are not exported as zero or infinity. The private context is independent of
the caller's active Decimal precision, exponent bounds and traps.

## Verification that matters

`tests/test_hot_gas_contract.py` covers:

- every gas scalar's malformed values, physical domains and fraction closure;
- native ownership and a real coupled coolant solution invariant to later mutation
  of a caller-owned specific-heat object;
- normal operating points and synthetic range attacks against a separately written
  100-digit mpmath log-space correlation and cube-root recovery expression;
- overflowing ratios/diameter, subnormal ratios/products, large Mach, the small
  kinetic contribution, and unrepresentable final outputs;
- failure before chemistry for invalid factory curvature.

The fault inventory includes deliberate removal of ownership/domain checks, rounded
subnormal intermediates, zero fallback, reduced precision, wrong pressure exponent,
cancelled recovery term, unchecked output and late factory admission. The mutation
runner requires a passing original control, an assertion failure on the fault and
a passing byte-restored control; crashes, skips and surviving faults remain open.

The measured Python 3.12 regression run passed 71 test methods and 722 subtests.
All 12 selected faults were killed, with 24 passing original/restored controls.
These counts describe this scoped run, not the entire package or hardware readiness.

Related regressions retain real-fluid coolant marching, variable geometry, fin/2-D/
axial conduction, API admission order and the published NASA hot-fire comparisons.
Those hot-fire comparisons' existing discrepancies remain model limitations; this
change does not tune Bartz to those measurements or remove their uncertainty.

## Limits and next evidence

Extreme synthetic cp, Pr, Mach and dimensions are arithmetic stress tests, not
physically validated operating points. Numerical representability is independent
of correlation validity. There is no newly certified pressure/temperature range,
burn duration, chamber lifetime or safe wall margin from this repair.
The actual motor/feed configuration, mission requirements, relevant environment,
instrumented brassboard, uncertainty, holdout comparison and scaling evidence are
still needed; see [TRL_ASSESSMENT.md](TRL_ASSESSMENT.md).
