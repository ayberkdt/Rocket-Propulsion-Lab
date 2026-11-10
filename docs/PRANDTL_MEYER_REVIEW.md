# Complete Prandtl–Meyer arithmetic and inverse admission

Parent: `fa037cf`. This is a numerical review of the calorically perfect-gas PM
forward and inverse calculations. It does not establish a real-gas operating
envelope, empirical expansion-fan accuracy, or motor/feed hardware readiness.

## Observed failures before repair

The exact supplied angles 0.08993033878768658 degrees at gamma=1000 and
0.0008988585731864637 degrees at gamma=100000 have independently calculated
roots below M=1000. Both were refused by the old native inverse. Other supplied
angles resolved to errors up to 129932 binary64 ULP in the exploratory matrix.
ULP measures arithmetic error here, not experimental uncertainty.

The pre-repair five-method independent test run retained 35 failing subtests.
It covers false refusals, quantized callback plateaus, endpoint classification,
forward cancellation and missing error fields. The historical raw observations
are retained in the preceding area/angle verification receipt; their hashes
identify the older source. They are not measurements of the repaired code.

## Calculation and admission

All forward PM angles now use a complete Decimal evaluation, including degree
conversion with pi=4 atan(1). The existing factorized near-sonic series is kept;
the general expression uses high-precision arctangents. Precision increases with
gamma's decimal exponent to retain the difference of nearly equal terms. Pi is
cached as an immutable Decimal keyed by precision, with at most 32 entries.
Cache misses explicitly own precision, rounding, exponent range, traps and flags.
Neither mutable Decimal defaults nor the caller's context configures that work.

All positive inverse PM angles use this complete expression. The physical
limit is evaluated without subtracting nearly equal square roots:

```text
nu_infinity = 180 / [(gamma-1)(sqrt((gamma+1)/(gamma-1)) + 1)]  [degrees]
d(nu)/dM = sqrt(M²-1) / [M(1+(gamma-1)M²/2)] × 180/pi
```

Negative angles and angles at or above the infinity limit receive a
`physical_domain` error on `value`. An angle below that limit but above the
complete value at M=1000 receives `outside_inverse_range`, also on `value`.
The endpoint is compared before binary64 export: rounding it first can admit a
mathematically out-of-bracket root. The existing inverse bracket is preserved.
Exactly zero degrees returns the sonic state. An extremely small positive input
can also round to M=1; the API does not claim recovery of an unrepresentable
Mach offset from a rounded input.

The solver updates a monotone bracket and uses the analytic Decimal derivative
for Newton steps. Out-of-bracket or stalled steps use the arithmetic midpoint.
A native calculation may propose a smaller initial bracket in a limited
ordinary-input range. Each proposed endpoint replaces the corresponding
original bound only after the complete Decimal expression confirms its sign.
Forcing the proposal to either end of the full bracket must preserve the root.
A computed zero or bounds exporting to the same binary64 value terminate the
iteration. Exhaustion raises `prandtl_meyer_not_converged`, never the last
iterate. Callback convergence failures retain their cause and identify the
supplied `value`. A nonzero forward angle that exports as zero is refused as
`numeric_range` on `prandtl_meyer_deg`.

These finite-precision signs are not rigorous interval enclosures or a proof
over every binary64 input. The legacy native Newton helper remains only for its
direct exhaustion regression and explicit historical fault restoration; the
public PM inverse no longer routes through it.

## Independent and hostile verification design

The new forward reference integrates the differential relation after
M=1+t² at 120 digits, scaling before integration for very large gamma. The
inverse reference uses separate mpmath transcendental arithmetic at 450 digits
and 160 bisections on the exact supplied float. It uses no production callback,
Decimal context or Newton iteration. A generating Mach is only a way to select
an input; its rounded angle is independently inverted before comparison.

Selected valid outputs must agree within two binary64 ULP. Boundary neighbors
test refusal classification before export. The matrix includes the old M=1.01,
angle=1e-4 and gamma=1e6 dispatch boundaries, near-unit gamma, gamma up to 1e308,
subnormal angles, positive forward underflow, the real API route and JSON output.
Hostile contexts are exercised with a cold pi cache. Deliberately shortened
inverse/series budgets and failed callbacks must produce typed failures.

The mutation fixture explicitly migrates the old precise-PM bypass to restore
native cancellation under its existing differential-quadrature selector. Other
historical results are not silently renewed. The selected new faults restore
rounded callbacks, premature endpoint export, lost range checks, incorrect
units/fields, inherited contexts, partial series and unconverged outputs.

## References

- [NASA Glenn PM relation](https://www.grc.nasa.gov/www/k-12/airplane/pranmyer.html).
- [NIST DLMF arctangent series and identities](https://dlmf.nist.gov/4.24).
- [Python 3.12 Decimal contexts](https://docs.python.org/3.12/library/decimal.html).

## Remaining requirements

The test-only high-precision reference is not a runtime dependency. Numerical
tests, a wheel build and hostile arithmetic checks cannot demonstrate the
unselected motor/feed hardware in its relevant environment. All nine hardware
evidence gates in `TRL_ASSESSMENT.md` remain separate, unmet requirements.

## Measured verification

The 10 November frozen PM candidate passed 301 tests and 12,703 subtests
across 23 modules. All 23 selected injected faults failed their respective
tests, with 46 passing original/restored controls. The copied-source build
and the installed wheel passed nine independently calculated PM cases with
no test-only arithmetic library installed in the smoke environment. The
receipts are under `output/assurance/pm_20261110` and
`output/assurance/pm_mutations_20261110`. These counts bind the candidate before
the subsequent normal-shock changes.

The original caller-context mutation survived because both public entry
points had already installed their own Decimal context before calling the pi
helper. The revised test seeds a cold cache directly under a five-digit
caller context and independently checks pi and subsequent public outputs.
The original/restored controls pass, while the mutation corrupts the cache.
