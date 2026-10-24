# Ideal rocket and steady thrust arithmetic review

This covers `propulsion.performance`, its positive log/exp helpers and two API
adapters. It verifies implementation consistency, not physical model or hardware
accuracy. The global review and all nine motor/feed gates remain incomplete.
No NASA consulting, endorsement or TRL certificate was obtained.

## Equations and reproduced defects

Prior public calls accepted boolean Isp, exported infinite fields, raised bare
inverse exp overflow, lost about 11.25% delta-v for adjacent large masses, and
rejected representable positive thrust after cancellation. Further checks exposed
rounded-intermediate reuse and unchecked result records. These are arithmetic
counterexamples, not measured operating points.

The ideal model is delta-v = ve*ln(m0/mf), with ve = Isp*g0 when Isp is supplied.
Steady thrust is F = mdot*V + A*(Pe-Pa), ve = F/mdot, Isp = F/(mdot*g0).
Signed pressure thrust is allowed while net thrust must be positive. Losses,
changing performance, separation and transient flow physics are not validated.

Exact rational ledgers retain binary inputs. Each required output rounds once
at its own boundary. Forward mass ratio and velocity exports are separate from
the exact ratio/coefficient used for delta-v. Thrust, velocity, Isp and impulse
use the unrounded momentum/pressure ledger, allowing cancellation before export.
Exact zero pressure thrust is valid; erased nonzero mandatory outputs fail typed.

## Directed log/exp export

Shared directed contexts enclose rational arithmetic. Python Decimal ln/exp use
correctly rounded nearest-even results; adjacent decimal values expand them
outward. Monotonicity encloses evaluation at rational input bounds. For positive
exponents exp(x)>1; negative exponents satisfy 0<exp(x)<1, permitting intersection
with those exact bounds.

Precision increases through 48, 96, 192, 384, 768, 1536 and 3072 digits. Output
requires both exact decimal bounds to round to the same nonzero finite binary64
value. The overflow midpoint is MAX+2**970 and the zero midpoint is 2**-1075.
An entirely impossible interval produces named DomainError; straddling or
otherwise unresolved cells produce ConvergenceError after the budget. Outputs
slightly beyond MAX may round to MAX, and nonzero subnormal cells remain valid.
This finite budget does not establish universal accepted-input convergence or
concurrent-service performance.

Inverse solves independently export exp(x), m0*exp(-x), m0*(1-exp(-x)) and
1-exp(-x), where x=delta-v/ve. Displayed masses may coincide while positive
consumption remains representable. If x>1024, e>2 proves mandatory exp(x) exceeds
the binary overflow midpoint. That early guard avoids huge exponentials without
excluding recoverable intermediate values.

## Records and compatibility

Strict scalar constructors use shared typed validation. Rocket results retain
solve_mode, performance_input and numerical_policy; thrust retains duration and
policy. Every derived field is replayed from retained source inputs. Inverse
records use their original delta-v/coefficient, not their rounded final mass.
This establishes internal consistency, not source/campaign authenticity.

Both API metadata versions advance to 2.0. Existing calculation arguments and
form aliases remain; result schemas add retained fields. Constructed impulse
records require duration, and inverse records require their mode/input tags.
Unsupported policies and inconsistent records fail closed. JSON round trips
contain finite quantities. Model outputs are neither clipped nor tolerance-repaired.

Tests use independent Decimal equations at two precisions, exact rational thrust
ledgers, scaling/rounding boundaries, malformed inputs, forged records and faults
between passing original-code controls. Baselines/later evidence remain locally
under ignored output/assurance. They do not replace calibrated integrated tests,
relevant-environment hardware evidence or independent readiness review.

The first development fault trial retained 17 kills and three survivors, all
between passing controls. Initial examples masked rounded-force reuse for Isp
and impulse; the uncertain-cell fixtures also missed a finite interior interval.
The production kernel and fault definitions were preserved while independent
rational/Decimal witnesses and that missing interval strengthened the tests.
Those incomplete development results are not reported as completed verification.

## References

- NASA Glenn, ideal rocket equation:
  <https://www1.grc.nasa.gov/beginners-guide-to-aeronautics/ideal-rocket-equation/>.
- NASA Glenn, rocket thrust:
  <https://www1.grc.nasa.gov/beginners-guide-to-aeronautics/rocket-thrust/>.
- NASA Glenn, specific impulse:
  <https://www1.grc.nasa.gov/beginners-guide-to-aeronautics/specific-impulse/>.
- Python 3.12, Decimal ln/exp and directed contexts:
  <https://docs.python.org/3.12/library/decimal.html>.
