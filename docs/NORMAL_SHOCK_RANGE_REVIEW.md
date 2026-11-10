# Normal-shock arithmetic and output admission

## Reproduced failures

The preceding implementation exported stagnation-pressure recovery
1.0000000000000002 at gamma 1.4 and the first binary64 Mach above one. In the
calorically perfect normal-shock model the complete recovery is below one
and rounds to one at this input. At the first binary64 gamma above one and
Mach 2, separate powers overflowed even though the reference recovery is
about 0.613419867379714. At gamma 1e308 and Mach 2, intermediate overflow
exported infinity and NaN even though all five derived outputs are finite.
Eleven assertions failed in the initial independent run, including refusal
metadata and the ordinary API dispatch path.

## Numerical method and contract

The normal-shock relations use complete intermediate expressions with a
private 80-digit Decimal context. The static ratios are evaluated before any
output is rounded. The logarithm of total-pressure recovery is

```text
log(p02/p01) = [gamma log(rho2/rho1) - log(p2/p1)] / (gamma - 1)
```

The solver exponentiates the combined logarithm and accepts positive finite
binary64 results, including subnormals. A recovery whose logarithm is below
-800 rounds to zero in binary64 and is rejected without exponentiation.
Other unrepresentable outputs raise `DomainError` with `numeric_range` and the
affected output field. Direct calls validate real, finite inputs, reject
booleans and name the actual invalid field. Mach and gamma retain the
existing greater-than-one domain. Near the sonic limit, recovery may round
to one.

## Independent and adversarial checks

The oracle derives the compressive specific-volume ratio from mass,
momentum and total enthalpy conservation. It obtains temperature and Mach
from that ratio and uses entropy computed from static temperature and pressure
for total-pressure recovery. The reference uses 520-digit mpmath arithmetic,
independent of production Decimal and with exact binary64 inputs. Cases
include adjacent sonic and gamma endpoints, finite overflow-balanced outputs,
representable subnormals, a truly underflowed recovery, a truly overflowed
pressure, malformed direct inputs and strict JSON from the API. The
caller-context case traps inexact Decimal arithmetic outside the solver.

Selected injected faults reverse the entropy sign, compute entropy from
prematurely rounded static outputs, accept zero recovery, report the wrong
overflow field, inherit the caller context and accept Boolean gamma.

This work concerns arithmetic and input/output contracts for the calorically
perfect model. The extreme numerical inputs are not the operating envelope
of a real shock or evidence of motor/feed readiness.

## References

- [NACA Report 1135, Equations, Tables, and Charts for Compressible Flow](https://ntrs.nasa.gov/citations/19930091059)
- [Python 3.12 Decimal contexts](https://docs.python.org/3.12/library/decimal.html)
