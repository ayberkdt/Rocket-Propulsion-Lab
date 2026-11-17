# Bounded first stage for Prandtl-Meyer and area calculations

Reviewed 2026-11-17 on top of commit `625d145`. The complete Decimal stage
introduced by the isentropic, area/Mach-angle and Prandtl-Meyer precision
reviews ([ISENTROPIC_FORWARD_RANGE_REVIEW.md](ISENTROPIC_FORWARD_RANGE_REVIEW.md),
[AREA_ANGLE_INVERSE_REVIEW.md](AREA_ANGLE_INVERSE_REVIEW.md),
[PRANDTL_MEYER_REVIEW.md](PRANDTL_MEYER_REVIEW.md)) evaluated every forward
Prandtl-Meyer angle and every area or Prandtl-Meyer inverse at 90 or more
digits. That made ordinary design work slow: a 20-characteristic MOC nozzle
took 5.42 s instead of 0.04 s at `ef00f3e`, and Rao and bell contours took
0.67 s and 0.72 s instead of 0.02 s and 0.01 s.

## What changed

`compressible/sonic_numerics.py` now tries a first stage before the complete
stage. The complete stage, its admissions, iterations and refusals are kept
unchanged and run whenever the first stage does not prove its answer.

- Inverses: a native estimate proposes a binary64 Mach number. Decimal
  residuals at a lower precision (the complete precision minus 60 digits) are
  evaluated at the exact decimal midpoints between that value and its two
  binary64 neighbours. Each residual carries an allowance of
  `10**(8 - precision)` times the magnitude of its terms. When both signs are
  certain and bracket the root, the value is the binary64 number nearest the
  root, which is what the complete stage returns. A proved wrong side moves
  the candidate by a Decimal Newton step and is checked again.
- Forward angles: the first-stage value is accepted only when the whole
  allowance interval rounds to one nonzero binary64 value.
- Admission: values near the physical limit, the inverse bracket ends or the
  extreme gamma ranges skip the first stage, so every refusal still comes from
  the complete stage.
- Any uncertain sign, exact tie, nonzero underflow, failed estimate or
  arithmetic exception returns `None`, and the complete stage runs.
- MOC mesh construction memoizes repeated Prandtl-Meyer inversions within one
  build (1213 inversions had only 209 distinct inputs).

## Evidence

Equality, not tolerance, is the acceptance rule.

- 2,314 cases per seeded run (624 forward angles, 624 Prandtl-Meyer inverses
  and 1,066 area inverses) over thirteen gamma values from 1.05 to 1e6, both
  area branches and near-sonic inputs, run with two seeds: no difference
  between the first and complete stages, refusals included. On the final
  code the first stage answered 590 and 598 of the 624 Prandtl-Meyer inverses
  (two seeds) and 916 and 919 of the 1,066 area inverses.
- 36 MOC, Rao and bell designs (four area ratios, three gamma values) and a
  1,596-point isentropic sweep: identical SHA-256 of the serialized outputs
  before and after.
- `tests/test_sonic_fast_path.py` keeps a bit-for-bit comparison, checks the
  allowance and midpoint logic, and asserts that ordinary nozzle states are
  answered by the first stage. `tests/test_sonic_complete_stage.py` reruns the
  three precision suites with the first stage disabled. Two budget-exhaustion
  tests now select the complete stage explicitly, because the budget belongs
  to it.
- All 35 existing numeric mutation cases on this code are still killed, and
  six new cases for the first stage (ignored allowance, unproved forward
  interval, wrong midpoint, reversed boundary sides, lost subsonic orientation)
  are killed.

Timing on the review machine, same inputs as above:

| Operation | `ef00f3e` | `625d145` | First stage |
|---|---:|---:|---:|
| MOC nozzle, 20 characteristics, area ratio 25 | 0.04 s | 5.42 s | 0.16 s |
| Rao contour, area ratio 25 | 0.02 s | 0.67 s | 0.09 s |
| Bell contour, area ratio 25 | 0.01 s | 0.72 s | 0.10 s |
| Area ratio inverse, supersonic | 20 us | 3786 us | 642 us |
| Prandtl-Meyer inverse | 16 us | 3442 us | 487 us |

The first stage is still slower than the original native code because each
proof needs Decimal logarithms or arctangents. It is a rounding proof for the
retained perfect-gas expressions, not a physical accuracy bound.

## Normal shocks

`compressible/normal_shock.py` runs its unchanged Decimal expressions first at
28 digits, for 1.05 <= gamma <= 100 and Mach up to 1e6. Each output is
accepted only when its relative interval of `10**(8 - 28)` rounds to one
positive binary64 value; otherwise everything reruns at 80 digits, where all
refusals are made. Within that range the outputs carry a few correctly rounded
operations, and the recovery adds two logarithms scaled by gamma/(gamma-1);
their estimated error stays about 5,000 times below the allowance. A call at
Mach 2.7 and gamma 1.22 drops from 609 us to 173 us.
`tests/test_normal_shock_first_stage.py` compares 612 gamma/Mach pairs bit for
bit, refusals included. The six existing normal-shock mutation cases and three
new first-stage cases (unproved interval, ignored allowance, unreset stage
marker) are killed; the unproved-interval case first survived until a 1e-3
allowance test was added.

## Area inverses: table logarithm and anchored boundaries

The area first stage no longer calls Decimal `ln` at each boundary. Its
logarithm reduces a value to `2**k * u` with `1 <= u < 2` and uses a cached
256-point table of `ln(1 + j/256)` (five guard digits) plus a short `log1p`
series. The first boundary point is evaluated in full; nearby points within a
relative distance of 1e-6 extrapolate from it with `log1p` series. These are
first-stage arithmetic only; the complete stage still uses `_log_area`, and
every result is still accepted only between proved half-ULP boundaries.

- Warm area inverse on 300 mixed cases: 529 us at `eee42c1`, 161 us now, none
  deferred to the complete stage.
- Off-design nozzle call (area ratio 8, eight ambient pressures): 6.3 ms at
  `ef00f3e`, 73.7 ms at `eee42c1`, 30.3 ms now.
- Bit-for-bit comparison with two seeds: no difference; the 36 designs and
  1,596-point sweep keep hash `9bbe8be8de885cc6`.
- The table logarithm is checked against 60-digit Decimal `ln`, including
  powers of two, table points and values beyond binary64. All 41 earlier sonic
  mutation cases stay killed; four new cases (dropped power of two, wrong table
  base, lost series sign, reversed anchor offset) are killed.

## Open

Stagnation-ratio inverses still use only their complete stage. Forward
Prandtl-Meyer angles in the first stage still use the Decimal arctangent, so an
off-design nozzle call remains about five times slower than at `ef00f3e`.
The allowances are conservative engineering bounds derived from correctly
rounded Decimal operations, not machine-checked proofs.

## References

- NACA-TR-1135, perfect-gas isentropic and Prandtl-Meyer relations.
  <https://ntrs.nasa.gov/citations/19930091059>
- PYTHON-3.12-DECIMAL, correctly rounded context arithmetic.
  <https://docs.python.org/3.12/library/decimal.html>
