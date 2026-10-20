# Native reserve arithmetic checkpoint

The Lunaris adapter now integrates supplied binary flow endpoints with an exact
rational mass ledger. This avoids intermediate overflow or underflow in quadratic
root formulas and preserves small consumption across sequential and overlapping
unprotected burns. Scaling occurs before the consumed mass is exported as a
binary64 scalar.

Protected, nonoverlapping burns select a guarded cutoff on the actual represented
epoch clock. The selected date does not exceed the margin-adjusted consumption
budget. Its outward-rounded mass shortfall is reported separately. A positive
cutoff that the clock cannot represent raises a domain error. Only the active
profile prefix must have distinct event dates; an unused tail cannot invalidate
a usable cutoff. The existing absolute inventory margin remains a numerical
model policy and is not a physical measurement uncertainty.

Window records copy their event dates and exact cutoff ratios. The dynamics
wrapper snapshots parameter overrides and uses the same prefix clock for thrust
evaluation and mass prediction. Invalid time or vehicle mass fails before the
wrapped right-hand side is called. Bridge identities include the exact inventory
clamp comparison policy.

Independent tests use exact rational mass integration, paired high-precision
Decimal calculations and adjacent represented-date inequalities. Cases cover
large and small flow ranges, rising and falling flow, absolute epoch resolution,
zero-flow plateaus, overlapping burns, sequential mass rounding, mutable inputs
and typed input rejection. The legacy clamped RK integration regression accepts
an explicit domain rejection or requires the original mass and velocity accuracy
when integration succeeds; it no longer requires a particular rounding failure.

This records the original `cf31f56` publication checkpoint. Follow-up arithmetic,
RHS retry, normalization and adversarial checks are documented in
[NATIVE_RESERVE_NUMERICS_REVIEW.md](NATIVE_RESERVE_NUMERICS_REVIEW.md).
Repository-wide review is partial. These numerical checks do not demonstrate
physical validation, hardware readiness or TRL 5. The hardware configuration,
mission requirements and nine hardware evidence gates remain unresolved.

## References

- Exact supplied-binary integration and epoch clock contracts:
  [TABULATED_TIME_CONTRACT_REVIEW.md](TABULATED_TIME_CONTRACT_REVIEW.md).
- Numerical implementation: `src/rocket_propulsion/integrations/native_mass.py`
  and `src/rocket_propulsion/integrations/lunaris.py`.
- Independent checks: `tests/test_native_reserve_numerics.py`.
- NASA model and simulation credibility scope:
  [NASA-STD-7009](https://standards.nasa.gov/standard/NASA/NASA-STD-7009).
