# Native reserve mass ledger and cutoff review

This review covers the Lunaris adapter's time-only mass consumption, represented
reserve cutoff dates, active prefix clocks, initial mass consistency and accepted
RHS start state. It does not establish correctness of the complete native solver,
all adapter behavior, physical validity or hardware readiness.

## Arithmetic contract

Public initial mass follows the existing finite-real-to-binary64 normalization
contract before entering the exact ledger. Parameter drivers and tabulated
endpoints have already passed their public normalization. Internal rational
coordinates remain unrounded until a scalar output is required.

For a segment of duration h, normalized elapsed time x in [0,1], and nonnegative
endpoint flow rates q0 and q1, consumed mass is

`h * (q0*x*(1-x/2) + q1*x*x/2)`.

Both endpoint weights are nonnegative. Exact supplied-binary durations and
endpoint values are integrated, summed and scaled before export. Intermediate
binary64 squares, duration products or per-burn mass subtraction cannot erase a
recoverable final result. Unrepresentable nonzero scalar exports fail explicitly.
This is exact arithmetic on supplied numbers, not a statement about real engine
flow accuracy.

The reserve target is the exact ledger difference between normalized vehicle
mass, protected dry mass and the existing declared inventory margin. That margin
is a numerical model policy. It is not a sensor uncertainty, hardware reserve
requirement or universal acceptance threshold. Exact state/limit comparison is
identified by the required hashed bridge policy
`exact_binary_state_and_limit_plus_declared_margin_v1`.

## Clock and cutoff contract

The inverse first finds the flow interval containing the first crossing, then
searches ordered binary64 epoch dates. The selected date's consumed mass is no
greater than the adjusted target. Where the crossing lies between dates, the
remaining numerical mass shortfall is exposed with an outward-rounded bound.
The bound covers this date selection only; it excludes physical input error,
solver integration error and model discrepancy.

Canonical dates name their exact profile coordinates. Other native dates use
exact epoch-minus-offset cancellation; native mass prediction and RHS evaluation
use the same canonical prefix clock. The stored integer ratio preserves the
selected cutoff coordinate. A protected cutoff truncates unused later events,
so collisions or overflow exclusively in that unused tail cannot reject an
otherwise representable prefix. Actual active-prefix collisions and positive
intervals that the clock cannot represent fail explicitly.

The binary lattice search takes at most 65 predicate calls: one upper-end probe
and at most 64 bisection probes. This bounds the number of probes, not wall time
or the cost of exact segment integration. A first crossing at a zero-flow
plateau's start is not moved to the plateau's end. The inverse makes a guarded
date choice and does not claim a universally exact analytical root.

Known initial mass is evaluated by summing all consumed mass before one export.
This retains small accumulated consumption and evaluates overlapping unprotected
burns at their actual ignition dates. Protected overlaps, state-dependent flow
surfaces and unmodeled jettisons are refused where they would invalidate this
precomputed mass model. Negative flow scales are rejected; a custom driver whose
declared domain permits zero retains that closed-flow control.

## Reproduced defects and repairs

- The pre-repair quadratic root returned zero for a positive approximately
  `5e-309 s` crossing and approximately twice the constant-flow reference time
  in a small-flow/large-duration case. Exact endpoint integration and guarded
  date inversion remove those range failures.
- Rejected initial RHS calls previously initialized the segment clock. Retrying
  an inconsistent mass could bypass its check, and the first accepted call at
  a later thrust jump could select the left side. Initial consistency now runs
  before the base RHS; the start date is saved only after a successful result.
  A base RHS failure also leaves the segment uninitialized.
- An accepted Fraction initial mass could be used exactly for cutoff construction
  but rounded in dynamics. A synthetic `2**53 + 3/4 kg` input produced predicted
  final mass one kilogram below its protected dry mass. Both paths now normalize
  the public input before exact arithmetic. This is a numerical counterexample,
  not a spacecraft mass recommendation.
- Event dictionaries were iterated as dates, silently dropping their values.
  Window records now reject mappings and snapshot accepted date/ratio inputs.

## Verification scope

Thirty methods in `tests/test_native_reserve_numerics.py` use independent exact
endpoint integration, paired Decimal precision (1200/1600), adjacent-date sign
inequalities and explicit malformed state/record cases. They include signed
epochs, subnormal outputs, rising and
falling flow, sequential and overlapping mass histories, a usable finite cutoff
with an overflowing nominal end, and an actually unrepresentable crossing.

The final methods are replayed unchanged against the verified pre-repair source
archive and the published `cf31f56` checkpoint. Historical failures caused by
missing new record/helper APIs are identified separately from assertion failures;
neither is counted as a distinct physical or software defect. Separate public
constant-flow reproductions establish the actual old arithmetic errors.

Fifteen additional curated faults cover unsafe clock neighbors, lost curvature,
premature scalar export, inward deficit rounding, rounded mass consistency,
premature RHS initialization, mutable driver aliases, public normalization,
record shape, scale sign and finite-clock rejection. Only assertion failures
count as killed faults. The finite set is not an exhaustive mutation analysis.

The existing RK integration regression compares against a closed-form rocket
velocity and retains its original mass and velocity tolerances.
With the corrected clamp arithmetic it can now succeed; the test accepts an
explicit domain rejection or requires those same physical invariants on success.
Collection and per-module execution budgets remain unchanged.

## Remaining evidence

Repository-wide function, catalogue and claim review remains partial. General
raw-document structure, unreviewed solver behavior and performance under broad
concurrent load remain separate work. No selected motor/feed design, approved
mission requirements, calibrated relevant-environment demonstration or independent
hardware validation exists. All nine hardware evidence gates remain unmet.
Software passes do not establish TRL 5 or certification.

## References

- [TABULATED_TIME_CONTRACT_REVIEW.md](TABULATED_TIME_CONTRACT_REVIEW.md): supplied
  binary clock, event ownership and public/native clock distinctions.
- [NATIVE_RESERVE_NUMERICS_CHECKPOINT.md](NATIVE_RESERVE_NUMERICS_CHECKPOINT.md):
  original publication checkpoint and its explicitly unfinished review scope.
- `src/rocket_propulsion/integrations/native_mass.py` and
  `src/rocket_propulsion/integrations/lunaris.py`: implementation and consumer.
- `tests/test_native_reserve_numerics.py`, `scripts/run_feed_mutations.py`:
  independent numerical oracles and configured semantic fault identities.
- [NASA-STD-7009](https://standards.nasa.gov/standard/NASA/NASA-STD-7009): external
  model credibility reference; this review does not assert full compliance.
