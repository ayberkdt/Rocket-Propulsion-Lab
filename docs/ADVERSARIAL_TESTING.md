# Adversarial testing

The `tests/test_adversarial_*.py` suites and `tests/test_constants_and_references.py`
are written to break the astrodynamics-facing code. They do not just confirm
that it works. Four rules apply:

1. **Independent oracles.** Every expected value is computed by code in the
   test file, not by the library: closed forms, exhaustive enumeration,
   independent integration with SciPy, or plain bisection. Each test module's
   docstring names the source of every equation it uses.
2. **Hostile inputs on purpose.** Large absolute epochs, boundaries hit from
   both sides, requests one ulp outside an achievable range, zero-flow
   plateaus, values of the wrong type, integers beyond the double range,
   NaN/Infinity tokens, and documents nested 50 000 levels deep.
3. **Seeded randomness.** Generators are seeded, so every failure reproduces.
4. **Mutation checks.** Each fix is checked by re-introducing the old
   behaviour (monkeypatching it back) and confirming that the suite fails. A
   test that cannot fail on the bug it targets was strengthened until it did.

## Suites

| Suite | What it attacks | Oracle |
| --- | --- | --- |
| `test_free_space_numerics.py`, `test_tabulated_ledger_numerics.py` | Intermediate range loss, false-zero diagnostics, directed tails, output rounding/range ties, scale/group order, HOLD and producer authority | Independent rational and paired-precision closed-form integrals/moments; direct enclosure inequalities; selected fault injections |
| `test_adversarial_epochs.py` | Left/right limits at every boundary of random jump profiles at epochs 1–3·10⁸ s; snapping; resolution limits; RK4 integration across biased burns at random epochs | Segment values set by the test; written-out linear interpolation; Δv = c ln(m0/mf) with exact trapezoid propellant |
| `test_adversarial_lunaris.py` | Coincident and back-to-back burns, bias that reorders burns, partial propagation windows, 400 random dry-mass cases (ramps up/down, plateaus, cuts on boundaries and at the stop), mass-check tolerance, predicted mass | Exact piecewise-quadratic mass and bisection; a two-burn lunar case integrated by SciPy DOP853 with a right-hand side written in the test |
| `test_adversarial_pulse_mode.py` | 120 random thrusters with thousands of requests, including interval end points ± 1 ulp; budget-limited trains; pulse limits; forged plans; tampered and hostile documents | Exhaustive enumeration of achievable Δv/impulse intervals per pulse count |
| `test_adversarial_json_readers.py` | All 14 persisted-document readers: more than 1 500 generated mutants plus raw hostile texts | Rule: refuse with `DomainError`/`InputError`, or return a model identical to the original |
| `test_constants_and_references.py` | Bare copies of physical constants, the 1976 standard atmosphere, docs without sources | Registered constants; independent RK4 integration of the hydrostatic equation |

## Defects these tests found

All of the following were fixed, and each has a regression test.

- **Clock ownership and event identity.** Earlier epoch-window repairs were
  superseded by exact relative ownership and once-rounded canonical epoch dates.
  Tiny interiors and distinct native markers are retained. Actual canonical-date
  collisions fail explicitly; ambiguous addition-order aliases do not override
  canonical events. See [the clock review](TABULATED_TIME_CONTRACT_REVIEW.md).
- **Dry-mass cut pushed past a zero-flow plateau.** A rounding difference in
  the sum moved the cut from 26.4 s to 48.5 s. The cut now uses the bridge
  clamp's margin (`INVENTORY_LIMIT_TOLERANCE`).
- **Near-stop cut handling.** Profile boundaries now use exact coordinates;
  distinct dates are retained. Reserve-root arithmetic remains a separate review.
- **Catastrophic cancellation in Δv for small burns.** `log(m0/(m0-u))` loses
  about 7 digits when u/m0 ≈ 10⁻⁷; pulse plans missed their target by
  1.3·10⁻⁹ relative. All five sites now use `-log1p(-u/m0)`.
- **Unphysical pulse plans.** Pulse plans could use the whole vehicle mass as
  propellant; `plan_delta_v` now requires a propellant budget.
- **Reader crashes.** In the first run, 58 sampled mutants or inputs broke the
  "refuse or unchanged" rule in 12 of the 14 readers, among them `OverflowError` from 10⁴⁰⁰,
  `ValueError`/`RecursionError` from NaN tokens and deep nesting, and
  `TypeError` on non-text input; bytes were silently accepted. All readers now
  go through `core.validation.read_json_document`.
- **U.S. Standard Atmosphere 1976 used the 2019 SI gas constant.** This put a
  pressure jump of about 3·10⁻⁵ at every layer base. It now uses the model's
  own R\* = 8.31432 J/(mol K).

## Known limits (deliberate)

- A relative profile has no minimum duration. Absolute epoch placement requires
  distinct canonical dates for distinct boundaries; actual collisions are refused.
- The study workspace has no content hash, so an edited file gives an edited
  workspace. It is checked for stability, not for immutability.

## References

- D. Goldberg, *What Every Computer Scientist Should Know About
  Floating-Point Arithmetic*, ACM Computing Surveys 23(1), 1991 — rounding of
  sums, the cancellation-free quadratic root and `log1p` (bibliographic
  citation).
- NIST DLMF §4.6.1, the series of ln(1+z) used to reason about `log1p`:
  <https://dlmf.nist.gov/4.6#E1>
- RFC 8259, *The JavaScript Object Notation (JSON) Data Interchange Format*:
  <https://www.rfc-editor.org/rfc/rfc8259>
- NASA Glenn, *Ideal Rocket Equation*:
  <https://www1.grc.nasa.gov/beginners-guide-to-aeronautics/ideal-rocket-equation/>
- U.S. Standard Atmosphere, 1976: <https://ntrs.nasa.gov/citations/19770009539>
- Orekit 13.1.5 `EventDetector` API (one-sided limits at events):
  <https://www.orekit.org/static/apidocs/org/orekit/propagation/events/EventDetector.html>
