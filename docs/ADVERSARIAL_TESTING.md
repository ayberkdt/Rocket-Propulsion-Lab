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
| `test_api_integer_budget.py` | Compact exponent expansion, allocated oversized integers, host-limit configuration and HTTP serialization | An independent Decimal conversion tripwire stops unsafe allocation; exact boundary controls and real HTTP requests verify refusal and preserved large seeds. |
| `test_blowdown_record_ownership.py`, `test_blowdown_evaluation_reuse.py` | Mutable numeric aliases, caller-owned sample lists and repeated pressure queries | Before/after caller mutation, owned scalar/type checks, single-pass rows and independent analytic pressure/time/impulse histories; custom-provider controls remain uncached. |
| `test_blowdown_direct_queries.py` | Direct-query malformed reals, false-zero/nonfinite volume export, volume overflow masking finite pressure, small growth erased before a large exponent | Independent 100-digit polytropic and log1p-series witnesses, forward/inverse closure and true final-range refusal. |
| `test_api_regen_admission_order.py`, `test_regen_input_records.py` | Invalid controls reaching chemistry, mutable dimensions, truthy nonboolean flow direction, false rejection at EOS/pitch bounds | Chemistry sentinels establish ordering only; exact EOS endpoints and positive land remain admitted; real thermal/geometry/conjugate suites retain physics and energy checks. |
| `test_free_space_numerics.py`, `test_tabulated_ledger_numerics.py` | Intermediate range loss, false-zero diagnostics, directed tails, output rounding/range ties, scale/group order, HOLD and producer authority | Independent rational and paired-precision closed-form integrals/moments; direct enclosure inequalities; selected fault injections |
| `test_adversarial_epochs.py` | Left/right limits at every boundary of random jump profiles at epochs 1–3·10⁸ s; snapping; resolution limits; RK4 integration across biased burns at random epochs | Segment values set by the test; written-out linear interpolation; Δv = c ln(m0/mf) with exact trapezoid propellant |
| `test_adversarial_lunaris.py` | Coincident and back-to-back burns, bias that reorders burns, partial propagation windows, 400 random dry-mass cases (ramps up/down, plateaus, cuts on boundaries and at the stop), mass-check tolerance, predicted mass | Exact piecewise-quadratic mass and bisection; a two-burn lunar case integrated by SciPy DOP853 with a right-hand side written in the test |
| `test_adversarial_pulse_mode.py` | 120 random thrusters with thousands of requests, including interval end points ± 1 ulp; budget-limited trains; pulse limits; forged plans; tampered and hostile documents | Exhaustive enumeration of achievable Δv/impulse intervals per pulse count |
| `test_adversarial_json_readers.py` | All 14 persisted-document readers: more than 1 500 generated mutants plus raw hostile texts | Rule: refuse with `DomainError`/`InputError`, or return a model identical to the original |
| `test_constants_and_references.py` | Bare copies of physical constants, the 1976 standard atmosphere, docs without sources | Registered constants; independent RK4 integration of the hydrostatic equation |
| `test_api_adversarial_fuzz.py`, `test_api_numeric_admission.py` | Declared top-level scalar readers, shared API helpers and direct getters, one field at a time; nested schemas are not comprehensively covered | Malformed declared numbers must be refused before model execution with the correct field; other fields and semantic extremes receive typed-error/strict-JSON smoke checks. Fault controls reject injected finite success and unrelated-field errors. See `API_NUMERIC_ADMISSION_REVIEW.md`. |
| `test_oracle_compressible.py` | Normal and oblique shocks, Prandtl–Meyer, Fanno and Rayleigh, every inverse of the isentropic relations, Taylor–Maccoll cones | Rankine–Hugoniot residuals and entropy rise; quadrature of the differential relations; the area-ratio mass-flux identity; SciPy `solve_ivp` on the Taylor–Maccoll equation |
| `test_oracle_conduction.py` | Bessel functions, erfc solutions, slab/cylinder/sphere eigenvalues and series, fins | SciPy Bessel, mpmath erfc, `brentq` eigenvalues to 1e−15, a BDF method-of-lines solution with a sparse Jacobian, `solve_bvp` fins |
| `test_oracle_real_fluids.py` | Every Helmholtz-energy derivative property, state round trips, saturation | Richardson-extrapolated finite differences of the EOS itself; Mayer, sound-speed, Gibbs, Joule–Thomson and Z identities; equal Gibbs energy and Clausius–Clapeyron on the saturation line |
| `test_oracle_equilibrium.py` | Gibbs-minimum composition | KKT element-potential certificate; SciPy SLSQP never finds a lower G; extensive scaling; Le Chatelier; HP/TP consistency |
| `test_oracle_kinetics.py` | Finite-rate H₂–O₂ reactor and Jacobian | SciPy Radau; atom conservation; Kc from species data; fourth-order finite-difference Jacobian; detailed balance at the Gibbs equilibrium (slow, `RPL_SLOW_TESTS=1`) |
| `test_oracle_aerodynamics.py` | Lifting line and the vortex-panel method | Elliptic-wing closed forms; `e ≤ 1`; Kutta–Joukowski against integrated pressure; d'Alembert; Kutta condition; refinement; thin-airfoil limit; an independent Hess–Smith code |

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
- **Vortex-panel lift wrong on cambered sections.** The closed-form source and
  vortex influence integrals used the unsigned control-point distance, which
  is valid only on convex contours. On the concave lower surface of a NACA
  4412 the zero-incidence lift came out as 1.44; an independent Hess–Smith
  code gives 0.5172. The integrals are now evaluated in each panel's own
  frame with signed coordinates (`atan2` and `log` of end-point distances).
  Convex sections are unchanged to 4e−15.
- **Isentropic inverses slow.** Bisection inverted the area ratio to 2e−11
  (worst of 300 random cases) in 2.6 ms per call. Stagnation and critical ratios are now inverted
  in closed form with `expm1`/`log`; area ratio and Prandtl–Meyer angle use a
  bracketed Newton iteration with analytic slopes. The round trip is within
  a few ulp of M plus the input ratio's own conditioning (8.3e−13 worst for
  M ≥ 0.05 over 2 000 random cases), and an area-ratio inversion takes
  0.017 ms.
- **Equilibrium hang at absurd pressure.** A 1e300 Pa chamber pressure ran
  the ideal-gas Gibbs solver without end. Pressures above 1e9 Pa are now
  refused (`pressure_outside_ideal_gas_range`).
- **41 contract breaks under extreme values.** The first extreme-value pass
  over 2 856 route/field/value probes found 33 bare exceptions (21 `OverflowError`, 10
  `ZeroDivisionError`, 2 `ValueError` from the square root of a negative
  throat area) and 8 runs over 20 s. Fixes: `core.validation.numeric_range_guard` converts
  `OverflowError`/`ZeroDivisionError` into `DomainError(code="numeric_range")`
  on the 22 functions where they escaped; positivity checks on route lengths
  and areas; a one-million-sample budget on response grids
  (`sample_budget_exceeded`); loss-chain efficiencies checked before the
  chemistry runs.
- **Thrust rounded twice.** With no nozzle exit area the propagator subtracted
  a float zero from an exact `Fraction`, turning the thrust into a float before
  it was projected onto the force and acceleration axes. The zero is now a
  `Fraction`.

## Known limits (deliberate)

- A relative profile has no minimum duration. Absolute epoch placement requires
  distinct canonical dates for distinct boundaries; actual collisions are refused.
- The study workspace has no content hash, so an edited file gives an edited
  workspace. It is checked for stability, not for immutability.

## Reproducing selected numeric fault checks

`scripts/run_numeric_mutations.py` and
`tests/fixtures/numeric_admission_mutations.json` keep the selected faults and
their test selectors in version control. Run the preflight without executing
tests, or execute the faults in a new output directory:

```text
python scripts/run_numeric_mutations.py --list
python scripts/run_numeric_mutations.py --output-dir output/assurance/numeric_mutations
```

`--case NAME` selects a particular fault and can be repeated. The current case
file contains 63 selected admission, integer-budget, ownership, range, reuse and
test-oracle faults. That count is an inventory, not a passing result or an
exhaustive mutation score. Each case loads a fresh source/test copy, runs its
original control, applies one compilable change, restores the exact original
bytes and runs the restored control. A kill requires an assertion failure;
timeouts, collection errors, unrelated crashes, skips and surviving faults
remain open. The before/after controls must pass the same test count. Results
include source origins, repository fingerprints, JUnit logs, restoration hashes
and a verified baseline source/test archive. The checkout itself is never
mutated. Only verified temporary source/test copies inside the output tree
are removed after a successful case.

A previous local run exposed two verification defects. Text-based restoration
changed LF bytes to CRLF on Windows, so restoration now uses original bytes.
After the tank's direct pressure implementation changed, a fault in the shared
pressure-scaling helper was still paired with a tank test that no longer called
it. The helper fault now targets the performance-law oracle, and a separate
fault targets the tank's new pressure fallback. The incomplete earlier run is
not evidence that all selected faults were killed. Any surviving fault must be
examined as a test gap, an invalid selector or a genuine equivalent change;
changing its label does not turn it into a pass.

The first rebinding trial reached the correct helper but its zero-output fault
triggered an unhandled `DomainError` before the numerical assertion. The strict
runner kept that result open instead of counting a crash as a kill. The two
valid pressure-scaling rescue cases now establish a finite positive independent
expected value before evaluation and explicitly fail on a false domain refusal.
Other exception types remain unexpected failures, not accepted assertion kills.

These checks establish sensitivity to particular equation and refusal defects.
They do not provide calibration data, relevant-environment hardware evidence or
TRL certification.

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
