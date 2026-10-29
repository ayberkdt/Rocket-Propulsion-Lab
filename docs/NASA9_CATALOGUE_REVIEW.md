# NASA9 species data and numerical evaluation

The ten bundled gas species retain their first two NASA9 intervals, 200--1000
and 1000--6000 K. This review covers the data, polynomial evaluator, species
construction, catalogue ownership and forward species property arithmetic.
It does not qualify mixture calculations, inverse accuracy, equilibrium,
finite-rate chemistry, real gases, motor/feed hardware or a physical domain.
The original full review and all nine hardware gates remain incomplete.

The subsequent [frozen-mixture forward review](MIXTURE_FORWARD_REVIEW.md)
extends software verification to mixture composition, entropy and retained
records. It does not add physical validation or certify the inverse solver.

## Primary source and corrected data

NASA/TP-2002-211556 equations (1)--(3), printed p.2/PDF p.10, define cp/R,
h/RT and s/R. Gas standard pressure is one bar. Appendix A, printed p.7/PDF
p.15, gives R=8.314510 J/(mol K). The prior implementation used the 2019 SI
constant with these fixed dimensionless fits, shifting dimensional properties
by approximately 5.7 parts per million. The fit now has its own registered
defining constant. The other models retain their selected constants.

The official NASA/CEA database was pinned to commit
`985fd33368b2334c4162a586c499e0ab41d12bd3`, SHA256
`fa7746572952d74e249e818a82a35c113829742fb421a308e167185528884363`.
Independent fixed-column extraction compared all 180 coefficients, both bounds
of twenty intervals, ten molecular masses, formulas and formation enthalpies.
The extraction is frozen in `assurance/nasa9_catalogue_reference.json`.

Two atomic-hydrogen coefficient discrepancies were found: low-temperature b1
was 25473.6599 instead of 25473.70801; high-temperature a2 was -0.181935441
instead of -0.1819354417. Both corrected values also appear in the original
2002 report, printed p.126/PDF p.134, visually inspected. They were not fitted
to current outputs. The low-temperature enthalpy reference correction is
about 0.40002 J/mol. The other 178 coefficients matched the pinned database.
No higher-temperature interval or additional species was implicitly enabled.

## Input ownership and fit arithmetic

Each interval owns exactly nine finite real coefficients. Boolean, textual,
nonfinite, unordered and wrong-length inputs fail at construction. Fixed-size
coefficient iteration reads at most ten items. Direct cp/h/s calls require a
finite positive numeric temperature. Merged onto current main (2026-10-29),
these raw polynomial primitives keep main's contract: they also evaluate
outside their fitted interval for reference-enthalpy checks, while
Species.properties refuses any state outside the species range. The shared
breakpoint belongs to the upper fit when accessed through Species; direct
interval calls permit either of their own endpoints.

Species own normalized positive molar mass, unique positive formula amounts
and newly constructed immutable interval records. Intervals must be ordered
and contiguous. This structural check does not prove thermodynamic adequacy
or global cp positivity of an arbitrary custom fit. Forward species evaluation
requires exact retained cp-R to be positive. The separate exported cp and R
may round equal even while a finite positive cv remains. Gamma and mass
properties are computed before their separate display rounding.

Rational polynomial products, powers and cancellation retain normalized
binary64 inputs exactly. Logarithmic terms use an isolated Decimal.ln context
with a conservative one-decimal-ULP enclosure. Exact rational combinations
of those bounds must round to the same binary64 value before export.
Precision is bounded to 48,96,192,384,768,1536,3072 digits; unresolved rounding
fails explicitly. Nonzero underflow and overflow of mandatory fields fail
with a typed domain error; exact zero is allowed. Caller Decimal precision,
exponent limits and traps cannot change these computations.

This resolves numerical rounding of the retained fit. It is not a fit
truncation, measurement, model discrepancy or physical uncertainty bound.
Temperature inversion still uses the shared scalar bisection algorithm;
only malformed-target normalization was added here. No global monotonicity
certificate, discontinuous-fit inverse policy or root enclosure is claimed.
Direct construction of SpeciesProperties is not source-bound replay evidence.

The cached catalogue stores immutable records. Each public call returns a new
dictionary, so replacing or deleting its entries cannot corrupt later lookups.
JSON duplicate members, nonfinite constants, malformed structures, unsupported
schema, coefficient order, standard pressure and defining R are refused.
This is a packaged-data cache; external edits require cache invalidation or a
new process. No live data update or automatic source authentication is claimed.

## Independent attacks and preserved evidence

The original production state and the first twenty attack methods are kept in
`output/assurance/history/nasa9_original_failures/manifest.json` (407 files).
All twenty methods have at least one failure in that state; subtests are not
distinct defects. Primary-data mismatch, wrong fit constant, shared catalogue
pollution, borrowed coefficient lists, extrapolation, malformed inputs,
intermediate overflow/underflow and hidden positive cv are independently
exercised. The original air-like-mixture assertion used the wrong constant;
only its expected R/M value was corrected, with its precision unchanged.

Forward equations use an independent 180-digit oracle over every interval;
110-digit independent differentiation checks dh/dT=cp and T ds/dT=cp.
Extreme finite temperatures, exact cancellation, one-ULP range escapes,
finite positive cv erased by display subtraction, malformed JSON and context
isolation are additional attacks. No source-table comparison is an engine
experiment or independent physical validation of the full software.

The first fault trial is preserved in
`history/nasa9_first_fault_trial/manifest.json` (406 files). It was incomplete:
one witness survived compensated floating summation, two experiments mixed
expected refusals with unexpected exceptions, and one whitespace anchor did
not inject. An additional non-dyadic cancellation witness strengthens the
original test. Two focused boundary methods supplement all original cases;
no production code or tolerances were changed to make these faults fail.
All 24 fault identities are retained; the gap anchor alone was narrowed to
its comparison expression. Passing original controls before and after each
fault are required. This is a curated fault assessment, not an exhaustive score.

A subsequent 1448-method/688-fault green run still missed an environment
dependency. Changing Decimal.DefaultContext's Inexact trap caused a valid
fit to raise raw decimal.Inexact because a new Context inherited that trap.
The green state is retained in
`history/nasa9_green_before_default_context_attack/manifest.json` (636 files).
The new valid-input regression assertion-fails on the same kernel; its exact
source and failure are retained in
`history/nasa9_default_context_failure/manifest.json` (409 files).
The repair specifies all context settings, traps and flags explicitly.
Polynomial equations, coefficients, log-enclosure rule and precision budget
are unchanged. The final suite has 23 new methods and 25 selected faults,
including a fault that restores this inheritance. All earlier methods and
fault identities remain present.

## Merge onto current main: exported records and solver internals

Merged on 2026-10-29. Exact evaluation applies to exported species records
(`Species.properties`, the interval cp/H/S methods and API output). Iterative
solvers (local equilibrium Gibbs loops and result mixing, transport, equilibrium
derivatives, condensed reactant enthalpy and the frozen rocket expansion) use a
separate `_binary64_properties` path that keeps the previous floating
expressions. It was compared with the pre-merge main over 14,400 random states
of all 36 gas species and over six full rocket-performance cases (three
propellant pairs, equilibrium and frozen): every field was bit-identical. The
same comparison found exported values differing from the old floating
evaluation by up to 140,687 binary64 ULPs in cancelling states, which the exact
path removes. A per-temperature cache shares each enclosed logarithm.

## References

- McBride, Zehe and Gordon, NASA/TP-2002-211556, September 2002, equations
  (1)--(3), Appendix A, Appendix C and gas data in Appendix D.
  <https://ntrs.nasa.gov/citations/20020085330>
- Official NASA/CEA thermodynamic database, pinned source commit above.
  <https://github.com/nasa/cea/blob/985fd33368b2334c4162a586c499e0ab41d12bd3/data/thermo.inp>

The primary PDF SHA256 is
`0ab4154b0fdeac29581f8d047a4cb6836c138c0b41da8f28990ef7d3e2756765`.
Sources and inspected equation/data pages are retained under
`output/assurance/sources/thermochemistry_review`. No NASA consulting,
endorsement, hardware qualification or certification was obtained.
