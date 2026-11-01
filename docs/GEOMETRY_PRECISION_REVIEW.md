# Near-sonic and subnormal sampled geometry

Baseline: `1b1028b6da3ba87f72718b17f5a23a609ff6fca6`. The six recorded finite
counterexamples are inputs, not accepted numerical limits: one-ULP-above-unity
exit area ratio and minimum positive binary64 throat area, for all three
general/Rao/planar-mapped-MOC generators. Hardware readiness is not assessed.

## Failures and numerical repair

`A_t/pi` could round to zero before a representable square root. Private
`pi*r*r` and segment areas could disappear before the complete area/sum.
`sqrt(epsilon)-1` and `r_e-r_t` erased real growth, falsely making Rao's tangent
intersection infeasible. Native area evaluation at M=1 could equal a target
one ULP above one; a broad station shortcut also flattened positive increments
to Mach one. PM subtraction and inverse termination lost small angles/offsets.
MOC discarded normalized wall height `1+growth`. Its fixed `1e-12 m` refinement
denominator falsely improved tiny-scale mesh indicators. The ordinary inverse
also silently returned a last iterate after exhausting 200 Newton iterations.

Throat radius now uses the existing checked rational square-root primitive on
exact `A_t/pi`, without exporting the private quotient. Pi is the binary64
`math.pi` basis used by geometry exchange. Growth is retained separately as
`r_t*(epsilon-1)/(sqrt(epsilon)+1)`; the denominator uses the native square root,
while its complete product/quotient remains rational until export. Lengths,
Rao control-point increments and translated Hermite/Bezier coordinates retain
growth before adding it to the throat radius.

Station area is complete `pi*r^2`; interior ratios divide that complete value
by the nominal input `A_t`, avoiding a second throat inferred from rounded rt.
The constructed throat/exit retain their defining ratios 1/epsilon, and exit
radius is rounded from complete `A_t*epsilon/pi`. These are design boundary
data, not clamps of caller geometry. A scale sweep exposed export-coherence
failures when rounded-radius ratios were substituted for those boundary data;
the exporter tolerance was retained. Each segment length is rounded from its
exact rational squared coordinate differences; frusta using those lengths
are summed exactly before final export. Individual segment areas need not be
representable. This is not globally correctly rounded continuous-wall quadrature.

Near-sonic area evaluation/inversion uses a logarithmic relation and private
Decimal context, with additional precision for large gamma. The new bracket
stops when both endpoints export to the same binary64 value; exhaustion refuses
a result. PM uses a factored small-root series and half-angle-reduced arctangent
for wider arguments. Degree conversion retains the binary64 pi basis. Public
flow inputs own finite reals/text; booleans and caller arithmetic are refused.

An additional environment probe found that `Decimal(float)` in precision
selection could trip a caller's `FloatOperation` trap before the private context
was entered. That conversion now uses exact `Decimal.from_float`. Every private
context setting, trap and flag is explicitly selected instead of inheriting
mutable `decimal.DefaultContext`. Tests poison the active context and process
defaults separately, check caller flags/settings remain unchanged, and restore
all global defaults. Earlier receipts are retained under
`pre_context_isolation/`; final source-bound controls are rerun after the repair.

The supported inverse Mach bracket is unchanged; it is an implementation limit,
not a physical validity claim. Ordinary inverse algorithms remain distinct;
their budget exhaustion now raises a typed error. These changes do not prove
every inverse, native callback plateau or floating-point input correct.

For small expansions MOC works in `q=y-1`: wall q=0, centreline q=-1.
Intersections recover wall height from the wall line, preserving small growth.
Dimensional mapping consumes that growth directly. Public coordinates remain
physical radii/positions. Refinement divides by the actual positive refined
length; repeating a capped characteristic count is not further refinement.

## Independent controls

`tests/test_geometry_precision.py` checks both area roots with independent
100-digit solves, forward area through mass flux, and PM through scaled
quadrature of its differential equation. Synthetic gamma extremes exercise
arithmetic, not physical-gas applicability. Hostile caller Decimal contexts,
legitimate rational inputs and malformed scalars are separate controls.

All six witnesses pass thermal geometry admission, actual mesh export and
strict JSON. Independent 120–160-digit calculations check throat radius,
frusta, retained spans and Rao control-point samples. A separately built
100-digit planar fan uses independent roots and normal-vector matrix
intersections. It shares the physical fan discretization; it is not an
independent physical model or measurement dataset. Its near-sonic length
tolerance explicitly scales with `ulp(M)/(M-1)` and count: a screening budget,
not a rigorous interval enclosure. Mesh/physical errors remain separate.

A bounded smooth exponential with known root demonstrates that 200 Newton
steps need not converge. Returning the last iterate is refused. Selected
faults reintroduce cancellation, premature rounding, offset loss, caller-context
dependence, false convergence and incorrect coordinate mapping. Originals,
assertion-only kills and byte-restored controls are required. Receipts are in
`output/assurance/geometry_precision_review/`; mutation anchor transitions are
explicit in `docs/assurance/geometry_precision_mutation_transitions.json`.
Historical receipts remain bound to their original sources.

## Measured validation receipt

The final source-bound run completed **200 test methods and 10067 subtests**
across 14 geometry, flow, export, cooling, HTTP and reference modules. Its JUnit
receipt has zero failures, errors and skips. The two mutation groups killed
20 distinct selected faults, with 40 passing original/byte-restored controls.
All mutant failures were assertions, not import errors, crashes or timeouts.
Their 455 archived source/test files still match the current tree. Four later
Sentinel research documents changed independently and are named explicitly in
the receipt, rather than silently treating the whole tree as byte-identical.

The first MOC offset-removal mutant survived a finite-output-only selector.
That OPEN receipt is retained. The existing independent high-precision matrix
intersection test detects its 4.61% length error; the selected test changed,
without weakening its tolerance or rewriting the physical reference. A separate
caller `FloatOperation` reproduction then led to the explicit private-context
repair. Final controls ran again after that repair. Earlier successes are not
substituted for those final-source controls.

A copied-source build and clean-venv wheel smoke completed on this Windows host,
including the new helpers, all six witnesses, actual cooling/export consumers,
nominal boundary ratios and hostile Decimal context. This is not a clean-machine,
cross-platform or hardware test. The source-bound final summary and byte hashes
are in `output/assurance/geometry_precision_review/summary.json`.

## Remaining scope

A separate `gamma=1e308, M=2` forward-state range failure is still observed.
Its complete-expression range is not repaired or qualified here. Ordinary
inverse rounding, all model domains, real-gas applicability, continuous-wall
accuracy, manufactured geometry and motor/feed relevant-environment evidence
remain separate work. These numerical controls do not establish hardware TRL.

## References

- NACA Report 1135: https://ntrs.nasa.gov/citations/19930091059
- NASA Glenn, isentropic equations and area branches:
  https://www.grc.nasa.gov/www/k-12/airplane/isentrop.html
- Aldo and Argrow, MOC context (not experimental validation of this code):
  https://ntrs.nasa.gov/citations/19940019183
- Python 3.12 rational arithmetic and isolated decimal contexts:
  https://docs.python.org/3.12/library/fractions.html and
  https://docs.python.org/3.12/library/decimal.html
