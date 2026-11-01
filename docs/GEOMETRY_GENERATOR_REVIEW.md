# Geometry generator input and sonic-limit review

Baseline: `f1df7bcf0105cf1d23d4503e478e5724432412db`. This review concerns
`generate_nozzle_contour`, `generate_rao_contour` and `generate_moc_nozzle`, their
generated sampled contour, and the unit-area-ratio comparison/API path. It does
not qualify a flow model, manufacturing geometry or motor/feed hardware.

## Reproduced problems and changes

The general generator accepted boolean throat area, leaked `TypeError` for a
fractional station count or a text angle, and leaked `AttributeError` for a
nontext family. The documented area ratio of one generated repeated x=0 samples.
Rao also lacked an early finite gamma check; MOC treated truthy strings and
integers as its axisymmetric-correction flag. The prior thermal-consumer guard
refused bad geometry, but did not repair its producers.

All physical scalars now use owned finite reals before geometry/flow evaluation;
booleans, text, erased nonzero inputs and nonfinite values are refused. Resolution
counts require actual bounded integers. MOC requires an actual boolean flag;
the general family name is owned plain text before case/alias processing. Existing
domain and count bounds are retained with stable offending-field labels. Known
native overflow/zero-division failures are reported as numerical-range errors.

## Preserve the physical limiting case

`A_e/A_t=1` remains supported by the general conical/preliminary generator. It
uses all requested samples on the cosine contraction, ending once at x=0 and
Mach one. Diverging length and angle are zero; the *divergence* factor is one,
not a statement of overall nozzle efficiency. Metadata identifies the actual
cosine contraction instead of claiming a divergent Hermite or cone was built.

Rao and the supersonic MOC construction still require area ratio greater than
one. An API comparison at one returns the two supported general-family limit
records and explicit Rao/MOC exclusion reasons. Its exit-angle screening metrics
use zero for the no-divergence case; ordinary comparisons retain their behavior.

## Output and scientific limits

Producer-owned contour output is checked for finite dimensions and station
values, positive mandatory geometry, nonnegative Mach and strictly increasing
x. Nonfinite totals, erased positive dimensions and unresolved axial samples
cannot be successful output. Caller-created `NozzleContour` constructors and
the distinct export/thermal snapshot contracts are not changed.

This is a postflight safeguard, not a guarantee that every mathematically
representable geometry can be evaluated by the current native-float formulas.
Near-sonic inversion/resolution, extreme geometric scales and continuous-wall
approximation need their own numerical studies. No rounding-cell tolerance or
gas/coolant property range was widened here. MOC remains the existing planar
net with radial mapping, not a new full axisymmetric solver; this change checks
its flag and generated contour, not every characteristic-field record.

## Verification requirements

`tests/test_geometry_generator_contract.py` distinguishes admission sentinels
from computed geometry. It covers malformed scalars for each generator,
count bounds/types, text/flag ownership, physical domains and mutable arithmetic.
Valid controls include count 401 and characteristic count 80, supported shape
endpoints, conical/Hermite/Rao/MOC scale similarity and sonic-limit sampling.

Independent 80-digit area-Mach and stagnation-ratio equations check the general
and Rao station fields; those one-dimensional equations are not imposed on a
mapped planar MOC wall. Independent tangent geometry checks nondefault cone
angles and a separate `hypot`/`fsum` calculation checks sampled wall area.
Parsed STL normals and actual API JSON exercise the unit-area-ratio deliverable.
Actual methane EOS/transport and coupled cooling in both directions independently
close coolant enthalpy rise for that limiting geometry.

A large converging-only witness has finite individual cross sections but cannot
have a finite total wall area: `S >= 2*pi*r_t*L_c = 3*A_t`, beyond binary64 for
`A_t=1e308`. It must yield a labelled range failure, not successful infinity.
The independent inequality avoids declaring an evaluation-order limitation to
be a physical impossibility for a balanced finite input.

Selected fault tests remove ownership/guards, reintroduce duplicated sonic
samples, falsify divergence metadata and alter geometry/flow equations. They
require passing originals, assertion failures on the changed code and passing
byte-restored controls. Execution receipts and exact input hashes are retained
under `output/assurance/geometry_generator_review/`. Counts belong to terminal
receipts, not this plan. Software checks do not supply a representative brassboard,
measurement uncertainty or relevant-environment demonstration for TRL 5.

## Measured results and open precision witnesses

The nine-module source-bound run passed **139 methods and 8849 subtests**, with
zero failures, errors or skips. All **16 selected faults were killed**, with
**32 passing original/restored controls** and byte-exact restoration. Production
inputs stayed unchanged during both runs. These are selected controls, not a
claim of exhaustive functional or mutation coverage.

A separate process loaded the prior commit from a Git archive and compared 24
non-unit-area-ratio geometries across three scales: all result values were
exactly equal to the current source on this host. That baseline also reproduced
the successful infinite wall-area output; the new version refuses it with
`numeric_range` and field `surface_area_m2`. This version comparison is separate
from the independent analytical checks above. A copied-source build and clean
local-venv wheel install verified actual generator admission, the sonic limit,
API comparison/exclusions, exports and a real-fluid cooling calculation.

Six finite-input witnesses remain open in
`open_native_precision_findings.json`: for each generator, an exit area ratio
one ULP above one and the minimum positive binary64 throat area. The current
arithmetic refuses them through resolution, inverse-flow or control-point
failures; they were **not** declared mathematically impossible or turned into
tests that demand refusal. Cancellation-safe span/radius arithmetic and
near-sonic inversion require a separate repair and independent numerical proof.

That subsequent repair and its scope are recorded in
[GEOMETRY_PRECISION_REVIEW.md](GEOMETRY_PRECISION_REVIEW.md). The six-witness
receipt above remains a record of this page's baseline, not a current refusal
policy or qualification of the full generator range.

## References

- NACA Report 1135, isentropic area-Mach and stagnation relations:
  https://ntrs.nasa.gov/citations/19930091059
- NASA Glenn, nozzle geometry/area-velocity context:
  https://www.grc.nasa.gov/www/k-12/airplane/nozzled.html
- Python 3.12, real and integral numeric protocols (physical type/domain
  exclusions are application policy):
  https://docs.python.org/3.12/library/numbers.html
