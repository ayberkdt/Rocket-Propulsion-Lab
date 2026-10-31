# Cooling contour admission and model reporting

Baseline: `6688d91560d6fb574e47ef595fa051ce5e02e15a`. The reviewed boundary is
the sampled geometry actually consumed by `solve_regen_chamber`, including its
axial-conduction path. This does not qualify a nozzle contour or cooling hardware.

## Reproduced defects

Before this change, empty/reversed station collections, repeated axial positions,
boolean coordinates, a negative throat radius and incompatible radius/area-ratio
pairs all reached fluid work and axial dispatch. A frozen `NozzleContour` could
also retain a caller-owned station list and mutable scalar subclasses. Checking
only its record class did not establish a valid, stable thermal geometry.

The converged axial result also retained the local march's statement that there
was no axial conduction, alongside a separate warning saying it was included.

## Contract and requirements

Both prescribed Bartz profiles and the coupled cooling march now use
`snapshot_thermal_contour`. The existing Bartz geometry validator is shared;
there is no second approximate-tolerance validator at the coupled boundary.
Its rounding cells describe binary64 representation, not manufacturing or
measurement uncertainty. Scalars are never clamped to make inconsistent inputs fit.

| Requirement | Behavior | Adversarial / physical control in `test_regen_contour_contract.py` |
|---|---|---|
| GEO-01 | Own finite native-float x, radius, A/At, Mach and throat radius | Malformed scalar sweep; dispatched types; mutation during a real EOS call |
| GEO-02 | At least two strictly increasing axial samples; ordered finite iterable | Empty/singleton, reversed/duplicate positions, wrong containers, one-pass iterator |
| GEO-03 | Positive radii; nonnegative Mach; coupled A/At at least one | Zero/negative physical inputs; one-ULP-subunity value refused before work |
| GEO-04 | Radius, throat radius and A/At have overlapping representation cells | Incompatible radii/ratios; existing independent prescribed Bartz geometry controls |
| GEO-05 | Admission precedes EOS and axial iteration; caller changes cannot alter a solve | Pre-work sentinels; mutated list/scalars during a real solve; real axial reuse |
| GEO-06 | Preserve legitimate prescribed profiles and existing generators | Coflow/counterflow, partial span without throat, Mach zero, Hermite/Rao/both MOC families |
| GEO-07 | Describe the model actually computed without implying hardware validation | Real axial result contains its coupling statement and retains hardware-validation limitation |

`solve_regen_chamber` accepts `NozzleContour` or the existing owned `BartzContour`.
The first is reduced to throat radius and station x/radius/A/At/Mach; the second
is copied and revalidated. The axial outer solve receives that owned profile;
inner marches revalidate owned inputs rather than revisiting a caller iterator.
Known malformed iteration protocol errors are typed; unrelated application callback
failures propagate. As elsewhere, variable-size Python iterables must be finite.

## Scientific and evidence limits

The Mach field is prescribed. Admission does not infer a gamma, impose an
isentropic relation, or verify the flow solver. A partial thermal span need not
contain x=0 or a sampled throat. `NozzleContour.area_m2`, reported lengths,
surface area and stagnation ratios are unused by this solver and are not claimed
to be checked. Its general-purpose constructors remain separate review work.

The coupling still uses its existing sampled explicit heat-load march and
constant-property gas assumptions. No numerical integration rule, coolant range,
wall conductivity, structural acceptance or hardware envelope was changed.
The model note now distinguishes local and axial marches; its warning that the
calculation has not been validated against hardware data is retained.

Real solutions use the existing LOX/CH4 equilibrium fixture, actual methane EOS,
transport and coupled wall equations. They compare the original contour with its
owned profile, attack ownership during EOS evaluation and independently recompute
coolant enthalpy rise. These are conservation and software controls, not a new
independent experimental validation dataset.

## Execution evidence

Execution receipts are retained under `output/assurance/regen_contour_review/`.
Selected deliberately broken implementations run in isolated copies with passing
originals and byte-restored controls. Crashes, errors, skips, timeouts and surviving
faults remain open results. Final counts are reported only from terminal receipts.
This scoped work does not close the motor/feed relevant-environment, test-article,
measurement, scaling or independent-review requirements for TRL 5.

The scoped 18-module run exercised **190 methods and 1284 subtests**: 189 methods
and all subtests passed; one documentation-index assertion failed because this
page's NASA overview URL was not yet registered. The URL was added to
`REFERENCES.md`; all 10 reference checks then passed. The red aggregate receipt
is retained rather than rewritten as a green command. Executed numerical source
and test-module bytes remained unchanged across this documentation repair.

All **12 selected faults were killed**, with **24 passing original/restored
controls**, byte-exact restoration and unchanged production inputs during the
final run. A first attempt at the axial-skip fault yielded an unexpected
`TypeError`; the harness correctly left it OPEN. Its selector was changed to the
existing reversed/duplicate-position test with valid containers, so the intended
missing axial guard is demonstrated by an assertion. No test assertion, solver
domain or kill-classification rule was weakened. The first attempt is retained.

The wheel and source distribution built from a copied snapshot; a fresh local
venv verified the installed shared contour API, a real-fluid solve, matching
owned-profile results and refusal before EOS/axial work. This is same-host
packaging evidence. It does not replace a supported-platform CI run or hardware
demonstration.

### Still-open generator findings

The generator itself is not repaired by protecting this consumer. A separate
probe reproduced a boolean throat area being accepted, a fractional station
count and text angle leaking `TypeError`, a nontext family leaking
`AttributeError`, and `area_ratio=1` producing repeated axial positions despite
being a documented accepted input. Receipts are in `open_generator_findings.json`.
Those construction contracts and the degenerate-span policy remain open work;
they must not be counted as qualified by these cooling tests.

## References

- NASA SP-125, printed pp.100–101, Bartz geometry and correlation context:
  https://ntrs.nasa.gov/api/citations/19710019929/downloads/19710019929.pdf
- Wadel, NASA/TM-1998-206313, sampled gas/wall/coolant marching method
  (bibliographic source; see `REGENERATIVE_COOLING.md`).
- Python 3.12, exact rational values and adjacent floating-point values used for
  representation-cell coherence, not physical uncertainty:
  https://docs.python.org/3.12/library/fractions.html and
  https://docs.python.org/3.12/library/math.html#math.nextafter
- NASA technology readiness levels, relevant-environment distinction:
  https://www.nasa.gov/directorates/somd/space-communications-navigation-program/technology-readiness-levels/
