# Sampled geometry exchange review

This review repairs preliminary text exchange of a prescribed sampled fluid
domain. It does not validate generic/Rao/MOC contour construction, ideal-flow
states, wall thickness, material, tolerances, manufacturing or motor/feed
hardware. No NASA consulting or endorsement was obtained. Source-bound results
are retained in `output/assurance/verification_summary.json`; all nine hardware
gates remain unmet.

## Contract and numerical policy

`ContourExportInputs` owns the plain contour name, throat area, nominal exit
area ratio, actual integer angular resolution and copies of all seven consumed
station fields. Scalars must be finite real values; booleans/text are rejected
at the Python boundary. Radius, area, area ratio and static/total ratios must
be positive, Mach nonnegative, and at least two axial stations must increase
strictly. The API retains its explicit numeric form parsing. Mapping reopening
snapshots nested data. Unused dimensions and claims on `NozzleContour` are not
authenticated by exporting its consumed station snapshot.

Area/radius/throat redundancies must overlap closed binary64 rounding-cell
enclosures under a **fixed binary64 math.pi basis**. These necessary
representation checks are not measurement uncertainty or flow validation.
CSV/DXF/STL numbers use 17 significant digits for binary64 round trips. DXF
declares metres (`$INSUNITS=6`), AC1015 and entity/polyline subclass markers.
STL side and cap winding is outward. Exact rational cross products of serialized
vertices determine normals; each component is rounded once from its exact
normalized squared value. Private facet areas/squared areas need not fit
binary64. Degenerate facets and unrepresentable mandatory nonzero normal
components fail with `DomainError`.

Cardinal ring axes are exact; other angles use platform sine/cosine. No
continuous-geometry tessellation error bound is established. STL has no
standardized unit metadata: preserve the bundle's `coordinate_unit='m'` and
choose metres in the receiver. `mesh_role='closed_fluid_domain_visualization'`;
`manufacturing_qualified` and `physical_model_qualified` must be actual False
booleans. This mesh has no wall thickness. SVG scales axes separately into a
viewport and is not a drawing to scale. Rational affine mapping avoids tiny
length floors and overflowing spans. XML metadata is escaped, literal carriage
returns preserved and forbidden characters rejected. The bundle normalizes
plain strings and replays all files from owned inputs; hostile equality cannot
authenticate forged text. The render cache retains at most one mesh.

## Independent checks and retained evidence

The 32 methods in `test_geometry_export_contract.py` inspect serialized text,
directed-edge pairing, nondegeneracy, normal direction and positive signed
volume. Cylinder and converging/diverging polygonal frustum volumes have
independent analytic oracles at several angular resolutions. Decimal 110/180
digit calculations check sampled normals. Tiny (r=L=1e-150) and large (r=1e150,
L=1e200) witnesses exercise private arithmetic range. Other checks cover exact
coordinate round trips, XML injection, extreme SVG spans, typed failures,
snapshot/mapping/tamper attacks and API modes. The existing DXF test now parses
all vertices instead of checking a rounded substring.

The committed failing baseline is retained in
`history/geometry_export_first_baseline/manifest.json` (389 verified files).
Its pytest summary `40 failed, 10 passed, 6 subtests passed` includes subtest
reports, not distinct defect counts. Some failures concern new ownership
requirements rather than previously implemented behavior. The first repair
exposed overflow in the test oracle itself: its cross product now converts
vertices to fractions **before** multiplication. That failure remains archived.

Twenty-six curated faults cover selected export/API decisions, not an exhaustive
mutation score. An accepted kill requires an assertion failure, no errors, and
passing original controls before and after injection. The first trial had 25
kills and one invalid mixed assertion/error result, preserved in
`history/geometry_export_first_fault_trial/manifest.json` (391 verified files).
A separate reversed-order witness isolates the intended assertion without
removing the original multi-case test. Production export kernels and the fault
operator did not change for that refinement. Browser checks execute the actual
scope renderer against an instrumented DOM; browser layout, commercial CAD
import and formal full-format conformance remain unassessed.

## Scientific geometry findings remain open

The CAD repair does not resolve six generic contour probes retained in
`geometry_followup.json`: boolean inputs, untyped float-count/None-selector
failures, forged dimensions and geometry overflow. Rao rounded-throat
construction, comparison proxies and scientific claims remain under review.

For At=0.01 m2, area ratio 25, gamma 1.22 and count 10, an independent diagnostic
checks 55 positive characteristic chords. Raw planar directions agree with
their midpoint characteristic angle to about 1.42e-14 degrees; returned mapped
coordinates deviate by up to 29.90374541382424 degrees. A tiny advertised wall
residual only checks its preserved boundary construction. Switching
`axisymmetric_correction` changes labels but returns identical geometry/Mach
and characteristic lines. At count 80 the refinement diagnostic compares the
same resolution to itself and reports zero change. These executable findings
are retained in `geometry_scientific_findings.json`; the solver/construction
review remains open. A multidimensional wall Mach need not obey the
quasi-one-dimensional global area-Mach relation; that relation is not used as
a universal MOC validity criterion.

The separate [interior characteristic unit-process review](MOC_UNIT_PROCESS_REVIEW.md)
adds a verified foundation with a genuine cylindrical source term. It has not
yet been integrated into a full nozzle solver and does not close these legacy
MOC findings.
The [axis/wall boundary review](MOC_BOUNDARY_REVIEW.md) supplies local
boundary evidence; the [axis departure review](MOC_AXIS_DEPARTURE_REVIEW.md)
retains a failed coupled stencil and its second-order regularity repair.
Neither qualifies the legacy mapped field.

## References

Autodesk documents vertex/subclass/closure codes and insertion units:
[LWPOLYLINE](https://help.autodesk.com/cloudhelp/2018/ENU/AutoCAD-DXF/files/GUID-748FC305-F3F2-4F74-825A-61F04D757A50.htm),
[INSUNITS](https://help.autodesk.com/cloudhelp/2025/ENU/AutoCAD-Core/files/GUID-A58A87BB-482B-4042-A00A-EEF55A2B4FD8.htm).
The Library of Congress describes oriented STL meshes and metadata limits:
[STL family](https://www.loc.gov/preservation/digital/formats/fdd/fdd000504.shtml),
[ASCII syntax](https://www.loc.gov/preservation/digital/formats/fdd/fdd000506.shtml).
Its account notes that the original specification was not recovered; it is not
treated as a formal current conformance standard. XML follows
[W3C XML 1.0](https://www.w3.org/TR/xml/).
