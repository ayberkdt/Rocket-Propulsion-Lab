# Interior characteristic unit-process review

This is a prerequisite for a genuine axisymmetric nozzle solver. The new
`solve_moc_interior_point` kernel solves one downstream interior intersection
of C+ and C- using constant-gamma, inviscid, irrotational flow. It is not a full
nozzle designer. The existing mapped `generate_moc_nozzle` remains unchanged;
its 29.9-degree mapped-net counterexample, ineffective axisymmetric flag and
maximum-resolution self-comparison remain open. The original project goal,
broader catalogue/claim review and nine motor/feed hardware gates are retained.

Separate [axis and prescribed-wall processes](MOC_BOUNDARY_REVIEW.md) now
provide software-checked local boundaries. Their assembly with the interior
kernel, a sonic starting field and complete nozzle validation remains open.

A [connected supersonic strip](MOC_FIELD_REVIEW.md) now assembles these
processes from an existing inlet to a finite straight wall, with complete
cell coverage and mass diagnostics. It still needs a precomputed supersonic
axis history and does not close sonic initialization or full-nozzle design.

## Equations and scope

With velocity components u,v and sound speed a, the characteristic slope is
lambda = tan(theta +/- mu). Each family obeys

```text
A du + B dv - S dx = 0
A = u*u - a*a
B = 2*u*v - A*lambda
S = sigma*a*a*v/r
sigma = 0 (planar), 1 (axisymmetric)
```

Velocities are nondimensionalized by sonic speed a*. The constant-gamma energy
relation gives `a*a = (gamma+1)/(2+(gamma-1)*M*M)` in these units. The two
geometric and two compatibility equations use arithmetic endpoint averages.
Velocity components and finite states are iterated together. Rational private
products, differences and line/linear solves avoid premature binary64
overflow/underflow; platform trigonometric results are binary inputs to these
ledgers. No rigorous enclosure of trigonometric or integration error is claimed.

The scanned Aldo/Argrow paper's sentence assigns sigma=0 to both geometries.
That text is inconsistent with cylindrical divergence. Sigma=1 is established
from cylindrical continuity, spherical radial-flow conservation and the
nonzero radius terms in TR R-33, not accepted on citation authority alone.
The Cartesian equations and endpoint-averaging principle are used here;
variable-composition results and published full-nozzle calculations are not
reproduced by this constant-gamma interior kernel.

Parents must be distinct, axially supersonic, and intersect downstream.
Axisymmetric radii must stay positive. The zero-radius regularity limit,
sonic starting line, corner expansion, wall streamline, characteristic
crossing/shocks, full field mass conservation and independent nozzle exit
profiles still need their own processes and validation. An x-parameterized
intersection cannot represent vertical characteristics. Unresolved iteration
or branch conditions fail with typed errors. No universal numerical-domain
or global convergence guarantee is established.

## Precision and retained records

The raw normalized discrete-equation residual remains visible. Its iteration
decision separately subtracts a conservative output-coordinate rounding-cell
allowance from geometric residuals. Compatibility residuals receive no such
allowance. The unchanged iteration threshold is 64*ulp(1); the default budget
is 128 iterations, with an explicit supported resource budget of 1..256.
`rounding_adjusted_residual_max`, `coordinate_rounding_allowance_max` and
`representation_limited` expose this distinction. All exported residuals are
rounded upward from their exact represented ledgers. Coordinate cells are not
measurement uncertainty or characteristic truncation bounds.
An exported point must still be strictly downstream of both parents;
coordinate-cell allowances cannot excuse coincident marching nodes. A
mathematically downstream intersection that rounds onto its parent x fails
with `moc_unit_unresolved_topology`.

`solution_kind` distinguishes a genuinely constant-state intersection from
an iterated one. A one-ULP parent Mach difference cannot be labelled exact
constant flow. The constant-state branch preserves high Mach directly rather
than reconstructing it from a rounded limiting velocity. Inputs/results own
their nested data and replay all fields. Qualification flags must be actual
False booleans; `discretization_error_bound` remains None.

## Independent tests and preserved failures

Twenty-five methods include closed-form constant flow, inclined planar flow,
exact radial/spherical flows, independent 110-digit Cartesian balances,
110/180-digit fixture stability, coordinate scaling by powers 2**(+/-1000),
M=1e200, and a translated planar intersection with overflowing private line
products. For exact radial flows, rho*V*R**d is constant (d=1 planar, d=2
axisymmetric). Their analytic Mach lines give polar-angle change
`+/- (nu_final-nu_parent)/d`. These oracles do not call the production MOC
equations, root solvers or mesh generator. Five spacings in four smooth
cases check third-order **local** error; they do not establish global nozzle
accuracy, an empirical domain or hardware maturity.

The first 20-method trial is archived at
`history/moc_unit_first_trial/manifest.json` (390 files). One fixed accuracy
bound was missed at the finest originally sampled planar step. A finer step
was added while retaining that bound and the order check. This exposed a real
coordinate representation limit in a different case: raw geometric residual
1.80638912423169e-14, with compatibility residuals below 9e-17. Its 128-iteration
failure is retained at `history/moc_unit_coordinate_rounding_failure/manifest.json`
(392 files). The coordinate-cell repair separates that effect without
relaxing the iteration budget, threshold or compatibility equations.

The subsequent 1320-method full pass and 568 selected fault kills still
missed a topology witness: at x=1e308 a valid private intersection rounded
onto its parents' x, while a huge coordinate allowance masked the geometric
residual. That green-but-inadequate state is retained in
`history/moc_unit_pre_topology_guard/manifest.json` (610 files). A separate
stored-topology guard and regression now reject it. The earlier green
results are not used as evidence that every numerical corner was covered.

Twenty-one curated faults cover source sign/magnitude, family directions,
coefficient averaging, private range, branch/budget decisions, coordinate
cells, visible residuals and retained records. Accepted kills require
assertions, no errors and passing unmodified controls before and after each
fault. This is selected semantic coverage, not an exhaustive mutation score.
After the topology guard was added, removal of the earlier upstream guard
survived the original generic exception assertion: the second guard still
refused the point, but reported output quantization instead of a mathematically
upstream intersection. That 20-kill/one-survivor trial is retained at
`history/moc_unit_post_topology_first_fault_trial/manifest.json` (391 files).
The test now checks the public machine error code as well, preserving the
distinction between invalid family ordering and unresolvable output topology.
The production kernel and injected operator were unchanged for this refinement.

A later full run passed 104 of 105 modules; the remaining module encountered
a disappearing root temporary directory during pytest collection. The failed
run is retained in `history/moc_unit_root_temp_collection_failure/manifest.json`
(604 files). Three existing workspace-scoped temporary fixtures now use separate
stable parents under `output/assurance/test_artifacts`, keeping transient peers
out of the project-root scan. Their assertions and project containment are
unchanged, as are collection/module deadlines, worker count, diagnostic plugin
and test selection.

TM X-1502's printed planar example was screened as a future full-nozzle
benchmark. Its finite-grid scheme differs from the legacy raw net; a
coordinate difference is not automatically a defect or a validation pass.
Source PDFs, relevant rendered pages and screening metadata are pinned under
`output/assurance/sources/geometry_review/moc_unit_source_manifest.json`.

## References

- Aldo and Argrow, 1993, Cartesian characteristic/compatibility equations
  (1a,b), PDF pp.3-4; the sigma text inconsistency is retained in screening.
  <https://ntrs.nasa.gov/citations/19940019183>
- Guentert and Neumann, NASA TR R-33, 1959, radius source and averaged
  characteristic coefficients, printed pp.3-4/PDF pp.5-6.
  <https://ntrs.nasa.gov/citations/19980232191>
- Vanco and Goldman, NASA TM X-1502, 1968, planar sharp-throat full-nozzle
  equations, grid sensitivity and Table II; future comparison only.
  <https://ntrs.nasa.gov/citations/19680005278>
