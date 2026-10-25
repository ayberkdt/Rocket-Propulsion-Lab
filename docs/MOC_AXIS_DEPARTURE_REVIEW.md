# Axis-to-interior characteristic process

`solve_moc_axis_departure` accepts a replayed `MocBoundaryResult` of kind
`axis` and a positive-radius C- parent. It owns the complete axis history;
gamma and planar/cylindrical mode come from that history. A bare axis point,
wall result, free gradient or arbitrary epsilon radius cannot replace it.
This supplies a missing local marching process, not a sonic starting line,
complete nozzle field, optimized wall or hardware validation.

The [supersonic field process](MOC_FIELD_REVIEW.md) now uses this departure
in a connected axis-to-straight-wall mesh. Its complete cell and mass audits
retain the unenforced C+ defect separately; sonic initialization and complete
nozzle design/validation remain open.

## Retained failure and the resulting method

The first implementation used endpoint-averaged Cartesian compatibility for
both families, carrying the nonzero quadratic-history source at the axis.
It passed isolated radial and spherical local-error tests. A new coupled
axis -> interior -> axis walk instead exposed instability: final relative
error over the same radial-flow interval was about 2.19e-4, 5.55e-5 and
3.70e-5 at 4, 8 and 16 cells, then grew to 2.25e-2 at 64 cells.
The original 3..6 refinement-ratio assertion failed and was retained.
The source, tests and failed run are preserved in
`history/moc_axis_departure_first_coupled_failure/manifest.json` (401 files).
Isolated green checks were insufficient evidence for this coupled process.

In frozen uniform cylindrical flow, equal axis spacing dx and exact outer
inflow, let e_j be a small axis axial-velocity error. The original closure
has amplification polynomial `26*z**3 + 27*z**2 - 6*z + 1`; one root lies
outside the unit circle. This is our local analysis, not a general nozzle
stability theorem or a result attributed to a NASA program.

The repaired cylindrical departure instead uses a regular near-axis Taylor
relation, alongside C- compatibility and both characteristic geometries:

```text
v(x,r) = r * [w + w_x*(x-x_axis)] + O(r*(x-x_axis)**2 + r**3)
w = c*u_x,  c = (M**2-1)/2
w_x = c*u_xx + c_u*u_x**2
c_u = (gamma+1)*u/(2*a**4)
```

Smooth symmetry gives `u_r=0` and `v=0` on the axis; cylindrical continuity
gives w, while irrotationality `u_r=v_x` supplies the regular expansion.
The constant-gamma energy relation gives c_u. Both axial derivatives come
from the same owned three-point quadratic (two history states and the
computed axis state). Velocities and sound speed use a* normalization.
Exact represented rational products retain the derivatives even when their
optional displayed source exceeds binary64 range.

The C+ **geometry** is still iterated. Its endpoint-averaged compatibility
equation is replaced by the Taylor relation in this cylindrical stencil.
`averaged_plus_compatibility_defect` exposes that unenforced discrete
equation separately; it is not subtracted, forced to zero, included in the
convergence claim or concealed as a two-family compatibility pass.
Iteration residuals cover both geometries, C- compatibility and the radial
Taylor relation. Coordinate-cell allowances apply only to geometry. They
are not characteristic truncation bounds or measurement uncertainty.

The same frozen uniform-flow analysis of this repaired stencil gives
`13*z**3 - 4*z**2 + 11*z - 4`, whose roots are inside the unit circle.
That supports this specific mode analysis only. Arbitrary meshes, nonlinear
fields, corners and complete-nozzle stability remain unproved.

Planar flow has no cylindrical singular source and retains the existing
two-family interior process, including C+ compatibility as an enforced
equation. Applying the Taylor replacement to planar flow missed an unchanged
2e-7 local accuracy bound, despite fixing the cylindrical coupled test.
That failed trial is retained at
`history/moc_axis_departure_planar_taylor_accuracy_failure/manifest.json`
(396 files); the regular planar process was restored, not the bound relaxed.

## Verification scope and records

Twenty-four methods cover signed spherical/planar expansion and compression,
110/180-digit independent fixtures, local refinement, uniform flow through
M=1e200, 2**(+/-900) scale changes, an unrepresentable source diagnostic at
2**-1026, budgets, private versus stored topology, owned replay, policy
callbacks, and explicit false qualification flags. Independent radial
oracles use mass conservation and analytic characteristic angles without
calling production compatibility equations. Coupled tests reuse every
computed axis/interior state, with analytic data supplied only on the outer
boundary and in the first two history states.

The original coupled refinement assertion now passes. A longer 64-cell
screen gives radial error about 8.83e-8, and a separate constant-flow walk
does not amplify a 1e-10 initial Mach perturbation in the sampled cases.
These finite checks do not prove a global stability or discretization bound.
The old kernel is loaded from its preserved snapshot in an isolated
screening process to reproduce the failure against the same oracles.

Twenty-one curated faults target axis gradients, curvature/nonlinear terms,
private source retention, family/source signs, planar routing, coefficients,
budgets, topology, visible defects and result contracts. Assertion failures
count only between passing original controls; errors are invalid runs.
This is selected semantic coverage, not an exhaustive mutation score.

A later 1366-method/609-fault green run still missed a branch witness:
quadratic axis history can have zero first derivative but nonzero curvature.
The original zero-source constant-state shortcut refused that valid local
input. The green run and the refusal are retained at
`history/moc_axis_departure_pre_curvature_guard/manifest.json` (642 files).
The constant branch now requires identical Mach throughout the owned axis
history, with a new curved-history regression and a fault that restores the
old shortcut. This does not relax an error bound or hide the earlier miss.

No iteration threshold, resource budget or original local/coupled accuracy
assertion was relaxed. Results replay from owned inputs and keep actual
False `physical_model_qualified` and `full_nozzle_solution` flags, with
`discretization_error_bound=None`. Platform trigonometry is not enclosed.
This process is not connected to the legacy HTTP nozzle generator. The
mapped-net counterexample, complete field/exit validation, sonic seed,
broader function/catalogue review and all nine hardware gates stay open.

## References

The governing Cartesian equations, cylindrical continuity and primary-page
screening are retained in the [interior](MOC_UNIT_PROCESS_REVIEW.md) and
[axis/wall](MOC_BOUNDARY_REVIEW.md) reviews. The Taylor departure closure and
the two amplification polynomials above are derived for this implementation;
they do not claim to reproduce NASA TR R-33's inverse prescribed-axis program
or to establish NASA consulting, endorsement or empirical validation.

- Cartesian governing equations (1a,b), screened with their printed sigma
  inconsistency in the interior review:
  <https://ntrs.nasa.gov/citations/19940019183>
- Cylindrical continuity, NASA TN D-4960, Appendix C, equation C1, printed
  page 30/PDF page 33:
  <https://ntrs.nasa.gov/citations/19690005234>
- NASA TR R-33, Appendix B, prescribed-axis inverse-program context, not the
  forward Taylor closure implemented here:
  <https://ntrs.nasa.gov/citations/19980232191>
