# Characteristic axis and prescribed-wall boundary review

`solve_moc_axis_point` and `solve_moc_wall_point` supply the next boundary
processes needed by the constant-gamma characteristic solver. The legacy
mapped nozzle generator remains unchanged and unqualified. These processes
do not supply a sonic/throat starting line, a complete characteristic field,
wall optimization, shock/crossing treatment, validated exit profiles or
hardware evidence. All nine motor/feed hardware gates remain unmet.

The [axis-to-interior review](MOC_AXIS_DEPARTURE_REVIEW.md) documents a new
departure process and a genuine instability found when naively combining
the locally passing axis process with a two-family endpoint average. The
repaired near-axis stencil does not turn these local tests into full-field
stability or nozzle validation.

## Axis regularity

At a smooth axis, v/r tends to v_r, not generally zero. Cylindrical
continuity is `(rho*u*r)_x + (rho*v*r)_r = 0`. Symmetry and isentropic
axial momentum give

```text
(1-M*M)*u_x + 2*v_r = 0
v_r = (M*M-1)*u_x/2
axis characteristic source = a*a*(M*M-1)*u_x/2
```

Here velocities and sound speed are normalized by a*, the sonic speed.
The two previous axial velocities and the new unknown velocity define a
quadratic interpolant. Its derivative at the new x supplies u_x. This is a
derived discretization, not a free user-supplied gradient and not a literal
implementation of TR R-33's prescribed-centerline inverse design program.
The new velocity appears in the derivative and is solved in the same
compatibility equation. Planar mode has no cylindrical source.

Axis history requires exactly two strictly ordered, supersonic points with
r=0 and theta=0. The incoming C- point must be above the axis. The new stored
point must resolve x greater than the incoming point and latest axis point;
an unresolved marching topology is rejected. The symmetry values themselves
are imposed conditions, not measured residuals. Validated supersonic history
must already exist: no arbitrary Mach epsilon stands in for throat analysis.

## Finite wall slip

The incoming C+ is intersected with a prescribed finite straight upper-wall
segment. Velocity is constrained by `v = wall_slope*u`, while the same
Cartesian compatibility equation determines its magnitude. Both private and
stored coordinates must stay in `(x0,x1]` and progress downstream; the
segment is not extrapolated. Positive radii avoid an unsupported wall/axis
singularity. A curved wall needs separately controlled segmentation and
convergence; an optimized wall needs a full design procedure.

The processes preserve signed source terms in expansion and compression.
Endpoint averages and exact binary-input rational intermediate ledgers are
used. Platform trigonometric outputs remain binary inputs to those ledgers;
no rigorous trigonometric or integration-error enclosure is claimed. Raw
equation residuals, coordinate rounding allowances and adjusted decisions are
reported separately. Compatibility/slip residuals are not reduced by
coordinate allowances. The iteration threshold and 1..256 resource-budget
contract follow the interior process. This does not establish physical
accuracy or a global nozzle error bound.

## Retained results and independent checks

`MocBoundaryResult` owns the incoming/history or wall data and replays all
states, decisions and diagnostics. `source_term_per_m` is the normalized
Cartesian source coefficient (dimension 1/m); an unrepresentable optional
diagnostic is None with `source_term_status='unrepresentable'`. This does not
discard its private value in the solve. The qualification flags remain
actual False booleans and `discretization_error_bound` remains None.

Twenty-two methods cover constant states through M=1e200, six-step local
refinement in planar/cylindrical expansion, five-step signed compression,
110/180-digit independent fixture controls, coordinate scaling by 2**(+/-900),
finite-wall containment, stored topology, budgets and hostile replay.
Radial/spherical exact solutions obey rho*V*R**d=constant and their analytic
characteristic angles; these oracles do not call production MOC equations.
The axis source has a separate second-order approximation check. A sequence
of computed axis reflections reuses computed history rather than reseeding
every result from the analytic solution. This is a boundary-marching check,
not a fully coupled characteristic field.

Twenty curated faults exercise source limits/signs, quadratic history,
characteristic family, slip, budgets, topology and records. The first trial
had 18 assertion kills and two invalid mixed assertion/harness errors,
retained in `history/moc_boundary_first_fault_trial/manifest.json` (393 files).
In those two cases, subTest swallowed a missing-exception assertion before
the test read its nonexistent exception attribute. Nesting the exception
assertion inside the subTest fixes the harness; the producer kernel and
injected operators remain unchanged. Invalid runs are not counted as kills.

## References

- Aldo/Argrow Cartesian characteristic equations (1a,b), as screened in the
  [interior review](MOC_UNIT_PROCESS_REVIEW.md).
  <https://ntrs.nasa.gov/citations/19940019183>
- NASA TN D-4960, Appendix C, equation C1 and streamline definitions, printed
  p.30/PDF p.33: cylindrical continuity. This source is not a hardware
  validation of this implementation.
  <https://ntrs.nasa.gov/citations/19690005234>
- NASA TR R-33, Appendix B, printed pp.16-17/PDF pp.18-19: characteristic
  averaging and prescribed centerline context. Its inverse-program boundary
  specification differs from this forward quadratic-history closure.
  <https://ntrs.nasa.gov/citations/19980232191>
