# Supersonic axis-to-wall characteristic field

`solve_moc_field` assembles the interior, axis-reflection, axis-departure and
finite straight-wall processes into a connected supersonic field. The inlet
is a vertical section containing an owned computed axis state and at least
two radial intervals. Its wall endpoint must match the origin and slip angle
of the prescribed finite wall. Gamma and planar/cylindrical mode are inherited
from the axis result, including its two earlier axis states.

A separate [Hall throat approximation](HALL_THROAT_REVIEW.md) now locates the
sonic curve of a retained near-throat polynomial. Its nonzero nonlinear/slip
defects and degenerate sonic state do not supply this process's required
compatible supersonic inlet or owned downstream axis history.

This is forward analysis of an already supersonic strip. Sonic/throat
initialization, a curved or optimized wall, corner fans, shocks, boundary
layers, exit design and physical validation are still missing. The legacy
`generate_moc_nozzle` and HTTP contour generator remain unchanged; their
mapped-net, ineffective-axisymmetric-toggle and self-refinement findings
remain open. The original catalogue/claim review is incomplete, and all nine
motor/feed hardware gates remain unmet.

## Connected marching and geometric coverage

For each old complete front A (axis through wall), every adjacent pair
produces one half-front point B. The first pair uses the owned axis-departure
process; the others use C+ and C- interior intersections. B supplies the new
axis reflection, interior intersections and finite-wall slip intersection,
forming complete front C. The new axis history includes the computed previous
axis point. No analytic downstream values are injected or used to reset a
computed state. Both complete and half fronts are retained.

At n inlet intervals the strip between A and C contains n inlet triangles,
n outlet triangles, n-1 characteristic diamonds and two boundary triangles:
3*n+1 cells. Exact represented cross products require every cell to be strictly
convex with positive orientation. Front radii must increase; off-axis nodes
must stay above the axis and inside the prescribed straight wall. Wall nodes
receive only their explicit coordinate-rounding-cell allowance. Cell folds
and collapsed areas receive no epsilon allowance. Lower-level causal,
iteration, energy-branch, finite-wall and stored-topology failures propagate.
These checks do not locate a shock or establish general nonlinear stability.

An early 26-method green suite checked only diamonds and boundary triangles.
That face list covered half of the strip. It is retained in
`history/moc_field_pre_complete_tessellation/manifest.json` (397 files).
An independent regression now requires full area coverage and exact
cancellation of every internal oriented edge. The old face list fails it;
the source and failed assertion are retained in
`history/moc_field_incomplete_tessellation_failure/manifest.json` (397 files).
The repair adds both Cauchy-front triangle fans to geometry and mass audits;
it does not remove a face or replace the whole-field requirement with a
selected-cell check.

## Mass transport diagnostic

Uniform stagnation enthalpy and entropy, constant gamma and inviscid
irrotational flow are assumed, as in the local kernels. Dimensionless u,v
use sonic speed a*. Density is normalized by the inlet axis density:

```text
rho/rho_ref = [(2+(gamma-1)*M_ref^2)/(2+(gamma-1)*M^2)]^(1/(gamma-1))
section flux = integral rho * (u dr - v dx)                 planar
section flux = integral r * rho * (u dr - v dx)             cylindrical
```

The common planar mirror factor 2 and cylindrical factor 2*pi cancel in the
reported ratios. These are fluxes through the actual generally slanted
polylines, not rho*u times a nominal circular exit area. Density-weighted
velocity is interpolated linearly on each edge. Integrating that interpolant
with linear radius gives the cylindrical endpoint weights `(2*r_a+r_b)/6`
and `(r_a+2*r_b)/6`. The code retains rational products, sums and differences
of the represented values before forming exported dimensionless ratios.
It never rescales Mach or density to force closure.

`relative_mass_flux` retains every complete section relative to the numerical
inlet integral. `max_relative_mass_drift` is the maximum absolute deviation
from one. `max_cell_relative_mass_imbalance` is the largest absolute closed-cell
net flux divided by that cell's sum of absolute edge fluxes. The accumulated
absolute wall-edge flux divided by inlet flux is reported separately as
`wall_mass_leakage_relative`. Section drift, local imbalance and wall leakage
are different diagnostics; none is a proven error bound or measured balance.

Density uses an isolated 80-digit Decimal context, independently of the
caller's rounding mode and traps. The relative logarithmic density magnitude
is limited to 16384 as an explicit private-ledger resource policy; exceeding
it fails with `moc_field_density_resource`. Uniform extreme Mach can remain
representable through density normalization even when absolute density would
underflow binary64. The Decimal power and platform trigonometry are not
interval-enclosed. This numerical domain is not a validated gas/hardware
operating envelope.

## Discrete residuals and retained records

The result separately retains the maximum raw equation residual, maximum
coordinate-rounding-adjusted residual, maximum coordinate allowance and
maximum unenforced C+ defect of the cylindrical axis-departure process. The
Taylor departure closure replaces that C+ compatibility equation; its defect
is not included in the solved-equation convergence claim or forced to zero.

A constant-flow example translated to x=1e308 with transverse scale 1e295
has raw residual about 0.007, zero adjusted residual and coordinate allowance
about 0.0087. Its mass ratios remain one. The regression preserves all three
quantities, rather than reporting an exact geometric intersection.

Inputs and results own and replay all nested records and every diagnostic.
`physical_model_qualified`, `sonic_initialization_complete` and
`full_nozzle_solution` must be actual False booleans;
`discretization_error_bound` remains None. Supported resource limits are
3..65 inlet states, 1..128 marching steps, at most 4096 newly computed points
and 1..256 iterations per local process (default 128). The local threshold
remains 64*ulp(1). These limits are not accuracy guarantees.

## Independent checks and retained failures

Radial planar and spherical expansion/compression fixtures use
`rho*V*R**d=constant` and independent analytic Mach-line angles. They supply
only the inlet and old axis history, then compare every computed complete
and half-front state. Inputs are independently rounded at 110 and 180 digits.
The sampled 4/8/16 interval fields meet unchanged finest-grid limits of
1e-5 relative Mach error and 1e-4 degree angle error, while local cell imbalance
and state errors decrease under refinement. Global mass drift can have
pre-asymptotic cancellation; the sampled compression drift was not monotone
at every refinement. No universal order or stability theorem is claimed.

After the missing triangles were added, all four 16-interval cases failed
the existing 2e-6 cell-imbalance limit (about 2.56e-6 to 5.83e-6 in three
cases). That failed trial is preserved in
`history/moc_field_complete_mesh_accuracy_failure/manifest.json` (397 files).
The original bound is retained and a 32-interval field is additionally
required to satisfy it. Every former state, mass-drift and refinement
assertion remains; no face is dropped. Heavy planar and cylindrical
refinements use separate modules while retaining the canonical 30-second
collection, 120-second module and four-worker budgets.

An earlier wall-end fixture accidentally changed the wall slope and correctly
failed inlet slip before reaching the intended finite-segment check.
`history/moc_field_first_target_wall_fixture_failure/manifest.json` retains
the 21-pass/one-failure trial (398 files). The corrected fixture truncates
both endpoint coordinates along the same wall and checks the public
`moc_wall_outside_segment` code. Production input/slip tolerances were not
relaxed.

Independent high-precision quadrature checks the represented section fluxes
and maximum closed-cell imbalance. Other checks cover exact constant flow,
M=1e200, 2**(+/-900) scaling, 2**(+/-1000) cell orientation, large translation,
noncausal inlet refusal, bounded iterator consumption, caller Decimal context,
owned replay, forged evidence and explicit resource/domain failures.
Twenty-four curated faults target stencil links, cylindrical source/gamma,
axis history, mesh coverage/orientation, mass transport, residual disclosure
and input/result contracts. Only assertion failures between passing original
controls count as kills; the curated set is not exhaustive.

The first fault trial produced 20 kills and four invalid experiments: a failed
`assertRaises` inside a `subTest` was followed by reading the absent exception
context. Those errors were not counted as kills. The exception-code assertions
now remain inside the outer subtest, so an absent required refusal produces
an assertion failure without a secondary harness error. Production kernels,
fault operators and expected machine codes are unchanged. That trial and an
unindexed-reference failure are retained in
`history/moc_field_first_fault_and_reference_trial/manifest.json` (404 files).

## References

- Governing Cartesian compatibility and the cylindrical source convention:
  [interior review](MOC_UNIT_PROCESS_REVIEW.md),
  <https://ntrs.nasa.gov/citations/19940019183>.
- Cylindrical continuity and the smooth-axis limit: [boundary review](MOC_BOUNDARY_REVIEW.md),
  <https://ntrs.nasa.gov/citations/19690005234>.
- Owned Taylor closure and its separately visible unenforced C+ equation:
  [axis-departure review](MOC_AXIS_DEPARTURE_REVIEW.md).
- Conservation and isentropic mass-flow context:
  <https://www1.grc.nasa.gov/beginners-guide-to-aeronautics/mass-flow-rate-equations/>.

The staggered assembly, exact cell bookkeeping and polyline quadrature above
are implementation choices derived here. They are not a reproduction of a
NASA nozzle-design program, NASA consulting, endorsement or empirical evidence.
