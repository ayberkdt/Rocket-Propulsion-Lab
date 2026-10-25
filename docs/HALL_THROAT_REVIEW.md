# Hall throat approximation and bracketed sonic curve

`HallThroatModel`, `evaluate_hall_throat` and `locate_hall_sonic_point` provide
a first/second-order constant-gamma throat approximation and the sonic curve
of its retained polynomial. They do not provide transonic CFD, a complete
sonic-to-MOC startup, physical validation or motor/feed TRL5. The original
function/catalogue/claim review and all nine hardware gates remain open.

## Equations and reference geometry

With rt the throat radius/planar half-height, R=rc/rt, eps=1/R, sigma=0 planar
or 1 cylindrical, the stretched-coordinate definitions are

```text
k = (sigma+1)*R/(gamma+1),  w = sqrt(k)*x/rt,  s = (radius/rt)^2
u/a* = 1 + eps*U1(w,s) + eps^2*U2(w,s)
v/a* = eps*sqrt(s/k) * [H1(w,s) + eps*H2(w,s)]
```

The polynomials were independently reconstructed from the full isentropic,
irrotational potential-flow equations. Planar coefficients agree with
McCabe's reproduced Hall results (Appendix A); Cartesian cylindrical
coefficients agree with Kliegel/Quan equations 4-1 through 4-4. The symbolic
reconstruction is retained in
`output/assurance/hall_second_order_symbolic_derivation.json`.

The explicit reference wall is hyperbolic:
`(radius_wall/rt)^2 = 1 + (x/rt)^2/R = 1 + eps*w^2/k`.
`radial_basis='throat'` uses radius divided by rt; `'wall'` uses the fraction
of the actual local wall radius. Inputs must lie inside this wall. No arbitrary
rounded or thrust-optimized contour is supplied. Prescribed wall curvature
and the truncated flow's streamline curvature are not asserted to agree exactly.

Cartesian velocities are retained to the selected order. The implementation
does not multiply untruncated speed/angle series and call that the same
approximation. Squared radial velocity, energy and polynomial coefficients
use exact rational ledgers of normalized binary64 inputs. Radical factors
remain symbolic until export. This avoids prematurely overflowing rc/rt or
erasing its reciprocal. Platform atan2 is not interval-enclosed.

## Visible model defects and display limits

The sound speed uses the full constant-gamma energy relation. The full
nonlinear potential equation is evaluated on the approximation; it is not
replaced with the vanishing formal-series coefficient residual.
`potential_residual` divides its magnitude by the leading Hall axial sonic
gradient, with velocities normalized by a* and lengths by rt.
`relative_term_balance` separately divides by the sum of absolute equation
terms. Near the sonic axis that latter ratio can approach one as its leading
terms vanish; the gradient-normalized defect can still be small. Both are
retained. `irrotational_residual` is the exact radial coefficient of the
derivative identity, including its axis limit.

At the actual wall, the magnitude of `v-u*wall_slope` is retained in sonic-speed
units. It is zero at the geometric throat, generally nonzero away from it.
Its status distinguishes `not_wall`, `finite` and `unrepresentable`.
Velocity/perturbation loss and an unrepresentable wall defect set
`representation_limited`; a displayed Mach of one cannot overwrite the private
subsonic/supersonic classification. Mandatory nonzero coordinates outside
output range fail. Optional unrepresentable quantities retain None/status;
nonzero residual upper bounds round upward even below the subnormal range.
Leaving the forward positive-energy branch fails without clipping.

A real missed witness was found after a 29-method green target: w=2**-1070,
rt=1, rc=5, gamma=1.4 at the wall retains finite nonzero radial velocity, but
its nonzero wall defect is below half the least binary64 subnormal. Its
status was correct while the overall representation flag was false. A
700-digit independent regression now catches it. Only flag aggregation changed;
coefficients, energy, residuals and root rules were retained. The earlier
green run and failing witness are preserved in
`history/hall_throat_hidden_wall_defect_failure/manifest.json` (406 files).

## Numerical sonic certificate

For wall fraction eta in [0,1], the condition is `u^2+v^2=1` with
`s=eta^2*(1+eps*w^2/k)`. Its rational polynomial has degree at most ten.
Eta=0 reaches the axis; eta=1 reaches the actual hyperbolic wall, not merely
the throat height. No 1.001 Mach or sound-speed offset is introduced.

On the supported search interval w in [-1,1], exact Bernstein coefficients
and at most eight subdivisions must certify strictly positive derivative.
The root must also be bracketed. This sufficient computational certificate
does not establish that every possible Hall branch can be found. Unbracketed
and uncertified cases fail with distinct codes. Bisection's default budget is
128, explicitly bounded to 1..256; convergence requires its internal endpoints
to round to the same binary64 parameter or an exact dyadic root.

Parameter and polynomial-function bounds are exported outward. Both x and
radius bounds enclose the retained curve point. Radius bounds include the
interior minimum when the parameter interval crosses zero. A 25-method
pre-extension run exposing only axial bounds is retained in
`history/hall_throat_pre_radial_bracket_extension/manifest.json` (403 files).
The root record was extended to enclose both spatial coordinates.

This brackets numerical root error inside the retained model, not Hall
truncation error, gas uncertainty or actual sonic position. Inputs/results own
and replay all fields. `physical_model_qualified` and
`sonic_initialization_complete` must be actual False booleans;
`approximation_error_bound` remains None.

## Primary comparisons and remaining work

`docs/assurance/hall_throat_theory_reference.json` retains twelve printed Hall
entries at R=5 and R=10. Eleven agree within their original four-decimal
rounding cells. Table 4-2's first-order gamma=1.4, 1/R=0.1 cell prints 0.0856,
whereas equation 4-1 gives 0.086602540378... . Full-page and magnified visual
inspection confirm the printed cell. The discrepancy is retained and flagged;
the code is not fitted, the tolerance not widened, and the cell is not counted
as a successful comparison. These are theoretical values, not measurements.

Independent high-precision equations check velocities, Mach and coordinates;
independent differentiation checks the nonlinear PDE diagnostic. At 180 digits,
sonic roots and polynomial-function values must lie within exported intervals.
Other checks cover 110/180-digit reference stability, 700-digit perturbations,
2**(+/-900) geometry scaling, 2**(+/-1000) radical bounds, exact zero roots,
resource limits, ownership, forged records and source anomaly retention.
Selected fault injection checks corrupted coefficients, equations, bounds,
input guards and evidence claims; it is not an exhaustive mutation score.

Kliegel/Quan, printed p.4-3, discusses Hall's local applicability and curvature
ratios of five or more. That theoretical comparison is not a universal
accuracy bound. Our sampled second-order PDE defect improves asymptotically
relative to first order; this does not validate arbitrary contours, tight
curvature, nonlinear startup, engine conditions or a complete nozzle.

The present sonic points retain nonlinear PDE and wall-slip defects. At sonic
speed the existing x-parameterized MOC process is degenerate. A nonlinear
startup/corrector, compatible curved-wall treatment and a verified connection
to the downstream field remain required. The legacy mapped/HTTP nozzle and
all prior open construction findings remain open.

## References

- McCabe, *Design of a Supersonic Nozzle*, ARC R&M3440, report March 1964,
  issued 1967; Appendix A, printed pp.12-13/PDF pp.13-14.
  <https://reports.aerade.cranfield.ac.uk/handle/1826.2/4023>
- Kliegel and Quan, *Convergent-Divergent Nozzle Flows*, 02874-6002-R000,
  30 December 1966; full potential equations, pp.2-1..2-3; Hall comparison,
  pp.4-1..4-6/PDF pp.30-35.
  <https://ntrs.nasa.gov/api/citations/19670009614/downloads/19670009614.pdf>
- Sauer, NACA TM1147, June 1947; historical curved-critical-line context,
  not the implemented second-order coefficients.
  <https://ntrs.nasa.gov/api/citations/20050019415/downloads/20050019415.pdf>

Primary PDFs, inspected pages, hashes and the anomalous-cell detail are pinned
under `output/assurance/sources/sonic_start`. No NASA consulting, endorsement,
engine experiment or TRL5 qualification was obtained.
