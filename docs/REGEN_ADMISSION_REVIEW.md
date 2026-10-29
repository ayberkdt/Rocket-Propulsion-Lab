# Regenerative cooling input ownership and evaluation order

This review covers the fixed-channel/coolant records and the public
`/api/v1/regen-chamber` request boundary. It does not qualify the physical
solver or hardware, and it does not claim every thermal record has been audited.

## Observed defects

Known-invalid controls reached gas equilibrium before rejection: nonpositive
channel dimensions, invalid count/helix/roughness, invalid coolant inlet/flow,
unknown fluid/material/wall model, incompatible closeout and an invalid axial
flag. A chemistry sentinel established the ordering defect; 17 failing test
reports include the subtest cases rather than 17 unrelated implementation bugs.

`CoolingChannels` and `Coolant` also kept caller numeric objects and accepted
nonboolean flow-direction values. In Python, the string `"false"` is truthy;
it cannot be treated as a library-level counterflow flag. The JSON route's
documented textual boolean conversion remains at the request boundary.

## Repair

Channel and coolant records normalize annotated scalars using the common
finite-real/exact-integer/text contract, validate their positive dimensions
and inlet/flow values, and require an actual boolean library flow direction.
The route resolves inexpensive controls before invoking gas chemistry:

- positive pressure/throat and exit-to-throat area ratio at least one;
- bounded API count, positive fixed dimensions, roughness and helix range;
- channel-normal throat pitch `2*pi*(r_throat+t_wall)*cos(beta)/N`, requiring
  positive land width, using the solver's hot-wall-thickness convention;
- coolant key and the reference EOS's declared inlet temperature/pressure
  bounds; the actual EOS state still has to be solved later;
- an explicit positive coolant flow, or a deferred fuel-flow default;
- material conductivity availability, wall model, compatible closeout and
  boolean parsing.

The default coolant flow is computed from the real gas characteristic velocity
after chemistry. No placeholder mass flow, arbitrary new physical upper
bound, shortened station mesh, looser residual or disabled physics is used.
Pressure/temperature within catalog bounds does not prove a single-phase
state; the existing solver still evaluates the EOS and boiling/choking limits.

## Evidence scope

`test_api_regen_admission_order.py` checks invalid controls, one-ULP outside
methane's declared EOS range, incompatible wall controls and both axial and
helical pitch limits. Valid controls must reach the chemistry boundary.
The sentinel tests prove ordering only. `test_regen_input_records.py` checks
caller arithmetic ownership, direct malformed reals and actual boolean flow
direction. The existing real thermal, geometry, conjugate-wall and contract
suites retain the heat/momentum calculations and independent energy balances.

The unchanged two regen fuzz methods were profiled in isolated source/test
copies before and after the repair. These methods accept typed refusal or
strict finite JSON for semantic extremes; they are not all baseline-valid
physical fixtures. Timings and complete-module outcomes belong to their
source snapshots; a selected timing result is not a full-suite pass.

The governing channel/energy equations and source ranges remain those cited
by `regen_chamber.py` (NASA/TM-1998-206313, NASA SP-125 and the real-fluid
catalog's primary equations of state). Relevant-environment hardware,
measurement uncertainty and independent held-out data are still required
for hardware maturity evidence.

## Owned variable channel profiles

The same ownership contract now covers `ChannelProfile`. Ordered finite
columns are captured once, with typed field errors for malformed scalars or
unordered containers. Roughness, helix angle and scalar/per-row counts become
owned native values. Direct geometry/count coordinates are normalized before
use, including the constant-count shortcut. Finite count queries retain their
existing endpoint clamping and right-continuous bifurcation convention.

Interpolation forms the complete convex combination from exact representations
of the supplied floats before one output conversion. A finite interpolated
dimension must not be lost because an axial span overflows or a barycentric
weight underflows. The existing normal-span extent tolerance is retained;
unrepresentable spans use an exact bound comparison.

The retained baseline showed raw errors from double-consumed iterators, accepted
boolean/text dimensions, caller-owned count/helix arithmetic and incorrect
interpolation at synthetic extreme scales. In a real coupled methane-cooling
calculation, mutating an accepted helix object changed total heat load from
7429071.108330447 W to 7450866.895891617 W without replacing the profile. The new
contract tests require the coupled solution to remain unchanged, preserve affine
interpolation and count transitions, and check malformed-field refusals. These
are equation/ownership checks; they are not physical heat-load validation or
hardware qualification. Hot-gas records and other unreviewed boundaries remain
outside this scoped repair.

## References

- M. F. Wadel, NASA/TM-1998-206313, *Comparison of High Aspect Ratio Cooling
  Channel Designs for a Rocket Combustion Chamber With Development of an
  Optimized Design*, 1998. Channel/thermal model context:
  <https://ntrs.nasa.gov/citations/19980017619>.
- NASA SP-125, *Design of Liquid Propellant Rocket Engines*. Chamber heat
  transfer context already cited by the solver:
  <https://ntrs.nasa.gov/citations/19710019929>.
- U. Setzmann and W. Wagner, JPCRD 20, 1061-1155 (1991), methane reference
  Helmholtz equation and stated range. The code's reference catalog supplies
  the actual inlet bounds checked by this route:
  <https://www.nist.gov/system/files/documents/srd/jpcrd425.pdf>.
