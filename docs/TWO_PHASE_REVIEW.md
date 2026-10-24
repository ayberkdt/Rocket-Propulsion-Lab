# Conditional two-phase exit-state review

This is prescribed exit-state accounting, not a coupled particle/gas solver,
validated physical bound, statistical interval, hardware selection or TRL proof.
The original full review, physical validation and hardware gates remain open.

## Reference state and coefficient meaning

The previous implementation scaled total thrust by a particle momentum factor.
With nonzero pressure thrust this also scales a term outside that mixture sum.
The NASA-hosted Drewniak et al. paper separates gas/particle momentum from pressure
thrust. Let F0,I0 be reference total thrust and Isp, P specified pressure force,
M=F0-P>0 reference momentum, chi condensed exhaust mass fraction and v particle/gas
axial exit velocity ratio. The conditional equations are:

eta_m = (1-chi) + chi*v

F = M*eta_m + P

Isp = I0*F/F0

Total mass flow and gas exit state are fixed. F0/I0 are reference values, not
asserted thermochemical maxima. P defaults to zero for old calls, with a visible
warning: zero is an assumption, not observed matched expansion. Negative P may
produce zero/negative conditional thrust; those endpoints are retained and
warned. Actual nozzle/backpressure compatibility remains unsolved. All-condensed
exhaust remains a warned mathematical limit, not a realizable carrier-gas claim.

`exit-velocity-ratio` is the default linear convention. In
`specific-kinetic-energy-ratio`, k=(vp^2/2)/(vg^2/2) **per unit mass**, so v=sqrt(k).
It is not total particle-phase/gas-phase energy. The legacy coupling-efficiency
argument names remain aliases for these ratios in [0,1]. Separate particle
velocity, momentum multiplier and total-thrust-ratio ranges are exported.
`performance_efficiency_range` is retained as a signed total-thrust-ratio alias.

## Arithmetic and retained records

Strict real inputs reject booleans, nonfinite/overflow values and erased nonzero
quantities. Numeric text is parsed explicitly by the API form adapter, while
direct Python quantity fields require real numbers. Coupling definitions are
validated text, not stringified controls. Source endpoints must be ordered.

Binary rational ledgers keep pressure and momentum terms exact until each
export. Private momentum overflow can therefore cancel to finite endpoints.
Velocity-mode arithmetic is rational. Energy mode rounds every final affine
square-root expression independently using exact squared comparisons on binary
rounding cells. Displayed velocity/force is not reused for dependent outputs.
The nonpositive-thrust flag uses an exact sign comparison. Genuine zeros and
signed values are supported; unrepresentable mandatory nonzero exports fail
typed. No absolute tolerance or saturation substitutes for the equations.

Frozen inputs, policies, assumptions and warnings replay every result field.
Fixed pairs read at most three entries; warning collections read at most five
because this schema emits at most four warnings. Warning text is normalized
before comparison, including mappings, so subclass callbacks cannot disguise
forged provenance. Known malformed mapping failures are typed; unrelated
programmer RuntimeError callbacks propagate. Replay is not authentication.

## Qualification and remaining physical work

`band_kind="conditional_exit_state_scenario_range"`,
`physical_bound_qualified=false` and `statistical_coverage_qualified=false`
are retained facts. Unknown particle distributions, thermal lag, phase change,
and gas reaction to particle acceleration are excluded. Kibbey's coupled-model
discussion identifies the gas momentum debit when particles accelerate. A
fixed-gas-state ledger cannot supply that missing conservation solution by
passing arithmetic tests. Coupled model and independent validation remain tasks.

Five published zero-gimbal air/water summary rows were visually checked in the
preceding stage. Their RMS oscillation amplitude is not measurement uncertainty.
They measure pressure response and cannot directly validate these fixed-exit-state
thrust scenarios. Gas-state/geometry data, uncertainty pedigree and an applicable
coupled model are needed. Calibration data must not become validation data merely
by reuse. The data are not evidence for a selected motor/feed system.

Metadata advances to 2.0. Ordinary calls retain their defaults; old manually
constructed result records lacking retained inputs fail typed. This migration
does not close physical validation or any original hardware requirement.

## References

- Drewniak, Boyer and Yetter, 2022, Eq.(3), separate momentum and pressure:
  <https://ntrs.nasa.gov/citations/20220006969>.
- Kibbey, 2012, coupled flow and gas momentum debit:
  <https://ntrs.nasa.gov/citations/20120015678>.
- Whitesides et al., 1995, air/water experiment, Table I:
  <https://ntrs.nasa.gov/citations/19970015335>.
