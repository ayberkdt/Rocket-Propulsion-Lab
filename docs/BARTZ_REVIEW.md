# Prescribed signed Bartz profile review

This module evaluates a conditional empirical formula with supplied constant
stagnation gas properties and a prescribed uniform wall temperature. It does
not solve combustion, reacting boundary layers, radiation, wall conduction,
coolant feedback, material temperatures, or hardware qualification. The original
function/catalogue/claim review and all motor/feed hardware gates remain open.

## Source equations and assumptions

NASA SP-125, printed pp.100-101 (PDF pages 109-110), gives the signed heat-flow
relation and the coefficient and property correction in equations 4-10, 4-13
and 4-14. In consistent SI units the force/mass conversion factor is 1. Writing
D=2*throat_radius, Ar=A/At, C=1+(gamma-1)*Mach^2/2 and
B=(Tw/Tc*C+1)/2, the implemented coefficient is:

    h = 0.026*mu^(1/5)*cp*Pr^(-3/5)*D^(-1/5)
        *(Pc/cstar)^(4/5)*(D/Rc)^(1/10)*Ar^(-9/10)
        *B^(-17/25)*C^(-3/25)

The 2012 NASA Glenn Mars hopper paper, printed p.8, equations 1 and 4, supplies
the selected turbulent recovery-temperature convention:

    recovery = Pr^(1/3)
    Taw = Tc*(1+recovery*(C-1))/C
    q = h*(Taw-Tw)

Positive q means heat into the wall. Tw below Tc does not imply Tw below local
Taw. Negative q is retained, rather than silently replaced by zero. The existing
Tw<Tc input condition remains an explicit restriction of this prescribed model;
it is not a universal physical law. Constant cp/mu/Pr and prescribed Mach/area
states are user assumptions, not a CEA calculation or an isentropic consistency
proof. The printed decimal coefficient and exponents are treated as exact
rationals for numerical evaluation, not as exact physical laws.

SP-125 explicitly discusses departures due to radiation, dissociation, mixing,
instabilities and deposits. The historical diameter and pressure warning
thresholds remain labelled legacy screening thresholds; this review has not
established an experimental validity domain for them.

## Surface reconstruction and numerical evaluation

For each straight meridional segment let length L, endpoint radii r0,r1 and
endpoint signed fluxes q0,q1 be given. Both q and radius are reconstructed
linearly in meridional position. The exact integral of that reconstructed field is:

    Qsegment = pi*L/3 * [(2*r0+r1)*q0 + (r0+2*r1)*q1]

The former mean-flux times frustum-area formula instead integrated the product
of means and omitted a cross term. Independent two-node Gaussian quadrature
of the reconstructed quadratic checks the new integral. This is not a bound on
the unsampled true profile. quadrature_error_bound_w is null.

All binary64 inputs enter exact rational ledgers. Positive fractional powers
use exact rational roots where available, otherwise directed logarithm and
exponential bounds. Geometry uses a declared fixed binary64 math.pi; bounds do
not include its difference from mathematical pi or physical model error.
Each mandatory export rounds independently after both bounds agree. Signed
nonzero underflow/overflow fails typed; exact zero flux and initial cumulative
zero are retained. Unresolved signs, peaks or rounding fail with ConvergenceError.

Private Pc/cstar, Mach^2 and segment-length squared may exceed binary64 without
making final outputs unrepresentable. Displayed Taw and prior cumulative heat
are not reused in dependent calculations. Peak location is selected from
unrounded flux comparisons; a tie in displayed values does not decide location.
Identical prescribed station states tie at the first station. Other unresolved
comparisons are not guessed.

## Inputs, records and user-facing scope

The thermal snapshot retains the consumed throat radius and station x, radius,
area ratio and Mach. Other NozzleContour metadata is not used to authenticate
or qualify this prescribed flow. Values are finite normalized reals; booleans
are rejected. A surface integral requires at least two axially ordered stations.
Positive radii and area ratios must agree within their exact binary64 rounding
cells with Ar=(r/rt)^2. This is a representation-coherence guard, not measurement
uncertainty, a geometry calibration or an isentropic solver; values are not clamped.

Owned inputs and complete results replay every station, integral, peak, warning,
policy and qualification flag. HeatFluxStation is finite station data; equation
and cumulative-integral replay belongs to its enclosing complete result. Mapping
round trips snapshot the source; warning text is normalized before comparison.
Replay is internal consistency, not experiment/campaign authentication.

physical_model_qualified is false. The legacy energy_balance_residual_w=0 is
labelled integral_ledger_identity_not_measured_energy_closure. No independent
energy measurement is part of that identity. The optional bulk coolant model
currently requires positive net absorbed heat and refuses nonpositive net loads
with unsupported_nonpositive_bulk_heat; it does not rectify signed heat. Its
own constant-property/hydraulic assumptions require separate review.

## Independent experimental evidence remains incomplete

The Mars hopper paper's Table VIII compares a particular GO2/GH2 heat-sink
experiment against a coupled finite-volume/CEA model. That is not a comparison
with this implementation: matching gas properties, complete geometry, boundary
conditions, raw measurement traces and measurement uncertainty have not been
recovered. The local evidence folder preserves the original PDFs, visually
inspected equation pages, a transcription of that table and these limitations.
The table cannot establish the current code's physical accuracy or transfer to
an unselected motor/feed design. Calibration and independent validation must
remain separate when sufficient data are available.

## References

- NASA SP-125, Design of Liquid Propellant Rocket Engines:
  https://ntrs.nasa.gov/citations/19710019929
- Original NASA SP-125 PDF, equations on PDF pages 109-110:
  https://ntrs.nasa.gov/api/citations/19710019929/downloads/19710019929.pdf
- Linne et al., Mars Surface Mobility Leading to Sustainable Exploration,
  AIAA 2012-0806, equation and comparison table on printed p.8:
  https://ntrs.nasa.gov/api/citations/20120004103/downloads/20120004103.pdf
- Python Decimal arithmetic:
  https://docs.python.org/3.12/library/decimal.html
