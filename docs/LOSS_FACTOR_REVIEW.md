# Loss-factor basis, arithmetic and covariance review

This is a preliminary coefficient model. It does not select hardware, solve a
coupled reacting nozzle, establish measurement coverage, or demonstrate TRL 5.
The original full function/catalogue/claim review and hardware gates remain open.

## Operating basis and coefficient meaning

The legacy model multiplied both thrust and Isp by every factor while calling
discharge a throat mass-flow coefficient. Since Isp=F/(mdot*g0), equal F/Isp
multipliers imply unchanged mass flow. NASA SP-8120 defines discharge as an
actual/inviscid mass-flow ratio; the two meanings require different assumptions.

Version 2 retains the default `fixed-mass-flow-lumped` mode for existing callers,
with a critical warning and explicit source labels: its discharge argument is
a legacy performance multiplier, not physical throat Cd. All four factors
multiply both F and Isp. Their product is the thrust/Isp ratio; flow ratio is 1.

The selectable `fixed-pressure-area` mode holds the reference Pc*At product
fixed. Define E as intrinsic thermochemical cstar ratio, D as mass-flow Cd,
and B,V as net thrust-coefficient factors. Relative to the supplied reference:

- thrust ratio = B*V;
- mass-flow ratio = D/E;
- Isp ratio = E*B*V/D.

These follow the declared mdot=Pc*At*Cd/cstar convention and Isp definition.
E is not measured Pc*At/mdot efficiency with discharge already included.
Disabled factors are replaced by 1. B and V scale net Cf, including its pressure
contribution; they are not independently derived momentum-only corrections.

Independent coefficient variation can yield Isp above the supplied reference.
The result keeps that signed change and a critical compatibility warning. It
is a coefficient sensitivity scenario, not a prediction that exceeds an ideal
physical limit. Coupled gas properties, Cd/Cf consistency, actual operating
state, overlapping empirical losses and independent experimental validation
still require evidence. Sequential telescoping does not prove nonoverlap.

## Arithmetic and records

All inputs are strict finite real quantities; switches are actual booleans.
Disabled coefficients still receive domain validation. Exact binary rational
ledgers retain every intermediate product, reciprocal and subtraction. Each
term, total and ratio rounds independently. Nonzero mandatory exports that
would erase to zero or overflow fail with a field-labelled domain error.

Losses use exact state differences before export. Subtracting two displayed
thrusts is not equivalent, especially for a small loss on a very large force.
Zero balance residuals describe the underlying algebra, not independently
measured closure; rounded terms need not sum exactly to a displayed total.

Frozen inputs, correlation matrix, ordered terms, basis, policies, warnings
and assumptions are retained. Both terms and complete results replay from
those inputs. JSON mappings are copied, including nested matrices. Internal
consistency is not source/campaign authentication.
An additional attack after a passing full suite exposed retained warning text
subclasses whose equality callbacks could disguise forged content. Plain UTF-8
text normalization now precedes warning comparison for record and mapping
inputs; the earlier passing suite and failing new witness remain archived.
Combined API/native tests also exposed partially initialized optional dependency
imports. A deterministic harness witness proved that whole `sys.modules`
fixture rollback erased modules loaded by another thread. Fixture cleanup now
restores only its owned names. Failed native runs and fresh import controls are
retained; this repair is not a universal startup-stability claim.

## Uncertainty model and limits

`relative_uncertainty` is one common relative input standard uncertainty for
each enabled factor. The 4x4 correlation matrix follows combustion, discharge,
boundary layer, divergence order. The default is identity; unknown experimental
correlations have not thereby been demonstrated zero. Exact symmetry, unit
diagonal, bounded entries and positive semidefiniteness are checked without
tolerance. Singular matrices and exact cancellation are permitted.

The first-order relative variance is u² e^T R e. Elasticities for physical-mode
thrust, Isp and flow are (0,0,1,1), (1,-1,1,1), (-1,1,0,0); disabled entries
are zero. Legacy elasticities are (1,1,1,1) for F/Isp and zero for flow.
No enabled factors means zero factor uncertainty; ideal-reference uncertainty
is separately unmodelled. Exact quadratic forms and algebraic square-root
rounding preserve subnormal uncertainties and final band endpoints.

`one_sigma_relative` is retained as the thrust proxy alias. Separate Isp/flow
relative proxies are exported. Bands are mean +/- the linearized standard
uncertainty, not probability quantiles or qualified coverage intervals.
`coverage_qualified=false` and an explicit kind accompany every result.
Input uncertainty evidence, covariance provenance, nonlinear adequacy and
distribution/design-specific coverage remain unqualified, as required by the
broader objective. NIST's propagation formula is first-order and includes
covariance; a symmetric band alone supplies no coverage evidence.

## References

- NASA Glenn, specific impulse definition:
  <https://www1.grc.nasa.gov/beginners-guide-to-aeronautics/specific-impulse/>.
- NASA SP-8120, printed pp.10,19-20,93, discharge and flow effects:
  <https://ntrs.nasa.gov/api/citations/19770009165/downloads/19770009165.pdf>.
- RocketIsp project documentation 0.1.13, its declared cstar/Cd convention;
  this is not NASA endorsement or independent tool validation:
  <https://rocketisp.readthedocs.io/en/latest/efficiencies.html>.
- NIST TN1297 Appendix A, measurement models and first-order covariance:
  <https://www.nist.gov/pml/nist-technical-note-1297/nist-tn-1297-appendix-law-propagation-uncertainty>.
