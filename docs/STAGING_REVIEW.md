# Conditional finite-grid staging review

This is a constant-Isp, prescribed-structure mass model. It does not size flight
hardware, model gravity/drag/steering losses, select a mission, or establish TRL 5.
The original catalogue/function/claim review and hardware gates remain open.

## Equations and finite-grid optimum

For each stage let U be the upper-stack mass, epsilon=dry/(dry+propellant),
c=g0*Isp, d the allocated delta-v, R=exp(d/c), and s=exp(-d/c). Solving the
ideal rocket equation gives:

- stage wet mass = U*(1-s)/(s-epsilon);
- initial stack mass = U*(1-epsilon)/(s-epsilon);
- dry mass = epsilon*wet; propellant mass = (1-epsilon)*wet;
- finite positive stage capacity requires s>epsilon.

Names and input quantities are normalized literal text and binary64 reals.
The grid is exact target*k/resolution, before any public delta-v rounding.
Its integer counts sum exactly to resolution. The achieved delta-v field
describes that exact ledger; separately rounded displayed stage values can
sum to a neighbouring float. No tolerance repairs the target or capacity.

The initial mass is payload times the product of stage mass factors. This
review derives the following optimizer property for this particular model;
it is not attributed to NASA or to experimental validation. With
f(d)=log((1-epsilon)/(exp(-d/c)-epsilon)):

    f'(d) = 1/[c*(1-epsilon*exp(d/c))]
    f''(d) = epsilon*exp(d/c)/[c^2*(1-epsilon*exp(d/c))^2] > 0

Thus each stage's successive equal-grid marginal costs increase strictly.
An exchange argument proves that selecting the least remaining marginal cost
at each unit minimizes the sum over every feasible integer allocation: a
solution using a larger available increment instead can be improved by swapping
that increment with a smaller one, while the prefix property is preserved.
The implementation compares the equivalent positive next/current mass-factor
ratios without rounded logarithms. Exactly equal cross-products are recognized
by collecting rational coefficients of identical rational exponents. Exact
ties allocate the next unit to the lower stage first. Other comparisons require
separated directed bounds; unresolved cases raise ConvergenceError.

This proves a conditional grid property. No continuous optimum or general
mission optimizer is claimed. A grid can be infeasible while a different grid
or the continuous model remains feasible; the error is infeasible_staging_grid.

## Arithmetic, records and reporting

Input products, grid fractions and exponents remain exact binary rationals.
Directed Decimal exponential bounds and interval arithmetic retain unrounded
upper-stack mass through every stage. Each mandatory field rounds independently
only when both bounds agree. Nonzero exports that cannot remain finite and
nonzero raise DomainError. A mass ratio rounding to 1 does not erase propellant,
and one rounding to 1/epsilon does not decide capacity. Precision budgets are
finite; unresolved rounding fails typed instead of claiming a result.

StagingInputs, stage solutions and the complete result retain owned inputs.
Results replay allocations, stage identity, every mass, target, grid diagnostic,
scope, numerical policy and actual boolean flags. Replay is internal consistency,
not campaign authentication or hardware validation. Maximum-size collections
are read with bounded cardinality checks.

The legacy convergence_delta_initial_mass_relative field is only the absolute
relative difference between the requested and comparison grids. A comparison
grid capped at 800 can be non-nested. At resolution=800 no second grid is solved:
the difference is null and refinement_performed=false. Zero from the same grid
would supply no convergence evidence. continuous_optimum_qualified is false.
Only zero delta-v allocations are labelled active constraints; proximity to a
strict capacity boundary is not an active finite-mass constraint.

The frontend records the solved inputs and full result artifact, even if a form
is later edited, and retains version 2.0 and the scope qualification.

## Verification scope

Independent tests use direct exponential equations at 1400 and 2200 decimal
digits for final exports. Separate 180-digit exhaustive enumeration checks every
allocation for two- and three-stage examples. Symmetric two/five-stage controls,
near-equal Isp, huge private products, the last float below capacity, an infeasible
odd grid with a feasible even-grid control, signed input attacks and forged
records exercise distinct failure modes. These are finite selected witnesses,
not an empirical qualification or an exhaustive test of every binary64 input.

Structural fractions, Isp dependence on operating state, separation dynamics,
flight constraints, independent measured comparisons and relevant-environment
hardware demonstrations remain unvalidated.

## References

- NASA ideal rocket equation, used only for the ideal per-stage relation:
  https://www1.grc.nasa.gov/beginners-guide-to-aeronautics/ideal-rocket-equation/
- Python Decimal arithmetic and correctly rounded exponential operations:
  https://docs.python.org/3.12/library/decimal.html
