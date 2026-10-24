# Prescribed bulk coolant balance review

This is a single-phase, constant-property bulk energy and major-pressure-loss
calculation with equal flow in identical channels. It does not predict wall
temperature, conjugate heat transfer, boiling, cavitation, material survival,
manifold distribution, or relevant-environment motor/feed hardware performance.
The original full review and hardware gates remain open.

## Energy and decision ledgers

With prescribed absorbed power H, total mass flow m and cp:

    capacity = m*cp
    rise = H/capacity
    outlet = inlet + rise

Every input quantity is a finite normalized real; bool is rejected. Channel
count is an actual positive integer. Positive heat and the existing restriction
maximum_temperature>inlet remain part of this model's interface. The numerical
ledgers retain exact binary input products and ratios before each export.

Temperature-limit exceedance compares the exact outlet expression with the
supplied limit. Comparing displayed temperatures can hide a real exceedance.
The critical warning and temperature_limit_exceeded flag use that exact decision.
The limit is a user requirement, not an inferred saturation or material limit.

temperature_limit_margin_k is maximum minus exact outlet, independently rounded.
An unrepresentable nonzero optional margin is null, with
temperature_limit_margin_unrepresentable=true; it does not invalidate other
representable mandatory outputs or suppress the exact decision.

energy_balance_residual_w=0 is labelled
bulk_energy_ledger_identity_not_measured_energy_closure. The separate
displayed_temperature_energy_residual_w is the absolute discrepancy obtained
by inserting the displayed outlet into the exact capacity ledger. It diagnoses
output rounding, not measurement/model agreement. Its nonzero unrepresentable
values also use null and an explicit flag. A positive rise can coexist with an
unchanged displayed outlet; no heat or rise is silently erased.

## Hydraulics and cross-section assumptions

Let N be channel count, A single-channel area, rho density, mu viscosity, D
hydraulic diameter and L length. No intermediate floating total-area product
or exported velocity is needed:

    velocity = m/(rho*N*A)
    Re = m*D/(N*A*mu)
    pressure_drop = f_D*(L/D)*rho*velocity^2/2

The laminar branch uses an explicit laminar_darcy_poiseuille_number K:

    f_D = K/Re

Default K=64 assumes a circular equivalent. It is not a universal consequence
of supplying a hydraulic diameter. NASA GFSSP v6, printed p.40 (PDF p.80),
equation 28, uses f_D=4*Po/Re with cross-section-dependent Po. This input K is
therefore Darcy f*Re, equal to 4*that Po; it is not a Fanning friction factor.
Area, hydraulic diameter and K are prescribed rather than independently derived
from a verified cross-section. No shape validity or flow distribution is inferred.

The legacy branch policy remains explicit: Re<2300 is laminar; [2300,4000) is
labelled transition and retains an unqualified Swamee-Jain extension; Re>=4000
uses the turbulent formula. Branches use exact Re even when its exported value
rounds to a threshold. This is not EPANET's transition policy: EPANET 2.2 uses
its own 2000/4000 switches and a transition interpolation. That implementation
is consulted here to verify the formula/logarithm convention, not to establish
equivalence or a qualified cooling-channel regime model.

    a = roughness/(3.7*D) + 5.74/Re^(9/10)
    f_D = 0.25/log10(a)^2

Decimal coefficients/exponents are exact rationals for numerical evaluation.
Directed positive power and logarithm bounds certify each finite output. Exact
rational roots and decimal-power logarithms handle recognized exact cases.
The a=1 singularity raises singular_swamee_jain_argument. Representable values
near the pole are not clamped or treated as singular based on a rounded a.
For a>1 the unsquared inverse-friction root is negative; the squared formula
extension is retained with a critical warning and no physical qualification.
Unresolved comparisons or rounding raise ConvergenceError. Pressure loss is
computed from unrounded f_D and dynamic head, rather than their exports.

friction_correlation_domain_qualified=false and physical_model_qualified=false
remain explicit. The original Swamee-Jain article's full text and an independently
accepted correlation-error envelope have not been recovered in this review.
Primary implementation checks do not establish experimental applicability to
these channels. The model includes major friction loss only; entrance, bend,
manifold, acceleration and changing-property losses remain absent.

## Records, API and evidence scope

Owned inputs and the complete result replay every output, optional diagnostic,
regime, warning, numerical policy and boolean qualification flag. Warning text
is normalized before comparison. Known malformed mappings fail typed, while
unrelated programmer RuntimeError callbacks propagate. Replay is internal
consistency, not authentication of measurements or hardware configuration.

The nested cooling result declares model_version 2.0. The parent nozzle-thermal
metadata continues to identify its gas-side model version 2.0. The API accepts
the explicit laminar coefficient, including numeric form text, and returns the
exact-limit flag even when displayed temperatures coincide. The UI shows that
decision and the displayed-temperature energy discrepancy separately.

Independent tests use exact Fraction energy and laminar equations and separate
180/360-digit Decimal turbulent references. Witnesses include hidden true limit
exceedance, equal-limit controls, exact Reynolds thresholds, subnormal heat,
optional diagnostic underflow, large exact channel counts, recoverable private
overflow, near/exact logarithm poles, unrounded-capacity/friction dependence,
cross-section factors and forged records. These finite checks do not qualify
the physical model, cover every possible input, or establish TRL 5.

## References

- NASA GFSSP v6 user manual, printed p.40, cross-section-dependent laminar factor:
  https://www.nasa.gov/wp-content/uploads/2024/04/gfssp-v6-usermanual.pdf
- EPANET 2.2 primary implementation, frictionFactor and A9=-2/ln(10):
  https://raw.githubusercontent.com/OpenWaterAnalytics/EPANET/v2.2/src/hydcoeffs.c
- Swamee and Jain, 1976, primary publisher abstract and bibliographic identity:
  https://doi.org/10.1061/JYCEAJ.0004542
- Python Decimal arithmetic:
  https://docs.python.org/3.12/library/decimal.html
