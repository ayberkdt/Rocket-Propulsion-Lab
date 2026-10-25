# Finite-rate chemistry in the nozzle

This work package adds one-dimensional finite-rate (kinetic) nozzle expansion,
the class of loss that NASA's ODK/TDK programs estimate, between the frozen and
shifting-equilibrium limits that the repository already computes. Everything is
pure Python (standard library) in three new modules:

- `rocket_propulsion.thermochemistry.kinetics`: elementary reactions, rates,
  equilibrium constants, production rates, analytic Jacobian and the RODAS stiff
  integrator.
- `rocket_propulsion.thermochemistry.kinetics_mechanisms`: the H2-O2 subset of
  GRI-Mech 3.0, transcribed line by line.
- `rocket_propulsion.thermochemistry.nozzle_kinetics`: steady quasi-1D
  finite-rate expansion through a prescribed area `A(x)` and Bray's
  sudden-freezing estimate.

The modules are not exported from `rocket_propulsion.thermochemistry`; import
them directly.

```python
from rocket_propulsion.thermochemistry import get_propellant_pair, reactants_for_mass_ratio
from rocket_propulsion.thermochemistry.nozzle_kinetics import (
    HyperbolicNozzle, bray_sudden_freezing, solve_nozzle_kinetics)

reactants = reactants_for_mass_ratio(get_propellant_pair("LOX_LH2"), 6.0)
nozzle = HyperbolicNozzle(throat_radius_m=0.05, contraction_area_ratio=3.0,
                          exit_area_ratio=40.0, convergent_half_angle_deg=30.0,
                          divergent_half_angle_deg=15.0)
result = solve_nozzle_kinetics(reactants=reactants, chamber_pressure_pa=6.89476e6,
                               nozzle=nozzle)
print(result.vacuum_specific_impulse_s, dict(result.exit.mole_fractions)["H"])
```

## Kinetics

Each reaction has a modified Arrhenius forward coefficient
`k_f = A T^b exp(-E/(R T))` (NASA RP-1328 eq. (2.4); MFiX eq. (4.8)). The
reverse coefficient is `k_r = k_f / K_c` with
`K_c = (p0/(R T))^dn exp(-dG0/(R T))` (RP-1328 eqs. (2.6)-(2.10); MFiX eqs.
(4.9)-(4.13)). `dG0` comes from the repository's NASA Glenn polynomials and
`p0 = 1 bar`, the standard state of those polynomials (RP-1328 writes the
1 atm form for its own data). Forward and reverse rates are therefore exactly
consistent with the repository's equilibrium solver. A test confirms this: at
three equilibrium states from `LocalEquilibriumProvider` the forward and
reverse rates of all 27 reactions agree to better than 1e-9.

Rates follow mass action. Third bodies `[M] = sum eps_l C_l` multiply both
directions (RP-1328 eqs. (2.29)-(2.32); MFiX eqs. (4.14)-(4.17)), and
collider-specific companion reactions follow Cantera's "Reaction Rates" page.
Falloff uses `k = k_inf (P_r/(1+P_r)) F` with the Troe `F` and base-10
logarithms (MFiX eqs. (4.18)-(4.20); Cantera "The Troe Falloff Function").
The net production rate is `omega_i = sum_j nu_ij (R_j - R_-j)` (RP-1328
eq. (2.33)).

`production_rate_jacobian` returns the analytic `d omega/d C`, including the
third-body and falloff terms. It agrees with central differences to better
than 1e-6 (row-scaled). Unit conversion from CHEMKIN defaults (cm, mol, s,
cal/mol, as stated in MFiX section "Units in the Mechanism") uses
1 cal_th = 4.184 J exactly (NIST SP 811 Appendix B.8).

## Mechanism

`gri30_h2o2_mechanism()` contains every reaction of GRI-Mech 3.0
`grimech30.dat` (SHA-256 `a7726208...cd763`) whose species all lie in
{H2, H, O, O2, OH, H2O, HO2, H2O2}: 27 reactions, including the duplicate
pairs and the Troe reaction `2OH(+M)<=>H2O2(+M)`. Each reaction records its
file line. A script parsed the downloaded file and compared every transcribed
value with it (A, b, E, efficiencies and DUPLICATE flags): 0 mismatches.

The following were omitted or adapted:

- Reverse rates use the NASA Glenn polynomials, not GRI's `thermo30.dat`.
- Efficiencies of CH4, CO, CO2, C2H6, N2 and AR are dropped because those
  species are absent. For the same reason the collider-specific lines 62
  (H+O2+N2) and 63 (H+O2+AR) are dropped.
- GRI-Mech 3.0 states that it is "designed to model natural gas combustion".
  Using it for rocket-nozzle recombination at 3500 K and 70 bar extrapolates
  beyond its optimization targets. Every nozzle result carries a warning that
  says so.
- No CH4/CO subset is provided.

## Stiff integrator

`solve_stiff` is the RODAS embedded Rosenbrock method of order 4(3). It
reproduces Hairer and Wanner's `rodas.f`: the `ROCOE` coefficients for
`METH=1`, the six-stage algebra of `ROSCOR`, the `DECOMR`/`SLVROD` case-1
linear algebra from `dc_decsol.f`, the error norm, step factors 1/6 to 5 with
safety 0.9, and the Gustafsson predictive controller. Finite-difference
Jacobian and `df/dt` increments are those of `rodas.f`. There is one
documented deviation: singular-stage failures count consecutively, not over
the whole call.

Verification against Hairer's stiff test set (`res_exact_pic` reference
solutions):

- ROBER (Robertson) at t = 1 ... 1e11, rtol 1e-8 with an analytic Jacobian:
  the largest relative error is 2.3e-8 (declared tolerance 1e-6). The error
  falls monotonically for rtol = 1e-4, 1e-6 and 1e-8.
- HIRES at t = 321.8122 and 421.8122, rtol 1e-8 with a finite-difference
  Jacobian: within the declared relative tolerance of 1e-5.

## Nozzle model

The model is the steady, inviscid, adiabatic assigned-area flow of RP-1328
section 2.4: continuity eq. (2.58), momentum eq. (2.60), energy eq. (2.61)
with no heat loss, and species `d sigma/dx = omega/(rho V)` (eqs. (2.33),
(2.57)).

Eliminating `V dV` gives eq. (2.65). With the mixture entropy of RP-1311
eqs. (2.16)-(2.17), this becomes `T ds = -sum mu_i d sigma_i`. The solver
therefore marches `(sigma_i, s)`. At each `x`, the temperature solves
`rho V = mdot/A(x)` on the current fixed-composition isentrope, with
`V = sqrt(2 (h0 - h))`. This is algebraically equivalent to RP-1328 eqs.
(2.80)-(2.85) and has the same singularity at frozen Mach 1 (a saddle point).

The Jacobian of this right-hand side is analytic, by implicit differentiation
through the closure. It agrees with central differences to better than 1e-4
(row-scaled). An earlier forward-difference version differed by up to 4e-3,
which is why the analytic form is used.

The nozzle geometry is a hyperbolic wall, `A/A_t = 1 + (x tan(theta)/r_t)^2`,
with separate convergent and divergent half angles. The throat radius sets the
residence time. Vacuum specific impulse is `(V_e + p_e A_e/mdot)/g0`
(RP-1311 eqs. (6.8)-(6.9)).

There are two start modes:

- `start="equilibrium_throat"` (default): shifting equilibrium from the
  chamber to `start_area_ratio` (default 1.1, where frozen Mach > 1), then
  finite rates. The equilibrium path uses the repository's equilibrium
  provider over the mechanism species, with the throat at maximum `rho V`.
  This mode works for realistic rates.
- `start="chamber"`: finite rates from the convergent section with the
  chamber composition. The choking mass flow is found by shooting (fixed point
  on `min A G_max`, then bisection), and the solver switches to the supersonic
  branch where the bracketing march choked. Zero rates reproduce
  frozen-at-chamber expansion exactly. With realistic (fast) throat chemistry
  the supersonic continuation leaves the saddle separatrix and re-chokes. This
  mode then fails closed with `nozzle_transonic_unresolved` after about 70 s.
  It succeeded for rate multipliers of 1e-3 (Isp 437.89 s) and 1e-2 (444.28 s)
  at 16-57 s per run.

### Bray sudden-freezing estimate

`bray_sudden_freezing` applies Bray's criterion `-(d alpha/dt)_e = K r_De`
(thesis eq. (36); K of order one, experimental range 0.42-1.25) along the
shifting-equilibrium path. It then expands with frozen composition after the
freezing point. For the multi-reaction gas, `alpha` is taken as total mol/kg,
which only `dn != 0` reactions change, and `r_De` as the one-way
mole-increasing rate. This generalization is a choice of this model; Bray's
thesis treats a single dissociation reaction.

## Results for the test case

Case: LOX/LH2 at O/F 6.0, 1000 psia (6.89476 MPa), A_e/A_t = 40, hyperbolic
nozzle with r_t = 5 cm and 30/15 degree half angles. Source: the run logs of
this work and `tests/test_finite_rate_kinetics.py`.

Bounds from independent NASA CEA 3.3.4 runs, with products restricted to the
eight species (`scripts/generate_finite_rate_cea_bounds.py`):

- frozen at the chamber: 430.4426 s
- frozen at the throat: 436.6738 s
- shifting equilibrium: 451.6687 s

The repository's own solver gives frozen 430.4413 s and shifting 451.6670 s.

Finite-rate results (default mode, rtol 1e-6):

- Vacuum Isp is 450.642 s, a kinetic loss of 1.03 s (0.23 %) against
  shifting equilibrium.
- Atomic hydrogen freezes near a mole fraction of 0.006 (it starts at 0.024).
  OH falls from 0.027 to 5e-4 and O from 1.5e-3 to 1.1e-5. The exit
  temperature is 1412.5 K, against 1440.8 K in equilibrium.
- Elements are conserved to 5e-10. This matches the equilibrium start
  composition's own residual; in chamber mode the residual is 1e-14.
- With r_t = 1 cm the vacuum Isp is 449.385 s, a loss of 2.28 s: smaller
  nozzles lose more, as expected from the shorter residence time.

Limits:

- Zero rates in chamber mode give 430.441325 s, against the repository's
  frozen 430.441314 s (2.6e-8 relative) and CEA's 430.4426 s. c* is
  2263.7533 m/s against CEA's 2263.7589 m/s.
- Rate multipliers 1, 1e2 and 1e4 leave gaps to shifting equilibrium of
  1.03, 0.060 and 0.0022 s. The gap decreases monotonically.
- Zero rates in the default mode give 440.42 s (frozen from A/A_t = 1.1),
  which lies between CEA's frozen-at-throat and shifting values.

Convergence:

- Changing rtol: 1e-5 gives 450.642021 s, 1e-6 gives 450.642034 s and 1e-7
  gives 450.642034 s.
- Changing the start station: A/A_t = 1.05 gives 450.641836 s, 1.1 gives
  450.642034 s and 1.2 gives 450.642521 s.

Bray estimate: 449.62, 450.34 and 450.48 s for K = 0.42, 1.0 and 1.25, with
freezing at A/A_t of 3.05, 3.76 and 3.97 (T_f of 2577, 2474 and 2448 K). The
sudden-freezing estimate falls 0.16-1.02 s below the finite-rate result. In
the finite-rate solution the bimolecular shuffle reactions keep converting OH
after three-body recombination has frozen, which a sudden freeze cannot
represent.

## What is not verified

- **No published ODK/TDK kinetic efficiency** for a fully specified case was
  read, so none is compared. NASA CR-393 (Lordi and Mates) was opened, but it
  studies catalysed hydrogen recombination for nuclear and electric rockets,
  not a LOX/LH2 engine case.
- The rate set's accuracy at rocket conditions is not validated against
  experiment (see the GRI-Mech caveat above).
- In chamber mode the transonic passage is not resolved for fast throat
  chemistry. In the default mode the assumption of equilibrium up to the
  start station is a model choice; its effect was measured as small for this
  case (2e-4 s for start ratios 1.05 to 1.1).
- No viscous, divergence, heat-loss, condensed-phase or two-dimensional
  effects are included.
- Runtime: the default mode takes 0.5-4 s per case and chamber mode tens of
  seconds.

## References

- NASA-RP-1328: K. Radhakrishnan, LSENS Part I, NASA RP-1328 (1994), chapter 2. <https://ntrs.nasa.gov/citations/19950008642>
- NASA RP-1311 Part I (Gordon and McBride, 1994), eqs. (2.14), (2.16)-(2.17), (6.8)-(6.9). <https://ntrs.nasa.gov/citations/19950013764>
- MFIX-25.2-CHEMKIN: NETL MFiX 25.2 documentation, section 4.14.4. <https://mfix.netl.doe.gov/doc/mfix/25.2/html/model_setup/reactions/chemkin.html>
- CANTERA-3.2-RATES: Cantera 3.2.0 "Rate Constant Parameterizations" and "Reaction Rates". <https://cantera.org/stable/reference/kinetics/rate-constants.html>, <https://cantera.org/stable/reference/kinetics/reaction-rates.html>
- GRI-MECH-3.0: Smith et al., GRI-Mech 3.0 documentation and grimech30.dat. <http://combustion.berkeley.edu/gri-mech/version30/text30.html>, <http://combustion.berkeley.edu/gri-mech/version30/files30/grimech30.dat>
- HAIRER-RODAS: Hairer and Wanner, rodas.f and dc_decsol.f. <https://www.unige.ch/~hairer/prog/stiff/rodas.f>, <https://www.unige.ch/~hairer/prog/stiff/dc_decsol.f>
- HAIRER-TESTSET: Hairer, testset of stiff ODEs (ROBER, HIRES). <https://www.unige.ch/~hairer/testset/testset.html>
- BRAY-1960-THESIS: K. N. C. Bray, PhD thesis, University of Southampton, 1960, section 5. <https://eprints.soton.ac.uk/458375/1/38084.pdf>
- NIST-SP811-B8: NIST Guide to the SI, Appendix B.8. <https://www.nist.gov/pml/special-publication-811/nist-guide-si-appendix-b-conversion-factors/nist-guide-si-appendix-b8>
- NASA-CEA-PY: NASA Glenn, cea 3.3.4 Python package (bounds generator, RP-1311 example 12 for `n_frz`). <https://pypi.org/project/cea/>
