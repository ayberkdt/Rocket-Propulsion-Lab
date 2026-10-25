# Combustion efficiency, product composition and chamber processes

Two modules cover the chamber, from the injector face to the throat:

- `rocket_propulsion.thermochemistry.finite_area_combustor` solves the NASA CEA
  finite-area combustor (FAC).
- `rocket_propulsion.thermochemistry.combustion_efficiency` covers stoichiometry,
  product composition, combustion-efficiency models, the JANNAF Isp loss chain
  and O/F optimisation.

Both use the local Gibbs solver and the shifting-equilibrium rocket solve. That
solve is checked against NASA CEA 3.3.4 in `tests/reference_data/cea_gate_points.json`.
Only the Python standard library is used.

## Finite-area combustor

NASA RP-1311 Part I, chapter 6, defines the FAC model [RP1311-FAC]. The chamber
has a constant area `A_c`, and propellants enter it with zero velocity. Heat
added in this duct raises the entropy and lowers the pressure, a Rayleigh-type
loss. The model has four stations:

1. `inj`: assigned-enthalpy equilibrium at the injector pressure.
2. `inf`: a fictitious infinite-area combustor. It has the same enthalpy, a
   lower pressure `p_inf`, and the entropy `s_inf` that the nozzle uses.
3. `c`: the combustor end. This is the subsonic point with `A/A_t = A_c/A_t`
   on the `s_inf` isentrope.
4. `t`: the throat on the same isentrope.

The solver iterates `p_inf` until the momentum balance (6.28) holds:
`p_inj = (p + rho u^2)_c`. It starts from the empirical estimate (6.30) and
updates with (6.32). From the second pass onwards it takes secant steps on
`ln p_inf`, which lead to the same fixed point. With option 2 the caller gives
`mdot/A_c` instead of `A_c/A_t`; the solver then sets `A_c/A_t` on every pass
from (6.33). The results use these definitions:

- `c* = p_inf A_t / mdot`, the CEA convention (6.11).
- `injector_characteristic_velocity_m_s = p_inj A_t / mdot`, reported
  separately and labelled as not a CEA output.
- Nozzle stations: the infinite-area results at `p_inf`.

NASA's CEA source implements the same loop in `RocketSolver_solve_fac`
[CEA-SRC-FAC].

```python
from rocket_propulsion.thermochemistry.finite_area_combustor import solve_finite_area_combustor
from rocket_propulsion.thermochemistry.rocket_equilibrium import PROPELLANT_PAIRS, reactants_for_mass_ratio

fac = solve_finite_area_combustor(
    reactants=reactants_for_mass_ratio(PROPELLANT_PAIRS["LOX_LH2"], 5.55157),
    injector_pressure_pa=5.33172e6, contraction_ratio=1.58, area_ratios=(25.0, 50.0, 75.0))
fac.stagnation_pressure_ratio, fac.combustor_end_mach_number, fac.characteristic_velocity_m_s
```

### Agreement with NASA CEA

The script `scripts/generate_cea_fac_points.py` records 6 FAC runs of CEA 3.3.4
in `tests/reference_data/cea_fac_points.json`:

- RP-1311 examples 9 and 10, LOX/LH2.
- LOX/LH2 with `A_c/A_t` = 1.3 and 3.
- LOX/CH4 and LOX/RP-1 with `A_c/A_t` = 2, at 68.95 bar.

The file also holds the NASA CEA2 regression outputs `example9.out` and
`example10.out` [CEA2-OUT]. Measured largest relative deviations of the local
solver from CEA 3.3.4, against the tolerances declared before comparison:

- `p_inf/p_inj`: 1.5e-5 (declared 2e-4).
- Combustor-end pressure: 4.0e-5 (5e-4).
- Combustor-end temperature: 3.3e-6 (1e-3).
- Combustor-end Mach number: 1.3e-4 (2e-3).
- Throat pressure: 1.7e-5 (5e-4).
- `c*`: 3.4e-6 (1e-3).
- Vacuum `I`: 3.9e-6 (1e-3).
- `A_c/A_t` with the mass-flux option: 1.1e-4 (2e-4). CEA itself stops at a
  contraction ratio of 1.58015 for an assigned 1.58.

The local `p_inj/p_inf` for example 9 is 1.084767, which equals the NASA CEA2
output value.

## Stoichiometry and mixture descriptors

The descriptors follow RP-1311 Part I, equations (9.13)-(9.19) [RP1311-EQ]:

- `V+` and `V-` are oxidation-state totals per kilogram. The valences are
  those of CEA's `atomic_data.f90` (H +1, C +4, N 0, O -2) [CEA-SRC-VALENCE].
- The chemical equivalence ratio is `r = -V+/V-`.
- The stoichiometric O/F is the O/F at `r = 1`.
- `phi = (o/f)_st/(o/f)`.
- The total-reactant density is `rho_0`, from equation (9.13).
- Density impulse is `I_s d`, from SP-125 equation (1-45) [SP125-DI]. SP-125
  equation (1-46) prints `(r_w - 1)` where RP-1311 (9.13) has `(o/f + 1)`. The
  RP-1311 form is used.

Stoichiometric O/F agrees with CEA's `chem_eq_ratio_to_of_ratio` to better than
1e-6 for LOX/LH2 (7.93668), LOX/CH4 (3.98926) and LOX/RP-1 (3.40567).

## Product composition and dissociation

`combustion_product_report` runs the shifting-equilibrium infinite-area solve.
It then re-solves each station (chamber, throat, exit) at its temperature and
pressure, so that every species above 1e-14 mol is kept. The report contains:

- Mole and mass fractions.
- A radical inventory. The list `DEFAULT_RADICALS` is a reporting choice.
- The element residual and the mass-fraction sum.
- An enthalpy split based on the assigned-enthalpy definition (4.2)
  [RP1311-H]: `h = h_f + h_s`, where `h_f = sum_j (n_j/m) Delta_f H_j(298.15)`.

The split gives these quantities:

- The available energy `h_c - h_ref`. The reference state is gas-phase
  equilibrium at 298.15 K and chamber pressure.
- The dissociation energy `h_f,c - h_f,ref` and its fraction of the available
  energy. The dissociation and sensible fractions sum to one by construction.
- The recombination heat release `h_f,c - h_f,station` at the throat and the
  exit.

For carbon fuels the 298 K reference contains CH4 and CO2 instead of CO. The
"dissociation energy" therefore includes the water-gas and methanation energy,
not only the energy held in radicals. Species whose NASA records start at 300 K
are absent at 298.15 K, as the solver prescribes.

The CEA 3.3.4 references cover 7 cases: LOX/LH2 at O/F 4, 6 and 8 and at 10,
68.95 and 200 bar, plus LOX/CH4 3.4 and LOX/RP-1 2.7. Mole fractions down to
1e-10 are compared at the chamber, the throat and an `A_e/A_t = 40` exit.
Measured largest relative deviations, against the declared tolerances:

- Mole fraction `x >= 1e-2`: 4.6e-6 (declared 2e-3).
- `1e-6 <= x < 1e-2`: 1.4e-4 (1e-2), for OH at the exit at 10 bar.
- `1e-10 <= x < 1e-6`: 2.5e-4 (5e-2), for O at 7.9e-10.
- Temperature: 7.2e-6 (1e-3).

CEA species missing from the bundled catalog appear only below 7e-10: C2O,
CH2OH and CH3O.

## Combustion-efficiency models

- Measured efficiencies come from NASA/TP-1999-208522, equations (6)-(11)
  [TP-208522]: `c* = p A_t/mdot`, `C_F = F/(p A_t)`, `I = F/(mdot g0)` and their
  ratios to ODE. The relation `eta_I = eta_c* eta_C_F` holds whenever
  `I g0 = c* C_F`.
- Vaporization-limited `c*` follows NACA TN 4219, equation (5) [TN4219]:
  `eta_c* = [c*(O w_o/F w_f)/c*(w_o/w_f)] (O w_o + F w_f)/(w_o + w_f)`.
  Equations (6) and (7) are the cases where one propellant is fully vaporized.
- Mixing and energy-release efficiency follow the two-streamtube Rupe model
  and the JANNAF simplified ERE bookkeeping implemented in ROCCID, NASA
  CR-187109 equations (2.4)-(2.11) [ROCCID].
  - The streamtubes are `MR_R = MR_M eta_m` and `MR_L = MR_M/eta_m`.
  - The mixing parameter is `eta_m = 0.449 - 0.399 E_m + 0.944 E_m^2`.
  - The efficiencies are `eta_c* = C*_Del/C*(MR)` and
    `eta_ERE = I_MZ/I(MR)`.

  `E_m` is evaluated as a fraction. This is an interpretation: ROCCID quotes
  `E_m` in percent in equation (2.6). At `E_m = 1` the fit gives
  `eta_m = 0.994`.
- The SP-8089 maximum-mixing orifice ratio uses equation (1) with the Table IV
  mixing factors [SP8089].
- Characteristic length is `L* = V_c/A_t`, SP-125 equation (4-4). The stay time
  is `t_s = V_c/(mdot v)`, equation (4-3). Results are checked against the
  Table 4-1 ranges [SP125-LSTAR].
- The chamber momentum pressure loss uses TP-1999-208522 equation (4). For the
  FAC solution with zero injection velocity it reproduces `p_inf/p_inj` to
  1e-6, which the tests check.

## Isp loss chain

`jannaf_isp_loss_chain` applies the losses in the order of TP-1999-208522,
printed p. 25 and Figures 20 and 24 [TP-208522]. The boundary-layer thrust
decrement is converted to an impulse decrement and subtracted from the Isp
[TM-89814]. The chain is:

1. ODE.
2. ODK: multiply by the kinetic efficiency.
3. Two-dimensional inviscid: multiply by the divergence and shock efficiency.
4. Boundary layer: subtract the impulse decrement.
5. Delivered: multiply by `eta_c*` for the energy-release loss.

The factors after kinetics go to `propulsion.losses.calculate_loss_budget`,
with `I_ODK` as the ideal value and the discharge factor disabled. That keeps
the existing loss-budget bookkeeping in one place. The tests check that the
delivered Isp from both paths agrees to 1e-12.

The kinetic efficiency is an input. `frozen_equilibrium_limits` gives the two
chemistry limits of RP-1311 sections 6.1 and 6.5: instantaneous and infinitely
slow reactions.

## O/F optimisation

`optimize_mixture_ratio` runs a golden-section search over a bracket. The
objective is either the vacuum `I` or the density impulse `rho_0 I_vac`, and
each evaluation is a CEA-equivalent solve. No propellant densities are bundled;
the caller supplies them.

The reference optima come from golden-section searches on CEA 3.3.4 itself, at
68.95 bar and `A_e/A_t = 40`:

- Maximum Isp: LOX/LH2 at O/F 4.727 (455.08 s), LOX/CH4 at 3.413 and LOX/RP-1
  at 2.755.
- Maximum density impulse, using declared test densities that are not property
  data: LOX/CH4 at 4.039, LOX/RP-1 at 2.918, and LOX/LH2 at the bracket edge.

The local solver reproduces each optimum O/F within 0.03 and each objective
within 2e-4. The fast suite runs only the LOX/LH2 Isp case; set
`RPL_SLOW_TESTS=1` for the rest.

## Not verified or not implemented

- The Priem-Heidmann generalized-length vaporization correlation of NASA
  TR R-67 is not implemented. The report could not be found on NTRS, and ROCCID
  only paraphrases it.
- ROCCID equations (3.17)-(3.19) are not used. As printed, they scale drop size
  and injection velocity with a sign that contradicts the stated physics.
- The original definition of the Rupe mixing efficiency `E_m` (JPL reports)
  was not opened. Only the ROCCID use of `E_m` is implemented.
- The JANNAF manual itself (CPIA 246) was not available. The loss chain follows
  its use in NASA TP-1999-208522 and TM-89814.
- No finite-rate kinetics (ODK) is computed here. Other modules provide it.

## References

- [RP1311-FAC] Gordon and McBride, NASA RP-1311 Part I (1994), chapter 6, eqs. (6.1)-(6.33). <https://ntrs.nasa.gov/citations/19950013764>
- [RP1311-EQ] Same report, chapter 9, eqs. (9.10)-(9.19), printed pp. 42-43. <https://ntrs.nasa.gov/api/citations/19950013764/downloads/19950013764.pdf>
- [RP1311-H] Same report, section 4.1, eqs. (4.1)-(4.2), printed p. 19. <https://ntrs.nasa.gov/citations/19950013764>
- [CEA-SRC-FAC] NASA CEA source `rocket.f90`, `RocketSolver_solve_fac`, commit 2f79a647. <https://github.com/nasa/cea/blob/2f79a647e85737742fb08a127da00d1ea1d069da/source/rocket.f90>
- [CEA-SRC-VALENCE] NASA CEA source `atomic_data.f90` and `mixture.f90`, commit 2f79a647. <https://github.com/nasa/cea/blob/2f79a647e85737742fb08a127da00d1ea1d069da/source/mixture.f90>
- [CEA2-OUT] NASA CEA2 regression output, RP-1311 examples 9 and 10. <https://github.com/nasa/cea/blob/2f79a647e85737742fb08a127da00d1ea1d069da/test/main_interface/reference_output/example9.out>
- [CEA-PY] NASA CEA Python package 3.3.4, used only to generate the reference files. <https://pypi.org/project/cea/>
- [TN4219] Heidmann and Priem, NACA TN 4219 (1958), eqs. (1)-(7). <https://ntrs.nasa.gov/citations/19930085136>
- [ROCCID] Muss, Nguyen and Johnson, NASA CR-187109 (1991), section 2.2, eqs. (2.4)-(2.11). <https://ntrs.nasa.gov/citations/19910014917>
- [SP8089] NASA SP-8089, Liquid Rocket Engine Injectors (1976), eq. (1) and Table IV. <https://ntrs.nasa.gov/citations/19760023196>
- [SP125-LSTAR] Huzel and Huang, NASA SP-125 (1967), eqs. (4-3)-(4-4) and Table 4-1. <https://ntrs.nasa.gov/citations/19710019929>
- [SP125-DI] Same book, eq. (1-45) and Table 1-7. <https://ntrs.nasa.gov/api/citations/19710019929/downloads/19710019929.pdf>
- [TP-208522] Jankovsky, Smith and Pavli, NASA/TP-1999-208522, eqs. (3)-(11) and Figs. 20 and 24. <https://ntrs.nasa.gov/citations/19990051032>
- [TM-89814] Davidian, NASA TM-89814 (1987). <https://ntrs.nasa.gov/citations/19870010946>

Full entries, with page locators and file hashes, are in
[REFERENCES.md](REFERENCES.md).
