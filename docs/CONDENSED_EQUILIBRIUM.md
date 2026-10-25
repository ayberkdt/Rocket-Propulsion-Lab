# Condensed-phase products in chemical equilibrium

`rocket_propulsion.thermochemistry.condensed_equilibrium` adds pure condensed
products (graphite, ice, liquid water, aluminium, alumina, aluminium chloride,
aluminium nitride, ...) to the local Gibbs-energy equilibrium solver and a
two-phase infinite-area-combustor (IAC) rocket calculation. It closes the
"Condensed-phase products (Al₂O₃, carbon)" gap without NASA CEA installed.
The existing gas-only `LocalEquilibriumProvider` is unchanged.

Status: experimental. Verified against NASA CEA 3.3.4 at 7 equilibrium points
and 4 rocket problems (17 stations), see [Verification](#verification-against-nasa-cea-334).

## Data

`scripts/sync_nasa_condensed.py` reads NASA's `thermo.inp` at commit
`2f79a647e85737742fb08a127da00d1ea1d069da` (SHA-256 checked,
`fa7746572952d74e249e818a82a35c113829742fb421a308e167185528884363`) with the
fixed-column reader of `scripts/sync_nasa_thermo.py` and writes
`src/rocket_propulsion/thermochemistry/data/nasa_cea_condensed.json`
[NASA-CEA-THERMO-CONDENSED]:

- 186 non-ionic gas products made only of Al, C, Cl, H, N, O that are not
  already in `nasa_glenn_2002.json` or `nasa_cea_thermo.json` (Al, ALO, AL2O,
  ALCL, ALCL2, ALCL3, ALOH, HALO2, HCL, CL, CL2, ..., and the C/H/O/N gases
  CEA also considers, such as C2H6);
- 16 condensed products: `AL(cr)`, `AL(L)`, `AL2O3(a)`, `AL2O3(L)`,
  `ALCL3(cr)`, `ALCL3(L)`, `ALH3(a)`, `ALN(cr)`, `ALN(L)`, `AL(OH)3(a)`,
  `AL4C3(cr)`, `C(gr)`, `H2O(cr)`, `H2O(L)`, `NH4CL(II)`, `NH4CL(III)`;
- 2 reactants with their 298.15 K enthalpy: `AL(cr)` (0 J/mol) and
  `NH4CLO4(I)` (-295767 J/mol). No binder (HTPB or similar) is bundled.

With the base catalog this is exactly the product set NASA CEA 3.3.4 uses
for these elements (`products_from_reactants=True`).

The script refuses to write unless these identities hold, and stores the
measured residuals in the file's `verification` block:

- heat of formation: `R T (H/RT)(298.15 K) = Delta_f H(298.15)`
  [NASA-RP1311-P1-CONDENSED, eqs. (4.1)-(4.2)] with R = 8314.51 J/(kg-mole K)
  [NASA-RP1311-P1-CONDENSED, section 2.1] for the 196 records whose first
  interval covers 298.15 K or starts at 300 K; worst residual 0.088 J/mol
  (`H2O(L)`), tolerance 0.5 J/mol. The reactant records reproduce their
  header values to 4.1e-5 J/mol (`AL(cr)`) and 2.8e-4 J/mol (`NH4CLO4(I)`);
- phase transitions: the Gibbs energies of consecutive phases agree at the
  shared temperature: alumina at 2327 K to 3.1e-8 RT (enthalpy of fusion
  111400 J/mol), aluminium at 933.61 K to 3.2e-8 RT, ice/water at 273.15 K to
  2.7e-4 RT (largest; a 0.03 K shift of the melting point);
- continuity of H/RT and S/R across interval boundaries: at most 1.24e-3
  (`ALOCL`, 1000 K). The source record `ALN(L)` carries a step of 3.03 in both
  H/RT and S/R at 2700 K; it is recorded, not altered.

## Method

The provider follows NASA RP-1311 Part I [NASA-RP1311-P1-CONDENSED]:

- Ideal gases plus pure condensed phases; the condensed phases occupy no
  volume, so `p/rho = n R T` with `n` the gas moles per unit mass and
  `M = 1/n` (section 2.2, eqs. (2.1)-(2.3)); `MW` over all species is eq.
  (2.4a).
- Chemical potentials (2.11): `mu_j/RT = mu°_j/RT + ln(n_j/n) + ln P` for gases,
  `mu°_j/RT` for condensed species; equilibrium condition (2.10)
  `mu_j/RT = sum_i pi_i a_ij` with `pi_i = -lambda_i/RT`.
- Reduced Newton-Raphson iteration (section 2.3.2): eqs. (2.24) element rows,
  (2.25) one row per included condensed species, (2.26) total gas moles,
  (2.27) for assigned enthalpy (HP) and (2.28) for assigned entropy (SP);
  gas corrections from (2.18). Unknowns `pi_i`, `Delta n_c`, `Delta ln n`
  and `Delta ln T`.
- Initial estimates of section 3.1 (`n = 0.1` kg-mol/kg spread over the
  gases, no condensed species, 3800 K for HP and SP) or a warm start from a
  neighbouring state; control factor eqs. (3.1)-(3.3) with SIZE = 18.420681;
  updates (3.4); convergence tests (3.5)-(3.6c) with tighter bounds (1e-10
  instead of 0.5e-5, 1e-6, 1e-4 and 0.5e-4).
- Singular Newton matrices (section 3.6, e.g. eq. (3.11): every significant
  gas in one element ratio, as for stoichiometric H2/O2 at 400 K) are not
  handled by CEA's reset and component switch but by a column-scaled
  Tikhonov least-squares step (weight 1e-12) that keeps the undetermined
  direction of the element potentials at its previous value. The major
  composition is unaffected; the split of the trace species in that
  direction is not resolved (for stoichiometric H2/O2 at 400 K the solver
  reports O2 at about 1e-11 instead of H2 and O2 in a 2:1 ratio near 1e-27).
- States without a gas phase (gas below 1e-8 of all moles, e.g.
  stoichiometric H2/O2 below its boiling point at 1 bar, Al-rich Al/O2 at
  1000 K) raise `DomainError` with code `no_gas_phase`: the ideal-gas
  formulation cannot carry the pressure. NASA CEA 3.3.4 reports
  non-convergence for the same two states.
- Condensed species (section 3.4): after each convergence, the candidate with
  the most negative `mu°_c/RT - sum_i pi_i a_ic` (eq. (3.7)) whose range
  covers T is included, at most one per convergence; a species whose amount
  becomes negative or zero is removed. Selection uses the per-mole value of
  eq. (3.7); RP-1311 says "largest negative change to Gibbs energy" without a
  normalization, so the path (not the converged state) may differ from CEA.
- Phase transitions (section 3.5): with T assigned only phases whose range
  covers T are used. With T variable, an included phase found outside its
  range by more than 50 K is switched for the neighbouring phase, within 50 K
  both phases are included, which pins T to the transition (for alumina the
  data meet within 1e-5 K of 2327 K). A phase with no neighbour is removed.
- Gas polynomials are evaluated with the nearest interval outside a species'
  range; a converged state in which such a species exceeds a mole fraction of
  1e-6 raises `DomainError`, smaller ones give the info warning
  `trace_species_extrapolated`.

Two-phase IAC rocket (chapter 6): chamber HP equilibrium; every nozzle state
is an SP equilibrium at the chamber entropy (no velocity or thermal lag,
section 6.1), warm-started from the nearest solved pressure; `u = sqrt(2(h_inf
- h))` (6.5); throat at maximum `rho u` (golden section on `ln p` between
`p_c/2.5` and `p_c/1.2`, equivalent to the sonic condition of section 6.3.3);
exits from continuity (6.1); vacuum impulse from (6.7) with zero ambient
pressure; `c* = p_c/(rho u)_t`. The function reuses `RocketStation`,
`RocketExit`, `RocketPerformance` and the root helpers of
`rocket_performance.py`.

`solve_rocket_performance(provider=CondensedEquilibriumProvider())` works only
when every chamber species above a mole fraction of 1e-10 is a base-catalog
gas, because that function freezes the chamber composition as a gas
`Mixture`; for LOX/LH2 it agrees with `solve_condensed_rocket_performance`
to 1e-7, and with graphite in the chamber it raises `InputError` (tested).

## Public API

- `CondensedEquilibriumProvider` (`name = "local-gibbs-condensed"`):
  `solve(problem) -> EquilibriumResult` (implements `EquilibriumProvider`;
  accepts `EquilibriumProblem` TP/HP and `CondensedEquilibriumProblem`
  TP/HP/SP; UV raises `InputError`) and
  `solve_detailed(problem, initial=None) -> CondensedEquilibriumSolution`.
- `CondensedEquilibriumProblem.tp/hp/sp(...)`: reactant keys of
  `get_condensed_reactant` (`"AL(cr)"`, `"NH4CLO4(I)"` plus every reactant of
  the base catalog).
- `CondensedEquilibriumSolution`: `result`, `gas_moles`, `condensed`
  (`CondensedPhaseAmount`), `condensed_mass_fraction`,
  `mixture_molecular_mass_kg_mol` (MW), `element_potentials`,
  `inclusion_tests` (eq. (3.7) values), `events`, `newton_iterations`,
  `state` (`CondensedEquilibriumState`, for warm starts).
- `solve_condensed_rocket_performance(reactants=..., chamber_pressure_pa=...,
  area_ratios=..., reactant_temperature_k=298.15, provider=None) ->
  CondensedRocketPerformance` (`performance: RocketPerformance`,
  `stations: tuple[CondensedStationDetail, ...]`).
- `condensed_catalog()`, `get_condensed_reactant(key)`, `CondensedSpecies`.

In `EquilibriumResult`, `total_moles` counts gas and condensed moles and
mole fractions refer to that total (as CEA prints them);
`molecular_mass_kg_mol` is RP-1311's `M` (mass per gas mole) so that
`p/((R/M) T)` is the two-phase density; `cp_j_kg_k` is the frozen
`sum n_j cp_j` over all species and `cv = cp - n R`.

## Verification against NASA CEA 3.3.4

`scripts/generate_cea_condensed_points.py` ran NASA CEA 3.3.4 [NASA-CEA-PY-CONDENSED]
(`thermo.lib` SHA-256 `8e5df1cc...ec52a`) with its full product set and wrote
`tests/reference_data/cea_condensed_points.json` (also storing the SHA-256 of
`thermo.inp` and of the bundled JSON). `tests/test_condensed_equilibrium.py`
(21 tests, about 4 s) needs no CEA install.

Points:

- rocket, IAC, `p_c` = 68.9476 bar: LOX/CH4 O/F 0.8 (graphite in chamber,
  throat and exits, area ratios 10 and 40), LOX/CH4 O/F 1.2 (graphite only at
  the exit, area ratio 40), LOX/RP-1 O/F 1.0 (area ratio 40),
  `AL(cr)`/`NH4CLO4(I)` 20/80 by mass (liquid alumina; area ratios 10, 40,
  60, 80, 100: liquid, liquid, liquid + alpha at 2327 K, liquid + alpha at
  2327 K, alpha);
- equilibrium: H2/O2 3:1 by moles at 1 bar and 250 K (ice), 300 K and 350 K
  (liquid water); Al/O2 50/50 by mass at 1 bar, 2000 K (alpha alumina) and
  3000 K (liquid), and HP at 20 bar (4743.9 K, liquid); CH4(L)/O2(L) O/F 1.0
  at 10 bar and 1000 K (graphite).

Declared tolerances (after one exploratory comparison, not before): 3e-5
relative for temperature, pressure, M, MW, c* and vacuum Isp; 3e-5 absolute
for condensed mass fractions and gas mole fractions at or above 1e-4; 2e-5
relative for entropy; 1 J/kg + 1e-5 relative for enthalpy. The set of
condensed phases must match CEA exactly at every point and station.

Largest measured deviations (this commit, all points):

- rocket: c* 5.1e-6 (LOX/CH4 O/F 1.2); temperature 2.0e-6 (LOX/CH4 O/F 0.8
  throat); throat pressure 8.9e-6, exit pressure 1.4e-6; M 1.8e-6; MW 1.7e-6;
  vacuum Isp 4.2e-6 (LOX/CH4 O/F 0.8, area ratio 40); condensed mass fraction
  per phase 5.8e-6 absolute (two-phase alumina at area ratio 80); gas mole
  fraction 1.7e-6 absolute;
- equilibrium: HP temperature 1.9e-13 relative; M 2.6e-9; MW 1.1e-9;
  condensed mass fraction 2.9e-9 absolute; gas mole fraction 4.0e-9 absolute;
  enthalpy below 0.001 J/kg absolute; entropy 5.7e-6 relative, which is the
  ratio of the gas constants (8.31451 in CEA, 8.314462618 J/(mol K) here).

The rocket deviations are dominated by CEA's own convergence: its tolerances
of eqs. (3.5)-(3.6c) and its area-ratio iteration (CEA reports 60.00106 for
a requested 60; the tests pass CEA's reported area ratios to the local
solver).

Robustness sweep (not a CEA comparison, run once while developing): LOX/CH4
and LOX/RP-1 rockets at O/F 0.4-3.6 (20 bar, area ratios 5, 50, 200),
Al/AP 5-40 % Al (50 bar, area ratios 3, 30, 300), H2/O2 TP 200-600 K at four
ratios and Al/O2 TP 400-6000 K at three ratios all converged except the
gas-free states above and Al/O2 90/10 HP at 1, 10 and 100 bar
(`ConvergenceError`; CEA 3.3.4 does not converge there either).

Also tested: element conservation below 1e-9; Gibbs phase rule (condensed
phases at most the number of elements at assigned T and p); eq. (2.10) for
included condensed phases (1e-8) and for gases above a mole fraction of 1e-8
(1e-6); no excluded candidate with a negative eq. (3.7) value; the Gibbs
energy with condensed products below the gas-only constrained minimum;
removal of graphite carried in by a warm start (equal to the cold start to
1e-10); the TP switch from alpha to liquid alumina; the SP "add at transition"
path; rejection of a gas-free state; the eq. (3.11) singular case; input
validation.

## Not verified / limitations

- No hardware or test-stand data; agreement with CEA verifies the
  implementation of the same model with the same data, not the model.
- Zero velocity and thermal lag between particles and gas; real aluminized
  motors lose impulse to two-phase flow lag, which is not modelled.
- No ionized species (CEA can include them); no binder data (HTPB) bundled.
- FAC (finite-area combustor), frozen two-phase expansion, UV/SV problems and
  equilibrium derivatives (`gamma_s`, sound speed) with condensed phases are
  not implemented; the reported `gamma` is the frozen ratio.
- CEA's exact choice among several eq. (3.7) candidates and its extrapolation
  rules for gas polynomials below 300 K were not read from the CEA source and
  are not claimed to be reproduced.
- Reactants other than `AL(cr)` and `NH4CLO4(I)` with assigned enthalpy are
  those of the base catalog; reactant temperatures other than the assigned
  ones are not supported for condensed reactants.

## References

- **NASA-RP1311-P1-CONDENSED**: Gordon and McBride, NASA RP-1311 Part I (1994),
  sections 2.1-2.3.2, 2.6, 3.1-3.6, 4.1 and chapter 6.
  <https://ntrs.nasa.gov/citations/19950013764>
- **NASA-CEA-THERMO-CONDENSED**: NASA CEA `thermo.inp`, commit
  2f79a647e85737742fb08a127da00d1ea1d069da.
  <https://github.com/nasa/cea/blob/2f79a647e85737742fb08a127da00d1ea1d069da/data/thermo.inp>
- **NASA-CEA-PY-CONDENSED**: NASA CEA Python package 3.3.4.
  <https://pypi.org/project/cea/>
- **NASA-CEA-REPO-CONDENSED**: NASA CEA repository. <https://github.com/nasa/cea>

Full entries: [REFERENCES.md](REFERENCES.md).
