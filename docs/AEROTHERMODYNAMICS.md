# Viscous flow, boundary layers and aerothermodynamics

Package: `rocket_propulsion.aerodynamics` — modules `boundary_layer`, `couette`,
`aerothermal` and `entry`. Python standard library only, SI units, frozen slotted
result records, typed `DomainError`/`InputError`/`ConvergenceError` on invalid
input (fail closed). Every equation and coefficient below was read from the cited
source in this work package; page locators are in the module docstrings and in
`tests/reference_data/aero_viscous.json`.

## Scope

**Included**

- Blasius laminar flat plate solved numerically (RK4 shooting on `f''(0)`), with
  skin friction, boundary-layer, displacement and momentum thicknesses
  (Anderson Eqs. 18.15-18.28).
- Compressible laminar and turbulent flat plate by the reference-temperature
  method: Eckert form `T*/Te = 1 + 0.032 Me^2 + 0.58 (Tw/Te - 1)` (Anderson
  Eq. 18.53) and the Meador-Smart laminar/turbulent forms printed by Anderson
  (Secs. 18.4.1, 19.2.2); heating by Reynolds analogy `CH = (cf/2) Pr^-2/3`
  (Anderson Eq. 18.50 = Tauber's Colburn form, Eq. 52).
- Turbulent incompressible flat plate: `delta = 0.37 x Re^-0.2`,
  `Cf = 0.074 Re^-0.2` (Anderson Eqs. 19.1-19.2), local `cf/2 = 0.0296 Re^-0.2`
  (Tauber Eq. 54) and Meador-Smart `0.02296/0.02667 Re^-0.139`.
- Van Driest (II) compressible turbulent skin friction exactly as written by
  Hopkins (NASA TN D-6945 Eqs. 1-18: Karman-Schoenherr incompressible law, Keyes or
  power-law viscosity, recovery factor 0.88 by default, `gamma = 1.4`).
- Recovery factors (laminar `Pr^1/2`, turbulent approx. `Pr^1/3`, Couette `Pr`) and
  adiabatic-wall temperature.
- Incompressible Couette flow (closed form, Anderson Eqs. 16.6-16.59) and
  compressible Couette flow for constant `Pr` and `cp` with any `mu(T)`: the energy
  integral (Eq. 16.80) gives `h(u)` exactly, the constant shear follows by Simpson
  quadrature of `D = (1/tau) int mu du` (Eq. 16.65), the heat flux by Eq. 16.95. This
  is an exact reduction of the book's nested shooting method under the same
  assumptions.
- Stagnation-point heating: Anderson's boundary-layer fits for a cylinder (0.57) and
  sphere (0.763) with the Newtonian velocity gradient (Eqs. 18.65, 18.70, 18.82);
  Sutton-Graves `qw = K (ps/R)^0.5 (hs - hw)` with `K` from their Table II or from the
  general mixture relation Eq. 42 and Table I (any mix of N2, O2, H2, He, Ne, Ar,
  CO2, NH3, CH4 — so CO2/N2/Ar planetary atmospheres are covered); Tauber's
  `1.83e-4 (rho/rn)^0.5 V^3 (1 - hw/Hs)` (Eq. 40 = Anderson Eq. 14.62).
- Tauber's air closed forms: swept cylinder (Eq. 41), laminar sharp cone (Eq. 46)
  and wedge/flat plate (Eq. 48), wing leading edge (Eq. 49), turbulent flat plate
  (Eqs. 56a/56b, valid 1.5-11 km/s, 3-60 km, ±15 %) and turbulent cone (`1.15 q_FP`,
  Eq. 59).
- Allen-Eggers ballistic entry closed form: velocity (Eq. 13), deceleration
  (Eq. 14), altitude/velocity/value of peak deceleration (Eqs. 15-17, sea-level case
  Eqs. 20-21), peak average heating (Eqs. 36-39, 41), peak stagnation heating for
  `qs = C (rho/sigma)^0.5 V^3` (Eqs. 45-48, 50), total heat input (Eqs. 33-34 =
  Anderson Eq. 14.60).
- Planar point-mass entry integrator (fixed-step RK4) of Chapman's polar equations
  of motion (NACA TN 4276 Eqs. 3-7, all terms kept) with constant `L/D`,
  inverse-square gravity (USSA 1976 Eq. 17), optional flat-earth mode (Allen-Eggers
  Eq. 1), exponential atmosphere, the U.S. Standard Atmosphere 1976 below 86 km, or a
  user density callable; reports peak deceleration, peak stagnation heating and the
  integrated heat load.
- U.S. Standard Atmosphere 1976 density below 86 km geometric altitude:
  `H = r0 Z/(r0 + Z)` (Eq. 18, `r0 = 6 356 766 m`) and `rho = P M0/(R* T_M)`
  (Eq. 42), reusing the repository's existing layer table and pressure function in
  `compressible/off_design_nozzle.py` rather than duplicating them.

**Not included — no verified source was opened, so nothing was implemented**

- Fay-Riddell in its original form (the 1958 J. Aeronaut. Sci. paper and the
  Avco report were not accessible); Detra-Kemp-Riddell; any radiative-heating
  correlation (Tauber-Sutton 1991 was not read). Radiative heating is therefore
  absent.
- A transition correlation. The transition Reynolds number is always an input;
  `ANDERSON_RULE_OF_THUMB_TRANSITION_REYNOLDS = 5e5` is exported only as the
  documented "rule of thumb" value (Anderson p. 935), never used as a default.
  Tauber's review shows measured cone transition Reynolds numbers spanning about
  1e6 to 3e7.
- The 1976 atmosphere above 86 km (species number-density model); the integrator
  fails closed if asked to start higher with `atmosphere="ussa1976"`.
- Real-gas/equilibrium-chemistry edge states, catalytic walls, ablation and mass
  addition, roughness, low-density (slip) effects, planet rotation and winds.

## Source notes and discrepancies found while reading

- Anderson Eq. 18.70 (sphere) is printed with `Pr^-0.65` and no exponent on
  `rho_e mu_e`, while the same page says the cylinder and sphere expressions differ
  only in the leading coefficient; the 1/2 power is also required dimensionally. The
  implementation uses `0.763 Pr^-0.6 (rho_e mu_e du_e/dx)^0.5 (haw - hw)`. Anderson
  attributes both fits to Van Driest (1956), not to Fay-Riddell.
- Anderson Example 18.3 prints `rho* = 0.500` in the `Re*` line (typo for 0.599);
  Example 19.3 prints `mu* = 2.05e-5` but uses 2.95e-5.
- Anderson Example 18.2 reads `T0/Te = 2.74` from a table for `Me = 2.94`; the
  isentropic formula gives 2.728, which shifts `T*` by 0.3 %.
- Anderson Example 16.2 gives `tau_w = 72 N/m^2` as data. The Sutherland-law
  quadrature in `compressible_couette` gives 83.1 N/m^2 for the same gap and plates,
  so the example is checked through `qw/tau_w = ue/2` (Eq. 16.95 with equal wall
  temperatures), which does not depend on the viscosity law.
- Hopkins Table 2 is reproduced only if the wall temperature is held at
  `(Tw/Taw) * Taw(r = 0.88)` when `r` is changed to 1.0 (documented in the JSON).
- Tauber Eq. 56 writes `0.9 - hw/H_c`; the nomenclature defines `H` as total
  enthalpy, so the function takes the inviscid-flow total enthalpy as input.

## Verification (measured in `tests/test_aero_viscous_aerothermal.py`)

- Blasius: `f''(0) = 0.332057` (book 0.332, +0.02 %); `eta_99 = 4.910` (book
  "approximately 5.0", -1.8 %); displacement coefficient 1.72079 (book 1.72,
  +0.05 %); momentum coefficient 0.66411 (book 0.664, +0.02 %); grid refinement
  2000→4000 steps changes `f''(0)` by less than 1e-9.
- Anderson Example 18.1a: drag 175.53 N (book 175.6, -0.04 %). Example 19.1a:
  `Cf = 2.770e-3`, drag 1351.9 N (book 1352). Example 18.2 (Eckert): 4981.3 N (book
  4976, +0.11 %), `T* = 610.9 K` (book 612.7, -0.29 %). Example 18.3
  (Meador-Smart): 5017.5 N (book 5008, +0.19 %). Example 19.2: 52 044 N (book
  51 890, +0.30 %). Example 19.3: 56 506 N (book 56 140, +0.65 %).
- Van Driest II vs Hopkins TN D-6945 Table 2 (8 rows × 4 ratio columns): worst
  deviation 1.6 % (row Me = 10, Tw/Taw = 1, Te = 222 K, power-law column at
  `Re_x = 1e5`: 0.919 vs 0.934); the next largest are 0.86 % and 0.77 % (Me = 10,
  Te = 222 K rows); mean absolute deviation over the 32 entries 0.30 %. For Anderson's Example 19.1b
  conditions (Me 2.94, Re 1.36e8, adiabatic) Van Driest II gives `Cf = 1.245e-3`
  against 1.34e-3 read by Anderson from Van Driest's 1951 chart (a different
  theory, Pr = 0.75): -7 %; this comparison is informational and not a test.
- Couette Example 16.1 (converted from English units): `u_mid = 100 ft/s`,
  `tau = 0.0897 lb/ft^2` (book 0.09), `Tmax = 519.59 °R` (519.6), `qw = 8.970
  ft·lb/(ft²·s)` (8.97), `Taw = 521.36 °R` (521.36). Example 16.2: `qw = 72 × ue/2 = 3.673e4 W/m^2`
  (book 3.68e4, -0.2 %), `Taw = 656.06 K` (book 656). The compressible solver reduces
  to `mu ue/D` at 1 m/s to 1e-6 and its quadrature is converged to 1e-9.
- Anderson Example 14.2 (Tauber Eq. 40): 45.784 W/cm² (book 45.78). For the same
  point, Sutton-Graves with Newtonian `ps = rho V^2` gives 43.7 W/cm²; Anderson
  quotes 45 W/cm² from Zoby's flight-data analysis.
- Sutton-Graves Eq. 42 vs their Table II computer correlation: air 0.1150 vs 0.1113
  (+3.3 %), CO2 0.1254 vs 0.1210 (+3.6 %), N2 0.1119 vs 0.1112 (+0.6 %). Table II's own
  average/maximum correlation errors for air are 3.3 %/9.8 %.
- U.S. Standard Atmosphere 1976 density vs Table I at 0, 38, 60, 70, 80 and 85.5 km:
  worst 1.6e-5 relative; geopotential altitudes within 1 m.
- Allen-Eggers: `e^-1/2 = 0.607`, `e^-1/3 = 0.717`, `e^-1/6 = 0.846` against the
  report's rounded 0.61, 0.72, 0.85. The integrator in flat-earth, no-gravity mode,
  exponential atmosphere (the report's own `rho0 = 0.0034 slug/ft^3`,
  `beta = 1/22000 ft^-1` fit), reproduces Eq. 13 to < 1e-6 along the path and the
  peak deceleration and peak stagnation heating to < 1e-3 for 90° and 30° entries.
- Measured effect of the Allen-Eggers assumptions (spherical earth, inverse-square
  gravity, exponential atmosphere with 7.2 km scale height, 7.5 km/s, 1000 kg/m²
  per unit `C_D A`): peak deceleration +3.9 % at 60°, +5.0 % at 20°, +25 % at 5°
  entry angle relative to the closed form. The closed form should not be used for
  shallow entries.

## Usage sketch

```python
from rocket_propulsion.aerodynamics.aerothermal import sutton_graves_heating
from rocket_propulsion.aerodynamics.entry import integrate_planar_entry

trajectory = integrate_planar_entry(
    entry_altitude_m=85_000.0, entry_velocity_m_s=7_000.0, entry_flight_path_angle_deg=-20.0,
    mass_kg=500.0, drag_coefficient=1.0, reference_area_m2=1.0,
    nose_radius_m=0.5, stagnation_heating_constant=1.83e-4, final_velocity_m_s=500.0,
)
peak = trajectory.peak_heating  # cold-wall Tauber stagnation heating along the path
```

## References

- **ANDERSON-FOA-2017** — J. D. Anderson Jr., *Fundamentals of Aerodynamics*, 6th ed.,
  McGraw-Hill Education, 2017, ISBN 978-1-259-12991-9 (bibliographic citation).
- **SUTTON-GRAVES-1971-TRR376** — NASA TR R-376, 1971:
  <https://ntrs.nasa.gov/citations/19720003329>
- **TAUBER-1989-TP2914** — NASA TP-2914, 1989:
  <https://ntrs.nasa.gov/citations/19890017745>
- **HOPKINS-1972-TND6945** — NASA TN D-6945, 1972:
  <https://ntrs.nasa.gov/citations/19730001588>
- **ALLEN-EGGERS-NACA-R1381** — NACA Report 1381:
  <https://ntrs.nasa.gov/citations/19930091020>
- **CHAPMAN-1958-TN4276** — NACA TN 4276, 1958:
  <https://ntrs.nasa.gov/citations/19930085059>
- **USSA-1976** — U.S. Standard Atmosphere, 1976:
  <https://ntrs.nasa.gov/citations/19770009539>
