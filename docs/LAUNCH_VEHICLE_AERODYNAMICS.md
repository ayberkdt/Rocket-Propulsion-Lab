# Launch-vehicle aerodynamics

Modules: `rocket_propulsion.aerodynamics.slender_body`, `rocket_propulsion.aerodynamics.launch_vehicle`,
`rocket_propulsion.aerodynamics.aero_database`. Tests: `tests/test_launch_vehicle_aero.py` with
reference values and page locators in `tests/reference_data/launch_vehicle_aero.json`.

Every equation and constant below was read from the cited source in this work package (2026-10-25).
The split follows what is public: (a) what is implemented and checked, (b) what industrial practice
documents but needs wind-tunnel or CFD data the user must supply, (c) what is not public.

## (a) Implemented and verified

### Body and fin normal force (`slender_body`)

- **Body geometry.** Pointed bodies made of a cone or tangent-ogive nose, cylinders and conical
  shoulders/boattails (Barrowman figure 3-7). Radius jumps are refused, as Barrowman requires.
  The tangent-ogive construction reproduces the NACA TR 1048 figure 6 bodies: generating radius
  25.000 in, nose 4.314 in, fineness ratios 21.1 and 13.1.
- **Slender-body (Munk) normal force.** `C_N_alpha = 2 [A(x1) - A(x0)] / A_r` per component and the
  centre of pressure `[x1 A1 - x0 A0 - V] / (A1 - A0)` (Barrowman equations 3-57 to 3-89). The tests
  check `C_N_alpha = 2` for a nose on its base area, a cone CP at 2/3 of its length, a negative
  boattail and the whole-body form `X = l - V/A(l)`.
- **Viscous cross flow (Allen and Perkins).** The TR 1048 cross-force distribution and equations 14
  and 15 give `C_N`, `C_L`, `Delta C_Dp`, `C_M` and the CP. For body I at M = 1.98 and alpha = 8 deg,
  with `c_d = 1.2 x 0.755` from the report's appendix, the code gives `C_L` = 0.699 (eq. 14) and
  0.710 (eq. 15). The report's calculated curve, read from figure 7(a) at 600 dpi, is 0.72 +- 0.03.
  At alpha = 7.8 deg the code is 22% (eq. 14) or 24% (eq. 15) above the measured point (0.55, read
  from the figure); the report's calculated curve also lies above that point. Integrating the line load numerically reproduces `C_N` and the CP.
- **Barrowman fins.** Subsonic `(C_N_alpha)_T` (equation 3-12) and the fin centre of pressure
  (equations 3-27, 3-30). Components combine with equations 3-106 and 3-107.
- **Fin-body interference.** `K_W(B)` and `K_B(W)` come from NACA TR 1307 equations 14 and 21. They
  match the report's figure 2 slender-body curve within reading (+-0.02) at r/s = 0.2, 0.4, 0.6 and
  0.8, and the stated limits 1 and 2, and 0 and 2. Barrowman's printed equation 3-96 omits the
  `+ pi/4` term; without it `K` would be negative at tau = 2, so it is not used. For r/s > 0.99999
  the function refuses the input, because cancellation makes the result inaccurate there.

### Drag, dynamic pressure and loads (`launch_vehicle`)

- **Dynamic pressure.** `q = gamma p M^2 / 2`, derived from the NASA Glenn relations
  `q = rho u^2 / 2`, `a = sqrt(gamma R T)` and `p = rho R T`. It is evaluated on the existing
  U.S. Standard Atmosphere 1976 pressure, `standard_atmosphere_pressure_pa`. No new atmosphere was
  needed, since density and temperature are not required. A helper picks the maximum q from sampled
  trajectory points.
- **Skin friction (Barrowman chapter 4).** Laminar `1.328/sqrt(R)`, Hama, transitional
  (`-1700/R`) and rough-wall (`51 (R_s/L)^-1.039`, `0.032 (R_s/L)^0.2`) relations, with Barrowman's
  subsonic and supersonic compressibility factors (equations 4-2 to 4-13). The tests confirm the
  thesis's claim that Hama agrees with Schoenherr within +-2% for Re = 1e6 to 1e9; the measured
  deviation is at most 1.7%.
- **Drag build-up.** Body friction drag with the `1 + 0.5/f_B` factor (equation 4-16) and subsonic
  forebody pressure drag (equations 4-41 to 4-46; the sonic value is reproduced to 0.03% because
  6/0.88 is printed rounded as 6.82). Base drag covers the subsonic range (equations 4-50 to 4-52)
  and the supersonic range (equations 4-54 to 4-56), with C*_DB = 0.185 and 0.30 at h/c = 0.1 as
  printed.
- **Loads.** Six body-axis coefficients are made dimensional as `q S C` and `q S l C`. A sectional
  line load is integrated to shear, bending and CP; this gives the aerodynamic part only.
- **SP-8035 gust and load combination.** The SP-8035 one-minus-cosine discrete gust follows
  Appendix A.2 (7.65 m/s, thickness 30-275 m). Random loads are root-sum-squared and superimposed on
  the quasi-steady load (section 4.2.2). The wind-induced incidence comes from adding the wind and
  vehicle velocity vectors.

### Database and uncertainty structure (`aero_database`)

- **Coefficient tables.** An N-dimensional rectilinear table, such as Mach, alpha and beta as in
  the SLS ascent look-up tables, uses multilinear interpolation. Queries outside the grid are
  refused.
- **Uncertainty bounds.** Uncertainty is either gridded or piecewise-constant over ranges of one
  axis. The SLS ascent bounds are separate for Mach 0.5-1.3 and 1.6-2.5, so a query at Mach 1.4 is
  refused. Each query returns `C_Q -+ U` (Favaregh et al. equation 1) and dispersed values
  `C_Q + k U` with `|k| <= 1`.
- **Uncertainty arithmetic.** The SLS papers' build-up is implemented:
  - the range method `sigma = R_bar / d2`, with d2 = 1.128, 1.693, 2.059, 2.326 (Lee et al.);
  - `U_EXP = 3 sqrt(n) sigma_rep` (Favaregh et al. equations 4-5);
  - 99.7% coverage bounds;
  - the net modelling term `sqrt(U_DBM^2 - U_rpt^2)`, or 0 (Shea et al. equations 2-3);
  - the factored root-sum-square `UF sqrt(sum w u^2)`, covering Shea et al. equation 4
    (UF = 1.05), Favaregh et al. equation 12 and Artemis I equations 5-7 (weight 2 for the doubled
    calibration term).
- **Hemsch fits.** Hemsch's least-squares slender-body-theory fits `C_N = C_N0 + C_N_alpha alpha`
  and `C_m = C_m0 + (x_CN - x_BMC) C_N` return the implied centre of pressure.

### Implementation choices not stated in the sources

- Multilinear interpolation of the tables. The papers say trajectory codes interpolate the
  delivered tables but do not name the interpolant.
- Interpolating gridded uncertainty bounds; the result carries an `info` warning.
- The order-statistic estimator for a "99.7% coverage" bound.
- The SP-8035 gust equation (1) slope. It is printed as `0.765 dH`, ten times the slope of
  equation (2) at the tangent point that the text defines; the tangent slope 0.0765 is used.
- Barrowman base drag. Equation 4-50 prints 0.29 where equations 4-48, 4-49 and 4-52 use 0.029, and
  equation 4-53 lacks a square; 0.029 and equation 4-52 are used. The two supersonic branches jump
  slightly at M_cr, and the result reports this jump as a warning.
- Barrowman equation 3-90 (tangent-ogive volume) is not used. As read, it gives one quarter of the
  integrated volume; the volume is integrated numerically instead.
- The fin aspect ratio is a required input. Barrowman's list of symbols defines it only in words.

## (b) Industrial practice documented, needs data the user must supply

- **SLS force-and-moment databases.** These are built from wind-tunnel tests (Ames UPWT, Boeing
  Polysonic, LaRC 14x22) and CFD (FUN3D, OVERFLOW, Cart3D, USM3D, Loci/Chem). The response surfaces
  are natural-neighbour interpolants of scattered data, and booster separation uses kriging over
  more than 20,000 Cart3D runs. This package takes the delivered tables and bounds as input. It does
  not build response surfaces and contains no SLS data; the published plots have their absolute
  values removed (SBU).
- **Ground-to-flight, tunnel-to-tunnel and code-to-code terms.** These come from paired CFD or test
  data. The package combines them but cannot produce them.
- **Total angle of attack and roll conversions.** The SLS liftoff/transition and ascent databases
  convert breakpoints between total angle of attack and roll and alpha-beta space. The papers read do
  not print the conversion formula, so it is not implemented.
- **Distributed line loads, buffet, protuberance loads and plume effects.** These come from CFD and
  tests. `integrate_line_load` accepts a user line-load table.
- **Structural ascent loads.** These need the mass and stiffness distribution, control-system and
  slosh models and a flight simulation (SP-8035 sections 2.2 and 4.2). Only the aerodynamic part and
  the combination rule are implemented.
- **Supersonic Barrowman methods** (second-order shock expansion, fin strip theory in appendix A)
  and the body-in-presence-of-tail CP curve (figure 3-10) are not implemented. Supply them as
  `NormalForceComponent` values.
- **Cross-flow drag versus cross-flow Mach and Reynolds number.** NACA TR 1048 figure 9 and the
  cylinder data it cites are not digitised. The user supplies `c_d`.

## (c) Not public

- **SpaceX launch-vehicle aerodynamic methods** are not published, and nothing about them is
  claimed or implemented. The only primary source read on SpaceX is Braun, Sforzo and Campbell
  (NTRS 20170008535). It states that SpaceX performed the first supersonic retropropulsion
  manoeuvre on a Falcon 9 first stage on 29 September 2013, and that NASA and SpaceX formed a
  three-year public-private partnership from FY2014 for flight reconstruction, CFD, imagery and Mars
  EDL design. It gives no implementable relations, so supersonic retropropulsion is out of scope.
- **SLS numerical database values and uncertainty magnitudes** are Sensitive But Unclassified, as
  stated in the notice of Shea et al. Only the uncertainty factor UF = 1.05 and the structure are
  printed.

## References

- NACA-TR-1048 — Allen and Perkins, NACA Report 1048, 1951. <https://ntrs.nasa.gov/citations/19930090962>
  (copy read: <https://www.pdas.com/refs/rep1048.pdf>)
- BARROWMAN-1967 — Barrowman, M.S. thesis, 1967. <https://ntrs.nasa.gov/citations/20010047838>
- NACA-TR-1307 — Pitts, Nielsen and Kaattari, NACA Report 1307, 1959. <https://ntrs.nasa.gov/citations/19930091008>
- NASA-SP-8035 — Wind Loads During Ascent, 1970. <https://ntrs.nasa.gov/citations/19700027637>
- SHEA-2020-SLS-ASCENT-DB — Shea et al., SLS ascent F&M database development. <https://ntrs.nasa.gov/citations/20200002382>
- FAVAREGH-2016-SLS-UQ — Favaregh, Houlden and Pinier, AIAA 2016-0795. <https://ntrs.nasa.gov/citations/20160007692>
- LEE-2024-ARTEMIS-UQ — Lee, Shea, Pinier and Chan, Artemis I UQ. <https://ntrs.nasa.gov/citations/20240000054>
- HEMSCH-2016-SLS-CFD-WT — Hemsch, AIAA 2016-0794. <https://ntrs.nasa.gov/citations/20160007695>
- BRAUN-2017-SRP — Braun, Sforzo and Campbell, SRP overview. <https://ntrs.nasa.gov/citations/20170008535>
- NASA-GRC-DYNAMIC-PRESSURE — <https://www1.grc.nasa.gov/beginners-guide-to-aeronautics/dynamic-pressure/>
- NASA-GRC-SPEED-OF-SOUND — <https://www1.grc.nasa.gov/beginners-guide-to-aeronautics/speed-of-sound-interactive/>
