# Advanced compressible (inviscid) flow

Five modules in `rocket_propulsion.compressible` extend the calorically perfect
one-dimensional relations (`isentropic`, `normal_shock`, `oblique_shock`) to
thermally perfect gases, conical flow, two-dimensional supersonic sections,
linearized compressibility corrections and inviscid hypersonic approximations.
They are stdlib-only, SI, fail closed with `DomainError`, and return frozen
dataclasses. Every equation was taken from a source opened in this work
package; locators below and in each docstring. Book pages are printed page
numbers of Anderson (6th ed.); the user's PDF is offset by +23 pages.

## Modules and public API

### `thermally_perfect.py` — calorically imperfect, thermally perfect gas

- `HarmonicVibratorGas(gamma_perfect=1.4, vibrational_temperature_k=5500 degR, gas_constant_j_kg_k)`:
  NACA 1135 diatomic gas with a harmonic-vibrator heat capacity, eqs. (176), (178),
  (180), (184), (189) (NACA 1135 pp. 627-628). Methods `cp_j_kg_k`, `gamma`,
  `enthalpy_j_kg`, `entropy_function_j_kg_k`.
- `PolynomialGas(source)`: frozen-composition wrapper around a repository
  `Species` or `Mixture` (NASA Glenn polynomials); fails closed outside the
  polynomial temperature range.
- `thermally_perfect_isentropic(mach, total_temperature_k, gas)`:
  energy equation `h(Tt) - h(T) = M^2 a^2/2` (eq. (186)), `p/pt = exp((phi(T)-phi(Tt))/R)`
  (eq. (189)), `A/A* = rho* a*/(rho V)` (eq. (191)).
- `thermally_perfect_normal_shock(upstream_mach, upstream_temperature_k, gas)`:
  continuity, momentum and energy (eqs. (192)-(194)) solved for `T2` by bisection.
- `thermally_perfect_oblique_shock(upstream_mach, turning_angle_deg, upstream_temperature_k, gas)`:
  normal-shock relations with `M1 sin(beta)` and the same `T1` (NACA 1135 pp. 629-630),
  `tan(beta - delta) = tan(beta)/(rho2/rho1)`, weak branch.

Dissociation is not modelled (NACA 1135 places its onset near 5000 degR).

### `cone_flow.py` — Taylor-Maccoll conical flow

- `cone_flow_from_shock_angle(mach, shock_angle_deg, gamma)`: inverse problem,
  Anderson section 13.6.3 (pp. 865-866): oblique-shock state behind the shock,
  `V'` from eq. (13.81), fourth-order Runge-Kutta integration of eq. (13.80)
  to `V'_theta = 0`.
- `solve_cone_flow(mach, cone_half_angle_deg, gamma)`: direct problem by
  bisection on the shock angle, bracketed by the Mach angle and the attached
  wedge wave angle (cone shock weaker than wedge shock, Anderson pp. 634-636).
- `maximum_cone_half_angle(mach, gamma)`: detachment limit (Anderson Fig. 13.20).

The integration step is `2e-4` rad, reduced to 1 % of `1 - (V_theta/a)^2`
near a weak shock where eq. (13.80) is nearly singular; results are
step-converged to about `1e-12`.

### `shock_expansion.py` — two-dimensional polygonal sections

- `AirfoilSection(upper, lower)`, builders `flat_plate()`, `diamond_airfoil(half_angle_deg)`,
  `wedge_section(half_angle_deg)`, `parabolic_biconvex(thickness_ratio, panels)`
  (the parabolic arc is this library's geometric choice; Anderson only sketches
  a bi-convex section, Fig. 12.3).
- `shock_expansion_coefficients(section, mach, angle_of_attack_deg, gamma, base_pressure_ratio, moment_reference)`:
  Anderson section 9.7 — weak oblique shocks at compression corners,
  Prandtl-Meyer fans (eq. (9.43)) at expansion corners, exact panel force
  integration (eqs. (1.15)-(1.19)); returns `c_l`, `c_d`, `c_n`, `c_a`,
  `c_m` (nose up, about the leading edge by default) and the panel solution.
- `integrate_section_forces(...)`: the shared exact integrator.

Base convention: `base_pressure_ratio = 1` puts free-stream pressure on a
blunt base (Anderson Example 9.6); Anderson Example 9.13b neglects the base
force, which corresponds to `base_pressure_ratio = 0`.

### `linearized_flow.py` — linear theory and corrections

- `prandtl_glauert`, `karman_tsien`, `laitone`, `corrected_pressure_coefficient`
  (Anderson eqs. (11.51), (11.54), (11.55)).
- `critical_pressure_coefficient(mach, gamma)` (eq. (11.60)),
  `critical_mach_number(cp0_min, rule, gamma)` (eq. (11.63) solved by bisection;
  for Karman-Tsien and Laitone the first intersection before the rule's own
  singularity is returned).
- `ackeret_pressure_coefficient` (eq. (12.15)), `ackeret_flat_plate`
  (eqs. (12.23)-(12.24)).
- `busemann_coefficients`, `busemann_pressure_coefficient`: `C1 = 2/sqrt(M^2-1)`,
  `C2 = ((gamma+1)M^4 - 4(M^2-1))/(2(M^2-1)^2)`, the first two terms of NACA 1135
  eq. (152); eq. (174) shows the expansion series shares them.
- `linearized_section_coefficients(section, mach, alpha, order=1|2)`: Ackeret or
  Busemann panel `Cp` integrated exactly over an `AirfoilSection`.
- `subsonic_validity_warnings`, `supersonic_validity_warnings`: Anderson's
  `M >= 0.7` Prandtl-Glauert caution (p. 751), the transonic band
  `0.8 <= M <= 1.2` (p. 757) and Van Dyke's `beta tau << 1` criterion (the
  0.3 flag level is this library's choice, not a sourced number).

### `hypersonic_inviscid.py` — Newtonian theory, limits, similarity

- `newtonian_pressure_coefficient` (eq. (14.4), shadow eq. (14.9)),
  `stagnation_pressure_coefficient` (`Cp,max`, p. 880),
  `modified_newtonian_pressure_coefficient` (eq. (14.7)).
- `newtonian_flat_plate` (eqs. (14.8)-(14.16)), `newtonian_maximum_lift` (p. 886),
  `newtonian_maximum_lift_to_drag` (eqs. (14.23)-(14.24), `2/3` printed as 0.67).
- `hypersonic_shock_limits` (eqs. (14.28)-(14.32), (14.38)-(14.39)),
  `limiting_wave_angle_ratio` (eq. (14.35)), `mach_independent_wedge`
  (the `M -> infinity` limit of eq. (9.23) with eq. (14.39), section 14.6).
- `hypersonic_similarity_parameter` (`M tau` or unified `beta tau`, Van Dyke p. 6),
  `small_disturbance_wedge` (Van Dyke eqs. (33)-(35), Linnell's thin wedge;
  `delta = tan(wedge angle)` as in Van Dyke Fig. 5).

Anderson's chapter 14 does not attribute the modified law to Lees by name;
the module uses the equation, not the attribution.

## Validation (measured against the sources)

All reference numbers live in `tests/reference_data/compressible_advanced.json`
with locators; `tests/test_compressible_advanced.py` runs in about 8 s.

- Ratio of specific heats, NACA 1135 p. 627 table (25 temperatures, 500-5000 degR):
  largest deviation from eq. (180) is 6.2e-4 (3000 degR); tolerance is one unit
  of the third printed decimal.
- Isentropic imperfect-gas charts 11, 12, 14 (pixel-digitized plateau values):
  deviations -0.0006 to +0.0023 (tolerance 0.003, 0.004 where the 5000 degR
  curve overlaps a grid line).
- Normal-shock imperfect-gas charts 15, 16, 18, 19, 20 at `M1 = 10` (visual
  readings, tolerance one minor grid division): T2/T1 -0.0085/-0.0082/-0.0047
  (division 0.01), rho2/rho1 up to +0.0032 (0.005), p2/p1 up to -0.0016 (0.002),
  M2 up to -0.0011 (0.002), pt2/pt1 up to -0.0047 (0.01). The chart-15 offsets
  are the largest; readings of charts 16 and 18 imply 0.958 for 2000 degR,
  closer to the model's 0.9545 than the chart-15 reading.
- Limits: `Theta/T = 30` reproduces the calorically perfect isentropic and
  normal-shock ratios to 1e-8; a vibrational temperature of 1e6 K reproduces
  `calculate_oblique_shock` to 1e-7 degrees.
- Sims cone tables (15 cases, 5-30 degrees, M = 1.5-10): shock angle within
  7.2e-6 rad, surface `Cp` within 2.2e-5, surface Mach within 1.6e-5 relative.
  These offsets exceed Sims' stated integration error (5e-8); Table 70 confirms
  the free-stream Mach number is exact, and this solver's post-shock state
  satisfies Sims eqs. (5)-(6) to machine precision, so the residual is not
  explained. Anderson Fig. 9.16 (M = 2, 20 degree cone, 37 degrees printed):
  37.80 degrees.
- Shock-expansion, Anderson examples (book values use chart-read wave angles):
  Example 9.6 `c_d` +1.6 %, Example 9.12 `c_l` -0.5 % and `c_d` -1.1 %,
  Example 9.13a `c_n` -1.8 %, Example 14.1(a) `c_l`, `c_d`, `c_n` about -1.5 %,
  `L/D` exact. Feeding the book's printed panel pressures to the force
  integrator reproduces its printed coefficients to the last digit.
- Linear theory: Examples 11.1-11.3, 12.1-12.3 reproduced to printed digits,
  including the `Cp,cr` table (p. 762), the eq. (11.63) iteration table
  (p. 763, four decimals) and `Mcr = 0.7371`. Busemann `Cp` error scales as
  `theta^3`, Ackeret's as `theta^2`, against the exact oblique shock.
- Hypersonic: Example 14.1(b) reproduced to printed digits; `c_l,max = 0.770`
  at 54.74 degrees; stagnation `Cp` = 1.2756 at M = 1 (printed 1.28).

## Discrepancies found in the sources

- Anderson Example 9.13b reads `beta = 23.5` degrees for M = 7, theta = 15 degrees;
  the theta-beta-M relation gives 21.60 degrees, so `p3/p1` is 7.58 rather than
  8.915 and the exact `c_l` is 0.203, not 0.241. Only the replay (printed
  pressures through the integrator) is tested for this example.
- Anderson Example 9.13a labels 0.126 as `c_l`, but its arithmetic omits the
  `cos(alpha)` of its own formula; 0.126 is `c_n` (`c_l = 0.1242` from the
  printed pressures).
- Anderson p. 879 gives the `M -> infinity` stagnation `Cp` as 1.86 for gamma = 1.4;
  the book's own normal-shock and isentropic relations give 1.8394. The test
  uses a 0.025 tolerance for this value and says so.
- Example 12.1 cites "Example 9.11" for the exact flat-plate result; in this
  edition that example is numbered 9.12.

## Not verified

- No published numbers were found for the moment coefficient; `c_m` is
  checked only by the flat-plate centre of pressure at mid-chord.
- The imperfect-gas oblique-shock chart (NACA 1135 chart 21-23) and
  Prandtl-Meyer chart 24 were not digitized; the oblique shock is verified
  through its normal-component construction and its perfect-gas limit.
- `PolynomialGas` is checked only near room temperature for N2 against the
  calorically perfect limit; no imperfect-gas air composition is shipped.
- Busemann third-order terms are not implemented: NACA 1135 shows they
  differ between compression (eq. (151)) and expansion (eq. (174)).
- Thickness and camber wave-drag closed forms (`g_c`, `g_t` in Anderson
  p. 818) are not given in the source and are not implemented; the panel
  integration is used instead.

## References

- NACA-TR-1135 — Ames Research Staff, Equations, Tables, and Charts for Compressible Flow, NACA Report 1135, 1953. <https://ntrs.nasa.gov/citations/19930091059>
- NASA-SP-3004 — J. L. Sims, Tables for Supersonic Flow Around Right Circular Cones at Zero Angle of Attack, 1964. <https://ntrs.nasa.gov/citations/19640009035>
- NACA-TR-1194 — M. D. Van Dyke, A Study of Hypersonic Small-Disturbance Theory, 1954. <https://ntrs.nasa.gov/citations/19930092209>
- NASA/TP-2002-211556 — NASA Glenn coefficients for species properties (used through `PolynomialGas`). <https://ntrs.nasa.gov/citations/20020085330>
- ANDERSON-FOA-2017 — J. D. Anderson Jr., Fundamentals of Aerodynamics, 6th ed., McGraw-Hill, 2017, ISBN 978-1-259-12991-9 (bibliographic citation).
