# Inviscid incompressible and lifting aerodynamics

Five modules in `rocket_propulsion.aerodynamics` cover two-dimensional potential
flow, NACA section geometry, thin airfoil theory, panel methods and Prandtl's
lifting-line theory. They are stdlib-only, SI, fail closed with `DomainError` /
`InputError`, and return frozen dataclasses. Every equation was taken from a
source opened in this work package; locators are below and in each docstring.
Anderson page numbers are printed book pages of the 6th edition (the user's
PDF is offset by +23). Reference values with their locators are in
`tests/reference_data/aero_subsonic.json`; the tests are in
`tests/test_aero_subsonic.py` (35 tests, about 3 s).

## Modules and public API

### `potential_flow.py` — elementary flows, cylinder, Joukowski sections

- `UniformFlow(speed, angle_rad)`, `SourceFlow(strength, x0, y0)`,
  `DoubletFlow(strength, x0, y0)`, `VortexFlow(circulation, x0, y0)`: each with
  `velocity`, `potential`, `stream_function` (Anderson Eqs. 3.53-3.57, 3.62,
  3.67, 3.72, 3.87-3.88, 3.105, 3.112, 3.114). Vortices are clockwise positive,
  as in Anderson, so `L' = rho V Gamma` is positive upward.
- `superpose_velocity`, `superpose_potential`, `superpose_stream_function`.
- `pressure_coefficient(local_speed, freestream_speed)` (Eq. 3.38) and
  `pressure_coefficient_from_pressure(...)` (Eq. 3.36).
- `kutta_joukowski_lift_per_span(density, speed, circulation)` (Eq. 3.140).
- `CircularCylinderFlow(freestream_speed, radius, circulation)`: field and
  surface velocity, surface `Cp`, `c_l = Gamma/(R V)`, stagnation points in all
  three regimes, and the surface angles where `p = p_inf`
  (Eqs. 3.92-3.138).
- `JoukowskiAirfoil(scale, center_x, center_y, angle_of_attack_rad, radius=None, freestream_speed=1)`:
  exact flow over the conformal image of a circle under
  `zeta = z' + a^2/z'` (NACA TR-411 Eq. I), with the circle flow of Anderson
  Eq. 3.125 rotated by the incidence, Kutta circulation (TR-411 Eq. VII) and
  the velocity transformation of TR-411 p. 232. `radius=None` gives the
  classical cusped section; a larger radius gives a smooth section with a
  rounded trailing edge, for which the rear stagnation point is placed at the
  rearmost point (TR-411 footnote 4, "or some equivalent assumption").
- `ellipse_surface_speed_ratio(psi0, alpha, theta)`: exact elliptic-section
  surface speed in the form of TR-411 Eq. XII.

### `airfoil_geometry.py` — NACA four- and five-digit sections

- `naca_section(designation) -> NacaSection` with `thickness(x)`,
  `mean_line(x)`, `upper_surface(x)`, `lower_surface(x)`,
  `leading_edge_radius`, `panel_nodes(n, close_trailing_edge=True)`.
- Thickness form and four-digit mean lines from NACA TR-460 pp. 4-5; surface
  construction from TR-824 Eqs. 1-4; five-digit simple (210-250) and reflexed
  (221-251) mean lines from TR-537 pp. 521-522 and Table I; first-digit scaling
  from TR-824 p. 5. The 211 line is refused (TR-537 withholds its data);
  reflexed lines are accepted only with first digit 2 (the only tabulated case).
- The TR-537 Table I `k2/k1` column is printed as 1e4 times the ratio given by
  the TR-537 formula `(3 (m - p)^2 - m^3)/(1 - m)^3`; the formula is used.
- `panel_nodes(..., close_trailing_edge=True)` closes the finite TR-460
  trailing-edge thickness by subtracting the ramp `x y_t(1)` from the
  half-thickness. This is a library modification, not part of the NACA
  definition (at most `0.0105 t` of chord, at the trailing edge).

### `thin_airfoil.py` — classical thin airfoil theory

- `thin_airfoil(camber_slope, alpha, slope_breakpoints=(), coefficient_count=4)`
  and `thin_airfoil_naca(section, alpha)` return `ThinAirfoilResult`:
  Fourier `A0` and `A1..AN` (Eqs. 4.50-4.51), `c_l` (4.56), lift slope `2 pi`
  (4.58), `alpha_L=0` (4.61), `c_m,le` (4.62), `c_m,c/4` (4.64), `x_cp/c` (4.66),
  the ideal angle (`A0 = 0`) and ideal lift coefficient `pi A1`.
- Warnings: thickness above 12 % (Anderson Sec. 4.10 states favourable
  agreement up to about 12 %) and `|alpha| > 0.15 rad` (library threshold for
  the small-angle assumption of Eq. 4.14).

### `panel_methods.py` — source and lifting panel methods

- `panel_geometry(nodes)`, `influence_matrices(geometry)` (Anderson
  Eqs. 3.163, 3.165), `trailing_edge_angle_deg(geometry)`.
- `solve_source_panels(nodes, freestream_speed, angle_of_attack_rad)`:
  Anderson Sec. 3.17, Eqs. 3.153-3.157.
- `solve_lifting_panels(nodes, alpha, freestream_speed, chord, moment_reference)`:
  constant-strength source panels plus a vortex sheet of uniform unknown
  strength on all panels, closed by the Kutta condition in the form of Eq. 4.81
  applied to the two trailing-edge panel velocities. The vortex influence
  integrals are those of Eqs. 4.72-4.77; with substitution 3.160-3.161 they
  reduce to the source integrals (`d theta/d n` has the integrand of
  `d ln r/d s`, and `d theta/d s` minus that of `d ln r/d n`); the test suite
  checks both closed forms against quadrature. Lift from `rho V Gamma`
  (Eqs. 4.82-4.83) and from pressure integration (Eqs. 1.15-1.19); moment
  about any point, nose-up positive.
- Nodes run clockwise from the trailing edge (lower surface first). Open
  contours, counterclockwise contours and trailing-edge angles below 2.5°
  (treated as cusps) are refused.

### `lifting_line.py` — Prandtl lifting line with Glauert's series

- `trapezoidal_wing(span, root_chord, taper_ratio, tip_twist_rad, zero_lift_angle_rad, section_lift_slope)`,
  `elliptic_wing(...)`, `custom_wing(span, chord, twist, zero_lift_angle, section_lift_slope, breakpoints)`.
- `solve_lifting_line(wing, alpha_root, terms=40)`: collocation of Anderson
  Eq. 5.51 in the general-slope form of NACA TR-572 Eq. 5; returns `A_n`,
  `C_L` (5.53), `C_D,i` (5.60), `delta` (5.61), `e` (5.62), and spanwise
  `Gamma/V`, `c_l` and `alpha_i` (5.57).
- `lift_slope(wing)`: `a`, `tau` (Eq. 5.70 solved for `tau`), `delta`, `e` for
  an untwisted wing with uniform section data.
- `finite_wing_lift_slope(a0, AR, tau)` (5.69-5.70) and
  `induced_drag_coefficient(C_L, AR, delta)` (5.43, 5.61).

## Verification and measured deviations

| Check | Reference | Result |
| --- | --- | --- |
| Lifting cylinder, `c_l = 5` | Anderson Ex. 3.16-3.17 (pp. 278-280) | `Cp_min = -6.816` (-6.82); stagnation 203.45°/336.55° (203.4/336.6); `p = p_inf` at 5.86/174.14/243.88/296.12° (5.85/174.1/243.8/296.23 — the printed 296.23 is inconsistent with 243.8 by symmetry); bottom `Cp = -0.450` (-0.45) |
| Lift per span | Anderson Ex. 3.18 (pp. 281-282) | `Gamma = 39.270` (39.27), `L' = 892.66 N/m` (892.7) |
| Source panels, 8-panel cylinder | Anderson Ex. 3.19 (pp. 290-293) | all printed `I_4j` and `lambda_j/(2 pi V)` reproduced to the 4 printed digits |
| Source panels, cylinder `Cp` | Eq. 3.101 | exact to 1e-11 for 8 and 32 panels (symmetric paneling) |
| Thin airfoil, flat plate 5° | Anderson Ex. 4.5 | `c_l = 0.5483`, `c_m,le = -0.1371` (0.5485, -0.137 with alpha rounded) |
| Thin airfoil, NACA 23012 | Anderson Ex. 4.6 | `alpha_L=0 = -1.094°` (-1.09), `c_l(4°) = 0.5586` (0.559), `A1 = 0.0955` (0.0954), `A2 = 0.0792`, `c_m,c/4 = -0.01284` (-0.0127), `x_cp = 0.2730` (0.273); vs experiment (-1.1°, 0.55, -0.01) |
| Thin airfoil, NACA 64 mean line | TR-824 p. 14 | `c_m,c/4 = -0.1594` (-0.157, +1.5 %), `c_li = 0.768` (0.76, +1.1 %), `alpha_i = 0.772°` (0.74°, +0.032° — measured deviation, origin not determined); NACA 2415 `alpha_L=0 = -2.077°` vs -2.06° from TR-824 Eq. 16 with its printed `alpha_i`, `c_li`; NACA 4415 `c_m,c/4 = -0.1062` (-0.105) |
| Thin airfoil, NACA 2412 | Anderson p. 330 (data, read from Fig. 4.10) | `alpha_L=0 = -2.077°` vs -2.1° |
| NACA geometry | TR-824 pp. 100-101 tables (1412, 0006, 23024) | ordinates within 0.0025 % chord of the printed values; LE radii 0.396/0.891/6.336 % (printed 0.40/0.89/6.33) |
| Five-digit mean lines | TR-537 Table I, TR-610 p. 2, Anderson Ex. 4.6 | 230 line: `k1/6 = 2.6595`, camber 1.839 % at 0.15c (1.8 %); simple and reflexed lines peak at the tabulated `p` (to 0.0005c) |
| Lifting panels vs exact ellipse (`t/c = 0.12`, 5°) | TR-411 Eqs. VII, XII | `c_l` error 0.031, 0.016, 0.008 for 40/80/160 panels (first order) |
| Lifting panels vs exact rounded-TE Joukowski section | TR-411 + Anderson Eq. 3.125 | `c_l` 1.1227, 1.0951, 1.0837, 1.0788 for 40/80/160/320 panels vs exact 1.0745 (0.4 % at 320) |
| Lifting panels, NACA 2412 (160 panels) | Anderson p. 330, Ex. 4.1 | `alpha_L=0 = -2.14°` (measured -2.1°); `c_l(4°) = 0.743` vs measured 0.65 (+14 %; inviscid slope 0.121/deg vs 0.107/deg measured) |
| Lifting panels, NACA 23012 (160 panels) | Anderson Ex. 4.6-4.7 | `alpha_L=0 = -1.17°` (-1.1°), `c_l(4°) = 0.625` (0.55, +14 %); panel `x_ac ≈ 0.26` vs measured 0.241 |
| Lifting line, elliptic wing | Eqs. 5.42-5.43, 5.69 | exact to round-off (`delta = 0`, constant `alpha_i`) |
| Lifting line helpers | Anderson Ex. 5.1, 5.2 | `a = 4.97/rad`, `C_L = 0.4335`, `C_D,i = 0.00789`; Ex. 5.2 values to the printed digits |
| Lifting line `delta` vs taper | Anderson Fig. 5.20 (read ±0.003) | trend and minimum (taper ≈ 0.36-0.4 vs read ≈ 0.3) agree, but the computed values are 0.0285/0.0483/0.0676/0.0859 at taper 1 for AR 4/6/8/10 against 0.042/0.062/0.085/0.105 read from the figure, and 0.0444 at AR 8, taper 0.8 against Anderson's 0.055 — about 70-80 % of the figure. Not resolved; pinned in the tests |

## Not verified / not provided

- Anderson's first-order vortex panel method with one strength per panel and a
  dropped control point (Eqs. 4.80-4.81) was tried; on thin or cusped trailing
  edges its lift did not converge with refinement, so it is not provided.
- The lifting panel method does not converge for cusped trailing edges
  (classical Joukowski sections) and refuses them; the exact cusped Joukowski
  solution itself is provided and tested (Kutta circulation, flat-plate limit).
- Kármán–Trefftz sections are not implemented (no primary source opened).
- The origin of the Fig. 5.20 discrepancy (McCormick's data, cited by
  Anderson) was not traced; NACA TR-572 Fig. 4 (rounded tips, `m0 = 5.79`)
  was read only qualitatively and is not used as an oracle.
- Experimental `c_m,ac` of the NACA 2412 (TR-824 p. 136 plot) was not
  digitised; its printed aerodynamic-centre positions (0.239-0.247) are quoted
  but not tested.
- No viscous, compressibility or stall effects; panel and lifting-line results
  are inviscid.

## References

Full records are in [REFERENCES.md](REFERENCES.md) (bibliographic
citation for the book; NTRS links for the NACA reports):

- ANDERSON-FOA-2017 — J. D. Anderson Jr., *Fundamentals of Aerodynamics*, 6th ed., McGraw-Hill Education, 2017, ISBN 978-1-259-12991-9 (bibliographic citation).
- NACA-TR-824 — Abbott, von Doenhoff, Stivers, Summary of Airfoil Data, 1945. <https://ntrs.nasa.gov/citations/19930090976>
- NACA-TR-460 — Jacobs, Ward, Pinkerton, 1933. <https://ntrs.nasa.gov/citations/19930091108>
- NACA-TR-537 — Jacobs, Pinkerton, 1936. <https://ntrs.nasa.gov/citations/19930091610>
- NACA-TR-610 — Jacobs, Pinkerton, Greenberg, 1937. <https://ntrs.nasa.gov/citations/19930091685>
- NACA-TR-411 — Theodorsen, Theory of Wing Sections of Arbitrary Shape, 1931. <https://ntrs.nasa.gov/citations/19930091485>
- NACA-TR-572 — R. F. Anderson, Determination of the Characteristics of Tapered Wings, 1937. <https://ntrs.nasa.gov/citations/19930091647>
- DLMF-3.5 — NIST DLMF Sec. 3.5(v), Gauss quadrature, Table 3.5.1. <https://dlmf.nist.gov/3.5>
