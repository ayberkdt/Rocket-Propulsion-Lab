# Heat conduction

Module: `src/rocket_propulsion/heat_transfer/conduction.py` (standard library only).
Tests: `tests/test_heat_conduction.py` with reference values in
`tests/reference_data/heat_conduction.json` (every value carries its source locator).

The module is a literature-verified conduction layer for propulsion hardware: chamber and
nozzle walls, liners, insulation, heat-sink sections and fins. Every equation is taken from
a source that was opened while the module was written; the locators below use the printed
page numbers of Lienhard and Lienhard, *A Heat Transfer Textbook*, 6th ed., v6.00
(`AHTT-V6`; printed page = PDF page - 14). All results are one-dimensional textbook models
and carry no hardware qualification.

## Public API

Steady walls

- `WallGeometry(kind, inner_radius_m=0, area_m2=1, length_m=1)` with
  `kind` in `plane`, `cylinder`, `sphere` (inner radius 0 = solid body).
- `WallLayer(thickness_m, thermal_conductivity, volumetric_heat_capacity=None,
  heat_generation_w_m3=0, name="")`; conductivity and `rho c` are a float or a
  `TabulatedProperty(temperatures_k, values)` (piecewise linear, never extrapolated).
- Surface conditions: `SurfaceTemperature`, `SurfaceHeatFlux` (positive into the wall),
  `SurfaceConvection(h, T_fluid, emissivity=0, surroundings_temperature_k=None)`.
- `solve_steady_wall(layers, geometry, outer_boundary, inner_boundary=None, *,
  contact_conductances_w_m2_k=None)` returns `SteadyWallResult` (interface temperatures on
  both sides of each contact, heat rates and fluxes, sampled profile, maximum temperature,
  energy and boundary residuals).
- Resistances: `plane_wall_resistance`, `cylindrical_wall_resistance`,
  `spherical_wall_resistance` (outer radius may be `math.inf`), `convection_resistance`,
  `contact_resistance`, `radiation_heat_transfer_coefficient`,
  `critical_insulation_radius`.

Transients

- `lumped_capacitance_response(...)` returns the temperature history, time constant and
  Biot number `h (V/A) / k`, with a warning above `biot_limit` (default 0.1).
- Semi-infinite solid: `semi_infinite_surface_temperature_step`,
  `semi_infinite_surface_flux_step`, `semi_infinite_convective_surface`,
  `semi_infinite_periodic_surface`, `contact_temperature`, and `integrated_erfc`.
- Series solutions: `transient_eigenvalues(shape, Bi, count)`,
  `one_term_coefficients(shape, Bi)` and `series_transient_solution(shape, Bi, Fo,
  position_ratio)` for `slab`, `cylinder`, `sphere` (`Bi = math.inf` allowed).
- `solve_transient_wall(layers, geometry, outer_boundary, inner_boundary=None, *,
  initial_temperature_k, time_step_s, step_count, cells_per_layer=20, theta=0.5,
  contact_conductances_w_m2_k=None, snapshot_every=0)` returns `TransientWallResult`
  (surface temperature and heat-rate histories, snapshots, final field, energy audit).

Fins and special functions

- `straight_fin`, `pin_fin`, `annular_fin`, `corrected_fin_length`; results are
  `FinResult` (heat rate, efficiency, effectiveness, resistance, tip temperature,
  transverse Biot number).
- `bessel_j0`, `bessel_j1`, `bessel_i0`, `bessel_i1`, `bessel_k0`, `bessel_k1` and the
  exponentially scaled `bessel_i0e`, `bessel_i1e`, `bessel_k0e`, `bessel_k1e`.
- `STEFAN_BOLTZMANN_W_M2_K4 = 5.670374419e-8` (CODATA 2022). It is defined locally because
  `rocket_propulsion.core.constants` does not register it yet.

## Equations and sources

Steady conduction. The conduction equation (2.10) with the Laplacians (2.12)-(2.14b)
(printed pp. 55-56). The Kirchhoff variable `phi = int k dT` turns the steady equation into
`div grad phi + q''' = 0` (ENCIT 2012 paper, Sec. 3, Eq. (5)); for a plane wall this is the
book's `q L = int k dT` (Problem 2.2, printed p. 86). With the energy balance
`Q(r) = Q(r_a) + q''' V(r_a, r)` and Fourier's law, the drop across a layer
`phi(r_a) - phi(r_b) = int Q/A dr` is integrated in closed form for plane, cylindrical and
spherical areas. For constant `k` and no generation it reduces to the resistances (2.16)
and (2.22) and the spherical-shell resistance that the book poses as Problem 2.21
(derived here and checked against the printed answer of Problem 2.15). Contact
resistances add `Q/(h_c A)` (Sec. 2.3, printed pp. 75-77). Convective and radiative
surfaces use `h (T_inf - T_s) + eps sigma (T_surr^4 - T_s^4)`, eqns. (1.34)-(1.35); the
transfer factor `F = eps` assumes large isothermal surroundings and is flagged on concave
inner surfaces. The unknown surface temperature or heat rate is found by bisection on the
outer-boundary residual, which is monotonic because `k > 0`.

Lumped capacity. Eqs. (1.19)-(1.23) (printed pp. 21-22); Biot criterion `Bi <~ 0.1`
(printed p. 24) with the length `V/A` of eqn. (5.7). The book's Example 1.4 uses `D/2`,
which for a sphere is three times `V/A`.

Semi-infinite solid. Temperature step (5.50), (5.48), (5.54); penetration depth (5.51);
convective surface (5.52)-(5.53), evaluated as `erf(z) + exp(-z^2) erfcx(z + beta)`;
flux step (5.55)-(5.56), with the interior profile `T - T_i = (2 q_w / k) sqrt(alpha t)
ierfc(z)` obtained by the same integration the book performs at the surface and the
repeated integral `ierfc` of DLMF 7.18.2 (closed form from 7.18.1 and 7.18.7; large-z tail
from DLMF 7.12.1); contact temperature (5.60); periodic surface (5.62).

Series solutions. Table 5.1 eigenvalue equations, `A_n` and `f_n`; eqns. (5.36), (5.38),
(5.40). `D_n` is the volume average of `A_n f_n`: `A_n sin(lambda)/lambda` (printed after
(5.38)) for the slab and, by the same averaging, `2 A_n J_1(lambda)/lambda` (cylinder) and
`3 A_n (sin lambda - lambda cos lambda)/lambda^3` (sphere). All three are checked against
every row of Table 5.2 and against a direct Simpson integration of the series.

Finite volumes. Cell-centred finite volumes with uniform cells per layer; theta time
integration (0.5 Crank-Nicolson, 1 backward Euler); stored energy `int rho c dT`; face flux
inside a layer from the Kirchhoff secant conductivity over the centre-to-centre shell
integral; series half-cell and contact resistances across interfaces; boundary half cells
solved together with the surface condition. Nonlinearities (k(T), rho c(T), radiation) are
iterated within each step with Newton-linearized storage and boundary terms. Each internal
face flux enters its two cells with opposite signs, so the scheme is conservative and the
result reports the energy residual.

Fins. One-dimensional straight fin, eqns. (4.27)-(4.58) (printed pp. 163-175): adiabatic
tip (4.41)-(4.45), convective tip (4.48)-(4.49), efficiency (4.52)-(4.53), effectiveness
(4.54)-(4.55) and resistance (4.56)-(4.58). Tip correction `L + t/2` from NACA Report 158
(printed p. 26, assumption (3); Eq. (5) on printed p. 21). Annular fin from Yovanovich,
Culham and Lemczyk (1988), Eqs. (26)-(41), in dimensional radius with `m = sqrt(h/(k t))`
(`t` = half thickness) and an optional edge coefficient. Their printed Eqs. (32) and (35)
contain a factor `m` on `Omega K_0(m alpha)` that is inconsistent with their Eq. (29) at the
root; the module applies the root condition to (29) directly, and the test confirms the
result against an independent finite-difference solution.

Bessel functions (DLMF Ch. 10). `I_0, I_1`: power series 10.25.2 for x <= 50, expansion
10.40.1 beyond. `K_0, K_1`: series 10.31.1-10.31.2 (digamma from 5.4.14, Euler's constant
5.2.3) for x <= 2, the integral 10.32.9 with the trapezoidal rule (DLMF 3.5.5) for
2 < x <= 50, expansion 10.40.2 (remainder bounded by the first omitted term) beyond.
`J_0, J_1`: Bessel's integral 10.9.2 with the periodic trapezoidal rule (DLMF 3.5.4).
Exponentially scaled forms avoid overflow in the annular-fin constants.

## Verification evidence (measured in this task)

Textbook tables

- Table 5.2: all 41 rows for slab, cylinder and sphere reproduce `lambda_1` within
  6e-6 and `A_1`, `D_1` within 6e-5 (the printed rounding).
- Table 5.3: all 26 rows of erf and erfc within 1.1e-5. Eqn. (5.33) printed terms
  reproduced within 6e-4.

Worked examples and printed answers (computed / printed)

- Example 1.5 black thermocouple 28.40 / 28.4 degC; Example 1.6 23.52 / 23.5 degC
  (exact radiative surface balance).
- Example 1.4 time constant 1.116 s; Example 2.8 53.1 s.
- Example 2.5: eqn. (2.24) profile and (2.25) heat rate to 1e-12.
- Example 2.6: critical radius 0.0037 m, ratio 1.48, break-even ratio 2.32.
- Example 2.7 resistor 72.32 / 72.3 degC with exact radiation. The printed
  `h_rad = 7.17` at 327 K is not reproduced by eqn. (2.31) (7.14, 0.45 percent lower).
- Example 2.9 contact share 0.091 / "about 10%"; Example 2.10 U = 192.08 / 192.1.
- Problem 2.15 Q = 16.89 / 16.9 W; Problem 2.20 centre 21.367 / 21.37 degC.
- Example 4.8 tip Theta 0.6889 / 0.6886 (the printed value is a rounding slip of
  1/(1.3986 + 0.0529) = 0.6889) and 0.7150 / 0.7150; Example 4.9 fin resistance
  2149.7 / 2150 K/W, equivalent 276.78 / 276.8 K/W, resistor 62.68 / 62.68 degC.
- Example 4.10 annular efficiency 0.908 / 0.89 (chart read-off), finned heat loss
  484.7 / 478 W/m.
- Example 5.2 centre Theta 0.860 / 0.84 (chart), Fo for Theta = 0.2 1.289 / 1.3 (chart),
  Phi 0.830 / 0.82 (chart).
- Example 5.4 0.045 s; Example 5.5 burn time 0.335 / 0.33 s; Example 5.7 depth
  2.783 / 2.783 m.

Bessel functions

- Wronskian DLMF 10.28.2 satisfied within 5e-15 for x from 1e-3 to 600.
- The series and integral representations of K agree within 5e-15 for x <= 2.
- J agrees with the independent power series DLMF 10.2.2 within 2e-15 for x <= 6;
  derivative identities 10.29.3 and 10.6.3 hold to finite-difference accuracy; the first
  zero of `J_0` matches Table 5.2 (2.40483).

Finite-volume convergence (observed orders)

- Spatial, against the series at Bi = 1, Fo = 0.5, 10/20/40 cells, Crank-Nicolson with
  dt proportional to dx: slab field 1.98 and 1.99, surface 2.00; cylinder field (max
  norm, dominated by the centre cell) 1.79 and 1.82, surface 2.01; sphere field 1.81 and
  1.83, surface 1.98 and 1.99.
- Temporal, from an eigenmode start on a fixed grid: Crank-Nicolson 2.000, backward
  Euler 0.977 then 0.989.
- Semi-infinite constant flux, surface temperature against eqn. (5.56): 2.00 (relative
  error 1.6e-4 at 200 cells).
- Nonlinear two-layer cylinder (tabulated k and rho c, generation, contact conductance,
  convective inner surface, convective and radiative outer surface) against the
  independent Kirchhoff steady solver: order 2.0 on 4/8/16 cells per layer, heat rate
  within 1e-5.
- A single layer with tabulated k reaches the Kirchhoff steady heat rate to 1e-12 at any
  grid size.
- Energy audit: residual below 1e-9 of the stored-energy change in every test run.
- Annular fin: agrees with a 4000-node finite-difference solution within 2e-6 (insulated
  and convective edge); tends to `tanh(mL)/mL` as the root radius grows.

## Limits and what is not verified

- All models are one-dimensional; ribs, channels, axial conduction in nozzle walls and
  two-dimensional spreading are not modelled. Boundary conditions of the transient solver
  are constant within a run; time-varying conditions are handled by chaining runs.
- The radiation transfer factor `F = eps` holds only for a surface that sees large
  isothermal surroundings.
- No rocket-specific temperature-dependent property data are built in; `k(T)` and
  `rho c(T)` are supplied by the caller. The GRCop-84 report NASA/CR-2000-210055 cited by
  `thermodynamics/materials.py` was opened, but its k(T) data were not tabulated here.
- No NASA transient thrust-chamber wall-heating report was used to validate the solver.
- Abramowitz and Stegun Bessel tables were not used: the copy obtained lacks the
  numerical tables of Ch. 9. Bessel values rest on DLMF identities, two independent DLMF
  representations and Table 5.2.
- The fin transverse-Biot warning level 0.1 is a library choice; the sources only require
  `Bi << 1`.
- Crank-Nicolson can ring after a sudden boundary change when `alpha dt / dx^2` is large;
  backward Euler (`theta = 1`) is monotone but first order in time.

## References

- AHTT-V6: <https://ahtt.mit.edu/wp-content/uploads/2024/04/AHTTv600.pdf> (landing page
  <https://ahtt.mit.edu/>)
- ENCIT-2012-KIRCHHOFF: <https://www.abcm.org.br/anais/encit/2012/links/pdf/ENCIT2012-352.pdf>
- DLMF-10-BESSEL: <https://dlmf.nist.gov/10>
- DLMF-7-ERFC: <https://dlmf.nist.gov/7.18>, <https://dlmf.nist.gov/7.12>
- DLMF-5-GAMMA: <https://dlmf.nist.gov/5.2>, <https://dlmf.nist.gov/5.4>
- DLMF-3.5-TRAPEZOID: <https://dlmf.nist.gov/3.5>
- YOVANOVICH-1988-ANNULAR-FIN: <https://uwaterloo.ca/microelectronics-heat-transfer-laboratory/sites/default/files/2024-03/mhtl88-1.pdf>
- NACA-TR-158: <https://ntrs.nasa.gov/api/citations/19930091223/downloads/19930091223.pdf>
- NIST-CODATA-2022-STEFAN-BOLTZMANN: <https://physics.nist.gov/cgi-bin/cuu/Value?sigma>
- PYTHON-3.12-MATH-ERF: <https://docs.python.org/3.12/library/math.html#math.erfc>
- AANDS-UBC-PDF (consulted, not used): <https://personal.math.ubc.ca/~cbm/aands/abramowitz_and_stegun.pdf>

Full entries: [REFERENCES.md](REFERENCES.md).
