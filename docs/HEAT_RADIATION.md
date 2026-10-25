# Thermal radiation

Package: `rocket_propulsion.heat_transfer.radiation`. Python standard library only,
SI units, frozen slotted result records, typed `DomainError`/`InputError`/
`ConvergenceError` on invalid input (fail closed), `WarningMessage` records for every
accepted extrapolation. Every equation and coefficient below was read from the cited
source in this work package; page and table locators are in the docstrings and in
`tests/reference_data/heat_radiation.json`.

## Scope

**Included**

- Black body: `e_b = sigma T^4` and Planck's law `e_lambda,b = 2 pi h c^2 /
  (lambda^5 [exp(hc/(k T lambda)) - 1])` (Lienhard Eqs. 1.28, 1.30); Wien peak from
  `d e_lambda,b / d lambda = 0`. The constants h, c and k are the exact SI values
  listed by NIST CODATA 2022; sigma is computed as `2 pi^5 k^4 / (15 h^3 c^2)`. sigma
  is registered in `core.constants` as `STEFAN_BOLTZMANN`, derived there from h, c, k.
- Band fractions `F(0 -> lambda T)`. With `z = c2/(lambda T)`,
  `F = (15/pi^4) int_z^inf x^3/(e^x - 1) dx`. For `z >= 1` the integral is
  `sum_n e^{-nz} (z^3/n + 3z^2/n^2 + 6z/n^3 + 6/n^4)` (geometric series of
  `1/(e^x - 1)` with DLMF 8.2.2, 8.4.8, 8.4.11; normalisation `pi^4/15 =
  Gamma(4) zeta(4)`, DLMF 25.5.1, 25.6.1). For `z < 1` the complement `int_0^z` is
  integrated term by term from the Bernoulli generating function (DLMF 24.2.1,
  B_0 to B_16 from Table 24.2.1). This derivation is the module's own; the tests check
  it against direct quadrature of Planck's law.
- View factors: Lienhard Table 10.2 (2-D strips, hinged strips, perpendicular strips,
  triangular duct by the crossed-string form, strip to cylinder, parallel cylinders,
  concentric cylinders), Table 10.3 (parallel rectangles = Howell C-11, perpendicular
  rectangles = C-14, coaxial disks = C-41, concentric spheres), Howell C-78 (cylinder
  interior to itself) and C-111 (frustum interior to base). Reciprocity, summation and
  a residual check for a full matrix (Eqs. 10.12, 10.15).
- Frustum (conical nozzle extension) enclosure: the openings' factor from C-41, the
  rest by summation and reciprocity, as in Lienhard Examples 10.3-10.4. The tests
  confirm the wall-to-end factors against the independent C-111 closed form and the
  cylinder limit against C-78.
- Gray diffuse surfaces: two-surface circuit (Eq. 10.23), parallel-plate and
  enclosed-body transfer factors (Eqs. 10.25, 10.27, 10.30), nested radiation shields
  (Eq. 10.31, generalised to several shields with different face emittances by adding
  the same series resistances) with shield temperatures, and the N-surface
  net-radiation method (Eqs. 10.39-10.44) with prescribed temperature or prescribed
  net flux (including adiabatic, re-radiating walls) on each surface.
- Gas radiation: mean beam length `L0 = 4V/A` (or 3.6V/A for strongly absorbing gas)
  and Table 10.4; weighted-sum-of-gray-gases (WSGG) total emissivity of H2O, CO2 and
  soot and their mixtures; total absorptivity by Lienhard Eq. 10.57,
  `alpha_g = eps_g(p_a L T_w/T_g, T_w) (T_g/T_w)^(1/2)`; gas-to-black-wall net flux
  `eps_g sigma T_g^4 - alpha_g sigma T_w^4` (Example 10.11); a bridge that reads a
  `RocketStation` (temperature, pressure, mole fractions) from
  `thermochemistry.rocket_performance`.
- Radiative-equilibrium wall temperature: the unique root of
  `q_imposed + h (T_r - T_w) = sigma sum_k F_k (T_w^4 - T_k^4)` (Lienhard
  Examples 1.5-1.6 and the transfer-factor form Eq. 10.2), any number of sinks (for
  example the inner and outer faces of a radiation-cooled nozzle extension).
- Interior of a gray, isothermal conical extension radiating through its two openings:
  transfer factors `e F_w-k / (1 - (1 - e) F_w-w)` (derived from Eqs. 10.39 and
  10.42 with black openings and checked against the N-surface solver).

**WSGG model and its limits**

- H2O and CO2: Bordbar, Fraga and Hostikka (2020) Tables 2-3, four gray gases plus a
  clear gas, `a_i = sum_j b_i,j (T/1200 K)^(j-1)`, `a_0 = 1 - sum a_i`. Fitted to
  line-by-line (HITEMP 2010) emissivities of the pure gas (mole fraction 1) at 1 atm
  total pressure for 300-2400 K and path lengths 0.01-60 m, i.e. pressure path
  lengths 0.01-60 atm m.
- Soot: Sadeghi et al. (2021) Table B9, Rayleigh soot with the Chang-Charalampopoulos
  refractive index, four gray gases without a clear gas, `a_j = max(sum b_j,n
  (T/1400 K)^n, 0)`, 300-3000 K. The paper does not state the `f_v L` range of its
  fit; a warning says so. As printed, the weights sum to 0.939-1.010 over 300-3000 K
  (not exactly one).
- Mixtures: superposition (Sadeghi et al. Eqs. 12-14): every combination of one gray
  gas per species, `kappa = sum kappa_i p_i`, `a = prod a_i`, assuming uncorrelated
  spectra. Molar-ratio limits do not apply because the pure-species tables are used
  (the reason Bordbar et al. built them). Bordbar et al. Table 4 reports, for this
  double-integration use of their Tables 2-3 in four 1-D LBL benchmarks,
  domain-average errors of 2.7-7.3 % and maxima up to 21.5 % (heat flux, case 1).
- Fail-closed rules: temperatures outside the fitted range always raise; a total
  pressure other than 1 atm raises unless `allow_pressure_extrapolation=True`
  (recorded as a critical warning: no pressure-broadening correction is applied); a
  pressure path length outside 0.01-60 atm m raises unless
  `allow_path_extrapolation=True`. The pure-gas fits are applied to diluted absorbers
  through `p_a L`, which neglects the difference between self- and foreign broadening.
- Consequence for rocket engines: combustion chambers (several MPa, about
  3000-3600 K) are outside these fits and are refused. Nozzle-exit and plume states
  below 2400 K can be evaluated, with an explicit pressure flag when the pressure is
  not 1 atm.

**Not included (no implementable relation opened and verified here)**

- Alumina (Al2O3) and other particle radiation in solid-motor plumes.
- Pressure-dependent WSGG models (the Bordbar et al. 2021 IJHMT paper could not be
  opened; the repository copy is behind a bot challenge).
- The CO WSGG coefficients of Sadeghi et al. Table B8: as printed (checked against the
  PDF image), the weight `a_2` becomes negative above about 1300 K and reaches
  -2.8 at 3000 K; with the paper's clipping that band vanishes. This could not be
  reconciled with the paper's text, so CO is not offered.
- Gas radiation inside enclosures with gray walls, non-isothermal gases, scattering,
  and the overlap correction for Hottel/Alberti charts (Lienhard only describes it).

## Source notes and discrepancies found while reading

- The free textbook site now serves only version 6.00 (2024) of Lienhard's book; all
  locators refer to it. Lienhard's own numerical examples use sigma = 5.67e-8 (and
  5.6704e-8, T(K) = T(C) + 273 in Chapter 1); the tests rescale or state the
  deviation.
- Lienhard Problem 10.2 gives "[41 % at 10,000 K]" for the visible fraction. With the
  visible band of his Table 1.2 (0.4-0.7 um) the band fraction is 32.7 %; 0.38-0.76 um
  gives 39.6 %. The problem's band is not stated, so this answer is not used as a test.
- Lienhard Example 10.2 reads 0.245 from the Fig. 10.9 chart; the Table 10.3 closed
  form (Howell C-14) gives 0.2406, so `F1-2` is 0.0406 rather than 0.045.
- Example 10.11 values (0.084, 0.073) come from Alberti et al. charts (CO2 at
  `p_a -> 0` in N2) and are a different model from the pure-CO2 WSGG data used here.

## Verification (measured in `tests/test_heat_radiation.py`)

- Constants: sigma, Wien and c2 agree with the 10 printed NIST digits (< 1e-9
  relative); Wien gives 2897.77 um K as printed by Lienhard. Simpson quadrature of
  Planck's law reproduces `sigma T^4` to 1e-9 at 300, 1000 and 3500 K.
- Band fractions agree with independent quadrature to 1e-10 on both sides of the
  series switch (a 40-digit check during development gave < 6e-16). Lienhard's
  statements: 25.005 % below the peak ("three-fourths" above); 90.06 % between 6.3 and
  42 um at 300 K and 89.96 % between 0.33 and 2.2 um at 5772 K (book: 90 %).
- Example 10.1: slit factors 0.060667, 0.939333, shield-to-room 0.084422 (book
  0.06067, 0.9393, 0.08441 with rounded intermediates), 611.95 W/m (611).
- Example 10.3/10.4: `F` heater-to-opening 0.19224 (0.192), heater-to-shield 0.80776
  (0.808), 1686.5 W (1687), shield area 0.097148 m2 (0.09715), shield-to-heater
  0.065304 (0.0653).
- Example 10.8: unshielded 0.6238 W/m (0.624), shielded 0.3283 W/m (0.328), shield
  213.4 K (213).
- Example 10.9: (a) B1 612.66, B2 379.52 W/m2 (612.1, 379.5), `Q1-2` 46.63 W/m
  (46.53); (b) B 988.5/657.8/823.1 (987.7/657.4/822.6), `Q` 198.5 W/m (198), reflector
  347.1 K (347). Book uses sigma = 5.67e-8; deviations <= 0.2 %.
- Example 10.10: -1296.9 W/m with the book's sigma (book -1294, which rounds the view
  factors to two digits): -0.22 %.
- Examples 1.5/1.6: 28.40 C and 23.52 C (book 28.4 C and 23.5 C).
- Example 10.11 (cross-model, not validation): WSGG pure-CO2 `eps_g` 0.0848 vs chart
  0.084 (+1.0 %), `alpha_g` 0.1351 vs 0.12 (+12.6 %), net 27.4 kW/m vs 27 kW/m.
- Identities: frustum algebra = Howell C-111 to 1e-12 and C-78 in the cylinder limit;
  C-40 = C-41 at equal radii; matrix residuals < 1e-14; N-surface solver = Eq. 10.27
  to 1e-12; closed-form extension transfer factors = solver to 1e-12; enclosure
  energy residual < 1e-9 W; WSGG superposition of two species equals
  `eps_w + eps_c - eps_w eps_c` to 1e-14; H2O and CO2 weights nonnegative and summing
  to one over 300-2400 K.

## Usage sketch

```python
from rocket_propulsion.heat_transfer.radiation import (
    RadiativeSink, frustum_interior_transfer_factors, radiative_equilibrium_wall_temperature,
)

# Conical radiation-cooled extension: 0.3 m attachment radius, 0.9 m exit, 1.5 m long.
to_upstream, to_exit, _ = frustum_interior_transfer_factors(
    small_radius_m=0.3, large_radius_m=0.9, length_m=1.5, wall_emissivity=0.85)
wall = radiative_equilibrium_wall_temperature(
    sinks=[RadiativeSink(to_exit, 3.0),          # inner face through the exit to space
           RadiativeSink(to_upstream, 1500.0),   # inner face toward the hot upstream nozzle
           RadiativeSink(0.85, 3.0)],            # outer face to space (Eq. 10.30)
    heat_transfer_coefficient_w_m2_k=250.0, recovery_temperature_k=2400.0)
```

The sink temperatures and the convective input above are placeholders chosen by the
user, not values taken from a source.

## References

- **LIENHARD-AHTT-6.00** — J. H. Lienhard IV and J. H. Lienhard V, *A Heat Transfer
  Textbook*, 6th ed., v6.00, 2024: <https://ahtt.mit.edu/>,
  <https://ahtt.mit.edu/wp-content/uploads/2024/04/AHTTv600.pdf>
- **NIST-CODATA-2022** — CODATA 2022 recommended values:
  <https://physics.nist.gov/cuu/Constants/index.html>,
  <https://physics.nist.gov/cgi-bin/cuu/Value?sigma>
- **DLMF-25.5** — NIST DLMF: <https://dlmf.nist.gov/25.5>, <https://dlmf.nist.gov/25.6>,
  <https://dlmf.nist.gov/8.2>, <https://dlmf.nist.gov/8.4>, <https://dlmf.nist.gov/24.2>
- **HOWELL-CATALOG-C** — J. R. Howell, configuration-factor catalog, Section C:
  <https://www.thermalradiation.net/indexCat.html>,
  <https://www.thermalradiation.net/calc/sectionc/C-41.html>,
  <https://www.thermalradiation.net/calc/sectionc/C-111.html>
- **BORDBAR-FRAGA-HOSTIKKA-2020** — Int. Commun. Heat Mass Transfer 110 (2020) 104400:
  <https://doi.org/10.1016/j.icheatmasstransfer.2019.104400>,
  <https://acris.aalto.fi/ws/portalfiles/portal/39751900/ENG_Bordbar_et_al_An_extended_weighted_sum_of_gray_gases_model_ICHMT.pdf>
- **SADEGHI-2021-FSJ** — Fire Safety Journal 125 (2021) 103420:
  <https://doi.org/10.1016/j.firesaf.2021.103420>,
  <https://acris.aalto.fi/ws/portalfiles/portal/66772583/1_s2.0_S0379711221001612_main_1.pdf>
