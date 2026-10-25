# Real-fluid thermodynamics for propellants

Package: `rocket_propulsion.thermodynamics.real_fluids` (standard library only).
Tests: `tests/test_real_fluids.py` (about 7 s). Reference data:
`tests/reference_data/nist_real_fluids.json`. All sources were opened on 2026-10-25;
full records are in [REFERENCES.md](REFERENCES.md).

## What it does

* **Reference equations of state** explicit in the reduced Helmholtz energy
  `alpha(tau, delta) = alpha0 + alphar`, with polynomial, exponential and Gaussian
  bell-shaped residual terms, for parahydrogen, normal hydrogen, oxygen, nitrogen, methane
  and helium.
* **Properties** from `alpha` and its derivatives (Span et al. 2000, Eqs. 56-71):
  p, rho, Z, u, h, s, g, cv, cp, speed of sound, Joule-Thomson coefficient,
  `(dp/drho)_T`, `(dp/dT)_rho`, and departure functions (`h - h_ig`, `s - s_ig`,
  `ln phi`, `cv - cv_ig`).
* **Density solve from (T, p)** on the stable branch: below `T_c` the vapor pressure
  picks liquid or vapor and the saturated densities bound the root, so metastable roots are
  never returned; above `T_c` the isotherm is searched upward from the ideal-gas density.
* **Saturation** from the Maxwell criterion (equal pressure and Gibbs energy), Newton
  iteration on `ln(delta')` and `ln(delta'')`. The source paper's vapor-pressure
  ancillary equation seeds the pressure and, where the paper gives one, its saturated-liquid
  density ancillary seeds the liquid root (oxygen Eq. 10, nitrogen Eq. 5, methane Eq. 3.4).
  Without ancillaries (helium; liquid side of hydrogen) the isotherm is scanned and the
  last spinodal loop is used. The result is checked to lie on the stable branches.
  `saturation_p` inverts `p_sat(T)` by regula falsi.
* **Cubic fallback**: Peng-Robinson and Soave-Redlich-Kwong, converted to residual
  Helmholtz form (Bell and Jäger 2016), so the same solvers and departure functions apply.
  Every cubic result carries the warning `cubic_eos_lower_accuracy`.

```python
from rocket_propulsion.thermodynamics.real_fluids import reference_fluid, cubic_counterpart

lox = reference_fluid("oxygen")
state = lox.state_tp(90.0, 1.0e6)          # FluidState, phase "liquid"
sat = lox.saturation_t(120.0)              # SaturationState with liquid and vapor
pr = cubic_counterpart("methane", "peng_robinson")
```

## Models, sources and validity

States outside the stated range raise `DomainError`; so do two-phase `(T, rho)` inputs, a
`(T, p)` exactly on the vapor-pressure curve, and saturation at `T >= T_c`. Saturation
within about `1e-5 T_c` of the critical point may raise `ConvergenceError`.

| Key | Equation (coefficients read from) | T range (K) | p max (MPa) |
|:---|:---|:---|---:|
| `parahydrogen` | Leachman et al. 2009, Tables 4-7 | 13.8033-1000 | 2000 |
| `normal_hydrogen` | Leachman et al. 2009, Tables 4-7 | 13.957-1000 | 2000 |
| `oxygen` | Schmidt and Wagner 1985, as printed in Stewart et al. 1991, Table 4 and Eq. (11) | 54.361-300 | 80 |
| `nitrogen` | Span et al. 2000, Tables 17-18 and Eq. (53) | 63.151-1000 | 2200 |
| `methane` | Setzmann and Wagner 1991, Tables 33 and 35 | 90.6941-625 | 1000 |
| `helium` | Ortiz Vega et al. 2023 (NIST IR 8474), Table 2 and Eq. (3) | 2.1768-1500 | 2000 |

Provenance and checks:

* The born-digital papers (hydrogen, nitrogen, helium) were read from the PDF text layer;
  the scanned papers (oxygen, methane) by comparing the OCR text with the rendered pages.
* Every residual coefficient with three or more significant digits was checked
  digit-for-digit against the CoolProp fluid files at commit `c00300e`; all matched.
  CoolProp served only as a second transcription; the numbers are those of the papers.
* The oxygen ideal-gas part is the analytic Helmholtz integral of Stewart et al. Eq. (11)
  (including the two-level electronic term), referenced to the ideal-gas enthalpy 8682 J/mol
  and entropy 205.037 J/(mol K) at 298.15 K and 0.101325 MPa (p. 919).
* Each fluid uses the gas constant printed in its paper (8.314472 for hydrogen and helium,
  8.31434 for oxygen, 8.31451 for nitrogen and methane).
* The hydrogen molar mass, 2.01588 g/mol, is not printed by Leachman et al.; it is the NIST
  WebBook species value.
* Setzmann and Wagner print `p_c = 4.5922 MPa` on p. 1064 but 4.5992 MPa in Eq. (2.3),
  Table 26 and p. 1091; 4.5992 MPa is used.
* Not modelled: the melting line (compressed-liquid results carry the info warning
  `melting_line_not_checked`), the helium lambda line, orthohydrogen, and transport
  properties.

## Measured agreement with NIST WebBook

The WebBook tables (SRD 69, `fluid.cgi`, 12 digits, fetched 2026-10-25, URL stored per
row) are computed from the same reference equations. The set has 106 single-phase states
covering liquid, vapor and supercritical, and 36 saturation pairs. Enthalpy and entropy are
compared as differences from the first state of each fluid, because the reference states
differ.

**Every relative deviation is below 11 ppm; for oxygen, nitrogen and methane the density and
saturation values agree to 0.04 ppm or better.** Values are maximum absolute deviations in
ppm (1e-6) relative, except: `h diff` is in ppm of `cp T`, `s diff` in ppm of `cp`, and
JT in 1e-6 K/MPa absolute. The Points column gives single-phase + saturation states.

| Fluid | Points | rho | cp | cv | w | h diff | s diff | p_sat | rho' | rho'' | dh_vap | JT |
|:---|:---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| Parahydrogen | 17 + 6 | 2.14 | 2.53 | 1.56 | 0.91 | 2.93 | 6.88 | 1.47 | 0.00 | 0.00 | 1.47 | 19.71 |
| Normal hydrogen | 17 + 6 | 2.16 | 2.56 | 1.57 | 0.92 | 2.93 | 6.55 | 1.47 | 0.00 | 0.00 | 1.47 | 19.97 |
| Oxygen | 17 + 6 | 0.00 | 3.15 | 4.52 | 1.02 | 0.43 | 0.54 | 0.00 | 0.00 | 0.00 | 0.00 | 58.29 |
| Nitrogen | 20 + 6 | 0.00 | 0.01 | 0.01 | 0.00 | 0.00 | 0.01 | 0.00 | 0.00 | 0.00 | 0.00 | 0.01 |
| Methane | 18 + 6 | 0.04 | 0.04 | 0.01 | 0.02 | 0.01 | 0.00 | 0.02 | 0.02 | 0.02 | 0.00 | 0.39 |
| Helium | 17 + 6 | 2.07 | 2.14 | 1.65 | 0.95 | 1.89 | 10.30 | 1.47 | 0.00 | 0.00 | 1.47 | 17.55 |

0.00 means below 0.005 ppm.

What explains the residuals:

* **Hydrogen and helium (about 1.5-2.6 ppm).** The WebBook evaluates these equations with
  R = 8.3144598 J/(mol K), whereas the papers print 8.314472. With the WebBook constant the
  density and sound speed agree to better than 1e-10, which a test checks. The
  papers' constant is kept.
* **Helium check values.** NIST IR 8474 Table 3 is reproduced to every printed digit with
  the SI-2019 gas constant, and to 1.1e-6 relative with the constant printed in its own
  Table 1. Table 4 (saturation) is reproduced to half a unit in the last printed digit.
* **Oxygen cp and cv (up to 4.5 ppm).** Only the printed ideal-gas heat capacity Eq. (11)
  is implemented; the WebBook ideal-gas representation differs slightly. Densities and
  saturation agree to 0.003 ppm.

Test tolerances are 10 ppm for density, cp, cv, w, p_sat, saturated densities and the
latent heat; 1e-4 K/MPa + 10 ppm for JT; 10 ppm of `cp T` for enthalpy differences; and
20 ppm of `cp` for entropy differences.

Thermodynamic-consistency tests, using central finite differences, check:

* `cp - cv = T (dp/dT)_rho^2 / (rho^2 (dp/drho)_T)`;
* the Maxwell relation `(ds/drho)_T = -(dp/dT)_rho / rho^2`;
* `cp = (dh/dT)_p`, `cv = (du/dT)_rho` and `mu_JT = -(dh/dp)_T / cp`;
* `g = h - T s`;
* equal Gibbs energy and pressure of the saturated phases (to 1e-10 RT/M);
* Clausius-Clapeyron, `dp_sat/dT = dh_vap / (T dv)`, to 1e-6;
* the `saturation_p` round trip;
* the departure enthalpy against `h(T, p) - h(T, 1 Pa)`;
* `ln phi` against the integral of `(Z - 1)/p dp`.

## Cubic fallback accuracy

`cubic_counterpart(key, kind)` uses the reference fluid's `T_c`, `p_c`, gas constant,
molar mass, ideal-gas part and validity range, with the WebBook acentric factor (hydrogen
-0.219, oxygen 0.0222, nitrogen 0.0372, methane 0.01142, helium -0.3836). The constants
(`Omega_a`, `Omega_b`, `m(omega)`, `Delta1`, `Delta2`) are those printed in Bell and Jäger
2016, Table 1. The Peng-Robinson 1976 and Soave 1972 papers were not opened (paywalled).
`cubic_model(...)` builds a cubic fluid from any sourced `T_c`, `p_c` and acentric factor,
together with an explicit validity range. It refuses temperatures where
`1 + m (1 - sqrt(T/Tc)) <= 0`.

**Cubic equations reproduce low-pressure vapor density to within about 6 %, but miss liquid
density by up to 41 % and cp by up to 157 %; use them only where no reference equation
exists.** Values are maximum absolute deviations in percent from the same WebBook states.

| EOS | Fluid | rho vapor <= 1 MPa | rho liquid | rho other | p_sat | cp | w |
|:---|:---|---:|---:|---:|---:|---:|---:|
| PR | Parahydrogen | 5.3 | 29.6 | 12.9 | 15.6 | 71.6 | 43.3 |
| PR | Normal hydrogen | 5.3 | 29.9 | 8.1 | 15.8 | 69.8 | 42.9 |
| PR | Oxygen | 1.0 | 12.8 | 11.9 | 9.1 | 24.0 | 19.5 |
| PR | Nitrogen | 1.4 | 13.3 | 11.6 | 3.8 | 15.0 | 72.8 |
| PR | Methane | 0.1 | 12.8 | 11.0 | 1.2 | 22.9 | 19.9 |
| PR | Helium | 5.7 | 40.9 | 16.6 | 9.5 | 131.0 | 82.2 |
| SRK | Parahydrogen | 1.6 | 15.4 | 2.3 | 30.7 | 90.2 | 47.1 |
| SRK | Normal hydrogen | 1.6 | 15.7 | 2.2 | 31.1 | 88.3 | 46.6 |
| SRK | Oxygen | 0.1 | 10.2 | 5.5 | 9.8 | 20.7 | 23.7 |
| SRK | Nitrogen | 0.1 | 6.6 | 12.0 | 2.9 | 11.6 | 86.2 |
| SRK | Methane | 0.0 | 11.8 | 8.7 | 6.6 | 20.9 | 19.4 |
| SRK | Helium | 1.6 | 26.3 | 24.5 | 28.4 | 156.7 | 95.0 |

The largest cp and w errors occur at near-critical and very-high-pressure states (for
example nitrogen at 1000 MPa). The quantum fluids, hydrogen and helium, are poorly
described by any corresponding-states cubic. The cubic tests check:

* the analytic derivatives against finite differences;
* the pressure against the pressure-explicit Eq. (9) of Bell and Jäger, to 1e-12;
* the critical point (inflection at `T_c` with `p = p_c` to 0.2 %, limited by the 5-digit
  `Omega` constants);
* nitrogen vapor density within 2 % and vapor pressure within 5 % of the WebBook.

## References

* LEACHMAN-2009-H2 — Leachman, Jacobsen, Penoncello and Lemmon, J. Phys. Chem. Ref. Data
  38, 721 (2009): <https://srd.nist.gov/jpcrdreprint/1.3160306.pdf>
* STEWART-1991-O2 — Stewart, Jacobsen and Wagner, J. Phys. Chem. Ref. Data 20, 917
  (1991): <https://srd.nist.gov/jpcrdreprint/1.555897.pdf>
* SCHMIDT-WAGNER-1985-O2 — Schmidt and Wagner, Fluid Phase Equilib. 19, 175 (1985), not
  opened: <https://doi.org/10.1016/0378-3812(85)87016-3>
* SPAN-2000-N2 — Span et al., J. Phys. Chem. Ref. Data 29, 1361 (2000):
  <https://srd.nist.gov/jpcrdreprint/1.1349047.pdf>
* SETZMANN-WAGNER-1991-CH4 — Setzmann and Wagner, J. Phys. Chem. Ref. Data 20, 1061
  (1991): <https://srd.nist.gov/jpcrdreprint/1.555898.pdf>
* ORTIZVEGA-2023-HE — Ortiz Vega et al., NIST IR 8474 (2023):
  <https://doi.org/10.6028/NIST.IR.8474>
* BELL-JAGER-2016-CUBIC — Bell and Jäger, J. Res. NIST 121, 238 (2016):
  <https://doi.org/10.6028/jres.121.011>
* PENG-ROBINSON-1976 — not opened: <https://doi.org/10.1021/i160057a011>
* SOAVE-1972 — not opened: <https://doi.org/10.1016/0009-2509(72)80096-4>
* COOLPROP-FLUIDS-C00300E — transcription cross-check only:
  <https://github.com/CoolProp/CoolProp/tree/c00300e97ad4cdb30032635d7379209d5e45a121/dev/fluids>
* NIST-WEBBOOK-FLUID-RF-2026 — NIST Chemistry WebBook SRD 69 fluid data:
  <https://webbook.nist.gov/chemistry/fluid/>, <https://webbook.nist.gov/cgi/fluid.cgi>
* NIST-WEBBOOK-H2-SPECIES — hydrogen molar mass:
  <https://webbook.nist.gov/cgi/cbook.cgi?ID=C1333740&Units=SI>
