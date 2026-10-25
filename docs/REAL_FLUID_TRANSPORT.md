# Real-fluid transport properties for propellants

Module: `rocket_propulsion.thermodynamics.real_fluids.transport` (standard library only).
Tests: `tests/test_real_fluid_transport.py` (about 3 s). Reference data:
`tests/reference_data/nist_real_fluid_transport.json`. All sources were opened on 2026-10-25,
except the two papers marked "not opened" below; full records are in
[REFERENCES.md](REFERENCES.md).

## What it does

Viscosity `mu`, thermal conductivity `k` and Prandtl number `Pr = cp mu / k` of cryogenic and
supercritical propellants and pressurants as functions of temperature and density, using the
correlations the NIST Chemistry WebBook uses. The density, `cp`, `cv`, `(drho/dp)_T` and
`(dp/dT)_rho` come from the reference equations of state of `real_fluids` (see
`REAL_FLUIDS.md`), which are the equations the transport correlations were fitted with.

```python
from rocket_propulsion.thermodynamics.real_fluids.transport import transport_tp

lox = transport_tp("oxygen", 90.0, 1.0e6)          # mu, k, Pr and their parts
lh2 = transport_tp("parahydrogen", 25.0, 2.0e6)    # no caller viscosity needed
he = transport_tp("helium", 300.0, 20.0e6)
```

Public API: `transport_tp(fluid, T, p)`, `transport_rho_t(fluid, rho, T)`,
`transport_from_state(fluid, FluidState)`, `viscosity_rho_t(fluid, rho, T)`,
`thermal_conductivity_rho_t(fluid, rho, T)`, `dilute_gas_transport(fluid, T)`,
`transport_validity(fluid)`, `TRANSPORT_FLUID_KEYS`, and the frozen dataclasses
`TransportProperties`, `TransportValidity`, `OlchowySengersParameters`. Every call accepts an
optional `viscosity_pa_s` override, which replaces the viscosity correlation and is also used in
the conductivity critical enhancement. Results carry the dilute, residual and critical parts and
`WarningMessage`s (`transport_near_critical` when the enhancement exceeds 10 % of `k`).

## Correlations, sources and validity

States outside the stated range raise `DomainError`; two-phase `(T, rho)` inputs are refused by
the equation of state. `transport_tp`, `transport_rho_t` and `transport_from_state` check both
ranges; `viscosity_rho_t` and `thermal_conductivity_rho_t` check only their own.

| Key | Viscosity (range) | Thermal conductivity (range) |
|:---|:---|:---|
| `nitrogen` | Lemmon and Jacobsen 2004, Eqs. (1)-(3) (63.151-1000 K, 2200 MPa) | Lemmon and Jacobsen 2004, Eqs. (4)-(11) (same) |
| `oxygen` | Lemmon and Jacobsen 2004, Eqs. (1)-(3) (54.361-300 K, 80 MPa) | Lemmon and Jacobsen 2004, Eqs. (4)-(11) (same) |
| `parahydrogen` | Muzny et al. 2013 via CoolProp, paper not opened (13.8033-1000 K, 200 MPa) | Assael et al. 2011, Table 5 (13.8033-1000 K, 100 MPa) |
| `normal_hydrogen` | Muzny et al. 2013 via CoolProp, paper not opened (13.957-1000 K, 200 MPa) | Assael et al. 2011, Table 2 (13.957-1000 K, 100 MPa) |
| `methane` | Friend, Ely and Ingham 1989, Eqs. (10a), (12), (15a) (90.6941-400 K, 55 MPa) | Friend, Ely and Ingham 1989, Eqs. (13a)-(22) (90.6941-625 K, 100 MPa) |
| `helium` | Arp, McCarty and Friend 1998 (NIST TN 1334), Eqs. (1)-(9) (2.5-1500 K, 100 MPa) | Hands and Arp 1981 via CoolProp, paper not opened (3.5-1500 K, 100 MPa) |

How the ranges were set: Lemmon and Jacobsen state validity "over all liquid and vapor states"
(abstract), so the equation-of-state range is used, with a warning above 100 MPa for nitrogen
(uncertainties unknown there, p. 62). Assael et al. give the range in the title; the Muzny et
al. range is the one in the NIST publication record. Friend et al. give "about 91 to 400 K for
pressures < 55 MPa" (viscosity) and "about 91 to 700 K for pressures below 100 MPa"
(conductivity) in the abstract; the conductivity is capped at 625 K by the Setzmann-Wagner
equation of state, and above the viscosity range its critical term uses the extrapolated
viscosity (info warning). TN 1334 gives 2.5-1500 K to 100 MPa for its helium I code (Sec. 1)
and an unknown conductivity validity below 3.5 K or above 100 MPa (Sec. 4.1).

Critical enhancement: the simplified Olchowy-Sengers model as printed by Lemmon and Jacobsen
(Eqs. 7-11, with their Boltzmann constant 1.380658e-23 J/K) and Assael et al. (Eqs. 4-7, with
the exact SI 2019 Boltzmann constant, which they do not print); for methane, Eqs. (18a)-(22) of
Friend et al. with the analytic equation of state everywhere (their scaled equation of state
for `|T*| < 0.03`, `|rho*| < 0.25` is not implemented; with the WebBook viscosity supplied,
four WebBook points fetched there were reproduced to 0.15 % without it, so only an info warning
is issued); for helium, the Hands-Arp term as transcribed by CoolProp (between 3.5 and 12 K).

## Provenance of the hydrogen viscosity and the helium conductivity

**The primary papers of these two correlations were not opened.** Muzny, Huber and Kazakov
(J. Chem. Eng. Data 58, 969, 2013) was not openly accessible (publisher page blocked; the NIST
publication record links only the DOI); Hands and Arp (Cryogenics 21, 697, 1981) is a
subscription article. Their coefficients are transcribed from CoolProp at commit
`c00300e97ad4cdb30032635d7379209d5e45a121` (`dev/fluids/Hydrogen.json`,
`dev/fluids/ParaHydrogen.json`, `src/Backends/Helmholtz/TransportRoutines.cpp`), and
**correctness rests on that transcription as verified against NIST WebBook rows**, which the
WebBook computes from the same correlations:

* Hydrogen viscosity: 16 WebBook rows per hydrogen form (liquid, vapor, supercritical to
  150 MPa, near-critical) agree to 4.6e-12. For hydrogen the WebBook uses the normal-hydrogen
  correlation for both forms, and so does this module. The erratum of the paper (opened) supplies
  the Avogadro constant 6.022137e23 /mol, the exponent `-i` of the initial-density sum and the
  scale density 90.909090909 kg/m3 (CoolProp agrees), and its three check values agree to
  6.1e-6 (they are printed to 5 significant digits).
* Helium conductivity: 15 WebBook rows (3.6-1500 K, to 90 MPa, liquid, vapor, supercritical,
  near-critical) agree to 9.0e-7. CoolProp notes that the critical term follows the REFPROP source
  rather than the paper, including two numerical factors whose origin NIST could not explain;
  they are reproduced as transcribed.

The helium viscosity is not in this category: TN 1334 was opened and Eqs. (1)-(9) were read from
the rendered pages (the OCR text layer is garbled in places). CoolProp's version truncates the
last coefficient of Eq. (2) (0.012520841 versus the printed 0.0125208416) and omits the linear
100-110 K blend of Sec. 4.2.2; with the printed values the WebBook is matched to 5.5e-10, while
CoolProp's version differs by up to 1.8 % between 100 and 110 K.

## Remaining gaps and deviations from the sources

* **Methane viscosity differs from the WebBook.** The WebBook uses an unpublished correlation
  (Quinones-Cisneros, Huber and Deiters, 2011). The published Friend-Ely-Ingham viscosity is
  implemented; it matches the paper's own Table A3 but differs from the WebBook by up to 4.4 %.
* **Methane conductivity, Eq. (16).** The reprint prints `delta*_sigma = 11` outside the
  saturated-vapor branch. With 11 the paper's own Table A3 is missed by 33 % in the liquid
  (100 K, 0.1 MPa: 134.8 instead of 199.8 mW/(m K)); with 1 the result is 199.7. The value 1 is
  used.
* **Argon and air**, also in Lemmon and Jacobsen 2004, are not included (no equation of state in
  `real_fluids`).

## Agreement with the NIST WebBook

76 WebBook rows (`Digits=12`), 9-16 per fluid, in the liquid, vapor and supercritical regions,
including 2-3 near-critical rows per fluid. One helium row (5.3 K, 0.235 MPa) was not used
because the WebBook returned "undefined" for `cp` and `k`.

Nitrogen, oxygen, hydrogen and helium match the WebBook to 2.2e-5 or better except in dense states
with a very small critical enhancement; methane differs at the percent level because of its
viscosity. Maximum |relative deviation|, evaluated at the WebBook `(T, rho)`:

| Fluid | Viscosity, all rows | k, vapor | k, supercritical | k, near-critical | k, dense, small enhancement |
|:---|---:|---:|---:|---:|---:|
| Nitrogen | 8.8e-10 | 2.6e-09 | 3.0e-07 | 2.3e-06 | 3.1e-04 |
| Oxygen | 8.6e-10 | 2.2e-09 | 2.0e-07 | 1.9e-06 | 4.7e-04 |
| Helium | 5.5e-10 | 2.7e-08 | 1.5e-07 | 3.3e-07 | 9.0e-07 |
| Normal hydrogen | 3.7e-12 | 3.1e-06 | 6.2e-06 | 5.8e-06 | 1.1e-03 |
| Parahydrogen | 4.6e-12 | 1.3e-05 | 2.2e-05 | 8.1e-06 | 5.8e-04 |
| Methane, own viscosity | 4.4e-02 | 6.9e-03 | 2.6e-03 | 1.8e-02 | 8.6e-04 |
| Methane, WebBook viscosity supplied | — | 7.1e-03 | 3.3e-03 | 3.2e-03 | 8.6e-04 |

Notes on the table:

* "Dense, small enhancement" are states with `rho > rho_c` whose enhancement is below 1 % of
  `k` (for methane, the liquid rows). There the WebBook differs from the Olchowy-Sengers model
  as printed by up to about 1.7e-4 W/(m K), in either direction, including states where the
  correlation-length bracket is negative and Lemmon and Jacobsen prescribe `lambda_c = 0`
  (p. 27). The papers do not explain this and the cause was not identified; the tests allow 2e-3
  there.
* Through `transport_tp` (density from the package equation of state) the hydrogen and helium
  deviations grow by the WebBook gas-constant effect described in `REAL_FLUIDS.md`: viscosity at
  most 2.5e-6, conductivity at most 2.2e-5.
* Methane: replacing the Friend-Ely-Ingham viscosity by the WebBook viscosity in Eq. (18a)
  removes most of the near-critical difference, and the remaining dilute-gas offset (300 K,
  0.1 MPa: +0.52 % in `k`) tracks the viscosity offset (+0.55 %). This suggests the WebBook
  uses its own viscosity in Eqs. (13a) and (18a); that is an inference, not documented by NIST.

## Checks against the papers' own values

* Lemmon and Jacobsen, Table V (6 digits): all 24 values within 2.4e-6, except the two
  near-critical rows (equation-of-state dependent), within 1.8e-5.
* Assael et al., Table 7 (5 digits): 12 values without critical enhancement within 4.6e-5.
  The 35 K, 30 kg/m3 row needs the background viscosity REFPROP used at the time, which is not
  printed; both fluids imply the same value (3.474 uPa s), which is a consistency check of the
  enhancement formula but is not used as a test.
* Muzny et al. erratum check values (40 K at 0, 50 and 100 kg/m3): within 6.1e-6.
* Friend et al., Table A1 (dilute gas, 3-4 digits): within 1.8e-4. Table A3 at its 4-digit
  densities: within 1.7e-3 (consistent with the density rounding).

## Dilute-gas limit against `thermochemistry/transport.py`

The zero-density values are compared with the NASA CEA fits over their overlap (200-1000 K,
oxygen to 300 K, methane to 625 K, helium to 1500 K). They are independent correlations of
overlapping data and should agree within their combined uncertainty; Lemmon and Jacobsen quote
2 % for the dilute-gas conductivity. Measured: viscosity within 1.0 % (normal hydrogen at
200 K), conductivity within 3.7 % (normal hydrogen at 1000 K, where Assael et al. quote 4-10 %
theory-versus-data differences above 600 K). Parahydrogen conductivity is not compared: its
rotational heat capacity differs from the normal hydrogen of the CEA data, and the
conductivities differ by 19 % at 200 K.

## References

* LEMMON-JACOBSEN-2004-TRANSPORT — Lemmon and Jacobsen, Int. J. Thermophys. 25, 21 (2004):
  <https://trc.nist.gov/refprop/FAQ/NAO.PDF>
* ASSAEL-2011-H2-CONDUCTIVITY — Assael et al., J. Phys. Chem. Ref. Data 40, 033101 (2011):
  <https://srd.nist.gov/jpcrdreprint/1.3606499.pdf>
* FRIEND-ELY-INGHAM-1989-CH4 — Friend, Ely and Ingham, J. Phys. Chem. Ref. Data 18, 583
  (1989): <https://srd.nist.gov/jpcrdreprint/1.555828.pdf>
* MUZNY-2022-H2-VISCOSITY-ERRATUM — erratum to Muzny, Huber and Kazakov (2013), opened; the
  2013 paper itself was not opened:
  <https://tsapps.nist.gov/publication/get_pdf.cfm?pub_id=934278>; NIST publication record
  (range): <https://www.nist.gov/node/574776>
* COOLPROP-TRANSPORT-C00300E — CoolProp commit c00300e, source of the hydrogen viscosity and
  helium conductivity coefficients:
  <https://github.com/CoolProp/CoolProp/blob/c00300e97ad4cdb30032635d7379209d5e45a121/dev/fluids/Hydrogen.json>,
  <https://github.com/CoolProp/CoolProp/blob/c00300e97ad4cdb30032635d7379209d5e45a121/dev/fluids/ParaHydrogen.json>,
  <https://github.com/CoolProp/CoolProp/blob/c00300e97ad4cdb30032635d7379209d5e45a121/src/Backends/Helmholtz/TransportRoutines.cpp>,
  <https://github.com/CoolProp/CoolProp/blob/c00300e97ad4cdb30032635d7379209d5e45a121/src/Backends/Helmholtz/HelmholtzEOSMixtureBackend.cpp>
* ARP-MCCARTY-FRIEND-1998-HE — NIST TN 1334 (revised), helium viscosity:
  <https://nvlpubs.nist.gov/nistpubs/Legacy/TN/nbstechnicalnote1334.pdf>
* NIST-WEBBOOK-FLUID-TRANSPORT-2026 — NIST Chemistry WebBook SRD 69 fluid data:
  <https://webbook.nist.gov/chemistry/fluid/>, <https://webbook.nist.gov/cgi/fluid.cgi>
* Equation-of-state sources (LEACHMAN-2009-H2, STEWART-1991-O2, SPAN-2000-N2,
  SETZMANN-WAGNER-1991-CH4, ORTIZVEGA-2023-HE) are listed in `REAL_FLUIDS.md` and
  `REFERENCES.md`.
