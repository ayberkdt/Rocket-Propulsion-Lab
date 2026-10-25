# Convective heat transfer (thrust-chamber cooling focus)

`rocket_propulsion.heat_transfer.convection` collects convection correlations
read equation by equation from the sources listed under References. It is
stdlib-only and SI, uses Darcy friction factors and kelvin temperature ratios,
and fails closed: outside a source's stated range a function raises
`DomainError` with code `outside_correlation_range`; with `extrapolate=True` it
returns the value with a `correlation_extrapolated` warning. Correlations with
no numerical range in the sources read carry a `range_not_stated` info warning.
Most functions return a frozen `CorrelationResult(value, quantity, correlation,
source, warnings)`. Lienhard page numbers are printed pages of the 6th edition
(v6.00); in the PDF opened, printed = PDF - 14 in Chapters 6-9.

Reference values with locators: `tests/reference_data/heat_convection.json`;
tests: `tests/test_heat_convection.py` (35 tests, about 3 s).

## Public API

### Groups, geometry, friction

| Function | Equation |
|---|---|
| `reynolds_number`, `mass_flux_reynolds_number`, `prandtl_number` (re-export), `peclet_number`, `stanton_number`, `nusselt_number`, `heat_transfer_coefficient` | definitions (Lienhard Ch. 6-7; Taylor symbols p. 3) |
| `graetz_number` | Lienhard Eq. (7.26) |
| `rayleigh_number` | Lienhard Example 8.1 |
| `hydraulic_diameter`, `rectangular_duct_hydraulic_diameter` | Eqs. (7.60), (7.61a) |
| `roughness_reynolds_number`, `roughness_regime` | Eqs. (7.48a)-(7.48d) |
| `laminar_darcy_friction_factor` | `64/Re`, Fig. 7.6 (Re <= 2300) |
| `filonenko_friction_factor` | Eq. (7.42) (range of Eq. 7.41 enforced) |
| `haaland_friction_factor` | Eq. (7.50), 4000 <= Re <= 1e8, 0 <= eps/D <= 0.05 |
| `chen_friction_factor` | explicit Colebrook used by NASA RTE, Wadel Eq. (12) |
| `ito_curvature_friction_factor` | `[Re (r/R)^2]^0.05`, CR-678 Eq. (4) = Wadel Eq. (13) |

### Internal flow

| Function | Equation and range |
|---|---|
| `laminar_fully_developed_nusselt` | Table 7.4 (circle, square, 2:1, 4:1, 8:1, plates), Eq. (7.63) one side insulated |
| `laminar_thermal_entry_length` | Eq. (7.24) |
| `laminar_entry_nusselt` | Eqs. (7.28) local / (7.29) average (isothermal), (7.30) local (uniform flux) |
| `gnielinski_nusselt` | Eq. (7.41), 2300 <= Re <= 5e6, 0.6 <= Pr <= 1e5 |
| `smooth_pipe_analogy_nusselt` | Eq. (7.35) (Petukhov-type), Pr >= 0.7 |
| `power_law_nusselt` | Colburn (7.39a, 0.023 Pr^1/3), Dittus-Boelter as printed by Lienhard (7.39b, 0.0243 Pr^0.4), McAdams (7.39c); Re >= 1e4, 0.67 <= Pr <= 100 |
| `sieder_tate_nusselt` | Eq. (7.40) |
| `gnielinski_power_law_nusselt` | Eqs. (7.43a)-(7.43b) |
| `property_ratio_nusselt_factor`, `property_ratio_friction_factor` | Eqs. (7.44)-(7.47) with their ratio ranges |
| `turbulent_entry_nusselt_ratio` | Eq. (7.31) and Table 7.2 (average Nu, Pr ~ 0.7, L/D > 3) |
| `rough_pipe_nusselt` | Bhatti-Shah Eq. (7.49), Re >= 1e4, 0.5 <= Pr <= 10, 0.002 <= eps/D <= 0.05, Re_eps > 70 |
| `liquid_metal_nusselt` | Skupinski Eq. (7.55) (uniform flux, Pe > 100), Seban-Shimazaki Eq. (7.54) |

### Rocket coolant side

| Function | Source |
|---|---|
| `taylor_temperature_ratio_exponent`, `taylor_hydrogen_nusselt` | Taylor TN D-4332 Eqs. (19)-(20); range Table I: 2 <= x/D <= 252, 1.1 <= Ts/Tb <= 27.6, 7500 <= Re_b <= 1.38e7 |
| `taylor_near_critical_region`, `pseudocritical_temperature` | Taylor p. 14 region 2 (25 K < T_i < T*, p_c < p_i < 3.65 MPa); T* = maximum of c_p from the package's Helmholtz EOS |
| `hess_kunz_nusselt`, `hess_kunz_cl_factor` | Hess-Kunz as reported by Taylor Eq. (14) and CR-678 Eq. (3) with the C_L table (p. 3; linear interpolation is this module's choice) |
| `sp125_sieder_tate_nusselt` | SP-125 Eq. (4-23), coolant-specific C1 supplied by the caller |
| `sp125_supercritical_film_coefficient` | SP-125 Eq. (4-25) as printed (dimensionally homogeneous; Nusselt form `0.029 Re^0.8 Pr^(1/3) (Tco/Twc)^0.55`) |
| `curvature_heat_transfer_factor` | RTE curvature multiplier `[Re (r/R)^2]^(+/-0.05)`, Wadel Eq. (11); SP-8087 p. 61 guidance (no factor for liquids) emitted as a warning |
| `rte_coolant_nusselt` | RTE coolant equation `C_c Re^0.8 Pr^0.4 phi_ent phi_cur`, Wadel Eq. (6), C_c = 0.023 as in the TM |
| `taylor_hydrogen_coefficient`, `hess_kunz_coefficient`, `sp125_supercritical_film_heat_transfer` | the above evaluated with real-fluid density and c_p (`thermodynamics.real_fluids`), returning `CoolantHeatTransfer` |
| `burnout_margin_screen` | SP-8087 p. 61: operate below 80 % of the mean burnout heat flux |

Viscosity and thermal conductivity are always caller inputs: the package has no
verified transport-property model for cryogenic coolants, and none was added.

### External, natural convection, phase change, film cooling

| Function | Source |
|---|---|
| `flat_plate_laminar_nusselt` | Eqs. (6.58)/(6.68), Pr >= 0.6 |
| `flat_plate_turbulent_nusselt` | Eqs. (6.102) + (6.111), Pr >= 0.5 |
| `flat_plate_turbulent_power_law_nusselt` | Eqs. (6.112) gases, (6.113) liquids (Pr <= 85) |
| `churchill_bernstein_cylinder_nusselt` | Eqs. (7.65), (7.66) Re < 4000, (7.68) 4e4 < Re < 4e5; Re Pr >= 0.2 |
| `churchill_chu_vertical_plate_nusselt` | Eqs. (8.13a) (Ra <= 1e9), (8.13b) all Ra |
| `rohsenow_heat_flux`, `rohsenow_wall_superheat` | Eq. (9.4) (typical errors 100 % in q, 25 % in dT) |
| `taylor_wavelength`, `peak_pool_boiling_heat_flux` | Eqs. (9.6a), (9.11) (C = 0.149, heater > ~3 lambda_d1 checked when given), (9.12) Zuber (C = 0.131) |
| `corrected_latent_heat`, `film_condensation_vertical_plate_nusselt` | Eqs. (8.61), (8.62b) |
| `hatch_papell_velocity_ratio_function`, `hatch_papell_effectiveness`, `hatch_papell_required_coolant_flow`, `film_cooled_wall_temperature` | Hatch-Papell TN D-130 Eqs. (10)-(13); data ranges of the summary enforced (slot 0.063-0.50 in, V_g 104-1040 ft/s, V_g/V_c 0.45-33.3, temperatures when given); warning for eta < 0.2 |

## Measured agreement with the sources

All worked examples reproduce the printed numbers to the printed digits:

| Check | Printed | Computed |
|---|---|---|
| Lienhard Ex. 7.2, Eq. (7.30), Gz = 43.56 | 5.06 | 5.057 |
| Ex. 7.3, Filonenko f / Gnielinski x (1.75)^0.11 / h | 0.0136 / 1570 / 8400 W/m²K | 0.01358 / 1569.7 / 8398 |
| Ex. 7.4, Haaland f / Re_eps / Eq. (7.49) | 0.02651 / 71.2 / 2959 | 0.026513 / 71.21 / 2959.0 |
| Ex. 7.5, Skupinski, Pe = 4470 | 24.2 | 24.14 |
| Ex. 7.6, Eq. (7.29), Gz = 11.8 | 4.27 | 4.266 |
| Ex. 7.7, square duct, f / Nu | 0.02646 / 49.82 | 0.026458 / 49.818 |
| Ex. 7.8, Eq. (7.66), Re = 457 | 10.69 | 10.694 |
| Ex. 6.5, Eq. (6.68) | 366.2 | 366.15 |
| Ex. 8.1, Eq. (8.13a) | 47.33 | 47.34 |
| Ex. 8.6, h'_fg / Nu_L (book via rounded delta) | 2281 kJ/kg / 3846 | 2280.6 / 3847.3 |
| Ex. 9.2, Rohsenow q at 22 K, Csf = 0.016 | ~800 kW/m² (read from figure) | 815 kW/m² (+1.9 %) |
| Ex. 9.5 / 9.6, Eq. (9.11) | 1.260 / 1.328 MW/m² | 1.2599 / 1.3282 |
| lambda_d1 for water (p. 488) | 2.72 cm | 2.726 cm |
| Taylor Fig. 2, C2 at x/D = 20, 40, 140 (read-off, +/-0.03) | 0.49, 0.53, 0.56 | 0.4905, 0.5303, 0.5586 |
| Hatch-Papell example: ln eta / eta | -0.326 / 0.722 | -0.3273 / 0.7208 |
| Hatch-Papell required flow, iterations 1 and 2 | 16.84 / 17.56 lb/s | within 0.5 % |

Notes: the Hatch-Papell example rounds (S V_g/alpha)^0.125 to 4.51 and f to 1.31
(exact 4.519, 1.314), and its main-gas velocity (2000 ft/s) lies above the
1040 ft/s top of its own data range, so the test runs with `extrapolate=True`.
The printed effectiveness 0.722 differs from (1472 - 800)/(1472 - 540) = 0.7210.

Cross-checks in overlap regions (measured here, used as regression guards):
Chen (RTE) vs Haaland friction within 1.6 % for 1e4 <= Re <= 1e7,
1e-5 <= eps/D <= 0.05; Gnielinski power laws (7.43a/b) vs Eq. (7.41) within
26 % over 2300 <= Re <= 5e6, 0.7 <= Pr <= 500 (largest at Pr = 1.5 near the
fit boundary and at Re = 2300, Pr = 500); Eq. (7.41) / Eq. (7.35) equals
(Re - 1000)/Re exactly; turbulent flat plate Eq. (6.111) vs Eq. (6.112) for
air within 4.9 % for 1e5 <= Re_x <= 1e7; Taylor's equation at C2 = 0 equals
the RTE form with C_c = 0.023; the RTE concave curvature factor equals Ito's
friction ratio for the same Re and r/R.

## Not implemented or not verified

- Hess and Kunz's original paper and validity range were not opened; the
  equation is used as reported by Taylor and CR-678.
- The RTE code itself (Naraghi; NTRS record 19950020934 has no full text) was
  not read; its coolant model is taken from Wadel's TM. Wadel's entrance factor
  Eq. (10) is printed with an undefined `+1` in the numerator and is not
  implemented (`rte_coolant_nusselt` takes `entrance_factor` as an input).
- Roughness enhancement for rocket coolant channels (Norris,
  Dipprey-Sabersky as reported by NASA) was not found in a readable NASA
  source; only Lienhard's Bhatti-Shah fully rough correlation is provided.
- Transpiration cooling: no source opened; not implemented. SP-8087 Table XI
  (liquid-film length and gaseous-film integration) was read but not
  implemented.
- Forced-convection burnout (CHF) correlations for propellants (SP-8087 ref.
  34) and SP-125 Eq. (4-24) (unspecified coolant constant C2) are not
  implemented; `burnout_margin_screen` takes the burnout flux as an input.
- The Miller-Seader-Trebes equation (CR-678 Eq. 2) is not implemented: its
  reference temperature is defined only in Taylor's symbol list.
- No NASA data points from CR-678 or Taylor's tables were digitised; the
  Taylor checks are the Fig. 2 exponent read-off and limits.
- Pseudocritical temperatures come from the package's Helmholtz equations; they
  were not compared with NIST tables in this work package.
- No hardware validation. The coupled gas/wall/coolant solve that uses these
  correlations is `heat_transfer.regen_chamber` ([REGENERATIVE_COOLING.md](REGENERATIVE_COOLING.md)).

## References

Full records are in [REFERENCES.md](REFERENCES.md):

- LIENHARD-AHTT-6.00 — J. H. Lienhard V and J. H. Lienhard IV, *A Heat Transfer Textbook*, 6th ed. (v6.00), Phlogiston Press, 2024. <https://ahtt.mit.edu/>, <https://ahtt.mit.edu/wp-content/uploads/2024/04/AHTTv600.pdf>
- NASA-TN-D-4332-TAYLOR-1968 — M. F. Taylor, NASA TN D-4332, 1968. <https://ntrs.nasa.gov/citations/19680004771>
- NASA-CR-678-AEROJET-1967 — Aerojet-General Corp., NASA CR-678, 1967. <https://ntrs.nasa.gov/citations/19670008827>
- NASA-TM-1998-206313-WADEL — M. F. Wadel, NASA/TM-1998-206313, 1998. <https://ntrs.nasa.gov/citations/19980017619>
- NASA-SP125-COOLANT-SIDE — Huzel and Huang, NASA SP-125, 2nd ed., 1971, pp. 105-106. <https://ntrs.nasa.gov/api/citations/19710019929/downloads/19710019929.pdf>
- NASA-SP-8087-1972 — NASA SP-8087, Liquid Rocket Engine Fluid-Cooled Combustion Chambers, 1972. <https://ntrs.nasa.gov/citations/19730022965>
- NASA-TN-D-130-HATCH-PAPELL-1959 — Hatch and Papell, NASA TN D-130, 1959. <https://ntrs.nasa.gov/citations/19890068390>
- LEACHMAN-2009-H2 and the other real-fluid equations of state are listed in `REFERENCES.md`.
