# Combustion instability

`rocket_propulsion.thermochemistry.combustion_instability` collects the
closed-form relations of NASA SP-194, *Liquid Propellant Rocket Combustion
Instability* (Harrje and Reardon, eds., 1972): chamber acoustic modes,
time-lag stability limits, and the sizing and measurement of damping. Every
equation was read from the scanned report; page numbers are in the
docstrings.

Tests: `tests/test_combustion_instability.py` (17 tests, reference values in
`tests/reference_data/combustion_instability.json`) and
`tests/test_feed_coupled_instability.py` (11 tests). HTTP route:
`POST /api/v1/combustion-stability`; workbench: **Yanma → Ürünler ve verim →
Kararlılık**.

## Public API

| Function | What it gives | SP-194 |
|---|---|---|
| `transverse_eigenvalue(nu, eta)` | `s_nu,eta`, root of `J_nu'(s) = 0` | Table 3.5.1, p. 143 |
| `acoustic_mode_frequency(...)` | `f = (a/2 pi) sqrt((s/r_c)^2 + (j pi/L)^2)`, with `(1 - M^2)` for longitudinal modes | Eq. (3.5.1-19), p. 144; p. 182 |
| `chamber_acoustic_modes(...)` | every T, R, L and combined mode below a frequency, sorted | same |
| `chugging_stability(...)` | gain `K`, unconditional stability, critical time lag, neutral frequency | Eq. (5.5.1-6), Example 5.5.1.1, p. 253 |
| `intrinsic_critical_time_lag(n, theta)` | Crocco-Cheng critical lag with rigid feeds (`n > 1/2`) | Eqs. (5.3.1-18)-(5.3.1-19), p. 239 |
| `neutral_interaction_index(omega_tau, G_R)` | `n = G_R/(1 - cos omega tau)`, minimum `G_R/2` | Eq. (4.2.2-23), Fig. 4.2.2a, p. 180 |
| `orifice_effective_length(...)` | `l_o + 0.85 d (1 - 0.7 sqrt(a))`, or the high-amplitude `0.375 delta` | Eqs. (8.3.2-1)-(8.3.2-2), pp. 410-411 |
| `helmholtz_resonator_frequency(...)` | `omega = a sqrt(A_o/(l_eff V))` | Eq. (8.3.1-2), p. 401 |
| `liner_resonant_frequency(...)` | flat-array form `a sqrt(a/(L_cav l_eff))` | Eq. (8.3.3-1), p. 415 |
| `response_factor(p, m)` | `N = sum(m'p')/sum(p'^2)`; positive drives the mode (Rayleigh) | p. 210, Eq. (4.4.1-5) |
| `quality_factor(f0, f1, f2)`, `damping_rate_from_decay(t, p)` | `Q = f0/(f2 - f1)`; decay rate from a ringing envelope | Sect. 9.6.6.2, pp. 523-524 |

## Feed-coupled chugging

`thermochemistry.feed_coupled_instability` couples the time-lag chamber to
each propellant's feed system through its injection admittance
`G = -(p̄/ṁ)(ṁ_i'/p')`, built from SP-194 Sect. 5.4: constant-area line
(Eq. 5.4.1-2) with wall-compliant sound speed (5.4.1-5), orifice with
inertance (5.4.1-9), manifold capacitance (5.4.1-16) and tank (Sect. 5.4.2).
`pressurized_tank_injection_admittance` is the closed form for a
tank-line-manifold-orifice feed; `feed_coupled_stability` applies the Nyquist
criterion of Sect. 5.5.1 to the single- or double-time-lag characteristic
equation (5.2.2-6 or 5.2.3-3), counting right-half-plane zeros from the
encirclements of the origin by `1 + f(iω)`.

The contour starts at zero frequency. A local derivative envelope bounds
each interval's image by a disk that excludes zero; unresolved intervals
are subdivided. The uncomputed tail is closed only when a global amplitude
envelope places all of it in `|f| < 1`, rather than testing one terminal
sample. The final phase increment includes that closure. `points` is now
the evaluation budget (1000–1000000); exhaustion or an unresolved boundary
root raises an error and supplies no stability verdict.

Physical `FeedSystem` inputs derive envelopes from their passive component
model. Custom callbacks now require `AdmittanceBounds(maximum_magnitude,
derivative_bound)` for both feeds: the first bounds `|G|` at every
nonnegative frequency; the second bounds `|dG/dω|` on each requested
interval, or is a global constant. A rigid feed uses its constant magnitude
and zero derivative. Finite sampling cannot prove these caller-supplied
bounds or absence of RHP poles. The result is conditional floating-point
numerical evidence, not an interval-arithmetic certificate or hardware
stability validation. `minimum_distance_to_minus_one` is a conservative
lower bound; reported gain crossovers are refined bracketed crossings and
do not enumerate unbracketed crossings or tangencies.

| Check (`tests/test_feed_coupled_instability.py`) | Result/criterion |
|---|---|
| Closed form against the composed component ratios | equal to 1e-12 |
| Large tank capacitance | tends to Eq. (5.4.2-2) to 1e-9 |
| Real-frequency form Eq. (5.4.2-3) | equal to 1e-12 |
| No manifold: quarter-wave and half-wave line resonances | admittance 0 and `1/Z_i`, as p. 244 states |
| Example 5.5.1.1 with rigid feeds | 0 zeros at K = 0.894; 2 (a conjugate pair, "two zeros" on p. 253) at K = 1.372 and 2.75 ms; the count changes across 1.78 ms; crossover at the neutral frequency |
| Short line, empty manifold | reduces to the rigid result |
| Unstable pair below the former first frequency | count 2, independently witnessed by a high-precision Lambert W root |
| Long lag, signed gain, rescaled residence times | counts agree with analytic delay crossing formulas; insufficient budget refuses a verdict |
| Finite tank, zero line length and zero combustion lag | counts agree with the independent quadratic characteristic roots |
| Imaginary-axis boundary; hidden high-frequency resonance | unresolved contour/tail is refused |

The scan of Eq. (5.4.2-1) shows the last denominator term as
`s² C_m Z_l Z_i`; that form does not reduce to Eq. (5.4.2-2) as `C_t → ∞`.
The composed form `s² C_m C_t Z_l Z_i` does, and is used.

## Example: the modes of a chamber

For a 0.30 m diameter, 0.40 m long chamber with a 1100 m/s sound speed,
`chamber_acoustic_modes(..., maximum_frequency_hz=4500)` returns:

| Mode | 1L | 1T | 1T1L | 2L | 1T2L | 2T | 2T1L | 3L | 1R |
|---|---|---|---|---|---|---|---|---|---|
| Frequency (Hz) | 1375 | 2149 | 2551 | 2750 | 3490 | 3565 | 3821 | 4125 | 4472 |

## Verification

| Check | Result |
|---|---|
| Table 3.5.1, eight eigenvalues | within 1.5e-4 of the printed values; every root satisfies `J_nu'(s) = 0` to 1e-12 |
| Table 3.5.1, 1T2R entry | printed 8.5263 is not a root (`J_1'` = 2.7e-3 there); the root is 8.5363 |
| Example 5.5.1.1, oxidizer drop 0.5 p_c | K = 0.894, unconditionally stable (printed K = 0.894) |
| Example 5.5.1.1, oxidizer drop 0.3 p_c | K = 1.372, critical lag 1.779 ms (printed 1.372 and 1.78 ms); 2.75 ms is unstable |
| Chugging limit | satisfies `K e^(-s tau)/(1 + s theta) = -1` at the neutral frequency to 1e-12 |
| Crocco-Cheng limit, n = 0.6 to 5 | satisfies Eq. (5.3.1-18) to 1e-12; n = 1 gives `tau = (pi/2) theta` |
| Fig. 4.2.2a | minimum `n/G_R = 0.5` at `tau/t_w = 0.5` |
| Q-factor example (Q = 55 at 550 Hz) | a 1 Hz error in one half-power frequency changes Q by 10 %, as p. 524 states |
| Helmholtz single resonator vs flat-array form | identical when `A_o/V = a/L_cav` |

![Chugging stability map of SP-194 Example 5.5.1.1: the neutral curve of Eq. 5.5.1-6 in the plane of injection gain K and time lag over residence time, with the engine at an oxidizer drop of 0.5 p_c (K = 0.89, stable) and 0.3 p_c (K = 1.37, unstable at its 2.75 ms lag)](media/validation-chugging-map.png)

*SP-194 Example 5.5.1.1 (O2/H2, r = 5.2, theta_c = 0.7 ms, lag 2.75 ms).
Shaded: unstable region of Eq. (5.5.1-6). Reducing the oxidizer injector drop
from 0.5 to 0.3 of chamber pressure raises the gain past 1 and the engine
chugs.*

## Limits

- Rigid-walled circular cylinder, nozzle end treated as a velocity node (the
  lowest-order approximation in Sect. 3.5.1). Transverse and mixed modes
  with a mean flow carry a `mach_correction_not_modelled` warning: SP-194
  gives the `(1 - M^2)` factor only for longitudinal modes.
- `chugging_stability` assumes purely resistive injectors and decoupled feed
  lines (the example's assumptions); `feed_coupled_stability` adds lines,
  manifolds, orifice inertance and tanks, but not pumps, branch lines or
  lines with area change, and assumes the feed adds no right-half-plane
  poles.
- `G_R` for the high-frequency neutral curve must be supplied; computing it
  needs the axial velocity and droplet profiles and the nozzle admittance.
- No baffles, liners in the chamber eigenvalue, annular chambers or
  triggered (nonlinear) instability.

## References

- **NASA-SP-194-1972**: <https://ntrs.nasa.gov/citations/19720026079>

Full record in [REFERENCES.md](REFERENCES.md).
