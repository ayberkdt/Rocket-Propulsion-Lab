# Nyquist and channel-wall boundary review

Review base: `bd449cbb95b8331c4ecd60f85d5f73fd7a944af1` on
`assurance-current-main`. The previous review turn made progress by
producing concrete counterexamples which changed the repair priorities.
This repair covers those two numerical paths; it does not complete the
whole function/catalogue audit or qualify motor/feed hardware.

## Counterexamples and corrected contracts

The old Nyquist grid began at `1e-3/theta`, leaving an unresolved jump from
zero. A rigid-feed example with `theta=0.0007 s`, lag `0.36 s`, and gain
`1.0001` returned stable with the default 40,000 samples. The characteristic
`1+theta*s+K*exp(-tau*s)` has a root with real part
`0.0002256996974 /s`, imaginary part `8.7097109235 /s`. An 80-digit Lambert W
regression substitutes the root into the characteristic and checks the
positive real part, independently of frequency-grid phase unwrapping.
Long lags also caused large undercounts, not merely inaccurate margins.

The new contour uses bounds on the response and its frequency derivative.
Every accepted interval lies in an origin-excluding disk about its
midpoint. This makes replacement of that curve segment by its chord
preserve winding. The analytic global amplitude envelope also excludes
tail encirclements. Terminal phase is explicitly closed before reflection
across the real axis; negative DC characteristics can have an odd count.
Budget exhaustion, boundary roots, invalid bounds and nonfinite responses
are errors with no verdict. Gain crossover output is restricted to
bracketed crossings; it does not claim to enumerate tangencies.

Custom admittances require global amplitude and interval derivative bounds,
plus the caller's real-system/no-RHP-poles assumptions. Observed violations
are rejected, but finite observations cannot validate an arbitrary callable's
global contract. Built-in passive feeds derive bounds automatically.
Floating-point enclosures with rounding slack are conditional numerical
evidence, not rigorous outward-rounded interval arithmetic.

For a physical feed on the imaginary axis, write
`G=A/(Zi*A +/- i*B)`, so `|G|<=1/Zi` and
`|dG/domega|=|A'B-AB'|/|D|^2`. For the constant-pressure tank,
`A=cos(omega*tl)-omega*Cm*Zl*sin(omega*tl)`, `B=Zl*sin(omega*tl)`;
the Wronskian magnitude is `Zl*(tl+Cm*Zl*sin(omega*tl)^2)`.
For a finite tank,
`A=omega*(Ct+Cm)*cos(omega*tl)+(1/Zl-omega^2*Cm*Ct*Zl)*sin(omega*tl)`
and `B=cos(omega*tl)-omega*Ct*Zl*sin(omega*tl)`; the Wronskian is
`tl*(1/Zl+omega^2*Ct^2*Zl)+Ct+Cm*B^2`. A determinant/Frobenius bound and
a midpoint Lipschitz bound lower-bound `|D|`. These formulas bound the
whole interval, including unsampled resonances. High-precision derivatives
are additional consistency checks, not a replacement for this derivation.

The old 2-D hot boundary added imposed heat directly to the cell RHS, while
reconstructing a source-free surface temperature. A 2 MW/m² source produced
a 2 MW/m² local Fourier/Robin imbalance even with a near-zero reported
global energy residual. The half-cell and film must be eliminated together:
`q*A=G*(Taw+q_source/h-Tcell)`, `Ts=Tcell+(dz/(2k))*q`.
The sparse reference now has separate surface unknowns, instead of repeating
the same eliminated source term. Local balances and the public heat result
are checked for radiation, signed imposed heat and combined sources.

Malformed scalar values are refused rather than coerced. Unrepresentable
or nonpositive absolute temperatures cannot become successful wall results.
A separate storage/work cap prevents a 20,000-cell wide-grid input from
allocating about 100 million entries in each band matrix. These budgets
are library resource policies, not model validity or runtime guarantees.

## Verification and open scope

The regressions were first executed on the faulty implementation and
produced the expected failures. Targeted thermal, conjugate-march,
geometry and published-comparison tests are included, together with API
integration, explicit bound-contract attacks, neutral modes, all analytic
delay crossings, signed gains and time-scale invariance. The complete
suite, lint/type checks and isolated mutation evidence are recorded with
source fingerprints under `output/assurance/oct27_repairs/`; only actual
completed checks may be reported as passes.

The CEA 3.3.4 vendor stub contains a syntax error in parameter ordering.
Only that optional third-party import boundary is treated as untyped by
mypy; the repository adapter and other project source remain checked.

Open findings from the base review include species alias mixing entropy,
mutable cached catalogues, NASA data/fit-constant consistency, and API fuzz
tests which do not distinguish mandatory refusal from finite successful
JSON. Passing the existing API fuzzer addresses unexpected exceptions and
JSON shape, not that missing refusal assertion. Other models, catalogues
and claims still require review. No full-system or TRL 5 completion claim
follows from these repairs; design/mission selection and relevant-environment
motor/feed test evidence remain absent.

## References

- Model equations and assumptions: NASA SP-194 (1972), Sects. 5.2, 5.4 and
  5.5.1, as located in [COMBUSTION_INSTABILITY.md](COMBUSTION_INSTABILITY.md).
  <https://ntrs.nasa.gov/citations/19720026079>
- Fourier conduction and surface energy balance: Lienhard and Lienhard,
  *A Heat Transfer Textbook*, v6.00 (2024), Chs. 1–2; see
  [REGENERATIVE_COOLING.md](REGENERATIVE_COOLING.md). <https://ahtt.mit.edu/>
- Lambert W definition, `W*exp(W)=z`: NIST DLMF §4.13, Eq. (4.13.1),
  checked 2026-10-27. The characteristic-root transformation and envelope
  derivations above are calculations for this review, not formulas claimed
  to be printed in that source. <https://dlmf.nist.gov/4.13>
