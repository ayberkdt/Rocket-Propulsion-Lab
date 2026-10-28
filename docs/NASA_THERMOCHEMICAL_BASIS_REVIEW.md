# NASA thermochemical data and dimensional basis

Review base: `e44158065454dab3338cc55ac646f32dc5142ea8`.

NASA/TP-2002-211556 defines its polynomial convention with
`R = 8.314510 J/(mol K)` in Appendix A, printed page 7. The pinned CEA
Python source exposes the equivalent `8314.51 J/(kmol K)` in
`source/bind/python/cea/constants.py`, Git blob
`d18bb89b0e65b6e0abbd576830037bb2be6c40c7`. The source report's H record
(Appendix D, printed page 126) and the pinned `thermo.inp` agree on two
digits that differed in the bundled manual subset:

| H coefficient | Before | Source |
|---|---:|---:|
| 200-1000 K, b1 | 25473.6599 | 25473.70801 |
| 1000-6000 K, a2 | -0.181935441 | -0.1819354417 |

The ten gas records' other 178 coefficients, interval bounds, formulae,
molar masses and assigned formation enthalpies match the pinned records.
The subset retains only its first two intervals through 6000 K, even where
the source contains additional higher-temperature intervals.

## Why the basis must be consistent

Previously, multiplying the dimensionless NASA coefficients by the general
2019 SI R changed dimensional cp, h and s by about -5.699 ppm. Mixture
entropy and the ideal EOS also used SI R, but assigned liquid enthalpies
remained in NASA's published J/mol. External CEA h/r input decks used the
SI denominator. This is distinct from an explicitly converted, fully
specified alternative thermochemical convention.

The repaired thermochemical model keeps NASA's defining R throughout
polynomial conversion, ideal-mixture entropy and caloric/EOS identities,
gas and condensed equilibrium, thermodynamic derivatives, reaction
transport, the CEA Eucken estimate, CEA energy input units, kinetic
concentrations, detailed balance and finite-rate nozzle energy/entropy/EOS
closure. The old thermochemistry export
`UNIVERSAL_GAS_CONSTANT_J_MOL_K` now denotes this NASA model value.

The general `core.constants.MOLAR_GAS_CONSTANT_J_MOL_K` retains the SI
value. Arrhenius forward-rate activation energies converted to J/mol retain
their existing SI R convention; reverse rates use the NASA thermodynamic
equilibrium constant. An independent high-precision Arrhenius control
guards this distinction. Empirical forward rates and a thermodynamic
equilibrium constant have different defining roles.

Real-fluid and atmosphere fits retain their own defining values. Span's
nitrogen EOS, Setzmann-Wagner methane EOS and Friend-Ely-Ingham methane
transport happened to contain the same numerical R as NASA; these are now
registered under their own source records, without changing their values
or tying them to the NASA alias. Registered `exact=True` marks an exact
choice of model convention, not an exact physical measurement of R.

The local gas and condensed providers carry new `nasa2002` version tags.
Comparisons with older result files must account for the earlier basis and
H coefficients. Historical CEA measurements in existing documentation are
labelled as historical; they are not re-presented as revised measurements.

## Independent controls and limits

`tests/reference_data/nasa_glenn_primary_records.json` preserves raw source
records and line locations. Its extraction selected source records without
reading bundled coefficient values or using the production parser. When
`RPLAB_THERMO_INP` is provided, tests verify the complete pinned input hash
and each raw block against it; the reference blocks also support offline
tests. The primary source hash is
`fa7746572952d74e249e818a82a35c113829742fb421a308e167185528884363`.

`test_nasa_primary_basis.py` checks the 180 literal coefficients, source
metadata and a 90-digit cp/h/s evaluation at each interval's two endpoints
and midpoint. Its roundoff allowance uses the absolute term sum, not a
possibly vanishing result. It is a sample-wise numerical check, not a
rigorous error enclosure over every representable temperature.

Further controls independently check ideal mixing and pressure entropy,
pure-gas UV energy/EOS recovery, assigned-liquid HP closure in source
dimensional units, frozen derivatives, concentration/detailed-balance
units, the unchanged SI Arrhenius convention, Eucken scaling and external
CEA enthalpy input. The finite-rate fixed-composition kernel is separately
checked against source energy, entropy and EOS values; this unit control
does not by itself verify its coupled nozzle ODE. Published/source
reference fixtures and numerical acceptance tolerances are not widened.
Existing independent oracles use the defining NASA value rather than
following the production alias.

The first full run found one remaining internal golden value computed with
SI R: the N2/O2 79/21 mixture expected R/M = 288.19294 J/(kg K). An independent
90-digit calculation with source molar masses and NASA R gives 288.19458
to five decimal places. That computed expectation is corrected while its
original four-decimal-place assertion is retained. This changes no published
NASA/NIST reference datum. The failed full run is retained, and the corrected
test tree receives a new full run. Production inputs and the eleven selected
mutation test methods are unchanged by this test/document correction.

The initial test draft tried a UV state below the provider's present
1800-6000 K search bracket and read a nonexistent molar-mixture output
field. Those observations are retained; the basis test uses an interior UV
state and converts the actual mass-specific entropy field. Broader inverse
reachability, polynomial conditioning, physical fit validity, other
catalogue/model contracts and mandatory-refusal API testing remain open.

Completed verification belongs in `output/assurance/oct28_nasa_basis/`.
Passing these checks does not qualify a chosen motor/feed assembly. No
hardware configuration, relevant environment or physical TRL 5 evidence
is established by this repair.

## Primary references

- NASA/TP-2002-211556, Appendices A and D:
  <https://ntrs.nasa.gov/citations/20020085330>.
- NASA CEA constants at the pinned commit:
  <https://github.com/nasa/cea/blob/2f79a647e85737742fb08a127da00d1ea1d069da/source/bind/python/cea/constants.py>.
- Pinned NASA thermodynamic input:
  <https://github.com/nasa/cea/blob/2f79a647e85737742fb08a127da00d1ea1d069da/data/thermo.inp>.
- Independently registered real-fluid model conventions: Span et al.
  (2000), p. 1364; Setzmann and Wagner (1991), p. 1064; Friend, Ely and
  Ingham (1989), Table 1, p. 590. Their source URLs are recorded with each
  constant in `core/constants.py`.
