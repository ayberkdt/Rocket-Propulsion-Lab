# Thermodynamic record and catalogue integrity

Review base: `025bb9e58239116c43f075b9e6f720dd9ad8b081` on
`assurance-current-main`. Both Python 3.11 and 3.12 GitHub jobs for that
previous Nyquist/wall checkpoint completed successfully. This is a later
repair batch; previous passing counts do not verify these edits.

## Defects reproduced before repair

- An external coefficient list could change a frozen NASA interval's heat
  capacity after construction. Clearing a mixture's external component list
  could change its molar mass to zero. NaN coefficients could reach a
  species' successful property result.
- Editing a returned cached dictionary could make lookup of H2 return N2.
  The extended gas/condensed/reactant catalogue had the same shared-state
  problem in all three tables.
- `{'H2': 0.4, ' h2 ': 0.6}` was treated as two mixing species: the
  spurious entropy at 300 K was about 2775.825 J/(kg K). Even a pure-gas
  fraction within the accepted sum tolerance retained non-unit weight.
- Invalid direct records bypassed factories: empty or overlapping interval
  sequences, nonfinite data, inconsistent mass/mole fractions and invalid
  property identities could be admitted. At a positive subnormal pressure,
  `log(p/p_ref)` could fail after the ratio underflowed to zero.

These are integrity defects, not measured errors in a selected engine.

## Contracts

NASA intervals own exactly nine finite real coefficients and finite,
positive, increasing bounds. The fixed-arity check reads at most ten
iterator entries; its adversarial stream fails safely if this budget is
exceeded. Formulae and ordered contiguous interval sequences are owned
tuples. Ordinary atom counts are finite and positive; signed electron
inventory `E` is retained. This does not add ion modelling to the existing
equilibrium provider, whose documented assumptions still apply.

Gas species, condensed phases and assigned-enthalpy reactants validate
their scalar fields and nested inputs. Catalogues use canonical lookup
keys; they preserve source spellings such as `C(gr)`, `AL2O3(a)` and
`C2H2,acetylene` in records and results. Private tables are immutable.
Public catalogue functions return independent, locally editable dictionary
snapshots of immutable records. Editing one snapshot cannot affect later
snapshots, lookup or provider calculations. Extended catalogue snapshots
still support `dataclasses.asdict` and strict JSON serialization.

The NASA polynomial methods remain numerical primitives: they accept a
positive finite temperature without declaring it physically inside the
fit. This preserves assigned-enthalpy checks at 298.15 K for records whose
first fitted interval begins at 300 K. `Species.properties` selects the
owning interval and refuses state extrapolation. Overflow or nonfinite
results cannot become successful primitive/property outputs. This is not
a proof of polynomial evaluation accuracy at every representable input.

Mixtures own their components and require one component per canonical
species identity. Factories validate positive real fractions before
coalescing aliases, preventing negative amounts from cancelling through an
alias. Totals must be within the existing 1e-9 input tolerance; accepted
roundoff is normalized before conversion and mixing. Direct records have
a tighter 1e-12 sum tolerance and must also satisfy the mass/mole relation.
Their accepted rounding is normalized; conflicting bases are rejected,
rather than silently choosing one. A common species temperature range is
required. No invalid component is silently dropped.

Property records check finite values and positive state quantities.
Separately rounded fields may have `cp == cv` and `gamma == 1` even when
the positive exact heat-capacity gap and ratio are valid: a 100-digit
constant-cp witness guards against refusing this representation. Rounded
ordering is `cp >= cv > 0`, `R > 0`, `gamma >= 1`; `cp-cv=R` has a 5e-13
relative tolerance plus the operands' binary64 ULP resolution. Other
identities (`gamma=cp/cv`, molar/mass conversions and `R_specific*M=R_molar`)
use 5e-13 relative roundoff tolerance. These are record-consistency checks,
not an outward-rounded proof of a fit or its evaluation. They do not establish physical
admissibility, source authenticity or accuracy. Pressure entropy uses
`log1p` near the reference and a difference of logarithms away from it,
avoiding an underflowed positive pressure ratio.

## Evidence and remaining scope

`tests/test_thermochemical_record_contract.py` attacks both constructors
and factories, mutates caller-owned inputs, pollutes snapshots, checks
source-name preservation and serialization, and compares mixing and
pressure entropy with independent 80-digit identities. A constant-cp
polynomial gives a separate analytic reference-enthalpy/entropy control.
The existing CEA gas/condensed equilibrium, nozzle, transport and kinetic
comparisons remain necessary integration checks; the source fixtures and
their physical tolerances were not loosened to accept these changes.

Profiling showed repeated scalar normalization dominated property-record
construction. An exact-native-float fast path preserves finite checks and
signed zero; float subclasses still take the conversion/error protocol.
The existing adversarial core tests guard this distinction. No value check
is memoized or dropped and no cross-campaign runtime guarantee is inferred.

Scalar normalization also avoids rewriting finite native floats that are
already normalized, while checking their finiteness on every call. Common
literal annotations take a direct classification path. Field metadata is
still read each time; no schema or field value is cached. Nonfinite floats,
float subclasses and other numeric types retain their original conversion
and structured-error path. Tests exercise a value becoming NaN on a later
call, failure without partial writes, changed field annotations, subtype
conversion, declared error identity and unexpected application exceptions.
This reduces normalization overhead; only completed end-to-end checks can
establish whether the concurrent API timeout has been resolved.

The first expanded run caught an implementation regression that uppercased
printed phase names, and an unnecessary restriction on raw reference
polynomial evaluation. Those failures were retained and the behavior was
corrected. The first full run retained a 600-second API-fuzz timeout; a
separate high-precision check also caught overly strict rounded-output
ordering, so that intermediate tree is not an accepted checkpoint.
Final source-bound full-suite, mutation, lint/type and package
results belong in `output/assurance/oct27_thermo_records/`; only completed
checks are passes.

Open work includes primary-source coefficient and fit-constant consistency,
numerical conditioning across the full polynomial domain, caloric jumps
at piecewise boundaries and inverse reachability, transport records and
other model/catalogue contracts, and the API fuzzer's mandatory-refusal
assertions. Structural interval contiguity does not prove monotonic or
gap-free enthalpy/entropy images. No formula-to-atomic-mass verification,
universal fit-admissibility proof, hardware qualification or TRL 5
completion follows from these checks. Motor/feed design, mission and
relevant-environment test evidence remain unselected or absent.

## References

- NASA/TP-2002-211556, *NASA Glenn Coefficients for Calculating
  Thermodynamic Properties of Individual Species*, NASA9 representation.
  <https://ntrs.nasa.gov/citations/20020085330>
- NASA RP-1311 Part I, ideal-mixture caloric and entropy relations.
  <https://ntrs.nasa.gov/citations/19950013764>
- Pinned CEA input data, also retaining signed electron inventory and
  original species spellings; its hash/provenance is recorded in the
  packaged source metadata. The present batch does not change its fitted
  coefficients. <https://github.com/nasa/cea>
