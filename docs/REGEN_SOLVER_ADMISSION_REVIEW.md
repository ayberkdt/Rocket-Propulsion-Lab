# Cooling solver controls and material applicability

Baseline: `f677e71cd072fad1cae3c12bab8e846aa471b86e`. This review covers the public
`solve_regen_chamber` option boundary and propagation of the material evidence it
actually uses. It does not qualify a liner or its operating envelope.

## Reproduced problems

- Boolean and numeric-string conductivity overrides could reach fluid calculations.
- A boolean closeout thickness could be accepted as zero or one metre; `None`
  leaked a bare `TypeError` instead of a field-labelled refusal.
- A malformed closeout conductivity, including NaN, was ignored when the closeout
  was inactive. Inactive optional values still need a coherent input contract.
- Axial dispatch occurred before the radiation option and some material inputs
  were checked; malformed options could reach the iterative solver first.
- Caller-owned conductivity scalar subclasses could be handed unchanged to axial
  iteration. Constructor normalization elsewhere did not protect this boundary.
- Results reported a numeric material temperature difference while dropping the
  criterion's source condition and whether conductivity came from a catalog or an
  explicit override.

## What changed

The entry point checks the required record types and owns text options before
membership tests. Physical options use the shared strict finite-real conversion:
closeout thickness is nonnegative; supplied liner/closeout conductivities are
positive, even if the closeout is absent. Zero thickness and optional `None`
conductivity remain supported. A supplied valid closeout conductivity may be
inactive at zero thickness; no fabricated heat-transfer role is assigned to it.

These checks and material selection precede axial dispatch or fluid state work.
The selected constant scalar values are native floats at the iterative boundary.
This does not constitute a complete audit of caller-created nozzle-contour
records, all source profiles, or the numerical iteration as a whole.

`RegenChamberResult.wall_property_use` preserves the selected constant conductivity,
its basis, the cited catalog property when actually used, and the separate
temperature criterion and basis. Catalog properties retain units, source locator,
reference condition and any published range. An explicit override has basis
`user_constant` and **no catalog conductivity citation**, even if its number equals
the catalog. Keeping a material for its temperature comparison does not turn a
user override into published property evidence. A run without a material leaves
the criterion and station margins unset.

Machine-readable warnings carry the constant-property assumption and criterion
condition through the existing JSON/UI warning channel. A melting-point fallback
retains basis `solidus_temperature`; it is never relabelled as a service limit.
The existing heat, pressure and temperature calculations remain unchanged for
valid constant-property inputs.

## Primary-source check

The GRCop-84 conductivity paragraph was visually checked at **printed p.8 / PDF
p.12** of [NASA/TM-2005-213566](https://ntrs.nasa.gov/api/citations/20050123582/downloads/20050123582.pdf):
305–320 W/(m K). The catalog's locator was correct; PDF text extraction placed
some footer text misleadingly. Its lower-end selection of 305 is retained as a
reference constant, not a newly fitted temperature law.

The 973 K comparison was visually checked at **printed p.13 / PDF p.17** of
[NASA/CR-2000-210055](https://ntrs.nasa.gov/api/citations/20000064095/downloads/20000064095.pdf).
The stated criterion is 100 MPa yield strength. It cannot supply the stress,
creep, fatigue or life acceptance of an unspecified liner. The two reports also
discuss lot/processing and impurity effects; their property descriptions should
not be transferred indiscriminately between manufacturing conditions.

Downloaded source hashes and page renders are preserved under
`output/assurance/regen_options_review/sources/`. No property value was changed to
make a thermal test pass, and no temperature-dependent constitutive law was
inferred from a constant catalog range.

## Verification scope

`tests/test_regen_solver_contract.py` separates sentinel checks of admission/order
from real coupled solutions. The latter use the established LOX/CH4 equilibrium
gas and fuel-flow regression fixture. A provisional low-flow manual fixture
correctly encountered the methane transport range guard; that limit was retained.

Checks include all optional scalar types/domains, inactive options, text subclasses,
required record types, pre-axial refusal, scalar ownership, valid endpoints and
optional None. Real solves compare catalog versus equal-valued user conductivity,
absence of invented citations/limits, GRCop-84's conditional criterion, copper's
melting basis, and strict JSON through the actual HTTP route.

The fault inventory deliberately removes checks and evidence, forges override
citations, mislabels melting as service capability and drops warnings. Selected
fault runs require passing originals, assertion failures on the changed code and
passing byte-restored controls. This is scoped numerical/interface/source
verification. Actual material lot identity, k(T), stress, cyclic life and the
engine/feed relevant-environment evidence remain open.

The measured scoped run passed 85 test methods and 778 subtests without failures,
errors or skips. All 12 selected mutations were killed, with 24 passing
original/restored controls. Source/test hashes were checked against the archived
mutation baseline. These counts do not establish package-wide correctness or
readiness of hardware.

## References

- Ellis, NASA/TM-2005-213566, printed p.8 / PDF p.12, conductivity range and
  manufacturing-condition discussion:
  https://ntrs.nasa.gov/api/citations/20050123582/downloads/20050123582.pdf
- Ellis and Keller, NASA/CR-2000-210055, printed p.13 / PDF p.17, conditional
  973 K comparison at the stated yield-strength criterion:
  https://ntrs.nasa.gov/api/citations/20000064095/downloads/20000064095.pdf
