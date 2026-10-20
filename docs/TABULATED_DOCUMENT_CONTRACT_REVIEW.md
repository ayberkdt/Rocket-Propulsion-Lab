# Tabulated JSON and CSV document boundary review

This review covers tabulated artifact raw JSON structure, the shared JSON text
parser used by its caller, UTF-8 metadata and the SI CSV producer. It does not
establish a closed schema for every other persisted model, source authenticity,
measurement validity or hardware readiness.

## Normalized model identity

Artifact hashes identify the normalized supported model. They are not signatures
and do not authenticate a source or prove engineering review. Caller-supplied
reference identifiers and source SHA-256 text remain literal metadata; uppercase
and lowercase source digest spellings retain distinct literal model identities.
Reference presence is not a validation or qualification claim.

The v1 artifact has closed root, segment and summary field sets. Unknown fields
are rejected rather than silently dropped during reconstruction. Writer-emitted
fields, including warnings, must be present. Only reference_ids is optional; an
explicit empty reference array has the same meaning as absence. Arrays at this
decoded-document boundary are JSON lists, with two-element list flow pairs.
Programmatic constructors retain their separate documented iterable interfaces.

Physical fields and summary values must have the declared numeric types;
booleans cannot masquerade as numerical agreement through Python's True == 1
comparison. Finite JSON integers can normalize to the model's binary64 fields.
Only the equivalent-specific-impulse summary can be null, and it must still
reproduce the model. Hash and derived-summary checks remain mandatory.

All artifact identifiers, phases, references, warnings and stream names must be
UTF-8 encodable before hashing. Unicode text remains supported. Accepted input
containers are copied; later changes to the raw document cannot alter the model.

## Shared JSON text parser

Duplicate member names are refused at every depth before last-wins information
loss. Non-finite constants, overflowed decimal literals and nonzero decimal
literals that would round to zero are refused before the model builder runs.
Actual zero remains valid, including large exponent spellings. Decoded keys and
values are checked for invalid Unicode; escaped lone surrogates cannot bypass a
check of the original source text.

These parser rules apply to paths using core.validation.read_json_document.
Readers with separate parsers are not implicitly covered. Large finite integers
remain subject to their field contracts and the runtime's parser limits; this
review does not introduce an arbitrary document-size or duration restriction.

## CSV contract

Required SI columns are time_s, delivered_thrust_n and total_mass_flow_kg_s.
Supported optional columns are axial_thrust_n, phase and stream:<identifier>.
Headers must be unique; unsupported columns are errors. Every record must match
the declared width. Quoting is parsed strictly. Declared numeric cells cannot be
blank or missing, and invalid phase cells cannot become invented phase labels.
Only absent axial, phase and stream columns receive their documented defaults.

Every physical sample is checked for finite nonnegative quantities, axial force
closure and named-flow closure, including the final row of a HOLD trace. The
existing numerical closure margins remain model policy, not measurement accuracy
or hardware acceptance criteria. A positive literal outside nonzero binary64
range fails instead of becoming a false zero flow or thrust.

LINEAR uses adjacent endpoint samples. HOLD uses row i on [t_i,t_(i+1)); the
terminal row sets the final time and has no following firing interval. Its
physical values are still validated and its literal bytes remain bound by the
source hash. This replaces the overly broad claim that every input sample was
represented unchanged in interval output. No smoothing or new measurement is
introduced. Source hashes bind the exact UTF-8 CSV text, including newline style.

## Compatibility and evidence

Previously ignored JSON extensions, Python-only decoded containers, missing
writer fields, duplicated CSV headers, unknown columns, blank declared cells
and invalid terminal HOLD samples now fail explicitly. Preserve unsupported
engineering information in a controlled external record until its schema is
defined; do not silently discard it. CSV warning text now describes validation
and interval semantics, so regenerated CSV artifact identities can change.

Twenty-five adversarial methods cover these contracts and valid controls:
Unicode, literal source hashes, empty references, detached containers, quoted
phases, exact LINEAR/HOLD impulse and mass totals, duplicates, unmodeled fields,
boolean summaries, malformed quoting, ragged records and numeric range loss.
The final tests are replayed unchanged against the verified pre-repair archive.
Public-only counterexamples separately record actual old field loss, invented
CSV defaults and nonzero-to-zero conversion. Historical failures are classified
by assertion and non-assertion exception; neither is a distinct-defect count.

At the document-review checkpoint, seventeen additional semantic faults targeted
schema omissions, JSON member loss,
numeric/Unicode checks, CSV headers/width/defaults, terminal physical validation
and strict quoting. Only assertion failures count as killed faults. This finite
set does not prove exhaustive parser correctness or denial-of-service resistance.
The later [core normalization review](CORE_VALIDATION_CONTRACT_REVIEW.md) retired
one local UTF-8 deletion after shared validation made it redundant on the tested
public paths. It was not counted as a killed fault in the new active set.

## Remaining work

Bridge, timeline and other exchange-model raw schemas need their own complete
field/metadata review. Programmatic Real normalization outside these text readers
also needs a separate underflow contract: a follow-up synthetic Fraction input
can still normalize to zero in the general real helper. This observation is not
a hardware operating recommendation. Broader catalogue, function and claim
review remains partial. All nine motor/feed hardware evidence gates remain unmet;
software checks do not establish TRL 5 or certification.

## References

- Repository implementation: core/validation.py, propulsion/burns/tabulated.py
  and propulsion/burns/tabulated_documents.py under src/rocket_propulsion.
- Independent checks: tests/test_tabulated_document_contract.py and the existing
  persisted-reader, ledger and time-contract regressions.
- [FREE_SPACE_TABULATED_NUMERICS_REVIEW.md](FREE_SPACE_TABULATED_NUMERICS_REVIEW.md) and
  [TABULATED_TIME_CONTRACT_REVIEW.md](TABULATED_TIME_CONTRACT_REVIEW.md): separate
  supplied-binary arithmetic and clock ownership contracts.
- [RFC 8259](https://www.rfc-editor.org/rfc/rfc8259): external JSON interoperability
  reference; repository model hashes are not an authentication protocol.
