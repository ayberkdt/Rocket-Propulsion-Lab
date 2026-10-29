# API numeric admission and refusal-test review

This change covers declared top-level scalar numbers at the dispatcher. It
does not qualify a physical model, nested document schema, motor or feed system.

## Request contract

For a registered calculation path, the request must be a JSON object. Every
supplied declared numeric field is checked before the adapter runs, even if
the selected branch will not use it. A malformed number raises `InputError`
with that field's name. Thus an invalid inactive control is no longer silently
ignored, and numeric errors can precede missing-field or selector errors.
Unknown paths retain their `KeyError` behavior.

The immutable `api.numeric_contracts` catalogue declares 413 binary64-number
fields and three exact-integer fields across 61 registered paths. These are
route/field pairs, not 416 unique names. Scalar readers reached through shared
API helpers are included. Runtime dispatch uses declarations, not AST parsing.

Binary64 admission reuses the existing `_number` reader: booleans, arrays,
objects, nonnumeric text, NaN, infinities and overflowing conversions fail.
Numeric form strings remain supported. Null and empty text retain their
adapter-specific optional/default/error semantics. Admission does not rewrite
the request or enforce all required fields, ranges, integrality, cross-field
constraints or selected-model validity.

The uncertainty route's `seed`, `sample_count` and `sobol_base_samples` instead
use their existing exact-integer reader. A JSON integer larger than binary64,
integral decimal form text and `1e309` are type-valid there; bounds remain with
the handler. The API also applies an explicit 4300-decimal-digit resource cap,
or a lower configured host integer/string limit, before decimal expansion.
A large exact seed within this budget is not forced through float conversion.
See [API_INTEGER_RESOURCE_REVIEW.md](API_INTEGER_RESOURCE_REVIEW.md).

## Defects and tests

The previous fuzz `probe` returned success for any finite JSON response to a
malformed numeric field, and for an error about any other field. Fault controls
now require both finite-success and unrelated-error injections to be reported
as failures. A correctly named refusal and a finite semantic-extreme result
remain passing controls.

The strengthened type-garbage family requires refusal for declared scalar
numbers and checks the error's field. It includes helper-read numbers in
addition to direct getters. Non-numeric fields still receive only a
typed-error/strict-JSON smoke check. Single-field requests are sufficient to
isolate the new admission contract; they are not baseline-valid full-model
fixtures. Baseline-valid semantic/branch and nested-schema audits remain open.

`test_api_numeric_admission.py` isolates dispatch with inert adapter bodies.
For each malformed scalar it checks refusal, the exact field, and that no
model body ran. Controls preserve request/value identity, optional placeholders,
numeric forms, large integers and integral exponent text. Additional checks
cover non-object payloads, conversion overflow, fractional exact integers,
immutable catalogues and a real bell-contour request with an invalid inactive
numeric control.

The static test audit follows positional and keyword payload aliases through
API helpers. Variable-key readers are explicitly bound to source default maps;
unresolved reader keys fail the catalogue guard. A computed-key fixture verifies
that unresolved sites are reported; another exercises keyword-only passage of
the payload and reader arguments. This is not a general Python dataflow proof:
nested library schemas, arbitrary payload transformations and unrecognized
reader idioms require additional review.

The first admission baseline recorded 2918 failed reports, six passing methods
and 2551 passing subtests. These reports include stub-contract cases and are
not a count of distinct production defects. Refusal-oracle baseline trials
retain the initial missing-keyword error separately from the two subsequent
assertion failures against the old behavior. A later test incorrectly treated
`1e309` as invalid exact-integer text; that failure was retained and the test
was corrected, with positive controls added. Production integer semantics were
not narrowed to make it pass.

Current full regression, selected fault injections and package evidence are
recorded against their exact source hashes under
`output/assurance/oct28_api_admission`. A passing software suite does not
establish complete malformed-input coverage or motor/feed TRL 5. Hardware
selection, relevant-environment measurements, uncertainty/holdout comparisons
and scaling remain separate open requirements.

## References

- RFC 8259, JSON values and numbers (sections 3 and 6):
  <https://www.rfc-editor.org/rfc/rfc8259>.
- NASA NPR 7123.1D, Appendix E, hardware and software readiness criteria:
  <https://nodis3.gsfc.nasa.gov/displayDir.cfm?Internal_ID=N_PR_7123_001D_&page_name=AppendixE>.
