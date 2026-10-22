# Ensemble document and retained-record review

Follow-up: PROBABILITY_NUMERICS_REVIEW.md records the subsequent statistical
numerics policy. The retained 1e-12 derived-value checks described in this earlier
stage are replaced there by deterministic display agreement under a new required
policy marker. Earlier evidence retains its historical scope and source bytes.

This bounded review covers hierarchical ensemble raw records, their future
Sidera handoff, scalar/collection normalization of the seven retained ensemble
record classes, and builder control arguments. It does not establish the
numerical adequacy of every uncertainty distribution, copula, weighted moment,
variance calculation, native orbit consumer or physical model. Mission and
hardware selection, applicable independent measurements and all nine motor/feed
hardware evidence gates remain unresolved. No certification is inferred.

## Reproduced failures

The previous readers silently truncated fractional seeds, sample counts and
realization indices, and converted JSON numeric strings and booleans into
numbers. Unknown root/nested fields and extra pair components could disappear
without changing the normalized document or its hash. Tuples could enter a raw
JSON dictionary as arrays. String-subclass behavior could affect schema checks.

Rehashed documents could contain duplicate metric/contribution/evidence
identities, missing/duplicate/undeclared or out-of-bounds realization inputs,
unlinked seed claims and unequal weights under an equally weighted sampling
method. A zero-output conditional sample exposes false weights even when its
moments remain unchanged. The old closeness checks also admitted small negative
standard-deviation/closure values. Handoff end-time metadata could disagree by
one output ULP, and indices could have gaps.

A public single-scenario ensemble exposed another defect: the writer retained
integer zero for the no-zero-thrust probability, but the reader converted it
to floating zero. Their serialized bytes and hashes then disagreed. The frozen
records also retained caller-owned lists without taking snapshots.

## Raw and programmatic contracts

Both document roots and every owned nested record have closed declared key sets.
Strings, objects, arrays, booleans, finite representable real numbers and integer
identity/count fields retain distinct JSON kinds. In raw dictionaries, integers
require built-in int values, excluding bool, numeric text and integral floats.
Real-number fields permit built-in int/float values but reject bool and text.
Input pairs require exactly two array entries. Matrix rows and directions use
arrays with validated numeric entries. Nested burn artifacts keep their own
reader and content identity checks; these are not bypassed by the outer schema.

The previously supported omitted nullable standard_deviation and correlation
fields remain supported; their null defaults normalize explicitly. No other
undeclared field is silently discarded. Key order and numeric int/float
representations of real fields can normalize. Content hashes still identify
the normalized representation; they are neither signatures nor source evidence.

Retained ensemble records normalize annotated scalars through the shared typed
boundary. Root record collections, conditional parameter lists, reference and
warning strings, and realization input pairs are copied to immutable tuples.
Builder controls reject malformed seeds, methods, tolerances and record types
before sampling. These changes preserve a public writer/reader roundtrip and
prevent later caller list edits from changing retained evidence. They do not
fully re-audit every nested uncertainty-model constructor or sampling function.

## Internal semantic evidence

Duplicate identities cannot be hidden by dictionary reconstruction. Realization
names must uniquely cover the definition and values must lie within its bounds.
Equal weighting is the builder's declared policy: each conditional weight must
equal the represented result of 1.0 / sample_count. It is not a generic
importance-weighted schema. Retained base/scenario/channel seed relationships
must match the documented deterministic derivation. Reordering complete records
and input pairs remains valid with a corresponding normalized content identity.

Standard deviations, variances, fractions and closure magnitudes cannot be negative.
Rounded variance fractions can lie just above one in a canonical deterministic
case; they are compared with the recomputed diagnostics under the existing
tolerance, rather than rejecting the writer's own value. Existing 1e-12 moment,
factorization and aggregate-probability comparison tolerances are retained.
Those tolerances are software consistency rules, not measurement uncertainty or
mission acceptance limits. Their full numerical adequacy is separate work.

Handoff job identities must be unique and complete within each scenario. Its
end time is a checked binary-input sum, and its supplied output must equal that
deterministic display. An arbitrary absolute 1e-12 allowance no longer permits
discarding a different end-time declaration. Handoff-owned literals and record
collections are normalized and copied. The source ensemble hash is a declared
link: without the source document it does not authenticate that source or prove
that a native weighted-ensemble consumer executed these jobs.

Seed linkage and finite/bounded input checks do not independently replay the
sampling campaign, regenerate the primary/audit histories, or validate reported
stability differences against measurements. The convergence flag is checked
against its reported diagnostics and tolerance. Source metadata remains supplied
data, and a hash cannot establish independent physical validation.

## Evidence discipline

Public tests retain canonical, correlated, large-seed Monte Carlo, reordered,
nullable-legacy and numeric-real representation controls. Negative cases use the
original hash where normalization previously hid information, and independent
JSON hashing where the inconsistent claim itself is included in the document.
Historical replay uses unchanged final tests against the verified previous
package archive; failed methods are not a count of distinct defects.

The initial test attempt had two fixture errors: its stochastic data did not
contain a zero-standard-deviation metric, and a numeric representation control
incorrectly hashed raw int/float spellings as identical bytes. Original source
and results are preserved in output/assurance/history/ensemble_initial_fixture_failures.
The zero-variance case now constructs deterministic histories, and the numeric
control retains the normalized model hash. Neither mistake is counted as a
software defect or a mutation kill. Final counts and source bindings belong to
the verification reports. Test workloads, seeds and execution limits are not
reduced to obtain a pass.

The first mutation campaign had two surviving variants. A test helper captured
the original reader as a default argument, bypassing a replaced callable;
runtime binding now uses the public reader. The negative-deviation test was also
being rejected by an overly strict new upper bound on a rounded fraction. Its
valid baseline now must reopen before the negative field is injected. The
incorrect new upper bound was removed; the pre-existing comparison tolerance
was not increased. The original failed campaign is retained, and removing the
actual nonnegative-deviation guard must now cause an assertion failure. A fixture
class import was changed to a module import to prevent pytest from collecting
five existing tests twice; all original test modules remain in the full run.

## References

- RFC 8259: distinct JSON value kinds, object members and arrays. Closed schemas
  and the application identity/count rules above are local contracts.
  <https://www.rfc-editor.org/rfc/rfc8259>
- Python 3.12 numbers: integer/real protocols; the raw JSON boundary deliberately
  distinguishes the built-in kinds before programmatic reconstruction.
  <https://docs.python.org/3.12/library/numbers.html>
- Python 3.12 fractions: exact values of supplied binary floats for the checked
  handoff end-time sum. <https://docs.python.org/3.12/library/fractions.html>
