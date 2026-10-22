# Shared programmatic normalization and literal record review

This review covers the shared scalar/container validators, supported dataclass
scalar normalization, and the related manual literal metadata boundaries in
tabulated artifacts and burn records. It does not establish validity of every
consumer model, physical input, arbitrary annotation or application callback.

## Scalar representation contract

Real physical fields normalize accepted Real values to binary64. Booleans are
not physical numbers. Overflow, non-finite values and conversion that erases a
supplied nonzero value to zero fail with a labeled DomainError. Actual zero,
signed zero and representable subnormals remain valid. Ordinary rounding to a
nonzero binary value remains part of normalization; this is not exact rational
storage or a measurement uncertainty bound. An input already rounded to zero
before reaching the helper cannot be distinguished from an intentional zero.

Integer fields accept Integral values excluding booleans and retain arbitrary
integer precision. They do not accept floats or Fractions merely because those
values happen to be integral. Numerical/application-specific count limits remain
the caller's responsibility. Private exact arithmetic paths such as _fraction
retain their rational intermediates; this review introduces no physical minimum
thrust, flow, mass or duration threshold.

Known numeric conversion errors and malformed iteration/mapping protocols become
field-labeled errors. Already typed domain errors retain their type and structured
context. Unrelated application failures such as RuntimeError propagate rather
than becoming misleading input errors. These are data/API contracts, not a
sandbox for arbitrary Python object code.

## Literal text and model identity

Text helpers validate the actual underlying Unicode characters, not overridable
strip/encode/string methods. Valid text normalizes to plain str without trimming,
case folding or Unicode normalization. Strings and mapping keys follow the same
rule. Numeric-text parsing reads the literal characters, so a str subclass cannot
turn a literal 1.25 into 99 through its __float__ implementation. JSON source text
also becomes a plain literal string before decoding.

Tabulated schemas, source identifiers/digests, interpolation policies and phases
store the normalized literals. Stream names, propellant/tank identifiers and
nullable summary binding labels also retain plain literals. Enum conversion uses
literal text before membership lookup; custom equality/hash methods cannot turn
an unsupported role or policy into a recognized member. Existing valid enum
members remain valid. Non-depleting streams still permit empty tank identifiers,
but their text must be valid Unicode and is normalized even on that branch.

These repairs prevent an accepted artifact from containing a forged digest or
schema spelling that fails its own JSON reopening. A source digest's valid format
and a reproducible model hash do not authenticate a source or prove engineering
review. Reference identifiers remain caller-supplied metadata.

## Ordered containers and dimensions

Numerical/component sequences require ordered finite iterables. Mapping and
Set-interface containers, including their set-like views, cannot implicitly
define component order. Text and raw binary buffers (bytes, bytearray, memoryview)
cannot become integer-valued numeric vectors. Lists, tuples, ordered iterators,
legacy indexed sequences and numeric arrays retain their order.

Three-component vectors inspect at most four items; inner name/value pairs inspect
at most three. Excess data is rejected rather than truncated. These bounded reads
also reject infinite iterators for fixed-cardinality inputs. Variable-size Python
iterables must be finite; no wall-time guarantee for arbitrary caller iterators is
claimed. Known protocol failures are labeled rather than leaking bare data errors.

Matrices are rectangular, with empty outer or zero-width dimensions still allowed
for callers that support them. Square shape, positivity, symmetry, covariance
properties and physical dimensions remain caller checks. Containers are detached
from input containers. Records retain the supplied typed record objects rather
than promising deep copies of arbitrary application objects.

Name/value sequences preserve repeated names as separate entries; a consumer that
requires uniqueness must check it. Mapping inputs instead have unique-key
semantics: a malformed custom mapping cannot silently lose duplicate items.
real_mapping(None) retains its existing empty-mapping meaning. allow_blank is a
real boolean option, not a truthiness shortcut.

## Supported dataclass fields

Supported primitive and optional-float spellings normalize without evaluating
annotation code. Actual Optional[float], PEP 604 unions and the documented quoted
scalar/typing spellings are recognized. Unsupported annotations and arbitrary
aliases remain model-owned; their values are not fetched or evaluated by this
helper. The helper is not a general runtime type system.

All recognized field validation completes before writes. A validation failure
does not partly normalize earlier fields of an existing dataclass. This does not
promise rollback of arbitrary custom descriptor side effects or later model-domain
checks. Call after supported fields are initialized, before model checks. Passing
a dataclass class instead of an instance, or malformed record-kind configuration,
is a programmer TypeError; malformed external values have labeled domain errors.

## Adversarial evidence and mutation accounting

Thirty-two final core methods cover the contracts and ordinary controls, including
public coupled-record cases. They are replayed unchanged against the verified
pre-repair package archive. Direct public probes separately reproduce nonzero
Fraction-to-zero conversion, bytearray vectors, indexed mapping pairs, ignored
optional field validation, partial dataclass writes, overridden numeric text and
accepted artifacts that fail their own JSON reopening. Historical failures include
deliberate contract extensions and are not a distinct-defect count.

The test oracle explicitly diagnoses recognized wrong data-protocol exception
types with assertions; unrelated test/programming errors remain errors. Direct
unittest failure/error outcomes are retained because JUnit can combine subtests.
The curated runner counts only assertion failures without unrelated errors.

Twenty-five new variants target real/count normalization, order/binary/cardinality
guards, error context, Unicode/literal conversion, mapping semantics, ragged
matrices, optional annotations, partial writes and coupled metadata/enum gates.
An exploratory run also exposed two older local Unicode deletions that survived
their selected tests. The branch audit found a missed non-depleting-stream path:
the strengthened public oracle killed that deletion before the manual stream
repair. The unmutated implementation also accepted forged encoders/roles.

After those real defects were repaired through shared literal normalization, the
strengthened public oracles still reject their inputs when either local Unicode
check is deleted. The two now redundant variants are retired from the active
set, not counted as kills. Original reports, the branch probe and the post-repair
outcomes are retained. This is evidence about the tested public paths, not a
private-helper equivalence theorem or exhaustive mutation analysis. Earlier review
counts describe their own checkpoints; they must not be added as current totals.

## Compatibility and remaining work

Previously accepted erased nonzero Real values, unordered/binary components,
ragged matrices, malformed pair protocols and truthiness-based blank options now
fail explicitly. Normalized text loses subclass behavior while preserving its
literal characters. Supported optional scalar fields now validate instead of
being silently skipped; regenerated identities can therefore change where older
objects stored unnormalized integer-versus-float metadata.

Further bridge/timeline/exchange schema work and the complete caller, catalogue,
function and technical-claim review remain open. These scalar/representation tests
do not establish physical uncertainty, calibration or the hardware acceptance
domain. Mission and motor/feed configuration remain undefined; all nine hardware
evidence gates are unmet. Software results do not establish TRL 5 or certification.

## References

- [Python 3.12 numeric abstract base classes](https://docs.python.org/3.12/library/numbers.html):
  Real/Integral protocol background. Boolean exclusions and accepted annotation
  spellings are repository contracts, not requirements imposed by that hierarchy.
- [Python 3.12 fractions](https://docs.python.org/3.12/library/fractions.html):
  exact rational representation of a supplied binary float, distinct from a
  decimal approximation or measurement uncertainty.

- src/rocket_propulsion/core/validation.py: shared normalization implementation.
- src/rocket_propulsion/propulsion/burns/models.py, tabulated.py and
  tabulated_documents.py: related manual scalar/metadata boundaries.
- tests/test_core_validation_contract.py and tests/test_burn_record_contracts.py:
  independent contract probes and the strengthened stream branch.
- scripts/run_feed_mutations.py: active source-bound semantic fault identities.
- [TABULATED_DOCUMENT_CONTRACT_REVIEW.md](TABULATED_DOCUMENT_CONTRACT_REVIEW.md)
  and [BURN_RECORD_CONTRACT_REVIEW.md](BURN_RECORD_CONTRACT_REVIEW.md): earlier
  document/record scopes and historical checkpoints.
