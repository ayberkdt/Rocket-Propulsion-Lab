# Exchange document schema and metadata review

This review covers the raw JSON/dictionary readers for propulsion bridges,
timelines, total-flow performance maps and stream-resolved performance maps.
It does not complete the other persistence, catalogue, solver or physical-model
reviews. Content hashes identify normalized contents; they are not signatures,
source authentication, calibration records or evidence of independent review.

## Defects reproduced through public readers

The previous readers discarded unknown fields at the root and in nested records.
Bridge reconstruction ignored declared artifact/surface hash links, operating-axis
names/units/limits, force/mass/time descriptions and reference identifiers. Timeline
reconstruction ignored each scheduled burn's declared bridge hash and unused bridge
registry entries. Those fields could disagree with the accepted normalized model
without changing its recomputed hash. This is inconsistent document acceptance,
not evidence that an authenticated signature was bypassed.

The final adversarial methods are replayed against the verified pre-repair package
archive. Their negative cases retain the original normalized hash so the oracle
distinguishes silent information loss from ordinary hash tampering. Controls cover
canonical reopening, evaluation equivalence, vacuum bridges, stream maps, legitimate
legacy omissions, reordered object keys and grid points, and independently rehashed
valid models. Historical failures are method outcomes, not a distinct-defect count.

## Reader contract

| Boundary | Enforced contract |
| --- | --- |
| All four document roots | Closed declared field sets and required metadata |
| Nested records | Closed axis, point, binding, constraint, parameter, scheduled-burn and jettison fields |
| JSON types | Native objects/arrays/text/numbers/booleans; no tuple, mapping-as-array, boolean-as-number or text subclass behavior |
| Numerical/text values | Finite representable numbers; UTF-8-encodable literal text; physical constraints still checked by model constructors |
| Bridge metadata | Artifact and surface hash links, operating axes, force/mass/time semantics and fixed implementation references agree with the rebuilt manifest |
| Timeline registry | Exactly the scheduled burn identifiers, with every supplied bridge validated |
| Scheduled burn | Declared bridge hash agrees with the bridge assigned to that burn |
| Nested artifact/maps | Their own raw schemas, numerical policies, model constraints and normalized content hashes remain enforced |

Schema checks run before reconstruction. Metadata consistency checks run after
reconstruction, in addition to existing model hashes and numerical-policy checks.
The comparison distinguishes boolean and number kinds while allowing equal integer
and floating representations. Each rejected field has a typed InputError. JSON
text still passes through the shared duplicate-member/non-finite/underflow decoder.

The bridge's optional dry-mass limit may explicitly be null, as written by its
writer. Absent legacy stream-surface and stream-hash fields are accepted only when
no stream surface is configured. A configured stream surface requires its correct
hash link. An explicit multilinear interpolation declaration remains accepted.
Nozzle exit area is omitted when unused; an explicitly present field must contain
a number. Point order may differ in a raw map: its sorted canonical identity and
evaluation are preserved, without promising identical storage tuple order.

## Implementation review

`exchange_documents.py` isolates structural object/array/text/number/boolean checks,
the kind-aware redundant-metadata comparison, and the three document-family
validators. The stream and total-flow map variants share the surface validator with
different point fields. The bridge validator checks its nested artifact and maps;
the timeline validator checks each registered bridge. Constructors and public hash
readers remain responsible for complete grids, unique identifiers, dimensional
constraints, physical ranges, ordering semantics and hash recomputation.

The generic checks are intentionally limited to finite built-in JSON trees at
these public document boundaries. They are not a sandbox for arbitrary Python
objects, a streaming parser, a resource-exhaustion defense or a general schema
language. Programmatic constructors retain their separate accepted input types.

Nine curated changes remove specific schema, array, metadata, link or registry
guards. Only test assertion failures without unrelated errors count as kills;
these variants are not exhaustive coverage of the schema or engine behavior.

## Verification status and retained failures

The first combined targeted run exceeded its unchanged 120-second deadline and
did not produce complete evidence. A subsequent isolated run exposed an overly
strict new dry-mass check; its source, unchanged test and failed report were
archived before allowing the legitimate null value. A later added control wrongly
required reordered grid points to have identical internal tuple order. That oracle
was corrected to require the documented canonical identity and evaluation instead;
the earlier test/report are retained. No physical acceptance tolerance was changed.

A separate public regression probe then found that validating the nested tabulated
artifact before the bridge reader's old exception wrapper leaked a DomainError
where the bridge boundary previously reported InputError. The new wrapper now
retains the input-error type, offending field, details and cause. The old boundary
also dropped that field context. The final unchanged historical replay includes
the added test; its non-assertion outcomes are reported separately from assertions.

Current targeted, historical and mutation reports live under `output/assurance/`.
Their source digests and terminal results determine which claims are current.
At the start of this review the full 87-module collection timed out at 30 seconds;
no test module ran in that attempt. Old regression JUnit files are not current
passes. A complete current regression result remains required separately from
the bounded checks above. Startup variability remains an open finding.

## Remaining requirements

The follow-up [power review](POWER_DOCUMENT_AND_NUMERICS_REVIEW.md) adds bounded
electric/throttle document, literal identity and arithmetic checks. It also
closes the explicit empty-driver-list gap in bridge metadata comparison. The
historical checkpoints above retain their original verification scopes.

Other exchange models still require their own complete raw document/metadata
reviews. All function, catalogue and technical-claim review is
still partial. The motor/feed design and mission requirements remain undefined;
all nine physical evidence gates remain unmet. These software checks do not
establish TRL 5, certification or NASA endorsement.

## References

- [RFC 8259](https://www.rfc-editor.org/rfc/rfc8259), sections 4-8: JSON object,
  array, number, string and interoperability background. The closed field sets,
  metadata agreement and default-driver rules here are repository decisions.
- [Python 3.12 fractions](https://docs.python.org/3.12/library/fractions.html):
  binary-input rational arithmetic background; not a source-authentication rule.
