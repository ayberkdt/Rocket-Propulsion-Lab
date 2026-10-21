# Electric/throttle documents, literal metadata and arithmetic review

This is a bounded continuation of the exchange review. It covers the electric
throttle table, discrete/fixed-efficiency model documents, throttleable model
documents, their owned literal metadata, electric level ratios, continuous
operating points and electric model output arithmetic. It does not establish
complete model adequacy, measured performance, calibration, hardware readiness
or exhaustive correctness of the remaining functions and catalogue claims.

## Reproduced defects and document contract

The prior readers discarded extra root, law, table-level, parameter, constraint
and throttle-curve fields. A fixed-efficiency document could silently discard an
entire supplied throttle table. Model reference identifiers and electric force,
mass and power semantics could disagree with the rebuilt manifest. A mapping
used as a stream pair could have its keys converted into physical components,
discarding its values. Empty explicit parameter lists could become invented
default driver sets in bridge, electric and throttleable readers.

The raw readers now require closed built-in JSON object/array/number/text fields.
The law discriminator selects an exact field set: discrete laws require their
complete table, while fixed-efficiency laws reject a supplied table or fields
from the other law. Parameters, constraints, levels, stream pairs and curves are
checked before reconstruction. Nullable fields remain nullable as emitted by
the writers. Fixed implementation metadata and explicit parameter definitions
must agree with the reconstructed model, in addition to existing hashes and
numerical-policy checks. An empty driver array cannot acquire new parameters.

This is not byte authentication: object key order and equal integer/float
representations may normalize, and valid changed models with their own content
hashes are accepted. Unknown extensions require an explicit schema contract.
Other persistence families still require their own reviews.

## Literal identity and corrected claims

Owned level, table, law, path and model labels normalize to plain literal strings,
preserving characters and whitespace. Schema comparisons and digest iteration
cannot be redirected by string-subclass callbacks. The scalar helpers already
validate Unicode; constructors now retain their returned literal values.

The source formerly called a supplied electric level qualified and asserted that
table values came from qualification/acceptance campaigns. Constructors cannot
establish either claim. Those statements now describe supplied values and source
metadata. A digest is a declared identifier, not proof that a campaign occurred
or that the data apply to a selected motor and mission.

## Exact arithmetic and checked outputs

Electric level efficiency uses the exact ratio of supplied binary scalars before
the consistency decision. The existing ceiling is the represented float value
of `1.0 + 1e-12`; it has not been enlarged or reinterpreted as a measurement
uncertainty. Decimal reference calculations check neighboring candidate thrusts
on both sides of that ceiling across wide binary exponents. Exact common scaling
also preserves the analytic efficiency of one half and constant specific impulse.

Efficiency and Isp exports are checked for finite, nonzero representation. A tiny
valid efficiency need not be exported during construction. Likewise, continuous
thrust/power gain remains exact inside an operating point: an unrepresentable
gain can recover to representable thrust and flow after multiplication by power
at saturation. Mandatory unrepresentable outputs, including unsaturated slopes,
fail with a typed domain error rather than becoming zero or infinity.

Continuous operating points validate public mode/power arguments. Finite negative
power still clips to zero for a held-on mode's integrator probes; both local
derivatives are zero below that clamp. At exactly zero the implementation uses
the right derivative; this is a numerical convention at a nonsmooth point.
Discrete operating points reject negative Python indices instead of selecting
the last level, and refuse malformed modes or non-finite power arguments.

Electric model duty/flow products, tank drain, power consumption, deficit, Isp
and flow partials retain exact intermediates until checked record exports. Tank
drain cannot silently vanish while a positive commanded thrust remains active.
Built-in laws share a private exact point ledger with the model, so a tiny
unscaled point slope is not exported before a driver can recover it. Public
operating-point calls still require their own complete float record to be
representable. Existing subclass callbacks retain their public point boundary;
arbitrary extension behavior and its physical validity are not assessed here.
These are representation/consistency guarantees for synthetic input cases, not
claims about physically realizable extreme motors or an empirical validity range.

## Adversarial evidence and unchanged workload

Public document mutations retain original hashes to distinguish ignored data
from ordinary content-hash rejection. The tests include both law variants,
nullable controls, reopening behavior/events, direct programmatic metadata and
independent Decimal arithmetic/scaling oracles. Final tests are replayed unchanged
against verified pre-review package bytes; assertion failures and unrelated
exception outcomes are recorded separately and may overlap within one method.
Historical method failures are not a count of independent defects.

Eighteen curated variants target specific raw-schema, metadata, literal-identity,
decision, index, gain, slope and drain failures. The mutation runner now preserves
property getter/setter/deleter descriptors while injecting getter variants. Only
assertion failures without unrelated errors count as kills.

A recorded full-suite timeout snapshot showed the JSON attack generator copying
entire documents. Its edit helpers now copy only the ancestors of one change;
readers receive serialized JSON strings, so shared untouched branches do not
escape. All 1,829 generated labels and JSON texts were compared byte-for-byte
against the archived generator, with identical random states and unchanged source
documents. All other module AST nodes, including fixtures, mutation count, test
selection and assertions, are unchanged. Observed generator timing is diagnostic;
it does not establish concurrent-load latency or resolve all earlier timeouts.

The full regression retains four workers, a 30-second collection budget and a
120-second budget per module. Final terminal reports and source bindings determine
whether it passes; earlier failed/cancelled attempts remain evidence of incomplete
verification. No test selection, workload or acceptance tolerance was reduced.

## Remaining requirements

Mode selection/hysteresis boundary policies and the rest of the model, catalogue,
native integration and technical-claim inventory still need complete review.
Independent hardware holdout applicability and physical uncertainty remain open.
Mission and motor/feed configuration are undefined and all nine hardware evidence
gates remain unmet. Software verification alone does not establish TRL 5,
certification or NASA endorsement.

## References

- [NASA/CR-2004-213099](https://ntrs.nasa.gov/api/citations/20040084644/downloads/20040084644.pdf),
  printed page 14, equations (2-9)-(2-10): jet power and system input power define
  total system efficiency. Algebra gives the ideal constant-efficiency law used
  here; this does not validate a supplied throttle map or select hardware.
- [NASA Glenn specific impulse](https://www1.grc.nasa.gov/beginners-guide-to-aeronautics/specific-impulse/):
  the thrust, propellant mass-flow and standard-gravity ratio defining Isp.
- [Python 3.12 fractions](https://docs.python.org/3.12/library/fractions.html):
  rational preservation of binary inputs. The retained ceiling, output-range
  checks, clamp derivative convention and test budgets are repository policies.
- [RFC 8259](https://www.rfc-editor.org/rfc/rfc8259): JSON type/interoperability
  background; it does not mandate this application's closed model schemas.
