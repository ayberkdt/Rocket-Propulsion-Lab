# Exact-integer expansion and response resource review

The uncertainty API previously converted compact decimal exponent text to an
integer before applying study bounds. A short string such as `1e1000000000`
therefore requested an integer with 1,000,000,001 decimal digits. The HTTP
body-size limit does not bound that expansion. Counts and Sobol budgets were
checked only after conversion; the seed was also returned in the JSON result.

## Contract and implementation

The API permits at most **4300 decimal digits** in each of `seed`,
`sample_count` and `sobol_base_samples`. A lower configured host integer/string
limit takes precedence. Disabling that host limit does not disable the API
cap. This is a resource and response-serialization policy, not a floating-point
precision limit or a physical range. Library uncertainty functions retain their
own seed semantics; the changed reader belongs to the JSON request boundary.

The policy matches CPython's default decimal integer/string conversion limit
and prevents compact exponent syntax from bypassing it. The reader never
changes the interpreter's global limit. Existing integers are compared with
`10**limit` without decimal stringification. Finite, integral Decimal values
are checked by their adjusted exponent before `int()` is called. Zero with
an enormous exponent remains zero. Fractional and nonfinite values remain
invalid integers. Sign and notation do not add digits to the value's budget.

An over-budget value raises `InputError` with code
`integer_budget_exceeded`, the offending field and
`details.maximum_decimal_digits`. The API still preserves `1e309`, a 4097-bit
integer seed, and values through `10**4300-1` on a default host. Study-specific
sample bounds still apply separately. A lower host setting may reduce the
accepted seed size; valid values are never coerced through binary64.

## Evidence and limits

An independent Decimal subclass intercepts `__int__` before unsafe allocation.
The first probe observed all three readers requesting the billion-digit
conversion for both signs. No such integer was constructed. A separate
tripwire in `test_api_integer_budget.py` makes the regression attacks safe
even if the production guard is removed. Assertions require refusal before
conversion and before the model runs, correct error code/field, exact boundary
values, host-limit behavior and a real HTTP 400 response. A positive HTTP
control runs the real uncertainty calculation and serializes the large seed.

The first regression baseline retained 23 failed reports, seven passing
reports and nine passing subtests; these are not counts of independent bugs.
The safety tripwire's own 10000-digit threshold is not the application policy.

This check does not establish complete request, memory, execution-time or
output-size coverage for all routes. In particular, passing integer admission
tests does not prove that all semantically invalid thermal controls are
rejected before costly chemistry or that physical model validity is complete.
Source-bound full-suite and distribution evidence remains separate.

## References

- Python 3.12, integer/string conversion limits, default policy and complexity:
  <https://docs.python.org/3.12/library/stdtypes.html#integer-string-conversion-length-limitation>.
- Python 3.12, Decimal adjusted exponent:
  <https://docs.python.org/3.12/library/decimal.html#decimal.Decimal.adjusted>.
