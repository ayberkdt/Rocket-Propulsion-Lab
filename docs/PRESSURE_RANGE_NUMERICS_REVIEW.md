# Pressure correlation, blowdown inverse and feed Jacobian range review

The expanded API helper-field probes exposed sixteen raw-exception reports
in four transient route methods. These were finite-input `OverflowError` or
`ZeroDivisionError` cases, not sixteen independent physical defects. The first
full run is retained under `output/assurance/oct28_api_admission/final`.

## Correctness contract

The equations and hardware-calibration requirements are unchanged. Numerical
admission must preserve a representable final result when an intermediate
power, product or ratio overflows or vanishes. A required positive quantity
that cannot be represented raises `DomainError`; it is not silently made zero.
An arbitrarily small or large finite exponent is not by itself invalid.

* The polytropic inverse compares logarithmic expansion with the empty-tank
  bound before exponentiation. Unreachable cutoffs return the documented zero
  remaining mass. Small expansions use cancellation-safe arithmetic.
* Forward pressure and correlated thrust, Isp and chamber pressure evaluate
  the complete scaled power. Ordinary cases retain the binary64 path. Near-one
  pressure ratios with bounded log contribution use log1p of the relative
  difference, avoiding amplification of a rounded ratio by a large exponent.
  Other large-exponent and range-sensitive cases use 80-digit Decimal log/exp.
  The loose logarithmic bounds only avoid exp of an astronomical argument;
  the final float conversion decides representability. Small log1p/expm1
  arguments below 1e-20 use their leading term (relative error below 5e-21).
* Stream flow, ideal thrust and chamber Isp use exact rational intermediates
  through their final output conversion. The stream scale itself need not fit
  in binary64. Nonzero flow is never discarded by rounding the scale to zero.
* An explicitly supplied maximum pressure of zero is refused; only `None`
  selects the reference-pressure default.
* Feed initial gas masses, recovered ideal-gas pressures and the three
  nonzero injector-pressure state derivatives use complete rational products
  and quotients. A huge ullage square cannot erase a representable Jacobian.

These changes do not make the entire feed energy equation, tank geometry,
performance interpolation or integration arithmetic range-complete. Record
validation, finite-input smoke tests and passing ordinary cases do not prove
that broader claim. The synthetic extreme geometries are numerical probes,
not physically endorsed hardware designs or calibrated operating points.

## Evidence and sensitivity

`test_blowdown_range_numerics.py` independently checks pressure and inverse
identities with 90-digit Decimal calculations. Positive controls cover scaled
powers whose intermediate factor is outside binary64, near-reference ratios
with amplified rounding error, and the finite flow from an overflowing g0*Isp
product. Refusal controls cover genuinely unrepresentable positive outputs;
reference evaluation still succeeds with exponents 1e9 and 1e300.

`test_feed_pressure_range.py` checks the analytic ideal-gas derivatives from
the state inputs, their signs and pressure elasticities, density compensation,
and a finite initial gas mass despite an overflowing pressure-volume product.
The existing open-system energy and independent RK4 conservation oracles remain
part of the regression; those tests were not replaced by finite-JSON checks.

The unchanged-source baseline recorded 19 failed reports, five passing methods
and three passing subtests in these two new modules; pytest subtest reports are
not a count of distinct defects. The initial repaired six-module selection
passed 54 methods and 69 subtests. Final source-bound regression, deliberate
fault controls and package results belong to their individual recorded runs;
an earlier passing build is not evidence for later source edits.

## References and qualification boundary

### Repeated operating-point evaluation

The repaired range arithmetic exposed avoidable repeated work: a source-bound
fixture made 762 built-in performance evaluations for only 257 distinct
pressures. Two real semantic-extreme fuzz methods completed in 234.78 seconds;
that was profiling evidence, not a passing full suite. Both the parallel and
isolated complete fuzz module reached their unchanged 600-second limit.

The exact built-in frozen `PressurePerformanceLaw` now reuses validated
operating points at repeated pressures within one simulation. An LRU retains
at most 2048 entries; each call to the simulator has its own cache and fixed
rated point. Other providers, including subclasses, retain their original
evaluation behavior. Validation still constructs each new pressure point.
No convergence tolerance, refinement budget, field probe or error assertion
is removed to improve runtime.

`test_blowdown_evaluation_reuse.py` checks full sampled curves against
independent analytic primitives. With expelled mass x, n=1 and constant Isp,
the pressure is p0*V0/(V0+x/rho). For F=F0*p/p0,
time(x)=g0*Isp/F0*(x+x*x/(2*rho*V0)); for constant F it is g0*Isp/F0*x.
Impulse is g0*Isp*x in both cases. Controls cover different thrust ratings,
different pressure exponents, custom forwarding providers and subclasses.
The original fixture failed the duplicate-work assertion while its analytic
curves and provider controls passed. Final timings and complete regression
results are retained separately against their actual source snapshots.

### Input ownership behind the frozen model

An accepted float subclass could change its exponent arithmetic after a
`PressurePerformanceLaw` had been constructed: the same rated point and
pressure changed from 50 N to 25 N. A caller-owned sample list likewise
changed an interpolated table point from 75 N to 80 N after row replacement.
Frozen dataclass assignment protection alone did not prevent these aliases.

The tank, correlation and map rows now use the established scalar coercion
contract: accepted `Real` values become native finite floats; booleans, text,
nonfinite values and nonzero inputs erased by float conversion fail with
typed field information. Map source identifiers are validated as text and
the sample iterable is captured as an ordered tuple. Integers and equivalent
real values consequently share canonical scalar serialization. This is an
intentional canonicalization at construction, not a new physical model.

`test_blowdown_record_ownership.py` retains the original mutation examples,
checks physical pressure/flow values, every numeric field's refusal, typed
container/source errors, and a single-pass sample iterator. The full sampled
analytic histories and custom-provider controls remain separate regression
checks. Construction checks do not establish purity of arbitrary callbacks
or prove uninspected model classes safe.

### Direct queries and final-volume range

Constructor normalization did not cover direct method arguments. A boolean
remaining mass was accepted and text leaked a raw `TypeError`. The direct
tank mass/cutoff and correlation/table pressure queries now normalize finite
real arguments with field-specific errors. Fraction inputs and signed-zero
remaining mass retain their valid meanings; nonpositive cutoff pressure is
refused with its pressure field.

Initial liquid, total tank and ullage volume outputs now form their complete
exact rational expressions before one checked output conversion. A positive
volume that would become zero or infinity is refused. Pressure evaluation has
its own range check and does not require an exported intermediate volume:
with p0=1e300, V0=1e200, m0=1e200, rho=1e-200 and n=1, the empty-tank
pressure is about 1e100 Pa even though its volume is outside binary64.
This synthetic scale check expands no physical hardware validity domain.

The complete pressure closure uses relative growth x=(m0-m)/(rho*V0).
It evaluates log(p)=log(p0)-n*log1p(x), with a bounded normal-float path and
an 80-digit fallback. Small growth is retained before a large exponent acts:
V0=1, expelled mass=1e-16, rho=1 and n=1e16 approaches p0/e, rather than
the unchanged p0 produced by rounding V0+dV first. The fallback uses the
leading log1p term only below 1e-20; its omitted terms are below the output
precision when the complete pressure remains representable.

`test_blowdown_direct_queries.py` uses independent 100-digit polytropic
closure and small-growth series, false-zero/nonfinite volume refusal,
representable-pressure rescue, genuine final underflow refusal, direct-query
typed errors and inverse/forward closure. These are implementation equation
checks, not measured pressure-performance calibration.

### Primary and implementation references

- NASA SP-8112, *Pressurization Systems for Liquid Rockets*:
  <https://ntrs.nasa.gov/citations/19760015212>.
- NASA CR-122347, *Monopropellant Hydrazine Resisto Jet*:
  <https://ntrs.nasa.gov/citations/19720011120>.
- Python 3.12, bounded LRU memoization:
  <https://docs.python.org/3.12/library/functools.html#functools.lru_cache>.

The governing equations remain the polytropic gas law, F=g0*Isp*mdot and
p=mRT/V. Source context is NASA SP-8112, NASA CR-122347 and NASA Glenn's thrust
and ideal-gas definitions, already registered in the model documentation.
High-precision calculations verify numerical implementation; they do not
replace test-derived pressure-performance maps, calibration uncertainties,
independent holdout measurements, relevant-environment tests or hardware
scaling evidence. Motor and feed-system TRL 5 remains an open hardware objective.
