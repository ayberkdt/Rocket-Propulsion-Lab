# Technology readiness: what the software shows and what TRL 5 still needs

This page maps the thermodynamics and aerodynamics evidence in this
repository to NASA's TRL definitions (NPR 7123.1D, Appendix E). It does not
assign a TRL to any hardware. A passing test suite is software evidence; it
is not a demonstration in a relevant environment.

## The TRL 5 criteria

NPR 7123.1D defines TRL 5 as a component or brassboard validated in a
relevant environment. For software it asks for end-to-end software elements
implemented and interfaced with existing systems or simulations, and tested in
a relevant environment against predicted performance. Its success criterion is
documented test performance that agrees with analytical predictions, plus
documented scaling requirements.

## Evidence that exists (software, independent references)

- Normal-shock arithmetic retains finite complete perfect-gas outputs and
  checks pressure recovery against independent mass, momentum and enthalpy
  relations (`NORMAL_SHOCK_RANGE_REVIEW.md`). This is numerical evidence for
  the gas model.

- Prandtl–Meyer forward and inverse arithmetic has differential-quadrature and
  exact-supplied-angle controls (`PRANDTL_MEYER_REVIEW.md`). Synthetic gamma
  ranges and numerical bracket screening do not establish a gas-model operating
  envelope or relevant-environment motor/feed hardware validation.

- Complete area/Mach-angle inverse arithmetic has independent entropy/mass-flux
  and trigonometric controls (`AREA_ANGLE_INVERSE_REVIEW.md`). Its numerical
  bracket and synthetic gamma matrix do not establish an operating envelope,
  continuous physical flow or relevant-environment hardware validation.

- Static/total and static/sonic scalar ratio inverses have exact-input
  entropy/energy controls and logarithmic bracket screening
  (`ISENTROPIC_RATIO_INVERSE_REVIEW.md`). This does not establish all inverse
  algorithms, a real-gas operating domain or hardware readiness.

- Complete forward isentropic ratios retain balanced extreme outputs and sonic
  normalization, with entropy/continuity oracles and checked positive exports
  (`ISENTROPIC_FORWARD_RANGE_REVIEW.md`). Synthetic arithmetic coverage does not
  establish perfect-gas applicability, physical-flow accuracy or hardware TRL.

- Near-sonic and subnormal geometry has complete-expression, coordinate-offset,
  high-precision and consumer controls (`GEOMETRY_PRECISION_REVIEW.md`). Synthetic
  ranges and a discrete planar-net comparison do not establish an experimentally
  validated flow or hardware environment.

- Geometry generators own physical inputs and discrete resolutions before flow
  work, preserve the converging-only sonic limit, and refuse invalid sampled
  output (`GEOMETRY_GENERATOR_REVIEW.md`). Independent equation, scale, export
  and real-cooling controls do not qualify manufactured geometry or hardware.

- Cooling contour admission now shares the prescribed Bartz geometry snapshot,
  with ordering/coherence and caller-mutation controls before EOS/axial work
  (`REGEN_CONTOUR_CONTRACT_REVIEW.md`). Conservation and generator regressions
  do not qualify a continuous wall, a prescribed flow field or physical hardware.

- Cooling controls and material use retain strict physical scalar admission,
  ownership before axial iteration and the actual conductivity/temperature-source
  conditions in results (`REGEN_SOLVER_ADMISSION_REVIEW.md`). A temperature
  difference to a conditional material criterion is not a structural or lifetime
  acceptance result.

- Frozen gas/cooling boundary: malformed and impossible states, caller-owned
  arithmetic, lost kinetic contributions and balanced extreme Bartz expressions
  receive independent numerical and coupled controls
  (`HOT_GAS_CONTRACT_REVIEW.md`). These synthetic range checks do not establish
  a relevant operating environment or certify thermal hardware margins.

- Pressure-fed numerical range: independent high-precision power, inverse,
  stream-flow and ideal-gas Jacobian checks preserve finite outputs through
  overflowing intermediates and refuse unrepresentable positive results
  (`PRESSURE_RANGE_NUMERICS_REVIEW.md`). Synthetic scale attacks are not
  a hardware performance envelope or relevant-environment validation.

- Dynamic liquid line balance and pressure cutoff: independent high-precision
  and closed-form checks cover signed near-equilibrium acceleration, the local
  Jacobian, mass closure against the exported pressure rate, and the stopping
  direction. An experimental water-line comparison in NASA CR-124308 limits
  the first-mode surrogate's usable bandwidth. See `DYNAMIC_FEED_LINE.md`.

- Combustion equilibrium and rocket performance: nine NASA CEA 3.3.4 points
  and the 0.5–25 MPa domain reproduced within 2e-5; RP-1311 example 8 at its
  published station (`VALIDATION_REPORT.md`).
- Transport properties and equilibrium derivatives: NASA CEA 3.3.4 transport
  output reproduced (species to 1e-15, rocket mixtures to 4e-4)
  (`TRANSPORT_AND_DERIVATIVES.md`).
- Real-fluid propellant and pressurant properties: H2, O2, N2, CH4 and He
  reference equations agree with NIST WebBook tables within 11 ppm
  (`REAL_FLUIDS.md`).
- Compressible flow: NACA Report 1135 tables and charts, NASA SP-3004 cone
  tables, Anderson's worked examples (`COMPRESSIBLE_ADVANCED.md`).
- Subsonic aerodynamics: exact potential-flow solutions, Anderson's worked
  examples and NACA Report 824 section data (`AERODYNAMICS_SUBSONIC.md`).
- Viscous flow and heating: Blasius, NASA TN D-6945 (Van Driest II), NASA
  TP-2914 and the 1976 Standard Atmosphere (`AEROTHERMODYNAMICS.md`).
- Launch-vehicle aerodynamics: NACA TR 1048 and TR 1307 data, Barrowman's
  method and NASA SP-8035 (`LAUNCH_VEHICLE_AERODYNAMICS.md`).
- Interfaces to existing simulations: the propagation bridge flies in
  Sidera's public propagator and in Lunaris (`SIDERA_NATIVE_WRAPPER.md`,
  `LUNARIS_INTEGRATION.md`).
- The optional SciPy reference trajectory closes a regulated single-tank
  feed against a motor pressure map, integrates translation and checks
  protective roots, command jumps and subsequent coast against independent
  analytic cases. This is an end-to-end software check for the stated model.
  Native Sidera additional-state support and validation against measured
  motor and feed responses remain separate requirements. See
  `REGULATED_FEED_DYNAMICS.md`.

- Independent-oracle tests (`ADVERSARIAL_TESTING.md`): compressible
  relations, conduction, real-fluid EOS derivatives, Gibbs equilibrium,
  kinetics, radiation view factors (Monte Carlo ray tracing), rocket
  performance identities, vortex panels, global sensitivity (Ishigami and
  rocket-equation closed forms) and the 2-D channel wall (an independent
  `scipy.sparse` solve of the same cells agrees to 1e-9); each fix was
  mutation-checked by re-introducing the defect.
- API numeric admission: declared top-level scalar numbers, including shared
  helper readers, must refuse malformed values before model execution with
  the correct field. Direct nonnumeric getters and numeric semantic extremes
  receive typed-error/strict-JSON smoke checks. Nested schemas and combinations
  of extreme fields are not comprehensively covered
  (`API_NUMERIC_ADMISSION_REVIEW.md`).
- Hot-fire comparisons of the regenerative-cooling gas side and heat load
  (NASA TN D-6396, TN D-3532; `REGENERATIVE_COOLING.md`).

Several comparisons are against wind-tunnel or flight-derived data (NACA TR
824 section lift, NACA TR 1048 body normal force). Where the model disagrees
with measurement the disagreement is reported, not tuned away: inviscid
panel lift is about 14 % above measured section lift, and the Allen–Perkins
normal force is 22–24 % above the TR 1048 measured point.

The NASA ICPTA LOX/LCH4 test report supplies a separate empirical feed-model
challenge. Three of eight ambient-thermal shutdown peaks in Table 1 exceed
the corresponding nominal instantaneous-closure estimates. The largest
reported pair is 91 versus 72 psid. Individual numerical uncertainties
are absent from the table. Cold testing and opening transients used distinct
conditions. The signed comparison is recorded under
`output/assurance/icpta_feed_benchmark_20261110` and does not identify a
single-pipe model or qualify the user's motor/feed configuration.

The Ares I upper-stage cold-flow article provides an independent check on
spatial resolution. Table 3.2 gives ten local peaks per test for nine
valve closures. The published EASY5 model locates the highest peak at a
different sensor in seven tests, although the absolute error in maximum
pressure anywhere in each run is at most 6.72%. Assigning one pressure to
every sensor would leave a worst local relative error of at least 22.93%
to 32.42%, depending on the run. The calculation is under
`output/assurance/ares_upper_stage_surge_benchmark_20261111`.
[Williams, Table 3.2, NASA NTRS](https://ntrs.nasa.gov/citations/20110007275)

For WH5b, the reported run maximum differs by 25 psia while the local
HfP-310 difference is 106 psia. This is a replay of published EASY5 outputs
against the article's measurements. The author identifies WH6r as the
baseline. The repeat-run statistics in §3.3 do not supply absolute error
bounds for each observation and prediction. The replay provides a
counterexample to relying on the run maximum alone.

Regulator catalogue limits need pressure-reference and configuration
traceability. For example, Moog lists a 260 psia single-stage lock-up for
the 250 psia family. The separate 50X1463 table lists 277 psig nominal for a
failed stage and 289 psia maximum lock-up without giving the applicable
stage for that limit. The failed-stage nominal converts to 291.70 psia at
standard atmospheric reference. The two conditions need supplier
clarification before they can define one envelope. The public NASA Robotic
Lunar Lander test abstract gives a different regulator's measured
performance bounds without run traces or uncertainty. Both sources are
audited under `output/assurance/regulator_catalogue_audit_20261111`.
[Moog variant table](https://www.moog.com/products/propulsion/space-propulsion/spacecraft-propulsion/propulsion-components/helium-regulator-pneumatic-pressure-regulator.html)
[NASA Robotic Lunar Lander regulator test](https://ntrs.nasa.gov/citations/20130014156)

## Thruster catalogue review

- Moog's 058E151 column lists 120 mN at a nominal 6.9 bar absolute nitrogen
  inlet and a nominal Isp of 65 s. The opening power ceiling and the holding
  power belong to distinct electrical phases. The mass covers the thruster,
  rather than the complete feed system.
  [Moog cold gas thrusters, 058E151 column](https://www.moog.com/content/dam/moog/literature/sdg/space/propulsion/moog-coldgasthrusters-datasheet.pdf)
- ArianeGroup lists 1 N, 220 s and 0.44 g/s in its 1 N nominal rows. The
  thrust and Isp imply about 0.4635 g/s, 5.34% above the nominal flow row
  and 3.69% above the upper end of the separately published 0.142–0.447 g/s
  range. The public tables leave the operating conditions and cause unresolved.
  The 10 N family gives 10 N, 292 s and 3.50 g/s. Its calculated flow is
  about 3.4922 g/s, inside the published 2.30–4.20 g/s interval. The S10-13
  catalogue record uses the 350 g single-seat valve mass. The 650 g dual-seat
  version is a different configuration. The pressure reference for the
  published inlet ranges is unspecified in the cited tables.
  [ArianeGroup 1 N specifications](https://www.space-propulsion.com/spacecraft-propulsion/hydrazine-thrusters/1n-hydrazine-thruster.html)
  [ArianeGroup 10 N specifications](https://www.space-propulsion.com/spacecraft-propulsion/bipropellant-thrusters/10-bipropellant-thrusters.html)
  [ArianeGroup portfolio, p. 6, both flow ranges](https://www.ariane.group/app/uploads/2024/06/Orbital_Propulsion_2019_07_PB_EN_Web.pdf)
- The ECAPS NASA survey record retains 204–235 s and an unscoped
  `survey_power` field. The separate Bradford record retains 204–231 s
  typical in steady vacuum operation and labels 8–10 W as regulated
  reactor preheat. Neither record supplies a paired thrust/Isp point or a
  complete firing bus power figure.
  [NASA in-space propulsion survey](https://www.nasa.gov/smallsat-institute/sst-soa/in-space_propulsion/)
  [Bradford ECAPS manufacturer sheet](https://satcatalog.s3.amazonaws.com/components/860/SatCatalog_-_Bradford_Space_-_1N_HPGP_-_Datasheet.pdf)
- Busek's August 2021 BIT-3 data sheet and product page agree on
  56–75 W, upper thrust of 1.1 mN, upper Isp of 2150 s and 1.40 kg with
  the gimbal. The RF family page gives 56–80 W, 1.25 mN and 2300 s.
  The catalogue uses the dated sheet for the specifications and records
  the unresolved difference in the family page.
  [Busek data sheet](https://www.busek.com/s/BIT3_v10.pdf)
  [Busek RF family page](https://www.busek.com/rf-ion-thrusters)

Sizing requires the thrust and Isp fields to share a labelled nominal
point and a source. Independent ranges remain unsized. The project
compares the calculated nominal mass flow with the reported range only
when both endpoints and the nominal pair share a source. The calculation
uses the exact rational values of the encoded inputs before rounding the
flow for display. A discrepancy keeps the preliminary budget available
and withholds the weighted ranking pending source review. A value inside
the range only establishes agreement with those published endpoints.
It supplies no pressure-dependent thrust, flow or efficiency map.
The project retains the catalogue snapshot and the separate nominal-flow
and range comparisons. The manufacturer nominal points support preliminary studies
with the documented limitations. Relevant-environment data are needed
for the integrated motor and feed system.

### Thruster inlet pressure in a preliminary study

A candidate may supply `inlet_pressures` for the same operating condition
at its thruster inlet ports. The pressure bounds use Pa and explicitly
declare an absolute, gauge or differential reference. A gauge window can
be converted with the lower and upper absolute reference pressures
and their review record. The reference bounds must describe the same
window. The code adds corresponding lower and upper bounds before
rounding for reporting. A gauge window without that review and a
differential window remain unassessed for this comparison.

The Moog 058E151 catalogue entry ties 6.9 bar absolute to its nominal
performance point. A variable window that includes that point leaves
other pressures outside its scope. The ArianeGroup 1 N table lists
5.5–22 bar without a pressure datum. The 10 N family table lists
10–23 bar without identifying the datum or providing separate inlet
ranges for fuel and oxidizer. Those catalogue values remain visible
with their limitations. If a bipropellant candidate supplies an inlet
window, the comparison lists both propellant ports.

An unresolved supplied pressure window withholds the weighted score and
retains the preliminary mass and time budgets. A matching nominal point
or coverage by a published interval is a catalogue screen. Relevant
operating conditions, dynamic response and component-to-system
applicability still require independent evidence.



The Moog 058E151 source also lists 27.6 bar absolute MEOP, 70.0 bar
proof pressure and 111.4 bar burst pressure. The catalogue stores all
three component ratings and uses only MEOP for the operating comparison.
The MEOP row prints 27.6 bar and 400 psia, which convert to 2,760,000 Pa
and approximately 2,757,902.917 Pa, respectively. The comparison retains
both printed values and uses the exact unit conversion at the boundary.
A window must be within both expressions to pass this screen. A window
entirely above both is an exceedance. A window above the lower expression
and entirely within the higher expression records a unit disagreement.
Other windows straddling either expression record a crossing. All three
unresolved outcomes withhold ranking independently of the nominal
performance screen. The span between the two printed expressions is
not a measured tolerance or a system pressure limit.

The component MEOP and the proof and burst test levels have different
uses. Below-MEOP pressure alone supplies no performance map or proof
that the regulator, lines and other components share that rating. System
evaluation needs the selected configuration and the pressure excursions
under relevant operating conditions.
[Moog cold-gas thruster specifications](https://www.moog.com/products/propulsion/space-propulsion/spacecraft-propulsion/thrusters/cold-gas-thruster.html)
[NBSIR 82-2401, Appendix 4, exact psi-to-Pa conversion](https://nvlpubs.nist.gov/nistpubs/Legacy/IR/nbsir82-2401.pdf)
[NASA-STD-7012A definitions](https://standards.nasa.gov/sites/default/files/standards/NASA/A/0/2023-02-22-NASA-STD-7012A-Approved.pdf)

The project compares the ArianeGroup 1 N references of 12 h for one
continuous burn, 50 h accumulated, 67 kg throughput and 135,000 N s
total impulse. The 10 N family publishes 8 h and 69 h for the two time
references. Valve and operating-profile applicability must be confirmed
for a specific S10-13 article. The tool records published cycle counts
without inferring them from the mission budget.

For accumulated firing, the lower bound on the maximum per-head demand is
mission on-time multiplied by the simultaneous firing count and divided
by the installed count. Consumed propellant and total scalar impulse are
already fleet totals. Each is divided by the installed count and compared
with the corresponding published value. These lower bounds allow ideal
sharing and can reveal a demand that exceeds a reference even with all
heads available. Below-reference results leave the real allocation open.
The single-burn check only uses ungrouped continuous occurrences with all
installed heads firing. A definite exceedance keeps the preliminary
budgets visible and removes the alternative from the weighted ranking
pending review. Catalogue references do not supply a qualification
envelope for the specific mission.
[ArianeGroup 1 N qualification table](https://www.space-propulsion.com/spacecraft-propulsion/hydrazine-thrusters/1n-hydrazine-thruster.html)
[ArianeGroup 10 N qualification table](https://www.space-propulsion.com/spacecraft-propulsion/bipropellant-thrusters/10-bipropellant-thrusters.html)

A candidate can supply an optional `command_history` in the SASS project
request. The same review is available directly through
`review_thruster_command_life(record, history)`. A history names its
source record, clock period in seconds, duration in integer ticks, all
installed head IDs and the heads already open just before tick zero.
Commands have exactly three fields: `tick`, `head_id` and the boolean
`opened`. Supply state transitions in time order. Heads at a common
tick change together. Commands at the final tick can close an observed
cycle. The reported peak counts intervals of positive duration.

The ledger records each head's observed open duration, longest open
interval and completed open-close commands. An interval open at either
window boundary is censored. Its observed duration and completed cycles
are lower bounds on exposure in the wider history. The duration comparison
uses the adjacent binary representations of the declared clock period
near a catalogue boundary. A clear exceedance of a published reference
or peak concurrency above the sizing assumption removes the preliminary
score pending review. The command window is supplementary evidence and
does not establish mission coverage, valve motion or physical cycle
equivalence. Supply telemetry, the full duty profile and the manufacturer's
qualification conditions for a hardware assessment.
[NASA JSC-67723, §§3.1.7–3.1.8](https://standards.nasa.gov/system/files/tmp/JSC-67723_Propulsion_Standard_Basic_Final%20%281%29_0.pdf)

For instrumented tests, `review_valve_responses(commanded, observed, limits)`
compares issued valve transitions with observed valve-state transitions.
Supply the same head IDs, common clock period and window in both histories.
`ObservedValveTrace` records which heads have reviewed continuous coverage,
the clock alignment record and a worst-case relative timestamp error in
ticks. `ValveLatencyLimits` records separately approved opening and
closing delays on that clock. The reviewer establishes coverage and
timing uncertainty from the instrument records.

Each observation can belong to a single command interval for its head.
An observation whose timing bound crosses an interval boundary remains
ambiguous. The tool records the lower and upper latency for a matched
transition. An upper bound within the approved delay supports a pass for
that transition. A lower bound above the delay supports failure. An
unobserved response supports failure only when the full deadline can be
observed before a superseding command. Incomplete coverage, ambiguous
prehistory and censored boundaries remain visible in the result.
This comparison addresses observed valve states and timing. Pressure
transients, structural response, chamber decay and thrust require their
own synchronized measurement channels.
[NASA NESC transient pressure guidance](https://ntrs.nasa.gov/citations/20210022275)

NASA Ames's discrete pneumatic-valve testbed paper gives opening and
closing requirements of 6 s and 3 s for its own ground-support hardware.
Its nominal cycle description gives approximate completion times, and its
prognostic comparisons use simulation. The related NASA Open Data records
DASHLINK_947 and DASHLINK_948 list publications as downloadable resources.
These sources motivate separate opening and closing limits. They do not
provide the synchronized command, valve-position, clock and coverage
records needed to run this assessment on the selected motor and feed
hardware.
[NASA Ames testbed paper, NTRS 20140008531](https://ntrs.nasa.gov/citations/20140008531)
[DASHLINK_947](https://data.nasa.gov/dataset/application-of-model-based-prognostics-to-pneumatic-valves-in-a-cryogenic-propellant-loadi)
[DASHLINK_948](https://data.nasa.gov/dataset/a-testbed-for-implementing-prognostic-methodologies-on-cryogenic-propellant-loading-system)

## Reducing measured pressure and shutdown thrust

Use `review_transient(trace, start_tick, end_tick, window_review_record=...)`
on each `TransientTrace` in Pa or N. Record the raw data identifier, channel
location or thrust direction, calibration, timestamp alignment, usable
bandwidth and accepted maximum gap. The pressure datum must match the model.
Check saturation, clipping, dropped records and the effect of signal filtering
before supplying a sample. The code uses an exact common tick period. Include
any timing error in the per-sample value bounds or resolve it in the acquisition
review. The stated absolute errors must hold simultaneously for the window.

Use a separate justification for the maximum physical slew in Pa/s or N/s.
Sampling frequency and observed sample-to-sample differences alone leave
between-sample extrema undetermined. The reducer propagates admissible
sample intervals across the window, bounds the physical signal between
samples and integrates the linear envelopes. It publishes a peak interval
and a signed integral. For a thrust channel it also publishes the integral
of positive thrust in N·s, which can remain nonzero when signed impulse
cancels. This last result applies to the declared direction. To examine
residual thrust in the reverse direction, use a separate reviewed channel
or reduce a trace with the sign convention reversed.

A missing boundary, oversized gap, missing slew or globally inconsistent
samples produces a specific status and no numeric enclosure. A bounded
result applies only within the selected window and at the named channel.
For a shutdown window, use synchronized closure observations and bracket
the event timing uncertainty through the observed end of the tail. Keep
the raw trace, calibration and independent window review with any result
transcribed into the holdout plan. Pressure envelopes at individual
sensors do not establish the spatial load history or structural response.

The NASA ICPTA publication distinguishes the 25 kHz dynamic manifold
pressure acquisition from the chamber-pressure instrumentation. The authors
state that they could not characterize the influence of the initial transient
on chamber-pressure oscillations with the installed sensors. The NTRS
entry provides the paper and figures. It provides no synchronized sample
files with the calibration, clock records and physical slew evidence needed
to exercise this reduction for that article.
[NASA ICPTA instrumentation and limitations, §§ III.A and IV](https://ntrs.nasa.gov/citations/20170005403)
[NASA NESC transient pressure guidance](https://ntrs.nasa.gov/citations/20210022275)

For an independently reviewed holdout target, attach a
`TransientValidationTarget` to the plan's `ValidationTarget`. Its channel
ID, exact start and end ticks, metric and window review record enter the
sample digest through
`validation_sample_sha256(evidence, sample, transient_target=recipe)`.
Identify the physical test run using the same condition record on the
plan target, the sample and the archived acquisition channel.

Group the reviewed traces in an `ArchivedTransientRun` and write
`transient_archive_bytes(run)` to the evidence's `raw_data_path`. This is
the interpreted acquisition export. Preserve the native DAQ files and the
reviewed conversion record as well. The archive retains the article ID,
run/channel pairing, source and acquisition review records, clock and
slew declarations and every sample value, error and tick. Binary64 values
are stored as hex strings. Use `transient_archive_from_bytes(data)` to
reopen the canonical export in a separate process.
In `assess_validation(evidence, plan=plan, transient_run=run)`, the on-disk
file and the reconstructed typed export must both match the declared
raw-data digest.

Use `transient_measurement(trace, recipe)` to obtain the exact submitted
binary64 center and outward absolute bound for the selected metric.
The calculation retains the original enclosure in the replay report.
The symmetric representation can widen the enclosure by rounding. For
a singleton, the smallest positive binary64 number is the radius
required by the positive-uncertainty validation interface. Peak uses
Pa or N. Signed integral uses Pa*s or N*s. Positive impulse applies
only to the N channel and uses N*s. Replay checks both the submitted
center and radius against the archive, even when the submitted row
already matches its independently reviewed digest. Missing acquisition
data or an assumed clock period or slew leaves the target unassessed.
The digest confirms record identity. The acquisition and independence
reviews establish whether those records represent the physical test.

### Calibration and holdout acquisition partition

For a reviewed model comparison with a transient target, attach an
`AcquisitionPartition` to the validation plan. Map each calibration
sample ID to its calibration condition. For every calibration and
holdout condition, list the contributing `NativeAcquisitionSlice`
records. Each slice identifies the physical firing, the SHA-256 of
the original acquisition file and its half-open range of native
record indices. Use indices in the original file's coordinate system
before trimming, clock rebasing or resampling. Use the same physical
firing ID across independently exported files and instrument streams.

Review the complete map against the native DAQ records and the actual
calibration inputs. Store the independently approved
`acquisition_partition_sha256(evidence, plan)` on the plan as
`reviewed_acquisition_partition_sha256`. It binds the map to the
article, interpreted acquisition bytes, model configuration and
calibration and holdout identities. The checker flags a shared firing
even if its files differ, or overlapping indices in the same native
file even if the firing IDs differ. Distinct firings may use adjacent
segments in a continuous DAQ file. A missing, incomplete or conflicting
partition withholds model agreement for affected evidence. Physical
acceptance still uses the independently replayed measurement.

Attach one `NativeAcquisitionFile` per SHA-256 referenced in the approved
map using the plan's `native_acquisitions` field. Each binding supplies
the expected digest and a local `Path` to the original file. The checker
hashes the complete file and reports missing, unreadable, unplanned or
modified sources under `native_acquisition_review`. A byte-identical
file may be relocated without changing the approved partition digest.
The interpreted export, model or prediction file cannot serve as
original DAQ evidence. An empty file cannot support a nonempty slice.

The reviewer substantiates the firing IDs, source indices and calibration
history against the native files. Byte verification establishes that the
reviewed source files are present.
[NASA-STD-7009B, §§4.2.1 and 4.2.5–4.2.6](https://standards.nasa.gov/standard/NASA/NASA-STD-7009)

### Native engineering-unit CSV replay

Use `NativeCsvTrace` in the partition's `native_csv_traces` field when the
original DAQ file stores force or pressure directly in engineering units
as ASCII CSV. The file has a unique header, a `tick` column of exact
integers and a named measurement column of finite decimal numbers.
Each mapping names the condition, interpreted channel, native file SHA-256,
measurement column, unit, half-open data-row range and integer tick offset.
The header is excluded from row indexing. The mapping also carries the
conversion review record and enters the partition approval digest.

For a declared conversion the checker requires a mapping for every planned
transient channel. Each row range must lie within a reviewed native slice
for the same condition and file. Replay compares every selected tick after
applying the offset and every binary64 sample center. The parser hashes
the bytes it reads, including rows outside the selected range, and checks
that digest against the approved source. A conversion mismatch withholds
model agreement and any physical acceptance that depends on the affected
force or pressure measurement. Valve timing retains its separate evidence.

An optional `absolute_error_column` names a symmetric native error bound in
the measurement unit. The approval digest binds its column name. Replay
requires the interpreted error to be at least the native error plus the
absolute difference between the exact decimal center and the binary64 center.
It compares rational values so rounding either native number cannot shrink
the interval. For this exact comparison, each decimal token is limited to
1,024 characters, with nonzero values restricted to adjusted decimal
exponents from −1075 to 308. Invalid, negative or narrowed errors cannot
support measurement acceptance.

An optional `tick_period_s` binds the reviewed native clock conversion in
seconds per tick to the partition approval digest. Replay requires the
same binary64 period on the interpreted trace. A mismatch withholds results
that depend on that channel even when the tick indices and values match.
The review record must identify the DAQ configuration or clock calibration
from which the period was derived. Clock drift and timing uncertainty
still require independent evidence.

Omitting an optional field retains the earlier digest representation.
Sensor placement, calibration, joint uncertainty coverage and slew remain
reviewed inputs. Other DAQ formats and valve histories still rely on
reviewed conversion records. When no CSV mapping is supplied, the native
trace review reports `not_declared`.

### Fixed-window model predictions

Use a reviewed `rocket-propulsion-prediction-run-v5` export for
`TransientValidationTarget` model comparisons. Each v5 row has
`pressure_referent`, `shutdown_scope` and `fixed_window` keys. The
applicable fixed window supplies `channel_id`, exact integer
`start_tick` and `end_tick`, `tick_period_s`, `window_review_record`
and `metric`. For example, forward impulse uses `positive_impulse`
and `N*s`, while signed impulse uses `signed_integral` and `N*s`.
The prediction scope distinguishes these statistics even when they
give the same number for a particular firing.

Use `null` for inapplicable scopes. A row can identify a fixed window
or a closure-relative shutdown. The v5 `fixed_window` replaces the
older `pressure_window` key and includes the acquisition channel ID.
Replay checks the scope and acquisition clock using the physical test
condition and channel as a pair. Earlier prediction archives remain
readable. Fixed-window model agreement stays unassessed when the
prediction lacks the channel and window binding. A valid independent
measurement remains available for physical acceptance.

### Pressure identity in acquisition and prediction archives

A `PressureReferent` names the physical pressure location and the reference.
Absolute pressure has no second location. Gauge pressure identifies its
reference environment. Differential pressure identifies the second location
and uses the order `location minus reference_location`. The validation
path compares these identities without applying an offset or a sign change.

Set `pressure_referent` on the Pa acquisition channel and its
`TransientValidationTarget`. The canonical acquisition export uses
`rocket-propulsion-transient-run-v3` when a pressure referent is present.
The independently reviewed sample digest includes the target referent.
If the acquired and planned referents differ, the pressure measurement
cannot enter the model comparison or physical acceptance.

A pressure target with an explicit referent requires a matching entry
from a reviewed v5 prediction archive. The row supplies the
`pressure_referent` with `location`, `reference` and `reference_location`
and the `fixed_window` above. The permitted metric and unit pairs are
`peak` with Pa and `signed_integral` with Pa*s.

The alignment review identifies the common time origin and alignment
of the prediction and acquisition windows. Replay checks the exact
integer ticks, the metric, the review record and the clock period from
the archived acquisition channel. A v2 archive can identify the
pressure referent. The v3 and v4 archives add the metric and time
window through `pressure_window`. Their pressure prediction metadata
has no acquisition channel ID.

Legacy numerical replays remain reproducible. An undeclared Pa channel
is marked `legacy_not_declared` and supplies no evidence of matching
pressure location or reference. The acquisition review must substantiate
that the archived channel was measured at the declared physical location.
[NASA-STD-7009B, §4.2.5](https://standards.nasa.gov/standard/NASA/NASA-STD-7009)
[NASA NESC pressure transient bulletin](https://ntrs.nasa.gov/citations/20220006583)



### Shutdown impulse anchored to valve position

`review_shutdown_impulse(commanded, observed, limits, thrust, close_commands,
window)` replays the issued and observed valve records. Select the issued
closing transition for each valve required for this shutdown. The closure
must be identifiable in the measured state history. A known late
closure can still anchor the impulse calculation. Its timing failure
remains visible in the returned valve review.

The force and command channels must share a nominal time origin and tick
scale. `ShutdownWindow` records an independent bound for force-to-command
alignment over the tail window, including offset, drift and timestamp
quantization. The reducer adds it to the valve-to-command bound. The
earliest and latest possible times of the final closure are the maxima
of the respective per-valve bounds. The window review establishes the
fixed end of the tail and the relationship of the measured thrust to the
selected valves, including whether other operating hardware contributes
to the force channel.

Forward impulse integrates max(thrust, 0) from closure to the fixed end.
Because the integrand is nonnegative, the latest possible closure gives
the lower bound and the earliest gives the upper bound. The integration
clips the force envelopes at the closure bounds, including when they fall
between samples. A later opening command or a possible observed reopening
inside the window prevents a shutdown-only result. The valve record must
cover the full tail with clock uncertainty. The force record must include
a sample at or before the earliest closure, a sample at the reviewed end
and no unapproved gaps.

The returned force review also contains the peak and signed integral over
the sampled window. Those quantities retain their full-window scope.
For closure-relative validation, use a `ShutdownValidationTarget` in the
plan's `transient_target` field. The reviewed digest includes the selected
closing commands, approved latency limits, force channel and tail window,
including the force clock-alignment review. Archive the issued and observed
histories in an `ArchivedValveAcquisition` for the same condition as the
force channel. The v2 acquisition archive covers both valve histories,
their observation coverage and alignment bound, and the force samples.
The v1 fixed-window archive stays readable for earlier records.

Reopen the export with `transient_archive_from_bytes` and pass it to
`assess_validation`. Use `shutdown_measurement` to derive the value and
absolute error for the submitted sample. The reviewer approves that
transcription and its recipe using `validation_sample_sha256`.
Replay resolves the closure from the archived valve state and checks the
submitted value and uncertainty. The report preserves the final closure
interval and the valve timing status.

A model comparison requires a reviewed
`rocket-propulsion-prediction-run-v4` or v5 export. Both carry
`shutdown_scope` on each row. Use `null` for other metrics.
A shutdown row uses unit `N*s` and supplies the following scope.

- `metric` is `post_closure_positive_impulse` and `channel_id` names the
  force channel.
- `closing_commands` lists the selected `head_id` and exact integer
  command `tick` for each valve.
- `window` supplies the four `ShutdownWindow` fields from the approved
  tail and force-clock review.
- `final_closure_ticks` gives the earliest and latest closure bounds in
  force-clock ticks. `tick_period_s` gives the acquisition clock scale.

Replay compares the scope with the reviewed recipe and the bounds
reconstructed from the acquisition. The reviewer must substantiate that
the model prediction bound covers the closure timing interval. An
identifiable late closure can support a numerical model comparison while
the timing assessment still reports failure. Missing prediction scope
leaves model agreement unassessed and preserves the physical acceptance
result from independent measurement replay.
[NASA-STD-7009B, §§4.2.5 and 4.2.7](https://standards.nasa.gov/standard/NASA/NASA-STD-7009)

NASA's AA-1 cold-gas developmental test gives an external example of
command and inferred valve behavior diverging. In Figs. 2–3, the measured
regulator and chamber pressures indicate different closure behavior for
short and long pulses and irregular response at lower supply pressure.
The public report publishes pressure and command plots for those tests.
The poppet-position plot in Fig. 11 is a model result. The public report
does not provide the synchronized force and directly observed valve-state
records needed to replay the closure-relative impulse.
[NASA/TM-2012-217271, §§2 and 4](https://ntrs.nasa.gov/citations/20120003741)



### Physical acceptance criteria for replayed measurements

A `ValidationTarget` can declare `PhysicalAcceptanceCriteria` with an
independent approval record, an inclusive lower bound, an inclusive upper
bound or both. Shutdown targets can also require a passing valve timing
review. Use the measurement units. Assumed limits cannot serve as approved
physical criteria. Include the criteria in the reviewed row digest by
passing `physical_criteria=criteria` to `validation_sample_sha256`.

`physical_acceptance_review` uses the original outward measurement
enclosure from acquisition replay. It records a pass when the full
enclosure satisfies the declared bounds and a failure when the enclosure
is disjoint from the acceptable range. An overlap at a limit leaves the
criterion indeterminate. Required valve timing is evaluated alongside
the measurement band. A confirmed failure remains visible if another
required measurement is missing. Passing requires every declared
criterion to be assessed and satisfied.

When valve timing is required, its review is replayed from the issued and
observed valve histories for the approved test condition. A missing force
channel, an unusable tail window or a mismatched impulse transcription leaves
the impulse unassessed. These conditions do not suppress a verified valve
timing failure. A plan with a timing-only criterion can report passing
valve timing with no force channel. A plan that also requires an impulse
band needs both results. The v2 archive can store the valve histories
without force samples. The target and its archive still require the
reviewed holdout identity and digest.



Physical acceptance is scoped to the named run and its approved limits.
It is evaluated even when the model comparison cannot run. The allowable
model error is a separate tolerance. Matching the measurement can coexist
with a physical limit exceedance. The archive and approvals remain subject
to the physical test review and the relevant environment and scaling
assessment.
[NASA-STD-7009B, §§4.1.1.5 and 4.2.7](https://standards.nasa.gov/standard/NASA/NASA-STD-7009)



## What TRL 5 still needs

- Demonstrate the integrated motor, feed and shutdown chain under the selected
  operating conditions. Compare commanded and observed closure times, pressure
  peak magnitudes and locations, chamber-pressure decay and residual thrust
  with uncertainty.
  The native propagator still requires integration of the additional feed
  states and protection transitions. The SciPy reference now detects the
  reproduced pressure excursion using dense output and pressure-rate
  reversals. Detection remains sensitive to solver tolerance and unresolved
  extrema between scan points. See [coupled trajectory reference](REGULATED_FEED_DYNAMICS.md#coupled-trajectory-reference)
  and [NASA NESC Technical Bulletin 22-03](https://ntrs.nasa.gov/citations/20220006583).
- A defined relevant environment for the user's engine and vehicle
  (pressure, temperature, Mach and Reynolds ranges, propellants).
- Test data from that environment: hot-fire chamber pressure, thrust and
  heat-flux traces; cold-flow feed data; wind-tunnel or flight aerodynamic
  coefficients for the actual vehicle.
- Documented comparison of those data with this software's predictions,
  with uncertainty, using the calibration and model-discrepancy tools already
  in the library (`HOT_FIRE_CALIBRATION.md`, `MODEL_DISCREPANCY_TRANSFER.md`,
  `CALIBRATION_FIELDS.md`).
- Use an independently reviewed `ValidationPlan` with `assess_validation`
  for planned holdout coverage. List the required observation IDs, test
  conditions, units and absolute error allowances before reviewing submitted
  comparisons. Include the conditions used for calibration so that a
  different sensor from the same run cannot enter the holdout. A missing
  sensor, substituted condition, changed unit or inflated allowance prevents
  a planned pass. During independent transcription review, compare each
  sample to the archived measurement and declared model output and record
  `validation_sample_sha256(approved_evidence, reviewed_sample)` in the
  target's `reviewed_sample_sha256`. Retain the approved plan separately
  from the evaluation submission. The digest covers numeric inputs, bounds,
  domain review, citation records, and the declared raw-data and model
  identities. At evaluation the archived file is rehashed. A row with a
  missing or mismatched review digest is unassessed. A matching failure
  remains visible when other targets are missing or unassessed. Omitting
  the plan leaves a comparison of submitted intervals with undeclared
  coverage and review. The reviewer still establishes the accuracy of the
  transcription and the suitability of the limits and uncertainties.
- For an archived model run, record the reviewed export's byte digest in
  `ValidationPlan.reviewed_predictions_sha256` and pass
  `ArchivedModelRun(model_path, predictions_path)` to `assess_validation`.
  The evaluator rehashes the model against `ValidationEvidence.model_sha256`
  and verifies and parses the same bytes of the approved export. The export
  is a UTF-8 JSON document up to 2 MiB, with exactly the keys shown below.
  Each prediction has a unique sample ID. Predictions and bounds use finite
  JSON numbers interpreted as binary64 values. A `null` bound denotes an
  unknown prediction bound. Planned samples must match by test condition,
  unit, prediction and bound. Additional archived predictions are listed
  for review and do not contribute to the comparison. An absent archive
  or mismatched row blocks a planned pass for this option. Model archives
  are read as bytes. The evaluator records file and row agreement, and a
  separate execution record establishes how the model produced the export.

  ```json
  {
    "schema": "rocket-propulsion-prediction-run-v1",
    "test_article_id": "article-A",
    "model_configuration_id": "configuration-7",
    "model_sha256": "<64 lowercase hexadecimal characters>",
    "predictions": [
      {
        "sample_id": "firing-12:injector-peak",
        "condition_record": "firing-12",
        "unit": "Pa",
        "predicted": 1800000,
        "prediction_bound": 2000
      }
    ]
  }
  ```
- Documented scaling requirements from that test article to flight.

No software change can supply these; they come from hardware tests. Until
they exist, no claim of zero modelling error is possible, and the passing
suite is software evidence, not a TRL 5 demonstration. Individual models state
validity ranges; the route fuzzer checks declared scalar admission and smoke
behavior, but a passing suite does not prove that every field, branch, nested
schema or combination of malformed inputs is covered. The current upstream reconciliation and remaining software
review are tracked in [CURRENT_UPSTREAM_REVIEW.md](CURRENT_UPSTREAM_REVIEW.md).

## References

- NASA-NPR-7123.1D, Appendix E, Technology Readiness Levels:
  https://nodis3.gsfc.nasa.gov/displayDir.cfm?Internal_ID=N_PR_7123_001D_&page_name=AppendixE
- NASA-7009B, model credibility and verification reporting:
  https://standards.nasa.gov/sites/default/files/standards/NASA/B/1/NASA-STD-7009B-Final-3-5-2024.pdf
- NASA NTRS 20170005403, ICPTA pressure-fed LOX/LCH4 RCS testing, §III.B and Table 1:
  https://ntrs.nasa.gov/citations/20170005403
