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

## What TRL 5 still needs

- A defined relevant environment for the user's engine and vehicle
  (pressure, temperature, Mach and Reynolds ranges, propellants).
- Test data from that environment: hot-fire chamber pressure, thrust and
  heat-flux traces; cold-flow feed data; wind-tunnel or flight aerodynamic
  coefficients for the actual vehicle.
- Documented comparison of those data with this software's predictions,
  with uncertainty, using the calibration and model-discrepancy tools already
  in the library (`HOT_FIRE_CALIBRATION.md`, `MODEL_DISCREPANCY_TRANSFER.md`,
  `CALIBRATION_FIELDS.md`).
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
