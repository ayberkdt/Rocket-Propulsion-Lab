# Tests

Reference-value tests cover each public engineering relation and its API adapter,
including duct flow, nozzles, Tsiolkovsky performance, thermodynamic connections,
oblique shocks, nozzle geometry, polytropic paths, sampled plot curves, and the
propellant catalog.
Values use common perfect-air textbook cases, especially `gamma = 1.4`, `M = 2`.

Run all tests with `scripts/test.ps1`. New modules should include:

1. a nominal reference case;
2. physical-domain rejection cases;
3. forward/inverse round trips where applicable;
4. route-level payload validation.

