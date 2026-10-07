# Documentation

- `ARCHITECTURE.md` explains dependency direction and extension points.
- `ROADMAP.md` records the staged engineering scope.
- `DEVELOPMENT_PLAN.md` turns the open gaps into a prioritized, test-gated
  implementation plan from version 0.3 to 1.0.
- `REFERENCES.md` records primary equation and screening-data sources.
- `GAP_ANALYSIS.md` distinguishes covered physics from the next fidelity limits.
- `USER_GUIDE.md` describes the end-to-end browser workflow and fidelity labels.
- `API.md` documents the JSON envelope and Python/PowerShell examples.
- `TOLERANCES.md` centralizes numerical acceptance thresholds.
- `VALIDATION_REPORT.md` records closed evidence and the remaining external CEA gate.
- `RELEASE_NOTES.md` summarizes packaged engineering-preview changes.
- `HOT_FIRE_CALIBRATION.md` documents scalar and paired test-data calibration.
- `MODEL_DISCREPANCY_TRANSFER.md` separates model-form and
  qualification-to-flight epistemic effects from hot-fire variability.
- `PROPAGATOR_PROPULSION_BRIDGE.md` defines the event-safe force, total/tank
  mass derivative, estimator-driver, and analytic-Jacobian boundary for Sidera.
- `OPERATING_CONDITION_SURFACES.md` documents bounded multidimensional
  supply-pressure/ambient/throttle/power performance maps and their gradients.
- `REGULATED_FEED_DYNAMICS.md` defines the pressurant bottle, dynamic
  regulator, thermal ullage, line-loss ODE and implicit feed/engine closure.
- `DYNAMIC_FEED_LINE.md` defines the propagated liquid-flow/manifold-pressure
  states, conservation law, cutoff roots and composed Jacobians used to add a
  calibrated first feed-line mode to an orbit propagator.

Equations live next to their implementation in Python docstrings so the code and
documentation cannot drift independently.

