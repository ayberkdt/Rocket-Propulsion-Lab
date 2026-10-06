# Compressible flow

Pure functions for one-dimensional, calorically perfect-gas flow.

- `isentropic.py` computes Mach angle, Prandtl-Meyer angle, stagnation ratios,
  critical ratios, area ratio, and inverse conversions.
- `normal_shock.py` computes downstream Mach number and static/stagnation jumps.
- `oblique_shock.py` solves weak and strong attached theta-beta-M branches and
  reports the detachment limit.
- `duct_flow.py` evaluates Fanno (friction) and Rayleigh (heat-transfer) lines.
- `nozzle.py` couples isentropic area relations to choked mass flow and thrust.
- `sweeps.py` supplies bounded sampled curves for plots and parameter studies.

All public result objects are immutable dataclasses. Real-gas and viscous-loss
models can be added with the same pattern.

