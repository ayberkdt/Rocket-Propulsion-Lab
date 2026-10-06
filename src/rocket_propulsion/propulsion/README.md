# Propulsion

Preliminary rocket-system performance tools, kept separate from gas-dynamic
nozzle equations.

- `performance.py` implements the Tsiolkovsky ideal rocket equation, thrust,
  effective exhaust velocity, specific impulse, and total impulse.
- `propellants.py` provides an immutable, filterable study catalog for common
  liquid bipropellant, solid, and hybrid systems.
- `losses.py` bridges ideal and delivered thrust/Isp through sequential,
  switchable combustion, discharge, boundary-layer, and divergence terms.
- `two_phase.py` provides an explicit condensed-phase coupling band rather than
  pretending to resolve particle size, residence time, slip, or phase change.

Catalog figures are representative study ranges, not design guarantees. Actual
performance depends on chamber pressure, mixture ratio, expansion ratio,
combustion efficiency, and hardware. Toxicity and handling notes are deliberately
high-level; this project does not contain preparation or manufacturing procedures.

