# Thermochemistry

Temperature-dependent ideal-gas properties based on the seven-term heat-capacity
representation documented in NASA/TP-2002-211556.

- `nasa_polynomials.py` implements the published `cp/R`, `h/RT`, and `s/R`
  equations.
- `species.py` loads source-labelled species data and prevents silent
  extrapolation.
- `mixtures.py` converts mass/mole fractions and evaluates ideal-mixture
  properties, including entropy of mixing and pressure dependence.
- `data/nasa_glenn_2002.json` is a machine-readable subset of the published NASA
  Glenn database for rocket-combustion species.

All dimensional results use SI units. Polynomial enthalpy includes the species'
heat of formation. This package calculates frozen-composition properties; it
does not yet solve chemical equilibrium.
