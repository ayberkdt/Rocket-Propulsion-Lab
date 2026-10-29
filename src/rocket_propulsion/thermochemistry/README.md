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

All dimensional results use SI units. These fits use their defining
R=8.314510 J/(mol K), from NASA/TP-2002-211556 Appendix A, rather than the
2019 SI constant used by other models. Polynomial enthalpy includes the
species' formation reference. The first two gas intervals of all ten bundled
species have been checked against the pinned official NASA/CEA database;
this subset remains limited to 200--6000 K.

Catalogue calls return a fresh dictionary over owned immutable records.
Both packaged catalogues are read under the strict JSON rules. Species
properties refuse extrapolation; raw polynomial primitives remain available
outside their interval for reference checks. Every exported property is
formed from exact rational terms with enclosed logarithms. See
`docs/NASA9_CATALOGUE_REVIEW.md` for numerical scope and retained failures.
Frozen mixtures and the experimental local equilibrium provider are separate
layers. Their existence does not establish independent equilibrium validation,
hardware operating limits or TRL5 readiness.

Frozen-mixture forward contracts are now reviewed separately in
`docs/MIXTURE_FORWARD_REVIEW.md`: aliases combined into one species identity,
retained raw composition and basis, explicit exact normalization of accepted fraction sums,
bounded logarithm exports and owned replay-checked property records. Inverse
root accuracy, equilibrium composition and physical applicability remain open.
