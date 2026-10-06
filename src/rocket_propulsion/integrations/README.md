# Integrations

Adapters in this package transform immutable Rocket Propulsion Lab results;
they do not participate in the propulsion calculation.

The Sidera adapter exports an exact constant finite burn using delivered
thrust, system-effective Isp, and unity target throttle. Unity is essential:
the source operating point has already resolved its commanded throttle, so
passing that throttle to Sidera again would scale thrust and mass flow twice.
For canted propulsion, the adapter sends net axial thrust and rocket-equivalent
Isp. This preserves both the reduced net force and the original total tank
drain through Sidera's current ideal `F/(g0 Isp)` mass-flow law.

Sidera is optional and is never discovered through a local checkout path.
Artifact export works without Sidera. Native object construction uses only the
public `sidera.core.propagation` facade and fails closed if its constructor
capability is incompatible.

`sidera_ensemble.py` implements the producer side of
`rocket_propulsion_sidera_ensemble_handoff_v1`. It binds each weighted
tabulated artifact to start time, normalized direction, and frame, fixes axial
thrust and total tank flow as independent channels, retains zero-thrust jobs,
and reports the exact public capabilities current Sidera still lacks. It does
not construct native objects or alter Sidera.
