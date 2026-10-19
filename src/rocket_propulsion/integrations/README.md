# Integrations

Adapters in this package define the propulsion/astrodynamics boundary without
moving orbit, frame, attitude, or event-location logic into this repository.

`sidera_propulsion.py` is the primary integration path. Its
`rocket_propulsion_sidera_provider_handoff_v1` document embeds the complete
propagator bridge, declares absolute/relative time mapping, additional-state
order, operating-condition inputs, estimable drivers, and every stopping or
restart event. An optional mounted engine cluster closes scalar engine thrust
to a body-frame force and moment with analytic gimbal derivatives. The hand-off
is executable and integrity-checkable without Sidera, so both repositories can
test the same contract before native integration.

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

`lunaris.py` is retained as a tested legacy compatibility adapter. New finite-
burn integration and capability growth target Sidera.
