# Sidera propulsion-provider contract

`SideraPropulsionProviderHandoff` is the primary interface between Rocket
Propulsion Lab and Sidera. It is intentionally independent of Sidera's Python
installation: a producer can build, hash, serialize, reopen, and evaluate the
complete propulsion model before the native Sidera consumer is present.

## Ownership

Rocket Propulsion Lab owns engine performance, named propellant drains, feed
and line dynamics, inventory limits, propulsion events, and analytic partials.
Sidera owns the absolute epoch, state vector, frames, attitude/guidance, event
location, integrator restarts, and translational/rotational propagation.

The hand-off converts Sidera epoch to bridge time as
`relative_time_s = epoch_s - ignition_epoch_s`. Its manifest records:

- the bridge and optional mounted-cluster hashes;
- vehicle-mass and additional-state layout;
- required operating-condition names and estimable drivers;
- burn, inventory, operating-domain, and interpolation-grid event IDs;
- force frame, mass/time semantics, and Sidera capability requirements.

## Chemical and bipropellant states

A stream-resolved surface independently interpolates force, oxidizer flow, and
fuel flow. `BipropellantDynamicFeedCoupling` combines that map with two
regulated tanks and two inertance/compliance feed lines. It exposes engine and
line mixture ratios, minimum/maximum O/F roots, first-component depletion,
component/line/vehicle mass-closure residuals, a four-state line Jacobian, and
acceleration derivatives through both manifold pressures.

Sidera should append every `additional_state_order` entry to its continuous
state and integrate the returned named derivative in exactly that order. It
must register every event in the manifest and restart on non-stopping grid
roots so variational equations do not step across derivative kinks.

## 6-DOF thrust vectoring

With a `MountedEngineCluster`, `evaluate_body_frame` composes ordered two-axis
body-fixed gimbal rotations, engine forces, and moments about Sidera's supplied
instantaneous centre of mass. The sum of enabled engine-command thrust
magnitudes must close to the propulsion bridge's scalar thrust; a mismatch is
refused. Sidera transforms the returned body force using its attitude state and
integrates the returned body moment. Rocket Propulsion Lab does not propagate
attitude or choose the gimbal commands.

## Sidera changes required

Current Sidera's constant `FiniteBurn` remains supported for eligible constant
results, but the full provider requires these public host capabilities:

1. a propulsion-provider callback returning force/acceleration and independent
   vehicle/named-state mass derivatives at every RHS call;
2. augmented-state registration and stable named-state ordering;
3. provider event-detector registration with left/right evaluation and restart;
4. analytic partial channels for STM/OD use;
5. optional body-frame force/moment ingestion for 6-DOF propagation.

The capability probe reports these fields separately. Until Sidera implements
them, the hand-off fails closed as an integration artifact instead of silently
reducing a transient, multi-stream engine to an equivalent constant Isp.
`missing_sidera_propulsion_capabilities(...)` reports the exact gaps for a
specific hand-off; `require_sidera_propulsion_capabilities(...)` is the native
consumer gate and raises a typed `FeatureUnavailableError` until every required
channel is present.

## References

- NASA-20040000363, transient liquid-propulsion system modeling:
  https://ntrs.nasa.gov/citations/20040000363
- NASA-20050206420, computer-controlled mixture-ratio study:
  https://ntrs.nasa.gov/citations/20050206420
- NASA-SP-8059, spacecraft attitude control during thrusting:
  https://ntrs.nasa.gov/citations/19710016722
- NASA/TM-20240012525, three-dimensional TVC gimbal geometry:
  https://ntrs.nasa.gov/citations/20240012525
- Orekit 13.1.5 `PropulsionModel` public contract:
  https://www.orekit.org/static/apidocs/org/orekit/forces/maneuvers/propulsion/PropulsionModel.html
- NASA-STD-7009B, model credibility and verification:
  https://standards.nasa.gov/standard/NASA/NASA-STD-7009
