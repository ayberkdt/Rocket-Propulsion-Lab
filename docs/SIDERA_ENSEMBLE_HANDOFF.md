# Weighted propulsion ensemble hand-off to Sidera

## Status

`rocket_propulsion_sidera_ensemble_handoff_v1` is the implemented producer-side
contract for a future Sidera weighted transient-burn consumer. It is generated
and verified entirely in Rocket Propulsion Lab. No Sidera source file is
modified, and current Sidera is correctly reported as missing the public
capabilities required for native consumption.

## Ownership boundary

Rocket Propulsion Lab supplies:

- one complete `rocket_propulsion_tabulated_burn_v1` artifact per realization;
- scenario, conditional, and joint probabilities;
- relative thrust and total tank-flow histories;
- propulsion source hashes and the parent ensemble hash.

The caller supplies Sidera execution inputs:

- solver-relative start time;
- unit direction and inertial/RIC/VNB frame.

Sidera remains responsible for trajectory state, epoch interpretation, frame
evaluation, attitude/guidance, gravity, events, numerical propagation, and
orbital results.

## Force and mass channels

The hand-off fixes two independent physical channels:

```text
force channel     = axial_thrust_n
spacecraft dm/dt  = -total_tank_mass_flow_kg_s
```

Delivered thrust is the engine/cluster production before hardware cant loss.
Axial thrust is the net force available along the commanded propulsion axis.
Total tank flow includes all declared chamber, gas-generator, coolant-dump,
vent, and other vehicle-depleting streams. Sidera must not reconstruct mass
flow from axial thrust and chamber/system Isp because that would either erase
cant loss or under-count propellant drain.

The existing exact-constant adapter now observes the same identity. It sends
constant axial thrust together with rocket-equivalent Isp, which lets Sidera's
current `F/(g0 Isp)` law reproduce the true total tank drain.

## Probability rules

Every job carries:

```text
scenario_probability
conditional_probability
joint_probability = scenario_probability * conditional_probability
```

The schema verifies conditional closure inside every scenario and joint closure
over the full hand-off. Zero-thrust jobs remain present. A consumer applies the
joint probability once, after propagating the corresponding history; it must
not renormalize successful runs after dropping failures.

## Integrity chain

The hand-off records:

- the parent `rocket_propulsion_hierarchical_ensemble_v1` hash;
- every nested burn artifact and artifact hash;
- absolute start/end times, frame, normalized direction, and weights;
- fixed physical-channel semantics and required capabilities;
- a SHA-256 over the complete hand-off document.

Reopening validates nested artifacts before validating the outer document.
Changing a force/mass channel label is rejected as a physical-semantic change,
not treated as optional metadata.

## Required future Sidera capabilities

Schema v1 requires:

1. a public tabulated finite-burn type;
2. mass flow independent of thrust/Isp;
3. propagation splitting at every artifact segment boundary;
4. a weighted maneuver-ensemble execution/result contract.

`missing_sidera_ensemble_capabilities()` reports the exact missing set.
`require_sidera_ensemble_capabilities()` fails closed until all public
capabilities exist. Constant `FiniteBurn` support alone is not considered
sufficient.

## Consumer algorithm

For each job, a future Sidera consumer should:

1. verify the hand-off and nested artifact hashes;
2. combine the job start time with all relative segment boundaries;
3. split propagation at every boundary and honor left/right discontinuities;
4. interpolate axial thrust and total tank flow independently;
5. evaluate direction in the declared frame from the current trajectory state;
6. store the source ensemble, hand-off, and artifact hashes with the result;
7. combine trajectory outputs using the declared joint weights;
8. retain invalid/terminated outcomes with explicit probability accounting.

## Primary references

- CCSDS 502.0-B-3, *Orbit Data Messages*:
  <https://ccsds.org/Pubs/502x0b3e1.pdf>
- NASA-STD-7009B, *Standard for Models and Simulations*:
  <https://standards.nasa.gov/standard/NASA/NASA-STD-7009>
- NASA/SP-2011-3421, *Probabilistic Risk Assessment Procedures Guide for NASA
  Managers and Practitioners*:
  <https://ntrs.nasa.gov/citations/20120001369>

