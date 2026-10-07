"""Evaluate a transient propulsion artifact at a propagator stage."""

from rocket_propulsion.propulsion.burns import (
    AdditionalStateConstraint,
    PropagatorPropulsionBridge,
    StreamStateBinding,
    TabulatedBurnArtifact,
    TabulatedBurnSegment,
)

artifact = TabulatedBurnArtifact(
    segments=(
        TabulatedBurnSegment(
            0.0,
            0.5,
            0.0,
            500.0,
            0.0,
            490.0,
            0.0,
            0.2,
            (("oxidizer", 0.0), ("fuel", 0.0)),
            (("oxidizer", 0.14), ("fuel", 0.06)),
            "ignition-rise",
        ),
        TabulatedBurnSegment(
            0.5,
            5.0,
            500.0,
            500.0,
            490.0,
            490.0,
            0.2,
            0.2,
            (("oxidizer", 0.14), ("fuel", 0.06)),
            (("oxidizer", 0.14), ("fuel", 0.06)),
            "steady",
        ),
    ),
    source_id="example-qualified-profile",
    source_sha256="1" * 64,
    reference_ids=("NASA-20040000363",),
)

bridge = PropagatorPropulsionBridge(
    artifact,
    stream_bindings=(
        StreamStateBinding("oxidizer", "oxidizer_mass_kg"),
        StreamStateBinding("fuel", "fuel_mass_kg"),
    ),
    state_constraints=(
        AdditionalStateConstraint("oxidizer_mass_kg", 0.5),
        AdditionalStateConstraint("fuel_mass_kg", 0.2),
    ),
    protected_dry_mass_kg=80.0,
)

evaluation = bridge.evaluate(
    relative_time_s=0.25,
    vehicle_mass_kg=100.0,
    direction=(1.0, 1.0, 0.0),
    additional_states={"oxidizer_mass_kg": 4.0, "fuel_mass_kg": 2.0},
)

print("active:", evaluation.active)
print("phase:", evaluation.phase)
print("force [N]:", evaluation.force_vector_n)
print("acceleration [m/s^2]:", evaluation.acceleration_m_s2)
print("vehicle mass rate [kg/s]:", evaluation.mass_derivative_kg_s)
print("tank rates [kg/s]:", dict(evaluation.additional_state_derivatives))
print("event count:", len(bridge.event_surfaces()))
print("model hash:", bridge.model_hash)
