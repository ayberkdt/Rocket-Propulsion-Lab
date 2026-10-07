"""Couple supply/ambient conditions to a transient propagator force call."""

from rocket_propulsion.propulsion.burns import (
    PerformanceAxis,
    PerformanceScalePoint,
    PropagatorPropulsionBridge,
    RectilinearPerformanceSurface,
    TabulatedBurnArtifact,
    TabulatedBurnSegment,
)

surface = RectilinearPerformanceSurface(
    axes=(
        PerformanceAxis("supply_pressure_pa", "Pa", (1.0e6, 2.0e6)),
        PerformanceAxis("ambient_pressure_pa", "Pa", (0.0, 100_000.0)),
    ),
    points=(
        PerformanceScalePoint((1.0e6, 0.0), 0.80, 0.90),
        PerformanceScalePoint((1.0e6, 100_000.0), 0.70, 0.95),
        PerformanceScalePoint((2.0e6, 0.0), 1.00, 1.10),
        PerformanceScalePoint((2.0e6, 100_000.0), 0.90, 1.15),
    ),
    source_id="example-qualified-operating-grid",
    source_sha256="2" * 64,
)

artifact = TabulatedBurnArtifact(
    segments=(
        TabulatedBurnSegment(
            0.0,
            10.0,
            500.0,
            500.0,
            490.0,
            490.0,
            0.2,
            0.2,
            phase="steady",
        ),
    ),
    source_id="example-time-profile",
    source_sha256="3" * 64,
)

bridge = PropagatorPropulsionBridge(
    artifact,
    performance_surface=surface,
)
evaluation = bridge.evaluate(
    relative_time_s=2.0,
    vehicle_mass_kg=100.0,
    direction=(1.0, 0.0, 0.0),
    operating_conditions={
        "supply_pressure_pa": 1.5e6,
        "ambient_pressure_pa": 50_000.0,
    },
)

print("force scale:", evaluation.performance_force_scale)
print("flow scale:", evaluation.performance_mass_flow_scale)
print("force [N]:", evaluation.force_vector_n)
print("mass rate [kg/s]:", evaluation.mass_derivative_kg_s)
print("d(acceleration)/d(condition):")
for name, derivative in evaluation.partials.acceleration_wrt_conditions:
    print(" ", name, derivative)
print("validity events:", len(bridge.operating_condition_surfaces()))
print("surface hash:", surface.surface_hash)
