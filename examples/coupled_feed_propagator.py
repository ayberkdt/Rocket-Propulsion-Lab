"""Resolve a regulated feed state and pressure-dependent finite burn together."""

from rocket_propulsion.propulsion.burns import (
    IdealPressurantGas,
    PerformanceAxis,
    PerformanceScalePoint,
    PropagatorPropulsionBridge,
    RectilinearPerformanceSurface,
    RegulatedFeedConfiguration,
    RegulatedFeedSystem,
    TabulatedBurnArtifact,
    TabulatedBurnSegment,
    evaluate_coupled_feed_propulsion,
)

surface = RectilinearPerformanceSurface(
    axes=(
        PerformanceAxis("supply_pressure_pa", "Pa", (1.0e6, 2.0e6)),
        PerformanceAxis("ambient_pressure_pa", "Pa", (0.0, 100_000.0)),
    ),
    points=(
        PerformanceScalePoint((1.0e6, 0.0), 0.8, 0.9),
        PerformanceScalePoint((1.0e6, 100_000.0), 0.7, 0.95),
        PerformanceScalePoint((2.0e6, 0.0), 1.0, 1.1),
        PerformanceScalePoint((2.0e6, 100_000.0), 0.9, 1.15),
    ),
    source_id="example-feed-performance-grid",
    source_sha256="4" * 64,
)
artifact = TabulatedBurnArtifact(
    segments=(
        TabulatedBurnSegment(
            0.0,
            20.0,
            500.0,
            500.0,
            490.0,
            490.0,
            0.2,
            0.2,
            phase="steady",
        ),
    ),
    source_id="example-feed-coupled-profile",
    source_sha256="5" * 64,
)
bridge = PropagatorPropulsionBridge(artifact, performance_surface=surface)

feed = RegulatedFeedSystem(
    RegulatedFeedConfiguration(
        gas=IdealPressurantGas(2_077.1, 1.667, "helium"),
        tank_internal_volume_m3=0.02,
        propellant_density_kg_m3=1_000.0,
        bottle_volume_m3=0.01,
        regulator_set_pressure_pa=2.1e6,
        regulator_full_open_error_pa=0.2e6,
        regulator_maximum_area_m2=1.0e-7,
        regulator_discharge_coefficient=0.8,
        valve_time_constant_s=0.05,
        feed_line_resistance_pa_s2_kg2=1_000.0,
        tank_heat_transfer_w_k=0.0,
        tank_wall_temperature_k=300.0,
        bottle_heat_transfer_w_k=0.0,
        bottle_wall_temperature_k=300.0,
        minimum_propellant_mass_kg=1.0,
        minimum_injector_pressure_pa=1.0e6,
        minimum_bottle_pressure_margin_pa=0.1e6,
        maximum_tank_pressure_pa=2.3e6,
    )
)
state = feed.initial_state(
    propellant_mass_kg=10.0,
    tank_pressure_pa=2.0e6,
    tank_gas_temperature_k=300.0,
    bottle_pressure_pa=20.0e6,
    bottle_gas_temperature_k=300.0,
    regulator_opening=0.5,
)

result = evaluate_coupled_feed_propulsion(
    bridge,
    feed,
    state,
    relative_time_s=2.0,
    vehicle_mass_kg=100.0,
    direction=(1.0, 0.0, 0.0),
    pressure_condition_name="supply_pressure_pa",
    other_operating_conditions={"ambient_pressure_pa": 50_000.0},
    initial_mass_flow_kg_s=0.2,
    relaxation=0.8,
)

print("injector pressure [Pa]:", result.feed.injector_pressure_pa)
print("force [N]:", result.propulsion.force_vector_n)
print("propellant flow [kg/s]:", -result.propulsion.mass_derivative_kg_s)
print("regulator flow [kg/s]:", result.feed.regulator_mass_flow_kg_s)
print("coupling iterations:", result.evidence.iterations)
print("flow residual [kg/s]:", result.evidence.mass_flow_residual_kg_s)
print("implicit Jacobian denominator:", result.partials.algebraic_denominator)
print("feed model hash:", feed.model_hash)
