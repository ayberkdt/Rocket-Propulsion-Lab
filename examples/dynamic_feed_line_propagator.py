"""Evaluate a pressure-dependent engine with propagated feed-line states."""

from rocket_propulsion.propulsion.burns import (
    DynamicFeedLine,
    DynamicFeedLineConfiguration,
    DynamicFeedLineState,
    IdealPressurantGas,
    PerformanceAxis,
    PerformanceScalePoint,
    PropagatorPropulsionBridge,
    RectilinearPerformanceSurface,
    RegulatedFeedConfiguration,
    RegulatedFeedSystem,
    TabulatedBurnArtifact,
    TabulatedBurnSegment,
    evaluate_dynamic_feed_propulsion,
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
    source_id="example-dynamic-feed-grid",
    source_sha256="6" * 64,
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
    source_id="example-dynamic-feed-profile",
    source_sha256="7" * 64,
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
        # DynamicFeedLine owns the line resistance in this composition.
        feed_line_resistance_pa_s2_kg2=0.0,
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
feed_state = feed.initial_state(
    propellant_mass_kg=10.0,
    tank_pressure_pa=2.0e6,
    tank_gas_temperature_k=300.0,
    bottle_pressure_pa=20.0e6,
    bottle_gas_temperature_k=300.0,
    regulator_opening=0.5,
)

# At 1.8 MPa supply and 50 kPa ambient the map requests 0.217 kg/s.
nominal_flow = 0.217
line = DynamicFeedLine(
    DynamicFeedLineConfiguration(
        inertance_pa_s2_kg=1.0e6,
        compliance_kg_pa=1.0e-8,
        resistance_pa_s2_kg2=(2.0e6 - 1.8e6) / nominal_flow**2,
        minimum_manifold_pressure_pa=1.0e6,
        maximum_manifold_pressure_pa=2.3e6,
        maximum_absolute_flow_kg_s=0.5,
    )
)
line_state = DynamicFeedLineState(
    mass_flow_kg_s=nominal_flow,
    manifold_pressure_pa=1.8e6,
)

result = evaluate_dynamic_feed_propulsion(
    bridge,
    feed,
    feed_state,
    line,
    line_state,
    relative_time_s=2.0,
    vehicle_mass_kg=100.0,
    direction=(1.0, 0.0, 0.0),
    pressure_condition_name="supply_pressure_pa",
    other_operating_conditions={"ambient_pressure_pa": 50_000.0},
)

print("tank pressure [Pa]:", result.regulated_feed.tank_pressure_pa)
print("manifold pressure [Pa]:", line_state.manifold_pressure_pa)
print("force [N]:", result.propulsion.force_vector_n)
print("engine demand [kg/s]:", -result.propulsion.mass_derivative_kg_s)
print("line flow derivative [kg/s^2]:", result.feed_line.derivative.mass_flow_kg_s2)
print("manifold derivative [Pa/s]:", result.feed_line.derivative.manifold_pressure_pa_s)
print("line storage rate [kg/s]:", result.feed_line.stored_liquid_mass_rate_kg_s)
print("mass closure residual [kg/s]:", result.feed_line.external_mass_closure_error_kg_s)
print("line state Jacobian:", result.partials.line_state_jacobian)
print("line model hash:", line.model_hash)
