"""Evaluate two regulated tanks and lines against a stream-resolved engine."""

from dataclasses import replace

from rocket_propulsion.propulsion.burns import (
    AdditionalStateConstraint,
    BipropellantDynamicFeedConfiguration,
    BipropellantDynamicFeedCoupling,
    BipropellantDynamicFeedState,
    DynamicFeedLine,
    DynamicFeedLineConfiguration,
    DynamicFeedLineState,
    IdealPressurantGas,
    PerformanceAxis,
    PropagatorPropulsionBridge,
    RectilinearStreamPerformanceSurface,
    RegulatedFeedConfiguration,
    RegulatedFeedSystem,
    StreamPerformanceScalePoint,
    StreamStateBinding,
    TabulatedBurnArtifact,
    TabulatedBurnSegment,
)

axes = (
    PerformanceAxis("oxidizer_pressure_pa", "Pa", (1.0e6, 2.0e6)),
    PerformanceAxis("fuel_pressure_pa", "Pa", (1.0e6, 2.0e6)),
)
surface = RectilinearStreamPerformanceSurface(
    axes=axes,
    points=tuple(
        StreamPerformanceScalePoint(
            (oxidizer_pressure, fuel_pressure),
            0.8 + 1.0e-7 * (oxidizer_pressure - 1.0e6) + 5.0e-8 * (fuel_pressure - 1.0e6),
            (
                ("oxidizer", 0.9 + 2.0e-7 * (oxidizer_pressure - 1.0e6)),
                ("fuel", 1.1 + 1.0e-7 * (fuel_pressure - 1.0e6)),
            ),
        )
        for oxidizer_pressure in axes[0].values
        for fuel_pressure in axes[1].values
    ),
    source_id="example-bipropellant-grid",
    source_sha256="8" * 64,
)
artifact = TabulatedBurnArtifact(
    segments=(
        TabulatedBurnSegment(
            0.0,
            10.0,
            150.0,
            150.0,
            150.0,
            150.0,
            3.0,
            3.0,
            (("oxidizer", 2.25), ("fuel", 0.75)),
            (("oxidizer", 2.25), ("fuel", 0.75)),
            "steady",
        ),
    ),
    source_id="example-bipropellant-profile",
    source_sha256="9" * 64,
)
bridge = PropagatorPropulsionBridge(
    artifact,
    stream_performance_surface=surface,
    stream_bindings=(
        StreamStateBinding("oxidizer", "oxidizer_mass_kg"),
        StreamStateBinding("fuel", "fuel_mass_kg"),
    ),
    state_constraints=(
        AdditionalStateConstraint("oxidizer_mass_kg", 1.0),
        AdditionalStateConstraint("fuel_mass_kg", 1.0),
    ),
)

base_configuration = RegulatedFeedConfiguration(
    gas=IdealPressurantGas(2_077.1, 1.667, "helium"),
    tank_internal_volume_m3=0.02,
    propellant_density_kg_m3=1_000.0,
    bottle_volume_m3=0.01,
    regulator_set_pressure_pa=2.1e6,
    regulator_full_open_error_pa=0.2e6,
    regulator_maximum_area_m2=1.0e-7,
    regulator_discharge_coefficient=0.8,
    valve_time_constant_s=0.05,
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
oxidizer_feed = RegulatedFeedSystem(replace(base_configuration, propellant_density_kg_m3=1_200.0))
fuel_feed = RegulatedFeedSystem(replace(base_configuration, propellant_density_kg_m3=800.0))
oxidizer_feed_state = oxidizer_feed.initial_state(
    propellant_mass_kg=10.0,
    tank_pressure_pa=2.0e6,
    tank_gas_temperature_k=300.0,
    bottle_pressure_pa=20.0e6,
    bottle_gas_temperature_k=300.0,
    regulator_opening=0.5,
)
fuel_feed_state = fuel_feed.initial_state(
    propellant_mass_kg=6.0,
    tank_pressure_pa=1.9e6,
    tank_gas_temperature_k=300.0,
    bottle_pressure_pa=20.0e6,
    bottle_gas_temperature_k=300.0,
    regulator_opening=0.5,
)

oxidizer_flow = 2.25 * 1.06
fuel_flow = 0.75 * 1.16
oxidizer_line = DynamicFeedLine(
    DynamicFeedLineConfiguration(
        1.0e6,
        1.0e-8,
        (2.0e6 - 1.8e6) / oxidizer_flow**2,
        1.0e6,
        2.2e6,
        5.0,
    )
)
fuel_line = DynamicFeedLine(
    DynamicFeedLineConfiguration(
        1.2e6,
        2.0e-8,
        (1.9e6 - 1.6e6) / fuel_flow**2,
        1.0e6,
        2.2e6,
        3.0,
    )
)
state = BipropellantDynamicFeedState(
    oxidizer_feed_state,
    DynamicFeedLineState(oxidizer_flow, 1.8e6),
    fuel_feed_state,
    DynamicFeedLineState(fuel_flow, 1.6e6),
)
coupling = BipropellantDynamicFeedCoupling(
    BipropellantDynamicFeedConfiguration(
        "oxidizer_mass_kg",
        "fuel_mass_kg",
        "oxidizer_pressure_pa",
        "fuel_pressure_pa",
        2.5,
        3.0,
    ),
    bridge,
    oxidizer_feed,
    oxidizer_line,
    fuel_feed,
    fuel_line,
)

result = coupling.evaluate(
    state,
    relative_time_s=0.5,
    vehicle_mass_kg=100.0,
    direction=(1.0, 0.0, 0.0),
)
print("force [N]:", result.propulsion.force_vector_n)
print("oxidizer engine flow [kg/s]:", result.oxidizer_engine_flow_kg_s)
print("fuel engine flow [kg/s]:", result.fuel_engine_flow_kg_s)
print("engine O/F:", result.engine_mixture_ratio)
print("line O/F:", result.line_mixture_ratio)
print("limiting component:", result.depletion.limiting_component)
print("vehicle mass closure [kg/s]:", result.vehicle_mass_closure_error_kg_s)
print("line-state Jacobian:", result.partials.line_state_jacobian)
print("coupling model hash:", coupling.model_hash)
