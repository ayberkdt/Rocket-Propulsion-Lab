"""Select an electric-thruster mode from available power and evaluate it.

The throttle table is synthetic and only illustrates the contract; replace it
with qualified data for a real thruster.
"""

from rocket_propulsion.propulsion.burns import (
    DiscreteThrottlePath,
    ElectricThrottleLevel,
    ElectricThrottleTable,
    PowerLimitedPropulsionModel,
)

table = ElectricThrottleTable(
    levels=(
        ElectricThrottleLevel("low", 1000.0, 0.040, 1.5e-6),
        ElectricThrottleLevel("mid", 2000.0, 0.080, 2.6e-6),
        ElectricThrottleLevel("high", 3000.0, 0.110, 3.2e-6),
    ),
    source_id="example-synthetic-throttle-table",
    source_sha256="0" * 64,
)
model = PowerLimitedPropulsionModel(
    DiscreteThrottlePath.from_table(table, objective="thrust", up_switch_hysteresis_w=50.0),
    propellant_reserve_kg=2.0,
    power_margin_w=100.0,
)

available_power_w = 2400.0
level = model.select_level(available_power_w)
evaluation = model.evaluate(
    level_index=level,
    available_power_w=available_power_w,
    vehicle_mass_kg=500.0,
    direction=(0.0, 1.0, 0.0),
    additional_states={"xenon_mass_kg": 20.0},
)

print("level:", evaluation.level_id)
print("thrust [N]:", evaluation.thrust_n)
print("Isp [s]:", round(evaluation.specific_impulse_s or 0.0, 1))
print("consumed power [W]:", evaluation.consumed_power_w)
print("xenon rate [kg/s]:", evaluation.mass_derivative_kg_s)
for surface in model.switching_surfaces(level):
    print("switch:", surface.event_id, surface.crossing, surface.threshold_w, "->", surface.to_level)
print("model hash:", model.model_hash)
