"""Build a paired thrust/Isp/duration residual model from hot-fire data."""

from rocket_propulsion.propulsion.burns import (
    HotFireVectorMeasurement,
    calibrate_hot_fire_multivariate,
)

PARAMETERS = ("thrust_scale", "specific_impulse_scale", "time_scale")
MEASUREMENT_COVARIANCE = (
    (1.00e-6, 0.20e-6, -0.10e-6),
    (0.20e-6, 0.64e-6, 0.08e-6),
    (-0.10e-6, 0.08e-6, 0.49e-6),
)
ENGINE_EFFECTS = (-0.030, -0.018, -0.006, 0.006, 0.018, 0.030)
SECONDARY_EFFECTS = (0.002, -0.001, 0.0015, -0.002, 0.0005, -0.001)
RUN_PATTERNS = (
    ((-0.004, 0.002, -0.001), (0.001, -0.003, 0.003), (0.003, 0.001, -0.002)),
    ((-0.003, 0.003, -0.002), (0.002, -0.002, 0.001), (0.001, -0.001, 0.001)),
)

measurements = tuple(
    HotFireVectorMeasurement(
        engine_id=f"E{engine_index + 1}",
        run_id=f"R{run_index + 1}",
        condition_id="rated-100pct",
        parameter_names=PARAMETERS,
        values=(
            1.010 + engine_effect + run_effect[0],
            0.997 + 0.55 * engine_effect + SECONDARY_EFFECTS[engine_index] + run_effect[1],
            1.004 - 0.30 * engine_effect
            + 0.5 * SECONDARY_EFFECTS[engine_index]
            + run_effect[2],
        ),
        measurement_covariance=MEASUREMENT_COVARIANCE,
        source="synthetic paired qualification campaign example",
    )
    for engine_index, engine_effect in enumerate(ENGINE_EFFECTS)
    for run_index, run_effect in enumerate(RUN_PATTERNS[engine_index % 2])
)

calibration = calibrate_hot_fire_multivariate(
    measurements,
    source="synthetic qualification campaign Q2",
    rationale="demonstrate paired thrust, Isp, and duration calibration",
)
parameters, correlation = calibration.residual_model()

print("Consensus corrections")
for marginal in calibration.marginals:
    print(f"  {marginal.parameter_name}: {marginal.consensus_scale:.9f}")
print("Population correlation")
for row in correlation.matrix:
    print("  " + " ".join(f"{value:+.6f}" for value in row))
print(f"Regularized: {calibration.population_correlation.regularized}")
print(f"Identity shrinkage: {calibration.population_correlation.identity_shrinkage:.9g}")
print(f"Residual parameters: {len(parameters)}")
print(f"Calibration hash: {calibration.calibration_hash}")
