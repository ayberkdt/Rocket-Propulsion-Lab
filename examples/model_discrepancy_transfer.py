"""Compose hot-fire variability with independent model/flight discrepancy."""

from rocket_propulsion.propulsion.burns import (
    DiscrepancyKind,
    DiscrepancyRatioObservation,
    HotFireVectorMeasurement,
    TabulatedBurnArtifact,
    TabulatedBurnSegment,
    build_layered_calibration_plan,
    calibrate_hot_fire_multivariate,
    calibrate_scale_discrepancy,
)

PARAMETERS = ("thrust_scale", "specific_impulse_scale")
ENGINE_EFFECTS = (-0.025, -0.015, -0.005, 0.005, 0.015, 0.025)
RUN_EFFECTS = ((-0.003, -0.001), (0.0, 0.002), (0.003, -0.001))
MEASUREMENT_COVARIANCE = ((1.0e-6, 0.1e-6), (0.1e-6, 0.64e-6))

hot_fire_measurements = tuple(
    HotFireVectorMeasurement(
        engine_id=f"E{engine_index + 1}",
        run_id=f"R{run_index + 1}",
        condition_id="qualification-rated",
        parameter_names=PARAMETERS,
        values=(
            1.01 + engine_effect + run_effect[0],
            0.995 + 0.5 * engine_effect + run_effect[1],
        ),
        measurement_covariance=MEASUREMENT_COVARIANCE,
        source="synthetic qualification hot-fire campaign HF-Q1",
    )
    for engine_index, engine_effect in enumerate(ENGINE_EFFECTS)
    for run_index, run_effect in enumerate(RUN_EFFECTS)
)
hot_fire = calibrate_hot_fire_multivariate(
    hot_fire_measurements,
    source="synthetic qualification hot-fire campaign HF-Q1",
    rationale="rated-condition engine and run variability",
)

model_form = calibrate_scale_discrepancy(
    "thrust_scale",
    "validation-residual",
    DiscrepancyKind.MODEL_FORM,
    tuple(
        DiscrepancyRatioObservation(
            case_id=f"V{index + 1}",
            applicability_id="flight-block-A-rated",
            ratio=ratio,
            standard_uncertainty=0.002,
            source="synthetic independent model-validation campaign MF-V1",
        )
        for index, ratio in enumerate((0.98, 0.99, 1.00, 1.01, 1.02, 1.00))
    ),
    evidence_id="MF-V1",
    source="synthetic independent model-validation campaign MF-V1",
    rationale="rated flight-block A prediction discrepancy",
)

plan = build_layered_calibration_plan(
    hot_fire,
    (model_form,),
    scenario_id="nominal",
    applicability_id="flight-block-A-rated",
    hot_fire_evidence_id="HF-Q1",
    sample_count=256,
    source="flight block A uncertainty qualification example",
    rationale="separate repeatability from model-form uncertainty",
)

base = TabulatedBurnArtifact(
    segments=(
        TabulatedBurnSegment(
            start_time_s=0.0,
            end_time_s=10.0,
            start_delivered_thrust_n=100.0,
            end_delivered_thrust_n=100.0,
            start_axial_thrust_n=95.0,
            end_axial_thrust_n=95.0,
            start_total_mass_flow_kg_s=0.05,
            end_total_mass_flow_kg_s=0.05,
            start_stream_mass_flows_kg_s=(("main", 0.05),),
            end_stream_mass_flows_kg_s=(("main", 0.05),),
        ),
    ),
    source_id="synthetic-base-artifact",
    source_sha256="a" * 64,
)
corrected = plan.corrected_artifact(base)

print("Deterministic consensus scales")
for name, value in plan.deterministic_scales:
    print(f"  {name}: {value:.9f}")
print("Residual components")
for parameter in plan.conditional_definition.parameters:
    print(f"  {parameter.name}: {parameter.uncertainty_class.value}")
print(f"Base impulse: {base.delivered_total_impulse_n_s:.6f} N*s")
print(f"Corrected impulse: {corrected.delivered_total_impulse_n_s:.6f} N*s")
print(f"Plan hash: {plan.plan_hash}")
print(f"Discrepancy hash: {model_form.calibration_hash}")
