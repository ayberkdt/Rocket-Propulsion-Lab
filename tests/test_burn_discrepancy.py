"""Model discrepancy, qualification transfer, and layered ensemble tests."""

from __future__ import annotations

import copy
import unittest

from rocket_propulsion.core.errors import DomainError, InputError
from rocket_propulsion.propulsion.burns import (
    DiscrepancyKind,
    DiscrepancyRatioObservation,
    DistributionKind,
    HotFireVectorMeasurement,
    PropulsionScenarioEnsemble,
    PropulsionScenarioRealization,
    ScaleDiscrepancyCalibration,
    ScenarioConditionalUncertainty,
    TabulatedBurnArtifact,
    TabulatedBurnSegment,
    UncertainParameter,
    UncertaintyClass,
    build_hierarchical_propulsion_ensemble,
    build_layered_calibration_plan,
    calibrate_hot_fire_multivariate,
    calibrate_scale_discrepancy,
)


def hot_fire_calibration():
    engine_effects = (-0.025, -0.015, -0.005, 0.005, 0.015, 0.025)
    run_effects = ((-0.003, -0.001), (0.0, 0.002), (0.003, -0.001))
    covariance = ((1.0e-6, 0.1e-6), (0.1e-6, 0.64e-6))
    measurements = tuple(
        HotFireVectorMeasurement(
            engine_id=f"E{engine_index + 1}",
            run_id=f"R{run_index + 1}",
            condition_id="qualification-rated",
            parameter_names=("thrust_scale", "specific_impulse_scale"),
            values=(
                1.01 + engine_effect + run_effect[0],
                0.995 + 0.5 * engine_effect + run_effect[1],
            ),
            measurement_covariance=covariance,
            source="qualification hot-fire evidence HF-Q1",
        )
        for engine_index, engine_effect in enumerate(engine_effects)
        for run_index, run_effect in enumerate(run_effects)
    )
    return calibrate_hot_fire_multivariate(
        measurements,
        source="qualification hot-fire campaign HF-Q1",
        rationale="engine and run variability at the qualification operating point",
    )


def discrepancy(
    target: str,
    component: str,
    kind: DiscrepancyKind,
    ratios: tuple[float, ...],
    evidence_id: str,
) -> ScaleDiscrepancyCalibration:
    observations = tuple(
        DiscrepancyRatioObservation(
            case_id=f"case-{index + 1}",
            applicability_id="flight-block-A-rated",
            ratio=ratio,
            standard_uncertainty=0.002,
            source=f"independent validation case {index + 1}",
        )
        for index, ratio in enumerate(ratios)
    )
    return calibrate_scale_discrepancy(
        target,
        component,
        kind,
        observations,
        evidence_id=evidence_id,
        source=f"reviewed {kind.value} campaign",
        rationale="predict performance inside the declared flight application domain",
    )


def base_artifact() -> TabulatedBurnArtifact:
    return TabulatedBurnArtifact(
        segments=(
            TabulatedBurnSegment(
                start_time_s=0.0,
                end_time_s=10.0,
                start_delivered_thrust_n=100.0,
                end_delivered_thrust_n=100.0,
                start_axial_thrust_n=90.0,
                end_axial_thrust_n=90.0,
                start_total_mass_flow_kg_s=0.05,
                end_total_mass_flow_kg_s=0.05,
                start_stream_mass_flows_kg_s=(("main", 0.05),),
                end_stream_mass_flows_kg_s=(("main", 0.05),),
            ),
        ),
        source_id="layered-base",
        source_sha256="a" * 64,
    )


class DiscrepancyCalibrationTests(unittest.TestCase):
    def setUp(self) -> None:
        self.model_form = discrepancy(
            "thrust_scale",
            "validation-residual",
            DiscrepancyKind.MODEL_FORM,
            (0.98, 0.99, 1.00, 1.01, 1.02, 1.00),
            "MF-V1",
        )
        self.transfer = discrepancy(
            "specific_impulse_scale",
            "qualification-flight-transfer",
            DiscrepancyKind.QUALIFICATION_TO_FLIGHT,
            (0.990, 0.995, 1.000, 1.005, 1.010, 0.998),
            "QF-T1",
        )

    def test_predictive_epistemic_variance_and_persistence(self) -> None:
        self.assertAlmostEqual(self.model_form.consensus_scale, 1.0, delta=0.002)
        self.assertGreater(
            self.model_form.predictive_epistemic_variance,
            self.model_form.between_case_variance,
        )
        parameter = self.model_form.residual_uncertain_parameter()
        self.assertEqual(parameter.uncertainty_class, UncertaintyClass.EPISTEMIC)
        self.assertTrue(parameter.name.startswith("thrust_scale::epistemic:model-form:"))
        self.assertIn(self.model_form.calibration_hash, parameter.source)

        reopened = ScaleDiscrepancyCalibration.from_json(self.model_form.to_json())
        self.assertEqual(reopened, self.model_form)
        changed = copy.deepcopy(self.model_form.to_dict())
        changed["between_case_variance"] += 0.01
        with self.assertRaisesRegex(InputError, "semantics"):
            ScaleDiscrepancyCalibration.from_dict(changed)

    def test_layered_plan_preserves_classes_and_physical_scale_products(self) -> None:
        hot_fire = hot_fire_calibration()
        plan = build_layered_calibration_plan(
            hot_fire,
            (self.model_form, self.transfer),
            scenario_id="nominal",
            applicability_id="flight-block-A-rated",
            hot_fire_evidence_id="HF-Q1",
            sample_count=32,
            source="flight block A uncertainty qualification",
            rationale="keep repeatability, model form, and transfer evidence separate",
        )
        definition = plan.conditional_definition
        classes = tuple(item.uncertainty_class for item in definition.parameters)
        self.assertEqual(classes.count(UncertaintyClass.ALEATORY), 2)
        self.assertEqual(classes.count(UncertaintyClass.EPISTEMIC), 2)
        self.assertEqual(len({item.name for item in definition.parameters}), 4)
        self.assertTrue(all("::" in item.name for item in definition.parameters))

        base = base_artifact()
        corrected = plan.corrected_artifact(base)
        deterministic = dict(plan.deterministic_scales)
        expected_impulse_ratio = deterministic["thrust_scale"]
        expected_mass_ratio = (
            deterministic["thrust_scale"] / deterministic["specific_impulse_scale"]
        )
        self.assertAlmostEqual(
            corrected.delivered_total_impulse_n_s / base.delivered_total_impulse_n_s,
            expected_impulse_ratio,
        )
        self.assertAlmostEqual(
            corrected.consumed_propellant_kg / base.consumed_propellant_kg,
            expected_mass_ratio,
        )

        scenarios = PropulsionScenarioEnsemble(
            objective="layered scale integration fixture",
            scenarios=(),
            realizations=(PropulsionScenarioRealization("nominal", 1.0, corrected, None),),
            statistics=(),
            zero_thrust_probability=0.0,
            reference_ids=(),
            warnings=(),
        )
        ensemble = build_hierarchical_propulsion_ensemble(
            scenarios,
            (definition,),
            seed=17,
        )
        self.assertEqual(len(ensemble.realizations), 32)
        for realization in ensemble.realizations:
            values = dict(realization.inputs)
            thrust_factor = 1.0
            isp_factor = 1.0
            for name, value in values.items():
                if name.startswith("thrust_scale::"):
                    thrust_factor *= value
                elif name.startswith("specific_impulse_scale::"):
                    isp_factor *= value
            self.assertAlmostEqual(
                realization.artifact.delivered_total_impulse_n_s,
                corrected.delivered_total_impulse_n_s * thrust_factor,
            )
            self.assertAlmostEqual(
                realization.artifact.consumed_propellant_kg,
                corrected.consumed_propellant_kg * thrust_factor / isp_factor,
            )

    def test_applicability_and_evidence_reuse_fail_closed(self) -> None:
        hot_fire = hot_fire_calibration()
        with self.assertRaisesRegex(DomainError, "reused"):
            build_layered_calibration_plan(
                hot_fire,
                (self.model_form,),
                scenario_id="nominal",
                applicability_id="flight-block-A-rated",
                hot_fire_evidence_id="MF-V1",
                sample_count=8,
                source="invalid double-count fixture",
                rationale="prove evidence reuse refusal",
            )

        shared_evidence = discrepancy(
            "specific_impulse_scale",
            "same-campaign-isp",
            DiscrepancyKind.MODEL_FORM,
            (0.99, 1.00, 1.01, 0.995, 1.005, 1.00),
            "MF-V1",
        )
        with self.assertRaisesRegex(DomainError, "requires an explicit correlation"):
            build_layered_calibration_plan(
                hot_fire,
                (self.model_form, shared_evidence),
                scenario_id="nominal",
                applicability_id="flight-block-A-rated",
                hot_fire_evidence_id="HF-Q1",
                sample_count=8,
                source="invalid shared-evidence fixture",
                rationale="prove shared evidence correlation refusal",
            )

        mixed = tuple(
            DiscrepancyRatioObservation(
                case_id=f"case-{index}",
                applicability_id=("domain-A" if index < 2 else "domain-B"),
                ratio=1.0,
                standard_uncertainty=0.01,
                source="mixed-domain fixture",
            )
            for index in range(3)
        )
        with self.assertRaisesRegex(DomainError, "mix applicability"):
            calibrate_scale_discrepancy(
                "thrust_scale",
                "invalid",
                DiscrepancyKind.MODEL_FORM,
                mixed,
                evidence_id="invalid",
                source="invalid fixture",
                rationale="prove domain pooling refusal",
            )

        mismatched_class = UncertainParameter(
            name="thrust_scale::aleatory:mislabelled",
            nominal=1.0,
            lower=0.99,
            upper=1.01,
            distribution=DistributionKind.UNIFORM,
            uncertainty_class=UncertaintyClass.EPISTEMIC,
            source="invalid qualified-name fixture",
            rationale="prove class/name mismatch refusal",
        )
        with self.assertRaisesRegex(DomainError, "encode their uncertainty class"):
            ScenarioConditionalUncertainty(
                scenario_id="invalid",
                parameters=(mismatched_class,),
                sample_count=8,
                source="invalid qualified-name fixture",
                rationale="prove qualified scale integrity",
            )


if __name__ == "__main__":
    unittest.main()
