"""Compose discrete engine outcomes with conditional performance variation.

The probabilities and uncertainty intervals are illustrative software fixtures,
not generic engine reliability or qualification values.

References
----------
NASA-STD-7009B: https://standards.nasa.gov/standard/NASA/NASA-STD-7009
NASA-SP-2011-3421: https://ntrs.nasa.gov/citations/20120001369
NASA-20050207429: https://ntrs.nasa.gov/citations/20050207429
NASA-20160007073: https://ntrs.nasa.gov/citations/20160007073
"""

from __future__ import annotations

from rocket_propulsion.integrations import export_sidera_ensemble_handoff
from rocket_propulsion.propulsion.burns import (
    ConstantPerformanceProvider,
    DistributionKind,
    EngineClusterMember,
    EngineTransientCommand,
    PropulsionScenario,
    PropulsionScenarioSet,
    ScenarioConditionalUncertainty,
    ScenarioEvent,
    ScenarioEventType,
    ThrottleSchedule,
    ThrottleSegment,
    UncertainParameter,
    UncertaintyClass,
    build_hierarchical_propulsion_ensemble,
    build_propulsion_scenario_ensemble,
    hierarchical_ensemble_hash,
    hierarchical_ensemble_to_json,
)


def scale(name: str, lower: float, upper: float) -> UncertainParameter:
    """Build one sourced bounded conditional scale.

    References
    ----------
    NASA-STD-7009B: https://standards.nasa.gov/standard/NASA/NASA-STD-7009
    """

    return UncertainParameter(
        name=name,
        nominal=1.0,
        lower=lower,
        upper=upper,
        distribution=DistributionKind.TRIANGULAR,
        uncertainty_class=UncertaintyClass.EPISTEMIC,
        source="illustrative qualification-data fixture",
        rationale="demonstrate conditional bounded performance uncertainty",
    )


def main() -> None:
    """Build and summarize a reproducible hierarchical ensemble.

    References
    ----------
    NASA-SP-2011-3421: https://ntrs.nasa.gov/citations/20120001369
    NASA-20160007073: https://ntrs.nasa.gov/citations/20160007073
    """

    point = ConstantPerformanceProvider.monopropellant(
        delivered_thrust_n=100.0,
        system_specific_impulse_s=220.0,
        chamber_pressure_pa=1.0e6,
        source="illustrative engine point",
    ).operating_point()
    members = (
        EngineClusterMember("engine-a", point),
        EngineClusterMember("engine-b", point),
    )
    schedule = ThrottleSchedule((ThrottleSegment(10.0, 1.0, 1.0, "steady"),))
    commands = tuple(EngineTransientCommand(item.engine_id, schedule) for item in members)
    scenarios = PropulsionScenarioSet(
        scenarios=(
            PropulsionScenario(
                "nominal",
                0.90,
                (),
                "illustrative event tree",
                "both engines complete the command",
            ),
            PropulsionScenario(
                "engine-a-out",
                0.08,
                (
                    ScenarioEvent(
                        ScenarioEventType.ENGINE_DISABLED,
                        ("engine-a",),
                        None,
                        "illustrative reliability fixture",
                        "engine-a does not ignite",
                    ),
                ),
                "illustrative event tree",
                "single-engine outcome",
            ),
            PropulsionScenario(
                "common-cause-no-ignition",
                0.02,
                (
                    ScenarioEvent(
                        ScenarioEventType.ENGINE_DISABLED,
                        ("engine-a", "engine-b"),
                        None,
                        "illustrative common-cause fixture",
                        "shared initiation system unavailable",
                    ),
                ),
                "illustrative event tree",
                "dependent loss of both engines",
            ),
        ),
        objective="illustrate joint propulsion uncertainty",
    )
    discrete = build_propulsion_scenario_ensemble(members, commands, scenarios)
    performance_scales = (
        scale("thrust_scale", 0.98, 1.02),
        scale("specific_impulse_scale", 0.99, 1.01),
        scale("time_scale", 0.99, 1.01),
    )
    definitions = (
        ScenarioConditionalUncertainty(
            "nominal",
            performance_scales,
            32,
            "illustrative nominal qualification fixture",
            "nominal engine-to-engine performance scatter",
        ),
        ScenarioConditionalUncertainty(
            "engine-a-out",
            performance_scales,
            32,
            "illustrative engine-out qualification fixture",
            "surviving engine performance scatter",
        ),
        ScenarioConditionalUncertainty(
            "common-cause-no-ignition",
            (),
            1,
            "illustrative deterministic failure fixture",
            "performance variation cannot restore disabled engines",
        ),
    )
    result = build_hierarchical_propulsion_ensemble(
        discrete,
        definitions,
        seed=20261006,
        replicate_tolerance=0.10,
    )
    impulse = next(
        item for item in result.statistics if item.metric == "delivered_impulse_n_s"
    )
    print(f"realizations: {result.realization_count}")
    print(f"joint probability: {result.joint_probability_sum:.12f}")
    print(f"expected delivered impulse: {impulse.expected_value:.6f} N*s")
    print(f"P05/P50/P95: {impulse.p05:.6f} / {impulse.p50:.6f} / {impulse.p95:.6f}")
    print(f"zero-thrust probability: {result.zero_thrust_probability:.6f}")
    decomposition = next(
        item
        for item in result.variance_decomposition
        if item.metric == "delivered_impulse_n_s"
    )
    print(
        "impulse variance fractions: "
        f"within={decomposition.within_fraction:.6f}, "
        f"between={decomposition.between_fraction:.6f}"
    )
    print(f"ensemble hash: {hierarchical_ensemble_hash(result)}")
    print(f"serialized bytes: {len(hierarchical_ensemble_to_json(result).encode('utf-8'))}")
    handoff = export_sidera_ensemble_handoff(
        result,
        t_start_s=120.0,
        direction=(1.0, 0.0, 0.0),
        frame="inertial",
    )
    print(f"Sidera hand-off jobs: {len(handoff.jobs)}")
    print(f"Sidera hand-off hash: {handoff.handoff_hash}")
    for item in result.sampling_evidence:
        print(
            f"{item.scenario_id}: N={item.sample_count}, "
            f"replicate-converged={item.converged}"
        )


if __name__ == "__main__":
    main()

