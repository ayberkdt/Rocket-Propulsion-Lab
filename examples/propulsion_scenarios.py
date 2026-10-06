"""Build nominal, engine-out, and common-cause propulsion histories.

References
----------
NASA-SP-2011-3421: https://ntrs.nasa.gov/citations/20120001369
NASA-20050207429: https://ntrs.nasa.gov/citations/20050207429
NASA-20160007073: https://ntrs.nasa.gov/citations/20160007073
"""

from __future__ import annotations

from rocket_propulsion.propulsion.burns import (
    ConstantPerformanceProvider,
    EngineClusterMember,
    EngineTransientCommand,
    PropulsionScenario,
    PropulsionScenarioSet,
    ScenarioEvent,
    ScenarioEventType,
    ThrottleSchedule,
    ThrottleSegment,
    build_propulsion_scenario_ensemble,
)


def main() -> None:
    """Print weighted impulse consequences for an illustrative event tree.

    The probabilities are software examples only, not generic engine failure
    rates or flight-reliability claims.

    References
    ----------
    NASA-SP-2011-3421: https://ntrs.nasa.gov/citations/20120001369
    NASA-20050207429: https://ntrs.nasa.gov/citations/20050207429
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
    commands = (
        EngineTransientCommand("engine-a", schedule),
        EngineTransientCommand("engine-b", schedule),
    )
    scenario_source = "illustrative event-tree fixture"
    event_source = "illustrative reliability fixture"
    scenario_set = PropulsionScenarioSet(
        scenarios=(
            PropulsionScenario(
                "nominal",
                0.90,
                (),
                scenario_source,
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
                        event_source,
                        "engine-a fails to ignite",
                    ),
                ),
                scenario_source,
                "single-engine contingency",
            ),
            PropulsionScenario(
                "common-cause-no-ignition",
                0.02,
                (
                    ScenarioEvent(
                        ScenarioEventType.ENGINE_DISABLED,
                        ("engine-a", "engine-b"),
                        None,
                        event_source,
                        "shared initiation system unavailable",
                    ),
                ),
                scenario_source,
                "dependent loss of both engines",
            ),
        ),
        objective="illustrate propulsion consequence weighting",
    )
    result = build_propulsion_scenario_ensemble(members, commands, scenario_set)
    impulse = next(
        statistic
        for statistic in result.statistics
        if statistic.metric == "delivered_impulse_n_s"
    )
    print(f"expected impulse: {impulse.expected_value:.9g} N*s")
    print(f"P05/P50/P95: {impulse.p05:.9g} / {impulse.p50:.9g} / {impulse.p95:.9g}")
    print(f"zero-thrust probability: {result.zero_thrust_probability:.9g}")


if __name__ == "__main__":
    main()
