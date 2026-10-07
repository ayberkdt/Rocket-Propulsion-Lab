"""Reduce a two-burn timeline to exact-mass impulses for Lunaris.

Works without Lunaris (prints the impulses); with Lunaris installed it also
builds the native ManeuverPlan. Numbers are illustrative.
"""

from rocket_propulsion.core.errors import FeatureUnavailableError
from rocket_propulsion.integrations import (
    lunaris_impulsive_equivalents,
    to_lunaris_maneuver_plan,
)
from rocket_propulsion.propulsion.burns import (
    PropagatorPropulsionBridge,
    PropulsionTimeline,
    ScheduledBurn,
    TabulatedBurnArtifact,
    TabulatedBurnSegment,
)


def burn(thrust_n: float, flow_kg_s: float, duration_s: float, name: str):
    artifact = TabulatedBurnArtifact(
        segments=(TabulatedBurnSegment(0.0, duration_s, thrust_n, thrust_n, thrust_n, thrust_n,
                                       flow_kg_s, flow_kg_s, phase="steady"),),
        source_id=name,
        source_sha256="0" * 64,
    )
    return PropagatorPropulsionBridge(artifact, require_complete_stream_binding=False)


timeline = PropulsionTimeline(
    burns=(
        ScheduledBurn("lower", burn(400.0, 0.14, 60.0, "deorbit-1"), 100.0),
        ScheduledBurn("trim", burn(400.0, 0.14, 15.0, "trim-1"), 3600.0),
    )
)
impulses = lunaris_impulsive_equivalents(
    timeline,
    initial_mass_kg=1000.0,
    directions={"lower": ((0.0, -1.0, 0.0), "ric"), "trim": ((0.0, 1.0, 0.0), "ric")},
)
for item in impulses:
    print(
        f"{item.burn_id:6s} t={item.t_burn_s:9.3f} s  dv={item.delta_v_m_s:8.4f} m/s  "
        f"isp={item.isp_s:7.2f} s  mass {item.initial_mass_kg:.3f} -> {item.final_mass_kg:.3f} kg  "
        f"(compresses a {item.finite_burn_duration_s:.0f} s burn)"
    )
try:
    plan = to_lunaris_maneuver_plan(impulses)
    print("Lunaris plan sha256:", plan.plan_sha256)
except FeatureUnavailableError as error:
    print("Lunaris not installed:", error)
