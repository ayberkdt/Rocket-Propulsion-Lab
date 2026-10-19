"""Inputs an optimizer needs: a(u), da/du, smooth Isp(u), impulsive equivalent.

Numbers are illustrative.
"""

from rocket_propulsion.propulsion.burns import (
    MonotoneCubicCurve,
    TabulatedBurnArtifact,
    TabulatedBurnSegment,
    ThrottleableEngineModel,
    free_space_burn_metrics,
)

engine = ThrottleableEngineModel(
    rated_thrust_n=30_000.0,
    rated_mass_flow_kg_s=9.5,
    minimum_throttle=0.4,
    maximum_throttle=1.05,
    stream_fractions=(("m_ox", 0.78), ("m_fu", 0.22)),
    # Thrust fraction below throttle at deep throttle: Isp falls.
    thrust_curve=MonotoneCubicCurve((0.4, 0.6, 0.8, 1.0, 1.05), (0.37, 0.58, 0.79, 1.0, 1.05)),
)
for throttle in (0.4, 0.55, 0.7, 0.85, 1.0):
    value = engine.evaluate(
        throttle=throttle,
        vehicle_mass_kg=1500.0,
        direction=(1.0, 0.0, 0.0),
        additional_states={"m_ox": 400.0, "m_fu": 120.0},
    )
    print(
        f"u={throttle:4.2f}  a={value.acceleration_m_s2[0]:7.3f} m/s^2  "
        f"da/du={value.partials.acceleration_wrt_throttle[0]:7.3f}  "
        f"Isp={value.vacuum_specific_impulse_s:6.1f} s"
    )

artifact = TabulatedBurnArtifact(
    segments=(
        TabulatedBurnSegment(0.0, 2.0, 0.0, 30_000.0, 0.0, 30_000.0, 0.0, 9.5, phase="start"),
        TabulatedBurnSegment(2.0, 120.0, 30_000.0, 30_000.0, 30_000.0, 30_000.0, 9.5, 9.5,
                             phase="steady"),
    ),
    source_id="example",
    source_sha256="0" * 64,
)
metrics = free_space_burn_metrics(artifact, initial_mass_kg=6000.0)
print(f"delta-v {metrics.delta_v_m_s:.4f} m/s, impulse at t* = "
      f"{metrics.impulsive_equivalent_time_s:.3f} s, error estimate "
      f"{metrics.delta_v_error_estimate_m_s:.1e} m/s")
