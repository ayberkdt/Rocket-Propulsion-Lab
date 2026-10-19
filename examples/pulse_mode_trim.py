"""Small trim maneuver with a pulse-mode thruster: MIB, quantization, artifact.

The impulse-bit table is synthetic and illustrative, not a real thruster.
"""

from rocket_propulsion.propulsion.burns import PulseModeThruster, free_space_burn_metrics
from rocket_propulsion.propulsion.performance import STANDARD_GRAVITY_M_S2

flow = 1.0 / (STANDARD_GRAVITY_M_S2 * 220.0)
on_times = (0.010, 0.020, 0.050, 0.100, 0.200)
thruster = PulseModeThruster(
    on_times_s=on_times,
    impulse_bits_n_s=tuple(1.0 * (t - 0.006) + 0.002 for t in on_times),
    propellant_per_pulse_kg=tuple(flow * (t + 0.003) for t in on_times),
    steady_thrust_n=1.0,
    steady_mass_flow_kg_s=flow,
    minimum_off_time_s=0.02,
    source_id="synthetic-1N",
    source_sha256="c" * 64,
    maximum_on_time_s=2.0,
)
print(f"MIB {thruster.minimum_impulse_bit_n_s * 1e3:.1f} mN*s, "
      f"Isp at MIB {thruster.pulse(0.010).specific_impulse_s:.0f} s")

mass = 12.0
for request in (0.0002, 0.0004, 0.05, 1.0):
    plan = thruster.plan_delta_v(request, mass, propellant_budget_kg=2.0)
    width = "-" if plan.pulse is None else f"{plan.pulse.on_time_s * 1e3:.2f} ms"
    print(f"dv {request:8.4f} m/s -> {plan.pulse_count:3d} pulses of {width:>9}, "
          f"delivered {plan.delivered:.6f} m/s, residual {plan.residual:+.2e}, "
          f"quantized={plan.quantized}")

plan = thruster.plan_delta_v(0.05, mass, propellant_budget_kg=2.0)
artifact = thruster.pulse_train_artifact(plan, period_s=plan.pulse.on_time_s + 0.05)
metrics = free_space_burn_metrics(artifact, mass)
print(f"pulse-train artifact: {len(artifact.segments)} segments, "
      f"dv {metrics.delta_v_m_s:.9f} m/s, final mass {metrics.final_mass_kg:.6f} kg")
