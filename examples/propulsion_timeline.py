"""Two upper-stage burns with an adapter jettison, integrated with RK4.

The burns are constant-thrust profiles with c = F / mdot = 3000 m/s, so the
integrated velocity can be checked against the rocket equation applied burn
by burn with the jettison in between. Numbers are illustrative.
"""

from math import log

from rocket_propulsion.propulsion.burns import (
    MassJettison,
    PropagatorPropulsionBridge,
    PropulsionTimeline,
    ScheduledBurn,
    StreamStateBinding,
    TabulatedBurnArtifact,
    TabulatedBurnSegment,
    TimelineEventKind,
)

EXHAUST = 3000.0


def burn(thrust_n: float, duration_s: float, name: str) -> PropagatorPropulsionBridge:
    flow = thrust_n / EXHAUST
    streams = (("ox", 0.78 * flow), ("fu", 0.22 * flow))
    artifact = TabulatedBurnArtifact(
        segments=(
            TabulatedBurnSegment(0.0, duration_s, thrust_n, thrust_n, thrust_n, thrust_n,
                                 flow, flow, streams, streams, "steady"),
        ),
        source_id=name,
        source_sha256="0" * 64,
        reference_ids=("NASA-20040000363",),
    )
    return PropagatorPropulsionBridge(
        artifact,
        stream_bindings=(StreamStateBinding("ox", "m_ox"), StreamStateBinding("fu", "m_fu")),
    )


timeline = PropulsionTimeline(
    burns=(
        ScheduledBurn("perigee", burn(30_000.0, 120.0, "perigee-burn"), 10.0),
        ScheduledBurn("apogee", burn(30_000.0, 40.0, "apogee-burn"), 400.0),
    ),
    jettisons=(MassJettison("adapter", 200.0, 150.0, "payload adapter"),),
)

time_s, mass, ox, fu, speed = 0.0, 6000.0, 2400.0, 700.0, 0.0
timed = [event for event in timeline.events() if event.time_s is not None]
end_s, step_s = 500.0, 1.0


def rates(t: float, m: float, o: float, f: float, side: str):
    value = timeline.evaluate(
        time_s=t, vehicle_mass_kg=m,
        directions={"perigee": (1.0, 0.0, 0.0), "apogee": (1.0, 0.0, 0.0)},
        additional_states={"m_ox": o, "m_fu": f}, side=side, apply_inventory_limits=False,
    )
    drains = dict(value.additional_state_derivatives)
    return value.acceleration_m_s2[0], value.mass_derivative_kg_s, drains["m_ox"], drains["m_fu"]


while time_s < end_s:
    target = min(end_s, time_s + step_s)
    landing = [e.time_s for e in timed if time_s + 1e-12 < e.time_s <= target]
    target = min(landing) if landing else target
    h = target - time_s
    y = (speed, mass, ox, fu)
    k1 = rates(time_s, *y[1:], "right")
    y2 = tuple(a + h / 2 * b for a, b in zip(y, k1, strict=True))
    k2 = rates(time_s + h / 2, *y2[1:], "right")
    y3 = tuple(a + h / 2 * b for a, b in zip(y, k2, strict=True))
    k3 = rates(time_s + h / 2, *y3[1:], "right")
    y4 = tuple(a + h * b for a, b in zip(y, k3, strict=True))
    k4 = rates(time_s + h, *y4[1:], "left")
    speed, mass, ox, fu = (
        a + h / 6 * (b + 2 * c + 2 * d + e) for a, b, c, d, e in zip(y, k1, k2, k3, k4, strict=True)
    )
    time_s = target
    for event in timed:
        if abs(event.time_s - time_s) < 1e-12 and event.kind is TimelineEventKind.JETTISON:
            mass = timeline.apply_jettison(event.owner_id, mass)

m1 = 6000.0 - 10.0 * 120.0
m2 = m1 - 150.0
m3 = m2 - 10.0 * 40.0
expected = EXHAUST * (log(6000.0 / m1) + log(m2 / m3))
print(f"integrated delta-v {speed:.6f} m/s, rocket equation {expected:.6f} m/s")
print(f"final mass {mass:.6f} kg (expected {m3:.6f} kg)")
print("timeline hash:", timeline.timeline_hash)
