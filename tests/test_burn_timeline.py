"""Multi-burn timelines: scheduling, jettison, summation, events, persistence."""

import copy
import unittest
from math import isclose, log

from rocket_propulsion.core.errors import DomainError, InputError
from rocket_propulsion.propulsion.burns import (
    AdditionalStateConstraint,
    MassJettison,
    PropagatorPropulsionBridge,
    PropulsionTimeline,
    ScheduledBurn,
    StreamStateBinding,
    TabulatedBurnArtifact,
    TabulatedBurnSegment,
    TimelineEventKind,
)

C = 50.0  # axial thrust / total flow in every profile below


def steady_bridge(thrust: float, duration: float, *, exit_area: float | None = None):
    flow = thrust / C
    streams = (("ox", 0.75 * flow), ("fu", 0.25 * flow))
    artifact = TabulatedBurnArtifact(
        segments=(
            TabulatedBurnSegment(0.0, duration, thrust, thrust, thrust, thrust, flow, flow,
                                 streams, streams, "steady"),
        ),
        source_id=f"steady-{thrust}-{duration}",
        source_sha256="4" * 64,
        reference_ids=("NASA-20040000363",),
    )
    return PropagatorPropulsionBridge(
        artifact,
        stream_bindings=(StreamStateBinding("ox", "m_ox"), StreamStateBinding("fu", "m_fu")),
        state_constraints=(AdditionalStateConstraint("m_ox", 1.0),),
        nozzle_exit_area_m2=exit_area,
    )


def two_burn_timeline() -> PropulsionTimeline:
    return PropulsionTimeline(
        burns=(
            ScheduledBurn("upper-2", steady_bridge(100.0, 4.0), 20.0),
            ScheduledBurn("upper-1", steady_bridge(150.0, 6.0), 2.0),
        ),
        jettisons=(MassJettison("interstage", 12.0, 30.0, "interstage adapter"),),
        timeline_id="two-burn-test",
    )


def integrate(timeline: PropulsionTimeline, end_s: float, step_s: float):
    """Event-landing RK4 over the timeline, applying jettisons."""

    time_s, mass, ox, fu = 0.0, 300.0, 60.0, 20.0
    velocity = 0.0
    events = timeline.events()
    timed = [event for event in events if event.time_s is not None]
    log_entries = []

    def rate(t, m, o, f, side):
        value = timeline.evaluate(
            time_s=t, vehicle_mass_kg=m,
            directions={burn.burn_id: (1.0, 0.0, 0.0) for burn in timeline.burns},
            additional_states={"m_ox": o, "m_fu": f}, side=side,
            apply_inventory_limits=False,
        )
        drains = dict(value.additional_state_derivatives)
        return value.acceleration_m_s2[0], value.mass_derivative_kg_s, drains["m_ox"], drains["m_fu"]

    while time_s < end_s - 1e-12:
        target = min(end_s, time_s + step_s)
        landing = [event for event in timed if time_s + 1e-12 < event.time_s <= target + 1e-12]
        if landing:
            target = min(event.time_s for event in landing)
        h = target - time_s
        y = (velocity, mass, ox, fu)

        def shifted(base, k, f):
            return tuple(a + f * b for a, b in zip(base, k, strict=True))

        def deriv(t, state, side):
            return rate(t, state[1], state[2], state[3], side)

        k1 = deriv(time_s, y, "right")
        k2 = deriv(time_s + h / 2, shifted(y, k1, h / 2), "right")
        k3 = deriv(time_s + h / 2, shifted(y, k2, h / 2), "right")
        k4 = deriv(time_s + h, shifted(y, k3, h), "left")
        velocity, mass, ox, fu = (
            a + h / 6 * (b + 2 * c + 2 * d + e)
            for a, b, c, d, e in zip(y, k1, k2, k3, k4, strict=True)
        )
        time_s = target
        for event in timed:
            if abs(event.time_s - time_s) <= 1e-12:
                log_entries.append(event.event_id)
                if event.kind is TimelineEventKind.JETTISON:
                    mass = timeline.apply_jettison(event.owner_id, mass)
    return velocity, mass, ox, fu, log_entries


class TimelineTests(unittest.TestCase):
    def test_two_burns_with_jettison_follow_the_rocket_equation(self) -> None:
        timeline = two_burn_timeline()
        self.assertEqual([burn.burn_id for burn in timeline.burns], ["upper-1", "upper-2"])
        velocity, mass, ox, fu, entries = integrate(timeline, 30.0, 0.25)
        m1_end = 300.0 - 3.0 * 6.0
        m2_start = m1_end - 30.0
        m2_end = m2_start - 2.0 * 4.0
        expected = C * log(300.0 / m1_end) + C * log(m2_start / m2_end)
        self.assertTrue(isclose(velocity, expected, rel_tol=1e-10))
        self.assertAlmostEqual(mass, m2_end, places=10)
        self.assertAlmostEqual(ox, 60.0 - 0.75 * (18.0 + 8.0), places=10)
        self.assertAlmostEqual(fu, 20.0 - 0.25 * (18.0 + 8.0), places=10)
        self.assertIn("jettison/interstage", entries)
        self.assertEqual(timeline.jettisoned_mass_kg(12.0), 30.0)
        self.assertEqual(timeline.jettisoned_mass_kg(12.0, include_at_time=False), 0.0)

    def test_events_are_namespaced_and_time_ordered(self) -> None:
        events = two_burn_timeline().events()
        timed = [event for event in events if event.time_s is not None]
        self.assertEqual([event.time_s for event in timed], sorted(e.time_s for e in timed))
        self.assertEqual(
            [event.event_id for event in timed],
            [
                "upper-1/burn:start", "upper-1/burn:stop", "jettison/interstage",
                "upper-2/burn:start", "upper-2/burn:stop",
            ],
        )
        state = [event for event in events if event.time_s is None]
        self.assertEqual(
            [event.event_id for event in state],
            ["upper-1/burn:reserve:m_ox", "upper-2/burn:reserve:m_ox"],
        )
        self.assertTrue(all(event.stops_burn for event in state))

    def test_overlap_is_refused_unless_declared_and_then_summed(self) -> None:
        with self.assertRaises(DomainError):
            PropulsionTimeline(
                burns=(
                    ScheduledBurn("a", steady_bridge(100.0, 4.0), 0.0),
                    ScheduledBurn("b", steady_bridge(100.0, 4.0), 3.0),
                )
            )
        touching = PropulsionTimeline(
            burns=(
                ScheduledBurn("a", steady_bridge(100.0, 4.0), 0.0),
                ScheduledBurn("b", steady_bridge(60.0, 4.0), 4.0),
            )
        )
        directions = {"a": (1.0, 0.0, 0.0), "b": (0.0, 1.0, 0.0)}
        states = {"m_ox": 50.0, "m_fu": 10.0}
        left = touching.evaluate(time_s=4.0, vehicle_mass_kg=100.0, directions=directions,
                                 additional_states=states, side="left")
        right = touching.evaluate(time_s=4.0, vehicle_mass_kg=100.0, directions=directions,
                                  additional_states=states, side="right")
        self.assertEqual(left.active_burns, ("a",))
        self.assertEqual(right.active_burns, ("b",))
        concurrent = PropulsionTimeline(
            burns=(
                ScheduledBurn("main", steady_bridge(100.0, 4.0), 0.0),
                ScheduledBurn("rcs", steady_bridge(10.0, 4.0), 1.0),
            ),
            allow_overlapping_burns=True,
        )
        value = concurrent.evaluate(
            time_s=2.0, vehicle_mass_kg=100.0,
            directions={"main": (1.0, 0.0, 0.0), "rcs": (0.0, 1.0, 0.0)},
            additional_states=states,
        )
        self.assertEqual(value.active_burns, ("main", "rcs"))
        self.assertEqual(value.acceleration_m_s2, (1.0, 0.1, 0.0))
        self.assertAlmostEqual(value.mass_derivative_kg_s, -(2.0 + 0.2))
        drains = dict(value.additional_state_derivatives)
        self.assertAlmostEqual(drains["m_ox"] + drains["m_fu"], value.mass_derivative_kg_s)

    def test_directions_stopped_burns_and_ambient_forwarding(self) -> None:
        timeline = PropulsionTimeline(
            burns=(
                ScheduledBurn("booster", steady_bridge(100.0, 4.0, exit_area=0.01), 0.0),
                ScheduledBurn("vacuum", steady_bridge(100.0, 4.0), 10.0),
            )
        )
        states = {"m_ox": 50.0, "m_fu": 10.0}
        coast = timeline.evaluate(time_s=6.0, vehicle_mass_kg=100.0, directions={},
                                  additional_states=states)
        self.assertEqual(coast.active_burns, ())
        with self.assertRaises(DomainError):
            timeline.evaluate(time_s=1.0, vehicle_mass_kg=100.0, directions={},
                              additional_states=states)
        booster = timeline.evaluate(
            time_s=1.0, vehicle_mass_kg=100.0, directions={"booster": (1.0, 0.0, 0.0)},
            additional_states=states, ambient_pressure_pa=1000.0,
        )
        self.assertAlmostEqual(booster.acceleration_m_s2[0], (100.0 - 10.0) / 100.0)
        vacuum = timeline.evaluate(
            time_s=11.0, vehicle_mass_kg=100.0, directions={"vacuum": (1.0, 0.0, 0.0)},
            additional_states=states, ambient_pressure_pa=1000.0,
        )
        self.assertAlmostEqual(vacuum.acceleration_m_s2[0], 1.0)
        stopped = timeline.evaluate(
            time_s=1.0, vehicle_mass_kg=100.0, directions={}, additional_states=states,
            stopped_burns={"booster"},
        )
        self.assertEqual(stopped.active_burns, ())
        with self.assertRaises(DomainError):
            timeline.evaluate(time_s=1.0, vehicle_mass_kg=100.0, directions={},
                              additional_states=states, stopped_burns={"unknown"})
        with self.assertRaises(DomainError):
            timeline.apply_jettison("none", 10.0)
        jettison = PropulsionTimeline(
            burns=(ScheduledBurn("a", steady_bridge(100.0, 4.0), 0.0),),
            jettisons=(MassJettison("big", 5.0, 200.0),),
        )
        with self.assertRaises(DomainError):
            jettison.apply_jettison("big", 150.0)
        with self.assertRaises(DomainError):
            PropulsionTimeline(burns=())
        with self.assertRaises(DomainError):
            PropulsionTimeline(
                burns=(ScheduledBurn("a", steady_bridge(100.0, 4.0), 0.0),),
                jettisons=(MassJettison("a", 5.0, 1.0),),
            )

    def test_persistence_round_trip_and_tamper_detection(self) -> None:
        timeline = two_burn_timeline()
        reopened = PropulsionTimeline.from_json(timeline.to_json())
        self.assertEqual(reopened, timeline)
        self.assertEqual(reopened.timeline_hash, timeline.timeline_hash)
        document = timeline.to_dict()
        for label, mutate in {
            "ignition": lambda d: d["burns"][0].__setitem__("ignition_time_s", 2.5),
            "jettison mass": lambda d: d["jettisons"][0].__setitem__("mass_kg", 31.0),
            "bridge artifact": lambda d: d["bridges"]["upper-1"]["artifact"]["segments"][0]
            .__setitem__("end_axial_thrust_n", 151.0),
            "missing bridge": lambda d: d["bridges"].pop("upper-2"),
        }.items():
            tampered = copy.deepcopy(document)
            mutate(tampered)
            with self.subTest(label), self.assertRaises(InputError):
                PropulsionTimeline.from_dict(tampered)


if __name__ == "__main__":
    unittest.main()
