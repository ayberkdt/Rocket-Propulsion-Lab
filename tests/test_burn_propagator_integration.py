"""End-to-end time integration of the propagator bridge with event handling.

The derivative tests in ``test_burn_propagator.py`` check one evaluation at a
time.  These tests integrate the bridge with a small event-aware RK4 loop that
behaves like a consuming orbit propagator: it lands on every time event, uses
the left limit at a step end and the right limit at a step start, locates
inventory roots by bisection and latches the burn off after a stopping event.
"""

import random
import unittest
from itertools import pairwise
from math import isclose, log, sqrt

from rocket_propulsion.core.errors import DomainError
from rocket_propulsion.propulsion.burns import (
    AdditionalStateConstraint,
    BurnEventKind,
    PerformanceAxis,
    PerformanceScalePoint,
    PropagatorPropulsionBridge,
    RectilinearPerformanceSurface,
    StreamStateBinding,
    TabulatedBurnArtifact,
    TabulatedBurnSegment,
)

STATES = ("oxidizer_mass_kg", "fuel_mass_kg")
EXHAUST_VELOCITY_M_S = 50.0  # axial thrust / total flow in every test profile


def ramp_steady_artifact() -> TabulatedBurnArtifact:
    return TabulatedBurnArtifact(
        segments=(
            TabulatedBurnSegment(
                0.0, 1.0, 120.0, 240.0, 100.0, 200.0, 2.0, 4.0,
                (("oxidizer", 1.5), ("fuel", 0.5)),
                (("oxidizer", 3.0), ("fuel", 1.0)),
                "ramp",
            ),
            TabulatedBurnSegment(
                1.0, 2.0, 180.0, 180.0, 150.0, 150.0, 3.0, 3.0,
                (("oxidizer", 2.25), ("fuel", 0.75)),
                (("oxidizer", 2.25), ("fuel", 0.75)),
                "steady",
            ),
        ),
        source_id="integration-test",
        source_sha256="b" * 64,
        reference_ids=("NASA-20040000363",),
    )


def steady_artifact() -> TabulatedBurnArtifact:
    return TabulatedBurnArtifact(
        segments=(
            TabulatedBurnSegment(
                0.0, 2.0, 180.0, 180.0, 150.0, 150.0, 3.0, 3.0,
                (("oxidizer", 2.25), ("fuel", 0.75)),
                (("oxidizer", 2.25), ("fuel", 0.75)),
                "steady",
            ),
        ),
        source_id="integration-steady",
        source_sha256="c" * 64,
        reference_ids=("NASA-20040000363",),
    )


def bridge_for(
    artifact: TabulatedBurnArtifact, *, oxidizer_reserve_kg: float = 2.0
) -> PropagatorPropulsionBridge:
    return PropagatorPropulsionBridge(
        artifact,
        stream_bindings=(
            StreamStateBinding("oxidizer", "oxidizer_mass_kg"),
            StreamStateBinding("fuel", "fuel_mass_kg"),
        ),
        state_constraints=(
            AdditionalStateConstraint("oxidizer_mass_kg", oxidizer_reserve_kg),
            AdditionalStateConstraint("fuel_mass_kg", 1.0),
        ),
        protected_dry_mass_kg=50.0,
    )


class EventAwareRk4:
    """Minimal fixed-step RK4 propagator with exact event landing."""

    def __init__(
        self,
        bridge: PropagatorPropulsionBridge,
        direction: tuple[float, float, float],
        *,
        ignition_time_bias_s: float = 0.0,
        apply_inventory_limits: bool = False,
    ) -> None:
        self.bridge = bridge
        self.apply_inventory_limits = apply_inventory_limits
        self.direction = direction
        self.overrides = {"ignition_time_bias_s": ignition_time_bias_s}
        self.events = bridge.event_surfaces()
        self.time_events = sorted(
            event.relative_time_s + ignition_time_bias_s
            for event in self.events
            if event.relative_time_s is not None
        )
        self.burn_latched_off = False
        self.stop_events: list[tuple[str, float]] = []

    def derivative(self, time_s: float, y: list[float], side: str) -> list[float]:
        velocity = y[3:6]
        if self.burn_latched_off:
            return [*velocity, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
        value = self.bridge.evaluate(
            relative_time_s=time_s,
            vehicle_mass_kg=y[6],
            direction=self.direction,
            additional_states=dict(zip(STATES, y[7:9], strict=True)),
            parameter_overrides=self.overrides,
            side=side,
            apply_inventory_limits=self.apply_inventory_limits,
        )
        drains = dict(value.additional_state_derivatives)
        return [
            *velocity,
            *value.acceleration_m_s2,
            value.mass_derivative_kg_s,
            drains[STATES[0]],
            drains[STATES[1]],
        ]

    def step(self, time_s: float, y: list[float], h: float) -> list[float]:
        def shifted(base: list[float], rate: list[float], factor: float) -> list[float]:
            return [a + factor * b for a, b in zip(base, rate, strict=True)]

        k1 = self.derivative(time_s, y, "right")
        k2 = self.derivative(time_s + h / 2.0, shifted(y, k1, h / 2.0), "right")
        k3 = self.derivative(time_s + h / 2.0, shifted(y, k2, h / 2.0), "right")
        k4 = self.derivative(time_s + h, shifted(y, k3, h), "left")
        return [
            value + h / 6.0 * (a + 2.0 * b + 2.0 * c + d)
            for value, a, b, c, d in zip(y, k1, k2, k3, k4, strict=True)
        ]

    def stopping_values(self, time_s: float, y: list[float]) -> dict[str, float]:
        states = dict(zip(STATES, y[7:9], strict=True))
        return {
            event.event_id: event.value(
                relative_time_s=time_s,
                vehicle_mass_kg=y[6],
                additional_states=states,
                ignition_time_bias_s=self.overrides["ignition_time_bias_s"],
            )
            for event in self.events
            if event.stops_burn and event.relative_time_s is None
        }

    def propagate(
        self, start_s: float, end_s: float, y0: list[float], max_step_s: float
    ) -> list[float]:
        time_s, y = start_s, list(y0)
        while time_s < end_s - 1e-15:
            target = min(end_s, time_s + max_step_s)
            for event_time in self.time_events:
                if time_s + 1e-15 < event_time < target:
                    target = event_time
                    break
            h = target - time_s
            candidate = self.step(time_s, y, h)
            if not self.burn_latched_off:
                crossed = [
                    name
                    for name, value in self.stopping_values(target, candidate).items()
                    if value <= 0.0
                ]
                if crossed:
                    low, high = 0.0, h
                    for _ in range(80):
                        middle = 0.5 * (low + high)
                        trial = self.step(time_s, y, middle)
                        values = self.stopping_values(time_s + middle, trial)
                        if min(values.values()) > 0.0:
                            low = middle
                        else:
                            high = middle
                    h = high
                    candidate = self.step(time_s, y, h)
                    self.burn_latched_off = True
                    self.stop_events.append((crossed[0], time_s + h))
            time_s, y = time_s + h, candidate
        return y


def initial_state(oxidizer_kg: float = 10.0, fuel_kg: float = 5.0) -> list[float]:
    return [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 100.0, oxidizer_kg, fuel_kg]


class BridgeTimeIntegrationTests(unittest.TestCase):
    direction = (3.0, 4.0, 0.0)
    unit = (0.6, 0.8, 0.0)

    def speed_along_axis(self, y: list[float]) -> float:
        return sum(component * axis for component, axis in zip(y[3:6], self.unit, strict=True))

    def test_full_burn_reproduces_rocket_equation_and_mass_closure(self) -> None:
        bridge = bridge_for(ramp_steady_artifact())
        propagator = EventAwareRk4(bridge, self.direction)
        final = propagator.propagate(-0.5, 3.0, initial_state(), 1e-2)

        consumed = bridge.artifact.consumed_propellant_kg
        self.assertAlmostEqual(consumed, 6.0)
        self.assertAlmostEqual(final[6], 100.0 - consumed, places=12)
        self.assertAlmostEqual(final[7], 10.0 - 4.5, places=12)
        self.assertAlmostEqual(final[8], 5.0 - 1.5, places=12)
        expected_dv = EXHAUST_VELOCITY_M_S * log(100.0 / final[6])
        self.assertTrue(isclose(self.speed_along_axis(final), expected_dv, rel_tol=1e-10))
        transverse = sqrt(
            sum(value * value for value in final[3:6]) - self.speed_along_axis(final) ** 2
        )
        self.assertLess(transverse, 1e-9)
        closure = bridge.mass_closure(
            vehicle_mass_kg=final[6],
            additional_states=dict(zip(STATES, final[7:9], strict=True)),
            non_propellant_mass_kg=85.0,
        )
        self.assertTrue(closure.is_closed(absolute_tolerance_kg=1e-11))
        self.assertEqual(propagator.stop_events, [])

    def test_steady_burn_position_matches_closed_form(self) -> None:
        bridge = bridge_for(steady_artifact())
        final = EventAwareRk4(bridge, self.direction).propagate(
            0.0, 2.0, initial_state(), 5e-3
        )
        flow, m0, duration = 3.0, 100.0, 2.0
        remaining = m0 - flow * duration
        expected_position = EXHAUST_VELOCITY_M_S * (
            duration + remaining / flow * log(remaining / m0)
        )
        position = sum(
            component * axis for component, axis in zip(final[0:3], self.unit, strict=True)
        )
        self.assertTrue(isclose(position, expected_position, rel_tol=1e-10))

    def test_reserve_root_cuts_burn_at_analytic_time_and_mass(self) -> None:
        # 3 kg of oxidizer above the reserve: 2.25 kg in the ramp, then 0.75 kg
        # at 2.25 kg/s, so the reserve is reached at t = 4/3 s.
        bridge = bridge_for(ramp_steady_artifact())
        propagator = EventAwareRk4(bridge, self.direction)
        final = propagator.propagate(0.0, 3.0, initial_state(oxidizer_kg=5.0), 1e-2)

        self.assertEqual(len(propagator.stop_events), 1)
        event_id, cut_time = propagator.stop_events[0]
        self.assertEqual(event_id, "burn:reserve:oxidizer_mass_kg")
        self.assertAlmostEqual(cut_time, 4.0 / 3.0, places=10)
        self.assertAlmostEqual(final[7], 2.0, places=10)
        self.assertAlmostEqual(final[8], 4.0, places=10)
        self.assertAlmostEqual(final[6], 96.0, places=10)
        expected_dv = EXHAUST_VELOCITY_M_S * log(100.0 / final[6])
        self.assertTrue(isclose(self.speed_along_axis(final), expected_dv, rel_tol=1e-10))

    def test_internal_inventory_inhibition_biases_rk_root_location(self) -> None:
        # Documents why event-handling consumers pass apply_inventory_limits=False:
        # the internal clamp zeroes the final RK stage past the root.
        bridge = bridge_for(ramp_steady_artifact())
        propagator = EventAwareRk4(bridge, self.direction, apply_inventory_limits=True)
        propagator.propagate(0.0, 3.0, initial_state(oxidizer_kg=5.0), 1e-2)
        _, cut_time = propagator.stop_events[0]
        self.assertGreater(abs(cut_time - 4.0 / 3.0), 1e-5)

    def test_ignition_bias_shifts_timeline_without_changing_the_burn(self) -> None:
        bridge = bridge_for(ramp_steady_artifact())
        nominal = EventAwareRk4(bridge, self.direction).propagate(
            -0.5, 3.0, initial_state(), 1e-2
        )
        biased_propagator = EventAwareRk4(bridge, self.direction, ignition_time_bias_s=0.25)
        self.assertEqual(biased_propagator.time_events, [0.25, 1.25, 2.25])
        biased = biased_propagator.propagate(-0.5, 3.0, initial_state(), 1e-2)
        for index in range(3, 9):
            self.assertAlmostEqual(biased[index], nominal[index], places=10)


class DirectionJacobianTests(unittest.TestCase):
    def evaluate(self, direction: tuple[float, float, float]):
        return bridge_for(ramp_steady_artifact()).evaluate(
            relative_time_s=0.4,
            vehicle_mass_kg=90.0,
            direction=direction,
            additional_states={"oxidizer_mass_kg": 10.0, "fuel_mass_kg": 5.0},
        )

    def test_direction_jacobian_matches_central_finite_differences(self) -> None:
        for direction in ((3.0, 4.0, 0.0), (0.3, -1.2, 2.0)):
            jacobian = self.evaluate(direction).partials.acceleration_wrt_direction
            step = 1e-6
            for column in range(3):
                plus = list(direction)
                minus = list(direction)
                plus[column] += step
                minus[column] -= step
                upper = self.evaluate(tuple(plus)).acceleration_m_s2  # type: ignore[arg-type]
                lower = self.evaluate(tuple(minus)).acceleration_m_s2  # type: ignore[arg-type]
                for row in range(3):
                    finite = (upper[row] - lower[row]) / (2.0 * step)
                    self.assertAlmostEqual(jacobian[row][column], finite, places=7)

    def test_direction_jacobian_annihilates_the_thrust_axis(self) -> None:
        direction = (0.3, -1.2, 2.0)
        jacobian = self.evaluate(direction).partials.acceleration_wrt_direction
        for row in jacobian:
            self.assertAlmostEqual(sum(a * b for a, b in zip(row, direction, strict=True)), 0.0)

    def test_inactive_burn_has_zero_direction_jacobian(self) -> None:
        inactive = bridge_for(ramp_steady_artifact()).evaluate(
            relative_time_s=5.0,
            vehicle_mass_kg=90.0,
            direction=(1.0, 0.0, 0.0),
            additional_states={"oxidizer_mass_kg": 10.0, "fuel_mass_kg": 5.0},
        )
        self.assertFalse(inactive.active)
        self.assertEqual(inactive.partials.acceleration_wrt_direction, ((0.0,) * 3,) * 3)


class MassClosureTests(unittest.TestCase):
    def test_residual_and_refusals(self) -> None:
        bridge = bridge_for(ramp_steady_artifact())
        closure = bridge.mass_closure(
            vehicle_mass_kg=100.0,
            additional_states={"oxidizer_mass_kg": 10.0, "fuel_mass_kg": 4.0},
            non_propellant_mass_kg=85.0,
        )
        self.assertAlmostEqual(closure.tracked_propellant_kg, 14.0)
        self.assertAlmostEqual(closure.residual_kg, 1.0)
        self.assertAlmostEqual(closure.relative_residual, 0.01)
        self.assertFalse(closure.is_closed(absolute_tolerance_kg=0.5))
        self.assertTrue(closure.is_closed(absolute_tolerance_kg=1.0))
        with self.assertRaises(DomainError):
            closure.is_closed(absolute_tolerance_kg=-1.0)
        with self.assertRaises(DomainError):
            bridge.mass_closure(
                vehicle_mass_kg=100.0,
                additional_states={"oxidizer_mass_kg": 10.0},
                non_propellant_mass_kg=85.0,
            )
        with self.assertRaises(DomainError):
            bridge.mass_closure(
                vehicle_mass_kg=0.0,
                additional_states={"oxidizer_mass_kg": 10.0, "fuel_mass_kg": 4.0},
                non_propellant_mass_kg=0.0,
            )


def _reference_artifact_segment(artifact: TabulatedBurnArtifact, time_s: float, side: str):
    """The original linear scan, kept as the oracle for the bisect lookup."""

    tolerance = max(1e-12, artifact.duration_s * 1e-12)
    for index, segment in enumerate(artifact.segments):
        if segment.start_time_s <= time_s < segment.end_time_s:
            if side == "left" and index > 0 and abs(time_s - segment.start_time_s) <= tolerance:
                return artifact.segments[index - 1]
            return segment
    return artifact.segments[-1]


def _reference_bridge_segment(artifact: TabulatedBurnArtifact, time_s: float, side: str):
    tolerance = max(1e-12, artifact.duration_s * 1e-12)
    at_boundary = any(
        abs(time_s - value) <= tolerance
        for value in (0.0, *(segment.end_time_s for segment in artifact.segments))
    )
    if time_s < 0.0 or time_s > artifact.duration_s:
        return None, at_boundary
    if abs(time_s) <= tolerance and side == "left":
        return None, True
    if abs(time_s - artifact.duration_s) <= tolerance and side == "right":
        return None, True
    for index, segment in enumerate(artifact.segments):
        if segment.start_time_s <= time_s < segment.end_time_s:
            if side == "left" and index > 0 and abs(time_s - segment.start_time_s) <= tolerance:
                return artifact.segments[index - 1], True
            return segment, at_boundary
    return artifact.segments[-1], at_boundary


class SegmentLookupTests(unittest.TestCase):
    def setUp(self) -> None:
        generator = random.Random(20261007)
        boundaries = [0.0]
        for _ in range(60):
            boundaries.append(boundaries[-1] + generator.uniform(0.01, 2.0))
        segments = []
        for index, (start, end) in enumerate(pairwise(boundaries)):
            thrust = 100.0 + index
            flow = thrust / EXHAUST_VELOCITY_M_S
            segments.append(
                TabulatedBurnSegment(
                    start, end, thrust, thrust + 1.0, thrust, thrust + 1.0, flow, flow,
                    phase=f"segment-{index}",
                )
            )
        self.artifact = TabulatedBurnArtifact(
            segments=tuple(segments),
            source_id="lookup-test",
            source_sha256="d" * 64,
            reference_ids=("NASA-20040000363",),
        )
        self.bridge = PropagatorPropulsionBridge(self.artifact)
        tolerance = self.artifact.boundary_tolerance_s
        self.times = [generator.uniform(-1.0, boundaries[-1] + 1.0) for _ in range(400)]
        for boundary in boundaries:
            for offset in (0.0, -0.5 * tolerance, 0.5 * tolerance, -2 * tolerance, 2 * tolerance):
                self.times.append(boundary + offset)

    def test_bridge_lookup_matches_linear_scan_oracle(self) -> None:
        for time_s in self.times:
            for side in ("left", "right"):
                expected = _reference_bridge_segment(self.artifact, time_s, side)
                self.assertEqual(self.bridge._segment_at(time_s, side), expected, (time_s, side))

    def test_artifact_lookup_matches_oracle_and_never_leaves_a_segment(self) -> None:
        for time_s in self.times:
            if not 0.0 <= time_s <= self.artifact.duration_s:
                continue
            for side in ("left", "right"):
                expected = _reference_artifact_segment(self.artifact, time_s, side)
                index = self.artifact.segment_index_at(time_s, side=side)
                self.assertIs(self.artifact.segments[index], expected)
                state = self.artifact.state_at(time_s, side=side)
                self.assertEqual(state["phase"], expected.phase)


class PerformanceGridLineTests(unittest.TestCase):
    def setUp(self) -> None:
        pressures = (1.0e6, 1.5e6, 2.0e6, 3.0e6)
        throttles = (0.4, 0.7, 1.0)
        points = tuple(
            PerformanceScalePoint(
                (pressure, throttle),
                force_scale=throttle**1.5 * (pressure / 2.0e6) ** (0.5 + 0.1 * i),
                mass_flow_scale=throttle * (pressure / 2.0e6) ** 0.5,
            )
            for i, pressure in enumerate(pressures)
            for throttle in throttles
        )
        self.surface = RectilinearPerformanceSurface(
            axes=(
                PerformanceAxis("supply_pressure_pa", "Pa", pressures),
                PerformanceAxis("throttle", "1", throttles),
            ),
            points=points,
            source_id="grid-line-test",
            source_sha256="e" * 64,
        )

    def test_interior_lines_are_non_stopping_signed_roots(self) -> None:
        lines = self.surface.grid_line_surfaces()
        self.assertEqual(
            [(line.condition_name, line.coordinate) for line in lines],
            [
                ("supply_pressure_pa", 1.5e6),
                ("supply_pressure_pa", 2.0e6),
                ("throttle", 0.7),
            ],
        )
        self.assertTrue(all(not line.stops_burn for line in lines))
        conditions = {"supply_pressure_pa": 1.8e6, "throttle": 0.7}
        self.assertAlmostEqual(lines[0].value(conditions), 0.3e6)
        self.assertAlmostEqual(lines[1].value(conditions), -0.2e6)
        self.assertEqual(lines[2].value(conditions), 0.0)
        with self.assertRaises(DomainError):
            lines[0].value({"throttle": 0.5})

    def test_gradient_jumps_across_each_registered_line(self) -> None:
        delta = 1e-3
        for line in self.surface.grid_line_surfaces():
            base = {"supply_pressure_pa": 2.5e6, "throttle": 0.85}
            if line.condition_name == "throttle":
                base["supply_pressure_pa"] = 1.25e6
            below = dict(base, **{line.condition_name: line.coordinate * (1.0 - delta)})
            above = dict(base, **{line.condition_name: line.coordinate * (1.0 + delta)})
            gradient_below = dict(self.surface.evaluate(below).force_scale_gradient)
            gradient_above = dict(self.surface.evaluate(above).force_scale_gradient)
            self.assertGreater(
                abs(gradient_above[line.condition_name] - gradient_below[line.condition_name]),
                1e-3 * abs(gradient_below[line.condition_name]),
            )

    def test_bridge_exposes_grid_lines_only_with_a_surface(self) -> None:
        self.assertEqual(bridge_for(ramp_steady_artifact()).operating_condition_grid_surfaces(), ())
        bridge = PropagatorPropulsionBridge(
            ramp_steady_artifact(),
            performance_surface=self.surface,
            stream_bindings=(
                StreamStateBinding("oxidizer", "oxidizer_mass_kg"),
                StreamStateBinding("fuel", "fuel_mass_kg"),
            ),
        )
        self.assertEqual(
            bridge.operating_condition_grid_surfaces(), self.surface.grid_line_surfaces()
        )
        stopping = [event.kind for event in bridge.event_surfaces() if event.stops_burn]
        self.assertEqual(stopping, [BurnEventKind.STOP])


if __name__ == "__main__":
    unittest.main()
