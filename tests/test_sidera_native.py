"""Native Sidera thrust-provider tests: fail-closed refusal and end-to-end burns.

Sidera is optional. Set ``RPLAB_SIDERA_SRC`` to a Sidera ``src`` directory to
put it on ``sys.path`` for these tests; without it the propagation tests skip
and only the refusal contract runs.
"""

from __future__ import annotations

import math
import os
import sys
import tempfile
import unittest

_SIDERA_SRC = os.environ.get("RPLAB_SIDERA_SRC", "").strip()
if _SIDERA_SRC and os.path.isdir(_SIDERA_SRC) and _SIDERA_SRC not in sys.path:
    # Never write bytecode or Numba JIT caches into the external checkout.
    sys.dont_write_bytecode = True
    os.environ.setdefault(
        "NUMBA_CACHE_DIR", os.path.join(tempfile.gettempdir(), "rplab_sidera_numba_cache")
    )
    sys.path.append(_SIDERA_SRC)

from rocket_propulsion.core.errors import DomainError, FeatureUnavailableError  # noqa: E402
from rocket_propulsion.integrations.sidera_native import (  # noqa: E402
    SideraNativeBurn,
    SideraNativeCapabilities,
    inspect_sidera_native_capabilities,
    missing_sidera_native_capabilities,
)
from rocket_propulsion.propulsion.burns import (  # noqa: E402
    PropagatorPropulsionBridge,
    free_space_burn_metrics,
)
from tests.test_burn_propagator import example_artifact, example_bridge  # noqa: E402
from tests.test_burn_stream_performance import example_stream_surface  # noqa: E402

CAPABILITIES = inspect_sidera_native_capabilities()
SIDERA_SUPPORTED = CAPABILITIES.status == "supported"

# Tolerances of the end-to-end checks (see docs/SIDERA_NATIVE_WRAPPER.md).
MASS_TOLERANCE_KG = 1e-10  # final mass vs exact bridge consumption (6 kg of 100 kg)
DELTA_V_TOLERANCE_M_S = 1e-8  # delivered |dv| vs free-space metric and rocket equation
EVENT_TIME_TOLERANCE_S = 1e-9  # located dry-mass root vs analytic root


def mass_only_bridge(**kwargs: object) -> PropagatorPropulsionBridge:
    return PropagatorPropulsionBridge(
        example_artifact(), require_complete_stream_binding=False, **kwargs
    )


def synthetic_capabilities(status: str) -> SideraNativeCapabilities:
    supported = status == "supported"
    return SideraNativeCapabilities(
        status=status,
        version="test" if supported else None,
        propagate_available=supported,
        dynamics_protocol_available=supported,
        extra_events_available=supported,
        frame_basis_available=supported,
    )


class SideraNativeRefusalTests(unittest.TestCase):
    def test_construction_refuses_cleanly_without_sidera(self) -> None:
        if CAPABILITIES.status != "unavailable":
            self.skipTest("Sidera is importable in this run.")
        self.assertIsNone(CAPABILITIES.version)
        with self.assertRaises(FeatureUnavailableError) as caught:
            SideraNativeBurn(
                object(), mass_only_bridge(), ignition_epoch_s=0.0, direction=(1.0, 0.0, 0.0)
            )
        self.assertIn("propagate_available", caught.exception.details["missing_capabilities"])

    def test_unavailable_capabilities_refuse_in_every_mode(self) -> None:
        with self.assertRaises(FeatureUnavailableError) as caught:
            SideraNativeBurn(
                object(),
                mass_only_bridge(),
                ignition_epoch_s=0.0,
                direction=(1.0, 0.0, 0.0),
                capabilities=synthetic_capabilities("unavailable"),
            )
        details = caught.exception.details
        self.assertEqual(details["sidera_status"], "unavailable")
        self.assertEqual(
            details["missing_capabilities"],
            [
                "propagate_available",
                "dynamics_protocol_available",
                "extra_events_available",
                "frame_basis_available",
            ],
        )

    def test_states_conditions_and_ambient_pressure_fail_closed(self) -> None:
        supported = synthetic_capabilities("supported")
        self.assertEqual(missing_sidera_native_capabilities(supported, mass_only_bridge()), ())
        self.assertEqual(
            missing_sidera_native_capabilities(supported, example_bridge()),
            ("additional_states_available",),
        )
        surface_bridge = mass_only_bridge(stream_performance_surface=example_stream_surface())
        self.assertEqual(
            missing_sidera_native_capabilities(supported, surface_bridge),
            ("operating_conditions_available",),
        )
        nozzle_bridge = mass_only_bridge(nozzle_exit_area_m2=0.01)
        self.assertEqual(
            missing_sidera_native_capabilities(supported, nozzle_bridge),
            ("ambient_pressure_available",),
        )
        with self.assertRaises(FeatureUnavailableError):
            SideraNativeBurn(
                object(),
                example_bridge(),
                ignition_epoch_s=0.0,
                direction=(1.0, 0.0, 0.0),
                capabilities=supported,
            )

    def test_non_bridge_is_a_domain_error(self) -> None:
        with self.assertRaises(DomainError):
            SideraNativeBurn(
                object(), object(), ignition_epoch_s=0.0, direction=(1.0, 0.0, 0.0)  # type: ignore[arg-type]
            )


@unittest.skipUnless(SIDERA_SUPPORTED, "Sidera is not on the explicit test path (RPLAB_SIDERA_SRC).")
class SideraNativePropagationTests(unittest.TestCase):
    """Short finite burns flown through Sidera's public ``propagate``.

    A circular Earth orbit at 350 000 km altitude keeps the gravity-gradient
    difference between the burn and a coast reference below 1e-10 m/s over
    the run, so the measured velocity change is the free-space delta-v.
    """

    INITIAL_MASS_KG = 100.0
    IGNITION_S = 10.0

    @classmethod
    def setUpClass(cls) -> None:
        import numpy
        from sidera.core.body_engine.orbits import build_body_dynamics
        from sidera.core.config import load_body_config
        from sidera.core.propagation import propagate
        from sidera.spacecraft.properties import SpacecraftProps

        cls.numpy = numpy
        cls.propagate = staticmethod(propagate)
        cls.config = load_body_config(
            "earth",
            periapsis_altitude_m=3.5e8,
            duration_s=15.0,
            output_dt_s=1.0,
            spacecraft=SpacecraftProps(mass_kg=cls.INITIAL_MASS_KG),
        )
        cls.engine = build_body_dynamics(cls.config)
        cls.y0 = numpy.append(cls.config.initial_state.to_array(), cls.INITIAL_MASS_KG)
        coast = propagate(cls.engine, cls.y0, cls.config.propagator, time_cfg=cls.config.time)
        cls.coast_final = numpy.asarray(coast.y[-1], dtype=numpy.float64)

    def fly(self, bridge: PropagatorPropulsionBridge, **kwargs: object):
        burn = SideraNativeBurn(
            self.engine,
            bridge,
            ignition_epoch_s=self.IGNITION_S,
            direction=kwargs.pop("direction", (1.0, 0.0, 0.0)),  # type: ignore[arg-type]
            **kwargs,  # type: ignore[arg-type]
        )
        run = burn.propagate(self.y0, self.config.propagator, time_cfg=self.config.time)
        return burn, run

    def test_finite_burn_closes_mass_delta_v_and_honours_events(self) -> None:
        numpy = self.numpy
        bridge = mass_only_bridge()
        burn, run = self.fly(bridge)
        from sidera.core.propagation.dynamics_protocol import PropagationDynamics

        self.assertIsInstance(burn.segment_dynamics(0), PropagationDynamics)
        metrics = free_space_burn_metrics(bridge.artifact, self.INITIAL_MASS_KG)
        times, states = run.stitched()

        # Mass: initial minus the bridge's exactly integrated consumption.
        expected_final = self.INITIAL_MASS_KG - metrics.consumed_propellant_kg
        self.assertTrue(run.completed)
        self.assertAlmostEqual(float(states[-1, 6]), expected_final, delta=MASS_TOLERANCE_KG)
        self.assertAlmostEqual(
            run.consumed_propellant_kg, metrics.consumed_propellant_kg, delta=MASS_TOLERANCE_KG
        )

        # Delta-v: burn minus coast reference, against both comparators.
        delta_v = states[-1, 3:6] - self.coast_final[3:6]
        exhaust = {
            segment.start_axial_thrust_n / segment.start_total_mass_flow_kg_s
            for segment in bridge.artifact.segments
        } | {
            segment.end_axial_thrust_n / segment.end_total_mass_flow_kg_s
            for segment in bridge.artifact.segments
        }
        self.assertEqual(len(exhaust), 1)  # constant F/q: the rocket equation is exact
        rocket = exhaust.pop() * math.log(self.INITIAL_MASS_KG / expected_final)
        self.assertAlmostEqual(float(numpy.linalg.norm(delta_v)), metrics.delta_v_m_s,
                               delta=DELTA_V_TOLERANCE_M_S)
        self.assertAlmostEqual(float(delta_v[0]), rocket, delta=DELTA_V_TOLERANCE_M_S)
        self.assertAlmostEqual(float(delta_v[1]), 0.0, delta=DELTA_V_TOLERANCE_M_S)
        self.assertAlmostEqual(float(delta_v[2]), 0.0, delta=DELTA_V_TOLERANCE_M_S)

        # Events: start, interior boundary and stop at the canonical dates.
        self.assertEqual(
            [(event.event_id, event.time_s) for event in run.events],
            [("burn:start", 10.0), ("burn:boundary:1", 11.0), ("burn:stop", 12.0)],
        )
        self.assertEqual(
            [label for label, _ in run.segments],
            ["coast:pre-burn", "burn:segment:1", "burn:segment:2", "coast:post-burn"],
        )
        before = times <= self.IGNITION_S
        after = times >= 12.0
        self.assertTrue(numpy.all(states[before, 6] == self.INITIAL_MASS_KG))
        self.assertTrue(numpy.all(states[after, 6] == states[-1, 6]))
        self.assertAlmostEqual(run.events[1].vehicle_mass_kg, 97.0, delta=MASS_TOLERANCE_KG)
        mapping = {item.event_id: item.mechanism for item in burn.event_mapping()}
        self.assertEqual(
            mapping,
            {
                "burn:start": "segment_boundary",
                "burn:boundary:1": "segment_boundary",
                "burn:stop": "segment_boundary",
            },
        )

    def test_ignition_bias_moves_every_time_event(self) -> None:
        _, run = self.fly(mass_only_bridge(), parameter_overrides={"ignition_time_bias_s": 0.5})
        self.assertEqual(
            [event.time_s for event in run.events], [10.5, 11.5, 12.5]
        )
        _, states = run.stitched()
        self.assertAlmostEqual(float(states[-1, 6]), 94.0, delta=MASS_TOLERANCE_KG)

    def test_dry_mass_event_stops_thrust_and_run_coasts_on(self) -> None:
        numpy = self.numpy
        bridge = mass_only_bridge(protected_dry_mass_kg=96.5)
        burn, run = self.fly(bridge)
        self.assertEqual(
            [item.mechanism for item in burn.event_mapping()][-1], "terminal_event"
        )
        # 3 kg by t = 1 s, then 3 kg/s: 0.5 kg more after 1/6 s.
        expected_cut = self.IGNITION_S + 1.0 + 0.5 / 3.0
        self.assertFalse(run.completed)
        self.assertFalse(run.stopped_early)
        self.assertEqual(run.events[-1].event_id, "burn:dry-mass")
        self.assertAlmostEqual(run.cutoff_date_s, expected_cut, delta=EVENT_TIME_TOLERANCE_S)
        self.assertAlmostEqual(run.mass_at_cutoff_kg, 96.5, delta=MASS_TOLERANCE_KG)
        times, states = run.stitched()
        self.assertAlmostEqual(float(times[-1]), 15.0)
        self.assertTrue(numpy.all(states[times >= expected_cut + 1e-6, 6] == states[-1, 6]))

    def test_vnb_direction_is_resolved_from_the_state(self) -> None:
        numpy = self.numpy
        _, run = self.fly(mass_only_bridge(), direction=(1.0, 0.0, 0.0), frame="vnb")
        _, states = run.stitched()
        self.assertAlmostEqual(float(states[-1, 6]), 94.0, delta=MASS_TOLERANCE_KG)
        delta_v = states[-1, 3:6] - self.coast_final[3:6]
        velocity = self.coast_final[3:6] / numpy.linalg.norm(self.coast_final[3:6])
        # Prograde: the velocity change lies along the velocity direction.
        self.assertGreater(float(numpy.dot(delta_v, velocity)), 0.999 * numpy.linalg.norm(delta_v))


if __name__ == "__main__":
    unittest.main()
