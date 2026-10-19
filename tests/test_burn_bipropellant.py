"""Dual-tank, dual-line bipropellant coupling and Jacobian tests."""

from __future__ import annotations

import unittest
from dataclasses import astuple, replace
from math import fsum

from rocket_propulsion.core.errors import DomainError
from rocket_propulsion.propulsion.burns import (
    AdditionalStateConstraint,
    BipropellantDynamicFeedConfiguration,
    BipropellantDynamicFeedCoupling,
    BipropellantDynamicFeedState,
    BipropellantMixtureRatioEventKind,
    BipropellantMixtureRatioEventSurface,
    DynamicFeedLine,
    DynamicFeedLineConfiguration,
    DynamicFeedLineState,
    PropagatorPropulsionBridge,
    RegulatedFeedSystem,
    RegulatedFeedState,
    StreamStateBinding,
)
from tests.test_burn_feed_system import example_feed_system
from tests.test_burn_propagator import example_artifact
from tests.test_burn_stream_performance import example_stream_surface


class BipropellantDynamicFeedTests(unittest.TestCase):
    def setUp(self) -> None:
        base_feed = example_feed_system()
        self.oxidizer_feed = RegulatedFeedSystem(
            replace(
                base_feed.configuration,
                propellant_density_kg_m3=1_200.0,
                feed_line_resistance_pa_s2_kg2=0.0,
            )
        )
        self.fuel_feed = RegulatedFeedSystem(
            replace(
                base_feed.configuration,
                propellant_density_kg_m3=800.0,
                feed_line_resistance_pa_s2_kg2=0.0,
            )
        )
        oxidizer_feed_state = self.oxidizer_feed.initial_state(
            propellant_mass_kg=10.0,
            tank_pressure_pa=2.0e6,
            tank_gas_temperature_k=300.0,
            bottle_pressure_pa=20.0e6,
            bottle_gas_temperature_k=300.0,
            regulator_opening=0.5,
        )
        fuel_feed_state = self.fuel_feed.initial_state(
            propellant_mass_kg=6.0,
            tank_pressure_pa=1.9e6,
            tank_gas_temperature_k=300.0,
            bottle_pressure_pa=20.0e6,
            bottle_gas_temperature_k=300.0,
            regulator_opening=0.5,
        )
        self.oxidizer_engine_flow = 2.25 * 1.06
        self.fuel_engine_flow = 0.75 * 1.16
        self.oxidizer_line = DynamicFeedLine(
            DynamicFeedLineConfiguration(
                inertance_pa_s2_kg=1.0e6,
                compliance_kg_pa=1.0e-8,
                resistance_pa_s2_kg2=(2.0e6 - 1.8e6) / self.oxidizer_engine_flow**2,
                minimum_manifold_pressure_pa=1.0e6,
                maximum_manifold_pressure_pa=2.2e6,
                maximum_absolute_flow_kg_s=5.0,
            )
        )
        self.fuel_line = DynamicFeedLine(
            DynamicFeedLineConfiguration(
                inertance_pa_s2_kg=1.2e6,
                compliance_kg_pa=2.0e-8,
                resistance_pa_s2_kg2=(1.9e6 - 1.6e6) / self.fuel_engine_flow**2,
                minimum_manifold_pressure_pa=1.0e6,
                maximum_manifold_pressure_pa=2.2e6,
                maximum_absolute_flow_kg_s=3.0,
            )
        )
        self.state = BipropellantDynamicFeedState(
            oxidizer_feed=oxidizer_feed_state,
            oxidizer_line=DynamicFeedLineState(self.oxidizer_engine_flow, 1.8e6),
            fuel_feed=fuel_feed_state,
            fuel_line=DynamicFeedLineState(self.fuel_engine_flow, 1.6e6),
        )
        self.bridge = PropagatorPropulsionBridge(
            example_artifact(),
            stream_performance_surface=example_stream_surface(),
            stream_bindings=(
                StreamStateBinding("oxidizer", "oxidizer_mass_kg"),
                StreamStateBinding("fuel", "fuel_mass_kg"),
            ),
            state_constraints=(
                AdditionalStateConstraint("oxidizer_mass_kg", 1.0),
                AdditionalStateConstraint("fuel_mass_kg", 1.0),
            ),
        )
        self.configuration = BipropellantDynamicFeedConfiguration(
            oxidizer_state_name="oxidizer_mass_kg",
            fuel_state_name="fuel_mass_kg",
            oxidizer_pressure_condition_name="oxidizer_pressure_pa",
            fuel_pressure_condition_name="fuel_pressure_pa",
            minimum_engine_mixture_ratio=2.5,
            maximum_engine_mixture_ratio=3.0,
        )
        self.coupling = BipropellantDynamicFeedCoupling(
            self.configuration,
            self.bridge,
            self.oxidizer_feed,
            self.oxidizer_line,
            self.fuel_feed,
            self.fuel_line,
        )

    def evaluate(
        self,
        state: BipropellantDynamicFeedState | None = None,
        *,
        relative_time_s: float = 0.5,
    ):
        return self.coupling.evaluate(
            self.state if state is None else state,
            relative_time_s=relative_time_s,
            vehicle_mass_kg=100.0,
            direction=(1.0, 0.0, 0.0),
        )

    def test_two_feed_legs_close_streams_storage_and_vehicle_mass(self) -> None:
        value = self.evaluate()
        self.assertTrue(value.propulsion.active)
        self.assertAlmostEqual(value.oxidizer_engine_flow_kg_s, self.oxidizer_engine_flow)
        self.assertAlmostEqual(value.fuel_engine_flow_kg_s, self.fuel_engine_flow)
        self.assertAlmostEqual(
            value.engine_mixture_ratio,
            self.oxidizer_engine_flow / self.fuel_engine_flow,
        )
        self.assertAlmostEqual(
            value.line_mixture_ratio,
            self.oxidizer_engine_flow / self.fuel_engine_flow,
        )
        self.assertAlmostEqual(value.oxidizer_line.derivative.mass_flow_kg_s2, 0.0)
        self.assertAlmostEqual(value.oxidizer_line.derivative.manifold_pressure_pa_s, 0.0)
        self.assertAlmostEqual(value.fuel_line.derivative.mass_flow_kg_s2, 0.0)
        self.assertAlmostEqual(value.fuel_line.derivative.manifold_pressure_pa_s, 0.0)
        self.assertAlmostEqual(value.stored_liquid_mass_rate_kg_s, 0.0)
        self.assertAlmostEqual(value.component_to_total_flow_closure_error_kg_s, 0.0)
        self.assertAlmostEqual(value.vehicle_mass_closure_error_kg_s, 0.0)
        self.assertEqual(value.depletion.limiting_component, "oxidizer")
        self.assertAlmostEqual(value.propulsion.acceleration_m_s2[0], 1.365)

    def test_mixture_ratio_events_are_positive_safe_and_inactive_safe(self) -> None:
        value = self.evaluate()
        surfaces = self.coupling.mixture_ratio_surfaces()
        self.assertEqual(
            tuple(surface.kind for surface in surfaces),
            (
                BipropellantMixtureRatioEventKind.MINIMUM,
                BipropellantMixtureRatioEventKind.MAXIMUM,
            ),
        )
        self.assertTrue(
            all(surface.value(value.engine_mixture_ratio) > 0.0 for surface in surfaces)
        )
        self.assertEqual(surfaces[0].value(2.5), 0.0)
        self.assertEqual(surfaces[1].value(3.0), 0.0)
        inactive = self.evaluate(relative_time_s=3.0)
        self.assertFalse(inactive.propulsion.active)
        self.assertIsNone(inactive.engine_mixture_ratio)
        self.assertTrue(all(surface.value(None) > 0.0 for surface in surfaces))

    def test_returned_tank_line_ledger_closes_across_reverse_event_trials(self) -> None:
        # Derive storage from the propagated pressure RHS, independently of
        # both the reported storage and the reported vehicle closure.
        flows = ((self.oxidizer_engine_flow, self.fuel_engine_flow), (0.0, 0.0),
                 (-0.1, 0.0), (0.0, -0.2), (-0.1, -0.2), (-0.1, 0.3),
                 (0.2, -0.1), (-1e-100, self.fuel_engine_flow))
        for time in (0.5, 3.0):
            for ox, fuel in flows:
                with self.subTest(time=time, oxidizer=ox, fuel=fuel):
                    state = replace(
                        self.state,
                        oxidizer_line=replace(self.state.oxidizer_line, mass_flow_kg_s=ox),
                        fuel_line=replace(self.state.fuel_line, mass_flow_kg_s=fuel),
                    )
                    value = self.evaluate(state, relative_time_s=time)
                    for flow, feed, line, model, engine in (
                        (ox, value.oxidizer_feed, value.oxidizer_line, self.oxidizer_line,
                         value.oxidizer_engine_flow_kg_s),
                        (fuel, value.fuel_feed, value.fuel_line, self.fuel_line,
                         value.fuel_engine_flow_kg_s),
                    ):
                        self.assertEqual(feed.derivative.propellant_mass_kg_s, -flow)
                        storage = (model.configuration.compliance_kg_pa
                                   * line.derivative.manifold_pressure_pa_s)
                        self.assertAlmostEqual(
                            fsum((feed.derivative.propellant_mass_kg_s, storage, engine)),
                            0.0, delta=1e-12)
                        self.assertEqual(feed.reverse_flow_event_extension, flow < 0.0)
                    actual = abs(fsum((
                        value.oxidizer_feed.derivative.propellant_mass_kg_s,
                        value.fuel_feed.derivative.propellant_mass_kg_s,
                        value.stored_liquid_mass_rate_kg_s,
                        -value.propulsion.mass_derivative_kg_s,
                    )))
                    self.assertEqual(value.vehicle_mass_closure_error_kg_s, actual)
                    self.assertLessEqual(actual, 1e-12)
                    self.assertEqual(value.reverse_flow_event_extension, ox < 0.0 or fuel < 0.0)

    def test_reverse_cutoff_root_remains_locatable_with_conservative_trial_rhs(self) -> None:
        from scipy.integrate import solve_ivp

        initial = replace(
            self.state,
            oxidizer_line=replace(self.state.oxidizer_line, mass_flow_kg_s=0.01,
                                  manifold_pressure_pa=2.1e6),
        )
        y0 = (*astuple(initial.oxidizer_feed), *astuple(initial.oxidizer_line),
              *astuple(initial.fuel_feed), *astuple(initial.fuel_line))
        reverse_trials = []

        def rhs(time, vector):
            state = BipropellantDynamicFeedState(
                RegulatedFeedState(*vector[:6]), DynamicFeedLineState(*vector[6:8]),
                RegulatedFeedState(*vector[8:14]), DynamicFeedLineState(*vector[14:16]))
            result = self.evaluate(state, relative_time_s=3.0 + time)
            reverse_trials.append(result.reverse_flow_event_extension)
            return (*result.oxidizer_feed.derivative.as_tuple(),
                    *result.oxidizer_line.derivative.as_tuple(),
                    *result.fuel_feed.derivative.as_tuple(),
                    *result.fuel_line.derivative.as_tuple())

        def cutoff(_time, vector):
            return vector[6]

        cutoff.terminal = True
        cutoff.direction = -1
        result = solve_ivp(rhs, (0.0, 0.3), y0, events=cutoff,
                           max_step=0.025, rtol=1e-8, atol=1e-10)
        self.assertTrue(result.success, result.message)
        self.assertEqual(result.status, 1)
        self.assertEqual(len(result.t_events[0]), 1)
        self.assertGreater(result.t_events[0][0], 0.0)
        self.assertLess(result.t_events[0][0], 0.3)
        self.assertAlmostEqual(result.y_events[0][0][6], 0.0, delta=1e-14)
        self.assertTrue(any(reverse_trials), "The event must bracket negative trial states.")
        # With the engine off, tank liquid + compliant line inventory is
        # constant through every accepted step, including the zero-flow root.
        compliance = self.oxidizer_line.configuration.compliance_kg_pa
        expected = y0[0] + compliance * y0[7]
        for column in result.y.T:
            self.assertAlmostEqual(column[0] + compliance * column[7], expected, delta=2e-12)

    def test_four_state_and_mixture_ratio_partials_match_finite_differences(self) -> None:
        nominal = self.evaluate()
        jacobian = nominal.partials.line_state_jacobian
        acceleration_partials = dict(nominal.partials.acceleration_wrt_line_state)
        perturbations = (
            ("oxidizer_line", "mass_flow_kg_s", 1.0e-5),
            ("oxidizer_line", "manifold_pressure_pa", 1.0),
            ("fuel_line", "mass_flow_kg_s", 1.0e-5),
            ("fuel_line", "manifold_pressure_pa", 1.0),
        )
        for column, (leg, field, step) in enumerate(perturbations):
            original_line = getattr(self.state, leg)
            plus_state = replace(
                self.state,
                **{leg: replace(original_line, **{field: getattr(original_line, field) + step})},
            )
            minus_state = replace(
                self.state,
                **{leg: replace(original_line, **{field: getattr(original_line, field) - step})},
            )
            plus = self.evaluate(plus_state)
            minus = self.evaluate(minus_state)
            plus_rhs = (
                *plus.oxidizer_line.derivative.as_tuple(),
                *plus.fuel_line.derivative.as_tuple(),
            )
            minus_rhs = (
                *minus.oxidizer_line.derivative.as_tuple(),
                *minus.fuel_line.derivative.as_tuple(),
            )
            finite_rhs = tuple(
                (right - left) / (2.0 * step)
                for right, left in zip(plus_rhs, minus_rhs, strict=True)
            )
            for row in range(4):
                self.assertAlmostEqual(
                    jacobian[row][column],
                    finite_rhs[row],
                    delta=max(1e-8, abs(finite_rhs[row]) * 2e-8),
                )
            finite_acceleration = tuple(
                (right - left) / (2.0 * step)
                for right, left in zip(
                    plus.propulsion.acceleration_m_s2,
                    minus.propulsion.acceleration_m_s2,
                    strict=True,
                )
            )
            state_name = nominal.partials.line_state_order[column]
            for analytic, numeric in zip(
                acceleration_partials[state_name], finite_acceleration, strict=True
            ):
                self.assertAlmostEqual(analytic, numeric, delta=1e-11)

        mixture_partials = dict(nominal.partials.engine_mixture_ratio_wrt_conditions)
        for leg, condition_name in (
            ("oxidizer_line", "oxidizer_pressure_pa"),
            ("fuel_line", "fuel_pressure_pa"),
        ):
            original_line = getattr(self.state, leg)
            step = 1.0
            plus = self.evaluate(
                replace(
                    self.state,
                    **{
                        leg: replace(
                            original_line,
                            manifold_pressure_pa=original_line.manifold_pressure_pa + step,
                        )
                    },
                )
            )
            minus = self.evaluate(
                replace(
                    self.state,
                    **{
                        leg: replace(
                            original_line,
                            manifold_pressure_pa=original_line.manifold_pressure_pa - step,
                        )
                    },
                )
            )
            assert plus.engine_mixture_ratio is not None
            assert minus.engine_mixture_ratio is not None
            finite = (plus.engine_mixture_ratio - minus.engine_mixture_ratio) / (2.0 * step)
            self.assertAlmostEqual(mixture_partials[condition_name], finite, places=11)

    def test_feed_state_partials_and_component_cutoff(self) -> None:
        nominal = self.evaluate()
        partials = dict(nominal.partials.line_derivative_wrt_feed_state)
        for leg, feed_field, row in (
            ("oxidizer_feed", "tank_pressurant_mass_kg", 0),
            ("fuel_feed", "tank_pressurant_mass_kg", 2),
        ):
            original_feed = getattr(self.state, leg)
            step = 1.0e-8
            plus = self.evaluate(
                replace(
                    self.state,
                    **{
                        leg: replace(
                            original_feed,
                            **{feed_field: getattr(original_feed, feed_field) + step},
                        )
                    },
                )
            )
            minus = self.evaluate(
                replace(
                    self.state,
                    **{
                        leg: replace(
                            original_feed,
                            **{feed_field: getattr(original_feed, feed_field) - step},
                        )
                    },
                )
            )
            plus_rhs = (
                *plus.oxidizer_line.derivative.as_tuple(),
                *plus.fuel_line.derivative.as_tuple(),
            )
            minus_rhs = (
                *minus.oxidizer_line.derivative.as_tuple(),
                *minus.fuel_line.derivative.as_tuple(),
            )
            finite = (plus_rhs[row] - minus_rhs[row]) / (2.0 * step)
            key = f"{'oxidizer' if row == 0 else 'fuel'}.{feed_field}"
            self.assertAlmostEqual(partials[key][row], finite, delta=max(1e-7, abs(finite) * 2e-8))

        cutoff = self.evaluate(
            replace(
                self.state,
                oxidizer_feed=replace(self.state.oxidizer_feed, propellant_mass_kg=1.0),
            )
        )
        self.assertFalse(cutoff.propulsion.active)
        self.assertIn("oxidizer_mass_kg", cutoff.propulsion.inhibited_by)
        self.assertEqual(cutoff.depletion.oxidizer_time_to_reserve_s, 0.0)
        self.assertEqual(cutoff.depletion.limiting_component, "oxidizer")

    def test_manifest_and_invalid_architectures_fail_closed(self) -> None:
        manifest = self.coupling.manifest()
        self.assertEqual(manifest["model_hash"], self.coupling.model_hash)
        self.assertIn("NASA-20050206420", manifest["reference_ids"])
        with self.assertRaisesRegex(DomainError, "Static feed resistance"):
            BipropellantDynamicFeedCoupling(
                self.configuration,
                self.bridge,
                example_feed_system(),
                self.oxidizer_line,
                self.fuel_feed,
                self.fuel_line,
            )
        with self.assertRaisesRegex(DomainError, "stream-resolved"):
            BipropellantDynamicFeedCoupling(
                self.configuration,
                PropagatorPropulsionBridge(
                    example_artifact(),
                    stream_bindings=self.bridge.stream_bindings,
                ),
                self.oxidizer_feed,
                self.oxidizer_line,
                self.fuel_feed,
                self.fuel_line,
            )

    def test_configuration_names_and_numeric_limits_have_typed_failures(self) -> None:
        for field in ("oxidizer_state_name", "fuel_state_name",
                      "oxidizer_pressure_condition_name", "fuel_pressure_condition_name", "model_id"):
            for value in (None, True, 0, [], {}, "", " "):
                with self.subTest(field=field, value=value), self.assertRaises(DomainError):
                    replace(self.configuration, **{field: value})
        for field in ("minimum_engine_mixture_ratio", "maximum_engine_mixture_ratio",
                      "mass_closure_absolute_tolerance_kg_s"):
            for value in (None, "2.0", [], {}, object(), float("nan"), float("inf"),
                          -float("inf"), 10**500):
                with self.subTest(field=field, value=value), self.assertRaises(DomainError):
                    replace(self.configuration, **{field: value})

    def test_boolean_limits_are_rejected_and_numeric_identity_is_canonical(self) -> None:
        base = replace(self.configuration, minimum_engine_mixture_ratio=0.5,
                       maximum_engine_mixture_ratio=2.0)
        for field in ("minimum_engine_mixture_ratio", "maximum_engine_mixture_ratio",
                      "mass_closure_absolute_tolerance_kg_s"):
            with self.subTest(field=field), self.assertRaises(DomainError):
                replace(base, **{field: True})
        integer_config = replace(self.configuration, minimum_engine_mixture_ratio=1,
                                 maximum_engine_mixture_ratio=3)
        float_config = replace(self.configuration, minimum_engine_mixture_ratio=1.0,
                               maximum_engine_mixture_ratio=3.0)
        self.assertEqual(replace(self.coupling, configuration=integer_config).model_hash,
                         replace(self.coupling, configuration=float_config).model_hash)

    def test_incomplete_state_and_coupling_records_are_rejected_at_construction(self) -> None:
        for field in ("oxidizer_feed", "oxidizer_line", "fuel_feed", "fuel_line"):
            for value in (None, True, object(), self.configuration):
                with self.subTest(state_field=field), self.assertRaises(DomainError):
                    replace(self.state, **{field: value})
        for field in ("configuration", "bridge", "oxidizer_feed", "oxidizer_line",
                      "fuel_feed", "fuel_line"):
            for value in (None, True, object(), self.state):
                with self.subTest(coupling_field=field), self.assertRaises(DomainError):
                    replace(self.coupling, **{field: value})

    def test_event_boundary_kind_and_flags_cannot_silently_change_cutoff_semantics(self) -> None:
        minimum = BipropellantMixtureRatioEventSurface(
            "minimum", "minimum_engine_mixture_ratio", 2.5)
        with self.subTest(case="registered string minimum"):
            self.assertEqual(minimum.value(3.0), 0.5)
        for kind in (None, True, 1, "", "typo", [], {}):
            with self.subTest(kind=kind), self.assertRaises(DomainError):
                BipropellantMixtureRatioEventSurface("event", kind, 2.5)
        for flag in (None, 1, "false", []):
            with self.subTest(flag=flag), self.assertRaises(DomainError):
                replace(self.coupling.mixture_ratio_surfaces()[0], stops_burn=flag)
        for boundary in (None, True, "2.5", [], float("nan"), float("inf"), 0.0, -1.0):
            with self.subTest(boundary=boundary), self.assertRaises(DomainError):
                replace(self.coupling.mixture_ratio_surfaces()[0], boundary_value=boundary)
        for ratio in (True, False, "3.0", [], float("nan"), float("inf"), 0.0, -1.0):
            with self.subTest(ratio=ratio), self.assertRaises(DomainError):
                minimum.value(ratio)
        self.assertEqual(minimum.value(None), 1.0)


if __name__ == "__main__":
    unittest.main()
