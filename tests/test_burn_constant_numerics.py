"""Public L0 burn inversion, event bounds and independently evaluated ledgers."""

import json
import unittest
from dataclasses import replace
from decimal import Decimal, localcontext
from math import nextafter, ulp

from rocket_propulsion.core.errors import DomainError
from rocket_propulsion.propulsion import STANDARD_GRAVITY_M_S2
from rocket_propulsion.propulsion.burns import (
    BurnDefinition, BurnTarget, BurnTargetKind, CutoffReason, EngineOperatingPoint,
    FlowDestination, PropellantFlow, PropellantRole, TankDefinition,
    burn_result_from_json, burn_result_to_json, simulate_constant_burn,
)
from tests.test_burns import example_definition, example_operating_point


def mono(*, flow=1., force=10., mass=100., loaded=10., reserve=0., dry=80.,
         kind=BurnTargetKind.DURATION, target=1., maximum=100., cant=1.):
    role = PropellantRole.MONOPROPELLANT
    point = EngineOperatingPoint(
        "numerical synthetic point", 1., 1., 1, force, force, 100.,
        force/(STANDARD_GRAVITY_M_S2*flow), flow, flow, 0., 0., None, 1e6, 0.,
        (PropellantFlow("main", "tank", "synthetic", role, FlowDestination.CHAMBER, flow),))
    definition = BurnDefinition("numerical synthetic burn", mass, dry,
        (TankDefinition("tank", "synthetic", role, loaded, reserve),),
        BurnTarget(kind, target), maximum, cant)
    return definition, point


def inverse_series(definition, point, precision):
    """Invert via the integral series sum (-x)^n/(n+1)!, without exp subtraction."""
    with localcontext() as context:
        context.prec = precision
        target, mass, force, flow, cant = map(Decimal.from_float, (
            definition.target.value, definition.initial_mass_kg, point.delivered_thrust_n,
            point.total_tank_flow_kg_s, definition.cant_efficiency))
        x = target*flow/(force*cant)
        assert x < Decimal('.1')
        term = total = Decimal(1)
        for index in range(1, 1000):
            term *= -x/Decimal(index+1)
            following = total+term
            if following == total:
                break
            total = following
        else:
            raise AssertionError("reference series did not settle")
        return float(target*mass/(force*cant)*total)


class ConstantBurnNumericalAttacks(unittest.TestCase):
    def close_ulp(self, actual, expected, count=4):
        self.assertLessEqual(abs(actual-expected), count*ulp(expected), (actual, expected))

    def test_small_delta_v_inverse_and_forward_result_do_not_flatten(self):
        point = example_operating_point()
        for target in (1e-10, 1e-12, 1e-14, 1e-100, 1e-300):
            with self.subTest(target=target):
                definition = example_definition(BurnTarget(BurnTargetKind.IDEAL_DELTA_V, target))
                reference = inverse_series(definition, point, 300)
                self.assertEqual(reference, inverse_series(definition, point, 200))
                result = simulate_constant_burn(definition, point)
                self.close_ulp(result.summary.duration_s, reference)
                self.close_ulp(result.summary.ideal_delta_v_m_s, target)
                self.assertTrue(result.summary.target_achieved)

    def test_small_burn_tank_consumption_is_not_subtracted_from_rounded_states(self):
        definition, point = mono(target=1e-20)
        summary = simulate_constant_burn(definition, point).summary
        self.assertEqual(summary.initial_mass_kg, summary.final_mass_kg)
        self.assertEqual(dict(summary.consumed_by_tank_kg)['tank'], 1e-20)
        self.assertEqual(summary.consumed_propellant_kg, 1e-20)
        self.assertEqual(summary.mass_closure_error_kg, 0.)

    def test_neighboring_maximum_duration_cannot_be_overridden_by_target_tolerance(self):
        for maximum, reached in ((nextafter(1., 0.), False), (1., True), (nextafter(1., 2.), True)):
            with self.subTest(maximum=maximum):
                definition, point = mono(maximum=maximum)
                result = simulate_constant_burn(definition, point).summary
                self.assertEqual(result.duration_s, min(1., maximum))
                self.assertEqual(result.target_achieved, reached)
                self.assertEqual(result.cutoff_reason,
                    CutoffReason.TARGET_REACHED if reached else CutoffReason.MAXIMUM_DURATION)

    def test_neighboring_tank_and_dry_limits_keep_binding_constraints(self):
        for mode in ('tank', 'dry'):
            # The existing definition tolerance accepts this one-ulp dry/load
            # budget mismatch; a target tolerance must not then cross the floor.
            dry = nextafter(98., 100.) if mode=='dry' else 80.
            boundary = 100.-dry if mode=='dry' else 2.
            for target, reached in ((nextafter(boundary, 0.), True), (boundary, True),
                                     (nextafter(boundary, 3.), False)):
                with self.subTest(mode=mode, target=target):
                    definition, point = mono(target=target, loaded=3. if mode=='tank' else 2.,
                                             reserve=1. if mode=='tank' else 0., dry=dry)
                    result = simulate_constant_burn(definition, point).summary
                    self.assertEqual(result.duration_s, min(boundary, target))
                    self.assertEqual(result.target_achieved, reached)
                    self.assertGreaterEqual(result.final_mass_kg, definition.protected_dry_mass_kg)
                    if not reached:
                        self.assertEqual(result.binding_constraint, 'tank' if mode=='tank' else 'protected_dry_mass_kg')

    def test_zero_available_inventory_returns_a_finite_unachieved_cutoff(self):
        definition, point = mono(loaded=1., reserve=1.)
        result = simulate_constant_burn(definition, point)
        summary = result.summary
        self.assertFalse(summary.target_achieved)
        self.assertEqual(summary.cutoff_reason, CutoffReason.TANK_RESERVE)
        self.assertEqual(summary.duration_s, 0.)
        self.assertEqual(summary.consumed_propellant_kg, 0.)
        self.assertEqual(summary.ideal_delta_v_m_s, 0.)
        self.assertEqual(summary.rocket_equivalent_specific_impulse_s,
                         point.delivered_thrust_n/(STANDARD_GRAVITY_M_S2*point.total_tank_flow_kg_s))
        self.assertEqual(burn_result_from_json(burn_result_to_json(result)), result)

    def test_small_mass_fraction_recovers_finite_forward_velocity(self):
        definition, point = mono(flow=1e-200, force=1e100, mass=1e200, loaded=1., dry=0.)
        summary = simulate_constant_burn(definition, point).summary
        self.close_ulp(summary.ideal_delta_v_m_s, 1e-100)
        self.close_ulp(summary.consumed_propellant_kg, 1e-200)
        self.assertEqual(summary.mass_closure_error_kg, 0.)

    def test_unrepresentable_exhaust_velocity_does_not_destroy_finite_delta_v(self):
        definition, point = mono(flow=1e-308, force=5.)
        summary = simulate_constant_burn(definition, point).summary
        self.close_ulp(summary.ideal_delta_v_m_s, .05)
        self.close_ulp(summary.rocket_equivalent_specific_impulse_s, point.system_specific_impulse_s)

    def test_nearly_depleted_burn_keeps_positive_reserve_and_logarithm(self):
        definition, point = mono(mass=1e16, loaded=1e16, reserve=1., dry=0., target=1e16, maximum=1e16)
        result = simulate_constant_burn(definition, point)
        self.assertFalse(result.summary.target_achieved)
        self.assertEqual(result.summary.cutoff_reason, CutoffReason.TANK_RESERVE)
        self.assertEqual(result.summary.final_mass_kg, 1.)
        self.assertEqual(dict(result.final_state.tank_masses_kg)['tank'], 1.)
        with localcontext() as context:
            context.prec=100
            reference=float(Decimal(10)*Decimal('1e16').ln())
        self.close_ulp(result.summary.ideal_delta_v_m_s, reference)

    def test_unreachable_mass_and_impulse_targets_retain_finite_maximum_cutoff(self):
        for kind in (BurnTargetKind.PROPELLANT_MASS, BurnTargetKind.TOTAL_IMPULSE):
            with self.subTest(kind=kind):
                definition, point = mono(flow=1e-200, force=1e-200, kind=kind,
                                         target=1e200, maximum=1.)
                result = simulate_constant_burn(definition, point).summary
                self.assertEqual(result.duration_s, 1.)
                self.assertFalse(result.target_achieved)
                self.assertEqual(result.cutoff_reason, CutoffReason.MAXIMUM_DURATION)

    def test_true_positive_output_underflow_is_a_typed_field_error(self):
        definition, point = mono(flow=1e-200, force=1e-200, target=1e-200)
        with self.assertRaises(DomainError) as caught:
            simulate_constant_burn(definition, point)
        self.assertIn('Consumed propellant', caught.exception.field or '')

    def test_true_impulse_overflow_is_a_typed_field_error(self):
        definition, point = mono(flow=1e100, force=1e308, mass=1e308, loaded=1e308,
                                 dry=0., target=10., maximum=10.)
        with self.assertRaises(DomainError) as caught:
            simulate_constant_burn(definition, point)
        self.assertIn('Delivered impulse', caught.exception.field or '')

    def test_declared_stream_mismatch_remains_visible_in_independent_mass_closure(self):
        definition, point = mono(flow=1., target=.125)
        modified = replace(point, streams=(replace(point.streams[0], mass_flow_kg_s=1.+2.**-34),))
        summary = simulate_constant_burn(definition, modified).summary
        self.assertEqual(dict(summary.consumed_by_tank_kg)['tank'], (1.+2.**-34)*.125)
        self.assertEqual(summary.mass_closure_error_kg, 2.**-37)

    def test_scaled_canted_inverse_and_forward_paths_match_independent_series(self):
        for scale in (1e-150, 1., 1e150):
            for cant in (.125, 1.):
                with self.subTest(scale=scale, cant=cant):
                    definition, point = mono(flow=scale, force=10.*scale, mass=100.*scale,
                        loaded=10.*scale, dry=80.*scale, kind=BurnTargetKind.IDEAL_DELTA_V,
                        target=1e-12, cant=cant)
                    reference = inverse_series(definition, point, 300)
                    self.assertEqual(reference, inverse_series(definition, point, 200))
                    result = simulate_constant_burn(definition, point)
                    self.close_ulp(result.summary.duration_s, reference)
                    self.close_ulp(result.summary.ideal_delta_v_m_s, 1e-12)
                    self.assertEqual(burn_result_from_json(burn_result_to_json(result)), result)
                    json.loads(burn_result_to_json(result))

    def test_forward_velocity_matches_atanh_integral_series(self):
        for time in (.125, .75, 2.):
            for cant in (.125, 1.):
                with self.subTest(time=time, cant=cant):
                    definition, point = mono(target=time, cant=cant)
                    values = []
                    for precision in (200, 300):
                        with localcontext() as context:
                            context.prec = precision
                            used = Decimal.from_float(time)
                            z = used/(Decimal(200)-used)
                            term = total = z
                            for index in range(1, 200):
                                term *= z*z
                                following = total+term/Decimal(2*index+1)
                                if following == total:
                                    break
                                total = following
                            values.append(float(Decimal(20)*Decimal.from_float(cant)*total))
                    self.assertEqual(*values)
                    self.close_ulp(simulate_constant_burn(definition, point).summary.ideal_delta_v_m_s, values[0])

    def test_huge_unreachable_delta_v_clips_before_exponential_output_limits(self):
        definition, point = mono(kind=BurnTargetKind.IDEAL_DELTA_V, target=1e308, maximum=2.)
        summary = simulate_constant_burn(definition, point).summary
        self.assertEqual(summary.duration_s, 2.)
        self.assertFalse(summary.target_achieved)
        self.assertEqual(summary.cutoff_reason, CutoffReason.MAXIMUM_DURATION)

    def test_ambient_decimal_precision_does_not_change_burn_results(self):
        definition, point = mono(kind=BurnTargetKind.IDEAL_DELTA_V, target=1e-100)
        reference = simulate_constant_burn(definition, point)
        with localcontext() as context:
            context.prec = 6
            self.assertEqual(simulate_constant_burn(definition, point), reference)
