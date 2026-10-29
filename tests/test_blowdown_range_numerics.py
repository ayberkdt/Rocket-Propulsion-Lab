"""Range attacks with independent high-precision pressure/stream oracles.

These synthetic cases exercise the declared polytropic and pressure-power
equations (NASA SP-8112 / CR-122347), not a calibrated hardware operating map.
"""

from __future__ import annotations

import math
import unittest
from decimal import Decimal, localcontext

from rocket_propulsion.core.errors import DomainError
from rocket_propulsion.propulsion.burns import (
    BlowdownTankModel,
    ConstantPerformanceProvider,
    PressurePerformanceLaw,
)


def decimal_power(value, pressure, reference, exponent):
    """Independent 90-digit log/exp oracle on the exact supplied floats."""
    with localcontext() as context:
        context.prec = 90
        v, p, r, n = map(Decimal, (value, pressure, reference, exponent))
        return float(v * ((p / r).ln() * n).exp())


def point(thrust=100., chamber_pressure=1.5e6):
    return ConstantPerformanceProvider.monopropellant(
        delivered_thrust_n=thrust, system_specific_impulse_s=220.,
        chamber_pressure_pa=chamber_pressure, source="synthetic range oracle",
    ).operating_point()


class BlowdownPressureRangeTests(unittest.TestCase):
    def test_tiny_positive_polytropic_exponent_clips_unreachable_cutoff_to_empty(self):
        tank = BlowdownTankModel(2e6, .01, 10., 1000., 1e-300)
        self.assertEqual(tank.remaining_mass_at_pressure(1e6), 0.)

    def test_inverse_preserves_finite_answer_when_unscaled_power_overflows(self):
        tank = BlowdownTankModel(1e300, 1e-300, 1e21, 1., .5)
        with localcontext() as context:
            context.prec = 90
            p0, v0, mass, rho, n, cutoff = map(Decimal,
                (tank.initial_pressure_pa, tank.initial_ullage_volume_m3,
                 tank.initial_propellant_mass_kg, tank.propellant_density_kg_m3,
                 tank.polytropic_exponent, 1e140))
            expected = mass - rho * (v0 * (p0 / cutoff) ** (1 / n) - v0)
        actual = tank.remaining_mass_at_pressure(1e140)
        self.assertTrue(math.isclose(actual, float(expected), rel_tol=3e-15), (actual, expected))
        self.assertGreater(actual, 0.)
        self.assertLess(actual, tank.initial_propellant_mass_kg)

    def test_inverse_resolves_near_initial_pressure_without_losing_expelled_volume(self):
        tank = BlowdownTankModel(2e6, .01, 1e-12, 1000., 1.2)
        cutoff = math.nextafter(tank.initial_pressure_pa, 0.)
        with localcontext() as context:
            context.prec = 90
            p0, v0, mass, rho, n, p = map(Decimal,
                (2e6, .01, 1e-12, 1000., 1.2, cutoff))
            expected = mass - rho * v0 * ((p0 / p) ** (1 / n) - 1)
        actual = tank.remaining_mass_at_pressure(cutoff)
        self.assertTrue(math.isclose(actual, float(expected), rel_tol=3e-15), (actual, expected))

    def test_forward_pressure_recovers_underflowed_power_after_scaling(self):
        tank = BlowdownTankModel(1e300, .01, 10., 1000., 1100.)
        expected = decimal_power(1e300, .01, .02, 1100.)
        actual = tank.pressure_pa(0.)
        self.assertGreater(actual, 0.)
        self.assertTrue(math.isclose(actual, expected, rel_tol=5e-14), (actual, expected))

    def test_unrepresentable_positive_pressure_is_refused_instead_of_returning_zero(self):
        tank = BlowdownTankModel(2e6, .01, 10., 1000., 1e300)
        with self.assertRaises(DomainError):
            tank.pressure_pa(0.)

    def test_cutoff_clipping_and_reference_pressure_remain_exact(self):
        tank = BlowdownTankModel(2e6, .01, 10., 1000., 1.2)
        self.assertEqual(tank.remaining_mass_at_pressure(2e6), 10.)
        self.assertEqual(tank.remaining_mass_at_pressure(3e6), 10.)
        self.assertEqual(tank.remaining_mass_at_pressure(1.), 0.)
        self.assertEqual(tank.pressure_pa(10.), 2e6)


class PressurePerformanceRangeTests(unittest.TestCase):
    def test_extreme_exponents_are_checked_at_evaluation_not_arbitrarily_bounded(self):
        for field in ('thrust_pressure_exponent', 'specific_impulse_pressure_exponent',
                      'chamber_pressure_exponent'):
            for exponent in (1e9, 1e300):
                law = PressurePerformanceLaw(2e6, 1e6, **{field: exponent})
                rated = point()
                with self.subTest(field=field, exponent=exponent):
                    at_reference = law.scale_operating_point(rated, 2e6)
                    self.assertEqual(at_reference.delivered_thrust_n, rated.delivered_thrust_n)
                    with self.assertRaises(DomainError):
                        law.scale_operating_point(rated, 1e6)

    def test_negative_isp_exponent_overflow_is_a_domain_error(self):
        law = PressurePerformanceLaw(2e6, 1e6, specific_impulse_pressure_exponent=-1e300)
        with self.assertRaises(DomainError):
            law.scale_operating_point(point(), 1e6)

    def test_overflowing_power_does_not_refuse_representable_scaled_thrust_or_streams(self):
        rated = point(1e-300, 1e-300)
        law = PressurePerformanceLaw(2e6, 1e6, thrust_pressure_exponent=1100.,
            chamber_pressure_exponent=1100., maximum_supply_pressure_pa=4e6)
        expected = decimal_power(1e-300, 4e6, 2e6, 1100.)
        self.assertTrue(math.isfinite(expected) and expected > 0.)
        try:
            result = law.scale_operating_point(rated, 4e6)
        except DomainError as error:
            self.fail(f'Finite scaled thrust {expected!r} was falsely refused: {error}')
        self.assertTrue(math.isclose(result.delivered_thrust_n, expected, rel_tol=5e-14))
        self.assertTrue(math.isclose(result.chamber_pressure_pa, expected, rel_tol=5e-14))
        self.assertGreater(result.streams[0].mass_flow_kg_s, 0.)
        self.assertEqual(result.total_tank_flow_kg_s, result.streams[0].mass_flow_kg_s)

    def test_underflowing_power_does_not_erase_representable_scaled_thrust_or_streams(self):
        rated = point(1e300, 1e300)
        law = PressurePerformanceLaw(2e6, 1e6, thrust_pressure_exponent=1100.,
            chamber_pressure_exponent=1100.)
        expected = decimal_power(1e300, 1e6, 2e6, 1100.)
        self.assertTrue(math.isfinite(expected) and expected > 0.)
        try:
            result = law.scale_operating_point(rated, 1e6)
        except DomainError as error:
            self.fail(f'Finite scaled thrust {expected!r} was falsely refused: {error}')
        self.assertTrue(math.isclose(result.delivered_thrust_n, expected, rel_tol=5e-14))
        self.assertTrue(math.isclose(result.chamber_pressure_pa, expected, rel_tol=5e-14))
        self.assertGreater(result.streams[0].mass_flow_kg_s, 0.)
        self.assertEqual(result.total_tank_flow_kg_s, result.streams[0].mass_flow_kg_s)

    def test_finite_flow_survives_overflowing_gravity_times_isp_intermediate(self):
        law = PressurePerformanceLaw(2e6, 1e6, specific_impulse_pressure_exponent=1015.,
            maximum_supply_pressure_pa=4e6)
        result = law.scale_operating_point(point(), 4e6)
        with localcontext() as context:
            context.prec = 90
            expected = Decimal(result.delivered_thrust_n) / (
                Decimal('9.80665') * Decimal(result.system_specific_impulse_s))
        self.assertGreater(result.total_tank_flow_kg_s, 0.)
        self.assertTrue(math.isclose(result.total_tank_flow_kg_s, float(expected), rel_tol=5e-14))

    def test_near_reference_ratio_rounding_is_not_amplified_by_large_exponent(self):
        law = PressurePerformanceLaw(2e6, 1e6, thrust_pressure_exponent=1e16)
        pressure = math.nextafter(2e6, 0.)
        result = law.scale_operating_point(point(), pressure)
        expected = decimal_power(100., pressure, 2e6, 1e16)
        self.assertTrue(math.isclose(result.delivered_thrust_n, expected, rel_tol=5e-14),
                        (result.delivered_thrust_n, expected))

    def test_explicit_zero_maximum_is_not_silently_replaced_by_reference(self):
        with self.assertRaises(DomainError):
            PressurePerformanceLaw(2e6, 1e6, maximum_supply_pressure_pa=0.)

    def test_zero_exponents_preserve_outputs_when_pressure_ratio_overflows(self):
        # The ratio is outside binary64, but r**0 is exactly one. This must
        # remain a valid constant correlation, including its complete streams.
        law = PressurePerformanceLaw(1e-300, 1e-310, thrust_pressure_exponent=0.,
            specific_impulse_pressure_exponent=0., chamber_pressure_exponent=0.,
            maximum_supply_pressure_pa=1e300)
        rated = point()
        result = law.scale_operating_point(rated, 1e300)
        for field in ('delivered_thrust_n', 'ideal_thrust_n', 'system_specific_impulse_s',
                      'chamber_pressure_pa', 'total_tank_flow_kg_s', 'chamber_mass_flow_kg_s'):
            self.assertTrue(math.isclose(getattr(result,field),getattr(rated,field),rel_tol=3e-15),field)
        self.assertEqual(len(result.streams),len(rated.streams))
        for actual,reference in zip(result.streams,rated.streams,strict=True):
            for field in ('stream_id','tank_id','propellant_key','role','destination',
                          'tank_depleting','contributes_to_delivered_thrust'):
                self.assertEqual(getattr(actual,field),getattr(reference,field),field)
            self.assertTrue(math.isclose(actual.mass_flow_kg_s,reference.mass_flow_kg_s,rel_tol=3e-15))


if __name__ == '__main__':
    unittest.main()
