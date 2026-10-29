"""Independent direct-query typing and complete polytropic closure oracles."""
import math
import unittest
from decimal import Decimal, localcontext
from fractions import Fraction

from rocket_propulsion.core.errors import DomainError
from rocket_propulsion.propulsion.burns import (
    BlowdownTankModel,
    ConstantPerformanceProvider,
    PressurePerformanceLaw,
    PressurePerformanceSample,
    PressurePerformanceTable,
)


def pressure_oracle(tank, remaining):
    with localcontext() as ctx:
        ctx.prec = 100
        p0, v0, m0, rho, power, mass = map(Decimal, (
            tank.initial_pressure_pa, tank.initial_ullage_volume_m3,
            tank.initial_propellant_mass_kg, tank.propellant_density_kg_m3,
            tank.polytropic_exponent, remaining))
        volume = v0+(m0-mass)/rho
        return float(p0*(v0/volume)**power)


class BlowdownDirectQueryTests(unittest.TestCase):
    def test_query_methods_refuse_malformed_and_erased_input_with_field(self):
        tank = BlowdownTankModel(2e6, .01, 10., 1000.)
        for method, field in ((tank.pressure_pa, 'remaining_propellant_mass_kg'),
                              (tank.ullage_volume_m3, 'remaining_propellant_mass_kg'),
                              (tank.remaining_mass_at_pressure, 'pressure_pa')):
            for value in (True, False, '5', None, [], complex(1), float('nan'),
                          float('inf'), Fraction(1, 10**500), Fraction(10**400)):
                with self.subTest(method=method.__name__, value=repr(value)[:40]):
                    with self.assertRaises(DomainError) as error:
                        method(value)
                    self.assertEqual(error.exception.field, field)

    def test_public_volume_outputs_do_not_export_false_zero_or_infinity(self):
        for tank, property_name in (
            (BlowdownTankModel(1., 1., 1e-300, 1e300), 'initial_propellant_volume_m3'),
            (BlowdownTankModel(1., 1e308, 1e308, 1.), 'tank_volume_m3'),
        ):
            with self.subTest(property=property_name), self.assertRaises(DomainError):
                getattr(tank, property_name)
        with self.assertRaises(DomainError):
            BlowdownTankModel(1., 1e308, 1e308, 1.).ullage_volume_m3(0.)

    def test_complete_pressure_survives_unrepresentable_intermediate_volume(self):
        for tank in (BlowdownTankModel(1e300, 1e200, 1e200, 1e-200),
                     BlowdownTankModel(1e300, 1e300, 1e300, 1e-300)):
            with self.subTest(tank=tank):
                expected = pressure_oracle(tank, 0.)
                try:
                    actual = tank.pressure_pa(0.)
                except DomainError as error:
                    self.fail(f'A finite complete pressure was refused: {error}')
                self.assertTrue(math.isclose(actual, expected, rel_tol=3e-13))
                self.assertGreater(actual, 0.)

    def test_small_volume_growth_is_not_erased_before_large_exponent(self):
        for mass, exponent in ((1e-16, 1e16), (1e-100, 1e100), (1e-300, 1e300)):
            tank = BlowdownTankModel(2e6, 1., mass, 1., exponent)
            # Independent log(1+x) series; direct (1+x) at even 100-digit
            # precision would erase the last case before exponentiation.
            with localcontext() as ctx:
                ctx.prec = 100
                x = Decimal(mass)
                logarithm = x-x*x/2+x*x*x/3
                expected = float(Decimal(2_000_000)*(-Decimal(exponent)*logarithm).exp())
            with self.subTest(mass=mass, exponent=exponent):
                self.assertTrue(math.isclose(tank.pressure_pa(0.), expected, rel_tol=3e-13))

    def test_genuinely_unrepresentable_pressure_still_refuses(self):
        tank = BlowdownTankModel(1e-300, 1e-200, 1e200, 1e-200)
        with self.assertRaises(DomainError):
            tank.pressure_pa(0.)

    def test_small_growth_fallback_handles_capacity_overflow_and_subnormal_ratio(self):
        cases = (BlowdownTankModel(2e6, 1e200, 1e100, 1e200, 1e300),
                 BlowdownTankModel(2e6, 1e300, 1e-24, 1., 1e308))
        for tank in cases:
            with localcontext() as ctx:
                ctx.prec = 100
                x = (Decimal(tank.initial_propellant_mass_kg)
                     /(Decimal(tank.propellant_density_kg_m3)*Decimal(tank.initial_ullage_volume_m3)))
                log_growth = x-x*x/2+x*x*x/3
                expected = float(Decimal(tank.initial_pressure_pa)
                                 *(-Decimal(tank.polytropic_exponent)*log_growth).exp())
            with self.subTest(tank=tank):
                self.assertEqual(tank.pressure_pa(0.), expected)
                self.assertLess(expected, tank.initial_pressure_pa)

    def test_nonpositive_cutoff_pressure_names_its_field(self):
        tank = BlowdownTankModel(2e6, .01, 10., 1000.)
        for value in (0., -1.):
            with self.subTest(value=value):
                with self.assertRaises(DomainError) as error:
                    tank.remaining_mass_at_pressure(value)
                self.assertEqual(error.exception.field, 'pressure_pa')

    def test_reference_and_fractional_queries_preserve_the_physical_closure(self):
        tank = BlowdownTankModel(2e6, .01, 10., 1000.)
        self.assertEqual(tank.pressure_pa(Fraction(10)), 2e6)
        expected = pressure_oracle(tank, 5.)
        self.assertTrue(math.isclose(tank.pressure_pa(Fraction(5)), expected, rel_tol=3e-14))
        self.assertTrue(math.isclose(tank.remaining_mass_at_pressure(expected), 5., abs_tol=2e-14))
        self.assertEqual(tank.remaining_mass_at_pressure(Fraction(1_000_000)), 0.)
        self.assertEqual(tank.pressure_pa(-0.), tank.pressure_pa(0.))

    def test_performance_queries_also_refuse_boolean_or_text_pressure(self):
        point = ConstantPerformanceProvider.monopropellant(
            delivered_thrust_n=100., system_specific_impulse_s=220., chamber_pressure_pa=1.5,
            source='synthetic low-pressure numeric boundary').operating_point()
        models = (PressurePerformanceLaw(2., .5), PressurePerformanceTable((
            PressurePerformanceSample(.5, 25., 220., .375),
            PressurePerformanceSample(2., 100., 220., 1.5)), 'synthetic map'))
        for model in models:
            for value in (True, False, '1', None, [], Fraction(1, 10**500)):
                with self.subTest(model=type(model).__name__, value=repr(value)[:35]):
                    with self.assertRaises(DomainError) as error:
                        model.scale_operating_point(point, value)
                    self.assertEqual(error.exception.field, 'supply_pressure_pa')
