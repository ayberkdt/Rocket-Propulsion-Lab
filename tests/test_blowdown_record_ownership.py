"""Input ownership controls for pressure-fed records and cached correlations.

These synthetic cases exercise mutable numeric subclasses, caller-owned maps,
typed refusal, and the F/(g0*Isp) identity. They are not hardware validation.
"""
from __future__ import annotations

import math
import unittest
from dataclasses import fields
from fractions import Fraction

from rocket_propulsion.core.errors import DomainError
from rocket_propulsion.propulsion.burns import (
    BlowdownTankModel,
    ConstantPerformanceProvider,
    PressurePerformanceLaw,
    PressurePerformanceSample,
    PressurePerformanceTable,
)


class MutableExponent(float):
    """Accepted Real whose arithmetic can change without changing its value."""

    def __new__(cls, value):
        result = super().__new__(cls, value)
        result.factor = 1.
        return result

    def __rpow__(self, base):
        return float(base) ** (float(self) * self.factor)


TANK = {'initial_pressure_pa':2_000_000, 'initial_ullage_volume_m3':.01,
        'initial_propellant_mass_kg':10, 'propellant_density_kg_m3':1000,
        'polytropic_exponent':1}
LAW = {'reference_supply_pressure_pa':2_000_000, 'minimum_supply_pressure_pa':1_000_000,
       'thrust_pressure_exponent':1, 'specific_impulse_pressure_exponent':0,
       'chamber_pressure_exponent':1, 'maximum_supply_pressure_pa':3_000_000}
ROW = {'supply_pressure_pa':1_000_000, 'delivered_thrust_n':50,
       'system_specific_impulse_s':220, 'chamber_pressure_pa':750_000}
CASES = ((BlowdownTankModel, TANK), (PressurePerformanceLaw, LAW),
         (PressurePerformanceSample, ROW))


def rated_point():
    return ConstantPerformanceProvider.monopropellant(
        delivered_thrust_n=100., system_specific_impulse_s=220.,
        chamber_pressure_pa=1.5e6, source='pressure-record ownership fixture',
    ).operating_point()


def rows():
    return [PressurePerformanceSample(**ROW), PressurePerformanceSample(2e6, 100., 220., 1.5e6)]


class PressureRecordOwnershipTests(unittest.TestCase):
    def test_caller_numeric_mutation_cannot_change_law_or_flow(self):
        exponent = MutableExponent(1.)
        law = PressurePerformanceLaw(2e6, 1e6, thrust_pressure_exponent=exponent)
        point = rated_point()
        before = law.scale_operating_point(point, 1e6)
        exponent.factor = 2.
        after = law.scale_operating_point(point, 1e6)
        self.assertEqual(before.delivered_thrust_n, 50.)
        self.assertEqual(after.delivered_thrust_n, 50., 'Caller-owned exponent changed a frozen law.')
        self.assertTrue(math.isclose(after.total_tank_flow_kg_s, 50/(9.80665*220), rel_tol=2e-14))
        self.assertIs(type(law.thrust_pressure_exponent), float)

    def test_list_mutation_cannot_change_table_range_or_interpolation(self):
        caller = rows()
        table = PressurePerformanceTable(caller, 'measured-map-synthetic-fixture')
        point = rated_point()
        self.assertEqual(table.scale_operating_point(point, 1.5e6).delivered_thrust_n, 75.)
        caller[0] = PressurePerformanceSample(1e6, 60., 220., 7.5e5)
        self.assertEqual(table.scale_operating_point(point, 1.5e6).delivered_thrust_n, 75.)
        caller.clear()
        self.assertEqual(table.minimum_supply_pressure_pa, 1e6)
        self.assertEqual(table.maximum_pressure_pa, 2e6)
        self.assertIs(type(table.samples), tuple)
        self.assertEqual(table.scale_operating_point(point, 1.5e6).delivered_thrust_n, 75.)

    def test_supported_real_values_are_owned_native_floats(self):
        for kind, defaults in CASES:
            for factory in (int, Fraction, MutableExponent):
                with self.subTest(record=kind.__name__, numeric=factory.__name__):
                    # Keep the non-integral ullage value when testing integer inputs.
                    values = {name: factory(value) if value == int(value) else Fraction(1, 100)
                              for name, value in defaults.items()}
                    record = kind(**values)
                    for item in fields(record):
                        self.assertIs(type(getattr(record, item.name)), float, item.name)
                    if kind is BlowdownTankModel:
                        self.assertTrue(math.isclose(record.pressure_pa(5.),
                                                     2e6*.01/.015, rel_tol=2e-14))
        self.assertIsNone(PressurePerformanceLaw(2e6, 1e6).maximum_supply_pressure_pa)

    def test_every_numeric_field_refuses_malformed_or_erased_input(self):
        invalid = (True, False, '1', [], {}, complex(1), float('nan'), float('inf'),
                   -float('inf'), Fraction(10**400), Fraction(1, 10**500))
        for kind, defaults in CASES:
            for name in defaults:
                for value in invalid:
                    with self.subTest(record=kind.__name__, field=name, bad=repr(value)[:55]):
                        with self.assertRaises(DomainError) as error:
                            kind(**(defaults | {name: value}))
                        self.assertIn(name, error.exception.field or '')

    def test_table_refuses_invalid_row_containers_with_typed_errors(self):
        for value in (None, 3, 'rows', {1: rows()[0]}, set(rows()), [object(), object()]):
            with self.subTest(value=repr(value)[:55]):
                with self.assertRaises(DomainError) as error:
                    PressurePerformanceTable(value, 'source')
                self.assertEqual(error.exception.field, 'samples')

    def test_table_source_is_validated_without_calling_foreign_strip(self):
        for source in (None, 1, True, [], '', '   ', '\ud800'):
            with self.subTest(source=repr(source)):
                with self.assertRaises(DomainError) as error:
                    PressurePerformanceTable(tuple(rows()), source)
                self.assertIn('source_id', error.exception.field or '')

    def test_ordered_iterable_rows_are_captured_once(self):
        table = PressurePerformanceTable(iter(rows()), 'single-pass-map')
        self.assertEqual(len(table.samples), 2)
        self.assertEqual(table.scale_operating_point(rated_point(), 1.5e6).delivered_thrust_n, 75.)


if __name__ == '__main__':
    unittest.main()
