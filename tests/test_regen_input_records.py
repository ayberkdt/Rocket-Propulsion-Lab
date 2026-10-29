"""Coolant direction and finite owned geometry at library entry points."""
import unittest
from dataclasses import fields
from fractions import Fraction

from rocket_propulsion.core.errors import DomainError, InputError
from rocket_propulsion.heat_transfer.regen_chamber import Coolant, CoolingChannels


class MutableLength(float):
    def __new__(cls, value):
        result = super().__new__(cls, value)
        result.scale = 1.
        return result

    def __mul__(self, other):
        return float(self)*self.scale*other


class RegenInputRecordTests(unittest.TestCase):
    def test_owned_width_cannot_be_changed_by_caller_arithmetic(self):
        width = MutableLength(.0015)
        channel = CoolingChannels(100, width, .005, .001)
        before = channel.width_m*channel.height_m
        width.scale = 1000.
        self.assertEqual(channel.width_m*channel.height_m, before)
        self.assertIs(type(channel.width_m), float)

    def test_coolant_direction_requires_an_actual_boolean(self):
        for value in ('false', 0, 1, None):
            with self.subTest(value=value), self.assertRaises(InputError) as error:
                Coolant('methane', 120., 1.5e7, 2., counterflow=value)
            self.assertEqual(error.exception.field, 'counterflow')

    def test_direct_physical_fields_refuse_boolean_text_and_erased_reals(self):
        configs = [(CoolingChannels, {'count':100,'width_m':.0015,'height_m':.005,
                                     'hot_wall_thickness_m':.001,'roughness_m':0.,'helix_angle_deg':0.}),
                   (Coolant, {'fluid':'methane','inlet_temperature_k':120.,
                              'inlet_pressure_pa':1.5e7,'mass_flow_kg_s':2.})]
        for kind, valid in configs:
            for name in valid.keys()-{'fluid'}:
                for value in (True, '1', float('nan'), float('inf'), Fraction(1, 10**500)):
                    with (self.subTest(record=kind.__name__, field=name, value=repr(value)[:30]),
                          self.assertRaises(DomainError)):
                        kind(**(valid|{name:value}))

    def test_fractional_geometry_and_native_direction_remain_supported(self):
        channel = CoolingChannels(1, Fraction(3,2000), Fraction(1,200), Fraction(1,1000))
        for item in fields(channel):
            self.assertIs(type(getattr(channel,item.name)), int if item.name=='count' else float)
        for direction in (True, False):
            fluid = Coolant('methane', Fraction(120), 15000000, Fraction(2), counterflow=direction)
            self.assertIs(fluid.counterflow, direction)
            self.assertEqual(fluid.mass_flow_kg_s, 2.)
