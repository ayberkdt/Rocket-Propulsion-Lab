"""Known-invalid regen controls must fail before expensive gas chemistry.

The sentinel proves evaluation order only. Real physics and energy-balance
controls remain in test_api_thermal and the regenerative solver suites.
"""
import math
import unittest
from unittest.mock import patch

from rocket_propulsion.api.routes import dispatch_calculation
from rocket_propulsion.core.errors import RocketPropulsionError


class RegenAdmissionOrderTests(unittest.TestCase):
    def test_invalid_controls_are_rejected_before_chemistry(self):
        cases = {
            'chamber_pressure_pa': 0., 'area_ratio': .5,
            'channel_count': 3, 'channel_width_m': 0., 'channel_height_m': 0.,
            'hot_wall_thickness_m': 0., 'roughness_m': -1.,
            'channel_helix_angle_deg': 75., 'coolant_inlet_temperature_k': 0.,
            'coolant_inlet_pressure_pa': 0., 'coolant_mass_flow_kg_s': 0.,
            'coolant': 'unknown-fluid', 'wall_material': 'unknown-wall',
            'wall_model': 'unknown-model', 'closeout_thickness_m': -.001,
            'axial_conduction': 'not-a-boolean',
        }
        for field, value in cases.items():
            with (self.subTest(field=field),
                  patch('rocket_propulsion.api.thermal_routes.hot_gas_from_propellant',
                        side_effect=AssertionError('Expensive chemistry reached for invalid controls.')) as gas):
                with self.assertRaises(RocketPropulsionError):
                    dispatch_calculation('/api/v1/regen-chamber', {field: value})
                gas.assert_not_called()

    def test_incompatible_closeout_is_rejected_before_chemistry(self):
        with patch('rocket_propulsion.api.thermal_routes.hot_gas_from_propellant',
                   side_effect=AssertionError('Chemistry reached for incompatible wall controls.')) as gas:
            with self.assertRaises(RocketPropulsionError) as error:
                dispatch_calculation('/api/v1/regen-chamber', {'wall_model':'fin', 'closeout_thickness_m':.001})
            self.assertEqual(error.exception.field, 'closeout_thickness_m')
            gas.assert_not_called()

    def test_declared_coolant_range_is_checked_before_chemistry(self):
        cases = [('coolant_inlet_temperature_k', math.nextafter(90.6941, 0.)),
                 ('coolant_inlet_temperature_k', math.nextafter(625., math.inf)),
                 ('coolant_inlet_pressure_pa', math.nextafter(1e9, math.inf))]
        for field, value in cases:
            with (self.subTest(field=field, value=value),
                  patch('rocket_propulsion.api.thermal_routes.hot_gas_from_propellant',
                        side_effect=AssertionError('Chemistry reached outside the methane EOS range.')) as gas):
                with self.assertRaises(RocketPropulsionError) as error:
                    dispatch_calculation('/api/v1/regen-chamber', {field:value})
                self.assertEqual(error.exception.field, field)
                gas.assert_not_called()

    def test_width_limit_uses_channel_normal_pitch_and_wall_thickness(self):
        for angle in (0., 60.):
            pitch = 2*math.pi*(math.sqrt(.01/math.pi)+.001)*math.cos(math.radians(angle))/100
            with (self.subTest(angle=angle),
                  patch('rocket_propulsion.api.thermal_routes.hot_gas_from_propellant',
                        side_effect=AssertionError('Chemistry reached for overlapping channels.')) as gas):
                with self.assertRaises(RocketPropulsionError) as error:
                    dispatch_calculation('/api/v1/regen-chamber', {
                        'channel_helix_angle_deg':angle,'channel_width_m':pitch})
                self.assertEqual(error.exception.field, 'channel_width_m')
                gas.assert_not_called()

    def test_valid_control_combinations_reach_the_real_chemistry_boundary(self):
        class ReachedChemistry(Exception):
            pass
        for payload in ({}, {'coolant_mass_flow_kg_s':None}, {'coolant_mass_flow_kg_s':''},
                        {'wall_model':'2d','closeout_thickness_m':.001},
                        {'axial_conduction':'false'}, {'channel_helix_angle_deg':60.,'channel_width_m':.001}):
            with (self.subTest(payload=payload),
                  patch('rocket_propulsion.api.thermal_routes.hot_gas_from_propellant',
                        side_effect=ReachedChemistry) as gas):
                with self.assertRaises(ReachedChemistry):
                    dispatch_calculation('/api/v1/regen-chamber', payload)
                gas.assert_called_once()

    def test_positive_land_and_declared_eos_endpoints_are_not_rejected(self):
        class ReachedChemistry(Exception):
            pass
        valid = [{'coolant_inlet_temperature_k':90.6941},
                 {'coolant_inlet_temperature_k':625.}, {'coolant_inlet_pressure_pa':1e9}]
        for angle in (0., 60.):
            pitch = 2*math.pi*(math.sqrt(.01/math.pi)+.001)*math.cos(math.radians(angle))/100
            valid.append({'channel_helix_angle_deg':angle,'channel_width_m':math.nextafter(pitch, 0.)})
        for payload in valid:
            with (self.subTest(payload=payload),
                  patch('rocket_propulsion.api.thermal_routes.hot_gas_from_propellant',
                        side_effect=ReachedChemistry) as gas):
                try:
                    dispatch_calculation('/api/v1/regen-chamber', payload)
                except ReachedChemistry:
                    pass
                except RocketPropulsionError as error:
                    self.fail(f'Positive land or a declared EOS endpoint was rejected: {error}')
                else:
                    self.fail('Valid preflight did not reach the chemistry boundary.')
                gas.assert_called_once()


if __name__ == '__main__':
    unittest.main()
