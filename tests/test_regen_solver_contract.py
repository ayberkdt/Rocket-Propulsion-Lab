"""Thermal option admission, ownership and material-evidence regression controls."""
import json
import math
import unittest
from dataclasses import asdict
from fractions import Fraction
from functools import cache
from unittest.mock import patch

from rocket_propulsion.core.errors import RocketPropulsionError
from rocket_propulsion.geometry import generate_nozzle_contour
from rocket_propulsion.heat_transfer.regen_chamber import (
    Coolant,
    CoolingChannels,
    HotGas,
    solve_regen_chamber,
)
from rocket_propulsion.thermodynamics.materials import get_wall_material


class MutableScalar(float):
    def __new__(cls, value):
        obj = super().__new__(cls, value)
        obj.factor = 1.
        return obj

    def __float__(self):
        return super().__float__() * self.factor


class PretendOption(str):
    def __eq__(self, other):
        return True

    __hash__ = str.__hash__


@cache
def fixture():
    gas = HotGas(7e6, 3500., 1800., 3200., 1e-4, .7, 1.2, .05)
    contour = generate_nozzle_contour(contour='conical', throat_area_m2=.01,
                                      area_ratio=2., gamma=gas.gamma, station_count=41)
    return {'gas': gas, 'contour': contour,
            'channels': CoolingChannels(100, .0015, .005, .001),
            'coolant': Coolant('methane', 120., 1.5e7, 2.)}


@cache
def evaluated(material='grcop-84', override=None):
    # Reuse the validated LOX/CH4 regression point, including its equilibrium
    # gas properties and fuel flow. Admission-only tests above stop before EOS work.
    from tests.test_regen_channel_geometry import CONSTANT, contour, coolant, gas

    return solve_regen_chamber(contour=contour(), gas=gas(), coolant=coolant(),
                               channels=CONSTANT, wall_material=material,
                               wall_thermal_conductivity_w_m_k=override)


class RegenSolverContractTests(unittest.TestCase):
    def test_optional_numbers_are_rejected_before_fluid_work_even_when_inactive(self):
        for field in ('wall_thermal_conductivity_w_m_k', 'closeout_thickness_m',
                      'closeout_conductivity_w_m_k'):
            values = [True, False, '0', [], float('nan'), float('inf'), -1., Fraction(1, 10**500)]
            values += [None] if field == 'closeout_thickness_m' else [0.]
            for value in values:
                with (self.subTest(field=field, value=repr(value)[:25]),
                      patch('rocket_propulsion.heat_transfer.regen_chamber.reference_fluid',
                            side_effect=AssertionError('Invalid option reached fluid work')) as fluid):
                    with self.assertRaises(RocketPropulsionError) as caught:
                        solve_regen_chamber(**fixture(), **{field: value})
                    self.assertEqual(caught.exception.field, field)
                    fluid.assert_not_called()

    def test_bad_options_are_rejected_before_axial_dispatch(self):
        for field, value in [('gas_radiation', 'invalid'), ('wall_model', 'invalid'),
                             ('coolant_correlation', 'invalid'), ('closeout_conductivity_w_m_k', True),
                             ('wall_thermal_conductivity_w_m_k', '305'),
                             ('closeout_thickness_m', False), ('wall_material', ''),
                             ('wall_material', 'unknown')]:
            with (self.subTest(field=field),
                  patch('rocket_propulsion.heat_transfer.regen_chamber._with_axial_conduction',
                        side_effect=AssertionError('Invalid option reached axial solve')) as axial):
                with self.assertRaises(RocketPropulsionError):
                    solve_regen_chamber(**fixture(), axial_conduction=True, **{field: value})
                axial.assert_not_called()

    def test_options_use_owned_text_not_subclass_equality(self):
        for field in ('coolant_correlation', 'gas_radiation', 'wall_model'):
            with (self.subTest(field=field),
                  patch('rocket_propulsion.heat_transfer.regen_chamber.reference_fluid',
                        side_effect=AssertionError('Pretended option reached fluid'))):
                with self.assertRaises(RocketPropulsionError) as caught:
                    solve_regen_chamber(**fixture(), **{field: PretendOption('unknown')})
                self.assertEqual(caught.exception.field, field)

    def test_required_records_fail_with_their_field_before_fluid_or_axial_work(self):
        for field in ('contour', 'gas', 'channels', 'coolant'):
            for value in (None, True, {}):
                with (self.subTest(field=field, value=value),
                      patch('rocket_propulsion.heat_transfer.regen_chamber.reference_fluid',
                            side_effect=AssertionError('Invalid record reached fluid')),
                      patch('rocket_propulsion.heat_transfer.regen_chamber._with_axial_conduction',
                            side_effect=AssertionError('Invalid record reached axial'))):
                    with self.assertRaises(RocketPropulsionError) as caught:
                        solve_regen_chamber(**(fixture() | {field: value}), axial_conduction=True)
                    self.assertEqual(caught.exception.field, field)

    def test_axial_dispatch_receives_owned_normalized_scalars(self):
        raw = {name: MutableScalar(value) for name, value in {
            'wall_thermal_conductivity_w_m_k': 305., 'closeout_thickness_m': .001,
            'closeout_conductivity_w_m_k': 200.}.items()}
        with patch('rocket_propulsion.heat_transfer.regen_chamber._with_axial_conduction',
                   side_effect=lambda **kw: kw) as axial:
            captured = solve_regen_chamber(**fixture(), wall_model='2d',
                                           axial_conduction=True, **raw)
        axial.assert_called_once()
        for field, value in raw.items():
            expected = float(value)
            value.factor = 10.
            self.assertIs(type(captured[field]), float)
            self.assertEqual(captured[field], expected)

    def test_valid_zero_closeout_and_optional_none_reach_fluid_work(self):
        class ReachedFluid(Exception):
            pass
        for extra in ({}, {'closeout_thickness_m': 0.}, {'closeout_conductivity_w_m_k': None},
                      {'wall_thermal_conductivity_w_m_k': None},
                      {'closeout_conductivity_w_m_k': 200.},
                      {'wall_material': None, 'wall_thermal_conductivity_w_m_k': 305.}):
            with (self.subTest(extra=extra),
                  patch('rocket_propulsion.heat_transfer.regen_chamber.reference_fluid',
                        side_effect=ReachedFluid) as fluid):
                with self.assertRaises(ReachedFluid):
                    solve_regen_chamber(**fixture(), **extra)
                fluid.assert_called_once()

    def test_grcop_source_and_conditional_limit_survive_the_real_solver(self):
        result = evaluated()
        use = result.wall_property_use
        self.assertIsNotNone(use)
        self.assertEqual(use.material_key, 'grcop-84')
        self.assertEqual(use.conductivity_basis, 'catalog_reference_constant')
        self.assertEqual(use.conductivity_w_m_k, 305.)
        self.assertEqual(use.conductivity_reference, get_wall_material('grcop-84').thermal_conductivity)
        self.assertEqual(use.conductivity_reference.upper_value, 320.)
        self.assertEqual(use.conductivity_reference.citation.source_id, 'NASA-TM-2005-213566')
        self.assertEqual(use.limit_basis, 'maximum_service_temperature')
        self.assertEqual(use.limit_reference.value, 973.)
        self.assertIn('100 MPa', use.limit_reference.condition)
        self.assertEqual(use.limit_reference.citation.source_id, 'NASA-CR-2000-210055')
        notes = {n.code: n for n in result.warnings}
        self.assertIn('wall_conductivity_constant', notes)
        self.assertIn('100 MPa', notes['wall_limit_reference_only'].message)
        self.assertGreater(result.total_heat_load_w, 0.)

    def test_override_does_not_claim_the_catalog_conductivity_source(self):
        implicit, explicit = evaluated(), evaluated(override=305.)
        use = explicit.wall_property_use
        self.assertEqual(use.conductivity_basis, 'user_constant')
        self.assertIsNone(use.conductivity_reference)
        self.assertEqual(use.limit_reference, implicit.wall_property_use.limit_reference)
        for field in ('total_heat_load_w', 'coolant_outlet_temperature_k',
                      'coolant_outlet_pressure_pa', 'maximum_hot_wall_temperature_k'):
            self.assertEqual(getattr(implicit, field), getattr(explicit, field))

    def test_no_material_does_not_invent_a_temperature_limit(self):
        result = evaluated(material=None, override=305.)
        use = result.wall_property_use
        self.assertIsNone(use.material_key)
        self.assertIsNone(use.limit_basis)
        self.assertIsNone(use.limit_reference)
        self.assertIsNone(result.wall_limit_temperature_k)
        self.assertTrue(all(s.wall_limit_margin_k is None for s in result.stations))
        self.assertIn('wall_limit_not_assessed', {n.code for n in result.warnings})
        self.assertTrue(math.isclose(result.total_heat_load_w, evaluated().total_heat_load_w,
                                     rel_tol=1e-13))

    def test_copper_melting_limit_is_not_relabelled_as_service_limit(self):
        result = evaluated(material='ofhc-copper')
        use = result.wall_property_use
        self.assertEqual(use.limit_basis, 'solidus_temperature')
        self.assertEqual(use.limit_reference.name, 'solidus_temperature')
        note = next(n for n in result.warnings if n.code == 'wall_limit_reference_only')
        self.assertIn('solidus_temperature', note.message)

    def test_material_evidence_is_strict_json_in_the_http_result(self):
        from rocket_propulsion.api.routes import dispatch_calculation
        from tests.test_regen_channel_geometry import gas

        with patch('rocket_propulsion.api.thermal_routes.hot_gas_from_propellant',
                   return_value=gas()):
            response = dispatch_calculation('/api/v1/regen-chamber', {})
        self.assertEqual(response['wall_property_use']['limit_reference']['value'], 973.)
        self.assertIn('100 MPa', response['wall_property_use']['limit_reference']['condition'])
        json.dumps(response, allow_nan=False)
        restored = json.loads(json.dumps(asdict(evaluated()), allow_nan=False))
        self.assertEqual(restored['wall_property_use']['conductivity_reference']['value'], 305.)


if __name__ == '__main__':
    unittest.main()
