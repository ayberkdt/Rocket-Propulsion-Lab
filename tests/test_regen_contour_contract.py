"""Adversarial cooling geometry admission and real-march ownership controls.

Mocks mark only the pre-EOS/pre-axial boundary. Physical controls use the actual
equilibrium fixture, methane EOS and coupled wall march. A geometric check does
not validate the prescribed Mach field or the accuracy of a continuous wall.
"""
import json
import math
import unittest
from dataclasses import asdict, replace
from fractions import Fraction
from functools import cache
from unittest.mock import patch

from rocket_propulsion.core.errors import DomainError
from rocket_propulsion.geometry import generate_moc_nozzle, generate_rao_contour
from rocket_propulsion.heat_transfer.regen_chamber import solve_regen_chamber
from rocket_propulsion.thermodynamics import snapshot_thermal_contour
from rocket_propulsion.thermodynamics.real_fluids import reference_fluid
from tests.test_regen_channel_geometry import CONSTANT, contour, coolant, gas
from tests.test_regen_solver_contract import MutableScalar, fixture


@cache
def nominal(counterflow=True):
    return solve_regen_chamber(contour=contour(), gas=gas(), channels=CONSTANT,
                               coolant=replace(coolant(), counterflow=counterflow))


def changed_station(source, **values):
    return replace(source, stations=(replace(source.stations[0], **values),
                                      *source.stations[1:]))


class ReachedFluid(Exception):
    pass


class RegenContourContractTests(unittest.TestCase):
    def refuse_before_work(self, source, axial):
        with (patch('rocket_propulsion.heat_transfer.regen_chamber.reference_fluid',
                    side_effect=AssertionError('Invalid contour reached EOS')) as fluid,
              patch('rocket_propulsion.heat_transfer.regen_chamber._with_axial_conduction',
                    side_effect=AssertionError('Invalid contour reached axial iteration')) as dispatch):
            with self.assertRaises(DomainError) as caught:
                solve_regen_chamber(**(fixture() | {'contour': source}), axial_conduction=axial)
            self.assertIsInstance(caught.exception.field, str)
            self.assertTrue(caught.exception.field.strip())
            fluid.assert_not_called()
            dispatch.assert_not_called()

    def test_malformed_consumed_quantities_fail_before_eos_and_axial_work(self):
        bad = (True, False, None, '1', [], math.nan, math.inf, -math.inf,
               Fraction(1, 10**500))
        for axial in (False, True):
            for field in ('x_m', 'radius_m', 'area_ratio', 'mach', 'throat_radius_m'):
                for value in bad:
                    with self.subTest(axial=axial, field=field, value=repr(value)[:25]):
                        c = fixture()['contour']
                        c = (replace(c, throat_radius_m=value) if field == 'throat_radius_m'
                             else changed_station(c, **{field: value}))
                        self.refuse_before_work(c, axial)

    def test_invalid_station_container_or_record_fails_before_work(self):
        c = fixture()['contour']
        for axial in (False, True):
            for points in (None, (), c.stations[:1], 'bad', b'bad', set(c.stations),
                           {0: c.stations[0]}, (c.stations[0], None),
                           (c.stations[0], asdict(c.stations[1]))):
                with self.subTest(axial=axial, points=repr(points)[:40]):
                    self.refuse_before_work(replace(c, stations=points), axial)

    def test_duplicate_or_reversed_positions_fail_before_work(self):
        c = fixture()['contour']
        for axial in (False, True):
            for points in (tuple(reversed(c.stations)),
                           (c.stations[0], replace(c.stations[1], x_m=c.stations[0].x_m),
                            *c.stations[2:])):
                with self.subTest(axial=axial, x=points[1].x_m):
                    self.refuse_before_work(replace(c, stations=points), axial)

    def test_incoherent_geometry_and_negative_physical_values_fail_before_work(self):
        c = fixture()['contour']
        cases = [changed_station(c, radius_m=c.stations[0].radius_m*2),
                 changed_station(c, area_ratio=100.), replace(c, throat_radius_m=-1.),
                 replace(c, throat_radius_m=c.throat_radius_m*2),
                 changed_station(c, radius_m=0.), changed_station(c, radius_m=-1.),
                 changed_station(c, area_ratio=0.), changed_station(c, mach=-1.)]
        for axial in (False, True):
            for n, source in enumerate(cases):
                with self.subTest(axial=axial, case=n):
                    self.refuse_before_work(source, axial)

    def test_subunity_area_ratio_is_not_clamped_even_with_overlapping_cells(self):
        # Representation overlap is weaker than the coupled Bartz A/At >= 1
        # domain. The one-ULP-under value must fail before either solver path.
        c = fixture()['contour']
        points = tuple(replace(s, radius_m=1., area_ratio=math.nextafter(1., 0.), mach=0.)
                       for s in c.stations[:2])
        c = replace(c, throat_radius_m=1., stations=points)
        snapshot_thermal_contour(c)  # geometric rounding cells overlap
        for axial in (False, True):
            with self.subTest(axial=axial):
                self.refuse_before_work(c, axial)

    def test_axial_dispatch_owns_container_stations_and_native_scalars(self):
        original = fixture()['contour']
        raw_scalars = []
        points = []
        for point in original.stations:
            values = {}
            for name in ('x_m', 'radius_m', 'area_ratio', 'mach'):
                values[name] = MutableScalar(getattr(point, name))
                raw_scalars.append(values[name])
            points.append(replace(point, **values))
        radius = MutableScalar(original.throat_radius_m)
        source = replace(original, stations=points, throat_radius_m=radius)
        with patch('rocket_propulsion.heat_transfer.regen_chamber._with_axial_conduction',
                   side_effect=lambda **kw: kw):
            received = solve_regen_chamber(**(fixture() | {'contour': source}),
                                            axial_conduction=True)['contour']
        points.clear()
        for scalar in [radius, *raw_scalars]:
            scalar.factor = 100.
        self.assertIs(type(received.stations), tuple)
        self.assertIs(type(received.throat_radius_m), float)
        self.assertEqual(received.throat_radius_m, original.throat_radius_m)
        self.assertEqual(len(received.stations), len(original.stations))
        for point, expected in zip(received.stations, original.stations, strict=True):
            for name in ('x_m', 'radius_m', 'area_ratio', 'mach'):
                with self.subTest(x=expected.x_m, field=name):
                    self.assertIs(type(getattr(point, name)), float)
                    self.assertEqual(getattr(point, name), getattr(expected, name))

    def test_ordered_iterator_is_snapshotted_once_before_fluid_work(self):
        source = fixture()['contour']
        calls = []

        class OnePass:
            def __iter__(self):
                calls.append(1)
                if len(calls) > 1:
                    raise AssertionError('Caller-owned stations were revisited')
                yield from source.stations

        with (patch('rocket_propulsion.heat_transfer.regen_chamber.reference_fluid',
                    side_effect=ReachedFluid), self.assertRaises(ReachedFluid)):
            solve_regen_chamber(**(fixture() | {'contour': replace(source, stations=OnePass())}))
        self.assertEqual(len(calls), 1)

    def test_malformed_iteration_is_typed_but_unrelated_callback_failure_propagates(self):
        for error_type in (TypeError, ValueError, OverflowError, RuntimeError):
            class BrokenIterator:
                def __init__(self, failure):
                    self.failure = failure

                def __iter__(self):
                    raise self.failure('deliberate iterator failure')
            source = replace(fixture()['contour'], stations=BrokenIterator(error_type))
            with self.subTest(error=error_type.__name__):
                if error_type is RuntimeError:
                    with self.assertRaisesRegex(RuntimeError, 'deliberate iterator failure'):
                        solve_regen_chamber(**(fixture() | {'contour': source}))
                else:
                    self.refuse_before_work(source, False)

    def test_real_march_survives_caller_mutation_at_the_eos_boundary(self):
        original = contour()
        raw_radius = MutableScalar(original.throat_radius_m)
        raw_mach = MutableScalar(original.stations[0].mach)
        points = [replace(original.stations[0], mach=raw_mach), *original.stations[1:]]
        source = replace(original, throat_radius_m=raw_radius, stations=points)

        def mutate_then_get_fluid(name):
            points.clear()
            raw_radius.factor = 3.
            raw_mach.factor = 3.
            return reference_fluid(name)

        with patch('rocket_propulsion.heat_transfer.regen_chamber.reference_fluid',
                   side_effect=mutate_then_get_fluid):
            actual = solve_regen_chamber(contour=source, gas=gas(), channels=CONSTANT,
                                         coolant=coolant())
        self.assertEqual(actual, nominal())
        json.dumps(asdict(actual), allow_nan=False)

    def check_energy(self, result, flow):
        fluid = reference_fluid(flow.fluid)
        inlet = fluid.state_tp(flow.inlet_temperature_k, flow.inlet_pressure_pa).enthalpy_j_kg
        outlet = fluid.state_tp(result.coolant_outlet_temperature_k,
                                result.coolant_outlet_pressure_pa).enthalpy_j_kg
        expected = flow.mass_flow_kg_s * (outlet-inlet)
        self.assertGreater(result.total_heat_load_w, 0.)
        self.assertLess(abs(expected-result.total_heat_load_w), 1e-8*expected)

    def test_real_counterflow_and_coflow_match_the_shared_owned_profile(self):
        for counterflow in (False, True):
            with self.subTest(counterflow=counterflow):
                flow = replace(coolant(), counterflow=counterflow)
                actual = solve_regen_chamber(contour=snapshot_thermal_contour(contour()),
                                             gas=gas(), channels=CONSTANT, coolant=flow)
                self.assertEqual(actual, nominal(counterflow))
                self.check_energy(actual, flow)

    def test_generated_hermite_rao_and_both_moc_families_run_with_real_fluids(self):
        from rocket_propulsion.geometry import generate_nozzle_contour

        families = [generate_nozzle_contour(throat_area_m2=.01, area_ratio=10., gamma=gas().gamma,
                                            station_count=41),
                    generate_rao_contour(throat_area_m2=.01, area_ratio=10., gamma=gas().gamma,
                                         station_count=41)]
        families.extend(generate_moc_nozzle(throat_area_m2=.01, area_ratio=10., gamma=gas().gamma,
                                           characteristic_count=10,
                                           axisymmetric_correction=axis).contour
                        for axis in (False, True))
        for source in families:
            with self.subTest(family=source.contour):
                result = solve_regen_chamber(contour=source, gas=gas(), channels=CONSTANT,
                                             coolant=coolant())
                self.assertEqual(tuple(s.x_m for s in result.stations),
                                 tuple(s.x_m for s in source.stations))
                self.check_energy(result, coolant())

    def test_partial_prescribed_profile_keeps_mach_zero_and_no_throat_sample(self):
        source = contour()
        # A thermal span may omit the throat and its Mach field is prescribed.
        # Do not impose an artificial isentropic solution or demand x=0.
        source = replace(source, stations=tuple(replace(p, mach=0.) for p in source.stations[-3:]))
        result = solve_regen_chamber(contour=source, gas=gas(), channels=CONSTANT, coolant=coolant())
        self.assertEqual(len(result.stations), 3)
        self.assertTrue(all(s.mach == 0. for s in result.stations))
        self.check_energy(result, coolant())

    def test_unused_metadata_is_not_silently_claimed_as_checked_physics(self):
        source = contour()
        source = replace(source, total_length_m=math.nan, throat_area_m2='unused',
                         surface_area_m2=-1., stations=tuple(replace(p, area_m2=math.nan,
                            pressure_total_ratio=False, temperature_total_ratio='unused')
                            for p in source.stations))
        actual = solve_regen_chamber(contour=source, gas=gas(), channels=CONSTANT, coolant=coolant())
        self.assertEqual(actual, nominal())

    def test_real_axial_solve_reuses_the_owned_profile_and_reports_its_model(self):
        from tests.test_regen_conjugate import CHANNELS, solved
        from tests.test_regen_conjugate import contour as axial_contour

        result = solve_regen_chamber(contour=snapshot_thermal_contour(axial_contour()),
                                     gas=gas(), channels=CHANNELS, coolant=coolant(),
                                     axial_conduction=True)
        self.assertEqual(result, solved(True))
        codes = {w.code: w.message for w in result.warnings}
        self.assertIn('axial_wall_conduction', codes)
        self.assertIn('axial wall conduction coupled', codes['rte_class_model'])
        self.assertNotIn('no axial wall conduction', codes['rte_class_model'])
        self.assertIn('not validated against hardware data', codes['rte_class_model'])
        self.check_energy(result, coolant())


if __name__ == '__main__':
    unittest.main()
