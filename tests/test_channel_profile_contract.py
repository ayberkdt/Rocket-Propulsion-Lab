"""Owned profile geometry, independent interpolation and real coupled controls."""
import math
import unittest
from fractions import Fraction

from rocket_propulsion.core.errors import DomainError
from rocket_propulsion.heat_transfer.regen_chamber import ChannelProfile, solve_regen_chamber

BASE = {'x_m': (0., .5, 1.), 'width_m': (.001, .002, .003),
        'height_m': (.005, .005, .005), 'hot_wall_thickness_m': (.001, .001, .001),
        'count': (100, 120, 60)}


class MutableReal(float):
    def __new__(cls, value):
        obj = super().__new__(cls, value)
        obj.factor = 1.
        return obj

    def __float__(self):
        return super().__float__() * self.factor


class MutableCount(int):
    def __new__(cls, value):
        obj = super().__new__(cls, value)
        obj.factor = 1.
        return obj

    def __rtruediv__(self, other):
        return other / int(self) * self.factor


class ChannelProfileContractTests(unittest.TestCase):
    def test_single_pass_columns_retain_interpolation_and_nodes(self):
        data = {k: iter(v) if k != 'count' else v for k, v in BASE.items()}
        try:
            profile = ChannelProfile(**data)
            answer = profile.geometry_at(.25)
        except (ValueError, IndexError) as exc:
            self.fail(f'Valid single-pass column lost its values: {exc}')
        self.assertEqual(profile.x_m, BASE['x_m'])
        self.assertAlmostEqual(answer[0], .0015, places=15)
        self.assertEqual(answer[1:], (.005, .001))

    def test_columns_refuse_invalid_scalars_with_correct_field(self):
        invalid = (True, False, '0.001', None, complex(1), float('nan'),
                   float('inf'), Fraction(10**400), Fraction(1, 10**500))
        for name in ('x_m', 'width_m', 'height_m', 'hot_wall_thickness_m'):
            for bad in invalid:
                column = (bad, *BASE[name][1:])
                with self.subTest(field=name, value=repr(bad)[:30]):
                    with self.assertRaises(DomainError) as caught:
                        ChannelProfile(**(BASE | {name: column}))
                    self.assertIn(name, caught.exception.field or '')

    def test_columns_refuse_unordered_or_noniterable_containers(self):
        for name in ('x_m', 'width_m', 'height_m', 'hot_wall_thickness_m'):
            for bad in ({0., .5, 1.}, {'a': .1}, '012', b'012', None, 3):
                with self.subTest(field=name, value=repr(bad)):
                    with self.assertRaises(DomainError) as caught:
                        ChannelProfile(**(BASE | {name: bad}))
                    self.assertIn(name, caught.exception.field or '')

    def test_fraction_columns_and_mutated_lists_keep_physical_geometry(self):
        x = [Fraction(0), Fraction(1, 2), Fraction(1)]
        widths = [Fraction(1, 1000), Fraction(2, 1000), Fraction(3, 1000)]
        profile = ChannelProfile(**(BASE | {'x_m': x, 'width_m': widths}))
        x.clear()
        widths[0] = 900.
        self.assertAlmostEqual(profile.geometry_at(.25)[0], .0015, places=15)
        self.assertTrue(all(type(v) is float for v in profile.x_m + profile.width_m))

    def test_scalar_roughness_and_helix_do_not_retain_caller_conversion(self):
        rough, helix = MutableReal(1e-6), MutableReal(10.)
        profile = ChannelProfile(**BASE, roughness_m=rough, helix_angle_deg=helix)
        rough.factor = 100.
        helix.factor = 2.
        self.assertEqual(float(profile.roughness_m), 1e-6)
        self.assertEqual(float(profile.helix_angle_deg), 10.)
        for name in ('roughness_m', 'helix_angle_deg'):
            for bad in (True, '1', None, float('nan'), Fraction(1, 10**500)):
                with self.subTest(field=name, value=repr(bad)[:20]):
                    with self.assertRaises(DomainError) as caught:
                        ChannelProfile(**BASE, **{name: bad})
                    self.assertIn(name, caught.exception.field or '')

    def test_native_counts_prevent_post_validation_flow_divisor_mutation(self):
        for piecewise in (False, True):
            count = MutableCount(100)
            raw = (count, count, count) if piecewise else count
            profile = ChannelProfile(**(BASE | {'count': raw}))
            count.factor = 2.
            self.assertEqual(10. / profile.count_at(.25), .1)
            self.assertIs(type(profile.count_at(.25)), int)
        for count in (True, 1.5, (100, True, 60), (100, 0, 60), (100, 60)):
            with self.assertRaises(DomainError):
                ChannelProfile(**(BASE | {'count': count}))

    def test_direct_queries_refuse_malformed_coordinates_even_for_constant_count(self):
        for count in (100, (100, 120, 60)):
            profile = ChannelProfile(**(BASE | {'count': count}))
            for method in (profile.geometry_at, profile.count_at):
                for bad in (True, '0.5', None, float('nan'), float('inf'),
                            Fraction(1, 10**500)):
                    with self.subTest(method=method.__name__, value=repr(bad)[:20]):
                        with self.assertRaises(DomainError) as caught:
                            method(bad)
                        self.assertEqual(caught.exception.field, 'x_m')

    def test_affine_geometry_survives_unrepresentable_axial_span(self):
        profile = ChannelProfile(x_m=(-1e308, 1e308), width_m=(.001, .003),
                                 height_m=(.003, .007), hot_wall_thickness_m=(.001, .005),
                                 count=100)
        for x, expected in ((-5e307, (.0015, .004, .002)),
                            (0., (.002, .005, .003)), (5e307, (.0025, .006, .004))):
            for value, reference in zip(profile.geometry_at(x), expected, strict=True):
                self.assertAlmostEqual(value, reference, places=15)

    def test_bifurcation_sides_nodes_and_endpoint_clamping(self):
        profile = ChannelProfile(**BASE)
        for x, expected in ((-1., 100), (0., 100), (math.nextafter(.5, 0.), 100),
                            (.5, 120), (1., 60), (2., 60)):
            self.assertEqual(profile.count_at(x), expected)
        for x in (0., .25, .5, .75, 1.):
            self.assertAlmostEqual(profile.geometry_at(x)[0], .001 + .002*x, places=15)
        self.assertEqual(profile.geometry_at(-1e-9)[0], .001)
        self.assertEqual(profile.geometry_at(1 + 1e-9)[0], .003)
        for x in (-2e-9, 1 + 2e-9):
            with self.assertRaises(DomainError):
                profile.geometry_at(x)

    def test_tiny_barycentric_terms_are_retained_before_geometry_scaling(self):
        for xs, widths, coordinate in (((0., 1e300), (1e-300, 1e300), 1e-300),
                                       ((-1e300, 0.), (1e300, 1e-300), -1e-300)):
            profile = ChannelProfile(x_m=xs, width_m=widths, height_m=(.005, .005),
                                     hot_wall_thickness_m=(.001, .001), count=100)
            # In both orientations the tiny barycentric contribution is
            # 1e-600 * 1e300 = 1e-300, in addition to the small endpoint.
            self.assertTrue(math.isclose(profile.geometry_at(coordinate)[0], 2e-300,
                                         rel_tol=3e-15))

    def test_real_coupled_solution_is_invariant_to_caller_helix_mutation(self):
        from tests.test_regen_channel_geometry import contour, coolant, gas

        nozzle = contour()
        xs = [s.x_m for s in nozzle.stations]
        helix = MutableReal(10.)
        profile = ChannelProfile(x_m=(min(xs), max(xs)), width_m=(.0015, .0015),
                                 height_m=(.005, .005), hot_wall_thickness_m=(.001, .001),
                                 count=100, roughness_m=1e-6, helix_angle_deg=helix)
        def evaluate():
            return solve_regen_chamber(contour=nozzle, gas=gas(), channels=profile,
                                       coolant=coolant())
        before = evaluate()
        helix.factor = 2.
        after = evaluate()
        self.assertGreater(before.total_heat_load_w, 0.)
        for name in ('total_heat_load_w', 'coolant_outlet_temperature_k',
                     'coolant_outlet_pressure_pa', 'maximum_hot_wall_temperature_k'):
            a, b = getattr(before, name), getattr(after, name)
            self.assertTrue(math.isclose(a, b, rel_tol=1e-12), (name, a, b))
