"""Generator admission, sonic limiting geometry and independent physics checks.

Sentinels verify ordering only. Valid constructions are checked with independent
area-Mach equations, analytic dimensions, parsed exports and real-fluid energy.
"""
import json
import math
import unittest
from dataclasses import asdict, fields, replace
from fractions import Fraction
from unittest.mock import patch

import mpmath as mp

from rocket_propulsion.core.errors import DomainError
from rocket_propulsion.geometry import (
    export_nozzle_contour,
    generate_moc_nozzle,
    generate_nozzle_contour,
    generate_rao_contour,
)
from rocket_propulsion.heat_transfer.regen_chamber import solve_regen_chamber
from rocket_propulsion.thermodynamics import snapshot_thermal_contour

GENERATORS = (generate_nozzle_contour, generate_rao_contour, generate_moc_nozzle)
BASE = {'throat_area_m2': .01, 'area_ratio': 10., 'gamma': 1.22}


class ArithmeticTrap(float):
    def __mul__(self, other):
        raise AssertionError('Unowned input arithmetic')

    __rmul__ = __mul__
    __truediv__ = __mul__
    __rtruediv__ = __mul__
    __add__ = __mul__
    __radd__ = __mul__
    __sub__ = __mul__
    __rsub__ = __mul__
    __pow__ = __mul__


class TextTrap(str):
    def strip(self, *args):
        raise AssertionError('Unowned contour name')


def inputs(generator):
    return BASE | ({'characteristic_count': 5} if generator is generate_moc_nozzle
                   else {'station_count': 41})


class GeneratorContractTests(unittest.TestCase):
    def valid(self, function, **kwargs):
        try:
            return function(**kwargs)
        except DomainError as error:
            self.fail(f'Valid construction was refused: {error.field}: {error}')

    def refuse_before_flow(self, generator, field, value):
        module = generator.__module__
        with (patch(module+'.mach_from_isentropic',
                    side_effect=AssertionError('Invalid input reached inverse flow')) as flow,
              patch('rocket_propulsion.geometry.moc._build_design',
                    side_effect=AssertionError('Invalid input reached characteristic build')) as build):
            try:
                generator(**(inputs(generator) | {field: value}))
            except DomainError as error:
                self.assertEqual(error.field, field)
            except (TypeError, ValueError, AttributeError, AssertionError, OverflowError) as error:
                self.fail(f'Expected labelled DomainError for {field}; got {type(error).__name__}: {error}')
            else:
                self.fail(f'Accepted invalid {field}: {value!r}')
            flow.assert_not_called()
            build.assert_not_called()

    def test_every_physical_scalar_rejects_malformed_values_before_flow(self):
        names = {generate_nozzle_contour: ('throat_area_m2', 'area_ratio', 'gamma', 'half_angle_deg',
                                         'chamber_radius_ratio', 'bell_length_fraction'),
                 generate_rao_contour: ('throat_area_m2', 'area_ratio', 'gamma', 'length_fraction',
                                       'throat_angle_deg', 'exit_angle_deg', 'chamber_radius_ratio'),
                 generate_moc_nozzle: ('throat_area_m2', 'area_ratio', 'gamma')}
        for generator, parameters in names.items():
            for field in parameters:
                for value in (True, False, '1', None, [], math.nan, math.inf, -math.inf,
                              Fraction(1, 10**500)):
                    with self.subTest(model=generator.__name__, field=field, value=repr(value)[:20]):
                        self.refuse_before_flow(generator, field, value)

    def test_counts_require_bounded_actual_integers(self):
        for generator in GENERATORS:
            field = 'characteristic_count' if generator is generate_moc_nozzle else 'station_count'
            low, high = (5, 80) if generator is generate_moc_nozzle else (41, 401)
            for value in (True, False, None, '41', 41., 41.5, low-1, high+1, 1 << 4096):
                with self.subTest(model=generator.__name__, value=repr(value)[:25]):
                    self.refuse_before_flow(generator, field, value)

    def test_family_name_is_owned_and_bad_names_are_typed(self):
        for value in (None, 1, True, [], '', ' ', 'unknown', '\ud800'):
            with self.subTest(value=repr(value)):
                self.refuse_before_flow(generate_nozzle_contour, 'contour', value)
        result = generate_nozzle_contour(**inputs(generate_nozzle_contour), contour=TextTrap(' Conical '))
        self.assertEqual(result.contour, 'conical')

    def test_moc_correction_is_an_actual_boolean_before_any_build(self):
        for value in ('false', 0, 1, None, []):
            with self.subTest(value=value):
                self.refuse_before_flow(generate_moc_nozzle, 'axisymmetric_correction', value)
        for value in (True, False):
            with self.subTest(value=value):
                result = generate_moc_nozzle(**inputs(generate_moc_nozzle), axisymmetric_correction=value)
                self.assertEqual(result.geometry_mode, 'axisymmetric-area-corrected' if value else 'planar')

    def test_nonphysical_domains_are_labelled_before_flow(self):
        cases = [(g, 'throat_area_m2', 0.) for g in GENERATORS]
        cases += [(g, 'gamma', value) for g in GENERATORS for value in (0., 1.)]
        cases += [(g, 'area_ratio', .99) for g in GENERATORS]
        cases += [(g, 'area_ratio', 1.) for g in (generate_rao_contour, generate_moc_nozzle)]
        cases += [(generate_nozzle_contour, 'half_angle_deg', 45.),
                  (generate_nozzle_contour, 'chamber_radius_ratio', 1.),
                  (generate_nozzle_contour, 'bell_length_fraction', .49),
                  (generate_rao_contour, 'chamber_radius_ratio', 1.),
                  (generate_rao_contour, 'length_fraction', .59),
                  (generate_rao_contour, 'throat_angle_deg', 14.),
                  (generate_rao_contour, 'exit_angle_deg', -1.)]
        for generator, field, value in cases:
            with self.subTest(model=generator.__name__, field=field, value=value):
                self.refuse_before_flow(generator, field, value)

    def test_scalar_subclasses_cannot_change_arithmetic_or_remain_in_results(self):
        for generator in GENERATORS:
            values = inputs(generator)
            raw = {name: ArithmeticTrap(value) if isinstance(value, float) else value
                   for name, value in values.items()}
            with self.subTest(model=generator.__name__):
                actual, expected = generator(**raw), generator(**values)
                self.assertEqual(actual, expected)
                c = actual.contour if generator is generate_moc_nozzle else actual
                for f in fields(c):
                    value = getattr(c, f.name)
                    if isinstance(value, float):
                        self.assertIs(type(value), float)
                json.dumps(asdict(actual), allow_nan=False)

    def test_unity_ratio_is_a_resolved_converging_only_sonic_exit(self):
        for family in ('conical', 'preliminary', 'bell'):
            for count in (41, 121, 401):
                with self.subTest(family=family, count=count):
                    c = self.valid(generate_nozzle_contour, throat_area_m2=.01, area_ratio=1.,
                                   contour=family, station_count=count)
                    self.assertEqual(len(c.stations), count)
                    self.assertEqual(c.diverging_length_m, 0.)
                    self.assertEqual(c.total_length_m, c.converging_length_m)
                    self.assertEqual(c.stations[-1].x_m, 0.)
                    self.assertEqual(c.stations[-1].radius_m, c.throat_radius_m)
                    self.assertEqual(c.stations[-1].mach, 1.)
                    self.assertEqual(c.exit_radius_m, c.throat_radius_m)
                    self.assertEqual(c.half_angle_deg, 0.)
                    self.assertEqual(c.estimated_divergence_efficiency, 1.)
                    self.assertIn('cosine contraction',c.model_fidelity)
                    self.assertIn('no diverging span',c.source)
                    self.assertTrue(all(b.x_m > a.x_m for a, b in zip(c.stations, c.stations[1:])))
                    self.assertTrue(all(p.x_m <= 0. and p.mach <= 1. for p in c.stations))
                    snapshot_thermal_contour(c)

    def test_unity_contour_exports_without_degenerate_faces_and_round_trips(self):
        c = generate_nozzle_contour(throat_area_m2=.01, area_ratio=1., station_count=41)
        export = export_nozzle_contour(c, angular_segments=8)
        self.assertEqual(len(export.csv.strip().splitlines()), 42)
        normal_lines = [s.split()[2:] for s in export.stl.splitlines() if 'facet normal' in s]
        self.assertTrue(normal_lines)
        for row in normal_lines:
            self.assertAlmostEqual(math.hypot(*map(float, row)), 1., delta=3e-15)
        self.assertFalse(export.physical_model_qualified)
        self.assertFalse(export.manufacturing_qualified)

    def test_real_cooling_accepts_unity_with_both_directions_and_closes_energy(self):
        from rocket_propulsion.thermodynamics.real_fluids import reference_fluid
        from tests.test_regen_channel_geometry import CONSTANT, coolant, gas

        c = generate_nozzle_contour(throat_area_m2=.01, area_ratio=1., gamma=gas().gamma,
                                    station_count=41)
        for counter in (True, False):
            with self.subTest(counterflow=counter):
                flow = replace(coolant(), counterflow=counter)
                result = solve_regen_chamber(contour=c, gas=gas(), channels=CONSTANT, coolant=flow)
                fluid = reference_fluid(flow.fluid)
                hi = fluid.state_tp(flow.inlet_temperature_k, flow.inlet_pressure_pa).enthalpy_j_kg
                ho = fluid.state_tp(result.coolant_outlet_temperature_k,
                                    result.coolant_outlet_pressure_pa).enthalpy_j_kg
                expected = flow.mass_flow_kg_s*(ho-hi)
                self.assertGreater(expected, 0.)
                self.assertLess(abs(result.total_heat_load_w-expected), 1e-8*expected)
                self.assertEqual(len(result.stations), 41)

    def test_independent_dimensions_area_mach_and_wall_area(self):
        for generator, extra in ((generate_nozzle_contour, {'contour': 'conical'}),
                                 (generate_nozzle_contour, {}), (generate_rao_contour, {})):
            c = generator(**(inputs(generator) | extra))
            with self.subTest(model=c.contour), mp.workdps(80):
                gamma = mp.mpf(BASE['gamma'])
                self.assertTrue(math.isclose(c.throat_radius_m, math.sqrt(.01/math.pi), rel_tol=3e-15))
                self.assertTrue(math.isclose(c.exit_radius_m/c.throat_radius_m, math.sqrt(10.), rel_tol=3e-15))
                area = math.fsum(math.pi*(a.radius_m+b.radius_m)*math.hypot(
                    b.x_m-a.x_m, b.radius_m-a.radius_m) for a,b in zip(c.stations,c.stations[1:]))
                self.assertTrue(math.isclose(area, c.surface_area_m2, rel_tol=3e-15))
                for p in c.stations:
                    mach = mp.mpf(p.mach)
                    compression = 1+(gamma-1)*mach**2/2
                    ratio = ((2/(gamma+1)*compression)**((gamma+1)/(2*(gamma-1))))/mach
                    with self.subTest(x=p.x_m):
                        self.assertLess(abs(float(ratio)/p.area_ratio-1), 3e-10)
                        self.assertLess(abs(float(1/compression)/p.temperature_total_ratio-1), 3e-15)
                        self.assertLess(abs(float(compression**(-gamma/(gamma-1)))/p.pressure_total_ratio-1), 3e-15)

    def test_similarity_across_large_and_small_representable_geometries(self):
        for generator in GENERATORS:
            base = generator(**inputs(generator))
            base = base.contour if generator is generate_moc_nozzle else base
            for scale in (1e-100, 1e100):
                with self.subTest(model=generator.__name__, scale=scale):
                    other = generator(**(inputs(generator) | {'throat_area_m2': BASE['throat_area_m2']*scale**2}))
                    other = other.contour if generator is generate_moc_nozzle else other
                    self.assertEqual(len(base.stations), len(other.stations))
                    for a,b in zip(base.stations,other.stations,strict=True):
                        self.assertTrue(math.isclose(a.radius_m*scale,b.radius_m,rel_tol=2e-14))
                        self.assertTrue(math.isclose(a.x_m*scale,b.x_m,rel_tol=2e-14,abs_tol=1e-115))
                    self.assertTrue(math.isclose(base.surface_area_m2*scale**2,other.surface_area_m2,rel_tol=2e-14))
                    json.dumps(asdict(other), allow_nan=False)

    def test_http_unity_result_and_exports_are_strict_json(self):
        from rocket_propulsion.api.routes import dispatch_calculation

        result = dispatch_calculation('/api/v1/nozzle-contour', {'throat_area_m2':.01,
            'area_ratio':1., 'station_count':41, 'include_exports':True, 'angular_segments':8})
        self.assertEqual(result['diverging_length_m'],0.)
        self.assertEqual(len(result['stations']),41)
        self.assertIn('exports',result)
        json.dumps(result,allow_nan=False)

    def test_http_unity_comparison_keeps_supported_models_and_explains_exclusions(self):
        from rocket_propulsion.api.routes import dispatch_calculation

        result = self.valid(dispatch_calculation, path='/api/v1/nozzle-contour',
            payload={'throat_area_m2':.01, 'area_ratio':1., 'station_count':41, 'include_comparison':True})
        self.assertEqual({r['contour'] for r in result['comparison']}, {'conical','preliminary'})
        self.assertIn('comparison_exclusions',result)
        self.assertEqual(set(result['comparison_exclusions']), {'rao','moc'})
        for item in result['comparison']:
            self.assertEqual(item['exit_flow_angle_deg'],0.)
            self.assertEqual(item['exit_uniformity_estimate'],1.)
            self.assertEqual(item['estimated_divergence_loss_fraction'],0.)
        json.dumps(result,allow_nan=False)

    def test_unrepresentable_required_exit_area_never_exports_success(self):
        for generator in GENERATORS:
            # 10 * 1e308 exceeds binary64: successful geometry cannot report
            # the required positive exit area. Do not forbid balanced inputs
            # merely because one particular evaluation order is fragile.
            with self.subTest(model=generator.__name__), self.assertRaises(DomainError):
                generator(**(inputs(generator) | {'throat_area_m2':1e308}))

    def test_valid_parameter_endpoints_are_not_rejected(self):
        cases = [(generate_nozzle_contour, {'half_angle_deg':2., 'station_count':401,
                                          'bell_length_fraction':.5}),
                 (generate_nozzle_contour, {'half_angle_deg':math.nextafter(45.,0.),
                                           'bell_length_fraction':1.2}),
                 (generate_rao_contour, {'length_fraction':.6, 'throat_angle_deg':40.,
                                         'exit_angle_deg':0., 'station_count':401}),
                 (generate_rao_contour, {'length_fraction':1., 'throat_angle_deg':40.,
                                         'exit_angle_deg':0.}),
                 (generate_moc_nozzle, {'characteristic_count':80})]
        for generator, change in cases:
            with self.subTest(model=generator.__name__, change=change):
                result = generator(**(inputs(generator) | change))
                json.dumps(asdict(result),allow_nan=False)

    def test_unrepresentable_total_wall_area_is_not_a_successful_infinite_result(self):
        # All cross sections fit: Ac/At = 1.1^2. But r >= rt and ds >= dx
        # imply S >= 2*pi*rt*Lc = 3*At, beyond the largest finite float.
        with self.assertRaises(DomainError) as caught:
            generate_nozzle_contour(throat_area_m2=1e308, area_ratio=1.,
                                    chamber_radius_ratio=1.1, station_count=41)
        self.assertEqual(caught.exception.code,'numeric_range')
        self.assertEqual(caught.exception.field,'surface_area_m2')

    def test_nondefault_conical_lengths_match_independent_tangent_geometry(self):
        with mp.workdps(80):
            for angle in (2.,15.,30.):
                c=generate_nozzle_contour(**inputs(generate_nozzle_contour),
                                          contour='conical',half_angle_deg=angle)
                expected=(mp.mpf(c.exit_radius_m)-mp.mpf(c.throat_radius_m))/mp.tan(mp.pi*angle/180)
                with self.subTest(angle=angle):
                    self.assertTrue(math.isclose(c.diverging_length_m,float(expected),rel_tol=3e-15))


if __name__ == '__main__':
    unittest.main()
