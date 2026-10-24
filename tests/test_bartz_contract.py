"""Independent signed Bartz equations, surface integration and hostile records."""
import unittest
from collections.abc import Mapping
from dataclasses import asdict, replace
from decimal import Decimal, localcontext
from math import inf, nextafter, pi

from rocket_propulsion.api import routes as api_routes
from rocket_propulsion.core.errors import ConvergenceError, DomainError
from rocket_propulsion.geometry import NozzleContour, NozzleStation, generate_rao_contour
from rocket_propulsion.thermodynamics import calculate_bartz_heat_transfer


def contour():
    # Prescribed states for arithmetic tests; no isentropic flow closure claim.
    points = (NozzleStation(0., 1., pi, 1., 1., 1., 1.),
              NozzleStation(1., 2., 4*pi, 4., 3., 1., 1.))
    return NozzleContour('prescribed-test', 1., 1., 2., 0., 1., 1., pi, 4*pi, 4., 0., points)


def inputs(**changes):
    return dict({'contour': contour(), 'chamber_pressure_pa': 5e6, 'chamber_temperature_k': 3500.,
        'wall_temperature_k': 800., 'characteristic_velocity_m_s': 1700., 'dynamic_viscosity_pa_s': 7e-5,
        'specific_heat_j_kg_k': 3000., 'prandtl_number': .7, 'gamma': 1.22, 'throat_radius_of_curvature_m': .03}, **changes)


def reference(x, precision=180):
    """Decimal equations plus two-node Gaussian integration of linear q and r.

    No production numerical helpers, output values or Simpson/trapezoid weights
    are used. pi is the declared fixed binary64 geometric constant.
    """
    with localcontext() as c:
        c.prec = precision
        d = lambda v: Decimal.from_float(float(v))
        pw = lambda v, n, k: (v.ln()*Decimal(n)/k).exp()
        tc, tw = d(x['chamber_temperature_k']), d(x['wall_temperature_k'])
        diameter = 2*d(x['contour'].throat_radius_m)
        pr = d(x['prandtl_number'])
        recovery = Decimal('.5') if pr == Decimal('.125') else Decimal(1) if pr == 1 else pw(pr, 1, 3)
        base = Decimal('.026')*pw(d(x['dynamic_viscosity_pa_s']), 1, 5)*d(x['specific_heat_j_kg_k'])
        base /= pw(pr, 3, 5)*pw(diameter, 1, 5)
        base *= pw(d(x['chamber_pressure_pa'])/d(x['characteristic_velocity_m_s']), 4, 5)
        base *= pw(diameter/d(x['throat_radius_of_curvature_m']), 1, 10)
        raw = []
        for station in x['contour'].stations:
            compression = 1+(d(x['gamma'])-1)*d(station.mach)**2/2
            sigma = pw((tw/tc*compression+1)/2, -17, 25)*pw(compression, -3, 25)
            h = base*pw(d(station.area_ratio), -9, 10)*sigma
            taw = tc*(1+recovery*(compression-1))/compression
            q = h*(taw-tw)
            raw.append({'x_m': d(station.x_m),'radius_m': d(station.radius_m),'mach': d(station.mach),
                'adiabatic_wall_temperature_k': taw,'gas_side_coefficient_w_m2_k': h,'heat_flux_w_m2': q})
        cumulative = Decimal(0); stations = []
        u0, u1 = Decimal('.5')-Decimal(3).sqrt()/6, Decimal('.5')+Decimal(3).sqrt()/6
        for index, station in enumerate(raw):
            if index:
                before = raw[index-1]
                length = ((station['x_m']-before['x_m'])**2+(station['radius_m']-before['radius_m'])**2).sqrt()
                value = Decimal(0)
                for u in (u0, u1):
                    q = before['heat_flux_w_m2']+u*(station['heat_flux_w_m2']-before['heat_flux_w_m2'])
                    r = before['radius_m']+u*(station['radius_m']-before['radius_m'])
                    value += q*r
                cumulative += d(pi)*length*value
            stations.append(dict(station, cumulative_heat_w=cumulative))
        peak = max(range(len(raw)), key=lambda i: raw[i]['heat_flux_w_m2'])
        return {'stations': tuple({k:float(v) for k,v in row.items()} for row in stations),
            'peak_heat_flux_w_m2': float(raw[peak]['heat_flux_w_m2']),'peak_location_m': float(raw[peak]['x_m']),
            'total_heat_load_w': float(cumulative),'contains_negative_heat_flux': any(row['heat_flux_w_m2'] < 0 for row in raw)}


class BartzAttacks(unittest.TestCase):
    def valid(self, x):
        try: return calculate_bartz_heat_transfer(**x)
        except (DomainError, ConvergenceError, ArithmeticError) as error:
            self.fail(f'Representable prescribed Bartz witness rejected: {error}')

    def check(self, x):
        expected = reference(x)
        self.assertEqual(expected, reference(x, 360), 'Independent reference did not stabilize')
        r = self.valid(x)
        for name in ('peak_heat_flux_w_m2', 'peak_location_m', 'total_heat_load_w', 'contains_negative_heat_flux'):
            self.assertEqual(getattr(r, name), expected[name], name)
        for actual, wanted in zip(r.stations, expected['stations']):
            for name, value in wanted.items(): self.assertEqual(getattr(actual, name), value, name)
        self.assertFalse(r.physical_model_qualified)
        self.assertIsNone(r.quadrature_error_bound_w)
        return r

    def test_boolean_pressure_is_not_one_pascal(self):
        self.valid(inputs(chamber_pressure_pa=1.))
        with self.assertRaises(DomainError): calculate_bartz_heat_transfer(**inputs(chamber_pressure_pa=True))

    def test_all_gas_quantities_reject_boolean_text_none_and_nonfinite(self):
        for name in inputs():
            if name == 'contour': continue
            for value in (True, None, '1', inf, float('nan')):
                with self.subTest(field=name, value=value), self.assertRaises(DomainError):
                    calculate_bartz_heat_transfer(**inputs(**{name:value}))

    def test_unrepresentable_mandatory_fields_fail_typed(self):
        with self.assertRaises(DomainError):
            calculate_bartz_heat_transfer(**inputs(specific_heat_j_kg_k=1e308))

    def test_empty_singleton_and_unordered_contours_fail_typed(self):
        for stations in ((), contour().stations[:1], set(), 'bad'):
            with self.subTest(stations=repr(stations)), self.assertRaises(DomainError):
                calculate_bartz_heat_transfer(**inputs(contour=replace(contour(), stations=stations)))

    def test_contour_rejects_zero_negative_and_nonfinite_consumed_fields(self):
        for field, value in (('area_ratio',0.), ('radius_m',-1.), ('mach',-1.), ('mach',inf), ('x_m',float('nan'))):
            points = (replace(contour().stations[0], **{field:value}),)+contour().stations[1:]
            with self.subTest(field=field, value=value), self.assertRaises(DomainError):
                calculate_bartz_heat_transfer(**inputs(contour=replace(contour(), stations=points)))

    def test_contour_requires_strict_axial_order(self):
        for points in (tuple(reversed(contour().stations)),
                       (contour().stations[0], replace(contour().stations[1], x_m=0.))):
            with self.subTest(points=points), self.assertRaises(DomainError):
                calculate_bartz_heat_transfer(**inputs(contour=replace(contour(), stations=points)))

    def test_contour_quantities_reject_booleans_before_they_become_valid_states(self):
        for name in ('x_m','radius_m','area_ratio','mach'):
            c=contour();points=(replace(c.stations[0],**{name:True}),)+c.stations[1:]
            with self.subTest(field=name),self.assertRaises(DomainError):
                calculate_bartz_heat_transfer(**inputs(contour=replace(c,stations=points)))

    def test_conflicting_area_ratio_and_radius_are_not_independent_geometry(self):
        points = (replace(contour().stations[0], area_ratio=100.),)+contour().stations[1:]
        with self.assertRaises(DomainError):
            calculate_bartz_heat_transfer(**inputs(contour=replace(contour(), stations=points)))

    def test_direct_profile_matches_independent_rational_power_equations(self):
        self.check(inputs())

    def test_negative_heat_flux_is_retained_when_wall_exceeds_local_recovery(self):
        r = self.check(inputs(chamber_temperature_k=4.,wall_temperature_k=3.5,prandtl_number=.125,gamma=1.5))
        self.assertLess(r.stations[-1].heat_flux_w_m2, 0.)
        self.assertTrue(r.contains_negative_heat_flux)

    def test_exact_zero_heat_flux_is_distinguished_from_underflow(self):
        points = (replace(contour().stations[0],mach=2.), replace(contour().stations[1],mach=3.))
        r = self.check(inputs(contour=replace(contour(),stations=points),chamber_temperature_k=4.,
                              wall_temperature_k=3.,prandtl_number=.125,gamma=1.5))
        self.assertEqual(r.stations[0].heat_flux_w_m2, 0.)
        self.assertLess(r.total_heat_load_w, 0.)

    def test_surface_integral_uses_radius_flux_product_not_product_of_means(self):
        x = inputs();r = self.check(x)
        a,b = reference(x)['stations']
        old = .5*(a['heat_flux_w_m2']+b['heat_flux_w_m2'])*pi*(a['radius_m']+b['radius_m'])*2**.5
        self.assertGreater(abs(r.total_heat_load_w-old)/abs(r.total_heat_load_w), .01)

    def test_shape_rescaling_preserves_coefficient_and_scales_surface_area(self):
        x=inputs();r=self.check(x);c=x['contour']
        scaled=replace(c,throat_radius_m=c.throat_radius_m*2,
            stations=tuple(replace(s,x_m=s.x_m*2,radius_m=s.radius_m*2) for s in c.stations))
        # cp compensates Dt**(-.2), and curvature diameter/radius is unchanged.
        y=inputs(contour=scaled,throat_radius_of_curvature_m=.06,specific_heat_j_kg_k=3000.*2**.2)
        s=self.check(y)
        self.assertAlmostEqual(s.total_heat_load_w/r.total_heat_load_w,4.,places=13)

    def test_private_pressure_velocity_ratio_overflow_can_cancel(self):
        x=inputs(chamber_pressure_pa=1e308,characteristic_velocity_m_s=1e-300,
                 dynamic_viscosity_pa_s=1e-300,specific_heat_j_kg_k=1e-300)
        self.check(x)

    def test_private_mach_squared_overflow_does_not_erase_finite_recovery(self):
        points=tuple(replace(s,mach=1e155) for s in contour().stations)
        self.check(inputs(contour=replace(contour(),stations=points)))

    def test_rounded_recovery_temperature_is_not_reused_for_heat_flux(self):
        x=inputs(chamber_temperature_k=4.,wall_temperature_k=3.,prandtl_number=.7,gamma=1.5)
        r=self.check(x)
        self.assertTrue(any(s.heat_flux_w_m2 != s.gas_side_coefficient_w_m2_k*(s.adiabatic_wall_temperature_k-3.) for s in r.stations))

    def test_heat_flux_survives_when_displayed_recovery_equals_wall(self):
        wall=reference(inputs())['stations'][-1]['adiabatic_wall_temperature_k']
        r=self.check(inputs(wall_temperature_k=wall))
        self.assertEqual(r.stations[-1].adiabatic_wall_temperature_k,wall)
        self.assertEqual(r.stations[-1].heat_flux_w_m2,5.040382872567959e-10)

    def test_peak_location_uses_unrounded_flux_when_displayed_values_tie(self):
        c=contour();a=replace(c.stations[0],mach=1.)
        b=replace(c.stations[1],radius_m=1.,area_ratio=1.,mach=nextafter(1.,-inf))
        r=self.check(inputs(contour=replace(c,stations=(a,b))))
        self.assertEqual(r.stations[0].heat_flux_w_m2,r.stations[1].heat_flux_w_m2)
        self.assertEqual(r.peak_location_m,1.)

    def test_private_segment_length_squared_overflow_is_recoverable(self):
        c=contour();points=(c.stations[0],replace(c.stations[1],x_m=1e200))
        self.check(inputs(contour=replace(c,stations=points)))

    def test_forged_complete_result_and_qualification_are_rejected(self):
        r=self.valid(inputs())
        for changes in ({'total_heat_load_w': 1e99},{'energy_balance_residual_w': 1.},
                        {'physical_model_qualified': True},{'physical_model_qualified': 0},{'contains_negative_heat_flux': 1},
                        {'numerical_policy': 'forged'},{'quadrature_error_bound_w': 0.}):
            with self.subTest(changes=changes),self.assertRaises(DomainError):replace(r,**changes)

    def test_forged_station_in_result_is_rejected(self):
        r=self.valid(inputs());points=(replace(r.stations[0],heat_flux_w_m2=0.),)+r.stations[1:]
        with self.assertRaises(DomainError):replace(r,stations=points)

    def test_mapping_round_trip_owns_source_and_station_records(self):
        r=self.valid(inputs());raw=asdict(r);restored=type(r)(**raw)
        self.assertEqual(restored,r)
        raw['inputs']['contour']['stations'][0]['radius_m']=999.
        self.assertEqual(restored.inputs.contour.stations[0].radius_m,1.)

    def test_literal_warning_and_policy_callbacks_cannot_disguise_forgery(self):
        class EqualText(str):
            def __eq__(self,other):return True
            def __ne__(self,other):return False
        r=self.valid(inputs())
        with self.assertRaises(DomainError):replace(r,numerical_policy=EqualText('forged'))
        warning=replace(r.warnings[0],message=EqualText('forged'))
        with self.assertRaises(DomainError):replace(r,warnings=(warning,)+r.warnings[1:])

    def test_warning_iterables_have_a_bounded_cardinality_read(self):
        r=self.valid(inputs())
        def endless():
            for _ in range(6):yield r.warnings[0]
            raise AssertionError('Read beyond thermal warning cardinality rejection')
        with self.assertRaises(DomainError):replace(r,warnings=endless())

    def test_known_mapping_failures_are_typed_and_programmer_errors_propagate(self):
        r=self.valid(inputs())
        class Broken(Mapping):
            def __init__(self,error):self.error=error
            def __len__(self):return 1
            def __getitem__(self,key):raise self.error('callback')
            def __iter__(self):raise self.error('callback')
        with self.assertRaises(DomainError):replace(r,inputs=Broken(KeyError))
        with self.assertRaises(RuntimeError):replace(r,inputs=Broken(RuntimeError))

    def test_zero_residual_is_labelled_ledger_identity_not_measured_energy_closure(self):
        r=self.valid(inputs())
        self.assertEqual(r.residual_kind,'integral_ledger_identity_not_measured_energy_closure')
        self.assertFalse(r.physical_model_qualified)

    def test_existing_rao_contour_remains_supported_without_rounding_clamps(self):
        c=generate_rao_contour(throat_area_m2=.01,area_ratio=25.,station_count=41)
        self.check(inputs(contour=c))


class BartzApiAttacks(unittest.TestCase):
    def payload(self):
        x=inputs();x.pop('contour')
        x.update(throat_area_m2=.01,area_ratio=25.,station_count=41,include_cooling=False)
        return x

    def test_api_retains_replay_inputs_and_explicit_unqualified_scope(self):
        r=api_routes.dispatch_calculation('/api/v1/nozzle-thermal',self.payload())
        self.assertIs(r['physical_model_qualified'],False)
        self.assertIsNone(r['quadrature_error_bound_w'])
        self.assertEqual(r['inputs']['chamber_pressure_pa'],5e6)
        self.assertEqual(r['residual_kind'],'integral_ledger_identity_not_measured_energy_closure')
        self.assertEqual(api_routes.calculation_metadata('/api/v1/nozzle-thermal')['version'],'2.0')

    def test_optional_cooling_does_not_silently_rectify_signed_heat(self):
        x=self.payload();x.update(chamber_temperature_k=4.,wall_temperature_k=3.9,prandtl_number=.125,gamma=1.5)
        r=api_routes.calculate_nozzle_thermal_route(x)
        self.assertLess(r['total_heat_load_w'],0.)
        x['include_cooling']=True
        with self.assertRaises(DomainError) as caught:api_routes.calculate_nozzle_thermal_route(x)
        self.assertEqual(caught.exception.code,'unsupported_nonpositive_bulk_heat')


if __name__=='__main__':unittest.main()
