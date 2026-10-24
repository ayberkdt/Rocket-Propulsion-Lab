"""Independent STL topology/orientation, exact interchange and hostile snapshots."""
import csv
import io
import unittest
import xml.etree.ElementTree as ET
from collections import Counter
from collections.abc import Mapping
from dataclasses import asdict, replace
from decimal import Decimal, localcontext
from fractions import Fraction
from math import isfinite, pi, sin, ulp

from rocket_propulsion.api import routes as api_routes
from rocket_propulsion.core.errors import DomainError
from rocket_propulsion.geometry import (
    NozzleContour,
    NozzleStation,
    export_nozzle_contour,
    generate_rao_contour,
)


def cylinder(radius=.125,length=.3,name='prescribed cylinder'):
    area=pi*radius**2
    stations=tuple(NozzleStation(x,radius,area,1.,1.,(5/6)**3.5,5/6) for x in (0.,length))
    return NozzleContour(name,radius,radius,radius,0.,length,length,area,area,1.,0.,stations)


def facets(text):
    """Parse actual serialized numbers, independent of production builders."""
    result=[];normal=None;points=[]
    for line in text.splitlines():
        row=line.split()
        if row[:2]==['facet','normal']:normal=tuple(float(v) for v in row[2:]);points=[]
        elif row[:1]==['vertex']:points.append(tuple(float(v) for v in row[1:]))
        elif row[:1]==['endfacet']:
            assert normal is not None and len(points)==3
            result.append((normal,tuple(points)))
    return result


def cross(a,b):
    # Convert BEFORE multiplication: even the oracle's private cross product
    # can exceed binary64 although every serialized vertex/normal is finite.
    a=tuple(map(Fraction,a));b=tuple(map(Fraction,b))
    return (a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0])
def vector(a,b):return tuple(Fraction(y)-Fraction(x) for x,y in zip(a,b))
def dot(a,b):return sum((Fraction(x)*Fraction(y) for x,y in zip(a,b)),Fraction(0))


def svg_tree(text):
    try:return ET.fromstring(text)
    except ET.ParseError as error:raise AssertionError(f'Exported SVG must be well-formed XML: {error}') from error


def dxf_points(text):
    rows=text.splitlines();assert len(rows)%2==0
    pairs=list(zip(rows[::2],rows[1::2]));xs=[float(v) for k,v in pairs if k=='10'];ys=[float(v) for k,v in pairs if k=='20']
    return tuple(zip(xs,ys)),pairs


class GeometryExportAttacks(unittest.TestCase):
    def valid(self,c,segments=8):
        try:return export_nozzle_contour(c,angular_segments=segments)
        except (DomainError,ArithmeticError) as error:self.fail(f'Representable export witness rejected: {error}')

    def mesh(self,c,segments=8):
        bundle=self.valid(c,segments);mesh=facets(bundle.stl)
        self.assertEqual(len(mesh),2*segments*len(c.stations))
        edges=Counter();volume=Fraction(0)
        for normal,(a,b,d) in mesh:
            n=cross(vector(a,b),vector(a,d))
            self.assertGreater(sum(v*v for v in n),0,'Serialized facet is degenerate')
            self.assertTrue(all(isfinite(v) for v in normal))
            self.assertLessEqual(abs(sum(v*v for v in normal)-1),8*ulp(1.))
            self.assertGreater(dot(n,normal),0,'Normal disagrees with serialized vertex order')
            for left,right in ((a,b),(b,d),(d,a)):edges[(left,right)]+=1
            volume+=dot(a,cross(b,d))/6
        for (a,b),count in edges.items():
            self.assertEqual(count,1,'Repeated directed mesh edge')
            self.assertEqual(edges[(b,a)],1,'Adjacent faces do not use opposite edge directions')
        self.assertGreater(volume,0,'Closed fluid-domain winding must be outward')
        return bundle,mesh,volume

    def test_closed_mesh_has_opposite_directed_edges_and_outward_normals(self):self.mesh(cylinder())

    def test_caps_and_sides_have_consistent_outward_orientation(self):
        c=cylinder();_,mesh,_=self.mesh(c)
        for normal,points in mesh:
            if all(p[0]==0 for p in points):self.assertLess(normal[0],0)
            elif all(p[0]==c.stations[-1].x_m for p in points):self.assertGreater(normal[0],0)
            else:
                radial=tuple(sum(Fraction(p[k]) for p in points) for k in (1,2))
                self.assertGreater(Fraction(normal[1])*radial[0]+Fraction(normal[2])*radial[1],0)

    def test_mesh_volume_matches_regular_polygon_cylinder_at_several_resolutions(self):
        c=cylinder()
        for count in (8,12,32):
            with self.subTest(segments=count):
                _,_,volume=self.mesh(c,count)
                expected=count/2*sin(2*pi/count)*c.throat_radius_m**2*c.total_length_m
                self.assertLess(abs(float(volume)/expected-1),1e-13)

    def test_tiny_facet_area_does_not_turn_valid_normals_into_zero(self):self.mesh(cylinder(1e-150,1e-150))

    def test_huge_private_cross_products_do_not_turn_valid_normals_into_nan(self):self.mesh(cylinder(1e150,1e200))

    def test_csv_round_trips_every_station_value_without_twelve_digit_loss(self):
        c=generate_rao_contour(throat_area_m2=.01,area_ratio=16.,station_count=41)
        rows=list(csv.DictReader(io.StringIO(self.valid(c).csv)))
        self.assertEqual(len(rows),len(c.stations))
        for row,station in zip(rows,c.stations):
            for name,value in row.items():self.assertEqual(float(value),getattr(station,name),name)

    def test_dxf_round_trips_vertices_and_declares_metres_and_entity_subclasses(self):
        c=generate_rao_contour(throat_area_m2=.01,area_ratio=16.,station_count=41)
        points,pairs=dxf_points(self.valid(c).dxf)
        wanted=tuple((s.x_m,s.radius_m) for s in c.stations)+tuple((s.x_m,-s.radius_m) for s in reversed(c.stations))
        self.assertEqual(points,wanted)
        self.assertIn(('100','AcDbEntity'),pairs);self.assertIn(('100','AcDbPolyline'),pairs)
        index=pairs.index(('9','$INSUNITS'));self.assertEqual(pairs[index+1],('70','6'))

    def test_stl_ring_axis_vertices_preserve_source_coordinates_exactly(self):
        c=generate_rao_contour(throat_area_m2=.01,area_ratio=16.,station_count=41)
        points={p for _,triangle in facets(self.valid(c).stl) for p in triangle}
        for s in c.stations:self.assertIn((s.x_m,s.radius_m,0.),points)

    def test_svg_tiny_geometry_uses_viewport_without_absolute_length_floor(self):
        tree=svg_tree(self.valid(cylinder(1e-150,1e-150)).svg)
        polygon=tree.find('{http://www.w3.org/2000/svg}polygon')
        points=[tuple(map(float,item.split(','))) for item in polygon.attrib['points'].split()]
        self.assertGreater(max(p[0] for p in points)-min(p[0] for p in points),1000.)
        self.assertGreater(max(p[1] for p in points)-min(p[1] for p in points),400.)

    def test_svg_metadata_is_literal_text_not_executable_xml(self):
        name='prescribed </metadata><script>alert(1)</script><metadata> & source'
        tree=svg_tree(self.valid(cylinder(name=name)).svg)
        self.assertFalse(any(e.tag.rsplit('}',1)[-1]=='script' for e in tree.iter()))
        self.assertIn(name,tree.find('{http://www.w3.org/2000/svg}metadata').text)

    def test_angular_resolution_rejects_bool_float_text_and_none_typed(self):
        for value in (True,8.,'8',None,7,129):
            with self.subTest(value=value),self.assertRaises(DomainError):export_nozzle_contour(cylinder(),angular_segments=value)

    def test_empty_singleton_and_unordered_station_sources_are_rejected(self):
        c=cylinder()
        for points in ((),c.stations[:1],set(),None,'bad'):
            with self.subTest(points=repr(points)),self.assertRaises(DomainError):export_nozzle_contour(replace(c,stations=points))

    def test_axial_order_and_positive_radii_are_required_before_serialization(self):
        c=cylinder()
        for points in (tuple(reversed(c.stations)),(c.stations[0],replace(c.stations[1],x_m=0.)),
                       (replace(c.stations[0],radius_m=-1.),c.stations[1])):
            with self.subTest(points=points),self.assertRaises(DomainError):export_nozzle_contour(replace(c,stations=points))

    def test_csv_quantities_cannot_hide_nonfinite_values_outside_xy(self):
        c=cylinder()
        for name in ('area_m2','area_ratio','mach','pressure_total_ratio','temperature_total_ratio'):
            points=(replace(c.stations[0],**{name:float('nan')}),)+c.stations[1:]
            with self.subTest(field=name),self.assertRaises(DomainError):export_nozzle_contour(replace(c,stations=points))

    def test_inconsistent_station_area_and_radius_are_rejected(self):
        c=cylinder();points=(replace(c.stations[0],area_m2=c.stations[0].area_m2*2),)+c.stations[1:]
        with self.assertRaises(DomainError):export_nozzle_contour(replace(c,stations=points))

    def test_input_snapshot_owns_mutable_station_collection(self):
        c=cylinder();points=list(c.stations);r=self.valid(replace(c,stations=points))
        points[0]=replace(points[0],radius_m=999.)
        self.assertEqual(r.inputs.stations[0].radius_m,c.stations[0].radius_m)

    def test_bundle_replay_rejects_forged_text_and_qualification_flags(self):
        r=self.valid(cylinder())
        for change in ({'csv': 'forged'},{'svg': 'forged'},{'dxf': 'forged'},{'stl': 'forged'},
                       {'manufacturing_qualified': True},{'manufacturing_qualified': 0},
                       {'physical_model_qualified': True},{'serialization_policy': 'forged'}):
            with self.subTest(change=change),self.assertRaises(DomainError):replace(r,**change)

    def test_mapping_round_trip_replays_and_snapshots_nested_station_records(self):
        r=self.valid(cylinder());raw=asdict(r);restored=type(r)(**raw)
        self.assertEqual(restored,r)
        raw['inputs']['stations'][0]['x_m']=999.
        self.assertEqual(restored.inputs.stations[0].x_m,0.)

    def test_text_comparison_callbacks_cannot_disguise_forged_files(self):
        class EqualText(str):
            def __eq__(self,other):return True
            def __ne__(self,other):return False
        r=self.valid(cylinder())
        with self.assertRaises(DomainError):replace(r,stl=EqualText('forged'))

    def test_malformed_mapping_is_typed_and_unrelated_callback_failure_propagates(self):
        class Broken(Mapping):
            def __init__(self,error):self.error=error
            def __len__(self):return 1
            def __getitem__(self,key):raise self.error('callback')
            def __iter__(self):raise self.error('callback')
        r=self.valid(cylinder())
        with self.assertRaises(DomainError):replace(r,inputs=Broken(KeyError))
        with self.assertRaises(RuntimeError):replace(r,inputs=Broken(RuntimeError))

    def test_mesh_role_and_units_do_not_claim_a_manufacturable_nozzle_wall(self):
        r=self.valid(cylinder())
        self.assertEqual(r.mesh_role,'closed_fluid_domain_visualization')
        self.assertEqual(r.coordinate_unit,'m')
        self.assertFalse(r.manufacturing_qualified);self.assertFalse(r.physical_model_qualified)

    def test_converging_and_diverging_frusta_reproduce_independent_volume(self):
        base=cylinder();radii=(.375,.125,.25);xs=(-.2,0.,.3)
        stations=tuple(NozzleStation(x,r,pi*r*r,(r/.125)**2,1.,.5,.8) for x,r in zip(xs,radii))
        c=replace(base,stations=stations,area_ratio=4.)
        for count in (9,11,128):
            with self.subTest(count=count):
                _,_,volume=self.mesh(c,count)
                expected=sum((xs[i+1]-xs[i])*(radii[i]**2+radii[i]*radii[i+1]+radii[i+1]**2)
                             for i in range(2))*count*sin(2*pi/count)/6
                self.assertLess(abs(float(volume)/expected-1),1e-13)

    def test_normals_match_independent_high_precision_vertex_cross_products(self):
        c=generate_rao_contour(throat_area_m2=.01,area_ratio=16.,station_count=41)
        mesh=facets(self.valid(c).stl)
        for _,(a,b,d) in mesh[::31]:
            exact=cross(vector(a,b),vector(a,d));expected=[]
            for precision in (110,180):
                with localcontext() as ctx:
                    ctx.prec=precision
                    values=[Decimal(v.numerator)/Decimal(v.denominator) for v in exact]
                    magnitude=sum(v*v for v in values).sqrt()
                    expected.append(tuple(float(v/magnitude) for v in values))
            self.assertEqual(expected[0],expected[1])
            actual=next(normal for normal,points in mesh if points==(a,b,d))
            self.assertEqual(actual,expected[1])

    def test_svg_rejects_xml_forbidden_characters_and_preserves_valid_literal_names(self):
        for name in ('bad\x00name','bad\x01name','bad\ufffename','bad\uffffname'):
            with self.subTest(name=repr(name)),self.assertRaises(DomainError):export_nozzle_contour(cylinder(name=name))
        name='A & B < C > D\r\n\t nozzle'
        tree=svg_tree(self.valid(cylinder(name=name)).svg)
        self.assertIn(name,tree.find('{http://www.w3.org/2000/svg}metadata').text)

    def test_export_boundary_types_all_consumed_station_scalars_and_contour(self):
        for c in (None,{},'contour'):
            with self.subTest(contour=c),self.assertRaises(DomainError):export_nozzle_contour(c)
        c=cylinder()
        for name in asdict(c.stations[0]):
            for value in (True,None,'1',float('inf')):
                with self.subTest(field=name,value=value),self.assertRaises(DomainError):
                    export_nozzle_contour(replace(c,stations=(replace(c.stations[0],**{name:value}),c.stations[1])))

    def test_stored_unit_role_and_nominal_geometry_cannot_be_tampered(self):
        r=self.valid(cylinder())
        for change in ({'coordinate_unit': 'mm'},{'mesh_role': 'manufacturable_wall'}):
            with self.subTest(change=change),self.assertRaises(DomainError):replace(r,**change)
        for change in ({'throat_area_m2': 2*r.inputs.throat_area_m2},{'area_ratio': 4.}):
            with self.subTest(change=change),self.assertRaises(DomainError):replace(r.inputs,**change)

    def test_large_opposite_axial_coordinates_do_not_overflow_svg_span(self):
        c=cylinder();points=(replace(c.stations[0],x_m=-1e308),replace(c.stations[1],x_m=1e308))
        r=self.valid(replace(c,stations=points))
        tree=svg_tree(r.svg)
        points=[tuple(map(float,p.split(','))) for p in tree.find('{http://www.w3.org/2000/svg}polygon').attrib['points'].split()]
        self.assertEqual((min(p[0] for p in points),max(p[0] for p in points)),(30.,1170.))

    def test_unrepresentable_mandatory_normal_component_fails_typed(self):
        c=cylinder(1e150,1e-200);end=replace(c.stations[1],radius_m=2e150,area_m2=pi*(2e150)**2,area_ratio=4.)
        with self.assertRaises(DomainError):export_nozzle_contour(replace(c,stations=(c.stations[0],end),area_ratio=4.))

    def test_boolean_mach_has_valid_numeric_positive_control(self):
        c=cylinder();self.valid(c)
        with self.assertRaises(DomainError):
            export_nozzle_contour(replace(c,stations=(replace(c.stations[0],mach=True),c.stations[1])))

    def test_float_angular_count_does_not_impersonate_an_actual_integer(self):
        self.valid(cylinder(),8);self.valid(cylinder(),128)
        with self.assertRaises(DomainError):export_nozzle_contour(cylinder(),angular_segments=8.)

    def test_reversed_axial_order_is_refused_before_rendering(self):
        c=cylinder();self.valid(c)
        with self.assertRaises(DomainError):export_nozzle_contour(replace(c,stations=tuple(reversed(c.stations))))


class GeometryExportApiAttacks(unittest.TestCase):
    def test_all_geometry_modes_return_owned_export_policy_and_matching_coordinates(self):
        for mode in ('conical','preliminary','rao','moc'):
            x={'throat_area_m2': .01,'area_ratio': 16.,'contour': mode,'station_count': 41,
                   'characteristic_count': 10,'include_exports': True,'angular_segments': 8}
            for invoke in (lambda*, x=x:api_routes.dispatch_calculation('/api/v1/nozzle-contour',x),
                           lambda*, x=x:api_routes.calculate_nozzle_contour_route(x)):
                with self.subTest(mode=mode,path=invoke):
                    r=invoke();export=r['exports']
                    self.assertIs(export['manufacturing_qualified'],False)
                    self.assertIs(export['physical_model_qualified'],False)
                    self.assertEqual(export['coordinate_unit'],'m')
                    self.assertEqual(export['mesh_role'],'closed_fluid_domain_visualization')
                    self.assertEqual(export['inputs']['angular_segments'],8)
                    points,_=dxf_points(export['dxf'])
                    self.assertEqual(points[:len(r['stations'])],tuple((s['x_m'],s['radius_m']) for s in r['stations']))


if __name__=='__main__':unittest.main()
