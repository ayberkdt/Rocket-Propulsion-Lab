"""Independent STL topology/orientation, exact interchange and hostile snapshots."""
from collections import Counter
from collections.abc import Mapping
import csv
from dataclasses import asdict, replace
from fractions import Fraction
import io
from math import isfinite, pi, sin, ulp
import unittest
import xml.etree.ElementTree as ET

from rocket_propulsion.core.errors import DomainError
from rocket_propulsion.geometry import NozzleContour, NozzleStation, export_nozzle_contour, generate_rao_contour


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


def cross(a,b):return (a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0])
def vector(a,b):return tuple(Fraction(y)-Fraction(x) for x,y in zip(a,b))
def dot(a,b):return sum((Fraction(x)*Fraction(y) for x,y in zip(a,b)),Fraction(0))


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
        tree=ET.fromstring(self.valid(cylinder(1e-150,1e-150)).svg)
        polygon=tree.find('{http://www.w3.org/2000/svg}polygon')
        points=[tuple(map(float,item.split(','))) for item in polygon.attrib['points'].split()]
        self.assertGreater(max(p[0] for p in points)-min(p[0] for p in points),1000.)
        self.assertGreater(max(p[1] for p in points)-min(p[1] for p in points),400.)

    def test_svg_metadata_is_literal_text_not_executable_xml(self):
        name='prescribed </metadata><script>alert(1)</script><metadata> & source'
        tree=ET.fromstring(self.valid(cylinder(name=name)).svg)
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
        for change in (dict(csv='forged'),dict(svg='forged'),dict(dxf='forged'),dict(stl='forged'),
                       dict(manufacturing_qualified=True),dict(manufacturing_qualified=0),
                       dict(physical_model_qualified=True),dict(serialization_policy='forged')):
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


if __name__=='__main__':unittest.main()
