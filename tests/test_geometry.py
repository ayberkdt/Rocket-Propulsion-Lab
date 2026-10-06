"""Tests for parametric nozzle contours and station flow."""

import unittest
from itertools import pairwise
from math import hypot, pi, sqrt

from rocket_propulsion.geometry import (
    export_nozzle_contour,
    generate_moc_nozzle,
    generate_nozzle_contour,
    generate_rao_contour,
)


class NozzleGeometryTests(unittest.TestCase):
    def test_conical_contour_connects_dimensions_and_flow(self) -> None:
        contour = generate_nozzle_contour(
            throat_area_m2=0.01,
            area_ratio=25.0,
            gamma=1.22,
            contour="conical",
            half_angle_deg=15.0,
            station_count=101,
        )
        self.assertEqual(len(contour.stations), 101)
        self.assertAlmostEqual(contour.throat_radius_m, sqrt(0.01 / pi))
        self.assertAlmostEqual(contour.exit_area_m2, 0.25)
        self.assertAlmostEqual(contour.exit_radius_m / contour.throat_radius_m, 5.0)
        self.assertLess(contour.stations[0].mach, 1.0)
        self.assertAlmostEqual(min(point.area_ratio for point in contour.stations), 1.0)
        self.assertGreater(contour.stations[-1].mach, 3.0)

    def test_bell_profile_is_positive_and_reaches_requested_exit(self) -> None:
        contour = generate_nozzle_contour(
            throat_area_m2=0.005,
            area_ratio=40.0,
            contour="bell",
            station_count=121,
        )
        self.assertTrue(all(point.radius_m > 0.0 for point in contour.stations))
        self.assertAlmostEqual(contour.stations[-1].radius_m, contour.exit_radius_m)
        throat_index = min(
            range(len(contour.stations)), key=lambda index: contour.stations[index].area_ratio
        )
        diverging_radii = [point.radius_m for point in contour.stations[throat_index:]]
        self.assertTrue(
            all(right >= left for left, right in pairwise(diverging_radii))
        )

    def test_rao_contour_honors_area_and_is_shorter_than_reference_cone(self) -> None:
        rao = generate_rao_contour(
            throat_area_m2=0.01,
            area_ratio=25.0,
            gamma=1.22,
            length_fraction=0.8,
        )
        cone = generate_nozzle_contour(
            throat_area_m2=0.01,
            area_ratio=25.0,
            gamma=1.22,
            contour="conical",
            half_angle_deg=15.0,
        )
        self.assertEqual(rao.contour, "rao")
        self.assertAlmostEqual(rao.exit_area_m2 / rao.throat_area_m2, 25.0)
        self.assertLess(rao.diverging_length_m, cone.diverging_length_m)
        self.assertGreater(rao.estimated_divergence_efficiency or 0.0, 0.95)

    def test_moc_meets_exit_state_and_boundary_compatibility(self) -> None:
        design = generate_moc_nozzle(
            throat_area_m2=0.01,
            area_ratio=25.0,
            gamma=1.22,
            characteristic_count=20,
        )
        contour = design.contour
        self.assertAlmostEqual(design.exit_mach_computed, design.exit_mach_target, places=8)
        self.assertLess(design.residuals.wall_tangency_max_deg, 1.0e-8)
        self.assertLess(design.residuals.centerline_symmetry_max_deg, 1.0e-12)
        self.assertLess(design.refined_length_change_relative, 0.01)
        self.assertAlmostEqual(contour.exit_area_m2 / contour.throat_area_m2, 25.0)
        self.assertTrue(
            all(
                right.x_m > left.x_m
                for left, right in zip(contour.stations, contour.stations[1:])
            )
        )
        for line in design.characteristic_lines:
            self.assertTrue(all(point.radius_m >= 0.0 for point in line))

    def test_exports_retain_contour_dimensions(self) -> None:
        contour = generate_rao_contour(throat_area_m2=0.01, area_ratio=16.0)
        bundle = export_nozzle_contour(contour, angular_segments=12)
        csv_lines = bundle.csv.strip().splitlines()
        self.assertEqual(len(csv_lines), len(contour.stations) + 1)
        self.assertIn("<svg", bundle.svg)
        self.assertIn("LWPOLYLINE", bundle.dxf)
        self.assertIn(f"{contour.stations[-1].x_m:.12g}", bundle.dxf)
        self.assertTrue(bundle.stl.startswith("solid rocket_propulsion_nozzle"))
        self.assertEqual(bundle.stl.count("facet normal"), 2 * 12 * len(contour.stations))
        vertices = []
        for line in bundle.stl.splitlines():
            tokens = line.strip().split()
            if tokens[:1] == ["vertex"]:
                vertices.append(tuple(float(value) for value in tokens[1:]))
        throat_vertices = [
            point for point in vertices if abs(point[0]) < 1.0e-10 and hypot(point[1], point[2]) > 0.0
        ]
        exit_vertices = [
            point
            for point in vertices
            if abs(point[0] - contour.stations[-1].x_m) < 1.0e-10
            and hypot(point[1], point[2]) > 0.0
        ]
        imported_throat_area = pi * min(hypot(y, z) for _, y, z in throat_vertices) ** 2
        imported_exit_area = pi * max(hypot(y, z) for _, y, z in exit_vertices) ** 2
        self.assertAlmostEqual(imported_throat_area, contour.throat_area_m2, places=9)
        self.assertAlmostEqual(imported_exit_area, contour.exit_area_m2, places=9)


if __name__ == "__main__":
    unittest.main()

