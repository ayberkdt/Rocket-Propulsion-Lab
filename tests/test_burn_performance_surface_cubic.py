"""C1 tensor-product cubic Hermite performance surfaces."""

import random
import unittest
from itertools import product

from rocket_propulsion.core.errors import DomainError
from rocket_propulsion.propulsion.burns import (
    PerformanceAxis,
    PerformanceScalePoint,
    PropagatorPropulsionBridge,
    RectilinearPerformanceSurface,
    TabulatedBurnArtifact,
    TabulatedBurnSegment,
)


def surface_from(function, axes_values, *, interpolation="cubic_hermite"):
    axes = tuple(
        PerformanceAxis(f"x{index}", "1", values) for index, values in enumerate(axes_values)
    )
    points = tuple(
        PerformanceScalePoint(coordinate, *function(*coordinate))
        for coordinate in product(*axes_values)
    )
    return RectilinearPerformanceSurface(
        axes=axes, points=points, source_id="cubic", source_sha256="9" * 64,
        interpolation=interpolation,
    )


def bilinear(x, y):
    return 1.0 + 0.1 * x + 0.2 * y + 0.05 * x * y, 2.0 - 0.1 * x + 0.03 * x * y


def smooth(x, y):
    return 1.0 + 0.3 * x**2 + 0.1 * x * y + 0.2 * y**3, 1.0 + 0.2 * x + 0.1 * y**2


GRID = ((0.0, 0.5, 1.5, 2.0, 3.0), (0.0, 1.0, 1.7, 3.0))


class CubicSurfaceTests(unittest.TestCase):
    def test_nodes_and_bilinear_functions_are_reproduced_exactly(self) -> None:
        surface = surface_from(bilinear, GRID)
        generator = random.Random(1)
        for x, y in product(*GRID):
            value = surface.evaluate({"x0": x, "x1": y})
            self.assertAlmostEqual(value.force_scale, bilinear(x, y)[0], places=13)
        for _ in range(200):
            x, y = generator.uniform(0, 3), generator.uniform(0, 3)
            value = surface.evaluate({"x0": x, "x1": y})
            force, flow = bilinear(x, y)
            self.assertAlmostEqual(value.force_scale, force, places=12)
            self.assertAlmostEqual(value.mass_flow_scale, flow, places=12)
            gradient = dict(value.force_scale_gradient)
            self.assertAlmostEqual(gradient["x0"], 0.1 + 0.05 * y, places=11)
            self.assertAlmostEqual(gradient["x1"], 0.2 + 0.05 * x, places=11)

    def test_gradient_is_continuous_across_grid_lines(self) -> None:
        cubic = surface_from(smooth, GRID)
        linear = surface_from(smooth, GRID, interpolation="multilinear")
        epsilon = 1e-9
        for axis_name, line, other in (("x0", 1.5, ("x1", 2.2)), ("x1", 1.7, ("x0", 0.8))):
            def gradient(surface, value, name=axis_name, other_axis=other):
                conditions = {name: value, other_axis[0]: other_axis[1]}
                return dict(surface.evaluate(conditions).force_scale_gradient)[name]

            self.assertAlmostEqual(gradient(cubic, line - epsilon), gradient(cubic, line + epsilon),
                                   places=6)
            # The multilinear surface really has a kink there.
            self.assertGreater(abs(gradient(linear, line - epsilon) - gradient(linear, line + epsilon)),
                               1e-3)
        self.assertEqual(cubic.grid_line_surfaces(), ())
        self.assertNotEqual(linear.grid_line_surfaces(), ())

    def test_gradients_match_central_finite_differences(self) -> None:
        surface = surface_from(smooth, GRID)
        generator = random.Random(2)
        h = 1e-6
        for _ in range(50):
            x, y = generator.uniform(0.01, 2.99), generator.uniform(0.01, 2.99)
            value = surface.evaluate({"x0": x, "x1": y})
            for name, gradients, field in (
                ("force", value.force_scale_gradient, "force_scale"),
                ("flow", value.mass_flow_scale_gradient, "mass_flow_scale"),
            ):
                analytic = dict(gradients)
                for axis, (dx, dy) in (("x0", (h, 0.0)), ("x1", (0.0, h))):
                    upper = getattr(surface.evaluate({"x0": x + dx, "x1": y + dy}), field)
                    lower = getattr(surface.evaluate({"x0": x - dx, "x1": y - dy}), field)
                    self.assertAlmostEqual(analytic[axis], (upper - lower) / (2 * h), places=6,
                                           msg=f"{name} d/{axis}")

    def test_three_dimensions_and_monotone_lines(self) -> None:
        def trilinear(x, y, z):
            return 1.0 + 0.1 * x + 0.2 * y + 0.3 * z + 0.01 * x * y * z, 1.0 + 0.05 * z

        surface = surface_from(trilinear, ((0.0, 1.0, 2.5), (0.0, 2.0, 3.0), (0.0, 0.5, 1.0)))
        generator = random.Random(3)
        for _ in range(100):
            x, y, z = (generator.uniform(0, 2.5), generator.uniform(0, 3), generator.uniform(0, 1))
            value = surface.evaluate({"x0": x, "x1": y, "x2": z})
            self.assertAlmostEqual(value.force_scale, trilinear(x, y, z)[0], places=12)
        # Along a grid line the cubic is monotone for monotone data (no overshoot).
        line = surface_from(lambda x, y: (1.0 + (x > 1.0) * 0.5 + 0.01 * y, 1.0), GRID)
        previous = None
        for k in range(301):
            value = line.evaluate({"x0": 3.0 * k / 300, "x1": 1.0}).force_scale
            if previous is not None:
                self.assertGreaterEqual(value, previous - 1e-12)
            self.assertLessEqual(value, 1.5 + 0.01 + 1e-12)
            previous = value

    def test_persistence_identity_and_bridge_use(self) -> None:
        cubic = surface_from(smooth, GRID)
        linear = surface_from(smooth, GRID, interpolation="multilinear")
        self.assertEqual(RectilinearPerformanceSurface.from_json(cubic.to_json()), cubic)
        self.assertNotIn("interpolation", linear.to_dict())
        self.assertNotEqual(cubic.surface_hash, linear.surface_hash)
        with self.assertRaises(DomainError):
            surface_from(smooth, GRID, interpolation="spline")
        artifact = TabulatedBurnArtifact(
            segments=(TabulatedBurnSegment(0.0, 5.0, 100.0, 100.0, 100.0, 100.0, 1.0, 1.0),),
            source_id="cubic-bridge", source_sha256="a" * 64,
        )
        bridge = PropagatorPropulsionBridge(artifact, performance_surface=cubic)
        value = bridge.evaluate(
            relative_time_s=1.0, vehicle_mass_kg=10.0, direction=(1.0, 0.0, 0.0),
            operating_conditions={"x0": 1.2, "x1": 0.4},
        )
        scale = cubic.evaluate({"x0": 1.2, "x1": 0.4}).force_scale
        self.assertAlmostEqual(value.axial_thrust_n, 100.0 * scale)
        self.assertEqual(bridge.operating_condition_grid_surfaces(), ())


if __name__ == "__main__":
    unittest.main()
