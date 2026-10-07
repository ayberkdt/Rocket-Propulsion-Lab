"""Thruster force/torque geometry and propellant-dependent mass properties."""

import unittest
from math import cos, radians, sin

from rocket_propulsion.core.errors import DomainError
from rocket_propulsion.propulsion.burns import (
    PropellantMassModel,
    PropellantStation,
    ThrusterLayout,
    ThrusterMount,
)


def canted_pair(cant_deg: float = 10.0) -> ThrusterLayout:
    c, s = cos(radians(cant_deg)), sin(radians(cant_deg))
    return ThrusterLayout(
        (
            ThrusterMount("A", (0.5, 0.0, -1.0), (-s, 0.0, c)),
            ThrusterMount("B", (-0.5, 0.0, -1.0), (s, 0.0, c)),
        )
    )


class ThrusterLayoutTests(unittest.TestCase):
    def test_symmetric_cant_cancels_torque_and_reports_cosine_loss(self) -> None:
        layout = canted_pair()
        result = layout.resultant({"A": 10.0, "B": 10.0}, (0.0, 0.0, 0.0))
        self.assertAlmostEqual(result.force_n[0], 0.0)
        self.assertAlmostEqual(result.force_n[2], 20.0 * cos(radians(10.0)))
        for value in result.torque_n_m:
            self.assertAlmostEqual(value, 0.0)
        self.assertAlmostEqual(
            layout.axial_efficiency({"A": 10.0, "B": 10.0}, (0.0, 0.0, 1.0)),
            cos(radians(10.0)),
        )

    def test_single_thruster_torque_and_center_of_mass_derivative(self) -> None:
        layout = canted_pair()
        center = (0.1, 0.2, 0.3)
        result = layout.resultant({"A": 10.0}, center)
        mount = layout.mounts[0]
        arm = tuple(p - c for p, c in zip(mount.position_m, center, strict=True))
        expected = (
            arm[1] * mount.thrust_axis[2] - arm[2] * mount.thrust_axis[1],
            arm[2] * mount.thrust_axis[0] - arm[0] * mount.thrust_axis[2],
            arm[0] * mount.thrust_axis[1] - arm[1] * mount.thrust_axis[0],
        )
        for actual, wanted in zip(result.torque_n_m, (10 * v for v in expected), strict=True):
            self.assertAlmostEqual(actual, wanted)
        step = 1e-6
        for column in range(3):
            plus = list(center)
            minus = list(center)
            plus[column] += step
            minus[column] -= step
            upper = layout.resultant({"A": 10.0}, tuple(plus)).torque_n_m  # type: ignore[arg-type]
            lower = layout.resultant({"A": 10.0}, tuple(minus)).torque_n_m  # type: ignore[arg-type]
            for row in range(3):
                self.assertAlmostEqual(
                    result.torque_wrt_center_of_mass[row][column],
                    (upper[row] - lower[row]) / (2 * step),
                    places=7,
                )

    def test_influence_matrix_reproduces_resultant(self) -> None:
        layout = canted_pair(25.0)
        center = (0.0, 0.1, -0.2)
        matrix = layout.influence_matrix(center)
        thrusts = {"A": 3.0, "B": 7.0}
        result = layout.resultant(thrusts, center)
        combined = [
            sum(matrix[row][col] * thrusts[name] for col, name in enumerate(layout.thruster_ids))
            for row in range(6)
        ]
        for actual, expected in zip(combined, (*result.force_n, *result.torque_n_m), strict=True):
            self.assertAlmostEqual(actual, expected)

    def test_refusals(self) -> None:
        layout = canted_pair()
        with self.assertRaises(DomainError):
            layout.resultant({"C": 1.0}, (0.0, 0.0, 0.0))
        with self.assertRaises(DomainError):
            layout.resultant({"A": -1.0}, (0.0, 0.0, 0.0))
        with self.assertRaises(DomainError):
            layout.axial_efficiency({}, (0.0, 0.0, 1.0))
        with self.assertRaises(DomainError):
            ThrusterMount("Z", (0.0, 0.0, 0.0), (0.0, 0.0, 0.0))
        with self.assertRaises(DomainError):
            ThrusterLayout((ThrusterMount("A", (0, 0, 0), (0, 0, 1)),) * 2)


class PropellantMassModelTests(unittest.TestCase):
    def setUp(self) -> None:
        self.model = PropellantMassModel(
            dry_mass_kg=400.0,
            dry_center_of_mass_m=(0.01, -0.02, 0.6),
            dry_inertia_kg_m2=((120.0, 1.0, -2.0), (1.0, 140.0, 0.5), (-2.0, 0.5, 90.0)),
            stations=(
                PropellantStation("oxidizer_mass_kg", (0.0, 0.3, 0.2)),
                PropellantStation("fuel_mass_kg", (0.1, -0.3, 0.25)),
            ),
        )
        self.states = {"oxidizer_mass_kg": 150.0, "fuel_mass_kg": 60.0}

    def test_center_of_mass_and_parallel_axis_inertia(self) -> None:
        value = self.model.evaluate(self.states)
        self.assertAlmostEqual(value.mass_kg, 610.0)
        expected_center = tuple(
            (400.0 * d + 150.0 * o + 60.0 * f) / 610.0
            for d, o, f in zip((0.01, -0.02, 0.6), (0.0, 0.3, 0.2), (0.1, -0.3, 0.25), strict=True)
        )
        for actual, expected in zip(value.center_of_mass_m, expected_center, strict=True):
            self.assertAlmostEqual(actual, expected)
        # Independent check: sum each body's inertia about the composite centre.
        def shifted(mass, position):
            r = tuple(p - c for p, c in zip(position, value.center_of_mass_m, strict=True))
            rr = sum(x * x for x in r)
            return [[mass * ((rr if i == j else 0.0) - r[i] * r[j]) for j in range(3)] for i in range(3)]
        total = [list(row) for row in self.model.dry_inertia_kg_m2]
        for mass, position in (
            (400.0, (0.01, -0.02, 0.6)),
            (150.0, (0.0, 0.3, 0.2)),
            (60.0, (0.1, -0.3, 0.25)),
        ):
            block = shifted(mass, position)
            for i in range(3):
                for j in range(3):
                    total[i][j] += block[i][j]
        for i in range(3):
            for j in range(3):
                self.assertAlmostEqual(value.inertia_kg_m2[i][j], total[i][j], places=9)

    def test_tank_mass_partials_match_finite_differences(self) -> None:
        nominal = self.model.evaluate(self.states)
        center_partials = dict(nominal.center_of_mass_wrt_state)
        inertia_partials = dict(nominal.inertia_wrt_state)
        step = 1e-4
        for name in self.states:
            plus = dict(self.states, **{name: self.states[name] + step})
            minus = dict(self.states, **{name: self.states[name] - step})
            upper = self.model.evaluate(plus)
            lower = self.model.evaluate(minus)
            for index in range(3):
                self.assertAlmostEqual(
                    center_partials[name][index],
                    (upper.center_of_mass_m[index] - lower.center_of_mass_m[index]) / (2 * step),
                    places=9,
                )
            for i in range(3):
                for j in range(3):
                    self.assertAlmostEqual(
                        inertia_partials[name][i][j],
                        (upper.inertia_kg_m2[i][j] - lower.inertia_kg_m2[i][j]) / (2 * step),
                        places=6,
                    )

    def test_refusals(self) -> None:
        with self.assertRaises(DomainError):
            self.model.evaluate({"oxidizer_mass_kg": 1.0})
        with self.assertRaises(DomainError):
            self.model.evaluate({"oxidizer_mass_kg": -1.0, "fuel_mass_kg": 1.0})
        with self.assertRaises(DomainError):
            PropellantMassModel(1.0, (0, 0, 0), ((1.0, 0, 0), (0, 1.0, 0), (0, 0, 5.0)), ())
        with self.assertRaises(DomainError):
            PropellantMassModel(1.0, (0, 0, 0), ((1.0, 0.2, 0), (0, 1.0, 0), (0, 0, 1.0)), ())


if __name__ == "__main__":
    unittest.main()
