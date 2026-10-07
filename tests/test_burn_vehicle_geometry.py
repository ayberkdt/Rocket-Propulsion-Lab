"""Thruster force/torque geometry and propellant-dependent mass properties."""

import random
import unittest
from itertools import combinations, product
from math import cos, radians, sin, sqrt

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


def rcs_cluster() -> ThrusterLayout:
    """Twelve-thruster cube-corner layout able to produce any force and torque."""

    mounts = []
    for axis in range(3):
        for sign in (1.0, -1.0):
            for offset in (1.0, -1.0):
                direction = [0.0, 0.0, 0.0]
                direction[axis] = sign
                position = [0.0, 0.0, 0.0]
                position[(axis + 1) % 3] = 0.5 * offset
                mounts.append(
                    ThrusterMount(
                        f"T{axis}{'p' if sign > 0 else 'm'}{'a' if offset > 0 else 'b'}",
                        tuple(position),
                        tuple(direction),
                    )
                )
    return ThrusterLayout(tuple(mounts))


def brute_force_optimum(layout, force, torque, center) -> float:
    """Minimum residual over every feasible passive set (exact NNLS optimum)."""

    matrix = layout.influence_matrix(center)
    target = (*force, *torque)
    best = sqrt(sum(value * value for value in target))
    count = len(layout.mounts)
    for size in range(1, 7):
        for subset in combinations(range(count), size):
            columns = [[matrix[row][j] for row in range(6)] for j in subset]
            gram = [[sum(a * b for a, b in zip(c1, c2, strict=True)) for c2 in columns]
                    for c1 in columns]
            rhs = [sum(a * b for a, b in zip(c, target, strict=True)) for c in columns]
            solution = _solve(gram, rhs)
            if solution is None or any(value < -1e-12 for value in solution):
                continue
            residual = [
                target[row] - sum(columns[k][row] * solution[k] for k in range(size))
                for row in range(6)
            ]
            best = min(best, sqrt(sum(value * value for value in residual)))
    return best


def brute_force_bounded(layout, force, torque, limits) -> float:
    """Exact box-constrained optimum: every variable at 0, at its limit, or free."""

    matrix = layout.influence_matrix((0.0, 0.0, 0.0))
    target = (*force, *torque)
    names = layout.thruster_ids
    count = len(names)
    best = float("inf")
    for states in product(("low", "high", "free"), repeat=count):
        fixed = [limits[names[j]] if states[j] == "high" else 0.0 for j in range(count)]
        free = [j for j in range(count) if states[j] == "free"]
        rhs_target = [
            target[row] - sum(matrix[row][j] * fixed[j] for j in range(count) if j not in free)
            for row in range(6)
        ]
        values = list(fixed)
        if free:
            columns = [[matrix[row][j] for row in range(6)] for j in free]
            gram = [[sum(a * b for a, b in zip(c1, c2, strict=True)) for c2 in columns]
                    for c1 in columns]
            rhs = [sum(a * b for a, b in zip(c, rhs_target, strict=True)) for c in columns]
            solution = _solve(gram, rhs)
            if solution is None:
                continue
            if any(v < -1e-12 or v > limits[names[j]] + 1e-12
                   for v, j in zip(solution, free, strict=True)):
                continue
            for value, j in zip(solution, free, strict=True):
                values[j] = value
        residual = [
            target[row] - sum(matrix[row][j] * values[j] for j in range(count))
            for row in range(6)
        ]
        best = min(best, sqrt(sum(value * value for value in residual)))
    return best


def _solve(matrix, rhs):
    size = len(rhs)
    a = [row[:] + [rhs[i]] for i, row in enumerate(matrix)]
    for column in range(size):
        pivot = max(range(column, size), key=lambda r: abs(a[r][column]))
        if abs(a[pivot][column]) < 1e-10:
            return None
        a[column], a[pivot] = a[pivot], a[column]
        for row in range(size):
            if row != column:
                factor = a[row][column] / a[column][column]
                for k in range(column, size + 1):
                    a[row][k] -= factor * a[column][k]
    return [a[i][size] / a[i][i] for i in range(size)]


class ThrustAllocationTests(unittest.TestCase):
    def test_full_authority_cluster_reproduces_any_demand_exactly(self) -> None:
        layout = rcs_cluster()
        generator = random.Random(8)
        for _ in range(20):
            force = tuple(generator.uniform(-5, 5) for _ in range(3))
            torque = tuple(generator.uniform(-2, 2) for _ in range(3))
            result = layout.allocate(force_n=force, torque_n_m=torque,
                                     center_of_mass_m=(0.0, 0.0, 0.0))
            self.assertTrue(result.exact)
            for actual, wanted in zip((*result.force_n, *result.torque_n_m),
                                      (*force, *torque), strict=True):
                self.assertAlmostEqual(actual, wanted, places=8)
            self.assertTrue(all(value >= 0.0 for _, value in result.thrusts_n))

    def test_matches_brute_force_optimum_on_underactuated_layouts(self) -> None:
        generator = random.Random(21)
        for _ in range(15):
            mounts = tuple(
                ThrusterMount(
                    f"T{i}",
                    tuple(generator.uniform(-1, 1) for _ in range(3)),
                    tuple(generator.uniform(-1, 1) for _ in range(3)),
                )
                for i in range(generator.randint(2, 7))
            )
            layout = ThrusterLayout(mounts)
            force = tuple(generator.uniform(-3, 3) for _ in range(3))
            torque = tuple(generator.uniform(-1, 1) for _ in range(3))
            result = layout.allocate(force_n=force, torque_n_m=torque,
                                     center_of_mass_m=(0.0, 0.0, 0.0))
            optimum = brute_force_optimum(layout, force, torque, (0.0, 0.0, 0.0))
            self.assertLessEqual(result.weighted_residual_norm, optimum + 1e-8)
            self.assertTrue(all(value >= 0.0 for _, value in result.thrusts_n))

    def test_unreachable_demand_returns_best_bounded_approximation(self) -> None:
        aft = ThrusterLayout((ThrusterMount("main", (0.0, 0.0, -1.0), (0.0, 0.0, 1.0)),))
        backwards = aft.allocate(force_n=(0.0, 0.0, -10.0), torque_n_m=(0.0, 0.0, 0.0),
                                 center_of_mass_m=(0.0, 0.0, 0.0))
        self.assertFalse(backwards.exact)
        self.assertEqual(dict(backwards.thrusts_n)["main"], 0.0)
        self.assertAlmostEqual(backwards.weighted_residual_norm, 10.0)
        strong = aft.allocate(force_n=(0.0, 0.0, 50.0), torque_n_m=(0.0, 0.0, 0.0),
                              center_of_mass_m=(0.0, 0.0, 0.0), maximum_thrust_n=20.0)
        self.assertFalse(strong.exact)
        self.assertEqual(strong.thrusters_at_maximum, ("main",))
        self.assertEqual(dict(strong.thrusts_n)["main"], 20.0)
        self.assertAlmostEqual(strong.weighted_residual_norm, 30.0)
        with self.assertRaises(DomainError):
            aft.allocate(force_n=(0, 0, 1), torque_n_m=(0, 0, 0), center_of_mass_m=(0, 0, 0),
                         maximum_thrust_n={"other": 1.0})
        with self.assertRaises(DomainError):
            aft.allocate(force_n=(0, 0, 1), torque_n_m=(0, 0, 0), center_of_mass_m=(0, 0, 0),
                         torque_length_scale_m=0.0)

    def test_bounds_are_imposed_so_a_feasible_exact_split_is_found(self) -> None:
        # Four aft thrusters; an unbounded sparse solution would need > 25 N on
        # two of them, but spreading the load meets the demand within 25 N.
        layout = ThrusterLayout(
            tuple(
                ThrusterMount(name, position, (0.0, 0.0, 1.0))
                for name, position in (
                    ("T1", (0.4, 0.0, -1.2)),
                    ("T2", (-0.4, 0.0, -1.2)),
                    ("T3", (0.0, 0.4, -1.2)),
                    ("T4", (0.0, -0.4, -1.2)),
                )
            )
        )
        result = layout.allocate(force_n=(0.0, 0.0, 60.0), torque_n_m=(0.0, 1.5, 0.0),
                                 center_of_mass_m=(0.0, 0.0, 0.3), maximum_thrust_n=25.0)
        self.assertTrue(result.exact)
        self.assertTrue(all(0.0 <= value <= 25.0 for _, value in result.thrusts_n))
        over = layout.allocate(force_n=(0.0, 0.0, 200.0), torque_n_m=(0.0, 0.0, 0.0),
                               center_of_mass_m=(0.0, 0.0, 0.3), maximum_thrust_n=25.0)
        self.assertFalse(over.exact)
        self.assertEqual(set(over.thrusters_at_maximum), {"T1", "T2", "T3", "T4"})
        self.assertAlmostEqual(over.weighted_residual_norm, 100.0)

    def test_matches_brute_force_bounded_optimum(self) -> None:
        generator = random.Random(34)
        for _ in range(12):
            count = generator.randint(2, 5)
            layout = ThrusterLayout(
                tuple(
                    ThrusterMount(
                        f"T{i}",
                        tuple(generator.uniform(-1, 1) for _ in range(3)),
                        tuple(generator.uniform(-1, 1) for _ in range(3)),
                    )
                    for i in range(count)
                )
            )
            limits = {f"T{i}": generator.uniform(0.5, 4.0) for i in range(count)}
            force = tuple(generator.uniform(-6, 6) for _ in range(3))
            torque = tuple(generator.uniform(-2, 2) for _ in range(3))
            result = layout.allocate(force_n=force, torque_n_m=torque,
                                     center_of_mass_m=(0.0, 0.0, 0.0),
                                     maximum_thrust_n=limits)
            for name, value in result.thrusts_n:
                self.assertGreaterEqual(value, 0.0)
                self.assertLessEqual(value, limits[name] + 1e-12)
            optimum = brute_force_bounded(layout, force, torque, limits)
            self.assertLessEqual(result.weighted_residual_norm, optimum + 1e-8)

    def test_torque_weighting_changes_the_compromise(self) -> None:
        layout = canted_pair(30.0)
        demand = {"force_n": (1.0, 0.0, 0.0), "torque_n_m": (0.0, 0.0, 0.0),
                  "center_of_mass_m": (0.0, 0.0, 0.0)}
        force_first = layout.allocate(**demand, torque_length_scale_m=100.0)
        torque_first = layout.allocate(**demand, torque_length_scale_m=0.01)
        torque_norm = lambda r: sqrt(sum(v * v for v in r.torque_residual_n_m))
        self.assertLess(torque_norm(torque_first), torque_norm(force_first))


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
