"""Tests for shared units, metadata, errors, and numerical evidence."""

import unittest

from rocket_propulsion.core import (
    ConvergenceError,
    InputError,
    from_si,
    get_unit,
    solve_bisect_monotonic,
    to_si,
)


class UnitTests(unittest.TestCase):
    def test_pressure_and_temperature_round_trip(self) -> None:
        self.assertAlmostEqual(to_si(1.0, "bar", quantity="pressure"), 100_000.0)
        temperature_k = to_si(68.0, "degF", quantity="temperature")
        self.assertAlmostEqual(temperature_k, 293.15, places=10)
        self.assertAlmostEqual(
            from_si(temperature_k, "degF", quantity="temperature"), 68.0
        )

    def test_quantity_mismatch_is_structured(self) -> None:
        with self.assertRaises(InputError) as context:
            get_unit("MPa", quantity="temperature")
        self.assertEqual(context.exception.code, "unit_quantity_mismatch")
        self.assertEqual(context.exception.field, "unit")


class NumericalEvidenceTests(unittest.TestCase):
    def test_bisection_returns_convergence_evidence(self) -> None:
        result = solve_bisect_monotonic(lambda value: value**2, 2.0, 0.0, 2.0)
        self.assertAlmostEqual(result.value**2, 2.0, places=11)
        self.assertGreater(result.iterations, 0)
        self.assertLessEqual(result.residual, result.tolerance * 2.0)
        self.assertEqual(result.method, "bisection")

    def test_iteration_limit_raises_typed_error_with_report(self) -> None:
        with self.assertRaises(ConvergenceError) as context:
            solve_bisect_monotonic(
                lambda value: value**2,
                2.0,
                0.0,
                2.0,
                tolerance=1.0e-30,
                max_iterations=1,
            )
        error = context.exception
        self.assertEqual(error.code, "non_convergence")
        self.assertEqual(error.details["iterations"], 1)
        self.assertIn("residual", error.details)


if __name__ == "__main__":
    unittest.main()
