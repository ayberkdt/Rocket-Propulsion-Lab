"""Global sensitivity analysis against independent code and closed-form indices.

Rank measures are recomputed here with numpy least squares and
``scipy.stats`` (Spearman, average ranks); Sobol indices are compared with
the analytic indices of the Ishigami function (Ishigami and Homma, 1990;
the formulas are written out below), of additive linear models
(``S_i = S_Ti = a_i^2 V_i / sum a_j^2 V_j``) and of the rocket equation
``dv = g0 Isp ln(MR)`` with independent uniform inputs, a product of two
independent factors ``X Y`` for which ``V = E[X^2]E[Y^2] - (E X E Y)^2``,
``S_X = V(X) (E Y)^2 / V``, ``S_Y = (E X)^2 V(Y) / V`` and the interaction
``V(X) V(Y) / V``.
"""

from __future__ import annotations

import math
import random
import unittest

import numpy as np
from scipy import stats

from rocket_propulsion.api.routes import dispatch_calculation
from rocket_propulsion.core.constants import STANDARD_GRAVITY_M_S2
from rocket_propulsion.core.errors import DomainError, InputError, RocketPropulsionError
from rocket_propulsion.studies.global_sensitivity import (
    RANK_R_SQUARED_ADEQUATE,
    rank_sensitivity,
    ranks,
    sobol_indices,
)
from rocket_propulsion.studies.sampling import ParameterRange


def numpy_rank_measures(columns: dict[str, np.ndarray], y: np.ndarray):
    names = list(columns)
    rx = {n: stats.rankdata(columns[n]) for n in names}
    ry = stats.rankdata(y)
    design = np.column_stack([np.ones(len(y))] + [rx[n] for n in names])
    beta, *_ = np.linalg.lstsq(design, ry, rcond=None)
    fitted = design @ beta
    r2 = 1 - np.sum((ry - fitted) ** 2) / np.sum((ry - ry.mean()) ** 2)
    out = {}
    for k, n in enumerate(names):
        others = np.column_stack([np.ones(len(y))] + [rx[m] for m in names if m != n])
        res_x = rx[n] - others @ np.linalg.lstsq(others, rx[n], rcond=None)[0]
        res_y = ry - others @ np.linalg.lstsq(others, ry, rcond=None)[0]
        prcc = np.corrcoef(res_x, res_y)[0, 1]
        srrc = beta[k + 1] * np.std(rx[n]) / np.std(ry)
        out[n] = (stats.spearmanr(columns[n], y)[0], prcc, srrc)
    return out, r2


class RankMeasuresAgainstNumpy(unittest.TestCase):
    def test_random_models_match_independent_least_squares(self) -> None:
        rng = np.random.default_rng(11)
        for trial in range(12):
            n, k = int(rng.integers(30, 400)), int(rng.integers(1, 6))
            columns = {f"x{i}": rng.uniform(-2, 3, n) for i in range(k)}
            weights = rng.normal(size=k)
            y = sum(w * np.tanh(c) for w, c in zip(weights, columns.values(), strict=True)) + 0.3 * rng.normal(size=n)
            if trial % 3 == 0:  # ties in an input and the output
                columns["x0"] = np.round(columns["x0"], 1)
                y = np.round(y, 1)
            expected, r2 = numpy_rank_measures(columns, y)
            result = rank_sensitivity({key: list(v) for key, v in columns.items()}, list(y))
            with self.subTest(trial=trial, n=n, k=k):
                self.assertAlmostEqual(result.rank_r_squared, r2, delta=1e-10)
                for record in result.coefficients:
                    spearman, prcc, srrc = expected[record.parameter]
                    self.assertAlmostEqual(record.spearman, spearman, delta=1e-12)
                    self.assertAlmostEqual(record.prcc, prcc, delta=1e-9)
                    self.assertAlmostEqual(record.srrc, srrc, delta=1e-9)
                prccs = [abs(r.prcc) for r in result.coefficients if r.prcc is not None]
                self.assertEqual(prccs, sorted(prccs, reverse=True))

    def test_ranks_match_scipy_average_ranks(self) -> None:
        rng = random.Random(5)
        for _ in range(50):
            values = [rng.choice([0.0, 1.5, -2.0, 3.25, 7.0]) + rng.choice([0.0, 0.0, 1e-9]) for _ in range(40)]
            self.assertEqual(list(ranks(values)), list(stats.rankdata(values)))

    def test_invariant_to_monotone_transforms_and_row_order(self) -> None:
        rng = random.Random(9)
        x1 = [rng.uniform(0.1, 3.0) for _ in range(200)]
        x2 = [rng.uniform(-1.0, 1.0) for _ in range(200)]
        y = [a * a + 0.5 * b + 0.1 * rng.random() for a, b in zip(x1, x2, strict=True)]
        base = rank_sensitivity({"x1": x1, "x2": x2}, y)
        transformed = rank_sensitivity({"x1": [math.log(a) for a in x1], "x2": [b ** 3 for b in x2]},
                                       [math.exp(v) for v in y])
        self.assertEqual(base.coefficients, transformed.coefficients)
        self.assertEqual(base.rank_r_squared, transformed.rank_r_squared)
        order = list(range(200))
        rng.shuffle(order)
        shuffled = rank_sensitivity({"x1": [x1[i] for i in order], "x2": [x2[i] for i in order]},
                                    [y[i] for i in order])
        for a, b in zip(base.coefficients, shuffled.coefficients, strict=True):
            self.assertEqual(a.parameter, b.parameter)
            for left, right in ((a.spearman, b.spearman), (a.prcc, b.prcc), (a.srrc, b.srrc)):
                assert left is not None and right is not None
                self.assertAlmostEqual(left, right, delta=1e-13)

    def test_independent_inputs_give_srrc_equal_to_rank_correlation(self) -> None:
        # Helton and Davis: RCC and SRRC coincide for orthogonal (independent) inputs
        rng = random.Random(2)
        n = 4000
        columns = {f"x{i}": [rng.random() for _ in range(n)] for i in range(3)}
        y = [3 * a + b - 0.5 * c + 0.2 * rng.random() for a, b, c in zip(*columns.values(), strict=True)]
        for record in rank_sensitivity(columns, y).coefficients:
            assert record.spearman is not None and record.srrc is not None
            self.assertAlmostEqual(record.spearman, record.srrc, delta=0.03)


class RankMeasuresKnowTheirLimits(unittest.TestCase):
    def test_non_monotone_output_is_flagged_and_sobol_still_sees_it(self) -> None:
        rng = random.Random(4)
        x1 = [rng.uniform(-1.0, 1.0) for _ in range(500)]
        x2 = [rng.uniform(-1.0, 1.0) for _ in range(500)]
        y = [a * a + 0.05 * b for a, b in zip(x1, x2, strict=True)]
        result = rank_sensitivity({"x1": x1, "x2": x2}, y)
        self.assertLess(result.rank_r_squared, RANK_R_SQUARED_ADEQUATE)
        self.assertIn("rank_regression_inadequate", {w.code for w in result.warnings})
        x1_record = next(r for r in result.coefficients if r.parameter == "x1")
        assert x1_record.spearman is not None
        self.assertLess(abs(x1_record.spearman), 0.15)
        sobol = sobol_indices(lambda p: p["x1"] ** 2 + 0.05 * p["x2"],
                              [ParameterRange("x1", -1, 1), ParameterRange("x2", -1, 1)],
                              base_samples=4096, bootstrap_resamples=0)
        self.assertGreater(sobol.indices[0].first_order, 0.9)

    def test_output_fixed_by_one_input_leaves_the_other_prcc_undefined(self) -> None:
        rng = random.Random(12)
        a = [rng.random() for _ in range(80)]
        b = [rng.random() for _ in range(80)]
        result = rank_sensitivity({"a": a, "b": b}, [math.exp(3 * v) for v in a])
        record = {r.parameter: r for r in result.coefficients}
        self.assertEqual(record["a"].spearman, 1.0)
        self.assertIsNone(record["b"].prcc)
        self.assertIn("output_determined_by_other_inputs", {w.code for w in result.warnings})

    def test_constant_collinear_and_hostile_inputs(self) -> None:
        rng = random.Random(1)
        x = [rng.random() for _ in range(50)]
        y = [2 * v + 0.1 * rng.random() for v in x]
        constant = rank_sensitivity({"x": x, "c": [1.0] * 50}, y)
        c = next(r for r in constant.coefficients if r.parameter == "c")
        self.assertEqual((c.spearman, c.prcc, c.srrc), (None, None, None))
        self.assertIn("constant_input", {w.code for w in constant.warnings})
        collinear = rank_sensitivity({"x": x, "x_twice": [2 * v for v in x]}, y)
        self.assertTrue(all(r.prcc is None and r.srrc is None for r in collinear.coefficients))
        self.assertIn("collinear_inputs", {w.code for w in collinear.warnings})
        flat = rank_sensitivity({"x": x}, [5.0] * 50)
        self.assertIsNone(flat.rank_r_squared)
        self.assertIn("constant_output", {w.code for w in flat.warnings})
        with self.assertRaises(DomainError):
            rank_sensitivity({"x": x[:3]}, y[:3])
        with self.assertRaises(DomainError):
            rank_sensitivity({"x": [*x[:-1], math.nan]}, y)
        with self.assertRaises(InputError):
            rank_sensitivity({"x": x[:-1]}, y)
        with self.assertRaises(InputError):
            rank_sensitivity({}, y)


def ishigami(p, a=7.0, b=0.1):
    return math.sin(p["x1"]) + a * math.sin(p["x2"]) ** 2 + b * p["x3"] ** 4 * math.sin(p["x1"])


class SobolAgainstClosedForms(unittest.TestCase):
    def test_ishigami(self) -> None:
        a, b = 7.0, 0.1
        pi4 = math.pi ** 4
        v = a * a / 8 + b * pi4 / 5 + b * b * pi4 * pi4 / 18 + 0.5
        v1, v2, v13 = 0.5 * (1 + b * pi4 / 5) ** 2, a * a / 8, b * b * pi4 * pi4 * (1 / 18 - 1 / 50)
        exact = {"x1": (v1 / v, (v1 + v13) / v), "x2": (v2 / v, v2 / v), "x3": (0.0, v13 / v)}
        result = sobol_indices(ishigami, [ParameterRange(n, -math.pi, math.pi) for n in ("x1", "x2", "x3")],
                               base_samples=8192, bootstrap_resamples=300, confidence=0.99, seed=17)
        self.assertEqual(result.model_evaluations, 8192 * 5)
        self.assertAlmostEqual(result.output_standard_deviation ** 2, v, delta=0.03 * v)
        for index in result.indices:
            first, total = exact[index.parameter]
            with self.subTest(parameter=index.parameter):
                self.assertLess(abs(index.first_order - first), 0.04)
                self.assertLess(abs(index.total - total), 0.03)
                lo, hi = index.first_order_interval
                self.assertTrue(lo <= first <= hi, (lo, first, hi))
                lo, hi = index.total_interval
                self.assertTrue(lo <= total <= hi, (lo, total, hi))
                self.assertGreaterEqual(index.total, index.first_order - 0.03)

    def test_additive_model_and_dummy_input(self) -> None:
        weights = {"a": 4.0, "b": -2.0, "c": 0.5}
        widths = {"a": 1.0, "b": 3.0, "c": 2.0}
        variances = {n: weights[n] ** 2 * widths[n] ** 2 / 12 for n in weights}
        total_variance = sum(variances.values())
        parameters = [ParameterRange(n, 0.0, widths[n]) for n in weights] + [ParameterRange("dummy", -1.0, 1.0)]
        result = sobol_indices(lambda p: sum(weights[n] * p[n] for n in weights), parameters,
                               base_samples=4096, bootstrap_resamples=200, seed=3)
        for index in result.indices:
            with self.subTest(parameter=index.parameter):
                if index.parameter == "dummy":
                    # f(A) == f(A_B^dummy) row by row: both estimators vanish exactly
                    self.assertEqual((index.first_order, index.total), (0.0, 0.0))
                    continue
                share = variances[index.parameter] / total_variance
                self.assertLess(abs(index.first_order - share), 0.03)
                self.assertLess(abs(index.total - share), 0.02)
        self.assertLess(abs(result.first_order_sum - 1.0), 0.04)

    def test_affine_output_change_and_reproducibility(self) -> None:
        params = [ParameterRange(n, -math.pi, math.pi) for n in ("x1", "x2", "x3")]
        base = sobol_indices(ishigami, params, base_samples=512, bootstrap_resamples=50, seed=8)
        shifted = sobol_indices(lambda p: 1.0e4 + 25.0 * ishigami(p), params, base_samples=512,
                                bootstrap_resamples=50, seed=8)
        for a, b in zip(base.indices, shifted.indices, strict=True):
            self.assertAlmostEqual(a.first_order, b.first_order, delta=1e-9)
            self.assertAlmostEqual(a.total, b.total, delta=1e-9)
        self.assertEqual(base, sobol_indices(ishigami, params, base_samples=512, bootstrap_resamples=50, seed=8))
        self.assertNotEqual(base.indices, sobol_indices(ishigami, params, base_samples=512,
                                                        bootstrap_resamples=50, seed=9).indices)

    def test_outputs_near_the_double_range_do_not_overflow(self) -> None:
        params = [ParameterRange("x", 1.0, 2.0), ParameterRange("y", 0.0, 1.0)]
        small = sobol_indices(lambda p: p["x"] + 0.3 * p["y"] ** 2, params, base_samples=256, bootstrap_resamples=50)
        for scale in (1e300, -1e300, 1e-300):
            with self.subTest(scale=scale):
                big = sobol_indices(lambda p, k=scale: k * (p["x"] + 0.3 * p["y"] ** 2), params, base_samples=256,
                                    bootstrap_resamples=50)
                for a, b in zip(small.indices, big.indices, strict=True):
                    self.assertAlmostEqual(a.first_order, b.first_order, delta=1e-9)
                    self.assertAlmostEqual(a.total, b.total, delta=1e-9)
                self.assertTrue(math.isfinite(big.output_standard_deviation))
                self.assertAlmostEqual(big.output_standard_deviation / abs(scale), small.output_standard_deviation,
                                       delta=1e-9 * small.output_standard_deviation)

    def test_hostile_arguments(self) -> None:
        params = [ParameterRange("x", 0.0, 1.0)]
        for kwargs in ({"base_samples": 7}, {"base_samples": True}, {"base_samples": 10 ** 7},
                       {"bootstrap_resamples": -1}, {"bootstrap_resamples": 10 ** 6}, {"confidence": 1.0},
                       {"confidence": math.nan}, {"seed": 1.5}):
            with self.subTest(kwargs=kwargs), self.assertRaises(InputError):
                sobol_indices(lambda p: p["x"], params, **kwargs)
        with self.assertRaises(InputError):
            sobol_indices(lambda p: p["x"], [*params, ParameterRange("x", 0.0, 2.0)])
        with self.assertRaises(DomainError):
            sobol_indices(lambda p: 3.0, params)
        with self.assertRaises(DomainError):
            sobol_indices(lambda p: math.nan if p["x"] > 0.5 else 1.0, params)


class BlowdownUncertaintyRoute(unittest.TestCase):
    """Rank measures on the blowdown ensemble, recomputed from the realizations the route returns."""

    def test_route_rank_measures_match_scipy_on_its_own_realizations(self) -> None:
        out = dispatch_calculation("/api/v1/transient/uncertainty", {"sample_count": 64, "seed": 5})
        series = out["series"]["realizations"]
        self.assertEqual(series["point_count"], series["source_point_count"])  # nothing thinned
        columns = series["columns"]
        inputs = {name.removeprefix("input_"): np.array(values, dtype=float)
                  for name, values in columns.items() if name.startswith("input_")}
        self.assertGreaterEqual(len(inputs), 2)
        for record in out["global_sensitivity"]:
            metric = record["output"]
            self.assertNotIn("unavailable", record, record)
            expected, r2 = numpy_rank_measures(inputs, np.array(columns[metric], dtype=float))
            with self.subTest(metric=metric):
                self.assertAlmostEqual(record["rank_r_squared"], r2, delta=1e-9)
                undefined = 0
                for coefficient in record["coefficients"]:
                    spearman, prcc, srrc = expected[coefficient["parameter"]]
                    self.assertAlmostEqual(coefficient["spearman"], spearman, delta=1e-12)
                    self.assertAlmostEqual(coefficient["srrc"], srrc, delta=1e-9)
                    if coefficient["prcc"] is None:
                        # only when another input fixes the output ranks exactly (|Spearman| = 1)
                        undefined += 1
                        self.assertTrue(any(abs(c["spearman"]) == 1.0 for c in record["coefficients"]
                                            if c["parameter"] != coefficient["parameter"]))
                    else:
                        self.assertAlmostEqual(coefficient["prcc"], prcc, delta=1e-9)
                codes = {w["code"] for w in record["warnings"]}
                self.assertEqual(undefined > 0, "output_determined_by_other_inputs" in codes)
                inadequate = "rank_regression_inadequate" in {w["code"] for w in record["warnings"]}
                self.assertEqual(inadequate, record["rank_r_squared"] < RANK_R_SQUARED_ADEQUATE)
        by_metric = {r["output"]: {c["parameter"]: c for c in r["coefficients"]} for r in out["global_sensitivity"]}
        # physics: with the tank and cutoff pressure fixed, the propellant expelled (and so the
        # impulse) is set by the initial tank pressure alone, monotonically
        for metric in ("delivered_impulse_n_s", "consumed_propellant_kg"):
            self.assertEqual(by_metric[metric]["initial_pressure_pa"]["spearman"], 1.0)
            self.assertIsNone(by_metric[metric]["rated_delivered_thrust_n"]["prcc"])


class RocketEquationThroughTheRoute(unittest.TestCase):
    def test_delta_v_indices_match_the_product_formula(self) -> None:
        isp_lo, isp_hi, mr_lo, mr_hi = 250.0, 330.0, 1.5, 5.0
        out = dispatch_calculation("/api/v1/uncertainty", {
            "minimum_specific_impulse_s": isp_lo, "maximum_specific_impulse_s": isp_hi,
            "minimum_mass_ratio": mr_lo, "maximum_mass_ratio": mr_hi, "sobol_base_samples": 8192})
        g0 = STANDARD_GRAVITY_M_S2
        ex, vx = g0 * (isp_lo + isp_hi) / 2, (g0 * (isp_hi - isp_lo)) ** 2 / 12

        def antiderivative(x: float, power: int) -> float:
            return x * math.log(x) - x if power == 1 else x * (math.log(x) ** 2 - 2 * math.log(x) + 2)
        ey = (antiderivative(mr_hi, 1) - antiderivative(mr_lo, 1)) / (mr_hi - mr_lo)
        ey2 = (antiderivative(mr_hi, 2) - antiderivative(mr_lo, 2)) / (mr_hi - mr_lo)
        vy = ey2 - ey * ey
        v = (vx + ex * ex) * ey2 - (ex * ey) ** 2
        exact = {"specific_impulse_s": (vx * ey * ey / v, (vx * ey * ey + vx * vy) / v),
                 "mass_ratio": (ex * ex * vy / v, (ex * ex * vy + vx * vy) / v)}
        sobol = out["global_sensitivity"]["sobol"]
        self.assertAlmostEqual(sobol["output_standard_deviation"] ** 2, v, delta=0.05 * v)
        for index in sobol["indices"]:
            first, total = exact[index["parameter"]]
            with self.subTest(parameter=index["parameter"]):
                self.assertLess(abs(index["first_order"] - first), 0.02)
                self.assertLess(abs(index["total"] - total), 0.01)
                self.assertTrue(index["first_order_interval"][0] <= first <= index["first_order_interval"][1])
                self.assertTrue(index["total_interval"][0] <= total <= index["total_interval"][1])
        rank = out["global_sensitivity"]["rank"]
        self.assertGreater(rank["rank_r_squared"], 0.95)
        self.assertEqual(rank["coefficients"][0]["parameter"], "mass_ratio")

    def test_route_refuses_bad_sobol_sizes(self) -> None:
        base = {"minimum_specific_impulse_s": 280, "maximum_specific_impulse_s": 320,
                "minimum_mass_ratio": 2.0, "maximum_mass_ratio": 4.0}
        for bad in (0, 7, 8193, True, "many", 1.5):
            with self.subTest(value=bad), self.assertRaises(RocketPropulsionError):
                dispatch_calculation("/api/v1/uncertainty", {**base, "sobol_base_samples": bad})


if __name__ == "__main__":
    unittest.main()
