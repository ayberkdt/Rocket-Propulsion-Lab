"""Global sensitivity analysis: rank-based measures and variance-based Sobol indices.

Two complementary families are provided:

* **Rank-based measures** from an existing sample [1]_: the rank correlation
  coefficient (Spearman), the partial rank correlation coefficient (PRCC)
  and the standardized rank regression coefficient (SRRC), with the
  coefficient of determination ``R^2`` of the rank regression. They capture
  any monotone relationship, are unchanged by monotone transformations of
  inputs or output, and PRCC removes the linear (in ranks) effect of the
  other inputs. A low rank ``R^2`` (below 0.7 here) says the monotone
  additive picture is inadequate and the coefficients must not be read as
  shares of the variance.
* **Sobol indices** [2]_: the first-order index ``S_i = V[E(Y|X_i)]/V(Y)``
  and the total index ``S_Ti = E[V(Y|X_~i)]/V(Y)`` by the Saltelli design
  (matrices ``A``, ``B`` and ``A_B^(i)``, ``N (k + 2)`` model runs) with the
  first-order estimator ``mean(f(B) (f(A_B^i) - f(A)))`` and the Jansen
  total estimator ``mean((f(A) - f(A_B^i))^2)/2``, both normalized by the
  variance of the pooled ``f(A)`` and ``f(B)`` outputs. These apportion
  variance with interactions included and need no monotonicity;
  ``S_Ti - S_i`` measures the interaction share. Percentile bootstrap
  intervals over the ``N`` base rows quantify the sampling error.

Inputs are independent uniform marginals on declared bounds. Correlated
inputs, dependent-input Sobol indices and quasi-random designs are not
provided.

References
----------
.. [1] SAND2001-0417. J. C. Helton and F. J. Davis, "Latin Hypercube
       Sampling and the Propagation of Uncertainty in Analyses of Complex
       Systems", SAND2001-0417, Sandia National Laboratories, 2002, Sect. 6
       (RCC, PRCC, SRRC, rank-regression R^2), PDF pp. 43-44;
       https://doi.org/10.2172/806696. Journal version: Reliability
       Engineering & System Safety 81(1), 23-69, 2003.
.. [2] SALTELLI-2010-CPC. A. Saltelli, P. Annoni, I. Azzini, F. Campolongo,
       M. Ratto and S. Tarantola, "Variance based sensitivity analysis of
       model output. Design and estimator for the total sensitivity index",
       Computer Physics Communications 181(2), 259-270, 2010,
       doi:10.1016/j.cpc.2009.09.018.
"""

from __future__ import annotations

from collections.abc import Callable, Mapping, Sequence
from dataclasses import dataclass
from math import fsum, isfinite, sqrt
from random import Random

from rocket_propulsion.core.errors import DomainError, InputError
from rocket_propulsion.core.metadata import WarningMessage

from .sampling import ParameterRange

__all__ = [
    "RANK_R_SQUARED_ADEQUATE",
    "RankSensitivity",
    "RankSensitivityResult",
    "SobolIndex",
    "SobolResult",
    "rank_sensitivity",
    "ranks",
    "sobol_indices",
]

#: Rank-regression R^2 below which the monotone additive reading is flagged (library choice).
RANK_R_SQUARED_ADEQUATE = 0.7
#: Upper bound on base sample size and bootstrap resamples (library choice: run time).
MAXIMUM_SOBOL_BASE_SAMPLES = 65_536
MAXIMUM_BOOTSTRAP_RESAMPLES = 2_000
#: Relative pivot below which a regression column counts as collinear.
_COLLINEAR = 1.0e-10


@dataclass(frozen=True, slots=True)
class RankSensitivity:
    """Rank-based importance of one input.

    Attributes
    ----------
    parameter : str
        Input name.
    spearman : float or None
        Rank correlation coefficient with the output.
    prcc : float or None
        Partial rank correlation coefficient.
    srrc : float or None
        Standardized rank regression coefficient.
    """

    parameter: str
    spearman: float | None
    prcc: float | None
    srrc: float | None


@dataclass(frozen=True, slots=True)
class RankSensitivityResult:
    """Rank-based sensitivity of one output to every input.

    Attributes
    ----------
    output : str
        Output name.
    sample_count : int
        Number of sample rows.
    rank_r_squared : float or None
        ``R^2`` of the regression of the output ranks on all input ranks.
    coefficients : tuple of RankSensitivity
        One record per input, ordered by decreasing ``|prcc|``.
    warnings : tuple of WarningMessage
        Constant inputs, collinearity and inadequate rank ``R^2``.
    """

    output: str
    sample_count: int
    rank_r_squared: float | None
    coefficients: tuple[RankSensitivity, ...]
    warnings: tuple[WarningMessage, ...]


@dataclass(frozen=True, slots=True)
class SobolIndex:
    """First-order and total Sobol index of one input with bootstrap intervals.

    Attributes
    ----------
    parameter : str
        Input name.
    first_order, total : float
        Point estimates of ``S_i`` and ``S_Ti``.
    first_order_interval, total_interval : tuple of float
        Percentile bootstrap interval at ``confidence``.
    """

    parameter: str
    first_order: float
    total: float
    first_order_interval: tuple[float, float]
    total_interval: tuple[float, float]


@dataclass(frozen=True, slots=True)
class SobolResult:
    """Variance-based sensitivity of a scalar model.

    Attributes
    ----------
    base_sample_count : int
        Rows ``N`` of each base matrix.
    model_evaluations : int
        ``N (k + 2)``.
    output_mean, output_standard_deviation : float
        Mean and standard deviation of the pooled ``f(A)``, ``f(B)`` outputs.
    confidence : float
        Bootstrap interval level.
    indices : tuple of SobolIndex
        One record per input, in input order.
    first_order_sum : float
        ``sum S_i``; ``1 - sum S_i`` estimates the interaction share.
    warnings : tuple of WarningMessage
        Sampling-error and estimator warnings.
    """

    base_sample_count: int
    model_evaluations: int
    output_mean: float
    output_standard_deviation: float
    confidence: float
    indices: tuple[SobolIndex, ...]
    first_order_sum: float
    warnings: tuple[WarningMessage, ...]


def ranks(values: Sequence[float]) -> tuple[float, ...]:
    """Return 1-based ranks with ties given their average rank.

    Parameters
    ----------
    values : sequence of float
        Finite values.

    Returns
    -------
    tuple of float
        Rank of each value in input order.
    """

    order = sorted(range(len(values)), key=lambda index: values[index])
    result = [0.0] * len(values)
    start = 0
    while start < len(order):
        stop = start
        while stop + 1 < len(order) and values[order[stop + 1]] == values[order[start]]:
            stop += 1
        average = 0.5 * (start + stop) + 1.0
        for position in range(start, stop + 1):
            result[order[position]] = average
        start = stop + 1
    return tuple(result)


def _centred(values: Sequence[float]) -> list[float]:
    mean = fsum(values) / len(values)
    return [value - mean for value in values]


def _correlation(left: Sequence[float], right: Sequence[float]) -> float | None:
    a, b = _centred(left), _centred(right)
    saa, sbb = fsum(x * x for x in a), fsum(y * y for y in b)
    if saa <= 0.0 or sbb <= 0.0:
        return None
    value = fsum(x * y for x, y in zip(a, b, strict=True)) / sqrt(saa * sbb)
    return max(-1.0, min(1.0, value))


def _orthonormal_basis(columns: Sequence[Sequence[float]]) -> list[list[float]]:
    """Modified Gram-Schmidt with one reorthogonalization; collinear columns are dropped."""
    basis: list[list[float]] = []
    for column in columns:
        vector = list(column)
        norm0 = sqrt(fsum(x * x for x in vector))
        for _ in range(2):
            for q in basis:
                projection = fsum(x * y for x, y in zip(q, vector, strict=True))
                vector = [x - projection * y for x, y in zip(vector, q, strict=True)]
        norm = sqrt(fsum(x * x for x in vector))
        if norm0 > 0.0 and norm > _COLLINEAR * norm0:
            basis.append([x / norm for x in vector])
    return basis


def _residual(vector: Sequence[float], basis: Sequence[Sequence[float]]) -> list[float]:
    result = list(vector)
    for _ in range(2):
        for q in basis:
            projection = fsum(x * y for x, y in zip(q, result, strict=True))
            result = [x - projection * y for x, y in zip(result, q, strict=True)]
    return result


def _solve_upper(matrix: list[list[float]], rhs: list[float]) -> list[float]:
    size = len(rhs)
    solution = [0.0] * size
    for row in range(size - 1, -1, -1):
        solution[row] = (rhs[row] - fsum(matrix[row][col] * solution[col] for col in range(row + 1, size))) \
            / matrix[row][row]
    return solution


def rank_sensitivity(inputs: Mapping[str, Sequence[float]], output: Sequence[float], *,
                     output_name: str = "output") -> RankSensitivityResult:
    """Return Spearman, PRCC and SRRC of one output to every input [1]_.

    Parameters
    ----------
    inputs : mapping of str to sequence of float
        Sample columns, one per input, all of the output's length.
    output : sequence of float
        Model output for each sample row.
    output_name : str, default "output"
        Label carried into the result.

    Returns
    -------
    RankSensitivityResult
        Coefficients ordered by decreasing ``|prcc|`` (undefined last).

    Raises
    ------
    InputError
        If there are no inputs or the columns differ in length.
    DomainError
        If a value is not finite, or there are fewer than ``k + 3`` rows
        for ``k`` inputs (PRCC needs residual degrees of freedom).

    Notes
    -----
    All columns are replaced by average ranks. PRCC of input ``j`` is the
    correlation of the residuals of ``rank(x_j)`` and ``rank(y)`` after
    least-squares regression on the other input ranks and a constant;
    SRRC is the coefficient of the regression of ``rank(y)`` on all input
    ranks scaled by ``sd(rank x_j)/sd(rank y)``. Regressions use modified
    Gram-Schmidt with reorthogonalization. A constant input gets ``None``
    coefficients; collinear inputs give ``None`` PRCC and SRRC.

    References
    ----------
    .. [1] SAND2001-0417. Helton and Davis, Sect. 6, PDF pp. 43-44.
    """

    if not inputs:
        raise InputError("At least one input column is required.", field="inputs")
    names = list(inputs)
    count = len(output)
    columns = []
    for name in names:
        column = [float(value) for value in inputs[name]]
        if len(column) != count:
            raise InputError(f"Input {name!r} has {len(column)} rows; the output has {count}.", field=name)
        columns.append(column)
    values = [float(value) for value in output]
    if not all(isfinite(value) for column in (*columns, values) for value in column):
        raise DomainError("Sensitivity samples must be finite.", field="output")
    if count < len(names) + 3:
        raise DomainError(f"Rank sensitivity of {len(names)} inputs needs at least {len(names) + 3} rows.",
                          field="output")
    rank_inputs = [_centred(ranks(column)) for column in columns]
    rank_output = _centred(ranks(values))
    warnings: list[WarningMessage] = []
    constant = {name for name, column in zip(names, rank_inputs, strict=True) if not any(column)}
    for name in sorted(constant):
        warnings.append(WarningMessage(code="constant_input", field=name,
                                       message=f"Input {name!r} is constant; its sensitivities are undefined."))
    output_constant = not any(rank_output)
    if output_constant:
        warnings.append(WarningMessage(code="constant_output",
                                       message="The output is constant; every sensitivity is undefined."))
    active = [index for index, name in enumerate(names) if name not in constant]

    # SRRC and R^2: regression of centred output ranks on all centred input ranks (QR by MGS).
    basis = _orthonormal_basis([rank_inputs[index] for index in active])
    full_rank = len(basis) == len(active)
    if not full_rank:
        warnings.append(WarningMessage(code="collinear_inputs",
                                       message="Input ranks are collinear; PRCC and SRRC are undefined."))
    srrc: dict[str, float | None] = {name: None for name in names}
    r_squared: float | None = None
    total = fsum(y * y for y in rank_output)
    if full_rank and active and not output_constant:
        residual = _residual(rank_output, basis)
        r_squared = max(0.0, min(1.0, 1.0 - fsum(r * r for r in residual) / total))
        # R and Q^T y from the basis: R[i][j] = q_i . a_j for j >= i
        r_matrix = [[fsum(q * a for q, a in zip(basis[i], rank_inputs[active[j]], strict=True)) if j >= i else 0.0
                     for j in range(len(active))] for i in range(len(active))]
        qty = [fsum(q * y for q, y in zip(basis[i], rank_output, strict=True)) for i in range(len(active))]
        beta = _solve_upper(r_matrix, qty)
        sd_y = sqrt(total)
        for coefficient, index in zip(beta, active, strict=True):
            sd_x = sqrt(fsum(x * x for x in rank_inputs[index]))
            srrc[names[index]] = coefficient * sd_x / sd_y
        if r_squared < RANK_R_SQUARED_ADEQUATE:
            warnings.append(WarningMessage(
                code="rank_regression_inadequate",
                message=(f"Rank-regression R^2 = {r_squared:.3f} < {RANK_R_SQUARED_ADEQUATE}: the output is not "
                         "a near-additive monotone function of the inputs; use Sobol indices.")))

    records = []
    for index, name in enumerate(names):
        spearman = None if name in constant or output_constant else _correlation(rank_inputs[index], rank_output)
        prcc = None
        if name not in constant and not output_constant and full_rank:
            others = _orthonormal_basis([rank_inputs[j] for j in active if j != index])
            prcc = _correlation(_residual(rank_inputs[index], others), _residual(rank_output, others))
        records.append(RankSensitivity(name, spearman, prcc, srrc[name]))
    records.sort(key=lambda item: (item.prcc is None, -abs(item.prcc or 0.0)))
    return RankSensitivityResult(output_name, count, r_squared, tuple(records), tuple(warnings))


def _percentile(sorted_values: Sequence[float], probability: float) -> float:
    position = probability * (len(sorted_values) - 1)
    lower = int(position)
    upper = min(lower + 1, len(sorted_values) - 1)
    return sorted_values[lower] + (position - lower) * (sorted_values[upper] - sorted_values[lower])


def _indices(f_a: Sequence[float], f_b: Sequence[float], f_ab: Sequence[float],
             rows: Sequence[int]) -> tuple[float, float, float] | None:
    pooled = [f_a[r] for r in rows] + [f_b[r] for r in rows]
    mean = fsum(pooled) / len(pooled)
    variance = fsum((y - mean) ** 2 for y in pooled) / len(pooled)
    if variance <= 0.0:
        return None
    first = fsum(f_b[r] * (f_ab[r] - f_a[r]) for r in rows) / len(rows) / variance
    total = 0.5 * fsum((f_a[r] - f_ab[r]) ** 2 for r in rows) / len(rows) / variance
    return first, total, variance


def sobol_indices(model: Callable[[Mapping[str, float]], float], parameters: Sequence[ParameterRange], *,
                  base_samples: int = 1024, seed: int = 42, bootstrap_resamples: int = 500,
                  confidence: float = 0.95) -> SobolResult:
    """Estimate first-order and total Sobol indices of ``model`` [1]_.

    Parameters
    ----------
    model : callable
        Maps ``{name: value}`` to a finite scalar.
    parameters : sequence of ParameterRange
        Independent uniform inputs with distinct names.
    base_samples : int, default 1024
        Rows ``N`` of each base matrix (8 to 65 536).
    seed : int, default 42
        Seed of the pseudo-random design and the bootstrap.
    bootstrap_resamples : int, default 500
        Percentile-bootstrap resamples of the ``N`` rows (0 to 2 000; 0 gives
        point intervals).
    confidence : float, default 0.95
        Interval level in (0, 1).

    Returns
    -------
    SobolResult
        Indices, intervals, pooled mean and standard deviation and warnings.

    Raises
    ------
    InputError
        If the parameters are empty or names repeat, or a count is outside
        its range.
    DomainError
        If the model returns a non-finite value or a constant output.

    Notes
    -----
    Rows of ``A`` and ``B`` are independent uniform draws; ``A_B^(i)`` is
    ``A`` with column ``i`` from ``B``. Estimators as in SALib's
    ``analyze.sobol`` (Saltelli et al. 2010): ``S_i = mean(f_B (f_ABi -
    f_A)) / V`` and ``S_Ti = mean((f_A - f_ABi)^2) / (2 V)``, ``V`` the
    variance of the pooled ``f_A`` and ``f_B``, after every output is
    centred on the mean of all runs (SALib does the same; without it the
    first-order estimator's variance grows with ``(mean/sd)^2``). A non-influential input
    gives exactly zero for both. Monte Carlo error decreases as
    ``N^(-1/2)``; the bootstrap intervals report it.

    References
    ----------
    .. [1] SALTELLI-2010-CPC. Saltelli et al., Computer Physics
           Communications 181(2), 259-270, 2010.
    """

    if not parameters:
        raise InputError("At least one parameter is required.", field="parameters")
    names = [parameter.name for parameter in parameters]
    if len(set(names)) != len(names):
        raise InputError("Parameter names must be distinct.", field="parameters")
    if isinstance(base_samples, bool) or not isinstance(base_samples, int) \
            or not 8 <= base_samples <= MAXIMUM_SOBOL_BASE_SAMPLES:
        raise InputError(f"base_samples must be an integer in [8, {MAXIMUM_SOBOL_BASE_SAMPLES}].",
                         field="base_samples")
    if isinstance(bootstrap_resamples, bool) or not isinstance(bootstrap_resamples, int) \
            or not 0 <= bootstrap_resamples <= MAXIMUM_BOOTSTRAP_RESAMPLES:
        raise InputError(f"bootstrap_resamples must be an integer in [0, {MAXIMUM_BOOTSTRAP_RESAMPLES}].",
                         field="bootstrap_resamples")
    if not (isfinite(confidence) and 0.0 < confidence < 1.0):
        raise InputError("confidence must be in (0, 1).", field="confidence")
    if isinstance(seed, bool) or not isinstance(seed, int):
        raise InputError("seed must be an integer.", field="seed")
    rng = Random(seed)

    def draw() -> list[float]:
        return [p.minimum + (p.maximum - p.minimum) * rng.random() for p in parameters]

    matrix_a = [draw() for _ in range(base_samples)]
    matrix_b = [draw() for _ in range(base_samples)]

    def evaluate(row: Sequence[float]) -> float:
        value = float(model(dict(zip(names, row, strict=True))))
        if not isfinite(value):
            raise DomainError("The model returned a non-finite value.", code="non_finite_model_output",
                              field="model")
        return value

    f_a = [evaluate(row) for row in matrix_a]
    f_b = [evaluate(row) for row in matrix_b]
    f_ab = [[evaluate([*a[:i], b[i], *a[i + 1:]]) for a, b in zip(matrix_a, matrix_b, strict=True)]
            for i in range(len(parameters))]
    rows = range(base_samples)
    # The indices are invariant to an affine change of the output. Divide by the
    # largest magnitude, centre every output on the mean of all N (k + 2) runs (as
    # SALib does: an uncentred mean inflates the first-order estimator's variance by
    # ~(mean/sd)^2) and divide again by the largest centred magnitude, so that no
    # square or sum can overflow even for outputs near the double range.
    every = (*f_a, *f_b, *(y for column in f_ab for y in column))
    magnitude = max(abs(y) for y in every)
    if magnitude == 0.0:
        raise DomainError("The model output is constant; Sobol indices are undefined.", code="constant_output")
    centre = fsum(y / magnitude for y in every) / len(every)
    spread = max(abs(y / magnitude - centre) for y in every)
    if spread == 0.0:
        raise DomainError("The model output is constant; Sobol indices are undefined.", code="constant_output")

    def normalized(values: Sequence[float]) -> list[float]:
        return [(y / magnitude - centre) / spread for y in values]
    pooled_mean = fsum(y / magnitude for y in (*f_a, *f_b)) / (2 * base_samples)
    f_a, f_b = normalized(f_a), normalized(f_b)
    f_ab = [normalized(column) for column in f_ab]
    pooled = f_a + f_b
    unit_mean = fsum(pooled) / len(pooled)
    variance = fsum((y - unit_mean) ** 2 for y in pooled) / len(pooled)
    mean = pooled_mean * magnitude
    standard_deviation = sqrt(variance) * spread * magnitude
    point = [_indices(f_a, f_b, column, rows) for column in f_ab]
    if any(value is None for value in point):
        raise DomainError("The model output is constant; Sobol indices are undefined.", code="constant_output")

    tail = 0.5 * (1.0 - confidence)
    intervals: list[tuple[tuple[float, float], tuple[float, float]]] = []
    for i, column in enumerate(f_ab):
        estimate = point[i]
        assert estimate is not None
        if bootstrap_resamples == 0:
            intervals.append(((estimate[0], estimate[0]), (estimate[1], estimate[1])))
            continue
        # per-row terms once; each resample is then five sums over drawn rows
        first_terms = [fb * (fab - fa) for fa, fb, fab in zip(f_a, f_b, column, strict=True)]
        total_terms = [0.5 * (fa - fab) ** 2 for fa, fab in zip(f_a, column, strict=True)]
        firsts, totals = [], []
        for _ in range(bootstrap_resamples):
            drawn = [int(base_samples * rng.random()) for _ in rows]
            pooled_sum = sum(f_a[r] for r in drawn) + sum(f_b[r] for r in drawn)
            pooled_square = sum(f_a[r] * f_a[r] for r in drawn) + sum(f_b[r] * f_b[r] for r in drawn)
            pooled_mean = pooled_sum / (2 * base_samples)
            pooled_variance = pooled_square / (2 * base_samples) - pooled_mean * pooled_mean
            if pooled_variance <= 1.0e-12 * variance:
                continue
            firsts.append(sum(first_terms[r] for r in drawn) / base_samples / pooled_variance)
            totals.append(sum(total_terms[r] for r in drawn) / base_samples / pooled_variance)
        firsts.sort()
        totals.sort()
        intervals.append(((_percentile(firsts, tail), _percentile(firsts, 1.0 - tail)),
                          (_percentile(totals, tail), _percentile(totals, 1.0 - tail))))

    indices = []
    for name, estimate, (first_ci, total_ci) in zip(names, point, intervals, strict=True):
        assert estimate is not None
        indices.append(SobolIndex(name, estimate[0], estimate[1], first_ci, total_ci))
    warnings: list[WarningMessage] = []
    widest = max(max(i.first_order_interval[1] - i.first_order_interval[0],
                     i.total_interval[1] - i.total_interval[0]) for i in indices)
    if bootstrap_resamples and widest > 0.1:
        warnings.append(WarningMessage(
            code="sobol_sampling_error",
            message=(f"Bootstrap intervals are up to {widest:.2f} wide; increase base_samples "
                     "before ranking inputs.")))
    if any(i.first_order < 0.0 for i in indices):
        warnings.append(WarningMessage(
            code="negative_first_order_estimate", severity="info",
            message="A first-order estimate is negative: the true index is near zero within sampling error."))
    return SobolResult(base_samples, base_samples * (len(parameters) + 2), mean, standard_deviation, confidence,
                       tuple(indices), fsum(i.first_order for i in indices), tuple(warnings))
