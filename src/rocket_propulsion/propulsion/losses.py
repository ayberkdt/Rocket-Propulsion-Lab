"""Declared performance-factor ledgers with first-order covariance proxy bands.

Combustion (c*), discharge, boundary-layer and divergence factors are applied in
that order to an ideal thrust and Isp, and each change is attributed to its factor.
Relative standard uncertainties are propagated to first order with a declared
correlation matrix. Coupled engine validity, double counting between factors and
coverage of the bands are not established.

Notes
-----
Two operating bases are supported. ``'fixed-mass-flow-lumped'`` treats all four
inputs as effective multipliers of thrust and Isp at fixed mass flow; there
``discharge_coefficient`` is a performance multiplier, not the actual/ideal throat
mass-flow ratio. ``'fixed-pressure-area'`` holds the reference chamber pressure
times throat area: the intrinsic c* efficiency multiplies Isp and divides mass
flow, the discharge coefficient multiplies mass flow and divides Isp, and the
boundary-layer and divergence factors multiply the net thrust coefficient and
hence thrust and Isp. [1]_ [2]_ [3]_ Uncertainty propagation follows the first-order law of
propagation with covariances [4]_ [5]_.

References
----------
.. [1] SUTTON-BIBLARZ-2017. G. P. Sutton and O. Biblarz, "Rocket Propulsion
       Elements", 9th ed., Wiley, 2017.
.. [2] NASA-SP8120-DISCHARGE. NASA SP-8120, "Liquid Rocket Engine Nozzles", 1976,
       discharge coefficient and flow effects, printed pp. 10, 19-20 and 93.
.. [3] ROCKETISP-CSTAR-CD-CONVENTION. RocketIsp 0.1.13 project documentation,
       declared intrinsic c*/discharge-coefficient relation; not an independent NASA
       qualification.
.. [4] NIST-TN1297-COVARIANCE. NIST Technical Note 1297, Appendix A: first-order
       propagation of uncertainty with sensitivity coefficients and covariances.
.. [5] JCGM 100:2008, "Evaluation of Measurement Data - Guide to the Expression of
       Uncertainty in Measurement", Joint Committee for Guides in Metrology, 2008.
"""
from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass
from fractions import Fraction

from rocket_propulsion.core.errors import DomainError
from rocket_propulsion.core.metadata import WarningMessage
from rocket_propulsion.core.validation import (
    coerce_fields,
    real_matrix,
    real_tuple,
    records,
    strings,
    text,
)

from .burns.numeric_outputs import _finite_rational
from .burns.probability_numerics import _sqrt_output
from .burns.sample_statistics import _affine_sqrt_output

_LOSS_POLICY = 'exact_binary_factor_ledgers_and_algebraic_covariance_bands_v1'
_UNCERTAINTY_KIND = 'first_order_relative_factor_covariance_proxy'
_IDENTITY = tuple(tuple(float(i == j) for j in range(4)) for i in range(4))
_KEYS = ('combustion', 'discharge', 'boundary_layer', 'divergence')
_EFFS = ('combustion_efficiency', 'discharge_coefficient', 'boundary_layer_efficiency', 'divergence_efficiency')
_FLAGS = ('enable_combustion_loss', 'enable_discharge_loss', 'enable_boundary_layer_loss', 'enable_divergence_loss')


def _checked_correlation(value):
    """Return an owned 4x4 correlation matrix after an exact semidefiniteness check.

    Parameters
    ----------
    value : sequence of sequences
        Candidate correlation matrix ordered combustion, discharge, boundary layer,
        divergence.

    Returns
    -------
    tuple of tuple of float

    Raises
    ------
    DomainError
        If the matrix is not 4x4, not symmetric, has a non-unit diagonal or an entry
        outside [-1, 1], or is not positive semidefinite in exact rational
        elimination.
    """
    rows = real_matrix(value, 'Loss correlation matrix')
    if len(rows) != 4 or any(len(row) != 4 for row in rows):
        raise DomainError('Loss correlations require a 4x4 matrix.', field='correlation_matrix')
    if any(rows[i][i] != 1 for i in range(4)) or any(
            abs(rows[i][j]) > 1 or rows[i][j] != rows[j][i] for i in range(4) for j in range(4)):
        raise DomainError('Correlations must be symmetric, unit diagonal and in [-1,1].', field='correlation_matrix')
    work = [[Fraction(v) for v in row] for row in rows]
    for i in range(4):
        pivot = work[i][i]
        if pivot < 0 or pivot == 0 and any(work[j][i] for j in range(i+1, 4)):
            raise DomainError('Correlation matrix is not positive semidefinite.', field='correlation_matrix')
        if pivot:
            for j in range(i+1, 4):
                for k in range(i+1, 4):
                    work[j][k] -= work[j][i]*work[i][k]/pivot
    return rows


@dataclass(frozen=True, slots=True)
class LossBudgetInputs:
    """Validated inputs for a declared performance-factor budget.

    Attributes
    ----------
    ideal_thrust_n : float
        Reference ideal thrust, N; positive.
    ideal_specific_impulse_s : float
        Reference ideal specific impulse, s; positive.
    discharge_coefficient : float, default 0.98
        Discharge factor, in (0, 1]; meaning depends on ``operating_basis``.
    divergence_efficiency : float, default 0.98
        Nozzle divergence factor of the thrust coefficient, in (0, 1].
    boundary_layer_efficiency : float, default 0.985
        Viscous (boundary-layer) factor of the thrust coefficient, in (0, 1].
    combustion_efficiency : float, default 0.97
        Combustion (c*) efficiency, in (0, 1].
    enable_discharge_loss, enable_divergence_loss : bool, default True
        Switches; a disabled factor is treated as one.
    enable_boundary_layer_loss, enable_combustion_loss : bool, default True
        Switches; a disabled factor is treated as one.
    relative_uncertainty : float, default 0.02
        Common relative standard uncertainty of each enabled factor, in [0, 0.25].
    operating_basis : {'fixed-mass-flow-lumped', 'fixed-pressure-area'}
        Interpretation of the factors.
    correlation_matrix : tuple of tuple of float
        4x4 factor correlation matrix; identity by default.

    Raises
    ------
    DomainError
        If a value, switch, basis or the correlation matrix is invalid.

    Notes
    -----
    Two operating bases are supported. ``'fixed-mass-flow-lumped'`` treats all four
    inputs as effective multipliers of thrust and Isp at fixed mass flow; there
    ``discharge_coefficient`` is a performance multiplier, not the actual/ideal throat
    mass-flow ratio. ``'fixed-pressure-area'`` holds the reference chamber pressure
    times throat area: the intrinsic c* efficiency multiplies Isp and divides mass
    flow, the discharge coefficient multiplies mass flow and divides Isp, and the
    boundary-layer and divergence factors multiply the net thrust coefficient and
    hence thrust and Isp. [1]_ [2]_ [3]_

    References
    ----------
    .. [1] SUTTON-BIBLARZ-2017. G. P. Sutton and O. Biblarz, "Rocket Propulsion
           Elements", 9th ed., Wiley, 2017.
    .. [2] NASA-SP8120-DISCHARGE. NASA SP-8120, "Liquid Rocket Engine Nozzles", 1976,
           discharge coefficient and flow effects, printed pp. 10, 19-20 and 93.
    .. [3] ROCKETISP-CSTAR-CD-CONVENTION. RocketIsp 0.1.13 project documentation,
           declared intrinsic c*/discharge-coefficient relation; not an independent NASA
           qualification.
    """
    ideal_thrust_n: float
    ideal_specific_impulse_s: float
    discharge_coefficient: float = .98
    divergence_efficiency: float = .98
    boundary_layer_efficiency: float = .985
    combustion_efficiency: float = .97
    enable_discharge_loss: bool = True
    enable_divergence_loss: bool = True
    enable_boundary_layer_loss: bool = True
    enable_combustion_loss: bool = True
    relative_uncertainty: float = .02
    operating_basis: str = 'fixed-mass-flow-lumped'
    correlation_matrix: tuple[tuple[float, ...], ...] = _IDENTITY

    def __post_init__(self):
        coerce_fields(self, 'Loss input')
        if self.ideal_thrust_n <= 0 or self.ideal_specific_impulse_s <= 0:
            raise DomainError('Ideal thrust and Isp must be positive.', field='Ideal performance')
        if self.operating_basis not in ('fixed-mass-flow-lumped', 'fixed-pressure-area'):
            raise DomainError('Unknown loss operating basis.', field='operating_basis')
        if not 0 <= self.relative_uncertainty <= .25:
            raise DomainError('Factor relative standard uncertainty must be in [0,.25].', field='relative_uncertainty')
        for name in _EFFS:
            if not 0 < getattr(self, name) <= 1:
                raise DomainError('Factor must be in (0,1].', field=name)
        for name in _FLAGS:
            if not isinstance(getattr(self, name), bool):
                raise DomainError('Loss switch must be boolean.', field=name)
        object.__setattr__(self, 'correlation_matrix', _checked_correlation(self.correlation_matrix))


def _owned(value, kind, label):
    """Return a typed record of ``kind`` from a record or mapping.

    Warning records also have their text fields and severity normalized.

    Parameters
    ----------
    value : object
        Record or mapping of its fields.
    kind : type
        Expected record type.
    label : str
        Label used in errors.

    Returns
    -------
    object
        Instance of ``kind``.

    Raises
    ------
    DomainError
        If the value is malformed, of the wrong type, or has an unknown severity.
    """
    if isinstance(value, Mapping):
        try:
            value = kind(**dict(value))
        except (TypeError, ValueError, KeyError, IndexError, OverflowError) as error:
            raise DomainError('Malformed retained ' + label + '.', field=label) from error
    if not isinstance(value, kind):
        raise DomainError('Retained ' + label + ' has the wrong type.', field=label)
    if kind is WarningMessage:
        code = text(value.code, 'Loss warning code')
        message = text(value.message, 'Loss warning message')
        severity = text(value.severity, 'Loss warning severity')
        field = None if value.field is None else text(value.field, 'Loss warning field')
        if severity not in ('info', 'warning', 'critical'):
            raise DomainError('Unknown loss warning severity.', field='Loss warning severity')
        value = WarningMessage(code, message, severity, field)
    return value


@dataclass(frozen=True, slots=True)
class LossTerm:
    """One ordered factor attribution, replayable from the retained inputs and its key.

    Attributes
    ----------
    key : {'combustion', 'discharge', 'boundary_layer', 'divergence'}
        Factor identifier.
    label : str
        Display label.
    enabled : bool
        Whether the factor is applied.
    efficiency : float
        Declared factor value.
    thrust_loss_n : float
        Signed thrust change attributed to this factor in sequence, N.
    specific_impulse_loss_s : float
        Signed Isp change attributed to this factor in sequence, s.
    source : str
        Physical meaning of the factor under the operating basis.
    inputs : LossBudgetInputs
        Retained inputs.
    numerical_policy : str
        Identifier of the arithmetic and rounding policy.

    Raises
    ------
    DomainError
        If the key, switch or policy is unknown or a field does not replay.

    References
    ----------
    .. [1] SUTTON-BIBLARZ-2017. G. P. Sutton and O. Biblarz, "Rocket Propulsion
           Elements", 9th ed., Wiley, 2017.
    """
    key: str
    label: str
    enabled: bool
    efficiency: float
    thrust_loss_n: float
    specific_impulse_loss_s: float
    source: str
    inputs: LossBudgetInputs | None = None
    numerical_policy: str = _LOSS_POLICY

    def __post_init__(self):
        coerce_fields(self, 'Loss term')
        source = _owned(self.inputs, LossBudgetInputs, 'Loss inputs')
        object.__setattr__(self, 'inputs', source)
        if self.key not in _KEYS or not isinstance(self.enabled, bool) or self.numerical_policy != _LOSS_POLICY:
            raise DomainError('Unknown loss term, switch or numerical policy.', field='Loss term')
        expected = _term_values(source)[_KEYS.index(self.key)]
        if any(getattr(self, name) != value for name, value in expected.items()):
            raise DomainError('Loss term does not replay retained source inputs.', field='Loss term')


@dataclass(frozen=True, slots=True)
class PerformanceLossBudget:
    """Ordered factor attribution with first-order uncertainty proxy bands.

    On construction the record evaluates the budget again from ``inputs`` and refuses
    any field that the evaluation does not reproduce. Zero balance residuals are
    algebraic identities of the sequential attribution, not measured validation, and
    ``coverage_qualified`` must be False.

    Attributes
    ----------
    ideal_thrust_n, delivered_thrust_n : float
        Reference and factored thrust, N.
    ideal_specific_impulse_s, delivered_specific_impulse_s : float
        Reference and factored Isp, s.
    total_efficiency : float
        Delivered/ideal thrust ratio.
    total_thrust_loss_n, total_specific_impulse_loss_s : float
        Ideal minus delivered thrust (N) and Isp (s).
    terms : tuple of LossTerm
        Ordered attributions.
    thrust_balance_residual_n, specific_impulse_balance_residual_s : float
        Zero by construction.
    one_sigma_relative : float
        First-order relative standard uncertainty of thrust.
    delivered_thrust_band_n, delivered_specific_impulse_band_s : tuple of float
        Delivered value times (1 -/+ relative standard uncertainty).
    assumptions : tuple of str
        Model assumptions.
    inputs : LossBudgetInputs
        Retained inputs.
    operating_basis : str
        Interpretation of the factors.
    specific_impulse_ratio, mass_flow_ratio : float
        Delivered/ideal Isp and mass-flow ratios.
    specific_impulse_one_sigma_relative, mass_flow_one_sigma_relative : float
        First-order relative standard uncertainties of Isp and mass flow.
    numerical_policy, uncertainty_kind : str
        Identifiers of the arithmetic and uncertainty policies.
    coverage_qualified : bool
        Always False.
    warnings : tuple of WarningMessage
        Scope warnings.

    Raises
    ------
    DomainError
        If a policy is unknown, a band is malformed, or a field does not replay.

    References
    ----------
    .. [1] SUTTON-BIBLARZ-2017. G. P. Sutton and O. Biblarz, "Rocket Propulsion
           Elements", 9th ed., Wiley, 2017.
    .. [2] NIST-TN1297-COVARIANCE. NIST Technical Note 1297, Appendix A: first-order
           propagation of uncertainty with sensitivity coefficients and covariances.
    .. [3] JCGM 100:2008, "Evaluation of Measurement Data - Guide to the Expression of
           Uncertainty in Measurement", Joint Committee for Guides in Metrology, 2008.
    """
    ideal_thrust_n: float
    delivered_thrust_n: float
    ideal_specific_impulse_s: float
    delivered_specific_impulse_s: float
    total_efficiency: float
    total_thrust_loss_n: float
    total_specific_impulse_loss_s: float
    terms: tuple[LossTerm, ...]
    thrust_balance_residual_n: float
    specific_impulse_balance_residual_s: float
    one_sigma_relative: float
    delivered_thrust_band_n: tuple[float, float]
    delivered_specific_impulse_band_s: tuple[float, float]
    assumptions: tuple[str, ...]
    inputs: LossBudgetInputs | None = None
    operating_basis: str = 'fixed-mass-flow-lumped'
    specific_impulse_ratio: float = 1.
    mass_flow_ratio: float = 1.
    specific_impulse_one_sigma_relative: float = 0.
    mass_flow_one_sigma_relative: float = 0.
    numerical_policy: str = _LOSS_POLICY
    uncertainty_kind: str = _UNCERTAINTY_KIND
    coverage_qualified: bool = False
    warnings: tuple[WarningMessage, ...] = ()

    def __post_init__(self):
        coerce_fields(self, 'Loss result')
        source = _owned(self.inputs, LossBudgetInputs, 'Loss inputs')
        object.__setattr__(self, 'inputs', source)
        if self.numerical_policy != _LOSS_POLICY or self.uncertainty_kind != _UNCERTAINTY_KIND:
            raise DomainError('Unknown loss arithmetic/uncertainty policy.', field='Loss policy')
        if self.coverage_qualified is not False:
            raise DomainError('Qualified uncertainty coverage is not established.', field='coverage_qualified')
        terms = tuple(_owned(v, LossTerm, 'Loss term') for v in records(self.terms, (LossTerm, Mapping), 'Loss terms'))
        warnings = tuple(_owned(v, WarningMessage, 'Loss warning') for v in records(
            self.warnings, (WarningMessage, Mapping), 'Loss warnings'))
        object.__setattr__(self, 'terms', terms)
        object.__setattr__(self, 'warnings', warnings)
        object.__setattr__(self, 'assumptions', strings(self.assumptions, 'Loss assumptions'))
        for name in ('delivered_thrust_band_n', 'delivered_specific_impulse_band_s'):
            value = real_tuple(getattr(self, name), name)
            if len(value) != 2:
                raise DomainError('Band must contain two endpoints.', field=name)
            object.__setattr__(self, name, value)
        expected = _budget_values(source)
        if any(getattr(self, name) != value for name, value in expected.items()):
            raise DomainError('Loss outputs do not replay retained inputs and basis.', field='Loss outputs')


def _ledgers(source):
    """Apply the factors in order and keep exact ratios until export.

    Parameters
    ----------
    source : LossBudgetInputs
        Validated inputs.

    Returns
    -------
    tuple
        ``(thrust_ratio, isp_ratio, flow_ratio, changes)`` as exact fractions, where
        ``changes`` holds the per-factor (thrust, Isp) ratio decrements.

    Notes
    -----
    Two operating bases are supported. ``'fixed-mass-flow-lumped'`` treats all four
    inputs as effective multipliers of thrust and Isp at fixed mass flow; there
    ``discharge_coefficient`` is a performance multiplier, not the actual/ideal throat
    mass-flow ratio. ``'fixed-pressure-area'`` holds the reference chamber pressure
    times throat area: the intrinsic c* efficiency multiplies Isp and divides mass
    flow, the discharge coefficient multiplies mass flow and divides Isp, and the
    boundary-layer and divergence factors multiply the net thrust coefficient and
    hence thrust and Isp. [1]_ [2]_ [3]_

    References
    ----------
    .. [1] SUTTON-BIBLARZ-2017. G. P. Sutton and O. Biblarz, "Rocket Propulsion
           Elements", 9th ed., Wiley, 2017.
    .. [2] NASA-SP8120-DISCHARGE. NASA SP-8120, "Liquid Rocket Engine Nozzles", 1976,
           discharge coefficient and flow effects, printed pp. 10, 19-20 and 93.
    .. [3] ROCKETISP-CSTAR-CD-CONVENTION. RocketIsp 0.1.13 project documentation,
           declared intrinsic c*/discharge-coefficient relation; not an independent NASA
           qualification.
    """
    force, impulse, flow = Fraction(1), Fraction(1), Fraction(1)
    physical = source.operating_basis == 'fixed-pressure-area'
    terms = []
    for index, (eff, flag) in enumerate(zip(_EFFS, _FLAGS)):
        eta = Fraction(getattr(source, eff)) if getattr(source, flag) else Fraction(1)
        before_f, before_i = force, impulse
        if physical:
            if index == 0:
                impulse *= eta
                flow /= eta
            elif index == 1:
                impulse /= eta
                flow *= eta
            else:
                force *= eta
                impulse *= eta
        else:
            force *= eta
            impulse *= eta
        terms.append((before_f-force, before_i-impulse))
    return force, impulse, flow, terms


def _term_values(source):
    """Return the field values of the four ordered :class:`LossTerm` records.

    Parameters
    ----------
    source : LossBudgetInputs
        Validated inputs.

    Returns
    -------
    list of dict

    References
    ----------
    .. [1] SUTTON-BIBLARZ-2017. G. P. Sutton and O. Biblarz, "Rocket Propulsion
           Elements", 9th ed., Wiley, 2017.
    .. [2] NASA-SP8120-DISCHARGE. NASA SP-8120, "Liquid Rocket Engine Nozzles", 1976,
           discharge coefficient and flow effects, printed pp. 10, 19-20 and 93.
    """
    _force, _impulse, _flow, changes = _ledgers(source)
    physical = source.operating_basis == 'fixed-pressure-area'
    labels = ('Yanma / c*', 'Deşarj / debi' if physical else 'Deşarj adlı eski performans çarpanı',
              'Sınır tabakası', 'Iraksama')
    meanings = (('intrinsic thermochemical cstar ratio; excludes mass-flow discharge correction',
                 'actual/reference gas-state mass-flow discharge coefficient',
                 'net thrust-coefficient viscous factor', 'net thrust-coefficient divergence factor') if physical else
                ('lumped performance factor at fixed mass flow',
                 'legacy performance multiplier; not a physical throat mass-flow coefficient',
                 'lumped performance factor at fixed mass flow', 'lumped performance factor at fixed mass flow'))
    return [{'key': key, 'label': labels[i], 'enabled': getattr(source, _FLAGS[i]),
        'efficiency': getattr(source, _EFFS[i]),
        'thrust_loss_n': _finite_rational(Fraction(source.ideal_thrust_n)*changes[i][0], 'Term thrust change'),
        'specific_impulse_loss_s': _finite_rational(Fraction(source.ideal_specific_impulse_s)*changes[i][1], 'Term Isp change'),
        'source': meanings[i], 'inputs': source, 'numerical_policy': _LOSS_POLICY} for i, key in enumerate(_KEYS)]


def _relative_variances(source):
    """Propagate the factor uncertainties to relative thrust, Isp and flow variances.

    Parameters
    ----------
    source : LossBudgetInputs
        Validated inputs.

    Returns
    -------
    tuple of Fraction
        Relative variances of thrust, Isp and mass flow.

    Notes
    -----
    For a product of power-law factors the first-order relative variance is
    u**2 * e^T R e, with u the common relative standard uncertainty, R the
    correlation matrix and e the log-elasticities of the output with respect to the
    enabled factors (zero for disabled ones), including signed covariance terms
    [1]_ [2]_.

    References
    ----------
    .. [1] NIST-TN1297-COVARIANCE. NIST Technical Note 1297, Appendix A: first-order
           propagation of uncertainty with sensitivity coefficients and covariances.
    .. [2] JCGM 100:2008, "Evaluation of Measurement Data - Guide to the Expression of
           Uncertainty in Measurement", Joint Committee for Guides in Metrology, 2008.
    """
    elasticities = ((0, 0, 1, 1), (1, -1, 1, 1), (-1, 1, 0, 0)) if (
        source.operating_basis == 'fixed-pressure-area') else ((1,)*4, (1,)*4, (0,)*4)
    variance = []
    for row in elasticities:
        e = [v if getattr(source, flag) else 0 for v, flag in zip(row, _FLAGS)]
        quadratic = sum((Fraction(source.correlation_matrix[i][j])*e[i]*e[j]
                         for i in range(4) for j in range(4)), Fraction(0))
        variance.append(Fraction(source.relative_uncertainty)**2*quadratic)
    return tuple(variance)


def _budget_values(source):
    """Evaluate the factor budget and proxy bands and return all result fields.

    Parameters
    ----------
    source : LossBudgetInputs
        Validated inputs.

    Returns
    -------
    dict
        Field values for :class:`PerformanceLossBudget`.

    Raises
    ------
    DomainError
        If an output is not representable.

    References
    ----------
    .. [1] SUTTON-BIBLARZ-2017. G. P. Sutton and O. Biblarz, "Rocket Propulsion
           Elements", 9th ed., Wiley, 2017.
    .. [2] NIST-TN1297-COVARIANCE. NIST Technical Note 1297, Appendix A: first-order
           propagation of uncertainty with sensitivity coefficients and covariances.
    """
    force, impulse, flow, _changes = _ledgers(source)
    f0, i0 = Fraction(source.ideal_thrust_n), Fraction(source.ideal_specific_impulse_s)
    f, isp = f0*force, i0*impulse
    vf, vi, vm = _relative_variances(source)
    physical = source.operating_basis == 'fixed-pressure-area'
    warnings = [WarningMessage('factor_model_not_empirically_validated',
        'Coefficients may share or overlap physical effects; sequential attribution does not prove absence of double counting.'),
        WarningMessage('first_order_uncertainty_unqualified',
        'Bands use first-order factor covariance only. Input uncertainties/correlations, nonlinear effects and coverage are not validated; ideal-reference uncertainty is excluded.')]
    if physical:
        warnings.append(WarningMessage('coefficient_compatibility_not_validated',
            'Fixed Pc*At coefficient ratios do not solve gas/nozzle compatibility. Isp may exceed the reference when coefficients are independently varied; this is not an achieved performance prediction.', 'critical'))
    else:
        warnings.append(WarningMessage('legacy_discharge_is_performance_multiplier',
            'In fixed-mass-flow-lumped mode, discharge_coefficient is the legacy F/Isp multiplier, not actual/ideal throat mass flow.', 'critical'))
    assumptions = ('ordered algebraic attribution; residual zero is not measured validation',
        'four factors ordered combustion, discharge, boundary_layer, divergence',
        'relative_uncertainty is the same relative input standard uncertainty for each enabled factor',
        'first-order covariance propagation; no probability coverage qualification; ideal inputs treated fixed',
        ('fixed reference chamber-pressure times throat-area; intrinsic cstar and independent net Cf factors; no absolute state or coupled flow solution'
         if physical else 'fixed mass flow; all four inputs are effective performance multipliers'))
    return {'ideal_thrust_n': source.ideal_thrust_n, 'delivered_thrust_n': _finite_rational(f, 'Delivered thrust'),
        'ideal_specific_impulse_s': source.ideal_specific_impulse_s,
        'delivered_specific_impulse_s': _finite_rational(isp, 'Delivered Isp'),
        'total_efficiency': _finite_rational(force, 'Thrust ratio'),
        'specific_impulse_ratio': _finite_rational(impulse, 'Isp ratio'), 'mass_flow_ratio': _finite_rational(flow, 'Mass-flow ratio'),
        'total_thrust_loss_n': _finite_rational(f0-f, 'Total thrust change'),
        'total_specific_impulse_loss_s': _finite_rational(i0-isp, 'Total Isp change'),
        'terms': tuple(LossTerm(**v) for v in _term_values(source)),
        'thrust_balance_residual_n': 0., 'specific_impulse_balance_residual_s': 0.,
        'one_sigma_relative': _sqrt_output(vf, 'Relative thrust uncertainty'),
        'specific_impulse_one_sigma_relative': _sqrt_output(vi, 'Relative Isp uncertainty'),
        'mass_flow_one_sigma_relative': _sqrt_output(vm, 'Relative flow uncertainty'),
        'delivered_thrust_band_n': tuple(_affine_sqrt_output(f, f*f*vf, direction, 'Thrust proxy band') for direction in (-1, 1)),
        'delivered_specific_impulse_band_s': tuple(_affine_sqrt_output(isp, isp*isp*vi, direction, 'Isp proxy band') for direction in (-1, 1)),
        'inputs': source, 'operating_basis': source.operating_basis, 'numerical_policy': _LOSS_POLICY,
        'uncertainty_kind': _UNCERTAINTY_KIND, 'coverage_qualified': False, 'warnings': tuple(warnings), 'assumptions': assumptions}


def calculate_loss_budget(
    *, ideal_thrust_n: float, ideal_specific_impulse_s: float,
    discharge_coefficient: float = .98, divergence_efficiency: float = .98,
    boundary_layer_efficiency: float = .985, combustion_efficiency: float = .97,
    enable_discharge_loss: bool = True, enable_divergence_loss: bool = True,
    enable_boundary_layer_loss: bool = True, enable_combustion_loss: bool = True,
    relative_uncertainty: float = .02, operating_basis: str = 'fixed-mass-flow-lumped',
    correlation_matrix: tuple[tuple[float, ...], ...] = _IDENTITY,
) -> PerformanceLossBudget:
    """Compute a declared performance-factor budget and unqualified covariance bands.

    Parameters
    ----------
    ideal_thrust_n : float
        Reference ideal thrust, N; positive.
    ideal_specific_impulse_s : float
        Reference ideal specific impulse, s; positive.
    discharge_coefficient, divergence_efficiency : float
        Discharge and divergence factors, in (0, 1].
    boundary_layer_efficiency, combustion_efficiency : float
        Boundary-layer and combustion (c*) factors, in (0, 1].
    enable_discharge_loss, enable_divergence_loss : bool
        Factor switches.
    enable_boundary_layer_loss, enable_combustion_loss : bool
        Factor switches.
    relative_uncertainty : float
        Common relative standard uncertainty of each enabled factor, in [0, 0.25].
    operating_basis : {'fixed-mass-flow-lumped', 'fixed-pressure-area'}
        Interpretation of the factors.
    correlation_matrix : tuple of tuple of float
        4x4 positive semidefinite correlation matrix.

    Returns
    -------
    PerformanceLossBudget
        Delivered thrust and Isp, ordered attributions, ratios, proxy bands and
        warnings.

    Raises
    ------
    DomainError
        If an input is invalid or an output is not representable.

    Notes
    -----
    Two operating bases are supported. ``'fixed-mass-flow-lumped'`` treats all four
    inputs as effective multipliers of thrust and Isp at fixed mass flow; there
    ``discharge_coefficient`` is a performance multiplier, not the actual/ideal throat
    mass-flow ratio. ``'fixed-pressure-area'`` holds the reference chamber pressure
    times throat area: the intrinsic c* efficiency multiplies Isp and divides mass
    flow, the discharge coefficient multiplies mass flow and divides Isp, and the
    boundary-layer and divergence factors multiply the net thrust coefficient and
    hence thrust and Isp. [1]_ [2]_ [3]_ The factors may share or overlap physical effects;
    sequential attribution does not prove the absence of double counting. Bands use
    first-order factor covariance only [4]_; nonlinear effects, the uncertainty of the
    ideal reference and probability coverage are not included.

    References
    ----------
    .. [1] SUTTON-BIBLARZ-2017. G. P. Sutton and O. Biblarz, "Rocket Propulsion
           Elements", 9th ed., Wiley, 2017.
    .. [2] NASA-SP8120-DISCHARGE. NASA SP-8120, "Liquid Rocket Engine Nozzles", 1976,
           discharge coefficient and flow effects, printed pp. 10, 19-20 and 93.
    .. [3] ROCKETISP-CSTAR-CD-CONVENTION. RocketIsp 0.1.13 project documentation,
           declared intrinsic c*/discharge-coefficient relation; not an independent NASA
           qualification.
    .. [4] NIST-TN1297-COVARIANCE. NIST Technical Note 1297, Appendix A: first-order
           propagation of uncertainty with sensitivity coefficients and covariances.
    """
    source = LossBudgetInputs(**locals())
    return PerformanceLossBudget(**_budget_values(source))
