"""Declared performance-factor ledgers; coupled engine validity remains open."""
from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass
from fractions import Fraction

from rocket_propulsion.core.errors import DomainError
from rocket_propulsion.core.metadata import WarningMessage
from rocket_propulsion.core.validation import coerce_fields, real_matrix, real_tuple, records, strings, text
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
    """Own a 4x4 correlation matrix and check exact positive semidefiniteness."""
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
    """One ordered attribution, replayed from retained inputs and its key."""
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
    """Owned component attribution; zero residuals are algebraic identities."""
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
    """Retain exact factor ratios until each separate public export."""
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
    force, impulse, flow, changes = _ledgers(source)
    physical = source.operating_basis == 'fixed-pressure-area'
    labels = ('Yanma / c*', 'Deşarj / debi' if physical else 'Deşarj adlı eski performans çarpanı',
              'Sınır tabakası', 'Iraksama')
    meanings = (('intrinsic thermochemical cstar ratio; excludes mass-flow discharge correction',
                 'actual/reference gas-state mass-flow discharge coefficient',
                 'net thrust-coefficient viscous factor', 'net thrust-coefficient divergence factor') if physical else
                ('lumped performance factor at fixed mass flow',
                 'legacy performance multiplier; not a physical throat mass-flow coefficient',
                 'lumped performance factor at fixed mass flow', 'lumped performance factor at fixed mass flow'))
    return [dict(key=key, label=labels[i], enabled=getattr(source, _FLAGS[i]),
        efficiency=getattr(source, _EFFS[i]),
        thrust_loss_n=_finite_rational(Fraction(source.ideal_thrust_n)*changes[i][0], 'Term thrust change'),
        specific_impulse_loss_s=_finite_rational(Fraction(source.ideal_specific_impulse_s)*changes[i][1], 'Term Isp change'),
        source=meanings[i], inputs=source, numerical_policy=_LOSS_POLICY) for i, key in enumerate(_KEYS)]


def _relative_variances(source):
    """First-order log-elasticity propagation, including signed covariance."""
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
    force, impulse, flow, changes = _ledgers(source)
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
    return dict(ideal_thrust_n=source.ideal_thrust_n, delivered_thrust_n=_finite_rational(f, 'Delivered thrust'),
        ideal_specific_impulse_s=source.ideal_specific_impulse_s,
        delivered_specific_impulse_s=_finite_rational(isp, 'Delivered Isp'),
        total_efficiency=_finite_rational(force, 'Thrust ratio'),
        specific_impulse_ratio=_finite_rational(impulse, 'Isp ratio'), mass_flow_ratio=_finite_rational(flow, 'Mass-flow ratio'),
        total_thrust_loss_n=_finite_rational(f0-f, 'Total thrust change'),
        total_specific_impulse_loss_s=_finite_rational(i0-isp, 'Total Isp change'),
        terms=tuple(LossTerm(**v) for v in _term_values(source)),
        thrust_balance_residual_n=0., specific_impulse_balance_residual_s=0.,
        one_sigma_relative=_sqrt_output(vf, 'Relative thrust uncertainty'),
        specific_impulse_one_sigma_relative=_sqrt_output(vi, 'Relative Isp uncertainty'),
        mass_flow_one_sigma_relative=_sqrt_output(vm, 'Relative flow uncertainty'),
        delivered_thrust_band_n=tuple(_affine_sqrt_output(f, f*f*vf, direction, 'Thrust proxy band') for direction in (-1, 1)),
        delivered_specific_impulse_band_s=tuple(_affine_sqrt_output(isp, isp*isp*vi, direction, 'Isp proxy band') for direction in (-1, 1)),
        inputs=source, operating_basis=source.operating_basis, numerical_policy=_LOSS_POLICY,
        uncertainty_kind=_UNCERTAINTY_KIND, coverage_qualified=False, warnings=tuple(warnings), assumptions=assumptions)


def calculate_loss_budget(
    *, ideal_thrust_n: float, ideal_specific_impulse_s: float,
    discharge_coefficient: float = .98, divergence_efficiency: float = .98,
    boundary_layer_efficiency: float = .985, combustion_efficiency: float = .97,
    enable_discharge_loss: bool = True, enable_divergence_loss: bool = True,
    enable_boundary_layer_loss: bool = True, enable_combustion_loss: bool = True,
    relative_uncertainty: float = .02, operating_basis: str = 'fixed-mass-flow-lumped',
    correlation_matrix: tuple[tuple[float, ...], ...] = _IDENTITY,
) -> PerformanceLossBudget:
    """Compute a declared factor budget and unqualified covariance proxy bands."""
    source = LossBudgetInputs(**locals())
    return PerformanceLossBudget(**_budget_values(source))
