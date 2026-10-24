"""Conditional two-phase exit-state accounting, not a coupled particle solver."""
from __future__ import annotations

from collections.abc import Mapping
from collections.abc import Set as AbstractSet
from dataclasses import dataclass
from fractions import Fraction
from itertools import islice

from rocket_propulsion.core.errors import DomainError
from rocket_propulsion.core.metadata import WarningMessage
from rocket_propulsion.core.validation import (
    _ordered_items,
    coerce_fields,
    real,
    records,
    strings,
    text,
)

from .burns.numeric_outputs import _finite_rational
from .burns.probability_numerics import _sqrt_output
from .burns.sample_statistics import _affine_sqrt_output

_VELOCITY = 'exit-velocity-ratio'
_ENERGY = 'specific-kinetic-energy-ratio'
_TWO_PHASE_POLICY = 'exact_binary_mixture_ledgers_and_algebraic_sqrt_exports_v1'
_BAND_KIND = 'conditional_exit_state_scenario_range'
_MODEL = 'fixed total mass flow and prescribed gas exit state; separated pressure thrust'
_PAIRS = ('particle_coupling_efficiency_range', 'particle_velocity_ratio_range',
          'momentum_efficiency_range', 'performance_efficiency_range', 'total_thrust_ratio_range',
          'thrust_range_n', 'specific_impulse_range_s')


@dataclass(frozen=True, slots=True)
class TwoPhaseInputs:
    ideal_thrust_n: float
    ideal_specific_impulse_s: float
    condensed_mass_fraction: float
    minimum_particle_coupling_efficiency: float
    maximum_particle_coupling_efficiency: float
    reference_pressure_thrust_n: float = 0.
    coupling_definition: str = _VELOCITY

    def __post_init__(self):
        coerce_fields(self, 'Two-phase input')
        if self.ideal_thrust_n <= 0 or self.ideal_specific_impulse_s <= 0:
            raise DomainError('Reference thrust and Isp must be positive.', field='Reference performance')
        if Fraction(self.ideal_thrust_n)-Fraction(self.reference_pressure_thrust_n) <= 0:
            raise DomainError('Reference momentum thrust must be positive.', field='reference_pressure_thrust_n')
        if not 0 <= self.condensed_mass_fraction <= 1:
            raise DomainError('Condensed exhaust mass fraction must be in [0,1].', field='condensed_mass_fraction')
        lower, upper = self.minimum_particle_coupling_efficiency, self.maximum_particle_coupling_efficiency
        if not 0 <= lower <= upper <= 1:
            raise DomainError('Ordered particle/gas ratios must be in [0,1].', field='Particle coupling range')
        if self.coupling_definition not in (_VELOCITY, _ENERGY):
            raise DomainError('Unknown coupling-ratio definition.', field='coupling_definition')


def _retained_inputs(value):
    if isinstance(value, Mapping):
        try:
            value = TwoPhaseInputs(**dict(value))
        except (TypeError, ValueError, KeyError, IndexError, OverflowError) as error:
            raise DomainError('Malformed retained two-phase inputs.', field='Two-phase inputs') from error
    if not isinstance(value, TwoPhaseInputs):
        raise DomainError('Typed retained two-phase inputs are required.', field='Two-phase inputs')
    return value


def _two_phase_warnings(value):
    """Bound the known warning collection and normalize literal text before replay."""
    if isinstance(value, (str, bytes, bytearray, memoryview, Mapping, AbstractSet)):
        raise DomainError('Warnings must be ordered records.', field='Two-phase warnings')
    try:
        items = records(islice(value, 5), (WarningMessage, Mapping), 'Two-phase warnings')
    except (TypeError, ValueError, KeyError, IndexError, OverflowError) as error:
        raise DomainError('Malformed two-phase warnings.', field='Two-phase warnings') from error
    if len(items) > 4:
        raise DomainError('Too many retained warnings.', field='Two-phase warnings')
    result = []
    for warning in items:
        if isinstance(warning, Mapping):
            try:
                warning = WarningMessage(**dict(warning))
            except (TypeError, ValueError, KeyError, IndexError, OverflowError) as error:
                raise DomainError('Malformed retained warning.', field='Two-phase warnings') from error
        code, message = text(warning.code, 'Warning code'), text(warning.message, 'Warning message')
        severity = text(warning.severity, 'Warning severity')
        field = None if warning.field is None else text(warning.field, 'Warning field')
        if severity not in ('info', 'warning', 'critical'):
            raise DomainError('Unknown warning severity.', field='Warning severity')
        result.append(WarningMessage(code, message, severity, field))
    return tuple(result)


@dataclass(frozen=True, slots=True)
class TwoPhasePerformanceBand:
    """Replayed conditional scenario endpoints, without qualification claims."""
    condensed_mass_fraction: float
    particle_coupling_efficiency_range: tuple[float, float]
    performance_efficiency_range: tuple[float, float]
    thrust_range_n: tuple[float, float]
    specific_impulse_range_s: tuple[float, float]
    warnings: tuple[WarningMessage, ...]
    model: str
    inputs: TwoPhaseInputs | None = None
    reference_pressure_thrust_n: float = 0.
    coupling_definition: str = _VELOCITY
    particle_velocity_ratio_range: tuple[float, float] = (0., 0.)
    momentum_efficiency_range: tuple[float, float] = (0., 0.)
    total_thrust_ratio_range: tuple[float, float] = (0., 0.)
    contains_nonpositive_thrust: bool = False
    physical_bound_qualified: bool = False
    statistical_coverage_qualified: bool = False
    band_kind: str = _BAND_KIND
    numerical_policy: str = _TWO_PHASE_POLICY
    assumptions: tuple[str, ...] = ()

    def __post_init__(self):
        coerce_fields(self, 'Two-phase result')
        source = _retained_inputs(self.inputs)
        object.__setattr__(self, 'inputs', source)
        if self.numerical_policy != _TWO_PHASE_POLICY or self.band_kind != _BAND_KIND:
            raise DomainError('Unknown two-phase policy or range kind.', field='Two-phase policy')
        if self.physical_bound_qualified is not False or self.statistical_coverage_qualified is not False:
            raise DomainError('Physical/statistical qualification has not been established.', field='Qualification')
        if not isinstance(self.contains_nonpositive_thrust, bool):
            raise DomainError('Conditional thrust flag must be boolean.', field='contains_nonpositive_thrust')
        for name in _PAIRS:
            pair = tuple(real(v, name) for v in _ordered_items(getattr(self, name), name, size=2))
            object.__setattr__(self, name, pair)
        object.__setattr__(self, 'warnings', _two_phase_warnings(self.warnings))
        object.__setattr__(self, 'assumptions', strings(self.assumptions, 'Two-phase assumptions'))
        expected = _band_values(source)
        if any(getattr(self, name) != value for name, value in expected.items()):
            raise DomainError('Two-phase outputs do not reproduce retained inputs.', field='Two-phase outputs')


def _nonpositive_algebraic(mean, radicand):
    """Compare mean+sqrt(radicand) with zero using exact squared quantities."""
    if mean > 0:
        return False
    return radicand <= mean*mean


def _endpoint(source, coupling):
    force, isp = Fraction(source.ideal_thrust_n), Fraction(source.ideal_specific_impulse_s)
    pressure, chi = Fraction(source.reference_pressure_thrust_n), Fraction(source.condensed_mass_fraction)
    momentum, k = force-pressure, Fraction(coupling)
    if source.coupling_definition == _VELOCITY:
        speed = k
        mixture = 1-chi+chi*speed
        thrust = momentum*mixture+pressure
        return {'velocity': _finite_rational(speed, 'Particle velocity ratio'),
            'mixture': _finite_rational(mixture, 'Momentum multiplier'),
            'thrust': _finite_rational(thrust, 'Conditional thrust'),
            'ratio': _finite_rational(thrust/force, 'Conditional thrust ratio'),
            'isp': _finite_rational(isp*thrust/force, 'Conditional Isp'), 'nonpositive': thrust <= 0}
    # Convert specific kinetic energy ratio to velocity without exporting the
    # root into dependent outputs. Each final affine square root rounds alone.
    scaled_particle = momentum*chi
    mean = force-scaled_particle
    radicand = scaled_particle*scaled_particle*k
    return {'velocity': _sqrt_output(k, 'Particle velocity ratio'),
        'mixture': _affine_sqrt_output(1-chi, chi*chi*k, 1, 'Momentum multiplier'),
        'thrust': _affine_sqrt_output(mean, radicand, 1, 'Conditional thrust'),
        'ratio': _affine_sqrt_output(mean/force, radicand/(force*force), 1, 'Conditional thrust ratio'),
        'isp': _affine_sqrt_output(isp*mean/force, isp*isp*radicand/(force*force), 1, 'Conditional Isp'),
        'nonpositive': _nonpositive_algebraic(mean, radicand)}


def _band_values(source):
    lower, upper = source.minimum_particle_coupling_efficiency, source.maximum_particle_coupling_efficiency
    endpoints = (_endpoint(source, lower), _endpoint(source, upper))
    nonpositive = any(row['nonpositive'] for row in endpoints)
    warnings = [WarningMessage('two_phase_band_only',
        'User-supplied exit-state scenarios, not a validated bound. Gas momentum/energy feedback, particle distributions, thermal lag and phase change are not solved.', 'critical')]
    if source.reference_pressure_thrust_n == 0:
        warnings.append(WarningMessage('zero_pressure_term_assumed',
            'Zero reference pressure thrust is an input assumption, not evidence of matched nozzle expansion.'))
    if source.condensed_mass_fraction == 1:
        warnings.append(WarningMessage('no_gas_phase_limit',
            'All-condensed exhaust is a mathematical limit; the prescribed gas-state reference does not establish a realizable nozzle flow.', 'critical'))
    if nonpositive:
        warnings.append(WarningMessage('nonpositive_conditional_thrust',
            'At least one prescribed-state endpoint has nonpositive thrust. Actual backpressure/nozzle compatibility is not solved.', 'critical'))
    ratios = tuple(row['ratio'] for row in endpoints)
    return {'condensed_mass_fraction': source.condensed_mass_fraction,
        'particle_coupling_efficiency_range': (lower, upper), 'performance_efficiency_range': ratios,
        'thrust_range_n': tuple(row['thrust'] for row in endpoints), 'specific_impulse_range_s': tuple(row['isp'] for row in endpoints),
        'warnings': tuple(warnings), 'model': _MODEL, 'inputs': source,
        'reference_pressure_thrust_n': source.reference_pressure_thrust_n, 'coupling_definition': source.coupling_definition,
        'particle_velocity_ratio_range': tuple(row['velocity'] for row in endpoints),
        'momentum_efficiency_range': tuple(row['mixture'] for row in endpoints), 'total_thrust_ratio_range': ratios,
        'contains_nonpositive_thrust': nonpositive, 'physical_bound_qualified': False, 'statistical_coverage_qualified': False,
        'band_kind': _BAND_KIND, 'numerical_policy': _TWO_PHASE_POLICY,
        'assumptions': ('fixed total mass flow and prescribed gas exit state; no coupled gas/particle solution',
            'condensed_mass_fraction is exhaust mass fraction, not initial propellant loading',
            'pressure thrust retained separately; positive reference momentum required',
            'specific kinetic energy ratio is particle/gas per-unit-mass energy, not total phase energy',
            'legacy performance_efficiency_range is a signed total-thrust-ratio alias',
            'scenario interval has neither validated physical-bound nor statistical-coverage qualification')}


def calculate_two_phase_band(
    *, ideal_thrust_n: float, ideal_specific_impulse_s: float, condensed_mass_fraction: float,
    minimum_particle_coupling_efficiency: float, maximum_particle_coupling_efficiency: float,
    reference_pressure_thrust_n: float = 0., coupling_definition: str = _VELOCITY,
) -> TwoPhasePerformanceBand:
    """Calculate a conditional exit-state range with separated pressure thrust."""
    source = TwoPhaseInputs(**locals())
    return TwoPhasePerformanceBand(**_band_values(source))
