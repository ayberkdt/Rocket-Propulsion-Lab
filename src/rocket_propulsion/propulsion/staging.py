"""Owned staging inputs and replayed conditional finite-grid optima."""
from __future__ import annotations
from collections.abc import Mapping, Set
from dataclasses import dataclass, fields, replace
from fractions import Fraction
from functools import lru_cache
from itertools import islice

from rocket_propulsion.core.errors import DomainError
from rocket_propulsion.core.validation import real, integer, text, _ordered_items
from .burns.numeric_outputs import _finite_rational
from .staging_numerics import (_allocate, _exports, _grid_difference, _STAGING_POLICY,
                               _feasible, _state, _Intervals, _positive_enclosed_output)

_SCOPE = 'global optimum on the declared finite delta-v grid under constant-Isp structural-fraction assumptions'
_TIE = 'lower-stage-first for exactly equal next-increment costs'


@dataclass(frozen=True, slots=True)
class StageDefinition:
    """Dry/(dry+propellant) fraction; payload and upper stages are excluded."""
    name: str
    specific_impulse_s: float
    structural_fraction: float

    def __post_init__(self):
        object.__setattr__(self, 'name', text(self.name, 'Stage name'))
        object.__setattr__(self, 'specific_impulse_s', real(self.specific_impulse_s, 'Stage Isp'))
        object.__setattr__(self, 'structural_fraction', real(self.structural_fraction, 'Structural fraction'))
        if self.specific_impulse_s <= 0 or not 0 < self.structural_fraction < 1:
            raise DomainError('Stage Isp must be positive and structural fraction in (0,1).', field='Stage definition')


def _stage(value):
    if isinstance(value, Mapping):
        try:
            return StageDefinition(**dict(value))
        except (TypeError, ValueError, KeyError, IndexError, OverflowError) as error:
            raise DomainError('Malformed stage mapping.', field='Stage definition') from error
    if isinstance(value, StageDefinition):
        return StageDefinition(value.name, value.specific_impulse_s, value.structural_fraction)
    raise DomainError('Stage definition record required.', field='Stage definition')


@dataclass(frozen=True, slots=True)
class StagingInputs:
    stages: tuple[StageDefinition, ...]
    target_delta_v_m_s: float
    payload_mass_kg: float
    resolution: int = 240

    def __post_init__(self):
        value = self.stages
        if isinstance(value, (str, bytes, bytearray, memoryview, Mapping, Set)):
            raise DomainError('Stages must be ordered records.', field='stages')
        try:
            values = tuple(islice(value, 6))
        except (TypeError, ValueError, KeyError, IndexError, OverflowError) as error:
            raise DomainError('Malformed stage iterable.', field='stages') from error
        if not 1 <= len(values) <= 5:
            raise DomainError('Stage count must be between 1 and 5.', field='stages')
        object.__setattr__(self, 'stages', tuple(_stage(x) for x in values))
        for name in ('target_delta_v_m_s', 'payload_mass_kg'):
            value = real(getattr(self, name), name)
            if value <= 0:
                raise DomainError('Target delta-v and payload must be positive.', field=name)
            object.__setattr__(self, name, value)
        value = integer(self.resolution, 'resolution')
        if not 40 <= value <= 800:
            raise DomainError('Resolution must be in [40,800].', field='resolution')
        object.__setattr__(self, 'resolution', value)


def _inputs(value):
    if isinstance(value, Mapping):
        try:
            return StagingInputs(**dict(value))
        except (TypeError, ValueError, KeyError, IndexError, OverflowError) as error:
            raise DomainError('Malformed retained staging inputs.', field='Staging inputs') from error
    if isinstance(value, StagingInputs):
        return StagingInputs(value.stages, value.target_delta_v_m_s, value.payload_mass_kg, value.resolution)
    raise DomainError('Retained staging inputs are required.', field='Staging inputs')


def _units(value, size):
    return tuple(integer(x, 'Allocation units') for x in _ordered_items(value, 'Allocation units', size=size))


def _literal_fields(value, names):
    for name in names:
        item = getattr(value, name)
        object.__setattr__(value, name, None if item is None and name == 'active_constraint' else text(item, name))


@dataclass(frozen=True, slots=True)
class StageSolution:
    name: str
    delta_v_m_s: float
    mass_ratio: float
    initial_mass_kg: float
    propellant_mass_kg: float
    dry_mass_kg: float
    upper_stack_mass_kg: float
    active_constraint: str | None
    inputs: StagingInputs | None = None
    stage_index: int = 0

    def __post_init__(self):
        source = _inputs(self.inputs); object.__setattr__(self, 'inputs', source)
        index = integer(self.stage_index, 'Stage index'); object.__setattr__(self, 'stage_index', index)
        if not 0 <= index < len(source.stages):
            raise DomainError('Stage index outside retained model.', field='Stage index')
        _literal_fields(self, ('name', 'active_constraint'))
        for name in ('delta_v_m_s', 'mass_ratio', 'initial_mass_kg', 'propellant_mass_kg', 'dry_mass_kg', 'upper_stack_mass_kg'):
            object.__setattr__(self, name, real(getattr(self, name), name))
        expected = _solve(source)['stage_values'][index]
        if any(getattr(self, name) != item for name, item in expected.items()):
            raise DomainError('Stage result does not replay from retained inputs.', field='Stage solution')


@dataclass(frozen=True, slots=True)
class StagingOptimizationResult:
    target_delta_v_m_s: float
    achieved_delta_v_m_s: float
    payload_mass_kg: float
    initial_mass_kg: float
    payload_fraction: float
    stages: tuple[StageSolution, ...]
    discretization_step_m_s: float
    convergence_delta_initial_mass_relative: float | None
    active_constraints: tuple[str, ...]
    feasible: bool
    inputs: StagingInputs | None = None
    allocation_units: tuple[int, ...] = ()
    refined_allocation_units: tuple[int, ...] = ()
    refined_resolution: int = 0
    refinement_performed: bool = False
    refinement_grid_nested: bool = False
    continuous_optimum_qualified: bool = False
    optimization_scope: str = _SCOPE
    tie_break_policy: str = _TIE
    numerical_policy: str = _STAGING_POLICY

    def __post_init__(self):
        source = _inputs(self.inputs); object.__setattr__(self, 'inputs', source)
        for name in ('target_delta_v_m_s', 'achieved_delta_v_m_s', 'payload_mass_kg', 'initial_mass_kg', 'payload_fraction', 'discretization_step_m_s'):
            object.__setattr__(self, name, real(getattr(self, name), name))
        if self.convergence_delta_initial_mass_relative is not None:
            object.__setattr__(self, 'convergence_delta_initial_mass_relative',
                real(self.convergence_delta_initial_mass_relative, 'Grid difference'))
        for name in ('feasible', 'refinement_performed', 'refinement_grid_nested', 'continuous_optimum_qualified'):
            if type(getattr(self, name)) is not bool:
                raise DomainError('Qualification and state flags must be booleans.', field=name)
        object.__setattr__(self, 'refined_resolution', integer(self.refined_resolution, 'Refined resolution'))
        _literal_fields(self, ('optimization_scope', 'tie_break_policy', 'numerical_policy'))
        for name in ('allocation_units', 'refined_allocation_units'):
            object.__setattr__(self, name, _units(getattr(self, name), len(source.stages)))
        constraints = tuple(text(x, 'Constraint') for x in _ordered_items(
            self.active_constraints, 'Constraints', size=len(_solve(source)['active_constraints'])))
        object.__setattr__(self, 'active_constraints', constraints)
        records = []
        for value in _ordered_items(self.stages, 'Stage solutions', size=len(source.stages)):
            if isinstance(value, Mapping):
                try:
                    value = StageSolution(**dict(value))
                except (TypeError, ValueError, KeyError, IndexError, OverflowError) as error:
                    raise DomainError('Malformed stage solution mapping.', field='Stage solutions') from error
            if not isinstance(value, StageSolution):
                raise DomainError('Typed stage solution required.', field='Stage solutions')
            records.append(StageSolution(**{f.name: getattr(value, f.name) for f in fields(StageSolution)}))
        object.__setattr__(self, 'stages', tuple(records))
        expected = _solve(source)
        if (any(getattr(self, name) != item for name, item in expected.items() if name != 'stage_values')
                or any(s.inputs != source or s.stage_index != i for i, s in enumerate(records))
                or self.optimization_scope != _SCOPE or self.tie_break_policy != _TIE or self.numerical_policy != _STAGING_POLICY):
            raise DomainError('Staging result does not replay from retained inputs.', field='Staging result')


@lru_cache(maxsize=64)
def _solve(inputs):
    allocations = _allocate(inputs.stages, inputs.target_delta_v_m_s, inputs.resolution)
    initial, fraction, numeric = _exports(inputs, allocations)
    step = Fraction(inputs.target_delta_v_m_s)/inputs.resolution
    stage_values = tuple(dict(value, name=stage.name,
        delta_v_m_s=_finite_rational(step*count, 'Stage delta-v'),
        active_constraint='zero delta-v allocation' if count == 0 else None)
        for stage, count, value in zip(inputs.stages, allocations, numeric))
    refined_resolution = min(800, inputs.resolution*2)
    performed = refined_resolution != inputs.resolution
    refined = replace(inputs, resolution=refined_resolution)
    refined_counts = _allocate(refined.stages, refined.target_delta_v_m_s, refined.resolution)
    return dict(target_delta_v_m_s=inputs.target_delta_v_m_s, achieved_delta_v_m_s=inputs.target_delta_v_m_s,
        payload_mass_kg=inputs.payload_mass_kg, initial_mass_kg=initial, payload_fraction=fraction,
        discretization_step_m_s=_finite_rational(step, 'Delta-v step'),
        convergence_delta_initial_mass_relative=_grid_difference(inputs, allocations, refined, refined_counts) if performed else None,
        active_constraints=tuple(f"{s['name']}: {s['active_constraint']}" for s in stage_values if s['active_constraint']),
        feasible=True, allocation_units=allocations, refined_allocation_units=refined_counts,
        refined_resolution=refined_resolution, refinement_performed=performed,
        refinement_grid_nested=performed and refined_resolution % inputs.resolution == 0,
        continuous_optimum_qualified=False, stage_values=stage_values)


def _mass_factor(stage, delta_v_m_s):
    """Compatibility helper; exported factor is never used to order allocations."""
    stage = _stage(stage); dv = real(delta_v_m_s, 'Stage delta-v')
    if dv < 0:
        raise DomainError('Stage delta-v must be nonnegative.', field='Stage delta-v')
    if not _feasible(stage, Fraction(dv)):
        return float('inf')
    def evaluate(p):
        math = _Intervals(p)
        denominator = _state(stage, Fraction(dv), p)[1]
        return math.div(math.rational(1-Fraction(stage.structural_fraction)), denominator)
    return _positive_enclosed_output(evaluate, 'Stage mass factor')


def optimize_staging(*, stages, target_delta_v_m_s, payload_mass_kg, resolution=240):
    """Minimize conditional initial mass on the requested finite allocation grid."""
    inputs = StagingInputs(stages, target_delta_v_m_s, payload_mass_kg, resolution)
    values = dict(_solve(inputs)); stage_values = values.pop('stage_values')
    solutions = tuple(StageSolution(**value, inputs=inputs, stage_index=i) for i, value in enumerate(stage_values))
    return StagingOptimizationResult(**values, stages=solutions, inputs=inputs)
