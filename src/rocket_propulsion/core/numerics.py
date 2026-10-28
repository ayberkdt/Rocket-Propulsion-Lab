"""Finite sampled-callback bisection with replayable acceptance evidence.

Continuity and monotonicity are caller preconditions, not proved by samples.
The residual rule is abs(f(x)-target) <= tol*max(1,abs(target)), evaluated on
exact normalized binary values. No root/physical qualification is conferred.
"""
from collections.abc import Callable, Mapping
from collections.abc import Set as AbstractSet
from dataclasses import dataclass
from fractions import Fraction as F
from itertools import islice

from .errors import ConvergenceError, DomainError, InputError
from .validation import _ordered_items, coerce_fields, integer, real, text

EVIDENCE_KIND = 'finite_binary_callback_bisection_trace_v1'


def _controls(target, lower, upper, tolerance, max_iterations):
    target=real(target,'target');lower=real(lower,'lower');upper=real(upper,'upper')
    tolerance=real(tolerance,'tolerance');max_iterations=integer(max_iterations,'max_iterations')
    if tolerance <= 0:
        raise InputError('Tolerance must be positive.',field='tolerance')
    if not 1 <= max_iterations <= 4096:
        raise InputError('Maximum iterations must be in 1..4096.',field='max_iterations')
    if lower > upper:
        raise InputError('Lower bound must not exceed upper bound.',code='invalid_bracket',
                         details={'lower':lower,'upper':upper})
    return target,lower,upper,tolerance,max_iterations


def _midpoint(lower, upper):
    return float((F(lower)+F(upper))/2)


def _residual(value, target):
    return abs(F(value)-F(target))


def _accepted(value, target, tolerance):
    return _residual(value,target) <= F(tolerance)*max(F(1),abs(F(target)))


def _pack(value, observed, iterations, lower, upper, low_value, high_value, target):
    residual=_residual(observed,target)
    return {'value':value,'value_function_value':observed,'iterations':iterations,
            'lower_bound':lower,'upper_bound':upper,
            'lower_function_value':low_value,'upper_function_value':high_value,
            'residual':float(residual),'residual_ratio':(str(residual.numerator),str(residual.denominator)),
            'scaled_residual':float(residual/max(F(1),abs(F(target)))),
            'termination':'exact_sample' if residual==0 else 'residual_tolerance'}


def _failure(code, iterations, lower, upper, low_value, high_value, target, tolerance, budget):
    value,observed=min(((lower,low_value),(upper,high_value)),key=lambda row:_residual(row[1],target))
    details=_pack(value,observed,iterations,lower,upper,low_value,high_value,target)
    details.update(target=target,tolerance=tolerance,max_iterations=budget,method='bisection',
                   root_enclosure_qualified=False,physical_model_qualified=False,
                   termination=code,acceptance_satisfied=False)
    raise ConvergenceError('Bisection cannot satisfy the requested residual criterion.',code=code,details=details)


def _run(evaluate, target, lower, upper, tolerance, budget):
    low_value=evaluate(lower)
    high_value=low_value if lower==upper else evaluate(upper)
    if not min(low_value,high_value) <= target <= max(low_value,high_value):
        raise DomainError('Target is outside the actual sampled bracket.',code='target_not_bracketed',
                          details={'target':target,'lower_bound':lower,'upper_bound':upper,
                                   'lower_function_value':low_value,'upper_function_value':high_value})
    increasing=high_value>=low_value
    value,observed=min(((lower,low_value),(upper,high_value)),key=lambda row:_residual(row[1],target))
    if _accepted(observed,target,tolerance):
        return _pack(value,observed,0,lower,upper,low_value,high_value,target)
    for iteration in range(1,budget+1):
        midpoint=_midpoint(lower,upper)
        if midpoint==lower or midpoint==upper:
            _failure('root_float_resolution',iteration-1,lower,upper,low_value,high_value,target,tolerance,budget)
        value=evaluate(midpoint)
        if not min(low_value,high_value) <= value <= max(low_value,high_value):
            raise DomainError('Callback observation violates sampled monotonic ordering.',
                              code='callback_not_monotonic',field='callback',
                              details={'x':midpoint,'value':value,'lower_function_value':low_value,
                                       'upper_function_value':high_value})
        if value != target:
            if (value < target)==increasing:
                lower,low_value=midpoint,value
            else:
                upper,high_value=midpoint,value
        if _accepted(value,target,tolerance):
            return _pack(midpoint,value,iteration,lower,upper,low_value,high_value,target)
    _failure('non_convergence',budget,lower,upper,low_value,high_value,target,tolerance,budget)


@dataclass(frozen=True,slots=True)
class RootEvaluation:
    """One finite normalized binary callback observation, without truth attestation."""
    x: float
    value: float

    def __post_init__(self):
        coerce_fields(self,'Root evaluation')


def _owned_trace(value, budget):
    if isinstance(value,(str,bytes,bytearray,memoryview,Mapping,AbstractSet)):
        raise DomainError('Root trace must be an ordered record iterable.',field='evaluations')
    try: rows=tuple(islice(iter(value),budget+3))
    except (TypeError,ValueError,KeyError,IndexError,OverflowError) as error:
        raise DomainError('Invalid root trace iterable.',field='evaluations') from error
    if not rows or len(rows)>budget+2 or any(not isinstance(row,RootEvaluation) for row in rows):
        raise DomainError('Root trace has invalid records or exceeds its budget.',field='evaluations')
    return tuple(RootEvaluation(row.x,row.value) for row in rows)


@dataclass(frozen=True,slots=True)
class RootResult:
    """Owned sampled observations and replay-checked residual acceptance.

    Bounds preserve sampled endpoint ordering. They are not certified
    enclosures of an unrounded mathematical or physical model's root.
    """
    value: float
    iterations: int
    residual: float
    lower_bound: float
    upper_bound: float
    tolerance: float
    target: float
    initial_lower_bound: float
    initial_upper_bound: float
    max_iterations: int
    value_function_value: float
    lower_function_value: float
    upper_function_value: float
    scaled_residual: float
    residual_ratio: tuple[str,str]
    termination: str
    evaluations: tuple[RootEvaluation,...]
    method: str = 'bisection'
    evidence_kind: str = EVIDENCE_KIND
    root_enclosure_qualified: bool = False
    physical_model_qualified: bool = False
    coordinate_error_bound: float | None = None

    def __post_init__(self):
        coerce_fields(self,'Root result')
        for name in ('root_enclosure_qualified','physical_model_qualified'):
            if type(getattr(self,name)) is not bool or getattr(self,name):
                raise DomainError('Sampled root evidence does not establish qualification.',field=name)
        if self.coordinate_error_bound is not None or self.method!='bisection' or self.evidence_kind!=EVIDENCE_KIND:
            raise DomainError('Unsupported root evidence or coordinate-bound claim.',field='evidence_kind')
        controls=_controls(self.target,self.initial_lower_bound,self.initial_upper_bound,self.tolerance,self.max_iterations)
        rows=_owned_trace(self.evaluations,self.max_iterations)
        ratio=tuple(text(x,'residual_ratio') for x in _ordered_items(self.residual_ratio,'residual_ratio',size=2))
        index=0
        def evaluate(x):
            nonlocal index
            if index>=len(rows) or rows[index].x!=x:
                raise DomainError('Root trace does not follow the declared bisection.',field='evaluations')
            value=rows[index].value;index+=1
            return value
        expected=_run(evaluate,*controls)
        if index!=len(rows) or any(getattr(self,name)!=value for name,value in expected.items() if name!='residual_ratio') or ratio!=expected['residual_ratio']:
            raise DomainError('Root result contradicts its owned trace.',field='result')
        object.__setattr__(self,'evaluations',rows)
        object.__setattr__(self,'residual_ratio',ratio)


def solve_bisect_monotonic(function: Callable[[float],float], target: float, lower: float, upper: float,
                           *, tolerance: float=1e-12,max_iterations: int=200) -> RootResult:
    """Solve a caller-supplied continuous monotonic relation with finite samples.

    Preserve exact sampled bracketing and the existing residual rule. An
    adjacent-float bracket that fails that rule is a typed resolution failure,
    not permission to loosen tolerance. Callback application errors propagate.
    """
    controls=_controls(target,lower,upper,tolerance,max_iterations)
    if not callable(function):
        raise InputError('Root function must be callable.',field='function')
    rows=[]
    def evaluate(x):
        value=real(function(x),'callback_value')
        rows.append(RootEvaluation(x,value))
        return value
    values=_run(evaluate,*controls)
    target,lower,upper,tolerance,budget=controls
    return RootResult(**values,target=target,initial_lower_bound=lower,initial_upper_bound=upper,
                      tolerance=tolerance,max_iterations=budget,evaluations=tuple(rows))


def bisect_monotonic(function: Callable[[float],float],target: float,lower: float,upper: float,
                     *,tolerance: float=1e-12,max_iterations: int=200) -> float:
    """Return the sampled-callback solution accepted by solve_bisect_monotonic."""
    return solve_bisect_monotonic(function,target,lower,upper,tolerance=tolerance,max_iterations=max_iterations).value
