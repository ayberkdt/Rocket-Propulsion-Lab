"""First off-axis intersection with an owned second-order axis expansion.

The C+ parent is a replayed axis-boundary result, including the quadratic
history that determines its finite cylindrical source. No epsilon radius or
free derivative is substituted for axis regularity. The radial velocity uses
the regular Taylor expansion, while C- compatibility and both characteristic
geometries are iterated. C+ endpoint-averaged compatibility is exposed as an
unenforced defect, not hidden inside a claimed two-family algebraic solve.
Planar flow retains the regular two-family interior process. This is a
constant-gamma local process, not a sonic starting solution or a complete field.
The governing equations and regularity derivation are documented in
MOC_UNIT_PROCESS_REVIEW.md and MOC_BOUNDARY_REVIEW.md.
"""
from dataclasses import dataclass
from fractions import Fraction as F
from functools import lru_cache

from rocket_propulsion.core.errors import ConvergenceError, DomainError
from rocket_propulsion.core.validation import integer, real, text
from rocket_propulsion.propulsion.burns.numeric_outputs import _finite_rational

from .moc_boundary import MocBoundaryResult, _axis_derivative_terms, _axis_rows, _optional_source
from .moc_unit_process import (
    _TOLERANCE,
    MocFlowPoint,
    _coefficients,
    _coordinate_cell_radius,
    _intersection,
    _owned,
    _point_from_velocity,
    _upper_residual,
    solve_moc_interior_point,
)

_POLICY='axis_departure_second_order_regularity_and_averaged_minus_v1'


@dataclass(frozen=True,slots=True)
class MocAxisDepartureInputs:
    axis_result:MocBoundaryResult
    minus_parent:MocFlowPoint
    max_iterations:int=128

    def __post_init__(self):
        axis=_owned(self.axis_result,MocBoundaryResult,'Axis boundary result')
        if axis.boundary_kind!='axis':raise DomainError('An axis result is required for the C+ parent.',field='axis_result')
        object.__setattr__(self,'axis_result',axis)
        minus=_owned(self.minus_parent,MocFlowPoint,'Off-axis C- parent')
        if minus.radius_m<=0:raise DomainError('The C- parent must be above the axis.',field='minus_parent')
        object.__setattr__(self,'minus_parent',minus)
        count=integer(self.max_iterations,'max_iterations')
        if not 1<=count<=256:raise DomainError('Iteration budget must be between 1 and 256.',field='max_iterations')
        object.__setattr__(self,'max_iterations',count)


def _parents(inputs):
    axis=inputs.axis_result;g=axis.inputs.gamma
    base=_coefficients(axis.point,g,1,False)
    source=_axis_rows(axis.inputs,axis.point)[-1]
    plus=(*base[:3],source,*base[4:])
    minus=_coefficients(inputs.minus_parent,g,-1,axis.inputs.axisymmetric)
    return plus,minus


def _axis_expansion(inputs):
    """v=r*(w+w_x*dx)+O(r*dx**2+r**3), from continuity and u_r=v_x.

    In cylindrical flow, w=(M**2-1)*u_x/2.
    Its derivative is w_x=c*u_xx+c_u*u_x**2. Both u derivatives
    come from the same recorded quadratic axis polynomial.
    """
    axis=inputs.axis_result;p=axis.point;g=F(axis.inputs.gamma);m=F(p.mach)
    u=_coefficients(p,axis.inputs.gamma,1,False)[4]
    constant,coefficient=_axis_derivative_terms(axis.inputs,F(p.x_m));first=constant+coefficient*u
    a,b=axis.inputs.axis_history
    u0=_coefficients(a,axis.inputs.gamma,1,False)[4];u1=_coefficients(b,axis.inputs.gamma,1,False)[4]
    second=2*((u-u1)/(F(p.x_m)-F(b.x_m))-(u1-u0)/(F(b.x_m)-F(a.x_m)))/(F(p.x_m)-F(a.x_m))
    a2=(g+1)/(2+(g-1)*m*m);c=(m*m-1)/2;c_u=(g+1)*u/(2*a2*a2)
    return c*first,c*second+c_u*first*first


def _rows(inputs,point):
    plus,minus=_parents(inputs);g=inputs.axis_result.inputs.gamma;mode=inputs.axis_result.inputs.axisymmetric
    at_plus=_coefficients(point,g,1,mode);at_minus=_coefficients(point,g,-1,mode)
    averaged=tuple(tuple((a+b)/2 for a,b in zip(parent[:4],candidate[:4]))
        for parent,candidate in ((plus,at_plus),(minus,at_minus)))
    return plus,minus,at_plus,averaged


def _metrics(inputs,point):
    plus,minus,at,averaged=_rows(inputs,point)
    a,b=inputs.axis_result.point,inputs.minus_parent
    scale=max(abs(F(a.x_m)-F(b.x_m)),F(b.radius_m))
    velocity=max(abs(v) for row in (plus,minus) for v in row[4:])
    cylindrical=inputs.axis_result.inputs.axisymmetric
    w,wx=_axis_expansion(inputs) if cylindrical else (F(0),F(0))
    regularity=abs(at[5]-F(point.radius_m)*(w+wx*(F(point.x_m)-F(a.x_m))))/velocity if cylindrical else F(0)
    raw=[regularity];adjusted=[regularity];allowances=[];plus_defect=F(0)
    for index,(parent,pc,(slope,aa,bb,source)) in enumerate(zip((a,b),(plus,minus),averaged)):
        dx=F(point.x_m)-F(parent.x_m);dy=F(point.radius_m)-F(parent.radius_m)
        geometry=abs(dy-slope*dx)/scale
        cell=(_coordinate_cell_radius(point.radius_m)+abs(slope)*_coordinate_cell_radius(point.x_m))/scale
        normalizer=max(abs(aa),abs(bb))*velocity
        if normalizer==0:raise ConvergenceError('Axis-departure compatibility is singular.',code='moc_departure_singular')
        compatibility=abs(aa*(at[4]-pc[4])+bb*(at[5]-pc[5])-source*dx)/normalizer
        raw.append(geometry);adjusted.append(max(F(0),geometry-cell));allowances.append(cell)
        if index==0:plus_defect=compatibility
        if index==1 or not cylindrical:raw.append(compatibility);adjusted.append(compatibility)
    return max(raw),max(adjusted),max(allowances),plus_defect


def _private_topology(inputs,x,y):
    if x<=max(F(inputs.axis_result.point.x_m),F(inputs.minus_parent.x_m)):
        raise ConvergenceError('Axis-departure parents do not intersect downstream.',code='moc_departure_noncausal')
    if y<=0:raise ConvergenceError('Axis departure must enter the upper half-plane.',code='moc_departure_axis_boundary')


def _stored_topology(inputs,point):
    if point.x_m<=max(inputs.axis_result.point.x_m,inputs.minus_parent.x_m) or point.radius_m<=0:
        raise ConvergenceError('Stored axis departure cannot resolve marching topology.',code='moc_departure_unresolved_topology')


def _pack(point,count,metrics,kind,source):
    raw,adjusted,cell,plus_defect=metrics;value,status=_optional_source(source)
    return {'point': point,'iterations': count,'equation_residual_max': _upper_residual(raw),
        'rounding_adjusted_residual_max': _upper_residual(adjusted),'coordinate_rounding_allowance_max': _upper_residual(cell),
        'representation_limited': raw>_TOLERANCE,'solution_kind': kind,'axis_source_term_per_m': value,'axis_source_term_status': status,
        'averaged_plus_compatibility_defect': _upper_residual(plus_defect)}


@lru_cache(maxsize=16)
def _solve(inputs):
    a,b=inputs.axis_result.point,inputs.minus_parent;g=inputs.axis_result.inputs.gamma
    plus,minus=_parents(inputs);x,y=_intersection(a,b,plus[0],minus[0]);_private_topology(inputs,x,y)
    if a.mach==b.mach and b.flow_angle_deg==0 and all(p.mach==a.mach for p in inputs.axis_result.inputs.axis_history):
        point=MocFlowPoint(_finite_rational(x,'Departure x'),_finite_rational(y,'Departure radius'),a.mach,0.)
        _stored_topology(inputs,point);metrics=_metrics(inputs,point)
        if metrics[1]<=_TOLERANCE:return _pack(point,1,metrics,'constant_axis_departure',plus[3])
        raise ConvergenceError('Constant axis departure cannot be represented.',code='moc_departure_unresolved')
    if not inputs.axis_result.inputs.axisymmetric:
        result=solve_moc_interior_point(plus_parent=a,minus_parent=b,gamma=g,axisymmetric=False,max_iterations=inputs.max_iterations)
        return _pack(result.point,result.iterations,_metrics(inputs,result.point),'iterated_planar_axis_departure',plus[3])
    point=MocFlowPoint(_finite_rational(x,'Departure x'),_finite_rational(y,'Departure radius'),
        _finite_rational((F(a.mach)+F(b.mach))/2,'Mean Mach'),b.flow_angle_deg/2)
    _stored_topology(inputs,point)
    w,wx=_axis_expansion(inputs)
    for count in range(1,inputs.max_iterations+1):
        plus,minus,_,((sp,_ap,_bp,_cp),(sm,am,bm,cm))=_rows(inputs,point)
        x,y=_intersection(a,b,sp,sm);_private_topology(inputs,x,y)
        rm=am*minus[4]+bm*minus[5]+cm*(x-F(b.x_m))
        if am==0:raise ConvergenceError('Axis-departure system is singular.',code='moc_departure_singular')
        v=y*(w+wx*(x-F(a.x_m)))
        u=(rm-bm*v)/am
        point=_point_from_velocity(x,y,u,v,g);_stored_topology(inputs,point);metrics=_metrics(inputs,point)
        if metrics[1]<=_TOLERANCE:return _pack(point,count,metrics,'iterated_axis_departure',plus[3])
    raise ConvergenceError('Axis-departure iteration budget exhausted.',code='moc_departure_not_converged')


@dataclass(frozen=True,slots=True)
class MocAxisDepartureResult:
    point:MocFlowPoint
    iterations:int
    equation_residual_max:float
    rounding_adjusted_residual_max:float
    coordinate_rounding_allowance_max:float
    representation_limited:bool
    solution_kind:str
    axis_source_term_per_m:float|None
    axis_source_term_status:str
    averaged_plus_compatibility_defect:float
    inputs:MocAxisDepartureInputs
    numerical_policy:str=_POLICY
    physical_model_qualified:bool=False
    full_nozzle_solution:bool=False
    discretization_error_bound:None=None

    def __post_init__(self):
        source=_owned(self.inputs,MocAxisDepartureInputs,'Axis-departure inputs');object.__setattr__(self,'inputs',source)
        point=_owned(self.point,MocFlowPoint,'Axis-departure point');object.__setattr__(self,'point',point);_stored_topology(source,point)
        object.__setattr__(self,'iterations',integer(self.iterations,'iterations'))
        for name in ('equation_residual_max','rounding_adjusted_residual_max','coordinate_rounding_allowance_max','averaged_plus_compatibility_defect'):
            object.__setattr__(self,name,real(getattr(self,name),name))
        if self.axis_source_term_per_m is not None:object.__setattr__(self,'axis_source_term_per_m',real(self.axis_source_term_per_m,'axis_source_term_per_m'))
        for name in ('solution_kind','axis_source_term_status','numerical_policy'):object.__setattr__(self,name,text(getattr(self,name),name))
        for name in ('representation_limited','physical_model_qualified','full_nozzle_solution'):
            if type(getattr(self,name)) is not bool:raise DomainError('Departure flags must be actual booleans.',field=name)
        if self.physical_model_qualified or self.full_nozzle_solution or self.discretization_error_bound is not None or self.numerical_policy!=_POLICY:
            raise DomainError('Axis departure cannot claim full-nozzle or hardware qualification.',field='Axis departure')
        if any(getattr(self,name)!=value for name,value in _solve(source).items()):
            raise DomainError('Axis-departure result does not replay from owned inputs.',field='Axis departure')


def solve_moc_axis_departure(*,axis_result,minus_parent,max_iterations=128):
    """March from a replayed axis using second-order radial regularity and C-.

    Gamma and planar/cylindrical mode come from the owned axis result. Its
    finite source and Taylor coefficients are recomputed privately from history,
    including when that source cannot be exported as a finite diagnostic.
    In cylindrical flow, the averaged C+ defect is visible and not an
    iteration condition; the radial Taylor relation replaces that equation.
    Planar flow retains both compatibility equations.
    This process supplies neither sonic initialization nor a global error bound.
    """
    inputs=MocAxisDepartureInputs(axis_result,minus_parent,max_iterations)
    return MocAxisDepartureResult(**_solve(inputs),inputs=inputs)
