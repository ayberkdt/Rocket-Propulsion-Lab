"""Constant-gamma characteristic axis reflection and prescribed straight-wall slip.

Axis regularity follows cylindrical continuity: v/r -> v_r and
v_r=(M**2-1)*u_x/2 at a smooth symmetry axis. A quadratic through two previous
axis velocities and the new velocity supplies u_x; it is a discretization,
not a prescribed free gradient or a transonic starting solution. The wall
process solves slip on one finite straight segment. Neither constructs a
complete nozzle or qualifies an attached/isentropic hardware flow.

Governing Cartesian equations: https://ntrs.nasa.gov/citations/19940019183
Axis/streamline continuity: https://ntrs.nasa.gov/citations/19690005234
Characteristic averaging context: https://ntrs.nasa.gov/citations/19980232191
"""
from dataclasses import dataclass, fields
from fractions import Fraction as F
from functools import lru_cache
from math import atan2, degrees

from rocket_propulsion.core.errors import DomainError, ConvergenceError
from rocket_propulsion.core.validation import _ordered_items, integer, real, text
from rocket_propulsion.propulsion.burns.numeric_outputs import _finite_rational
from .moc_unit_process import (MocFlowPoint, _owned, _coefficients, _point_from_velocity,
    _coordinate_cell_radius, _upper_residual, _TOLERANCE)

_POLICY='averaged_cartesian_axis_quadratic_regularity_or_finite_straight_wall_v1'


def _common(inputs):
    g=real(inputs.gamma,'gamma')
    if g<=1:raise DomainError('Gamma must exceed one.',field='gamma')
    object.__setattr__(inputs,'gamma',g)
    if type(inputs.axisymmetric) is not bool:raise DomainError('Axisymmetric selector must be boolean.',field='axisymmetric')
    count=integer(inputs.max_iterations,'max_iterations')
    if not 1<=count<=256:raise DomainError('Iteration budget must be between 1 and 256.',field='max_iterations')
    object.__setattr__(inputs,'max_iterations',count)
    parent=_owned(inputs.incoming,MocFlowPoint,'Incoming characteristic')
    if parent.radius_m<=0:raise DomainError('Incoming point must be above the axis.',field='radius_m')
    object.__setattr__(inputs,'incoming',parent)


@dataclass(frozen=True,slots=True)
class MocAxisInputs:
    incoming:MocFlowPoint
    axis_history:tuple[MocFlowPoint,MocFlowPoint]
    gamma:float=1.4
    axisymmetric:bool=True
    max_iterations:int=128

    def __post_init__(self):
        _common(self)
        history=tuple(_owned(p,MocFlowPoint,'Axis history') for p in _ordered_items(self.axis_history,'Axis history',size=2))
        if any(p.radius_m!=0 or p.flow_angle_deg!=0 for p in history):
            raise DomainError('Axis history requires r=0 and zero flow angle.',field='axis_history')
        if history[1].x_m<=history[0].x_m:raise DomainError('Axis history must increase strictly in x.',field='axis_history')
        object.__setattr__(self,'axis_history',history)


@dataclass(frozen=True,slots=True)
class MocWallSegment:
    x0_m:float
    radius0_m:float
    x1_m:float
    radius1_m:float

    def __post_init__(self):
        for f in fields(self):object.__setattr__(self,f.name,real(getattr(self,f.name),f.name))
        if self.x1_m<=self.x0_m or min(self.radius0_m,self.radius1_m)<=0:
            raise DomainError('An upper wall segment requires increasing x and positive radii.',field='wall_segment')


@dataclass(frozen=True,slots=True)
class MocWallInputs:
    incoming:MocFlowPoint
    wall_segment:MocWallSegment
    gamma:float=1.4
    axisymmetric:bool=True
    max_iterations:int=128

    def __post_init__(self):
        _common(self)
        object.__setattr__(self,'wall_segment',_owned(self.wall_segment,MocWallSegment,'Wall segment'))


def _axis_derivative_terms(inputs,x):
    """Lagrange derivative at the new point: u_x=constant+coefficient*u_new."""
    a,b=inputs.axis_history;x0,x1=F(a.x_m),F(b.x_m)
    if x<=x1:raise ConvergenceError('Axis intersection does not follow its history.',code='moc_axis_noncausal')
    u0=_coefficients(a,inputs.gamma,-1,False)[4]
    u1=_coefficients(b,inputs.gamma,-1,False)[4]
    c0=(x-x1)/((x0-x1)*(x0-x))
    c1=(x-x0)/((x1-x0)*(x1-x))
    cp=1/(x-x0)+1/(x-x1)
    return c0*u0+c1*u1,cp


def _axis_rows(inputs,point):
    parent=_coefficients(inputs.incoming,inputs.gamma,-1,inputs.axisymmetric)
    at_point=_coefficients(point,inputs.gamma,-1,False)
    m=F(point.mach);g=F(inputs.gamma)
    a2=(g+1)/(2+(g-1)*m*m)
    k=a2*(m*m-1)/2 if inputs.axisymmetric else F(0)
    constant,coefficient=_axis_derivative_terms(inputs,F(point.x_m))
    source=k*(constant+coefficient*at_point[4])
    return parent,at_point,tuple((a+b)/2 for a,b in zip(parent[:3],at_point[:3])),k,source


def _axis_metrics(inputs,point):
    parent,current,(slope,aa,bb),_,source=_axis_rows(inputs,point)
    p=inputs.incoming;a,b=inputs.axis_history
    dx=F(point.x_m)-F(p.x_m)
    scale=max(F(p.radius_m),F(point.x_m)-F(b.x_m),F(b.x_m)-F(a.x_m))
    geometry=abs(F(p.radius_m)+slope*dx)/scale
    allowance=abs(slope)*_coordinate_cell_radius(point.x_m)/scale
    normalizer=max(abs(aa),abs(bb))*max(abs(parent[4]),abs(parent[5]),abs(current[4]))
    if normalizer==0:raise ConvergenceError('Axis compatibility equation is singular.',code='moc_axis_singular')
    compatibility=abs(aa*(current[4]-parent[4])-bb*parent[5]-(parent[3]+source)*dx/2)/normalizer
    return max(geometry,compatibility),max(max(F(0),geometry-allowance),compatibility),allowance,source


def _axis_topology(inputs,point):
    if point.radius_m!=0 or point.flow_angle_deg!=0:raise DomainError('Axis point must preserve symmetry exactly.',field='Axis point')
    if point.x_m<=max(inputs.incoming.x_m,inputs.axis_history[-1].x_m):
        raise ConvergenceError('Stored axis point cannot resolve downstream topology.',code='moc_axis_unresolved_topology')


def _optional_source(value):
    try:return _finite_rational(value,'Boundary source term'),'finite'
    except DomainError:return None,'unrepresentable'


def _pack(point,iteration,metrics,kind):
    raw,adjusted,allowance,source=metrics
    source_value,status=_optional_source(source)
    return dict(point=point,iterations=iteration,equation_residual_max=_upper_residual(raw),
        rounding_adjusted_residual_max=_upper_residual(adjusted),coordinate_rounding_allowance_max=_upper_residual(allowance),
        representation_limited=raw>_TOLERANCE,source_term_per_m=source_value,source_term_status=status,solution_kind=kind)


@lru_cache(maxsize=16)
def _solve_axis(inputs):
    p=inputs.incoming
    parent=_coefficients(p,inputs.gamma,-1,inputs.axisymmetric)
    guess=MocFlowPoint(p.x_m,0.,p.mach,0.)
    slope=(parent[0]+_coefficients(guess,inputs.gamma,-1,False)[0])/2
    if slope==0:raise ConvergenceError('Axis characteristic is parallel to the axis.',code='moc_axis_parallel')
    x=F(p.x_m)-F(p.radius_m)/slope
    if x<=max(F(p.x_m),F(inputs.axis_history[-1].x_m)):
        raise ConvergenceError('Axis intersection is not downstream.',code='moc_axis_noncausal')
    point=MocFlowPoint(_finite_rational(x,'Axis x'),0.,p.mach,0.)
    _axis_topology(inputs,point)
    if p.flow_angle_deg==0 and all(a.mach==p.mach for a in inputs.axis_history):
        metrics=_axis_metrics(inputs,point)
        if metrics[1]<=_TOLERANCE:return _pack(point,1,metrics,'uniform_axis_intersection')
        raise ConvergenceError('Uniform axis point cannot be represented.',code='moc_axis_unresolved')
    for iteration in range(1,inputs.max_iterations+1):
        incoming,_,(slope,aa,bb),k,_=_axis_rows(inputs,point)
        if slope==0:raise ConvergenceError('Axis characteristic is parallel.',code='moc_axis_parallel')
        x=F(p.x_m)-F(p.radius_m)/slope;dx=x-F(p.x_m)
        if x<=F(p.x_m):raise ConvergenceError('Axis iteration is upstream.',code='moc_axis_noncausal')
        constant,coefficient=_axis_derivative_terms(inputs,x)
        denominator=aa-dx*k*coefficient/2
        if denominator==0:raise ConvergenceError('Axis regularity system is singular.',code='moc_axis_singular')
        u=(aa*incoming[4]+bb*incoming[5]+dx*(incoming[3]+k*constant)/2)/denominator
        point=_point_from_velocity(x,F(0),u,F(0),inputs.gamma)
        _axis_topology(inputs,point)
        metrics=_axis_metrics(inputs,point)
        if metrics[1]<=_TOLERANCE:return _pack(point,iteration,metrics,'iterated_quadratic_axis_regularity')
    raise ConvergenceError('Axis iteration budget exhausted.',code='moc_axis_not_converged')


def _wall_slope(segment):
    return (F(segment.radius1_m)-F(segment.radius0_m))/(F(segment.x1_m)-F(segment.x0_m))


def _wall_angle(segment):
    dx=F(segment.x1_m)-F(segment.x0_m);dy=F(segment.radius1_m)-F(segment.radius0_m)
    scale=max(abs(dx),abs(dy))
    return degrees(atan2(_finite_rational(dy/scale,'Wall radial direction'),_finite_rational(dx/scale,'Wall axial direction')))


def _wall_position(inputs,slope):
    p=inputs.incoming;w=inputs.wall_segment;wall=_wall_slope(w)
    if slope==wall:raise ConvergenceError('Characteristic and wall are parallel.',code='moc_wall_parallel')
    x=(F(w.radius0_m)-F(p.radius_m)+slope*F(p.x_m)-wall*F(w.x0_m))/(slope-wall)
    y=F(w.radius0_m)+wall*(x-F(w.x0_m))
    if x<=max(F(p.x_m),F(w.x0_m)):raise ConvergenceError('Wall intersection is not downstream.',code='moc_wall_noncausal')
    if x>F(w.x1_m):raise ConvergenceError('Wall intersection lies outside the prescribed segment.',code='moc_wall_outside_segment')
    return x,y


def _wall_topology(inputs,point):
    if point.x_m<=max(inputs.incoming.x_m,inputs.wall_segment.x0_m):
        raise ConvergenceError('Stored wall point cannot resolve downstream topology.',code='moc_wall_unresolved_topology')
    if point.x_m>inputs.wall_segment.x1_m:
        raise ConvergenceError('Stored wall point exceeds the prescribed segment.',code='moc_wall_outside_segment')
    if point.radius_m<=0:raise DomainError('Wall point must remain above the axis.',field='Wall point')


def _wall_rows(inputs,point):
    parent=_coefficients(inputs.incoming,inputs.gamma,1,inputs.axisymmetric)
    current=_coefficients(point,inputs.gamma,1,inputs.axisymmetric)
    return parent,current,tuple((a+b)/2 for a,b in zip(parent[:4],current[:4]))


def _wall_metrics(inputs,point):
    parent,current,(slope,aa,bb,source)=_wall_rows(inputs,point)
    p=inputs.incoming;w=inputs.wall_segment;wall=_wall_slope(w)
    dx=F(point.x_m)-F(p.x_m);dy=F(point.radius_m)-F(p.radius_m)
    scale=max(F(w.x1_m)-F(w.x0_m),abs(F(w.radius1_m)-F(w.radius0_m)),abs(dy))
    g1=abs(dy-slope*dx)/scale
    g2=abs(F(point.radius_m)-F(w.radius0_m)-wall*(F(point.x_m)-F(w.x0_m)))/scale
    a1=(_coordinate_cell_radius(point.radius_m)+abs(slope)*_coordinate_cell_radius(point.x_m))/scale
    a2=(_coordinate_cell_radius(point.radius_m)+abs(wall)*_coordinate_cell_radius(point.x_m))/scale
    velocity=max(abs(v) for row in (parent,current) for v in row[4:])
    normalizer=max(abs(aa),abs(bb))*velocity
    if normalizer==0:raise ConvergenceError('Wall compatibility equation is singular.',code='moc_wall_singular')
    compatibility=abs(aa*(current[4]-parent[4])+bb*(current[5]-parent[5])-source*dx)/normalizer
    slip=abs(current[5]-wall*current[4])/velocity
    return max(g1,g2,compatibility,slip),max(max(F(0),g1-a1),max(F(0),g2-a2),compatibility,slip),max(a1,a2),current[3]


@lru_cache(maxsize=16)
def _solve_wall(inputs):
    p=inputs.incoming;wall=_wall_slope(inputs.wall_segment)
    parent=_coefficients(p,inputs.gamma,1,inputs.axisymmetric)
    if wall==0 and p.flow_angle_deg==0:
        x,y=_wall_position(inputs,parent[0])
        point=MocFlowPoint(_finite_rational(x,'Wall x'),_finite_rational(y,'Wall radius'),p.mach,0.)
        _wall_topology(inputs,point);metrics=_wall_metrics(inputs,point)
        if metrics[1]<=_TOLERANCE:return _pack(point,1,metrics,'uniform_horizontal_wall')
        raise ConvergenceError('Uniform wall point cannot be represented.',code='moc_wall_unresolved')
    angle=_wall_angle(inputs.wall_segment)
    guess=MocFlowPoint(p.x_m,p.radius_m,p.mach,angle)
    slope=(parent[0]+_coefficients(guess,inputs.gamma,1,inputs.axisymmetric)[0])/2
    x,y=_wall_position(inputs,slope)
    point=MocFlowPoint(_finite_rational(x,'Wall x'),_finite_rational(y,'Wall radius'),p.mach,angle)
    _wall_topology(inputs,point)
    for iteration in range(1,inputs.max_iterations+1):
        incoming,_,(slope,aa,bb,source)=_wall_rows(inputs,point)
        x,y=_wall_position(inputs,slope);dx=x-F(p.x_m)
        denominator=aa+bb*wall
        if denominator==0:raise ConvergenceError('Wall slip system is singular.',code='moc_wall_singular')
        u=(aa*incoming[4]+bb*incoming[5]+source*dx)/denominator
        point=_point_from_velocity(x,y,u,wall*u,inputs.gamma)
        _wall_topology(inputs,point);metrics=_wall_metrics(inputs,point)
        if metrics[1]<=_TOLERANCE:return _pack(point,iteration,metrics,'iterated_straight_wall_slip')
    raise ConvergenceError('Wall iteration budget exhausted.',code='moc_wall_not_converged')


@dataclass(frozen=True,slots=True)
class MocBoundaryResult:
    point:MocFlowPoint
    iterations:int
    equation_residual_max:float
    rounding_adjusted_residual_max:float
    coordinate_rounding_allowance_max:float
    representation_limited:bool
    source_term_per_m:float|None
    source_term_status:str
    solution_kind:str
    boundary_kind:str
    inputs:MocAxisInputs|MocWallInputs
    numerical_policy:str=_POLICY
    physical_model_qualified:bool=False
    full_nozzle_solution:bool=False
    discretization_error_bound:None=None

    def __post_init__(self):
        kind=text(self.boundary_kind,'boundary_kind');object.__setattr__(self,'boundary_kind',kind)
        if kind not in ('axis','wall'):raise DomainError('Unknown characteristic boundary.',field='boundary_kind')
        source=_owned(self.inputs,MocAxisInputs if kind=='axis' else MocWallInputs,'Boundary inputs')
        object.__setattr__(self,'inputs',source)
        point=_owned(self.point,MocFlowPoint,'Boundary point');object.__setattr__(self,'point',point)
        (_axis_topology if kind=='axis' else _wall_topology)(source,point)
        object.__setattr__(self,'iterations',integer(self.iterations,'iterations'))
        for name in ('equation_residual_max','rounding_adjusted_residual_max','coordinate_rounding_allowance_max'):
            object.__setattr__(self,name,real(getattr(self,name),name))
        if self.source_term_per_m is not None:object.__setattr__(self,'source_term_per_m',real(self.source_term_per_m,'source_term_per_m'))
        for name in ('source_term_status','solution_kind','numerical_policy'):
            object.__setattr__(self,name,text(getattr(self,name),name))
        for name in ('representation_limited','physical_model_qualified','full_nozzle_solution'):
            if type(getattr(self,name)) is not bool:raise DomainError('Boundary flags must be actual booleans.',field=name)
        if self.physical_model_qualified or self.full_nozzle_solution or self.discretization_error_bound is not None or self.numerical_policy!=_POLICY:
            raise DomainError('Boundary process cannot claim full-nozzle or physical qualification.',field='Boundary result')
        expected=(_solve_axis if kind=='axis' else _solve_wall)(source)
        if any(getattr(self,name)!=value for name,value in expected.items()):
            raise DomainError('Boundary result does not replay from owned inputs.',field='Boundary result')


def solve_moc_axis_point(*,incoming,axis_history,gamma=1.4,axisymmetric=True,max_iterations=128):
    """Reflect an incoming C- onto r=0 using two prior axis velocities.

    The quadratic regularity approximation is not a sonic starting condition
    or a bound on discretization error. History points must already exist.
    """
    inputs=MocAxisInputs(incoming,axis_history,gamma,axisymmetric,max_iterations)
    return MocBoundaryResult(**_solve_axis(inputs),boundary_kind='axis',inputs=inputs)


def solve_moc_wall_point(*,incoming,wall_segment,gamma=1.4,axisymmetric=True,max_iterations=128):
    """Intersect an incoming C+ with a prescribed finite upper-wall segment.

    The accepted x lies in (x0,x1]; no boundary extrapolation is accepted.
    This imposes slip on given geometry, not a thrust-optimal wall design.
    """
    inputs=MocWallInputs(incoming,wall_segment,gamma,axisymmetric,max_iterations)
    return MocBoundaryResult(**_solve_wall(inputs),boundary_kind='wall',inputs=inputs)
