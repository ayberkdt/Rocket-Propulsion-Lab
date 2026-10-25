"""Interior characteristic intersection for a constant-gamma, irrotational gas.

This module provides one unit process, not a complete nozzle designer. Axis
and wall boundary processes are in :mod:`.moc_boundary`; sonic/throat
starting lines, shock handling and wall design are outside its scope. Neither
this process nor the mapped minimum-length contour of :mod:`.moc` is
validated against hardware.

Cartesian compatibility follows Aldo and Argrow, equations (1a,b) [1]_. The
cylindrical continuity term has ``sigma = 1`` (planar ``sigma = 0``); the
sentence assigning ``sigma = 0`` to both cases in that scanned report is a
typographical inconsistency, and the cylindrical divergence form of
continuity [3]_ determines the convention used here. Finite-difference
coefficient averaging follows NASA TR R-33 [2]_.

References
----------
.. [1] NASA-MOC-CARTESIAN-1993. Aldo and Argrow, Cartesian characteristic and
       compatibility equations (1a,b). https://ntrs.nasa.gov/citations/19940019183
.. [2] NASA-TRR33-AXISYMMETRIC. Guentert and Neumann, 1959, radius source and
       averaged characteristic coefficients.
       https://ntrs.nasa.gov/citations/19980232191
.. [3] NASA-TND4960-AXIS-CONTINUITY. Anderson, 1969, Appendix C equation C1 and
       streamline definitions. https://ntrs.nasa.gov/citations/19690005234
"""
from collections.abc import Mapping
from dataclasses import dataclass, fields
from fractions import Fraction
from functools import lru_cache
from math import asin, atan2, cos, degrees, inf, isfinite, nextafter, pi, sin, tan, ulp

from rocket_propulsion.core.errors import ConvergenceError, DomainError
from rocket_propulsion.core.validation import integer, real, text
from rocket_propulsion.propulsion.burns.numeric_outputs import _finite_rational
from rocket_propulsion.propulsion.burns.probability_numerics import _sqrt_output

_POLICY = 'averaged_cartesian_interior_characteristics_exact_binary_ledgers_v1'
_TOLERANCE = Fraction(64 * ulp(1.0))


def _owned(value, kind, label):
    """Return a fresh ``kind`` instance copied from a mapping or a ``kind`` record.

    Raises
    ------
    DomainError
        If ``value`` is neither a mapping nor a ``kind`` instance, or if
        construction from it fails with a type, value, key, index, attribute
        or overflow error.
    """
    try:
        if isinstance(value, Mapping):return kind(**dict(value))
        if isinstance(value, kind):return kind(**{f.name:getattr(value,f.name) for f in fields(kind)})
    except (TypeError, ValueError, KeyError, IndexError, AttributeError, OverflowError) as error:
        raise DomainError('Malformed characteristic record.',field=label) from error
    raise DomainError('Typed characteristic record required.',field=label)


@dataclass(frozen=True, slots=True)
class MocFlowPoint:
    """A supersonic flow state at one point; positions in m, angle in degrees.

    Attributes
    ----------
    x_m, radius_m : float
        Axial and radial (or transverse, planar) coordinates in m.
    mach : float
        Mach number, greater than one.
    flow_angle_deg : float
        Flow angle from the axis in degrees, ``|theta| < 90``.

    Raises
    ------
    DomainError
        If the state is not supersonic, not axially directed, or its axial
        Mach component ``M cos(theta)`` does not exceed one (the unit
        processes march in ``x``).

    Notes
    -----
    The characteristic equations used by the unit processes [1]_ are written
    with ``x`` as the marching coordinate, which requires a supersonic axial
    velocity component.

    References
    ----------
    .. [1] NASA-MOC-CARTESIAN-1993. Aldo and Argrow, Cartesian characteristic and
           compatibility equations (1a,b). https://ntrs.nasa.gov/citations/19940019183
    """
    x_m:float
    radius_m:float
    mach:float
    flow_angle_deg:float

    def __post_init__(self):
        for f in fields(self):object.__setattr__(self,f.name,real(getattr(self,f.name),f.name))
        if self.mach<=1 or abs(self.flow_angle_deg)>=90:
            raise DomainError('An axially directed supersonic state is required.',field='Characteristic point')
        theta=_finite_rational(Fraction(self.flow_angle_deg)*Fraction(pi)/180,'Flow angle radians')
        if Fraction(self.mach)*Fraction(cos(theta))<=1:
            raise DomainError('This x-parameterized unit process requires supersonic axial Mach.',field='Characteristic point')


@dataclass(frozen=True, slots=True)
class MocUnitInputs:
    """Validated inputs of :func:`solve_moc_interior_point`.

    Attributes
    ----------
    plus_parent : MocFlowPoint
        Upstream point on the ``C+`` (left-running) characteristic.
    minus_parent : MocFlowPoint
        Upstream point on the ``C-`` (right-running) characteristic.
    gamma : float, default 1.4
        Ratio of specific heats, greater than one.
    axisymmetric : bool, default True
        Include the cylindrical source term.
    max_iterations : int, default 128
        Corrector budget, between 1 and 256.

    Raises
    ------
    DomainError
        If a field is malformed, the parents coincide, or an axisymmetric
        parent lies on or below the axis.
    """

    plus_parent:MocFlowPoint
    minus_parent:MocFlowPoint
    gamma:float=1.4
    axisymmetric:bool=True
    max_iterations:int=128

    def __post_init__(self):
        for name in ('plus_parent','minus_parent'):
            object.__setattr__(self,name,_owned(getattr(self,name),MocFlowPoint,name))
        g=real(self.gamma,'gamma')
        if g<=1:raise DomainError('Gamma must exceed one.',field='gamma')
        object.__setattr__(self,'gamma',g)
        if type(self.axisymmetric) is not bool:raise DomainError('Axisymmetric selector must be a boolean.',field='axisymmetric')
        count=integer(self.max_iterations,'max_iterations')
        if not 1<=count<=256:raise DomainError('Iteration budget must be between 1 and 256.',field='max_iterations')
        object.__setattr__(self,'max_iterations',count)
        a,b=self.plus_parent,self.minus_parent
        if (a.x_m,a.radius_m)==(b.x_m,b.radius_m):raise DomainError('Distinct parents are required.',field='Characteristic parents')
        if self.axisymmetric and min(a.radius_m,b.radius_m)<=0:
            raise DomainError('Axisymmetric interior parents require positive radii; the axis needs a separate boundary process.',field='radius_m')


def _coefficients(point, gamma, family, axisymmetric):
    """Return characteristic coefficients of ``point`` for one family.

    Velocities are normalized by the critical speed of sound ``a*`` and
    ``(a/a*)^2 = (gamma + 1)/(2 + (gamma - 1) M^2)`` is kept exact.

    Parameters
    ----------
    point : MocFlowPoint
        Flow state.
    gamma : float
        Ratio of specific heats.
    family : int
        ``+1`` for ``C+`` or ``-1`` for ``C-``.
    axisymmetric : bool
        Include the cylindrical source term.

    Returns
    -------
    tuple
        ``(lambda, Q, R, S, u, v)`` with characteristic slope
        ``lambda = tan(theta + family mu)``, ``Q = u^2 - a^2``,
        ``R = 2 u v - Q lambda`` and ``S = sigma a^2 v / r``, such that the
        compatibility relation along the characteristic is
        ``Q du + R dv - S dx = 0`` [1]_, [2]_.

    References
    ----------
    .. [1] NASA-MOC-CARTESIAN-1993. Aldo and Argrow, Cartesian characteristic and
           compatibility equations (1a,b). https://ntrs.nasa.gov/citations/19940019183
    .. [2] NASA-TRR33-AXISYMMETRIC. Guentert and Neumann, 1959, radius source and
           averaged characteristic coefficients.
           https://ntrs.nasa.gov/citations/19980232191
    """
    m=Fraction(point.mach);g=Fraction(gamma)
    a2=(g+1)/(2+(g-1)*m*m)
    speed=_sqrt_output(m*m*a2,'Characteristic speed / a*')
    theta=_finite_rational(Fraction(point.flow_angle_deg)*Fraction(pi)/180,'Flow angle radians')
    u=Fraction(_finite_rational(Fraction(speed)*Fraction(cos(theta)),'Axial velocity / a*'))
    v=Fraction(_finite_rational(Fraction(speed)*Fraction(sin(theta)),'Radial velocity / a*'))
    mu=asin(_finite_rational(1/m,'Inverse Mach'))
    slope=Fraction(tan(theta+family*mu))
    aa=u*u-a2
    bb=2*u*v-aa*slope
    source=a2*v/Fraction(point.radius_m) if axisymmetric else Fraction(0)
    return slope,aa,bb,source,u,v


def _intersection(a,b,lp,lm):
    """Return the exact intersection of the lines through ``a`` and ``b``.

    The slopes ``lp`` and ``lm`` are used as given (binary64 values); no
    parallelism tolerance is applied.

    Raises
    ------
    ConvergenceError
        If the slopes are exactly equal.
    """
    denominator=lp-lm
    if denominator==0:raise ConvergenceError('Characteristic slopes coincide.',code='moc_unit_parallel')
    x=(Fraction(b.radius_m)-Fraction(a.radius_m)+lp*Fraction(a.x_m)-lm*Fraction(b.x_m))/denominator
    y=Fraction(a.radius_m)+lp*(x-Fraction(a.x_m))
    return x,y


def _midpoint_equations(inputs, point):
    """Return parent and candidate coefficients with family-wise averages.

    The averaged ``(lambda, Q, R, S)`` of each family is the mean of its
    parent and candidate-point values, the averaged-coefficient corrector of
    [1]_.

    References
    ----------
    .. [1] NASA-TRR33-AXISYMMETRIC. Guentert and Neumann, 1959, radius source and
           averaged characteristic coefficients.
           https://ntrs.nasa.gov/citations/19980232191
    """
    plus=_coefficients(inputs.plus_parent,inputs.gamma,1,inputs.axisymmetric)
    minus=_coefficients(inputs.minus_parent,inputs.gamma,-1,inputs.axisymmetric)
    pplus=_coefficients(point,inputs.gamma,1,inputs.axisymmetric)
    pminus=_coefficients(point,inputs.gamma,-1,inputs.axisymmetric)
    averaged=tuple(tuple((a+b)/2 for a,b in zip(parent[:4],candidate[:4]))
                   for parent,candidate in ((plus,pplus),(minus,pminus)))
    return plus,minus,pplus,averaged


def _coordinate_cell_radius(value):
    """Return the larger half-gap from ``value`` to its binary64 neighbours.

    It describes coordinate representation only, not a physical error.
    """
    return max(abs(Fraction(n)-Fraction(value))/2 for n in (nextafter(value,-inf),nextafter(value,inf)) if isfinite(n))


def _residual(inputs,point):
    """Return the residuals of the discrete interior equations at ``point``.

    Returns
    -------
    tuple of Fraction
        ``(raw, adjusted, allowance)``: the largest normalized geometric
        (characteristic line) or compatibility residual; the same with each
        geometric residual reduced by its coordinate rounding allowance; and
        the largest allowance. They measure the represented discrete algebra,
        not the physical or truncation error.
    """
    plus,minus,at_point,averaged=_midpoint_equations(inputs,point)
    a,b=inputs.plus_parent,inputs.minus_parent
    scale=max(abs(Fraction(a.x_m)-Fraction(b.x_m)),abs(Fraction(a.radius_m)-Fraction(b.radius_m)))
    velocity=max(abs(v) for row in (plus,minus) for v in row[4:])
    values=[];adjusted=[];allowances=[]
    for parent,coefficients,avg in ((a,plus,averaged[0]),(b,minus,averaged[1])):
        slope,aa,bb,source=avg
        dx=Fraction(point.x_m)-Fraction(parent.x_m)
        dy=Fraction(point.radius_m)-Fraction(parent.radius_m)
        geometry=abs(dy-slope*dx)/scale
        allowance=(_coordinate_cell_radius(point.radius_m)+abs(slope)*_coordinate_cell_radius(point.x_m))/scale
        values.append(geometry);adjusted.append(max(Fraction(0),geometry-allowance));allowances.append(allowance)
        normalizer=max(abs(aa),abs(bb))*velocity
        if normalizer==0:raise ConvergenceError('Degenerate characteristic equation.',code='moc_unit_singular')
        compatibility=abs(aa*(at_point[4]-coefficients[4])+bb*(at_point[5]-coefficients[5])-source*dx)/normalizer
        values.append(compatibility);adjusted.append(compatibility)
    return max(values),max(adjusted),max(allowances)


def _upper_residual(value):
    """Return the smallest binary64 value not below the exact rational ``value``."""
    number=float(value)
    return nextafter(number,inf) if Fraction(number)<value else number


def _point_from_velocity(x,y,u,v,gamma):
    """Return the flow point at ``(x, y)`` for velocity ``(u, v)`` in units of ``a*``.

    ``(a/a*)^2 = ((gamma + 1) - (gamma - 1)(u^2 + v^2))/2`` from the energy
    equation; Mach number is ``|V|/a``.

    Raises
    ------
    ConvergenceError
        If the state is not on the supersonic, forward-flowing energy branch.
    """
    g=Fraction(gamma);speed2=u*u+v*v
    a2=((g+1)-(g-1)*speed2)/2
    if a2<=0 or speed2<=a2 or u<=0:
        raise ConvergenceError('Interior iteration left the supersonic energy branch.',code='moc_unit_invalid_iteration')
    mach=_sqrt_output(speed2/a2,'Characteristic Mach')
    angle=degrees(atan2(_finite_rational(v,'Radial velocity'),_finite_rational(u,'Axial velocity')))
    return MocFlowPoint(_finite_rational(x,'Characteristic x'),_finite_rational(y,'Characteristic radius'),mach,angle)


def _check_output_topology(point,inputs):
    """Return ``point`` if it lies strictly downstream of both parents.

    Raises
    ------
    ConvergenceError
        If binary64 output collapses the downstream ordering.
    """
    if point.x_m<=max(inputs.plus_parent.x_m,inputs.minus_parent.x_m):
        raise ConvergenceError('The exported point cannot resolve downstream topology.',code='moc_unit_unresolved_topology')
    return point


@lru_cache(maxsize=16)
def _solve(inputs):
    """Solve the interior point; return the fields of :class:`MocUnitResult`.

    A predictor intersects the parent characteristics; the corrector repeats
    the intersection and the 2x2 compatibility solve with averaged
    coefficients [1]_, [2]_ until the rounding-adjusted residual is within 64
    units in the last place of one. Constant planar flow and parallel axial
    axisymmetric flow are returned directly.

    Raises
    ------
    ConvergenceError
        If the intersection is not downstream, reaches the axis, becomes
        singular, leaves the supersonic branch, or the budget is exhausted.

    References
    ----------
    .. [1] NASA-MOC-CARTESIAN-1993. Aldo and Argrow, Cartesian characteristic and
           compatibility equations (1a,b). https://ntrs.nasa.gov/citations/19940019183
    .. [2] NASA-TRR33-AXISYMMETRIC. Guentert and Neumann, 1959, radius source and
           averaged characteristic coefficients.
           https://ntrs.nasa.gov/citations/19980232191
    """
    a,b=inputs.plus_parent,inputs.minus_parent
    left=_coefficients(a,inputs.gamma,1,inputs.axisymmetric)
    right=_coefficients(b,inputs.gamma,-1,inputs.axisymmetric)
    x,y=_intersection(a,b,left[0],right[0])
    if x<=max(Fraction(a.x_m),Fraction(b.x_m)):
        raise ConvergenceError('The parent families do not intersect downstream.',code='moc_unit_noncausal')
    # Constant planar flow, or parallel axial axisymmetric flow, is an exact
    # solution. Preserve Mach directly instead of inverting a rounded limiting
    # velocity, which cannot distinguish very high Mach states.
    constant=(a.mach==b.mach and a.flow_angle_deg==b.flow_angle_deg
              and (not inputs.axisymmetric or a.flow_angle_deg==0))
    if constant:
        point=MocFlowPoint(_finite_rational(x,'Characteristic x'),_finite_rational(y,'Characteristic radius'),a.mach,a.flow_angle_deg)
        _check_output_topology(point,inputs)
        if inputs.axisymmetric and point.radius_m<=0:
            raise ConvergenceError('Interior intersection reached the axis; a boundary solve is required.',code='moc_unit_axis_boundary')
        residual,adjusted,allowance=_residual(inputs,point)
        if adjusted>_TOLERANCE:
            raise ConvergenceError('Constant-state intersection could not be represented.',code='moc_unit_unresolved')
        return point,1,_upper_residual(residual),'constant_state_intersection',_upper_residual(adjusted),_upper_residual(allowance),residual>_TOLERANCE
    point=MocFlowPoint(_finite_rational(x,'Characteristic x'),_finite_rational(y,'Characteristic radius'),
                       _finite_rational((Fraction(a.mach)+Fraction(b.mach))/2,'Mean Mach'),(a.flow_angle_deg+b.flow_angle_deg)/2)
    for iteration in range(1,inputs.max_iterations+1):
        if inputs.axisymmetric and point.radius_m<=0:
            raise ConvergenceError('Interior intersection reached the axis; a boundary solve is required.',code='moc_unit_axis_boundary')
        plus,minus,_,averaged=_midpoint_equations(inputs,point)
        lp,ap,bp,sp=averaged[0];lm,am,bm,sm=averaged[1]
        x,y=_intersection(a,b,lp,lm)
        if x<=max(Fraction(a.x_m),Fraction(b.x_m)):
            raise ConvergenceError('Interior iteration left the downstream branch.',code='moc_unit_noncausal')
        cp=ap*plus[4]+bp*plus[5]+sp*(x-Fraction(a.x_m))
        cm=am*minus[4]+bm*minus[5]+sm*(x-Fraction(b.x_m))
        determinant=ap*bm-am*bp
        if determinant==0:raise ConvergenceError('Interior compatibility equations are singular.',code='moc_unit_singular')
        u=(cp*bm-cm*bp)/determinant;v=(ap*cm-am*cp)/determinant
        point=_point_from_velocity(x,y,u,v,inputs.gamma)
        _check_output_topology(point,inputs)
        if inputs.axisymmetric and point.radius_m<=0:
            raise ConvergenceError('Interior intersection reached the axis; a boundary solve is required.',code='moc_unit_axis_boundary')
        residual,adjusted,allowance=_residual(inputs,point)
        if adjusted<=_TOLERANCE:return point,iteration,_upper_residual(residual),'iterated_midpoint_intersection',_upper_residual(adjusted),_upper_residual(allowance),residual>_TOLERANCE
    raise ConvergenceError('Interior characteristic iteration budget exhausted.',code='moc_unit_not_converged')


@dataclass(frozen=True, slots=True)
class MocUnitResult:
    """Interior characteristic point and its discrete-equation diagnostics.

    Construction re-solves the point from ``inputs`` and rejects any field
    that does not match.

    Attributes
    ----------
    point : MocFlowPoint
        Solved point.
    iterations : int
        Corrector iterations used.
    equation_residual_max : float
        Upper bound of the largest normalized discrete residual.
    inputs : MocUnitInputs
        Owned inputs.
    solution_kind : str
        ``"constant_state_intersection"`` or
        ``"iterated_midpoint_intersection"``.
    rounding_adjusted_residual_max : float
        Residual after subtracting coordinate rounding allowances.
    coordinate_rounding_allowance_max : float
        Largest coordinate rounding allowance.
    representation_limited : bool
        ``True`` when the raw residual exceeds the acceptance tolerance only
        because of coordinate rounding.
    numerical_policy : str
        Identifier of the numerical policy.
    physical_model_qualified, full_nozzle_solution : bool
        Always ``False``.
    discretization_error_bound : None
        Always ``None``; no truncation-error bound is available.

    Raises
    ------
    DomainError
        If a field is malformed, a qualification claim is made, or the record
        does not reproduce from ``inputs``.

    Notes
    -----
    The residuals refer to the averaged characteristic and compatibility
    equations of [1]_ and [2]_.

    References
    ----------
    .. [1] NASA-MOC-CARTESIAN-1993. Aldo and Argrow, Cartesian characteristic and
           compatibility equations (1a,b). https://ntrs.nasa.gov/citations/19940019183
    .. [2] NASA-TRR33-AXISYMMETRIC. Guentert and Neumann, 1959, radius source and
           averaged characteristic coefficients.
           https://ntrs.nasa.gov/citations/19980232191
    """

    point:MocFlowPoint
    iterations:int
    equation_residual_max:float
    inputs:MocUnitInputs
    solution_kind:str
    rounding_adjusted_residual_max:float
    coordinate_rounding_allowance_max:float
    representation_limited:bool
    numerical_policy:str=_POLICY
    physical_model_qualified:bool=False
    full_nozzle_solution:bool=False
    discretization_error_bound:None=None

    def __post_init__(self):
        source=_owned(self.inputs,MocUnitInputs,'Unit process inputs');object.__setattr__(self,'inputs',source)
        point=_owned(self.point,MocFlowPoint,'Unit process point');object.__setattr__(self,'point',point)
        _check_output_topology(point,source)
        count=integer(self.iterations,'iterations');object.__setattr__(self,'iterations',count)
        value=real(self.equation_residual_max,'equation_residual_max');object.__setattr__(self,'equation_residual_max',value)
        for name in ('rounding_adjusted_residual_max','coordinate_rounding_allowance_max'):
            object.__setattr__(self,name,real(getattr(self,name),name))
        if type(self.representation_limited) is not bool:raise DomainError('Representation flag must be a boolean.',field='representation_limited')
        policy=text(self.numerical_policy,'numerical_policy');object.__setattr__(self,'numerical_policy',policy)
        kind=text(self.solution_kind,'solution_kind');object.__setattr__(self,'solution_kind',kind)
        for name in ('physical_model_qualified','full_nozzle_solution'):
            if type(getattr(self,name)) is not bool or getattr(self,name):raise DomainError('Interior unit process cannot claim full-flow or physical qualification.',field=name)
        if self.discretization_error_bound is not None or policy!=_POLICY or (
                point,count,value,kind,self.rounding_adjusted_residual_max,self.coordinate_rounding_allowance_max,self.representation_limited)!=_solve(source):
            raise DomainError('Unit process result does not replay from its inputs.',field='MOC unit process')


def solve_moc_interior_point(*,plus_parent,minus_parent,gamma=1.4,axisymmetric=True,max_iterations=128):
    """Intersect C+ and C- with averaged Cartesian compatibility coefficients.

    Parameters
    ----------
    plus_parent, minus_parent : MocFlowPoint or mapping
        Upstream points on the ``C+`` and ``C-`` characteristics.
    gamma : float, default 1.4
        Ratio of specific heats, greater than one.
    axisymmetric : bool, default True
        Include the cylindrical source term ``a^2 v / r``.
    max_iterations : int, default 128
        Corrector budget, between 1 and 256.

    Returns
    -------
    MocUnitResult
        Solved point and its discrete-equation diagnostics.

    Raises
    ------
    DomainError
        If the inputs are invalid.
    ConvergenceError
        If no downstream supersonic intersection is found within the budget,
        or the point reaches the axis.

    Notes
    -----
    Along each characteristic ``dr/dx = tan(theta +- mu)`` and
    ``Q du + R dv - S dx = 0`` [1]_, with coefficients averaged between the
    parent and the new point [2]_. The residual bounds only the represented
    discrete algebra, not the characteristic integration/truncation error.
    Axisymmetric radii must stay positive: no zero-radius regularity limit or
    sonic boundary is assumed. Coordinate rounding cells are reported
    separately and subtracted only from geometric residuals for the iteration
    decision; they are not measurement uncertainty or a continuous
    characteristic error estimate.

    References
    ----------
    .. [1] NASA-MOC-CARTESIAN-1993. Aldo and Argrow, Cartesian characteristic and
           compatibility equations (1a,b). https://ntrs.nasa.gov/citations/19940019183
    .. [2] NASA-TRR33-AXISYMMETRIC. Guentert and Neumann, 1959, radius source and
           averaged characteristic coefficients.
           https://ntrs.nasa.gov/citations/19980232191
    """
    inputs=MocUnitInputs(plus_parent,minus_parent,gamma,axisymmetric,max_iterations)
    point,count,residual,kind,adjusted,allowance,limited=_solve(inputs)
    return MocUnitResult(point,count,residual,inputs,kind,adjusted,allowance,limited)
