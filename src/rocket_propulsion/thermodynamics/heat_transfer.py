"""Owned prescribed-state Bartz profiles; no coupled wall or hardware validation."""
from __future__ import annotations
from collections.abc import Mapping
from dataclasses import dataclass, fields
from fractions import Fraction
from functools import lru_cache
from math import inf, isfinite, nextafter

from rocket_propulsion.core.errors import DomainError
from rocket_propulsion.core.metadata import WarningMessage
from rocket_propulsion.core.validation import real, text, _ordered_items
from rocket_propulsion.geometry import NozzleContour, NozzleStation
from .bartz_numerics import _profile_exports, _BARTZ_POLICY

_INTEGRATION='piecewise_linear_heat_flux_and_radius_on_straight_meridional_segments'
_RESIDUAL='integral_ledger_identity_not_measured_energy_closure'
_CORRELATION='Bartz prescribed constant-property gas-side correlation'
_PROPERTIES='user-supplied constant stagnation bulk gas properties'
_SOURCE='NASA SP-125 equations 4-13/4-14; Linne et al. 2012 recovery and signed flux equations 1/4'
_URL='https://ntrs.nasa.gov/citations/19710019929'
_VALIDITY=('prescribed gas state and uniform wall temperature; no conjugate wall/coolant solve',
           'empirical operating-domain applicability is unqualified',
           'no radiation, reacting boundary-layer, film cooling or deposition model',
           'linear sampled-profile reconstruction; no continuous-profile quadrature error bound')


def _owned(value,kind,label):
    if isinstance(value,Mapping):
        try:return kind(**dict(value))
        except (TypeError,ValueError,KeyError,IndexError,OverflowError) as error:
            raise DomainError('Malformed retained mapping.',field=label) from error
    if isinstance(value,kind):
        return kind(**{f.name:getattr(value,f.name) for f in fields(kind)})
    raise DomainError('Typed retained record required.',field=label)


@dataclass(frozen=True,slots=True)
class BartzContourStation:
    x_m:float
    radius_m:float
    area_ratio:float
    mach:float

    def __post_init__(self):
        for f in fields(self):object.__setattr__(self,f.name,real(getattr(self,f.name),f.name))
        if self.radius_m<=0 or self.area_ratio<=0 or self.mach<0:
            raise DomainError('Radius/area ratio must be positive; Mach nonnegative.',field='Thermal contour')


def _rounding_cell(value):
    """Exact binary64 rounding-cell endpoints; these are not measurement errors."""
    center=Fraction(value);before=Fraction(nextafter(value,-inf));after=nextafter(value,inf)
    return (center+before)/2, (center+Fraction(after))/2 if isfinite(after) else center+(center-before)/2


@dataclass(frozen=True,slots=True)
class BartzContour:
    throat_radius_m:float
    stations:tuple[BartzContourStation,...]

    def __post_init__(self):
        radius=real(self.throat_radius_m,'Throat radius')
        if radius<=0:raise DomainError('Throat radius must be positive.',field='Throat radius')
        object.__setattr__(self,'throat_radius_m',radius)
        points=tuple(_owned(s,BartzContourStation,'Contour station') for s in _ordered_items(self.stations,'Contour stations'))
        if len(points)<2:raise DomainError('A surface profile needs at least two stations.',field='Contour stations')
        rlo,rhi=_rounding_cell(radius)
        for index,point in enumerate(points):
            if index and point.x_m<=points[index-1].x_m:
                raise DomainError('Contour x must increase strictly.',field='Contour stations')
            lo,hi=_rounding_cell(point.radius_m);alo,ahi=_rounding_cell(point.area_ratio)
            if hi<rlo or ahi<(lo/rhi)**2 or alo>(hi/rlo)**2:
                raise DomainError('Area ratio and radius disagree beyond binary64 representation cells.',field='Contour geometry')
        object.__setattr__(self,'stations',points)


def _contour(value):
    if isinstance(value,NozzleContour):
        points=[]
        for s in _ordered_items(value.stations,'Contour stations'):
            if not isinstance(s,NozzleStation):raise DomainError('Nozzle station required.',field='Contour stations')
            points.append(BartzContourStation(s.x_m,s.radius_m,s.area_ratio,s.mach))
        return BartzContour(value.throat_radius_m,tuple(points))
    return _owned(value,BartzContour,'Thermal contour')


@dataclass(frozen=True,slots=True)
class BartzInputs:
    contour:BartzContour
    chamber_pressure_pa:float
    chamber_temperature_k:float
    wall_temperature_k:float
    characteristic_velocity_m_s:float
    dynamic_viscosity_pa_s:float
    specific_heat_j_kg_k:float
    prandtl_number:float
    gamma:float
    throat_radius_of_curvature_m:float

    def __post_init__(self):
        object.__setattr__(self,'contour',_contour(self.contour))
        for f in fields(self):
            if f.name=='contour':continue
            value=real(getattr(self,f.name),f.name)
            if value<=0:raise DomainError('Thermal quantities must be positive.',field=f.name)
            object.__setattr__(self,f.name,value)
        if self.gamma<=1:raise DomainError('Gamma must exceed one.',field='gamma')
        if self.wall_temperature_k>=self.chamber_temperature_k:
            raise DomainError('This prescribed model requires wall temperature below chamber temperature.',field='wall_temperature_k')


@dataclass(frozen=True,slots=True)
class HeatFluxStation:
    """Finite station data; equation replay belongs to the enclosing profile."""
    x_m:float
    radius_m:float
    mach:float
    adiabatic_wall_temperature_k:float
    gas_side_coefficient_w_m2_k:float
    heat_flux_w_m2:float
    cumulative_heat_w:float

    def __post_init__(self):
        for f in fields(self):object.__setattr__(self,f.name,real(getattr(self,f.name),f.name))
        if self.radius_m<=0 or self.mach<0 or self.adiabatic_wall_temperature_k<=0 or self.gas_side_coefficient_w_m2_k<=0:
            raise DomainError('Invalid finite thermal station.',field='Heat-flux station')


def _warnings(value):
    # At most two base, two legacy-screening and one signed-flux warning.
    from itertools import islice
    from collections.abc import Set
    if isinstance(value,(str,bytes,bytearray,memoryview,Mapping,Set)):
        raise DomainError('Warnings must be ordered records.',field='Thermal warnings')
    try:items=tuple(islice(value,6))
    except (TypeError,ValueError,KeyError,IndexError,OverflowError) as error:
        raise DomainError('Malformed thermal warnings.',field='Thermal warnings') from error
    if len(items)>5:raise DomainError('Too many thermal warnings.',field='Thermal warnings')
    result=[]
    for item in items:
        item=_owned(item,WarningMessage,'Thermal warning')
        severity=text(item.severity,'Warning severity')
        if severity not in ('info','warning','critical'):raise DomainError('Unknown warning severity.',field='Warning severity')
        result.append(WarningMessage(text(item.code,'Warning code'),text(item.message,'Warning message'),
            severity,None if item.field is None else text(item.field,'Warning field')))
    return tuple(result)


@dataclass(frozen=True,slots=True)
class BartzHeatTransferResult:
    stations:tuple[HeatFluxStation,...]
    peak_heat_flux_w_m2:float
    peak_location_m:float
    total_heat_load_w:float
    energy_balance_residual_w:float
    correlation:str
    property_model:str
    validity:tuple[str,...]
    warnings:tuple[WarningMessage,...]
    source:str
    source_url:str
    inputs:BartzInputs|None=None
    contains_negative_heat_flux:bool=False
    physical_model_qualified:bool=False
    quadrature_error_bound_w:float|None=None
    integration_model:str=_INTEGRATION
    residual_kind:str=_RESIDUAL
    numerical_policy:str=_BARTZ_POLICY
    pi_basis:str='fixed_binary64_math_pi'

    def __post_init__(self):
        source=_owned(self.inputs,BartzInputs,'Bartz inputs');object.__setattr__(self,'inputs',source)
        for name in ('peak_heat_flux_w_m2','peak_location_m','total_heat_load_w','energy_balance_residual_w'):
            object.__setattr__(self,name,real(getattr(self,name),name))
        for name in ('contains_negative_heat_flux','physical_model_qualified'):
            if type(getattr(self,name)) is not bool:raise DomainError('Thermal flags must be actual booleans.',field=name)
        if self.quadrature_error_bound_w is not None:
            object.__setattr__(self,'quadrature_error_bound_w',real(self.quadrature_error_bound_w,'Quadrature error bound'))
        for name in ('correlation','property_model','source','source_url','integration_model','residual_kind','numerical_policy','pi_basis'):
            object.__setattr__(self,name,text(getattr(self,name),name))
        object.__setattr__(self,'validity',tuple(text(s,'Thermal scope') for s in _ordered_items(self.validity,'Thermal scope',size=len(_VALIDITY))))
        object.__setattr__(self,'warnings',_warnings(self.warnings))
        stations=tuple(_owned(s,HeatFluxStation,'Heat-flux station') for s in _ordered_items(self.stations,'Heat-flux stations',size=len(source.contour.stations)))
        object.__setattr__(self,'stations',stations)
        expected=_solve(source)
        if any(getattr(self,name)!=value for name,value in expected.items() if name!='station_values') or any(
                any(getattr(s,name)!=value for name,value in wanted.items()) for s,wanted in zip(stations,expected['station_values'])):
            raise DomainError('Thermal result does not replay from owned inputs.',field='Bartz result')


@lru_cache(maxsize=64)
def _solve(inputs):
    values=_profile_exports(inputs)
    warnings=[WarningMessage('bartz_prescribed_model_unqualified',
        'Prescribed constant gas properties and wall temperature; physical applicability, wall/coolant closure and hardware are unqualified.'),
        WarningMessage('thermal_integral_ledger_identity','Zero residual is an integral ledger identity, not measured energy closure.','info')]
    diameter=2*Fraction(inputs.contour.throat_radius_m)
    if not Fraction(.01)<=diameter<=1:
        warnings.append(WarningMessage('bartz_diameter_extrapolation','Outside legacy diameter screening thresholds; these are not an independently validated correlation domain.'))
    if inputs.chamber_pressure_pa>30e6:
        warnings.append(WarningMessage('bartz_pressure_extrapolation','Above the legacy pressure screening threshold; it does not establish physical applicability.'))
    if values['contains_negative_heat_flux']:
        warnings.append(WarningMessage('negative_prescribed_heat_flux','Some prescribed wall locations transfer heat toward the gas; signed flux and integrated load are retained.'))
    values.update(energy_balance_residual_w=0.,correlation=_CORRELATION,property_model=_PROPERTIES,
        validity=_VALIDITY,warnings=tuple(warnings),source=_SOURCE,source_url=_URL,
        physical_model_qualified=False,quadrature_error_bound_w=None,integration_model=_INTEGRATION,
        residual_kind=_RESIDUAL,numerical_policy=_BARTZ_POLICY,pi_basis='fixed_binary64_math_pi')
    return values


def calculate_bartz_heat_transfer(*,contour,chamber_pressure_pa,chamber_temperature_k,wall_temperature_k,
        characteristic_velocity_m_s,dynamic_viscosity_pa_s,specific_heat_j_kg_k,prandtl_number,gamma,throat_radius_of_curvature_m):
    """Evaluate a replayed prescribed profile; signed heat is positive into the wall."""
    inputs=BartzInputs(contour,chamber_pressure_pa,chamber_temperature_k,wall_temperature_k,
        characteristic_velocity_m_s,dynamic_viscosity_pa_s,specific_heat_j_kg_k,prandtl_number,gamma,throat_radius_of_curvature_m)
    values=dict(_solve(inputs));stations=tuple(HeatFluxStation(**s) for s in values.pop('station_values'))
    return BartzHeatTransferResult(**values,stations=stations,inputs=inputs)
