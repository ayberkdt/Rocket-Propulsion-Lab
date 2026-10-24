"""Owned constant-property coolant balances with explicit unqualified friction basis."""
from __future__ import annotations
from collections.abc import Mapping, Set
from dataclasses import dataclass, fields
from functools import lru_cache
from itertools import islice

from rocket_propulsion.core.errors import DomainError
from rocket_propulsion.core.metadata import WarningMessage
from rocket_propulsion.core.validation import real, integer, text, _ordered_items
from .heat_transfer import _owned
from .regenerative_numerics import _values, _COOLANT_POLICY

_RESIDUAL='bulk_energy_ledger_identity_not_measured_energy_closure'
_ASSUMPTIONS=('positive absorbed heat; prescribed single-phase constant coolant properties',
 'equal flow in identical parallel channels; supplied area and hydraulic diameter',
 'laminar Darcy f*Re is a supplied coefficient; default 64 assumes a circular equivalent',
 'legacy Re<2300 laminar switch; [2300,4000) uses an explicitly unqualified Swamee-Jain extension',
 'Darcy major loss only; no entrance, bend, manifold or acceleration losses',
 'no boiling, wall conduction, conjugate cooling or physical qualification',
 'zero energy residual is an exact ledger identity; displayed-temperature discrepancy is separate')
_FLOAT_FIELDS=('inlet_temperature_k','outlet_temperature_k','temperature_rise_k','absorbed_heat_w',
 'coolant_capacity_w_k','channel_velocity_m_s','reynolds_number','darcy_friction_factor','pressure_drop_pa','energy_balance_residual_w')
_FLAGS=('experimental','physical_model_qualified','friction_correlation_domain_qualified','temperature_limit_exceeded',
 'temperature_limit_margin_unrepresentable','displayed_temperature_energy_residual_unrepresentable')


@dataclass(frozen=True,slots=True)
class RegenerativeInputs:
    absorbed_heat_w:float
    coolant_mass_flow_kg_s:float
    coolant_specific_heat_j_kg_k:float
    coolant_inlet_temperature_k:float
    coolant_density_kg_m3:float
    coolant_dynamic_viscosity_pa_s:float
    channel_count:int
    channel_flow_area_m2:float
    hydraulic_diameter_m:float
    channel_length_m:float
    roughness_m:float
    maximum_coolant_temperature_k:float
    laminar_darcy_poiseuille_number:float=64.

    def __post_init__(self):
        for field in fields(self):
            name=field.name
            if name=='channel_count':
                value=integer(self.channel_count,name)
                if value<1:raise DomainError('Channel count must be positive.',field=name)
            else:
                value=real(getattr(self,name),name)
                if value<0 or value==0 and name!='roughness_m':
                    raise DomainError('Positive coolant quantity required; roughness may be zero.',field=name)
            object.__setattr__(self,name,value)
        if self.maximum_coolant_temperature_k<=self.coolant_inlet_temperature_k:
            raise DomainError('Maximum coolant temperature must exceed inlet temperature.',field='maximum_coolant_temperature_k')


def _warnings(value):
    if isinstance(value,(str,bytes,bytearray,memoryview,Mapping,Set)):
        raise DomainError('Warnings must be ordered records.',field='Coolant warnings')
    try:items=tuple(islice(value,7))
    except (TypeError,ValueError,KeyError,IndexError,OverflowError) as error:
        raise DomainError('Malformed coolant warning iterable.',field='Coolant warnings') from error
    if len(items)>6:raise DomainError('Too many coolant warnings.',field='Coolant warnings')
    result=[]
    for value in items:
        value=_owned(value,WarningMessage,'Coolant warning')
        severity=text(value.severity,'Warning severity')
        if severity not in ('info','warning','critical'):raise DomainError('Unknown warning severity.',field='Warning severity')
        result.append(WarningMessage(text(value.code,'Warning code'),text(value.message,'Warning message'),severity,
            None if value.field is None else text(value.field,'Warning field')))
    return tuple(result)


@dataclass(frozen=True,slots=True)
class RegenerativeCoolingResult:
    inlet_temperature_k:float
    outlet_temperature_k:float
    temperature_rise_k:float
    absorbed_heat_w:float
    coolant_capacity_w_k:float
    channel_velocity_m_s:float
    reynolds_number:float
    darcy_friction_factor:float
    pressure_drop_pa:float
    energy_balance_residual_w:float
    experimental:bool
    warnings:tuple[WarningMessage,...]
    assumptions:tuple[str,...]
    inputs:RegenerativeInputs|None=None
    flow_regime:str=''
    friction_model:str=''
    physical_model_qualified:bool=False
    friction_correlation_domain_qualified:bool=False
    temperature_limit_exceeded:bool=False
    temperature_limit_margin_k:float|None=None
    temperature_limit_margin_unrepresentable:bool=False
    displayed_temperature_energy_residual_w:float|None=None
    displayed_temperature_energy_residual_unrepresentable:bool=False
    residual_kind:str=_RESIDUAL
    numerical_policy:str=_COOLANT_POLICY
    model_version:str='2.0'
    model_fidelity:str='prescribed_single_phase_bulk_channel'

    def __post_init__(self):
        source=_owned(self.inputs,RegenerativeInputs,'Coolant inputs');object.__setattr__(self,'inputs',source)
        for name in _FLOAT_FIELDS:object.__setattr__(self,name,real(getattr(self,name),name))
        for name in ('temperature_limit_margin_k','displayed_temperature_energy_residual_w'):
            value=getattr(self,name)
            if value is not None:object.__setattr__(self,name,real(value,name))
        for name in _FLAGS:
            if type(getattr(self,name)) is not bool:raise DomainError('Coolant flags must be actual booleans.',field=name)
        for name in ('flow_regime','friction_model','residual_kind','numerical_policy','model_version','model_fidelity'):
            object.__setattr__(self,name,text(getattr(self,name),name))
        object.__setattr__(self,'warnings',_warnings(self.warnings))
        object.__setattr__(self,'assumptions',tuple(text(x,'Coolant assumption') for x in _ordered_items(self.assumptions,'Coolant assumptions',size=len(_ASSUMPTIONS))))
        if any(getattr(self,name)!=value for name,value in _solve(source).items()):
            raise DomainError('Coolant result does not replay from owned inputs.',field='Regenerative result')


@lru_cache(maxsize=64)
def _solve(inputs):
    result=_values(inputs);argument_side=result.pop('swamee_argument_side')
    warnings=[WarningMessage('coolant_model_unqualified',
        'Prescribed single-phase constant properties and equal channel flows. Supplied area, hydraulic diameter and laminar coefficient do not establish shape validity, boiling limits or physical qualification.'),
        WarningMessage('coolant_energy_ledger_identity','Zero energy residual is a ledger identity, not measured heat closure; displayed-temperature discrepancy is separate.','info')]
    if result['temperature_limit_exceeded']:
        warnings.append(WarningMessage('coolant_temperature_limit_exceeded',
            'Exact bulk outlet temperature exceeds the supplied limit; rounded displayed temperatures may coincide.','critical'))
    if result['flow_regime']=='transition':
        warnings.append(WarningMessage('coolant_transition_regime',
            'Legacy [2300,4000) Reynolds band uses an unqualified Swamee-Jain extension; this is not a validated transition correlation.'))
    if argument_side>0:
        warnings.append(WarningMessage('swamee_jain_nonphysical_inverse_root',
            'The logarithm argument exceeds one: its unsquared inverse-friction root is negative. The squared formula extension is retained but is not physically qualified.','critical'))
    if result['outlet_temperature_k']==inputs.coolant_inlet_temperature_k:
        warnings.append(WarningMessage('coolant_temperature_rise_hidden',
            'Positive exact temperature rise is hidden by outlet-temperature display rounding; retain the rise and displayed-energy discrepancy.'))
    regime=result['flow_regime']
    model='prescribed_laminar_poiseuille' if regime=='laminar' else 'unqualified_swamee_jain_transition_extension' if regime=='transition' else 'swamee_jain_turbulent'
    result.update(experimental=True,warnings=tuple(warnings),assumptions=_ASSUMPTIONS,friction_model=model,
        physical_model_qualified=False,friction_correlation_domain_qualified=False,residual_kind=_RESIDUAL,
        numerical_policy=_COOLANT_POLICY,model_version='2.0',model_fidelity='prescribed_single_phase_bulk_channel')
    return result


def calculate_regenerative_cooling(*,absorbed_heat_w,coolant_mass_flow_kg_s,coolant_specific_heat_j_kg_k,
        coolant_inlet_temperature_k,coolant_density_kg_m3,coolant_dynamic_viscosity_pa_s,channel_count,
        channel_flow_area_m2,hydraulic_diameter_m,channel_length_m,roughness_m,maximum_coolant_temperature_k,
        laminar_darcy_poiseuille_number=64.):
    """Return a replayed prescribed bulk balance, not a coupled cooling design."""
    source=RegenerativeInputs(absorbed_heat_w,coolant_mass_flow_kg_s,coolant_specific_heat_j_kg_k,
        coolant_inlet_temperature_k,coolant_density_kg_m3,coolant_dynamic_viscosity_pa_s,channel_count,
        channel_flow_area_m2,hydraulic_diameter_m,channel_length_m,roughness_m,maximum_coolant_temperature_k,
        laminar_darcy_poiseuille_number)
    return RegenerativeCoolingResult(**_solve(source),inputs=source)
