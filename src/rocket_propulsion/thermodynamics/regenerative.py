"""Prescribed constant-property bulk coolant balance for regenerative channels.

The model computes the bulk coolant temperature rise for a prescribed
absorbed heat load and the Darcy major pressure loss in identical parallel
channels with equal flow. It does not predict conjugate heat transfer,
boiling, cavitation, manifold distribution or hardware performance, and
results carry ``physical_model_qualified=False``.

Two optional extensions are available through keyword-only arguments. A
catalog coolant [4]_ can fill the constant property inputs; the record's
reference condition is then reported as a warning on
:class:`CoolantSelection`. A one-dimensional steady wall-conduction step
[5]_ estimates the hot-gas-side wall temperature from a prescribed local
heat flux and checks it against the wall material's sourced limit; see
:mod:`rocket_propulsion.thermodynamics.wall_conduction`. Without these
arguments the result is unchanged.

The laminar friction factor uses a supplied Darcy ``f Re`` coefficient
(cross-section dependent, [1]_); transition and turbulent friction use the
Swamee-Jain explicit formula [2]_ in the base-ten convention also used by
EPANET [3]_. The friction correlation's domain for cooling channels is not
qualified (``friction_correlation_domain_qualified=False``).

References
----------
.. [1] NASA-GFSSP6-NONCIRCULAR. NASA GFSSP v6 user manual, printed p. 40
       (PDF p. 80), equation 28,
       https://www.nasa.gov/wp-content/uploads/2024/04/gfssp-v6-usermanual.pdf
.. [2] SWAMEE-JAIN-1976-IDENTITY. Swamee and Jain, "Explicit Equations
       for Pipe-Flow Problems", J. Hydraulics Division 102(5),
       657-664, 1976, https://doi.org/10.1061/JYCEAJ.0004542
.. [3] EPANET22-FRICTION-SOURCE. EPANET v2.2, ``src/hydcoeffs.c``,
       ``frictionFactor`` and base-ten logarithm conversion,
       https://raw.githubusercontent.com/OpenWaterAnalytics/EPANET/v2.2/src/hydcoeffs.c
.. [4] NIST-WEBBOOK-FLUID. E. W. Lemmon, M. O. McLinden and D. G. Friend,
       "Thermophysical Properties of Fluid Systems", NIST Chemistry WebBook,
       SRD 69, https://webbook.nist.gov/chemistry/fluid/
.. [5] NASA-SP125-REGEN-WALL. D. K. Huzel and D. H. Huang, "Design of
       Liquid Propellant Rocket Engines", NASA SP-125, 2nd ed., 1971,
       printed p. 104 (PDF p. 113), equations (4-19) to (4-22),
       https://ntrs.nasa.gov/api/citations/19710019929/downloads/19710019929.pdf
"""
from __future__ import annotations

from collections.abc import Mapping
from collections.abc import Set as AbstractSet
from dataclasses import dataclass, fields
from functools import lru_cache
from itertools import islice

from rocket_propulsion.core.errors import DomainError, InputError
from rocket_propulsion.core.metadata import WarningMessage
from rocket_propulsion.core.validation import _ordered_items, integer, real, strings, text

from .heat_transfer import _owned
from .materials import get_coolant
from .regenerative_numerics import _COOLANT_POLICY, _values
from .wall_conduction import WallConductionResult, calculate_wall_conduction

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
    """Owned, validated inputs of the bulk coolant balance.

    Attributes
    ----------
    absorbed_heat_w : float
        Prescribed net heat absorbed by the coolant, W; must be positive.
    coolant_mass_flow_kg_s : float
        Total coolant mass flow, kg/s.
    coolant_specific_heat_j_kg_k : float
        Constant coolant specific heat, J/(kg K).
    coolant_inlet_temperature_k : float
        Bulk inlet temperature, K.
    coolant_density_kg_m3 : float
        Constant coolant density, kg/m^3.
    coolant_dynamic_viscosity_pa_s : float
        Constant coolant dynamic viscosity, Pa s.
    channel_count : int
        Number of identical parallel channels; at least one.
    channel_flow_area_m2 : float
        Flow area of one channel, m^2.
    hydraulic_diameter_m : float
        Prescribed hydraulic diameter, m.
    channel_length_m : float
        Channel length, m.
    roughness_m : float
        Absolute wall roughness, m; may be zero.
    maximum_coolant_temperature_k : float
        User temperature limit, K; must exceed the inlet temperature.
    laminar_darcy_poiseuille_number : float, optional
        Laminar Darcy ``f Re`` coefficient ``K``. The default 64 assumes a
        circular equivalent; for other cross-sections ``K`` equals four times
        the cross-section-dependent Poiseuille number [1]_. It is not a
        Fanning factor.

    Raises
    ------
    DomainError
        If a quantity is not a finite real (or, for ``channel_count``, an
        actual integer), is non-positive (zero roughness excepted), or the
        temperature limit does not exceed the inlet temperature.

    Notes
    -----
    The temperature limit is a user requirement, not an inferred saturation
    or material limit.

    References
    ----------
    .. [1] NASA-GFSSP6-NONCIRCULAR. NASA GFSSP v6 user manual, printed p. 40
           (PDF p. 80), equation 28,
           https://www.nasa.gov/wp-content/uploads/2024/04/gfssp-v6-usermanual.pdf
    """
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
    """Normalize coolant warning records; at most six are accepted."""
    if isinstance(value,(str,bytes,bytearray,memoryview,Mapping,AbstractSet)):
        raise DomainError('Warnings must be ordered records.',field='Coolant warnings')
    try:items=tuple(islice(value,7))
    except (TypeError,ValueError,KeyError,IndexError,OverflowError) as error:
        raise DomainError('Malformed coolant warning iterable.',field='Coolant warnings') from error
    if len(items)>6:raise DomainError('Too many coolant warnings.',field='Coolant warnings')
    result=[]
    for warning in items:
        warning=_owned(warning,WarningMessage,'Coolant warning')
        severity=text(warning.severity,'Warning severity')
        if severity not in ('info','warning','critical'):raise DomainError('Unknown warning severity.',field='Warning severity')
        result.append(WarningMessage(text(warning.code,'Warning code'),text(warning.message,'Warning message'),severity,
            None if warning.field is None else text(warning.field,'Warning field')))
    return tuple(result)


_CATALOG_FIELDS=(('coolant_specific_heat_j_kg_k','specific_heat'),('coolant_density_kg_m3','density'),
 ('coolant_dynamic_viscosity_pa_s','dynamic_viscosity'))


def _selection_warnings(record,filled,inlet,outlet):
    """Return the reference-condition and critical-temperature warnings of a catalog fill.

    References
    ----------
    .. [1] NIST-WEBBOOK-FLUID. E. W. Lemmon, M. O. McLinden and D. G. Friend,
           "Thermophysical Properties of Fluid Systems", NIST Chemistry
           WebBook, SRD 69, https://webbook.nist.gov/chemistry/fluid/
    """
    warnings=[]
    if filled:
        parts=[]
        for name,prop in _CATALOG_FIELDS:
            if name in filled:
                item=getattr(record,prop)
                where='unstated temperature' if item.temperature_k is None else f'{item.temperature_k:.6g} K'
                if item.pressure_pa is not None:
                    where+=f', {item.pressure_pa:.6g} Pa'
                parts.append(f'{prop} {item.value:.6g} {item.unit} at {where}')
        warnings.append(WarningMessage('coolant_catalog_reference_condition',
            f'Catalog {record.key} properties are constants at their reference conditions ({"; ".join(parts)}); '
            f'they are not evaluated at the bulk inlet ({inlet:.6g} K) or outlet ({outlet:.6g} K) temperature.',
            'warning','coolant'))
    else:
        warnings.append(WarningMessage('coolant_catalog_properties_not_used',
            'All coolant property inputs were supplied by the caller; the catalog record is used only for its critical-temperature check.',
            'info','coolant'))
    critical=record.critical_temperature
    if critical is not None and outlet>critical.value:
        warnings.append(WarningMessage('coolant_outlet_above_critical_temperature',
            f'Bulk outlet temperature exceeds the {record.key} critical temperature ({critical.value:.6g} K); '
            'the constant single-phase liquid property assumption does not hold.','critical','coolant'))
    return tuple(warnings)


@dataclass(frozen=True,slots=True)
class CoolantSelection:
    """Record of a catalog coolant used to fill constant coolant properties.

    Attributes
    ----------
    coolant : str
        Catalog key of :func:`~rocket_propulsion.thermodynamics.materials.get_coolant`;
        normalized on construction.
    filled_fields : tuple of str
        Input names filled from the catalog: a subset of
        ``coolant_specific_heat_j_kg_k``, ``coolant_density_kg_m3`` and
        ``coolant_dynamic_viscosity_pa_s``. Inputs given explicitly by the
        caller are not listed.
    inlet_temperature_k, outlet_temperature_k : float
        Bulk temperatures of the balance, K, used for the warnings.
    warnings : tuple of WarningMessage
        Reference-condition caveat (or an info record when nothing was
        filled) and, when the outlet exceeds the catalog critical
        temperature, a critical warning.

    Raises
    ------
    InputError
        If the coolant key is unknown.
    DomainError
        If a filled field is unknown, repeated or absent from the catalog
        record, a temperature is not positive, or the warnings do not replay.

    References
    ----------
    .. [1] NIST-WEBBOOK-FLUID. E. W. Lemmon, M. O. McLinden and D. G. Friend,
           "Thermophysical Properties of Fluid Systems", NIST Chemistry
           WebBook, SRD 69, https://webbook.nist.gov/chemistry/fluid/
    .. [2] NISTIR-6646-RP1. J. W. Magee et al., "Thermophysical Properties
           Measurements and Models for Rocket Propellant RP-1: Phase I",
           NISTIR 6646, 2007,
           https://nvlpubs.nist.gov/nistpubs/Legacy/IR/nistir6646.pdf
    """
    coolant:str
    filled_fields:tuple[str,...]
    inlet_temperature_k:float
    outlet_temperature_k:float
    warnings:tuple[WarningMessage,...]=()

    def __post_init__(self):
        record=get_coolant(self.coolant)
        object.__setattr__(self,'coolant',record.key)
        filled=strings(self.filled_fields,'Filled coolant fields')
        known=dict(_CATALOG_FIELDS)
        if len(set(filled))!=len(filled) or any(name not in known or getattr(record,known[name]) is None for name in filled):
            raise DomainError('Filled coolant fields must be distinct catalog-backed inputs.',field='filled_fields')
        object.__setattr__(self,'filled_fields',tuple(name for name,_ in _CATALOG_FIELDS if name in filled))
        for name in ('inlet_temperature_k','outlet_temperature_k'):
            value=real(getattr(self,name),name)
            if value<=0:
                raise DomainError('Coolant selection temperatures must be positive.',field=name)
            object.__setattr__(self,name,value)
        object.__setattr__(self,'warnings',_warnings(self.warnings))
        if self.warnings!=_selection_warnings(record,self.filled_fields,self.inlet_temperature_k,self.outlet_temperature_k):
            raise DomainError('Coolant selection warnings do not replay.',field='Coolant selection')

    def filled_values(self):
        """Return the catalog values of the filled inputs.

        Returns
        -------
        dict of str to float
            Input name mapped to the SI catalog value.

        References
        ----------
        .. [1] NIST-WEBBOOK-FLUID. E. W. Lemmon, M. O. McLinden and D. G.
               Friend, "Thermophysical Properties of Fluid Systems", NIST
               Chemistry WebBook, SRD 69, https://webbook.nist.gov/chemistry/fluid/
        """
        record=get_coolant(self.coolant)
        known=dict(_CATALOG_FIELDS)
        return {name:getattr(record,known[name]).value for name in self.filled_fields}


@dataclass(frozen=True,slots=True)
class RegenerativeCoolingResult:
    """Replayed bulk coolant balance with explicit qualification flags.

    Attributes
    ----------
    inlet_temperature_k, outlet_temperature_k : float
        Bulk coolant temperatures, K.
    temperature_rise_k : float
        ``H / (m cp)``, K.
    absorbed_heat_w : float
        Prescribed absorbed heat, W.
    coolant_capacity_w_k : float
        ``m cp``, W/K.
    channel_velocity_m_s : float
        ``m / (rho N A)``, m/s.
    reynolds_number : float
        ``m D / (N A mu)``.
    darcy_friction_factor : float
        Darcy friction factor of the selected regime.
    pressure_drop_pa : float
        Darcy major loss ``f (L/D) rho V^2 / 2``, Pa.
    energy_balance_residual_w : float
        Always zero; labelled by ``residual_kind`` as a ledger identity, not
        a measured heat closure.
    experimental : bool
        Always True for this model.
    warnings : tuple of WarningMessage
        Model, ledger, limit, regime and rounding warnings.
    assumptions : tuple of str
        Declared modelling assumptions.
    inputs : RegenerativeInputs
        Owned inputs from which the result is replayed.
    flow_regime : {"laminar", "transition", "turbulent"}
        Regime selected from the exact Reynolds number.
    friction_model : str
        Identifier of the friction branch used.
    physical_model_qualified, friction_correlation_domain_qualified : bool
        Always False.
    temperature_limit_exceeded : bool
        Exact comparison of the outlet temperature with the user limit.
    temperature_limit_margin_k : float or None
        Limit minus exact outlet temperature; ``None`` when that nonzero
        value is not representable as a finite float.
    temperature_limit_margin_unrepresentable : bool
        True when the margin is ``None`` for that reason.
    displayed_temperature_energy_residual_w : float or None
        Absolute discrepancy obtained by inserting the rounded outlet
        temperature into the exact capacity ledger; it diagnoses output
        rounding only.
    displayed_temperature_energy_residual_unrepresentable : bool
        True when that discrepancy is nonzero but not representable.
    residual_kind, numerical_policy, model_version, model_fidelity : str
        Identifiers of the residual meaning, numerical policy and model.
    coolant_selection : CoolantSelection or None
        Catalog coolant used to fill property inputs; ``None`` when no
        catalog coolant was selected.
    wall_conduction : WallConductionResult or None
        Optional wall-conduction step; ``None`` when not requested.

    Raises
    ------
    DomainError
        If a field is malformed, a flag is not an actual ``bool``, any
        value differs from the replay of ``inputs``, a filled input differs
        from its catalog value, or the extension records disagree with the
        bulk balance.

    References
    ----------
    .. [1] NASA-GFSSP6-NONCIRCULAR. NASA GFSSP v6 user manual, printed p. 40
           (PDF p. 80), equation 28,
           https://www.nasa.gov/wp-content/uploads/2024/04/gfssp-v6-usermanual.pdf
    .. [2] SWAMEE-JAIN-1976-IDENTITY. Swamee and Jain, "Explicit Equations
           for Pipe-Flow Problems", J. Hydraulics Division 102(5),
           657-664, 1976, https://doi.org/10.1061/JYCEAJ.0004542
    .. [3] EPANET22-FRICTION-SOURCE. EPANET v2.2, ``src/hydcoeffs.c``,
           ``frictionFactor`` and base-ten logarithm conversion,
           https://raw.githubusercontent.com/OpenWaterAnalytics/EPANET/v2.2/src/hydcoeffs.c
    """
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
    coolant_selection:CoolantSelection|None=None
    wall_conduction:WallConductionResult|None=None

    @property
    def all_warnings(self):
        """Return the bulk-balance warnings followed by those of the extensions.

        Returns
        -------
        tuple of WarningMessage
            ``warnings``, then ``coolant_selection.warnings`` and
            ``wall_conduction.warnings`` when present.

        References
        ----------
        .. [1] NASA-SP125-REGEN-WALL. D. K. Huzel and D. H. Huang, "Design of
               Liquid Propellant Rocket Engines", NASA SP-125, 2nd ed., 1971,
               printed p. 104 (PDF p. 113), equations (4-19) to (4-22),
               https://ntrs.nasa.gov/api/citations/19710019929/downloads/19710019929.pdf
        """
        extra=()
        if self.coolant_selection is not None:
            extra+=self.coolant_selection.warnings
        if self.wall_conduction is not None:
            extra+=self.wall_conduction.warnings
        return self.warnings+extra

    def __post_init__(self):
        source=_owned(self.inputs,RegenerativeInputs,'Coolant inputs');object.__setattr__(self,'inputs',source)
        self._own_extensions()
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
        selection=self.coolant_selection
        if selection is not None:
            if (selection.inlet_temperature_k!=source.coolant_inlet_temperature_k
                    or selection.outlet_temperature_k!=self.outlet_temperature_k
                    or any(getattr(source,name)!=value for name,value in selection.filled_values().items())):
                raise DomainError('Coolant selection does not match the bulk balance inputs.',field='coolant_selection')
        wall=self.wall_conduction
        if wall is not None and wall.inputs.coolant_bulk_temperature_k not in (None,self.outlet_temperature_k):
            raise DomainError('Wall step must use the bulk outlet temperature.',field='wall_conduction')

    def _own_extensions(self):
        """Replace extension records (or their mappings) by owned, replayed copies."""
        if self.coolant_selection is not None:
            object.__setattr__(self,'coolant_selection',_owned(self.coolant_selection,CoolantSelection,'Coolant selection'))
        if self.wall_conduction is not None:
            object.__setattr__(self,'wall_conduction',_owned(self.wall_conduction,WallConductionResult,'Wall conduction'))


@lru_cache(maxsize=64)
def _solve(inputs):
    """Return every result field of the bulk coolant balance for ``inputs``.

    References
    ----------
    .. [1] NASA-GFSSP6-NONCIRCULAR. NASA GFSSP v6 user manual, printed p. 40
           (PDF p. 80), equation 28,
           https://www.nasa.gov/wp-content/uploads/2024/04/gfssp-v6-usermanual.pdf
    .. [2] SWAMEE-JAIN-1976-IDENTITY. Swamee and Jain, "Explicit Equations
           for Pipe-Flow Problems", J. Hydraulics Division 102(5),
           657-664, 1976, https://doi.org/10.1061/JYCEAJ.0004542
    .. [3] EPANET22-FRICTION-SOURCE. EPANET v2.2, ``src/hydcoeffs.c``,
           ``frictionFactor`` and base-ten logarithm conversion,
           https://raw.githubusercontent.com/OpenWaterAnalytics/EPANET/v2.2/src/hydcoeffs.c
    """
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


def calculate_regenerative_cooling(*,absorbed_heat_w,coolant_mass_flow_kg_s,coolant_specific_heat_j_kg_k=None,
        coolant_inlet_temperature_k,coolant_density_kg_m3=None,coolant_dynamic_viscosity_pa_s=None,channel_count,
        channel_flow_area_m2,hydraulic_diameter_m,channel_length_m,roughness_m,maximum_coolant_temperature_k,
        laminar_darcy_poiseuille_number=64.,coolant=None,wall_material=None,wall_thickness_m=None,
        wall_heat_flux_w_m2=None,coolant_side_wall_temperature_k=None,
        coolant_side_heat_transfer_coefficient_w_m2_k=None,wall_thermal_conductivity_w_m_k=None):
    """Evaluate a prescribed bulk coolant energy balance and channel pressure loss.

    Parameters
    ----------
    absorbed_heat_w : float
        Net heat absorbed by the coolant, W; must be positive.
    coolant_mass_flow_kg_s : float
        Total coolant mass flow, kg/s.
    coolant_specific_heat_j_kg_k : float, optional
        Constant coolant specific heat, J/(kg K). Required unless filled
        from ``coolant``.
    coolant_inlet_temperature_k : float
        Bulk inlet temperature, K.
    coolant_density_kg_m3 : float, optional
        Constant coolant density, kg/m^3. Required unless filled from
        ``coolant``.
    coolant_dynamic_viscosity_pa_s : float, optional
        Constant coolant dynamic viscosity, Pa s. Required unless filled
        from ``coolant``.
    channel_count : int
        Number of identical parallel channels.
    channel_flow_area_m2 : float
        Flow area of one channel, m^2.
    hydraulic_diameter_m : float
        Hydraulic diameter, m.
    channel_length_m : float
        Channel length, m.
    roughness_m : float
        Absolute wall roughness, m; may be zero.
    maximum_coolant_temperature_k : float
        User temperature limit, K.
    laminar_darcy_poiseuille_number : float, optional
        Laminar Darcy ``f Re`` coefficient; 64 for a circular equivalent.
    coolant : str, optional
        Key of :func:`~rocket_propulsion.thermodynamics.materials.get_coolant`.
        Each of ``coolant_specific_heat_j_kg_k``, ``coolant_density_kg_m3``
        and ``coolant_dynamic_viscosity_pa_s`` left as ``None`` is filled
        with the catalog value at its reference condition [4]_; explicit
        values are kept. The fill is recorded in ``coolant_selection``.
    wall_material : str, optional
        Wall-material catalog key; requests the wall-conduction step.
    wall_thickness_m : float, optional
        Hot-gas to coolant wall thickness, m.
    wall_heat_flux_w_m2 : float, optional
        Local heat flux through the wall, W/m^2 (for example a peak Bartz
        station value; it is independent of ``absorbed_heat_w``).
    coolant_side_wall_temperature_k : float, optional
        Prescribed coolant-side wall temperature, K.
    coolant_side_heat_transfer_coefficient_w_m2_k : float, optional
        Coolant film coefficient, W/(m^2 K); the coolant-side wall
        temperature is then ``T_out + q / h_c`` with the bulk outlet
        temperature ``T_out`` [5]_. Give exactly one of this and
        ``coolant_side_wall_temperature_k``.
    wall_thermal_conductivity_w_m_k : float, optional
        Conductivity replacing the catalog value, W/(m K).

    Returns
    -------
    RegenerativeCoolingResult
        Temperatures, hydraulics, warnings and qualification flags. This is
        a prescribed bulk balance, not a coupled cooling design. Extension
        warnings are on ``coolant_selection`` and ``wall_conduction`` and are
        collected by ``all_warnings``.

    Raises
    ------
    InputError
        If the coolant or wall material is unknown, a property to be filled
        is absent from the catalog (code ``"coolant_property_unavailable"``),
        the wall step is incomplete (code ``"incomplete_wall_conduction"``)
        or its coolant-side specification is ambiguous.
    DomainError
        If an input fails the checks of :class:`RegenerativeInputs`, an
        exported value is outside the nonzero finite range, or the
        Swamee-Jain logarithm argument is exactly one (code
        ``"singular_swamee_jain_argument"``).
    ConvergenceError
        If a rounding or branch decision cannot be resolved within the
        numerical precision budget.

    Notes
    -----
    With ``N`` channels of area ``A``:

    * ``rise = H / (m cp)`` and ``outlet = inlet + rise``;
    * ``V = m / (rho N A)`` and ``Re = m D / (N A mu)``;
    * laminar (``Re < 2300``): ``f = K / Re`` with the supplied ``K`` [1]_;
    * ``2300 <= Re < 4000`` is labelled transition and keeps the
      Swamee-Jain formula as an unqualified extension; ``Re >= 4000`` uses
      it as the turbulent branch: ``a = e/(3.7 D) + 5.74 / Re^(9/10)`` and
      ``f = 0.25 / log10(a)^2`` [2]_, [3]_;
    * ``dp = f (L/D) rho V^2 / 2`` (major loss only; no entrance, bend,
      manifold or acceleration losses).

    Branches use the exact Reynolds number, even when its exported value
    rounds to a threshold. For ``a > 1`` the unsquared inverse-friction root
    is negative; the squared extension is retained with a critical warning.
    The temperature-limit decision compares the exact outlet expression,
    because displayed temperatures can hide an exceedance.

    The optional wall step uses ``T_wg = T_wc + q t / k`` [5]_ and flags
    ``wall_conduction.wall_limit_exceeded`` against the material's sourced
    maximum service temperature, or its melting point when no service
    temperature is catalogued.

    References
    ----------
    .. [1] NASA-GFSSP6-NONCIRCULAR. NASA GFSSP v6 user manual, printed p. 40
           (PDF p. 80), equation 28,
           https://www.nasa.gov/wp-content/uploads/2024/04/gfssp-v6-usermanual.pdf
    .. [2] SWAMEE-JAIN-1976-IDENTITY. Swamee and Jain, "Explicit Equations
           for Pipe-Flow Problems", J. Hydraulics Division 102(5),
           657-664, 1976, https://doi.org/10.1061/JYCEAJ.0004542
    .. [3] EPANET22-FRICTION-SOURCE. EPANET v2.2, ``src/hydcoeffs.c``,
           ``frictionFactor`` and base-ten logarithm conversion,
           https://raw.githubusercontent.com/OpenWaterAnalytics/EPANET/v2.2/src/hydcoeffs.c
    .. [4] NIST-WEBBOOK-FLUID. E. W. Lemmon, M. O. McLinden and D. G. Friend,
           "Thermophysical Properties of Fluid Systems", NIST Chemistry
           WebBook, SRD 69, https://webbook.nist.gov/chemistry/fluid/
    .. [5] NASA-SP125-REGEN-WALL. D. K. Huzel and D. H. Huang, "Design of
           Liquid Propellant Rocket Engines", NASA SP-125, 2nd ed., 1971,
           printed p. 104 (PDF p. 113), equations (4-19) to (4-22),
           https://ntrs.nasa.gov/api/citations/19710019929/downloads/19710019929.pdf
    """
    properties={'coolant_specific_heat_j_kg_k':coolant_specific_heat_j_kg_k,'coolant_density_kg_m3':coolant_density_kg_m3,
        'coolant_dynamic_viscosity_pa_s':coolant_dynamic_viscosity_pa_s}
    record=None
    filled=[]
    if coolant is not None:
        record=get_coolant(coolant)
        for name,prop in _CATALOG_FIELDS:
            if properties[name] is not None:
                continue
            item=getattr(record,prop)
            if item is None:
                raise InputError(f'Catalog coolant {record.key} has no verified {prop}; supply {name}.',
                    code='coolant_property_unavailable',field=name)
            properties[name]=item.value
            filled.append(name)
    wall_arguments=(wall_material,wall_thickness_m,wall_heat_flux_w_m2,coolant_side_wall_temperature_k,
        coolant_side_heat_transfer_coefficient_w_m2_k,wall_thermal_conductivity_w_m_k)
    wall_requested=any(value is not None for value in wall_arguments)
    if wall_requested and None in (wall_material,wall_thickness_m,wall_heat_flux_w_m2):
        raise InputError('The wall step needs wall_material, wall_thickness_m and wall_heat_flux_w_m2.',
            code='incomplete_wall_conduction',field='wall_material')
    source=RegenerativeInputs(absorbed_heat_w,coolant_mass_flow_kg_s,properties['coolant_specific_heat_j_kg_k'],
        coolant_inlet_temperature_k,properties['coolant_density_kg_m3'],properties['coolant_dynamic_viscosity_pa_s'],
        channel_count,channel_flow_area_m2,hydraulic_diameter_m,channel_length_m,roughness_m,
        maximum_coolant_temperature_k,laminar_darcy_poiseuille_number)
    values=_solve(source)
    outlet=values['outlet_temperature_k']
    selection=None
    if record is not None:
        selection=CoolantSelection(record.key,tuple(filled),source.coolant_inlet_temperature_k,outlet,
            _selection_warnings(record,tuple(filled),source.coolant_inlet_temperature_k,outlet))
    wall=None
    if wall_requested:
        wall=calculate_wall_conduction(wall_material=wall_material,wall_thickness_m=wall_thickness_m,
            heat_flux_w_m2=wall_heat_flux_w_m2,coolant_side_wall_temperature_k=coolant_side_wall_temperature_k,
            coolant_side_heat_transfer_coefficient_w_m2_k=coolant_side_heat_transfer_coefficient_w_m2_k,
            coolant_bulk_temperature_k=None if coolant_side_heat_transfer_coefficient_w_m2_k is None else outlet,
            wall_thermal_conductivity_w_m_k=wall_thermal_conductivity_w_m_k)
    return RegenerativeCoolingResult(**values,inputs=source,coolant_selection=selection,wall_conduction=wall)
