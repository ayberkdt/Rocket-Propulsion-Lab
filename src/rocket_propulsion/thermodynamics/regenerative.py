"""Prescribed constant-property bulk coolant balance for regenerative channels.

The model computes the bulk coolant temperature rise for a prescribed
absorbed heat load and the Darcy major pressure loss in identical parallel
channels with equal flow. It does not predict wall temperature, conjugate
heat transfer, boiling, cavitation, manifold distribution or hardware
performance, and results carry ``physical_model_qualified=False``.

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
"""
from __future__ import annotations

from collections.abc import Mapping
from collections.abc import Set as AbstractSet
from dataclasses import dataclass, fields
from functools import lru_cache
from itertools import islice

from rocket_propulsion.core.errors import DomainError
from rocket_propulsion.core.metadata import WarningMessage
from rocket_propulsion.core.validation import _ordered_items, integer, real, text

from .heat_transfer import _owned
from .regenerative_numerics import _COOLANT_POLICY, _values

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

    Raises
    ------
    DomainError
        If a field is malformed, a flag is not an actual ``bool``, or any
        value differs from the replay of ``inputs``.

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


def calculate_regenerative_cooling(*,absorbed_heat_w,coolant_mass_flow_kg_s,coolant_specific_heat_j_kg_k,
        coolant_inlet_temperature_k,coolant_density_kg_m3,coolant_dynamic_viscosity_pa_s,channel_count,
        channel_flow_area_m2,hydraulic_diameter_m,channel_length_m,roughness_m,maximum_coolant_temperature_k,
        laminar_darcy_poiseuille_number=64.):
    """Evaluate a prescribed bulk coolant energy balance and channel pressure loss.

    Parameters
    ----------
    absorbed_heat_w : float
        Net heat absorbed by the coolant, W; must be positive.
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

    Returns
    -------
    RegenerativeCoolingResult
        Temperatures, hydraulics, warnings and qualification flags. This is
        a prescribed bulk balance, not a coupled cooling design.

    Raises
    ------
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
    source=RegenerativeInputs(absorbed_heat_w,coolant_mass_flow_kg_s,coolant_specific_heat_j_kg_k,
        coolant_inlet_temperature_k,coolant_density_kg_m3,coolant_dynamic_viscosity_pa_s,channel_count,
        channel_flow_area_m2,hydraulic_diameter_m,channel_length_m,roughness_m,maximum_coolant_temperature_k,
        laminar_darcy_poiseuille_number)
    return RegenerativeCoolingResult(**_solve(source),inputs=source)
