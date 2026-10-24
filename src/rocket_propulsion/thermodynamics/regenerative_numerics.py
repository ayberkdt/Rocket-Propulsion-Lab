"""Exact coolant ledgers and directed Swamee-Jain formula exports.

These implement a prescribed constant-property bulk model. They do not bound
correlation error, establish cross-section validity or predict boiling.
"""
from decimal import Decimal
from fractions import Fraction
from functools import lru_cache

from rocket_propulsion.core.errors import DomainError, ConvergenceError
from rocket_propulsion.propulsion.burns.numeric_outputs import _finite_rational
from rocket_propulsion.propulsion.burns.marginal_numerics import _Intervals, _log_interval, _PRECISIONS
from rocket_propulsion.propulsion.cycle_numerics import _integer_root
from rocket_propulsion.propulsion.performance_numerics import _positive_enclosed_output
from .bartz_numerics import _power

_COOLANT_POLICY='exact_binary_coolant_ledgers_and_enclosed_swamee_jain_exports_v1'


def _ledgers(inputs):
    f=Fraction;m=f(inputs.coolant_mass_flow_kg_s);cp=f(inputs.coolant_specific_heat_j_kg_k)
    area=inputs.channel_count*f(inputs.channel_flow_area_m2)
    density=f(inputs.coolant_density_kg_m3);diameter=f(inputs.hydraulic_diameter_m)
    capacity=m*cp;rise=f(inputs.absorbed_heat_w)/capacity
    outlet=f(inputs.coolant_inlet_temperature_k)+rise
    velocity=m/(density*area)
    reynolds=m*diameter/(area*f(inputs.coolant_dynamic_viscosity_pa_s))
    pressure_scale=f(inputs.channel_length_m)/diameter*density*velocity**2/2
    return dict(coolant_capacity_w_k=capacity,temperature_rise_k=rise,outlet_temperature_k=outlet,
                channel_velocity_m_s=velocity,reynolds_number=reynolds,pressure_scale=pressure_scale)


def _exact_argument(reynolds,relative_roughness):
    n,d=_integer_root(reynolds.numerator,10),_integer_root(reynolds.denominator,10)
    if n is None or d is None:return None
    return relative_roughness/Fraction(37,10)+Fraction(287,50)*Fraction(d,n)**9


def _decimal_log_integer(value):
    n,d=value.numerator,value.denominator
    if n!=1 and d!=1:return None
    item=n if d==1 else d;power=0
    while item>1 and item%10==0:item//=10;power+=1
    return (power if d==1 else -power) if item==1 else None


@lru_cache(maxsize=2048)
def _argument_bounds(reynolds,relative_roughness,precision):
    math=_Intervals(precision)
    exact=_exact_argument(reynolds,relative_roughness)
    if exact is not None:return math.rational(exact)
    return math.add(math.rational(relative_roughness/Fraction(37,10)),
        math.mul(math.rational(Fraction(287,50)),_power(reynolds,Fraction(-9,10),precision)))


def _friction_basis(reynolds,relative_roughness):
    exact=_exact_argument(reynolds,relative_roughness)
    if exact==1:
        raise DomainError('Swamee-Jain logarithm is exactly zero.',code='singular_swamee_jain_argument',field='roughness_m')
    if exact is not None:
        power=_decimal_log_integer(exact)
        return (None if power is None else Fraction(1,4*power**2)),(-1 if exact<1 else 1)
    for precision in _PRECISIONS:
        lo,hi=_argument_bounds(reynolds,relative_roughness,precision)
        if hi<1:return None,-1
        if lo>1:return None,1
    raise ConvergenceError('Swamee-Jain singularity-side comparison unresolved.',field='Darcy friction factor')


@lru_cache(maxsize=1024)
def _friction_bounds(reynolds,relative_roughness,precision):
    for p in _PRECISIONS:
        if p<precision:continue
        math=_Intervals(p);lo,hi=_argument_bounds(reynolds,relative_roughness,p)
        if lo<=0 or lo<=1<=hi:continue
        logarithm=(math.near.next_minus(math.near.ln(lo)),math.near.next_plus(math.near.ln(hi)))
        if logarithm[0]>0:
            squared=(math.down.multiply(logarithm[0],logarithm[0]),math.up.multiply(logarithm[1],logarithm[1]))
        elif logarithm[1]<0:
            squared=(math.down.multiply(logarithm[1],logarithm[1]),math.up.multiply(logarithm[0],logarithm[0]))
        else:continue
        log_ten=_log_interval(Fraction(10),p)
        numerator=math.mul(log_ten,log_ten)
        return math.div(numerator,math.mul(math.rational(Fraction(4)),squared))
    raise ConvergenceError('Swamee-Jain logarithm denominator unresolved.',field='Darcy friction factor')


def _optional_output(value,label):
    try:return _finite_rational(value,label)
    except DomainError:return None


def _values(inputs):
    ledgers=_ledgers(inputs);reynolds=ledgers['reynolds_number'];pressure_scale=ledgers.pop('pressure_scale')
    result={name:_finite_rational(value,name) for name,value in ledgers.items()}
    regime='laminar' if reynolds<2300 else 'transition' if reynolds<4000 else 'turbulent'
    if regime=='laminar':
        exact_friction=Fraction(inputs.laminar_darcy_poiseuille_number)/reynolds
        argument_side=0
    else:
        roughness=Fraction(inputs.roughness_m)/Fraction(inputs.hydraulic_diameter_m)
        exact_friction,argument_side=_friction_basis(reynolds,roughness)
    if exact_friction is not None:
        friction=_finite_rational(exact_friction,'Darcy friction factor')
        pressure=_finite_rational(exact_friction*pressure_scale,'Pressure drop')
    else:
        friction=_positive_enclosed_output(lambda p:_friction_bounds(reynolds,roughness,p),'Darcy friction factor')
        def pressure_bounds(p):
            math=_Intervals(p)
            return math.mul(math.rational(pressure_scale),_friction_bounds(reynolds,roughness,p))
        pressure=_positive_enclosed_output(pressure_bounds,'Pressure drop')
    margin=Fraction(inputs.maximum_coolant_temperature_k)-ledgers['outlet_temperature_k']
    display_residual=abs(Fraction(inputs.absorbed_heat_w)-ledgers['coolant_capacity_w_k']*(
        Fraction(result['outlet_temperature_k'])-Fraction(inputs.coolant_inlet_temperature_k)))
    exported_margin=_optional_output(margin,'Temperature limit margin')
    exported_display_residual=_optional_output(display_residual,'Displayed-temperature energy discrepancy')
    result.update(inlet_temperature_k=inputs.coolant_inlet_temperature_k,absorbed_heat_w=inputs.absorbed_heat_w,
        darcy_friction_factor=friction,pressure_drop_pa=pressure,flow_regime=regime,
        temperature_limit_exceeded=margin<0,temperature_limit_margin_k=exported_margin,
        temperature_limit_margin_unrepresentable=exported_margin is None,
        displayed_temperature_energy_residual_w=exported_display_residual,
        displayed_temperature_energy_residual_unrepresentable=exported_display_residual is None,
        energy_balance_residual_w=0.,swamee_argument_side=argument_side)
    return result
