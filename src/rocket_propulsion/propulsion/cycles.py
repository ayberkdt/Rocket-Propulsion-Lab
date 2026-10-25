"""Preliminary engine-cycle pressure, routed mass-flow and shaft-power balances.

Pressure-fed cycles are sized for the injector feed-pressure requirement only.
Gas-generator and staged-combustion cycles additionally close an ideal pump/turbine
power balance for an explicitly declared pump routing. These are sizing
equations, not an integrated engine model or a hardware qualification.

References
----------
.. [1] NASA-TMX64749-CYCLE. R. H. Schmucker, 1973, equivalent cycle power balances
       and whole pumped mass flows; theory source, not hardware validation.
.. [2] NASA-ENGINESIM-THEORY. NASA Glenn Research Center, EngineSim theory notes:
       ideal gas turbine temperature/pressure and shaft-work relations.
.. [3] SUTTON-BIBLARZ-2017. G. P. Sutton and O. Biblarz, "Rocket Propulsion
       Elements", 9th ed., Wiley, 2017.
"""
from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass
from fractions import Fraction

from rocket_propulsion.core.errors import DomainError
from rocket_propulsion.core.metadata import WarningMessage
from rocket_propulsion.core.validation import coerce_fields, records, strings

from .burns.numeric_outputs import _finite_rational
from .cycle_numerics import _CYCLE_POLICY, component_power_outputs

_PUMP_FIELDS = ('propellant_density_kg_m3','pump_efficiency','turbine_efficiency',
    'turbine_inlet_temperature_k','turbine_specific_heat_j_kg_k','turbine_gamma','turbine_pressure_ratio')


@dataclass(frozen=True, slots=True)
class EngineCycleInputs:
    """Validated inputs for a preliminary engine-cycle balance.

    Attributes
    ----------
    cycle : {'pressure-fed', 'gas-generator', 'staged-combustion'}
        Cycle architecture; case and underscores are normalized on construction.
    chamber_pressure_pa : float
        Chamber pressure Pc, Pa; positive.
    chamber_mass_flow_kg_s : float
        Propellant mass flow delivered to the main chamber, kg/s; positive.
    inlet_pressure_pa : float
        Supply (tank or pump-inlet) pressure, Pa; non-negative.
    injector_pressure_drop_fraction : float
        Injector pressure drop as a fraction of Pc, in (0, 1].
    propellant_density_kg_m3 : float or None
        Single equivalent liquid density of the pumped propellant, kg/m^3.
    pump_efficiency : float or None
        Pump efficiency, in (0, 1].
    turbine_efficiency : float or None
        Turbine isentropic efficiency, in (0, 1].
    turbine_inlet_temperature_k : float or None
        Turbine total inlet temperature, K.
    turbine_specific_heat_j_kg_k : float or None
        Turbine-gas specific heat at constant pressure, J/(kg K).
    turbine_gamma : float or None
        Turbine-gas ratio of specific heats; greater than one.
    turbine_pressure_ratio : float or None
        Turbine total inlet/outlet pressure ratio; greater than one.
    pump_feed_policy : {'shared-pump', 'external-gg-feed'}
        ``'shared-pump'`` routes generator propellant through the main pump;
        ``'external-gg-feed'`` (gas-generator cycle only) supplies it from a separate,
        unmodelled source.

    Raises
    ------
    DomainError
        If the cycle or routing policy is unknown, a quantity is out of domain, or a
        pumped cycle lacks a pump/turbine input. Pump and turbine inputs may be None
        only for the pressure-fed cycle.

    References
    ----------
    .. [1] NASA-TMX64749-CYCLE. R. H. Schmucker, 1973, equivalent cycle power balances
           and whole pumped mass flows; theory source, not hardware validation.
    .. [2] SUTTON-BIBLARZ-2017. G. P. Sutton and O. Biblarz, "Rocket Propulsion
           Elements", 9th ed., Wiley, 2017.
    """
    cycle: str
    chamber_pressure_pa: float
    chamber_mass_flow_kg_s: float
    inlet_pressure_pa: float
    injector_pressure_drop_fraction: float
    propellant_density_kg_m3: float | None = None
    pump_efficiency: float | None = None
    turbine_efficiency: float | None = None
    turbine_inlet_temperature_k: float | None = None
    turbine_specific_heat_j_kg_k: float | None = None
    turbine_gamma: float | None = None
    turbine_pressure_ratio: float | None = None
    pump_feed_policy: str = 'shared-pump'

    def __post_init__(self):
        coerce_fields(self,'Cycle input')
        cycle = self.cycle.strip().lower().replace('_','-')
        if cycle not in ('pressure-fed','gas-generator','staged-combustion'):
            raise DomainError('Unknown cycle.',field='cycle')
        object.__setattr__(self,'cycle',cycle)
        if self.pump_feed_policy not in ('shared-pump','external-gg-feed'):
            raise DomainError('Unknown pump routing policy.',field='pump_feed_policy')
        if cycle != 'gas-generator' and self.pump_feed_policy != 'shared-pump':
            raise DomainError('External generator feed applies only to gas-generator mode.',field='pump_feed_policy')
        if self.chamber_pressure_pa <= 0 or self.chamber_mass_flow_kg_s <= 0:
            raise DomainError('Chamber pressure and flow must be positive.',field='Chamber inputs')
        if self.inlet_pressure_pa < 0 or not 0 < self.injector_pressure_drop_fraction <= 1:
            raise DomainError('Inlet pressure must be nonnegative and injector fraction in (0,1].',field='Feed inputs')
        for name in _PUMP_FIELDS:
            value = getattr(self,name)
            if value is None:
                if cycle != 'pressure-fed':
                    raise DomainError('Pumped cycles require active pump/turbine inputs.',field=name)
                continue
            if value <= 0 or (name.endswith('efficiency') and value > 1):
                raise DomainError('Invalid positive component quantity or efficiency.',field=name)
            if name in ('turbine_gamma','turbine_pressure_ratio') and value <= 1:
                raise DomainError('Turbine gamma and inlet/outlet pressure ratio must exceed one.',field=name)


@dataclass(frozen=True, slots=True)
class EngineCycleResult:
    """Engine-cycle balance, replayable from its retained typed inputs.

    On construction the record evaluates the balance again from ``inputs`` and refuses
    any output that the evaluation does not reproduce. Mapping inputs and warnings are
    converted to their typed records.

    Attributes
    ----------
    cycle : str
        Normalized cycle name.
    chamber_pressure_pa : float
        Chamber pressure, Pa.
    required_feed_pressure_pa : float
        Pc * (1 + injector_pressure_drop_fraction), Pa.
    pump_pressure_rise_pa : float
        max(0, required feed pressure - inlet pressure), Pa; zero when pressure-fed.
    pump_power_w : float
        Hydraulic pump shaft power, W.
    turbine_specific_work_j_kg : float
        Ideal turbine specific work with isentropic efficiency, J/kg.
    turbine_mass_flow_kg_s : float
        Turbine drive-gas mass flow, kg/s.
    turbine_flow_fraction : float
        Turbine flow divided by total propellant flow.
    turbine_power_w : float
        Turbine shaft power, W; equal to pump power by construction.
    power_balance_residual_w : float
        Zero; an algebraic identity of the model, not measured closure.
    chamber_mass_flow_kg_s : float
        Main-chamber propellant flow, kg/s.
    total_propellant_flow_kg_s : float
        Total propellant drawn from the tanks, kg/s.
    delivered_mass_fraction : float
        Chamber flow divided by total flow.
    warnings : tuple of WarningMessage
        Model-scope and advisory warnings.
    assumptions : tuple of str
        Model assumptions for the evaluated cycle.
    inputs : EngineCycleInputs
        Retained typed inputs.
    pump_mass_flow_kg_s : float
        Flow through the main pump, kg/s.
    inlet_pressure_margin_pa : float or None
        Inlet pressure minus required feed pressure, Pa; None if not representable.
    inlet_pressure_margin_representable : bool
        Whether the signed margin is a finite binary64 value.
    pressure_fed_supply_sufficient : bool or None
        For pressure-fed cycles, whether the declared inlet pressure meets the
        requirement; None for pumped cycles.
    energy_balance_scope : str
        Identifier of what the power balance covers.
    numerical_policy : str
        Identifier of the arithmetic and rounding policy.

    Raises
    ------
    DomainError
        If inputs are missing or malformed, a flag has the wrong type, the numerical
        policy is unknown, or an output does not replay from the inputs.

    References
    ----------
    .. [1] NASA-TMX64749-CYCLE. R. H. Schmucker, 1973, equivalent cycle power balances
           and whole pumped mass flows; theory source, not hardware validation.
    .. [2] NASA-ENGINESIM-THEORY. NASA Glenn Research Center, EngineSim theory notes:
           ideal gas turbine temperature/pressure and shaft-work relations.
    """
    cycle: str
    chamber_pressure_pa: float
    required_feed_pressure_pa: float
    pump_pressure_rise_pa: float
    pump_power_w: float
    turbine_specific_work_j_kg: float
    turbine_mass_flow_kg_s: float
    turbine_flow_fraction: float
    turbine_power_w: float
    power_balance_residual_w: float
    chamber_mass_flow_kg_s: float
    total_propellant_flow_kg_s: float
    delivered_mass_fraction: float
    warnings: tuple[WarningMessage, ...]
    assumptions: tuple[str, ...]
    inputs: EngineCycleInputs | None = None
    pump_mass_flow_kg_s: float = 0.
    inlet_pressure_margin_pa: float | None = 0.
    inlet_pressure_margin_representable: bool = True
    pressure_fed_supply_sufficient: bool | None = None
    energy_balance_scope: str = 'unspecified'
    numerical_policy: str = _CYCLE_POLICY

    def __post_init__(self):
        coerce_fields(self,'Cycle result')
        source = self.inputs
        if isinstance(source,Mapping):
            try:source = EngineCycleInputs(**dict(source))
            except (TypeError,ValueError,KeyError,IndexError,OverflowError) as error:
                raise DomainError('Invalid retained cycle inputs.',field='Cycle inputs') from error
        if not isinstance(source,EngineCycleInputs):
            raise DomainError('Cycle result requires retained typed inputs.',field='Cycle inputs')
        if self.numerical_policy != _CYCLE_POLICY:
            raise DomainError('Unknown cycle numerical policy.',field='Numerical policy')
        if self.pressure_fed_supply_sufficient is not None and not isinstance(self.pressure_fed_supply_sufficient,bool):
            raise DomainError('Supply sufficiency must be boolean or null.',field='Supply sufficiency')
        if not isinstance(self.inlet_pressure_margin_representable,bool):
            raise DomainError('Margin representability must be boolean.',field='Inlet pressure margin')
        object.__setattr__(self,'inputs',source)
        object.__setattr__(self,'warnings',_cycle_warnings(self.warnings))
        object.__setattr__(self,'assumptions',strings(self.assumptions,'Cycle assumptions'))
        expected = _cycle_values(source)
        if any(getattr(self,name) != value for name,value in expected.items()):
            raise DomainError('Cycle outputs do not reproduce retained inputs and routing.',field='Cycle outputs')


def _cycle_warnings(value):
    """Return retained cycle warnings as typed ``WarningMessage`` records.

    Parameters
    ----------
    value : iterable
        ``WarningMessage`` records or mappings of their fields.

    Returns
    -------
    tuple of WarningMessage

    Raises
    ------
    DomainError
        If a warning record is malformed.
    """
    owned = []
    for warning in records(value,(WarningMessage,Mapping),'Cycle warnings'):
        if isinstance(warning,Mapping):
            try:warning = WarningMessage(**dict(warning))
            except (TypeError,ValueError,KeyError,IndexError,OverflowError) as error:
                raise DomainError('Invalid retained cycle warning.',field='Cycle warnings') from error
        owned.append(warning)
    return tuple(owned)


def _cycle_values(inputs):
    """Evaluate the cycle balance on exact rationals and return all result fields.

    Parameters
    ----------
    inputs : EngineCycleInputs
        Validated cycle inputs.

    Returns
    -------
    dict
        Field values for :class:`EngineCycleResult` (excluding ``inputs``).

    Raises
    ------
    DomainError
        If the routing cannot close the power balance or an output is not
        representable.
    ConvergenceError
        If an enclosed turbine expression is unresolved within the precision budget.

    Notes
    -----
    Required feed pressure p_req = Pc * (1 + f_inj); pump pressure rise
    Delta_p = max(0, p_req - p_in); pump specific work s = Delta_p / (rho * eta_p)
    [1]_ [2]_. Turbine specific work
    w_t = eta_t * cp * T_in * (1 - PR**(-(gamma - 1) / gamma)) [3]_. Pressure-fed
    cycles report requirement sizing only: feed-line loss, transient tank pressure,
    cavitation and hardware capability are not modelled.

    References
    ----------
    .. [1] SUTTON-BIBLARZ-2017. G. P. Sutton and O. Biblarz, "Rocket Propulsion
           Elements", 9th ed., Wiley, 2017.
    .. [2] NASA-TMX64749-CYCLE. R. H. Schmucker, 1973, equivalent cycle power balances
           and whole pumped mass flows; theory source, not hardware validation.
    .. [3] NASA-ENGINESIM-THEORY. NASA Glenn Research Center, EngineSim theory notes:
           ideal gas turbine temperature/pressure and shaft-work relations.
    """
    pressure, flow, inlet, fraction = map(Fraction,(inputs.chamber_pressure_pa,
        inputs.chamber_mass_flow_kg_s,inputs.inlet_pressure_pa,inputs.injector_pressure_drop_fraction))
    required = pressure*(1+fraction)
    margin = inlet-required
    required_output = _finite_rational(required,'Required feed pressure')
    try:
        margin_output = _finite_rational(margin,'Inlet pressure margin')
    except DomainError:
        margin_output = None
    values = {'cycle': inputs.cycle,'chamber_pressure_pa': inputs.chamber_pressure_pa,
        'required_feed_pressure_pa': required_output,
        'inlet_pressure_margin_pa': margin_output,'inlet_pressure_margin_representable': margin_output is not None,
        'chamber_mass_flow_kg_s': inputs.chamber_mass_flow_kg_s}
    if inputs.cycle == 'pressure-fed':
        warnings = (() if margin >= 0 else (WarningMessage('insufficient_pressure_fed_supply',
            'Declared inlet pressure is below the injector-feed requirement; this is sizing, not achieved delivery.',
            severity='critical',field='inlet_pressure_pa'),))
        values.update(pump_pressure_rise_pa=0.,pump_power_w=0.,turbine_specific_work_j_kg=0.,
            turbine_mass_flow_kg_s=0.,turbine_flow_fraction=0.,turbine_power_w=0.,power_balance_residual_w=0.,
            total_propellant_flow_kg_s=inputs.chamber_mass_flow_kg_s,delivered_mass_fraction=1.,
            pump_mass_flow_kg_s=0.,pressure_fed_supply_sufficient=margin>=0,
            energy_balance_scope='pressure_fed_requirements_sizing_only',warnings=warnings,
            assumptions=('pressure-fed requirements sizing; sufficiency covers declared injector loss only',
                         'feed-line loss, transient tank pressure, cavitation and hardware capability not modelled'))
        return values
    rise = max(Fraction(0),-margin)
    specific = rise/(Fraction(inputs.propellant_density_kg_m3)*Fraction(inputs.pump_efficiency))
    coefficient = Fraction(inputs.turbine_efficiency)*Fraction(inputs.turbine_specific_heat_j_kg_k)*Fraction(inputs.turbine_inlet_temperature_k)
    gamma = Fraction(inputs.turbine_gamma)
    exponent = (gamma-1)/gamma
    outputs, high_fraction = component_power_outputs(inputs.cycle,inputs.pump_feed_policy,flow,specific,
        coefficient,exponent,Fraction(inputs.turbine_pressure_ratio))
    warnings = [WarningMessage('hot_gas_pressure_not_closed',
        'Hot-gas supply/exhaust pressure compatibility and integrated cycle feasibility are not established.')]
    if inputs.pump_feed_policy == 'external-gg-feed':
        scope = 'ideal_main_pump_turbine_balance_excludes_external_gg_supply'
        warnings.append(WarningMessage('external_gg_feed_not_modelled',
            'Generator propellant is supplied separately; its pressure source and energy duty are excluded.',severity='critical'))
    else:
        scope = 'ideal_shared_pump_and_turbine_balance'
    if high_fraction and inputs.cycle == 'gas-generator':
        warnings.append(WarningMessage('high_gas_generator_flow_fraction',
            'Computed generator fraction triggers the heuristic 8% advisory; this is not an acceptance limit.'))
    values.update(outputs,pump_pressure_rise_pa=_finite_rational(rise,'Pump pressure rise'),
        pressure_fed_supply_sufficient=None,energy_balance_scope=scope,warnings=tuple(warnings),
        assumptions=('single equivalent liquid density and common cold-side pressure rise',
            'constant efficiencies and calorically perfect turbine gas; total inlet temperature and total inlet/outlet pressure ratio',
            'no pump map, cavitation, gearbox, bearing or mechanical loss model',
            'generator turbine flow is dumped; staged turbine flow returns to the chamber',
            'zero power residual is an algebraic model identity, not measured closure'))
    return values


def calculate_engine_cycle(
    *, cycle: str, chamber_pressure_pa: float, chamber_mass_flow_kg_s: float,
    inlet_pressure_pa: float, injector_pressure_drop_fraction: float,
    propellant_density_kg_m3: float | None = None, pump_efficiency: float | None = None,
    turbine_efficiency: float | None = None, turbine_inlet_temperature_k: float | None = None,
    turbine_specific_heat_j_kg_k: float | None = None, turbine_gamma: float | None = None,
    turbine_pressure_ratio: float | None = None, pump_feed_policy: str = 'shared-pump',
) -> EngineCycleResult:
    """Size pressure-fed feed requirements or close an explicitly routed shaft balance.

    Parameters
    ----------
    cycle : {'pressure-fed', 'gas-generator', 'staged-combustion'}
        Cycle architecture.
    chamber_pressure_pa : float
        Chamber pressure, Pa; positive.
    chamber_mass_flow_kg_s : float
        Main-chamber propellant flow, kg/s; positive.
    inlet_pressure_pa : float
        Supply pressure, Pa; non-negative.
    injector_pressure_drop_fraction : float
        Injector pressure drop as a fraction of chamber pressure, in (0, 1].
    propellant_density_kg_m3, pump_efficiency, turbine_efficiency : float, optional
        Pump and turbine inputs; required for pumped cycles.
    turbine_inlet_temperature_k, turbine_specific_heat_j_kg_k : float, optional
        Turbine total inlet temperature (K) and specific heat (J/(kg K)).
    turbine_gamma, turbine_pressure_ratio : float, optional
        Turbine-gas gamma and total inlet/outlet pressure ratio; both above one.
    pump_feed_policy : {'shared-pump', 'external-gg-feed'}, default 'shared-pump'
        Generator propellant routing; ``'external-gg-feed'`` applies only to the
        gas-generator cycle.

    Returns
    -------
    EngineCycleResult
        Feed-pressure requirement, pump/turbine powers, routed flows, warnings and
        assumptions.

    Raises
    ------
    DomainError
        If an input is out of domain, the routing cannot close the power balance
        (code ``infeasible_cycle_power_balance``), or an output is not representable.
    ConvergenceError
        If an enclosed turbine expression is unresolved within the precision budget.

    Notes
    -----
    In the shared-pump gas-generator routing the generator flow is included in the
    pump duty; ``'external-gg-feed'`` keeps the separate, unmodelled generator supply
    and its energy duty explicit [1]_. Generator turbine flow is dumped; staged
    turbine flow returns to the chamber. Efficiencies are constant and the turbine gas
    is calorically perfect [2]_. Neither routing establishes an operable integrated
    engine or closes hot-gas pressure compatibility, pump maps, cavitation or hardware
    capability.

    References
    ----------
    .. [1] NASA-TMX64749-CYCLE. R. H. Schmucker, 1973, equivalent cycle power balances
           and whole pumped mass flows; theory source, not hardware validation.
    .. [2] NASA-ENGINESIM-THEORY. NASA Glenn Research Center, EngineSim theory notes:
           ideal gas turbine temperature/pressure and shaft-work relations.
    .. [3] SUTTON-BIBLARZ-2017. G. P. Sutton and O. Biblarz, "Rocket Propulsion
           Elements", 9th ed., Wiley, 2017.
    """
    source = EngineCycleInputs(cycle,chamber_pressure_pa,chamber_mass_flow_kg_s,inlet_pressure_pa,
        injector_pressure_drop_fraction,propellant_density_kg_m3,pump_efficiency,turbine_efficiency,
        turbine_inlet_temperature_k,turbine_specific_heat_j_kg_k,turbine_gamma,turbine_pressure_ratio,pump_feed_policy)
    return EngineCycleResult(**_cycle_values(source),inputs=source)
