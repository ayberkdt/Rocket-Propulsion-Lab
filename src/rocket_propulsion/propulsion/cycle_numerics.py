"""Exact head ledgers and finite-budget turbine/shaft component expressions."""
from decimal import Decimal
from fractions import Fraction
from functools import lru_cache

from rocket_propulsion.core.errors import ConvergenceError, DomainError

from .burns.marginal_numerics import _PRECISIONS, _Intervals, _log_interval
from .burns.numeric_outputs import _finite_rational
from .performance_numerics import _positive_enclosed_output

_CYCLE_POLICY = 'exact_binary_pressure_ledgers_and_enclosed_component_cycle_v1'


def _integer_root(value, degree):
    if value == 1 or degree == 1:
        return value
    if degree > value.bit_length():
        return None
    low, high = 1, 1 << ((value.bit_length()+degree-1)//degree)
    while low <= high:
        middle = (low+high)//2
        power = middle**degree
        if power == value:return middle
        if power < value:low = middle+1
        else:high = middle-1
    return None


def _exact_turbine_work(coefficient, exponent, ratio):
    # For coprime p/q, rational r**(p/q) requires both rational-base
    # integers to be perfect q-th powers. Large q is rejected by bit length.
    degree = exponent.denominator
    numerator = _integer_root(ratio.numerator, degree)
    denominator = _integer_root(ratio.denominator, degree)
    if numerator is None or denominator is None:return None
    survival = Fraction(denominator, numerator)**exponent.numerator
    return coefficient*(1-survival)


@lru_cache(maxsize=2048)
def _work_bounds(coefficient, exponent, ratio, precision):
    math = _Intervals(precision)
    power = math.mul(math.rational(exponent), _log_interval(ratio, precision))
    negative = math.sub((Decimal(0), Decimal(0)), power)
    survival = (math.near.next_minus(math.near.exp(negative[0])),
                min(Decimal(1), math.near.next_plus(math.near.exp(negative[1]))))
    complement = math.sub((Decimal(1), Decimal(1)), survival)
    return math.mul(math.rational(coefficient), complement)


def _work_relation(specific, coefficient, exponent, ratio):
    if not specific:return 1, _PRECISIONS[0]
    for precision in _PRECISIONS:
        math = _Intervals(precision)
        difference = math.sub(_work_bounds(coefficient, exponent, ratio, precision),
                              math.rational(specific))
        if difference[0] > 0:return 1, precision
        if difference[1] < 0:return -1, precision
    raise ConvergenceError('Turbine/pump capacity comparison is unresolved.', field='Cycle power balance')


def _power_feasible(cycle, policy, relation):
    if ((cycle == 'gas-generator' and policy == 'shared-pump' and relation <= 0)
            or (cycle == 'staged-combustion' and relation < 0)):
        raise DomainError('Available turbine work cannot close the declared pump routing.',
                          code='infeasible_cycle_power_balance', field='Cycle power balance')


def _fraction_outputs(cycle, policy, flow, specific, work):
    if not specific:
        return {'turbine_mass_flow_kg_s': Fraction(0),'total_propellant_flow_kg_s': flow,
                    'pump_mass_flow_kg_s': flow,'turbine_flow_fraction': Fraction(0),
                    'delivered_mass_fraction': Fraction(1),'pump_power_w': Fraction(0)}
    if cycle == 'gas-generator' and policy == 'shared-pump':
        turbine = flow*specific/(work-specific)
        total = flow*work/(work-specific)
        pumped = total
        fraction = specific/work
        delivered = 1-fraction
    else:
        turbine = flow*specific/work
        total = flow+turbine if cycle == 'gas-generator' else flow
        pumped = flow
        fraction = turbine/total
        delivered = flow/total
    return {'turbine_mass_flow_kg_s': turbine,'total_propellant_flow_kg_s': total,
                'pump_mass_flow_kg_s': pumped,'turbine_flow_fraction': fraction,
                'delivered_mass_fraction': delivered,'pump_power_w': pumped*specific}


def _dependent_bounds(cycle, policy, flow, specific, work, precision):
    math = _Intervals(precision)
    f, s = math.rational(flow), math.rational(specific)
    one = (Decimal(1), Decimal(1))
    if cycle == 'gas-generator' and policy == 'shared-pump':
        denominator = math.sub(work, s)
        turbine = math.div(math.mul(f, s), denominator)
        total = math.div(math.mul(f, work), denominator)
        pumped = total
        fraction = math.div(s, work)
        delivered = math.sub(one, fraction)
    else:
        turbine = math.div(math.mul(f, s), work)
        total = math.add(f, turbine) if cycle == 'gas-generator' else f
        pumped = f
        fraction = math.div(turbine, total)
        delivered = math.div(f, total)
    return {'turbine_mass_flow_kg_s': turbine,'total_propellant_flow_kg_s': total,
                'pump_mass_flow_kg_s': pumped,'turbine_flow_fraction': fraction,
                'delivered_mass_fraction': delivered,'pump_power_w': math.mul(pumped, s)}


def component_power_outputs(cycle, policy, flow, specific, coefficient, exponent, ratio):
    exact_work = _exact_turbine_work(coefficient, exponent, ratio)
    if exact_work is not None:
        relation = (exact_work > specific)-(exact_work < specific)
        _power_feasible(cycle, policy, relation)
        work = _finite_rational(exact_work, 'Turbine specific work')
        values = {name:_finite_rational(value, name) for name,value in
                  _fraction_outputs(cycle,policy,flow,specific,exact_work).items()}
        high_fraction = values['turbine_flow_fraction'] > .08
    else:
        relation, minimum_precision = _work_relation(specific,coefficient,exponent,ratio)
        _power_feasible(cycle,policy,relation)
        work = _positive_enclosed_output(lambda p:_work_bounds(coefficient,exponent,ratio,p),
                                         'Turbine specific work')
        if not specific:
            values = {name:float(value) for name,value in _fraction_outputs(cycle,policy,flow,specific,Fraction(1)).items()}
            high_fraction = False
        else:
            def bounds(precision):
                precision = max(precision, minimum_precision)
                return _dependent_bounds(cycle,policy,flow,specific,
                    _work_bounds(coefficient,exponent,ratio,precision),precision)
            names = ('turbine_mass_flow_kg_s','total_propellant_flow_kg_s','pump_mass_flow_kg_s',
                     'turbine_flow_fraction','delivered_mass_fraction','pump_power_w')
            values = {name:_positive_enclosed_output(lambda p, *, name=name:bounds(p)[name],name) for name in names}
            # The warning is an explicitly heuristic display threshold. The
            # rounded coefficient does not enter any mass/power equation.
            high_fraction = values['turbine_flow_fraction'] > .08
    values.update(turbine_specific_work_j_kg=work,turbine_power_w=values['pump_power_w'],
                  power_balance_residual_w=0.)
    return values, high_fraction
