"""Exact pump-head ledgers and enclosed turbine/shaft component expressions.

Pump specific work and the routed mass flows are rational in the inputs and are
carried exactly. Turbine specific work contains a rational power of the pressure
ratio; it is evaluated exactly when that power is rational and otherwise enclosed
with directed decimal bounds at increasing precision.

References
----------
.. [1] NASA-TMX64749-CYCLE. R. H. Schmucker, 1973, equivalent cycle power balances
       and whole pumped mass flows; theory source, not hardware validation.
.. [2] NASA-ENGINESIM-THEORY. NASA Glenn Research Center, EngineSim theory notes:
       ideal gas turbine temperature/pressure and shaft-work relations.
.. [3] PYTHON-3.12-DECIMAL. Python 3.12 documentation, decimal module: precision
       contexts and exact conversion of binary floats.
"""
from decimal import Decimal
from fractions import Fraction
from functools import lru_cache

from rocket_propulsion.core.errors import ConvergenceError, DomainError

from .burns.marginal_numerics import _PRECISIONS, _Intervals, _log_interval
from .burns.numeric_outputs import _finite_rational
from .performance_numerics import _positive_enclosed_output

_CYCLE_POLICY = 'exact_binary_pressure_ledgers_and_enclosed_component_cycle_v1'


def _integer_root(value, degree):
    """Return the exact integer ``degree``-th root of ``value``, or None.

    Parameters
    ----------
    value : int
        Positive integer.
    degree : int
        Positive root degree.

    Returns
    -------
    int or None
        ``r`` with ``r**degree == value``, or None if no integer root exists or the
        degree exceeds the bit length of ``value``.
    """
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
    """Return the turbine specific work exactly when it is rational, else None.

    Parameters
    ----------
    coefficient : Fraction
        eta_t * cp * T_in, J/kg.
    exponent : Fraction
        (gamma - 1) / gamma.
    ratio : Fraction
        Total inlet/outlet pressure ratio, greater than one.

    Returns
    -------
    Fraction or None
        w_t = coefficient * (1 - ratio**(-exponent)) when ``ratio**exponent`` is
        rational, otherwise None.

    Notes
    -----
    Ideal-gas turbine work with an isentropic efficiency [1]_. For coprime
    exponent p/q the power is rational only if numerator and denominator of the ratio
    are perfect q-th powers.

    References
    ----------
    .. [1] NASA-ENGINESIM-THEORY. NASA Glenn Research Center, EngineSim theory notes:
           ideal gas turbine temperature/pressure and shaft-work relations.
    """
    degree = exponent.denominator
    numerator = _integer_root(ratio.numerator, degree)
    denominator = _integer_root(ratio.denominator, degree)
    if numerator is None or denominator is None:return None
    survival = Fraction(denominator, numerator)**exponent.numerator
    return coefficient*(1-survival)


@lru_cache(maxsize=2048)
def _work_bounds(coefficient, exponent, ratio, precision):
    """Enclose the turbine specific work at the given decimal precision.

    Parameters
    ----------
    coefficient, exponent, ratio : Fraction
        As in :func:`_exact_turbine_work`.
    precision : int
        Decimal working precision.

    Returns
    -------
    tuple of Decimal
        Lower and upper bounds of coefficient * (1 - exp(-exponent * ln(ratio))).

    References
    ----------
    .. [1] NASA-ENGINESIM-THEORY. NASA Glenn Research Center, EngineSim theory notes:
           ideal gas turbine temperature/pressure and shaft-work relations.
    .. [2] PYTHON-3.12-DECIMAL. Python 3.12 documentation, decimal module: precision
           contexts and exact conversion of binary floats.
    """
    math = _Intervals(precision)
    power = math.mul(math.rational(exponent), _log_interval(ratio, precision))
    negative = math.sub((Decimal(0), Decimal(0)), power)
    survival = (math.near.next_minus(math.near.exp(negative[0])),
                min(Decimal(1), math.near.next_plus(math.near.exp(negative[1]))))
    complement = math.sub((Decimal(1), Decimal(1)), survival)
    return math.mul(math.rational(coefficient), complement)


def _work_relation(specific, coefficient, exponent, ratio):
    """Compare turbine specific work with pump specific work.

    Parameters
    ----------
    specific : Fraction
        Pump specific work Delta_p / (rho * eta_p), J/kg.
    coefficient, exponent, ratio : Fraction
        Turbine work parameters as in :func:`_exact_turbine_work`.

    Returns
    -------
    tuple
        ``(sign, precision)``: +1 if turbine work exceeds pump work, -1 if it is
        smaller, and the first precision that resolved the comparison. Zero pump
        work returns +1.

    Raises
    ------
    ConvergenceError
        If the comparison is unresolved within the precision budget.

    References
    ----------
    .. [1] NASA-TMX64749-CYCLE. R. H. Schmucker, 1973, equivalent cycle power balances
           and whole pumped mass flows; theory source, not hardware validation.
    .. [2] PYTHON-3.12-DECIMAL. Python 3.12 documentation, decimal module: precision
           contexts and exact conversion of binary floats.
    """
    if not specific:return 1, _PRECISIONS[0]
    for precision in _PRECISIONS:
        math = _Intervals(precision)
        difference = math.sub(_work_bounds(coefficient, exponent, ratio, precision),
                              math.rational(specific))
        if difference[0] > 0:return 1, precision
        if difference[1] < 0:return -1, precision
    raise ConvergenceError('Turbine/pump capacity comparison is unresolved.', field='Cycle power balance')


def _power_feasible(cycle, policy, relation):
    """Refuse routings whose turbine work cannot close the shaft power balance.

    Parameters
    ----------
    cycle : str
        Normalized cycle name.
    policy : str
        Pump feed policy.
    relation : int
        Sign of turbine work minus pump specific work.

    Raises
    ------
    DomainError
        ``infeasible_cycle_power_balance`` when a shared-pump gas generator has
        turbine work not exceeding pump work (the routed flow would be unbounded), or
        a staged-combustion cycle has turbine work below pump work.

    References
    ----------
    .. [1] NASA-TMX64749-CYCLE. R. H. Schmucker, 1973, equivalent cycle power balances
           and whole pumped mass flows; theory source, not hardware validation.
    """
    if ((cycle == 'gas-generator' and policy == 'shared-pump' and relation <= 0)
            or (cycle == 'staged-combustion' and relation < 0)):
        raise DomainError('Available turbine work cannot close the declared pump routing.',
                          code='infeasible_cycle_power_balance', field='Cycle power balance')


def _fraction_outputs(cycle, policy, flow, specific, work):
    """Return exact routed mass flows and pump power for the shaft balance.

    Parameters
    ----------
    cycle, policy : str
        Normalized cycle name and pump feed policy.
    flow : Fraction
        Chamber mass flow m_c, kg/s.
    specific : Fraction
        Pump specific work s, J/kg.
    work : Fraction
        Turbine specific work w, J/kg.

    Returns
    -------
    dict of Fraction
        Turbine, total, pumped flows, turbine flow fraction, delivered fraction and
        pump power.

    Notes
    -----
    The shaft balance m_t * w = m_pump * s is closed for the declared routing [1]_.
    Shared-pump gas generator: the pump also feeds the generator, so
    m_t = m_c * s / (w - s) and m_total = m_pump = m_c * w / (w - s). External
    generator feed: m_t = m_c * s / w, m_total = m_c + m_t, m_pump = m_c. Staged
    combustion: the turbine flow returns to the chamber, m_t = m_c * s / w and
    m_total = m_pump = m_c.

    References
    ----------
    .. [1] NASA-TMX64749-CYCLE. R. H. Schmucker, 1973, equivalent cycle power balances
           and whole pumped mass flows; theory source, not hardware validation.
    .. [2] SUTTON-BIBLARZ-2017. G. P. Sutton and O. Biblarz, "Rocket Propulsion
           Elements", 9th ed., Wiley, 2017.
    """
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
    """Enclose the routed flows and pump power when turbine work is an interval.

    Parameters
    ----------
    cycle, policy : str
        Normalized cycle name and pump feed policy.
    flow, specific : Fraction
        Chamber mass flow and pump specific work.
    work : tuple of Decimal
        Enclosure of the turbine specific work.
    precision : int
        Decimal working precision.

    Returns
    -------
    dict
        Decimal enclosures of the quantities returned by :func:`_fraction_outputs`.

    References
    ----------
    .. [1] NASA-TMX64749-CYCLE. R. H. Schmucker, 1973, equivalent cycle power balances
           and whole pumped mass flows; theory source, not hardware validation.
    .. [2] PYTHON-3.12-DECIMAL. Python 3.12 documentation, decimal module: precision
           contexts and exact conversion of binary floats.
    """
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
    """Close the pump/turbine shaft balance for a pumped cycle routing.

    Parameters
    ----------
    cycle : {'gas-generator', 'staged-combustion'}
        Normalized cycle name.
    policy : {'shared-pump', 'external-gg-feed'}
        Pump feed policy.
    flow : Fraction
        Chamber mass flow, kg/s.
    specific : Fraction
        Pump specific work Delta_p / (rho * eta_p), J/kg.
    coefficient : Fraction
        eta_t * cp * T_in, J/kg.
    exponent : Fraction
        (gamma - 1) / gamma.
    ratio : Fraction
        Turbine total inlet/outlet pressure ratio.

    Returns
    -------
    values : dict
        Turbine specific work, turbine/total/pumped mass flows, turbine flow
        fraction, delivered fraction, pump power, turbine power and the zero power
        residual, each rounded once to binary64.
    high_fraction : bool
        True when the turbine flow fraction exceeds the heuristic advisory threshold
        of 0.08. The threshold is a display advisory, not an acceptance limit, and does
        not enter any mass or power equation.

    Raises
    ------
    DomainError
        If the routing cannot close the power balance or an output is not
        representable.
    ConvergenceError
        If an enclosed output or comparison is unresolved within the precision
        budget.

    Notes
    -----
    Pump power P = m_pump * s with s = Delta_p / (rho * eta_p) [3]_. Turbine power
    equals pump power by construction, so the zero residual is an algebraic identity,
    not measured closure [1]_. Turbine specific work
    w_t = eta_t * cp * T_in * (1 - PR**(-(gamma - 1) / gamma)) [2]_. Mechanical,
    bearing and gearbox losses are not modelled.

    References
    ----------
    .. [1] NASA-TMX64749-CYCLE. R. H. Schmucker, 1973, equivalent cycle power balances
           and whole pumped mass flows; theory source, not hardware validation.
    .. [2] NASA-ENGINESIM-THEORY. NASA Glenn Research Center, EngineSim theory notes:
           ideal gas turbine temperature/pressure and shaft-work relations.
    .. [3] SUTTON-BIBLARZ-2017. G. P. Sutton and O. Biblarz, "Rocket Propulsion
           Elements", 9th ed., Wiley, 2017.
    """
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
