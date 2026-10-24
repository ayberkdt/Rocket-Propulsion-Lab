"""Certified finite-grid allocations and independently rounded stage ledgers.

The log of each stage mass factor is strictly convex in delta-v. Selecting the
least next grid increment therefore minimizes the separable objective over all
integer allocations. Comparisons use unrounded expressions and directed bounds;
an unresolved comparison fails typed rather than guessing a winner.
"""
from decimal import Decimal
from fractions import Fraction
from functools import lru_cache

from rocket_propulsion.core.constants import STANDARD_GRAVITY_M_S2
from rocket_propulsion.core.errors import ConvergenceError, DomainError

from .burns.marginal_numerics import _PRECISIONS, _Intervals
from .performance_numerics import _exp_interval, _positive_enclosed_output

_ONE = (Decimal(1), Decimal(1))
_ZERO = (Decimal(0), Decimal(0))
_G = Fraction(STANDARD_GRAVITY_M_S2)
_STAGING_POLICY = 'exact_grid_convex_increment_comparisons_and_enclosed_mass_ledgers_v1'


def _exponent(stage, dv):
    return dv/(_G*Fraction(stage.specific_impulse_s))


@lru_cache(maxsize=16384)
def _feasible(stage, dv):
    x = _exponent(stage, dv)
    # exp(-x) < 2**-2048 for x > 2048, below every positive binary64 eps.
    if x > 2048:
        return False
    for precision in _PRECISIONS:
        lo, hi = _exp_interval(-x, precision)
        eps = Fraction(stage.structural_fraction)
        if Fraction(lo) > eps:
            return True
        if Fraction(hi) <= eps:
            return False
    raise ConvergenceError('Stage capacity comparison unresolved.', field='Stage capacity')


def _state(stage, dv, precision):
    """Survival and positive denominator; refine before any interval division."""
    x = _exponent(stage, dv)
    for p in _PRECISIONS:
        if p < precision:
            continue
        math = _Intervals(p)
        survival = _exp_interval(-x, p)
        denominator = math.sub(survival, math.rational(Fraction(stage.structural_fraction)))
        if denominator[0] > 0:
            return survival, denominator
    raise ConvergenceError('Positive stage denominator unresolved.', field='Stage denominator')


def _increment_identity(a, ak, b, bk, step):
    """Recognize exact equal cross-products by collecting rational exponents."""
    ax, bx = _exponent(a, step), _exponent(b, step)
    ae, be = Fraction(a.structural_fraction), Fraction(b.structural_fraction)
    terms = {}
    for power, coefficient in ((-ak*ax-(bk+1)*bx, Fraction(1)),
            (-ak*ax, -be), (-(bk+1)*bx, -ae),
            (-bk*bx-(ak+1)*ax, Fraction(-1)),
            (-bk*bx, ae), (-(ak+1)*ax, be)):
        terms[power] = terms.get(power, Fraction(0))+coefficient
    return not any(terms.values())


def _compare_increment(a, ak, b, bk, step):
    # The common (1-eps) numerator cancels in each stage's next/current ratio.
    if _increment_identity(a, ak, b, bk, step):
        return 0
    for precision in _PRECISIONS:
        math = _Intervals(precision)
        an = math.sub(_exp_interval(-_exponent(a, step*ak), precision),
                      math.rational(Fraction(a.structural_fraction)))
        ad = math.sub(_exp_interval(-_exponent(a, step*(ak+1)), precision),
                      math.rational(Fraction(a.structural_fraction)))
        bn = math.sub(_exp_interval(-_exponent(b, step*bk), precision),
                      math.rational(Fraction(b.structural_fraction)))
        bd = math.sub(_exp_interval(-_exponent(b, step*(bk+1)), precision),
                      math.rational(Fraction(b.structural_fraction)))
        if min(ad[0], bd[0]) <= 0:
            continue
        difference = math.sub(math.mul(an, bd), math.mul(bn, ad))
        if difference[1] < 0:
            return -1
        if difference[0] > 0:
            return 1
    raise ConvergenceError('Grid increment ordering unresolved.', field='Grid optimum')


@lru_cache(maxsize=64)
def _allocate(stages, target, resolution):
    step = Fraction(target)/resolution
    counts = [0]*len(stages)
    if len(stages) == 1:
        if _feasible(stages[0], Fraction(target)):
            return (resolution,)
    else:
        for _ in range(resolution):
            chosen = None
            for index, stage in enumerate(stages):
                if not _feasible(stage, step*(counts[index]+1)):
                    continue
                if chosen is None or _compare_increment(stage, counts[index],
                        stages[chosen], counts[chosen], step) < 0:
                    chosen = index
            if chosen is None:
                break
            counts[chosen] += 1
        else:
            return tuple(counts)
    raise DomainError('No finite-mass allocation exists on the requested delta-v grid; '
                      'this does not by itself prove continuous-model infeasibility.',
                      code='infeasible_staging_grid', field='resolution')


@lru_cache(maxsize=256)
def _ledger(inputs, allocations, precision):
    math = _Intervals(precision)
    upper = math.rational(Fraction(inputs.payload_mass_kg)); stages = []
    step = Fraction(inputs.target_delta_v_m_s)/inputs.resolution
    for stage, count in zip(reversed(inputs.stages), reversed(allocations)):
        eps = math.rational(Fraction(stage.structural_fraction))
        survival, denominator = _state(stage, step*count, precision)
        wet = math.div(math.mul(upper, math.sub(_ONE, survival)), denominator)
        initial = math.div(math.mul(upper, math.sub(_ONE, eps)), denominator)
        stages.append({'initial_mass_kg': initial, 'upper_stack_mass_kg': upper,
            'dry_mass_kg': math.mul(eps, wet), 'propellant_mass_kg': math.mul(math.sub(_ONE, eps), wet),
            'mass_ratio': _exp_interval(_exponent(stage, step*count), precision)})
        upper = initial
    return upper, tuple(reversed(stages))


def _objective_polynomial(inputs, allocations):
    """Normalize the product denominator; identical polynomials imply equal mass."""
    terms = {Fraction(0): Fraction(1)}
    step = Fraction(inputs.target_delta_v_m_s)/inputs.resolution
    for stage, count in zip(inputs.stages, allocations):
        eps = Fraction(stage.structural_fraction)
        x = _exponent(stage, step*count); following = {}
        for power, coefficient in terms.items():
            following[power-x] = following.get(power-x, Fraction(0))+coefficient/(1-eps)
            following[power] = following.get(power, Fraction(0))-coefficient*eps/(1-eps)
        terms = {power: coefficient for power, coefficient in following.items() if coefficient}
    return terms


def _grid_difference(coarse, counts, refined, refined_counts):
    if _objective_polynomial(coarse, counts) == _objective_polynomial(refined, refined_counts):
        return 0.
    def evaluate(precision):
        math = _Intervals(precision)
        a = _ledger(coarse, counts, precision)[0]
        b = _ledger(refined, refined_counts, precision)[0]
        lo, hi = math.sub(a, b)
        absolute = (lo, hi) if lo >= 0 else (math.down.minus(hi), math.up.minus(lo)) if hi <= 0 else (Decimal(0), max(math.up.minus(lo), hi))
        return math.div(absolute, b)
    return _positive_enclosed_output(evaluate, 'Grid mass difference')


def _exports(inputs, allocations):
    result = []
    for index, (stage, count) in enumerate(zip(inputs.stages, allocations)):
        fields = {}
        for name in ('initial_mass_kg', 'upper_stack_mass_kg', 'dry_mass_kg',
                     'propellant_mass_kg', 'mass_ratio'):
            fields[name] = (0. if count == 0 and name in ('dry_mass_kg', 'propellant_mass_kg') else
                _positive_enclosed_output(lambda p, *, index=index, name=name: _ledger(inputs, allocations, p)[1][index][name], name))
        result.append(fields)
    initial = _positive_enclosed_output(lambda p: _ledger(inputs, allocations, p)[0], 'Initial mass')
    def fraction(p):
        math = _Intervals(p)
        return math.div(math.rational(Fraction(inputs.payload_mass_kg)), _ledger(inputs, allocations, p)[0])
    return initial, _positive_enclosed_output(fraction, 'Payload fraction'), tuple(result)
