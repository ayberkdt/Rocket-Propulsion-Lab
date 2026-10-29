"""Exact composition ledgers and bounded logarithmic exports for frozen mixtures.

These are numerical bounds on the supplied ideal-gas fits, not physical error
bounds. Log bases are made pairwise coprime by GCD splitting, without factoring.
This preserves exact cancellations such as 3*ln(2)-ln(8).
"""
from decimal import MAX_EMAX, MIN_EMIN, ROUND_HALF_EVEN, Context, Decimal
from fractions import Fraction as F
from functools import lru_cache
from math import gcd, inf

from rocket_propulsion.core.errors import ConvergenceError, DomainError

from .nasa_polynomials import UNIVERSAL_GAS_CONSTANT_J_MOL_K, _rounded


def _rational(value, label):
    result = _rounded(value)
    if result in (inf, -inf) or (value != 0 and result == 0):
        raise DomainError('Mixture quantity is outside the nonzero finite output range.',
                          code='mixture_output_range', field=label)
    return result


@lru_cache(maxsize=1024)
def _composition(species, values, basis):
    """Normalize accepted raw values once; never use component display fractions.

    Arguments are immutable tuples, so equal requests share one exact ledger.
    """
    raw = tuple(map(F, values))
    total = sum(raw, F(0))
    if any(x <= 0 for x in raw) or abs(total-1) > F(1e-9):
        raise DomainError('Composition fractions must sum to one within 1e-9.',
                          code='composition_sum', field='input_fractions')
    amounts = raw if basis == 'mole' else tuple(v/F(s.molecular_mass_kg_mol)
                                               for s,v in zip(species,raw,strict=True))
    amount = sum(amounts)
    mole = tuple(n/amount for n in amounts)
    mean = sum(x*F(s.molecular_mass_kg_mol) for s,x in zip(species,mole,strict=True))
    mass = tuple(x*F(s.molecular_mass_kg_mol)/mean for s,x in zip(species,mole,strict=True))
    return mole, mass, mean


def _coalesced(species, values):
    """Combine aliases of one species exactly, in first-occurrence order."""
    order, amounts = [], {}
    for item, value in zip(species, values, strict=True):
        key = item.key.upper()
        if key not in amounts:
            order.append(item)
            amounts[key] = F(0)
        amounts[key] += F(value)
    return tuple(order), tuple(amounts[item.key.upper()] for item in order)


def _ledger(mixture):
    """Return owned species and the exact composition of retained raw inputs."""
    species = tuple(c.species for c in mixture.components)
    index = {item.key.upper(): item for item in species}
    order, values = _coalesced(tuple(index[key.upper()] for key, _ in mixture.input_fractions),
                               [v for _, v in mixture.input_fractions])
    return order, _composition(order, values, mixture.composition_basis)


def _log_basis(terms):
    """Reduce rational log arguments to pairwise-coprime integer bases."""
    ledger: dict[int, F] = {}
    def add(base, coefficient):
        if base > 1 and coefficient:
            ledger[base] = ledger.get(base,F(0))+coefficient
            if not ledger[base]: del ledger[base]
    for coefficient,argument in terms:
        add(argument.numerator,coefficient)
        add(argument.denominator,-coefficient)
    for _ in range(4096):
        bases = sorted(ledger)
        split = None
        for i,a in enumerate(bases):
            for b in bases[i+1:]:
                common = gcd(a,b)
                if common > 1:
                    split = (a,b,common)
                    break
            if split: break
        if split is None: return tuple(sorted(ledger.items()))
        a,b,common = split
        ca,cb = ledger.pop(a),ledger.pop(b)
        add(common,ca+cb);add(a//common,ca);add(b//common,cb)
    raise ConvergenceError('Mixture logarithm basis exceeded its split budget.',
                           code='mixture_log_resource')


@lru_cache(maxsize=512)
def _integer_log(base, precision):
    context = Context(prec=precision, Emax=MAX_EMAX, Emin=MIN_EMIN,
                      rounding=ROUND_HALF_EVEN, traps=[], flags=[], capitals=1, clamp=0)
    value = context.ln(Decimal(base))
    unit = F(Decimal((0,(1,),value.adjusted()-precision+1)))
    return F(value),unit


def _log_export(base, terms, label):
    basis = _log_basis(terms)
    if not basis: return _rational(base,label)
    for precision in (48,96,192,384,768,1536,3072):
        center,radius = base,F(0)
        for integer,coefficient in basis:
            logarithm,bound = _integer_log(integer,precision)
            center += coefficient*logarithm
            radius += abs(coefficient)*bound
        lower,upper = _rounded(center-radius),_rounded(center+radius)
        if lower == upper:
            # Nonempty coprime integer log bases with rational coefficients
            # cannot cancel a rational base to zero. Erased nonzero is refused.
            if lower in (inf,-inf) or lower == 0:
                raise DomainError('Mixture quantity is outside the nonzero finite output range.',
                                  code='mixture_output_range',field=label)
            return lower
    raise ConvergenceError('Mixture logarithm rounding is unresolved.',
                           code='mixture_rounding_unresolved',field=label)


def _forward(mixture, temperature, pressure):
    species,(mole,_,mean) = _ledger(mixture)
    r = F(UNIVERSAL_GAS_CONSTANT_J_MOL_K)
    cp,h,hlog,s,slog = F(0),F(0),F(0),F(0),F(0)
    for item,x in zip(species,mole,strict=True):
        (c,_),(a,b),(d,e) = item._interval(temperature)._terms(temperature)
        if c <= r:
            raise DomainError('Mixture component has nonpositive cv.',
                              code='invalid_thermodynamic_property',field=item.key)
        cp += x*c;h += x*a;hlog += x*b;s += x*d;slog += x*e
    entropy_logs = [(slog/mean,F(temperature)),(-r/mean,F(pressure)/100_000)]
    entropy_logs.extend((-r*x/mean,x) for x in mole)
    return {'molecular_mass_kg_mol': _rational(mean,'molecular_mass_kg_mol'),
                'gas_constant_j_kg_k': _rational(r/mean,'gas_constant_j_kg_k'),
                'cp_j_kg_k': _rational(cp/mean,'cp_j_kg_k'),
                'cv_j_kg_k': _rational((cp-r)/mean,'cv_j_kg_k'),
                'gamma': _rational(cp/(cp-r),'gamma'),
                'enthalpy_j_kg': _log_export(h/mean,[(hlog/mean,F(temperature))],'enthalpy_j_kg'),
                'entropy_j_kg_k': _log_export(s/mean,entropy_logs,'entropy_j_kg_k')}
