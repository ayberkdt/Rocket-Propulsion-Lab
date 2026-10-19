"""Exact rational phase ledger with a single final binary conversion.

References
----------
BIPM-SI: disjoint phase energy and watt=joule/second.
NASA-7009B: numerical implementation limits, distinct from physical uncertainty.
"""

from fractions import Fraction

from .rational_numerics import binary_ledger


def solve_power_budget(firing: float, on_time: float, duty: float,
                       standby: float, warmup: float, warmup_time: float) -> tuple:
    """Evaluate validated binary inputs exactly before converting final channels.

    Standby duration is t*(1-duty)/duty. Subtracting two rounded elapsed times
    can erase or double this phase near continuous duty. Rational evaluation
    also preserves sums whose individual binary products would underflow.
    Physical inputs remain analyst assumptions; exact arithmetic does not
    establish their validity or quantify their uncertainty.

    References
    ----------
    BIPM-SI: E=sum(P_i*t_i), average power=E/total_time.
    NASA-7009B: numerical domain and finite representation boundaries.
    """
    p, t, d, idle, warm, tw = map(Fraction, (firing, on_time, duty, standby, warmup, warmup_time))
    standby_time = t * (1 - d) / d
    elapsed = t + standby_time + tw
    energy = p * t + idle * standby_time + warm * tw
    peak = max(p if t else 0, idle if standby_time else 0, warm if tw else 0)
    rows = (("energy", energy, "J"), ("elapsed_time", elapsed, "s"),
            ("peak_bus_power", peak, "W"),
            ("average_bus_power", energy / elapsed if elapsed else 0, "W"))
    return binary_ledger(rows, "Power budget")
