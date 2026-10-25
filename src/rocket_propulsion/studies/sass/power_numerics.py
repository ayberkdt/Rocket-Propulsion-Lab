"""Exact rational phase ledger with a single final binary conversion.

References
----------
.. [1] BIPM-SI. BIPM, "SI Brochure, 9th edition"; disjoint phase energy and
       watt=joule/second.
.. [2] NASA-7009B. NASA, "Standard for Models and Simulations, NASA-STD-7009B";
       numerical implementation limits, distinct from physical uncertainty.
"""

from fractions import Fraction

from .rational_numerics import binary_ledger


def solve_power_budget(firing: float, on_time: float, duty: float,
                       standby: float, warmup: float, warmup_time: float) -> tuple:
    """Evaluate the phase energy ledger exactly before converting the final channels.

    Parameters
    ----------
    firing, standby, warmup : float
        Phase bus powers, W.
    on_time, warmup_time : float
        Firing on-time and warm-up duration, s.
    duty : float
        Duty fraction in (0, 1].

    Returns
    -------
    tuple of (str, float, str)
        ``energy``, ``elapsed_time``, ``peak_bus_power`` and ``average_bus_power``.

    Raises
    ------
    DomainError
        If an output is not representable.

    Notes
    -----
    Standby duration is t_on * (1 - duty) / duty; E = sum(P_i * t_i) and average power
    is E / elapsed time [1]_. Subtracting two rounded elapsed times could erase or
    double the standby phase near continuous duty, and exact rational evaluation also
    keeps sums whose individual binary products would underflow [2]_. Exact
    arithmetic does not validate the inputs or quantify their uncertainty.

    References
    ----------
    .. [1] BIPM-SI. BIPM, "SI Brochure, 9th edition"; E=sum(P_i*t_i), average
           power=E/total_time.
    .. [2] NASA-7009B. NASA, "Standard for Models and Simulations, NASA-STD-7009B";
           numerical domain and finite representation boundaries.
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
