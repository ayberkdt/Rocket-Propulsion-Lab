"""Direct-acting pressure-regulator characteristic with lockup, droop and hysteresis.

The steady force balance of a direct-acting, pressure-unbalanced regulator
[1]_ sets the regulated pressure as a linear function of inlet-valve stroke
and supply pressure.  Larger flow demand needs a larger stroke, so the
regulated pressure droops with flow; zero stroke defines the lockup pressure
at which all flow stops.  Coulomb friction on the moving assembly [1]_ makes
the stroke stick until the force imbalance exceeds the friction force, which
gives separate opening and closing branches and a memory state: the stroke.

Gas flow through the inlet valve uses the isentropic compressible-flow
relation, choked at Mach one [2]_, [3]_.  The resulting regulated pressure can
replace the fixed set pressure of the existing single-node regulated feed
model (`feed_configuration_with_regulated_pressure`).

References
----------
.. [1] NASA-SP8080-FORCE-BALANCE. NASA SP-8080, "Liquid Rocket Pressure
       Regulators, Relief Valves, Check Valves, Burst Disks, and Explosive
       Valves", 1973. https://ntrs.nasa.gov/citations/19740002611
.. [2] NASA-GRC-MASS-FLOW-EQ10. NASA Glenn Research Center, "Mass Flow
       Choking". https://www.grc.nasa.gov/www/k-12/BGP/mflchk.html
.. [3] NASA-GRC-ISENTROPIC. NASA Glenn Research Center, "Isentropic Flow".
       https://www.grc.nasa.gov/www/k-12/airplane/isentrop.html
"""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass, field, replace
from enum import Enum
from math import isfinite, sqrt

from rocket_propulsion.core.errors import DomainError
from rocket_propulsion.core.validation import coerce_fields, real

from .feed_system import IdealPressurantGas, RegulatedFeedConfiguration

REGULATOR_HYSTERESIS_SCHEMA = "rocket_propulsion_regulator_hysteresis_v1"
REGULATOR_HYSTERESIS_REFERENCE_IDS = (
    "NASA-SP8080-FORCE-BALANCE",
    "NASA-GRC-MASS-FLOW-EQ10",
    "NASA-GRC-ISENTROPIC",
)
#: Orifice discharge coefficient recommended by NASA SP-8080 (printed p. 76,
#: section 3.1.2.1) for preliminary inlet-valve sizing "unless a different
#: value is known to be appropriate".
SP8080_DEFAULT_DISCHARGE_COEFFICIENT = 0.6
_BISECTION_ITERATIONS = 200


def _positive(value: float, label: str) -> float:
    """Return ``value`` as float, raising `DomainError` unless finite and positive."""
    number = real(value, label)
    if number <= 0.0:
        raise DomainError(f"{label} must be finite and greater than zero.", field=label)
    return number


def _nonnegative(value: float, label: str) -> float:
    """Return ``value`` as float, raising `DomainError` unless finite and non-negative."""
    number = real(value, label)
    if number < 0.0:
        raise DomainError(f"{label} must be finite and non-negative.", field=label)
    return number


def compressible_valve_mass_flow_kg_s(
    *,
    gas: IdealPressurantGas,
    effective_area_m2: float,
    upstream_pressure_pa: float,
    upstream_temperature_k: float,
    downstream_pressure_pa: float,
) -> float:
    """Return the isentropic mass flow through a restriction, choked when sonic.

    Parameters
    ----------
    gas : IdealPressurantGas
        Calorically perfect gas.
    effective_area_m2 : float
        Non-negative effective flow area ``Cd A``, m^2.
    upstream_pressure_pa : float
        Positive upstream total pressure ``pt``, Pa.
    upstream_temperature_k : float
        Positive upstream total temperature ``Tt``, K.
    downstream_pressure_pa : float
        Positive downstream static pressure ``p``, Pa.

    Returns
    -------
    float
        Mass flow, kg/s; zero for zero area or non-forward pressure drop.

    Raises
    ------
    DomainError
        If an input is outside its domain.

    Notes
    -----
    The throat Mach number follows from ``p / pt = (1 + (gam-1)/2 M**2)**(-gam/(gam-1))``
    [2]_, limited to ``M <= 1``, and the flow is [1]_

    ``mdot = A pt / sqrt(Tt) sqrt(gam/R) M (1 + (gam-1)/2 M**2)**(-(gam+1)/(2(gam-1)))``.

    At or below the critical pressure ratio the choked value ``M = 1`` is
    used.

    References
    ----------
    .. [1] NASA-GRC-MASS-FLOW-EQ10. NASA Glenn Research Center, "Mass Flow
           Choking", eq. 10 and the choked form.
           https://www.grc.nasa.gov/www/k-12/BGP/mflchk.html
    .. [2] NASA-GRC-ISENTROPIC. NASA Glenn Research Center, "Isentropic Flow",
           eq. 6. https://www.grc.nasa.gov/www/k-12/airplane/isentrop.html
    """

    if not isinstance(gas, IdealPressurantGas):
        raise DomainError("Gas must be an IdealPressurantGas.", field="gas")
    area = _nonnegative(effective_area_m2, "Effective valve area")
    upstream = _positive(upstream_pressure_pa, "Upstream pressure")
    temperature = _positive(upstream_temperature_k, "Upstream temperature")
    downstream = _positive(downstream_pressure_pa, "Downstream pressure")
    if area == 0.0 or downstream >= upstream:
        return 0.0
    gamma = gas.gamma
    exponent = (gamma - 1.0) / gamma
    # Isentropic eq. 6 inverted: (pt/p)**((gam-1)/gam) = 1 + (gam-1)/2 M**2.
    ratio_power = (upstream / downstream) ** exponent
    mach_squared = min(1.0, 2.0 / (gamma - 1.0) * (ratio_power - 1.0))
    mach = sqrt(mach_squared)
    bracket = 1.0 + 0.5 * (gamma - 1.0) * mach_squared
    flow = (area * upstream / sqrt(temperature) * sqrt(gamma / gas.specific_gas_constant_j_kg_k)
            * mach * bracket ** (-(gamma + 1.0) / (2.0 * (gamma - 1.0))))
    return _nonnegative(flow, "Valve mass flow")


@dataclass(frozen=True, slots=True)
class DirectActingRegulator:
    """Geometry, springs and friction of a direct-acting regulator.

    Parameters
    ----------
    gas : IdealPressurantGas
        Regulated gas.
    actuator_area_m2 : float
        Positive sensor/actuator area ``A_A``, m^2.
    inlet_valve_area_m2 : float
        Positive unbalanced inlet-valve area ``A_i`` on which supply pressure
        acts, m^2; must be smaller than ``A_A``.
    reference_spring_load_n : float
        Non-negative reference-spring installed load ``L``, N.
    reference_spring_rate_n_m : float
        Non-negative reference-spring rate ``Y``, N/m.
    bias_spring_load_n : float
        Non-negative bias-spring installed load ``B``, N.
    bias_spring_rate_n_m : float
        Non-negative bias-spring rate ``Z``, N/m; ``Y + Z`` must be positive.
    reference_pressure_pa : float
        Non-negative pressure ``P_A`` on the reference side of the actuator
        (ambient or reference cavity), Pa.
    maximum_stroke_m : float
        Positive inlet-valve stroke limit, m.
    flow_area_per_stroke_m : float
        Positive geometric flow-area gradient ``dA/dX`` of the inlet valve,
        m^2/m.
    coulomb_friction_n : float, default 0.0
        Non-negative dry-friction force magnitude ``F_c`` on the moving
        assembly, N.
    maximum_flow_area_m2 : float or None, default None
        Optional positive geometric flow-area cap, m^2.
    discharge_coefficient : float, default 0.6
        Discharge coefficient in ``(0, 1]``; the default is the SP-8080
        preliminary-sizing value.
    model_id : str, default "direct-acting-regulator"
        Non-blank model identifier.
    schema : str, default REGULATOR_HYSTERESIS_SCHEMA
        Only the current schema is accepted.

    Raises
    ------
    DomainError
        If a value is outside its domain, ``A_i >= A_A``, ``Y + Z = 0`` or the
        identifier is blank.

    Notes
    -----
    SP-8080 writes the steady force balance on the inlet valve, actuator and
    springs of a direct-acting, pressure-unbalanced regulator as [1]_

    ``P_A A_A + L - Y X - P_R A_A + P_R A_i - B - Z X - P_S A_i = 0``

    so that the regulated-pressure error is
    ``dP_R = -dP_S A_i / (A_A - A_i) - dX (Y + Z) / (A_A - A_i)``.  The
    left-hand side is the net force in the opening direction.  With Coulomb
    damping, a force constant in magnitude and opposed to the motion [1]_,
    the valve opens only when that force exceeds ``F_c`` and closes only
    when it falls below ``-F_c``.  SP-8080 lists mechanical hysteresis and
    friction among the regulated-pressure error sources and recommends
    keeping the hysteresis and friction error below 25 percent of the
    allowable tolerance [1]_.  The effective flow area is
    ``Cd min(X dA/dX, A_max)``; SP-8080 recommends ``Cd = 0.6`` unless a
    better value is known [1]_.

    References
    ----------
    .. [1] NASA-SP8080-FORCE-BALANCE. NASA SP-8080, "Liquid Rocket Pressure
           Regulators, Relief Valves, Check Valves, Burst Disks, and
           Explosive Valves", 1973, printed pp. 11, 72-73 (eqs. 2-3), 76 and
           104 (Glossary). https://ntrs.nasa.gov/citations/19740002611
    """

    gas: IdealPressurantGas
    actuator_area_m2: float
    inlet_valve_area_m2: float
    reference_spring_load_n: float
    reference_spring_rate_n_m: float
    bias_spring_load_n: float
    bias_spring_rate_n_m: float
    reference_pressure_pa: float
    maximum_stroke_m: float
    flow_area_per_stroke_m: float
    coulomb_friction_n: float = 0.0
    maximum_flow_area_m2: float | None = None
    discharge_coefficient: float = SP8080_DEFAULT_DISCHARGE_COEFFICIENT
    model_id: str = "direct-acting-regulator"
    schema: str = REGULATOR_HYSTERESIS_SCHEMA

    def __post_init__(self) -> None:
        coerce_fields(self, "Regulator")
        if self.schema != REGULATOR_HYSTERESIS_SCHEMA:
            raise DomainError("Unsupported regulator schema.", field="schema")
        if not isinstance(self.gas, IdealPressurantGas):
            raise DomainError("Regulator gas must be an IdealPressurantGas.", field="gas")
        if not self.model_id.strip():
            raise DomainError("Regulator model identifier must not be empty.", field="model_id")
        for value, label in (
            (self.actuator_area_m2, "Actuator area"),
            (self.inlet_valve_area_m2, "Inlet valve area"),
            (self.maximum_stroke_m, "Maximum stroke"),
            (self.flow_area_per_stroke_m, "Flow area per stroke"),
        ):
            _positive(value, label)
        for value, label in (
            (self.reference_spring_load_n, "Reference spring load"),
            (self.reference_spring_rate_n_m, "Reference spring rate"),
            (self.bias_spring_load_n, "Bias spring load"),
            (self.bias_spring_rate_n_m, "Bias spring rate"),
            (self.reference_pressure_pa, "Reference pressure"),
            (self.coulomb_friction_n, "Coulomb friction"),
        ):
            _nonnegative(value, label)
        if self.maximum_flow_area_m2 is not None:
            _positive(self.maximum_flow_area_m2, "Maximum flow area")
        if self.inlet_valve_area_m2 >= self.actuator_area_m2:
            raise DomainError("Inlet valve area must be smaller than actuator area.",
                              field="inlet_valve_area_m2")
        if self.spring_rate_n_m <= 0.0:
            raise DomainError("Reference plus bias spring rate must be positive.",
                              field="reference_spring_rate_n_m")
        if not 0.0 < self.discharge_coefficient <= 1.0:
            raise DomainError("Discharge coefficient must be in (0, 1].",
                              field="discharge_coefficient")

    @property
    def net_area_m2(self) -> float:
        """Net sensing area ``A_A - A_i``, m^2 [1]_.

        References
        ----------
        .. [1] NASA-SP8080-FORCE-BALANCE. NASA SP-8080, printed p. 73, eq. 3.
               https://ntrs.nasa.gov/citations/19740002611
        """

        return self.actuator_area_m2 - self.inlet_valve_area_m2

    @property
    def spring_rate_n_m(self) -> float:
        """Combined spring rate ``Y + Z``, N/m [1]_.

        References
        ----------
        .. [1] NASA-SP8080-FORCE-BALANCE. NASA SP-8080, printed p. 73, eq. 3.
               https://ntrs.nasa.gov/citations/19740002611
        """

        return self.reference_spring_rate_n_m + self.bias_spring_rate_n_m

    @property
    def hysteresis_band_pa(self) -> float:
        """Closing-minus-opening branch pressure ``2 F_c / (A_A - A_i)``, Pa [1]_.

        References
        ----------
        .. [1] NASA-SP8080-FORCE-BALANCE. NASA SP-8080, printed pp. 72-73
               (eq. 2) and p. 104 (Coulomb damping).
               https://ntrs.nasa.gov/citations/19740002611
        """

        return 2.0 * self.coulomb_friction_n / self.net_area_m2

    def opening_force_n(self, *, stroke_m: float, regulated_pressure_pa: float,
                        supply_pressure_pa: float) -> float:
        """Return the net opening force, the left side of SP-8080 eq. (2) [1]_.

        Parameters
        ----------
        stroke_m : float
            Inlet-valve stroke ``X``, m.
        regulated_pressure_pa : float
            Regulated (outlet) pressure ``P_R``, Pa.
        supply_pressure_pa : float
            Supply (inlet) pressure ``P_S``, Pa.

        Returns
        -------
        float
            ``P_A A_A + L - Y X - P_R A_A + P_R A_i - B - Z X - P_S A_i``, N.

        References
        ----------
        .. [1] NASA-SP8080-FORCE-BALANCE. NASA SP-8080, printed p. 72, eq. 2.
               https://ntrs.nasa.gov/citations/19740002611
        """

        return (self.reference_pressure_pa * self.actuator_area_m2
                + self.reference_spring_load_n - self.bias_spring_load_n
                - self.spring_rate_n_m * stroke_m
                - regulated_pressure_pa * self.net_area_m2
                - supply_pressure_pa * self.inlet_valve_area_m2)

    def branch_pressure_pa(self, *, stroke_m: float, supply_pressure_pa: float,
                           branch: RegulatorBranch) -> float:
        """Return the regulated pressure on one branch at a given stroke.

        Parameters
        ----------
        stroke_m : float
            Inlet-valve stroke ``X``, m.
        supply_pressure_pa : float
            Supply pressure ``P_S``, Pa.
        branch : RegulatorBranch
            ``OPENING`` (force equals ``+F_c``), ``CLOSING`` (``-F_c``) or
            ``FRICTIONLESS`` (zero).

        Returns
        -------
        float
            ``P_R`` that makes the net opening force equal the branch
            friction, Pa.

        Notes
        -----
        Solving SP-8080 eq. (2) for ``P_R`` with friction ``s F_c`` gives [1]_

        ``P_R = (P_A A_A + L - B - (Y + Z) X - P_S A_i - s F_c) / (A_A - A_i)``

        with ``s = +1`` on the opening and ``-1`` on the closing branch.

        References
        ----------
        .. [1] NASA-SP8080-FORCE-BALANCE. NASA SP-8080, printed p. 72, eq. 2,
               and p. 104 (Coulomb damping).
               https://ntrs.nasa.gov/citations/19740002611
        """

        sign = {RegulatorBranch.OPENING: 1.0, RegulatorBranch.CLOSING: -1.0,
                RegulatorBranch.FRICTIONLESS: 0.0}[RegulatorBranch(branch)]
        force = self.opening_force_n(stroke_m=stroke_m, regulated_pressure_pa=0.0,
                                     supply_pressure_pa=supply_pressure_pa)
        return (force - sign * self.coulomb_friction_n) / self.net_area_m2

    def lockup_pressure_pa(self, supply_pressure_pa: float) -> float:
        """Return the pressure at which the closing valve seats and all flow stops.

        Parameters
        ----------
        supply_pressure_pa : float
            Positive supply pressure, Pa.

        Returns
        -------
        float
            Closing-branch pressure at zero stroke, Pa.

        Notes
        -----
        SP-8080 defines lockup as the stop of all gas flow when the regulated
        pressure exceeds the specified value, and requires the lockup
        pressure to lie a minimum above the maximum regulated pressure [1]_.
        In the force-balance model it is the closing-branch pressure at
        ``X = 0``.

        References
        ----------
        .. [1] NASA-SP8080-FORCE-BALANCE. NASA SP-8080, printed p. 9 and
               p. 104 (Glossary, "lockup").
               https://ntrs.nasa.gov/citations/19740002611
        """

        supply = _positive(supply_pressure_pa, "Supply pressure")
        return self.branch_pressure_pa(stroke_m=0.0, supply_pressure_pa=supply,
                                       branch=RegulatorBranch.CLOSING)

    def mass_flow_kg_s(self, *, stroke_m: float, supply_pressure_pa: float,
                       supply_temperature_k: float, regulated_pressure_pa: float) -> float:
        """Return inlet-valve mass flow at a given stroke and pressures [1]_, [2]_.

        Parameters
        ----------
        stroke_m : float
            Stroke in ``[0, X_max]``, m.
        supply_pressure_pa : float
            Positive supply pressure, Pa.
        supply_temperature_k : float
            Positive supply temperature, K.
        regulated_pressure_pa : float
            Positive regulated pressure, Pa.

        Returns
        -------
        float
            Mass flow, kg/s.

        References
        ----------
        .. [1] NASA-SP8080-FORCE-BALANCE. NASA SP-8080, printed p. 76,
               section 3.1.2.1. https://ntrs.nasa.gov/citations/19740002611
        .. [2] NASA-GRC-MASS-FLOW-EQ10. NASA Glenn Research Center, "Mass Flow
               Choking". https://www.grc.nasa.gov/www/k-12/BGP/mflchk.html
        """

        stroke = _nonnegative(stroke_m, "Regulator stroke")
        if stroke > self.maximum_stroke_m:
            raise DomainError("Regulator stroke exceeds its limit.", field="stroke_m")
        geometric = stroke * self.flow_area_per_stroke_m
        if self.maximum_flow_area_m2 is not None:
            geometric = min(geometric, self.maximum_flow_area_m2)
        return compressible_valve_mass_flow_kg_s(
            gas=self.gas,
            effective_area_m2=self.discharge_coefficient * geometric,
            upstream_pressure_pa=supply_pressure_pa,
            upstream_temperature_k=supply_temperature_k,
            downstream_pressure_pa=regulated_pressure_pa,
        )


class RegulatorBranch(str, Enum):
    """Branch of the regulator characteristic on which an operating point lies.

    ``OPENING`` and ``CLOSING`` points have the friction force at its limit
    in the respective direction; ``STUCK`` points lie inside the band with
    the stroke held by friction; ``FRICTIONLESS`` is the zero-friction
    curve; ``LOCKUP`` is the seated valve with zero flow [1]_.

    References
    ----------
    .. [1] NASA-SP8080-FORCE-BALANCE. NASA SP-8080, printed p. 72 (eq. 2) and
           p. 104 (Glossary). https://ntrs.nasa.gov/citations/19740002611
    """

    OPENING = "opening"
    CLOSING = "closing"
    STUCK = "stuck"
    FRICTIONLESS = "frictionless"
    LOCKUP = "lockup"


@dataclass(frozen=True, slots=True)
class RegulatorOperatingPoint:
    """Quasi-steady regulator state for one flow demand.

    Attributes
    ----------
    mass_flow_kg_s : float
        Delivered mass flow, kg/s.
    supply_pressure_pa : float
        Supply pressure, Pa.
    regulated_pressure_pa : float
        Outlet pressure delivered at that flow, Pa.
    stroke_m : float
        Inlet-valve stroke (the memory state), m.
    branch : RegulatorBranch
        Characteristic branch of the point.
    choked : bool
        Whether the inlet valve is sonic.
    saturated : bool
        Whether the stroke is at its limit, so the demand fixes the outlet
        pressure rather than the regulator.
    droop_pa : float
        Outlet pressure below the lockup pressure at the same supply, Pa.

    References
    ----------
    .. [1] NASA-SP8080-FORCE-BALANCE. NASA SP-8080, printed pp. 72-73.
           https://ntrs.nasa.gov/citations/19740002611
    """

    mass_flow_kg_s: float
    supply_pressure_pa: float
    regulated_pressure_pa: float
    stroke_m: float
    branch: RegulatorBranch
    choked: bool
    saturated: bool
    droop_pa: float
    reference_ids: tuple[str, ...] = field(default=REGULATOR_HYSTERESIS_REFERENCE_IDS)


def _is_choked(gas: IdealPressurantGas, supply: float, regulated: float) -> bool:
    """Whether ``regulated / supply`` is at or below the critical ratio."""
    gamma = gas.gamma
    critical = (2.0 / (gamma + 1.0)) ** (gamma / (gamma - 1.0))
    return regulated <= critical * supply


def _bisect(function, low: float, high: float) -> float:
    """Root of a monotone function with a sign change on ``[low, high]``."""
    f_low = function(low)
    for _ in range(_BISECTION_ITERATIONS):
        middle = 0.5 * (low + high)
        if middle in (low, high):
            break
        f_middle = function(middle)
        if (f_middle > 0.0) == (f_low > 0.0):
            low, f_low = middle, f_middle
        else:
            high = middle
    return 0.5 * (low + high)


def regulator_operating_point(
    regulator: DirectActingRegulator,
    *,
    mass_flow_kg_s: float,
    supply_pressure_pa: float,
    supply_temperature_k: float,
    previous_stroke_m: float = 0.0,
) -> RegulatorOperatingPoint:
    """Return the quasi-steady outlet pressure for a flow demand, with memory.

    Parameters
    ----------
    regulator : DirectActingRegulator
        Regulator description.
    mass_flow_kg_s : float
        Non-negative demanded mass flow, kg/s.
    supply_pressure_pa : float
        Positive supply pressure, Pa.
    supply_temperature_k : float
        Positive supply temperature, K.
    previous_stroke_m : float, default 0.0
        Stroke at the previous operating point, in ``[0, X_max]``, m.

    Returns
    -------
    RegulatorOperatingPoint
        Operating point; its ``stroke_m`` is the memory for the next call.

    Raises
    ------
    DomainError
        If an input is outside its domain, the lockup pressure is not below
        the supply pressure, or the demanded flow cannot pass the fully open
        valve at any positive outlet pressure.

    Notes
    -----
    Zero demand seats the valve at the lockup pressure (`lockup_pressure_pa`)
    [1]_.  For positive demand the previous stroke is held while the outlet
    pressure that passes the demand through that stroke lies inside the
    friction band; otherwise the stroke moves along the opening branch (a
    pressure below the band) or the closing branch (above it) until the
    valve passes the demand.  Along both branches the flow is monotone in
    stroke because the outlet pressure falls as the stroke rises [1]_, so the
    root is found by bisection.  The flow law is the isentropic restriction
    relation [2]_.

    References
    ----------
    .. [1] NASA-SP8080-FORCE-BALANCE. NASA SP-8080, printed pp. 9, 72-73 and
           104. https://ntrs.nasa.gov/citations/19740002611
    .. [2] NASA-GRC-MASS-FLOW-EQ10. NASA Glenn Research Center, "Mass Flow
           Choking". https://www.grc.nasa.gov/www/k-12/BGP/mflchk.html
    """

    if not isinstance(regulator, DirectActingRegulator):
        raise DomainError("Regulator must be a DirectActingRegulator.", field="regulator")
    demand = _nonnegative(mass_flow_kg_s, "Demanded mass flow")
    supply = _positive(supply_pressure_pa, "Supply pressure")
    temperature = _positive(supply_temperature_k, "Supply temperature")
    previous = _nonnegative(previous_stroke_m, "Previous stroke")
    if previous > regulator.maximum_stroke_m:
        raise DomainError("Previous stroke exceeds the stroke limit.", field="previous_stroke_m")
    lockup = regulator.lockup_pressure_pa(supply)
    if not 0.0 < lockup < supply:
        raise DomainError("Lockup pressure must be positive and below the supply pressure.")
    frictionless = regulator.coulomb_friction_n == 0.0

    def flow_at(stroke: float, pressure: float) -> float:
        if pressure <= 0.0:
            pressure = supply * 1e-12
        return regulator.mass_flow_kg_s(stroke_m=stroke, supply_pressure_pa=supply,
                                        supply_temperature_k=temperature,
                                        regulated_pressure_pa=min(pressure, supply))

    def point(pressure: float, stroke: float, branch: RegulatorBranch,
              saturated: bool = False) -> RegulatorOperatingPoint:
        if not isfinite(pressure) or pressure <= 0.0:
            raise DomainError("Demand cannot be delivered at a positive outlet pressure.")
        return RegulatorOperatingPoint(
            mass_flow_kg_s=demand, supply_pressure_pa=supply,
            regulated_pressure_pa=pressure, stroke_m=stroke, branch=branch,
            choked=_is_choked(regulator.gas, supply, pressure) and demand > 0.0,
            saturated=saturated, droop_pa=lockup - pressure)

    if demand == 0.0:
        return point(lockup, 0.0, RegulatorBranch.LOCKUP)

    def pressure_through(stroke: float) -> float:
        """Outlet pressure that passes the demand through a fixed stroke."""
        if flow_at(stroke, 0.0) < demand:
            return float("-inf")
        return _bisect(lambda p: flow_at(stroke, p) - demand, 0.0, supply)

    if previous > 0.0 and not frictionless:
        stuck = pressure_through(previous)
        opening_limit = regulator.branch_pressure_pa(
            stroke_m=previous, supply_pressure_pa=supply, branch=RegulatorBranch.OPENING)
        closing_limit = regulator.branch_pressure_pa(
            stroke_m=previous, supply_pressure_pa=supply, branch=RegulatorBranch.CLOSING)
        if opening_limit <= stuck <= closing_limit:
            return point(stuck, previous, RegulatorBranch.STUCK)
        moving_open = stuck < opening_limit
    else:
        moving_open = True
    if frictionless:
        branch, low, high = RegulatorBranch.FRICTIONLESS, 0.0, regulator.maximum_stroke_m
    elif moving_open:
        branch, low, high = RegulatorBranch.OPENING, previous, regulator.maximum_stroke_m
    else:
        branch, low, high = RegulatorBranch.CLOSING, 0.0, previous

    def branch_pressure(stroke: float) -> float:
        return regulator.branch_pressure_pa(stroke_m=stroke, supply_pressure_pa=supply,
                                            branch=branch)

    def excess(stroke: float) -> float:
        return flow_at(stroke, branch_pressure(stroke)) - demand

    if excess(high) < 0.0:
        if branch is RegulatorBranch.CLOSING:  # pragma: no cover - closing keeps flow above demand
            raise DomainError("Closing-branch solution not bracketed.")
        stroke = regulator.maximum_stroke_m
        return point(pressure_through(stroke), stroke, branch, saturated=True)
    if excess(low) >= 0.0:
        stroke = low
    else:
        stroke = _bisect(excess, low, high)
    return point(branch_pressure(stroke), stroke, branch)


def sweep_regulator(
    regulator: DirectActingRegulator,
    *,
    mass_flows_kg_s: Sequence[float],
    supply_pressures_pa: Sequence[float] | float,
    supply_temperature_k: float,
    initial_stroke_m: float = 0.0,
) -> tuple[RegulatorOperatingPoint, ...]:
    """Evaluate a sequence of demands, carrying the stroke memory forward.

    Parameters
    ----------
    regulator : DirectActingRegulator
        Regulator description.
    mass_flows_kg_s : sequence of float
        Non-empty sequence of non-negative demands, kg/s.
    supply_pressures_pa : sequence of float or float
        Supply pressure per demand, or one value for all, Pa.
    supply_temperature_k : float
        Positive supply temperature, K.
    initial_stroke_m : float, default 0.0
        Stroke before the first demand, m.

    Returns
    -------
    tuple of RegulatorOperatingPoint
        One point per demand; an increasing then decreasing sweep traces the
        opening and closing branches of the hysteresis loop [1]_.

    Raises
    ------
    DomainError
        If the sequences are empty or of different lengths, or a point
        fails.

    References
    ----------
    .. [1] NASA-SP8080-FORCE-BALANCE. NASA SP-8080, printed pp. 11 and 72-73.
           https://ntrs.nasa.gov/citations/19740002611
    """

    flows = tuple(mass_flows_kg_s)
    if not flows:
        raise DomainError("At least one demand is required.", field="mass_flows_kg_s")
    if isinstance(supply_pressures_pa, (int, float)) and not isinstance(supply_pressures_pa, bool):
        supplies = (supply_pressures_pa,) * len(flows)
    else:
        supplies = tuple(supply_pressures_pa)
    if len(supplies) != len(flows):
        raise DomainError("Supply pressures and demands must have equal length.",
                          field="supply_pressures_pa")
    stroke = initial_stroke_m
    points = []
    for flow, supply in zip(flows, supplies, strict=True):
        result = regulator_operating_point(
            regulator, mass_flow_kg_s=flow, supply_pressure_pa=supply,
            supply_temperature_k=supply_temperature_k, previous_stroke_m=stroke)
        points.append(result)
        stroke = result.stroke_m
    return tuple(points)


def feed_configuration_with_regulated_pressure(
    configuration: RegulatedFeedConfiguration, operating_point: RegulatorOperatingPoint,
) -> RegulatedFeedConfiguration:
    """Return a feed configuration whose set pressure is the regulator outlet pressure.

    Parameters
    ----------
    configuration : RegulatedFeedConfiguration
        Existing single-node regulated feed configuration.
    operating_point : RegulatorOperatingPoint
        Regulator operating point, for example at the current flow and
        bottle pressure.

    Returns
    -------
    RegulatedFeedConfiguration
        Copy with ``regulator_set_pressure_pa`` replaced; all validation of
        the feed configuration is repeated.

    Raises
    ------
    DomainError
        If an argument has the wrong type or the new set pressure violates
        the feed configuration (for example the overpressure limit).

    Notes
    -----
    The existing feed model regulates toward a fixed set pressure [1]_.  This
    helper lets a caller update that set point between integration steps
    from the droop and hysteresis characteristic of SP-8080 [2]_; it does not
    change the feed model's own regulator dynamics.

    References
    ----------
    .. [1] NASA SP-8112. NASA, "Pressurization Systems for Liquid Rockets",
           1975. https://ntrs.nasa.gov/citations/19760015212
    .. [2] NASA-SP8080-FORCE-BALANCE. NASA SP-8080, printed pp. 72-73.
           https://ntrs.nasa.gov/citations/19740002611
    """

    if not isinstance(configuration, RegulatedFeedConfiguration):
        raise DomainError("Configuration must be a RegulatedFeedConfiguration.")
    if not isinstance(operating_point, RegulatorOperatingPoint):
        raise DomainError("Operating point must be a RegulatorOperatingPoint.")
    return replace(configuration, regulator_set_pressure_pa=operating_point.regulated_pressure_pa)


def regulated_pressure_sensitivities(regulator: DirectActingRegulator) -> tuple[float, float]:
    """Return ``(dP_R/dP_S, dP_R/dX)`` at fixed branch from SP-8080 eq. (3).

    Parameters
    ----------
    regulator : DirectActingRegulator
        Regulator description.

    Returns
    -------
    tuple of float
        ``(-A_i / (A_A - A_i), -(Y + Z) / (A_A - A_i))``: the supply-pressure
        error per pascal and the position (flow-demand) error per metre of
        stroke, in Pa/Pa and Pa/m [1]_.

    References
    ----------
    .. [1] NASA-SP8080-FORCE-BALANCE. NASA SP-8080, printed p. 73, eq. 3.
           https://ntrs.nasa.gov/citations/19740002611
    """

    if not isinstance(regulator, DirectActingRegulator):
        raise DomainError("Regulator must be a DirectActingRegulator.", field="regulator")
    net = regulator.net_area_m2
    return (-regulator.inlet_valve_area_m2 / net, -regulator.spring_rate_n_m / net)

