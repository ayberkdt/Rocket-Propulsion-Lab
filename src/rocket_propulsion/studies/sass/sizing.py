"""Mission-level propellant, impulse and resource budgets with complete evidence.

References
----------
.. [1] NASA-ROCKET. NASA Glenn Research Center, "Ideal Rocket Equation"; inverse
       ideal rocket equation.
.. [2] NASA-ISP. NASA Glenn Research Center, "Specific Impulse"; thrust/flow/impulse
       identities.
.. [3] NASA-SE-6. NASA, "Systems Engineering Handbook: Crosscutting Technical
       Management"; analyst-defined margins, constraints and subsystem interfaces.
"""

from __future__ import annotations

from dataclasses import dataclass
from math import fsum, isfinite

from rocket_propulsion.core.constants import STANDARD_GRAVITY_M_S2
from rocket_propulsion.core.errors import DomainError
from rocket_propulsion.core.validation import integer, real, records, text

from .evidence import AnalysisResult, Citation, CitedScalar, derived, quantity
from .mass_numerics import solve_mass_budget
from .power_numerics import solve_power_budget


@dataclass(frozen=True, slots=True)
class DeltaVItem:
    """One maneuver allocation with an explicit margin and occurrence count.

    Attributes
    ----------
    name : str
        Unique allocation name.
    nominal_delta_v : CitedScalar
        Nominal delta-v per occurrence, m/s; non-negative.
    margin_fraction : CitedScalar
        Analyst-chosen margin as a fraction of the nominal value; non-negative.
    count : int, default 1
        Number of occurrences; positive.

    Raises
    ------
    DomainError
        If a field is missing, has the wrong unit, or is out of domain.

    Notes
    -----
    Allocations and margins are explicit project inputs [1]_; no default margin is
    supplied.

    References
    ----------
    .. [1] NASA-SE-6. NASA, "Systems Engineering Handbook: Crosscutting Technical
           Management"; sections 6.2 and 6.8, allocated requirements and assumptions.
    """

    name: str
    nominal_delta_v: CitedScalar
    margin_fraction: CitedScalar
    count: int = 1

    def __post_init__(self) -> None:
        """Validate allocations without supplying an unsupported default margin.

        References
        ----------
        .. [1] NASA-SE-6. NASA, "Systems Engineering Handbook: Crosscutting Technical
               Management"; sections 6.2 and 6.8, controlled requirements and margins.
        """
        text(self.name, "maneuver name")
        quantity(self.nominal_delta_v, "m/s", minimum=0)
        quantity(self.margin_fraction, "1", minimum=0)
        count = integer(self.count, "maneuver count")
        real(count, "finite maneuver count")
        if count < 1:
            raise DomainError("Maneuver count must be positive.")
        object.__setattr__(self, "count", count)


def delta_v_budget(items: tuple[DeltaVItem, ...]) -> AnalysisResult:
    """Sum a delta-v allocation without applying any margin twice.

    Parameters
    ----------
    items : tuple of DeltaVItem
        Non-empty, uniquely named allocations.

    Returns
    -------
    AnalysisResult
        ``nominal_delta_v`` = sum(count * nominal), ``allocated_delta_v`` =
        sum(count * nominal * (1 + margin)) and one margined entry per item, m/s.

    Raises
    ------
    DomainError
        If the collection is empty, of the wrong type, or has duplicate names.

    Notes
    -----
    The margin is a project choice [1]_, not a fixed percentage prescribed by NASA.
    Allocations are inputs, not a computed station-keeping or collision-avoidance
    history.

    References
    ----------
    .. [1] NASA-SE-6. NASA, "Systems Engineering Handbook: Crosscutting Technical
           Management"; sections 6.2 and 6.8, allocated budgets and explicit
           assumptions.
    """
    entries = records(items, DeltaVItem, "delta-v items")
    if not entries or len({x.name for x in entries}) != len(entries):
        raise DomainError("Delta-v budget needs uniquely named entries.")
    inputs = tuple(v for x in entries for v in (x.nominal_delta_v, x.margin_fraction))
    methods = (Citation("NASA-SE-6", "6.2/6.8: explicit allocated contingency arithmetic"),)
    nominal = fsum(x.count * x.nominal_delta_v.value for x in entries)
    margined = fsum(x.count * x.nominal_delta_v.value * (1 + x.margin_fraction.value)
                   for x in entries)
    values = [("nominal_delta_v", nominal), ("allocated_delta_v", margined)]
    values.extend((x.name, x.count * x.nominal_delta_v.value * (1 + x.margin_fraction.value))
                  for x in entries)
    return AnalysisResult("explicit_delta_v_allocation", tuple(
        (name, derived(value, "m/s", inputs, methods)) for name, value in values
    ), ("Allocations are inputs, not a computed station-keeping or collision-avoidance history.",))


@dataclass(frozen=True, slots=True)
class Maneuver:
    """One executed ideal burn at a constant system operating point.

    Attributes
    ----------
    name : str
        Unique maneuver name.
    delta_v : CitedScalar
        Useful (axial) delta-v of the burn, m/s; non-negative.
    scalar_thrust : CitedScalar
        Total scalar thrust of the simultaneously firing thrusters, N; positive.
    system_isp : CitedScalar
        Total-system specific impulse at the same operating point, s; positive.
    axial_efficiency : CitedScalar
        Fraction of scalar thrust acting along the velocity change, in (0, 1].
        It reduces useful force only; tank mass flow is unchanged.
    duty_fraction : CitedScalar
        On-time divided by elapsed time, in (0, 1].

    Raises
    ------
    DomainError
        If a field is out of domain or the derived exhaust velocity or mass flow is
        not a positive finite number.

    Notes
    -----
    Mass flow mdot = F / (g0 * Isp) [1]_; the useful acceleration is
    eta_axial * F / m [2]_.

    References
    ----------
    .. [1] NASA-ISP. NASA Glenn Research Center, "Specific Impulse"; system Isp from
           total exhaust flow.
    .. [2] MIT-522-L5. MIT OCW, "16.522 Session 5: Review of Classical Astrodynamics";
           equation (19), useful force/vehicle mass.
    """

    name: str
    delta_v: CitedScalar
    scalar_thrust: CitedScalar
    system_isp: CitedScalar
    axial_efficiency: CitedScalar
    duty_fraction: CitedScalar

    def __post_init__(self) -> None:
        """Reject zero flow, invalid efficiency or duty fraction.

        References
        ----------
        .. [1] NASA-ISP. NASA Glenn Research Center, "Specific Impulse"; positive thrust and
               Isp.
        .. [2] NASA-SE-6. NASA, "Systems Engineering Handbook: Crosscutting Technical
               Management"; explicit operational inputs.
        """
        text(self.name, "maneuver name")
        quantity(self.delta_v, "m/s", minimum=0)
        quantity(self.scalar_thrust, "N", strictly_positive=True)
        quantity(self.system_isp, "s", strictly_positive=True)
        for item in (self.axial_efficiency, self.duty_fraction):
            if quantity(item, "1", strictly_positive=True) > 1:
                raise DomainError("Axial efficiency and duty fraction must not exceed one.")
        exhaust = STANDARD_GRAVITY_M_S2 * self.system_isp.value
        effective = exhaust * self.axial_efficiency.value
        if not isfinite(exhaust) or not isfinite(effective) or effective <= 0:
            raise DomainError("Effective exhaust velocity is outside the finite numerical domain.")
        flow = self.scalar_thrust.value / exhaust
        if not isfinite(flow) or flow <= 0:
            raise DomainError("Mass flow is outside the positive finite numerical domain.")


@dataclass(frozen=True, slots=True)
class MissionSizing:
    """Mission mass summary and the ordered per-burn mass history.

    Attributes
    ----------
    summary : AnalysisResult
        Dry, wet, loaded, consumed, retained and final mass, total delta-v,
        impulse, on-time, elapsed time and closure residuals.
    burns : tuple of (str, AnalysisResult)
        Per-burn ledgers in execution order.
    numerical_evidence : tuple of (str, str)
        Decimal working-precision diagnostics of the mass closure.

    Notes
    -----
    Dry mass includes all dry hardware; burns follow the sequential rocket-equation
    mass ratios [1]_. Results and limitations are preserved for the trade report
    [2]_.

    References
    ----------
    .. [1] NASA-ROCKET. NASA Glenn Research Center, "Ideal Rocket Equation"; sequential
           rocket-equation mass ratios.
    .. [2] NASA-SE-APP. NASA, "Systems Engineering Handbook: Appendix"; preserve
           computational results and limitations in trade reports.
    """

    summary: AnalysisResult
    burns: tuple[tuple[str, AnalysisResult], ...]
    numerical_evidence: tuple[tuple[str, str], ...] = ()

    def to_dict(self) -> dict[str, object]:
        """Preserve evidence in the mission summary and every executed maneuver.

        References
        ----------
        .. [1] NASA-SE-APP. NASA, "Systems Engineering Handbook: Appendix"; Trade Study
               Report, results and all data sources.
        """
        return {"summary": self.summary.to_dict(), "numerics": dict(self.numerical_evidence),
                "burns": [{"name": name, **value.to_dict()} for name, value in self.burns]}


def size_mission(dry_mass: CitedScalar, maneuvers: tuple[Maneuver, ...], *,
                 residual_fraction: CitedScalar, fixed_reserve: CitedScalar) -> MissionSizing:
    """Close the wet mass when the residual is a fraction of LOADED propellant.

    Parameters
    ----------
    dry_mass : CitedScalar
        Spacecraft dry mass including thrusters, tanks, pressurant hardware and all
        other dry hardware, kg; positive.
    maneuvers : tuple of Maneuver
        Non-empty, uniquely named burns in execution order.
    residual_fraction : CitedScalar
        Unusable residual as a fraction of loaded propellant, in [0, 1).
    fixed_reserve : CitedScalar
        Additional retained propellant reserve, kg; non-negative.

    Returns
    -------
    MissionSizing
        Mission summary, per-burn ledgers and numerical closure evidence.

    Raises
    ------
    DomainError
        If an input is out of domain, no finite wet mass closes the residual and
        delta-v, an output is not representable, or the closure does not resolve
        within the decimal precision budget.

    Notes
    -----
    For S = sum(dv / (g0 * Isp * eta_axial)) and A = exp(S), the loaded propellant is
    [dry * (A - 1) + A * fixed_reserve] / (1 - A * residual_fraction) and the final
    retained propellant is residual_fraction * loaded + fixed_reserve. This algebraic
    rearrangement of the rocket equation [1]_ is derived in this module; neither the
    residual fraction nor the fixed reserve is prescribed by the source. Flow and
    on-time follow mdot = F / (g0 * Isp) and t_on = m_consumed / mdot [2]_. Evaluation
    uses the equivalent inverse-ratio form with isolated Decimal precision, separate
    total and propellant closure and agreement at two precisions; g0 is the
    registered conventional standard gravity [3]_. Margins and allocations are
    project assumptions [4]_.

    References
    ----------
    .. [1] NASA-ROCKET. NASA Glenn Research Center, "Ideal Rocket Equation"; inverse
           mass ratio, applied to each ordered burn.
    .. [2] NASA-ISP. NASA Glenn Research Center, "Specific Impulse"; mdot=F/(g0 Isp),
           on-time=consumed/mdot.
    .. [3] BIPM-SI. BIPM, "SI Brochure, 9th edition"; conventional standard gravity from
           core.constants.
    .. [4] NASA-SE-6. NASA, "Systems Engineering Handbook: Crosscutting Technical
           Management"; section 6.8, project assumptions and resource allocations.
    """
    dry = quantity(dry_mass, "kg", strictly_positive=True)
    residual = quantity(residual_fraction, "1", minimum=0)
    reserve = quantity(fixed_reserve, "kg", minimum=0)
    burns = records(maneuvers, Maneuver, "mission maneuvers")
    if residual >= 1 or not burns or len({x.name for x in burns}) != len(burns):
        raise DomainError("Use residual fraction <1 and a non-empty unique maneuver sequence.")
    inputs = (dry_mass, residual_fraction, fixed_reserve, *(
        item for burn in burns for item in
        (burn.delta_v, burn.scalar_thrust, burn.system_isp,
         burn.axial_efficiency, burn.duty_fraction)
    ))
    methods = (Citation("NASA-ROCKET", "Inverse rocket equation; residual closure derived here"),
               Citation("NASA-ISP", "Isp=F/(mdot*g0); impulse=exhausted_mass*g0*Isp"),
               Citation("BIPM-SI", "Conventional g0 in core.constants"),
               Citation("NASA-7009B", "4.2.3: documented numerical verification and calculation limits"))
    summary_values, burn_values, numerical_evidence = solve_mass_budget(dry, residual, reserve, tuple(
        (b.delta_v.value, b.scalar_thrust.value, b.system_isp.value,
         b.axial_efficiency.value, b.duty_fraction.value) for b in burns))
    results = []
    for burn, values in zip(burns, burn_values, strict=True):
        results.append((burn.name, AnalysisResult("constant_system_isp_burn", tuple(
            (name, derived(value, unit, inputs, methods)) for name, value, unit in values
        ), ("Ideal useful delta-v; no gravity/drag/steering history or duty-cycle trajectory.",))))
    summary = AnalysisResult("ordered_mission_residual_closed_mass_budget", tuple(
        (name, derived(value, unit, inputs, methods)) for name, value, unit in summary_values
    ), ("Dry mass must include thrusters, tanks, pressurant and all other dry spacecraft hardware.",
        "Isp is total-system Isp at the supplied operating point; no catalogue-envelope inference.",
        "Fixed reserve is retained in addition to the loaded-propellant residual fraction."))
    return MissionSizing(summary, tuple(results), numerical_evidence)


def bipropellant_split(total_propellant: CitedScalar, oxidizer_fuel_ratio: CitedScalar
                       ) -> AnalysisResult:
    """Split a total propellant mass into oxidizer and fuel at a constant mixture ratio.

    Parameters
    ----------
    total_propellant : CitedScalar
        Total propellant mass, kg; non-negative.
    oxidizer_fuel_ratio : CitedScalar
        Mixture ratio O/F = m_ox / m_fuel; positive.

    Returns
    -------
    AnalysisResult
        ``fuel_mass`` = m / (1 + O/F) and ``oxidizer_mass`` = m * O/F / (1 + O/F), kg.

    Raises
    ------
    DomainError
        If a quantity has the wrong unit or is out of domain.

    Notes
    -----
    The mixture ratio is oxidizer over fuel [1]_ and the total propellant mass is the
    sum of both components [2]_. A constant total-drain mixture ratio is assumed;
    component residuals require separate allocations.

    References
    ----------
    .. [1] ARIANE-10N. ArianeGroup, "10 N Bipropellant Thruster"; mixture ratio is
           oxidizer/fuel.
    .. [2] NASA-ISP. NASA Glenn Research Center, "Specific Impulse"; total propellant
           mass.
    """
    total = quantity(total_propellant, "kg", minimum=0)
    ratio = quantity(oxidizer_fuel_ratio, "1", strictly_positive=True)
    fuel = total / (1 + ratio)
    methods = (Citation("ARIANE-10N", "Definition of nominal oxidizer/fuel mixture ratio"),)
    inputs = (total_propellant, oxidizer_fuel_ratio)
    return AnalysisResult("constant_mixture_mass_split", (
        ("fuel_mass", derived(fuel, "kg", inputs, methods)),
        ("oxidizer_mass", derived(total * (ratio / (1 + ratio)), "kg", inputs, methods)),
    ), ("Constant total-drain mixture ratio; component residuals require separate allocations.",))


def power_resources(firing_power: CitedScalar, on_time: CitedScalar, duty_fraction: CitedScalar,
                    standby_power: CitedScalar, warmup_power: CitedScalar,
                    warmup_time: CitedScalar) -> AnalysisResult:
    """Compute propulsion bus energy over disjoint warm-up, firing and standby phases.

    Parameters
    ----------
    firing_power : CitedScalar
        Complete bus load while firing, W; non-negative.
    on_time : CitedScalar
        Total on-time, s; non-negative.
    duty_fraction : CitedScalar
        On-time divided by firing-period elapsed time, in (0, 1].
    standby_power : CitedScalar
        Bus load during the off part of the duty cycle, W; non-negative.
    warmup_power : CitedScalar
        Bus load during warm-up, W; non-negative.
    warmup_time : CitedScalar
        Warm-up duration preceding firing, s; non-negative.

    Returns
    -------
    AnalysisResult
        ``energy`` (J), ``elapsed_time`` (s), ``peak_bus_power`` (W) and
        ``average_bus_power`` (W).

    Raises
    ------
    DomainError
        If a quantity is out of domain or an output is not representable.

    Notes
    -----
    E = P_fire * t_on + P_idle * (t_on / duty - t_on) + P_warm * t_warm, summing
    energy over disjoint intervals with W = J/s [1]_. Each input power is the complete
    bus load of its own phase and must not already include another phase. Solar
    array, battery, PPU efficiency and thermal effects are outside this propulsion/EPS
    interface budget [2]_.

    References
    ----------
    .. [1] BIPM-SI. BIPM, "SI Brochure, 9th edition"; watt=joule/second, energy
           integrated over disjoint intervals.
    .. [2] NASA-SE-6. NASA, "Systems Engineering Handbook: Crosscutting Technical
           Management"; section 6.3, propulsion/EPS resource interfaces.
    """
    p = quantity(firing_power, "W", minimum=0)
    t = quantity(on_time, "s", minimum=0)
    duty = quantity(duty_fraction, "1", strictly_positive=True)
    idle = quantity(standby_power, "W", minimum=0)
    warm = quantity(warmup_power, "W", minimum=0)
    tw = quantity(warmup_time, "s", minimum=0)
    if duty > 1:
        raise DomainError("Duty fraction must not exceed one.")
    values = solve_power_budget(p, t, duty, idle, warm, tw)
    inputs = (firing_power, on_time, duty_fraction, standby_power, warmup_power, warmup_time)
    methods = (Citation("BIPM-SI", "W=J/s; sum of energy over disjoint operating intervals"),
               Citation("NASA-SE-6", "6.3: explicit propulsion-to-EPS operational interfaces"))
    return AnalysisResult("disjoint_phase_bus_energy_budget", tuple(
        (name, derived(value, unit, inputs, methods)) for name, value, unit in values
    ), (("Phase ledger uses exact rational input arithmetic with final binary conversion; "
        "physical uncertainty is not quantified."),
        "No solar-array, battery, PPU efficiency or spacecraft thermal-model inference."))
