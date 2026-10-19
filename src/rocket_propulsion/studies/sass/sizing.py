"""Mission-level propellant, impulse and resource budgets with complete evidence.

References
----------
NASA-ROCKET: inverse ideal rocket equation.
NASA-ISP: thrust/flow/impulse identities.
NASA-SE-6: analyst-defined margins, constraints and subsystem interfaces.
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
    """A maneuver allocation, explicit margin and occurrence count.

    References
    ----------
    NASA-SE-6: sections 6.2 and 6.8, allocated requirements and assumptions.
    """

    name: str
    nominal_delta_v: CitedScalar
    margin_fraction: CitedScalar
    count: int = 1

    def __post_init__(self) -> None:
        """Validate allocations without supplying an unsupported default margin.

        References
        ----------
        NASA-SE-6: sections 6.2 and 6.8, controlled requirements and margins.
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
    """Sum count*nominal and count*nominal*(1+margin); do not add margin twice.

    Margin is a project choice, not a fixed percentage prescribed by NASA.

    References
    ----------
    NASA-SE-6: sections 6.2 and 6.8, allocated budgets and explicit assumptions.
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
    """One executed ideal burn with explicit system performance and duty fraction.

    ``scalar_thrust`` and ``system_isp`` refer to the same operating point.
    ``axial_efficiency`` reduces useful force only; tank flow remains unchanged.

    References
    ----------
    NASA-ISP: system Isp from total exhaust flow.
    MIT-522-L5: equation (19), useful force/vehicle mass.
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
        NASA-ISP: positive thrust and Isp; NASA-SE-6: explicit operational inputs.
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
    """Summary and ordered mass history; dry mass includes all dry hardware.

    References
    ----------
    NASA-ROCKET: sequential rocket-equation mass ratios.
    NASA-SE-APP: preserve computational results and limitations in trade reports.
    """

    summary: AnalysisResult
    burns: tuple[tuple[str, AnalysisResult], ...]
    numerical_evidence: tuple[tuple[str, str], ...] = ()

    def to_dict(self) -> dict[str, object]:
        """Preserve evidence in the mission summary and every executed maneuver.

        References
        ----------
        NASA-SE-APP: Trade Study Report, results and all data sources.
        """
        return {"summary": self.summary.to_dict(), "numerics": dict(self.numerical_evidence),
                "burns": [{"name": name, **value.to_dict()} for name, value in self.burns]}


def size_mission(dry_mass: CitedScalar, maneuvers: tuple[Maneuver, ...], *,
                 residual_fraction: CitedScalar, fixed_reserve: CitedScalar) -> MissionSizing:
    """Close wet mass when residual is a fraction of LOADED propellant.

    For S=sum(dv/(g0 Isp axial_efficiency)), A=exp(S), loaded propellant is
    [dry*(A-1) + A*fixed_reserve] / (1-A*residual_fraction). Final retained
    propellant is residual_fraction*loaded + fixed_reserve. This algebraic
    rearrangement is derived here from the cited rocket equation; neither the
    residual fraction nor the fixed reserve is prescribed by the source.
    Evaluation uses the equivalent inverse-ratio form with isolated Decimal
    precision, separate total/propellant closure and two-precision agreement.

    References
    ----------
    NASA-ROCKET: inverse mass ratio, applied to each ordered burn.
    NASA-ISP: mdot=F/(g0 Isp), on-time=consumed/mdot.
    NASA-SE-6: section 6.8, project assumptions and resource allocations.
    BIPM-SI: conventional standard gravity from core.constants.
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
    """Split total mass using m_ox/m_fuel=O/F and exact component mass closure.

    References
    ----------
    ARIANE-10N: mixture ratio is oxidizer/fuel; NASA-ISP: total propellant mass.
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
    """Compute E=P_fire*t_on+P_idle*(t_on/duty-t_on)+P_warm*t_warm.

    Warmup precedes firing. Input powers are complete bus loads in their own
    phases: they must not already include another supplied phase's power.

    References
    ----------
    BIPM-SI: watt=joule/second, energy integrated over disjoint intervals.
    NASA-SE-6: section 6.3, propulsion/EPS resource interfaces.
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
    ), ("Phase ledger uses exact rational input arithmetic with final binary conversion; "
        "physical uncertainty is not quantified.",
        "No solar-array, battery, PPU efficiency or spacecraft thermal-model inference."))
