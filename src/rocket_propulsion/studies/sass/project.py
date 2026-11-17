"""A reproducible provisional course-project study assembled from cited modules.

References
----------
.. [1] ITU-UZM508. Istanbul Technical University, "UZM 508E - Spacecraft
       Architecture and Subsystems Design"; course identity only; the scenario is
       not an official assignment.
.. [2] NASA-SE-APP. NASA, "Systems Engineering Handbook: Appendix"; Trade Study
       Report.
.. [3] NASA-SE-6. NASA, "Systems Engineering Handbook: Crosscutting Technical
       Management"; interface and decision data.
"""

from __future__ import annotations

import json
from collections.abc import Mapping
from dataclasses import asdict, dataclass, replace
from fractions import Fraction
from math import fsum, nextafter
from typing import Literal

from rocket_propulsion.core.errors import DomainError
from rocket_propulsion.core.validation import integer, real, records, text
from rocket_propulsion.propulsion.burns.serialization import canonical_sha256

from .catalog import (
    CatalogueFlowCheck,
    ThrusterRecord,
    _catalogue_flow_check,
    _nominal_study_point,
    thruster,
)
from .command_life import (
    CommandLifeReview,
    ThrusterCommandHistory,
    ValveCommand,
    review_thruster_command_life,
)
from .evidence import AnalysisResult, Citation, CitedScalar, derived, quantity
from .pressure_interface import (
    CataloguePressureCheck,
    InletPressureWindow,
    _catalogue_pressure_checks,
    inlet_ports,
)
from .references import REFERENCES
from .sizing import DeltaVItem, Maneuver, MissionSizing, delta_v_budget, size_mission
from .trade import AlternativeScore, Criterion, Requirement, weighted_trade

STUDY_SCHEMA = "rocket_propulsion_sass_study_v1"


@dataclass(frozen=True, slots=True)
class Candidate:
    """One propulsion alternative: installed hardware and operating assumptions.

    ``other_propulsion_dry_mass`` excludes the catalogue head mass, which is added
    once per installed head [1]_. ``firing_bus_power`` is optional because
    manufacturer valve power is not complete propulsion bus power; a missing value
    remains unknown in the screens [2]_.

    Attributes
    ----------
    identifier : str
        Unique alternative name.
    catalog_key : str
        Catalogue variant identifier.
    other_propulsion_dry_mass : CitedScalar
        Tanks, feed, pressurant inventory, electronics and mounting, kg.
    installed_count : int
        Installed heads, including non-firing redundant heads.
    firing_count : int
        Simultaneously firing heads; 1 <= firing_count <= installed_count.
    axial_efficiency : CitedScalar
        Useful fraction of scalar thrust, in (0, 1].
    duty_fraction : CitedScalar
        On-time divided by elapsed time, in (0, 1].
    firing_bus_power : CitedScalar or None
        Complete bus power while firing, W.

    Raises
    ------
    DomainError
        If the catalogue key is unknown, counts are inconsistent or a quantity is out
        of domain.

    References
    ----------
    .. [1] NASA-SOA26. NASA, "State-of-the-Art of Small Spacecraft Technology: In-Space
           Propulsion"; component versus integrated-system boundary.
    .. [2] NASA-SE-6. NASA, "Systems Engineering Handbook: Crosscutting Technical
           Management"; sections 6.3 and 6.8, interface budgets and explicit
           assumptions.
    """

    identifier: str
    catalog_key: str
    other_propulsion_dry_mass: CitedScalar
    installed_count: int
    firing_count: int
    axial_efficiency: CitedScalar
    duty_fraction: CitedScalar
    firing_bus_power: CitedScalar | None = None
    command_history: ThrusterCommandHistory | None = None
    inlet_pressures: tuple[InletPressureWindow, ...] = ()

    def __post_init__(self) -> None:
        """Check configuration identity and counts without filling missing hardware mass.

        References
        ----------
        .. [1] NASA-SE-6. NASA, "Systems Engineering Handbook: Crosscutting Technical
               Management"; sections 6.3, 6.5 and 6.8, interfaces and configuration.
        """
        text(self.identifier, "candidate identifier")
        record = thruster(self.catalog_key)
        installed = integer(self.installed_count, "installed engine count")
        firing = integer(self.firing_count, "firing engine count")
        real(installed, "finite installed count")
        real(firing, "finite firing count")
        if not installed >= firing >= 1:
            raise DomainError("Require installed engine count >= firing engine count >= 1.")
        object.__setattr__(self, "installed_count", installed)
        object.__setattr__(self, "firing_count", firing)
        quantity(self.other_propulsion_dry_mass, "kg", minimum=0)
        for value in (self.axial_efficiency, self.duty_fraction):
            if quantity(value, "1", strictly_positive=True) > 1:
                raise DomainError("Candidate axial efficiency and duty fraction must be <=1.")
        if self.firing_bus_power is not None:
            quantity(self.firing_bus_power, "W", minimum=0)
        history = self.command_history
        if history is not None and (not isinstance(history, ThrusterCommandHistory)
                                    or len(history.head_ids) != installed):
            raise DomainError("Command history must identify every installed head in the candidate.")
        pressures: tuple[InletPressureWindow, ...] = records(
            self.inlet_pressures, InletPressureWindow, "inlet pressure windows"
        )
        ports = {window.port for window in pressures}
        if len(ports) != len(pressures) or not ports.issubset(inlet_ports(record)):
            raise DomainError("Inlet pressure windows must identify unique ports on the selected thruster.")
        if len({window.condition_record for window in pressures}) > 1:
            raise DomainError("Inlet pressure windows must refer to the same operating condition.")
        object.__setattr__(self, "inlet_pressures", pressures)

    def input_record(self) -> dict[str, object]:
        """Keep an optional command history in the immutable request snapshot.

        References
        ----------
        .. [1] NASA-SE-6. NASA Systems Engineering Handbook, 6.6.
        """
        record = asdict(self)
        if self.command_history is None:
            del record["command_history"]
        if not self.inlet_pressures:
            del record["inlet_pressures"]
        return record


@dataclass(frozen=True, slots=True)
class CatalogueLifeCheck:
    """A demand comparison with a published per-thruster reference.

    A cumulative demand is a lower bound on the worst head under ideal
    sharing across all installed heads. A value within the reference
    leaves the actual duty schedule and qualification applicability open.

    References
    ----------
    .. [1] ARIANE-1N. ArianeGroup, "1 N Monopropellant Hydrazine Thruster";
           published burn, throughput and impulse reference values.
    .. [2] NASA-SE-6. NASA, "Systems Engineering Handbook"; preliminary screening.
    """

    metric: str
    demand: CitedScalar | None
    published_limit: CitedScalar | None
    status: Literal["exceeds_reference", "no_exceedance_identified", "unassessed"]
    basis: str


def _catalogue_life_checks(
    config: Candidate, product: ThrusterRecord, sizing: MissionSizing | None,
    items: tuple[DeltaVItem, ...],
) -> tuple[CatalogueLifeCheck, ...]:
    """Screen stated operating demands against applicable catalogue references.

    References
    ----------
    .. [1] ARIANE-1N. ArianeGroup, "1 N Monopropellant Hydrazine Thruster";
           qualification table.
    .. [2] ARIANE-10N. ArianeGroup, "10 N Bipropellant Thruster";
           qualified single, accumulated and cycle life.
    .. [3] NASA-SE-6. NASA, "Systems Engineering Handbook"; section 6.8,
           retain limitations of decision-analysis inputs.
    """
    limits = {field.name: field for field in product.fields}
    summary = None if sizing is None else sizing.summary
    # Individual firing durations are available for a single, continuous
    # occurrence with every installed head firing for the full duration.
    known_single_burns: list[CitedScalar] = []
    singles_complete = False
    if (sizing is not None and config.duty_fraction.value == 1.0
            and config.installed_count == config.firing_count):
        known_single_burns = [
            burn.get("on_time") for item, (_, burn) in zip(items, sizing.burns, strict=True)
            if item.count == 1
        ]
        singles_complete = bool(known_single_burns) and len(known_single_burns) == len(items)
    longest = max(known_single_burns, key=lambda value: value.value, default=None)
    installed = config.installed_count
    specs = (
        ("single_burn_life", "s", longest, Fraction(1), singles_complete,
         "Longest resolved single continuous burn with all installed heads firing."),
        ("accumulated_burn_life", "s", None if summary is None else summary.get("on_time"),
         Fraction(config.firing_count, installed), True,
         "Mission on-time × simultaneous firing heads ÷ installed heads."),
        ("propellant_throughput", "kg", None if summary is None else summary.get("consumed_propellant"),
         Fraction(1, installed), True, "Consumed propellant ÷ installed heads."),
        ("total_impulse_life", "N*s", None if summary is None else summary.get("scalar_impulse"),
         Fraction(1, installed), True, "Total scalar impulse ÷ installed heads."),
        ("cycle_life", "1", None, Fraction(1), False,
         "Cycle count requires a per-head command history."),
    )
    checks = []
    for name, unit, source, scale, complete, basis in specs:
        field = limits.get(name)
        limit = None if field is None else field.upper
        if field is not None:
            if field.kind != "upper_bound":
                raise DomainError(f"Catalogue {name} must be a published upper reference.")
            quantity(field.upper, unit, strictly_positive=True)
        demand = None
        status: Literal["exceeds_reference", "no_exceedance_identified", "unassessed"] = "unassessed"
        if source is not None:
            value = quantity(source, unit, minimum=0)
            exact_demand = Fraction(value) * scale
            demand = derived(real(exact_demand, name + " demand"), unit, (source,), (
                Citation("NASA-SE-6", "6.8: screening with declared operating assumptions"),
            ), basis)
            if limit is not None:
                # A result at the binary rounding boundary does not establish
                # exceedance from the exported mass and time ledger.
                if Fraction(nextafter(value, 0.0)) * scale > Fraction(limit.value):
                    status = "exceeds_reference"
                elif complete:
                    status = "no_exceedance_identified"
        if name == "single_burn_life" and not complete:
            basis += " Some individual durations require a command and head assignment schedule."
        checks.append(CatalogueLifeCheck(name, demand, limit, status, basis))
    return tuple(checks)


@dataclass(frozen=True, slots=True)
class CandidateStudy:
    """Computed budgets of one alternative, or the reasons it cannot yet be sized.

    Attributes
    ----------
    candidate : Candidate
        The evaluated alternative.
    sizing : MissionSizing or None
        Mission mass budget; None when no nominal operating point exists.
    operating_point : AnalysisResult or None
        Published nominal study point.
    metrics : tuple of (str, CitedScalar)
        Metrics offered to the trade study.
    unresolved : tuple of str
        Catalogue limitations and missing evidence.
    catalog_record : ThrusterRecord
        Immutable catalogue record used for this analysis.
    catalogue_life_checks : tuple of CatalogueLifeCheck
        Demand comparisons with published per-thruster references.
    command_life_review : CommandLifeReview or None
        Per-head observations and comparisons when a command history is supplied.
    catalogue_flow_check : CatalogueFlowCheck or None
        Calculated nominal flow compared with the reported range when available.

    References
    ----------
    .. [1] NASA-SE-APP. NASA, "Systems Engineering Handbook: Appendix"; Trade Study
           Report, retain all considered alternatives.
    """

    candidate: Candidate
    sizing: MissionSizing | None
    operating_point: AnalysisResult | None
    metrics: tuple[tuple[str, CitedScalar], ...]
    unresolved: tuple[str, ...]
    catalog_record: ThrusterRecord
    catalogue_life_checks: tuple[CatalogueLifeCheck, ...] = ()
    command_life_review: CommandLifeReview | None = None
    catalogue_flow_check: CatalogueFlowCheck | None = None
    catalogue_pressure_checks: tuple[CataloguePressureCheck, ...] = ()


@dataclass(frozen=True, slots=True)
class ProjectStudy:
    """A dated, source-complete preparation study; not a graded course submission.

    Attributes
    ----------
    title : str
        Study title.
    budget : AnalysisResult
        Mission delta-v allocation.
    candidates : tuple of CandidateStudy
        Per-alternative results in input order.
    scores : tuple of AlternativeScore
        Trade results in input order.
    input_json : str
        Canonical JSON snapshot of the request.

    Notes
    -----
    The course reference identifies the course only [1]_; the study follows the trade
    study report structure [2]_.

    References
    ----------
    .. [1] ITU-UZM508. Istanbul Technical University, "UZM 508E - Spacecraft
           Architecture and Subsystems Design"; course identity.
    .. [2] NASA-SE-APP. NASA, "Systems Engineering Handbook: Appendix"; Trade Study
           Report.
    """

    title: str
    budget: AnalysisResult
    candidates: tuple[CandidateStudy, ...]
    scores: tuple[AlternativeScore, ...]
    input_json: str

    @property
    def input_data(self) -> dict[str, object]:
        """Return a detached copy of the inputs so callers cannot rewrite provenance.

        Returns
        -------
        dict

        References
        ----------
        .. [1] NASA-SE-6. NASA, "Systems Engineering Handbook: Crosscutting Technical
               Management"; sections 6.5/6.6, configuration and technical data integrity.
        """
        return json.loads(self.input_json)

    def to_dict(self) -> dict[str, object]:
        """Export inputs, catalogue variants, sources, results and a snapshot hash.

        Returns
        -------
        dict
            JSON-native study document with an ``assurance`` block stating its permitted
            use and a ``study_hash`` (canonical SHA-256 of the remaining payload).

        Notes
        -----
        The hash identifies this output snapshot [2]_; no verified-document reader or
        independent qualification certificate is implied. All alternatives, results and
        data sources are retained [1]_.

        References
        ----------
        .. [1] NASA-SE-APP. NASA, "Systems Engineering Handbook: Appendix"; Trade Study
               Report, alternatives, computational results and data sources.
        .. [2] NASA-SE-6. NASA, "Systems Engineering Handbook: Crosscutting Technical
               Management"; section 6.5, configuration identification.
        """
        candidate_rows = []
        for item in self.candidates:
            candidate_rows.append({
                "configuration": item.candidate.input_record(),
                "catalog_record": asdict(item.catalog_record),
                "operating_point": None if item.operating_point is None
                else item.operating_point.to_dict(),
                "sizing": None if item.sizing is None else item.sizing.to_dict(),
                "metrics": {name: asdict(value) for name, value in item.metrics},
                "unresolved": list(item.unresolved),
                "catalogue_life_checks": [asdict(check) for check in item.catalogue_life_checks],
                **({} if item.catalogue_flow_check is None else
                   {"catalogue_flow_check": asdict(item.catalogue_flow_check)}),
                **({} if not item.catalogue_pressure_checks else
                   {"catalogue_pressure_checks": [asdict(check) for check in item.catalogue_pressure_checks]}),
                **({} if item.command_life_review is None else
                   {"command_life_review": asdict(item.command_life_review)}),
            })
        payload = {
            "schema": STUDY_SCHEMA, "title": self.title,
            "course": "UZM 508E - provisional propulsion preparation",
            "assignment_status": "unconfirmed_scope_not_an_official_course_brief",
            "assurance": {
                "permitted_use": "preliminary_design_screening_only",
                "hardware_scope": "motor_and_feed_system",
                "hardware_readiness": "not_assessed",
                "experimental_validation": "not_performed",
                "qualification_evidence": False,
                "uncertainty": "not_quantified; scores and passes are deterministic screens only",
                "ranking_robustness": "not_assessed; numerical precision does not resolve uncertain inputs",
                "catalog_evidence": "published nominal values; not validated joint operating maps",
                "pulse_life_evidence": (
                    "commanded cycles compared for supplied observation windows; "
                    "actuator response and mission coverage unassessed"
                    if any(item.command_life_review is not None for item in self.candidates)
                    else "not_assessed; budget occurrences are aggregated"
                ),
                "catalogue_life_evidence": "preliminary comparison with published per-thruster references",
            },
            "inputs": dict(self.input_data), "delta_v_budget": self.budget.to_dict(),
            "candidates": candidate_rows, "trade": [asdict(x) for x in self.scores],
            "references": [asdict(x) for x in REFERENCES],
            "extrapolation_policy": "disabled; no inferred catalogue performance points",
            "selection_rule": "hard screens and catalogue life review precede ranking; "
                              "published life exceedances require review; otherwise exact affine "
                              "weighted utilities use original metrics, anchors and weights",
        }
        # JSON-native lists make Python and persisted output contracts identical.
        payload = json.loads(json.dumps(payload, allow_nan=False))
        return {**payload, "study_hash": canonical_sha256(payload)}


def analyze_project(title: str, bus_dry_mass: CitedScalar, items: tuple[DeltaVItem, ...],
                     candidates: tuple[Candidate, ...], residual_fraction: CitedScalar,
                     fixed_reserve: CitedScalar, requirements: tuple[Requirement, ...],
                     criteria: tuple[Criterion, ...]) -> ProjectStudy:
    """Build comparable mission budgets for each alternative and score the trade.

    Parameters
    ----------
    title : str
        Study title.
    bus_dry_mass : CitedScalar
        Spacecraft dry mass excluding propulsion, kg; positive.
    items : tuple of DeltaVItem
        Delta-v allocations, executed in order.
    candidates : tuple of Candidate
        Non-empty, uniquely identified alternatives.
    residual_fraction : CitedScalar
        Residual fraction of loaded propellant.
    fixed_reserve : CitedScalar
        Additional retained propellant, kg.
    requirements : tuple of Requirement
        Hard screens.
    criteria : tuple of Criterion
        Weighted trade criteria.

    Returns
    -------
    ProjectStudy

    Raises
    ------
    DomainError
        If an input is invalid, a sized candidate lacks a nominal hardware mass, or a
        mass budget does not close.

    Notes
    -----
    Only catalogue records with a jointly labelled nominal thrust/Isp pair are sized;
    others stay unresolved and unscored. Each allocated budget line is executed at
    the candidate's constant study point using the rocket-equation and impulse
    identities [1]_ [2]_. Dry mass is bus mass plus other propulsion dry mass plus
    installed heads times the catalogue head mass, so the head mass is added exactly
    once [3]_ [4]_.

    References
    ----------
    .. [1] NASA-ROCKET. NASA Glenn Research Center, "Ideal Rocket Equation"; mission
           mass identities.
    .. [2] NASA-ISP. NASA Glenn Research Center, "Specific Impulse"; mission impulse
           identities.
    .. [3] NASA-SOA26. NASA, "State-of-the-Art of Small Spacecraft Technology: In-Space
           Propulsion"; component versus complete system scope.
    .. [4] NASA-SE-6. NASA, "Systems Engineering Handbook: Crosscutting Technical
           Management"; sections 6.3 and 6.8, interfaces and explicit decision criteria.
    """
    text(title, "study title")
    bus = quantity(bus_dry_mass, "kg", strictly_positive=True)
    entries = records(items, DeltaVItem, "project delta-v items")
    configurations = records(candidates, Candidate, "project candidates")
    if not configurations or len({x.identifier for x in configurations}) != len(configurations):
        raise DomainError("Project candidates must be non-empty and uniquely identified.")
    budget = delta_v_budget(entries)
    studies = []
    alternatives = {}
    catalogue_review_needed: set[str] = set()
    for config in configurations:
        product = thruster(config.catalog_key)
        operating = _nominal_study_point(product)
        flow_check = _catalogue_flow_check(product, operating)
        pressure_checks = _catalogue_pressure_checks(product, config.inlet_pressures)
        sizing = None
        metrics = {}
        unresolved = list(product.limitations)
        if flow_check is not None and flow_check.status in ("above_published_range", "below_published_range"):
            unresolved.append("Calculated nominal mass flow falls outside the published range. Review before ranking.")
            catalogue_review_needed.add(config.identifier)
        for check in pressure_checks:
            if check.status not in ("matches_nominal", "within_published_range"):
                unresolved.append(f"Review required for {check.port}: {check.reason}")
                catalogue_review_needed.add(config.identifier)
            if check.meop is not None and check.meop.status != "within_meop_reference":
                unresolved.append(f"Operating pressure review required for {check.port}. {check.meop.reason}")
                catalogue_review_needed.add(config.identifier)
        # A record containing bounds/ranges is useful evidence but not a model.
        if operating is None:
            unresolved.append("Missing jointly labelled nominal operating point; no sizing or score.")
        else:
            component_mass = product.field("dry_mass").upper
            if product.field("dry_mass").kind != "nominal":
                raise DomainError("A nominal installed hardware mass is required for project sizing.")
            installed_mass = quantity(component_mass, "kg", strictly_positive=True)
            dry = derived(bus + config.other_propulsion_dry_mass.value
                          + config.installed_count * installed_mass, "kg",
                          (bus_dry_mass, config.other_propulsion_dry_mass, component_mass), (
                              Citation("NASA-SE-6", "6.3: non-overlapping propulsion/spacecraft mass ledger"),
                          ), "Installed count includes non-firing redundant heads.")
            force = derived(config.firing_count * operating.get("thrust").value, "N",
                            (operating.get("thrust"),), (
                                Citation("NASA-ISP", "Add scalar thrust of simultaneously firing heads"),
                            ), f"{config.firing_count} simultaneously firing heads.")
            burns = tuple(Maneuver(x.name, budget.get(x.name), force, operating.get("isp"),
                                   config.axial_efficiency, config.duty_fraction) for x in entries)
            sizing = size_mission(dry, burns, residual_fraction=residual_fraction,
                                 fixed_reserve=fixed_reserve)
            metrics = dict(sizing.summary.values)
            if config.firing_bus_power is not None:
                metrics["firing_bus_power"] = config.firing_bus_power
            else:
                unresolved.append("Complete firing bus power is missing; valve-only power is not substituted.")
        life_checks = _catalogue_life_checks(config, product, sizing, entries)
        exceeded = [check.metric for check in life_checks if check.status == "exceeds_reference"]
        if exceeded:
            unresolved.append("Review required for published life exceedance: " + ", ".join(exceeded) + ".")
            catalogue_review_needed.add(config.identifier)
        command_review = None
        if config.command_history is not None:
            command_review = review_thruster_command_life(product, config.command_history)
            if command_review.peak_simultaneous_open > config.firing_count:
                unresolved.append("Command history exceeds the sizing assumption for simultaneous heads.")
                catalogue_review_needed.add(config.identifier)
            if any(check.status == "exceeds_reference" for check in command_review.comparisons):
                unresolved.append("Review required for per-head command history life exceedance.")
                catalogue_review_needed.add(config.identifier)
        alternatives[config.identifier] = metrics
        studies.append(CandidateStudy(config, sizing, operating, tuple(metrics.items()),
                                      tuple(unresolved), product, life_checks, command_review, flow_check, pressure_checks))
    scores = tuple(
        replace(score, eligibility="incomplete", score=None)
        if score.eligibility == "eligible" and score.name in catalogue_review_needed
        else score
        for score in weighted_trade(alternatives, criteria, requirements)
    )
    inputs = {
        "schema": STUDY_SCHEMA, "title": title,
        "bus_dry_mass_excluding_propulsion": asdict(bus_dry_mass),
        "items": [asdict(x) for x in entries], "candidates": [x.input_record() for x in configurations],
        "residual_fraction": asdict(residual_fraction), "fixed_reserve": asdict(fixed_reserve),
        "requirements": [asdict(x) for x in requirements], "criteria": [asdict(x) for x in criteria],
    }
    return ProjectStudy(title, budget, tuple(studies), scores,
                        json.dumps(inputs, sort_keys=True, allow_nan=False))


def _object(value: object, fields: set[str], label: str,
             optional: set[str] | None = None) -> Mapping:
    """Require an exact set of JSON fields so misspelled assumptions cannot disappear.

    References
    ----------
    .. [1] NASA-SE-6. NASA, "Systems Engineering Handbook: Crosscutting Technical
           Management"; section 6.6, controlled and interpretable technical input data.
    """
    if not isinstance(value, Mapping) or set(value) - fields - (optional or set()):
        raise DomainError(f"{label}: expected an object with known fields.")
    if fields - set(value):
        raise DomainError(f"{label}: missing required fields {sorted(fields - set(value))}.")
    return value


def scalar_from_dict(value: object) -> CitedScalar:
    """Read a cited-scalar envelope; a bare number is deliberately not accepted.

    Parameters
    ----------
    value : Mapping
        Object with ``value``, ``unit``, ``citations``, ``status`` and ``note``, and
        optionally ``assumption_notes``.

    Returns
    -------
    CitedScalar

    Raises
    ------
    DomainError
        If fields are missing or unknown, citations are not an array, or the
        quantity is invalid.

    References
    ----------
    .. [1] NASA-SE-6. NASA, "Systems Engineering Handbook: Crosscutting Technical
           Management"; section 6.6, evidence-bearing technical data.
    """
    raw = _object(value, {"value", "unit", "citations", "status", "note"}, "cited scalar",
                  {"assumption_notes"})
    if not isinstance(raw["citations"], (list, tuple)):
        raise DomainError("Quantity citations must be a JSON array.")
    citations = []
    for item in raw["citations"]:
        source = _object(item, {"source_id", "locator"}, "citation")
        citations.append(Citation(source["source_id"], source["locator"]))
    return CitedScalar(raw["value"], raw["unit"], tuple(citations), raw["status"], raw["note"],
                        raw.get("assumption_notes", ()))


def project_from_dict(payload: object) -> ProjectStudy:
    """Execute a strict, citation-bearing v1 project request with no hidden defaults.

    Parameters
    ----------
    payload : Mapping
        Request with exactly the fields ``schema``, ``title``,
        ``bus_dry_mass_excluding_propulsion``, ``items``, ``candidates``,
        ``residual_fraction``, ``fixed_reserve``, ``requirements`` and ``criteria``.

    Returns
    -------
    ProjectStudy

    Raises
    ------
    DomainError
        If the schema is unsupported, a field is missing, unknown or malformed, no
        positive trade weight is given, or the analysis fails.

    References
    ----------
    .. [1] NASA-SE-APP. NASA, "Systems Engineering Handbook: Appendix"; reproducible
           trade study inputs, selection rule and sources.
    """
    raw = _object(payload, {"schema", "title", "bus_dry_mass_excluding_propulsion", "items",
                            "candidates", "residual_fraction", "fixed_reserve", "requirements",
                            "criteria"}, "SASS project")
    if raw["schema"] != STUDY_SCHEMA:
        raise DomainError("Unsupported SASS study request schema.")
    for label in ("items", "candidates", "requirements", "criteria"):
        if not isinstance(raw[label], list):
            raise DomainError(f"{label} must be a JSON array.")
    items = []
    for value in raw["items"]:
        item = _object(value, {"name", "nominal_delta_v", "margin_fraction", "count"}, "budget item")
        items.append(DeltaVItem(item["name"], scalar_from_dict(item["nominal_delta_v"]),
                               scalar_from_dict(item["margin_fraction"]), item["count"]))
    candidates = []
    for value in raw["candidates"]:
        item = _object(value, {"identifier", "catalog_key", "other_propulsion_dry_mass",
                               "installed_count", "firing_count", "axial_efficiency", "duty_fraction",
                               "firing_bus_power"}, "candidate", {"command_history", "inlet_pressures"})
        candidates.append(Candidate(
            item["identifier"], item["catalog_key"], scalar_from_dict(item["other_propulsion_dry_mass"]),
            item["installed_count"], item["firing_count"], scalar_from_dict(item["axial_efficiency"]),
            scalar_from_dict(item["duty_fraction"]), None if item["firing_bus_power"] is None
            else scalar_from_dict(item["firing_bus_power"]),
            None if item.get("command_history") is None else _command_history_from_dict(item["command_history"]),
            _inlet_pressures_from_dict(item.get("inlet_pressures", [])),
        ))
    requirements = []
    for value in raw["requirements"]:
        item = _object(value, {"identifier", "description", "metric", "relation", "limit",
                               "verification"}, "requirement")
        requirements.append(Requirement(item["identifier"], item["description"], item["metric"],
                                        item["relation"], scalar_from_dict(item["limit"]),
                                        item["verification"]))
    criteria = []
    for value in raw["criteria"]:
        item = _object(value, {"metric", "best", "worst", "weight"}, "criterion")
        criteria.append(Criterion(item["metric"], scalar_from_dict(item["best"]),
                                  scalar_from_dict(item["worst"]), scalar_from_dict(item["weight"])))
    # fsum is used by the trade module; reject an empty request before any result is emitted.
    if fsum(x.weight.value for x in criteria) <= 0:
        raise DomainError("At least one positive trade weight is required.")
    return analyze_project(raw["title"], scalar_from_dict(raw["bus_dry_mass_excluding_propulsion"]),
                           tuple(items), tuple(candidates), scalar_from_dict(raw["residual_fraction"]),
                           scalar_from_dict(raw["fixed_reserve"]), tuple(requirements), tuple(criteria))


def _inlet_pressures_from_dict(value: object) -> tuple[InletPressureWindow, ...]:
    """Read cited pressure windows with no implicit datum or ambient reference.

    References
    ----------
    .. [1] NASA-SE-6, 6.3 and 6.6, traceable interface data.
    """
    if not isinstance(value, list):
        raise DomainError("Inlet pressure windows must be a JSON array.")
    windows = []
    for row in value:
        data = _object(row, {"port", "reference", "lower", "upper", "condition_record"},
                       "inlet pressure window",
                       {"gauge_reference_lower", "gauge_reference_upper", "gauge_reference_record"})
        ref_low, ref_high = data.get("gauge_reference_lower"), data.get("gauge_reference_upper")
        windows.append(InletPressureWindow(
            data["port"], data["reference"], scalar_from_dict(data["lower"]),
            scalar_from_dict(data["upper"]), data["condition_record"],
            None if ref_low is None else scalar_from_dict(ref_low),
            None if ref_high is None else scalar_from_dict(ref_high), data.get("gauge_reference_record"),
        ))
    return tuple(windows)


def _command_history_from_dict(value: object) -> ThrusterCommandHistory:
    """Read the optional per-head command history from a project request.

    References
    ----------
    .. [1] NASA-SE-6. NASA Systems Engineering Handbook, 6.6.
    """
    fields = {"source_record", "tick_period", "duration_ticks", "head_ids", "initially_open", "commands"}
    raw = _object(value, fields, "thruster command history")
    for key in ("head_ids", "initially_open", "commands"):
        if not isinstance(raw[key], list):
            raise DomainError(f"Thruster command history {key} must be a JSON array.")
    commands = []
    for command in raw["commands"]:
        event = _object(command, {"tick", "head_id", "opened"}, "valve command")
        commands.append(ValveCommand(event["tick"], event["head_id"], event["opened"]))
    return ThrusterCommandHistory(
        raw["source_record"], scalar_from_dict(raw["tick_period"]), raw["duration_ticks"],
        tuple(raw["head_ids"]), tuple(raw["initially_open"]), tuple(commands),
    )


def project_from_json(value: str) -> ProjectStudy:
    """Parse a bounded JSON request with the strict reader and execute it.

    Parameters
    ----------
    value : str
        JSON request text.

    Returns
    -------
    ProjectStudy

    Raises
    ------
    InputError
        If the text is not strict JSON (duplicate fields, non-finite constants,
        oversize or malformed input).
    DomainError
        If the request content is invalid.

    References
    ----------
    .. [1] NASA-SE-6. NASA, "Systems Engineering Handbook: Crosscutting Technical
           Management"; section 6.6, controlled technical inputs and explicit error
           handling.
    """
    from .strict_json import read_project_document
    return read_project_document(value, project_from_dict)
