"""A reproducible provisional course-project study assembled from cited modules.

References
----------
ITU-UZM508: course identity only; the scenario is not an official assignment.
NASA-SE-APP: Trade Study Report; NASA-SE-6: interface and decision data.
"""

from __future__ import annotations

import json
from collections.abc import Mapping
from dataclasses import asdict, dataclass
from math import fsum

from rocket_propulsion.core.errors import DomainError
from rocket_propulsion.core.validation import integer, real, records, text
from rocket_propulsion.propulsion.burns.serialization import canonical_sha256

from .catalog import nominal_operating_point, thruster
from .evidence import AnalysisResult, Citation, CitedScalar, derived, quantity
from .references import REFERENCES
from .sizing import DeltaVItem, Maneuver, MissionSizing, delta_v_budget, size_mission
from .trade import AlternativeScore, Criterion, Requirement, weighted_trade

STUDY_SCHEMA = "rocket_propulsion_sass_study_v1"


@dataclass(frozen=True, slots=True)
class Candidate:
    """Installed hardware and operations; other dry mass excludes catalogue head mass.

    ``firing_bus_power`` is optional because manufacturer valve power is not
    complete propulsion bus power. Missing values remain unknown in screens.

    References
    ----------
    NASA-SOA26: component versus integrated-system boundary.
    NASA-SE-6: sections 6.3 and 6.8, interface budgets and explicit assumptions.
    """

    identifier: str
    catalog_key: str
    other_propulsion_dry_mass: CitedScalar
    installed_count: int
    firing_count: int
    axial_efficiency: CitedScalar
    duty_fraction: CitedScalar
    firing_bus_power: CitedScalar | None = None

    def __post_init__(self) -> None:
        """Check configuration identity and counts without filling missing hardware mass.

        References
        ----------
        NASA-SE-6: sections 6.3, 6.5 and 6.8, interfaces and configuration.
        """
        text(self.identifier, "candidate identifier")
        thruster(self.catalog_key)
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


@dataclass(frozen=True, slots=True)
class CandidateStudy:
    """Computed budgets or an explicit reason why a candidate cannot yet be sized.

    References
    ----------
    NASA-SE-APP: Trade Study Report, retain all considered alternatives.
    """

    candidate: Candidate
    sizing: MissionSizing | None
    operating_point: AnalysisResult | None
    metrics: tuple[tuple[str, CitedScalar], ...]
    unresolved: tuple[str, ...]


@dataclass(frozen=True, slots=True)
class ProjectStudy:
    """A dated, source-complete preparation artifact, not a graded course submission.

    References
    ----------
    ITU-UZM508: course identity; NASA-SE-APP: Trade Study Report.
    """

    title: str
    budget: AnalysisResult
    candidates: tuple[CandidateStudy, ...]
    scores: tuple[AlternativeScore, ...]
    input_json: str

    @property
    def input_data(self) -> dict[str, object]:
        """Return a detached input copy so a caller cannot rewrite calculation provenance.

        References
        ----------
        NASA-SE-6: sections 6.5/6.6, configuration and technical data integrity.
        """
        return json.loads(self.input_json)

    def to_dict(self) -> dict[str, object]:
        """Export complete inputs, catalogue variants, sources and computed results.

        The hash identifies this output snapshot; no verified-document reader
        or independent qualification certificate is claimed.

        References
        ----------
        NASA-SE-APP: Trade Study Report, alternatives, computational results and data sources.
        NASA-SE-6: section 6.5, configuration identification.
        """
        candidate_rows = []
        for item in self.candidates:
            candidate_rows.append({
                "configuration": asdict(item.candidate),
                "catalog_record": asdict(thruster(item.candidate.catalog_key)),
                "operating_point": None if item.operating_point is None
                else item.operating_point.to_dict(),
                "sizing": None if item.sizing is None else item.sizing.to_dict(),
                "metrics": {name: asdict(value) for name, value in item.metrics},
                "unresolved": list(item.unresolved),
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
                "pulse_life_evidence": "not_assessed; budget occurrences are aggregated",
            },
            "inputs": dict(self.input_data), "delta_v_budget": self.budget.to_dict(),
            "candidates": candidate_rows, "trade": [asdict(x) for x in self.scores],
            "references": [asdict(x) for x in REFERENCES],
            "extrapolation_policy": "disabled; no inferred catalogue performance points",
            "selection_rule": "hard screens then exact affine weighted utilities from original "
                              "metrics/anchors/weights; final score rounded once; incomplete rows unranked",
        }
        # JSON-native lists make Python and persisted output contracts identical.
        payload = json.loads(json.dumps(payload, allow_nan=False))
        return {**payload, "study_hash": canonical_sha256(payload)}


def analyze_project(title: str, bus_dry_mass: CitedScalar, items: tuple[DeltaVItem, ...],
                     candidates: tuple[Candidate, ...], residual_fraction: CitedScalar,
                     fixed_reserve: CitedScalar, requirements: tuple[Requirement, ...],
                     criteria: tuple[Criterion, ...]) -> ProjectStudy:
    """Build comparable mission budgets from actual nominal catalog pairs only.

    Every allocated budget line is executed in order at one candidate's
    constant study point. Other propulsion dry mass includes tanks, feed,
    pressurant inventory, electronics and mounting as an explicit allocation;
    the catalogue thruster-head mass is added exactly once.

    References
    ----------
    NASA-ROCKET / NASA-ISP: mission mass and impulse identities.
    NASA-SOA26: component versus complete system scope.
    NASA-SE-6: sections 6.3 and 6.8, interfaces and explicit decision criteria.
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
    for config in configurations:
        product = thruster(config.catalog_key)
        operating = None
        sizing = None
        metrics = {}
        unresolved = list(product.limitations)
        # A record containing bounds/ranges is useful evidence but not a model.
        if product.field("thrust").kind != "nominal" or product.field("isp").kind != "nominal":
            unresolved.append("Missing jointly labelled nominal operating point; no sizing or score.")
        else:
            operating = nominal_operating_point(config.catalog_key)
            component_mass = product.field("dry_mass").upper
            if product.field("dry_mass").kind != "nominal":
                raise DomainError("A nominal installed hardware mass is required for project sizing.")
            dry = derived(bus + config.other_propulsion_dry_mass.value
                          + config.installed_count * component_mass.value, "kg",
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
        alternatives[config.identifier] = metrics
        studies.append(CandidateStudy(config, sizing, operating, tuple(metrics.items()),
                                      tuple(unresolved)))
    scores = weighted_trade(alternatives, criteria, requirements)
    inputs = {
        "schema": STUDY_SCHEMA, "title": title,
        "bus_dry_mass_excluding_propulsion": asdict(bus_dry_mass),
        "items": [asdict(x) for x in entries], "candidates": [asdict(x) for x in configurations],
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
    NASA-SE-6: section 6.6, controlled and interpretable technical input data.
    """
    if not isinstance(value, Mapping) or set(value) - fields - (optional or set()):
        raise DomainError(f"{label}: expected an object with known fields.")
    if fields - set(value):
        raise DomainError(f"{label}: missing required fields {sorted(fields - set(value))}.")
    return value


def scalar_from_dict(value: object) -> CitedScalar:
    """Read the evidence envelope; a bare number is deliberately not accepted.

    References
    ----------
    NASA-SE-6: section 6.6, evidence-bearing technical data.
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
    """Execute a strict, citation-bearing v1 project request; no hidden defaults.

    References
    ----------
    NASA-SE-APP: reproducible trade study inputs, selection rule and sources.
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
                               "firing_bus_power"}, "candidate")
        candidates.append(Candidate(
            item["identifier"], item["catalog_key"], scalar_from_dict(item["other_propulsion_dry_mass"]),
            item["installed_count"], item["firing_count"], scalar_from_dict(item["axial_efficiency"]),
            scalar_from_dict(item["duty_fraction"]), None if item["firing_bus_power"] is None
            else scalar_from_dict(item["firing_bus_power"]),
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


def project_from_json(value: str) -> ProjectStudy:
    """Parse a bounded JSON document through the shared strict reader.

    References
    ----------
    NASA-SE-6: section 6.6, controlled technical inputs and explicit error handling.
    """
    from .strict_json import read_project_document
    return read_project_document(value, project_from_dict)
