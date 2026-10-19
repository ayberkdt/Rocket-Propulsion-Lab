"""Small, field-cited spacecraft thruster catalog; no envelope-to-point inference.

References
----------
MOOG-COLD, ARIANE-1N, ARIANE-10N: manufacturer nominal characteristics.
NASA-SOA26: ECAPS row, ranges only. BUSEK-BIT3 / BUSEK-RF: unresolved page discrepancy.
NASA-SE-6: keep data provenance and system/component boundaries visible.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Literal

from rocket_propulsion.core.constants import STANDARD_GRAVITY_M_S2
from rocket_propulsion.core.errors import DomainError
from rocket_propulsion.core.validation import records, strings, text

from .evidence import AnalysisResult, Citation, CitedScalar, derived


@dataclass(frozen=True, slots=True)
class CatalogField:
    """A scalar, interval or one-sided bound with its original interpretation.

    A maximum is NOT a nominal operating point. Different fields in the same
    table need not be simultaneously attainable.

    References
    ----------
    NASA-SOA26: section 4.5, survey envelopes are not technology limits or coupled maps.
    """

    name: str
    kind: Literal["nominal", "range", "upper_bound"]
    lower: CitedScalar | None
    upper: CitedScalar

    def __post_init__(self) -> None:
        """Keep interval bounds, source labels and units consistent.

        References
        ----------
        NASA-SE-6: section 6.6, preserve source-data semantics.
        """
        text(self.name, "catalog field")
        if not isinstance(self.kind, str) or self.kind not in {"nominal", "range", "upper_bound"}:
            raise DomainError("Unsupported catalogue field interpretation.")
        if not isinstance(self.upper, CitedScalar) or self.upper.status != "published":
            raise DomainError("Catalogue specifications must have published evidence.")
        if self.kind == "upper_bound":
            if self.lower is not None:
                raise DomainError("An upper bound does not invent a lower bound.")
        elif not isinstance(self.lower, CitedScalar):
            raise DomainError("A nominal value or interval needs its stated lower value.")
        else:
            if self.lower.status != "published" or self.lower.unit != self.upper.unit:
                raise DomainError("Catalogue endpoints must have matching published units.")
            if self.lower.value > self.upper.value:
                raise DomainError("Catalogue interval is reversed.")
            if self.kind == "nominal" and self.lower.value != self.upper.value:
                raise DomainError("Nominal catalogue fields contain one value.")


@dataclass(frozen=True, slots=True)
class ThrusterRecord:
    """Published hardware record; mass_scope prevents component/system confusion.

    References
    ----------
    NASA-SOA26: section 4.5, integrated systems versus thruster heads.
    NASA-SE-6: section 6.6, source-specific technical data.
    """

    key: str
    name: str
    family: str
    propellant: str
    mass_scope: str
    fields: tuple[CatalogField, ...]
    citations: tuple[Citation, ...]
    limitations: tuple[str, ...]

    def __post_init__(self) -> None:
        """Require references for textual catalog facts as well as each number.

        References
        ----------
        NASA-SE-6: section 6.6, source identification and metadata.
        """
        for label in ("key", "name", "family", "propellant", "mass_scope"):
            text(getattr(self, label), label)
        entries = records(self.fields, CatalogField, "catalog fields")
        if not entries or len({f.name for f in entries}) != len(entries):
            raise DomainError("Catalog fields must be non-empty and uniquely named.")
        citations = records(self.citations, Citation, "catalog text citations")
        if not citations:
            raise DomainError("Catalog textual information needs citations.")
        limitations = strings(self.limitations, "catalog limitations")
        if not limitations:
            raise DomainError("Catalog records must state their evidence limitations.")
        object.__setattr__(self, "fields", entries)
        object.__setattr__(self, "citations", citations)
        object.__setattr__(self, "limitations", limitations)

    def field(self, name: str) -> CatalogField:
        """Retrieve the sourced specification without replacing missing data by zero.

        References
        ----------
        NASA-SE-6: section 6.6, retrievable source data.
        """
        for item in self.fields:
            if item.name == name:
                return item
        raise DomainError(f"{self.key} has no published {name!r} field in this catalog.")


def _field(name: str, kind: str, low: float | None, high: float, unit: str,
           source_id: str, locator: str) -> CatalogField:
    """Encode a specific table row and its direct SI unit conversion.

    References
    ----------
    BIPM-SI: SI prefixes; NASA-SE-6: section 6.6, preserve the original row locator.
    """
    citations = (Citation(source_id, locator),)
    return CatalogField(name, kind,  # type: ignore[arg-type]
                        None if low is None else CitedScalar(low, unit, citations),
                        CitedScalar(high, unit, citations))


THRUSTERS = (
    ThrusterRecord(
        "moog-058e151", "Moog 058E151", "cold_gas", "Gaseous nitrogen",
        "thruster head and valve only; tank, regulator and feed hardware excluded",
        (
            _field("thrust", "nominal", 0.12, 0.12, "N", "MOOG-COLD",
                   "058E151: 120 mN at 6.9 bar nitrogen inlet"),
            _field("isp", "nominal", 65, 65, "s", "MOOG-COLD", "058E151: nominal Isp"),
            _field("inlet_pressure", "nominal", 690000, 690000, "Pa", "MOOG-COLD",
                   "058E151 nominal thrust condition: 6.9 bar absolute"),
            _field("dry_mass", "nominal", 0.07, 0.07, "kg", "MOOG-COLD",
                   "058E151 mass: 70 grams, component only"),
            _field("opening_power", "upper_bound", None, 10.5, "W", "MOOG-COLD",
                   "058E151: <10.5 W at 28 Vdc, 20 degC; strict bound represented as a ceiling"),
            _field("holding_power", "nominal", 1.5, 1.5, "W", "MOOG-COLD",
                   "058E151: 1.5 W holding at 10 Vdc"),
        ), (Citation("MOOG-COLD", "058E151 column, propellant and mass rows"),),
        ("Fixed nominal point assumes maintained inlet condition; not a blowdown map.",
         "Opening and holding powers are different electrical operating phases."),
    ),
    ThrusterRecord(
        "ariane-1n", "ArianeGroup 1 N hydrazine thruster", "monopropellant", "Hydrazine",
        "thruster with valve only; not complete propulsion-system dry mass",
        (
            _field("thrust", "nominal", 1, 1, "N", "ARIANE-1N", "Nominal thrust row"),
            _field("isp", "nominal", 220, 220, "s", "ARIANE-1N", "Nominal Isp row"),
            _field("reported_mass_flow", "nominal", 0.00044, 0.00044, "kg/s",
                   "ARIANE-1N", "Nominal mass flow: 0.44 g/s; retained independently of F/Isp"),
            _field("dry_mass", "nominal", 0.29, 0.29, "kg", "ARIANE-1N",
                   "Mass thruster with valve: 290 g"),
        ), (Citation("ARIANE-1N", "1N characteristics and propellant rows"),),
        ("Published rounded nominal thrust/Isp/flow do not close exactly; discrepancy is reported.",
         "No complete-system mass or firing/heater power has been inferred."),
    ),
    ThrusterRecord(
        "ariane-s10-13", "ArianeGroup S10-13, single-seat valve", "bipropellant", "MMH / NTO",
        "thruster with single-seat valve only; no tanks, feed or pressurant",
        (
            _field("thrust", "nominal", 10, 10, "N", "ARIANE-10N", "Nominal thrust row"),
            _field("isp", "nominal", 292, 292, "s", "ARIANE-10N", "Isp at nominal point"),
            _field("reported_mass_flow", "nominal", 0.0035, 0.0035, "kg/s",
                   "ARIANE-10N", "Nominal flow: 3.50 g/s"),
            _field("oxidizer_fuel_ratio", "range", 1.6, 1.65, "1", "ARIANE-10N",
                   "Nominal mixture-ratio interval; no preferred midpoint inferred"),
            _field("dry_mass", "nominal", 0.35, 0.35, "kg", "ARIANE-10N",
                   "Single-seat valve version: 350 g; dual-seat version is different"),
        ), (Citation("ARIANE-10N", "S10-13 identification, propellants, single-seat mass row"),),
        ("No inlet-pressure performance interpolation or full-system mass model.",),
    ),
    ThrusterRecord(
        "ecaps-1n-hpgp", "ECAPS 1 N HPGP", "green_monopropellant", "LMP-103S",
        "thruster only; NASA table lists 0.38 kg component mass",
        (
            _field("thrust", "range", 0.25, 1, "N", "NASA-SOA26", "ECAPS 1 N HPGP row"),
            _field("isp", "range", 204, 235, "s", "NASA-SOA26", "ECAPS 1 N HPGP row"),
            _field("dry_mass", "nominal", 0.38, 0.38, "kg", "NASA-SOA26",
                   "ECAPS 1 N HPGP row, thruster mass"),
            _field("power", "range", 8, 10, "W", "NASA-SOA26", "ECAPS 1 N HPGP power column"),
        ), (Citation("NASA-SOA26", "ECAPS 1 N HPGP device row, 2026 edition"),),
        ("Independent reported intervals; no paired thrust-Isp-power point is asserted.",),
    ),
    ThrusterRecord(
        "busek-bit3-page", "Busek BIT-3 (/bit3 page snapshot)", "gridded_ion", "Iodine",
        "manufacturer states integrated system dry mass with gimbal; no spacecraft EPS",
        (
            _field("thrust", "upper_bound", None, 0.0011, "N", "BUSEK-BIT3",
                   "Performance panel: up to 1.1 mN, not a paired operating point"),
            _field("isp", "upper_bound", None, 2150, "s", "BUSEK-BIT3",
                   "Performance panel: up to 2150 s, not a paired operating point"),
            _field("power", "range", 56, 75, "W", "BUSEK-BIT3", "System input power panel"),
            _field("dry_mass", "nominal", 1.4, 1.4, "kg", "BUSEK-BIT3",
                   "System dry mass: 1.40 kg with gimbal"),
        ), (Citation("BUSEK-BIT3", "BIT-3 performance panel"),
            Citation("BUSEK-RF", "Conflicting RF-family panel: 56-80 W, 1.25 mN, 2300 s")),
        ("Manufacturer pages disagree; this dated record preserves /bit3 values only.",
         "Upper bounds cannot be turned into a nominal propulsion point.",
         "No flight-readiness inference from old launch wording on the product page."),
    ),
)


def list_thrusters() -> tuple[ThrusterRecord, ...]:
    """Return the immutable, dated records with all missing-data distinctions intact.

    References
    ----------
    NASA-SE-6: section 6.6, provenance-preserving technical data.
    """
    return THRUSTERS


def thruster(key: str) -> ThrusterRecord:
    """Resolve a catalog record by its exact variant identifier.

    References
    ----------
    NASA-SE-6: sections 6.5/6.6, configuration identity and technical data.
    """
    for item in THRUSTERS:
        if item.key == key:
            return item
    raise DomainError(f"Unknown SASS thruster: {key!r}.")


def nominal_operating_point(key: str, *, require_validated: bool = False) -> AnalysisResult:
    """Use only published nominal thrust/Isp and expose rounded flow inconsistency.

    No point is constructed for ECAPS interval-only data or BIT-3 maxima.
    Derived flow closes F/(g0*Isp); any independently published flow remains a
    separate channel rather than silently changing thrust or Isp.

    References
    ----------
    NASA-ISP: mdot=F/(g0*Isp).
    MOOG-COLD, ARIANE-1N, ARIANE-10N: labelled nominal characteristics.
    NASA-SOA26: unpaired survey ranges cannot define a performance map.
    """
    item = thruster(key)
    if not isinstance(require_validated, bool):
        raise DomainError("require_validated must be boolean.")
    if require_validated:
        raise DomainError("No catalog record contains independently validated joint operating data; "
                          "nominal specifications cannot qualify motor and feed hardware.")
    force, isp = item.field("thrust"), item.field("isp")
    if force.kind != "nominal" or isp.kind != "nominal":
        raise DomainError("No paired nominal thrust/Isp is published in this catalog record.")
    flow = force.upper.value / (STANDARD_GRAVITY_M_S2 * isp.upper.value)
    methods = (Citation("NASA-ISP", "Isp=F/(mdot*g0); derive stream-closed study flow"),
               Citation("BIPM-SI", "Conventional standard gravity"))
    values = [("thrust", force.upper), ("isp", isp.upper),
              ("closed_mass_flow", derived(flow, "kg/s", (force.upper, isp.upper), methods))]
    if any(field.name == "reported_mass_flow" for field in item.fields):
        reported = item.field("reported_mass_flow").upper
        values.extend((
            ("reported_mass_flow", reported),
            ("reported_flow_relative_difference", derived(
                flow / reported.value - 1, "1", (force.upper, isp.upper, reported), methods,
                "Difference is retained; not interpreted as measurement uncertainty.",
            )),
        ))
    return AnalysisResult("manufacturer_nominal_constant_study_point", tuple(values), (
        *item.limitations,
        "Nominal characteristics are preliminary study evidence, not an engine qualification map.",
        "Joint operating conditions, measurement uncertainty and applicable feed-system tests are missing.",
    ))
