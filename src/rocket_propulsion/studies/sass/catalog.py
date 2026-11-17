"""Small, field-cited spacecraft thruster catalog; no envelope-to-point inference.

References
----------
.. [1] MOOG-COLD. Moog, "Cold Gas Thrusters: Performance Characteristics";
       manufacturer nominal characteristics.
.. [2] ARIANE-1N. ArianeGroup, "1 N Monopropellant Hydrazine Thruster"; manufacturer
       nominal characteristics.
.. [3] ARIANE-10N. ArianeGroup, "10 N Bipropellant Thruster"; manufacturer nominal
       characteristics.
.. [4] NASA-SOA26. NASA, "State-of-the-Art of Small Spacecraft Technology: In-Space
       Propulsion"; ECAPS row, ranges only.
.. [5] BUSEK-BIT3. Busek, "BIT-3 RF Ion Thruster product page"; unresolved page
       discrepancy.
.. [6] BUSEK-RF. Busek, "RF Ion Thrusters product family page"; unresolved page
       discrepancy.
.. [7] NASA-SE-6. NASA, "Systems Engineering Handbook: Crosscutting Technical
       Management"; keep data provenance and system/component boundaries visible.
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from typing import Literal

from rocket_propulsion.core.constants import STANDARD_GRAVITY_M_S2
from rocket_propulsion.core.errors import DomainError
from rocket_propulsion.core.validation import real, records, strings, text

from .evidence import AnalysisResult, Citation, CitedScalar, derived, quantity


@dataclass(frozen=True, slots=True)
class CatalogField:
    """A published scalar, interval or one-sided bound with its original interpretation.

    A maximum is NOT a nominal operating point, and different fields of the same
    table need not be simultaneously attainable [1]_.

    Attributes
    ----------
    name : str
        Field name, for example ``'thrust'`` or ``'isp'``.
    kind : {'nominal', 'range', 'upper_bound'}
        Interpretation of the published value.
    lower : CitedScalar or None
        Published lower value; equal to ``upper`` for a nominal field and None for an
        upper bound.
    upper : CitedScalar
        Published upper or nominal value.
    nominal_point_id : str or None
        Identifier of the source's labelled nominal study point. Absent for
        fields without evidence of that pairing.
    reported_psia : CitedScalar or None
        Alternate printed value for the same absolute inlet pressure ceiling.

    Raises
    ------
    DomainError
        If the kind is unknown, a bound is not published evidence, units differ, an
        interval is reversed, or a nominal field has two values.

    References
    ----------
    .. [1] NASA-SOA26. NASA, "State-of-the-Art of Small Spacecraft Technology: In-Space
           Propulsion"; section 4.5, survey envelopes are not technology limits or
           coupled maps.
    """

    name: str
    kind: Literal["nominal", "range", "upper_bound"]
    lower: CitedScalar | None
    upper: CitedScalar
    nominal_point_id: str | None = None
    pressure_reference: Literal["absolute", "unspecified"] | None = None
    pressure_port: Literal["inlet", "fuel_inlet", "oxidizer_inlet"] | None = None
    reported_psia: CitedScalar | None = None

    def __post_init__(self) -> None:
        """Keep interval bounds, source labels and units consistent.

        References
        ----------
        .. [1] NASA-SE-6. NASA, "Systems Engineering Handbook: Crosscutting Technical
               Management"; section 6.6, preserve source-data semantics.
        """
        text(self.name, "catalog field")
        if not isinstance(self.kind, str) or self.kind not in {"nominal", "range", "upper_bound"}:
            raise DomainError("Unsupported catalogue field interpretation.")
        if self.nominal_point_id is not None:
            text(self.nominal_point_id, "catalogue nominal point ID")
            if self.kind != "nominal":
                raise DomainError("Only nominal fields can belong to a labelled nominal study point.")
        if not isinstance(self.upper, CitedScalar) or self.upper.status != "published":
            raise DomainError("Catalogue specifications must have published evidence.")
        if self.pressure_reference not in (None, "absolute", "unspecified"):
            raise DomainError("Unsupported catalogue pressure reference.")
        if self.pressure_port not in (None, "inlet", "fuel_inlet", "oxidizer_inlet"):
            raise DomainError("Unsupported catalogue pressure port.")
        if self.pressure_port is not None or self.pressure_reference is not None:
            quantity(self.upper, "Pa", minimum=0)
            if self.pressure_port is not None and self.pressure_reference is None:
                raise DomainError("A catalogue pressure port requires explicit reference metadata.")
        if self.reported_psia is not None:
            quantity(self.reported_psia, "psia", strictly_positive=True)
            if (self.kind != "upper_bound" or self.pressure_reference != "absolute"
                    or self.pressure_port is None or self.reported_psia.status != "published"
                    or not set(self.upper.citations).intersection(self.reported_psia.citations)):
                raise DomainError("The alternate psia value requires a published absolute inlet ceiling from the same row.")
            quantity(self.upper, "Pa", strictly_positive=True)
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
    """Published hardware record with an explicit mass scope.

    ``mass_scope`` states whether the dry mass covers a thruster head, a thruster with
    valve or an integrated system, so component and system masses are not confused
    [1]_. Every textual fact and every number carries its source [2]_.

    Attributes
    ----------
    key : str
        Exact variant identifier.
    name, family, propellant : str
        Display name, technology family and propellant.
    mass_scope : str
        Hardware covered by the published dry mass.
    fields : tuple of CatalogField
        Uniquely named published specifications.
    citations : tuple of Citation
        Sources of the textual information.
    limitations : tuple of str
        Evidence limitations of the record.

    Raises
    ------
    DomainError
        If a text field is empty, fields are missing or duplicated, or citations or
        limitations are missing.

    References
    ----------
    .. [1] NASA-SOA26. NASA, "State-of-the-Art of Small Spacecraft Technology: In-Space
           Propulsion"; section 4.5, integrated systems versus thruster heads.
    .. [2] NASA-SE-6. NASA, "Systems Engineering Handbook: Crosscutting Technical
           Management"; section 6.6, source-specific technical data.
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
        .. [1] NASA-SE-6. NASA, "Systems Engineering Handbook: Crosscutting Technical
               Management"; section 6.6, source identification and metadata.
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
        """Return a published specification without replacing missing data by zero.

        Parameters
        ----------
        name : str
            Field name.

        Returns
        -------
        CatalogField

        Raises
        ------
        DomainError
            If the record has no published field of that name.

        References
        ----------
        .. [1] NASA-SE-6. NASA, "Systems Engineering Handbook: Crosscutting Technical
               Management"; section 6.6, retrievable source data.
        """
        for item in self.fields:
            if item.name == name:
                return item
        raise DomainError(f"{self.key} has no published {name!r} field in this catalog.")

    def nominal_pair(self) -> tuple[CitedScalar, CitedScalar] | None:
        """Return thrust and Isp when their source labels a shared nominal point.

        References
        ----------
        .. [1] NASA-SE-6. NASA, "Systems Engineering Handbook"; section 6.6,
               controlled technical data and explicit conditions.
        """
        fields = {field.name: field for field in self.fields}
        force, isp = fields.get("thrust"), fields.get("isp")
        if (force is None or isp is None or force.nominal_point_id is None
                or force.nominal_point_id != isp.nominal_point_id):
            return None
        if not {citation.source_id for citation in force.upper.citations}.intersection(
            citation.source_id for citation in isp.upper.citations
        ):
            return None
        return force.upper, isp.upper


def _field(name: str, kind: str, low: float | None, high: float, unit: str,
           source_id: str, locator: str, *, nominal_point_id: str | None = None,
           pressure_reference: Literal["absolute", "unspecified"] | None = None,
           pressure_port: Literal["inlet", "fuel_inlet", "oxidizer_inlet"] | None = None,
           reported_psia: float | None = None) -> CatalogField:
    """Encode a specific table row and its direct SI unit conversion.

    References
    ----------
    .. [1] BIPM-SI. BIPM, "SI Brochure, 9th edition"; SI prefixes.
    .. [2] NASA-SE-6. NASA, "Systems Engineering Handbook: Crosscutting Technical
           Management"; section 6.6, preserve the original row locator.
    """
    citations = (Citation(source_id, locator),)
    return CatalogField(name, kind,  # type: ignore[arg-type]
                        None if low is None else CitedScalar(low, unit, citations),
                        CitedScalar(high, unit, citations), nominal_point_id, pressure_reference, pressure_port,
                        None if reported_psia is None else CitedScalar(reported_psia, "psia", citations))


THRUSTERS = (
    ThrusterRecord(
        "moog-058e151", "Moog 058E151", "cold_gas", "Gaseous nitrogen",
        "thruster head and valve only; tank, regulator and feed hardware excluded",
        (
            _field("thrust", "nominal", 0.12, 0.12, "N", "MOOG-COLD",
                   "058E151: 120 mN at 6.9 bar nitrogen inlet", nominal_point_id="058E151-GN2-nominal"),
            _field("isp", "nominal", 65, 65, "s", "MOOG-COLD", "058E151: nominal Isp",
                   nominal_point_id="058E151-GN2-nominal"),
            _field("inlet_pressure", "nominal", 690000, 690000, "Pa", "MOOG-COLD",
                    "058E151 nominal thrust condition: 6.9 bar absolute",
                   nominal_point_id="058E151-GN2-nominal", pressure_reference="absolute", pressure_port="inlet"),
            _field("maximum_expected_operating_pressure", "upper_bound", None, 2_760_000, "Pa", "MOOG-COLD",
                   "058E151 column: MEOP 27.6 bar (400 psia)",
                   pressure_reference="absolute", pressure_port="inlet", reported_psia=400),
            _field("proof_pressure", "nominal", 7_000_000, 7_000_000, "Pa", "MOOG-COLD",
                   "058E151 column: proof 70.0 bar (1015 psia)",
                   pressure_reference="absolute", pressure_port="inlet"),
            _field("burst_pressure", "nominal", 11_140_000, 11_140_000, "Pa", "MOOG-COLD",
                   "058E151 column: burst 111.4 bar (1615 psia)",
                   pressure_reference="absolute", pressure_port="inlet"),
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
            _field("thrust", "nominal", 1, 1, "N", "ARIANE-1N", "Nominal thrust row",
                   nominal_point_id="ArianeGroup-1N-nominal"),
            _field("isp", "nominal", 220, 220, "s", "ARIANE-1N", "Nominal Isp row",
                   nominal_point_id="ArianeGroup-1N-nominal"),
            _field("reported_mass_flow", "nominal", 0.00044, 0.00044, "kg/s",
                   "ARIANE-1N", "Nominal mass flow: 0.44 g/s; retained independently of F/Isp"),
            _field("reported_mass_flow_range", "range", 0.000142, 0.000447, "kg/s",
                   "ARIANE-1N", "Mass Flow Rate, Range: 0.142 to 0.447 g/s"),
            _field("thrust_range", "range", 0.32, 1.1, "N", "ARIANE-1N",
                   "Thrust Range: 0.320 to 1.1 N; separate from the nominal row"),
            _field("inlet_pressure_range", "range", 550000, 2200000, "Pa", "ARIANE-1N",
                    "Inlet Pressure Range: 5.5 to 22 bar. Pressure reference unspecified.",
                   pressure_reference="unspecified", pressure_port="inlet"),
            _field("dry_mass", "nominal", 0.29, 0.29, "kg", "ARIANE-1N",
                   "Mass thruster with valve: 290 g"),
            _field("single_burn_life", "upper_bound", None, 12 * 3600, "s", "ARIANE-1N",
                   "Single Burn Life: 12 h per thruster"),
            _field("accumulated_burn_life", "upper_bound", None, 50 * 3600, "s", "ARIANE-1N",
                   "Accumulated Burn Life: 50 h per thruster"),
            _field("propellant_throughput", "upper_bound", None, 67, "kg", "ARIANE-1N",
                   "Propellant Throughput: 67 kg per thruster"),
            _field("total_impulse_life", "upper_bound", None, 135_000, "N*s", "ARIANE-1N",
                   "Total Impulse: 135,000 N s per thruster"),
            _field("cycle_life", "upper_bound", None, 59_000, "1", "ARIANE-1N",
                   "Cycle Life: 59,000 cycles per thruster"),
        ), (Citation("ARIANE-1N", "1N characteristics and propellant rows"),),
        ("Published nominal thrust, Isp and mass flow do not close; rounding alone is unconfirmed.",
         "No complete-system mass or firing/heater power has been inferred."),
    ),
    ThrusterRecord(
        "ariane-s10-13", "ArianeGroup S10-13, single-seat valve", "bipropellant", "MMH / NTO",
        "thruster with single-seat valve only; no tanks, feed or pressurant",
        (
            _field("thrust", "nominal", 10, 10, "N", "ARIANE-10N", "Nominal thrust row",
                   nominal_point_id="ArianeGroup-10N-nominal"),
            _field("isp", "nominal", 292, 292, "s", "ARIANE-10N", "Isp at nominal point",
                   nominal_point_id="ArianeGroup-10N-nominal"),
            _field("reported_mass_flow", "nominal", 0.0035, 0.0035, "kg/s",
                   "ARIANE-10N", "Nominal flow: 3.50 g/s"),
            _field("reported_mass_flow_range", "range", 0.0023, 0.0042, "kg/s",
                   "ARIANE-10N", "Flow Rate, Range: 2.30 to 4.20 g/s; 10 N family table"),
            _field("thrust_range", "range", 6, 12.5, "N", "ARIANE-10N",
                   "Thrust Range: 6.0 to 12.5 N; 10 N family table"),
            _field("inlet_pressure_range", "range", 1000000, 2300000, "Pa", "ARIANE-10N",
                    "Inlet Pressure Range: 10 to 23 bar. Port and pressure reference unspecified.",
                   pressure_reference="unspecified"),
            _field("oxidizer_fuel_ratio", "range", 1.6, 1.65, "1", "ARIANE-10N",
                   "Nominal mixture-ratio interval; no preferred midpoint inferred"),
            _field("dry_mass", "nominal", 0.35, 0.35, "kg", "ARIANE-10N",
                   "Single-seat valve version: 350 g; dual-seat version is different"),
            _field("single_burn_life", "upper_bound", None, 8 * 3600, "s", "ARIANE-10N",
                   "10 N family: qualified longest single burn, 8 h"),
            _field("accumulated_burn_life", "upper_bound", None, 69 * 3600, "s", "ARIANE-10N",
                   "10 N family: qualified accumulated burn life, 69 h"),
            _field("cycle_life", "upper_bound", None, 1_000_000, "1", "ARIANE-10N",
                   "10 N family: qualified cycle life, 1,000,000 cycles"),
        ), (Citation("ARIANE-10N", "S10-13 identification, propellants, single-seat mass row"),),
        ("No inlet-pressure performance interpolation or full-system mass model.",
         "Published family life figures require confirmation for the selected valve and duty profile."),
    ),
    ThrusterRecord(
        "ecaps-1n-hpgp", "ECAPS 1 N HPGP", "green_monopropellant", "LMP-103S",
        "thruster only; NASA table lists 0.38 kg component mass",
        (
            _field("thrust", "range", 0.25, 1, "N", "NASA-SOA26", "ECAPS 1 N HPGP row"),
            _field("isp", "range", 204, 235, "s", "NASA-SOA26", "ECAPS 1 N HPGP row"),
            _field("dry_mass", "nominal", 0.38, 0.38, "kg", "NASA-SOA26",
                   "ECAPS 1 N HPGP row, thruster mass"),
            _field("survey_power", "range", 8, 10, "W", "NASA-SOA26", "ECAPS 1 N HPGP power column"),
        ), (Citation("NASA-SOA26", "ECAPS 1 N HPGP device row, 2026 edition"),),
        ("Survey power has no phase or bus scope in the device row.",
         "Independent survey intervals cannot establish a paired operating point."),
    ),
    ThrusterRecord(
        "ecaps-1n-hpgp-bradford", "Bradford ECAPS 1 N HPGP (manufacturer sheet)",
        "green_monopropellant", "LMP-103S", "published thruster mass; complete-system mass unspecified",
        (
            _field("thrust", "range", 0.25, 1, "N", "ECAPS-1N-BRADFORD", "Specifications: thrust range"),
            _field("isp", "range", 204, 231, "s", "ECAPS-1N-BRADFORD",
                   "Specifications: typical steady-state vacuum specific impulse"),
            _field("dry_mass", "nominal", 0.38, 0.38, "kg", "ECAPS-1N-BRADFORD", "Specifications: mass"),
            _field("reactor_preheat_power", "range", 8, 10, "W", "ECAPS-1N-BRADFORD",
                   "Specifications: regulated reactor preheating power"),
        ), (Citation("ECAPS-1N-BRADFORD", "Specifications: 1 N HPGP hardware and propellant"),),
        ("Preheat power is a reactor load. Complete firing bus power is unspecified.",
         "The typical steady-state Isp range differs from the NASA survey and lacks a paired map."),
    ),
    ThrusterRecord(
        "busek-bit3-page", "Busek BIT-3 (/bit3 page snapshot)", "gridded_ion", "Iodine",
        "manufacturer states integrated system dry mass with gimbal; no spacecraft EPS",
        (
            _field("thrust", "upper_bound", None, 0.0011, "N", "BUSEK-BIT3-PDF",
                   "Performance panel: up to 1.1 mN, not a paired operating point"),
            _field("isp", "upper_bound", None, 2150, "s", "BUSEK-BIT3-PDF",
                   "Performance panel: up to 2150 s, not a paired operating point"),
            _field("power", "range", 56, 75, "W", "BUSEK-BIT3-PDF", "Standard Specifications: system input power"),
            _field("dry_mass", "nominal", 1.4, 1.4, "kg", "BUSEK-BIT3-PDF",
                   "System dry mass: 1.40 kg with gimbal"),
        ), (Citation("BUSEK-BIT3", "BIT-3 performance panel"),
            Citation("BUSEK-RF", "Conflicting RF-family panel: 56-80 W, 1.25 mN, 2300 s")),
        ("Manufacturer pages disagree; this dated record preserves /bit3 values only.",
         "Upper bounds cannot be turned into a nominal propulsion point.",
         "The August 2021 version 1.0 sheet confirms these values for the gimbaled mass scope.",
         "The older launch wording does not establish current flight-readiness."),
    ),
)


def list_thrusters() -> tuple[ThrusterRecord, ...]:
    """Return the immutable, dated catalogue records.

    Returns
    -------
    tuple of ThrusterRecord
        All records, with nominal, range and upper-bound distinctions intact [1]_.

    References
    ----------
    .. [1] NASA-SE-6. NASA, "Systems Engineering Handbook: Crosscutting Technical
           Management"; section 6.6, provenance-preserving technical data.
    """
    return THRUSTERS


def thruster(key: str) -> ThrusterRecord:
    """Resolve a catalogue record by its exact variant identifier.

    Parameters
    ----------
    key : str
        Variant identifier, for example ``'moog-058e151'``.

    Returns
    -------
    ThrusterRecord

    Raises
    ------
    DomainError
        If no record has that identifier.

    References
    ----------
    .. [1] NASA-SE-6. NASA, "Systems Engineering Handbook: Crosscutting Technical
           Management"; sections 6.5/6.6, configuration identity and technical data.
    """
    for item in THRUSTERS:
        if item.key == key:
            return item
    raise DomainError(f"Unknown SASS thruster: {key!r}.")


def nominal_operating_point(key: str, *, require_validated: bool = False) -> AnalysisResult:
    """Build a constant study point from published nominal thrust and Isp only.

    Parameters
    ----------
    key : str
        Catalogue variant identifier.
    require_validated : bool, default False
        Request independently validated joint operating data. No record contains
        such data, so True always raises.

    Returns
    -------
    AnalysisResult
        ``thrust``, ``isp`` and ``closed_mass_flow``; when a mass flow is published,
        also ``reported_mass_flow`` and ``reported_flow_relative_difference``.

    Raises
    ------
    DomainError
        If the key is unknown, ``require_validated`` is not a bool or is True, or the
        record has no paired nominal thrust and Isp.

    Notes
    -----
    The closed flow is mdot = F / (g0 * Isp) [1]_ from the labelled nominal
    characteristics [2]_ [3]_ [4]_. An independently published flow remains a separate
    channel and its relative difference is reported, not interpreted as measurement
    uncertainty. No point is constructed from interval-only survey data or from
    maxima [5]_.

    References
    ----------
    .. [1] NASA-ISP. NASA Glenn Research Center, "Specific Impulse"; mdot=F/(g0*Isp).
    .. [2] MOOG-COLD. Moog, "Cold Gas Thrusters: Performance Characteristics"; labelled
           nominal characteristics.
    .. [3] ARIANE-1N. ArianeGroup, "1 N Monopropellant Hydrazine Thruster"; labelled
           nominal characteristics.
    .. [4] ARIANE-10N. ArianeGroup, "10 N Bipropellant Thruster"; labelled nominal
           characteristics.
    .. [5] NASA-SOA26. NASA, "State-of-the-Art of Small Spacecraft Technology: In-Space
           Propulsion"; unpaired survey ranges cannot define a performance map.
    """
    if not isinstance(require_validated, bool):
        raise DomainError("require_validated must be boolean.")
    item = thruster(key)
    if require_validated:
        raise DomainError("No catalog record contains independently validated joint operating data; "
                          "nominal specifications cannot qualify motor and feed hardware.")
    operating = _nominal_study_point(item)
    if operating is None:
        raise DomainError("No paired nominal thrust/Isp is published in this catalog record.")
    return operating


@dataclass(frozen=True, slots=True)
class CatalogueFlowCheck:
    """Comparison of a calculated nominal flow against a published component reference.

    References
    ----------
    .. [1] NASA-ISP, specific impulse and mass flow identity.
    .. [2] NASA-SE-6, 6.6, scope and pedigree of manufacturer data.
    """

    status: Literal["below_published_range", "above_published_range", "within_published_range", "unassessed"]
    nominal_flow: CitedScalar | None
    reported_range: CatalogField
    basis: str


def _closed_flow(force_value: float, isp_value: float) -> Fraction:
    """Compute the nominal flow without rounding intermediate products.

    References
    ----------
    .. [1] NASA-ISP, specific impulse identity.
    """
    return Fraction(force_value) / (Fraction(STANDARD_GRAVITY_M_S2) * Fraction(isp_value))


def _catalogue_flow_check(
    record: ThrusterRecord, operating: AnalysisResult | None,
) -> CatalogueFlowCheck | None:
    """Screen nominal flow against a comparable source's published flow range.

    A discrepancy triggers source review. Published ranges alone do not
    establish a coupled performance map or a qualification envelope.

    References
    ----------
    .. [1] NASA-SE-6, 6.6, comparison scope and source traceability.
    .. [2] ARIANE-1N, nominal values and independently listed mass flow range.
    .. [3] ARIANE-10N, family nominal values and separately listed mass flow range.
    """
    bounds = next((field for field in record.fields if field.name == "reported_mass_flow_range"), None)
    if bounds is None:
        return None
    if operating is None:
        return CatalogueFlowCheck("unassessed", None, bounds, "No paired nominal thrust and Isp.")
    nominal = operating.get("closed_mass_flow")
    if bounds.kind != "range" or bounds.lower is None:
        return CatalogueFlowCheck("unassessed", nominal, bounds, "The source did not provide two range endpoints.")
    force, isp = operating.get("thrust"), operating.get("isp")
    cited = [{citation.source_id for citation in value.citations}
             for value in (force, isp, bounds.lower, bounds.upper)]
    if not set.intersection(*cited):
        return CatalogueFlowCheck("unassessed", nominal, bounds, "No common source for both endpoints and the nominal pair.")
    lo = quantity(bounds.lower, "kg/s", strictly_positive=True)
    hi = quantity(bounds.upper, "kg/s", strictly_positive=True)
    flow = _closed_flow(quantity(force, "N", strictly_positive=True),
                        quantity(isp, "s", strictly_positive=True))
    status: Literal["below_published_range", "above_published_range", "within_published_range"] = (
        "below_published_range" if flow < Fraction(lo) else
        "above_published_range" if flow > Fraction(hi) else "within_published_range"
    )
    return CatalogueFlowCheck(status, nominal, bounds,
                              "Same-source published range screened against thrust divided by standard gravity and Isp.")


def _nominal_study_point(item: ThrusterRecord) -> AnalysisResult | None:
    """Calculate a preliminary point from the already resolved catalogue record.

    References
    ----------
    .. [1] NASA-ISP. NASA Glenn Research Center, "Specific Impulse"; mass-flow identity.
    .. [2] NASA-SE-6. NASA, "Systems Engineering Handbook"; technical-data provenance.
    """
    pair = item.nominal_pair()
    if pair is None:
        return None
    force, isp = pair
    force_value = quantity(force, "N", strictly_positive=True)
    isp_value = quantity(isp, "s", strictly_positive=True)
    flow = _closed_flow(force_value, isp_value)
    methods = (Citation("NASA-ISP", "Isp=F/(mdot*g0); derive stream-closed study flow"),
               Citation("BIPM-SI", "Conventional standard gravity"))
    values = [("thrust", force), ("isp", isp),
              ("closed_mass_flow", derived(real(flow, "closed nominal flow"), "kg/s", pair, methods))]
    reported_field = next((field for field in item.fields if field.name == "reported_mass_flow"), None)
    if reported_field is not None and reported_field.kind == "nominal":
        reported = reported_field.upper
        reported_flow = quantity(reported, "kg/s", strictly_positive=True)
        values.extend((
            ("reported_mass_flow", reported),
            ("reported_flow_relative_difference", derived(
                real(flow / Fraction(reported_flow) - 1, "reported flow relative difference"),
                "1", (*pair, reported), methods,
                "Difference is retained; not interpreted as measurement uncertainty.",
            )),
        ))
    return AnalysisResult("manufacturer_nominal_constant_study_point", tuple(values), (
        *item.limitations,
        "Nominal characteristics are preliminary study evidence, not an engine qualification map.",
        "Joint operating conditions, measurement uncertainty and applicable feed-system tests are missing.",
    ))
