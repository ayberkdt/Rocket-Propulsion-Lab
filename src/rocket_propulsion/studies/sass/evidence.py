"""Mandatory field-level evidence and explicit design assumptions.

References
----------
NASA-SE-APP: Appendix D and Trade Study Report definition.
NASA-SE-6: sections 6.6 and 6.8; the Python representation is our implementation.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass
from typing import Literal

from rocket_propulsion.core.errors import DomainError
from rocket_propulsion.core.validation import real, records, strings, text

from .references import reference


@dataclass(frozen=True, slots=True)
class Citation:
    """A registered source with a claim-specific locator.

    References
    ----------
    NASA-SE-APP: Appendix D, source and unique identification.
    """

    source_id: str
    locator: str

    def __post_init__(self) -> None:
        """Reject unreviewed sources or locators missing the relevant claim.

        References
        ----------
        NASA-SE-APP: Appendix D, definitive source requirement.
        """
        reference(self.source_id)
        text(self.locator, "citation locator")


@dataclass(frozen=True, slots=True)
class CitedScalar:
    """An SI quantity with evidence or an explicitly labelled analyst assumption.

    ``assumption`` requires a rationale; its citation documents the method for
    recording assumptions, not independent evidence for the chosen number.
    There is deliberately no extrapolated-data status.

    References
    ----------
    NASA-SE-6: section 6.8, preserve assumptions and decision-analysis inputs.
    """

    value: float
    unit: str
    citations: tuple[Citation, ...]
    status: Literal["published", "assumption", "derived"] = "published"
    note: str = ""
    assumption_notes: tuple[str, ...] = ()

    def __post_init__(self) -> None:
        """Require finite values, known sources, units and assumption rationale.

        References
        ----------
        NASA-SE-6: sections 6.6 and 6.8, technical data and decision assumptions.
        """
        object.__setattr__(self, "value", real(self.value, "cited value"))
        text(self.unit, "quantity unit")
        citations = records(self.citations, Citation, "quantity citations")
        if not citations:
            raise DomainError("Every SASS quantity requires at least one citation.")
        if not isinstance(self.status, str) or self.status not in {"published", "assumption", "derived"}:
            raise DomainError("Unsupported evidence status; extrapolation is disabled.")
        if not isinstance(self.note, str):
            raise DomainError("Quantity note must be text.")
        if self.status == "assumption":
            text(self.note, "analyst assumption rationale")
        notes = strings(self.assumption_notes, "assumption lineage")
        if self.status == "assumption":
            notes = tuple(dict.fromkeys((*notes, self.note)))
        if self.status == "published" and (notes or all(
            c.source_id in {"NASA-SE-6", "NASA-SE-APP"} for c in citations
        )):
            raise DomainError("Published data cannot be an assumption or cite only bookkeeping methods.")
        object.__setattr__(self, "assumption_notes", notes)
        object.__setattr__(self, "citations", citations)


def assumption(value: float, unit: str, rationale: str) -> CitedScalar:
    """Label an analyst input without presenting it as a published measurement.

    References
    ----------
    NASA-SE-6: section 6.8, document assumptions in decision analysis.
    """
    return CitedScalar(value, unit, (Citation("NASA-SE-6", "6.8: analyst assumption"),),
                       "assumption", rationale)


def quantity(value: CitedScalar, unit: str, *, minimum: float | None = None,
             strictly_positive: bool = False) -> float:
    """Check quantity units and the domain without stripping evidence silently.

    References
    ----------
    BIPM-SI: dimensional units; NASA-SE-6: section 6.6, controlled technical data.
    """
    if not isinstance(value, CitedScalar) or value.unit != unit:
        raise DomainError(f"Expected an evidence-bearing quantity in {unit}.")
    if strictly_positive and value.value <= 0:
        raise DomainError(f"Quantity in {unit} must be positive.")
    if minimum is not None and value.value < minimum:
        raise DomainError(f"Quantity in {unit} must be at least {minimum}.")
    return value.value


def derived(value: float, unit: str, inputs: tuple[CitedScalar, ...],
            methods: tuple[Citation, ...], note: str = "") -> CitedScalar:
    """Carry both method and original input sources into a computed quantity.

    References
    ----------
    NASA-SE-APP: Trade Study Report, models and all data sources.
    """
    input_records = records(inputs, CitedScalar, "derived inputs")
    citations = tuple(dict.fromkeys(
        (*records(methods, Citation, "calculation methods"),
         *(citation for item in input_records for citation in item.citations))
    ))
    if not methods:
        raise DomainError("A derived quantity needs a method citation.")
    assumptions = tuple(dict.fromkeys(n for item in input_records for n in item.assumption_notes))
    return CitedScalar(value, unit, citations, "derived", note, assumptions)


@dataclass(frozen=True, slots=True)
class AnalysisResult:
    """Named outputs plus the declared model boundary and propagated evidence.

    References
    ----------
    NASA-SE-APP: Trade Study Report, model, results, uncertainty and sources.
    """

    model: str
    values: tuple[tuple[str, CitedScalar], ...]
    limitations: tuple[str, ...]

    def __post_init__(self) -> None:
        """Prevent duplicate output channels and citation-free result values.

        References
        ----------
        NASA-SE-6: section 6.6, coherent technical data.
        """
        text(self.model, "analysis model")
        if not isinstance(self.values, (tuple, list)) or not self.values:
            raise DomainError("Analysis needs non-empty, uniquely named outputs.")
        values = []
        for pair in self.values:
            if not isinstance(pair, (tuple, list)) or len(pair) != 2:
                raise DomainError("Analysis outputs must be name/quantity pairs.")
            name, value = pair
            text(name, "output name")
            if not isinstance(value, CitedScalar):
                raise DomainError("Analysis outputs must retain citations.")
            values.append((name, value))
        if len({name for name, _ in values}) != len(values):
            raise DomainError("Analysis needs uniquely named outputs.")
        limitations = strings(self.limitations, "model limitations")
        if not limitations:
            raise DomainError("Analysis must state its model limitations.")
        object.__setattr__(self, "values", tuple(values))
        object.__setattr__(self, "limitations", limitations)

    def get(self, name: str) -> CitedScalar:
        """Return a named output together with its evidence.

        References
        ----------
        NASA-SE-6: section 6.6, retrievable technical data.
        """
        for key, value in self.values:
            if key == name:
                return value
        raise DomainError(f"Unknown analysis output: {name!r}.")

    def to_dict(self) -> dict[str, object]:
        """Serialize values and complete source details for independent review.

        References
        ----------
        NASA-SE-APP: Trade Study Report, preserve computational results and sources.
        """
        source_ids = sorted({c.source_id for _, v in self.values for c in v.citations})
        return {
            "model": self.model,
            "values": {name: asdict(value) for name, value in self.values},
            "limitations": list(self.limitations),
            "references": [asdict(reference(key)) for key in source_ids],
        }
