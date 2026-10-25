"""Mandatory field-level evidence and explicit design assumptions.

Every quantity used by the SASS tools carries registered citations with
claim-specific locators; analyst assumptions are labelled as such and derived
values carry both their method and their input sources.

References
----------
.. [1] NASA-SE-APP. NASA, "Systems Engineering Handbook: Appendix"; Appendix D and
       Trade Study Report definition.
.. [2] NASA-SE-6. NASA, "Systems Engineering Handbook: Crosscutting Technical
       Management"; sections 6.6 and 6.8; the Python representation is this
       library's implementation.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass
from typing import Literal

from rocket_propulsion.core.errors import DomainError
from rocket_propulsion.core.validation import real, records, strings, text

from .references import reference


@dataclass(frozen=True, slots=True)
class Citation:
    """A registered source together with a claim-specific locator.

    Attributes
    ----------
    source_id : str
        Key of an entry in :data:`REFERENCES`.
    locator : str
        Section, table, equation or row supporting the claim.

    References
    ----------
    .. [1] NASA-SE-APP. NASA, "Systems Engineering Handbook: Appendix"; Appendix D,
           source and unique identification.
    """

    source_id: str
    locator: str

    def __post_init__(self) -> None:
        """Reject unreviewed sources or locators missing the relevant claim.

        References
        ----------
        .. [1] NASA-SE-APP. NASA, "Systems Engineering Handbook: Appendix"; Appendix D,
               definitive source requirement.
        """
        reference(self.source_id)
        text(self.locator, "citation locator")


@dataclass(frozen=True, slots=True)
class CitedScalar:
    """An SI quantity with published evidence, a derivation or a labelled assumption.

    An ``'assumption'`` requires a rationale in ``note``; its citation documents the
    method of recording assumptions [1]_, not independent evidence for the chosen
    number. There is deliberately no extrapolated-data status.

    Attributes
    ----------
    value : float
        Finite value in ``unit``.
    unit : str
        Unit symbol.
    citations : tuple of Citation
        At least one source; published data must cite more than bookkeeping methods.
    status : {'published', 'assumption', 'derived'}, default 'published'
        Evidence kind.
    note : str
        Rationale (required for assumptions) or remark.
    assumption_notes : tuple of str
        Assumption lineage propagated from inputs.

    Raises
    ------
    DomainError
        If the value is not finite, a field is malformed, the status is unknown, or
        published data carries assumptions or only bookkeeping citations.

    References
    ----------
    .. [1] NASA-SE-6. NASA, "Systems Engineering Handbook: Crosscutting Technical
           Management"; section 6.8, preserve assumptions and decision-analysis inputs.
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
        .. [1] NASA-SE-6. NASA, "Systems Engineering Handbook: Crosscutting Technical
               Management"; sections 6.6 and 6.8, technical data and decision assumptions.
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

    Parameters
    ----------
    value : float
        Assumed value.
    unit : str
        Unit symbol.
    rationale : str
        Non-empty justification.

    Returns
    -------
    CitedScalar
        Quantity with status ``'assumption'``.

    Raises
    ------
    DomainError
        If the value or rationale is invalid.

    References
    ----------
    .. [1] NASA-SE-6. NASA, "Systems Engineering Handbook: Crosscutting Technical
           Management"; section 6.8, document assumptions in decision analysis.
    """
    return CitedScalar(value, unit, (Citation("NASA-SE-6", "6.8: analyst assumption"),),
                       "assumption", rationale)


def quantity(value: CitedScalar, unit: str, *, minimum: float | None = None,
             strictly_positive: bool = False) -> float:
    """Return the value of an evidence-bearing quantity after unit and domain checks.

    Parameters
    ----------
    value : CitedScalar
        Quantity to check.
    unit : str
        Required unit symbol.
    minimum : float, optional
        Inclusive lower limit.
    strictly_positive : bool, default False
        Require a value greater than zero.

    Returns
    -------
    float

    Raises
    ------
    DomainError
        If the value is not a ``CitedScalar`` in ``unit`` or violates the domain.

    References
    ----------
    .. [1] BIPM-SI. BIPM, "SI Brochure, 9th edition"; dimensional units.
    .. [2] NASA-SE-6. NASA, "Systems Engineering Handbook: Crosscutting Technical
           Management"; section 6.6, controlled technical data.
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
    """Create a computed quantity carrying its method and all original input sources.

    Parameters
    ----------
    value : float
        Computed value.
    unit : str
        Unit symbol.
    inputs : tuple of CitedScalar
        Quantities the value was computed from.
    methods : tuple of Citation
        At least one citation of the calculation method.
    note : str, default ''
        Remark.

    Returns
    -------
    CitedScalar
        Quantity with status ``'derived'`` and merged citations and assumption lineage.

    Raises
    ------
    DomainError
        If no method citation is given or an input is malformed.

    References
    ----------
    .. [1] NASA-SE-APP. NASA, "Systems Engineering Handbook: Appendix"; Trade Study
           Report, models and all data sources.
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
    """Named outputs with the declared model boundary and propagated evidence.

    Attributes
    ----------
    model : str
        Model identifier.
    values : tuple of (str, CitedScalar)
        Uniquely named outputs.
    limitations : tuple of str
        Non-empty statement of the model boundary.

    Raises
    ------
    DomainError
        If outputs are empty, duplicated or lack citations, or limitations are
        missing.

    References
    ----------
    .. [1] NASA-SE-APP. NASA, "Systems Engineering Handbook: Appendix"; Trade Study
           Report, model, results, uncertainty and sources.
    """

    model: str
    values: tuple[tuple[str, CitedScalar], ...]
    limitations: tuple[str, ...]

    def __post_init__(self) -> None:
        """Prevent duplicate output channels and citation-free result values.

        References
        ----------
        .. [1] NASA-SE-6. NASA, "Systems Engineering Handbook: Crosscutting Technical
               Management"; section 6.6, coherent technical data.
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

        Parameters
        ----------
        name : str
            Output name.

        Returns
        -------
        CitedScalar

        Raises
        ------
        DomainError
            If no output has that name.

        References
        ----------
        .. [1] NASA-SE-6. NASA, "Systems Engineering Handbook: Crosscutting Technical
               Management"; section 6.6, retrievable technical data.
        """
        for key, value in self.values:
            if key == name:
                return value
        raise DomainError(f"Unknown analysis output: {name!r}.")

    def to_dict(self) -> dict[str, object]:
        """Serialize values, limitations and full source records for independent review.

        Returns
        -------
        dict
            ``model``, ``values``, ``limitations`` and ``references`` (registry entries
            of every cited source, sorted by key).

        References
        ----------
        .. [1] NASA-SE-APP. NASA, "Systems Engineering Handbook: Appendix"; Trade Study
               Report, preserve computational results and sources.
        """
        source_ids = sorted({c.source_id for _, v in self.values for c in v.citations})
        return {
            "model": self.model,
            "values": {name: asdict(value) for name, value in self.values},
            "limitations": list(self.limitations),
            "references": [asdict(reference(key)) for key in source_ids],
        }
