"""Printable, self-contained study report with citations beside every result.

References
----------
.. [1] NASA-SE-APP. NASA, "Systems Engineering Handbook: Appendix"; Trade Study
       Report, computational results, selection rule and data sources.
"""

from __future__ import annotations

import json
from html import escape

from .evidence import AnalysisResult, CitedScalar
from .project import ProjectStudy
from .references import REFERENCES


def _citations(value: CitedScalar) -> str:
    """Link every quantity to its source and retain the claim-specific locator.

    References
    ----------
    .. [1] NASA-SE-APP. NASA, "Systems Engineering Handbook: Appendix"; definitive
           evidence sources in a technical study.
    """
    return "; ".join(
        f'<a href="#ref-{escape(c.source_id)}" title="{escape(c.locator, quote=True)}">'
        f'{escape(c.source_id)}</a>' for c in value.citations
    )


def _analysis_table(result: AnalysisResult) -> str:
    """Render value, unit, assumption lineage and citations without recalculation.

    References
    ----------
    .. [1] NASA-SE-APP. NASA, "Systems Engineering Handbook: Appendix"; Trade Study
           Report, models, assumptions and computational results.
    """
    rows = []
    for name, q in result.values:
        rows.append(f"<tr><th>{escape(name)}</th><td>{q.value:.9g}</td>"
                    f"<td>{escape(q.unit)}</td><td>{_citations(q)}"
                    f"<br><small>{escape(q.note)}</small></td></tr>")
    boundary = "".join(f"<li>{escape(line)}</li>" for line in result.limitations)
    return ("<table><thead><tr><th>Metric</th><th>Value</th><th>Unit</th>"
            "<th>Evidence / assumptions</th></tr></thead><tbody>" + "".join(rows)
            + "</tbody></table><ul>" + boundary + "</ul>")


def build_project_html(study: ProjectStudy) -> str:
    """Produce a self-contained offline HTML report of a project study.

    Parameters
    ----------
    study : ProjectStudy
        Completed study.

    Returns
    -------
    str
        HTML document with the delta-v allocation, trade matrix, per-candidate
        budgets and screens, bibliography, inputs and an embedded JSON snapshot.

    Notes
    -----
    All values come from stored result objects; nothing is recalculated. The embedded
    JSON escapes ``<`` so input text cannot terminate the script element. The report
    keeps all alternatives, sources and the selection rule [1]_ and the decision
    rationale and uncertainties [2]_.

    References
    ----------
    .. [1] NASA-SE-APP. NASA, "Systems Engineering Handbook: Appendix"; Trade Study
           Report, all alternatives, sources and selection rule.
    .. [2] NASA-SE-6. NASA, "Systems Engineering Handbook: Crosscutting Technical
           Management"; 6.8.1.2.6, retain decision rationale and uncertainties.
    """
    document = study.to_dict()
    sections = ["<h2>Mission delta-v allocation</h2>", _analysis_table(study.budget)]
    sections.append("<h2>Conditional trade matrix</h2><p>Weights and anchors are analyst choices. "
                    "A numerical pass is preliminary analysis only. Unknown evidence remains "
                    "unranked. Scores use the original metrics, anchors and weights before final "
                    "rounding; rounded utility displays are not scoring inputs. "
                    "Ranking robustness under input and model uncertainty is not assessed. "
                    "No final course-project architecture is selected.</p>"
                    "<table><tr><th>Alternative</th><th>Eligibility</th><th>Score</th><th>Evidence</th></tr>")
    for row in study.scores:
        score = "unknown" if row.score is None else f"{row.score.value:.6g}"
        cites = "" if row.score is None else _citations(row.score)
        sections.append(f"<tr><th>{escape(row.name)}</th><td>{row.eligibility}</td>"
                        f"<td>{score}</td><td>{cites}</td></tr>")
    sections.append("</table>")
    for item, score in zip(study.candidates, study.scores, strict=True):
        sections.append(f"<h2>{escape(item.candidate.identifier)}</h2>")
        sections.append("<p>Catalog variant: " + escape(item.candidate.catalog_key) + "</p>")
        sections.append("<ul>" + "".join(f"<li>{escape(x)}</li>" for x in item.unresolved) + "</ul>")
        if item.operating_point is not None:
            sections.append("<h3>Published nominal study point</h3>" + _analysis_table(item.operating_point))
        if item.sizing is not None:
            sections.append("<h3>Mass and time budget</h3>" + _analysis_table(item.sizing.summary))
            sections.append("<details><summary>Numerical closure evidence</summary><p>"
                            "Closure residuals above use the reported binary values. Decimal "
                            "working residuals, precision agreement and denominator conditioning "
                            "are recorded below. These are numerical diagnostics; they do not "
                            "quantify input or physical model uncertainty.</p><pre>"
                            + escape(json.dumps(dict(item.sizing.numerical_evidence), indent=2))
                            + "</pre></details>")
            for name, burn in item.sizing.burns:
                sections.append(f"<details><summary>{escape(name)}</summary>"
                                + _analysis_table(burn) + "</details>")
        sections.append("<h3>Requirements screening</h3><table><tr><th>Requirement</th>"
                        "<th>Status</th><th>Margin</th><th>Evidence</th></tr>")
        for check in score.checks:
            margin = "unknown" if check.margin is None else f"{check.margin.value:.9g} {check.margin.unit}"
            evidence = _citations(check.requirement.limit)
            if check.margin is not None:
                evidence += "; " + _citations(check.margin)
            sections.append(f"<tr><th>{escape(check.requirement.identifier)}: "
                            f"{escape(check.requirement.description)}</th><td>{check.status}</td>"
                            f"<td>{escape(margin)}</td><td>{evidence}</td></tr>")
        sections.append("</table>")
    sections.append("<h2>References</h2><ol>")
    for source in REFERENCES:
        sections.append(f'<li id="ref-{escape(source.key)}"><strong>{escape(source.key)}</strong>: '
                        f'<a href="{escape(source.url, quote=True)}">{escape(source.title)}</a>. '
                        f'{escape(source.publisher)}. Accessed {source.accessed}. '
                        f'{escape(source.locator)}</li>')
    sections.append("</ol><details><summary>Complete reproducibility inputs</summary><pre>"
                    + escape(json.dumps(document["inputs"], indent=2, ensure_ascii=False))
                    + "</pre></details>")
    snapshot = json.dumps(document, sort_keys=True, allow_nan=False).replace("<", "\\u003c")
    style = """
    :root{font-family:system-ui,sans-serif;color:#13263a;background:#edf2f6}
    body{max-width:1200px;margin:2rem auto;padding:2rem;background:white;line-height:1.5}
    h1,h2,h3{line-height:1.2} h2{margin-top:2.2rem;color:#154d6e}
    .notice{padding:1rem;border-left:4px solid #b97818;background:#fff5df}
    table{border-collapse:collapse;width:100%;font-size:.88rem;margin:1rem 0}
    td,th{border:1px solid #cad5df;padding:.55rem;text-align:left;vertical-align:top}
    thead,th{background:#f1f5f8} a{color:#125f8c} small{color:#42576a}
    pre{white-space:pre-wrap;overflow-wrap:anywhere;font-size:.8rem}
    li{margin:.35rem 0;overflow-wrap:anywhere} details{margin:1rem 0}
    @media print{body{margin:0;padding:0;background:white;max-width:none}
    :root{background:white} h2,h3{break-after:avoid} tr{break-inside:avoid} a{color:inherit}}
    """
    return ("<!doctype html><html lang=en><meta charset=utf-8><meta name=viewport "
            "content='width=device-width,initial-scale=1'><title>" + escape(study.title)
            + "</title><style>" + style + "</style><body><h1>" + escape(study.title)
            + "</h1><p class=notice>UZM 508E propulsion preparation. Assignment scope is unconfirmed. "
            "The example inputs are explicitly labelled assumptions. Extrapolation is disabled. "
            "This report is preliminary engineering analysis. Motor and feed hardware readiness "
            "is NOT ASSESSED. No experimental validation, quantified result uncertainty or "
            "qualification evidence is provided by these calculations.</p><p>Snapshot hash: <code>"
            + str(document["study_hash"]) + "</code></p>" + "".join(sections)
            + '<script type="application/json" id="study-document">' + snapshot + "</script></body></html>")
