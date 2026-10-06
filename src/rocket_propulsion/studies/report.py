"""Dependency-free report packages for portable study review."""

from __future__ import annotations

import csv
import io
from dataclasses import dataclass
from html import escape
from math import isfinite

from .workspace import Workspace, workspace_to_json


@dataclass(frozen=True, slots=True)
class ReportBundle:
    """Machine-readable and printable representations of one workspace."""

    json: str
    csv: str
    svg: str
    html: str


def _numeric_rows(workspace: Workspace) -> list[tuple[str, str, float]]:
    rows: list[tuple[str, str, float]] = []
    for case in workspace.cases:
        for metric, raw_value in sorted(case.results.items()):
            if isinstance(raw_value, bool):
                continue
            try:
                value = float(raw_value)
            except (TypeError, ValueError):
                continue
            if isfinite(value):
                rows.append((case.case_id, metric, value))
    return rows


def _csv_report(workspace: Workspace) -> str:
    output = io.StringIO(newline="")
    writer = csv.writer(output, lineterminator="\n")
    writer.writerow(("case_id", "case_name", "model", "metric", "value"))
    names = {case.case_id: (case.name, case.model) for case in workspace.cases}
    for case_id, metric, value in _numeric_rows(workspace):
        name, model = names[case_id]
        writer.writerow((case_id, name, model, metric, f"{value:.12g}"))
    return output.getvalue()


def _svg_report(workspace: Workspace) -> str:
    rows = _numeric_rows(workspace)
    if not rows:
        return '<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 800 300"><text x="20" y="40">No scalar results</text></svg>\n'
    first_metric = rows[0][1]
    values = [(case_id, value) for case_id, metric, value in rows if metric == first_metric]
    maximum = max(abs(value) for _, value in values) or 1.0
    bars = []
    for index, (case_id, value) in enumerate(values):
        y = 45 + index * 45
        width = 600.0 * abs(value) / maximum
        bars.append(
            f'<text x="10" y="{y + 18}">{escape(case_id)}</text>'
            f'<rect x="150" y="{y}" width="{width:.3f}" height="24" fill="#0b615d"/>'
            f'<text x="{160 + width:.3f}" y="{y + 18}">{value:.6g}</text>'
        )
    height = 90 + len(values) * 45
    return (
        f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 800 {height}" role="img" '
        f'aria-label="{escape(first_metric)} case comparison">'
        f'<text x="10" y="25" font-weight="bold">{escape(first_metric)}</text>'
        + "".join(bars)
        + "</svg>\n"
    )


def _html_report(workspace: Workspace, svg: str) -> str:
    cases = []
    for case in workspace.cases:
        input_rows = "".join(
            f"<tr><th>{escape(str(key))}</th><td>{escape(str(value))}</td></tr>"
            for key, value in sorted(case.inputs.items())
        )
        result_rows = "".join(
            f"<tr><th>{escape(str(key))}</th><td>{escape(str(value))}</td></tr>"
            for key, value in sorted(case.results.items())
        )
        versions = ", ".join(
            f"{escape(str(key))}: {escape(str(value))}"
            for key, value in sorted(case.model_versions.items())
        )
        warnings = "".join(f"<li>{escape(item)}</li>" for item in case.warnings)
        cases.append(
            f"<section><h2>{escape(case.name)}</h2><p><b>Model:</b> {escape(case.model)}"
            f" · <b>Versions:</b> {versions or 'not recorded'}</p>"
            f"<h3>Inputs</h3><table>{input_rows}</table>"
            f"<h3>Results</h3><table>{result_rows}</table>"
            f"<h3>Warnings</h3><ul>{warnings or '<li>None</li>'}</ul></section>"
        )
    return (
        "<!doctype html><html lang=\"en\"><head><meta charset=\"utf-8\">"
        f"<title>{escape(workspace.name)}</title><style>"
        "body{font-family:system-ui;max-width:1100px;margin:auto;padding:2rem;color:#173332}"
        "table{border-collapse:collapse;width:100%}th,td{border:1px solid #ccd7d3;padding:.45rem;text-align:left}"
        "section{break-inside:avoid;margin:2rem 0}@media print{body{padding:0}}</style></head><body>"
        f"<h1>{escape(workspace.name)}</h1><p>Schema {workspace.schema_version}</p>{svg}"
        + "".join(cases)
        + "</body></html>\n"
    )


def build_report_bundle(workspace: Workspace) -> ReportBundle:
    """Build JSON, tidy CSV, comparison SVG, and a printable HTML summary."""

    svg = _svg_report(workspace)
    return ReportBundle(
        json=workspace_to_json(workspace),
        csv=_csv_report(workspace),
        svg=svg,
        html=_html_report(workspace, svg),
    )
