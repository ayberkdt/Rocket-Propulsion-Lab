"""Dependency-free HTML reporting for standalone burn results."""

from __future__ import annotations

from html import escape

from .models import BurnResult, ProfiledBurnResult


def build_burn_html_report(result: BurnResult | ProfiledBurnResult) -> str:
    """Render a self-contained engineering summary without recomputing physics."""

    summary = result.summary
    status = "Target reached" if summary.target_achieved else "Constraint limited"
    warning_items = "".join(f"<li>{escape(item)}</li>" for item in result.warnings)
    if not warning_items:
        warning_items = "<li>None</li>"
    tank_rows = "".join(
        "<tr>"
        f"<td>{escape(tank_id)}</td>"
        f"<td>{consumed:.9g}</td>"
        "</tr>"
        for tank_id, consumed in summary.consumed_by_tank_kg
    )
    schedule_section = ""
    if isinstance(result, ProfiledBurnResult):
        schedule_rows = "".join(
            "<tr>"
            f"<td>{index}</td>"
            f"<td>{escape(segment.phase)}</td>"
            f"<td>{segment.duration_s:.9g}</td>"
            f"<td>{segment.start_throttle:.9g}</td>"
            f"<td>{segment.end_throttle:.9g}</td>"
            "</tr>"
            for index, segment in enumerate(result.schedule.segments, start=1)
        )
        schedule_section = (
            "<h2>Throttle schedule</h2>"
            "<table><thead><tr><th>#</th><th>Phase</th><th>Duration (s)</th>"
            "<th>Start</th><th>End</th></tr></thead>"
            f"<tbody>{schedule_rows}</tbody></table>"
        )
    return f"""<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width, initial-scale=1">
  <title>{escape(result.definition.name)} — burn report</title>
  <style>
    body {{ font: 16px/1.45 system-ui, sans-serif; max-width: 900px; margin: 2rem auto; padding: 0 1rem; color: #18202a; }}
    table {{ border-collapse: collapse; width: 100%; margin: 1rem 0; }}
    th, td {{ border-bottom: 1px solid #d7dde5; padding: .45rem; text-align: left; }}
    code {{ overflow-wrap: anywhere; }}
    .status {{ font-weight: 700; }}
    .note {{ background: #f3f6f9; padding: .8rem; border-left: 4px solid #607d9b; }}
  </style>
</head>
<body>
  <h1>{escape(result.definition.name)}</h1>
  <p class="status">{status}: {escape(summary.cutoff_reason.value)}</p>
  <table>
    <tbody>
      <tr><th>Duration</th><td>{summary.duration_s:.9g} s</td></tr>
      <tr><th>Initial / final mass</th><td>{summary.initial_mass_kg:.9g} / {summary.final_mass_kg:.9g} kg</td></tr>
      <tr><th>Consumed propellant</th><td>{summary.consumed_propellant_kg:.9g} kg</td></tr>
      <tr><th>Delivered impulse</th><td>{summary.delivered_total_impulse_n_s:.9g} N s</td></tr>
      <tr><th>System equivalent Isp</th><td>{summary.system_equivalent_specific_impulse_s:.9g} s</td></tr>
      <tr><th>Ideal 1-D delta-v</th><td>{summary.ideal_delta_v_m_s:.9g} m/s</td></tr>
      <tr><th>Binding constraint</th><td>{escape(summary.binding_constraint or "None")}</td></tr>
    </tbody>
  </table>
  {schedule_section}
  <h2>Tank consumption</h2>
  <table><thead><tr><th>Tank</th><th>Consumed (kg)</th></tr></thead><tbody>{tank_rows}</tbody></table>
  <h2>Closure evidence</h2>
  <table>
    <tbody>
      <tr><th>Method</th><td>{escape(result.evidence.method)}</td></tr>
      <tr><th>Mass residual</th><td>{summary.mass_closure_error_kg:.6e} kg</td></tr>
      <tr><th>Impulse residual</th><td>{summary.impulse_closure_error_n_s:.6e} N s</td></tr>
      <tr><th>Delta-v residual</th><td>{summary.delta_v_closure_error_m_s:.6e} m/s</td></tr>
    </tbody>
  </table>
  <h2>Warnings</h2><ul>{warning_items}</ul>
  <p class="note">Propulsion-only preliminary result. Gravity, drag, steering, attitude, and orbital propagation are outside this model.</p>
  <p><strong>Schema:</strong> {escape(result.schema)}<br>
     <strong>Input hash:</strong> <code>{escape(result.input_hash)}</code><br>
     <strong>Result hash:</strong> <code>{escape(result.result_hash)}</code></p>
</body>
</html>
"""

