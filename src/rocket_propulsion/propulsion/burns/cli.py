"""Command-line surface for standalone propulsion burns."""

from __future__ import annotations

import argparse
import json
from collections.abc import Sequence
from pathlib import Path

from rocket_propulsion.api.routes import export_sidera_burn_route, simulate_burn_route
from rocket_propulsion.core.errors import InputError

from .report import build_burn_html_report
from .serialization import burn_result_from_dict, burn_result_from_json, burn_result_to_json
from .tabulated import InterpolationPolicy, tabulated_artifact_from_csv


def build_parser() -> argparse.ArgumentParser:
    """Build the ``rocket-propulsion burn`` command tree."""

    parser = argparse.ArgumentParser(prog="rocket-propulsion burn")
    commands = parser.add_subparsers(dest="command", required=True)

    simulate = commands.add_parser("simulate", help="Run an L0 constant or L1 profiled burn.")
    simulate.add_argument("input", type=Path, help="JSON file using the burn API request shape.")
    simulate.add_argument("--output", type=Path, help="Write the canonical result JSON here.")
    simulate.add_argument("--report", type=Path, help="Write a standalone HTML report here.")
    simulate.add_argument("--json", action="store_true", help="Print the full result JSON.")

    validate = commands.add_parser("validate", help="Validate and reproduce a burn result.")
    validate.add_argument("result", type=Path, help="Canonical burn-result JSON file.")
    validate.add_argument("--json", action="store_true", help="Print a machine-readable status.")

    trace = commands.add_parser(
        "tabulate-trace",
        help="Convert an SI thrust/flow CSV into the isolated transient artifact.",
    )
    trace.add_argument("input", type=Path, help="CSV thrust/flow trace.")
    trace.add_argument("--output", type=Path, required=True, help="Artifact JSON path.")
    trace.add_argument(
        "--source-id",
        help="Stable source identifier; defaults to the input filename.",
    )
    trace.add_argument(
        "--interpolation",
        choices=tuple(item.value for item in InterpolationPolicy),
        default=InterpolationPolicy.LINEAR.value,
    )

    export = commands.add_parser("export", help="Export an exact or reduced constant burn.")
    export.add_argument("result", type=Path, help="Canonical burn-result JSON file.")
    export.add_argument("--target", choices=("sidera",), required=True)
    export.add_argument("--output", type=Path, required=True)
    export.add_argument("--start", type=float, required=True, help="Sidera start time in seconds.")
    export.add_argument(
        "--direction", type=float, nargs=3, required=True, metavar=("X", "Y", "Z")
    )
    export.add_argument("--frame", choices=("inertial", "ric", "vnb"), default="inertial")
    export.add_argument(
        "--mode",
        choices=("exact", "equivalent_constant"),
        default="exact",
        help="Require exact mapping or explicitly permit impulse/mass-preserving reduction.",
    )
    return parser


def _read_text(path: Path) -> str:
    try:
        return path.read_text(encoding="utf-8")
    except OSError as error:
        raise InputError(f"Unable to read {path}.", field="input") from error


def _write_text(path: Path, text: str) -> None:
    try:
        path.write_text(text, encoding="utf-8")
    except OSError as error:
        raise InputError(f"Unable to write {path}.", field="output") from error


def main(arguments: Sequence[str] | None = None) -> int:
    """Run a burn command and return a process-style status code."""

    options = build_parser().parse_args(arguments)
    if options.command == "simulate":
        try:
            payload = json.loads(_read_text(options.input))
        except json.JSONDecodeError as error:
            raise InputError("Burn input is not valid JSON.", field="input") from error
        if not isinstance(payload, dict):
            raise InputError("Burn input root must be an object.", field="input")
        result_payload = simulate_burn_route(payload)
        result = burn_result_from_dict(result_payload)
        result_json = burn_result_to_json(result)
        if options.output is not None:
            _write_text(options.output, result_json + "\n")
        if options.report is not None:
            _write_text(options.report, build_burn_html_report(result))
        if options.json:
            print(result_json)
        else:
            summary = result.summary
            print(
                f"{result.definition.name}: {summary.cutoff_reason.value}; "
                f"duration={summary.duration_s:.9g} s; "
                f"propellant={summary.consumed_propellant_kg:.9g} kg; "
                f"ideal_delta_v={summary.ideal_delta_v_m_s:.9g} m/s"
            )
        return 0 if result.summary.target_achieved else 2

    if options.command == "tabulate-trace":
        artifact = tabulated_artifact_from_csv(
            _read_text(options.input),
            source_id=options.source_id or options.input.name,
            interpolation=InterpolationPolicy(options.interpolation),
        )
        _write_text(options.output, artifact.to_json() + "\n")
        print(
            f"Tabulated {len(artifact.segments)} segments; "
            f"impulse={artifact.delivered_total_impulse_n_s:.9g} N*s; "
            f"propellant={artifact.consumed_propellant_kg:.9g} kg; "
            f"hash={artifact.artifact_hash}"
        )
        return 0

    result = burn_result_from_json(_read_text(options.result))
    if options.command == "validate":
        if options.json:
            print(json.dumps({"valid": True, "result_hash": result.result_hash}, sort_keys=True))
        else:
            print(f"Valid {result.schema}: {result.result_hash}")
        return 0

    artifact = export_sidera_burn_route(
        {
            "burn_result": json.loads(burn_result_to_json(result, indent=None)),
            "t_start_s": options.start,
            "direction": options.direction,
            "frame": options.frame,
            "mode": options.mode,
        }
    )
    _write_text(
        options.output,
        json.dumps(artifact, ensure_ascii=False, allow_nan=False, sort_keys=True, indent=2)
        + "\n",
    )
    print(f"Exported {options.mode} Sidera burn: {options.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

