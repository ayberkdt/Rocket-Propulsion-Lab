"""Run a cited provisional propulsion study from JSON and produce reviewable outputs.

References
----------
NASA-SE-APP: Trade Study Report, repeatable inputs, results and data sources.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

from rocket_propulsion.core.errors import DomainError, InputError
from .artifacts import distinct_artifact_paths, write_artifact_pair
from .project import project_from_json
from .report import build_project_html


def main(argv: list[str] | None = None) -> int:
    """Read strict input evidence and write JSON/HTML; return nonzero on invalid data.

    References
    ----------
    NASA-SE-APP: reproducible Trade Study Report and complete data sources.
    """
    parser = argparse.ArgumentParser(description="Cited preliminary UZM 508E propulsion study")
    parser.add_argument("request", type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--report", required=True, type=Path)
    options = parser.parse_args(argv)
    try:
        distinct_artifact_paths((options.request, options.output, options.report))
        study = project_from_json(options.request.read_text(encoding="utf-8"))
        output = json.dumps(study.to_dict(), indent=2, allow_nan=False, ensure_ascii=False)
        html = build_project_html(study)
        write_artifact_pair(((options.output, output + "\n"), (options.report, html)))
    except (DomainError, InputError, OSError) as error:
        print(str(error), file=sys.stderr)
        return 2
    print(f"Study JSON: {options.output.resolve()}")
    print(f"Study report: {options.report.resolve()}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
