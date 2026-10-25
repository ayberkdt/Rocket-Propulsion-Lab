"""Run a cited provisional propulsion study from JSON and produce reviewable outputs.

References
----------
.. [1] NASA-SE-APP. NASA, "Systems Engineering Handbook: Appendix"; Trade Study
       Report, repeatable inputs, results and data sources.
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
    """Run a study from a JSON request and write the JSON result and HTML report.

    Parameters
    ----------
    argv : list of str, optional
        Command-line arguments: ``request --output PATH --report PATH``. Defaults to
        ``sys.argv[1:]``.

    Returns
    -------
    int
        0 on success; 2 if the input is invalid or a file operation fails, with the
        message written to standard error.

    References
    ----------
    .. [1] NASA-SE-APP. NASA, "Systems Engineering Handbook: Appendix"; reproducible
           Trade Study Report and complete data sources.
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
