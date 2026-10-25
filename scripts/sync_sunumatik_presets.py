"""Vendor the Sunumatik presets that the browser workbench uses.

The web UI borrows its palette tokens, chart theme, table preset, focus-layout
pattern, Turkish number formatting, typefaces and duotone icons from the
Sunumatik preset library (https://github.com/ayberkdt/sunumatik). Copies live
under ``src/rocket_propulsion/web/vendor/sunumatik`` so the server works offline
and without the source checkout.

    python scripts/sync_sunumatik_presets.py --source "<sunumatik checkout>"
    python scripts/sync_sunumatik_presets.py --check

``--source`` copies and regenerates every file and rewrites PROVENANCE.json with
the source commit and SHA-256 hashes. ``--check`` re-hashes the vendored files
and exits with status 1 when any file is missing or differs from its record.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import re
import subprocess
import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[1]
VENDOR_ROOT = PROJECT_ROOT / "src" / "rocket_propulsion" / "web" / "vendor" / "sunumatik"
PROVENANCE = VENDOR_ROOT / "PROVENANCE.json"
SOURCE_REPOSITORY = "https://github.com/ayberkdt/sunumatik"

# Stylesheets are vendored as rule subsets: only the selectors the workbench uses are
# kept, so unused deck rules (violin, tooltip, gradient grounds, 1920 px stage type)
# do not enter the page's style budget. A rule survives when every selector in its
# list matches one of the patterns; @media blocks are filtered recursively.
CSS_SUBSETS = {
    "presets/color_themes/palette-library.css": (
        "palette-library.css", [r":root", r'\[data-palette="[^"]+"\]']),
    "presets/charts_icons/chart-theme.css": (
        "chart-theme.css", [
            r":root", r"\.sci-chart", r"\.sci-chart-legend",
            r'\.sci-chart \[data-series-slot="\d"\]',
            r"\.sci-chart \.series-line", r"\.sci-chart \.series-marker",
        ]),
    "presets/color_themes/components/table-presets.css": (
        "table-presets.css", [
            r"\.sci-table", r"\.sci-table (th|td)",
            r"\.sci-table thead th", r"\.sci-table tbody (td|th)",
            r"\.sci-table tbody tr:last-child (td|th)",
            r"\.sci-table--data (td|thead th):not\(:first-child\)",
        ]),
    "presets/core/pano-odak.css": (
        "pano-odak.css", [
            r"\.po-sekmeler", r'\.po-sekmeler \[role="tab"\]',
            r'\.po-sekmeler \[role="tab"\]\[aria-selected="true"\]',
            r'\.po-sekmeler \[role="tab"\]:focus-visible',
            r"\.po-panel", r"\.po-panel\[hidden\]",
        ]),
}

# (source path inside the Sunumatik checkout, vendored path)
COPIED = [
    ("presets/core/pano-odak.mjs", "pano-odak.mjs"),
    ("presets/core/tr-bicim.mjs", "tr-bicim.mjs"),
    ("presets/color_themes/typefaces/source-sans-3/source-sans-3.woff2",
     "typefaces/source-sans-3.woff2"),
    ("presets/color_themes/typefaces/source-sans-3/LICENSE", "typefaces/source-sans-3.LICENSE"),
    ("presets/color_themes/typefaces/jetbrains-mono/jetbrains-mono.woff2",
     "typefaces/jetbrains-mono.woff2"),
    ("presets/color_themes/typefaces/jetbrains-mono/LICENSE", "typefaces/jetbrains-mono.LICENSE"),
]

PALETTE_SOURCE = "presets/color_themes/palette-library.json"
ICON_SOURCES = [
    "presets/charts_icons/icons/science-icons.svg",
    "presets/charts_icons/domain-icons/domain-icons.svg",
]
# One icon per workbench view plus the brand mark; all from the duotone v2 set.
ICONS = [
    "i-rocket", "i-fluid-flow", "i-shock-diamonds", "i-friction", "i-pressure",
    "i-work-energy", "i-entropy", "i-mixture-ratio", "i-propellant-tank",
    "i-nozzle-expansion", "i-nozzle-throat", "i-regenerative-cooling", "i-delta-v-burn",
    "i-thrust-curve", "i-rocket-stages",
]
SERIES_TOKENS = ("data1", "data2", "data3", "data4", "data5", "data6")


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def source_commit(source: Path) -> dict[str, object]:
    """Return the checkout commit and whether any vendored source path is modified."""

    def git(*arguments: str) -> str:
        return subprocess.run(
            ["git", "-C", str(source), *arguments],
            check=True, capture_output=True, text=True,
        ).stdout.strip()

    paths = [item[0] for item in COPIED] + list(CSS_SUBSETS) + [PALETTE_SOURCE, *ICON_SOURCES]
    return {
        "commit": git("rev-parse", "HEAD"),
        "modified_source_paths": [line[3:] for line in git("status", "--porcelain", "--", *paths)
                                  .splitlines() if line],
    }


def relative_luminance(colour: str) -> float:
    """WCAG 2.2 relative luminance of a #RRGGBB colour."""

    def channel(value: int) -> float:
        c = value / 255
        return c / 12.92 if c <= 0.04045 else ((c + 0.055) / 1.055) ** 2.4

    red, green, blue = (int(colour[index:index + 2], 16) for index in (1, 3, 5))
    return 0.2126 * channel(red) + 0.7152 * channel(green) + 0.0722 * channel(blue)


def slim_palettes(source: Path) -> bytes:
    """Keep only palettes with six data series so every chart slot has a palette colour."""

    library = json.loads(source.read_text(encoding="utf-8"))
    palettes = []
    for palette in library["palettes"]:
        tokens = palette.get("tokens", {})
        if not all(token in tokens for token in SERIES_TOKENS):
            continue
        # Mode follows the measured canvas, not the descriptive tags (graticule-slate
        # is tagged "dark" yet has a light canvas).
        dark = relative_luminance(tokens["canvas"]) < 0.18
        palettes.append({
            "id": palette["id"],
            "ad": palette.get("ad") or palette.get("name") or palette["id"],
            "dark": dark,
            "canvas": tokens["canvas"],
            "accent": tokens["accent"],
        })
    body = {"source": PALETTE_SOURCE, "palettes": palettes}
    return (json.dumps(body, ensure_ascii=False, indent=1) + "\n").encode("utf-8")


def css_blocks(text: str) -> list[tuple[str, str]]:
    """Split CSS into top-level (prelude, body) pairs; comments are dropped."""

    text = re.sub(r"/\*.*?\*/", "", text, flags=re.DOTALL)
    blocks, index = [], 0
    while True:
        start = text.find("{", index)
        if start < 0:
            return blocks
        depth, cursor = 1, start + 1
        while depth:
            depth += {"{": 1, "}": -1}.get(text[cursor], 0)
            cursor += 1
        blocks.append((text[index:start].strip(), text[start + 1:cursor - 1]))
        index = cursor


def css_subset(text: str, patterns: list[str]) -> tuple[str, int, int]:
    """Keep rules whose every selector fully matches a pattern; return (css, kept, dropped)."""

    allowed = [re.compile(pattern) for pattern in patterns]
    kept, dropped, output = 0, 0, []
    for prelude, body in css_blocks(text):
        if prelude.startswith("@media"):
            inner, inner_kept, inner_dropped = css_subset(body, patterns)
            kept, dropped = kept + inner_kept, dropped + inner_dropped
            if inner_kept:
                output.append(f"{prelude} {{\n{inner}}}\n")
            continue
        selectors = [" ".join(part.split()) for part in prelude.split(",")]
        if prelude.startswith("@") or not all(
            any(rule.fullmatch(selector) for rule in allowed) for selector in selectors
        ):
            dropped += 1
            continue
        kept += 1
        declarations = "\n".join(f"  {line.strip()}" for line in body.strip().splitlines()
                                 if line.strip())
        output.append(f"{', '.join(selectors)} {{\n{declarations}\n}}\n")
    return "".join(output), kept, dropped


def icon_subset(sources: list[Path]) -> bytes:
    """Extract the listed <symbol> elements into one hidden sprite."""

    symbols: dict[str, str] = {}
    pattern = re.compile(r'<symbol id="([^"]+)".*?</symbol>', re.DOTALL)
    for source in sources:
        for match in pattern.finditer(source.read_text(encoding="utf-8")):
            symbols.setdefault(match.group(1), match.group(0))
    missing = [name for name in ICONS if name not in symbols]
    if missing:
        raise SystemExit(f"Icons missing from the Sunumatik sprites: {', '.join(missing)}")
    lines = [
        '<svg xmlns="http://www.w3.org/2000/svg" style="display:none">',
        "  <!-- Sunumatik duotone v2 icons (science-icons + domain-icons), generated subset. -->",
        *(f"  {symbols[name]}" for name in ICONS),
        "</svg>",
    ]
    return ("\n".join(lines) + "\n").encode("utf-8")


def sync(source: Path) -> None:
    if not (source / "presets").is_dir():
        raise SystemExit(f"Not a Sunumatik checkout: {source}")
    VENDOR_ROOT.mkdir(parents=True, exist_ok=True)
    records = []
    for relative, target in COPIED:
        origin = source / relative
        destination = VENDOR_ROOT / target
        destination.parent.mkdir(parents=True, exist_ok=True)
        destination.write_bytes(origin.read_bytes())
        records.append({"path": target, "source": relative, "source_sha256": sha256(origin),
                        "sha256": sha256(destination), "transform": "copy"})
    for relative, (target, patterns) in CSS_SUBSETS.items():
        origin = source / relative
        css, kept, dropped = css_subset(origin.read_text(encoding="utf-8"), patterns)
        header = (f"/* Sunumatik {relative} — subset vendored by scripts/sync_sunumatik_presets.py.\n"
                  f"   {kept} rules kept, {dropped} unused rules dropped; source hash in PROVENANCE.json. */\n")
        destination = VENDOR_ROOT / target
        destination.write_text(header + css, encoding="utf-8")
        records.append({"path": target, "source": relative, "source_sha256": sha256(origin),
                        "sha256": sha256(destination),
                        "transform": f"rule subset: {kept} kept, {dropped} dropped"})
    generated = [
        ("palettes.json", [PALETTE_SOURCE], slim_palettes(source / PALETTE_SOURCE),
         "palettes with six data series: id, Turkish name, mode, canvas and accent"),
        ("icons.svg", ICON_SOURCES, icon_subset([source / path for path in ICON_SOURCES]),
         f"symbol subset: {', '.join(ICONS)}"),
    ]
    for target, origins, payload, transform in generated:
        destination = VENDOR_ROOT / target
        destination.write_bytes(payload)
        records.append({"path": target, "source": origins,
                        "source_sha256": [sha256(source / path) for path in origins],
                        "sha256": sha256(destination), "transform": transform})
    provenance = {
        "repository": SOURCE_REPOSITORY,
        **source_commit(source),
        "licenses": {
            "typefaces": "SIL OFL 1.1 (license text beside each font)",
            "presets": "Sunumatik, same author as this project",
        },
        "files": records,
    }
    PROVENANCE.write_text(json.dumps(provenance, ensure_ascii=False, indent=2) + "\n",
                          encoding="utf-8")
    print(f"Vendored {len(records)} files from {provenance['commit']} into {VENDOR_ROOT}")
    if provenance["modified_source_paths"]:
        print("WARNING: source paths have uncommitted changes:",
              ", ".join(provenance["modified_source_paths"]))


def check() -> int:
    if not PROVENANCE.is_file():
        print(f"Missing {PROVENANCE}")
        return 1
    records = json.loads(PROVENANCE.read_text(encoding="utf-8"))["files"]
    failures = []
    for record in records:
        path = VENDOR_ROOT / record["path"]
        if not path.is_file():
            failures.append(f"missing {record['path']}")
        elif sha256(path) != record["sha256"]:
            failures.append(f"changed {record['path']}")
    for failure in failures:
        print(failure)
    print(f"{len(records) - len(failures)}/{len(records)} vendored files match PROVENANCE.json")
    return 1 if failures else 0


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--source", type=Path, help="Sunumatik checkout to vendor from.")
    group.add_argument("--check", action="store_true", help="Verify vendored file hashes.")
    arguments = parser.parse_args()
    if arguments.check:
        return check()
    sync(arguments.source.resolve())
    return 0


if __name__ == "__main__":
    sys.exit(main())
