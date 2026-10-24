# Geometry export review: work in progress

The geometry export regression tests are an unfinished review checkpoint.
Production geometry/export code has not yet been repaired. The new tests
include both reproduced defects and proposed input/result ownership contracts;
failures of the latter describe missing behavior, rather than regressions from
an earlier implemented contract.

Observed export defects include inconsistent side/cap mesh orientation,
lossy CSV/DXF/STL coordinate serialization, tiny-profile SVG scaling, and
unescaped SVG metadata. Tests also exercise extreme representable geometry,
input validation and retained result replay. The existing DXF test now parses
every vertex and checks exact round trips instead of matching a rounded string.

The targeted baseline on 2026-10-24 was run with:

```text
python -m pytest tests/test_geometry.py tests/test_geometry_export_contract.py -q --tb=short
```

It completed with failures. Pytest reported `40 failed, 10 passed, 6 subtests
passed`; the failure total includes subtest reports and is not a count of
independent defects. The local source-bound runner confirmed unchanged source
during this run. Its report is retained locally as
`output/assurance/power_contract_geometry_export_commit_baseline.json`.

The prior sealed full-suite result of 1264 passing methods across 103 modules
predates these geometry test changes and does not establish a green result for
this work-in-progress checkpoint. Rao/MOC construction and scientific-model
validation also remain open. No CAD manufacturing qualification or motor/feed
hardware maturity is established by these checks.
