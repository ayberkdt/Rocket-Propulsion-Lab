# Studies

This package owns reproducibility and engineering decision support, not primary
physics equations.

- `workspace.py` defines the versioned `.rplab.json` schema and v1-to-v2
  migration. Serialization is deterministic and rejects non-finite JSON.
- `comparison.py` produces baseline-relative scalar result differences.
- `sampling.py` provides parameter sweeps, seeded Monte Carlo and Latin
  hypercube sampling, percentile propagation, and correlation sensitivity.
- `report.py` packages a workspace as JSON, tidy CSV, SVG, and printable HTML
  while retaining model versions and warnings.
