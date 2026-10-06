# Architecture

The project follows a simple ports-and-adapters boundary:

```text
Browser UI -> JSON API -> domain packages -> core utilities

domain packages = compressible / thermochemistry / thermodynamics / geometry /
                  propulsion / studies
```

The domain packages never import the API or web layers. They accept plain Python
values and return immutable dataclasses. This keeps calculations reusable from a
notebook, test suite, command-line tool, or future desktop application.

## Design rules

1. Keep equations pure and deterministic.
2. Validate physical domains at module boundaries.
3. Put inverse-relation root finding in `core.numerics`.
4. Use SI units for dimensional inputs and explicit ratio names for dimensionless values.
5. Add reference-value tests with every new relation.
6. Keep the HTTP API versioned under `/api/v1`.
7. Keep propellant records immutable and label their values as screening ranges.
8. Generate plot data in Python; JavaScript only scales and renders returned values.
9. Keep external CEA behind the equilibrium provider protocol; never bundle or import it in ideal-gas modules.
10. Version saved workspaces and supply explicit migrations before changing their shape.

