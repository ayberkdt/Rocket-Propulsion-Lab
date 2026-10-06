# Core

Shared infrastructure used by the engineering domains.

- `errors.py`: stable error codes with field and diagnostic context.
- `metadata.py`: model identity, assumptions, validity ranges, and warnings.
- `units.py`: explicit affine conversions at application boundaries; domain
  equations remain SI-only.
- `numerics.py`: dependency-free monotonic bisection with convergence evidence.

Domain equations should not duplicate numerical solvers or HTTP-specific errors.

