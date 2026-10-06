# Numerical tolerance matrix

| Model / evidence | Quantity | Acceptance | Rationale |
|---|---|---:|---|
| Perfect-gas isentropic | tabulated ratios | `1e-12` relative | closed-form double precision |
| Normal shock | tabulated ratios | `1e-12` relative | closed-form double precision |
| NASA species | published 298 K properties | per fixture | rounded source coefficients |
| Thermochemical inverse | `h→T→h`, `s→T→s` | `1e-8` relative | bracketed numerical inverse |
| Local equilibrium | element residual | `1e-8` relative | experimental Gibbs solve |
| External CEA release gate | `Tc`, `c*`, Isp | declared per fixture; pending | must be measured against installed CEA |
| Internal-shock nozzle | mass-flow residual | `1e-8` relative | quasi-1D conservation |
| MOC | exit Mach | `1e-8` absolute | characteristic compatibility |
| MOC | refined contour length | `1%` relative | 2× characteristic count |
| STL round trip | throat/exit area | `1e-9` absolute in tests | generated vertices retain radii |
| Loss budget | thrust/Isp residual | `1e-10` absolute | exact sequential accounting |
| Cooling | bulk energy residual | `1e-8` W | same heat load on both sides |
| Cycle | pump/turbine power residual | `1e-8` W | solved turbine flow |
| Workspace | save/open | exact object and JSON equality | sorted finite serialization |

Tests use the stricter of the tabulated value and the physical invariant. A
tolerance is never interpreted as permission to extrapolate a model outside its
declared validity range.

