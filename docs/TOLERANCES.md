# Numerical tolerance matrix

| Model / evidence | Quantity | Acceptance | Rationale |
|---|---|---:|---|
| Perfect-gas isentropic | tabulated ratios | `1e-12` relative | closed-form double precision |
| Normal shock | tabulated ratios | `1e-12` relative | closed-form double precision |
| NASA species | published 298 K properties | per fixture | rounded source coefficients |
| Thermochemical inverse | `h→T→h`, `s→T→s` | `1e-8` relative | bracketed numerical inverse |
| Local equilibrium | element residual | `1e-8` relative | experimental Gibbs solve |
| External CEA release gate | `Tc`, `c*`, vacuum Isp | `1%` relative declared; measured `≤ 2e-5` | NASA CEA 3.3.4 reference, `cea_gate_points.json` |
| CEA pressure domain | `Tc`, `c*`, vacuum Isp | `2e-5` relative | four pairs at 0.5 and 25 MPa |
| Chemistry-calibrated nozzle | mass flow, design exit pressure | `1e-9` relative | calibrated perfect gas |
| Internal-shock nozzle | mass-flow residual | `1e-8` relative | quasi-1D conservation |
| MOC | exit Mach | `1e-8` absolute | characteristic compatibility |
| MOC | refined contour length | `1%` relative | 2× characteristic count |
| STL round trip | throat/exit area | `1e-9` absolute in tests | generated vertices retain radii |
| Loss budget | thrust/Isp residual | `1e-10` absolute | exact sequential accounting |
| Cooling | bulk energy residual | `1e-8` W | same heat load on both sides |
| Cycle | pump/turbine power residual | `1e-8` W | solved turbine flow |
| Workspace | save/open | exact object and JSON equality | sorted finite serialization |
| Finite-area combustor | `p∞/p_inj`, end Mach, Isp vs CEA | `2e-4`, `2e-3`, `1e-3` relative | declared before comparison; measured 1.5e-5, 1.3e-4, 3.9e-6 |
| Product composition | mole fractions vs CEA | `2e-3` (x ≥ 1e-2), `1e-2` (1e-6 to 1e-2), `5e-2` (1e-10 to 1e-6) relative | declared before comparison; measured 4.6e-6, 1.4e-4, 2.5e-4 |
| Condensed equilibrium | T, p, M, MW, c*, Isp vs CEA | `3e-5` relative | declared after one exploratory comparison |
| Finite-rate nozzle | element conservation | equal to the start composition's residual (5e-10) | equilibrium start |
| Real-fluid transport | NIST WebBook viscosity, conductivity | per region; `2e-3` in dense states with small enhancement | WebBook vs Olchowy–Sengers as printed |
| Transient conduction | energy residual | `1e-9` of the stored-energy change | finite-volume conservation |
| Gray enclosures | summation, reciprocity; energy | `1e-12`; `1e-9` W | exact linear algebra |
| Regenerative cooling | coolant energy balance; wall network | `1e-9` of the heat load; `1e-9` relative | station bisection to 1e-9 |
| Loss chain | delivered Isp by two paths | `1e-12` absolute | same bookkeeping |

Tests use the stricter of the tabulated value and the physical invariant. A
tolerance is never interpreted as permission to extrapolate a model outside its
declared validity range.

## References

Each tolerance row inherits its source from the model it checks; the
sources are listed in [`REFERENCES.md`](REFERENCES.md). In particular:
NASA Glenn coefficients, NASA/TP-2002-211556,
<https://ntrs.nasa.gov/citations/20020085330>; U.S. Standard Atmosphere
1976, <https://ntrs.nasa.gov/citations/19770009539>.
