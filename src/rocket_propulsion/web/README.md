# Web interface

The browser client is intentionally framework-free. It renders the complete
engineering workbench and sends inputs to the Python JSON API; no engineering
equation is duplicated in JavaScript.

- `index.html` defines semantic forms and result grids.
- `styles.css` provides the responsive visual system.
- `app.js` handles API requests, formatting, and validation feedback.

The standalone Burn workbench resolves an engine operating point through
`/api/v1/engine-operating-point`, then sends the immutable point and burn
definition and optional piecewise-linear throttle schedule to
`/api/v1/burn/simulate`. It displays duration, impulse, ideal one-dimensional
delta-v, tank inventory, mass/thrust histories and closure evidence.
Canonical JSON, printable HTML, workspace persistence and optional Sidera
export all consume the returned artifact rather than recalculating results in
the browser. Variable profiles require the user to explicitly choose the
clearly labelled equivalent-constant approximation; exact mode fails closed.

The interface can later move to a component framework without touching domain
calculations or route contracts.

